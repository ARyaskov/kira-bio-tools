use super::*;
use std::io::Write;

fn write_tmp_vcf(name: &str, body: &str) -> std::path::PathBuf {
    let mut path = std::env::temp_dir();
    path.push(format!(
        "kira_ktile_writer_test_{}_{}.vcf",
        name,
        std::process::id()
    ));
    let mut f = std::fs::File::create(&path).unwrap();
    f.write_all(body.as_bytes()).unwrap();
    path
}

#[test]
fn roundtrip_tiny_vcf() {
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\trs1\tA\tT\t.\t.\t.\n1\t200\trs2\tC\tG,A\t.\t.\t.\n2\t50\t.\tA\tAT\t.\t.\t.\n";
    let in_path = write_tmp_vcf("roundtrip", vcf);
    let out_path = in_path.with_extension("ktile");

    let stats = write_ktile_from_vcf(&in_path, &out_path).unwrap();
    assert_eq!(stats.n_records, 3);
    assert!(stats.bytes_written > KtileHeader::SIZE as u64);

    let reader = super::super::reader::KtileReader::open(&out_path).unwrap();
    assert_eq!(reader.n_records(), 3);
    assert_eq!(reader.header().version, KTILE_LEGACY_VERSION);
    assert!(!reader.header().has_elias_fano_offsets());
    assert_eq!(reader.line_owned(0), "1\t100\trs1\tA\tT\t.\t.\t.");
    assert_eq!(reader.line_owned(1), "1\t200\trs2\tC\tG,A\t.\t.\t.");
    assert_eq!(reader.line_owned(2), "2\t50\t.\tA\tAT\t.\t.\t.");
    assert_eq!(reader.position(0), 100);
    assert_eq!(reader.position(2), 50);

    assert!(reader.has_ref_alt_columns());
    assert_eq!(reader.ref_slice(0).as_deref(), Some(b"A".as_slice()));
    assert_eq!(reader.alt_slice(0).as_deref(), Some(b"T".as_slice()));
    assert_eq!(reader.ref_slice(1).as_deref(), Some(b"C".as_slice()));
    assert_eq!(reader.alt_slice(1).as_deref(), Some(b"G,A".as_slice()));
    assert_eq!(reader.ref_slice(2).as_deref(), Some(b"A".as_slice()));
    assert_eq!(reader.alt_slice(2).as_deref(), Some(b"AT".as_slice()));

    let _ = std::fs::remove_file(&in_path);
    let _ = std::fs::remove_file(&out_path);
}

#[test]
fn compact_offsets_round_trip_both_pool_modes_against_input_lines() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("records.vcf");
    let headers = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let lines: Vec<_> = (0..2051)
        .map(|i| {
            let chrom = i % 2 + 1;
            let (reference, alternate) = if i % 3 == 0 {
                ("AT", "A,C")
            } else {
                ("A", "G")
            };
            format!(
                "{chrom}\t{}\tmarker_λ_{i}\t{reference}\t{alternate}\t.\t.\tNOTE={}",
                i + 1,
                "x".repeat(i % 37)
            )
        })
        .collect();
    std::fs::write(&input, format!("{headers}\n{}\n", lines.join("\n"))).unwrap();
    for compressed in [false, true] {
        let output = dir.path().join(format!("records-{compressed}.ktile"));
        let stats = write_ktile_from_vcf_with(
            &input,
            &output,
            KtileWriteOptions {
                compressed,
                lines_per_chunk: 1024,
                deflate_level: 1,
            },
        )
        .unwrap();
        let reader = super::super::reader::KtileReader::open(&output).unwrap();
        let header = reader.header();
        assert_eq!(header.version, KTILE_VERSION);
        assert!(header.has_elias_fano_offsets());
        assert_eq!(header.has_compressed_pool(), compressed);
        assert_eq!(reader.headers_block(), headers);
        assert_eq!(reader.n_records(), lines.len());
        assert_eq!(stats.n_records, lines.len());
        assert_eq!(
            stats.bytes_written,
            std::fs::metadata(&output).unwrap().len()
        );
        assert!(header.chr_ids_off - header.line_offsets_off < ((lines.len() + 1) * 8) as u64);
        let mut batch = crate::annotate::cpu_v2::ReadBatch::with_capacity(0, lines.len());
        for (i, line) in lines.iter().enumerate() {
            assert_eq!(
                reader.line_owned(i),
                *line,
                "compressed={compressed}, row={i}"
            );
            assert_eq!(reader.chr_id(i), (i % 2 + 1) as u32);
            assert_eq!(reader.position(i), (i + 1) as u32);
            reader.push_line_into_batch(i, &mut batch);
        }
        assert_eq!(
            batch.iter().collect::<Vec<_>>(),
            lines.iter().map(String::as_str).collect::<Vec<_>>()
        );
        for i in [2050, 1024, 1023, 0, 2048, 1] {
            assert_eq!(reader.line_owned(i), lines[i]);
            let expected_ref = if i % 3 == 0 {
                b"AT".as_slice()
            } else {
                b"A".as_slice()
            };
            let expected_alt = if i % 3 == 0 {
                b"A,C".as_slice()
            } else {
                b"G".as_slice()
            };
            assert_eq!(reader.ref_slice(i).as_deref(), Some(expected_ref));
            assert_eq!(reader.alt_slice(i).as_deref(), Some(expected_alt));
        }
    }
}

#[test]
fn compact_offset_selection_keeps_tiny_and_invalid_sequences_raw() {
    assert!(encode_compact_offsets(&[0, 100, 200]).is_none());
    assert!(encode_compact_offsets(&[0, u64::MAX]).is_none());
    assert!(encode_compact_offsets(&vec![1; 2048]).is_none());
    let offsets: Vec<_> = (0..2048).map(|i| i * 71).collect();
    let bytes = encode_compact_offsets(&offsets).unwrap();
    assert!(bytes.len() < std::mem::size_of_val(offsets.as_slice()));
    let mut cursor = 0;
    let decoded = EliasFano::read_from(&bytes, &mut cursor).unwrap();
    assert_eq!(cursor, bytes.len());
    assert_eq!(decoded.len(), offsets.len());
    for (i, expected) in offsets.into_iter().enumerate() {
        assert_eq!(decoded.get(i), expected);
    }
}
