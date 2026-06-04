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
