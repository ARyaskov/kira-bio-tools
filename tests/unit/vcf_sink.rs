use super::*;

#[test]
fn regression_pass_header_is_first_filter_and_preserves_description() {
    let dir = tempfile::tempdir().unwrap();
    for pass in [None, Some("##FILTER=<ID=PASS,Description=\"Accepted\">")] {
        let mut headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FILTER=<ID=LowQual,Description=\"Low quality\">".to_string(),
        ];
        if let Some(pass) = pass {
            headers.push(pass.to_string());
        }
        headers.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string());
        for (kind, extension) in [
            (OutputKind::Vcf, "vcf"),
            (OutputKind::VcfGz(1), "vcf.gz"),
            (OutputKind::Bcf(1), "bcf"),
        ] {
            let path = dir.path().join(extension);
            let mut sink = VcfSink::open(Some(&path), kind, &headers).unwrap();
            sink.write_header(&headers).unwrap();
            sink.finish().unwrap();
            let mut reader = crate::vcf::VcfReader::open(&path).unwrap();
            let header = reader.header().unwrap();
            let passes: Vec<_> = header
                .iter()
                .filter(|h| h.starts_with("##FILTER=<ID=PASS,"))
                .collect();
            assert_eq!(passes.len(), 1);
            assert_eq!(
                passes[0].as_str(),
                pass.unwrap_or("##FILTER=<ID=PASS,Description=\"All filters passed\">")
            );
            if !matches!(kind, OutputKind::Bcf(_)) {
                assert_eq!(header[1], *passes[0]);
            }
        }
    }
}

#[test]
fn regression_streamed_headers_and_records_survive_all_output_formats() {
    let dir = tempfile::tempdir().unwrap();
    let text = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tC\t.\tPASS\tEND=250\n";
    for (kind, extension) in [
        (OutputKind::Vcf, "vcf"),
        (OutputKind::VcfGz(1), "vcf.gz"),
        (OutputKind::Bcf(0), "u.bcf"),
        (OutputKind::Bcf(1), "bcf"),
    ] {
        let path = dir.path().join(extension);
        let mut sink = StreamingVcfSink::new(Some(&path), kind);
        for chunk in text.as_bytes().chunks(7) {
            sink.write_all(chunk).unwrap();
        }
        sink.finish().unwrap();
        let mut reader = crate::vcf::VcfReader::open(&path).unwrap();
        assert!(
            reader
                .header()
                .unwrap()
                .iter()
                .any(|h| h.starts_with("#CHROM"))
        );
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.info, "END=250");
        assert!(reader.next_record().unwrap().is_none());
        if matches!(kind, OutputKind::VcfGz(_) | OutputKind::Bcf(1..)) {
            write_output_index(Some(&path), "csi").unwrap();
            let mut indexed = crate::csi::IndexedVcfReader::open(&path).unwrap();
            let mut hits = 0;
            indexed
                .query("chr1", 150, 150, |_| {
                    hits += 1;
                    Ok(true)
                })
                .unwrap();
            assert_eq!(hits, 1);
        }
    }
}

#[test]
fn regression_missing_header_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("invalid.vcf");
    let mut sink = StreamingVcfSink::new(Some(&path), OutputKind::Vcf);
    assert!(sink.write_all(b"chr1\t10\t.\tA\tC\t.\t.\t.\n").is_err());
    assert!(sink.finish().is_err());
    assert!(write_output_index(None, "csi").is_err());
}
