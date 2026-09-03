    use super::*;
    use crate::csi::IndexedVcfReader;

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("kira_csi_{}_{}_{name}", std::process::id(), std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()))
    }

    fn write_bgzf_vcf(path: &Path, lines: &[String]) {
        let mut w = crate::bgzf::BgzfWriter::create(path).unwrap();
        for l in lines {
            w.write_all(l.as_bytes()).unwrap();
            w.write_all(b"\n").unwrap();
        }
        w.finish().unwrap();
    }

    fn header() -> Vec<String> {
        vec![
            "##fileformat=VCFv4.2".into(),
            "##contig=<ID=scaffold_9,length=500000>".into(),
            "##contig=<ID=scaffold_10,length=500000>".into(),
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"x\">".into(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".into(),
        ]
    }

    #[test]
    fn interval_from_line() {
        assert_eq!(vcf_line_interval("1\t100\t.\tA\tT\t.\t.\t."), Some(("1", 99, 100)));
        assert_eq!(vcf_line_interval("1\t100\t.\tACGT\tA\t.\t.\tDP=3"), Some(("1", 99, 103)));
        assert_eq!(vcf_line_interval("1\t100\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=250"), Some(("1", 99, 250)));
        assert_eq!(vcf_line_interval("1\t100\t.\tA\tT"), Some(("1", 99, 100)));
    }

    #[test]
    fn index_then_query_deletion_span_and_non_human_contigs() {
        let path = tmp("a.vcf.gz");
        let mut lines = header();
        // scaffold_10 first: index order follows the data, not the header.
        lines.push("scaffold_10\t50\t.\tA\tT\t.\t.\t.".into());
        lines.push("scaffold_10\t100\t.\tACGTACGTAC\tA\t.\t.\t.".into()); // spans 100..109
        lines.push("scaffold_10\t200\t.\tA\t<DEL>\t.\t.\tEND=40000".into());
        for p in (1..=30000u32).step_by(7) {
            lines.push(format!("scaffold_9\t{p}\t.\tG\tC\t.\t.\t."));
        }
        write_bgzf_vcf(&path, &lines);
        for kind in [IndexKind::Csi, IndexKind::Tbi] {
            let ip = tmp(if kind == IndexKind::Csi { "a.csi" } else { "a.tbi" });
            let idx = build_index(&path, &ip, kind, None).unwrap();
            assert_eq!(idx.names(), &["scaffold_10", "scaffold_9"]);
            assert_eq!(idx.n_records(0), Some(3));
            assert_eq!(idx.n_records(1), Some((1..=30000u32).step_by(7).count() as u64));

            let mut r = IndexedVcfReader::open_with_index(&path, &ip).unwrap();
            let mut got = Vec::new();
            r.query("scaffold_10", 105, 106, |l| { got.push(l.to_string()); Ok(true) }).unwrap();
            assert_eq!(got.len(), 1, "{kind:?}: deletion spanning the query must be found");
            assert!(got[0].starts_with("scaffold_10\t100\t"));

            got.clear();
            r.query("scaffold_10", 30000, 30001, |l| { got.push(l.to_string()); Ok(true) }).unwrap();
            assert_eq!(got.len(), 1, "{kind:?}: INFO/END span must be found");
            assert!(got[0].contains("<DEL>"));

            got.clear();
            r.query("scaffold_9", 20000, 20100, |l| { got.push(l.to_string()); Ok(true) }).unwrap();
            let expect = (1..=30000u32).step_by(7).filter(|p| (20000..=20100).contains(p)).count();
            assert_eq!(got.len(), expect);

            got.clear();
            r.query("chr1", 1, 10, |l| { got.push(l.to_string()); Ok(true) }).unwrap();
            assert!(got.is_empty());
            let _ = std::fs::remove_file(&ip);
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn unsorted_input_is_rejected() {
        let path = tmp("u.vcf.gz");
        let mut lines = header();
        lines.push("scaffold_9\t500\t.\tA\tT\t.\t.\t.".into());
        lines.push("scaffold_9\t100\t.\tA\tT\t.\t.\t.".into());
        write_bgzf_vcf(&path, &lines);
        assert!(build_index_in_memory(&path, IndexKind::Csi, None).is_err());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn plain_text_cannot_be_indexed() {
        let path = tmp("p.vcf");
        std::fs::write(&path, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tT\t.\t.\t.\n").unwrap();
        assert!(build_index_in_memory(&path, IndexKind::Csi, None).is_err());
        let _ = std::fs::remove_file(&path);
    }
