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

    #[test]
    fn tabix_presets_parse_intervals() {
        assert_eq!(tabix_line_interval(&TabixConf::BED, "chr1\t100\t200\tx"), Some(("chr1", 100, 200)));
        assert_eq!(tabix_line_interval(&TabixConf::GFF, "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g"), Some(("chr1", 100, 200)));
        assert_eq!(tabix_line_interval(&TabixConf::SAM, "r1\t0\tchr1\t100\t60\t10M2D5M3I\t*\t0\t0\tACGT\t*"), Some(("chr1", 99, 116)));
        assert_eq!(tabix_line_interval(&TabixConf::VCF, "1\t100\t.\tACGT\tA\t.\t.\tEND=250"), Some(("1", 99, 250)));
        let custom = TabixConf { preset: TBX_GENERIC, sc: 2, bc: 3, ec: 4, meta: b'#', skip: 1 };
        assert_eq!(tabix_line_interval(&custom, "x\tchrZ\t5\t9"), Some(("chrZ", 4, 9)));
        assert_eq!(tabix_line_interval(&TabixConf::BED, "chr1\tabc\t200"), None);
        assert_eq!(TabixConf::from_extension(Path::new("a/b.gff3.gz")), Some(TabixConf::GFF));
        assert_eq!(TabixConf::from_extension(Path::new("x.bed")), Some(TabixConf::BED));
        assert_eq!(TabixConf::from_extension(Path::new("x.txt.gz")), None);
        // CSI depth grows with the declared sequence length, 100 Gbp when none is declared.
        assert_eq!(adjust_n_lvls(14, 6, 100 * 1024 * 1024 * 1024), 8);
        assert_eq!(adjust_n_lvls(14, 6, 249_250_621), 6);
        assert_eq!(header_ref_len(&TabixConf::VCF, "##contig=<ID=1,length=249250621>"), Some(249250621));
        assert_eq!(header_ref_len(&TabixConf::SAM, "@SQ\tSN:chr1\tLN:5000"), Some(5000));
    }

    /// A BED file indexed with the bed preset answers overlap queries by its
    /// own columns, and the index header records the layout.
    #[test]
    fn bed_index_then_query() {
        let bed = tmp("regions.bed.gz");
        write_bgzf_vcf(&bed, &[
            "#chrom\tstart\tend\tname".into(),
            "chr1\t100\t200\ta".into(),
            "chr1\t150\t400\tb".into(),
            "chr1\t1000\t1001\tc".into(),
            "chr2\t5\t10\td".into(),
        ]);
        let idx = build_tabix_index(&bed, &TabixConf::BED, IndexKind::Tbi, None).unwrap();
        let h = idx.tabix_header().unwrap();
        assert_eq!((h.format, h.col_seq, h.col_beg, h.col_end), (TBX_UCSC, 1, 2, 3));
        assert_eq!(h.names, vec!["chr1".to_string(), "chr2".to_string()]);
        let mut ip = bed.as_os_str().to_os_string();
        ip.push(".tbi");
        idx.save(&ip).unwrap();
        let mut r = IndexedVcfReader::open(&bed).unwrap();
        assert_eq!(r.headers(), ["#chrom\tstart\tend\tname".to_string()]);
        let mut hits = Vec::new();
        r.query("chr1", 190, 210, |l| { hits.push(l.to_string()); Ok(true) }).unwrap();
        assert_eq!(hits, vec!["chr1\t100\t200\ta".to_string(), "chr1\t150\t400\tb".to_string()]);
        hits.clear();
        r.query("chr1", 1001, 1001, |l| { hits.push(l.to_string()); Ok(true) }).unwrap();
        assert_eq!(hits, vec!["chr1\t1000\t1001\tc".to_string()]);
        hits.clear();
        r.query("chr2", 1, 5, |l| { hits.push(l.to_string()); Ok(true) }).unwrap();
        assert!(hits.is_empty());
        let _ = std::fs::remove_file(&bed);
        let _ = std::fs::remove_file(&ip);
    }

    /// Overlapping and adjacent regions read the file once: every record is
    /// emitted a single time, in file order, whichever regions cover it.
    #[test]
    fn multi_region_query_reads_each_record_once() {
        let vcf = tmp("multi.vcf.gz");
        let mut lines = header();
        for pos in (100..=100_000).step_by(100) {
            lines.push(format!("scaffold_9\t{pos}\t.\tACGT\tA\t.\t.\t."));
        }
        write_bgzf_vcf(&vcf, &lines);
        let idx = build_index_in_memory(&vcf, IndexKind::Csi, None).unwrap();
        let mut ip = vcf.as_os_str().to_os_string();
        ip.push(".csi");
        idx.save(&ip).unwrap();
        let regions = crate::regions::RegionSet::from_cli(
            "scaffold_9:150-260,scaffold_9:261-400,scaffold_9:390-450,scaffold_9:50000-50003,scaffold_9:99990-200000",
        )
        .unwrap();
        let mut seen: Vec<u32> = Vec::new();
        regions.stream_with_index(&vcf, 1, |l| {
            seen.push(l.split('\t').nth(1).unwrap().parse().unwrap());
            Ok(())
        }).unwrap();
        // 200 overlaps 150-260 by its 4-base span (200..203); 100 (100..103) does not reach 150.
        assert_eq!(seen, vec![200, 300, 400, 50000, 100000]);
        let mut r = IndexedVcfReader::open(&vcf).unwrap();
        let mut n = 0;
        r.query_multi("scaffold_9", &[(1, 1000), (1001, 100_000)], |_| { n += 1; Ok(true) }).unwrap();
        assert_eq!(n, 1000);
        let _ = std::fs::remove_file(&vcf);
        let _ = std::fs::remove_file(&ip);
    }
