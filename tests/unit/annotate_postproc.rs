    use super::*;

    #[test]
    fn parse_remove_basic() {
        let r = PostProcessor::parse_remove("ID,INFO/DP,FORMAT/GQ").unwrap();
        assert!(r.drop_id);
        assert_eq!(r.info_tags, vec!["DP"]);
        assert_eq!(r.format_tags, vec!["GQ"]);
        assert!(!r.inverse);
    }

    #[test]
    fn parse_remove_inverse() {
        let r = PostProcessor::parse_remove("^INFO/AC,INFO/AN").unwrap();
        assert!(r.inverse);
        assert_eq!(r.info_tags, vec!["AC", "AN"]);
    }

    #[test]
    fn parse_set_id_template_chrom_pos() {
        let t = PostProcessor::parse_set_id("%CHROM:%POS:%REF:%ALT").unwrap();
        assert_eq!(t.tokens.len(), 7);
        assert!(!t.fill_only_missing);
    }

    #[test]
    fn parse_set_id_plus_prefix() {
        let t = PostProcessor::parse_set_id("+chr%CHROM_%POS").unwrap();
        assert!(t.fill_only_missing);
    }

    #[test]
    fn region_filter_cli() {
        let r = RegionFilter::from_cli("1:100-200,1:150-300,2:500").unwrap();
        assert!(r.contains("1", 175));
        assert!(r.contains("1", 250));
        assert!(!r.contains("1", 400));
        assert!(r.contains("2", 500));
        assert!(!r.contains("2", 501));
        assert_eq!(r.by_chr.get("1").unwrap().len(), 1);
    }

    #[test]
    fn region_line_passes_fast() {
        let r = RegionFilter::from_cli("1:100-200").unwrap();
        assert!(r.line_passes("1\t150\trs1\tA\tT\t.\t.\t."));
        assert!(!r.line_passes("1\t300\trs2\tA\tT\t.\t.\t."));
        assert!(!r.line_passes("2\t150\trs3\tA\tT\t.\t.\t."));
    }

    #[test]
    fn parse_output_type_all() {
        assert_eq!(parse_output_type("v").unwrap(), OutputKind::Vcf);
        assert_eq!(parse_output_type("z").unwrap(), OutputKind::VcfGz(6));
        assert_eq!(parse_output_type("z3").unwrap(), OutputKind::VcfGz(3));
        assert_eq!(parse_output_type("b9").unwrap(), OutputKind::Bcf(9));
        assert_eq!(parse_output_type("u").unwrap(), OutputKind::Bcf(0));
    }

    #[test]
    fn pair_logic_parse() {
        assert_eq!(PairLogic::parse("some").unwrap(), PairLogic::Some_);
        assert_eq!(PairLogic::parse("exact").unwrap(), PairLogic::Exact);
        assert!(PairLogic::parse("foo").is_err());
    }

    #[test]
    fn render_id_template_basic() {
        let t = PostProcessor::parse_set_id("%CHROM_%POS_%REF/%ALT").unwrap();
        let out = render_id_template(&t, "1", "100", "A", "T", ".", "DP=10");
        assert_eq!(out, "1_100_A/T");
    }

    #[test]
    fn render_id_info_tag() {
        let t = PostProcessor::parse_set_id("%CHROM:%POS:%INFO/RSID").unwrap();
        let out = render_id_template(&t, "1", "100", "A", "T", ".", "RSID=rs9999;DP=10");
        assert_eq!(out, "1:100:rs9999");
    }

    #[test]
    fn filter_info_drops_listed() {
        let out = filter_info("DP=10;AC=2;AN=4", &["DP".into(), "AN".into()], false, false);
        assert_eq!(out, "AC=2");
    }

    #[test]
    fn rename_info_keys_simple() {
        let mut m = FxHashMap::default();
        m.insert("DP".to_string(), "DEPTH".to_string());
        assert_eq!(rename_info_keys("DP=10;AC=2", &m), "DEPTH=10;AC=2");
    }

    #[test]
    fn header_remove_info_dp() {
        let rm = PostProcessor::parse_remove("INFO/DP").unwrap();
        let opts = HeaderOptions {
            no_version: true, extra_header_lines: &[], remove: Some(&rm),
            rename_chrs: None, rename_annots: None, mark_sites: None, set_id: false,
            samples_keep: None, version_line: None,
        };
        let h = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"AC\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let out = apply_to_header(h, &opts);
        assert!(!out.iter().any(|l| l.contains("ID=DP,")));
        assert!(out.iter().any(|l| l.contains("ID=AC,")));
    }

    #[test]
    fn rename_chrs_in_contig() {
        let mut m = FxHashMap::default();
        m.insert("1".to_string(), "chr1".to_string());
        let opts = HeaderOptions {
            no_version: true, extra_header_lines: &[], remove: None,
            rename_chrs: Some(&m), rename_annots: None, mark_sites: None, set_id: false,
            samples_keep: None, version_line: None,
        };
        let h = vec![
            "##contig=<ID=1,length=249250621>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let out = apply_to_header(h, &opts);
        assert!(out.iter().any(|l| l.contains("ID=chr1,")));
    }
