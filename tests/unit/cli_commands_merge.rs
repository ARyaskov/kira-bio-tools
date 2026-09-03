    use super::*;

    #[test]
    fn info_rules_defaults_and_parse() {
        let m = parse_info_rules(None).unwrap();
        assert!(matches!(m["DP"], InfoRule::Sum));
        assert!(matches!(m["DP4"], InfoRule::Sum));
        let m = parse_info_rules(Some("AN:sum,AC:sum,AF:join")).unwrap();
        assert!(matches!(m["AF"], InfoRule::Join));
        assert!(parse_info_rules(Some("-")).unwrap().is_empty());
        assert!(parse_info_rules(Some("DP:bogus")).is_err());
    }

    #[test]
    fn filters_union_and_exclude() {
        assert_eq!(merge_filters(["PASS", "PASS"].into_iter(), FilterLogic::Union), "PASS");
        assert_eq!(merge_filters(["PASS", "q10"].into_iter(), FilterLogic::Union), "q10");
        assert_eq!(merge_filters(["PASS", "q10"].into_iter(), FilterLogic::Exclude), "PASS");
        assert_eq!(merge_filters(["a", "b;c"].into_iter(), FilterLogic::Exclude), "a;b;c");
        assert_eq!(merge_filters([".", "PASS"].into_iter(), FilterLogic::Union), "PASS");
    }

    #[test]
    fn merge_mode_uses_shared_variant_types() {
        let mut d = ContigDict::new();
        let snp = parse_record("1\t10\t.\tA\tC\t.\t.\t.".into(), &mut d).unwrap();
        let snp2 = parse_record("1\t10\t.\tA\tG\t.\t.\t.".into(), &mut d).unwrap();
        let indel = parse_record("1\t10\t.\tA\tAT\t.\t.\t.".into(), &mut d).unwrap();
        assert!(parse_record("1\tx\t.\tA\tC\t.\t.\t.".into(), &mut d).is_err());
        assert!(matches_merge_mode(&snp2, &snp, MergeMode::Snps));
        assert!(!matches_merge_mode(&indel, &snp, MergeMode::Snps));
        assert!(!matches_merge_mode(&indel, &snp, MergeMode::Both));
        assert!(matches_merge_mode(&indel, &snp, MergeMode::All));
        assert!(!matches_merge_mode(&snp2, &snp, MergeMode::None_));
    }

    #[test]
    fn vector_rules() {
        assert_eq!(vector_fold(&["1,2", "3,4"], |a, b| a + b), "4,6");
        assert_eq!(vector_fold(&["5", "3"], f64::min), "3");
        assert_eq!(normalize_num(2.5), "2.5");
        assert_eq!(normalize_num(3.0), "3");
    }
