    use super::*;

    fn rec(s: &str) -> Rec {
        Rec::from_line(s.to_string()).unwrap()
    }

    #[test]
    fn collapse_modes() {
        let a = rec("1\t10\trs1\tA\tC\t.\t.\t.");
        let b = rec("1\t10\trs2\tA\tG\t.\t.\t.");
        let c = rec("1\t10\trs1\tA\tAT\t.\t.\t.");
        let d = rec("1\t10\t.\tA\tC,G\t.\t.\t.");
        assert!(!same_site(&a, &b, Collapse::None));
        assert!(same_site(&a, &b, Collapse::Snps));
        assert!(same_site(&a, &b, Collapse::All));
        assert!(!same_site(&a, &c, Collapse::Snps));
        assert!(!same_site(&a, &c, Collapse::Both));
        assert!(same_site(&a, &c, Collapse::Id));
        assert!(same_site(&a, &d, Collapse::Some_));
        assert!(!same_site(&c, &d, Collapse::Some_));
    }

    #[test]
    fn nfiles_specs() {
        assert!(matches_nfiles(&[true, true], &parse_nfiles(Some("=2"), 2, false).unwrap()));
        assert!(!matches_nfiles(&[true, false], &parse_nfiles(Some("=2"), 2, false).unwrap()));
        assert!(matches_nfiles(&[true, false], &parse_nfiles(Some("+1"), 2, false).unwrap()));
        assert!(matches_nfiles(&[true, false], &parse_nfiles(Some("-1"), 2, false).unwrap()));
        assert!(!matches_nfiles(&[true, true], &parse_nfiles(Some("-1"), 2, false).unwrap()));
        assert!(matches_nfiles(&[true, false, true], &parse_nfiles(Some("~101"), 3, false).unwrap()));
        assert!(!matches_nfiles(&[true, true, true], &parse_nfiles(Some("~101"), 3, false).unwrap()));
        assert!(matches_nfiles(&[true, false], &parse_nfiles(None, 2, true).unwrap()));
        assert!(!matches_nfiles(&[true, true], &parse_nfiles(None, 2, true).unwrap()));
        assert!(parse_nfiles(Some("~1"), 2, false).is_err());
    }
