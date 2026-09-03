    use super::*;

    #[test]
    fn dup_window_modes() {
        let mut w = DupWindow::default();
        assert!(!w.check("1\t10\t.\tA\tC\t.\t.\t.", RmDups::Exact));
        assert!(w.check("1\t10\t.\tA\tC\t.\t.\t.", RmDups::Exact));
        assert!(!w.check("1\t10\t.\tA\tG\t.\t.\t.", RmDups::Exact));
        let mut w = DupWindow::default();
        assert!(!w.check("1\t10\t.\tA\tC\t.\t.\t.", RmDups::Snps));
        assert!(w.check("1\t10\t.\tA\tG\t.\t.\t.", RmDups::Snps));
        assert!(!w.check("1\t10\t.\tA\tAT\t.\t.\t.", RmDups::Snps));
        assert!(!w.check("1\t11\t.\tA\tC\t.\t.\t.", RmDups::All));
        assert!(w.check("1\t11\t.\tA\tAT\t.\t.\t.", RmDups::All));
    }

    #[test]
    fn rm_dups_parse() {
        assert_eq!(RmDups::parse(None, true).unwrap(), RmDups::Exact);
        assert_eq!(RmDups::parse(None, false).unwrap(), RmDups::None);
        assert_eq!(RmDups::parse(Some("both"), false).unwrap(), RmDups::Both);
        assert!(RmDups::parse(Some("x"), false).is_err());
    }
