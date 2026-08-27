    use super::*;

    #[test]
    fn parse_annotate_default() {
        let s = AnnotateSpec::parse(None).unwrap();
        assert!(s.fmt_ad);
        assert!(s.fmt_dp);
        assert!(s.fmt_pl);
        assert!(!s.fmt_qs);
    }

    #[test]
    fn parse_annotate_selective() {
        let s = AnnotateSpec::parse(Some("FORMAT/QS,FORMAT/SP,INFO/AD")).unwrap();
        assert!(s.fmt_qs);
        assert!(s.fmt_sp);
        assert!(s.info_ad);
        assert!(!s.fmt_ad);
        assert!(s.fmt_pl);
    }

    #[test]
    fn preset_ont_baq_disabled() {
        let p = PresetConfig::parse("ont").unwrap();
        assert_eq!(p.no_baq, Some(true));
        assert_eq!(p.min_mq, Some(7));
    }

    #[test]
    fn flag_filter_require_paired() {
        let f = FlagFilters::from_args(None, None, Some(0x1), None);
        assert!(f.passes(0x1 | 0x40));
        assert!(!f.passes(0x40));
    }

    #[test]
    fn flag_filter_exclude_dup() {
        let f = FlagFilters::from_args(Some(0x400), None, None, None);
        assert!(f.passes(0x1));
        assert!(!f.passes(0x400 | 0x1));
    }

    #[test]
    fn parse_samples_filter_cli() {
        let r = parse_samples_filter(Some("S1,S2,S3"), None).unwrap().unwrap();
        assert_eq!(r, vec!["S1", "S2", "S3"]);
    }
