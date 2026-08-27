    use super::*;

    #[test]
    fn test_parse_modes() {
        let (mode, rest) = AnnotateMode::parse("TAG");
        assert!(mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse("+TAG");
        assert!(mode.replace_missing);
        assert!(!mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse("-TAG");
        assert!(mode.replace_non_missing);
        assert!(!mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse("=TAG");
        assert!(mode.set_or_append);
        assert!(!mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse(".TAG");
        assert!(mode.carry_over_missing);
        assert!(mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse(".+TAG");
        assert!(mode.carry_over_missing);
        assert!(mode.replace_missing);
        assert!(!mode.replace_all);
        assert_eq!(rest, "TAG");

        let (mode, rest) = AnnotateMode::parse("~ID");
        assert!(mode.match_value);
        assert_eq!(rest, "ID");
    }

    #[test]
    fn test_should_transfer() {
        let mode_tag = AnnotateMode::default_mode();
        assert!(mode_tag.should_transfer(false, false, false));
        assert!(mode_tag.should_transfer(false, true, false));
        assert!(mode_tag.should_transfer(false, true, true));
        assert!(!mode_tag.should_transfer(true, false, false));

        let (mode_plus, _) = AnnotateMode::parse("+TAG");
        assert!(mode_plus.should_transfer(false, false, false));
        assert!(!mode_plus.should_transfer(false, true, false));
        assert!(mode_plus.should_transfer(false, true, true));

        let (mode_minus, _) = AnnotateMode::parse("-TAG");
        assert!(!mode_minus.should_transfer(false, false, false));
        assert!(mode_minus.should_transfer(false, true, false));
        assert!(!mode_minus.should_transfer(false, true, true));

        let (mode_dot, _) = AnnotateMode::parse(".TAG");
        assert!(mode_dot.should_transfer(true, false, false));
        assert!(mode_dot.should_transfer(true, true, false));
    }
