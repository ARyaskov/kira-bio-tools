    use super::*;

    const SPECS: &[OptSpec] = &[
        opt(Some('O'), "output-type", true),
        opt(Some('G'), "GTs-only", true),
        opt(Some('i'), "ignore-homref", false),
        opt(Some('I'), "skip-indels", false),
        opt(None, "AF-dflt", true),
        opt(None, "AF-file", true),
    ];

    fn toks(s: &str) -> Vec<String> {
        s.split_whitespace().map(str::to_string).collect()
    }

    #[test]
    fn attached_separate_and_long_forms() {
        let p = parse(
            &toks("-Or -G 30 --AF-dflt 0.4 --AF-f=x.tab in.vcf -iI"),
            SPECS,
        )
        .unwrap();
        assert_eq!(
            p.opts,
            vec![
                ("output-type", Some("r".into())),
                ("GTs-only", Some("30".into())),
                ("AF-dflt", Some("0.4".into())),
                ("AF-file", Some("x.tab".into())),
                ("ignore-homref", None),
                ("skip-indels", None),
            ]
        );
        assert_eq!(p.positional, vec!["in.vcf"]);
    }

    #[test]
    fn errors_and_dashdash() {
        assert!(parse(&toks("-Z"), SPECS).is_err());
        assert!(parse(&toks("--AF"), SPECS).is_err());
        assert!(parse(&toks("-G"), SPECS).is_err());
        let p = parse(&toks("-- -Or"), SPECS).unwrap();
        assert_eq!(p.positional, vec!["-Or"]);
        let p = parse(&toks("-"), SPECS).unwrap();
        assert_eq!(p.positional, vec!["-"]);
    }
