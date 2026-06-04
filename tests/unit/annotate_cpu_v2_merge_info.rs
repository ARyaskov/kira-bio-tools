    use super::*;

    fn parse_to_pairs(info: &str) -> Vec<(String, String)> {
        parse_existing_info(info)
            .into_iter()
            .map(|(k, v)| (k, v))
            .collect()
    }

    #[test]
    fn parse_info_empty_and_dot() {
        assert!(parse_existing_info("").is_empty());
        assert!(parse_existing_info(".").is_empty());
    }

    #[test]
    fn parse_info_single_kv() {
        assert_eq!(parse_to_pairs("DP=42"), vec![("DP".into(), "42".into())]);
    }

    #[test]
    fn parse_info_flag() {
        assert_eq!(parse_to_pairs("SOMATIC"), vec![("SOMATIC".into(), "".into())]);
    }

    #[test]
    fn parse_info_mixed() {
        assert_eq!(
            parse_to_pairs("DP=42;AF=0.5;SOMATIC;CLN=A,B"),
            vec![
                ("DP".into(), "42".into()),
                ("AF".into(), "0.5".into()),
                ("SOMATIC".into(), "".into()),
                ("CLN".into(), "A,B".into()),
            ]
        );
    }

    #[test]
    fn parse_info_value_contains_equals() {
        // Split on the FIRST `=` only.
        assert_eq!(
            parse_to_pairs("EXTRA=k=v"),
            vec![("EXTRA".into(), "k=v".into())]
        );
    }

    #[test]
    fn parse_info_preserves_insertion_order() {
        let parsed = parse_existing_info("ZZZ=1;AAA=2;MMM=3");
        let keys: Vec<&str> = parsed.keys().map(|k| k.as_str()).collect();
        assert_eq!(keys, vec!["ZZZ", "AAA", "MMM"]);
    }
