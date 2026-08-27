    use super::build_sample_map;

    #[test]
    fn test_build_sample_map_by_name() {
        let input = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let db = vec!["B".to_string(), "A".to_string()];
        let map = build_sample_map(&input, &db);
        assert_eq!(map, vec![Some(1), Some(0), None]);
    }
