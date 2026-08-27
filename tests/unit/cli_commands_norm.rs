    use super::join_adjacent_snps;

    #[test]
    fn joins_three_consecutive_snps() {
        let mut buf = vec![
            "1\t100\t.\tA\tT\t.\t.\t.".to_string(),
            "1\t101\t.\tC\tG\t.\t.\t.".to_string(),
            "1\t102\t.\tG\tA\t.\t.\t.".to_string(),
        ];
        let out = join_adjacent_snps(&mut buf);
        assert_eq!(out.len(), 1);
        let cols: Vec<&str> = out[0].split('\t').collect();
        assert_eq!(cols[3], "ACG");
        assert_eq!(cols[4], "TGA");
    }

    #[test]
    fn skips_non_snp() {
        let mut buf = vec![
            "1\t100\t.\tA\tT\t.\t.\t.".to_string(),
            "1\t101\t.\tCG\tC\t.\t.\t.".to_string(),
            "1\t102\t.\tG\tA\t.\t.\t.".to_string(),
        ];
        let out = join_adjacent_snps(&mut buf);
        assert_eq!(out.len(), 3);
    }

    #[test]
    fn does_not_join_non_adjacent() {
        let mut buf = vec![
            "1\t100\t.\tA\tT\t.\t.\t.".to_string(),
            "1\t105\t.\tC\tG\t.\t.\t.".to_string(),
        ];
        let out = join_adjacent_snps(&mut buf);
        assert_eq!(out.len(), 2);
    }
