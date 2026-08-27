    use super::*;

    #[test]
    fn parse_prior_freqs_af_tag() {
        let s = PriorFreqsSpec::Af("AF".into());
        let af = s.extract_af("DP=10;AF=0.3;AN=100", 2).unwrap();
        assert!((af[0] - 0.7).abs() < 1e-6);
        assert!((af[1] - 0.3).abs() < 1e-6);
    }

    #[test]
    fn parse_prior_freqs_an_ac() {
        let s = PriorFreqsSpec::AnAc { an: "AN".into(), ac: "AC".into() };
        let af = s.extract_af("DP=10;AN=100;AC=20", 2).unwrap();
        assert!((af[0] - 0.8).abs() < 1e-6);
        assert!((af[1] - 0.2).abs() < 1e-6);
    }

    #[test]
    fn ploidy_at_haploid_region() {
        let regions = vec![PloidyRegion {
            chrom: "X".into(), beg: 1, end: 1000, sex: "M".into(), ploidy: 1,
        }];
        let mut sex = HashMap::new();
        sex.insert("S1".into(), "M".to_string());
        sex.insert("S2".into(), "F".to_string());
        let samples = vec!["S1".to_string(), "S2".to_string()];
        let p = ploidy_at_site("X", 500, &samples, &sex, &regions, 2);
        assert_eq!(p, vec![1, 2]);
    }
