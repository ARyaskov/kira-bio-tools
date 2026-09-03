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

    #[test]
    fn ploidy_file_accepts_wildcards() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ploidy.txt");
        std::fs::write(&p, "X\t1\t1000\tM\t1\n*\t*\t*\tM\t2\n*\t*\t*\tF\t2\n").unwrap();
        let regions = parse_ploidy_file(&p).unwrap();
        assert_eq!(regions.len(), 3);
        assert_eq!(regions[1].beg, 0);
        assert_eq!(regions[1].end, u32::MAX);

        let mut sex = HashMap::new();
        sex.insert("S1".to_string(), "M".to_string());
        let samples = vec!["S1".to_string()];
        // Inside the named region the specific line wins.
        assert_eq!(ploidy_at_site("X", 500, &samples, &sex, &regions, 9), vec![1]);
        // Outside it the catch-all applies, not the caller's default.
        assert_eq!(ploidy_at_site("X", 5000, &samples, &sex, &regions, 9), vec![2]);
        assert_eq!(ploidy_at_site("1", 5000, &samples, &sex, &regions, 9), vec![2]);
    }

    #[test]
    fn sex_file_reads_both_shapes() {
        let dir = tempfile::tempdir().unwrap();
        let two_col = dir.path().join("samples.txt");
        std::fs::write(&two_col, "HG00100\tF\nHG00101\tM\n").unwrap();
        let m = parse_sex_file(&two_col).unwrap();
        assert_eq!(m.get("HG00101").map(String::as_str), Some("M"));

        let ped = dir.path().join("fam.ped");
        std::fs::write(&ped, "Fam1\tsmpl1\t0\t0\t1\nFam1\tsmpl2\t0\t0\t2\n").unwrap();
        // PED sex codes become the M/F labels a ploidy file is keyed on.
        let m = parse_sex_file(&ped).unwrap();
        assert_eq!(m.get("smpl1").map(String::as_str), Some("M"));
        assert_eq!(m.get("smpl2").map(String::as_str), Some("F"));
    }
