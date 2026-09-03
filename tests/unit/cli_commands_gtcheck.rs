    use super::*;

    #[test]
    fn scientific_format_matches_c() {
        assert_eq!(sci6(0.0), "0.000000e+00");
        assert_eq!(sci6(0.1234567891), "1.234568e-01");
        assert_eq!(sci6(12345.678), "1.234568e+04");
    }

    #[test]
    fn dosage_masks() {
        assert_eq!(gt_to_dsg(Some([0, 0])), 1);
        assert_eq!(gt_to_dsg(Some([0, 1])), 2);
        assert_eq!(gt_to_dsg(Some([1, 1])), 4);
        assert_eq!(gt_to_dsg(None), 0);
        assert_eq!(pl_to_dsg(Some([0, 10, 20])), 1);
        assert_eq!(pl_to_dsg(Some([0, 0, 20])), 3);
        assert_eq!(pl_to_dsg(None), 0);
    }

    #[test]
    fn error_model_scores_like_bcftools() {
        // -E 10: e = 0.1, -ln e = 2.302585; 0/0 vs 1/1 costs two errors.
        let e = -(10f64.powf(-1.0)).ln();
        let a = gt_to_prob(1, e);
        let d = gt_to_prob(2, e);
        let h = gt_to_prob(4, e);
        let min = |x: [f64; 3], y: [f64; 3]| (0..3).map(|k| x[k] + y[k]).fold(f64::INFINITY, f64::min);
        let ln10 = std::f64::consts::LN_10;
        assert!((min(a, d) - ln10).abs() < 1e-9);
        assert!((min(a, h) - 2.0 * ln10).abs() < 1e-9);
    }

    #[test]
    fn lrand48_matches_the_c_library() {
        // First draw of lrand48() after srand48(0): ((a*0x330E + c) mod 2^48) >> 17.
        let mut r = Lrand48::new(0);
        assert_eq!(r.next(), 366850414);
        assert_ne!(r.next(), 366850414);
    }

    #[test]
    fn record_genotypes_require_diploid_max_ploidy() {
        let r = Rec::parse("1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1".to_string()).unwrap();
        let g = r.genotypes().unwrap().unwrap();
        assert_eq!(g[0], Some([0, 1]));
        assert_eq!(g[1], None);
        let r = Rec::parse("1\t10\t.\tA\tC\t.\t.\t.\tGT\t0\t1".to_string()).unwrap();
        assert!(r.genotypes().unwrap().is_none());
        let r = Rec::parse("1\t10\t.\tA\tC\t.\t.\tAN=4;AC=1\tGT\t0/1\t0/0".to_string()).unwrap();
        assert_eq!(r.allele_counts(), Some((3, 1)));
    }
