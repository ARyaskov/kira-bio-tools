    use super::*;

    #[test]
    fn transition_matrices_are_column_stochastic() {
        let m = init_tprob_matrix(N_STATES, 1e-9, 0.5).unwrap();
        for j in 0..N_STATES {
            let s: f64 = (0..N_STATES).map(|i| m[N_STATES * i + j]).sum();
            assert!((s - 1.0).abs() < 1e-15);
        }
        assert!((m[0] - (1.0 - 3e-9)).abs() < 1e-18);
        let n = N_STATES * N_STATES;
        let m2 = init_tprob_matrix(n, 1e-9, 0.5).unwrap();
        for j in 0..n {
            let s: f64 = (0..n).map(|i| m2[n * i + j]).sum();
            assert!((s - 1.0).abs() < 1e-12, "col {j}: {s}");
        }
        // Same-state pairs are favoured over different-state pairs.
        assert!(m2[n * 10 + 10] > m2[n * 9 + 10]);
        assert!(init_tprob_matrix(N_STATES, 0.3, 0.5).is_err());
    }

    #[test]
    fn initial_probabilities() {
        let p = init_iprobs(N_STATES, 0.5);
        assert_eq!(p, vec![0.5 / 3.0, 0.5 / 3.0, 0.5, 0.5 / 3.0]);
        let p2 = init_iprobs(16, 0.5);
        assert!((p2.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        // The bcftools quirk: the first cn line prints these marginals.
        let cn2: f64 = (0..4).map(|j| p2[CN2 * 4 + j]).sum();
        assert!((cn2 - 0.5625).abs() < 1e-6, "{cn2}");
    }

    #[test]
    fn gaussian_helpers() {
        assert!((norm_cdf(0.0, 0.04) - 0.5).abs() < 1e-9);
        assert!((norm_cdf(0.5, 0.04) - 1.0).abs() < 1e-9);
        let mut s = Sample {
            baf_dev2: 0.04 * 0.04,
            cell_frac: 1.0,
            ..Default::default()
        };
        set_gauss_params(&mut s);
        assert!((s.gauss[CN3_PK_RRA].mean - 1.0 / 3.0).abs() < 1e-6);
        assert!((s.gauss[CN3_PK_RAA].mean - 2.0 / 3.0).abs() < 1e-6);
        // Density at the mean is the peak.
        let p = &s.gauss[CN2_PK_RA];
        assert!(norm_prob(0.5, p) > norm_prob(0.45, p));
        assert_eq!(copy_number_state(false, 3, 0), '3');
        assert_eq!(copy_number_state(true, 4 + 2, 0), '1');
        assert_eq!(copy_number_state(true, 4 + 2, 1), '2');
    }

    #[test]
    fn moving_average_smoothing() {
        let mut d = vec![0.0f32, 0.0, 3.0, 0.0, 0.0];
        smooth_data(&mut d, 3);
        assert_eq!(d, vec![0.0, 1.0, 1.0, 1.0, 0.0]);
        // Short arrays: the window is primed past the end with zeros.
        let mut short = vec![1.0f32, 1.0, 1.0, 1.0];
        smooth_data(&mut short, 10);
        assert_eq!(short, vec![0.8, 0.8, 0.8, 0.8]);
        let mut e = vec![1.0f32, 2.0, 3.0];
        smooth_data(&mut e, 1);
        assert_eq!(e, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn option_parsing() {
        let a: Vec<String> = "-s Q -c C -d 0.05,0.06 -k 0.3 -a 0.7,0.8 -O 0.8 -o out in.vcf"
            .split(' ')
            .map(String::from)
            .collect();
        let c = parse_args(&a).unwrap();
        assert_eq!(c.query_sample.as_deref(), Some("Q"));
        assert_eq!(c.control_sample.as_deref(), Some("C"));
        assert!((c.baf_dev2.0 - 0.0025).abs() < 1e-6 && (c.baf_dev2.1 - 0.0036).abs() < 1e-6);
        assert!((c.lrr_dev2.1 - 0.09).abs() < 1e-6);
        assert_eq!(c.cell_frac, (0.7, 0.8));
        assert_eq!(c.optimize_frac, 0.8);
        assert_eq!(c.input, PathBuf::from("in.vcf"));
        let a: Vec<String> = vec!["in.vcf".into()];
        assert!(
            parse_args(&a)
                .unwrap_err()
                .to_string()
                .contains("Expected -o option")
        );
    }
