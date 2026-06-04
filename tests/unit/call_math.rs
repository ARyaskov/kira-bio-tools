    use super::*;

    #[test]
    fn pl_to_prob_table_matches_formula() {
        let t = init_pl2p();
        assert!((pl_to_prob(0, &t) - 1.0).abs() < 1e-12);
        assert!((pl_to_prob(10, &t) - 0.1).abs() < 1e-12);
        assert!((pl_to_prob(20, &t) - 0.01).abs() < 1e-12);
        assert!((pl_to_prob(30, &t) - 0.001).abs() < 1e-12);
    }

    #[test]
    fn log10_sumexp_basic() {
        let r = log10_sum_exp((0.5f64).log10(), (0.5f64).log10());
        assert!((r - 0f64).abs() < 1e-12, "got {r}");
    }

    #[test]
    fn gt_index_canonical() {
        assert_eq!(gt_index(0, 0), 0);
        assert_eq!(gt_index(0, 1), 1);
        assert_eq!(gt_index(1, 1), 2);
        assert_eq!(gt_index(0, 2), 3);
        assert_eq!(gt_index(1, 2), 4);
        assert_eq!(gt_index(2, 2), 5);
        assert_eq!(gt_index(2, 1), 4);
    }

    #[test]
    fn n_genotypes_correct() {
        assert_eq!(n_genotypes(1), 1);
        assert_eq!(n_genotypes(2), 3);
        assert_eq!(n_genotypes(3), 6);
        assert_eq!(n_genotypes(4), 10);
    }

    #[test]
    fn watterson_2n_alleles() {
        let w2 = watterson_factor(2);
        assert!((w2 - 1.0).abs() < 1e-12);
        let w4 = watterson_factor(4);
        assert!((w4 - (1.0 + 0.5 + 1.0/3.0)).abs() < 1e-12);
    }
