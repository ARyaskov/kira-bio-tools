    use super::*;

    #[test]
    fn pl_pure_ref() {
        let em = ErrorModel::new();
        let gl = em.likelihoods(2, &[20, 0], &[600, 0]);
        assert_eq!(gl.pl[0], 0);
        assert!(gl.pl[1] > gl.pl[0]);
        assert!(gl.pl[2] > gl.pl[1]);
    }

    #[test]
    fn pl_pure_alt() {
        let em = ErrorModel::new();
        let gl = em.likelihoods(2, &[0, 20], &[0, 600]);
        assert_eq!(gl.pl[2], 0);
        assert!(gl.pl[1] > gl.pl[2]);
        assert!(gl.pl[0] > gl.pl[1]);
    }

    #[test]
    fn pl_balanced_het() {
        let em = ErrorModel::new();
        let gl = em.likelihoods(2, &[10, 10], &[300, 300]);
        assert_eq!(gl.pl[1], 0);
        assert!(gl.pl[0] > gl.pl[1]);
        assert!(gl.pl[2] > gl.pl[1]);
    }
