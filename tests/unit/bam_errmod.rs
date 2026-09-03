    use super::*;

    fn bases(n_ref: usize, n_alt: usize, q: u8) -> Vec<u16> {
        let mut v = Vec::new();
        for i in 0..n_ref {
            v.push(pack_base(q, i % 2 == 1, 0));
        }
        for i in 0..n_alt {
            v.push(pack_base(q, i % 2 == 1, 1));
        }
        v
    }

    #[test]
    fn pl_pure_ref_matches_htslib_het_term() {
        let em = ErrorModel::new();
        let mut b = bases(10, 0, 30);
        let pl = em.pls(&mut b, &[0, 1]);
        // hom-ref is the best; the het term is -4.343 * lhet[10, 0] = 30.1 phred.
        assert_eq!(pl[0], 0);
        assert_eq!(pl[1], 30);
        assert!(pl[2] > pl[1]);
    }

    #[test]
    fn pl_pure_alt() {
        let em = ErrorModel::new();
        let mut b = bases(0, 20, 30);
        let pl = em.pls(&mut b, &[0, 1]);
        assert_eq!(pl[2], 0);
        assert!(pl[1] > pl[2]);
        assert!(pl[0] > pl[1]);
    }

    #[test]
    fn pl_balanced_het() {
        let em = ErrorModel::new();
        let mut b = bases(10, 10, 30);
        let pl = em.pls(&mut b, &[0, 1]);
        assert_eq!(pl[1], 0);
        assert!(pl[0] > pl[1]);
        assert!(pl[2] > pl[1]);
    }

    #[test]
    fn dependent_errors_weigh_repeats_less() {
        // Ten alt bases all on one strand carry less evidence than five per strand.
        let em = ErrorModel::new();
        let mut same: Vec<u16> = (0..10).map(|_| pack_base(30, false, 1)).collect();
        same.extend((0..10).map(|_| pack_base(30, false, 0)));
        let mut split = bases(10, 10, 30);
        let pl_same = em.pls(&mut same, &[0, 1]);
        let pl_split = em.pls(&mut split, &[0, 1]);
        assert!(pl_same[0] < pl_split[0], "{pl_same:?} vs {pl_split:?}");
    }

    #[test]
    fn deep_sites_are_subsampled_deterministically() {
        let em = ErrorModel::new();
        let mut a = bases(300, 300, 30);
        let mut b = a.clone();
        let pa = em.pls(&mut a, &[0, 1]);
        let pb = em.pls(&mut b, &[0, 1]);
        assert_eq!(pa[1], 0);
        assert_eq!(pa.len(), pb.len());
    }
