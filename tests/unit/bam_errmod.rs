
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
    let pl = em.pls(&mut b, &[0, 1], &mut Drand48::new(0));
    // hom-ref is the best; the het term is -4.343 * lhet[10, 0] = 30.1 phred.
    assert_eq!(pl[0], 0);
    assert_eq!(pl[1], 30);
    assert!(pl[2] > pl[1]);
}

#[test]
fn pl_pure_alt() {
    let em = ErrorModel::new();
    let mut b = bases(0, 20, 30);
    let pl = em.pls(&mut b, &[0, 1], &mut Drand48::new(0));
    assert_eq!(pl[2], 0);
    assert!(pl[1] > pl[2]);
    assert!(pl[0] > pl[1]);
}

#[test]
fn pl_balanced_het() {
    let em = ErrorModel::new();
    let mut b = bases(10, 10, 30);
    let pl = em.pls(&mut b, &[0, 1], &mut Drand48::new(0));
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
    let pl_same = em.pls(&mut same, &[0, 1], &mut Drand48::new(0));
    let pl_split = em.pls(&mut split, &[0, 1], &mut Drand48::new(0));
    assert!(pl_same[0] < pl_split[0], "{pl_same:?} vs {pl_split:?}");
}

#[test]
fn deep_sites_are_subsampled_deterministically() {
    let em = ErrorModel::new();
    let mut a = bases(300, 300, 30);
    let mut b = a.clone();
    let pa = em.pls(&mut a, &[0, 1], &mut Drand48::new(0));
    let pb = em.pls(&mut b, &[0, 1], &mut Drand48::new(0));
    assert_eq!(pa[1], 0);
    assert_eq!(pa.len(), pb.len());
}

#[test]
fn drand48_matches_hts_srand48() {
    // srand48(0); drand48() -> 0.170828...
    let mut r = Drand48::new(0);
    assert!((r.next() - 0.170828).abs() < 1e-6);
    // Same seed, same stream; a site key gives a stable stream too.
    let (mut a, mut b) = (
        Drand48::for_site(7, 1, 100, 0),
        Drand48::for_site(7, 1, 100, 0),
    );
    assert_eq!(a.next(), b.next());
    assert_ne!(
        Drand48::for_site(7, 1, 100, 0).next(),
        Drand48::for_site(8, 1, 100, 0).next()
    );
}

#[test]
fn deep_sites_subsample_deterministically() {
    let em = ErrorModel::new();
    let mut a = bases(300, 200, 30);
    let mut b = a.clone();
    let pa = em.pls(&mut a, &[0, 1], &mut Drand48::new(3));
    let pb = em.pls(&mut b, &[0, 1], &mut Drand48::new(3));
    assert_eq!(pa, pb);
    assert_eq!(a, b);
    // A different seed keeps a different subset.
    let mut c = bases(300, 200, 30);
    em.pls(&mut c, &[0, 1], &mut Drand48::new(4));
    assert_ne!(a[..MAX_BASES], c[..MAX_BASES]);
}
