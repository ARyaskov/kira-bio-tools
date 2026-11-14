use super::*;

#[test]
fn discordant_genotypes_do_not_increment_hwe_match_count() {
    let mut check = Gtcheck {
        use_pls: 40,
        neg_log_e: 4.0,
        hom_only: false,
        calc_hwe: true,
        pairs: None,
        cross: false,
        acc: vec![PairAcc::default()],
        cnt: Counters::default(),
        ds: None,
        rng: Lrand48::new(0),
    };
    check.score(
        0,
        &(1, gt_to_prob(1, 4.0)),
        &(4, gt_to_prob(4, 4.0)),
        &[0.0; 8],
        &mut Vec::new(),
    );
    assert_eq!(check.acc[0].ncnt, 1);
    assert_eq!(check.acc[0].nmatch, 0);
    check.score(
        0,
        &(1, gt_to_prob(1, 4.0)),
        &(1, gt_to_prob(1, 4.0)),
        &[0.0; 8],
        &mut Vec::new(),
    );
    assert_eq!(check.acc[0].nmatch, 1);
}

#[test]
fn overlapping_reference_blocks_compare_with_variant_sites() {
    let side = || Side {
        samples: vec!["A".into()],
        sel: vec![0],
        use_gt: true,
        include: None,
        exclude: None,
    };
    let mut check = Gtcheck {
        use_pls: 0,
        neg_log_e: 4.0,
        hom_only: false,
        calc_hwe: true,
        pairs: None,
        cross: false,
        acc: vec![PairAcc::default()],
        cnt: Counters::default(),
        ds: None,
        rng: Lrand48::new(0),
    };
    let mut query = Rec::parse("1\t2\t.\tA\t<*>\t.\t.\tEND=5\tGT\t0/0".into());
    let mut panel = Rec::parse("1\t3\t.\tC\tT\t.\t.\t.\tGT\t0/0".into());
    assert!(compare_active(
        &mut check,
        &mut query,
        &mut panel,
        &side(),
        &side()
    ));
    panel = Rec::parse("1\t4\t.\tG\tA\t.\t.\t.\tGT\t1/1".into());
    assert!(compare_active(
        &mut check,
        &mut query,
        &mut panel,
        &side(),
        &side()
    ));
    assert_eq!(check.acc[0].ncnt, 2);
    assert_eq!(check.acc[0].ndiff, 1);
    panel = Rec::parse("1\t6\t.\tG\tA\t.\t.\t.\tGT\t0/0".into());
    assert!(!compare_active(
        &mut check,
        &mut query,
        &mut panel,
        &side(),
        &side()
    ));
    assert!(query.is_none());
}

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
