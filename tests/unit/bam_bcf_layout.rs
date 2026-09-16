use super::*;

#[test]
fn erfc_matches_libm() {
    assert!((kf_erfc(0.0) - 1.0).abs() < 1e-12);
    assert!((kf_erfc(1.0) - 0.15729920705).abs() < 1e-9);
    assert!((kf_erfc(-1.0) - 1.84270079295).abs() < 1e-9);
    assert_eq!(kf_erfc(40.0), 0.0);
    assert_eq!(kf_erfc(-40.0), 2.0);
}

#[test]
fn fisher_two_tails() {
    let (_, _, two) = kt_fisher_exact(3, 1, 1, 3);
    assert!((two - 0.4857142857).abs() < 1e-6, "{two}");
    let (_, _, two) = kt_fisher_exact(5, 5, 5, 5);
    assert!((two - 1.0).abs() < 1e-9);
    let (_, _, two) = kt_fisher_exact(20, 0, 0, 20);
    assert!(two < 1e-9);
    assert_eq!(strand_bias_sp(&[1, 1, 1, 1]), 0);
    assert!(strand_bias_sp(&[20, 0, 0, 20]) > 60);
}

#[test]
fn mwu_z_is_zero_for_identical_histograms() {
    let a = [0, 3, 5, 2, 0];
    assert!(mwu_z(&a, &a).unwrap().abs() < 1e-12);
    assert!(mwu_z(&a, &[0; 5]).is_none());
    // U counts a-below-b pairs: all of b above a is strongly positive.
    assert!(mwu_z(&[5, 0, 0, 0], &[0, 0, 0, 5]).unwrap() > 2.0);
    assert!(mwu_z(&[0, 0, 0, 5], &[5, 0, 0, 0]).unwrap() < -2.0);
}

#[test]
fn vdb_two_reads_at_one_position() {
    let mut pos = [0i32; NPOS];
    pos[40] = 2;
    // dp == 2, mean_diff 0: (2*100-2-1)*1/99 = 1 (integer division), /50.
    assert!((calc_vdb(&pos).unwrap() - 0.02).abs() < 1e-12);
    pos[40] = 1;
    assert!(calc_vdb(&pos).is_none());
}

fn call_with_qs(qs: [i64; 4]) -> SampleCall {
    let mut c = new_call();
    c.qs[..4].copy_from_slice(&qs);
    c.p = vec![0.0; 25];
    c
}

#[test]
fn allele_order_follows_bcftools() {
    // Only the reference seen: the unseen allele is the last zero-qsum base.
    let comb = combine(&[call_with_qs([10, 0, 0, 0])], Some(0));
    assert_eq!(comb.n_alleles, 2);
    assert_eq!(comb.unseen, 1);
    assert_eq!(comb.a[1], 3);
    assert!(!comb.has_alt);
    // Ref C with G seen: alleles C, G, then <*> standing in for T.
    let comb = combine(&[call_with_qs([0, 9, 1, 0])], Some(1));
    assert_eq!(comb.n_alleles, 3);
    assert_eq!(&comb.a[..3], &[1, 2, 3]);
    assert_eq!(comb.unseen, 2);
    assert!(comb.has_alt);
    assert!((comb.qsum[0] - 0.9).abs() < 1e-6 && (comb.qsum[1] - 0.1).abs() < 1e-6);
    // Reference N: N first, then the seen base, then <*>.
    let comb = combine(&[call_with_qs([5, 0, 0, 0])], Some(4));
    assert_eq!(comb.n_alleles, 3);
    assert_eq!(&comb.a[..3], &[4, 0, 3]);
    assert_eq!(comb.unseen, 2);
    // Indel site: no unseen allele, one alt.
    let comb = combine(&[call_with_qs([5, 2, 0, 0])], None);
    assert_eq!(comb.n_alleles, 2);
    assert_eq!(comb.unseen, -1);
    assert!(comb.has_alt);
}

#[test]
fn default_header_carries_bias_tags() {
    let tags = AnnotateSpec::bcftools_defaults();
    let lines = header_lines(&tags);
    assert!(lines.iter().any(|l| l.starts_with("##INFO=<ID=VDB,")));
    assert!(lines.iter().any(|l| l.starts_with("##INFO=<ID=MQ0F,")));
    assert!(!lines.iter().any(|l| l.starts_with("##FORMAT=<ID=DP,")));
    let spec = AnnotateSpec::parse_bcftools(Some("DP,DV,-AD,INFO/SCR,FMT/SCR")).unwrap();
    assert!(spec.fmt_dp && spec.fmt_dv && !spec.fmt_ad && spec.info_scr && spec.fmt_scr);
    assert!(AnnotateSpec::parse_bcftools(Some("VDB")).is_err());
}
