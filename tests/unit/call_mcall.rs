    use super::*;

    fn site_ref_only(n_smpl: usize) -> CallSite {
        let n_gt = 3;
        let mut pls = vec![0i32; n_smpl * n_gt];
        for s in 0..n_smpl {
            pls[s * n_gt + 0] = 0;
            pls[s * n_gt + 1] = 60;
            pls[s * n_gt + 2] = 100;
        }
        CallSite::new(n_smpl, 2, pls)
    }

    fn site_homozygous_alt(n_smpl: usize) -> CallSite {
        let n_gt = 3;
        let mut pls = vec![0i32; n_smpl * n_gt];
        for s in 0..n_smpl {
            pls[s * n_gt + 0] = 100;
            pls[s * n_gt + 1] = 60;
            pls[s * n_gt + 2] = 0;
        }
        CallSite::new(n_smpl, 2, pls)
    }

    fn site_het(n_smpl: usize) -> CallSite {
        let n_gt = 3;
        let mut pls = vec![0i32; n_smpl * n_gt];
        for s in 0..n_smpl {
            pls[s * n_gt + 0] = 30;
            pls[s * n_gt + 1] = 0;
            pls[s * n_gt + 2] = 30;
        }
        CallSite::new(n_smpl, 2, pls)
    }

    #[test]
    fn pure_ref_no_variant() {
        let caller = Caller::new(CallerOpts::default(), 3);
        let mut s = site_ref_only(3);
        let r = caller.call_site(&mut s);
        match r {
            CallResult::Called { alleles_kept, gts, .. } => {
                assert_eq!(alleles_kept, vec![0]);
                for &(a, b) in &gts { assert_eq!((a, b), (0, 0)); }
            }
            _ => panic!(),
        }
    }

    #[test]
    fn pure_alt_homozygous() {
        let caller = Caller::new(CallerOpts::default(), 3);
        let mut s = site_homozygous_alt(3);
        let r = caller.call_site(&mut s);
        match r {
            CallResult::Called { alleles_kept, gts, qual, .. } => {
                assert_eq!(alleles_kept, vec![0, 1]);
                for &(a, b) in &gts { assert_eq!((a, b), (1, 1)); }
                assert!(qual > 30.0, "expected high QUAL for clear variant, got {qual}");
            }
            _ => panic!(),
        }
    }

    #[test]
    fn het_calls_correctly() {
        let caller = Caller::new(CallerOpts::default(), 5);
        let mut s = site_het(5);
        let r = caller.call_site(&mut s);
        match r {
            CallResult::Called { alleles_kept, gts, .. } => {
                assert_eq!(alleles_kept, vec![0, 1]);
                for &(a, b) in &gts {
                    assert!((a, b) == (0, 1) || (a, b) == (1, 0), "expected het, got {a}/{b}");
                }
            }
            _ => panic!(),
        }
    }

    #[test]
    fn variants_only_skips_pure_ref() {
        let opts = CallerOpts { variants_only: true, ..Default::default() };
        let caller = Caller::new(opts, 3);
        let mut s = site_ref_only(3);
        assert!(matches!(caller.call_site(&mut s), CallResult::Skip));
    }

    #[test]
    fn indel_theta_used_for_indel_site() {
        let opts_a = CallerOpts { theta: 0.5, indel_theta: 0.5, ..Default::default() };
        let opts_b = CallerOpts { theta: 0.5, indel_theta: 1e-30, ..Default::default() };
        let mk = || {
            let mut s = site_het(3);
            s.is_indel = true;
            s
        };
        let r_a = Caller::new(opts_a, 3).call_site(&mut mk());
        let r_b = Caller::new(opts_b, 3).call_site(&mut mk());
        let n_a = match r_a { CallResult::Called { alleles_kept, .. } => alleles_kept.len(), _ => 0 };
        let n_b = match r_b { CallResult::Called { alleles_kept, .. } => alleles_kept.len(), _ => 0 };
        assert!(n_a >= n_b, "high indel_theta should keep more alts, n_a={n_a} n_b={n_b}");
    }

    #[test]
    fn keep_alts_preserves_input_alts() {
        let mut s = site_ref_only(3);
        s.n_alleles = 3;
        let n_gt = 6;
        let mut pls = vec![0i32; 3 * n_gt];
        for i in 0..3 {
            pls[i * n_gt + 0] = 0;
            pls[i * n_gt + 1] = 60;
            pls[i * n_gt + 2] = 100;
            pls[i * n_gt + 3] = 60;
            pls[i * n_gt + 4] = 100;
            pls[i * n_gt + 5] = 200;
        }
        s.pls = pls;
        let opts = CallerOpts { keep_alts: true, ..Default::default() };
        let caller = Caller::new(opts, 3);
        let r = caller.call_site(&mut s);
        match r {
            CallResult::Called { alleles_kept, .. } => {
                assert_eq!(alleles_kept, vec![0, 1, 2], "keep_alts must preserve all input alts");
            }
            _ => panic!(),
        }
    }

    #[test]
    fn per_sample_ploidy_haploid_sample() {
        let opts = CallerOpts {
            ploidy: 2,
            per_sample_ploidy: Some(vec![2, 1, 2]),
            ..Default::default()
        };
        let caller = Caller::new(opts, 3);
        let mut s = site_homozygous_alt(3);
        let r = caller.call_site(&mut s);
        match r {
            CallResult::Called { gts, .. } => {
                assert_eq!(gts.len(), 3);
                assert_eq!(gts[0], (1, 1));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn prior_af_used_when_set() {
        let opts = CallerOpts { prior_af: Some(vec![0.9, 0.1]), ..Default::default() };
        let caller = Caller::new(opts, 3);
        let mut s = site_het(3);
        let r = caller.call_site(&mut s);
        match r {
            CallResult::Called { alleles_kept, .. } => assert!(!alleles_kept.is_empty()),
            _ => panic!(),
        }
    }

    #[test]
    fn grouped_calling_independent() {
        let groups = vec![
            SampleGroup { name: "g1".into(), sample_idxs: vec![0, 1] },
            SampleGroup { name: "g2".into(), sample_idxs: vec![2, 3] },
        ];
        let opts = CallerOpts { groups: Some(groups), ..Default::default() };
        let caller = Caller::new(opts, 4);
        let n_gt = 3;
        let mut pls = vec![0i32; 4 * n_gt];
        for s in 0..4 {
            pls[s * n_gt + 0] = if s < 2 { 0 } else { 255 };
            pls[s * n_gt + 1] = 60;
            pls[s * n_gt + 2] = if s < 2 { 255 } else { 0 };
        }
        let mut site = CallSite::new(4, 2, pls);
        let r = caller.call_site(&mut site);
        match r {
            CallResult::Called { gts, .. } => {
                assert_eq!(gts[0], (0, 0));
                assert_eq!(gts[1], (0, 0));
                assert_eq!(gts[2], (1, 1));
                assert_eq!(gts[3], (1, 1));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn trio_constraint_fixes_mendelian_violation() {
        let f = (0, 0);
        let m = (0, 0);
        let c_violation = (1, 1);
        assert!(!mendelian_ok(f, m, c_violation));
        let mut gts = vec![f, m, c_violation];
        let mut gqs = vec![50u32, 50, 50];
        let fam = TrioFamily { father: Some(0), mother: Some(1), child: Some(2), is_son: false };
        apply_trio_constraint(&mut gts, &mut gqs, &fam, &[0, 1, 2]);
        assert_eq!(gts[2], (0, 0));
        assert!(gqs[2] < 50);
    }

    #[test]
    fn mendelian_consistent_no_change() {
        let f = (0, 1);
        let m = (0, 0);
        let c = (0, 1);
        assert!(mendelian_ok(f, m, c));
        let mut gts = vec![f, m, c];
        let mut gqs = vec![50u32, 50, 50];
        let fam = TrioFamily { father: Some(0), mother: Some(1), child: Some(2), is_son: false };
        apply_trio_constraint(&mut gts, &mut gqs, &fam, &[0, 1, 2]);
        assert_eq!(gts[2], (0, 1));
        assert_eq!(gqs[2], 50);
    }
