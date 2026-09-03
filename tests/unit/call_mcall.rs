    use super::*;

    fn site_ref_only(n_smpl: usize) -> CallSite {
        let n_gt = 3;
        let mut pls = vec![0i32; n_smpl * n_gt];
        for s in 0..n_smpl {
            pls[s * n_gt] = 0;
            pls[s * n_gt + 1] = 60;
            pls[s * n_gt + 2] = 100;
        }
        CallSite::new(n_smpl, 2, pls)
    }

    fn site_homozygous_alt(n_smpl: usize) -> CallSite {
        let n_gt = 3;
        let mut pls = vec![0i32; n_smpl * n_gt];
        for s in 0..n_smpl {
            pls[s * n_gt] = 100;
            pls[s * n_gt + 1] = 60;
            pls[s * n_gt + 2] = 0;
        }
        CallSite::new(n_smpl, 2, pls)
    }

    fn site_het(n_smpl: usize) -> CallSite {
        let n_gt = 3;
        let mut pls = vec![0i32; n_smpl * n_gt];
        for s in 0..n_smpl {
            pls[s * n_gt] = 30;
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
            CallResult::Called { alleles_kept, gts, qual, pls, an, .. } => {
                assert_eq!(alleles_kept, vec![0]);
                for gt in &gts {
                    assert_eq!(*gt, Some((0, 0)));
                }
                assert_eq!(an, 6);
                // Ref site: QUAL is phred(P(variant)), small but positive; PL is dropped.
                let q = qual.expect("ref-site QUAL");
                assert!(q > 0.0 && q < 40.0, "ref-site QUAL {q}");
                assert!(pls[0].is_empty());
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
            CallResult::Called { alleles_kept, gts, qual, gqs, pls, ac, an, .. } => {
                assert_eq!(alleles_kept, vec![0, 1]);
                for gt in &gts {
                    assert_eq!(*gt, Some((1, 1)));
                }
                let q = qual.unwrap();
                assert!(q > 30.0, "expected high QUAL for clear variant, got {q}");
                assert!(gqs.iter().all(|&g| g >= 50), "{gqs:?}");
                // The input PLs are kept verbatim for the retained alleles.
                assert_eq!(pls[0], vec![100, 60, 0]);
                assert_eq!((ac, an), (vec![0, 6], 6));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn qual_matches_bcftools_formula() {
        // Two samples, QS given: QUAL = -4.343*(ref_lk - logsumexp(lk_sum, ref_lk)).
        // ref_lk = ln(6.66e-7) + ln(0.666) = -14.63; the 0/1 hypothesis at
        // qsum 0.5/0.5 scores -2.196 + ln(theta*a_m) = -8.40, so QUAL = 27.06.
        let pls = vec![60, 3, 0, 0, 3, 60];
        let mut s = CallSite::new(2, 2, pls);
        s.qs = Some(vec![1.0, 1.0]);
        let r = Caller::new(CallerOpts::default(), 2).call_site(&mut s);
        let CallResult::Called { qual, gts, gps, .. } = r else { panic!() };
        // Equal frequencies give the 1:2:1 HWE prior: 0/1 (2*0.501*0.25) just
        // beats 1/1 (1*0.25).
        assert_eq!(gts[0], Some((0, 1)));
        assert!((gps[0][1] - 0.5006).abs() < 1e-3, "{:?}", gps[0]);
        assert!((gps[0][2] - 0.4994).abs() < 1e-3, "{:?}", gps[0]);
        let q = qual.unwrap();
        assert!((q - 27.06).abs() < 0.2, "QUAL {q}");
    }

    #[test]
    fn het_calls_correctly() {
        let caller = Caller::new(CallerOpts::default(), 5);
        let mut s = site_het(5);
        let r = caller.call_site(&mut s);
        match r {
            CallResult::Called { alleles_kept, gts, gps, .. } => {
                assert_eq!(alleles_kept, vec![0, 1]);
                for gt in &gts {
                    assert!(*gt == Some((0, 1)) || *gt == Some((1, 0)), "expected het, got {gt:?}");
                }
                let gp = &gps[0];
                assert_eq!(gp.len(), 3);
                assert!((gp.iter().sum::<f64>() - 1.0).abs() < 1e-9);
                assert!(gp[1] > gp[0] && gp[1] > gp[2]);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn gq_is_phred_of_one_minus_best_posterior() {
        // Balanced weights between two genotypes give GQ = 3 (phred 0.5).
        let (gq, gp) = gq_gp(&[0.0, 1.0, 1.0]);
        assert_eq!(gq, 3);
        assert!((gp[1] - 0.5).abs() < 1e-12);
        let (gq, _) = gq_gp(&[1.0, 1e-30, 0.0]);
        assert_eq!(gq, 127);
    }

    #[test]
    fn variants_only_skips_pure_ref() {
        let opts = CallerOpts { variants_only: true, ..Default::default() };
        let caller = Caller::new(opts, 3);
        let mut s = site_ref_only(3);
        assert!(matches!(caller.call_site(&mut s), CallResult::Skip));
    }

    #[test]
    fn missing_samples_get_missing_genotypes() {
        let mut pls = vec![100, 60, 0];
        pls.extend([PL_MISSING; 3]);
        let mut s = CallSite::new(2, 2, pls);
        let r = Caller::new(CallerOpts::default(), 2).call_site(&mut s);
        let CallResult::Called { gts, an, gqs, .. } = r else { panic!() };
        assert_eq!(gts[0], Some((1, 1)));
        assert_eq!(gts[1], None);
        assert_eq!(an, 2);
        assert_eq!(gqs[1], 0);
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
            pls[i * n_gt] = 0;
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
            CallResult::Called { alleles_kept, pls, .. } => {
                assert_eq!(alleles_kept, vec![0, 1, 2], "keep_alts must preserve all input alts");
                assert_eq!(pls[0].len(), 6);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn unseen_allele_is_never_emitted() {
        // Alleles: REF, ALT, <*>; the caller must not keep the unseen one.
        let n_gt = 6;
        let mut pls = Vec::new();
        for _ in 0..3 {
            pls.extend([100, 60, 0, 100, 60, 100]);
        }
        let mut s = CallSite::new(3, 3, pls);
        s.qs = Some(vec![0.0, 27.0, 0.0]);
        s.unseen = Some(2);
        let r = Caller::new(CallerOpts::default(), 3).call_site(&mut s);
        let CallResult::Called { alleles_kept, gts, pls, .. } = r else { panic!() };
        assert_eq!(alleles_kept, vec![0, 1]);
        assert_eq!(gts[0], Some((1, 1)));
        assert_eq!(pls[0], vec![100, 60, 0]);
        assert_eq!(pls[0].len(), n_gt / 2);
    }

    #[test]
    fn per_sample_ploidy_haploid_sample() {
        let opts = CallerOpts { ploidy: 2, per_sample_ploidy: Some(vec![2, 1, 2]), ..Default::default() };
        let caller = Caller::new(opts, 3);
        let mut s = site_homozygous_alt(3);
        let r = caller.call_site(&mut s);
        match r {
            CallResult::Called { gts, pls, an, .. } => {
                assert_eq!(gts.len(), 3);
                assert_eq!(gts[0], Some((1, 1)));
                assert_eq!(gts[1], Some((1, 1)));
                assert_eq!(pls[1], vec![100, 0], "haploid PL has one value per allele");
                assert_eq!(an, 5);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn panel_counts_shrink_the_frequency_estimate() {
        // A panel that only ever saw the ALT pulls the frequency estimate
        // towards it: the 1/1 posterior rises, though PL 30 keeps the het call.
        let gp11 = |prior: Option<(f64, Vec<f64>)>| {
            let mut s = site_het(3);
            s.qs = Some(vec![10.0, 10.0]);
            s.prior_an_ac = prior;
            let r = Caller::new(CallerOpts::default(), 3).call_site(&mut s);
            let CallResult::Called { alleles_kept, gts, gps, .. } = r else { panic!() };
            assert_eq!(alleles_kept, vec![0, 1]);
            assert!(gts.iter().all(|g| *g == Some((0, 1))), "{gts:?}");
            gps[0][2]
        };
        let flat = gp11(None);
        let shrunk = gp11(Some((1000.0, vec![1000.0])));
        assert!(shrunk > 10.0 * flat, "flat {flat} shrunk {shrunk}");
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
            pls[s * n_gt] = if s < 2 { 0 } else { 255 };
            pls[s * n_gt + 1] = 60;
            pls[s * n_gt + 2] = if s < 2 { 255 } else { 0 };
        }
        let mut site = CallSite::new(4, 2, pls);
        // Group frequencies from per-sample AD fractions.
        site.sample_af = Some(vec![vec![1.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0], vec![0.0, 1.0]]);
        let r = caller.call_site(&mut site);
        match r {
            CallResult::Called { gts, alleles_kept, .. } => {
                assert_eq!(alleles_kept, vec![0, 1]);
                assert_eq!(gts[0], Some((0, 0)));
                assert_eq!(gts[1], Some((0, 0)));
                assert_eq!(gts[2], Some((1, 1)));
                assert_eq!(gts[3], Some((1, 1)));
            }
            _ => panic!(),
        }
    }

    fn trio_opts() -> CallerOpts {
        CallerOpts {
            constrain: ConstrainMode::Trio,
            families: vec![TrioFamily { father: Some(0), mother: Some(1), child: Some(2), is_son: false }],
            ..Default::default()
        }
    }

    #[test]
    fn trio_joint_likelihood_overrides_weak_de_novo() {
        // Parents clearly 0/0; the child's data favour 1/1 by 60 phred, enough
        // for the site to be called variant but not to beat the 1e-8 de-novo
        // rate, so the joint model keeps the child at 0/0.
        let pls = vec![0, 80, 120, 0, 80, 120, 60, 30, 0];
        let mut s = CallSite::new(3, 2, pls);
        let r = Caller::new(trio_opts(), 3).call_site(&mut s);
        let CallResult::Called { alleles_kept, gts, gqs, .. } = r else { panic!() };
        assert_eq!(alleles_kept, vec![0, 1]);
        assert_eq!(gts[0], Some((0, 0)));
        assert_eq!(gts[1], Some((0, 0)));
        assert_eq!(gts[2], Some((0, 0)), "child should follow Mendelian inheritance");
        assert!(gqs[2] > 0 && gqs[2] < 127, "GQ {}", gqs[2]);
    }

    #[test]
    fn trio_keeps_consistent_genotypes() {
        let pls = vec![40, 0, 40, 0, 80, 120, 40, 0, 40];
        let mut s = CallSite::new(3, 2, pls);
        let r = Caller::new(trio_opts(), 3).call_site(&mut s);
        let CallResult::Called { gts, .. } = r else { panic!() };
        assert_eq!(gts[0], Some((0, 1)));
        assert_eq!(gts[1], Some((0, 0)));
        assert_eq!(gts[2], Some((0, 1)));
        assert_eq!(transmissions((0, 1), (0, 0), (0, 1)), 2);
        assert_eq!(transmissions((0, 0), (0, 0), (1, 1)), 0);
    }

    #[test]
    fn strong_de_novo_survives_the_prior() {
        // Overwhelming child evidence beats the 1e-8 de-novo rate (80 phred > 1e-8).
        let pls = vec![0, 80, 120, 0, 80, 120, 120, 100, 0];
        let mut s = CallSite::new(3, 2, pls);
        let r = Caller::new(trio_opts(), 3).call_site(&mut s);
        let CallResult::Called { gts, .. } = r else { panic!() };
        assert_eq!(gts[2], Some((1, 1)));
    }
