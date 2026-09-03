    use super::*;

    #[test]
    fn fisher_two_tailed_matches_known_values() {
        // Balanced table: no bias.
        assert!((fisher_exact_two_tailed(5, 5, 5, 5) - 1.0).abs() < 1e-9);
        // Classic tea-tasting table: p = 0.4857.
        let p = fisher_exact_two_tailed(3, 1, 1, 3);
        assert!((p - 0.4857142857).abs() < 1e-6, "{p}");
        // Strongly biased.
        let p = fisher_exact_two_tailed(20, 0, 0, 20);
        assert!(p < 1e-9);
        assert_eq!(strand_bias_phred(10, 10, 5, 5), 0);
        assert!(strand_bias_phred(20, 0, 0, 20) > 60);
        assert_eq!(strand_bias_phred(0, 0, 0, 0), 0);
    }

    fn mk(start: u32, flags: u16, cigar: Vec<(Kind, u32)>) -> LiveRead {
        LiveRead::new(b"ACGTACGT", &[30; 8], cigar.into_iter().collect(), start, 0, 60, 0, flags)
    }

    #[test]
    fn strand_counts_and_softclip_are_recorded() {
        let mut fwd = mk(0, 0, vec![(Kind::Match, 8)]);
        let mut rev = mk(0, 0x10, vec![(Kind::Match, 8)]);
        let mut clipped = mk(0, 0, vec![(Kind::SoftClip, 2), (Kind::Match, 6)]);
        let mut site = PileupSite::default();
        site.reset(0, 1, 1);
        observe(&mut fwd, 1, 13, true, &mut site);
        observe(&mut rev, 1, 13, true, &mut site);
        observe(&mut clipped, 1, 13, true, &mut site);
        let s = &site.per_sample[0];
        // fwd/rev read base at pos 1 is 'C' (index 1); clipped read: query index 3 = 'T'.
        assert_eq!(s.base_counts[1], 2);
        assert_eq!(s.fwd_counts[1], 1);
        assert_eq!(s.rev_count(1), 1);
        assert_eq!(s.base_counts[3], 1);
        assert_eq!(s.n_softclip, 1);
        assert_eq!(s.depth, 3);
        assert_eq!(fwd.seq(), b"ACGTACGT");
        assert_eq!(fwd.qual(), &[30; 8]);
    }

    #[test]
    fn missing_qualities_are_padded() {
        let lr = LiveRead::new(b"ACGT", &[], CigarOps::from_slice(&[(Kind::Match, 4)]), 5, 0, 10, 0, 0);
        assert_eq!(lr.qual(), &[30; 4]);
        assert_eq!(lr.ref_end(), 9);
    }

    #[test]
    fn cursor_walks_cigar_once() {
        // 3M2I3M1D2M over AAACCGGGTT: ref 0-2 = AAA, "CC" inserted after ref 2,
        // ref 3-5 = GGG, ref 6 deleted, ref 7-8 = TT.
        let cigar: CigarOps = [(Kind::Match, 3), (Kind::Insertion, 2), (Kind::Match, 3), (Kind::Deletion, 1), (Kind::Match, 2)]
            .into_iter()
            .collect();
        let mut lr = LiveRead::new(b"AAACCGGGTT", &[20; 10], cigar, 0, 0, 60, 0, 0);
        assert_eq!(lr.ref_end(), 9);
        let mut site = PileupSite::default();
        let mut got = Vec::new();
        for p in 0..9u32 {
            site.reset(0, p, 1);
            observe(&mut lr, p, 0, false, &mut site);
            let s = &site.per_sample[0];
            let base = (0..4).find(|&i| s.base_counts[i] == 1).map(|i| b"ACGT"[i] as char);
            got.push((base, s.ins_alleles.first().map(|a| a.0.clone()), s.del_alleles.first().map(|d| d.0)));
        }
        assert_eq!(got[0], (Some('A'), None, None));
        assert_eq!(got[2], (Some('A'), Some("CC".to_string()), None));
        assert_eq!(got[3], (Some('G'), None, None));
        assert_eq!(got[5], (Some('G'), None, Some(1)));
        assert_eq!(got[6], (None, None, None));
        assert_eq!(got[7], (Some('T'), None, None));
        assert_eq!(got[8], (Some('T'), None, None));
        // The cursor never rewinds: it ends on the last op.
        assert_eq!(lr.cur_op, 4);
    }

    #[test]
    fn engine_emits_only_covered_sites_and_reuses_buffers() {
        let reads = vec![mk(2, 0, vec![(Kind::Match, 4)]), mk(4, 0, vec![(Kind::Match, 4)])];
        let mut seen = Vec::new();
        mpileup_engine_from_records(vec![reads], 0, 0, true, false, None, &mut |site, live| {
            seen.push((site.pos, site.total_depth(), live.len()));
        })
        .unwrap();
        assert_eq!(seen, vec![(2, 1, 1), (3, 1, 1), (4, 2, 2), (5, 2, 2), (6, 1, 1), (7, 1, 1)]);
    }

    #[test]
    fn overlapping_mates_count_once() {
        // Mates of one fragment overlap at ref 4..8 with equal bases: the
        // second mate's qualities drop to 0 there, so depth is 1 not 2.
        let a = mk(0, 0x1 | 0x2, vec![(Kind::Match, 8)]).with_mate(8, 4);
        let b = mk(4, 0x1 | 0x2 | 0x10, vec![(Kind::Match, 8)]).with_mate(8, 0);
        let mut seen = Vec::new();
        mpileup_engine_from_records(vec![vec![a.clone(), b.clone()]], 0, 1, true, true, None, &mut |site, _| {
            seen.push((site.pos, site.total_depth()));
        })
        .unwrap();
        assert_eq!(seen[4], (4, 1));
        assert_eq!(seen[7], (7, 1));
        assert_eq!(seen[8], (8, 1));
        // Without overlap detection both mates count.
        let mut seen2 = Vec::new();
        mpileup_engine_from_records(vec![vec![a, b]], 0, 1, true, false, None, &mut |site, _| {
            seen2.push((site.pos, site.total_depth()));
        })
        .unwrap();
        assert_eq!(seen2[4], (4, 2));
    }

    #[test]
    fn tweak_overlap_quality_pools_agreeing_bases() {
        // qname hash bit 1 set: `a` keeps the pooled quality.
        let mut a = mk(0, 0x1, vec![(Kind::Match, 8)]).with_mate(1, 4);
        let mut b = LiveRead::new(b"ACGTAAAA", &[20; 8], CigarOps::from_slice(&[(Kind::Match, 8)]), 4, 0, 60, 0, 0x1).with_mate(1, 0);
        tweak_overlap_quality(&mut a, &mut b);
        // a[4..8] = "ACGT" vs b[0..4] = "ACGT": agree -> a gets 30+20, b gets 0.
        assert_eq!(&a.qual()[4..8], &[50, 50, 50, 50]);
        assert_eq!(&b.qual()[..4], &[0, 0, 0, 0]);
        assert_eq!(&b.qual()[4..], &[20, 20, 20, 20]);
        // Hash bit clear: `b` keeps it instead.
        let mut a = mk(0, 0x1, vec![(Kind::Match, 8)]).with_mate(2, 4);
        let mut b = LiveRead::new(b"ACGTAAAA", &[20; 8], CigarOps::from_slice(&[(Kind::Match, 8)]), 4, 0, 60, 0, 0x1).with_mate(2, 0);
        tweak_overlap_quality(&mut a, &mut b);
        assert_eq!(&a.qual()[4..8], &[0, 0, 0, 0]);
        assert_eq!(&b.qual()[..4], &[50, 50, 50, 50]);
    }
