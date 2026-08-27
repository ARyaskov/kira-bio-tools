    use super::*;
    use noodles_sam::alignment::record::cigar::op::Kind;

    fn mk_read(seq: &[u8], cigar: Vec<(Kind, u32)>, start: u32) -> LiveRead {
        let ref_end_cached = {
            let mut e = start;
            for &(k, l) in &cigar {
                if matches!(k, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip) {
                    e += l;
                }
            }
            e
        };
        LiveRead {
            seq: seq.to_vec(),
            qual: vec![30; seq.len()],
            cigar_pairs: cigar,
            ref_start: start,
            ref_end_cached,
            ref_id: 0, mapq: 60, sample_idx: 0, flags: 0,
        }
    }

    #[test]
    fn insertion_apply() {
        let r = b"AAAAATTTTT";
        let out = apply_insertion(r, 5, b"CC");
        assert_eq!(&out, b"AAAAACCTTTTT");
    }

    #[test]
    fn deletion_apply() {
        let r = b"AAAAATTTTT";
        let out = apply_deletion(r, 5, 2);
        assert_eq!(&out, b"AAAAATTT");
    }

    #[test]
    fn perfect_match_score_high() {
        let s = align_score(b"ACGT", b"ACGT");
        assert!(s >= MATCH * 4);
    }

    #[test]
    fn alt_haplotype_wins_for_ins_read() {
        let read = mk_read(b"AAAAACCTTTTT", vec![(Kind::Match, 5), (Kind::Insertion, 2), (Kind::Match, 5)], 100);
        let refw = b"AAAAATTTTT";
        let cand = IndelCandidate { chr: "1".into(), pos: 104, ref_base: b'A',
            kind: IndelKind::Insertion, seq: b"CC".to_vec(), length: 2 };
        let scores = realign_reads_at_indel(&cand, &[&read], refw, 100);
        assert_eq!(scores.len(), 2);
        let alt_score = scores[1].log_prob;
        let ref_score = scores[0].log_prob;
        assert!(alt_score >= ref_score, "ALT should match better, alt={alt_score}, ref={ref_score}");
    }
