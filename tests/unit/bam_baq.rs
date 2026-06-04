    use super::*;

    #[test]
    fn capping_no_panic_empty() {
        let mut q: Vec<u8> = Vec::new();
        apply_baq_capping(&[], &mut q, &[(Kind::Match, 0)], 30);
    }

    #[test]
    fn caps_near_insertion() {
        let mut qual = vec![60u8; 20];
        let cigar = vec![(Kind::Match, 10), (Kind::Insertion, 2), (Kind::Match, 8)];
        apply_baq_capping(&[], &mut qual, &cigar, 30);
        for q in &qual[5..15] { assert!(*q <= 30); }
    }

    #[test]
    fn hmm_perfect_match_high_posterior() {
        let seq = b"ACGTACGTAC";
        let mut qual = vec![40u8; seq.len()];
        let cigar = vec![(Kind::Match, seq.len() as u32)];
        let refseq = b"ACGTACGTAC";
        apply_baq_hmm(seq, &mut qual, &cigar, refseq, 5);
        for q in &qual { assert!(*q >= 25, "expected high quality for perfect match, got {q}"); }
    }

    #[test]
    fn hmm_fallback_no_ref() {
        let seq = b"ACGT";
        let mut qual = vec![40u8; seq.len()];
        let cigar = vec![(Kind::Match, 4)];
        apply_baq_hmm(seq, &mut qual, &cigar, &[], 5);
    }
