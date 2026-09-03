    use super::*;

    const REF: &[u8] = b"TTGACCGTAGACGTACGTACGTACGTACGTACGGATCCATTG";

    #[test]
    fn perfect_match_keeps_most_of_the_quality() {
        let read = &REF[10..30];
        let mut qual = vec![40u8; read.len()];
        let cigar = vec![(Kind::Match, read.len() as u32)];
        assert!(apply_baq_hmm(read, &mut qual, &cigar, REF, 0, 10));
        assert!(qual.iter().all(|&q| q > 0), "{qual:?}");
        let high = qual.iter().filter(|&&q| q >= 20).count();
        assert!(high * 2 >= qual.len(), "{qual:?}");
    }

    #[test]
    fn inserted_bases_get_zero_quality() {
        let mut read = Vec::new();
        read.extend_from_slice(&REF[10..18]);
        read.extend_from_slice(b"TT");
        read.extend_from_slice(&REF[18..26]);
        let mut qual = vec![40u8; read.len()];
        let cigar = vec![(Kind::Match, 8), (Kind::Insertion, 2), (Kind::Match, 8)];
        assert!(apply_baq_hmm(&read, &mut qual, &cigar, REF, 0, 10));
        assert_eq!(qual[8], 0);
        assert_eq!(qual[9], 0);
        assert!(qual[..8].iter().any(|&q| q > 0));
    }

    #[test]
    fn reference_skip_and_missing_window_are_no_ops() {
        let read = &REF[10..20];
        let mut qual = vec![40u8; read.len()];
        assert!(!apply_baq_hmm(read, &mut qual, &[(Kind::Match, 5), (Kind::Skip, 3), (Kind::Match, 5)], REF, 0, 10));
        assert!(!apply_baq_hmm(read, &mut qual, &[(Kind::Match, 10)], &[], 0, 10));
        assert!(qual.iter().all(|&q| q == 40));
    }
