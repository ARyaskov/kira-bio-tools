
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
fn inserted_and_clipped_bases_keep_their_quality() {
    // `sam_prob_realn` only caps aligned bases; inserted and soft-clipped
    // ones are untouched.
    let mut read = Vec::new();
    read.extend_from_slice(&REF[10..18]);
    read.extend_from_slice(b"TT");
    read.extend_from_slice(&REF[18..26]);
    let mut qual = vec![40u8; read.len()];
    let cigar = vec![(Kind::Match, 8), (Kind::Insertion, 2), (Kind::Match, 8)];
    assert!(apply_baq_hmm(&read, &mut qual, &cigar, REF, 0, 10));
    assert_eq!(qual[8], 40);
    assert_eq!(qual[9], 40);
    assert!(qual[..8].iter().any(|&q| q > 0));

    let mut clipped = b"GGGG".to_vec();
    clipped.extend_from_slice(&REF[10..26]);
    let mut qual = vec![35u8; clipped.len()];
    let cigar = vec![(Kind::SoftClip, 4), (Kind::Match, 16)];
    assert!(apply_baq_hmm(&clipped, &mut qual, &cigar, REF, 0, 10));
    assert!(qual[..4].iter().all(|&q| q == 35), "{qual:?}");
}

#[test]
fn bq_tag_applies_stored_offsets() {
    let mut qual = vec![40u8, 30, 20, 10];
    // 64 = no change, 74 = minus 10, 200 = floors at 0.
    assert!(apply_bq_tag(&mut qual, &[64, 74, 64, 200]));
    assert_eq!(qual, vec![40, 20, 20, 0]);
    assert!(!apply_bq_tag(&mut qual, &[64, 64]));
}

#[test]
fn reference_skip_and_missing_window_are_no_ops() {
    let read = &REF[10..20];
    let mut qual = vec![40u8; read.len()];
    assert!(!apply_baq_hmm(
        read,
        &mut qual,
        &[(Kind::Match, 5), (Kind::Skip, 3), (Kind::Match, 5)],
        REF,
        0,
        10
    ));
    assert!(!apply_baq_hmm(
        read,
        &mut qual,
        &[(Kind::Match, 10)],
        &[],
        0,
        10
    ));
    assert!(qual.iter().all(|&q| q == 40));
}
