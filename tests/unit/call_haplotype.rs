use super::*;

const REFW: &[u8] = b"ACGTGGCCTTAAGGCCTTAACCGGTTACGTAC";

#[test]
fn discovers_deletion() {
    let mut read = Vec::new();
    read.extend_from_slice(&REFW[..12]);
    read.extend_from_slice(&REFW[14..]); // 2bp deleted at offset 12
    let (ind, mm) = discover_indel(&read, REFW).expect("should find a deletion");
    assert!(!ind.is_ins);
    assert_eq!(ind.len, 2);
    assert_eq!(ind.ref_off, 12);
    assert_eq!(mm, 0, "clean deletion read has no residual mismatches");
    assert_eq!(apply_indel(REFW, &ind), read, "apply should reconstruct the read");
}

#[test]
fn discovers_insertion() {
    let mut read = Vec::new();
    read.extend_from_slice(&REFW[..12]);
    read.extend_from_slice(b"TT");
    read.extend_from_slice(&REFW[12..]);
    let (ind, mm) = discover_indel(&read, REFW).expect("should find an insertion");
    assert!(ind.is_ins);
    assert_eq!(ind.len, 2);
    assert_eq!(ind.ref_off, 12);
    assert_eq!(ind.bases, b"TT");
    assert_eq!(mm, 0);
    assert_eq!(apply_indel(REFW, &ind), read);
}

#[test]
fn residual_mismatch_is_counted() {
    // a deletion read carrying one extra SNP -> residual mm == 1 (the clean-fit signal)
    let mut read = Vec::new();
    read.extend_from_slice(&REFW[..12]);
    read.extend_from_slice(&REFW[14..]);
    read[18] = if read[18] == b'A' { b'C' } else { b'A' };
    let (ind, mm) = discover_indel(&read, REFW).expect("should still find the deletion");
    assert!(!ind.is_ins);
    assert_eq!(ind.len, 2);
    assert_eq!(mm, 1, "one SNP -> one residual mismatch");
}

#[test]
fn clean_read_finds_no_indel() {
    assert!(discover_indel(REFW, REFW).is_none(), "exact read should yield no indel");
}

#[test]
fn haplotype_pls_separate_carriers_from_reference_reads() {
    use crate::bam::pileup::CigarOps;
    use noodles_sam::alignment::record::cigar::op::Kind;
    let mut alt = Vec::new();
    alt.extend_from_slice(&REFW[..12]);
    alt.extend_from_slice(&REFW[14..]);
    let call = AssembledCall {
        pos1: 12,
        ref_str: "GG".into(),
        alt_str: "G".into(),
        support: 2,
        total: 4,
        win_lo: 0,
        win_hi: REFW.len() as u32,
        hap_ref: REFW.to_vec(),
        hap_alt: alt.clone(),
    };
    let mk = |seq: &[u8], cigar: Vec<(Kind, u32)>, sample: usize| {
        LiveRead::new(seq, &vec![35; seq.len()], cigar.into_iter().collect::<CigarOps>(), 0, 0, 60, sample, 0)
    };
    let reads = vec![
        mk(REFW, vec![(Kind::Match, REFW.len() as u32)], 0),
        mk(REFW, vec![(Kind::Match, REFW.len() as u32)], 0),
        mk(&alt, vec![(Kind::Match, 12), (Kind::Deletion, 2), (Kind::Match, 18)], 1),
        mk(&alt, vec![(Kind::Match, 12), (Kind::Deletion, 2), (Kind::Match, 18)], 1),
    ];
    let hs = haplotype_pls(&reads, 2, &call);
    assert_eq!(hs[0].pl[0], 0, "{:?}", hs[0]);
    assert!(hs[0].pl[2] > 20, "{:?}", hs[0]);
    assert_eq!((hs[0].n_ref, hs[0].n_alt), (2, 0));
    assert_eq!(hs[1].pl[2], 0);
    assert!(hs[1].pl[0] > 20, "{:?}", hs[1]);
    assert_eq!((hs[1].n_ref, hs[1].n_alt), (0, 2));
}
