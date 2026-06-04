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
