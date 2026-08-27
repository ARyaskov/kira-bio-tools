use super::*;

fn q(n: usize) -> Vec<u8> {
    vec![40u8; n]
}

#[test]
fn exact_match_beats_mismatches() {
    let hap = b"ACGTGGCCTTAAGGCCTTAACCGGTTACGTAC";
    let exact = &hap[..24];
    let mut mism = exact.to_vec();
    mism[5] = b'A';
    mism[12] = b'T';
    mism[18] = b'G';
    let le = read_vs_hap_loglik(exact, &q(24), hap);
    let lm = read_vs_hap_loglik(&mism, &q(24), hap);
    assert!(le > lm, "exact {le} should beat 3-mismatch {lm}");
}

#[test]
fn deletion_read_prefers_alt_haplotype() {
    let refh = b"ACGTGGCCTTAAGGCCTTAACCGGTTACGTAC";
    // read = ref with the 2bp at [12,14) deleted
    let mut read = Vec::new();
    read.extend_from_slice(&refh[..12]);
    read.extend_from_slice(&refh[14..]);
    let alt = read.clone(); // haplotype carrying the deletion
    let lr = loglik_ratio(&read, &q(read.len()), refh, &alt);
    assert!(lr > 5.0, "deletion read should strongly prefer alt hap, LR={lr}");
}

#[test]
fn insertion_read_prefers_alt_haplotype() {
    let refh = b"ACGTGGCCTTAAGGCCTTAACCGGTTACGTAC";
    let mut read = Vec::new();
    read.extend_from_slice(&refh[..12]);
    read.extend_from_slice(b"TT"); // 2bp insertion
    read.extend_from_slice(&refh[12..28]);
    let alt = read.clone();
    let lr = loglik_ratio(&read, &q(read.len()), refh, &alt);
    assert!(lr > 5.0, "insertion read should prefer alt hap, LR={lr}");
}

#[test]
fn clean_ref_read_rejects_spurious_indel() {
    // Specificity: a read matching the reference must NOT prefer a spurious-indel haplotype.
    let refh = b"ACGTGGCCTTAAGGCCTTAACCGGTTACGTAC";
    let read = &refh[2..30]; // 28bp, matches ref
    let mut alt = Vec::new(); // spurious: ref with a 1bp deletion at 15
    alt.extend_from_slice(&refh[..15]);
    alt.extend_from_slice(&refh[16..]);
    let lr = loglik_ratio(read, &q(read.len()), refh, &alt);
    assert!(lr < 0.0, "clean ref read must reject spurious indel, LR={lr}");
}
