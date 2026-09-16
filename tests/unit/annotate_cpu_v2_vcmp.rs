use super::*;

fn d<'a>(a: &'a str, b: &'a str) -> Option<RefDiff<'a>> {
    diff_refs(a.as_bytes(), b.as_bytes())
}
fn m(a: &str, b: &str, diff: &RefDiff<'_>) -> bool {
    matches_allele(a.as_bytes(), b.as_bytes(), diff)
}

#[test]
fn equal_refs_collapse_to_exact_match() {
    let diff = d("A", "A").unwrap();
    assert_eq!(diff.ndref, 0);
    assert!(m("T", "T", &diff));
    assert!(m("T", "t", &diff));
    assert!(!m("T", "C", &diff));
    assert!(!m("AT", "T", &diff));
}

#[test]
fn incompatible_refs_return_none() {
    assert!(d("ACG", "ATG").is_none());
    assert!(d("A", "T").is_none());
}

#[test]
fn left_padded_insert_matches() {
    let diff = d("A", "ATC").unwrap();
    assert!(diff.ndref < 0);
    assert_eq!(diff.tail, b"TC");
    assert!(m("AT", "ATTC", &diff));
    assert!(!m("AT", "AGTC", &diff));
    let diff2 = d("ATC", "A").unwrap();
    assert!(diff2.ndref > 0);
    assert!(m("ATTC", "AT", &diff2));
}

#[test]
fn left_padded_delete_matches() {
    let diff = d("AGT", "AG").unwrap();
    assert!(diff.ndref > 0);
    assert_eq!(diff.tail, b"T");
    assert!(!m("A", "A", &diff));
}

#[test]
fn snp_with_padding_difference() {
    let diff = d("A", "AT").unwrap();
    assert_eq!(diff.ndref, -1);
    assert_eq!(diff.tail, b"T");
    assert!(m("C", "CT", &diff));
    assert!(!m("C", "GT", &diff));
}

#[test]
fn find_allele_returns_first_match_index() {
    let diff = d("A", "AT").unwrap();
    let als1: Vec<&[u8]> = vec![b"G", b"C", b"T"];
    assert_eq!(find_allele(&als1, b"CT", &diff), Some(1));
    assert_eq!(find_allele(&als1, b"AT", &diff), None);
}

#[test]
fn case_insensitive_throughout() {
    let diff = d("acg", "ACG").unwrap();
    assert_eq!(diff.ndref, 0);
    assert!(m("att", "ATT", &diff));
}

#[test]
fn canonical_cores_ignore_padding_and_anchor_side() {
    assert_eq!(
        canonical_allele(100, b"GAT", b"GCT"),
        (101, &b"A"[..], &b"C"[..])
    );
    assert_eq!(
        canonical_allele(100, b"GATC", b"GC"),
        (101, &b"AT"[..], &b""[..])
    );
    assert_eq!(
        canonical_allele(101, b"ATC", b"C"),
        (101, &b"AT"[..], &b""[..])
    );
    assert_eq!(
        canonical_allele(100, b"A", b"AT"),
        (101, &b""[..], &b"T"[..])
    );
    assert_eq!(
        canonical_allele(100, b"AC", b"ATC"),
        (101, &b""[..], &b"T"[..])
    );
    assert_eq!(canonical_allele(5, b"A", b"C"), (5, &b"A"[..], &b"C"[..]));
    assert_eq!(canonical_allele(5, b"A", b"A"), (5, &b""[..], &b""[..]));
    // Re-expressed over another record's REF.
    assert_eq!(
        reframe_allele(100, b"GAT", 101, b"A", b"C").as_deref(),
        Some(&b"GCT"[..])
    );
    assert_eq!(reframe_allele(101, b"ATC", 200, b"GATC", b"GC"), None);
    assert_eq!(
        reframe_allele(101, b"ATC", 100, b"GATC", b"GC").as_deref(),
        Some(&b"C"[..])
    );
    assert_eq!(
        reframe_allele(100, b"A", 100, b"AC", b"ATC").as_deref(),
        Some(&b"AT"[..])
    );
}
