use super::*;

#[test]
fn repeat_units_are_the_smallest_period() {
    assert_eq!(repeat_unit_len(b"A"), 1);
    assert_eq!(repeat_unit_len(b"AAAA"), 1);
    assert_eq!(repeat_unit_len(b"CACA"), 2);
    assert_eq!(repeat_unit_len(b"CAG"), 3);
    assert_eq!(repeat_unit_len(b"CAGCAGCAG"), 3);
    assert_eq!(repeat_unit_len(b"ACGT"), 4);
}

#[test]
fn copies_count_both_sides_of_the_event() {
    // Deleting one A from a run of six: five copies remain around/inside the tract.
    let left = b"GGCAAA";
    let right = b"AAAGGC";
    let c = str_context(b"A", left, right, false);
    assert_eq!((c.unit, c.copies, c.insertion), (1, 6, false));
    // Inserting CA next to three CA copies.
    let c = str_context(b"CA", b"TTCACACA", b"GGT", true);
    assert_eq!((c.unit, c.copies), (2, 3));
    // Non-repetitive insertion.
    let c = str_context(b"GATTACA", b"CCCC", b"GGGG", true);
    assert_eq!((c.unit, c.copies), (7, 0));
}

#[test]
fn error_grows_with_copies_and_is_steepest_for_homopolymers() {
    let hp = |n| StrContext {
        unit: 1,
        copies: n,
        insertion: false,
    };
    let di = |n| StrContext {
        unit: 2,
        copies: n,
        insertion: false,
    };
    assert!(error_rate(&hp(1)) < error_rate(&hp(4)));
    assert!(error_rate(&hp(8)) > error_rate(&di(8)));
    assert_eq!(indel_quality(&hp(0)), 38);
    assert_eq!(indel_quality(&hp(1)), 38);
    assert!(indel_quality(&hp(12)) < indel_quality(&hp(6)));
    assert_eq!(indel_quality(&hp(40)), 8);
    assert!(!error_prone(&hp(2)));
    assert!(error_prone(&hp(8)));
    assert!(error_prone(&di(9)));
}
