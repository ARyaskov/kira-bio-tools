use super::*;

#[test]
fn reference_coordinates_resolve_deletion_boundaries() {
    let read = LiveRead::new(
        b"AAAAAA",
        &[30; 6],
        vec![(Kind::Match, 3), (Kind::Deletion, 2), (Kind::Match, 3)].into(),
        10,
        0,
        60,
        0,
        0,
    );
    assert_eq!(query_position(&read, 14, true), (3, 13));
    assert_eq!(query_position(&read, 14, false), (3, 15));
    assert_eq!(query_position(&read, 18, false), (6, 18));
}

#[test]
fn homopolymer_sequence_quality_scales_with_gap_length() {
    let options = IndelOptions {
        min_support: 2,
        min_fraction: 0.05,
        per_sample: false,
        window: 110,
        open_quality: 40,
        extension_quality: 20,
        tandem_quality: 500,
        bias: 1.0,
        minimum_quality: 1,
        ambiguous: 0,
    };
    assert_eq!(sequence_quality(&options, 1, 20), 25);
    assert_eq!(sequence_quality(&options, -2, 20), 50);
    assert_eq!(sequence_quality(&options, 1, 2), 40);
}
