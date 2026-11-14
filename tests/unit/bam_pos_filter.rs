use super::*;

fn map(positions: impl IntoIterator<Item = (usize, Vec<u32>)>) -> InterestingMap {
    InterestingMap {
        per_ref: positions
            .into_iter()
            .map(|(rid, positions)| (rid, Positions::Plain(positions)))
            .collect(),
    }
}

#[test]
fn frozen_candidates_match_membership_and_successor_oracles() {
    let positions: Vec<_> = (0..4096u32).map(|i| i * 3).collect();
    let mut candidates = map([(2, positions.clone()), (7, vec![0, u32::MAX])]);
    candidates.freeze();
    assert!(matches!(
        candidates.per_ref.get(&2),
        Some(Positions::Compact(_))
    ));
    assert!(matches!(
        candidates.per_ref.get(&7),
        Some(Positions::Plain(_))
    ));
    assert_eq!(candidates.total(), positions.len() + 2);
    for pos in (0..13000).chain([u32::MAX - 1, u32::MAX]) {
        assert_eq!(
            candidates.contains(2, pos),
            positions.binary_search(&pos).is_ok(),
            "{pos}"
        );
        assert_eq!(
            candidates.next_at_or_after(2, pos),
            positions.iter().copied().find(|&p| p >= pos),
            "{pos}"
        );
    }
    assert_eq!(candidates.next_at_or_after(7, 1), Some(u32::MAX));
    assert!(candidates.contains(7, u32::MAX));
    assert_eq!(candidates.next_at_or_after(9, 0), None);
    assert!(!candidates.contains(9, 0));
}

#[test]
fn freezing_requires_savings_against_four_byte_positions() {
    let mut candidates = map([
        (0, vec![1, 2, 3]),
        (1, (0..256).map(|i| u32::MAX - 255 + i).collect()),
        (2, (0..4096).collect()),
    ]);
    candidates.freeze();
    assert!(matches!(
        candidates.per_ref.get(&0),
        Some(Positions::Plain(_))
    ));
    assert!(matches!(
        candidates.per_ref.get(&1),
        Some(Positions::Plain(_))
    ));
    let Some(Positions::Compact(compact)) = candidates.per_ref.get(&2) else {
        panic!("dense positions should compact");
    };
    assert!(compact.memory_usage() + std::mem::size_of::<Positions>() < 4096 * 4);
    candidates.freeze();
    assert_eq!(candidates.total(), 3 + 256 + 4096);
}

#[test]
fn merging_frozen_candidates_preserves_union_and_can_freeze_again() {
    let even: Vec<_> = (0..4096).map(|i| i * 2).collect();
    let odd: Vec<_> = (0..4096).map(|i| i * 2 + 1).collect();
    let mut left = map([(1, even), (2, vec![u32::MAX])]);
    let mut right = map([(1, odd), (2, vec![u32::MAX]), (3, vec![9])]);
    left.freeze();
    right.freeze();
    left.merge(right);
    left.freeze();
    assert_eq!(left.total(), 8194);
    for pos in 0..8192 {
        assert!(left.contains(1, pos));
        assert_eq!(left.next_at_or_after(1, pos), Some(pos));
    }
    assert_eq!(left.next_at_or_after(1, 8192), None);
    assert!(left.contains(2, u32::MAX));
    assert_eq!(left.next_at_or_after(3, 0), Some(9));
    assert!(!left.contains(4, 9));
}

#[test]
fn known_site_set_keeps_public_construction_and_coordinate_semantics() {
    let known = KnownSites::Set(crate::fx::FxHashSet::from_iter([
        7,
        (2u64 << 32) | u32::MAX as u64,
    ]));
    assert_eq!(known.len(), 2);
    assert!(known.contains(0, 7));
    assert!(!known.contains(0, 8));
    assert!(known.contains(2, u32::MAX));
    assert!(!known.contains(1, u32::MAX));
}
