use fxhash::FxHashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::annotate::structs::ani::{AniEntry, make_variant_key};

/// Compose a `u64` variant key for the MPH builder. Thin wrapper over
/// [`make_variant_key`] — kept under the historic name `make_position_key`
/// so existing call sites read naturally.
#[inline]
pub fn make_position_key(chr_id: u32, pos: u32, ref_allele: &str, alt: &str) -> u64 {
    make_variant_key(chr_id, pos, ref_allele.as_bytes(), alt.as_bytes())
}

/// Map value: the variant entry plus its insertion order plus a verification
/// triplet (chr_id, pos, ref+alt offset) — the last lets us distinguish a
/// *real* duplicate (same variant) from a *hash collision* (different
/// variants colliding on the u64 MPH key, which would silently corrupt
/// output).
pub struct EntryEntry {
    pub entry: AniEntry,
    pub order: usize,
    pub ref_alt: Vec<u8>, // ref + 0xFF + alt — small (~10 B for SNPs)
}

pub fn insert_or_update_entry(
    key: u64,
    entry: AniEntry,
    entries_map: &mut FxHashMap<u64, EntryEntry>,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    collisions_detected: &AtomicUsize,
    debug: bool,
    chr: &str,
    pos: u32,
    rf: &str,
    alt: &str,
) {
    let mut ref_alt = Vec::with_capacity(rf.len() + 1 + alt.len());
    ref_alt.extend_from_slice(rf.as_bytes());
    ref_alt.push(0xFF);
    ref_alt.extend_from_slice(alt.as_bytes());

    if let Some(existing) = entries_map.get(&key) {
        let same_variant = existing.entry.chr_id == entry.chr_id
            && existing.entry.pos == entry.pos
            && existing.ref_alt == ref_alt;
        if same_variant {
            duplicates_skipped.fetch_add(1, Ordering::Relaxed);
            if debug {
                eprintln!("[ani-build] True duplicate skipped: {chr}:{pos} {rf}>{alt}");
            }
            return;
        } else {
            collisions_detected.fetch_add(1, Ordering::Relaxed);
            eprintln!(
                "[ani-build] WARNING: hash collision on {key:016x}: \
                 keeping previous entry, dropping ({chr}:{pos} {rf}>{alt})"
            );
            return;
        }
    }

    entries_map.insert(
        key,
        EntryEntry {
            entry,
            order: *insertion_order,
            ref_alt,
        },
    );
    *insertion_order += 1;
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_builder_v2_entry_processing.rs"]
mod tests;
