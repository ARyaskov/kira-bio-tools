use fxhash::FxHashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::annotate::structs::ani::AniEntry;
use crate::util::chr_name_to_id;

pub fn make_position_key(chr_id: u8, pos: u32, ref_allele: &str, alt: &str) -> u64 {
    use crate::util::fast_hash64;

    let mut h = (chr_id as u64) << 32 | (pos as u64);
    h ^= fast_hash64(ref_allele.as_bytes());
    h ^= fast_hash64(alt.as_bytes());
    h
}

pub fn insert_or_update_entry(
    key: u64,
    entry: AniEntry,
    entries_map: &mut FxHashMap<u64, (AniEntry, usize)>,
    insertion_order: &mut usize,
    duplicates_skipped: &AtomicUsize,
    debug: bool,
    chr: &str,
    pos: u32,
    rf: &str,
    alt: &str,
) {
    if entries_map.contains_key(&key) {
        duplicates_skipped.fetch_add(1, Ordering::Relaxed);

        if debug {
            eprintln!(
                "[ani-build] Overwriting duplicate: {}:{} {} {}",
                chr, pos, rf, alt
            );
        }
    }

    entries_map.insert(key, (entry, *insertion_order));
    *insertion_order += 1;
}

pub fn parse_chrom_and_pos(chr: &str, pos_str: &str) -> Option<(u8, u32)> {
    let chr_id = chr_name_to_id(chr)?;
    let pos = pos_str.parse::<u32>().ok()?;
    Some((chr_id, pos))
}
