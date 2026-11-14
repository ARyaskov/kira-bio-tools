use super::*;
use crate::kbi::KbiBuilder;
use std::collections::BTreeMap;
use std::io::Write;

fn build_index(entries: &[(u32, u32, u64)]) -> KbiIndex {
    let mut contigs = ContigDict::new();
    for name in ["chr1", "chr2", "empty"] {
        contigs.insert(name);
    }
    let mut builder = KbiBuilder::new();
    builder.set_contigs(contigs);
    for &(chr, pos, offset) in entries {
        builder.add(GenomicKey::new(chr, pos), offset);
    }
    builder.build().unwrap()
}

fn saved_bytes(index: &KbiIndex) -> Vec<u8> {
    let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    index.save(&path).unwrap();
    std::fs::read(path).unwrap()
}

fn load_bytes(bytes: &[u8]) -> Result<KbiIndex> {
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(bytes).unwrap();
    file.flush().unwrap();
    KbiIndex::load(file.path())
}

fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

fn u64_at(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

#[test]
fn kbi_public_queries_match_first_inserted_btree_oracle() {
    let mut entries = vec![
        (1, 0, 900),
        (0, u32::MAX, 400),
        (0, 10, 100),
        (1, u32::MAX, 999),
    ];
    let mut state = 17u64;
    for i in 0..4096u64 {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        entries.push((
            (state >> 63) as u32,
            ((state >> 32) as u32) % 1400,
            1_000_000 - i,
        ));
    }
    entries.extend_from_slice(&[(0, 10, 777), (0, u32::MAX, 888), (1, 0, 1)]);
    let mut oracle = BTreeMap::new();
    for &(chr, pos, offset) in &entries {
        oracle.entry((chr, pos)).or_insert(offset);
    }
    let built = build_index(&entries);
    let loaded = load_bytes(&saved_bytes(&built)).unwrap();
    for index in [&built, &loaded] {
        assert_eq!(index.len(), oracle.len());
        assert!(!index.is_empty());
        let query_keys: Vec<_> = (0..3)
            .flat_map(|chr| (0..1500).map(move |pos| GenomicKey::new(chr, pos)))
            .collect();
        let expected: Vec<_> = query_keys
            .iter()
            .map(|k| oracle.get(&(k.chr(), k.position())).copied())
            .collect();
        assert_eq!(index.get_batch(&query_keys), expected);
        for (key, expected) in query_keys.into_iter().zip(expected) {
            assert_eq!(index.get(key), expected);
            assert_eq!(index.contains(key), expected.is_some());
        }
        for chr in 0..3 {
            for (start, end) in [
                (0, 0),
                (0, 100),
                (10, 10),
                (11, 13),
                (1400, u32::MAX),
                (u32::MAX, u32::MAX),
                (0, u32::MAX),
                (100, 10),
            ] {
                let expected: Vec<_> = oracle
                    .iter()
                    .filter_map(|(&(c, pos), &offset)| {
                        (c == chr && start <= pos && pos <= end).then_some((pos, offset))
                    })
                    .collect();
                assert_eq!(index.range(chr, start, end), expected);
                assert_eq!(index.has_range(chr, start, end), !expected.is_empty());
                assert_eq!(
                    index.range_by_name(["chr1", "chr2", "empty"][chr as usize], start, end),
                    expected
                );
            }
        }
        assert!(index.range_by_name("absent", 0, u32::MAX).is_empty());
        assert_eq!(index.contigs_with_records(), ["chr1", "chr2"]);
        assert_eq!(index.get(GenomicKey::new(0, 10)), Some(100));
        assert_eq!(index.get(GenomicKey::new(0, u32::MAX)), Some(400));
        assert_eq!(index.get(GenomicKey::new(1, 0)), Some(900));
    }
}

#[test]
fn kbi_empty_and_singleton_boundaries() {
    assert!(matches!(
        KbiBuilder::new().build(),
        Err(KbiError::EmptyDataset)
    ));
    let mut builder = KbiBuilder::new();
    builder.add(GenomicKey::from_u64(u64::MAX), u64::MAX);
    let built = builder.build().unwrap();
    let loaded = load_bytes(&saved_bytes(&built)).unwrap();
    for index in [&built, &loaded] {
        assert_eq!(index.get(GenomicKey::from_u64(u64::MAX)), Some(u64::MAX));
        assert_eq!(index.get(GenomicKey::from_u64(u64::MAX - 1)), None);
        assert_eq!(
            index.range(u32::MAX, u32::MAX, u32::MAX),
            [(u32::MAX, u64::MAX)]
        );
        assert!(index.has_range(u32::MAX, 0, u32::MAX));
        assert!(index.range(u32::MAX, 2, 1).is_empty());
        assert!(!index.has_range(u32::MAX, 2, 1));
    }
}

#[test]
fn kbi_compact_keys_remove_duplicate_keys_and_slot_table() {
    let entries: Vec<_> = (0..20_000)
        .map(|pos| (0, pos * 2, u64::MAX - pos as u64))
        .collect();
    let index = build_index(&entries);
    let bytes = saved_bytes(&index);
    assert!(matches!(&index.keys, KbiKeys::Compact(_)));
    assert_eq!(bytes[80], 1);
    assert_eq!(&bytes[..8], b"KBIV0005");
    assert_eq!(u64_at(&bytes, 40), u64_at(&bytes, 48));
    assert_eq!(u64_at(&bytes, 56), u64_at(&bytes, 72));
    assert_eq!(u64_at(&bytes, 64), 0);
    assert!(index.memory_usage() < entries.len() * 16);
    assert!(bytes.len() < entries.len() * 20);
    let loaded = load_bytes(&bytes).unwrap();
    assert!(loaded.memory_usage() < entries.len() * 16);
    assert!(matches!(&loaded.keys, KbiKeys::Compact(_)));
    let expected: Vec<_> = entries
        .iter()
        .map(|&(_, pos, offset)| (pos, offset))
        .collect();
    assert_eq!(loaded.range(0, 0, u32::MAX), expected);
    assert_eq!(
        loaded.range(0, 10, 16),
        [
            (10, u64::MAX - 5),
            (12, u64::MAX - 6),
            (14, u64::MAX - 7),
            (16, u64::MAX - 8)
        ]
    );
}

#[test]
fn kbi_v5_rejects_malformed_header_and_metadata() {
    let valid = saved_bytes(&build_index(&[(0, 10, 100), (1, 20, 3)]));
    for (offset, value) in [
        (16, 0),
        (16, u64::MAX),
        (24, 0),
        (24, u64::MAX),
        (32, 0),
        (32, u64::MAX),
        (40, 0),
        (40, u64::MAX),
        (48, 0),
        (48, u64::MAX),
        (56, 0),
        (56, u64::MAX),
        (64, 1),
        (72, 0),
        (72, u64::MAX),
    ] {
        let mut malformed = valid.clone();
        put_u64(&mut malformed, offset, value);
        assert!(
            load_bytes(&malformed).is_err(),
            "offset {offset}, value {value}"
        );
    }
    for offset in [0usize, 8, 12] {
        let mut malformed = valid.clone();
        malformed[offset] ^= 1;
        assert!(load_bytes(&malformed).is_err(), "offset {offset}");
    }
    for len in [0, 1, 79, 80, valid.len() - 1] {
        assert!(load_bytes(&valid[..len]).is_err(), "truncation {len}");
    }
    let names_start = u64_at(&valid, 72) as usize;
    let mut malformed = valid.clone();
    malformed[names_start..names_start + 4].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.clone();
    malformed[names_start + 4..names_start + 6].copy_from_slice(&u16::MAX.to_le_bytes());
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.clone();
    malformed[names_start + 6] = 0xff;
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.clone();
    malformed[names_start + 12..names_start + 16].copy_from_slice(b"chr1");
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.clone();
    malformed.push(0);
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.clone();
    malformed[80] ^= 1;
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.clone();
    put_u64(&mut malformed, 16, 1);
    let offsets_start = u64_at(&valid, 48);
    malformed.drain(offsets_start as usize..offsets_start as usize + 8);
    put_u64(&mut malformed, 56, offsets_start + 8);
    put_u64(&mut malformed, 72, offsets_start + 8);
    assert!(load_bytes(&malformed).is_err());
}

#[test]
fn kbi_frozen_v4_preserves_exact_offsets_ranges_and_contigs() {
    let bytes = include_bytes!("../fixtures/kbi/legacy_v4.kbi");
    assert_eq!(&bytes[..8], b"KBIV0004");
    let slots_start = u64_at(bytes, 56) as usize;
    let slots: Vec<_> = bytes[slots_start..slots_start + 20]
        .as_chunks::<4>()
        .0
        .iter()
        .map(|chunk| u32::from_le_bytes(*chunk))
        .collect();
    assert_eq!(slots, [4, 2, 0, 1, 3]);
    let index = load_bytes(bytes).unwrap();
    assert_eq!(index.len(), 5);
    assert_eq!(index.contigs().names(), ["chr1", "chr2", "empty"]);
    assert_eq!(
        index.range_by_name("chr1", 0, u32::MAX),
        [(10, 100), (20, 200), (u32::MAX, 400)]
    );
    assert_eq!(
        index.range_by_name("chr2", 0, u32::MAX),
        [(0, 900), (u32::MAX, 999)]
    );
    assert_eq!(index.get(GenomicKey::new(0, 10)), Some(100));
    assert_eq!(index.get(GenomicKey::new(0, 11)), None);
    assert_eq!(index.get(GenomicKey::new(1, u32::MAX)), Some(999));
    assert_eq!(index.contigs_with_records(), ["chr1", "chr2"]);
    let back = load_bytes(&saved_bytes(&index)).unwrap();
    assert_eq!(back.range(0, 0, u32::MAX), index.range(0, 0, u32::MAX));
    assert_eq!(back.range(1, 0, u32::MAX), index.range(1, 0, u32::MAX));
}

#[test]
fn kbi_v4_rejects_malformed_sections_keys_and_exact_slot_mapping() {
    let valid = include_bytes!("../fixtures/kbi/legacy_v4.kbi");
    for (offset, value) in [
        (16, 0),
        (16, u64::MAX),
        (24, u64::MAX),
        (32, 0),
        (40, u64::MAX),
        (48, 0),
        (56, 0),
        (64, 0),
        (64, u64::MAX),
        (72, u64::MAX),
    ] {
        let mut malformed = valid.to_vec();
        put_u64(&mut malformed, offset, value);
        assert!(
            load_bytes(&malformed).is_err(),
            "offset {offset}, value {value}"
        );
    }
    let keys_start = u64_at(valid, 40) as usize;
    let slots_start = u64_at(valid, 56) as usize;
    let mut malformed = valid.to_vec();
    put_u64(&mut malformed, keys_start + 8, u64_at(valid, keys_start));
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.to_vec();
    malformed[keys_start..keys_start + 8]
        .copy_from_slice(&u64_at(valid, keys_start + 8).to_le_bytes());
    malformed[keys_start + 8..keys_start + 16]
        .copy_from_slice(&u64_at(valid, keys_start).to_le_bytes());
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.to_vec();
    malformed[slots_start..slots_start + 4].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.to_vec();
    malformed[slots_start..slots_start + 4].copy_from_slice(&5u32.to_le_bytes());
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.to_vec();
    let first = malformed[slots_start..slots_start + 4].to_vec();
    let second = malformed[slots_start + 4..slots_start + 8].to_vec();
    malformed[slots_start..slots_start + 4].copy_from_slice(&second);
    malformed[slots_start + 4..slots_start + 8].copy_from_slice(&first);
    assert!(load_bytes(&malformed).is_err());
    for len in [79, 80, keys_start, slots_start, valid.len() - 1] {
        assert!(load_bytes(&valid[..len]).is_err(), "truncation {len}");
    }
    let mut malformed = valid.to_vec();
    malformed.push(0);
    assert!(load_bytes(&malformed).is_err());
}

#[test]
fn kbi_v5_rejects_invalid_raw_and_compact_key_payloads() {
    let valid = saved_bytes(&build_index(&[(0, 10, 100), (1, 20, 3)]));
    assert_eq!(valid[80], 0);
    let mut malformed = valid.clone();
    malformed[80] = 2;
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.clone();
    put_u64(&mut malformed, 89, u64_at(&valid, 81));
    assert!(load_bytes(&malformed).is_err());
    let mut malformed = valid.clone();
    put_u64(&mut malformed, 81, u64_at(&valid, 89));
    put_u64(&mut malformed, 89, u64_at(&valid, 81));
    assert!(load_bytes(&malformed).is_err());
    let entries: Vec<_> = (0..4096).map(|pos| (0, pos, pos as u64)).collect();
    let valid = saved_bytes(&build_index(&entries));
    assert_eq!(valid[80], 1);
    for (offset, value) in [(81, 0), (81, 4095), (81, u64::MAX), (89, 0), (98, u64::MAX)] {
        let mut malformed = valid.clone();
        put_u64(&mut malformed, offset, value);
        assert!(
            load_bytes(&malformed).is_err(),
            "offset {offset}, value {value}"
        );
    }
    let mut malformed = valid.clone();
    malformed[97] = 57;
    assert!(load_bytes(&malformed).is_err());
    let index_end = u64_at(&valid, 48) as usize;
    let mut malformed = valid.clone();
    malformed.insert(index_end, 0);
    for offset in [24, 40, 48, 56, 72] {
        put_u64(&mut malformed, offset, u64_at(&valid, offset) + 1);
    }
    assert!(load_bytes(&malformed).is_err());
}

#[test]
fn kbi_rejects_oversized_contig_names_before_replacing_output() {
    let mut contigs = ContigDict::new();
    contigs.insert(&"x".repeat(u16::MAX as usize + 1));
    let mut builder = KbiBuilder::new();
    builder.set_contigs(contigs);
    builder.add(GenomicKey::new(0, 1), 0);
    let index = builder.build().unwrap();
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(b"preserved").unwrap();
    file.flush().unwrap();
    assert!(index.save(file.path()).is_err());
    assert_eq!(std::fs::read(file.path()).unwrap(), b"preserved");
}

#[test]
fn kbi_public_legacy_header_constructor_preserves_v4_layout() {
    let header = KbiHeader::new(3, 7, 3);
    header.validate().unwrap();
    let KbiHeader {
        magic,
        version,
        endian,
        n_entries,
        index_len,
        off_index,
        off_keys,
        off_offsets,
        off_slots,
        n_slots,
        off_names,
    } = header;
    assert_eq!(magic, *b"KBIV0004");
    assert_eq!(version, 4);
    assert_eq!(endian, crate::kbi::structs::ENDIAN_TAG);
    assert_eq!(n_entries, 3);
    assert_eq!(index_len, 7);
    assert_eq!(off_index, 80);
    assert_eq!(off_keys, 87);
    assert_eq!(off_offsets, 111);
    assert_eq!(off_slots, 135);
    assert_eq!(n_slots, 3);
    assert_eq!(off_names, 147);
    assert!(
        KbiHeader::new(usize::MAX, usize::MAX, usize::MAX)
            .validate()
            .is_err()
    );
}
