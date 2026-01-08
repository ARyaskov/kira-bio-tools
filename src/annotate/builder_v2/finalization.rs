use anyhow::Result;
use kira_kv_engine::{BuildConfig, Builder};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;

use crate::annotate::builder_v2::StringPool;
use crate::annotate::structs::ani::{AniBlockEntry, AniEntry, AniHeaderV4, ANI_MAGIC, ANI_VERSION};

pub fn finalize_ani_index(
    rows: Vec<(u64, AniEntry)>,
    mut pool: StringPool,
    output: &Path,
    timing: bool,
) -> Result<()> {
    if rows.is_empty() {
        anyhow::bail!("No valid entries to index");
    }

    let n = rows.len();
    let mph_start = std::time::Instant::now();

    let (gamma, rehash_limit) = calculate_mph_params(n);

    if timing {
        print_mph_config(gamma, rehash_limit, n, &rows);
    }

    let keys_bytes: Vec<[u8; 8]> = rows.iter().map(|(k, _)| k.to_le_bytes()).collect();

    if timing {
        verify_key_uniqueness(&keys_bytes, &rows);
    }

    let mph = Builder::new()
        .with_config(BuildConfig {
            gamma,
            rehash_limit: rehash_limit as u32,
            salt: 0x9E3779B185EBCA87,
        })
        .build(keys_bytes.iter().map(|b| b.as_slice()))?;

    if timing {
        eprintln!(
            "[ani-build] MPH construction: {:.3}s",
            mph_start.elapsed().as_secs_f64()
        );
    }

    let entries: Vec<AniEntry> = rows.into_iter().map(|(_, entry)| entry).collect();

    if timing {
        verify_mph_correctness(&mph, &keys_bytes, &entries);
    }

    write_ani_file(output, &mph, &entries, &mut pool, timing)?;
    pool.cleanup();

    Ok(())
}

fn calculate_mph_params(n: usize) -> (f64, usize) {
    match n {
        0..=100_000 => (1.2, 16),
        100_001..=1_000_000 => (1.5, 32),
        1_000_001..=10_000_000 => (2.0, 64),
        _ => (2.5, 100),
    }
}

fn print_mph_config(gamma: f64, rehash_limit: usize, n: usize, rows: &[(u64, AniEntry)]) {
    eprintln!(
        "[ani-build] MPHF config: gamma={:.1}, rehash_limit={}, entries={}",
        gamma, rehash_limit, n
    );
    eprintln!("[ani-build] First 5 keys for verification:");
    for (i, (k, e)) in rows.iter().enumerate().take(5) {
        eprintln!("  [{}] key={:016x} chr={} pos={}", i, k, e.chr_id, e.pos);
    }
}

fn verify_key_uniqueness(keys_bytes: &[[u8; 8]], rows: &[(u64, AniEntry)]) {
    let mut key_set = std::collections::HashSet::new();
    let mut dup_count = 0;

    for (i, key) in keys_bytes.iter().enumerate() {
        let k = u64::from_le_bytes(*key);
        if !key_set.insert(k) {
            dup_count += 1;
            if dup_count <= 5 {
                let (_, e) = &rows[i];
                eprintln!(
                    "[ani-build] WARNING: Duplicate key {:016x} for chr={} pos={}",
                    k, e.chr_id, e.pos
                );
            }
        }
    }

    if dup_count > 0 {
        eprintln!(
            "[ani-build] ERROR: Found {} duplicate keys! MPH will not work correctly!",
            dup_count
        );
    }
}

fn verify_mph_correctness(
    mph: &kira_kv_engine::Mphf,
    keys_bytes: &[[u8; 8]],
    entries: &[AniEntry],
) {
    eprintln!("[ani-build] Verifying MPH lookups (checking 100 random entries):");
    let mut errors = 0;
    let check_count = 100.min(entries.len());
    let step = entries.len() / check_count;

    for i in (0..entries.len()).step_by(step.max(1)).take(check_count) {
        let key = u64::from_le_bytes(keys_bytes[i]);
        let idx = mph.index(&keys_bytes[i]) as usize;
        let retrieved = &entries[idx];
        let expected = &entries[i];

        if retrieved.chr_id != expected.chr_id || retrieved.pos != expected.pos {
            errors += 1;
            if errors <= 5 {
                eprintln!("  [{}] ERROR: key={:016x} -> mph_idx={} -> chr={} pos={} (expected chr={} pos={})", i, key, idx, retrieved.chr_id, retrieved.pos, expected.chr_id, expected.pos);
            }
        } else if i < 3 {
            eprintln!(
                "  [{}] OK: key={:016x} -> mph_idx={} -> chr={} pos={}",
                i, key, idx, retrieved.chr_id, retrieved.pos
            );
        }
    }

    if errors > 0 {
        eprintln!(
            "[ani-build] ERROR: MPH verification FAILED with {} errors out of {} checks!",
            errors, check_count
        );
        eprintln!("[ani-build] This means the index will NOT work correctly!");
    } else {
        eprintln!(
            "[ani-build] MPH verification: All {} checks passed ✓",
            check_count
        );
    }
}

fn write_ani_file(
    output: &Path,
    mph: &kira_kv_engine::Mphf,
    entries: &[AniEntry],
    pool: &mut StringPool,
    timing: bool,
) -> Result<()> {
    let n = entries.len();
    let hdr_size = mem::size_of::<AniHeaderV4>();
    let g_size = mph.g.len() * 4;
    let ent_size = n * mem::size_of::<AniEntry>();
    let block_size = pool.block_size();
    let blocks = pool.blocks();
    let block_index_size = blocks.len() * mem::size_of::<AniBlockEntry>();
    let mut blocks_size = 0usize;
    for b in blocks {
        blocks_size += b.data.len();
    }

    let block_index_off = (hdr_size + g_size + ent_size) as u64;
    let block_data_off = block_index_off + block_index_size as u64;

    let header = AniHeaderV4 {
        magic: ANI_MAGIC,
        version: ANI_VERSION,
        n_entries: n as u64,
        mph_m: mph.m as u64,
        mph_salt: mph.salt,
        off_mph_g: hdr_size as u64,
        off_entries: (hdr_size + g_size) as u64,
        off_strings: block_index_off,
        off_block_index: block_index_off,
        n_blocks: blocks.len() as u64,
        block_size: block_size as u32,
        _pad: 0,
    };

    let write_start = std::time::Instant::now();
    let file = File::create(output)?;
    let mut file = BufWriter::with_capacity(8 * 1024 * 1024, file);

    let hdr_bytes: &[u8] =
        unsafe { std::slice::from_raw_parts(&header as *const _ as *const u8, hdr_size) };
    file.write_all(hdr_bytes)?;

    let g_bytes: &[u8] = unsafe { std::slice::from_raw_parts(mph.g.as_ptr() as *const u8, g_size) };
    file.write_all(g_bytes)?;

    let entries_bytes: &[u8] =
        unsafe { std::slice::from_raw_parts(entries.as_ptr() as *const u8, ent_size) };
    file.write_all(entries_bytes)?;

    let mut block_entries = Vec::with_capacity(blocks.len());
    let mut cur_off = block_data_off;
    for b in blocks {
        block_entries.push(AniBlockEntry {
            raw_start: b.raw_start,
            raw_len: b.raw_len,
            data_len: b.data.len() as u32,
            data_off: cur_off,
        });
        cur_off += b.data.len() as u64;
    }

    let block_entries_bytes: &[u8] = unsafe {
        std::slice::from_raw_parts(
            block_entries.as_ptr() as *const u8,
            block_entries.len() * mem::size_of::<AniBlockEntry>(),
        )
    };
    file.write_all(block_entries_bytes)?;

    for b in blocks {
        file.write_all(&b.data)?;
    }

    file.flush()?;

    if timing {
        eprintln!(
            "[ani-build] Write to disk: {:.3}s",
            write_start.elapsed().as_secs_f64()
        );
        eprintln!(
            "[ani-build] Total ANI size: {} bytes",
            hdr_size + g_size + ent_size + block_index_size + blocks_size
        );
    }

    Ok(())
}
