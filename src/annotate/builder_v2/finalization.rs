use anyhow::Result;
use bytemuck;
use kira_kv_engine::{BackendKind, Index, IndexBuilder, IndexConfig};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem;
use std::path::Path;

use crate::annotate::builder_v2::StringPool;
use crate::annotate::structs::ani::{
    ANI_MAGIC, ANI_STR_NONE, ANI_VERSION, AniBlockEntry, AniEntry, AniHeaderV6, AniInfoBlobHeader,
    AniInfoPair, AniPosBlock, AniPosContig, AniPosIndexHeader,
};
use crate::util::url_decode_info_value;

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

    let mut config = IndexConfig::default();
    config.auto_detect_numeric = false;
    config.backend = BackendKind::CHD;

    if timing {
        print_mph_config(&config, n, &rows);
    }

    let keys_bytes: Vec<[u8; 8]> = rows.iter().map(|(k, _)| k.to_le_bytes()).collect();

    if timing {
        verify_key_uniqueness(&keys_bytes, &rows);
    }

    let index = IndexBuilder::new()
        .with_config(config)
        .build_index(keys_bytes.clone())?;

    if timing {
        eprintln!(
            "[ani-build] Index construction: {:.3}s",
            mph_start.elapsed().as_secs_f64()
        );
    }

    let mut reordered_entries: Vec<Option<AniEntry>> = vec![None; n];
    let mut reordered_keys: Vec<Option<u64>> = vec![None; n];
    for (i, (key, entry)) in rows.into_iter().enumerate() {
        let idx = index.get(&keys_bytes[i])?;
        if idx >= n {
            anyhow::bail!("Index returned out-of-range idx {}", idx);
        }
        if reordered_entries[idx].is_some() {
            anyhow::bail!("Index collision at idx {}", idx);
        }
        reordered_entries[idx] = Some(entry);
        reordered_keys[idx] = Some(key);
    }
    let entries: Vec<AniEntry> = reordered_entries
        .into_iter()
        .map(|v| v.ok_or_else(|| anyhow::anyhow!("Missing entry after reordering")))
        .collect::<Result<Vec<_>>>()?;
    let keys_ordered: Vec<u64> = reordered_keys
        .into_iter()
        .map(|v| v.ok_or_else(|| anyhow::anyhow!("Missing key after reordering")))
        .collect::<Result<Vec<_>>>()?;

    if timing {
        let stats = index.stats();
        eprintln!(
            "[ani-build] Index stats: engine={}, total_keys={}, total_memory={}",
            stats.engine, stats.total_keys, stats.total_memory
        );
        eprintln!("[ani-build] Sanity (first 5 keys):");
        for (i, key_bytes) in keys_bytes.iter().take(5).enumerate() {
            let key = u64::from_le_bytes(*key_bytes);
            let contains = index.contains(key_bytes);
            let lookup = index.get(key_bytes).ok();
            let lookup_u64 = index.lookup_u64(key).ok();
            eprintln!(
                "  [{}] key={:016x} len={} contains={} lookup={:?} lookup_u64={:?}",
                i,
                key,
                key_bytes.len(),
                contains,
                lookup,
                lookup_u64
            );
        }
        verify_index_correctness(&index, &keys_bytes, &keys_ordered);
    }

    let index_bytes = index.serialize()?;
    let pos_index_bytes = build_pos_index(&entries);
    let info_blob_bytes = build_info_blob(&entries, &mut pool)?;
    write_ani_file(
        output,
        &index_bytes,
        &pos_index_bytes,
        &info_blob_bytes,
        &entries,
        &mut pool,
        timing,
    )?;
    pool.cleanup();

    Ok(())
}

fn print_mph_config(config: &IndexConfig, n: usize, rows: &[(u64, AniEntry)]) {
    eprintln!(
        "[ani-build] MPHF config: gamma={:.2}, rehash_limit={}, entries={}",
        config.mph_config.gamma, config.mph_config.rehash_limit, n
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

fn verify_index_correctness(index: &Index, keys_bytes: &[[u8; 8]], keys_ordered: &[u64]) {
    eprintln!("[ani-build] Verifying index lookups (checking 100 random entries):");
    let mut errors = 0;
    let check_count = 100.min(keys_bytes.len());
    let step = keys_bytes.len() / check_count;

    for i in (0..keys_bytes.len()).step_by(step.max(1)).take(check_count) {
        let key_bytes = &keys_bytes[i];
        let key = u64::from_le_bytes(*key_bytes);
        let idx = match index.get(key_bytes) {
            Ok(v) => v,
            Err(err) => {
                errors += 1;
                if errors <= 5 {
                    eprintln!(
                        "  [{}] ERROR: key={:016x} -> index lookup failed: {:?}",
                        i, key, err
                    );
                }
                continue;
            }
        };
        if idx >= keys_ordered.len() || keys_ordered[idx] != key {
            errors += 1;
            if errors <= 5 {
                let got = keys_ordered.get(idx).copied().unwrap_or(0);
                eprintln!(
                    "  [{}] ERROR: key={:016x} -> mph_idx={} -> key={:016x}",
                    i, key, idx, got
                );
            }
        } else if i < 3 {
            eprintln!("  [{}] OK: key={:016x} -> mph_idx={}", i, key, idx);
        }
    }

    if errors > 0 {
        eprintln!(
            "[ani-build] ERROR: Index verification FAILED with {} errors out of {} checks!",
            errors, check_count
        );
        eprintln!("[ani-build] This means the index will NOT work correctly!");
    } else {
        eprintln!(
            "[ani-build] Index verification: All {} checks passed ✓",
            check_count
        );
    }
}

fn write_ani_file(
    output: &Path,
    index_bytes: &[u8],
    pos_index_bytes: &[u8],
    info_blob_bytes: &[u8],
    entries: &[AniEntry],
    pool: &mut StringPool,
    timing: bool,
) -> Result<()> {
    let n = entries.len();
    let hdr_size = mem::size_of::<AniHeaderV6>();
    let index_size = index_bytes.len();
    let ent_size = n * mem::size_of::<AniEntry>();
    let block_size = pool.block_size();
    let blocks = pool.blocks();
    let block_index_size = blocks.len() * mem::size_of::<AniBlockEntry>();
    let mut blocks_size = 0usize;
    for b in blocks {
        blocks_size += b.data.len();
    }

    let block_index_off = (hdr_size + index_size + ent_size) as u64;
    let block_data_off = block_index_off + block_index_size as u64;
    let pos_index_off = block_data_off + blocks_size as u64;
    let blob_off = pos_index_off + pos_index_bytes.len() as u64;

    let header = AniHeaderV6 {
        magic: ANI_MAGIC,
        version: ANI_VERSION,
        n_entries: n as u64,
        index_len: index_size as u64,
        off_index: hdr_size as u64,
        off_entries: (hdr_size + index_size) as u64,
        off_strings: block_index_off,
        off_block_index: block_index_off,
        n_blocks: blocks.len() as u64,
        block_size: block_size as u32,
        _pad: 0,
        off_pos_index: pos_index_off,
        pos_index_len: pos_index_bytes.len() as u64,
        off_blob: blob_off,
        blob_len: info_blob_bytes.len() as u64,
    };

    let write_start = std::time::Instant::now();
    let file = File::create(output)?;
    let mut file = BufWriter::with_capacity(8 * 1024 * 1024, file);

    let hdr_bytes: &[u8] =
        unsafe { std::slice::from_raw_parts(&header as *const _ as *const u8, hdr_size) };
    file.write_all(hdr_bytes)?;

    file.write_all(index_bytes)?;

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

    file.write_all(pos_index_bytes)?;
    file.write_all(info_blob_bytes)?;

    file.flush()?;

    if timing {
        eprintln!(
            "[ani-build] Write to disk: {:.3}s",
            write_start.elapsed().as_secs_f64()
        );
        eprintln!(
            "[ani-build] Total ANI size: {} bytes",
            hdr_size
                + index_size
                + ent_size
                + block_index_size
                + blocks_size
                + pos_index_bytes.len()
                + info_blob_bytes.len()
        );
    }

    Ok(())
}

fn build_pos_index(entries: &[AniEntry]) -> Vec<u8> {
    let mut per_chr: Vec<Vec<(u32, u32)>> = vec![Vec::new(); 256];
    for (idx, entry) in entries.iter().enumerate() {
        per_chr[entry.chr_id as usize].push((entry.pos, idx as u32));
    }

    let mut contigs: Vec<AniPosContig> = Vec::new();
    let mut blocks: Vec<AniPosBlock> = Vec::new();
    let mut pos_offsets: Vec<u32> = Vec::new();
    let mut pos_counts: Vec<u16> = Vec::new();
    let mut entry_indices: Vec<u32> = Vec::new();

    for (chr_id, list) in per_chr.iter_mut().enumerate() {
        if list.is_empty() {
            continue;
        }
        list.sort_by_key(|(pos, _)| *pos);

        let min_pos = list.first().unwrap().0;
        let max_pos = list.last().unwrap().0;
        let block_start = blocks.len() as u32;

        let mut current_base: Option<u32> = None;
        let mut current_block = AniPosBlock {
            base_pos: 0,
            _pad: 0,
            masks: [0u64; 8],
            offsets_start: 0,
            _pad2: 0,
        };

        let mut i = 0usize;
        while i < list.len() {
            let pos = list[i].0;
            let base = (pos / 512) * 512;
            if current_base != Some(base) {
                if current_base.is_some() {
                    blocks.push(current_block);
                }
                current_base = Some(base);
                current_block = AniPosBlock {
                    base_pos: base,
                    _pad: 0,
                    masks: [0u64; 8],
                    offsets_start: pos_offsets.len() as u32,
                    _pad2: 0,
                };
            }

            let mut count = 0u16;
            let start_idx = i;
            while i < list.len() && list[i].0 == pos {
                let entry_idx = list[i].1;
                entry_indices.push(entry_idx);
                count = count.wrapping_add(1);
                i += 1;
            }

            let bit = (pos - base) as usize;
            let word = bit / 64;
            let bit_in_word = bit % 64;
            current_block.masks[word] |= 1u64 << bit_in_word;
            pos_offsets.push((entry_indices.len() as u32) - (count as u32));
            pos_counts.push(count);

            if i == start_idx {
                i += 1;
            }
        }

        if current_base.is_some() {
            blocks.push(current_block);
        }

        let block_count = blocks.len() as u32 - block_start;
        contigs.push(AniPosContig {
            chr_id: chr_id as u16,
            _pad: 0,
            min_pos,
            max_pos,
            block_start,
            block_count,
        });
    }

    let header_size = mem::size_of::<AniPosIndexHeader>();
    let mut out = vec![0u8; header_size];
    let mut offset = header_size;

    offset = align8(&mut out, offset);
    let off_contigs = offset;
    let contig_bytes = bytemuck::cast_slice(&contigs);
    out.extend_from_slice(contig_bytes);
    offset += contig_bytes.len();

    offset = align8(&mut out, offset);
    let off_blocks = offset;
    let block_bytes = bytemuck::cast_slice(&blocks);
    out.extend_from_slice(block_bytes);
    offset += block_bytes.len();

    offset = align8(&mut out, offset);
    let off_pos_offsets = offset;
    let pos_offsets_bytes = bytemuck::cast_slice(&pos_offsets);
    out.extend_from_slice(pos_offsets_bytes);
    offset += pos_offsets_bytes.len();

    offset = align8(&mut out, offset);
    let off_pos_counts = offset;
    let pos_counts_bytes = bytemuck::cast_slice(&pos_counts);
    out.extend_from_slice(pos_counts_bytes);
    offset += pos_counts_bytes.len();

    offset = align8(&mut out, offset);
    let off_entry_indices = offset;
    let entry_indices_bytes = bytemuck::cast_slice(&entry_indices);
    out.extend_from_slice(entry_indices_bytes);

    let header = AniPosIndexHeader {
        contig_count: contigs.len() as u32,
        block_count: blocks.len() as u32,
        pos_count: pos_offsets.len() as u32,
        entry_index_count: entry_indices.len() as u32,
        off_contigs: off_contigs as u32,
        off_blocks: off_blocks as u32,
        off_pos_offsets: off_pos_offsets as u32,
        off_pos_counts: off_pos_counts as u32,
        off_entry_indices: off_entry_indices as u32,
    };

    let hdr_bytes: &[u8] =
        unsafe { std::slice::from_raw_parts(&header as *const _ as *const u8, header_size) };
    out[..header_size].copy_from_slice(hdr_bytes);
    out
}

fn build_info_blob(entries: &[AniEntry], pool: &mut StringPool) -> Result<Vec<u8>> {
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
    let raw = pool.materialize()?;
    let mut dict_map: std::collections::HashMap<String, u32> = std::collections::HashMap::new();
    let mut dict_list: Vec<String> = Vec::new();
    let mut entry_offsets: Vec<u32> = Vec::with_capacity(entries.len());
    let mut entry_counts: Vec<u16> = Vec::with_capacity(entries.len());
    let mut pairs: Vec<AniInfoPair> = Vec::new();
    let mut values: Vec<u8> = Vec::new();

    for entry in entries {
        entry_offsets.push(pairs.len() as u32);
        let mut count: u16 = 0;
        if entry.info_ofs == ANI_STR_NONE || entry.info_len == 0 {
            entry_counts.push(0);
            continue;
        }
        let ofs = entry.info_ofs as usize;
        let len = entry.info_len as usize;
        if ofs + len > raw.len() {
            entry_counts.push(0);
            continue;
        }
        let info = &raw[ofs..ofs + len];
        let info_str = std::str::from_utf8(info).unwrap_or("");
        let decoded_info = url_decode_info_value(info_str);
        if decoded_info.is_empty() || decoded_info == "." {
            entry_counts.push(0);
            continue;
        }
        for token in decoded_info.split(';') {
            if token.is_empty() {
                continue;
            }
            let (key, value_opt) = if let Some((k, v)) = token.split_once('=') {
                (k, Some(v))
            } else {
                (token, None)
            };
            if key.is_empty() {
                continue;
            }
            let tag_id = *dict_map.entry(key.to_string()).or_insert_with(|| {
                let id = dict_list.len() as u32;
                dict_list.push(key.to_string());
                id
            });
            if let Some(value) = value_opt {
                let value_off = values.len() as u32;
                let value_len = value.len() as u32;
                values.extend_from_slice(value.as_bytes());
                pairs.push(AniInfoPair {
                    tag_id,
                    value_off,
                    value_len,
                });
            } else {
                pairs.push(AniInfoPair {
                    tag_id,
                    value_off: 0,
                    value_len: 0,
                });
            }
            count = count.wrapping_add(1);
        }
        entry_counts.push(count);
    }

    if debug {
        eprintln!(
            "[ANI-BLOB] entries={}, pairs={}, dict={}, values={}",
            entries.len(),
            pairs.len(),
            dict_list.len(),
            values.len()
        );
    }

    let mut dict_offsets: Vec<u32> = Vec::with_capacity(dict_list.len());
    let mut dict_data: Vec<u8> = Vec::new();
    for key in &dict_list {
        dict_offsets.push(dict_data.len() as u32);
        dict_data.extend_from_slice(key.as_bytes());
        dict_data.push(0);
    }

    let header_size = mem::size_of::<AniInfoBlobHeader>();
    let mut out = vec![0u8; header_size];
    let mut offset = header_size;

    offset = align8(&mut out, offset);
    let off_dict_offsets = offset;
    let dict_offsets_bytes = bytemuck::cast_slice(&dict_offsets);
    out.extend_from_slice(dict_offsets_bytes);
    offset += dict_offsets_bytes.len();

    offset = align8(&mut out, offset);
    let off_dict_data = offset;
    out.extend_from_slice(&dict_data);
    offset += dict_data.len();

    offset = align8(&mut out, offset);
    let off_entry_offsets = offset;
    let entry_offsets_bytes = bytemuck::cast_slice(&entry_offsets);
    out.extend_from_slice(entry_offsets_bytes);
    offset += entry_offsets_bytes.len();

    offset = align8(&mut out, offset);
    let off_entry_counts = offset;
    let entry_counts_bytes = bytemuck::cast_slice(&entry_counts);
    out.extend_from_slice(entry_counts_bytes);
    offset += entry_counts_bytes.len();

    offset = align8(&mut out, offset);
    let off_pairs = offset;
    let pairs_bytes = bytemuck::cast_slice(&pairs);
    out.extend_from_slice(pairs_bytes);
    offset += pairs_bytes.len();

    offset = align8(&mut out, offset);
    let off_values = offset;
    out.extend_from_slice(&values);

    let header = AniInfoBlobHeader {
        n_entries: entries.len() as u64,
        dict_count: dict_list.len() as u32,
        pair_count: pairs.len() as u32,
        off_dict_offsets: off_dict_offsets as u64,
        off_dict_data: off_dict_data as u64,
        off_entry_offsets: off_entry_offsets as u64,
        off_entry_counts: off_entry_counts as u64,
        off_pairs: off_pairs as u64,
        off_values: off_values as u64,
    };

    let hdr_bytes: &[u8] =
        unsafe { std::slice::from_raw_parts(&header as *const _ as *const u8, header_size) };
    out[..header_size].copy_from_slice(hdr_bytes);
    Ok(out)
}

fn align8(out: &mut Vec<u8>, mut offset: usize) -> usize {
    let pad = (8 - (offset % 8)) % 8;
    for _ in 0..pad {
        out.push(0u8);
    }
    offset += pad;
    offset
}
