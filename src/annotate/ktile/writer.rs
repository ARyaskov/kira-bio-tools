//! `.ktile` writer — streaming.
//!
//! Reads the VCF/BGZF input once and streams the line-pool bytes straight
//! to the output as they arrive. File layout matches [`super::format`]:
//! header (finalised via `seek(0)` after data is on disk), headers_blob,
//! line_pool, then the typed-column tables.

use std::fs::OpenOptions;
use std::io::{Seek, SeekFrom, Write};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::{Compress, Compression, FlushCompress, Status};

use super::format::{
    CompressedChunkEntry, KTILE_MAGIC, KTILE_VERSION, KtileHeader, flags as ktile_flags,
};
use crate::util::chr_name_to_id;
use crate::vcf::UnifiedVcfReader;

/// Build-time options for `.ktile` writers.
#[derive(Debug, Clone, Copy)]
pub struct KtileWriteOptions {
    /// When true, `line_pool` is stored as raw-deflate-compressed
    /// fixed-line chunks.
    pub compressed: bool,
    /// Lines per compressed chunk.
    pub lines_per_chunk: u32,
    /// Deflate level (0-9).
    pub deflate_level: u32,
}

impl Default for KtileWriteOptions {
    fn default() -> Self {
        Self {
            compressed: true,
            lines_per_chunk: 4096,
            deflate_level: 1,
        }
    }
}

/// Build a `.ktile` sidecar from a VCF/BGZF input.
pub fn write_ktile_from_vcf<I: AsRef<Path>, O: AsRef<Path>>(
    input: I,
    output: O,
) -> Result<KtileBuildStats> {
    write_ktile_from_vcf_with(input, output, KtileWriteOptions::default())
}

/// Variant of [`write_ktile_from_vcf`] with explicit build options.
pub fn write_ktile_from_vcf_with<I: AsRef<Path>, O: AsRef<Path>>(
    input: I,
    output: O,
    opts: KtileWriteOptions,
) -> Result<KtileBuildStats> {
    let input_ref = input.as_ref();
    let (source_size, source_mtime_unix) = read_source_metadata(input_ref);

    let mut reader = UnifiedVcfReader::open(input_ref)
        .with_context(|| format!("opening input {:?}", input_ref))?;
    let headers = reader
        .header()
        .context("reading VCF header block")?;

    // Open output + write placeholder header (finalised after data is on disk).
    let out_file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(output.as_ref())
        .with_context(|| format!("creating output {:?}", output.as_ref()))?;
    let mut out = std::io::BufWriter::with_capacity(16 * 1024 * 1024, out_file);

    let mut cursor: u64 = 0;

    let placeholder = vec![0u8; KtileHeader::SIZE];
    out.write_all(&placeholder)?;
    cursor += KtileHeader::SIZE as u64;

    let headers_blob: String = headers.join("\n");
    let headers_bytes = headers_blob.as_bytes();
    let headers_off = cursor;
    out.write_all(headers_bytes)?;
    cursor += headers_bytes.len() as u64;

    // Pad so the following line_offsets section is 8-aligned.
    cursor = pad_with_zeros(&mut out, cursor, 8)?;

    // Stream line_pool + collect typed columns in RAM.
    let line_pool_off = cursor;
    let mut line_offsets: Vec<u64> = Vec::with_capacity(8_192);
    let mut chr_ids: Vec<u32> = Vec::with_capacity(8_192);
    let mut positions: Vec<u32> = Vec::with_capacity(8_192);
    let mut ref_offsets: Vec<u32> = Vec::with_capacity(8_192);
    let mut ref_lens: Vec<u32> = Vec::with_capacity(8_192);
    let mut alt_offsets: Vec<u32> = Vec::with_capacity(8_192);
    let mut alt_lens: Vec<u32> = Vec::with_capacity(8_192);
    line_offsets.push(0);
    let mut pool_cursor: u64 = 0;

    let mut chunk_uncompressed_buf: Vec<u8> = if opts.compressed {
        Vec::with_capacity(32 * 1024 * 1024)
    } else {
        Vec::new()
    };
    let mut chunk_scratch: Vec<u8> = Vec::new();
    let mut compressor = Compress::new(Compression::new(opts.deflate_level), false);
    let mut chunks: Vec<CompressedChunkEntry> = Vec::new();
    let mut chunk_uncompressed_off: u64 = 0;
    let mut chunk_compressed_off: u64 = 0;
    let mut line_in_chunk: u32 = 0;
    let mut compressed_pool_written: u64 = 0;

    while let Some(line) = reader.read_line().context("reading next VCF line")? {
        if line.is_empty() {
            continue;
        }
        let bytes = line.as_bytes();
        let mut tabs = [0u32; 5];
        let mut tab_count = 0usize;
        for (i, b) in bytes.iter().enumerate() {
            if *b == b'\t' {
                tabs[tab_count] = i as u32;
                tab_count += 1;
                if tab_count == 5 {
                    break;
                }
            }
        }
        if tab_count < 5 {
            continue;
        }
        let chrom = unsafe { std::str::from_utf8_unchecked(&bytes[..tabs[0] as usize]) };
        let pos_str = unsafe {
            std::str::from_utf8_unchecked(&bytes[(tabs[0] as usize + 1)..tabs[1] as usize])
        };
        let Ok(pos) = pos_str.parse::<u32>() else {
            continue;
        };
        let chr_id = chr_name_to_id(chrom).map(u32::from).unwrap_or(u32::MAX);
        let ref_off = tabs[2] + 1;
        let ref_end = tabs[3];
        let alt_off = tabs[3] + 1;
        let alt_end = tabs[4];

        chr_ids.push(chr_id);
        positions.push(pos);
        ref_offsets.push(ref_off);
        ref_lens.push(ref_end - ref_off);
        alt_offsets.push(alt_off);
        alt_lens.push(alt_end - alt_off);

        if opts.compressed {
            chunk_uncompressed_buf.extend_from_slice(bytes);
            line_in_chunk += 1;
            if line_in_chunk >= opts.lines_per_chunk {
                let written = flush_compressed_chunk(
                    &mut out,
                    &mut compressor,
                    &chunk_uncompressed_buf,
                    &mut chunk_scratch,
                )?;
                chunks.push(CompressedChunkEntry {
                    compressed_off: chunk_compressed_off,
                    compressed_size: written as u32,
                    uncompressed_off: chunk_uncompressed_off,
                    uncompressed_size: chunk_uncompressed_buf.len() as u32,
                });
                chunk_uncompressed_off += chunk_uncompressed_buf.len() as u64;
                chunk_compressed_off += written as u64;
                compressed_pool_written += written as u64;
                chunk_uncompressed_buf.clear();
                line_in_chunk = 0;
            }
        } else {
            out.write_all(bytes)?;
        }
        pool_cursor += line.len() as u64;
        line_offsets.push(pool_cursor);
    }

    if opts.compressed && !chunk_uncompressed_buf.is_empty() {
        let written = flush_compressed_chunk(
            &mut out,
            &mut compressor,
            &chunk_uncompressed_buf,
            &mut chunk_scratch,
        )?;
        chunks.push(CompressedChunkEntry {
            compressed_off: chunk_compressed_off,
            compressed_size: written as u32,
            uncompressed_off: chunk_uncompressed_off,
            uncompressed_size: chunk_uncompressed_buf.len() as u32,
        });
        compressed_pool_written += written as u64;
        chunk_uncompressed_buf.clear();
    }

    let n_records = chr_ids.len();
    let line_pool_len = if opts.compressed {
        compressed_pool_written
    } else {
        pool_cursor
    };
    cursor += line_pool_len;

    cursor = pad_with_zeros(&mut out, cursor, 8)?;
    let line_offsets_off = cursor;
    out.write_all(bytemuck::cast_slice(&line_offsets))?;
    cursor += (line_offsets.len() * std::mem::size_of::<u64>()) as u64;

    let chr_ids_off = cursor;
    out.write_all(bytemuck::cast_slice(&chr_ids))?;
    cursor += (chr_ids.len() * std::mem::size_of::<u32>()) as u64;
    cursor = pad_with_zeros(&mut out, cursor, 4)?;
    let positions_off = cursor;
    out.write_all(bytemuck::cast_slice(&positions))?;
    cursor += (positions.len() * std::mem::size_of::<u32>()) as u64;

    let off_ref_offsets = cursor;
    out.write_all(bytemuck::cast_slice(&ref_offsets))?;
    cursor += (ref_offsets.len() * std::mem::size_of::<u32>()) as u64;
    let off_ref_lens = cursor;
    out.write_all(bytemuck::cast_slice(&ref_lens))?;
    cursor += (ref_lens.len() * std::mem::size_of::<u32>()) as u64;
    let off_alt_offsets = cursor;
    out.write_all(bytemuck::cast_slice(&alt_offsets))?;
    cursor += (alt_offsets.len() * std::mem::size_of::<u32>()) as u64;
    let off_alt_lens = cursor;
    out.write_all(bytemuck::cast_slice(&alt_lens))?;
    cursor += (alt_lens.len() * std::mem::size_of::<u32>()) as u64;

    let (off_chunk_index, n_chunks) = if opts.compressed {
        cursor = pad_with_zeros(&mut out, cursor, 8)?;
        let off = cursor;
        out.write_all(bytemuck::cast_slice(&chunks))?;
        cursor += (chunks.len() * std::mem::size_of::<CompressedChunkEntry>()) as u64;
        (off, chunks.len() as u32)
    } else {
        (0, 0)
    };

    let mut flag_bits = ktile_flags::HAS_REF_ALT_COLUMNS;
    if opts.compressed {
        flag_bits |= ktile_flags::HAS_COMPRESSED_POOL;
    }
    let header = KtileHeader {
        magic: KTILE_MAGIC,
        version: KTILE_VERSION,
        flags: flag_bits,
        n_records: n_records as u64,
        headers_off,
        headers_len: headers_bytes.len() as u64,
        line_offsets_off,
        chr_ids_off,
        positions_off,
        line_pool_off,
        line_pool_len,
        source_size,
        source_mtime_unix,
        off_ref_offsets,
        off_ref_lens,
        off_alt_offsets,
        off_alt_lens,
        lines_per_chunk: if opts.compressed { opts.lines_per_chunk } else { 0 },
        n_chunks,
        off_chunk_index,
    };
    out.flush()?;
    let mut file = out.into_inner().map_err(|e| {
        let io_err: std::io::Error = e.into_error();
        anyhow::anyhow!("BufWriter into_inner: {io_err}")
    })?;
    file.seek(SeekFrom::Start(0))?;
    file.write_all(bytemuck::bytes_of(&header))?;
    file.flush()?;

    Ok(KtileBuildStats {
        n_records,
        bytes_written: cursor,
    })
}

/// Raw-deflate compresses `data` as an independent stream, writes it to
/// `out`, returns the number of bytes written.
fn flush_compressed_chunk<W: Write>(
    out: &mut W,
    compressor: &mut Compress,
    data: &[u8],
    scratch: &mut Vec<u8>,
) -> Result<usize> {
    let bound = data.len() + (data.len() >> 12) + (data.len() >> 14) + (data.len() >> 25) + 13;
    if scratch.len() < bound {
        scratch.resize(bound, 0);
    }
    compressor.reset();
    let before_out = compressor.total_out();
    let before_in = compressor.total_in();
    let status = compressor
        .compress(data, scratch, FlushCompress::Finish)
        .map_err(|e| anyhow::anyhow!("flate2 compress failed: {e}"))?;
    if status != Status::StreamEnd || (compressor.total_in() - before_in) as usize != data.len() {
        anyhow::bail!(
            "ktile chunk compress: scratch undersized or input not consumed"
        );
    }
    let n = (compressor.total_out() - before_out) as usize;
    out.write_all(&scratch[..n])?;
    Ok(n)
}

/// Pads `out` with zero bytes up to the next multiple of `align`.
fn pad_with_zeros<W: Write>(
    out: &mut W,
    cursor: u64,
    align: u64,
) -> Result<u64> {
    debug_assert!(align.is_power_of_two());
    let aligned = (cursor + align - 1) & !(align - 1);
    let pad = (aligned - cursor) as usize;
    if pad > 0 {
        let zeros = [0u8; 8];
        out.write_all(&zeros[..pad])?;
    }
    Ok(aligned)
}

#[derive(Debug, Clone, Copy)]
pub struct KtileBuildStats {
    pub n_records: usize,
    pub bytes_written: u64,
}

/// Reads `path`'s byte size + mtime (Unix seconds), returning `(0, 0)` if
/// the file can't be stat'd.
fn read_source_metadata(path: &Path) -> (u64, u64) {
    let Ok(meta) = std::fs::metadata(path) else {
        return (0, 0);
    };
    let size = meta.len();
    let mtime = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_secs())
        .unwrap_or(0);
    (size, mtime)
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_ktile_writer.rs"]
mod tests;
