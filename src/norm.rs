use anyhow::Result;
use memmap2::{Mmap, MmapOptions};
use rayon::prelude::*;
use std::cell::RefCell;
use std::fs::{File, OpenOptions};
use std::path::Path;
use thread_local::ThreadLocal;

#[inline]
pub fn normalize(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let rb = ref_allele.as_bytes();
    let ab = alt_allele.as_bytes();

    let mut prefix = 0;
    while prefix < rb.len() && prefix < ab.len() && rb[prefix] == ab[prefix] {
        prefix += 1;
    }

    let mut suffix = 0;
    let mut ri = rb.len();
    let mut ai = ab.len();
    while ri > prefix && ai > prefix && rb[ri - 1] == ab[ai - 1] {
        ri -= 1;
        ai -= 1;
        suffix += 1;
    }

    (prefix, suffix)
}

#[inline]
fn extract_tabs(line: &[u8]) -> Option<(usize, usize, usize, usize)> {
    let mut t = [0usize; 5];
    let mut k = 0;

    for i in 0..line.len() {
        if line[i] == b'\t' {
            t[k] = i;
            k += 1;
            if k == 5 {
                return Some((t[2] + 1, t[3], t[3] + 1, t[4]));
            }
        }
    }
    None
}

pub fn turbo_norm_vcf(input: &Path, output: &Path) -> Result<()> {
    let file = File::open(input)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data = &mmap[..];

    let mut line_positions = Vec::with_capacity(data.len() / 80);
    let mut pos = 0;

    while let Some(off) = memchr::memchr(b'\n', &data[pos..]) {
        line_positions.push(pos + off);
        pos += off + 1;
    }

    let arenas: ThreadLocal<RefCell<Vec<u8>>> = ThreadLocal::new();

    line_positions.par_iter().enumerate().for_each(|(i, &end)| {
        let arena_cell = arenas.get_or(|| RefCell::new(Vec::with_capacity(1024 * 1024)));
        let mut arena = arena_cell.borrow_mut();

        let start = if i == 0 { 0 } else { line_positions[i - 1] + 1 };
        let line = &data[start..end];

        if line.is_empty() || line[0] == b'#' {
            arena.extend_from_slice(line);
            arena.push(b'\n');
            return;
        }

        let Some((rs, re, as_, ae)) = extract_tabs(line) else {
            arena.extend_from_slice(line);
            arena.push(b'\n');
            return;
        };

        let ref_allele = &line[rs..re];
        let alt_allele = &line[as_..ae];

        let ref_str = unsafe { std::str::from_utf8_unchecked(ref_allele) };
        let alt_str = unsafe { std::str::from_utf8_unchecked(alt_allele) };

        let (prefix, suffix) = normalize(ref_str, alt_str);

        let nr = &ref_allele[prefix..ref_allele.len() - suffix];
        let na = &alt_allele[prefix..alt_allele.len() - suffix];

        arena.extend_from_slice(&line[..rs]);
        arena.extend_from_slice(nr);
        arena.push(b'\t');
        arena.extend_from_slice(na);
        arena.extend_from_slice(&line[alt_allele.len() + rs..]);
        arena.push(b'\n');
    });

    let mut final_buf = Vec::with_capacity(data.len() + data.len() / 3);
    arenas
        .into_iter()
        .for_each(|buf_cell| final_buf.extend_from_slice(&buf_cell.into_inner()));

    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    {
        let mut writer = BufWriter::with_capacity(64 * 1024 * 1024, File::create(output)?);
        writer.write_all(&final_buf)?;
        writer.flush()?;
        return Ok(());
    }

    #[cfg(not(all(target_os = "macos", target_arch = "aarch64")))]
    {
        let out_file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(output)?;
        out_file.set_len(final_buf.len() as u64)?;
        let mut mmap_out = unsafe { MmapOptions::new().map_mut(&out_file)? };
        unsafe {
            std::ptr::copy_nonoverlapping(
                final_buf.as_ptr(),
                mmap_out.as_mut_ptr(),
                final_buf.len(),
            );
        }
        mmap_out.flush()?;
        return Ok(());
    }
}
