//! Pre-scan: walk reads vs FASTA once, producing both
//! - `InterestingMap`: positions with ≥2 reads of non-ref evidence (engine skip-list)
//! - `needs_baq[read_idx]`: per-read flag — `false` means the read is fully ref-match
//!   with no indels in CIGAR, so BAQ HMM is a no-op for it and can be skipped.

use fxhash::FxHashMap;
use noodles_sam::alignment::record::cigar::op::Kind;
use rayon::prelude::*;

use crate::bam::pileup::LiveRead;
use crate::bam::reader::FastaLike;

pub struct InterestingMap {
    per_ref: FxHashMap<usize, Vec<u32>>,
}

impl InterestingMap {
    #[inline]
    pub fn next_at_or_after(&self, ref_id: usize, pos: u32) -> Option<u32> {
        let v = self.per_ref.get(&ref_id)?;
        let i = v.partition_point(|&p| p < pos);
        v.get(i).copied()
    }

    pub fn total(&self) -> usize {
        self.per_ref.values().map(|v| v.len()).sum()
    }

    #[inline]
    pub fn contains(&self, ref_id: usize, pos: u32) -> bool {
        self.per_ref.get(&ref_id).is_some_and(|v| v.binary_search(&pos).is_ok())
    }

    /// Union with another map (positions from all samples are interesting).
    pub fn merge(&mut self, other: InterestingMap) {
        for (rid, mut v) in other.per_ref {
            let e = self.per_ref.entry(rid).or_default();
            e.append(&mut v);
            e.sort_unstable();
            e.dedup();
        }
    }
}

pub struct PreScan {
    pub pos_filter: InterestingMap,
    pub needs_baq: Vec<bool>,
    pub skipped_baq: usize,
}

/// Walk all reads once vs FASTA. Returns positions with ≥`min_alt_reads` non-ref
/// supporting reads (strict filter — drops single-read errors) AND a per-read
/// flag of whether BAQ is needed (mismatch or indel present).
pub fn pre_scan(
    records_per_sample: &[std::sync::Arc<Vec<LiveRead>>],
    fasta: &(impl FastaLike + Sync),
    ref_names: &[String],
    min_alt_reads: u32,
) -> PreScan {
    use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

    // Flatten records, but keep stable per-read indexing per sample so the
    // returned needs_baq aligns with the original Vec<LiveRead> slot.
    let mut total_reads = 0usize;
    let mut sample_offsets: Vec<usize> = Vec::with_capacity(records_per_sample.len() + 1);
    for sample in records_per_sample {
        sample_offsets.push(total_reads);
        total_reads += sample.len();
    }
    sample_offsets.push(total_reads);

    if total_reads == 0 {
        return PreScan {
            pos_filter: InterestingMap { per_ref: FxHashMap::default() },
            needs_baq: Vec::new(),
            skipped_baq: 0,
        };
    }

    let no_pb = std::env::var("KIRA_BT_NO_PROGRESS").is_ok();
    let pb = ProgressBar::new(total_reads as u64);
    pb.set_style(
        ProgressStyle::with_template("[SCAN] {bar:40.green/blue} {pos}/{len} reads ({per_sec}, ETA {eta})")
            .unwrap_or_else(|_| ProgressStyle::default_bar()),
    );
    pb.set_draw_target(if no_pb { ProgressDrawTarget::hidden() } else { ProgressDrawTarget::stderr_with_hz(2) });
    let counter = Arc::new(AtomicU64::new(0));
    let stop = Arc::new(AtomicBool::new(false));
    let ticker = {
        let counter = Arc::clone(&counter);
        let stop = Arc::clone(&stop);
        let pb_clone = pb.clone();
        std::thread::spawn(move || {
            while !stop.load(Ordering::Relaxed) {
                pb_clone.set_position(counter.load(Ordering::Relaxed));
                std::thread::sleep(std::time::Duration::from_millis(500));
            }
            pb_clone.set_position(counter.load(Ordering::Relaxed));
        })
    };

    // Per-thread partial results: (positions per ref, per-read needs_baq for a range of indices).
    struct Partial {
        positions: FxHashMap<usize, Vec<u32>>,
        needs_baq_range: Vec<(usize, bool)>, // (global_read_idx, needs_baq)
    }

    let chunk_size = (total_reads / rayon::current_num_threads()).max(1024);

    // Build chunk specs: (sample_idx, local_lo, local_hi, global_offset)
    let mut chunks: Vec<(usize, usize, usize, usize)> = Vec::new();
    for (si, sample) in records_per_sample.iter().enumerate() {
        let base = sample_offsets[si];
        let mut lo = 0usize;
        while lo < sample.len() {
            let hi = (lo + chunk_size).min(sample.len());
            chunks.push((si, lo, hi, base + lo));
            lo = hi;
        }
    }

    let partials: Vec<Partial> = chunks
        .par_iter()
        .map(|&(si, lo, hi, gbase)| {
            let sample = &records_per_sample[si];
            let mut positions: FxHashMap<usize, Vec<u32>> = FxHashMap::default();
            let mut needs_baq_range: Vec<(usize, bool)> = Vec::with_capacity(hi - lo);
            for (k, lr) in sample[lo..hi].iter().enumerate() {
                let nb = scan_read(lr, fasta, ref_names, &mut positions);
                needs_baq_range.push((gbase + k, nb));
            }
            counter.fetch_add((hi - lo) as u64, Ordering::Relaxed);
            Partial { positions, needs_baq_range }
        })
        .collect();

    stop.store(true, Ordering::Relaxed);
    let _ = ticker.join();
    pb.finish_with_message("SCAN done");

    // Merge needs_baq
    let mut needs_baq = vec![false; total_reads];
    let mut all_positions: FxHashMap<usize, Vec<u32>> = FxHashMap::default();
    for p in partials {
        for (idx, nb) in p.needs_baq_range {
            needs_baq[idx] = nb;
        }
        for (rid, mut v) in p.positions {
            all_positions.entry(rid).or_default().append(&mut v);
        }
    }

    // Strict filter: dedup + count, keep positions with ≥ min_alt_reads.
    all_positions.par_iter_mut().for_each(|(_, v)| {
        v.sort_unstable();
        if min_alt_reads <= 1 {
            v.dedup();
            return;
        }
        let mut filtered: Vec<u32> = Vec::with_capacity(v.len() / 2);
        let mut i = 0;
        while i < v.len() {
            let p = v[i];
            let mut j = i + 1;
            while j < v.len() && v[j] == p { j += 1; }
            if (j - i) as u32 >= min_alt_reads {
                filtered.push(p);
            }
            i = j;
        }
        *v = filtered;
    });

    let skipped_baq = needs_baq.iter().filter(|&&b| !b).count();
    PreScan {
        pos_filter: InterestingMap { per_ref: all_positions },
        needs_baq,
        skipped_baq,
    }
}

/// Returns true if the read has any mismatch (in Match ops) OR any indel
/// (Insertion/Deletion) — i.e. anything BAQ would change. Pushes mismatch
/// and indel positions into `out` (for skip-list construction).
fn scan_read(
    lr: &LiveRead,
    fasta: &impl FastaLike,
    ref_names: &[String],
    out: &mut FxHashMap<usize, Vec<u32>>,
) -> bool {
    let Some(chr) = ref_names.get(lr.ref_id) else { return false };
    let mut span: u32 = 0;
    let mut has_indel = false;
    for &(k, l) in &lr.cigar_pairs {
        match k {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip => {
                span += l;
            }
            _ => {}
        }
        if matches!(k, Kind::Insertion | Kind::Deletion) { has_indel = true; }
    }
    if span == 0 { return has_indel; }
    let Some(ref_slice) = fasta.slice(chr, lr.ref_start + 1, span as usize) else {
        return has_indel;
    };

    let entry = out.entry(lr.ref_id).or_default();
    let mut r_off: usize = 0;
    let mut q_off: usize = 0;
    let mut has_mismatch = false;
    for &(kind, len) in &lr.cigar_pairs {
        let l = len as usize;
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..l {
                    let qb = lr.seq().get(q_off + i).copied().unwrap_or(b'N').to_ascii_uppercase();
                    let rb = ref_slice.get(r_off + i).copied().unwrap_or(b'N').to_ascii_uppercase();
                    if rb != b'N' && qb != b'N' && qb != rb {
                        entry.push(lr.ref_start + (r_off + i) as u32);
                        has_mismatch = true;
                    }
                }
                r_off += l;
                q_off += l;
            }
            Kind::Insertion | Kind::SoftClip => {
                if matches!(kind, Kind::Insertion) && r_off > 0 {
                    // Anchor = last match base before the insertion.
                    entry.push(lr.ref_start + (r_off - 1) as u32);
                }
                q_off += l;
            }
            Kind::Deletion | Kind::Skip => {
                if matches!(kind, Kind::Deletion) && r_off > 0 {
                    // Anchor = last match base before the deletion. indel_after()
                    // returns the deletion event when cur_pos == anchor.
                    entry.push(lr.ref_start + (r_off - 1) as u32);
                }
                r_off += l;
            }
            _ => {}
        }
    }
    has_mismatch || has_indel
}

// Legacy entry point — kept for compatibility, delegates to pre_scan with min=2.
pub fn build(
    records_per_sample: &[std::sync::Arc<Vec<LiveRead>>],
    fasta: &(impl FastaLike + Sync),
    ref_names: &[String],
) -> InterestingMap {
    pre_scan(records_per_sample, fasta, ref_names, 2).pos_filter
}

const RECAL_Q: usize = 64;
const RECAL_CTX: usize = 64;
const RECAL_CYCLES: usize = 8;
/// Cells with fewer observations keep the reported quality.
const RECAL_MIN_OBS: u32 = 200;

/// Empirical base-quality table (`mpileup --recal`): (reported Q, reference
/// trinucleotide, strand, cycle bin) → observed bases and mismatches,
/// collected at sites the pre-scan did not flag as candidate variants. A
/// first-pass recalibration without a known-sites database.
pub struct RecalTable {
    n: Vec<u32>,
    mism: Vec<u32>,
}

#[inline]
fn recal_cell(q: u8, ctx: usize, rev: bool, cycle: usize) -> usize {
    (((q.min(63) as usize) * RECAL_CTX + ctx) * 2 + rev as usize) * RECAL_CYCLES + cycle
}

/// Trinucleotide code of the reference around `i`, or `None` near N/edges.
#[inline]
fn tri_ctx(r: &[u8], i: usize) -> Option<usize> {
    if i == 0 || i + 1 >= r.len() {
        return None;
    }
    let code = |b: u8| match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    };
    Some(code(r[i - 1])? * 16 + code(r[i])? * 4 + code(r[i + 1])?)
}

/// Visit every aligned base of `lr` as (query index, reference offset in `ref_slice`).
fn for_each_aligned(lr: &LiveRead, mut f: impl FnMut(usize, usize)) {
    let (mut r_off, mut q_off) = (0usize, 0usize);
    for &(kind, len) in &lr.cigar_pairs {
        let l = len as usize;
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..l {
                    f(q_off + i, r_off + i);
                }
                r_off += l;
                q_off += l;
            }
            Kind::Insertion | Kind::SoftClip => q_off += l,
            Kind::Deletion | Kind::Skip => r_off += l,
            _ => {}
        }
    }
}

fn read_ref_span(lr: &LiveRead) -> u32 {
    lr.cigar_pairs
        .iter()
        .filter(|(k, _)| matches!(k, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip))
        .map(|(_, l)| *l)
        .sum()
}

impl RecalTable {
    fn new() -> Self {
        let cells = RECAL_Q * RECAL_CTX * 2 * RECAL_CYCLES;
        Self { n: vec![0; cells], mism: vec![0; cells] }
    }

    fn merge(&mut self, other: &RecalTable) {
        for (a, b) in self.n.iter_mut().zip(&other.n) {
            *a += b;
        }
        for (a, b) in self.mism.iter_mut().zip(&other.mism) {
            *a += b;
        }
    }

    /// Collect counts from every read, skipping positions in `candidates`.
    pub fn build(
        records_per_sample: &[std::sync::Arc<Vec<LiveRead>>],
        fasta: &(impl FastaLike + Sync),
        ref_names: &[String],
        candidates: &InterestingMap,
    ) -> Self {
        let chunk = 4096;
        let partials: Vec<RecalTable> = records_per_sample
            .iter()
            .flat_map(|s| s.chunks(chunk))
            .collect::<Vec<_>>()
            .par_iter()
            .map(|reads| {
                let mut t = RecalTable::new();
                for lr in reads.iter() {
                    let Some(chr) = ref_names.get(lr.ref_id) else { continue };
                    let span = read_ref_span(lr);
                    if span == 0 {
                        continue;
                    }
                    let Some(rs) = fasta.slice(chr, lr.ref_start + 1, span as usize) else { continue };
                    let (seq, qual) = (lr.seq(), lr.qual());
                    let len = seq.len().max(1);
                    let rev = lr.is_reverse();
                    for_each_aligned(lr, |qi, ri| {
                        if qi >= seq.len() || ri >= rs.len() {
                            return;
                        }
                        let Some(ctx) = tri_ctx(rs, ri) else { return };
                        if candidates.contains(lr.ref_id, lr.ref_start + ri as u32) {
                            return;
                        }
                        let qb = seq[qi].to_ascii_uppercase();
                        if qb == b'N' {
                            return;
                        }
                        let cycle = if rev { len - 1 - qi } else { qi } * RECAL_CYCLES / len;
                        let cell = recal_cell(qual[qi], ctx, rev, cycle.min(RECAL_CYCLES - 1));
                        t.n[cell] += 1;
                        if qb != rs[ri].to_ascii_uppercase() {
                            t.mism[cell] += 1;
                        }
                    });
                }
                t
            })
            .collect();
        let mut out = RecalTable::new();
        for p in &partials {
            out.merge(p);
        }
        out
    }

    /// Empirical quality for a cell, or `None` when it has too few observations.
    #[inline]
    pub fn empirical(&self, q: u8, ctx: usize, rev: bool, cycle: usize) -> Option<u8> {
        let cell = recal_cell(q, ctx, rev, cycle);
        let n = self.n[cell];
        if n < RECAL_MIN_OBS {
            return None;
        }
        let e = (self.mism[cell] as f64 + 1.0) / (n as f64 + 2.0);
        Some((-10.0 * e.log10()).round().clamp(2.0, 60.0) as u8)
    }

    /// Number of cells with enough observations to recalibrate.
    pub fn n_calibrated(&self) -> usize {
        self.n.iter().filter(|&&n| n >= RECAL_MIN_OBS).count()
    }

    /// Replace reported qualities by the empirical ones in place.
    pub fn apply(&self, reads: &mut [LiveRead], fasta: &(impl FastaLike + Sync), ref_names: &[String]) {
        reads.par_chunks_mut(4096).for_each(|chunk| {
            for lr in chunk.iter_mut() {
                let Some(chr) = ref_names.get(lr.ref_id) else { continue };
                let span = read_ref_span(lr);
                if span == 0 {
                    continue;
                }
                let Some(rs) = fasta.slice(chr, lr.ref_start + 1, span as usize) else { continue };
                let rev = lr.is_reverse();
                let mut updates: Vec<(usize, u8)> = Vec::new();
                {
                    let seq = lr.seq();
                    let qual = lr.qual();
                    let len = seq.len().max(1);
                    for_each_aligned(lr, |qi, ri| {
                        if qi >= seq.len() || ri >= rs.len() {
                            return;
                        }
                        let Some(ctx) = tri_ctx(rs, ri) else { return };
                        let cycle = if rev { len - 1 - qi } else { qi } * RECAL_CYCLES / len;
                        if let Some(nq) = self.empirical(qual[qi], ctx, rev, cycle.min(RECAL_CYCLES - 1)) {
                            updates.push((qi, nq));
                        }
                    });
                }
                let (_, qual) = lr.seq_qual_mut();
                for (qi, nq) in updates {
                    qual[qi] = nq;
                }
            }
        });
    }
}
