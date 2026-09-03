//! Pileup engine: walks sorted reads per sample position by position. Each
//! read keeps a CIGAR cursor that only moves forward, so a pileup pass costs
//! O(read length) per read instead of a CIGAR walk from the start at every
//! position. Per-site buffers are reused across positions.

use anyhow::Result;
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use smallvec::SmallVec;
use std::collections::VecDeque;

use crate::bam::pos_filter::InterestingMap;

/// CIGAR ops, inline for typical short reads (one heap block for long ones).
pub type CigarOps = SmallVec<[(Kind, u32); 8]>;

/// Per-sample observations at one site. Counts are split by strand so
/// ADF/ADR and the strand-bias test can be derived.
#[derive(Default, Clone)]
pub struct SampleSiteCounts {
    pub depth: u32,
    pub base_counts: [u32; 4],
    pub base_quals: [u32; 4],
    /// Forward-strand part of `base_counts`.
    pub fwd_counts: [u32; 4],
    pub mq_sum: u64,
    /// Reads contributing a base here that carry a soft clip (bcftools SCR).
    pub n_softclip: u32,
    /// Per-read `(ACGT index | strand << 4, MAPQ-capped quality)` for the
    /// genotype-likelihood model. `base_counts`/`base_quals` stay raw (for
    /// AD/QS); only this carries the `min(BQ, MAPQ)` cap.
    pub obs: Vec<(u8, u8)>,
    /// Insertions after this position: (sequence, count, forward count).
    pub ins_alleles: Vec<(String, u32, u32)>,
    /// Deletions after this position: (length, count, forward count).
    pub del_alleles: Vec<(u32, u32, u32)>,
}

impl SampleSiteCounts {
    #[inline]
    pub fn rev_count(&self, i: usize) -> u32 {
        self.base_counts[i].saturating_sub(self.fwd_counts[i])
    }

    /// Reset for the next site, keeping the allocations.
    fn clear(&mut self) {
        self.depth = 0;
        self.base_counts = [0; 4];
        self.base_quals = [0; 4];
        self.fwd_counts = [0; 4];
        self.mq_sum = 0;
        self.n_softclip = 0;
        self.obs.clear();
        self.ins_alleles.clear();
        self.del_alleles.clear();
    }
}

#[derive(Default, Clone)]
pub struct PileupSite {
    pub ref_id: usize,
    pub pos: u32,
    pub per_sample: Vec<SampleSiteCounts>,
}

impl PileupSite {
    pub fn total_depth(&self) -> u32 { self.per_sample.iter().map(|s| s.depth).sum() }

    pub fn aggregated(&self) -> SampleSiteCounts {
        let mut a = SampleSiteCounts::default();
        for s in &self.per_sample {
            a.depth += s.depth;
            for i in 0..4 {
                a.base_counts[i] += s.base_counts[i];
                a.base_quals[i] += s.base_quals[i];
                a.fwd_counts[i] += s.fwd_counts[i];
            }
            a.mq_sum += s.mq_sum;
            a.n_softclip += s.n_softclip;
            for (seq, c, f) in &s.ins_alleles {
                if let Some(e) = a.ins_alleles.iter_mut().find(|(s, _, _)| s == seq) { e.1 += c; e.2 += f; } else { a.ins_alleles.push((seq.clone(), *c, *f)); }
            }
            for (len, c, f) in &s.del_alleles {
                if let Some(e) = a.del_alleles.iter_mut().find(|(l, _, _)| l == len) { e.1 += c; e.2 += f; } else { a.del_alleles.push((*len, *c, *f)); }
            }
        }
        a
    }

    /// Prepare for a new site, reusing the per-sample buffers.
    fn reset(&mut self, ref_id: usize, pos: u32, n_samples: usize) {
        self.ref_id = ref_id;
        self.pos = pos;
        self.per_sample.resize_with(n_samples, SampleSiteCounts::default);
        for s in &mut self.per_sample {
            s.clear();
        }
    }
}

#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b { b'A' | b'a' => Some(0), b'C' | b'c' => Some(1), b'G' | b'g' => Some(2), b'T' | b't' => Some(3), _ => None }
}

#[inline]
fn is_aligned(k: Kind) -> bool {
    matches!(k, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch)
}

/// One aligned read. Bases and qualities share one allocation; the CIGAR
/// cursor (`cur_*`) only moves forward, so a pileup walk visits every op once.
#[derive(Clone)]
pub struct LiveRead {
    /// Query bases followed by their qualities.
    data: Vec<u8>,
    len: u32,
    pub cigar_pairs: CigarOps,
    pub ref_start: u32,
    pub ref_end_cached: u32,
    pub ref_id: usize,
    pub mapq: u8,
    pub sample_idx: usize,
    pub flags: u16,
    pub has_softclip: bool,
    /// Hash of the read name, for pairing overlapping mates.
    pub qname_hash: u64,
    /// 0-based start of a mapped mate on the same contig, else `u32::MAX`.
    pub mate_pos: u32,
    cur_op: u32,
    cur_r: u32,
    cur_q: u32,
}

impl LiveRead {
    #[allow(clippy::too_many_arguments)]
    pub fn new(seq: &[u8], qual: &[u8], cigar: CigarOps, ref_start: u32, ref_id: usize, mapq: u8, sample_idx: usize, flags: u16) -> Self {
        let mut data = Vec::with_capacity(seq.len() * 2);
        data.extend_from_slice(seq);
        push_quals(&mut data, qual.iter().copied(), seq.len());
        Self::from_parts(data, seq.len(), cigar, ref_start, ref_id, mapq, sample_idx, flags, 0, u32::MAX)
    }

    /// Attach mate information (see [`LiveRead::qname_hash`], [`LiveRead::mate_pos`]).
    pub fn with_mate(mut self, qname_hash: u64, mate_pos: u32) -> Self {
        self.qname_hash = qname_hash;
        self.mate_pos = mate_pos;
        self
    }

    #[allow(clippy::too_many_arguments)]
    fn from_parts(data: Vec<u8>, len: usize, cigar: CigarOps, ref_start: u32, ref_id: usize, mapq: u8, sample_idx: usize, flags: u16, qname_hash: u64, mate_pos: u32) -> Self {
        let has_softclip = cigar.iter().any(|(k, _)| *k == Kind::SoftClip);
        let ref_end_cached = compute_ref_end(ref_start, &cigar);
        Self {
            data,
            len: len as u32,
            cigar_pairs: cigar,
            ref_start,
            ref_end_cached,
            ref_id,
            mapq,
            sample_idx,
            flags,
            has_softclip,
            qname_hash,
            mate_pos,
            cur_op: 0,
            cur_r: ref_start,
            cur_q: 0,
        }
    }

    /// `(reference position, query index)` of every aligned (M/=/X) base.
    pub fn aligned_positions(&self) -> Vec<(u32, u32)> {
        let mut out = Vec::with_capacity(self.len as usize);
        let (mut r, mut q) = (self.ref_start, 0u32);
        for &(kind, len) in &self.cigar_pairs {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for i in 0..len {
                        out.push((r + i, q + i));
                    }
                    r += len;
                    q += len;
                }
                Kind::Insertion | Kind::SoftClip => q += len,
                Kind::Deletion | Kind::Skip => r += len,
                _ => {}
            }
        }
        out
    }

    #[inline]
    pub fn seq(&self) -> &[u8] { &self.data[..self.len as usize] }

    #[inline]
    pub fn qual(&self) -> &[u8] { &self.data[self.len as usize..] }

    /// Bases and mutable qualities (BAQ rewrites qualities in place).
    #[inline]
    pub fn seq_qual_mut(&mut self) -> (&[u8], &mut [u8]) {
        let (s, q) = self.data.split_at_mut(self.len as usize);
        (s, q)
    }

    pub fn skip_by_flags(&self) -> bool {
        let f = self.flags;
        let unmapped = f & 0x4 != 0;
        let secondary = f & 0x100 != 0;
        let supplementary = f & 0x800 != 0;
        let duplicate = f & 0x400 != 0;
        let qcfail = f & 0x200 != 0;
        unmapped || secondary || supplementary || duplicate || qcfail
    }

    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.flags & 0x10 != 0
    }

    #[inline]
    pub fn ref_end(&self) -> u32 {
        self.ref_end_cached
    }

    /// Move the cursor to the op covering `ref_pos`; positions only increase.
    #[inline]
    fn seek(&mut self, ref_pos: u32) {
        while let Some(&(kind, len)) = self.cigar_pairs.get(self.cur_op as usize) {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if ref_pos < self.cur_r + len { return; }
                    self.cur_r += len;
                    self.cur_q += len;
                }
                Kind::Deletion | Kind::Skip => {
                    if ref_pos < self.cur_r + len { return; }
                    self.cur_r += len;
                }
                Kind::Insertion | Kind::SoftClip => self.cur_q += len,
                Kind::HardClip | Kind::Pad => {}
            }
            self.cur_op += 1;
        }
    }

    /// Base and quality at `ref_pos` after [`seek`]; `None` inside a
    /// deletion/skip or past the read.
    #[inline]
    fn base_at_cursor(&self, ref_pos: u32) -> Option<(u8, u8)> {
        let &(kind, len) = self.cigar_pairs.get(self.cur_op as usize)?;
        if !is_aligned(kind) || ref_pos < self.cur_r || ref_pos >= self.cur_r + len {
            return None;
        }
        let qi = (self.cur_q + (ref_pos - self.cur_r)) as usize;
        if qi >= self.len as usize {
            return None;
        }
        Some((self.data[qi], self.data[self.len as usize + qi]))
    }

    /// The indel that follows `ref_pos` when it is the last aligned base
    /// before an insertion/deletion op (after [`seek`]).
    #[inline]
    fn indel_after_cursor(&self, ref_pos: u32) -> Option<IndelEvent<'_>> {
        let &(kind, len) = self.cigar_pairs.get(self.cur_op as usize)?;
        if !is_aligned(kind) || ref_pos + 1 != self.cur_r + len {
            return None;
        }
        let &(next, nlen) = self.cigar_pairs.get(self.cur_op as usize + 1)?;
        match next {
            Kind::Insertion => {
                let start = (self.cur_q + len) as usize;
                let end = (start + nlen as usize).min(self.len as usize);
                (start < end).then(|| IndelEvent::Ins(&self.data[start..end]))
            }
            Kind::Deletion => Some(IndelEvent::Del(nlen)),
            _ => None,
        }
    }

    /// Query index range `[s, e)` of the bases spanning reference window
    /// `[w_lo, w_hi)`, including bases the read inserts inside it. `None`
    /// unless the read fully covers the window.
    fn query_span(&self, w_lo: u32, w_hi: u32) -> Option<(usize, usize)> {
        if self.ref_start > w_lo || self.ref_end() < w_hi {
            return None;
        }
        let mut r = self.ref_start;
        let mut q: u32 = 0;
        // Query offset after the last reference-consuming op, so a window
        // ending exactly at the alignment end resolves.
        let mut q_ref_end: u32 = 0;
        let mut qs: Option<u32> = None;
        let mut qe: Option<u32> = None;
        for &(kind, len) in &self.cigar_pairs {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let r2 = r + len;
                    if qs.is_none() && w_lo >= r && w_lo < r2 {
                        qs = Some(q + (w_lo - r));
                    }
                    if qe.is_none() && w_hi >= r && w_hi < r2 {
                        qe = Some(q + (w_hi - r));
                    }
                    r = r2;
                    q += len;
                    q_ref_end = q;
                }
                Kind::Insertion | Kind::SoftClip => {
                    q += len;
                }
                Kind::Deletion | Kind::Skip => {
                    let r2 = r + len;
                    if qs.is_none() && w_lo >= r && w_lo < r2 {
                        qs = Some(q);
                    }
                    if qe.is_none() && w_hi >= r && w_hi < r2 {
                        qe = Some(q);
                    }
                    r = r2;
                    q_ref_end = q;
                }
                _ => {}
            }
            if qs.is_some() && qe.is_some() {
                break;
            }
        }
        if qe.is_none() && w_hi == self.ref_end() {
            qe = Some(q_ref_end);
        }
        let (s, e) = (qs? as usize, qe? as usize);
        if e <= s || e > self.len as usize {
            return None;
        }
        Some((s, e))
    }

    /// The read's query bases spanning reference window [w_lo, w_hi) (0-based ref coords),
    /// INCLUDING any bases the read inserts inside the window (so a read carrying an
    /// insertion yields a longer substring than a ref-matching read). Returns None unless
    /// the read fully spans the window. Used for local indel realignment.
    pub fn query_window(&self, w_lo: u32, w_hi: u32) -> Option<Vec<u8>> {
        let (s, e) = self.query_span(w_lo, w_hi)?;
        Some(self.seq()[s..e].to_vec())
    }

    /// [`query_window`] with the matching base qualities.
    pub fn query_window_qual(&self, w_lo: u32, w_hi: u32) -> Option<(Vec<u8>, Vec<u8>)> {
        let (s, e) = self.query_span(w_lo, w_hi)?;
        Some((self.seq()[s..e].to_vec(), self.qual()[s..e].to_vec()))
    }
}

/// htslib `tweak_overlap_quality`: where two mates of one fragment cover the
/// same reference base, agreeing bases pool their quality (capped at 200)
/// on one mate and zero the other; disagreeing bases keep 80% of the better
/// one and zero the other. Which mate keeps the pooled quality is decided
/// by a bit of the read-name hash, as htslib does. The fragment then counts
/// once in depth and likelihoods.
pub fn tweak_overlap_quality(a: &mut LiveRead, b: &mut LiveRead) {
    let pa = a.aligned_positions();
    let pb = b.aligned_positions();
    let (la, lb) = (a.len as usize, b.len as usize);
    let a_keeps = a.qname_hash & 1 == 1;
    let (mut i, mut j) = (0usize, 0usize);
    while i < pa.len() && j < pb.len() {
        let (ra, qa) = pa[i];
        let (rb, qb) = pb[j];
        if ra < rb {
            i += 1;
            continue;
        }
        if rb < ra {
            j += 1;
            continue;
        }
        let (qa, qb) = (qa as usize, qb as usize);
        if qa >= la || qb >= lb {
            return;
        }
        let (ba, bb) = (a.data[qa], b.data[qb]);
        let (qual_a, qual_b) = (a.data[la + qa], b.data[lb + qb]);
        if ba.eq_ignore_ascii_case(&bb) {
            let pooled = (qual_a as u32 + qual_b as u32).min(200) as u8;
            a.data[la + qa] = if a_keeps { pooled } else { 0 };
            b.data[lb + qb] = if a_keeps { 0 } else { pooled };
        } else if qual_a > qual_b {
            a.data[la + qa] = (0.8 * qual_a as f64) as u8;
            b.data[lb + qb] = 0;
        } else if qual_a < qual_b {
            b.data[lb + qb] = (0.8 * qual_b as f64) as u8;
            a.data[la + qa] = 0;
        } else {
            a.data[la + qa] = if a_keeps { (0.8 * qual_a as f64) as u8 } else { 0 };
            b.data[lb + qb] = if a_keeps { 0 } else { (0.8 * qual_b as f64) as u8 };
        }
        i += 1;
        j += 1;
    }
}

/// Add a read to the live set, first folding it into an overlapping mate
/// that is already there (`bcftools mpileup` default; `-x` disables).
/// Like htslib's `overlap_push`, only proper pairs with a mapped mate on
/// the same contig take part.
fn push_live(live: &mut VecDeque<LiveRead>, mut lr: LiveRead, overlaps: bool) {
    if overlaps && lr.flags & 0x2 != 0 && lr.mate_pos != u32::MAX && lr.mate_pos <= lr.ref_start {
        let mate = live
            .iter_mut()
            .rev()
            .find(|m| m.sample_idx == lr.sample_idx && m.qname_hash == lr.qname_hash && m.ref_start == lr.mate_pos && m.mate_pos == lr.ref_start);
        if let Some(mate) = mate {
            tweak_overlap_quality(mate, &mut lr);
        }
    }
    live.push_back(lr);
}

/// Append `n` qualities, padding with 30 when the record carries fewer.
fn push_quals(data: &mut Vec<u8>, quals: impl Iterator<Item = u8>, n: usize) {
    let before = data.len();
    data.extend(quals.take(n));
    if data.len() - before < n {
        data.resize(before + n, 30);
    }
}

enum IndelEvent<'a> {
    Ins(&'a [u8]),
    Del(u32),
}

/// Record one read's observation at `cur_pos` into `site`.
#[inline]
fn observe(lr: &mut LiveRead, cur_pos: u32, min_bq: u8, skip_indels: bool, site: &mut PileupSite) {
    lr.seek(cur_pos);
    if let Some((b, q)) = lr.base_at_cursor(cur_pos) {
        if q >= min_bq {
            if let Some(i) = base_to_idx(b) {
                let s = &mut site.per_sample[lr.sample_idx];
                s.base_counts[i] += 1;
                s.base_quals[i] += q as u32;
                if !lr.is_reverse() { s.fwd_counts[i] += 1; }
                if lr.has_softclip { s.n_softclip += 1; }
                s.mq_sum += lr.mapq as u64;
                // Base index in the low nibble, strand in bit 4 (see
                // `errmod::pack_base`); quality capped by MAPQ and `capQ` (60).
                s.obs.push(((i as u8) | ((lr.is_reverse() as u8) << 4), q.min(lr.mapq.min(60))));
                s.depth += 1;
            }
        }
    }
    if skip_indels {
        return;
    }
    let f = u32::from(!lr.is_reverse());
    let s = &mut site.per_sample[lr.sample_idx];
    match lr.indel_after_cursor(cur_pos) {
        Some(IndelEvent::Ins(seq)) => {
            if let Some(e) = s.ins_alleles.iter_mut().find(|(k, _, _)| k.as_bytes() == seq) {
                e.1 += 1;
                e.2 += f;
            } else {
                s.ins_alleles.push((String::from_utf8_lossy(seq).into_owned(), 1, f));
            }
        }
        Some(IndelEvent::Del(l)) => {
            if let Some(e) = s.del_alleles.iter_mut().find(|(k, _, _)| *k == l) {
                e.1 += 1;
                e.2 += f;
            } else {
                s.del_alleles.push((l, 1, f));
            }
        }
        None => {}
    }
}

/// A sorted per-sample stream of reads.
pub trait ReadSource {
    /// Next read passing the flag and MAPQ filters, tagged with `sample_idx`.
    fn next_read(&mut self, min_mq: u8, sample_idx: usize) -> Option<LiveRead>;
}

impl ReadSource for std::vec::IntoIter<LiveRead> {
    fn next_read(&mut self, min_mq: u8, sample_idx: usize) -> Option<LiveRead> {
        for mut lr in self.by_ref() {
            if lr.skip_by_flags() || lr.mapq < min_mq { continue; }
            lr.sample_idx = sample_idx;
            return Some(lr);
        }
        None
    }
}

impl ReadSource for crossbeam_channel::IntoIter<LiveRead> {
    fn next_read(&mut self, min_mq: u8, sample_idx: usize) -> Option<LiveRead> {
        for mut lr in self.by_ref() {
            if lr.skip_by_flags() || lr.mapq < min_mq { continue; }
            lr.sample_idx = sample_idx;
            return Some(lr);
        }
        None
    }
}

/// Position-by-position walk over `sources` (one per sample). `emit` sees
/// every site with depth and the reads overlapping it. `pos_filter` skips
/// to the next listed position; `overlaps` enables mate-overlap merging.
#[allow(clippy::too_many_arguments)]
pub fn run_engine<S: ReadSource, F>(
    mut sources: Vec<S>,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    overlaps: bool,
    pos_filter: Option<&InterestingMap>,
    emit: &mut F,
) -> Result<()>
where
    F: FnMut(&PileupSite, &[LiveRead]),
{
    let n_samples = sources.len();
    let mut next: Vec<Option<LiveRead>> = sources.iter_mut().enumerate().map(|(i, s)| s.next_read(min_mq, i)).collect();
    let mut live: VecDeque<LiveRead> = VecDeque::new();
    let mut site = PileupSite::default();
    let mut cur_ref: i32 = -1;
    let mut cur_pos: u32 = 0;

    loop {
        if live.is_empty() {
            // Jump to the earliest pending read.
            let mut pick: Option<(usize, i32, u32)> = None;
            for (i, r) in next.iter().enumerate() {
                if let Some(lr) = r {
                    let cand = (i, lr.ref_id as i32, lr.ref_start);
                    pick = Some(match pick {
                        None => cand,
                        Some(p) => if (cand.1, cand.2) < (p.1, p.2) { cand } else { p },
                    });
                }
            }
            let Some((idx, rid, p)) = pick else { break };
            cur_ref = rid;
            cur_pos = p;
            push_live(&mut live, next[idx].take().unwrap(), overlaps);
            next[idx] = sources[idx].next_read(min_mq, idx);
            continue;
        }

        if let Some(filter) = pos_filter {
            match filter.next_at_or_after(cur_ref as usize, cur_pos) {
                Some(next_pos) => {
                    if next_pos > cur_pos {
                        cur_pos = next_pos;
                        live.retain(|lr| lr.ref_end() > cur_pos);
                    }
                }
                None => {
                    live.clear();
                    continue;
                }
            }
        }

        for i in 0..n_samples {
            while let Some(lr) = &next[i] {
                if lr.ref_id as i32 == cur_ref && lr.ref_start <= cur_pos {
                    push_live(&mut live, next[i].take().unwrap(), overlaps);
                    next[i] = sources[i].next_read(min_mq, i);
                } else {
                    break;
                }
            }
        }

        live.retain(|lr| lr.ref_end() > cur_pos);
        site.reset(cur_ref as usize, cur_pos, n_samples);
        for lr in live.iter_mut() {
            observe(lr, cur_pos, min_bq, skip_indels, &mut site);
        }
        if site.total_depth() > 0 {
            // Every live read starts at or before `cur_pos` and ends after it.
            emit(&site, live.make_contiguous());
        }
        cur_pos += 1;
    }
    Ok(())
}

pub fn mpileup_engine<F>(bam: &mut crate::bam::BamReader, min_mq: u8, min_bq: u8, skip_indels: bool, overlaps: bool, mut emit: F) -> Result<()>
where
    F: FnMut(&PileupSite, &[LiveRead]),
{
    mpileup_engine_multi(std::slice::from_mut(bam), min_mq, min_bq, skip_indels, overlaps, &mut emit)
}

/// Engine over decoder threads: nothing but the channel buffers is in RAM.
/// Decoder errors surface after the walk.
pub fn mpileup_engine_streaming<F>(streams: Vec<crate::bam::StreamingBam>, min_mq: u8, min_bq: u8, skip_indels: bool, overlaps: bool, emit: &mut F) -> Result<()>
where
    F: FnMut(&PileupSite, &[LiveRead]),
{
    let mut handles = Vec::with_capacity(streams.len());
    let sources: Vec<crossbeam_channel::IntoIter<LiveRead>> = streams
        .into_iter()
        .map(|s| {
            handles.push(s.handle);
            s.rx.into_iter()
        })
        .collect();
    run_engine(sources, min_mq, min_bq, skip_indels, overlaps, None, emit)?;
    for h in handles.into_iter().flatten() {
        h.join().map_err(|_| anyhow::anyhow!("BAM decoder thread panicked"))??;
    }
    Ok(())
}

pub fn mpileup_engine_multi<F>(bams: &mut [crate::bam::BamReader], min_mq: u8, min_bq: u8, skip_indels: bool, overlaps: bool, emit: &mut F) -> Result<()>
where
    F: FnMut(&PileupSite, &[LiveRead]),
{
    let records_per_sample: Vec<Vec<LiveRead>> = bams.iter_mut().map(|b| std::mem::take(&mut b.records_buf)).collect();
    mpileup_engine_from_records(records_per_sample, min_mq, min_bq, skip_indels, overlaps, None, emit)
}

/// Engine entry point that takes owned records directly. Pass `pos_filter`
/// to visit only listed positions (the engine jumps to the next one).
#[allow(clippy::too_many_arguments)]
pub fn mpileup_engine_from_records<F>(
    records_per_sample: Vec<Vec<LiveRead>>,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    overlaps: bool,
    pos_filter: Option<&InterestingMap>,
    emit: &mut F,
) -> Result<()>
where
    F: FnMut(&PileupSite, &[LiveRead]),
{
    let sources: Vec<std::vec::IntoIter<LiveRead>> = records_per_sample.into_iter().map(|v| v.into_iter()).collect();
    run_engine(sources, min_mq, min_bq, skip_indels, overlaps, pos_filter, emit)
}

/// Mate start for overlap detection: paired, mate mapped on the same contig.
fn mate_pos_of(flags: u16, ref_id: usize, mate_ref: Option<usize>, mate_start: Option<usize>) -> u32 {
    let paired = flags & 0x1 != 0;
    let mate_unmapped = flags & 0x8 != 0;
    match (paired && !mate_unmapped, mate_ref, mate_start) {
        (true, Some(mr), Some(ms)) if mr == ref_id && ms > 0 => (ms - 1) as u32,
        _ => u32::MAX,
    }
}

pub fn build_live_from_bam(rec: &bam::Record, sample_idx: usize) -> Option<LiveRead> {
    let start = usize::from(rec.alignment_start().transpose().ok().flatten()?) as u32 - 1;
    let ref_id = rec.reference_sequence_id().transpose().ok().flatten()?;
    // MAPQ 255 (unavailable) counts as 20, like `bcf_call_glfgen`.
    let mapq: u8 = rec.mapping_quality().map(u8::from).unwrap_or(20);
    let flags = u16::from(rec.flags());
    let seq = rec.sequence();
    let n = seq.len();
    let mut data = Vec::with_capacity(n * 2);
    data.extend(seq.iter());
    push_quals(&mut data, rec.quality_scores().iter(), n);
    let mut cigar = CigarOps::new();
    for op in rec.cigar().iter() {
        let op = op.ok()?;
        cigar.push((op.kind(), op.len() as u32));
    }
    let qname_hash = rec.name().map(|n| fxhash::hash64(n.as_ref() as &[u8])).unwrap_or(0);
    let mate_ref = rec.mate_reference_sequence_id().transpose().ok().flatten();
    let mate_start = rec.mate_alignment_start().transpose().ok().flatten().map(usize::from);
    let mate_pos = mate_pos_of(flags, ref_id, mate_ref, mate_start);
    Some(LiveRead::from_parts(data, n, cigar, start, ref_id, mapq, sample_idx, flags, qname_hash, mate_pos))
}

pub fn build_live_from_cram(rb: &noodles_sam::alignment::RecordBuf, sample_idx: usize) -> Option<LiveRead> {
    use noodles_sam::alignment::record::cigar::Op;
    let start = usize::from(rb.alignment_start()?) as u32 - 1;
    let ref_id = rb.reference_sequence_id()?;
    let mapq: u8 = rb.mapping_quality().map(u8::from).unwrap_or(20);
    let flags = u16::from(rb.flags());
    let seq: &[u8] = rb.sequence().as_ref();
    let qual: &[u8] = rb.quality_scores().as_ref();
    let cigar_ops: &[Op] = rb.cigar().as_ref();
    if cigar_ops.is_empty() { return None; }
    let cigar: CigarOps = cigar_ops.iter().map(|op| (op.kind(), op.len() as u32)).collect();
    let qname_hash = rb.name().map(|n| fxhash::hash64(n.as_ref() as &[u8])).unwrap_or(0);
    let mate_pos = mate_pos_of(flags, ref_id, rb.mate_reference_sequence_id(), rb.mate_alignment_start().map(usize::from));
    Some(LiveRead::new(seq, qual, cigar, start, ref_id, mapq, sample_idx, flags).with_mate(qname_hash, mate_pos))
}

#[inline]
fn compute_ref_end(start: u32, cigar_pairs: &[(Kind, u32)]) -> u32 {
    let mut r = start;
    for &(k, l) in cigar_pairs {
        if matches!(k, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip) {
            r += l;
        }
    }
    r
}

/// Two-tailed Fisher exact test p-value for a 2x2 table (htslib `kt_fisher_exact`).
pub fn fisher_exact_two_tailed(n11: u32, n12: u32, n21: u32, n22: u32) -> f64 {
    let n = n11 + n12 + n21 + n22;
    if n == 0 { return 1.0; }
    let r1 = n11 + n12;
    let c1 = n11 + n21;
    let lnf = |k: u32| -> f64 {
        let mut s = 0.0;
        for i in 2..=k { s += (i as f64).ln(); }
        s
    };
    let ln_const = lnf(r1) + lnf(n - r1) + lnf(c1) + lnf(n - c1) - lnf(n);
    let prob = |a: u32| -> f64 {
        let b = r1 - a;
        let c = c1 - a;
        let d = n - r1 - c;
        (ln_const - lnf(a) - lnf(b) - lnf(c) - lnf(d)).exp()
    };
    let lo = r1.saturating_sub(n - c1);
    let hi = r1.min(c1);
    let p_obs = prob(n11);
    let mut p = 0.0;
    for a in lo..=hi {
        let pa = prob(a);
        if pa <= p_obs * (1.0 + 1e-7) { p += pa; }
    }
    p.min(1.0)
}

/// Phred-scaled strand bias (bcftools SP): Fisher test on ref/alt by strand.
pub fn strand_bias_phred(ref_fwd: u32, ref_rev: u32, alt_fwd: u32, alt_rev: u32) -> u32 {
    let p = fisher_exact_two_tailed(ref_fwd, ref_rev, alt_fwd, alt_rev);
    if p <= 0.0 { return 255; }
    ((-10.0 * p.log10()).round() as i64).clamp(0, 255) as u32
}

#[cfg(test)]
#[path = "../../tests/unit/bam_pileup.rs"]
mod tests;
