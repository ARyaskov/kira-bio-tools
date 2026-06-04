use anyhow::Result;
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use std::collections::VecDeque;

#[derive(Default, Clone)]
pub struct SampleSiteCounts {
    pub depth: u32,
    pub base_counts: [u32; 4],
    pub base_quals: [u32; 4],
    pub mq_sum: u64,
    pub ins_alleles: Vec<(String, u32)>,
    pub del_alleles: Vec<(u32, u32)>,
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
            for i in 0..4 { a.base_counts[i] += s.base_counts[i]; a.base_quals[i] += s.base_quals[i]; }
            a.mq_sum += s.mq_sum;
            for (seq, c) in &s.ins_alleles {
                if let Some(e) = a.ins_alleles.iter_mut().find(|(s, _)| s == seq) { e.1 += c; } else { a.ins_alleles.push((seq.clone(), *c)); }
            }
            for (len, c) in &s.del_alleles {
                if let Some(e) = a.del_alleles.iter_mut().find(|(l, _)| l == len) { e.1 += c; } else { a.del_alleles.push((*len, *c)); }
            }
        }
        a
    }
}

#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b { b'A' | b'a' => Some(0), b'C' | b'c' => Some(1), b'G' | b'g' => Some(2), b'T' | b't' => Some(3), _ => None }
}

#[derive(Clone)]
pub struct LiveRead {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub cigar_pairs: Vec<(Kind, u32)>,
    pub ref_start: u32,
    pub ref_end_cached: u32,
    pub ref_id: usize,
    pub mapq: u8,
    pub sample_idx: usize,
    pub flags: u16,
}

impl LiveRead {
    pub fn skip_by_flags(&self) -> bool {
        let f = self.flags;
        let unmapped = f & 0x4 != 0;
        let secondary = f & 0x100 != 0;
        let supplementary = f & 0x800 != 0;
        let duplicate = f & 0x400 != 0;
        let qcfail = f & 0x200 != 0;
        unmapped || secondary || supplementary || duplicate || qcfail
    }
}

impl LiveRead {
    fn base_at(&self, ref_pos: u32) -> Option<(u8, u8)> {
        let mut r_pos = self.ref_start;
        let mut q_pos: u32 = 0;
        for &(kind, len) in &self.cigar_pairs {
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if ref_pos >= r_pos && ref_pos < r_pos + len {
                        let off = (ref_pos - r_pos) as usize;
                        let qi = (q_pos as usize) + off;
                        if qi < self.seq.len() {
                            let b = self.seq[qi];
                            let q = if qi < self.qual.len() { self.qual[qi] } else { 30 };
                            return Some((b, q));
                        }
                        return None;
                    }
                    r_pos += len;
                    q_pos += len;
                }
                Kind::Insertion | Kind::SoftClip => q_pos += len,
                Kind::Deletion | Kind::Skip => {
                    if ref_pos >= r_pos && ref_pos < r_pos + len { return None; }
                    r_pos += len;
                }
                Kind::HardClip | Kind::Pad => {}
            }
            if r_pos > ref_pos { return None; }
        }
        None
    }

    #[inline]
    pub fn ref_end(&self) -> u32 {
        self.ref_end_cached
    }

    /// The read's query bases spanning reference window [w_lo, w_hi) (0-based ref coords),
    /// INCLUDING any bases the read inserts inside the window (so a read carrying an
    /// insertion yields a longer substring than a ref-matching read). Returns None unless
    /// the read fully spans the window. Used for local indel realignment.
    pub fn query_window(&self, w_lo: u32, w_hi: u32) -> Option<Vec<u8>> {
        if self.ref_start > w_lo || self.ref_end() < w_hi {
            return None;
        }
        let mut r = self.ref_start;
        let mut q: u32 = 0;
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
                }
                _ => {}
            }
            if qs.is_some() && qe.is_some() {
                break;
            }
        }
        let (s, e) = (qs?, qe?);
        if e <= s || e as usize > self.seq.len() {
            return None;
        }
        Some(self.seq[s as usize..e as usize].to_vec())
    }

    fn indel_after(&self, ref_pos: u32) -> Option<IndelEvent> {
        let mut r_pos = self.ref_start;
        let mut q_pos: u32 = 0;
        let mut found_anchor = false;
        for &(kind, len) in &self.cigar_pairs {
            if found_anchor {
                return match kind {
                    Kind::Insertion => {
                        let start = q_pos as usize;
                        let end = (start + len as usize).min(self.seq.len());
                        if start >= self.seq.len() { return None; }
                        Some(IndelEvent::Ins(self.seq[start..end].to_vec()))
                    }
                    Kind::Deletion => Some(IndelEvent::Del(len)),
                    _ => None,
                };
            }
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let end = r_pos + len;
                    if ref_pos >= r_pos && ref_pos < end && ref_pos + 1 == end {
                        found_anchor = true;
                    }
                    r_pos = end;
                    q_pos += len;
                }
                Kind::Insertion | Kind::SoftClip => { q_pos += len; }
                Kind::Deletion | Kind::Skip => { r_pos += len; }
                _ => {}
            }
            if r_pos > ref_pos + 1 && !found_anchor { break; }
        }
        None
    }
}

enum IndelEvent { Ins(Vec<u8>), Del(u32) }

pub fn mpileup_engine<F>(
    bam: &mut crate::bam::BamReader,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    mut emit: F,
) -> Result<()>
where F: FnMut(&PileupSite, &[&LiveRead]),
{
    mpileup_engine_multi(std::slice::from_mut(bam), min_mq, min_bq, skip_indels, &mut emit)
}

pub fn mpileup_engine_streaming<F>(
    streams: Vec<crate::bam::StreamingBam>,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    emit: &mut F,
) -> Result<()>
where F: FnMut(&PileupSite, &[&LiveRead]),
{
    let n_samples = streams.len();
    let mut iters: Vec<crossbeam_channel::IntoIter<LiveRead>> = streams.into_iter()
        .map(|s| s.rx.into_iter()).collect();
    let mut next: Vec<Option<LiveRead>> = iters.iter_mut()
        .enumerate().map(|(i, it)| next_aligned_chan(it, min_mq, i)).collect();

    let mut live: std::collections::VecDeque<LiveRead> = std::collections::VecDeque::new();
    let mut cur_ref: i32 = -1;
    let mut cur_pos: u32 = 0;

    loop {
        let mut any = !live.is_empty();
        for r in &next { if r.is_some() { any = true; break; } }
        if !any { break; }

        if live.is_empty() {
            let mut pick: Option<(usize, i32, u32)> = None;
            for (i, r) in next.iter().enumerate() {
                if let Some(lr) = r {
                    let cand = (i, lr.ref_id as i32, lr.ref_start);
                    pick = Some(match pick { None => cand, Some(p) => if (cand.1, cand.2) < (p.1, p.2) { cand } else { p } });
                }
            }
            let Some((idx, rid, p)) = pick else { break; };
            cur_ref = rid; cur_pos = p;
            let lr = next[idx].take().unwrap();
            live.push_back(lr);
            next[idx] = next_aligned_chan(&mut iters[idx], min_mq, idx);
            continue;
        }

        for i in 0..n_samples {
            while let Some(lr) = &next[i] {
                if lr.ref_id as i32 == cur_ref && lr.ref_start <= cur_pos {
                    let taken = next[i].take().unwrap();
                    live.push_back(taken);
                    next[i] = next_aligned_chan(&mut iters[i], min_mq, i);
                } else { break; }
            }
        }

        let mut site = PileupSite { ref_id: cur_ref as usize, pos: cur_pos, per_sample: vec![SampleSiteCounts::default(); n_samples] };
        let mut keep: std::collections::VecDeque<LiveRead> = std::collections::VecDeque::with_capacity(live.len());
        while let Some(lr) = live.pop_front() {
            if lr.ref_end() <= cur_pos { continue; }
            if let Some((b, q)) = lr.base_at(cur_pos) {
                if q >= min_bq {
                    if let Some(i) = base_to_idx(b) {
                        let s = &mut site.per_sample[lr.sample_idx];
                        s.base_counts[i] += 1;
                        s.base_quals[i] += q as u32;
                        s.mq_sum += lr.mapq as u64;
                        s.depth += 1;
                    }
                }
            }
            if !skip_indels {
                if let Some(ev) = lr.indel_after(cur_pos) {
                    let s = &mut site.per_sample[lr.sample_idx];
                    match ev {
                        IndelEvent::Ins(seq) => {
                            let key = String::from_utf8_lossy(&seq).into_owned();
                            if let Some(e) = s.ins_alleles.iter_mut().find(|(k, _)| *k == key) { e.1 += 1; } else { s.ins_alleles.push((key, 1)); }
                        }
                        IndelEvent::Del(l) => {
                            if let Some(e) = s.del_alleles.iter_mut().find(|(k, _)| *k == l) { e.1 += 1; } else { s.del_alleles.push((l, 1)); }
                        }
                    }
                }
            }
            keep.push_back(lr);
        }
        live = keep;
        if site.total_depth() > 0 {
            let overlapping: Vec<&LiveRead> = live.iter().filter(|lr| lr.ref_start <= cur_pos && lr.ref_end() > cur_pos).collect();
            emit(&site, &overlapping);
        }
        cur_pos += 1;
    }
    Ok(())
}

fn next_aligned_chan(iter: &mut crossbeam_channel::IntoIter<LiveRead>, min_mq: u8, sample_idx: usize) -> Option<LiveRead> {
    for mut lr in iter.by_ref() {
        if lr.skip_by_flags() || lr.mapq < min_mq { continue; }
        lr.sample_idx = sample_idx;
        return Some(lr);
    }
    None
}

pub fn mpileup_engine_multi<F>(
    bams: &mut [crate::bam::BamReader],
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    emit: &mut F,
) -> Result<()>
where F: FnMut(&PileupSite, &[&LiveRead]),
{
    let records_per_sample: Vec<Vec<LiveRead>> = bams
        .iter_mut()
        .enumerate()
        .map(|(i, b)| {
            let mut v = std::mem::take(&mut b.records_buf);
            for lr in v.iter_mut() { lr.sample_idx = i; }
            v
        })
        .collect();
    mpileup_engine_from_records(records_per_sample, min_mq, min_bq, skip_indels, None, emit)
}

/// Engine entry point that takes owned records directly — used by the
/// chunk-parallel runner so each rayon worker can feed its own slice without
/// constructing a synthetic BamReader. Pass `pos_filter` to skip ref-only
/// positions; the engine then jumps cur_pos to the next interesting site.
pub fn mpileup_engine_from_records<F>(
    records_per_sample: Vec<Vec<LiveRead>>,
    min_mq: u8,
    min_bq: u8,
    skip_indels: bool,
    pos_filter: Option<&crate::bam::pos_filter::InterestingMap>,
    emit: &mut F,
) -> Result<()>
where F: FnMut(&PileupSite, &[&LiveRead]),
{
    let n_samples = records_per_sample.len();
    let mut iters: Vec<std::vec::IntoIter<LiveRead>> = records_per_sample
        .into_iter()
        .enumerate()
        .map(|(i, mut v)| {
            for lr in v.iter_mut() { lr.sample_idx = i; }
            v.into_iter()
        })
        .collect();
    let mut next: Vec<Option<LiveRead>> = iters.iter_mut().map(|it| next_aligned_live(it, min_mq)).collect();

    let mut live: VecDeque<LiveRead> = VecDeque::new();
    let mut cur_ref: i32 = -1;
    let mut cur_pos: u32 = 0;

    loop {
        let mut any_active = !live.is_empty();
        for r in &next { if r.is_some() { any_active = true; break; } }
        if !any_active { break; }

        if live.is_empty() {
            let mut pick: Option<(usize, i32, u32)> = None;
            for (i, r) in next.iter().enumerate() {
                if let Some(lr) = r {
                    let cand = (i, lr.ref_id as i32, lr.ref_start);
                    pick = Some(match pick {
                        None => cand,
                        Some(prev) => if (cand.1, cand.2) < (prev.1, prev.2) { cand } else { prev },
                    });
                }
            }
            let Some((idx, rid, p)) = pick else { break; };
            cur_ref = rid;
            cur_pos = p;
            let lr = next[idx].take().unwrap();
            live.push_back(lr);
            next[idx] = next_aligned_live(&mut iters[idx], min_mq);
            continue;
        }

        // Skip ref-only positions when a filter is provided.
        if let Some(filter) = pos_filter {
            if cur_ref >= 0 {
                if let Some(next_pos) = filter.next_at_or_after(cur_ref as usize, cur_pos) {
                    if next_pos > cur_pos {
                        cur_pos = next_pos;
                        // Drop reads that no longer cover cur_pos after the jump.
                        live.retain(|lr| lr.ref_end() > cur_pos);
                    }
                } else {
                    // No more interesting positions on this ref — drain live and move on.
                    live.clear();
                    continue;
                }
            }
        }

        for i in 0..n_samples {
            while let Some(lr) = &next[i] {
                if lr.ref_id as i32 == cur_ref && lr.ref_start <= cur_pos {
                    let taken = next[i].take().unwrap();
                    live.push_back(taken);
                    next[i] = next_aligned_live(&mut iters[i], min_mq);
                } else { break; }
            }
        }

        let mut site = PileupSite { ref_id: cur_ref as usize, pos: cur_pos, per_sample: vec![SampleSiteCounts::default(); n_samples] };
        let mut keep: VecDeque<LiveRead> = VecDeque::with_capacity(live.len());
        while let Some(lr) = live.pop_front() {
            if lr.ref_end() <= cur_pos { continue; }
            if let Some((b, q)) = lr.base_at(cur_pos) {
                if q >= min_bq {
                    if let Some(i) = base_to_idx(b) {
                        let s = &mut site.per_sample[lr.sample_idx];
                        s.base_counts[i] += 1;
                        s.base_quals[i] += q as u32;
                        s.mq_sum += lr.mapq as u64;
                        s.depth += 1;
                    }
                }
            }
            if !skip_indels {
                if let Some(ev) = lr.indel_after(cur_pos) {
                    let s = &mut site.per_sample[lr.sample_idx];
                    match ev {
                        IndelEvent::Ins(seq) => {
                            let key = String::from_utf8_lossy(&seq).into_owned();
                            if let Some(e) = s.ins_alleles.iter_mut().find(|(k, _)| *k == key) { e.1 += 1; } else { s.ins_alleles.push((key, 1)); }
                        }
                        IndelEvent::Del(l) => {
                            if let Some(e) = s.del_alleles.iter_mut().find(|(k, _)| *k == l) { e.1 += 1; } else { s.del_alleles.push((l, 1)); }
                        }
                    }
                }
            }
            keep.push_back(lr);
        }
        live = keep;
        if site.total_depth() > 0 {
            let overlapping: Vec<&LiveRead> = live.iter().filter(|lr| lr.ref_start <= cur_pos && lr.ref_end() > cur_pos).collect();
            emit(&site, &overlapping);
        }
        cur_pos += 1;
    }
    Ok(())
}

fn next_aligned_live(iter: &mut std::vec::IntoIter<LiveRead>, min_mq: u8) -> Option<LiveRead> {
    for lr in iter.by_ref() {
        if lr.skip_by_flags() || lr.mapq < min_mq { continue; }
        return Some(lr);
    }
    None
}

pub fn build_live_from_bam(rec: &bam::Record, sample_idx: usize) -> Option<LiveRead> {
    let start = usize::from(rec.alignment_start().transpose().ok().flatten()?) as u32 - 1;
    let mapq: u8 = rec.mapping_quality().map(|m| u8::from(m)).unwrap_or(0);
    let ref_id = rec.reference_sequence_id().transpose().ok().flatten().unwrap_or(0);
    let flags = u16::from(rec.flags());
    let seq: Vec<u8> = rec.sequence().iter().collect();
    let mut qual: Vec<u8> = rec.quality_scores().iter().collect();
    let mut cigar_pairs: Vec<(Kind, u32)> = Vec::with_capacity(8);
    for op in rec.cigar().iter() {
        let op = op.ok()?;
        cigar_pairs.push((op.kind(), op.len() as u32));
    }
    crate::bam::baq::apply_baq_capping(&seq, &mut qual, &cigar_pairs, 30);
    let ref_end_cached = compute_ref_end(start, &cigar_pairs);
    Some(LiveRead { seq, qual, cigar_pairs, ref_start: start, ref_end_cached, ref_id, mapq, sample_idx, flags })
}

pub fn build_live_from_cram(rb: &noodles_sam::alignment::RecordBuf, sample_idx: usize) -> Option<LiveRead> {
    use noodles_sam::alignment::record::cigar::Op;
    let start = usize::from(rb.alignment_start()?) as u32 - 1;
    let mapq: u8 = rb.mapping_quality().map(|m| u8::from(m)).unwrap_or(0);
    let ref_id = rb.reference_sequence_id().unwrap_or(0);
    let flags = u16::from(rb.flags());
    let seq: Vec<u8> = rb.sequence().as_ref().to_vec();
    let mut qual: Vec<u8> = rb.quality_scores().as_ref().to_vec();
    let cigar_ops: &[Op] = rb.cigar().as_ref();
    let cigar_pairs: Vec<(Kind, u32)> = cigar_ops.iter().map(|op| (op.kind(), op.len() as u32)).collect();
    if cigar_pairs.is_empty() { return None; }
    crate::bam::baq::apply_baq_capping(&seq, &mut qual, &cigar_pairs, 30);
    let ref_end_cached = compute_ref_end(start, &cigar_pairs);
    Some(LiveRead { seq, qual, cigar_pairs, ref_start: start, ref_end_cached, ref_id, mapq, sample_idx, flags })
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

pub struct Pileup;
