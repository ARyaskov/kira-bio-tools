use noodles_sam::alignment::record::cigar::op::Kind;

use super::pileup::{IndelEvent, LiveRead};
use crate::align::glocal::glocal_loglik_with;
use crate::align::{GlocalParams, GlocalScratch, encode_nt};

pub struct IndelOptions {
    pub min_support: usize,
    pub min_fraction: f64,
    pub per_sample: bool,
    pub window: usize,
    pub open_quality: i32,
    pub extension_quality: i32,
    pub tandem_quality: i32,
    pub bias: f64,
    pub minimum_quality: i32,
    pub ambiguous: u8,
}

pub struct IndelRead {
    pub read: usize,
    pub query: usize,
    pub gap: i32,
    pub allele: usize,
    pub quality: i32,
    pub sequence_quality: i32,
}

pub struct IndelCalls {
    pub types: Vec<i32>,
    pub insertions: Vec<Vec<u8>>,
    pub observations: Vec<IndelRead>,
    pub extent: usize,
    pub max_support: i32,
    pub max_fraction: f32,
    pub minimum_quality: i32,
    pub ambiguous: u8,
}

fn aligned(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

pub(super) fn bias_position(read: &LiveRead, query: usize, gap: i32) -> (usize, usize) {
    let (mut clip, mut distance, mut at_left) = (0i64, -1i64, true);
    let (mut position, mut length) = (query as i64, read.seq().len() as i64);
    for &(kind, size) in &read.cigar_pairs {
        if kind == Kind::SoftClip {
            length -= size as i64;
            if at_left {
                clip += size as i64;
                position -= clip;
                distance = position;
            } else {
                let right = read.seq().len() as i64 - size as i64 - query as i64;
                if distance < 0 || distance > right {
                    clip = size as i64;
                    distance = right;
                }
            }
        } else if kind != Kind::HardClip {
            at_left = false;
        }
    }
    if gap > 0 && length - (position + gap as i64) < position {
        position += gap as i64 - 1;
    }
    let position = ((position as f64 / (length + 1) as f64) * 100.0) as usize;
    let clip = if clip > 0 {
        ((15.0 * clip as f64 / (distance + 1) as f64) as usize).min(99)
    } else {
        0
    };
    (position.min(99), clip)
}

fn query_position(read: &LiveRead, position: usize, left: bool) -> (usize, usize) {
    let (mut reference, mut query, mut last) = (read.ref_start as usize, 0, 0);
    for &(kind, length) in &read.cigar_pairs {
        let length = length as usize;
        if aligned(kind) {
            if read.ref_start as usize > position {
                return (query, read.ref_start as usize);
            }
            if reference + length > position {
                return (query + position - reference, position);
            }
            reference += length;
            query += length;
            last = query;
        } else if matches!(kind, Kind::Insertion | Kind::SoftClip) {
            query += length;
        } else if matches!(kind, Kind::Deletion | Kind::Skip) {
            if reference + length > position {
                return (query, if left { reference } else { reference + length });
            }
            reference += length;
        }
    }
    (last, reference)
}

fn gap(read: &LiveRead, position: usize) -> i32 {
    match read.indel_after(position as u32) {
        Some(IndelEvent::Ins(sequence)) => sequence.len() as i32,
        Some(IndelEvent::Del(length)) => -(length as i32),
        None => 0,
    }
}

fn consensus(
    reference: &[u8],
    reads: &[LiveRead],
    indices: &[usize],
    left: usize,
    right: usize,
) -> Vec<u8> {
    let mut counts = vec![(0u32, 0u32); right - left];
    for &index in indices {
        let read = &reads[index];
        let (mut r, mut q) = (read.ref_start as usize, 0);
        for &(kind, length) in &read.cigar_pairs {
            let length = length as usize;
            if aligned(kind) {
                for at in left.max(r)..right.min(r + length) {
                    let query = q + at - r;
                    if let Some(&base) = read.seq().get(query) {
                        if base.eq_ignore_ascii_case(&reference[at]) {
                            counts[at - left].0 += 1;
                        } else {
                            counts[at - left].1 += 1;
                        }
                    }
                }
                r += length;
                q += length;
            } else if matches!(kind, Kind::Insertion | Kind::SoftClip) {
                q += length;
            } else if matches!(kind, Kind::Deletion | Kind::Skip) {
                r += length;
            }
            if r > right {
                break;
            }
        }
    }
    let mut result: Vec<_> = reference[left..right]
        .iter()
        .map(|&base| encode_nt(base))
        .collect();
    let (mut first, mut second) = (None, None);
    for i in 0..counts.len() {
        if first.is_none_or(|k: usize| counts[i].1 >= counts[k].1) {
            second = first;
            first = Some(i);
        } else if second.is_none_or(|k: usize| counts[i].1 >= counts[k].1) {
            second = Some(i);
        }
    }
    for index in [first, second].into_iter().flatten() {
        let (ref_count, alt_count) = counts[index];
        if (ref_count as f64 / (ref_count + alt_count) as f64) < 0.7 || ref_count + alt_count == 0 {
            result[index] = 4;
        }
    }
    result
}

fn repeat_regions(sequence: &[u8]) -> Vec<(usize, usize, usize)> {
    let mut result: Vec<(usize, usize, usize)> = Vec::new();
    let mut word = 0u32;
    for (position, &base) in sequence.iter().enumerate() {
        word = (word << 2) | base as u32;
        for step in 0..8 {
            let length = if position < 15 { step + 1 } else { 8 - step };
            if position < 15 && (length > 7 || position < 2 * length - 1) {
                continue;
            }
            let mask = (1u32 << (2 * length)) - 1;
            if word & mask != (word >> (2 * length)) & mask {
                continue;
            }
            let start = position + 1 - 2 * length;
            if result
                .last()
                .is_some_and(|&(s, e, _)| s <= start && e >= position)
            {
                if position >= 15 {
                    break;
                } else {
                    continue;
                }
            }
            let mut end = position + 1;
            while end < sequence.len() && sequence[end - length] == sequence[end] {
                end += 1;
            }
            while result
                .last()
                .is_some_and(|&(s, e, _)| e >= start && s >= start)
            {
                result.pop();
            }
            result.push((start, end - 1, length));
            if position >= 15 {
                break;
            }
        }
    }
    result
}

fn sequence_quality(options: &IndelOptions, length: i32, run: usize) -> i32 {
    let length = length.abs();
    let open = options.open_quality + options.extension_quality * (length - 1);
    let homopolymer = if run >= 3 {
        (options.tandem_quality as f64 * length as f64 / run as f64 + 0.499) as i32
    } else {
        1000
    };
    open.min(homopolymer)
}

pub fn prepare(
    reads: &[LiveRead],
    samples: &[usize],
    position: usize,
    reference: &[u8],
    options: &IndelOptions,
) -> Option<IndelCalls> {
    if position + 1 >= reference.len() {
        return None;
    }
    let groups: Vec<Vec<usize>> = samples
        .iter()
        .map(|sample| {
            reads
                .iter()
                .enumerate()
                .filter(|(_, read)| read.sample_idx == *sample)
                .map(|(i, _)| i)
                .collect()
        })
        .collect();
    let indices: Vec<_> = groups.iter().flatten().copied().collect();
    let total = indices.len();
    if total == 0 {
        return None;
    }
    let mut types = vec![0];
    let (mut support, mut max_support, mut max_fraction, mut sample_support) =
        (0, 0, 0.0f64, false);
    for group in &groups {
        let mut count = 0;
        for &index in group {
            let length = gap(&reads[index], position);
            if length != 0 {
                count += 1;
                types.push(length);
            }
        }
        let fraction = count as f64 / group.len() as f64;
        sample_support |= count >= options.min_support && fraction >= options.min_fraction;
        if count > max_support && fraction > 0.0 {
            max_support = count;
            max_fraction = fraction;
        }
        support += count;
    }
    types.sort_unstable();
    types.dedup();
    if types.len() < 2
        || types.len() >= 64
        || if options.per_sample {
            !sample_support
        } else {
            support < options.min_support || (support as f64 / total as f64) < options.min_fraction
        }
    {
        return None;
    }
    let max_read = indices
        .iter()
        .map(|&i| reads[i].seq().len())
        .max()
        .unwrap_or(0);
    let end = (position + (2 * options.window).min(max_read)).min(reference.len());
    if reference[position..end]
        .iter()
        .filter(|&&base| base == b'N')
        .count()
        * 2
        > end - position
    {
        return None;
    }
    let ref_type = types.iter().position(|&length| length == 0).unwrap();
    let left = position.saturating_sub(options.window);
    let right = (position + options.window + (-types[0]).max(0) as usize).min(reference.len());
    let mut run_left = position + 1;
    let mut run_right = position + 2;
    let run_base = reference[position + 1];
    if run_base != b'N' {
        while run_left > 0 && reference[run_left - 1] == run_base {
            run_left -= 1;
        }
        while run_right < reference.len() && reference[run_right] == run_base {
            run_right += 1;
        }
    }
    let run = run_right - run_left;
    let mut insertions = vec![Vec::new(); types.len()];
    for (t, &length) in types.iter().enumerate() {
        if length <= 0 {
            continue;
        }
        let mut counts = vec![[0u32; 5]; length as usize];
        for &index in &indices {
            if let Some(IndelEvent::Ins(sequence)) = reads[index].indel_after(position as u32)
                && sequence.len() == length as usize
            {
                for (j, &base) in sequence.iter().enumerate() {
                    counts[j][encode_nt(base) as usize] += 1;
                }
            }
        }
        for row in counts {
            let best = (0..5)
                .max_by_key(|&i| (row[i], std::cmp::Reverse(i)))
                .unwrap();
            insertions[t].push(if row[best] > 0 { best as u8 } else { 4 });
        }
    }
    if insertions.iter().any(|sequence| sequence.contains(&4)) {
        return None;
    }
    let mut extent = 0;
    for (t, &length) in types.iter().enumerate() {
        if length == 0 {
            continue;
        }
        let length = length.unsigned_abs() as usize;
        let (mut score, mut best) = (0, 0);
        for at in position + 1..reference.len() {
            let offset = at - position - 1;
            let base = if types[t] > 0 {
                b"ACGTN"[insertions[t][offset % length] as usize]
            } else {
                reference[position + 1 + offset % length]
            };
            score += if base == reference[at] { 1 } else { -10 };
            if score < 0 {
                break;
            }
            if score > best {
                best = score;
                extent = extent.max(at - position);
            }
        }
    }
    let consensuses: Vec<_> = groups
        .iter()
        .map(|group| consensus(reference, reads, group, left, right))
        .collect();
    let mut scores = vec![vec![0i32; types.len()]; total];
    let mut scratch = GlocalScratch::default();
    for (t, &length) in types.iter().enumerate() {
        let mut row = 0;
        for (sample, group) in groups.iter().enumerate() {
            let offset = position + 1 - left;
            let mut haplotype = consensuses[sample][..offset].to_vec();
            if length > 0 {
                haplotype.extend_from_slice(&insertions[t]);
            }
            let after = offset + (-length).max(0) as usize;
            if after < consensuses[sample].len() {
                haplotype.extend_from_slice(&consensuses[sample][after..]);
            }
            haplotype.resize(
                right - left + 2 + 2 * types.last().copied().unwrap_or(0).max(-types[0]) as usize,
                4,
            );
            for &index in group {
                let read = &reads[index];
                if read.cigar_pairs.iter().any(|&(kind, _)| kind == Kind::Skip) {
                    row += 1;
                    continue;
                }
                let long = read.seq().len() > 1000;
                let left2 = left
                    + if long && position - left >= options.window {
                        options.window / 2
                    } else {
                        0
                    };
                let right2 = right
                    - if long && right - position >= options.window {
                        options.window / 2
                    } else {
                        0
                    };
                let (qbeg, mut tbeg) = query_position(read, left2, false);
                let (qposition, _) = query_position(read, position, false);
                let (qend, tend) = query_position(read, right2, true);
                if length < 0 {
                    tbeg = tbeg.saturating_sub((-length) as usize).max(left);
                }
                if tend <= tbeg || qend <= qbeg || qend > read.seq().len() {
                    scores[row][t] = 0xffffff;
                    row += 1;
                    continue;
                }
                let query: Vec<_> = read.seq()[qbeg..qend]
                    .iter()
                    .map(|&base| encode_nt(base))
                    .collect();
                let quality: Vec<_> = (qbeg..qend)
                    .map(|q| {
                        read.qual()[q]
                            .saturating_add(read.baq_delta.as_ref().map_or(0, |delta| delta[q]))
                            .clamp(7, 30)
                    })
                    .collect();
                let start = tbeg - left;
                let size = tend - tbeg + length.unsigned_abs() as usize;
                if start + size > haplotype.len() {
                    scores[row][t] = 0xffffff;
                    row += 1;
                    continue;
                }
                let segment = &haplotype[start..start + size];
                let params = GlocalParams {
                    d: if long { 1e-3 } else { 1e-4 },
                    e: if long { 1e-1 } else { 1e-2 },
                    bw: length.unsigned_abs() as usize + 3,
                };
                let likelihood =
                    glocal_loglik_with(&mut scratch, segment, &query, Some(&quality), &params);
                let phred = (-4.343 * likelihood + 0.499) as i32;
                let normalized = ((100.0 * phred as f64 / query.len() as f64 + 0.499) as i32 as f64
                    * options.bias) as i32;
                let mut repeat_score = 0;
                let qpos = qposition.saturating_sub(qbeg);
                for (start, end, unit) in repeat_regions(segment) {
                    if start <= qpos && end >= qpos {
                        repeat_score += (end - start) / unit;
                        if start + tbeg <= read.ref_start as usize
                            || end + tbeg >= read.ref_end() as usize - 1
                        {
                            repeat_score += 2 * (end - start);
                        }
                    }
                }
                let adjusted =
                    (normalized.min(255) as f64 * 0.8 + (repeat_score * 2) as f64) as i32;
                scores[row][t] = (phred << 8) | adjusted.min(255);
                row += 1;
            }
        }
    }
    let mut sums = vec![0i32; types.len()];
    let mut observations = Vec::new();
    for (row, &index) in indices.iter().enumerate() {
        let mut order: Vec<_> = (0..types.len()).collect();
        order.sort_by_key(|&t| (scores[row][t], t));
        let best = order[0];
        let other = if best == ref_type { order[1] } else { ref_type };
        let raw = (scores[row][other] >> 8) - (scores[row][best] >> 8);
        let sequence_quality = sequence_quality(
            options,
            types[if best == ref_type { other } else { best }],
            run,
        )
        .min(255);
        let normalized = scores[row][best] & 255;
        let quality = if normalized > 111 {
            0
        } else {
            ((1.0 - normalized as f64 / 111.0) * raw as f64 + 0.499) as i32
        }
        .min(sequence_quality)
        .min(255);
        sums[best] += quality.min(sequence_quality);
        let (query, _) = query_position(&reads[index], position, false);
        observations.push(IndelRead {
            read: index,
            query: query.min(reads[index].seq().len().saturating_sub(1)),
            gap: gap(&reads[index], position),
            allele: best,
            quality,
            sequence_quality,
        });
    }
    let mut order: Vec<_> = (0..types.len()).collect();
    order.sort_by_key(|&t| std::cmp::Reverse((sums[t], t)));
    order.retain(|&t| t != ref_type);
    order.insert(0, ref_type);
    order.truncate(4);
    for observation in &mut observations {
        observation.allele = order
            .iter()
            .position(|&t| t == observation.allele)
            .unwrap_or(4);
        if observation.allele == 4 {
            observation.quality = 0;
            observation.sequence_quality = 0;
        }
    }
    if !observations
        .iter()
        .any(|observation| observation.allele > 0)
    {
        return None;
    }
    Some(IndelCalls {
        types: order.iter().map(|&t| types[t]).collect(),
        insertions: order.iter().map(|&t| insertions[t].clone()).collect(),
        observations,
        extent,
        max_support: max_support as i32,
        max_fraction: max_fraction as f32,
        minimum_quality: options.minimum_quality,
        ambiguous: options.ambiguous,
    })
}

#[cfg(test)]
#[path = "../../tests/unit/bam_indel.rs"]
mod tests;
