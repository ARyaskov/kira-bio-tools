//! Local assembly of an active region (`mpileup --assemble`).
//!
//! The reads spanning the window are threaded through a de Bruijn graph
//! together with the reference; the best-supported source-to-sink paths are
//! the candidate haplotypes. Each haplotype is aligned globally to the
//! reference window (Gotoh, gaps pushed left) to read off its indel events,
//! and every spanning read is scored against every haplotype with the
//! pair-HMM. Genotype likelihoods then follow a diploid haplotype model: for
//! an event, a read's likelihood under the best haplotype carrying it and
//! under the best one without it stand in for the two alleles, so SNPs and
//! other indels on the same haplotype no longer count against a carrier.

use crate::bam::pileup::LiveRead;
use crate::call::haplotype::{AssembledCall, Indel, indel_alleles};
use crate::call::pairhmm::read_vs_hap_loglik;
use crate::fx::FxHashMap;

#[derive(Clone, Debug)]
pub struct AsmParams {
    /// k-mer sizes tried in order until the reference k-mers are unique and
    /// an alternative haplotype appears (each at most 31).
    pub ks: Vec<usize>,
    /// Reads that must support a non-reference edge.
    pub min_kmer: u32,
    /// Base quality below which a k-mer is not counted.
    pub min_bq: u8,
    /// Alternative haplotypes kept per window.
    pub max_haps: usize,
    /// Reads whose best haplotype must carry an event for it to be called.
    pub min_support: u32,
    /// Longest indel event reported.
    pub max_indel: usize,
}

impl Default for AsmParams {
    fn default() -> Self {
        Self {
            ks: vec![21, 27, 31],
            min_kmer: 2,
            min_bq: 10,
            max_haps: 8,
            min_support: 3,
            max_indel: 50,
        }
    }
}

/// Bound on graph walks per window, against bubble blow-up in repeats.
const MAX_STEPS: usize = 20_000;

fn enc(b: u8) -> Option<u64> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

struct Edge {
    to: usize,
    count: u32,
    on_ref: bool,
}

struct Node {
    out: Vec<Edge>,
    last: u8,
}

/// Read-threading de Bruijn graph of one window.
pub struct Graph {
    nodes: Vec<Node>,
    index: FxHashMap<u64, usize>,
    source: usize,
    sink: usize,
    /// Bases of the source k-mer.
    head: Vec<u8>,
}

/// 2-bit k-mers of `seq` (`None` where a base is not ACGT or a quality is
/// below `min_bq`).
fn kmers(seq: &[u8], quals: Option<&[u8]>, k: usize, min_bq: u8) -> Vec<Option<u64>> {
    if seq.len() < k {
        return Vec::new();
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut out = Vec::with_capacity(seq.len() - k + 1);
    let mut code = 0u64;
    // Length of the valid run ending at the current base.
    let mut run = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        let ok = quals.is_none_or(|q| q.get(i).copied().unwrap_or(0) >= min_bq);
        match enc(b) {
            Some(c) if ok => {
                code = ((code << 2) | c) & mask;
                run += 1;
            }
            _ => run = 0,
        }
        if i + 1 >= k {
            out.push((run >= k).then_some(code));
        }
    }
    out
}

impl Graph {
    fn node(&mut self, kmer: u64, last: u8) -> usize {
        if let Some(&i) = self.index.get(&kmer) {
            return i;
        }
        let i = self.nodes.len();
        self.nodes.push(Node {
            out: Vec::new(),
            last,
        });
        self.index.insert(kmer, i);
        i
    }

    fn edge(&mut self, a: usize, b: usize, count: u32, on_ref: bool) {
        let out = &mut self.nodes[a].out;
        match out.iter_mut().find(|e| e.to == b) {
            Some(e) => {
                e.count += count;
                e.on_ref |= on_ref;
            }
            None => out.push(Edge {
                to: b,
                count,
                on_ref,
            }),
        }
    }

    /// The graph of `refw` and the read windows; `None` when the reference
    /// path repeats a k-mer or contains a non-ACGT base.
    pub fn build(k: usize, refw: &[u8], reads: &[(&[u8], &[u8])], min_bq: u8) -> Option<Self> {
        if !(8..=31).contains(&k) || refw.len() < k + 2 {
            return None;
        }
        let ref_kmers = kmers(refw, None, k, 0);
        if ref_kmers.iter().any(Option::is_none) {
            return None;
        }
        let mut g = Self {
            nodes: Vec::new(),
            index: FxHashMap::default(),
            source: 0,
            sink: 0,
            head: refw[..k].to_vec(),
        };
        let mut prev: Option<usize> = None;
        for (i, km) in ref_kmers.iter().enumerate() {
            let km = km.unwrap();
            if g.index.contains_key(&km) {
                return None;
            }
            let n = g.node(km, refw[i + k - 1]);
            if let Some(p) = prev {
                g.edge(p, n, 0, true);
            }
            prev = Some(n);
        }
        g.sink = prev?;
        for (seq, quals) in reads {
            let ks = kmers(seq, Some(quals), k, min_bq);
            let mut prev: Option<usize> = None;
            for (i, km) in ks.iter().enumerate() {
                let Some(km) = km else {
                    prev = None;
                    continue;
                };
                let n = g.node(*km, seq[i + k - 1]);
                if let Some(p) = prev {
                    g.edge(p, n, 1, false);
                }
                prev = Some(n);
            }
        }
        Some(g)
    }

    /// Source-to-sink paths as sequences with their support (the smallest
    /// read count on a non-reference edge; `u32::MAX` for the reference).
    /// Reference edges are always walkable, other edges need `min_kmer`
    /// reads; a k-mer is used at most once per path.
    pub fn haplotypes(&self, min_kmer: u32, max_len: usize) -> Vec<(Vec<u8>, u32)> {
        let mut out: Vec<(Vec<u8>, u32)> = Vec::new();
        let mut seq = self.head.clone();
        let mut visited = vec![false; self.nodes.len()];
        let mut steps = 0usize;
        // Explicit stack of (node, next edge index, support on entry).
        let mut stack: Vec<(usize, usize, u32)> = vec![(self.source, 0, u32::MAX)];
        visited[self.source] = true;
        if self.source == self.sink {
            out.push((seq.clone(), u32::MAX));
        }
        while let Some(top) = stack.last_mut() {
            let (node, ei, support) = *top;
            let edges = &self.nodes[node].out;
            if ei >= edges.len() || steps >= MAX_STEPS || seq.len() >= max_len {
                visited[node] = false;
                stack.pop();
                if !stack.is_empty() {
                    seq.pop();
                }
                continue;
            }
            top.1 += 1;
            let e = &edges[ei];
            if visited[e.to] || (!e.on_ref && e.count < min_kmer) {
                continue;
            }
            steps += 1;
            let sup = if e.on_ref {
                support
            } else {
                support.min(e.count)
            };
            visited[e.to] = true;
            seq.push(self.nodes[e.to].last);
            stack.push((e.to, 0, sup));
            if e.to == self.sink {
                out.push((seq.clone(), sup));
            }
        }
        out
    }
}

/// Indel events of `hap` relative to `refw` from a global affine-gap
/// alignment (mismatch 4, gap open 6, extend 1), gaps left-normalised,
/// anchored inside the window and at most `max_indel` long.
pub fn indel_events(refw: &[u8], hap: &[u8], max_indel: usize) -> Vec<Indel> {
    const MM: i32 = 4;
    const GO: i32 = 6;
    const GE: i32 = 1;
    const INF: i32 = i32::MAX / 4;
    let (n, m) = (hap.len(), refw.len());
    let w = m + 1;
    let idx = |i: usize, j: usize| i * w + j;
    let mut mm = vec![INF; (n + 1) * w];
    let mut x = vec![INF; (n + 1) * w]; // gap in reference: insertion
    let mut y = vec![INF; (n + 1) * w]; // gap in haplotype: deletion
    mm[0] = 0;
    for i in 1..=n {
        x[idx(i, 0)] = GO + (i as i32 - 1) * GE;
    }
    for j in 1..=m {
        y[idx(0, j)] = GO + (j as i32 - 1) * GE;
    }
    for i in 1..=n {
        for j in 1..=m {
            let d = idx(i - 1, j - 1);
            let best = mm[d].min(x[d]).min(y[d]);
            let sub = if hap[i - 1].eq_ignore_ascii_case(&refw[j - 1]) {
                0
            } else {
                MM
            };
            mm[idx(i, j)] = best.saturating_add(sub);
            let u = idx(i - 1, j);
            x[idx(i, j)] = (mm[u].saturating_add(GO))
                .min(x[u].saturating_add(GE))
                .min(y[u].saturating_add(GO));
            let l = idx(i, j - 1);
            y[idx(i, j)] = (mm[l].saturating_add(GO))
                .min(y[l].saturating_add(GE))
                .min(x[l].saturating_add(GO));
        }
    }
    // Traceback; ties prefer the diagonal, then insertion, so gaps end up
    // as far left as the scores allow.
    let (mut i, mut j) = (n, m);
    let e = idx(n, m);
    let mut state = if mm[e] <= x[e] && mm[e] <= y[e] {
        0
    } else if x[e] <= y[e] {
        1
    } else {
        2
    };
    let mut ins: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut del: Vec<(usize, usize)> = Vec::new();
    while i > 0 || j > 0 {
        match state {
            0 => {
                let d = idx(i - 1, j - 1);
                let here = mm[idx(i, j)]
                    - if hap[i - 1].eq_ignore_ascii_case(&refw[j - 1]) {
                        0
                    } else {
                        MM
                    };
                state = if mm[d] == here {
                    0
                } else if x[d] == here {
                    1
                } else {
                    2
                };
                i -= 1;
                j -= 1;
            }
            1 => {
                // hap[i-1] inserted before refw[j].
                match ins.last_mut() {
                    Some((o, b)) if *o == j => b.insert(0, hap[i - 1]),
                    _ => ins.push((j, vec![hap[i - 1]])),
                }
                let here = x[idx(i, j)];
                let u = idx(i - 1, j);
                state = if i > 1 && x[u] + GE == here {
                    1
                } else if mm[u] + GO == here || i == 1 {
                    0
                } else {
                    2
                };
                i -= 1;
            }
            _ => {
                // refw[j-1] deleted.
                match del.last_mut() {
                    Some((o, len)) if *o == j => {
                        *o -= 1;
                        *len += 1;
                    }
                    _ => del.push((j - 1, 1)),
                }
                let here = y[idx(i, j)];
                let l = idx(i, j - 1);
                state = if j > 1 && y[l] + GE == here {
                    2
                } else if mm[l] + GO == here || j == 1 {
                    0
                } else {
                    1
                };
                j -= 1;
            }
        }
    }
    let mut out = Vec::new();
    for (mut o, mut bases) in ins {
        while o > 0 && refw[o - 1].eq_ignore_ascii_case(bases.last().unwrap()) {
            bases.rotate_right(1);
            o -= 1;
        }
        if o >= 1 && o < m && bases.len() <= max_indel {
            out.push(Indel {
                is_ins: true,
                ref_off: o,
                len: bases.len(),
                bases,
            });
        }
    }
    for (mut o, len) in del {
        while o > 0 && refw[o - 1].eq_ignore_ascii_case(&refw[o + len - 1]) {
            o -= 1;
        }
        if o >= 1 && o + len < m && len <= max_indel {
            out.push(Indel {
                is_ins: false,
                ref_off: o,
                len,
                bases: Vec::new(),
            });
        }
    }
    out.sort_by_key(|e| (e.ref_off, e.is_ins, e.len));
    out
}

/// Alternative haplotypes of the window that carry at least one indel,
/// best supported first, at most `max_haps`.
pub fn candidate_haplotypes(refw: &[u8], reads: &[(&[u8], &[u8])], p: &AsmParams) -> Vec<Vec<u8>> {
    for &k in &p.ks {
        let Some(g) = Graph::build(k, refw, reads, p.min_bq) else {
            continue;
        };
        let mut haps: Vec<(Vec<u8>, u32)> = g
            .haplotypes(p.min_kmer, refw.len() + p.max_indel)
            .into_iter()
            .filter(|(h, s)| {
                *s != u32::MAX && h != refw && !indel_events(refw, h, p.max_indel).is_empty()
            })
            .collect();
        if haps.is_empty() {
            continue;
        }
        haps.sort_by(|a, b| {
            b.1.cmp(&a.1)
                .then(
                    a.0.len()
                        .abs_diff(refw.len())
                        .cmp(&b.0.len().abs_diff(refw.len())),
                )
                .then(a.0.cmp(&b.0))
        });
        haps.dedup_by(|a, b| a.0 == b.0);
        haps.truncate(p.max_haps);
        return haps.into_iter().map(|(h, _)| h).collect();
    }
    Vec::new()
}

/// One call per distinct indel event of the haplotypes, with every read's
/// likelihood under the best carrier and the best non-carrier haplotype.
pub fn calls_from_haplotypes(
    refw: &[u8],
    win_lo: u32,
    alts: &[Vec<u8>],
    reads: &[(Vec<u8>, Vec<u8>, usize)],
    p: &AsmParams,
) -> Vec<AssembledCall> {
    let mut haps: Vec<&[u8]> = Vec::with_capacity(alts.len() + 1);
    haps.push(refw);
    haps.extend(alts.iter().map(Vec::as_slice));
    // Distinct events and the haplotypes carrying them.
    let mut events: Vec<(Indel, Vec<bool>)> = Vec::new();
    for (h, hap) in haps.iter().enumerate().skip(1) {
        for ev in indel_events(refw, hap, p.max_indel) {
            match events.iter_mut().find(|(e, _)| *e == ev) {
                Some((_, c)) => c[h] = true,
                None => {
                    let mut c = vec![false; haps.len()];
                    c[h] = true;
                    events.push((ev, c));
                }
            }
        }
    }
    if events.is_empty() {
        return Vec::new();
    }
    // Read-vs-haplotype log-likelihoods.
    let ll: Vec<Vec<f64>> = reads
        .iter()
        .map(|(b, q, _)| haps.iter().map(|h| read_vs_hap_loglik(b, q, h)).collect())
        .collect();
    let total = reads.len() as u32;
    let mut out = Vec::new();
    for (ev, carriers) in &events {
        let Some((anchor, ref_str, alt_str)) = indel_alleles(refw, ev) else {
            continue;
        };
        let mut read_ll = Vec::with_capacity(reads.len());
        let mut support = 0u32;
        for (r, (_, _, sample)) in reads.iter().enumerate() {
            let best = |carry: bool| -> f64 {
                ll[r]
                    .iter()
                    .zip(carriers)
                    .filter(|&(_, &c)| c == carry)
                    .map(|(&l, _)| l)
                    .fold(f64::NEG_INFINITY, f64::max)
            };
            let (lref, lalt) = (best(false), best(true));
            if !lref.is_finite() || !lalt.is_finite() {
                continue;
            }
            if lalt > lref {
                support += 1;
            }
            read_ll.push((*sample, lref, lalt));
        }
        if support < p.min_support {
            continue;
        }
        out.push(AssembledCall {
            pos1: win_lo as u64 + anchor as u64 + 1,
            ref_str,
            alt_str,
            support,
            total,
            read_ll,
        });
    }
    out
}

/// Assemble the window `[win_lo, win_lo + refw.len())` from the reads that
/// span it; empty when no supported indel haplotype exists.
pub fn assemble_window(
    reads: &[LiveRead],
    win_lo: u32,
    refw: &[u8],
    p: &AsmParams,
) -> Vec<AssembledCall> {
    let hi = win_lo + refw.len() as u32;
    let sub: Vec<(Vec<u8>, Vec<u8>, usize)> = reads
        .iter()
        .filter_map(|lr| {
            let (b, q) = lr.query_window_qual(win_lo, hi)?;
            (b.len() >= 12).then_some((b, q, lr.sample_idx))
        })
        .collect();
    if (sub.len() as u32) < p.min_support {
        return Vec::new();
    }
    let views: Vec<(&[u8], &[u8])> = sub
        .iter()
        .map(|(b, q, _)| (b.as_slice(), q.as_slice()))
        .collect();
    let alts = candidate_haplotypes(refw, &views, p);
    if alts.is_empty() {
        return Vec::new();
    }
    calls_from_haplotypes(refw, win_lo, &alts, &sub, p)
}

#[cfg(test)]
#[path = "../../tests/unit/call_assembly.rs"]
mod tests;
