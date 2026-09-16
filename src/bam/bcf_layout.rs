//! bcftools `mpileup` record layout (bam2bcf.c): alleles ordered by quality
//! sum with the unseen `<*>`, PL over that order, I16/QS/MQ0F and the bias
//! statistics (VDB, SGB, Mann-Whitney U-z tests, Fisher strand bias),
//! written in bcftools' INFO/FORMAT order so `bcftools call` and filtering
//! recipes consume the output unchanged. Float arithmetic follows the C
//! types (single precision where bcftools uses `float`) so values agree to
//! the printed digit.

use std::io::Write;

use noodles_sam::alignment::record::cigar::op::Kind;

use crate::bam::errmod::{Drand48, ErrorModel, lgamma};
use crate::bam::mpileup_opts::AnnotateSpec;
use crate::bam::pileup::{IndelEvent, LiveRead, PileupSite, adjusted_baseq};
use crate::util::fmt_g;

const NPOS: usize = 100;
const NQUAL: usize = 60;
const CAP_DIST: i64 = 25;
const CAP_Q: i64 = 60;
const MAX_ALLELES: usize = 5;
const N_NM: usize = 32;

/// Per-site knobs shared by every record of a run.
pub struct LayoutCtx<'a> {
    pub tags: &'a AnnotateSpec,
    pub em: &'a ErrorModel,
    /// The run's `hts_drand48` stream (`--seed`), shared by every site.
    pub rng: &'a std::sync::Mutex<Drand48>,
    pub delta_bq: Option<u8>,
    pub max_bq: Option<u8>,
}

/// One indel allele as kira's candidate detection reports it.
pub struct IndelSpec<'a> {
    pub ref_str: &'a str,
    pub alt_str: &'a str,
    pub is_ins: bool,
    pub ins_seq: &'a [u8],
    pub del_len: u32,
    pub indel_qual: u8,
}

/// `bcf_callret1_t`: one sample's contribution to a site.
struct SampleCall {
    anno: [f64; 16],
    p: Vec<f32>,
    qs: [i64; MAX_ALLELES],
    adf: [i32; MAX_ALLELES],
    adr: [i32; MAX_ALLELES],
    ori_depth: i32,
    mq0: i32,
    scr: i32,
    ref_nm: [i32; N_NM],
    alt_nm: [i32; N_NM],
}

/// The bias-test histograms of `bcf_callaux_t`, pooled over samples.
struct Bias {
    ref_nm: [i32; N_NM],
    alt_nm: [i32; N_NM],
    /// Mismatch observations and their sum for ref/alt reads (INFO/NM).
    nnm: [u32; 2],
    nm: [f32; 2],
    want_nm: bool,
    ref_pos: [i32; NPOS],
    alt_pos: [i32; NPOS],
    ref_bq: [i32; NQUAL],
    alt_bq: [i32; NQUAL],
    ref_mq: [i32; NQUAL],
    alt_mq: [i32; NQUAL],
    fwd_mqs: [i32; NQUAL],
    rev_mqs: [i32; NQUAL],
    ref_scl: [i32; NPOS],
    alt_scl: [i32; NPOS],
}

impl Bias {
    fn new(want_nm: bool) -> Self {
        Self {
            ref_nm: [0; N_NM],
            alt_nm: [0; N_NM],
            nnm: [0; 2],
            nm: [0.0; 2],
            want_nm,
            ref_pos: [0; NPOS],
            alt_pos: [0; NPOS],
            ref_bq: [0; NQUAL],
            alt_bq: [0; NQUAL],
            ref_mq: [0; NQUAL],
            alt_mq: [0; NQUAL],
            fwd_mqs: [0; NQUAL],
            rev_mqs: [0; NQUAL],
            ref_scl: [0; NPOS],
            alt_scl: [0; NPOS],
        }
    }
}

/// `seq_nt16_int[seq_nt16_table[c]]`: A/C/G/T to 0..3, anything else 4.
#[inline]
pub fn ref4_of(c: u8) -> usize {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// `get_position`: distance of query index `qi` from the aligned start and
/// the aligned length, plus the nearest soft clip's length and distance
/// (with bcftools' quirk that a nearer right clip reports a zero length).
fn read_position(lr: &LiveRead, qi: usize) -> (i64, i64, i64, i64) {
    let cig = &lr.cigar_pairs;
    let l_qseq = lr.seq().len() as i64;
    let mut i = 0usize;
    let mut sc_left = 0i64;
    while i < cig.len() {
        match cig[i].0 {
            Kind::HardClip => {}
            Kind::SoftClip => sc_left += cig[i].1 as i64,
            _ => break,
        }
        i += 1;
    }
    let mut sc_right = 0i64;
    let mut j = cig.len() as i64 - 1;
    while j >= i as i64 {
        match cig[j as usize].0 {
            Kind::HardClip => {}
            Kind::SoftClip => sc_right += cig[j as usize].1 as i64,
            _ => break,
        }
        j -= 1;
    }
    let qpos = qi as i64;
    let sc_left_dist = if sc_left > 0 { qpos + 1 - sc_left } else { -1 };
    let sc_right_dist = if sc_right > 0 {
        l_qseq - sc_right - qpos
    } else {
        -1
    };
    let (sc_len, sc_dist) = if sc_left_dist >= 0 {
        if sc_right_dist < 0 || sc_left_dist < sc_right_dist {
            (sc_left, sc_left_dist)
        } else {
            (0, 0)
        }
    } else if sc_right_dist >= 0 {
        (sc_right, sc_right_dist)
    } else {
        (0, 0)
    };
    (
        qpos + 1 - sc_left,
        l_qseq - sc_left - sc_right,
        sc_len,
        sc_dist,
    )
}

/// One observation folded into a sample's call (`bcf_call_glfgen` body).
#[allow(clippy::too_many_arguments)]
fn add_observation(
    call: &mut SampleCall,
    bias: &mut Bias,
    bases: &mut Vec<u16>,
    lr: &LiveRead,
    qi: usize,
    b: usize,
    q_in: i64,
    is_diff: bool,
) {
    let rev = lr.is_reverse();
    let base_q = q_in;
    let inm = if bias.want_nm {
        aux_nm(lr, !is_diff)
    } else {
        None
    };
    if let Some(v) = inm {
        bias.nnm[usize::from(is_diff)] += 1;
        bias.nm[usize::from(is_diff)] += v as f32;
    }
    let mapq0 = lr.mapq as i64;
    if mapq0 == 0 {
        call.mq0 += 1;
    }
    let mapq = mapq0.min(CAP_Q);
    let q = q_in.min(mapq).clamp(4, 63);
    bases.push(((q as u16) << 5) | ((rev as u16) << 4) | (b as u16));
    if b < 4 {
        call.qs[b] += q;
        if rev {
            call.adr[b] += 1
        } else {
            call.adf[b] += 1
        }
    }
    let d = usize::from(is_diff);
    call.anno[(d << 1) | usize::from(rev)] += 1.0;
    let l_qseq = lr.seq().len() as i64;
    let min_dist = (l_qseq - 1 - qi as i64).min(qi as i64).min(CAP_DIST);
    call.anno[4 + (d << 1)] += base_q as f64;
    call.anno[5 + (d << 1)] += (base_q * base_q) as f64;
    call.anno[8 + (d << 1)] += mapq as f64;
    call.anno[9 + (d << 1)] += (mapq * mapq) as f64;
    call.anno[12 + (d << 1)] += min_dist as f64;
    call.anno[13 + (d << 1)] += (min_dist * min_dist) as f64;

    let base_q = base_q.min(59) as usize;
    let mapq = mapq.min(59) as usize;
    let (pos, len, sc_len, sc_dist) = read_position(lr, qi);
    let epos = ((pos as f64 / (len + 1) as f64) * (NPOS - 1) as f64) as i64;
    let epos = epos.clamp(0, NPOS as i64 - 1) as usize;
    let scl = if sc_len > 0 {
        ((15.0 * sc_len as f64 / (sc_dist + 1) as f64) as i64).clamp(0, 99) as usize
    } else {
        0
    };
    if rev {
        bias.rev_mqs[mapq] += 1
    } else {
        bias.fwd_mqs[mapq] += 1
    }
    if is_diff {
        bias.alt_pos[epos] += 1;
        bias.alt_bq[base_q] += 1;
        bias.alt_mq[mapq] += 1;
        bias.alt_scl[scl] += 1;
        if let Some(v) = inm {
            bias.alt_nm[v] += 1;
            call.alt_nm[v] += 1;
        }
    } else {
        bias.ref_pos[epos] += 1;
        bias.ref_bq[base_q] += 1;
        bias.ref_mq[mapq] += 1;
        bias.ref_scl[scl] += 1;
        if let Some(v) = inm {
            bias.ref_nm[v] += 1;
            call.ref_nm[v] += 1;
        }
    }
}

fn new_call() -> SampleCall {
    SampleCall {
        anno: [0.0; 16],
        p: Vec::new(),
        qs: [0; MAX_ALLELES],
        adf: [0; MAX_ALLELES],
        adr: [0; MAX_ALLELES],
        ori_depth: 0,
        mq0: 0,
        scr: 0,
        ref_nm: [0; N_NM],
        alt_nm: [0; N_NM],
    }
}

/// `get_aux_nm`: the read's mismatch bin, counting indels as single events,
/// soft clips as mismatches, and discounting the base itself (two for an
/// alternate base, since de novo SNVs cluster).
fn aux_nm(lr: &LiveRead, is_ref: bool) -> Option<usize> {
    if lr.nm < 0 {
        return None;
    }
    let mut nm = lr.nm as i64;
    for &(k, l) in lr.cigar_pairs.iter() {
        match k {
            Kind::SoftClip => nm += l as i64,
            Kind::Insertion | Kind::Deletion if l > 1 => {
                nm -= l as i64 - 1;
            }
            _ => {}
        }
    }
    nm -= if is_ref { 1 } else { 2 };
    Some(nm.clamp(0, N_NM as i64 - 1) as usize)
}

fn softclip_reads(reads: &[LiveRead], sample: usize) -> i32 {
    reads
        .iter()
        .filter(|r| r.sample_idx == sample && r.has_softclip)
        .count() as i32
}

/// `bcf_call_combine` output.
struct Combined {
    a: [i64; MAX_ALLELES],
    n_alleles: usize,
    unseen: i64,
    qsum: [f32; MAX_ALLELES],
    pl: Vec<i32>,
    anno: [f64; 16],
    ori_depth: i32,
    mq0: i32,
    has_alt: bool,
}

fn combine(calls: &[SampleCall], ref4: Option<usize>) -> Combined {
    let mut qsum = [0f32; MAX_ALLELES];
    for c in calls {
        let mut sum = 0f32;
        for j in 0..4 {
            sum += c.qs[j] as f32;
        }
        if sum != 0.0 {
            for j in 0..4 {
                qsum[j] += c.qs[j] as f32 / sum;
            }
        }
    }
    // Insertion sort of the four base indices by ascending qsum, stable.
    let mut ptr = [0usize, 1, 2, 3, 4];
    for i in 1..4 {
        let mut j = i;
        while j > 0 && qsum[ptr[j]] < qsum[ptr[j - 1]] {
            ptr.swap(j, j - 1);
            j -= 1;
        }
    }
    let mut a = [-1i64; MAX_ALLELES];
    let mut out_qsum = [0f32; MAX_ALLELES];
    let mut unseen = -1i64;
    let r4 = ref4.unwrap_or(0);
    a[0] = r4 as i64;
    let mut j = 1usize;
    let mut i: i64 = 3;
    while i >= 0 {
        let ipos = ptr[i as usize];
        if ipos == r4 {
            out_qsum[0] = qsum[ipos];
        } else {
            if qsum[ipos] == 0.0 {
                break;
            }
            out_qsum[j] = qsum[ipos];
            a[j] = ipos as i64;
            j += 1;
        }
        i -= 1;
    }
    if ref4.is_some() && ((r4 < 4 && j < 4) || (r4 == 4 && j < 5)) && i >= 0 {
        unseen = j as i64;
        a[j] = ptr[i as usize] as i64;
        j += 1;
    }
    let n_alleles = j;
    let has_alt = !(n_alleles == 2 && unseen != -1);

    let x = n_alleles * (n_alleles + 1) / 2;
    let mut g = Vec::with_capacity(x);
    for i in 0..n_alleles {
        for jj in 0..=i {
            g.push((a[jj] * 5 + a[i]) as usize);
        }
    }
    let mut pl = Vec::with_capacity(x * calls.len());
    for c in calls {
        let mut min = f32::MAX;
        for &gi in &g {
            let v = c.p.get(gi).copied().unwrap_or(0.0);
            if min > v {
                min = v;
            }
        }
        for &gi in &g {
            let v = c.p.get(gi).copied().unwrap_or(0.0);
            let y = ((v - min) as f64 + 0.499) as i32;
            pl.push(y.min(255));
        }
    }
    let mut anno = [0f64; 16];
    let mut ori_depth = 0;
    let mut mq0 = 0;
    for c in calls {
        ori_depth += c.ori_depth;
        mq0 += c.mq0;
        for k in 0..16 {
            anno[k] += c.anno[k];
        }
    }
    Combined {
        a,
        n_alleles,
        unseen,
        qsum: out_qsum,
        pl,
        anno,
        ori_depth,
        mq0,
        has_alt,
    }
}

/// `kf_erfc` (AS66).
pub fn kf_erfc(x: f64) -> f64 {
    const P0: f64 = 220.2068679123761;
    const P1: f64 = 221.2135961699311;
    const P2: f64 = 112.0792914978709;
    const P3: f64 = 33.912866078383;
    const P4: f64 = 6.37396220353165;
    const P5: f64 = 0.7003830644436881;
    const P6: f64 = 0.03526249659989109;
    const Q0: f64 = 440.4137358247522;
    const Q1: f64 = 793.8265125199484;
    const Q2: f64 = 637.3336333788311;
    const Q3: f64 = 296.5642487796737;
    const Q4: f64 = 86.78073220294608;
    const Q5: f64 = 16.06417757920695;
    const Q6: f64 = 1.755667163182642;
    const Q7: f64 = 0.08838834764831844;
    let z = x.abs() * std::f64::consts::SQRT_2;
    if z > 37.0 {
        return if x > 0.0 { 0.0 } else { 2.0 };
    }
    let expntl = (z * z * -0.5).exp();
    let p = if z < 10.0 / std::f64::consts::SQRT_2 {
        expntl * ((((((P6 * z + P5) * z + P4) * z + P3) * z + P2) * z + P1) * z + P0)
            / (((((((Q7 * z + Q6) * z + Q5) * z + Q4) * z + Q3) * z + Q2) * z + Q1) * z + Q0)
    } else {
        expntl / 2.506628274631001 / (z + 1.0 / (z + 2.0 / (z + 3.0 / (z + 4.0 / (z + 0.65)))))
    };
    if x > 0.0 { 2.0 * p } else { 2.0 * (1.0 - p) }
}

/// `calc_vdb`: Variant Distance Bias from the alt-base position histogram;
/// `None` below 2x.
pub fn calc_vdb(pos: &[i32; NPOS]) -> Option<f64> {
    const PARAM: [[f32; 3]; 15] = [
        [3.0, 0.079, 18.0],
        [4.0, 0.09, 19.8],
        [5.0, 0.1, 20.5],
        [6.0, 0.11, 21.5],
        [7.0, 0.125, 21.6],
        [8.0, 0.135, 22.0],
        [9.0, 0.14, 22.2],
        [10.0, 0.153, 22.3],
        [15.0, 0.19, 22.8],
        [20.0, 0.22, 23.2],
        [30.0, 0.26, 23.4],
        [40.0, 0.29, 23.5],
        [50.0, 0.35, 23.65],
        [100.0, 0.5, 23.7],
        [200.0, 0.7, 23.7],
    ];
    const READLEN: i64 = 100;
    let mut dp = 0i32;
    let mut mean_pos = 0f32;
    for (i, &c) in pos.iter().enumerate() {
        if c == 0 {
            continue;
        }
        dp += c;
        mean_pos += (c * i as i32) as f32;
    }
    if dp < 2 {
        return None;
    }
    mean_pos /= dp as f32;
    let mut mean_diff = 0f32;
    for (i, &c) in pos.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let diff = (i as f32 - mean_pos).abs() as f64;
        mean_diff = (mean_diff as f64 + c as f64 * diff) as f32;
    }
    mean_diff /= dp as f32;
    let ipos = mean_diff as i64;
    if dp == 2 {
        let v = (2 * READLEN - 2 * (ipos + 1) - 1) * (ipos + 1) / (READLEN - 1);
        return Some(v as f64 / (READLEN as f64 * 0.5));
    }
    let n = PARAM.len();
    let mut i = n;
    if dp < 200 {
        i = 0;
        while i < n {
            if PARAM[i][0] >= dp as f32 {
                break;
            }
            i += 1;
        }
    }
    let (pscale, pshift): (f32, f32) = if i == n {
        (PARAM[n - 1][1], PARAM[n - 1][2])
    } else if i > 0 && PARAM[i][0] != dp as f32 {
        (
            ((PARAM[i - 1][1] + PARAM[i][1]) as f64 * 0.5) as f32,
            ((PARAM[i - 1][2] + PARAM[i][2]) as f64 * 0.5) as f32,
        )
    } else {
        (PARAM[i][1], PARAM[i][2])
    };
    let arg = -((mean_diff - pshift) * pscale);
    Some(0.5 * kf_erfc(arg as f64))
}

/// `calc_mwu_biasZ` with `do_Z`: tie-corrected Mann-Whitney U as a z-score;
/// `None` when one side is empty.
pub fn mwu_z(a: &[i32], b: &[i32]) -> Option<f64> {
    let n = a.len().min(b.len());
    let b_empty = b.iter().all(|&x| x == 0);
    let (mut e, mut l, mut na, mut nb, mut t) = (0i64, 0i64, 0i64, 0i64, 0i64);
    if b_empty {
        for i in (0..n).rev() {
            let ai = a[i] as i64;
            na += ai;
            t += (ai * ai - 1) * ai;
        }
    } else {
        for i in (0..n).rev() {
            let (ai, bi) = (a[i] as i64, b[i] as i64);
            e += ai * bi;
            l += ai * nb;
            na += ai;
            nb += bi;
            let p = ai + bi;
            t += (p * p - 1) * p;
        }
    }
    if na == 0 || nb == 0 {
        return None;
    }
    let u = l as f64 + e as f64 * 0.5;
    let m = (na * nb) as f64 / 2.0;
    let var2 = (na * nb) as f64 / 12.0
        * ((na + nb + 1) as f64 - t as f64 / ((na + nb) * (na + nb - 1)) as f64);
    if var2 <= 0.0 {
        return Some(0.0);
    }
    Some((u - m) / var2.sqrt())
}

fn logsumexp2(a: f64, b: f64) -> f64 {
    if a > b {
        (1.0 + (b - a).exp()).ln() + a
    } else {
        (1.0 + (a - b).exp()).ln() + b
    }
}

/// `calc_SegBias`.
fn calc_seg_bias(calls: &[SampleCall], anno: &[f64; 16]) -> Option<f64> {
    let n = calls.len();
    let nr = (anno[2] + anno[3]) as i64;
    if nr == 0 {
        return None;
    }
    let avg_dp = ((anno[0] + anno[1] + nr as f64) / n as f64) as i64;
    let mut m = (nr as f64 / avg_dp as f64 + 0.5).floor();
    if m > n as f64 {
        m = n as f64;
    } else if m == 0.0 {
        m = 1.0;
    }
    let f = m / 2.0 / n as f64;
    let p = nr as f64 / n as f64;
    let q = nr as f64 / m;
    let log2 = 2f64.ln();
    let mut sum = 0.0;
    for c in calls {
        let oi = (c.anno[2] + c.anno[3]) as i64;
        let tmp = if oi != 0 {
            let mut t = logsumexp2((2.0 * (1.0 - f)).ln(), f.ln() + oi as f64 * log2 - q);
            t += f.ln() + oi as f64 * (q / p).ln() - q + p;
            t
        } else {
            (2.0 * f * (1.0 - f) * (-q).exp() + f * f * (-2.0 * q).exp() + (1.0 - f) * (1.0 - f))
                .ln()
                + p
        };
        sum += tmp;
    }
    Some(sum)
}

fn lbinom(n: i64, k: i64) -> f64 {
    if k == 0 || n == k {
        return 0.0;
    }
    lgamma((n + 1) as f64) - lgamma((k + 1) as f64) - lgamma((n - k + 1) as f64)
}

fn hypergeo(n11: i64, n1_: i64, n_1: i64, n: i64) -> f64 {
    (lbinom(n1_, n11) + lbinom(n - n1_, n_1 - n11) - lbinom(n, n_1)).exp()
}

struct HgAcc {
    n11: i64,
    n1_: i64,
    n_1: i64,
    n: i64,
    p: f64,
}

fn hypergeo_acc(n11: i64, n1_: i64, n_1: i64, n: i64, aux: &mut HgAcc) -> f64 {
    if n1_ != 0 || n_1 != 0 || n != 0 {
        aux.n11 = n11;
        aux.n1_ = n1_;
        aux.n_1 = n_1;
        aux.n = n;
    } else {
        if n11 % 11 != 0 && n11 + aux.n - aux.n1_ - aux.n_1 != 0 {
            if n11 == aux.n11 + 1 {
                aux.p *= (aux.n1_ - aux.n11) as f64 / n11 as f64 * (aux.n_1 - aux.n11) as f64
                    / (n11 + aux.n - aux.n1_ - aux.n_1) as f64;
                aux.n11 = n11;
                return aux.p;
            }
            if n11 == aux.n11 - 1 {
                aux.p *= aux.n11 as f64 / (aux.n1_ - n11) as f64
                    * (aux.n11 + aux.n - aux.n1_ - aux.n_1) as f64
                    / (aux.n_1 - n11) as f64;
                aux.n11 = n11;
                return aux.p;
            }
        }
        aux.n11 = n11;
    }
    aux.p = hypergeo(aux.n11, aux.n1_, aux.n_1, aux.n);
    aux.p
}

/// `kt_fisher_exact`: `(left, right, two)` tail p-values.
pub fn kt_fisher_exact(n11: i64, n12: i64, n21: i64, n22: i64) -> (f64, f64, f64) {
    let n1_ = n11 + n12;
    let n_1 = n11 + n21;
    let n = n11 + n12 + n21 + n22;
    let max = n_1.min(n1_);
    let min = (n1_ + n_1 - n).max(0);
    if min == max {
        return (1.0, 1.0, 1.0);
    }
    let mut aux = HgAcc {
        n11: 0,
        n1_: 0,
        n_1: 0,
        n: 0,
        p: 0.0,
    };
    let q = hypergeo_acc(n11, n1_, n_1, n, &mut aux);
    if q == 0.0 {
        return if n11 * (n + 2) < (n_1 + 1) * (n1_ + 1) {
            (0.0, 1.0, 0.0)
        } else {
            (1.0, 0.0, 0.0)
        };
    }
    let mut p = hypergeo_acc(min, 0, 0, 0, &mut aux);
    let mut left = 0.0;
    let mut i = min + 1;
    while p < 0.99999999 * q && i <= max {
        left += p;
        p = hypergeo_acc(i, 0, 0, 0, &mut aux);
        i += 1;
    }
    i -= 1;
    if p < 1.00000001 * q {
        left += p;
    } else {
        i -= 1;
    }
    p = hypergeo_acc(max, 0, 0, 0, &mut aux);
    let mut right = 0.0;
    let mut j = max - 1;
    while p < 0.99999999 * q && j >= 0 {
        right += p;
        p = hypergeo_acc(j, 0, 0, 0, &mut aux);
        j -= 1;
    }
    j += 1;
    if p < 1.00000001 * q {
        right += p;
    } else {
        j += 1;
    }
    let mut two = left + right;
    if two > 1.0 {
        two = 1.0;
    }
    if (i - n11).abs() < (j - n11).abs() {
        right = 1.0 - left + q;
    } else {
        left = 1.0 - right + q;
    }
    (left, right, two)
}

/// FORMAT/SP: phred-scaled Fisher strand-bias p-value of one sample's DP4.
fn strand_bias_sp(dp4: &[i32; 4]) -> i32 {
    let (fr, rr, fa, ra) = (dp4[0], dp4[1], dp4[2], dp4[3]);
    if fr + rr < 2 || fa + ra < 2 || fr + fa < 2 || rr + ra < 2 {
        return 0;
    }
    let (_, _, two) = kt_fisher_exact(fr as i64, rr as i64, fa as i64, ra as i64);
    if two <= 0.0 {
        return 255;
    }
    let x = (-4.343 * two.ln() + 0.499) as i32;
    x.min(255)
}

fn fmt_f32(v: f64) -> String {
    fmt_g(v as f32 as f64)
}

fn push_list<T: std::fmt::Display>(out: &mut String, vals: impl Iterator<Item = T>) {
    let mut first = true;
    for v in vals {
        if !first {
            out.push(',');
        }
        first = false;
        out.push_str(&v.to_string());
    }
}

/// Render the record once every sample has been folded in.
#[allow(clippy::too_many_arguments)]
fn write_record(
    out: &mut Vec<u8>,
    chr: &str,
    pos1: u64,
    ref_str: &str,
    alts: &[String],
    calls: &[SampleCall],
    comb: &Combined,
    bias: &Bias,
    ctx: &LayoutCtx,
    indel: Option<(i32, f32)>,
) {
    let tags = ctx.tags;
    let n = comb.n_alleles;
    let nsmpl = calls.len();
    let mut info = String::new();
    let sep = |s: &mut String| {
        if !s.is_empty() {
            s.push(';');
        }
    };
    if let Some((idv, imf)) = indel {
        info.push_str("INDEL");
        if tags.info_idv {
            sep(&mut info);
            info.push_str(&format!("IDV={idv}"));
        }
        if tags.info_imf {
            sep(&mut info);
            info.push_str(&format!("IMF={}", fmt_f32(imf as f64)));
        }
    }
    sep(&mut info);
    info.push_str(&format!("DP={}", comb.ori_depth));
    // Per-sample allele counts in output allele order.
    let adf: Vec<Vec<i32>> = calls
        .iter()
        .map(|c| (0..n).map(|j| c.adf[comb.a[j] as usize]).collect())
        .collect();
    let adr: Vec<Vec<i32>> = calls
        .iter()
        .map(|c| (0..n).map(|j| c.adr[comb.a[j] as usize]).collect())
        .collect();
    let adf_tot: Vec<i32> = (0..n).map(|j| adf.iter().map(|v| v[j]).sum()).collect();
    let adr_tot: Vec<i32> = (0..n).map(|j| adr.iter().map(|v| v[j]).sum()).collect();
    let ad_tot: Vec<i32> = (0..n).map(|j| adf_tot[j] + adr_tot[j]).collect();
    if tags.info_adf {
        sep(&mut info);
        info.push_str("ADF=");
        push_list(&mut info, adf_tot.iter());
    }
    if tags.info_adr {
        sep(&mut info);
        info.push_str("ADR=");
        push_list(&mut info, adr_tot.iter());
    }
    if tags.info_ad {
        sep(&mut info);
        info.push_str("AD=");
        push_list(&mut info, ad_tot.iter());
    }
    if tags.info_dpr {
        sep(&mut info);
        info.push_str("DPR=");
        push_list(&mut info, ad_tot.iter());
    }
    if tags.info_scr {
        let scr: i32 = calls.iter().map(|c| c.scr).sum();
        sep(&mut info);
        info.push_str(&format!("SCR={scr}"));
    }
    sep(&mut info);
    info.push_str("I16=");
    push_list(&mut info, comb.anno.iter().map(|&v| fmt_f32(v)));
    sep(&mut info);
    info.push_str("QS=");
    push_list(&mut info, (0..n).map(|j| fmt_f32(comb.qsum[j] as f64)));
    if comb.has_alt {
        if tags.info_vdb
            && let Some(v) = calc_vdb(&bias.alt_pos)
        {
            sep(&mut info);
            info.push_str(&format!("VDB={}", fmt_f32(v)));
        }
        if tags.info_sgb
            && let Some(v) = calc_seg_bias(calls, &comb.anno)
        {
            sep(&mut info);
            info.push_str(&format!("SGB={}", fmt_f32(v)));
        }
        if tags.info_nm && (bias.nnm[0] > 0 || bias.nnm[1] > 0) {
            let nm: Vec<String> = (0..2)
                .map(|i| {
                    fmt_g(if bias.nnm[i] > 0 {
                        (bias.nm[i] / bias.nnm[i] as f32) as f64
                    } else {
                        0.0
                    })
                })
                .collect();
            sep(&mut info);
            info.push_str(&format!("NM={}", nm.join(",")));
        }
        let z_tags: [(bool, &str, &[i32], &[i32]); 6] = [
            (tags.info_rpbz, "RPBZ", &bias.ref_pos, &bias.alt_pos),
            (tags.info_mqbz, "MQBZ", &bias.ref_mq, &bias.alt_mq),
            (tags.info_mqsbz, "MQSBZ", &bias.fwd_mqs, &bias.rev_mqs),
            (tags.info_bqbz, "BQBZ", &bias.ref_bq, &bias.alt_bq),
            (tags.info_nmbz, "NMBZ", &bias.ref_nm, &bias.alt_nm),
            (tags.info_scbz, "SCBZ", &bias.ref_scl, &bias.alt_scl),
        ];
        for (on, name, a, b) in z_tags {
            if !on {
                continue;
            }
            if let Some(z) = mwu_z(a, b) {
                sep(&mut info);
                info.push_str(&format!("{name}={}", fmt_f32(z)));
            }
        }
        if tags.info_fs {
            let (_, _, two) = kt_fisher_exact(
                comb.anno[0] as i64,
                comb.anno[1] as i64,
                comb.anno[2] as i64,
                comb.anno[3] as i64,
            );
            sep(&mut info);
            info.push_str(&format!("FS={}", fmt_f32(two)));
        }
    }
    if tags.info_mq0f {
        let v = if comb.ori_depth > 0 {
            comb.mq0 as f32 / comb.ori_depth as f32
        } else {
            0.0
        };
        sep(&mut info);
        info.push_str(&format!("MQ0F={}", fmt_g(v as f64)));
    }

    let mut fmt_keys = vec!["PL"];
    if tags.fmt_dp {
        fmt_keys.push("DP");
    }
    if tags.fmt_dv {
        fmt_keys.push("DV");
    }
    if tags.fmt_sp {
        fmt_keys.push("SP");
    }
    if tags.fmt_dp4 {
        fmt_keys.push("DP4");
    }
    if tags.fmt_adf {
        fmt_keys.push("ADF");
    }
    if tags.fmt_adr {
        fmt_keys.push("ADR");
    }
    if tags.fmt_ad {
        fmt_keys.push("AD");
    }
    if tags.fmt_dpr {
        fmt_keys.push("DPR");
    }
    if tags.fmt_scr {
        fmt_keys.push("SCR");
    }
    if tags.fmt_qs {
        fmt_keys.push("QS");
    }
    if tags.fmt_nmbz && comb.has_alt {
        fmt_keys.push("NMBZ");
    }

    let _ = write!(
        out,
        "{chr}\t{pos1}\t.\t{ref_str}\t{}\t0\t.\t{info}\t{}",
        alts.join(","),
        fmt_keys.join(":")
    );
    let x = n * (n + 1) / 2;
    for (si, c) in calls.iter().enumerate() {
        let mut col = String::new();
        push_list(&mut col, comb.pl[si * x..(si + 1) * x].iter());
        let dp4 = [
            c.anno[0] as i32,
            c.anno[1] as i32,
            c.anno[2] as i32,
            c.anno[3] as i32,
        ];
        if tags.fmt_dp {
            col.push_str(&format!(":{}", dp4.iter().sum::<i32>()));
        }
        if tags.fmt_dv {
            col.push_str(&format!(":{}", dp4[2] + dp4[3]));
        }
        if tags.fmt_sp {
            col.push_str(&format!(":{}", strand_bias_sp(&dp4)));
        }
        if tags.fmt_dp4 {
            col.push(':');
            push_list(&mut col, dp4.iter());
        }
        if tags.fmt_adf {
            col.push(':');
            push_list(&mut col, adf[si].iter());
        }
        if tags.fmt_adr {
            col.push(':');
            push_list(&mut col, adr[si].iter());
        }
        if tags.fmt_ad || tags.fmt_dpr {
            let ad: Vec<i32> = (0..n).map(|j| adf[si][j] + adr[si][j]).collect();
            if tags.fmt_ad {
                col.push(':');
                push_list(&mut col, ad.iter());
            }
            if tags.fmt_dpr {
                col.push(':');
                push_list(&mut col, ad.iter());
            }
        }
        if tags.fmt_scr {
            col.push_str(&format!(":{}", c.scr));
        }
        if tags.fmt_qs {
            col.push(':');
            push_list(
                &mut col,
                (0..n).map(|j| c.qs[comb.a[j] as usize].min(i32::MAX as i64)),
            );
        }
        if tags.fmt_nmbz && comb.has_alt {
            col.push_str(&format!(
                ":{}",
                fmt_f32(mwu_z(&c.ref_nm, &c.alt_nm).unwrap_or(0.0))
            ));
        }
        let _ = write!(out, "\t{col}");
    }
    let _ = writeln!(out);
    let _ = nsmpl;
}

/// One SNP record: every sample in `kept`, alleles ordered by quality sum,
/// the unseen `<*>` appended.
#[allow(clippy::too_many_arguments)]
pub fn render_snp(
    site: &PileupSite,
    reads: &[LiveRead],
    kept: &[usize],
    ref_char: u8,
    chr: &str,
    pos1: u64,
    ctx: &LayoutCtx,
    out: &mut Vec<u8>,
) {
    let ref4 = ref4_of(ref_char);
    let want_nm = ctx.tags.info_nm || ctx.tags.info_nmbz || ctx.tags.fmt_nmbz;
    let mut bias = Bias::new(want_nm);
    let mut calls: Vec<SampleCall> = Vec::with_capacity(kept.len());
    let mut bases: Vec<u16> = Vec::new();
    for &si in kept {
        let s = &site.per_sample[si];
        let mut call = new_call();
        call.ori_depth = s.raw_depth as i32;
        call.scr = softclip_reads(reads, si);
        bases.clear();
        for &(ri, qi) in &s.obs_src {
            let lr = &reads[ri as usize];
            let qi = qi as usize;
            let b = ref4_of(lr.seq()[qi]);
            let q = adjusted_baseq(lr.qual(), qi, ctx.delta_bq, ctx.max_bq) as i64;
            let is_diff = !(ref4 < 4 && b == ref4);
            add_observation(&mut call, &mut bias, &mut bases, lr, qi, b, q, is_diff);
        }
        call.p = ctx
            .em
            .cal(&mut bases, MAX_ALLELES, &mut ctx.rng.lock().unwrap());
        calls.push(call);
    }
    let comb = combine(&calls, Some(ref4));
    let ref_str = (b"ACGTN"[ref4] as char).to_string();
    let alts: Vec<String> = (1..comb.n_alleles)
        .map(|i| {
            if comb.unseen == i as i64 {
                "<*>".to_string()
            } else {
                (b"ACGT"[comb.a[i] as usize] as char).to_string()
            }
        })
        .collect();
    write_record(
        out, chr, pos1, &ref_str, &alts, &calls, &comb, &bias, ctx, None,
    );
}

/// One biallelic indel record for a kira indel candidate: reads carrying
/// the indel are the alternate allele at the indel quality, every other
/// read the reference at quality 30. The realignment-based indel model of
/// bcftools is not reproduced, only its layout.
#[allow(clippy::too_many_arguments)]
pub fn render_indel(
    site: &PileupSite,
    reads: &[LiveRead],
    kept: &[usize],
    chr: &str,
    pos1: u64,
    spec: &IndelSpec,
    ctx: &LayoutCtx,
    out: &mut Vec<u8>,
) {
    let pos0 = (pos1 - 1) as u32;
    let want_nm = ctx.tags.info_nm || ctx.tags.info_nmbz || ctx.tags.fmt_nmbz;
    let mut bias = Bias::new(want_nm);
    let mut calls: Vec<SampleCall> = Vec::with_capacity(kept.len());
    let mut bases: Vec<u16> = Vec::new();
    let mut idv = 0i32;
    let mut imf = 0f32;
    for &si in kept {
        let s = &site.per_sample[si];
        let mut call = new_call();
        call.ori_depth = s.raw_depth as i32;
        call.scr = softclip_reads(reads, si);
        bases.clear();
        let mut support = 0i32;
        for &(ri, qi) in &s.obs_src {
            let lr = &reads[ri as usize];
            let qi = qi as usize;
            let carries = match lr.indel_after(pos0) {
                Some(IndelEvent::Ins(seq)) => spec.is_ins && seq == spec.ins_seq,
                Some(IndelEvent::Del(l)) => !spec.is_ins && l == spec.del_len,
                None => false,
            };
            let (b, q) = if carries {
                (1usize, spec.indel_qual as i64)
            } else {
                (0usize, 30i64)
            };
            if carries {
                support += 1;
            }
            add_observation(&mut call, &mut bias, &mut bases, lr, qi, b, q, carries);
        }
        call.p = ctx
            .em
            .cal(&mut bases, MAX_ALLELES, &mut ctx.rng.lock().unwrap());
        if support > idv {
            idv = support;
        }
        if s.raw_depth > 0 {
            let f = support as f32 / s.raw_depth as f32;
            if f > imf {
                imf = f;
            }
        }
        calls.push(call);
    }
    let comb = combine(&calls, None);
    if comb.n_alleles < 2 {
        return;
    }
    let alts = vec![spec.alt_str.to_string()];
    write_record(
        out,
        chr,
        pos1,
        spec.ref_str,
        &alts,
        &calls,
        &comb,
        &bias,
        ctx,
        Some((idv, imf)),
    );
}

/// Header lines bcftools `mpileup` writes for the requested annotations, in
/// its order; contigs and the reference line are the caller's.
pub fn header_lines(tags: &AnnotateSpec) -> Vec<&'static str> {
    let mut v = vec![
        "##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">",
        "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">",
    ];
    if tags.info_idv {
        v.push("##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of raw reads supporting an indel\">");
    }
    if tags.info_imf {
        v.push("##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of raw reads supporting an indel\">");
    }
    v.push("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
    if tags.info_vdb {
        v.push("##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">");
    }
    if tags.info_rpbz {
        v.push("##INFO=<ID=RPBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)\">");
    }
    if tags.info_mqbz {
        v.push("##INFO=<ID=MQBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)\">");
    }
    if tags.info_bqbz {
        v.push("##INFO=<ID=BQBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)\">");
    }
    if tags.info_mqsbz {
        v.push("##INFO=<ID=MQSBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)\">");
    }
    if tags.info_nm {
        v.push("##INFO=<ID=NM,Number=2,Type=Float,Description=\"Average number of mismatches in ref and alt reads (approximate, experimental, make me localized?)\">");
    }
    if tags.info_nmbz {
        v.push("##INFO=<ID=NMBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Number of Mismatches within supporting reads (closer to 0 is better; approximate, experimental, make me localized?)\">");
    }
    if tags.fmt_nmbz {
        v.push("##FORMAT=<ID=NMBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Number of Mismatches within supporting reads (closer to 0 is better)\">");
    }
    if tags.info_scbz {
        v.push("##INFO=<ID=SCBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)\">");
    }
    if tags.info_fs {
        v.push("##INFO=<ID=FS,Number=1,Type=Float,Description=\"Fisher's exact test P-value to detect strand bias\">");
    }
    if tags.info_sgb {
        v.push("##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf\">");
    }
    if tags.info_mq0f {
        v.push("##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">");
    }
    v.push("##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">");
    v.push("##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">");
    v.push("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
    if tags.fmt_dp {
        v.push(
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">",
        );
    }
    if tags.fmt_dv {
        v.push("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">");
    }
    if tags.fmt_dpr {
        v.push("##FORMAT=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    }
    if tags.info_dpr {
        v.push("##INFO=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    }
    if tags.fmt_dp4 {
        v.push("##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases\">");
    }
    if tags.fmt_sp {
        v.push("##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">");
    }
    if tags.fmt_ad {
        v.push("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths (high-quality bases)\">");
    }
    if tags.fmt_adf {
        v.push("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand (high-quality bases)\">");
    }
    if tags.fmt_adr {
        v.push("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand (high-quality bases)\">");
    }
    if tags.fmt_qs {
        v.push("##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Phred-score allele quality sum used by `call -mG` and `+trio-dnm`\">");
    }
    if tags.info_ad {
        v.push("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths (high-quality bases)\">");
    }
    if tags.info_adf {
        v.push("##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand (high-quality bases)\">");
    }
    if tags.info_scr {
        v.push("##INFO=<ID=SCR,Number=1,Type=Integer,Description=\"Number of soft-clipped reads (at high-quality bases)\">");
    }
    if tags.fmt_scr {
        v.push("##FORMAT=<ID=SCR,Number=1,Type=Integer,Description=\"Per-sample number of soft-clipped reads (at high-quality bases)\">");
    }
    if tags.info_adr {
        v.push("##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand (high-quality bases)\">");
    }
    v
}

#[cfg(test)]
#[path = "../../tests/unit/bam_bcf_layout.rs"]
mod tests;
