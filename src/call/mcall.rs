//! Multi-allelic Bayesian variant caller (bcftools `call -m`, vcfcall/mcall.c).
//!
//! 1. FORMAT/PL → P(D|G) per sample, normalised per sample (`set_pdg`).
//! 2. Allele frequencies come from INFO/QS (per-group FORMAT/AD fractions
//!    with `-G`), optionally shrunk towards `-F AN,AC` panel counts. They are
//!    kept in single precision, as bcftools does, so QUAL agrees to the digit.
//! 3. Enumerate allele subsets (singletons, pairs, triples) per group; each
//!    gets the HWE likelihood at those frequencies plus a θ prior per
//!    non-ref allele. The union of the best subsets is the call; `lk_sum`
//!    is the log-sum over the non-reference hypotheses.
//! 4. QUAL = -4.343·(ref_lk − logsumexp(lk_sum, ref_lk)) for variant sites
//!    (the best group), -4.343·(lk_sum − logsumexp(lk_sum, ref_lk)) for
//!    reference sites.
//! 5. Per-sample GT = argmax of P(D|G)·HWE(qsum) over the kept alleles;
//!    GQ = phred(1 − P(best)) capped at 127; GP the normalised posteriors;
//!    PL is the input PL subset.
//! 6. `-C trio`: the trio genotypes maximise the joint likelihood with a
//!    Mendelian transmission prior and a de-novo rate.

use super::math::*;

#[derive(Clone, Debug)]
pub struct CallerOpts {
    pub theta: f64,
    pub indel_theta: f64,
    pub keep_alts: bool,
    pub variants_only: bool,
    pub min_ac: u32,
    pub ploidy: u8,
    pub per_sample_ploidy: Option<Vec<u8>>,
    /// Chromosomes behind the Watterson factor; bcftools counts every
    /// sample at the maximum ploidy of the ploidy file.
    pub prior_alleles: Option<usize>,
    pub gvcf: Option<GvcfOpts>,
    pub groups: Option<Vec<SampleGroup>>,
    pub constrain: ConstrainMode,
    pub families: Vec<TrioFamily>,
    /// De-novo rates for SNPs, deletions and insertions (`-X`).
    pub novel_rate: [f64; 3],
}

#[derive(Clone, Debug)]
pub struct GvcfOpts {
    pub min_dp: u32,
    pub min_qual: f64,
    pub blocks: Vec<u32>,
}

#[derive(Clone, Debug)]
pub struct SampleGroup {
    pub name: String,
    pub sample_idxs: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConstrainMode {
    None,
    Alleles,
    Trio,
}

#[derive(Clone, Debug)]
pub struct TrioFamily {
    pub father: Option<usize>,
    pub mother: Option<usize>,
    pub child: Option<usize>,
    pub is_son: bool,
}

impl Default for CallerOpts {
    fn default() -> Self {
        Self {
            theta: 1.1e-3,
            indel_theta: 1.5e-4,
            keep_alts: false,
            variants_only: false,
            min_ac: 0,
            ploidy: 2,
            per_sample_ploidy: None,
            prior_alleles: None,
            gvcf: None,
            groups: None,
            constrain: ConstrainMode::None,
            families: Vec::new(),
            novel_rate: [1e-8, 1e-9, 1e-9],
        }
    }
}

pub struct Caller {
    pub opts: CallerOpts,
    pl2p: [f64; PL2P_SIZE],
    log_theta_snv: f64,
    log_theta_indel: f64,
}

/// GQ is stored in a signed byte by bcftools.
const GQ_MAX: f64 = 127.0;
/// The `-4.343` of mcall.c.
const PHRED: f64 = -4.343;

/// One site's input to the caller.
pub struct CallSite {
    pub n_samples: usize,
    /// Alleles in the record, the unseen `<*>` included.
    pub n_alleles: usize,
    /// Diploid PLs, `n_genotypes(n_alleles)` per sample; [`PL_MISSING`] for
    /// `.`, [`PL_VECTOR_END`] past a short vector. Missing values are
    /// filled in place, as `set_pdg` does.
    pub pls: Vec<i32>,
    pub is_indel: bool,
    pub depths: Option<Vec<u32>>,
    /// INFO/QS per allele (bcftools' frequency estimate); PL-derived when absent.
    pub qs: Option<Vec<f64>>,
    /// Index of the unseen allele (`<*>`), never emitted.
    pub unseen: Option<usize>,
    /// Per-sample allele fractions (FORMAT/AD) for `-G` group calling.
    pub sample_af: Option<Vec<Vec<f64>>>,
    /// `-F AN,AC` panel counts: (AN, AC per ALT).
    pub prior_an_ac: Option<(f64, Vec<f64>)>,
}

impl CallSite {
    pub fn new(n_samples: usize, n_alleles: usize, pls: Vec<i32>) -> Self {
        Self { n_samples, n_alleles, pls, is_indel: false, depths: None, qs: None, unseen: None, sample_af: None, prior_an_ac: None }
    }
}

/// Missing PL value in [`CallSite::pls`] and [`CallResult::pls`].
pub const PL_MISSING: i32 = i32::MIN;
/// End of a PL vector shorter than the diploid genotype count.
pub const PL_VECTOR_END: i32 = i32::MIN + 1;

pub enum CallResult {
    Called {
        alleles_kept: Vec<u32>,
        /// A non-reference allele set was the best hypothesis (bcftools'
        /// `is_variant`; the genotypes may still all be 0/0).
        is_variant: bool,
        /// `None` prints as `.`: no read support and no reference call.
        qual: Option<f64>,
        /// Original allele ids; `None` is a missing genotype.
        gts: Vec<Option<(u32, u32)>>,
        gqs: Vec<u32>,
        /// Posterior genotype probabilities over the kept genotypes (GP).
        gps: Vec<Vec<f64>>,
        /// Input PLs restricted to the kept alleles (`n_genotypes(kept)`
        /// values, `kept` for haploids, empty when PL is dropped).
        pls: Vec<Vec<i32>>,
        ac: Vec<u32>,
        an: u32,
    },
    Skip,
}

struct BestAlleles {
    als: u32,
    max_lk: f64,
    ref_lk: f64,
    lk_sum: f64,
}

impl Caller {
    pub fn new(opts: CallerOpts, n_samples: usize) -> Self {
        let n_total = opts.prior_alleles.unwrap_or_else(|| total_alleles(n_samples, &opts));
        let a_m = watterson_factor(n_total);
        let scaled_snv = (opts.theta * a_m).min(0.99);
        let scaled_indel = (opts.indel_theta * a_m).min(0.99);
        let log_theta_snv = if scaled_snv > 0.0 { scaled_snv.ln() } else { f64::NEG_INFINITY };
        let log_theta_indel = if scaled_indel > 0.0 { scaled_indel.ln() } else { f64::NEG_INFINITY };
        Self { opts, pl2p: init_pl2p(), log_theta_snv, log_theta_indel }
    }

    fn ploidy_for(&self, sample_idx: usize) -> u8 {
        if let Some(p) = &self.opts.per_sample_ploidy {
            return p.get(sample_idx).copied().unwrap_or(self.opts.ploidy);
        }
        self.opts.ploidy
    }

    fn log_theta_for(&self, is_indel: bool) -> f64 {
        if is_indel { self.log_theta_indel } else { self.log_theta_snv }
    }

    /// Allele frequency estimate of one sample group, as `mcall()` builds
    /// it, in bcftools' single precision.
    fn group_qsum(&self, site: &CallSite, group: &[usize], pdg: &[f64]) -> Vec<f32> {
        let n_als = site.n_alleles;
        let mut qsum: Vec<f32> = if let Some(af) = &site.sample_af {
            let mut q = vec![0.0f32; n_als];
            for &si in group {
                if let Some(row) = af.get(si) {
                    for (j, v) in row.iter().enumerate().take(n_als) {
                        q[j] += *v as f32;
                    }
                }
            }
            q
        } else if let Some(qs) = &site.qs {
            let mut q: Vec<f32> = qs.iter().map(|v| *v as f32).collect();
            q.resize(n_als, 0.0);
            q
        } else {
            compute_qsum(pdg, group, n_als).into_iter().map(|v| v as f32).collect()
        };
        if let Some((an, ac)) = &site.prior_an_ac {
            if *an > 0.0 && ac.len() == n_als - 1 {
                let n = group.len() as f64;
                let mut ac0 = *an;
                for (i, &c) in ac.iter().enumerate() {
                    ac0 -= c;
                    qsum[i + 1] = ((qsum[i + 1] as f64 + 0.5 * c) / (n + 0.5 * an)) as f32;
                }
                qsum[0] = ((qsum[0] as f64 + 0.5 * ac0.max(0.0)) / (n + 0.5 * an)) as f32;
            }
        }
        let total: f32 = qsum.iter().sum();
        if total > 0.0 {
            for q in qsum.iter_mut() {
                *q /= total;
            }
        }
        qsum
    }

    /// Call a single site given per-sample PLs.
    pub fn call_site(&self, site: &mut CallSite) -> CallResult {
        let n_smpl = site.n_samples;
        let n_als = site.n_alleles;
        let n_gt = n_genotypes(n_als);
        let groups: Vec<Vec<usize>> = match &self.opts.groups {
            Some(g) => g.iter().map(|x| x.sample_idxs.iter().copied().filter(|&i| i < n_smpl).collect()).collect(),
            None => vec![(0..n_smpl).collect()],
        };
        let ploidies: Vec<u8> = (0..n_smpl).map(|si| self.ploidy_for(si)).collect();

        let mut pdg: Vec<f64> = vec![0.0; n_smpl * n_gt];
        for si in 0..n_smpl {
            set_pdg_row(&mut site.pls[si * n_gt..(si + 1) * n_gt], &mut pdg[si * n_gt..(si + 1) * n_gt], &self.pl2p, n_als, site.unseen);
        }
        let log_theta = self.log_theta_for(site.is_indel);

        // Best allele set per group; the union is the call.
        let mut als_new: u32 = 0;
        let mut best: Option<(f64, f64, f64)> = None; // (qual, lk_sum, ref_lk) of the best group
        let mut group_qsums: Vec<Vec<f32>> = Vec::with_capacity(groups.len());
        let mut group_als: Vec<u32> = Vec::with_capacity(groups.len());
        for g in &groups {
            let qsum = self.group_qsum(site, g, &pdg);
            let r = find_best_alleles(&pdg, g, n_als, &qsum, &ploidies, log_theta);
            als_new |= r.als;
            group_als.push(r.als);
            if r.max_lk.is_finite() {
                let qual = PHRED * (r.ref_lk - ln_sum_exp(r.lk_sum, r.ref_lk));
                if best.is_none_or(|b| qual > b.0) {
                    best = Some((qual, r.lk_sum, r.ref_lk));
                }
            }
            group_qsums.push(qsum);
        }
        als_new |= 1;
        if let Some(u) = site.unseen {
            if u > 0 {
                als_new &= !(1u32 << u);
            }
        }
        let is_variant = als_new != 1;
        if self.opts.variants_only && !is_variant {
            return CallResult::Skip;
        }
        if self.opts.keep_alts {
            for i in 1..n_als {
                if Some(i) != site.unseen {
                    als_new |= 1 << i;
                }
            }
        }
        let alleles_kept: Vec<u32> = (0..n_als as u32).filter(|&i| als_new & (1 << i) != 0).collect();
        let nk = alleles_kept.len();

        let mut gts: Vec<Option<(u32, u32)>> = vec![None; n_smpl];
        let mut gqs: Vec<u32> = vec![0; n_smpl];
        let mut gps: Vec<Vec<f64>> = vec![Vec::new(); n_smpl];
        let mut ac = vec![0u32; nk];
        let drop_pl;
        if !is_variant {
            // Reference site: 0/0 wherever there is any likelihood.
            drop_pl = nk == 1;
            for si in 0..n_smpl {
                let row = &pdg[si * n_gt..(si + 1) * n_gt];
                if ploidies[si] == 0 || row.iter().all(|&p| p == 0.0) {
                    continue;
                }
                gts[si] = Some((0, 0));
                ac[0] += ploidies[si].min(2) as u32;
            }
        } else {
            drop_pl = false;
            for (gi, g) in groups.iter().enumerate() {
                // Genotypes range over the group's own best allele set, as
                // `mcall_call_genotypes` does with `grp->als`.
                let mask = group_als[gi] & als_new;
                for &si in g {
                    let row = &pdg[si * n_gt..(si + 1) * n_gt];
                    let (gt, gq, gp) = call_genotype(row, &alleles_kept, mask, &group_qsums[gi], ploidies[si]);
                    gts[si] = gt;
                    gqs[si] = gq;
                    gps[si] = gp;
                    if let Some((a, b)) = gt {
                        ac[pos_of(&alleles_kept, a)] += 1;
                        if ploidies[si] == 2 {
                            ac[pos_of(&alleles_kept, b)] += 1;
                        }
                    }
                }
            }
            if self.opts.constrain == ConstrainMode::Trio {
                for fam in &self.opts.families {
                    self.call_trio(fam, &pdg, n_gt, &alleles_kept, &group_qsums[0], &ploidies, site.is_indel, &mut gts, &mut gqs, &mut gps);
                }
                ac.iter_mut().for_each(|c| *c = 0);
                for (si, gt) in gts.iter().enumerate() {
                    if let Some((a, b)) = gt {
                        ac[pos_of(&alleles_kept, *a)] += 1;
                        if ploidies[si] == 2 {
                            ac[pos_of(&alleles_kept, *b)] += 1;
                        }
                    }
                }
            }
        }
        let n_ac: u32 = ac.iter().skip(1).sum();
        if n_ac == 0 && self.opts.variants_only {
            return CallResult::Skip;
        }
        let an: u32 = ac.iter().sum();

        let qual = if n_ac > 0 {
            best.map(|b| b.0)
        } else {
            match best {
                Some((_, lk_sum, ref_lk)) if lk_sum.is_finite() => Some(PHRED * (lk_sum - ln_sum_exp(lk_sum, ref_lk))),
                _ if ac[0] > 0 => Some(if log_theta.is_finite() { PHRED * log_theta } else { 0.0 }),
                _ => None,
            }
        };

        let pls = if drop_pl { vec![Vec::new(); n_smpl] } else { subset_pls(site, &alleles_kept, &ploidies) };
        CallResult::Called { alleles_kept, is_variant, qual, gts, gqs, gps, pls, ac, an }
    }

    /// Joint trio call: the (father, mother, child) genotype triple maximising
    /// P(D_f|G_f)P(G_f) · P(D_m|G_m)P(G_m) · P(D_c|G_c)P(G_c) · T(G_c | G_f, G_m),
    /// with T = ν + (1 − ν)·P_Mendel and ν the de-novo rate. GQ and GP
    /// come from the joint marginals.
    #[allow(clippy::too_many_arguments)]
    fn call_trio(
        &self,
        fam: &TrioFamily,
        pdg: &[f64],
        n_gt: usize,
        kept: &[u32],
        qsum: &[f32],
        ploidies: &[u8],
        is_indel: bool,
        gts: &mut [Option<(u32, u32)>],
        gqs: &mut [u32],
        gps: &mut [Vec<f64>],
    ) {
        let (Some(fi), Some(mi), Some(ci)) = (fam.father, fam.mother, fam.child) else { return };
        if fi >= ploidies.len() || mi >= ploidies.len() || ci >= ploidies.len() {
            return;
        }
        if ploidies[fi] != 2 || ploidies[mi] != 2 || ploidies[ci] != 2 {
            return;
        }
        let nk = kept.len();
        let n_new = n_genotypes(nk);
        let pairs: Vec<(usize, usize)> = (0..n_new).map(idx_to_pair).collect();
        let novel = if is_indel { self.opts.novel_rate[1] } else { self.opts.novel_rate[0] }.clamp(1e-300, 1.0);
        let logs = |si: usize| -> Option<Vec<f64>> {
            let row = &pdg[si * n_gt..(si + 1) * n_gt];
            if row.iter().all(|&p| p == 0.0) {
                return None;
            }
            Some(genotype_weights(row, kept, u32::MAX, qsum, 2).iter().map(|w| if *w > 0.0 { w.ln() } else { f64::NEG_INFINITY }).collect())
        };
        let (Some(lf), Some(lm), Some(lc)) = (logs(fi), logs(mi), logs(ci)) else { return };
        let mut best = f64::NEG_INFINITY;
        let mut best_triple = (0usize, 0usize, 0usize);
        let mut marg = [vec![f64::NEG_INFINITY; n_new], vec![f64::NEG_INFINITY; n_new], vec![f64::NEG_INFINITY; n_new]];
        for (a, &pf) in lf.iter().enumerate() {
            if pf == f64::NEG_INFINITY {
                continue;
            }
            for (b, &pm) in lm.iter().enumerate() {
                if pm == f64::NEG_INFINITY {
                    continue;
                }
                for (c, &pc) in lc.iter().enumerate() {
                    if pc == f64::NEG_INFINITY {
                        continue;
                    }
                    let mendel = transmissions(pairs[a], pairs[b], pairs[c]) as f64 / 4.0;
                    let joint = pf + pm + pc + (novel + (1.0 - novel) * mendel).ln();
                    if joint > best {
                        best = joint;
                        best_triple = (a, b, c);
                    }
                    marg[0][a] = ln_sum_exp(marg[0][a], joint);
                    marg[1][b] = ln_sum_exp(marg[1][b], joint);
                    marg[2][c] = ln_sum_exp(marg[2][c], joint);
                }
            }
        }
        if best == f64::NEG_INFINITY {
            return;
        }
        for (si, m, gi) in [(fi, &marg[0], best_triple.0), (mi, &marg[1], best_triple.1), (ci, &marg[2], best_triple.2)] {
            let (i, j) = pairs[gi];
            gts[si] = Some((kept[i], kept[j]));
            let weights: Vec<f64> = {
                let max = m.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                m.iter().map(|v| if v.is_finite() { (v - max).exp() } else { 0.0 }).collect()
            };
            let (gq, gp) = gq_gp(&weights);
            gqs[si] = gq;
            gps[si] = gp;
        }
    }
}

/// `mcall_find_best_alleles`: the most likely allele subset of one group.
fn find_best_alleles(pdg: &[f64], group: &[usize], n_als: usize, qsum: &[f32], ploidies: &[u8], log_theta: f64) -> BestAlleles {
    let n_gt = n_genotypes(n_als);
    let mut max_lk = f64::NEG_INFINITY;
    let mut max_als = 0u32;
    let mut ref_lk = f64::NEG_INFINITY;
    let mut lk_sum = f64::NEG_INFINITY;

    for ia in 0..n_als {
        let iaa = gt_index(ia, ia);
        let mut lk_tot = 0.0;
        let mut set = false;
        for &si in group {
            let v = pdg[si * n_gt + iaa];
            if v > 0.0 {
                lk_tot += v.ln();
                set = true;
            }
        }
        if ia == 0 {
            ref_lk = lk_tot;
        } else {
            lk_tot += log_theta;
        }
        if set && lk_tot > max_lk {
            max_lk = lk_tot;
            max_als = 1 << ia;
        }
        if ia > 0 && set {
            lk_sum = ln_sum_exp(lk_sum, lk_tot);
        }
    }

    if n_als > 1 {
        for ia in 0..n_als {
            if qsum[ia] == 0.0 {
                continue;
            }
            for ib in 0..ia {
                if qsum[ib] == 0.0 {
                    continue;
                }
                // Single precision, as in mcall.c.
                let fa = (qsum[ia] / (qsum[ia] + qsum[ib])) as f64;
                let fb = (qsum[ib] / (qsum[ia] + qsum[ib])) as f64;
                let (fa2, fb2, fab) = (fa * fa, fb * fb, 2.0 * fa * fb);
                let (iaa, ibb, iab) = (gt_index(ia, ia), gt_index(ib, ib), gt_index(ia, ib));
                let mut lk_tot = 0.0;
                let mut set = false;
                for &si in group {
                    let row = &pdg[si * n_gt..];
                    let val = match ploidies[si] {
                        2 => fa2 * row[iaa] + fb2 * row[ibb] + fab * row[iab],
                        1 => fa * row[iaa] + fb * row[ibb],
                        _ => 0.0,
                    };
                    if val > 0.0 {
                        lk_tot += val.ln();
                        set = true;
                    }
                }
                if ia != 0 {
                    lk_tot += log_theta;
                }
                if ib != 0 {
                    lk_tot += log_theta;
                }
                if set && lk_tot > max_lk {
                    max_lk = lk_tot;
                    max_als = (1 << ia) | (1 << ib);
                }
                if set {
                    lk_sum = ln_sum_exp(lk_sum, lk_tot);
                }
            }
        }
    }

    if n_als > 2 {
        for ia in 0..n_als {
            if qsum[ia] == 0.0 {
                continue;
            }
            for ib in 0..ia {
                if qsum[ib] == 0.0 {
                    continue;
                }
                for ic in 0..ib {
                    if qsum[ic] == 0.0 {
                        continue;
                    }
                    let total = qsum[ia] + qsum[ib] + qsum[ic];
                    let (fa, fb, fc) = ((qsum[ia] / total) as f64, (qsum[ib] / total) as f64, (qsum[ic] / total) as f64);
                    let (fa2, fb2, fc2) = (fa * fa, fb * fb, fc * fc);
                    let (fab, fac, fbc) = (2.0 * fa * fb, 2.0 * fa * fc, 2.0 * fb * fc);
                    let (iaa, ibb, icc) = (gt_index(ia, ia), gt_index(ib, ib), gt_index(ic, ic));
                    let (iab, iac, ibc) = (gt_index(ia, ib), gt_index(ia, ic), gt_index(ib, ic));
                    let mut lk_tot = 0.0;
                    let mut set = false;
                    for &si in group {
                        let row = &pdg[si * n_gt..];
                        let val = match ploidies[si] {
                            2 => fa2 * row[iaa] + fb2 * row[ibb] + fc2 * row[icc] + fab * row[iab] + fac * row[iac] + fbc * row[ibc],
                            1 => fa * row[iaa] + fb * row[ibb] + fc * row[icc],
                            _ => 0.0,
                        };
                        if val > 0.0 {
                            lk_tot += val.ln();
                            set = true;
                        }
                    }
                    for a in [ia, ib, ic] {
                        if a != 0 {
                            lk_tot += log_theta;
                        }
                    }
                    if set && lk_tot > max_lk {
                        max_lk = lk_tot;
                        max_als = (1 << ia) | (1 << ib) | (1 << ic);
                    }
                    if set {
                        lk_sum = ln_sum_exp(lk_sum, lk_tot);
                    }
                }
            }
        }
    }
    BestAlleles { als: max_als, max_lk, ref_lk, lk_sum }
}

/// P(D|G)·HWE(qsum) over the kept genotypes (VCF order), as
/// `mcall_call_genotypes` weighs them; genotypes with an allele outside
/// `mask` (bits over original allele ids) weigh 0. Haploids get one
/// weight per allele.
fn genotype_weights(row: &[f64], kept: &[u32], mask: u32, qsum: &[f32], ploidy: u8) -> Vec<f64> {
    let nk = kept.len();
    let allowed = |a: u32| mask & (1 << a) != 0;
    if ploidy == 1 {
        return kept.iter().map(|&a| if allowed(a) { row[gt_index(a as usize, a as usize)] * qsum[a as usize] as f64 } else { 0.0 }).collect();
    }
    let mut w = vec![0.0; n_genotypes(nk)];
    for j in 0..nk {
        for i in 0..=j {
            if !allowed(kept[i]) || !allowed(kept[j]) {
                continue;
            }
            let (a, b) = (kept[i] as usize, kept[j] as usize);
            let (qa, qb) = (qsum[a] as f64, qsum[b] as f64);
            w[gt_index(i, j)] = if i == j { row[gt_index(a, a)] * qa * qa } else { 2.0 * row[gt_index(a, b)] * qa * qb };
        }
    }
    w
}

/// One sample's GT (original allele ids), GQ and GP. Ties keep the earlier
/// genotype in bcftools' iteration order: homozygotes by allele, then
/// heterozygotes; 0/0 when every weight is zero.
fn call_genotype(row: &[f64], kept: &[u32], mask: u32, qsum: &[f32], ploidy: u8) -> (Option<(u32, u32)>, u32, Vec<f64>) {
    if ploidy == 0 || row.iter().all(|&p| p == 0.0) {
        return (None, 0, Vec::new());
    }
    let w = genotype_weights(row, kept, mask, qsum, ploidy);
    let nk = kept.len();
    let mut best_lk = 0.0;
    let mut best = (kept[0], kept[0]);
    if ploidy == 1 {
        for (i, &v) in w.iter().enumerate() {
            if best_lk < v {
                best_lk = v;
                best = (kept[i], kept[i]);
            }
        }
    } else {
        for i in 0..nk {
            let v = w[gt_index(i, i)];
            if best_lk < v {
                best_lk = v;
                best = (kept[i], kept[i]);
            }
        }
        for j in 0..nk {
            for i in 0..j {
                let v = w[gt_index(i, j)];
                if best_lk < v {
                    best_lk = v;
                    best = (kept[i], kept[j]);
                }
            }
        }
    }
    let (gq, gp) = gq_gp(&w);
    (Some(best), gq, gp)
}

/// GQ = phred(1 − w_best/Σw) truncated and capped at 127; GP = w/Σw.
fn gq_gp(w: &[f64]) -> (u32, Vec<f64>) {
    let sum: f64 = w.iter().sum();
    if sum <= 0.0 {
        return (0, vec![0.0; w.len()]);
    }
    let max = w.iter().copied().fold(0.0, f64::max);
    let p = max / sum;
    let gq = if p >= 1.0 { GQ_MAX } else { (-4.34294 * (1.0 - p).ln()).min(GQ_MAX) };
    (gq.max(0.0) as u32, w.iter().map(|v| v / sum).collect())
}

fn pos_of(kept: &[u32], allele: u32) -> usize {
    kept.iter().position(|&k| k == allele).unwrap_or(0)
}

/// Number of the four parental transmissions (one allele from each parent)
/// that produce the child's genotype.
fn transmissions(f: (usize, usize), m: (usize, usize), c: (usize, usize)) -> u32 {
    let mut n = 0;
    for &a in &[f.0, f.1] {
        for &b in &[m.0, m.1] {
            if (a.min(b), a.max(b)) == (c.0.min(c.1), c.0.max(c.1)) {
                n += 1;
            }
        }
    }
    n
}

fn total_alleles(n_samples: usize, opts: &CallerOpts) -> usize {
    if let Some(p) = &opts.per_sample_ploidy {
        p.iter().take(n_samples).map(|x| *x as usize).sum()
    } else {
        n_samples * opts.ploidy as usize
    }
}

#[inline]
fn gt_idx(a: usize, b: usize) -> usize {
    gt_index(a.min(b), a.max(b))
}

/// `set_pdg`: PL → probabilities normalised per sample. A leading missing
/// value or a short vector makes the sample missing; inner missing values
/// are filled from the unseen allele's PLs (255 without one); a flat
/// vector (every PL 0) is missing too. The fills stay in `pls`.
fn set_pdg_row(pls: &mut [i32], pdg: &mut [f64], pl2p: &[f64; PL2P_SIZE], n_als: usize, unseen: Option<usize>) {
    let n_gt = pdg.len();
    let mut sum = 0.0;
    let mut j = 0;
    while j < n_gt {
        let v = pls[j];
        if v == PL_VECTOR_END {
            j = 0;
            break;
        }
        if v == PL_MISSING {
            break;
        }
        pdg[j] = pl_to_prob(v, pl2p);
        sum += pdg[j];
        j += 1;
    }
    if j == 0 {
        pdg.fill(0.0);
        return;
    }
    if j < n_gt {
        sum = 0.0;
        match unseen {
            None => {
                for k in 0..n_gt {
                    if pls[k] == PL_MISSING {
                        pls[k] = 255;
                    }
                    pdg[k] = pl_to_prob(pls[k], pl2p);
                    sum += pdg[k];
                }
            }
            Some(u) => {
                let mut k = 0;
                for ia in 0..n_als {
                    for ib in 0..=ia {
                        if pls[k] == PL_MISSING {
                            let mut src = gt_idx(ia, u);
                            if pls[src] == PL_MISSING {
                                src = gt_idx(ib, u);
                            }
                            if pls[src] == PL_MISSING {
                                src = gt_idx(u, u);
                            }
                            pls[k] = if pls[src] == PL_MISSING { 255 } else { pls[src] };
                        }
                        pdg[k] = pl_to_prob(pls[k], pl2p);
                        sum += pdg[k];
                        k += 1;
                    }
                }
            }
        }
    }
    if n_gt > 1 && sum == n_gt as f64 {
        pdg.fill(0.0);
    } else {
        for p in pdg.iter_mut() {
            *p /= sum;
        }
    }
}

/// Frequency estimate from the PLs when INFO/QS is absent (kira fallback).
fn compute_qsum(pdg: &[f64], group: &[usize], n_als: usize) -> Vec<f64> {
    let n_gt = n_genotypes(n_als);
    let mut qsum = vec![0.0f64; n_als];
    for &s in group {
        let row = &pdg[s * n_gt..(s + 1) * n_gt];
        for a in 0..n_als {
            qsum[a] += row[gt_index(a, a)];
            for b in 0..n_als {
                if b != a {
                    qsum[a] += 0.5 * row[gt_idx(a, b)];
                }
            }
        }
    }
    qsum
}

fn idx_to_pair(idx: usize) -> (usize, usize) {
    let mut k = idx;
    let mut j = 0;
    while k > j {
        k -= j + 1;
        j += 1;
    }
    (k, j)
}

/// Input PLs restricted to the kept alleles (`mcall_trim_and_update_PLs`):
/// bcftools keeps the original values. Haploids get one value per allele.
fn subset_pls(site: &CallSite, kept: &[u32], ploidies: &[u8]) -> Vec<Vec<i32>> {
    let n_gt = n_genotypes(site.n_alleles);
    let nk = kept.len();
    (0..site.n_samples)
        .map(|si| {
            let row = &site.pls[si * n_gt..(si + 1) * n_gt];
            let missing = |v: i32| if v == PL_MISSING || v == PL_VECTOR_END { PL_MISSING } else { v };
            match ploidies.get(si).copied().unwrap_or(2) {
                1 => kept.iter().map(|&a| missing(row[gt_index(a as usize, a as usize)])).collect(),
                0 => vec![PL_MISSING],
                _ => {
                    let mut out = Vec::with_capacity(n_genotypes(nk));
                    for j in 0..nk {
                        for i in 0..=j {
                            out.push(missing(row[gt_index(kept[i] as usize, kept[j] as usize)]));
                        }
                    }
                    out
                }
            }
        })
        .collect()
}

fn ln_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

#[cfg(test)]
#[path = "../../tests/unit/call_mcall.rs"]
mod tests;
