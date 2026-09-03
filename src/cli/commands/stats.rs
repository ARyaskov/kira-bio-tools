//! `bcftools stats` port: same sections and column layout, so plot-vcfstats
//! and MultiQC can consume the output.

use crate::annotate::postproc::RegionFilter;
use crate::cli::args::StatsArgs;
use crate::filter::FilterEngine;
use crate::vcf::alleles::gt_alleles;
use crate::vcf::header::{ContigDict, extract_samples};
use crate::vcf::variant_type::{VT_INDEL, VT_MNP, VT_OTHER, VT_REF, VT_SNP, allele_type, record_type};
use crate::vcf::{UnifiedVcfReader, parse_vcf_line};
use anyhow::{Context, Result, bail};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

const AF_BINS: usize = 100;
const VAF_BINS: usize = 21;

/// One parsed record, independent of the input file.
struct Site {
    rank: usize,
    pos: u32,
    id: String,
    refa: String,
    alts: Vec<String>,
    alt_types: Vec<u32>,
    rtype: u32,
    qual: Option<f64>,
    info: String,
    format: Vec<String>,
    samples: Vec<Vec<String>>,
}

impl Site {
    fn parse(line: &str, contigs: &mut ContigDict) -> Option<Self> {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            return None;
        }
        let pos: u32 = cols[1].parse().ok()?;
        let refa = cols[3].to_string();
        let alts: Vec<String> = if cols[4] == "." || cols[4].is_empty() {
            Vec::new()
        } else {
            cols[4].split(',').map(|s| s.to_string()).collect()
        };
        let alt_types: Vec<u32> = alts.iter().map(|a| allele_type(&refa, a).ty).collect();
        let rtype = record_type(&refa, cols[4]);
        let qual = if cols[5] == "." { None } else { cols[5].parse().ok() };
        let format: Vec<String> = if cols.len() > 8 { cols[8].split(':').map(|s| s.to_string()).collect() } else { Vec::new() };
        let samples: Vec<Vec<String>> = if cols.len() > 9 {
            cols[9..].iter().map(|s| s.split(':').map(|v| v.to_string()).collect()).collect()
        } else {
            Vec::new()
        };
        let rank = contigs.insert(cols[0]) as usize;
        Some(Self {
            rank,
            pos,
            id: cols[2].to_string(),
            refa,
            alts,
            alt_types,
            rtype,
            qual,
            info: cols[7].to_string(),
            format,
            samples,
        })
    }

    fn fmt_idx(&self, key: &str) -> Option<usize> {
        self.format.iter().position(|k| k == key)
    }

    fn indel_len(&self, alt: usize) -> i64 {
        self.alts[alt - 1].len() as i64 - self.refa.len() as i64
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Gt {
    HomRR,
    HetRA,
    HomAA,
    HetAA,
    HaplR,
    HaplA,
    Unknown,
}

/// `bcf_gt_type`: genotype class plus the two allele indices (haploid: ial only).
fn gt_type(gt: &str) -> (Gt, Option<usize>, Option<usize>) {
    let alleles = gt_alleles(gt);
    let called: Vec<usize> = alleles.iter().flatten().copied().collect();
    if called.len() != alleles.len() || alleles.is_empty() {
        return (Gt::Unknown, None, None);
    }
    if called.len() == 1 {
        return if called[0] == 0 { (Gt::HaplR, Some(0), None) } else { (Gt::HaplA, Some(called[0]), None) };
    }
    let (a, b) = (called[0], called[1]);
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    let g = if lo == hi {
        if lo == 0 { Gt::HomRR } else { Gt::HomAA }
    } else if lo == 0 {
        Gt::HetRA
    } else {
        Gt::HetAA
    };
    (g, Some(lo), Some(hi))
}

fn is_ts(r: u8, a: u8) -> bool {
    matches!((r, a), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C'))
}

/// REF/ALT bases of a SNP allele the way bcftools sees them: the first base of
/// each allele; `None` when they are equal or not ACGT.
fn snp_bases(refa: &str, alt: &str) -> Option<(u8, u8)> {
    let r = refa.bytes().next()?.to_ascii_uppercase();
    let a = alt.bytes().next()?.to_ascii_uppercase();
    if !b"ACGT".contains(&r) || !b"ACGT".contains(&a) || r == a {
        return None;
    }
    Some((r, a))
}

/// AF bin index as in bcftools `init_iaf`: singletons and sites without
/// allele counts land in bin 0.
fn af_bin(ac: u32, an: u32) -> usize {
    if an == 0 || ac == 1 {
        return 0;
    }
    ((ac as f64 * (AF_BINS as f64 - 1.0) / an as f64) as usize).min(AF_BINS - 1)
}

/// Bin of ALT allele `i` (0-based), from `--af-tag` when given.
fn alt_bin(s: &Site, ac: &[u32], an: u32, af_tag: Option<&str>, i: usize) -> usize {
    match af_tag {
        Some(tag) => info_af_bin(&s.info, tag, i),
        None => af_bin(ac.get(i + 1).copied().unwrap_or(0), an),
    }
}

/// HWE het-fraction bin; single precision reproduces bcftools' rounding.
fn frac_bin(num: u32, den: u32) -> usize {
    if den == 0 {
        return 0;
    }
    let f = num as f32 / den as f32;
    ((f * (AF_BINS as f32 - 1.0)) as usize).min(AF_BINS - 1)
}

fn vaf_bin(vaf: f64) -> usize {
    ((vaf * (VAF_BINS as f64 - 1.0) + 0.5) as usize).min(VAF_BINS - 1)
}

struct Cfg {
    dp: DpBins,
    af_tag: Option<String>,
    first_only: bool,
    has_ad: bool,
}

#[derive(Default, Clone)]
struct SampleStats {
    ref_hom: u64,
    nonref_hom: u64,
    hets: u64,
    ts: u64,
    tv: u64,
    indels: u64,
    dp_sum: f64,
    dp_n: u64,
    singletons: u64,
    hap_ref: u64,
    hap_alt: u64,
    missing: u64,
    ins_hets: u64,
    del_hets: u64,
    ins_homs: u64,
    del_homs: u64,
    vaf_snv: [u64; VAF_BINS],
    vaf_indel: [u64; VAF_BINS],
}

#[derive(Default, Clone)]
struct TypeCounts {
    snps: u64,
    ts: u64,
    tv: u64,
    indels: u64,
    repeat_consistent: u64,
    repeat_inconsistent: u64,
    not_applicable: u64,
}

impl TypeCounts {
    fn add_snp(&mut self, ts: bool) {
        self.snps += 1;
        if ts { self.ts += 1 } else { self.tv += 1 }
    }

    fn add_indel(&mut self) {
        self.indels += 1;
        self.not_applicable += 1;
    }
}

#[derive(Default)]
struct Stats {
    n_samples: usize,
    /// Samples in the file header; printed for per-file sets.
    n_header_samples: usize,
    n_records: u64,
    n_no_alts: u64,
    n_snps: u64,
    n_mnps: u64,
    n_indels: u64,
    n_others: u64,
    n_multi: u64,
    n_multi_snps: u64,
    ts: u64,
    tv: u64,
    ts_alt1: u64,
    tv_alt1: u64,
    subst: BTreeMap<(u8, u8), u64>,
    singletons: TypeCounts,
    af: BTreeMap<usize, TypeCounts>,
    /// QUAL bin -> (snps, ts alt1, tv alt1, indels); None = missing QUAL
    qual: BTreeMap<Option<i64>, (u64, u64, u64, u64)>,
    /// indel length -> (sites, genotypes, vaf sum)
    idd: BTreeMap<i64, (u64, u64, f64)>,
    dp_gt: BTreeMap<u32, u64>,
    dp_site: BTreeMap<u32, u64>,
    dp_gt_total: u64,
    dp_site_total: u64,
    hwe: BTreeMap<usize, Vec<u64>>,
    samples: Vec<SampleStats>,
    sample_names: Vec<String>,
}

impl Stats {
    fn new(names: &[String], n_header_samples: usize) -> Self {
        Self {
            n_samples: names.len(),
            n_header_samples,
            samples: vec![SampleStats::default(); names.len()],
            sample_names: names.to_vec(),
            ..Default::default()
        }
    }

    fn add(&mut self, s: &Site, selected: &[usize], cfg: &Cfg) {
        self.n_records += 1;
        let rt = s.rtype;
        if s.alts.is_empty() || rt == VT_REF {
            self.n_no_alts += 1;
        }
        if rt & VT_SNP != 0 { self.n_snps += 1; }
        if rt & VT_MNP != 0 { self.n_mnps += 1; }
        if rt & VT_INDEL != 0 { self.n_indels += 1; }
        if rt & VT_OTHER != 0 { self.n_others += 1; }
        if s.alts.len() > 1 {
            self.n_multi += 1;
            if s.alt_types.iter().all(|t| *t == VT_SNP) { self.n_multi_snps += 1; }
        }

        let n_al = s.alts.len() + 1;
        let (ac, an) = allele_counts(s, selected, n_al);
        let af_tag = cfg.af_tag.as_deref();

        // Per-allele site stats; `-1` looks at the first ALT only.
        let mut alt1_ts: Option<bool> = None;
        for (i, a) in s.alts.iter().enumerate() {
            if cfg.first_only && i > 0 {
                break;
            }
            let t = s.alt_types[i];
            let ac_i = ac.get(i + 1).copied().unwrap_or(0);
            let singleton = an > 0 && ac_i == 1;
            let bin = alt_bin(s, &ac, an, af_tag, i);
            if t == VT_SNP {
                let Some((r, b)) = snp_bases(&s.refa, a) else { continue };
                let ts = is_ts(r, b);
                *self.subst.entry((r, b)).or_default() += 1;
                if ts { self.ts += 1 } else { self.tv += 1 }
                if i == 0 {
                    alt1_ts = Some(ts);
                    if ts { self.ts_alt1 += 1 } else { self.tv_alt1 += 1 }
                }
                self.af.entry(bin).or_default().add_snp(ts);
                if singleton { self.singletons.add_snp(ts) }
            } else if t == VT_INDEL {
                self.af.entry(bin).or_default().add_indel();
                if singleton { self.singletons.add_indel() }
                self.idd.entry(s.indel_len(i + 1)).or_insert((0, 0, 0.0)).0 += 1;
            }
        }

        let q = self.qual.entry(s.qual.map(|q| q.floor() as i64)).or_default();
        if let Some(ts) = alt1_ts {
            q.0 += 1;
            if ts { q.1 += 1 } else { q.2 += 1 }
        }
        if rt & VT_INDEL != 0 { q.3 += 1; }

        if let Some(dp) = s.info.split(';').find_map(|kv| kv.strip_prefix("DP=")).and_then(|v| v.parse::<u32>().ok()) {
            *self.dp_site.entry(cfg.dp.bin(dp)).or_default() += 1;
            self.dp_site_total += 1;
        }

        if selected.is_empty() {
            return;
        }
        let gt_idx = s.fmt_idx("GT");
        let dp_idx = s.fmt_idx("DP");
        let ad_idx = s.fmt_idx("AD");
        let (mut n_rr, mut n_ra, mut n_aa) = (0u32, 0u32, 0u32);
        let mut nonref: Vec<usize> = Vec::new();
        for (k, &si) in selected.iter().enumerate() {
            let Some(smp) = s.samples.get(si) else { continue };
            let st = &mut self.samples[k];
            let ad: Option<Vec<f64>> = ad_idx.and_then(|i| smp.get(i)).map(|v| {
                v.split(',').map(|x| x.parse::<f64>().unwrap_or(0.0)).collect()
            });
            let dp = dp_idx
                .and_then(|i| smp.get(i))
                .and_then(|v| v.parse::<u32>().ok())
                .or_else(|| ad.as_ref().map(|a| a.iter().sum::<f64>() as u32));
            if let Some(d) = dp {
                st.dp_sum += d as f64;
                st.dp_n += 1;
                *self.dp_gt.entry(cfg.dp.bin(d)).or_default() += 1;
                self.dp_gt_total += 1;
            }
            let Some(gi) = gt_idx else { continue };
            let (g, ial, jal) = gt_type(smp.get(gi).map(String::as_str).unwrap_or("."));
            match g {
                Gt::Unknown => { st.missing += 1; continue; }
                Gt::HaplR => { st.hap_ref += 1; continue; }
                Gt::HaplA => { st.hap_alt += 1; continue; }
                Gt::HomRR => n_rr += 1,
                Gt::HetRA => n_ra += 1,
                Gt::HomAA => n_aa += 1,
                Gt::HetAA => {}
            }
            if g != Gt::HomRR { nonref.push(k); }
            let mut alleles: Vec<usize> = Vec::with_capacity(2);
            for a in [ial, jal].into_iter().flatten() {
                if a > 0 && a <= s.alts.len() && !alleles.contains(&a) { alleles.push(a); }
            }
            // Class counters follow the sample's own allele types (VT_REF for 0/0).
            let var_type = alleles.iter().fold(VT_REF, |m, &a| m | s.alt_types[a - 1]);
            if var_type & VT_SNP != 0 || var_type == VT_REF {
                match g {
                    Gt::HomRR => st.ref_hom += 1,
                    Gt::HomAA => st.nonref_hom += 1,
                    Gt::HetRA | Gt::HetAA => st.hets += 1,
                    _ => {}
                }
            }
            let is_het = matches!(g, Gt::HetRA | Gt::HetAA);
            let (mut has_ins, mut has_del) = (false, false);
            for &a in &alleles {
                match s.alt_types[a - 1] {
                    VT_SNP => {
                        if let Some((r, b)) = snp_bases(&s.refa, &s.alts[a - 1]) {
                            if is_ts(r, b) { st.ts += 1 } else { st.tv += 1 }
                        }
                    }
                    VT_INDEL => {
                        if s.indel_len(a) > 0 { has_ins = true } else { has_del = true }
                    }
                    _ => {}
                }
            }
            if has_ins || has_del {
                st.indels += 1;
                if is_het {
                    if has_ins { st.ins_hets += 1 }
                    if has_del { st.del_hets += 1 }
                } else {
                    if has_ins { st.ins_homs += 1 }
                    if has_del { st.del_homs += 1 }
                }
            }
            // VAF from AD: per carried allele, SNVs and indels binned separately.
            if let Some(ad) = &ad {
                let tot: f64 = ad.iter().sum();
                if tot > 0.0 {
                    for &a in &alleles {
                        let x = ad.get(a).copied().unwrap_or(0.0);
                        if x <= 0.0 { continue; }
                        let vaf = x / tot;
                        match s.alt_types[a - 1] {
                            VT_SNP => st.vaf_snv[vaf_bin(vaf)] += 1,
                            VT_INDEL => {
                                st.vaf_indel[vaf_bin(vaf)] += 1;
                                let e = self.idd.entry(s.indel_len(a)).or_insert((0, 0, 0.0));
                                e.1 += 1;
                                e.2 += vaf;
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
        if let [k] = nonref[..] {
            self.samples[k].singletons += 1;
        }
        let n_hwe = n_rr + n_ra + n_aa;
        if !s.alts.is_empty() && n_hwe > 0 {
            let bin = alt_bin(s, &ac, an, af_tag, 0);
            let hb = frac_bin(n_ra, n_hwe);
            self.hwe.entry(bin).or_insert_with(|| vec![0; AF_BINS])[hb] += 1;
        }
    }
}

/// Per-allele counts (REF first) and AN: INFO AC/AN when both are present,
/// otherwise every genotype of the record once samples are requested
/// (`bcf_calc_ac` with `BCF_UN_INFO|BCF_UN_FMT`).
fn allele_counts(s: &Site, selected: &[usize], n_al: usize) -> (Vec<u32>, u32) {
    let mut ac = vec![0u32; n_al];
    let mut info_ac: Option<Vec<u32>> = None;
    let mut info_an: Option<u32> = None;
    for kv in s.info.split(';') {
        if let Some(v) = kv.strip_prefix("AC=") {
            info_ac = Some(v.split(',').map(|x| x.parse::<u32>().unwrap_or(0)).collect());
        } else if let Some(v) = kv.strip_prefix("AN=") {
            info_an = v.parse().ok();
        }
    }
    if let (Some(iac), Some(an)) = (info_ac, info_an) {
        for (i, c) in iac.iter().enumerate() {
            if i + 1 < ac.len() { ac[i + 1] = *c; }
        }
        ac[0] = an.saturating_sub(iac.iter().sum());
        return (ac, an);
    }
    let mut an = 0u32;
    if let (Some(gi), false) = (s.fmt_idx("GT"), selected.is_empty()) {
        for smp in &s.samples {
            let gt = smp.get(gi).map(String::as_str).unwrap_or(".");
            for a in gt_alleles(gt).into_iter().flatten() {
                an += 1;
                if a < ac.len() { ac[a] += 1; }
            }
        }
    }
    (ac, an)
}

fn info_af_bin(info: &str, tag: &str, alt_idx: usize) -> usize {
    for kv in info.split(';') {
        if let Some(v) = kv.strip_prefix(tag).and_then(|r| r.strip_prefix('=')) {
            if let Some(af) = v.split(',').nth(alt_idx).and_then(|x| x.parse::<f64>().ok()) {
                return ((af * (AF_BINS as f64 - 1.0)) as usize).min(AF_BINS - 1);
            }
        }
    }
    0
}

struct DpBins {
    min: u32,
    max: u32,
    step: u32,
}

impl DpBins {
    fn parse(s: &str) -> Result<Self> {
        let v: Vec<u32> = s.split(',').map(|x| x.trim().parse::<u32>()).collect::<Result<_, _>>().context("-d/--depth")?;
        if v.len() != 3 || v[2] == 0 { bail!("--depth expects min,max,step"); }
        Ok(Self { min: v[0], max: v[1], step: v[2] })
    }

    /// Bin value printed in the DP section; depths beyond `max` collapse into `max+step`.
    fn bin(&self, dp: u32) -> u32 {
        if dp > self.max {
            return self.max + self.step;
        }
        let d = dp.max(self.min);
        self.min + ((d - self.min) / self.step) * self.step
    }
}

/// Genotype concordance between two files at shared sites.
#[derive(Default)]
struct Concordance {
    /// per AF bin: [rr_m, ra_m, aa_m, rr_x, ra_x, aa_x, n]
    by_af: BTreeMap<usize, [u64; 7]>,
    by_af_dosage: BTreeMap<usize, (f64, f64, f64, f64, f64, u64)>,
    /// per sample: [rr_m, ra_m, aa_m, rr_x, ra_x, aa_x]
    by_sample: Vec<[u64; 6]>,
    by_sample_dosage: Vec<(f64, f64, f64, f64, f64, u64)>,
    /// per sample 5x5 transition table: [from][to] with classes RR, RA, AA, AAhet, missing
    table: Vec<[[u64; 5]; 5]>,
}

impl Concordance {
    fn new(n: usize) -> Self {
        Self {
            by_sample: vec![[0; 6]; n],
            by_sample_dosage: vec![(0.0, 0.0, 0.0, 0.0, 0.0, 0); n],
            table: vec![[[0; 5]; 5]; n],
            ..Default::default()
        }
    }

    fn class(g: Gt) -> usize {
        match g {
            Gt::HomRR => 0,
            Gt::HetRA => 1,
            Gt::HomAA => 2,
            Gt::HetAA => 3,
            _ => 4,
        }
    }

    fn dosage(g: Gt) -> Option<f64> {
        match g {
            Gt::HomRR => Some(0.0),
            Gt::HetRA => Some(1.0),
            Gt::HomAA | Gt::HetAA => Some(2.0),
            _ => None,
        }
    }

    fn add(&mut self, a: &Site, b: &Site, sel_a: &[usize], sel_b: &[usize], af_bin: usize) {
        let (Some(ga), Some(gb)) = (a.fmt_idx("GT"), b.fmt_idx("GT")) else { return };
        for (k, (&sa, &sb)) in sel_a.iter().zip(sel_b.iter()).enumerate() {
            let gta = a.samples.get(sa).and_then(|s| s.get(ga)).map(String::as_str).unwrap_or(".");
            let gtb = b.samples.get(sb).and_then(|s| s.get(gb)).map(String::as_str).unwrap_or(".");
            let (ta, ..) = gt_type(gta);
            let (tb, ..) = gt_type(gtb);
            self.table[k][Self::class(ta)][Self::class(tb)] += 1;
            if matches!(ta, Gt::Unknown | Gt::HaplR | Gt::HaplA) || matches!(tb, Gt::Unknown | Gt::HaplR | Gt::HaplA) {
                continue;
            }
            let ca = Self::class(ta).min(2);
            let matched = ta == tb;
            let e = self.by_af.entry(af_bin).or_default();
            e[6] += 1;
            let idx = if matched { ca } else { 3 + ca };
            e[idx] += 1;
            self.by_sample[k][idx] += 1;
            if let (Some(da), Some(db)) = (Self::dosage(ta), Self::dosage(tb)) {
                let d = self.by_af_dosage.entry(af_bin).or_insert((0.0, 0.0, 0.0, 0.0, 0.0, 0));
                d.0 += da; d.1 += db; d.2 += da * da; d.3 += db * db; d.4 += da * db; d.5 += 1;
                let s = &mut self.by_sample_dosage[k];
                s.0 += da; s.1 += db; s.2 += da * da; s.3 += db * db; s.4 += da * db; s.5 += 1;
            }
        }
    }
}

/// Dosage r²; `None` when either file has no dosage variance.
fn r_squared(d: &(f64, f64, f64, f64, f64, u64)) -> Option<f64> {
    let n = d.5 as f64;
    if n < 2.0 { return None; }
    let cov = d.4 / n - (d.0 / n) * (d.1 / n);
    let va = d.2 / n - (d.0 / n).powi(2);
    let vb = d.3 / n - (d.1 / n).powi(2);
    if va <= 1e-12 || vb <= 1e-12 { return None; }
    Some((cov * cov / (va * vb)).clamp(0.0, 1.0))
}

fn fmt_r2(r2: Option<f64>) -> String {
    match r2 {
        Some(v) => format!("{:.6}", v),
        None => "0".to_string(),
    }
}

struct Input {
    reader: UnifiedVcfReader,
    next: Option<Site>,
    selected: Vec<usize>,
    region: Option<RegionFilter>,
    regions_overlap: u8,
    target: Option<RegionFilter>,
    target_inverse: bool,
    apply_filters: Option<Vec<String>>,
    include: Option<FilterEngine>,
    exclude: Option<FilterEngine>,
}

impl Input {
    fn advance(&mut self, contigs: &mut ContigDict) -> Result<()> {
        loop {
            let Some(line) = self.reader.read_line()? else { self.next = None; return Ok(()); };
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            if let Some(rf) = &self.region {
                if !rf.line_passes_mode(&line, self.regions_overlap) { continue; }
            }
            if let Some(tf) = &self.target {
                if tf.line_passes_mode(&line, 0) == self.target_inverse { continue; }
            }
            if let Some(af) = &self.apply_filters {
                let f = line.splitn(8, '\t').nth(6).unwrap_or(".");
                let pass = if f == "." || f.is_empty() { af.iter().any(|a| a == ".") } else { f.split(';').any(|t| af.iter().any(|a| a == t)) };
                if !pass { continue; }
            }
            if self.include.is_some() || self.exclude.is_some() {
                if let Some(rec) = parse_vcf_line(&line) {
                    if let Some(e) = &self.include {
                        if !e.eval(&rec).map(|r| r.pass_site).unwrap_or(true) { continue; }
                    }
                    if let Some(e) = &self.exclude {
                        if e.eval(&rec).map(|r| r.pass_site).unwrap_or(false) { continue; }
                    }
                }
            }
            if let Some(site) = Site::parse(&line, contigs) {
                self.next = Some(site);
                return Ok(());
            }
        }
    }

    fn take(&mut self, contigs: &mut ContigDict) -> Result<Option<Site>> {
        let s = self.next.take();
        if s.is_some() {
            self.advance(contigs)?;
        }
        Ok(s)
    }
}

pub fn cmd_stats(args: StatsArgs) -> Result<()> {
    if args.inputs.is_empty() { bail!("stats: at least one input required"); }
    if args.inputs.len() > 2 { bail!("stats: at most two input files can be compared"); }
    let region = if let Some(s) = &args.regions {
        Some(RegionFilter::from_cli(s)?)
    } else if let Some(p) = &args.regions_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };
    let target = if let Some(s) = &args.targets {
        Some(RegionFilter::from_cli(s.trim_start_matches('^'))?)
    } else if let Some(p) = &args.targets_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };
    let target_inverse = args.targets.as_deref().is_some_and(|s| s.starts_with('^'));
    let apply_filters: Option<Vec<String>> = args.apply_filters.as_deref().map(|s| s.split(',').map(|t| t.trim().to_string()).collect());

    let mut contigs = ContigDict::new();
    let mut inputs: Vec<Input> = Vec::new();
    let mut names: Vec<Vec<String>> = Vec::new();
    let mut header_n: Vec<usize> = Vec::new();
    let mut has_ad = false;
    for p in &args.inputs {
        let r = UnifiedVcfReader::open(p).with_context(|| format!("open {}", p.display()))?;
        let headers = r.header()?;
        has_ad |= headers.iter().any(|h| h.starts_with("##FORMAT=<ID=AD,"));
        let all = extract_samples(&headers);
        header_n.push(all.len());
        let selected = select_samples(&all, args.samples.as_deref(), args.samples_file.as_deref())?;
        let include = args.include.as_deref().map(|e| FilterEngine::new(&headers, Some(e), false)).transpose().context("-i")?;
        let exclude = args.exclude.as_deref().map(|e| FilterEngine::new(&headers, Some(e), false)).transpose().context("-e")?;
        names.push(selected.iter().map(|&i| all[i].clone()).collect());
        let mut inp = Input {
            reader: r,
            next: None,
            selected,
            region: region.clone(),
            regions_overlap: args.regions_overlap,
            target: target.clone(),
            target_inverse,
            apply_filters: apply_filters.clone(),
            include,
            exclude,
        };
        inp.advance(&mut contigs)?;
        inputs.push(inp);
    }
    let cfg = Cfg {
        dp: DpBins::parse(&args.depth)?,
        af_tag: args.af_tag.clone(),
        first_only: args.first_allele_only,
        has_ad,
    };

    // Sets: one per file, plus (two files) the intersection, or (-I) novel/known.
    let two = inputs.len() == 2;
    let split_id = args.split_by_id && !two;
    let file_name = |i: usize| args.inputs[i].display().to_string();
    let mut sets: Vec<Stats> = Vec::new();
    let mut set_names: Vec<Vec<String>> = Vec::new();
    if two {
        // The shared set uses the samples common to both files.
        let common: Vec<String> = names[0].iter().filter(|n| names[1].contains(n)).cloned().collect();
        sets.push(Stats::new(&names[0], header_n[0]));
        sets.push(Stats::new(&names[1], header_n[1]));
        sets.push(Stats::new(&common, common.len()));
        set_names.push(vec![file_name(0)]);
        set_names.push(vec![file_name(1)]);
        set_names.push(vec![file_name(0), file_name(1)]);
    } else if split_id {
        sets.push(Stats::new(&names[0], header_n[0]));
        sets.push(Stats::new(&names[0], header_n[0]));
        set_names.push(vec![file_name(0)]);
        set_names.push(vec![file_name(0)]);
    } else {
        sets.push(Stats::new(&names[0], header_n[0]));
        set_names.push(vec![file_name(0)]);
    }
    let n_common = sets.last().map(|s| s.n_samples).unwrap_or(0);
    let mut conc_snps = Concordance::new(n_common);
    let mut conc_indels = Concordance::new(n_common);
    let (sel_a, sel_b): (Vec<usize>, Vec<usize>) = if two {
        let common = &sets[2].sample_names;
        let pick = |f: usize, n: &String| inputs[f].selected[names[f].iter().position(|x| x == n).unwrap()];
        (common.iter().map(|n| pick(0, n)).collect(), common.iter().map(|n| pick(1, n)).collect())
    } else {
        (Vec::new(), Vec::new())
    };
    let sel0 = inputs[0].selected.clone();
    let sel1 = inputs.get(1).map(|i| i.selected.clone()).unwrap_or_default();

    if !two {
        while let Some(site) = inputs[0].take(&mut contigs)? {
            let set = if split_id && site.id != "." { 1 } else { 0 };
            sets[set].add(&site, &sel0, &cfg);
        }
    } else {
        loop {
            let ka = inputs[0].next.as_ref().map(|s| (s.rank, s.pos));
            let kb = inputs[1].next.as_ref().map(|s| (s.rank, s.pos));
            let (Some(a), Some(b)) = (ka, kb) else {
                if ka.is_some() {
                    let s = inputs[0].take(&mut contigs)?.unwrap();
                    sets[0].add(&s, &sel0, &cfg);
                    continue;
                }
                if kb.is_some() {
                    let s = inputs[1].take(&mut contigs)?.unwrap();
                    sets[1].add(&s, &sel1, &cfg);
                    continue;
                }
                break;
            };
            if a < b {
                let s = inputs[0].take(&mut contigs)?.unwrap();
                sets[0].add(&s, &sel0, &cfg);
                continue;
            }
            if b < a {
                let s = inputs[1].take(&mut contigs)?.unwrap();
                sets[1].add(&s, &sel1, &cfg);
                continue;
            }
            // Same position: pair records by alleles, the rest stay file-private.
            let mut ra: Vec<Site> = Vec::new();
            let mut rb: Vec<Site> = Vec::new();
            while inputs[0].next.as_ref().is_some_and(|s| (s.rank, s.pos) == a) {
                ra.push(inputs[0].take(&mut contigs)?.unwrap());
            }
            while inputs[1].next.as_ref().is_some_and(|s| (s.rank, s.pos) == b) {
                rb.push(inputs[1].take(&mut contigs)?.unwrap());
            }
            let mut used_b = vec![false; rb.len()];
            for sa in &ra {
                let m = rb.iter().enumerate().position(|(j, sb)| {
                    !used_b[j]
                        && sa.refa.eq_ignore_ascii_case(&sb.refa)
                        && sa.alts.len() == sb.alts.len()
                        && sa.alts.iter().zip(sb.alts.iter()).all(|(x, y)| x.eq_ignore_ascii_case(y))
                });
                match m {
                    Some(j) => {
                        used_b[j] = true;
                        let sb = &rb[j];
                        sets[2].add(sa, &sel_a, &cfg);
                        let (ac, an) = allele_counts(sa, &sel_a, sa.alts.len() + 1);
                        let bin = alt_bin(sa, &ac, an, cfg.af_tag.as_deref(), 0);
                        if sa.rtype & VT_SNP != 0 {
                            conc_snps.add(sa, sb, &sel_a, &sel_b, bin);
                        } else if sa.rtype & VT_INDEL != 0 {
                            conc_indels.add(sa, sb, &sel_a, &sel_b, bin);
                        }
                    }
                    None => sets[0].add(sa, &sel0, &cfg),
                }
            }
            for (j, sb) in rb.iter().enumerate() {
                if !used_b[j] {
                    sets[1].add(sb, &sel1, &cfg);
                }
            }
        }
    }

    print_stats(&args, &sets, &set_names, &cfg, if two { Some((&conc_snps, &conc_indels)) } else { None })
}

fn select_samples(all: &[String], cli: Option<&str>, file: Option<&Path>) -> Result<Vec<usize>> {
    let mut names: Vec<String> = Vec::new();
    match cli {
        Some("-") => return Ok((0..all.len()).collect()),
        Some(s) => {
            let (inv, body) = match s.strip_prefix('^') { Some(b) => (true, b), None => (false, s) };
            for t in body.split(',') {
                let t = t.trim();
                if !t.is_empty() { names.push(t.to_string()); }
            }
            if inv {
                return Ok(all.iter().enumerate().filter(|(_, n)| !names.contains(n)).map(|(i, _)| i).collect());
            }
        }
        None => {}
    }
    if let Some(p) = file {
        let text = std::fs::read_to_string(p).with_context(|| format!("-S {}", p.display()))?;
        for l in text.lines() {
            let t = l.trim();
            if !t.is_empty() && !t.starts_with('#') { names.push(t.split_whitespace().next().unwrap_or("").to_string()); }
        }
    }
    if names.is_empty() && file.is_none() && cli.is_none() {
        return Ok(Vec::new());
    }
    let mut out = Vec::new();
    for n in &names {
        match all.iter().position(|s| s == n) {
            Some(i) => out.push(i),
            None => bail!("sample {n:?} not found"),
        }
    }
    Ok(out)
}

fn print_stats(args: &StatsArgs, sets: &[Stats], set_names: &[Vec<String>], cfg: &Cfg, conc: Option<(&Concordance, &Concordance)>) -> Result<()> {
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::with_capacity(64 * 1024, stdout.lock());
    let cmd: Vec<String> = std::env::args().skip(2).filter(|a| a != "--").collect();

    writeln!(out, "# This file was produced by kira-bt stats ({}) and can be plotted using plot-vcfstats.", env!("CARGO_PKG_VERSION"))?;
    writeln!(out, "# The command line was:\tkira-bt stats {}", cmd.join(" "))?;
    writeln!(out, "#")?;
    writeln!(out, "# Definition of sets:")?;
    writeln!(out, "# ID\t[2]id\t[3]tab-separated file names")?;
    for (i, n) in set_names.iter().enumerate() {
        writeln!(out, "ID\t{}\t{}", i, n.join("\t"))?;
    }
    writeln!(out, "# SN, Summary numbers:")?;
    writeln!(out, "#   number of records   .. number of data rows in the VCF")?;
    writeln!(out, "#   number of no-ALTs   .. reference-only sites, ALT is either \".\" or identical to REF")?;
    writeln!(out, "#   number of SNPs      .. number of rows with a SNP")?;
    writeln!(out, "#   number of MNPs      .. number of rows with a MNP, such as CC>TT")?;
    writeln!(out, "#   number of indels    .. number of rows with an indel")?;
    writeln!(out, "#   number of others    .. number of rows with other type, for example a symbolic allele or")?;
    writeln!(out, "#                          a complex substitution, such as ACT>TCGA")?;
    writeln!(out, "#   number of multiallelic sites     .. number of rows with multiple alternate alleles")?;
    writeln!(out, "#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs")?;
    writeln!(out, "# ")?;
    writeln!(out, "#   Note that rows containing multiple types will be counted multiple times, in each")?;
    writeln!(out, "#   counter. For example, a row with a SNP and an indel increments both the SNP and")?;
    writeln!(out, "#   the indel counter.")?;
    writeln!(out, "# ")?;
    writeln!(out, "# SN\t[2]id\t[3]key\t[4]value")?;
    let has_samples = args.samples.is_some() || args.samples_file.is_some();
    if has_samples {
        for (i, s) in sets.iter().enumerate() {
            if sets.len() == 3 && i == 2 { continue; }
            writeln!(out, "SN\t{}\tnumber of samples:\t{}", i, s.n_header_samples)?;
        }
    }
    for (i, s) in sets.iter().enumerate() {
        writeln!(out, "SN\t{}\tnumber of records:\t{}", i, s.n_records)?;
        writeln!(out, "SN\t{}\tnumber of no-ALTs:\t{}", i, s.n_no_alts)?;
        writeln!(out, "SN\t{}\tnumber of SNPs:\t{}", i, s.n_snps)?;
        writeln!(out, "SN\t{}\tnumber of MNPs:\t{}", i, s.n_mnps)?;
        writeln!(out, "SN\t{}\tnumber of indels:\t{}", i, s.n_indels)?;
        writeln!(out, "SN\t{}\tnumber of others:\t{}", i, s.n_others)?;
        writeln!(out, "SN\t{}\tnumber of multiallelic sites:\t{}", i, s.n_multi)?;
        writeln!(out, "SN\t{}\tnumber of multiallelic SNP sites:\t{}", i, s.n_multi_snps)?;
    }

    writeln!(out, "# TSTV, transitions/transversions:")?;
    writeln!(out, "# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)")?;
    for (i, s) in sets.iter().enumerate() {
        writeln!(out, "TSTV\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i, s.ts, s.tv, ratio(s.ts, s.tv), s.ts_alt1, s.tv_alt1, ratio(s.ts_alt1, s.tv_alt1))?;
    }

    writeln!(out, "# SiS, Singleton stats:")?;
    writeln!(out, "# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")?;
    for (i, s) in sets.iter().enumerate() {
        let c = &s.singletons;
        writeln!(out, "SiS\t{}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i, c.snps, c.ts, c.tv, c.indels, c.repeat_consistent, c.repeat_inconsistent, c.not_applicable)?;
    }

    writeln!(out, "# AF, Stats by non-reference allele frequency:")?;
    writeln!(out, "# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")?;
    for (i, s) in sets.iter().enumerate() {
        for (bin, c) in &s.af {
            if c.snps == 0 && c.indels == 0 { continue; }
            writeln!(out, "AF\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", i, *bin as f64 / AF_BINS as f64, c.snps, c.ts, c.tv, c.indels, c.repeat_consistent, c.repeat_inconsistent, c.not_applicable)?;
        }
    }

    writeln!(out, "# QUAL, Stats by quality")?;
    writeln!(out, "# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels")?;
    for (i, s) in sets.iter().enumerate() {
        for (q, v) in &s.qual {
            if v.0 == 0 && v.3 == 0 { continue; }
            let qs = q.map(|x| x.to_string()).unwrap_or_else(|| ".".into());
            writeln!(out, "QUAL\t{}\t{}\t{}\t{}\t{}\t{}", i, qs, v.0, v.1, v.2, v.3)?;
        }
    }

    writeln!(out, "# IDD, InDel distribution:")?;
    writeln!(out, "# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF")?;
    for (i, s) in sets.iter().enumerate() {
        for (len, (sites, gts, vaf)) in &s.idd {
            let mean = if *gts > 0 { format!("{:.2}", vaf / *gts as f64) } else { ".".to_string() };
            writeln!(out, "IDD\t{}\t{}\t{}\t{}\t{}", i, len, sites, gts, mean)?;
        }
    }

    writeln!(out, "# ST, Substitution types:")?;
    writeln!(out, "# ST\t[2]id\t[3]type\t[4]count")?;
    let types = [
        (b'A', b'C'), (b'A', b'G'), (b'A', b'T'),
        (b'C', b'A'), (b'C', b'G'), (b'C', b'T'),
        (b'G', b'A'), (b'G', b'C'), (b'G', b'T'),
        (b'T', b'A'), (b'T', b'C'), (b'T', b'G'),
    ];
    for (i, s) in sets.iter().enumerate() {
        for (r, a) in &types {
            let c = s.subst.get(&(*r, *a)).copied().unwrap_or(0);
            writeln!(out, "ST\t{}\t{}>{}\t{}", i, *r as char, *a as char, c)?;
        }
    }

    if let Some((cs, ci)) = conc {
        if has_samples {
            writeln!(out, "SN\t2\tnumber of samples:\t{}", sets[2].n_samples)?;
        }
        print_concordance(&mut out, cs, ci, &sets[2].sample_names)?;
    }

    writeln!(out, "# DP, Depth distribution")?;
    writeln!(out, "# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)")?;
    for (i, s) in sets.iter().enumerate() {
        let mut bins: Vec<u32> = s.dp_gt.keys().chain(s.dp_site.keys()).copied().collect();
        bins.sort_unstable();
        bins.dedup();
        for bin in bins {
            let g = s.dp_gt.get(&bin).copied().unwrap_or(0);
            let st = s.dp_site.get(&bin).copied().unwrap_or(0);
            let gf = if s.dp_gt_total > 0 { 100.0 * g as f64 / s.dp_gt_total as f64 } else { 0.0 };
            let sf = if s.dp_site_total > 0 { 100.0 * st as f64 / s.dp_site_total as f64 } else { 0.0 };
            let label = if bin > cfg.dp.max { format!(">{}", cfg.dp.max) } else { bin.to_string() };
            writeln!(out, "DP\t{}\t{}\t{}\t{:.6}\t{}\t{:.6}", i, label, g, gf, st, sf)?;
        }
    }

    writeln!(out, "# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.")?;
    writeln!(out, "# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing")?;
    for (i, s) in sets.iter().enumerate() {
        for (j, name) in s.sample_names.iter().enumerate() {
            let st = &s.samples[j];
            let avg = if st.dp_n > 0 { st.dp_sum / st.dp_n as f64 } else { 0.0 };
            writeln!(out, "PSC\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
                i, name, st.ref_hom, st.nonref_hom, st.hets, st.ts, st.tv, st.indels, avg, st.singletons, st.hap_ref, st.hap_alt, st.missing)?;
        }
    }

    writeln!(out, "# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del allele are counted twice, in both nInsHets and nDelHets.")?;
    writeln!(out, "# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms")?;
    for (i, s) in sets.iter().enumerate() {
        for (j, name) in s.sample_names.iter().enumerate() {
            let st = &s.samples[j];
            writeln!(out, "PSI\t{}\t{}\t0\t0\t0\t0.00\t{}\t{}\t{}\t{}", i, name, st.ins_hets, st.del_hets, st.ins_homs, st.del_homs)?;
        }
    }

    writeln!(out, "# HWE")?;
    writeln!(out, "# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile")?;
    for (i, s) in sets.iter().enumerate() {
        for (bin, v) in &s.hwe {
            let n: u64 = v.iter().sum();
            if n == 0 { continue; }
            let pct = |p: f64| -> f64 {
                let target = (n as f64 * p).ceil().max(1.0) as u64;
                let mut acc = 0u64;
                for (k, c) in v.iter().enumerate() {
                    acc += c;
                    if acc >= target { return k as f64 / AF_BINS as f64; }
                }
                0.0
            };
            writeln!(out, "HWE\t{}\t{:.6}\t{}\t{:.6}\t{:.6}\t{:.6}", i, *bin as f64 / AF_BINS as f64, n, pct(0.25), pct(0.5), pct(0.75))?;
        }
    }

    if cfg.has_ad {
        writeln!(out, "# VAF, Variant Allele Frequency determined as fraction of alternate reads in FORMAT/AD")?;
        writeln!(out, "# VAF\t[2]id\t[3]sample\t[4]SNV VAF distribution\t[5]indel VAF distribution")?;
        for (i, s) in sets.iter().enumerate() {
            for (j, name) in s.sample_names.iter().enumerate() {
                let st = &s.samples[j];
                let snv: Vec<String> = st.vaf_snv.iter().map(u64::to_string).collect();
                let ind: Vec<String> = st.vaf_indel.iter().map(u64::to_string).collect();
                writeln!(out, "VAF\t{}\t{}\t{}\t{}", i, name, snv.join(","), ind.join(","))?;
            }
        }
    }

    out.flush()?;
    Ok(())
}

fn ratio(a: u64, b: u64) -> String {
    if b == 0 { "0.00".to_string() } else { format!("{:.2}", a as f64 / b as f64) }
}

/// NRD and the three per-class discordance rates, in percent.
fn nrd(m: &[u64; 6]) -> (f64, f64, f64, f64) {
    let pct = |x: u64, tot: u64| if tot > 0 { 100.0 * x as f64 / tot as f64 } else { 0.0 };
    let nrd = pct(m[3] + m[4] + m[5], m[3] + m[4] + m[5] + m[1] + m[2]);
    (nrd, pct(m[3], m[3] + m[0]), pct(m[4], m[4] + m[1]), pct(m[5], m[5] + m[2]))
}

/// Two-file concordance sections in bcftools order: GC*AF/NRD*, GC*S, GCT*.
fn print_concordance<W: Write>(out: &mut W, cs: &Concordance, ci: &Concordance, samples: &[String]) -> Result<()> {
    let both = [(cs, "s", "SNPs"), (ci, "i", "indels")];
    for (c, tag, what) in both {
        writeln!(out, "# GC{tag}AF, Genotype concordance by non-reference allele frequency ({what})")?;
        writeln!(out, "# GC{tag}AF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t[11]number of genotypes")?;
        for (bin, v) in &c.by_af {
            // Bin 0 (singletons) is accumulated but never printed here.
            if *bin == 0 { continue; }
            let r2 = c.by_af_dosage.get(bin).and_then(r_squared);
            writeln!(out, "GC{tag}AF\t2\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", *bin as f64 / AF_BINS as f64, v[0], v[1], v[2], v[3], v[4], v[5], fmt_r2(r2), v[6])?;
        }
        if tag == "s" {
            writeln!(out, "# NRD and discordance is calculated as follows:")?;
            writeln!(out, "#   m .. number of matches")?;
            writeln!(out, "#   x .. number of mismatches")?;
            writeln!(out, "#   NRD = 100 * (xRR + xRA + xAA) / (xRR + xRA + xAA + mRA + mAA)")?;
            writeln!(out, "#   RR discordance = 100 * xRR / (xRR + mRR)")?;
            writeln!(out, "#   RA discordance = 100 * xRA / (xRA + mRA)")?;
            writeln!(out, "#   AA discordance = 100 * xAA / (xAA + mAA)")?;
        }
        writeln!(out, "# Non-Reference Discordance (NRD), {what}")?;
        writeln!(out, "# NRD{tag}\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance")?;
        let sum: [u64; 6] = c.by_sample.iter().fold([0; 6], |mut a, s| { for k in 0..6 { a[k] += s[k]; } a });
        let (a, b, cc, d) = nrd(&sum);
        writeln!(out, "NRD{tag}\t2\t{:.6}\t{:.6}\t{:.6}\t{:.6}", a, b, cc, d)?;
    }
    for (c, tag, what) in both {
        writeln!(out, "# GC{tag}S, Genotype concordance by sample ({what})")?;
        writeln!(out, "# GC{tag}S\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\t[11]dosage r-squared")?;
        for (k, name) in samples.iter().enumerate() {
            let m = &c.by_sample[k];
            let (n, ..) = nrd(m);
            let r2 = r_squared(&c.by_sample_dosage[k]);
            writeln!(out, "GC{tag}S\t2\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", name, n, m[0], m[1], m[2], m[3], m[4], m[5], fmt_r2(r2))?;
        }
    }
    let cls = ["RR Hom", "RA Het", "AA Hom", "AA Het", "missing"];
    for (c, tag, what) in both {
        let mut hdr = format!("# GCT{tag}\t[2]sample");
        let mut col = 3;
        for from in cls {
            for to in cls {
                hdr.push_str(&format!("\t[{col}]{from} -> {to}"));
                col += 1;
            }
        }
        writeln!(out, "# GCT{tag}, Genotype concordance table ({what})")?;
        writeln!(out, "{hdr}")?;
        for (k, name) in samples.iter().enumerate() {
            let t = &c.table[k];
            let cells: Vec<String> = t.iter().flat_map(|row| row.iter().map(|v| v.to_string())).collect();
            writeln!(out, "GCT{tag}\t{}\t{}", name, cells.join("\t"))?;
        }
    }
    Ok(())
}
