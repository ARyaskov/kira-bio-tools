use crate::cli::args::StatsArgs;
use crate::vcf::UnifiedVcfReader;
use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::path::Path;

pub fn cmd_stats(args: StatsArgs) -> Result<()> {
    if args.inputs.is_empty() { anyhow::bail!("stats: at least one input required"); }
    let files: Vec<&Path> = args.inputs.iter().map(|p| p.as_path()).collect();
    let mut sets: Vec<Stats> = Vec::with_capacity(files.len());
    for f in &files {
        sets.push(compute_stats(f).with_context(|| format!("stats {:?}", f))?);
    }
    print_stats(&files, &sets, &args)
}

#[derive(Default)]
struct Stats {
    n_records: u64,
    n_no_alts: u64,
    n_snps: u64,
    n_mnps: u64,
    n_indels: u64,
    n_other: u64,
    n_multi: u64,
    n_multi_snps: u64,
    n_singleton: u64,
    samples: Vec<String>,
    qual_dist: BTreeMap<u32, u64>,
    af_dist: BTreeMap<u32, u64>,
    sub_counts: BTreeMap<(u8, u8), u64>,
    indel_dist: BTreeMap<i32, u64>,
    psc_n_snps: Vec<u64>,
    psc_n_indels: Vec<u64>,
    psc_n_het: Vec<u64>,
    psc_n_hom: Vec<u64>,
    psc_n_ref_hom: Vec<u64>,
    psc_n_miss: Vec<u64>,
    psc_n_singletons: Vec<u64>,
    psc_dp_sum: Vec<u64>,
    psc_dp_count: Vec<u64>,
    psi_n_in_frame: Vec<u64>,
    psi_n_out_frame: Vec<u64>,
    psi_n_not_applicable: Vec<u64>,
    psi_n_indels: Vec<u64>,
    psi_n_indel_hets: Vec<u64>,
    psi_n_indel_alts: Vec<u64>,
    sis_n_snps: u64,
    sis_n_ts: u64,
    sis_n_tv: u64,
    sis_n_indels: u64,
    sis_n_multi_snps: u64,
    sis_n_repeat_consistent: u64,
    sis_n_repeat_inconsistent: u64,
    sis_n_not_applicable: u64,
    hwe_bins: BTreeMap<u32, (u64, u64)>,
    dp_bins: BTreeMap<u32, u64>,
}

fn is_transition(r: u8, a: u8) -> bool {
    matches!((r, a), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C'))
}

fn compute_stats(path: &Path) -> Result<Stats> {
    let mut reader = UnifiedVcfReader::open(path)?;
    let headers = reader.header()?;
    let mut s = Stats::default();
    s.samples = extract_samples(&headers);
    let n = s.samples.len();
    s.psc_n_snps = vec![0; n];
    s.psc_n_indels = vec![0; n];
    s.psc_n_het = vec![0; n];
    s.psc_n_hom = vec![0; n];
    s.psc_n_ref_hom = vec![0; n];
    s.psc_n_miss = vec![0; n];
    s.psc_n_singletons = vec![0; n];
    s.psc_dp_sum = vec![0; n];
    s.psc_dp_count = vec![0; n];
    s.psi_n_in_frame = vec![0; n];
    s.psi_n_out_frame = vec![0; n];
    s.psi_n_not_applicable = vec![0; n];
    s.psi_n_indels = vec![0; n];
    s.psi_n_indel_hets = vec![0; n];
    s.psi_n_indel_alts = vec![0; n];

    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 { continue; }
        s.n_records += 1;

        let refa = cols[3]; let alt = cols[4];
        let alts: Vec<&str> = alt.split(',').filter(|a| !a.is_empty() && *a != ".").collect();
        if alts.is_empty() { s.n_no_alts += 1; continue; }
        if alts.len() > 1 { s.n_multi += 1; }
        let mut had_snp = false; let mut had_indel = false; let mut had_mnp = false; let mut had_other = false;
        for a in &alts {
            if refa.len() == 1 && a.len() == 1 {
                had_snp = true;
                let r = refa.as_bytes()[0].to_ascii_uppercase();
                let aa = a.as_bytes()[0].to_ascii_uppercase();
                *s.sub_counts.entry((r, aa)).or_default() += 1;
            } else if refa.len() == a.len() && refa.len() > 1 { had_mnp = true; }
            else if refa.len() != a.len() {
                had_indel = true;
                let d = a.len() as i32 - refa.len() as i32;
                *s.indel_dist.entry(d).or_default() += 1;
            } else { had_other = true; }
        }
        if had_snp { s.n_snps += 1; if alts.len() > 1 { s.n_multi_snps += 1; } }
        if had_mnp { s.n_mnps += 1; }
        if had_indel { s.n_indels += 1; }
        if had_other { s.n_other += 1; }

        if let Ok(q) = cols[5].parse::<f64>() {
            let bucket = if q.is_finite() { (q.max(0.0).min(999.0)).round() as u32 } else { 0 };
            *s.qual_dist.entry(bucket).or_default() += 1;
        }
        for kv in cols[7].split(';') {
            if let Some(v) = kv.strip_prefix("AF=") {
                if let Ok(af) = v.split(',').next().unwrap_or("0").parse::<f64>() {
                    let bucket = ((af * 1000.0).round() as u32).min(1000);
                    *s.af_dist.entry(bucket).or_default() += 1;
                }
                break;
            }
        }
        if had_snp {
            s.sis_n_snps += 1;
            if alts.len() > 1 { s.sis_n_multi_snps += 1; }
            for a in &alts {
                if refa.len() == 1 && a.len() == 1 {
                    let r = refa.as_bytes()[0].to_ascii_uppercase();
                    let aa = a.as_bytes()[0].to_ascii_uppercase();
                    if is_transition(r, aa) { s.sis_n_ts += 1; } else { s.sis_n_tv += 1; }
                }
            }
        }
        if had_indel { s.sis_n_indels += 1; }

        if cols.len() > 8 && !s.samples.is_empty() {
            let fmt = cols[8];
            let fmt_keys: Vec<&str> = fmt.split(':').collect();
            let gt_idx = fmt_keys.iter().position(|k| *k == "GT");
            let dp_idx = fmt_keys.iter().position(|k| *k == "DP");
            let mut alt_alleles_total = 0u32;
            let mut alt_carrier: Vec<usize> = Vec::new();
            let mut gts_collected: Vec<(usize, Vec<u32>)> = Vec::new();
            let mut n_ref_alleles = 0u32;
            let mut n_alt_alleles = 0u32;
            let mut record_dp_sum = 0u64;
            if let Some(gi) = gt_idx {
                let is_snp = had_snp; let is_indel = had_indel;
                for (si, raw) in cols[9..].iter().enumerate() {
                    if si >= s.samples.len() { break; }
                    let parts: Vec<&str> = raw.split(':').collect();
                    let gt = parts.get(gi).copied().unwrap_or(".");
                    if let Some(di) = dp_idx {
                        if let Some(d) = parts.get(di).and_then(|s| s.parse::<u32>().ok()) {
                            s.psc_dp_sum[si] += d as u64;
                            s.psc_dp_count[si] += 1;
                            record_dp_sum += d as u64;
                        }
                    }
                    let alleles: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();
                    if alleles.iter().any(|a| *a == "." || a.is_empty()) {
                        s.psc_n_miss[si] += 1;
                        continue;
                    }
                    if is_snp { s.psc_n_snps[si] += 1; }
                    if is_indel { s.psc_n_indels[si] += 1; }
                    let nums: Vec<u32> = alleles.iter().filter_map(|a| a.parse().ok()).collect();
                    let any_nonref = nums.iter().any(|n| *n > 0);
                    let all_ref = nums.iter().all(|n| *n == 0);
                    let all_same = nums.iter().all(|n| *n == nums[0]);
                    if all_ref { s.psc_n_ref_hom[si] += 1; }
                    else if all_same { s.psc_n_hom[si] += 1; } else { s.psc_n_het[si] += 1; }
                    if any_nonref {
                        let alt_n: u32 = nums.iter().filter(|n| **n > 0).count() as u32;
                        alt_alleles_total += alt_n;
                        alt_carrier.push(si);
                    }
                    if is_indel && any_nonref {
                        s.psi_n_indel_alts[si] += 1;
                        if !all_same { s.psi_n_indel_hets[si] += 1; }
                        s.psi_n_indels[si] += 1;
                        for a in &alts {
                            let d = a.len() as i32 - refa.len() as i32;
                            if d % 3 == 0 { s.psi_n_in_frame[si] += 1; } else { s.psi_n_out_frame[si] += 1; }
                        }
                    }
                    n_ref_alleles += nums.iter().filter(|n| **n == 0).count() as u32;
                    n_alt_alleles += nums.iter().filter(|n| **n > 0).count() as u32;
                    gts_collected.push((si, nums));
                }
            }
            if alt_alleles_total == 1 && alt_carrier.len() == 1 {
                s.psc_n_singletons[alt_carrier[0]] += 1;
                s.n_singleton += 1;
            }
            if had_snp {
                let n = (n_ref_alleles + n_alt_alleles) as f64;
                if n > 0.0 {
                    let p = n_ref_alleles as f64 / n;
                    let bin = ((p * 100.0).round() as u32).min(100);
                    let entry = s.hwe_bins.entry(bin).or_default();
                    entry.0 += 1;
                    let mut chi2 = 0.0;
                    let mut n_hom_ref = 0u32; let mut n_het = 0u32; let mut n_hom_alt = 0u32;
                    for (_, nums) in &gts_collected {
                        if nums.len() != 2 { continue; }
                        match (nums[0] == 0, nums[1] == 0) {
                            (true, true) => n_hom_ref += 1,
                            (false, false) => n_hom_alt += 1,
                            _ => n_het += 1,
                        }
                    }
                    let total = (n_hom_ref + n_het + n_hom_alt) as f64;
                    if total > 0.0 {
                        let q = 1.0 - p;
                        let exp_hom_ref = p * p * total;
                        let exp_het = 2.0 * p * q * total;
                        let exp_hom_alt = q * q * total;
                        for (obs, exp) in [(n_hom_ref as f64, exp_hom_ref), (n_het as f64, exp_het), (n_hom_alt as f64, exp_hom_alt)] {
                            if exp > 0.0 { chi2 += (obs - exp).powi(2) / exp; }
                        }
                        if chi2 < 3.84 { entry.1 += 1; }
                    }
                }
            }
            if !s.samples.is_empty() {
                let avg_dp = record_dp_sum / s.samples.len() as u64;
                let bin = (avg_dp.min(500)) as u32;
                *s.dp_bins.entry(bin).or_default() += 1;
            }
        }
    }
    Ok(s)
}

fn extract_samples(h: &[String]) -> Vec<String> {
    for line in h {
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 9 { return cols[9..].iter().map(|s| s.to_string()).collect(); }
        }
    }
    Vec::new()
}

fn print_stats(files: &[&Path], sets: &[Stats], _args: &StatsArgs) -> Result<()> {
    use std::io::Write;
    let mut out = std::io::BufWriter::with_capacity(64 * 1024, std::io::stdout());

    writeln!(out, "# This file was produced by kira_bt stats (compatible with bcftools stats output)")?;
    writeln!(out, "# Definition of sets:")?;
    for (i, f) in files.iter().enumerate() {
        writeln!(out, "# ID\t{}\t{}", i, f.display())?;
    }

    writeln!(out, "# SN, Summary numbers:")?;
    writeln!(out, "# SN\t[2]id\t[3]key\t[4]value")?;
    for (i, s) in sets.iter().enumerate() {
        writeln!(out, "SN\t{}\tnumber of samples:\t{}", i, s.samples.len())?;
        writeln!(out, "SN\t{}\tnumber of records:\t{}", i, s.n_records)?;
        writeln!(out, "SN\t{}\tnumber of no-ALTs:\t{}", i, s.n_no_alts)?;
        writeln!(out, "SN\t{}\tnumber of SNPs:\t{}", i, s.n_snps)?;
        writeln!(out, "SN\t{}\tnumber of MNPs:\t{}", i, s.n_mnps)?;
        writeln!(out, "SN\t{}\tnumber of indels:\t{}", i, s.n_indels)?;
        writeln!(out, "SN\t{}\tnumber of others:\t{}", i, s.n_other)?;
        writeln!(out, "SN\t{}\tnumber of multiallelic sites:\t{}", i, s.n_multi)?;
        writeln!(out, "SN\t{}\tnumber of multiallelic SNP sites:\t{}", i, s.n_multi_snps)?;
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
            let c = s.sub_counts.get(&(*r, *a)).copied().unwrap_or(0);
            writeln!(out, "ST\t{}\t{}>{}\t{}", i, *r as char, *a as char, c)?;
        }
    }

    writeln!(out, "# IDD, InDel distribution:")?;
    writeln!(out, "# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF")?;
    for (i, s) in sets.iter().enumerate() {
        for (d, n) in &s.indel_dist {
            writeln!(out, "IDD\t{}\t{}\t{}\t0\t0", i, d, n)?;
        }
    }

    writeln!(out, "# QUAL, Stats by quality:")?;
    writeln!(out, "# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels")?;
    for (i, s) in sets.iter().enumerate() {
        for (q, n) in &s.qual_dist {
            writeln!(out, "QUAL\t{}\t{}\t{}\t0\t0\t0", i, q, n)?;
        }
    }

    writeln!(out, "# AF, Stats by non-reference allele frequency:")?;
    writeln!(out, "# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")?;
    for (i, s) in sets.iter().enumerate() {
        for (bin, n) in &s.af_dist {
            let af = *bin as f64 / 1000.0;
            writeln!(out, "AF\t{}\t{:.6}\t{}\t0\t0\t0\t0\t0\t0", i, af, n)?;
        }
    }

    writeln!(out, "# SiS, Singleton stats:")?;
    writeln!(out, "# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")?;
    for (i, s) in sets.iter().enumerate() {
        writeln!(out, "SiS\t{}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            i, s.sis_n_snps, s.sis_n_ts, s.sis_n_tv, s.sis_n_indels,
            s.sis_n_repeat_consistent, s.sis_n_repeat_inconsistent, s.sis_n_not_applicable)?;
    }

    writeln!(out, "# HWE, hardy weinberg equilibrium (PLINK definition):")?;
    writeln!(out, "# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25% to 75% percentile")?;
    for (i, s) in sets.iter().enumerate() {
        for (bin, (obs, in_eq)) in &s.hwe_bins {
            let af = *bin as f64 / 100.0;
            writeln!(out, "HWE\t{}\t{:.6}\t{}\t{}", i, af, obs, in_eq)?;
        }
    }

    writeln!(out, "# DP, Depth distribution")?;
    writeln!(out, "# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)")?;
    for (i, s) in sets.iter().enumerate() {
        let total: u64 = s.dp_bins.values().sum();
        for (bin, n) in &s.dp_bins {
            let pct = if total > 0 { 100.0 * (*n as f64) / (total as f64) } else { 0.0 };
            writeln!(out, "DP\t{}\t{}\t0\t0.000000\t{}\t{:.6}", i, bin, n, pct)?;
        }
    }

    writeln!(out, "# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.")?;
    writeln!(out, "# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing")?;
    for (i, s) in sets.iter().enumerate() {
        for (j, name) in s.samples.iter().enumerate() {
            let avg_dp = if s.psc_dp_count[j] > 0 { s.psc_dp_sum[j] as f64 / s.psc_dp_count[j] as f64 } else { 0.0 };
            writeln!(out, "PSC\t{}\t{}\t{}\t{}\t{}\t0\t0\t{}\t{:.1}\t{}\t0\t0\t{}",
                i, name, s.psc_n_ref_hom[j], s.psc_n_hom[j], s.psc_n_het[j],
                s.psc_n_indels[j], avg_dp, s.psc_n_singletons[j], s.psc_n_miss[j])?;
        }
    }

    writeln!(out, "# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del are counted twice, in both nInsHets and nDelHets.")?;
    writeln!(out, "# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms")?;
    for (i, s) in sets.iter().enumerate() {
        for (j, name) in s.samples.iter().enumerate() {
            let denom = (s.psi_n_in_frame[j] + s.psi_n_out_frame[j]) as f64;
            let ratio = if denom > 0.0 { s.psi_n_out_frame[j] as f64 / denom } else { 0.0 };
            writeln!(out, "PSI\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t0\t{}\t0",
                i, name, s.psi_n_in_frame[j], s.psi_n_out_frame[j],
                s.psi_n_not_applicable[j], ratio, s.psi_n_indel_hets[j], s.psi_n_indel_alts[j])?;
        }
    }

    out.flush()?;
    Ok(())
}
