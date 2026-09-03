//! `bcftools merge` port: records at a position are grouped by the `-m`
//! rules (SNVs first, then indels), INFO follows the bcftools ordering
//! (plain tags, rule tags, A/R/G tags, AN/AC) and AN/AC are recomputed from
//! the merged genotypes.

use crate::annotate::postproc::{RegionFilter, version_header_line};
use crate::cli::args::MergeArgs;
use crate::vcf::alleles::{gt_alleles, join_info, remap_samples, remap_value, split_info};
use crate::vcf::header::{ContigDict, FieldNumber, FieldType, HeaderInfo, extract_samples};
use crate::vcf::sink::{OutputKind, parse_output_type};
use crate::vcf::variant_type::{VT_INDEL, VT_REF, VT_SNP, allele_type};
#[cfg(test)]
use crate::vcf::variant_type::record_type;
use crate::vcf::{UnifiedVcfReader, VcfSink};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

pub fn cmd_merge(args: MergeArgs) -> Result<()> {
    let mut inputs: Vec<PathBuf> = args.inputs.clone();
    if let Some(fl) = &args.file_list {
        for line in BufReader::new(File::open(fl)?).lines() {
            let l = line?;
            let t = l.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            inputs.push(PathBuf::from(t));
        }
    }
    if inputs.is_empty() { bail!("merge: no inputs"); }

    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);
    let do_gvcf = args.gvcf.is_some();
    let merge_mode = MergeMode::parse(&args.merge)?;
    let filter_logic = match args.filter_logic.as_str() {
        "+" => FilterLogic::Union,
        "x" => FilterLogic::Exclude,
        o => bail!("--filter-logic: expected + or x, got {o:?}"),
    };
    let apply_filters: Option<Vec<String>> = args.apply_filters.as_deref().map(|s| s.split(',').map(|t| t.trim().to_string()).collect());
    let region_filter = if let Some(s) = &args.regions {
        Some(RegionFilter::from_cli(s)?)
    } else if let Some(p) = &args.regions_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };

    let mut readers: Vec<Source> = Vec::with_capacity(inputs.len());
    let mut all_samples: Vec<String> = Vec::new();
    let mut all_meta_lines: Vec<String> = Vec::new();
    let mut seen_meta: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut contigs = ContigDict::new();

    for (fi, p) in inputs.iter().enumerate() {
        let r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let headers = r.header()?;
        let samples = extract_samples(&headers);
        let mut local_idx_in_global: Vec<usize> = Vec::with_capacity(samples.len());
        for s in &samples {
            let mut name = s.clone();
            if all_samples.contains(&name) {
                if !args.force_samples {
                    bail!("Error: Duplicate sample names ({s}), use --force-samples to proceed anyway.");
                }
                // bcftools prefixes the file number until the name is unique.
                while all_samples.contains(&name) {
                    name = format!("{}:{}", fi + 1, name);
                }
            }
            local_idx_in_global.push(all_samples.len());
            all_samples.push(name);
        }
        for h in &headers {
            if h.starts_with("##") && seen_meta.insert(h.clone()) {
                all_meta_lines.push(h.clone());
            }
            if let Some((id, len)) = crate::vcf::header::parse_contig_line(h) {
                contigs.insert_with_length(&id, len);
            }
        }
        let mut src = Source {
            reader: r,
            next: None,
            local_to_global: local_idx_in_global,
            region_filter: region_filter.clone(),
            regions_overlap: args.regions_overlap,
            apply_filters: apply_filters.clone(),
            gvcf_block: None,
        };
        src.advance(&mut contigs)?;
        readers.push(src);
    }

    if let Some(p) = &args.use_header {
        let text = std::fs::read_to_string(p).with_context(|| format!("--use-header {}", p.display()))?;
        all_meta_lines = text.lines().filter(|l| l.starts_with("##")).map(|l| l.to_string()).collect();
    }

    let merged_hdr = HeaderInfo::parse(&all_meta_lines);
    let rules = InfoRules::build(args.info_rules.as_deref(), &merged_hdr, all_samples.is_empty(), do_gvcf)?;
    // `join` turns the tag into a variable-length vector in the header.
    for (k, r) in &rules.rules {
        if matches!(r, InfoRule::Join) {
            for l in all_meta_lines.iter_mut() {
                if l.starts_with(&format!("##INFO=<ID={k},")) && !l.contains("Number=.") {
                    if let Some(s) = l.find("Number=") {
                        let e = l[s..].find(',').map(|x| s + x).unwrap_or(l.len());
                        l.replace_range(s..e, "Number=.");
                    }
                }
            }
        }
    }
    let mut full_headers: Vec<String> = all_meta_lines.clone();
    if !args.no_version { full_headers.push(version_header_line()); }
    let mut chrom_line = String::from("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !all_samples.is_empty() {
        chrom_line.push_str("\tFORMAT");
        for s in &all_samples { chrom_line.push('\t'); chrom_line.push_str(s); }
    }
    full_headers.push(chrom_line);

    let mut sink = VcfSink::open(args.output.as_deref(), kind, &full_headers)?;
    sink.write_header(&full_headers)?;
    if args.print_header { sink.finish()?; return Ok(()); }

    let ctx = MergeCtx {
        hdr: &merged_hdr,
        all_samples: &all_samples,
        missing_to_ref: args.missing_to_ref,
        rules: &rules,
        filter_logic,
        do_gvcf,
    };

    loop {
        // Next site: smallest (contig rank, pos) across readers.
        let Some((key_rank, key_pos)) = readers.iter().filter_map(|s| s.next.as_ref().map(|r| (r.rank, r.pos))).min() else { break };

        // Every record of every reader at this position.
        let mut bufs: Vec<Vec<Record>> = Vec::with_capacity(readers.len());
        for r in readers.iter_mut() {
            let mut v = Vec::new();
            while r.next.as_ref().is_some_and(|n| n.rank == key_rank && n.pos == key_pos) {
                v.push(r.next.take().unwrap());
                r.advance(&mut contigs)?;
            }
            bufs.push(v);
        }
        let mut done: Vec<Vec<bool>> = bufs.iter().map(|v| vec![false; v.len()]).collect();
        loop {
            let staged = stage_round(&bufs, &done, merge_mode);
            if staged.iter().all(|s| s.is_none()) { break; }
            let group: Vec<(usize, Record)> = staged
                .iter()
                .enumerate()
                .filter_map(|(i, s)| s.map(|j| (i, bufs[i][j].clone())))
                .collect();
            for (i, s) in staged.iter().enumerate() {
                if let Some(j) = s { done[i][*j] = true; }
            }
            let merged = merge_group(&group, &readers, &ctx, key_rank, key_pos);
            sink.write_line(&merged)?;
        }
    }
    sink.finish()?;
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FilterLogic { Union, Exclude }

#[derive(Clone, Debug)]
enum InfoRule { Sum, Avg, Min, Max, Join, First }

fn parse_info_rules(spec: Option<&str>) -> Result<HashMap<String, InfoRule>> {
    let mut m: HashMap<String, InfoRule> = HashMap::new();
    // bcftools defaults.
    m.insert("DP".into(), InfoRule::Sum);
    m.insert("DP4".into(), InfoRule::Sum);
    let Some(s) = spec else { return Ok(m); };
    if s == "-" {
        m.clear();
        return Ok(m);
    }
    for tok in s.split(',') {
        let tok = tok.trim();
        if tok.is_empty() { continue; }
        let Some((key, rule)) = tok.split_once(':') else { bail!("--info-rules: expected TAG:METHOD, got {tok:?}") };
        let r = match rule.trim() {
            "sum" => InfoRule::Sum,
            "avg" => InfoRule::Avg,
            "min" => InfoRule::Min,
            "max" => InfoRule::Max,
            "join" => InfoRule::Join,
            "first" => InfoRule::First,
            o => bail!("--info-rules: unknown method {o:?} for {key}"),
        };
        m.insert(key.trim().to_string(), r);
    }
    Ok(m)
}

/// Active INFO rules: the user's list, or the bcftools defaults (`DP:sum,DP4:sum`,
/// plus `AN:sum,AC:sum` for sites-only output and the gVCF tags).
struct InfoRules {
    rules: BTreeMap<String, InfoRule>,
    /// AN/AC were given explicitly and are not recomputed from genotypes.
    keep_ac_an: bool,
}

impl InfoRules {
    fn build(spec: Option<&str>, hdr: &HeaderInfo, no_samples: bool, do_gvcf: bool) -> Result<Self> {
        let mut rules: BTreeMap<String, InfoRule> = BTreeMap::new();
        let mut keep_ac_an = false;
        match spec {
            Some("-") => {}
            Some(s) => {
                for (k, r) in parse_info_rules(Some(s))? {
                    if (k == "AC" || k == "AN") && s.split(',').any(|t| t.trim().starts_with(&format!("{k}:"))) {
                        keep_ac_an = true;
                    }
                    if s.split(',').any(|t| t.trim().starts_with(&format!("{k}:"))) {
                        rules.insert(k, r);
                    }
                }
            }
            None => {
                if hdr.info.contains_key("DP") { rules.insert("DP".into(), InfoRule::Sum); }
                if hdr.info.contains_key("DP4") { rules.insert("DP4".into(), InfoRule::Sum); }
                if do_gvcf {
                    for (k, r) in [("QS", InfoRule::Sum), ("MIN_DP", InfoRule::Min), ("MinDP", InfoRule::Min), ("I16", InfoRule::Sum), ("IDV", InfoRule::Max), ("IMF", InfoRule::Max)] {
                        if hdr.info.contains_key(k) { rules.insert(k.into(), r); }
                    }
                }
                if no_samples {
                    if hdr.info.contains_key("AN") { rules.insert("AN".into(), InfoRule::Sum); }
                    if hdr.info.contains_key("AC") { rules.insert("AC".into(), InfoRule::Sum); }
                }
            }
        }
        Ok(Self { rules, keep_ac_an })
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum MergeMode { Snps, Indels, Both, All, None_, Id, SnpInsDel }

impl MergeMode {
    fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "snps" => Self::Snps,
            "indels" => Self::Indels,
            "both" => Self::Both,
            "all" | "any" => Self::All,
            "none" => Self::None_,
            "id" => Self::Id,
            "snp-ins-del" => Self::SnpInsDel,
            o => bail!("--merge: unknown {o:?}"),
        })
    }

    fn merges_snps(self) -> bool { matches!(self, Self::Snps | Self::Both | Self::SnpInsDel) }
    fn merges_indels(self) -> bool { matches!(self, Self::Indels | Self::Both) }
}

// Variant-type bits of the bcftools merge buffer (htslib types shifted by one).
const M_REF: u32 = 1;
const M_SNP: u32 = 1 << 1;
const M_MNP: u32 = 1 << 2;
const M_INDEL: u32 = 1 << 3;
const M_OTHER: u32 = 1 << 4;
const M_INS: u32 = 1 << 7;
const M_DEL: u32 = 1 << 8;

/// Type mask of one ALT allele.
fn allele_mask(refa: &str, alt: &str) -> u32 {
    let t = allele_type(refa, alt);
    let mut m = 0;
    if t.ty & VT_SNP != 0 { m |= M_SNP; }
    if t.ty & crate::vcf::variant_type::VT_MNP != 0 { m |= M_MNP; }
    if t.ty & VT_INDEL != 0 {
        m |= M_INDEL;
        m |= if t.n > 0 { M_INS } else { M_DEL };
    }
    if t.ty & crate::vcf::variant_type::VT_OTHER != 0 { m |= M_OTHER; }
    if t.ty & crate::vcf::variant_type::VT_BND != 0 { m |= M_OTHER; }
    m
}

/// Type mask of a record (`ref_mask` when it has no variant alleles).
fn record_mask(r: &Record, mode: MergeMode) -> u32 {
    let mut m = 0;
    if r.alt != "." && !r.alt.is_empty() {
        for a in r.alt.split(',') {
            m |= allele_mask(&r.refa, a);
        }
    }
    if mode == MergeMode::SnpInsDel { m &= !M_INDEL; }
    let mut m = if m == 0 { M_REF } else { m };
    if is_gvcf_block(r) { m |= M_REF; }
    m
}

/// Pad `alt` of a record with REF `refa` to the longer reference `long_ref`.
fn pad_alt(alt: &str, refa: &str, long_ref: &str) -> String {
    if alt.starts_with('<') || alt == "*" || alt.starts_with('.') || long_ref.len() <= refa.len() {
        return alt.to_uppercase();
    }
    if long_ref.as_bytes()[..refa.len()].eq_ignore_ascii_case(refa.as_bytes()) {
        format!("{}{}", alt.to_uppercase(), long_ref[refa.len()..].to_uppercase())
    } else {
        alt.to_uppercase()
    }
}

fn refs_compatible(a: &str, b: &str) -> bool {
    let n = a.len().min(b.len());
    a.as_bytes()[..n].eq_ignore_ascii_case(&b.as_bytes()[..n])
}

/// bcftools `types_compatible`: may `rec` join the records selected so far?
fn types_compatible(mode: MergeMode, sel_types: u32, sel_ref: &str, sel_alts: &[String], rec: &Record, rec_types: u32) -> bool {
    if mode == MergeMode::All { return true; }
    if sel_types & M_REF != 0 && sel_types & !M_REF == 0 { return true; }
    if rec_types & M_REF != 0 && rec_types & !M_REF == 0 { return true; }
    if mode != MergeMode::None_ {
        if mode.merges_snps() && rec_types & M_SNP != 0 && sel_types & M_SNP != 0 { return true; }
        if mode.merges_indels() && rec_types & M_INDEL != 0 && sel_types & M_INDEL != 0 { return true; }
        if mode == MergeMode::SnpInsDel {
            if rec_types & M_INS != 0 && sel_types & M_INS != 0 { return true; }
            if rec_types & M_DEL != 0 && sel_types & M_DEL != 0 { return true; }
        }
    }
    // Exact matching: same variant types, compatible REF and a shared allele.
    let (mut x, mut y) = (sel_types >> 1, rec_types >> 1);
    while x != 0 && y != 0 { x >>= 1; y >>= 1; }
    if x != 0 || y != 0 { return false; }
    if !refs_compatible(sel_ref, &rec.refa) { return false; }
    let long_ref = if sel_ref.len() >= rec.refa.len() { sel_ref.to_string() } else { rec.refa.to_uppercase() };
    let sel_padded: Vec<String> = sel_alts.iter().map(|a| pad_alt(a, sel_ref, &long_ref)).collect();
    if rec.alt == "." || rec.alt.is_empty() { return false; }
    for a in rec.alt.split(',') {
        if allele_type(&rec.refa, a).ty == VT_REF { continue; }
        let p = pad_alt(a, &rec.refa, &long_ref);
        if sel_padded.iter().any(|s| *s == p) { return true; }
    }
    false
}

/// One merge round at a position (`can_merge` + `stage_line`): at most one
/// record per reader is staged; the rest wait for the next round.
fn stage_round(bufs: &[Vec<Record>], done: &[Vec<bool>], mode: MergeMode) -> Vec<Option<usize>> {
    let n = bufs.len();
    let mut staged: Vec<Option<usize>> = vec![None; n];
    let types: Vec<Vec<u32>> = bufs.iter().map(|v| v.iter().map(|r| record_mask(r, mode)).collect()).collect();
    let mut var_all = 0u32;
    let mut reader_types = vec![0u32; n];
    let mut ntodo = 0;
    let mut first_id: Option<&str> = None;
    for i in 0..n {
        for j in 0..bufs[i].len() {
            if done[i][j] { continue; }
            ntodo += 1;
            if mode == MergeMode::Id && first_id.is_none() { first_id = Some(&bufs[i][j].id); continue; }
            var_all |= types[i][j];
            reader_types[i] |= types[i][j];
        }
    }
    if ntodo == 0 { return staged; }

    // Candidate records compatible with the growing selection.
    let mut sel_types = 0u32;
    let mut sel_ref = String::new();
    let mut sel_alts: Vec<String> = Vec::new();
    let mut cand: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..bufs[i].len() {
            if done[i][j] { continue; }
            let rec = &bufs[i][j];
            if mode == MergeMode::Id {
                if first_id != Some(rec.id.as_str()) { continue; }
            } else if sel_types != 0 && !types_compatible(mode, sel_types, &sel_ref, &sel_alts, rec, types[i][j]) {
                continue;
            } else if (mode.merges_snps() || mode == MergeMode::None_) && var_all & M_SNP != 0 && types[i][j] & (M_SNP | M_REF) == 0 {
                // SNVs go first when present.
                continue;
            }
            sel_types |= types[i][j];
            // Extend the selected allele set (padding to the longest REF).
            let up = rec.refa.to_uppercase();
            if up.len() > sel_ref.len() {
                let old = sel_ref.clone();
                sel_ref = up.clone();
                if !old.is_empty() {
                    sel_alts = sel_alts.iter().map(|a| pad_alt(a, &old, &sel_ref)).collect();
                }
            }
            if rec.alt != "." && !rec.alt.is_empty() {
                for a in rec.alt.split(',') {
                    let p = pad_alt(a, &rec.refa, &sel_ref);
                    if !sel_alts.contains(&p) { sel_alts.push(p); }
                }
            }
            cand[i].push(j);
        }
    }

    // Allele counts across the candidates; the most frequent non-ref allele leads.
    let mut cnt: Vec<u32> = vec![0; sel_alts.len()];
    let mut al_types: Vec<u32> = vec![0; sel_alts.len()];
    let rec_alleles = |i: usize, j: usize| -> Vec<usize> {
        let rec = &bufs[i][j];
        let mut out = Vec::new();
        if rec.alt != "." && !rec.alt.is_empty() {
            for a in rec.alt.split(',') {
                let p = pad_alt(a, &rec.refa, &sel_ref);
                if let Some(k) = sel_alts.iter().position(|s| *s == p) { out.push(k); }
            }
        }
        out
    };
    for i in 0..n {
        for &j in &cand[i] {
            for k in rec_alleles(i, j) {
                cnt[k] += 1;
                let t = allele_mask(&sel_ref, &sel_alts[k]);
                al_types[k] = if t == 0 { M_REF } else if mode == MergeMode::SnpInsDel { t & !M_INDEL } else { t };
            }
        }
    }
    let mut icnt: Option<usize> = None;
    for k in 0..sel_alts.len() {
        if al_types[k] & M_REF != 0 { continue; }
        if icnt.is_none_or(|c| cnt[c] < cnt[k]) { icnt = Some(k); }
    }
    let selected_type = icnt.map(|k| al_types[k]).unwrap_or(M_REF);

    for i in 0..n {
        let mut cur: Option<usize> = None;
        for &j in &cand[i] {
            if mode == MergeMode::Id { cur = Some(j); break; }
            let t = types[i][j];
            if selected_type & reader_types[i] != 0 && selected_type & t == 0 { continue; }
            if selected_type & reader_types[i] == 0 && t & M_REF != 0 { cur = Some(j); break; }
            if selected_type == M_REF { cur = Some(j); break; }
            if let Some(k) = icnt {
                if rec_alleles(i, j).contains(&k) { cur = Some(j); break; }
            }
        }
        if cur.is_none() && mode != MergeMode::None_ && mode != MergeMode::Id {
            for &j in &cand[i] {
                let t = types[i][j];
                if mode == MergeMode::All { cur = Some(j); break; }
                if var_all & M_SNP != 0 && t & M_SNP != 0 && mode.merges_snps() { cur = Some(j); break; }
                if var_all & M_INDEL != 0 && t & M_INDEL != 0 && mode.merges_indels() { cur = Some(j); break; }
                if var_all & M_INS != 0 && t & M_INS != 0 && mode == MergeMode::SnpInsDel { cur = Some(j); break; }
                if var_all & M_DEL != 0 && t & M_DEL != 0 && mode == MergeMode::SnpInsDel { cur = Some(j); break; }
                if t & M_REF != 0 {
                    if var_all & M_SNP != 0 && mode.merges_snps() { cur = Some(j); break; }
                    if var_all & M_INDEL != 0 && mode.merges_indels() { cur = Some(j); break; }
                    if var_all & (M_INS | M_DEL) != 0 && mode == MergeMode::SnpInsDel { cur = Some(j); break; }
                    if var_all & M_REF != 0 { cur = Some(j); break; }
                } else if var_all & M_REF != 0 {
                    if t & M_SNP != 0 && mode.merges_snps() { cur = Some(j); break; }
                    if t & M_INDEL != 0 && mode.merges_indels() { cur = Some(j); break; }
                    if t & (M_INS | M_DEL) != 0 && mode == MergeMode::SnpInsDel { cur = Some(j); break; }
                }
            }
        }
        staged[i] = cur;
    }
    if staged.iter().all(|s| s.is_none()) {
        // Nothing staged although work remains: take the first unprocessed record.
        for i in 0..n {
            if let Some(j) = (0..bufs[i].len()).find(|&j| !done[i][j]) {
                staged[i] = Some(j);
                break;
            }
        }
    }
    staged
}

struct Source {
    reader: UnifiedVcfReader,
    next: Option<Record>,
    local_to_global: Vec<usize>,
    region_filter: Option<RegionFilter>,
    regions_overlap: u8,
    apply_filters: Option<Vec<String>>,
    /// Open gVCF reference block (`<NON_REF>`/`<*>` with END) that later sites may fall into.
    gvcf_block: Option<Record>,
}

impl Source {
    fn advance(&mut self, contigs: &mut ContigDict) -> Result<()> {
        loop {
            let Some(line) = self.reader.read_line()? else {
                self.next = None;
                return Ok(());
            };
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            if let Some(rf) = &self.region_filter {
                if !rf.line_passes_mode(&line, self.regions_overlap) { continue; }
            }
            let rec = parse_record(line, contigs)?;
            if is_gvcf_block(&rec) {
                self.gvcf_block = Some(rec.clone());
            }
            if let Some(af) = &self.apply_filters {
                let pass = if rec.filter == "." || rec.filter.is_empty() {
                    af.iter().any(|a| a == ".")
                } else {
                    rec.filter.split(';').any(|t| af.iter().any(|a| a == t))
                };
                if !pass { continue; }
            }
            self.next = Some(rec);
            return Ok(());
        }
    }
}

#[derive(Clone, Default)]
struct Record {
    chrom: String,
    rank: usize,
    pos: u32,
    id: String,
    refa: String,
    alt: String,
    qual: String,
    filter: String,
    info: String,
    format: String,
    samples: Vec<String>,
}

fn parse_record(line: String, contigs: &mut ContigDict) -> Result<Record> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 {
        bail!("merge: malformed record (fewer than 8 columns): {line:?}");
    }
    let pos: u32 = cols[1].parse().with_context(|| format!("merge: invalid POS {:?}", cols[1]))?;
    let format = if cols.len() > 8 { cols[8].to_string() } else { String::new() };
    let samples = if cols.len() > 9 { cols[9..].iter().map(|s| s.to_string()).collect() } else { Vec::new() };
    let chrom = cols.first().copied().unwrap_or("").to_string();
    let rank = contigs.insert(&chrom) as usize;
    Ok(Record {
        chrom,
        rank,
        pos,
        id: cols.get(2).copied().unwrap_or(".").to_string(),
        refa: cols.get(3).copied().unwrap_or(".").to_string(),
        alt: cols.get(4).copied().unwrap_or(".").to_string(),
        qual: cols.get(5).copied().unwrap_or(".").to_string(),
        filter: cols.get(6).copied().unwrap_or(".").to_string(),
        info: cols.get(7).copied().unwrap_or(".").to_string(),
        format,
        samples,
    })
}

fn is_gvcf_block(r: &Record) -> bool {
    (r.alt == "<NON_REF>" || r.alt == "<*>") && r.info.split(';').any(|kv| kv.starts_with("END="))
}

/// Pairwise compatibility of two records under `-m` (kept for the unit tests;
/// the merge itself uses `stage_round`).
#[cfg(test)]
fn matches_merge_mode(rec: &Record, anchor: &Record, mode: MergeMode) -> bool {
    let rt = record_type(&rec.refa, &rec.alt);
    let at = record_type(&anchor.refa, &anchor.alt);
    let same_ref = rec.refa.eq_ignore_ascii_case(&anchor.refa);
    match mode {
        MergeMode::All => true,
        MergeMode::None_ => same_ref && rec.alt.eq_ignore_ascii_case(&anchor.alt),
        MergeMode::Id => rec.id == anchor.id,
        MergeMode::Both => {
            let snp = rt & VT_SNP != 0 && at & VT_SNP != 0;
            let indel = rt & VT_INDEL != 0 && at & VT_INDEL != 0;
            snp || indel || (rt == at && same_ref)
        }
        MergeMode::Snps => (rt & VT_SNP != 0 && at & VT_SNP != 0) || (rt & VT_INDEL == 0 && at & VT_INDEL == 0 && same_ref && rec.alt == anchor.alt),
        MergeMode::Indels => (rt & VT_INDEL != 0 && at & VT_INDEL != 0) || (rt & VT_SNP == 0 && at & VT_SNP == 0 && same_ref && rec.alt == anchor.alt),
        MergeMode::SnpInsDel => {
            let cls = |r: &Record| -> u8 {
                let t = record_type(&r.refa, &r.alt);
                if t & VT_SNP != 0 { 1 } else if t & VT_INDEL != 0 { if r.alt.len() > r.refa.len() { 2 } else { 3 } } else { 4 }
            };
            cls(rec) == cls(anchor)
        }
    }
}

struct MergeCtx<'a> {
    hdr: &'a HeaderInfo,
    all_samples: &'a [String],
    missing_to_ref: bool,
    rules: &'a InfoRules,
    filter_logic: FilterLogic,
    do_gvcf: bool,
}

/// htslib drops the phase of a missing allele (`0|.` is written `0/.`).
fn normalize_gt(gt: &str) -> String {
    if !gt.contains('|') { return gt.to_string(); }
    let mut out = String::with_capacity(gt.len());
    let mut cur = String::new();
    let mut seps: Vec<char> = Vec::new();
    let mut alleles: Vec<String> = Vec::new();
    for c in gt.chars() {
        if c == '|' || c == '/' {
            seps.push(c);
            alleles.push(std::mem::take(&mut cur));
        } else {
            cur.push(c);
        }
    }
    alleles.push(cur);
    for (i, a) in alleles.iter().enumerate() {
        if i > 0 {
            out.push(if a == "." { '/' } else { seps[i - 1] });
        }
        out.push_str(a);
    }
    out
}

/// C `%g` with six significant digits (htslib float output).
fn fmt_g(v: f64) -> String {
    if v == 0.0 { return "0".into(); }
    if v.is_nan() { return "nan".into(); }
    if v.is_infinite() { return if v > 0.0 { "inf".into() } else { "-inf".into() }; }
    let sci = format!("{:.5e}", v);
    let (mant, exp) = sci.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap_or(0);
    if !(-4..6).contains(&exp) {
        let m = mant.trim_end_matches('0').trim_end_matches('.');
        return format!("{}e{}{:02}", m, if exp < 0 { '-' } else { '+' }, exp.abs());
    }
    let decimals = (5 - exp).max(0) as usize;
    let s = format!("{:.*}", decimals, v);
    if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_string() } else { s }
}

fn fmt_float_list(v: &str) -> String {
    v.split(',')
        .map(|x| match x.parse::<f64>() {
            Ok(f) if x != "." => fmt_g(f),
            _ => x.to_string(),
        })
        .collect::<Vec<_>>()
        .join(",")
}

/// Merge the records of one site. `readers` carries the per-file sample maps
/// and any open gVCF block used to fill absent samples.
fn merge_group(group: &[(usize, Record)], readers: &[Source], ctx: &MergeCtx<'_>, rank: usize, pos: u32) -> String {
    let anchor = &group[0].1;

    // Union of alleles with REF padding (shorter REFs are prefixes of the longest).
    let longest = group.iter().map(|(_, r)| r.refa.len()).max().unwrap_or(0);
    let long_ref = group.iter().find(|(_, r)| r.refa.len() == longest).map(|(_, r)| r.refa.clone()).unwrap_or_default();
    let mut all_alts: Vec<String> = Vec::new();
    let mut maps: Vec<Vec<Option<usize>>> = Vec::with_capacity(group.len());
    for (_, r) in group {
        let tail = if long_ref.len() > r.refa.len() && long_ref.as_bytes()[..r.refa.len()].eq_ignore_ascii_case(r.refa.as_bytes()) {
            &long_ref[r.refa.len()..]
        } else {
            ""
        };
        let mut map = vec![Some(0)];
        if r.alt != "." && !r.alt.is_empty() {
            for a in r.alt.split(',') {
                let padded = if a.starts_with('<') || a == "*" || a.starts_with('.') || tail.is_empty() { a.to_string() } else { format!("{a}{tail}") };
                let idx = match all_alts.iter().position(|x| x.eq_ignore_ascii_case(&padded)) {
                    Some(i) => i,
                    None => {
                        all_alts.push(padded);
                        all_alts.len() - 1
                    }
                };
                map.push(Some(idx + 1));
            }
        }
        maps.push(map);
    }
    let n_new = all_alts.len() + 1;
    let long_ref = normalize_alleles(long_ref, &mut all_alts);

    let merged_id = {
        let mut ids: Vec<&str> = Vec::new();
        for (_, r) in group {
            for id in r.id.split(';') {
                if id != "." && !ids.contains(&id) { ids.push(id); }
            }
        }
        if ids.is_empty() { ".".to_string() } else { ids.join(";") }
    };
    let merged_qual = group
        .iter()
        .filter_map(|(_, r)| r.qual.parse::<f64>().ok())
        .fold(None, |m: Option<f64>, q| Some(m.map_or(q, |v| v.max(q))))
        .map(fmt_g)
        .unwrap_or_else(|| ".".to_string());
    let merged_filter = merge_filters(group.iter().map(|(_, r)| r.filter.as_str()), ctx.filter_logic);

    let merged_format = find_common_format(group);
    let format_keys: Vec<&str> = merged_format.split(':').collect();
    let gt_idx = if merged_format.is_empty() { None } else { format_keys.iter().position(|k| *k == "GT") };
    let float_keys: Vec<bool> = format_keys.iter().map(|k| matches!(ctx.hdr.format_type(k), Some(FieldType::Float))).collect();

    // Per-sample columns, default missing.
    let missing_col = |keys: &[&str]| -> String {
        keys.iter().map(|k| if *k == "GT" { "./." } else { "." }).collect::<Vec<_>>().join(":")
    };
    let mut sample_cols: Vec<String> = vec![missing_col(&format_keys); ctx.all_samples.len()];
    let mut sample_filled: Vec<bool> = vec![false; ctx.all_samples.len()];

    for (gi, (src_i, r)) in group.iter().enumerate() {
        if r.samples.is_empty() { continue; }
        let n_old = maps[gi].len();
        let samples: Vec<&str> = r.samples.iter().map(String::as_str).collect();
        let remapped = remap_samples(&r.format, &samples, ctx.hdr, n_old, n_new, &maps[gi]);
        let local_keys: Vec<&str> = r.format.split(':').collect();
        for (li, sval) in remapped.iter().enumerate() {
            let Some(&global) = readers[*src_i].local_to_global.get(li) else { continue };
            if global >= sample_cols.len() { continue; }
            let local_vals: Vec<&str> = sval.split(':').collect();
            let out_vals: Vec<String> = format_keys
                .iter()
                .enumerate()
                .map(|(ki, k)| {
                    let v = local_keys
                        .iter()
                        .position(|lk| lk == k)
                        .and_then(|i| local_vals.get(i).copied())
                        .unwrap_or(if *k == "GT" { "./." } else { "." });
                    if *k == "GT" { normalize_gt(v) } else if float_keys[ki] { fmt_float_list(v) } else { v.to_string() }
                })
                .collect();
            sample_cols[global] = out_vals.join(":");
            sample_filled[global] = true;
        }
    }

    // Samples absent at this site: gVCF reference blocks or --missing-to-ref.
    for (src_i, src) in readers.iter().enumerate() {
        let in_group = group.iter().any(|(g, _)| *g == src_i);
        if in_group { continue; }
        let block = if ctx.do_gvcf {
            src.gvcf_block.as_ref().filter(|b| b.rank == rank && b.pos <= pos && block_end(b) >= pos)
        } else {
            None
        };
        for (li, &global) in src.local_to_global.iter().enumerate() {
            if global >= sample_cols.len() || sample_filled[global] { continue; }
            if let Some(b) = block {
                let bkeys: Vec<&str> = b.format.split(':').collect();
                let bvals: Vec<&str> = b.samples.get(li).map(|s| s.split(':').collect()).unwrap_or_default();
                let out_vals: Vec<String> = format_keys
                    .iter()
                    .map(|k| {
                        if *k == "GT" {
                            let gt = bkeys.iter().position(|x| *x == "GT").and_then(|i| bvals.get(i).copied()).unwrap_or("./.");
                            // Any non-ref index inside a ref block is meaningless; keep ploidy.
                            gt_alleles(gt).iter().map(|a| if a.is_some() { "0" } else { "." }).collect::<Vec<_>>().join(if gt.contains('|') { "|" } else { "/" })
                        } else if let Some(i) = bkeys.iter().position(|x| x == k) {
                            let v = bvals.get(i).copied().unwrap_or(".");
                            match ctx.hdr.format_number(k) {
                                FieldNumber::A => vec!["."; n_new - 1].join(","),
                                FieldNumber::R => {
                                    let first = v.split(',').next().unwrap_or(".");
                                    let mut vv = vec![first]; vv.extend(std::iter::repeat_n(".", n_new - 1)); vv.join(",")
                                }
                                FieldNumber::G => {
                                    let first = v.split(',').next().unwrap_or(".");
                                    let n = crate::vcf::alleles::n_diploid_genotypes(n_new);
                                    let mut vv = vec![first]; vv.extend(std::iter::repeat_n(".", n - 1)); vv.join(",")
                                }
                                _ => v.to_string(),
                            }
                        } else {
                            ".".to_string()
                        }
                    })
                    .collect();
                sample_cols[global] = out_vals.join(":");
                sample_filled[global] = true;
            } else if ctx.missing_to_ref {
                if let Some(gi) = gt_idx {
                    let mut parts: Vec<&str> = sample_cols[global].split(':').map(|_| ".").collect();
                    parts[gi] = "0/0";
                    sample_cols[global] = parts.join(":");
                }
            }
        }
    }

    let merged_info = merge_info(group, &maps, n_new, ctx, gt_idx, &sample_cols);

    let mut out = String::new();
    out.push_str(&anchor.chrom); out.push('\t');
    out.push_str(&pos.to_string()); out.push('\t');
    out.push_str(&merged_id); out.push('\t');
    out.push_str(&long_ref); out.push('\t');
    out.push_str(&if all_alts.is_empty() { ".".to_string() } else { all_alts.join(",") }); out.push('\t');
    out.push_str(&merged_qual); out.push('\t');
    out.push_str(&merged_filter); out.push('\t');
    out.push_str(&merged_info);
    if !merged_format.is_empty() && !ctx.all_samples.is_empty() {
        out.push('\t'); out.push_str(&merged_format);
        for s in &sample_cols { out.push('\t'); out.push_str(s); }
    }
    out
}

/// bcftools `normalize_alleles`: drop trailing bases shared by every allele.
fn normalize_alleles(refa: String, alts: &mut [String]) -> String {
    if refa.len() < 2 { return refa; }
    let rb = refa.as_bytes();
    let mut i = 1usize;
    let mut done = false;
    while i < rb.len() {
        for a in alts.iter() {
            let ab = a.as_bytes();
            if i >= ab.len() { done = true; }
            if ab.get(ab.len().wrapping_sub(i)) != Some(&rb[rb.len() - i]) { done = true; break; }
        }
        if done { break; }
        i += 1;
    }
    if i <= 1 { return refa; }
    let cut = i - 1;
    for a in alts.iter_mut() {
        let n = a.len() - cut;
        a.truncate(n);
    }
    refa[..refa.len() - cut].to_string()
}

fn block_end(r: &Record) -> u32 {
    r.info
        .split(';')
        .find_map(|kv| kv.strip_prefix("END=").and_then(|v| v.parse::<u32>().ok()))
        .unwrap_or(r.pos)
}

fn merge_filters<'a, I: Iterator<Item = &'a str>>(filters: I, logic: FilterLogic) -> String {
    let mut set: Vec<&str> = Vec::new();
    let mut any = false;
    let mut any_pass = false;
    let mut all_pass = true;
    for f in filters {
        any = true;
        if f == "PASS" {
            any_pass = true;
            continue;
        }
        all_pass = false;
        if f == "." || f.is_empty() { continue; }
        for t in f.split(';') {
            if !set.contains(&t) { set.push(t); }
        }
    }
    if !any { return ".".into(); }
    match logic {
        FilterLogic::Union => {
            if all_pass { "PASS".into() } else if set.is_empty() { if any_pass { "PASS".into() } else { ".".into() } } else { set.join(";") }
        }
        FilterLogic::Exclude => {
            // PASS if any record passes (an unfiltered "." counts), else the union of failing filters.
            if any_pass || set.is_empty() { "PASS".into() } else { set.join(";") }
        }
    }
}

fn find_common_format(group: &[(usize, Record)]) -> String {
    let mut keys: Vec<String> = Vec::new();
    for (_, r) in group {
        if r.format.is_empty() { continue; }
        for k in r.format.split(':') {
            if !keys.iter().any(|x| x == k) { keys.push(k.to_string()); }
        }
    }
    if let Some(p) = keys.iter().position(|k| k == "GT") {
        if p != 0 {
            let gt = keys.remove(p);
            keys.insert(0, gt);
        }
    }
    keys.join(":")
}

/// INFO in bcftools order: plain tags as first seen, then rule tags
/// (alphabetical), then A/R/G tags, then AN/AC recomputed from the genotypes.
fn merge_info(
    group: &[(usize, Record)],
    maps: &[Vec<Option<usize>>],
    n_new: usize,
    ctx: &MergeCtx<'_>,
    gt_idx: Option<usize>,
    sample_cols: &[String],
) -> String {
    let hdr = ctx.hdr;
    let is_float = |k: &str| matches!(hdr.info_type(k), Some(FieldType::Float));
    let mut plain: Vec<(String, Option<String>)> = Vec::new();
    let mut rule_vals: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut agr: Vec<(String, String)> = Vec::new();
    for (gi, (_, r)) in group.iter().enumerate() {
        for (k, v) in split_info(&r.info) {
            if !ctx.rules.keep_ac_an && (k == "AC" || k == "AN") && gt_idx.is_some() { continue; }
            let num = hdr.info_number(k);
            if let Some(rule) = ctx.rules.rules.get(k) {
                if let Some(v) = v {
                    // A/R/G vectors are remapped before folding; join keeps the raw text.
                    let val = if num.is_per_allele() && !matches!(rule, InfoRule::Join) {
                        remap_value(v, num, maps[gi].len(), n_new, &maps[gi]).unwrap_or_else(|| v.to_string())
                    } else {
                        v.to_string()
                    };
                    rule_vals.entry(k.to_string()).or_default().push(val);
                }
                continue;
            }
            if num.is_per_allele() {
                let val = match v {
                    Some(v) => remap_value(v, num, maps[gi].len(), n_new, &maps[gi]).unwrap_or_else(|| v.to_string()),
                    None => continue,
                };
                match agr.iter_mut().find(|(kk, _)| kk == k) {
                    Some(e) => e.1 = fill_missing(&e.1, &val),
                    None => agr.push((k.to_string(), val)),
                }
                continue;
            }
            if !plain.iter().any(|(kk, _)| kk == k) {
                let val = v.map(|v| if is_float(k) { fmt_float_list(v) } else { v.to_string() });
                plain.push((k.to_string(), val));
            }
        }
    }

    let mut items: Vec<(String, Option<String>)> = plain;
    for (k, vals) in &rule_vals {
        let rule = &ctx.rules.rules[k];
        let present: Vec<&str> = vals.iter().map(String::as_str).collect();
        let merged = match rule {
            InfoRule::Sum => vector_fold_fill(&present, |a, b| a + b, Some(0.0)),
            InfoRule::Avg => {
                let n = present.len() as f64;
                let s = vector_fold_fill(&present, |a, b| a + b, Some(0.0));
                s.split(',').map(|x| x.parse::<f64>().map(|v| normalize_num(v / n)).unwrap_or_else(|_| x.to_string())).collect::<Vec<_>>().join(",")
            }
            InfoRule::Min => vector_fold_fill(&present, f64::min, None),
            InfoRule::Max => vector_fold_fill(&present, f64::max, None),
            InfoRule::Join => present.join(","),
            InfoRule::First => present[0].to_string(),
        };
        let merged = if is_float(k) { fmt_float_list(&merged) } else { merged };
        items.push((k.clone(), Some(merged)));
    }
    for (k, v) in agr {
        let v = if is_float(&k) { fmt_float_list(&v) } else { v };
        items.push((k, Some(v)));
    }

    // AN/AC follow the merged genotypes.
    if let (Some(gi), false) = (gt_idx, ctx.rules.keep_ac_an) {
        if hdr.info.contains_key("AC") || hdr.info.contains_key("AN") {
            let mut ac = vec![0u32; n_new];
            let mut an = 0u32;
            let mut any_gt = false;
            for s in sample_cols {
                let gt = s.split(':').nth(gi).unwrap_or(".");
                for a in gt_alleles(gt).into_iter().flatten() {
                    any_gt = true;
                    an += 1;
                    if a < ac.len() { ac[a] += 1; }
                }
            }
            if any_gt || !sample_cols.is_empty() {
                items.retain(|(k, _)| k != "AN" && k != "AC");
                if hdr.info.contains_key("AN") { items.push(("AN".into(), Some(an.to_string()))); }
                if hdr.info.contains_key("AC") && n_new > 1 {
                    items.push(("AC".into(), Some(ac[1..].iter().map(u32::to_string).collect::<Vec<_>>().join(","))));
                }
            }
        }
    }
    join_info(&items)
}

fn fill_missing(cur: &str, new: &str) -> String {
    let a: Vec<&str> = cur.split(',').collect();
    let b: Vec<&str> = new.split(',').collect();
    if a.len() != b.len() { return cur.to_string(); }
    a.iter().zip(b.iter()).map(|(x, y)| if *x == "." { *y } else { *x }).collect::<Vec<_>>().join(",")
}

/// Element-wise fold over comma vectors (scalars are length-1 vectors).
#[cfg(test)]
fn vector_fold(vals: &[&str], f: fn(f64, f64) -> f64) -> String {
    vector_fold_fill(vals, f, None)
}

/// Element-wise fold of comma-separated vectors; slots with no numeric value
/// print `fill` when given.
fn vector_fold_fill(vals: &[&str], f: fn(f64, f64) -> f64, fill: Option<f64>) -> String {
    let mut acc: Vec<Option<f64>> = Vec::new();
    for v in vals {
        for (i, x) in v.split(',').enumerate() {
            if acc.len() <= i { acc.push(None); }
            if let Ok(n) = x.parse::<f64>() {
                acc[i] = Some(match acc[i] { Some(c) => f(c, n), None => n });
            }
        }
    }
    acc.iter().map(|x| x.or(fill).map(normalize_num).unwrap_or_else(|| ".".into())).collect::<Vec<_>>().join(",")
}

fn normalize_num(v: f64) -> String {
    if v.fract() == 0.0 && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        format!("{:.6}", v).trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_merge.rs"]
mod tests;
