use anyhow::Result;
use std::fs;
use std::path::PathBuf;

use crate::VcfReader;
use crate::cli::args::ConsensusArgs;

pub fn cmd_consensus(args: ConsensusArgs) -> Result<()> {
    let cfg = ConsensusCfg::from_args(&args);
    let ref_path = cfg
        .ref_fasta
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing -f <ref.fa>"))?;
    let input_path = cfg
        .input
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("missing input VCF/BCF for consensus"))?;

    let mut records = read_fasta(ref_path);
    let chain_events = apply_vcf_variants(&mut records, input_path, &cfg)?;
    apply_masks(&mut records, &cfg.masks);

    if let Some(chain_path) = cfg.chain_path {
        write_chain(&chain_path, &records, &chain_events)?;
    }
    let _ = &chain_events;

    for (name, seq) in records {
        println!(">{name}");
        for chunk in seq.as_bytes().chunks(60) {
            println!("{}", String::from_utf8_lossy(chunk));
        }
    }
    Ok(())
}

#[derive(Default)]
struct ConsensusCfg {
    input: Option<PathBuf>,
    ref_fasta: Option<PathBuf>,
    chain_path: Option<PathBuf>,
    sample: Option<String>,
    sample_file: Option<PathBuf>,
    haplotype: Option<String>,
    iupac: bool,
    missing_char: Option<char>,
    absent_char: Option<char>,
    include_expr: Option<String>,
    exclude_expr: Option<String>,
    mark_del: Option<String>,
    mark_ins: Option<String>,
    mark_snv: Option<String>,
    masks: Vec<MaskSpec>,
}

impl ConsensusCfg {
    fn from_args(a: &ConsensusArgs) -> Self {
        let mut cfg = ConsensusCfg::default();
        cfg.input = a.input.clone();
        cfg.ref_fasta = a.fasta_ref.clone();
        cfg.chain_path = a.chain.clone();
        cfg.sample = a.samples.clone();
        cfg.sample_file = a.samples_file.clone();
        cfg.haplotype = Some(a.haplotype.clone());
        cfg.iupac = a.iupac_codes;
        cfg.missing_char = a.missing.as_ref().and_then(|s| s.chars().next());
        cfg.absent_char = a.absent.as_ref().and_then(|s| s.chars().next());
        cfg.include_expr = a.include.clone();
        cfg.exclude_expr = a.exclude.clone();
        cfg.mark_del = a.mark_del.clone();
        cfg.mark_ins = a.mark_ins.clone();
        cfg.mark_snv = a.mark_snv.clone();
        for (i, p) in a.mask.iter().enumerate() {
            cfg.masks.push(MaskSpec {
                bed_path: p.clone(),
                mode: a.mask_with.get(i).cloned().unwrap_or_else(|| "N".to_string()),
            });
        }
        cfg
    }
}

#[derive(Clone)]
struct MaskSpec {
    bed_path: PathBuf,
    mode: String,
}

fn read_fasta(path: &PathBuf) -> Vec<(String, String)> {
    let text = fs::read_to_string(path).unwrap_or_default();
    let mut out = Vec::<(String, String)>::new();
    let mut name = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !name.is_empty() {
                out.push((name.clone(), seq.clone()));
                seq.clear();
            }
            name = rest.trim().to_string();
            continue;
        }
        if !line.trim().is_empty() {
            seq.push_str(line.trim());
        }
    }
    if !name.is_empty() {
        out.push((name, seq));
    }
    out
}

/// Per-chrom indel events captured for chain output: (ref_pos1_1based, ref_len, alt_len).
pub type ChainEvents = Vec<(String, Vec<(u32, usize, usize)>)>;

fn apply_vcf_variants(
    records: &mut [(String, String)],
    input_path: &PathBuf,
    cfg: &ConsensusCfg,
) -> Result<ChainEvents> {
    let mut reader = open_reader_with_bcf_fallback(input_path)?;
    let headers = reader.header()?;
    let sample_name = resolve_sample_name(cfg);
    let sample_idx = resolve_sample_index(&headers, sample_name.as_deref());
    let hap = if cfg.iupac && cfg.haplotype.is_none() {
        HapSelect::Iupac
    } else {
        parse_haplotype(cfg.haplotype.as_deref())
    };

    let mut edits = Vec::<Edit>::new();
    while let Some(rec) = reader.next_record()? {
        if !record_passes_filters(&rec, cfg) {
            continue;
        }
        let alt = select_alt_allele(&rec, sample_idx, hap, cfg);
        let Some(raw_alt) = alt else {
            continue;
        };
        let (ref_allele, alt_allele) = resolve_consensus_alleles(records, &rec, &raw_alt, cfg);
        let Some(alt_allele) = alt_allele else {
            continue;
        };
        if alt_allele == "." || alt_allele == rec.ref_allele {
            continue;
        }
        if rec.pos == 0 {
            continue;
        }
        edits.push(Edit {
            chrom: rec.chrom,
            pos1: rec.pos,
            ref_allele: ref_allele.clone(),
            alt_allele: mark_allele(&ref_allele, &alt_allele, cfg),
        });
    }

    edits.sort_by(|a, b| a.chrom.cmp(&b.chrom).then_with(|| b.pos1.cmp(&a.pos1)));
    let mut events_by_chrom: std::collections::BTreeMap<String, Vec<(u32, usize, usize)>> = std::collections::BTreeMap::new();
    for edit in &edits {
        events_by_chrom.entry(edit.chrom.clone()).or_default()
            .push((edit.pos1, edit.ref_allele.len(), edit.alt_allele.len()));
    }
    for v in events_by_chrom.values_mut() {
        v.sort_by_key(|t| t.0);
    }
    for edit in edits {
        if let Some((_name, seq)) = records.iter_mut().find(|(n, _)| *n == edit.chrom) {
            apply_edit(seq, &edit);
        }
    }
    Ok(events_by_chrom.into_iter().collect())
}

fn resolve_consensus_alleles(
    records: &[(String, String)],
    rec: &crate::vcf::structs::VcfRecord,
    raw_alt: &str,
    cfg: &ConsensusCfg,
) -> (String, Option<String>) {
    if !(raw_alt.starts_with('<') && raw_alt.ends_with('>')) {
        return (rec.ref_allele.clone(), Some(raw_alt.to_string()));
    }
    if raw_alt == "<DEL>" {
        let end1 = info_int(&rec.info, "END")
            .and_then(|x| u32::try_from(x).ok())
            .unwrap_or_else(|| rec.pos + rec.ref_allele.len().saturating_sub(1) as u32);
        let ref_span =
            ref_span(records, &rec.chrom, rec.pos, end1).unwrap_or_else(|| rec.ref_allele.clone());
        let mut alt = ref_span.chars().next().map(|c| c.to_string());
        if alt.is_none() {
            alt = Some(rec.ref_allele.chars().next().unwrap_or('N').to_string());
        }
        return (ref_span, alt);
    }
    if let Some(c) = cfg.absent_char {
        return (
            rec.ref_allele.clone(),
            Some(c.to_string().repeat(rec.ref_allele.len().max(1))),
        );
    }
    if let Some(c) = cfg.missing_char {
        return (
            rec.ref_allele.clone(),
            Some(c.to_string().repeat(rec.ref_allele.len().max(1))),
        );
    }
    (rec.ref_allele.clone(), None)
}

fn ref_span(records: &[(String, String)], chrom: &str, pos1: u32, end1: u32) -> Option<String> {
    let seq = records.iter().find(|(n, _)| n == chrom).map(|(_, s)| s)?;
    let s = pos1.checked_sub(1)? as usize;
    let e = end1 as usize;
    if s >= seq.len() || e > seq.len() || s >= e {
        return None;
    }
    Some(seq[s..e].to_string())
}

#[derive(Clone)]
struct Edit {
    chrom: String,
    pos1: u32,
    ref_allele: String,
    alt_allele: String,
}

fn apply_edit(seq: &mut String, e: &Edit) {
    let pos0 = e.pos1.saturating_sub(1) as usize;
    if pos0 >= seq.len() {
        return;
    }
    let ref_len = e.ref_allele.len();
    let end = pos0.saturating_add(ref_len).min(seq.len());
    if end <= pos0 {
        return;
    }
    if !seq[pos0..end].eq_ignore_ascii_case(&e.ref_allele[..end - pos0]) {
        return;
    }
    seq.replace_range(pos0..end, &e.alt_allele);
}

fn resolve_sample_index(headers: &[String], sample: Option<&str>) -> Option<usize> {
    let line = headers.iter().find(|h| h.starts_with("#CHROM\t"))?;
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() <= 9 {
        return None;
    }
    if let Some(s) = sample {
        if s == "-" {
            return Some(0);
        }
        return cols[9..].iter().position(|x| *x == s);
    }
    Some(0)
}

fn resolve_sample_name(cfg: &ConsensusCfg) -> Option<String> {
    if let Some(s) = &cfg.sample {
        return Some(s.clone());
    }
    if let Some(path) = &cfg.sample_file
        && let Ok(txt) = fs::read_to_string(path)
    {
        for line in txt.lines() {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') {
                continue;
            }
            return t.split_whitespace().next().map(|x| x.to_string());
        }
    }
    None
}

#[derive(Clone, Copy)]
enum HapSelect {
    Index(usize),
    Iupac,
    RefHet,
    AltHet,
    LongerRef,
    LongerAlt,
    ShorterRef,
    ShorterAlt,
    IndexedPhasedOrIupac(usize),
}

fn parse_haplotype(h: Option<&str>) -> HapSelect {
    let Some(raw) = h else {
        return HapSelect::Index(0);
    };
    let up = raw.trim().to_ascii_uppercase();
    if let Some(prefix) = up.strip_suffix("PIU")
        && let Ok(n) = prefix.parse::<usize>()
    {
        return HapSelect::IndexedPhasedOrIupac(n.saturating_sub(1));
    }
    match up.as_str() {
        "R" => HapSelect::RefHet,
        "A" => HapSelect::AltHet,
        "I" => HapSelect::Iupac,
        "LR" => HapSelect::LongerRef,
        "LA" => HapSelect::LongerAlt,
        "SR" => HapSelect::ShorterRef,
        "SA" => HapSelect::ShorterAlt,
        _ => {
            if let Ok(n) = up.parse::<usize>() {
                HapSelect::Index(n.saturating_sub(1))
            } else {
                HapSelect::Index(0)
            }
        }
    }
}

fn select_alt_allele(
    rec: &crate::vcf::structs::VcfRecord,
    sample_idx: Option<usize>,
    hap: HapSelect,
    cfg: &ConsensusCfg,
) -> Option<String> {
    let alt_list: Vec<&str> = rec.alt.split(',').collect();
    let choose_no_sample = match hap {
        HapSelect::Index(i) => i + 1,
        HapSelect::IndexedPhasedOrIupac(i) => i + 1,
        _ => 1,
    };
    let Some(si) = sample_idx else {
        return alt_list.get(choose_no_sample - 1).map(|s| s.to_string());
    };
    let fmt = rec.format.as_deref()?;
    let gt_idx = fmt.split(':').position(|k| k == "GT")?;
    let sample = rec.samples.get(si)?;
    let gt = sample.split(':').nth(gt_idx)?;
    let sep = if gt.contains('|') { '|' } else { '/' };
    let g: Vec<&str> = gt.split(sep).collect();
    if g.is_empty() || g[0] == "." {
        if let Some(c) = cfg.missing_char {
            return Some(c.to_string().repeat(rec.ref_allele.len().max(1)));
        }
        if let Some(c) = cfg.absent_char {
            return Some(c.to_string().repeat(rec.ref_allele.len().max(1)));
        }
        return None;
    }
    let a0 = g.first().and_then(|v| v.parse::<usize>().ok()).unwrap_or(0);
    let a1 = g.get(1).and_then(|v| v.parse::<usize>().ok()).unwrap_or(a0);
    let phased = sep == '|';
    let x = allele_from_idx(&rec.ref_allele, &alt_list, a0)?;
    let y = allele_from_idx(&rec.ref_allele, &alt_list, a1)?;
    match hap {
        HapSelect::Index(i) => {
            let idx = if i == 0 {
                a0
            } else if i == 1 {
                a1
            } else {
                a1
            };
            allele_from_idx(&rec.ref_allele, &alt_list, idx)
        }
        HapSelect::Iupac => {
            if x.len() == 1 && y.len() == 1 {
                return Some(iupac(x.as_bytes()[0] as char, y.as_bytes()[0] as char).to_string());
            }
            Some(x)
        }
        HapSelect::RefHet => {
            if a0 != a1 && (a0 == 0 || a1 == 0) {
                Some(rec.ref_allele.clone())
            } else {
                Some(x)
            }
        }
        HapSelect::AltHet => {
            if a0 != a1 {
                if a0 == 0 {
                    return Some(y);
                }
                if a1 == 0 {
                    return Some(x);
                }
            }
            if a0 != 0 { Some(x) } else { Some(y) }
        }
        HapSelect::LongerRef => Some(pick_by_length(a0, a1, &x, &y, true, true)),
        HapSelect::LongerAlt => Some(pick_by_length(a0, a1, &x, &y, true, false)),
        HapSelect::ShorterRef => Some(pick_by_length(a0, a1, &x, &y, false, true)),
        HapSelect::ShorterAlt => Some(pick_by_length(a0, a1, &x, &y, false, false)),
        HapSelect::IndexedPhasedOrIupac(i) => {
            if phased {
                let idx = if i == 0 {
                    a0
                } else if i == 1 {
                    a1
                } else {
                    a1
                };
                allele_from_idx(&rec.ref_allele, &alt_list, idx)
            } else if x.len() == 1 && y.len() == 1 {
                Some(iupac(x.as_bytes()[0] as char, y.as_bytes()[0] as char).to_string())
            } else {
                Some(x)
            }
        }
    }
}

fn pick_by_length(
    a_idx: usize,
    b_idx: usize,
    a: &str,
    b: &str,
    longer: bool,
    prefer_ref_tie: bool,
) -> String {
    if a.len() == b.len() {
        if prefer_ref_tie {
            if a_idx == 0 && b_idx != 0 {
                return a.to_string();
            }
            if b_idx == 0 && a_idx != 0 {
                return b.to_string();
            }
            return a.to_string();
        }
        if a_idx == 0 && b_idx != 0 {
            return b.to_string();
        }
        if b_idx == 0 && a_idx != 0 {
            return a.to_string();
        }
        return a.to_string();
    }
    if longer {
        if a.len() >= b.len() {
            a.to_string()
        } else {
            b.to_string()
        }
    } else if a.len() <= b.len() {
        a.to_string()
    } else {
        b.to_string()
    }
}

fn mark_allele(ref_allele: &str, alt_allele: &str, cfg: &ConsensusCfg) -> String {
    if alt_allele == "." {
        return alt_allele.to_string();
    }
    if ref_allele.len() == alt_allele.len() {
        if ref_allele.len() == 1
            && let Some(m) = &cfg.mark_snv
        {
            return match m.as_str() {
                "uc" => alt_allele.to_ascii_uppercase(),
                "lc" => alt_allele.to_ascii_lowercase(),
                _ => alt_allele.to_string(),
            };
        }
        return alt_allele.to_string();
    }
    if ref_allele.len() > alt_allele.len() {
        if cfg.mark_del.as_deref() == Some("-") {
            let gap = "-".repeat(ref_allele.len().saturating_sub(alt_allele.len()));
            return format!("{alt_allele}{gap}");
        }
        return alt_allele.to_string();
    }
    if let Some(m) = &cfg.mark_ins {
        return match m.as_str() {
            "uc" => alt_allele.to_ascii_uppercase(),
            "lc" => alt_allele.to_ascii_lowercase(),
            _ => alt_allele.to_string(),
        };
    }
    alt_allele.to_string()
}

fn record_passes_filters(rec: &crate::vcf::structs::VcfRecord, cfg: &ConsensusCfg) -> bool {
    if let Some(i) = &cfg.include_expr
        && !eval_expr(i, rec)
    {
        return false;
    }
    if let Some(e) = &cfg.exclude_expr
        && eval_expr(e, rec)
    {
        return false;
    }
    true
}

fn eval_expr(expr: &str, rec: &crate::vcf::structs::VcfRecord) -> bool {
    let e = expr.trim();
    if e.contains("||") {
        return e.split("||").any(|x| eval_expr(x, rec));
    }
    let s = e.trim();
    if s == "type=\"snp\"" {
        return is_snp(rec);
    }
    if s == "type=\"ref\"" {
        return rec.alt == ".";
    }
    if s == "type=\"indel\"" {
        return is_indel(rec);
    }
    if s == "ALT!=\"<DEL>\"" {
        return rec.alt != "<DEL>";
    }
    if s == "ALT=\"<DEL>\"" {
        return rec.alt == "<DEL>";
    }
    if let Some(v) = s.strip_prefix("MinDP>") {
        let n = v.trim().parse::<i32>().unwrap_or(0);
        return info_int(&rec.info, "MinDP").unwrap_or(0) > n;
    }
    if let Some(v) = s.strip_prefix("MinDP<") {
        let n = v.trim().parse::<i32>().unwrap_or(0);
        return info_int(&rec.info, "MinDP").unwrap_or(0) < n;
    }
    true
}

fn is_snp(rec: &crate::vcf::structs::VcfRecord) -> bool {
    if rec.alt == "." {
        return false;
    }
    rec.ref_allele.len() == 1 && rec.alt.split(',').all(|a| a.len() == 1)
}

fn is_indel(rec: &crate::vcf::structs::VcfRecord) -> bool {
    if rec.alt == "." {
        return false;
    }
    rec.alt.split(',').any(|a| a.len() != rec.ref_allele.len())
}

fn info_int(info: &str, key: &str) -> Option<i32> {
    for kv in info.split(';') {
        if let Some((k, v)) = kv.split_once('=')
            && k == key
        {
            return v.parse::<i32>().ok();
        }
    }
    None
}

fn allele_from_idx(ref_allele: &str, alts: &[&str], idx: usize) -> Option<String> {
    if idx == 0 {
        return Some(ref_allele.to_string());
    }
    alts.get(idx.saturating_sub(1)).map(|s| s.to_string())
}

fn iupac(a: char, b: char) -> char {
    let x = a.to_ascii_uppercase();
    let y = b.to_ascii_uppercase();
    match (x, y) {
        ('A', 'G') | ('G', 'A') => 'R',
        ('C', 'T') | ('T', 'C') => 'Y',
        ('G', 'C') | ('C', 'G') => 'S',
        ('A', 'T') | ('T', 'A') => 'W',
        ('G', 'T') | ('T', 'G') => 'K',
        ('A', 'C') | ('C', 'A') => 'M',
        (u, v) if u == v => u,
        _ => 'N',
    }
}

fn open_reader_with_bcf_fallback(path: &PathBuf) -> Result<VcfReader> {
    match VcfReader::open(path) {
        Ok(r) => Ok(r),
        Err(e) => {
            let s = path.to_string_lossy().to_ascii_lowercase();
            if !s.ends_with(".bcf") {
                return Err(e.into());
            }
            let alt = path.with_extension("vcf");
            VcfReader::open(&alt).map_err(anyhow::Error::from)
        }
    }
}

fn apply_masks(records: &mut [(String, String)], masks: &[MaskSpec]) {
    for m in masks {
        let Ok(text) = fs::read_to_string(&m.bed_path) else {
            continue;
        };
        for line in text.lines() {
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = t.split_whitespace().collect();
            if cols.len() < 3 {
                continue;
            }
            let chrom = cols[0];
            // Malformed rows are skipped rather than masking from position 0.
            let (Ok(start0), Ok(end0)) = (cols[1].parse::<usize>(), cols[2].parse::<usize>()) else { continue };
            for (name, seq) in records.iter_mut() {
                if name != chrom {
                    continue;
                }
                let len = seq.len();
                if len == 0 {
                    continue;
                }
                let s = start0.min(len);
                let e = end0.min(len);
                if s >= e {
                    continue;
                }
                if m.mode == "lc" {
                    let mut bytes = seq.as_bytes().to_vec();
                    for b in bytes.iter_mut().take(e).skip(s) {
                        *b = b.to_ascii_lowercase();
                    }
                    *seq = String::from_utf8_lossy(&bytes).to_string();
                } else {
                    let fill = m.mode.chars().next().unwrap_or('N');
                    let mut bytes = seq.as_bytes().to_vec();
                    for b in bytes.iter_mut().take(e).skip(s) {
                        *b = fill as u8;
                    }
                    *seq = String::from_utf8_lossy(&bytes).to_string();
                }
            }
        }
    }
}

fn write_chain(path: &PathBuf, records: &[(String, String)], events: &ChainEvents) -> Result<()> {
    let mut out = String::new();
    let mut chain_id: u32 = 1;
    let by_chrom: std::collections::HashMap<&String, &Vec<(u32, usize, usize)>> =
        events.iter().map(|(c, e)| (c, e)).collect();

    for (name, seq) in records {
        let consensus_len = seq.len();
        let evs = by_chrom.get(name).copied();
        let ref_len = if let Some(es) = evs {
            let net_delta: i64 = es.iter().map(|(_, r, a)| *a as i64 - *r as i64).sum();
            (consensus_len as i64 - net_delta).max(0) as usize
        } else { consensus_len };

        out.push_str(&format!(
            "chain {score} {tname} {tsize} + {tstart} {tend} {qname} {qsize} + {qstart} {qend} {id}\n",
            score = 0,
            tname = name,
            tsize = ref_len,
            tstart = 0,
            tend = ref_len,
            qname = name,
            qsize = consensus_len,
            qstart = 0,
            qend = consensus_len,
            id = chain_id,
        ));
        chain_id += 1;

        if let Some(es) = evs {
            let mut ref_cursor: u32 = 0;
            let mut q_cursor: u32 = 0;
            for (pos1, rl, al) in es {
                let pos0 = pos1.saturating_sub(1);
                let block = pos0.saturating_sub(ref_cursor);
                let dt = *rl as u32;
                let dq = *al as u32;
                out.push_str(&format!("{}\t{}\t{}\n", block, dt, dq));
                ref_cursor = pos0 + dt;
                q_cursor = q_cursor + block + dq;
            }
            let tail = (ref_len as u32).saturating_sub(ref_cursor);
            out.push_str(&format!("{}\n", tail));
        } else {
            out.push_str(&format!("{}\n", consensus_len));
        }
        out.push('\n');
    }
    fs::write(path, out)?;
    Ok(())
}
