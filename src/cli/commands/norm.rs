use crate::annotate::postproc::{OutputKind, parse_output_type, version_header_line};
use crate::cli::args::NormArgs;
use crate::vcf::UnifiedVcfReader;
use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum MultiMode { None, SplitAll, SplitSnps, SplitIndels, JoinAll, JoinSnps, JoinIndels }

impl MultiMode {
    fn parse(s: Option<&str>) -> Result<Self> {
        let Some(s) = s else { return Ok(Self::None); };
        Ok(match s {
            "-" | "-any" | "-both" => Self::SplitAll,
            "-snps" => Self::SplitSnps,
            "-indels" => Self::SplitIndels,
            "+" | "+any" | "+both" => Self::JoinAll,
            "+snps" => Self::JoinSnps,
            "+indels" => Self::JoinIndels,
            other => bail!("--multiallelics: unknown {other:?}"),
        })
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum CheckRef { Exit, Warn, Exclude, Set }

impl CheckRef {
    fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "e" => Self::Exit, "w" => Self::Warn, "x" => Self::Exclude, "s" => Self::Set,
            o => bail!("-c, --check-ref: unknown {o:?} (expected e|w|x|s)"),
        })
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum RmDup { None, All, Exact, Snps, Indels, Both, AnyAllele, AnyId }

impl RmDup {
    fn parse(s: Option<&str>) -> Result<Self> {
        let Some(s) = s else { return Ok(Self::None); };
        Ok(match s {
            "none" => Self::None,
            "exact" => Self::Exact,
            "snps" => Self::Snps,
            "indels" => Self::Indels,
            "both" => Self::Both,
            "any" | "all" => Self::All,
            "id" => Self::AnyId,
            o => bail!("--rm-dup: unknown {o:?} (none|exact|snps|indels|both|all|id)"),
        })
    }
}

#[derive(Clone, Debug)]
struct SplitOpts {
    multi_overlaps: String,
    keep_sum_keys: Vec<String>,
    old_rec_tag: Option<String>,
}

pub fn cmd_norm(args: NormArgs) -> Result<()> {
    let multi = MultiMode::parse(args.multiallelics.as_deref())?;
    let check_ref = CheckRef::parse(&args.check_ref)?;
    let rm_dup = RmDup::parse(args.rm_dup.as_deref())?;
    let split_opts = SplitOpts {
        multi_overlaps: args.multi_overlaps.clone(),
        keep_sum_keys: args.keep_sum.as_deref()
            .map(|s| s.split(',').map(|t| t.trim().to_string()).collect())
            .unwrap_or_default(),
        old_rec_tag: args.old_rec_tag.clone(),
    };

    let fasta = args.fasta_ref.as_ref().map(|p| load_fasta(p)).transpose()?;
    if matches!(check_ref, CheckRef::Exit | CheckRef::Warn | CheckRef::Exclude | CheckRef::Set) && fasta.is_none() && !args.do_not_normalize {
        if args.fasta_ref.is_none() && !matches!(multi, MultiMode::SplitAll | MultiMode::SplitSnps | MultiMode::SplitIndels | MultiMode::JoinAll | MultiMode::JoinSnps | MultiMode::JoinIndels) && !args.atomize {
            // ok: no fasta needed when only filtering
        }
    }

    let out_path = args.output.clone().unwrap_or_else(|| {
        let mut p = args.input.clone(); p.set_extension("norm.vcf"); p
    });
    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);
    let mut sink = open_sink(&out_path, kind)?;

    let mut reader = UnifiedVcfReader::open(&args.input).context("open input")?;
    let headers = reader.header()?;

    // VCF spec: `#CHROM` MUST be the final header line. Insert the kira version
    // line just before it; otherwise downstream `bcftools` parses it as a record
    // and bails on the malformed chrom name.
    let version = version_header_line();
    let mut wrote_version = false;
    for h in &headers {
        if h.starts_with("#CHROM") && !args.no_version && !wrote_version {
            sink.write_all(version.as_bytes())?; sink.write_all(b"\n")?;
            wrote_version = true;
        }
        sink.write_all(h.as_bytes())?; sink.write_all(b"\n")?;
    }
    if !args.no_version && !wrote_version {
        // No `#CHROM` in input (degenerate) — append at end as best-effort.
        sink.write_all(version.as_bytes())?; sink.write_all(b"\n")?;
    }

    let mut dup_window: Vec<(u32, Vec<String>)> = Vec::new();
    let mut current_chrom = String::new();

    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
        let records = expand_record(&line, multi, args.atomize, &split_opts)?;
        for rec in records {
            let cols: Vec<&str> = rec.split('\t').collect();
            if cols.len() < 8 { continue; }
            if cols[0] != current_chrom { current_chrom = cols[0].to_string(); dup_window.clear(); }
            let pos: u32 = cols[1].parse().unwrap_or(0);

            if let Some(fa) = &fasta {
                match verify_ref(&rec, fa, check_ref) {
                    RefAction::Keep => {}
                    RefAction::Skip => continue,
                    RefAction::Fix(new) => {
                        let mut new_cols = cols.clone();
                        new_cols[3] = unsafe { &*(new.as_str() as *const str) };
                        // Build manually so the borrow stays alive
                        let fixed = new_cols.join("\t");
                        if !args.do_not_normalize && fa.has(cols[0]) {
                            if let Some(left) = left_align(&fixed, fa) {
                                if !is_dup(&dup_window, pos, &left, rm_dup) {
                                    sink.write_all(left.as_bytes())?; sink.write_all(b"\n")?;
                                    push_dup(&mut dup_window, pos, left);
                                }
                                continue;
                            }
                        }
                        if !is_dup(&dup_window, pos, &fixed, rm_dup) {
                            sink.write_all(fixed.as_bytes())?; sink.write_all(b"\n")?;
                            push_dup(&mut dup_window, pos, fixed);
                        }
                        continue;
                    }
                    RefAction::Fail(msg) => bail!("-c e: {msg}"),
                    RefAction::Warn(msg) => { eprintln!("[norm] warn: {msg}"); }
                }
            }

            let final_rec = if !args.do_not_normalize {
                if let Some(fa) = &fasta {
                    if fa.has(cols[0]) { left_align(&rec, fa).unwrap_or(rec) } else { rec }
                } else { rec }
            } else { rec };

            if !is_dup(&dup_window, pos, &final_rec, rm_dup) {
                sink.write_all(final_rec.as_bytes())?; sink.write_all(b"\n")?;
                push_dup(&mut dup_window, pos, final_rec);
            }
        }
    }
    sink.flush()?;
    Ok(())
}

fn open_sink(p: &Path, kind: OutputKind) -> Result<Box<dyn Write>> {
    match kind {
        OutputKind::Vcf => Ok(Box::new(BufWriter::with_capacity(1 << 20, File::create(p)?))),
        OutputKind::VcfGz(lvl) => {
            let w = crate::bgzf::BgzfWriter::with_compression(p, flate2::Compression::new(lvl))?;
            Ok(Box::new(w))
        }
        OutputKind::Bcf(_) => bail!("-O u|b (BCF) not yet supported in norm"),
    }
}

fn expand_record(line: &str, multi: MultiMode, atomize: bool, opts: &SplitOpts) -> Result<Vec<String>> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 { return Ok(vec![line.to_string()]); }
    let refa = cols[3]; let alt = cols[4];
    let alts: Vec<&str> = alt.split(',').collect();

    let split_now = match multi {
        MultiMode::SplitAll => alts.len() > 1,
        MultiMode::SplitSnps => alts.len() > 1 && alts.iter().all(|a| a.len() == 1 && refa.len() == 1),
        MultiMode::SplitIndels => alts.len() > 1 && alts.iter().any(|a| a.len() != refa.len()),
        _ => false,
    };

    let mut result: Vec<String> = if split_now {
        let info = cols[7];
        let orig_pos = cols[1];
        let orig_ref = cols[3];
        let orig_alt = cols[4];
        let mut out_recs: Vec<String> = Vec::with_capacity(alts.len());
        for i in 0..alts.len() {
            if alts[i] == "*" {
                match opts.multi_overlaps.as_str() {
                    "0" => continue,
                    "." => {}
                    _ => {}
                }
            }
            let mut c = cols.clone();
            let alt_eff: String = if alts[i] == "*" && opts.multi_overlaps == "." {
                ".".to_string()
            } else { alts[i].to_string() };
            let alt_static: &str = unsafe { &*(alt_eff.as_str() as *const str) };
            c[4] = alt_static;
            let mut new_info = split_info_per_allele_keep_sum(info, i, alts.len(), &opts.keep_sum_keys);
            if let Some(tag) = &opts.old_rec_tag {
                let extra = format!("{}={}|{}|{}", tag, orig_pos, orig_ref, orig_alt);
                if new_info == "." || new_info.is_empty() { new_info = extra; }
                else { new_info.push(';'); new_info.push_str(&extra); }
            }
            let mut s = String::new();
            for (idx, col) in c.iter().enumerate() {
                if idx > 0 { s.push('\t'); }
                if idx == 7 { s.push_str(&new_info); } else { s.push_str(col); }
            }
            if cols.len() > 8 {
                let new_samples = split_samples_per_allele(&cols[8..], i, alts.len());
                for col in &new_samples {
                    s.push('\t'); s.push_str(col);
                }
            }
            out_recs.push(s);
        }
        out_recs
    } else { vec![line.to_string()] };

    if atomize {
        let mut atomic = Vec::new();
        for rec in result {
            atomic.extend(atomize_record(&rec));
        }
        result = atomic;
    }

    Ok(result)
}

/// Stream-level MNP joining: collapse runs of adjacent SNVs at positions
/// p, p+1, p+2, ... (same REF length=1, same ALT length=1, same FILTER)
/// into single MNP record at position p. Used when `multi` is JoinSnps.
pub fn join_adjacent_snps(buffer: &mut Vec<String>) -> Vec<String> {
    if buffer.len() < 2 { return std::mem::take(buffer); }
    let mut out: Vec<String> = Vec::new();
    let mut group: Vec<String> = Vec::new();
    let mut group_last_pos: u32 = 0;

    let flush = |group: &mut Vec<String>, out: &mut Vec<String>| {
        if group.len() <= 1 {
            out.extend(group.drain(..));
            return;
        }
        let first_cols: Vec<&str> = group[0].split('\t').collect();
        let mut new_ref = String::new();
        let mut new_alts: Vec<String> = first_cols[4].split(',').map(|s| s.to_string()).collect();
        for rec in group.iter() {
            let c: Vec<&str> = rec.split('\t').collect();
            new_ref.push_str(c[3]);
            let alts: Vec<&str> = c[4].split(',').collect();
            for (i, a) in alts.iter().enumerate() {
                if i >= new_alts.len() { new_alts.push(String::new()); }
                if rec.as_ptr() == group[0].as_ptr() { continue; }
                new_alts[i].push_str(a);
            }
        }
        let mut joined: Vec<String> = first_cols.iter().map(|s| s.to_string()).collect();
        joined[3] = new_ref;
        joined[4] = new_alts.join(",");
        out.push(joined.join("\t"));
        group.clear();
    };

    for rec in buffer.drain(..) {
        let cols: Vec<&str> = rec.split('\t').collect();
        if cols.len() < 8 { out.push(rec); continue; }
        let pos: u32 = cols[1].parse().unwrap_or(0);
        let refa = cols[3]; let alt = cols[4];
        let is_snp = refa.len() == 1 && alt.split(',').all(|a| a.len() == 1 && a != ".");
        if !is_snp {
            if !group.is_empty() {
                let mut g = std::mem::take(&mut group);
                flush(&mut g, &mut out);
            }
            out.push(rec);
            continue;
        }
        let chrom_match = group.first().map_or(true, |first| {
            first.split('\t').next() == cols.first().copied()
        });
        if group.is_empty() || (chrom_match && pos == group_last_pos + 1) {
            group_last_pos = pos;
            group.push(rec);
        } else {
            let mut g = std::mem::take(&mut group);
            flush(&mut g, &mut out);
            group_last_pos = pos;
            group.push(rec);
        }
    }
    if !group.is_empty() {
        flush(&mut group, &mut out);
    }
    out
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_norm.rs"]
mod tests_join;

fn split_info_per_allele_keep_sum(info: &str, i: usize, n_alts: usize, keep_sum_keys: &[String]) -> String {
    if info == "." || info.is_empty() { return ".".into(); }
    let mut out = String::with_capacity(info.len());
    let mut first = true;
    for kv in info.split(';') {
        let (k, v) = match kv.split_once('=') { Some(p) => (p.0, Some(p.1)), None => (kv, None) };
        if let Some(val) = v {
            let vals: Vec<&str> = val.split(',').collect();
            let new_val = if keep_sum_keys.iter().any(|t| t == k) {
                let sum: f64 = vals.iter().filter_map(|s| s.parse::<f64>().ok()).sum();
                if vals.len() == n_alts {
                    if i == 0 { sum.to_string() } else { "0".to_string() }
                } else if vals.len() == n_alts + 1 {
                    if i == 0 { sum.to_string() } else { "0".to_string() }
                } else { val.to_string() }
            } else if vals.len() == n_alts {
                vals[i].to_string()
            } else if vals.len() == n_alts + 1 {
                format!("{},{}", vals[0], vals[i + 1])
            } else { val.to_string() };
            if !first { out.push(';'); }
            out.push_str(k); out.push('='); out.push_str(&new_val);
        } else {
            if !first { out.push(';'); }
            out.push_str(kv);
        }
        first = false;
    }
    out
}

fn split_info_per_allele(info: &str, i: usize, n_alts: usize) -> String {
    if info == "." || info.is_empty() { return ".".into(); }
    let mut out = String::with_capacity(info.len());
    let mut first = true;
    for kv in info.split(';') {
        let (k, v) = match kv.split_once('=') { Some(p) => (p.0, Some(p.1)), None => (kv, None) };
        if let Some(val) = v {
            let vals: Vec<&str> = val.split(',').collect();
            let new_val = if vals.len() == n_alts { vals[i].to_string() } else { val.to_string() };
            if !first { out.push(';'); }
            out.push_str(k); out.push('='); out.push_str(&new_val);
        } else {
            if !first { out.push(';'); }
            out.push_str(kv);
        }
        first = false;
    }
    out
}

fn split_samples_per_allele(samples: &[&str], allele_i: usize, n_alts: usize) -> Vec<String> {
    if samples.is_empty() { return Vec::new(); }
    let fmt_keys: Vec<&str> = samples[0].split(':').collect();
    let _ = fmt_keys; let _ = n_alts;
    let mut out = Vec::with_capacity(samples.len());
    out.push(samples[0].to_string());
    for s in &samples[1..] {
        let vals: Vec<&str> = s.split(':').collect();
        let new: Vec<String> = vals.iter().enumerate().map(|(j, v)| {
            if j == 0 { remap_gt(v, allele_i) } else { v.to_string() }
        }).collect();
        out.push(new.join(":"));
    }
    out
}

fn remap_gt(gt: &str, keep_allele: usize) -> String {
    let sep = if gt.contains('|') { '|' } else { '/' };
    let mapped: Vec<String> = gt.split(|c| c == '/' || c == '|').map(|a| {
        if a == "." || a.is_empty() { ".".into() }
        else if let Ok(n) = a.parse::<usize>() {
            if n == 0 { "0".into() } else if n == keep_allele + 1 { "1".into() } else { "0".into() }
        } else { a.to_string() }
    }).collect();
    mapped.join(&sep.to_string())
}

fn atomize_record(line: &str) -> Vec<String> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 { return vec![line.to_string()]; }
    let refa = cols[3].as_bytes(); let alt = cols[4];
    if alt.contains(',') { return vec![line.to_string()]; }
    let alta = alt.as_bytes();
    let pos: u32 = cols[1].parse().unwrap_or(0);
    if refa.len() < 2 || alta.len() < 2 || refa.len() != alta.len() { return vec![line.to_string()]; }

    let mut out = Vec::new();
    for (k, (&r, &a)) in refa.iter().zip(alta.iter()).enumerate() {
        if r == a { continue; }
        let mut c = cols.clone();
        let p = pos + k as u32;
        let r_str = (r as char).to_string();
        let a_str = (a as char).to_string();
        let p_str = p.to_string();
        c[1] = unsafe { &*(p_str.as_str() as *const str) };
        c[3] = unsafe { &*(r_str.as_str() as *const str) };
        c[4] = unsafe { &*(a_str.as_str() as *const str) };
        out.push(c.join("\t"));
    }
    if out.is_empty() { out.push(line.to_string()); }
    out
}

fn is_dup(window: &[(u32, Vec<String>)], pos: u32, line: &str, mode: RmDup) -> bool {
    if matches!(mode, RmDup::None) { return false; }
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 5 { return false; }
    for (p, lines) in window.iter().rev() {
        if *p != pos { continue; }
        for prev in lines {
            let pcols: Vec<&str> = prev.split('\t').collect();
            if pcols.len() < 5 { continue; }
            match mode {
                RmDup::Exact => if cols[3] == pcols[3] && cols[4] == pcols[4] { return true; },
                RmDup::All => return true,
                RmDup::AnyAllele => {
                    let new_alts: Vec<&str> = cols[4].split(',').collect();
                    let prev_alts: Vec<&str> = pcols[4].split(',').collect();
                    if new_alts.iter().any(|a| prev_alts.contains(a)) && cols[3] == pcols[3] { return true; }
                }
                RmDup::Snps => if cols[3].len() == 1 && cols[4].len() == 1 && pcols[3].len() == 1 && pcols[4].len() == 1 { return true; },
                RmDup::Indels => if cols[3].len() != cols[4].len() && pcols[3].len() != pcols[4].len() { return true; },
                RmDup::Both => if cols[3] == pcols[3] && cols[4] == pcols[4] { return true; },
                RmDup::AnyId => if cols[2] == pcols[2] && cols[2] != "." { return true; },
                RmDup::None => {}
            }
        }
    }
    false
}

fn push_dup(window: &mut Vec<(u32, Vec<String>)>, pos: u32, line: String) {
    if let Some((p, lines)) = window.last_mut() { if *p == pos { lines.push(line); return; } }
    window.push((pos, vec![line]));
    if window.len() > 64 { window.drain(..32); }
}

enum RefAction { Keep, Skip, Fix(String), Fail(String), Warn(String) }

struct Fasta { seqs: fxhash::FxHashMap<String, Vec<u8>> }
impl Fasta {
    fn has(&self, chr: &str) -> bool { self.seqs.contains_key(chr) }
    fn base(&self, chr: &str, pos: u32) -> Option<u8> {
        self.seqs.get(chr).and_then(|s| s.get((pos as usize).saturating_sub(1)).copied())
    }
    fn slice(&self, chr: &str, pos: u32, len: usize) -> Option<&[u8]> {
        let s = self.seqs.get(chr)?;
        let start = (pos as usize).saturating_sub(1);
        s.get(start..start + len)
    }
}

fn load_fasta(p: &Path) -> Result<Fasta> {
    let mut seqs: fxhash::FxHashMap<String, Vec<u8>> = fxhash::FxHashMap::default();
    let data = std::fs::read(p).with_context(|| format!("open fasta {:?}", p))?;
    let mut name: Option<String> = None;
    let mut cur: Vec<u8> = Vec::new();
    for line in data.split(|&b| b == b'\n') {
        if line.is_empty() { continue; }
        if line[0] == b'>' {
            if let Some(n) = name.take() { seqs.insert(n, std::mem::take(&mut cur)); }
            let rest = &line[1..];
            let end = rest.iter().position(|&b| b == b' ' || b == b'\t' || b == b'\r').unwrap_or(rest.len());
            name = Some(std::str::from_utf8(&rest[..end])?.to_string());
        } else {
            for &b in line { if b != b'\r' { cur.push(b.to_ascii_uppercase()); } }
        }
    }
    if let Some(n) = name { seqs.insert(n, cur); }
    Ok(Fasta { seqs })
}

fn verify_ref(line: &str, fa: &Fasta, mode: CheckRef) -> RefAction {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 5 { return RefAction::Keep; }
    let pos: u32 = match cols[1].parse() { Ok(v) => v, Err(_) => return RefAction::Keep };
    let refa = cols[3].as_bytes();
    let Some(seq) = fa.slice(cols[0], pos, refa.len()) else {
        let msg = format!("REF mismatch at {}:{} (no fasta for chromosome)", cols[0], pos);
        return match mode { CheckRef::Exit => RefAction::Fail(msg), CheckRef::Warn => RefAction::Warn(msg), CheckRef::Exclude => RefAction::Skip, CheckRef::Set => RefAction::Keep };
    };
    let ru: Vec<u8> = refa.iter().map(|b| b.to_ascii_uppercase()).collect();
    if seq == ru.as_slice() { return RefAction::Keep; }
    let new_ref = std::str::from_utf8(seq).unwrap_or("N").to_string();
    let msg = format!("REF mismatch at {}:{} (vcf={}, fasta={})", cols[0], pos, cols[3], new_ref);
    match mode {
        CheckRef::Exit => RefAction::Fail(msg),
        CheckRef::Warn => RefAction::Warn(msg),
        CheckRef::Exclude => RefAction::Skip,
        CheckRef::Set => RefAction::Fix(new_ref),
    }
}

fn is_symbolic_alt(a: &str) -> bool {
    a.starts_with('<') || a == "*" || a == "."
}

/// Largest symbolic SV span left-alignable without materialising the whole
/// reference allele. Beyond this the record passes through unshifted.
const MAX_SV_SPAN: u32 = 1_000_000;

fn left_align(line: &str, fa: &Fasta) -> Option<String> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 { return None; }
    let chr = cols[0];
    let mut pos: u32 = cols[1].parse().ok()?;
    let alt_raw: Vec<&str> = cols[4].split(',').collect();
    if alt_raw.is_empty() {
        return None;
    }

    // Symbolic alleles (<DEL>, <DUP>, ...) carry their span in INFO/END and are
    // realigned alongside any sequence alleles, shifting POS and END together.
    if alt_raw.iter().any(|a| is_symbolic_alt(a)) {
        return left_align_symbolic(&cols, chr, pos, &alt_raw, fa);
    }

    let mut r: Vec<u8> = cols[3].as_bytes().to_ascii_uppercase();
    let mut alts: Vec<Vec<u8>> = alt_raw.iter().map(|a| a.as_bytes().to_ascii_uppercase()).collect();
    if r.is_empty() || alts.iter().any(|a| a.is_empty()) {
        return None;
    }

    left_shift(&mut r, &mut alts, &mut pos, chr, fa)?;
    // Left-trim a shared leading base while every allele keeps length >= 2 (drop the redundant
    // anchor on complex/MNP-style records).
    while r.len() >= 2 && alts.iter().all(|a| a.len() >= 2 && a[0] == r[0]) {
        r.remove(0);
        for a in alts.iter_mut() { a.remove(0); }
        pos += 1;
    }

    let new_ref = String::from_utf8(r).ok()?;
    let new_alt = alts
        .into_iter()
        .map(|a| String::from_utf8(a).unwrap_or_default())
        .collect::<Vec<_>>()
        .join(",");
    let mut out: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
    out[1] = pos.to_string();
    out[3] = new_ref;
    out[4] = new_alt;
    Some(out.join("\t"))
}

/// Canonical left-alignment loop (Tan/Abecasis/Kang 2015; bcftools/vt `norm`):
/// while ref and all alts share a trailing base, truncate it; if that would
/// empty an allele, extend one reference base to the left (shifting `pos`) so
/// the indel slides to its leftmost equivalent position.
fn left_shift(r: &mut Vec<u8>, alts: &mut [Vec<u8>], pos: &mut u32, chr: &str, fa: &Fasta) -> Option<()> {
    loop {
        let last = *r.last().unwrap();
        if !alts.iter().all(|a| *a.last().unwrap() == last) {
            break;
        }
        if r.len() == 1 || alts.iter().any(|a| a.len() == 1) {
            if *pos <= 1 { break; }
            *pos -= 1;
            let prev = fa.base(chr, *pos)?.to_ascii_uppercase();
            r.insert(0, prev);
            for a in alts.iter_mut() { a.insert(0, prev); }
        } else {
            r.pop();
            for a in alts.iter_mut() { a.pop(); }
        }
    }
    Some(())
}

/// Left-align a record carrying at least one symbolic ALT. The shift is driven
/// by the sequence representation (either the explicit REF + sequence ALTs, or
/// — for a pure symbolic deletion — the reference span POS..END), then POS, END
/// and the REF anchor are moved left by the same amount. Symbolic ALTs are
/// preserved verbatim.
fn left_align_symbolic(
    cols: &[&str],
    chr: &str,
    orig_pos: u32,
    alt_raw: &[&str],
    fa: &Fasta,
) -> Option<String> {
    let end = info_get_u32(cols[7], "END")?;
    if end < orig_pos || end - orig_pos > MAX_SV_SPAN {
        return None;
    }
    let seq_idx: Vec<usize> = alt_raw
        .iter()
        .enumerate()
        .filter(|(_, a)| !is_symbolic_alt(a))
        .map(|(i, _)| i)
        .collect();

    let mut pos = orig_pos;
    let symbolic_only = seq_idx.is_empty();
    let (mut r, mut seq_alts): (Vec<u8>, Vec<Vec<u8>>) = if symbolic_only {
        // Implied sequence form of a symbolic deletion: REF = ref[POS..=END],
        // ALT = the single anchor base.
        let mut refseq = Vec::with_capacity((end - orig_pos + 1) as usize);
        for p in orig_pos..=end {
            refseq.push(fa.base(chr, p)?.to_ascii_uppercase());
        }
        let anchor = vec![refseq[0]];
        (refseq, vec![anchor])
    } else {
        let r = cols[3].as_bytes().to_ascii_uppercase();
        let alts = seq_idx
            .iter()
            .map(|&i| alt_raw[i].as_bytes().to_ascii_uppercase())
            .collect();
        (r, alts)
    };
    if r.is_empty() || seq_alts.iter().any(|a| a.is_empty()) {
        return None;
    }

    left_shift(&mut r, &mut seq_alts, &mut pos, chr, fa)?;

    let delta = orig_pos - pos;
    if delta == 0 {
        return None; // nothing moved; leave the record untouched
    }
    let new_end = end - delta;

    let new_ref = if symbolic_only {
        (r[0] as char).to_string()
    } else {
        String::from_utf8(r).ok()?
    };
    let mut seq_iter = seq_alts.into_iter();
    let new_alt = alt_raw
        .iter()
        .map(|a| {
            if is_symbolic_alt(a) {
                a.to_string()
            } else {
                seq_iter
                    .next()
                    .map(|v| String::from_utf8(v).unwrap_or_default())
                    .unwrap_or_default()
            }
        })
        .collect::<Vec<_>>()
        .join(",");

    let mut out: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
    out[1] = pos.to_string();
    out[3] = new_ref;
    out[4] = new_alt;
    out[7] = info_set_u32(cols[7], "END", new_end);
    Some(out.join("\t"))
}

fn info_get_u32(info: &str, key: &str) -> Option<u32> {
    for f in info.split(';') {
        if let Some((k, v)) = f.split_once('=') {
            if k == key {
                return v.parse().ok();
            }
        }
    }
    None
}

fn info_set_u32(info: &str, key: &str, val: u32) -> String {
    info.split(';')
        .map(|f| match f.split_once('=') {
            Some((k, _)) if k == key => format!("{key}={val}"),
            _ => f.to_string(),
        })
        .collect::<Vec<_>>()
        .join(";")
}
