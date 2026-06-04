use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{Cursor, Read, Write};
use std::path::{Path, PathBuf};

use crate::cli::args::ReheaderArgs;

pub fn cmd_reheader(args: ReheaderArgs) -> Result<()> {
    let mut argv: Vec<String> = Vec::new();
    if let Some(p) = &args.header { argv.push("-h".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(p) = &args.samples_file { argv.push("-s".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(p) = &args.fai { argv.push("-f".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(p) = &args.output { argv.push("-o".into()); argv.push(p.to_string_lossy().into_owned()); }
    argv.extend(args.passthrough.iter().cloned());
    let cfg = parse_reheader_args(&argv)?;
    let (input_text, input_is_bcf) = read_input_text(&args.input)?;
    let (orig_meta, orig_chrom, data_lines) = split_vcf_text(&input_text);

    let (mut meta, mut chrom_line) = if let Some(hdr_path) = &cfg.header_path {
        let hdr_text = fs::read_to_string(hdr_path)?;
        let (h_meta, h_chrom, _) = split_vcf_text(&hdr_text);
        if input_is_bcf {
            merge_headers_for_bcf(&orig_meta, &h_meta, &h_chrom)
        } else {
            (h_meta, h_chrom)
        }
    } else {
        (orig_meta, orig_chrom)
    };

    if let Some(samples_path) = &cfg.samples_path {
        chrom_line = apply_sample_renames(&chrom_line, samples_path)?;
    }

    if let Some(fai_path) = &cfg.fai_path {
        let fai = read_fai(fai_path)?;
        apply_fai_to_contigs(&mut meta, &fai);
    }

    let mut out = String::new();
    for line in &meta {
        out.push_str(line);
        out.push('\n');
    }
    out.push_str(&chrom_line);
    out.push('\n');
    for line in &data_lines {
        out.push_str(line);
        out.push('\n');
    }

    if let Some(path) = cfg.output_path {
        let mut f = fs::File::create(path)?;
        f.write_all(out.as_bytes())?;
    } else {
        print!("{out}");
    }
    Ok(())
}

#[derive(Default)]
struct ReheaderCfg {
    header_path: Option<PathBuf>,
    samples_path: Option<PathBuf>,
    fai_path: Option<PathBuf>,
    output_path: Option<PathBuf>,
}

fn parse_reheader_args(args: &[String]) -> Result<ReheaderCfg> {
    let mut cfg = ReheaderCfg::default();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-h" => {
                i += 1;
                cfg.header_path = Some(PathBuf::from(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -h"))?,
                ));
            }
            "-s" => {
                i += 1;
                cfg.samples_path = Some(PathBuf::from(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -s"))?,
                ));
            }
            "-f" => {
                i += 1;
                cfg.fai_path = Some(PathBuf::from(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -f"))?,
                ));
            }
            "-o" => {
                i += 1;
                cfg.output_path = Some(PathBuf::from(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -o"))?,
                ));
            }
            _ => {}
        }
        i += 1;
    }
    Ok(cfg)
}

fn read_input_text(input: &Path) -> Result<(String, bool)> {
    if input == Path::new("-") {
        let mut bytes = Vec::new();
        std::io::stdin().read_to_end(&mut bytes)?;
        return Ok((decode_input_bytes(&bytes, "")?, false));
    }

    let is_bcf_ext = input.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("bcf"))
        .unwrap_or(false);

    if is_bcf_ext || is_bcf_magic(input)? {
        let mut sibling = input.to_path_buf();
        sibling.set_extension("vcf");
        if sibling.exists() {
            return Ok((fs::read_to_string(sibling)?, true));
        }
        return decode_bcf_to_vcf(input).map(|s| (s, true));
    }

    let bytes = fs::read(input)?;
    Ok((decode_input_bytes(&bytes, &input.to_string_lossy())?, false))
}

fn is_bcf_magic(input: &Path) -> Result<bool> {
    use std::io::Read as _;
    let mut f = fs::File::open(input)?;
    let mut probe = [0u8; 5];
    let n = f.read(&mut probe)?;
    if n < 5 { return Ok(false); }
    if probe == *crate::bcf::BCF_MAGIC { return Ok(true); }
    if probe[0] == 0x1F && probe[1] == 0x8B {
        let f2 = fs::File::open(input)?;
        let mut rd = noodles_bgzf::io::Reader::new(f2);
        let mut probe2 = [0u8; 5];
        if rd.read_exact(&mut probe2).is_ok() && probe2 == *crate::bcf::BCF_MAGIC {
            return Ok(true);
        }
    }
    Ok(false)
}

fn decode_bcf_to_vcf(input: &Path) -> Result<String> {
    let mut r = crate::bcf::BcfReader::open(input)?;
    let mut out = String::new();
    for h in &r.header_lines {
        out.push_str(h);
        out.push('\n');
    }
    while let Some(line) = r.read_record_line()? {
        out.push_str(&line);
        if !line.ends_with('\n') { out.push('\n'); }
    }
    Ok(out)
}

fn decode_input_bytes(bytes: &[u8], path_hint: &str) -> Result<String> {
    let is_gz_magic = bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b;
    let is_gz_ext = path_hint.ends_with(".gz") || path_hint.ends_with(".bgz");

    if is_gz_magic || is_gz_ext {
        if bytes.len() >= 7 && bytes[0] == 0x1F && bytes[1] == 0x8B {
            let mut probe = MultiGzDecoder::new(Cursor::new(bytes));
            let mut head = [0u8; 5];
            if probe.read_exact(&mut head).is_ok() && head == *crate::bcf::BCF_MAGIC {
                let tmp = std::env::temp_dir().join(format!("kira-bt-reh-{}.bcf", std::process::id()));
                fs::write(&tmp, bytes)?;
                let s = decode_bcf_to_vcf(&tmp);
                let _ = fs::remove_file(&tmp);
                return s;
            }
        }
        let mut s = String::new();
        MultiGzDecoder::new(Cursor::new(bytes)).read_to_string(&mut s)?;
        return Ok(s);
    }

    if bytes.len() >= 5 && &bytes[..5] == crate::bcf::BCF_MAGIC {
        let tmp = std::env::temp_dir().join(format!("kira-bt-reh-{}.bcf", std::process::id()));
        fs::write(&tmp, bytes)?;
        let s = decode_bcf_to_vcf(&tmp);
        let _ = fs::remove_file(&tmp);
        return s;
    }

    String::from_utf8(bytes.to_vec())
        .map_err(|_| anyhow::anyhow!("input is binary but not BCF; cannot decode"))
}

fn split_vcf_text(text: &str) -> (Vec<String>, String, Vec<String>) {
    let mut meta = Vec::new();
    let mut chrom = String::new();
    let mut data = Vec::new();

    for line in text.lines() {
        if line.starts_with("##") {
            meta.push(line.to_string());
        } else if line.starts_with('#') {
            chrom = line.to_string();
        } else if !line.trim().is_empty() {
            data.push(line.to_string());
        }
    }

    (meta, chrom, data)
}

fn parse_samples_map_line(line: &str) -> Option<(String, String)> {
    let t = line.trim();
    if t.is_empty() {
        return None;
    }
    let mut split_at = None;
    let mut escaped = false;
    for (i, c) in t.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        if c == '\\' {
            escaped = true;
            continue;
        }
        if c.is_whitespace() {
            split_at = Some(i);
            break;
        }
    }

    if let Some(i) = split_at {
        let (lhs, rhs0) = t.split_at(i);
        let rhs = rhs0.trim_start();
        if rhs.is_empty() {
            let v = unescape_backslash(lhs);
            Some((String::new(), v))
        } else {
            Some((unescape_backslash(lhs), unescape_backslash(rhs)))
        }
    } else {
        Some((String::new(), unescape_backslash(t)))
    }
}

fn unescape_backslash(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut esc = false;
    for c in s.chars() {
        if esc {
            out.push(c);
            esc = false;
        } else if c == '\\' {
            esc = true;
        } else {
            out.push(c);
        }
    }
    if esc {
        out.push('\\');
    }
    out
}

fn apply_sample_renames(chrom_line: &str, samples_path: &Path) -> Result<String> {
    let map_text = fs::read_to_string(samples_path)?;
    let mut entries = Vec::<(String, String)>::new();
    for line in map_text.lines() {
        if let Some(v) = parse_samples_map_line(line) {
            entries.push(v);
        }
    }

    let mut cols: Vec<String> = chrom_line.split('\t').map(|s| s.to_string()).collect();
    if cols.len() <= 9 {
        return Ok(chrom_line.to_string());
    }

    let is_pair_map = entries.iter().any(|(k, _)| !k.is_empty());
    if is_pair_map {
        let mut m = HashMap::<String, String>::new();
        for (k, v) in entries {
            if !k.is_empty() {
                m.insert(k, v);
            }
        }
        for sample in cols.iter_mut().skip(9) {
            if let Some(new_name) = m.get(sample) {
                *sample = new_name.clone();
            }
        }
    } else {
        for (idx, (_, v)) in entries.iter().enumerate() {
            let col = 9 + idx;
            if col < cols.len() {
                cols[col] = v.clone();
            }
        }
    }

    Ok(cols.join("\t"))
}

fn read_fai(path: &Path) -> Result<Vec<(String, u64)>> {
    let mut out = Vec::<(String, u64)>::new();
    let s = fs::read_to_string(path)?;
    for line in s.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.split('\t');
        if let (Some(name), Some(len)) = (it.next(), it.next()) {
            if let Ok(v) = len.parse::<u64>() {
                out.push((name.to_string(), v));
            }
        }
    }
    Ok(out)
}

fn apply_fai_to_contigs(meta: &mut Vec<String>, fai: &[(String, u64)]) {
    let fai_map: HashMap<String, u64> = fai.iter().map(|(k, v)| (k.clone(), *v)).collect();
    let mut seen_ids = HashSet::<String>::new();
    let mut out = Vec::<String>::with_capacity(meta.len() + fai.len());

    for line in meta.iter() {
        if !line.starts_with("##contig=<") {
            out.push(line.clone());
            continue;
        }

        let Some(id) = extract_contig_id(line) else {
            continue;
        };
        let Some(new_len) = fai_map.get(&id) else {
            continue;
        };

        let updated = upsert_contig_length(line, *new_len);
        out.push(updated);
        seen_ids.insert(id);
    }

    for (id, len) in fai {
        if !seen_ids.contains(id) {
            out.push(format!("##contig=<ID={id},length={len}>"));
        }
    }

    *meta = out;
}

fn extract_contig_id(line: &str) -> Option<String> {
    let body = line.strip_prefix("##contig=<")?.strip_suffix('>')?;
    for part in body.split(',') {
        if let Some(v) = part.strip_prefix("ID=") {
            return Some(v.to_string());
        }
    }
    None
}

fn upsert_contig_length(line: &str, len: u64) -> String {
    let Some(body) = line
        .strip_prefix("##contig=<")
        .and_then(|s| s.strip_suffix('>'))
    else {
        return line.to_string();
    };
    let fields = split_header_fields(body);
    if fields.is_empty() {
        return line.to_string();
    }

    let mut id_val = None::<String>;
    let mut rest = Vec::<String>::new();
    for f in fields {
        if let Some(v) = f.strip_prefix("ID=") {
            id_val = Some(v.to_string());
        } else if !f.starts_with("length=") {
            rest.push(f);
        }
    }
    let Some(id) = id_val else {
        return line.to_string();
    };
    let mut out = format!("##contig=<ID={id},length={len}");
    for f in rest {
        out.push(',');
        out.push_str(&f);
    }
    out.push('>');
    out
}

fn split_header_fields(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut escaped = false;
    for ch in s.chars() {
        if escaped {
            cur.push(ch);
            escaped = false;
            continue;
        }
        if ch == '\\' {
            cur.push(ch);
            escaped = true;
            continue;
        }
        if ch == '"' {
            in_quotes = !in_quotes;
            cur.push(ch);
            continue;
        }
        if ch == ',' && !in_quotes {
            out.push(cur);
            cur = String::new();
            continue;
        }
        cur.push(ch);
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

fn extract_meta_key(line: &str) -> Option<String> {
    if line.starts_with("##INFO=<")
        || line.starts_with("##FORMAT=<")
        || line.starts_with("##FILTER=<")
        || line.starts_with("##contig=<")
    {
        let id = extract_contig_id_like(line)?;
        if line.starts_with("##INFO=<") {
            return Some(format!("INFO:{id}"));
        }
        if line.starts_with("##FORMAT=<") {
            return Some(format!("FORMAT:{id}"));
        }
        if line.starts_with("##FILTER=<") {
            return Some(format!("FILTER:{id}"));
        }
        return Some(format!("CONTIG:{id}"));
    }
    None
}

fn extract_contig_id_like(line: &str) -> Option<String> {
    let body = line.split_once('<')?.1.strip_suffix('>')?;
    for p in split_header_fields(body) {
        if let Some(v) = p.strip_prefix("ID=") {
            return Some(v.to_string());
        }
    }
    None
}

fn merge_headers_for_bcf(
    old_meta: &[String],
    new_meta: &[String],
    new_chrom: &str,
) -> (Vec<String>, String) {
    let mut new_by_key = HashMap::<String, String>::new();
    let mut new_fileformat = None::<String>;
    let mut new_non_key = HashSet::<String>::new();
    for line in new_meta {
        if line.starts_with("##fileformat=") {
            new_fileformat = Some(line.clone());
            continue;
        }
        if let Some(k) = extract_meta_key(line) {
            new_by_key.insert(k, line.clone());
        } else {
            new_non_key.insert(line.clone());
        }
    }

    let mut out = Vec::<String>::new();
    let mut seen = HashSet::<String>::new();
    let mut fileformat_done = false;

    for line in old_meta {
        if line.starts_with("##fileformat=") {
            if let Some(ff) = &new_fileformat {
                out.push(ff.clone());
            } else {
                out.push(line.clone());
            }
            fileformat_done = true;
            continue;
        }
        if let Some(k) = extract_meta_key(line) {
            if let Some(repl) = new_by_key.get(&k) {
                out.push(repl.clone());
                seen.insert(k);
            }
        } else {
            if new_non_key.contains(line) {
                out.push(line.clone());
            }
        }
    }

    if !fileformat_done {
        if let Some(ff) = new_fileformat {
            out.insert(0, ff);
        }
    }

    for line in new_meta {
        if let Some(k) = extract_meta_key(line) {
            if !seen.contains(&k) {
                out.push(line.clone());
            }
        } else if !out.iter().any(|l| l == line) && !line.starts_with("##fileformat=") {
            out.push(line.clone());
        }
    }

    (out, new_chrom.to_string())
}
