use anyhow::{Context, Result, anyhow, bail};
use fxhash::FxHashMap;
use memchr::memchr;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Default)]
pub struct PostProcessor {
    pub remove: Option<RemoveSpec>,
    pub set_id: Option<IdTemplate>,
    pub mark_sites: Option<MarkSpec>,
    pub include: Option<Predicate>,
    pub exclude: Option<Predicate>,
    pub keep_sites: bool,
    pub rename_chrs: Option<FxHashMap<String, String>>,
    pub rename_annots: Option<RenameAnnots>,
    pub samples_keep: Option<Vec<usize>>,
    pub no_version: bool,
    pub extra_header_lines: Vec<String>,
    pub match_tags: Vec<String>,
}

pub struct RenameAnnots {
    pub info: FxHashMap<String, String>,
    pub format: FxHashMap<String, String>,
    pub filter: FxHashMap<String, String>,
}

pub struct RemoveSpec {
    pub drop_id: bool,
    pub drop_qual: bool,
    pub drop_filter: bool,
    pub drop_all_info: bool,
    pub drop_all_format: bool,
    pub info_tags: Vec<String>,
    pub format_tags: Vec<String>,
    pub inverse: bool,
}

pub enum MarkOp { Present, Absent }

pub struct MarkSpec {
    pub tag: String,
    pub op: MarkOp,
}

pub struct IdTemplate {
    pub tokens: Vec<IdTok>,
    pub fill_only_missing: bool,
}

pub enum IdTok {
    Lit(String),
    Chrom,
    Pos,
    Ref,
    Alt,
    Id,
    InfoTag(String),
}

pub struct Predicate {
    pub engine: crate::filter::FilterEngine,
}

impl PostProcessor {
    pub fn parse_remove(spec: &str) -> Result<RemoveSpec> {
        let (inverse, body) = match spec.strip_prefix('^') {
            Some(rest) => (true, rest),
            None => (false, spec),
        };
        let mut out = RemoveSpec {
            drop_id: false, drop_qual: false, drop_filter: false,
            drop_all_info: false, drop_all_format: false,
            info_tags: Vec::new(), format_tags: Vec::new(), inverse,
        };
        for raw in body.split(',') {
            let item = raw.trim();
            if item.is_empty() { continue; }
            match item {
                "ID" => out.drop_id = true,
                "QUAL" => out.drop_qual = true,
                "FILTER" => out.drop_filter = true,
                "INFO" => out.drop_all_info = true,
                "FORMAT" | "FMT" => out.drop_all_format = true,
                _ => {
                    if let Some(tag) = item.strip_prefix("INFO/") {
                        out.info_tags.push(tag.to_string());
                    } else if let Some(tag) = item.strip_prefix("FORMAT/").or_else(|| item.strip_prefix("FMT/")) {
                        out.format_tags.push(tag.to_string());
                    } else {
                        out.info_tags.push(item.to_string());
                    }
                }
            }
        }
        Ok(out)
    }

    pub fn parse_set_id(spec: &str) -> Result<IdTemplate> {
        let (fill_only_missing, body) = match spec.strip_prefix('+') {
            Some(rest) => (true, rest),
            None => (false, spec),
        };
        let mut tokens = Vec::new();
        let bytes = body.as_bytes();
        let mut i = 0;
        let mut lit = String::new();
        while i < bytes.len() {
            if bytes[i] == b'%' && i + 1 < bytes.len() {
                let rest = &body[i + 1..];
                let (name, advance) = parse_id_token(rest);
                if !name.is_empty() {
                    if !lit.is_empty() {
                        tokens.push(IdTok::Lit(std::mem::take(&mut lit)));
                    }
                    let tok = match name.as_str() {
                        "CHROM" => IdTok::Chrom,
                        "POS" => IdTok::Pos,
                        "REF" => IdTok::Ref,
                        "ALT" => IdTok::Alt,
                        "ID" => IdTok::Id,
                        other => {
                            let tag = other.strip_prefix("INFO/").unwrap_or(other).to_string();
                            IdTok::InfoTag(tag)
                        }
                    };
                    tokens.push(tok);
                    i += 1 + advance;
                    continue;
                }
            }
            lit.push(bytes[i] as char);
            i += 1;
        }
        if !lit.is_empty() { tokens.push(IdTok::Lit(lit)); }
        Ok(IdTemplate { tokens, fill_only_missing })
    }

    pub fn parse_mark_sites(spec: &str) -> Result<MarkSpec> {
        let (op, tag) = match spec.as_bytes().first() {
            Some(b'+') => (MarkOp::Present, &spec[1..]),
            Some(b'-') => (MarkOp::Absent, &spec[1..]),
            _ => (MarkOp::Present, spec),
        };
        if tag.is_empty() { bail!("--mark-sites: empty tag"); }
        Ok(MarkSpec { tag: tag.to_string(), op })
    }

    pub fn read_rename_chrs<P: AsRef<Path>>(p: P) -> Result<FxHashMap<String, String>> {
        let mut m = FxHashMap::default();
        for line in BufReader::new(File::open(p.as_ref())?).lines() {
            let line = line?;
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            let mut it = t.split('\t');
            let (Some(a), Some(b)) = (it.next(), it.next()) else { continue };
            m.insert(a.to_string(), b.to_string());
        }
        Ok(m)
    }

    pub fn read_rename_annots<P: AsRef<Path>>(p: P) -> Result<RenameAnnots> {
        let mut r = RenameAnnots {
            info: FxHashMap::default(),
            format: FxHashMap::default(),
            filter: FxHashMap::default(),
        };
        for line in BufReader::new(File::open(p.as_ref())?).lines() {
            let line = line?;
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            let mut it = t.split('\t');
            let (Some(lhs), Some(new)) = (it.next(), it.next()) else { continue };
            let (kind, old) = lhs.split_once('/').ok_or_else(|| anyhow!("rename-annots: expect TYPE/OLD"))?;
            match kind {
                "INFO" => { r.info.insert(old.to_string(), new.to_string()); }
                "FORMAT" | "FMT" => { r.format.insert(old.to_string(), new.to_string()); }
                "FILTER" => { r.filter.insert(old.to_string(), new.to_string()); }
                _ => bail!("rename-annots: unknown TYPE '{}'", kind),
            }
        }
        Ok(r)
    }
}

fn parse_id_token(s: &str) -> (String, usize) {
    let bytes = s.as_bytes();
    let mut end = 0;
    while end < bytes.len() && bytes[end].is_ascii_uppercase() {
        end += 1;
    }
    if end == 0 { return (String::new(), 0); }
    if bytes.get(end) == Some(&b'/')
        && bytes.get(end + 1).map_or(false, |b| b.is_ascii_alphanumeric() || *b == b'_')
    {
        let mut p = end + 1;
        while p < bytes.len() && (bytes[p].is_ascii_alphanumeric() || bytes[p] == b'_') {
            p += 1;
        }
        return (s[..p].to_string(), p);
    }
    (s[..end].to_string(), end)
}

pub struct HeaderOptions<'a> {
    pub no_version: bool,
    pub extra_header_lines: &'a [String],
    pub remove: Option<&'a RemoveSpec>,
    pub rename_chrs: Option<&'a FxHashMap<String, String>>,
    pub rename_annots: Option<&'a RenameAnnots>,
    pub mark_sites: Option<&'a MarkSpec>,
    pub set_id: bool,
    pub samples_keep: Option<&'a [usize]>,
    pub version_line: Option<&'a str>,
}

pub fn apply_to_header(headers: Vec<String>, opts: &HeaderOptions<'_>) -> Vec<String> {
    let mut out = Vec::with_capacity(headers.len() + opts.extra_header_lines.len() + 4);
    let mut chrom_line: Option<String> = None;
    for h in headers {
        if h.starts_with("#CHROM") {
            chrom_line = Some(h);
            continue;
        }
        if let Some(rest) = h.strip_prefix("##contig=") {
            if let Some(rmap) = opts.rename_chrs {
                out.push(format!("##contig={}", rewrite_contig_id(rest, rmap)));
                continue;
            }
        }
        if let Some(line) = filter_header_line(&h, opts) {
            out.push(line);
        }
    }
    for extra in opts.extra_header_lines {
        if !out.iter().any(|h| h == extra) { out.push(extra.clone()); }
    }
    if !opts.no_version {
        if let Some(v) = opts.version_line {
            out.push(v.to_string());
        }
    }
    if let Some(m) = opts.mark_sites {
        let line = format!(
            "##INFO=<ID={},Number=0,Type=Flag,Description=\"Site is {} in -a file\">",
            m.tag,
            match m.op { MarkOp::Present => "present", MarkOp::Absent => "absent" },
        );
        if !out.iter().any(|h| h.starts_with(&format!("##INFO=<ID={},", m.tag))) {
            out.push(line);
        }
    }
    if let Some(line) = chrom_line {
        if let Some(keep) = opts.samples_keep {
            out.push(filter_chrom_samples(&line, keep));
        } else {
            out.push(line);
        }
    }
    out
}

fn rewrite_contig_id(rest: &str, rmap: &FxHashMap<String, String>) -> String {
    if let Some(body) = rest.strip_prefix('<').and_then(|s| s.strip_suffix('>')) {
        let mut parts: Vec<String> = Vec::new();
        for kv in split_top_commas(body) {
            if let Some(rest) = kv.strip_prefix("ID=") {
                if let Some(new) = rmap.get(rest) {
                    parts.push(format!("ID={}", new));
                    continue;
                }
            }
            parts.push(kv.to_string());
        }
        format!("<{}>", parts.join(","))
    } else {
        rest.to_string()
    }
}

fn split_top_commas(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let bytes = s.as_bytes();
    let mut start = 0;
    let mut depth = 0i32;
    let mut in_quote = false;
    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'"' => in_quote = !in_quote,
            b'<' if !in_quote => depth += 1,
            b'>' if !in_quote => depth -= 1,
            b',' if !in_quote && depth == 0 => {
                out.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    out.push(&s[start..]);
    out
}

fn filter_header_line(h: &str, opts: &HeaderOptions<'_>) -> Option<String> {
    let kind_id = parse_struct_header(h);
    if let Some((kind, id)) = kind_id {
        if let Some(rm) = opts.remove {
            if remove_drops_header(rm, kind, id) { return None; }
        }
        if let Some(ren) = opts.rename_annots {
            let map = match kind { "INFO" => &ren.info, "FORMAT" => &ren.format, "FILTER" => &ren.filter, _ => return Some(h.to_string()) };
            if let Some(new) = map.get(id) {
                return Some(rewrite_struct_id(h, new));
            }
        }
    }
    Some(h.to_string())
}

fn parse_struct_header(h: &str) -> Option<(&str, &str)> {
    let kind = if let Some(rest) = h.strip_prefix("##INFO=<") { ("INFO", rest) }
        else if let Some(rest) = h.strip_prefix("##FORMAT=<") { ("FORMAT", rest) }
        else if let Some(rest) = h.strip_prefix("##FILTER=<") { ("FILTER", rest) }
        else { return None; };
    let body = kind.1.strip_suffix('>')?;
    let id_part = body.strip_prefix("ID=")?;
    let end = id_part.find(',').unwrap_or(id_part.len());
    Some((kind.0, &id_part[..end]))
}

fn rewrite_struct_id(h: &str, new_id: &str) -> String {
    let (prefix_end, _) = h.match_indices("ID=").next().expect("must contain ID=");
    let after_id = prefix_end + 3;
    let tail_start = h[after_id..].find(',').map(|p| after_id + p).unwrap_or_else(|| h.len() - 1);
    let mut s = String::with_capacity(h.len() + new_id.len());
    s.push_str(&h[..after_id]);
    s.push_str(new_id);
    s.push_str(&h[tail_start..]);
    s
}

fn remove_drops_header(rm: &RemoveSpec, kind: &str, id: &str) -> bool {
    let listed = match kind {
        "INFO" => rm.drop_all_info || rm.info_tags.iter().any(|t| t == id),
        "FORMAT" => rm.drop_all_format || rm.format_tags.iter().any(|t| t == id),
        "FILTER" => rm.drop_filter,
        _ => false,
    };
    if rm.inverse { !listed } else { listed }
}

fn filter_chrom_samples(line: &str, keep: &[usize]) -> String {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() <= 9 { return line.to_string(); }
    let mut out = String::with_capacity(line.len());
    out.push_str(&cols[..9].join("\t"));
    for &i in keep {
        let sample_idx = 9 + i;
        if sample_idx < cols.len() {
            out.push('\t');
            out.push_str(cols[sample_idx]);
        }
    }
    out
}

pub enum LineAction {
    Keep,
    Drop,
    Replace(String),
}

pub fn process_record_line(
    line: &str,
    pp: &PostProcessor,
    matched_db: bool,
) -> LineAction {
    if line.is_empty() || line.as_bytes()[0] == b'#' {
        return LineAction::Keep;
    }

    let cols_borrow: Vec<&str> = line.split('\t').collect();
    if cols_borrow.len() < 8 { return LineAction::Keep; }

    let chrom_renamed: Option<String> = pp.rename_chrs.as_ref().and_then(|m| m.get(cols_borrow[0]).cloned());
    let chrom: &str = chrom_renamed.as_deref().unwrap_or(cols_borrow[0]);

    if let Some(p) = pp.include.as_ref().or(pp.exclude.as_ref()) {
        let include_mode = pp.include.is_some();
        let pass = eval_predicate(p, line);
        let keep = if include_mode { pass } else { !pass };
        if !keep {
            return if pp.keep_sites { LineAction::Keep } else { LineAction::Drop };
        }
    }

    let mut modified = chrom_renamed.is_some();
    let mut id = cols_borrow[2].to_string();
    let mut qual = cols_borrow[5].to_string();
    let mut filter = cols_borrow[6].to_string();
    let mut info = cols_borrow[7].to_string();
    let mut format = if cols_borrow.len() > 8 { cols_borrow[8].to_string() } else { String::new() };
    let mut samples: Vec<String> = if cols_borrow.len() > 9 {
        cols_borrow[9..].iter().map(|s| s.to_string()).collect()
    } else { Vec::new() };

    if let Some(rm) = &pp.remove {
        if rm.drop_id && !rm.inverse { id = ".".into(); modified = true; }
        if rm.drop_qual && !rm.inverse { qual = ".".into(); modified = true; }
        if rm.drop_filter && !rm.inverse { filter = ".".into(); modified = true; }
        if rm.drop_all_info && !rm.inverse {
            info = ".".into(); modified = true;
        } else if !rm.info_tags.is_empty() || rm.inverse {
            let new = filter_info(&info, &rm.info_tags, rm.inverse, rm.drop_all_info);
            if new != info { info = new; modified = true; }
        }
        if !format.is_empty() && (!rm.format_tags.is_empty() || rm.drop_all_format || rm.inverse) {
            let (new_format, removed_idx) = filter_format_keys(&format, &rm.format_tags, rm.inverse, rm.drop_all_format);
            if !removed_idx.is_empty() || new_format != format {
                if new_format.is_empty() {
                    samples.clear();
                    format.clear();
                } else {
                    for s in samples.iter_mut() {
                        *s = strip_format_indices(s, &removed_idx);
                    }
                    format = new_format;
                }
                modified = true;
            }
        }
    }

    if let Some(t) = &pp.set_id {
        if !t.fill_only_missing || id == "." || id.is_empty() {
            let new = render_id_template(t, chrom, cols_borrow[1], cols_borrow[3], cols_borrow[4], &id, &info);
            if new != id { id = new; modified = true; }
        }
    }

    if let Some(m) = &pp.mark_sites {
        let detected = if !pp.match_tags.is_empty() && !(info == "." || info.is_empty()) {
            info.split(';').any(|kv| {
                let k = match memchr(b'=', kv.as_bytes()) { Some(p) => &kv[..p], None => kv };
                pp.match_tags.iter().any(|t| t == k)
            })
        } else {
            matched_db
        };
        let present = match m.op { MarkOp::Present => detected, MarkOp::Absent => !detected };
        if present {
            if info == "." || info.is_empty() {
                info = m.tag.clone();
            } else if !info.split(';').any(|kv| kv == m.tag || kv.starts_with(&format!("{}=", m.tag))) {
                info.push(';');
                info.push_str(&m.tag);
            }
            modified = true;
        }
    }

    if let Some(keep) = &pp.samples_keep {
        if !samples.is_empty() {
            let filtered: Vec<String> = keep.iter().filter_map(|&i| samples.get(i).cloned()).collect();
            samples = filtered;
            modified = true;
        }
    }

    if let Some(ren) = &pp.rename_annots {
        if !ren.info.is_empty() {
            let new = rename_info_keys(&info, &ren.info);
            if new != info { info = new; modified = true; }
        }
        if !ren.format.is_empty() && !format.is_empty() {
            let new = rename_format_keys(&format, &ren.format);
            if new != format { format = new; modified = true; }
        }
        if !ren.filter.is_empty() && !filter.is_empty() && filter != "." {
            let new = rename_filter(&filter, &ren.filter);
            if new != filter { filter = new; modified = true; }
        }
    }

    if !modified { return LineAction::Keep; }

    let mut out = String::with_capacity(line.len() + 32);
    out.push_str(chrom); out.push('\t');
    out.push_str(cols_borrow[1]); out.push('\t');
    out.push_str(&id); out.push('\t');
    out.push_str(cols_borrow[3]); out.push('\t');
    out.push_str(cols_borrow[4]); out.push('\t');
    out.push_str(&qual); out.push('\t');
    out.push_str(&filter); out.push('\t');
    out.push_str(&info);
    if !format.is_empty() {
        out.push('\t');
        out.push_str(&format);
        for s in &samples {
            out.push('\t');
            out.push_str(s);
        }
    }
    LineAction::Replace(out)
}

fn filter_info(info: &str, tags: &[String], inverse: bool, drop_all: bool) -> String {
    if info == "." || info.is_empty() { return info.to_string(); }
    let mut out = String::with_capacity(info.len());
    let mut first = true;
    for kv in info.split(';') {
        let key = match memchr(b'=', kv.as_bytes()) { Some(p) => &kv[..p], None => kv };
        let listed = tags.iter().any(|t| t == key);
        let drop = if inverse { !(listed || (drop_all && false)) } else { listed || drop_all };
        if drop { continue; }
        if !first { out.push(';'); }
        out.push_str(kv);
        first = false;
    }
    if out.is_empty() { ".".into() } else { out }
}

fn filter_format_keys(format: &str, tags: &[String], inverse: bool, drop_all: bool) -> (String, Vec<usize>) {
    let mut new_keys = Vec::new();
    let mut removed = Vec::new();
    for (idx, k) in format.split(':').enumerate() {
        let listed = tags.iter().any(|t| t == k);
        let drop = if inverse { !listed && !drop_all } else { listed || drop_all };
        if drop && k != "GT" { removed.push(idx); } else { new_keys.push(k); }
    }
    (new_keys.join(":"), removed)
}

fn strip_format_indices(sample: &str, drop: &[usize]) -> String {
    if drop.is_empty() { return sample.to_string(); }
    sample.split(':').enumerate()
        .filter(|(i, _)| !drop.contains(i))
        .map(|(_, v)| v).collect::<Vec<_>>().join(":")
}

fn rename_info_keys(info: &str, map: &FxHashMap<String, String>) -> String {
    if info == "." || info.is_empty() { return info.to_string(); }
    let mut out = String::with_capacity(info.len());
    let mut first = true;
    for kv in info.split(';') {
        if !first { out.push(';'); }
        first = false;
        if let Some(eq) = memchr(b'=', kv.as_bytes()) {
            let key = &kv[..eq];
            if let Some(new) = map.get(key) {
                out.push_str(new);
                out.push_str(&kv[eq..]);
                continue;
            }
        } else if let Some(new) = map.get(kv) {
            out.push_str(new);
            continue;
        }
        out.push_str(kv);
    }
    out
}

fn rename_format_keys(format: &str, map: &FxHashMap<String, String>) -> String {
    format.split(':').map(|k| map.get(k).map(String::as_str).unwrap_or(k)).collect::<Vec<_>>().join(":")
}

fn rename_filter(filter: &str, map: &FxHashMap<String, String>) -> String {
    filter.split(';').map(|k| map.get(k).map(String::as_str).unwrap_or(k)).collect::<Vec<_>>().join(";")
}

fn render_id_template(t: &IdTemplate, chrom: &str, pos: &str, ref_: &str, alt: &str, id: &str, info: &str) -> String {
    let mut s = String::new();
    for tok in &t.tokens {
        match tok {
            IdTok::Lit(v) => s.push_str(v),
            IdTok::Chrom => s.push_str(chrom),
            IdTok::Pos => s.push_str(pos),
            IdTok::Ref => s.push_str(ref_),
            IdTok::Alt => s.push_str(alt),
            IdTok::Id => s.push_str(id),
            IdTok::InfoTag(tag) => {
                let v = info.split(';').find_map(|kv| {
                    let (k, val) = kv.split_once('=')?;
                    if k == tag { Some(val) } else { None }
                }).unwrap_or(".");
                s.push_str(v);
            }
        }
    }
    s
}

fn eval_predicate(p: &Predicate, line: &str) -> bool {
    let Some(rec) = parse_record_for_filter(line) else { return true; };
    p.engine.eval(&rec).map(|r| r.pass_site).unwrap_or(true)
}

fn parse_record_for_filter(line: &str) -> Option<crate::vcf::VcfRecord> {
    let cols: Vec<&str> = line.trim_end().split('\t').collect();
    if cols.len() < 8 { return None; }
    let pos: u32 = cols[1].parse().ok()?;
    let format = if cols.len() > 8 { Some(cols[8].to_string()) } else { None };
    let samples = if cols.len() > 9 { cols[9..].iter().map(|s| s.to_string()).collect() } else { Vec::new() };
    Some(crate::vcf::VcfRecord {
        chrom: cols[0].to_string(),
        pos,
        id: cols[2].to_string(),
        ref_allele: cols[3].to_string(),
        alt: cols[4].to_string(),
        qual: cols[5].to_string(),
        filter: cols[6].to_string(),
        info: cols[7].to_string(),
        format,
        samples,
        chr_id: 0,
        position: pos,
        offset: 0,
    })
}

pub fn parse_output_type(s: &str) -> Result<OutputKind> {
    let bytes = s.as_bytes();
    if bytes.is_empty() { bail!("-O: empty"); }
    let (kind_byte, lvl_str) = bytes.split_first().unwrap();
    let kind = match *kind_byte {
        b'v' => OutputKind::Vcf,
        b'z' => OutputKind::VcfGz(parse_level(lvl_str).unwrap_or(6)),
        b'u' => OutputKind::Bcf(0),
        b'b' => OutputKind::Bcf(parse_level(lvl_str).unwrap_or(6)),
        _ => bail!("-O: unknown type '{}', expected v|z|u|b", s),
    };
    Ok(kind)
}

fn parse_level(b: &[u8]) -> Option<u32> {
    if b.is_empty() { return None; }
    std::str::from_utf8(b).ok().and_then(|s| s.parse().ok()).filter(|&n: &u32| n <= 9)
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum OutputKind {
    Vcf,
    VcfGz(u32),
    Bcf(u32),
}

pub fn read_samples_file<P: AsRef<Path>>(p: P) -> Result<(Vec<String>, bool)> {
    let mut names: Vec<String> = Vec::new();
    let mut inverse = false;
    for (i, line) in BufReader::new(File::open(p.as_ref()).with_context(|| format!("open {:?}", p.as_ref()))?).lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') { continue; }
        if i == 0 {
            if let Some(rest) = t.strip_prefix('^') { inverse = true; names.push(rest.to_string()); continue; }
        }
        names.push(t.to_string());
    }
    Ok((names, inverse))
}

pub fn parse_samples_cli(spec: &str) -> (Vec<String>, bool) {
    let (inverse, body) = match spec.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, spec),
    };
    (body.split(',').filter(|s| !s.is_empty()).map(str::to_string).collect(), inverse)
}

pub fn resolve_samples_keep(input_samples: &[String], requested: &[String], inverse: bool) -> Vec<usize> {
    let req_set: std::collections::HashSet<&str> = requested.iter().map(|s| s.as_str()).collect();
    input_samples.iter().enumerate()
        .filter_map(|(i, s)| {
            let listed = req_set.contains(s.as_str());
            if inverse { (!listed).then_some(i) } else { listed.then_some(i) }
        })
        .collect()
}

pub fn read_columns_file<P: AsRef<Path>>(p: P) -> Result<(Vec<String>, Vec<Option<String>>)> {
    let mut cols = Vec::new();
    let mut types = Vec::new();
    for line in BufReader::new(File::open(p.as_ref())?).lines() {
        let line = line?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') { continue; }
        let mut it = t.splitn(2, char::is_whitespace);
        let name = it.next().unwrap_or("").to_string();
        let typ = it.next().map(|s| s.trim().to_string());
        if name.is_empty() { continue; }
        cols.push(name);
        types.push(typ);
    }
    Ok((cols, types))
}

#[derive(Clone, Default)]
pub struct RegionFilter {
    pub by_chr: FxHashMap<String, Vec<(u32, u32)>>,
}

impl RegionFilter {
    pub fn has_index_for<P: AsRef<Path>>(input: P) -> Option<std::path::PathBuf> {
        let p = input.as_ref();
        let s = p.to_string_lossy();
        for ext in &[".tbi", ".csi"] {
            let cand = std::path::PathBuf::from(format!("{}{}", s, ext));
            if cand.exists() { return Some(cand); }
        }
        None
    }

    /// Stream lines from a BGZF file using CSI index to skip irrelevant blocks.
    /// Returns lines that pass the region filter. If index not present,
    /// falls back to streaming everything through `line_passes`.
    pub fn stream_with_index<P: AsRef<Path>>(
        &self,
        input: P,
        contigs: &[String],
        mut cb: impl FnMut(&str),
    ) -> anyhow::Result<()> {
        use crate::bgzf::BgzfSeekReader;
        use crate::csi::CsiQuery;
        let path = input.as_ref();
        let Some(idx_path) = Self::has_index_for(path) else {
            anyhow::bail!("no .tbi/.csi index for {:?}", path);
        };
        let q = CsiQuery::open(&idx_path)?;
        let mut chunks: Vec<(u64, u64)> = Vec::new();
        for (rid, contig) in contigs.iter().enumerate() {
            if let Some(ranges) = self.by_chr.get(contig) {
                for _ in ranges {
                    let cs = q.query(rid, 0, u32::MAX);
                    chunks.extend(cs);
                }
            }
        }
        chunks.sort_unstable_by_key(|c| c.0);
        chunks.dedup();

        let mut reader = BgzfSeekReader::open(path)?;
        for (start, end) in chunks {
            reader.seek_to(start)?;
            while let Some(line) = reader.read_line()? {
                if line.is_empty() || line.starts_with('#') { continue; }
                if !self.line_passes(&line) { continue; }
                cb(&line);
                if reader_past(&reader, end) { break; }
            }
        }
        Ok(())
    }

    pub fn from_cli(spec: &str) -> Result<Self> {
        let mut f = RegionFilter::default();
        for raw in spec.split(',') {
            let item = raw.trim();
            if item.is_empty() { continue; }
            let (chr, beg, end) = parse_region_item(item)?;
            f.by_chr.entry(chr).or_default().push((beg, end));
        }
        f.finalize();
        Ok(f)
    }

    pub fn from_file<P: AsRef<Path>>(p: P) -> Result<Self> {
        let mut f = RegionFilter::default();
        let path = p.as_ref();
        let is_bed = matches!(path.extension().and_then(|e| e.to_str()), Some("bed") | Some("BED"));
        for line in BufReader::new(File::open(path)?).lines() {
            let line = line?;
            let t = line.trim();
            if t.is_empty() || t.starts_with('#') || t.starts_with("track") || t.starts_with("browser") { continue; }
            let parts: Vec<&str> = t.split('\t').collect();
            if parts.is_empty() { continue; }
            let chr = parts[0].to_string();
            let (beg, end) = match parts.len() {
                1 => (1, u32::MAX),
                2 => {
                    let b: u32 = parts[1].parse().context("region beg")?;
                    (if is_bed { b + 1 } else { b }, u32::MAX)
                }
                _ => {
                    let b: u32 = parts[1].parse().context("region beg")?;
                    let e: u32 = parts[2].parse().context("region end")?;
                    (if is_bed { b + 1 } else { b }, e)
                }
            };
            f.by_chr.entry(chr).or_default().push((beg, end));
        }
        f.finalize();
        Ok(f)
    }

    fn finalize(&mut self) {
        for v in self.by_chr.values_mut() {
            v.sort_unstable_by_key(|r| r.0);
            let mut merged: Vec<(u32, u32)> = Vec::with_capacity(v.len());
            for r in v.drain(..) {
                if let Some(last) = merged.last_mut() {
                    if r.0 <= last.1.saturating_add(1) { last.1 = last.1.max(r.1); continue; }
                }
                merged.push(r);
            }
            *v = merged;
        }
    }

    pub fn contains(&self, chr: &str, pos: u32) -> bool {
        let Some(ranges) = self.by_chr.get(chr) else { return false; };
        let idx = ranges.partition_point(|r| r.1 < pos);
        idx < ranges.len() && ranges[idx].0 <= pos
    }

    pub fn line_passes(&self, line: &str) -> bool {
        self.line_passes_mode(line, 1)
    }

    /// Overlap mode: 0 = POS in region; 1 = record overlaps (REF span);
    /// 2 = variant overlaps (REF + ALT span). Matches bcftools `--regions-overlap`.
    pub fn line_passes_mode(&self, line: &str, mode: u8) -> bool {
        let bytes = line.as_bytes();
        let Some(t1) = memchr(b'\t', bytes) else { return false; };
        let chr = &line[..t1];
        let rest = &line[t1 + 1..];
        let Some(t2) = memchr(b'\t', rest.as_bytes()) else { return false; };
        let pos_str = &rest[..t2];
        let Ok(pos) = pos_str.parse::<u32>() else { return false; };
        if mode == 0 { return self.contains(chr, pos); }
        let after = &rest[t2 + 1..];
        let cols: Vec<&str> = after.splitn(4, '\t').collect();
        if cols.len() < 3 { return self.contains(chr, pos); }
        let refa = cols[1];
        let alt = cols[2];
        let end = if mode == 2 {
            let max_alt: usize = alt.split(',').filter_map(|a| if a.starts_with('<') { None } else { Some(a.len()) }).max().unwrap_or(refa.len());
            pos + (refa.len().max(max_alt) as u32) - 1
        } else {
            pos + (refa.len() as u32) - 1
        };
        self.overlaps_range(chr, pos, end)
    }

    pub fn overlaps_range(&self, chr: &str, beg: u32, end: u32) -> bool {
        let Some(ranges) = self.by_chr.get(chr) else { return false; };
        let idx = ranges.partition_point(|r| r.1 < beg);
        idx < ranges.len() && ranges[idx].0 <= end
    }
}

fn reader_past(_r: &crate::bgzf::BgzfSeekReader, _end: u64) -> bool { false }

fn parse_region_item(s: &str) -> Result<(String, u32, u32)> {
    let (chr, range) = match s.split_once(':') {
        Some((c, r)) => (c.to_string(), Some(r)),
        None => (s.to_string(), None),
    };
    match range {
        None => Ok((chr, 1, u32::MAX)),
        Some(r) => {
            if let Some((b, e)) = r.split_once('-') {
                let beg: u32 = b.parse().context("region beg")?;
                let end: u32 = if e.is_empty() { u32::MAX } else { e.parse().context("region end")? };
                Ok((chr, beg, end))
            } else {
                let beg: u32 = r.parse().context("region beg")?;
                Ok((chr, beg, beg))
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum PairLogic { Snps, Indels, Both, All, Some_, Exact, Id }

impl PairLogic {
    pub fn parse(s: &str) -> Result<Self> {
        Ok(match s {
            "snps" => Self::Snps,
            "indels" => Self::Indels,
            "both" => Self::Both,
            "all" => Self::All,
            "some" => Self::Some_,
            "exact" => Self::Exact,
            "id" => Self::Id,
            _ => bail!("--pair-logic: unknown '{}', expected snps|indels|both|all|some|exact|id", s),
        })
    }
}

pub fn version_header_line() -> String {
    format!(
        "##kira_bt_annotateVersion={}+htslib-compat",
        env!("CARGO_PKG_VERSION")
    )
}

#[cfg(test)]
#[path = "../../tests/unit/annotate_postproc.rs"]
mod tests;
