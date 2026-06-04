use crate::cli::args::IsecArgs;
use crate::vcf::UnifiedVcfReader;
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Clone, Copy, Debug)]
enum Collapse { None, All, Snps, Indels, Both, Some_, Id }

fn parse_collapse(s: Option<&str>) -> Result<Collapse> {
    let Some(s) = s else { return Ok(Collapse::None); };
    match s {
        "none" => Ok(Collapse::None),
        "all" => Ok(Collapse::All),
        "snps" => Ok(Collapse::Snps),
        "indels" => Ok(Collapse::Indels),
        "both" => Ok(Collapse::Both),
        "some" => Ok(Collapse::Some_),
        "id" => Ok(Collapse::Id),
        _ => bail!("-c: unknown collapse mode '{}', expected none|all|snps|indels|both|some|id", s),
    }
}

fn parse_apply_filters(s: Option<&str>) -> Option<Vec<String>> {
    s.map(|s| s.split(',').map(|t| t.trim().to_string()).collect())
}

fn filter_passes(line: &str, allow: &Option<Vec<String>>) -> bool {
    let Some(allow) = allow else { return true; };
    let cols: Vec<&str> = line.splitn(8, '\t').collect();
    if cols.len() < 7 { return true; }
    let f = cols[6];
    f.split(';').any(|tok| allow.iter().any(|a| a == tok))
}

fn is_snp(refa: &str, alt: &str) -> bool {
    if refa.len() != 1 { return false; }
    alt.split(',').all(|a| a.len() == 1 && a != "*" && a != ".")
}

fn is_indel(refa: &str, alt: &str) -> bool {
    alt.split(',').any(|a| a.len() != refa.len() && a != "*" && a != ".")
}

fn alt_set(alt: &str) -> Vec<&str> {
    alt.split(',').filter(|a| *a != "*" && *a != ".").collect()
}

fn site_key_for(line: &str, mode: Collapse) -> Option<String> {
    let cols: Vec<&str> = line.splitn(8, '\t').collect();
    if cols.len() < 5 { return None; }
    let chrom = cols[0];
    let pos = cols[1];
    let id = cols[2];
    let refa = cols[3];
    let alt = cols[4];
    let k = match mode {
        Collapse::None => format!("{}\t{}\t{}\t{}", chrom, pos, refa, alt),
        Collapse::All => format!("{}\t{}", chrom, pos),
        Collapse::Snps => if is_snp(refa, alt) { format!("snp\t{}\t{}", chrom, pos) } else { format!("{}\t{}\t{}\t{}", chrom, pos, refa, alt) },
        Collapse::Indels => if is_indel(refa, alt) { format!("indel\t{}\t{}", chrom, pos) } else { format!("{}\t{}\t{}\t{}", chrom, pos, refa, alt) },
        Collapse::Both => {
            if is_snp(refa, alt) { format!("snp\t{}\t{}", chrom, pos) }
            else if is_indel(refa, alt) { format!("indel\t{}\t{}", chrom, pos) }
            else { format!("{}\t{}\t{}\t{}", chrom, pos, refa, alt) }
        }
        Collapse::Some_ => {
            let mut alts = alt_set(alt);
            alts.sort();
            format!("some\t{}\t{}\t{}\t{}", chrom, pos, refa, alts.join("|"))
        }
        Collapse::Id => format!("id\t{}", id),
    };
    Some(k)
}

fn key_matches_some(a: &str, b: &str) -> bool {
    let pa: Vec<&str> = a.splitn(5, '\t').collect();
    let pb: Vec<&str> = b.splitn(5, '\t').collect();
    if pa.len() < 5 || pb.len() < 5 { return a == b; }
    if pa[0] != "some" || pa[..4] != pb[..4] { return a == b; }
    let sa: std::collections::HashSet<&str> = pa[4].split('|').collect();
    let sb: std::collections::HashSet<&str> = pb[4].split('|').collect();
    sa.intersection(&sb).next().is_some()
}

/// BED/region set for `-R`/`-r` filtering. Intervals are stored 0-based half-open per chromosome,
/// sorted by start; lookup is by 1-based VCF position.
struct RegionSet {
    by_chr: FxHashMap<String, Vec<(u32, u32)>>,
}
impl RegionSet {
    fn contains(&self, chrom: &str, pos1: u32) -> bool {
        let Some(v) = self.by_chr.get(chrom) else { return false };
        let p0 = pos1.saturating_sub(1);
        let idx = v.partition_point(|&(s, _)| s <= p0);
        idx > 0 && {
            let (s, e) = v[idx - 1];
            s <= p0 && p0 < e
        }
    }
}
fn load_regions(file: Option<&std::path::Path>, regions: Option<&str>) -> Result<Option<RegionSet>> {
    let mut by_chr: FxHashMap<String, Vec<(u32, u32)>> = FxHashMap::default();
    if let Some(f) = file {
        for line in BufReader::new(File::open(f).with_context(|| format!("open regions {:?}", f))?).lines() {
            let l = line?;
            let t = l.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            let mut it = t.split('\t');
            let (Some(c), Some(s), Some(e)) = (it.next(), it.next(), it.next()) else { continue };
            if let (Ok(s), Ok(e)) = (s.parse::<u32>(), e.parse::<u32>()) {
                by_chr.entry(c.to_string()).or_default().push((s, e));
            }
        }
    } else if let Some(r) = regions {
        // -r chr[:start-end][,...]  (1-based inclusive)
        for tok in r.split(',') {
            let tok = tok.trim();
            if tok.is_empty() { continue; }
            match tok.split_once(':') {
                None => { by_chr.entry(tok.to_string()).or_default().push((0, u32::MAX)); }
                Some((c, range)) => {
                    let (a, b) = range.split_once('-').unwrap_or((range, range));
                    if let (Ok(a), Ok(b)) = (a.parse::<u32>(), b.parse::<u32>()) {
                        by_chr.entry(c.to_string()).or_default().push((a.saturating_sub(1), b));
                    }
                }
            }
        }
    } else {
        return Ok(None);
    }
    for v in by_chr.values_mut() { v.sort_unstable(); }
    Ok(Some(RegionSet { by_chr }))
}

pub fn cmd_isec(args: IsecArgs) -> Result<()> {
    let mut inputs: Vec<PathBuf> = args.inputs.clone();
    if let Some(fl) = &args.file_list {
        for line in BufReader::new(File::open(fl)?).lines() {
            let l = line?;
            let t = l.trim();
            if t.is_empty() || t.starts_with('#') { continue; }
            inputs.push(PathBuf::from(t));
        }
    }
    if inputs.len() < 2 { bail!("isec: need at least 2 inputs"); }

    let nfiles_spec = parse_nfiles(args.nfiles.as_deref(), inputs.len(), args.complement)?;
    let write_files = parse_write(args.write.as_deref(), inputs.len())?;
    let collapse = parse_collapse(args.collapse.as_deref())?;
    let apply_filters = parse_apply_filters(args.apply_filters.as_deref());
    let region_set = load_regions(args.regions_file.as_deref(), args.regions.as_deref())?;

    let mut sites: FxHashMap<String, Vec<bool>> = FxHashMap::default();
    let mut headers: Vec<Vec<String>> = Vec::with_capacity(inputs.len());
    let mut records: Vec<Vec<(String, String)>> = Vec::with_capacity(inputs.len());

    for (i, p) in inputs.iter().enumerate() {
        let mut r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        headers.push(r.header()?);
        let mut recs: Vec<(String, String)> = Vec::new();
        while let Some(line) = r.read_line()? {
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            if !filter_passes(&line, &apply_filters) { continue; }
            if let Some(rs) = &region_set {
                let mut c = line.splitn(3, '\t');
                let chrom = c.next().unwrap_or("");
                let pos: u32 = c.next().and_then(|p| p.parse().ok()).unwrap_or(0);
                if !rs.contains(chrom, pos) { continue; }
            }
            let Some(k) = site_key_for(&line, collapse) else { continue; };
            sites.entry(k.clone()).or_insert_with(|| vec![false; inputs.len()])[i] = true;
            recs.push((k, line));
        }
        records.push(recs);
    }

    if matches!(collapse, Collapse::Some_) {
        let keys: Vec<String> = sites.keys().cloned().collect();
        let mut merged: FxHashMap<String, Vec<bool>> = FxHashMap::default();
        let mut used = vec![false; keys.len()];
        for i in 0..keys.len() {
            if used[i] { continue; }
            used[i] = true;
            let mut acc = sites.get(&keys[i]).cloned().unwrap_or_default();
            for j in (i + 1)..keys.len() {
                if used[j] { continue; }
                if key_matches_some(&keys[i], &keys[j]) {
                    used[j] = true;
                    let other = sites.get(&keys[j]).cloned().unwrap_or_default();
                    for (a, b) in acc.iter_mut().zip(other.iter()) { *a = *a || *b; }
                }
            }
            merged.insert(keys[i].clone(), acc);
        }
        sites = merged;
    }

    let dst = args.prefix.clone();
    if let Some(d) = &dst { fs::create_dir_all(d).with_context(|| format!("create dir {:?}", d))?; }

    let mut writers: Vec<Option<BufWriter<File>>> = Vec::with_capacity(inputs.len());
    if let Some(dir) = &dst {
        for i in 0..inputs.len() {
            if write_files.is_some() && !write_files.as_ref().unwrap().contains(&i) {
                writers.push(None); continue;
            }
            let path = dir.join(format!("{:04}.vcf", i));
            let mut w = BufWriter::with_capacity(1 << 20, File::create(&path)?);
            for h in &headers[i] { writeln!(w, "{}", h)?; }
            writers.push(Some(w));
        }
    }

    for (i, recs) in records.iter().enumerate() {
        if let Some(d) = &dst {
            if write_files.as_ref().map_or(false, |w| !w.contains(&i)) { continue; }
            let path = d.join(format!("sites_{:04}.txt", i));
            let mut sw = BufWriter::with_capacity(1 << 20, File::create(&path)?);
            for (k, line) in recs {
                let presence = sites.get(k).cloned().unwrap_or_default();
                if !matches_nfiles(&presence, &nfiles_spec) { continue; }
                if let Some(w) = writers[i].as_mut() { writeln!(w, "{}", line)?; }
                let bits: String = presence.iter().map(|b| if *b { '1' } else { '0' }).collect();
                let cols: Vec<&str> = line.splitn(8, '\t').collect();
                if cols.len() >= 5 {
                    writeln!(sw, "{}\t{}\t{}\t{}\t{}", cols[0], cols[1], cols[3], cols[4], bits)?;
                }
            }
            sw.flush()?;
        } else {
            let mut stdout = BufWriter::with_capacity(1 << 20, std::io::stdout());
            if i == 0 {
                for h in &headers[0] { writeln!(stdout, "{}", h)?; }
            }
            for (k, line) in recs {
                let presence = sites.get(k).cloned().unwrap_or_default();
                if !matches_nfiles(&presence, &nfiles_spec) { continue; }
                if write_files.as_ref().map_or(true, |w| w.contains(&i)) {
                    writeln!(stdout, "{}", line)?;
                }
            }
            stdout.flush()?;
            break;
        }
    }

    for w in writers.iter_mut() {
        if let Some(w) = w { w.flush()?; }
    }
    Ok(())
}

#[derive(Debug)]
enum NSpec { Exact(usize), AtLeast(usize), Mask(Vec<bool>) }

fn parse_nfiles(s: Option<&str>, n: usize, complement: bool) -> Result<NSpec> {
    if complement { return Ok(NSpec::Mask(std::iter::once(true).chain(std::iter::repeat(false).take(n - 1)).collect())); }
    let Some(s) = s else { return Ok(NSpec::AtLeast(1)); };
    if let Some(rest) = s.strip_prefix('=') {
        return Ok(NSpec::Exact(rest.parse().context("-n =N: parse N")?));
    }
    if let Some(rest) = s.strip_prefix('+') {
        return Ok(NSpec::AtLeast(rest.parse().context("-n +N: parse N")?));
    }
    if let Some(rest) = s.strip_prefix('~') {
        let bits: Vec<bool> = rest.chars().map(|c| c == '1').collect();
        if bits.len() != n { bail!("-n ~MASK: mask length {} != #files {}", bits.len(), n); }
        return Ok(NSpec::Mask(bits));
    }
    Ok(NSpec::Exact(s.parse().context("-n N: parse N")?))
}

fn matches_nfiles(presence: &[bool], spec: &NSpec) -> bool {
    let n_present = presence.iter().filter(|b| **b).count();
    match spec {
        NSpec::Exact(k) => n_present == *k,
        NSpec::AtLeast(k) => n_present >= *k,
        NSpec::Mask(mask) => presence.iter().zip(mask.iter()).all(|(p, m)| !*m || *p),
    }
}

fn parse_write(s: Option<&str>, n: usize) -> Result<Option<Vec<usize>>> {
    let Some(s) = s else { return Ok(None); };
    let mut out = Vec::new();
    for tok in s.split(',') {
        let i: usize = tok.parse().context("-w: parse index")?;
        if i == 0 || i > n { bail!("-w: index {i} out of range 1..={n}"); }
        out.push(i - 1);
    }
    Ok(Some(out))
}
