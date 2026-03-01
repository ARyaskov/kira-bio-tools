use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io::{Cursor, Read};

use crate::VcfReader;
use crate::cli::args::IsecArgs;

pub fn cmd_isec(args: IsecArgs) -> Result<()> {
    let cfg = parse_isec_args(&args.bcftools_args)?;
    let mut inputs = Vec::<InputData>::new();

    for path in &args.inputs {
        let mut r = VcfReader::open(path)?;
        let headers = r.header()?;
        let mut recs = Vec::<Rec>::new();
        while let Some(rec) = r.next_record()? {
            if !cfg.matches_regions(&rec.chrom, rec.pos) {
                continue;
            }
            if !cfg.matches_expr(&rec.ref_allele) {
                continue;
            }
            let key = format!(
                "{}\t{}\t{}\t{}",
                rec.chrom, rec.pos, rec.ref_allele, rec.alt
            );
            let full_line = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}",
                rec.chrom,
                rec.pos,
                rec.id,
                rec.ref_allele,
                rec.alt,
                rec.qual,
                rec.filter,
                rec.info,
                match &rec.format {
                    Some(f) => {
                        let mut s = String::new();
                        s.push('\t');
                        s.push_str(f);
                        for sample in &rec.samples {
                            s.push('\t');
                            s.push_str(sample);
                        }
                        s
                    }
                    None => String::new(),
                }
            );
            recs.push(Rec {
                chrom: rec.chrom,
                pos: rec.pos,
                ref_allele: rec.ref_allele,
                alt: rec.alt,
                key,
                full_line,
            });
        }
        inputs.push(InputData {
            headers,
            records: recs,
        });
    }

    let n_inputs = inputs.len();
    let mut map = BTreeMap::<String, Vec<bool>>::new();
    let mut selected_lines = HashMap::<(usize, String), String>::new();

    for (i, inp) in inputs.iter().enumerate() {
        for r in &inp.records {
            map.entry(r.key.clone())
                .or_insert_with(|| vec![false; n_inputs])[i] = true;
            selected_lines.insert((i, r.key.clone()), r.full_line.clone());
        }
    }

    let mut selected_keys = Vec::<(String, Vec<bool>)>::new();
    for (k, bits) in map {
        if cfg.pass_bits(&bits) {
            selected_keys.push((k, bits));
        }
    }

    if let Some(w) = cfg.write_index {
        let idx = w.saturating_sub(1);
        if let Some(inp) = inputs.get(idx) {
            for h in &inp.headers {
                println!("{h}");
            }
            for (k, _bits) in &selected_keys {
                if let Some(line) = selected_lines.get(&(idx, k.clone())) {
                    println!("{line}");
                }
            }
        }
        return Ok(());
    }

    for (k, bits) in selected_keys {
        let mut cols = k.split('\t');
        let chrom = cols.next().unwrap_or(".");
        let pos = cols.next().unwrap_or("0");
        let r = cols.next().unwrap_or(".");
        let a = cols.next().unwrap_or(".");
        let bitmap: String = bits.iter().map(|b| if *b { '1' } else { '0' }).collect();
        println!("{chrom}\t{pos}\t{r}\t{a}\t{bitmap}");
    }

    Ok(())
}

#[derive(Default)]
struct IsecCfg {
    n_mode: Option<NMode>,
    complement: bool,
    write_index: Option<usize>,
    regions: Vec<Region>,
    include_strlen_ref_eq_2: bool,
}

#[derive(Clone, Copy)]
enum NMode {
    Eq(usize),
    Ge(usize),
    Le(usize),
}

#[derive(Clone)]
struct Region {
    chrom: String,
    start: Option<u32>,
    end: Option<u32>,
}

impl IsecCfg {
    fn matches_regions(&self, chrom: &str, pos: u32) -> bool {
        if self.regions.is_empty() {
            return true;
        }
        self.regions.iter().any(|r| {
            if r.chrom != chrom {
                return false;
            }
            if let Some(s) = r.start {
                if pos < s {
                    return false;
                }
            }
            if let Some(e) = r.end {
                if pos > e {
                    return false;
                }
            }
            true
        })
    }

    fn matches_expr(&self, ref_allele: &str) -> bool {
        if self.include_strlen_ref_eq_2 {
            return ref_allele.len() == 2;
        }
        true
    }

    fn pass_bits(&self, bits: &[bool]) -> bool {
        let c = bits.iter().filter(|b| **b).count();
        if self.complement {
            return bits.first().copied().unwrap_or(false) && c == 1;
        }
        if let Some(n) = self.n_mode {
            return match n {
                NMode::Eq(v) => c == v,
                NMode::Ge(v) => c >= v,
                NMode::Le(v) => c <= v,
            };
        }
        true
    }
}

fn parse_isec_args(args: &[String]) -> Result<IsecCfg> {
    let mut cfg = IsecCfg::default();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-n" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -n"))?;
                cfg.n_mode = Some(parse_n(v)?);
            }
            "-C" => cfg.complement = true,
            "-w" => {
                i += 1;
                cfg.write_index = Some(
                    args.get(i)
                        .ok_or_else(|| anyhow::anyhow!("missing value for -w"))?
                        .parse::<usize>()?,
                );
            }
            "-r" => {
                i += 1;
                let v = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -r"))?;
                for part in v.split(',') {
                    if let Some(r) = parse_region(part) {
                        cfg.regions.push(r);
                    }
                }
            }
            "-R" | "-T" => {
                i += 1;
                let p = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -R/-T"))?;
                for r in load_regions_from_file(p)? {
                    cfg.regions.push(r);
                }
            }
            "-i" => {
                i += 1;
                let e = args
                    .get(i)
                    .ok_or_else(|| anyhow::anyhow!("missing value for -i"))?;
                if e.trim() == "STRLEN(REF)==2" {
                    cfg.include_strlen_ref_eq_2 = true;
                }
            }
            _ => {}
        }
        i += 1;
    }
    Ok(cfg)
}

fn parse_n(s: &str) -> Result<NMode> {
    if let Some(v) = s.strip_prefix('=') {
        return Ok(NMode::Eq(v.parse::<usize>()?));
    }
    if let Some(v) = s.strip_prefix('+') {
        return Ok(NMode::Ge(v.parse::<usize>()?));
    }
    if let Some(v) = s.strip_prefix('-') {
        return Ok(NMode::Le(v.parse::<usize>()?));
    }
    Ok(NMode::Eq(s.parse::<usize>()?))
}

fn parse_region(s: &str) -> Option<Region> {
    if let Some((chrom, rest)) = s.split_once(':') {
        if let Some((a, b)) = rest.split_once('-') {
            return Some(Region {
                chrom: chrom.to_string(),
                start: a.parse::<u32>().ok(),
                end: b.parse::<u32>().ok(),
            });
        }
        return Some(Region {
            chrom: chrom.to_string(),
            start: rest.parse::<u32>().ok(),
            end: rest.parse::<u32>().ok(),
        });
    }
    Some(Region {
        chrom: s.to_string(),
        start: None,
        end: None,
    })
}

fn load_regions_from_file(path: &str) -> Result<Vec<Region>> {
    let bytes = fs::read(path)?;
    let text = decode_bytes(&bytes, path)?;
    let mut out = Vec::new();
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        if t.contains(':') {
            if let Some(r) = parse_region(t) {
                out.push(r);
            }
            continue;
        }
        let cols: Vec<&str> = t.split_whitespace().collect();
        if cols.len() >= 3 {
            let start0 = cols[1].parse::<u32>().ok();
            let end = cols[2].parse::<u32>().ok();
            out.push(Region {
                chrom: cols[0].to_string(),
                start: start0.map(|v| v + 1),
                end,
            });
        } else if cols.len() == 2 {
            let p = cols[1].parse::<u32>().ok();
            out.push(Region {
                chrom: cols[0].to_string(),
                start: p,
                end: p,
            });
        }
    }
    Ok(out)
}

fn decode_bytes(bytes: &[u8], path_hint: &str) -> Result<String> {
    let is_gz = (bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b)
        || path_hint.ends_with(".gz")
        || path_hint.ends_with(".bgz");
    if is_gz {
        let mut s = String::new();
        MultiGzDecoder::new(Cursor::new(bytes)).read_to_string(&mut s)?;
        return Ok(s);
    }
    String::from_utf8(bytes.to_vec()).map_err(|_| anyhow::anyhow!("non-UTF8 region file"))
}

struct InputData {
    headers: Vec<String>,
    records: Vec<Rec>,
}

struct Rec {
    chrom: String,
    pos: u32,
    ref_allele: String,
    alt: String,
    key: String,
    full_line: String,
}
