//! Sample and read-group selection of `bcftools mpileup` (bam_sample.c):
//! `-s/-S` with `^` exclusion and `OLD NEW` renaming, `-G` read-group rows
//! (`ID [FILE] [SAMPLE]`, `^` to exclude), and the per-read sample lookup
//! through the `RG` tag. Samples are numbered in order of first appearance
//! across the input headers, as bcftools does.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, Result, bail};
use noodles_sam as sam;

#[derive(Clone, Debug, Default)]
pub struct SampleOpts {
    /// `-s/-S`: keep (`true`) or drop (`false`) the listed samples; the
    /// value is the output name.
    pub samples: Option<(bool, HashMap<String, String>)>,
    /// `-G`: keys are `ID` or `ID\tFILE`; the value is the new sample name,
    /// or `\t` to keep the header's.
    pub read_groups: Option<(bool, HashMap<String, String>)>,
    /// `--ignore-RG`: one sample per file, named after the file.
    pub ignore_rg: bool,
}

/// Whitespace-separated fields with backslash escapes (`hts_readlist` rows).
fn split_fields(row: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut escaped = false;
    let mut in_field = false;
    for c in row.chars() {
        if c == '\\' && !escaped {
            escaped = true;
            continue;
        }
        if c.is_whitespace() && !escaped {
            if in_field {
                out.push(std::mem::take(&mut cur));
                in_field = false;
            }
            continue;
        }
        cur.push(c);
        in_field = true;
        escaped = false;
    }
    if in_field {
        out.push(cur);
    }
    out
}

/// `hts_readlist`: a file (one item per line) or a comma list.
fn read_list(spec: &str, is_file: bool) -> Result<Vec<String>> {
    if is_file {
        let text = std::fs::read_to_string(spec).with_context(|| format!("read {spec}"))?;
        Ok(text.lines().map(|l| l.trim_end_matches('\r').to_string()).filter(|l| !l.is_empty()).collect())
    } else {
        Ok(spec.split(',').filter(|s| !s.is_empty()).map(|s| s.to_string()).collect())
    }
}

impl SampleOpts {
    pub fn parse(samples: Option<&str>, samples_file: Option<&Path>, read_groups: Option<&Path>, ignore_rg: bool) -> Result<Self> {
        let mut o = Self { ignore_rg, ..Default::default() };
        let mut add_samples = |spec: &str, is_file: bool| -> Result<()> {
            let (logic, spec) = match spec.strip_prefix('^') {
                Some(rest) => (false, rest),
                None => (true, spec),
            };
            let (_, map) = o.samples.get_or_insert_with(|| (logic, HashMap::new()));
            for row in read_list(spec, is_file)? {
                let f = split_fields(&row);
                if f.is_empty() {
                    continue;
                }
                let new = f.get(1).cloned().unwrap_or_else(|| f[0].clone());
                map.insert(f[0].clone(), new);
            }
            if let Some(s) = o.samples.as_mut() {
                s.0 = logic;
            }
            Ok(())
        };
        if let Some(s) = samples {
            add_samples(s, false)?;
        }
        if let Some(p) = samples_file {
            add_samples(&p.to_string_lossy(), true)?;
        }
        if let Some(p) = read_groups {
            let spec = p.to_string_lossy().into_owned();
            let (logic, spec) = match spec.strip_prefix('^') {
                Some(rest) => (false, rest.to_string()),
                None => (true, spec),
            };
            let mut map: HashMap<String, String> = HashMap::new();
            for row in read_list(&spec, true)? {
                let f = split_fields(&row);
                if f.is_empty() {
                    continue;
                }
                let (key, value) = match f.len() {
                    1 => (f[0].clone(), "\t".to_string()),
                    2 => (f[0].clone(), f[1].clone()),
                    _ => (format!("{}\t{}", f[0], f[1]), f[2].clone()),
                };
                match map.get(&key) {
                    None => {
                        map.insert(key, value);
                    }
                    Some(v) if *v != value => bail!(
                        "The read group \"{key}\" was assigned to two different samples: \"{v}\" and \"{value}\""
                    ),
                    _ => {}
                }
            }
            o.read_groups = Some((logic, map));
        }
        Ok(o)
    }
}

#[derive(Default, Clone, Debug)]
struct FileMap {
    default_idx: Option<usize>,
    /// Read group id -> sample; `None` when the group was seen but excluded.
    rg2idx: HashMap<String, Option<usize>>,
}

/// bcftools' `bam_smpl_t`.
#[derive(Default, Clone, Debug)]
pub struct SampleMap {
    pub names: Vec<String>,
    name2idx: HashMap<String, usize>,
    files: Vec<FileMap>,
}

impl SampleMap {
    fn sample_idx(&mut self, name: &str) -> usize {
        if let Some(&i) = self.name2idx.get(name) {
            return i;
        }
        let i = self.names.len();
        self.names.push(name.to_string());
        self.name2idx.insert(name.to_string(), i);
        i
    }

    fn add_readgroup(&mut self, file: usize, rg_id: &str, smpl_name: Option<&str>) {
        let ismpl = smpl_name.map(|n| self.sample_idx(n));
        let f = &mut self.files[file];
        if rg_id == "*" {
            f.default_idx = ismpl;
            return;
        }
        f.rg2idx.entry(rg_id.to_string()).or_insert(ismpl);
    }

    /// `bsmpl_keep_readgroup`: whether `-G` keeps this read group, renaming
    /// the sample when the row says so.
    fn keep_readgroup<'a>(rgs: &'a (bool, HashMap<String, String>), fname: &str, rg_id: &str, smpl_name: &mut Option<&'a str>) -> bool {
        let (logic, map) = rgs;
        let rg_smpl = map
            .get(rg_id)
            .or_else(|| map.get(&format!("{rg_id}\t{fname}")))
            .or_else(|| map.get(&format!("*\t{fname}")));
        if rg_smpl.is_none() && *logic {
            return false;
        }
        if rg_smpl.is_some() && !*logic {
            return false;
        }
        if let Some(s) = rg_smpl {
            if !s.starts_with('\t') {
                *smpl_name = Some(s.as_str());
            }
        }
        true
    }

    /// `bam_smpl_add_bam`: register one input; `None` when none of its read
    /// groups is selected (bcftools then drops the file).
    pub fn add_bam(&mut self, header: &sam::Header, fname: &str, opts: &SampleOpts) -> Option<usize> {
        let file = self.files.len();
        self.files.push(FileMap::default());
        if opts.ignore_rg {
            self.add_readgroup(file, "*", Some(fname));
            return Some(file);
        }
        let mut bam_smpls: HashSet<String> = HashSet::new();
        let mut first_smpl: Option<usize> = None;
        let mut nskipped = 0usize;
        for (id, rg) in header.read_groups() {
            let id = String::from_utf8_lossy(id.as_ref()).into_owned();
            let sm = rg
                .other_fields()
                .iter()
                .find(|(k, _)| k.as_ref() == b"SM")
                .map(|(_, v)| String::from_utf8_lossy(v.as_ref()).into_owned());
            // bcftools stops at the first @RG without an SM field.
            let Some(sm) = sm else { break };
            let mut r: &str = &sm;
            let mut accept = true;
            if let Some((logic, list)) = &opts.samples {
                match (list.get(r), *logic) {
                    (Some(_), false) => accept = false,
                    (None, false) => {}
                    (None, true) => accept = false,
                    (Some(new), true) => r = new,
                }
            }
            if accept {
                if let Some(rgs) = &opts.read_groups {
                    let mut renamed: Option<&str> = None;
                    accept = Self::keep_readgroup(rgs, fname, &id, &mut renamed);
                    if let Some(n) = renamed {
                        r = n;
                    }
                }
            }
            let r = r.to_string();
            if accept {
                self.add_readgroup(file, &id, Some(&r));
            } else {
                self.add_readgroup(file, &id, None);
                nskipped += 1;
            }
            if first_smpl.is_none() {
                first_smpl = self.name2idx.get(&r).copied();
            }
            bam_smpls.insert(r);
        }
        let nsmpls = bam_smpls.len();
        let mut smpl_name: Option<&str> = None;
        let mut accept_null_rg = true;
        if let Some(rgs) = &opts.read_groups {
            if !Self::keep_readgroup(rgs, fname, "?", &mut smpl_name) {
                accept_null_rg = false;
            }
        }
        if opts.samples.is_some() && first_smpl.is_none() {
            accept_null_rg = false;
        }
        if !accept_null_rg && first_smpl.is_none() {
            self.files.pop();
            return None;
        }
        if !accept_null_rg {
            return Some(file);
        }
        if nsmpls == 1 && nskipped == 0 {
            self.files[file].default_idx = first_smpl;
            return Some(file);
        }
        let name: String = match smpl_name {
            Some(n) => n.to_string(),
            None => match first_smpl {
                Some(i) => self.names[i].clone(),
                None => fname.to_string(),
            },
        };
        self.add_readgroup(file, "?", Some(&name));
        Some(file)
    }

    /// `bam_smpl_get_sample_id`: the sample of a read of file `bam_id` with
    /// read group `rg`; `None` drops the read.
    pub fn sample_id(&self, bam_id: usize, rg: Option<&[u8]>) -> Option<usize> {
        let f = self.files.get(bam_id)?;
        if let Some(d) = f.default_idx {
            return Some(d);
        }
        let key: std::borrow::Cow<'_, str> = match rg {
            Some(b) => String::from_utf8_lossy(b),
            None => "?".into(),
        };
        if let Some(v) = f.rg2idx.get(key.as_ref()) {
            return *v;
        }
        f.rg2idx.get("?").copied().flatten()
    }

    pub fn n_files(&self) -> usize {
        self.files.len()
    }

    /// Samples a file contributes to.
    pub fn file_samples(&self, bam_id: usize) -> Vec<usize> {
        let Some(f) = self.files.get(bam_id) else { return Vec::new() };
        let mut v: Vec<usize> = f.default_idx.into_iter().chain(f.rg2idx.values().copied().flatten()).collect();
        v.sort_unstable();
        v.dedup();
        v
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bam_samples.rs"]
mod tests;
