//! Allele-indexed field rewriting for Number=A/R/G tags and GT when the set of
//! alleles of a record changes (norm split/join, merge, view --trim-alt-alleles).

use crate::vcf::header::{FieldNumber, HeaderInfo};

/// Mapping from old allele index (0 = REF) to new allele index, or `None`
/// when the allele is dropped.
pub type AlleleMap = [Option<usize>];

#[inline]
pub fn diploid_gt_index(a: usize, b: usize) -> usize {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    hi * (hi + 1) / 2 + lo
}

#[inline]
pub fn n_diploid_genotypes(n_alleles: usize) -> usize {
    n_alleles * (n_alleles + 1) / 2
}

fn inverse_map(n_new: usize, map: &AlleleMap) -> Vec<Option<usize>> {
    let mut inv = vec![None; n_new];
    for (old, new) in map.iter().enumerate() {
        if let Some(j) = *new {
            if j < n_new && inv[j].is_none() {
                inv[j] = Some(old);
            }
        }
    }
    inv
}

/// Rewrite a comma-separated value according to `number`. Returns `None` when
/// the value does not have the declared cardinality (left unchanged by callers).
pub fn remap_value(
    value: &str,
    number: FieldNumber,
    n_old: usize,
    n_new: usize,
    map: &AlleleMap,
) -> Option<String> {
    if value == "." || value.is_empty() {
        return Some(value.to_string());
    }
    let vals: Vec<&str> = value.split(',').collect();
    let inv = inverse_map(n_new, map);
    match number {
        FieldNumber::A => {
            if vals.len() != n_old.saturating_sub(1) {
                return None;
            }
            let out: Vec<&str> = (1..n_new)
                .map(|j| inv[j].filter(|&i| i >= 1).map(|i| vals[i - 1]).unwrap_or("."))
                .collect();
            Some(out.join(","))
        }
        FieldNumber::R => {
            if vals.len() != n_old {
                return None;
            }
            let out: Vec<&str> = (0..n_new).map(|j| inv[j].map(|i| vals[i]).unwrap_or(".")).collect();
            Some(out.join(","))
        }
        FieldNumber::G => {
            if vals.len() == n_old {
                let out: Vec<&str> = (0..n_new).map(|j| inv[j].map(|i| vals[i]).unwrap_or(".")).collect();
                return Some(out.join(","));
            }
            if vals.len() != n_diploid_genotypes(n_old) {
                return None;
            }
            let mut out: Vec<&str> = Vec::with_capacity(n_diploid_genotypes(n_new));
            for b in 0..n_new {
                for a in 0..=b {
                    let v = match (inv[a], inv[b]) {
                        (Some(oa), Some(ob)) => vals[diploid_gt_index(oa, ob)],
                        _ => ".",
                    };
                    out.push(v);
                }
            }
            Some(out.join(","))
        }
        _ => None,
    }
}

/// Remap GT allele indices, preserving phasing separators. Dropped alleles
/// become missing.
pub fn remap_gt(gt: &str, map: &AlleleMap) -> String {
    let mut out = String::with_capacity(gt.len());
    let mut cur = String::new();
    let flush = |cur: &mut String, out: &mut String| {
        if cur.is_empty() {
            return;
        }
        if let Ok(n) = cur.parse::<usize>() {
            match map.get(n).copied().flatten() {
                Some(j) => out.push_str(&j.to_string()),
                None => out.push('.'),
            }
        } else {
            out.push_str(cur);
        }
        cur.clear();
    };
    for c in gt.chars() {
        if c == '/' || c == '|' {
            flush(&mut cur, &mut out);
            out.push(c);
        } else {
            cur.push(c);
        }
    }
    flush(&mut cur, &mut out);
    out
}

/// Allele indices of a GT, `None` for missing.
pub fn gt_alleles(gt: &str) -> Vec<Option<usize>> {
    gt.split(['/', '|'])
        .map(|a| if a == "." || a.is_empty() { None } else { a.parse::<usize>().ok() })
        .collect()
}

pub fn gt_is_phased(gt: &str) -> bool {
    gt.contains('|')
}

/// Rewrite the INFO column for a new allele set.
pub fn remap_info(info: &str, hdr: &HeaderInfo, n_old: usize, n_new: usize, map: &AlleleMap) -> String {
    if info == "." || info.is_empty() {
        return info.to_string();
    }
    let mut out = String::with_capacity(info.len());
    let mut first = true;
    for kv in info.split(';') {
        if !first {
            out.push(';');
        }
        first = false;
        match kv.split_once('=') {
            Some((k, v)) => {
                let num = hdr.info_number(k);
                if num.is_per_allele() {
                    if let Some(nv) = remap_value(v, num, n_old, n_new, map) {
                        if nv.is_empty() {
                            // No alleles left for a per-ALT tag: drop it.
                            if out.ends_with(';') {
                                out.pop();
                            }
                            first = out.is_empty();
                            continue;
                        }
                        out.push_str(k);
                        out.push('=');
                        out.push_str(&nv);
                        continue;
                    }
                }
                out.push_str(kv);
            }
            None => out.push_str(kv),
        }
    }
    if out.is_empty() { ".".to_string() } else { out }
}

/// Rewrite FORMAT sample columns for a new allele set. GT is remapped, A/R/G
/// tags are subset/expanded, everything else is copied.
pub fn remap_samples(
    format: &str,
    samples: &[&str],
    hdr: &HeaderInfo,
    n_old: usize,
    n_new: usize,
    map: &AlleleMap,
) -> Vec<String> {
    let keys: Vec<&str> = format.split(':').collect();
    let numbers: Vec<FieldNumber> = keys.iter().map(|k| hdr.format_number(k)).collect();
    let mut out = Vec::with_capacity(samples.len());
    for s in samples {
        let mut parts: Vec<String> = Vec::with_capacity(keys.len());
        for (i, v) in s.split(':').enumerate() {
            let Some(k) = keys.get(i) else {
                parts.push(v.to_string());
                continue;
            };
            if *k == "GT" {
                parts.push(remap_gt(v, map));
            } else if numbers[i].is_per_allele() {
                let nv = remap_value(v, numbers[i], n_old, n_new, map).unwrap_or_else(|| v.to_string());
                parts.push(if nv.is_empty() { ".".to_string() } else { nv });
            } else {
                parts.push(v.to_string());
            }
        }
        out.push(parts.join(":"));
    }
    out
}

/// Split a `k=v;k2;k3=v3` INFO string.
pub fn split_info(info: &str) -> Vec<(&str, Option<&str>)> {
    if info == "." || info.is_empty() {
        return Vec::new();
    }
    info.split(';')
        .map(|kv| match kv.split_once('=') {
            Some((k, v)) => (k, Some(v)),
            None => (kv, None),
        })
        .collect()
}

pub fn join_info(items: &[(String, Option<String>)]) -> String {
    if items.is_empty() {
        return ".".to_string();
    }
    let mut out = String::new();
    for (i, (k, v)) in items.iter().enumerate() {
        if i > 0 {
            out.push(';');
        }
        out.push_str(k);
        if let Some(v) = v {
            out.push('=');
            out.push_str(v);
        }
    }
    out
}

#[cfg(test)]
#[path = "../../tests/unit/vcf_alleles.rs"]
mod tests;
