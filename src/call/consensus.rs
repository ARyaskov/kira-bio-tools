use anyhow::{Context, Result, ensure};
use indexmap::IndexMap;

use crate::bam::bcf_layout::kt_fisher_exact;
use crate::bam::errmod::lgamma;
use crate::util::{erfc, fmt_g};
use crate::vcf::VcfRecord;
use crate::vcf::header::{FieldNumber, HeaderInfo};

use super::stream::CallConfig;

struct Spectrum {
    reference: f64,
    variant: f64,
    folded_reference: f64,
    folded_variant: f64,
    frequency: f64,
    count: usize,
}

fn spectrum(likelihoods: &[[f64; 3]], indel: bool) -> Spectrum {
    let copies = 2 * likelihoods.len();
    let mut previous = vec![0.0; copies + 1];
    let mut next = vec![0.0; copies + 1];
    previous[0] = 1.0;
    let (mut low, mut high) = (0usize, 0usize);
    for (sample, row) in likelihoods.iter().enumerate() {
        let total = 2 * sample;
        while low < high && previous[low] < 1e-20 {
            previous[low] = 0.0;
            next[low] = 0.0;
            low += 1;
        }
        while high > low && previous[high] < 1e-20 {
            previous[high] = 0.0;
            next[high] = 0.0;
            high -= 1;
        }
        high += 2;
        let mut sum = 0.0;
        for k in low..=high {
            let remaining = total as i64 - k as i64;
            let a = (remaining + 1) as f64 * (remaining + 2) as f64 * row[0] * previous[k];
            let b = if k > 0 {
                k as f64 * (remaining + 2) as f64 * (2.0 * row[1]) * previous[k - 1]
            } else {
                0.0
            };
            let c = if k > 1 {
                k as f64 * (k - 1) as f64 * row[2] * previous[k - 2]
            } else {
                0.0
            };
            next[k] = a + b + c;
            sum += next[k];
        }
        for value in &mut next[low..=high] {
            *value /= sum;
        }
        if low > 0 {
            next[low - 1] = 0.0;
        }
        if low > 1 {
            next[low - 2] = 0.0;
        }
        if sample + 1 < likelihoods.len() {
            next[high + 1] = 0.0;
            next[high + 2] = 0.0;
        }
        std::mem::swap(&mut previous, &mut next);
    }
    let mut prior: Vec<f64> = (0..copies).map(|k| 1e-3 / (copies - k) as f64).collect();
    prior.push(1.0 - prior.iter().sum::<f64>());
    if indel {
        for value in &mut prior[..copies] {
            *value *= 0.15;
        }
        prior[copies] = 1.0 - prior[copies] * 0.15;
    }
    let normalization: f64 = prior
        .iter()
        .zip(&previous)
        .map(|(prior, likelihood)| prior * likelihood)
        .sum();
    let posterior: Vec<_> = prior
        .iter()
        .zip(&previous)
        .map(|(prior, likelihood)| prior * likelihood / normalization)
        .collect();
    let folded_normalization: f64 = (0..=copies)
        .map(|k| (prior[k] + prior[copies - k]) / 2.0 * previous[k])
        .sum();
    let folded_variant = (1..copies)
        .map(|k| (prior[k] + prior[copies - k]) / 2.0 * previous[k])
        .sum::<f64>()
        / folded_normalization;
    let mut best = 0;
    for k in 1..=copies {
        if previous[k] > previous[best] {
            best = k;
        }
    }
    Spectrum {
        reference: posterior[copies],
        variant: posterior[..copies].iter().sum(),
        folded_reference: (prior[copies] + prior[0]) / 2.0 * (previous[copies] + previous[0])
            / folded_normalization,
        folded_variant,
        frequency: posterior
            .iter()
            .enumerate()
            .map(|(k, probability)| k as f64 * probability)
            .sum::<f64>()
            / copies as f64,
        count: copies - best,
    }
}

fn minimize_frequency(mut a: f64, mut b: f64, likelihoods: &[[f64; 3]]) -> f64 {
    let objective = |frequency: f64| {
        if !(0.0..=1.0).contains(&frequency) {
            return 1e300;
        }
        let weights = [
            (1.0 - frequency).powi(2),
            2.0 * frequency * (1.0 - frequency),
            frequency.powi(2),
        ];
        let (mut product, mut logarithm) = (1.0f64, 0.0);
        for row in likelihoods {
            product *= row[0] * weights[0] + row[1] * weights[1] + row[2] * weights[2];
            if product < 1e-200 {
                logarithm -= product.ln();
                product = 1.0;
            }
        }
        logarithm - product.ln()
    };
    let gold1 = (std::f64::consts::GOLDEN_RATIO * 1e10).round() / 1e10;
    let (gold2, tiny) = (0.381_966_011_3, 1e-20);
    let (mut fa, mut fb) = (objective(a), objective(b));
    if fb > fa {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = b + gold1 * (b - a);
    let mut fc = objective(c);
    while fb > fc {
        let bound = b + 100.0 * (c - b);
        let r = (b - a) * (fb - fc);
        let q = (b - c) * (fb - fa);
        let denominator = if (q - r).abs() < tiny {
            if q > r { tiny } else { -tiny }
        } else {
            q - r
        };
        let mut u = b - ((b - c) * q - (b - a) * r) / (2.0 * denominator);
        let mut fu;
        if (b > u && u > c) || (b < u && u < c) {
            fu = objective(u);
            if fu < fc {
                a = b;
                b = u;
                fb = fu;
                break;
            } else if fu > fb {
                c = u;
                break;
            }
            u = c + gold1 * (c - b);
            fu = objective(u);
        } else if (c > u && u > bound) || (c < u && u < bound) {
            fu = objective(u);
            if fu < fc {
                b = c;
                c = u;
                u = c + gold1 * (c - b);
                fb = fc;
                fc = fu;
                fu = objective(u);
            } else {
                a = b;
                b = c;
                c = u;
                fb = fc;
                break;
            }
        } else {
            u = if (u > bound && bound > c) || (u < bound && bound < c) {
                bound
            } else {
                c + gold1 * (c - b)
            };
            fu = objective(u);
        }
        a = b;
        b = c;
        c = u;
        fa = fb;
        fb = fc;
        fc = fu;
    }
    if a > c {
        std::mem::swap(&mut a, &mut c);
    }
    let (mut e, mut d) = (0.0f64, 0.0f64);
    let (mut w, mut v, mut fv, mut fw) = (b, b, fb, fb);
    for _ in 0..100 {
        let middle = 0.5 * (a + c);
        let tolerance = 1e-5 * b.abs() + tiny;
        let twice = 2.0 * tolerance;
        if (b - middle).abs() <= twice - 0.5 * (c - a) {
            break;
        }
        if e.abs() > tolerance {
            let r = (b - w) * (fb - fv);
            let q = (b - v) * (fb - fw);
            let mut p = (b - v) * q - (b - w) * r;
            let q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            let q = q.abs();
            let old = e;
            e = d;
            if p.abs() >= (0.5 * q * old).abs() || p <= q * (a - b) || p >= q * (c - b) {
                e = if b >= middle { a - b } else { c - b };
                d = gold2 * e;
            } else {
                d = p / q;
                let u = b + d;
                if u - a < twice || c - u < twice {
                    d = if middle > b { tolerance } else { -tolerance };
                }
            }
        } else {
            e = if b >= middle { a - b } else { c - b };
            d = gold2 * e;
        }
        let u = b + if d.abs() >= tolerance {
            d
        } else if d > 0.0 {
            tolerance
        } else {
            -tolerance
        };
        let fu = objective(u);
        if fu <= fb {
            if u >= b {
                a = b;
            } else {
                c = b;
            }
            v = w;
            w = b;
            b = u;
            fv = fw;
            fw = fb;
            fb = fu;
        } else {
            if u < b {
                a = u;
            } else {
                c = u;
            }
            if fu <= fw || w == b {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if fu <= fv || v == b || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    b
}

fn estimate_frequency(likelihoods: &[[f64; 3]]) -> Option<f64> {
    let mut counts = [0usize; 3];
    for row in likelihoods {
        if row == &[1.0; 3] {
            continue;
        }
        let best = if row[0] > row[1] { 0 } else { 1 };
        let best = if row[best] > row[2] { best } else { 2 };
        counts[best] += 1;
    }
    let total = counts.iter().sum::<usize>();
    if total == 0 {
        return None;
    }
    let mut frequency = (0.5 * counts[1] as f64 + counts[2] as f64) / total as f64;
    let initial = frequency;
    for _ in 0..10 {
        let weights = [
            (1.0 - frequency) * (1.0 - frequency),
            2.0 * frequency * (1.0 - frequency),
            frequency * frequency,
        ];
        let updated = likelihoods
            .iter()
            .map(|row| {
                (row[1] * weights[1] + 2.0 * row[2] * weights[2])
                    / (row[0] * weights[0] + row[1] * weights[1] + row[2] * weights[2])
            })
            .sum::<f64>()
            / (likelihoods.len() * 2) as f64;
        let difference = (updated - frequency).abs();
        frequency = updated;
        if difference < 1e-5 {
            return Some(frequency);
        }
    }
    Some(minimize_frequency(
        if initial == frequency {
            0.5 * initial
        } else {
            initial
        },
        frequency,
        likelihoods,
    ))
}

fn genotype_frequencies(likelihoods: &[[f64; 3]], frequency: f64) -> ([f64; 3], f64) {
    let hwe = [
        (1.0 - frequency) * (1.0 - frequency),
        2.0 * frequency * (1.0 - frequency),
        frequency * frequency,
    ];
    let mut frequencies = hwe;
    for _ in 0..50 {
        let mut updated = [0.0; 3];
        for row in likelihoods {
            let values = [
                row[0] * frequencies[0],
                row[1] * frequencies[1],
                row[2] * frequencies[2],
            ];
            let normalization = values.iter().sum::<f64>() * likelihoods.len() as f64;
            for i in 0..3 {
                updated[i] += values[i] / normalization;
            }
        }
        let error = (0..3)
            .map(|i| (updated[i] - frequencies[i]).abs())
            .fold(0.0, f64::max);
        frequencies = updated;
        if error < 1e-5 {
            break;
        }
    }
    let ratio = likelihoods
        .iter()
        .map(|row| {
            row.iter()
                .zip(frequencies)
                .map(|(value, weight)| value * weight)
                .sum::<f64>()
                / row
                    .iter()
                    .zip(hwe)
                    .map(|(value, weight)| value * weight)
                    .sum::<f64>()
        })
        .product::<f64>();
    (frequencies, erfc(ratio.ln().max(0.0).sqrt()))
}

fn beta_fraction(a: f64, b: f64, x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    if x == 1.0 {
        return 1.0;
    }
    let (mut c, mut d, mut fraction) = (1.0f64, 0.0f64, 1.0f64);
    for j in 1..200 {
        let m = (j / 2) as f64;
        let coefficient = if j & 1 == 1 {
            -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
        } else {
            m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
        };
        d = (1.0 + coefficient * d).max(1e-290);
        c = (1.0 + coefficient / c).max(1e-290);
        d = 1.0 / d;
        let change = c * d;
        fraction *= change;
        if (change - 1.0).abs() < 1e-14 {
            break;
        }
    }
    (lgamma(a + b) - lgamma(a) - lgamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp() / a / fraction
}

fn beta(a: f64, b: f64, x: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        beta_fraction(a, b, x)
    } else {
        1.0 - beta_fraction(b, a, 1.0 - x)
    }
}

fn t_test(n1: f64, n2: f64, values: &[f64]) -> f64 {
    if n1 == 0.0 || n2 == 0.0 || n1 + n2 < 3.0 {
        return 1.0;
    }
    let (mean1, mean2) = (values[0] / n1, values[2] / n2);
    if mean1 <= mean2 {
        return 1.0;
    }
    let degrees = n1 + n2 - 2.0;
    let t = (mean1 - mean2)
        / (((values[1] - n1 * mean1 * mean1 + values[3] - n2 * mean2 * mean2) / degrees)
            * (1.0 / n1 + 1.0 / n2))
            .sqrt();
    if t < 0.0 {
        1.0
    } else {
        0.5 * beta(0.5 * degrees, 0.5, degrees / (degrees + t * t))
    }
}

fn float(value: f64) -> String {
    fmt_g(value as f32 as f64)
}

pub(super) fn call_record(
    rec: &VcfRecord,
    selected: &[usize],
    cfg: &CallConfig,
    hdr: &HeaderInfo,
) -> Result<Option<String>> {
    let keys: Vec<_> = rec.format.as_deref().unwrap_or("").split(':').collect();
    let pl_index = keys.iter().position(|key| *key == "PL");
    ensure!(pl_index.is_some(), "call -c requires FORMAT/PL");
    let pl_index = pl_index.unwrap();
    let ploidies = selected.len();
    ensure!(ploidies > 0, "call -c requires samples");
    let alleles: Vec<_> = std::iter::once(rec.ref_allele.as_str())
        .chain(rec.alt.split(',').filter(|allele| *allele != "."))
        .collect();
    ensure!(alleles.len() >= 2, "call -c requires at least two alleles");
    let samples: Vec<Vec<&str>> = selected
        .iter()
        .map(|sample| {
            rec.samples
                .get(*sample)
                .map(String::as_str)
                .unwrap_or(".")
                .split(':')
                .collect()
        })
        .collect();
    let genotypes = alleles.len() * (alleles.len() + 1) / 2;
    let pls: Vec<Vec<u32>> = samples
        .iter()
        .map(|sample| {
            let value = sample.get(pl_index).copied().unwrap_or(".");
            if value == "." {
                return Ok(vec![0; genotypes]);
            }
            let values = value
                .split(',')
                .map(|value| {
                    value
                        .parse::<u32>()
                        .with_context(|| format!("call -c: invalid PL value {value:?}"))
                })
                .collect::<Result<Vec<_>>>()?;
            ensure!(
                values.len() == genotypes,
                "call -c requires complete diploid PL vectors"
            );
            Ok(values)
        })
        .collect::<Result<_>>()?;
    let likelihoods: Vec<_> = pls
        .iter()
        .map(|row| {
            let minimum = *row[..3].iter().min().unwrap();
            let minimum = if minimum > 3000 { minimum } else { 0 };
            [
                10f64.powf(-((row[2] - minimum) as f64) / 10.0),
                10f64.powf(-((row[1] - minimum) as f64) / 10.0),
                10f64.powf(-((row[0] - minimum) as f64) / 10.0),
            ]
        })
        .collect();
    let indel = rec.ref_allele.len() != 1
        || alleles
            .iter()
            .skip(1)
            .any(|allele| !matches!(*allele, "<*>" | "<X>" | "<NON_REF>") && allele.len() != 1);
    let prediction = spectrum(&likelihoods, indel);
    let is_variant = prediction.reference < cfg.consensus_threshold;
    if !is_variant && cfg.variants_only {
        return Ok(None);
    }
    let mut info: IndexMap<String, String> = rec
        .info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .map(|field| {
            field
                .split_once('=')
                .map(|(key, value)| (key.to_string(), value.to_string()))
                .unwrap_or_else(|| (field.to_string(), String::new()))
        })
        .collect();
    let annotation = info
        .get("I16")
        .map(|value| {
            value
                .split(',')
                .map(|value| value.parse::<f32>().unwrap_or(0.0) as f64)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    if let Some(frequency) = estimate_frequency(&likelihoods) {
        info.insert("AF1".into(), float(1.0 - frequency));
        let (genotypes, pvalue) = genotype_frequencies(&likelihoods, frequency);
        if (0.0..=0.05).contains(&pvalue) {
            info.insert(
                "G3".into(),
                [genotypes[2], genotypes[1], genotypes[0]]
                    .map(float)
                    .join(","),
            );
            info.insert("HWE".into(), float(pvalue));
        }
    }
    info.insert("AC1".into(), prediction.count.to_string());
    if annotation.len() == 16 {
        let depth = annotation[..4].iter().sum::<f64>();
        info.insert(
            "DP4".into(),
            annotation[..4]
                .iter()
                .map(|value| (*value as i64).to_string())
                .collect::<Vec<_>>()
                .join(","),
        );
        let mq = if depth > 0.0 {
            (((annotation[9] + annotation[11]) / depth).sqrt() + 0.499) as i32
        } else {
            0
        };
        info.insert("MQ".into(), mq.to_string());
    }
    let fq = if prediction.folded_reference < 0.5 {
        -4.343 * prediction.folded_reference.ln()
    } else {
        4.343 * prediction.folded_variant.ln()
    };
    info.insert("FQ".into(), float(fq.clamp(-999.0, 999.0)));
    if annotation.len() == 16 {
        let (n1, n2) = (annotation[0] + annotation[1], annotation[2] + annotation[3]);
        if n1 > 0.0 && n2 > 0.0 {
            let mut pvalues = [1.0; 4];
            pvalues[0] = kt_fisher_exact(
                annotation[0] as i64,
                annotation[1] as i64,
                annotation[2] as i64,
                annotation[3] as i64,
            )
            .2;
            for i in 1..4 {
                pvalues[i] = t_test(n1, n2, &annotation[4 * i..4 * i + 4]);
            }
            info.insert("PV4".into(), pvalues.map(float).join(","));
        }
    }
    info.shift_remove("I16");
    info.shift_remove("QS");
    let error = if is_variant {
        prediction.reference
    } else {
        prediction.variant
    } as f32 as f64;
    let quality = if error < 1e-100 {
        999.0
    } else {
        (-4.343 * error.ln()).min(999.0)
    };
    let mut ranks: Vec<_> = (0..alleles.len())
        .map(|allele| {
            (
                pls.iter()
                    .map(|row| {
                        row.get((allele + 1) * (allele + 2) / 2 - 1)
                            .copied()
                            .unwrap_or(255) as u64
                    })
                    .sum::<u64>(),
                allele,
            )
        })
        .collect();
    ranks.sort_unstable();
    let reference_rank = ranks.iter().position(|(_, allele)| *allele == 0).unwrap();
    let mut kept = if !is_variant && !cfg.keep_alts {
        1
    } else {
        reference_rank.max(1) + 1
    };
    if cfg.keep_alts
        && alleles
            .get(kept - 1)
            .is_some_and(|allele| matches!(*allele, "<*>" | "<X>" | "<NON_REF>"))
    {
        kept -= 1;
    }
    kept = kept.min(alleles.len());
    let mut output_keys = vec!["GT"];
    output_keys.extend(
        keys.iter()
            .copied()
            .filter(|key| *key != "GT" && (!cfg.emit_gq || *key != "GQ")),
    );
    if cfg.emit_gq {
        output_keys.push("GQ");
    }
    let mut output_samples = Vec::new();
    let frequency = prediction.frequency;
    let weights = [
        (1.0 - frequency) * (1.0 - frequency),
        2.0 * frequency * (1.0 - frequency),
        frequency * frequency,
    ];
    for (i, row) in likelihoods.iter().enumerate() {
        let mut probabilities = [
            row[0] * weights[0],
            row[1] * weights[1],
            row[2] * weights[2],
        ];
        let normalization = probabilities.iter().sum::<f64>();
        for value in &mut probabilities {
            *value /= normalization;
        }
        let mut best = 0;
        for genotype in 1..3 {
            if probabilities[genotype] > probabilities[best] {
                best = genotype;
            }
        }
        if !is_variant {
            best = 2;
        }
        let genotype = ["1/1", "0/1", "0/0"][best];
        let gq = ((-4.343 * (1.0 - probabilities[best]).max(1e-308).ln() + 0.499) as i32).min(99);
        let mut values = vec![genotype.to_string()];
        for key in output_keys.iter().skip(1) {
            if *key == "GQ" && cfg.emit_gq {
                values.push(gq.to_string());
                continue;
            }
            let index = keys.iter().position(|source| source == key).unwrap();
            let value = samples[i].get(index).copied().unwrap_or(".");
            if *key == "PL" && kept < alleles.len() {
                values.push(
                    pls[i]
                        .iter()
                        .take(kept * (kept + 1) / 2)
                        .map(u32::to_string)
                        .collect::<Vec<_>>()
                        .join(","),
                );
            } else if kept < alleles.len() && hdr.format_number(key) == FieldNumber::R {
                values.push(value.split(',').take(kept).collect::<Vec<_>>().join(","));
            } else {
                values.push(value.to_string());
            }
        }
        output_samples.push(values.join(":"));
    }
    let info = info
        .iter()
        .map(|(key, value)| {
            if value.is_empty() {
                key.clone()
            } else if kept < alleles.len() && hdr.info_number(key) == FieldNumber::R {
                format!(
                    "{key}={}",
                    value.split(',').take(kept).collect::<Vec<_>>().join(",")
                )
            } else {
                format!("{key}={value}")
            }
        })
        .collect::<Vec<_>>()
        .join(";");
    Ok(Some(format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.chrom,
        rec.pos,
        rec.id,
        rec.ref_allele,
        if kept == 1 {
            ".".into()
        } else {
            alleles[1..kept].join(",")
        },
        float(quality),
        rec.filter,
        info,
        output_keys.join(":"),
        output_samples.join("\t")
    )))
}

pub(super) fn header_lines() -> &'static [&'static str] {
    &[
        "##INFO=<ID=AF1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first ALT allele frequency (assuming HWE)\">",
        "##INFO=<ID=AF2,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)\">",
        "##INFO=<ID=AC1,Number=1,Type=Float,Description=\"Max-likelihood estimate of the first ALT allele count (no HWE assumption)\">",
        "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Root-mean-square mapping quality of covering reads\">",
        "##INFO=<ID=FQ,Number=1,Type=Float,Description=\"Phred probability of all samples being the same\">",
        "##INFO=<ID=PV4,Number=4,Type=Float,Description=\"P-values for strand bias, baseQ bias, mapQ bias and tail distance bias\">",
        "##INFO=<ID=G3,Number=3,Type=Float,Description=\"ML estimate of genotype frequencies\">",
        "##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Chi^2 based HWE test P-value based on G3\">",
        "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases\">",
    ]
}

#[cfg(test)]
#[path = "../../tests/unit/call_consensus.rs"]
mod tests;
