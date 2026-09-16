//! `getopt_long`-style parsing of bcftools option tokens passed through
//! verbatim (after `--` or as a raw argument list): `-Xval`, `-X val`,
//! bundled flags, `--long val`, `--long=val` and unambiguous long-option
//! prefixes.

use anyhow::{Result, bail};

#[derive(Clone, Copy, Debug)]
pub struct OptSpec {
    pub short: Option<char>,
    pub long: &'static str,
    pub has_arg: bool,
}

pub const fn opt(short: Option<char>, long: &'static str, has_arg: bool) -> OptSpec {
    OptSpec {
        short,
        long,
        has_arg,
    }
}

/// Parsed options in order of appearance (by long name) and the positionals.
pub struct Parsed {
    pub opts: Vec<(&'static str, Option<String>)>,
    pub positional: Vec<String>,
}

fn find_long<'a>(specs: &'a [OptSpec], name: &str) -> Result<&'a OptSpec> {
    if let Some(s) = specs.iter().find(|s| s.long == name) {
        return Ok(s);
    }
    let cands: Vec<&OptSpec> = specs.iter().filter(|s| s.long.starts_with(name)).collect();
    match cands.len() {
        1 => Ok(cands[0]),
        0 => bail!("Unknown option: --{name}"),
        _ => bail!("Option --{name} is ambiguous"),
    }
}

pub fn parse(tokens: &[String], specs: &[OptSpec]) -> Result<Parsed> {
    let mut out = Parsed {
        opts: Vec::new(),
        positional: Vec::new(),
    };
    let mut i = 0;
    let mut only_positional = false;
    while i < tokens.len() {
        let t = &tokens[i];
        i += 1;
        if only_positional || t == "-" || !t.starts_with('-') {
            out.positional.push(t.clone());
            continue;
        }
        if t == "--" {
            only_positional = true;
            continue;
        }
        if let Some(rest) = t.strip_prefix("--") {
            let (name, inline_val) = match rest.split_once('=') {
                Some((n, v)) => (n, Some(v.to_string())),
                None => (rest, None),
            };
            let spec = find_long(specs, name)?;
            if spec.has_arg {
                let v = match inline_val {
                    Some(v) => v,
                    None => {
                        if i >= tokens.len() {
                            bail!("Option --{} requires an argument", spec.long);
                        }
                        i += 1;
                        tokens[i - 1].clone()
                    }
                };
                out.opts.push((spec.long, Some(v)));
            } else {
                if inline_val.is_some() {
                    bail!("Option --{} does not take an argument", spec.long);
                }
                out.opts.push((spec.long, None));
            }
            continue;
        }
        let chars: Vec<char> = t[1..].chars().collect();
        let mut k = 0;
        while k < chars.len() {
            let c = chars[k];
            let Some(spec) = specs.iter().find(|s| s.short == Some(c)) else {
                bail!("Unknown option: -{c}");
            };
            k += 1;
            if spec.has_arg {
                let v = if k < chars.len() {
                    let s: String = chars[k..].iter().collect();
                    k = chars.len();
                    s
                } else {
                    if i >= tokens.len() {
                        bail!("Option -{c} requires an argument");
                    }
                    i += 1;
                    tokens[i - 1].clone()
                };
                out.opts.push((spec.long, Some(v)));
            } else {
                out.opts.push((spec.long, None));
            }
        }
    }
    Ok(out)
}

/// `strtod`-style parse of a whole token.
pub fn parse_f64(v: &str, what: &str) -> Result<f64> {
    v.trim()
        .parse::<f64>()
        .map_err(|_| anyhow::anyhow!("Could not parse: {what} {v}"))
}

pub fn parse_i64(v: &str, what: &str) -> Result<i64> {
    v.trim()
        .parse::<i64>()
        .map_err(|_| anyhow::anyhow!("Could not parse: {what} {v}"))
}

#[cfg(test)]
#[path = "../../tests/unit/cli_getopt.rs"]
mod tests;
