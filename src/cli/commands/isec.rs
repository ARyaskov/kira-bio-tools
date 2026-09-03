use crate::annotate::postproc::RegionFilter;
use crate::cli::args::IsecArgs;
use crate::filter::FilterEngine;
use crate::vcf::header::ContigDict;
use crate::vcf::sink::{OutputKind, parse_output_type};
use crate::vcf::variant_type::{VT_INDEL, VT_SNP, record_type};
use crate::vcf::{UnifiedVcfReader, VcfSink, parse_vcf_line};
use anyhow::{Context, Result, bail};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Collapse { None, All, Snps, Indels, Both, Some_, Id }

fn parse_collapse(s: Option<&str>) -> Result<Collapse> {
    let Some(s) = s else { return Ok(Collapse::None); };
    match s {
        "none" => Ok(Collapse::None),
        "all" | "any" => Ok(Collapse::All),
        "snps" => Ok(Collapse::Snps),
        "indels" => Ok(Collapse::Indels),
        "both" => Ok(Collapse::Both),
        "some" => Ok(Collapse::Some_),
        "id" => Ok(Collapse::Id),
        _ => bail!("-c: unknown collapse mode '{}', expected none|all|snps|indels|both|some|id", s),
    }
}

struct Rec {
    line: String,
    id: String,
    refa: String,
    alt: String,
    vt: u32,
}

impl Rec {
    fn from_line(line: String) -> Option<Self> {
        let mut it = line.splitn(6, '\t');
        let _chrom = it.next()?;
        let _pos = it.next()?;
        let id = it.next()?.to_string();
        let refa = it.next()?.to_string();
        let alt = it.next()?.to_string();
        let vt = record_type(&refa, &alt);
        Some(Self { line, id, refa, alt, vt })
    }

    fn alts(&self) -> Vec<&str> {
        self.alt.split(',').filter(|a| *a != "*" && *a != ".").collect()
    }
}

/// Do two records at the same position count as the same site?
fn same_site(a: &Rec, b: &Rec, mode: Collapse) -> bool {
    match mode {
        Collapse::All => true,
        Collapse::None => a.refa.eq_ignore_ascii_case(&b.refa) && a.alt.eq_ignore_ascii_case(&b.alt),
        Collapse::Snps => {
            let s = a.vt & VT_SNP != 0 && b.vt & VT_SNP != 0;
            s || (a.refa.eq_ignore_ascii_case(&b.refa) && a.alt.eq_ignore_ascii_case(&b.alt))
        }
        Collapse::Indels => {
            let i = a.vt & VT_INDEL != 0 && b.vt & VT_INDEL != 0;
            i || (a.refa.eq_ignore_ascii_case(&b.refa) && a.alt.eq_ignore_ascii_case(&b.alt))
        }
        Collapse::Both => {
            (a.vt & VT_SNP != 0 && b.vt & VT_SNP != 0)
                || (a.vt & VT_INDEL != 0 && b.vt & VT_INDEL != 0)
                || (a.refa.eq_ignore_ascii_case(&b.refa) && a.alt.eq_ignore_ascii_case(&b.alt))
        }
        Collapse::Some_ => {
            a.refa.eq_ignore_ascii_case(&b.refa) && a.alts().iter().any(|x| b.alts().iter().any(|y| x.eq_ignore_ascii_case(y)))
        }
        Collapse::Id => a.id != "." && a.id == b.id,
    }
}

struct Source {
    reader: UnifiedVcfReader,
    headers: Vec<String>,
    next: Option<(usize, u32, String, Rec)>,
    region: Option<RegionFilter>,
    regions_overlap: u8,
    target: Option<RegionFilter>,
    target_inverse: bool,
    apply_filters: Option<Vec<String>>,
    include: Option<FilterEngine>,
    exclude: Option<FilterEngine>,
}

impl Source {
    fn advance(&mut self, contigs: &mut ContigDict) -> Result<()> {
        loop {
            let Some(line) = self.reader.read_line()? else {
                self.next = None;
                return Ok(());
            };
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            if let Some(rf) = &self.region {
                if !rf.line_passes_mode(&line, self.regions_overlap) { continue; }
            }
            if let Some(tf) = &self.target {
                if tf.line_passes_mode(&line, 0) == self.target_inverse { continue; }
            }
            if let Some(af) = &self.apply_filters {
                let f = line.splitn(8, '\t').nth(6).unwrap_or(".");
                let pass = if f == "." || f.is_empty() { af.iter().any(|a| a == ".") } else { f.split(';').any(|t| af.iter().any(|a| a == t)) };
                if !pass { continue; }
            }
            if self.include.is_some() || self.exclude.is_some() {
                if let Some(rec) = parse_vcf_line(&line) {
                    if let Some(e) = &self.include {
                        if !e.eval(&rec).map(|r| r.pass_site).unwrap_or(true) { continue; }
                    }
                    if let Some(e) = &self.exclude {
                        if e.eval(&rec).map(|r| r.pass_site).unwrap_or(false) { continue; }
                    }
                }
            }
            let mut it = line.splitn(3, '\t');
            let chrom = it.next().unwrap_or("").to_string();
            let pos: u32 = it
                .next()
                .and_then(|p| p.parse().ok())
                .ok_or_else(|| anyhow::anyhow!("isec: malformed POS in record {line:?}"))?;
            let rank = contigs.insert(&chrom) as usize;
            let Some(rec) = Rec::from_line(line) else { continue };
            self.next = Some((rank, pos, chrom, rec));
            return Ok(());
        }
    }
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
    let has_targets = args.targets.is_some() || args.targets_file.is_some();
    if inputs.is_empty() || (inputs.len() < 2 && !has_targets) {
        bail!("isec: need at least 2 inputs (or one input with -t/-T)");
    }
    // A single file with targets behaves like intersecting with the target set.
    let n_files = inputs.len();

    let nfiles_spec = parse_nfiles(args.nfiles.as_deref(), n_files, args.complement)?;
    let mut write_files = parse_write(args.write.as_deref(), n_files)?;
    // A single input without -p streams its matching records as VCF (like -w 1).
    if write_files.is_none() && n_files == 1 && args.prefix.is_none() {
        write_files = Some(vec![0]);
    }
    let collapse = parse_collapse(args.collapse.as_deref())?;
    let apply_filters: Option<Vec<String>> = args.apply_filters.as_deref().map(|s| s.split(',').map(|t| t.trim().to_string()).collect());
    let region = if let Some(s) = &args.regions {
        Some(RegionFilter::from_cli(s)?)
    } else if let Some(p) = &args.regions_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };
    let target = if let Some(s) = &args.targets {
        Some(RegionFilter::from_cli(s.trim_start_matches('^'))?)
    } else if let Some(p) = &args.targets_file {
        Some(RegionFilter::from_file(p)?)
    } else {
        None
    };
    let target_inverse = args.targets.as_deref().is_some_and(|s| s.starts_with('^'));
    let kind = args.output_type.as_deref().map(parse_output_type).transpose()?.unwrap_or(OutputKind::Vcf);

    let mut contigs = ContigDict::new();
    let mut sources: Vec<Source> = Vec::with_capacity(n_files);
    for p in &inputs {
        let r = UnifiedVcfReader::open(p).with_context(|| format!("open {:?}", p))?;
        let headers = r.header()?;
        for h in &headers {
            if let Some((id, len)) = crate::vcf::header::parse_contig_line(h) {
                contigs.insert_with_length(&id, len);
            }
        }
        let include = args.include.as_deref().map(|e| FilterEngine::new(&headers, Some(e), false)).transpose().context("-i")?;
        let exclude = args.exclude.as_deref().map(|e| FilterEngine::new(&headers, Some(e), false)).transpose().context("-e")?;
        let mut s = Source {
            reader: r,
            headers,
            next: None,
            region: region.clone(),
            regions_overlap: args.regions_overlap,
            target: target.clone(),
            target_inverse,
            apply_filters: apply_filters.clone(),
            include,
            exclude,
        };
        s.advance(&mut contigs)?;
        sources.push(s);
    }

    // Output destinations.
    let dst = args.prefix.clone();
    let mut writers: Vec<Option<VcfSink>> = Vec::with_capacity(n_files);
    let mut sites_out: Option<BufWriter<File>> = None;
    let mut stdout_sink: Option<VcfSink> = None;
    let mut stdout_sites: Option<BufWriter<std::io::Stdout>> = None;
    if let Some(d) = &dst {
        fs::create_dir_all(d).with_context(|| format!("create dir {:?}", d))?;
        let ext = match kind {
            OutputKind::Vcf => "vcf",
            OutputKind::VcfGz(_) => "vcf.gz",
            OutputKind::Bcf(_) => "bcf",
        };
        let mut readme = BufWriter::new(File::create(d.join("README.txt"))?);
        writeln!(readme, "This file was produced by vcfisec.")?;
        let cmd: Vec<String> = std::env::args().collect();
        writeln!(readme, "The command line was:\t{}\n", cmd.join(" "))?;
        writeln!(readme, "Using the following file names:")?;
        for (i, p) in inputs.iter().enumerate() {
            let selected = write_files.as_ref().is_none_or(|w| w.contains(&i));
            if !selected {
                writers.push(None);
                continue;
            }
            let path = d.join(format!("{:04}.{ext}", i));
            let what = if args.complement && i == 0 { "for records private to".to_string() } else { "for stripped".to_string() };
            writeln!(readme, "{}\t{}\t{}", path.display(), what, p.display())?;
            let mut w = VcfSink::open(Some(&path), kind, &sources[i].headers)?;
            w.write_header(&sources[i].headers)?;
            writers.push(Some(w));
        }
        readme.flush()?;
        sites_out = Some(BufWriter::with_capacity(1 << 20, File::create(d.join("sites.txt"))?));
    } else if let Some(w) = &write_files {
        // Without -p, -w selects one file whose records are printed.
        let i = *w.first().unwrap_or(&0);
        let mut s = VcfSink::open(args.output.as_deref(), kind, &sources[i].headers)?;
        s.write_header(&sources[i].headers)?;
        stdout_sink = Some(s);
    } else {
        stdout_sites = Some(BufWriter::with_capacity(1 << 20, std::io::stdout()));
    }
    let stdout_file = write_files.as_ref().and_then(|w| w.first().copied()).unwrap_or(0);

    // k-way walk over the sorted inputs.
    loop {
        let Some((rank, pos)) = sources.iter().filter_map(|s| s.next.as_ref().map(|n| (n.0, n.1))).min() else { break };
        // All records at this site from every file.
        let mut here: Vec<Vec<Rec>> = (0..n_files).map(|_| Vec::new()).collect();
        let mut chrom = String::new();
        for (i, s) in sources.iter_mut().enumerate() {
            while s.next.as_ref().is_some_and(|n| n.0 == rank && n.1 == pos) {
                let (_, _, c, rec) = s.next.take().unwrap();
                chrom = c;
                here[i].push(rec);
                s.advance(&mut contigs)?;
            }
        }
        // Presence mask per record: which files have a matching record.
        let mut masks: Vec<Vec<Vec<bool>>> = vec![Vec::new(); n_files];
        for i in 0..n_files {
            for r in &here[i] {
                let mut m = vec![false; n_files];
                for j in 0..n_files {
                    m[j] = here[j].iter().any(|o| same_site(r, o, collapse));
                }
                masks[i].push(m);
            }
        }
        // With targets and a single file, the target hit counts as the second set.
        // sites.txt: one line per distinct site, from the first file carrying it.
        let mut emitted: Vec<(usize, &Rec)> = Vec::new();
        for i in 0..n_files {
            for (k, r) in here[i].iter().enumerate() {
                let m = &masks[i][k];
                let pass = matches_nfiles(m, &nfiles_spec);
                if pass {
                    if let Some(w) = writers.get_mut(i).and_then(|w| w.as_mut()) {
                        w.write_line(&r.line)?;
                    }
                    if i == stdout_file {
                        if let Some(s) = stdout_sink.as_mut() {
                            s.write_line(&r.line)?;
                        }
                    }
                }
                // Every record of a file is listed; a later file's record is
                // skipped only when an earlier file already carried that site.
                if !pass || emitted.iter().any(|(f, e)| *f != i && same_site(e, r, collapse)) {
                    continue;
                }
                emitted.push((i, r));
                let bits: String = m.iter().map(|b| if *b { '1' } else { '0' }).collect();
                let line = format!("{}\t{}\t{}\t{}\t{}", chrom, pos, r.refa, r.alt, bits);
                if let Some(w) = sites_out.as_mut() {
                    writeln!(w, "{line}")?;
                }
                if let Some(w) = stdout_sites.as_mut() {
                    writeln!(w, "{line}")?;
                }
            }
        }
    }

    for w in writers.into_iter().flatten() {
        w.finish()?;
    }
    if let Some(w) = sites_out.as_mut() { w.flush()?; }
    if let Some(s) = stdout_sink { s.finish()?; }
    if let Some(w) = stdout_sites.as_mut() { w.flush()?; }
    Ok(())
}

#[derive(Debug)]
enum NSpec { Exact(usize), AtLeast(usize), AtMost(usize), Mask(Vec<Option<bool>>) }

fn parse_nfiles(s: Option<&str>, n: usize, complement: bool) -> Result<NSpec> {
    if complement {
        // -C: records private to the first file.
        let mut m = vec![Some(false); n];
        m[0] = Some(true);
        return Ok(NSpec::Mask(m));
    }
    let Some(s) = s else { return Ok(NSpec::AtLeast(1)); };
    if let Some(rest) = s.strip_prefix('=') {
        return Ok(NSpec::Exact(rest.parse().context("-n =N: parse N")?));
    }
    if let Some(rest) = s.strip_prefix('+') {
        return Ok(NSpec::AtLeast(rest.parse().context("-n +N: parse N")?));
    }
    if let Some(rest) = s.strip_prefix('-') {
        return Ok(NSpec::AtMost(rest.parse().context("-n -N: parse N")?));
    }
    if let Some(rest) = s.strip_prefix('~') {
        let bits: Vec<Option<bool>> = rest.chars().map(|c| match c { '1' => Some(true), '0' => Some(false), _ => None }).collect();
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
        NSpec::AtMost(k) => n_present <= *k,
        NSpec::Mask(mask) => presence.iter().zip(mask.iter()).all(|(p, m)| m.is_none_or(|m| m == *p)),
    }
}

fn parse_write(s: Option<&str>, n: usize) -> Result<Option<Vec<usize>>> {
    let Some(s) = s else { return Ok(None); };
    let mut out = Vec::new();
    for tok in s.split(',') {
        let i: usize = tok.trim().parse().context("-w: parse index")?;
        if i == 0 || i > n { bail!("-w: index {i} out of range 1..={n}"); }
        out.push(i - 1);
    }
    Ok(Some(out))
}

#[cfg(test)]
#[path = "../../../tests/unit/cli_commands_isec.rs"]
mod tests;
