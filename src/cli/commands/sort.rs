use anyhow::Result;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::PathBuf;

use crate::cli::args::SortArgs;

pub fn cmd_sort(args: SortArgs) -> Result<()> {
    let _ = &args.bcftools_args;
    let out_path = args
        .output
        .clone()
        .unwrap_or_else(|| PathBuf::from("out.sort.vcf"));

    let mut header = Vec::new();
    let mut records = Vec::new();
    let mut contig_order = std::collections::HashMap::<String, usize>::new();
    let mut has_pass_filter_header = false;

    let reader: Box<dyn Read> = if matches!(
        args.input.extension().and_then(|x| x.to_str()),
        Some("gz" | "bgz" | "bgzf")
    ) {
        Box::new(MultiGzDecoder::new(File::open(&args.input)?))
    } else {
        Box::new(File::open(&args.input)?)
    };

    for (idx, line) in BufReader::new(reader).lines().enumerate() {
        let line = line?;
        if line.starts_with('#') {
            if let Some(id) = parse_contig_id(&line) {
                let next = contig_order.len();
                contig_order.entry(id).or_insert(next);
            }
            if line.starts_with("##FILTER=<ID=PASS,") {
                has_pass_filter_header = true;
            }
            header.push(line);
            continue;
        }
        let mut parts = line.split('\t');
        let chrom = parts.next().unwrap_or("").to_string();
        let pos = parts
            .next()
            .and_then(|x| x.parse::<u32>().ok())
            .unwrap_or(u32::MAX);
        let _id = parts.next().unwrap_or("").to_string();
        let r = parts.next().unwrap_or("").to_string();
        let a = parts.next().unwrap_or("").to_string();
        records.push((chrom, pos, r, a, idx, line));
    }

    if !has_pass_filter_header {
        let insert_at = header
            .iter()
            .position(|h| h.starts_with("##reference="))
            .unwrap_or(1);
        header.insert(
            insert_at,
            "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
        );
    }

    records.sort_by(|a, b| {
        let ao = contig_order
            .get(&a.0)
            .copied()
            .unwrap_or(usize::MAX.saturating_sub(1));
        let bo = contig_order
            .get(&b.0)
            .copied()
            .unwrap_or(usize::MAX.saturating_sub(1));
        ao.cmp(&bo)
            .then_with(|| a.0.cmp(&b.0))
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
            .then_with(|| a.3.cmp(&b.3))
            .then_with(|| a.4.cmp(&b.4))
    });

    let mut out = File::create(out_path)?;
    for h in &header {
        writeln!(out, "{h}")?;
    }
    for (_, _, _, _, _, line) in &records {
        writeln!(out, "{line}")?;
    }
    Ok(())
}

fn parse_contig_id(header_line: &str) -> Option<String> {
    if !header_line.starts_with("##contig=<") {
        return None;
    }
    let body = header_line
        .trim_start_matches("##contig=<")
        .trim_end_matches('>');
    for kv in body.split(',') {
        if let Some(id) = kv.strip_prefix("ID=") {
            return Some(id.to_string());
        }
    }
    None
}
