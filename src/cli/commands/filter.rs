use anyhow::Result;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::cli::args::FilterArgs;
use crate::filter::{eval_ast, parse_expr};
use crate::vcf::{parse_vcf_full_line, VcfReader};

pub fn cmd_filter(args: &FilterArgs) -> Result<()> {
    let ast = parse_expr(&args.expr)?;
    let mut reader = VcfReader::open(&args.input)?;
    let mut out = BufWriter::new(File::create(&args.output)?);

    for h in &reader.header()? {
        writeln!(out, "{}", h)?;
    }

    let soft = args.soft_filter.clone();
    let pass_only = args.pass_only;

    while let Some((line, _offset)) = reader.next_raw_line()? {
        if let Some(mut rec) = parse_vcf_full_line(&line) {
            let pass = eval_ast(&ast, &rec);

            match (pass, &soft, pass_only) {
                (true, _, false) => rec.filter = "PASS".to_string(),
                (true, _, true) => rec.filter = "PASS".to_string(),
                (false, None, _) => continue,
                (false, Some(name), false) => rec.filter = name.clone(),
                (false, Some(_), true) => continue,
            }

            writeln!(out, "{}", rec.to_line())?;
        }
    }

    Ok(())
}
