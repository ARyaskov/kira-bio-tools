#![cfg(feature = "cli")]

use std::path::Path;
use std::process::Command;

fn run(input: &Path, output: &Path, extra: &[&str]) {
    let result = Command::new(env!("CARGO_BIN_EXE_kira-bt"))
        .args(["cnv", "--", "-s", "Q", "-c", "C", "-o"])
        .arg(output)
        .args(extra)
        .arg(input)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

fn records(directory: &Path, chromosome: Option<&str>) -> Vec<String> {
    let mut result = Vec::new();
    for name in [
        "summary.tab",
        "summary.Q.tab",
        "summary.C.tab",
        "cn.Q.tab",
        "cn.C.tab",
        "dat.Q.tab",
        "dat.C.tab",
    ] {
        result.push(format!("###{name}"));
        let text = std::fs::read_to_string(directory.join(name)).unwrap();
        result.extend(
            text.lines()
                .filter(|line| {
                    !line.starts_with('#')
                        && chromosome.is_none_or(|chr| line.split('\t').next() == Some(chr))
                })
                .map(str::to_string),
        );
    }
    result
}

#[test]
fn short_default_window_matches_the_bounded_bcftools_reference() {
    let work = tempfile::tempdir().unwrap();
    let input = work.path().join("input.vcf");
    std::fs::write(&input, include_str!("cnv/test6/in.vcf")).unwrap();
    let output = work.path().join("out");
    run(&input, &output, &["-t", "2:100-400"]);
    assert_eq!(
        records(&output, None),
        include_str!("fixtures/cnv/short-default.records.tsv")
            .lines()
            .map(str::to_string)
            .collect::<Vec<_>>()
    );
}

#[test]
fn short_chromosomes_are_independent_of_previous_buffers() {
    let work = tempfile::tempdir().unwrap();
    let source = include_str!("cnv/test6/in.vcf");
    let single = work.path().join("single.vcf");
    std::fs::write(
        &single,
        source
            .lines()
            .filter(|line| line.starts_with('#') || line.starts_with("2\t"))
            .collect::<Vec<_>>()
            .join("\n")
            + "\n",
    )
    .unwrap();
    let prefixed = work.path().join("prefixed.vcf");
    let rows: Vec<_> = source
        .lines()
        .map(|line| {
            if !line.starts_with("1\t") {
                return line.to_string();
            }
            let mut fields: Vec<_> = line.split('\t').map(str::to_string).collect();
            fields[9] = "0/1:0.50:1000".to_string();
            fields.join("\t")
        })
        .collect();
    std::fs::write(&prefixed, rows.join("\n") + "\n").unwrap();
    for (index, extra) in [vec![], vec!["-O", "0.8"]].iter().enumerate() {
        let alone = work.path().join(format!("alone-{index}"));
        let after_prefix = work.path().join(format!("prefix-{index}"));
        run(&single, &alone, extra);
        run(&prefixed, &after_prefix, extra);
        assert_eq!(
            records(&alone, Some("2")),
            records(&after_prefix, Some("2"))
        );
    }
}

#[test]
fn a_single_site_accepts_an_oversized_window() {
    let work = tempfile::tempdir().unwrap();
    let source = include_str!("cnv/test6/in.vcf");
    let input = work.path().join("single.vcf");
    std::fs::write(
        &input,
        source
            .lines()
            .filter(|line| line.starts_with('#') || line.starts_with("2\t100\t"))
            .collect::<Vec<_>>()
            .join("\n")
            + "\n",
    )
    .unwrap();
    let output = work.path().join("out");
    run(&input, &output, &["-L", "9223372036854775807"]);
    let data = std::fs::read_to_string(output.join("cn.Q.tab")).unwrap();
    let rows: Vec<_> = data.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0], "2\t100\t2\t0.145833\t0.145833\t0.562500\t0.145833");
}
