#![cfg(feature = "cli")]

use std::path::Path;
use std::process::Command;

fn check(case: &str, files: &[&str], flags: &[&str], expected: &str, all_records: bool) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/mpileup")
        .join(case);
    let work = tempfile::tempdir().unwrap();
    let reference = work.path().join("ref.fa");
    std::fs::copy(root.join("ref.fa"), &reference).unwrap();
    let mut command = Command::new(env!("CARGO_BIN_EXE_kira-bt"));
    command.arg("mpileup");
    for file in files {
        command.arg(root.join(file));
    }
    let output = command
        .arg("--")
        .arg("-f")
        .arg(&reference)
        .args(flags)
        .env("KIRA_BT_NO_PROGRESS", "1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let actual = String::from_utf8(output.stdout).unwrap();
    let records: Vec<_> = actual
        .lines()
        .filter(|line| !line.starts_with('#') && (all_records || line.contains("INDEL;")))
        .collect();
    assert_eq!(
        records,
        expected.lines().collect::<Vec<_>>(),
        "{case} {flags:?}"
    );
}

#[test]
fn ambiguous_indel_reads_match_bcftools_depth_modes() {
    for (mode, expected) in [
        ("drop", include_str!("fixtures/mpileup/test17.records.tsv")),
        ("incAD", include_str!("fixtures/mpileup/test18.records.tsv")),
        (
            "incAD0",
            include_str!("fixtures/mpileup/test19.records.tsv"),
        ),
    ] {
        check(
            "test17",
            &["indel-AD.2.bam"],
            &["-a", "AD", "-r", "11:75", "--ambig-reads", mode],
            expected,
            false,
        );
    }
}

#[test]
fn mixed_indel_likelihoods_and_persistent_biases_match_bcftools() {
    check(
        "test4",
        &["mpileup.1.bam", "mpileup.2.bam", "mpileup.3.bam"],
        &["-r", "17:100-600", "-a", "DP,DV,-AD"],
        include_str!("fixtures/mpileup/test4.records.tsv"),
        false,
    );
    check(
        "test16",
        &["indel-AD.1.bam"],
        &["-a", "AD"],
        include_str!("fixtures/mpileup/test16.records.tsv"),
        false,
    );
}

#[test]
fn indexed_repeat_indels_keep_next_mate_context_and_multiple_alts() {
    check(
        "test25",
        &["annot-NMBZ.3.1.bam", "annot-NMBZ.3.2.bam"],
        &["-r", "chr16:75", "-a", "-AD,INFO/NMBZ"],
        include_str!("fixtures/mpileup/test25.records.tsv"),
        true,
    );
}
