#![cfg(feature = "cli")]

use std::process::Command;

fn normalized(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|line| !line.starts_with("##bcftools") && !line.starts_with("##kira_bt"))
        .collect()
}

#[test]
fn consensus_calls_match_reference_posteriors_genotypes_and_quality() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    std::fs::write(&input, include_str!("fixtures/call-consensus/input.vcf")).unwrap();
    for (flags, expected) in [
        (
            vec!["-c"],
            include_str!("fixtures/call-consensus/reference.vcf"),
        ),
        (
            vec!["-c", "-f", "GQ"],
            include_str!("fixtures/call-consensus/gq.vcf"),
        ),
        (
            vec!["-cA"],
            include_str!("fixtures/call-consensus/keep-alts.vcf"),
        ),
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_kira-bt"))
            .arg("call")
            .arg(&input)
            .arg("--")
            .args(&flags)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let actual = String::from_utf8(output.stdout).unwrap();
        assert_eq!(normalized(&actual), normalized(expected), "{flags:?}");
    }
}

#[test]
fn consensus_rejects_invalid_likelihoods() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("invalid.vcf");
    let source = include_str!("fixtures/call-consensus/input.vcf");
    let mut lines: Vec<_> = source
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(str::to_string)
        .collect();
    let row = source.lines().find(|line| !line.starts_with('#')).unwrap();
    let mut fields: Vec<_> = row.split('\t').map(str::to_string).collect();
    fields[9] = "1,invalid,0:5:0".to_string();
    lines.push(fields.join("\t"));
    std::fs::write(&input, lines.join("\n")).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bt"))
        .arg("call")
        .arg(&input)
        .args(["--", "-c"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid PL"));
}
