#![cfg(feature = "cli")]

use std::process::Command;

#[test]
fn indexed_mpileup_merges_repeated_overlapping_regions() {
    let bam =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/mpileup/test2/mpileup.1.bam");
    let run = |regions: &str| {
        let output = Command::new(env!("CARGO_BIN_EXE_kira-bt"))
            .arg("mpileup")
            .arg(&bam)
            .args(["--", "-B", "-I", "-r", regions])
            .env("KIRA_BT_NO_PROGRESS", "1")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_string)
            .collect::<Vec<_>>()
    };
    let expected = run("17:100-105");
    assert_eq!(expected.len(), 6);
    let actual = run("17:104-105,17:102-103,17:103-104,17:100-102,17:100-105");
    assert_eq!(actual, expected);
}
