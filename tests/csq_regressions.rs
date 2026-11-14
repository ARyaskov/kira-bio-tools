#![cfg(feature = "cli")]

use std::path::Path;
use std::process::Command;

#[test]
fn complex_alleles_have_one_feature_consequence_and_correct_sample_bits() {
    let dir = tempfile::tempdir().unwrap();
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/csq/test16");
    for name in ["ref.fa", "ann.gff", "in.vcf"] {
        std::fs::copy(fixture.join(name), dir.path().join(name)).unwrap();
    }
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bt"))
        .current_dir(dir.path())
        .args([
            "csq",
            "--",
            "-f",
            "ref.fa",
            "-g",
            "ann.gff",
            "--unify-chr-names",
            "1",
            "in.vcf",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
        .collect();
    let complex = rows
        .iter()
        .find(|row| row[0] == "1" && row[1] == "140")
        .unwrap();
    assert_eq!(
        complex[7]
            .matches("splice_region&feature_elongation")
            .count(),
        1
    );
    assert!(
        complex[7].contains("feature_elongation|XYZ|ENST00000000001|protein_coding|+|8LR>7QRR")
    );
    assert_eq!(complex[9], "1|0:5");
    for pos in ["113", "152"] {
        let row = rows
            .iter()
            .find(|row| row[0] == "3" && row[1] == pos)
            .unwrap();
        assert_eq!(row[9], "1|0:1");
    }
}
