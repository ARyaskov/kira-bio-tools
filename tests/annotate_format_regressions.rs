#![cfg(feature = "cli")]

use std::path::Path;
use std::process::Command;

#[test]
fn annotate_format_retains_cardinality_ploidy_and_implicit_contigs() {
    let dir = tempfile::tempdir().unwrap();
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/annotate/test18");
    for name in ["in.vcf.gz", "in.vcf.gz.tbi", "db.vcf.gz", "db.vcf.gz.tbi"] {
        std::fs::copy(fixture.join(name), dir.path().join(name)).unwrap();
    }
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bt"))
        .current_dir(dir.path())
        .args([
            "annotate",
            "-a",
            "db.vcf.gz",
            "-c",
            "FMT",
            "in.vcf.gz",
            "-o",
            "out.vcf",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(dir.path().join("out.vcf")).unwrap();
    assert!(text.lines().any(|line| line == "##contig=<ID=X>"));
    let rows: Vec<Vec<&str>> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 6);
    assert_eq!(rows[2][9], "0/1:4,.,5:4,.,5:4,.,.,5,.,6:4,.,.,5,.,6");
    assert_eq!(rows[4][10], "0:1,3:1,3:1,3:1,3");
    assert_eq!(rows[5][10], "0:1,.,2:1,.,2:1,.,2:1,.,2");
}
