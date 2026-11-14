#![cfg(feature = "cli")]

use std::path::Path;
use std::process::Command;

fn write_vcf(path: &Path, version: &str, sample: &str, row: &str) {
    std::fs::write(path, format!("##fileformat=VCFv{version}\n##contig=<ID=1,length=1000>\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n{row}\n")).unwrap();
}

#[test]
fn merge_keeps_structural_extents_vcf_version_and_missing_phase() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.vcf");
    let b = dir.path().join("b.vcf");
    write_vcf(
        &a,
        "4.1",
        "A",
        "1\t10\tSV1\tN\t<DEL>\t.\t.\tEND=100\tGT\t0|.",
    );
    write_vcf(
        &b,
        "4.4",
        "B",
        "1\t10\tSV2\tN\t<DEL>\t.\t.\tEND=200\tGT\t1/1",
    );
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bt"))
        .args([
            "merge",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--",
            "--no-version",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.starts_with("##fileformat=VCFv4.4\n"));
    assert_eq!(text.matches("##fileformat=").count(), 1);
    let records: Vec<_> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(records.len(), 2);
    assert!(records[0].ends_with("END=100\tGT\t0|.\t./."));
    assert!(records[1].ends_with("END=200\tGT\t./.\t1/1"));
}

#[test]
fn merge_none_pairs_mixed_records_with_a_shared_allele() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.vcf");
    let b = dir.path().join("b.vcf");
    write_vcf(&a, "4.2", "A", "1\t10\t.\tA\tC,<NON_REF>\t.\t.\t.\tGT\t0|1");
    write_vcf(
        &b,
        "4.2",
        "B",
        "1\t10\t.\tA\tACCC,C,<NON_REF>\t.\t.\t.\tGT\t0/2",
    );
    let output = Command::new(env!("CARGO_BIN_EXE_kira-bt"))
        .args([
            "merge",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--",
            "--merge",
            "none",
            "--no-version",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    let records: Vec<_> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(
        records,
        ["1\t10\t.\tA\tC,<NON_REF>,ACCC\t.\t.\t.\tGT\t0|1\t0/1"]
    );
}
