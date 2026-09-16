use super::*;
use crate::vcf::{VcfRecordRef, parse_vcf_line};

fn headers() -> Vec<String> {
    vec![
        "##fileformat=VCFv4.2".into(),
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">".into(),
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"\">".into(),
        "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"\">".into(),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">".into(),
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"\">".into(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2".into(),
    ]
}

const LINES: [&str; 4] = [
    "1\t100\trs1\tA\tG\t50\tPASS\tDP=10;AF=0.5;DB\tGT:DP\t0/1:7\t1/1:3",
    "1\t200\t.\tAT\tA\t3\tq10\tDP=2;AF=0.01\tGT:DP\t0/0:1\t./.:.",
    "1\t300\t.\tC\tT,G\t.\t.\tDP=40;AF=0.2,0.3\tGT:DP\t1/2:20\t0/1:20",
    "2\t5\t.\tG\tC\t99\tPASS\tDP=100\tGT:DP\t0/1:50\t0/1:50",
];

const EXPRS: [&str; 10] = [
    "QUAL>10",
    "INFO/DP>5 && AF[0]<0.4",
    "DB",
    "TYPE=\"indel\"",
    "FILTER=\"PASS\"",
    "GT[0]=\"het\"",
    "FMT/DP[1]>10",
    "N_ALT>1",
    "CHROM=\"2\" || POS<150",
    "MAX(FORMAT/DP)>=50",
];

/// `eval_ref` on the borrowed record gives the same answer as `eval` on
/// the owned one.
#[test]
fn borrowed_and_owned_records_evaluate_alike() {
    let hdr = headers();
    for expr in EXPRS {
        let engine =
            FilterEngine::new(&hdr, Some(expr), false).unwrap_or_else(|e| panic!("{expr}: {e}"));
        for line in LINES {
            let owned = parse_vcf_line(line).unwrap();
            let borrowed = VcfRecordRef::parse(line).unwrap();
            let a = engine.eval(&owned).unwrap();
            let b = engine.eval_ref(&borrowed).unwrap();
            assert_eq!(a.pass_site, b.pass_site, "{expr} on {line}");
            assert_eq!(a.pass_samples, b.pass_samples, "{expr} on {line}");
        }
    }
}

#[test]
fn exclude_inverts_the_site_result() {
    let hdr = headers();
    let inc = FilterEngine::new(&hdr, Some("QUAL>10"), false).unwrap();
    let exc = FilterEngine::new(&hdr, Some("QUAL>10"), true).unwrap();
    let r = VcfRecordRef::parse(LINES[0]).unwrap();
    assert!(inc.eval_ref(&r).unwrap().pass_site);
    assert!(!exc.eval_ref(&r).unwrap().pass_site);
}
