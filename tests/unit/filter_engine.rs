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

#[test]
fn file_lists_preserve_exact_scalar_membership_after_loading() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("whitelist.txt");
    let mut listed = vec![
        "hit".to_string(),
        "hit suffix".to_string(),
        "marker_λ".to_string(),
    ];
    listed.extend((0..4096).map(|i| format!("entry_{i:04}")));
    let content = format!("\r\n  hit  \r\nhit\r\n{}\r\n\t\r\n", listed.join("\r\n"));
    std::fs::write(&path, content).unwrap();
    let hdr = headers();
    for field in ["ID", "CHROM", "FILTER", "ALT"] {
        for (op, exclude) in [("=", false), ("!=", false), ("=", true)] {
            let expression = format!("{field}{op}@{}", path.display());
            let engine = FilterEngine::new(&hdr, Some(&expression), exclude).unwrap();
            for value in [
                "hit",
                "hit suffix",
                "marker_λ",
                "entry_0000",
                "entry_4095",
                "entry_4096",
                "entry_0000x",
                "marker_Λ",
                "Hit",
                ".",
                "missing",
            ] {
                let line = match field {
                    "ID" => format!("1\t1\t{value}\tA\tC\t.\t.\t."),
                    "CHROM" => format!("{value}\t1\t.\tA\tC\t.\t.\t."),
                    "FILTER" => format!("1\t1\t.\tA\tC\t.\t{value}\t."),
                    _ => format!("1\t1\t.\tA\t{value}\t.\t.\t."),
                };
                let contains = listed.iter().any(|entry| entry == value);
                let expected = if value == "." && field == "FILTER" {
                    exclude
                } else {
                    ((op == "=") == contains) != exclude
                };
                let owned = parse_vcf_line(&line).unwrap();
                let borrowed = VcfRecordRef::parse(&line).unwrap();
                for result in [
                    engine.eval(&owned).unwrap(),
                    engine.eval_ref(&borrowed).unwrap(),
                ] {
                    assert_eq!(result.pass_site, expected, "{expression}: {value}");
                    assert_eq!(result.pass_samples, None);
                }
            }
        }
    }
    let expression = format!("ID=@{}", path.display());
    let engine = FilterEngine::new(&hdr, Some(&expression), false).unwrap();
    std::fs::write(&path, "replacement\n").unwrap();
    let record = VcfRecordRef::parse("1\t1\thit\tA\tC\t.\t.\t.").unwrap();
    assert!(engine.eval_ref(&record).unwrap().pass_site);
}

#[test]
fn file_lists_preserve_sample_masks_and_missing_values() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("sample-list.txt");
    std::fs::write(&path, " hit \nhit\nmarker_λ\n\n").unwrap();
    let hdr = vec![
        "##fileformat=VCFv4.2".into(),
        "##FORMAT=<ID=LABEL,Number=.,Type=String,Description=\"\">".into(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3".into(),
    ];
    let lines = [
        ("hit\tmiss\t.", [true, false, false], [false, true, false]),
        (
            "miss,hit\tmarker_λ\tother",
            [true, true, false],
            [true, false, true],
        ),
        (".\t.\t.", [false; 3], [false; 3]),
    ];
    for (samples, equal, unequal) in lines {
        let line = format!("1\t1\t.\tA\tC\t.\t.\t.\tLABEL\t{samples}");
        let owned = parse_vcf_line(&line).unwrap();
        let borrowed = VcfRecordRef::parse(&line).unwrap();
        for (op, expected) in [("=", equal), ("!=", unequal)] {
            for reverse in [false, true] {
                let list = format!("@{}", path.display());
                let expression = if reverse {
                    format!("{list}{op}FMT/LABEL")
                } else {
                    format!("FMT/LABEL{op}{list}")
                };
                let engine = FilterEngine::new(&hdr, Some(&expression), false).unwrap();
                for result in [
                    engine.eval(&owned).unwrap(),
                    engine.eval_ref(&borrowed).unwrap(),
                ] {
                    assert_eq!(
                        result.pass_site,
                        expected.into_iter().any(|v| v),
                        "{expression}"
                    );
                    assert_eq!(result.pass_samples, Some(expected.to_vec()), "{expression}");
                }
            }
        }
        let expression = format!("FMT/LABEL[1]=@{}", path.display());
        let engine = FilterEngine::new(&hdr, Some(&expression), false).unwrap();
        let result = engine.eval_ref(&borrowed).unwrap();
        assert_eq!(result.pass_samples, Some(vec![false, equal[1], false]));
        let engine = FilterEngine::new(&hdr, Some(&expression), true).unwrap();
        let result = engine.eval_ref(&borrowed).unwrap();
        assert_eq!(result.pass_site, !equal[1]);
        assert_eq!(result.pass_samples, Some(vec![true, !equal[1], true]));
    }
}
