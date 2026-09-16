use super::*;

fn headers() -> Vec<String> {
    vec![
        "##fileformat=VCFv4.2".into(),
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">".into(),
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"\">".into(),
        "##INFO=<ID=XRF,Number=R,Type=Float,Description=\"\">".into(),
        "##INFO=<ID=TXT,Number=1,Type=String,Description=\"\">".into(),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">".into(),
        "##FORMAT=<ID=FGF,Number=G,Type=Float,Description=\"\">".into(),
        "##FORMAT=<ID=FGS,Number=G,Type=String,Description=\"\">".into(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB".into(),
    ]
}

#[test]
fn scalar_values_render_like_kputd() {
    for (input, want) in [
        ("1e-1", "0.1"),
        ("1e6", "1e+06"),
        ("50.0", "50"),
        ("0.500", "0.5"),
        ("3.1415927", "3.14159"),
        ("+1.5", "1.5"),
        ("007", "7"),
        ("NaN", "nan"),
        ("Infinity", "inf"),
        ("-0.0", "-0"),
        ("0.00001", "1e-05"),
        ("1234567", "1.23457e+06"),
    ] {
        assert_eq!(canonical_float(input), want, "{input}");
    }
    for keep in [
        ".", "0", "-0", "0.1", "0.0001", "999999", "5e-05", "1e+06", "1.5e+100", "abc", "", "inf",
        "-inf", "nan", "0.123457",
    ] {
        assert!(
            matches!(canonical_float(keep), Cow::Borrowed(_)),
            "{keep} should be left alone"
        );
    }
}

#[test]
fn vectors_keep_missing_entries() {
    assert_eq!(canonical_floats("1e-1,.,2e-2"), "0.1,.,0.02");
    assert!(matches!(canonical_floats("0.1,.,0.02"), Cow::Borrowed(_)));
    assert_eq!(canonical_floats("0.1,0.20"), "0.1,0.2");
}

#[test]
fn line_rewrites_only_float_typed_fields() {
    let mut c = RecordCanon::from_headers(&headers());
    let line = "1\t100\t.\tA\tC,G\t238.0\tPASS\tDP=1e1;XRF=1e6,2e6,3e6;TXT=1e6;AF=1e-1,5e-2\tGT:FGS:FGF\t0/1:1e-1,2,3:1e-1,2e-2,3e-3\t.";
    let want = "1\t100\t.\tA\tC,G\t238\tPASS\tDP=1e1;XRF=1e+06,2e+06,3e+06;TXT=1e6;AF=0.1,0.05\tGT:FGS:FGF\t0/1:1e-1,2,3:0.1,0.02,0.003\t.:.:.";
    assert_eq!(c.canon_line(line).as_deref(), Some(want));
    assert_eq!(c.canon_line(want), None);
    assert_eq!(c.canon_line("#CHROM\tPOS"), None);
    assert_eq!(c.canon_line("1\t100\t.\tA\tC\t.\t.\t."), None);
    assert_eq!(
        c.canon_line("1\t100\t.\tA\tC\t1e2\t.\t.").as_deref(),
        Some("1\t100\t.\tA\tC\t100\t.\t.")
    );
    // A FORMAT column without Float keys leaves values alone; a short
    // sample column is padded with missing fields as htslib prints it.
    assert_eq!(
        c.canon_line("1\t100\t.\tA\tC\t.\t.\t.\tGT:FGS\t0/1:1e-1\t0/1")
            .as_deref(),
        Some("1\t100\t.\tA\tC\t.\t.\t.\tGT:FGS\t0/1:1e-1\t0/1:.")
    );
    assert_eq!(
        c.canon_line("1\t100\t.\tA\tC\t.\t.\t.\tDP:AD:FR\t.\t2:2,11:2,11\t3:3")
            .as_deref(),
        Some("1\t100\t.\tA\tC\t.\t.\t.\tDP:AD:FR\t.:.:.\t2:2,11:2,11\t3:3:.")
    );
    assert_eq!(c.canon_line("1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1\t."), None);
    assert_eq!(
        c.canon_line("1\t100\t.\tA\tC\t.\t.\t.\tFGF:GT\t1e-1:0/1\t.")
            .as_deref(),
        Some("1\t100\t.\tA\tC\t.\t.\t.\tFGF:GT\t0.1:0/1\t.:.")
    );
    assert_eq!(c.canon_info("AF=1e-1;DP=1e1"), "AF=0.1;DP=1e1");
    assert!(c.info_is_float("AF") && !c.info_is_float("DP") && c.format_is_float("FGF"));
}
