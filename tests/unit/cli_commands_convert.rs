use super::*;

#[test]
fn pl_conversion_uses_reference_accumulation_precision() {
    assert_eq!(
        sample_to_gen_probs("0|0:0,12,62", Some("GT:PL"), "PL"),
        ("0.940649".into(), "0.059351".into(), "0.000001".into())
    );
}
