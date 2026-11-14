use super::*;

#[test]
fn default_layout_is_compatible_and_native_layout_is_explicit() {
    assert!(
        MpileupArgs::try_parse_from(["mpileup", "in.bam"])
            .unwrap()
            .bcftools_layout()
    );
    assert!(
        !MpileupArgs::try_parse_from(["mpileup", "--native-layout", "in.bam"])
            .unwrap()
            .bcftools_layout()
    );
    assert!(
        MpileupArgs::try_parse_from(["mpileup", "--native-layout", "--bcftools-layout", "in.bam"])
            .unwrap()
            .bcftools_layout()
    );
    assert!(
        !MpileupArgs::try_parse_from(["mpileup", "--variants-only", "in.bam"])
            .unwrap()
            .bcftools_layout()
    );
}
