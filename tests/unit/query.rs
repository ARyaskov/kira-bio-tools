use super::*;

#[test]
fn repeated_header_flag_removes_column_numbers() {
    let samples = vec!["A".to_string(), "B".to_string()];
    let context = QueryFormatContext::new(
        &[],
        &samples,
        vec![1, 0],
        "%CHROM\t[%SAMPLE\t]\n".to_string(),
    )
    .unwrap();
    assert_eq!(context.header_line(false), "#[1]CHROM\t[2]SAMPLE\t\n");
    assert_eq!(context.header_line(true), "#CHROM\tB:SAMPLE\tA:SAMPLE\t\n");
    let args = ["-H", "-H", "-f", "%POS", "in.vcf"].map(String::from);
    assert_eq!(parse_query_args(&args).unwrap().header_mode, 2);
}
