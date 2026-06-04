    use super::*;
    use crate::bcf::header::parse_header_to_dict;

    fn dict_with_dp() -> BcfHeaderDict {
        let h = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "##contig=<ID=1>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2".to_string(),
        ];
        parse_header_to_dict(&h)
    }

    #[test]
    fn record_roundtrip_basic() {
        let dict = dict_with_dp();
        let line = "1\t100\t.\tA\tT\t60\tPASS\tDP=10\tGT\t0/1\t1|1";
        let mut buf = Vec::new();
        encode_record(&mut buf, line, &dict).unwrap();
        let mut c = Cursor::new(buf);
        let out = decode_record_to_vcf(&mut c, &dict).unwrap().unwrap();
        eprintln!("OUT: {out}");
        assert!(out.starts_with("1\t100\t.\tA\tT"), "got {out:?}");
        assert!(out.contains("DP=10"), "got {out:?}");
        assert!(out.contains("GT"), "got {out:?}");
    }
