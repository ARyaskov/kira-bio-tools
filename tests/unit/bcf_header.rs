    use super::*;

    #[test]
    fn parse_basic_header() {
        let h = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"AC\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "##FILTER=<ID=q10,Description=\"low q\">".to_string(),
            "##contig=<ID=1,length=249250621>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2".to_string(),
        ];
        let d = parse_header_to_dict(&h);
        assert_eq!(d.info.len(), 2);
        assert_eq!(d.format.len(), 1);
        assert_eq!(d.filter.len(), 2); // PASS prepended + q10
        assert!(d.filter_idx.contains_key("PASS"));
        assert_eq!(d.samples, vec!["NA1", "NA2"]);
        assert_eq!(d.contigs, vec!["1"]);
    }

    #[test]
    fn serialize_round_trip() {
        let h = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"D\">".to_string(),
            "##contig=<ID=1>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let d = parse_header_to_dict(&h);
        let s = serialize_header(&d);
        assert!(s.contains("##INFO=<ID=DP"));
        assert!(s.contains("IDX="));
        assert!(s.contains("#CHROM"));
    }
