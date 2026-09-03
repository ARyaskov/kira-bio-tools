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
        assert_eq!(d.filter_idx["PASS"], 0);
        assert_eq!(d.info_field(d.info_idx["AC"]).unwrap().id, "AC");
        assert_eq!(d.contig_name(0), Some("1"));
    }

    #[test]
    fn shared_id_space_like_htslib() {
        let h = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "##FILTER=<ID=q10,Description=\"low q\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let d = parse_header_to_dict(&h);
        // Same tag name in INFO and FORMAT shares one IDX; PASS is 0.
        assert_eq!(d.info_idx["DP"], d.format_idx["DP"]);
        assert_eq!(d.info_idx["DP"], 1);
        assert_eq!(d.format_idx["GT"], 2);
        assert_eq!(d.filter_idx["q10"], 3);
        let s = serialize_header(&d);
        assert!(s.contains("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\",IDX=1>"));
        assert!(s.contains("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\",IDX=1>"));
        assert!(s.contains("##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>"));
    }

    #[test]
    fn explicit_idx_is_honoured() {
        let h = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>".to_string(),
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"AC\",IDX=5>".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"D\",IDX=2>".to_string(),
            "##contig=<ID=chr2,length=10,IDX=1>".to_string(),
            "##contig=<ID=chr1,length=10,IDX=0>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let d = parse_header_to_dict(&h);
        assert_eq!(d.info_idx["AC"], 5);
        assert_eq!(d.info_idx["DP"], 2);
        assert_eq!(d.contig_name(0), Some("chr1"));
        assert_eq!(d.contig_name(1), Some("chr2"));
        assert_eq!(d.info_field(5).unwrap().id, "AC");
        assert!(d.info_field(3).is_none());
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
        let lines: Vec<String> = s.lines().map(|l| l.to_string()).collect();
        let d2 = parse_header_to_dict(&lines);
        assert_eq!(d2.info_idx["DP"], d.info_idx["DP"]);
    }
