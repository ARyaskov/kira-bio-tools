    use super::*;

    #[test]
    fn contig_dict_from_header() {
        let lines = [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=248956422>",
            "##contig=<ID=chr2>",
            "##contig=<ID=chr1_KI270706v1_random,length=175055,assembly=\"a,b\">",
        ];
        let d = ContigDict::from_header_lines(lines.iter().copied());
        assert_eq!(d.len(), 3);
        assert_eq!(d.id("chr1"), Some(0));
        assert_eq!(d.id("chr2"), Some(1));
        assert_eq!(d.id("chr1_KI270706v1_random"), Some(2));
        assert_eq!(d.length(0), Some(248956422));
        assert_eq!(d.length(1), None);
        assert_eq!(d.length(2), Some(175055));
        assert_eq!(d.name(2), Some("chr1_KI270706v1_random"));
        assert_eq!(d.id("chrX"), None);
    }

    #[test]
    fn contig_dict_insert_appends() {
        let mut d = ContigDict::new();
        assert_eq!(d.insert("scaffold_1"), 0);
        assert_eq!(d.insert("scaffold_2"), 1);
        assert_eq!(d.insert("scaffold_1"), 0);
        assert_eq!(d.len(), 2);
    }

    #[test]
    fn header_info_numbers() {
        let lines = vec![
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count, per ALT\">".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">".to_string(),
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">".to_string(),
            "##FILTER=<ID=q10,Description=\"low\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2".to_string(),
        ];
        let h = HeaderInfo::parse(&lines);
        assert_eq!(h.info_number("AC"), FieldNumber::A);
        assert_eq!(h.info_number("DP"), FieldNumber::Fixed(1));
        assert_eq!(h.info_number("NOPE"), FieldNumber::Dot);
        assert_eq!(h.format_number("PL"), FieldNumber::G);
        assert_eq!(h.format_number("AD"), FieldNumber::R);
        assert_eq!(h.info.get("AC").unwrap().description, "Allele count, per ALT");
        assert_eq!(h.filters, vec!["q10"]);
        assert_eq!(h.samples, vec!["S1", "S2"]);
    }
