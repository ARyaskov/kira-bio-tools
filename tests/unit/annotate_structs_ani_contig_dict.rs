    use super::*;

    #[test]
    fn roundtrip_empty() {
        let dict = ContigDict::default();
        let mut buf = Vec::new();
        dict.write_to(&mut buf).unwrap();
        let restored = ContigDict::parse_bytes(&buf).unwrap();
        assert_eq!(restored.len(), 0);
    }

    #[test]
    fn roundtrip_named() {
        let dict = ContigDict::from_names(["chr1", "chr2", "chrX", "chrM"]).unwrap();
        let mut buf = Vec::new();
        dict.write_to(&mut buf).unwrap();
        assert_eq!(buf.len(), dict.serialized_len());
        let restored = ContigDict::parse_bytes(&buf).unwrap();
        assert_eq!(restored.len(), 4);
        assert_eq!(restored.id("chr1"), Some(0));
        assert_eq!(restored.id("chr2"), Some(1));
        assert_eq!(restored.id("chrX"), Some(2));
        assert_eq!(restored.id("chrM"), Some(3));
        assert_eq!(restored.id("chrUnknown"), None);
        assert_eq!(restored.name(0), Some("chr1"));
    }

    #[test]
    fn parses_vcf_header_lines() {
        let lines = [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=249250621>",
            "##contig=<ID=chrM,length=16571>",
            "##contig=<ID=GL000207.1,length=4262>",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"...\">",
        ];
        let dict = ContigDict::from_header_lines(lines.iter().copied());
        assert_eq!(dict.len(), 3);
        assert_eq!(dict.id("chr1"), Some(0));
        assert_eq!(dict.id("chrM"), Some(1));
        assert_eq!(dict.id("GL000207.1"), Some(2));
    }

    #[test]
    fn insert_dedupes() {
        let mut dict = ContigDict::default();
        let a1 = dict.insert("chr1");
        let a2 = dict.insert("chr2");
        let a3 = dict.insert("chr1"); // duplicate
        assert_eq!(a1, 0);
        assert_eq!(a2, 1);
        assert_eq!(a3, 0); // same id reused
        assert_eq!(dict.len(), 2);
    }

    #[test]
    fn duplicate_in_from_names_is_rejected() {
        assert!(ContigDict::from_names(["chr1", "chr1"]).is_none());
    }
