    use super::*;

    #[test]
    fn region_parsing() {
        let r = Region::parse("chr1:1,000-2,000").unwrap();
        assert_eq!((r.chr.as_str(), r.start, r.end), ("chr1", Some(1000), Some(2000)));
        let r = Region::parse("chr1").unwrap();
        assert_eq!(r.bounds(), (1, u32::MAX));
        let r = Region::parse("chr1:100").unwrap();
        assert_eq!(r.bounds(), (100, 100));
        let r = Region::parse_with("chr1:100", false).unwrap();
        assert_eq!(r.bounds(), (100, u32::MAX));
        let r = Region::parse("chr1:100-").unwrap();
        assert_eq!(r.bounds(), (100, u32::MAX));
        let r = Region::parse("HLA-A*01:01:01:01:5-10").unwrap();
        assert_eq!(r.chr, "HLA-A*01:01:01:01");
        assert_eq!(r.bounds(), (5, 10));
        assert_eq!(parse_coordinate("1.5k"), None);
        assert_eq!(parse_coordinate("2k"), Some(2000));
        assert_eq!(parse_coordinate("3M"), Some(3_000_000));
    }

    #[test]
    fn genomic_key_holds_large_contig_ids() {
        let k = GenomicKey::new(70_000, 123);
        assert_eq!(k.chr(), 70_000);
        assert_eq!(k.position(), 123);
    }
