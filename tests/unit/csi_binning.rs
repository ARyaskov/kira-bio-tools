    use super::*;

    const B: u64 = 1 << 16; // one BGZF block in virtual-position units

    fn sample() -> BinIndex {
        let mut bins = BTreeMap::new();
        bins.insert(4681, BinEntry { loffset: B, chunks: vec![(B, 2 * B), (5 * B, 6 * B)] });
        bins.insert(4682, BinEntry { loffset: 6 * B, chunks: vec![(6 * B, 7 * B)] });
        bins.insert(585, BinEntry { loffset: B, chunks: vec![(10, 3 * B)] });
        let r = RefIndex {
            bins,
            linear: vec![B, 6 * B],
            meta: Some(RefMeta { beg_off: 10, end_off: 7 * B, n_mapped: 7, n_unmapped: 0 }),
        };
        BinIndex {
            kind: IndexKind::Csi,
            min_shift: 14,
            depth: 5,
            header: Some(TabixHeader::vcf(vec!["chr1".into(), "chr2".into()])),
            refs: vec![r, RefIndex::default()],
            n_no_coor: Some(0),
        }
    }

    #[test]
    fn csi_roundtrip() {
        let idx = sample();
        let bytes = idx.to_bytes();
        let back = BinIndex::parse(&bytes).unwrap();
        assert_eq!(back.kind, IndexKind::Csi);
        assert_eq!(back.min_shift, 14);
        assert_eq!(back.depth, 5);
        assert_eq!(back.names(), &["chr1", "chr2"]);
        assert_eq!(back.refs.len(), 2);
        assert_eq!(back.refs[0].bins[&4681].chunks, vec![(B, 2 * B), (5 * B, 6 * B)]);
        assert_eq!(back.refs[0].bins[&4681].loffset, B);
        assert_eq!(back.n_records(0), Some(7));
        assert_eq!(back.n_records(1), None);
        assert_eq!(back.n_no_coor, Some(0));
    }

    #[test]
    fn tbi_roundtrip() {
        let mut idx = sample();
        idx.kind = IndexKind::Tbi;
        let bytes = idx.to_bytes();
        assert_eq!(&bytes[..4], b"TBI\x01");
        let back = BinIndex::parse(&bytes).unwrap();
        assert_eq!(back.kind, IndexKind::Tbi);
        assert_eq!(back.refs[0].linear, vec![B, 6 * B]);
        assert_eq!(back.refs[0].bins[&4682].chunks, vec![(6 * B, 7 * B)]);
        assert_eq!(back.n_records(0), Some(7));
        assert_eq!(back.ref_id("chr2"), Some(1));
    }

    #[test]
    fn query_uses_bins_and_min_offset() {
        let idx = sample();
        // First 16 kb window: parent-bin chunk clipped to loffset of leaf 4681,
        // merged with the leaf's first chunk; the second leaf chunk is a
        // separate block and stays separate.
        let c = idx.query(0, 0, 10);
        assert_eq!(c, vec![(B, 3 * B), (5 * B, 6 * B)]);
        // Second window: min_off = loffset(4682) drops the earlier chunks.
        let c = idx.query(0, 16384, 16390);
        assert_eq!(c, vec![(6 * B, 7 * B)]);
        assert!(idx.query(1, 0, 10).is_empty());
        assert!(idx.query(5, 0, 10).is_empty());
        assert!(idx.query(0, 10, 10).is_empty());

        // TBI uses the linear index for the minimum offset.
        let mut t = sample();
        t.kind = IndexKind::Tbi;
        assert_eq!(t.query(0, 16384, 16390), vec![(6 * B, 7 * B)]);
    }

    #[test]
    fn rejects_garbage() {
        assert!(BinIndex::parse(b"XYZ\x01").is_err());
        assert!(BinIndex::parse(b"CSI\x01\x0e\x00\x00").is_err());
    }
