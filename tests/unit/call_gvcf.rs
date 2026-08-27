    use super::*;

    #[test]
    fn collapses_adjacent_ref_sites_same_bin() {
        let mut b = GvcfBlocker::new(vec![0, 10, 20]);
        let mut out: Vec<u8> = Vec::new();
        b.add_ref_site("1", 100, "A", 15, 50.0, &mut out).unwrap();
        b.add_ref_site("1", 101, "C", 16, 60.0, &mut out).unwrap();
        b.add_ref_site("1", 102, "G", 17, 70.0, &mut out).unwrap();
        b.flush(&mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert!(s.contains("END=102"));
        assert!(s.contains("MIN_DP=15"));
        assert!(s.contains("<NON_REF>"));
    }

    #[test]
    fn flushes_on_bin_change() {
        let mut b = GvcfBlocker::new(vec![0, 10, 20]);
        let mut out: Vec<u8> = Vec::new();
        b.add_ref_site("1", 100, "A", 15, 50.0, &mut out).unwrap();
        b.add_ref_site("1", 101, "C", 25, 60.0, &mut out).unwrap();
        b.flush(&mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        let blocks: Vec<&str> = s.lines().collect();
        assert_eq!(blocks.len(), 2, "expected 2 blocks, got {:?}", blocks);
    }

    #[test]
    fn flushes_on_chrom_change() {
        let mut b = GvcfBlocker::new(vec![0]);
        let mut out: Vec<u8> = Vec::new();
        b.add_ref_site("1", 100, "A", 10, 50.0, &mut out).unwrap();
        b.add_ref_site("2", 100, "A", 10, 50.0, &mut out).unwrap();
        b.flush(&mut out).unwrap();
        let s = String::from_utf8(out).unwrap();
        assert_eq!(s.lines().count(), 2);
    }
