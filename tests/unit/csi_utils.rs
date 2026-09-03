    use super::*;

    #[test]
    fn bin_arithmetic_matches_htslib() {
        assert_eq!(bin_first(0), 0);
        assert_eq!(bin_first(1), 1);
        assert_eq!(bin_first(2), 9);
        assert_eq!(bin_first(5), 4681);
        assert_eq!(n_bins(5), 37449);
        assert_eq!(metadata_bin(5), 37450);
        assert_eq!(bin_parent(4681), 585);
        assert_eq!(bin_level(4681, 5), 5);
        assert_eq!(bin_level(0, 5), 0);
        assert_eq!(bin_bot(4681, 5), 0);
        assert_eq!(bin_bot(1, 5), 0);
        assert_eq!(bin_bot(2, 5), 4096);
    }

    #[test]
    fn reg2bin_examples() {
        // 16 kb leaves: [0,1) -> first leaf bin.
        assert_eq!(reg2bin(0, 1, 14, 5), 4681);
        assert_eq!(reg2bin(16384, 16385, 14, 5), 4682);
        // Spanning two leaves goes up a level.
        assert_eq!(reg2bin(16000, 17000, 14, 5), 585);
        let mut bins = Vec::new();
        reg2bins(0, 1, 14, 5, &mut bins);
        assert_eq!(bins, vec![0, 1, 9, 73, 585, 4681]);
        reg2bins(0, 0, 14, 5, &mut bins);
        assert!(bins.is_empty());
    }

    #[test]
    fn depth_from_contig_length() {
        assert_eq!(depth_for(14, 249_250_621), 5);
        assert_eq!(depth_for(14, 2_147_483_647), 6);
        assert_eq!(depth_for(14, 10_000), 0);
        assert_eq!(tabix_depth_for(14), 6);
    }
