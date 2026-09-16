    use super::*;

    /// The vendored hasher matches the `fxhash` crate bit for bit; the
    /// values also pin the algorithm for `.ani` keys written so far.
    #[test]
    fn matches_the_fxhash_crate() {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for len in 0..40usize {
            for _ in 0..8 {
                let bytes: Vec<u8> = (0..len).map(|_| next() as u8).collect();
                assert_eq!(hash64(&bytes[..]), fxhash::hash64(&bytes[..]), "{bytes:?}");
                assert_eq!(hash64(&bytes), fxhash::hash64(&bytes), "vec {bytes:?}");
            }
        }
        assert_eq!(hash64("ACGT"), fxhash::hash64("ACGT"));
        assert_eq!(hash64(&123_456_789u64), fxhash::hash64(&123_456_789u64));
        assert_eq!(hash64(&(1u32, "x")), fxhash::hash64(&(1u32, "x")));
        let mut m: FxHashMap<&str, u32> = FxHashMap::default();
        m.insert("a", 1);
        assert_eq!(m.get("a"), Some(&1));
    }
