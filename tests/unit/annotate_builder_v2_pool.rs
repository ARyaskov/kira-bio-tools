    use super::*;

    #[test]
    fn interns_repeated_medium_strings_to_same_offset() {
        let mut pool = StringPool::new();
        let rs1 = pool.append_cstr("rs148327885");
        let snv1 = pool.append_cstr("single_nucleotide_variant");
        let rs2 = pool.append_cstr("rs148327885");
        let snv2 = pool.append_cstr("single_nucleotide_variant");
        assert_eq!(rs1, rs2);
        assert_eq!(snv1, snv2);
        assert_eq!(pool.intern_hits(), 2);
        assert_eq!(pool.intern_bytes_saved(), 12 + 26);
    }

    #[test]
    fn tiny_strings_below_intern_min_are_not_interned() {
        let mut pool = StringPool::new();
        let a = pool.append_cstr(".");
        let b = pool.append_cstr(".");
        let c = pool.append_cstr("PASS");
        let d = pool.append_cstr("PASS");
        assert_ne!(a, b);
        assert_ne!(c, d);
        assert_eq!(pool.intern_hits(), 0);
    }

    #[test]
    fn long_strings_bypass_interning_and_always_advance_offset() {
        let mut pool = StringPool::new();
        let long_a = "X".repeat(INTERN_MAX_LEN + 1);
        let long_b = "X".repeat(INTERN_MAX_LEN + 1);
        let ofs_a = pool.append_cstr(&long_a);
        let ofs_b = pool.append_cstr(&long_b);
        assert_ne!(ofs_a, ofs_b, "long strings must not intern");
        assert_eq!(pool.intern_hits(), 0);
    }

    #[test]
    fn interned_strings_decode_to_original_content_on_lookup() {
        let mut pool = StringPool::new();
        let ofs_chr = pool.append_cstr("chr1_KZ208906v1");
        let ofs_rs = pool.append_cstr("rs148327885");
        let _ofs_chr2 = pool.append_cstr("chr1_KZ208906v1");
        let raw = pool.materialize().unwrap();
        let n_chr = "chr1_KZ208906v1".len();
        let n_rs = "rs148327885".len();
        assert_eq!(&raw[ofs_chr..ofs_chr + n_chr], b"chr1_KZ208906v1");
        assert_eq!(raw[ofs_chr + n_chr], 0);
        assert_eq!(&raw[ofs_rs..ofs_rs + n_rs], b"rs148327885");
        assert_eq!(raw[ofs_rs + n_rs], 0);
    }

    #[test]
    fn cleanup_releases_intern_map() {
        let mut pool = StringPool::new();
        for i in 0..1000 {
            pool.append_cstr(&format!("variant_id_{i:08}"));
        }
        assert!(pool.intern.as_ref().is_some_and(|m| m.len() >= 1000));
        pool.cleanup();
        assert!(pool.intern.is_none());
    }
