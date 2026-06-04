    use super::process_vcf_line_multiallelic_simd;
    use crate::annotate::builder_v2::StringPool;
    use crate::annotate::builder_v2::entry_processing::{EntryEntry, make_position_key};
    use crate::annotate::structs::ani::{ANI_STR_NONE, ContigDict};
    use crate::annotate::structs::bundle::FieldNumber;
    use crate::util::read_cstring;
    use fxhash::FxHashMap;
    use std::collections::HashMap;
    use std::sync::atomic::AtomicUsize;

    #[test]
    fn test_samples_are_not_truncated() {
        let line = b"1\t3000002\tid\tC\tT\t99\tq99\tFLAG;IINT=88,99;IFLT=8.8,9.9;ISTR=888,999\tGT:FINT:FFLT:FSTR\t1|1:88,99:8.8,9.9:888,999\t0|1:77:7.7:77";
        let mut contigs = ContigDict::default();
        let mut entries_map: FxHashMap<u64, EntryEntry> = FxHashMap::default();
        let mut pool = StringPool::new();
        let mut insertion_order = 0usize;
        let duplicates_skipped = AtomicUsize::new(0);
        let collisions_detected = AtomicUsize::new(0);
        let multiallelic_count = AtomicUsize::new(0);
        let field_meta: HashMap<String, FieldNumber> = HashMap::new();
        let format_meta: HashMap<String, FieldNumber> = HashMap::new();

        let processed = process_vcf_line_multiallelic_simd(
            line,
            &mut contigs,
            &mut entries_map,
            &mut pool,
            &mut insertion_order,
            &duplicates_skipped,
            &collisions_detected,
            &multiallelic_count,
            &field_meta,
            &format_meta,
            2,
            false,
        )
        .unwrap();

        assert_eq!(processed, 1);
        let entry = entries_map.values().next().unwrap().entry;
        assert_ne!(entry.samples_ofs, ANI_STR_NONE);
        let pool_bytes = pool.materialize().unwrap();
        let samples = read_cstring(&pool_bytes, entry.samples_ofs as usize);
        assert_eq!(samples, "1|1:88,99:8.8,9.9:888,999\t0|1:77:7.7:77");
    }

    #[test]
    fn test_multiallelic_format_a_r_g_are_split_per_alt() {
        let line = b"1\t10\t.\tC\tA,T\t.\t.\t.\tGT:AD:AF:PL\t0/2:10,1,9:0.1,0.9:0,1,2,3,4,5";
        let mut contigs = ContigDict::default();
        let mut entries_map: FxHashMap<u64, EntryEntry> = FxHashMap::default();
        let mut pool = StringPool::new();
        let mut insertion_order = 0usize;
        let duplicates_skipped = AtomicUsize::new(0);
        let collisions_detected = AtomicUsize::new(0);
        let multiallelic_count = AtomicUsize::new(0);
        let field_meta: HashMap<String, FieldNumber> = HashMap::new();
        let mut format_meta: HashMap<String, FieldNumber> = HashMap::new();
        format_meta.insert("AD".to_string(), FieldNumber::R);
        format_meta.insert("AF".to_string(), FieldNumber::A);
        format_meta.insert("PL".to_string(), FieldNumber::G);

        process_vcf_line_multiallelic_simd(
            line,
            &mut contigs,
            &mut entries_map,
            &mut pool,
            &mut insertion_order,
            &duplicates_skipped,
            &collisions_detected,
            &multiallelic_count,
            &field_meta,
            &format_meta,
            1,
            false,
        )
        .unwrap();

        let key = make_position_key(0, 10, "C", "T");
        let entry = entries_map.get(&key).unwrap().entry;
        let pool_bytes = pool.materialize().unwrap();
        let samples = read_cstring(&pool_bytes, entry.samples_ofs as usize);
        assert_eq!(samples, "0/2:10,9:0.9:0,3,5");
    }
