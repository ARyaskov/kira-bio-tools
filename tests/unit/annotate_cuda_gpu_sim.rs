    use super::*;
    use kira_kv_engine::{IndexBuilder, IndexConfig};

    fn build(keys: &[u64]) -> kira_kv_engine::Index {
        let key_bytes: Vec<[u8; 8]> = keys.iter().map(|k| k.to_le_bytes()).collect();
        IndexBuilder::new()
            .with_config(IndexConfig::default())
            .build_index_ref(&key_bytes)
            .expect("build index")
    }

    /// Cross-validation: simulator's `Found(idx)` must equal `Index::lookup_u64`.
    #[test]
    fn simulator_matches_cpu_index_on_built_keys() {
        let mut keys: Vec<u64> = (0..1024u64).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        keys.sort_unstable();
        keys.dedup();
        let index = build(&keys);

        let export = match index.gpu_export() {
            Some(e) => e,
            None => return,
        };

        for &key in &keys {
            let cpu = index.lookup_u64(key).expect("cpu lookup");
            let sim = lookup_u64(&export, key);
            assert_eq!(
                sim,
                GpuLookup::Found(cpu as u32),
                "simulator/CPU divergence for key 0x{key:016x}: cpu={cpu}, sim={sim:?}"
            );
        }
    }

    /// kira_kv_engine 0.6.3 splits above ~32768 keys, so a real-sized index takes
    /// the partitioned path: part selector, per-part salt, bucket/slot offsets.
    #[test]
    fn simulator_matches_cpu_index_across_partitions() {
        let mut keys: Vec<u64> = (0..200_000u64).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        keys.sort_unstable();
        keys.dedup();
        let index = build(&keys);

        let export = match index.gpu_export() {
            Some(e) => e,
            None => return,
        };
        assert!(
            export.parts.len() > 1,
            "expected a multi-part index, got {} part(s)",
            export.parts.len()
        );

        for &key in &keys {
            let cpu = index.lookup_u64(key).expect("cpu lookup");
            let sim = lookup_u64(&export, key);
            assert_eq!(
                sim,
                GpuLookup::Found(cpu as u32),
                "simulator/CPU divergence for key 0x{key:016x}: cpu={cpu}, sim={sim:?}"
            );
        }
    }

    /// Negative-lookup correctness: a key NOT in the build set must be
    /// rejected by either bloom or fingerprint.
    #[test]
    fn simulator_rejects_foreign_keys() {
        let keys: Vec<u64> = (0..512u64).map(|i| i * 2).collect();
        let index = build(&keys);
        let Some(export) = index.gpu_export() else { return; };

        let mut found = 0;
        for k in (1..512u64).step_by(2).take(256) {
            if matches!(lookup_u64(&export, k), GpuLookup::Found(_)) {
                found += 1;
            }
        }
        assert!(
            found <= 1,
            "expected ≈0 foreign-key false positives, got {found}"
        );
    }

    /// The bloom lane shift comes from the export: 26 for filters built here.
    /// Hardcoding the legacy 27 would drop real keys at the prefilter.
    #[test]
    fn bloom_uses_the_exported_bit_shift() {
        let keys: Vec<u64> = (0..4096u64).map(|i| i.wrapping_mul(0xD6E8_FEB8_6659_FD93)).collect();
        let index = build(&keys);
        let Some(export) = index.gpu_export() else { return; };
        let Some(bloom) = export.bloom.as_ref() else { return };
        assert_eq!(bloom.bit_shift, 26, "0.6.3 builds a 6-bit lane index");
        for &key in &keys {
            assert!(
                bloom_contains(bloom, canonical_u64(key, export.prehash_seed)),
                "bloom rejected a key that was built into it: 0x{key:016x}"
            );
        }
    }
