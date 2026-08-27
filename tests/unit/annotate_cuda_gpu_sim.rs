    use super::*;
    use kira_kv_engine::{IndexBuilder, IndexConfig};

    /// Cross-validation: simulator's `Found(idx)` must equal `Index::lookup_u64`.
    #[test]
    fn simulator_matches_cpu_index_on_built_keys() {
        let mut keys: Vec<u64> = (0..1024u64).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        keys.sort_unstable();
        keys.dedup();
        let key_bytes: Vec<[u8; 8]> = keys.iter().map(|k| k.to_le_bytes()).collect();

        let cfg = IndexConfig::default();
        let index = IndexBuilder::new()
            .with_config(cfg)
            .build_index(key_bytes.clone())
            .expect("build small index");

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

    /// Negative-lookup correctness: a key NOT in the build set must be
    /// rejected by either bloom or fingerprint.
    #[test]
    fn simulator_rejects_foreign_keys() {
        let keys: Vec<u64> = (0..512u64).map(|i| i * 2).collect();
        let key_bytes: Vec<[u8; 8]> = keys.iter().map(|k| k.to_le_bytes()).collect();
        let index = IndexBuilder::new()
            .with_config(IndexConfig::default())
            .build_index(key_bytes)
            .expect("build");
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
