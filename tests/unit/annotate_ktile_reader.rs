    use super::*;
    use crate::annotate::ktile::writer::write_ktile_from_vcf;
    use std::io::Write;

    fn write_tmp_vcf(name: &str, body: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "kira_ktile_reader_strategy_{}_{}.vcf",
            name,
            std::process::id()
        ));
        let mut f = File::create(&path).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        path
    }

    fn tiny_vcf() -> &'static str {
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t100\trs1\tA\tT\t.\t.\t.\n\
         1\t200\trs2\tC\tG\t.\t.\t.\n\
         2\t50\t.\tA\tAT\t.\t.\t.\n"
    }

    /// Combined test — env vars are process-wide and `cargo test` runs
    /// tests in parallel by default, so we sequence all three strategies
    /// (Whole, Sliding, Compressed) in a single test to avoid the race.
    #[test]
    fn all_three_strategies_round_trip() {
        use crate::annotate::ktile::writer::{
            KtileWriteOptions, write_ktile_from_vcf_with,
        };

        let vcf = write_tmp_vcf("strategies", tiny_vcf());
        let ktile_uncompressed = vcf.with_extension("ktile.uncompressed");
        let ktile_compressed = vcf.with_extension("ktile.compressed");

        // Build both flavours from the same VCF.
        write_ktile_from_vcf_with(
            &vcf,
            &ktile_uncompressed,
            KtileWriteOptions {
                compressed: false,
                lines_per_chunk: 0,
                deflate_level: 0,
            },
        )
        .unwrap();
        write_ktile_from_vcf_with(
            &vcf,
            &ktile_compressed,
            KtileWriteOptions {
                compressed: true,
                lines_per_chunk: 2, // forces multiple chunks for 3 records
                deflate_level: 1,
            },
        )
        .unwrap();

        // ---- Phase 1: uncompressed, default threshold → Whole ----
        let reader = KtileReader::open(&ktile_uncompressed).unwrap();
        assert_eq!(reader.n_records(), 3);
        assert_eq!(reader.line_owned(1), "1\t200\trs2\tC\tG\t.\t.\t.");
        drop(reader);

        // ---- Phase 2: uncompressed, force Sliding via env var ----
        unsafe {
            std::env::set_var("KIRA_BT_KTILE_MMAP_MAX_MB", "0");
            std::env::set_var("KIRA_BT_KTILE_CHUNK_MB", "0");
        }
        let reader = KtileReader::open(&ktile_uncompressed).unwrap();
        let sliding_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            assert!(reader.is_sliding(), "sliding mode must engage at 0 MB threshold");
            assert_eq!(reader.n_records(), 3);
            assert_eq!(reader.line_owned(0), "1\t100\trs1\tA\tT\t.\t.\t.");
            assert_eq!(reader.line_owned(1), "1\t200\trs2\tC\tG\t.\t.\t.");
            assert_eq!(reader.line_owned(2), "2\t50\t.\tA\tAT\t.\t.\t.");
            assert_eq!(reader.line_owned(0), "1\t100\trs1\tA\tT\t.\t.\t.");
        }));
        drop(reader);
        unsafe {
            std::env::remove_var("KIRA_BT_KTILE_MMAP_MAX_MB");
            std::env::remove_var("KIRA_BT_KTILE_CHUNK_MB");
        }

        // ---- Phase 3: compressed mode (LRU-cached decompression) ----
        let reader = KtileReader::open(&ktile_compressed).unwrap();
        let compressed_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            assert!(reader.is_compressed(), "compressed pool must auto-detect");
            assert_eq!(reader.n_records(), 3);
            // Sequential reads — first 2 in chunk 0, third in chunk 1.
            assert_eq!(reader.line_owned(0), "1\t100\trs1\tA\tT\t.\t.\t.");
            assert_eq!(reader.line_owned(1), "1\t200\trs2\tC\tG\t.\t.\t.");
            assert_eq!(reader.line_owned(2), "2\t50\t.\tA\tAT\t.\t.\t.");
            // Out-of-order — re-decompress chunk 0 from the LRU.
            assert_eq!(reader.line_owned(0), "1\t100\trs1\tA\tT\t.\t.\t.");
            // Phase 3 columns still work through the decompressed buffer.
            assert!(reader.has_ref_alt_columns());
            assert_eq!(reader.ref_slice(2).as_deref(), Some(b"A".as_slice()));
            assert_eq!(reader.alt_slice(2).as_deref(), Some(b"AT".as_slice()));
        }));

        let _ = std::fs::remove_file(&vcf);
        let _ = std::fs::remove_file(&ktile_uncompressed);
        let _ = std::fs::remove_file(&ktile_compressed);
        if let Err(p) = sliding_result {
            std::panic::resume_unwind(p);
        }
        if let Err(p) = compressed_result {
            std::panic::resume_unwind(p);
        }
    }
