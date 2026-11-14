use super::*;
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
    use crate::annotate::ktile::writer::{KtileWriteOptions, write_ktile_from_vcf_with};

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
        assert!(
            reader.is_sliding(),
            "sliding mode must engage at 0 MB threshold"
        );
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

fn legacy_fixture() -> (Vec<u8>, KtileHeader, Vec<&'static str>) {
    use crate::annotate::ktile::format::{KTILE_LEGACY_VERSION, KTILE_MAGIC};
    use bytemuck::Zeroable;
    let lines = vec![
        "1\t7\tfirst\tA\tG\t.\t.\t.",
        "2\t19\tsecond\tAT\tA,C\t.\t.\t.",
    ];
    let headers = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
    let mut bytes = vec![0; KtileHeader::SIZE];
    let mut header = KtileHeader::zeroed();
    header.magic = KTILE_MAGIC;
    header.version = KTILE_LEGACY_VERSION;
    header.n_records = lines.len() as u64;
    header.headers_off = bytes.len() as u64;
    header.headers_len = headers.len() as u64;
    bytes.extend_from_slice(headers);
    bytes.resize(bytes.len().next_multiple_of(8), 0);
    header.line_pool_off = bytes.len() as u64;
    for line in &lines {
        bytes.extend_from_slice(line.as_bytes());
    }
    header.line_pool_len = lines.iter().map(|line| line.len() as u64).sum();
    bytes.resize(bytes.len().next_multiple_of(8), 0);
    header.line_offsets_off = bytes.len() as u64;
    for offset in [0, lines[0].len() as u64, header.line_pool_len] {
        bytes.extend_from_slice(&offset.to_le_bytes());
    }
    header.chr_ids_off = bytes.len() as u64;
    for chrom in [1u32, 2] {
        bytes.extend_from_slice(&chrom.to_le_bytes());
    }
    header.positions_off = bytes.len() as u64;
    for position in [7u32, 19] {
        bytes.extend_from_slice(&position.to_le_bytes());
    }
    bytes[..KtileHeader::SIZE].copy_from_slice(bytemuck::bytes_of(&header));
    (bytes, header, lines)
}

fn assert_invalid_ktile(path: &Path, bytes: &[u8]) {
    std::fs::write(path, bytes).unwrap();
    let result = std::panic::catch_unwind(|| KtileReader::open(path));
    assert!(
        matches!(result, Ok(Err(_))),
        "malformed file was accepted or panicked"
    );
}

#[test]
fn legacy_fixture_reads_without_using_the_current_writer() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("legacy.ktile");
    let (bytes, _, lines) = legacy_fixture();
    std::fs::write(&path, bytes).unwrap();
    let reader = KtileReader::open(&path).unwrap();
    assert_eq!(reader.header().version, 1);
    assert_eq!(reader.n_records(), lines.len());
    assert!(!reader.has_ref_alt_columns());
    for (i, line) in lines.into_iter().enumerate() {
        assert_eq!(reader.line_owned(i), line);
    }
    assert_eq!(reader.chr_id(1), 2);
    assert_eq!(reader.position(1), 19);
}

#[test]
fn malformed_legacy_sections_are_rejected_before_typed_access() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("invalid.ktile");
    let (bytes, header, _) = legacy_fixture();
    for mutate in [
        (|h: &mut KtileHeader| h.n_records = u64::MAX) as fn(&mut KtileHeader),
        |h| h.n_records = u64::MAX / 8,
        |h| h.headers_off = u64::MAX,
        |h| h.headers_len = u64::MAX,
        |h| h.line_pool_off = u64::MAX,
        |h| h.line_pool_len = u64::MAX,
        |h| h.line_offsets_off = u64::MAX - 1,
        |h| h.line_offsets_off += 1,
        |h| h.chr_ids_off += 1,
        |h| h.positions_off = u64::MAX,
    ] {
        let mut invalid_header = header;
        mutate(&mut invalid_header);
        let mut invalid = bytes.clone();
        invalid[..KtileHeader::SIZE].copy_from_slice(bytemuck::bytes_of(&invalid_header));
        assert_invalid_ktile(&path, &invalid);
    }
    for (index, offset) in [(0, 1u64), (1, 0), (2, header.line_pool_len + 1)] {
        let mut invalid = bytes.clone();
        let start = header.line_offsets_off as usize + index * 8;
        invalid[start..start + 8].copy_from_slice(&offset.to_le_bytes());
        assert_invalid_ktile(&path, &invalid);
    }
    for length in [
        KtileHeader::SIZE - 1,
        header.line_offsets_off as usize + 7,
        bytes.len() - 1,
    ] {
        assert_invalid_ktile(&path, &bytes[..length]);
    }
}

#[test]
fn malformed_compact_offsets_and_chunk_ranges_are_rejected() {
    use crate::annotate::ktile::writer::{KtileWriteOptions, write_ktile_from_vcf_with};
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("records.vcf");
    let output = dir.path().join("records.ktile");
    let invalid_path = dir.path().join("invalid.ktile");
    let mut content =
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n".to_string();
    for i in 1..=1024 {
        content.push_str(&format!("1\t{i}\t.\tA\tC\t.\t.\t.\n"));
    }
    std::fs::write(&input, content).unwrap();
    write_ktile_from_vcf_with(
        &input,
        &output,
        KtileWriteOptions {
            compressed: true,
            lines_per_chunk: 512,
            deflate_level: 1,
        },
    )
    .unwrap();
    let bytes = std::fs::read(&output).unwrap();
    let header: KtileHeader = bytemuck::pod_read_unaligned(&bytes[..KtileHeader::SIZE]);
    assert!(header.has_elias_fano_offsets());
    for mutate in [
        (|h: &mut KtileHeader| h.version = 1) as fn(&mut KtileHeader),
        |h| h.chr_ids_off = h.line_offsets_off + 3,
        |h| h.line_offsets_off = u64::MAX,
        |h| h.n_records += 1,
        |h| h.lines_per_chunk = 0,
        |h| h.n_chunks += 1,
        |h| h.off_chunk_index += 1,
        |h| h.off_chunk_index = u64::MAX,
    ] {
        let mut invalid_header = header;
        mutate(&mut invalid_header);
        let mut invalid = bytes.clone();
        invalid[..KtileHeader::SIZE].copy_from_slice(bytemuck::bytes_of(&invalid_header));
        assert_invalid_ktile(&invalid_path, &invalid);
    }
    let mut invalid = bytes.clone();
    let offset = header.line_offsets_off as usize;
    invalid[offset..offset + 8].copy_from_slice(&u64::MAX.to_le_bytes());
    assert_invalid_ktile(&invalid_path, &invalid);
    for (relative_offset, value) in [(0, u64::MAX), (8, u64::MAX), (16, u64::MAX)] {
        let mut invalid = bytes.clone();
        let start = header.off_chunk_index as usize + relative_offset;
        invalid[start..start + 8].copy_from_slice(&value.to_le_bytes());
        assert_invalid_ktile(&invalid_path, &invalid);
    }
    assert_invalid_ktile(&invalid_path, &bytes[..header.chr_ids_off as usize - 1]);
}

#[test]
fn compressed_cache_disconnects_blocked_output_before_joining_workers() {
    use crossbeam_channel::SendTimeoutError;
    use std::time::Duration;
    let (output, prefetched) = bounded(0);
    let (finished_tx, finished_rx) = bounded(1);
    let worker = std::thread::spawn(move || {
        let result = output.send_timeout((0, Vec::new()), Duration::from_secs(2));
        finished_tx.send(result).unwrap();
    });
    let cache = CompressedCache {
        entries: Vec::new(),
        decompressor: Decompress::new(false),
        prefetched,
        request_tx: None,
        handles: vec![worker],
        next_to_prefetch: 0,
        n_chunks: 1,
        pending: FxHashMap::default(),
    };
    drop(cache);
    assert!(matches!(
        finished_rx.recv().unwrap(),
        Err(SendTimeoutError::Disconnected(_))
    ));
}

#[test]
fn synchronous_compressed_fallback_keeps_prefetch_progress() {
    use crate::annotate::ktile::writer::{KtileWriteOptions, write_ktile_from_vcf_with};
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("fallback.vcf");
    let output = dir.path().join("fallback.ktile");
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let rows = (1..=64)
        .map(|pos| format!("1\t{pos}\t.\tA\tG\t.\t.\t.\n"))
        .collect::<String>();
    std::fs::write(&input, format!("{header}{rows}")).unwrap();
    write_ktile_from_vcf_with(
        &input,
        &output,
        KtileWriteOptions {
            compressed: true,
            lines_per_chunk: 1,
            deflate_level: 1,
        },
    )
    .unwrap();
    let reader = KtileReader::open(&output).unwrap();
    let PoolStrategy::Compressed { cache, .. } = &reader.strategy else {
        panic!("expected compressed reader");
    };
    let initial = {
        let mut cache = cache.borrow_mut();
        let initial = cache.next_to_prefetch;
        cache.request_tx.take();
        cache.prefetched = crossbeam_channel::never();
        for worker in cache.handles.drain(..) {
            worker.join().unwrap();
        }
        initial
    };
    assert!(initial > 1);
    assert_eq!(reader.line_owned(0), "1\t1\t.\tA\tG\t.\t.\t.");
    assert_eq!(cache.borrow().next_to_prefetch, initial);
}
