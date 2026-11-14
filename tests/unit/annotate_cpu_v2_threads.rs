use super::*;

#[test]
fn reused_batch_scratch_matches_scalar_after_growth_shrink_and_empty() {
    let dir = tempfile::tempdir().unwrap();
    let db = dir.path().join("db.tab");
    let path = dir.path().join("db.ani");
    std::fs::write(&db, "1\t10\tA\tC,G\t0.1,0.2\n1\t20\tT\tA\t0.3\n").unwrap();
    crate::annotate::builder_v2::build_ani_index_from_tab(&db, &path, Some("CHROM,POS,REF,ALT,FA"))
        .unwrap();
    let ani = Arc::new(AniIndex::open(&path).unwrap());
    let meta = Arc::new(HashMap::from([("FA".to_string(), FieldNumber::A)]));
    let modes = Arc::new(vec![("FA".to_string(), AnnotateMode::default_mode())]);
    let samples = Arc::new(Vec::new());
    let lines = [
        "1\t10\t.\tA\tC,G\t.\t.\t.",
        "1\t20\t.\tT\tA\t.\t.\t.",
        "1\t21\t.\tT\tA\t.\t.\t.",
        "unknown\t10\t.\tA\tC\t.\t.\t.",
        "broken",
    ];
    let mut scratch = BatchScratch::default();
    let mut max_capacity = 0;
    for count in [37, 1, 0, 19, 0, 2] {
        let mut batch = ReadBatch::with_capacity(0, count);
        for i in 0..count {
            batch.push_line(lines[i % lines.len()]);
        }
        let result = annotate_batch_simd(
            &batch,
            &ani,
            &meta,
            &modes,
            &samples,
            false,
            false,
            false,
            &mut scratch,
        );
        for (i, actual) in result.iter().enumerate() {
            let expected =
                annotate_line(batch.line(i), &ani, &meta, &modes, &samples, false, false);
            assert_eq!(actual.as_deref().unwrap_or(batch.line(i)), expected);
        }
        max_capacity = max_capacity.max(scratch.flat_keys.capacity());
        assert_eq!(scratch.flat_keys.capacity(), max_capacity);
    }
}
