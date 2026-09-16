use super::*;
use std::collections::BTreeSet;

fn lcg(seed: &mut u64) -> u32 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 33) as u32
}

#[test]
fn membership_and_ranges_match_a_set_oracle() {
    let mut seed = 7u64;
    let mut oracle = BTreeSet::new();
    let mut c = ContigBitmap::default();
    // Sparse positions across many chunks plus one dense chunk that must
    // turn into a bitset.
    for _ in 0..3000 {
        let p = lcg(&mut seed) % 5_000_000;
        oracle.insert(p);
        c.insert(p);
    }
    for _ in 0..6000 {
        let p = 1_000_000 + lcg(&mut seed) % 60_000;
        oracle.insert(p);
        c.insert(p);
    }
    assert_eq!(c.len(), oracle.len());
    assert!(matches!(
        c.chunks.get(&((1_000_000u32 >> 16) as u16)),
        Some(Container::Bitset(_))
    ));
    for _ in 0..20_000 {
        let p = lcg(&mut seed) % 5_100_000;
        assert_eq!(c.contains(p), oracle.contains(&p), "contains {p}");
    }
    for _ in 0..20_000 {
        let a = lcg(&mut seed) % 5_100_000;
        let b = a + lcg(&mut seed) % 3000;
        let want = oracle.range(a..=b).next().is_some();
        assert_eq!(c.any_in_range(a, b), want, "range {a}..={b}");
    }
    assert!(!c.any_in_range(10, 5));
    let got: Vec<u32> = c.iter().collect();
    let want: Vec<u32> = oracle.iter().copied().collect();
    assert_eq!(got, want);
}

#[test]
fn set_algebra() {
    let mut a = ContigBitmap::default();
    let mut b = ContigBitmap::default();
    for p in [1u32, 2, 3, 70_000, 70_001, 200_000] {
        a.insert(p);
    }
    for p in [2u32, 3, 4, 70_001, 300_000] {
        b.insert(p);
    }
    let and: Vec<u32> = a.and(&b).iter().collect();
    assert_eq!(and, vec![2, 3, 70_001]);
    let or: Vec<u32> = a.or(&b).iter().collect();
    assert_eq!(or, vec![1, 2, 3, 4, 70_000, 70_001, 200_000, 300_000]);
    let diff: Vec<u32> = a.and_not(&b).iter().collect();
    assert_eq!(diff, vec![1, 70_000, 200_000]);
    // Empty results drop the chunk entirely.
    assert!(a.and(&ContigBitmap::default()).is_empty());
}

#[test]
fn per_contig_operations_and_roundtrip() {
    let mut x = PosBitmap::default();
    let mut y = PosBitmap::default();
    for (c, p) in [("chr1", 10u32), ("chr1", 20), ("chr2", 5), ("chrX", 1)] {
        x.insert(c, p);
    }
    for (c, p) in [("chr1", 20u32), ("chr2", 6), ("chrY", 9)] {
        y.insert(c, p);
    }
    assert_eq!(x.len(), 4);
    assert!(x.contains("chr1", 20));
    assert!(!x.contains("chr3", 20));
    assert!(x.any_in_range("chr2", 1, 5));
    assert!(!x.any_in_range("chr2", 6, 100));
    let and = x.and(&y);
    assert_eq!(and.len(), 1);
    assert!(and.contains("chr1", 20));
    assert!(and.contig("chr2").is_none());
    assert_eq!(x.or(&y).len(), 6);
    let d = x.and_not(&y);
    assert_eq!(d.len(), 3);
    assert!(d.contains("chrX", 1) && !d.contains("chr1", 20));

    // REF spans: a 4-base deletion at chr2:100 covers 101..=103 as well.
    x.insert_span("chr2", 100, 4);
    assert_eq!(x.tail_len(), 3);
    assert!(x.any_span_in_range("chr2", 102, 110));
    assert!(!x.any_in_range("chr2", 102, 110));
    assert!(!x.any_span_in_range("chr2", 104, 110));
    // Spans survive a union but not an intersection.
    assert_eq!(x.or(&y).tail_len(), 3);
    assert_eq!(x.and(&x).tail_len(), 0);

    // A dense contig survives the file format too.
    let mut dense = ContigBitmap::default();
    for p in 0..5000u32 {
        dense.insert(p * 3);
    }
    x.contigs.insert("dense".to_string(), dense);
    let dir = std::env::temp_dir().join(format!("kbm_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("x.kbm");
    x.save(&path).unwrap();
    let back = PosBitmap::load(&path).unwrap();
    assert_eq!(back, x);
    assert_eq!(
        PosBitmap::sidecar_path(std::path::Path::new("a/b.vcf.gz")),
        std::path::PathBuf::from("a/b.vcf.kbm")
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn from_vcf_reads_every_record_position_and_span() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let path = root.join("tests/filter/test31/in.vcf");
    if !path.exists() {
        return;
    }
    let bm = PosBitmap::from_vcf(&path).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let mut want = BTreeSet::new();
    let mut tails = BTreeSet::new();
    for line in text.lines() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let p: u32 = f[1].parse().unwrap();
        want.insert((f[0].to_string(), p));
        for q in p + 1..p + f[3].len() as u32 {
            tails.insert((f[0].to_string(), q));
        }
    }
    let got: BTreeSet<(String, u32)> = bm
        .contig_names()
        .flat_map(|c| {
            bm.contig(c)
                .unwrap()
                .iter()
                .map(move |p| (c.to_string(), p))
        })
        .collect();
    assert_eq!(got, want);
    assert!(!tails.is_empty());
    for (c, q) in &tails {
        assert!(bm.any_span_in_range(c, *q, *q), "{c}:{q} spanned");
    }
    assert_eq!(bm.tail_len(), tails.len());
}
