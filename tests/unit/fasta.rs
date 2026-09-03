use super::*;

#[test]
fn indexed_and_scanned_access_agree() {
    let dir = tempfile::tempdir().unwrap();
    let fa = dir.path().join("ref.fa");
    std::fs::write(&fa, ">chr1 desc\nACGTAC\nGTACGT\nAC\n>chr2\nTTTT\nGG\n").unwrap();
    let scanned = IndexedFasta::open(&fa).unwrap();
    assert_eq!(scanned.contig("chr1").unwrap(), b"ACGTACGTACGTAC");
    assert_eq!(scanned.contig("chr2").unwrap(), b"TTTTGG");
    assert_eq!(scanned.base("chr1", 7), Some(b'G'));
    assert_eq!(scanned.slice("chr1", 13, 2), Some(&b"AC"[..]));
    assert_eq!(scanned.slice("chr1", 13, 3), None);
    assert_eq!(scanned.slice_bytes("chr1", 13, 3), Some(&b"AC"[..]));
    assert!(scanned.has("chr2") && !scanned.has("chr3"));

    // The same file with a .fai: contigs are read by offset.
    std::fs::write(fai_path(&fa), "chr1\t14\t11\t6\t7\nchr2\t6\t34\t4\t5\n").unwrap();
    let mut indexed = IndexedFasta::open(&fa).unwrap();
    assert_eq!(indexed.contig("chr1").unwrap(), b"ACGTACGTACGTAC");
    assert_eq!(indexed.contig("chr2").unwrap(), b"TTTTGG");
    assert_eq!(indexed.length("chr2"), Some(6));
    indexed.evict_except("chr2");
    assert!(indexed.seqs[0].get().is_none());
    assert_eq!(indexed.contig("chr1").unwrap(), b"ACGTACGTACGTAC");
}

#[test]
fn irregular_layout_is_kept_in_memory() {
    let dir = tempfile::tempdir().unwrap();
    let fa = dir.path().join("odd.fa");
    std::fs::write(&fa, ">c\nACG\nTACGTA\nC\n").unwrap();
    let f = IndexedFasta::open(&fa).unwrap();
    assert_eq!(f.contig("c").unwrap(), b"ACGTACGTAC");
}

#[test]
fn gzip_input_loads_whole() {
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let fa = dir.path().join("ref.fa.gz");
    let mut enc = flate2::write::GzEncoder::new(std::fs::File::create(&fa).unwrap(), flate2::Compression::default());
    enc.write_all(b">x\nacgt\n>y\nGG\n").unwrap();
    enc.finish().unwrap();
    let f = IndexedFasta::open(&fa).unwrap();
    assert_eq!(f.contig("x").unwrap(), b"ACGT");
    assert_eq!(f.names().collect::<Vec<_>>(), vec!["x", "y"]);
}
