use super::*;

#[test]
fn region_parsing() {
    let r = Region::parse("chr1:1,000-2,000").unwrap();
    assert_eq!(
        (r.chr.as_str(), r.start, r.end),
        ("chr1", Some(1000), Some(2000))
    );
    let r = Region::parse("chr1").unwrap();
    assert_eq!(r.bounds(), (1, u32::MAX));
    let r = Region::parse("chr1:100").unwrap();
    assert_eq!(r.bounds(), (100, 100));
    let r = Region::parse_with("chr1:100", false).unwrap();
    assert_eq!(r.bounds(), (100, u32::MAX));
    let r = Region::parse("chr1:100-").unwrap();
    assert_eq!(r.bounds(), (100, u32::MAX));
    let r = Region::parse("HLA-A*01:01:01:01:5-10").unwrap();
    assert_eq!(r.chr, "HLA-A*01:01:01:01");
    assert_eq!(r.bounds(), (5, 10));
    assert_eq!(parse_coordinate("1.5k"), None);
    assert_eq!(parse_coordinate("2k"), Some(2000));
    assert_eq!(parse_coordinate("3M"), Some(3_000_000));
}

#[test]
fn genomic_key_holds_large_contig_ids() {
    let k = GenomicKey::new(70_000, 123);
    assert_eq!(k.chr(), 70_000);
    assert_eq!(k.position(), 123);
}

/// `fmt_g` against htslib's `kputd` (libhts 1.19) on `fixtures/kputd.tsv`:
/// input text, then what bcftools prints for it after the float round trip.
#[test]
fn fmt_g_matches_htslib_kputd() {
    let mut n = 0;
    for line in include_str!("fixtures/kputd.tsv").lines() {
        let line = line.trim_end_matches('\r');
        let Some((input, want)) = line.split_once('\t') else {
            continue;
        };
        let Ok(d) = input.parse::<f64>() else {
            panic!("unparsable oracle input {input}")
        };
        assert_eq!(fmt_g(d as f32 as f64), want, "{input}");
        n += 1;
    }
    assert!(n > 1000, "{n} oracle rows");
    assert_eq!(fmt_g(f64::NAN), "nan");
    assert_eq!(fmt_g(f64::NEG_INFINITY), "-inf");
    assert_eq!(fmt_g(-0.0), "-0");
    assert_eq!(fmt_g(1e100), "1e+100");
    assert_eq!(fmt_g(0.5e-6), "5e-07");
}
