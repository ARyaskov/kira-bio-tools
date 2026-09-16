//! Property tests (proptest) for the parsers and formatters that face
//! untrusted bytes: they must never panic, round trips must hold, and
//! canonical forms must be fixed points.

use proptest::prelude::*;

use crate::bcf::header::parse_header_to_dict;
use crate::bcf::record::{Undefined, decode_record_to_vcf, encode_record_opts};
use crate::csi::{TabixConf, tabix_line_interval};
use crate::util::fmt_g;
use crate::vcf::{RecordCanon, canonical_float, canonical_floats};

fn dict() -> crate::bcf::BcfHeaderDict {
    parse_header_to_dict(&[
        "##fileformat=VCFv4.2".to_string(),
        "##contig=<ID=1>".to_string(),
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">".to_string(),
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">".to_string(),
        "##FILTER=<ID=q10,Description=\"q\">".to_string(),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">".to_string(),
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"r\">".to_string(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2".to_string(),
    ])
}

proptest! {
    /// `fmt_g` keeps six significant digits and its output is already canonical.
    #[test]
    fn fmt_g_keeps_six_digits(v in any::<f32>().prop_filter("finite", |v| v.is_finite())) {
        let s = fmt_g(v as f64);
        let back: f64 = s.parse().unwrap();
        let v = v as f64;
        prop_assert!((v == 0.0 && back == 0.0) || (back - v).abs() <= v.abs() * 5e-6, "{v} -> {s}");
        prop_assert!(matches!(canonical_float(&s), std::borrow::Cow::Borrowed(_)), "{s} not canonical");
    }

    /// Canonicalising float text is idempotent, for scalars and vectors.
    #[test]
    fn canonical_float_is_a_fixed_point(s in "[+-]?[0-9]{0,8}(\\.[0-9]{0,8})?([eE][+-]?[0-9]{1,3})?") {
        let once = canonical_float(&s).into_owned();
        prop_assert_eq!(canonical_float(&once).into_owned(), once.clone());
        let list = format!("{s},.,{s}");
        let once = canonical_floats(&list).into_owned();
        prop_assert_eq!(canonical_floats(&once).into_owned(), once);
    }

    /// A record survives the BCF encode/decode round trip byte for byte.
    #[test]
    fn bcf_record_round_trip(
        pos in 1u32..100_000_000,
        refa in "[ACGT]{1,6}",
        alt in "[ACGT]{1,6}",
        dp in 0i32..100_000,
        af in 0.0f32..1.0,
        gt in "[01][/|][01]",
        ad in (0u32..500, 0u32..500),
    ) {
        let dict = dict();
        let af = fmt_g(af as f64);
        let line = format!("1\t{pos}\t.\t{refa}\t{alt}\t50\tq10\tDP={dp};AF={af}\tGT:AD\t{gt}:{},{}\t0/0:.", ad.0, ad.1);
        let mut buf = Vec::new();
        encode_record_opts(&mut buf, &line, &dict, Undefined::Error).unwrap();
        let out = decode_record_to_vcf(&mut std::io::Cursor::new(buf), &dict).unwrap().unwrap();
        prop_assert_eq!(out, line);
    }

    /// The text canonicaliser never panics and is a fixed point.
    #[test]
    fn record_canon_never_panics(line in "\\PC{0,80}") {
        let mut c = RecordCanon::from_headers(&[
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">".to_string(),
            "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"g\">".to_string(),
        ]);
        if let Some(out) = c.canon_line(&line) {
            prop_assert_eq!(c.canon_line(&out), None);
        }
    }

    /// Tabix line parsing and genotype parsing accept any text.
    #[test]
    fn line_parsers_never_panic(line in "\\PC{0,60}") {
        for conf in [TabixConf::VCF, TabixConf::BED, TabixConf::GFF, TabixConf::SAM] {
            let _ = tabix_line_interval(&conf, &line);
        }
        let _ = crate::consensus::parse_gt(&line);
        let _ = crate::vcf::parse_vcf_line(&line);
    }
}
