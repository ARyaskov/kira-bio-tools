use super::build_ani_index_from_tab;
use crate::annotate::structs::ani::AniIndex;
use std::fs;

#[test]
fn test_tab_infers_numbers_for_split_info() {
    let dir = std::env::temp_dir();
    let tab_path = dir.join("kira_tab_infer_test.tab");
    let ani_path = dir.join("kira_tab_infer_test.ani");

    let tab = "1\t1\tC\tA,T\t0,1.1\t1.1,0,2.2\n";
    fs::write(&tab_path, tab).unwrap();

    build_ani_index_from_tab(&tab_path, &ani_path, Some("CHROM,POS,REF,ALT,FA,FR")).unwrap();

    let ani = AniIndex::open(&ani_path).unwrap();
    let bundle = ani.lookup_exact("1", 1, "C", "T").unwrap();

    let fa = bundle
        .info
        .iter()
        .find(|f| f.key == "FA")
        .map(|f| f.values.clone())
        .unwrap();
    let fr = bundle
        .info
        .iter()
        .find(|f| f.key == "FR")
        .map(|f| f.values.clone())
        .unwrap();

    assert_eq!(fa, vec!["1.1"]);
    assert_eq!(fr, vec!["1.1", "2.2"]);

    let _ = fs::remove_file(&tab_path);
    let _ = fs::remove_file(&ani_path);
}

/// A record written with extra padding, or anchored on the other side of an
/// indel, still finds its annotation, and the annotation's alleles come back
/// in the target's own frame.
#[test]
fn test_canonical_allele_lookup() {
    let dir = std::env::temp_dir();
    let tab_path = dir.join(format!("kira_canon_{}.tab", std::process::id()));
    let ani_path = dir.join(format!("kira_canon_{}.ani", std::process::id()));
    let tab = "1	101	A	C	0.5
1	200	GATC	GC	0.25
1	300	A	AT,G	0.1,0.2
";
    fs::write(&tab_path, tab).unwrap();
    build_ani_index_from_tab(&tab_path, &ani_path, Some("CHROM,POS,REF,ALT,FA")).unwrap();
    let ani = AniIndex::open(&ani_path).unwrap();
    let chr = ani.contig_id("1").unwrap();
    let meta = std::collections::HashMap::new();
    let hash = |s: &str| crate::util::fast_hash64(s.as_bytes());
    // GAT>GCT at 100 is A>C at 101.
    let b = ani
        .lookup_exact_by_chr_id_pos_index_opts(
            chr,
            100,
            "GAT",
            hash("GAT"),
            "GCT",
            &meta,
            true,
            false,
        )
        .expect("padded SNP matches");
    assert_eq!(b.alt, "GCT");
    assert_eq!(b.db_ref, "GAT");
    assert_eq!(
        b.info
            .iter()
            .find(|f| f.key == "FA")
            .map(|f| f.values.clone()),
        Some(vec!["0.5".to_string()])
    );
    // ATC>C at 201 is the deletion GATC>GC at 200 anchored on the right.
    let b = ani
        .lookup_exact_by_chr_id_pos_index_opts(
            chr,
            201,
            "ATC",
            hash("ATC"),
            "C",
            &meta,
            true,
            false,
        )
        .expect("right-anchored deletion matches");
    assert_eq!(b.alt, "C");
    // GA>GAT at 299 is the insertion A>AT at 300 with a base of left padding
    // (REF-suffix padding alone is vcmp's job); the tab builder stores one
    // entry per ALT, so the matched allele comes back alone, in the target frame.
    let b = ani
        .lookup_exact_by_chr_id_pos_index_opts(
            chr,
            299,
            "GA",
            hash("GA"),
            "GAT",
            &meta,
            true,
            false,
        )
        .expect("left-padded insertion matches");
    assert_eq!(b.alt, "GAT");
    assert_eq!(b.db_ref, "GA");
    // Unrelated alleles do not match.
    assert!(
        ani.lookup_exact_by_chr_id_pos_index_opts(
            chr,
            100,
            "GAT",
            hash("GAT"),
            "GTT",
            &meta,
            true,
            false
        )
        .is_none()
    );
    let _ = fs::remove_file(&tab_path);
    let _ = fs::remove_file(&ani_path);
}
