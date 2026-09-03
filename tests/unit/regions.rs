use super::*;

#[test]
fn cli_regions_merge_and_lookup() {
    // Comma lists cannot carry thousands separators (bcftools splits on ',' too).
    let r = RegionSet::from_cli("1:100-200,1:150-300,2,1:1000-1050").unwrap();
    assert_eq!(r.intervals("1").unwrap(), &[(100, 300), (1000, 1050)]);
    assert_eq!(r.intervals("2").unwrap(), &[(1, u32::MAX)]);
    assert!(r.contains("1", 250));
    assert!(!r.contains("1", 301));
    assert!(!r.contains("3", 1));
    assert!(r.overlaps_range("1", 90, 100));
    assert!(!r.overlaps_range("1", 301, 999));
    assert_eq!(r.next_interval("1", 400), Some((1000, 1050)));
    assert_eq!(r.next_interval("1", 1051), None);
    assert_eq!(r.contigs(), &["1".to_string(), "2".to_string()]);
    assert_eq!(r.iter().collect::<Vec<_>>(), vec![("1", 100, 300), ("1", 1000, 1050), ("2", 1, u32::MAX)]);
    assert!(RegionSet::from_cli("1:100-50").is_err());
    assert!(RegionSet::from_cli("1:abc").is_err());
}

#[test]
fn overlap_modes_on_lines() {
    let r = RegionSet::from_cli("1:105-110").unwrap();
    let del = "1\t100\t.\tACGTACGT\tA\t.\t.\t.";
    assert!(!r.line_passes_mode(del, 0));
    assert!(r.line_passes_mode(del, 1));
    assert!(r.line_passes_mode(del, 2));
    let ins = "1\t100\t.\tA\tACGTACGTAC\t.\t.\t.";
    assert!(!r.line_passes_mode(ins, 0));
    assert!(!r.line_passes_mode(ins, 1));
    assert!(r.line_passes_mode(ins, 2));
    let sym = "1\t100\t.\tA\t<DEL>\t.\t.\t.";
    assert!(!r.line_passes_mode(sym, 2));
    assert!(!r.line_passes("garbage"));
}

#[test]
fn bed_and_tab_files() {
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("r.bed");
    std::fs::write(&bed, "track name=x\n1\t99\t200\n1\t150\t300\n#c\n2\t0\t10\n").unwrap();
    let r = RegionSet::from_file(&bed).unwrap();
    assert_eq!(r.intervals("1").unwrap(), &[(100, 300)]);
    assert_eq!(r.intervals("2").unwrap(), &[(1, 10)]);

    let tab = dir.path().join("r.txt");
    std::fs::write(&tab, "1\t100\n1\t200\t250\n3\n").unwrap();
    let r = RegionSet::from_file(&tab).unwrap();
    assert_eq!(r.intervals("1").unwrap(), &[(100, 100), (200, 250)]);
    assert_eq!(r.intervals("3").unwrap(), &[(1, u32::MAX)]);

    let bad = dir.path().join("bad.txt");
    std::fs::write(&bad, "1\tx\n").unwrap();
    assert!(RegionSet::from_file(&bad).is_err());

    let both = RegionSet::from_args(Some("1:1-50"), Some(&tab)).unwrap().unwrap();
    assert_eq!(both.intervals("1").unwrap(), &[(1, 50), (100, 100), (200, 250)]);
    assert!(RegionSet::from_args(None, None).unwrap().is_none());
}
