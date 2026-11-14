use super::*;

#[test]
fn likelihood_blocks_preserve_samples_and_select_pl_by_heterozygous_likelihood() {
    let mut blocker = GvcfBlocker::from_spec("0,2,5").unwrap();
    let mut output = Vec::new();
    blocker.write_likelihoods("1\t100\t.\tA\t<*>\t0\t.\tQS=2,0\tPL:DP:DV\t0,30,90:10:0\t0,20,80:3:0\n1\t101\t.\tC\t<*>\t0\t.\tQS=2,0\tPL:DP:DV\t0,12,120:4:0\t0,20,60:2:0\n1\t102\t.\tG\tT,<*>\t0\t.\tDP=8\tPL:DP:DV\t10,0,20:4:2\t0,10,20:4:0\n", &mut output).unwrap();
    blocker.flush(&mut output).unwrap();
    assert_eq!(
        String::from_utf8(output).unwrap(),
        "1\t100\t.\tA\t<*>\t.\t.\tEND=101;MIN_DP=2;QS=2,0\tPL:DP\t0,12,120:4\t0,20,60:2\n1\t102\t.\tG\tT,<*>\t0\t.\tDP=8\tPL:DP:DV\t10,0,20:4:2\t0,10,20:4:0\n"
    );
}

#[test]
fn malformed_depth_thresholds_are_rejected() {
    for spec in ["", "1,x", "5,2", "2,2", "1,"] {
        assert!(GvcfBlocker::from_spec(spec).is_err(), "{spec}");
    }
}

#[test]
fn collapses_adjacent_ref_sites_same_bin() {
    let mut b = GvcfBlocker::new(vec![0, 10, 20]);
    let mut out: Vec<u8> = Vec::new();
    b.add_ref_site("1", 100, "A", 15, 50.0, &mut out).unwrap();
    b.add_ref_site("1", 101, "C", 16, 60.0, &mut out).unwrap();
    b.add_ref_site("1", 102, "G", 17, 70.0, &mut out).unwrap();
    b.flush(&mut out).unwrap();
    let s = String::from_utf8(out).unwrap();
    assert!(s.contains("END=102"));
    assert!(s.contains("MIN_DP=15"));
    assert!(s.contains("<NON_REF>"));
}

#[test]
fn flushes_on_bin_change() {
    let mut b = GvcfBlocker::new(vec![0, 10, 20]);
    let mut out: Vec<u8> = Vec::new();
    b.add_ref_site("1", 100, "A", 15, 50.0, &mut out).unwrap();
    b.add_ref_site("1", 101, "C", 25, 60.0, &mut out).unwrap();
    b.flush(&mut out).unwrap();
    let s = String::from_utf8(out).unwrap();
    let blocks: Vec<&str> = s.lines().collect();
    assert_eq!(blocks.len(), 2, "expected 2 blocks, got {:?}", blocks);
}

#[test]
fn flushes_on_chrom_change() {
    let mut b = GvcfBlocker::new(vec![0]);
    let mut out: Vec<u8> = Vec::new();
    b.add_ref_site("1", 100, "A", 10, 50.0, &mut out).unwrap();
    b.add_ref_site("2", 100, "A", 10, 50.0, &mut out).unwrap();
    b.flush(&mut out).unwrap();
    let s = String::from_utf8(out).unwrap();
    assert_eq!(s.lines().count(), 2);
}
