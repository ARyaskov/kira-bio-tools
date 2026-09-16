use super::*;

#[test]
fn borrowed_record_splits_the_line() {
    let line = "chr1\t12345\trs1\tA\tG,T\t50\tPASS\tDP=10;AF=0.5\tGT:DP\t0/1:7\t1/1:3\r\n";
    let r = VcfRecordRef::parse(line).unwrap();
    assert_eq!(r.chrom, "chr1");
    assert_eq!(r.pos, 12345);
    assert_eq!(r.id, "rs1");
    assert_eq!(r.ref_allele, "A");
    assert_eq!(r.alt, "G,T");
    assert_eq!(r.qual, "50");
    assert_eq!(r.filter, "PASS");
    assert_eq!(r.info, "DP=10;AF=0.5");
    assert_eq!(r.format, Some("GT:DP"));
    assert_eq!(r.samples, vec!["0/1:7", "1/1:3"]);

    let sites_only = VcfRecordRef::parse("1\t5\t.\tA\tC\t.\t.\t.").unwrap();
    assert_eq!(sites_only.format, None);
    assert!(sites_only.samples.is_empty());
    assert!(VcfRecordRef::parse("#CHROM\tPOS").is_none());
    assert!(VcfRecordRef::parse("1\tx\t.\tA\tC\t.\t.\t.").is_none());
    assert!(VcfRecordRef::parse("1\t5\t.\tA\tC").is_none());
}

#[test]
fn from_record_and_from_cols_agree() {
    let line = "2\t10\t.\tC\tT\t.\tq10\tAC=1\tGT\t0|1";
    let owned = VcfRecord {
        chrom: "2".into(),
        pos: 10,
        id: ".".into(),
        ref_allele: "C".into(),
        alt: "T".into(),
        qual: ".".into(),
        filter: "q10".into(),
        info: "AC=1".into(),
        format: Some("GT".into()),
        samples: vec!["0|1".into()],
        chr_id: 0,
        position: 10,
        offset: 0,
    };
    let a = VcfRecordRef::from_record(&owned);
    let cols: Vec<&str> = line.split('\t').collect();
    let b = VcfRecordRef::from_cols(&cols).unwrap();
    let c = VcfRecordRef::parse(line).unwrap();
    for r in [&a, &b, &c] {
        assert_eq!(
            (
                r.chrom,
                r.pos,
                r.filter,
                r.info,
                r.format,
                r.samples.as_slice()
            ),
            ("2", 10, "q10", "AC=1", Some("GT"), &["0|1"][..])
        );
    }
    assert!(VcfRecordRef::from_cols(&cols[..7]).is_none());
}
