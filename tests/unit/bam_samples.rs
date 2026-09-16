use super::*;
use noodles_sam::header::record::value::{Map, map::ReadGroup};

fn header(rgs: &[(&str, Option<&str>)]) -> sam::Header {
    let mut b = sam::Header::builder();
    for (id, sm) in rgs {
        let mut rg = Map::<ReadGroup>::default();
        if let Some(s) = sm {
            rg.other_fields_mut()
                .insert((*b"SM").try_into().unwrap(), s.to_string().into());
        }
        b = b.add_read_group(id.to_string(), rg);
    }
    b.build()
}

#[test]
fn one_sample_per_file_uses_default_index() {
    let mut m = SampleMap::default();
    let o = SampleOpts::default();
    let a = m
        .add_bam(&header(&[("rg1", Some("S1"))]), "a.bam", &o)
        .unwrap();
    let b = m
        .add_bam(&header(&[("rg2", Some("S2"))]), "b.bam", &o)
        .unwrap();
    assert_eq!(m.names, vec!["S1", "S2"]);
    assert_eq!(m.sample_id(a, Some(b"anything")), Some(0));
    assert_eq!(m.sample_id(b, None), Some(1));
}

#[test]
fn no_read_groups_names_the_file() {
    let mut m = SampleMap::default();
    let id = m
        .add_bam(&header(&[]), "x.bam", &SampleOpts::default())
        .unwrap();
    assert_eq!(m.names, vec!["x.bam"]);
    assert_eq!(m.sample_id(id, None), Some(0));
    // With -s the file has no usable read group and is dropped.
    let o = SampleOpts::parse(Some("S1"), None, None, false).unwrap();
    let mut m = SampleMap::default();
    assert!(m.add_bam(&header(&[]), "x.bam", &o).is_none());
}

#[test]
fn exclusion_and_renaming() {
    let o = SampleOpts::parse(Some("^S2"), None, None, false).unwrap();
    let mut m = SampleMap::default();
    let f = m
        .add_bam(
            &header(&[("rg1", Some("S1")), ("rg2", Some("S2"))]),
            "a.bam",
            &o,
        )
        .unwrap();
    assert_eq!(m.names, vec!["S1"]);
    assert_eq!(m.sample_id(f, Some(b"rg1")), Some(0));
    assert_eq!(m.sample_id(f, Some(b"rg2")), None);

    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("s.txt");
    std::fs::write(&p, "S2 NEW2\n").unwrap();
    let o = SampleOpts::parse(None, Some(&p), None, false).unwrap();
    let mut m = SampleMap::default();
    let f = m
        .add_bam(
            &header(&[("rg1", Some("S1")), ("rg2", Some("S2"))]),
            "a.bam",
            &o,
        )
        .unwrap();
    assert_eq!(m.names, vec!["NEW2"]);
    assert_eq!(m.sample_id(f, Some(b"rg2")), Some(0));
    assert_eq!(m.sample_id(f, Some(b"rg1")), None);
}

#[test]
fn read_group_rows_rename_and_order_by_first_appearance() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("rg.txt");
    std::fs::write(&p, "rgA\nrgB\tSB\nrgC\ta.bam\tSC\n").unwrap();
    let o = SampleOpts::parse(None, None, Some(&p), false).unwrap();
    let mut m = SampleMap::default();
    let f = m
        .add_bam(
            &header(&[
                ("rgB", Some("S1")),
                ("rgA", Some("S1")),
                ("rgC", Some("S3")),
                ("rgD", Some("S4")),
            ]),
            "a.bam",
            &o,
        )
        .unwrap();
    assert_eq!(m.names, vec!["SB", "S1", "SC"]);
    assert_eq!(m.sample_id(f, Some(b"rgA")), Some(1));
    assert_eq!(m.sample_id(f, Some(b"rgB")), Some(0));
    assert_eq!(m.sample_id(f, Some(b"rgC")), Some(2));
    assert_eq!(m.sample_id(f, Some(b"rgD")), None);
    assert_eq!(m.sample_id(f, None), None);
}
