    use super::*;

    fn hdr() -> HeaderInfo {
        HeaderInfo::parse(&[
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">",
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"x\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"x\">",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"x\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"x\">",
        ])
    }

    fn opts() -> SplitOpts {
        SplitOpts { missing_for_overlap: false, keep_sum_keys: Vec::new(), old_rec_tag: None }
    }

    #[test]
    fn split_rewrites_pl_ad_and_gt() {
        let line = "1\t100\t.\tA\tC,G\t50\tPASS\tAC=1,2;AN=4;AF=0.25,0.5\tGT:AD:PL:DP\t0/1:5,3,0:10,0,20,30,40,50:8\t1/2:0,2,4:60,50,40,30,20,0:6";
        let recs = expand_record(line, MultiMode::SplitAll, false, &opts(), &hdr()).unwrap();
        assert_eq!(recs.len(), 2);
        // 1/2 becomes 1/0 for the first ALT: allele order is kept, like bcftools.
        assert_eq!(recs[0], "1\t100\t.\tA\tC\t50\tPASS\tAC=1;AN=4;AF=0.25\tGT:AD:PL:DP\t0/1:5,3:10,0,20:8\t1/0:0,2:60,50,40:6");
        assert_eq!(recs[1], "1\t100\t.\tA\tG\t50\tPASS\tAC=2;AN=4;AF=0.5\tGT:AD:PL:DP\t0/0:5,0:10,30,50:8\t0/1:0,4:60,30,0:6");
    }

    #[test]
    fn split_keep_sum_and_old_rec_tag() {
        let mut o = opts();
        o.keep_sum_keys = vec!["AD".into()];
        o.old_rec_tag = Some("OLD".into());
        let line = "1\t100\t.\tA\tC,G\t.\t.\t.\tGT:AD\t1/2:1,2,4";
        let recs = expand_record(line, MultiMode::SplitAll, false, &o, &hdr()).unwrap();
        assert_eq!(recs[0], "1\t100\t.\tA\tC\t.\t.\tOLD=1|100|A|C,G|1\tGT:AD\t1/0:5,2");
        assert_eq!(recs[1], "1\t100\t.\tA\tG\t.\t.\tOLD=1|100|A|C,G|2\tGT:AD\t0/1:3,4");
    }

    #[test]
    fn split_overlap_allele_to_missing() {
        let mut o = opts();
        o.missing_for_overlap = true;
        let line = "1\t100\t.\tA\tC,*\t.\t.\t.\tGT\t1/2";
        let recs = expand_record(line, MultiMode::SplitAll, false, &o, &hdr()).unwrap();
        assert_eq!(recs[0], "1\t100\t.\tA\tC\t.\t.\t.\tGT\t1/.");
        assert_eq!(recs[1], "1\t100\t.\tA\t*\t.\t.\t.\tGT\t0/1");
    }

    #[test]
    fn join_biallelic_records() {
        let h = hdr();
        let mut j = Joiner::new(MultiMode::JoinAll);
        let a: Vec<String> = "1\t100\trs1\tA\tC\t50\tPASS\tAC=1;AN=4\tGT:AD:PL\t0/1:5,3:10,0,20\t0/0:6,0:0,10,20".split('\t').map(String::from).collect();
        let b: Vec<String> = "1\t100\t.\tA\tG\t70\tq10\tAC=1;AN=4\tGT:AD:PL\t0/0:5,0:0,10,20\t0/1:6,4:10,0,20".split('\t').map(String::from).collect();
        assert!(j.push(a, &h).is_empty());
        assert!(j.push(b, &h).is_empty());
        let out = j.finish(&h);
        assert_eq!(out.len(), 1);
        assert_eq!(
            out[0].join("\t"),
            "1\t100\trs1\tA\tC,G\t70\tq10\tAC=1,1;AN=4\tGT:AD:PL\t0/1:5,3,0:10,0,20,10,.,20\t0/2:6,0,4:0,10,20,0,.,20"
        );
    }

    #[test]
    fn join_pads_shorter_ref() {
        let h = hdr();
        let mut j = Joiner::new(MultiMode::JoinAll);
        let a: Vec<String> = "1\t100\t.\tAT\tA\t.\t.\t.\tGT\t0/1".split('\t').map(String::from).collect();
        let b: Vec<String> = "1\t100\t.\tA\tC\t.\t.\t.\tGT\t0/1".split('\t').map(String::from).collect();
        j.push(a, &h);
        j.push(b, &h);
        let out = j.finish(&h);
        assert_eq!(out[0].join("\t"), "1\t100\t.\tAT\tA,CT\t.\t.\t.\tGT\t1/2");
    }

    #[test]
    fn join_only_snps_leaves_indels() {
        let h = hdr();
        let mut j = Joiner::new(MultiMode::JoinSnps);
        let a: Vec<String> = "1\t100\t.\tA\tC\t.\t.\t.".split('\t').map(String::from).collect();
        let b: Vec<String> = "1\t100\t.\tA\tAT\t.\t.\t.".split('\t').map(String::from).collect();
        let c: Vec<String> = "1\t100\t.\tA\tG\t.\t.\t.".split('\t').map(String::from).collect();
        let mut out = Vec::new();
        out.extend(j.push(a, &h));
        out.extend(j.push(b, &h));
        out.extend(j.push(c, &h));
        out.extend(j.finish(&h));
        let lines: Vec<String> = out.iter().map(|c| c.join("\t")).collect();
        assert_eq!(lines, vec!["1\t100\t.\tA\tC\t.\t.\t.", "1\t100\t.\tA\tAT\t.\t.\t.", "1\t100\t.\tA\tG\t.\t.\t."]);
    }

    #[test]
    fn fix_ref_swaps_alleles() {
        let h = hdr();
        let mut cols: Vec<String> = "1\t100\t.\tA\tC,G\t.\t.\tAC=1,2;AN=4\tGT:AD:PL\t1/2:1,2,4:0,1,2,3,4,5".split('\t').map(String::from).collect();
        fix_ref(&mut cols, "G", &h);
        // New genotype order over [G,C,A]: 00=old 22, 01=old 12, 11=old 11, 02=old 02, 12=old 01, 22=old 00.
        assert_eq!(cols.join("\t"), "1\t100\t.\tG\tC,A\t.\t.\tAC=1,.;AN=4\tGT:AD:PL\t1/0:4,2,1:5,4,2,3,1,0");
    }

    #[test]
    fn dedup_by_type() {
        let mut w = DedupWindow::default();
        let a: Vec<String> = "1\t100\t.\tA\tC\t.\t.\t.".split('\t').map(String::from).collect();
        let b: Vec<String> = "1\t100\t.\tA\tG\t.\t.\t.".split('\t').map(String::from).collect();
        let c: Vec<String> = "1\t100\t.\tA\tAT\t.\t.\t.".split('\t').map(String::from).collect();
        w.push(a.clone());
        assert!(w.is_dup(&b, RmDup::Snps));
        assert!(!w.is_dup(&b, RmDup::Exact));
        assert!(w.is_dup(&a, RmDup::Exact));
        assert!(!w.is_dup(&c, RmDup::Snps));
        assert!(w.is_dup(&c, RmDup::All));
    }

    #[test]
    fn filters_merge_like_bcftools() {
        assert_eq!(merge_filters(["PASS", "PASS"].into_iter()), "PASS");
        assert_eq!(merge_filters(["PASS", "q10"].into_iter()), "q10");
        assert_eq!(merge_filters([".", "."].into_iter()), ".");
        assert_eq!(merge_filters(["a;b", "b;c"].into_iter()), "a;b;c");
    }
