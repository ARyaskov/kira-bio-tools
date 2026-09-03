    use super::*;

    fn hdr() -> HeaderInfo {
        HeaderInfo::parse(&[
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"x\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"x\">",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"x\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"x\">",
        ])
    }

    #[test]
    fn split_second_alt() {
        // alleles [REF, A1, A2] -> keep [REF, A2]
        let map = [Some(0), None, Some(1)];
        assert_eq!(remap_value("5,7", FieldNumber::A, 3, 2, &map).unwrap(), "7");
        assert_eq!(remap_value("10,5,7", FieldNumber::R, 3, 2, &map).unwrap(), "10,7");
        // diploid G over 3 alleles: 00,01,11,02,12,22 -> 00,02,22
        assert_eq!(remap_value("0,10,20,30,40,50", FieldNumber::G, 3, 2, &map).unwrap(), "0,30,50");
        assert_eq!(remap_gt("0/2", &map), "0/1");
        assert_eq!(remap_gt("1|2", &map), ".|1");
        assert_eq!(remap_gt("./.", &map), "./.");
        let h = hdr();
        assert_eq!(remap_info("AC=5,7;DP=20;AF=0.1,0.2", &h, 3, 2, &map), "AC=7;DP=20;AF=0.2");
        let s = remap_samples("GT:AD:PL:DP", &["1/2:1,2,3:0,1,2,3,4,5:9"], &h, 3, 2, &map);
        assert_eq!(s, vec!["./1:1,3:0,3,5:9"]);
    }

    #[test]
    fn expand_to_union() {
        // old [REF, T] -> new [REF, G, T]
        let map = [Some(0), Some(2)];
        assert_eq!(remap_value("3", FieldNumber::A, 2, 3, &map).unwrap(), ".,3");
        assert_eq!(remap_value("0,10,20", FieldNumber::G, 2, 3, &map).unwrap(), "0,.,.,10,.,20");
        assert_eq!(remap_gt("0/1", &map), "0/2");
    }

    #[test]
    fn wrong_cardinality_is_left_alone() {
        let map = [Some(0), Some(1)];
        assert!(remap_value("1,2,3", FieldNumber::A, 2, 2, &map).is_none());
        assert_eq!(remap_value(".", FieldNumber::A, 2, 2, &map).unwrap(), ".");
        assert_eq!(remap_value("0", FieldNumber::G, 1, 1, &[Some(0)]).unwrap(), "0");
    }
