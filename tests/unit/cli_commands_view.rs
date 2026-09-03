    use super::*;

    fn hdr() -> HeaderInfo {
        HeaderInfo::parse(&[
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">",
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"x\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"x\">",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"x\">",
        ])
    }

    #[test]
    fn trim_rewrites_allele_indexed_fields() {
        let h = hdr();
        let line = "1\t10\t.\tA\tC,G\t.\t.\tAC=1,0;AN=4;AF=0.25,0\tGT:AD:PL\t0/1:5,3,0:10,0,20,30,40,50\t0/0:6,0,0:0,10,20,30,40,50";
        let out = trim_alt_alleles(line, &h, true, false);
        assert_eq!(
            out,
            "1\t10\t.\tA\tC\t.\t.\tAC=1;AN=4;AF=0.25\tGT:AD:PL\t0/1:5,3:10,0,20\t0/0:6,0:0,10,20"
        );
        // Nothing observed but REF: ALT becomes '.'
        let line2 = "1\t10\t.\tA\tC\t.\t.\tAC=0;AN=2\tGT\t0/0";
        assert_eq!(trim_alt_alleles(line2, &h, true, false), "1\t10\t.\tA\t.\t.\t.\tAN=2\tGT\t0/0");
        // --trim-unseen-allele only touches <*>.
        let line3 = "1\t10\t.\tA\tC,<*>\t.\t.\tAC=1,0\tGT\t0/1";
        assert_eq!(trim_alt_alleles(line3, &h, false, true), "1\t10\t.\tA\tC\t.\t.\tAC=1\tGT\t0/1");
    }

    #[test]
    fn ac_an_recomputed_after_subsetting() {
        let h = hdr();
        let line = "1\t10\t.\tA\tC\t.\t.\tAC=3;AN=6\tGT\t0/1\t1/1";
        assert_eq!(update_ac_an(line, &h), "1\t10\t.\tA\tC\t.\t.\tAC=3;AN=4\tGT\t0/1\t1/1");
    }

    #[test]
    fn type_filter_uses_shared_classification() {
        let t = TypeFilter::parse(Some("snps"), None).unwrap().unwrap();
        assert!(t.passes("A", "C,AT"));
        assert!(!t.passes("AT", "A"));
        let t = TypeFilter::parse(None, Some("indels")).unwrap().unwrap();
        assert!(!t.passes("A", "C,AT"));
        assert!(t.passes("A", "C"));
        let t = TypeFilter::parse(Some("ref"), None).unwrap().unwrap();
        assert!(t.passes("A", "."));
        assert!(!t.passes("A", "C"));
    }

    #[test]
    fn genotype_filter_negation() {
        let f = GenotypeFilter::parse(Some("^miss")).unwrap().unwrap();
        assert!(f.passes("GT", &["0/1", "1/1"], None));
        assert!(!f.passes("GT", &["0/1", "./."], None));
        let f = GenotypeFilter::parse(Some("het")).unwrap().unwrap();
        assert!(f.passes("GT", &["0/1"], None));
        assert!(!f.passes("GT", &["1/1"], None));
    }

    #[test]
    fn allele_count_modes() {
        let (ac, an) = AcAfFilter::counts("GT", &["0/1", "1/1", "0/2"], None, ".", 3);
        assert_eq!((ac.clone(), an), (vec![2, 3, 1], 6));
        assert_eq!(AcAfFilter::select(&ac, AcMode::Nref), 4);
        assert_eq!(AcAfFilter::select(&ac, AcMode::Alt1), 3);
        assert_eq!(AcAfFilter::select(&ac, AcMode::Minor), 2);
        assert_eq!(AcAfFilter::select(&ac, AcMode::Major), 3);
        assert_eq!(AcAfFilter::select(&ac, AcMode::Nonmajor), 3);
        let (ac, an) = AcAfFilter::counts("", &[], None, "AC=2,1;AN=10", 3);
        assert_eq!((ac, an), (vec![7, 2, 1], 10));
    }
