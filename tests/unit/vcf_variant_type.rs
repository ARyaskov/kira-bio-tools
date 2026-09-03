    use super::*;

    #[test]
    fn classifies_like_htslib() {
        assert_eq!(allele_type("A", "C").ty, VT_SNP);
        assert_eq!(allele_type("A", "A").ty, VT_REF);
        assert_eq!(allele_type("A", ".").ty, VT_REF);
        assert_eq!(allele_type("A", "*").ty, VT_OVERLAP);
        assert_eq!(allele_type("A", "AT"), AlleleType { ty: VT_INDEL, n: 1 });
        assert_eq!(allele_type("AT", "A"), AlleleType { ty: VT_INDEL, n: -1 });
        assert_eq!(allele_type("ATT", "AGT").ty, VT_SNP);
        assert_eq!(allele_type("AC", "TG").ty, VT_MNP);
        assert_eq!(allele_type("ACT", "TCGA").ty, VT_OTHER);
        assert_eq!(allele_type("A", "<DEL>").ty, VT_OTHER);
        assert_eq!(allele_type("A", "<*>").ty, VT_REF);
        assert_eq!(allele_type("A", "<NON_REF>").ty, VT_REF);
        assert_eq!(allele_type("G", "G]17:198982]").ty, VT_BND);
        assert_eq!(allele_type("T", "]13:123456]T").ty, VT_BND);
        assert_eq!(allele_type("a", "T").ty, VT_SNP);
        assert_eq!(allele_type("a", "A").ty, VT_REF);
    }

    #[test]
    fn record_union() {
        assert_eq!(record_type("A", "C,AT"), VT_SNP | VT_INDEL);
        assert_eq!(record_type("A", "."), VT_REF);
        assert!(has_snp("A", "C,AT"));
        assert!(has_indel("A", "C,AT"));
        assert!(!is_pure_snp("A", "C,AT"));
        assert!(is_pure_snp("A", "C,G"));
        assert_eq!(type_name(record_type("AC", "A")), "indel");
        assert_eq!(parse_type_mask("snps,indels"), Some(VT_SNP | VT_INDEL));
        assert_eq!(parse_type_mask("bogus"), None);
    }
