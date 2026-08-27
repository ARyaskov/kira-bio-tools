    use super::*;
    use std::path::PathBuf;

    #[test]
    fn ktile_path_appends_extension() {
        let p = ktile_path_for(&PathBuf::from("data/9.vcf.gz"));
        assert!(p.to_string_lossy().ends_with("9.vcf.gz.ktile"));
    }

    #[test]
    fn freshness_returns_io_error_on_missing_source() {
        let tmp_dir = std::env::temp_dir();
        let ktile = tmp_dir.join(format!(
            "kira_freshness_test_{}.ktile",
            std::process::id()
        ));
        let vcf = tmp_dir.join(format!(
            "kira_freshness_test_{}.vcf",
            std::process::id()
        ));
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tT\t.\t.\t.\n",
        )
        .unwrap();
        super::super::writer::write_ktile_from_vcf(&vcf, &ktile).unwrap();
        std::fs::remove_file(&vcf).unwrap();
        let res = check_ktile_freshness(&ktile, &vcf);
        assert!(res.is_err());
        let _ = std::fs::remove_file(&ktile);
    }
