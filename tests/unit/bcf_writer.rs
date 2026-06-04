    use super::*;
    use crate::bcf::reader::BcfReader;

    #[test]
    fn write_then_read_back() {
        let h = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "##contig=<ID=1>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1".to_string(),
        ];
        let mut tmp = std::env::temp_dir();
        tmp.push(format!("kira_bt_bcf_test_{}.bcf", std::process::id()));
        let _ = std::fs::remove_file(&tmp);
        {
            let mut w = BcfWriter::create(&tmp, false, 0, &h).unwrap();
            w.write_vcf_line("1\t100\t.\tA\tT\t60\tPASS\tDP=42\tGT\t0/1").unwrap();
            w.write_vcf_line("1\t200\trs1\tC\tG\t.\t.\t.\tGT\t1|1").unwrap();
            w.finish().unwrap();
        }
        let mut r = BcfReader::open(&tmp).unwrap();
        let l1 = r.read_record_line().unwrap().unwrap();
        let l2 = r.read_record_line().unwrap().unwrap();
        let l3 = r.read_record_line().unwrap();
        assert!(l1.starts_with("1\t100\t."));
        assert!(l1.contains("DP=42"));
        assert!(l2.contains("rs1"));
        assert!(l3.is_none());
        let _ = std::fs::remove_file(&tmp);
    }
