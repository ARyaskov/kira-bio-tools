    use super::build_ani_index_from_tab;
    use crate::annotate::structs::ani::AniIndex;
    use std::fs;

    #[test]
    fn test_tab_infers_numbers_for_split_info() {
        let dir = std::env::temp_dir();
        let tab_path = dir.join("kira_tab_infer_test.tab");
        let ani_path = dir.join("kira_tab_infer_test.ani");

        let tab = "1\t1\tC\tA,T\t0,1.1\t1.1,0,2.2\n";
        fs::write(&tab_path, tab).unwrap();

        build_ani_index_from_tab(&tab_path, &ani_path, Some("CHROM,POS,REF,ALT,FA,FR")).unwrap();

        let ani = AniIndex::open(&ani_path).unwrap();
        let bundle = ani.lookup_exact("1", 1, "C", "T").unwrap();

        let fa = bundle
            .info
            .iter()
            .find(|f| f.key == "FA")
            .map(|f| f.values.clone())
            .unwrap();
        let fr = bundle
            .info
            .iter()
            .find(|f| f.key == "FR")
            .map(|f| f.values.clone())
            .unwrap();

        assert_eq!(fa, vec!["1.1"]);
        assert_eq!(fr, vec!["1.1", "2.2"]);

        let _ = fs::remove_file(&tab_path);
        let _ = fs::remove_file(&ani_path);
    }
