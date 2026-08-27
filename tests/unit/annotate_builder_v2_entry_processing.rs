    use super::*;
    use std::sync::atomic::Ordering;

    fn entry(id_ofs: u32) -> AniEntry {
        AniEntry {
            chr_id: 0,
            pos: 10,
            ref_ofs: 1,
            alt_ofs: 2,
            id_ofs,
            qual_ofs: 0,
            filter_ofs: 0,
            info_ofs: 0,
            info_len: 0,
            format_ofs: u32::MAX,
            samples_ofs: u32::MAX,
        }
    }

    #[test]
    fn duplicate_variant_keeps_first_entry() {
        let key = make_position_key(0, 10, "A", "C");
        let mut map = FxHashMap::default();
        let mut order = 0;
        let duplicates = AtomicUsize::new(0);
        let collisions = AtomicUsize::new(0);
        insert_or_update_entry(
            key,
            entry(11),
            &mut map,
            &mut order,
            &duplicates,
            &collisions,
            false,
            "1",
            10,
            "A",
            "C",
        );
        insert_or_update_entry(
            key,
            entry(22),
            &mut map,
            &mut order,
            &duplicates,
            &collisions,
            false,
            "1",
            10,
            "A",
            "C",
        );
        assert_eq!(map.get(&key).unwrap().entry.id_ofs, 11);
        assert_eq!(duplicates.load(Ordering::Relaxed), 1);
        assert_eq!(order, 1);
    }
