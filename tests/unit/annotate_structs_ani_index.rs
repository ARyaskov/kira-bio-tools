    use super::*;
    use crate::annotate::structs::ani::make_variant_key;

    /// Verify the `cached_entry_keys` slice matches what `make_variant_key`
    /// would compute for every entry.
    #[test]
    fn cached_entry_keys_match_recomputed() {
        use crate::annotate::builder_v2::build_ani_index_auto_v2;
        let tmp_dir = std::env::temp_dir();
        let vcf = tmp_dir.join(format!(
            "kira_entry_keys_test_{}.vcf",
            std::process::id()
        ));
        let ani = tmp_dir.join(format!(
            "kira_entry_keys_test_{}.ani",
            std::process::id()
        ));
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##contig=<ID=1>\n\
             ##contig=<ID=2>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t100\trs1\tA\tT\t.\t.\tDP=10\n\
             1\t200\trs2\tC\tG\t.\t.\tDP=20\n\
             1\t300\trs3\tA\tAT\t.\t.\tDP=30\n\
             2\t50\trs4\tG\tA\t.\t.\tDP=40\n",
        )
        .unwrap();

        build_ani_index_auto_v2(&vcf, &ani).expect("ani build");
        let index = AniIndex::open(&ani).expect("ani open");

        let cached = index
            .cached_entry_keys()
            .expect("HAS_ENTRY_KEYS flag should be set by the current builder");
        assert_eq!(cached.len(), index.entries.len());

        let mut checked = 0usize;
        for (i, entry) in index.entries.iter().enumerate() {
            if entry.chr_id == crate::annotate::structs::ani::ANI_SENTINEL_CHR_ID {
                continue;
            }
            let ref_str = index.read_cstring(entry.ref_ofs as usize);
            let alt_str = index.read_cstring(entry.alt_ofs as usize);
            let expected = make_variant_key(
                entry.chr_id,
                entry.pos,
                ref_str.as_ref().as_bytes(),
                alt_str.as_ref().as_bytes(),
            );
            assert_eq!(
                cached[i], expected,
                "entry_keys[{i}] mismatch: cached={:016x} expected={:016x}",
                cached[i], expected
            );
            checked += 1;
        }
        assert_eq!(checked, 4, "expected to verify 4 real entries");

        let _ = std::fs::remove_file(&vcf);
        let _ = std::fs::remove_file(&ani);
    }
