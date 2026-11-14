use super::*;
use crate::annotate::structs::ani::make_variant_key;

#[test]
fn unknown_contig_cannot_alias_a_dictionary_id() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("in.vcf");
    let ani = dir.path().join("in.ani");
    std::fs::write(&vcf, "##fileformat=VCFv4.2\n##contig=<ID=scaf0>\n##contig=<ID=scaf1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nscaf1\t10\t.\tA\tC\t.\t.\t.\n").unwrap();
    crate::annotate::builder_v2::build_ani_index_auto_v2(&vcf, &ani).unwrap();
    let index = AniIndex::open(&ani).unwrap();
    assert_eq!(index.contig_id("scaf1"), Some(1));
    assert!(index.lookup_exact("scaf1", 10, "A", "C").is_some());
    assert!(index.lookup_exact("chr1", 10, "A", "C").is_none());
    assert!(
        crate::annotate::cpu_v2::extract_line_keys("chr1\t10\t.\tA\tC\t.\t.\t.", &index, false)
            .is_none()
    );
}

#[test]
fn interval_queries_match_the_linear_oracle() {
    let mut entries = vec![IntervalEntry {
        start: 0,
        end: 4000,
        entry_idx: 0,
    }];
    for i in 1..1000usize {
        entries.push(IntervalEntry {
            start: i as u32 * 2,
            end: i as u32 * 2 + (i % 7) as u32,
            entry_idx: i,
        });
    }
    let index = IntervalIndex::new(entries.clone());
    for pos in 0..5000 {
        let expected = entries
            .iter()
            .rev()
            .find(|e| e.start <= pos && pos <= e.end)
            .map(|e| e.entry_idx);
        assert_eq!(index.find(pos), expected, "{pos}");
    }
    assert_eq!(IntervalIndex::new(Vec::new()).find(0), None);
}

/// Verify the `cached_entry_keys` slice matches what `make_variant_key`
/// would compute for every entry.
#[test]
fn cached_entry_keys_match_recomputed() {
    use crate::annotate::builder_v2::build_ani_index_auto_v2;
    let tmp_dir = std::env::temp_dir();
    let vcf = tmp_dir.join(format!("kira_entry_keys_test_{}.vcf", std::process::id()));
    let ani = tmp_dir.join(format!("kira_entry_keys_test_{}.ani", std::process::id()));
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

#[test]
fn sparse_position_routing_matches_postings_and_range_oracles() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("sparse.vcf");
    let ani = dir.path().join("sparse.ani");
    let mut source = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=scaf0>\n##contig=<ID=scaf1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    let mut positions = Vec::new();
    for i in 0..4096u32 {
        let pos = i * 1024 + 511;
        positions.push(pos);
        source.push_str(&format!("scaf0\t{pos}\t.\tA\tC\t.\t.\t.\n"));
        if i % 257 == 0 {
            source.push_str(&format!("scaf0\t{pos}\t.\tA\tG\t.\t.\t.\n"));
        }
    }
    source.push_str("scaf1\t511\t.\tA\tT\t.\t.\t.\n");
    std::fs::write(&vcf, source).unwrap();
    crate::annotate::builder_v2::build_ani_index_auto_v2(&vcf, &ani).unwrap();
    let index = AniIndex::open(&ani).unwrap();
    let chr = index.contig_id("scaf0").unwrap();
    let contig_idx = index.pos_contig_lut.get(chr).unwrap();
    assert!(index.pos_block_lut.contains_key(&contig_idx));
    assert_eq!(mem::size_of::<AniPosBlock>(), 80);
    assert_eq!(mem::size_of::<AniPosIndexHeader>(), 36);
    for (rank, &pos) in positions.iter().enumerate() {
        let postings = index.lookup_pos_index(chr, pos).unwrap();
        assert_eq!(postings.len(), if rank % 257 == 0 { 2 } else { 1 });
        assert!(postings.iter().all(|&entry| {
            let entry = &index.entries[entry as usize];
            entry.chr_id == chr && entry.pos == pos
        }));
        assert!(index.has_site(chr, pos));
        assert!(!index.has_site(chr, pos - 1));
        assert!(!index.has_site(chr, pos + 1));
    }
    for (lo, hi) in [
        (0, 0),
        (0, 512),
        (500, 1540),
        (511, 511),
        (512, 1535),
        (u32::MAX - 1, u32::MAX),
    ] {
        let expected: Vec<_> = positions
            .iter()
            .copied()
            .filter(|&p| lo <= p && p <= hi)
            .collect();
        let actual: Vec<_> = index
            .pos_index_range(chr, lo, hi)
            .map(|(p, list)| {
                assert_eq!(Some(list), index.lookup_pos_index(chr, p));
                p
            })
            .collect();
        assert_eq!(actual, expected, "{lo}..={hi}");
    }
    assert!(index.pos_index_range(chr, 1024, 512).next().is_none());
    assert!(
        index
            .pos_index_range(u32::MAX, 0, u32::MAX)
            .next()
            .is_none()
    );
    assert_eq!(
        index
            .lookup_pos_index(index.contig_id("scaf1").unwrap(), 511)
            .unwrap()
            .len(),
        1
    );
}

#[test]
fn canonical_range_keeps_first_hit_end_and_interval_matching() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("matching.vcf");
    let ani = dir.path().join("matching.ani");
    let padding = "A".repeat(32);
    let source = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=scaf0>\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nscaf0\t480\tfirst\t{padding}C\t{padding}T\t.\t.\t.\nscaf0\t511\tsecond\tAC\tAT\t.\t.\t.\nscaf0\t700\tshort\tA\t<DEL>\t.\t.\tEND=710\nscaf0\t700\tlong\tA\t<DEL>\t.\t.\tEND=720\nscaf0\t1000\tinterval\t.\t2000\t.\t.\t.\n"
    );
    std::fs::write(&vcf, source).unwrap();
    crate::annotate::builder_v2::build_ani_index_auto_v2(&vcf, &ani).unwrap();
    let index = AniIndex::open(&ani).unwrap();
    let chr = index.contig_id("scaf0").unwrap();
    let meta = HashMap::new();
    let lookup = |pos, rf: &str, alt: &str, end| {
        index.lookup_exact_by_chr_id_pos_index_matching(
            chr,
            pos,
            rf,
            crate::util::fast_hash64(rf.as_bytes()),
            alt,
            &meta,
            true,
            false,
            None,
            end,
        )
    };
    let canonical = lookup(512, "C", "T", None).unwrap();
    assert_eq!(canonical.id.as_deref(), Some("first"));
    assert_eq!(canonical.db_ref, "C");
    assert_eq!(canonical.alt, "T");
    assert!(lookup(512, "C", "T", Some(512)).is_none());
    assert_eq!(
        lookup(700, "A", "<DEL>", Some(720)).unwrap().id.as_deref(),
        Some("long")
    );
    assert!(lookup(700, "A", "<DEL>", Some(730)).is_none());
    assert_eq!(
        lookup(1500, "C", "G", None).unwrap().id.as_deref(),
        Some("interval")
    );
    let populated: Vec<_> = index
        .pos_index_range(chr, 480, 512)
        .map(|(p, _)| p)
        .collect();
    assert_eq!(populated, [480, 511]);
}

#[test]
fn block_directory_skips_dense_and_small_contigs() {
    let block = |base_pos| AniPosBlock {
        base_pos,
        _pad: 0,
        masks: [0; 8],
        offsets_start: 0,
        _pad2: 0,
    };
    let mut blocks: Vec<_> = (0..2048u32).map(|i| block(i * 512)).collect();
    blocks.extend((0..2047u32).map(|i| block(i * 1024)));
    let index = AniPosIndex {
        contigs: vec![
            AniPosContig {
                chr_id: 0,
                min_pos: 1,
                max_pos: 2047 * 512,
                block_start: 0,
                block_count: 2048,
            },
            AniPosContig {
                chr_id: 1,
                min_pos: 1,
                max_pos: 2046 * 1024,
                block_start: 2048,
                block_count: 2047,
            },
        ],
        blocks,
        pos_offsets: Vec::new(),
        pos_counts: Vec::new(),
        entry_indices: Vec::new(),
    };
    assert!(build_pos_block_lut(&index).is_empty());
}

#[test]
fn populated_range_preserves_word_and_coordinate_boundaries() {
    let positions = [0, 63, 64, 511, 512, 575, 576, 1023, u32::MAX];
    let mut blocks: Vec<AniPosBlock> = Vec::new();
    for (rank, &pos) in positions.iter().enumerate() {
        let base = pos / 512 * 512;
        if blocks.last().is_none_or(|b| b.base_pos != base) {
            blocks.push(AniPosBlock {
                base_pos: base,
                _pad: 0,
                masks: [0; 8],
                offsets_start: rank as u32,
                _pad2: 0,
            });
        }
        let bit = (pos - base) as usize;
        blocks.last_mut().unwrap().masks[bit / 64] |= 1u64 << (bit % 64);
    }
    let index = AniPosIndex {
        contigs: Vec::new(),
        blocks,
        pos_offsets: (0..positions.len() as u32).collect(),
        pos_counts: vec![1; positions.len()],
        entry_indices: (0..positions.len() as u32).collect(),
    };
    for (lo, hi) in [
        (0, u32::MAX),
        (63, 64),
        (511, 512),
        (575, 576),
        (65, 510),
        (u32::MAX, u32::MAX),
    ] {
        let range = AniPosRange {
            pos_index: Some(&index),
            blocks: &index.blocks,
            block_idx: index
                .blocks
                .partition_point(|b| b.base_pos < lo / 512 * 512),
            word_idx: 0,
            word_base: 0,
            rank_base: 0,
            original_mask: 0,
            remaining: 0,
            lo,
            hi,
        };
        let actual: Vec<_> = range.map(|(p, list)| (p, list[0])).collect();
        let expected: Vec<_> = positions
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, p)| lo <= p && p <= hi)
            .map(|(rank, p)| (p, rank as u32))
            .collect();
        assert_eq!(actual, expected, "{lo}..={hi}");
    }
}

#[test]
fn malformed_block_geometry_does_not_overflow_point_or_range_routing() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("geometry.vcf");
    let ani = dir.path().join("geometry.ani");
    std::fs::write(&vcf, "##fileformat=VCFv4.2\n##contig=<ID=scaf0>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nscaf0\t512\t.\tA\tC\t.\t.\t.\nscaf0\t1024\t.\tA\tG\t.\t.\t.\n").unwrap();
    crate::annotate::builder_v2::build_ani_index_auto_v2(&vcf, &ani).unwrap();
    let mut index = AniIndex::open(&ani).unwrap();
    let chr = index.contig_id("scaf0").unwrap();
    let contig_idx = index.pos_contig_lut.get(chr).unwrap();
    assert!(index.pos_block_lut.is_empty());
    index.pos_index.as_mut().unwrap().contigs[contig_idx].min_pos = 0;
    assert!(index.lookup_pos_index(chr, 1).is_none());
    {
        let positions = index.pos_index.as_mut().unwrap();
        let contig = &positions.contigs[contig_idx];
        let start = contig.block_start as usize;
        let end = start + contig.block_count as usize;
        positions.blocks[start..end].reverse();
    }
    assert!(index.lookup_pos_index(chr, 1).is_none());
    {
        let positions = index.pos_index.as_mut().unwrap();
        let contig = &mut positions.contigs[contig_idx];
        contig.block_start = 0;
        contig.block_count = 1;
        contig.max_pos = u32::MAX;
        positions.blocks = vec![AniPosBlock {
            base_pos: u32::MAX,
            _pad: 0,
            masks: [0; 8],
            offsets_start: 0,
            _pad2: 0,
        }];
    }
    assert!(
        index
            .pos_index_range(chr, u32::MAX, u32::MAX)
            .next()
            .is_none()
    );
    assert!(index.pos_index_range(chr, 0, u32::MAX).next().is_none());
}

#[test]
fn block_directory_rejects_misaligned_and_unordered_geometry() {
    let blocks: Vec<_> = (0..4096u32)
        .map(|rank| AniPosBlock {
            base_pos: rank * 1024,
            _pad: 0,
            masks: [0; 8],
            offsets_start: 0,
            _pad2: 0,
        })
        .collect();
    let mut index = AniPosIndex {
        contigs: vec![AniPosContig {
            chr_id: 0,
            min_pos: 0,
            max_pos: 4095 * 1024,
            block_start: 0,
            block_count: 4096,
        }],
        blocks,
        pos_offsets: Vec::new(),
        pos_counts: Vec::new(),
        entry_indices: Vec::new(),
    };
    index.blocks[100].base_pos += 1;
    assert!(build_pos_block_lut(&index).is_empty());
    index.blocks[100].base_pos -= 1;
    index.blocks.swap(100, 101);
    assert!(build_pos_block_lut(&index).is_empty());
}
