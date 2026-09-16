use super::*;

fn write_gff(dir: &std::path::Path, text: &str) -> std::path::PathBuf {
    let p = dir.join("ann.gff");
    std::fs::write(&p, text).unwrap();
    p
}

#[test]
fn biotype_matching_mirrors_strncmp() {
    assert_eq!(
        parse_biotype(Some(b"protein_coding;extra")),
        GF_PROTEIN_CODING as i32
    );
    assert_eq!(parse_biotype(Some(b"mRNA")), GF_PROTEIN_CODING as i32);
    assert_eq!(parse_biotype(Some(b"Mt_rRNA")), GF_MT_TRNA as i32);
    // strncmp(line,"lncRNA",7) only matches at the end of the line.
    assert_eq!(parse_biotype(Some(b"lncRNA")), GF_LNCRNA as i32);
    assert_eq!(parse_biotype(Some(b"lncRNA;Name=x")), 0);
    assert_eq!(
        parse_biotype(Some(b"nonsense_mediated_decay")),
        GF_NMD as i32
    );
    assert_eq!(parse_biotype(Some(b"unknown_thing")), 0);
    assert_eq!(parse_biotype(None), -1);
    assert_eq!(gf_type2gff_string(GF_PROTEIN_CODING), "protein_coding");
    assert_eq!(gf_type2gff_string(GF_LINCRNA), "lincRNA");
    assert_eq!(gf_type2gff_string(GF_NMD), "NMD");
    assert!(gf_is_coding(GF_IG_V));
    assert!(!gf_is_coding(GF_PSEUDOGENE));
}

#[test]
fn regidx_returns_overlaps_sorted_longest_first() {
    let mut idx = RegIdx::<u32>::default();
    idx.push("1", 10, 20, 1);
    idx.push("1", 5, 30, 2);
    idx.push("1", 10, 25, 3);
    idx.push("1", 40, 50, 4);
    idx.push("2", 0, 100, 5);
    idx.finalize();
    assert_eq!(idx.overlap("1", 12, 12), vec![2, 3, 1]);
    assert_eq!(idx.overlap("1", 0, 4), Vec::<u32>::new());
    assert_eq!(idx.overlap("1", 30, 40), vec![2, 4]);
    assert_eq!(idx.overlap("3", 0, 1), Vec::<u32>::new());
    assert_eq!(idx.nregs(), 5);
}

#[test]
fn parses_ensembl_style_gff_and_fixes_phase() {
    let dir = tempfile::tempdir().unwrap();
    let p = write_gff(
        dir.path(),
        "##gff-version 3\n\
                 chr1\t.\tgene\t90\t200\t.\t+\t.\tID=gene:G1;Name=XYZ;biotype=protein_coding\n\
                 chr1\t.\ttranscript\t90\t200\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n\
                 chr1\t.\texon\t90\t110\t.\t+\t.\tParent=transcript:T1\n\
                 chr1\t.\tfive_prime_UTR\t90\t98\t.\t+\t.\tParent=transcript:T1\n\
                 chr1\t.\tCDS\t99\t110\t.\t+\t1\tParent=transcript:T1\n\
                 chr1\t.\tCDS\t120\t130\t.\t+\t1\tParent=transcript:T1\n\
                 chr1\t.\tthree_prime_UTR\t131\t200\t.\t+\t.\tParent=transcript:T1\n\
                 chr1\t.\tbiological_region\t1\t5\t.\t+\t.\tID=skip\n\
                 chr1\t.\ttranscript\t300\t400\t.\t-\t.\tID=transcript:T2;Parent=gene:G2;biotype=TCE\n",
    );
    let gff = parse(
        &p,
        &GffOpts {
            strip_chr_names: true,
            force: false,
            verbosity: 0,
        },
    )
    .unwrap();
    assert_eq!(gff.seqs, vec!["1"]);
    assert_eq!(gff.tscripts.len(), 1);
    let tr = &gff.tscripts[0];
    assert_eq!(gff.tscript_name(tr.id), "T1");
    assert_eq!(gff.gene_name(tr.gene), Some("XYZ"));
    assert_eq!(tr.strand, STRAND_FWD);
    assert_eq!(tr.trim, TRIM_5PRIME | TRIM_3PRIME);
    // Phase 1 drops one base from the first CDS: 99..110 -> 100..110 (0-based 99, len 11).
    let c0 = &gff.cds[tr.cds[0]];
    assert_eq!((c0.beg, c0.len, c0.pos, c0.phase), (99, 11, 0, 0));
    // 11 + 11 = 22 bases, one trimmed from the 3' end.
    let c1 = &gff.cds[tr.cds[1]];
    assert_eq!((c1.beg, c1.len, c1.pos), (119, 10, 11));
    assert_eq!(gff.idx_cds.overlap("1", 105, 105), vec![tr.cds[0]]);
    assert_eq!(gff.idx_utr.overlap("1", 95, 95).len(), 1);
    // Exons are indexed with the intronic splice-region margin.
    assert_eq!(gff.idx_exon.overlap("1", 115, 115).len(), 1);
    assert_eq!(gff.idx_tscript.overlap("1", 150, 150), vec![0]);
}

#[test]
fn wrong_phase_is_an_error_unless_forced() {
    let dir = tempfile::tempdir().unwrap();
    let text = "1\t.\tgene\t1\t100\t.\t+\t.\tID=gene:G1;biotype=protein_coding\n\
                        1\t.\ttranscript\t1\t100\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n\
                        1\t.\tCDS\t1\t10\t.\t+\t0\tParent=transcript:T1\n\
                        1\t.\tCDS\t20\t30\t.\t+\t0\tParent=transcript:T1\n";
    let p = write_gff(dir.path(), text);
    let err = match parse(&p, &GffOpts::default()) {
        Ok(_) => panic!("expected a phase error"),
        Err(e) => e.to_string(),
    };
    assert!(err.contains("GFF3 assumption failed"), "{err}");
    let gff = parse(
        &p,
        &GffOpts {
            force: true,
            ..Default::default()
        },
    )
    .unwrap();
    // The transcript is kept but its CDS are not indexed.
    assert_eq!(gff.idx_tscript.nregs(), 1);
    assert_eq!(gff.idx_cds.nregs(), 0);
}
