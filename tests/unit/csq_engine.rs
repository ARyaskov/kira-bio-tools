    use super::*;

    #[test]
    fn codon_tables() {
        assert_eq!(dna2aa(b"ATG"), b'M');
        assert_eq!(dna2aa(b"TAA"), b'*');
        assert_eq!(dna2aa(b"tgg"), b'W');
        // Reverse strand: CAT on the reference is ATG on the coding strand.
        assert_eq!(cdna2aa(b"CAT"), b'M');
        assert_eq!(cdna2aa(b"TTA"), b'*');
        assert_eq!(dna2aa(b"NNN"), b'X');
    }

    #[test]
    fn genotype_parsing() {
        let s = vec![
            "0|1".to_string(),
            "1/2:3".to_string(),
            ".".to_string(),
            "1".to_string(),
        ];
        let (g, n) = parse_genotypes(Some("GT:DP"), &s);
        assert_eq!(n, 2);
        assert_eq!(
            g[0],
            vec![
                Gt {
                    allele: 0,
                    phased: false
                },
                Gt {
                    allele: 1,
                    phased: true
                }
            ]
        );
        assert_eq!(
            g[1],
            vec![
                Gt {
                    allele: 1,
                    phased: false
                },
                Gt {
                    allele: 2,
                    phased: false
                }
            ]
        );
        assert_eq!(
            g[2],
            vec![Gt {
                allele: -1,
                phased: false
            }]
        );
        assert_eq!(g[3].len(), 1);
        assert_eq!(parse_genotypes(Some("DP"), &s).1, 0);
        assert_eq!(parse_genotypes(None, &s).1, 0);
    }

    #[test]
    fn bitmask_layout() {
        assert_eq!(ncsq2_to_nfmt(30), 1);
        assert_eq!(ncsq2_to_nfmt(31), 2);
        assert_eq!(ncsq2_to_nfmt(64), 3);
        assert_eq!(icsq2_to_bit(31), (1, 1));
    }

    #[test]
    fn translate_forward_and_reverse() {
        // Spliced reference: 10 N pad + ATGGCCTAA + 10 N pad.
        let mut sref = vec![b'N'; N_REF_PAD];
        sref.extend_from_slice(b"ATGGCCTAA");
        sref.extend(vec![b'N'; N_REF_PAD]);
        // Whole transcript, forward.
        let aa = cds_translate(&sref, b"ATGGCCTAA", 9, 0, 0, 9, STRAND_FWD, false).unwrap();
        assert_eq!(aa, b"MA*");
        // One base in the middle of the second codon: padding from the reference.
        let aa = cds_translate(&sref, b"T", 9, 4, 4, 5, STRAND_FWD, false).unwrap();
        assert_eq!(aa, b"V");
        // Frameshift fill runs to the transcript end.
        let aa = cds_translate(&sref, b"GT", 10, 3, 3, 4, STRAND_FWD, true).unwrap();
        assert_eq!(aa, b"VL");
        // Reverse strand: TTAGGCCAT is ATGGCCTAA read backwards.
        let mut sref = vec![b'N'; N_REF_PAD];
        sref.extend_from_slice(b"TTAGGCCAT");
        sref.extend(vec![b'N'; N_REF_PAD]);
        let aa = cds_translate(&sref, b"TTAGGCCAT", 9, 0, 0, 9, gff::STRAND_REV, false).unwrap();
        assert_eq!(aa, b"MA*");
        let aa = cds_translate(&sref, b"A", 9, 4, 4, 5, gff::STRAND_REV, false).unwrap();
        assert_eq!(aa, b"V");
        assert_eq!(
            cds_translate(&sref, b"", 9, 0, 0, 0, STRAND_FWD, false).unwrap(),
            b"?"
        );
    }

    #[test]
    fn consequence_string_rendering() {
        let mut gff = Gff::default();
        gff.tscript_ids.push("ENST1".into());
        gff.genes.push(gff::Gene {
            name: Some("XYZ".into()),
            ..Default::default()
        });
        gff.tscripts.push(gff::Tscript {
            id: 0,
            beg: 0,
            end: 10,
            strand: STRAND_FWD,
            used: true,
            cds: Vec::new(),
            trim: 0,
            ty: gff::GF_PROTEIN_CODING,
            gene: 0,
        });
        let mut v = Vcsq::new(&gff, 0, CSQ_MISSENSE_VARIANT, 1);
        v.vstr = Some("|2V>2I|103G>A".into());
        let mut s = String::new();
        kput_vcsq(&gff, &mut v, &mut s);
        assert_eq!(s, "missense|XYZ|ENST1|protein_coding|+|2V>2I|103G>A");

        let mut v = Vcsq::new(&gff, 0, CSQ_SPLICE_REGION | CSQ_UTR5, 1);
        let mut s = String::new();
        kput_vcsq(&gff, &mut v, &mut s);
        assert_eq!(s, "splice_region&5_prime_utr|XYZ|ENST1|protein_coding");

        let mut v = Vcsq::new(&gff, 0, CSQ_INTRON, 1);
        let mut s = String::new();
        kput_vcsq(&gff, &mut v, &mut s);
        assert_eq!(s, "intron|XYZ||protein_coding");

        let mut v = Vcsq::new(&gff, 0, CSQ_PRINTED_UPSTREAM | CSQ_SPLICE_REGION, 1);
        v.ref_ = Some((7, 106));
        let mut s = String::new();
        kput_vcsq(&gff, &mut v, &mut s);
        assert_eq!(s, "@107");

        // Stop/start flags drop the accompanying missense; upstream stop adds `*`.
        let mut v = Vcsq::new(
            &gff,
            0,
            CSQ_STOP_GAINED | CSQ_MISSENSE_VARIANT | CSQ_UPSTREAM_STOP,
            1,
        );
        v.vstr = Some("|x".into());
        let mut s = String::new();
        kput_vcsq(&gff, &mut v, &mut s);
        assert_eq!(s, "*stop_gained|XYZ|ENST1|protein_coding|+|x");
    }

    #[test]
    fn brief_protein_predictions() {
        let mut s = String::new();
        kprint_aa_prediction(0, 1, b"MAKV*", &mut s);
        assert_eq!(s, "MAKV*");
        let mut s = String::new();
        kprint_aa_prediction(2, 12, b"YVRTLKQ*", &mut s);
        assert_eq!(s, "YV..19");
    }

    #[test]
    fn info_and_format_updates() {
        assert_eq!(caller::update_info(".", "BCSQ", "x"), "BCSQ=x");
        assert_eq!(
            caller::update_info("EXP=a;DP=3", "BCSQ", "x"),
            "EXP=a;DP=3;BCSQ=x"
        );
        assert_eq!(
            caller::update_info("BCSQ=old;DP=3", "BCSQ", "x"),
            "BCSQ=x;DP=3"
        );
        let mut rec = crate::vcf::VcfRecord {
            chrom: "1".into(),
            pos: 5,
            id: ".".into(),
            ref_allele: "A".into(),
            alt: "T".into(),
            qual: ".".into(),
            filter: ".".into(),
            info: ".".into(),
            format: Some("GT:DP".into()),
            samples: vec!["0|1:3".into(), "1|1".into()],
            chr_id: Default::default(),
            position: 5,
            offset: 0,
        };
        caller::update_format(&mut rec, "BCSQ", Some(vec!["2".into(), "3,0".into()]));
        assert_eq!(rec.format.as_deref(), Some("GT:DP:BCSQ"));
        assert_eq!(rec.samples, vec!["0|1:3:2", "1|1:.:3,0"]);
        caller::update_format(&mut rec, "BCSQ", None);
        assert_eq!(rec.format.as_deref(), Some("GT:DP"));
        assert_eq!(rec.samples, vec!["0|1:3", "1|1:."]);
    }
