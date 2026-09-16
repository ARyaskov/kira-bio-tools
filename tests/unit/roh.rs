    use super::*;

    fn rec(chrom: &str, pos: u32, alt: &str, info: &str, fmt: &str, samples: &[&str]) -> VcfRecord {
        VcfRecord {
            chrom: chrom.into(),
            pos,
            id: ".".into(),
            ref_allele: "A".into(),
            alt: alt.into(),
            qual: ".".into(),
            filter: ".".into(),
            info: info.into(),
            format: Some(fmt.into()),
            samples: samples.iter().map(|s| s.to_string()).collect(),
            chr_id: Default::default(),
            position: pos,
            offset: 0,
        }
    }

    fn headers() -> Vec<String> {
        vec![
            "##fileformat=VCFv4.2".into(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">".into(),
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"\">".into(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1".into(),
        ]
    }

    fn run(opts: RohOpts, recs: &[VcfRecord]) -> (String, RohStats) {
        let buf: std::sync::Arc<std::sync::Mutex<Vec<u8>>> = Default::default();
        struct W(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl Write for W {
            fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(b);
                Ok(b.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let mut c = RohCaller::new(
            opts,
            &headers(),
            vec!["S1".into()],
            vec![0],
            Box::new(W(buf.clone())),
        )
        .unwrap();
        for r in recs {
            c.process(r, true).unwrap();
        }
        c.finish().unwrap();
        let stats = c.stats.clone();
        let s = String::from_utf8(buf.lock().unwrap().clone()).unwrap();
        (s, stats)
    }

    #[test]
    fn genotype_and_pl_parsing() {
        let r = rec("1", 10, "G", ".", "GT:PL", &["0/1:20,0,30", "1:5", "./.:."]);
        let g = parse_gt_pairs(&r).unwrap();
        assert_eq!(g[0], [Some(0), Some(1)]);
        assert_eq!(g[1], [Some(1), Some(-65)]);
        assert_eq!(g[2], [None, None]);
        let (p, n) = parse_pls(&r).unwrap();
        assert_eq!(n, 3);
        assert_eq!(p[0], vec![20, 0, 30]);
        assert_eq!(p[1], vec![5, -2, -2]);
        assert_eq!(p[2], vec![-1, -2, -2]);
        // Haploid-only sites are not diploid.
        assert!(parse_gt_pairs(&rec("1", 10, "G", ".", "GT", &["1", "0"])).is_none());
        assert!(bcf_is_snp(&["A".into(), "G".into(), "<*>".into()]));
        assert!(!bcf_is_snp(&["AT".into(), "A".into()]));
        assert_eq!(alleles2gt(0, 1), 1);
        assert_eq!(alleles2gt(1, 1), 2);
        assert_eq!(alleles2gt(0, 2), 3);
    }

    #[test]
    fn homozygous_stretch_becomes_an_rg_region() {
        let opts = RohOpts {
            gts_only: Some(30.0),
            dflt_af: 0.4,
            ..Default::default()
        };
        let mut recs = Vec::new();
        for i in 0..2000u32 {
            let gt = if i < 500 || i >= 1500 { "0/1" } else { "1/1" };
            recs.push(rec("1", 1000 + i * 1000, "G", ".", "GT", &[gt]));
        }
        let (out, stats) = run(opts, &recs);
        assert_eq!(stats.ntot, 2000);
        let rg: Vec<&str> = out.lines().filter(|l| l.starts_with("RG\t")).collect();
        assert_eq!(rg.len(), 1, "{out}");
        let f: Vec<&str> = rg[0].split('\t').collect();
        assert_eq!(f[1], "S1");
        assert_eq!(f[2], "1");
        let beg: u32 = f[3].parse().unwrap();
        let end: u32 = f[4].parse().unwrap();
        assert!(
            beg >= 1000 + 480 * 1000 && beg <= 1000 + 520 * 1000,
            "{rg:?}"
        );
        assert!(
            end >= 1000 + 1480 * 1000 && end <= 1000 + 1520 * 1000,
            "{rg:?}"
        );
        let st: Vec<&str> = out.lines().filter(|l| l.starts_with("ST\t")).collect();
        assert_eq!(st.len(), 2000);
        // bcftools quirk: the first site reports the initial state probability.
        assert!(st[0].ends_with("\t0\t3.0"), "{}", st[0]);
    }

    #[test]
    fn site_filters_are_counted() {
        let opts = RohOpts {
            gts_only: Some(30.0),
            dflt_af: 0.4,
            ..Default::default()
        };
        let recs = vec![
            rec("1", 10, "G", ".", "GT", &["0/1"]),
            rec("1", 10, "G", ".", "GT", &["0/1"]),
            rec("1", 20, "G,T", ".", "GT", &["0/1"]),
            rec("1", 30, ".", ".", "GT", &["0/0"]),
            rec("1", 40, "<*>", ".", "GT", &["0/0"]),
            rec("1", 50, "G", ".", "GT", &["./."]),
            rec("1", 60, "G", ".", "GT", &["1/1"]),
        ];
        let (_, stats) = run(opts, &recs);
        assert_eq!(
            (stats.ntot, stats.ndup, stats.nmultiallelic, stats.nno_alt),
            (7, 1, 1, 2)
        );
        // Without any AF source the sites are dropped.
        let (_, stats) = run(
            RohOpts {
                gts_only: Some(30.0),
                ..Default::default()
            },
            &recs[..1],
        );
        assert_eq!(stats.nno_af, 1);
    }

    #[test]
    fn genetic_map_rate() {
        let mut m = GenMap {
            map: vec![(0, 0.0), (1000, 0.01), (3000, 0.05)],
            igenmap: 0,
        };
        assert!((m.rate(0, 1000) - 0.01).abs() < 1e-12);
        assert!((m.rate(1000, 2000) - 0.02).abs() < 1e-12);
        assert_eq!(m.rate(5000, 6000), 0.0);
        let mut t = [0.5, 0.1, 0.2, 0.5];
        RohTransition::RecRate(1e-3).apply(100, 200, &mut t);
        assert!((t[1] - 0.01).abs() < 1e-12 && (t[2] - 0.02).abs() < 1e-12);
        assert!((t[3] - 0.99).abs() < 1e-12 && (t[0] - 0.98).abs() < 1e-12);
    }
