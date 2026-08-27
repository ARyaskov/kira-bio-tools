    use super::*;

    fn mk_site(chrom: &str, pos: u32, gt: GtClass, af: f64) -> RohSite {
        RohSite { chrom: chrom.into(), pos, genetic_pos: None, gt, af }
    }

    #[test]
    fn pure_homozygous_run_called_az() {
        let mut sites = Vec::new();
        for i in 0..100 {
            let gt = if i % 5 == 0 { GtClass::HomAlt } else { GtClass::HomRef };
            sites.push(mk_site("1", 1000 + i * 1000, gt, 0.3));
        }
        let opts = RohOpts::default();
        let path = viterbi(&sites, &opts);
        let az_count = path.iter().filter(|s| **s == State::AZ).count();
        assert!(az_count >= 50, "expected majority AZ in homozygous block, got {}", az_count);
    }

    #[test]
    fn het_heavy_region_called_hw() {
        let mut sites = Vec::new();
        for i in 0..100 {
            let gt = if i % 2 == 0 { GtClass::Het } else { GtClass::HomRef };
            sites.push(mk_site("1", 1000 + i * 1000, gt, 0.3));
        }
        let opts = RohOpts::default();
        let path = viterbi(&sites, &opts);
        let hw_count = path.iter().filter(|s| **s == State::HW).count();
        assert!(hw_count >= 80, "expected majority HW with hets, got {}", hw_count);
    }

    #[test]
    fn segments_collapse_correctly() {
        let sites = vec![
            mk_site("1", 100, GtClass::HomRef, 0.3),
            mk_site("1", 200, GtClass::HomRef, 0.3),
            mk_site("1", 300, GtClass::Het, 0.3),
            mk_site("1", 400, GtClass::Het, 0.3),
        ];
        let path = vec![State::AZ, State::AZ, State::HW, State::HW];
        let segs = segments(&sites, &path);
        assert_eq!(segs.len(), 2);
        assert_eq!(segs[0].state, State::AZ);
        assert_eq!(segs[0].n_sites, 2);
        assert_eq!(segs[1].state, State::HW);
        assert_eq!(segs[1].n_sites, 2);
    }

    #[test]
    fn estimate_af_balanced_window() {
        let gts = vec![GtClass::HomRef, GtClass::Het, GtClass::HomAlt, GtClass::Het, GtClass::HomRef];
        let af = estimate_af(&gts, 5);
        for f in af { assert!(f >= 0.01 && f <= 0.99); }
    }

    #[test]
    fn baum_welch_converges_on_pure_az() {
        let mut sites = Vec::new();
        for i in 0..200 {
            sites.push(mk_site("1", 1000 + i * 1000, GtClass::HomRef, 0.3));
        }
        let mut opts = RohOpts::default();
        let init_hw_az = opts.hw_to_az;
        baum_welch_train(&sites, &mut opts, 3);
        assert!(opts.hw_to_az != init_hw_az || opts.az_to_hw != RohOpts::default().az_to_hw,
                "expected param update after training");
    }
