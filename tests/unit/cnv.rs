    use super::*;

    fn pt(chr: &str, pos: u32, baf: f64, lrr: f64) -> SamplePoint {
        SamplePoint {
            chr: chr.to_string(),
            pos,
            baf,
            lrr,
            cn: 2,
            probs: [0.0; N_STATES],
            is_het: true,
        }
    }

    #[test]
    fn viterbi_segments_diploid_run() {
        let mut points: Vec<SamplePoint> = (0..50)
            .map(|i| pt("1", i as u32 * 100 + 1, 0.5, 0.0))
            .collect();
        run_viterbi_per_chrom(&mut points, 0.9999, 1.0, 0.2, 0.04, 0.2);
        let cn2_count = points.iter().filter(|p| p.cn == 2).count();
        assert!(cn2_count >= 40, "expected most CN=2 in diploid run, got {}", cn2_count);
    }

    #[test]
    fn viterbi_detects_duplication_block() {
        let mut points: Vec<SamplePoint> = Vec::new();
        for i in 0..30 {
            points.push(pt("1", i as u32 * 100 + 1, 0.5, 0.0));
        }
        for i in 0..30 {
            // BAF cluster around 1/3 or 2/3 — duplication signature
            let b = if i % 2 == 0 { 1.0 / 3.0 } else { 2.0 / 3.0 };
            points.push(pt("1", 3000 + i as u32 * 100, b, 0.35));
        }
        run_viterbi_per_chrom(&mut points, 0.999, 1.0, 0.2, 0.04, 0.2);
        let dup_count = points[30..].iter().filter(|p| p.cn == 3).count();
        assert!(dup_count >= 20, "expected >=20 CN3 in dup block, got {}", dup_count);
    }
