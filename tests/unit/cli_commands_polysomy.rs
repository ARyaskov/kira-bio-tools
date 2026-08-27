    use super::*;

    #[test]
    fn gmm_single_peak_diploid() {
        let mut bafs = Vec::new();
        for k in 0..200 {
            let x = 0.5 + ((k as f64 / 100.0) - 1.0) * 0.05;
            bafs.push(x.clamp(0.05, 0.95));
        }
        let fit = fit_gmm_baf(&bafs);
        assert!(!fit.peaks.is_empty());
        let predicted = infer_cn_from_peaks(&fit.peaks);
        assert_eq!(predicted, 2);
    }

    #[test]
    fn gmm_two_peaks_trisomy() {
        let mut bafs = Vec::new();
        for k in 0..100 {
            let drift = (k as f64 / 100.0 - 0.5) * 0.03;
            bafs.push((1.0 / 3.0 + drift).clamp(0.05, 0.95));
            bafs.push((2.0 / 3.0 + drift).clamp(0.05, 0.95));
        }
        let fit = fit_gmm_baf(&bafs);
        let predicted = infer_cn_from_peaks(&fit.peaks);
        assert_eq!(predicted, 3, "fit peaks={:?}", fit.peaks);
    }
