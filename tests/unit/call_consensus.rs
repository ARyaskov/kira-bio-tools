use super::*;

#[test]
fn allele_count_posterior_matches_exhaustive_genotypes() {
    let likelihoods = [[0.2, 0.7, 0.1], [0.1, 0.3, 0.8]];
    let mut counts = [0.0; 5];
    for a in 0..3 {
        for b in 0..3 {
            let heterozygotes = (a == 1) as u32 + (b == 1) as u32;
            counts[a + b] +=
                likelihoods[0][a] * likelihoods[1][b] * 2f64.powi(heterozygotes as i32);
        }
    }
    let binomial = [1.0, 4.0, 6.0, 4.0, 1.0];
    let mut prior = [0.0; 5];
    for k in 0..4 {
        prior[k] = 1e-3 / (4 - k) as f64;
    }
    prior[4] = 1.0 - prior[..4].iter().sum::<f64>();
    let posterior: Vec<_> = (0..5).map(|k| counts[k] / binomial[k] * prior[k]).collect();
    let normalization = posterior.iter().sum::<f64>();
    let result = spectrum(&likelihoods, false);
    assert!((result.reference - posterior[4] / normalization).abs() < 1e-12);
    let frequency = posterior
        .iter()
        .enumerate()
        .map(|(k, probability)| k as f64 * probability)
        .sum::<f64>()
        / normalization
        / 4.0;
    assert!((result.frequency - frequency).abs() < 1e-12);
}

#[test]
fn incomplete_beta_has_uniform_distribution_limits() {
    for x in [0.0, 0.1, 0.5, 0.9, 1.0] {
        assert!((beta(1.0, 1.0, x) - x).abs() < 1e-12);
    }
}
