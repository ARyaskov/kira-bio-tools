use super::*;

fn two_state() -> Hmm {
    // P(i|j): stay 0.9, switch 0.1.
    Hmm::new(2, &[0.9, 0.1, 0.1, 0.9], 4, ())
}

#[test]
fn viterbi_and_marginals_match_exhaustive_state_paths() {
    let transitions = [0.8, 0.3, 0.2, 0.7];
    let emissions = [0.9, 0.1, 0.25, 0.75, 0.1, 0.9, 0.85, 0.15];
    let sites = [1, 2, 3, 4];
    let prior = [0.8, 0.2];
    let mut total = 0.0;
    let mut marginals = [[0.0; 2]; 4];
    let mut best = (0.0, 0usize);
    for path in 0..32usize {
        let mut probability = prior[path & 1];
        for i in 0..4 {
            let previous = (path >> i) & 1;
            let state = (path >> (i + 1)) & 1;
            probability *= transitions[state * 2 + previous] * emissions[i * 2 + state];
        }
        total += probability;
        if probability > best.0 {
            best = (probability, path);
        }
        for (i, marginal) in marginals.iter_mut().enumerate() {
            marginal[(path >> (i + 1)) & 1] += probability;
        }
    }
    let mut hmm = Hmm::new(2, &transitions, 4, ());
    hmm.init_states(Some(&prior));
    hmm.run_viterbi(&emissions, &sites);
    hmm.run_fwd_bwd(&emissions, &sites);
    for (i, marginal) in marginals.iter().enumerate() {
        assert_eq!(hmm.vstate(i) as usize, (best.1 >> (i + 1)) & 1);
        for (state, expected) in marginal.iter().enumerate() {
            assert!((hmm.fwd_bwd_prob()[(i + 1) * 2 + state] - expected / total).abs() < 1e-12);
        }
    }
}

#[test]
fn tprob_powers_are_precomputed() {
    let h = two_state();
    assert_eq!(h.tprob(), &[0.9, 0.1, 0.1, 0.9]);
    // M^2 = [[0.82, 0.18], [0.18, 0.82]].
    let m2 = &h.tprob_arr[4..8];
    assert!((m2[0] - 0.82).abs() < 1e-12 && (m2[1] - 0.18).abs() < 1e-12);
    let mut h = two_state();
    h.set_curr_tprob(5);
    let expect = &h.tprob_arr[..4];
    // M^6 by direct multiplication.
    let mut m = expect.to_vec();
    let mut tmp = vec![0.0; 4];
    for _ in 0..5 {
        mat_mul(2, expect, &m, &mut tmp);
        m.copy_from_slice(&tmp);
    }
    for i in 0..4 {
        assert!(
            (h.curr_tprob[i] - m[i]).abs() < 1e-12,
            "{:?} vs {:?}",
            h.curr_tprob,
            m
        );
    }
}

#[test]
fn viterbi_and_posteriors_follow_the_emissions() {
    let mut h = two_state();
    // Emissions strongly favour state 0 then state 1.
    let e = [0.9, 0.1, 0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9];
    let sites = [10u32, 11, 12, 13, 14, 15];
    h.run_viterbi(&e, &sites);
    let path: Vec<u8> = (0..6).map(|i| h.vstate(i)).collect();
    assert_eq!((path[0], path[5]), (0, 1), "{path:?}");
    assert!(
        path.windows(2).filter(|w| w[0] != w[1]).count() == 1,
        "{path:?}"
    );
    h.run_fwd_bwd(&e, &sites);
    let fb = h.fwd_bwd_prob();
    // Row 0 is the initial distribution, row i+1 the posterior of site i.
    assert_eq!(&fb[..2], &[0.5, 0.5]);
    assert!(fb[2] > 0.9);
    assert!(fb[6 * 2 + 1] > 0.9);
    for i in 1..=6 {
        assert!((fb[i * 2] + fb[i * 2 + 1] - 1.0).abs() < 1e-12);
    }
}

#[test]
fn snapshot_resumes_a_run() {
    let e = [0.9, 0.1, 0.9, 0.1, 0.9, 0.1, 0.1, 0.9];
    let sites = [1u32, 2, 3, 4];
    let mut full = two_state();
    full.run_viterbi(&e, &sites);
    full.run_fwd_bwd(&e, &sites);
    let want = full.fwd_bwd_prob()[4 * 2..5 * 2].to_vec();

    let mut h = two_state();
    h.request_snapshot(None, 2);
    h.run_viterbi(&e[..4], &sites[..2]);
    h.run_fwd_bwd(&e[..4], &sites[..2]);
    let snap = h.take_snapshot().unwrap();
    assert_eq!(snap.snap_at_pos, 2);
    h.restore(Some(&snap));
    h.run_viterbi(&e[4..], &sites[2..]);
    h.run_fwd_bwd(&e[4..], &sites[2..]);
    let got = h.fwd_bwd_prob()[2 * 2..3 * 2].to_vec();
    // The forward part agrees; the backward part is truncated so only
    // the last site of both runs is comparable.
    assert!((got[0] - want[0]).abs() < 1e-9, "{got:?} vs {want:?}");
    let mut s = snap;
    h.reset(Some(&mut s));
    assert_eq!(s.snap_at_pos, 0);
}

#[test]
fn baum_welch_returns_a_stochastic_matrix() {
    let mut h = two_state();
    let e = [0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.9, 0.1];
    let sites = [1u32, 2, 3, 4, 5];
    let t = h.run_baum_welch(&e, &sites).to_vec();
    for j in 0..2 {
        let col: f64 = (0..2).map(|i| t[2 * i + j]).sum();
        assert!((col - 1.0).abs() < 1e-12, "{t:?}");
    }
}
