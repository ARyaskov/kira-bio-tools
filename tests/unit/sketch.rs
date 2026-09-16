use super::*;

fn lcg(seed: &mut u64) -> u32 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 33) as u32
}

/// `founders` unrelated samples plus one noisy copy of each (5 % genotype
/// errors, 10 % missing) at `sites` sites with random allele frequencies.
fn cohort(founders: usize, sites: usize, seed: u64) -> Vec<Vec<i8>> {
    let mut rng = seed;
    let n = 2 * founders;
    let mut mat = vec![vec![-1i8; n]; sites];
    for row in mat.iter_mut() {
        let af = 0.05 + (lcg(&mut rng) % 900) as f64 / 1000.0;
        for f in 0..founders {
            let mut d = 0i8;
            for _ in 0..2 {
                if (lcg(&mut rng) % 10_000) as f64 / 10_000.0 < af {
                    d += 1;
                }
            }
            row[f] = d;
            let r = lcg(&mut rng) % 100;
            row[founders + f] = if r < 10 {
                -1
            } else if r < 15 {
                (d + 1) % 3
            } else {
                d
            };
        }
    }
    mat
}

#[test]
fn copies_are_each_others_nearest_sketch() {
    let founders = 25;
    let mat = cohort(founders, 3000, 11);
    let mut sk = DosageSketch::new(2 * founders, 256);
    for row in &mat {
        sk.add_site(row);
    }
    assert_eq!(sk.bits(), 256);
    assert!(sk.sites() > 2500, "{} sites used", sk.sites());
    let sketches = sk.finish();
    assert_eq!(sketches.len(), 2 * founders);
    assert_eq!(sketches[0].len(), 4);
    let pairs = nearest_pairs(&sketches, 1);
    let want: Vec<(usize, usize)> = (0..founders).map(|f| (founders + f, f)).collect();
    assert_eq!(pairs, want);
    // A copy is far closer to its founder than to anyone else.
    for f in 0..founders {
        let own = hamming(&sketches[f], &sketches[founders + f]);
        for g in 0..2 * founders {
            if g != f && g != founders + f {
                assert!(
                    own < hamming(&sketches[f], &sketches[g]),
                    "founder {f} vs {g}"
                );
            }
        }
    }
    // Deterministic and symmetric.
    assert_eq!(sk.finish(), sketches);
    assert_eq!(
        hamming(&sketches[1], &sketches[2]),
        hamming(&sketches[2], &sketches[1])
    );
}

#[test]
fn uninformative_sites_are_skipped_and_k_widens_the_candidate_set() {
    let mut sk = DosageSketch::new(3, 1);
    assert_eq!(sk.bits(), 64);
    sk.add_site(&[1, 1, 1]);
    sk.add_site(&[-1, 2, -1]);
    sk.add_site(&[0, -1, -1]);
    assert_eq!(sk.sites(), 0);
    sk.add_site(&[0, 2, -1]);
    assert_eq!(sk.sites(), 1);

    let mat = cohort(6, 500, 3);
    let mut sk = DosageSketch::new(12, 128);
    for row in &mat {
        sk.add_site(row);
    }
    let s = sk.finish();
    let one = nearest_pairs(&s, 1);
    let three = nearest_pairs(&s, 3);
    assert!(one.len() <= three.len());
    assert!(one.iter().all(|p| three.contains(p)));
    assert!(three.iter().all(|&(a, b)| b < a));
    assert_eq!(nearest_pairs(&s, 0), Vec::<(usize, usize)>::new());
    assert_eq!(nearest_pairs(&s, 100).len(), 12 * 11 / 2);
}
