//! Pair-HMM forward-backward, before and after the structure-of-arrays
//! rewrite, on 150-base reads against a 200-base reference window with
//! mismatches and one indel, i.e. the BAQ workload of `mpileup`.

use criterion::{Criterion, black_box, criterion_group, criterion_main};

#[path = "glocal_old.rs"]
mod glocal_old;

use kira_bio_tools::align::{GlocalParams, GlocalScratch, encode_nt, glocal_with};

fn synth() -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    let bases = *b"ACGT";
    (0..64)
        .map(|_| {
            let refw: Vec<u8> = (0..200).map(|_| bases[(next() % 4) as usize]).collect();
            let start = 20 + (next() % 20) as usize;
            let mut read: Vec<u8> = refw[start..start + 150].to_vec();
            for _ in 0..3 {
                let p = (next() % 150) as usize;
                read[p] = bases[(next() % 4) as usize];
            }
            let p = 40 + (next() % 60) as usize;
            if next() % 2 == 0 {
                read.remove(p);
            } else {
                read.insert(p, b'A');
            }
            let qual: Vec<u8> = (0..read.len()).map(|_| 20 + (next() % 20) as u8).collect();
            let r: Vec<u8> = refw.iter().map(|&b| encode_nt(b)).collect();
            let q: Vec<u8> = read.iter().map(|&b| encode_nt(b)).collect();
            (r, q, qual)
        })
        .collect()
}

fn bench(c: &mut Criterion) {
    let data = synth();
    let conf = GlocalParams {
        d: 0.001,
        e: 0.1,
        bw: 10,
    };
    let old_conf = glocal_old::GlocalParams {
        d: 0.001,
        e: 0.1,
        bw: 10,
    };
    let mut g = c.benchmark_group("glocal_fwd_bwd_64_reads");
    g.bench_function("old", |b| {
        let mut sc = glocal_old::GlocalScratch::default();
        b.iter(|| {
            for (r, q, qual) in &data {
                black_box(glocal_old::glocal_with(
                    &mut sc,
                    r,
                    q,
                    Some(qual),
                    &old_conf,
                    true,
                ));
            }
        })
    });
    g.bench_function("new", |b| {
        let mut sc = GlocalScratch::default();
        b.iter(|| {
            for (r, q, qual) in &data {
                black_box(glocal_with(&mut sc, r, q, Some(qual), &conf, true));
            }
        })
    });
    g.finish();
    // The two must agree bit for bit.
    let mut sc = GlocalScratch::default();
    let mut so = glocal_old::GlocalScratch::default();
    for (r, q, qual) in &data {
        let a = glocal_with(&mut sc, r, q, Some(qual), &conf, true).unwrap();
        let b = glocal_old::glocal_with(&mut so, r, q, Some(qual), &old_conf, true).unwrap();
        assert_eq!(a.loglik.to_bits(), b.loglik.to_bits());
        assert_eq!(a.state, b.state);
        assert_eq!(a.q, b.q);
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
