use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use kira_bio_tools::{GenomicKey, KbiBuilder};
use rand::prelude::*;

fn bench_kbi_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("kbi_build");

    for size in [1_000, 10_000, 100_000, 1_000_000].iter() {
        group.throughput(Throughput::Elements(*size as u64));

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            let mut rng = StdRng::seed_from_u64(42);
            let entries: Vec<(GenomicKey, u64)> = (0..size)
                .map(|i| {
                    let chr = rng.gen_range(1..=22);
                    let pos = rng.gen_range(1..=250_000_000);
                    (GenomicKey::new(chr, pos), i as u64 * 100)
                })
                .collect();

            b.iter(|| {
                let mut builder = KbiBuilder::with_capacity(size);
                for (key, offset) in &entries {
                    builder.add(*key, *offset);
                }
                black_box(builder.build().unwrap())
            });
        });
    }

    group.finish();
}

fn bench_kbi_build_parallel(c: &mut Criterion) {
    let mut group = c.benchmark_group("kbi_build_parallel");

    for size in [100_000, 1_000_000].iter() {
        group.throughput(Throughput::Elements(*size as u64));

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            let mut rng = StdRng::seed_from_u64(42);
            let entries: Vec<(GenomicKey, u64)> = (0..size)
                .map(|i| {
                    let chr = rng.gen_range(1..=22);
                    let pos = rng.gen_range(1..=250_000_000);
                    (GenomicKey::new(chr, pos), i as u64 * 100)
                })
                .collect();

            b.iter(|| {
                let mut builder = KbiBuilder::with_capacity(size);
                for (key, offset) in &entries {
                    builder.add(*key, *offset);
                }
                black_box(builder.build().unwrap())
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_kbi_build, bench_kbi_build_parallel);
criterion_main!(benches);
