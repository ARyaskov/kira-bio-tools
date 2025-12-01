use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use kira_bio_tools::{GenomicKey, KbiBuilder, KbiIndex};
use rand::prelude::*;

fn create_test_index(size: usize) -> KbiIndex {
    let mut rng = StdRng::seed_from_u64(42);
    let mut builder = KbiBuilder::with_capacity(size);

    for i in 0..size {
        let chr = rng.gen_range(1..=22);
        let pos = rng.gen_range(1..=250_000_000);
        builder.add(GenomicKey::new(chr, pos), i as u64 * 100);
    }

    builder.build().unwrap()
}

fn bench_kbi_point_query(c: &mut Criterion) {
    let mut group = c.benchmark_group("kbi_point_query");

    for size in [10_000, 100_000, 1_000_000].iter() {
        let index = create_test_index(*size);

        let mut rng = StdRng::seed_from_u64(123);
        let queries: Vec<GenomicKey> = (0..10_000)
            .map(|_| {
                let chr = rng.gen_range(1..=22);
                let pos = rng.gen_range(1..=250_000_000);
                GenomicKey::new(chr, pos)
            })
            .collect();

        group.throughput(Throughput::Elements(10_000));

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("size_{}", size)),
            &(&index, &queries),
            |b, (index, queries)| {
                b.iter(|| {
                    for key in *queries {
                        black_box(index.get(*key));
                    }
                });
            },
        );
    }

    group.finish();
}

fn bench_kbi_range_query(c: &mut Criterion) {
    let mut group = c.benchmark_group("kbi_range_query");

    for size in [10_000, 100_000, 1_000_000].iter() {
        let index = create_test_index(*size);

        group.throughput(Throughput::Elements(1000));

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("size_{}", size)),
            &index,
            |b, index| {
                let mut rng = StdRng::seed_from_u64(456);

                b.iter(|| {
                    for _ in 0..1000 {
                        let chr = rng.gen_range(1..=22);
                        let start = rng.gen_range(1..=240_000_000);
                        let end = start + rng.gen_range(1000..100_000);
                        black_box(index.range(chr, start, end));
                    }
                });
            },
        );
    }

    group.finish();
}

fn bench_kbi_batch_query(c: &mut Criterion) {
    let mut group = c.benchmark_group("kbi_batch_query");

    let size = 1_000_000;
    let index = create_test_index(size);

    for batch_size in [100, 1000, 10000].iter() {
        let mut rng = StdRng::seed_from_u64(789);
        let queries: Vec<GenomicKey> = (0..*batch_size)
            .map(|_| {
                let chr = rng.gen_range(1..=22);
                let pos = rng.gen_range(1..=250_000_000);
                GenomicKey::new(chr, pos)
            })
            .collect();

        group.throughput(Throughput::Elements(*batch_size as u64));

        group.bench_with_input(
            BenchmarkId::from_parameter(batch_size),
            &(&index, &queries),
            |b, (index, queries)| {
                b.iter(|| black_box(index.get_batch(queries)));
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_kbi_point_query,
    bench_kbi_range_query,
    bench_kbi_batch_query
);
criterion_main!(benches);
