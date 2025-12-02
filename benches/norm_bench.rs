use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use kira_bio_tools::norm::{normalize, normalize_scalar, NormContext};

#[cfg(target_arch = "x86_64")]
use kira_bio_tools::norm::normalize_avx2;

#[cfg(target_arch = "aarch64")]
use kira_bio_tools::norm::normalize_neon;

fn bench_normalization_methods(c: &mut Criterion) {
    let test_cases = vec![
        ("small_snp", "A", "T"),
        ("indel_prefix", "ATGC", "ATT"),
        ("indel_suffix", "GCTA", "TTA"),
        ("both_prefix_suffix", "ATGCAT", "ATTCAT"),
        ("long_indel", "ATATATATATATATAT", "ATAT"),
        (
            "very_long",
            "ATGCATGCATGCATGCATGCATGCATGC",
            "ATGCATGCTTTATGCATGCATGCATGC",
        ),
    ];

    let mut group = c.benchmark_group("normalization_comparison");

    for (name, ref_allele, alt_allele) in test_cases.iter() {
        // Scalar baseline
        group.bench_with_input(
            BenchmarkId::new("scalar", name),
            &(ref_allele, alt_allele),
            |b, (r, a)| {
                b.iter(|| {
                    let result = normalize_scalar(black_box(r), black_box(a));
                    black_box(result);
                });
            },
        );

        // Auto-detected (runtime)
        group.bench_with_input(
            BenchmarkId::new("auto", name),
            &(ref_allele, alt_allele),
            |b, (r, a)| {
                b.iter(|| {
                    let result = normalize(black_box(r), black_box(a));
                    black_box(result);
                });
            },
        );

        // AVX2 (x86_64 only)
        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("avx2") {
                group.bench_with_input(
                    BenchmarkId::new("avx2", name),
                    &(ref_allele, alt_allele),
                    |b, (r, a)| {
                        b.iter(|| {
                            let result = normalize_avx2(black_box(r), black_box(a));
                            black_box(result);
                        });
                    },
                );
            }
        }

        // NEON (aarch64 only)
        #[cfg(target_arch = "aarch64")]
        {
            group.bench_with_input(
                BenchmarkId::new("neon", name),
                &(ref_allele, alt_allele),
                |b, (r, a)| {
                    b.iter(|| {
                        let result = normalize_neon(black_box(r), black_box(a));
                        black_box(result);
                    });
                },
            );
        }
    }

    group.finish();
}

fn bench_batch_normalization(c: &mut Criterion) {
    let batch_sizes = vec![100, 1_000, 10_000];
    let ctx = NormContext::detect();

    let mut group = c.benchmark_group("batch_normalization");

    for size in batch_sizes {
        group.throughput(Throughput::Elements(size as u64));

        let refs: Vec<String> = (0..size).map(|i| format!("ATGC{}", i % 10)).collect();

        let alts: Vec<String> = (0..size).map(|i| format!("ATT{}", i % 10)).collect();

        group.bench_with_input(
            BenchmarkId::from_parameter(size),
            &(refs, alts),
            |b, (r, a)| {
                b.iter(|| {
                    for i in 0..r.len() {
                        let result = ctx.normalize(black_box(&r[i]), black_box(&a[i]));
                        black_box(result);
                    }
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_normalization_methods,
    bench_batch_normalization
);
criterion_main!(benches);
