use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use kira_bio_tools::annotate::{annotate_vcf_ani_v2, build_ani_index_auto_v2, AniIndex};
use std::fs;
use tempfile::TempDir;

fn bench_bgzf_annotation(c: &mut Criterion) {
    let temp_dir = TempDir::new().unwrap();

    let plain_vcf = "annotate.vcf";
    let plain_db = "annotate.db.vcf";

    let bgzf_vcf = temp_dir.path().join("input.vcf.gz");
    let bgzf_db = temp_dir.path().join("db.vcf.gz");

    std::process::Command::new("bgzip")
        .args(&["-c", plain_vcf])
        .stdout(fs::File::create(&bgzf_vcf).unwrap())
        .status()
        .expect("bgzip failed");

    std::process::Command::new("bgzip")
        .args(&["-c", plain_db])
        .stdout(fs::File::create(&bgzf_db).unwrap())
        .status()
        .expect("bgzip failed");

    let ani_path = temp_dir.path().join("db.ani");
    build_ani_index_auto_v2(&bgzf_db, &ani_path).unwrap();

    let mut group = c.benchmark_group("annotation");
    group.throughput(Throughput::Elements(16)); // 16 variants in test file

    group.bench_function(BenchmarkId::new("plain", "vcf"), |b| {
        b.iter(|| {
            let output = temp_dir.path().join("output_plain.vcf");
            annotate_vcf_ani_v2(
                black_box(&ani_path),
                black_box(&plain_vcf.as_ref()),
                black_box(&output),
            )
            .unwrap();
        });
    });

    group.bench_function(BenchmarkId::new("bgzf", "vcf"), |b| {
        b.iter(|| {
            let output = temp_dir.path().join("output_bgzf.vcf");
            annotate_vcf_ani_v2(
                black_box(&ani_path),
                black_box(&bgzf_vcf),
                black_box(&output),
            )
            .unwrap();
        });
    });

    group.finish();
}

fn bench_ani_index_building(c: &mut Criterion) {
    let temp_dir = TempDir::new().unwrap();
    let plain_db = "annotate.db.vcf";
    let bgzf_db = temp_dir.path().join("db.vcf.gz");

    std::process::Command::new("bgzip")
        .args(&["-c", plain_db])
        .stdout(fs::File::create(&bgzf_db).unwrap())
        .status()
        .expect("bgzip failed");

    let mut group = c.benchmark_group("ani_index_build");
    group.throughput(Throughput::Elements(18)); // 18 variants in DB

    group.bench_function(BenchmarkId::new("plain", "vcf"), |b| {
        b.iter(|| {
            let output = temp_dir.path().join("ani_plain.ani");
            build_ani_index_auto_v2(black_box(&plain_db.as_ref()), black_box(&output)).unwrap();
        });
    });

    group.bench_function(BenchmarkId::new("bgzf", "vcf"), |b| {
        b.iter(|| {
            let output = temp_dir.path().join("ani_bgzf.ani");
            build_ani_index_auto_v2(black_box(&bgzf_db), black_box(&output)).unwrap();
        });
    });

    group.finish();
}

criterion_group!(benches, bench_bgzf_annotation, bench_ani_index_building);
criterion_main!(benches);
