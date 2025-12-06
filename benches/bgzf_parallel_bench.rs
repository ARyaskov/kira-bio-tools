use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};
use kira_bio_tools::bgzf_parallel::*;
use std::io::Write;
use tempfile::NamedTempFile;

fn bench_parallel_decompression(c: &mut Criterion) {
    let mut tmp = NamedTempFile::new().unwrap();
    let path = tmp.path().to_owned();

    // Create 100MB test file
    {
        let mut writer = ParallelBgzfWriter::create(&path).unwrap();
        let line = "chr1\t1000\trs123\tA\tT\t30\tPASS\tDP=10;AF=0.5\tGT\t0/1\n";
        for _ in 0..1_000_000 {
            writer.write_all(line.as_bytes()).unwrap();
        }
        writer.finish().unwrap();
    }

    let file_size = std::fs::metadata(&path).unwrap().len();

    let mut group = c.benchmark_group("bgzf_parallel");
    group.throughput(Throughput::Bytes(file_size));

    group.bench_function("parallel_read", |b| {
        b.iter(|| {
            let mut reader = ParallelBgzfReader::open(&path).unwrap();
            let mut total = 0usize;

            loop {
                let blocks = reader.read_batch().unwrap();
                if blocks.is_empty() {
                    break;
                }
                total += blocks.iter().map(|b| b.uncompressed.len()).sum::<usize>();
            }

            black_box(total);
        });
    });

    group.finish();
}

criterion_group!(benches, bench_parallel_decompression);
criterion_main!(benches);
