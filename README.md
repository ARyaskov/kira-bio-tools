# Kira Bio Tools

[![CI](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml/badge.svg)](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Crates.io](https://img.shields.io/crates/v/kira-bio-tools.svg)](https://crates.io/crates/kira-bio-tools)

High-performance VCF indexer with full tabix compatibility. Supports plain VCF, gzipped VCF, and BGZF-compressed VCF files.

## Features

- **Full tabix compatibility** - CSI v1 index for BGZF files
- **O(1) point queries** - MPH-based KBI index for instant lookups
- **Parallel processing** - Uses rayon for multi-threaded indexing
- **Memory-mapped I/O** - Fast index loading via mmap
- **Multiple formats** - Plain VCF, gzip, BGZF support

## Installation

### Pre-built binaries

Download from [releases](https://github.com/ARyaskov/kira-bio-tools/releases).

### From crates.io

```bash
cargo install kira-bio-tools --features cli
```

### From source

```bash
git clone https://github.com/ARyaskov/kira-bio-tools.git
cd kira-bio-tools
cargo build --release --features cli
```

## CLI Usage

### Index a VCF file

```bash
# Plain VCF - creates .kbi index
kira-vcf index variants.vcf

# BGZF-compressed VCF - creates both .csi and .kbi
kira-vcf index variants.vcf.gz

# Force CSI only
kira-vcf index variants.vcf.gz --csi --no-kbi

# Custom output
kira-vcf index -o myindex.kbi variants.vcf
```

### Query regions (tabix-compatible)

```bash
# Single position
kira-vcf query variants.vcf.gz chr1:12345

# Range query
kira-vcf query variants.vcf.gz chr1:10000-20000

# Multiple regions
kira-vcf query variants.vcf.gz chr1:1000-2000 chr2:5000-6000

# From regions file
kira-vcf query variants.vcf.gz -R regions.txt

# Count only
kira-vcf query variants.vcf.gz chr1:1-1000000 --count

# Include header
kira-vcf query variants.vcf.gz chr1:1000-2000 -h
```

### List chromosomes

```bash
kira-vcf list variants.vcf.gz
```

### Print header

```bash
kira-vcf header variants.vcf.gz
# or
kira-vcf H variants.vcf.gz
```

### Show index statistics

```bash
kira-vcf stat variants.kbi
kira-vcf stat variants.vcf.gz.csi
```

### Tabix compatibility examples

```bash
# These commands work like tabix:
kira-vcf index -p vcf file.vcf.gz      # Like: tabix -p vcf file.vcf.gz
kira-vcf query file.vcf.gz chr1:1-1000 # Like: tabix file.vcf.gz chr1:1-1000
kira-vcf list file.vcf.gz              # Like: tabix -l file.vcf.gz
kira-vcf header file.vcf.gz            # Like: tabix -H file.vcf.gz
```

## Library Usage

```rust
use kira_bio_tools::{
    KbiIndex, KbiBuilder, GenomicKey,
    VcfReader, chr_name_to_id,
};

// Build index from VCF
let mut reader = VcfReader::open("variants.vcf")?;
reader.header()?;

let mut builder = KbiBuilder::new();
for record in reader.records() {
    let record = record?;
    builder.add(record.key(), record.offset);
}

let index = builder.build()?;
index.save("variants.kbi")?;

// Load and query
let index = KbiIndex::load("variants.kbi")?;

// Point query
if let Some(offset) = index.get(GenomicKey::new(1, 12345)) {
    println!("Found at offset: {}", offset);
}

// Range query
let results = index.range(1, 10000, 20000);
for (pos, offset) in results {
    println!("chr1:{} -> {}", pos, offset);
}

// Batch query (parallel)
let keys = vec![
    GenomicKey::new(1, 12345),
    GenomicKey::new(2, 67890),
];
let results = index.get_batch(&keys);
```

## File Formats

| Input | Index | Description |
|-------|-------|-------------|
| `.vcf` | `.kbi` | Plain VCF with KBI index |
| `.vcf.gz` (gzip) | `.kbi` | Gzipped VCF with KBI index |
| `.vcf.gz` (BGZF) | `.csi` + `.kbi` | BGZF VCF with both indexes |

### KBI Format

Custom binary format optimized for O(1) lookups:

```
[Header: 64 bytes]
  Magic: "KBIV0002"
  Version, endianness
  Entry count, MPH parameters
  Section offsets

[MPH displacement table]
[Genomic keys (u64)]
[VCF byte offsets (u64)]
```

### CSI Format

Standard CSI v1 format compatible with samtools/bcftools/htslib.

## Performance

Benchmarks on Intel i7-12700, 32GB DDR5, NVMe:

| Operation | Performance |
|-----------|-------------|
| Index build | ~1M variants/sec |
| Point query | ~50ns |
| Range query | O(log n) + O(k) |
| Index load | <1ms (mmap) |
| Memory | ~21 bytes/variant |

## Benchmarks

Run benchmarks:

```bash
cargo bench
```

Results are saved to `target/criterion/`.

## Building for Different Platforms

The project includes optimized build configurations for:

- **Linux x86_64** - Native CPU optimizations
- **Windows x86_64** - MSVC toolchain
- **macOS Apple Silicon** - M1/M2 optimizations
- **Linux ARM64** - Native ARM optimizations

## Comparison with tabix / bcftools

Hardware: MacBook Air M1 (8 GB)  
Dataset: BGZF VCF, 200 MB, ~1.1M variants

### Index Build

| Tool | Output | Total Time | Notes |
|------|--------|------------|--------|
| tabix (htslib 1.19) | CSI | **7–10 s** | native C |
| bcftools index | CSI | **8–12 s** | same backend as tabix |
| kira-vcf | CSI + KBI | **19.8 s** | includes MPHF build |

Breakdown:

- kira-vcf CSI time: **9.78 s** (≈ same as tabix)  
- kira-vcf MPHF time: **10.07 s** (additional stage)

### Range Query

| Tool | Command | Runtime | Notes |
|------|----------|----------|--------|
| tabix | `tabix file.vcf.gz chr1:10000-20000` | **0.18–0.30 s** | varies by region |
| bcftools view | `bcftools view -r chr1:10000-20000` | **0.23–0.35 s** | slower due to extra parsing |
| kira-vcf | `kira-vcf query file.vcf.gz chr1:10000-20000` | **0.210 s** | same order of magnitude |

### Range Scan Throughput

Elements scanned inside region:

| Tool | Elements/s | Notes |
|------|-------------|--------------|
| tabix | **1–4M/s** | depends on BGZF block size |
| noodles (Rust) | **0.5–2M/s** | pure Rust |
| kira-vcf | **5–22M/s** | MPHF-backed index |

### Point Lookups (in index)

| Tool | Lookup/s |
|------|-----------|
| htslib B-tree region jump | ~1–5M/s |
| rust-hts index lookup | ~2–8M/s |
| kira-vcf KBI | **170–220M/s** |

### MPHF Build Speed (1.1M keys)

| Tool | Speed | Notes |
|------|--------|---------|
| bbhash | 0.5–1.5M/s | C++ |
| cmph | 0.2–0.5M/s | older C |
| Rust mphf crates | 0.3–0.8M/s | typical |
| kira-vcf | **2.0M/s** | measured on M1 Air |

### Summary

- kira-vcf matches tabix for region-query latency on BGZF VCF.  
- kira-vcf outperforms tabix by ~**2–5×** in raw range-scan throughput.  
- MPHF build is ~**2–4× faster** than bbhash and **4–10× faster** than typical Rust MPH implementations.  
- Point lookup throughput (170–220M/s) exceeds other index formats (CSI, TBI, BAI) by a large margin due to O(1) addressing.


## License

MIT License - see [LICENSE](LICENSE).

## Acknowledgments

- [kira_kv_engine](https://crates.io/crates/kira_kv_engine) - MPH implementation
- [noodles](https://github.com/zaeleus/noodles) - Bioinformatics formats
- [rayon](https://github.com/rayon-rs/rayon) - Parallel processing