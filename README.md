# Kira Bio Tools

[![CI](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml/badge.svg)](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

High-performance VCF indexer for bioinformatics using learned indexes.

## Features

- **O(1) lookup** - Minimal perfect hashing guarantees constant-time queries
- **100M+ variants** - Scales to large cohort studies
- **Memory efficient** - ~10-12 bytes per variant
- **Instant load** - Memory-mapped index files
- **Cross-platform** - Linux, macOS, Windows

## Quick Start

### Library Usage

```rust
use kira_bio_tools::{VcfIndex, VcfIndexBuilder, GenomicKey, chr_name_to_id};

// Build index
let mut builder = VcfIndexBuilder::new();
builder.add(GenomicKey::new(1, 12345), 1024)?;  // chr1:12345 at byte 1024
builder.add(GenomicKey::new(1, 67890), 2048)?;  // chr1:67890 at byte 2048
let index = builder.build()?;

// Query
if let Some(offset) = index.get(GenomicKey::new(1, 12345)) {
    println!("Found at byte offset: {}", offset);
}

// Save/Load
index.save("variants.kbi")?;
let index = VcfIndex::load_mmap("variants.kbi")?;
```

### CLI Usage

```bash
# Build index from VCF
kira-vcf index -i variants.vcf -o variants.kbi

# Query positions
kira-vcf query -i variants.kbi -p "chr1:12345,chr2:67890"

# Query with VCF line retrieval
kira-vcf query -i variants.kbi -v variants.vcf -p "chr1:12345"

# Range query
kira-vcf range -i variants.kbi -c chr1 -s 10000 -e 20000

# Statistics
kira-vcf stat -i variants.kbi
```

## Installation

### From Source

```bash
# Library only
cargo build --release

# With CLI
cargo build --release --features cli
```

### As Dependency

```toml
[dependencies]
kira-bio-tools = { git = "https://github.com/ARyaskov/kira-bio-tools" }
```

## Performance

Benchmarks on Intel i7-12700, 32GB DDR5, NVMe:

| Metric | Value |
|--------|-------|
| Build rate | ~2M variants/sec |
| Query time | ~50ns |
| Memory | ~10 bytes/variant |
| Load time (mmap) | <1ms |

## File Format

Index files (`.kbi`) use a simple binary format:

```
[Header: 64 bytes]
  - Magic: "KBIV0001"
  - Version, endianness
  - Entry count, MPH parameters
  - Section offsets

[MPH displacement table]
[Genomic keys (u64)]
[VCF byte offsets (u64)]
```

## Architecture

```
┌─────────────────────────────────────────────┐
│                VcfIndex                      │
├─────────────────────────────────────────────┤
│  GenomicKey                                  │
│  ┌─────────────────────────────────────┐    │
│  │ u64 = (chr_id << 32) | position     │    │
│  └─────────────────────────────────────┘    │
│                    │                         │
│                    ▼                         │
│  ┌─────────────────────────────────────┐    │
│  │ kira_kv_engine::Mphf (BDZ MPH)      │    │
│  │ - O(1) key → index mapping          │    │
│  │ - ~2.5 bytes/key overhead           │    │
│  └─────────────────────────────────────┘    │
│                    │                         │
│                    ▼                         │
│  ┌─────────────────────────────────────┐    │
│  │ offsets[index] → VCF byte offset    │    │
│  └─────────────────────────────────────┘    │
└─────────────────────────────────────────────┘
```

## License

MIT License - see [LICENSE](LICENSE) for details.