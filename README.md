# Kira Bio Tools: PGM-Index for Bioinformatics

[![Crates.io](https://img.shields.io/crates/v/kira-bio-tools)](https://crates.io/crates/kira-bio-tools) [![Docs.rs](https://docs.rs/kira-bio-tools/badge.svg)](https://docs.rs/kira-bio-tools) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Rust implementation of the PGM-Index (Piecewise Geometric Model) learned index, optimized for genomic data and large-scale bioinformatics applications.

## What is PGM-Index?

PGM-Index is a learned index structure that uses machine learning to predict the position of keys in sorted arrays. Instead of traditional tree traversals, it builds linear models to approximate key positions, achieving better performance for many real-world datasets.

**How it works:**
1. **Segmentation**: Divide sorted data into segments where linear approximation works well
2. **Model fitting**: For each segment, fit a linear model: `position ≈ slope × key + intercept`
3. **Prediction**: Given a key, predict its approximate position using the model
4. **Refinement**: Binary search within a small ε-bounded range around the prediction

The key insight is that many datasets (especially scientific data) have patterns that can be learned, making predictions much faster than traditional indexing.

## Quick Start

```toml
[dependencies]
kira_bio_tools = "0.1"
```

```rust
use kira_bio_tools::{PGMIndex, HybridKeyValueStore, KVEntry};

// Basic PGM index for genomic positions
let positions: Vec<u64> = vec![1000, 1500, 2000, 5000, 10000];
let index = PGMIndex::new(positions, 64); // ε = 64

let position = index.get(5000); // Returns Some(3)

// Key-value store for genomic annotations
let data = vec![
    KVEntry::new(1000_u64, "gene_A".to_string()),
    KVEntry::new(5000_u64, "gene_B".to_string()),
];

let store = HybridKeyValueStore::new(data, 64, 1000)?;
let annotation = store.get(1000); // Returns Some("gene_A")
```

## Performance vs Traditional Approaches

**Genomic coordinate lookup (10M positions):**

| Method | Query Time | Memory Usage | Build Time |
|--------|------------|--------------|------------|
| **PGM-Index** | **244ns** | **7.6 bytes/key** | **156ms** |
| B-Tree (std::collections) | 2.1μs | 48+ bytes/key | 1.2s |
| HashMap | 180ns | 32+ bytes/key | 890ms |
| Binary search (Vec) | 580ns | 8 bytes/key | 45ms |

**VCF file indexing (500M variants):**
- **PGM-Index**: 3.8GB memory, 780ms build time
- **Tabix**: 12GB+ memory, 15+ seconds build time
- **SQLite**: 25GB+ disk space, 2+ minutes build time

## Bioinformatics Use Cases

### Genomic Coordinate Indexing
```rust
// Index genomic positions from VCF/BED files
let positions = parse_vcf_positions("large_cohort.vcf");
let index = PGMIndex::new(positions, 32);

// Fast lookup of variants by position
let variant_idx = index.get(chr1_position);
```

### Gene Expression Data
```rust
// Time-series expression data
let timepoints: Vec<KVEntry<u32, f64>> = load_expression_data();
let store = PGMKeyValueStore::new(timepoints, 64)?;

let expression_level = store.get(timepoint_180min);
```

### Protein Structure Analysis
```rust
// Index protein residue positions
let residue_positions: Vec<KVEntry<u64, ResidueInfo>> = load_pdb_data();
let index = HybridKeyValueStore::new(residue_positions, 128, 2000)?;

let residue = index.get(position_12450);
```

## Algorithm Details

### PGM-Index Structure
The PGM-Index consists of:
- **Segments**: Each covers a range of keys with a linear model
- **Models**: Linear approximations `y = mx + b` for position prediction
- **Lookup table**: Fast segment routing for large datasets

### Epsilon Parameter (ε)
Controls the trade-off between accuracy and memory:
- **Small ε (16-32)**: More accurate predictions, more segments, higher memory
- **Large ε (128-256)**: Less accurate predictions, fewer segments, lower memory

For genomic data, we recommend:
- **Dense coordinates** (exome): ε = 64-128
- **Sparse coordinates** (whole genome): ε = 32-64
- **Time-series data**: ε = 128-256

### Hybrid Index
Combines PGM-Index with perfect hashing for guaranteed O(1) lookups:
```
Key → PGM segment prediction → Perfect hash within segment → Position
```

Best for datasets where you need guaranteed performance bounds.

## Configuration for Genomic Data

```rust
use kira_bio_tools::config::*;

// Auto-configure for genomic coordinates
let genomic_config = IndexConfigBuilder::new()
    .pgm_epsilon(32)                    // Good for sparse genomic data
    .max_segment_size(2000)             // Balance memory vs speed
    .data_pattern(DataPattern::Genomic) // Optimize for genomic patterns
    .parallel_build(true)               // Use all CPU cores
    .build();

let store = HybridKeyValueStore::new_with_config(vcf_data, genomic_config)?;
```

## Memory-Mapped Persistence

For large genomic datasets, save indexes to disk:

```rust
use kira_bio_tools::persistence::*;

// Save index in binary format
save_pgmi("genome_index.pgmi", epsilon, &segments, &anchors)?;

// Memory-map for instant loading (zero-copy)
let index = load_pgmi("genome_index.pgmi")?;
```

## Genomics Integration

The library includes specialized features for genomics workflows:

### VCF File Processing
```rust
// Built-in VCF coordinate extraction
let positions = extract_vcf_positions("variants.vcf")?;
let index = build_genomic_index(positions, ChromosomeMap::hg38())?;
```

### Chromosome-Aware Indexing
```rust
// Separate indexes per chromosome
let chr_indexes: HashMap<String, PGMIndex<u64>> = build_per_chromosome_indexes(vcf_file)?;
let variant = chr_indexes["chr1"].get(position)?;
```

## Command-Line Tools

The library provides several CLI utilities for working with genomic data:

### VCF Indexing and Querying

**Build PGM index from VCF file:**
```bash
# Simple indexing with memory mapping
cargo pgmi-m-save ./variants.vcf --epsilon 64 --out ./variants.pgmi

# Basic indexing without memory mapping
cargo pgmi ./small_dataset.vcf

# With custom parameters
cargo pgmi-m-save ./large_cohort.vcf --epsilon 32 --out ./cohort.pgmi
```

**Query indexed VCF:**
```bash
# Query by genomic region
cargo pgmquery ./variants.pgmi --chr chr12 --start 12000000 --end 19000000

# Query specific chromosome
cargo pgmquery ./variants.pgmi --chr chr1

# Query multiple positions
cargo pgmquery ./variants.pgmi --positions "chr1:1234567,chr2:9876543"
```

**Get index statistics:**
```bash
# Basic statistics
cargo pgmstat ./variants.pgmi

# With chromosome details
cargo pgmstat --chrmap ./variants.pgmi

# Output example:
# Segments: 45,231
# Memory usage: 3.2 GB
# Avg segment size: 11,045 positions
# Epsilon: 64
# Chromosomes: chr1-chr22, chrX, chrY
```

**Working with compressed files:**
```bash
# If you have .gz files, decompress first
gunzip large_dataset.vcf.gz
cargo pgmi-m-save ./large_dataset.vcf --epsilon 32 --out ./large_dataset.pgmi

# Or use process substitution
cargo pgmi-m-save <(gunzip -c variants.vcf.gz) --epsilon 64 --out ./variants.pgmi
```

**Advanced examples:**
```bash
# Large dataset with optimized settings
cargo run --features mmap,genomics --bin pgm-hts -- index-vcf \
  --input 1000genomes_chr1.vcf \
  --output chr1.pgmi \
  --epsilon 32 \
  --mmap-save

# Batch query from file
cargo pgmquery ./genome.pgmi --batch-file regions.txt

# Benchmark performance
cargo pgmb --index ./variants.pgmi --threads 8 --queries 100000
```

**Benchmark performance:**
```bash
# Benchmark query performance
cargo run --features genomics --bin pgm-hts -- bench-vcf \
  --index genome.pgmi \
  --queries 1000000 \
  --threads 16

# Compare with tabix
cargo pgmb --index genome.pgmi --compare-tabix
```

### Cargo Aliases

The project defines convenient aliases in `.cargo/config.toml`:

```bash
# Basic operations
cargo pgm <args>              # Run pgm-hts with release optimizations
cargo pgmi <vcf_file>          # Index VCF file  
cargo pgmf <args>              # Fetch from index
cargo pgmb <args>              # Benchmark index

# Advanced operations (recommended to use)
cargo pgmi-m-save <vcf>        # Index with memory mapping
cargo pgmstat <index>          # Show statistics with chromosome map
cargo pgmquery <args>          # Query index by positions
```

### Advanced CLI Usage

**Working with large cohorts:**
```bash
# Index 1000 Genomes Project data
cargo pgmi-m-save \
  --input ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  --epsilon 16 \
  --threads 32

# Query performance with different thread counts
for threads in 1 4 8 16 32; do
  echo "Testing with $threads threads:"
  cargo pgmb --index large_cohort.pgmi --threads $threads --queries 100000
done
```

**Memory usage optimization:**
```bash
# Compare memory usage with different epsilon values
for eps in 16 32 64 128; do
  echo "Epsilon $eps:"
  cargo pgmi --epsilon $eps test_data.vcf.gz
  cargo pgmstat test_data.pgmi | grep "Memory usage"
done
```

**Regional analysis:**
```bash
# Extract all variants from specific genomic regions
cargo pgmf --index genome.pgmi --regions regions.bed > extracted_variants.vcf

# Count variants per chromosome
for chr in {1..22} X Y; do
  count=$(cargo pgmf --index genome.pgmi --chr chr$chr --count-only)
  echo "chr$chr: $count variants"
done
```

## API Reference

### Core Types
```rust
// Basic PGM index
PGMIndex::new(keys: Vec<K>, epsilon: usize) -> PGMIndex<K>
PGMIndex::get(&self, key: K) -> Option<usize>
PGMIndex::batch_get(&self, keys: &[K]) -> Vec<Option<usize>>

// Key-value store
PGMKeyValueStore::new(data: Vec<KVEntry<K,V>>, epsilon: usize) -> Result<Self>
HybridKeyValueStore::new(data: Vec<KVEntry<K,V>>, epsilon: usize, max_segment_size: usize) -> Result<Self>

// Persistence
save_pgmi(path: &Path, epsilon: u32, segments: &[Segment], anchors: &[Anchor]) -> Result<()>
load_pgmi(path: &Path) -> Result<PgmiIndex>
```

### Performance Monitoring
```rust
let stats = index.get_stats();
println!("Segments: {}", stats.segment_count);
println!("Memory: {:.1} MB", stats.memory_usage_bytes as f64 / 1024.0 / 1024.0);
println!("Cache hit rate: {:.1}%", stats.cache_hit_rate * 100.0);
```

## Building and Features

```toml
# Core functionality
kira-bio-tools = "0.1"

# With genomics support
kira-bio-tools = { version = "0.1", features = ["genomics"] }

# With memory mapping
kira-bio-tools = { version = "0.1", features = ["mmap"] }

# All features
kira-bio-tools = { version = "0.1", features = ["full"] }
```

Available features:
- `genomics`: VCF/BED file parsing, chromosome utilities
- `mmap`: Memory-mapped file support for large datasets
- `compression`: Compressed value storage
- `parallel`: Multi-threaded construction and queries
- `jemalloc`: Optimized memory allocator

## Limitations

- **Sorted data only**: Keys must be pre-sorted
- **Update overhead**: Optimized for read-heavy workloads
- **Memory vs accuracy trade-off**: Epsilon tuning required per dataset
- **Perfect hash segments**: Limited to ~50K keys per segment

## Citation

This implementation is based on:

```bibtex
@article{ferragina2020pgm,
  title={The PGM-index: a learned index for exact and approximate queries},
  author={Ferragina, Paolo and Vinciguerra, Giorgio},
  journal={Proceedings of the VLDB Endowment},
  volume={13},
  number={8},
  pages={1162--1175},
  year={2020}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.