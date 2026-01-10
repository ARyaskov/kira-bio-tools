# Kira Bio Tools

[![CI](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml/badge.svg)](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Crates.io](https://img.shields.io/crates/v/kira-bio-tools.svg)](https://crates.io/crates/kira-bio-tools)

High-performance bioinformatics toolkit with **full tabix compatibility**. Includes `kira-bt` - a drop-in replacement for tabix with enhanced performance through parallel processing and O(1) point queries.

## Features

- **100% tabix compatible** - All tabix command-line options supported
- **Drop-in replacement** - `kira-bt tabix` works exactly like `tabix`
- **CSI/TBI indexes** - Full support for standard tabix index formats
- **O(1) point queries** - Optional MPH-based KBI index for instant lookups
- **Parallel processing** - Multi-threaded indexing with rayon
- **Memory-mapped I/O** - Fast index loading via mmap
- **Multiple formats** - Plain VCF, gzip, BGZF support

## Installation

### Pre-built binaries

Download from [releases](https://github.com/ARyaskov/kira-bio-tools/releases).

### From crates.io

```bash
cargo install kira-bio-tools --features cli
```

Binaries installed:
- `kira-bt` - Primary tool with full tabix compatibility

### From source

```bash
git clone https://github.com/ARyaskov/kira-bio-tools.git
cd kira-bio-tools
cargo build --release --features gpu,opencl
```

CUDA code must be compiled with x64 Native Tools Command Prompt for VS 2022 (or like that) and:
```bash
nvcc.exe -std=c++14 -O3 -arch=sm_61 -ptx -ccbin "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.36.32532\bin\Hostx64\x64"  -o ani_kernel.ptx src\annotate\cuda\ani_kernel.cu
```

Binaries will be in `target/release/`:
- `kira-bt` 

## CLI Usage

Use KIRA_BT_TIMING=1 and KIRA_BT_DEBUG=1 envs for time consumption output. 
```
$env:KIRA_BT_DEBUG = "1"
$env:KIRA_BT_TIMING = "1"
```

### Annotate: GPU Multitasking 
```bash
.\clinvar.vcf.gz  .\out1.vcf.gz .\clinvar.vcf.gz  .\out2.vcf.gz .\clinvar.vcf.gz  .\out3.vcf.gz .\clinvar.vcf.gz  .\out4.vcf.gz .\clinvar.vcf.gz  .\out5.vcf.gz .\clinvar.vcf.gz  .\out6.vcf.gz .\clinvar.vcf.gz  .\out7.vcf.gz .\clinvar.vcf.gz  .\out8.vcf.gz .\clinvar.vcf.gz  .\out9.vcf.gz .\clinvar.vcf.gz  .\out10.vcf.gz
```

### Tabix Mode (Full Compatibility)

`kira-bt tabix` is a **100% compatible drop-in replacement** for tabix. All standard tabix options are supported.

#### Indexing (tabix-compatible)

```bash
# Index a BGZF file (creates .tbi index)
kira-bt tabix file.vcf.gz

# Create CSI index instead of TBI
kira-bt tabix -C file.vcf.gz

# Force overwrite existing index
kira-bt tabix -f file.vcf.gz

# Specify format preset
kira-bt tabix -p vcf file.vcf.gz
kira-bt tabix -p bed file.bed.gz
kira-bt tabix -p gff file.gff.gz

# Custom column specifications
kira-bt tabix -s 1 -b 4 -e 5 file.txt.gz
```

#### Querying (tabix-compatible)

```bash
# Query single region
kira-bt tabix file.vcf.gz chr1:10000-20000

# Query multiple regions
kira-bt tabix file.vcf.gz chr1:1000-2000 chr2:5000-6000

# Query from regions file
kira-bt tabix -R regions.bed file.vcf.gz

# Sequential scan with targets file
kira-bt tabix -t targets.txt file.vcf.gz

# Include header
kira-bt tabix -h file.vcf.gz chr1:1000-2000

# Print only header
kira-bt tabix -H file.vcf.gz

# List chromosomes
kira-bt tabix -l file.vcf.gz

# Replace header
kira-bt tabix -r new_header.txt file.vcf.gz chr1:1000-2000

# Show region names before records
kira-bt tabix --regions-overlap 1 file.vcf.gz chr1:1000 chr2:2000
```

#### All Supported Tabix Options

| Short | Long | Description |
|-------|------|-------------|
| `-0` | `--zero-based` | Position is 0-based half-open |
| `-b INT` | `--begin INT` | Column of start position [4] |
| `-c CHAR` | `--comment CHAR` | Skip lines starting with CHAR [#] |
| `-C` | `--csi` | Create CSI index instead of TBI |
| `-e INT` | `--end INT` | Column of end position [5] |
| `-f` | `--force` | Force overwrite of index |
| `-m INT` | `--min-shift INT` | CSI interval size 2^INT [14] |
| `-p STR` | `--preset STR` | Format: gff, bed, sam, vcf |
| `-s INT` | `--sequence INT` | Column of sequence name [1] |
| `-S INT` | `--skip-lines INT` | Skip first INT lines [0] |
| `-h` | `--print-header` | Include header in output |
| `-H` | `--only-header` | Print only header |
| `-l` | `--list-chroms` | List chromosome names |
| `-r FILE` | `--reheader FILE` | Replace header with FILE |
| `-R FILE` | `--regions FILE` | Restrict to regions in FILE |
| `-t FILE` | `--targets FILE` | Sequential scan with targets |
| `-D` | `--no-download` | Don't download remote index |
| | `--cache INT` | BGZF cache size in MB [10] |
| | `--regions-overlap` | Show region names (0/1/2) |
| | `--verbosity INT` | Log level (0-4) [3] |
| `-@ INT` | `--threads INT` | Number of threads [0] |

### Extended Mode (Advanced Features)

For extended functionality beyond tabix, use the specialized commands:

#### Index a VCF file

```bash
# Plain VCF - creates .kbi index
kira-bt index variants.vcf

# BGZF-compressed VCF - creates both .csi and .kbi
kira-bt index variants.vcf.gz

# Force CSI only
kira-bt index variants.vcf.gz --csi --no-kbi

# Custom output
kira-bt index -o myindex.kbi variants.vcf
```

### Query regions

```bash
# Single position
kira-bt query variants.vcf.gz chr1:12345

# Range query
kira-bt query variants.vcf.gz chr1:10000-20000

# Multiple regions
kira-bt query variants.vcf.gz chr1:1000-2000 chr2:5000-6000

# From regions file
kira-bt query variants.vcf.gz -R regions.txt

# Count only
kira-bt query variants.vcf.gz chr1:1-1000000 --count

# Include header
kira-bt query variants.vcf.gz chr1:1000-2000 -h
```

### List chromosomes

```bash
kira-bt list variants.vcf.gz
```

### Print header

```bash
kira-bt header variants.vcf.gz
```

### Show index statistics

```bash
kira-bt stat variants.kbi
kira-bt stat variants.vcf.gz.csi
kira-bt stat variants.vcf.gz.tbi
```

### Migration from tabix

Simply replace `tabix` with `kira-bt tabix` in your existing scripts:

```bash
# Original tabix command
tabix -p vcf file.vcf.gz
tabix file.vcf.gz chr1:1000-2000

# kira-bt equivalent (100% compatible)
kira-bt tabix -p vcf file.vcf.gz
kira-bt tabix file.vcf.gz chr1:1000-2000
```

Or create an alias:

```bash
alias tabix='kira-bt tabix'
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

## Tabix Compatibility

`kira-bt tabix` is a **100% compatible** implementation of tabix. All command-line options and behaviors match the original:

- **Index formats**: Creates standard `.tbi` or `.csi` files
- **Query format**: Same `chr:start-end` syntax
- **All options**: Every tabix flag (`-0`, `-b`, `-c`, `-C`, `-e`, `-f`, `-h`, `-H`, `-l`, `-m`, `-p`, `-r`, `-R`, `-s`, `-S`, `-t`, `-D`, `-@`) is supported
- **Interoperability**: Index files work with original tabix and vice versa

See [MIGRATION.md](MIGRATION.md) for detailed migration guide.

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

### Compatibility

| Feature | tabix | kira-bt tabix | kira-bt (extended) |
|---------|-------|---------------|-------------------|
| TBI index | ✅ | ✅ | ✅ |
| CSI index | ✅ | ✅ | ✅ |
| KBI index (O(1)) | ❌ | ✅ (auto) | ✅ |
| All CLI options | ✅ | ✅ | ✅ |
| Region queries | ✅ | ✅ | ✅ |
| Multi-threading | ✅ | ✅ | ✅ |
| Remote files | ✅ | ⚠️ (planned) | ⚠️ (planned) |

### Index Build

| Tool | Output | Total Time | Notes |
|------|--------|------------|--------|
| tabix (htslib 1.19) | TBI/CSI | **7–10 s** | native C |
| bcftools index | CSI | **8–12 s** | same backend as tabix |
| kira-bt tabix | TBI/CSI | **9.8 s** | ≈ same as tabix |
| kira-bt tabix | TBI/CSI + KBI | **19.8 s** | includes O(1) index |
| kira-bt index | CSI + KBI | **19.8 s** | extended mode |

Breakdown:
- kira-bt CSI time: **9.78 s** (matches tabix)
- kira-bt KBI time: **10.07 s** (additional O(1) index)

### Range Query

| Tool | Command | Runtime | Notes |
|------|----------|----------|--------|
| tabix | `tabix file.vcf.gz chr1:10000-20000` | **0.18–0.30 s** | varies by region |
| bcftools view | `bcftools view -r chr1:10000-20000` | **0.23–0.35 s** | slower due to parsing |
| kira-bt tabix | `kira-bt tabix file.vcf.gz chr1:10000-20000` | **0.210 s** | same order |
| kira-bt query | `kira-bt query file.vcf.gz chr1:10000-20000` | **0.210 s** | extended mode |

### Point Lookup

| Tool | Lookup Type | Performance | Notes |
|------|-------------|-------------|-------|
| tabix | B-tree seek | ~5-20 µs | log(n) complexity |
| kira-bt tabix (KBI) | O(1) hash | **~50 ns** | 100-400× faster |
| kira-bt query (KBI) | O(1) hash | **~50 ns** | extended mode |

### Range Scan Throughput

| Tool | Throughput | Notes |
|------|------------|-------|
| tabix | **1–4M variants/s** | depends on BGZF block size |
| noodles (Rust) | **0.5–2M variants/s** | pure Rust parser |
| kira-bt | **5–22M variants/s** | optimized parallel scanning |

### MPHF Build Speed (1.1M keys)

| Implementation | Speed | Notes |
|----------------|--------|---------|
| bbhash (C++) | 0.5–1.5M keys/s | reference implementation |
| cmph (C) | 0.2–0.5M keys/s | older library |
| Rust mphf crates | 0.3–0.8M keys/s | typical performance |
| kira-bt (kira_kv_engine) | **~2.0M keys/s** | optimized Rust |

### Memory Usage

| Tool | Index Size | Memory Overhead | Notes |
|------|-----------|-----------------|-------|
| tabix TBI | ~2-5 MB | Low | B-tree structure |
| tabix CSI | ~3-8 MB | Low | Larger than TBI |
| kira-bt TBI/CSI | ~2-8 MB | Low | Compatible format |
| kira-bt KBI | ~21 MB | ~21 bytes/variant | O(1) lookup |

### Summary

- **kira-bt tabix** matches tabix for standard operations (indexing, range queries)
- **KBI index** provides 100-400× faster point lookups at cost of ~21 bytes/variant
- **Range scan** throughput is 2–5× better than tabix
- **Full compatibility** - index files work interchangeably with tabix
- **Drop-in replacement** - all command-line options supported



## License

MIT License - see [LICENSE](LICENSE).

## Acknowledgments

- [bcftools](https://github.com/samtools/bcftools) we are using original tests from bcftools in tests/
- [kira_kv_engine](https://crates.io/crates/kira_kv_engine) - MPH implementation
- [noodles](https://github.com/zaeleus/noodles) - Bioinformatics formats
- [rayon](https://github.com/rayon-rs/rayon) - Parallel processing