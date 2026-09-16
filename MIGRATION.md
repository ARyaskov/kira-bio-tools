# Migration Guide: tabix → kira-bt

## Overview

`kira-bt tabix` is a **100% compatible drop-in replacement** for tabix. This guide helps you migrate existing workflows.

## Quick Start

### Installation

```bash
# Install from crates.io
cargo install kira-bio-tools --features cli

# Or build from source
git clone https://github.com/ARyaskov/kira-bio-tools.git
cd kira-bio-tools
cargo build --release --features cli
```

### Simple Replacement

Replace `tabix` with `kira-bt tabix` in all your scripts:

```bash
# Before
tabix -p vcf file.vcf.gz
tabix file.vcf.gz chr1:1000-2000

# After
kira-bt tabix -p vcf file.vcf.gz
kira-bt tabix file.vcf.gz chr1:1000-2000
```

### Create an Alias

Add to your `.bashrc` or `.zshrc`:

```bash
alias tabix='kira-bt tabix'
```

## Index Formats

`kira-bt tabix` writes `.tbi` unless `-C` asks for `.csi`; `kira-bt index`
writes `.csi` unless `-t` asks for `.tbi`. Both formats are read by htslib
tools (`tabix`, `bcftools view -r`, `bcftools index -n/-s`) and kira-bt reads
indexes made by htslib. Plain gzip and uncompressed inputs are accepted for
streaming; indexing needs BGZF.

## Command Mapping

All tabix commands work identically with `kira-bt tabix`:

| Operation | tabix | kira-bt tabix |
|-----------|-------|---------------|
| Index file | `tabix -p vcf file.vcf.gz` | `kira-bt tabix -p vcf file.vcf.gz` |
| Query region | `tabix file.vcf.gz chr1:1-1000` | `kira-bt tabix file.vcf.gz chr1:1-1000` |
| Multiple regions | `tabix file.vcf.gz chr1:1-1000 chr2:1-1000` | `kira-bt tabix file.vcf.gz chr1:1-1000 chr2:1-1000` |
| With header | `tabix -h file.vcf.gz chr1:1-1000` | `kira-bt tabix -h file.vcf.gz chr1:1-1000` |
| Only header | `tabix -H file.vcf.gz` | `kira-bt tabix -H file.vcf.gz` |
| List chroms | `tabix -l file.vcf.gz` | `kira-bt tabix -l file.vcf.gz` |
| Regions file | `tabix -R regions.bed file.vcf.gz` | `kira-bt tabix -R regions.bed file.vcf.gz` |
| CSI index | `tabix -C file.vcf.gz` | `kira-bt tabix -C file.vcf.gz` |
| Force overwrite | `tabix -f file.vcf.gz` | `kira-bt tabix -f file.vcf.gz` |

## Option Compatibility

### All Short Options Supported

```
-0  --zero-based       Position is 0-based half-open
-b  --begin           Column of start position [4]
-c  --comment         Skip lines starting with character [#]
-C  --csi             Create CSI index instead of TBI
-e  --end             Column of end position [5]
-f  --force           Force overwrite of index
-m  --min-shift       CSI interval size 2^INT [14]
-p  --preset          Format: gff, bed, sam, vcf
-s  --sequence        Column of sequence name [1]
-S  --skip-lines      Skip first INT lines [0]
-h  --print-header    Include header in output
-H  --only-header     Print only header
-l  --list-chroms     List chromosome names
-r  --reheader        Replace header with FILE
-R  --regions         Restrict to regions in FILE
-t  --targets         Sequential scan with targets
-D  --no-download     Don't download remote index
-@  --threads         Number of threads [0]
```

### Extended Options

```
    --cache            BGZF cache size in MB [10]
    --regions-overlap  Show region names (0/1/2)
    --verbosity        Log level (0-4) [3]
```

## Script Migration Examples

### Shell Script

```bash
#!/bin/bash
# Before
for file in *.vcf.gz; do
    tabix -p vcf "$file"
    tabix "$file" chr1:1000000-2000000 > "${file%.vcf.gz}.region.vcf"
done

# After (just add 'kira-bt' before 'tabix')
for file in *.vcf.gz; do
    kira-bt tabix -p vcf "$file"
    kira-bt tabix "$file" chr1:1000000-2000000 > "${file%.vcf.gz}.region.vcf"
done
```

### Python Script

```python
# Before
import subprocess
subprocess.run(['tabix', '-p', 'vcf', 'file.vcf.gz'])
result = subprocess.run(['tabix', 'file.vcf.gz', 'chr1:1000-2000'], 
                       capture_output=True, text=True)

# After
import subprocess
subprocess.run(['kira-bt', 'tabix', '-p', 'vcf', 'file.vcf.gz'])
result = subprocess.run(['kira-bt', 'tabix', 'file.vcf.gz', 'chr1:1000-2000'], 
                       capture_output=True, text=True)
```

### Snakemake Workflow

```python
# Before
rule index_vcf:
    input: "data/{sample}.vcf.gz"
    output: "data/{sample}.vcf.gz.tbi"
    shell: "tabix -p vcf {input}"

# After
rule index_vcf:
    input: "data/{sample}.vcf.gz"
    output: "data/{sample}.vcf.gz.tbi"
    shell: "kira-bt tabix -p vcf {input}"
```

### Nextflow Pipeline

```groovy
// Before
process index_vcf {
    input:
    path vcf
    
    output:
    path "${vcf}.tbi"
    
    script:
    """
    tabix -p vcf ${vcf}
    """
}

// After
process index_vcf {
    input:
    path vcf
    
    output:
    path "${vcf}.tbi"
    
    script:
    """
    kira-bt tabix -p vcf ${vcf}
    """
}
```

## Presets and Custom Layouts

`kira-bt tabix` indexes the same text layouts as tabix:

```bash
kira-bt tabix -p vcf calls.vcf.gz      # VCF: POS and REF length, or INFO/END
kira-bt tabix -p bed regions.bed.gz    # BED: 0-based half-open columns 1-3
kira-bt tabix -p gff genes.gff3.gz     # GFF/GTF: 1-based columns 1, 4, 5
kira-bt tabix -p sam reads.sam.gz      # SAM: RNAME/POS and the CIGAR span
kira-bt tabix -s 2 -b 3 -e 4 table.tsv.gz   # any tab-delimited layout
```

Without `-p` the layout is guessed from the file name (`.vcf`, `.bed`,
`.gff`/`.gff3`/`.gtf`, `.sam`, with or without `.gz`); `-s/-b/-e` select a
generic layout (defaults 1, 4, 5 as in tabix), and `-0`, `-c` and `-S` adjust
whichever applies. The `.tbi` records the layout, so later queries, `-h`,
`-H` and `-l` read the file with the same columns and comment character.

## Performance

Indexing is single-threaded, as in tabix, and runs at a comparable speed;
`-@/--threads` is accepted for command-line compatibility and currently has
no effect. Queries read the same `.tbi`/`.csi` chunks tabix reads. No extra
index files are created.

## Backward Compatibility

### Index Files

kira-bt reads and writes standard `.tbi` and `.csi` files:

- **`.tbi`** files created by kira-bt work with tabix
- **`.csi`** files created by kira-bt work with tabix
- **`.tbi/.csi`** files created by tabix work with kira-bt

## Troubleshooting

### Command Not Found

```bash
# Ensure kira-bt is in PATH
export PATH="$HOME/.cargo/bin:$PATH"

# Or use full path
/path/to/kira-bt tabix file.vcf.gz chr1:1000-2000
```

### Index Already Exists

```bash
# Use -f to force overwrite
kira-bt tabix -f -p vcf file.vcf.gz
```

### File Not BGZF Compressed

```bash
# Compress with bgzip first
bgzip input.vcf
kira-bt tabix -p vcf input.vcf.gz
```

## Extended Features

Beyond tabix compatibility, kira-bt offers extended commands:

```bash
# CSI by default, .tbi with -t; -n/-s report record counts from the index
kira-bt index file.vcf.gz

# bcftools-style region queries and formatted output
kira-bt query -r chr1:1-1000000 -f '%CHROM\t%POS\n' file.vcf.gz

# bcftools view on an indexed file
kira-bt view -r chr1:1-1000000 file.vcf.gz
```

## Testing Your Migration

### Verify Installation

```bash
kira-bt --version
kira-bt tabix --help
```

### Test Basic Functionality

```bash
# Create test file
echo -e "chr1\t1000\t1001\tA\tT" | bgzip > test.vcf.gz

# Index
kira-bt tabix -p vcf test.vcf.gz

# Query
kira-bt tabix test.vcf.gz chr1:1000-1001

# Verify index exists
ls -lh test.vcf.gz.tbi
```

### Compare Output

```bash
# Run both tools on same file
tabix file.vcf.gz chr1:1000-2000 > tabix.out
kira-bt tabix file.vcf.gz chr1:1000-2000 > kira.out

# Compare (should be identical)
diff tabix.out kira.out
```

## Support

For issues or questions:

- GitHub Issues: https://github.com/ARyaskov/kira-bio-tools/issues
- Documentation: https://github.com/ARyaskov/kira-bio-tools#readme

## License

kira-bt is released under the MIT License