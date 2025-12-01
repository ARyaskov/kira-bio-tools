# kira-bt Examples

## Basic Usage

### Indexing

```bash
# Index a BGZF-compressed VCF
kira-bt tabix -p vcf variants.vcf.gz

# Force overwrite existing index
kira-bt tabix -f -p vcf variants.vcf.gz

# Create CSI index instead of TBI
kira-bt tabix -C -p vcf variants.vcf.gz

# Index with custom parameters
kira-bt tabix -s 1 -b 2 -e 3 data.txt.gz
```

### Querying

```bash
# Query single region
kira-bt tabix variants.vcf.gz chr1:1000-2000

# Query multiple regions
kira-bt tabix variants.vcf.gz chr1:1000-2000 chr2:3000-4000 chr3:5000-6000

# Query with header
kira-bt tabix -h variants.vcf.gz chr1:1000-2000

# Query specific position
kira-bt tabix variants.vcf.gz chr1:12345

# Query from file of regions
kira-bt tabix -R regions.bed variants.vcf.gz
```

### Information

```bash
# List all chromosomes
kira-bt tabix -l variants.vcf.gz

# Print only header
kira-bt tabix -H variants.vcf.gz

# Show statistics
kira-bt stat variants.vcf.gz.tbi
kira-bt stat variants.kbi
```

## Advanced Usage

### Working with Different Formats

```bash
# BED files
kira-bt tabix -p bed regions.bed.gz
kira-bt tabix regions.bed.gz chr1:1000-2000

# GFF files
kira-bt tabix -p gff annotations.gff.gz
kira-bt tabix annotations.gff.gz chr1:5000-10000

# SAM files
kira-bt tabix -p sam alignments.sam.gz
kira-bt tabix alignments.sam.gz chr1:1000-2000
```

### Custom Column Specifications

```bash
# Custom tab-delimited file with:
# - Column 1: chromosome
# - Column 4: start position
# - Column 5: end position
kira-bt tabix -s 1 -b 4 -e 5 custom.txt.gz

# 0-based coordinates (UCSC style)
kira-bt tabix -0 -p bed regions.bed.gz

# Skip first 5 lines
kira-bt tabix -S 5 -p vcf variants.vcf.gz

# Use different comment character
kira-bt tabix -c '@' -p sam alignments.sam.gz
```

### Region Files

```bash
# BED format regions file
cat > regions.bed <<EOF
chr1    1000    2000
chr2    5000    6000
chr3    10000   15000
EOF

kira-bt tabix -R regions.bed variants.vcf.gz

# Simple text format (chr:start-end)
cat > regions.txt <<EOF
chr1:1000-2000
chr2:5000-6000
chr3:10000-15000
EOF

kira-bt tabix -R regions.txt variants.vcf.gz
```

### Sequential Scan with Targets

```bash
# Targets file (does not require sorting)
cat > targets.txt <<EOF
chr3:10000-15000
chr1:1000-2000
chr2:5000-6000
EOF

kira-bt tabix -t targets.txt variants.vcf.gz
```

### Header Replacement

```bash
# Create new header
cat > new_header.txt <<EOF
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
EOF

# Use new header when querying
kira-bt tabix -r new_header.txt variants.vcf.gz chr1:1000-2000
```

### Region Overlap Markers

```bash
# Mode 1: Print region as comment before records
kira-bt tabix --regions-overlap 1 variants.vcf.gz chr1:1000 chr2:2000

# Output:
# #chr1:1000
# chr1	1000	.	A	T	...
# #chr2:2000
# chr2	2000	.	G	C	...

# Mode 2: Prepend region to each line (tab-separated)
kira-bt tabix --regions-overlap 2 variants.vcf.gz chr1:1000 chr2:2000

# Output:
# chr1:1000	chr1	1000	.	A	T	...
# chr2:2000	chr2	2000	.	G	C	...
```

### Multi-threading

```bash
# Use 8 threads for operations
kira-bt tabix -@ 8 -p vcf large_file.vcf.gz

# Query with multiple threads
kira-bt tabix -@ 4 large_file.vcf.gz chr1:1-1000000
```

## Pipeline Integration

### Simple Shell Pipeline

```bash
#!/bin/bash
set -euo pipefail

VCF="variants.vcf.gz"

# Index
echo "Indexing..."
kira-bt tabix -p vcf "$VCF"

# Extract regions
echo "Extracting regions..."
for chr in chr{1..22}; do
    kira-bt tabix "$VCF" "${chr}:1-1000000" > "${chr}_region.vcf"
done

# Count variants per chromosome
echo "Counting variants..."
for chr in chr{1..22}; do
    count=$(kira-bt tabix "$VCF" "$chr" | wc -l)
    echo "$chr: $count variants"
done
```

### Processing Multiple Files

```bash
#!/bin/bash
set -euo pipefail

# Index all VCF files
for vcf in *.vcf.gz; do
    echo "Indexing $vcf..."
    kira-bt tabix -f -p vcf "$vcf"
done

# Query all files for same region
REGION="chr1:1000000-2000000"
for vcf in *.vcf.gz; do
    echo "=== $vcf ==="
    kira-bt tabix "$vcf" "$REGION"
done > combined_region.vcf
```

### Parallel Processing

```bash
#!/bin/bash
set -euo pipefail

VCF="large_variants.vcf.gz"
REGIONS="regions.txt"

# Index with 8 threads
kira-bt tabix -@ 8 -p vcf "$VCF"

# Process regions in parallel
cat "$REGIONS" | parallel -j 4 \
    "kira-bt tabix $VCF {} > {/.}.vcf"
```

### Integration with bcftools

```bash
#!/bin/bash
set -euo pipefail

VCF="variants.vcf.gz"

# Index
kira-bt tabix -p vcf "$VCF"

# Extract region and pipe to bcftools
kira-bt tabix "$VCF" chr1:1000000-2000000 | \
    bcftools view -i 'QUAL>30' | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n'
```

### Integration with awk/grep

```bash
#!/bin/bash
set -euo pipefail

VCF="variants.vcf.gz"

# Extract region and filter with awk
kira-bt tabix -h "$VCF" chr1:1000000-2000000 | \
    awk '$6 > 30' > high_quality.vcf

# Count variants by type
kira-bt tabix "$VCF" chr1 | \
    awk '{print length($4), length($5)}' | \
    awk '{if($1==1 && $2==1) print "SNP"; else print "INDEL"}' | \
    sort | uniq -c
```

## Performance Optimization

### Cache Size Tuning

```bash
# Increase cache for large queries
kira-bt tabix --cache 100 large_file.vcf.gz chr1:1-100000000

# Disable cache
kira-bt tabix --cache 0 file.vcf.gz chr1:1000-2000
```

### Batch Queries

```bash
#!/bin/bash
set -euo pipefail

VCF="variants.vcf.gz"

# Single query for multiple regions (faster)
kira-bt tabix "$VCF" \
    chr1:1000-2000 \
    chr1:5000-6000 \
    chr1:10000-11000 \
    > regions.vcf

# Instead of multiple separate queries (slower)
# for region in chr1:1000-2000 chr1:5000-6000 chr1:10000-11000; do
#     kira-bt tabix "$VCF" "$region"
# done
```

### Extended Mode for O(1) Lookups

```bash
# Create both TBI and KBI indexes
kira-bt index variants.vcf.gz

# Point queries use O(1) KBI index automatically
time kira-bt query variants.vcf.gz chr1:12345

# Range queries still use TBI/CSI
time kira-bt query variants.vcf.gz chr1:10000-20000
```

## Error Handling

### Check if Index Exists

```bash
#!/bin/bash
set -euo pipefail

VCF="variants.vcf.gz"
INDEX="${VCF}.tbi"

if [ ! -f "$INDEX" ]; then
    echo "Index not found, creating..."
    kira-bt tabix -p vcf "$VCF"
fi

kira-bt tabix "$VCF" chr1:1000-2000
```

### Handle Missing Files

```bash
#!/bin/bash
set -euo pipefail

VCF="variants.vcf.gz"

if [ ! -f "$VCF" ]; then
    echo "Error: $VCF not found"
    exit 1
fi

kira-bt tabix -p vcf "$VCF" || {
    echo "Error: Failed to index $VCF"
    exit 1
}
```

### Verbose Logging

```bash
# Enable debug output
kira-bt tabix --verbosity 4 -p vcf variants.vcf.gz

# Only errors
kira-bt tabix --verbosity 1 -p vcf variants.vcf.gz

# Silent mode
kira-bt tabix --verbosity 0 -p vcf variants.vcf.gz
```

## Comparison with Original tabix

All these commands work identically with both tools:

```bash
# Original tabix
tabix -p vcf file.vcf.gz
tabix file.vcf.gz chr1:1000-2000
tabix -h file.vcf.gz chr1:1000-2000
tabix -l file.vcf.gz
tabix -R regions.bed file.vcf.gz

# kira-bt (identical behavior)
kira-bt tabix -p vcf file.vcf.gz
kira-bt tabix file.vcf.gz chr1:1000-2000
kira-bt tabix -h file.vcf.gz chr1:1000-2000
kira-bt tabix -l file.vcf.gz
kira-bt tabix -R regions.bed file.vcf.gz
```

The index files are fully compatible:

```bash
# Create index with tabix
tabix -p vcf file.vcf.gz

# Query with kira-bt (uses tabix index)
kira-bt tabix file.vcf.gz chr1:1000-2000

# Create index with kira-bt
kira-bt tabix -p vcf file2.vcf.gz

# Query with original tabix (uses kira-bt index)
tabix file2.vcf.gz chr1:1000-2000
```