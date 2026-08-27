## Annotation Database Formats

`kira-bt` supports annotation databases in multiple formats:

### 1. VCF Format (Plain or BGZF)
Standard VCF files with header:
```vcf
##fileformat=VCFv4.1
##INFO=<ID=T_INT,Number=.,Type=Integer,Description="Test Integer">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	3000150	snp_3000150	C	T	999	PASS	T_INT=1,2;T_FLOAT=1e-10,2e-10
```

Both plain `.vcf` and BGZF-compressed `.vcf.gz` are supported.

### 2. TAB Format (Plain)
Tab-delimited format (9 columns, no header):
```
CHROM	POS	REF	ALT	ID	QUAL	T_INT	T_FLOAT	T_STR/FLAG
1	3000150	C	T	snp_3000150	999	1,2	1e-10,2e-10	.
```

## Usage

### Build ANI Index (Auto-detection)
```bash
# Works with plain VCF, BGZF VCF, or TAB
kira-bt annotate-index database.vcf -o db.ani
kira-bt annotate-index database.vcf.gz -o db.ani
kira-bt annotate-index database.tab -o db.ani
```

### Annotation
```bash
# Both input and database support plain/BGZF auto-detection
kira-bt annotate -a db.ani input.vcf -o output.vcf
kira-bt annotate -a db.ani input.vcf.gz -o output.vcf
kira-bt annotate -a database.vcf.gz input.vcf.gz -o output.vcf
```

### Performance Comparison

| Tool | Format | Throughput | Notes |
|------|--------|------------|-------|
| bcftools annotate | BGZF | 1-4M var/s | Standard baseline |
| kira-bt annotate (plain) | Plain VCF | 15-20M var/s | No compression overhead |
| kira-bt annotate (BGZF) | BGZF | 5-15M var/s | 3-5× faster than bcftools |

### BGZF Optimizations

kira-bt uses several optimizations for BGZF files:

1. **Batched parallel processing**: 200K variants per batch
2. **Lock-free output pipeline**: Minimal contention
3. **Optimized INFO merging**: Reduced allocations
4. **Native BGZF support**: No external dependencies

### Example: bcftools Compatible Workflow
```bash
# Index annotation database (BGZF)
kira-bt annotate-index annotate.db.vcf.gz -o annotate.ani

# Annotate input (BGZF)
kira-bt annotate -a annotate.ani input.vcf.gz -o output.vcf

# Performance timing
KIRA_BT_TIMING=1 kira-bt annotate -a annotate.ani input.vcf.gz -o output.vcf
```

### Memory Usage

| Component | Size | Notes |
|-----------|------|-------|
| ANI Index | ~21 bytes/variant | Loaded into memory |
| Input Buffer | 64MB | For output writing |
| Batch Buffer | ~50MB | For parallel processing |
| Total | <150MB | For 1M variant file |

### Compatibility

- **Input formats**: Plain VCF, BGZF VCF (auto-detected)
- **Database formats**: Plain VCF, BGZF VCF, TAB (auto-detected)
- **Output format**: Same as input (plain or BGZF)
- **bcftools compatibility**: Full command-line compatibility