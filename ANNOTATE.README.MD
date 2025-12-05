## Annotation Database Formats

`kira-bt` supports two annotation database formats:

### 1. VCF Format (Standard)
Standard VCF files with header:
```vcf
##fileformat=VCFv4.1
##INFO=<ID=T_INT,Number=.,Type=Integer,Description="Test Integer">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	3000150	snp_3000150	C	T	999	PASS	T_INT=1,2;T_FLOAT=1e-10,2e-10
```

### 2. TAB Format (bcftools-compatible)
Tab-delimited format (9 columns, no header):
```
CHROM	POS	REF	ALT	ID	QUAL	T_INT	T_FLOAT	T_STR/FLAG
1	3000150	C	T	snp_3000150	999	1,2	1e-10,2e-10	.
```

**Column mapping:**
- Cols 1-6: Standard VCF fields
- Col 7: `T_INT` values (comma-separated or `.`)
- Col 8: `T_FLOAT` values (comma-separated or `.`)
- Col 9: `T_STR` value or `1` for INDEL flag

## Usage

### Auto-detection (Recommended)
```bash
# Works with both VCF and TAB
kira-bt annotate-index database.tab -o db.ani
kira-bt annotate-index database.vcf -o db.ani
```

### Annotation
```bash
kira-bt annotate -a db.ani input.vcf -o output.vcf
# Or let it build ANI on-the-fly:
kira-bt annotate -a database.tab input.vcf -o output.vcf
```

### Examples from bcftools test suite
```bash
# Download bcftools test data
git clone https://github.com/samtools/bcftools.git
cd bcftools/test

# Index TAB database
kira-bt annotate-index annotate.tab -o annotate.ani

# Annotate
kira-bt annotate -a annotate.ani annotate.vcf -o result.vcf

# Verify compatibility
diff result.vcf annotate.out
```

## Performance
TAB format parsing is ~10% faster than VCF due to simpler structure.