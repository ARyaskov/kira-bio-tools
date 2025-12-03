# Annotate


## Annotate — High‑Speed INFO Field Annotation

`kira-bt annotate` provides a **drop‑in replacement for ****bcftools annotate -a**, but with a dramatically faster architecture based on:

* **ANI (Annotation Index)** — custom, ultra‑compact MPH‑based annotation index
* **O(1) lookup per variant** (hash‑based, not binary‑search‑based)
* **Zero decompression** and **zero parsing** during annotation
* **Streaming VCF processing** with negligible RAM usage

### Why an Annotation Index (ANI)?

Traditional `bcftools annotate -a db.vcf` performs:

1. Text parsing of the annotation database
2. Position lookup via B‑tree or CSI
3. Per‑record string comparisons
4. INFO merge

This becomes slow for large VCFs.

`kira-bt` instead pre‑builds an annotation-only MPH index:

```text
ANI = MPH(key = CHR + POS + REF + ALT) → INFO
```

This gives:

* **Instant lookups** (50–80 ns each)
* **No database parsing at annotation time**
* **Small index** (≈ 1.5–3× smaller than equivalent KBI)
* **100% deterministic matching** (REF/ALT aware)

### Building the Annotation Index

```bash
kira-bt annotate-index big_db.vcf
```

Output:

* `big_db.ani` — annotation index

You can also specify output file manually:

```bash
kira-bt annotate-index -o my_index.ani big_db.vcf
```

Example (real timings on Apple M1 8GB):

```
[annotate-index] Input  = "big_db.vcf"
[annotate-index] Output = "big_db.ani"
[ani] reading db VCF...
[ani] building MPH...
[ani] DONE: 10000 variants
Total time: 0.011s
```

### Performing Annotation

```bash
kira-bt annotate -a big_db.vcf input.vcf -o output.vcf
```

`kira-bt` automatically resolves:

* if `big_db.ani` exists → use instant ANI lookup
* otherwise → full VCF-to-VCF slow annotation

Explicit ANI usage:

```bash
kira-bt annotate -a big_db.ani input.vcf -o output.vcf
```

Annotation progress is reported:

```
[annotate] ANI = "big_db.ani"
[annotate] Loading ANI index...
[annotate] Annotating...
[annotate] processed=100000 annotated=2314 0.040s
...
[annotate] DONE: processed=636106 annotated=15422 total=0.223s
```

### Performance Characteristics

#### Annotation Speed

Dataset: VCF 27 MB → 636k records

| Tool                            | Annotation DB  | Time             |
| ------------------------------- | -------------- | ---------------- |
| bcftools annotate               | 433 KB VCF     | ~4–8 seconds     |
| kira-bt annotate (VCF mode)     | 433 KB VCF     | ~1.8 seconds     |
| **kira-bt annotate (ANI mode)** | **big_db.ani** | **0.22 seconds** |

Speedup vs bcftools: **≈20–40× faster**

#### Annotation of Large VCFs (11 GB)

```
[annotate] ... processed=1103547 ... total=67.48s
```

Scaling: **~15–17 million variants/minute** on Apple M1 Air 2020.

This is extraordinarily fast — comparable to optimized C++ pipelines.

### Example Pipeline

```bash
# 1. Normalize (optional)
kira-bt norm sample.vcf -o sample.norm.vcf

# 2. Index annotation DB
kira-bt annotate-index db.vcf -o db.ani

# 3. Annotate
kira-bt annotate -a db.ani sample.norm.vcf -o sample.annot.vcf

# 4. (Optional) tabix-index again if needed
kira-bt tabix sample.annot.vcf.gz
```

### INFO Merge Logic

`kira-bt annotate` strictly follows **bcftools default behavior**:

* INFO fields from annotation DB are appended to the input record
* If the key already exists → **value is replaced**
* Merging does **not** touch FORMAT, QUAL, FILTER
* ALT-aware fields (Number=A/R/G) match exactly on REF and ALT

### Limitations

* Does not yet support FORMAT-level annotation
* Does not yet support gVCFs
* Does not support multi-allelic records in DB (planned)

### Summary

`kira-bt annotate` with ANI indexing delivers:

* **Superb performance** (10–40× faster than bcftools)
* **Near-zero memory usage**
* **Deterministic REF/ALT-aware merging**
* **A compact MPH annotation index format (ANI)**
* **Perfect reproducibility for large-scale pipelines**

