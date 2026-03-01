# Kira Bio Tools

[![CI](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml/badge.svg)](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Crates.io](https://img.shields.io/crates/v/kira-bio-tools.svg)](https://crates.io/crates/kira-bio-tools)

`kira-bt` is a high-performance toolkit for VCF/BCF processing with bcftools-oriented command-line compatibility.

## Installation

From crates.io:

```bash
cargo install kira-bio-tools --features cli
```

From source:

```bash
git clone https://github.com/ARyaskov/kira-bio-tools.git
cd kira-bio-tools
cargo build --release
```

Binary name: `kira-bt`

## Runtime Flags

- `KIRA_BT_TIMING=1` enables timing output.
- `KIRA_BT_DEBUG=1` enables debug output.

## CLI Modes and Commands

Top-level entrypoint:

```bash
kira-bt <command> [arguments]
```

For exact command syntax:

```bash
kira-bt --help
kira-bt <command> --help
```

### 1. Annotation and Indexing

- `annotate-index`  
  Build `.ani` index from annotation VCF.  
  Args: `<input>`, `-o/--output`.

- `annotate`  
  Annotate input VCF using annotation source.  
  Required source group: `-a/--annotations <vcf|tab|ani>` or `--ani <ani>`.  
  Common args: `<input>`, `-o/--output`, `-c/--columns`, `--gpu`, `--opencl`, `--bgzf-level`, `--cache-plain`, `--bgzf-after`, `--mmap-output`, `--mmap-no-flush`, `--ram-output`, `--ram-max-mb`.

- `annotate-serve`  
  Server mode for annotate pipeline.  
  Uses the same annotation source group and backend/performance flags as `annotate`.

- `db-build`  
  Build ANI database file from input annotation data.  
  Args: `<input>`, `-o/--output`.

### 2. Tabix and Region Utilities

- `tabix`  
  Tabix-compatible indexing and querying.  
  Args: `<input> [regions...]` and options such as `--csi`, `--preset`, `--begin`, `--end`, `--sequence`, `--regions`, `--targets`, `--print-header`, `--only-header`, `--list-chroms`, `--threads`.

- `region-index` (alias: `kindex`)  
  Extended regional index builder.  
  Args: `<input>`, `-o/--output`, `--preset`, `--force`, `--min-shift`, `--depth`, `--csi`, `--no-kbi`.

- `region-query` (alias: `rquery`)  
  Query indexed files by region list.  
  Args: `<file> [regions...]`, `-R/--regions-file`, `-c/--count`, `-h/--print-header`, `-H/--only-header`.

- `stat`  
  Show index statistics.  
  Args: `<index>`.

- `list`  
  List chromosome names from indexed file.  
  Args: `<file>`.

- `header` (alias: `H`)  
  Print VCF header.  
  Args: `<file>`.

### 3. Variant Processing Commands

- `filter`  
  Native filtering mode with expression and output controls.  
  Args include: `-i/--include`, `-e/--exclude`, `--expr`, `-s/--soft-filter`, `-m/--mode`, `-S/--set-GTs`, `-g/--SnpGap`, `-G/--IndelGap`, `-M/--mask-file`, `-o/--output`, `-O/--output-type`, `--threads`, `--gpu`, `--opencl`, `<input>`.

- `norm`  
  Normalization mode.  
  Args: `<input>`, `-f/--fasta-ref`, `-o/--output`.

- `index` (alias: `idx`)  
  VCF/BCF index command with bcftools-style tail arguments.  
  Pattern: `kira-bt index <input> -- [bcftools-like options]`.

- `query`  
  Transform VCF/BCF into custom text output (bcftools-style compatible args).  
  Pattern: `kira-bt query -- [bcftools-like options]`.

- `head`  
  View VCF/BCF headers (bcftools-style args passthrough).  
  Pattern: `kira-bt head -- [bcftools-like options]`.

- `sort`  
  Sort VCF/BCF input.  
  Args: `<input>`, `-o/--output`, and additional compatibility options after `--`.

- `stats`  
  Produce VCF/BCF stats.  
  Args: `<inputs...>` and additional options after `--`.

### 4. Compatibility Command Surface

The following commands are available in the CLI and use bcftools-style argument patterns (usually with extra args after `--`):

- `call` (input + optional output + trailing compatibility args)
- `mpileup` (inputs + optional output + trailing compatibility args)
- `gtcheck`
- `csq`
- `convert`
- `consensus`
- `concat`
- `roh` (input + trailing compatibility args)
- `reheader` (input + trailing compatibility args)
- `polysomy` (input + trailing compatibility args)
- `merge` (inputs + trailing compatibility args)
- `isec` (inputs + trailing compatibility args)
- `cnv`

## Command-Line Conventions

- Commands with `trailing_var_arg` accept additional compatibility arguments after `--`.
- Typical pattern:

```bash
kira-bt <command> <required-positional-args> -- <extra-options>
```

Examples:

```bash
kira-bt query -- -f '%CHROM\t%POS\n' in.vcf.gz
kira-bt index in.vcf.gz -- --csi -f
kira-bt call in.vcf.gz -o out.vcf -- -mv
```

## Tests

Test suites are grouped per command in `tests/<mode>/testN`.

Current scale: approximately **400+ test cases** (about **405** at the moment).

Each suite contains:
- command scripts (`bcftools.sh`, `kira.sh`)
- input assets
- reference outputs (`out.bcf.ref.vcf`, etc.)
- per-suite documentation in `tests/<mode>/README.md`

## License

MIT, see [LICENSE](LICENSE).
