# Kira Bio Tools

[![CI](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml/badge.svg)](https://github.com/ARyaskov/kira-bio-tools/actions/workflows/ci.yml)
[![License: Apache](https://img.shields.io/badge/License-Apache-violet.svg)](LICENSE)
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
- `KIRA_BT_NO_PROGRESS=1` suppresses BAM-reading progress bars.
- `KIRA_BT_BGZF_THREADS=<n>` sets the BGZF decompression worker count.
- `KIRA_SOLID_DUMP_BAM=<path>` writes the `solid` pipeline's post-markdup BAM.
- `KIRA_NM_WEIGHT` / `KIRA_INDEL_REALIGN` supply the default for `mpileup`'s
  `--nm-weight` / `--indel-realign` (and `solid`'s `--mpileup-` forms).

Other `KIRA_*` variables in the source are internal tuning knobs and are not
part of the supported interface.

## CLI Modes and Commands

Top-level entrypoint:

```bash
kira-bt <command> [arguments]
```

The lists below name the positional arguments and the options you are most
likely to need. For the complete, authoritative option set of any command:

```bash
kira-bt --help
kira-bt <command> --help
```

### Argument conventions

Commands come in three shapes.

**1. Native flags.** Most commands parse the bcftools options directly. A
bare `--` in the command line is ignored for these commands, so bcftools
invocations written as `kira-bt cmd ... -- <options>` keep working:

```bash
kira-bt call in.vcf.gz -o out.vcf -m -v
kira-bt view in.bcf -r chr1:1-10000 -O z -o region.vcf.gz
kira-bt stats -- -s - in.vcf.gz
```

Commands: `annotate-index`, `annotate`, `annotate-serve`, `db-build`, `tabix`,
`region-index`, `region-query`, `stat`, `list`, `header`, `head`, `norm`,
`filter`, `gtcheck`, `csq`, `consensus`, `concat`, `call`, `mpileup`, `stats`,
`merge`, `isec`, `view`, `samview`, `solid`, `tile`.

**2. Native flags plus a passthrough tail.** These accept extra
bcftools-compatible arguments after `--`:

```bash
kira-bt index in.vcf.gz -- --tbi -f
kira-bt sort in.vcf -o out.vcf -- --temp-dir /scratch
```

Commands: `index`, `convert`, `sort`, `roh`, `reheader`, `polysomy`.

Streaming: every reader accepts `-` (or no path where the input is optional)
for stdin and detects plain, gzip, BGZF, VCF and BCF by content; `-O z`, `-O b`
and `-O u` write to stdout when no `-o` is given.

**3. Free-form bcftools argv.** These take the whole bcftools command line as
positional arguments, with no native options of their own:

```bash
kira-bt query -f '%CHROM\t%POS\n' in.vcf.gz
kira-bt cnv --output-dir outdir in.vcf
```

Commands: `query`, `cnv`.

### 1. Fused Pipeline

- `solid`  
  FASTQ → VCF in one process: align → coordinate-sort → markdup → mpileup →
  call, entirely in memory. The output VCF is the only file written to disk.

  Required: `--aligner-ref <REF>`, `--aligner-r1 <R1>`, `-o/--output <VCF>`.

  Aligner stage: `--aligner-r2`, `--aligner-index`, `--aligner-rg`,
  `--aligner-insert-size`, `--aligner-batch`.

  In-memory BAM stage: `--bam-markdup` (off by default — PCR duplicates
  otherwise count as independent support).

  Pileup stage: `--mpileup-annotate` (default `AD,DP,SP`; `PL` is always
  emitted), `--mpileup-max-depth`, `--mpileup-min-mq` (default 0),
  `--mpileup-min-bq` (default 13), `--mpileup-nm-weight` and
  `--mpileup-indel-realign` (the `mpileup` extensions below, off by default).
  `--mpileup-variants-only` is on in this build and cannot currently be
  switched off, so ref-only positions are always skipped.

  Call stage: `--call` enables it (off by default, since it changes what the
  output VCF means), `--call-ploidy <0|1|2>`, `--call-ploidy-file`,
  `--call-consensus`, `--call-annotate <GQ[,GP]>`, `--call-prior`.

  Memory: `--window-mb <MB>` processes the genome in reference windows instead
  of holding every alignment at once, with `--window-tmpdir` for the spill
  directory. The default (0) is fastest but peak memory scales with the input.

  ```bash
  kira-bt solid \
    --aligner-ref ref.fa --aligner-r1 r1.fq.gz --aligner-r2 r2.fq.gz \
    --aligner-rg 'ID:s1\tSM:sample' \
    --bam-markdup --mpileup-min-mq 20 --mpileup-min-bq 20 \
    --call --call-ploidy 1 \
    -t 16 -o variants.vcf
  ```

  Without `--call` the output carries the pileup's per-site maximum-likelihood
  genotypes, not a called set. `solid` covers alignment through calling only:
  normalization, filtering, annotation and QC remain separate steps.

  Debug aid: `KIRA_SOLID_DUMP_BAM=<path>` writes the post-markdup BAM so a fused
  run can be compared against a step-by-step one. Non-windowed runs only.

### 2. Annotation and Indexing

- `annotate-index`  
  Build a `.ani` index from an annotation VCF.  
  Args: `<input>`, `-o/--output`.

- `annotate`  
  Annotate an input VCF from an annotation source.  
  Required source: `-a/--annotations <vcf|tab|ani>` or `--ani <ani>`.  
  Common: `<input>`, `-o/--output`, `-O/--output-type`, `-c/--columns`,
  `-C/--columns-file`, `-h/--header-lines`, `-H/--header-line`, `-i/--include`,
  `-e/--exclude`, `-I/--set-id`, `-m/--mark-sites`, `-x/--remove`,
  `-r/--regions`, `-R/--regions-file`, `-s/--samples`, `-S/--samples-file`,
  `--pair-logic`, `--rename-chrs`, `--rename-annots`, `--threads`,
  `-W/--write-index`.  
  Backend and I/O tuning: `--gpu`, `--bgzf-level`, `--cache-plain`,
  `--bgzf-after`, `--mmap-output`, `--mmap-no-flush`, `--ram-output`,
  `--ram-max-mb`, `--no-ktile`, `--no-build-ktile`, `--force-build-ktile`.

- `annotate-serve`  
  Server mode driven by stdin commands. Same annotation source group and
  backend flags as `annotate`, plus `-c/--columns` as the default column spec.

- `db-build`  
  Build an ANI database from annotation data.  
  Args: `<input>`, `-o/--output`.

- `tile`  
  Build and manage `.ktile` sidecars that let `annotate` skip the BGZF decode
  and per-line parse pass.  
  Subcommand: `kira-bt tile build -i/--input <input> [-o/--output <ktile>]`.

### 3. Tabix, Regions and Headers

- `tabix`  
  Tabix-compatible indexing and querying. Writes `.tbi` by default and `.csi`
  with `-C`; both are readable by htslib, and htslib-made indexes are used for
  queries. Regions use tabix semantics (`chr:beg-end`, `-0`, `-R/-T` region
  files).  
  Args: `<input> [regions...]`.  
  Options: `-p/--preset`, `-C/--csi`, `-m/--min-shift`, `-s/--sequence`,
  `-b/--begin`, `-e/--end`, `-c/--comment`, `-S/--skip-lines`, `-0/--zero-based`,
  `-f/--force`, `-r/--reheader`, `-R/--regions`, `-t/--targets`,
  `-h/--print-header`, `-H/--only-header`, `-l/--list-chroms`, `-D/--no-download`,
  `--cache`, `-@/--threads`.

- `index` (alias `idx`)  
  VCF/BCF indexing. Writes an htslib-compatible `.csi` by default (`-t` for
  `.tbi`); `-n`/`-s` read record counts from an existing index. Index entries
  span `POS..POS+len(REF)-1` (or `INFO/END`), so overlap queries behave like
  bcftools.  
  Args: `<input> [-- <passthrough>...]`.  
  Options: `-t/--tbi`, `-c/--csi`, `-m/--min-shift`, `-f/--force`,
  `-n/--nrecords`, `-s/--stats`, `--all`, `-o/--output`, `--threads`.

- `region-index` (alias `kindex`)  
  Extended regional index builder.  
  Args: `<input>`, `-o/--output`, `-p/--preset`, `-f/--force`, `-s/--min-shift`,
  `-d/--depth`, `-C/--csi`, `--no-kbi`.

- `region-query` (alias `rquery`)  
  Query an indexed file by region list.  
  Args: `<file> [regions...]`, `-R/--regions-file`, `-c/--count`,
  `-h/--print-header`, `-H/--only-header`.

- `stat`  
  Index statistics.  
  Args: `<index>` (`.kbi` or `.csi`).

- `list`  
  Chromosome names from an indexed file.  
  Args: `<file>`.

- `header` (alias `H`)  
  Print the VCF header.  
  Args: `<file>`.

- `head`  
  View leading header and record lines.  
  Args: `<input>`, `-h/--headers`, `-n/--records`, `-s/--samples`.

### 4. Variant Processing

- `view`  
  Convert, view, subset and filter VCF/BCF.  
  Args: `<input>`, `-o/--output`, `-O/--output-type`, `-l/--compression-level`,
  `-r/--regions`, `-R/--regions-file`, `-t/--targets`, `-T/--targets-file`,
  `-s/--samples`, `-S/--samples-file`, `--force-samples`, `-i/--include`,
  `-e/--exclude`, `-a/--trim-alt-alleles`, `-I/--no-update`,
  `-G/--drop-genotypes`, `-h/--header-only`, `-H/--no-header`,
  `-f/--apply-filters`, `-v/--types`, `-V/--exclude-types`, `--known`,
  `--novel`, `-m/--min-alleles`, `-M/--max-alleles`, `-g/--genotype`,
  `-p/--phased`, `-u/--uncalled`, `-x/--private`, `-c/--min-ac`, `-C/--max-ac`,
  `-q/--min-af`, `-Q/--max-af`, `--threads`, `-W/--write-index`.

- `filter`  
  Expression filtering with soft-filter and gap controls.  
  Args: `<input>`, `-i/--include`, `-e/--exclude`, `--expr`, `-s/--soft-filter`,
  `-m/--mode`, `-S/--set-GTs`, `-g/--SnpGap`, `-G/--IndelGap`, `--mask`,
  `-M/--mask-file`, `--mask-overlap`, `-o/--output`, `-O/--output-type`,
  `--threads`, `--gpu`, `-W/--write-index`.

- `norm` (alias `N`)  
  Left-align, split/join multiallelics, check REF, remove duplicates.  
  Args: `<input>`, `-f/--fasta-ref`, `--fasta-ref-fai`, `-o/--output`,
  `-O/--output-type`, `-m/--multiallelics`, `--atomize`, `--atom-overlaps`,
  `-c/--check-ref`, `-N/--do-not-normalize`, `-d/--rm-dup`, `--site-win`,
  `--strict-filter`, `--keep-sum`, `--old-rec-tag`, `-S/--sort`, plus the usual
  region/target/include/exclude selectors.

- `sort`  
  Sort VCF/BCF.  
  Args: `<input> [-- <passthrough>...]`, `-o/--output`, `-O/--output-type`,
  `-m/--max-mem`, `-T/--temp-dir`, `--threads`, `-W/--write-index`.

- `concat`  
  Concatenate files sharing a sample set. `-a` walks overlapping inputs the
  way the bcftools synced reader does (chromosomes in order of first
  appearance, identical records grouped per input, `-d`/`-D` collapsing);
  `-l` ligates phased chunks with bcftools' half-overlap rule, `PQ` at the
  switch site and `PS` phase sets (`--ligate-warn`, `--ligate-force`);
  `-n` copies BGZF blocks without decoding.  
  Args: `<inputs...>`, `-f/--file-list`, `-o/--output`, `-O/--output-type`,
  `-a/--allow-overlaps`, `-l/--ligate`, `-D/--remove-duplicates`,
  `-d/--rm-dups`, `-n/--naive`, `-G/--drop-genotypes`, `-c/--compact-PS`,
  `-q/--min-PQ`, region selectors, `--threads`, `-W/--write-index`.

- `merge`  
  Merge files with non-overlapping sample sets. Records at one position are
  grouped by the `-m` rules (SNVs first, then indels), contigs follow header
  order, QUAL is the maximum, FILTER follows `-F`, INFO rules default to
  `DP:sum,DP4:sum` (plus `AN:sum,AC:sum` for sites-only output), Number=A/R/G
  tags are remapped to the merged alleles and AN/AC are recomputed from the
  genotypes; `--gvcf` fills absent samples from reference blocks.  
  Args: `<inputs...>`, `-l/--file-list`, `-o/--output`, `-O/--output-type`,
  `-m/--merge`, `-f/--apply-filters`, `-F/--filter-logic`, `--force-samples`,
  `-i/--info-rules`, `--missing-rules`, `-0/--missing-to-ref`, `-g/--gvcf`,
  `-L/--local-alleles`, region selectors, `--threads`, `-W/--write-index`.

- `isec`  
  Intersections and complements with bcftools' file layout: `-p` writes
  `0000.vcf`, `0001.vcf`, ..., `sites.txt` and `README.txt`; without `-p` the
  site list (or the `-w` file) goes to stdout.  
  Args: `<inputs...>`, `-p/--prefix`, `-n/--nfiles`, `-w/--write`,
  `-c/--collapse`, `-C/--complement`, `-f/--apply-filters`, `-l/--file-list`,
  `-o/--output`, `-O/--output-type`, region/target/include/exclude selectors.

- `query`  
  Render VCF/BCF as user-defined text. Takes a bcftools `query` command line
  verbatim.  
  Example: `kira-bt query -f '%CHROM\t%POS\t%REF\t%ALT\n' in.vcf.gz`.

- `stats`  
  VCF/BCF statistics in the `bcftools stats` layout (SN, TSTV, SiS, AF, QUAL,
  IDD, ST, DP, PSC, PSI, HWE, VAF, and the two-file concordance sections), so
  `plot-vcfstats` and MultiQC read the output.  
  Args: `<inputs...>`, `-s/--samples`, `-S/--samples-file`, `-f/--apply-filters`,
  `-F/--fasta-ref`, `-d/--depth`, `--af-bins`, `--af-tag`, `-1/--1st-allele-only`,
  `-I/--split-by-ID`, `-c/--collapse`, `-E/--exons`, `-u/--user-tstv`,
  `-i/--include`, `-e/--exclude`, region/target selectors, `--threads`.

- `convert`  
  Convert to and from GEN/HAP/legend/sample and TSV formats.  
  Args: `[input] [-- <passthrough>...]`, `-o/--output`, `-O/--output-type`,
  `--gvcf2vcf`, `--tsv2vcf`, `-G/--gensample`, `-g/--gensample2vcf`,
  `--haplegendsample[2vcf]`, `-H/--hapsample`, `--hapsample2vcf`, `--haploid`,
  `-c/--columns`, `-f/--fasta-ref`, `--vcf-ids`, `--sex`, `--tag`, selectors.

- `consensus`  
  Apply variants to a reference to build a consensus sequence.  
  Args: `[input]`, `-f/--fasta-ref`, `-c/--chain`, `-H/--haplotype`,
  `-I/--iupac-codes`, `-a/--absent`, `--mark-del`, `--mark-ins`, `--mark-snv`,
  `-m/--mask`, `--mask-with`, `-M/--missing`, `-s/--samples`,
  `-S/--samples-file`, `-p/--prefix`, `-o/--output`, selectors.

- `reheader`  
  Replace the header or rename samples.  
  Args: `<input> [-- <passthrough>...]`, `-h/--header`, `-s/--samples`,
  `--samples-list`, `-f/--fai`, `-o/--output`, `--temp-prefix`, `--threads`.

### 5. Calling and Analysis

- `mpileup`  
  Multi-way pileup producing genotype likelihoods. `FORMAT/PL` is always
  emitted regardless of `-a/--annotate`; `-a` also accepts `AD`, `ADF`, `ADR`,
  `DP`, `QS`, `SP`, `SCR`, `GQ` and their `INFO/` forms. `##contig` lines carry
  the reference lengths from the BAM header, results do not depend on
  `--threads`, and CRAM inputs are decoded against `-f/--fasta-ref` (region
  queries need the `.crai`).  
  Args: `<inputs...>`, `-o/--output`, `-O/--output-type`,
  `-f/--fasta-ref`, `-b/--bam-list`, `-g/--gvcf`, `-a/--annotate`,
  `-d/--max-depth`, `-D/--max-idepth`, `-q/--min-MQ`, `-Q/--min-BQ`,
  `-m/--min-ireads`, `-I/--skip-indels`, `-B/--no-BAQ`, `-E/--redo-BAQ`,
  `-x/--ignore-overlaps`, `-X/--config`, `-F/--gap-frac`, `--indel-size`,
  `--threads`, `--stream`, flag filters (`--ef`, `--df`, `--if`, `--nf`,
  `--rf`, `--ff`), and selectors.
  Fused calling shortcuts: `-v/--variants-only`, `--prior`, `--min-alt-reads`,
  `--min-af`, `--min-qual`.

  The likelihood model follows htslib/bcftools: base qualities are capped by
  MAPQ and 60 (MAPQ 255 counts as 20); BAQ is the `probaln_glocal` pair-HMM
  applied as in `sam_prob_realn` (extended mode; `-E` recomputes it, `-B`
  skips it); overlapping mates are merged so a fragment counts once (`-x`
  disables); genotype likelihoods come from the dependent-error model of
  `errmod.c`, the site prior is a Watterson-scaled θ, and `INFO/QS` is the
  per-sample-normalised quality sum per allele. Kira extensions, all off by
  default and deviations from bcftools when enabled:
  `--nm-weight off|auto|FULL,SLOPE` (down-weight high-mismatch reads),
  `--hp-indel` (stricter indels in homopolymer/STR tracts),
  `--indel-realign off|ins|all` (pair-HMM realignment of indel candidates),
  `--assemble` (local assembly in active regions with haplotype PLs from the
  same pair-HMM kernel) and `--recal` (empirical base-quality recalibration
  by reported Q × trinucleotide × strand × cycle, learned at non-variant
  sites; needs `-f`, not available with `--stream`).
  `--nm-weight` and `--indel-realign` also read `KIRA_NM_WEIGHT` and
  `KIRA_INDEL_REALIGN` when the flag is absent, so scripts written against
  the pre-0.7 environment-variable interface keep working; an explicit flag
  always wins.

- `call`  
  SNP/indel calling from genotype likelihoods.  
  Args: `<input>`, `-o/--output`, `-O/--output-type`, `-c/--consensus-caller`,
  `-m/--multiallelic-caller` (default), `-v/--variants-only`, `-A/--keep-alts`,
  `-M/--keep-masked-ref`, `--keep-unseen-allele`, `-p/--pval-threshold`,
  `-P/--prior`, `-F/--prior-freqs <AN,AC>`, `--prior-af <INFO/TAG>`,
  `-X/--novel-rate`, `-C/--constrain <trio|alleles>`, `-g/--gvcf`,
  `-a/--annotate <GQ[,GP]>`, `--ploidy <0|1|2>`, `--ploidy-file`,
  `-G/--group-samples`, `--group-samples-tag`, `-s/--samples`,
  `-S/--samples-file`, `-V/--skip-variants`, selectors, `--threads`,
  `-W/--write-index`.

  The multiallelic caller reproduces `bcftools call -m`: allele frequencies
  come from `INFO/QS` (per-group `FORMAT/AD` fractions with `-G`, `-G -` for
  one group per sample), optionally shrunk towards `-F AN,AC` panel counts;
  `--prior-af INFO/TAG` does the same from an allele-frequency annotation,
  weighted like one diploid pseudo-sample per called sample. QUAL is the
  marginal `-4.343·(ref_lk − logsumexp(lk_sum, ref_lk))`, genotypes are the
  posterior mode under the HWE prior at those frequencies, GQ is
  phred(1 − P(best)) capped at 127, PL is dropped at reference-only sites,
  and `I16` is folded into `DP4`/`MQ`. `-C trio` with a PED given through
  `-S` is a kira extension (bcftools 1.19 has it disabled): the trio
  genotypes maximise the joint likelihood under a Mendelian transmission
  prior with the `-X` de-novo rate.

  `--ploidy` takes a plain number; region-dependent ploidy goes through
  `--ploidy-file` (`chrom beg end sex ploidy`, `*` allowed in the first three
  columns), with sample sex read from `-S/--samples-file` as either a
  `sample sex` table or a PED. Haploid samples are emitted as single-allele
  genotypes and contribute one allele each to `AC`/`AN`. The consensus caller
  (`-c`) is diploid-only and rejects other ploidies rather than silently
  calling diploid.

  ```bash
  kira-bt call pileup.vcf -o calls.vcf -m -v --ploidy 1
  kira-bt call pileup.vcf -o calls.vcf -m -v --ploidy-file ploidy.txt -S samples.txt
  ```

- `csq`  
  Haplotype-aware consequence caller.  
  Args: `[input]`, `-f/--fasta-ref`, `-g/--gff-annot`, `-l/--local-csq`,
  `-p/--phase`, `-c/--custom-tag`, `-B/--brief-predictions`, `-H/--haplotypes`,
  `--ncsq`, `-o/--output`, `-O/--output-type`, selectors, `--threads`.

- `roh`  
  Runs of homozygosity / autozygosity.  
  Args: `<input> [-- <passthrough>...]`, `--AF-tag`, `--AF-file`, `--AF-dflt`,
  `-G/--GTs-only`, `-I/--ignore-homref`, `--include-noalt`, `-X/--skip-indels`,
  `-V/--viterbi-training`, `-M/--rec-rate`, `-a/--hw-to-az`, `-H/--az-to-hw`,
  `-b/--buffer-size`, `-m/--genetic-map`, `-E/--estimate-AF`, `-o/--output`,
  `-O/--output-type`, selectors.

- `cnv`  
  Copy-number variation caller. Takes a bcftools `cnv` command line verbatim.

- `polysomy`  
  Detect contamination and whole-chromosome aberrations.  
  Args: `<input> [-- <bcftools args>...]`.

- `gtcheck`  
  Sample concordance, swap and contamination detection.  
  Args: `[input]`, `-g/--genotypes`, `-p/--pairs`, `-P/--pairs-file`,
  `-s/--samples`, `-S/--samples-file`, `-e/--error-probability`,
  `-n/--n-matches`, `-H/--homs-only`, `-u/--use`, `--no-HWE-prob`,
  `--keep-refs`, `--distinctive-sites`, `--dry-run`, `-o/--output`,
  `-O/--output-type`, selectors.

### 6. Alignment Files

- `samview`  
  BAM viewer and filter, the `samtools view -b` analogue.  
  Args: `<input>`, `-o/--output`, `-b/--bam`, `-C/--cram`, `-h/--with-header`,
  `-H/--header-only`, `-c/--count`, `-f/--require-flags`, `-F/--exclude-flags`,
  `-q/--min-MQ`, `-r/--region`, `-L/--regions-file`, `-s/--subsample`,
  `--threads`.


## Performance notes

- One thread budget per process: `--threads N` (or `-@ N`) caps every pool
  (BGZF compression and decompression, rayon, BAM decoding). Without it the
  budget is the machine's core count.
- BGZF output is compressed with libdeflate and BGZF input is decompressed
  with libdeflate (through noodles); plain gzip goes through zlib-rs. Only a
  C compiler is needed for the build, no cmake.
- `-r`/`-R` on an indexed VCF (`.csi`/`.tbi`) reads only the blocks that
  overlap the regions (`view`, `query`); every other selector filters through
  one sorted-interval set with binary search, so a 200k-exon BED costs a
  lookup per record, not a scan.
- `-f/--fasta-ref` is read through the `.fai` index one contig at a time
  (the index is built in memory when missing, and a `.fai` whose line layout
  does not match the file is ignored). Gzipped FASTA is loaded whole.
- `mpileup --stream` decodes each BAM on its own thread and pileups from the
  channel, so no BAM is held in memory; the reference pre-scan skip-list is
  not available on this path.
- `gtcheck -g` walks the panel in lock-step with the query; both inputs must
  be position-sorted with the same chromosome order.

## Tests

```bash
cargo test
```

Rust unit tests live beside the code (included from `tests/unit/`), with
end-to-end pipeline tests in `tests/*.rs`.

Alongside them are **404** bcftools comparison fixtures, grouped per command in
`tests/<mode>/testN`. Each fixture contains:
- command scripts (`bcftools.sh`, `kira.sh`)
- input assets
- reference outputs (`out.bcf.ref.vcf`, `out.kira.ref.vcf`)
- per-suite documentation in `tests/<mode>/README.md`

## License

Apache 2.0, see [LICENSE](LICENSE).
