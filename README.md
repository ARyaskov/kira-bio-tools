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
- `KIRA_BT_BCF_DROP_UNDEFINED=1` drops FILTER/INFO/FORMAT tags the header does
  not declare when writing BCF (warning once per tag) instead of failing; BCF
  stores tags by header index, so an undeclared tag cannot be written and
  htslib would produce a file its own reader rejects.
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

Record text: output floats are rendered as htslib renders them (`kputd`,
i.e. `%g` with six significant digits, after the value is narrowed to
`float`), so `1e-1` prints as `0.1`, `60.00` as `60` and `1e6` as `1e+06`,
and a sample column shorter than FORMAT is padded with missing fields (`.`
becomes `.:.:.`), as bcftools prints it. Text sinks apply this to QUAL, to
Float-typed INFO/FORMAT values and to short sample columns only, leaving every
other byte as read; BCF decoding, `query` and `merge` use the same formatter.

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
  Column semantics follow bcftools: `-c INFO` / `-c FMT` transfer every tag
  of the annotation header one by one (target tags the source lacks or
  leaves missing are kept; Flag tags take the source record's state),
  `DST:=SRC` renames on transfer (`FMT/newGT:=GT`, `FMT/GT:=newGT`, with
  genotypes copied whole and `+FMT/GT` filling only `./.`),
  `INFO/TAG:=./FILTER` copies the target's own FILTER column (filters joined
  by `,`), `-c ALT` pairs records on POS/REF and replaces ALT, and a tag
  defined differently in the two headers keeps the input definition with a
  warning. Records are matched first as written, then with bcftools' REF
  padding rule (vcmp), then by their representation-free core (REF/ALT
  stripped of shared suffix and prefix, POS moved accordingly), so
  `GAT>GCT` at 100 finds `A>C` at 101 and a deletion anchored on either side
  finds the other; a core match is re-expressed over the target's REF before
  values are mapped. Left-alignment inside repeats needs a reference and is
  not attempted. Not yet mirrored: INFO/END matching for symbolic alleles,
  allele-aware FORMAT vectors under `-c FMT`, and bcftools' `+`-mode quirk
  that replaces an Integer vector whole when its first value is missing.  

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
  queries. Layouts: `-p vcf|bed|gff|sam` (guessed from the file name when
  `-p` is absent) or custom `-s/-b/-e` columns with `-0`, `-c`, `-S`; the
  index records the layout, so queries, `-h`/`-H` and `-l` follow it for any
  text file. Regions use tabix semantics (`chr:beg-end`, `-0`, `-R/-T` region
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
  Next to the `.kbi` it writes a `.kbm` position bitmap of the same pass:
  per contig a roaring-style set (65536-position chunks stored as sorted
  arrays or bitsets) of the record positions, plus the bases REF alleles
  span beyond POS. Membership, "any record in a range" and set algebra
  (and/or/and-not) run on the containers without touching the data file;
  `filter -M`, `isec` and `gtcheck -g` use it (see Performance notes).

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
  BCF in, BCF out (`-O b`/`-O u`) with only `-r/-R/-t/-T` (overlap modes 0
  and 1) copies the encoded records through without decoding them; any
  option that reads or rewrites a record takes the text path.

- `filter`  
  Expression filtering with soft-filter and gap controls. Without `-s`,
  `-S`, `--mask` or gap options the expression is evaluated on the record
  borrowed from its line and the line is written back unchanged (FILTER
  `.` becomes `PASS`, as bcftools does).  
  The engine follows bcftools for per-sample subscripts (`AD[:1]`,
  `FR[0:1]`), Flag tags compared as 0/1 (`INDEL=0`), aggregates against
  FORMAT fields (`sum(AD[*]) > FORMAT/DP`) and `-g` type lists
  (`-g2:mnp,indel,other`, spans from INFO/END); `-S .` with `-e` masks the
  samples the expression matches and recomputes AN/AC.  
  Args: `<input>`, `-i/--include`, `-e/--exclude`, `--expr`, `-s/--soft-filter`,
  `-m/--mode`, `-S/--set-GTs`, `-g/--SnpGap`, `-G/--IndelGap`, `--mask`,
  `-M/--mask-file`, `--mask-overlap`, `-o/--output`, `-O/--output-type`,
  `--threads`, `--gpu`, `-W/--write-index`.  
  `-M` takes a BED/VCF mask, a `.kbi` (record lines fetched by offset) or a
  `.kbm` position bitmap from `region-index`: the mask records' POS and REF
  spans are then answered by a set lookup, with the same result as the
  VCF mask for every `--mask-overlap` mode.

- `norm` (alias `N`)  
  Left-align, split/join multiallelics, check REF, remove duplicates.  
  Args: `<input>`, `-f/--fasta-ref`, `--fasta-ref-fai`, `-o/--output`,
  `-O/--output-type`, `-m/--multiallelics`, `--atomize`, `--atom-overlaps`,
  `-c/--check-ref`, `-N/--do-not-normalize`, `-d/--rm-dup`, `--site-win`,
  `--strict-filter`, `--keep-sum`, `--old-rec-tag`, `-S/--sort`, plus the usual
  region/target/include/exclude selectors.

- `sort`  
  Sort VCF/BCF by contig rank and position, ties broken as bcftools does
  (alleles compared one by one, case-insensitively, then input order).
  Records are kept as text lines, or as encoded BCF blocks when the input is
  BCF, so nothing is decoded unless the output type differs; runs beyond
  `-m` spill to BGZF-compressed temporary files and merge back k-way.  
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
  A port of bcftools' `consensus.c`: the FASTA is streamed line by line with
  the same buffer bookkeeping, output wraps at 60 columns, and sample
  selection is bcftools 1.19's (without `-s` every sample is used and
  genotypes become IUPAC codes; `-s -` applies every ALT; `-H` needs exactly
  one sample and takes `N`, `R`, `A`, `L[R|A]`, `S[R|A]`, `I` or `NpIu`).
  Masks are BED (0-based) or TAB (1-based) regions; a character mask also
  skips the variants it overlaps, `uc`/`lc` masks only change case.
  `--mask-with` pairs with `-m` by position on the command line.  
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
  `--rf`, `--ff`), selectors, `--bcftools-layout`, `-G/--read-groups`,
  `--ignore-RG`, `--max-read-len`, `-E/--redo-BAQ`, `--seed` (deep-site
  subsampling: one `hts_drand48` stream like bcftools in the bcftools layout,
  a per-site stream in the kira layout so results do not depend on
  `--threads`).
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
  `--hp-indel` (an STR-aware indel error model: each candidate's quality
  and required support follow the unit length and copy number of the repeat
  tract it sits in, a log-linear slippage model in the spirit of Dindel and
  Strelka2, instead of a fixed homopolymer penalty),
  `--indel-realign off|ins|all` (pair-HMM realignment of indel candidates),
  `--assemble` (local assembly in active regions: a read-threading de Bruijn
  graph over the spanning reads and the reference yields candidate
  haplotypes, a global alignment of each haplotype reads off its indels,
  and per-sample PLs follow a diploid haplotype model on read-vs-haplotype
  pair-HMM likelihoods, the best haplotype carrying an event against the
  best without it, so a SNP or second indel on the same haplotype no longer
  penalises a carrier read; the pair-HMM realignment vote remains the
  fallback when the graph finds nothing) and `--recal` (empirical base-quality recalibration
  by reported Q × trinucleotide × strand × cycle, learned at non-variant
  sites; `--recal-known <sites.ani|vcf>` also masks known variant sites, a
  `.ani` index being queried through its position index without decoding a
  record; needs `-f`, not available with `--stream`).
  `--nm-weight` and `--indel-realign` also read `KIRA_NM_WEIGHT` and
  `KIRA_INDEL_REALIGN` when the flag is absent, so scripts written against
  the pre-0.7 environment-variable interface keep working; an explicit flag
  always wins.

  `--bcftools-layout` writes records exactly the way `bcftools mpileup`
  does, so `bcftools call`, `kira-bt call` and the usual hard-filter recipes
  read them unchanged: alleles ordered by quality sum with the unseen `<*>`,
  PL only (no GT), `INFO/DP` as raw depth, `I16`, `QS`, `MQ0F` and the bias
  statistics (`VDB`, `SGB`, `RPBZ`, `MQBZ`, `BQBZ`, `MQSBZ`, `SCBZ`, `NMBZ`,
  `FS`) in bcftools' tag order, `-a` with bcftools' tag names and defaults
  (`-a DP,AD,INFO/AD,...`), bcftools' read and base filters (`-Q 1`,
  `--delta-BQ 30`, `--max-BQ 60`, `--ff UNMAP,SECONDARY,QCFAIL,DUP` as the
  only flag exclusion, orphan pairs dropped unless `-A`, `-d` per-file depth
  cap, `--max-read-len` for BAQ), lazy BAQ applied per site only where a
  pileup read carries an indel (bcftools 1.19's partial realignment), stored
  `BQ` tags honoured unless `-E`, and `-r` regions walked one at a time in
  the order given. Sample selection follows bcftools in both layouts:
  `-s ^A,B` excludes, `-S FILE` rows may rename (`OLD NEW`), `-G FILE` lists
  read groups (`ID [FILE] [SAMPLE]`, `^FILE` to exclude) and `--ignore-RG`
  names one sample per file. The indel records keep kira's candidate model
  behind bcftools' layout, and `-g`, `-v` and `--assemble` need the kira
  layout. On the bcftools test suite the SNP records match bcftools 1.19
  byte for byte on 18 of 28 cases; the rest differ in a few base-quality sums
  where BAQ posteriors disagree by one or two phred at read ends.

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
  Haplotype-aware consequence caller, a port of `bcftools csq`: the same GFF3
  model (Ensembl-style `gene`/`transcript`/`exon`/`CDS`/UTR lines, phase
  checks, incomplete-CDS trimming), per-sample haplotype trees combining
  phased variants that share a codon or a frameshift, `INFO/BCSQ` strings
  (`@pos` back references, `*` upstream stop, `&`-joined effects) and the
  `FORMAT/BCSQ` bitmask decoded by `bcftools query -f '%TBCSQ'`.  
  Args: `[input]`, `-f/--fasta-ref`, `-g/--gff-annot`, `-l/--local-csq`,
  `-p/--phase a|m|r|R|s` (default `r`), `-c/--custom-tag`, `-B/--trim-protein-seq INT`,
  `-n/--ncsq INT` (default 15), `--force`, `--unify-chr-names 0|1` (default 1),
  `-s/--samples` (`-` ignores genotypes), `-S/--samples-file`, `-i/-e`,
  `-r/-R/--regions-overlap`, `-t/-T/--targets-overlap`, `-o/--output`,
  `-O/--output-type b|u|z|v|t` (`t` is the tab-delimited `CSQ` listing),
  `--no-version`, `--write-index`, `--threads`, `-v/--verbose 0-2`.  
  Output matches bcftools 1.19 byte for byte on `tests/csq`, including its
  quirks (stale `FORMAT/BCSQ` zeros, tie order of transcripts ending at the
  same position). One deliberate difference: with a `-s` subset bcftools
  writes the bits of UTR/splice/intron consequences into the columns of the
  first samples in the file; kira writes them into the selected samples'
  own columns.

- `roh`  
  Runs of autozygosity, a port of `bcftools roh` on the same two-state HMM
  (`HW`/`AZ`, precomputed transition powers, genetic-map or constant
  recombination rate, Viterbi + forward-backward with bcftools' quality
  bookkeeping, `-b` buffering with snapshots, `-V` Viterbi training).
  Output is bcftools' `# RG`/`# ST` tables (`-Oz` bgzips them) and matches
  bcftools 1.19 byte for byte on `tests/roh`.  
  Args: `[input] [-- <bcftools options>...]`, `--AF-dflt`, `--AF-tag`,
  `--AF-file`, `-b/--buffer-size`, `-e/--estimate-AF [GT,|PL,]FILE|-`,
  `-G/--GTs-only`, `-i/--ignore-homref`, `-I/--skip-indels`,
  `--include-noalt`, `-m/--genetic-map`, `-M/--rec-rate`, `-a/--hw-to-az`,
  `-H/--az-to-hw`, `-V/--viterbi-training`, `--include/--exclude`,
  `-o/--output`, `-O/--output-type s|r|z`, `-s/-S`, `-r/-R/--regions-overlap`,
  `-t/-T/--targets-overlap`, `--threads`. The short letters follow bcftools
  (`-e` estimates AF, `-i` ignores hom-ref, `-I` skips indels).

- `cnv`  
  Copy-number variation from Illumina BAF/LRR, a port of `bcftools cnv`:
  CN0..CN3 states (16 paired states with `-c`), per-genotype Gaussian BAF
  mixtures, LRR smoothing, `-O` aberrant-cell-fraction estimation and the
  `dat.*.tab`/`cn.*.tab`/`summary*.tab` outputs (plus the plotting
  scripts). Matches bcftools 1.19 byte for byte on `tests/cnv`.  
  Args: `[-- ]<bcftools cnv options> FILE.vcf`: `-s/--query-sample`,
  `-c/--control-sample`, `-f/--AF-file`, `-o/--output-dir`,
  `-p/--plot-threshold`, `-a/--aberrant`, `-b/--BAF-weight`, `-d/--BAF-dev`,
  `-e/--err-prob`, `-k/--LRR-dev`, `-l/--LRR-weight`, `-L/--LRR-smooth-win`,
  `-O/--optimize`, `-P/--same-prob`, `-x/--xy-prob`, `-r/-R/--regions-overlap`,
  `-t/-T/--targets-overlap`.

- `polysomy`  
  Detect contamination and whole-chromosome aberrations.  
  Args: `<input> [-- <bcftools args>...]`.

- `gtcheck`  
  Sample concordance, swap and contamination detection.  
  `--sketch K` (all-vs-all only) first sketches every sample's dosage vector
  with `--sketch-bits` (256) SimHash projections in one pass, then scores
  exactly only the K nearest samples of each sample: duplicate and relative
  search in a large cohort costs O(n·m·bits) instead of O(n²·m).  
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


### 7. Columnar Export

- `zarr`
  Export a VCF/BCF to a VCF Zarr store (specification 0.2, the layout
  sgkit and `vcf2zarr` read) as a Zarr v2 directory: `variant_contig`,
  `variant_position`, `variant_id`, `variant_allele`, `variant_quality`,
  `variant_filter`, `call_genotype` with `call_genotype_phased` and
  `call_genotype_mask`, the `contig_id`/`contig_length`/`filter_id`/
  `sample_id` tables, and every Number=1 Integer/Float INFO/FORMAT field
  plus INFO flags as `variant_<TAG>`/`call_<TAG>` (or the `--fields`
  selection). Arrays are chunked along variants and samples and
  zlib-compressed (numcodecs `zlib`); strings are fixed-width `S`; the
  header is kept in the root attributes. Two passes over the input: the
  first sizes the arrays (variants, alleles, ploidy, string widths), the
  second writes the chunks.
  Args: `<input>`, `-o/--output`, `--chunk-variants`, `--chunk-samples`,
  `-l/--level`, `--fields`, `-f/--force`.

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
  not available on this path. The in-memory path hands each parallel chunk
  reads that share their sequence buffers (copied only when mate overlap
  adjusts qualities), BAQ reuses its pair-HMM matrices per thread, and sample
  columns are rendered straight into the output buffer. The pair-HMM keeps
  its three states in separate planes so the match and insertion cells of a
  row are computed in one vectorisable loop (bit-identical to htslib's
  probaln), and identical reads (same placement, bases and qualities) reuse
  one BAQ result through a per-thread cache (`KIRA_BT_TIMING=1` reports the
  hit count).
- `gtcheck -g` walks the panel in lock-step with the query; both inputs must
  be position-sorted with the same chromosome order. With a `.kbm` next to
  the panel, query sites the panel has no record at are counted as
  unmatched without parsing their genotypes.
- `isec` with `.kbm` bitmaps next to every input a record must be present in
  (`-n=N`/`+N` equal to the number of files, or the `1`s of `-n ~MASK`)
  intersects the bitmaps first and drops records outside the intersection
  before parsing them.

## Tests

```bash
cargo test
```

Rust unit tests live beside the code (included from `tests/unit/`), with
end-to-end pipeline tests in `tests/*.rs`. `tests/unit/props.rs` holds
proptest property tests for the parsers and formatters that face untrusted
input (float rendering, record canonicalisation, BCF round trips, tabix and
genotype parsing never panic). The crate builds on stable Rust
(`rust-version` in `Cargo.toml`); CI runs `cargo clippy --all-targets -D
warnings` under the default lint set, with the few structural lints the code
does not follow listed at the top of `src/lib.rs`.

Alongside them are **404** bcftools comparison fixtures, grouped per command in
`tests/<mode>/testN`. Each fixture contains:
- command scripts (`bcftools.sh`, `kira.sh`)
- input assets
- reference outputs (`out.bcf.ref.vcf`, `out.kira.ref.vcf`)
- per-suite documentation in `tests/<mode>/README.md`

## License

Apache 2.0, see [LICENSE](LICENSE).
