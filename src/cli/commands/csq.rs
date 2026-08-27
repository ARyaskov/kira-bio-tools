use anyhow::{Result, anyhow};
use std::collections::HashMap;
use std::fs;

use crate::VcfReader;
use crate::cli::args::CsqArgs;

pub fn cmd_csq(args: CsqArgs) -> Result<()> {
    let mut argv: Vec<String> = Vec::new();
    if let Some(p) = &args.input { argv.push(p.to_string_lossy().into_owned()); }
    if let Some(p) = &args.fasta_ref { argv.push("-f".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(p) = &args.gff { argv.push("-g".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(s) = &args.samples { argv.push("-s".into()); argv.push(s.clone()); }
    if let Some(p) = &args.samples_file { argv.push("-S".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(s) = &args.regions { argv.push("-r".into()); argv.push(s.clone()); }
    if let Some(p) = &args.regions_file { argv.push("-R".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(p) = &args.output { argv.push("-o".into()); argv.push(p.to_string_lossy().into_owned()); }
    if let Some(s) = &args.include { argv.push("-i".into()); argv.push(s.clone()); }
    if let Some(s) = &args.exclude { argv.push("-e".into()); argv.push(s.clone()); }
    let cfg = parse_args(&argv)?;
    let gff_path = cfg
        .gff
        .as_ref()
        .ok_or_else(|| anyhow!("missing -g <ann.gff> for csq"))?;
    let input = cfg
        .input
        .as_ref()
        .ok_or_else(|| anyhow!("missing input VCF for csq"))?;

    let anno = AnnotationDb::from_gff(gff_path)?;
    let fasta = if let Some(f) = &cfg.fasta {
        Some(FastaDb::from_path(f)?)
    } else {
        None
    };

    let mut reader = VcfReader::open(input)?;
    let mut headers = reader.header()?;
    headers.retain(|h| !h.starts_with("##INFO=<ID=BCSQ,"));
    headers.retain(|h| !h.starts_with("##FORMAT=<ID=BCSQ,"));

    let mut inserted = false;
    for h in &headers {
        if !inserted && h.starts_with("#CHROM\t") {
            println!(
                "##INFO=<ID=BCSQ,Number=.,Type=String,Description=\"Consequence annotation compatible with bcftools csq\">"
            );
            println!(
                "##FORMAT=<ID=BCSQ,Number=.,Type=Integer,Description=\"Per-sample csq indexes\">"
            );
            inserted = true;
        }
        println!("{h}");
    }

    // Buffer per chromosome so the haplotype-aware pass can combine phased
    // variants that share a codon before any record is emitted.
    let mut buf: Vec<crate::vcf::structs::VcfRecord> = Vec::new();
    let mut cur_chrom = String::new();
    while let Some(rec) = reader.next_record()? {
        if !buf.is_empty() && rec.chrom != cur_chrom {
            emit_chunk(&buf, &anno, fasta.as_ref());
            buf.clear();
        }
        cur_chrom = rec.chrom.clone();
        buf.push(rec);
    }
    if !buf.is_empty() {
        emit_chunk(&buf, &anno, fasta.as_ref());
    }

    Ok(())
}

/// Haplotype combination for one record: replaces the coding consequence of the
/// combined transcript with the merged annotation (on the leader) or a back
/// reference `@<leader pos>` (on the follower).
struct Combo {
    leader: bool,
    leader_pos: u32,
    effect: String,
    aa: String,
    dna: String,
    gene: String,
    tx: String,
    biotype: String,
    strand: char,
}

/// Annotate and print one chromosome's worth of records, merging phased
/// same-codon variants into haplotype-aware consequences (bcftools `csq`).
fn emit_chunk(
    buf: &[crate::vcf::structs::VcfRecord],
    anno: &AnnotationDb,
    fasta: Option<&FastaDb>,
) {
    let combos = build_combinations(buf, anno, fasta);
    for (ri, rec0) in buf.iter().enumerate() {
        let mut rec = rec0.clone();
        let alts = rec
            .alt
            .split(',')
            .map(str::trim)
            .filter(|a| !a.is_empty() && *a != ".")
            .map(|a| a.to_string())
            .collect::<Vec<_>>();

        let mut consequences = Vec::<String>::new();
        let mut alt_to_csq_idxs = Vec::<Vec<usize>>::new();
        for alt in &alts {
            let mut per_alt = annotate_variant(&rec, alt, anno, fasta);
            if let Some(combo) = combos.get(&ri) {
                apply_combo(&mut per_alt, combo);
            }
            let mut idxs = Vec::<usize>::with_capacity(per_alt.len());
            for c in per_alt.drain(..) {
                consequences.push(c);
                idxs.push(consequences.len());
            }
            alt_to_csq_idxs.push(idxs);
        }

        if consequences.is_empty() {
            consequences.push(".".to_string());
        }
        let bcsq_val = consequences.join(",");

        let mut info = parse_info_map(&rec.info);
        info.insert("BCSQ".to_string(), bcsq_val);
        rec.info = render_info_map(&info);

        if let Some(fmt) = rec.format.clone() {
            let keys = fmt.split(':').collect::<Vec<_>>();
            let gt_idx = keys.iter().position(|k| *k == "GT");
            rec.format = Some(format!("{fmt}:BCSQ"));
            for s in &mut rec.samples {
                let parts = s.split(':').collect::<Vec<_>>();
                let gt = gt_idx.and_then(|i| parts.get(i).copied()).unwrap_or("./.");
                let bcsq = bcsq_for_gt(gt, &alt_to_csq_idxs);
                s.push(':');
                s.push_str(&bcsq);
            }
        }

        print_record(&rec);
    }
}

/// Coding effects carry strand|aa|dna fields and are the elements the
/// haplotype merge replaces.
fn is_coding_effect(effect: &str) -> bool {
    matches!(
        effect,
        "missense"
            | "synonymous"
            | "stop_gained"
            | "stop_lost"
            | "start_lost"
            | "inframe_deletion"
            | "inframe_insertion"
            | "frameshift"
            | "complex_substitution"
    )
}

/// True if a consequence string names a coding effect on transcript `tx` —
/// i.e. the element the haplotype merge should replace.
fn is_coding_for_tx(s: &str, tx: &str) -> bool {
    let f: Vec<&str> = s.split('|').collect();
    f.len() > 2 && f[2] == tx && is_coding_effect(f[0])
}

fn apply_combo(per_alt: &mut Vec<String>, combo: &Combo) {
    if combo.leader {
        let combined = format!(
            "{}|{}|{}|{}|{}|{}|{}",
            combo.effect, combo.gene, combo.tx, combo.biotype, combo.strand, combo.aa, combo.dna
        );
        match per_alt.iter().position(|s| is_coding_for_tx(s, &combo.tx)) {
            Some(i) => per_alt[i] = combined,
            None => per_alt.push(combined),
        }
    } else {
        per_alt.retain(|s| !is_coding_for_tx(s, &combo.tx));
        per_alt.push(format!("@{}", combo.leader_pos));
    }
}

fn csq_is_snv(rec: &crate::vcf::structs::VcfRecord) -> bool {
    rec.ref_allele.len() == 1 && {
        let a = rec.alt.trim();
        a.len() == 1 && a.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    }
}

fn csq_ref_char(rec: &crate::vcf::structs::VcfRecord) -> char {
    rec.ref_allele.chars().next().unwrap_or('N').to_ascii_uppercase()
}

fn csq_alt_char(rec: &crate::vcf::structs::VcfRecord) -> char {
    rec.alt.trim().chars().next().unwrap_or('N').to_ascii_uppercase()
}

/// Two records are on the same haplotype if some sample carries the ALT (allele
/// "1") at the same phased haplotype index in both.
fn phased_together(a: &crate::vcf::structs::VcfRecord, b: &crate::vcf::structs::VcfRecord) -> bool {
    let (Some(fa), Some(fb)) = (a.format.as_ref(), b.format.as_ref()) else {
        return false;
    };
    let (Some(gi_a), Some(gi_b)) = (
        fa.split(':').position(|k| k == "GT"),
        fb.split(':').position(|k| k == "GT"),
    ) else {
        return false;
    };
    let n = a.samples.len().min(b.samples.len());
    for si in 0..n {
        let ga = a.samples[si].split(':').nth(gi_a).unwrap_or(".");
        let gb = b.samples[si].split(':').nth(gi_b).unwrap_or(".");
        if !ga.contains('|') || !gb.contains('|') {
            continue;
        }
        let ha: Vec<&str> = ga.split('|').collect();
        let hb: Vec<&str> = gb.split('|').collect();
        for h in 0..ha.len().min(hb.len()) {
            if ha[h] == "1" && hb[h] == "1" {
                return true;
            }
        }
    }
    false
}

/// Find phased SNV pairs that share a codon in the same protein-coding
/// transcript and build the merged haplotype consequence for each.
fn build_combinations(
    buf: &[crate::vcf::structs::VcfRecord],
    anno: &AnnotationDb,
    fasta: Option<&FastaDb>,
) -> HashMap<usize, Combo> {
    let mut out = HashMap::new();
    let Some(fasta) = fasta else {
        return out;
    };
    // Group coding-SNV record indices by (transcript id, codon index).
    let mut groups: HashMap<(String, usize), Vec<usize>> = HashMap::new();
    for (ri, rec) in buf.iter().enumerate() {
        if !csq_is_snv(rec) {
            continue;
        }
        let Some(txs) = anno.by_chrom.get(&rec.chrom) else {
            continue;
        };
        for tx in txs {
            if tx.biotype != "protein_coding" {
                continue;
            }
            if !tx.cds.iter().any(|r| rec.pos >= r.start && rec.pos <= r.end) {
                continue;
            }
            let g2c = build_genome_to_cds_map(tx);
            if let Some(&cds) = g2c.get(&rec.pos) {
                groups
                    .entry((tx.id.clone(), (cds / 3) as usize))
                    .or_default()
                    .push(ri);
            }
        }
    }

    for ((tx_id, _codon), mut idxs) in groups {
        if idxs.len() < 2 {
            continue;
        }
        idxs.sort_by_key(|&i| buf[i].pos);
        let (lo, hi) = (idxs[0], idxs[1]);
        if buf[lo].pos == buf[hi].pos || !phased_together(&buf[lo], &buf[hi]) {
            continue;
        }
        let Some(txs) = anno.by_chrom.get(&buf[lo].chrom) else {
            continue;
        };
        let Some(tx) = txs.iter().find(|t| t.id == tx_id) else {
            continue;
        };
        let var_lo = (buf[lo].pos, csq_ref_char(&buf[lo]), csq_alt_char(&buf[lo]));
        let var_hi = (buf[hi].pos, csq_ref_char(&buf[hi]), csq_alt_char(&buf[hi]));
        let Some((aapos, aa_ref, aa_alt)) = combine_phased_codon(var_lo, var_hi, tx, fasta) else {
            continue;
        };
        let effect = if aa_ref == aa_alt {
            "synonymous"
        } else if aa_alt == '*' {
            "stop_gained"
        } else {
            "missense"
        };
        let aa = format!("{}{}>{}{}", aapos, aa_ref, aapos, aa_alt);
        // DNA change is always genomic-order; the leader is the variant first in
        // transcription order (higher genomic position on the minus strand).
        let dna = format!(
            "{}{}>{}+{}{}>{}",
            var_lo.0, var_lo.1, var_lo.2, var_hi.0, var_hi.1, var_hi.2
        );
        let (leader_i, follower_i) = if tx.strand == '-' { (hi, lo) } else { (lo, hi) };
        let leader_pos = buf[leader_i].pos;
        out.insert(
            leader_i,
            Combo {
                leader: true,
                leader_pos,
                effect: effect.to_string(),
                aa,
                dna,
                gene: tx.gene_name.clone(),
                tx: tx.id.clone(),
                biotype: tx.biotype.clone(),
                strand: tx.strand,
            },
        );
        out.insert(
            follower_i,
            Combo {
                leader: false,
                leader_pos,
                effect: String::new(),
                aa: String::new(),
                dna: String::new(),
                gene: tx.gene_name.clone(),
                tx: tx.id.clone(),
                biotype: tx.biotype.clone(),
                strand: tx.strand,
            },
        );
    }
    out
}

#[derive(Default)]
struct CsqCfg {
    fasta: Option<String>,
    gff: Option<String>,
    input: Option<String>,
}

fn parse_args(args: &[String]) -> Result<CsqCfg> {
    let mut cfg = CsqCfg::default();
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-f" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.fasta = Some(v.clone());
                }
            }
            "-g" => {
                i += 1;
                if let Some(v) = args.get(i) {
                    cfg.gff = Some(v.clone());
                }
            }
            "-p" | "--ncsq" | "--unify-chr-names" => {
                i += 1;
            }
            a => {
                if !a.starts_with('-') || a == "-" {
                    cfg.input = Some(a.to_string());
                }
            }
        }
        i += 1;
    }
    Ok(cfg)
}

#[derive(Clone)]
struct Region {
    start: u32,
    end: u32,
    /// GFF3 CDS phase (column 8): bases to skip to reach the first codon.
    /// 0 for non-CDS features.
    phase: u8,
}

#[derive(Clone)]
struct Transcript {
    id: String,
    chrom: String,
    start: u32,
    end: u32,
    strand: char,
    gene_id: String,
    gene_name: String,
    biotype: String,
    exons: Vec<Region>,
    cds: Vec<Region>,
    utr5: Vec<Region>,
    utr3: Vec<Region>,
}

#[derive(Default)]
struct AnnotationDb {
    by_chrom: HashMap<String, Vec<Transcript>>,
}

#[derive(Clone)]
struct GeneMeta {
    name: String,
    biotype: String,
}

#[derive(Default)]
struct TxBuild {
    tx: Option<Transcript>,
    parent: String,
}

impl AnnotationDb {
    fn from_gff(path: &str) -> Result<Self> {
        let text = fs::read_to_string(path)?;
        let mut genes = HashMap::<String, GeneMeta>::new();
        let mut tx_map = HashMap::<String, TxBuild>::new();

        for line in text.lines() {
            let l = line.trim();
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            let cols = l.split('\t').collect::<Vec<_>>();
            if cols.len() < 9 {
                continue;
            }
            let chrom = cols[0].to_string();
            let feature = cols[2];
            let start = cols[3].parse::<u32>().unwrap_or(0);
            let end = cols[4].parse::<u32>().unwrap_or(0);
            let strand = cols[6].chars().next().unwrap_or('+');
            let attrs = parse_attrs(cols[8]);

            if feature == "gene" || feature == "lincRNA_gene" {
                if let Some(id_raw) = attrs.get("ID") {
                    let id = clean_id(id_raw);
                    let name = attrs.get("Name").cloned().unwrap_or_else(|| id.clone());
                    let biotype = attrs
                        .get("biotype")
                        .or_else(|| attrs.get("gene_biotype"))
                        .cloned()
                        .unwrap_or_else(|| {
                            if feature == "lincRNA_gene" {
                                "lincRNA".to_string()
                            } else {
                                "protein_coding".to_string()
                            }
                        });
                    genes.insert(id, GeneMeta { name, biotype });
                }
                continue;
            }

            if feature == "transcript" || feature == "lincRNA" {
                let tx_id = attrs
                    .get("ID")
                    .map(|x| clean_id(x))
                    .unwrap_or_else(|| format!("tx_{chrom}_{start}_{end}"));
                let parent = attrs.get("Parent").map(|x| clean_id(x)).unwrap_or_default();

                let gm = genes.get(&parent).cloned().unwrap_or(GeneMeta {
                    name: parent.clone(),
                    biotype: attrs
                        .get("biotype")
                        .or_else(|| attrs.get("gene_biotype"))
                        .cloned()
                        .unwrap_or_else(|| {
                            if feature == "lincRNA" {
                                "lincRNA".to_string()
                            } else {
                                "protein_coding".to_string()
                            }
                        }),
                });

                tx_map.insert(
                    tx_id.clone(),
                    TxBuild {
                        tx: Some(Transcript {
                            id: tx_id,
                            chrom,
                            start,
                            end,
                            strand,
                            gene_id: parent.clone(),
                            gene_name: gm.name,
                            biotype: gm.biotype,
                            exons: Vec::new(),
                            cds: Vec::new(),
                            utr5: Vec::new(),
                            utr3: Vec::new(),
                        }),
                        parent,
                    },
                );
                continue;
            }

            if feature == "exon"
                || feature == "CDS"
                || feature == "five_prime_UTR"
                || feature == "three_prime_UTR"
            {
                let parent = attrs.get("Parent").map(|x| clean_id(x)).unwrap_or_default();
                let build = tx_map.entry(parent.clone()).or_default();
                if build.tx.is_none() {
                    let gm = genes.get(&build.parent).cloned().unwrap_or(GeneMeta {
                        name: build.parent.clone(),
                        biotype: "protein_coding".to_string(),
                    });
                    build.tx = Some(Transcript {
                        id: parent.clone(),
                        chrom: chrom.clone(),
                        start,
                        end,
                        strand,
                        gene_id: build.parent.clone(),
                        gene_name: gm.name,
                        biotype: gm.biotype,
                        exons: Vec::new(),
                        cds: Vec::new(),
                        utr5: Vec::new(),
                        utr3: Vec::new(),
                    });
                }
                let tx = build.tx.as_mut().expect("tx exists");
                if start < tx.start {
                    tx.start = start;
                }
                if end > tx.end {
                    tx.end = end;
                }
                let phase = if feature == "CDS" {
                    cols[7].parse::<u8>().unwrap_or(0)
                } else {
                    0
                };
                let r = Region { start, end, phase };
                match feature {
                    "exon" => tx.exons.push(r),
                    "CDS" => tx.cds.push(r),
                    "five_prime_UTR" => tx.utr5.push(r),
                    "three_prime_UTR" => tx.utr3.push(r),
                    _ => {}
                }
            }
        }

        let mut by_chrom = HashMap::<String, Vec<Transcript>>::new();
        for (_, build) in tx_map {
            if let Some(mut tx) = build.tx {
                tx.exons.sort_by_key(|r| r.start);
                tx.cds.sort_by_key(|r| r.start);
                tx.utr5.sort_by_key(|r| r.start);
                tx.utr3.sort_by_key(|r| r.start);
                by_chrom.entry(tx.chrom.clone()).or_default().push(tx);
            }
        }
        for v in by_chrom.values_mut() {
            v.sort_by(|a, b| a.id.cmp(&b.id));
        }

        Ok(Self { by_chrom })
    }
}

fn parse_attrs(s: &str) -> HashMap<String, String> {
    let mut out = HashMap::<String, String>::new();
    for part in s.split(';') {
        if let Some((k, v)) = part.split_once('=') {
            out.insert(k.trim().to_string(), v.trim().to_string());
        }
    }
    out
}

fn clean_id(v: &str) -> String {
    v.rsplit(':').next().unwrap_or(v).to_string()
}

struct FastaDb {
    seqs: HashMap<String, String>,
}

impl FastaDb {
    fn from_path(path: &str) -> Result<Self> {
        let text = fs::read_to_string(path)?;
        let mut seqs = HashMap::<String, String>::new();
        let mut name = String::new();
        let mut seq = String::new();
        for line in text.lines() {
            if let Some(h) = line.strip_prefix('>') {
                if !name.is_empty() {
                    seqs.insert(name.clone(), seq.clone());
                    seq.clear();
                }
                name = h.split_whitespace().next().unwrap_or(h).to_string();
            } else if !line.trim().is_empty() {
                seq.push_str(&line.trim().to_ascii_uppercase());
            }
        }
        if !name.is_empty() {
            seqs.insert(name, seq);
        }
        Ok(Self { seqs })
    }

    fn base_at(&self, chrom: &str, pos1: u32) -> Option<char> {
        let idx = pos1.checked_sub(1)? as usize;
        self.seqs
            .get(chrom)
            .and_then(|s| s.as_bytes().get(idx).copied())
            .map(|b| b as char)
    }

    fn slice(&self, chrom: &str, start1: u32, end1: u32) -> Option<String> {
        let s = self.seqs.get(chrom)?;
        let a = start1.checked_sub(1)? as usize;
        let b = end1 as usize;
        if a >= s.len() || b > s.len() || a >= b {
            return None;
        }
        Some(s[a..b].to_string())
    }
}

fn annotate_variant(
    rec: &crate::vcf::structs::VcfRecord,
    alt: &str,
    anno: &AnnotationDb,
    fasta: Option<&FastaDb>,
) -> Vec<String> {
    let start = rec.pos;
    let end = rec.pos + rec.ref_allele.len().saturating_sub(1) as u32;

    let Some(txs) = anno.by_chrom.get(&rec.chrom) else {
        return vec![format!(
            "intergenic|.|.|.|.|{}{}>{}|{}{}>{}",
            rec.pos, rec.ref_allele, alt, rec.pos, rec.ref_allele, alt
        )];
    };

    let mut out = Vec::<String>::new();
    for tx in txs {
        if end < tx.start || start > tx.end {
            continue;
        }
        let effect = infer_effect(rec, alt, tx, fasta);
        // bcftools prints strand|aa|dna only for coding consequences; UTR /
        // intron / splice / non-coding get just effect|gene|tx|biotype.
        let consequence = if is_coding_effect(&effect) && tx.biotype == "protein_coding" {
            let aa_field = coding_aa_change(rec, alt, tx, fasta)
                .map(|(ar, aa)| {
                    let cds_pos = build_genome_to_cds_map(tx).get(&rec.pos).copied().unwrap_or(0);
                    let aa_pos = (cds_pos / 3) + 1;
                    if ar == aa {
                        format!("{}{}", aa_pos, ar)
                    } else {
                        format!("{}{}>{}{}", aa_pos, ar, aa_pos, aa)
                    }
                })
                .unwrap_or_else(|| ".".to_string());
            let dna_change = format!("{}{}>{}", rec.pos, rec.ref_allele, alt);
            format!(
                "{}|{}|{}|{}|{}|{}|{}",
                effect, tx.gene_name, tx.id, tx.biotype, tx.strand, aa_field, dna_change
            )
        } else {
            format!("{}|{}|{}|{}", effect, tx.gene_name, tx.id, tx.biotype)
        };
        out.push(consequence);
    }

    if out.is_empty() {
        out.push(format!(
            "intergenic|.|.|.|.|.|{}{}>{}",
            rec.pos, rec.ref_allele, alt
        ));
    }

    out
}

/// Combine two adjacent phased SNVs that fall within the same codon
/// into a single haplotype-aware aa change. Returns Some(combined_aa)
/// when both variants on the same haplotype hit the same codon.
fn combine_phased_codon(
    var_a: (u32, char, char),
    var_b: (u32, char, char),
    tx: &Transcript,
    fasta: &FastaDb,
) -> Option<(u32, char, char)> {
    let cds = build_cds_sequence(tx, fasta)?;
    let g2c = build_genome_to_cds_map(tx);
    let a_cds = *g2c.get(&var_a.0)? as usize;
    let b_cds = *g2c.get(&var_b.0)? as usize;
    let codon_a = a_cds / 3;
    let codon_b = b_cds / 3;
    if codon_a != codon_b { return None; }
    let codon_start = codon_a * 3;
    if codon_start + 3 > cds.len() { return None; }
    let mut codon = cds[codon_start..codon_start + 3].as_bytes().to_vec();
    let aa_ref = translate_codon(&cds[codon_start..codon_start + 3]);
    codon[a_cds % 3] = if tx.strand == '-' { complement(var_a.2) as u8 } else { var_a.2 as u8 };
    codon[b_cds % 3] = if tx.strand == '-' { complement(var_b.2) as u8 } else { var_b.2 as u8 };
    let new_aa = translate_codon(std::str::from_utf8(&codon).ok()?);
    Some((codon_a as u32 + 1, aa_ref, new_aa))
}

fn infer_effect(
    rec: &crate::vcf::structs::VcfRecord,
    alt: &str,
    tx: &Transcript,
    fasta: Option<&FastaDb>,
) -> String {
    if tx.biotype != "protein_coding" {
        return "non_coding".to_string();
    }

    let start = rec.pos;
    let end = rec.pos + rec.ref_allele.len().saturating_sub(1) as u32;

    let overlaps_cds = tx.cds.iter().any(|r| overlap(start, end, r.start, r.end));
    if overlaps_cds {
        if rec.ref_allele.len() == alt.len() {
            if rec.ref_allele.len() == 1 {
                if let Some((aa_ref, aa_alt)) = coding_aa_change(rec, alt, tx, fasta) {
                    if aa_ref == aa_alt {
                        return "synonymous".to_string();
                    }
                    if aa_alt == '*' {
                        return "stop_gained".to_string();
                    }
                    return "missense".to_string();
                }
                return "missense".to_string();
            }
            return "complex_substitution".to_string();
        }

        let diff = alt.len() as i32 - rec.ref_allele.len() as i32;
        if diff % 3 == 0 {
            if diff < 0 {
                return "inframe_deletion".to_string();
            }
            if diff > 0 {
                return "inframe_insertion".to_string();
            }
            return "complex_substitution".to_string();
        }
        return "frameshift".to_string();
    }

    if tx.utr5.iter().any(|r| overlap(start, end, r.start, r.end)) {
        return "5_prime_utr".to_string();
    }
    if tx.utr3.iter().any(|r| overlap(start, end, r.start, r.end)) {
        return "3_prime_utr".to_string();
    }

    if is_splice_site(start, tx) {
        return splice_kind(start, tx).to_string();
    }

    if tx.exons.iter().any(|r| overlap(start, end, r.start, r.end)) {
        return "splice_region".to_string();
    }

    if start >= tx.start && start <= tx.end {
        return "intron".to_string();
    }

    "intergenic".to_string()
}

fn overlap(a1: u32, a2: u32, b1: u32, b2: u32) -> bool {
    !(a2 < b1 || b2 < a1)
}

fn is_splice_site(pos: u32, tx: &Transcript) -> bool {
    for ex in &tx.exons {
        if pos.abs_diff(ex.start) <= 1 || pos.abs_diff(ex.end) <= 1 {
            return true;
        }
    }
    false
}

fn splice_kind(pos: u32, tx: &Transcript) -> &'static str {
    for ex in &tx.exons {
        if tx.strand == '+' {
            if pos.abs_diff(ex.end) <= 1 {
                return "splice_donor";
            }
            if pos.abs_diff(ex.start) <= 1 {
                return "splice_acceptor";
            }
        } else {
            if pos.abs_diff(ex.start) <= 1 {
                return "splice_donor";
            }
            if pos.abs_diff(ex.end) <= 1 {
                return "splice_acceptor";
            }
        }
    }
    "splice_region"
}

fn coding_aa_change(
    rec: &crate::vcf::structs::VcfRecord,
    alt: &str,
    tx: &Transcript,
    fasta: Option<&FastaDb>,
) -> Option<(char, char)> {
    let fasta = fasta?;
    let cds = build_cds_sequence(tx, fasta)?;
    let g2c = build_genome_to_cds_map(tx);
    let cds_idx = *g2c.get(&rec.pos)? as usize;

    if cds_idx >= cds.len() {
        return None;
    }

    let codon_start = (cds_idx / 3) * 3;
    if codon_start + 3 > cds.len() {
        return None;
    }
    let mut codon = cds[codon_start..codon_start + 3].as_bytes().to_vec();
    let alt_base = alt.chars().next()?.to_ascii_uppercase();
    let mut ref_base = rec.ref_allele.chars().next()?.to_ascii_uppercase();
    if tx.strand == '-' {
        ref_base = complement(ref_base);
        codon[cds_idx % 3] = complement(alt_base) as u8;
    } else {
        codon[cds_idx % 3] = alt_base as u8;
    }

    let ref_codon = cds[codon_start..codon_start + 3].to_string();
    let alt_codon = String::from_utf8(codon).ok()?;
    let aa_ref = translate_codon(&ref_codon);
    let aa_alt = translate_codon(&alt_codon);
    let _ = ref_base;
    Some((aa_ref, aa_alt))
}

/// Start phase of a transcript: the GFF phase of its first CDS segment in
/// transcription order (lowest start on `+`, highest end on `-`). This many
/// bases are dropped from the 5' end before translation begins.
fn cds_start_phase(tx: &Transcript) -> usize {
    let first = if tx.strand == '-' {
        tx.cds.iter().max_by_key(|r| r.end)
    } else {
        tx.cds.iter().min_by_key(|r| r.start)
    };
    first.map(|r| r.phase as usize).unwrap_or(0)
}

fn build_cds_sequence(tx: &Transcript, fasta: &FastaDb) -> Option<String> {
    let mut parts = Vec::<String>::new();
    let mut cds = tx.cds.clone();
    cds.sort_by_key(|r| r.start);
    for r in &cds {
        parts.push(fasta.slice(&tx.chrom, r.start, r.end)?);
    }
    let mut seq = parts.join("");
    if tx.strand == '-' {
        seq = revcomp(&seq);
    }
    // Drop the leading partial codon indicated by the start phase.
    let p = cds_start_phase(tx);
    if p >= seq.len() {
        return Some(String::new());
    }
    Some(seq[p..].to_string())
}

fn build_genome_to_cds_map(tx: &Transcript) -> HashMap<u32, u32> {
    let mut out = HashMap::<u32, u32>::new();
    let mut cds = tx.cds.clone();
    cds.sort_by_key(|r| r.start);
    // Start the index at -phase so the first translated codon base maps to 0;
    // the dropped leading bases get negative indices and are skipped.
    let mut idx: i64 = -(cds_start_phase(tx) as i64);

    let push = |out: &mut HashMap<u32, u32>, pos: u32, idx: &mut i64| {
        if *idx >= 0 {
            out.insert(pos, *idx as u32);
        }
        *idx += 1;
    };

    if tx.strand == '+' {
        for r in cds {
            for p in r.start..=r.end {
                push(&mut out, p, &mut idx);
            }
        }
    } else {
        cds.reverse();
        for r in cds {
            let mut p = r.end;
            loop {
                push(&mut out, p, &mut idx);
                if p == r.start {
                    break;
                }
                p -= 1;
            }
        }
    }

    out
}

fn translate_codon(c: &str) -> char {
    match c {
        "TTT" | "TTC" => 'F',
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
        "ATT" | "ATC" | "ATA" => 'I',
        "ATG" => 'M',
        "GTT" | "GTC" | "GTA" | "GTG" => 'V',
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
        "CCT" | "CCC" | "CCA" | "CCG" => 'P',
        "ACT" | "ACC" | "ACA" | "ACG" => 'T',
        "GCT" | "GCC" | "GCA" | "GCG" => 'A',
        "TAT" | "TAC" => 'Y',
        "CAT" | "CAC" => 'H',
        "CAA" | "CAG" => 'Q',
        "AAT" | "AAC" => 'N',
        "AAA" | "AAG" => 'K',
        "GAT" | "GAC" => 'D',
        "GAA" | "GAG" => 'E',
        "TGT" | "TGC" => 'C',
        "TGG" => 'W',
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
        "GGT" | "GGC" | "GGA" | "GGG" => 'G',
        "TAA" | "TAG" | "TGA" => '*',
        _ => 'X',
    }
}

fn complement(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => 'N',
    }
}

fn revcomp(s: &str) -> String {
    s.chars().rev().map(complement).collect()
}

fn bcsq_for_gt(gt: &str, alt_to_csq_idxs: &[Vec<usize>]) -> String {
    let mut out = Vec::<usize>::new();
    for a in gt.split(['/', '|']) {
        if a == "." {
            continue;
        }
        let Ok(v) = a.parse::<usize>() else {
            continue;
        };
        if v == 0 {
            continue;
        }
        let alt_idx = v - 1;
        let Some(csq_idxs) = alt_to_csq_idxs.get(alt_idx) else {
            continue;
        };
        for idx in csq_idxs {
            if !out.contains(idx) {
                out.push(*idx);
            }
        }
    }
    if out.is_empty() {
        ".".to_string()
    } else {
        out.into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn parse_info_map(info: &str) -> HashMap<String, String> {
    let mut out = HashMap::<String, String>::new();
    if info == "." || info.is_empty() {
        return out;
    }
    for kv in info.split(';') {
        if kv.is_empty() {
            continue;
        }
        if let Some((k, v)) = kv.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        } else {
            out.insert(kv.to_string(), String::new());
        }
    }
    out
}

fn render_info_map(info: &HashMap<String, String>) -> String {
    if info.is_empty() {
        return ".".to_string();
    }
    let mut keys = info.keys().cloned().collect::<Vec<_>>();
    keys.sort();
    keys.into_iter()
        .map(|k| {
            let v = info.get(&k).map(|s| s.as_str()).unwrap_or("");
            if v.is_empty() { k } else { format!("{k}={v}") }
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn print_record(rec: &crate::vcf::structs::VcfRecord) {
    print!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        rec.chrom, rec.pos, rec.id, rec.ref_allele, rec.alt, rec.qual, rec.filter, rec.info
    );
    if let Some(fmt) = &rec.format {
        print!("\t{fmt}");
        for s in &rec.samples {
            print!("\t{s}");
        }
    }
    println!();
}
