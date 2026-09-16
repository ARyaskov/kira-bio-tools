//! GFF3 model for `csq`: a port of bcftools `gff.c`/`gff.h` and the parts of
//! `regidx.c` it relies on. Coordinates are 0-based inclusive and kept as
//! `u32` with wrapping arithmetic where the C code mixes signed and unsigned
//! operands, so out-of-range edge cases behave the same way.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Result, bail};
use fxhash::FxHashMap;

pub const STRAND_REV: u8 = 0;
pub const STRAND_FWD: u8 = 1;
pub const TRIM_5PRIME: u8 = 1;
pub const TRIM_3PRIME: u8 = 2;
pub const CDS_PHASE_UNKN: u8 = 3;

pub const N_SPLICE_DONOR: u32 = 2;
pub const N_SPLICE_REGION_EXON: u32 = 3;
pub const N_SPLICE_REGION_INTRON: u32 = 8;

const GF_CODING_BIT: u32 = 6;
#[inline]
pub fn gf_is_coding(t: u32) -> bool {
    t & (1 << GF_CODING_BIT) != 0
}
const fn coding(n: u32) -> u32 {
    n | (1 << GF_CODING_BIT)
}
const fn special(n: u32) -> u32 {
    (1 << (GF_CODING_BIT + 1)) + n
}

pub const GF_MT_RRNA: u32 = 1;
pub const GF_MT_TRNA: u32 = 2;
pub const GF_LINCRNA: u32 = 3;
pub const GF_MIRNA: u32 = 4;
pub const GF_MISC_RNA: u32 = 5;
pub const GF_RRNA: u32 = 6;
pub const GF_SNRNA: u32 = 7;
pub const GF_SNORNA: u32 = 8;
pub const GF_PROCESSED_TRANSCRIPT: u32 = 9;
pub const GF_ANTISENSE: u32 = 10;
pub const GF_MACRO_LNCRNA: u32 = 11;
pub const GF_RIBOZYME_LC: u32 = 12;
pub const GF_SRNA: u32 = 13;
pub const GF_SCRNA: u32 = 14;
pub const GF_SCARNA: u32 = 15;
pub const GF_SENSE_INTRONIC: u32 = 16;
pub const GF_SENSE_OVERLAPPING: u32 = 17;
pub const GF_PSEUDOGENE: u32 = 18;
pub const GF_PROCESSED_PSEUDOGENE: u32 = 19;
pub const GF_ARTIFACT: u32 = 20;
pub const GF_IG_PSEUDOGENE: u32 = 21;
pub const GF_IG_C_PSEUDOGENE: u32 = 22;
pub const GF_IG_J_PSEUDOGENE: u32 = 23;
pub const GF_IG_V_PSEUDOGENE: u32 = 24;
pub const GF_TR_V_PSEUDOGENE: u32 = 25;
pub const GF_TR_J_PSEUDOGENE: u32 = 26;
pub const GF_MT_TRNA_PSEUDOGENE: u32 = 27;
pub const GF_MISC_RNA_PSEUDOGENE: u32 = 28;
pub const GF_MIRNA_PSEUDOGENE: u32 = 29;
pub const GF_RIBOZYME: u32 = 30;
pub const GF_RETAINED_INTRON: u32 = 31;
pub const GF_RETROTRANSPOSED: u32 = 32;
pub const GF_TRNA_PSEUDOGENE: u32 = 33;
pub const GF_TRANSCRIBED_PROCESSED_PSEUDOGENE: u32 = 34;
pub const GF_TRANSCRIBED_UNPROCESSED_PSEUDOGENE: u32 = 35;
pub const GF_TRANSCRIBED_UNITARY_PSEUDOGENE: u32 = 36;
pub const GF_TRANSLATED_UNPROCESSED_PSEUDOGENE: u32 = 37;
pub const GF_TRANSLATED_PROCESSED_PSEUDOGENE: u32 = 38;
pub const GF_KNOWN_NCRNA: u32 = 39;
pub const GF_UNITARY_PSEUDOGENE: u32 = 40;
pub const GF_UNPROCESSED_PSEUDOGENE: u32 = 41;
pub const GF_LRG_GENE: u32 = 42;
pub const GF_3PRIME_OVERLAPPING_NCRNA: u32 = 43;
pub const GF_DISRUPTED_DOMAIN: u32 = 44;
pub const GF_VAULTRNA: u32 = 45;
pub const GF_BIDIRECTIONAL_PROMOTER_LNCRNA: u32 = 46;
pub const GF_AMBIGUOUS_ORF: u32 = 47;
pub const GF_LNCRNA: u32 = 48;

pub const GF_PROTEIN_CODING: u32 = coding(1);
pub const GF_POLYMORPHIC_PSEUDOGENE: u32 = coding(2);
pub const GF_IG_C: u32 = coding(3);
pub const GF_IG_D: u32 = coding(4);
pub const GF_IG_J: u32 = coding(5);
pub const GF_IG_LV: u32 = coding(6);
pub const GF_IG_V: u32 = coding(7);
pub const GF_TR_C: u32 = coding(8);
pub const GF_TR_D: u32 = coding(9);
pub const GF_TR_J: u32 = coding(10);
pub const GF_TR_V: u32 = coding(11);
pub const GF_NMD: u32 = coding(12);
pub const GF_NON_STOP_DECAY: u32 = coding(13);

pub const GF_CDS: u32 = special(1);
pub const GF_EXON: u32 = special(2);
pub const GF_UTR3: u32 = special(3);
pub const GF_UTR5: u32 = special(4);

const GF_STRINGS_NONCODING: [&str; 48] = [
    "MT_rRNA",
    "MT_tRNA",
    "lincRNA",
    "miRNA",
    "misc_RNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "processed_transcript",
    "antisense",
    "macro_lncRNA",
    "ribozyme",
    "sRNA",
    "scRNA",
    "scaRNA",
    "sense_intronic",
    "sense_overlapping",
    "pseudogene",
    "processed_pseudogene",
    "artifact",
    "IG_pseudogene",
    "IG_C_pseudogene",
    "IG_J_pseudogene",
    "IG_V_pseudogene",
    "TR_V_pseudogene",
    "TR_J_pseudogene",
    "MT_tRNA_pseudogene",
    "misc_RNA_pseudogene",
    "miRNA_pseudogene",
    "ribozyme",
    "retained_intron",
    "retrotransposed",
    "Trna_pseudogene",
    "transcribed_processed_pseudogene",
    "transcribed_unprocessed_pseudogene",
    "transcribed_unitary_pseudogene",
    "translated_unprocessed_pseudogene",
    "translated_processed_pseudogene",
    "known_ncRNA",
    "unitary_pseudogene",
    "unprocessed_pseudogene",
    "LRG_gene",
    "3_prime_overlapping_ncRNA",
    "disrupted_domain",
    "vaultRNA",
    "bidirectional_promoter_lncRNA",
    "ambiguous_orf",
    "lncRNA",
];
const GF_STRINGS_CODING: [&str; 13] = [
    "protein_coding",
    "polymorphic_pseudogene",
    "IG_C",
    "IG_D",
    "IG_J",
    "IG_LV",
    "IG_V",
    "TR_C",
    "TR_D",
    "TR_J",
    "TR_V",
    "NMD",
    "non_stop_decay",
];
const GF_STRINGS_SPECIAL: [&str; 4] = ["CDS", "exon", "3_prime_UTR", "5_prime_UTR"];

/// `gf_type2gff_string`.
pub fn gf_type2gff_string(ty: u32) -> &'static str {
    if !gf_is_coding(ty) {
        if ty < (1 << GF_CODING_BIT) {
            return GF_STRINGS_NONCODING
                .get(ty.wrapping_sub(1) as usize)
                .copied()
                .unwrap_or("");
        }
        let t = ty & ((1 << (GF_CODING_BIT + 1)) - 1);
        return GF_STRINGS_SPECIAL
            .get(t.wrapping_sub(1) as usize)
            .copied()
            .unwrap_or("");
    }
    let t = ty & ((1 << GF_CODING_BIT) - 1);
    GF_STRINGS_CODING
        .get(t.wrapping_sub(1) as usize)
        .copied()
        .unwrap_or("")
}

/// `strncmp(s, pat, n) == 0` on a C string: bytes past the end of `s` read as
/// NUL, so a pattern longer than the remaining text never matches.
fn strncmp_eq(s: &[u8], pat: &[u8], n: usize, fold: bool) -> bool {
    for i in 0..n {
        let a = s.get(i).copied().unwrap_or(0);
        let b = pat.get(i).copied().unwrap_or(0);
        let (a, b) = if fold {
            (a.to_ascii_lowercase(), b.to_ascii_lowercase())
        } else {
            (a, b)
        };
        if a != b {
            return false;
        }
        if a == 0 {
            return true;
        }
    }
    true
}

/// `gff_parse_biotype`: `-1` for a missing string, `0` when unknown. `s`
/// runs to the end of the line, exactly like the C pointer it mirrors.
pub fn parse_biotype(s: Option<&[u8]>) -> i32 {
    let Some(s) = s else { return -1 };
    let eq = |pat: &str| strncmp_eq(s, pat.as_bytes(), pat.len(), false);
    let ieq = |pat: &str, n: usize| strncmp_eq(s, pat.as_bytes(), n, true);
    let t = match s.first().copied().unwrap_or(0) {
        b'p' => {
            if eq("protein_coding") {
                GF_PROTEIN_CODING
            } else if eq("pseudogene") {
                GF_PSEUDOGENE
            } else if eq("processed_transcript") {
                GF_PROCESSED_TRANSCRIPT
            } else if eq("processed_pseudogene") {
                GF_PROCESSED_PSEUDOGENE
            } else if eq("polymorphic_pseudogene") {
                GF_POLYMORPHIC_PSEUDOGENE
            } else {
                0
            }
        }
        b'a' => {
            if eq("artifact") {
                GF_ARTIFACT
            } else if eq("antisense") {
                GF_ANTISENSE
            } else if eq("ambiguous_orf") {
                GF_AMBIGUOUS_ORF
            } else {
                0
            }
        }
        b'I' => {
            if eq("IG_pseudogene") {
                GF_IG_PSEUDOGENE
            } else if eq("IG_C_pseudogene") {
                GF_IG_C_PSEUDOGENE
            } else if eq("IG_J_pseudogene") {
                GF_IG_J_PSEUDOGENE
            } else if eq("IG_V_pseudogene") {
                GF_IG_V_PSEUDOGENE
            } else if eq("IG_C") {
                GF_IG_C
            } else if eq("IG_D") {
                GF_IG_D
            } else if eq("IG_J") {
                GF_IG_J
            } else if eq("IG_V") {
                GF_IG_V
            } else if eq("IG_LV") {
                GF_IG_LV
            } else {
                0
            }
        }
        b'T' => {
            if eq("TR_V_pseudogene") {
                GF_TR_V_PSEUDOGENE
            } else if eq("TR_J_pseudogene") {
                GF_TR_J_PSEUDOGENE
            } else if eq("TR_C") {
                GF_TR_C
            } else if eq("TR_D") {
                GF_TR_D
            } else if eq("TR_J") {
                GF_TR_J
            } else if eq("TR_V") {
                GF_TR_V
            } else {
                0
            }
        }
        b'M' => {
            if eq("Mt_tRNA_pseudogene") {
                GF_MT_TRNA_PSEUDOGENE
            } else if ieq("Mt_tRNA", 7) {
                GF_MT_TRNA
            } else if ieq("Mt_rRNA", 7) {
                GF_MT_TRNA
            } else if ieq("MRNA", 4) {
                GF_PROTEIN_CODING
            } else {
                0
            }
        }
        b'l' => {
            if eq("lincRNA") {
                GF_LINCRNA
            } else if strncmp_eq(s, b"lncRNA", 7, false) {
                GF_LNCRNA
            } else {
                0
            }
        }
        b'm' => {
            if eq("macro_lncRNA") {
                GF_MACRO_LNCRNA
            } else if eq("misc_RNA_pseudogene") {
                GF_MISC_RNA_PSEUDOGENE
            } else if eq("miRNA_pseudogene") {
                GF_MIRNA_PSEUDOGENE
            } else if eq("miRNA") {
                GF_MIRNA
            } else if eq("misc_RNA") {
                GF_MISC_RNA
            } else if ieq("mRNA", 4) {
                GF_PROTEIN_CODING
            } else {
                0
            }
        }
        b'r' => {
            if eq("rRNA") {
                GF_RRNA
            } else if eq("ribozyme") {
                GF_RIBOZYME
            } else if eq("retained_intron") {
                GF_RETAINED_INTRON
            } else if eq("retrotransposed") {
                GF_RETROTRANSPOSED
            } else {
                0
            }
        }
        b's' => {
            if eq("snRNA") {
                GF_SNRNA
            } else if eq("sRNA") {
                GF_SRNA
            } else if eq("scRNA") {
                GF_SCRNA
            } else if eq("scaRNA") {
                GF_SCARNA
            } else if eq("snoRNA") {
                GF_SNORNA
            } else if eq("sense_intronic") {
                GF_SENSE_INTRONIC
            } else if eq("sense_overlapping") {
                GF_SENSE_OVERLAPPING
            } else {
                0
            }
        }
        b't' => {
            if eq("tRNA_pseudogene") {
                GF_TRNA_PSEUDOGENE
            } else if eq("transcribed_processed_pseudogene") {
                GF_TRANSCRIBED_PROCESSED_PSEUDOGENE
            } else if eq("transcribed_unprocessed_pseudogene") {
                GF_TRANSCRIBED_UNPROCESSED_PSEUDOGENE
            } else if eq("transcribed_unitary_pseudogene") {
                GF_TRANSCRIBED_UNITARY_PSEUDOGENE
            } else if eq("translated_unprocessed_pseudogene") {
                GF_TRANSLATED_UNPROCESSED_PSEUDOGENE
            } else if eq("translated_processed_pseudogene") {
                GF_TRANSLATED_PROCESSED_PSEUDOGENE
            } else {
                0
            }
        }
        b'n' => {
            if eq("nonsense_mediated_decay") {
                GF_NMD
            } else if eq("non_stop_decay") {
                GF_NON_STOP_DECAY
            } else {
                0
            }
        }
        b'N' => {
            if eq("NMD") {
                GF_NMD
            } else {
                0
            }
        }
        b'k' => {
            if eq("known_ncrna") {
                GF_KNOWN_NCRNA
            } else {
                0
            }
        }
        b'u' => {
            if eq("unitary_pseudogene") {
                GF_UNITARY_PSEUDOGENE
            } else if eq("unprocessed_pseudogene") {
                GF_UNPROCESSED_PSEUDOGENE
            } else {
                0
            }
        }
        b'L' => {
            if eq("LRG_gene") {
                GF_LRG_GENE
            } else {
                0
            }
        }
        b'3' => {
            if ieq("3prime_overlapping_ncRNA", 24) {
                GF_3PRIME_OVERLAPPING_NCRNA
            } else if ieq("3_prime_overlapping_ncRNA", 25) {
                GF_3PRIME_OVERLAPPING_NCRNA
            } else {
                0
            }
        }
        b'd' => {
            if eq("disrupted_domain") {
                GF_DISRUPTED_DOMAIN
            } else {
                0
            }
        }
        b'v' => {
            if eq("vaultRNA") {
                GF_VAULTRNA
            } else {
                0
            }
        }
        b'b' => {
            if eq("bidirectional_promoter_lncRNA") {
                GF_BIDIRECTIONAL_PROMOTER_LNCRNA
            } else {
                0
            }
        }
        _ => 0,
    };
    t as i32
}

/// Interval index with `regidx` semantics: regions of one sequence sorted by
/// start, longer first; a query returns every overlapping region in that
/// order.
pub struct RegIdx<T: Copy> {
    seqs: FxHashMap<String, usize>,
    lists: Vec<RegList<T>>,
}

struct RegList<T: Copy> {
    regs: Vec<(u32, u32, T)>,
    /// Running maximum of `end` over the sorted list, for the first-hit search.
    max_end: Vec<u32>,
    sorted: bool,
}

impl<T: Copy> Default for RegIdx<T> {
    fn default() -> Self {
        Self {
            seqs: FxHashMap::default(),
            lists: Vec::new(),
        }
    }
}

impl<T: Copy> RegIdx<T> {
    pub fn push(&mut self, chr: &str, beg: u32, end: u32, payload: T) {
        let i = match self.seqs.get(chr) {
            Some(&i) => i,
            None => {
                self.seqs.insert(chr.to_string(), self.lists.len());
                self.lists.push(RegList {
                    regs: Vec::new(),
                    max_end: Vec::new(),
                    sorted: false,
                });
                self.lists.len() - 1
            }
        };
        let l = &mut self.lists[i];
        l.regs.push((beg, end, payload));
        l.sorted = false;
    }

    pub fn nregs(&self) -> usize {
        self.lists.iter().map(|l| l.regs.len()).sum()
    }

    fn sort(l: &mut RegList<T>) {
        l.regs.sort_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
        l.max_end.clear();
        let mut m = 0u32;
        for r in &l.regs {
            m = m.max(r.1);
            l.max_end.push(m);
        }
        l.sorted = true;
    }

    /// Sort every list; queries on an unsorted index sort lazily otherwise.
    pub fn finalize(&mut self) {
        for l in &mut self.lists {
            if !l.sorted {
                Self::sort(l);
            }
        }
    }

    /// Payloads of the regions overlapping `[beg, end]`.
    pub fn overlap(&self, chr: &str, beg: u32, end: u32) -> Vec<T> {
        let Some(&i) = self.seqs.get(chr) else {
            return Vec::new();
        };
        let l = &self.lists[i];
        debug_assert!(l.sorted, "RegIdx::finalize must run before queries");
        let first = l.max_end.partition_point(|&m| m < beg);
        let mut out = Vec::new();
        for r in &l.regs[first..] {
            if r.0 > end {
                break;
            }
            if r.1 >= beg && r.0 <= end {
                out.push(r.2);
            }
        }
        out
    }
}

#[derive(Clone, Debug)]
pub struct Cds {
    pub tr: usize,
    pub beg: u32,
    /// Offset of the first base within the spliced transcript.
    pub pos: u32,
    pub len: u32,
    pub icds: u32,
    pub phase: u8,
}

#[derive(Clone, Debug, Default)]
pub struct Gene {
    pub name: Option<String>,
    pub iseq: usize,
    pub beg: u32,
    pub end: u32,
    pub strand: u8,
    pub used: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Prime {
    P3,
    P5,
}

#[derive(Clone, Debug)]
pub struct Utr {
    pub which: Prime,
    pub beg: u32,
    pub end: u32,
    pub tr: usize,
}

#[derive(Clone, Debug)]
pub struct Exon {
    pub beg: u32,
    pub end: u32,
    pub tr: usize,
}

#[derive(Clone, Debug)]
pub struct Tscript {
    /// Index into [`Gff::tscript_ids`].
    pub id: u32,
    pub beg: u32,
    pub end: u32,
    pub strand: u8,
    pub used: bool,
    /// Indices into [`Gff::cds`], sorted by start after [`parse`].
    pub cds: Vec<usize>,
    pub trim: u8,
    pub ty: u32,
    pub gene: usize,
}

#[derive(Clone, Debug, Default)]
pub struct GffOpts {
    pub strip_chr_names: bool,
    pub force: bool,
    pub verbosity: u8,
}

#[derive(Default)]
pub struct Gff {
    pub tscript_ids: Vec<String>,
    pub genes: Vec<Gene>,
    pub tscripts: Vec<Tscript>,
    pub cds: Vec<Cds>,
    pub utrs: Vec<Utr>,
    pub exons: Vec<Exon>,
    /// Sequence names in order of first appearance (`aux->seq`).
    pub seqs: Vec<String>,
    pub idx_cds: RegIdx<usize>,
    pub idx_utr: RegIdx<usize>,
    pub idx_exon: RegIdx<usize>,
    pub idx_tscript: RegIdx<usize>,
}

impl Gff {
    pub fn tscript_name(&self, id: u32) -> &str {
        &self.tscript_ids[id as usize]
    }
    pub fn gene_name(&self, gene: usize) -> Option<&str> {
        self.genes[gene].name.as_deref()
    }
}

struct Ftr {
    ty: u32,
    beg: u32,
    end: u32,
    trid: u32,
    strand: u8,
    phase: u8,
}

#[derive(Default)]
struct Warned {
    unknown_tscript_biotype: u32,
    unknown_strand: u32,
    unknown_phase: u32,
    duplicate_id: u32,
    unknown_cds_phase: u32,
    incomplete_cds: u32,
    wrong_phase: u32,
    overlapping_cds: u32,
}

struct Parser<'a> {
    opts: &'a GffOpts,
    gff: Gff,
    ftrs: Vec<Ftr>,
    seq2int: FxHashMap<String, usize>,
    gene_ids: Vec<String>,
    gene_str2id: FxHashMap<String, u32>,
    gid2gene: HashMap<u32, usize>,
    tscript_str2id: FxHashMap<String, u32>,
    id2tr: HashMap<u32, usize>,
    /// Transcripts in registration order, for deterministic indexing.
    tr_order: Vec<u32>,
    ignored_biotypes: Vec<(String, u32)>,
    warned: Warned,
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    crate::util::open_text_maybe_gz(path)
}

fn warn_once(count: &mut u32, verbosity: u8, msg: impl FnOnce() -> String) {
    if verbosity > 0 && (*count == 0 || verbosity > 1) {
        eprintln!("{}", msg());
    }
    *count += 1;
}

impl<'a> Parser<'a> {
    fn feature_set_seq(&mut self, chr: &str) -> usize {
        if let Some(&i) = self.seq2int.get(chr) {
            return i;
        }
        let i = self.gff.seqs.len();
        self.gff.seqs.push(chr.to_string());
        self.seq2int.insert(chr.to_string(), i);
        i
    }

    fn register_id(tbl: &mut Vec<String>, map: &mut FxHashMap<String, u32>, s: &str) -> u32 {
        if let Some(&id) = map.get(s) {
            return id;
        }
        let id = tbl.len() as u32;
        tbl.push(s.to_string());
        map.insert(s.to_string(), id);
        id
    }

    fn gene_init(&mut self, gene_id: u32) -> usize {
        if let Some(&g) = self.gid2gene.get(&gene_id) {
            return g;
        }
        self.gff.genes.push(Gene::default());
        let g = self.gff.genes.len() - 1;
        self.gid2gene.insert(gene_id, g);
        g
    }

    /// `gff_parse_line`: returns the feature to keep for the second pass.
    fn parse_line(&mut self, line: &str) -> Result<Option<Ftr>> {
        let b = line.as_bytes();
        if b.is_empty() || b[0] == b'#' {
            return Ok(None);
        }
        let bad = || anyhow::anyhow!("Could not parse the line: {line}");
        let mut tabs = line.match_indices('\t').map(|(i, _)| i);
        let t1 = tabs.next().ok_or_else(bad)?;
        let t2 = tabs.next().ok_or_else(bad)?;
        let t3 = tabs.next().ok_or_else(bad)?;
        let t4 = tabs.next().ok_or_else(bad)?;
        let t5 = tabs.next().ok_or_else(bad)?;
        let t6 = tabs.next().ok_or_else(bad)?;
        let mut chr = &line[..t1];
        if self.opts.strip_chr_names && chr.len() >= 3 && chr[..3].eq_ignore_ascii_case("chr") {
            chr = &chr[3..];
        }
        let type_col = &line[t2 + 1..t3];
        let mut ty = 0u32;
        let mut is_gene_line = false;
        match type_col {
            "exon" => ty = GF_EXON,
            "CDS" => ty = GF_CDS,
            "three_prime_UTR" => ty = GF_UTR3,
            "five_prime_UTR" => ty = GF_UTR5,
            "biological_region" => return Ok(None),
            "gene" => is_gene_line = true,
            _ => {}
        }
        let beg = strtol_prefix(&line[t3 + 1..t4])
            .ok_or_else(bad)?
            .wrapping_sub(1) as u32;
        let end = strtol_prefix(&line[t4 + 1..t5])
            .ok_or_else(bad)?
            .wrapping_sub(1) as u32;
        let rest = &line[t6 + 1..];
        // Strand and phase are read as single characters, then skipped by two.
        let rb = rest.as_bytes();
        let mut strand = -1i32;
        match rb.first() {
            Some(b'+') => strand = STRAND_FWD as i32,
            Some(b'-') => strand = STRAND_REV as i32,
            _ => {}
        }
        let rest = rest.get(2..).unwrap_or("");
        let rb = rest.as_bytes();
        let mut phase = -1i32;
        match rb.first() {
            Some(b'0') => phase = 0,
            Some(b'1') => phase = 1,
            Some(b'2') => phase = 2,
            Some(b'.') => phase = CDS_PHASE_UNKN as i32,
            _ => {}
        }
        let attrs = rest.get(2..).unwrap_or("");
        let mut id: Option<&str> = None;
        let mut name: Option<&str> = None;
        let mut parent: Option<&str> = None;
        let mut biotype_rest: Option<&str> = None;
        let mut ss = attrs;
        loop {
            let es = ss.find(';').unwrap_or(ss.len());
            let tok = &ss[..es];
            if let Some(v) = tok.strip_prefix("ID=") {
                let v = if let Some(x) = v.strip_prefix("gene:") {
                    is_gene_line = true;
                    x
                } else if let Some(x) = v.strip_prefix("transcript:") {
                    x
                } else {
                    v
                };
                id = Some(v);
            } else if let Some(v) = tok.strip_prefix("Name=") {
                name = Some(v);
            } else if let Some(v) = tok.strip_prefix("Parent=") {
                let v = v
                    .strip_prefix("gene:")
                    .or_else(|| v.strip_prefix("transcript:"))
                    .unwrap_or(v);
                parent = Some(v);
            } else if let Some(v) = tok.strip_prefix("biotype=") {
                biotype_rest = Some(&ss[(es - v.len())..]);
            } else if let Some(v) = tok.strip_prefix("gene_biotype=") {
                biotype_rest = Some(&ss[(es - v.len())..]);
            }
            if es == ss.len() {
                break;
            }
            ss = &ss[es + 1..];
        }
        let strand_u = if strand == STRAND_REV as i32 {
            STRAND_REV
        } else {
            STRAND_FWD
        };
        let phase_u = if phase < 0 {
            CDS_PHASE_UNKN
        } else {
            phase as u8
        };

        if is_gene_line || parent.is_none() {
            // gff_parse_gene
            let Some(id) = id else { return Ok(None) };
            let gene_id = Self::register_id(&mut self.gene_ids, &mut self.gene_str2id, id);
            let g = self.gene_init(gene_id);
            if self.gff.genes[g].name.is_some() {
                let v = self.opts.verbosity;
                warn_once(&mut self.warned.duplicate_id, v, || {
                    format!("Warning: The GFF contains features with duplicate id .. {line}")
                });
                return Ok(None);
            }
            let iseq = self.feature_set_seq(chr);
            let gene = &mut self.gff.genes[g];
            gene.iseq = iseq;
            gene.beg = beg;
            gene.end = end;
            gene.strand = strand_u;
            gene.name = Some(
                name.map(str::to_string)
                    .unwrap_or_else(|| self.gene_ids[gene_id as usize].clone()),
            );
            return Ok(None);
        }
        if ty != 0 {
            // gff_parse_exon
            let Some(parent) = parent else {
                bail!(
                    "Could not parse the line, neither \"Parent=transcript:\" nor \"Parent=\" substring found: {line}"
                );
            };
            let trid =
                Self::register_id(&mut self.gff.tscript_ids, &mut self.tscript_str2id, parent);
            let v = self.opts.verbosity;
            if strand == -1 {
                warn_once(&mut self.warned.unknown_strand, v, || {
                    format!("Warning: Ignoring GFF feature with unknown strand .. {line}")
                });
            }
            if phase == -1 {
                warn_once(&mut self.warned.unknown_phase, v, || {
                    format!("Warning: Ignoring GFF feature with unknown phase .. {line}")
                });
            }
            self.feature_set_seq(chr);
            return Ok(Some(Ftr {
                ty,
                beg,
                end,
                trid,
                strand: strand_u,
                phase: phase_u,
            }));
        }
        // gff_parse_transcript
        let mut t = parse_biotype(biotype_rest.map(str::as_bytes));
        if t <= 0 {
            t = parse_biotype(Some(type_col.as_bytes()));
        }
        if t <= 0 {
            if let Some(bt) = biotype_rest {
                let key = bt.split(';').next().unwrap_or("").to_string();
                match self.ignored_biotypes.iter_mut().find(|(k, _)| *k == key) {
                    Some(e) => e.1 += 1,
                    None => self.ignored_biotypes.push((key, 1)),
                }
            } else {
                let v = self.opts.verbosity;
                warn_once(&mut self.warned.unknown_tscript_biotype, v, || {
                    format!("Warning: Ignoring transcript with unknown biotype .. {line}")
                });
            }
            return Ok(None);
        }
        let Some(id) = id else {
            bail!(
                "Could not parse the line, neither \"ID=transcript:\" nor \"ID=\" substring is present: {line}"
            );
        };
        let Some(parent) = parent else {
            bail!(
                "Could not parse the line, neither \"Parent=gene:\" nor \"Parent=\" substring is present: {line}"
            );
        };
        let trid = Self::register_id(&mut self.gff.tscript_ids, &mut self.tscript_str2id, id);
        let gene_id = Self::register_id(&mut self.gene_ids, &mut self.gene_str2id, parent);
        let gene = self.gene_init(gene_id);
        let tr = Tscript {
            id: trid,
            beg,
            end,
            strand: strand_u,
            used: false,
            cds: Vec::new(),
            trim: 0,
            ty: t as u32,
            gene,
        };
        match self.id2tr.get(&trid) {
            Some(&i) => self.gff.tscripts[i] = tr,
            None => {
                self.gff.tscripts.push(tr);
                self.id2tr.insert(trid, self.gff.tscripts.len() - 1);
                self.tr_order.push(trid);
            }
        }
        Ok(None)
    }

    fn chr_of_gene(&self, gene: usize) -> &str {
        &self.gff.seqs[self.gff.genes[gene].iseq]
    }

    fn register(&mut self, ftr: &Ftr) -> Result<()> {
        let Some(&ti) = self.id2tr.get(&ftr.trid) else {
            return Ok(());
        };
        let gene = self.gff.tscripts[ti].gene;
        self.gff.tscripts[ti].used = true;
        self.gff.genes[gene].used = true;
        let chr = self.chr_of_gene(gene).to_string();
        match ftr.ty {
            GF_CDS => {
                if self.gff.tscripts[ti].strand != ftr.strand {
                    bail!(
                        "Conflicting strand in transcript {} .. {} vs {}",
                        ftr.trid,
                        self.gff.tscripts[ti].strand,
                        ftr.strand
                    );
                }
                self.gff.cds.push(Cds {
                    tr: ti,
                    beg: ftr.beg,
                    pos: 0,
                    len: ftr.end.wrapping_sub(ftr.beg).wrapping_add(1),
                    icds: 0,
                    phase: ftr.phase,
                });
                let ci = self.gff.cds.len() - 1;
                self.gff.tscripts[ti].cds.push(ci);
            }
            GF_EXON => {
                self.gff.exons.push(Exon {
                    beg: ftr.beg,
                    end: ftr.end,
                    tr: ti,
                });
                let ei = self.gff.exons.len() - 1;
                self.gff.idx_exon.push(
                    &chr,
                    ftr.beg.wrapping_sub(N_SPLICE_REGION_INTRON),
                    ftr.end.wrapping_add(N_SPLICE_REGION_INTRON),
                    ei,
                );
            }
            GF_UTR3 | GF_UTR5 => {
                let which = if ftr.ty == GF_UTR3 {
                    Prime::P3
                } else {
                    Prime::P5
                };
                self.gff.utrs.push(Utr {
                    which,
                    beg: ftr.beg,
                    end: ftr.end,
                    tr: ti,
                });
                let ui = self.gff.utrs.len() - 1;
                self.gff.idx_utr.push(&chr, ftr.beg, ftr.end, ui);
            }
            _ => bail!("unexpected feature type {}", ftr.ty),
        }
        Ok(())
    }

    /// `tscript_init_cds`: sort CDS, verify and fix phase, trim incomplete
    /// ends, set offsets and index.
    fn init_cds(&mut self) -> Result<()> {
        let v = self.opts.verbosity;
        for k in 0..self.tr_order.len() {
            let ti = self.id2tr[&self.tr_order[k]];
            let (gene, beg, end) = {
                let t = &self.gff.tscripts[ti];
                (t.gene, t.beg, t.end)
            };
            let chr = self.chr_of_gene(gene).to_string();
            self.gff.idx_tscript.push(&chr, beg, end, ti);
            if self.gff.tscripts[ti].cds.is_empty() {
                continue;
            }
            let name = self.gff.tscript_ids[self.gff.tscripts[ti].id as usize].clone();
            let strand = self.gff.tscripts[ti].strand;
            let mut order = self.gff.tscripts[ti].cds.clone();
            let cds = &mut self.gff.cds;
            order.sort_by_key(|&c| cds[c].beg);
            self.gff.tscripts[ti].cds = order.clone();
            let ncds = order.len();
            let mut len: i64 = 0;
            let phase_check = |cds: &mut Vec<Cds>,
                               warned: &mut Warned,
                               i: usize,
                               len: &mut i64,
                               force: bool|
             -> Result<bool> {
                let c = &cds[order[i]];
                if c.phase == CDS_PHASE_UNKN {
                    warn_once(&mut warned.unknown_cds_phase, v, || {
                        format!(
                            "Warning: CDS with unknown phase, could not verify reading frame in transcript {name}"
                        )
                    });
                    *len += c.len as i64;
                    return Ok(true);
                }
                let phase = if c.phase != 0 { 3 - c.phase as i64 } else { 0 };
                if phase != *len % 3 {
                    if !force {
                        bail!(
                            "Error: GFF3 assumption failed for transcript {}, CDS={}: phase!=len%3 (phase={}, len={}). Use the --force option to proceed anyway (at your own risk).",
                            name,
                            c.beg.wrapping_add(1),
                            phase,
                            len
                        );
                    }
                    warn_once(&mut warned.wrong_phase, v, || {
                        format!(
                            "Warning: The GFF has inconsistent phase column in transcript {}, skipping. CDS pos={}: phase!=len%3 (phase={}, len={})",
                            name,
                            c.beg.wrapping_add(1),
                            phase,
                            len
                        )
                    });
                    return Ok(false);
                }
                *len += c.len as i64;
                Ok(true)
            };
            let mut ok = true;
            if strand == STRAND_FWD {
                let c0 = order[0];
                if cds[c0].phase != CDS_PHASE_UNKN {
                    if cds[c0].phase != 0 {
                        self.gff.tscripts[ti].trim |= TRIM_5PRIME;
                    }
                    let p = cds[c0].phase as u32;
                    cds[c0].beg = cds[c0].beg.wrapping_add(p);
                    cds[c0].len = cds[c0].len.wrapping_sub(p);
                    cds[c0].phase = 0;
                }
                for i in 0..ncds {
                    if !phase_check(cds, &mut self.warned, i, &mut len, self.opts.force)? {
                        ok = false;
                        break;
                    }
                }
            } else {
                let last = order[ncds - 1];
                if cds[last].phase != CDS_PHASE_UNKN {
                    let mut i = ncds as i64 - 1;
                    let mut phase = cds[last].phase as i64;
                    if phase != 0 {
                        self.gff.tscripts[ti].trim |= TRIM_5PRIME;
                    }
                    while i >= 0 && phase > cds[order[i as usize]].len as i64 {
                        phase -= cds[order[i as usize]].len as i64;
                        cds[order[i as usize]].phase = 0;
                        cds[order[i as usize]].len = 0;
                        i -= 1;
                    }
                    // `i` can only go negative with a corrupt phase; C would
                    // read before the array, so stop at the first exon.
                    let ci = order[i.max(0) as usize];
                    if cds[ci].phase != 0 {
                        warn_once(&mut self.warned.incomplete_cds, v, || {
                            format!(
                                "Note: truncated transcript {name} with incomplete CDS (this is very common)"
                            )
                        });
                    }
                    let p = cds[ci].phase as u32;
                    cds[ci].len = cds[ci].len.wrapping_sub(p);
                    cds[ci].phase = 0;
                }
                for i in (0..ncds).rev() {
                    if !phase_check(cds, &mut self.warned, i, &mut len, self.opts.force)? {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                continue;
            }
            let mut len: i64 = 0;
            for i in 0..ncds {
                cds[order[i]].icds = i as u32;
                len += cds[order[i]].len as i64;
                if i == 0 {
                    continue;
                }
                let a = &cds[order[i - 1]];
                let b = &cds[order[i]];
                if a.beg.wrapping_add(a.len).wrapping_sub(1) >= b.beg {
                    let msg = format!(
                        "Warning: GFF contains overlapping CDS {}, {}-{} and {}-{} (ribosomal slippage?)",
                        name,
                        a.beg.wrapping_add(1),
                        a.beg.wrapping_add(a.len),
                        b.beg.wrapping_add(1),
                        b.beg.wrapping_add(b.len)
                    );
                    warn_once(&mut self.warned.overlapping_cds, v, || msg);
                }
            }
            if len % 3 != 0 {
                warn_once(&mut self.warned.incomplete_cds, v, || {
                    format!(
                        "Note: truncated transcript {name} with incomplete CDS (this is very common)"
                    )
                });
                self.gff.tscripts[ti].trim |= TRIM_3PRIME;
                if strand == STRAND_FWD {
                    let mut i = ncds as i64 - 1;
                    while i >= 0 && len % 3 != 0 {
                        let c = &mut cds[order[i as usize]];
                        let dlen = if c.len as i64 >= len % 3 {
                            len % 3
                        } else {
                            c.len as i64
                        };
                        c.len -= dlen as u32;
                        len -= dlen;
                        i -= 1;
                    }
                } else {
                    let mut i = 0usize;
                    while i < ncds && len % 3 != 0 {
                        let c = &mut cds[order[i]];
                        let dlen = if c.len as i64 >= len % 3 {
                            len % 3
                        } else {
                            c.len as i64
                        };
                        c.len -= dlen as u32;
                        c.beg = c.beg.wrapping_add(dlen as u32);
                        len -= dlen;
                        i += 1;
                    }
                }
            }
            let mut len: u32 = 0;
            for &ci in &order {
                cds[ci].pos = len;
                len = len.wrapping_add(cds[ci].len);
                let (b, l) = (cds[ci].beg, cds[ci].len);
                self.gff
                    .idx_cds
                    .push(&chr, b, b.wrapping_add(l).wrapping_sub(1), ci);
            }
        }
        Ok(())
    }
}

/// `strtol` on the leading digits (optional sign), `None` when none.
fn strtol_prefix(s: &str) -> Option<i64> {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && b[i].is_ascii_whitespace() {
        i += 1;
    }
    let neg = i < b.len() && b[i] == b'-';
    if i < b.len() && (b[i] == b'-' || b[i] == b'+') {
        i += 1;
    }
    let start = i;
    let mut v: i64 = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        v = v.saturating_mul(10).saturating_add((b[i] - b'0') as i64);
        i += 1;
    }
    if i == start {
        return None;
    }
    Some(if neg { -v } else { v })
}

/// `gff_parse`: read the file, connect exons/CDS/UTRs to transcripts and
/// build the four interval indexes.
pub fn parse(path: &Path, opts: &GffOpts) -> Result<Gff> {
    if opts.verbosity > 0 {
        eprintln!("Parsing {} ...", path.display());
    }
    let mut p = Parser {
        opts,
        gff: Gff::default(),
        ftrs: Vec::new(),
        seq2int: FxHashMap::default(),
        gene_ids: Vec::new(),
        gene_str2id: FxHashMap::default(),
        gid2gene: HashMap::new(),
        tscript_str2id: FxHashMap::default(),
        id2tr: HashMap::new(),
        tr_order: Vec::new(),
        ignored_biotypes: Vec::new(),
        warned: Warned::default(),
    };
    let mut rd = open_text(path)?;
    let mut line = String::new();
    loop {
        line.clear();
        if rd.read_line(&mut line)? == 0 {
            break;
        }
        let l = line.trim_end_matches(['\n', '\r']);
        if let Some(f) = p.parse_line(l)? {
            p.ftrs.push(f);
        }
    }
    let ftrs = std::mem::take(&mut p.ftrs);
    for f in &ftrs {
        p.register(f)?;
    }
    p.init_cds()?;
    p.gff.idx_cds.finalize();
    p.gff.idx_utr.finalize();
    p.gff.idx_exon.finalize();
    p.gff.idx_tscript.finalize();
    if opts.verbosity > 0 {
        eprintln!(
            "Indexed {} transcripts, {} exons, {} CDSs, {} UTRs",
            p.gff.idx_tscript.nregs(),
            p.gff.idx_exon.nregs(),
            p.gff.idx_cds.nregs(),
            p.gff.idx_utr.nregs()
        );
        if !p.ignored_biotypes.is_empty() {
            eprintln!("Ignored the following biotypes:");
            for (bt, n) in &p.ignored_biotypes {
                let bt = if bt == "TCE" {
                    "TCE (\"To be Experimentally Confirmed\")"
                } else {
                    bt.as_str()
                };
                eprintln!("\t{n}x\t.. {bt}");
            }
        }
        let w = &p.warned;
        let mut nwarn = 0u32;
        for c in [
            w.unknown_tscript_biotype,
            w.unknown_strand,
            w.unknown_phase,
            w.duplicate_id,
            w.unknown_cds_phase,
            w.incomplete_cds,
            w.wrong_phase,
            w.overlapping_cds,
        ] {
            if c > 0 && opts.verbosity <= 1 {
                nwarn += c - 1;
            }
        }
        if nwarn > 0 {
            eprintln!(
                "Warning: {nwarn} warnings were supressed, run with `--verbose 2` to see them all"
            );
        }
    }
    if p.gff.idx_tscript.nregs() == 0 {
        bail!(
            "Error: No usable transcripts found, likely a failure to parse a non-standard GFF file. Please check if the misc/gff2gff\n       or misc/gff2gff.py script can fix the problem (both do different things). See also the man page for the description\n       of the expected format http://samtools.github.io/bcftools/bcftools-man.html#csq"
        );
    }
    Ok(p.gff)
}

#[cfg(test)]
#[path = "../../tests/unit/csq_gff.rs"]
mod tests;
