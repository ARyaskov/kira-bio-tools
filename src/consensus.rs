//! Port of bcftools 1.19 `consensus.c`: apply VCF variants to a FASTA read
//! line by line, with the same buffer bookkeeping (`fa_ori_pos`,
//! `fa_mod_off`, `fa_frz_pos`), sample/haplotype selection, IUPAC codes,
//! masks, `--mark-*`, `-a`/`-M` characters and chain output.

use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};

pub const PICK_REF: u32 = 1;
pub const PICK_ALT: u32 = 2;
pub const PICK_LONG: u32 = 4;
pub const PICK_SHORT: u32 = 8;
pub const PICK_IUPAC: u32 = 16;

const VCF_REF: u32 = 0;
const VCF_SNP: u32 = 1;
const VCF_MNP: u32 = 2;
const VCF_INDEL: u32 = 4;
const VCF_OTHER: u32 = 8;
const VCF_BND: u32 = 16;
const VCF_OVERLAP: u32 = 32;

/// `--mark-ins`/`--mark-snv`: `uc`, `lc` or one character.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Mark {
    #[default]
    None,
    Upper,
    Lower,
    Char(u8),
}

impl Mark {
    pub fn parse(s: &str, opt: &str) -> Result<Self> {
        if s.eq_ignore_ascii_case("uc") {
            return Ok(Self::Upper);
        }
        if s.eq_ignore_ascii_case("lc") {
            return Ok(Self::Lower);
        }
        let b = s.as_bytes();
        if b.len() == 1 && b[0] > 32 && b[0] < 127 {
            return Ok(Self::Char(b[0]));
        }
        bail!("The argument is not recognised: --{opt} {s}")
    }
}

/// `--mask-with`: change case, or replace with a character (which also
/// skips variants overlapping the mask).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaskWith {
    Upper,
    Lower,
    Char(u8),
}

impl MaskWith {
    pub fn parse(s: &str) -> Result<Self> {
        if s.eq_ignore_ascii_case("uc") {
            return Ok(Self::Upper);
        }
        if s.eq_ignore_ascii_case("lc") {
            return Ok(Self::Lower);
        }
        if s.len() != 1 {
            bail!("Expected \"lc\", \"uc\", or a single character with the --mask-with option");
        }
        Ok(Self::Char(s.as_bytes()[0]))
    }
}

/// Regions of one `-m` file (regidx semantics: `.bed` is 0-based half-open,
/// anything else 1-based inclusive), stored 0-based inclusive.
#[derive(Debug)]
pub struct Mask {
    pub with: MaskWith,
    regions: HashMap<String, Vec<(i64, i64)>>,
}

impl Mask {
    pub fn load(path: &Path, with: MaskWith) -> Result<Self> {
        let name = path.to_string_lossy().to_ascii_lowercase();
        let is_bed = name.ends_with(".bed") || name.ends_with(".bed.gz");
        let reader = crate::util::open_text_maybe_gz(path)?;
        let mut regions: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for line in reader.lines() {
            let line = line?;
            let line = line.trim_end_matches('\r');
            if line.is_empty() || line.starts_with('#') || (is_bed && (line.starts_with("track") || line.starts_with("browser"))) {
                continue;
            }
            let mut it = line.split(['\t', ' ']).filter(|s| !s.is_empty());
            let Some(chr) = it.next() else { continue };
            let Some(b) = it.next().and_then(|s| s.parse::<i64>().ok()) else {
                bail!("Could not parse the region: {line}")
            };
            let e = it.next().and_then(|s| s.parse::<i64>().ok());
            let (beg, end) = if is_bed {
                let Some(e) = e else { bail!("Could not parse the region: {line}") };
                (b, e - 1)
            } else {
                (b - 1, e.unwrap_or(b) - 1)
            };
            regions.entry(chr.to_string()).or_default().push((beg, end));
        }
        for v in regions.values_mut() {
            v.sort_unstable();
        }
        Ok(Self { with, regions })
    }

    /// A character mask skips the variants it overlaps; case masks do not.
    fn skips(&self) -> bool {
        matches!(self.with, MaskWith::Char(_))
    }

    fn overlapping(&self, chr: &str, beg: i64, end: i64) -> impl Iterator<Item = (i64, i64)> + '_ {
        self.regions
            .get(chr)
            .into_iter()
            .flat_map(|v| v.iter().copied())
            .filter(move |&(rb, re)| rb <= end && re >= beg)
    }

    fn overlaps(&self, chr: &str, beg: i64, end: i64) -> bool {
        self.overlapping(chr, beg, end).next().is_some()
    }
}

/// One GT value as htslib stores it (`bcf_gt_*`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GtVal {
    Allele { idx: i32, phased: bool },
    Missing { phased: bool },
    /// `bcf_int8_vector_end` padding of a shorter genotype.
    End,
}

impl GtVal {
    fn is_missing(self) -> bool {
        matches!(self, Self::Missing { .. })
    }

    /// `bcf_gt_is_phased`; vector_end (0x81) has the bit set.
    fn is_phased(self) -> bool {
        match self {
            Self::Allele { phased, .. } | Self::Missing { phased } => phased,
            Self::End => true,
        }
    }

    /// `bcf_gt_allele`: -1 for missing.
    fn allele(self) -> i32 {
        match self {
            Self::Allele { idx, .. } => idx,
            Self::Missing { .. } => -1,
            Self::End => 63,
        }
    }
}

/// Parse a text genotype (`0/1`, `.|1`, `.`); phase marks the allele it precedes.
pub fn parse_gt(s: &str) -> Vec<GtVal> {
    let mut out = Vec::with_capacity(2);
    let mut phased = s.starts_with('|');
    for tok in s.trim_start_matches('|').split(['/', '|']) {
        let v = if tok == "." || tok.is_empty() {
            GtVal::Missing { phased }
        } else {
            match tok.parse::<i32>() {
                Ok(i) => GtVal::Allele { idx: i, phased },
                Err(_) => GtVal::Missing { phased },
            }
        };
        out.push(v);
        phased = false;
    }
    let bytes = s.as_bytes();
    let mut k = 1usize;
    let mut i = 0usize;
    while i < bytes.len() && k < out.len() {
        if bytes[i] == b'|' || bytes[i] == b'/' {
            if let GtVal::Allele { phased, .. } | GtVal::Missing { phased } = &mut out[k] {
                *phased = bytes[i] == b'|';
            }
            k += 1;
        }
        i += 1;
    }
    out
}

/// One VCF record as consensus needs it.
#[derive(Clone, Debug)]
pub struct Rec {
    /// 0-based.
    pub pos: i64,
    /// htslib `rlen`: INFO/END when present, else the REF length.
    pub rlen: i64,
    pub has_end: bool,
    pub alleles: Vec<Vec<u8>>,
    /// GT per sample, padded to the record's ploidy with [`GtVal::End`];
    /// `None` when the record has no GT field.
    pub gts: Option<Vec<Vec<GtVal>>>,
}

impl Rec {
    pub fn from_record(rec: &crate::vcf::VcfRecord) -> Self {
        let mut alleles: Vec<Vec<u8>> = vec![rec.ref_allele.as_bytes().to_vec()];
        if rec.alt != "." {
            alleles.extend(rec.alt.split(',').map(|a| a.as_bytes().to_vec()));
        }
        let end = rec
            .info
            .split(';')
            .find_map(|kv| kv.strip_prefix("END=").and_then(|v| v.parse::<i64>().ok()));
        let pos = rec.pos as i64 - 1;
        let (rlen, has_end) = match end {
            Some(e) if e > pos => (e - pos, true),
            _ => (rec.ref_allele.len() as i64, false),
        };
        let gts = rec.format.as_deref().and_then(|f| f.split(':').position(|k| k == "GT")).map(|gi| {
            let mut gts: Vec<Vec<GtVal>> = rec
                .samples
                .iter()
                .map(|s| parse_gt(s.split(':').nth(gi).unwrap_or(".")))
                .collect();
            let n = gts.iter().map(Vec::len).max().unwrap_or(0);
            for g in &mut gts {
                g.resize(n, GtVal::End);
            }
            gts
        });
        Self { pos, rlen, has_end, alleles, gts }
    }
}

/// `bcf_set_variant_type`: the type of `alt` against `ref` and its length
/// difference (`bcf_variant_t.n`, zero unless an indel/MNP).
pub fn variant_type(refa: &[u8], alt: &[u8]) -> (u32, i64) {
    let up = |c: u8| c.to_ascii_uppercase();
    if alt == b"*" {
        return (VCF_OVERLAP, 0);
    }
    if refa.len() == 1 && alt.len() == 1 {
        if alt[0] == b'.' || refa[0] == alt[0] || alt[0] == b'X' {
            return (VCF_REF, 0);
        }
        return (VCF_SNP, 0);
    }
    if alt.first() == Some(&b'<') {
        if alt == b"<X>" || alt == b"<*>" || alt == b"<NON_REF>" {
            return (VCF_REF, 0);
        }
        return (VCF_OTHER, 0);
    }
    if alt.first() == Some(&b']') || alt.first() == Some(&b'[') {
        return (VCF_BND, 0);
    }
    let (mut r, mut a) = (0usize, 0usize);
    while r < refa.len() && a < alt.len() && up(refa[r]) == up(alt[a]) {
        r += 1;
        a += 1;
    }
    if a < alt.len() && r == refa.len() {
        if alt[a] == b']' || alt[a] == b'[' {
            return (VCF_BND, 0);
        }
        return (VCF_INDEL, alt.len() as i64 - a as i64 - (refa.len() as i64 - r as i64));
    } else if r < refa.len() && a == alt.len() {
        return (VCF_INDEL, (alt.len() as i64 - a as i64) - (refa.len() as i64 - r as i64));
    } else if r == refa.len() && a == alt.len() {
        return (VCF_REF, 0);
    }
    let (mut re, mut ae) = (refa.len() - 1, alt.len() - 1);
    while re > r && ae > a && up(refa[re]) == up(alt[ae]) {
        re -= 1;
        ae -= 1;
    }
    if ae == a {
        if re == r {
            return (VCF_SNP, 0);
        }
        let n = -((re - r) as i64);
        return (if up(refa[re]) == up(alt[ae]) { VCF_INDEL } else { VCF_OTHER }, n);
    } else if re == r {
        let n = (ae - a) as i64;
        return (if up(refa[re]) == up(alt[ae]) { VCF_INDEL } else { VCF_OTHER }, n);
    }
    let ty = if re - r == ae - a { VCF_MNP } else { VCF_OTHER };
    let n = if re - r > ae - a { -((re - r + 1) as i64) } else { (ae - a + 1) as i64 };
    (ty, n)
}

fn iupac2bitmask(c: u8) -> i32 {
    match c.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        b'M' => 3,
        b'R' => 5,
        b'W' => 9,
        b'S' => 6,
        b'Y' => 10,
        b'K' => 12,
        b'V' => 7,
        b'H' => 11,
        b'D' => 13,
        b'B' => 14,
        b'N' => 15,
        _ => -1,
    }
}

fn bitmask2iupac(m: u8) -> u8 {
    const IUPAC: &[u8; 16] = b".ACMGRSVTWYHKDBN";
    if m == 0 || m > 15 { 0 } else { IUPAC[m as usize] }
}

/// Ungapped blocks between the reference and the consensus (`-c`).
struct Chain {
    block_lengths: Vec<i64>,
    ref_gaps: Vec<i64>,
    alt_gaps: Vec<i64>,
    ori_pos: i64,
    ref_last_block_ori: i64,
    alt_last_block_ori: i64,
}

impl Chain {
    fn new(ref_ori_pos: i64) -> Self {
        Self {
            block_lengths: Vec::new(),
            ref_gaps: Vec::new(),
            alt_gaps: Vec::new(),
            ori_pos: ref_ori_pos,
            ref_last_block_ori: ref_ori_pos,
            alt_last_block_ori: ref_ori_pos,
        }
    }

    fn push_gap(&mut self, ref_start: i64, ref_len: i64, alt_start: i64, alt_len: i64) {
        let num = self.block_lengths.len();
        if num > 0 && ref_start <= self.ref_last_block_ori {
            // Back-to-back with the previous variant.
            self.ref_last_block_ori = ref_start + ref_len;
            self.alt_last_block_ori = alt_start + alt_len;
            self.ref_gaps[num - 1] += ref_len;
            self.alt_gaps[num - 1] += alt_len;
        } else {
            self.block_lengths.push(ref_start - self.ref_last_block_ori);
            self.ref_gaps.push(ref_len);
            self.alt_gaps.push(alt_len);
            self.ref_last_block_ori = ref_start + ref_len;
            self.alt_last_block_ori = alt_start + alt_len;
        }
    }
}

/// Everything that is not the VCF and the FASTA.
pub struct Opts {
    pub output_iupac: bool,
    /// `-H N`; 0 when not given.
    pub haplotype: i64,
    /// `PICK_*` bits from `-H`.
    pub allele: u32,
    /// The single sample used with `-H`.
    pub isample: Option<usize>,
    /// `-s`/`-S` without `-H`: IUPAC codes from these samples' genotypes.
    pub iupac_gts: Option<Vec<usize>>,
    pub missing_allele: Option<u8>,
    pub absent_allele: Option<u8>,
    pub mark_del: Option<u8>,
    pub mark_ins: Mark,
    pub mark_snv: Mark,
    pub chr_prefix: Option<String>,
    pub masks: Vec<Mask>,
}

/// The FASTA writer state of consensus.c's `args_t`.
pub struct Consensus<'a> {
    opts: &'a Opts,
    vcf: &'a HashMap<String, Vec<Rec>>,
    /// Contigs the VCF header declares (bcf_hdr_name2id).
    contigs: &'a HashMap<String, ()>,
    out: Box<dyn Write>,
    chain_out: Option<Box<dyn Write>>,
    fa_buf: Vec<u8>,
    fa_ori_pos: i64,
    fa_frz_pos: i64,
    fa_mod_off: i64,
    fa_frz_mod: i64,
    fa_end_pos: i64,
    fa_length: i64,
    fa_src_pos: i64,
    prev_base: u8,
    prev_base_pos: i64,
    prev_is_insert: bool,
    /// Sequence name of the current region (FASTA header up to `:` or space).
    chr: String,
    /// The region's records and the read cursor; `None` when the VCF does
    /// not know the sequence.
    recs: Option<&'a [Rec]>,
    cursor: usize,
    chain: Option<Chain>,
    chain_id: u32,
    iupac_als: Vec<bool>,
    iupac_bitmask: Vec<u8>,
    napplied: u64,
    warned_haplotype: bool,
    multibase_ref_del_warned: bool,
}

impl<'a> Consensus<'a> {
    pub fn new(
        opts: &'a Opts,
        vcf: &'a HashMap<String, Vec<Rec>>,
        contigs: &'a HashMap<String, ()>,
        out: Box<dyn Write>,
        chain_out: Option<Box<dyn Write>>,
    ) -> Self {
        if opts.isample.is_none() && opts.iupac_gts.is_none() {
            eprintln!("Note: the --samples option not given, applying all records regardless of the genotype");
        }
        Self {
            opts,
            vcf,
            contigs,
            out,
            chain_out,
            fa_buf: Vec::new(),
            fa_ori_pos: 0,
            fa_frz_pos: -1,
            fa_mod_off: 0,
            fa_frz_mod: -1,
            fa_end_pos: 0,
            fa_length: 0,
            fa_src_pos: 0,
            prev_base: 0,
            prev_base_pos: -1,
            prev_is_insert: false,
            chr: String::new(),
            recs: None,
            cursor: 0,
            chain: None,
            chain_id: 0,
            iupac_als: Vec::new(),
            iupac_bitmask: Vec::new(),
            napplied: 0,
            warned_haplotype: false,
            multibase_ref_del_warned: false,
        }
    }

    /// Number of variants applied.
    pub fn applied(&self) -> u64 {
        self.napplied
    }

    /// Run over a FASTA (`consensus()` in C).
    pub fn run<R: BufRead>(&mut self, fasta: R) -> Result<()> {
        for line in fasta.lines() {
            let line = line.context("read FASTA")?;
            let line = line.trim_end_matches('\r');
            if line.is_empty() {
                continue;
            }
            if let Some(name) = line.strip_prefix('>') {
                self.print_chain()?;
                // Apply the cached variants and those a short FASTA may have missed.
                while self.recs.is_some() {
                    let Some(i) = self.next_vcf_line() else { break };
                    let rec = &self.recs.unwrap()[i];
                    if self.fa_end_pos != 0 && rec.pos > self.fa_end_pos {
                        break;
                    }
                    self.apply_variant(rec)?;
                }
                if self.opts.absent_allele.is_some() {
                    let pos = self.peek_pos().unwrap_or(0);
                    self.apply_absent(pos);
                }
                self.flush_fa_buffer(0)?;
                self.init_region(name)?;
                continue;
            }
            let mut seq = line.as_bytes().to_vec();
            self.fa_length += seq.len() as i64;
            self.fa_src_pos += seq.len() as i64;
            if !self.opts.masks.is_empty() {
                self.mask_region(&mut seq);
            }
            self.fa_buf.extend_from_slice(&seq);
            let mut pending = false;
            while self.recs.is_some() {
                let Some(i) = self.next_vcf_line() else { break };
                let rec = &self.recs.unwrap()[i];
                // Still the same region? If not, keep the record and flush.
                if self.fa_end_pos != 0 && rec.pos > self.fa_end_pos {
                    self.unread_vcf_line();
                    break;
                }
                // Record well beyond the cached FASTA: flush.
                if self.fa_ori_pos + self.fa_buf.len() as i64 - self.fa_mod_off <= rec.pos {
                    self.unread_vcf_line();
                    break;
                }
                // Cached FASTA too short for the record: read more, no flushing.
                if self.fa_ori_pos + self.fa_buf.len() as i64 - self.fa_mod_off < rec.pos + rec.rlen {
                    self.unread_vcf_line();
                    pending = true;
                    break;
                }
                self.apply_variant(rec)?;
            }
            if !pending {
                if self.opts.absent_allele.is_some() {
                    let pos = self.fa_ori_pos - self.fa_mod_off + self.fa_buf.len() as i64;
                    self.apply_absent(pos);
                }
                self.flush_fa_buffer(60)?;
            }
        }
        while self.recs.is_some() {
            let Some(i) = self.next_vcf_line() else { break };
            let rec = &self.recs.unwrap()[i];
            if self.fa_end_pos != 0 && rec.pos > self.fa_end_pos {
                break;
            }
            if self.fa_ori_pos + self.fa_buf.len() as i64 - self.fa_mod_off <= rec.pos {
                break;
            }
            self.apply_variant(rec)?;
        }
        self.print_chain()?;
        if self.opts.absent_allele.is_some() {
            self.apply_absent(i64::MAX);
        }
        self.flush_fa_buffer(0)?;
        self.out.flush()?;
        if let Some(c) = self.chain_out.as_mut() {
            c.flush()?;
        }
        Ok(())
    }

    fn next_vcf_line(&mut self) -> Option<usize> {
        let recs = self.recs?;
        if self.cursor < recs.len() {
            self.cursor += 1;
            Some(self.cursor - 1)
        } else {
            None
        }
    }

    fn unread_vcf_line(&mut self) {
        self.cursor -= 1;
    }

    /// Position of the buffered (unread) record, as `args->vcf_buf[f]->pos`.
    fn peek_pos(&self) -> Option<i64> {
        self.recs.and_then(|r| r.get(self.cursor)).map(|r| r.pos)
    }

    fn init_region(&mut self, header: &str) -> Result<()> {
        let name_end = header.find(|c: char| c.is_whitespace() || c == ':').unwrap_or(header.len());
        let name = &header[..name_end];
        let (mut from, mut to) = (0i64, 0i64);
        if name_end < header.len() && header.as_bytes()[name_end] == b':' {
            let rest = &header[name_end + 1..];
            let digits = |s: &str| s.bytes().take_while(u8::is_ascii_digit).count();
            let nf = digits(rest);
            if nf == 0 || rest[nf..].as_bytes().first() != Some(&b'-') {
                from = 0;
            } else {
                from = rest[..nf].parse::<i64>().unwrap_or(0) - 1;
                let tail = &rest[nf + 1..];
                let nt = digits(tail);
                if nt == 0 || tail[nt..].chars().next().is_some_and(|c| !c.is_whitespace()) {
                    from = 0;
                    to = 0;
                } else {
                    to = tail[..nt].parse::<i64>().unwrap_or(0) - 1;
                }
            }
        }
        self.chr = name.to_string();
        let known = self.vcf.contains_key(name) || self.contigs.contains_key(name);
        if !known {
            eprintln!("Warning: Sequence \"{name}\" not in the VCF");
        }
        self.recs = if known { Some(self.vcf.get(name).map(Vec::as_slice).unwrap_or(&[])) } else { None };
        self.cursor = 0;
        self.prev_base_pos = -1;
        self.fa_buf.clear();
        self.fa_length = 0;
        self.fa_end_pos = to;
        self.fa_ori_pos = from;
        self.fa_src_pos = from;
        self.fa_mod_off = 0;
        self.fa_frz_pos = -1;
        self.fa_frz_mod = -1;
        writeln!(self.out, ">{}{}", self.opts.chr_prefix.as_deref().unwrap_or(""), header)?;
        if self.chain_out.is_some() {
            self.chain = Some(Chain::new(self.fa_ori_pos));
        }
        Ok(())
    }

    fn print_chain(&mut self) -> Result<()> {
        let (Some(chain), Some(out)) = (self.chain.as_ref(), self.chain_out.as_mut()) else {
            return Ok(());
        };
        let ref_end_pos = self.fa_length + chain.ori_pos;
        let last_block_size = ref_end_pos - chain.ref_last_block_ori;
        let alt_end_pos = chain.alt_last_block_ori + last_block_size;
        let score: i64 = chain.block_lengths.iter().sum::<i64>() + last_block_size;
        self.chain_id += 1;
        writeln!(
            out,
            "chain {} {} {} + {} {} {} {} + {} {} {}",
            score, self.chr, ref_end_pos, chain.ori_pos, ref_end_pos, self.chr, alt_end_pos, chain.ori_pos, alt_end_pos, self.chain_id
        )?;
        for n in 0..chain.block_lengths.len() {
            writeln!(out, "{} {} {}", chain.block_lengths[n], chain.ref_gaps[n], chain.alt_gaps[n])?;
        }
        writeln!(out, "{}\n", last_block_size)?;
        Ok(())
    }

    fn flush_fa_buffer(&mut self, len: usize) -> Result<()> {
        if self.fa_buf.is_empty() {
            return Ok(());
        }
        let mut nwr = 0usize;
        while nwr + 60 <= self.fa_buf.len() {
            self.out.write_all(&self.fa_buf[nwr..nwr + 60])?;
            self.out.write_all(b"\n")?;
            nwr += 60;
        }
        if nwr > 0 {
            self.fa_ori_pos += nwr as i64;
        }
        self.fa_frz_mod -= nwr as i64;
        if len > 0 {
            // Not finished on this sequence yet: keep the tail.
            self.fa_buf.drain(..nwr);
            return Ok(());
        }
        if nwr == self.fa_buf.len() {
            self.fa_buf.clear();
            return Ok(());
        }
        self.out.write_all(&self.fa_buf[nwr..])?;
        self.out.write_all(b"\n")?;
        self.fa_ori_pos += (self.fa_buf.len() - nwr) as i64 - self.fa_mod_off;
        self.fa_mod_off = 0;
        self.fa_buf.clear();
        Ok(())
    }

    fn apply_absent(&mut self, pos: i64) {
        let Some(absent) = self.opts.absent_allele else { return };
        if self.fa_buf.is_empty() || pos <= self.fa_frz_pos + 1 || pos <= self.fa_ori_pos {
            return;
        }
        let l = self.fa_buf.len() as i64;
        let ie = if pos != 0 && pos - self.fa_ori_pos + self.fa_mod_off < l { pos - self.fa_ori_pos + self.fa_mod_off } else { l };
        let ib = if self.fa_frz_mod < 0 { 0 } else { self.fa_frz_mod };
        let mut i = ib;
        while i < ie {
            self.fa_buf[i as usize] = absent;
            i += 1;
        }
    }

    fn freeze_ref(&mut self, pos: i64, rlen: i64) {
        if self.fa_frz_pos >= pos + rlen - 1 {
            return;
        }
        self.fa_frz_pos = pos + rlen - 1;
        self.fa_frz_mod = pos - self.fa_ori_pos + self.fa_mod_off + rlen;
    }

    fn mask_region(&mut self, seq: &mut [u8]) {
        let len = seq.len() as i64;
        let start = self.fa_src_pos - len;
        let end = self.fa_src_pos;
        for mask in &self.opts.masks {
            for (rb, re) in mask.overlapping(&self.chr, start, end) {
                let idx_start = (rb - start).max(0);
                let idx_end = (re - start).min(len - 1);
                let mut j = idx_start;
                while j <= idx_end {
                    let c = &mut seq[j as usize];
                    match mask.with {
                        MaskWith::Upper => *c = c.to_ascii_uppercase(),
                        MaskWith::Lower => *c = c.to_ascii_lowercase(),
                        MaskWith::Char(m) => *c = m,
                    }
                    j += 1;
                }
            }
        }
    }

    fn iupac_init(&mut self, n_allele: usize) {
        self.iupac_als.clear();
        self.iupac_als.resize(n_allele, false);
    }

    fn iupac_add_gt(&mut self, gt: &[GtVal], n_allele: usize) -> Result<bool> {
        let mut is_set = false;
        for &g in gt {
            match g {
                GtVal::Missing { .. } => continue,
                GtVal::End => break,
                GtVal::Allele { idx, .. } => {
                    if idx < 0 || idx as usize >= n_allele {
                        bail!("Invalid VCF, too few ALT alleles at {}:{}", self.chr, self.cur_pos1());
                    }
                    self.iupac_als[idx as usize] = true;
                    is_set = true;
                }
            }
        }
        Ok(is_set)
    }

    fn iupac_set_allele(&mut self, alleles: &mut [Vec<u8>]) -> i32 {
        let (mut max_len, mut alt_len, mut ialt, mut fallback_alt) = (0usize, 0usize, -1i32, -1i32);
        for (i, a) in alleles.iter().enumerate() {
            if !self.iupac_als[i] {
                continue;
            }
            if fallback_alt <= 0 {
                fallback_alt = i as i32;
            }
            if a.iter().any(|&c| iupac2bitmask(c) < 0) {
                continue;
            }
            if a.len() > max_len {
                self.iupac_bitmask.resize(a.len(), 0);
                for j in max_len..a.len() {
                    self.iupac_bitmask[j] = 0;
                }
                max_len = a.len();
            }
            if i > 0 && a.len() > alt_len {
                alt_len = a.len();
                ialt = i as i32;
            }
            for (j, &c) in a.iter().enumerate() {
                self.iupac_bitmask[j] |= iupac2bitmask(c) as u8;
            }
        }
        if alt_len > 0 {
            for j in 0..alt_len {
                alleles[ialt as usize][j] = bitmask2iupac(self.iupac_bitmask[j]);
            }
            return ialt;
        }
        if fallback_alt >= 0 {
            return fallback_alt;
        }
        ialt
    }

    /// 1-based position of the record being applied, for messages.
    fn cur_pos1(&self) -> i64 {
        self.recs.and_then(|r| r.get(self.cursor.wrapping_sub(1))).map_or(0, |r| r.pos + 1)
    }

    fn apply_variant(&mut self, rec: &Rec) -> Result<()> {
        let pos1 = rec.pos + 1;
        if self.opts.absent_allele.is_some() {
            self.apply_absent(rec.pos);
        }
        if rec.alleles.len() == 1 && self.opts.missing_allele.is_none() && self.opts.absent_allele.is_none() {
            return Ok(());
        }
        if !self.opts.masks.is_empty() {
            let (start, end) = (rec.pos, rec.pos + rec.rlen - 1);
            if self.opts.masks.iter().any(|m| m.skips() && m.overlaps(&self.chr, start, end)) {
                return Ok(());
            }
        }
        let mut alleles = rec.alleles.clone();
        let mut rlen = rec.rlen;
        let n_allele = alleles.len();
        let mut ialt: i32 = 1;
        if let Some(smpl) = &self.opts.iupac_gts {
            let Some(gts) = &rec.gts else { return Ok(()) };
            ialt = -1;
            self.iupac_init(n_allele);
            let mut is_set = false;
            for &i in smpl {
                is_set |= self.iupac_add_gt(&gts[i], n_allele)?;
            }
            if !is_set && self.opts.missing_allele.is_none() {
                return Ok(());
            }
            if is_set {
                ialt = self.iupac_set_allele(&mut alleles);
            }
        } else if let Some(isample) = self.opts.isample {
            let Some(gts) = &rec.gts else { return Ok(()) };
            let ptr = &gts[isample];
            let n = ptr.len();
            #[derive(PartialEq)]
            enum Action {
                UseHap,
                UseIupac,
                PickOne,
            }
            let mut action = Action::UseHap;
            if self.opts.allele == PICK_IUPAC {
                if self.opts.haplotype == 0 {
                    action = Action::UseIupac;
                }
                if !ptr[0].is_phased() && !ptr[n - 1].is_phased() {
                    action = Action::UseIupac;
                }
            } else if self.opts.output_iupac {
                action = Action::UseIupac;
            } else if self.opts.haplotype == 0 {
                action = Action::PickOne;
            }
            match action {
                Action::UseHap => {
                    if self.opts.haplotype > n as i64 {
                        if ptr[n - 1].is_missing() || ptr[0].is_missing() {
                            if self.opts.missing_allele.is_none() {
                                return Ok(());
                            }
                            ialt = -1;
                        } else {
                            if !self.warned_haplotype {
                                eprintln!(
                                    "Can't apply {}-th haplotype at {}:{}. (This warning is printed only once.)",
                                    self.opts.haplotype, self.chr, pos1
                                );
                                self.warned_haplotype = true;
                            }
                            return Ok(());
                        }
                    } else {
                        match ptr[self.opts.haplotype as usize - 1] {
                            GtVal::Allele { idx, .. } => ialt = idx,
                            _ => {
                                if self.opts.missing_allele.is_none() {
                                    return Ok(());
                                }
                                ialt = -1;
                            }
                        }
                    }
                }
                Action::UseIupac => {
                    ialt = -1;
                    self.iupac_init(n_allele);
                    let is_set = self.iupac_add_gt(ptr, n_allele)?;
                    if !is_set && self.opts.missing_allele.is_none() {
                        return Ok(());
                    }
                    if is_set {
                        ialt = self.iupac_set_allele(&mut alleles);
                    }
                }
                Action::PickOne => {
                    let mut is_hom = true;
                    for i in 0..n {
                        match ptr[i] {
                            GtVal::Missing { .. } => {
                                if self.opts.missing_allele.is_none() {
                                    return Ok(());
                                }
                                ialt = -1;
                                break;
                            }
                            GtVal::End => break,
                            GtVal::Allele { idx, .. } => {
                                ialt = idx;
                                if i > 0 && idx != ptr[i - 1].allele() {
                                    is_hom = false;
                                    break;
                                }
                            }
                        }
                    }
                    if !is_hom {
                        let mut prev_len = 0i64;
                        for i in 0..n {
                            if ptr[i] == GtVal::End {
                                break;
                            }
                            let jalt = ptr[i].allele();
                            if jalt < 0 || n_allele as i32 <= jalt {
                                bail!("Broken VCF, too few alts at {}:{}", self.chr, pos1);
                            }
                            let a = self.opts.allele;
                            if a & (PICK_LONG | PICK_SHORT) != 0 {
                                let len = if jalt == 0 { rlen } else { alleles[jalt as usize].len() as i64 };
                                if i == 0 {
                                    ialt = jalt;
                                    prev_len = len;
                                } else if len == prev_len {
                                    if a & PICK_REF != 0 && jalt == 0 {
                                        ialt = jalt;
                                        prev_len = len;
                                    } else if a & PICK_ALT != 0 && ialt == 0 {
                                        ialt = jalt;
                                        prev_len = len;
                                    }
                                } else if a & PICK_LONG != 0 && len > prev_len {
                                    ialt = jalt;
                                    prev_len = len;
                                } else if a & PICK_SHORT != 0 && len < prev_len {
                                    ialt = jalt;
                                    prev_len = len;
                                }
                            } else if a & PICK_REF != 0 && jalt == 0 {
                                ialt = jalt;
                            } else if a & PICK_ALT != 0 && ialt == 0 {
                                ialt = jalt;
                            }
                        }
                    }
                }
            }
            if ialt == 0 {
                if self.opts.absent_allele.is_some() {
                    self.freeze_ref(rec.pos, rlen);
                }
                return Ok(());
            }
            if n_allele as i32 <= ialt {
                bail!("Broken VCF, too few alts at {}:{}", self.chr, pos1);
            }
        } else if self.opts.output_iupac && n_allele > 1 {
            let (mut alen, mut mlen, mut ialt_l) = (0usize, 0usize, 0usize);
            for i in 0..n_allele {
                let a = &alleles[i];
                if a.iter().any(|&c| iupac2bitmask(c) < 0) {
                    continue;
                }
                if a.len() > mlen {
                    self.iupac_bitmask.resize(a.len(), 0);
                    for j in mlen..a.len() {
                        self.iupac_bitmask[j] = 0;
                    }
                    mlen = a.len();
                }
                if i > 0 && a.len() > alen {
                    alen = a.len();
                    ialt_l = i;
                }
                for (j, &c) in a.iter().enumerate() {
                    self.iupac_bitmask[j] |= iupac2bitmask(c) as u8;
                }
            }
            if alen > 0 {
                for j in 0..alen {
                    alleles[ialt_l][j] = bitmask2iupac(self.iupac_bitmask[j]);
                }
                ialt = ialt_l as i32;
            } else {
                ialt = 1;
            }
        }
        if n_allele == 1 && ialt != -1 {
            if self.opts.absent_allele.is_some() {
                self.freeze_ref(rec.pos, rlen);
            }
            return Ok(());
        }
        if ialt == -1 {
            let missing = self.opts.missing_allele.unwrap_or(b'.');
            alleles = vec![vec![alleles[0][0]], vec![missing]];
            if !rec.has_end {
                rlen = 1;
            }
            ialt = 1;
        }
        let ialt = ialt as usize;
        // For some types POS+REF refer to the base before the event.
        let mut trim_beg = 0i64;
        let (var_type, mut var_len) = variant_type(&alleles[0], &alleles[ialt]);
        if var_type & VCF_INDEL != 0 {
            trim_beg = if alleles[0][0] == alleles[ialt][0] { 1 } else { 0 };
        } else if var_type & VCF_OTHER != 0 && alleles[ialt].eq_ignore_ascii_case(b"<DEL>") {
            trim_beg = 1;
            var_len = 1 - rlen;
        } else if var_type & VCF_OTHER != 0 && alleles[ialt].len() >= 4 && alleles[ialt][..4].eq_ignore_ascii_case(b"<INS") {
            trim_beg = 1;
        }
        if rec.pos <= self.fa_frz_pos {
            let overlap = rec.pos < self.fa_frz_pos || trim_beg == 0 || var_len == 0 || self.prev_is_insert;
            if overlap {
                eprintln!("The site {}:{} overlaps with another variant, skipping...", self.chr, pos1);
                return Ok(());
            }
        }
        let mut alt_allele = alleles[ialt].clone();
        let mut len_diff: i64;
        let mut alen: i64;
        let idx = rec.pos - self.fa_ori_pos + self.fa_mod_off;
        if idx < 0 {
            eprintln!("Warning: ignoring overlapping variant starting at {}:{}", self.chr, pos1);
            return Ok(());
        }
        let buf_len = self.fa_buf.len() as i64;
        if rlen > buf_len - idx {
            rlen = buf_len - idx;
            if alt_allele[0] != b'<' {
                alen = alt_allele.len() as i64;
                if alen > rlen {
                    eprintln!(
                        "Warning: trimming variant \"{}\" starting at {}:{}",
                        String::from_utf8_lossy(&alt_allele),
                        self.chr,
                        pos1
                    );
                    alt_allele.truncate(rlen as usize);
                }
            }
        }
        if idx >= buf_len {
            bail!(
                "FIXME: {}:{} .. idx={}, ori_pos={}, len={}, off={}",
                self.chr,
                pos1,
                idx,
                self.fa_ori_pos,
                buf_len,
                self.fa_mod_off
            );
        }
        let idx = idx as usize;
        let refa = alleles[0].clone();
        if alt_allele[0] == b'<' {
            if !alt_allele.eq_ignore_ascii_case(b"<DEL>")
                && !alt_allele.eq_ignore_ascii_case(b"<*>")
                && !alt_allele.eq_ignore_ascii_case(b"<NON_REF>")
            {
                bail!(
                    "Symbolic alleles other than <DEL>, <*> or <NON_REF> are currently not supported, e.g. \"{}\" at {}:{}.\nPlease use filtering expressions to exclude such sites, for example by running with: -e 'ALT~\"<.*>\"'",
                    String::from_utf8_lossy(&alt_allele),
                    self.chr,
                    pos1
                );
            }
            if alt_allele.eq_ignore_ascii_case(b"<DEL>") {
                if refa.len() > 1 && !self.multibase_ref_del_warned {
                    eprintln!(
                        "Warning: one REF base is expected with <DEL>, assuming the actual deletion starts at POS+1 at {}:{}.\n         (This warning is printed only once.)",
                        self.chr, pos1
                    );
                    self.multibase_ref_del_warned = true;
                }
                if let Some(m) = self.opts.mark_del {
                    alt_allele = mark_del(&refa, rlen as usize, None, m);
                    alen = rlen;
                    len_diff = 0;
                } else {
                    len_diff = 1 - rlen;
                    alt_allele = refa.clone();
                    alen = 1;
                }
            } else {
                // <*> or <NON_REF>: evidence for the reference over the block.
                self.freeze_ref(rec.pos, rlen);
                return Ok(());
            }
        } else if !strncasecmp_eq(&refa, &self.fa_buf[idx..], rlen as usize) {
            // A SNP or insert right after a deletion loses its reference base
            // in the buffer; allow that one-base overlap.
            let mut fail = true;
            if self.prev_base_pos == rec.pos && refa[0].to_ascii_uppercase() == self.prev_base.to_ascii_uppercase() {
                if rlen == 1 {
                    fail = false;
                } else if strncasecmp_eq(&refa[1..], &self.fa_buf[idx + 1..], rlen as usize - 1) {
                    fail = false;
                }
            }
            if fail {
                let shown = &self.fa_buf[idx..(idx + rlen as usize).min(self.fa_buf.len())];
                bail!(
                    "The fasta sequence does not match the REF allele at {}:{}:\n   REF .vcf: [{}]\n   ALT .vcf: [{}]\n   REF .fa : [{}]",
                    self.chr,
                    pos1,
                    String::from_utf8_lossy(&refa),
                    String::from_utf8_lossy(&alt_allele),
                    String::from_utf8_lossy(shown)
                );
            }
            alen = alt_allele.len() as i64;
            len_diff = alen - rlen;
            if let (Some(m), true) = (self.opts.mark_del, len_diff < 0) {
                alt_allele = mark_del(&refa, rlen as usize, Some(&alt_allele), m);
                alen = rlen;
                len_diff = 0;
            }
        } else {
            alen = alt_allele.len() as i64;
            len_diff = alen - rlen;
            if let (Some(m), true) = (self.opts.mark_del, len_diff < 0) {
                alt_allele = mark_del(&refa, rlen as usize, Some(&alt_allele), m);
                alen = rlen;
                len_diff = 0;
            }
        }
        let alen_u = alen as usize;
        let rlen_u = rlen as usize;
        let upper = self.fa_buf[idx].to_ascii_uppercase() == self.fa_buf[idx];
        for c in alt_allele.iter_mut().take(alen_u) {
            *c = if upper { c.to_ascii_uppercase() } else { c.to_ascii_lowercase() };
        }
        if self.opts.mark_ins != Mark::None && len_diff > 0 {
            mark_ins(&refa, &mut alt_allele, self.opts.mark_ins);
        }
        if self.opts.mark_snv != Mark::None {
            mark_snv(&refa, &mut alt_allele, self.opts.mark_snv);
        }
        if len_diff <= 0 {
            // Deletion or same-size event.
            assert!(self.fa_buf.len() >= idx + rlen_u);
            self.prev_base = self.fa_buf[idx + rlen_u - 1];
            self.prev_base_pos = rec.pos + rlen - 1;
            self.prev_is_insert = false;
            self.fa_frz_mod = idx as i64 + alen;
            for i in trim_beg as usize..alen_u {
                self.fa_buf[idx + i] = alt_allele[i];
            }
            if len_diff != 0 {
                self.fa_buf.drain(idx + alen_u..idx + rlen_u);
            }
        } else {
            self.prev_is_insert = true;
            self.prev_base_pos = rec.pos;
            let ins = std::iter::repeat_n(0u8, len_diff as usize);
            self.fa_buf.splice(idx + rlen_u..idx + rlen_u, ins);
            // Bases unchanged by the insertion must not overwrite preceding
            // variants (C>T then C>CAA at one position gives TAA).
            let mut ibeg = 0usize;
            while ibeg < alen_u && refa.get(ibeg) == alt_allele.get(ibeg) && rec.pos + ibeg as i64 <= self.prev_base_pos {
                ibeg += 1;
            }
            for i in ibeg..alen_u {
                self.fa_buf[idx + i] = alt_allele[i];
            }
            self.fa_frz_mod = idx as i64 + alen - ibeg as i64 + 1;
        }
        if let (Some(chain), true) = (self.chain.as_mut(), len_diff != 0) {
            if refa[0].eq_ignore_ascii_case(&alt_allele[0]) {
                // Indels carry the base before the event: block is 1 bp longer, alleles 1 bp shorter.
                chain.push_gap(rec.pos + 1, rlen - 1, rec.pos + 1 + self.fa_mod_off, alen - 1);
            } else {
                chain.push_gap(rec.pos, rlen, rec.pos + self.fa_mod_off, alen);
            }
        }
        self.fa_mod_off += len_diff;
        self.fa_frz_pos = rec.pos + rlen - 1;
        self.napplied += 1;
        Ok(())
    }
}

/// C `strncasecmp(a, b, n) == 0` with `a` NUL-terminated at its end.
fn strncasecmp_eq(a: &[u8], b: &[u8], n: usize) -> bool {
    for i in 0..n {
        match (a.get(i), b.get(i)) {
            (None, None) => return true,
            (Some(x), Some(y)) if x.eq_ignore_ascii_case(y) => continue,
            _ => return false,
        }
    }
    true
}

fn mark_del(refa: &[u8], rlen: usize, alt: Option<&[u8]>, mark: u8) -> Vec<u8> {
    let mut out = Vec::with_capacity(rlen);
    out.extend_from_slice(alt.unwrap_or(refa));
    out.truncate(rlen);
    out.resize(rlen, mark);
    out
}

fn mark_ins(refa: &[u8], alt: &mut [u8], mark: Mark) {
    for c in alt.iter_mut().skip(refa.len()) {
        match mark {
            Mark::Lower => *c = c.to_ascii_lowercase(),
            Mark::Upper => *c = c.to_ascii_uppercase(),
            Mark::Char(m) => *c = m,
            Mark::None => {}
        }
    }
}

fn mark_snv(refa: &[u8], alt: &mut [u8], mark: Mark) {
    let n = refa.len().min(alt.len());
    for i in 0..n {
        if refa[i].to_ascii_lowercase() != alt[i].to_ascii_lowercase() {
            match mark {
                Mark::Lower => alt[i] = alt[i].to_ascii_lowercase(),
                Mark::Upper => alt[i] = alt[i].to_ascii_uppercase(),
                Mark::Char(m) => alt[i] = m,
                Mark::None => {}
            }
        }
    }
}

#[cfg(test)]
#[path = "../tests/unit/consensus.rs"]
mod tests;
