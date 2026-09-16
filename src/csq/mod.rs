//! Haplotype-aware consequence calling (`csq`): a port of bcftools `csq.c`.
//!
//! Records are buffered per position while transcripts overlapping them are
//! still open; every sample haplotype grows a path through a per-transcript
//! tree of variant nodes, and when the transcript closes the paths are
//! spliced, translated and broken into independent consequences. Coding
//! consequences are printed once (`@pos` references the record carrying the
//! text) and FORMAT/BCSQ carries a per-haplotype bitmask into INFO/BCSQ.

mod caller;
pub mod gff;
mod hap;
mod splice;

use std::collections::VecDeque;
use std::rc::Rc;

use anyhow::{Result, bail};
use crate::fx::FxHashMap;

pub use caller::{Caller, CsqOpts, CsqOutput, Phase};
use gff::{Gff, STRAND_FWD, gf_type2gff_string};

pub const N_REF_PAD: usize = 10;

pub const CSQ_PRINTED_UPSTREAM: u32 = 1 << 0;
pub const CSQ_SYNONYMOUS_VARIANT: u32 = 1 << 1;
pub const CSQ_MISSENSE_VARIANT: u32 = 1 << 2;
pub const CSQ_STOP_LOST: u32 = 1 << 3;
pub const CSQ_STOP_GAINED: u32 = 1 << 4;
pub const CSQ_INFRAME_DELETION: u32 = 1 << 5;
pub const CSQ_INFRAME_INSERTION: u32 = 1 << 6;
pub const CSQ_FRAMESHIFT_VARIANT: u32 = 1 << 7;
pub const CSQ_SPLICE_ACCEPTOR: u32 = 1 << 8;
pub const CSQ_SPLICE_DONOR: u32 = 1 << 9;
pub const CSQ_START_LOST: u32 = 1 << 10;
pub const CSQ_SPLICE_REGION: u32 = 1 << 11;
pub const CSQ_STOP_RETAINED: u32 = 1 << 12;
pub const CSQ_UTR5: u32 = 1 << 13;
pub const CSQ_UTR3: u32 = 1 << 14;
pub const CSQ_NON_CODING: u32 = 1 << 15;
pub const CSQ_INTRON: u32 = 1 << 16;
pub const CSQ_INFRAME_ALTERING: u32 = 1 << 18;
pub const CSQ_UPSTREAM_STOP: u32 = 1 << 19;
pub const CSQ_INCOMPLETE_CDS: u32 = 1 << 20;
pub const CSQ_CODING_SEQUENCE: u32 = 1 << 21;
pub const CSQ_ELONGATION: u32 = 1 << 22;
pub const CSQ_START_RETAINED: u32 = 1 << 23;

/// Consequences printed on one record only; the others reference it.
pub const CSQ_COMPOUND: u32 = CSQ_SYNONYMOUS_VARIANT
    | CSQ_MISSENSE_VARIANT
    | CSQ_STOP_LOST
    | CSQ_STOP_GAINED
    | CSQ_INFRAME_DELETION
    | CSQ_INFRAME_INSERTION
    | CSQ_FRAMESHIFT_VARIANT
    | CSQ_START_LOST
    | CSQ_STOP_RETAINED
    | CSQ_INFRAME_ALTERING
    | CSQ_INCOMPLETE_CDS
    | CSQ_UPSTREAM_STOP
    | CSQ_START_RETAINED;
pub const CSQ_START_STOP: u32 =
    CSQ_STOP_LOST | CSQ_STOP_GAINED | CSQ_STOP_RETAINED | CSQ_START_LOST | CSQ_START_RETAINED;
pub(crate) const CSQ_PRN_TSCRIPT: u32 = !(CSQ_INTRON | CSQ_NON_CODING);
pub(crate) const CSQ_PRN_NMD: u32 = !(CSQ_INTRON | CSQ_NON_CODING);

#[inline]
pub(crate) fn csq_prn_strand(c: u32) -> bool {
    c & CSQ_COMPOUND != 0 && c & (CSQ_SPLICE_ACCEPTOR | CSQ_SPLICE_DONOR | CSQ_SPLICE_REGION) == 0
}

const CSQ_STRINGS: [Option<&str>; 24] = [
    None,
    Some("synonymous"),
    Some("missense"),
    Some("stop_lost"),
    Some("stop_gained"),
    Some("inframe_deletion"),
    Some("inframe_insertion"),
    Some("frameshift"),
    Some("splice_acceptor"),
    Some("splice_donor"),
    Some("start_lost"),
    Some("splice_region"),
    Some("stop_retained"),
    Some("5_prime_utr"),
    Some("3_prime_utr"),
    Some("non_coding"),
    Some("intron"),
    Some("intergenic"),
    Some("inframe_altering"),
    None,
    None,
    Some("coding_sequence"),
    Some("feature_elongation"),
    Some("start_retained"),
];

pub(crate) const HAP_CDS: u8 = 0;
pub(crate) const HAP_ROOT: u8 = 1;
pub(crate) const HAP_SSS: u8 = 2;

/// Standard genetic code indexed by `nt4(a)<<4 | nt4(b)<<2 | nt4(c)`.
const GENCODE: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

#[inline]
fn nt4(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}
#[inline]
fn cnt4(b: u8) -> usize {
    match b {
        b'A' | b'a' => 3,
        b'C' | b'c' => 2,
        b'G' | b'g' => 1,
        b'T' | b't' => 0,
        _ => 4,
    }
}

/// Amino acid of a forward-strand codon; `X` when a base is not ACGT and the
/// lookup falls outside the table.
#[inline]
pub fn dna2aa(x: &[u8]) -> u8 {
    let i = nt4(x[0]) << 4 | nt4(x[1]) << 2 | nt4(x[2]);
    GENCODE.get(i).copied().unwrap_or(b'X')
}

/// Amino acid of a reverse-strand codon given on the reference strand.
#[inline]
pub fn cdna2aa(x: &[u8]) -> u8 {
    let i = cnt4(x[2]) << 4 | cnt4(x[1]) << 2 | cnt4(x[0]);
    GENCODE.get(i).copied().unwrap_or(b'X')
}

/// `ncsq2_to_nfmt`: integers per sample needed for `ncsq2` bits (one bit of
/// each `int32` is left unused so the value never equals the BCF missing
/// marker).
pub fn ncsq2_to_nfmt(ncsq2: usize) -> usize {
    1 + (ncsq2 - 1) / 30
}

#[inline]
pub fn icsq2_to_bit(icsq2: usize) -> (usize, usize) {
    (icsq2 / 30, icsq2 % 30)
}

pub type RecId = u64;

/// One genotype allele; `allele < 0` is missing, `phased` mirrors the `|`
/// separator in front of it.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Gt {
    pub allele: i32,
    pub phased: bool,
}

/// The parts of a VCF record the caller works with. Shared through `Rc` so a
/// record can sit in the output buffer while haplotype nodes refer to it.
#[derive(Clone, Debug)]
pub struct RecInfo {
    pub id: RecId,
    pub chrom: String,
    /// 0-based.
    pub pos: u32,
    pub rlen: u32,
    /// REF first.
    pub alleles: Rc<Vec<String>>,
    /// Per header sample; shorter vectors are padded with vector-end.
    pub gts: Rc<Vec<Vec<Gt>>>,
    /// Maximum ploidy (`bcf_get_genotypes` / nsamples); 0 without GT.
    pub ngts: usize,
}

impl RecInfo {
    pub fn n_allele(&self) -> usize {
        self.alleles.len()
    }
    /// `None` is the vector-end padding of a shorter genotype.
    #[inline]
    pub fn gt(&self, ismpl: usize, ihap: usize) -> Option<Gt> {
        self.gts.get(ismpl).and_then(|g| g.get(ihap)).copied()
    }
    pub fn allele_is_symbolic(&self, ial: usize) -> bool {
        let a = self.alleles[ial].as_bytes();
        matches!(a.first(), Some(b'<') | Some(b'*'))
    }
}

/// Parse the GT column of every sample: `(gts, max ploidy)`.
pub fn parse_genotypes(format: Option<&str>, samples: &[String]) -> (Vec<Vec<Gt>>, usize) {
    let Some(fmt) = format else {
        return (Vec::new(), 0);
    };
    let Some(igt) = fmt.split(':').position(|k| k == "GT") else {
        return (Vec::new(), 0);
    };
    let mut out = Vec::with_capacity(samples.len());
    let mut ngts = 0usize;
    for s in samples {
        let field = s.split(':').nth(igt).unwrap_or(".");
        let mut v = Vec::with_capacity(2);
        let mut phased = false;
        let mut tok = String::new();
        let push = |tok: &str, phased: bool, v: &mut Vec<Gt>| {
            let allele = if tok == "." || tok.is_empty() {
                -1
            } else {
                tok.parse::<i32>().unwrap_or(-1)
            };
            v.push(Gt { allele, phased });
        };
        for c in field.chars() {
            if c == '/' || c == '|' {
                push(&tok, phased, &mut v);
                tok.clear();
                phased = c == '|';
            } else {
                tok.push(c);
            }
        }
        push(&tok, phased, &mut v);
        ngts = ngts.max(v.len());
        out.push(v);
    }
    (out, ngts)
}

/// Consequence text fields (`vcsq_t`).
#[derive(Clone, Debug)]
pub struct Vcsq {
    pub strand: u8,
    pub ty: u32,
    pub trid: u32,
    pub vcf_ial: u32,
    pub biotype: u32,
    pub gene: usize,
    /// With `CSQ_PRINTED_UPSTREAM`: the record carrying the text.
    pub ref_: Option<(RecId, u32)>,
    /// `None` mirrors an unallocated kstring.
    pub vstr: Option<String>,
}

impl Vcsq {
    pub(crate) fn new(gff: &Gff, tr: usize, ty: u32, vcf_ial: u32) -> Self {
        let t = &gff.tscripts[tr];
        Self {
            strand: t.strand,
            ty,
            trid: t.id,
            vcf_ial,
            biotype: t.ty,
            gene: t.gene,
            ref_: None,
            vstr: None,
        }
    }
    fn vstr_len(&self) -> usize {
        self.vstr.as_ref().map_or(0, |s| s.len())
    }
}

/// A consequence tied to a haplotype (`csq_t`).
#[derive(Clone, Debug)]
pub struct Csq {
    pub pos: u32,
    /// `(vbuf id, index within the vbuf)` once pushed.
    pub vrec: Option<(usize, usize)>,
    /// 0-based index of the consequence at the record, for FORMAT/BCSQ.
    pub idx: usize,
    pub ty: Vcsq,
}

impl Csq {
    pub(crate) fn new(pos: u32, ty: Vcsq) -> Self {
        Self {
            pos,
            vrec: None,
            idx: 0,
            ty,
        }
    }
}

/// A buffered VCF record with the consequences attached to it.
pub struct Vrec {
    pub info: RecInfo,
    pub line: crate::vcf::VcfRecord,
    /// `nsmpl * nfmt_bcsq` bits, first/second haplotype interleaved.
    pub fmt_bm: Vec<u32>,
    pub nfmt: usize,
    pub vcsq: Vec<Vcsq>,
}

/// Records sharing one position.
pub struct Vbuf {
    pub pos: u32,
    pub vrec: Vec<Vrec>,
    /// Largest end of an open transcript overlapping these records.
    pub keep_until: u32,
    /// Slot in the bcftools-style ring of position buffers.
    pub slot: usize,
}

/// FIFO of position buffers; ids grow monotonically so `id - base` indexes
/// the deque.
#[derive(Default)]
pub struct VcfBuf {
    pub vbufs: VecDeque<Vbuf>,
    pub base: usize,
    pub pos2vbuf: FxHashMap<u32, usize>,
}

impl VcfBuf {
    pub fn get_mut(&mut self, id: usize) -> Option<&mut Vbuf> {
        id.checked_sub(self.base)
            .and_then(|i| self.vbufs.get_mut(i))
    }

    /// `csq_push`: merge the consequence into the record's list following the
    /// bcftools de-duplication rules. Returns `false` when it already existed.
    pub fn csq_push(&mut self, gff: &Gff, csq: &mut Csq, rec: RecId, chrom: &str) -> Result<bool> {
        let Some(&id) = self.pos2vbuf.get(&csq.pos) else {
            bail!(
                "This should not happen. {}:{}  {}",
                chrom,
                csq.pos + 1,
                csq.ty.vstr.as_deref().unwrap_or("")
            );
        };
        let vbuf = self.get_mut(id).expect("live vbuf");
        let Some(vi) = vbuf.vrec.iter().position(|v| v.info.id == rec) else {
            bail!(
                "This should not happen.. {}:{}  {}",
                chrom,
                csq.pos + 1,
                csq.ty.vstr.as_deref().unwrap_or("")
            );
        };
        let vrec = &mut vbuf.vrec[vi];
        if csq.ty.ty & CSQ_SPLICE_REGION != 0
            && csq.ty.ty & (CSQ_SPLICE_DONOR | CSQ_SPLICE_ACCEPTOR) != 0
        {
            csq.ty.ty &= !CSQ_SPLICE_REGION;
        }
        let gene_eq = |a: usize, b: usize| {
            a == b || (gff.genes[a].name.is_none() && gff.genes[b].name.is_none())
        };
        let mut dup: Option<usize> = None;
        if csq.ty.ty & CSQ_PRINTED_UPSTREAM != 0 {
            for (i, v) in vrec.vcsq.iter_mut().enumerate() {
                if csq.ty.ty & CSQ_START_STOP != 0 && v.ty & CSQ_START_STOP != 0 {
                    *v = csq.ty.clone();
                    dup = Some(i);
                    break;
                }
                if v.ty & CSQ_PRINTED_UPSTREAM == 0 {
                    continue;
                }
                if csq.ty.ref_.map(|r| r.0) != v.ref_.map(|r| r.0) {
                    continue;
                }
                dup = Some(i);
                break;
            }
        } else if csq.ty.ty & CSQ_COMPOUND != 0 {
            for (i, v) in vrec.vcsq.iter_mut().enumerate() {
                if csq.ty.trid != v.trid && (csq.ty.ty | v.ty) & CSQ_PRN_TSCRIPT != 0 {
                    continue;
                }
                if csq.ty.biotype != v.biotype
                    || !gene_eq(csq.ty.gene, v.gene)
                    || csq.ty.vcf_ial != v.vcf_ial
                {
                    continue;
                }
                if (csq.ty.ty & CSQ_UPSTREAM_STOP) ^ (v.ty & CSQ_UPSTREAM_STOP) != 0 {
                    continue;
                }
                if csq.ty.vstr.is_some() || v.vstr.is_some() {
                    if csq.ty.vstr.is_none() || v.vstr.is_none() {
                        if csq.ty.ty & CSQ_START_STOP != 0 && v.ty & CSQ_START_STOP != 0 {
                            v.ty |= csq.ty.ty;
                            if v.ty & CSQ_STOP_RETAINED != 0 {
                                v.ty &= !(CSQ_STOP_LOST | CSQ_SYNONYMOUS_VARIANT);
                            }
                            if v.vstr.is_none() {
                                v.vstr = csq.ty.vstr.clone();
                            }
                            dup = Some(i);
                            break;
                        }
                        continue;
                    }
                    if csq.ty.vstr != v.vstr {
                        continue;
                    }
                }
                v.ty |= csq.ty.ty;
                dup = Some(i);
                break;
            }
        } else {
            for (i, v) in vrec.vcsq.iter_mut().enumerate() {
                if csq.ty.trid != v.trid && (csq.ty.ty | v.ty) & CSQ_PRN_TSCRIPT != 0 {
                    continue;
                }
                if csq.ty.biotype != v.biotype {
                    continue;
                }
                if v.ty & CSQ_COMPOUND == 0 {
                    v.ty |= csq.ty.ty;
                    dup = Some(i);
                    break;
                }
                if v.ty == (v.ty | csq.ty.ty) {
                    dup = Some(i);
                    break;
                }
            }
        }
        csq.vrec = Some((id, vi));
        match dup {
            Some(i) => {
                csq.idx = i;
                Ok(false)
            }
            None => {
                csq.idx = vrec.vcsq.len();
                vrec.vcsq.push(csq.ty.clone());
                Ok(true)
            }
        }
    }
}

/// `kput_vcsq`: render one consequence. Start/stop flags of incomplete CDS
/// and missense accompanying a start/stop are dropped here.
pub fn kput_vcsq(gff: &Gff, csq: &mut Vcsq, out: &mut String) {
    if csq.ty & CSQ_INCOMPLETE_CDS != 0
        && csq.ty & !(CSQ_START_STOP | CSQ_INCOMPLETE_CDS | CSQ_UPSTREAM_STOP) != 0
    {
        csq.ty &= !(CSQ_START_STOP | CSQ_INCOMPLETE_CDS);
    }
    if csq.ty & CSQ_START_STOP != 0 && csq.ty & CSQ_MISSENSE_VARIANT != 0 {
        csq.ty &= !CSQ_MISSENSE_VARIANT;
    }
    if csq.ty & CSQ_PRINTED_UPSTREAM != 0
        && let Some((_, pos)) = csq.ref_ {
            out.push('@');
            out.push_str(&(pos + 1).to_string());
            return;
        }
    if csq.ty & CSQ_UPSTREAM_STOP != 0 {
        out.push('*');
    }
    let mut has_csq = false;
    let mut first = true;
    for (i, s) in CSQ_STRINGS.iter().enumerate().skip(1) {
        let Some(s) = s else { continue };
        if csq.ty & (1 << i) == 0 {
            continue;
        }
        has_csq = true;
        if !first {
            out.push('&');
        }
        first = false;
        out.push_str(s);
    }
    if csq.biotype == gff::GF_NMD && csq.ty & CSQ_PRN_NMD != 0 {
        if has_csq {
            out.push('&');
        }
        out.push_str("NMD_transcript");
    }
    out.push('|');
    if let Some(g) = gff.gene_name(csq.gene) {
        out.push_str(g);
    }
    out.push('|');
    if csq.ty & CSQ_PRN_TSCRIPT != 0 {
        out.push_str(gff.tscript_name(csq.trid));
    }
    out.push('|');
    out.push_str(gf_type2gff_string(csq.biotype));
    if csq_prn_strand(csq.ty) || csq.vstr_len() > 0 {
        out.push_str(if csq.strand == STRAND_FWD { "|+" } else { "|-" });
    }
    if let Some(v) = &csq.vstr {
        out.push_str(v);
    }
}

/// `kprint_aa_prediction`: the protein string, abbreviated to
/// `brief` residues (`-B`) with the end position appended.
pub(crate) fn kprint_aa_prediction(brief: usize, beg: usize, aa: &[u8], out: &mut String) {
    if brief == 0 || (aa.len() as i64) - (brief as i64) < 3 {
        out.push_str(&String::from_utf8_lossy(aa));
        return;
    }
    let mut len = aa.len();
    if aa[len - 1] == b'*' {
        len -= 1;
    }
    for &c in aa.iter().take(len.min(brief)) {
        out.push(c as char);
    }
    out.push_str("..");
    out.push_str(&(beg + len).to_string());
}

/// `cds_translate`: translate the variable part of a spliced transcript.
///
/// `sref`: spliced reference padded with `N_REF_PAD` on both ends; `seq`: the
/// part of the query transcript to translate (reference strand) and `seq_m`
/// the full query transcript length; `sbeg`: offset of `seq` in the query
/// transcript; `rbeg`/`rend`: offsets of the corresponding reference part;
/// `fill`: frameshift, keep translating to the end of the transcript.
pub fn cds_translate(
    sref: &[u8],
    seq: &[u8],
    seq_m: usize,
    sbeg: usize,
    rbeg: usize,
    rend: usize,
    strand: u8,
    fill: bool,
) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    if seq.is_empty() {
        out.push(b'?');
        return Ok(out);
    }
    let pad = N_REF_PAD;
    let mut tmp = [0u8; 3];
    let rf = |i: isize| -> u8 {
        if i < 0 {
            b'N'
        } else {
            sref.get(i as usize).copied().unwrap_or(b'N')
        }
    };
    if strand == STRAND_FWD {
        let npad = sbeg % 3;
        if npad > rbeg {
            bail!("cds_translate: npad {npad} > rbeg {rbeg}");
        }
        let mut i = 0usize;
        while i < npad {
            tmp[i] = rf((rbeg + i + pad) as isize - npad as isize);
            i += 1;
        }
        while i < 3 && i - npad < seq.len() {
            tmp[i] = seq[i - npad];
            i += 1;
        }
        let len = seq.len() + npad - i;
        let mut codon_ref = rend + pad;
        if i == 3 {
            out.push(dna2aa(&tmp));
            let mut codon = 3 - npad;
            let end = codon + len - 1 - (len % 3);
            while codon < end {
                out.push(dna2aa(&seq[codon..codon + 3]));
                codon += 3;
            }
            i = 0;
            while codon + i < seq.len() {
                tmp[i] = seq[codon + i];
                i += 1;
            }
        }
        if i > 0 {
            while i < 3 {
                tmp[i] = rf(codon_ref as isize);
                codon_ref += 1;
                i += 1;
            }
            out.push(dna2aa(&tmp));
        }
        if fill {
            let end = sref.len().saturating_sub(pad);
            while codon_ref + 3 <= end {
                out.push(dna2aa(&sref[codon_ref..codon_ref + 3]));
                codon_ref += 3;
            }
        }
    } else {
        let rem = seq_m as i64 - (sbeg + seq.len()) as i64;
        if rem < 0 {
            bail!(
                "cds_translate: sbeg={sbeg} seq.l={} seq.m={seq_m}",
                seq.len()
            );
        }
        let npad = (rem % 3) as usize;
        let mut i: i32;
        if npad == 2 {
            tmp[1] = rf((rend + pad) as isize);
            tmp[2] = rf((rend + pad + 1) as isize);
            i = 0;
        } else if npad == 1 {
            tmp[2] = rf((rend + pad) as isize);
            i = 1;
        } else {
            i = 2;
        }
        let mut end = seq.len();
        while i >= 0 && end > 0 {
            end -= 1;
            tmp[i as usize] = seq[end];
            i -= 1;
        }
        if i == -1 {
            out.push(cdna2aa(&tmp));
            let mut codon = end as isize - 3;
            while codon >= 0 {
                out.push(cdna2aa(&seq[codon as usize..codon as usize + 3]));
                codon -= 3;
            }
            if -codon == 2 {
                tmp[2] = seq[0];
                i = 1;
            } else if -codon == 1 {
                tmp[1] = seq[0];
                tmp[2] = seq[1];
                i = 0;
            } else {
                i = -1;
            }
        }
        let mut end_ref = (pad + rbeg) as isize;
        if i >= 0 {
            while i >= 0 && end_ref >= 0 {
                end_ref -= 1;
                tmp[i as usize] = rf(end_ref);
                i -= 1;
            }
            out.push(cdna2aa(&tmp));
        }
        if fill {
            let mut codon = end_ref - 3;
            while codon >= pad as isize {
                out.push(cdna2aa(&sref[codon as usize..codon as usize + 3]));
                codon -= 3;
            }
        }
    }
    Ok(out)
}

/// `len` bytes of `s` from `beg`; `N` where the request leaves the buffer
/// (the C code would read padding or adjacent memory).
pub(crate) fn sub(s: &[u8], beg: i64, len: i64) -> Vec<u8> {
    let mut out = Vec::with_capacity(len.max(0) as usize);
    for i in 0..len.max(0) {
        let p = beg + i;
        out.push(if p < 0 {
            b'N'
        } else {
            s.get(p as usize).copied().unwrap_or(b'N')
        });
    }
    out
}

/// Byte `i` of a C string: NUL past either end.
pub(crate) fn c_at(s: &[u8], i: i64) -> u8 {
    if i < 0 {
        0
    } else {
        s.get(i as usize).copied().unwrap_or(0)
    }
}

/// Truncate a translation at its first stop codon; true when one was found.
pub(crate) fn truncate_at_stop(aa: &mut Vec<u8>) -> bool {
    if let Some(i) = aa.iter().position(|&c| c == b'*') {
        aa.truncate(i + 1);
        true
    } else {
        false
    }
}

#[cfg(test)]
#[path = "../../tests/unit/csq_engine.rs"]
mod tests;
