//! Driver: buffers records, opens transcripts, runs the per-record tests and
//! flushes annotated records (`process`, `test_*`, `hap_flush`, `vbuf_flush`
//! in csq.c).

use std::io::Write;

use anyhow::{Result, bail};

use super::gff::{Gff, Prime, STRAND_FWD, gf_is_coding};
use super::hap::{HapNode, HapState, TscriptAux, hap_finalize, hap_init, tscript_splice_ref};
use super::splice::{SPLICE_INSIDE, SPLICE_OVERLAP, Splice, splice_csq};
use super::*;
use crate::fasta::IndexedFasta;
use crate::vcf::{VcfRecord, VcfSink};

/// `--phase`: how unphased heterozygous genotypes are handled.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Phase {
    /// `r`: error on an unphased het.
    Require,
    /// `m`: merge all genotypes into one haplotype.
    Merge,
    /// `a`: take genotypes as they are.
    AsIs,
    /// `s`: skip unphased hets.
    Skip,
    /// `R`: build non-reference haplotypes when possible.
    NonRef,
    /// `-s -`: ignore samples, apply every variant.
    DropGt,
}

pub struct CsqOpts {
    pub phase: Phase,
    /// `2 * --ncsq`.
    pub ncsq2_max: usize,
    pub nfmt_bcsq: usize,
    pub local_csq: bool,
    /// `-B`; 0 prints the full protein change.
    pub brief_predictions: usize,
    pub unify_chr_names: bool,
    pub verbosity: u8,
    pub bcsq_tag: String,
    /// Header samples, 0 with `--phase`-less `-s -`.
    pub hdr_nsmpl: usize,
    pub fa_fname: String,
}

pub enum CsqOutput {
    Vcf(VcfSink),
    Text(Box<dyn Write>),
}

impl CsqOutput {
    fn is_text(&self) -> bool {
        matches!(self, Self::Text(_))
    }
    fn text(&mut self, s: &str) -> Result<()> {
        if let Self::Text(w) = self {
            w.write_all(s.as_bytes())?;
        }
        Ok(())
    }
}

#[derive(Default)]
pub(crate) struct Warn {
    pub overlaps_warned: bool,
    pub multiploid_warned: bool,
    pub ncsq2_small_warned: usize,
    pub symbolic_warned: bool,
    pub small_ref_padding_warned: bool,
}

/// Borrowed working set handed to the staging functions.
pub(crate) struct Ctx<'a> {
    pub gff: &'a Gff,
    pub opts: &'a CsqOpts,
    pub buf: &'a mut VcfBuf,
    /// Header indexes of the samples to process.
    pub smpl: &'a [usize],
    pub samples: &'a [String],
    pub warn: &'a mut Warn,
    pub out: &'a mut CsqOutput,
    pub chrom: &'a str,
}

fn drop_chr(unify: bool, chr: &str) -> &str {
    if unify && chr.len() >= 3 && chr[..3].eq_ignore_ascii_case("chr") {
        &chr[3..]
    } else {
        chr
    }
}

impl Ctx<'_> {
    pub fn drop_chr_prefix(&self, chr: &str) -> String {
        drop_chr(self.opts.unify_chr_names, chr).to_string()
    }

    pub fn csq_push(&mut self, csq: &mut Csq, rec: RecId) -> Result<bool> {
        self.buf.csq_push(self.gff, csq, rec, self.chrom)
    }

    fn csq_print_text(&mut self, csq: &mut Csq, ismpl: Option<usize>, ihap: usize) -> Result<()> {
        if csq.ty.ty & CSQ_PRINTED_UPSTREAM != 0 {
            return Ok(());
        }
        let smpl = ismpl.map(|i| self.samples[i].as_str()).unwrap_or("-");
        let mut s = String::new();
        kput_vcsq(self.gff, &mut csq.ty, &mut s);
        let hap = if ihap > 0 {
            ihap.to_string()
        } else {
            "-".to_string()
        };
        let line = format!(
            "CSQ\t{}\t{}\t{}\t{}\t{}\n",
            smpl,
            hap,
            self.chrom,
            csq.pos + 1,
            s
        );
        self.out.text(&line)
    }

    /// `csq_stage`: register the consequence on its record and set the
    /// FORMAT bits of every haplotype carrying the allele.
    pub fn csq_stage(&mut self, csq: &mut Csq, rec: &RecInfo) -> Result<()> {
        let added = self.csq_push(csq, rec.id)?;
        if !added && self.opts.phase == Phase::DropGt {
            return Ok(());
        }
        let ngt = if self.opts.phase != Phase::DropGt {
            rec.ngts
        } else {
            0
        };
        if ngt == 0 {
            if self.out.is_text() {
                self.csq_print_text(csq, None, 0)?;
            }
            return Ok(());
        }
        if ngt > 2 {
            bail!(
                "csq: genotypes with ploidy {} are not supported at {}:{}",
                ngt,
                rec.chrom,
                rec.pos + 1
            );
        }
        if self.out.is_text() {
            for i in 0..self.smpl.len() {
                let ismpl = self.smpl[i];
                for j in 0..ngt {
                    let Some(g) = rec.gt(ismpl, j) else { continue };
                    if g.allele < 0 {
                        continue;
                    }
                    let ial = g.allele as u32;
                    if ial == 0 || ial != csq.ty.vcf_ial {
                        continue;
                    }
                    self.csq_print_text(csq, Some(ismpl), j + 1)?;
                }
            }
            return Ok(());
        }
        let Some((vb, vi)) = csq.vrec else {
            return Ok(());
        };
        let nfmt_bcsq = self.opts.nfmt_bcsq;
        let ncsq2_max = self.opts.ncsq2_max;
        let verbosity = self.opts.verbosity;
        let vrec = &mut self.buf.get_mut(vb).expect("live vbuf").vrec[vi];
        for i in 0..self.smpl.len() {
            let ismpl = self.smpl[i];
            for j in 0..ngt {
                let Some(g) = rec.gt(ismpl, j) else { continue };
                if g.allele < 0 {
                    continue;
                }
                let ial = g.allele as u32;
                if ial == 0 || ial != csq.ty.vcf_ial {
                    continue;
                }
                let icsq2 = 2 * csq.idx + j;
                if icsq2 >= ncsq2_max {
                    if verbosity > 0 && (self.warn.ncsq2_small_warned == 0 || verbosity > 1) {
                        eprintln!(
                            "Warning: Too many consequences for sample {} at {}:{}, keeping the first {} and skipping the rest.",
                            self.samples[ismpl],
                            self.chrom,
                            vrec.info.pos + 1,
                            icsq2 + 1
                        );
                        if self.warn.ncsq2_small_warned == 0 {
                            eprintln!(
                                "         The limit can be increased by setting the --ncsq parameter. This warning is printed only once."
                            );
                        }
                        self.warn.ncsq2_small_warned = 1;
                    }
                    if self.warn.ncsq2_small_warned < icsq2 {
                        self.warn.ncsq2_small_warned = icsq2;
                    }
                    break;
                }
                let (ival, ibit) = icsq2_to_bit(icsq2);
                if vrec.nfmt < 1 + ival {
                    vrec.nfmt = 1 + ival;
                }
                vrec.fmt_bm[ismpl * nfmt_bcsq + ival] |= 1 << ibit;
            }
        }
        Ok(())
    }
}

/// `hap_stage_vcf`: FORMAT bits for every consequence on the haplotype
/// ending at `node`.
fn hap_stage_vcf(ctx: &mut Ctx, aux: &TscriptAux, ismpl: usize, ihap: usize, node: Option<usize>) {
    let Some(node) = node else { return };
    let nfmt_bcsq = ctx.opts.nfmt_bcsq;
    for csq in &aux.nodes[node].csq_list {
        let Some((vb, vi)) = csq.vrec else { continue };
        let icsq2 = 2 * csq.idx + ihap;
        let vrec = &mut ctx.buf.get_mut(vb).expect("live vbuf").vrec[vi];
        if icsq2 >= ctx.opts.ncsq2_max {
            if ctx.opts.verbosity > 0
                && (ctx.warn.ncsq2_small_warned == 0 || ctx.opts.verbosity > 1)
            {
                eprintln!(
                    "Warning: Too many consequences for sample {} at {}:{}, keeping the first {} and skipping the rest.",
                    ctx.samples[ismpl],
                    ctx.chrom,
                    vrec.info.pos + 1,
                    csq.idx
                );
                if ctx.warn.ncsq2_small_warned == 0 {
                    eprintln!(
                        "         The limit can be increased by setting the --ncsq parameter. This warning is printed only once."
                    );
                }
            }
            if ctx.warn.ncsq2_small_warned < icsq2 {
                ctx.warn.ncsq2_small_warned = icsq2;
            }
            break;
        }
        let (ival, ibit) = icsq2_to_bit(icsq2);
        if vrec.nfmt < 1 + ival {
            vrec.nfmt = 1 + ival;
        }
        vrec.fmt_bm[ismpl * nfmt_bcsq + ival] |= 1 << ibit;
    }
}

fn hap_print_text(
    ctx: &mut Ctx,
    aux: &mut TscriptAux,
    ismpl: Option<usize>,
    ihap: usize,
    node: Option<usize>,
) -> Result<()> {
    let Some(node) = node else { return Ok(()) };
    let smpl = ismpl
        .map(|i| ctx.samples[i].clone())
        .unwrap_or_else(|| "-".to_string());
    let hap = if ihap > 0 {
        ihap.to_string()
    } else {
        "-".to_string()
    };
    for csq in &mut aux.nodes[node].csq_list {
        if csq.ty.ty & CSQ_PRINTED_UPSTREAM != 0 {
            continue;
        }
        let mut s = String::new();
        kput_vcsq(ctx.gff, &mut csq.ty, &mut s);
        let line = format!(
            "CSQ\t{}\t{}\t{}\t{}\t{}\n",
            smpl,
            hap,
            ctx.chrom,
            csq.pos + 1,
            s
        );
        ctx.out.text(&line)?;
    }
    Ok(())
}

/// `tscript_init_ref`: the transcript sequence with `N_REF_PAD` flanks,
/// `N`-padded where the contig ends.
fn tscript_init_ref(
    fasta: &IndexedFasta,
    opts: &CsqOpts,
    gff: &Gff,
    tr_i: usize,
    chr: &str,
) -> Result<Vec<u8>> {
    let tr = &gff.tscripts[tr_i];
    let pad = N_REF_PAD as i64;
    let pad_beg = if tr.beg as i64 >= pad {
        pad
    } else {
        tr.beg as i64
    };
    let mut name = chr.to_string();
    if !fasta.has(&name) {
        name = drop_chr(opts.unify_chr_names, chr).to_string();
        if !fasta.has(&name) {
            name = if opts.unify_chr_names {
                format!("chr{chr}")
            } else {
                chr.to_string()
            };
        }
    }
    let Some(clen) = fasta.length(&name) else {
        bail!(
            "faidx_fetch_seq failed {}:{}-{}",
            chr,
            tr.beg + 1,
            tr.end + 1
        );
    };
    let beg0 = tr.beg as i64 - pad_beg;
    let mut end0 = tr.end as i64 + pad;
    if end0 >= clen as i64 {
        end0 = clen as i64 - 1;
    }
    if beg0 > end0 || beg0 < 0 {
        bail!(
            "faidx_fetch_seq failed {}:{}-{}",
            chr,
            tr.beg + 1,
            tr.end + 1
        );
    }
    let len = (end0 - beg0 + 1) as usize;
    let Some(seq) = fasta.slice(&name, beg0 as u32 + 1, len) else {
        bail!(
            "faidx_fetch_seq failed {}:{}-{}",
            chr,
            tr.beg + 1,
            tr.end + 1
        );
    };
    let tr_len = tr.end as i64 - tr.beg as i64 + 1;
    let pad_end = len as i64 - (tr_len + pad_beg);
    if pad_beg + pad_end == 2 * pad {
        return Ok(seq.to_vec());
    }
    let mut out = Vec::with_capacity((tr_len + 2 * pad) as usize);
    out.extend(std::iter::repeat_n(b'N', (pad - pad_beg).max(0) as usize));
    out.extend_from_slice(seq);
    out.extend(std::iter::repeat_n(b'N', (pad - pad_end).max(0) as usize));
    Ok(out)
}

/// `sanity_check_ref`: the VCF REF allele must match the reference.
fn sanity_check_ref(aux: &TscriptAux, tr: &gff::Tscript, rec: &RecInfo) -> Result<()> {
    let mut rbeg = rec.pos as i64 - tr.beg as i64 + N_REF_PAD as i64;
    let mut vbeg = 0i64;
    if rbeg < 0 {
        vbeg += -rbeg;
        rbeg = 0;
    }
    let vcf = rec.alleles[0].as_bytes();
    let mut i = 0usize;
    loop {
        let r = c_at(&aux.ref_, rbeg + i as i64);
        let v = c_at(vcf, vbeg + i as i64);
        if r == 0 || v == 0 {
            break;
        }
        if r != v && !r.eq_ignore_ascii_case(&v) {
            bail!(
                "Error: the fasta reference does not match the VCF REF allele at {}:{} .. fasta={} vcf={}",
                rec.chrom,
                rec.pos as i64 + vbeg + 1,
                r as char,
                v as char
            );
        }
        i += 1;
    }
    Ok(())
}

/// Min-heap of open transcripts keyed by end position, with the exact
/// sift-up/sift-down of bcftools `kheap.h` so ties flush in the same order.
#[derive(Default)]
struct TrHeap {
    dat: Vec<(u32, usize)>,
}

impl TrHeap {
    fn is_empty(&self) -> bool {
        self.dat.is_empty()
    }
    fn top(&self) -> Option<(u32, usize)> {
        self.dat.first().copied()
    }
    fn insert(&mut self, end: u32, tr: usize) {
        self.dat.push((end, tr));
        let mut i = self.dat.len() - 1;
        while i > 0 && end < self.dat[(i - 1) / 2].0 {
            self.dat[i] = self.dat[(i - 1) / 2];
            i = (i - 1) / 2;
        }
        self.dat[i] = (end, tr);
    }
    fn delete(&mut self) {
        let Some(last) = self.dat.pop() else { return };
        if self.dat.is_empty() {
            return;
        }
        self.dat[0] = last;
        let mut i = 0usize;
        loop {
            let n = self.dat.len();
            let (l, r) = (2 * i + 1, 2 * i + 2);
            let mut extreme = if l < n && self.dat[l].0 < self.dat[i].0 {
                l
            } else {
                i
            };
            if r < n && self.dat[r].0 < self.dat[extreme].0 {
                extreme = r;
            }
            if extreme == i {
                break;
            }
            self.dat.swap(i, extreme);
            i = extreme;
        }
    }
}

/// bcftools keeps `vrec_t` objects in a ring of position buffers and never
/// resets their `nfmt`, so a record that gains consequences but no sample
/// bits still carries a zero FORMAT/BCSQ when its slot was used before.
/// This mirrors that ring (`rbuf.h`) closely enough to reproduce the output.
#[derive(Default)]
struct VrecRing {
    /// `nfmt` remembered per vrec index of each slot's `vbuf_t`; `None` is
    /// a slot bcftools has not allocated yet.
    slots: Vec<Option<Vec<usize>>>,
    m: usize,
    n: usize,
    f: usize,
}

impl VrecRing {
    /// `rbuf_expand0` + `rbuf_append`: the slot of a new position buffer.
    fn append(&mut self) -> usize {
        let need = self.n + 1;
        if need > self.m {
            let mut m = need + self.f;
            m -= 1;
            m |= m >> 1;
            m |= m >> 2;
            m |= m >> 4;
            m |= m >> 8;
            m |= m >> 16;
            m += 1;
            self.slots.resize(m, None);
            for k in 0..self.f {
                self.slots[self.m + k] = self.slots[k].take();
            }
            self.m = m;
        }
        self.n += 1;
        let i = self.f + self.n;
        let slot = if i <= self.m { i - 1 } else { i - self.m - 1 };
        if self.slots[slot].is_none() {
            self.slots[slot] = Some(Vec::new());
        }
        slot
    }
    fn stale_nfmt(&self, slot: usize, k: usize) -> usize {
        self.slots[slot]
            .as_ref()
            .and_then(|v| v.get(k))
            .copied()
            .unwrap_or(0)
    }
    fn remember(&mut self, slot: usize, k: usize, nfmt: usize) {
        let v = self.slots[slot].get_or_insert_with(Vec::new);
        if v.len() <= k {
            v.resize(k + 1, 0);
        }
        v[k] = nfmt;
    }
    /// `rbuf_shift`.
    fn shift(&mut self) {
        if self.n == 0 {
            return;
        }
        self.f += 1;
        if self.f >= self.m {
            self.f = 0;
        }
        self.n -= 1;
    }
}

/// The consequence caller: feed sorted records with [`Caller::process`],
/// then [`Caller::finish`].
pub struct Caller {
    gff: Gff,
    fasta: IndexedFasta,
    opts: CsqOpts,
    smpl: Vec<usize>,
    samples: Vec<String>,
    buf: VcfBuf,
    ring: VrecRing,
    aux: Vec<Option<TscriptAux>>,
    /// Open transcripts by end position.
    active_tr: TrHeap,
    hap: HapState,
    rm_tr: Vec<usize>,
    /// Chromosome of the buffered records.
    rid: Option<String>,
    /// First record of the current chromosome, for the sortedness check.
    prev: Option<(String, u32)>,
    next_id: RecId,
    warn: Warn,
    out: CsqOutput,
}

impl Caller {
    pub fn new(
        gff: Gff,
        fasta: IndexedFasta,
        opts: CsqOpts,
        smpl: Vec<usize>,
        samples: Vec<String>,
        out: CsqOutput,
    ) -> Self {
        let ntr = gff.tscripts.len();
        let mut aux = Vec::with_capacity(ntr);
        aux.resize_with(ntr, || None);
        Self {
            gff,
            fasta,
            opts,
            smpl,
            samples,
            buf: VcfBuf::default(),
            ring: VrecRing::default(),
            aux,
            active_tr: TrHeap::default(),
            hap: HapState::default(),
            rm_tr: Vec::new(),
            rid: None,
            prev: None,
            next_id: 0,
            warn: Warn::default(),
            out,
        }
    }

    pub fn into_output(self) -> CsqOutput {
        self.out
    }

    fn rec_info(&mut self, line: &VcfRecord) -> RecInfo {
        let mut alleles = vec![line.ref_allele.clone()];
        if line.alt != "." && !line.alt.is_empty() {
            alleles.extend(line.alt.split(',').map(str::to_string));
        }
        let (gts, ngts) = if self.opts.phase != Phase::DropGt && self.opts.hdr_nsmpl > 0 {
            parse_genotypes(line.format.as_deref(), &line.samples)
        } else {
            (Vec::new(), 0)
        };
        let id = self.next_id;
        self.next_id += 1;
        RecInfo {
            id,
            chrom: line.chrom.clone(),
            pos: line.pos.saturating_sub(1),
            rlen: line.ref_allele.len() as u32,
            alleles: Rc::new(alleles),
            gts: Rc::new(gts),
            ngts,
        }
    }

    /// `process`: annotate one record. `filter_pass` is the result of
    /// `-i`/`-e`; a failing record is written through unannotated.
    pub fn process(&mut self, line: VcfRecord, filter_pass: Option<bool>) -> Result<()> {
        let rec = self.rec_info(&line);
        if self.prev.as_ref().map(|p| p.0.as_str()) != Some(rec.chrom.as_str()) {
            self.prev = Some((rec.chrom.clone(), rec.pos));
            let c = &rec.chrom;
            let dropped = drop_chr(self.opts.unify_chr_names, c);
            let added = if self.opts.unify_chr_names {
                format!("chr{c}")
            } else {
                c.clone()
            };
            if !self.fasta.has(c) && !self.fasta.has(dropped) && !self.fasta.has(&added) {
                bail!(
                    "Error: the chromosome \"{}\" is not present in {}",
                    c,
                    self.opts.fa_fname
                );
            }
        }
        if let Some((_, ppos)) = &self.prev
            && *ppos > rec.pos
        {
            bail!(
                "Error: The file is not sorted, {}:{} comes before {}:{}",
                rec.chrom,
                ppos + 1,
                rec.chrom,
                rec.pos + 1
            );
        }
        let mut call_csq = true;
        if rec.n_allele() < 2 {
            call_csq = false;
        } else if rec.n_allele() == 2 && {
            let a = rec.alleles[1].as_bytes();
            a.first() == Some(&b'*') || a.get(1) == Some(&b'*')
        } {
            call_csq = false;
        } else if rec.alleles[1].starts_with('<') && !rec.alleles[1].starts_with("<INS") {
            call_csq = false;
        }
        if call_csq && let Some(p) = filter_pass {
            call_csq = p;
        }
        if !call_csq {
            if self.out.is_text() {
                return Ok(());
            }
            self.vbuf_push(rec.clone(), line);
            let p = rec.pos.wrapping_sub(1);
            self.hap_flush(p)?;
            self.vbuf_flush(p)?;
            return Ok(());
        }
        if self.rid.as_deref() != Some(rec.chrom.as_str()) {
            self.hap_flush(u32::MAX)?;
            self.vbuf_flush(u32::MAX)?;
        }
        self.rid = Some(rec.chrom.clone());
        let vb = self.vbuf_push(rec.clone(), line);
        if !rec.alleles[1].starts_with('<') {
            let mut hit = if self.opts.local_csq {
                self.test_cds_local(&rec)?
            } else {
                self.test_cds(&rec, vb)?
            };
            hit |= self.test_utr(&rec)?;
            hit |= self.test_splice(&rec)?;
            if !hit {
                self.test_tscript(&rec)?;
            }
        } else {
            self.test_symbolic_alt(&rec)?;
        }
        if rec.pos > 0 {
            self.hap_flush(rec.pos - 1)?;
            self.vbuf_flush(rec.pos - 1)?;
        }
        Ok(())
    }

    /// Flush everything; call once after the last record.
    pub fn finish(&mut self) -> Result<()> {
        self.hap_flush(u32::MAX)?;
        self.vbuf_flush(u32::MAX)?;
        if self.warn.ncsq2_small_warned > 0 {
            eprintln!(
                "Note: Some samples had too many consequences to be represented in {} bytes. If you need to record them all,\n      the limit can be increased by running with `--ncsq {}`.",
                ncsq2_to_nfmt(self.opts.ncsq2_max) / 8,
                1 + self.warn.ncsq2_small_warned / 2
            );
        }
        if let CsqOutput::Text(w) = &mut self.out {
            w.flush()?;
        }
        Ok(())
    }

    fn ctx(&mut self) -> Ctx<'_> {
        Ctx {
            gff: &self.gff,
            opts: &self.opts,
            buf: &mut self.buf,
            smpl: &self.smpl,
            samples: &self.samples,
            warn: &mut self.warn,
            out: &mut self.out,
            chrom: self.rid.as_deref().unwrap_or(""),
        }
    }

    fn warn_overlap(
        &mut self,
        chr: &str,
        rec: &RecInfo,
        sample: Option<usize>,
        ial: usize,
    ) -> Result<()> {
        let smpl = sample
            .map(|s| format!(", sample {}", self.samples[s]))
            .unwrap_or_default();
        if self.opts.verbosity > 0 && (!self.warn.overlaps_warned || self.opts.verbosity > 1) {
            eprintln!(
                "Warning: Skipping overlapping variants at {}:{}{}\t{}>{}.",
                chr,
                rec.pos + 1,
                smpl,
                rec.alleles[0],
                rec.alleles[ial]
            );
            if !self.warn.overlaps_warned {
                eprintln!(
                    "         This message is printed only once, the verbosity can be increased with `--verbose 2`"
                );
            }
            self.warn.overlaps_warned = true;
        }
        let line = format!(
            "LOG\tWarning: Skipping overlapping variants at {}:{}{}\t{}>{}\n",
            chr,
            rec.pos + 1,
            smpl,
            rec.alleles[0],
            rec.alleles[ial]
        );
        self.out.text(&line)
    }

    /// `test_cds`: extend the haplotype trees of the coding transcripts
    /// whose CDS overlap the record.
    fn test_cds(&mut self, rec: &RecInfo, vb: usize) -> Result<bool> {
        let chr = drop_chr(self.opts.unify_chr_names, &rec.chrom).to_string();
        let hits = self
            .gff
            .idx_cds
            .overlap(&chr, rec.pos, rec.pos.wrapping_add(rec.rlen));
        let mut ret = false;
        for ci in hits {
            let tr_i = self.gff.cds[ci].tr;
            if !gf_is_coding(self.gff.tscripts[tr_i].ty) {
                continue;
            }
            let tr_end = self.gff.tscripts[tr_i].end;
            if let Some(vbuf) = self.buf.get_mut(vb)
                && vbuf.keep_until < tr_end
            {
                vbuf.keep_until = tr_end;
            }
            ret = true;
            if self.aux[tr_i].is_none() {
                let r = tscript_init_ref(&self.fasta, &self.opts, &self.gff, tr_i, &chr)?;
                let nhap = if self.opts.phase == Phase::DropGt {
                    1
                } else {
                    2 * self.smpl.len()
                };
                let root = HapNode {
                    nend: nhap as u32,
                    ty: HAP_ROOT,
                    ..Default::default()
                };
                self.aux[tr_i] = Some(TscriptAux {
                    ref_: r,
                    sref: Vec::new(),
                    nodes: vec![root],
                    hap: vec![None; nhap],
                });
                self.active_tr.insert(tr_end, tr_i);
            }
            sanity_check_ref(
                self.aux[tr_i].as_ref().unwrap(),
                &self.gff.tscripts[tr_i],
                rec,
            )?;
            let phase = self.opts.phase;
            if phase == Phase::DropGt {
                if rec.allele_is_symbolic(1) {
                    continue;
                }
                let mut ctx = Ctx {
                    gff: &self.gff,
                    opts: &self.opts,
                    buf: &mut self.buf,
                    smpl: &self.smpl,
                    samples: &self.samples,
                    warn: &mut self.warn,
                    out: &mut self.out,
                    chrom: self.rid.as_deref().unwrap_or(""),
                };
                let aux = self.aux[tr_i].as_mut().unwrap();
                let parent = aux.hap[0].unwrap_or(0);
                let mut child = HapNode::default();
                let hap_ret =
                    hap_init(&mut ctx, &*aux, &aux.nodes[parent], &mut child, ci, rec, 1)?;
                if hap_ret != 0 {
                    if hap_ret == 1 {
                        self.warn_overlap(&chr, rec, None, 1)?;
                    }
                    continue;
                }
                if child.ty == HAP_SSS {
                    let mut csq = Csq::new(rec.pos, Vcsq::new(&self.gff, tr_i, child.csq, 1));
                    ctx.csq_stage(&mut csq, rec)?;
                    continue;
                }
                aux.nodes[parent].nend -= 1;
                let idx = aux.nodes.len();
                aux.nodes.push(child);
                aux.nodes[parent].children = vec![idx];
                aux.hap[0] = Some(idx);
                aux.nodes[idx].nend = 1;
                continue;
            }
            let ngts = rec.ngts;
            if ngts != 1 && ngts != 2 {
                if self.opts.verbosity > 0
                    && (!self.warn.multiploid_warned || self.opts.verbosity > 1)
                {
                    eprintln!(
                        "Warning: Skipping site with non-diploid/non-haploid genotypes at {}:{}\t{}>{}.",
                        chr,
                        rec.pos + 1,
                        rec.alleles[0],
                        rec.alleles[1]
                    );
                    if !self.warn.multiploid_warned {
                        eprintln!(
                            "         This message is printed only once, the verbosity can be increased with `--verbose 2`"
                        );
                    }
                    self.warn.multiploid_warned = true;
                }
                let line = format!(
                    "LOG\tWarning: Skipping site with non-diploid/non-haploid genotypes at {}:{}\t{}>{}\n",
                    chr,
                    rec.pos + 1,
                    rec.alleles[0],
                    rec.alleles[1]
                );
                self.out.text(&line)?;
                continue;
            }
            for ismpl_i in 0..self.smpl.len() {
                let ismpl = self.smpl[ismpl_i];
                let Some(g0) = rec.gt(ismpl, 0) else { continue };
                if g0.allele < 0 {
                    continue;
                }
                let mut gt = [Some(g0), rec.gt(ismpl, 1)];
                if ngts > 1
                    && let Some(g1) = gt[1]
                    && g1.allele >= 0
                    && g0.allele != g1.allele
                {
                    if phase == Phase::Merge && g0.allele == 0 {
                        gt[0] = Some(g1);
                    }
                    if !g0.phased && !g1.phased {
                        match phase {
                            Phase::Require => bail!(
                                "Unphased heterozygous genotype at {}:{}, sample {}. See the --phase option.",
                                chr,
                                rec.pos + 1,
                                self.samples[ismpl]
                            ),
                            Phase::Skip => continue,
                            Phase::NonRef => {
                                if gt[0].unwrap().allele == 0 {
                                    gt[0] = gt[1];
                                } else if g1.allele == 0 {
                                    gt[1] = gt[0];
                                }
                            }
                            _ => {}
                        }
                    }
                }
                for ihap in 0..ngts {
                    let Some(g) = gt[ihap] else { continue };
                    if g.allele < 0 {
                        continue;
                    }
                    let i = 2 * ismpl_i + ihap;
                    let ial = g.allele as usize;
                    if ial == 0 {
                        continue;
                    }
                    if ial >= rec.n_allele() {
                        bail!(
                            "csq: genotype allele {} out of range at {}:{}",
                            ial,
                            chr,
                            rec.pos + 1
                        );
                    }
                    if rec.allele_is_symbolic(ial) {
                        continue;
                    }
                    let mut ctx = Ctx {
                        gff: &self.gff,
                        opts: &self.opts,
                        buf: &mut self.buf,
                        smpl: &self.smpl,
                        samples: &self.samples,
                        warn: &mut self.warn,
                        out: &mut self.out,
                        chrom: self.rid.as_deref().unwrap_or(""),
                    };
                    let aux = self.aux[tr_i].as_mut().unwrap();
                    let parent = aux.hap[i].unwrap_or(0);
                    if aux.nodes[parent].cur_rec == Some(rec.id)
                        && aux.nodes[parent].cur_child.get(ial).copied().unwrap_or(-1) >= 0
                    {
                        let c =
                            aux.nodes[parent].children[aux.nodes[parent].cur_child[ial] as usize];
                        aux.hap[i] = Some(c);
                        aux.nodes[c].nend += 1;
                        aux.nodes[parent].nend -= 1;
                        continue;
                    }
                    let mut child = HapNode::default();
                    let hap_ret = hap_init(
                        &mut ctx,
                        &*aux,
                        &aux.nodes[parent],
                        &mut child,
                        ci,
                        rec,
                        ial,
                    )?;
                    if hap_ret != 0 {
                        if hap_ret == 1 {
                            self.warn_overlap(&chr, rec, Some(ismpl), ial)?;
                        }
                        continue;
                    }
                    if child.ty == HAP_SSS {
                        let mut csq =
                            Csq::new(rec.pos, Vcsq::new(&self.gff, tr_i, child.csq, ial as u32));
                        ctx.csq_stage(&mut csq, rec)?;
                        continue;
                    }
                    if aux.nodes[parent].cur_rec != Some(rec.id) {
                        aux.nodes[parent].cur_child = vec![-1; rec.n_allele()];
                        aux.nodes[parent].cur_rec = Some(rec.id);
                    }
                    let j = aux.nodes[parent].children.len();
                    let idx = aux.nodes.len();
                    aux.nodes.push(child);
                    aux.nodes[parent].cur_child[ial] = j as i32;
                    aux.nodes[parent].children.push(idx);
                    aux.hap[i] = Some(idx);
                    aux.nodes[idx].nend += 1;
                    aux.nodes[parent].nend -= 1;
                }
            }
        }
        Ok(ret)
    }

    /// `test_cds_local` (`-l`): each record on its own, no haplotypes.
    fn test_cds_local(&mut self, rec: &RecInfo) -> Result<bool> {
        let chr = drop_chr(self.opts.unify_chr_names, &rec.chrom).to_string();
        let hits = self
            .gff
            .idx_cds
            .overlap(&chr, rec.pos, rec.pos.wrapping_add(rec.rlen));
        let mut ret = false;
        let root = HapNode {
            ty: HAP_ROOT,
            ..Default::default()
        };
        let pad = N_REF_PAD;
        for ci in hits {
            let tr_i = self.gff.cds[ci].tr;
            let tr = &self.gff.tscripts[tr_i];
            if !gf_is_coding(tr.ty) {
                continue;
            }
            ret = true;
            if self.aux[tr_i].is_none() {
                let r = tscript_init_ref(&self.fasta, &self.opts, &self.gff, tr_i, &chr)?;
                let mut aux = TscriptAux {
                    ref_: r,
                    sref: Vec::new(),
                    nodes: Vec::new(),
                    hap: Vec::new(),
                };
                tscript_splice_ref(&mut aux, tr, &self.gff.cds);
                self.aux[tr_i] = Some(aux);
                self.active_tr.insert(tr.end, tr_i);
            }
            let aux = self.aux[tr_i].as_ref().unwrap();
            sanity_check_ref(aux, tr, rec)?;
            let strand = tr.strand;
            let nsref = aux.sref.len() as i64;
            let mut ctx = Ctx {
                gff: &self.gff,
                opts: &self.opts,
                buf: &mut self.buf,
                smpl: &self.smpl,
                samples: &self.samples,
                warn: &mut self.warn,
                out: &mut self.out,
                chrom: self.rid.as_deref().unwrap_or(""),
            };
            for ial in 1..rec.n_allele() {
                if rec.allele_is_symbolic(ial) {
                    continue;
                }
                let mut node = HapNode::default();
                if hap_init(&mut ctx, aux, &root, &mut node, ci, rec, ial)? != 0 {
                    continue;
                }
                let mut csq = Csq::new(rec.pos, Vcsq::new(&self.gff, tr_i, 0, ial as u32));
                let mut csq_type = node.csq;
                if node.ty == HAP_SSS {
                    csq.ty.ty = csq_type;
                    ctx.csq_stage(&mut csq, rec)?;
                    continue;
                }
                let sref = &aux.sref;
                let alen = node.seq.len();
                let fill = node.dlen % 3 != 0 && alen > 0;
                let sbeg = node.sbeg as usize;
                let rlen = node.rlen.max(0) as usize;
                let mut tseq = cds_translate(
                    sref,
                    &node.seq,
                    (nsref - 2 * pad as i64 + node.dlen as i64).max(0) as usize,
                    sbeg,
                    sbeg,
                    sbeg + rlen,
                    strand,
                    fill,
                )?;
                let rseq = sub(sref, (pad + sbeg) as i64, rlen as i64);
                let mut tref = cds_translate(
                    sref,
                    &rseq,
                    (nsref - 2 * pad as i64) as usize,
                    sbeg,
                    sbeg,
                    sbeg + rlen,
                    strand,
                    fill,
                )?;
                truncate_at_stop(&mut tref);
                truncate_at_stop(&mut tseq);
                let tref_last = tref.last().copied().unwrap_or(0);
                let tseq_last = tseq.last().copied().unwrap_or(0);
                if csq_type & CSQ_STOP_LOST != 0 {
                    if tref_last == b'*' && tref_last == tseq_last {
                        csq_type &= !CSQ_STOP_LOST;
                        csq_type |= CSQ_STOP_RETAINED;
                    } else if tref_last != b'*' {
                        if tseq_last == b'*' {
                            csq_type &= !CSQ_STOP_GAINED;
                            csq_type |= CSQ_STOP_RETAINED;
                        } else {
                            csq_type |= CSQ_INCOMPLETE_CDS;
                        }
                    }
                }
                if csq_type & CSQ_START_LOST != 0 && tref.first().copied() != Some(b'M') {
                    csq_type &= !CSQ_START_LOST;
                }
                if node.dlen != 0 {
                    if node.dlen % 3 != 0 {
                        csq_type |= CSQ_FRAMESHIFT_VARIANT;
                    } else if node.dlen < 0 {
                        csq_type |= CSQ_INFRAME_DELETION;
                    } else {
                        csq_type |= CSQ_INFRAME_INSERTION;
                    }
                    if tref_last != b'*' && tseq_last == b'*' {
                        csq_type |= CSQ_STOP_GAINED;
                    }
                } else {
                    let mut aa_change = false;
                    for j in 0..tref.len() {
                        let t = tref[j];
                        let s = tseq.get(j).copied().unwrap_or(0);
                        if t == s {
                            continue;
                        }
                        aa_change = true;
                        if t == b'*' {
                            csq_type |= CSQ_STOP_LOST;
                        } else if s == b'*' {
                            csq_type |= CSQ_STOP_GAINED;
                        } else {
                            csq_type |= CSQ_MISSENSE_VARIANT;
                        }
                    }
                    if !aa_change {
                        csq_type |= CSQ_SYNONYMOUS_VARIANT;
                    }
                }
                let mut compound_ran = false;
                if csq_type & CSQ_COMPOUND != 0 {
                    let mut s = String::new();
                    let aa_rbeg = if strand == STRAND_FWD {
                        sbeg as i64 / 3 + 1
                    } else {
                        (nsref - 2 * pad as i64 - sbeg as i64 - rlen as i64) / 3 + 1
                    };
                    let aa_sbeg = if strand == STRAND_FWD {
                        sbeg as i64 / 3 + 1
                    } else {
                        (nsref - 2 * pad as i64 + node.dlen as i64 - sbeg as i64 - alen as i64) / 3
                            + 1
                    };
                    s.push('|');
                    s.push_str(&aa_rbeg.to_string());
                    kprint_aa_prediction(
                        self.opts.brief_predictions,
                        aa_rbeg.max(0) as usize,
                        &tref,
                        &mut s,
                    );
                    if csq_type & CSQ_SYNONYMOUS_VARIANT == 0 {
                        s.push('>');
                        s.push_str(&aa_sbeg.to_string());
                        kprint_aa_prediction(
                            self.opts.brief_predictions,
                            aa_sbeg.max(0) as usize,
                            &tseq,
                            &mut s,
                        );
                    }
                    s.push('|');
                    s.push_str(&(rec.pos + 1).to_string());
                    s.push_str(&node.var);
                    csq.ty.vstr = Some(s);
                    csq.ty.ty = csq_type & CSQ_COMPOUND;
                    ctx.csq_stage(&mut csq, rec)?;
                    compound_ran = true;
                }
                if csq_type & !CSQ_COMPOUND != 0 {
                    csq.ty.ty = csq_type & !CSQ_COMPOUND;
                    csq.ty.vstr = if compound_ran {
                        Some(String::new())
                    } else {
                        None
                    };
                    ctx.csq_stage(&mut csq, rec)?;
                }
            }
        }
        Ok(ret)
    }

    fn test_utr(&mut self, rec: &RecInfo) -> Result<bool> {
        let chr = drop_chr(self.opts.unify_chr_names, &rec.chrom).to_string();
        let hits = self
            .gff
            .idx_utr
            .overlap(&chr, rec.pos, rec.pos.wrapping_add(rec.rlen));
        if hits.is_empty() {
            return Ok(false);
        }
        let mut sp = Splice::init(rec, chr);
        let mut ret = false;
        let mut ctx = self.ctx();
        for ui in hits {
            let (tr_i, which, beg, end) = {
                let u = &ctx.gff.utrs[ui];
                (u.tr, u.which, u.beg, u.end)
            };
            sp.tr = tr_i;
            for ial in 1..rec.n_allele() {
                if rec.allele_is_symbolic(ial) {
                    continue;
                }
                sp.vcf.alt = rec.alleles[ial].as_bytes().to_vec();
                sp.csq = 0;
                let sret = splice_csq(&mut ctx, &mut sp, &[], rec, beg, end)?;
                if sret != SPLICE_INSIDE && sret != SPLICE_OVERLAP {
                    continue;
                }
                let ty = if which == Prime::P5 {
                    CSQ_UTR5
                } else {
                    CSQ_UTR3
                };
                let mut csq = Csq::new(rec.pos, Vcsq::new(ctx.gff, tr_i, ty, ial as u32));
                ctx.csq_stage(&mut csq, rec)?;
                ret = true;
            }
        }
        Ok(ret)
    }

    fn test_splice(&mut self, rec: &RecInfo) -> Result<bool> {
        let chr = drop_chr(self.opts.unify_chr_names, &rec.chrom).to_string();
        let hits = self
            .gff
            .idx_exon
            .overlap(&chr, rec.pos, rec.pos.wrapping_add(rec.rlen));
        if hits.is_empty() {
            return Ok(false);
        }
        let mut sp = Splice::init(rec, chr);
        sp.check_acceptor = true;
        sp.check_donor = true;
        let mut ret = false;
        let mut ctx = self.ctx();
        for ei in hits {
            let (tr_i, ex_beg, ex_end) = {
                let e = &ctx.gff.exons[ei];
                (e.tr, e.beg, e.end)
            };
            let tr = &ctx.gff.tscripts[tr_i];
            if tr.cds.is_empty() {
                continue;
            }
            sp.tr = tr_i;
            sp.check_region_beg = tr.beg != ex_beg;
            sp.check_region_end = tr.end != ex_end;
            for ial in 1..rec.n_allele() {
                // bcftools checks the first ALT here, not allele `ial`.
                if rec.allele_is_symbolic(1) {
                    continue;
                }
                sp.vcf.alt = rec.alleles[ial].as_bytes().to_vec();
                sp.csq = 0;
                splice_csq(&mut ctx, &mut sp, &[], rec, ex_beg, ex_end)?;
                if sp.csq != 0 {
                    ret = true;
                }
            }
        }
        Ok(ret)
    }

    fn test_tscript(&mut self, rec: &RecInfo) -> Result<bool> {
        let chr = drop_chr(self.opts.unify_chr_names, &rec.chrom).to_string();
        let hits = self
            .gff
            .idx_tscript
            .overlap(&chr, rec.pos, rec.pos.wrapping_add(rec.rlen));
        if hits.is_empty() {
            return Ok(false);
        }
        let mut sp = Splice::init(rec, chr);
        let mut ret = false;
        let mut ctx = self.ctx();
        for tr_i in hits {
            let (beg, end, coding) = {
                let t = &ctx.gff.tscripts[tr_i];
                (t.beg, t.end, gf_is_coding(t.ty))
            };
            sp.tr = tr_i;
            for ial in 1..rec.n_allele() {
                if rec.allele_is_symbolic(ial) {
                    continue;
                }
                sp.vcf.alt = rec.alleles[ial].as_bytes().to_vec();
                sp.csq = 0;
                let sret = splice_csq(&mut ctx, &mut sp, &[], rec, beg, end)?;
                if sret != SPLICE_INSIDE && sret != SPLICE_OVERLAP {
                    continue;
                }
                let ty = if coding { CSQ_INTRON } else { CSQ_NON_CODING };
                let mut csq = Csq::new(rec.pos, Vcsq::new(ctx.gff, tr_i, ty, 0));
                ctx.csq_stage(&mut csq, rec)?;
                ret = true;
            }
        }
        Ok(ret)
    }

    /// `test_symbolic_alt`: `<INS...>` alleles, reported as feature elongation.
    fn test_symbolic_alt(&mut self, rec: &RecInfo) -> Result<()> {
        if self.opts.verbosity > 0 && !self.warn.symbolic_warned {
            eprintln!("Warning: The support for symbolic ALT insertions is experimental.");
            self.warn.symbolic_warned = true;
        }
        let chr = drop_chr(self.opts.unify_chr_names, &rec.chrom).to_string();
        let beg = rec.pos + 1;
        let end = beg;
        let csq_class = CSQ_ELONGATION;
        let mut hit = false;
        let mut ctx = self.ctx();
        for ci in ctx.gff.idx_cds.overlap(&chr, beg, end) {
            let tr_i = ctx.gff.cds[ci].tr;
            let coding = gf_is_coding(ctx.gff.tscripts[tr_i].ty);
            let ty = (if coding {
                CSQ_CODING_SEQUENCE
            } else {
                CSQ_NON_CODING
            }) | csq_class;
            let mut csq = Csq::new(rec.pos, Vcsq::new(ctx.gff, tr_i, ty, 0));
            ctx.csq_stage(&mut csq, rec)?;
            hit = true;
        }
        for ui in ctx.gff.idx_utr.overlap(&chr, beg, end) {
            let (tr_i, which) = {
                let u = &ctx.gff.utrs[ui];
                (u.tr, u.which)
            };
            let ty = (if which == Prime::P5 {
                CSQ_UTR5
            } else {
                CSQ_UTR3
            }) | csq_class;
            let mut csq = Csq::new(rec.pos, Vcsq::new(ctx.gff, tr_i, ty, 0));
            ctx.csq_stage(&mut csq, rec)?;
            hit = true;
        }
        let exons = ctx.gff.idx_exon.overlap(&chr, beg, end);
        if !exons.is_empty() {
            let mut sp = Splice::init(rec, chr.clone());
            sp.check_acceptor = true;
            sp.check_donor = true;
            for ei in exons {
                let (tr_i, ex_beg, ex_end) = {
                    let e = &ctx.gff.exons[ei];
                    (e.tr, e.beg, e.end)
                };
                let tr = &ctx.gff.tscripts[tr_i];
                if tr.cds.is_empty() {
                    continue;
                }
                sp.tr = tr_i;
                sp.check_region_beg = tr.beg != ex_beg;
                sp.check_region_end = tr.end != ex_end;
                sp.vcf.alt = rec.alleles[1].as_bytes().to_vec();
                sp.csq = csq_class;
                splice_csq(&mut ctx, &mut sp, &[], rec, ex_beg, ex_end)?;
                if sp.csq != 0 {
                    hit = true;
                }
            }
        }
        if !hit {
            let trs = ctx.gff.idx_tscript.overlap(&chr, beg, end);
            let mut sp = Splice::init(rec, chr.clone());
            for tr_i in trs {
                let (tbeg, tend, coding) = {
                    let t = &ctx.gff.tscripts[tr_i];
                    (t.beg, t.end, gf_is_coding(t.ty))
                };
                sp.tr = tr_i;
                sp.vcf.alt = rec.alleles[1].as_bytes().to_vec();
                sp.csq = csq_class;
                let sret = splice_csq(&mut ctx, &mut sp, &[], rec, tbeg, tend)?;
                if sret != SPLICE_INSIDE && sret != SPLICE_OVERLAP {
                    continue;
                }
                let ty = (if coding { CSQ_INTRON } else { CSQ_NON_CODING }) | csq_class;
                let mut csq = Csq::new(rec.pos, Vcsq::new(ctx.gff, tr_i, ty, 0));
                ctx.csq_stage(&mut csq, rec)?;
            }
        }
        Ok(())
    }

    /// `hap_flush`: finalize every transcript ending at or before `pos`.
    fn hap_flush(&mut self, pos: u32) -> Result<()> {
        while let Some((end, tr_i)) = self.active_tr.top() {
            if end > pos {
                break;
            }
            self.active_tr.delete();
            let has_tree = self.aux[tr_i]
                .as_ref()
                .is_some_and(|a| a.root_has_children());
            if has_tree {
                let mut ctx = Ctx {
                    gff: &self.gff,
                    opts: &self.opts,
                    buf: &mut self.buf,
                    smpl: &self.smpl,
                    samples: &self.samples,
                    warn: &mut self.warn,
                    out: &mut self.out,
                    chrom: self.rid.as_deref().unwrap_or(""),
                };
                let aux = self.aux[tr_i].as_mut().unwrap();
                hap_finalize(&mut ctx, &mut self.hap, aux, tr_i)?;
                if ctx.out.is_text() {
                    if ctx.opts.phase == Phase::DropGt {
                        hap_print_text(&mut ctx, aux, None, 0, aux.hap[0])?;
                    } else {
                        for i in 0..ctx.smpl.len() {
                            let ismpl = ctx.smpl[i];
                            for j in 0..2 {
                                let h = aux.hap[i * 2 + j];
                                hap_print_text(&mut ctx, aux, Some(ismpl), j + 1, h)?;
                            }
                        }
                    }
                } else if ctx.opts.phase != Phase::DropGt {
                    for i in 0..ctx.smpl.len() {
                        let ismpl = ctx.smpl[i];
                        for j in 0..2 {
                            hap_stage_vcf(&mut ctx, aux, ismpl, j, aux.hap[i * 2 + j]);
                        }
                    }
                }
            }
            self.rm_tr.push(tr_i);
        }
        Ok(())
    }

    fn vbuf_push(&mut self, rec: RecInfo, line: VcfRecord) -> usize {
        let pos = rec.pos;
        let new = self.buf.vbufs.back().is_none_or(|v| v.pos != pos);
        if new {
            let slot = self.ring.append();
            self.buf.vbufs.push_back(Vbuf {
                pos,
                vrec: Vec::new(),
                keep_until: 0,
                slot,
            });
        }
        let id = self.buf.base + self.buf.vbufs.len() - 1;
        let fmt_bm = if self.opts.phase != Phase::DropGt && !self.smpl.is_empty() {
            vec![0u32; self.opts.hdr_nsmpl * self.opts.nfmt_bcsq]
        } else {
            Vec::new()
        };
        let vbuf = self.buf.vbufs.back_mut().unwrap();
        let nfmt = self.ring.stale_nfmt(vbuf.slot, vbuf.vrec.len());
        vbuf.vrec.push(Vrec {
            info: rec,
            line,
            fmt_bm,
            nfmt,
            vcsq: Vec::new(),
        });
        self.buf.pos2vbuf.insert(pos, id);
        id
    }

    /// `vbuf_flush`: write the position buffers no open transcript needs.
    fn vbuf_flush(&mut self, pos: u32) -> Result<()> {
        while !self.buf.vbufs.is_empty() {
            if !self.opts.local_csq
                && !self.active_tr.is_empty()
                && self.buf.vbufs[0].keep_until > pos
            {
                break;
            }
            let mut vbuf = self.buf.vbufs.pop_front().unwrap();
            self.buf.base += 1;
            self.ring.shift();
            for (k, vrec) in vbuf.vrec.iter_mut().enumerate() {
                self.ring.remember(vbuf.slot, k, vrec.nfmt);
                if self.out.is_text() {
                    continue;
                }
                let line = if vrec.vcsq.is_empty() {
                    format_line(&vrec.line)
                } else {
                    let mut s = String::new();
                    for (j, v) in vrec.vcsq.iter_mut().enumerate() {
                        if j > 0 {
                            s.push(',');
                        }
                        kput_vcsq(&self.gff, v, &mut s);
                    }
                    let mut line = vrec.line.clone();
                    line.info = update_info(&line.info, &self.opts.bcsq_tag, &s);
                    if self.opts.hdr_nsmpl > 0 {
                        let vals = (vrec.nfmt > 0).then(|| {
                            (0..self.opts.hdr_nsmpl)
                                .map(|i| {
                                    (0..vrec.nfmt)
                                        .map(|k| {
                                            vrec.fmt_bm[i * self.opts.nfmt_bcsq + k].to_string()
                                        })
                                        .collect::<Vec<_>>()
                                        .join(",")
                                })
                                .collect::<Vec<_>>()
                        });
                        update_format(&mut line, &self.opts.bcsq_tag, vals);
                    }
                    format_line(&line)
                };
                if let CsqOutput::Vcf(sink) = &mut self.out {
                    sink.write_line(&line)?;
                }
            }
            self.buf.pos2vbuf.remove(&vbuf.pos);
        }
        if !self.active_tr.is_empty() {
            return Ok(());
        }
        for tr_i in self.rm_tr.drain(..) {
            self.aux[tr_i] = None;
        }
        Ok(())
    }
}

/// Set or replace `tag=value` in an INFO column.
pub(crate) fn update_info(info: &str, tag: &str, value: &str) -> String {
    let mut parts: Vec<String> = if info == "." || info.is_empty() {
        Vec::new()
    } else {
        info.split(';').map(str::to_string).collect()
    };
    let new = format!("{tag}={value}");
    match parts
        .iter()
        .position(|p| p == tag || p.starts_with(tag) && p.as_bytes().get(tag.len()) == Some(&b'='))
    {
        Some(i) => parts[i] = new,
        None => parts.push(new),
    }
    parts.join(";")
}

/// Set (or with `None` remove) a per-sample FORMAT tag.
pub(crate) fn update_format(line: &mut VcfRecord, tag: &str, vals: Option<Vec<String>>) {
    let Some(fmt) = line.format.clone() else {
        return;
    };
    let keys: Vec<&str> = fmt.split(':').collect();
    let pos = keys.iter().position(|k| *k == tag);
    match (pos, vals) {
        (Some(k), Some(vals)) => {
            for (i, s) in line.samples.iter_mut().enumerate() {
                let mut f: Vec<String> = s.split(':').map(str::to_string).collect();
                while f.len() <= k {
                    f.push(".".to_string());
                }
                f[k] = vals.get(i).cloned().unwrap_or_else(|| "0".to_string());
                *s = f.join(":");
            }
        }
        (None, Some(vals)) => {
            line.format = Some(format!("{fmt}:{tag}"));
            for (i, s) in line.samples.iter_mut().enumerate() {
                let mut f: Vec<String> = s.split(':').map(str::to_string).collect();
                while f.len() < keys.len() {
                    f.push(".".to_string());
                }
                f.push(vals.get(i).cloned().unwrap_or_else(|| "0".to_string()));
                *s = f.join(":");
            }
        }
        (Some(k), None) => {
            let new_keys: Vec<&str> = keys
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != k)
                .map(|(_, v)| *v)
                .collect();
            line.format = Some(new_keys.join(":"));
            for s in line.samples.iter_mut() {
                let f: Vec<&str> = s
                    .split(':')
                    .enumerate()
                    .filter(|(i, _)| *i != k)
                    .map(|(_, v)| v)
                    .collect();
                *s = f.join(":");
            }
        }
        (None, None) => {}
    }
}

fn format_line(rec: &VcfRecord) -> String {
    let mut s = String::with_capacity(128);
    s.push_str(&rec.chrom);
    s.push('\t');
    s.push_str(&rec.pos.to_string());
    s.push('\t');
    s.push_str(&rec.id);
    s.push('\t');
    s.push_str(&rec.ref_allele);
    s.push('\t');
    s.push_str(&rec.alt);
    s.push('\t');
    s.push_str(&rec.qual);
    s.push('\t');
    s.push_str(&rec.filter);
    s.push('\t');
    s.push_str(&rec.info);
    if let Some(f) = &rec.format {
        s.push('\t');
        s.push_str(f);
        for smp in &rec.samples {
            s.push('\t');
            s.push_str(smp);
        }
    }
    s
}
