//! Consequences around exon boundaries: splice sites, start/stop overlap and
//! allele trimming (`splice_*` and `shifted_del_synonymous` in csq.c).
//! `ref_beg`/`ref_end` are `u32` like the C fields, so comparisons against
//! `ex_beg - N` wrap the same way near position 0.

use anyhow::{Result, bail};

use super::caller::Ctx;
use super::gff::{N_SPLICE_DONOR, N_SPLICE_REGION_INTRON, STRAND_FWD, STRAND_REV, Tscript};
use super::*;

pub(crate) const SPLICE_VAR_REF: i32 = 0;
pub(crate) const SPLICE_OUTSIDE: i32 = 1;
pub(crate) const SPLICE_INSIDE: i32 = 2;
pub(crate) const SPLICE_OVERLAP: i32 = 3;

pub(crate) struct SpliceVcf {
    pub pos: i32,
    pub rlen: i32,
    pub alen: i32,
    pub ial: i32,
    pub ref_: Vec<u8>,
    pub alt: Vec<u8>,
}

pub(crate) struct Splice {
    pub tr: usize,
    /// Record chromosome with the `chr` prefix unified.
    pub chr: String,
    pub vcf: SpliceVcf,
    pub check_acceptor: bool,
    pub check_start: bool,
    pub check_stop: bool,
    pub check_donor: bool,
    pub check_region_beg: bool,
    pub check_region_end: bool,
    pub check_utr: bool,
    pub set_refalt: bool,
    pub csq: u32,
    pub tbeg: i32,
    pub tend: i32,
    pub ref_beg: u32,
    pub ref_end: u32,
    pub kref: Vec<u8>,
    pub kalt: Vec<u8>,
}

impl Splice {
    pub fn init(rec: &RecInfo, chr: String) -> Self {
        Self {
            tr: 0,
            chr,
            vcf: SpliceVcf {
                pos: rec.pos as i32,
                rlen: rec.rlen as i32,
                alen: 0,
                ial: 0,
                ref_: rec.alleles[0].as_bytes().to_vec(),
                alt: Vec::new(),
            },
            check_acceptor: false,
            check_start: false,
            check_stop: false,
            check_donor: false,
            check_region_beg: false,
            check_region_end: false,
            check_utr: false,
            set_refalt: false,
            csq: 0,
            tbeg: 0,
            tend: 0,
            ref_beg: 0,
            ref_end: 0,
            kref: Vec::new(),
            kalt: Vec::new(),
        }
    }
}

/// `strncmp(a, b, n) == 0` with bytes past either end reading as NUL.
fn strncmp_eq(a: &[u8], b: &[u8], n: usize) -> bool {
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(0);
        let y = b.get(i).copied().unwrap_or(0);
        if x != y {
            return false;
        }
        if x == 0 {
            return true;
        }
    }
    true
}

/// `splice_build_hap`: the reference and alternate sequence of `len` bases
/// from `beg` (`len < 0`: `beg` is the last base, filled from the left).
pub(crate) fn splice_build_hap(sp: &mut Splice, tref: &[u8], tr: &Tscript, beg: u32, len: i32) {
    let pad = N_REF_PAD as i64;
    let tr_beg = tr.beg as i64;
    let tr_end = tr.end as i64;
    let (mut rlen, mut alen, rbeg, abeg): (i64, i64, i64, i64);
    if len < 0 {
        rlen = -len as i64;
        alen = rlen;
        rbeg = beg as i64 - rlen + 1;
        let mut dlen = (sp.vcf.alen - sp.vcf.rlen) as i64;
        if dlen < 0 && beg < sp.ref_end {
            dlen += (sp.ref_end - beg) as i64;
        }
        abeg = rbeg + dlen;
    } else {
        rbeg = beg as i64;
        abeg = rbeg;
        rlen = len as i64;
        alen = rlen;
    }
    sp.kref.clear();
    sp.kalt.clear();
    let vpos = sp.vcf.pos as i64;
    let vrlen = sp.vcf.rlen as i64;
    let valen = sp.vcf.alen as i64;

    let roff: i64 = if rbeg < vpos {
        sp.kref.extend(sub(tref, pad + rbeg - tr_beg, vpos - rbeg));
        0
    } else {
        rbeg - vpos
    };
    if roff < vrlen && (sp.kref.len() as i64) < rlen {
        let mut l = vrlen - roff;
        if l > rlen - sp.kref.len() as i64 {
            l = rlen - sp.kref.len() as i64;
        }
        sp.kref.extend(sub(&sp.vcf.ref_, roff, l));
    }
    let end = vpos + vrlen;
    if (sp.kref.len() as i64) < rlen {
        if end + rlen - sp.kref.len() as i64 - 1 > tr_end {
            rlen -= end + rlen - sp.kref.len() as i64 - 1 - tr_end;
        }
        if (sp.kref.len() as i64) < rlen {
            let l = rlen - sp.kref.len() as i64;
            sp.kref.extend(sub(tref, pad + end - tr_beg, l));
        }
    }
    let mut aoff: i64;
    if abeg < vpos {
        sp.kalt.extend(sub(tref, pad + abeg - tr_beg, vpos - abeg));
        aoff = 0;
    } else {
        aoff = abeg - vpos;
    }
    if aoff < valen && (sp.kalt.len() as i64) < alen {
        let mut l = valen - aoff;
        if l > alen - sp.kalt.len() as i64 {
            l = alen - sp.kalt.len() as i64;
        }
        sp.kalt.extend(sub(&sp.vcf.alt, aoff, l));
        aoff -= l;
    }
    if aoff < 0 {
        aoff = 0;
    } else {
        aoff -= 1;
    }
    let end = vpos + vrlen;
    if (sp.kalt.len() as i64) < alen {
        if end + alen + aoff - sp.kalt.len() as i64 - 1 > tr_end {
            alen -= end + alen + aoff - sp.kalt.len() as i64 - 1 - tr_end;
        }
        if alen > 0 && alen > sp.kalt.len() as i64 {
            let l = alen - sp.kalt.len() as i64;
            sp.kalt.extend(sub(tref, aoff + pad + end - tr_beg, l));
        }
    }
}

impl Ctx<'_> {
    /// `csq_stage_utr`: stage the first UTR of transcript `tr` among `hits`.
    pub(crate) fn csq_stage_utr(
        &mut self,
        hits: &[usize],
        rec: &RecInfo,
        tr: usize,
        ty: u32,
        ial: i32,
    ) -> Result<u32> {
        for &ui in hits {
            let utr = &self.gff.utrs[ui];
            if utr.tr != tr {
                continue;
            }
            let t = (if utr.which == gff::Prime::P5 {
                CSQ_UTR5
            } else {
                CSQ_UTR3
            }) | ty;
            let mut csq = Csq::new(rec.pos, Vcsq::new(self.gff, tr, t, ial as u32));
            self.csq_stage(&mut csq, rec)?;
            return Ok(csq.ty.ty);
        }
        Ok(0)
    }

    pub(crate) fn csq_stage_splice(
        &mut self,
        rec: &RecInfo,
        tr: usize,
        ty: u32,
        ial: i32,
    ) -> Result<()> {
        if ty == 0 {
            return Ok(());
        }
        let mut csq = Csq::new(rec.pos, Vcsq::new(self.gff, tr, ty, ial as u32));
        self.csq_stage(&mut csq, rec)
    }
}

fn region_and_donor_after(
    sp: &mut Splice,
    strand: u8,
    ex_end: u32,
    refalt: Option<(&[u8], &[u8])>,
    strict_region: bool,
    syn_off: usize,
) {
    let n8 = N_SPLICE_REGION_INTRON;
    let n2 = N_SPLICE_DONOR;
    let region_beg = if strict_region {
        sp.ref_beg < ex_end.wrapping_add(n8)
    } else {
        sp.ref_beg <= ex_end.wrapping_add(n8)
    };
    if region_beg && sp.ref_end > ex_end.wrapping_add(n2) {
        sp.csq |= CSQ_SPLICE_REGION;
        if let Some((r, a)) = refalt
            && strncmp_eq(r, a, n8 as usize)
        {
            sp.csq |= CSQ_SYNONYMOUS_VARIANT;
        }
    }
    let donor = if strict_region {
        sp.ref_beg < ex_end.wrapping_add(n2)
    } else {
        sp.ref_beg <= ex_end.wrapping_add(n2)
    };
    if donor {
        if sp.check_donor && strand == STRAND_FWD {
            sp.csq |= CSQ_SPLICE_DONOR;
        }
        if sp.check_acceptor && strand == STRAND_REV {
            sp.csq |= CSQ_SPLICE_ACCEPTOR;
        }
        if let Some((r, a)) = refalt {
            let r = r.get(syn_off..).unwrap_or(&[]);
            let a = a.get(syn_off..).unwrap_or(&[]);
            if strncmp_eq(r, a, n2 as usize) {
                sp.csq |= CSQ_SYNONYMOUS_VARIANT;
            }
        }
    }
}

fn region_and_donor_before(
    sp: &mut Splice,
    strand: u8,
    ex_beg: u32,
    refalt: Option<(&[u8], &[u8])>,
    inclusive: bool,
) {
    let n8 = N_SPLICE_REGION_INTRON;
    let n2 = N_SPLICE_DONOR;
    let region_end = if inclusive {
        sp.ref_end >= ex_beg.wrapping_sub(n8)
    } else {
        sp.ref_end > ex_beg.wrapping_sub(n8)
    };
    if region_end && sp.ref_beg < ex_beg.wrapping_sub(n2) {
        sp.csq |= CSQ_SPLICE_REGION;
        if let Some((r, a)) = refalt
            && strncmp_eq(r, a, n8 as usize)
        {
            sp.csq |= CSQ_SYNONYMOUS_VARIANT;
        }
    }
    let donor = if inclusive {
        sp.ref_end >= ex_beg.wrapping_sub(n2)
    } else {
        sp.ref_end > ex_beg.wrapping_sub(n2)
    };
    if donor {
        if sp.check_donor && strand == STRAND_REV {
            sp.csq |= CSQ_SPLICE_DONOR;
        }
        if sp.check_acceptor && strand == STRAND_FWD {
            sp.csq |= CSQ_SPLICE_ACCEPTOR;
        }
        if let Some((r, a)) = refalt {
            let off = (n8 - n2) as usize;
            let r = r.get(off..).unwrap_or(&[]);
            let a = a.get(off..).unwrap_or(&[]);
            if strncmp_eq(r, a, n2 as usize) {
                sp.csq |= CSQ_SYNONYMOUS_VARIANT;
            }
        }
    }
}

fn splice_csq_ins(
    ctx: &mut Ctx,
    sp: &mut Splice,
    tref: &[u8],
    rec: &RecInfo,
    ex_beg: u32,
    ex_end: u32,
) -> Result<i32> {
    let tr = &ctx.gff.tscripts[sp.tr];
    let strand = tr.strand;
    if sp.tbeg != 0 || sp.vcf.ref_.first() != sp.vcf.alt.first() {
        sp.ref_beg = (sp.vcf.pos + sp.tbeg - 1) as u32;
        sp.ref_end = (sp.vcf.pos + sp.vcf.rlen - sp.tend) as u32;
    } else {
        if sp.tend != 0 {
            sp.tend -= 1;
        }
        sp.ref_beg = sp.vcf.pos as u32;
        sp.ref_end = (sp.vcf.pos + sp.vcf.rlen - sp.tend) as u32;
    }
    if sp.ref_beg >= ex_end {
        if sp.check_utr {
            let p = sp.ref_beg.wrapping_add(1);
            let hits = ctx.gff.idx_utr.overlap(&sp.chr, p, p);
            if !hits.is_empty() {
                let ret = ctx.csq_stage_utr(&hits, rec, sp.tr, sp.csq, sp.vcf.ial)?;
                if ret != 0 {
                    return Ok(SPLICE_OUTSIDE);
                }
            }
        }
        if !sp.check_region_end {
            return Ok(SPLICE_OUTSIDE);
        }
        let mut refalt = None;
        let (kr, ka);
        if sp.set_refalt {
            splice_build_hap(
                sp,
                tref,
                tr,
                ex_end.wrapping_add(1),
                N_SPLICE_REGION_INTRON as i32,
            );
            kr = sp.kref.clone();
            ka = sp.kalt.clone();
            refalt = Some((kr.as_slice(), ka.as_slice()));
        }
        region_and_donor_after(sp, strand, ex_end, refalt, true, 0);
        ctx.csq_stage_splice(rec, sp.tr, sp.csq, sp.vcf.ial)?;
        return Ok(SPLICE_OUTSIDE);
    }
    if sp.ref_end < ex_beg || (sp.ref_end == ex_beg && !sp.check_region_beg) {
        if sp.check_utr {
            let p = sp.ref_end.wrapping_sub(1);
            let hits = ctx.gff.idx_utr.overlap(&sp.chr, p, p);
            if !hits.is_empty() {
                let ret = ctx.csq_stage_utr(&hits, rec, sp.tr, sp.csq, sp.vcf.ial)?;
                if ret != 0 {
                    return Ok(SPLICE_OUTSIDE);
                }
            }
        }
        if !sp.check_region_beg {
            return Ok(SPLICE_OUTSIDE);
        }
        let mut refalt = None;
        let (kr, ka);
        if sp.set_refalt {
            splice_build_hap(
                sp,
                tref,
                tr,
                ex_beg.wrapping_sub(N_SPLICE_REGION_INTRON),
                N_SPLICE_REGION_INTRON as i32,
            );
            kr = sp.kref.clone();
            ka = sp.kalt.clone();
            refalt = Some((kr.as_slice(), ka.as_slice()));
        }
        region_and_donor_before(sp, strand, ex_beg, refalt, false);
        ctx.csq_stage_splice(rec, sp.tr, sp.csq, sp.vcf.ial)?;
        return Ok(SPLICE_OUTSIDE);
    }
    if sp.ref_beg <= ex_beg.wrapping_add(2) {
        if sp.check_region_beg {
            sp.csq |= CSQ_SPLICE_REGION;
        }
        if strand == STRAND_FWD {
            if sp.check_start {
                sp.csq |= CSQ_START_LOST;
            }
        } else if sp.check_stop {
            sp.csq |= CSQ_STOP_LOST;
        }
    }
    if sp.ref_end > ex_end.wrapping_sub(2) {
        if sp.check_region_end {
            sp.csq |= CSQ_SPLICE_REGION;
        }
        if strand == STRAND_REV {
            if sp.check_start {
                sp.csq |= CSQ_START_LOST;
            }
        } else if sp.check_stop {
            sp.csq |= CSQ_STOP_LOST;
        }
    }
    if sp.set_refalt {
        if sp.ref_beg < sp.vcf.pos as u32 {
            let dlen = sp.vcf.pos - sp.ref_beg as i32;
            if dlen != 1 {
                bail!(
                    "csq: unexpected insertion trimming at {}:{}",
                    sp.chr,
                    sp.vcf.pos + 1
                );
            }
            sp.tbeg += dlen;
            if sp.tbeg + sp.tend == sp.vcf.rlen {
                sp.tend -= dlen;
            }
            sp.ref_beg = sp.vcf.pos as u32;
        }
        if sp.ref_end == ex_beg {
            sp.tend -= 1;
        }
        let l = sp.vcf.alen - sp.tend - sp.tbeg + 1;
        splice_build_hap(sp, tref, tr, sp.ref_beg, l);
        sp.vcf.rlen -= sp.tbeg + sp.tend - 1;
        if sp.vcf.rlen >= 0 && sp.kref.len() > sp.vcf.rlen as usize {
            sp.kref.truncate(sp.vcf.rlen as usize);
        }
    }
    ctx.csq_stage_splice(rec, sp.tr, sp.csq, sp.vcf.ial)?;
    Ok(SPLICE_INSIDE)
}

/// `shifted_del_synonymous`: can the deleted bases be replaced by the
/// reference sequence that follows, leaving the start codon intact?
fn shifted_del_synonymous(
    ctx: &mut Ctx,
    sp: &Splice,
    tref: &[u8],
    ex_beg: u32,
    ex_end: u32,
) -> bool {
    let tr = &ctx.gff.tscripts[sp.tr];
    if tr.strand == STRAND_REV && ((sp.vcf.pos + sp.vcf.rlen + 2) as u32) <= ex_end {
        return false;
    }
    if tr.strand == STRAND_FWD && sp.vcf.pos as u32 >= ex_beg.wrapping_add(3) {
        return false;
    }
    let ref_len = sp.vcf.ref_.len() as i64;
    let alt_len = sp.vcf.alt.len() as i64;
    if ref_len <= alt_len {
        return false;
    }
    let ndel = ref_len - alt_len;
    let pad = N_REF_PAD as i64;
    let pos = sp.vcf.pos as i64;
    let ptr_vcf = &sp.vcf.ref_[alt_len as usize..];

    let ref_off: i64 = if tr.strand == STRAND_REV {
        let vcf_ref_end = pos + ref_len - 1;
        let tr_ref_end = tr.end as i64 + pad;
        if vcf_ref_end + ndel > tr_ref_end {
            if !ctx.warn.small_ref_padding_warned {
                eprintln!(
                    "Warning: Could not verify synonymous start/stop at {}:{} due to small N_REF_PAD. (Improve me?)",
                    ctx.chrom,
                    sp.vcf.pos + 1
                );
                ctx.warn.small_ref_padding_warned = true;
            }
            return false;
        }
        pad + (vcf_ref_end + 1 - tr.beg as i64)
    } else {
        let vcf_block_beg = pos + ref_len - 2 * ndel;
        if vcf_block_beg < 0 {
            return false;
        }
        if ((pad + vcf_block_beg) as u32) < ex_beg {
            if !ctx.warn.small_ref_padding_warned {
                eprintln!(
                    "Warning: Could not verify synonymous start/stop at {}:{} due to small N_REF_PAD. (Improve me?)",
                    ctx.chrom,
                    sp.vcf.pos + 1
                );
                ctx.warn.small_ref_padding_warned = true;
            }
            return false;
        }
        pad + vcf_block_beg - tr.beg as i64
    };
    let mut i = 0i64;
    while c_at(ptr_vcf, i) != 0 && c_at(ptr_vcf, i) == c_at(tref, ref_off + i) {
        i += 1;
    }
    c_at(ptr_vcf, i) == 0
}

fn splice_csq_del(
    ctx: &mut Ctx,
    sp: &mut Splice,
    tref: &[u8],
    rec: &RecInfo,
    ex_beg: u32,
    ex_end: u32,
) -> Result<i32> {
    if sp.check_start && shifted_del_synonymous(ctx, sp, tref, ex_beg, ex_end) {
        sp.csq |= CSQ_START_RETAINED;
        return Ok(SPLICE_OVERLAP);
    }
    let strand = ctx.gff.tscripts[sp.tr].strand;
    sp.ref_beg = (sp.vcf.pos + sp.tbeg - 1) as u32;
    sp.ref_end = (sp.vcf.pos + sp.vcf.rlen - sp.tend - 1) as u32;
    if sp.ref_beg.wrapping_add(1) < ex_beg {
        if sp.check_region_beg {
            let mut csq = 0;
            if sp.check_utr {
                let hits = ctx
                    .gff
                    .idx_utr
                    .overlap(&sp.chr, sp.ref_beg, ex_beg.wrapping_sub(1));
                if !hits.is_empty() {
                    csq = ctx.csq_stage_utr(&hits, rec, sp.tr, sp.csq, sp.vcf.ial)?;
                }
            }
            if csq == 0 {
                let mut refalt = None;
                let (kr, ka);
                if sp.set_refalt {
                    let tr = &ctx.gff.tscripts[sp.tr];
                    splice_build_hap(
                        sp,
                        tref,
                        tr,
                        ex_beg.wrapping_sub(N_SPLICE_REGION_INTRON),
                        N_SPLICE_REGION_INTRON as i32,
                    );
                    kr = sp.kref.clone();
                    ka = sp.kalt.clone();
                    refalt = Some((kr.as_slice(), ka.as_slice()));
                }
                region_and_donor_before(sp, strand, ex_beg, refalt, true);
            }
        }
        if sp.ref_end >= ex_beg {
            sp.tbeg = sp.ref_beg.wrapping_sub(sp.vcf.pos as u32).wrapping_add(1) as i32;
            sp.ref_beg = ex_beg.wrapping_sub(1);
            if sp.tbeg + sp.tend == sp.vcf.alen {
                if sp.tend == 0 {
                    sp.csq |= CSQ_CODING_SEQUENCE;
                    return Ok(SPLICE_OVERLAP);
                }
                sp.tend -= 1;
            }
        }
    }
    if ex_end < sp.ref_end {
        if sp.check_region_end {
            let mut csq = 0;
            if sp.check_utr {
                let hits = ctx
                    .gff
                    .idx_utr
                    .overlap(&sp.chr, ex_end.wrapping_add(1), sp.ref_end);
                if !hits.is_empty() {
                    csq = ctx.csq_stage_utr(&hits, rec, sp.tr, sp.csq, sp.vcf.ial)?;
                }
            }
            if csq == 0 {
                let mut refalt = None;
                let (kr, ka);
                if sp.set_refalt {
                    let tr = &ctx.gff.tscripts[sp.tr];
                    splice_build_hap(
                        sp,
                        tref,
                        tr,
                        ex_end.wrapping_add(1),
                        N_SPLICE_REGION_INTRON as i32,
                    );
                    kr = sp.kref.clone();
                    ka = sp.kalt.clone();
                    refalt = Some((kr.as_slice(), ka.as_slice()));
                }
                region_and_donor_after(
                    sp,
                    strand,
                    ex_end,
                    refalt,
                    true,
                    N_SPLICE_REGION_INTRON as usize - N_SPLICE_DONOR as usize,
                );
            }
        }
        if sp.ref_beg < ex_end {
            sp.tend = (sp.vcf.rlen as u32)
                .wrapping_sub(sp.ref_end.wrapping_sub(sp.vcf.pos as u32).wrapping_add(1))
                as i32;
            sp.ref_end = ex_end;
        }
    }
    if sp.ref_end < ex_beg || sp.ref_beg >= ex_end {
        ctx.csq_stage_splice(rec, sp.tr, sp.csq, sp.vcf.ial)?;
        return Ok(SPLICE_OUTSIDE);
    }
    if sp.ref_beg < ex_beg.wrapping_add(2) {
        if sp.check_region_beg {
            sp.csq |= CSQ_SPLICE_REGION;
        }
        if strand == STRAND_FWD {
            if sp.check_start {
                sp.csq |= CSQ_START_LOST;
            }
        } else if sp.check_stop {
            sp.csq |= CSQ_STOP_LOST;
        }
    }
    if sp.ref_end > ex_end.wrapping_sub(3) {
        if sp.check_region_end {
            sp.csq |= CSQ_SPLICE_REGION;
        }
        if strand == STRAND_REV {
            if sp.check_start {
                sp.csq |= CSQ_START_LOST;
            }
        } else if sp.check_stop {
            sp.csq |= CSQ_STOP_LOST;
        }
    }
    if sp.set_refalt {
        if sp.tbeg > 0 {
            sp.tbeg -= 1;
        }
        if sp.vcf.rlen > sp.tbeg + sp.tend && sp.vcf.alen > sp.tbeg + sp.tend {
            sp.vcf.rlen -= sp.tbeg + sp.tend;
            sp.vcf.alen -= sp.tbeg + sp.tend;
        }
        sp.kref = sub(&sp.vcf.ref_, sp.tbeg as i64, sp.vcf.rlen as i64);
        sp.kalt = sub(&sp.vcf.alt, sp.tbeg as i64, sp.vcf.alen as i64);
        if (sp.ref_beg.wrapping_add(1) < ex_beg && sp.ref_end >= ex_beg)
            || (sp.ref_beg.wrapping_add(1) < ex_end && sp.ref_end >= ex_end)
        {
            sp.csq |= if !sp.ref_end.wrapping_sub(sp.ref_beg).is_multiple_of(3) {
                CSQ_FRAMESHIFT_VARIANT
            } else {
                CSQ_INFRAME_DELETION
            };
            return Ok(SPLICE_OVERLAP);
        }
    }
    ctx.csq_stage_splice(rec, sp.tr, sp.csq, sp.vcf.ial)?;
    Ok(SPLICE_INSIDE)
}

fn splice_csq_mnp(
    ctx: &mut Ctx,
    sp: &mut Splice,
    rec: &RecInfo,
    ex_beg: u32,
    ex_end: u32,
) -> Result<i32> {
    if sp.tbeg + sp.tend == sp.vcf.rlen {
        return Ok(SPLICE_VAR_REF);
    }
    let strand = ctx.gff.tscripts[sp.tr].strand;
    sp.ref_beg = (sp.vcf.pos + sp.tbeg) as u32;
    sp.ref_end = (sp.vcf.pos + sp.vcf.rlen - sp.tend - 1) as u32;
    if sp.ref_beg < ex_beg {
        if sp.check_region_beg {
            let mut csq = 0;
            if sp.check_utr {
                let hits = ctx
                    .gff
                    .idx_utr
                    .overlap(&sp.chr, sp.ref_beg, ex_beg.wrapping_sub(1));
                if !hits.is_empty() {
                    csq = ctx.csq_stage_utr(&hits, rec, sp.tr, sp.csq, sp.vcf.ial)?;
                }
            }
            if csq == 0 {
                region_and_donor_before(sp, strand, ex_beg, None, true);
            }
        }
        if sp.ref_end >= ex_beg {
            sp.tbeg = sp.ref_beg.wrapping_sub(sp.vcf.pos as u32) as i32;
            sp.ref_beg = ex_beg;
        }
    }
    if ex_end < sp.ref_end {
        if sp.check_region_end {
            let mut csq = 0;
            if sp.check_utr {
                let hits = ctx
                    .gff
                    .idx_utr
                    .overlap(&sp.chr, ex_end.wrapping_add(1), sp.ref_end);
                if !hits.is_empty() {
                    csq = ctx.csq_stage_utr(&hits, rec, sp.tr, sp.csq, sp.vcf.ial)?;
                }
            }
            if csq == 0 {
                region_and_donor_after(sp, strand, ex_end, None, false, 0);
            }
        }
        if sp.ref_beg <= ex_end {
            sp.tend = (sp.vcf.rlen as u32)
                .wrapping_sub(sp.ref_end.wrapping_sub(sp.vcf.pos as u32).wrapping_add(1))
                as i32;
            sp.ref_end = ex_end;
        }
    }
    if sp.ref_end < ex_beg || sp.ref_beg > ex_end {
        ctx.csq_stage_splice(rec, sp.tr, sp.csq, sp.vcf.ial)?;
        return Ok(SPLICE_OUTSIDE);
    }
    if sp.ref_beg < ex_beg.wrapping_add(3) {
        if sp.check_region_beg {
            sp.csq |= CSQ_SPLICE_REGION;
        }
        if strand == STRAND_FWD {
            if sp.check_start {
                sp.csq |= CSQ_START_LOST;
            }
        } else if sp.check_stop {
            sp.csq |= CSQ_STOP_LOST;
        }
    }
    if sp.ref_end > ex_end.wrapping_sub(3) {
        if sp.check_region_end {
            sp.csq |= CSQ_SPLICE_REGION;
        }
        if strand == STRAND_REV {
            if sp.check_start {
                sp.csq |= CSQ_START_LOST;
            }
        } else if sp.check_stop {
            sp.csq |= CSQ_STOP_LOST;
        }
    }
    if sp.set_refalt {
        sp.vcf.rlen -= sp.tbeg + sp.tend;
        sp.kref = sub(&sp.vcf.ref_, sp.tbeg as i64, sp.vcf.rlen as i64);
        sp.kalt = sub(&sp.vcf.alt, sp.tbeg as i64, sp.vcf.rlen as i64);
    }
    ctx.csq_stage_splice(rec, sp.tr, sp.csq, sp.vcf.ial)?;
    Ok(SPLICE_INSIDE)
}

/// `splice_csq`: trim shared bases from both ends of REF/ALT and classify
/// the variant against the exon `[ex_beg, ex_end]`.
pub(crate) fn splice_csq(
    ctx: &mut Ctx,
    sp: &mut Splice,
    tref: &[u8],
    rec: &RecInfo,
    ex_beg: u32,
    ex_end: u32,
) -> Result<i32> {
    sp.vcf.alen = sp.vcf.alt.len() as i32;
    let mut rlen1 = sp.vcf.rlen - 1;
    let mut alen1 = sp.vcf.alen - 1;
    let mut i = 0i32;
    sp.tbeg = 0;
    sp.tend = 0;
    while i <= rlen1 && i <= alen1 {
        if c_at(&sp.vcf.ref_, (rlen1 - i) as i64) != c_at(&sp.vcf.alt, (alen1 - i) as i64) {
            break;
        }
        i += 1;
    }
    sp.tend = i;
    rlen1 -= i;
    alen1 -= i;
    i = 0;
    while i <= rlen1 && i <= alen1 {
        if c_at(&sp.vcf.ref_, i as i64) != c_at(&sp.vcf.alt, i as i64) {
            break;
        }
        i += 1;
    }
    sp.tbeg = i;
    if sp.vcf.rlen == sp.vcf.alen {
        return splice_csq_mnp(ctx, sp, rec, ex_beg, ex_end);
    }
    if sp.vcf.rlen < sp.vcf.alen {
        return splice_csq_ins(ctx, sp, tref, rec, ex_beg, ex_end);
    }
    splice_csq_del(ctx, sp, tref, rec, ex_beg, ex_end)
}
