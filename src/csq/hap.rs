//! Per-transcript haplotype tree and the assembly of coding consequences
//! (`hap_init`, `hap_finalize`, `hap_add_csq`, `tscript_splice_ref` in csq.c).

use anyhow::{Result, bail};

use super::caller::Ctx;
use super::gff::{STRAND_FWD, TRIM_3PRIME, TRIM_5PRIME, Tscript};
use super::splice::{SPLICE_OUTSIDE, SPLICE_OVERLAP, SPLICE_VAR_REF, Splice, splice_csq};
use super::*;

/// Node of a haplotype tree: one applied variant and the CDS segment
/// between the previous node and this one.
#[derive(Default, Debug)]
pub(crate) struct HapNode {
    /// CDS segment `[parent, this)` on the reference strand, ALT included.
    pub seq: Vec<u8>,
    /// `ref>alt`.
    pub var: String,
    pub ty: u8,
    pub csq: u32,
    /// ALT minus REF length.
    pub dlen: i32,
    /// Variant start (0-based) after trimming.
    pub rbeg: u32,
    pub rlen: i32,
    /// Position on the spliced reference transcript.
    pub sbeg: u32,
    pub icds: u32,
    pub children: Vec<usize>,
    pub cur_rec: Option<RecId>,
    /// Record id and position.
    pub rec: Option<(RecId, u32)>,
    pub vcf_ial: i32,
    /// Haplotypes ending here.
    pub nend: u32,
    /// Allele index to the active child for `cur_rec`.
    pub cur_child: Vec<i32>,
    pub csq_list: Vec<Csq>,
}

/// Runtime data of an open transcript (`tscript_t`).
pub(crate) struct TscriptAux {
    /// Reference padded with `N_REF_PAD` on both ends.
    pub ref_: Vec<u8>,
    /// Spliced reference, padded; empty until first needed.
    pub sref: Vec<u8>,
    /// Node arena; index 0 is the root in haplotype mode.
    pub nodes: Vec<HapNode>,
    /// Leaf of each sample haplotype.
    pub hap: Vec<Option<usize>>,
}

impl TscriptAux {
    pub fn root_has_children(&self) -> bool {
        self.nodes.first().is_some_and(|r| !r.children.is_empty())
    }
}

/// `tscript_splice_ref`: concatenate the CDS with `N_REF_PAD` flanks. The
/// right flank copies the first bases of the last CDS, as in bcftools.
pub(crate) fn tscript_splice_ref(aux: &mut TscriptAux, tr: &Tscript, cds: &[gff::Cds]) {
    let pad = N_REF_PAD as i64;
    let mut s = Vec::new();
    let first = &cds[tr.cds[0]];
    s.extend(sub(&aux.ref_, first.beg as i64 - tr.beg as i64, pad));
    for &ci in &tr.cds {
        let c = &cds[ci];
        s.extend(sub(
            &aux.ref_,
            pad + c.beg as i64 - tr.beg as i64,
            c.len as i64,
        ));
    }
    let last = &cds[*tr.cds.last().unwrap()];
    s.extend(sub(&aux.ref_, pad + last.beg as i64 - tr.beg as i64, pad));
    aux.sref = s;
}

/// `hap_init`: apply allele `ial` of `rec` to `parent`, filling `child`.
/// Returns 0 when added, 1 for an overlapping variant, 2 when silently
/// discarded (intronic or REF-identical).
pub(crate) fn hap_init(
    ctx: &mut Ctx,
    aux: &TscriptAux,
    parent: &HapNode,
    child: &mut HapNode,
    ci: usize,
    rec: &RecInfo,
    ial: usize,
) -> Result<i32> {
    let gff = ctx.gff;
    let cds = &gff.cds[ci];
    let tr_i = cds.tr;
    let tr = &gff.tscripts[tr_i];
    let pad = N_REF_PAD as i64;
    child.icds = cds.icds;
    child.vcf_ial = ial as i32;
    let chr = ctx.drop_chr_prefix(&rec.chrom);
    let mut sp = Splice::init(rec, chr);
    sp.tr = tr_i;
    sp.vcf.ial = ial as i32;
    sp.vcf.alt = rec.alleles[ial].as_bytes().to_vec();
    sp.check_acceptor = true;
    sp.check_donor = true;
    sp.set_refalt = true;
    sp.check_utr = true;
    let ncds = tr.cds.len() as u32;
    if tr.trim & TRIM_5PRIME == 0 {
        if tr.strand == STRAND_FWD {
            if child.icds == 0 {
                sp.check_start = true;
            }
        } else if child.icds == ncds - 1 {
            sp.check_start = true;
        }
    }
    if tr.trim & TRIM_3PRIME == 0 {
        if tr.strand == STRAND_FWD {
            if child.icds == ncds - 1 {
                sp.check_stop = true;
            }
        } else if child.icds == 0 {
            sp.check_stop = true;
        }
    }
    if sp.check_start {
        let off = pad + cds.beg as i64 - tr.beg as i64;
        if tr.strand == STRAND_FWD {
            if dna2aa(&sub(&aux.ref_, off, 3)) != b'M' {
                sp.check_start = false;
            }
        } else if cdna2aa(&sub(&aux.ref_, off + cds.len as i64 - 3, 3)) != b'M' {
            sp.check_start = false;
        }
    }
    if child.icds != 0 {
        sp.check_region_beg = true;
    }
    if child.icds != ncds - 1 {
        sp.check_region_end = true;
    }
    let ret = splice_csq(
        ctx,
        &mut sp,
        &aux.ref_,
        rec,
        cds.beg,
        cds.beg.wrapping_add(cds.len).wrapping_sub(1),
    )?;
    if ret == SPLICE_VAR_REF {
        return Ok(2);
    }
    if ret == SPLICE_OUTSIDE || ret == SPLICE_OVERLAP || sp.csq == CSQ_START_LOST {
        if sp.csq == 0 {
            return Ok(2);
        }
        child.seq = Vec::new();
        child.sbeg = 0;
        child.rbeg = rec.pos;
        child.rlen = 0;
        child.dlen = 0;
        child.var = format!("{}>{}", rec.alleles[0], rec.alleles[ial]);
        child.ty = HAP_SSS;
        child.csq = sp.csq;
        child.rec = Some((rec.id, rec.pos));
        return Ok(0);
    }
    if sp.csq & CSQ_SYNONYMOUS_VARIANT != 0 {
        sp.csq &= !CSQ_SYNONYMOUS_VARIANT;
    }
    let mut dbeg = 0usize;
    if sp.ref_beg < cds.beg {
        dbeg = (cds.beg - sp.ref_beg) as usize;
        let l = sp.kref.len().saturating_sub(dbeg);
        sp.kref.truncate(l);
        sp.ref_beg = cds.beg;
        if dbeg > sp.kalt.len() {
            bail!(
                "csq: variant at {}:{} extends past the exon in an unexpected way",
                rec.chrom,
                rec.pos + 1
            );
        }
    }
    if parent.ty == HAP_SSS {
        bail!("csq: internal error, splice node used as a haplotype parent");
    }
    let mut s: Vec<u8> = Vec::new();
    let tr_beg = tr.beg as i64;
    if parent.ty == HAP_CDS {
        let mut i = parent.icds as usize;
        if i != cds.icds as usize {
            let c = &gff.cds[tr.cds[i]];
            let len = c
                .len
                .wrapping_sub(parent.rbeg)
                .wrapping_sub(parent.rlen as u32)
                .wrapping_add(c.beg) as i32;
            if len > 0 {
                s.extend(sub(
                    &aux.ref_,
                    pad + parent.rbeg as i64 + parent.rlen as i64 - tr_beg,
                    len as i64,
                ));
            }
        }
        i += 1;
        while i < cds.icds as usize {
            let c = &gff.cds[tr.cds[i]];
            s.extend(sub(&aux.ref_, pad + c.beg as i64 - tr_beg, c.len as i64));
            i += 1;
        }
        if parent.icds == child.icds {
            let len = sp
                .ref_beg
                .wrapping_sub(parent.rbeg)
                .wrapping_sub(parent.rlen as u32) as i32;
            if len < 0 {
                return Ok(1);
            }
            s.extend(sub(
                &aux.ref_,
                pad + parent.rbeg as i64 + parent.rlen as i64 - tr_beg,
                len as i64,
            ));
        } else {
            s.extend(sub(
                &aux.ref_,
                pad + cds.beg as i64 - tr_beg,
                sp.ref_beg as i64 - cds.beg as i64,
            ));
        }
    }
    s.extend_from_slice(&sp.kalt[dbeg..]);
    child.seq = s;
    child.sbeg = cds.pos.wrapping_add(sp.ref_beg.wrapping_sub(cds.beg));
    child.rbeg = sp.ref_beg;
    child.rlen = sp.kref.len() as i32;
    child.ty = HAP_CDS;
    child.rec = Some((rec.id, rec.pos));
    child.csq = sp.csq;
    child.dlen = rec.alleles[ial].len() as i32 - rec.alleles[0].len() as i32;
    child.var = format!("{}>{}", rec.alleles[0], rec.alleles[ial]);
    if child.rbeg.wrapping_add(child.rlen as u32) > cds.beg.wrapping_add(cds.len) {
        child.ty = HAP_SSS;
        if child.csq == 0 {
            child.csq |= CSQ_CODING_SEQUENCE;
        }
    }
    Ok(0)
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct HStack {
    pub node: usize,
    pub ichild: i64,
    pub dlen: i32,
    pub slen: usize,
}

/// Traversal state shared by all transcripts (`hap_t`).
#[derive(Default)]
pub(crate) struct HapState {
    pub stack: Vec<HStack>,
    /// Spliced haplotype sequence on the reference strand.
    pub sseq: Vec<u8>,
    /// Translated variable part of the haplotype / of the reference.
    pub tseq: Vec<u8>,
    pub tref: Vec<u8>,
    pub sbeg: u32,
    pub upstream_stop: bool,
}

impl HapState {
    fn node2soff(&self, aux: &TscriptAux, i: usize) -> i64 {
        let n = &aux.nodes[self.stack[i].node];
        self.stack[i].slen as i64 - (n.rlen as i64 + n.dlen as i64)
    }
    fn node2sbeg(&self, aux: &TscriptAux, i: usize) -> i64 {
        self.sbeg as i64 + self.node2soff(aux, i)
    }
    fn node2send(&self, i: usize) -> i64 {
        self.sbeg as i64 + self.stack[i].slen as i64
    }
    fn node2rbeg(&self, aux: &TscriptAux, i: usize) -> i64 {
        aux.nodes[self.stack[i].node].sbeg as i64
    }
    fn node2rend(&self, aux: &TscriptAux, i: usize) -> i64 {
        let n = &aux.nodes[self.stack[i].node];
        n.sbeg as i64 + n.rlen as i64
    }
    fn node2rpos(&self, aux: &TscriptAux, i: usize) -> u32 {
        aux.nodes[self.stack[i].node].rec.map(|r| r.1).unwrap_or(0)
    }
}

/// `hap_finalize`: walk the tree, build each haplotype's spliced sequence and
/// translate it piecewise into independent consequences.
pub(crate) fn hap_finalize(
    ctx: &mut Ctx,
    hap: &mut HapState,
    aux: &mut TscriptAux,
    tr_i: usize,
) -> Result<()> {
    let tr = &ctx.gff.tscripts[tr_i];
    if aux.sref.is_empty() {
        tscript_splice_ref(aux, tr, &ctx.gff.cds);
    }
    let sref = std::mem::take(&mut aux.sref);
    let nsref = sref.len() as i64;
    let pad = N_REF_PAD as i64;
    let strand = tr.strand;
    hap.stack.clear();
    hap.stack.push(HStack {
        node: 0,
        ichild: -1,
        dlen: 0,
        slen: 0,
    });
    hap.sseq.clear();
    hap.tseq.clear();
    let mut istack: i64 = 0;
    let res: Result<()> = (|| {
        while istack >= 0 {
            let is = istack as usize;
            let node_i = hap.stack[is].node;
            let nchild = aux.nodes[node_i].children.len() as i64;
            hap.stack[is].ichild += 1;
            if hap.stack[is].ichild >= nchild {
                istack -= 1;
                continue;
            }
            let child_i = aux.nodes[node_i].children[hap.stack[is].ichild as usize];
            istack += 1;
            hap.stack.truncate(istack as usize);
            let parent_slen = hap.stack[is].slen;
            let parent_dlen = hap.stack[is].dlen;
            hap.sseq.truncate(parent_slen);
            if aux.nodes[child_i].ty == HAP_CDS {
                hap.sseq.extend_from_slice(&aux.nodes[child_i].seq);
            }
            hap.stack.push(HStack {
                node: child_i,
                ichild: -1,
                dlen: parent_dlen + aux.nodes[child_i].dlen,
                slen: hap.sseq.len(),
            });
            if aux.nodes[child_i].nend == 0 {
                continue;
            }
            let is = istack as usize;
            let sseq_m = nsref - 2 * pad + hap.stack[is].dlen as i64;
            hap.upstream_stop = false;
            hap.sbeg = aux.nodes[hap.stack[1].node].sbeg;
            if aux.nodes[hap.stack[is].node].ty == HAP_SSS {
                bail!("csq: internal error, splice node as a haplotype leaf");
            }
            let mut dlen: i32 = 0;
            let mut indel = false;
            if strand == STRAND_FWD {
                let mut ibeg: i64 = -1;
                let mut i = 0usize;
                while i < is {
                    i += 1;
                    let n_dlen = aux.nodes[hap.stack[i].node].dlen;
                    dlen += n_dlen;
                    if n_dlen != 0 {
                        indel = true;
                    }
                    if i < is {
                        if dlen % 3 != 0 {
                            if ibeg == -1 {
                                ibeg = i as i64;
                            }
                            continue;
                        }
                        let mut icur = hap.node2sbeg(aux, i);
                        let inext = hap.node2sbeg(aux, i + 1);
                        if n_dlen > 0 {
                            icur += n_dlen as i64;
                        } else if n_dlen < 0 {
                            icur += 1;
                        }
                        if icur / 3 == inext / 3 {
                            if ibeg == -1 {
                                ibeg = i as i64;
                            }
                            continue;
                        }
                    }
                    if ibeg < 0 {
                        ibeg = i as i64;
                    }
                    let ib = ibeg as usize;
                    let ioff = hap.node2soff(aux, ib);
                    let icur = hap.node2sbeg(aux, ib);
                    let rbeg = hap.node2rbeg(aux, ib);
                    let rend = hap.node2rend(aux, i);
                    let mut fill = dlen % 3 != 0;
                    let alt: Vec<u8> = if !hap.sseq.is_empty() {
                        sub(&hap.sseq, ioff, hap.stack[i].slen as i64 - ioff)
                    } else {
                        fill = false;
                        Vec::new()
                    };
                    hap.tseq = cds_translate(
                        &sref,
                        &alt,
                        sseq_m as usize,
                        icur as usize,
                        rbeg as usize,
                        rend as usize,
                        strand,
                        fill,
                    )?;
                    let rseq = sub(&sref, pad + rbeg, hap.node2rend(aux, i) - rbeg);
                    hap.tref = cds_translate(
                        &sref,
                        &rseq,
                        (nsref - 2 * pad) as usize,
                        rbeg as usize,
                        rbeg as usize,
                        rend as usize,
                        strand,
                        fill,
                    )?;
                    hap_add_csq(ctx, hap, aux, tr_i, child_i, nsref, 0, ib, i, dlen, indel)?;
                    ibeg = -1;
                    dlen = 0;
                    indel = false;
                }
            } else {
                let mut ibeg: i64 = -1;
                let mut i = is + 1;
                while i > 1 {
                    i -= 1;
                    let n_dlen = aux.nodes[hap.stack[i].node].dlen;
                    dlen += n_dlen;
                    if n_dlen != 0 {
                        indel = true;
                    }
                    if i > 1 {
                        if dlen % 3 != 0 {
                            if ibeg == -1 {
                                ibeg = i as i64;
                            }
                            continue;
                        }
                        let mut icur = sseq_m - 1 - hap.node2sbeg(aux, i);
                        let mut inext = sseq_m - 1 - hap.node2sbeg(aux, i - 1);
                        if n_dlen > 0 {
                            icur += n_dlen as i64 - 1;
                        } else if n_dlen < 0 {
                            icur -= n_dlen as i64;
                        }
                        let p_dlen = aux.nodes[hap.stack[i - 1].node].dlen;
                        if p_dlen > 0 {
                            inext -= p_dlen as i64;
                        }
                        if icur / 3 == inext / 3 {
                            if ibeg == -1 {
                                ibeg = i as i64;
                            }
                            continue;
                        }
                    }
                    if ibeg < 0 {
                        ibeg = i as i64;
                    }
                    let ib = ibeg as usize;
                    let ioff = hap.node2soff(aux, i);
                    let icur = hap.node2sbeg(aux, i);
                    let rbeg = hap.node2rbeg(aux, i);
                    let rend = hap.node2rend(aux, ib);
                    let mut fill = dlen % 3 != 0;
                    let alt: Vec<u8> = if !hap.sseq.is_empty() {
                        sub(&hap.sseq, ioff, hap.stack[ib].slen as i64 - ioff)
                    } else {
                        fill = false;
                        Vec::new()
                    };
                    hap.tseq = cds_translate(
                        &sref,
                        &alt,
                        sseq_m as usize,
                        icur as usize,
                        rbeg as usize,
                        rend as usize,
                        strand,
                        fill,
                    )?;
                    let rseq = sub(&sref, pad + rbeg, hap.node2rend(aux, ib) - rbeg);
                    hap.tref = cds_translate(
                        &sref,
                        &rseq,
                        (nsref - 2 * pad) as usize,
                        rbeg as usize,
                        rbeg as usize,
                        rend as usize,
                        strand,
                        fill,
                    )?;
                    hap_add_csq(
                        ctx, hap, aux, tr_i, child_i, nsref, sseq_m, i, ib, dlen, indel,
                    )?;
                    ibeg = -1;
                    dlen = 0;
                    indel = false;
                }
            }
        }
        Ok(())
    })();
    aux.sref = sref;
    res
}

/// `hap_add_csq`: classify the translated window `[ibeg, iend]` of the
/// current haplotype and attach the consequence to `node`.
#[allow(clippy::too_many_arguments)]
fn hap_add_csq(
    ctx: &mut Ctx,
    hap: &mut HapState,
    aux: &mut TscriptAux,
    tr_i: usize,
    node: usize,
    nsref: i64,
    tlen: i64,
    ibeg: usize,
    iend: usize,
    dlen: i32,
    indel: bool,
) -> Result<()> {
    let gff = ctx.gff;
    let tr = &gff.tscripts[tr_i];
    let pad = N_REF_PAD as i64;
    let ref_node = if tr.strand == STRAND_FWD { ibeg } else { iend };
    let ref_rec = aux.nodes[hap.stack[ref_node].node]
        .rec
        .expect("variant node has a record");
    let mut csq = Csq::new(
        ref_rec.1,
        Vcsq::new(gff, tr_i, 0, aux.nodes[node].vcf_ial as u32),
    );
    let mut rm_csq = 0u32;
    for i in ibeg..=iend {
        csq.ty.ty |= aux.nodes[hap.stack[i].node].csq & CSQ_COMPOUND;
    }
    if dlen == 0 && indel {
        csq.ty.ty |= CSQ_INFRAME_ALTERING;
    }
    let has_upstream_stop = hap.upstream_stop;
    let beg_is_sss = aux.nodes[hap.stack[ibeg].node].ty == HAP_SSS;
    if !beg_is_sss {
        truncate_at_stop(&mut hap.tref);
        if truncate_at_stop(&mut hap.tseq) {
            hap.upstream_stop = true;
        }
        let tref_last = hap.tref.last().copied().unwrap_or(0);
        let tseq_last = hap.tseq.last().copied().unwrap_or(0);
        if csq.ty.ty & CSQ_STOP_LOST != 0 {
            if tref_last == b'*' && tref_last == tseq_last {
                rm_csq |= CSQ_STOP_LOST;
                csq.ty.ty |= CSQ_STOP_RETAINED;
            } else if tref_last != b'*' {
                if tseq_last == b'*' {
                    rm_csq |= CSQ_STOP_GAINED;
                    csq.ty.ty |= CSQ_STOP_RETAINED;
                } else {
                    csq.ty.ty |= CSQ_INCOMPLETE_CDS;
                }
            }
        }
        if csq.ty.ty & CSQ_START_LOST != 0 && hap.tref.first().copied() != Some(b'M') {
            rm_csq |= CSQ_START_LOST;
            csq.ty.ty &= !CSQ_START_LOST;
        }
        if dlen != 0 {
            if dlen % 3 != 0 {
                csq.ty.ty |= CSQ_FRAMESHIFT_VARIANT;
            } else if dlen < 0 {
                csq.ty.ty |= CSQ_INFRAME_DELETION;
            } else {
                csq.ty.ty |= CSQ_INFRAME_INSERTION;
            }
            if tref_last != b'*' && tseq_last == b'*' {
                csq.ty.ty |= CSQ_STOP_GAINED;
            }
        } else {
            let mut aa_change = false;
            for i in 0..hap.tref.len() {
                let t = hap.tref[i];
                let s = hap.tseq.get(i).copied().unwrap_or(0);
                if t == s {
                    continue;
                }
                aa_change = true;
                if t == b'*' {
                    csq.ty.ty |= CSQ_STOP_LOST;
                } else if s == b'*' {
                    csq.ty.ty |= CSQ_STOP_GAINED;
                } else {
                    csq.ty.ty |= CSQ_MISSENSE_VARIANT;
                }
            }
            if !aa_change {
                csq.ty.ty |= CSQ_SYNONYMOUS_VARIANT;
            }
        }
    }
    if ibeg != iend
        && csq.ty.ty & (CSQ_INFRAME_DELETION | CSQ_INFRAME_INSERTION | CSQ_INFRAME_ALTERING) != 0
        && hap.tseq.last().copied() == Some(b'*')
    {
        rm_csq |= CSQ_INFRAME_DELETION | CSQ_INFRAME_INSERTION | CSQ_INFRAME_ALTERING;
        csq.ty.ty |= CSQ_FRAMESHIFT_VARIANT | CSQ_STOP_GAINED;
    }
    if has_upstream_stop {
        csq.ty.ty |= CSQ_UPSTREAM_STOP;
    }
    csq.ty.ty &= !rm_csq;
    if beg_is_sss {
        let n = &aux.nodes[hap.stack[ibeg].node];
        csq.ty.ty |= n.csq & !rm_csq;
        csq.ty.ref_ = n.rec;
        csq.ty.biotype = tr.ty;
        let rec_id = n.rec.map(|r| r.0).unwrap_or(0);
        ctx.csq_push(&mut csq, rec_id)?;
        aux.nodes[node].csq_list.push(csq);
        return Ok(());
    }
    let mut s = String::new();
    let aa_rbeg = if tr.strand == STRAND_FWD {
        hap.node2rbeg(aux, ibeg) / 3 + 1
    } else {
        (nsref - 2 * pad - hap.node2rend(aux, iend)) / 3 + 1
    };
    let aa_sbeg = if tr.strand == STRAND_FWD {
        hap.node2sbeg(aux, ibeg) / 3 + 1
    } else {
        (tlen - hap.node2send(iend)) / 3 + 1
    };
    s.push('|');
    s.push_str(&aa_rbeg.to_string());
    kprint_aa_prediction(
        ctx.opts.brief_predictions,
        aa_rbeg.max(0) as usize,
        &hap.tref,
        &mut s,
    );
    if csq.ty.ty & CSQ_SYNONYMOUS_VARIANT == 0 {
        s.push('>');
        s.push_str(&aa_sbeg.to_string());
        kprint_aa_prediction(
            ctx.opts.brief_predictions,
            aa_sbeg.max(0) as usize,
            &hap.tseq,
            &mut s,
        );
    }
    s.push('|');
    for i in ibeg..=iend {
        if i > ibeg {
            s.push('+');
        }
        s.push_str(&(hap.node2rpos(aux, i) + 1).to_string());
        s.push_str(&aux.nodes[hap.stack[i].node].var);
    }
    csq.ty.vstr = Some(s.clone());
    ctx.csq_push(&mut csq, ref_rec.0)?;
    let main_ty = csq.ty.ty;
    aux.nodes[node].csq_list.push(csq);
    for i in ibeg..=iend {
        let (n_csq, n_rec) = {
            let n = &aux.nodes[hap.stack[i].node];
            (n.csq, n.rec.expect("variant node has a record"))
        };
        if n_csq & !CSQ_COMPOUND != 0 {
            let mut tmp = Csq::new(
                n_rec.1,
                Vcsq::new(gff, tr_i, n_csq & !CSQ_COMPOUND & !rm_csq, 0),
            );
            tmp.ty.vstr = Some(s.clone());
            ctx.csq_push(&mut tmp, n_rec.0)?;
            aux.nodes[node].csq_list.push(tmp);
        }
        if i != ref_node && (main_ty & CSQ_COMPOUND != 0 || n_csq & !CSQ_COMPOUND == 0) {
            let mut tmp = Csq::new(
                n_rec.1,
                Vcsq::new(gff, tr_i, CSQ_PRINTED_UPSTREAM | n_csq, 0),
            );
            tmp.ty.ref_ = Some(ref_rec);
            ctx.csq_push(&mut tmp, n_rec.0)?;
            aux.nodes[node].csq_list.push(tmp);
        }
    }
    Ok(())
}
