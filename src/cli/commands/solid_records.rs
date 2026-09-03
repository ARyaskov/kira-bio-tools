//! Build noodles `RecordBuf`s straight from the aligner's scored batches,
//! bypassing the SAM text round trip.
//!
//! Must agree with the aligner's SAM emitter (`kira_ls_aligner::io::SamFormatter`)
//! field for field; the two paths are cross-checked in
//! `tests/solid_record_bridge.rs`.
//!
//! Only what the fused consumers read is materialised: name, flags, coordinates,
//! CIGAR and mate fields for sorting and markdup; sequence, qualities and MAPQ
//! for mpileup; `NM` for the NM-aware quality weighting. `MD`, `AS`, `XS` and
//! `RG` are emitted by the text path but read by nothing here, and `MD` is a
//! per-record `String` — skipping them is most of this module's memory win.
//!
//! Indel left-normalization is not repeated here: kira-ls-aligner 0.4.6 applies
//! it inside the pipeline, before pairing, so the scored batches already carry
//! canonical CIGARs.

use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record::{Flags, MappingQuality};
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence};

use kira_ls_aligner::pipeline::stage5_scoring::ScoredBatch;
use kira_ls_aligner::types::{Alignment, CigarKind, MateInfo, ReadRecord};

/// Phred+33 ASCII → raw Phred score, the encoding `RecordBuf` holds.
#[inline]
fn phred_from_ascii(c: u8) -> u8 {
    c.saturating_sub(33)
}

/// Complement of an IUPAC base, mirroring the aligner's emitter.
#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// SAM FLAG for a mapped alignment. Mirrors `kira_ls_aligner::io::sam_flag`.
fn mapped_flags(aln: &Alignment) -> u16 {
    let mut flag = 0u16;
    if aln.mate.is_paired {
        flag |= 0x1;
        if aln.mate.is_proper_pair {
            flag |= 0x2;
        }
        if aln.mate.mate_is_unmapped {
            flag |= 0x8;
        }
        if aln.mate.mate_is_rev {
            flag |= 0x20;
        }
        if aln.mate.is_first_in_pair {
            flag |= 0x40;
        }
        if aln.mate.is_second_in_pair {
            flag |= 0x80;
        }
    }
    if aln.is_rev {
        flag |= 0x10;
    }
    if aln.is_secondary {
        flag |= 0x100;
    }
    if aln.is_supplementary {
        flag |= 0x800;
    }
    flag
}

/// SAM FLAG for an unmapped read. Mirrors `SamFormatter::append_unmapped_with_mate`.
fn unmapped_flags(mate: Option<&MateInfo>) -> u16 {
    let mut flag: u16 = 0x4;
    if let Some(m) = mate {
        if m.is_paired {
            flag |= 0x1;
        }
        if m.mate_is_unmapped {
            flag |= 0x8;
        } else if m.mate_is_rev {
            flag |= 0x20;
        }
        if m.is_first_in_pair {
            flag |= 0x40;
        }
        if m.is_second_in_pair {
            flag |= 0x80;
        }
    }
    flag
}

fn cigar_kind(k: CigarKind) -> Kind {
    match k {
        CigarKind::Match => Kind::Match,
        CigarKind::Ins => Kind::Insertion,
        CigarKind::Del => Kind::Deletion,
        CigarKind::SoftClip => Kind::SoftClip,
        CigarKind::Skipped => Kind::Skip,
    }
}

/// Lengths of the leading and trailing soft clips, mirroring the emitter's
/// `terminal_softclips`.
fn terminal_softclips(ops: &[kira_ls_aligner::types::CigarOp]) -> (u32, u32) {
    let lead = match ops.first() {
        Some(op) if op.op == CigarKind::SoftClip => op.len,
        _ => 0,
    };
    let trail = match ops.last() {
        Some(op) if op.op == CigarKind::SoftClip && ops.len() > 1 => op.len,
        _ => 0,
    };
    (lead, trail)
}

/// Query bases consumed by a CIGAR — the emitter's `cigar_query_consumed`.
fn cigar_query_consumed(ops: &[kira_ls_aligner::types::CigarOp]) -> u32 {
    ops.iter()
        .filter(|op| {
            matches!(
                op.op,
                CigarKind::Match | CigarKind::Ins | CigarKind::SoftClip
            )
        })
        .fold(0u32, |n, op| n.saturating_add(op.len))
}

/// Convert one mapped alignment. Returns `None` when the CIGAR does not account
/// for exactly the read's bases — the text emitter writes `*` there, which the
/// sorter and pileup cannot use, so the record is dropped instead of being
/// carried through in a shape that silently misplaces bases.
fn mapped_record(read: &ReadRecord, aln: &Alignment) -> Option<RecordBuf> {
    let seq_len = read.seq.len() as u32;
    // The consumed-length check runs against the soft-clipped form, since that
    // is what the aligner produced.
    if cigar_query_consumed(&aln.cigar) != seq_len {
        return None;
    }

    // bwa-mem hard-clips supplementary segments unless `-Y`, which the fused
    // aligner leaves at its default: the terminal clips become `H` and their
    // bases leave SEQ/QUAL, so a chimeric read's sequence is stored once, on
    // its primary record.
    let hard_clip = aln.is_supplementary;
    let (clip_lead, clip_trail) = if hard_clip {
        terminal_softclips(&aln.cigar)
    } else {
        (0, 0)
    };
    let keep_end = (seq_len as usize).saturating_sub(clip_trail as usize);
    let keep_start = (clip_lead as usize).min(keep_end);

    // SAM stores the sequence as it aligns to the forward strand; the clip
    // range indexes that oriented sequence, so it is applied after flipping.
    let seq: Vec<u8> = if aln.is_rev {
        read.seq.iter().rev().map(|&b| complement(b)).collect()
    } else {
        read.seq.clone()
    };
    let seq = seq[keep_start..keep_end].to_vec();
    // FASTQ carries Phred+33 ASCII; `RecordBuf` stores raw Phred values. Passing
    // the ASCII through unchanged silently inflates every base quality.
    let qual: Vec<u8> = match read.qual.as_ref() {
        Some(q) if aln.is_rev => q.iter().rev().map(|&c| phred_from_ascii(c)).collect(),
        Some(q) => q.iter().map(|&c| phred_from_ascii(c)).collect(),
        None => Vec::new(),
    };
    let qual = if qual.is_empty() { qual } else { qual[keep_start..keep_end].to_vec() };

    let last_op = aln.cigar.len().saturating_sub(1);
    let cigar_ops: Vec<Op> = aln
        .cigar
        .iter()
        .enumerate()
        .map(|(i, op)| {
            let terminal_clip = op.op == CigarKind::SoftClip && (i == 0 || i == last_op);
            let kind = if hard_clip && terminal_clip { Kind::HardClip } else { cigar_kind(op.op) };
            Op::new(kind, op.len as usize)
        })
        .collect();
    let cigar = Cigar::from(cigar_ops);

    let mut builder = RecordBuf::builder()
        .set_name(read.id.as_bytes())
        .set_flags(Flags::from(mapped_flags(aln)))
        .set_reference_sequence_id(aln.ref_id as usize)
        .set_cigar(cigar)
        .set_sequence(Sequence::from(seq))
        .set_quality_scores(QualityScores::from(qual));

    // SAM is 1-based; `Position` rejects 0, the unmapped case handled elsewhere.
    if let Some(pos) = noodles_core::Position::new(aln.ref_start as usize + 1) {
        builder = builder.set_alignment_start(pos);
    }
    if let Some(mq) = MappingQuality::new(aln.mapq) {
        builder = builder.set_mapping_quality(mq);
    }
    if aln.mate.is_paired
        && let Some(mate_ref) = aln.mate.mate_ref_id
        && !aln.mate.mate_is_unmapped
    {
        builder = builder.set_mate_reference_sequence_id(mate_ref as usize);
        if let Some(mpos) = noodles_core::Position::new(aln.mate.mate_pos as usize + 1) {
            builder = builder.set_mate_alignment_start(mpos);
        }
        builder = builder.set_template_length(aln.mate.tlen);
    }

    // NM only — see the module docs.
    let mut data = Data::default();
    data.insert(Tag::EDIT_DISTANCE, Value::from(aln.nm as i32));
    builder = builder.set_data(data);

    Some(builder.build())
}

/// Convert one unmapped read.
fn unmapped_record(read: &ReadRecord, mate: Option<&MateInfo>) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_name(read.id.as_bytes())
        .set_flags(Flags::from(unmapped_flags(mate)))
        .set_sequence(Sequence::from(read.seq.clone()))
        .set_quality_scores(QualityScores::from(
            read.qual
                .as_ref()
                .map(|q| q.iter().map(|&c| phred_from_ascii(c)).collect::<Vec<u8>>())
                .unwrap_or_default(),
        ));
    if let Some(m) = mate
        && let Some(mate_ref) = m.mate_ref_id
    {
        builder = builder.set_mate_reference_sequence_id(mate_ref as usize);
        if let Some(mpos) = noodles_core::Position::new(m.mate_pos as usize + 1) {
            builder = builder.set_mate_alignment_start(mpos);
        }
    }
    builder.build()
}

/// Append one scored batch's records to `out`.
///
/// `max_alignments` caps reported primary/secondary records per read exactly as
/// the SAM path's `retain_reported_alignments` does; supplementary records are
/// never dropped.
pub fn append_batch(batch: ScoredBatch, max_alignments: usize, out: &mut Vec<RecordBuf>) {
    let ScoredBatch {
        reads,
        alignments,
        unmapped_mate_info,
        ..
    } = batch;

    for ((read, alns), mate) in reads
        .iter()
        .zip(alignments.iter())
        .zip(unmapped_mate_info.iter())
    {
        if alns.is_empty() {
            out.push(unmapped_record(read, mate.as_ref()));
            continue;
        }
        let mut reported = 0usize;
        for aln in alns {
            if !aln.is_supplementary {
                if max_alignments > 0 && reported >= max_alignments {
                    continue;
                }
                reported += 1;
            }
            if let Some(rec) = mapped_record(read, aln) {
                out.push(rec);
            }
        }
    }
}
