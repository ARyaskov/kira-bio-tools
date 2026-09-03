//! The fused `solid` pipeline builds `RecordBuf`s directly from the aligner's
//! scored batches instead of round-tripping through SAM text. That bridge must
//! agree with the aligner's SAM emitter field for field: a mismatch is silent,
//! since record counts still match and only the contents drift. (The first
//! version passed Phred+33 ASCII where `RecordBuf` expects raw Phred, halving the
//! chr20 variant count with no error.) Both paths run the same reads here and the
//! records are compared field by field.

use std::io::Write;

use noodles_sam::alignment::RecordBuf;
use noodles_sam::{self as sam};

use kira_bio_tools::cli::commands::solid_records::append_batch;
use kira_ls_aligner::aligner_core::Aligner;
use kira_ls_aligner::cli::commands::mem::{FusedAlignerParams, build_short_pe_aligner};
use kira_ls_aligner::index::{Index, IndexConfig};
use kira_ls_aligner::io::read_reference;

/// Deterministic pseudo-random DNA.
fn dna(len: usize, mut state: u64) -> Vec<u8> {
    let alphabet = b"ACGT";
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            alphabet[(state >> 60) as usize & 3]
        })
        .collect()
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Build a small reference plus paired reads sampled from it, with substitutions
/// and a short indel so mapped, clipped and gapped records all appear.
fn make_fixture(dir: &std::path::Path) -> (std::path::PathBuf, Vec<std::path::PathBuf>) {
    let reference = dna(200_000, 0xC0FFEE);
    let ref_path = dir.join("ref.fa");
    let mut f = std::fs::File::create(&ref_path).unwrap();
    writeln!(f, ">chrT").unwrap();
    for chunk in reference.chunks(60) {
        f.write_all(chunk).unwrap();
        f.write_all(b"\n").unwrap();
    }
    drop(f);

    let r1_path = dir.join("r1.fq");
    let r2_path = dir.join("r2.fq");
    let mut f1 = std::fs::File::create(&r1_path).unwrap();
    let mut f2 = std::fs::File::create(&r2_path).unwrap();
    let mut state = 0x1234_5678u64;
    for i in 0..600 {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let frag = 300 + (state >> 55) as usize % 100;
        let start = (state >> 20) as usize % (reference.len() - frag - 1);
        let mut r1 = reference[start..start + 150].to_vec();
        let mut r2: Vec<u8> = reference[start + frag - 150..start + frag]
            .iter()
            .rev()
            .map(|&b| complement(b))
            .collect();
        // A substitution in most reads, a deletion in some.
        if i % 3 == 0 {
            r1[40] = if r1[40] == b'A' { b'T' } else { b'A' };
        }
        if i % 7 == 0 {
            r2.remove(70);
            r2.push(b'A');
        }
        // Vary the qualities so an encoding mistake cannot hide behind a constant.
        let q1: String = (0..r1.len())
            .map(|j| (b'#' + ((j + i) % 40) as u8) as char)
            .collect();
        let q2: String = (0..r2.len())
            .map(|j| (b'#' + ((j * 3 + i) % 40) as u8) as char)
            .collect();
        writeln!(f1, "@rd{i}/1\n{}\n+\n{}", String::from_utf8_lossy(&r1), q1).unwrap();
        writeln!(f2, "@rd{i}/2\n{}\n+\n{}", String::from_utf8_lossy(&r2), q2).unwrap();
    }
    drop(f1);
    drop(f2);
    (ref_path, vec![r1_path, r2_path])
}

fn build_aligner(ref_path: &std::path::Path) -> Aligner {
    let params = FusedAlignerParams {
        reference: ref_path.to_path_buf(),
        index: None,
        threads: 2,
        num_p_threads: None,
        num_e_threads: None,
        batch_bases: 100_000, // small, so several batches are exercised
        read_group: None,
        paired: true,
        interleaved: false,
        insert_size: "0,1000,350,50".to_string(),
        n_read_files: 2,
    };
    build_short_pe_aligner(&params).unwrap().0
}

fn describe(r: &RecordBuf) -> String {
    format!(
        "name={:?} flags={:#06x} ref={:?} pos={:?} mapq={:?} cigar={:?} \
         mate_ref={:?} mate_pos={:?} tlen={} seq={:?} qual={:?}",
        r.name(),
        u16::from(r.flags()),
        r.reference_sequence_id(),
        r.alignment_start().map(usize::from),
        r.mapping_quality().map(u8::from),
        r.cigar(),
        r.mate_reference_sequence_id(),
        r.mate_alignment_start().map(usize::from),
        r.template_length(),
        r.sequence().as_ref(),
        r.quality_scores().as_ref(),
    )
}

#[test]
fn streaming_records_match_the_sam_text_path() {
    let dir = tempfile::tempdir().unwrap();
    let (ref_path, reads) = make_fixture(dir.path());

    // Path 1: SAM text, then parse back.
    let aligner = build_aligner(&ref_path);
    let index = Index::build(read_reference(&ref_path).unwrap(), IndexConfig::default());
    let sam_bytes = aligner.align_to_sam_bytes(index, &reads).unwrap();
    let mut reader = sam::io::Reader::new(std::io::Cursor::new(&sam_bytes[..]));
    let header = reader.read_header().unwrap();
    let via_text: Vec<RecordBuf> = reader
        .record_bufs(&header)
        .map(|r| r.unwrap())
        .collect();

    // Path 2: structured batches straight to records.
    let aligner = build_aligner(&ref_path);
    let index = Index::build(read_reference(&ref_path).unwrap(), IndexConfig::default());
    let mut via_records: Vec<RecordBuf> = Vec::new();
    aligner
        .align_streaming(index, &reads, |batch| {
            append_batch(batch, 1, &mut via_records);
            Ok(())
        })
        .unwrap();

    assert!(!via_text.is_empty(), "fixture produced no alignments");
    assert_eq!(
        via_text.len(),
        via_records.len(),
        "record counts differ: text={} records={}",
        via_text.len(),
        via_records.len()
    );

    for (i, (t, r)) in via_text.iter().zip(via_records.iter()).enumerate() {
        assert_eq!(t.name(), r.name(), "record {i}: name");
        assert_eq!(t.flags(), r.flags(), "record {i}: flags\n{}\n{}", describe(t), describe(r));
        assert_eq!(
            t.reference_sequence_id(),
            r.reference_sequence_id(),
            "record {i}: reference id"
        );
        assert_eq!(t.alignment_start(), r.alignment_start(), "record {i}: position");
        assert_eq!(
            t.mapping_quality(),
            r.mapping_quality(),
            "record {i}: mapping quality"
        );
        assert_eq!(t.cigar(), r.cigar(), "record {i}: cigar");
        assert_eq!(
            t.mate_reference_sequence_id(),
            r.mate_reference_sequence_id(),
            "record {i}: mate reference id"
        );
        assert_eq!(
            t.mate_alignment_start(),
            r.mate_alignment_start(),
            "record {i}: mate position"
        );
        assert_eq!(
            t.template_length(),
            r.template_length(),
            "record {i}: template length"
        );
        assert_eq!(
            t.sequence().as_ref(),
            r.sequence().as_ref(),
            "record {i}: sequence"
        );
        assert_eq!(
            t.quality_scores().as_ref(),
            r.quality_scores().as_ref(),
            "record {i}: quality scores\n{}\n{}",
            describe(t),
            describe(r)
        );
    }
}

/// The bridge deliberately drops tags nothing downstream reads, but `NM` feeds
/// the NM-aware quality weighting and must survive.
#[test]
fn nm_tag_is_preserved() {
    use noodles_sam::alignment::record::data::field::Tag;

    let dir = tempfile::tempdir().unwrap();
    let (ref_path, reads) = make_fixture(dir.path());
    let aligner = build_aligner(&ref_path);
    let index = Index::build(read_reference(&ref_path).unwrap(), IndexConfig::default());
    let mut recs: Vec<RecordBuf> = Vec::new();
    aligner
        .align_streaming(index, &reads, |batch| {
            append_batch(batch, 1, &mut recs);
            Ok(())
        })
        .unwrap();

    let mapped: Vec<&RecordBuf> = recs
        .iter()
        .filter(|r| u16::from(r.flags()) & 0x4 == 0)
        .collect();
    assert!(!mapped.is_empty(), "no mapped records");
    assert!(
        mapped.iter().all(|r| r.data().get(&Tag::EDIT_DISTANCE).is_some()),
        "every mapped record must carry NM"
    );
    // The fixture plants mismatches, so at least one NM must be non-zero —
    // otherwise a bridge that always wrote 0 would pass the check above.
    assert!(
        mapped.iter().any(|r| matches!(
            r.data().get(&Tag::EDIT_DISTANCE),
            Some(v) if v.as_int().unwrap_or(0) > 0
        )),
        "expected some non-zero NM values"
    );
}
