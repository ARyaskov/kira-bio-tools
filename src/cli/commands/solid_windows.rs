//! Windowed execution for the fused `solid` pipeline.
//!
//! Alignments are spilled, as they are produced, into one temporary BAM per
//! reference window; each window is then loaded, sorted, deduplicated and called
//! on its own, so peak memory is one window's depth instead of the whole run.
//!
//! Two boundary problems have to be handled explicitly, because both fail
//! silently:
//!
//! * **Pileup coverage.** A read starting just before a window boundary still
//!   covers positions after it. Every record is therefore written to each window
//!   its reference span touches, and a window only *emits* calls inside its own
//!   coordinate range — so the extra copies add coverage without double-counting
//!   sites.
//!
//! * **Duplicate marking.** `mark_duplicates_in_memory` picks one winning
//!   template per position key and flags the rest. Without an `ms` tag the score
//!   is the read's own base-quality sum, which differs between R1 and R2, so if
//!   the two mates of a template landed in different windows each would pick its
//!   winner from a different score population and could disagree. Templates are
//!   therefore also written to the window of their *leftmost* mate, keeping a
//!   template's records together for the dedup pass.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::{self as sam, Header};

/// Identifies one reference window: `(reference id, window index)`.
pub type WindowId = (usize, u32);

/// Writes alignment records into per-window temporary BAMs. Writers are opened
/// lazily, so only windows that receive reads get a file.
pub struct WindowSpiller {
    dir: PathBuf,
    window_len: u32,
    writers: BTreeMap<WindowId, BamFileWriter>,
    counts: BTreeMap<WindowId, u64>,
    unplaced: Option<BamFileWriter>,
    unplaced_count: u64,
}

/// `bam::io::Writer` wraps the sink in a BGZF writer, so the concrete type is
/// one layer deeper than the file it is built from.
type BamFileWriter = bam::io::Writer<noodles_bgzf::io::Writer<BufWriter<File>>>;

fn window_path(dir: &Path, id: WindowId) -> PathBuf {
    dir.join(format!("w_{}_{}.bam", id.0, id.1))
}

impl WindowSpiller {
    pub fn new(dir: PathBuf, window_len: u32) -> Self {
        Self {
            dir,
            window_len,
            writers: BTreeMap::new(),
            counts: BTreeMap::new(),
            unplaced: None,
            unplaced_count: 0,
        }
    }

    /// Reference span of a record as `(ref_id, start, end)`, 0-based half-open.
    /// `None` for records with no placement.
    fn span(rec: &RecordBuf) -> Option<(usize, usize, usize)> {
        use noodles_sam::alignment::record::cigar::op::Kind;
        let ref_id = rec.reference_sequence_id()?;
        let start = usize::from(rec.alignment_start()?) - 1;
        let span: usize = rec
            .cigar()
            .as_ref()
            .iter()
            .filter(|op| {
                matches!(
                    op.kind(),
                    Kind::Match
                        | Kind::Deletion
                        | Kind::Skip
                        | Kind::SequenceMatch
                        | Kind::SequenceMismatch
                )
            })
            .map(|op| op.len())
            .sum();
        Some((ref_id, start, start + span.max(1)))
    }

    /// Spill one record into every window that needs it.
    pub fn push(&mut self, header: &Header, rec: &RecordBuf) -> Result<()> {
        let Some((ref_id, start, end)) = Self::span(rec) else {
            return self.push_unplaced(header, rec);
        };

        let wlen = self.window_len as usize;
        let first_covered = start / wlen;
        // `end` is exclusive, so a record ending exactly on a boundary does not
        // reach into the next window.
        let last_covered = (end - 1) / wlen;

        // Keep both mates of a template together so the dedup pass ranks them
        // against the same population.
        let mut template = start;
        if rec.mate_reference_sequence_id() == Some(ref_id)
            && let Some(mp) = rec.mate_alignment_start()
        {
            template = template.min(usize::from(mp) - 1);
        }
        let template_window = template / wlen;

        let lo = first_covered.min(template_window);
        let hi = last_covered.max(template_window);
        for w in lo..=hi {
            let id = (ref_id, w as u32);
            // Split the borrows: `dir` is read while `writers` is mutated.
            let Self { dir, writers, .. } = self;
            if !writers.contains_key(&id) {
                let path = window_path(dir, id);
                let file = File::create(&path)
                    .with_context(|| format!("create window spill {}", path.display()))?;
                let mut w = bam::io::Writer::new(BufWriter::with_capacity(1 << 20, file));
                w.write_header(header).context("write spill header")?;
                writers.insert(id, w);
            }
            writers
                .get_mut(&id)
                .expect("writer just inserted")
                .write_alignment_record(header, rec as &dyn sam::alignment::Record)
                .context("write window record")?;
            *self.counts.entry(id).or_insert(0) += 1;
        }
        Ok(())
    }

    /// Unmapped/unplaced records take no part in any pileup, but are kept so a
    /// run's record count stays complete.
    fn push_unplaced(&mut self, header: &Header, rec: &RecordBuf) -> Result<()> {
        if self.unplaced.is_none() {
            let path = self.dir.join("w_unplaced.bam");
            let file =
                File::create(&path).with_context(|| format!("create {}", path.display()))?;
            let mut w = bam::io::Writer::new(BufWriter::with_capacity(1 << 16, file));
            w.write_header(header).context("write spill header")?;
            self.unplaced = Some(w);
        }
        self.unplaced
            .as_mut()
            .expect("writer just inserted")
            .write_alignment_record(header, rec as &dyn sam::alignment::Record)
            .context("write unplaced record")?;
        self.unplaced_count += 1;
        Ok(())
    }

    pub fn unplaced_records(&self) -> u64 {
        self.unplaced_count
    }

    /// Flush every writer and return the populated windows in coordinate order.
    pub fn finish(mut self) -> Result<Vec<SpilledWindow>> {
        for (_, w) in std::mem::take(&mut self.writers) {
            let mut w = w;
            w.try_finish().context("finish window spill")?;
        }
        if let Some(mut w) = self.unplaced.take() {
            w.try_finish().context("finish unplaced spill")?;
        }
        let dir = self.dir.clone();
        let window_len = self.window_len;
        Ok(self
            .counts
            .iter()
            .map(|(&id, &records)| SpilledWindow {
                id,
                path: window_path(&dir, id),
                records,
                window_len,
            })
            .collect())
    }
}

/// One populated window on disk.
#[derive(Clone, Debug)]
pub struct SpilledWindow {
    pub id: WindowId,
    pub path: PathBuf,
    pub records: u64,
    window_len: u32,
}

impl SpilledWindow {
    /// Window width in bases.
    pub fn window_len(&self) -> usize {
        self.window_len as usize
    }

    /// The reference range this window is responsible for emitting, formatted as
    /// mpileup's `-r` expects (1-based, inclusive).
    pub fn region_spec(&self, header: &Header) -> Option<String> {
        let (ref_id, widx) = self.id;
        let (name, seq) = header.reference_sequences().get_index(ref_id)?;
        let len = usize::from(seq.length());
        let start = widx as usize * self.window_len as usize;
        if start >= len {
            return None;
        }
        let end = ((widx as usize + 1) * self.window_len as usize).min(len);
        Some(format!("{name}:{}-{end}", start + 1))
    }

    /// Read this window's records back.
    pub fn load(&self) -> Result<(Header, Vec<RecordBuf>)> {
        let file = File::open(&self.path)
            .with_context(|| format!("open window spill {}", self.path.display()))?;
        let mut reader = bam::io::Reader::new(file);
        let header = reader.read_header().context("read spill header")?;
        let mut out = Vec::with_capacity(self.records as usize);
        for rec in reader.record_bufs(&header) {
            out.push(rec.context("read spill record")?);
        }
        Ok((header, out))
    }

    pub fn remove(&self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

/// Choose a directory for the spill files.
pub fn resolve_tmpdir(explicit: Option<&Path>, output: &Path) -> PathBuf {
    match explicit {
        Some(d) => d.to_path_buf(),
        None => output
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from(".")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record::cigar::{Op, op::Kind};
    use noodles_sam::alignment::record_buf::Cigar;
    use noodles_sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;

    fn header() -> Header {
        Header::builder()
            .add_reference_sequence(
                "chrA",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .build()
    }

    fn rec(start: usize, len: usize, mate: Option<usize>) -> RecordBuf {
        let mut b = RecordBuf::builder()
            .set_name(format!("r{start}").as_bytes())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::new(start + 1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, len)]));
        if let Some(m) = mate {
            b = b
                .set_mate_reference_sequence_id(0)
                .set_mate_alignment_start(noodles_core::Position::new(m + 1).unwrap());
        }
        b.build()
    }

    fn windows_touched(window_len: u32, r: &RecordBuf) -> Vec<u32> {
        let dir = tempfile::tempdir().unwrap();
        let mut s = WindowSpiller::new(dir.path().to_path_buf(), window_len);
        s.push(&header(), r).unwrap();
        s.finish().unwrap().iter().map(|w| w.id.1).collect()
    }

    #[test]
    fn a_read_inside_one_window_is_written_once() {
        assert_eq!(windows_touched(100, &rec(10, 20, None)), [0]);
    }

    /// A read overlapping the boundary must reach the next window, or coverage
    /// there silently drops.
    #[test]
    fn a_read_spanning_a_boundary_reaches_both_windows() {
        assert_eq!(windows_touched(100, &rec(90, 20, None)), [0, 1]);
    }

    #[test]
    fn a_read_ending_exactly_on_the_boundary_stays_in_one_window() {
        // [80, 100) touches window 0 only; `end` is exclusive.
        assert_eq!(windows_touched(100, &rec(80, 20, None)), [0]);
    }

    /// Mates that straddle a boundary must both be present in the leftmost
    /// window so duplicate marking ranks them against the same population.
    #[test]
    fn mates_straddling_a_boundary_are_kept_together() {
        // Read at 250 whose mate starts at 90: template window is 0, covered
        // window is 2, so it is written to 0, 1 and 2.
        assert_eq!(windows_touched(100, &rec(250, 20, Some(90))), [0, 1, 2]);
        // The mate itself sits in window 0 already.
        assert_eq!(windows_touched(100, &rec(90, 20, Some(250))), [0, 1]);
    }

    #[test]
    fn unplaced_records_go_to_the_unplaced_file_and_no_window() {
        let dir = tempfile::tempdir().unwrap();
        let mut s = WindowSpiller::new(dir.path().to_path_buf(), 100);
        let unmapped = RecordBuf::builder().set_name(b"u").build();
        s.push(&header(), &unmapped).unwrap();
        assert_eq!(s.unplaced_records(), 1);
        assert!(s.finish().unwrap().is_empty());
    }

    #[test]
    fn region_spec_is_one_based_and_clamped_to_the_contig() {
        let h = header();
        let dir = tempfile::tempdir().unwrap();
        let mut s = WindowSpiller::new(dir.path().to_path_buf(), 400);
        s.push(&h, &rec(10, 20, None)).unwrap();
        s.push(&h, &rec(950, 20, None)).unwrap();
        let ws = s.finish().unwrap();
        assert_eq!(ws[0].region_spec(&h).unwrap(), "chrA:1-400");
        // Last window is clipped to the contig length, not 1200.
        assert_eq!(ws.last().unwrap().region_spec(&h).unwrap(), "chrA:801-1000");
    }

    #[test]
    fn records_round_trip_through_the_spill() {
        let h = header();
        let dir = tempfile::tempdir().unwrap();
        let mut s = WindowSpiller::new(dir.path().to_path_buf(), 1000);
        for start in [10usize, 40, 700] {
            s.push(&h, &rec(start, 30, None)).unwrap();
        }
        let ws = s.finish().unwrap();
        assert_eq!(ws.len(), 1);
        let (_, loaded) = ws[0].load().unwrap();
        assert_eq!(loaded.len(), 3);
        let starts: Vec<usize> = loaded
            .iter()
            .map(|r| usize::from(r.alignment_start().unwrap()) - 1)
            .collect();
        assert_eq!(starts, [10, 40, 700]);
    }
}
