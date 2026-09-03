//! Streaming reader: one decoder thread per BAM pushes `LiveRead`s into a
//! bounded channel; the pileup engine pulls, so no BAM is ever fully in RAM.

use crate::bam::pileup::{LiveRead, build_live_from_bam};
use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, bounded};
use noodles_bam as bam;
use noodles_sam as sam;
use std::fs::File;
use std::num::NonZeroUsize;
use std::path::Path;
use std::thread::JoinHandle;

const CHANNEL_DEPTH: usize = 4096;

/// Per-read hook run on the decoder thread (flag filtering, BAQ); returns
/// `false` to drop the read. Receives the contig names of the file.
pub type ReadHook = Box<dyn FnMut(&mut LiveRead, &[String]) -> bool + Send>;

pub struct StreamingBam {
    pub header: sam::Header,
    pub samples: Vec<String>,
    pub ref_names: Vec<String>,
    pub ref_lengths: Vec<u64>,
    pub rx: Receiver<LiveRead>,
    /// Decoder thread; its error (a truncated or corrupt BAM) is reported
    /// when joined after the walk.
    pub handle: Option<JoinHandle<Result<()>>>,
}

impl StreamingBam {
    pub fn open<P: AsRef<Path>>(p: P, sample_idx: usize) -> Result<Self> {
        Self::open_with(p, sample_idx, 1, None)
    }

    /// `workers` BGZF decompression threads feed the decoder.
    pub fn open_with<P: AsRef<Path>>(p: P, sample_idx: usize, workers: usize, mut hook: Option<ReadHook>) -> Result<Self> {
        let path = p.as_ref().to_path_buf();
        let f = File::open(&path).with_context(|| format!("open BAM {:?}", path))?;
        let bgzf = noodles_bgzf::io::MultithreadedReader::with_worker_count(NonZeroUsize::new(workers.max(1)).unwrap(), f);
        let mut reader = bam::io::Reader::from(bgzf);
        let header = reader.read_header().context("read BAM header")?;
        let samples = crate::bam::reader::extract_samples_helper(&header);
        let ref_names: Vec<String> = header
            .reference_sequences()
            .keys()
            .map(|k| String::from_utf8_lossy(k).into_owned())
            .collect();
        let ref_lengths: Vec<u64> = header.reference_sequences().values().map(|s| usize::from(s.length()) as u64).collect();
        let names = ref_names.clone();

        let (tx, rx) = bounded::<LiveRead>(CHANNEL_DEPTH);
        let handle = std::thread::Builder::new()
            .name(format!("bam-stream-{}", path.display()))
            .spawn(move || -> Result<()> {
                for r in reader.records() {
                    let rec = r.with_context(|| format!("read BAM record from {}", path.display()))?;
                    let Some(mut lr) = build_live_from_bam(&rec, sample_idx) else { continue };
                    if let Some(h) = hook.as_mut() {
                        if !h(&mut lr, &names) {
                            continue;
                        }
                    }
                    if tx.send(lr).is_err() {
                        break;
                    }
                }
                Ok(())
            })
            .context("spawn bam stream thread")?;

        Ok(Self { header, samples, ref_names, ref_lengths, rx, handle: Some(handle) })
    }
}
