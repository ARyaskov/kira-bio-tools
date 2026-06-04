//! Streaming reader: spawns one thread per BAM that decodes records and
//! pushes LiveRead into a bounded channel. Pileup engine pulls — no full
//! Vec<Record> in RAM. For huge BAMs (multi-GB).

use crate::bam::pileup::{LiveRead, build_live_from_bam};
use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, bounded};
use noodles_bam as bam;
use noodles_sam as sam;
use std::fs::File;
use std::path::Path;
use std::thread::JoinHandle;

const CHANNEL_DEPTH: usize = 4096;

pub struct StreamingBam {
    pub header: sam::Header,
    pub samples: Vec<String>,
    pub ref_names: Vec<String>,
    pub rx: Receiver<LiveRead>,
    _handle: JoinHandle<()>,
}

impl StreamingBam {
    pub fn open<P: AsRef<Path>>(p: P, sample_idx: usize) -> Result<Self> {
        let path = p.as_ref().to_path_buf();
        let f = File::open(&path).with_context(|| format!("open BAM {:?}", path))?;
        let mut reader = bam::io::Reader::new(f);
        let header = reader.read_header().context("read BAM header")?;
        let samples = crate::bam::reader::extract_samples_helper(&header);
        let ref_names: Vec<String> = header.reference_sequences().keys()
            .map(|k| std::str::from_utf8(k).unwrap_or("").to_string()).collect();

        let (tx, rx) = bounded::<LiveRead>(CHANNEL_DEPTH);
        let handle = std::thread::Builder::new()
            .name(format!("bam-stream-{}", path.display()))
            .spawn(move || {
                for r in reader.records() {
                    match r {
                        Ok(rec) => {
                            if let Some(lr) = build_live_from_bam(&rec, sample_idx) {
                                if tx.send(lr).is_err() { break; }
                            }
                        }
                        Err(_) => break,
                    }
                }
            })
            .context("spawn bam stream thread")?;

        Ok(Self { header, samples, ref_names, rx, _handle: handle })
    }
}
