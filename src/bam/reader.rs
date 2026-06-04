use crate::bam::pileup::{LiveRead, build_live_from_bam, build_live_from_cram};
use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam as sam;
use std::fs::File;
use std::path::{Path, PathBuf};

pub struct BamReader {
    pub records_buf: Vec<LiveRead>,
    pub header: sam::Header,
    pub samples: Vec<String>,
    pub ref_names: Vec<String>,
    pub source_path: PathBuf,
}

impl BamReader {
    /// Construct a reader from in-memory records (no file). For fused in-process
    /// pipelines that hand sorted records straight to mpileup.
    pub fn from_parts(
        records_buf: Vec<LiveRead>,
        header: sam::Header,
        samples: Vec<String>,
        ref_names: Vec<String>,
    ) -> Self {
        Self {
            records_buf,
            header,
            samples,
            ref_names,
            source_path: PathBuf::from("<memory>"),
        }
    }
}

pub fn apply_hmm_baq_to_reads(
    reads: &mut [LiveRead],
    ref_names: &[String],
    fasta: &(impl FastaLike + Sync),
) {
    apply_hmm_baq_to_reads_masked(reads, ref_names, fasta, None);
}

/// Same as `apply_hmm_baq_to_reads`, but skips reads where `needs_baq[i] == false`.
/// A read with no mismatches and no indels is unchanged by BAQ — skipping is a
/// pure-correctness no-op that saves 30-50% of HMM cost on typical WGS data.
pub fn apply_hmm_baq_to_reads_masked(
    reads: &mut [LiveRead],
    ref_names: &[String],
    fasta: &(impl FastaLike + Sync),
    needs_baq: Option<&[bool]>,
) {
    use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
    use rayon::prelude::*;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

    let total = reads.len() as u64;
    if total == 0 { return; }

    // Auto-select the NM down-weighting profile from the observed read length (platform proxy),
    // once, before the parallel pass. KIRA_NM_WEIGHT overrides ("auto"|"off"|"F,S").
    {
        let n = reads.len().min(1000);
        let sum: usize = reads.iter().take(n).map(|r| r.seq.len()).sum();
        crate::bam::baq::init_nm_profile(if n > 0 { sum / n } else { 0 });
    }

    let no_pb = std::env::var("KIRA_BT_NO_PROGRESS").is_ok();
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::with_template("[BAQ] {bar:40.cyan/blue} {pos}/{len} reads ({per_sec}, ETA {eta})")
            .unwrap_or_else(|_| ProgressStyle::default_bar()),
    );
    pb.set_draw_target(if no_pb { ProgressDrawTarget::hidden() } else { ProgressDrawTarget::stderr_with_hz(2) });
    let counter = Arc::new(AtomicU64::new(0));
    let stop = Arc::new(AtomicBool::new(false));
    let ticker = {
        let counter = Arc::clone(&counter);
        let stop = Arc::clone(&stop);
        let pb_clone = pb.clone();
        std::thread::spawn(move || {
            while !stop.load(Ordering::Relaxed) {
                pb_clone.set_position(counter.load(Ordering::Relaxed));
                std::thread::sleep(std::time::Duration::from_millis(500));
            }
            pb_clone.set_position(counter.load(Ordering::Relaxed));
        })
    };

    const BATCH: usize = 256;
    let needs = needs_baq;

    reads
        .par_chunks_mut(BATCH)
        .enumerate()
        .for_each(|(batch_idx, batch)| {
            let batch_off = batch_idx * BATCH;
            for (i, lr) in batch.iter_mut().enumerate() {
                let global = batch_off + i;
                if let Some(mask) = needs {
                    if !mask.get(global).copied().unwrap_or(true) {
                        continue;
                    }
                }
                let Some(chr) = ref_names.get(lr.ref_id) else { continue };
                let span = ref_span(&lr.cigar_pairs);
                if let Some(slice) = fasta.slice(chr, lr.ref_start + 1, span as usize) {
                    let LiveRead { seq, qual, cigar_pairs, .. } = lr;
                    crate::bam::baq::apply_baq_hmm(&seq[..], &mut qual[..], &cigar_pairs[..], slice, 7);
                    crate::bam::baq::nm_weight_qual(&seq[..], &mut qual[..], &cigar_pairs[..], slice);
                }
            }
            counter.fetch_add(batch.len() as u64, Ordering::Relaxed);
        });

    stop.store(true, Ordering::Relaxed);
    let _ = ticker.join();
    pb.finish_with_message("BAQ done");
}

fn ref_span(cigar: &[(noodles_sam::alignment::record::cigar::op::Kind, u32)]) -> u32 {
    use noodles_sam::alignment::record::cigar::op::Kind;
    cigar.iter().filter(|(k, _)| matches!(k, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip))
        .map(|(_, l)| *l).sum()
}

pub trait FastaLike {
    fn slice(&self, chr: &str, pos: u32, len: usize) -> Option<&[u8]>;
}

impl BamReader {
    pub fn open<P: AsRef<Path>>(p: P) -> Result<Self> {
        let path = p.as_ref();
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
        if ext == "cram" { return Self::open_cram(path); }
        let f = File::open(path).with_context(|| format!("open BAM {:?}", path))?;
        // MultithreadedReader scales BGZF decompression across cores (~3-5× vs single-threaded).
        let workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).min(16);
        let bgzf_r = noodles_bgzf::io::MultithreadedReader::with_worker_count(
            std::num::NonZeroUsize::new(workers).unwrap(),
            f,
        );
        let mut inner = bam::io::Reader::from(bgzf_r);
        let header = inner.read_header().context("read BAM header")?;
        let mut records_buf: Vec<LiveRead> = Vec::with_capacity(64 * 1024);
        for r in inner.records() {
            let rec = r.context("read BAM record")?;
            if let Some(lr) = build_live_from_bam(&rec, 0) { records_buf.push(lr); }
        }
        Ok(Self::build(records_buf, header, path))
    }

    pub fn open_with_region<P: AsRef<Path>>(p: P, region: &str) -> Result<Self> {
        let path = p.as_ref();
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
        if ext == "cram" {
            eprintln!("[mpileup] note: CRAM + region query not yet supported; reading all records");
            return Self::open_cram(path);
        }
        let mut ir = bam::io::indexed_reader::Builder::default()
            .build_from_path(path).with_context(|| format!("open indexed BAM {:?}", path))?;
        let header = ir.read_header().context("read header")?;
        let reg: noodles_core::Region = region.parse()
            .map_err(|e| anyhow::anyhow!("parse region {region:?}: {e}"))?;
        let query = ir.query(&header, &reg).context("query region")?;
        let mut records_buf: Vec<LiveRead> = Vec::new();
        for r in query.records() {
            let rec = r.context("read indexed BAM record")?;
            if let Some(lr) = build_live_from_bam(&rec, 0) { records_buf.push(lr); }
        }
        Ok(Self::build(records_buf, header, path))
    }

    fn open_cram(path: &Path) -> Result<Self> {
        let f = File::open(path).with_context(|| format!("open CRAM {:?}", path))?;
        let mut inner = noodles_cram::io::Reader::new(f);
        let header = inner.read_header().context("read CRAM header")?;
        let mut records_buf: Vec<LiveRead> = Vec::new();
        for r in inner.records(&header) {
            let rb = r.context("read CRAM record")?;
            if let Some(lr) = build_live_from_cram(&rb, 0) { records_buf.push(lr); }
        }
        Ok(Self::build(records_buf, header, path))
    }

    fn build(records_buf: Vec<LiveRead>, header: sam::Header, p: &Path) -> Self {
        let samples = extract_samples(&header);
        let ref_names: Vec<String> = header.reference_sequences().keys()
            .map(|k| std::str::from_utf8(k).unwrap_or("").to_string()).collect();
        Self { records_buf, header, samples, ref_names, source_path: p.to_path_buf() }
    }
}

pub fn extract_samples_helper(h: &sam::Header) -> Vec<String> { extract_samples(h) }

fn extract_samples(h: &sam::Header) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for (_id, rg) in h.read_groups() {
        for (k, v) in rg.other_fields().iter() {
            if k.as_ref() == b"SM" {
                let s = std::str::from_utf8(v).unwrap_or("").to_string();
                if !out.contains(&s) { out.push(s); }
            }
        }
    }
    out
}
