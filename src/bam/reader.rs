use crate::bam::pileup::{LiveRead, build_live_from_bam, build_live_from_cram};
use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;
use std::fs::File;
use std::path::{Path, PathBuf};

pub struct BamReader {
    pub records_buf: Vec<LiveRead>,
    pub header: sam::Header,
    pub samples: Vec<String>,
    pub ref_names: Vec<String>,
    /// Contig lengths from the `@SQ` lines, parallel to `ref_names`.
    pub ref_lengths: Vec<u64>,
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
        let ref_lengths = header_lengths(&header);
        Self {
            records_buf,
            header,
            samples,
            ref_names,
            ref_lengths,
            source_path: PathBuf::from("<memory>"),
        }
    }
}

fn header_lengths(h: &sam::Header) -> Vec<u64> {
    h.reference_sequences().values().map(|s| usize::from(s.length()) as u64).collect()
}

pub fn apply_hmm_baq_to_reads(
    reads: &mut [LiveRead],
    ref_names: &[String],
    fasta: &(impl FastaLike + Sync),
) {
    apply_hmm_baq_to_reads_masked(reads, ref_names, fasta, None, "off");
}

/// Same as `apply_hmm_baq_to_reads`, but skips reads where `needs_baq[i] == false`.
/// A read with no mismatches and no indels is unchanged by BAQ — skipping is a
/// pure-correctness no-op that saves 30-50% of HMM cost on typical WGS data.
/// `nm_weight` is the `--nm-weight` spec.
pub fn apply_hmm_baq_to_reads_masked(
    reads: &mut [LiveRead],
    ref_names: &[String],
    fasta: &(impl FastaLike + Sync),
    needs_baq: Option<&[bool]>,
    nm_weight: &str,
) {
    use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
    use rayon::prelude::*;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

    let total = reads.len() as u64;
    if total == 0 { return; }

    // The NM profile is chosen once from the observed read length (a platform proxy).
    {
        let n = reads.len().min(1000);
        let sum: usize = reads.iter().take(n).map(|r| r.seq().len()).sum();
        crate::bam::baq::init_nm_profile(nm_weight, if n > 0 { sum / n } else { 0 });
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
                baq_read(lr, ref_names, fasta);
            }
            counter.fetch_add(batch.len() as u64, Ordering::Relaxed);
        });

    stop.store(true, Ordering::Relaxed);
    let _ = ticker.join();
    pb.finish_with_message("BAQ done");
}

/// BAQ and NM down-weighting for one read; a no-op when its contig is not in
/// the reference. The reference window reaches a read length past both
/// ends of the alignment, as `sam_prob_realn` may look that far.
pub fn baq_read(lr: &mut LiveRead, ref_names: &[String], fasta: &impl FastaLike) {
    let Some(chr) = ref_names.get(lr.ref_id) else { return };
    let span = ref_span(&lr.cigar_pairs);
    let ref_start = lr.ref_start;
    let pad = lr.seq().len() as u32 + 16;
    let lo = ref_start.saturating_sub(pad);
    let len = (ref_start - lo) as usize + span as usize + pad as usize;
    let cigar = lr.cigar_pairs.clone();
    if let Some(window) = fasta.slice(chr, lo + 1, len) {
        let (seq, qual) = lr.seq_qual_mut();
        crate::bam::baq::apply_baq_hmm(seq, qual, &cigar, window, lo, ref_start);
        let off = (ref_start - lo) as usize;
        if off < window.len() {
            crate::bam::baq::nm_weight_qual(seq, qual, &cigar, &window[off..]);
        }
    }
}

fn ref_span(cigar: &[(noodles_sam::alignment::record::cigar::op::Kind, u32)]) -> u32 {
    use noodles_sam::alignment::record::cigar::op::Kind;
    cigar.iter().filter(|(k, _)| matches!(k, Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip))
        .map(|(_, l)| *l).sum()
}

pub trait FastaLike {
    fn slice(&self, chr: &str, pos: u32, len: usize) -> Option<&[u8]>;
}

/// Reference sequences for CRAM decoding, backed by an indexed FASTA. The
/// `.fai` is read when present and computed otherwise.
pub fn fasta_repository(fasta_path: &Path) -> Result<fasta::Repository> {
    let mut fai_path = fasta_path.as_os_str().to_os_string();
    fai_path.push(".fai");
    let fai_path = PathBuf::from(fai_path);
    let index = if fai_path.exists() {
        fasta::fai::fs::read(&fai_path).with_context(|| format!("read {}", fai_path.display()))?
    } else {
        fasta::fs::index(fasta_path).with_context(|| format!("index {}", fasta_path.display()))?
    };
    let reader = fasta::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(fasta_path)
        .with_context(|| format!("open reference {}", fasta_path.display()))?;
    Ok(fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader)))
}

impl BamReader {
    pub fn open<P: AsRef<Path>>(p: P) -> Result<Self> {
        Self::open_with_reference(p, None)
    }

    /// Open a BAM or CRAM. CRAM decoding needs the reference the file was
    /// encoded against (`--fasta-ref`) unless it embeds its sequences.
    pub fn open_with_reference<P: AsRef<Path>>(p: P, reference: Option<&Path>) -> Result<Self> {
        let path = p.as_ref();
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
        if ext == "cram" { return Self::open_cram(path, reference, None); }
        if ext == "sam" { return Self::open_sam(path, None); }
        let f = File::open(path).with_context(|| format!("open BAM {:?}", path))?;
        // MultithreadedReader scales BGZF decompression across cores (~3-5× vs single-threaded).
        let bgzf_r = noodles_bgzf::io::MultithreadedReader::with_worker_count(
            std::num::NonZeroUsize::new(crate::threads::bam_workers()).unwrap(),
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
        Self::open_with_region_and_reference(p, region, None)
    }

    /// `region` is `chr`, `chr:pos` or `chr:beg-end` (1-based, inclusive).
    pub fn open_with_region_and_reference<P: AsRef<Path>>(p: P, region: &str, reference: Option<&Path>) -> Result<Self> {
        let (chrom, beg, end) = parse_region(region)?;
        Self::open_with_regions_and_reference(p, &[(chrom, beg, end)], reference)
    }

    /// Reads overlapping any of the (sorted, non-overlapping) regions, in position order.
    pub fn open_with_regions_and_reference<P: AsRef<Path>>(p: P, regions: &[(String, u32, u32)], reference: Option<&Path>) -> Result<Self> {
        let path = p.as_ref();
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
        if ext == "cram" {
            return Self::open_cram(path, reference, Some(regions));
        }
        if ext == "sam" {
            return Self::open_sam(path, Some(regions));
        }
        let mut ir = bam::io::indexed_reader::Builder::default()
            .build_from_path(path).with_context(|| format!("open indexed BAM {:?}", path))?;
        let header = ir.read_header().context("read header")?;
        let mut records_buf: Vec<LiveRead> = Vec::new();
        for (chrom, beg, end) in regions {
            let reg = make_region(chrom, *beg, *end);
            let query = ir.query(&header, &reg).with_context(|| format!("query region {chrom}:{beg}-{end}"))?;
            for r in query.records() {
                let rec = r.context("read indexed BAM record")?;
                if let Some(lr) = build_live_from_bam(&rec, 0) { records_buf.push(lr); }
            }
        }
        records_buf.sort_by_key(|lr| (lr.ref_id, lr.ref_start));
        Ok(Self::build(records_buf, header, path))
    }

    fn open_cram(path: &Path, reference: Option<&Path>, regions: Option<&[(String, u32, u32)]>) -> Result<Self> {
        let repo = match reference {
            Some(r) => fasta_repository(r)?,
            None => fasta::Repository::default(),
        };
        let mut records_buf: Vec<LiveRead> = Vec::new();
        let header = match regions {
            Some(regions) => {
                let mut ir = cram::io::indexed_reader::Builder::default()
                    .set_reference_sequence_repository(repo)
                    .build_from_path(path)
                    .with_context(|| format!("open indexed CRAM {:?} (needs {:?}.crai)", path, path))?;
                let header = ir.read_header().context("read CRAM header")?;
                for (chrom, beg, end) in regions {
                    let reg = make_region(chrom, *beg, *end);
                    let query = ir.query(&header, &reg).with_context(|| format!("query CRAM region {chrom}:{beg}-{end}"))?;
                    for r in query {
                        let rb = r.context("read CRAM record")?;
                        if let Some(lr) = build_live_from_cram(&rb, 0) { records_buf.push(lr); }
                    }
                }
                records_buf.sort_by_key(|lr| (lr.ref_id, lr.ref_start));
                header
            }
            None => {
                let mut inner = cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repo)
                    .build_from_path(path)
                    .with_context(|| format!("open CRAM {:?}", path))?;
                let header = inner.read_header().context("read CRAM header")?;
                for r in inner.records(&header) {
                    let rb = match r {
                        Ok(rb) => rb,
                        Err(e) if reference.is_none() => bail!(
                            "read CRAM record: {e}; reference-encoded CRAM needs the reference FASTA (-f/--fasta-ref)"
                        ),
                        Err(e) => return Err(e).context("read CRAM record"),
                    };
                    if let Some(lr) = build_live_from_cram(&rb, 0) { records_buf.push(lr); }
                }
                header
            }
        };
        Ok(Self::build(records_buf, header, path))
    }

    /// Plain SAM: no index, so regions are applied to the decoded records.
    fn open_sam(path: &Path, regions: Option<&[(String, u32, u32)]>) -> Result<Self> {
        let f = File::open(path).with_context(|| format!("open SAM {:?}", path))?;
        let mut reader = sam::io::Reader::new(std::io::BufReader::with_capacity(1 << 20, f));
        let header = reader.read_header().context("read SAM header")?;
        let names: Vec<String> = header.reference_sequences().keys().map(|k| std::str::from_utf8(k).unwrap_or("").to_string()).collect();
        let mut records_buf: Vec<LiveRead> = Vec::new();
        for r in reader.records() {
            let rec = r.context("read SAM record")?;
            let rb = sam::alignment::RecordBuf::try_from_alignment_record(&header, &rec).context("decode SAM record")?;
            let Some(lr) = build_live_from_cram(&rb, 0) else { continue };
            if let Some(regs) = regions {
                let chrom = names.get(lr.ref_id).map(String::as_str).unwrap_or("");
                let keep = regs.iter().any(|(c, b, e)| c == chrom && lr.ref_start < *e && lr.ref_end_cached >= *b);
                if !keep { continue; }
            }
            records_buf.push(lr);
        }
        records_buf.sort_by_key(|lr| (lr.ref_id, lr.ref_start));
        Ok(Self::build(records_buf, header, path))
    }

    fn build(records_buf: Vec<LiveRead>, header: sam::Header, p: &Path) -> Self {
        let samples = extract_samples(&header);
        let ref_names: Vec<String> = header.reference_sequences().keys()
            .map(|k| std::str::from_utf8(k).unwrap_or("").to_string()).collect();
        let ref_lengths = header_lengths(&header);
        Self { records_buf, header, samples, ref_names, ref_lengths, source_path: p.to_path_buf() }
    }
}

/// `chr`, `chr:pos` or `chr:beg-end` (1-based, inclusive; open end allowed).
fn parse_region(s: &str) -> Result<(String, u32, u32)> {
    let s = s.trim();
    if let Some((c, r)) = s.rsplit_once(':') {
        if r.chars().next().is_some_and(|ch| ch.is_ascii_digit()) {
            let num = |t: &str| crate::util::parse_coordinate(t).ok_or_else(|| anyhow::anyhow!("bad region {s:?}"));
            return match r.split_once('-') {
                Some((b, e)) => Ok((c.to_string(), (num(b)? as u32).max(1), if e.is_empty() { u32::MAX } else { num(e)? as u32 })),
                None => { let p = num(r)? as u32; Ok((c.to_string(), p.max(1), p)) }
            };
        }
    }
    Ok((s.to_string(), 1, u32::MAX))
}

fn make_region(chrom: &str, beg: u32, end: u32) -> noodles_core::Region {
    use noodles_core::Position;
    let start = Position::try_from(beg.max(1) as usize).unwrap_or(Position::MIN);
    if end == u32::MAX {
        noodles_core::Region::new(chrom, start..)
    } else {
        let stop = Position::try_from((end.max(beg)) as usize).unwrap_or(Position::MIN);
        noodles_core::Region::new(chrom, start..=stop)
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
