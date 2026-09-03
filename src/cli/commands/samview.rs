use crate::bam::{BamWriter, writer::OutputKind};
use crate::cli::args::SamViewArgs;
use anyhow::{Context, Result, bail};
use noodles_bam as bam;
use noodles_sam as sam;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn cmd_samview(args: SamViewArgs) -> Result<()> {
    let ext = args.input.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    if ext == "cram" { bail!("samview: CRAM input not yet supported (use samtools view)"); }

    let use_index = args.region.is_some() && bai_path(&args.input).is_some();

    let header;
    let mut iter: Box<dyn Iterator<Item = std::io::Result<bam::Record>>>;
    if use_index {
        let mut ir = bam::io::indexed_reader::Builder::default()
            .build_from_path(&args.input).with_context(|| format!("open indexed BAM {:?}", args.input))?;
        header = ir.read_header().context("read BAM header")?;
        let reg_str = args.region.as_deref().unwrap();
        let reg: noodles_core::Region = reg_str.parse()
            .map_err(|e| anyhow::anyhow!("parse region {reg_str:?}: {e}"))?;
        let query = ir.query(&header, &reg).context("query region")?;
        let records: Vec<bam::Record> = query.records()
            .collect::<std::io::Result<Vec<_>>>()
            .context("read indexed BAM records")?;
        iter = Box::new(records.into_iter().map(Ok));
        eprintln!("[samview] using index for region '{}'", reg_str);
    } else {
        let f = File::open(&args.input).with_context(|| format!("open {:?}", args.input))?;
        let mut reader = bam::io::Reader::new(f);
        header = reader.read_header().context("read BAM header")?;
        let records: Vec<bam::Record> = reader.records()
            .collect::<std::io::Result<Vec<_>>>()
            .context("read BAM records")?;
        iter = Box::new(records.into_iter().map(Ok));
    }

    let mut bam_writer: Option<BamWriter> = None;
    let mut text_writer: Option<Box<dyn Write>> = None;

    if args.bam || args.cram {
        if args.cram {
            eprintln!("[samview] note: CRAM output uses embedded reference sequences (no -f required); records will be encoded with reference-free strategy where supported");
        }
        let (default_name, kind) = if args.cram { ("out.cram", OutputKind::Cram) } else { ("out.bam", OutputKind::Bam) };
        let out_path = args.output.clone().unwrap_or_else(|| std::path::PathBuf::from(default_name));
        bam_writer = Some(BamWriter::create_kind(&out_path, &header, kind)?);
    } else {
        let writer: Box<dyn Write> = match args.output.as_deref() {
            Some(p) => Box::new(BufWriter::with_capacity(1 << 20, File::create(p)?)),
            None => Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout())),
        };
        text_writer = Some(writer);
    }

    if args.with_header || args.header_only {
        if let Some(w) = text_writer.as_mut() {
            write_text_header(w, &header)?;
        }
    }
    if args.header_only { return finish(bam_writer, text_writer); }

    let region_filter = if use_index { None } else { parse_region(args.region.as_deref(), &header)? };
    let mut n_count: u64 = 0;
    let mut rng_state: u64 = 0x9E3779B97F4A7C15;

    for r in iter.by_ref() {
        let rec = r.context("read BAM record")?;
        let flags = u16::from(rec.flags());
        if args.require_flags != 0 && (flags & args.require_flags) != args.require_flags { continue; }
        if args.exclude_flags != 0 && (flags & args.exclude_flags) != 0 { continue; }
        let mq: u32 = rec.mapping_quality().map(|m| u8::from(m) as u32).unwrap_or(0);
        if mq < args.min_mq { continue; }

        if let Some((rid, beg, end)) = region_filter {
            let r_rid = rec.reference_sequence_id().transpose().ok().flatten();
            if r_rid != Some(rid) { continue; }
            let pos = rec.alignment_start().transpose().ok().flatten().map(|p| usize::from(p) as u32);
            let Some(pos) = pos else { continue; };
            if pos > end { continue; }
            let read_end = read_end_pos(&rec).unwrap_or(pos);
            if read_end < beg { continue; }
        }

        if args.subsample < 1.0 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let frac = ((rng_state >> 11) as f64) / ((1u64 << 53) as f64);
            if frac > args.subsample { continue; }
        }

        if args.count { n_count += 1; continue; }

        if let Some(bw) = bam_writer.as_mut() {
            bw.write_record(&rec)?;
        } else if let Some(tw) = text_writer.as_mut() {
            write_sam_record(tw, &header, &rec)?;
        }
    }

    if args.count {
        let mut out: Box<dyn Write> = match args.output.as_deref() {
            Some(p) => Box::new(BufWriter::new(File::create(p)?)),
            None => Box::new(BufWriter::new(std::io::stdout())),
        };
        writeln!(out, "{}", n_count)?;
    }
    finish(bam_writer, text_writer)
}

fn finish(bam_writer: Option<BamWriter>, text_writer: Option<Box<dyn Write>>) -> Result<()> {
    if let Some(b) = bam_writer { b.finish()?; }
    if let Some(mut t) = text_writer { t.flush()?; }
    Ok(())
}

fn write_text_header<W: Write>(w: &mut W, header: &sam::Header) -> Result<()> {
    let mut s = sam::io::Writer::new(Vec::new());
    s.write_header(header).context("serialize SAM header")?;
    let bytes = s.into_inner();
    w.write_all(&bytes)?;
    Ok(())
}

fn write_sam_record<W: Write>(w: &mut W, header: &sam::Header, rec: &bam::Record) -> Result<()> {
    let qname = rec.name().map(|n| std::str::from_utf8(n.as_ref()).unwrap_or("*").to_string()).unwrap_or_else(|| "*".to_string());
    let flag = u16::from(rec.flags());
    let rname = match rec.reference_sequence_id().transpose().ok().flatten() {
        Some(rid) => header.reference_sequences().get_index(rid)
            .map(|(k, _)| std::str::from_utf8(k.as_ref()).unwrap_or("*").to_string())
            .unwrap_or_else(|| "*".to_string()),
        None => "*".to_string(),
    };
    let pos = rec.alignment_start().transpose().ok().flatten().map(|p| usize::from(p)).unwrap_or(0);
    let mapq = rec.mapping_quality().map(|m| u8::from(m)).unwrap_or(255);
    let cigar = {
        let mut s = String::new();
        use noodles_sam::alignment::record::cigar::op::Kind;
        for op in rec.cigar().iter() {
            let op = op.map_err(|e| anyhow::anyhow!("cigar parse: {e}"))?;
            let c = match op.kind() {
                Kind::Match => 'M', Kind::Insertion => 'I', Kind::Deletion => 'D',
                Kind::Skip => 'N', Kind::SoftClip => 'S', Kind::HardClip => 'H',
                Kind::Pad => 'P', Kind::SequenceMatch => '=', Kind::SequenceMismatch => 'X',
            };
            s.push_str(&format!("{}{}", op.len(), c));
        }
        if s.is_empty() { "*".to_string() } else { s }
    };
    let mrnm = match rec.mate_reference_sequence_id().transpose().ok().flatten() {
        Some(rid) => header.reference_sequences().get_index(rid)
            .map(|(k, _)| std::str::from_utf8(k.as_ref()).unwrap_or("*").to_string())
            .unwrap_or_else(|| "*".to_string()),
        None => "*".to_string(),
    };
    let mpos = rec.mate_alignment_start().transpose().ok().flatten().map(|p| usize::from(p)).unwrap_or(0);
    let isize_v = rec.template_length();
    let seq: Vec<u8> = rec.sequence().iter().collect();
    let qual_raw: Vec<u8> = rec.quality_scores().iter().collect();
    let seq_str = if seq.is_empty() { "*".to_string() } else { String::from_utf8_lossy(&seq).into_owned() };
    let qual_str = if qual_raw.is_empty() { "*".to_string() } else {
        qual_raw.iter().map(|q| (q + 33) as char).collect()
    };
    writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize_v, seq_str, qual_str)?;
    Ok(())
}

fn parse_region(s: Option<&str>, header: &sam::Header) -> Result<Option<(usize, u32, u32)>> {
    let Some(s) = s else { return Ok(None); };
    let (chr, beg, end) = crate::regions::parse_region_spec(s)?;
    let rid = header.reference_sequences().get_index_of(chr.as_bytes())
        .ok_or_else(|| anyhow::anyhow!("unknown contig {chr:?}"))?;
    Ok(Some((rid, beg, end)))
}

fn read_end_pos(rec: &bam::Record) -> Option<u32> {
    use noodles_sam::alignment::record::cigar::op::Kind;
    let start = rec.alignment_start().transpose().ok().flatten().map(|p| usize::from(p) as u32)?;
    let mut span = 0u32;
    for op in rec.cigar().iter() {
        let op = op.ok()?;
        if matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip) {
            span += op.len() as u32;
        }
    }
    Some(start + span - 1)
}

fn bai_path(p: &Path) -> Option<std::path::PathBuf> {
    let s = p.to_string_lossy();
    for ext in &[".bai", ".csi"] {
        let cand = std::path::PathBuf::from(format!("{}{}", s, ext));
        if cand.exists() { return Some(cand); }
    }
    None
}
