pub mod annotation;
pub mod column_spec;
pub mod field_metadata;
pub mod lookup;
pub mod merge_info;
pub mod merge_info_helpers;
pub mod threads;
pub mod vcf_parsing;

pub use annotation::*;
pub use column_spec::*;
pub use field_metadata::*;
pub use lookup::*;
pub use merge_info::*;
pub use threads::*;
pub use vcf_parsing::*;

use anyhow::Result;
use crossbeam_channel::bounded;
use std::path::Path;
use std::sync::Arc;
use std::thread;
use std::time::Instant;

use super::constants::*;
use super::reader::{StreamingVcfReader, VcfAnnotationReader};
use crate::annotate::structs::ani::AniIndex;
use crate::util::{detect_format, VcfFormat};

pub fn annotate_vcf_ani_v2(
    db: &Path,
    input: &Path,
    output: &Path,
    columns: &[String],
) -> Result<()> {
    let timing = std::env::var("KIRA_BT_TIMING").is_ok();
    let debug = std::env::var("KIRA_BT_DEBUG").is_ok() || timing;
    let start = Instant::now();

    let column_specs = ColumnSpec::parse_all(columns);
    let field_order: Vec<String> = column_specs.iter().map(|c| c.dst_key.clone()).collect();

    let num_threads = rayon::current_num_threads() / 2;
    if debug {
        eprintln!("[annotate] Using {} CPU threads", num_threads);
        eprintln!("[annotate] Batch size: {} lines", BATCH_SIZE);
        eprintln!(
            "[annotate] Column specs: {:?}",
            column_specs
                .iter()
                .map(|c| format!("{}{}", c.mode, c.key))
                .collect::<Vec<_>>()
        );
    }

    let load_start = Instant::now();
    let ani = AniIndex::open(db)?;
    if debug {
        eprintln!(
            "[annotate] ANI load: {:.3}s",
            load_start.elapsed().as_secs_f64()
        );
    }

    let field_meta = load_and_infer_metadata(&ani, debug)?;

    let input_format = detect_format(input)?;
    let use_bgzf = matches!(input_format, VcfFormat::Bgzf);
    if debug {
        eprintln!(
            "[annotate] Input: {:?}, Output BGZF: {}",
            input_format, use_bgzf
        );
    }

    let reader_start = Instant::now();
    let input_reader = VcfAnnotationReader::open(input)?;
    let streaming_reader = StreamingVcfReader::new(input_reader);
    let (headers, mut reader) = streaming_reader.into_headers_and_self()?;
    if debug {
        eprintln!(
            "[annotate] Reader init: {:.3}s, headers: {}",
            reader_start.elapsed().as_secs_f64(),
            headers.len()
        );
    }

    let merged_headers = merge_annotation_headers(&headers, &ani)?;

    let (read_tx, read_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);
    let (work_tx, work_rx) = bounded::<Vec<String>>(CHANNEL_DEPTH);

    let ani_clone = Arc::new(ani);
    let field_meta_clone = Arc::new(field_meta);
    let field_order_arc = Arc::new(field_order);
    let column_specs_arc = Arc::new(column_specs);

    let worker = thread::spawn(move || {
        worker_thread(
            read_rx,
            work_tx,
            ani_clone.clone(),
            field_meta_clone.clone(),
            field_order_arc.clone(),
            column_specs_arc.clone(),
            num_threads,
        )
    });

    let output_clone = output.to_path_buf();
    let writer = thread::spawn(move || {
        writer_thread(work_rx, merged_headers, &output_clone, use_bgzf, timing)
    });

    read_batches(&mut reader, read_tx, timing)?;

    worker.join().unwrap()?;
    writer.join().unwrap()?;

    if debug {
        eprintln!(
            "[annotate] Total time: {:.3}s",
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

fn read_batches(
    reader: &mut StreamingVcfReader,
    read_tx: crossbeam_channel::Sender<Vec<String>>,
    timing: bool,
) -> Result<()> {
    let mut batch = Vec::with_capacity(BATCH_SIZE);
    let mut total_lines = 0usize;
    let mut last_report = Instant::now();

    while let Some(line) = reader.read_line()? {
        batch.push(line);

        if batch.len() >= BATCH_SIZE {
            total_lines += batch.len();
            if read_tx
                .send(std::mem::replace(
                    &mut batch,
                    Vec::with_capacity(BATCH_SIZE),
                ))
                .is_err()
            {
                break;
            }

            if timing && last_report.elapsed().as_secs() >= 2 {
                eprintln!("[annotate] Progress: {} lines read", total_lines);
                last_report = Instant::now();
            }
        }
    }

    if !batch.is_empty() {
        total_lines += batch.len();
        let _ = read_tx.send(batch);
    }

    Ok(())
}

fn merge_annotation_headers(vcf_headers: &[String], _ani: &AniIndex) -> Result<Vec<String>> {
    Ok(vcf_headers.to_vec())
}
