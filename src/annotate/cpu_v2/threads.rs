use anyhow::Result;
use crossbeam_channel::{Receiver, Sender};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use crate::annotate::constants::OUTPUT_BUFFER_SIZE;
use crate::annotate::cpu_v2::annotation::annotate_line;
use crate::annotate::cpu_v2::column_spec::ColumnSpec;
use crate::annotate::structs::ani::AniIndex;
use crate::annotate::structs::annotate_mode::AnnotateMode;
use crate::annotate::structs::bundle::FieldNumber;
use crate::bgzf::BgzfWriter;

pub fn worker_thread(
    rx: Receiver<Vec<String>>,
    tx: Sender<Vec<String>>,
    ani: Arc<AniIndex>,
    field_meta: Arc<HashMap<String, FieldNumber>>,
    column_specs: Arc<Vec<ColumnSpec>>,
    sample_map: Arc<Vec<Option<usize>>>,
    info_overwrite_all: bool,
    format_overwrite_all: bool,
    num_threads: usize,
) -> Result<()> {
    let column_modes: Arc<Vec<(String, AnnotateMode)>> = Arc::new(
        column_specs
            .iter()
            .map(|c| (c.key.clone(), c.mode))
            .collect(),
    );

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    while let Ok(batch) = rx.recv() {
        let annotated: Vec<String> = pool.install(|| {
            batch
                .par_iter()
                .map(|line| {
                    annotate_line(
                        line,
                        &ani,
                        &field_meta,
                        &column_modes,
                        &sample_map,
                        info_overwrite_all,
                        format_overwrite_all,
                    )
                })
                .collect()
        });

        if tx.send(annotated).is_err() {
            break;
        }
    }

    Ok(())
}

pub fn writer_thread(
    rx: Receiver<Vec<String>>,
    headers: Vec<String>,
    output: &Path,
    use_bgzf: bool,
    timing: bool,
) -> Result<()> {
    let start = Instant::now();
    let mut lines_written = 0usize;
    let mut bytes_written = 0usize;

    if use_bgzf {
        let mut writer = BgzfWriter::create(output)?;

        for h in headers {
            writeln!(writer, "{}", h)?;
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for line in batch {
                writeln!(writer, "{}", line)?;
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }

        writer.finish()?;
    } else {
        let file = File::create(output)?;
        let mut writer = BufWriter::with_capacity(OUTPUT_BUFFER_SIZE, file);

        for h in headers {
            writeln!(writer, "{}", h)?;
            bytes_written += h.len() + 1;
        }

        while let Ok(batch) = rx.recv() {
            for line in batch {
                writeln!(writer, "{}", line)?;
                bytes_written += line.len() + 1;
                lines_written += 1;
            }
        }

        writer.flush()?;
    }

    if timing {
        let elapsed = start.elapsed().as_secs_f64();
        let mb_sec = (bytes_written as f64 / 1_048_576.0) / elapsed;
        eprintln!(
            "[annotate] Write complete: {} lines, {:.1} MB/s",
            lines_written, mb_sec
        );
    }

    Ok(())
}
