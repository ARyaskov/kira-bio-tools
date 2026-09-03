//! One thread budget per process. `--threads N` (or `-@ N`) sets it; every
//! pool (BGZF compression/decompression, rayon, BAM decoding) derives its
//! worker count from it instead of `num_cpus` so the pools do not oversubscribe.

use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};

static BUDGET: AtomicUsize = AtomicUsize::new(0);
static RAYON: OnceLock<()> = OnceLock::new();

fn hardware_threads() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

/// Set the budget (0 restores the hardware default). Call before any pool is built.
pub fn set_budget(n: usize) {
    BUDGET.store(n, Ordering::Relaxed);
}

/// Threads this process may use in total.
pub fn budget() -> usize {
    match BUDGET.load(Ordering::Relaxed) {
        0 => hardware_threads(),
        n => n,
    }
}

/// Workers for a BGZF compression pool (the writer thread is extra).
pub fn compress_workers() -> usize {
    budget().clamp(1, 16)
}

/// Workers for a BGZF decompression pool; decompression is cheap relative to
/// parsing, so a quarter of the budget keeps the reader ahead of the consumer.
pub fn decompress_workers() -> usize {
    (budget() / 4).clamp(1, 8)
}

/// Workers for BAM decoding, where decompression dominates the load phase.
pub fn bam_workers() -> usize {
    budget().clamp(1, 16)
}

/// Build the global rayon pool from the budget (first call wins; later calls
/// and pools already built elsewhere are left alone).
pub fn init_rayon() {
    RAYON.get_or_init(|| {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(budget()).build_global();
    });
}

/// `--threads N`, `--threads=N` or `-@ N` from a raw argv, if present.
pub fn budget_from_argv<S: AsRef<str>>(argv: &[S]) -> Option<usize> {
    let mut it = argv.iter().map(|s| s.as_ref());
    while let Some(a) = it.next() {
        if a == "--threads" || a == "-@" {
            return it.next().and_then(|v| v.parse().ok());
        }
        if let Some(v) = a.strip_prefix("--threads=") {
            return v.parse().ok();
        }
    }
    None
}
