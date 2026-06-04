//! Pipeline sizing constants.

pub const BATCH_SIZE: usize = 1000_000;
pub const CHANNEL_DEPTH: usize = 8;
pub const OUTPUT_BUFFER_SIZE: usize = 512 * 1024 * 1024;
pub const ESTIMATE_BYTES_PER_LINE: usize = 200;

/// Target raw bytes per batch in the byte-bounded reader. Override with
/// `KIRA_BT_BATCH_BYTES=<bytes>`.
pub const BATCH_TARGET_BYTES: usize = 256 * 1024 * 1024;
/// Floor on lines/batch.
pub const BATCH_MIN_LINES: usize = 1_000;
/// Hard ceiling on lines/batch regardless of byte budget.
pub const BATCH_MAX_LINES: usize = 400_000;

/// Resolve the per-batch byte budget at runtime, honouring the env override.
#[inline]
pub fn batch_target_bytes() -> usize {
    std::env::var("KIRA_BT_BATCH_BYTES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&n| n >= 1024)
        .unwrap_or(BATCH_TARGET_BYTES)
}
