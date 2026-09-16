//! Short-tandem-repeat aware indel error model (`mpileup --hp-indel`).
//!
//! Polymerase slippage makes indel errors in a repeat tract far more likely
//! than elsewhere, and the rate grows with the number of copies of the
//! repeat unit. Instead of a fixed quality penalty per homopolymer run, the
//! model gives each indel candidate a per-read error probability from the
//! unit length and copy number of the tract it sits in, in the spirit of
//! Dindel's and Strelka2's log-linear models: a base rate for non-repetitive
//! sequence, multiplied by a per-copy growth factor that is steepest for
//! homopolymers. The probability becomes the quality of the synthetic indel
//! observations fed to the genotype likelihoods, and the tracts where the
//! model expects errors decide how much support a call needs.

/// The repeat context of one indel candidate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StrContext {
    /// Length of the repeat unit the indel adds or removes (1 = homopolymer).
    pub unit: usize,
    /// Copies of that unit in the reference tract the indel sits in.
    pub copies: u32,
    /// Whether the indel is an insertion.
    pub insertion: bool,
}

/// Per-read probability of a spurious indel observation in this context.
pub fn error_rate(ctx: &StrContext) -> f64 {
    // Non-repetitive base rates: deletions are somewhat more frequent in
    // Illumina data than insertions.
    let base: f64 = if ctx.insertion { 1.0e-4 } else { 1.5e-4 };
    // Growth per additional copy of the unit; homopolymers slip most.
    let growth: f64 = match ctx.unit {
        0 | 1 => 1.5,
        2 => 1.35,
        _ => 1.25,
    };
    let extra = ctx.copies.saturating_sub(1) as f64;
    (base * growth.powf(extra)).min(0.5)
}

/// Phred quality of an indel observation in this context, clamped to
/// `[8, 40]` so a single tract never dominates and clean sequence keeps the
/// default indel quality.
pub fn indel_quality(ctx: &StrContext) -> u8 {
    let q = -10.0 * error_rate(ctx).log10();
    q.round().clamp(8.0, 40.0) as u8
}

/// Tracts where errors are expected: the model's rate is above the
/// per-read rate a default caller tolerates without extra support.
pub fn error_prone(ctx: &StrContext) -> bool {
    error_rate(ctx) > 1.0e-3
}

/// Smallest period of `seq` (its repeat unit length), `seq.len()` when it
/// does not repeat.
pub fn repeat_unit_len(seq: &[u8]) -> usize {
    let n = seq.len();
    for p in 1..n {
        if n.is_multiple_of(p) && (p..n).all(|i| seq[i].eq_ignore_ascii_case(&seq[i - p])) {
            return p;
        }
    }
    n.max(1)
}

/// Copies of `unit` tiled through `tract`, forwards from its start.
fn copies_forward(tract: &[u8], unit: &[u8]) -> u32 {
    let u = unit.len();
    if u == 0 {
        return 0;
    }
    let mut n = 0u32;
    let mut i = 0usize;
    while i + u <= tract.len() && tract[i..i + u].eq_ignore_ascii_case(unit) {
        n += 1;
        i += u;
    }
    n
}

/// Copies of `unit` ending at the end of `tract`, walking backwards.
fn copies_backward(tract: &[u8], unit: &[u8]) -> u32 {
    let u = unit.len();
    if u == 0 {
        return 0;
    }
    let mut n = 0u32;
    let mut end = tract.len();
    while end >= u && tract[end - u..end].eq_ignore_ascii_case(unit) {
        n += 1;
        end -= u;
    }
    n
}

/// The context of an indel whose inserted or deleted sequence is `seq`,
/// anchored so that `seq` sits right after `left` in the reference and
/// `right` follows it. `left` is the reference before the anchor point
/// (its end adjoins the event), `right` the reference from the event on.
pub fn str_context(seq: &[u8], left: &[u8], right: &[u8], insertion: bool) -> StrContext {
    if seq.is_empty() {
        return StrContext {
            unit: 1,
            copies: 0,
            insertion,
        };
    }
    let u = repeat_unit_len(seq);
    let unit = &seq[..u];
    // Copies of the unit in the reference on both sides of the event; the
    // deleted bases themselves are part of the reference tract.
    let copies = copies_backward(left, unit) + copies_forward(right, unit);
    StrContext {
        unit: u,
        copies,
        insertion,
    }
}

#[cfg(test)]
#[path = "../../tests/unit/bam_str_model.rs"]
mod tests;
