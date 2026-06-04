//! Variant-key composition for the ANI MPH.
//!
//! Single canonical key constructor used by both builder and lookup sites.
//! Previously each call site (5 of them: finalize, tab_processing,
//! structs/ani/lookup, cuda/lookup, opencl/v2) hand-rolled the same XOR
//! recipe — that was both duplication and a correctness trap:
//!
//! **Legacy bug:** the old formula `(chr<<32 | pos) ^ fxhash(ref) ^ fxhash(alt)`
//! is XOR-commutative. `A>T` and `T>A` at the same position produce the
//! *same* key, silently overwriting each other in the MPH and making one of
//! the two variants unreachable at lookup time. Real-world impact: clinical
//! databases that carry both strand polarities lose annotations on the
//! second-written variant.
//!
//! **Fix:** stream-hash (chr, pos, ref_len, ref, separator, alt_len, alt) so
//! the algebra is fully order-dependent. Length prefixes + the 0xFF separator
//! also defeat synthetic collisions from refs/alts that share a prefix
//! (e.g. `REF=AA ALT=A` vs `REF=A ALT=AA`).

use std::hash::Hasher;

use fxhash::FxHasher;

/// Composes a deterministic, non-commutative `u64` key from a variant
/// coordinate. Used identically at build and lookup time — any change here
/// must be made in lock-step on both sides.
#[inline]
pub fn make_variant_key(chr_id: u32, pos: u32, ref_bytes: &[u8], alt_bytes: &[u8]) -> u64 {
    let mut h = FxHasher::default();
    h.write_u32(chr_id);
    h.write_u32(pos);
    h.write_u16(ref_bytes.len() as u16);
    h.write(ref_bytes);
    // Hard separator: guards against ref/alt boundary ambiguity.
    h.write_u8(0xFF);
    h.write_u16(alt_bytes.len() as u16);
    h.write(alt_bytes);
    h.finish()
}

#[cfg(test)]
#[path = "../../../../tests/unit/annotate_structs_ani_key.rs"]
mod tests;
