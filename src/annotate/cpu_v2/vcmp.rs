//! Bcftools-compatible REF/ALT comparison for differently-padded indels.
//!
//! Port of `bcftools/vcmp.c` (Petr Danecek, MIT) restructured for our
//! parallel/stateless use.

/// Resolved REF padding difference between two annotation sources.
#[derive(Debug, Clone, Copy)]
pub struct RefDiff<'a> {
    pub ndref: i32,
    pub tail: &'a [u8],
}

impl<'a> RefDiff<'a> {
    pub const fn equal() -> Self {
        Self {
            ndref: 0,
            tail: &[],
        }
    }
}

/// Case-insensitive ASCII byte compare.
#[inline]
fn eq_ascii_lower(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

/// Walk two slices side-by-side, returning the length of the common
/// case-insensitive prefix.
#[inline]
fn common_prefix_len(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut i = 0;
    while i < n && eq_ascii_lower(a[i], b[i]) {
        i += 1;
    }
    i
}

/// Returns the REF-padding diff between two annotation sources, or `None`
/// if the refs are incompatible.
#[inline]
pub fn diff_refs<'a>(ref1: &'a [u8], ref2: &'a [u8]) -> Option<RefDiff<'a>> {
    let n_match = common_prefix_len(ref1, ref2);
    let a_rest = &ref1[n_match..];
    let b_rest = &ref2[n_match..];
    match (a_rest.is_empty(), b_rest.is_empty()) {
        (true, true) => Some(RefDiff::equal()),
        (false, false) => None,
        (false, true) => Some(RefDiff {
            ndref: a_rest.len() as i32,
            tail: a_rest,
        }),
        (true, false) => Some(RefDiff {
            ndref: -(b_rest.len() as i32),
            tail: b_rest,
        }),
    }
}

/// Tests whether two alleles describe the same biological variant given the
/// REF-padding diff established by [`diff_refs`].
#[inline]
pub fn matches_allele(al1: &[u8], al2: &[u8], diff: &RefDiff<'_>) -> bool {
    let n_match = common_prefix_len(al1, al2);
    let a_rest = &al1[n_match..];
    let b_rest = &al2[n_match..];

    if !a_rest.is_empty() && !b_rest.is_empty() {
        return false;
    }

    if diff.ndref == 0 {
        return a_rest.is_empty() && b_rest.is_empty();
    }

    if !a_rest.is_empty() {
        if diff.ndref <= 0 {
            return false;
        }
        let n = diff.ndref as usize;
        a_rest.len() == n && suffix_matches(a_rest, diff.tail, n)
    } else {
        if diff.ndref >= 0 {
            return false;
        }
        let n = (-diff.ndref) as usize;
        b_rest.len() == n && suffix_matches(b_rest, diff.tail, n)
    }
}

#[inline]
fn suffix_matches(rest: &[u8], tail: &[u8], n: usize) -> bool {
    for j in 0..n {
        if !eq_ascii_lower(rest[j], tail[j]) {
            return false;
        }
    }
    true
}

/// Returns the 0-based index in `als1` of the allele matching `al2` under
/// the given REF diff, or `None` if no such allele exists.
#[inline]
pub fn find_allele(als1: &[&[u8]], al2: &[u8], diff: &RefDiff<'_>) -> Option<usize> {
    als1.iter().position(|a| matches_allele(a, al2, diff))
}

#[cfg(test)]
#[path = "../../../tests/unit/annotate_cpu_v2_vcmp.rs"]
mod tests;
