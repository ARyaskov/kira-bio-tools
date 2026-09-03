//! Variant classification. Port of htslib `bcf_set_variant_type` /
//! `bcf_get_variant_types` so every command agrees on what a SNP or indel is.

pub const VT_REF: u32 = 0;
pub const VT_SNP: u32 = 1;
pub const VT_MNP: u32 = 2;
pub const VT_INDEL: u32 = 4;
pub const VT_OTHER: u32 = 8;
pub const VT_BND: u32 = 16;
pub const VT_OVERLAP: u32 = 32;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AlleleType {
    pub ty: u32,
    /// Length change: >0 insertion, <0 deletion, 1 for SNP, 0 for REF.
    pub n: i32,
}

#[inline]
fn eq_ic(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

/// Type of one ALT allele relative to REF.
pub fn allele_type(r: &str, a: &str) -> AlleleType {
    let rb = r.as_bytes();
    let ab = a.as_bytes();

    if ab == b"*" {
        return AlleleType { ty: VT_OVERLAP, n: 0 };
    }
    if rb.len() == 1 && ab.len() == 1 {
        if ab[0] == b'.' || eq_ic(rb[0], ab[0]) || ab[0] == b'X' {
            return AlleleType { ty: VT_REF, n: 0 };
        }
        return AlleleType { ty: VT_SNP, n: 1 };
    }
    if ab.is_empty() || ab == b"." {
        return AlleleType { ty: VT_REF, n: 0 };
    }
    if ab[0] == b'<' {
        if ab == b"<X>" || ab == b"<*>" || ab == b"<NON_REF>" {
            return AlleleType { ty: VT_REF, n: 0 };
        }
        return AlleleType { ty: VT_OTHER, n: 0 };
    }
    if ab[0] == b'[' || ab[0] == b']' {
        return AlleleType { ty: VT_BND, n: 0 };
    }

    let mut p = 0usize;
    while p < rb.len() && p < ab.len() && eq_ic(rb[p], ab[p]) {
        p += 1;
    }
    let (rr, aa) = (&rb[p..], &ab[p..]);
    if !aa.is_empty() && rr.is_empty() {
        if aa[0] == b']' || aa[0] == b'[' {
            return AlleleType { ty: VT_BND, n: 0 };
        }
        return AlleleType { ty: VT_INDEL, n: aa.len() as i32 };
    }
    if !rr.is_empty() && aa.is_empty() {
        return AlleleType { ty: VT_INDEL, n: -(rr.len() as i32) };
    }
    if rr.is_empty() && aa.is_empty() {
        return AlleleType { ty: VT_REF, n: 0 };
    }
    if aa.iter().any(|&b| b == b'[' || b == b']') {
        return AlleleType { ty: VT_BND, n: 0 };
    }

    let mut re = rr.len() - 1;
    let mut ae = aa.len() - 1;
    while re > 0 && ae > 0 && eq_ic(rr[re], aa[ae]) {
        re -= 1;
        ae -= 1;
    }
    if ae == 0 {
        if re == 0 {
            return AlleleType { ty: VT_SNP, n: 1 };
        }
        let ty = if eq_ic(rr[re], aa[ae]) { VT_INDEL } else { VT_OTHER };
        return AlleleType { ty, n: -(re as i32) };
    }
    if re == 0 {
        let ty = if eq_ic(rr[re], aa[ae]) { VT_INDEL } else { VT_OTHER };
        return AlleleType { ty, n: ae as i32 };
    }
    let ty = if re == ae { VT_MNP } else { VT_OTHER };
    let n = if re > ae { -(re as i32 + 1) } else { ae as i32 + 1 };
    AlleleType { ty, n }
}

/// Union of allele types of a record (`bcf_get_variant_types`).
pub fn record_type(r: &str, alt_field: &str) -> u32 {
    if alt_field.is_empty() || alt_field == "." {
        return VT_REF;
    }
    let mut t = VT_REF;
    for a in alt_field.split(',') {
        t |= allele_type(r, a).ty;
    }
    t
}

/// Per-allele types in ALT order.
pub fn allele_types(r: &str, alt_field: &str) -> Vec<AlleleType> {
    if alt_field.is_empty() || alt_field == "." {
        return Vec::new();
    }
    alt_field.split(',').map(|a| allele_type(r, a)).collect()
}

#[inline]
pub fn has_snp(r: &str, alt: &str) -> bool {
    record_type(r, alt) & VT_SNP != 0
}

#[inline]
pub fn has_indel(r: &str, alt: &str) -> bool {
    record_type(r, alt) & VT_INDEL != 0
}

#[inline]
pub fn has_mnp(r: &str, alt: &str) -> bool {
    record_type(r, alt) & VT_MNP != 0
}

/// True when every ALT is a SNP (and there is at least one).
pub fn is_pure_snp(r: &str, alt: &str) -> bool {
    let ts = allele_types(r, alt);
    !ts.is_empty() && ts.iter().all(|t| t.ty == VT_SNP)
}

/// Name used by the filter language `TYPE` field.
pub fn type_name(ty: u32) -> &'static str {
    if ty & VT_BND != 0 {
        "bnd"
    } else if ty & VT_OTHER != 0 {
        "other"
    } else if ty & VT_INDEL != 0 {
        "indel"
    } else if ty & VT_MNP != 0 {
        "mnp"
    } else if ty & VT_SNP != 0 {
        "snp"
    } else if ty & VT_OVERLAP != 0 {
        "overlap"
    } else {
        "ref"
    }
}

/// Parse a `-v/-V` type list (`snps,indels,mnps,other,bnd,ref`) into a mask.
pub fn parse_type_mask(spec: &str) -> Option<u32> {
    let mut m = 0u32;
    for t in spec.split(',') {
        m |= match t.trim() {
            "snps" | "snp" => VT_SNP,
            "indels" | "indel" => VT_INDEL,
            "mnps" | "mnp" => VT_MNP,
            "other" => VT_OTHER,
            "bnd" | "bnds" => VT_BND,
            "ref" => VT_REF,
            "overlap" => VT_OVERLAP,
            "" => 0,
            _ => return None,
        };
    }
    Some(m)
}

#[cfg(test)]
#[path = "../../tests/unit/vcf_variant_type.rs"]
mod tests;
