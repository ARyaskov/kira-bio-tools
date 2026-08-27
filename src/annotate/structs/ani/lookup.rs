use crate::annotate::structs::ani::header::{ANI_STR_NONE, AniEntry};
use crate::annotate::structs::ani::index::{AniIndex, CStrRef};
use crate::annotate::structs::ani::key::make_variant_key;
use crate::annotate::structs::bundle::FieldNumber;
use crate::annotate::structs::bundle::{AnnotationBundle, parse_info_field};
use crate::util::chr_name_to_id;
use std::collections::HashMap;
use std::time::Instant;

pub struct BundleTiming {
    pub read_s: f64,
    pub info_s: f64,
    pub optional_s: f64,
    pub samples_s: f64,
}

impl AniIndex {
    pub fn build_bundle_from_entry(&self, e: &AniEntry) -> AnnotationBundle {
        self.build_bundle_from_entry_opts(e, true, true)
    }

    pub fn build_bundle_from_entry_timed(&self, e: &AniEntry) -> (AnnotationBundle, BundleTiming) {
        self.build_bundle_from_entry_timed_opts(e, true, true)
    }

    pub fn build_bundle_from_entry_opts(
        &self,
        e: &AniEntry,
        need_info: bool,
        need_format: bool,
    ) -> AnnotationBundle {
        let ref_str = self.read_cstring(e.ref_ofs as usize);
        let alt_str = self.read_cstring(e.alt_ofs as usize);
        let id_str = self.read_cstring(e.id_ofs as usize);
        let qual_str = self.read_cstring(e.qual_ofs as usize);
        let filter_str = self.read_cstring(e.filter_ofs as usize);
        let info = if need_info {
            let info_str = self.read_cstring(e.info_ofs as usize);
            parse_info_field(info_str.as_ref())
        } else {
            Vec::new()
        };

        let (format_opt, samples) = if need_format && e.format_ofs != ANI_STR_NONE {
            let format_str = self.read_cstring(e.format_ofs as usize);
            let samples_str = if e.samples_ofs != ANI_STR_NONE {
                self.read_cstring(e.samples_ofs as usize)
            } else {
                CStrRef::empty()
            };
            let format_opt = parse_optional(format_str.as_ref());
            let samples = if format_opt.is_some() && !samples_str.as_ref().is_empty() {
                let s = samples_str.as_ref();
                let mut count = 1usize;
                for b in s.as_bytes() {
                    if *b == b'\t' {
                        count += 1;
                    }
                }
                let mut out = Vec::with_capacity(count);
                for v in s.split('\t') {
                    out.push(v.to_string());
                }
                out
            } else {
                Vec::new()
            };
            (format_opt, samples)
        } else {
            (None, Vec::new())
        };

        AnnotationBundle {
            alt: alt_str.to_string(),
            id: parse_optional(id_str.as_ref()),
            qual: parse_optional(qual_str.as_ref()),
            filter: parse_optional(filter_str.as_ref()),
            info,
            format_str: format_opt,
            format_samples: samples,
            db_ref: ref_str.to_string(),
        }
    }

    pub fn build_bundle_from_entry_timed_opts(
        &self,
        e: &AniEntry,
        need_info: bool,
        need_format: bool,
    ) -> (AnnotationBundle, BundleTiming) {
        let read_start = Instant::now();
        let ref_str = self.read_cstring(e.ref_ofs as usize);
        let alt_str = self.read_cstring(e.alt_ofs as usize);
        let id_str = self.read_cstring(e.id_ofs as usize);
        let qual_str = self.read_cstring(e.qual_ofs as usize);
        let filter_str = self.read_cstring(e.filter_ofs as usize);
        let info_str = if need_info {
            self.read_cstring(e.info_ofs as usize)
        } else {
            CStrRef::empty()
        };
        let (format_str, samples_str) = if need_format && e.format_ofs != ANI_STR_NONE {
            let format_str = self.read_cstring(e.format_ofs as usize);
            let samples_str = if e.samples_ofs != ANI_STR_NONE {
                self.read_cstring(e.samples_ofs as usize)
            } else {
                CStrRef::empty()
            };
            (format_str, samples_str)
        } else {
            (CStrRef::empty(), CStrRef::empty())
        };
        let read_s = read_start.elapsed().as_secs_f64();

        let info_start = Instant::now();
        let info = if need_info {
            parse_info_field(info_str.as_ref())
        } else {
            Vec::new()
        };
        let info_s = info_start.elapsed().as_secs_f64();

        let opt_start = Instant::now();
        let format_opt = if need_format {
            parse_optional(format_str.as_ref())
        } else {
            None
        };
        let id_opt = parse_optional(id_str.as_ref());
        let qual_opt = parse_optional(qual_str.as_ref());
        let filter_opt = parse_optional(filter_str.as_ref());
        let optional_s = opt_start.elapsed().as_secs_f64();

        let samples_start = Instant::now();
        let samples = if format_opt.is_some() && !samples_str.as_ref().is_empty() {
            let s = samples_str.as_ref();
            let mut count = 1usize;
            for b in s.as_bytes() {
                if *b == b'\t' {
                    count += 1;
                }
            }
            let mut out = Vec::with_capacity(count);
            for v in s.split('\t') {
                out.push(v.to_string());
            }
            out
        } else {
            Vec::new()
        };
        let samples_s = samples_start.elapsed().as_secs_f64();

        (
            AnnotationBundle {
                alt: alt_str.to_string(),
                id: id_opt,
                qual: qual_opt,
                filter: filter_opt,
                info,
                format_str: format_opt,
                format_samples: samples,
                db_ref: ref_str.to_string(),
            },
            BundleTiming {
                read_s,
                info_s,
                optional_s,
                samples_s,
            },
        )
    }

    pub fn lookup_exact(
        &self,
        chr: &str,
        pos: u32,
        rf: &str,
        alt: &str,
    ) -> Option<AnnotationBundle> {
        // Header-derived contig dict first, fall back to the static name table
        // for legacy `.ani` files without an embedded dict (contig_count==0).
        let chr_id = self
            .contig_id(chr)
            .or_else(|| chr_name_to_id(chr).map(u32::from))?;
        // `rf_hash` was a micro-opt for the legacy XOR key (compute once,
        // reuse across multiple ALTs). The new non-commutative `make_variant_key`
        // hashes REF inline; the parameter is kept for API stability with the
        // existing GPU/OpenCL call sites and is ignored.
        self.lookup_exact_by_chr_id(chr_id, pos, rf, 0, alt)
    }

    pub fn lookup_exact_by_chr_id(
        &self,
        chr_id: u32,
        pos: u32,
        rf: &str,
        _rf_hash: u64,
        alt: &str,
    ) -> Option<AnnotationBundle> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

        let h = make_variant_key(chr_id, pos, rf.as_bytes(), alt.as_bytes());

        if debug {
            eprintln!(
                "[LOOKUP] Searching: {}:{} {}>{} key={:016x}",
                chr_id, pos, rf, alt, h
            );
        }

        // kira_kv_engine 0.6.0: lookup_u64_fast skips the BackendDispatch enum
        // match for PtrHash25-backed indexes (~5-10 ns / lookup faster than
        // lookup_u64). Returns None only for non-PtrHash25 engines, which 0.6.0
        // doesn't have yet — fall back to lookup_u64 to stay future-proof.
        let idx = match self
            .index
            .lookup_u64_fast(h)
            .unwrap_or_else(|| self.index.lookup_u64(h))
        {
            Ok(v) => v,
            Err(_) => {
                if debug {
                    eprintln!("[LOOKUP] Index miss for key={:016x}", h);
                }
                if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                    return Some(self.build_bundle_from_entry(&self.entries[idx]));
                }
                return None;
            }
        };

        if idx >= self.entries.len() {
            if debug {
                eprintln!(
                    "[LOOKUP] Index returned idx {} >= entries.len() {}",
                    idx,
                    self.entries.len()
                );
            }
            return None;
        }

        let e = &self.entries[idx];

        if e.chr_id != chr_id || e.pos != pos {
            if debug {
                eprintln!(
                    "[LOOKUP] Position mismatch: entry chr={} pos={}",
                    e.chr_id, e.pos
                );
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                return Some(self.build_bundle_from_entry(&self.entries[idx]));
            }
            return None;
        }

        let rf_str = self.read_cstring(e.ref_ofs as usize);
        if rf_str.as_ref() != rf {
            if debug {
                eprintln!("[LOOKUP] REF mismatch: expected {}, got {}", rf, rf_str);
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                return Some(self.build_bundle_from_entry(&self.entries[idx]));
            }
            return None;
        }

        let alt_str = self.read_cstring(e.alt_ofs as usize);
        if alt_str.as_ref() != alt {
            if debug {
                eprintln!("[LOOKUP] ALT mismatch: expected {}, got {}", alt, alt_str);
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                return Some(self.build_bundle_from_entry(&self.entries[idx]));
            }
            return None;
        }

        if debug {
            eprintln!("[LOOKUP] Found entry at idx={}", idx);
        }

        Some(self.build_bundle_from_entry(e))
    }

    pub fn lookup_exact_by_chr_id_opts(
        &self,
        chr_id: u32,
        pos: u32,
        rf: &str,
        _rf_hash: u64,
        alt: &str,
        need_info: bool,
        need_format: bool,
    ) -> Option<AnnotationBundle> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

        let h = make_variant_key(chr_id, pos, rf.as_bytes(), alt.as_bytes());

        if debug {
            eprintln!(
                "[LOOKUP] Searching: {}:{} {}>{} key={:016x}",
                chr_id, pos, rf, alt, h
            );
        }

        // kira_kv_engine 0.6.0: lookup_u64_fast skips the BackendDispatch enum
        // match for PtrHash25-backed indexes (~5-10 ns / lookup faster than
        // lookup_u64). Returns None only for non-PtrHash25 engines, which 0.6.0
        // doesn't have yet — fall back to lookup_u64 to stay future-proof.
        let idx = match self
            .index
            .lookup_u64_fast(h)
            .unwrap_or_else(|| self.index.lookup_u64(h))
        {
            Ok(v) => v,
            Err(_) => {
                if debug {
                    eprintln!("[LOOKUP] Index miss for key={:016x}", h);
                }
                if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                    return Some(self.build_bundle_from_entry_opts(
                        &self.entries[idx],
                        need_info,
                        need_format,
                    ));
                }
                return None;
            }
        };

        if idx >= self.entries.len() {
            if debug {
                eprintln!(
                    "[LOOKUP] Index returned idx {} >= entries.len() {}",
                    idx,
                    self.entries.len()
                );
            }
            return None;
        }

        let e = &self.entries[idx];

        if e.chr_id != chr_id || e.pos != pos {
            if debug {
                eprintln!(
                    "[LOOKUP] Position mismatch: entry chr={} pos={}",
                    e.chr_id, e.pos
                );
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                return Some(self.build_bundle_from_entry_opts(
                    &self.entries[idx],
                    need_info,
                    need_format,
                ));
            }
            return None;
        }

        let rf_str = self.read_cstring(e.ref_ofs as usize);
        if rf_str.as_ref() != rf {
            if debug {
                eprintln!("[LOOKUP] REF mismatch: expected {}, got {}", rf, rf_str);
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                return Some(self.build_bundle_from_entry_opts(
                    &self.entries[idx],
                    need_info,
                    need_format,
                ));
            }
            return None;
        }

        let alt_str = self.read_cstring(e.alt_ofs as usize);
        if alt_str.as_ref() != alt {
            if debug {
                eprintln!("[LOOKUP] ALT mismatch: expected {}, got {}", alt, alt_str);
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                return Some(self.build_bundle_from_entry_opts(
                    &self.entries[idx],
                    need_info,
                    need_format,
                ));
            }
            return None;
        }

        if debug {
            eprintln!("[LOOKUP] Found entry at idx={}", idx);
        }

        Some(self.build_bundle_from_entry_opts(e, need_info, need_format))
    }

    pub fn lookup_exact_by_chr_id_pos_index_opts(
        &self,
        chr_id: u32,
        pos: u32,
        rf: &str,
        rf_hash: u64,
        alt: &str,
        field_meta: &HashMap<String, FieldNumber>,
        need_info: bool,
        need_format: bool,
    ) -> Option<AnnotationBundle> {
        self.lookup_exact_by_chr_id_pos_index_opts_filtered(
            chr_id, pos, rf, rf_hash, alt, field_meta, need_info, need_format, None,
        )
    }

    /// Selective-info variant. When `info_filter` is `Some`, the bundle's
    /// `info` field only carries the matching tags — saves ~85% of per-record
    /// `String` allocations on selective annotation runs
    /// (`-c INFO/A,INFO/B,INFO/C` vs full INFO).
    #[allow(clippy::too_many_arguments)]
    pub fn lookup_exact_by_chr_id_pos_index_opts_filtered(
        &self,
        chr_id: u32,
        pos: u32,
        rf: &str,
        rf_hash: u64,
        alt: &str,
        field_meta: &HashMap<String, FieldNumber>,
        need_info: bool,
        need_format: bool,
        info_filter: Option<&std::collections::HashSet<&str>>,
    ) -> Option<AnnotationBundle> {
        use crate::annotate::cpu_v2::vcmp;

        if let Some(list) = self.lookup_pos_index(chr_id, pos) {
            // Pass 1: exact REF + exact ALT. The fastest and overwhelmingly
            // common case — same source/target naming conventions.
            for &idx in list {
                let e = &self.entries[idx as usize];
                if e.chr_id != chr_id || e.pos != pos {
                    continue;
                }
                let rf_str = self.read_cstring(e.ref_ofs as usize);
                if rf_str.as_ref() != rf {
                    continue;
                }
                let alt_str = self.read_cstring(e.alt_ofs as usize);
                if alt_str.as_ref() != alt {
                    continue;
                }
                return Some(self.build_bundle_from_entry_idx_opts_with_meta_filtered(
                    idx as usize,
                    field_meta,
                    need_info,
                    need_format,
                    info_filter,
                ));
            }
            // Pass 2: vcmp REF-padding-aware match. Catches the bcftools
            // case where source and target encode the same indel with
            // different REF padding (e.g. user has `REF=A ALT=AT` but the
            // database has `REF=ATC ALT=ATTC`). Pass 1's exact-match would
            // have hit any candidate where rf is byte-identical to the
            // stored REF, so anything reaching here has a REF length
            // mismatch — exactly the vcmp's territory.
            for &idx in list {
                let e = &self.entries[idx as usize];
                if e.chr_id != chr_id || e.pos != pos {
                    continue;
                }
                let rf_str = self.read_cstring(e.ref_ofs as usize);
                let Some(diff) = vcmp::diff_refs(rf.as_bytes(), rf_str.as_ref().as_bytes())
                else {
                    continue;
                };
                let alt_str = self.read_cstring(e.alt_ofs as usize);
                if !vcmp::matches_allele(
                    alt.as_bytes(),
                    alt_str.as_ref().as_bytes(),
                    &diff,
                ) {
                    continue;
                }
                if std::env::var("KIRA_BT_DEBUG").is_ok() {
                    eprintln!(
                        "[LOOKUP] vcmp hit: {chr_id}:{pos} ({rf}>{alt}) ↔ ({}>{})",
                        rf_str.as_ref(),
                        alt_str.as_ref()
                    );
                }
                return Some(self.build_bundle_from_entry_idx_opts_with_meta_filtered(
                    idx as usize,
                    field_meta,
                    need_info,
                    need_format,
                    info_filter,
                ));
            }
        }

        self.lookup_exact_by_chr_id_opts(chr_id, pos, rf, rf_hash, alt, need_info, need_format)
    }

    pub fn lookup_exact_by_chr_id_timed(
        &self,
        chr_id: u32,
        pos: u32,
        rf: &str,
        rf_hash: u64,
        alt: &str,
    ) -> Option<(AnnotationBundle, BundleTiming)> {
        self.lookup_exact_by_chr_id_timed_opts(chr_id, pos, rf, rf_hash, alt, true, true)
    }

    pub fn lookup_exact_by_chr_id_timed_opts(
        &self,
        chr_id: u32,
        pos: u32,
        rf: &str,
        _rf_hash: u64,
        alt: &str,
        need_info: bool,
        need_format: bool,
    ) -> Option<(AnnotationBundle, BundleTiming)> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

        let h = make_variant_key(chr_id, pos, rf.as_bytes(), alt.as_bytes());

        if debug {
            eprintln!(
                "[LOOKUP] Searching: {}:{} {}>{} key={:016x}",
                chr_id, pos, rf, alt, h
            );
        }

        // kira_kv_engine 0.6.0: lookup_u64_fast skips the BackendDispatch enum
        // match for PtrHash25-backed indexes (~5-10 ns / lookup faster than
        // lookup_u64). Returns None only for non-PtrHash25 engines, which 0.6.0
        // doesn't have yet — fall back to lookup_u64 to stay future-proof.
        let idx = match self
            .index
            .lookup_u64_fast(h)
            .unwrap_or_else(|| self.index.lookup_u64(h))
        {
            Ok(v) => v,
            Err(_) => {
                if debug {
                    eprintln!("[LOOKUP] Index miss for key={:016x}", h);
                }
                if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                    let (bundle, t) = self.build_bundle_from_entry_timed_opts(
                        &self.entries[idx],
                        need_info,
                        need_format,
                    );
                    return Some((bundle, t));
                }
                return None;
            }
        };

        if idx >= self.entries.len() {
            if debug {
                eprintln!(
                    "[LOOKUP] Index returned idx {} >= entries.len() {}",
                    idx,
                    self.entries.len()
                );
            }
            return None;
        }

        let e = &self.entries[idx];

        if e.chr_id != chr_id || e.pos != pos {
            if debug {
                eprintln!(
                    "[LOOKUP] Position mismatch: entry chr={} pos={}",
                    e.chr_id, e.pos
                );
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                let (bundle, t) = self.build_bundle_from_entry_timed_opts(
                    &self.entries[idx],
                    need_info,
                    need_format,
                );
                return Some((bundle, t));
            }
            return None;
        }

        let rf_str = self.read_cstring(e.ref_ofs as usize);
        if rf_str.as_ref() != rf {
            if debug {
                eprintln!("[LOOKUP] REF mismatch: expected {}, got {}", rf, rf_str);
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                let (bundle, t) = self.build_bundle_from_entry_timed_opts(
                    &self.entries[idx],
                    need_info,
                    need_format,
                );
                return Some((bundle, t));
            }
            return None;
        }

        let alt_str = self.read_cstring(e.alt_ofs as usize);
        if alt_str.as_ref() != alt {
            if debug {
                eprintln!("[LOOKUP] ALT mismatch: expected {}, got {}", alt, alt_str);
            }
            if let Some(idx) = self.find_interval_entry(chr_id, pos) {
                let (bundle, t) = self.build_bundle_from_entry_timed_opts(
                    &self.entries[idx],
                    need_info,
                    need_format,
                );
                return Some((bundle, t));
            }
            return None;
        }

        if debug {
            eprintln!("[LOOKUP] Found entry at idx={}", idx);
        }

        let (bundle, t) = self.build_bundle_from_entry_timed_opts(e, need_info, need_format);
        Some((bundle, t))
    }

    pub fn lookup_any_alt(&self, chr: &str, pos: u32, rf: &str) -> Option<AnnotationBundle> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
        let chr_id = self
            .contig_id(chr)
            .or_else(|| chr_name_to_id(chr).map(u32::from))?;

        // Use the pos-index to bound candidates to the exact position window
        // instead of scanning all 4M entries. Falls back to None (rather than
        // a linear scan) when no pos-index is present — legacy `.ani` users
        // can rebuild.
        let candidates = self.lookup_pos_index(chr_id, pos)?;

        let mut found: Option<&AniEntry> = None;
        for &entry_idx in candidates {
            let e = &self.entries[entry_idx as usize];
            if e.chr_id != chr_id || e.pos != pos {
                continue;
            }
            let rf_str = self.read_cstring(e.ref_ofs as usize);
            if rf_str.as_ref() != rf {
                continue;
            }
            if found.is_some() {
                if debug {
                    eprintln!(
                        "[LOOKUP] Multiple ALT candidates found for {}:{} {}",
                        chr, pos, rf
                    );
                }
                return None;
            }
            found = Some(e);
        }

        let entry = found?;

        if debug {
            let alt = self.read_cstring(entry.alt_ofs as usize);
            eprintln!(
                "[LOOKUP] Using single-alt fallback: {}:{} {}>{}",
                chr, pos, rf, alt
            );
        }

        Some(self.build_bundle_from_entry(entry))
    }
}

fn parse_optional(s: &str) -> Option<String> {
    if s.is_empty() || s == "." {
        None
    } else {
        Some(s.to_string())
    }
}
