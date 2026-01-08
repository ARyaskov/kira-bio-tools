use crate::util::fast_hash64;

use crate::annotate::structs::ani::header::{AniEntry, ANI_STR_NONE};
use crate::annotate::structs::ani::index::{AniIndex, CStrRef};
use crate::annotate::structs::bundle::{parse_info_field, AnnotationBundle};
use crate::util::chr_name_to_id;
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
        }
    }

    pub fn build_bundle_from_entry_timed_opts(
        &self,
        e: &AniEntry,
        need_info: bool,
        need_format: bool,
    ) -> (AnnotationBundle, BundleTiming) {
        let read_start = Instant::now();
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
        let chr_id = chr_name_to_id(chr)? as u8;
        let rf_hash = fast_hash64(rf.as_bytes());

        self.lookup_exact_by_chr_id(chr_id, pos, rf, rf_hash, alt)
    }

    pub fn lookup_exact_by_chr_id(
        &self,
        chr_id: u8,
        pos: u32,
        rf: &str,
        rf_hash: u64,
        alt: &str,
    ) -> Option<AnnotationBundle> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

        let mut h = (chr_id as u64) << 32 | pos as u64;
        h ^= rf_hash;
        h ^= fast_hash64(alt.as_bytes());

        if debug {
            eprintln!(
                "[LOOKUP] Searching: {}:{} {}>{} key={:016x}",
                chr_id, pos, rf, alt, h
            );
        }

        let idx = self.mph.index(&h.to_le_bytes()) as usize;

        if idx >= self.entries.len() {
            if debug {
                eprintln!(
                    "[LOOKUP] MPH returned idx {} >= entries.len() {}",
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
        chr_id: u8,
        pos: u32,
        rf: &str,
        rf_hash: u64,
        alt: &str,
        need_info: bool,
        need_format: bool,
    ) -> Option<AnnotationBundle> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

        let mut h = (chr_id as u64) << 32 | pos as u64;
        h ^= rf_hash;
        h ^= fast_hash64(alt.as_bytes());

        if debug {
            eprintln!(
                "[LOOKUP] Searching: {}:{} {}>{} key={:016x}",
                chr_id, pos, rf, alt, h
            );
        }

        let idx = self.mph.index(&h.to_le_bytes()) as usize;

        if idx >= self.entries.len() {
            if debug {
                eprintln!(
                    "[LOOKUP] MPH returned idx {} >= entries.len() {}",
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

    pub fn lookup_exact_by_chr_id_timed(
        &self,
        chr_id: u8,
        pos: u32,
        rf: &str,
        rf_hash: u64,
        alt: &str,
    ) -> Option<(AnnotationBundle, BundleTiming)> {
        self.lookup_exact_by_chr_id_timed_opts(chr_id, pos, rf, rf_hash, alt, true, true)
    }

    pub fn lookup_exact_by_chr_id_timed_opts(
        &self,
        chr_id: u8,
        pos: u32,
        rf: &str,
        rf_hash: u64,
        alt: &str,
        need_info: bool,
        need_format: bool,
    ) -> Option<(AnnotationBundle, BundleTiming)> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

        let mut h = (chr_id as u64) << 32 | pos as u64;
        h ^= rf_hash;
        h ^= fast_hash64(alt.as_bytes());

        if debug {
            eprintln!(
                "[LOOKUP] Searching: {}:{} {}>{} key={:016x}",
                chr_id, pos, rf, alt, h
            );
        }

        let idx = self.mph.index(&h.to_le_bytes()) as usize;

        if idx >= self.entries.len() {
            if debug {
                eprintln!(
                    "[LOOKUP] MPH returned idx {} >= entries.len() {}",
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
        let chr_id = chr_name_to_id(chr)? as u8;

        let mut found: Option<&AniEntry> = None;

        for e in &self.entries {
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
