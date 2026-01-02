use fxhash::hash64;

use crate::annotate::structs::ani::header::{AniEntry, ANI_STR_NONE};
use crate::annotate::structs::ani::index::AniIndex;
use crate::annotate::structs::bundle::{parse_info_field, AnnotationBundle};
use crate::util::{chr_name_to_id, read_cstring};

impl AniIndex {
    pub fn lookup_exact(
        &self,
        chr: &str,
        pos: u32,
        rf: &str,
        alt: &str,
    ) -> Option<AnnotationBundle> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();

        let chr_id = chr_name_to_id(chr)? as u8;

        let mut h = (chr_id as u64) << 32 | pos as u64;
        h ^= hash64(rf.as_bytes());
        h ^= hash64(alt.as_bytes());

        if debug {
            eprintln!(
                "[LOOKUP] Searching: {}:{} {}>{} key={:016x}",
                chr, pos, rf, alt, h
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
            return None;
        }

        let rf_str = read_cstring(&self.strings, e.ref_ofs as usize);
        if rf_str != rf {
            if debug {
                eprintln!("[LOOKUP] REF mismatch: expected {}, got {}", rf, rf_str);
            }
            return None;
        }

        let alt_str = read_cstring(&self.strings, e.alt_ofs as usize);
        if alt_str != alt {
            if debug {
                eprintln!("[LOOKUP] ALT mismatch: expected {}, got {}", alt, alt_str);
            }
            return None;
        }

        if debug {
            eprintln!("[LOOKUP] Found entry at idx={}", idx);
        }

        self.build_bundle(e)
    }

    pub fn lookup_any_alt(&self, chr: &str, pos: u32, rf: &str) -> Option<AnnotationBundle> {
        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
        let chr_id = chr_name_to_id(chr)? as u8;

        let mut found: Option<&AniEntry> = None;

        for e in &self.entries {
            if e.chr_id != chr_id || e.pos != pos {
                continue;
            }

            let rf_str = read_cstring(&self.strings, e.ref_ofs as usize);
            if rf_str != rf {
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
            let alt = read_cstring(&self.strings, entry.alt_ofs as usize);
            eprintln!(
                "[LOOKUP] Using single-alt fallback: {}:{} {}>{}",
                chr, pos, rf, alt
            );
        }

        self.build_bundle(entry)
    }

    fn build_bundle(&self, e: &AniEntry) -> Option<AnnotationBundle> {
        let alt_str = read_cstring(&self.strings, e.alt_ofs as usize);
        let id_str = read_cstring(&self.strings, e.id_ofs as usize);
        let qual_str = read_cstring(&self.strings, e.qual_ofs as usize);
        let filter_str = read_cstring(&self.strings, e.filter_ofs as usize);
        let info_str = read_cstring(&self.strings, e.info_ofs as usize);
        let format_str = if e.format_ofs != ANI_STR_NONE {
            read_cstring(&self.strings, e.format_ofs as usize)
        } else {
            ""
        };
        let samples_str = if e.samples_ofs != ANI_STR_NONE {
            read_cstring(&self.strings, e.samples_ofs as usize)
        } else {
            ""
        };

        let info = parse_info_field(info_str);

        let format_opt = parse_optional(format_str);
        let samples = if format_opt.is_some() && !samples_str.is_empty() {
            samples_str.split('\t').map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        Some(AnnotationBundle {
            alt: alt_str.to_string(),
            id: parse_optional(id_str),
            qual: parse_optional(qual_str),
            filter: parse_optional(filter_str),
            info,
            format_str: format_opt,
            format_samples: samples,
        })
    }
}

fn parse_optional(s: &str) -> Option<String> {
    if s.is_empty() || s == "." {
        None
    } else {
        Some(s.to_string())
    }
}
