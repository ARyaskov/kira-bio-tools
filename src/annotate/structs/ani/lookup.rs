use super::header::AniEntry;
use super::index::AniIndex;
use crate::annotate::structs::bundle::parse_info_field;
use crate::annotate::AnnotationBundle;
use crate::util::read_cstring;

impl AniIndex {
    pub fn lookup(&self, chr: &str, pos: u32, rf: &str, alt: &str) -> Option<AnnotationBundle> {
        let entry = self.find_entry(chr, pos, rf, alt)?;
        self.build_bundle(&entry)
    }

    pub fn lookup_any_alt(&self, chr: &str, pos: u32, rf: &str) -> Option<AnnotationBundle> {
        use crate::chr_name_to_id;

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
                return None;
            }
            found = Some(e);
        }

        let entry = found?;
        self.build_bundle(entry)
    }

    fn find_entry(&self, chr: &str, pos: u32, rf: &str, alt: &str) -> Option<&AniEntry> {
        use crate::chr_name_to_id;
        use fxhash::hash64;

        let debug = std::env::var("KIRA_BT_DEBUG").is_ok();
        let chr_id = chr_name_to_id(chr)? as u8;

        let mut h = (chr_id as u64) << 32 | (pos as u64);
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
                    "[LOOKUP] Position mismatch: entry has chr={} pos={}",
                    e.chr_id, e.pos
                );
            }
            return None;
        }

        if !self.verify_ref_alt(e, rf, alt, debug) {
            return None;
        }

        if debug {
            eprintln!("[LOOKUP] Found entry at idx={}", idx);
        }

        Some(e)
    }

    fn verify_ref_alt(&self, e: &AniEntry, rf: &str, alt: &str, debug: bool) -> bool {
        let rf_str = read_cstring(&self.strings, e.ref_ofs as usize);
        if rf_str != rf {
            if debug {
                eprintln!(
                    "[DEBUG-LOOKUP] REF mismatch: expected {}, got {}",
                    rf, rf_str
                );
            }
            return false;
        }

        let alt_str = read_cstring(&self.strings, e.alt_ofs as usize);
        if alt_str != alt {
            if debug {
                eprintln!(
                    "[DEBUG-LOOKUP] ALT mismatch: expected {}, got {}",
                    alt, alt_str
                );
            }
            return false;
        }

        true
    }

    fn build_bundle(&self, e: &AniEntry) -> Option<AnnotationBundle> {
        let alt_str = read_cstring(&self.strings, e.alt_ofs as usize);
        let id_str = read_cstring(&self.strings, e.id_ofs as usize);
        let qual_str = read_cstring(&self.strings, e.qual_ofs as usize);
        let filter_str = read_cstring(&self.strings, e.filter_ofs as usize);
        let info_str = read_cstring(&self.strings, e.info_ofs as usize);

        let info_fields = parse_info_field(info_str);

        Some(AnnotationBundle {
            alt: alt_str.to_string(),
            id: parse_optional_field(id_str),
            qual: parse_optional_field(qual_str),
            filter: parse_optional_field(filter_str),
            info: info_fields,
        })
    }
}

fn parse_optional_field(s: &str) -> Option<String> {
    if s == "." || s.is_empty() {
        None
    } else {
        Some(s.to_string())
    }
}
