use super::header::BcfHeaderDict;
use super::typed::*;
use anyhow::{Context, Result, bail};
use std::io::{Cursor, Read, Write};

/// Sanity cap for one record's payload (`l_shared + l_indiv`).
pub const MAX_RECORD_BYTES: u32 = 1 << 30;

#[derive(Default, Clone)]
pub struct BcfRecord {
    pub rid: i32,
    pub pos: i32,
    pub rlen: i32,
    pub qual: f32,
    pub id: String,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    pub filters: Vec<u32>,
    pub info: Vec<(u32, BcfValue)>,
    pub n_sample: u32,
    pub fmt: Vec<(u32, BcfFmtValues)>,
}

#[derive(Clone)]
pub enum BcfFmtValues {
    Ints { n_per_sample: usize, vals: Vec<i32> },
    Floats { n_per_sample: usize, vals: Vec<f32> },
    Strings { lens: Vec<usize>, data: Vec<u8> },
    Gt { ploidy: usize, vals: Vec<i32> },
}

/// `(rid, pos0, rlen)` of a record from its shared block.
#[derive(Clone, Copy, Debug)]
pub struct BcfRecordMeta {
    pub rid: i32,
    pub pos: i32,
    pub rlen: i32,
}

pub fn record_meta(shared: &[u8]) -> Option<BcfRecordMeta> {
    if shared.len() < 12 {
        return None;
    }
    let rid = i32::from_le_bytes([shared[0], shared[1], shared[2], shared[3]]);
    let pos = i32::from_le_bytes([shared[4], shared[5], shared[6], shared[7]]);
    let rlen = i32::from_le_bytes([shared[8], shared[9], shared[10], shared[11]]);
    Some(BcfRecordMeta { rid, pos, rlen })
}

/// INFO/END for a record, or `pos + max(rlen, 1)` when absent, as the 0-based
/// exclusive end used for indexing (htslib `rlen`).
pub fn record_end0(meta: &BcfRecordMeta) -> u64 {
    meta.pos as u64 + (meta.rlen.max(1)) as u64
}

pub fn encode_record<W: Write>(w: &mut W, line: &str, dict: &BcfHeaderDict) -> Result<()> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() < 8 { bail!("VCF line has <8 fields"); }
    let chrom = cols[0];
    let rid = *dict.contig_idx.get(chrom)
        .ok_or_else(|| anyhow::anyhow!("contig {chrom:?} not in header"))?;
    let pos: i32 = cols[1].parse::<i32>().context("POS parse")? - 1;
    let id = if cols[2] == "." { "" } else { cols[2] };
    let ref_a = cols[3];
    let alt_str = cols[4];
    let alts: Vec<&str> = if alt_str == "." || alt_str.is_empty() {
        Vec::new()
    } else { alt_str.split(',').collect() };
    let qual: f32 = if cols[5] == "." { f32::from_bits(FLOAT_MISSING_BITS) } else { cols[5].parse().context("QUAL parse")? };
    let filter_str = cols[6];
    let filters: Vec<u32> = if filter_str == "." || filter_str.is_empty() {
        Vec::new()
    } else {
        filter_str.split(';').filter_map(|f| dict.filter_idx.get(f).copied()).collect()
    };
    let info_str = cols[7];

    let mut info_entries: Vec<(u32, &str, Option<&str>)> = Vec::new();
    let mut info_end: Option<i64> = None;
    if info_str != "." && !info_str.is_empty() {
        for kv in info_str.split(';') {
            if let Some((k, v)) = kv.split_once('=') {
                if k == "END" {
                    info_end = v.parse::<i64>().ok();
                }
                if let Some(idx) = dict.info_idx.get(k).copied() {
                    info_entries.push((idx, k, Some(v)));
                }
            } else if let Some(idx) = dict.info_idx.get(kv).copied() {
                info_entries.push((idx, kv, None));
            }
        }
    }
    // htslib: rlen = max(len(REF), END - POS + 1) for symbolic alleles.
    let mut rlen = ref_a.len() as i32;
    if let Some(end) = info_end {
        let span = end - (pos as i64 + 1) + 1;
        if span > rlen as i64 && alts.iter().any(|a| a.starts_with('<')) {
            rlen = span.min(i32::MAX as i64) as i32;
        }
    }
    let n_allele = (alts.len() + 1) as u32;
    let n_info = info_entries.len() as u32;

    let (format_str, samples) = if cols.len() > 8 {
        (cols[8], &cols[9..])
    } else { ("", &[][..]) };
    let fmt_keys: Vec<&str> = if format_str.is_empty() { Vec::new() } else { format_str.split(':').collect() };
    let n_fmt = fmt_keys.len() as u32;
    let n_sample = samples.len() as u32;

    let mut shared: Vec<u8> = Vec::with_capacity(128);
    shared.extend_from_slice(&rid.to_le_bytes());
    shared.extend_from_slice(&pos.to_le_bytes());
    shared.extend_from_slice(&rlen.to_le_bytes());
    shared.extend_from_slice(&qual.to_le_bytes());
    // BCF spec: `n_allele<<16 | n_info`.
    let n_info_allele = (n_allele << 16) | n_info;
    shared.extend_from_slice(&n_info_allele.to_le_bytes());
    let n_fmt_sample = (n_fmt << 24) | (n_sample & 0x00FFFFFF);
    shared.extend_from_slice(&n_fmt_sample.to_le_bytes());

    write_typed_string(&mut shared, id.as_bytes())?;
    write_typed_string(&mut shared, ref_a.as_bytes())?;
    for a in &alts { write_typed_string(&mut shared, a.as_bytes())?; }
    let filter_ints: Vec<i32> = filters.iter().map(|&v| v as i32).collect();
    if filter_ints.is_empty() {
        shared.push(BT_NULL);
    } else {
        write_typed_ints(&mut shared, &filter_ints)?;
    }

    for (idx, k, v) in &info_entries {
        write_typed_int_one(&mut shared, *idx as i64)?;
        let typ = dict.info_field(*idx).map(|f| f.typ.as_str()).unwrap_or("String");
        encode_info_value(&mut shared, typ, *v, k)?;
    }

    let mut indiv: Vec<u8> = Vec::new();
    if n_fmt > 0 {
        let sample_parts: Vec<Vec<&str>> = samples.iter().map(|s| s.split(':').collect()).collect();
        for (ki, key) in fmt_keys.iter().enumerate() {
            let key_idx = match dict.format_idx.get(*key) {
                Some(i) => *i,
                None => continue,
            };
            let key_typ = dict.format_field(key_idx).map(|f| f.typ.as_str()).unwrap_or("String");
            let per_sample_vals: Vec<&str> = sample_parts
                .iter()
                .map(|p| p.get(ki).copied().unwrap_or("."))
                .collect();
            write_typed_int_one(&mut indiv, key_idx as i64)?;
            if *key == "GT" {
                encode_fmt_gt(&mut indiv, &per_sample_vals)?;
            } else {
                encode_fmt_value(&mut indiv, key_typ, &per_sample_vals)?;
            }
        }
    }

    let l_shared = shared.len() as u32;
    let l_indiv = indiv.len() as u32;
    w.write_all(&l_shared.to_le_bytes())?;
    w.write_all(&l_indiv.to_le_bytes())?;
    w.write_all(&shared)?;
    w.write_all(&indiv)?;
    Ok(())
}

fn encode_info_value<W: Write>(w: &mut W, typ: &str, v: Option<&str>, _key: &str) -> Result<()> {
    let Some(val) = v else { w.write_all(&[(1u8 << 4) | BT_INT8, 1])?; return Ok(()); };
    match typ {
        "Integer" => {
            let vals: Vec<i32> = val.split(',').map(|s| if s == "." { INT32_MISSING } else { s.parse::<i32>().unwrap_or(INT32_MISSING) }).collect();
            write_typed_ints(w, &vals)?;
        }
        "Float" => {
            let vals: Vec<f32> = val.split(',').map(|s| if s == "." { f32::from_bits(FLOAT_MISSING_BITS) } else { s.parse::<f32>().unwrap_or(f32::from_bits(FLOAT_MISSING_BITS)) }).collect();
            write_typed_floats(w, &vals)?;
        }
        "Flag" => { w.write_all(&[(1u8 << 4) | BT_INT8, 1])?; }
        _ => write_typed_string(w, val.as_bytes())?,
    }
    Ok(())
}

fn encode_fmt_gt<W: Write>(w: &mut W, samples: &[&str]) -> Result<()> {
    let mut max_ploidy = 1usize;
    for s in samples {
        let p = s.split(['/', '|']).count();
        if p > max_ploidy { max_ploidy = p; }
    }
    let mut vals: Vec<i32> = Vec::with_capacity(samples.len() * max_ploidy);
    for s in samples {
        let mut cur = String::new();
        let mut alleles: Vec<(Option<u32>, bool)> = Vec::new();
        let mut next_phased = false;
        for c in s.chars() {
            if c == '/' || c == '|' {
                alleles.push((parse_gt_allele(&cur), next_phased));
                cur.clear();
                next_phased = c == '|';
            } else {
                cur.push(c);
            }
        }
        alleles.push((parse_gt_allele(&cur), next_phased));
        for i in 0..max_ploidy {
            if let Some(&(a, p)) = alleles.get(i) {
                vals.push(encode_gt(a, p));
            } else {
                vals.push(INT32_VECTOR_END);
            }
        }
    }
    let typ = min_int_type(&vals);
    write_type_desc(w, max_ploidy, typ)?;
    for &v in &vals { write_int_as(w, v, typ)?; }
    Ok(())
}

fn parse_gt_allele(s: &str) -> Option<u32> {
    if s == "." || s.is_empty() { None } else { s.parse().ok() }
}

fn encode_fmt_value<W: Write>(w: &mut W, typ: &str, samples: &[&str]) -> Result<()> {
    match typ {
        "Integer" => {
            let mut max_n = 1usize;
            for s in samples { let n = if *s == "." { 1 } else { s.split(',').count() }; if n > max_n { max_n = n; } }
            let mut vals: Vec<i32> = Vec::with_capacity(samples.len() * max_n);
            for s in samples {
                let parts: Vec<&str> = if *s == "." { vec!["."] } else { s.split(',').collect() };
                for i in 0..max_n {
                    if i < parts.len() {
                        let v = parts[i];
                        vals.push(if v == "." { INT32_MISSING } else { v.parse().unwrap_or(INT32_MISSING) });
                    } else {
                        vals.push(INT32_VECTOR_END);
                    }
                }
            }
            let t = min_int_type(&vals);
            write_type_desc(w, max_n, t)?;
            for &v in &vals { write_int_as(w, v, t)?; }
        }
        "Float" => {
            let mut max_n = 1usize;
            for s in samples { let n = if *s == "." { 1 } else { s.split(',').count() }; if n > max_n { max_n = n; } }
            write_type_desc(w, max_n, BT_FLOAT)?;
            let missing = f32::from_bits(FLOAT_MISSING_BITS);
            let vend = f32::from_bits(FLOAT_VECTOR_END_BITS);
            for s in samples {
                let parts: Vec<&str> = if *s == "." { vec!["."] } else { s.split(',').collect() };
                for i in 0..max_n {
                    let f = if i < parts.len() {
                        let v = parts[i];
                        if v == "." { missing } else { v.parse().unwrap_or(missing) }
                    } else {
                        vend
                    };
                    w.write_all(&f.to_le_bytes())?;
                }
            }
        }
        _ => {
            let mut max_len = 1usize;
            for s in samples { if s.len() > max_len { max_len = s.len(); } }
            write_type_desc(w, max_len, BT_CHAR)?;
            for s in samples {
                let bytes = if *s == "." { b"." as &[u8] } else { s.as_bytes() };
                w.write_all(bytes)?;
                for _ in bytes.len()..max_len { w.write_all(&[0])?; }
            }
        }
    }
    Ok(())
}

/// Read the raw `(shared, indiv)` blocks of the next record, `None` at EOF.
pub fn read_record_raw<R: Read>(r: &mut R) -> Result<Option<(Vec<u8>, Vec<u8>)>> {
    let mut buf = [0u8; 8];
    match r.read_exact(&mut buf) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e.into()),
    }
    let l_shared = u32::from_le_bytes([buf[0], buf[1], buf[2], buf[3]]);
    let l_indiv = u32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]);
    if l_shared < 24 || l_shared.saturating_add(l_indiv) > MAX_RECORD_BYTES {
        bail!("corrupt BCF record: l_shared={l_shared} l_indiv={l_indiv}");
    }
    let mut shared = vec![0u8; l_shared as usize];
    r.read_exact(&mut shared).context("read shared")?;
    let mut indiv = vec![0u8; l_indiv as usize];
    r.read_exact(&mut indiv).context("read indiv")?;
    Ok(Some((shared, indiv)))
}

pub fn decode_record_to_vcf<R: Read>(r: &mut R, dict: &BcfHeaderDict) -> Result<Option<String>> {
    match read_record_raw(r)? {
        Some((shared, indiv)) => decode_blocks_to_vcf(&shared, &indiv, dict).map(Some),
        None => Ok(None),
    }
}

#[inline]
fn remaining(c: &Cursor<&[u8]>) -> usize {
    c.get_ref().len().saturating_sub(c.position() as usize)
}

pub fn decode_blocks_to_vcf(shared: &[u8], indiv: &[u8], dict: &BcfHeaderDict) -> Result<String> {
    let mut sc = Cursor::new(shared);
    let mut buf4 = [0u8; 4];
    sc.read_exact(&mut buf4)?; let rid = i32::from_le_bytes(buf4);
    sc.read_exact(&mut buf4)?; let pos = i32::from_le_bytes(buf4);
    sc.read_exact(&mut buf4)?; let _rlen = i32::from_le_bytes(buf4);
    sc.read_exact(&mut buf4)?; let qual = f32::from_le_bytes(buf4);
    sc.read_exact(&mut buf4)?; let n_info_allele = u32::from_le_bytes(buf4);
    sc.read_exact(&mut buf4)?; let n_fmt_sample = u32::from_le_bytes(buf4);
    let n_allele = (n_info_allele >> 16) & 0xFFFF;
    let n_info = n_info_allele & 0xFFFF;
    let n_fmt = (n_fmt_sample >> 24) & 0xFF;
    let n_sample = n_fmt_sample & 0x00FFFFFF;

    let id = read_typed_string(&mut sc)?;
    let ref_a = read_typed_string(&mut sc)?;
    let mut alts: Vec<String> = Vec::with_capacity((n_allele as usize).saturating_sub(1));
    for _ in 1..n_allele { alts.push(read_typed_string(&mut sc)?); }

    let lim = remaining(&sc);
    let filt_val = read_typed_limited(&mut sc, lim)?;
    let filter = match filt_val {
        BcfValue::Null => ".".to_string(),
        BcfValue::Ints(vs) => {
            let names: Vec<&str> = vs.iter().filter_map(|&v| {
                if v < 0 { return None; }
                dict.filter_field(v as u32).map(|f| f.id.as_str())
            }).collect();
            if names.is_empty() { ".".into() } else { names.join(";") }
        }
        _ => ".".into(),
    };

    let mut info_str = String::new();
    let mut first = true;
    for _ in 0..n_info {
        let key_idx = read_typed_int_one(&mut sc)? as u32;
        let (key_name, key_typ): (&str, &str) = match dict.info_field(key_idx) {
            Some(f) => (f.id.as_str(), f.typ.as_str()),
            None => ("", "String"),
        };
        let unk;
        let key_name = if key_name.is_empty() { unk = format!("UNK{}", key_idx); unk.as_str() } else { key_name };
        let lim = remaining(&sc);
        let val = read_typed_limited(&mut sc, lim)?;
        if !first { info_str.push(';'); }
        first = false;
        info_str.push_str(key_name);
        match val {
            BcfValue::Null => {}
            BcfValue::Ints(vs) => {
                if key_typ == "Flag" { continue; }
                info_str.push('=');
                let strs: Vec<String> = vs.iter().take_while(|&&v| v != INT32_VECTOR_END).map(|&v| if v == INT32_MISSING { ".".into() } else { v.to_string() }).collect();
                info_str.push_str(&strs.join(","));
            }
            BcfValue::Floats(vs) => {
                info_str.push('=');
                let strs: Vec<String> = vs.iter().take_while(|&&f| !float_is_vector_end(f)).map(|&f| if float_is_missing(f) { ".".into() } else { format_float(f) }).collect();
                info_str.push_str(&strs.join(","));
            }
            BcfValue::Str(bs) => {
                if key_typ != "Flag" {
                    info_str.push('=');
                    info_str.push_str(std::str::from_utf8(&bs).unwrap_or("."));
                }
            }
        }
    }
    if info_str.is_empty() { info_str.push('.'); }

    let chrom = dict.contig_name(rid as u32).map(|s| s.to_string()).unwrap_or_else(|| rid.to_string());
    let id_str = if id.is_empty() { ".".to_string() } else { id };
    let qual_str = if float_is_missing(qual) { ".".to_string() } else { format_float(qual) };
    let alt_str = if alts.is_empty() { ".".to_string() } else { alts.join(",") };

    let mut line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        chrom, pos + 1, id_str, ref_a, alt_str, qual_str, filter, info_str);

    if n_fmt > 0 && n_sample > 0 {
        let mut ic = Cursor::new(indiv);
        let mut fmt_keys: Vec<String> = Vec::with_capacity(n_fmt as usize);
        let mut per_key_vals: Vec<Vec<String>> = Vec::with_capacity(n_fmt as usize);
        for _ in 0..n_fmt {
            let key_idx = read_typed_int_one(&mut ic)? as u32;
            let key_name = dict.format_field(key_idx).map(|f| f.id.clone()).unwrap_or_else(|| format!("UNK{}", key_idx));
            let (n_per, typ) = read_type_desc(&mut ic)?;
            let need = n_per.saturating_mul(type_size(typ)).saturating_mul(n_sample as usize);
            if need > remaining(&ic) {
                bail!("corrupt BCF record: FORMAT {key_name} needs {need} bytes");
            }
            let mut col_vals: Vec<String> = Vec::with_capacity(n_sample as usize);
            for _ in 0..n_sample {
                let s = decode_format_cell(&mut ic, &key_name, n_per, typ)?;
                col_vals.push(s);
            }
            fmt_keys.push(key_name);
            per_key_vals.push(col_vals);
        }
        line.push('\t');
        line.push_str(&fmt_keys.join(":"));
        for si in 0..n_sample as usize {
            line.push('\t');
            let parts: Vec<&str> = per_key_vals.iter().map(|v| v[si].as_str()).collect();
            line.push_str(&parts.join(":"));
        }
    }
    Ok(line)
}

fn read_typed_string(r: &mut Cursor<&[u8]>) -> Result<String> {
    let lim = remaining(r);
    match read_typed_limited(r, lim)? {
        BcfValue::Null => Ok(String::new()),
        BcfValue::Str(b) => Ok(String::from_utf8_lossy(&b).into_owned()),
        _ => bail!("expected typed string"),
    }
}

fn decode_format_cell<R: Read>(r: &mut R, key: &str, n_per: usize, typ: u8) -> Result<String> {
    let mut buf = vec![0u8; n_per * type_size(typ)];
    r.read_exact(&mut buf)?;
    if key == "GT" {
        let mut out = String::with_capacity(n_per * 2);
        for i in 0..n_per {
            let v: i32 = match typ {
                BT_INT8 => buf[i] as i8 as i32,
                BT_INT16 => i16::from_le_bytes([buf[i*2], buf[i*2+1]]) as i32,
                BT_INT32 => i32::from_le_bytes([buf[i*4], buf[i*4+1], buf[i*4+2], buf[i*4+3]]),
                _ => return Ok(".".into()),
            };
            let is_end = (typ == BT_INT8 && v == INT8_VECTOR_END as i32)
                || (typ == BT_INT16 && v == INT16_VECTOR_END as i32)
                || (typ == BT_INT32 && v == INT32_VECTOR_END);
            if is_end { break; }
            let is_missing = (typ == BT_INT8 && v == INT8_MISSING as i32)
                || (typ == BT_INT16 && v == INT16_MISSING as i32)
                || (typ == BT_INT32 && v == INT32_MISSING);
            let (allele, phased) = if is_missing { (None, false) } else { decode_gt(v) };
            if i > 0 { out.push(if phased { '|' } else { '/' }); }
            match allele {
                Some(a) => out.push_str(&a.to_string()),
                None => out.push('.'),
            }
        }
        if out.is_empty() { out.push('.'); }
        return Ok(out);
    }
    let vs: Vec<String> = match typ {
        BT_INT8 => (0..n_per).map(|i| {
            let b = buf[i] as i8;
            if b == INT8_MISSING { ".".into() } else if b == INT8_VECTOR_END { String::new() } else { b.to_string() }
        }).take_while(|s| !s.is_empty()).collect(),
        BT_INT16 => (0..n_per).map(|i| {
            let s = i16::from_le_bytes([buf[i*2], buf[i*2+1]]);
            if s == INT16_MISSING { ".".into() } else if s == INT16_VECTOR_END { String::new() } else { s.to_string() }
        }).take_while(|s| !s.is_empty()).collect(),
        BT_INT32 => (0..n_per).map(|i| {
            let v = i32::from_le_bytes([buf[i*4], buf[i*4+1], buf[i*4+2], buf[i*4+3]]);
            if v == INT32_MISSING { ".".into() } else if v == INT32_VECTOR_END { String::new() } else { v.to_string() }
        }).take_while(|s| !s.is_empty()).collect(),
        BT_FLOAT => (0..n_per).map(|i| {
            let f = f32::from_le_bytes([buf[i*4], buf[i*4+1], buf[i*4+2], buf[i*4+3]]);
            if float_is_missing(f) { ".".into() } else if float_is_vector_end(f) { String::new() } else { format_float(f) }
        }).take_while(|s| !s.is_empty()).collect(),
        BT_CHAR => {
            let end = buf.iter().position(|&b| b == 0).unwrap_or(buf.len());
            return Ok(String::from_utf8_lossy(&buf[..end]).into_owned());
        }
        _ => vec![],
    };
    Ok(if vs.is_empty() { ".".into() } else { vs.join(",") })
}

/// htslib prints floats with `%g`; integers-valued floats print without a
/// fraction.
pub fn format_float(f: f32) -> String {
    if f.fract() == 0.0 && f.abs() < 1e10 { format!("{}", f as i64) } else { format!("{}", f) }
}

#[cfg(test)]
#[path = "../../tests/unit/bcf_record.rs"]
mod tests;
