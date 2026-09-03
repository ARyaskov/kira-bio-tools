use std::io::{Error, ErrorKind, Read, Result, Write};

pub const BT_NULL: u8 = 0;
pub const BT_INT8: u8 = 1;
pub const BT_INT16: u8 = 2;
pub const BT_INT32: u8 = 3;
pub const BT_FLOAT: u8 = 5;
pub const BT_CHAR: u8 = 7;

pub const INT8_MISSING: i8 = i8::MIN;
pub const INT8_VECTOR_END: i8 = i8::MIN + 1;
pub const INT16_MISSING: i16 = i16::MIN;
pub const INT16_VECTOR_END: i16 = i16::MIN + 1;
pub const INT32_MISSING: i32 = i32::MIN;
pub const INT32_VECTOR_END: i32 = i32::MIN + 1;
pub const FLOAT_MISSING_BITS: u32 = 0x7F800001;
pub const FLOAT_VECTOR_END_BITS: u32 = 0x7F800002;

/// Largest typed vector accepted from a file (bytes). Guards against corrupt
/// length fields allocating gigabytes.
pub const MAX_TYPED_BYTES: usize = 1 << 28;

#[inline]
pub fn float_is_missing(f: f32) -> bool { f.to_bits() == FLOAT_MISSING_BITS }

#[inline]
pub fn float_is_vector_end(f: f32) -> bool { f.to_bits() == FLOAT_VECTOR_END_BITS }

pub fn write_type_desc<W: Write>(w: &mut W, n: usize, typ: u8) -> Result<()> {
    if n < 15 {
        w.write_all(&[((n as u8) << 4) | typ])?;
    } else {
        w.write_all(&[(15u8 << 4) | typ])?;
        write_typed_int_one(w, n as i64)?;
    }
    Ok(())
}

pub fn write_typed_int_one<W: Write>(w: &mut W, x: i64) -> Result<()> {
    if x >= i8::MIN as i64 + 2 && x <= i8::MAX as i64 {
        w.write_all(&[(1u8 << 4) | BT_INT8])?;
        w.write_all(&[(x as i8) as u8])?;
    } else if x >= i16::MIN as i64 + 2 && x <= i16::MAX as i64 {
        w.write_all(&[(1u8 << 4) | BT_INT16])?;
        w.write_all(&(x as i16).to_le_bytes())?;
    } else if x >= i32::MIN as i64 + 2 && x <= i32::MAX as i64 {
        w.write_all(&[(1u8 << 4) | BT_INT32])?;
        w.write_all(&(x as i32).to_le_bytes())?;
    } else {
        return Err(Error::new(ErrorKind::InvalidData, "int out of range"));
    }
    Ok(())
}

pub fn min_int_type(vals: &[i32]) -> u8 {
    let mut t = BT_INT8;
    for &v in vals {
        if v == INT32_VECTOR_END || v == INT32_MISSING { continue; }
        if v < i8::MIN as i32 + 2 || v > i8::MAX as i32 {
            if v < i16::MIN as i32 + 2 || v > i16::MAX as i32 { return BT_INT32; }
            if t < BT_INT16 { t = BT_INT16; }
        }
    }
    t
}

/// Write one integer of width `typ`, translating the missing / vector-end sentinels.
#[inline]
pub fn write_int_as<W: Write>(w: &mut W, v: i32, typ: u8) -> Result<()> {
    match typ {
        BT_INT8 => {
            let b = if v == INT32_MISSING { INT8_MISSING as u8 }
                else if v == INT32_VECTOR_END { INT8_VECTOR_END as u8 }
                else { (v as i8) as u8 };
            w.write_all(&[b])
        }
        BT_INT16 => {
            let s: i16 = if v == INT32_MISSING { INT16_MISSING }
                else if v == INT32_VECTOR_END { INT16_VECTOR_END }
                else { v as i16 };
            w.write_all(&s.to_le_bytes())
        }
        _ => w.write_all(&v.to_le_bytes()),
    }
}

pub fn write_typed_ints<W: Write>(w: &mut W, vals: &[i32]) -> Result<()> {
    if vals.is_empty() { w.write_all(&[BT_NULL])?; return Ok(()); }
    let typ = min_int_type(vals);
    write_type_desc(w, vals.len(), typ)?;
    for &v in vals { write_int_as(w, v, typ)?; }
    Ok(())
}

pub fn write_typed_floats<W: Write>(w: &mut W, vals: &[f32]) -> Result<()> {
    if vals.is_empty() { w.write_all(&[BT_NULL])?; return Ok(()); }
    write_type_desc(w, vals.len(), BT_FLOAT)?;
    for v in vals { w.write_all(&v.to_le_bytes())?; }
    Ok(())
}

/// Typed character string; an empty string is a zero-length CHAR vector (htslib's missing ID).
pub fn write_typed_string<W: Write>(w: &mut W, s: &[u8]) -> Result<()> {
    write_type_desc(w, s.len(), BT_CHAR)?;
    w.write_all(s)?;
    Ok(())
}

pub fn read_type_desc<R: Read>(r: &mut R) -> Result<(usize, u8)> {
    let mut b = [0u8; 1];
    r.read_exact(&mut b)?;
    let typ = b[0] & 0x0F;
    let mut n = (b[0] >> 4) as usize;
    if n == 15 {
        let v = read_typed_int_one(r)?;
        if v < 0 { return Err(Error::new(ErrorKind::InvalidData, "negative typed length")); }
        n = v as usize;
    }
    Ok((n, typ))
}

pub fn read_typed_int_one<R: Read>(r: &mut R) -> Result<i64> {
    let (n, typ) = read_type_desc(r)?;
    if n != 1 { return Err(Error::new(ErrorKind::InvalidData, "expected typed int size 1")); }
    let mut buf4 = [0u8; 4];
    Ok(match typ {
        BT_INT8 => { let mut b=[0u8;1]; r.read_exact(&mut b)?; b[0] as i8 as i64 }
        BT_INT16 => { let mut b=[0u8;2]; r.read_exact(&mut b)?; i16::from_le_bytes(b) as i64 }
        BT_INT32 => { r.read_exact(&mut buf4)?; i32::from_le_bytes(buf4) as i64 }
        _ => return Err(Error::new(ErrorKind::InvalidData, "non-int typed value")),
    })
}

#[derive(Clone)]
pub enum BcfValue { Null, Ints(Vec<i32>), Floats(Vec<f32>), Str(Vec<u8>) }

pub fn type_size(t: u8) -> usize {
    match t { BT_INT8 | BT_CHAR => 1, BT_INT16 => 2, BT_INT32 | BT_FLOAT => 4, _ => 1 }
}

/// Read a typed value, refusing lengths beyond `limit` bytes (the remaining
/// payload of the record) so a corrupt length field cannot trigger a huge
/// allocation.
pub fn read_typed_limited<R: Read>(r: &mut R, limit: usize) -> Result<BcfValue> {
    let (n, typ) = read_type_desc(r)?;
    let bytes = n.saturating_mul(type_size(typ));
    if bytes > limit.min(MAX_TYPED_BYTES) {
        return Err(Error::new(ErrorKind::InvalidData, format!("typed vector of {bytes} bytes exceeds record payload")));
    }
    match typ {
        BT_NULL => Ok(BcfValue::Null),
        BT_INT8 => {
            let mut buf = vec![0u8; n];
            r.read_exact(&mut buf)?;
            Ok(BcfValue::Ints(buf.iter().map(|&b| {
                let s = b as i8;
                if s == INT8_MISSING { INT32_MISSING } else if s == INT8_VECTOR_END { INT32_VECTOR_END } else { s as i32 }
            }).collect()))
        }
        BT_INT16 => {
            let mut buf = vec![0u8; n * 2];
            r.read_exact(&mut buf)?;
            let vals = (0..n).map(|i| {
                let s = i16::from_le_bytes([buf[i*2], buf[i*2+1]]);
                if s == INT16_MISSING { INT32_MISSING } else if s == INT16_VECTOR_END { INT32_VECTOR_END } else { s as i32 }
            }).collect();
            Ok(BcfValue::Ints(vals))
        }
        BT_INT32 => {
            let mut buf = vec![0u8; n * 4];
            r.read_exact(&mut buf)?;
            let vals = (0..n).map(|i| i32::from_le_bytes([buf[i*4], buf[i*4+1], buf[i*4+2], buf[i*4+3]])).collect();
            Ok(BcfValue::Ints(vals))
        }
        BT_FLOAT => {
            let mut buf = vec![0u8; n * 4];
            r.read_exact(&mut buf)?;
            let vals = (0..n).map(|i| f32::from_le_bytes([buf[i*4], buf[i*4+1], buf[i*4+2], buf[i*4+3]])).collect();
            Ok(BcfValue::Floats(vals))
        }
        BT_CHAR => {
            let mut buf = vec![0u8; n];
            r.read_exact(&mut buf)?;
            Ok(BcfValue::Str(buf))
        }
        _ => Err(Error::new(ErrorKind::InvalidData, format!("unknown BCF type {typ}"))),
    }
}

pub fn read_typed<R: Read>(r: &mut R) -> Result<BcfValue> {
    read_typed_limited(r, MAX_TYPED_BYTES)
}

pub fn write_typed_string_with_terminator<W: Write>(w: &mut W, s: &[u8]) -> Result<()> {
    write_typed_string(w, s)
}

/// GT allele encoding: `(allele + 1) << 1 | phased`, missing allele = -1.
pub fn encode_gt(allele_idx: Option<u32>, phased: bool) -> i32 {
    let a = allele_idx.map(|a| a as i32).unwrap_or(-1);
    ((a + 1) << 1) | (phased as i32)
}

pub fn decode_gt(v: i32) -> (Option<u32>, bool) {
    if v == INT32_MISSING || v == INT32_VECTOR_END { return (None, false); }
    let phased = (v & 1) != 0;
    let a = (v >> 1) - 1;
    if a < 0 { (None, phased) } else { (Some(a as u32), phased) }
}

#[cfg(test)]
#[path = "../../tests/unit/bcf_typed.rs"]
mod tests;
