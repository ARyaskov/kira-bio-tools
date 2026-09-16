//! VCF Zarr export (`kira-bt zarr`): the columnar layout sgkit and
//! `vcf2zarr` read (VCF Zarr specification 0.2) as a Zarr v2 directory
//! store. Every field becomes one array chunked along variants (and
//! samples), so a downstream analysis touches only the columns and the
//! genomic slab it needs instead of decoding whole records: positions,
//! alleles, filters, `call_genotype` with its phase and missingness masks,
//! and the Number=1 numeric INFO/FORMAT fields (or the ones selected).
//! Chunks are zlib-compressed (numcodecs `zlib`), strings fixed-width `S`.

use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::vcf::UnifiedVcfReader;
use crate::vcf::header::{FieldNumber, FieldType, HeaderInfo};

#[derive(Clone, Debug)]
pub struct ZarrOpts {
    pub chunk_variants: usize,
    pub chunk_samples: usize,
    /// zlib level 0-9.
    pub level: u32,
    /// `INFO/X` and `FORMAT/Y` names to export besides GT; `None` selects
    /// every Number=1 Integer/Float/Flag field.
    pub fields: Option<Vec<String>>,
}

impl Default for ZarrOpts {
    fn default() -> Self {
        Self {
            chunk_variants: 10_000,
            chunk_samples: 1_000,
            level: 6,
            fields: None,
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct ZarrSummary {
    pub variants: usize,
    pub samples: usize,
    pub arrays: Vec<String>,
}

/// Element type of an exported field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Int,
    Float,
    Bool,
}

#[derive(Clone, Debug)]
struct Field {
    id: String,
    kind: Kind,
    /// `variant_<id>` (INFO) or `call_<id>` (FORMAT).
    per_sample: bool,
}

/// Shape and dtype of one array; the store writes its chunks and metadata.
struct ArrayMeta {
    name: String,
    shape: Vec<usize>,
    chunks: Vec<usize>,
    dtype: String,
    fill: String,
    dims: Vec<&'static str>,
}

fn json_str(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn json_list<T: std::fmt::Display>(v: &[T]) -> String {
    let parts: Vec<String> = v.iter().map(|x| x.to_string()).collect();
    format!("[{}]", parts.join(", "))
}

struct Store {
    root: PathBuf,
    level: u32,
}

impl Store {
    fn create(root: &Path) -> Result<Self> {
        std::fs::create_dir_all(root).with_context(|| format!("create {}", root.display()))?;
        std::fs::write(root.join(".zgroup"), "{\n    \"zarr_format\": 2\n}")?;
        Ok(Self {
            root: root.to_path_buf(),
            level: 6,
        })
    }

    fn write_root_attrs(&self, header: &str) -> Result<()> {
        let attrs = format!(
            "{{\n    \"source\": {},\n    \"vcf_header\": {},\n    \"vcf_zarr_version\": \"0.2\"\n}}",
            json_str(&format!("kira-bt {}", env!("CARGO_PKG_VERSION"))),
            json_str(header)
        );
        std::fs::write(self.root.join(".zattrs"), attrs)?;
        Ok(())
    }

    fn array_dir(&self, name: &str) -> Result<PathBuf> {
        let d = self.root.join(name);
        std::fs::create_dir_all(&d)?;
        Ok(d)
    }

    fn write_meta(&self, m: &ArrayMeta) -> Result<()> {
        let d = self.array_dir(&m.name)?;
        let zarray = format!(
            "{{\n    \"chunks\": {},\n    \"compressor\": {{\n        \"id\": \"zlib\",\n        \"level\": {}\n    }},\n    \"dimension_separator\": \".\",\n    \"dtype\": \"{}\",\n    \"fill_value\": {},\n    \"filters\": null,\n    \"order\": \"C\",\n    \"shape\": {},\n    \"zarr_format\": 2\n}}",
            json_list(&m.chunks),
            self.level,
            m.dtype,
            m.fill,
            json_list(&m.shape)
        );
        std::fs::write(d.join(".zarray"), zarray)?;
        let dims: Vec<String> = m.dims.iter().map(|d| json_str(d)).collect();
        std::fs::write(
            d.join(".zattrs"),
            format!("{{\n    \"_ARRAY_DIMENSIONS\": [{}]\n}}", dims.join(", ")),
        )?;
        Ok(())
    }

    fn write_chunk(&self, name: &str, key: &[usize], raw: &[u8]) -> Result<()> {
        let d = self.array_dir(name)?;
        let key: Vec<String> = key.iter().map(|k| k.to_string()).collect();
        let f = std::fs::File::create(d.join(key.join(".")))?;
        let mut enc =
            flate2::write::ZlibEncoder::new(f, flate2::Compression::new(self.level.min(9)));
        enc.write_all(raw)?;
        enc.finish()?;
        Ok(())
    }
}

fn s_dtype(width: usize) -> String {
    format!("|S{}", width.max(1))
}

fn fixed_strings(items: &[String], width: usize) -> Vec<u8> {
    let w = width.max(1);
    let mut out = vec![0u8; items.len() * w];
    for (i, s) in items.iter().enumerate() {
        let b = s.as_bytes();
        let n = b.len().min(w);
        out[i * w..i * w + n].copy_from_slice(&b[..n]);
    }
    out
}

/// One string array with one dimension.
fn write_string_array(
    store: &Store,
    name: &str,
    dim: &'static str,
    items: &[String],
) -> Result<()> {
    let width = items.iter().map(|s| s.len()).max().unwrap_or(0).max(1);
    let m = ArrayMeta {
        name: name.to_string(),
        shape: vec![items.len()],
        chunks: vec![items.len().max(1)],
        dtype: s_dtype(width),
        fill: "\"\"".to_string(),
        dims: vec![dim],
    };
    store.write_meta(&m)?;
    // A chunk always holds a full chunk's worth of elements, even when the
    // array is empty.
    let mut raw = fixed_strings(items, width);
    raw.resize(items.len().max(1) * width, 0);
    store.write_chunk(name, &[0], &raw)?;
    Ok(())
}

/// The fields to export: the requested `INFO/X` / `FORMAT/Y` names, or
/// every Number=1 Integer/Float/Flag field of the header.
fn select_fields(h: &HeaderInfo, wanted: Option<&[String]>) -> Result<Vec<Field>> {
    let kind_of = |ty: FieldType, number: FieldNumber, per_sample: bool| -> Option<Kind> {
        if number != FieldNumber::Fixed(1) && !(ty == FieldType::Flag && !per_sample) {
            return None;
        }
        match ty {
            FieldType::Integer => Some(Kind::Int),
            FieldType::Float => Some(Kind::Float),
            FieldType::Flag if !per_sample => Some(Kind::Bool),
            _ => None,
        }
    };
    let mut out = Vec::new();
    match wanted {
        Some(list) => {
            for w in list {
                let (per_sample, id) = if let Some(id) = w.strip_prefix("INFO/") {
                    (false, id)
                } else if let Some(id) = w.strip_prefix("FORMAT/") {
                    (true, id)
                } else {
                    bail!("--fields: expected INFO/TAG or FORMAT/TAG, got {w:?}");
                };
                if id == "GT" {
                    continue;
                }
                let map = if per_sample { &h.format } else { &h.info };
                let Some(def) = map.get(id) else {
                    bail!("--fields: {w} is not defined in the header");
                };
                let Some(kind) = kind_of(def.ty, def.number, per_sample) else {
                    bail!("--fields: {w} is not a Number=1 Integer/Float (or INFO Flag) field");
                };
                out.push(Field {
                    id: id.to_string(),
                    kind,
                    per_sample,
                });
            }
        }
        None => {
            let mut info: Vec<&String> = h.info.keys().collect();
            info.sort();
            for id in info {
                let def = &h.info[id];
                if let Some(kind) = kind_of(def.ty, def.number, false) {
                    out.push(Field {
                        id: id.clone(),
                        kind,
                        per_sample: false,
                    });
                }
            }
            let mut fmt: Vec<&String> = h.format.keys().collect();
            fmt.sort();
            for id in fmt {
                if id == "GT" {
                    continue;
                }
                let def = &h.format[id];
                if let Some(kind) = kind_of(def.ty, def.number, true) {
                    out.push(Field {
                        id: id.clone(),
                        kind,
                        per_sample: true,
                    });
                }
            }
        }
    }
    Ok(out)
}

/// What the first pass learns about the file.
struct Scan {
    n_variants: usize,
    max_alleles: usize,
    allele_width: usize,
    id_width: usize,
    ploidy: usize,
    contigs: Vec<String>,
    contig_lengths: Vec<i64>,
    filters: Vec<String>,
}

fn scan(path: &Path, h: &HeaderInfo) -> Result<Scan> {
    let mut reader = UnifiedVcfReader::open(path)?;
    let _ = reader.header()?;
    let mut contigs: Vec<String> = h.contigs.names().to_vec();
    let mut contig_lengths: Vec<i64> = (0..contigs.len() as u32)
        .map(|i| h.contigs.length(i).map(|l| l as i64).unwrap_or(-1))
        .collect();
    let mut contig_idx: BTreeMap<String, usize> = contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), i))
        .collect();
    let mut filters: Vec<String> = vec!["PASS".to_string()];
    for f in &h.filters {
        if !filters.contains(f) {
            filters.push(f.clone());
        }
    }
    let mut s = Scan {
        n_variants: 0,
        max_alleles: 1,
        allele_width: 1,
        id_width: 1,
        ploidy: 1,
        contigs: Vec::new(),
        contig_lengths: Vec::new(),
        filters: Vec::new(),
    };
    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' {
            continue;
        }
        let mut it = line.split('\t');
        let chrom = it.next().unwrap_or("");
        if !contig_idx.contains_key(chrom) {
            contig_idx.insert(chrom.to_string(), contigs.len());
            contigs.push(chrom.to_string());
            contig_lengths.push(-1);
        }
        let _pos = it.next();
        let id = it.next().unwrap_or(".");
        s.id_width = s.id_width.max(id.len());
        let refa = it.next().unwrap_or("");
        let alt = it.next().unwrap_or(".");
        s.allele_width = s.allele_width.max(refa.len());
        let mut n_alleles = 1;
        if alt != "." && !alt.is_empty() {
            for a in alt.split(',') {
                s.allele_width = s.allele_width.max(a.len());
                n_alleles += 1;
            }
        }
        s.max_alleles = s.max_alleles.max(n_alleles);
        let _qual = it.next();
        let filter = it.next().unwrap_or(".");
        if filter != "." && filter != "PASS" {
            for f in filter.split(';') {
                if !filters.iter().any(|x| x == f) {
                    filters.push(f.to_string());
                }
            }
        }
        let _info = it.next();
        if let Some(format) = it.next()
            && let Some(gi) = format.split(':').position(|k| k == "GT")
        {
            for sample in it {
                if let Some(gt) = sample.split(':').nth(gi) {
                    let p = gt.split(['/', '|']).count();
                    s.ploidy = s.ploidy.max(p);
                }
            }
        }
        s.n_variants += 1;
    }
    s.contigs = contigs;
    s.contig_lengths = contig_lengths;
    s.filters = filters;
    Ok(s)
}

/// Column buffers of one variant chunk.
struct Chunk {
    contig: Vec<i32>,
    position: Vec<i32>,
    id: Vec<String>,
    id_mask: Vec<u8>,
    allele: Vec<String>,
    quality: Vec<f32>,
    filter: Vec<u8>,
    gt: Vec<i8>,
    phased: Vec<u8>,
    gt_mask: Vec<u8>,
    info_int: Vec<Vec<i32>>,
    info_float: Vec<Vec<f32>>,
    info_bool: Vec<Vec<u8>>,
    call_int: Vec<Vec<i32>>,
    call_float: Vec<Vec<f32>>,
    n: usize,
}

impl Chunk {
    fn new(fields: &[Field]) -> Self {
        let count = |per_sample: bool, kind: Kind| {
            fields
                .iter()
                .filter(|f| f.per_sample == per_sample && f.kind == kind)
                .count()
        };
        Self {
            contig: Vec::new(),
            position: Vec::new(),
            id: Vec::new(),
            id_mask: Vec::new(),
            allele: Vec::new(),
            quality: Vec::new(),
            filter: Vec::new(),
            gt: Vec::new(),
            phased: Vec::new(),
            gt_mask: Vec::new(),
            info_int: vec![Vec::new(); count(false, Kind::Int)],
            info_float: vec![Vec::new(); count(false, Kind::Float)],
            info_bool: vec![Vec::new(); count(false, Kind::Bool)],
            call_int: vec![Vec::new(); count(true, Kind::Int)],
            call_float: vec![Vec::new(); count(true, Kind::Float)],
            n: 0,
        }
    }

    fn clear(&mut self) {
        self.contig.clear();
        self.position.clear();
        self.id.clear();
        self.id_mask.clear();
        self.allele.clear();
        self.quality.clear();
        self.filter.clear();
        self.gt.clear();
        self.phased.clear();
        self.gt_mask.clear();
        for v in self.info_int.iter_mut().chain(self.call_int.iter_mut()) {
            v.clear();
        }
        for v in self.info_float.iter_mut().chain(self.call_float.iter_mut()) {
            v.clear();
        }
        for v in self.info_bool.iter_mut() {
            v.clear();
        }
        self.n = 0;
    }
}

fn le_bytes_i32(v: &[i32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

fn le_bytes_f32(v: &[f32]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

fn le_bytes_i64(v: &[i64]) -> Vec<u8> {
    v.iter().flat_map(|x| x.to_le_bytes()).collect()
}

/// Pad a row-major buffer of `rows` rows × `row_len` elements to `total_rows` rows.
fn padded<T: Copy>(v: &[T], row_len: usize, total_rows: usize, fill: T) -> Vec<T> {
    let mut out = Vec::with_capacity(row_len * total_rows);
    out.extend_from_slice(v);
    out.resize(row_len * total_rows, fill);
    out
}

/// Split a `[chunk_variants, samples, inner]` buffer into sample chunks.
fn sample_blocks<T: Copy>(
    v: &[T],
    n_samples: usize,
    inner: usize,
    chunk_samples: usize,
    fill: T,
) -> Vec<Vec<T>> {
    let rows = v.len() / (n_samples * inner).max(1);
    let n_blocks = n_samples.div_ceil(chunk_samples.max(1)).max(1);
    let mut out = Vec::with_capacity(n_blocks);
    for b in 0..n_blocks {
        let s0 = b * chunk_samples;
        let s1 = (s0 + chunk_samples).min(n_samples);
        let mut blk = Vec::with_capacity(rows * chunk_samples * inner);
        for r in 0..rows {
            let base = r * n_samples * inner;
            blk.extend_from_slice(&v[base + s0 * inner..base + s1 * inner]);
            // Pad the sample dimension of a short last block.
            for _ in (s1 - s0)..chunk_samples {
                blk.extend(std::iter::repeat_n(fill, inner));
            }
        }
        out.push(blk);
    }
    out
}

/// Export `input` to the Zarr store at `out` (created; must not exist
/// unless `force`).
pub fn export_vcf_zarr(
    input: &Path,
    out: &Path,
    opts: &ZarrOpts,
    force: bool,
) -> Result<ZarrSummary> {
    if out.exists() {
        if !force {
            bail!("{}: exists (use --force to overwrite)", out.display());
        }
        std::fs::remove_dir_all(out)?;
    }
    let reader =
        UnifiedVcfReader::open(input).with_context(|| format!("open {}", input.display()))?;
    let headers = reader.header()?;
    let h = HeaderInfo::parse(&headers);
    let fields = select_fields(&h, opts.fields.as_deref())?;
    let sc = scan(input, &h)?;
    let n_samples = h.samples.len();
    let cv = opts.chunk_variants.max(1);
    let cs = if n_samples == 0 {
        1
    } else {
        opts.chunk_samples.max(1).min(n_samples)
    };
    let ploidy = sc.ploidy.max(1);

    let mut store = Store::create(out)?;
    store.level = opts.level.min(9);
    store.write_root_attrs(&headers.join("\n"))?;
    let mut summary = ZarrSummary {
        variants: sc.n_variants,
        samples: n_samples,
        arrays: Vec::new(),
    };

    // Dimension tables.
    write_string_array(&store, "contig_id", "contigs", &sc.contigs)?;
    store.write_meta(&ArrayMeta {
        name: "contig_length".into(),
        shape: vec![sc.contigs.len()],
        chunks: vec![sc.contigs.len().max(1)],
        dtype: "<i8".into(),
        fill: "-1".into(),
        dims: vec!["contigs"],
    })?;
    let mut lengths = sc.contig_lengths.clone();
    lengths.resize(lengths.len().max(1), -1);
    store.write_chunk("contig_length", &[0], &le_bytes_i64(&lengths))?;
    write_string_array(&store, "filter_id", "filters", &sc.filters)?;
    write_string_array(&store, "sample_id", "samples", &h.samples)?;
    summary
        .arrays
        .extend(["contig_id", "contig_length", "filter_id", "sample_id"].map(String::from));

    // Variant arrays.
    let nv = sc.n_variants;
    let meta = |name: &str,
                shape: Vec<usize>,
                chunks: Vec<usize>,
                dtype: &str,
                fill: &str,
                dims: Vec<&'static str>| ArrayMeta {
        name: name.to_string(),
        shape,
        chunks,
        dtype: dtype.to_string(),
        fill: fill.to_string(),
        dims,
    };
    let mut metas = vec![
        meta(
            "variant_contig",
            vec![nv],
            vec![cv],
            "<i4",
            "-1",
            vec!["variants"],
        ),
        meta(
            "variant_position",
            vec![nv],
            vec![cv],
            "<i4",
            "-1",
            vec!["variants"],
        ),
        meta(
            "variant_id",
            vec![nv],
            vec![cv],
            &s_dtype(sc.id_width),
            "\"\"",
            vec!["variants"],
        ),
        meta(
            "variant_id_mask",
            vec![nv],
            vec![cv],
            "|b1",
            "false",
            vec!["variants"],
        ),
        meta(
            "variant_allele",
            vec![nv, sc.max_alleles],
            vec![cv, sc.max_alleles],
            &s_dtype(sc.allele_width),
            "\"\"",
            vec!["variants", "alleles"],
        ),
        meta(
            "variant_quality",
            vec![nv],
            vec![cv],
            "<f4",
            "\"NaN\"",
            vec!["variants"],
        ),
        meta(
            "variant_filter",
            vec![nv, sc.filters.len()],
            vec![cv, sc.filters.len().max(1)],
            "|b1",
            "false",
            vec!["variants", "filters"],
        ),
    ];
    let has_gt = n_samples > 0 && h.format.contains_key("GT");
    if has_gt {
        metas.push(meta(
            "call_genotype",
            vec![nv, n_samples, ploidy],
            vec![cv, cs, ploidy],
            "|i1",
            "-1",
            vec!["variants", "samples", "ploidy"],
        ));
        metas.push(meta(
            "call_genotype_phased",
            vec![nv, n_samples],
            vec![cv, cs],
            "|b1",
            "false",
            vec!["variants", "samples"],
        ));
        metas.push(meta(
            "call_genotype_mask",
            vec![nv, n_samples, ploidy],
            vec![cv, cs, ploidy],
            "|b1",
            "false",
            vec!["variants", "samples", "ploidy"],
        ));
    }
    for f in &fields {
        let (dtype, fill) = match f.kind {
            Kind::Int => ("<i4", "-1"),
            Kind::Float => ("<f4", "\"NaN\""),
            Kind::Bool => ("|b1", "false"),
        };
        if f.per_sample {
            metas.push(meta(
                &format!("call_{}", f.id),
                vec![nv, n_samples],
                vec![cv, cs],
                dtype,
                fill,
                vec!["variants", "samples"],
            ));
        } else {
            metas.push(meta(
                &format!("variant_{}", f.id),
                vec![nv],
                vec![cv],
                dtype,
                fill,
                vec!["variants"],
            ));
        }
    }
    for m in &metas {
        store.write_meta(m)?;
        summary.arrays.push(m.name.clone());
    }

    // Second pass: fill and flush chunks.
    let contig_idx: BTreeMap<&str, i32> = sc
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_str(), i as i32))
        .collect();
    let filter_idx: BTreeMap<&str, usize> = sc
        .filters
        .iter()
        .enumerate()
        .map(|(i, f)| (f.as_str(), i))
        .collect();
    let info_fields: Vec<&Field> = fields.iter().filter(|f| !f.per_sample).collect();
    let call_fields: Vec<&Field> = fields.iter().filter(|f| f.per_sample).collect();
    let mut chunk = Chunk::new(&fields);
    let mut chunk_no = 0usize;
    let mut reader = UnifiedVcfReader::open(input)?;
    let _ = reader.header()?;
    let flush = |chunk: &mut Chunk, chunk_no: usize| -> Result<()> {
        if chunk.n == 0 {
            return Ok(());
        }
        let key = [chunk_no];
        store.write_chunk(
            "variant_contig",
            &key,
            &le_bytes_i32(&padded(&chunk.contig, 1, cv, -1)),
        )?;
        store.write_chunk(
            "variant_position",
            &key,
            &le_bytes_i32(&padded(&chunk.position, 1, cv, -1)),
        )?;
        let mut ids = chunk.id.clone();
        ids.resize(cv, String::new());
        store.write_chunk("variant_id", &key, &fixed_strings(&ids, sc.id_width))?;
        store.write_chunk("variant_id_mask", &key, &padded(&chunk.id_mask, 1, cv, 0))?;
        let mut alleles = chunk.allele.clone();
        alleles.resize(cv * sc.max_alleles, String::new());
        store.write_chunk(
            "variant_allele",
            &[chunk_no, 0],
            &fixed_strings(&alleles, sc.allele_width),
        )?;
        store.write_chunk(
            "variant_quality",
            &key,
            &le_bytes_f32(&padded(&chunk.quality, 1, cv, f32::NAN)),
        )?;
        store.write_chunk(
            "variant_filter",
            &[chunk_no, 0],
            &padded(&chunk.filter, sc.filters.len().max(1), cv, 0),
        )?;
        if has_gt {
            let gt = padded(&chunk.gt, n_samples * ploidy, cv, -1i8);
            for (b, blk) in sample_blocks(&gt, n_samples, ploidy, cs, -1i8)
                .into_iter()
                .enumerate()
            {
                let bytes: Vec<u8> = blk.iter().map(|&x| x as u8).collect();
                store.write_chunk("call_genotype", &[chunk_no, b, 0], &bytes)?;
            }
            let ph = padded(&chunk.phased, n_samples, cv, 0u8);
            for (b, blk) in sample_blocks(&ph, n_samples, 1, cs, 0u8)
                .into_iter()
                .enumerate()
            {
                store.write_chunk("call_genotype_phased", &[chunk_no, b], &blk)?;
            }
            let mask = padded(&chunk.gt_mask, n_samples * ploidy, cv, 1u8);
            for (b, blk) in sample_blocks(&mask, n_samples, ploidy, cs, 1u8)
                .into_iter()
                .enumerate()
            {
                store.write_chunk("call_genotype_mask", &[chunk_no, b, 0], &blk)?;
            }
        }
        let (mut ii, mut fi, mut bi) = (0, 0, 0);
        for f in &info_fields {
            let name = format!("variant_{}", f.id);
            match f.kind {
                Kind::Int => {
                    store.write_chunk(
                        &name,
                        &key,
                        &le_bytes_i32(&padded(&chunk.info_int[ii], 1, cv, -1)),
                    )?;
                    ii += 1;
                }
                Kind::Float => {
                    store.write_chunk(
                        &name,
                        &key,
                        &le_bytes_f32(&padded(&chunk.info_float[fi], 1, cv, f32::NAN)),
                    )?;
                    fi += 1;
                }
                Kind::Bool => {
                    store.write_chunk(&name, &key, &padded(&chunk.info_bool[bi], 1, cv, 0))?;
                    bi += 1;
                }
            }
        }
        let (mut ci, mut cf) = (0, 0);
        for f in &call_fields {
            let name = format!("call_{}", f.id);
            match f.kind {
                Kind::Int => {
                    let v = padded(&chunk.call_int[ci], n_samples, cv, -1);
                    for (b, blk) in sample_blocks(&v, n_samples, 1, cs, -1)
                        .into_iter()
                        .enumerate()
                    {
                        store.write_chunk(&name, &[chunk_no, b], &le_bytes_i32(&blk))?;
                    }
                    ci += 1;
                }
                _ => {
                    let v = padded(&chunk.call_float[cf], n_samples, cv, f32::NAN);
                    for (b, blk) in sample_blocks(&v, n_samples, 1, cs, f32::NAN)
                        .into_iter()
                        .enumerate()
                    {
                        store.write_chunk(&name, &[chunk_no, b], &le_bytes_f32(&blk))?;
                    }
                    cf += 1;
                }
            }
        }
        chunk.clear();
        Ok(())
    };

    while let Some(line) = reader.read_line()? {
        if line.is_empty() || line.as_bytes()[0] == b'#' {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            continue;
        }
        chunk
            .contig
            .push(contig_idx.get(cols[0]).copied().unwrap_or(-1));
        chunk.position.push(cols[1].parse().unwrap_or(-1));
        chunk.id.push(if cols[2] == "." {
            String::new()
        } else {
            cols[2].to_string()
        });
        chunk.id_mask.push(u8::from(cols[2] == "."));
        chunk.allele.push(cols[3].to_string());
        let mut n_alt = 0;
        if cols[4] != "." && !cols[4].is_empty() {
            for a in cols[4].split(',') {
                chunk.allele.push(a.to_string());
                n_alt += 1;
            }
        }
        for _ in n_alt + 1..sc.max_alleles {
            chunk.allele.push(String::new());
        }
        chunk.quality.push(if cols[5] == "." {
            f32::NAN
        } else {
            cols[5].parse().unwrap_or(f32::NAN)
        });
        let mut fl = vec![0u8; sc.filters.len().max(1)];
        if cols[6] != "." {
            for f in cols[6].split(';') {
                if let Some(&i) = filter_idx.get(f) {
                    fl[i] = 1;
                }
            }
        }
        chunk.filter.extend_from_slice(&fl);
        // INFO.
        let info: Vec<(&str, Option<&str>)> = if cols[7] == "." {
            Vec::new()
        } else {
            cols[7]
                .split(';')
                .map(|kv| match kv.split_once('=') {
                    Some((k, v)) => (k, Some(v)),
                    None => (kv, None),
                })
                .collect()
        };
        let (mut ii, mut fi, mut bi) = (0, 0, 0);
        for f in &info_fields {
            let val = info.iter().find(|(k, _)| *k == f.id);
            match f.kind {
                Kind::Int => {
                    chunk.info_int[ii].push(
                        val.and_then(|(_, v)| *v)
                            .and_then(|v| v.parse().ok())
                            .unwrap_or(-1),
                    );
                    ii += 1;
                }
                Kind::Float => {
                    chunk.info_float[fi].push(
                        val.and_then(|(_, v)| *v)
                            .and_then(|v| v.parse().ok())
                            .unwrap_or(f32::NAN),
                    );
                    fi += 1;
                }
                Kind::Bool => {
                    chunk.info_bool[bi].push(u8::from(val.is_some()));
                    bi += 1;
                }
            }
        }
        // Samples.
        if n_samples > 0 {
            let format: Vec<&str> = cols
                .get(8)
                .map(|f| f.split(':').collect())
                .unwrap_or_default();
            let gi = format.iter().position(|k| *k == "GT");
            let call_pos: Vec<Option<usize>> = call_fields
                .iter()
                .map(|f| format.iter().position(|k| *k == f.id))
                .collect();
            for s in 0..n_samples {
                let sample = cols.get(9 + s).copied().unwrap_or("");
                let parts: Vec<&str> = sample.split(':').collect();
                if has_gt {
                    let gt = gi.and_then(|i| parts.get(i).copied()).unwrap_or(".");
                    chunk.phased.push(u8::from(gt.contains('|')));
                    let mut k = 0;
                    for a in gt.split(['/', '|']) {
                        if k >= ploidy {
                            break;
                        }
                        let v: i8 = if a == "." || a.is_empty() {
                            -1
                        } else {
                            a.parse().unwrap_or(-1)
                        };
                        chunk.gt.push(v);
                        chunk.gt_mask.push(u8::from(v < 0));
                        k += 1;
                    }
                    for _ in k..ploidy {
                        chunk.gt.push(-2);
                        chunk.gt_mask.push(1);
                    }
                }
                let (mut ci, mut cf) = (0, 0);
                for (f, pos) in call_fields.iter().zip(&call_pos) {
                    let raw = pos.and_then(|i| parts.get(i).copied()).unwrap_or(".");
                    match f.kind {
                        Kind::Int => {
                            chunk.call_int[ci].push(if raw == "." {
                                -1
                            } else {
                                raw.parse().unwrap_or(-1)
                            });
                            ci += 1;
                        }
                        _ => {
                            chunk.call_float[cf].push(if raw == "." {
                                f32::NAN
                            } else {
                                raw.parse().unwrap_or(f32::NAN)
                            });
                            cf += 1;
                        }
                    }
                }
            }
        }
        chunk.n += 1;
        if chunk.n == cv {
            flush(&mut chunk, chunk_no)?;
            chunk_no += 1;
        }
    }
    flush(&mut chunk, chunk_no)?;
    Ok(summary)
}

#[cfg(test)]
#[path = "../tests/unit/zarr.rs"]
mod tests;
