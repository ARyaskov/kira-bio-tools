//! persistence.rs — flat binary PGMI v2.4 format + mmap read/write.
//!
//! File format v2.4 (section order is strictly fixed):
//!   [PgmiHeader]
//!   [DiskSegment; n_segments]                // PGM segments
//!   [u64 anchors; n_anchors]                 // (optional; can be 0)
//!   [ChrIndexEntry; n_chr]                   // chr -> {name_off, name_len, id}
//!   [u8 chr_names_blob; chr_names_len]       // concatenation of chr names
//!   [ChrIdMapEntry; n_idmap]                 // id -> {name_off, name_len} (id — 0..max_id); n_idmap = max_id+1
//!   [u8  offsets_comp; offsets_comp_len]     // varint-encoded deltas of VCF line offsets
//!   [OffsCheckpoint; n_ckpts]                // checkpoints for varint stream (every stride records)
//!
//! Offset compression:
//!   - store absolute off[0] separately in first checkpoint (abs_offset),
//!     then the stream contains varint-encoded deltas (off[i] - off[i-1]) >= 0.
//!   - Checkpoint: { index, abs_offset, blob_pos }, where blob_pos is position in offsets_comp.
//!   - Default stride = 1024 (defined by CKPT_STRIDE constant).

use std::{
    fs::{File, OpenOptions},
    io::{self, BufWriter, Read, Seek, Write},
    mem,
    path::Path,
    slice,
    sync::Arc,
};

#[cfg(feature = "mmap")]
use memmap2::{Mmap, MmapMut, MmapOptions};

pub type Anchor = u64;

pub const CKPT_STRIDE: usize = 1024;

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct DiskSegment {
    pub key_lo: u64,
    pub key_hi: u64,
    pub slope: f32,
    pub intercept: f32,
    pub base_rank: u64,
}
pub type Segment = DiskSegment;

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ChrIndexEntry {
    pub name_off: u64,
    pub name_len: u32,
    pub id: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct ChrIdMapEntry {
    pub name_off: u64,
    pub name_len: u32,
    pub _pad: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct OffsCheckpoint {
    pub index: u64,
    pub abs_offset: u64,
    pub blob_pos: u64,
}

const PGMI_MAGIC: &[u8; 8] = b"PGMIv2.4";
const PGMI_ENDIAN_TAG: u32 = 0x01020304;
const PGMI_VERSION: u32 = 0x0002_0004;

#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct PgmiHeader {
    pub magic: [u8; 8],
    pub version: u32,
    pub endian: u32,
    pub epsilon: u32,
    pub ckpt_stride: u32,

    pub n_segments: u64,
    pub off_segments: u64,

    pub n_anchors: u64,
    pub off_anchors: u64,

    pub n_chr: u64,
    pub off_chrindex: u64,

    pub chr_names_len: u64,
    pub off_chrnames: u64,

    pub n_idmap: u64,
    pub off_idmap: u64,

    pub offsets_comp_len: u64,
    pub off_offsets_comp: u64,

    pub n_ckpts: u64,
    pub off_ckpts: u64,
}

pub struct PgmiOwned {
    pub header: PgmiHeader,
    pub segments: Vec<Segment>,
    pub anchors: Vec<Anchor>,
    pub chr_index: Vec<ChrIndexEntry>,
    pub chr_names: Vec<u8>,
    pub idmap: Vec<ChrIdMapEntry>,
    pub offsets_comp: Vec<u8>,
    pub ckpts: Vec<OffsCheckpoint>,
}

#[cfg(feature = "mmap")]
pub struct PgmiMapped {
    _mmap: Arc<Mmap>,
    pub header: &'static PgmiHeader,
    pub segments: &'static [Segment],
    pub anchors: &'static [Anchor],
    pub chr_index: &'static [ChrIndexEntry],
    pub chr_names: &'static [u8],
    pub idmap: &'static [ChrIdMapEntry],
    pub offsets_comp: &'static [u8],
    pub ckpts: &'static [OffsCheckpoint],
}

#[inline]
unsafe fn as_u8_slice<T: Sized>(v: &[T]) -> &[u8] {
    slice::from_raw_parts(v.as_ptr() as *const u8, v.len() * mem::size_of::<T>())
}
#[inline]
unsafe fn as_u8_slice_mut<T: Sized>(v: &mut [T]) -> &mut [u8] {
    slice::from_raw_parts_mut(v.as_mut_ptr() as *mut u8, v.len() * mem::size_of::<T>())
}

#[inline]
fn check_header(h: &PgmiHeader) -> io::Result<()> {
    if &h.magic != PGMI_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "pgmi: bad magic (expect v2.4)",
        ));
    }
    if h.version != PGMI_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "pgmi: version mismatch (expect v2.4)",
        ));
    }
    if h.endian != PGMI_ENDIAN_TAG {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "pgmi: endian mismatch",
        ));
    }
    if h.ckpt_stride as usize != CKPT_STRIDE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "pgmi: ckpt stride mismatch",
        ));
    }
    Ok(())
}

pub fn compute_pgmi_size(
    n_segments: usize,
    n_anchors: usize,
    n_chr: usize,
    chr_names_len: usize,
    n_idmap: usize,
    offsets_comp_len: usize,
    n_ckpts: usize,
) -> u64 {
    let header = mem::size_of::<PgmiHeader>() as u64;
    let seg = (n_segments * mem::size_of::<Segment>()) as u64;
    let anch = (n_anchors * mem::size_of::<Anchor>()) as u64;
    let chridx = (n_chr * mem::size_of::<ChrIndexEntry>()) as u64;
    let chrnames = chr_names_len as u64;
    let idmap = (n_idmap * mem::size_of::<ChrIdMapEntry>()) as u64;
    let offs = offsets_comp_len as u64;
    let ckpts = (n_ckpts * mem::size_of::<OffsCheckpoint>()) as u64;
    header + seg + anch + chridx + chrnames + idmap + offs + ckpts
}

pub fn save_pgmi(
    path: &Path,
    epsilon: u32,
    segments: &[Segment],
    anchors: &[Anchor],
) -> io::Result<()> {
    save_pgmi_v24(path, epsilon, segments, anchors, &[], &[], &[], &[], &[])
}


pub fn save_pgmi_v24(
    path: &Path,
    epsilon: u32,
    segments: &[Segment],
    anchors: &[Anchor],
    chr_index: &[ChrIndexEntry],
    chr_names: &[u8],
    idmap: &[ChrIdMapEntry],
    offsets_comp: &[u8],
    ckpts: &[OffsCheckpoint],
) -> io::Result<()> {
    save_pgmi_flat(
        path,
        epsilon,
        segments,
        anchors,
        chr_index,
        chr_names,
        idmap,
        offsets_comp,
        ckpts,
    )
}

pub fn save_pgmi_flat(
    path: &Path,
    epsilon: u32,
    segments: &[Segment],
    anchors: &[Anchor],
    chr_index: &[ChrIndexEntry],
    chr_names: &[u8],
    idmap: &[ChrIdMapEntry],
    offsets_comp: &[u8],
    ckpts: &[OffsCheckpoint],
) -> io::Result<()> {
    let mut file = File::create(path)?;
    let mut bw = BufWriter::with_capacity(32 * 1024 * 1024, &mut file);

    let header_size = mem::size_of::<PgmiHeader>() as u64;
    let seg_bytes = (segments.len() * mem::size_of::<Segment>()) as u64;
    let anch_bytes = (anchors.len() * mem::size_of::<Anchor>()) as u64;
    let chridx_bytes = (chr_index.len() * mem::size_of::<ChrIndexEntry>()) as u64;
    let chrnames_bytes = chr_names.len() as u64;
    let idmap_bytes = (idmap.len() * mem::size_of::<ChrIdMapEntry>()) as u64;
    let offs_bytes = offsets_comp.len() as u64;
    let ckpt_bytes = (ckpts.len() * mem::size_of::<OffsCheckpoint>()) as u64;

    let off_segments = header_size;
    let off_anchors = off_segments + seg_bytes;
    let off_chrindex = off_anchors + anch_bytes;
    let off_chrnames = off_chrindex + chridx_bytes;
    let off_idmap = off_chrnames + chrnames_bytes;
    let off_offsets_comp = off_idmap + idmap_bytes;
    let off_ckpts = off_offsets_comp + offs_bytes;

    let header = PgmiHeader {
        magic: *PGMI_MAGIC,
        version: PGMI_VERSION,
        endian: PGMI_ENDIAN_TAG,
        epsilon,
        ckpt_stride: CKPT_STRIDE as u32,

        n_segments: segments.len() as u64,
        off_segments,

        n_anchors: anchors.len() as u64,
        off_anchors,

        n_chr: chr_index.len() as u64,
        off_chrindex,

        chr_names_len: chr_names.len() as u64,
        off_chrnames,

        n_idmap: idmap.len() as u64,
        off_idmap,

        offsets_comp_len: offsets_comp.len() as u64,
        off_offsets_comp,

        n_ckpts: ckpts.len() as u64,
        off_ckpts,
    };

    unsafe {
        bw.write_all(as_u8_slice(slice::from_ref(&header)))?;
    }
    unsafe {
        bw.write_all(as_u8_slice(segments))?;
    }
    unsafe {
        bw.write_all(as_u8_slice(anchors))?;
    }
    unsafe {
        bw.write_all(as_u8_slice(chr_index))?;
    }
    bw.write_all(chr_names)?;
    unsafe {
        bw.write_all(as_u8_slice(idmap))?;
    }
    bw.write_all(offsets_comp)?;
    unsafe {
        bw.write_all(as_u8_slice(ckpts))?;
    }
    bw.flush()?;
    Ok(())
}

#[cfg(feature = "mmap")]
pub fn save_pgmi_mmap(
    path: &Path,
    epsilon: u32,
    segments: &[Segment],
    anchors: &[Anchor],
    chr_index: &[ChrIndexEntry],
    chr_names: &[u8],
    idmap: &[ChrIdMapEntry],
    offsets_comp: &[u8],
    ckpts: &[OffsCheckpoint],
) -> io::Result<()> {
    let total = compute_pgmi_size(
        segments.len(),
        anchors.len(),
        chr_index.len(),
        chr_names.len(),
        idmap.len(),
        offsets_comp.len(),
        ckpts.len(),
    );

    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(path)?;
    file.set_len(total)?;
    let mut mmap = unsafe { MmapMut::map_mut(&file)? };

    let header_size = mem::size_of::<PgmiHeader>() as u64;
    let seg_bytes = (segments.len() * mem::size_of::<Segment>()) as u64;
    let anch_bytes = (anchors.len() * mem::size_of::<Anchor>()) as u64;
    let chridx_bytes = (chr_index.len() * mem::size_of::<ChrIndexEntry>()) as u64;
    let chrnames_bytes = chr_names.len() as u64;
    let idmap_bytes = (idmap.len() * mem::size_of::<ChrIdMapEntry>()) as u64;
    let offs_bytes = offsets_comp.len() as u64;
    let ckpt_bytes = (ckpts.len() * mem::size_of::<OffsCheckpoint>()) as u64;

    let off_segments = header_size;
    let off_anchors = off_segments + seg_bytes;
    let off_chrindex = off_anchors + anch_bytes;
    let off_chrnames = off_chrindex + chridx_bytes;
    let off_idmap = off_chrnames + chrnames_bytes;
    let off_offsets_comp = off_idmap + idmap_bytes;
    let off_ckpts = off_offsets_comp + offs_bytes;

    let header = PgmiHeader {
        magic: *PGMI_MAGIC,
        version: PGMI_VERSION,
        endian: PGMI_ENDIAN_TAG,
        epsilon,
        ckpt_stride: CKPT_STRIDE as u32,

        n_segments: segments.len() as u64,
        off_segments,

        n_anchors: anchors.len() as u64,
        off_anchors,

        n_chr: chr_index.len() as u64,
        off_chrindex,

        chr_names_len: chr_names.len() as u64,
        off_chrnames,

        n_idmap: idmap.len() as u64,
        off_idmap,

        offsets_comp_len: offsets_comp.len() as u64,
        off_offsets_comp,

        n_ckpts: ckpts.len() as u64,
        off_ckpts,
    };

    unsafe {
        // header
        mmap[..mem::size_of::<PgmiHeader>()].copy_from_slice(as_u8_slice(slice::from_ref(&header)));
        // segments
        let mut p = off_segments as usize;
        mmap[p..p + seg_bytes as usize].copy_from_slice(as_u8_slice(segments));
        p += seg_bytes as usize;
        // anchors
        mmap[p..p + anch_bytes as usize].copy_from_slice(as_u8_slice(anchors));
        p += anch_bytes as usize;
        // chrindex
        mmap[p..p + chridx_bytes as usize].copy_from_slice(as_u8_slice(chr_index));
        p += chridx_bytes as usize;
        // chrnames
        mmap[p..p + chrnames_bytes as usize].copy_from_slice(chr_names);
        p += chrnames_bytes as usize;
        // idmap
        mmap[p..p + idmap_bytes as usize].copy_from_slice(as_u8_slice(idmap));
        p += idmap_bytes as usize;
        // offsets comp blob
        mmap[p..p + offs_bytes as usize].copy_from_slice(offsets_comp);
        p += offs_bytes as usize;
        // checkpoints
        mmap[p..p + ckpt_bytes as usize].copy_from_slice(as_u8_slice(ckpts));
    }

    mmap.flush()?;
    Ok(())
}

trait SeekExt {
    fn seek_io(&mut self, pos: std::io::SeekFrom) -> io::Result<u64>;
}
impl SeekExt for File {
    fn seek_io(&mut self, pos: std::io::SeekFrom) -> io::Result<u64> {
        self.seek(pos)
    }
}

pub fn load_pgmi_owned(path: &Path) -> io::Result<PgmiOwned> {
    let mut file = File::open(path)?;
    if (file.metadata()?.len() as usize) < mem::size_of::<PgmiHeader>() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "pgmi: too small",
        ));
    }
    let mut header = PgmiHeader {
        magic: *PGMI_MAGIC,
        version: 0,
        endian: 0,
        epsilon: 0,
        ckpt_stride: CKPT_STRIDE as u32,
        n_segments: 0,
        off_segments: 0,
        n_anchors: 0,
        off_anchors: 0,
        n_chr: 0,
        off_chrindex: 0,
        chr_names_len: 0,
        off_chrnames: 0,
        n_idmap: 0,
        off_idmap: 0,
        offsets_comp_len: 0,
        off_offsets_comp: 0,
        n_ckpts: 0,
        off_ckpts: 0,
    };
    unsafe {
        let hdr = as_u8_slice_mut(slice::from_mut(&mut header));
        file.read_exact(hdr)?;
    }
    check_header(&header)?;

    let segments =
        read_vec_from_file::<Segment>(path, header.off_segments, header.n_segments as usize)?;
    let anchors =
        read_vec_from_file::<Anchor>(path, header.off_anchors, header.n_anchors as usize)?;
    let chr_index =
        read_vec_from_file::<ChrIndexEntry>(path, header.off_chrindex, header.n_chr as usize)?;
    let chr_names = read_bytes_from_file(path, header.off_chrnames, header.chr_names_len as usize)?;
    let idmap =
        read_vec_from_file::<ChrIdMapEntry>(path, header.off_idmap, header.n_idmap as usize)?;
    let offsets_comp = read_bytes_from_file(
        path,
        header.off_offsets_comp,
        header.offsets_comp_len as usize,
    )?;
    let ckpts =
        read_vec_from_file::<OffsCheckpoint>(path, header.off_ckpts, header.n_ckpts as usize)?;

    Ok(PgmiOwned {
        header,
        segments,
        anchors,
        chr_index,
        chr_names,
        idmap,
        offsets_comp,
        ckpts,
    })
}

fn read_vec_from_file<T: Sized>(path: &Path, offset: u64, count: usize) -> io::Result<Vec<T>> {
    if count == 0 {
        return Ok(Vec::new());
    }
    let mut f = File::open(path)?;
    f.seek_io(std::io::SeekFrom::Start(offset))?;
    let mut v: Vec<T> = Vec::with_capacity(count);
    unsafe {
        v.set_len(count);
    }
    let bytes = unsafe { as_u8_slice_mut(&mut v) };
    f.read_exact(bytes)?;
    Ok(v)
}
fn read_bytes_from_file(path: &Path, offset: u64, len: usize) -> io::Result<Vec<u8>> {
    if len == 0 {
        return Ok(Vec::new());
    }
    let mut f = File::open(path)?;
    f.seek_io(std::io::SeekFrom::Start(offset))?;
    let mut v = vec![0u8; len];
    f.read_exact(&mut v)?;
    Ok(v)
}

#[cfg(feature = "mmap")]
pub fn load_pgmi_mmap(path: &Path) -> io::Result<PgmiMapped> {
    let file = File::open(path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    let arc = Arc::new(mmap);
    if arc.len() < mem::size_of::<PgmiHeader>() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "pgmi: too small for header",
        ));
    }
    let header = unsafe { &*(arc.as_ptr() as *const PgmiHeader) };
    check_header(header)?;

    let base = arc.as_ptr() as usize;

    let segments: &'static [Segment] = if header.n_segments > 0 {
        let ptr = (base + header.off_segments as usize) as *const Segment;
        unsafe { slice::from_raw_parts(ptr, header.n_segments as usize) }
    } else {
        &[]
    };

    let anchors: &'static [Anchor] = if header.n_anchors > 0 {
        let ptr = (base + header.off_anchors as usize) as *const Anchor;
        unsafe { slice::from_raw_parts(ptr, header.n_anchors as usize) }
    } else {
        &[]
    };

    let chr_index: &'static [ChrIndexEntry] = if header.n_chr > 0 {
        let ptr = (base + header.off_chrindex as usize) as *const ChrIndexEntry;
        unsafe { slice::from_raw_parts(ptr, header.n_chr as usize) }
    } else {
        &[]
    };

    let chr_names: &'static [u8] = if header.chr_names_len > 0 {
        let ptr = (base + header.off_chrnames as usize) as *const u8;
        unsafe { slice::from_raw_parts(ptr, header.chr_names_len as usize) }
    } else {
        &[]
    };

    let idmap: &'static [ChrIdMapEntry] = if header.n_idmap > 0 {
        let ptr = (base + header.off_idmap as usize) as *const ChrIdMapEntry;
        unsafe { slice::from_raw_parts(ptr, header.n_idmap as usize) }
    } else {
        &[]
    };

    let offsets_comp: &'static [u8] = if header.offsets_comp_len > 0 {
        let ptr = (base + header.off_offsets_comp as usize) as *const u8;
        unsafe { slice::from_raw_parts(ptr, header.offsets_comp_len as usize) }
    } else {
        &[]
    };

    let ckpts: &'static [OffsCheckpoint] = if header.n_ckpts > 0 {
        let ptr = (base + header.off_ckpts as usize) as *const OffsCheckpoint;
        unsafe { slice::from_raw_parts(ptr, header.n_ckpts as usize) }
    } else {
        &[]
    };

    Ok(PgmiMapped {
        _mmap: arc,
        header,
        segments,
        anchors,
        chr_index,
        chr_names,
        idmap,
        offsets_comp,
        ckpts,
    })
}

pub enum PgmiIndex {
    #[cfg(feature = "mmap")]
    Mapped(PgmiMapped),
    Owned(PgmiOwned),
}

impl PgmiIndex {
    #[inline]
    pub fn segments(&self) -> &[Segment] {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => m.segments,
            PgmiIndex::Owned(o) => &o.segments,
        }
    }
    #[inline]
    pub fn anchors(&self) -> &[Anchor] {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => m.anchors,
            PgmiIndex::Owned(o) => &o.anchors,
        }
    }
    #[inline]
    pub fn chr_index(&self) -> &[ChrIndexEntry] {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => m.chr_index,
            PgmiIndex::Owned(o) => &o.chr_index,
        }
    }
    #[inline]
    pub fn chr_names(&self) -> &[u8] {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => m.chr_names,
            PgmiIndex::Owned(o) => &o.chr_names,
        }
    }
    #[inline]
    pub fn idmap(&self) -> &[ChrIdMapEntry] {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => m.idmap,
            PgmiIndex::Owned(o) => &o.idmap,
        }
    }
    #[inline]
    pub fn offsets_comp(&self) -> &[u8] {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => m.offsets_comp,
            PgmiIndex::Owned(o) => &o.offsets_comp,
        }
    }
    #[inline]
    pub fn ckpts(&self) -> &[OffsCheckpoint] {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => m.ckpts,
            PgmiIndex::Owned(o) => &o.ckpts,
        }
    }
    #[inline]
    pub fn epsilon(&self) -> u32 {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => m.header.epsilon,
            PgmiIndex::Owned(o) => o.header.epsilon,
        }
    }
    #[inline]
    pub fn header(&self) -> PgmiHeader {
        match self {
            #[cfg(feature = "mmap")]
            PgmiIndex::Mapped(m) => *m.header,
            PgmiIndex::Owned(o) => o.header,
        }
    }
}

pub fn load_pgmi(path: &Path) -> io::Result<PgmiIndex> {
    #[cfg(feature = "mmap")]
    {
        return load_pgmi_mmap(path).map(PgmiIndex::Mapped);
    }
    #[cfg(not(feature = "mmap"))]
    {
        return load_pgmi_owned(path).map(PgmiIndex::Owned);
    }
}
