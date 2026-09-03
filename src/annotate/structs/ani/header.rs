use anyhow::{Result, anyhow};

pub const ANI_MAGIC: u64 = 0x494E4149524B4256;
pub const ANI_VERSION: u64 = 6;
pub const ANI_HEADER_END: &str = "##KIRA_BT_ANI_HEADER_END";
pub const ANI_STR_NONE: u32 = u32::MAX;

/// Sentinel chr_id marking an empty MPH slot.
pub const ANI_SENTINEL_CHR_ID: u32 = u32::MAX;

/// Default-state `AniEntry` used to fill near-MPH slot holes.
pub const fn ani_sentinel_entry() -> AniEntry {
    AniEntry {
        chr_id: ANI_SENTINEL_CHR_ID,
        pos: u32::MAX,
        ref_ofs: ANI_STR_NONE,
        alt_ofs: ANI_STR_NONE,
        id_ofs: ANI_STR_NONE,
        qual_ofs: ANI_STR_NONE,
        filter_ofs: ANI_STR_NONE,
        info_ofs: ANI_STR_NONE,
        info_len: 0,
        format_ofs: ANI_STR_NONE,
        samples_ofs: ANI_STR_NONE,
    }
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
pub struct AniHeaderV3 {
    pub magic: u64,
    pub version: u64,
    pub n_entries: u64,
    pub mph_m: u64,
    pub mph_salt: u64,
    pub off_mph_g: u64,
    pub off_entries: u64,
    pub off_strings: u64,
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct AniHeaderV4 {
    pub magic: u64,
    pub version: u64,
    pub n_entries: u64,
    pub mph_m: u64,
    pub mph_salt: u64,
    pub off_mph_g: u64,
    pub off_entries: u64,
    pub off_strings: u64,
    pub off_block_index: u64,
    pub n_blocks: u64,
    pub block_size: u32,
    pub _pad: u32,
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct AniHeaderV5 {
    pub magic: u64,
    pub version: u64,
    pub n_entries: u64,
    pub index_len: u64,
    pub off_index: u64,
    pub off_entries: u64,
    pub off_strings: u64,
    pub off_block_index: u64,
    pub n_blocks: u64,
    pub block_size: u32,
    pub _pad: u32,
}

/// On-disk v6 header.
#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
pub struct AniHeaderV6 {
    pub magic: u64,
    pub version: u64,
    pub n_entries: u64,
    pub index_len: u64,
    pub off_index: u64,
    pub off_entries: u64,
    pub off_strings: u64,
    pub off_block_index: u64,
    pub n_blocks: u64,
    pub block_size: u32,
    pub flags: u32,
    pub off_pos_index: u64,
    pub pos_index_len: u64,
    pub off_blob: u64,
    pub blob_len: u64,
    pub off_contigs: u64,
    pub contigs_len: u64,
    /// Offset of the `entry_keys` section (`u64[n_entries]`). Zero means absent.
    pub off_entry_keys: u64,
    /// Length of the `entry_keys` section in bytes (`n_entries * 8` when present).
    pub entry_keys_len: u64,
}

/// `flags` bits in [`AniHeaderV6`].
pub mod ani_flags {
    /// At least one variant in the index is an interval-style annotation
    /// (REF=`.`, ALT=`<end_pos>`).
    pub const HAS_INTERVALS: u32 = 1 << 0;
    /// Per-entry MPH verification keys are cached in the `entry_keys` section.
    pub const HAS_ENTRY_KEYS: u32 = 1 << 1;
}

pub struct AniHeader {
    pub magic: u64,
    pub version: u64,
    pub n_entries: u64,
    pub index_len: u64,
    pub off_index: u64,
    pub off_entries: u64,
    pub off_strings: u64,
    pub off_block_index: u64,
    pub n_blocks: u64,
    pub block_size: u32,
    pub flags: u32,
    pub off_pos_index: u64,
    pub pos_index_len: u64,
    pub off_blob: u64,
    pub blob_len: u64,
    pub off_contigs: u64,
    pub contigs_len: u64,
    pub off_entry_keys: u64,
    pub entry_keys_len: u64,
}

impl AniHeader {
    #[inline]
    pub fn has_intervals(&self) -> bool {
        (self.flags & ani_flags::HAS_INTERVALS) != 0
    }

    #[inline]
    pub fn has_entry_keys(&self) -> bool {
        (self.flags & ani_flags::HAS_ENTRY_KEYS) != 0
    }
}

impl AniHeader {
    pub fn validate(&self) -> Result<()> {
        if self.magic != ANI_MAGIC {
            return Err(anyhow!("Bad ANI magic"));
        }
        if self.version != ANI_VERSION {
            return Err(anyhow!("ANI version mismatch"));
        }
        Ok(())
    }
}

#[repr(C)]
#[derive(Clone, Copy, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub struct AniBlockEntry {
    pub raw_start: u64,
    pub raw_len: u32,
    pub data_len: u32,
    pub data_off: u64,
}

#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct AniEntryV2 {
    pub chr_id: u8,
    pub pos: u32,
    pub ref_ofs: u32,
    pub alt_ofs: u32,
    pub id_ofs: u32,
    pub qual_ofs: u32,
    pub filter_ofs: u32,
    pub info_ofs: u32,
    pub info_len: u32,
}

/// On-disk variant record.
#[repr(C)]
#[derive(Clone, Copy, Default, bytemuck::Pod, bytemuck::Zeroable)]
pub struct AniEntry {
    pub chr_id: u32,
    pub pos: u32,
    pub ref_ofs: u32,
    pub alt_ofs: u32,
    pub id_ofs: u32,
    pub qual_ofs: u32,
    pub filter_ofs: u32,
    pub info_ofs: u32,
    pub info_len: u32,
    pub format_ofs: u32,
    pub samples_ofs: u32,
}
