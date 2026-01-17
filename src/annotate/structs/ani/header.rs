use anyhow::{anyhow, Result};

pub const ANI_MAGIC: u64 = 0x494E4149524B4256;
pub const ANI_VERSION: u64 = 5;
pub const ANI_HEADER_END: &str = "##KIRA_BT_ANI_HEADER_END";
pub const ANI_STR_NONE: u32 = u32::MAX;

#[repr(C)]
#[derive(Clone, Copy)]
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
#[derive(Clone, Copy, Default)]
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

#[repr(C)]
#[derive(Clone, Copy, Default)]
pub struct AniEntry {
    pub chr_id: u8,
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
