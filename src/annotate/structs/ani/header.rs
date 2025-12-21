use anyhow::{anyhow, Result};

pub const ANI_MAGIC: u64 = 0x494E4149524B4256;
pub const ANI_VERSION: u64 = 2;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct AniHeader {
    pub magic: u64,
    pub version: u64,
    pub n_entries: u64,
    pub mph_m: u64,
    pub mph_salt: u64,
    pub off_mph_g: u64,
    pub off_entries: u64,
    pub off_strings: u64,
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
}
