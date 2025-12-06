use super::BGZF_BLOCK_SIZE;
use flate2::read::DeflateDecoder;
use std::io::{self, Read};

const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];
const BGZF_EXTRA_LEN: u16 = 6;
const BGZF_SI1: u8 = 0x42; // 'B'
const BGZF_SI2: u8 = 0x43; // 'C'

#[derive(Clone)]
pub struct BgzfBlock {
    pub compressed: Vec<u8>,
    pub uncompressed: Vec<u8>,
    pub cdata_size: u16,
    pub file_offset: u64,
}

impl BgzfBlock {
    pub fn new(file_offset: u64) -> Self {
        Self {
            compressed: Vec::with_capacity(BGZF_BLOCK_SIZE),
            uncompressed: Vec::with_capacity(BGZF_BLOCK_SIZE),
            cdata_size: 0,
            file_offset,
        }
    }

    pub fn virtual_offset(&self) -> u64 {
        self.file_offset << 16
    }
}

pub struct BlockDecoder;

impl BlockDecoder {
    #[inline]
    pub fn decode(block: &mut BgzfBlock) -> io::Result<usize> {
        if block.compressed.len() < 18 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Block too small",
            ));
        }

        // Validate BGZF header
        if &block.compressed[0..4] != &BGZF_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid BGZF magic",
            ));
        }

        // Extract BSIZE from extra field
        let xlen = u16::from_le_bytes([block.compressed[10], block.compressed[11]]);
        if xlen < BGZF_EXTRA_LEN {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid XLEN"));
        }

        // Verify BGZF signature
        if block.compressed[12] != BGZF_SI1 || block.compressed[13] != BGZF_SI2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid BGZF signature",
            ));
        }

        let bsize = u16::from_le_bytes([block.compressed[16], block.compressed[17]]);
        block.cdata_size = bsize;

        // DEFLATE payload starts after header
        let header_size = 12 + xlen as usize;
        let trailer_size = 8; // CRC32 + ISIZE
        let deflate_size = (bsize as usize + 1) - header_size - trailer_size;

        if block.compressed.len() < header_size + deflate_size + trailer_size {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Incomplete block",
            ));
        }

        // Decompress DEFLATE stream
        block.uncompressed.clear();
        let deflate_data = &block.compressed[header_size..header_size + deflate_size];

        let mut decoder = DeflateDecoder::new(deflate_data);
        decoder.read_to_end(&mut block.uncompressed)?;

        Ok(block.uncompressed.len())
    }

    #[inline]
    pub fn is_eof_marker(data: &[u8]) -> bool {
        data.len() >= 28
            && &data[0..28]
                == &[
                    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42,
                    0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                    0x00, 0x00,
                ]
    }
}
