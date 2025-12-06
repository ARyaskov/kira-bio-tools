mod reader;
pub mod structs;
mod utils;
mod writer;

pub use reader::{BgzfLineReader, BgzfReader};
pub use structs::{BgzfError, Result, VirtualPosition};
pub use utils::is_bgzf;
pub use writer::BgzfWriter;
