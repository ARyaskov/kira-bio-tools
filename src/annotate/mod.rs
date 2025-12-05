pub mod builder;
pub mod builder_v2;
pub mod cpu;
pub mod cpu_v2;
pub mod reader;
pub mod structs;

#[cfg(feature = "gpu")]
pub mod cuda;

#[cfg(feature = "opencl")]
pub mod opencl;

pub use builder::*;
pub use builder_v2::*;
pub use cpu::annotate_vcf_ani;
pub use cpu_v2::annotate_vcf_ani_v2;
pub use reader::*;
pub use structs::*;

#[cfg(feature = "gpu")]
pub use cuda::*;

#[cfg(feature = "opencl")]
pub use opencl::*;
