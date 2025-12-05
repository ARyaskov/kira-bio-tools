pub mod builder;
pub mod cpu;
pub mod structs;

#[cfg(feature = "gpu")]
pub mod cuda;

#[cfg(feature = "opencl")]
pub mod opencl;

pub use builder::*;
pub use cpu::annotate_vcf_ani;
pub use structs::*;

#[cfg(feature = "gpu")]
pub use cuda::*;

#[cfg(feature = "opencl")]
pub use opencl::*;
