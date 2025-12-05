#[cfg(feature = "opencl")]
mod v1;

#[cfg(feature = "opencl")]
mod v2;

#[cfg(feature = "opencl")]
pub use v1::*;

#[cfg(feature = "opencl")]
pub use v2::*;
