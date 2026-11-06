// The simulator is plain Rust; building it in tests keeps the kernel's
// reference formula checked against the CPU index even without CUDA.
#[cfg(any(feature = "gpu", test))]
pub mod gpu_sim;
#[cfg(feature = "gpu")]
mod lookup;

#[cfg(feature = "gpu")]
pub use lookup::*;
