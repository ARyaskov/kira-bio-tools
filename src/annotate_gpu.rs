#![cfg(feature = "gpu")]

use anyhow::{Context, Result};
use cust::prelude::*;

use crate::annotate_index::{AniEntry, AniIndex};
use crate::chr_name_to_id;

/// GPU ANI context: holds device buffers and loaded PTX module.
pub struct GpuAni {
    _ctx: cust::context::Context, // must be kept alive
    module: Module,
    stream: Stream,
    d_g: DeviceBuffer<u32>,
    d_entries: DeviceBuffer<AniEntry>,
    m: u32,
    n_entries: usize,
    string_block: Vec<u8>,
}

impl GpuAni {
    /// Load ANI index into GPU buffers.
    pub fn load(ani: &AniIndex) -> Result<Self> {
        // 1) Initialize CUDA driver API
        let _ctx = cust::quick_init().context("Failed CUDA init")?;

        // 2) Load our compiled PTX (ani_kernel.ptx)
        let ptx =
            std::fs::read_to_string("ani_kernel.ptx").context("Failed to read ani_kernel.ptx")?;
        let module = Module::from_ptx(&ptx, &[]).context("Failed to load PTX")?;

        // 3) Create stream
        let stream =
            Stream::new(StreamFlags::NON_BLOCKING, None).context("Failed to create CUDA stream")?;

        // 4) Upload MPH array g[] to GPU
        let d_g = DeviceBuffer::from_slice(&ani.mph.g).context("Failed to upload g[]")?;

        // 5) Upload entries[] to GPU
        let d_entries =
            DeviceBuffer::from_slice(&ani.entries).context("Failed to upload ANI entries")?;

        Ok(Self {
            _ctx,
            module,
            stream,
            d_g,
            d_entries,
            m: ani.mph.m,
            n_entries: ani.entries.len(),
            string_block: ani.string_block.clone(),
        })
    }

    /// GPU batch-lookup for ANI keys
    pub fn lookup_batch(&self, keys: &[u64]) -> Result<Vec<i64>> {
        let n = keys.len();

        // 1) Upload keys
        let d_keys = DeviceBuffer::from_slice(keys).context("Failed to upload lookup keys")?;

        // 2) Output buffer
        let mut d_out =
            DeviceBuffer::<i64>::zeroed(n).context("Failed to allocate output buffer")?;

        let threads = 256;
        let blocks = ((n as u32) + threads - 1) / threads;

        // 3) Lookup kernel
        let func = self
            .module
            .get_function("ani_lookup_kernel")
            .context("Failed to get ani_lookup_kernel")?;

        unsafe {
            launch!(
                func<<<blocks, threads, 0, self.stream>>>(
                    d_keys.as_device_ptr(),
                    self.d_g.as_device_ptr(),
                    self.m,
                    self.d_entries.as_device_ptr(),
                    d_out.as_device_ptr(),
                    n as i32
                )
            )?;
        }

        // 4) Sync
        self.stream.synchronize()?;

        // 5) Download results
        let mut out = vec![0i64; n];
        d_out.copy_to(&mut out)?;

        Ok(out)
    }
}

pub fn annotate_vcf_ani_gpu(
    gpu: &GpuAni,
    ani: &AniIndex,
    input_vcf: &Path,
    output_vcf: &Path,
) -> Result<()> {
    eprintln!("[gpu] annotate_vcf_ani_gpu() is not yet implemented");
    anyhow::bail!("GPU annotation is not yet implemented");
}
