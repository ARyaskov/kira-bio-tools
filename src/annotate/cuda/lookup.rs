use anyhow::{Context, Result};
use cust::prelude::*;
use std::path::Path;

use crate::annotate::structs::*;

pub struct GpuAni {
    _ctx: cust::context::Context,
    module: Module,
    stream: Stream,
    d_g: DeviceBuffer<u32>,
    d_entries: DeviceBuffer<AniEntry>,
    m: u32,
    n_entries: usize,
    string_block: Vec<u8>,
}

impl GpuAni {
    pub fn load(ani: &AniIndex) -> Result<Self> {
        let _ctx = cust::quick_init().context("Failed CUDA init")?;

        let ptx =
            std::fs::read_to_string("ani_kernel.ptx").context("Failed to read ani_kernel.ptx")?;
        let module = Module::from_ptx(&ptx, &[]).context("Failed to load PTX")?;

        let stream =
            Stream::new(StreamFlags::NON_BLOCKING, None).context("Failed to create CUDA stream")?;

        let d_g = DeviceBuffer::from_slice(&ani.mph.g).context("Failed to upload g[]")?;

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
            string_block: ani.strings_owned(),
        })
    }

    pub fn lookup_batch(&self, keys: &[u64]) -> Result<Vec<i64>> {
        let n = keys.len();

        let d_keys = DeviceBuffer::from_slice(keys).context("Failed to upload lookup keys")?;

        let mut d_out =
            DeviceBuffer::<i64>::zeroed(n).context("Failed to allocate output buffer")?;

        let threads = 256;
        let blocks = ((n as u32) + threads - 1) / threads;

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

        self.stream.synchronize()?;

        let mut out = vec![0i64; n];
        d_out.copy_to(&mut out)?;

        Ok(out)
    }
}

pub fn annotate_vcf_ani_gpu(
    _gpu: &GpuAni,
    _ani: &AniIndex,
    _input_vcf: &Path,
    _output_vcf: &Path,
) -> Result<()> {
    eprintln!("[gpu] annotate_vcf_ani_gpu() is not yet implemented");
    anyhow::bail!("GPU annotation is not yet implemented");
}
