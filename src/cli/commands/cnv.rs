use anyhow::Result;

use crate::cli::args::CnvArgs;

pub fn cmd_cnv(args: CnvArgs) -> Result<()> {
    crate::cnv::run_from_args(&args.bcftools_args)
}
