use clap::Parser;

#[derive(Parser)]
pub struct GtcheckArgs {
    #[arg(allow_hyphen_values = true, trailing_var_arg = true)]
    pub bcftools_args: Vec<String>,
}
