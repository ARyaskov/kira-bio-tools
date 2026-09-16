//! Alignment kernels shared by BAQ, haplotype scoring and indel discovery.

pub mod glocal;

pub use glocal::{GlocalParams, GlocalResult, GlocalScratch, encode_nt, glocal, glocal_loglik, glocal_with};
