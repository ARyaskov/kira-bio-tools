//! `solid --call` fuses the caller into the pileup with no VCF in between, so
//! only running both paths on the same reads proves they compute the same thing.

use std::io::Write;
use std::path::{Path, PathBuf};

use clap::Parser;

use kira_bio_tools::VcfReader;
use kira_bio_tools::call::stream::{CallConfig, call_stream};
use kira_bio_tools::cli::args::SolidArgs;
use kira_bio_tools::cli::commands::solid::cmd_solid;

/// Deterministic pseudo-random DNA.
fn dna(len: usize, mut state: u64) -> Vec<u8> {
    let alphabet = b"ACGT";
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            alphabet[(state >> 60) as usize & 3]
        })
        .collect()
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

fn substitute(b: u8) -> u8 {
    match b {
        b'A' => b'G',
        b'G' => b'A',
        b'C' => b'T',
        _ => b'C',
    }
}

/// Reads sampled from a donor copy of the reference carrying substitutions at
/// fixed positions, so the pileup sees real variant sites.
fn make_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let reference = dna(50_000, 0xC0FFEE);
    let ref_path = dir.join("ref.fa");
    let mut f = std::fs::File::create(&ref_path).unwrap();
    writeln!(f, ">chrT").unwrap();
    for chunk in reference.chunks(60) {
        f.write_all(chunk).unwrap();
        f.write_all(b"\n").unwrap();
    }
    drop(f);

    let mut donor = reference.clone();
    for pos in [5_000usize, 12_500, 20_000, 31_000, 44_000] {
        donor[pos] = substitute(donor[pos]);
    }

    let r1_path = dir.join("r1.fq");
    let r2_path = dir.join("r2.fq");
    let mut f1 = std::fs::File::create(&r1_path).unwrap();
    let mut f2 = std::fs::File::create(&r2_path).unwrap();
    let qual: String = "I".repeat(150);
    let mut state = 0x1234_5678u64;
    for i in 0..2_000 {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let frag = 300 + (state >> 55) as usize % 100;
        let start = (state >> 20) as usize % (donor.len() - frag - 1);
        let r1 = donor[start..start + 150].to_vec();
        let r2: Vec<u8> = donor[start + frag - 150..start + frag]
            .iter()
            .rev()
            .map(|&b| complement(b))
            .collect();
        writeln!(f1, "@rd{i}/1\n{}\n+\n{qual}", String::from_utf8_lossy(&r1)).unwrap();
        writeln!(f2, "@rd{i}/2\n{}\n+\n{qual}", String::from_utf8_lossy(&r2)).unwrap();
    }
    drop(f1);
    drop(f2);
    (ref_path, r1_path, r2_path)
}

fn run_solid(ref_path: &Path, r1: &Path, r2: &Path, out: &Path, extra: &[&str]) {
    let mut argv: Vec<String> = vec![
        "solid".into(),
        "--aligner-ref".into(),
        ref_path.display().to_string(),
        "--aligner-r1".into(),
        r1.display().to_string(),
        "--aligner-r2".into(),
        r2.display().to_string(),
        "--aligner-rg".into(),
        "ID:t\\tSM:testsample".into(),
        "--threads".into(),
        "2".into(),
        "--output".into(),
        out.display().to_string(),
    ];
    argv.extend(extra.iter().map(|s| s.to_string()));
    let args = SolidArgs::try_parse_from(&argv).unwrap();
    cmd_solid(args).unwrap();
}

fn data_lines(vcf: &str) -> Vec<&str> {
    vcf.lines().filter(|l| !l.starts_with('#')).collect()
}

fn sample_col(line: &str) -> &str {
    line.split('\t').nth(9).expect("sample column")
}

/// Apply the caller to an already-written pileup VCF, the step-by-step way.
fn call_separately(pileup: &Path, cfg: CallConfig) -> String {
    let mut reader = VcfReader::open(pileup).unwrap();
    let headers = reader.header().unwrap();
    let mut records = Vec::new();
    while let Some(rec) = reader.next_record().unwrap() {
        records.push(rec);
    }
    let mut out: Vec<u8> = Vec::new();
    call_stream(records, &headers, cfg, &mut out).unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn fused_call_matches_the_separate_call_step() {
    let dir = tempfile::tempdir().unwrap();
    let (ref_path, r1, r2) = make_fixture(dir.path());

    // Pileup only — what the pipeline emitted before the call stage existed.
    let raw = dir.path().join("raw.vcf");
    run_solid(&ref_path, &r1, &r2, &raw, &[]);
    let raw_text = std::fs::read_to_string(&raw).unwrap();
    assert!(
        !data_lines(&raw_text).is_empty(),
        "fixture produced no pileup sites"
    );

    // Fused: the same pileup piped straight into the caller.
    let fused = dir.path().join("called.vcf");
    run_solid(&ref_path, &r1, &r2, &fused, &["--call"]);
    let fused_text = std::fs::read_to_string(&fused).unwrap();

    let expected = call_separately(
        &raw,
        CallConfig {
            variants_only: true,
            ..CallConfig::default()
        },
    );
    assert_eq!(
        fused_text, expected,
        "fused call differs from the separate call step"
    );
    assert!(!data_lines(&fused_text).is_empty(), "no sites called");
}

#[test]
fn haploid_call_stage_emits_single_allele_genotypes() {
    let dir = tempfile::tempdir().unwrap();
    let (ref_path, r1, r2) = make_fixture(dir.path());

    let diploid = dir.path().join("diploid.vcf");
    run_solid(&ref_path, &r1, &r2, &diploid, &["--call"]);
    let diploid_text = std::fs::read_to_string(&diploid).unwrap();
    let diploid_calls = data_lines(&diploid_text);
    assert!(!diploid_calls.is_empty());
    assert!(
        diploid_calls.iter().all(|l| sample_col(l).contains('/')),
        "diploid run should emit paired genotypes"
    );

    let haploid = dir.path().join("haploid.vcf");
    run_solid(
        &ref_path,
        &r1,
        &r2,
        &haploid,
        &["--call", "--call-ploidy", "1"],
    );
    let haploid_text = std::fs::read_to_string(&haploid).unwrap();
    let haploid_calls = data_lines(&haploid_text);
    assert!(!haploid_calls.is_empty());
    for line in &haploid_calls {
        let gt = sample_col(line);
        assert!(
            !gt.contains('/'),
            "haploid run emitted a diploid genotype: {gt} in {line}"
        );
    }
    // The two runs see the same reads, so a ploidy change must not silently
    // change which sites are called.
    assert_eq!(haploid_calls.len(), diploid_calls.len());
}

#[test]
fn windowed_mode_runs_the_call_stage_too() {
    let dir = tempfile::tempdir().unwrap();
    let (ref_path, r1, r2) = make_fixture(dir.path());

    let whole = dir.path().join("whole.vcf");
    run_solid(&ref_path, &r1, &r2, &whole, &["--call", "--call-ploidy", "1"]);
    let whole_text = std::fs::read_to_string(&whole).unwrap();

    let windowed = dir.path().join("windowed.vcf");
    run_solid(
        &ref_path,
        &r1,
        &r2,
        &windowed,
        &["--call", "--call-ploidy", "1", "--window-mb", "1"],
    );
    let windowed_text = std::fs::read_to_string(&windowed).unwrap();

    let sites = |t: &str| -> Vec<String> {
        data_lines(t)
            .iter()
            .map(|l| {
                let mut f = l.split('\t');
                format!(
                    "{}:{}",
                    f.next().unwrap_or(""),
                    f.next().unwrap_or("")
                )
            })
            .collect()
    };
    assert!(!sites(&whole_text).is_empty());
    assert_eq!(sites(&windowed_text), sites(&whole_text));
}
