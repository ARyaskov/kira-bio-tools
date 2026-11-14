#![cfg(feature = "cli")]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

fn kira() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kira-bt"))
}

fn tool(name: &str) -> PathBuf {
    std::env::var_os(name.to_ascii_uppercase())
        .map(PathBuf::from)
        .unwrap_or_else(|| name.into())
}

fn run(exe: &Path, args: &[&str]) -> Output {
    let output = Command::new(exe)
        .args(args)
        .env("KIRA_BT_NO_PROGRESS", "1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{} {:?}: {}",
        exe.display(),
        args,
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn data(bytes: &[u8]) -> Vec<String> {
    String::from_utf8(bytes.to_vec())
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(str::to_string)
        .collect()
}

fn decode(bytes: &[u8]) -> Vec<String> {
    let mut child = Command::new(tool("bcftools"))
        .args(["view", "-H", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(bytes).unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    data(&output.stdout)
}

fn vcf(path: &Path, contigs: &[&str], rows: &str) {
    let mut text = String::from("##fileformat=VCFv4.2\n");
    for name in contigs {
        text.push_str(&format!("##contig=<ID={name},length=1000>\n"));
    }
    text.push_str("##INFO=<ID=TAG,Number=1,Type=String,Description=\"Probe\">\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    text.push_str(rows);
    std::fs::write(path, text).unwrap();
}

fn packed(path: &Path) -> PathBuf {
    let output = path.with_extension("vcf.gz");
    run(
        &tool("bcftools"),
        &[
            "view",
            "-Oz",
            "-o",
            output.to_str().unwrap(),
            path.to_str().unwrap(),
        ],
    );
    run(&tool("bcftools"), &["index", output.to_str().unwrap()]);
    output
}

#[test]
#[ignore = "requires bcftools"]
fn annotation_contigs_stdout_and_all_formats_match_bcftools() {
    let dir = tempfile::tempdir().unwrap();
    for (i, names, db_name, input_name) in [
        (0, vec!["scaf0", "scaf1"], "scaf1", "chr1"),
        (1, vec!["chr1"], "chr1", "1"),
        (2, vec!["chr1"], "chr1", "chr1"),
        (3, vec!["chrM"], "chrM", "MT"),
    ] {
        let db = dir.path().join(format!("db{i}.vcf"));
        let input = dir.path().join(format!("in{i}.vcf"));
        vcf(
            &db,
            &names,
            &format!("{db_name}\t10\t.\tA\tC\t.\tPASS\tTAG=HIT\n"),
        );
        vcf(
            &input,
            &[input_name],
            &format!("{input_name}\t10\t.\tA\tC\t.\tPASS\t.\n"),
        );
        let db_gz = packed(&db);
        let input_gz = packed(&input);
        let expected = run(
            &tool("bcftools"),
            &[
                "annotate",
                "--no-version",
                "-a",
                db_gz.to_str().unwrap(),
                "-c",
                "INFO/TAG",
                input_gz.to_str().unwrap(),
            ],
        );
        for kind in ["v", "z", "u", "b"] {
            let actual = run(
                &kira(),
                &[
                    "annotate",
                    "--no-ktile",
                    "--no-version",
                    "--threads",
                    "1",
                    "-a",
                    db.to_str().unwrap(),
                    "-c",
                    "INFO/TAG",
                    "-O",
                    kind,
                    input.to_str().unwrap(),
                ],
            );
            assert_eq!(
                decode(&actual.stdout),
                data(&expected.stdout),
                "case={i} format={kind}"
            );
        }
        if i == 2 {
            let mut child = Command::new(kira())
                .args([
                    "annotate",
                    "--no-version",
                    "--threads",
                    "1",
                    "-a",
                    db.to_str().unwrap(),
                    "-c",
                    "INFO/TAG",
                    "-Ou",
                    "-",
                ])
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .unwrap();
            child
                .stdin
                .take()
                .unwrap()
                .write_all(&std::fs::read(&input).unwrap())
                .unwrap();
            let output = child.wait_with_output().unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert_eq!(decode(&output.stdout), data(&expected.stdout));
        }
        let output = dir.path().join(format!("out{i}.bcf"));
        run(
            &kira(),
            &[
                "annotate",
                "--no-ktile",
                "--threads",
                "1",
                "-a",
                db.to_str().unwrap(),
                "-c",
                "INFO/TAG",
                "-Ob",
                "-o",
                output.to_str().unwrap(),
                "-W=csi",
                input.to_str().unwrap(),
            ],
        );
        let queried = run(
            &tool("bcftools"),
            &[
                "view",
                "-H",
                "-r",
                &format!("{input_name}:10"),
                output.to_str().unwrap(),
            ],
        );
        assert_eq!(data(&queried.stdout), data(&expected.stdout));
    }
}

#[test]
#[ignore = "requires bcftools"]
fn end_spans_survive_cross_tool_bcf_index_queries() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    vcf(
        &input,
        &["chr1"],
        "chr1\t100\t.\tA\tC\t.\tPASS\tEND=150\nchr1\t200\t.\tA\t.\t.\tPASS\tEND=250\nchr1\t300\t.\tA\t<NON_REF>\t.\tPASS\tEND=350\nchr1\t400\t.\tA\t<DEL>\t.\tPASS\tEND=450\n",
    );
    let regions = "chr1:120,chr1:220,chr1:320,chr1:420";
    for writer in [kira(), tool("bcftools")] {
        let bcf = dir.path().join(if writer == kira() {
            "kira.bcf"
        } else {
            "oracle.bcf"
        });
        run(
            &writer,
            &[
                "view",
                "-Ob",
                "-o",
                bcf.to_str().unwrap(),
                input.to_str().unwrap(),
            ],
        );
        for indexer in [kira(), tool("bcftools")] {
            run(&indexer, &["index", "-f", bcf.to_str().unwrap()]);
            for reader in [kira(), tool("bcftools")] {
                let rows = run(
                    &reader,
                    &["view", "-H", "-r", regions, bcf.to_str().unwrap()],
                );
                assert_eq!(
                    data(&rows.stdout).len(),
                    4,
                    "writer={} indexer={} reader={}",
                    writer.display(),
                    indexer.display(),
                    reader.display()
                );
            }
        }
    }
}

fn bam_fixture(dir: &Path, name: &str, reverse_dictionary: bool) -> PathBuf {
    let sam = dir.join(format!("{name}.sam"));
    let sq = if reverse_dictionary {
        "@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr1\tLN:1000\n"
    } else {
        "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n"
    };
    let mut text = format!("@HD\tVN:1.6\tSO:coordinate\n{sq}@RG\tID:{name}\tSM:{name}\n");
    for i in 0..12 {
        let seq = if i < 6 {
            "AAAAAAAAAACAAAAAAAAA"
        } else {
            "AAAAAAAAAAAAAAAAAAAA"
        };
        text.push_str(&format!(
            "r{i}\t0\tchr1\t100\t60\t20M\t*\t0\t0\t{seq}\tIIIIIIIIIIIIIIIIIIII\tRG:Z:{name}\n"
        ));
    }
    std::fs::write(&sam, text).unwrap();
    let bam = dir.join(format!("{name}.bam"));
    run(
        &tool("samtools"),
        &[
            "view",
            "-b",
            "-o",
            bam.to_str().unwrap(),
            sam.to_str().unwrap(),
        ],
    );
    bam
}

fn genotype_fields(bytes: &[u8]) -> Vec<Vec<String>> {
    data(bytes)
        .iter()
        .map(|row| {
            let cols: Vec<&str> = row.split('\t').collect();
            let mut out: Vec<String> = cols[..5].iter().map(|s| s.to_string()).collect();
            let keys: Vec<&str> = cols[8].split(':').collect();
            for sample in &cols[9..] {
                let fields: Vec<&str> = sample.split(':').collect();
                for tag in ["PL", "AD", "DP"] {
                    out.push(fields[keys.iter().position(|k| *k == tag).unwrap()].to_string());
                }
            }
            out
        })
        .collect()
}

#[test]
#[ignore = "requires bcftools and samtools"]
fn pileup_stdout_formats_parallelism_and_bam_validation() {
    let dir = tempfile::tempdir().unwrap();
    let reference = dir.path().join("ref.fa");
    std::fs::write(
        &reference,
        format!(">chr1\n{}\n>chr2\n{}\n", "A".repeat(1000), "C".repeat(1000)),
    )
    .unwrap();
    run(&tool("samtools"), &["faidx", reference.to_str().unwrap()]);
    let bam = bam_fixture(dir.path(), "s1", false);
    let expected = run(
        &tool("bcftools"),
        &[
            "mpileup",
            "-B",
            "-f",
            reference.to_str().unwrap(),
            "-a",
            "FORMAT/AD,FORMAT/DP",
            bam.to_str().unwrap(),
        ],
    );
    let expected_fields = genotype_fields(&expected.stdout);
    assert!(!expected_fields.is_empty());
    for (threads, mode, kind) in [
        ("1", false, "v"),
        ("2", false, "b"),
        ("1", true, "u"),
        ("2", false, "z"),
    ] {
        let mut args = vec![
            "mpileup",
            "--bcftools-layout",
            "-B",
            "-f",
            reference.to_str().unwrap(),
            "-a",
            "FORMAT/AD,FORMAT/DP",
            "--threads",
            threads,
            "-O",
            kind,
            bam.to_str().unwrap(),
        ];
        if mode {
            args.push("--stream");
        }
        let actual = run(&kira(), &args);
        let decoded = decode(&actual.stdout).join("\n");
        assert_eq!(
            genotype_fields(decoded.as_bytes()),
            expected_fields,
            "threads={threads} stream={mode} format={kind}"
        );
    }
    let actual = run(&kira(), &["samview", "-c", bam.to_str().unwrap()]);
    let expected = run(&tool("samtools"), &["view", "-c", bam.to_str().unwrap()]);
    assert_eq!(actual.stdout, expected.stdout);
    let output = dir.path().join("pileup.bcf");
    run(
        &kira(),
        &[
            "mpileup",
            "--bcftools-layout",
            "-B",
            "-f",
            reference.to_str().unwrap(),
            "-a",
            "FORMAT/AD,FORMAT/DP",
            "-Ob",
            "-o",
            output.to_str().unwrap(),
            "-W=csi",
            bam.to_str().unwrap(),
        ],
    );
    let queried = run(
        &tool("bcftools"),
        &["view", "-H", "-r", "chr1:110-110", output.to_str().unwrap()],
    );
    assert_eq!(
        genotype_fields(&queried.stdout),
        expected_fields
            .iter()
            .filter(|row| row[1] == "110")
            .cloned()
            .collect::<Vec<_>>()
    );
    let other = bam_fixture(dir.path(), "s2", true);
    let failed = Command::new(kira())
        .args([
            "mpileup",
            "-B",
            "-f",
            reference.to_str().unwrap(),
            bam.to_str().unwrap(),
            other.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(
        String::from_utf8_lossy(&failed.stderr).contains("incompatible reference dictionaries")
    );
}

#[test]
#[ignore = "requires bcftools"]
fn unsorted_overlapping_masks_match_bcftools() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    vcf(
        &input,
        &["chr1"],
        "chr1\t5\t.\tAA\tA\t.\tPASS\t.\nchr1\t10\t.\tA\tC\t.\tPASS\t.\nchr1\t20\t.\tA\tC\t.\tPASS\t.\nchr1\t30\t.\tA\tC\t.\tPASS\t.\n",
    );
    for overlap in ["0", "1", "2"] {
        let args = [
            "filter",
            "--no-version",
            "--soft-filter",
            "MASK",
            "--mask",
            "chr1:19-21,chr1:4-9,chr1:8-11",
            "--mask-overlap",
            overlap,
            input.to_str().unwrap(),
        ];
        assert_eq!(
            data(&run(&kira(), &args).stdout),
            data(&run(&tool("bcftools"), &args).stdout),
            "overlap={overlap}"
        );
    }
    vcf(&input, &["chr1"], "chr1\t9\t.\tA\tAT\t.\tPASS\t.\n");
    for mask in ["chr1:9-9,chr1:10-10", "chr1:9-10"] {
        let args = [
            "filter",
            "--no-version",
            "--soft-filter",
            "MASK",
            "--mask",
            mask,
            "--mask-overlap",
            "2",
            input.to_str().unwrap(),
        ];
        assert_eq!(
            data(&run(&kira(), &args).stdout),
            data(&run(&tool("bcftools"), &args).stdout)
        );
    }
}

#[test]
#[ignore = "requires bcftools"]
fn indel_flags_and_extents_survive_cross_tool_indexed_queries() {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/mpileup/test25");
    let dir = tempfile::tempdir().unwrap();
    let reference = dir.path().join("ref.fa");
    std::fs::copy(fixture.join("ref.fa"), &reference).unwrap();
    let first = fixture.join("annot-NMBZ.3.1.bam");
    let second = fixture.join("annot-NMBZ.3.2.bam");
    let expected: Vec<String> = include_str!("fixtures/mpileup/test25.records.tsv")
        .lines()
        .map(str::to_string)
        .collect();
    for (kind, name) in [("b", "indels.bcf"), ("z", "indels.vcf.gz")] {
        let output = dir.path().join(name);
        run(
            &kira(),
            &[
                "mpileup",
                first.to_str().unwrap(),
                second.to_str().unwrap(),
                "-f",
                reference.to_str().unwrap(),
                "-r",
                "chr16:75",
                "-a",
                "-AD,INFO/NMBZ",
                "-O",
                kind,
                "-o",
                output.to_str().unwrap(),
                "-W=csi",
            ],
        );
        let records = run(&tool("bcftools"), &["view", "-H", output.to_str().unwrap()]);
        assert_eq!(data(&records.stdout), expected, "{kind}");
        let queried = run(
            &tool("bcftools"),
            &["view", "-H", "-r", "chr16:80", output.to_str().unwrap()],
        );
        assert_eq!(data(&queried.stdout), expected[1..], "{kind}");
    }
}

#[test]
#[ignore = "requires bcftools"]
fn cnv_smoothing_matches_defined_bcftools_windows() {
    let work = tempfile::tempdir().unwrap();
    let summarize = |directory: &Path| {
        [
            "summary.tab",
            "summary.Q.tab",
            "summary.C.tab",
            "cn.Q.tab",
            "cn.C.tab",
            "dat.Q.tab",
            "dat.C.tab",
        ]
        .into_iter()
        .map(|name| {
            let bytes = std::fs::read(directory.join(name)).unwrap();
            data(&bytes)
        })
        .collect::<Vec<_>>()
    };
    for count in [1usize, 2, 3, 4, 5, 12, 32] {
        let input = work.path().join(format!("n{count}.vcf"));
        let mut text = "##fileformat=VCFv4.2\n##contig=<ID=1,length=10000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=BAF,Number=1,Type=Float,Description=\"B allele frequency\">\n##FORMAT=<ID=LRR,Number=1,Type=Float,Description=\"Log R ratio\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQ\tC\n".to_string();
        for index in 0..count {
            let baf = ["0.4", "0.5", "0.6"][index % 3];
            let lrr = ["-0.1", "0.3", "0.1"][index % 3];
            text.push_str(&format!(
                "1\t{}\t.\tA\tG\t.\tPASS\t.\tGT:BAF:LRR\t0/1:{baf}:{lrr}\t0/1:0.5:0\n",
                100 + index * 10
            ));
        }
        std::fs::write(&input, text).unwrap();
        for window in [3usize, 10] {
            let reference_window = window.min(count * 2).to_string();
            let window_text = window.to_string();
            let expected = work.path().join(format!("reference-{count}-{window}"));
            let actual = work.path().join(format!("kira-{count}-{window}"));
            run(
                &tool("bcftools"),
                &[
                    "cnv",
                    "-s",
                    "Q",
                    "-c",
                    "C",
                    "-L",
                    &reference_window,
                    "-l",
                    "0.8",
                    "-o",
                    expected.to_str().unwrap(),
                    input.to_str().unwrap(),
                ],
            );
            let mut args = vec![
                "cnv",
                "-s",
                "Q",
                "-c",
                "C",
                "-l",
                "0.8",
                "-o",
                actual.to_str().unwrap(),
                input.to_str().unwrap(),
            ];
            if window != 10 {
                args.extend(["-L", window_text.as_str()]);
            }
            run(&kira(), &args);
            assert_eq!(
                summarize(&actual),
                summarize(&expected),
                "sites={count}, window={window}, reference_window={reference_window}"
            );
        }
    }
}
