use super::*;
use std::io::Read;

fn tmp(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("kira_zarr_{}_{}", std::process::id(), name))
}

fn chunk_bytes(store: &Path, array: &str, key: &str) -> Vec<u8> {
    let f = std::fs::File::open(store.join(array).join(key)).expect(key);
    let mut dec = flate2::read::ZlibDecoder::new(f);
    let mut out = Vec::new();
    dec.read_to_end(&mut out).unwrap();
    out
}

fn i32s(b: &[u8]) -> Vec<i32> {
    b.chunks(4)
        .map(|c| i32::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

fn f32s(b: &[u8]) -> Vec<f32> {
    b.chunks(4)
        .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

fn strings(b: &[u8], width: usize) -> Vec<String> {
    b.chunks(width)
        .map(|c| {
            String::from_utf8_lossy(c)
                .trim_end_matches('\0')
                .to_string()
        })
        .collect()
}

const VCF: &str = "##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2>
##FILTER=<ID=q10,Description=\"low qual\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbsnp\">
##INFO=<ID=MQ,Number=1,Type=Float,Description=\"mq\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"dp\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"gq\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t10\trs1\tA\tG\t50\tPASS\tDP=20;DB;MQ=59.5\tGT:DP:GQ\t0/1:10:30\t1|1:5:.\t./.:.:3
chr1\t20\t.\tC\tT,CTT\t.\tq10\tDP=7;AF=0.1,0.2\tGT:DP\t0/2:4\t1/2:2\t0:1
chr2\t5\t.\tGA\tG\t3.5\t.\tMQ=10\tGT:GQ\t0|1:12\t.\t1/1:44
chr9\t7\trs9\tT\tA\t9\tq10;other\tDP=1\tGT\t0/0\t0/1\t1/1
";

#[test]
fn store_layout_and_values_roundtrip() {
    let dir = tmp("basic");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let vcf = dir.join("in.vcf");
    std::fs::write(&vcf, VCF).unwrap();
    let out = dir.join("out.zarr");
    let opts = ZarrOpts {
        chunk_variants: 3,
        chunk_samples: 2,
        level: 6,
        fields: None,
    };
    let s = export_vcf_zarr(&vcf, &out, &opts, false).unwrap();
    assert_eq!((s.variants, s.samples), (4, 3));
    for name in [
        "contig_id",
        "contig_length",
        "filter_id",
        "sample_id",
        "variant_contig",
        "variant_position",
        "variant_allele",
        "variant_quality",
        "variant_filter",
        "call_genotype",
        "call_genotype_phased",
        "call_genotype_mask",
        "variant_DP",
        "variant_MQ",
        "variant_DB",
        "call_DP",
        "call_GQ",
    ] {
        assert!(s.arrays.contains(&name.to_string()), "{name}");
        assert!(out.join(name).join(".zarray").exists(), "{name}/.zarray");
    }
    assert!(
        !s.arrays.contains(&"variant_AF".to_string()),
        "Number=A is not exported"
    );
    let zgroup = std::fs::read_to_string(out.join(".zgroup")).unwrap();
    assert!(zgroup.contains("\"zarr_format\": 2"));
    let attrs = std::fs::read_to_string(out.join(".zattrs")).unwrap();
    assert!(attrs.contains("\"vcf_zarr_version\": \"0.2\""));
    assert!(attrs.contains("##fileformat=VCFv4.2\\n##contig"));

    let za = std::fs::read_to_string(out.join("variant_allele").join(".zarray")).unwrap();
    assert!(za.contains("\"shape\": [4, 3]"), "{za}");
    assert!(za.contains("\"chunks\": [3, 3]"), "{za}");
    assert!(za.contains("\"dtype\": \"|S3\""), "{za}");
    assert!(za.contains("\"dimension_separator\": \".\""));
    let zg = std::fs::read_to_string(out.join("call_genotype").join(".zarray")).unwrap();
    assert!(zg.contains("\"shape\": [4, 3, 2]"), "{zg}");
    assert!(zg.contains("\"chunks\": [3, 2, 2]"), "{zg}");
    assert!(zg.contains("\"dtype\": \"|i1\""), "{zg}");
    let dims = std::fs::read_to_string(out.join("call_genotype").join(".zattrs")).unwrap();
    assert!(
        dims.contains("[\"variants\", \"samples\", \"ploidy\"]"),
        "{dims}"
    );

    // Dimension tables: header contigs first, then the unseen chr9.
    assert_eq!(
        strings(&chunk_bytes(&out, "contig_id", "0"), 4),
        ["chr1", "chr2", "chr9"]
    );
    let cl: Vec<i64> = chunk_bytes(&out, "contig_length", "0")
        .chunks(8)
        .map(|c| i64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    assert_eq!(cl, [1000, -1, -1]);
    assert_eq!(
        strings(&chunk_bytes(&out, "filter_id", "0"), 5),
        ["PASS", "q10", "other"]
    );
    assert_eq!(
        strings(&chunk_bytes(&out, "sample_id", "0"), 2),
        ["S1", "S2", "S3"]
    );

    // Variant chunk 0 holds three records, chunk 1 one record padded to three.
    assert_eq!(
        i32s(&chunk_bytes(&out, "variant_position", "0")),
        [10, 20, 5]
    );
    assert_eq!(
        i32s(&chunk_bytes(&out, "variant_position", "1")),
        [7, -1, -1]
    );
    assert_eq!(i32s(&chunk_bytes(&out, "variant_contig", "0")), [0, 0, 1]);
    assert_eq!(i32s(&chunk_bytes(&out, "variant_contig", "1")), [2, -1, -1]);
    assert_eq!(
        strings(&chunk_bytes(&out, "variant_allele", "0.0"), 3),
        ["A", "G", "", "C", "T", "CTT", "GA", "G", ""]
    );
    assert_eq!(
        strings(&chunk_bytes(&out, "variant_id", "0"), 3),
        ["rs1", "", ""]
    );
    assert_eq!(chunk_bytes(&out, "variant_id_mask", "0"), [0, 1, 1]);
    let q = f32s(&chunk_bytes(&out, "variant_quality", "0"));
    assert_eq!(q[0], 50.0);
    assert!(q[1].is_nan());
    assert_eq!(q[2], 3.5);
    // filters: PASS, q10, other per record.
    assert_eq!(
        chunk_bytes(&out, "variant_filter", "0.0"),
        [1, 0, 0, 0, 1, 0, 0, 0, 0]
    );
    assert_eq!(
        chunk_bytes(&out, "variant_filter", "1.0"),
        [0, 1, 1, 0, 0, 0, 0, 0, 0]
    );

    // INFO fields.
    assert_eq!(i32s(&chunk_bytes(&out, "variant_DP", "0")), [20, 7, -1]);
    assert_eq!(chunk_bytes(&out, "variant_DB", "0"), [1, 0, 0]);
    let mq = f32s(&chunk_bytes(&out, "variant_MQ", "0"));
    assert_eq!((mq[0], mq[2]), (59.5, 10.0));
    assert!(mq[1].is_nan());

    // Genotypes: sample blocks of two; block 1 holds S3 plus padding.
    let g00: Vec<i8> = chunk_bytes(&out, "call_genotype", "0.0.0")
        .iter()
        .map(|&b| b as i8)
        .collect();
    // A bare `.` is one missing allele; the unused ploidy slot is -2.
    assert_eq!(g00, [0, 1, 1, 1, 0, 2, 1, 2, 0, 1, -1, -2]);
    let g01: Vec<i8> = chunk_bytes(&out, "call_genotype", "0.1.0")
        .iter()
        .map(|&b| b as i8)
        .collect();
    assert_eq!(g01[..2], [-1, -1], "./. is missing");
    assert_eq!(g01[4..6], [0, -2], "haploid call padded with -2");
    assert_eq!(g01[8..10], [1, 1]);
    assert_eq!(
        chunk_bytes(&out, "call_genotype_phased", "0.0"),
        [0, 1, 0, 0, 1, 0]
    );
    assert_eq!(
        chunk_bytes(&out, "call_genotype_mask", "0.1.0")[..6],
        [1, 1, 1, 1, 0, 1]
    );
    let g10: Vec<i8> = chunk_bytes(&out, "call_genotype", "1.0.0")
        .iter()
        .map(|&b| b as i8)
        .collect();
    assert_eq!(g10[..4], [0, 0, 0, 1]);
    assert_eq!(g10[4..], vec![-1; 8], "padded variants are missing");

    // FORMAT fields.
    assert_eq!(
        i32s(&chunk_bytes(&out, "call_DP", "0.0")),
        [10, 5, 4, 2, -1, -1]
    );
    assert_eq!(i32s(&chunk_bytes(&out, "call_GQ", "0.1"))[..2], [3, -1]);

    // Explicit field selection and overwrite protection.
    let opts = ZarrOpts {
        fields: Some(vec!["INFO/MQ".into(), "FORMAT/GQ".into()]),
        ..ZarrOpts::default()
    };
    assert!(export_vcf_zarr(&vcf, &out, &opts, false).is_err());
    let s = export_vcf_zarr(&vcf, &out, &opts, true).unwrap();
    assert!(s.arrays.contains(&"variant_MQ".to_string()));
    assert!(!s.arrays.contains(&"variant_DP".to_string()));
    assert!(s.arrays.contains(&"call_GQ".to_string()));
    assert!(!out.join("call_DP").exists());
    assert_eq!(
        i32s(&chunk_bytes(&out, "variant_position", "0"))[..4],
        [10, 20, 5, 7]
    );
    let bad = ZarrOpts {
        fields: Some(vec!["INFO/AF".into()]),
        ..ZarrOpts::default()
    };
    assert!(export_vcf_zarr(&vcf, &out, &bad, true).is_err());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn sites_only_file_has_no_call_arrays() {
    let dir = tmp("sites");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let vcf = dir.join("s.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t5\t.\tA\tC\t1\t.\tDP=3\n1\t9\t.\tA\tC\t1\t.\t.\n",
    )
    .unwrap();
    let out = dir.join("s.zarr");
    let s = export_vcf_zarr(&vcf, &out, &ZarrOpts::default(), false).unwrap();
    assert_eq!((s.variants, s.samples), (2, 0));
    assert!(!s.arrays.iter().any(|a| a.starts_with("call_")));
    assert_eq!(i32s(&chunk_bytes(&out, "variant_DP", "0"))[..2], [3, -1]);
    let sid = std::fs::read_to_string(out.join("sample_id").join(".zarray")).unwrap();
    assert!(sid.contains("\"shape\": [0]"));
    let _ = std::fs::remove_dir_all(&dir);
}
