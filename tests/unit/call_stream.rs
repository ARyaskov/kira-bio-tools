    use super::*;
    use crate::vcf::parse_vcf_line;

    const HEADER: &[&str] = &[
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled likelihoods\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    ];

    /// One hom-alt site: PL puts all the weight on 1/1.
    const SITE: &str = "chr1\t100\t.\tA\tG\t.\t.\tDP=30\tGT:PL\t0/0:255,255,0";

    fn headers() -> Vec<String> {
        HEADER.iter().map(|s| s.to_string()).collect()
    }

    fn call_one(cfg: CallConfig, line: &str) -> String {
        let rec = parse_vcf_line(line).expect("parse");
        let mut out: Vec<u8> = Vec::new();
        call_stream(vec![rec], &headers(), cfg, &mut out).expect("call");
        String::from_utf8(out).expect("utf8")
    }

    fn data_lines(vcf: &str) -> Vec<&str> {
        vcf.lines().filter(|l| !l.starts_with('#')).collect()
    }

    fn sample_col(line: &str) -> &str {
        line.split('\t').nth(9).expect("sample column")
    }

    fn info_of(line: &str) -> &str {
        line.split('\t').nth(7).expect("info column")
    }

    #[test]
    fn diploid_site_keeps_the_pair() {
        let vcf = call_one(CallConfig::default(), SITE);
        let calls = data_lines(&vcf);
        assert_eq!(calls.len(), 1);
        assert_eq!(sample_col(calls[0]), "1/1");
        assert!(info_of(calls[0]).contains("AN=2"), "{}", info_of(calls[0]));
        assert!(info_of(calls[0]).contains("AC=2"), "{}", info_of(calls[0]));
    }

    #[test]
    fn haploid_site_emits_one_allele() {
        let cfg = CallConfig { ploidy: 1, ..CallConfig::default() };
        let vcf = call_one(cfg, SITE);
        let calls = data_lines(&vcf);
        assert_eq!(calls.len(), 1);
        assert_eq!(sample_col(calls[0]), "1");
        // A haploid sample contributes exactly one allele to AC/AN.
        assert!(info_of(calls[0]).contains("AN=1"), "{}", info_of(calls[0]));
        assert!(info_of(calls[0]).contains("AC=1"), "{}", info_of(calls[0]));
    }

    #[test]
    fn ploidy_regions_override_the_default() {
        let cfg = CallConfig {
            ploidy: 2,
            ploidy_regions: vec![PloidyRegion {
                chrom: "chr1".into(),
                beg: 1,
                end: 1000,
                sex: "*".into(),
                ploidy: 1,
            }],
            ..CallConfig::default()
        };
        let vcf = call_one(cfg.clone(), SITE);
        assert_eq!(sample_col(data_lines(&vcf)[0]), "1");

        // Outside the region the uniform ploidy applies again.
        let far = SITE.replace("\t100\t", "\t5000\t");
        let vcf = call_one(cfg, &far);
        assert_eq!(sample_col(data_lines(&vcf)[0]), "1/1");
    }

    #[test]
    fn consensus_caller_rejects_non_diploid() {
        let cfg = CallConfig {
            mode: CallMode::Consensus,
            ploidy: 1,
            ..CallConfig::default()
        };
        assert!(cfg.validate().is_err());
    }

    #[test]
    fn ploidy_above_two_is_rejected() {
        let cfg = CallConfig { ploidy: 3, ..CallConfig::default() };
        assert!(cfg.validate().is_err());
    }

    /// A writer splits wherever its buffer ends, mid-line included.
    #[test]
    fn sink_reassembles_lines_split_across_writes() {
        let mut vcf_text = HEADER.join("\n");
        vcf_text.push('\n');
        vcf_text.push_str(SITE);
        vcf_text.push('\n');

        let whole = {
            let mut out: Vec<u8> = Vec::new();
            let mut sink = CallSink::new(&mut out, CallConfig::default()).unwrap();
            sink.write_all(vcf_text.as_bytes()).unwrap();
            sink.finish().unwrap();
            String::from_utf8(out).unwrap()
        };

        for chunk in [1usize, 3, 7, 64] {
            let mut out: Vec<u8> = Vec::new();
            let mut sink = CallSink::new(&mut out, CallConfig::default()).unwrap();
            for part in vcf_text.as_bytes().chunks(chunk) {
                sink.write_all(part).unwrap();
            }
            sink.finish().unwrap();
            let got = String::from_utf8(out).unwrap();
            assert_eq!(got, whole, "chunk size {chunk}");
        }
    }

    /// A producer that does not end its last line must not lose that site.
    #[test]
    fn sink_flushes_an_unterminated_last_line() {
        let mut vcf_text = HEADER.join("\n");
        vcf_text.push('\n');
        vcf_text.push_str(SITE);

        let mut out: Vec<u8> = Vec::new();
        let mut sink = CallSink::new(&mut out, CallConfig::default()).unwrap();
        sink.write_all(vcf_text.as_bytes()).unwrap();
        sink.finish().unwrap();
        let got = String::from_utf8(out).unwrap();
        assert_eq!(data_lines(&got).len(), 1);
    }

    /// The sink emits a usable header even when the producer found no sites.
    #[test]
    fn sink_writes_header_for_an_empty_call_set() {
        let mut out: Vec<u8> = Vec::new();
        let mut sink = CallSink::new(&mut out, CallConfig::default()).unwrap();
        for h in HEADER {
            sink.write_all(h.as_bytes()).unwrap();
            sink.write_all(b"\n").unwrap();
        }
        sink.finish().unwrap();
        let got = String::from_utf8(out).unwrap();
        assert!(got.starts_with("##fileformat=VCFv4.2"));
        assert!(got.lines().any(|l| l.starts_with("#CHROM")));
        assert!(data_lines(&got).is_empty());
    }

    #[test]
    fn variants_only_drops_reference_sites() {
        let ref_site = "chr1\t100\t.\tA\t.\t.\t.\tDP=30\tGT:PL\t0/0:0";
        let cfg = CallConfig { variants_only: true, ..CallConfig::default() };
        let vcf = call_one(cfg, ref_site);
        assert!(data_lines(&vcf).is_empty(), "{vcf}");
    }

    #[test]
    fn push_before_header_is_an_error() {
        let mut out: Vec<u8> = Vec::new();
        let mut stream = CallStream::new(&mut out, CallConfig::default()).unwrap();
        let rec = parse_vcf_line(SITE).unwrap();
        assert!(stream.push(&rec).is_err());
    }
