use super::*;

fn opts() -> Opts {
    Opts {
        output_iupac: false,
        haplotype: 0,
        allele: 0,
        isample: None,
        iupac_gts: None,
        missing_allele: None,
        absent_allele: None,
        mark_del: None,
        mark_ins: Mark::None,
        mark_snv: Mark::None,
        chr_prefix: None,
        masks: Vec::new(),
    }
}

fn rec(pos1: i64, refa: &str, alts: &[&str], gt: Option<&str>) -> Rec {
    let mut alleles = vec![refa.as_bytes().to_vec()];
    alleles.extend(alts.iter().map(|a| a.as_bytes().to_vec()));
    let gts = gt.map(|g| {
        let mut v = parse_gt(g);
        v.resize(2, GtVal::End);
        vec![v]
    });
    Rec {
        pos: pos1 - 1,
        rlen: refa.len() as i64,
        has_end: false,
        alleles,
        gts,
    }
}

fn run(opts: &Opts, fasta: &str, recs: Vec<Rec>) -> (String, String) {
    let mut vcf = HashMap::new();
    vcf.insert("1".to_string(), recs);
    let contigs = HashMap::new();
    let out = std::sync::Arc::new(std::sync::Mutex::new(Vec::<u8>::new()));
    let chain = std::sync::Arc::new(std::sync::Mutex::new(Vec::<u8>::new()));
    struct W(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
    impl Write for W {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.0.lock().unwrap().extend_from_slice(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut c = Consensus::new(
        opts,
        &vcf,
        &contigs,
        Box::new(W(out.clone())),
        Some(Box::new(W(chain.clone()))),
    );
    c.run(std::io::Cursor::new(fasta.as_bytes())).unwrap();
    let o = String::from_utf8(out.lock().unwrap().clone()).unwrap();
    let ch = String::from_utf8(chain.lock().unwrap().clone()).unwrap();
    (o, ch)
}

#[test]
fn variant_types_follow_htslib() {
    assert_eq!(variant_type(b"A", b"G"), (VCF_SNP, 0));
    assert_eq!(variant_type(b"A", b"A"), (VCF_REF, 0));
    assert_eq!(variant_type(b"AC", b"A"), (VCF_INDEL, -1));
    assert_eq!(variant_type(b"A", b"ACG"), (VCF_INDEL, 2));
    assert_eq!(variant_type(b"AC", b"GT"), (VCF_MNP, 2));
    assert_eq!(variant_type(b"A", b"<DEL>"), (VCF_OTHER, 0));
    assert_eq!(variant_type(b"A", b"*"), (VCF_OVERLAP, 0));
}

#[test]
fn genotype_text_parses_phase_and_missing() {
    assert_eq!(
        parse_gt("0/1"),
        vec![
            GtVal::Allele {
                idx: 0,
                phased: false
            },
            GtVal::Allele {
                idx: 1,
                phased: false
            }
        ]
    );
    assert_eq!(
        parse_gt(".|1"),
        vec![
            GtVal::Missing { phased: false },
            GtVal::Allele {
                idx: 1,
                phased: true
            }
        ]
    );
    assert_eq!(parse_gt("."), vec![GtVal::Missing { phased: false }]);
}

/// SNP, insertion and deletion applied to a 70-base sequence: output
/// wraps at 60 columns and the chain records the two indels.
#[test]
fn applies_snp_indel_and_writes_chain() {
    let seq: String = "ACGT".repeat(17) + "AC";
    let fasta = format!(">1\n{}\n{}\n", &seq[..40], &seq[40..]);
    let recs = vec![
        rec(2, "C", &["T"], None),
        rec(10, "C", &["CAA"], None),
        rec(30, "CGT", &["C"], None),
    ];
    let (out, chain) = run(&opts(), &fasta, recs);
    let mut want = seq.clone();
    want.replace_range(1..2, "T");
    want.replace_range(9..10, "CAA");
    want.replace_range(31..34, "C");
    let lines: Vec<&str> = out.lines().collect();
    assert_eq!(lines[0], ">1");
    assert_eq!(lines[1..].concat(), want);
    assert_eq!(lines[1].len(), 60);
    assert_eq!(
        chain,
        "chain 68 1 70 + 0 70 1 70 + 0 70 1\n10 0 2\n20 2 0\n38\n\n"
    );
}

#[test]
fn samples_and_absent_characters() {
    let fasta = ">1\nACGTACGT\n";
    // Only the sample's ALT calls apply; called REF bases are kept; absent bases become N.
    let mut o = opts();
    o.isample = Some(0);
    o.haplotype = 1;
    o.absent_allele = Some(b'N');
    let recs = vec![
        rec(2, "C", &["T"], Some("1/1")),
        rec(5, "A", &["G"], Some("0/0")),
        rec(7, "G", &["A"], Some("0/1")),
    ];
    let (out, _) = run(&o, fasta, recs);
    assert_eq!(out, ">1\nNTNNANGN\n");
    // IUPAC codes from genotypes without -H.
    let mut o = opts();
    o.iupac_gts = Some(vec![0]);
    let (out, _) = run(
        &o,
        fasta,
        vec![
            rec(2, "C", &["T"], Some("0/1")),
            rec(5, "A", &["G"], Some("1/1")),
        ],
    );
    assert_eq!(out, ">1\nAYGTGCGT\n");
}
