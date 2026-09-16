use super::*;

const BASES: [u8; 4] = *b"ACGT";

fn lcg(seed: &mut u64) -> u32 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 33) as u32
}

fn random_ref(len: usize, seed: u64) -> Vec<u8> {
    let mut s = seed;
    (0..len)
        .map(|_| BASES[(lcg(&mut s) % 4) as usize])
        .collect()
}

fn with_deletion(refw: &[u8], at: usize, len: usize) -> Vec<u8> {
    let mut h = refw[..at].to_vec();
    h.extend_from_slice(&refw[at + len..]);
    h
}

fn with_insertion(refw: &[u8], at: usize, bases: &[u8]) -> Vec<u8> {
    let mut h = refw[..at].to_vec();
    h.extend_from_slice(bases);
    h.extend_from_slice(&refw[at..]);
    h
}

fn q(len: usize) -> Vec<u8> {
    vec![35; len]
}

#[test]
fn graph_paths_recover_the_haplotypes_present_in_the_reads() {
    let refw = random_ref(80, 5);
    let del = with_deletion(&refw, 40, 3);
    let ins = with_insertion(&refw, 25, b"TTAG");
    let mut reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    for _ in 0..5 {
        reads.push((refw.clone(), q(refw.len())));
        reads.push((del.clone(), q(del.len())));
        reads.push((ins.clone(), q(ins.len())));
    }
    // One read with a sequencing error and one with a low-quality stretch.
    let mut err = del.clone();
    err[10] = if err[10] == b'A' { b'C' } else { b'A' };
    reads.push((err, q(del.len())));
    let mut lowq = q(refw.len());
    for x in lowq[30..50].iter_mut() {
        *x = 2;
    }
    reads.push((refw.clone(), lowq));
    let views: Vec<(&[u8], &[u8])> = reads
        .iter()
        .map(|(b, q)| (b.as_slice(), q.as_slice()))
        .collect();

    let g = Graph::build(21, &refw, &views, 10).expect("unique reference k-mers");
    let haps = g.haplotypes(2, refw.len() + 50);
    let seqs: Vec<&[u8]> = haps.iter().map(|(h, _)| h.as_slice()).collect();
    assert!(seqs.contains(&refw.as_slice()));
    assert!(seqs.contains(&del.as_slice()));
    assert!(seqs.contains(&ins.as_slice()));
    let sup = |s: &[u8]| haps.iter().find(|(h, _)| h == s).unwrap().1;
    assert_eq!(sup(&refw), u32::MAX);
    assert_eq!(
        sup(&del),
        6,
        "five clean copies plus the read with an error elsewhere"
    );
    assert_eq!(sup(&ins), 5);
    // The lone error k-mers are below min_kmer and open no path of their own.
    assert_eq!(
        haps.len(),
        3,
        "{:?}",
        haps.iter().map(|(h, s)| (h.len(), *s)).collect::<Vec<_>>()
    );

    let p = AsmParams::default();
    let alts = candidate_haplotypes(&refw, &views, &p);
    assert_eq!(alts.len(), 2);
    assert_eq!(alts[0], del, "better supported haplotype first");
    assert_eq!(alts[1], ins);
}

#[test]
fn events_are_left_normalised_and_anchored() {
    // Homopolymer: deleting any A of the run is the same event, reported at
    // the leftmost position.
    let refw = b"CCGTAAAAAGTCCG".to_vec();
    let hap = with_deletion(&refw, 6, 2);
    let ev = indel_events(&refw, &hap, 50);
    assert_eq!(ev.len(), 1);
    assert!(!ev[0].is_ins);
    assert_eq!((ev[0].ref_off, ev[0].len), (4, 2));
    let (anchor, r, a) = indel_alleles(&refw, &ev[0]).unwrap();
    assert_eq!((anchor, r.as_str(), a.as_str()), (3, "TAA", "T"));

    // Insertion of a repeat unit inside a dinucleotide run shifts left too.
    let refw = b"GGCATATATGGC".to_vec();
    let hap = with_insertion(&refw, 7, b"AT");
    let ev = indel_events(&refw, &hap, 50);
    assert_eq!(ev.len(), 1);
    assert!(ev[0].is_ins);
    assert_eq!(ev[0].ref_off, 3);
    assert_eq!(ev[0].bases, b"AT");

    // A haplotype differing by SNPs only has no events; a deletion flush
    // with the window end is not anchored and dropped.
    let refw = random_ref(60, 9);
    let mut snp = refw.clone();
    snp[20] = if snp[20] == b'G' { b'T' } else { b'G' };
    assert!(indel_events(&refw, &snp, 50).is_empty());
    // (The last five bases could shift one step left, so cut four.)
    let tail = refw[..56].to_vec();
    assert!(indel_events(&refw, &tail, 50).is_empty());
    // Two events on one haplotype are both reported, in order.
    let two = with_insertion(&with_deletion(&refw, 45, 4), 15, b"GCA");
    let ev = indel_events(&refw, &two, 50);
    assert_eq!(ev.len(), 2);
    assert!(ev[0].is_ins && ev[0].ref_off <= 15);
    assert!(!ev[1].is_ins && ev[1].len == 4);
}

#[test]
fn haplotype_pair_genotyping_separates_carriers() {
    let refw = random_ref(70, 21);
    let del = with_deletion(&refw, 33, 2);
    // A carrier haplotype also has a SNP: reads carrying both still score
    // as clean deletion carriers.
    let mut del_snp = del.clone();
    del_snp[10] = if del_snp[10] == b'A' { b'C' } else { b'A' };
    let mut reads: Vec<(Vec<u8>, Vec<u8>, usize)> = Vec::new();
    for _ in 0..6 {
        reads.push((refw.clone(), q(70), 0)); // sample 0: hom ref
        reads.push((refw.clone(), q(70), 1)); // sample 1: het
        reads.push((del_snp.clone(), q(68), 1));
        reads.push((del_snp.clone(), q(68), 2)); // sample 2: hom alt
    }
    let alts = vec![del_snp.clone()];
    let calls = calls_from_haplotypes(&refw, 1000, &alts, &reads, &AsmParams::default());
    assert_eq!(calls.len(), 1);
    let c = &calls[0];
    assert_eq!(c.pos1, 1000 + 33);
    assert_eq!(c.ref_str.len(), 3);
    assert_eq!(c.alt_str.len(), 1);
    assert_eq!((c.support, c.total), (12, 24));
    let pls = crate::call::haplotype::haplotype_pls(3, c);
    assert_eq!(pls[0].pl[0], 0, "{:?}", pls[0]);
    assert!(pls[0].pl[1] > 10 && pls[0].pl[2] > 30, "{:?}", pls[0]);
    assert_eq!(pls[1].pl[1], 0, "{:?}", pls[1]);
    assert!(pls[1].pl[0] > 10 && pls[1].pl[2] > 10, "{:?}", pls[1]);
    assert_eq!(pls[2].pl[2], 0, "{:?}", pls[2]);
    assert!(pls[2].pl[0] > 30, "{:?}", pls[2]);
    assert_eq!((pls[1].n_ref, pls[1].n_alt), (6, 6));
}

#[test]
fn assemble_window_end_to_end_with_live_reads() {
    use crate::bam::pileup::CigarOps;
    use noodles_sam::alignment::record::cigar::op::Kind;
    let refw = random_ref(60, 33);
    let ins = with_insertion(&refw, 28, b"GATC");
    let mk = |seq: &[u8], cigar: Vec<(Kind, u32)>, sample: usize| {
        LiveRead::new(
            seq,
            &vec![35; seq.len()],
            cigar.into_iter().collect::<CigarOps>(),
            0,
            0,
            60,
            sample,
            0,
        )
    };
    let mut reads = Vec::new();
    for _ in 0..4 {
        reads.push(mk(&refw, vec![(Kind::Match, 60)], 0));
        reads.push(mk(
            &ins,
            vec![(Kind::Match, 28), (Kind::Insertion, 4), (Kind::Match, 32)],
            1,
        ));
        reads.push(mk(
            &ins,
            vec![(Kind::Match, 28), (Kind::Insertion, 4), (Kind::Match, 32)],
            0,
        ));
    }
    let calls = assemble_window(&reads, 0, &refw, &AsmParams::default());
    assert_eq!(calls.len(), 1);
    let c = &calls[0];
    assert_eq!(c.pos1, 28);
    assert_eq!(c.ref_str.len(), 1);
    assert_eq!(c.alt_str.len(), 5);
    assert_eq!(&c.alt_str.as_bytes()[1..], b"GATC");
    assert_eq!((c.support, c.total), (8, 12));
    let pls = crate::call::haplotype::haplotype_pls(2, c);
    assert_eq!(pls[0].pl[1], 0);
    assert_eq!(pls[1].pl[2], 0);
    // Too few spanning reads: nothing is assembled.
    assert!(assemble_window(&reads[..2], 0, &refw, &AsmParams::default()).is_empty());
}
