use super::*;

fn enc(s: &[u8]) -> Vec<u8> {
    s.iter().map(|&b| encode_nt(b)).collect()
}

const REF: &[u8] = b"TTGACCGTAGACGTTCAGGCTAGCGATCCATTGCAAGTCGA";

#[test]
fn perfect_match_has_confident_states() {
    let r = enc(REF);
    let q = enc(&REF[10..26]);
    let res = glocal(&r, &q, Some(&[30; 16]), &GlocalParams { d: 0.001, e: 0.1, bw: 30 }, true).unwrap();
    assert!(res.loglik.is_finite() && res.loglik < 0.0);
    for (i, st) in res.state.iter().enumerate() {
        assert_eq!(st & 3, 0, "base {i} should be a match");
        assert_eq!((st >> 2) as usize, 10 + i, "base {i} should sit on its own column");
    }
    let interior = &res.q[2..14];
    assert!(interior.iter().all(|&q| q >= 20), "posterior qualities too low: {:?}", res.q);
}

#[test]
fn likelihood_prefers_the_matching_haplotype() {
    let hap = enc(REF);
    let read = enc(&REF[8..28]);
    let mut mism = read.clone();
    mism[5] ^= 1;
    mism[12] ^= 2;
    let p = GlocalParams { d: 0.001, e: 0.1, bw: 30 };
    let l_ok = glocal_loglik(&hap, &read, Some(&[40; 20]), &p);
    let l_bad = glocal_loglik(&hap, &mism, Some(&[40; 20]), &p);
    assert!(l_ok > l_bad + 5.0, "{l_ok} vs {l_bad}");
}

#[test]
fn reference_read_rejects_a_deletion_haplotype() {
    let read = enc(REF);
    let mut del = REF.to_vec();
    del.drain(12..14);
    let del = enc(&del);
    let p = GlocalParams { d: 0.001, e: 0.1, bw: 12 };
    let q_read = vec![35u8; read.len()];
    let q_del = vec![35u8; del.len()];
    let l_ref = glocal_loglik(&read, &read, Some(&q_read), &p);
    let l_del = glocal_loglik(&del, &read, Some(&q_read), &p);
    assert!(l_ref > l_del + 3.0, "ref {l_ref} vs deletion haplotype {l_del}");
    // And the carrier read prefers the deletion haplotype.
    let l_c_del = glocal_loglik(&del, &del, Some(&q_del), &p);
    let l_c_ref = glocal_loglik(&read, &del, Some(&q_del), &p);
    assert!(l_c_del > l_c_ref + 3.0, "carrier: del {l_c_del} vs ref {l_c_ref}");
}

#[test]
fn insertion_shows_in_map_states() {
    let r = enc(REF);
    let mut q = Vec::new();
    q.extend_from_slice(&r[..12]);
    q.extend_from_slice(&enc(b"TT"));
    q.extend_from_slice(&r[12..28]);
    let res = glocal(&r, &q, Some(&[40; 30]), &GlocalParams { d: 0.001, e: 0.1, bw: 12 }, true).unwrap();
    let ins: Vec<usize> = res.state.iter().enumerate().filter(|(_, s)| (*s & 3) == 1).map(|(i, _)| i).collect();
    assert_eq!(ins, vec![12, 13], "states {:?}", res.state);
    // Bases after the insertion map back onto their reference columns.
    assert_eq!((res.state[14] >> 2) as usize, 12);
}

#[test]
fn empty_inputs_are_rejected() {
    assert!(glocal(&[], &[0, 1], None, &GlocalParams::default(), true).is_none());
    assert!(glocal(&[0, 1], &[], None, &GlocalParams::default(), false).is_none());
}
