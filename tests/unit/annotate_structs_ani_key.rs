    use super::*;

    #[test]
    fn xor_inversion_no_longer_collides() {
        // The classical XOR-commutative bug: A>T and T>A used to collide.
        let k_at = make_variant_key(1, 1000, b"A", b"T");
        let k_ta = make_variant_key(1, 1000, b"T", b"A");
        assert_ne!(k_at, k_ta, "A>T and T>A must hash to distinct keys");
    }

    #[test]
    fn length_prefix_separates_boundary_ambiguity() {
        // Without the length prefix + 0xFF separator these would collide:
        // ref="AA"+alt="A" vs ref="A"+alt="AA".
        let k1 = make_variant_key(1, 1000, b"AA", b"A");
        let k2 = make_variant_key(1, 1000, b"A", b"AA");
        assert_ne!(k1, k2);
    }

    #[test]
    fn deterministic_across_runs() {
        // Same input must always produce the same hash — critical because
        // build-time and lookup-time use the same function but on different
        // process invocations.
        let a = make_variant_key(7, 5_000_001, b"ACGT", b"A");
        let b = make_variant_key(7, 5_000_001, b"ACGT", b"A");
        assert_eq!(a, b);
    }

    #[test]
    fn different_chroms_distinct() {
        let k1 = make_variant_key(1, 1000, b"A", b"T");
        let k2 = make_variant_key(2, 1000, b"A", b"T");
        assert_ne!(k1, k2);
    }
