    use super::*;

    #[test]
    fn keys_follow_header_contig_order_not_lexicographic() {
        let mut order = HashMap::new();
        order.insert("chr2".to_string(), 0);
        order.insert("chr10".to_string(), 1);
        let a = sort_key("chr10\t5\t.\tA\tT", 0, &order);
        let b = sort_key("chr2\t500\t.\tA\tT", 1, &order);
        assert!(b < a, "chr2 must sort before chr10 when the header says so");
        let c = sort_key("chrUn\t1\t.\tA\tT", 2, &order);
        assert!(a < c, "unknown contigs go last");
    }

    #[test]
    fn heap_pops_smallest_first() {
        let mut order = HashMap::new();
        order.insert("1".to_string(), 0);
        let mut heap = BinaryHeap::new();
        for (i, p) in [300u32, 100, 200].iter().enumerate() {
            let line = format!("1\t{p}\t.\tA\tT");
            heap.push(HeapEntry { key: sort_key(&line, i, &order), line, src: i });
        }
        let got: Vec<u32> = std::iter::from_fn(|| heap.pop().map(|e| e.key.pos)).collect();
        assert_eq!(got, vec![100, 200, 300]);
    }

    #[test]
    fn max_mem_units() {
        assert_eq!(parse_max_mem("10K").unwrap(), 10 * 1024);
        assert_eq!(parse_max_mem("1M").unwrap(), 1024 * 1024);
        assert_eq!(parse_max_mem("768M").unwrap(), 768 * 1024 * 1024);
        assert!(parse_max_mem("3X").is_err());
    }
