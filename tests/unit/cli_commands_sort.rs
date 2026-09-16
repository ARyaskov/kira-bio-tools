use super::*;

#[test]
fn keys_follow_header_contig_order_not_lexicographic() {
    let mut order = HashMap::new();
    order.insert("chr2".to_string(), 0);
    order.insert("chr10".to_string(), 1);
    let a = sort_key("chr10\t5\t.\tA\tT", &order);
    let b = sort_key("chr2\t500\t.\tA\tT", &order);
    assert!(b < a, "chr2 must sort before chr10 when the header says so");
    let c = sort_key("chrUn\t1\t.\tA\tT", &order);
    assert!(a < c, "unknown contigs go last");
    assert_eq!(pack_key(3, 7), (3 << 32) | 7);
}

/// Same position: alleles compared one by one, case-insensitively, shorter
/// allele lists first (bcftools' cmp_bcf_pos); then input order.
#[test]
fn ties_break_like_bcftools() {
    let item = |idx: u32, line: &str| Item {
        key: 1,
        idx,
        data: Data::Text(line.to_string()),
    };
    let a = item(0, "1\t5\t.\tA\tT,G");
    let b = item(1, "1\t5\t.\ta\tt");
    let c = item(2, "1\t5\t.\tA\tC");
    assert_eq!(cmp_items(&b, &a), Ordering::Less, "shorter ALT list first");
    assert_eq!(cmp_items(&c, &a), Ordering::Less, "C before T");
    let d = item(3, "1\t5\t.\tA\tT,G");
    assert_eq!(
        cmp_items(&a, &d),
        Ordering::Less,
        "input order on full ties"
    );
}

#[test]
fn heap_pops_smallest_first() {
    let mut order = HashMap::new();
    order.insert("1".to_string(), 0);
    let mut heap = BinaryHeap::new();
    for (i, p) in [300u32, 100, 200].iter().enumerate() {
        let line = format!("1\t{p}\t.\tA\tT");
        heap.push(HeapEntry {
            item: Item {
                key: sort_key(&line, &order),
                idx: i as u32,
                data: Data::Text(line),
            },
            src: i,
        });
    }
    let got: Vec<u64> =
        std::iter::from_fn(|| heap.pop().map(|e| e.item.key & 0xffff_ffff)).collect();
    assert_eq!(got, vec![100, 200, 300]);
}

/// Records spilled to compressed chunks come back merged in order, both as
/// text and as encoded BCF blocks.
#[test]
fn spills_and_merges_chunks() {
    let tmp = std::env::temp_dir();
    let mut order = HashMap::new();
    order.insert("1".to_string(), 0);
    let mut s = Sorter::new(false, tmp.clone(), 1);
    for p in [50u32, 10, 40, 20, 30, 20] {
        let line = format!("1\t{p}\t.\tA\tT");
        s.push(sort_key(&line, &order), Data::Text(line)).unwrap();
    }
    assert!(s.chunks.len() >= 5);
    let mut out = Vec::new();
    s.finish(
        &|d| match d {
            Data::Text(l) => sort_key(l, &order),
            _ => unreachable!(),
        },
        &mut |d| {
            if let Data::Text(l) = d {
                out.push(l.split('\t').nth(1).unwrap().parse::<u32>().unwrap());
            }
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(out, vec![10, 20, 20, 30, 40, 50]);

    let dict = crate::bcf::header::parse_header_to_dict(&[
        "##fileformat=VCFv4.2".to_string(),
        "##contig=<ID=1>".to_string(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
    ]);
    let mut s = Sorter::new(true, tmp, 1);
    for (p, alt) in [(7u32, "T"), (3, "G"), (7, "C")] {
        let mut buf = Vec::new();
        crate::bcf::record::encode_record(
            &mut buf,
            &format!("1\t{p}\t.\tA\t{alt}\t.\t.\t."),
            &dict,
        )
        .unwrap();
        let (shared, indiv) = crate::bcf::record::read_record_raw(&mut std::io::Cursor::new(buf))
            .unwrap()
            .unwrap();
        s.push(
            bcf_key(&Data::bcf(shared.clone(), indiv.clone())),
            Data::bcf(shared, indiv),
        )
        .unwrap();
    }
    let mut out = Vec::new();
    s.finish(&bcf_key, &mut |d| {
        if let Data::Bcf(f) = d {
            let line = decode_record_to_vcf(&mut std::io::Cursor::new(f.as_slice()), &dict)
                .unwrap()
                .unwrap();
            out.push(line);
        }
        Ok(())
    })
    .unwrap();
    assert_eq!(
        out,
        vec![
            "1\t3\t.\tA\tG\t.\t.\t.",
            "1\t7\t.\tA\tC\t.\t.\t.",
            "1\t7\t.\tA\tT\t.\t.\t."
        ]
    );
}

#[test]
fn max_mem_units() {
    assert_eq!(parse_max_mem("10K").unwrap(), 10 * 1024);
    assert_eq!(parse_max_mem("1M").unwrap(), 1024 * 1024);
    assert_eq!(parse_max_mem("768M").unwrap(), 768 * 1024 * 1024);
    assert!(parse_max_mem("3X").is_err());
}
