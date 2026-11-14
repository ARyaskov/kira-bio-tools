use super::AfFile;

fn alleles(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

#[test]
fn indexed_sites_preserve_first_rows_allele_order_and_deferred_errors() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("targets.tsv");
    let content = concat!(
        "# targets\r\n\r\n",
        "contig_λ\t9\tC,T,G\t0.25\r\n",
        "chr1\t0\tA,C\t.\r\n",
        "chr1\t17\tA,C\tinvalid\r\n",
        "chr1\t17\tC,T\t0.5\r\n",
        "chr1\t19\tA\t0.4\r\n",
        "chr1\t23\tA,C\t1.1\r\n",
        "chr1\t25\tA,C\r\n",
        "chr1\t27\tA,C\t0.5 trailing\r\n",
        "chr1\t4294967295\tA,C\t0.75\r\n",
        "contig_λ\t1\tA,C\t0\r\n",
    );
    std::fs::write(&path, content).unwrap();
    let file = AfFile::load(&path).unwrap();
    for (chrom, pos) in [
        ("chr1", 0),
        ("chr1", 17),
        ("chr1", u32::MAX),
        ("contig_λ", 9),
    ] {
        assert!(file.has_site(chrom, pos));
    }
    for (chrom, pos) in [("1", 17), ("chr1", 18), ("contig_Λ", 9), ("contig_λ", 17)] {
        assert!(!file.has_site(chrom, pos));
        assert_eq!(
            file.read_af(chrom, pos, &alleles(&["A", "C"])).unwrap(),
            None
        );
    }
    assert_eq!(
        file.read_af("contig_λ", 9, &alleles(&["C", "T", "G"]))
            .unwrap(),
        Some(0.25)
    );
    assert_eq!(
        file.read_af("contig_λ", 9, &alleles(&["C", "G", "T"]))
            .unwrap(),
        None
    );
    assert_eq!(
        file.read_af("contig_λ", 9, &alleles(&["C", "T"])).unwrap(),
        None
    );
    assert_eq!(
        file.read_af("chr1", 0, &alleles(&["A", "C"])).unwrap(),
        None
    );
    assert_eq!(
        file.read_af("chr1", u32::MAX, &alleles(&["A", "C"]))
            .unwrap(),
        Some(0.75)
    );
    assert_eq!(
        file.read_af("contig_λ", 1, &alleles(&["A", "C"])).unwrap(),
        Some(0.0)
    );
    assert_eq!(
        file.read_af("chr1", 27, &alleles(&["A", "C"])).unwrap(),
        Some(0.5)
    );
    assert_eq!(
        file.read_af("chr1", 17, &alleles(&["C", "T"])).unwrap(),
        None
    );
    assert_eq!(
        file.read_af("chr1", 23, &alleles(&["A", "G"])).unwrap(),
        None
    );
    assert_eq!(
        file.read_af("chr1", 25, &alleles(&["C", "A"])).unwrap(),
        None
    );
    let error = file.read_af("chr1", 17, &alleles(&["A", "C"])).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Could not parse: [chr1\t17\tA,C\tinvalid]"
    );
    let error = file.read_af("chr1", 23, &alleles(&["A", "C"])).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Could not parse AF: [chr1\t23\tA,C\t1.1]"
    );
    let error = file.read_af("chr1", 25, &alleles(&["A", "C"])).unwrap_err();
    assert_eq!(error.to_string(), "Could not parse: [chr1\t25\tA,C]");
    let error = file.read_af("chr1", 19, &alleles(&["T", "G"])).unwrap_err();
    assert_eq!(
        error.to_string(),
        format!(
            "Expected two comma-separated alleles (REF,ALT) in the third column of {}, found:\n\tchr1\t19\tA\t0.4",
            path.display()
        )
    );
}

#[test]
fn indexed_site_queries_match_an_independent_first_row_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("many-targets.tsv");
    let contigs = ["chr1", "chr01", "chrX", "contig_λ"];
    let mut rows = Vec::new();
    for (chrom_index, chrom) in contigs.iter().enumerate() {
        for i in 0..257u32 {
            let pos = if i == 256 { u32::MAX } else { i * 17 };
            let pair = if i % 2 == 0 { ["A", "C"] } else { ["G", "T"] };
            let af = [0.0, 0.25, 0.75, 1.0][(i as usize + chrom_index) % 4];
            rows.push(((*chrom).to_string(), pos, alleles(&pair), af));
        }
    }
    rows.reverse();
    for (chrom, pos, _, _) in rows.clone().into_iter().step_by(31) {
        rows.push((chrom, pos, alleles(&["C", "G"]), 0.5));
    }
    let mut content = String::new();
    for (chrom, pos, pair, af) in &rows {
        content.push_str(&format!("{chrom}\t{pos}\t{}\t{af}\n", pair.join(",")));
    }
    std::fs::write(&path, content).unwrap();
    let file = AfFile::load(&path).unwrap();
    for chrom in ["chr1", "chr01", "chrX", "contig_λ", "unknown"] {
        for pos in (0..=4370).chain([u32::MAX - 1, u32::MAX]) {
            let expected = rows.iter().find(|row| row.0 == chrom && row.1 == pos);
            assert_eq!(
                file.has_site(chrom, pos),
                expected.is_some(),
                "{chrom}:{pos}"
            );
            for query in [
                alleles(&["A", "C"]),
                alleles(&["G", "T"]),
                alleles(&["C", "G"]),
            ] {
                let expected_af = expected.filter(|row| row.2 == query).map(|row| row.3);
                assert_eq!(
                    file.read_af(chrom, pos, &query).unwrap(),
                    expected_af,
                    "{chrom}:{pos} {query:?}"
                );
            }
        }
    }
}

#[test]
fn indexed_site_files_accept_empty_and_reject_invalid_coordinates() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("empty.tsv");
    std::fs::write(&path, "# no targets\n\n").unwrap();
    let file = AfFile::load(&path).unwrap();
    assert!(!file.has_site("chr1", 1));
    assert_eq!(
        file.read_af("chr1", 1, &alleles(&["A", "C"])).unwrap(),
        None
    );
    for line in [
        "chr1\t1",
        "chr1\t-1\tA,C\t0.5",
        "chr1\t4294967296\tA,C\t0.5",
    ] {
        std::fs::write(&path, line).unwrap();
        let error = match AfFile::load(&path) {
            Ok(_) => panic!("invalid target accepted: {line}"),
            Err(error) => error,
        };
        assert_eq!(
            error.to_string(),
            format!(
                "Could not parse the targets file {}: {line}",
                path.display()
            )
        );
    }
}
