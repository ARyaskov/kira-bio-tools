    use super::*;
    use crate::annotate::cpu_v2::{ParsedFormat, ParsedSample};
    use std::collections::HashMap;

    #[test]
    fn test_empty_sample_column_is_dot() {
        let vals = vec![String::new()];
        assert_eq!(normalize_sample_values(&vals), ".");
    }

    #[test]
    fn test_missing_gt_with_format_fields() {
        let vals = vec![String::new(), String::new(), String::new(), String::new()];
        assert_eq!(normalize_sample_values(&vals), ".:.:.:.");
    }

    #[test]
    fn test_missing_subfields_are_dots() {
        let vals = vec![
            "0/0".to_string(),
            "".to_string(),
            "1.1".to_string(),
            "".to_string(),
        ];
        assert_eq!(normalize_sample_values(&vals), "0/0:.:1.1:.");
    }

    #[test]
    fn test_format_samples_mapped_by_name() {
        let parsed = ParsedVcfRecord {
            chrom: "1",
            pos: 1,
            id: ".",
            ref_allele: "A",
            alt: "C",
            qual: ".",
            filter: ".",
            info: ".",
            format: Some(ParsedFormat {
                raw: "GT:FINT:FFLT:FSTR",
            }),
            samples: vec![
                ParsedSample {
                    raw: "0/0:11:1.1:AAA",
                },
                ParsedSample {
                    raw: "0/1:22:2.2:BBB",
                },
                ParsedSample {
                    raw: "0/0:33:3.3:CCC",
                },
            ],
        };

        let bundle = AnnotationBundle {
            alt: "C".to_string(),
            id: None,
            qual: None,
            filter: None,
            info: Vec::new(),
            format_str: Some("GT:FINT:FFLT:FSTR".to_string()),
            format_samples: vec![
                "1/1:88:8.8:BBB_DB".to_string(),
                "0/1:77:7.7:AAA_DB".to_string(),
            ],
            db_ref: "A".to_string(),
        };

        let sample_map = vec![Some(1), Some(0), None];

        let (fmt, samples) = merge_all_format(
            &parsed,
            &bundle,
            AnnotateMode::default_mode(),
            &sample_map,
            None,
            true,
            false,
        );

        assert_eq!(fmt, Some("GT:FINT:FFLT:FSTR".to_string()));
        assert_eq!(samples[0], "0/1:77:7.7:AAA_DB");
        assert_eq!(samples[1], "1/1:88:8.8:BBB_DB");
        assert_eq!(samples[2], "0/0:33:3.3:CCC");
    }

    #[test]
    fn test_alt_column_replaces_alt() {
        let parsed = ParsedVcfRecord {
            chrom: "1",
            pos: 10,
            id: ".",
            ref_allele: "A",
            alt: "C",
            qual: ".",
            filter: ".",
            info: ".",
            format: None,
            samples: Vec::new(),
        };
        let bundle = AnnotationBundle {
            alt: "G".to_string(),
            id: None,
            qual: None,
            filter: None,
            info: Vec::new(),
            format_str: None,
            format_samples: Vec::new(),
            db_ref: "A".to_string(),
        };
        let out = annotate_record_with_bundles(
            &parsed,
            &[(0, bundle)],
            &HashMap::new(),
            &[("ALT".to_string(), AnnotateMode::default_mode())],
            &[],
            false,
            false,
            false,
        );
        assert_eq!(out, "1\t10\t.\tA\tG\t.\t.\t.");
    }

    #[test]
    fn test_plus_alt_does_not_replace_existing_alt() {
        let (mode, key) = AnnotateMode::parse("+ALT");
        let parsed = ParsedVcfRecord {
            chrom: "1",
            pos: 10,
            id: ".",
            ref_allele: "A",
            alt: "C",
            qual: ".",
            filter: ".",
            info: ".",
            format: None,
            samples: Vec::new(),
        };
        let bundle = AnnotationBundle {
            alt: "G".to_string(),
            id: None,
            qual: None,
            filter: None,
            info: Vec::new(),
            format_str: None,
            format_samples: Vec::new(),
            db_ref: "A".to_string(),
        };
        let out = annotate_record_with_bundles(
            &parsed,
            &[(0, bundle)],
            &HashMap::new(),
            &[(key.to_string(), mode)],
            &[],
            false,
            false,
            false,
        );
        assert_eq!(out, "1\t10\t.\tA\tC\t.\t.\t.");
    }
