    use super::*;
    use crate::bcf::header::parse_header_to_dict;

    fn dict_full() -> BcfHeaderDict {
        let h = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">".to_string(),
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">".to_string(),
            "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"GP\">".to_string(),
            "##contig=<ID=1>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2".to_string(),
        ];
        parse_header_to_dict(&h)
    }

    fn roundtrip(line: &str) -> String {
        let dict = dict_full();
        let mut buf = Vec::new();
        encode_record(&mut buf, line, &dict).unwrap();
        let mut c = Cursor::new(buf);
        decode_record_to_vcf(&mut c, &dict).unwrap().unwrap()
    }

    #[test]
    fn record_roundtrip_basic() {
        let out = roundtrip("1\t100\t.\tA\tT\t60\tPASS\tDP=10\tGT\t0/1\t1|1");
        assert_eq!(out, "1\t100\t.\tA\tT\t60\tPASS\tDP=10\tGT\t0/1\t1|1");
    }

    #[test]
    fn phasing_and_missing_alleles_survive() {
        let out = roundtrip("1\t100\t.\tA\tT\t.\t.\t.\tGT\t0|1\t.|.");
        assert!(out.ends_with("\tGT\t0|1\t.|."), "got {out:?}");
        let out = roundtrip("1\t100\t.\tA\tT\t.\t.\t.\tGT\t./.\t1");
        assert!(out.ends_with("\tGT\t./.\t1"), "got {out:?}");
    }

    #[test]
    fn ragged_float_vectors_use_vector_end() {
        let out = roundtrip("1\t100\t.\tA\tT,G\t.\t.\t.\tGT:GP\t0/1:0.1,0.2,0.7\t0/2:0.5,0.5");
        assert!(out.ends_with("\tGT:GP\t0/1:0.1,0.2,0.7\t0/2:0.5,0.5"), "got {out:?}");
        let out = roundtrip("1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:3,4\t0/0:.");
        assert!(out.ends_with("\tGT:AD\t0/1:3,4\t0/0:."), "got {out:?}");
    }

    #[test]
    fn symbolic_end_sets_rlen() {
        let dict = dict_full();
        let mut buf = Vec::new();
        encode_record(&mut buf, "1\t100\t.\tA\t<DEL>\t.\t.\tEND=200", &dict).unwrap();
        let (shared, _) = read_record_raw(&mut Cursor::new(buf)).unwrap().unwrap();
        let m = record_meta(&shared).unwrap();
        assert_eq!(m.pos, 99);
        assert_eq!(m.rlen, 101);
        assert_eq!(record_end0(&m), 200);
    }

    #[test]
    fn corrupt_length_is_rejected() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&u32::MAX.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        assert!(read_record_raw(&mut Cursor::new(buf)).is_err());
    }

    /// Decompressed BCF written by bcftools 1.19 (shared IDX dictionary, five records).
    const BCFTOOLS_BCF_HEX: &[&str] = &[
        "424346020201030000232366696c65666f726d61743d56434676342e300a232346494c5445523d3c49443d504153532c",
        "4465736372697074696f6e3d22416c6c2066696c7465727320706173736564222c4944583d303e0a2323494e464f3d3c",
        "49443d58582c4e756d6265723d312c547970653d496e74656765722c4465736372697074696f6e3d2254657374222c49",
        "44583d313e0a2323494e464f3d3c49443d44502c4e756d6265723d312c547970653d496e74656765722c446573637269",
        "7074696f6e3d22546f74616c204465707468222c4944583d323e0a2323464f524d41543d3c49443d47512c4e756d6265",
        "723d312c547970653d496e74656765722c4465736372697074696f6e3d2247656e6f74797065205175616c697479222c",
        "4944583d333e0a2323464f524d41543d3c49443d44502c4e756d6265723d312c547970653d496e74656765722c446573",
        "6372697074696f6e3d2252656164204465707468222c4944583d323e0a2323464f524d41543d3c49443d47542c4e756d",
        "6265723d312c547970653d537472696e672c4465736372697074696f6e3d2247656e6f74797065222c4944583d343e0a",
        "232346494c5445523d3c49443d4661696c2c4465736372697074696f6e3d224661696c222c4944583d353e0a2323636f",
        "6e7469673d3c49443d312c6c656e6774683d36323433353936342c4944583d303e0a2323636f6e7469673d3c49443d32",
        "2c6c656e6774683d36323433353936342c4944583d313e0a2323626366746f6f6c735f7669657756657273696f6e3d31",
        "2e31392b6874736c69622d312e31390a2323626366746f6f6c735f76696577436f6d6d616e643d76696577202d4f6220",
        "2d6f2074657374732f636f6e6361742f74657374352f636f6e6361742e322e612e6263662074657374732f636f6e6361",
        "742f74657374352f636f6e6361742e322e612e7663662e677a3b20446174653d547565204665622032342030363a3034",
        "3a343120323032360a234348524f4d09504f530949440952454609414c54095155414c0946494c54455209494e464f09",
        "464f524d415409410a00230000000e000000010000008b0000000100000000c035440100020001000003071741174711",
        "001102111e110421020411031296001102111e2d0000000d000000010000009f00000005000000000076430100040001",
        "00000307575441414141275441275443175411001102110a11042102061103110c1102110a290000000e000000000000",
        "006d000000010000000000e04402000300010000030717431754174711051101110b110211201104210204110312f500",
        "11021120260000000e00000000000000810000000300000000007e440100020001000003073747414127474711001102",
        "11161104210204110312d40011021116230000000e00000000000000810000000100000000007e440100020001000003",
        "07174717541100110211161104210204110312d40011021116",
    ];

    fn bcftools_bcf_bytes() -> Vec<u8> {
        let joined: String = BCFTOOLS_BCF_HEX.concat();
        (0..joined.len()).step_by(2).map(|i| u8::from_str_radix(&joined[i..i + 2], 16).unwrap()).collect()
    }

    #[test]
    fn decodes_bcftools_written_bcf_stream() {
        let bytes = bcftools_bcf_bytes();
        let mut r = crate::bcf::BcfReader::from_bufread(Box::new(std::io::Cursor::new(bytes))).unwrap();
        assert!(r.header_lines.iter().any(|l| l.starts_with("#CHROM")));
        let mut n = 0;
        loop {
            match r.read_record_line() {
                Ok(Some(line)) => {
                    n += 1;
                    if n == 1 {
                        assert_eq!(line, "2	140	.	A	G	727	PASS	DP=30	GT:GQ:DP	0/1:150:30");
                    }
                }
                Ok(None) => break,
                Err(e) => panic!("record {}: {:#}", n + 1, e),
            }
        }
        assert_eq!(n, 5);
    }
