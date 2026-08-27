    use super::*;
    use flate2::read::MultiGzDecoder;
    use std::io::Read;

    #[test]
    fn test_bgzf_crc_decompresses_as_gzip() {
        let path = std::env::temp_dir().join(format!(
            "kira_bgzf_crc_{}_{}.vcf.gz",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let payload = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tC\t.\t.\t.\n";
        {
            let mut writer = BgzfWriter::create(&path).unwrap();
            writer.write_all(&payload[..8]).unwrap();
            writer.write_all(&payload[8..]).unwrap();
            writer.finish().unwrap();
        }
        let mut decoded = Vec::new();
        let file = File::open(&path).unwrap();
        MultiGzDecoder::new(file).read_to_end(&mut decoded).unwrap();
        assert_eq!(decoded, payload);
        let _ = std::fs::remove_file(path);
    }
