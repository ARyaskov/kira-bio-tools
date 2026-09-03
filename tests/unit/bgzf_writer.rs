    use super::*;
    use flate2::read::MultiGzDecoder;
    use std::io::Read;

    fn tmp_path(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "kira_bgzf_{tag}_{}_{}.vcf.gz",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    const PAYLOAD: &[u8] = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tC\t.\t.\t.\n";

    fn decode(path: &std::path::Path) -> Vec<u8> {
        let mut decoded = Vec::new();
        MultiGzDecoder::new(File::open(path).unwrap()).read_to_end(&mut decoded).unwrap();
        decoded
    }

    #[test]
    fn test_bgzf_crc_decompresses_as_gzip() {
        let path = tmp_path("crc");
        {
            let mut writer = BgzfWriter::create(&path).unwrap();
            writer.write_all(&PAYLOAD[..8]).unwrap();
            writer.write_all(&PAYLOAD[8..]).unwrap();
            writer.finish().unwrap();
        }
        assert_eq!(decode(&path), PAYLOAD);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn drop_without_finish_still_writes_everything() {
        let path = tmp_path("drop");
        let big: Vec<u8> = (0..2_000_000u32).map(|i| b"ACGT\n"[(i % 5) as usize]).collect();
        {
            let mut writer = BgzfWriter::create(&path).unwrap();
            writer.write_all(&big).unwrap();
            // no finish(): Drop must finalize
        }
        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.ends_with(&BGZF_EOF), "missing BGZF EOF marker");
        assert_eq!(decode(&path), big);
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn writes_to_arbitrary_writer() {
        let sink = std::sync::Arc::new(std::sync::Mutex::new(Vec::<u8>::new()));
        struct Shared(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl Write for Shared {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let mut w = BgzfWriter::from_writer_buffered(Shared(sink.clone()), Compression::new(3), 1 << 16).unwrap();
        w.write_all(PAYLOAD).unwrap();
        w.finish().unwrap();
        let bytes = sink.lock().unwrap().clone();
        assert!(bytes.ends_with(&BGZF_EOF));
        let mut decoded = Vec::new();
        MultiGzDecoder::new(&bytes[..]).read_to_end(&mut decoded).unwrap();
        assert_eq!(decoded, PAYLOAD);
    }
