#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parallel_read_write_roundtrip() {
        let data = b"##fileformat=VCFv4.3\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                     chr1\t1000\t.\tA\tT\t30\tPASS\tDP=10\n";

        let mut tmp = NamedTempFile::new().unwrap();
        let path = tmp.path().to_owned();

        // Write
        {
            let mut writer = ParallelBgzfWriter::create(&path).unwrap();
            writer.write_all(data).unwrap();
            writer.finish().unwrap();
        }

        // Read
        let mut reader = ParallelBgzfReader::open(&path).unwrap();
        let blocks = reader.read_batch().unwrap();

        let mut output = Vec::new();
        for block in blocks {
            output.extend_from_slice(&block.uncompressed);
        }

        assert_eq!(&output[..], data);
    }

    #[test]
    fn test_batched_line_reader() {
        // Create test BGZF file
        let mut tmp = NamedTempFile::new().unwrap();
        let path = tmp.path().to_owned();

        {
            let mut writer = ParallelBgzfWriter::create(&path).unwrap();
            for i in 0..1000 {
                writeln!(writer, "chr1\t{}\t.\tA\tT\t30\tPASS\tDP={}", i, i).unwrap();
            }
            writer.finish().unwrap();
        }

        let reader = ParallelBgzfReader::open(&path).unwrap();
        let mut line_reader = BatchedLineReader::new(reader, 100);

        let batch = line_reader.read_batch().unwrap();
        assert!(!batch.is_empty());
        assert!(batch[0].0.starts_with("chr1"));
    }
}
