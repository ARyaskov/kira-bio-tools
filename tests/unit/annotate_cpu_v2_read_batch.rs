    use super::*;

    #[test]
    fn push_and_read_roundtrip() {
        let mut b = ReadBatch::with_capacity(64, 8);
        b.push_line("1\t100\trs1\tA\tT\t.\t.\t.");
        b.push_line("1\t200\trs2\tC\tG\t.\t.\t.\n");
        b.push_line("2\t1\t.\tA\tAT\t.\t.\t.\r\n");
        assert_eq!(b.len(), 3);
        assert_eq!(b.line(0), "1\t100\trs1\tA\tT\t.\t.\t.");
        assert_eq!(b.line(1), "1\t200\trs2\tC\tG\t.\t.\t.");
        assert_eq!(b.line(2), "2\t1\t.\tA\tAT\t.\t.\t.");
    }

    #[test]
    fn empty_batch_iter() {
        let b = ReadBatch::with_capacity(0, 0);
        assert_eq!(b.iter().count(), 0);
        assert_eq!(b.line(0), "");
    }

    #[test]
    fn iter_visits_in_order() {
        let mut b = ReadBatch::with_capacity(64, 4);
        b.push_line("a");
        b.push_line("bb");
        b.push_line("ccc");
        let lines: Vec<&str> = b.iter().collect();
        assert_eq!(lines, vec!["a", "bb", "ccc"]);
    }
