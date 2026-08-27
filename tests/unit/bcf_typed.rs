    use super::*;
    use std::io::Cursor;

    #[test]
    fn roundtrip_typed_ints_small() {
        let vals = vec![1i32, 2, 3];
        let mut buf = Vec::new();
        write_typed_ints(&mut buf, &vals).unwrap();
        let mut c = Cursor::new(buf);
        match read_typed(&mut c).unwrap() {
            BcfValue::Ints(v) => assert_eq!(v, vals),
            _ => panic!(),
        }
    }

    #[test]
    fn roundtrip_typed_ints_large() {
        let vals: Vec<i32> = (1..=20).collect();
        let mut buf = Vec::new();
        write_typed_ints(&mut buf, &vals).unwrap();
        let mut c = Cursor::new(buf);
        match read_typed(&mut c).unwrap() {
            BcfValue::Ints(v) => assert_eq!(v, vals),
            _ => panic!(),
        }
    }

    #[test]
    fn roundtrip_int16_range() {
        let vals = vec![0i32, 1000, -500, 32000];
        let mut buf = Vec::new();
        write_typed_ints(&mut buf, &vals).unwrap();
        let mut c = Cursor::new(buf);
        match read_typed(&mut c).unwrap() {
            BcfValue::Ints(v) => assert_eq!(v, vals),
            _ => panic!(),
        }
    }

    #[test]
    fn roundtrip_floats() {
        let vals = vec![1.5f32, 2.25, -3.125];
        let mut buf = Vec::new();
        write_typed_floats(&mut buf, &vals).unwrap();
        let mut c = Cursor::new(buf);
        match read_typed(&mut c).unwrap() {
            BcfValue::Floats(v) => assert_eq!(v, vals),
            _ => panic!(),
        }
    }

    #[test]
    fn roundtrip_string() {
        let s = b"hello world";
        let mut buf = Vec::new();
        write_typed_string(&mut buf, s).unwrap();
        let mut c = Cursor::new(buf);
        match read_typed(&mut c).unwrap() {
            BcfValue::Str(v) => assert_eq!(v, s.to_vec()),
            _ => panic!(),
        }
    }

    #[test]
    fn gt_encode_decode() {
        assert_eq!(encode_gt(Some(0), false), 2);
        assert_eq!(encode_gt(Some(1), true), 5);
        assert_eq!(encode_gt(None, false), 0);
        assert_eq!(decode_gt(2), (Some(0), false));
        assert_eq!(decode_gt(5), (Some(1), true));
        assert_eq!(decode_gt(0), (None, false));
    }
