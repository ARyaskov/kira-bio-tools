    use super::*;

    #[test]
    fn header_size_matches_constant() {
        assert_eq!(KtileHeader::SIZE, 144);
    }

    #[test]
    fn header_validate_rejects_bad_magic() {
        let mut h = KtileHeader::zeroed();
        h.version = KTILE_VERSION;
        assert!(matches!(h.validate(), Err(KtileError::BadMagic)));
    }

    #[test]
    fn header_validate_rejects_old_version() {
        let mut h = KtileHeader::zeroed();
        h.magic = KTILE_MAGIC;
        h.version = 0;
        assert!(matches!(
            h.validate(),
            Err(KtileError::UnsupportedVersion { .. })
        ));
    }

    #[test]
    fn header_validate_accepts_current() {
        let mut h = KtileHeader::zeroed();
        h.magic = KTILE_MAGIC;
        h.version = KTILE_VERSION;
        assert!(h.validate().is_ok());
    }
