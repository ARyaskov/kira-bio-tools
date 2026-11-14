use super::*;

#[test]
fn header_size_matches_constant() {
    assert_eq!(KtileHeader::SIZE, 144);
}

#[test]
fn header_validate_rejects_bad_magic() {
    let mut h = KtileHeader::zeroed();
    h.version = KTILE_VERSION;
    std::assert_matches!(h.validate(), Err(KtileError::BadMagic));
}

#[test]
fn header_validate_rejects_old_version() {
    let mut h = KtileHeader::zeroed();
    h.magic = KTILE_MAGIC;
    h.version = 0;
    std::assert_matches!(h.validate(), Err(KtileError::UnsupportedVersion { .. }));
}

#[test]
fn header_validate_accepts_current() {
    let mut h = KtileHeader::zeroed();
    h.magic = KTILE_MAGIC;
    h.version = KTILE_VERSION;
    assert!(h.validate().is_ok());
}

#[test]
fn header_validate_accepts_legacy_and_rejects_legacy_compact_offsets() {
    let mut header = KtileHeader::zeroed();
    header.magic = KTILE_MAGIC;
    header.version = KTILE_LEGACY_VERSION;
    assert!(header.validate().is_ok());
    header.flags = flags::HAS_ELIAS_FANO_OFFSETS;
    std::assert_matches!(
        header.validate(),
        Err(KtileError::Invalid {
            section: "line_offsets"
        })
    );
    header.version = KTILE_VERSION;
    assert!(header.validate().is_ok());
}
