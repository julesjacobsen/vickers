use std::borrow::Cow;
use crate::variant_key::{decode_variant_key, encode_variant_key, Variant};

#[cfg(test)]
#[test]
fn test_encode_variant_key() {
    assert_eq!(0x0807728e88e80000, encode_variant_key("1", 976157, "T", "C"));
}

#[test]
fn test_decode_variant_key() {
    assert_eq!(
        decode_variant_key(&0x0806b567a0fee000),
        Variant {
            chrom: "1",
            pos: 879311,
            reference: Cow::Borrowed("TTTC"),
            alternate: Cow::Borrowed("T"),
        }
    );
    assert_eq!(
        decode_variant_key(&0x0807728e88e80000),
        Variant {
            chrom: "1",
            pos: 976157,
            reference: Cow::Borrowed("T"),
            alternate: Cow::Borrowed("C"),
        }
    );
}