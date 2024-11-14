use std::borrow::Cow;
use crate::variant_key::{decode_variant_key, encode_variant_key, Variant, encode_ref_alt};

#[cfg(test)]
#[test]
fn test_encode_variant_key() {
    assert_eq!(0x0807728e88e80000, encode_variant_key("1", 976157, "T", "C"));
    assert_eq!(0x0806ac8d5179a93f, encode_variant_key("1", 874778, "GCCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCCCCCTCATCA", "G"));
    assert_eq!(0x7973c63e408cbc69, encode_variant_key("15", 48729212, "AACCATCTGTATTGATGCACTGTCCATGTTT", "A"));
    assert_eq!(0xc003670d08900000, encode_variant_key("Y", 445978, "A", "G"));
    assert_eq!(0xc003670d6842610d, encode_variant_key("Y", 445978, "AAAGAAAGAAAGAAAGAAAG", "A"));
}

#[test]
fn test_encode_ref_alt() {
    println!("{:x}", encode_ref_alt("T", "C"));
    println!("{:x}", encode_ref_alt("ACG", "ACGTACGTAC"));
    println!("{:x}", encode_ref_alt("GCCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCCCCCTCATCA", "G"));
    println!("{:x}", encode_ref_alt("AACCATCTGTATTGATGCACTGTCCATGTTT", "A"));
    println!("{:x}", encode_ref_alt("A", "ACGTacgtACGT"));
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