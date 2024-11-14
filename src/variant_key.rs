use crate::hash::encode_refalt_hash;
use lazy_static::lazy_static;
use regex::Regex;
use std::borrow::Cow;
use std::fmt::{Display, Formatter};

lazy_static! {
    static ref VARIANT_REGEX: Regex = Regex::new(r"(?<chrom>1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y|M|MT)[-:](?<pos>[0-9]+)[-:](?<ref>[ACGT]+)[-:](?<alt>[ACGT]+)").unwrap();
}

#[derive(Debug, PartialOrd, PartialEq)]
pub struct Variant<'a> {
    pub chrom: &'a str,
    pub pos: u32,
    pub reference: Cow<'a, str>,
    pub alternate: Cow<'a, str>,
}

impl Display for Variant<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.pad(&format!("{}-{}-{}-{}", self.chrom, self.pos, self.reference, self.alternate))
    }
}

impl Variant<'_> {
    pub fn to_variant_key(&self) -> VariantKey {
        encode_variant_key(self.chrom, self.pos, self.reference.as_ref(), self.alternate.as_ref())
    }
}

pub type VariantKey = u64;

pub fn parse_variant(variant_string: &str) -> Option<Variant> {
    VARIANT_REGEX.captures(&variant_string).map(|caps| {
        let chrom = caps.name("chrom").unwrap().as_str();
        let pos = caps.name("pos").unwrap().as_str().parse::<u32>().unwrap();
        let reference = Cow::Borrowed(caps.name("ref").unwrap().as_str());
        let alternate = Cow::Borrowed(caps.name("alt").unwrap().as_str());
        Variant { chrom, pos, reference, alternate }
    })
}


const VKMASK_CHROM: u64 = 0xF800000000000000;  // VariantKey binary mask for CHROM     [ 11111000 00000000 00000000 00000000 00000000 00000000 00000000 00000000 ]
const VKMASK_POS: u64 = 0x07FFFFFF80000000;  // VariantKey binary mask for POS       [ 00000111 11111111 11111111 11111111 10000000 00000000 00000000 00000000 ]
const VKSHIFT_CHROM: u64 = 59; // CHROM LSB position from the VariantKey LSB
const VKSHIFT_POS: u64 = 31; //  POS LSB position from the VariantKey LSB
const ALLELE_LEN_MASK: u64 = 0b1111;
const BASE_MASK: u64 = 0b11;


pub fn encode_variant_key(chrom: &str, pos: u32, reference: &str, alternate: &str) -> VariantKey {
    if pos > (1 << 28) {
        panic!("Position {} too large for variant key! Max size is {}", pos, 1 << 28);
    }
    // if (reference.len() + alternate.len()) > 11 {
    //     panic!("Total length of ref and alt alleles must be max 11 bases");
    // }
    encode_chrom(chrom) << VKSHIFT_CHROM | (pos as u64) << VKSHIFT_POS | encode_ref_alt(reference, alternate)
}

fn encode_chrom(chrom: &str) -> u64 {
    match chrom {
        n if chrom.parse::<u64>().is_ok() => n.parse::<u64>().unwrap(),
        "X" => 23,
        "Y" => 24,
        "M" => 25,
        "MT" => 25,
        _ => 0,
        // _ => panic!("Illegal chromosome string '{}'. Should be one of 1-22,X,Y,M.", chrom),
    }
}

pub fn encode_ref_alt(reference: &str, alternate: &str) -> u64 {
    let ref_alt = if reference.len() + alternate.len() <= 11 && is_just_acgt(reference) && is_just_acgt(alternate) {
        encode_ref_alt_rev(reference, alternate)
    } else {
        encode_refalt_hash(reference.as_ref(), alternate.as_ref()).into()
    };
    ref_alt
}

fn is_just_acgt(allele: &str) -> bool {
    for c in allele.chars() {
        if c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'a' && c != 'c' && c != 'g' && c != 't' {
            return false;
        }
    }
    true
}

fn encode_ref_alt_rev(reference: &str, alternate: &str) -> u64 {
    let mut bits: u64 = 0;

    let ref_length = reference.len() as u64;
    bits |= ref_length << 27;

    let alt_length = alternate.len() as u64;
    bits |= alt_length << 23;

    let ref_encoded = encode_allele(reference);
    bits |= ref_encoded << (23 - (ref_length << 1));

    let alt_encoded = encode_allele(alternate);
    bits |= alt_encoded << (23 - ((ref_length + alt_length) << 1));

    bits
}

fn encode_allele(allele: &str) -> u64 {
    let mut encoded: u64 = 0;
    let allele_length = allele.len();
    for (i, c) in allele.chars().enumerate() {
        let enc = encode_base(c);
        encoded |= (enc as u64) << ((allele_length - 1 - i) << 1);
    }
    encoded
}

fn encode_base(c: char) -> u8 {
    match c {
        'A' | 'a' => 0b00,
        'T' | 't' => 0b11,
        'C' | 'c' => 0b01,
        'G' | 'g' => 0b10,
        _ => panic!("Unable to encode base {} must be [AaTtGgCc]", c),
    }
}

pub fn decode_variant_key(variant_key: &VariantKey) -> Variant {
    let chrom = decode_chrom(variant_key);
    let pos = decode_pos(variant_key);
    let reference = decode_ref(variant_key);
    let alternate = decode_alt(variant_key);
    Variant { chrom, pos, reference, alternate }
}

const BASE_CHARS: [char; 4] = ['A', 'C', 'G', 'T'];

fn decode_base(base: u8) -> char {
    BASE_CHARS[base as usize]
}

const BASE_STR: [&str; 4] = ["A", "C", "G", "T"];

fn decode_base_as_str(base: u64) -> &'static str {
    BASE_STR[base as usize]
}

const CHROMS: [&str; 26] = ["NA",
    "1", "2", "3", "4", "5",
    "6", "7", "8", "9", "10",
    "11", "12", "13", "14", "15",
    "16", "17", "18", "19", "20",
    "21", "22", "X", "Y", "MT"];

fn decode_chrom<'a>(variant_key: &VariantKey) -> &'static str {
    let chr_bits = ((variant_key & VKMASK_CHROM) >> VKSHIFT_CHROM) as usize;
    match chr_bits {
        chr @ 1..=25 => CHROMS[chr],
        _ => CHROMS[0]
    }
}

fn decode_pos(variant_key: &VariantKey) -> u32 {
    ((variant_key & VKMASK_POS) >> VKSHIFT_POS) as u32
}

fn ref_length(variant_key: &VariantKey) -> u64 {
    (variant_key >> 27) & ALLELE_LEN_MASK
}

fn alt_length(variant_key: &VariantKey) -> u64 {
    (variant_key >> 23) & ALLELE_LEN_MASK
}

fn decode_ref(variant_key: &VariantKey) -> Cow<'_, str> {
    let ref_length = ref_length(variant_key);
    decode_allele(ref_length, 21, variant_key)
}

fn decode_alt(variant_key: &VariantKey) -> Cow<'_, str> {
    let alt_length = alt_length(variant_key);
    decode_allele(alt_length, 21 - (ref_length(variant_key) << 1), variant_key)
}

fn decode_allele(allele_length: u64, allele_offset: u64, variant_key: &VariantKey) -> Cow<'_, str> {
    match allele_length {
        0 => Cow::Borrowed(""),
        1 => Cow::Borrowed(decode_base_as_str((variant_key >> allele_offset) & BASE_MASK)),
        2..=10 => Cow::Owned(new_allele_string(allele_length, allele_offset, variant_key)),
        _ => Cow::Borrowed("?"),
    }
}

fn new_allele_string(allele_length: u64, allele_offset: u64, variant_key: &VariantKey) -> String {
    let mut bases = String::with_capacity(allele_length as usize);
    for i in 0..allele_length {
        let base = decode_base(((variant_key >> (allele_offset - (i << 1))) & BASE_MASK) as u8);
        bases.push(base);
    }
    bases
}