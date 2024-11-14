// based on https://github.com/Genomicsplc/variantkey/blob/master/c/src/variantkey/variantkey.h
// LICENSE
//
// Copyright (c) 2017-2018 GENOMICS plc
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

use std::u32;

// hash
fn muxhash(mut k: u32, mut h: u32) -> u32 {
    // println!("k={}, h={}", k, h);
    k = k.wrapping_mul(0xcc9e2d51);
    k = (k >> 17) | (k << 15);
    k = k.wrapping_mul(0x1b873593);
    // println!("K={}", k);
    h ^= k;
    h = (h >> 19) | (h << 13);
    let mux = h.wrapping_mul(5).wrapping_add(0xe6546b64);
    // println!("mux={}", mux);
    mux
}

fn encode_packchar(c: u8) -> u32 {
    if c < b'A' {
        return 27;
    }
    if c >= b'a' {
        return (c - b'a' + 1).into();
    }
    (c - b'A' + 1).into()
}

// pack blocks of 6 characters in 32 bit (6 x 5 bit + 2 spare bit) [ 01111122 22233333 44444555 55666660 ]
fn pack_chars_tail(value: &[u8]) -> u32 {
    let mut h: u32 = 0;

    if value.len() >= 5 {
        h ^= encode_packchar(value[4]) << (1 + 5);
    }
    if value.len() >= 4 {
        h ^= encode_packchar(value[3]) << (1 + (5 * 2));
    }
    if value.len() >= 3 {
        h ^= encode_packchar(value[2]) << (1 + (5 * 3));
    }
    if value.len() >= 2 {
        h ^= encode_packchar(value[1]) << (1 + (5 * 4));
    }
    if !value.is_empty() {
        h ^= encode_packchar(value[0]) << (1 + (5 * 5));
    }
    // println!("packed tail={:x}", h);
    h
}

fn pack_chars(value: &[u8]) -> u32 {
    let packed = (encode_packchar(value[5]) << 1)
        ^ (encode_packchar(value[4]) << (1 + 5))
        ^ (encode_packchar(value[3]) << (1 + (5 * 2)))
        ^ (encode_packchar(value[2]) << (1 + (5 * 3)))
        ^ (encode_packchar(value[1]) << (1 + (5 * 4)))
        ^ (encode_packchar(value[0]) << (1 + (5 * 5)));
    // println!("packed={}", packed);
    packed
}

// Return a 32 bit hash of a nucleotide string
fn hash32(mut value: &[u8]) -> u32 {
    let mut h: u32 = 0;
    let len = 6;
    while value.len() >= len {
        h = muxhash(pack_chars(value), h);
        value = &value[len..];
        // println!("{:?}", String::from_utf8_lossy(value))
    }
    if !value.is_empty() {
        h = muxhash(pack_chars_tail(value), h);
    }
    // println!("h={}", h);
    h
}

pub fn encode_refalt_hash(reference: &str, alternative: &str) -> u32 {
    // 0x3 is the separator character between REF and ALT [00000000 00000000 00000000 00000011]
    let mut h = muxhash(hash32(alternative.as_ref()), muxhash(0x3, hash32(reference.as_ref())));
    // MurmurHash3 finalization mix - force all bits of a hash block to avalanche
    h ^= h >> 16;
    h = h.wrapping_mul(0x85ebca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2ae35);
    h ^= h >> 16;
    h >> 1 | 0x1 // 0x1 is the set bit to indicate HASH mode [00000000 00000000 00000000 00000001]
}

#[test]
fn test_encode_variant_key() {
    println!("{:x}", encode_refalt_hash("T", "C"));
    println!("{:x}", encode_refalt_hash("GCCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCCCCCTCATCA", "G"));
    println!("{:x}", encode_refalt_hash("AACCATCTGTATTGATGCACTGTCCATGTTT", "A"));
    // assert_eq!(0x0807728e88e80000, encode_refalt_hash("T", "C"));
    // assert_eq!(0x0806ac8d5179a93f, encode_refalt_hash("GCCTCCCCAGCCACGGTGAGGACCCACCCTGGCATGATCCCCCTCATCA", "G"));
    // assert_eq!(0x7973c63e408cbc69, encode_refalt_hash("AACCATCTGTATTGATGCACTGTCCATGTTT", "A"));
}