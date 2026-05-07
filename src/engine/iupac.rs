//! IUPAC nucleotide-code utilities.
//!
//! The four canonical bases are encoded as bits in a 4-bit mask:
//! A=0001, C=0010, G=0100, T=1000. IUPAC ambiguity codes are unions of these
//! (e.g. R = A|G = 0101, N = 1111). The bitmask form makes consensus building
//! and consensus-vs-base matching one bitwise op each.

pub const MASK_A: u8 = 0b0001;
pub const MASK_C: u8 = 0b0010;
pub const MASK_G: u8 = 0b0100;
pub const MASK_T: u8 = 0b1000;

#[inline]
pub const fn base_mask(b: u8) -> u8 {
    match b {
        b'A' | b'a' => MASK_A,
        b'C' | b'c' => MASK_C,
        b'G' | b'g' => MASK_G,
        b'T' | b't' | b'U' | b'u' => MASK_T,
        b'R' | b'r' => MASK_A | MASK_G,
        b'Y' | b'y' => MASK_C | MASK_T,
        b'S' | b's' => MASK_C | MASK_G,
        b'W' | b'w' => MASK_A | MASK_T,
        b'K' | b'k' => MASK_G | MASK_T,
        b'M' | b'm' => MASK_A | MASK_C,
        b'B' | b'b' => MASK_C | MASK_G | MASK_T,
        b'D' | b'd' => MASK_A | MASK_G | MASK_T,
        b'H' | b'h' => MASK_A | MASK_C | MASK_T,
        b'V' | b'v' => MASK_A | MASK_C | MASK_G,
        b'N' | b'n' => 0b1111,
        _ => 0,
    }
}

#[inline]
pub const fn mask_to_iupac(mask: u8) -> u8 {
    match mask {
        MASK_A => b'A',
        MASK_C => b'C',
        MASK_G => b'G',
        MASK_T => b'T',
        0b0101 => b'R',
        0b1010 => b'Y',
        0b0110 => b'S',
        0b1001 => b'W',
        0b1100 => b'K',
        0b0011 => b'M',
        0b1110 => b'B',
        0b1101 => b'D',
        0b1011 => b'H',
        0b0111 => b'V',
        0b1111 => b'N',
        _ => b'N',
    }
}

#[inline]
pub fn is_ambiguous(b: u8) -> bool {
    matches!(
        b,
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N'
    )
}

#[inline]
pub const fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        other => other,
    }
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}
