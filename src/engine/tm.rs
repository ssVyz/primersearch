//! Nearest-neighbor thermodynamic Tm calculation.
//!
//! Uses the SantaLucia (1998) unified nearest-neighbor parameters for ΔH
//! and ΔS, with two-end initiation (G/C vs A/T) and a monovalent-equivalent
//! salt correction that folds in Mg2+ via [Na⁺]_eq = [Na⁺] + 120·√(Mg_free)
//! (von Ahsen 2001 / Owczarzy 2008 method 4), where Mg_free = Mg − dNTP
//! (one Mg sequestered per dNTP). For oligo sequences containing IUPAC
//! ambiguity codes, every A/C/G/T expansion is enumerated and the median
//! Tm is returned.
//!
//! Concentration convention for the CT/x denominator:
//! the user-supplied oligo concentration is interpreted as a single
//! strand's concentration in a non-self-complementary duplex with both
//! strands at the same concentration. Total CT = 2·oligo_conc; x = 4;
//! so effectively the divisor inside the log is oligo_conc / 2.
//!
//! References:
//!   - SantaLucia J. (1998) PNAS 95(4): 1460–1465.
//!   - von Ahsen N. et al. (2001) Clin Chem 47(11): 1956–1961.
//!   - Owczarzy R. et al. (2008) Biochemistry 47(19): 5336–5353.

use crate::engine::iupac::base_mask;

/// Thermodynamic conditions for the Tm calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TmParams {
    pub oligo_conc_um: f64,
    pub na_conc_mm: f64,
    pub mg_conc_mm: f64,
    pub dntp_conc_mm: f64,
}

impl TmParams {
    /// Equivalent monovalent concentration in M.
    pub fn na_eq_molar(&self) -> f64 {
        let mg_free_mm = (self.mg_conc_mm - self.dntp_conc_mm).max(0.0);
        let na_eq_mm = self.na_conc_mm + 120.0 * mg_free_mm.sqrt();
        (na_eq_mm * 1e-3).max(1e-10)
    }

    /// Denominator term inside `ln(...)` of the Tm formula, in M.
    /// `2·oligo_conc / 4 = oligo_conc / 2` — see module docs.
    pub fn ct_eff_molar(&self) -> f64 {
        ((self.oligo_conc_um * 1e-6) * 0.5).max(1e-30)
    }
}

/// Universal gas constant, cal/(K·mol).
const R_GAS: f64 = 1.987;

// SantaLucia 1998 unified NN parameters. Indexed by `(top << 2) | bot`
// where A=0, C=1, G=2, T=3. Top is the 5' base of the dinucleotide, bot
// the 3' base. Order: AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT.
#[rustfmt::skip]
const NN_DH: [f64; 16] = [
    -7.9, -8.4, -7.8, -7.2, // AA AC AG AT
    -8.5, -8.0, -10.6, -7.8, // CA CC CG CT
    -8.2, -9.8, -8.0, -8.4, // GA GC GG GT
    -7.2, -8.2, -8.5, -7.9, // TA TC TG TT
];
#[rustfmt::skip]
const NN_DS: [f64; 16] = [
    -22.2, -22.4, -21.0, -20.4, // AA AC AG AT
    -22.7, -19.9, -27.2, -21.0, // CA CC CG CT
    -22.2, -24.4, -19.9, -22.4, // GA GC GG GT
    -21.3, -22.2, -22.7, -22.2, // TA TC TG TT
];

/// Initiation contribution (ΔH in kcal/mol, ΔS in cal/(K·mol)) for one
/// end of the duplex. SantaLucia 1998 distinguishes G/C-terminated from
/// A/T-terminated ends.
#[inline]
fn init_end(b: u8) -> (f64, f64) {
    match b {
        b'C' | b'c' | b'G' | b'g' => (0.1, -2.8),
        _ => (2.3, 4.1),
    }
}

#[inline]
fn base_to_idx(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Tm in °C for an exact A/C/G/T sequence of length ≥ 2. Returns `f64::NAN`
/// for inputs that are too short or contain non-ACGT characters.
fn tm_acgt(seq: &[u8], params: &TmParams) -> f64 {
    if seq.len() < 2 {
        return f64::NAN;
    }
    let mut dh = 0.0_f64;
    let mut ds = 0.0_f64;
    let Some(mut prev) = base_to_idx(seq[0]) else {
        return f64::NAN;
    };
    for &b in &seq[1..] {
        let Some(cur) = base_to_idx(b) else {
            return f64::NAN;
        };
        let k = (prev << 2) | cur;
        dh += NN_DH[k];
        ds += NN_DS[k];
        prev = cur;
    }
    let (h5, s5) = init_end(seq[0]);
    let (h3, s3) = init_end(*seq.last().expect("len ≥ 2"));
    dh += h5 + h3;
    ds += s5 + s3;

    let n = seq.len() as f64;
    let ds_salt = ds + 0.368 * (n - 1.0) * params.na_eq_molar().ln();
    let denom = ds_salt + R_GAS * params.ct_eff_molar().ln();
    (dh * 1000.0) / denom - 273.15
}

/// Tm in °C for a sequence that may contain IUPAC ambiguity codes. The
/// median over all A/C/G/T expansions is returned. For pure A/C/G/T input
/// this is identical to `tm_acgt`.
pub fn calculate_tm(sequence: &[u8], params: &TmParams) -> f64 {
    if sequence.len() < 2 {
        return 0.0;
    }
    let mut amb_positions: Vec<usize> = Vec::new();
    for (i, &b) in sequence.iter().enumerate() {
        let m = base_mask(b);
        if m == 0 {
            return 0.0;
        }
        if m.count_ones() > 1 {
            amb_positions.push(i);
        }
    }
    if amb_positions.is_empty() {
        let tm = tm_acgt(sequence, params);
        return if tm.is_nan() { 0.0 } else { tm };
    }
    let mut variant: Vec<u8> = sequence.to_vec();
    let mut tms: Vec<f64> = Vec::new();
    enumerate_variants(&mut variant, &amb_positions, 0, params, &mut tms);
    if tms.is_empty() {
        return 0.0;
    }
    tms.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = tms.len();
    if n % 2 == 1 {
        tms[n / 2]
    } else {
        (tms[n / 2 - 1] + tms[n / 2]) / 2.0
    }
}

fn enumerate_variants(
    variant: &mut Vec<u8>,
    positions: &[usize],
    depth: usize,
    params: &TmParams,
    tms: &mut Vec<f64>,
) {
    if depth == positions.len() {
        let tm = tm_acgt(variant, params);
        if !tm.is_nan() {
            tms.push(tm);
        }
        return;
    }
    let pos = positions[depth];
    let original = variant[pos];
    let mask = base_mask(original);
    const BASES: [(u8, u8); 4] = [
        (0b0001, b'A'),
        (0b0010, b'C'),
        (0b0100, b'G'),
        (0b1000, b'T'),
    ];
    for (bit, base) in BASES {
        if mask & bit != 0 {
            variant[pos] = base;
            enumerate_variants(variant, positions, depth + 1, params, tms);
        }
    }
    variant[pos] = original;
}

/// Smallest length `L ≥ 2` of the substring starting at `start` whose Tm
/// reaches `tm_threshold`, or `None` if extending to the end of `seq`
/// does not. Computed incrementally — O(L) per call. Input is assumed
/// A/C/G/T (the engine's quality filter rejects sequences with anything
/// else); a non-ACGT base aborts and returns `None`.
pub fn determine_oligo_length(
    seq: &[u8],
    start: usize,
    tm_threshold: f64,
    params: &TmParams,
) -> Option<usize> {
    if start >= seq.len() {
        return None;
    }
    let max_len = seq.len() - start;
    if max_len < 2 {
        return None;
    }
    let salt_ln = params.na_eq_molar().ln();
    let ct_term = R_GAS * params.ct_eff_molar().ln();

    let first_b = seq[start];
    let first_idx = base_to_idx(first_b)?;
    let (h5, s5) = init_end(first_b);
    let mut dh_nn = 0.0_f64;
    let mut ds_nn = 0.0_f64;
    let mut prev = first_idx;

    for length in 2..=max_len {
        let cur_b = seq[start + length - 1];
        let cur_idx = base_to_idx(cur_b)?;
        let k = (prev << 2) | cur_idx;
        dh_nn += NN_DH[k];
        ds_nn += NN_DS[k];
        prev = cur_idx;

        let (h3, s3) = init_end(cur_b);
        let dh = dh_nn + h5 + h3;
        let ds = ds_nn + s5 + s3;
        let ds_salt = ds + 0.368 * (length as f64 - 1.0) * salt_ln;
        let denom = ds_salt + ct_term;
        let tm = (dh * 1000.0) / denom - 273.15;
        if tm >= tm_threshold {
            return Some(length);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Operating conditions per `program_instructions.md` / user spec.
    fn ref_params() -> TmParams {
        TmParams {
            oligo_conc_um: 0.2,
            na_conc_mm: 50.0,
            mg_conc_mm: 3.0,
            dntp_conc_mm: 0.8,
        }
    }

    /// All five reference values from example_sets/example_oligo_tm.csv
    /// must agree to within 2 °C.
    #[test]
    fn matches_reference_oligo_tm_csv() {
        let p = ref_params();
        let cases = [
            ("TCTGGACTGGAT", 44.4),
            ("GGGTCCCATATAAACCTG", 56.2),
            ("AAATGTGGCTAGGTCTTTGCCCACGCACCT", 73.9),
            ("AAYCGGTCCTCGTTTRCTTT", 62.1),
            ("TYAGGACMCGATGSACYA", 60.9),
        ];
        for (seq, expected) in cases {
            let got = calculate_tm(seq.as_bytes(), &p);
            let delta = (got - expected).abs();
            assert!(
                delta <= 2.0,
                "Tm for {seq}: predicted {got:.2}, expected {expected:.1}, off by {delta:.2}°C"
            );
        }
    }

    /// `determine_oligo_length` agrees with `calculate_tm` at the chosen
    /// length: the Tm at L meets the threshold and at L-1 does not.
    #[test]
    fn determine_oligo_length_consistent_with_tm() {
        let p = ref_params();
        let seq = b"GGGTCCCATATAAACCTGAAAGGCCTTACG"; // 30 nt
        let threshold = 60.0;
        let len = determine_oligo_length(seq, 0, threshold, &p)
            .expect("should reach 60°C within 30 nt");
        let tm_at = calculate_tm(&seq[..len], &p);
        assert!(tm_at >= threshold, "Tm at L={len} is {tm_at:.2}, below threshold");
        if len > 2 {
            let tm_before = calculate_tm(&seq[..len - 1], &p);
            assert!(
                tm_before < threshold,
                "Tm at L-1={} is {:.2}, already >= threshold",
                len - 1,
                tm_before
            );
        }
    }

    #[test]
    fn ambiguous_tm_is_median_of_variants() {
        let p = ref_params();
        // R = {A, G}. Sequence has exactly one R, so median = average of
        // the two variant Tms.
        let amb = b"ACGTACGTRCGTACGT";
        let v1 = b"ACGTACGTACGTACGT";
        let v2 = b"ACGTACGTGCGTACGT";
        let tm_a = calculate_tm(amb, &p);
        let tm_1 = calculate_tm(v1, &p);
        let tm_2 = calculate_tm(v2, &p);
        let mid = (tm_1 + tm_2) / 2.0;
        assert!(
            (tm_a - mid).abs() < 1e-9,
            "ambiguous Tm {tm_a:.4} != median {mid:.4}"
        );
    }
}
