//! Public types for the engine — kept free of CLI / serde concerns so the
//! engine remains self-contained.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchMode {
    NoAmbiguities,
    Incremental,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    Reverse,
}

#[derive(Debug, Clone)]
pub struct SearchSettings {
    pub tm_threshold: f64,
    pub oligo_concentration_um: f64,
    pub na_concentration_mm: f64,
    pub mg_concentration_mm: f64,
    pub dntp_concentration_mm: f64,
    pub mode: SearchMode,
    pub target_coverage_pct: f64,
    pub max_ambiguities: usize,
    pub exclude_n: bool,
    pub only_twofold: bool,
    pub orientation: Orientation,
    pub three_prime_match: usize,
    /// Per-range cap on the number of unique subsequences tried as
    /// seeds in incremental mode. `0` disables the cap (try every unique
    /// subsequence).
    pub max_seeds: usize,
    /// Fixed-slice mode. When `true`, the search is skipped: the entire
    /// input alignment is treated as a single user-chosen slice and the
    /// engine only generates the variant(s) needed to cover every input
    /// sequence at that slice. The configured `mode` still selects how
    /// variants are formed (exact vs. IUPAC consensus), but `tm_threshold`
    /// is not used as a gate — see [`crate::engine::find_primers_fixed`].
    pub fixed: bool,
}

impl SearchSettings {
    pub fn tm_params(&self) -> crate::engine::tm::TmParams {
        crate::engine::tm::TmParams {
            oligo_conc_um: self.oligo_concentration_um,
            na_conc_mm: self.na_concentration_mm,
            mg_conc_mm: self.mg_concentration_mm,
            dntp_conc_mm: self.dntp_concentration_mm,
        }
    }
    pub fn is_reverse(&self) -> bool {
        matches!(self.orientation, Orientation::Reverse)
    }
}

#[derive(Debug, Clone)]
pub struct PrimerCandidate {
    /// The primer sequence as it should be displayed. Already reverse-
    /// complemented when `orientation == Reverse`.
    pub sequence: Vec<u8>,
    pub coverage_count: usize,
    pub coverage_pct: f64,
    pub tm: f64,
    /// Alignment range [start, end), zero-based, in the *forward* coordinates
    /// of the input alignment.
    pub align_start: usize,
    pub align_end: usize,
    pub ambiguity_count: usize,
    /// `true` if this primer was supplied by the user via `--inject` rather
    /// than discovered by the search. Injected oligos are obligatory: they
    /// are placed before the search runs and emitted regardless of Tm or the
    /// variant-generation constraints.
    pub injected: bool,
}

#[derive(Debug, Clone, Default)]
pub struct PrimerSearchResult {
    pub primers: Vec<PrimerCandidate>,
    pub total_sequences: usize,
    pub message: String,
}

#[derive(Debug, Clone, Default)]
pub struct QualityReport {
    pub original_count: usize,
    pub valid_count: usize,
    pub removed_ambiguous: usize,
    pub removed_gaps: usize,
    pub removed_wrong_length: usize,
    pub removed_invalid: usize,
    pub majority_length: usize,
    pub valid_sequences: Vec<Vec<u8>>,
}

/// Progress sink. Implementors must be `Sync + Send` so the engine can call
/// `report` from any worker thread.
pub trait Progress: Sync + Send {
    fn report(&self, message: &str, pct: f64);
    fn cancelled(&self) -> bool {
        false
    }
}

pub struct NoProgress;
impl Progress for NoProgress {
    fn report(&self, _: &str, _: f64) {}
}
