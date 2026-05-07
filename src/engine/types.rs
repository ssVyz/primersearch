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
    pub na_concentration_mm: f64,
    pub mode: SearchMode,
    pub target_coverage_pct: f64,
    pub max_ambiguities: usize,
    pub exclude_n: bool,
    pub only_twofold: bool,
    pub orientation: Orientation,
    pub three_prime_match: usize,
}

impl SearchSettings {
    pub fn na_conc_molar(&self) -> f64 {
        self.na_concentration_mm / 1000.0
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
