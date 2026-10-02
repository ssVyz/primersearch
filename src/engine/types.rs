//! Public types for the engine — kept free of CLI / serde concerns so the
//! engine remains self-contained.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchMode {
    NoAmbiguities,
    Incremental,
    /// Exhaustive search for the best set of `n` oligos tolerating mismatches.
    /// Handled by [`crate::engine::find_primers_by_mismatch`], not by the
    /// greedy round loops; see [`MismatchSettings`].
    OptimizeByMismatch,
}

/// Coverage criterion of [`SearchMode::OptimizeByMismatch`]. In both cases a
/// sequence is scored by its *best-matching* oligo in the set (fewest
/// mismatches).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MismatchOp {
    /// A sequence counts when its best match has at most `mismatches`
    /// mismatches.
    LowerOrEqual,
    /// A sequence counts only when its best match has exactly `mismatches`
    /// mismatches; sequences matched better than that do not count.
    Exact,
}

/// Parameters of [`SearchMode::OptimizeByMismatch`]. Ignored by the other
/// modes.
#[derive(Debug, Clone)]
pub struct MismatchSettings {
    pub op: MismatchOp,
    /// Number of oligos in the optimised set (`n`). Injected oligos count
    /// toward it.
    pub oligo_count: usize,
    /// Mismatch count `x` of the coverage criterion.
    pub mismatches: usize,
    /// Ambiguity codes per candidate oligo (`y`); fewer only where a window
    /// has fewer eligible variable positions.
    pub ambiguities: usize,
    /// Upper bound on the candidate oligos enumerated, checked before any
    /// work starts. `0` = no limit.
    pub max_candidates: u64,
    /// Upper bound on candidate evaluations during the set search. `0` = no
    /// limit.
    pub max_work: u64,
}

impl Default for MismatchSettings {
    fn default() -> Self {
        Self {
            op: MismatchOp::LowerOrEqual,
            oligo_count: 1,
            mismatches: 0,
            ambiguities: 0,
            max_candidates: DEFAULT_MAX_CANDIDATES,
            max_work: DEFAULT_MAX_WORK,
        }
    }
}

pub const DEFAULT_MAX_CANDIDATES: u64 = 500_000_000;
pub const DEFAULT_MAX_WORK: u64 = 2_000_000_000;

/// Set-level summary produced by [`SearchMode::OptimizeByMismatch`]. All
/// counts are numbers of input sequences.
#[derive(Debug, Clone, Default)]
pub struct MismatchReport {
    /// `level_counts[j]` = sequences whose best-matching oligo in the set has
    /// exactly `j` mismatches, for `j = 0..=mismatches`.
    pub level_counts: Vec<usize>,
    /// Sequences with more than `mismatches` mismatches to every oligo, or a
    /// mismatch in the protected 3' region.
    pub not_covered: usize,
    /// Sequences meeting the coverage criterion (the optimised objective).
    pub counted: usize,
    /// Alignment windows searched (1 in fixed-slice mode).
    pub windows: usize,
    /// Distinct candidate oligos generated across all windows.
    pub candidates_generated: u64,
    /// Candidates left after dropping duplicates and dominated ones.
    pub candidates_reduced: usize,
    /// Candidate evaluations spent in the set search. With several threads
    /// it can vary slightly from run to run; the chosen set does not.
    pub evaluations: u64,
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
    /// Parameters of [`SearchMode::OptimizeByMismatch`].
    pub mismatch: MismatchSettings,
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
    /// Only set by [`SearchMode::OptimizeByMismatch`].
    pub mismatch: Option<MismatchReport>,
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

/// Stage of a run, carried by [`ProgressEvent`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProgressPhase {
    /// Placing injected oligos (any mode).
    Inject,
    /// A greedy search round: collecting and evaluating ranges, round done.
    Round,
    /// Fixed-slice variant generation (greedy modes).
    Fixed,
    /// Optimize-by-mismatch: collecting windows.
    Windows,
    /// Optimize-by-mismatch: per-window candidate generation.
    Candidates,
    /// Optimize-by-mismatch: reducing the candidate pool.
    Reduce,
    /// Optimize-by-mismatch: branch-and-bound set search.
    SetSearch,
}

impl ProgressPhase {
    /// Short machine-readable name, e.g. `"set_search"`.
    pub fn as_str(self) -> &'static str {
        match self {
            ProgressPhase::Inject => "inject",
            ProgressPhase::Round => "round",
            ProgressPhase::Fixed => "fixed",
            ProgressPhase::Windows => "windows",
            ProgressPhase::Candidates => "candidates",
            ProgressPhase::Reduce => "reduce",
            ProgressPhase::SetSearch => "set_search",
        }
    }
}

/// A structured progress report: the message and percentage of
/// [`Progress::report`] plus the phase and the counters known at that point,
/// as `(name, value)` pairs (e.g. `[("done", 12), ("total", 40)]`).
#[derive(Debug, Clone, Copy)]
pub struct ProgressEvent<'a> {
    pub message: &'a str,
    pub pct: f64,
    pub phase: Option<ProgressPhase>,
    pub counters: &'a [(&'static str, u64)],
}

impl<'a> ProgressEvent<'a> {
    pub fn new(
        phase: ProgressPhase,
        message: &'a str,
        pct: f64,
        counters: &'a [(&'static str, u64)],
    ) -> Self {
        Self {
            message,
            pct,
            phase: Some(phase),
            counters,
        }
    }
}

/// Progress sink. Implementors must be `Sync + Send` so the engine can call
/// `report` from any worker thread.
pub trait Progress: Sync + Send {
    fn report(&self, message: &str, pct: f64);
    /// Structured form of [`report`](Self::report), used by the engine for
    /// every report. The default forwards the message and percentage, so a
    /// sink that only shows text implements `report` alone.
    fn report_event(&self, event: &ProgressEvent<'_>) {
        self.report(event.message, event.pct);
    }
    fn cancelled(&self) -> bool {
        false
    }
}

pub struct NoProgress;
impl Progress for NoProgress {
    fn report(&self, _: &str, _: f64) {}
}
