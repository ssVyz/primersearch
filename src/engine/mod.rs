//! Primer-search analysis engine.
//!
//! This module is intentionally self-contained: it depends only on `std` and
//! `rayon` and does not touch the filesystem, CLI, or any serialization
//! framework. To lift it into another project, copy the `engine/` directory.
//!
//! Some re-exports below are unused by the bundled CLI but are part of the
//! engine's public surface for external embedders.

#![allow(dead_code, unused_imports)]

mod fasta;
mod iupac;
mod search;
pub mod tm;
mod types;

pub use fasta::{parse_fasta, quality_filter};
pub use iupac::{is_ambiguous, reverse_complement};
pub use search::{find_primers, find_primers_fixed};
pub use tm::{calculate_tm, determine_oligo_length, TmParams};
pub use types::{
    NoProgress, Orientation, PrimerCandidate, PrimerSearchResult, Progress, QualityReport,
    SearchMode, SearchSettings,
};
