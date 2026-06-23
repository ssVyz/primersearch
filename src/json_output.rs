//! JSON rendering of the search result for programmatic consumers — e.g. a
//! GUI front-end that drives primersearch as a subprocess and parses its
//! output. Opt-in via `--format json`; the default text output is unchanged.
//!
//! The engine types are intentionally serde-free (see `engine/mod.rs`), so
//! this module defines its own `Serialize` DTOs that copy from them rather
//! than deriving `Serialize` on the engine structs.

use std::io::{self, Write};

use serde::Serialize;

use crate::engine::{
    Orientation, PrimerSearchResult, QualityReport, SearchMode, SearchSettings,
};

/// Bump this when the JSON shape changes in a backward-incompatible way so
/// consumers can guard against schema drift.
const FORMAT_VERSION: u32 = 1;

#[derive(Serialize)]
struct Document<'a> {
    format_version: u32,
    preprocessing: Preprocessing,
    settings: Settings,
    result: ResultBlock<'a>,
}

#[derive(Serialize)]
struct Preprocessing {
    original_count: usize,
    valid_count: usize,
    majority_length: usize,
    removed: Removed,
}

#[derive(Serialize)]
struct Removed {
    gaps: usize,
    ambiguous: usize,
    invalid: usize,
    wrong_length: usize,
    total: usize,
}

#[derive(Serialize)]
struct Settings {
    /// `"no_ambiguities"` or `"incremental"`.
    mode: &'static str,
    fixed: bool,
    /// `"forward"` or `"reverse"`.
    orientation: &'static str,
    tm_threshold: f64,
    /// `false` in fixed mode, where the Tm threshold is reported but not used
    /// as a gate (injected oligos also bypass it regardless of this flag).
    tm_threshold_enforced: bool,
    oligo_concentration_um: f64,
    na_concentration_mm: f64,
    mg_concentration_mm: f64,
    dntp_concentration_mm: f64,
    target_coverage_pct: f64,
    max_ambiguities: usize,
    exclude_n: bool,
    only_twofold: bool,
    three_prime_match: usize,
    max_seeds: usize,
}

#[derive(Serialize)]
struct ResultBlock<'a> {
    total_sequences: usize,
    primer_count: usize,
    injected_count: usize,
    message: &'a str,
    primers: Vec<Primer>,
}

#[derive(Serialize)]
struct Primer {
    /// 1-based rank in the result list (matches the `#` column of the text
    /// output).
    index: usize,
    /// The primer as it should be ordered: already reverse-complemented for a
    /// reverse-orientation run, no display spacers.
    sequence: String,
    coverage_count: usize,
    coverage_pct: f64,
    /// Running total of `coverage_pct` through this primer (the `Total%`
    /// column of the text output).
    cumulative_pct: f64,
    tm: f64,
    ambiguity_count: usize,
    injected: bool,
    /// Zero-based inclusive start of the alignment range.
    align_start: usize,
    /// Zero-based exclusive end of the alignment range (half-open
    /// `[align_start, align_end)`). The primer spans `align_end - align_start`
    /// bases.
    align_end: usize,
    /// 1-based inclusive range exactly as shown in the text output, e.g.
    /// `"12-30"`.
    position_label: String,
}

/// Serialize the full run (preprocessing report + echoed settings + primers)
/// as a single pretty-printed JSON document followed by a trailing newline.
pub fn write_json<W: Write>(
    report: &QualityReport,
    result: &PrimerSearchResult,
    settings: &SearchSettings,
    w: &mut W,
) -> io::Result<()> {
    let removed = Removed {
        gaps: report.removed_gaps,
        ambiguous: report.removed_ambiguous,
        invalid: report.removed_invalid,
        wrong_length: report.removed_wrong_length,
        total: report.removed_gaps
            + report.removed_ambiguous
            + report.removed_invalid
            + report.removed_wrong_length,
    };

    let mut cumulative = 0.0_f64;
    let primers: Vec<Primer> = result
        .primers
        .iter()
        .enumerate()
        .map(|(i, p)| {
            cumulative += p.coverage_pct;
            Primer {
                index: i + 1,
                sequence: String::from_utf8_lossy(&p.sequence).into_owned(),
                coverage_count: p.coverage_count,
                coverage_pct: p.coverage_pct,
                cumulative_pct: cumulative,
                tm: p.tm,
                ambiguity_count: p.ambiguity_count,
                injected: p.injected,
                align_start: p.align_start,
                align_end: p.align_end,
                position_label: format!("{}-{}", p.align_start + 1, p.align_end),
            }
        })
        .collect();

    let injected_count = result.primers.iter().filter(|p| p.injected).count();

    let doc = Document {
        format_version: FORMAT_VERSION,
        preprocessing: Preprocessing {
            original_count: report.original_count,
            valid_count: report.valid_count,
            majority_length: report.majority_length,
            removed,
        },
        settings: Settings {
            mode: match settings.mode {
                SearchMode::NoAmbiguities => "no_ambiguities",
                SearchMode::Incremental => "incremental",
            },
            fixed: settings.fixed,
            orientation: match settings.orientation {
                Orientation::Forward => "forward",
                Orientation::Reverse => "reverse",
            },
            tm_threshold: settings.tm_threshold,
            tm_threshold_enforced: !settings.fixed,
            oligo_concentration_um: settings.oligo_concentration_um,
            na_concentration_mm: settings.na_concentration_mm,
            mg_concentration_mm: settings.mg_concentration_mm,
            dntp_concentration_mm: settings.dntp_concentration_mm,
            target_coverage_pct: settings.target_coverage_pct,
            max_ambiguities: settings.max_ambiguities,
            exclude_n: settings.exclude_n,
            only_twofold: settings.only_twofold,
            three_prime_match: settings.three_prime_match,
            max_seeds: settings.max_seeds,
        },
        result: ResultBlock {
            total_sequences: result.total_sequences,
            primer_count: result.primers.len(),
            injected_count,
            message: &result.message,
            primers,
        },
    };

    let s = serde_json::to_string_pretty(&doc)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
    writeln!(w, "{s}")
}
