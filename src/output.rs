//! Text rendering of the quality report and the primer search results. The
//! layout intentionally matches the reference Python output (semantically,
//! not byte-for-byte): same headers, same columns, same numeric formatting.

use std::io::Write;

use crate::engine::{Orientation, PrimerSearchResult, QualityReport, SearchMode, SearchSettings};

pub fn format_quality_report<W: Write>(r: &QualityReport, w: &mut W) -> std::io::Result<()> {
    writeln!(w, "PRIMER SEARCH PREPROCESSING")?;
    writeln!(w, "{}", "=".repeat(28))?;
    writeln!(
        w,
        "Original sequences:         {}",
        thousands(r.original_count)
    )?;
    writeln!(
        w,
        "Valid sequences:            {}",
        thousands(r.valid_count)
    )?;
    writeln!(w, "Majority sequence length:   {} bp", r.majority_length)?;
    writeln!(w)?;
    writeln!(w, "Removed due to quality issues:")?;
    writeln!(
        w,
        "  - Gaps (- or .):          {}",
        thousands(r.removed_gaps)
    )?;
    writeln!(
        w,
        "  - Ambiguous bases:        {}",
        thousands(r.removed_ambiguous)
    )?;
    writeln!(
        w,
        "  - Invalid characters:     {}",
        thousands(r.removed_invalid)
    )?;
    writeln!(
        w,
        "  - Wrong length:           {}",
        thousands(r.removed_wrong_length)
    )?;
    let total =
        r.removed_gaps + r.removed_ambiguous + r.removed_invalid + r.removed_wrong_length;
    writeln!(w, "  - Total removed:          {}", thousands(total))?;
    writeln!(w)?;
    Ok(())
}

pub fn format_results<W: Write>(
    result: &PrimerSearchResult,
    settings: &SearchSettings,
    add_spacer: bool,
    w: &mut W,
) -> std::io::Result<()> {
    let bar = "=".repeat(90);
    let dash = "-".repeat(90);
    writeln!(w, "{bar}")?;
    writeln!(w, "PRIMER SEARCH RESULTS")?;
    writeln!(w, "{bar}")?;
    writeln!(w)?;

    let mode_label = match settings.mode {
        SearchMode::NoAmbiguities => "No Ambiguities (exact match)".to_string(),
        SearchMode::Incremental => format!(
            "Incremental (target {:.0}%, max {} ambiguities)",
            settings.target_coverage_pct, settings.max_ambiguities,
        ),
    };
    let orient_label = match settings.orientation {
        Orientation::Forward => "Forward (sense)",
        Orientation::Reverse => "Reverse (anti-sense)",
    };
    writeln!(w, "Search Mode:      {mode_label}")?;
    writeln!(w, "Orientation:      {orient_label}")?;
    writeln!(
        w,
        "Total Sequences:  {}",
        thousands(result.total_sequences)
    )?;
    writeln!(w, "Primers Found:    {}", result.primers.len())?;
    writeln!(w, "Tm Threshold:     {:.1}°C", settings.tm_threshold)?;
    writeln!(w, "Oligo Conc:       {:.3} µM", settings.oligo_concentration_um)?;
    writeln!(w, "Na+ Conc:         {:.1} mM", settings.na_concentration_mm)?;
    writeln!(w, "Mg2+ Conc:        {:.2} mM", settings.mg_concentration_mm)?;
    writeln!(w, "dNTP Conc:        {:.2} mM", settings.dntp_concentration_mm)?;
    if settings.three_prime_match > 0 {
        writeln!(w, "3' Perfect Match: {} bases", settings.three_prime_match)?;
    }
    if matches!(settings.mode, SearchMode::Incremental) {
        if settings.exclude_n {
            writeln!(w, "Exclude N:        Yes")?;
        }
        if settings.only_twofold {
            writeln!(w, "Only 2-fold:      Yes")?;
        }
    }
    if !result.message.is_empty() {
        writeln!(w, "Note:             {}", result.message)?;
    }
    writeln!(w)?;
    writeln!(w, "{dash}")?;
    writeln!(w)?;

    if result.primers.is_empty() {
        writeln!(w, "No primers found.")?;
        writeln!(w)?;
        writeln!(w, "{bar}")?;
        return Ok(());
    }

    let formatted_seqs: Vec<String> = result
        .primers
        .iter()
        .map(|p| {
            let s = std::str::from_utf8(&p.sequence)
                .unwrap_or("?")
                .to_string();
            if add_spacer { add_spaces(&s, 3) } else { s }
        })
        .collect();
    let num_w = result.primers.len().to_string().len().max(1);
    let max_seq_len = formatted_seqs
        .iter()
        .map(|s| s.len())
        .max()
        .unwrap_or(0)
        .max("Sequence".len());
    let max_count_len = result
        .primers
        .iter()
        .map(|p| thousands(p.coverage_count).len())
        .max()
        .unwrap_or(0)
        .max("Count".len());

    writeln!(
        w,
        "{:>w1$}   {:<w2$}   {:>w3$}      %    Total%   Tm(°C)   Position",
        "#",
        "Sequence",
        "Count",
        w1 = num_w,
        w2 = max_seq_len,
        w3 = max_count_len,
    )?;
    writeln!(w, "{dash}")?;

    let mut cumulative = 0.0_f64;
    for (i, p) in result.primers.iter().enumerate() {
        cumulative += p.coverage_pct;
        let pct_s = format!("{:.1}%", p.coverage_pct);
        let cum_s = format!("{:.1}%", cumulative);
        let tm_s = format!("{:.1}", p.tm);
        writeln!(
            w,
            "{:>w1$}   {:<w2$}   {:>w3$}   {:>6}   {:>7}   {:>6}   {}-{}",
            i + 1,
            formatted_seqs[i],
            thousands(p.coverage_count),
            pct_s,
            cum_s,
            tm_s,
            p.align_start + 1,
            p.align_end,
            w1 = num_w,
            w2 = max_seq_len,
            w3 = max_count_len,
        )?;
    }

    writeln!(w)?;
    writeln!(w, "{bar}")?;
    Ok(())
}

fn thousands(n: usize) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

fn add_spaces(s: &str, group: usize) -> String {
    let mut out = String::with_capacity(s.len() + s.len() / group.max(1));
    for (i, c) in s.chars().enumerate() {
        if i > 0 && i % group == 0 {
            out.push(' ');
        }
        out.push(c);
    }
    out
}
