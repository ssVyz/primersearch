//! FASTA parsing and quality filtering.

use std::collections::HashMap;

use crate::engine::iupac::is_ambiguous;
use crate::engine::types::QualityReport;

/// Parse FASTA-format text into raw uppercase byte sequences. Whitespace
/// inside sequence lines is stripped. Records are separated by `>` headers.
pub fn parse_fasta(text: &str) -> Vec<Vec<u8>> {
    let mut sequences: Vec<Vec<u8>> = Vec::new();
    let mut current: Vec<u8> = Vec::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if !current.is_empty() {
                sequences.push(std::mem::take(&mut current));
            }
        } else {
            for &b in line.as_bytes() {
                if !b.is_ascii_whitespace() {
                    current.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    if !current.is_empty() {
        sequences.push(current);
    }
    sequences
}

/// Filter sequences, keeping only those matching the majority alignment
/// length and containing only A / C / G / T. Tracks why each rejection
/// happened.
pub fn quality_filter(sequences: Vec<Vec<u8>>) -> QualityReport {
    let mut report = QualityReport {
        original_count: sequences.len(),
        ..Default::default()
    };
    if sequences.is_empty() {
        return report;
    }

    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    for s in &sequences {
        *length_counts.entry(s.len()).or_insert(0) += 1;
    }
    // Pick the most common length deterministically: highest count, then
    // smallest length on tie.
    let mut entries: Vec<(usize, usize)> = length_counts.into_iter().collect();
    entries.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    report.majority_length = entries.first().map(|(l, _)| *l).unwrap_or(0);

    let mut valid: Vec<Vec<u8>> = Vec::new();
    for seq in sequences {
        if seq.iter().any(|&b| b == b'-' || b == b'.') {
            report.removed_gaps += 1;
            continue;
        }
        if seq.iter().any(|&b| is_ambiguous(b)) {
            report.removed_ambiguous += 1;
            continue;
        }
        if !seq
            .iter()
            .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        {
            report.removed_invalid += 1;
            continue;
        }
        if seq.len() != report.majority_length {
            report.removed_wrong_length += 1;
            continue;
        }
        valid.push(seq);
    }
    report.valid_count = valid.len();
    report.valid_sequences = valid;
    report
}
