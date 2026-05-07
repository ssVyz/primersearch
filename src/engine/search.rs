//! Greedy round-based primer search.
//!
//! Each round:
//!  1. (single-threaded) For every (sequence, start_pos) pair among the
//!     remaining sequences, walk forward to find the smallest length whose
//!     subsequence reaches the Tm threshold. Collect the unique
//!     `(start, end)` alignment ranges seen.
//!  2. (parallel) Evaluate every unique range against every remaining
//!     sequence.
//!  3. (single-threaded) Pick the best candidate with a stable tie-break,
//!     remove its covered sequences from the pool, repeat.
//!
//! The tie-break (max coverage_count, then lowest start, lowest end, fewest
//! ambiguities, lex-smallest sequence) is deterministic and identical in
//! single-threaded and parallel execution — running with `--threads 1` and
//! `--threads 16` produces the same primer set.

use std::collections::{BTreeSet, HashMap, HashSet};

use rayon::prelude::*;

use crate::engine::iupac::{base_mask, is_ambiguous, mask_to_iupac, reverse_complement, MASK_C, MASK_G};
use crate::engine::types::{
    Orientation, PrimerCandidate, PrimerSearchResult, Progress, SearchMode, SearchSettings,
};

// ---------------------------------------------------------------------------
// Tm
// ---------------------------------------------------------------------------

/// Salt-adjusted Tm:  81.5 + 16.6·log10([Na+]) + 0.41·%GC − 675/n.
///
/// IUPAC ambiguity codes contribute their *fractional* GC content (e.g. an
/// `R = {A, G}` position contributes 0.5 GC). For pure A/C/G/T inputs this
/// is identical to the integer count.
pub fn calculate_tm(sequence: &[u8], na_conc_molar: f64) -> f64 {
    let n = sequence.len();
    if n == 0 {
        return 0.0;
    }
    let mut gc_score = 0.0_f64;
    for &b in sequence {
        let mask = base_mask(b);
        let total = mask.count_ones();
        if total == 0 {
            continue;
        }
        let gc = (mask & (MASK_G | MASK_C)).count_ones();
        gc_score += gc as f64 / total as f64;
    }
    let gc_pct = (gc_score / n as f64) * 100.0;
    let na = na_conc_molar.max(1e-10);
    81.5 + 16.6 * na.log10() + 0.41 * gc_pct - 675.0 / (n as f64)
}

/// Smallest oligo length starting at `start` whose Tm reaches the threshold,
/// or `None` if extending to the end of `seq` does not. Computed
/// incrementally — O(L) per call, not O(L²) like the naïve form.
pub fn determine_oligo_length(
    seq: &[u8],
    start: usize,
    tm_threshold: f64,
    na_conc_molar: f64,
) -> Option<usize> {
    if start >= seq.len() {
        return None;
    }
    let max_len = seq.len() - start;
    let na = na_conc_molar.max(1e-10);
    let na_log = na.log10();
    let mut gc_count: usize = 0;
    for length in 1..=max_len {
        let added = seq[start + length - 1];
        if added == b'G' || added == b'C' {
            gc_count += 1;
        }
        let n = length as f64;
        let gc_pct = (gc_count as f64 / n) * 100.0;
        let tm = 81.5 + 16.6 * na_log + 0.41 * gc_pct - 675.0 / n;
        if tm >= tm_threshold {
            return Some(length);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn find_primers(
    sequences: &[Vec<u8>],
    settings: &SearchSettings,
    progress: &dyn Progress,
) -> PrimerSearchResult {
    let total = sequences.len();
    let mut result = PrimerSearchResult {
        total_sequences: total,
        ..Default::default()
    };
    if total == 0 {
        result.message = "No sequences to analyze".to_string();
        return result;
    }
    let na_molar = settings.na_conc_molar();
    let is_reverse = matches!(settings.orientation, Orientation::Reverse);
    let seq_len = sequences[0].len();

    let mut remaining: BTreeSet<usize> = (0..total).collect();
    let mut round = 0;

    while !remaining.is_empty() {
        if progress.cancelled() {
            result.message = "Search cancelled".to_string();
            break;
        }
        round += 1;
        let remaining_list: Vec<usize> = remaining.iter().copied().collect();
        let remaining_seqs: Vec<&[u8]> = remaining_list
            .iter()
            .map(|&i| sequences[i].as_slice())
            .collect();

        progress.report(
            &format!(
                "Round {}: collecting candidate ranges ({} sequences remaining)",
                round,
                remaining_list.len()
            ),
            0.0,
        );

        // Phase 1: collect unique alignment ranges (start, end). Each range
        // came from at least one remaining sequence whose slice reaches the
        // Tm threshold; we don't need to remember which one.
        let mut seen: HashSet<(usize, usize)> = HashSet::new();
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        for seq in remaining_seqs.iter() {
            for start in 0..seq_len {
                if let Some(len) =
                    determine_oligo_length(seq, start, settings.tm_threshold, na_molar)
                {
                    let end = start + len;
                    if end > seq_len {
                        continue;
                    }
                    if seen.insert((start, end)) {
                        ranges.push((start, end));
                    }
                }
            }
        }

        if ranges.is_empty() {
            result.message = "Could not find any more primers".to_string();
            break;
        }

        progress.report(
            &format!("Round {}: evaluating {} unique ranges", round, ranges.len()),
            10.0,
        );

        // Phase 2: evaluate ranges in parallel, collect candidates.
        let candidates: Vec<EvalResult> = ranges
            .par_iter()
            .filter_map(|&(start, end)| {
                if progress.cancelled() {
                    return None;
                }
                evaluate_range(start, end, &remaining_seqs, settings, na_molar, is_reverse)
            })
            .collect();

        if candidates.is_empty() {
            result.message = "Could not find any more primers".to_string();
            break;
        }

        // Phase 3: reduce. Match the reference Python (`if count >
        // best_coverage:`) — strict greater wins, ties keep the first
        // candidate in iteration order. `rayon::par_iter().collect()`
        // preserves order, so `candidates` is in the same order as
        // `ranges`, which is in (seq_i, start_pos) iteration order.
        let mut best_idx = 0usize;
        let mut best_count = candidates[0].count;
        for (i, c) in candidates.iter().enumerate().skip(1) {
            if c.count > best_count {
                best_count = c.count;
                best_idx = i;
            }
        }
        let best = candidates.into_iter().nth(best_idx).expect("non-empty");

        // Map local-to-this-round indices back to global indices.
        let global_covered: BTreeSet<usize> = best
            .local_indices
            .iter()
            .map(|&j| remaining_list[j])
            .collect();

        let display_seq = if is_reverse {
            reverse_complement(&best.sequence)
        } else {
            best.sequence.clone()
        };
        let coverage_pct = (best.count as f64 / total as f64) * 100.0;
        let primer = PrimerCandidate {
            sequence: display_seq,
            coverage_count: best.count,
            coverage_pct,
            tm: best.tm,
            align_start: best.start,
            align_end: best.end,
            ambiguity_count: best.amb_count,
        };

        for idx in &global_covered {
            remaining.remove(idx);
        }
        result.primers.push(primer);
        let covered_so_far = total - remaining.len();
        progress.report(
            &format!(
                "Round {} done: primer covers {} sequences ({:.1}%). Total: {}/{} ({:.1}%)",
                round,
                best.count,
                coverage_pct,
                covered_so_far,
                total,
                covered_so_far as f64 / total as f64 * 100.0
            ),
            100.0,
        );
    }

    result
}

// ---------------------------------------------------------------------------
// Per-range evaluation
// ---------------------------------------------------------------------------

struct EvalResult {
    start: usize,
    end: usize,
    sequence: Vec<u8>,
    count: usize,
    local_indices: Vec<usize>,
    amb_count: usize,
    tm: f64,
}

fn evaluate_range(
    start: usize,
    end: usize,
    remaining_seqs: &[&[u8]],
    settings: &SearchSettings,
    na_molar: f64,
    is_reverse: bool,
) -> Option<EvalResult> {
    let subsequences: Vec<&[u8]> = remaining_seqs.iter().map(|s| &s[start..end]).collect();

    let (variant, count, local_indices, amb_count) = match settings.mode {
        SearchMode::NoAmbiguities => {
            let (v, c, l) = eval_no_ambiguities(&subsequences);
            if c == 0 {
                return None;
            }
            (v, c, l, 0)
        }
        SearchMode::Incremental => eval_incremental(
            &subsequences,
            settings.target_coverage_pct,
            settings.max_ambiguities,
            settings.exclude_n,
            settings.only_twofold,
            settings.three_prime_match,
            is_reverse,
        )?,
    };

    // Tm of the actual displayed primer (the consensus / chosen variant).
    let tm = calculate_tm(&variant, na_molar);

    Some(EvalResult {
        start,
        end,
        sequence: variant,
        count,
        local_indices,
        amb_count,
        tm,
    })
}

/// No-ambiguity mode: the most frequent exact subsequence wins. Ties broken
/// by *first-appearance order* in `subsequences` — matches the reference
/// Python implementation, which iterates a dict in insertion order with a
/// `>` (strict-greater) update.
fn eval_no_ambiguities(subsequences: &[&[u8]]) -> (Vec<u8>, usize, Vec<usize>) {
    // First-appearance order: build a Vec in insertion order with a HashMap
    // for fast lookup of which entry to extend.
    let mut index_of: HashMap<&[u8], usize> = HashMap::new();
    let mut entries: Vec<(&[u8], Vec<usize>)> = Vec::new();
    for (i, s) in subsequences.iter().enumerate() {
        if let Some(&idx) = index_of.get(s) {
            entries[idx].1.push(i);
        } else {
            index_of.insert(*s, entries.len());
            entries.push((*s, vec![i]));
        }
    }
    // Strict `>` update — first to reach a given count wins.
    let mut best_seq: &[u8] = b"";
    let mut best: Vec<usize> = Vec::new();
    let mut best_count: usize = 0;
    for (seq, idx) in entries.into_iter() {
        if idx.len() > best_count {
            best_count = idx.len();
            best_seq = seq;
            best = idx;
        }
    }
    (best_seq.to_vec(), best_count, best)
}

/// Incremental mode: try increasing ambiguity budgets; for each budget try
/// up to 50 distinct subsequences as seeds. For each seed, greedily absorb
/// other unique subsequences as long as the resulting consensus stays within
/// the budget. Pick the best consensus across all seeds.
fn eval_incremental(
    subsequences: &[&[u8]],
    target_pct: f64,
    max_ambiguities: usize,
    exclude_n: bool,
    only_twofold: bool,
    three_prime_len: usize,
    is_reverse: bool,
) -> Option<(Vec<u8>, usize, Vec<usize>, usize)> {
    if subsequences.is_empty() {
        return None;
    }

    // Build `unique` in first-appearance order — required to match the
    // Python iteration order (dict insertion order). The greedy expansion
    // loop below depends on this order.
    let mut index_of: HashMap<&[u8], usize> = HashMap::new();
    let mut unique: Vec<(&[u8], Vec<usize>)> = Vec::new();
    for (i, s) in subsequences.iter().enumerate() {
        if let Some(&idx) = index_of.get(s) {
            unique[idx].1.push(i);
        } else {
            index_of.insert(*s, unique.len());
            unique.push((*s, vec![i]));
        }
    }

    let total = subsequences.len();
    let target_count = (target_pct / 100.0 * total as f64).ceil() as usize;

    // Seed iteration order: by count descending, ties broken by first-
    // appearance order (Python's `sorted()` is stable, so `sorted(seqs,
    // key=-count)` preserves insertion order on ties).
    let mut seed_order: Vec<usize> = (0..unique.len()).collect();
    seed_order.sort_by(|&a, &b| unique[b].1.len().cmp(&unique[a].1.len()));
    let max_seeds = seed_order.len().min(50);

    let mut best_consensus: Option<Vec<u8>> = None;
    let mut best_count: usize = 0;
    let mut best_amb: usize = usize::MAX;
    let mut best_local: Vec<usize> = Vec::new();
    let mut found_target = false;

    for amb_level in 0..=max_ambiguities {
        if found_target {
            break;
        }
        for &seed_idx in seed_order.iter().take(max_seeds) {
            let seed_seq = unique[seed_idx].0;

            let mut group: Vec<&[u8]> = Vec::with_capacity(unique.len());
            group.push(seed_seq);
            // Inner expansion: iterate `unique` in first-appearance order
            // (NOT seed_order). Matches Python's
            // `for other_seq in unique_seqs:` over insertion-ordered keys.
            for (oi, (other_seq, _)) in unique.iter().enumerate() {
                if oi == seed_idx {
                    continue;
                }
                group.push(*other_seq);
                let ok = match create_consensus(&group, exclude_n, only_twofold) {
                    Some((_, amb_count)) => amb_count <= amb_level,
                    None => false,
                };
                if !ok {
                    group.pop();
                }
            }

            let (consensus, amb_count) = match create_consensus(&group, exclude_n, only_twofold) {
                Some((c, a)) if a <= amb_level => (c, a),
                _ => continue,
            };

            if !check_three_prime_ok(&consensus, three_prime_len, is_reverse) {
                continue;
            }

            // Total coverage = sum of bucket counts for every unique
            // subseq that this consensus matches.
            let mut coverage_count: usize = 0;
            let mut local_indices: Vec<usize> = Vec::new();
            for (uniq_seq, indices) in unique.iter() {
                if sequence_matches_consensus(uniq_seq, &consensus) {
                    coverage_count += indices.len();
                    local_indices.extend_from_slice(indices);
                }
            }

            let is_better = if found_target {
                coverage_count > best_count
                    || (coverage_count == best_count && amb_count < best_amb)
            } else if coverage_count >= target_count {
                true
            } else {
                coverage_count > best_count
                    || (coverage_count == best_count && amb_count < best_amb)
            };

            if is_better {
                best_consensus = Some(consensus);
                best_count = coverage_count;
                best_amb = amb_count;
                best_local = local_indices;
                if best_count >= target_count {
                    found_target = true;
                }
            }
        }
    }

    if best_consensus.is_none() || best_count == 0 {
        // Fallback: most frequent exact subsequence.
        let (seq, indices) = unique.into_iter().next()?;
        let count = indices.len();
        return Some((seq.to_vec(), count, indices, 0));
    }

    best_local.sort();
    Some((best_consensus.expect("checked above"), best_count, best_local, best_amb))
}

/// Build an IUPAC consensus over a group of (filtered, A/C/G/T-only)
/// subsequences. Returns `(consensus, ambiguity_count)`, or `None` if a
/// forbidden code (N when `exclude_n`, or a 3-/4-fold code when
/// `only_twofold`) would be required.
fn create_consensus(seqs: &[&[u8]], exclude_n: bool, only_twofold: bool) -> Option<(Vec<u8>, usize)> {
    if seqs.is_empty() {
        return Some((Vec::new(), 0));
    }
    let len = seqs[0].len();
    let mut consensus = Vec::with_capacity(len);
    let mut amb_count = 0;
    for pos in 0..len {
        let mut mask = 0u8;
        for s in seqs {
            mask |= base_mask(s[pos]);
        }
        let popcount = mask.count_ones();
        if popcount == 1 {
            consensus.push(mask_to_iupac(mask));
        } else {
            let code = mask_to_iupac(mask);
            if exclude_n && code == b'N' {
                return None;
            }
            if only_twofold && popcount > 2 {
                return None;
            }
            consensus.push(code);
            amb_count += 1;
        }
    }
    Some((consensus, amb_count))
}

/// True iff every base of `seq` (assumed A/C/G/T) is in the IUPAC code set
/// of the corresponding consensus position.
fn sequence_matches_consensus(seq: &[u8], consensus: &[u8]) -> bool {
    if seq.len() != consensus.len() {
        return false;
    }
    seq.iter()
        .zip(consensus.iter())
        .all(|(s, c)| base_mask(*c) & base_mask(*s) != 0)
}

/// Verify that the 3'-end region is free of ambiguity codes. For forward
/// orientation that's the last `n` bases; for reverse orientation it's the
/// first `n` bases (since the displayed sequence will be reverse-
/// complemented before output).
fn check_three_prime_ok(consensus: &[u8], three_prime_len: usize, is_reverse: bool) -> bool {
    if three_prime_len == 0 {
        return true;
    }
    let n = three_prime_len.min(consensus.len());
    let region = if is_reverse {
        &consensus[..n]
    } else {
        &consensus[consensus.len() - n..]
    };
    !region.iter().any(|&c| is_ambiguous(c))
}
