//! Greedy round-based primer search.
//!
//! Each round:
//!  1. (single-threaded) For every (sequence, start_pos) pair among the
//!     remaining sequences, walk forward to find the smallest length whose
//!     subsequence reaches the Tm threshold. Collect the unique
//!     `(start, end)` alignment ranges seen.
//!  2. (parallel) Evaluate every unique range against every remaining
//!     sequence. Only variants whose own Tm meets the threshold may be
//!     chosen; in incremental mode this constraint applies to the *seed*
//!     (an exact A/C/G/T slice from one of the input rows) — the
//!     resulting consensus may have a lower median Tm once ambiguity
//!     codes resolve to A/T-rich expansions, which is accepted.
//!  3. (single-threaded) Pick the best candidate with a stable tie-break,
//!     remove its covered sequences from the pool, repeat.
//!
//! The tie-break (max coverage_count, then lowest start, lowest end, fewest
//! ambiguities, lex-smallest sequence) is deterministic and identical in
//! single-threaded and parallel execution — running with `--threads 1` and
//! `--threads 16` produces the same primer set.

use std::collections::{BTreeSet, HashMap, HashSet};

use rayon::prelude::*;

use crate::engine::iupac::{base_mask, is_ambiguous, mask_to_iupac, reverse_complement};
use crate::engine::tm::{calculate_tm, determine_oligo_length, TmParams};
use crate::engine::types::{
    Orientation, PrimerCandidate, PrimerSearchResult, Progress, SearchMode, SearchSettings,
};

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn find_primers(
    sequences: &[Vec<u8>],
    settings: &SearchSettings,
    injected: &[Vec<u8>],
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
    let tm_params = settings.tm_params();
    let is_reverse = matches!(settings.orientation, Orientation::Reverse);
    let seq_len = sequences[0].len();

    let mut remaining: BTreeSet<usize> = (0..total).collect();

    // Obligatory user-supplied oligos are placed first; the sequences they
    // cover are removed from the pool before the search begins.
    apply_injected_oligos(
        sequences,
        injected,
        &mut remaining,
        &mut result,
        &tm_params,
        is_reverse,
        progress,
    );

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
                    determine_oligo_length(seq, start, settings.tm_threshold, &tm_params)
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
                evaluate_range(start, end, &remaining_seqs, settings, &tm_params, is_reverse)
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

        let primer = make_primer(&best, total, is_reverse);

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
                best.count as f64 / total as f64 * 100.0,
                covered_so_far,
                total,
                covered_so_far as f64 / total as f64 * 100.0
            ),
            100.0,
        );
    }

    result
}

/// Fixed-slice analysis — the search is skipped entirely.
///
/// The whole input alignment is treated as a single, user-chosen slice
/// `[0, seq_len)`: the caller has already trimmed the FASTA down to exactly
/// the region they want the primer at. The engine simply generates the
/// oligo variant(s) required to cover *every* input sequence at that slice,
/// using the same greedy round loop as [`find_primers`] but with the single
/// fixed range instead of the dynamically-discovered ranges. Each round
/// emits the highest-coverage variant over the still-uncovered sequences and
/// removes those; the loop runs until coverage is 100%.
///
/// Variant generation obeys the configured [`SearchMode`] (exact vs. IUPAC
/// consensus) plus the IUPAC / 3'-conservation / orientation settings.
/// **The Tm threshold is not enforced**: the user picked the region, so
/// every variant needed for full coverage is emitted regardless of its Tm.
/// The per-variant Tm is still computed (from the configured concentrations)
/// and reported on each primer, so a cold slice is reported, not silently
/// dropped.
pub fn find_primers_fixed(
    sequences: &[Vec<u8>],
    settings: &SearchSettings,
    injected: &[Vec<u8>],
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
    let tm_params = settings.tm_params();
    let is_reverse = matches!(settings.orientation, Orientation::Reverse);
    let seq_len = sequences[0].len();
    let (start, end) = (0usize, seq_len);

    // The slice spans the whole alignment. Neutralise Tm gating inside the
    // variant evaluators (a fixed slice must yield whatever covers the input,
    // hot or cold). The real Tm is still computed via `tm_params`, which is
    // independent of `tm_threshold`.
    let mut eval_settings = settings.clone();
    eval_settings.tm_threshold = f64::NEG_INFINITY;

    let mut remaining: BTreeSet<usize> = (0..total).collect();

    // Obligatory user-supplied oligos are placed first; the sequences they
    // cover are removed before variants are generated for the fixed slice.
    apply_injected_oligos(
        sequences,
        injected,
        &mut remaining,
        &mut result,
        &tm_params,
        is_reverse,
        progress,
    );

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
                "Fixed slice: generating variant {} ({} sequences uncovered)",
                round,
                remaining_list.len()
            ),
            0.0,
        );

        // One range, evaluated directly — no phase-1 discovery, no rayon.
        let Some(best) = evaluate_range(
            start,
            end,
            &remaining_seqs,
            &eval_settings,
            &tm_params,
            is_reverse,
        ) else {
            // Unreachable for non-degenerate input (with Tm gating off the
            // evaluators always yield at least the most-frequent exact
            // slice), but guard rather than loop forever.
            result.message = "Could not generate a variant for the fixed slice".to_string();
            break;
        };

        let global_covered: BTreeSet<usize> = best
            .local_indices
            .iter()
            .map(|&j| remaining_list[j])
            .collect();

        let primer = make_primer(&best, total, is_reverse);
        for idx in &global_covered {
            remaining.remove(idx);
        }
        result.primers.push(primer);

        let covered_so_far = total - remaining.len();
        progress.report(
            &format!(
                "Fixed slice: variant {} covers {} sequences ({:.1}%). Total: {}/{} ({:.1}%)",
                round,
                best.count,
                best.count as f64 / total as f64 * 100.0,
                covered_so_far,
                total,
                covered_so_far as f64 / total as f64 * 100.0
            ),
            100.0,
        );
    }

    result
}

/// Build the displayed [`PrimerCandidate`] from a per-range evaluation
/// result: reverse-complement the sequence when the orientation is reverse
/// and fill in the coverage percentage. Shared by the search and fixed-slice
/// loops.
fn make_primer(best: &EvalResult, total: usize, is_reverse: bool) -> PrimerCandidate {
    let display_seq = if is_reverse {
        reverse_complement(&best.sequence)
    } else {
        best.sequence.clone()
    };
    PrimerCandidate {
        sequence: display_seq,
        coverage_count: best.count,
        coverage_pct: (best.count as f64 / total as f64) * 100.0,
        tm: best.tm,
        align_start: best.start,
        align_end: best.end,
        ambiguity_count: best.amb_count,
        injected: false,
    }
}

// ---------------------------------------------------------------------------
// Injected (obligatory) oligos
// ---------------------------------------------------------------------------

/// Place the user-supplied obligatory oligos before the search runs.
///
/// Each `oligo` is given in *display* orientation (exactly as the user typed
/// it). For a reverse-orientation run that means it is the reverse complement
/// of the alignment slice it binds, so it is reverse-complemented back into
/// forward alignment coordinates before matching — mirroring how the search
/// stores forward consensuses and only reverse-complements for display.
///
/// For each oligo we:
///   1. choose the alignment start that maximises coverage over the *full*
///      input set (tie-break: smallest start) — the oligo's intrinsic
///      best-fit position, independent of processing order;
///   2. emit it as a [`PrimerCandidate`] (`injected = true`) whose reported
///      coverage is the number of *still-uncovered* sequences it matches at
///      that position — the same greedy accounting the search loops use;
///   3. remove those sequences from `remaining`.
///
/// Injected oligos are obligatory: the Tm threshold and the
/// ambiguity/3'/IUPAC constraints are **not** enforced. Each oligo's Tm and
/// ambiguity count are still computed and reported. Oligos may contain IUPAC
/// ambiguity codes (e.g. an existing degenerate primer).
fn apply_injected_oligos(
    sequences: &[Vec<u8>],
    injected: &[Vec<u8>],
    remaining: &mut BTreeSet<usize>,
    result: &mut PrimerSearchResult,
    tm_params: &TmParams,
    is_reverse: bool,
    progress: &dyn Progress,
) {
    if injected.is_empty() {
        return;
    }
    let total = sequences.len();
    let seq_len = sequences.first().map(|s| s.len()).unwrap_or(0);

    for (i, oligo) in injected.iter().enumerate() {
        // Forward alignment form: undo the display-time reverse complement.
        let forward = if is_reverse {
            reverse_complement(oligo)
        } else {
            oligo.clone()
        };
        let l = forward.len();
        // An empty or over-long oligo cannot be placed at any position; the
        // CLI validates these away, so this is only a defensive guard.
        if l == 0 || l > seq_len {
            continue;
        }

        progress.report(
            &format!("Injecting oligo {}/{}", i + 1, injected.len()),
            0.0,
        );

        // 1. Intrinsic position: the start covering the most input sequences.
        let mut best_start = 0usize;
        let mut best_cov = 0usize;
        for start in 0..=(seq_len - l) {
            let cov = sequences
                .iter()
                .filter(|s| sequence_matches_consensus(&s[start..start + l], &forward))
                .count();
            if cov > best_cov {
                best_cov = cov;
                best_start = start;
            }
        }
        let (start, end) = (best_start, best_start + l);

        // 2. Credit only the still-uncovered sequences, and remove them.
        let covered: Vec<usize> = remaining
            .iter()
            .copied()
            .filter(|&j| sequence_matches_consensus(&sequences[j][start..end], &forward))
            .collect();
        for j in &covered {
            remaining.remove(j);
        }

        let amb_count = forward.iter().filter(|&&b| is_ambiguous(b)).count();
        let eval = EvalResult {
            start,
            end,
            sequence: forward,
            count: covered.len(),
            local_indices: covered,
            amb_count,
            tm: calculate_tm(oligo, tm_params),
        };
        let mut primer = make_primer(&eval, total, is_reverse);
        primer.injected = true;
        result.primers.push(primer);
    }
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
    tm_params: &TmParams,
    is_reverse: bool,
) -> Option<EvalResult> {
    let subsequences: Vec<&[u8]> = remaining_seqs.iter().map(|s| &s[start..end]).collect();

    let (variant, count, local_indices, amb_count) = match settings.mode {
        SearchMode::NoAmbiguities => {
            let (v, c, l) =
                eval_no_ambiguities(&subsequences, settings.tm_threshold, tm_params)?;
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
            settings.max_seeds,
            settings.tm_threshold,
            tm_params,
        )?,
    };

    // Tm of the actual displayed primer (the consensus / chosen variant).
    let tm = calculate_tm(&variant, tm_params);

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

/// No-ambiguity mode: the most frequent exact subsequence *whose own Tm
/// meets `tm_threshold`* wins. Ties broken by *first-appearance order* in
/// `subsequences` — matches the reference Python implementation, which
/// iterates a dict in insertion order with a `>` (strict-greater) update.
/// Returns `None` if no variant at this range meets the threshold.
fn eval_no_ambiguities(
    subsequences: &[&[u8]],
    tm_threshold: f64,
    tm_params: &TmParams,
) -> Option<(Vec<u8>, usize, Vec<usize>)> {
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
    let mut best_seq: &[u8] = b"";
    let mut best: Vec<usize> = Vec::new();
    let mut best_count: usize = 0;
    for (seq, idx) in entries.into_iter() {
        if calculate_tm(seq, tm_params) < tm_threshold {
            continue;
        }
        if idx.len() > best_count {
            best_count = idx.len();
            best_seq = seq;
            best = idx;
        }
    }
    if best_count == 0 {
        return None;
    }
    Some((best_seq.to_vec(), best_count, best))
}

/// Incremental mode: try increasing ambiguity budgets; for each budget try
/// up to `max_seeds` distinct subsequences as seeds (`0` = no cap). For
/// each seed, build a candidate group with greedy expansion in *several*
/// orderings and keep whichever produces the highest-coverage consensus.
/// Pick the best consensus across all seeds.
///
/// The orderings tried per seed:
///   1. First-appearance order (input order — matches the legacy Python
///      tool's behaviour).
///   2. Count descending, ties broken lexicographically (popular slices
///      absorbed first).
///   3. Count ascending, ties broken lexicographically (rare slices
///      absorbed first — leaves budget for high-frequency slices later).
///
/// Different orderings produce *structurally different* consensuses
/// (different IUPAC codes at different positions), not just tie-broken
/// equivalents. Trying multiple orderings and keeping the best is a cheap
/// way to approximate the per-range max-coverage objective more closely.
fn eval_incremental(
    subsequences: &[&[u8]],
    target_pct: f64,
    max_ambiguities: usize,
    exclude_n: bool,
    only_twofold: bool,
    three_prime_len: usize,
    is_reverse: bool,
    max_seeds: usize,
    tm_threshold: f64,
    tm_params: &TmParams,
) -> Option<(Vec<u8>, usize, Vec<usize>, usize)> {
    if subsequences.is_empty() {
        return None;
    }

    // Build `unique` in first-appearance order.
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
    let n_unique = unique.len();

    // Per-unique-variant Tm. A seed is a real A/C/G/T slice from one of
    // the input rows, so its own Tm must meet the user's threshold.
    // Consensuses grown from a passing seed are accepted even if
    // ambiguity expansion drops the median Tm below the threshold.
    let unique_tm: Vec<f64> = unique
        .iter()
        .map(|(seq, _)| calculate_tm(seq, tm_params))
        .collect();

    let total = subsequences.len();
    let target_count = (target_pct / 100.0 * total as f64).ceil() as usize;

    // Seed iteration order: count desc, ties broken by first-appearance.
    // Variants below the Tm threshold are filtered out before sorting so
    // they can never be chosen as a seed. `sort_by` is stable, so within
    // a count bucket first-appearance order is preserved.
    let mut seed_order: Vec<usize> = (0..n_unique)
        .filter(|&i| unique_tm[i] >= tm_threshold)
        .collect();
    seed_order.sort_by(|&a, &b| unique[b].1.len().cmp(&unique[a].1.len()));
    let seed_limit = if max_seeds == 0 {
        seed_order.len()
    } else {
        seed_order.len().min(max_seeds)
    };

    // Pre-compute the three expansion orderings — same for every seed.
    let order_first_app: Vec<usize> = (0..n_unique).collect();
    let mut order_count_desc: Vec<usize> = (0..n_unique).collect();
    order_count_desc.sort_by(|&a, &b| {
        unique[b]
            .1
            .len()
            .cmp(&unique[a].1.len())
            .then(unique[a].0.cmp(unique[b].0))
    });
    let mut order_count_asc: Vec<usize> = (0..n_unique).collect();
    order_count_asc.sort_by(|&a, &b| {
        unique[a]
            .1
            .len()
            .cmp(&unique[b].1.len())
            .then(unique[a].0.cmp(unique[b].0))
    });
    let orderings: [&[usize]; 3] = [
        &order_first_app,
        &order_count_desc,
        &order_count_asc,
    ];

    let mut best_consensus: Option<Vec<u8>> = None;
    let mut best_count: usize = 0;
    let mut best_amb: usize = usize::MAX;
    let mut best_local: Vec<usize> = Vec::new();
    let mut found_target = false;

    for amb_level in 0..=max_ambiguities {
        if found_target {
            break;
        }
        for &seed_idx in seed_order.iter().take(seed_limit) {
            // Try each ordering for this seed; remember the highest-
            // coverage valid consensus.
            let mut seed_best: Option<(Vec<u8>, usize, usize, Vec<usize>)> = None;

            for order in orderings.iter() {
                let Some((consensus, amb_count)) = greedy_expand(
                    seed_idx,
                    &unique,
                    order,
                    amb_level,
                    exclude_n,
                    only_twofold,
                ) else {
                    continue;
                };
                if !check_three_prime_ok(&consensus, three_prime_len, is_reverse) {
                    continue;
                }

                let mut coverage_count: usize = 0;
                let mut local_indices: Vec<usize> = Vec::new();
                for (uniq_seq, indices) in unique.iter() {
                    if sequence_matches_consensus(uniq_seq, &consensus) {
                        coverage_count += indices.len();
                        local_indices.extend_from_slice(indices);
                    }
                }

                let take = match &seed_best {
                    None => true,
                    Some((_, sb_amb, sb_count, _)) => {
                        coverage_count > *sb_count
                            || (coverage_count == *sb_count && amb_count < *sb_amb)
                    }
                };
                if take {
                    seed_best = Some((consensus, amb_count, coverage_count, local_indices));
                }
            }

            let Some((consensus, amb_count, coverage_count, local_indices)) = seed_best else {
                continue;
            };

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
        // Fallback: most-frequent threshold-passing exact subsequence.
        // `seed_order` is already filtered by Tm and sorted by count desc,
        // so its first element (if any) is the right pick. If no seed
        // passes, this range cannot yield a valid primer.
        let &i = seed_order.first()?;
        let seq = unique[i].0.to_vec();
        let indices = unique[i].1.clone();
        let count = indices.len();
        return Some((seq, count, indices, 0));
    }

    best_local.sort();
    Some((best_consensus.expect("checked above"), best_count, best_local, best_amb))
}

/// One greedy expansion pass: start from `seed_idx`, walk `unique` in the
/// given order, and absorb each candidate slice if the resulting consensus
/// stays within `amb_level` ambiguities and respects exclude_n /
/// only_twofold. Returns `(consensus, amb_count)` if the final consensus is
/// valid and within budget, `None` otherwise.
fn greedy_expand(
    seed_idx: usize,
    unique: &[(&[u8], Vec<usize>)],
    order: &[usize],
    amb_level: usize,
    exclude_n: bool,
    only_twofold: bool,
) -> Option<(Vec<u8>, usize)> {
    let mut group: Vec<&[u8]> = Vec::with_capacity(unique.len());
    group.push(unique[seed_idx].0);
    for &oi in order {
        if oi == seed_idx {
            continue;
        }
        group.push(unique[oi].0);
        let ok = match create_consensus(&group, exclude_n, only_twofold) {
            Some((_, amb_count)) => amb_count <= amb_level,
            None => false,
        };
        if !ok {
            group.pop();
        }
    }
    create_consensus(&group, exclude_n, only_twofold).filter(|(_, a)| *a <= amb_level)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::types::NoProgress;

    fn params() -> TmParams {
        TmParams {
            oligo_conc_um: 0.2,
            na_conc_mm: 50.0,
            mg_conc_mm: 3.0,
            dntp_conc_mm: 0.8,
        }
    }

    /// A `SearchSettings` for fixed-slice tests. Concentrations match the
    /// reference operating conditions; `target_coverage_pct = 100` so
    /// incremental mode absorbs as much as the ambiguity budget allows.
    fn fixed_settings(mode: SearchMode) -> SearchSettings {
        SearchSettings {
            tm_threshold: 60.0,
            oligo_concentration_um: 0.2,
            na_concentration_mm: 50.0,
            mg_concentration_mm: 3.0,
            dntp_concentration_mm: 0.8,
            mode,
            target_coverage_pct: 100.0,
            max_ambiguities: 2,
            exclude_n: false,
            only_twofold: false,
            orientation: Orientation::Forward,
            three_prime_match: 0,
            max_seeds: 0,
            fixed: true,
        }
    }

    /// The most-frequent variant at a range is AT-rich (sub-threshold); a
    /// minority variant is GC-rich and passes. The Tm filter must reject
    /// the popular-but-cold variant and pick the rare-but-hot one.
    #[test]
    fn no_ambiguities_skips_low_tm_variant() {
        let p = params();
        let low: &[u8] = b"AAATAAATAAATAAAT";
        let high: &[u8] = b"GCGCGCGCGCGCGCGC";
        assert!(calculate_tm(low, &p) < 60.0, "test setup: low Tm");
        assert!(calculate_tm(high, &p) >= 60.0, "test setup: high Tm");

        let subs: Vec<&[u8]> = vec![low, low, low, low, low, high, high];
        let (picked, count, _) =
            eval_no_ambiguities(&subs, 60.0, &p).expect("at least one variant passes");
        assert_eq!(picked, high.to_vec());
        assert_eq!(count, 2);
    }

    /// When no variant at the range meets threshold, the range yields no
    /// candidate — `evaluate_range` must propagate `None`.
    #[test]
    fn no_ambiguities_returns_none_when_all_below_threshold() {
        let p = params();
        let low: &[u8] = b"AAATAAATAAATAAAT";
        let subs: Vec<&[u8]> = vec![low, low, low];
        assert!(eval_no_ambiguities(&subs, 80.0, &p).is_none());
    }

    /// In incremental mode, a sub-threshold high-count seed must be
    /// skipped, leaving a passing minority seed to drive the consensus.
    /// With `max_ambiguities = 1` the two slices cannot be merged into one
    /// consensus, so the result is the high-Tm slice alone.
    #[test]
    fn incremental_skips_low_tm_seed() {
        let p = params();
        let low: &[u8] = b"AAATAAATAAATAAAT";
        let high: &[u8] = b"GCGCGCGCGCGCGCGC";
        let subs: Vec<&[u8]> = vec![low, low, low, low, low, high, high];
        let (consensus, count, _local, _amb) = eval_incremental(
            &subs, 100.0, 1, false, false, 0, false, 0, 60.0, &p,
        )
        .expect("high-Tm seed exists");
        assert_eq!(consensus, high.to_vec());
        assert_eq!(count, 2);
    }

    /// Incremental fallback: if no seed passes threshold at all, return
    /// `None` so the range is skipped rather than emitting a sub-threshold
    /// primer.
    #[test]
    fn incremental_returns_none_when_no_seed_passes() {
        let p = params();
        let low: &[u8] = b"AAATAAATAAATAAAT";
        let subs: Vec<&[u8]> = vec![low, low, low];
        assert!(
            eval_incremental(&subs, 100.0, 1, false, false, 0, false, 0, 80.0, &p).is_none()
        );
    }

    /// Fixed no-ambiguities: every distinct exact slice becomes its own
    /// variant, ordered by coverage descending, and together they cover the
    /// whole input. The reported range is always the full slice.
    #[test]
    fn fixed_no_ambiguities_emits_all_distinct_variants() {
        let settings = fixed_settings(SearchMode::NoAmbiguities);
        let a = b"GCGCGCGCGCGCGCGC".to_vec(); // count 3
        let b = b"ATATATATATATATAT".to_vec(); // count 2
        let c = b"TTTTGGGGCCCCAAAA".to_vec(); // count 1
        let seqs = vec![a.clone(), a.clone(), a.clone(), b.clone(), b.clone(), c.clone()];

        let res = find_primers_fixed(&seqs, &settings, &[], &NoProgress);

        assert_eq!(res.primers.len(), 3, "one variant per distinct slice");
        assert_eq!(res.primers[0].coverage_count, 3);
        assert_eq!(res.primers[1].coverage_count, 2);
        assert_eq!(res.primers[2].coverage_count, 1);
        for p in &res.primers {
            assert_eq!((p.align_start, p.align_end), (0, 16), "full fixed slice");
        }
        let covered: usize = res.primers.iter().map(|p| p.coverage_count).sum();
        assert_eq!(covered, seqs.len(), "100% coverage");
    }

    /// Fixed mode ignores the Tm threshold: a slice whose Tm is far below an
    /// (unreachable) threshold is still emitted. The equivalent search would
    /// find nothing here.
    #[test]
    fn fixed_ignores_tm_threshold() {
        let mut settings = fixed_settings(SearchMode::NoAmbiguities);
        settings.tm_threshold = 90.0; // unreachable for this AT-rich slice
        let low = b"AAATAAATAAATAAAT".to_vec();
        let seqs = vec![low.clone(), low.clone(), low.clone()];

        let res = find_primers_fixed(&seqs, &settings, &[], &NoProgress);

        assert_eq!(res.primers.len(), 1);
        assert_eq!(res.primers[0].coverage_count, 3);
        // The cold Tm is still reported, not hidden.
        assert!(res.primers[0].tm < 90.0);
    }

    /// Fixed incremental: slices differing at a single position are merged
    /// into one IUPAC consensus when the ambiguity budget allows, covering
    /// all input with a single variant.
    #[test]
    fn fixed_incremental_merges_into_consensus() {
        let mut settings = fixed_settings(SearchMode::Incremental);
        settings.max_ambiguities = 1;
        let a = b"ACGTACGTACGTACGA".to_vec();
        let b = b"ACGTACGTACGTACGG".to_vec(); // differs only at the last base
        let seqs = vec![a.clone(), a.clone(), b.clone()];

        let res = find_primers_fixed(&seqs, &settings, &[], &NoProgress);

        assert_eq!(res.primers.len(), 1, "merged into a single consensus");
        assert_eq!(res.primers[0].coverage_count, 3);
        assert_eq!(res.primers[0].ambiguity_count, 1);
        assert_eq!(res.primers[0].sequence, b"ACGTACGTACGTACGR".to_vec());
    }

    /// Fixed reverse orientation: the displayed sequence is the reverse
    /// complement of the alignment slice, while the reported range stays in
    /// forward alignment coordinates.
    #[test]
    fn fixed_reverse_orientation_reverse_complements() {
        let mut settings = fixed_settings(SearchMode::NoAmbiguities);
        settings.orientation = Orientation::Reverse;
        let s = b"AAAACCCCGGGGTTTT".to_vec();
        let seqs = vec![s.clone(), s.clone()];

        let res = find_primers_fixed(&seqs, &settings, &[], &NoProgress);

        assert_eq!(res.primers.len(), 1);
        assert_eq!(res.primers[0].sequence, reverse_complement(&s));
        assert_eq!((res.primers[0].align_start, res.primers[0].align_end), (0, 16));
    }

    /// A `SearchSettings` for the regular (non-fixed) search. `tm_threshold`
    /// is 0 so the search itself never gates anything out — these tests are
    /// about the injection step, not the Tm filter.
    fn search_settings(mode: SearchMode) -> SearchSettings {
        SearchSettings {
            tm_threshold: 0.0,
            oligo_concentration_um: 0.2,
            na_concentration_mm: 50.0,
            mg_concentration_mm: 3.0,
            dntp_concentration_mm: 0.8,
            mode,
            target_coverage_pct: 100.0,
            max_ambiguities: 2,
            exclude_n: false,
            only_twofold: false,
            orientation: Orientation::Forward,
            three_prime_match: 0,
            max_seeds: 0,
            fixed: false,
        }
    }

    /// Injected oligos are emitted first, in input order, each marked
    /// `injected`, and the sequences they cover are removed so the search
    /// adds nothing when injection already reaches 100%.
    #[test]
    fn injected_oligos_listed_first_and_remove_coverage() {
        let settings = search_settings(SearchMode::NoAmbiguities);
        let a = b"ACGTACGTACGTACGT".to_vec();
        let b = b"TTTTGGGGCCCCAAAA".to_vec();
        let seqs = vec![a.clone(), a.clone(), a.clone(), b.clone(), b.clone()];
        // Inject both exact slices, B before A, to check order is preserved.
        let injected = vec![b.clone(), a.clone()];

        let res = find_primers(&seqs, &settings, &injected, &NoProgress);

        assert_eq!(res.primers.len(), 2, "only the injected oligos were needed");
        assert!(res.primers.iter().all(|p| p.injected));
        assert_eq!(res.primers[0].sequence, b);
        assert_eq!(res.primers[0].coverage_count, 2);
        assert_eq!(res.primers[1].sequence, a);
        assert_eq!(res.primers[1].coverage_count, 3);
        let covered: usize = res.primers.iter().map(|p| p.coverage_count).sum();
        assert_eq!(covered, seqs.len());
    }

    /// An oligo shorter than the alignment is placed at the offset where it
    /// covers the most input sequences, and the reported range reflects that
    /// offset (1-based for display is handled by the formatter).
    #[test]
    fn injected_oligo_positioned_at_best_offset() {
        let settings = search_settings(SearchMode::NoAmbiguities);
        let s1 = b"AAAAGGGGTTTTCCCC".to_vec();
        let s2 = b"CCCCGGGGAAAATTTT".to_vec();
        let seqs = vec![s1, s2];
        let injected = vec![b"GGGG".to_vec()]; // common motif at offset 4

        let res = find_primers(&seqs, &settings, &injected, &NoProgress);

        let inj = &res.primers[0];
        assert!(inj.injected);
        assert_eq!((inj.align_start, inj.align_end), (4, 8));
        assert_eq!(inj.coverage_count, 2);
    }

    /// In a reverse-orientation run the user supplies the oligo in display
    /// (reverse-complement) form; it is matched in forward coordinates and
    /// echoed back in the form provided.
    #[test]
    fn injected_oligo_reverse_orientation_matches_rc() {
        let mut settings = search_settings(SearchMode::NoAmbiguities);
        settings.orientation = Orientation::Reverse;
        // Non-palindromic, so reverse_complement(s) != s — this genuinely
        // exercises the display-form -> forward-coordinate round trip.
        let s = b"AAAACCCCGGGGTTTA".to_vec();
        let seqs = vec![s.clone(), s.clone()];
        let rc = reverse_complement(&s);
        assert_ne!(rc, s, "test setup: sequence must not be self-complementary");
        let injected = vec![rc.clone()];

        let res = find_primers(&seqs, &settings, &injected, &NoProgress);

        assert_eq!(res.primers.len(), 1);
        let p = &res.primers[0];
        assert!(p.injected);
        assert_eq!(p.coverage_count, 2);
        assert_eq!(p.sequence, rc);
        assert_eq!((p.align_start, p.align_end), (0, 16));
    }

    /// An injected oligo may carry IUPAC ambiguity codes; it then covers
    /// every input variant the code admits, and its ambiguity count is
    /// reported.
    #[test]
    fn injected_ambiguous_oligo_covers_variants() {
        let settings = search_settings(SearchMode::NoAmbiguities);
        let a = b"ACGTACGTACGTACGA".to_vec();
        let g = b"ACGTACGTACGTACGG".to_vec(); // differs only at the last base
        let seqs = vec![a.clone(), a.clone(), g.clone()];
        let injected = vec![b"ACGTACGTACGTACGR".to_vec()]; // R = {A, G}

        let res = find_primers(&seqs, &settings, &injected, &NoProgress);

        assert_eq!(res.primers.len(), 1);
        let p = &res.primers[0];
        assert!(p.injected);
        assert_eq!(p.coverage_count, 3);
        assert_eq!(p.ambiguity_count, 1);
        assert_eq!(p.sequence, b"ACGTACGTACGTACGR".to_vec());
    }

    /// Injection runs first; the search then covers whatever the injected
    /// oligos left behind, and those primers are not marked injected.
    #[test]
    fn injection_then_search_covers_remainder() {
        let settings = search_settings(SearchMode::NoAmbiguities);
        let a = b"ACGTACGTACGTACGT".to_vec();
        let b = b"TTTTTTTTGGGGGGGG".to_vec();
        let seqs = vec![a.clone(), a.clone(), a.clone(), b.clone(), b.clone()];
        let injected = vec![a.clone()]; // covers the 3 A's; B's left for search

        let res = find_primers(&seqs, &settings, &injected, &NoProgress);

        assert!(res.primers[0].injected);
        assert_eq!(res.primers[0].coverage_count, 3);
        assert!(res.primers[1..].iter().all(|p| !p.injected));
        let covered: usize = res.primers.iter().map(|p| p.coverage_count).sum();
        assert_eq!(covered, seqs.len(), "100% coverage overall");
    }

    /// Injection composes with fixed-slice mode: injected oligos are placed
    /// first, then the fixed-slice loop generates variants for the rest.
    #[test]
    fn injected_oligo_in_fixed_mode() {
        let settings = fixed_settings(SearchMode::NoAmbiguities);
        let a = b"GCGCGCGCGCGCGCGC".to_vec();
        let b = b"ATATATATATATATAT".to_vec();
        let seqs = vec![a.clone(), a.clone(), b.clone()];
        let injected = vec![a.clone()];

        let res = find_primers_fixed(&seqs, &settings, &injected, &NoProgress);

        assert_eq!(res.primers.len(), 2);
        assert!(res.primers[0].injected);
        assert_eq!(res.primers[0].coverage_count, 2);
        assert!(!res.primers[1].injected);
        assert_eq!(res.primers[1].coverage_count, 1);
        assert_eq!(res.primers[1].sequence, b);
    }
}
