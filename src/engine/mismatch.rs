//! Optimize-by-mismatch: the best set of `n` oligos when every oligo carries
//! exactly `y` IUPAC ambiguity codes and target sequences may be bound with
//! mismatches.
//!
//! Unlike the greedy round loops in `search.rs`, this mode optimises the whole
//! set at once and never emits more than `n` oligos.
//!
//! Pipeline:
//!
//! 1. **Windows.** In fixed-slice mode the whole alignment is the only window.
//!    Otherwise the windows are the Tm-derived `(start, end)` ranges that the
//!    regular search discovers in its first phase, collected over every input
//!    sequence. Candidates from all windows compete in one pool, so the oligos
//!    of a set may sit at different positions.
//! 2. **Candidates.** At each window the IUPAC consensus is deconstructed:
//!    `y` of its variable positions keep their consensus code (or, where
//!    `--exclude-n` / `--only-twofold` forbid that code, the widest allowed
//!    sub-code) and every other variable position is resolved to one base
//!    observed there. A variant more than `x` mismatches away from every input
//!    slice covers nothing, so only variants within reach of at least one
//!    slice are enumerated — the outcome is the same as enumerating all of
//!    them. Candidates matching an excluded 3' signature and, outside fixed
//!    mode, candidates whose own Tm misses the threshold are dropped.
//! 3. **Reduction.** Each candidate is reduced to its coverage profile, the
//!    mismatch level of every sequence it covers. A candidate that another
//!    valid candidate matches or beats on every sequence can be swapped for it
//!    in any set without making the set worse, so it is dropped. The surviving
//!    pool does not depend on processing order.
//! 4. **Set search.** Exact depth-first branch-and-bound over sets of up to
//!    `n` pool candidates, pruned with submodular upper bounds. The first
//!    strictly-best set in a fixed visiting order wins, so results do not
//!    depend on the thread count.
//!
//! Scoring: `lower_or_equal` maximises the number of sequences whose best
//! match has at most `x` mismatches, ties broken by the fewest total
//! mismatches over those sequences. `exact` maximises the number of sequences
//! whose best match has exactly `x` mismatches.
//!
//! A sequence with a mismatch in the protected 3' region (`three_prime_match`
//! bases) counts as not covered, and candidates never place ambiguity codes
//! there.
//!
//! Both the enumeration and the set search are exponential in the worst case.
//! The enumeration is estimated before any work starts and the search counts
//! its work; exceeding either limit aborts with an error instead of reporting
//! a set that is not proven optimal.
//!
//! Profiles are bitsets with one plane per mismatch level, over distinct
//! window slices while generating and over distinct full-length sequences in
//! the pool, so unions, marginal gains and dominance tests are word-wide bit
//! operations weighted by sequence multiplicity.

use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use crate::engine::iupac::{base_mask, is_ambiguous, mask_to_iupac, reverse_complement};
use crate::engine::search::ExcludeSet;
use crate::engine::tm::{calculate_tm, determine_oligo_length, TmParams};
use crate::engine::types::{
    MismatchOp, MismatchReport, PrimerCandidate, PrimerSearchResult, Progress, SearchSettings,
};

const BASES: [u8; 4] = *b"ACGT";

/// Lists at least this long are scanned / evaluated in parallel.
const PAR_MIN_ITEMS: usize = 2048;

/// Progress is reported (and cancellation polled) every this many
/// evaluations during the set search.
const REPORT_EVERY: u64 = 1 << 22;

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Find the best set of up to `settings.mismatch.oligo_count` oligos under the
/// mismatch-tolerant coverage criterion of [`MismatchOp`].
///
/// `injected` oligos (display orientation) are obligatory members of the set
/// and count toward its size; each is placed at the alignment offset where it
/// alone scores best. `excluded` 3' signatures are applied to candidates in
/// both search and fixed-slice mode, since this mode chooses among candidates
/// rather than having to cover every sequence.
///
/// Returns an error when a work or memory limit is hit, or when the
/// parameters are inconsistent.
pub fn find_primers_by_mismatch(
    sequences: &[Vec<u8>],
    settings: &SearchSettings,
    injected: &[Vec<u8>],
    excluded: &[Vec<u8>],
    progress: &dyn Progress,
) -> Result<PrimerSearchResult, String> {
    let ms = &settings.mismatch;
    let total = sequences.len();
    let mut result = PrimerSearchResult {
        total_sequences: total,
        ..Default::default()
    };
    if total == 0 {
        result.message = "No sequences to analyze".to_string();
        return Ok(result);
    }
    let n = ms.oligo_count;
    if n == 0 {
        return Err("--n-oligos must be at least 1".to_string());
    }
    if injected.len() > n {
        return Err(format!(
            "{} oligo(s) were injected but the set holds only {} (--n-oligos); \
             injected oligos count toward the set",
            injected.len(),
            n
        ));
    }

    let x = ms.mismatches;
    let is_reverse = settings.is_reverse();
    let tm_params = settings.tm_params();
    let three_prime = settings.three_prime_match;
    let uni = Universe::new(sequences);
    let obj = Objective::new(ms.op, x, uni.total);
    let exclude = ExcludeSet::new(excluded);

    // 1. Windows.
    progress.report("Optimize by mismatch: collecting windows", 0.0);
    let ranges: Vec<(usize, usize)> = if settings.fixed {
        vec![(0, uni.len)]
    } else {
        discover_windows(&uni, settings.tm_threshold, &tm_params)
    };
    let windows: Vec<Window> = ranges
        .par_iter()
        .map(|&(s, e)| Window::build(&uni, s, e, three_prime, is_reverse))
        .collect();

    let gen_cfg = GenConfig {
        obj,
        y: ms.ambiguities,
        exclude_n: settings.exclude_n,
        only_twofold: settings.only_twofold,
        is_reverse,
        tm_gate: (!settings.fixed).then_some((settings.tm_threshold, tm_params)),
        exclude: &exclude,
    };

    // Refuse up front when the enumeration alone would be too large.
    let estimate: f64 = windows.iter().map(|w| gen_cfg.estimate(w)).sum();
    if ms.max_candidates > 0 && estimate > ms.max_candidates as f64 {
        return Err(format!(
            "the candidate space is too large: about {:.2e} candidate oligos would be \
             enumerated across {} window(s), above the limit of {} (--max-candidates; \
             0 = no limit). Reduce --mismatches and/or --ambiguities, run --fixed on a \
             narrower slice, or raise --max-candidates.",
            estimate,
            windows.len(),
            ms.max_candidates
        ));
    }

    // 2 + 3. Candidates per window, reduced within the window, then across
    // windows (profiles re-expressed over the full-length sequences).
    let done = AtomicUsize::new(0);
    let per_window: Vec<Result<WindowOutput, String>> = windows
        .par_iter()
        .enumerate()
        .map(|(wi, w)| {
            let out = gen_cfg.window_candidates(wi as u32, w, progress);
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            progress.report(
                &format!(
                    "Optimize by mismatch: generated candidates for {}/{} window(s)",
                    d,
                    windows.len()
                ),
                10.0 + 50.0 * d as f64 / windows.len() as f64,
            );
            out
        })
        .collect();

    let mut generated = 0u64;
    let mut window_cands: Vec<Vec<Cand>> = Vec::with_capacity(windows.len());
    for out in per_window {
        let out = out?;
        generated += out.generated;
        window_cands.push(out.cands);
    }

    progress.report("Optimize by mismatch: reducing the candidate pool", 60.0);
    let mut pool: Vec<Cand> = if window_cands.len() == 1 {
        // A single window's survivors are already mutually non-dominated.
        let cands = window_cands.pop().expect("one window");
        cands
            .into_par_iter()
            .map(|c| expand(c, &windows[0], &uni, obj.planes))
            .collect::<Result<_, _>>()?
    } else {
        let mut sky = Skyline::new(obj, uni.words);
        for (wi, cands) in window_cands.into_iter().enumerate() {
            if progress.cancelled() {
                return Err("search cancelled".to_string());
            }
            let expanded: Vec<Cand> = cands
                .into_par_iter()
                .map(|c| expand(c, &windows[wi], &uni, obj.planes))
                .collect::<Result<_, _>>()?;
            for c in expanded {
                sky.insert(c)?;
            }
        }
        sky.into_items()
    };
    pool.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.sec.cmp(&b.sec))
            .then(a.key.cmp(&b.key))
    });

    // Injected oligos: fixed members, placed at their best stand-alone offset.
    let n3 = three_prime;
    let mut members: Vec<Member> = Vec::new();
    let mut state = try_zeroed(obj.planes * uni.words, "allocating the search state")?;
    for oligo in injected {
        let forward = if is_reverse {
            reverse_complement(oligo)
        } else {
            oligo.clone()
        };
        let l = forward.len();
        // The CLI validates these away; this is only a defensive guard.
        if l == 0 || l > uni.len {
            continue;
        }
        let mut best_start = 0usize;
        let mut best_score = i64::MIN;
        for start in 0..=(uni.len - l) {
            let score: i64 = uni
                .seqs
                .iter()
                .zip(&uni.weight)
                .map(|(s, &w)| {
                    let lv = level_of(&forward, &s[start..start + l], x, n3, is_reverse);
                    w as i64 * obj.unit(lv)
                })
                .sum();
            if score > best_score {
                best_score = score;
                best_start = start;
            }
        }
        let levels: Vec<Option<usize>> = uni
            .seqs
            .iter()
            .map(|s| level_of(&forward, &s[best_start..best_start + l], x, n3, is_reverse))
            .collect();
        let planes = planes_from_levels(&levels, &uni, obj)?;
        for (s, p) in state.iter_mut().zip(&planes) {
            *s |= p;
        }
        members.push(Member {
            forward,
            display: oligo.clone(),
            start: best_start,
            end: best_start + l,
            injected: true,
            pool_index: usize::MAX,
        });
    }

    // 4. Set search.
    let slots = n - members.len();
    if slots > 0 && !pool.is_empty() {
        progress.report(
            &format!(
                "Optimize by mismatch: searching sets over {} candidates",
                pool.len()
            ),
            65.0,
        );
    }
    let limits = SearchLimits {
        max_work: ms.max_work,
        generated,
        set_size: n,
    };
    let (best_score, best_set, evaluations) =
        search_sets(obj, &uni.wt, uni.words, &pool, state, slots, &limits, progress)?;

    for &i in &best_set {
        let c = &pool[i];
        let forward = c.oligo.to_vec();
        let display = if is_reverse {
            reverse_complement(&forward)
        } else {
            forward.clone()
        };
        members.push(Member {
            forward,
            display,
            start: c.start as usize,
            end: c.end as usize,
            injected: false,
            pool_index: i,
        });
    }

    // Report: per-sequence best level, set breakdown, per-oligo credit.
    let levels: Vec<Vec<Option<usize>>> = members
        .iter()
        .map(|m| {
            uni.seqs
                .iter()
                .map(|s| level_of(&m.forward, &s[m.start..m.end], x, n3, is_reverse))
                .collect()
        })
        .collect();
    let best_level: Vec<Option<usize>> = (0..uni.seqs.len())
        .map(|u| levels.iter().filter_map(|l| l[u]).min())
        .collect();
    let counted: Vec<bool> = best_level
        .iter()
        .map(|&b| match ms.op {
            MismatchOp::LowerOrEqual => b.is_some(),
            MismatchOp::Exact => b == Some(x),
        })
        .collect();
    debug_assert_eq!(
        best_level
            .iter()
            .zip(&uni.weight)
            .map(|(&b, &w)| w as i64 * obj.unit(b))
            .sum::<i64>(),
        best_score,
        "search score must match the directly recomputed score"
    );

    let mut level_counts = vec![0usize; x + 1];
    let mut not_covered = 0usize;
    for (b, &w) in best_level.iter().zip(&uni.weight) {
        match b {
            Some(l) => level_counts[*l] += w as usize,
            None => not_covered += w as usize,
        }
    }

    // Each counted sequence is credited to its best-matching oligo, ties to
    // the oligo listed first, so per-oligo counts add up to the set's
    // coverage. Injected oligos are listed first (input order), then the
    // optimised ones by decreasing credit (ties: pool order).
    let is_best = |m: usize, u: usize| counted[u] && levels[m][u] == best_level[u];
    let mut order: Vec<usize> = (0..members.len()).filter(|&m| members[m].injected).collect();
    let mut credited = vec![false; uni.seqs.len()];
    let credit_of = |m: usize, credited: &[bool]| -> usize {
        (0..uni.seqs.len())
            .filter(|&u| !credited[u] && is_best(m, u))
            .map(|u| uni.weight[u] as usize)
            .sum()
    };
    let mut credit = vec![0usize; members.len()];
    for &m in &order {
        credit[m] = credit_of(m, &credited);
        for (u, c) in credited.iter_mut().enumerate() {
            if is_best(m, u) {
                *c = true;
            }
        }
    }
    let mut rest: Vec<usize> = (0..members.len()).filter(|&m| !members[m].injected).collect();
    while !rest.is_empty() {
        let (pos, gain) = rest
            .iter()
            .enumerate()
            .map(|(p, &m)| (p, credit_of(m, &credited)))
            .max_by(|a, b| {
                a.1.cmp(&b.1).then(
                    members[rest[b.0]]
                        .pool_index
                        .cmp(&members[rest[a.0]].pool_index),
                )
            })
            .expect("non-empty");
        let m = rest.remove(pos);
        credit[m] = gain;
        for (u, c) in credited.iter_mut().enumerate() {
            if is_best(m, u) {
                *c = true;
            }
        }
        order.push(m);
    }

    for &m in &order {
        let mb = &members[m];
        result.primers.push(PrimerCandidate {
            tm: calculate_tm(&mb.display, &tm_params),
            ambiguity_count: mb.forward.iter().filter(|&&b| is_ambiguous(b)).count(),
            sequence: mb.display.clone(),
            coverage_count: credit[m],
            coverage_pct: credit[m] as f64 / total as f64 * 100.0,
            align_start: mb.start,
            align_end: mb.end,
            injected: mb.injected,
        });
    }

    if result.primers.is_empty() {
        result.message = if pool.is_empty() {
            "No candidate oligo satisfies the constraints".to_string()
        } else {
            "Could not find any primers".to_string()
        };
    } else if result.primers.len() < n {
        result.message = if pool.is_empty() {
            "No candidate oligo satisfies the constraints".to_string()
        } else {
            format!(
                "{} of {} requested oligos used; further oligos would not improve the result",
                result.primers.len(),
                n
            )
        };
    }

    result.mismatch = Some(MismatchReport {
        counted: counted
            .iter()
            .zip(&uni.weight)
            .filter(|(c, _)| **c)
            .map(|(_, &w)| w as usize)
            .sum(),
        level_counts,
        not_covered,
        windows: windows.len(),
        candidates_generated: generated,
        candidates_reduced: pool.len(),
        evaluations,
    });
    Ok(result)
}

/// One oligo of the final set.
struct Member {
    forward: Vec<u8>,
    display: Vec<u8>,
    start: usize,
    end: usize,
    injected: bool,
    /// Index into the sorted pool; `usize::MAX` for injected oligos.
    pool_index: usize,
}

// ---------------------------------------------------------------------------
// Sequences, weights, objective
// ---------------------------------------------------------------------------

/// Distinct full-length input sequences with their multiplicities.
struct Universe<'u> {
    /// Distinct sequences in first-appearance order.
    seqs: Vec<&'u [u8]>,
    weight: Vec<u64>,
    /// Bit position of each sequence in pool-level bitsets (heaviest first,
    /// so that equal weights share words).
    bit: Vec<u32>,
    words: usize,
    wt: Weights,
    total: u64,
    /// Alignment length.
    len: usize,
}

impl<'u> Universe<'u> {
    fn new(sequences: &'u [Vec<u8>]) -> Self {
        let mut index: HashMap<&'u [u8], usize> = HashMap::new();
        let mut seqs: Vec<&'u [u8]> = Vec::new();
        let mut weight: Vec<u64> = Vec::new();
        for s in sequences {
            match index.get(s.as_slice()) {
                Some(&i) => weight[i] += 1,
                None => {
                    index.insert(s.as_slice(), seqs.len());
                    seqs.push(s.as_slice());
                    weight.push(1);
                }
            }
        }
        let bit = rank_by_weight(&weight);
        let mut by_bit = vec![0u64; seqs.len()];
        for (u, &b) in bit.iter().enumerate() {
            by_bit[b as usize] = weight[u];
        }
        Universe {
            words: words_for(seqs.len()),
            wt: Weights::new(&by_bit),
            total: sequences.len() as u64,
            len: sequences[0].len(),
            seqs,
            weight,
            bit,
        }
    }
}

fn words_for(n: usize) -> usize {
    n.div_ceil(64)
}

/// Rank of each element when sorted by weight descending, ties by index.
fn rank_by_weight(weight: &[u64]) -> Vec<u32> {
    let mut order: Vec<usize> = (0..weight.len()).collect();
    order.sort_by(|&a, &b| weight[b].cmp(&weight[a]).then(a.cmp(&b)));
    let mut rank = vec![0u32; weight.len()];
    for (r, &i) in order.iter().enumerate() {
        rank[i] = r as u32;
    }
    rank
}

/// Per-bit element weights for weighted popcounts. Elements are laid out
/// heaviest first, so most words hold elements of a single weight and are
/// counted with one popcount.
struct Weights {
    bit_w: Vec<u64>,
    /// The common weight of all elements in the word, or `u64::MAX` if mixed.
    uniform: Vec<u64>,
}

impl Weights {
    fn new(by_bit: &[u64]) -> Self {
        let n = by_bit.len();
        let words = words_for(n);
        let mut bit_w = vec![0u64; words * 64];
        bit_w[..n].copy_from_slice(by_bit);
        let uniform = (0..words)
            .map(|i| {
                let span = &bit_w[i * 64..((i + 1) * 64).min(n)];
                if span.iter().all(|&v| v == span[0]) {
                    span[0]
                } else {
                    u64::MAX
                }
            })
            .collect();
        Weights { bit_w, uniform }
    }

    /// Total weight of the set bits of word `i`.
    #[inline]
    fn word(&self, i: usize, w: u64) -> u64 {
        if w == 0 {
            return 0;
        }
        let u = self.uniform[i];
        if u != u64::MAX {
            return u * w.count_ones() as u64;
        }
        let base = i * 64;
        let mut s = 0;
        let mut w = w;
        while w != 0 {
            s += self.bit_w[base + w.trailing_zeros() as usize];
            w &= w - 1;
        }
        s
    }

    fn pop(&self, bits: &[u64]) -> u64 {
        bits.iter().enumerate().map(|(i, &w)| self.word(i, w)).sum()
    }
}

/// The optimisation objective and the layout of coverage profiles.
///
/// `LowerOrEqual` profiles have `x + 1` cumulative planes, plane `j` holding
/// the elements bound with at most `j` mismatches. A set scores
/// `big * |covered| + Σ_{j<x} |≤ j|`, which ranks by coverage within `x`
/// mismatches, then by fewest total mismatches (`Σ_{j<x} |≤ j|` equals
/// `Σ (x - level)` over covered elements). The score is monotone submodular.
///
/// `Exact` profiles have two planes: `X` (exactly `x` mismatches) and `Lo`
/// (fewer). A set counts the elements in `∪X` and not in `∪Lo`.
#[derive(Clone, Copy)]
struct Objective {
    op: MismatchOp,
    x: usize,
    planes: usize,
    big: i64,
}

impl Objective {
    fn new(op: MismatchOp, x: usize, total: u64) -> Self {
        let planes = match op {
            MismatchOp::LowerOrEqual => x + 1,
            MismatchOp::Exact => 2,
        };
        Objective {
            op,
            x,
            planes,
            big: x as i64 * total as i64 + 1,
        }
    }

    #[inline]
    fn coef(&self, j: usize) -> i64 {
        if j == self.x {
            self.big
        } else {
            1
        }
    }

    /// Contribution of one sequence (per unit weight) given its best level.
    fn unit(&self, level: Option<usize>) -> i64 {
        match (self.op, level) {
            (_, None) => 0,
            (MismatchOp::LowerOrEqual, Some(l)) => self.big + (self.x - l) as i64,
            (MismatchOp::Exact, Some(l)) => (l == self.x) as i64,
        }
    }

    /// Score of a set from the union of its profiles.
    fn union_score(&self, planes: &[u64], words: usize, wt: &Weights) -> i64 {
        match self.op {
            MismatchOp::LowerOrEqual => (0..self.planes)
                .map(|j| self.coef(j) * wt.pop(&planes[j * words..(j + 1) * words]) as i64)
                .sum(),
            MismatchOp::Exact => (0..words)
                .map(|i| wt.word(i, planes[i] & !planes[words + i]) as i64)
                .sum(),
        }
    }

    /// Stand-alone score of one candidate, plus a secondary sort key that
    /// places dominating candidates before the ones they dominate.
    fn standalone(&self, planes: &[u64], words: usize, wt: &Weights) -> (i64, i64) {
        match self.op {
            MismatchOp::LowerOrEqual => (self.union_score(planes, words, wt), 0),
            MismatchOp::Exact => (
                wt.pop(&planes[..words]) as i64,
                wt.pop(&planes[words..2 * words]) as i64,
            ),
        }
    }

    /// Whether `a` is at least as good as `b` in every set: `a` binds every
    /// element at least as well (`LowerOrEqual`), or `a`'s exact hits contain
    /// `b`'s and `a`'s below-`x` hits are contained in `b`'s (`Exact`).
    fn dominates(&self, a: &[u64], b: &[u64], words: usize) -> bool {
        match self.op {
            MismatchOp::LowerOrEqual => a.iter().zip(b).all(|(&a, &b)| b & !a == 0),
            MismatchOp::Exact => {
                let (ax, alo) = a.split_at(words);
                let (bx, blo) = b.split_at(words);
                ax.iter().zip(bx).all(|(&a, &b)| b & !a == 0)
                    && alo.iter().zip(blo).all(|(&a, &b)| a & !b == 0)
            }
        }
    }
}

/// Mismatches between `oligo` (forward-orientation IUPAC) and the equally
/// long `slice`, or `None` when they exceed `x` or a mismatch falls in the
/// protected 3' region of `n3` bases.
fn level_of(oligo: &[u8], slice: &[u8], x: usize, n3: usize, is_reverse: bool) -> Option<usize> {
    let l = oligo.len();
    let n3 = n3.min(l);
    let mut mm = 0;
    for (k, (&o, &s)) in oligo.iter().zip(slice).enumerate() {
        if base_mask(o) & base_mask(s) == 0 {
            let in3 = if is_reverse { k < n3 } else { k >= l - n3 };
            if in3 {
                return None;
            }
            mm += 1;
            if mm > x {
                return None;
            }
        }
    }
    Some(mm)
}

/// Pool-level profile from per-sequence levels.
fn planes_from_levels(
    levels: &[Option<usize>],
    uni: &Universe,
    obj: Objective,
) -> Result<Vec<u64>, String> {
    let words = uni.words;
    let mut planes = try_zeroed(obj.planes * words, "allocating an oligo profile")?;
    for (u, lv) in levels.iter().enumerate() {
        let Some(l) = *lv else { continue };
        let b = uni.bit[u] as usize;
        let (w, m) = (b / 64, 1u64 << (b % 64));
        match obj.op {
            MismatchOp::LowerOrEqual => {
                for j in l..=obj.x {
                    planes[j * words + w] |= m;
                }
            }
            MismatchOp::Exact => {
                let plane = if l == obj.x { 0 } else { 1 };
                planes[plane * words + w] |= m;
            }
        }
    }
    Ok(planes)
}

// ---------------------------------------------------------------------------
// Windows and candidate generation
// ---------------------------------------------------------------------------

/// Tm-derived windows, exactly as phase 1 of the regular search collects them
/// (distinct sequences in first-appearance order yield the same windows in
/// the same order as all sequences would).
fn discover_windows(uni: &Universe, tm_threshold: f64, tm: &TmParams) -> Vec<(usize, usize)> {
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    let mut out = Vec::new();
    for seq in &uni.seqs {
        for start in 0..uni.len {
            if let Some(len) = determine_oligo_length(seq, start, tm_threshold, tm) {
                let end = start + len;
                if end <= uni.len && seen.insert((start, end)) {
                    out.push((start, end));
                }
            }
        }
    }
    out
}

/// One alignment window: its distinct slices ("groups") and, for every
/// variable position, which groups carry which base.
struct Window<'u> {
    start: usize,
    end: usize,
    /// Distinct slices, heaviest first; a group's index is its bit position.
    slices: Vec<&'u [u8]>,
    /// Distinct full-length sequences (universe indices) behind each group.
    members: Vec<Vec<u32>>,
    words: usize,
    wt: Weights,
    valid: Vec<u64>,
    /// Window offsets of the variable positions.
    var: Vec<usize>,
    /// Observed base mask at each variable position.
    obs: Vec<u8>,
    /// Whether each variable position lies in the protected 3' region.
    in3: Vec<bool>,
    /// Group bitset per (variable position, base index):
    /// `gbits[(vi * 4 + b) * words..][..words]`.
    gbits: Vec<u64>,
    /// Base index of group `g` at variable position `vi`:
    /// `gbase[g * var.len() + vi]`.
    gbase: Vec<u8>,
}

#[inline]
fn base_index(b: u8) -> u8 {
    // Input is quality-filtered to A/C/G/T.
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    }
}

impl<'u> Window<'u> {
    fn build(uni: &Universe<'u>, start: usize, end: usize, three_prime: usize, is_reverse: bool) -> Self {
        let len = end - start;
        let mut index: HashMap<&'u [u8], usize> = HashMap::new();
        let mut first_slices: Vec<&'u [u8]> = Vec::new();
        let mut first_members: Vec<Vec<u32>> = Vec::new();
        let mut first_weight: Vec<u64> = Vec::new();
        for (u, s) in uni.seqs.iter().enumerate() {
            let sl: &'u [u8] = &s[start..end];
            match index.get(sl) {
                Some(&g) => {
                    first_members[g].push(u as u32);
                    first_weight[g] += uni.weight[u];
                }
                None => {
                    index.insert(sl, first_slices.len());
                    first_slices.push(sl);
                    first_members.push(vec![u as u32]);
                    first_weight.push(uni.weight[u]);
                }
            }
        }
        let ng = first_slices.len();
        let rank = rank_by_weight(&first_weight);
        let mut order = vec![0usize; ng];
        for (g, &r) in rank.iter().enumerate() {
            order[r as usize] = g;
        }
        let slices: Vec<&'u [u8]> = order.iter().map(|&g| first_slices[g]).collect();
        let members: Vec<Vec<u32>> = order
            .iter()
            .map(|&g| std::mem::take(&mut first_members[g]))
            .collect();
        let weight: Vec<u64> = order.iter().map(|&g| first_weight[g]).collect();

        let words = words_for(ng);
        let mut valid = vec![!0u64; words];
        if !ng.is_multiple_of(64) {
            valid[words - 1] = (1u64 << (ng % 64)) - 1;
        }

        let mut var = Vec::new();
        let mut obs = Vec::new();
        for off in 0..len {
            let m = slices.iter().fold(0u8, |m, s| m | base_mask(s[off]));
            if m.count_ones() >= 2 {
                var.push(off);
                obs.push(m);
            }
        }
        let n3 = three_prime.min(len);
        let in3: Vec<bool> = var
            .iter()
            .map(|&off| if is_reverse { off < n3 } else { off >= len - n3 })
            .collect();

        let nvar = var.len();
        let mut gbits = vec![0u64; nvar * 4 * words];
        let mut gbase = vec![0u8; ng * nvar];
        for (g, s) in slices.iter().enumerate() {
            for (vi, &off) in var.iter().enumerate() {
                let b = base_index(s[off]);
                gbase[g * nvar + vi] = b;
                gbits[(vi * 4 + b as usize) * words + g / 64] |= 1u64 << (g % 64);
            }
        }

        Window {
            start,
            end,
            wt: Weights::new(&weight),
            slices,
            members,
            words,
            valid,
            var,
            obs,
            in3,
            gbits,
            gbase,
        }
    }

    #[inline]
    fn group_bits(&self, vi: usize, b: usize) -> &[u64] {
        let o = (vi * 4 + b) * self.words;
        &self.gbits[o..o + self.words]
    }
}

/// Ambiguity codes allowed at a variable position with observed bases `obs`:
/// the full consensus code, or — where `--exclude-n` / `--only-twofold`
/// forbid it — the widest allowed sub-codes. Narrower codes are never
/// needed: a wider code binds a superset of the sequences.
fn allowed_codes(obs: u8, exclude_n: bool, only_twofold: bool) -> Vec<u8> {
    let forbidden = |m: u8| (exclude_n && m == 0b1111) || (only_twofold && m.count_ones() > 2);
    if !forbidden(obs) {
        return vec![obs];
    }
    let subs: Vec<u8> = (1u8..16)
        .filter(|&m| m & !obs == 0 && m.count_ones() >= 2 && !forbidden(m))
        .collect();
    let widest = subs.iter().map(|m| m.count_ones()).max().unwrap_or(0);
    subs.into_iter().filter(|m| m.count_ones() == widest).collect()
}

/// All `k`-subsets of `0..n` in lexicographic order.
fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    if k > n {
        return out;
    }
    let mut c: Vec<usize> = (0..k).collect();
    loop {
        out.push(c.clone());
        let mut i = k;
        loop {
            if i == 0 {
                return out;
            }
            i -= 1;
            if c[i] != i + n - k {
                break;
            }
        }
        c[i] += 1;
        for j in i + 1..k {
            c[j] = c[j - 1] + 1;
        }
    }
}

/// A candidate oligo and its coverage profile.
struct Cand {
    /// Generation order `(window, ambiguity-position set, local)`; the
    /// deterministic final tie-break.
    key: (u32, u32, u32),
    score: i64,
    sec: i64,
    /// `obj.planes` bitsets, `words` each, over the window's groups during
    /// generation and over the distinct sequences once in the pool.
    planes: Box<[u64]>,
    /// Forward-orientation IUPAC sequence.
    oligo: Box<[u8]>,
    start: u32,
    end: u32,
}

struct WindowOutput {
    cands: Vec<Cand>,
    generated: u64,
}

struct GenConfig<'a> {
    obj: Objective,
    y: usize,
    exclude_n: bool,
    only_twofold: bool,
    is_reverse: bool,
    /// `(threshold, params)` when each candidate's own Tm must reach the
    /// threshold (search mode); `None` in fixed-slice mode.
    tm_gate: Option<(f64, TmParams)>,
    exclude: &'a ExcludeSet,
}

impl GenConfig<'_> {
    fn eligible(&self, w: &Window) -> Vec<usize> {
        (0..w.var.len()).filter(|&vi| !w.in3[vi]).collect()
    }

    /// Upper bound on the enumeration work at a window: per group and per
    /// ambiguity layout, every resolution within `x` changes. Computed with
    /// a generating polynomial over the eligible positions (`z` marks an
    /// ambiguity position, `t` a changed base).
    fn estimate(&self, w: &Window) -> f64 {
        let x = self.obj.x;
        let eligible = self.eligible(w);
        let y = self.y.min(eligible.len());
        let mut poly = vec![vec![0f64; x + 1]; y + 1];
        poly[0][0] = 1.0;
        for &vi in &eligible {
            let alt = (w.obs[vi].count_ones() - 1) as f64;
            let codes = allowed_codes(w.obs[vi], self.exclude_n, self.only_twofold).len() as f64;
            let prev = poly.clone();
            for z in 0..=y {
                for t in 0..=x {
                    let mut v = prev[z][t];
                    if t > 0 {
                        v += prev[z][t - 1] * alt;
                    }
                    if z > 0 {
                        v += prev[z - 1][t] * codes;
                    }
                    poly[z][t] = v;
                }
            }
        }
        let layouts = poly[y][0];
        let per_group = match self.obj.op {
            MismatchOp::LowerOrEqual => poly[y].iter().sum(),
            MismatchOp::Exact if x == 0 => layouts,
            MismatchOp::Exact => poly[y][x] + layouts,
        };
        per_group * w.slices.len() as f64
    }

    /// All valid, mutually non-dominated candidates of one window.
    fn window_candidates(
        &self,
        wi: u32,
        w: &Window,
        progress: &dyn Progress,
    ) -> Result<WindowOutput, String> {
        let eligible = self.eligible(w);
        let y = self.y.min(eligible.len());
        let codes: Vec<Vec<u8>> = eligible
            .iter()
            .map(|&vi| allowed_codes(w.obs[vi], self.exclude_n, self.only_twofold))
            .collect();
        let combos = combinations(eligible.len(), y);
        let parts: Vec<Result<(Vec<Cand>, u64), String>> = combos
            .par_iter()
            .enumerate()
            .map(|(pi, combo)| {
                let amb: Vec<usize> = combo.iter().map(|&e| eligible[e]).collect();
                let amb_codes: Vec<&[u8]> = combo.iter().map(|&e| codes[e].as_slice()).collect();
                self.layout_candidates(wi, pi as u32, w, &amb, &amb_codes, progress)
            })
            .collect();

        let mut generated = 0u64;
        let mut lists = Vec::with_capacity(parts.len());
        for p in parts {
            let (cands, g) = p?;
            generated += g;
            lists.push(cands);
        }
        let cands = if lists.len() == 1 {
            lists.pop().expect("one list")
        } else {
            let mut sky = Skyline::new(self.obj, w.words);
            for list in lists {
                for c in list {
                    sky.insert(c)?;
                }
            }
            sky.into_items()
        };
        Ok(WindowOutput { cands, generated })
    }

    /// Candidates of one window with ambiguity codes at the variable
    /// positions `amb`, for every combination of their allowed codes.
    fn layout_candidates(
        &self,
        wi: u32,
        pi: u32,
        w: &Window,
        amb: &[usize],
        amb_codes: &[&[u8]],
        progress: &dyn Progress,
    ) -> Result<(Vec<Cand>, u64), String> {
        let nvar = w.var.len();
        let mut is_amb = vec![false; nvar];
        for &vi in amb {
            is_amb[vi] = true;
        }
        let resolved: Vec<usize> = (0..nvar).filter(|&vi| !is_amb[vi]).collect();
        let changeable: Vec<usize> = (0..resolved.len())
            .filter(|&r| !w.in3[resolved[r]])
            .collect();
        let x = self.obj.x;
        let mut st = LayoutState {
            cfg: self,
            w,
            wi,
            pi,
            amb,
            resolved: &resolved,
            changeable: &changeable,
            codes: vec![0u8; amb.len()],
            amb_match: vec![0u64; amb.len() * w.words],
            template: w.slices[0].to_vec(),
            seen: KeySet::default(),
            projected: KeySet::default(),
            therm: vec![0u64; (x + 1) * w.words],
            uncov: vec![0u64; w.words],
            buf: vec![0u64; self.obj.planes * w.words],
            sky: Skyline::new(self.obj, w.words),
            local: 0,
            generated: 0,
            error: None,
        };

        // Odometer over the allowed codes, last position fastest.
        let mut idx = vec![0usize; amb.len()];
        'layouts: loop {
            if progress.cancelled() {
                return Err("search cancelled".to_string());
            }
            for k in 0..amb.len() {
                st.codes[k] = amb_codes[k][idx[k]];
            }
            st.run()?;
            let mut k = amb.len();
            loop {
                if k == 0 {
                    break 'layouts;
                }
                k -= 1;
                idx[k] += 1;
                if idx[k] < amb_codes[k].len() {
                    break;
                }
                idx[k] = 0;
            }
        }
        Ok((st.sky.into_items(), st.generated))
    }
}

/// Enumeration state for one ambiguity layout (positions + codes).
struct LayoutState<'a> {
    cfg: &'a GenConfig<'a>,
    w: &'a Window<'a>,
    wi: u32,
    pi: u32,
    /// Variable-position indices carrying an ambiguity code.
    amb: &'a [usize],
    /// Variable-position indices resolved to a single base.
    resolved: &'a [usize],
    /// Indices into `resolved` that may differ from the source slice (all
    /// but the protected 3' region).
    changeable: &'a [usize],
    codes: Vec<u8>,
    /// Groups matching each ambiguity code: `amb.len()` bitsets.
    amb_match: Vec<u64>,
    template: Vec<u8>,
    /// Resolutions already visited for the current codes.
    seen: KeySet,
    /// Source projections (with their budget) already expanded.
    projected: KeySet,
    therm: Vec<u64>,
    uncov: Vec<u64>,
    /// Profile of the resolution being visited.
    buf: Vec<u64>,
    sky: Skyline,
    local: u32,
    generated: u64,
    error: Option<String>,
}

impl LayoutState<'_> {
    /// Enumerate every resolution within reach of some group for the current
    /// codes, inserting the valid ones into the skyline.
    fn run(&mut self) -> Result<(), String> {
        let w = self.w;
        let words = w.words;
        let nvar = w.var.len();
        for (k, &vi) in self.amb.iter().enumerate() {
            let dst = &mut self.amb_match[k * words..(k + 1) * words];
            dst.fill(0);
            for b in 0..4 {
                if self.codes[k] & (1 << b) != 0 {
                    for (d, s) in dst.iter_mut().zip(w.group_bits(vi, b)) {
                        *d |= s;
                    }
                }
            }
        }
        self.seen.clear();
        self.projected.clear();
        let x = self.cfg.obj.x;
        let mut key: Vec<u8> = vec![0; self.resolved.len()];
        for g in 0..w.slices.len() {
            let gb = &w.gbase[g * nvar..(g + 1) * nvar];
            let m0 = self
                .amb
                .iter()
                .enumerate()
                .filter(|&(k, &vi)| self.codes[k] & (1 << gb[vi]) == 0)
                .count();
            if m0 > x {
                continue;
            }
            let budget = x - m0;
            for (k, &vi) in key.iter_mut().zip(self.resolved) {
                *k = gb[vi];
            }
            if !self.projected.insert(&key, budget) {
                continue;
            }
            self.neighbours(&mut key, 0, budget);
            if let Some(e) = self.error.take() {
                return Err(e);
            }
        }
        Ok(())
    }

    /// Visit `key` and every variant with up to `budget` (exactly `budget`
    /// in `Exact` mode) further changes at changeable positions `from..`.
    fn neighbours(&mut self, key: &mut [u8], from: usize, budget: usize) {
        if self.error.is_some() {
            return;
        }
        if budget == 0 || self.cfg.obj.op == MismatchOp::LowerOrEqual {
            self.visit(key);
        }
        if budget == 0 {
            return;
        }
        for ci in from..self.changeable.len() {
            let r = self.changeable[ci];
            let orig = key[r];
            let obs = self.w.obs[self.resolved[r]];
            for b in 0..4u8 {
                if b != orig && obs & (1 << b) != 0 {
                    key[r] = b;
                    self.neighbours(key, ci + 1, budget - 1);
                }
            }
            key[r] = orig;
        }
    }

    fn visit(&mut self, key: &[u8]) {
        if !self.seen.insert(key, 0) {
            return;
        }
        self.generated += 1;
        self.profile(key);
        let w = self.w;
        let (obj, words) = (self.cfg.obj, w.words);
        let covered = match obj.op {
            MismatchOp::LowerOrEqual => &self.buf[obj.x * words..],
            MismatchOp::Exact => &self.buf[..words],
        };
        if covered.iter().all(|&b| b == 0) {
            return;
        }
        let (score, sec) = obj.standalone(&self.buf, words, &w.wt);
        // Nearly every resolution is dominated; only survivors pay for the
        // oligo, the exclusion and the Tm check. (A dominated candidate is
        // dropped regardless of its validity: every skyline item is valid.)
        let Some((hash, kill)) = self.sky.check(&self.buf, score) else {
            return;
        };
        let mut oligo = self.template.clone();
        for (k, &vi) in self.amb.iter().enumerate() {
            oligo[w.var[vi]] = mask_to_iupac(self.codes[k]);
        }
        for (r, &vi) in self.resolved.iter().enumerate() {
            oligo[w.var[vi]] = BASES[key[r] as usize];
        }
        if self.cfg.exclude.excludes(&oligo, self.cfg.is_reverse) {
            return;
        }
        if let Some((threshold, params)) = self.cfg.tm_gate
            && calculate_tm(&oligo, &params) < threshold
        {
            return;
        }
        let mut planes = Vec::new();
        if planes.try_reserve_exact(self.buf.len()).is_err() {
            self.error = Some(out_of_memory("building candidate profiles", self.buf.len() * 8));
            return;
        }
        planes.extend_from_slice(&self.buf);
        let cand = Cand {
            key: (self.wi, self.pi, self.local),
            score,
            sec,
            planes: planes.into_boxed_slice(),
            oligo: oligo.into_boxed_slice(),
            start: w.start as u32,
            end: w.end as u32,
        };
        self.local += 1;
        if let Err(e) = self.sky.commit(cand, hash, kill) {
            self.error = Some(e);
        }
    }

    /// Group-level profile of the current codes with resolution `key`, into
    /// `buf`: mismatches per group counted in thermometer planes (`therm[j]`
    /// = groups with more than `j` mismatches), then turned into the
    /// objective's planes.
    fn profile(&mut self, key: &[u8]) {
        let w = self.w;
        let words = w.words;
        let x = self.cfg.obj.x;
        self.therm.fill(0);
        self.uncov.fill(0);
        for k in 0..self.amb.len() {
            let m = &self.amb_match[k * words..(k + 1) * words];
            add_mismatches(&mut self.therm, &mut self.uncov, m, &w.valid, x, false);
        }
        for (r, &vi) in self.resolved.iter().enumerate() {
            let m = w.group_bits(vi, key[r] as usize);
            add_mismatches(&mut self.therm, &mut self.uncov, m, &w.valid, x, w.in3[vi]);
        }
        let (therm, uncov, buf) = (&self.therm, &self.uncov, &mut self.buf);
        let le = |j: usize, i: usize| !therm[j * words + i] & !uncov[i] & w.valid[i];
        match self.cfg.obj.op {
            MismatchOp::LowerOrEqual => {
                for j in 0..=x {
                    for i in 0..words {
                        buf[j * words + i] = le(j, i);
                    }
                }
            }
            MismatchOp::Exact => {
                for i in 0..words {
                    let below = if x > 0 { le(x - 1, i) } else { 0 };
                    buf[i] = le(x, i) & !below;
                    buf[words + i] = below;
                }
            }
        }
    }
}

/// Set of resolutions (base indices 0..4, all of one length) with a small
/// tag, packed two bits per base when they fit in 128 bits.
#[derive(Default)]
struct KeySet {
    short: HashSet<(u128, usize)>,
    long: HashSet<(Vec<u8>, usize)>,
}

impl KeySet {
    /// Insert `(key, tag)`; `true` if it was not present.
    fn insert(&mut self, key: &[u8], tag: usize) -> bool {
        if key.len() <= 64 {
            let packed = key.iter().fold(0u128, |v, &b| (v << 2) | b as u128);
            self.short.insert((packed, tag))
        } else {
            self.long.insert((key.to_vec(), tag))
        }
    }

    fn clear(&mut self) {
        self.short.clear();
        self.long.clear();
    }
}

/// Count one position: groups outside `matched` gain a mismatch, or become
/// uncovered if the position is in the protected 3' region.
#[inline]
fn add_mismatches(
    therm: &mut [u64],
    uncov: &mut [u64],
    matched: &[u64],
    valid: &[u64],
    x: usize,
    in3: bool,
) {
    let words = valid.len();
    for i in 0..words {
        let mism = !matched[i] & valid[i];
        if mism == 0 {
            continue;
        }
        if in3 {
            uncov[i] |= mism;
            continue;
        }
        for j in (1..=x).rev() {
            therm[j * words + i] |= therm[(j - 1) * words + i] & mism;
        }
        therm[i] |= mism;
    }
}

/// Re-express a window-level profile over the distinct full-length
/// sequences.
fn expand(c: Cand, w: &Window, uni: &Universe, planes: usize) -> Result<Cand, String> {
    let words = uni.words;
    let mut out = try_zeroed(planes * words, "building the candidate pool")?;
    for j in 0..planes {
        let src = &c.planes[j * w.words..(j + 1) * w.words];
        let dst = &mut out[j * words..(j + 1) * words];
        for (i, &word) in src.iter().enumerate() {
            let mut bits = word;
            while bits != 0 {
                let g = i * 64 + bits.trailing_zeros() as usize;
                for &u in &w.members[g] {
                    let b = uni.bit[u as usize] as usize;
                    dst[b / 64] |= 1u64 << (b % 64);
                }
                bits &= bits - 1;
            }
        }
    }
    Ok(Cand {
        planes: out.into_boxed_slice(),
        ..c
    })
}

// ---------------------------------------------------------------------------
// Dominance reduction
// ---------------------------------------------------------------------------

/// The non-dominated candidates seen so far. Insertion keeps the set an
/// antichain: a newcomer is discarded if an item dominates it (identical
/// profiles keep the earlier item), otherwise it displaces every item it
/// dominates. Only valid candidates are committed. The final set is the
/// non-dominated subset of everything inserted, independent of insertion order up to identical profiles, where
/// the first inserted wins.
struct Skyline {
    obj: Objective,
    words: usize,
    items: Vec<Cand>,
    alive: Vec<bool>,
    n_alive: usize,
    by_hash: HashMap<u64, Vec<u32>>,
}

impl Skyline {
    fn new(obj: Objective, words: usize) -> Self {
        Skyline {
            obj,
            words,
            items: Vec::new(),
            alive: Vec::new(),
            n_alive: 0,
            by_hash: HashMap::new(),
        }
    }

    /// Insert `c` unless an item matches or beats it everywhere.
    fn insert(&mut self, c: Cand) -> Result<(), String> {
        match self.check(&c.planes, c.score) {
            Some((hash, kill)) => self.commit(c, hash, kill),
            None => Ok(()),
        }
    }

    /// `None` when an item has an identical profile or dominates `planes`;
    /// otherwise the profile hash and the items `planes` dominates, to pass
    /// to [`Skyline::commit`] if the candidate turns out valid.
    fn check(&self, planes: &[u64], score: i64) -> Option<(u64, Vec<usize>)> {
        let h = hash_words(planes);
        if let Some(ids) = self.by_hash.get(&h)
            && ids
                .iter()
                .any(|&i| self.alive[i as usize] && *self.items[i as usize].planes == *planes)
        {
            return None;
        }
        let (obj, words) = (self.obj, self.words);
        let beaten_by = |k: &Cand| k.score >= score && obj.dominates(&k.planes, planes, words);
        let beats = |k: &Cand| score >= k.score && obj.dominates(planes, &k.planes, words);
        // The items form an antichain, so the newcomer either is dominated by
        // some item or dominates some items, never both.
        if self.items.len() >= PAR_MIN_ITEMS {
            let (items, alive) = (&self.items, &self.alive);
            if (0..items.len())
                .into_par_iter()
                .any(|i| alive[i] && beaten_by(&items[i]))
            {
                return None;
            }
            let kill = (0..items.len())
                .into_par_iter()
                .filter(|&i| alive[i] && beats(&items[i]))
                .collect();
            Some((h, kill))
        } else {
            let mut kill = Vec::new();
            for (i, k) in self.items.iter().enumerate() {
                if !self.alive[i] {
                    continue;
                }
                if beaten_by(k) {
                    return None;
                }
                if beats(k) {
                    kill.push(i);
                }
            }
            Some((h, kill))
        }
    }

    /// Add `c`, displacing `kill`, as returned by [`Skyline::check`] for its
    /// profile with no insertion in between.
    fn commit(&mut self, c: Cand, hash: u64, kill: Vec<usize>) -> Result<(), String> {
        for i in kill {
            self.alive[i] = false;
            self.n_alive -= 1;
        }
        self.items
            .try_reserve(1)
            .map_err(|_| out_of_memory("reducing the candidate pool", 0))?;
        self.by_hash
            .entry(hash)
            .or_default()
            .push(self.items.len() as u32);
        self.items.push(c);
        self.alive.push(true);
        self.n_alive += 1;
        let dead = self.items.len() - self.n_alive;
        if dead > 4096 && dead > self.n_alive {
            self.compact();
        }
        Ok(())
    }

    fn compact(&mut self) {
        let items = std::mem::take(&mut self.items);
        self.items = items
            .into_iter()
            .zip(&self.alive)
            .filter(|(_, a)| **a)
            .map(|(c, _)| c)
            .collect();
        self.alive = vec![true; self.items.len()];
        self.by_hash.clear();
        for (i, c) in self.items.iter().enumerate() {
            self.by_hash
                .entry(hash_words(&c.planes))
                .or_default()
                .push(i as u32);
        }
    }

    fn into_items(self) -> Vec<Cand> {
        self.items
            .into_iter()
            .zip(self.alive)
            .filter(|(_, a)| *a)
            .map(|(c, _)| c)
            .collect()
    }
}

fn hash_words(w: &[u64]) -> u64 {
    let mut h: u64 = 0x9E37_79B9_7F4A_7C15;
    for &x in w {
        h = (h ^ x).wrapping_mul(0x0100_0000_01B3).rotate_left(23);
    }
    h
}

// ---------------------------------------------------------------------------
// Set search
// ---------------------------------------------------------------------------

/// Work limit of the set search, plus context for its error message.
struct SearchLimits {
    max_work: u64,
    generated: u64,
    set_size: usize,
}

/// Best way to add up to `slots` candidates of `pool` (sorted by stand-alone
/// score descending) to the set whose union profile is `state`. Returns the
/// best score, the chosen pool indices (empty when nothing improves on
/// `state`) and the number of evaluations spent.
#[allow(clippy::too_many_arguments)]
fn search_sets(
    obj: Objective,
    wt: &Weights,
    words: usize,
    pool: &[Cand],
    state: Vec<u64>,
    slots: usize,
    limits: &SearchLimits,
    progress: &dyn Progress,
) -> Result<(i64, Vec<usize>, u64), String> {
    let base = obj.union_score(&state, words, wt);
    let mut search = SetSearch {
        obj,
        wt,
        words,
        pool,
        score: base,
        best: base,
        best_set: Vec::new(),
        chosen: Vec::new(),
        saved: vec![vec![0u64; state.len()]; slots],
        state,
        evals: 0,
        max_work: limits.max_work,
        next_report: REPORT_EVERY,
        progress,
        generated: limits.generated,
        set_size: limits.set_size,
    };
    if slots > 0 && !pool.is_empty() {
        let root: Vec<(u32, i64)> = pool
            .iter()
            .enumerate()
            .map(|(i, c)| (i as u32, c.score))
            .collect();
        search.dfs(&root, slots, 0)?;
    }
    Ok((search.best, search.best_set, search.evals))
}

/// Depth-first branch-and-bound over sets of pool candidates.
///
/// A node holds the union profile of the chosen candidates and a list of
/// remaining candidates with upper bounds on their marginal gain. Bounds
/// computed at a node stay valid below it (gains only shrink as the set
/// grows), so a subtree is skipped when the current score plus the best
/// `k` remaining bounds — or the score of adding *all* remaining candidates
/// — cannot beat the incumbent. Children are visited in decreasing bound
/// order, so the first leaf reached is the greedy solution.
struct SetSearch<'a> {
    obj: Objective,
    wt: &'a Weights,
    words: usize,
    pool: &'a [Cand],
    /// Union profile of the chosen candidates (plus injected oligos).
    state: Vec<u64>,
    /// Per-depth copies of `state` for backtracking.
    saved: Vec<Vec<u64>>,
    score: i64,
    best: i64,
    best_set: Vec<usize>,
    chosen: Vec<usize>,
    evals: u64,
    max_work: u64,
    next_report: u64,
    progress: &'a dyn Progress,
    generated: u64,
    set_size: usize,
}

impl SetSearch<'_> {
    /// Marginal gain of `planes` on the current set, and an upper bound on
    /// its gain on any superset (`Exact` gains can turn negative later; only
    /// newly hit, previously uncovered sequences can add).
    fn delta(&self, planes: &[u64]) -> (i64, i64) {
        let (w, st) = (self.words, &self.state);
        match self.obj.op {
            MismatchOp::LowerOrEqual => {
                let mut d = 0i64;
                for j in 0..self.obj.planes {
                    let o = j * w;
                    let s: u64 = (0..w)
                        .map(|i| self.wt.word(i, planes[o + i] & !st[o + i]))
                        .sum();
                    d += self.obj.coef(j) * s as i64;
                }
                (d, d)
            }
            MismatchOp::Exact => {
                let (mut gain, mut loss) = (0u64, 0u64);
                for i in 0..w {
                    let (sx, slo) = (st[i], st[w + i]);
                    gain += self.wt.word(i, planes[i] & !sx & !slo);
                    loss += self.wt.word(i, planes[w + i] & sx & !slo);
                }
                (gain as i64 - loss as i64, gain as i64)
            }
        }
    }

    /// `suf[t]` = upper bound on the score after adding any subset of
    /// `list[t..]`: the score of adding all of them (`LowerOrEqual`), or the
    /// count plus every uncovered sequence some of them hit exactly (`Exact`).
    fn suffix_bounds(&self, list: &[(u32, i64, i64)]) -> Vec<i64> {
        let w = self.words;
        let mut acc = self.state.clone();
        let mut s = self.score;
        let mut suf = vec![0i64; list.len()];
        for t in (0..list.len()).rev() {
            let o = &self.pool[list[t].0 as usize].planes;
            match self.obj.op {
                MismatchOp::LowerOrEqual => {
                    for j in 0..self.obj.planes {
                        let b = j * w;
                        let mut g = 0u64;
                        for i in 0..w {
                            g += self.wt.word(i, o[b + i] & !acc[b + i]);
                            acc[b + i] |= o[b + i];
                        }
                        s += self.obj.coef(j) * g as i64;
                    }
                }
                MismatchOp::Exact => {
                    for i in 0..w {
                        s += self.wt.word(i, o[i] & !acc[i] & !self.state[w + i]) as i64;
                        acc[i] |= o[i];
                    }
                }
            }
            suf[t] = s;
        }
        suf
    }

    fn tick(&mut self, n: u64) -> Result<(), String> {
        self.evals += n;
        if self.max_work > 0 && self.evals > self.max_work {
            return Err(format!(
                "the set search exceeded the work limit of {} evaluations (--max-work; \
                 0 = no limit) before the best set could be proven. Pool: {} candidate \
                 oligos after reduction ({} generated), set size {}. Reduce --n-oligos, \
                 --mismatches or --ambiguities, narrow the search with --fixed, or raise \
                 --max-work.",
                self.max_work,
                self.pool.len(),
                self.generated,
                self.set_size
            ));
        }
        if self.evals >= self.next_report {
            self.next_report = self.evals + REPORT_EVERY;
            if self.progress.cancelled() {
                return Err("search cancelled".to_string());
            }
            self.progress.report(
                &format!(
                    "Optimize by mismatch: searching sets over {} candidates ({} evaluations)",
                    self.pool.len(),
                    self.evals
                ),
                70.0,
            );
        }
        Ok(())
    }

    fn record(&mut self, extra: Option<usize>) {
        self.best_set.clear();
        self.best_set.extend_from_slice(&self.chosen);
        self.best_set.extend(extra);
    }

    /// Explore every way of adding up to `k` candidates from `list` (sorted
    /// by bound descending, then pool index), each with an upper bound on its
    /// gain on the current set.
    fn dfs(&mut self, list: &[(u32, i64)], k: usize, depth: usize) -> Result<(), String> {
        if k == 1 {
            for &(i, bound) in list {
                if self.score + bound <= self.best {
                    break;
                }
                let (d, _) = self.delta(&self.pool[i as usize].planes);
                self.tick(1)?;
                if self.score + d > self.best {
                    self.best = self.score + d;
                    self.record(Some(i as usize));
                }
            }
            return Ok(());
        }

        let fresh: Vec<(u32, i64, i64)> = {
            let this = &*self;
            let eval = |&(i, _): &(u32, i64)| {
                let (d, b) = this.delta(&this.pool[i as usize].planes);
                (i, d, b)
            };
            if list.len() >= PAR_MIN_ITEMS {
                list.par_iter().map(eval).collect()
            } else {
                list.iter().map(eval).collect()
            }
        };
        self.tick(list.len() as u64)?;
        let mut fresh: Vec<(u32, i64, i64)> = fresh.into_iter().filter(|f| f.2 > 0).collect();
        fresh.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)));
        let len = fresh.len();
        if len == 0 {
            return Ok(());
        }
        let mut psum = vec![0i64; len + 1];
        for t in 0..len {
            psum[t + 1] = psum[t] + fresh[t].2;
        }
        let suf = self.suffix_bounds(&fresh);
        self.tick(len as u64)?;
        let pairs: Vec<(u32, i64)> = fresh.iter().map(|&(i, _, b)| (i, b)).collect();

        for t in 0..len {
            let top = psum[(t + k).min(len)] - psum[t];
            if (self.score + top).min(suf[t]) <= self.best {
                break;
            }
            let (i, d, _) = fresh[t];
            self.saved[depth].copy_from_slice(&self.state);
            let saved_score = self.score;
            for (s, o) in self.state.iter_mut().zip(self.pool[i as usize].planes.iter()) {
                *s |= o;
            }
            self.score += d;
            self.chosen.push(i as usize);
            if self.score > self.best {
                self.best = self.score;
                self.record(None);
            }
            let rest = &pairs[t + 1..];
            if !rest.is_empty() {
                let child_top = psum[(t + k).min(len)] - psum[t + 1];
                if (self.score + child_top).min(suf[t]) > self.best {
                    self.dfs(rest, k - 1, depth + 1)?;
                }
            }
            self.chosen.pop();
            self.state.copy_from_slice(&self.saved[depth]);
            self.score = saved_score;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Memory
// ---------------------------------------------------------------------------

fn try_zeroed(len: usize, what: &str) -> Result<Vec<u64>, String> {
    let mut v: Vec<u64> = Vec::new();
    v.try_reserve_exact(len)
        .map_err(|_| out_of_memory(what, len.saturating_mul(8)))?;
    v.resize(len, 0);
    Ok(v)
}

fn out_of_memory(what: &str, bytes: usize) -> String {
    let size = if bytes > 0 {
        format!(" (failed to allocate {:.1} MB)", bytes as f64 / 1e6)
    } else {
        String::new()
    };
    format!(
        "out of memory while {what}{size}. The search space is too large for this \
         machine: reduce --mismatches, --ambiguities or --n-oligos, or run --fixed on \
         a narrower slice."
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::types::{MismatchSettings, NoProgress, Orientation, SearchMode};

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    fn settings(fixed: bool, op: MismatchOp, n: usize, x: usize, y: usize) -> SearchSettings {
        SearchSettings {
            tm_threshold: 0.0,
            oligo_concentration_um: 0.2,
            na_concentration_mm: 50.0,
            mg_concentration_mm: 3.0,
            dntp_concentration_mm: 0.8,
            mode: SearchMode::OptimizeByMismatch,
            target_coverage_pct: 100.0,
            max_ambiguities: 0,
            exclude_n: false,
            only_twofold: false,
            orientation: Orientation::Forward,
            three_prime_match: 0,
            max_seeds: 0,
            fixed,
            mismatch: MismatchSettings {
                op,
                oligo_count: n,
                mismatches: x,
                ambiguities: y,
                max_candidates: 0,
                max_work: 0,
            },
        }
    }

    fn seqs(list: &[&str]) -> Vec<Vec<u8>> {
        list.iter().map(|s| s.as_bytes().to_vec()).collect()
    }

    /// Clades of a random reference with extra per-sequence noise.
    fn random_alignment(rng: &mut Lcg, nseq: usize, len: usize, rate: u64) -> Vec<Vec<u8>> {
        let reference: Vec<u8> = (0..len).map(|_| BASES[rng.below(4) as usize]).collect();
        let clades: Vec<Vec<u8>> = (0..3)
            .map(|_| {
                reference
                    .iter()
                    .map(|&b| if rng.below(100) < rate { BASES[rng.below(4) as usize] } else { b })
                    .collect()
            })
            .collect();
        (0..nseq)
            .map(|_| {
                let c = &clades[rng.below(3) as usize];
                c.iter()
                    .map(|&b| if rng.below(100) < rate / 2 { BASES[rng.below(4) as usize] } else { b })
                    .collect()
            })
            .collect()
    }

    /// Score of a reported set, recomputed from scratch.
    fn reported_score(sequences: &[Vec<u8>], s: &SearchSettings, res: &PrimerSearchResult) -> i64 {
        let x = s.mismatch.mismatches;
        let obj = Objective::new(s.mismatch.op, x, sequences.len() as u64);
        sequences
            .iter()
            .map(|q| {
                let best = res
                    .primers
                    .iter()
                    .filter_map(|p| {
                        let fwd = if s.is_reverse() {
                            reverse_complement(&p.sequence)
                        } else {
                            p.sequence.clone()
                        };
                        level_of(&fwd, &q[p.align_start..p.align_end], x, s.three_prime_match, s.is_reverse())
                    })
                    .min();
                obj.unit(best)
            })
            .sum()
    }

    /// The literal method: deconstruct every window's consensus into *every*
    /// variant (full product), filter, then try every set of up to `n`.
    /// Identical profiles are merged first, which cannot change the optimum.
    fn brute_force(sequences: &[Vec<u8>], s: &SearchSettings, excluded: &[Vec<u8>]) -> i64 {
        let ms = &s.mismatch;
        let x = ms.mismatches;
        let rev = s.is_reverse();
        let tm = s.tm_params();
        let obj = Objective::new(ms.op, x, sequences.len() as u64);
        let exclude = ExcludeSet::new(excluded);
        let len = sequences[0].len();
        let windows: Vec<(usize, usize)> = if s.fixed {
            vec![(0, len)]
        } else {
            let mut seen = HashSet::new();
            let mut out = Vec::new();
            for q in sequences {
                for st in 0..len {
                    if let Some(l) = determine_oligo_length(q, st, s.tm_threshold, &tm)
                        && st + l <= len
                        && seen.insert((st, st + l))
                    {
                        out.push((st, st + l));
                    }
                }
            }
            out
        };
        let mut profiles: HashSet<Vec<Option<usize>>> = HashSet::new();
        for &(st, en) in &windows {
            let wl = en - st;
            let obs: Vec<u8> = (0..wl)
                .map(|o| sequences.iter().fold(0u8, |m, q| m | base_mask(q[st + o])))
                .collect();
            let n3 = s.three_prime_match.min(wl);
            let in3 = |o: usize| if rev { o < n3 } else { o >= wl - n3 };
            let var: Vec<usize> = (0..wl).filter(|&o| obs[o].count_ones() >= 2).collect();
            let elig: Vec<usize> = var.iter().copied().filter(|&o| !in3(o)).collect();
            let y = ms.ambiguities.min(elig.len());
            for combo in combinations(elig.len(), y) {
                let amb: Vec<usize> = combo.iter().map(|&e| elig[e]).collect();
                // Per variable position: the options (code masks or bases).
                let options: Vec<(usize, Vec<u8>)> = var
                    .iter()
                    .map(|&o| {
                        if amb.contains(&o) {
                            (o, allowed_codes(obs[o], s.exclude_n, s.only_twofold))
                        } else {
                            (o, (0..4).filter(|b| obs[o] & (1 << b) != 0).map(|b| 1u8 << b).collect())
                        }
                    })
                    .collect();
                let mut idx = vec![0usize; options.len()];
                loop {
                    let mut oligo: Vec<u8> = sequences[0][st..en].to_vec();
                    for (k, (o, opts)) in options.iter().enumerate() {
                        oligo[*o] = mask_to_iupac(opts[idx[k]]);
                    }
                    let ok = !exclude.excludes(&oligo, rev)
                        && (s.fixed || calculate_tm(&oligo, &tm) >= s.tm_threshold);
                    if ok {
                        let prof: Vec<Option<usize>> = sequences
                            .iter()
                            .map(|q| level_of(&oligo, &q[st..en], x, s.three_prime_match, rev))
                            .collect();
                        profiles.insert(prof);
                    }
                    let mut k = options.len();
                    let mut done = true;
                    while k > 0 {
                        k -= 1;
                        idx[k] += 1;
                        if idx[k] < options[k].1.len() {
                            done = false;
                            break;
                        }
                        idx[k] = 0;
                    }
                    if done {
                        break;
                    }
                }
            }
        }
        let profiles: Vec<Vec<Option<usize>>> = profiles.into_iter().collect();
        let mut best = 0i64;
        for size in 1..=ms.oligo_count.min(profiles.len()) {
            for set in combinations(profiles.len(), size) {
                let score: i64 = (0..sequences.len())
                    .map(|q| obj.unit(set.iter().filter_map(|&p| profiles[p][q]).min()))
                    .sum();
                best = best.max(score);
            }
        }
        best
    }

    /// Randomised comparison against the literal brute force: fixed slices
    /// and searched windows, both criteria, IUPAC restrictions, 3' rule and
    /// both orientations.
    #[test]
    fn matches_brute_force() {
        let mut rng = Lcg(12345);
        for case in 0..160 {
            let fixed = case % 2 == 0;
            let op = if case % 4 < 2 { MismatchOp::LowerOrEqual } else { MismatchOp::Exact };
            let len = if fixed { 5 + rng.below(4) as usize } else { 9 + rng.below(3) as usize };
            let nseq = 4 + rng.below(6) as usize;
            let seqsv = random_alignment(&mut rng, nseq, len, 30);
            let n = 1 + rng.below(if fixed { 3 } else { 2 }) as usize;
            let x = rng.below(3) as usize;
            let y = rng.below(3) as usize;
            let mut s = settings(fixed, op, n, x, y);
            s.exclude_n = rng.below(3) == 0;
            s.only_twofold = rng.below(3) == 0;
            s.three_prime_match = [0, 0, 2][rng.below(3) as usize];
            if rng.below(2) == 0 {
                s.orientation = Orientation::Reverse;
            }
            if !fixed {
                // Short windows: a few bases reach this threshold.
                s.tm_threshold = -20.0 + rng.below(25) as f64;
            }
            let res = find_primers_by_mismatch(&seqsv, &s, &[], &[], &NoProgress)
                .unwrap_or_else(|e| panic!("case {case}: {e}"));
            assert!(res.primers.len() <= n, "case {case}: at most n oligos");
            let got = reported_score(&seqsv, &s, &res);
            let want = brute_force(&seqsv, &s, &[]);
            assert_eq!(got, want, "case {case}: fixed={fixed} op={op:?} n={n} x={x} y={y}");
            let rep = res.mismatch.as_ref().expect("report");
            assert_eq!(rep.counted, res.primers.iter().map(|p| p.coverage_count).sum::<usize>());
            assert_eq!(rep.level_counts.iter().sum::<usize>() + rep.not_covered, nseq);
        }
    }

    /// The same search on one thread and on many gives the same set.
    #[test]
    fn thread_count_does_not_change_result() {
        let mut rng = Lcg(99);
        let seqsv = random_alignment(&mut rng, 120, 30, 8);
        let s = {
            let mut s = settings(false, MismatchOp::LowerOrEqual, 3, 1, 1);
            s.tm_threshold = 30.0;
            s
        };
        let run = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            pool.install(|| find_primers_by_mismatch(&seqsv, &s, &[], &[], &NoProgress).unwrap())
        };
        let a = run(1);
        let b = run(8);
        let key = |r: &PrimerSearchResult| {
            r.primers
                .iter()
                .map(|p| (p.sequence.clone(), p.align_start, p.align_end, p.coverage_count))
                .collect::<Vec<_>>()
        };
        assert_eq!(key(&a), key(&b));
        assert_eq!(a.mismatch.unwrap().candidates_reduced, b.mismatch.unwrap().candidates_reduced);
    }

    /// Fixed slice, one mismatch allowed: a single exact oligo reaches every
    /// sequence that differs from it at one non-3' position.
    #[test]
    fn one_mismatch_covers_neighbours() {
        let s = settings(true, MismatchOp::LowerOrEqual, 1, 1, 0);
        let input = seqs(&["ACGTACGTAC", "ACGTACGTAC", "ACCTACGTAC", "ACGTACCTAC", "TTTTTTTTTT"]);
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        assert_eq!(res.primers.len(), 1);
        assert_eq!(res.primers[0].sequence, b"ACGTACGTAC".to_vec());
        assert_eq!(res.primers[0].coverage_count, 4);
        let rep = res.mismatch.unwrap();
        assert_eq!(rep.level_counts, vec![2, 2]);
        assert_eq!(rep.not_covered, 1);
    }

    /// Equal coverage is broken by fewer total mismatches, so a second oligo
    /// that only turns 1-mismatch hits into perfect ones is still used, and
    /// each sequence is credited to its best-matching oligo.
    #[test]
    fn credit_goes_to_best_match() {
        let s = settings(true, MismatchOp::LowerOrEqual, 2, 1, 0);
        let input = seqs(&["AAAAAAAA", "AAAAAAAA", "AAAAAAAA", "AAAACAAA", "AAAACAAA"]);
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        let got: Vec<(Vec<u8>, usize)> = res
            .primers
            .iter()
            .map(|p| (p.sequence.clone(), p.coverage_count))
            .collect();
        assert_eq!(got, vec![(b"AAAAAAAA".to_vec(), 3), (b"AAAACAAA".to_vec(), 2)]);
        assert_eq!(res.mismatch.unwrap().level_counts, vec![5, 0]);
    }

    /// With the 3' rule a mismatch in the protected bases is not covered,
    /// even when the total mismatch count is within the limit.
    #[test]
    fn three_prime_mismatch_is_not_covered() {
        let mut s = settings(true, MismatchOp::LowerOrEqual, 1, 1, 0);
        s.three_prime_match = 3;
        let input = seqs(&["ACGTACGTAC", "ACGTACGTAC", "ACGTACGTAG"]);
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        assert_eq!(res.primers[0].sequence, b"ACGTACGTAC".to_vec());
        assert_eq!(res.primers[0].coverage_count, 2);
        // Reverse orientation protects the alignment start instead.
        s.orientation = Orientation::Reverse;
        let input = seqs(&["ACGTACGTAC", "ACGTACGTAC", "TCGTACGTAC"]);
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        assert_eq!(res.primers[0].coverage_count, 2);
    }

    /// `exact` counts only sequences whose best match has exactly x
    /// mismatches: perfectly matched sequences do not count.
    #[test]
    fn exact_mode_ignores_better_matches() {
        let s = settings(true, MismatchOp::Exact, 1, 1, 0);
        // Three copies of A and one variant, one mismatch apart.
        let input = seqs(&["AAAAAAAA", "AAAAAAAA", "AAAAAAAA", "AAAACAAA"]);
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        // The perfect match for the copies would count only the variant; the
        // variant's own sequence hits the three copies with one mismatch.
        assert_eq!(res.primers[0].sequence, b"AAAACAAA".to_vec());
        let rep = res.mismatch.unwrap();
        assert_eq!(rep.counted, 3);
        assert_eq!(rep.level_counts, vec![1, 3]);
    }

    /// Ambiguity codes: with y = 1 the consensus code absorbs the variable
    /// position, and the 3' region never carries a code.
    #[test]
    fn ambiguity_code_from_consensus() {
        let mut s = settings(true, MismatchOp::LowerOrEqual, 1, 0, 1);
        let input = seqs(&["ACGTAAGTAC", "ACGTAGGTAC", "ACGTACGTAC"]);
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        assert_eq!(res.primers[0].sequence, b"ACGTAVGTAC".to_vec());
        assert_eq!(res.primers[0].coverage_count, 3);
        assert_eq!(res.primers[0].ambiguity_count, 1);
        s.exclude_n = true;
        s.only_twofold = true;
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        assert_eq!(res.primers[0].coverage_count, 2, "a 2-fold code covers two of three");
        s.only_twofold = false;
        s.three_prime_match = 5;
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        assert_eq!(res.primers[0].ambiguity_count, 0, "no codes in the 3' region");
    }

    /// Injected oligos are fixed members and count toward n; exclusion also
    /// applies in fixed-slice mode.
    #[test]
    fn injection_and_exclusion() {
        let s = settings(true, MismatchOp::LowerOrEqual, 2, 0, 0);
        let input = seqs(&["ACGTACGT", "ACGTACGT", "ACGTACGT", "TTGTACGA", "TTGTACGA", "GGGTACGC"]);
        let res = find_primers_by_mismatch(&input, &s, &[b"GGGTACGC".to_vec()], &[], &NoProgress).unwrap();
        assert_eq!(res.primers.len(), 2);
        assert!(res.primers[0].injected);
        assert_eq!(res.primers[0].coverage_count, 1);
        assert_eq!(res.primers[1].sequence, b"ACGTACGT".to_vec());
        let three = vec![b"ACGT".to_vec(); 3];
        assert!(find_primers_by_mismatch(&input, &s, &three, &[], &NoProgress).is_err());

        let res = find_primers_by_mismatch(&input, &s, &[], &[b"ACGTACGT".to_vec()], &NoProgress).unwrap();
        assert!(res.primers.iter().all(|p| p.sequence != b"ACGTACGT".to_vec()));
        assert_eq!(res.primers[0].sequence, b"TTGTACGA".to_vec());
    }

    /// Search mode gates each candidate on its own Tm.
    #[test]
    fn tm_gate_applies_to_candidates() {
        let mut s = settings(false, MismatchOp::LowerOrEqual, 1, 0, 0);
        s.tm_threshold = 40.0;
        let input = seqs(&["GCGCGCGCGCGCGCGCAAAA", "GCGCGCGCGCGCGCGCAAAA"]);
        let res = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap();
        assert!(res.primers.iter().all(|p| p.tm >= 40.0));
    }

    /// Both work limits abort with an explanatory error.
    #[test]
    fn limits_abort_with_error() {
        let mut rng = Lcg(7);
        let input = random_alignment(&mut rng, 60, 14, 25);
        let mut s = settings(true, MismatchOp::LowerOrEqual, 3, 1, 1);
        s.mismatch.max_candidates = 10;
        let e = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap_err();
        assert!(e.contains("--max-candidates"), "{e}");
        s.mismatch.max_candidates = 0;
        s.mismatch.max_work = 1;
        let e = find_primers_by_mismatch(&input, &s, &[], &[], &NoProgress).unwrap_err();
        assert!(e.contains("--max-work"), "{e}");
    }

    /// Random profile over `elems` weighted elements: a few clustered runs
    /// of covered elements with random levels, so candidates overlap the way
    /// neighbouring windows do.
    fn random_planes(rng: &mut Lcg, obj: Objective, elems: usize) -> Vec<Option<usize>> {
        let mut levels = vec![None; elems];
        for _ in 0..1 + rng.below(3) {
            let start = rng.below(elems as u64) as usize;
            let run = 1 + rng.below(elems as u64 / 3) as usize;
            for lv in &mut levels[start..(start + run).min(elems)] {
                *lv = Some(rng.below(obj.x as u64 + 1) as usize);
            }
        }
        levels
    }

    fn planes_of(levels: &[Option<usize>], obj: Objective, words: usize) -> Vec<u64> {
        let mut planes = vec![0u64; obj.planes * words];
        for (e, lv) in levels.iter().enumerate() {
            let Some(l) = *lv else { continue };
            let m = 1u64 << (e % 64);
            match obj.op {
                MismatchOp::LowerOrEqual => {
                    for j in l..=obj.x {
                        planes[j * words + e / 64] |= m;
                    }
                }
                MismatchOp::Exact => planes[(l != obj.x) as usize * words + e / 64] |= m,
            }
        }
        planes
    }

    fn dummy_cand(i: u32, planes: Vec<u64>, obj: Objective, words: usize, wt: &Weights) -> Cand {
        let (score, sec) = obj.standalone(&planes, words, wt);
        Cand {
            key: (0, 0, i),
            score,
            sec,
            planes: planes.into_boxed_slice(),
            oligo: Box::new([]),
            start: 0,
            end: 0,
        }
    }

    /// The branch-and-bound finds the same optimum as trying every set, on
    /// pools far larger than the end-to-end brute force can reach, with and
    /// without a pre-filled (injected) state.
    #[test]
    fn set_search_matches_exhaustive() {
        let mut rng = Lcg(2024);
        for case in 0..60 {
            let op = if case % 2 == 0 { MismatchOp::LowerOrEqual } else { MismatchOp::Exact };
            let x = rng.below(3) as usize;
            let elems = 70 + rng.below(130) as usize;
            let words = words_for(elems);
            let weights: Vec<u64> = (0..elems).map(|_| 1 + rng.below(4)).collect();
            let wt = Weights::new(&weights);
            let obj = Objective::new(op, x, weights.iter().sum());
            let size = 20 + rng.below(50) as usize;
            let mut pool: Vec<Cand> = (0..size)
                .map(|i| {
                    let lv = random_planes(&mut rng, obj, elems);
                    dummy_cand(i as u32, planes_of(&lv, obj, words), obj, words, &wt)
                })
                .collect();
            pool.sort_by(|a, b| b.score.cmp(&a.score).then(a.sec.cmp(&b.sec)).then(a.key.cmp(&b.key)));
            let state = if case % 3 == 0 {
                planes_of(&random_planes(&mut rng, obj, elems), obj, words)
            } else {
                vec![0u64; obj.planes * words]
            };
            let slots = 1 + rng.below(3) as usize;
            let limits = SearchLimits { max_work: 0, generated: 0, set_size: slots };
            let (best, set, _) =
                search_sets(obj, &wt, words, &pool, state.clone(), slots, &limits, &NoProgress).unwrap();

            let score_of = |members: &[usize]| {
                let mut u = state.clone();
                for &m in members {
                    for (a, b) in u.iter_mut().zip(pool[m].planes.iter()) {
                        *a |= b;
                    }
                }
                obj.union_score(&u, words, &wt)
            };
            let mut want = score_of(&[]);
            for k in 1..=slots.min(pool.len()) {
                for c in combinations(pool.len(), k) {
                    want = want.max(score_of(&c));
                }
            }
            assert_eq!(best, want, "case {case}: op={op:?} x={x} slots={slots} pool={size}");
            assert_eq!(score_of(&set), best, "case {case}: reported set reaches the score");
            assert!(set.len() <= slots);
        }
    }

    /// The skyline keeps exactly the valid candidates no other valid
    /// candidate dominates (identical profiles: the first one).
    #[test]
    fn skyline_keeps_exactly_the_non_dominated() {
        let mut rng = Lcg(77);
        for case in 0..40 {
            let op = if case % 2 == 0 { MismatchOp::LowerOrEqual } else { MismatchOp::Exact };
            let x = rng.below(3) as usize;
            let elems = 40 + rng.below(60) as usize;
            let words = words_for(elems);
            let weights: Vec<u64> = (0..elems).map(|_| 1 + rng.below(3)).collect();
            let wt = Weights::new(&weights);
            let obj = Objective::new(op, x, weights.iter().sum());
            let mut cands: Vec<(Vec<u64>, bool)> = Vec::new();
            for _ in 0..150 {
                let planes = if !cands.is_empty() && rng.below(4) == 0 {
                    // A duplicate or a weakened copy of an earlier profile.
                    let mut p = cands[rng.below(cands.len() as u64) as usize].0.clone();
                    if rng.below(2) == 0 {
                        let e = rng.below(elems as u64) as usize;
                        for j in 0..obj.planes {
                            p[j * words + e / 64] &= !(1u64 << (e % 64));
                        }
                    }
                    p
                } else {
                    planes_of(&random_planes(&mut rng, obj, elems), obj, words)
                };
                if planes[..words].iter().all(|&w| w == 0) && op == MismatchOp::Exact {
                    continue;
                }
                cands.push((planes, rng.below(5) != 0));
            }
            let mut sky = Skyline::new(obj, words);
            for (i, (p, valid)) in cands.iter().enumerate() {
                // An invalid candidate is never committed (and so displaces
                // nothing), exactly as when generation rejects it.
                if *valid {
                    sky.insert(dummy_cand(i as u32, p.clone(), obj, words, &wt)).unwrap();
                }
            }
            let mut got: Vec<u32> = sky.into_items().iter().map(|c| c.key.2).collect();
            got.sort();

            let beats = |a: &[u64], b: &[u64]| obj.dominates(a, b, words);
            let want: Vec<u32> = (0..cands.len())
                .filter(|&i| cands[i].1)
                .filter(|&i| {
                    !(0..cands.len()).any(|j| {
                        j != i
                            && cands[j].1
                            && beats(&cands[j].0, &cands[i].0)
                            && (cands[j].0 != cands[i].0 || j < i)
                    })
                })
                .map(|i| i as u32)
                .collect();
            assert_eq!(got, want, "case {case}: op={op:?} x={x}");
        }
    }
}
