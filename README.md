# primersearch

Find an efficient set of PCR primers covering an aligned DNA sequence set.

A Rust CLI port of an earlier Tkinter-based Python tool (see
`reference_program/`). Same analysis logic, multi-threaded, no GUI.

## What it does

Given a FASTA file of aligned DNA sequences (same length, A/C/G/T only —
gaps and ambiguous bases are filtered out as a preprocessing step), the
program returns a small set of primers that together match every input
sequence. Each round, it picks the primer that covers the largest number
of *remaining* sequences, removes those, and repeats until coverage is
100%.

Two greedy search modes:

- **`no-ambiguities`** — exact-match primers only.
- **`incremental`** — allow IUPAC ambiguity codes (R, Y, …, up to a
  user-set budget) to absorb closely-related variants into one primer.
  Falls back to exact-match if no ambiguity expansion qualifies.

And one exhaustive mode:

- **`optimize-by-mismatch`** — the best set of (at most) `n` oligos, each
  with a given number of ambiguity codes, when sequences may be bound with
  mismatches. See "Optimize-by-mismatch mode" below.

There is also a **fixed-slice** mode (`--fixed`) for when you already know
*where* the primer should sit. See "Fixed-slice mode" below.

Constraints applied to each candidate primer: minimum Tm (nearest-neighbor
thermodynamic calculation — see "Tm calculation" below), maximum ambiguity
count, optional 3' conserved tail, optional ban on N / non-2-fold codes,
forward or reverse orientation.

## Build

Requires Rust 1.85+ (edition 2024).

```
cargo build --release
```

The binary lands at `target/release/primersearch[.exe]`. For day-to-day
use, copy it to a directory of your choice — `settings.ini` will live
next to the binary.

## Usage

Simplest:

```
primersearch input.fasta
```

Writes `output.txt` in the current directory.

Common invocations:

```
primersearch input.fasta -o results.txt
primersearch input.fasta -o results.txt --rev
primersearch input.fasta --tm 62 --na 200 --mode incremental --target 60 --max-amb 3 --exclude-n --three-prime 5
primersearch slice.fasta --fixed --mode incremental --max-amb 3   # generate variants for a fixed slice
primersearch input.fasta --inject ACGTTGCACGTACGTACGT   # force one or more obligatory oligos
primersearch input.fasta --exclude CTAAATCYCGTG          # forbid candidates matching a 3' signature
primersearch input.fasta --mode optimize-by-mismatch --n-oligos 3 --mismatches 1 --ambiguities 2
primersearch --mkini                 # write defaults to settings.ini
primersearch input.fasta -j 8        # 8 worker threads
primersearch input.fasta --silent    # no progress / info output
```

### Fixed-slice mode

`--fixed` skips the search entirely. Instead of hunting for the alignment
positions that let one primer cover the most sequences, you hand the tool a
FASTA that has already been trimmed to exactly the region you want the
primer at — the whole alignment is treated as a single slice — and it just
generates the oligo variant(s) needed to cover every input sequence there.

It still runs the same per-round greedy coverage loop, but over that one
fixed slice: each round emits the highest-coverage variant over the
sequences not yet covered and removes them, repeating until coverage is
100%. The configured `--mode` chooses how variants are formed:

- `--fixed --mode no-ambiguities` — one exact primer per distinct sequence
  in the slice, ordered by coverage.
- `--fixed --mode incremental` — IUPAC consensus primers that absorb
  related variants up to `--max-amb` (and respect `--exclude-n`,
  `--only-twofold`, `--three-prime`, `--target`).

Orientation (`--rev`/`--fwd`) applies as usual. The Tm threshold (`--tm`)
is **not** enforced in fixed mode — you chose the region, so every variant
required for full coverage is emitted regardless of its Tm; each variant's
Tm is still computed and reported so you can judge the choice. The input is
still quality-filtered (same length, A/C/G/T only), and every primer's
reported position spans the full slice.

### Injecting obligatory oligos

`--inject` lets you hand the tool one or more oligos that it **must** include
in the primer set — for example primers you have already validated in the lab
and want to keep, while letting the search fill in coverage for whatever they
miss.

Injected oligos are processed *before* the search runs (before any new
positions are searched or new variants generated):

1. Each oligo is positioned at the alignment offset where it covers the most
   input sequences (its intrinsic best-fit position).
2. It is emitted as a primer (listed first, marked `[injected]` in the output)
   whose reported coverage is the number of sequences it matches.
3. The sequences it covers are removed from the pool, so the search only has
   to cover what the injected oligos left behind.

Supply several oligos by repeating the flag or with a comma-separated list:

```
primersearch input.fasta --inject ACGTTGCA... --inject GGCATTAC...
primersearch input.fasta --inject ACGTTGCA...,GGCATTAC...
```

Details:

- Oligos may contain IUPAC ambiguity codes (e.g. an existing degenerate
  primer); each one then covers every A/C/G/T variant the codes admit.
- Provide oligos in the **same orientation as the run**. For a `--rev` run
  that means the reverse-complement form you would actually order; it is
  matched in alignment coordinates and echoed back in the form you supplied.
- Injected oligos are **obligatory**: the Tm threshold and the
  ambiguity / 3' / IUPAC restrictions are *not* enforced on them. They are
  always emitted, and their Tm and ambiguity count are still computed and
  reported so you can judge them. (The regular search applied to the leftover
  sequences still respects all of those constraints.)
- An oligo longer than the alignment, empty, or containing a non-IUPAC
  character is rejected with an error before the run starts.
- `--inject` works in both the regular search and `--fixed` mode.

### Excluding primers by 3' signature

`--exclude` lets you hand the tool one or more primer sequences that the search
must **not** reproduce — for example primers already used elsewhere in a
multiplex, whose 3' end you want to keep clear of new candidates. Because primer
interactions are dominated by the 3' end, the match is **anchored at the 3'
end**:

1. The candidate and the excluded primer are aligned at their 3' ends (right
   edge).
2. Over the shorter of the two lengths, every position must share at least one
   base under IUPAC codes (sets must *intersect*). A single mismatch anywhere in
   that overlap — even far from the 3' end — means it is a different signature
   and the candidate is kept.
3. Whenever a candidate would be the best pick for a region but matches an
   excluded signature, it is dropped and the search takes the best
   **non-excluded** candidate for that region instead.

With `Excluded: CTAAATCYCGTG`:

```
AAACTAAATCCCGTG   -> excluded (last 12 bases match under Y={C,T})
AAACTAAATCTCRTG   -> excluded (R={A,G} still shares G; T shares with Y)
TAAATCYCGTG       -> excluded (shorter; its whole length matches the 3' tail)
TTTAAACTAAATCCCGT -> kept (3'-most base differs: ...GT vs ...TG)
AAACTAAATCCCGTGA  -> kept (extra 3' base shifts the anchor)
GATAATCYCGTG      -> kept (matches the last 10 bases, differs at the 5' end)
```

Supply several by repeating the flag or with a comma-separated list:

```
primersearch input.fasta --exclude CTAAATCYCGTG --exclude ACGT...
primersearch input.fasta --exclude CTAAATCYCGTG,ACGT...
```

Details:

- Excluded primers may contain IUPAC ambiguity codes, and may be longer than the
  alignment (a candidate is then matched against the excluded primer's 3' tail
  over the candidate's own length).
- Provide them in the **same orientation as the run**, exactly like `--inject`:
  for a `--rev` run that means the reverse-complement form you would order.
  Candidates are compared in their display form.
- `--inject`ed oligos are **exempt** — they are obligatory and user-forced, so
  they are emitted even if they match an excluded signature.
- Exclusion applies to the regular search only. It is **ignored in `--fixed`
  mode** (which must emit whatever covers the fixed slice); combining
  `--exclude` with `--fixed` prints a warning.
- Running without `--exclude` produces exactly the same results as before the
  feature existed — the check is skipped entirely when no signatures are given.
- In `optimize-by-mismatch` mode exclusion also applies with `--fixed` (see
  below).

### Optimize-by-mismatch mode

`--mode optimize-by-mismatch` replaces the greedy round loop with an exhaustive
search for the **best set of `n` oligos** (`--n-oligos`), where each oligo
carries exactly `y` IUPAC ambiguity codes (`--ambiguities`) and sequences may
be bound with mismatches (`--mismatches`). It never lists more than `n`
oligos, even if sequences remain uncovered.

```
primersearch input.fasta --mode optimize-by-mismatch --n-oligos 3 --mismatches 1 --ambiguities 2
primersearch slice.fasta --fixed --mode optimize-by-mismatch --n-oligos 2 --mismatches 1 --mismatch-mode exact
```

Every sequence is scored by its **best-matching oligo** in the set (fewest
mismatches). The coverage criterion is chosen with `--mismatch-mode`:

- `lower-or-equal` (default) — a sequence counts when its best match has at
  most `--mismatches` mismatches. Among sets with equal coverage, the one with
  the fewest total mismatches over the covered sequences wins.
- `exact` — a sequence counts only when its best match has *exactly*
  `--mismatches` mismatches; sequences matched better than that do not count.

Where the candidates come from:

- With `--fixed` the whole alignment is one slice. Its IUPAC consensus is
  deconstructed into variants: `y` of the variable positions keep their
  consensus code and every other variable position is resolved to one of the
  bases observed there (exactly `y` codes, or fewer when the slice has fewer
  variable positions). The Tm threshold is not enforced.
- Otherwise the windows are the same Tm-derived ranges the regular search
  uses (from every start position of every sequence, the shortest length
  reaching `--tm`), each deconstructed the same way, and every candidate's own
  Tm (median over its variants) must reach `--tm`. Candidates from all
  windows compete, so the oligos of a set may sit at different positions.

A variant more than `--mismatches` away from every input sequence can never
cover anything, so those are skipped. The result is the same as enumerating
every variant.

How the other options apply:

- `--three-prime N` — no ambiguity codes in the 3'-most `N` bases, **and** a
  sequence with a mismatch in those bases counts as not covered.
- `--exclude-n` / `--only-twofold` — where the consensus code at a position is
  forbidden, the widest allowed sub-codes are used instead.
- `--exclude` — candidates matching an excluded 3' signature are dropped, in
  this mode also with `--fixed` (it chooses among candidates rather than
  having to cover every sequence).
- `--inject` — injected oligos are fixed members of the set and **count
  toward `--n-oligos`**. Each is placed where it alone scores best, under the
  same mismatch rules.
- `--target`, `--max-amb` and `--max-seeds` are not used.

Output: the table lists the set, injected oligos first, then by decreasing
Count. Each counted sequence is credited to its best-matching oligo (ties to
the one listed first), so Total% adds up to the set's coverage. A mismatch
breakdown follows the table: how many sequences the set binds with 0, 1, …
mismatches, and how many it does not cover. If further oligos cannot improve
the result, fewer than `n` are listed and a note says so.

#### Search space and limits

The reported set is optimal for the criterion, which makes the method
exponential in the worst case. Candidates with an identical coverage profile,
or that another candidate matches at least as well on every sequence, are
dropped first (this cannot change the optimum), and sets are searched with
branch-and-bound. Two limits keep runs bounded; exceeding either aborts with
an error explaining how to shrink the problem, and no partial result is
reported:

- `--max-candidates N` (default 500,000,000) — candidate variants that would
  be enumerated, estimated before any work starts.
- `--max-work N` (default 2,000,000,000) — candidate evaluations in the set
  search.

`0` disables a limit. The largest allocations are checked, so running out of
memory normally ends with an error message rather than a crash. The biggest
levers are `--mismatches` and `--ambiguities` (number of candidates), then
`--n-oligos` (set search); `--fixed` on a narrow slice is far cheaper than
searching the whole alignment, which pools candidates from every window.

The mode can be selected in `settings.ini` (`search_mode =
"optimize_by_mismatch"`); its parameters are CLI-only (defaults: 1 oligo,
0 mismatches, 0 ambiguities, `lower-or-equal`).

### Settings precedence

`built-in defaults` < `settings.ini` < `CLI flags`. The settings file
lives next to the executable and uses TOML syntax (the `.ini` extension is
just a hint to the user). Run `primersearch --mkini` to write a fresh
defaults file, or pass `--config <path>` to point at an alternative.

### Selected flags

| flag | meaning |
|------|---------|
| `-o, --output FILE` | Output destination. Text mode defaults to `output.txt`; JSON mode defaults to stdout. `-o -` forces stdout in either mode |
| `--format {text,json}` | Output format (default `text`). `json` emits a structured document for programmatic consumers; see "Programmatic / subprocess use" |
| `--no-config` | Ignore the settings file entirely: built-in defaults overlaid by CLI flags only, no `settings.ini` read or created |
| `--rev` / `--fwd` | Search reverse / forward orientation |
| `--tm C` | Minimum Tm in °C |
| `--oligo UM` | Oligo (primer) concentration in µM |
| `--na MM` | Na⁺ concentration in mM |
| `--mg MM` | Mg²⁺ concentration in mM |
| `--dntp MM` | dNTP concentration in mM (one Mg²⁺ sequestered per dNTP) |
| `--mode {no-ambiguities,incremental,optimize-by-mismatch}` | Variant-generation mode |
| `--fixed` | Skip the search; treat the whole input as one slice and generate the variants needed to cover it (Tm threshold not enforced) |
| `--target PCT` | Coverage target % at which the ambiguity counter is allowed to increase (incremental) |
| `--max-amb N` | Maximum ambiguity codes per primer (incremental) |
| `--exclude-n` / `--only-twofold` | IUPAC restrictions (incremental) |
| `--three-prime N` | Number of 3' bases that must be perfectly conserved |
| `--inject OLIGO` | Obligatory oligo placed before the search (repeatable / comma-separated). See "Injecting obligatory oligos" |
| `--exclude OLIGO` | Primer whose 3' signature must not be reproduced by the search (repeatable / comma-separated; ignored in `--fixed` except in optimize-by-mismatch). See "Excluding primers by 3' signature" |
| `--max-seeds N` | Per-range seed cap in incremental mode (0 = no cap, default 50) |
| `--n-oligos N` | Oligos in the optimized set (optimize-by-mismatch; injected oligos count toward it) |
| `--mismatches X` | Mismatch count of the coverage criterion (optimize-by-mismatch) |
| `--ambiguities Y` | Ambiguity codes per oligo (optimize-by-mismatch) |
| `--mismatch-mode {lower-or-equal,exact}` | Count sequences matched with at most / exactly `X` mismatches (optimize-by-mismatch) |
| `--max-candidates N` / `--max-work N` | Search-space limits of optimize-by-mismatch (0 = no limit) |
| `-j, --threads N` | Worker threads (0 = all logical cores) |
| `-s, --silent` | Suppress progress / info |
| `--mkini` | Write a default `settings.ini` and exit |

`primersearch --help` for the full list.

## Output

The output file contains:

1. A **preprocessing report** — original / valid sequence counts,
   majority alignment length, and a breakdown of why sequences were
   removed (gaps, ambiguous bases, invalid characters, wrong length).
2. A **results block** — search settings echoed back, then a table of
   primers with: position in the alignment, sequence (with optional
   spacing every 3 bases), coverage count, per-primer and cumulative
   coverage %, and Tm of the displayed primer.

The same preprocessing report is also written to stdout during the run
so you can abort early if the input data looks bad. Progress messages
go to stderr (and are suppressed by `--silent`).

## Programmatic / subprocess use

primersearch is designed to be driven as a child process by another tool
(e.g. a GUI front-end) that passes inputs as CLI flags and captures the
result. For that use, three opt-in flags make a run self-describing,
deterministic, and easy to capture:

```
primersearch input.fasta --format json --no-config --silent [other flags...]
```

- `--format json` emits a single JSON document instead of the text table.
  By default it goes to **stdout** (use `-o FILE` to write it to a file
  instead); in JSON mode the human-readable preprocessing echo on stdout is
  suppressed so stdout carries *only* the JSON.
- `--no-config` makes the run depend solely on built-in defaults plus the
  flags you pass — no `settings.ini` is read or created, so the result does
  not depend on hidden file state or the binary's location/permissions.
- `--silent` suppresses the progress spinner and the `primersearch: …`
  informational lines on stderr. Errors are *not* suppressed: any failure
  still prints `error: …` to stderr and exits non-zero (exit code `1`;
  argument-parsing errors exit `2`).

The JSON document (`format_version: 1`) has three top-level objects:

- `preprocessing` — `original_count`, `valid_count`, `majority_length`, and
  a `removed` breakdown (`gaps`, `ambiguous`, `invalid`, `wrong_length`,
  `total`).
- `settings` — the fully resolved run settings (`mode`, `fixed`,
  `orientation`, Tm and concentrations, `tm_threshold_enforced` — `false`
  in fixed mode — etc.).
- `result` — `total_sequences`, `primer_count`, `injected_count`,
  `message`, and a `primers` array. Each primer has `index`, `sequence`
  (display orientation, no spacers), `coverage_count`, `coverage_pct`,
  `cumulative_pct`, `tm`, `ambiguity_count`, `injected`, the zero-based
  half-open range `align_start` / `align_end`, and a 1-based
  `position_label` (e.g. `"9-20"`) matching the text output.

In `optimize-by-mismatch` mode, `settings` also carries `mismatch_mode`
(`"lower_or_equal"` / `"exact"`), `n_oligos`, `mismatches`, `ambiguities`,
`max_candidates` and `max_work`, and `result` carries `mismatch_breakdown`:
`counted` / `counted_pct` (sequences meeting the criterion), `levels` (one
`{mismatches, count, pct}` entry per mismatch count `0..=mismatches`, each
sequence scored by its best-matching oligo), `not_covered` /
`not_covered_pct`, and the search statistics `windows`,
`candidates_generated`, `candidates_after_reduction` and `evaluations`. These
keys are absent in the other modes, whose documents are unchanged.

Bump-guard on `format_version` before parsing; it increments only on a
backward-incompatible shape change.

## How it works

```
parse + quality-filter FASTA
while remaining sequences exist:
    Phase 1 (single-thread): walk every (sequence, start_pos),
        compute the smallest length whose subsequence reaches the Tm
        threshold, and collect the unique (start, end) alignment ranges.
    Phase 2 (parallel via rayon): for each unique range, evaluate
        the highest-coverage primer that fits the constraints.
    Phase 3 (single-thread reduce): pick the highest-coverage candidate
        and remove its covered sequences.
output table
```

Per-range evaluation depends on mode:

- **No-ambiguities**: bucket the slices, return the most frequent.
- **Incremental**: for each ambiguity budget from 0 up to `max_amb`, try
  the top `max_seeds` unique slices as seeds. For each seed, run greedy
  expansion in three different orderings (first-appearance, count
  descending, count ascending) and keep the highest-coverage consensus.
  Stop early once the coverage target is met. Both the seed cap and the
  multi-ordering expansion are heuristics — see `claude_evaluation.md`
  for a detailed discussion of where they may fall short of the true
  per-round optimum.

Output is byte-identical across thread counts: phase 1 produces ranges in
a fixed order, `rayon::par_iter().collect()` preserves that order, and
the phase 3 reduce uses a deterministic tie-break.

Fixed-slice mode (`--fixed`) reuses the same round loop and the same
per-range evaluators, but skips phase 1 (range discovery) and phase 2's
parallelism: there is exactly one range — the whole slice `[0, len)` — so
each round evaluates just that, with the Tm threshold disabled as a gate.

Optimize-by-mismatch (`src/engine/mismatch.rs`) does not use the round loop:

```
windows: the fixed slice, or phase-1 ranges over all sequences
per window (parallel): enumerate consensus variants within reach of
    some slice, compute each one's coverage profile (mismatch level per
    sequence, as bitsets), keep only valid, non-dominated candidates
pool: merge the windows' candidates, again keeping non-dominated ones
branch-and-bound over sets of up to n pool candidates (parallel branches)
```

Results are identical across thread counts: the reduced pool does not depend
on processing order, and the set search, although it explores branches in
parallel, breaks ties by position in a fixed visiting order and so keeps the
same best set as a sequential search. Only the evaluation count can vary a
little between multi-threaded runs, since a branch started in parallel may
not yet know the best set found by an earlier one.

## Tm calculation

Tm is computed from a full nearest-neighbor thermodynamic model
(`src/engine/tm.rs`), not a salt-adjusted GC fraction:

- ΔH / ΔS parameters from SantaLucia (1998) unified table, with
  end-dependent initiation (G/C vs. A/T) for both 5' and 3' ends.
- Salt correction: `ΔS_salt = ΔS + 0.368·(N-1)·ln([Na⁺]_eq)` where
  `[Na⁺]_eq = [Na⁺] + 120·√(Mg_free)` and `Mg_free = max(Mg − dNTP, 0)`
  (one Mg²⁺ sequestered per dNTP — von Ahsen 2001 / Owczarzy 2008).
- Strand-concentration term: `Tm = ΔH·1000 / (ΔS_salt + R·ln(C_T/4)) − 273.15`
  with the user-supplied oligo concentration interpreted as a single
  strand of a non-self-complementary duplex (`C_T = 2·oligo_conc`).
- IUPAC ambiguity codes in a primer expand to all A/C/G/T variants;
  the reported Tm is the median across variants.

For the operating conditions used by this tool (oligo = 0.2 µM,
Na⁺ = 50 mM, Mg²⁺ = 3 mM, dNTP = 0.8 mM), the predicted Tms agree with
a commercial NN calculator (see `example_sets/example_oligo_tm.csv`) to
within 2 °C across the reference panel.

## Repository layout

```
src/
  main.rs          CLI entry, argument resolution, top-level pipeline
  cli.rs           clap-based argument definitions
  config.rs        settings.ini load / write (TOML syntax)
  output.rs        text rendering of the preprocessing + results blocks
  progress.rs      indicatif-based progress sink
  engine/          self-contained analysis logic — depends only on std
                   and rayon, no CLI / serde / I/O. Drop the directory
                   into another project to reuse:
    mod.rs           module re-exports + the public API surface
    types.rs         SearchSettings / PrimerCandidate / etc., Progress trait
    iupac.rs         bitmask-based IUPAC code utilities
    fasta.rs         FASTA parser + quality filter
    tm.rs            nearest-neighbor thermodynamic Tm
    search.rs        greedy round loop, per-range evaluators
    mismatch.rs      optimize-by-mismatch: candidate enumeration,
                     dominance reduction, branch-and-bound set search
example_sets/      reference inputs and Python-tool outputs for testing
reference_program/ original Python implementation (kept for reference)
program_instructions.md  initial spec / Q&A
claude_evaluation.md     algorithmic evaluation and improvement notes
```

## Example reference data

```
example_sets/example1_fasta.fasta   2741 sequences, length 123
example_sets/example1_result1.txt   Python output for comparison

example_sets/example2_fasta.fasta   608 sequences, length 80
example_sets/example2_result1.txt   Python output for comparison
```

To reproduce the example1 reference run:

```
primersearch example_sets/example1_fasta.fasta -o ex1.txt --tm 62 --na 200 --mode no-ambiguities --three-prime 5
```

For example2:

```
primersearch example_sets/example2_fasta.fasta -o ex2.txt --rev --tm 62 --na 200 --mode incremental --target 60 --max-amb 3 --exclude-n --three-prime 5
```
