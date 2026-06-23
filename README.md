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

Two search modes:

- **`no-ambiguities`** — exact-match primers only.
- **`incremental`** — allow IUPAC ambiguity codes (R, Y, …, up to a
  user-set budget) to absorb closely-related variants into one primer.
  Falls back to exact-match if no ambiguity expansion qualifies.

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

### Settings precedence

`built-in defaults` < `settings.ini` < `CLI flags`. The settings file
lives next to the executable and uses TOML syntax (the `.ini` extension is
just a hint to the user). Run `primersearch --mkini` to write a fresh
defaults file, or pass `--config <path>` to point at an alternative.

### Selected flags

| flag | meaning |
|------|---------|
| `-o, --output FILE` | Output file (default `output.txt`) |
| `--rev` / `--fwd` | Search reverse / forward orientation |
| `--tm C` | Minimum Tm in °C |
| `--oligo UM` | Oligo (primer) concentration in µM |
| `--na MM` | Na⁺ concentration in mM |
| `--mg MM` | Mg²⁺ concentration in mM |
| `--dntp MM` | dNTP concentration in mM (one Mg²⁺ sequestered per dNTP) |
| `--mode {no-ambiguities,incremental}` | Variant-generation mode |
| `--fixed` | Skip the search; treat the whole input as one slice and generate the variants needed to cover it (Tm threshold not enforced) |
| `--target PCT` | Coverage target % at which the ambiguity counter is allowed to increase (incremental) |
| `--max-amb N` | Maximum ambiguity codes per primer (incremental) |
| `--exclude-n` / `--only-twofold` | IUPAC restrictions (incremental) |
| `--three-prime N` | Number of 3' bases that must be perfectly conserved |
| `--inject OLIGO` | Obligatory oligo placed before the search (repeatable / comma-separated). See "Injecting obligatory oligos" |
| `--max-seeds N` | Per-range seed cap in incremental mode (0 = no cap, default 50) |
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
