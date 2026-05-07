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

Constraints applied to each candidate primer: minimum Tm (salt-adjusted
formula), maximum ambiguity count, optional 3' conserved tail, optional
ban on N / non-2-fold codes, forward or reverse orientation.

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
primersearch --mkini                 # write defaults to settings.ini
primersearch input.fasta -j 8        # 8 worker threads
primersearch input.fasta --silent    # no progress / info output
```

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
| `--na MM` | Na⁺ concentration in mM |
| `--mode {no-ambiguities,incremental}` | Search mode |
| `--target PCT` | Coverage target % at which the ambiguity counter is allowed to increase (incremental) |
| `--max-amb N` | Maximum ambiguity codes per primer (incremental) |
| `--exclude-n` / `--only-twofold` | IUPAC restrictions (incremental) |
| `--three-prime N` | Number of 3' bases that must be perfectly conserved |
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
    search.rs        greedy round loop, per-range evaluators, Tm calc
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
