# Changelog

All notable changes to this project are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning: [SemVer](https://semver.org/).

## [Unreleased]

## [0.1.1] - 2026-10-02

### Added
- `--progress {spinner,jsonl}`. `jsonl` writes throttled, machine-readable
  progress lines (phase and counters) to stderr, also with `--silent`, and
  ends a failed run with a `{"type":"error",…}` line. Default output is
  unchanged.
- Engine: `ProgressEvent`, `ProgressPhase` and `Progress::report_event`
  (the default forwards to `report`).
- Optimize-by-mismatch progress during the cross-window pool reduction, when
  placing injected oligos, and when candidate generation starts.
- End-to-end tests in `tests/progress_cli.rs`.

### Changed
- The optimize-by-mismatch set search reports progress about every 250 ms
  instead of every 2^22 evaluations, so reports keep coming on large inputs.

## [0.1.0] - 2026-10-02

### Added
- Initial versioned baseline. Existing features: greedy search modes
  (`no-ambiguities`, `incremental`), exhaustive `optimize-by-mismatch` mode,
  fixed-slice mode, nearest-neighbor Tm, oligo injection (`--inject`),
  primer exclusion (`--exclude`), JSON output, `--no-config`, parallel
  evaluation.
