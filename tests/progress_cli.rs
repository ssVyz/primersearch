//! End-to-end checks of `--progress jsonl` against the built binary.

use std::path::PathBuf;
use std::process::{Command, Output};

use serde_json::Value;

/// A 24 bp slice in three clades of 12, 8 and 5 sequences.
fn fixture() -> PathBuf {
    let clades = [
        ("ACGTACGGTCAGTTGACCATGGCA", 12),
        ("ACGTACGGTCTGTTGACCATCGCA", 8),
        ("ACGAACGGTCAGTAGACCATGGCA", 5),
    ];
    let mut fasta = String::new();
    for (ci, (seq, n)) in clades.iter().enumerate() {
        for i in 0..*n {
            fasta.push_str(&format!(">c{ci}_{i}\n{seq}\n"));
        }
    }
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("progress_cli.fasta");
    std::fs::write(&path, fasta).unwrap();
    path
}

fn run(extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_primersearch"))
        .arg(fixture())
        .args(["--format", "json", "--no-config", "--silent", "--fixed"])
        .args(["--mode", "optimize-by-mismatch", "--n-oligos", "2"])
        .args(extra)
        .output()
        .unwrap()
}

fn stderr_lines(out: &Output) -> Vec<Value> {
    String::from_utf8(out.stderr.clone())
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("not JSON: {l:?}: {e}")))
        .collect()
}

#[test]
fn jsonl_leaves_stdout_unchanged_and_reports_phases() {
    let plain = run(&[]);
    let jsonl = run(&["--progress", "jsonl"]);
    assert!(plain.status.success() && jsonl.status.success());
    assert!(plain.stderr.is_empty());
    assert_eq!(plain.stdout, jsonl.stdout);

    let lines = stderr_lines(&jsonl);
    for l in &lines {
        assert_eq!(l["type"], "progress");
        assert!(l["message"].is_string() && l["pct"].is_number());
    }
    let phases: Vec<&str> = lines.iter().filter_map(|l| l["phase"].as_str()).collect();
    for phase in ["windows", "candidates", "reduce", "set_search"] {
        assert!(phases.contains(&phase), "no {phase} line in {phases:?}");
    }
    let search = lines.iter().find(|l| l["phase"] == "set_search").unwrap();
    assert_eq!(search["candidates"], 3);
    assert!(search["max_work"].is_u64());
}

#[test]
fn jsonl_ends_a_failed_run_with_an_error_line() {
    let out = run(&["--progress", "jsonl", "--max-work", "1"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(out.stdout.is_empty());
    let stderr = String::from_utf8(out.stderr).unwrap();
    let last = stderr.lines().last().unwrap();
    let error: Value = serde_json::from_str(last).unwrap();
    assert_eq!(error["type"], "error");
    assert!(error["message"].as_str().unwrap().contains("work limit"));
    // The plain error line stays, right before the JSON one.
    assert!(stderr.contains("\nerror: the set search exceeded the work limit"));
}
