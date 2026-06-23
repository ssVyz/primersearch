mod cli;
mod config;
mod engine;
mod json_output;
mod output;
mod progress;

use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use clap::Parser;

use crate::cli::{Args, CliOutputFormat};
use crate::engine::{
    find_primers, find_primers_fixed, parse_fasta, quality_filter, PrimerSearchResult,
    QualityReport, SearchSettings,
};
use crate::progress::CliProgress;

fn main() {
    let args = Args::parse();
    if let Err(e) = run(args) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    if args.threads > 0 {
        // build_global may fail if the pool is already initialised — for a
        // single-process CLI this only happens if main runs twice, which it
        // doesn't. Suppress the (harmless) error either way.
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let config_path = args
        .config
        .clone()
        .unwrap_or_else(config::default_config_path);

    if args.mkini {
        config::write_defaults(&config_path)?;
        if !args.silent {
            eprintln!(
                "primersearch: wrote defaults to {}",
                config_path.display()
            );
        }
        return Ok(());
    }

    // Resolve configuration. With `--no-config` no settings file is read or
    // created: built-in defaults are overlaid by CLI flags only, so the run is
    // fully deterministic and free of filesystem side effects (handy when
    // driving primersearch as a subprocess).
    let cfg = if args.no_config {
        if !args.silent {
            eprintln!("primersearch: ignoring settings file (--no-config)");
        }
        config::ConfigFile::default()
    } else {
        let (cfg, created) = config::load_or_create(&config_path)?;
        if !args.silent {
            if created {
                eprintln!(
                    "primersearch: created defaults at {}",
                    config_path.display()
                );
            } else {
                eprintln!("primersearch: using config {}", config_path.display());
            }
        }
        cfg
    };
    let resolved = cli::resolve(&args, &cfg);
    let settings = resolved.settings;

    // `input` is required unless `--mkini` is used, which we already handled.
    let input_path = args
        .input
        .as_ref()
        .expect("input is required_unless_present mkini");
    let input_text = fs::read_to_string(input_path)
        .map_err(|e| format!("failed to read {}: {e}", input_path.display()))?;

    let raw = parse_fasta(&input_text);
    let report = quality_filter(raw);

    // Decide where output goes. Text mode defaults to a file (`output.txt`);
    // JSON mode defaults to stdout. `-o -` forces stdout in either mode.
    let to_stdout = match args.output.as_deref() {
        Some(p) => p == Path::new("-"),
        None => matches!(args.format, CliOutputFormat::Json),
    };

    // The human-readable preprocessing echo to stdout (so an interactive user
    // can abort early on bad input) only applies to text mode writing to a
    // file. In JSON mode stdout is reserved for the JSON document; when text
    // output itself goes to stdout the report is already part of it.
    if !args.silent && matches!(args.format, CliOutputFormat::Text) && !to_stdout {
        let stdout = std::io::stdout();
        let mut h = stdout.lock();
        output::format_quality_report(&report, &mut h)?;
        h.flush()?;
    }

    if report.valid_count == 0 {
        return Err("no valid sequences after quality filtering".into());
    }

    // Validate any user-supplied obligatory oligos against the alignment
    // length now that we know the valid sequences share one length.
    let align_len = report.valid_sequences[0].len();
    let injected = cli::prepare_injected(&args.inject, align_len)?;
    if !args.silent && !injected.is_empty() {
        eprintln!("primersearch: injecting {} obligatory oligo(s)", injected.len());
    }

    let progress = CliProgress::new(args.silent);
    let result = if settings.fixed {
        find_primers_fixed(&report.valid_sequences, &settings, &injected, &progress)
    } else {
        find_primers(&report.valid_sequences, &settings, &injected, &progress)
    };
    progress.finish();

    let dest_label = if to_stdout {
        let stdout = std::io::stdout();
        let mut w = stdout.lock();
        write_output(&args, &report, &result, &settings, resolved.add_spacer, &mut w)?;
        w.flush()?;
        "<stdout>".to_string()
    } else {
        // In text mode the default destination is `output.txt`; JSON mode only
        // reaches this branch when an explicit `-o <file>` was given.
        let path = args.output.as_deref().unwrap_or_else(|| Path::new("output.txt"));
        let out_file = fs::File::create(path)
            .map_err(|e| format!("failed to open {} for writing: {e}", path.display()))?;
        let mut w = BufWriter::new(out_file);
        write_output(&args, &report, &result, &settings, resolved.add_spacer, &mut w)?;
        w.flush()?;
        path.display().to_string()
    };

    if !args.silent {
        eprintln!(
            "primersearch: {} primer(s) found, wrote {}",
            result.primers.len(),
            dest_label
        );
    }

    Ok(())
}

/// Render the run to `w` in the format selected by `--format`.
fn write_output<W: Write>(
    args: &Args,
    report: &QualityReport,
    result: &PrimerSearchResult,
    settings: &SearchSettings,
    add_spacer: bool,
    w: &mut W,
) -> std::io::Result<()> {
    match args.format {
        CliOutputFormat::Text => {
            output::format_quality_report(report, w)?;
            output::format_results(result, settings, add_spacer, w)?;
        }
        CliOutputFormat::Json => {
            json_output::write_json(report, result, settings, w)?;
        }
    }
    Ok(())
}
