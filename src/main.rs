mod cli;
mod config;
mod engine;
mod output;
mod progress;

use std::fs;
use std::io::{BufWriter, Write};

use clap::Parser;

use crate::cli::Args;
use crate::engine::{find_primers, parse_fasta, quality_filter};
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

    let (cfg, created) = config::load_or_create(&config_path)?;
    let resolved = cli::resolve(&args, &cfg);
    let settings = resolved.settings;

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

    // `input` is required unless `--mkini` is used, which we already handled.
    let input_path = args
        .input
        .as_ref()
        .expect("input is required_unless_present mkini");
    let input_text = fs::read_to_string(input_path)
        .map_err(|e| format!("failed to read {}: {e}", input_path.display()))?;

    let raw = parse_fasta(&input_text);
    let report = quality_filter(raw);

    if !args.silent {
        let stdout = std::io::stdout();
        let mut h = stdout.lock();
        output::format_quality_report(&report, &mut h)?;
        h.flush()?;
    }

    if report.valid_count == 0 {
        return Err("no valid sequences after quality filtering".into());
    }

    let progress = CliProgress::new(args.silent);
    let result = find_primers(&report.valid_sequences, &settings, &progress);
    progress.finish();

    let out_file = fs::File::create(&args.output)
        .map_err(|e| format!("failed to open {} for writing: {e}", args.output.display()))?;
    let mut w = BufWriter::new(out_file);
    output::format_quality_report(&report, &mut w)?;
    output::format_results(&result, &settings, resolved.add_spacer, &mut w)?;
    w.flush()?;

    if !args.silent {
        eprintln!(
            "primersearch: {} primer(s) found, wrote {}",
            result.primers.len(),
            args.output.display()
        );
    }

    Ok(())
}
