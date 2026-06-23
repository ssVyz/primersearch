//! Command-line interface. Defaults < settings.ini < CLI flags.

use std::path::PathBuf;

use clap::{Parser, ValueEnum};

use crate::config::ConfigFile;
use crate::engine::{base_mask, Orientation, SearchMode, SearchSettings};

#[derive(Parser, Debug)]
#[command(
    name = "primersearch",
    about = "Find optimal PCR primer regions in aligned DNA sequences.",
    version,
    after_help = "Defaults come from settings.ini next to the executable. \
                  CLI flags always override settings.ini."
)]
pub struct Args {
    /// Input FASTA file (aligned, gap- and ambiguity-free).
    #[arg(required_unless_present = "mkini")]
    pub input: Option<PathBuf>,

    /// Output file.
    #[arg(short = 'o', long, default_value = "output.txt")]
    pub output: PathBuf,

    /// Path to settings.ini (default: next to the executable).
    #[arg(long)]
    pub config: Option<PathBuf>,

    /// Write a default settings.ini (next to the executable, or to --config
    /// if given) and exit. Overwrites any existing file at that path.
    #[arg(long)]
    pub mkini: bool,

    /// Worker threads (0 = use all logical cores).
    #[arg(short = 'j', long, default_value_t = 0)]
    pub threads: usize,

    /// Suppress progress and informational output.
    #[arg(short = 's', long)]
    pub silent: bool,

    // ----- Search parameters -----
    /// Tm threshold (°C).
    #[arg(long, value_name = "C")]
    pub tm: Option<f64>,

    /// Na+ concentration (mM).
    #[arg(long, value_name = "MM")]
    pub na: Option<f64>,

    /// Mg2+ concentration (mM).
    #[arg(long, value_name = "MM")]
    pub mg: Option<f64>,

    /// dNTP concentration (mM). Each dNTP is assumed to sequester one Mg2+.
    #[arg(long, value_name = "MM")]
    pub dntp: Option<f64>,

    /// Total oligo (primer) concentration (µM).
    #[arg(long = "oligo", value_name = "UM")]
    pub oligo_conc: Option<f64>,

    /// Search mode.
    #[arg(long, value_enum)]
    pub mode: Option<CliSearchMode>,

    /// Fixed-slice mode: skip the search and treat the entire input alignment
    /// as one user-chosen slice, generating only the variant(s) needed to
    /// cover every input sequence at that slice. Combine with --mode to pick
    /// exact (no-ambiguities) or IUPAC (incremental) variant generation. The
    /// Tm threshold is not enforced in this mode; each variant's Tm is still
    /// reported.
    #[arg(long)]
    pub fixed: bool,

    /// Target coverage per variant, % (incremental mode).
    #[arg(long, value_name = "PCT")]
    pub target: Option<f64>,

    /// Max ambiguities per variant (incremental mode).
    #[arg(long, value_name = "N")]
    pub max_amb: Option<usize>,

    /// Per-range cap on seeds tried in incremental mode (0 = no cap).
    /// Higher = closer to optimal per-round coverage but slower.
    #[arg(long, value_name = "N")]
    pub max_seeds: Option<usize>,

    /// Forbid N (4-fold) consensus codes (incremental mode).
    #[arg(long)]
    pub exclude_n: bool,
    /// Allow N consensus codes (incremental mode); overrides settings.ini.
    #[arg(long, conflicts_with = "exclude_n")]
    pub no_exclude_n: bool,

    /// Only allow 2-fold (R/Y/S/W/K/M) consensus codes (incremental mode).
    #[arg(long)]
    pub only_twofold: bool,
    /// Allow all-fold consensus codes (incremental mode); overrides settings.ini.
    #[arg(long, conflicts_with = "only_twofold")]
    pub no_only_twofold: bool,

    /// Use reverse (anti-sense) orientation.
    #[arg(long)]
    pub rev: bool,
    /// Use forward (sense) orientation; overrides settings.ini.
    #[arg(long, conflicts_with = "rev")]
    pub fwd: bool,

    /// Number of 3' bases that must be perfectly conserved (0 disables).
    #[arg(long, value_name = "N")]
    pub three_prime: Option<usize>,

    /// Add a space every 3 bases in displayed primer sequences.
    #[arg(long)]
    pub spacer: bool,
    /// Render primer sequences without spacers; overrides settings.ini.
    #[arg(long, conflicts_with = "spacer")]
    pub no_spacer: bool,

    /// Obligatory oligo(s) to use. Each given oligo is placed before the
    /// search runs: it is positioned at the alignment offset where it covers
    /// the most input sequences, emitted as a primer, and the sequences it
    /// covers are removed before any new variants are generated. Repeat the
    /// flag or pass a comma-separated list for several oligos (e.g.
    /// `--inject ACGT... --inject TTGC...` or `--inject ACGT...,TTGC...`).
    /// Oligos may contain IUPAC ambiguity codes. Provide them in the same
    /// orientation as the run (so for `--rev`, the reverse-complement form you
    /// would order). The Tm threshold and ambiguity/3' constraints are not
    /// enforced on injected oligos.
    #[arg(long = "inject", value_name = "OLIGO", value_delimiter = ',')]
    pub inject: Vec<String>,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CliSearchMode {
    NoAmbiguities,
    Incremental,
}

pub struct ResolvedConfig {
    pub settings: SearchSettings,
    pub add_spacer: bool,
}

pub fn resolve(args: &Args, cfg: &ConfigFile) -> ResolvedConfig {
    let mut settings = cfg.to_settings();
    let mut add_spacer = cfg.add_spacer;

    if let Some(v) = args.tm {
        settings.tm_threshold = v;
    }
    if let Some(v) = args.na {
        settings.na_concentration_mm = v;
    }
    if let Some(v) = args.mg {
        settings.mg_concentration_mm = v;
    }
    if let Some(v) = args.dntp {
        settings.dntp_concentration_mm = v;
    }
    if let Some(v) = args.oligo_conc {
        settings.oligo_concentration_um = v;
    }
    if let Some(m) = args.mode {
        settings.mode = match m {
            CliSearchMode::NoAmbiguities => SearchMode::NoAmbiguities,
            CliSearchMode::Incremental => SearchMode::Incremental,
        };
    }
    settings.fixed = args.fixed;
    if let Some(v) = args.target {
        settings.target_coverage_pct = v;
    }
    if let Some(v) = args.max_amb {
        settings.max_ambiguities = v;
    }
    if let Some(v) = args.max_seeds {
        settings.max_seeds = v;
    }
    if args.exclude_n {
        settings.exclude_n = true;
    } else if args.no_exclude_n {
        settings.exclude_n = false;
    }
    if args.only_twofold {
        settings.only_twofold = true;
    } else if args.no_only_twofold {
        settings.only_twofold = false;
    }
    if args.rev {
        settings.orientation = Orientation::Reverse;
    } else if args.fwd {
        settings.orientation = Orientation::Forward;
    }
    if let Some(v) = args.three_prime {
        settings.three_prime_match = v;
    }
    if args.spacer {
        add_spacer = true;
    } else if args.no_spacer {
        add_spacer = false;
    }

    ResolvedConfig { settings, add_spacer }
}

/// Validate and normalise the `--inject` oligo strings against the alignment
/// length. Returns the cleaned oligos (whitespace-stripped, uppercased, in
/// display orientation) ready to hand to the engine, or an error describing
/// the first problem found. Empty input yields an empty vec.
pub fn prepare_injected(raw: &[String], align_len: usize) -> Result<Vec<Vec<u8>>, String> {
    let mut out = Vec::with_capacity(raw.len());
    for (i, s) in raw.iter().enumerate() {
        let oligo: Vec<u8> = s
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .map(|b| b.to_ascii_uppercase())
            .collect();
        if oligo.is_empty() {
            return Err(format!("injected oligo #{} is empty", i + 1));
        }
        if let Some(&bad) = oligo.iter().find(|&&b| base_mask(b) == 0) {
            return Err(format!(
                "injected oligo #{} ({}) contains invalid base '{}': only \
                 A/C/G/T and IUPAC ambiguity codes are allowed",
                i + 1,
                String::from_utf8_lossy(&oligo),
                bad as char,
            ));
        }
        if oligo.len() > align_len {
            return Err(format!(
                "injected oligo #{} (length {}) is longer than the alignment \
                 (length {})",
                i + 1,
                oligo.len(),
                align_len,
            ));
        }
        out.push(oligo);
    }
    Ok(out)
}
