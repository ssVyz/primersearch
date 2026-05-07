//! Command-line interface. Defaults < settings.ini < CLI flags.

use std::path::PathBuf;

use clap::{Parser, ValueEnum};

use crate::config::ConfigFile;
use crate::engine::{Orientation, SearchMode, SearchSettings};

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

    /// Search mode.
    #[arg(long, value_enum)]
    pub mode: Option<CliSearchMode>,

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
    if let Some(m) = args.mode {
        settings.mode = match m {
            CliSearchMode::NoAmbiguities => SearchMode::NoAmbiguities,
            CliSearchMode::Incremental => SearchMode::Incremental,
        };
    }
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
