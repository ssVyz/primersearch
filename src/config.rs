//! Settings file (TOML syntax, `.ini` extension). Lives next to the
//! executable by default; auto-created with built-in defaults on first run.

use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::engine::{Orientation, SearchMode, SearchSettings};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfigFile {
    pub tm_threshold: f64,
    #[serde(default = "default_oligo_concentration_um")]
    pub oligo_concentration_um: f64,
    pub na_concentration_mm: f64,
    #[serde(default = "default_mg_concentration_mm")]
    pub mg_concentration_mm: f64,
    #[serde(default = "default_dntp_concentration_mm")]
    pub dntp_concentration_mm: f64,
    pub search_mode: ConfigSearchMode,
    pub target_coverage_pct: f64,
    pub max_ambiguities: usize,
    pub exclude_n: bool,
    pub only_twofold: bool,
    pub orientation: ConfigOrientation,
    pub three_prime_match: usize,
    pub add_spacer: bool,
    /// Per-range cap on seeds tried in incremental mode. 0 = no cap.
    #[serde(default = "default_max_seeds")]
    pub max_seeds: usize,
}

fn default_max_seeds() -> usize {
    50
}
fn default_oligo_concentration_um() -> f64 {
    0.2
}
fn default_mg_concentration_mm() -> f64 {
    3.0
}
fn default_dntp_concentration_mm() -> f64 {
    0.8
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConfigSearchMode {
    NoAmbiguities,
    Incremental,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ConfigOrientation {
    Forward,
    Reverse,
}

impl Default for ConfigFile {
    fn default() -> Self {
        Self {
            tm_threshold: 55.0,
            oligo_concentration_um: 0.2,
            na_concentration_mm: 50.0,
            mg_concentration_mm: 3.0,
            dntp_concentration_mm: 0.8,
            search_mode: ConfigSearchMode::NoAmbiguities,
            target_coverage_pct: 50.0,
            max_ambiguities: 1,
            exclude_n: false,
            only_twofold: false,
            orientation: ConfigOrientation::Forward,
            three_prime_match: 0,
            add_spacer: true,
            max_seeds: 50,
        }
    }
}

impl From<ConfigSearchMode> for SearchMode {
    fn from(m: ConfigSearchMode) -> Self {
        match m {
            ConfigSearchMode::NoAmbiguities => SearchMode::NoAmbiguities,
            ConfigSearchMode::Incremental => SearchMode::Incremental,
        }
    }
}

impl From<ConfigOrientation> for Orientation {
    fn from(o: ConfigOrientation) -> Self {
        match o {
            ConfigOrientation::Forward => Orientation::Forward,
            ConfigOrientation::Reverse => Orientation::Reverse,
        }
    }
}

impl ConfigFile {
    pub fn to_settings(&self) -> SearchSettings {
        SearchSettings {
            tm_threshold: self.tm_threshold,
            oligo_concentration_um: self.oligo_concentration_um,
            na_concentration_mm: self.na_concentration_mm,
            mg_concentration_mm: self.mg_concentration_mm,
            dntp_concentration_mm: self.dntp_concentration_mm,
            mode: self.search_mode.into(),
            target_coverage_pct: self.target_coverage_pct,
            max_ambiguities: self.max_ambiguities,
            exclude_n: self.exclude_n,
            only_twofold: self.only_twofold,
            orientation: self.orientation.into(),
            three_prime_match: self.three_prime_match,
            max_seeds: self.max_seeds,
        }
    }
}

pub fn default_config_path() -> PathBuf {
    let exe = std::env::current_exe().unwrap_or_else(|_| PathBuf::from("./primersearch"));
    exe.parent()
        .map(|p| p.join("settings.ini"))
        .unwrap_or_else(|| PathBuf::from("settings.ini"))
}

/// Load `settings.ini` from `path`. If absent, write a fully-populated
/// defaults file there and return those defaults.
pub fn load_or_create(path: &Path) -> std::io::Result<(ConfigFile, bool)> {
    if path.exists() {
        let text = fs::read_to_string(path)?;
        let cfg: ConfigFile = toml::from_str(&text).map_err(|e| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("failed to parse {}: {}", path.display(), e),
            )
        })?;
        Ok((cfg, false))
    } else {
        write_defaults(path)?;
        Ok((ConfigFile::default(), true))
    }
}

/// Write a fully-populated defaults file at `path`, overwriting any
/// existing file. Used by `load_or_create` and by the `--mkini` CLI flag.
pub fn write_defaults(path: &Path) -> std::io::Result<()> {
    let cfg = ConfigFile::default();
    let body = toml::to_string_pretty(&cfg).map_err(|e| {
        std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("failed to serialize defaults: {e}"),
        )
    })?;
    let banner = "# primersearch settings (TOML syntax). Edit values to change\n\
                  # the defaults used when a CLI flag is not supplied. CLI flags\n\
                  # always override values from this file.\n\n";
    fs::write(path, format!("{banner}{body}"))
}
