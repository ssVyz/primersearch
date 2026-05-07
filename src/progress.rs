//! CLI progress sink. Wraps an indicatif spinner; in `--silent` mode it's
//! a no-op. Implements the engine's `Progress` trait.

use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};

use crate::engine::Progress;

pub struct CliProgress {
    bar: Option<ProgressBar>,
    cancelled: AtomicBool,
}

impl CliProgress {
    pub fn new(silent: bool) -> Self {
        let bar = if silent {
            None
        } else {
            let style = ProgressStyle::with_template("{spinner:.cyan} {msg}")
                .unwrap_or_else(|_| ProgressStyle::default_spinner());
            let pb = ProgressBar::new_spinner();
            pb.set_style(style);
            pb.enable_steady_tick(Duration::from_millis(120));
            Some(pb)
        };
        Self {
            bar,
            cancelled: AtomicBool::new(false),
        }
    }

    pub fn finish(&self) {
        if let Some(bar) = &self.bar {
            bar.finish_and_clear();
        }
    }
}

impl Progress for CliProgress {
    fn report(&self, message: &str, _pct: f64) {
        if let Some(bar) = &self.bar {
            bar.set_message(message.to_string());
        }
    }
    fn cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }
}
