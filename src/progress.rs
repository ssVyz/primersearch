//! CLI progress sinks implementing the engine's `Progress` trait: an
//! indicatif spinner (`--progress spinner`, a no-op in `--silent` mode) or a
//! throttled JSON-lines stream on stderr (`--progress jsonl`).

use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};
use std::thread;
use std::time::{Duration, Instant};

use indicatif::{ProgressBar, ProgressStyle};
use serde::ser::{SerializeMap, Serializer};
use serde::Serialize;

use crate::cli::CliProgressMode;
use crate::engine::{Progress, ProgressEvent, ProgressPhase};

/// Minimum interval between two `progress` lines of the same phase.
const MIN_INTERVAL: Duration = Duration::from_millis(200);

/// How often the JSON-lines sink checks for a held-back line that is due.
const FLUSH_POLL: Duration = Duration::from_millis(50);

pub struct CliProgress {
    sink: Sink,
    cancelled: AtomicBool,
}

enum Sink {
    Silent,
    Spinner(ProgressBar),
    Jsonl(Arc<Jsonl>),
}

impl CliProgress {
    pub fn new(mode: CliProgressMode, silent: bool) -> Self {
        let sink = match mode {
            CliProgressMode::Jsonl => Sink::Jsonl(Jsonl::start()),
            CliProgressMode::Spinner if silent => Sink::Silent,
            CliProgressMode::Spinner => {
                let style = ProgressStyle::with_template("{spinner:.cyan} {msg}")
                    .unwrap_or_else(|_| ProgressStyle::default_spinner());
                let pb = ProgressBar::new_spinner();
                pb.set_style(style);
                pb.enable_steady_tick(Duration::from_millis(120));
                Sink::Spinner(pb)
            }
        };
        Self {
            sink,
            cancelled: AtomicBool::new(false),
        }
    }

    pub fn finish(&self) {
        match &self.sink {
            Sink::Silent => {}
            Sink::Spinner(bar) => bar.finish_and_clear(),
            Sink::Jsonl(jsonl) => jsonl.finish(),
        }
    }
}

impl Progress for CliProgress {
    fn report(&self, message: &str, pct: f64) {
        self.report_event(&ProgressEvent {
            message,
            pct,
            phase: None,
            counters: &[],
        });
    }
    fn report_event(&self, event: &ProgressEvent<'_>) {
        match &self.sink {
            Sink::Silent => {}
            Sink::Spinner(bar) => bar.set_message(event.message.to_string()),
            Sink::Jsonl(jsonl) => jsonl.offer(event),
        }
    }
    fn cancelled(&self) -> bool {
        self.cancelled.load(Ordering::Relaxed)
    }
}

/// Write the `error` line of `--progress jsonl` for a failed run.
pub fn emit_error(message: &str) {
    if let Some(line) = error_line(message) {
        write_line(line);
    }
}

/// The JSON-lines sink. Lines are written while holding `state`, so they
/// never interleave and appear in the order the throttle released them. A
/// background thread writes a held-back line once it is due, so the latest
/// report is never withheld for long, however long the next one takes.
struct Jsonl {
    state: Mutex<JsonlState>,
}

struct JsonlState {
    throttle: Throttle,
    finished: bool,
}

impl Jsonl {
    fn start() -> Arc<Self> {
        let jsonl = Arc::new(Jsonl {
            state: Mutex::new(JsonlState {
                throttle: Throttle::default(),
                finished: false,
            }),
        });
        let flusher = Arc::clone(&jsonl);
        thread::spawn(move || {
            loop {
                thread::sleep(FLUSH_POLL);
                let mut st = flusher.lock();
                if st.finished {
                    break;
                }
                if let Some(line) = st.throttle.due(Instant::now()) {
                    write_line(line);
                }
            }
        });
        jsonl
    }

    fn lock(&self) -> MutexGuard<'_, JsonlState> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    fn offer(&self, event: &ProgressEvent<'_>) {
        let Some(line) = progress_line(event) else {
            return;
        };
        let mut st = self.lock();
        if st.finished {
            return;
        }
        if let Some(line) = st.throttle.offer(line, event.phase, Instant::now()) {
            write_line(line);
        }
    }

    /// Write any held-back line and stop, so nothing follows on stderr from
    /// this sink (the `error` line of a failed run comes last).
    fn finish(&self) {
        let mut st = self.lock();
        if let Some(line) = st.throttle.take_pending() {
            write_line(line);
        }
        st.finished = true;
    }
}

/// Rate limit of the JSON-lines sink.
#[derive(Default)]
struct Throttle {
    /// Time and phase of the last line written.
    last: Option<(Instant, Option<ProgressPhase>)>,
    /// The latest line held back, with its phase.
    pending: Option<(String, Option<ProgressPhase>)>,
}

impl Throttle {
    /// Offer `line` of `phase` at `now` and return it if it is to be written
    /// now: always the first line and on a phase change, otherwise at most
    /// one per [`MIN_INTERVAL`]. A line held back replaces the pending one.
    fn offer(
        &mut self,
        line: String,
        phase: Option<ProgressPhase>,
        now: Instant,
    ) -> Option<String> {
        let write = match self.last {
            None => true,
            Some((at, last_phase)) => {
                last_phase != phase || now.saturating_duration_since(at) >= MIN_INTERVAL
            }
        };
        if write {
            self.last = Some((now, phase));
            self.pending = None;
            Some(line)
        } else {
            self.pending = Some((line, phase));
            None
        }
    }

    /// The pending line, once [`MIN_INTERVAL`] has passed since the last
    /// line written.
    fn due(&mut self, now: Instant) -> Option<String> {
        let (at, _) = self.last?;
        if now.saturating_duration_since(at) < MIN_INTERVAL {
            return None;
        }
        let (line, phase) = self.pending.take()?;
        self.last = Some((now, phase));
        Some(line)
    }

    fn take_pending(&mut self) -> Option<String> {
        self.pending.take().map(|(line, _)| line)
    }
}

/// A `progress` line: `type`, `message`, `pct`, then `phase` and the
/// counters where present.
struct ProgressLine<'a>(&'a ProgressEvent<'a>);

impl Serialize for ProgressLine<'_> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let event = self.0;
        let mut map = serializer.serialize_map(None)?;
        map.serialize_entry("type", "progress")?;
        map.serialize_entry("message", event.message)?;
        map.serialize_entry("pct", &event.pct)?;
        if let Some(phase) = event.phase {
            map.serialize_entry("phase", phase.as_str())?;
        }
        for (name, value) in event.counters {
            map.serialize_entry(name, value)?;
        }
        map.end()
    }
}

#[derive(Serialize)]
struct ErrorLine<'a> {
    #[serde(rename = "type")]
    kind: &'static str,
    message: &'a str,
}

fn progress_line(event: &ProgressEvent<'_>) -> Option<String> {
    serde_json::to_string(&ProgressLine(event)).ok()
}

fn error_line(message: &str) -> Option<String> {
    serde_json::to_string(&ErrorLine {
        kind: "error",
        message,
    })
    .ok()
}

/// Write `line` and a newline to stderr in one call and flush. Write errors
/// (e.g. a closed pipe) are ignored: progress must never fail the run.
fn write_line(mut line: String) {
    line.push('\n');
    let mut stderr = std::io::stderr().lock();
    let _ = stderr.write_all(line.as_bytes());
    let _ = stderr.flush();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progress_line_has_documented_keys_in_order() {
        let event = ProgressEvent::new(
            ProgressPhase::Candidates,
            "generated 3/7",
            32.5,
            &[("done", 3), ("total", 7)],
        );
        assert_eq!(
            progress_line(&event).unwrap(),
            r#"{"type":"progress","message":"generated 3/7","pct":32.5,"phase":"candidates","done":3,"total":7}"#
        );
    }

    #[test]
    fn progress_line_omits_absent_keys() {
        let event = ProgressEvent {
            message: "m",
            pct: 0.0,
            phase: None,
            counters: &[],
        };
        assert_eq!(
            progress_line(&event).unwrap(),
            r#"{"type":"progress","message":"m","pct":0.0}"#
        );
    }

    #[test]
    fn lines_stay_on_one_line() {
        let event = ProgressEvent::new(ProgressPhase::SetSearch, "a \"b\"\nc", 70.0, &[]);
        assert!(!progress_line(&event).unwrap().contains('\n'));
        assert_eq!(
            error_line("bad\nthing").unwrap(),
            r#"{"type":"error","message":"bad\nthing"}"#
        );
    }

    #[test]
    fn throttle_writes_phase_changes_and_limits_rate() {
        let t0 = Instant::now();
        let ms = |n| t0 + Duration::from_millis(n);
        let mut t = Throttle::default();
        let mut offer = |s: &str, phase, at| t.offer(s.to_string(), phase, at).is_some();
        assert!(offer("a", Some(ProgressPhase::Windows), ms(0)));
        assert!(offer("b", Some(ProgressPhase::Candidates), ms(1)));
        assert!(!offer("c", Some(ProgressPhase::Candidates), ms(100)));
        assert!(!offer("d", Some(ProgressPhase::Candidates), ms(200)));
        assert!(offer("e", Some(ProgressPhase::Candidates), ms(201)));
        assert!(offer("f", Some(ProgressPhase::Reduce), ms(202)));
        assert!(offer("g", None, ms(203)));
        assert!(!offer("h", None, ms(300)));
    }

    #[test]
    fn throttle_releases_the_latest_held_line_when_due() {
        let t0 = Instant::now();
        let ms = |n| t0 + Duration::from_millis(n);
        let phase = Some(ProgressPhase::Round);
        let mut t = Throttle::default();
        assert!(t.offer("a".into(), phase, ms(0)).is_some());
        assert!(t.offer("b".into(), phase, ms(10)).is_none());
        assert!(t.offer("c".into(), phase, ms(20)).is_none());
        assert_eq!(t.due(ms(150)), None);
        assert_eq!(t.due(ms(200)).as_deref(), Some("c"));
        assert_eq!(t.due(ms(500)), None);
        // The released line restarts the interval.
        assert!(t.offer("d".into(), phase, ms(300)).is_none());
        // A phase change supersedes a held line.
        assert!(t.offer("e".into(), Some(ProgressPhase::Fixed), ms(310)).is_some());
        assert_eq!(t.take_pending(), None);
        assert!(t.offer("f".into(), Some(ProgressPhase::Fixed), ms(320)).is_none());
        assert_eq!(t.take_pending().as_deref(), Some("f"));
    }
}
