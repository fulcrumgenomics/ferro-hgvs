//! Extract the hermetic reference windows for the biocommons conformance gate
//! (#325, #478 pillar 4).
//!
//! The biocommons `axis_normalized` test needs real reference bases, which live
//! only in a multi-GB manifest absent from CI. This generator captures *exactly*
//! the bases that pass touches into a small committed fixture
//! (`tests/fixtures/biocommons-normalize/reference-windows.json`), so the gate
//! can run hermetically on every PR via [`WindowProvider`].
//!
//! It does not guess window sizes. A [`RecordingProvider`] wraps the real
//! manifest provider, the same per-case normalize loop the test runs is
//! executed, and every transcript resolved and every genomic range read is
//! recorded. Transcripts are stored whole; genomic accesses are clustered into
//! disjoint windows, each padded by a safety margin, and their bases are
//! captured from the manifest. The realized access set already reflects each
//! case's actual shuffle distance, so the margin is pure safety.
//!
//! Unlike `generate_conformance_summary` (which derives from `cases.json` and so
//! runs in CI), this generator **requires the manifest**, so its `--check` is a
//! local / nightly guard — the per-PR gate consumes the committed fixture and
//! needs no manifest.
//!
//! Run:   `cargo run --features dev --example extract_biocommons_windows [-- --manifest PATH]`
//! Check: `cargo run --features dev --example extract_biocommons_windows -- --check`

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::PathBuf;
use std::process::ExitCode;
use std::rc::Rc;
use std::sync::Arc;

use clap::Parser;
use ferro_hgvs::conformance::biocommons::{Case, Fixture};
use ferro_hgvs::conformance::reference_window::{GenomicWindow, WindowFixture};
use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{
    parse_hgvs, FerroError, MultiFastaProvider, NormalizeConfig, Normalizer, ReferenceProvider,
    ShuffleDirection,
};

const CASES_JSON: &str = "tests/fixtures/biocommons-normalize/cases.json";
const WINDOWS_JSON: &str = "tests/fixtures/biocommons-normalize/reference-windows.json";

/// Safety margin (bp) added to each side of a captured genomic window beyond the
/// realized access range. Generous relative to any HGVS shuffle distance.
const PAD: u64 = 64;
/// Recorded genomic ranges on one contig closer than this gap are merged into a
/// single window (keeps the fixture small when one contig has several loci while
/// still splitting genuinely distant loci into separate windows).
const MERGE_GAP: u64 = 2 * PAD;

#[derive(Parser, Debug)]
#[command(about = "Extract hermetic reference windows for the biocommons conformance gate")]
struct Cli {
    /// Manifest path. Falls back to $FERRO_MANIFEST then well-known locations.
    #[arg(long)]
    manifest: Option<PathBuf>,
    /// Compare the committed fixture to a fresh extraction; exit non-zero on
    /// drift and write nothing. Requires the manifest.
    #[arg(long)]
    check: bool,
}

// ----------------------------------------------------------------------------
// Recording provider
// ----------------------------------------------------------------------------

#[derive(Default)]
struct Recorder {
    /// Every transcript resolved during the pass, keyed by id (deduped).
    transcripts: BTreeMap<String, Transcript>,
    /// Recorded genomic ranges per contig, as `(start, end)` 0-based half-open.
    contig_ranges: BTreeMap<String, Vec<(u64, u64)>>,
    /// `get_sequence` ranges whose id is not (yet) known to be a transcript;
    /// resolved against `transcripts` after the pass — a leftover id is a contig.
    seq_ranges: Vec<(String, u64, u64)>,
}

/// Wraps the real manifest provider, forwarding every call and recording which
/// transcripts and genomic ranges the normalize pass touches. `Clone` is cheap
/// and shares the recorder (`Rc<RefCell<_>>`), so a fresh clone can be handed to
/// each per-case `Normalizer` (which takes its provider by value) while all
/// accesses accumulate into one `Recorder`.
#[derive(Clone)]
struct RecordingProvider {
    inner: Arc<MultiFastaProvider>,
    rec: Rc<RefCell<Recorder>>,
}

impl RecordingProvider {
    fn new(inner: Arc<MultiFastaProvider>) -> Self {
        Self {
            inner,
            rec: Rc::new(RefCell::new(Recorder::default())),
        }
    }

    fn capture(&self, tx: &Arc<Transcript>) {
        self.rec
            .borrow_mut()
            .transcripts
            .entry(tx.id.clone())
            .or_insert_with(|| (**tx).clone());
    }
}

impl ReferenceProvider for RecordingProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        let tx = self.inner.get_transcript(id)?;
        self.capture(&tx);
        Ok(tx)
    }

    fn get_transcript_for_variant(
        &self,
        variant: &ferro_hgvs::HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
        // Delegate to the inner provider's build-aware resolution (#332) so
        // NG/NC-anchored inputs capture the correct transcript.
        let tx = self.inner.get_transcript_for_variant(variant)?;
        self.capture(&tx);
        Ok(tx)
    }

    fn get_transcript_for_accession(
        &self,
        accession: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
        let tx = self.inner.get_transcript_for_accession(accession)?;
        self.capture(&tx);
        Ok(tx)
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        let out = self.inner.get_sequence(id, start, end)?;
        self.rec
            .borrow_mut()
            .seq_ranges
            .push((id.to_string(), start, end));
        Ok(out)
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let out = self.inner.get_genomic_sequence(contig, start, end)?;
        self.rec
            .borrow_mut()
            .contig_ranges
            .entry(contig.to_string())
            .or_default()
            .push((start, end));
        Ok(out)
    }

    fn has_transcript(&self, id: &str) -> bool {
        self.inner.has_transcript(id)
    }

    fn has_transcript_version_exact(&self, id: &str) -> bool {
        self.inner.has_transcript_version_exact(id)
    }

    fn has_genomic_data(&self) -> bool {
        self.inner.has_genomic_data()
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.inner.get_protein_sequence(accession, start, end)
    }

    fn has_protein_data(&self) -> bool {
        self.inner.has_protein_data()
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        self.inner.get_sequence_length(id)
    }
}

// ----------------------------------------------------------------------------
// Normalize loop (mirrors tests/biocommons_normalize_tests.rs::axis_normalized)
// ----------------------------------------------------------------------------

fn direction_from_str(s: &str) -> ShuffleDirection {
    match s {
        "3prime" => ShuffleDirection::ThreePrime,
        "5prime" => ShuffleDirection::FivePrime,
        other => panic!("unknown shuffle_direction in fixture: {other}"),
    }
}

fn build_config(case: &Case) -> NormalizeConfig {
    let mut cfg =
        NormalizeConfig::default().with_direction(direction_from_str(&case.shuffle_direction));
    if case.cross_boundaries {
        cfg = cfg.allow_crossing_boundaries();
    }
    cfg
}

/// Run one case through the recording provider, ignoring the outcome — we only
/// need the reference accesses, which are recorded as they happen. Panics and
/// errors still leave their accesses recorded.
fn run_case(provider: &RecordingProvider, case: &Case) {
    let cfg = if case.expects_error {
        build_config(case).with_error_mode(ErrorMode::Strict)
    } else {
        build_config(case)
    };
    let _ = catch_unwind(AssertUnwindSafe(|| {
        let normalizer = Normalizer::with_config(provider.clone(), cfg);
        if let Ok(v) = parse_hgvs(&case.input) {
            let _ = normalizer.normalize(&v);
        }
    }));
}

// ----------------------------------------------------------------------------
// Window assembly
// ----------------------------------------------------------------------------

/// Merge `(start, end)` ranges that overlap or sit within `MERGE_GAP` of each
/// other into padded, disjoint intervals (each clamped to `[0, contig_len]`).
fn merged_windows(mut ranges: Vec<(u64, u64)>, contig_len: Option<u64>) -> Vec<(u64, u64)> {
    ranges.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    for (s, e) in ranges {
        let s = s.saturating_sub(PAD);
        let e = match contig_len {
            Some(len) => (e + PAD).min(len),
            None => e + PAD,
        };
        match out.last_mut() {
            Some(last) if s <= last.1 + MERGE_GAP => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

fn build_fixture(provider: &RecordingProvider, captured_from: String) -> WindowFixture {
    let rec = provider.rec.borrow();

    // Resolve leftover get_sequence ids: any id not known as a transcript is a
    // contig whose bases were read via get_sequence rather than get_genomic_sequence.
    let mut contig_ranges = rec.contig_ranges.clone();
    for (id, start, end) in &rec.seq_ranges {
        if !rec.transcripts.contains_key(id) {
            contig_ranges
                .entry(id.clone())
                .or_default()
                .push((*start, *end));
        }
    }

    let transcripts: Vec<Transcript> = rec.transcripts.values().cloned().collect();

    let mut contig_lengths: std::collections::BTreeMap<String, u64> =
        std::collections::BTreeMap::new();
    let mut genomic: Vec<GenomicWindow> = Vec::new();
    for (contig, ranges) in &contig_ranges {
        let contig_len = provider.inner.get_sequence_length(contig).ok();
        if let Some(len) = contig_len {
            contig_lengths.insert(contig.clone(), len);
        }
        for (start, end) in merged_windows(ranges.clone(), contig_len) {
            let bases = provider
                .inner
                .get_genomic_sequence(contig, start, end)
                .unwrap_or_else(|e| {
                    panic!("capturing window {contig}:{start}-{end}: {e}");
                });
            genomic.push(GenomicWindow {
                contig: contig.clone(),
                start,
                bases,
            });
        }
    }
    // Deterministic order for stable --check.
    genomic.sort_by(|a, b| (a.contig.as_str(), a.start).cmp(&(b.contig.as_str(), b.start)));

    WindowFixture {
        description: format!(
            "GENERATED by `cargo run --features dev --example extract_biocommons_windows`. \
             Hermetic reference windows for the biocommons conformance gate (#325, #478 \
             pillar 4): every transcript and padded genomic window the normalize pass over \
             cases.json touches. Do not hand-edit; regenerate from the manifest. \
             Margin: {PAD} bp."
        ),
        captured_from,
        contig_lengths,
        transcripts,
        genomic,
        // The biocommons corpus references no NG_/LRG_ parents, so it captures
        // no version-independent placements.
        placements: Vec::new(),
    }
}

// ----------------------------------------------------------------------------
// Manifest resolution + main
// ----------------------------------------------------------------------------

fn resolve_manifest(cli_manifest: Option<PathBuf>) -> Option<PathBuf> {
    if let Some(p) = cli_manifest {
        return p.exists().then_some(p);
    }
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let p = PathBuf::from(path);
        return p.exists().then_some(p);
    }
    // Relative, machine-independent fallback only. Machine-specific locations
    // belong in `--manifest` or `FERRO_MANIFEST`, never hardcoded here.
    let fallback = PathBuf::from("benchmark-output/manifest.json");
    fallback.exists().then_some(fallback)
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    let Some(manifest_path) = resolve_manifest(cli.manifest) else {
        eprintln!(
            "extract_biocommons_windows: no manifest (pass --manifest, set FERRO_MANIFEST, or \
             place it at benchmark-output/manifest.json). This generator requires the full reference manifest."
        );
        return ExitCode::FAILURE;
    };

    let provider = match MultiFastaProvider::from_manifest(&manifest_path) {
        Ok(p) => Arc::new(p),
        Err(e) => {
            eprintln!("from_manifest({}) failed: {e}", manifest_path.display());
            return ExitCode::FAILURE;
        }
    };

    let cases_content = match std::fs::read_to_string(CASES_JSON) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("read {CASES_JSON}: {e}");
            return ExitCode::FAILURE;
        }
    };
    let fixture: Fixture = match serde_json::from_str(&cases_content) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("parse {CASES_JSON}: {e}");
            return ExitCode::FAILURE;
        }
    };

    let recording = RecordingProvider::new(provider);
    for case in &fixture.cases {
        if case.to_test {
            run_case(&recording, case);
        }
    }

    let captured_from = read_prepared_at(&manifest_path).unwrap_or_default();
    let window_fixture = build_fixture(&recording, captured_from);
    let rendered = match window_fixture.to_json() {
        Ok(s) => s,
        Err(e) => {
            eprintln!("serialize window fixture: {e}");
            return ExitCode::FAILURE;
        }
    };

    println!(
        "extract_biocommons_windows: {} transcripts, {} genomic windows ({} bp captured)",
        window_fixture.transcripts.len(),
        window_fixture.genomic.len(),
        window_fixture
            .genomic
            .iter()
            .map(|w| w.bases.len())
            .sum::<usize>(),
    );

    if cli.check {
        let committed = std::fs::read_to_string(WINDOWS_JSON).unwrap_or_default();
        if committed == rendered {
            println!("{WINDOWS_JSON} is up to date");
            return ExitCode::SUCCESS;
        }
        eprintln!(
            "{WINDOWS_JSON} is out of date; regenerate with \
             `cargo run --features dev --example extract_biocommons_windows`"
        );
        return ExitCode::FAILURE;
    }

    if let Err(e) = std::fs::write(WINDOWS_JSON, &rendered) {
        eprintln!("write {WINDOWS_JSON}: {e}");
        return ExitCode::FAILURE;
    }
    println!("wrote {WINDOWS_JSON}");
    ExitCode::SUCCESS
}

/// Read the manifest's `prepared_at` field for provenance, if present.
fn read_prepared_at(manifest_path: &std::path::Path) -> Option<String> {
    let content = std::fs::read_to_string(manifest_path).ok()?;
    let value: serde_json::Value = serde_json::from_str(&content).ok()?;
    value
        .get("prepared_at")
        .and_then(|v| v.as_str())
        .map(str::to_string)
}
