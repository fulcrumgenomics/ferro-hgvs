//! Extract the hermetic reference windows for the mutalyzer **genomic-axis**
//! conformance gate (#737).
//!
//! The mutalyzer `axis_genomic` test projects each transcript-coordinate input
//! onto its genomic (`g.`) frame and normalizes the result — which needs real
//! reference bases (transcript + chromosome) plus the version-independent
//! `NG_`/`LRG_` → chromosome placements (#728) that live only in a multi-GB
//! manifest absent from CI. This generator captures *exactly* what that pass
//! touches into a small committed fixture
//! (`tests/fixtures/mutalyzer-normalize/genomic-windows.json`) so the
//! `gate_genomic_snapshot` test can run hermetically on every PR via
//! [`WindowProvider`].
//!
//! Mirrors `extract_biocommons_windows.rs`: a [`RecordingProvider`] wraps the
//! real manifest provider, the same per-case `project_to_genomic` + `normalize`
//! loop the gate runs is executed, and every transcript resolved, every genomic
//! range read (including the `NG_` bases the #737 re-normalization synthesizes
//! from the chromosome), and every `genomic_placement` consulted is recorded.
//! Transcripts are stored whole; genomic accesses are clustered into disjoint
//! padded windows; placements are stored as `DerivedPlacement` records.
//!
//! Run:   `cargo run --features dev --example extract_mutalyzer_genomic_windows [-- --manifest PATH]`
//! Check: `cargo run --features dev --example extract_mutalyzer_genomic_windows -- --check`

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::PathBuf;
use std::process::ExitCode;
use std::rc::Rc;
use std::sync::Arc;

use clap::Parser;
use ferro_hgvs::conformance::mutalyzer::Fixture;
use ferro_hgvs::conformance::reference_window::{GenomicWindow, WindowFixture};
use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::project::VariantProjector;
use ferro_hgvs::reference::derived_placement::DerivedPlacement;
use ferro_hgvs::reference::provider::GenomicPlacement;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{
    parse_hgvs, FerroError, HgvsVariant, MultiFastaProvider, Normalizer, ReferenceProvider,
};

const CASES_JSON: &str = "tests/fixtures/mutalyzer-normalize/cases.json";
const WINDOWS_JSON: &str = "tests/fixtures/mutalyzer-normalize/genomic-windows.json";

/// Safety margin (bp) added to each side of a captured genomic window beyond the
/// realized access range. Generous relative to any HGVS shuffle distance.
const PAD: u64 = 64;
/// Recorded ranges on one accession closer than this gap are merged into a
/// single window (keeps the fixture small while still splitting distant loci).
const MERGE_GAP: u64 = 2 * PAD;

#[derive(Parser, Debug)]
#[command(about = "Extract hermetic reference windows for the mutalyzer genomic gate (#737)")]
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
    /// Recorded genomic ranges per contig (`get_genomic_sequence`), 0-based half-open.
    contig_ranges: BTreeMap<String, Vec<(u64, u64)>>,
    /// `get_sequence` ranges per id; an id that is neither a transcript nor a
    /// known contig (e.g. a synthesized `NG_`) becomes its own window accession.
    seq_ranges: Vec<(String, u64, u64)>,
    /// Every `NG_`/`LRG_` placement consulted, keyed by parent accession.
    placements: BTreeMap<String, GenomicPlacement>,
}

/// Wraps the real manifest provider, forwarding every call and recording which
/// transcripts, genomic ranges, and placements the pass touches. `Clone` shares
/// the recorder (`Rc<RefCell<_>>`) so a fresh clone handed to each per-case
/// `VariantProjector`/`Normalizer` (taken by value) accumulates into one record.
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
        variant: &HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
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

    fn genomic_placement(
        &self,
        parent: &ferro_hgvs::hgvs::variant::Accession,
    ) -> Option<GenomicPlacement> {
        let p = self.inner.genomic_placement(parent)?;
        self.rec
            .borrow_mut()
            .placements
            .entry(parent.full())
            .or_insert_with(|| p.clone());
        Some(p)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &ferro_hgvs::hgvs::variant::Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        let p = self.inner.genomic_placement_on_build(parent, build)?;
        self.rec
            .borrow_mut()
            .placements
            .entry(parent.full())
            .or_insert_with(|| p.clone());
        Some(p)
    }

    fn resolve_legacy_gene_selector(&self, selector: &str) -> Option<String> {
        self.inner.resolve_legacy_gene_selector(selector)
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
// Genomic pass (mirrors tests/it/mutalyzer_normalize_tests.rs::axis_genomic)
// ----------------------------------------------------------------------------

/// Run one case through the recording provider, ignoring the outcome — we only
/// need the reference accesses, which are recorded as they happen. Panics and
/// errors still leave their accesses recorded.
fn run_case(provider: &RecordingProvider, cdot: &CdotMapper, input: &str) {
    let _ = catch_unwind(AssertUnwindSafe(|| {
        let Ok(v) = parse_hgvs(input) else { return };
        // `project_to_genomic` maps c.→NC through the cdot mapper directly, not
        // `get_transcript`, so capture the resolved transcript explicitly — the
        // gate rebuilds its `CdotMapper` from the fixture's transcripts, which
        // must therefore carry every transcript the gateable cases project on.
        let _ = provider.get_transcript_for_variant(&v);
        let projector = Projector::new(cdot.clone());
        let vp = VariantProjector::new(projector, provider.clone());
        let g = if matches!(
            v,
            HgvsVariant::Cds(_) | HgvsVariant::Tx(_) | HgvsVariant::Rna(_)
        ) {
            match vp.project_to_genomic(&v) {
                Ok(g) => g,
                Err(_) => return,
            }
        } else {
            v
        };
        let normalizer = Normalizer::new(provider.clone());
        let _ = normalizer.normalize(&g);
    }));
}

// ----------------------------------------------------------------------------
// Window assembly
// ----------------------------------------------------------------------------

/// Merge `(start, end)` ranges that overlap or sit within `MERGE_GAP` into
/// padded, disjoint intervals (each clamped to `[0, len]`).
fn merged_windows(mut ranges: Vec<(u64, u64)>, len: Option<u64>) -> Vec<(u64, u64)> {
    ranges.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    for (s, e) in ranges {
        let s = s.saturating_sub(PAD);
        let e = match len {
            Some(l) => (e + PAD).min(l),
            None => e + PAD,
        };
        match out.last_mut() {
            Some(last) if s <= last.1 + MERGE_GAP => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

fn strand_str(s: Strand) -> &'static str {
    match s {
        Strand::Plus => "+",
        Strand::Minus => "-",
        _ => "?",
    }
}

fn build_fixture(provider: &RecordingProvider, captured_from: String) -> WindowFixture {
    let rec = provider.rec.borrow();

    // Cluster every read range by accession. A `get_sequence` id that is not a
    // resolved transcript is itself a windowed accession (a chromosome read via
    // get_sequence, or a synthesized `NG_`): its bases are served back through
    // `get_sequence` too, so capture them the same way.
    let mut acc_ranges: BTreeMap<String, Vec<(u64, u64)>> = rec.contig_ranges.clone();
    for (id, start, end) in &rec.seq_ranges {
        if !rec.transcripts.contains_key(id) {
            acc_ranges
                .entry(id.clone())
                .or_default()
                .push((*start, *end));
        }
    }

    // Keep only the `NM_`/`NR_` transcripts the coverable c. rows project on.
    // An `NG_`/`LRG_` resolved here as a "transcript" is spurious (its placement
    // is captured separately, and the genomic gate maps c.→NC through the cdot
    // exons, not these); dropping them — chiefly the multi-hundred-KB RefSeqGene
    // genomic sequences — keeps the committed fixture small.
    let transcripts: Vec<Transcript> = rec
        .transcripts
        .values()
        .filter(|tx| tx.id.starts_with("NM_") || tx.id.starts_with("NR_"))
        .cloned()
        .collect();

    let mut contig_lengths: BTreeMap<String, u64> = BTreeMap::new();
    let mut genomic: Vec<GenomicWindow> = Vec::new();
    for (acc, ranges) in &acc_ranges {
        let len = provider.inner.get_sequence_length(acc).ok();
        if let Some(l) = len {
            contig_lengths.insert(acc.clone(), l);
        }
        for (start, end) in merged_windows(ranges.clone(), len) {
            // `get_sequence` serves both real chromosome FASTA records and
            // synthesized `NG_`/`LRG_` bases (the #737 placement synthesis),
            // matching exactly how the gate's pass will read them back.
            let bases = provider
                .inner
                .get_sequence(acc, start, end)
                .unwrap_or_else(|e| panic!("capturing window {acc}:{start}-{end}: {e}"));
            genomic.push(GenomicWindow {
                contig: acc.clone(),
                start,
                bases,
            });
        }
    }
    genomic.sort_by(|a, b| (a.contig.as_str(), a.start).cmp(&(b.contig.as_str(), b.start)));

    let mut placements: Vec<DerivedPlacement> = rec
        .placements
        .iter()
        .map(|(parent, p)| DerivedPlacement {
            parent: parent.clone(),
            nc: p.nc.full(),
            nc_start: p.nc_start,
            nc_end: p.nc_end,
            strand: strand_str(p.strand).to_string(),
            anchored_by: String::new(),
            mismatch_fraction: 0.0,
        })
        .collect();
    placements.sort_by(|a, b| a.parent.cmp(&b.parent));

    WindowFixture {
        description: format!(
            "GENERATED by `cargo run --features dev --example \
             extract_mutalyzer_genomic_windows`. Hermetic reference windows for the \
             mutalyzer genomic-axis conformance gate (#737): every transcript, padded \
             genomic window, and NG_/LRG_ placement the project_to_genomic + normalize \
             pass over cases.json touches. Do not hand-edit; regenerate from the \
             manifest. Margin: {PAD} bp."
        ),
        captured_from,
        contig_lengths,
        transcripts,
        genomic,
        placements,
    }
}

// ----------------------------------------------------------------------------
// Manifest resolution + main
// ----------------------------------------------------------------------------

fn resolve_manifest(cli: &Cli) -> Option<PathBuf> {
    if let Some(p) = &cli.manifest {
        return Some(p.clone());
    }
    if let Ok(p) = std::env::var("FERRO_MANIFEST") {
        return Some(PathBuf::from(p));
    }
    None
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> anyhow::Result<ExitCode> {
    let Some(manifest) = resolve_manifest(cli) else {
        eprintln!(
            "no manifest (pass --manifest or set FERRO_MANIFEST); cannot extract genomic windows"
        );
        return Ok(ExitCode::FAILURE);
    };
    let provider = Arc::new(
        MultiFastaProvider::from_manifest(&manifest)
            .map_err(|e| anyhow::anyhow!("load manifest {}: {e}", manifest.display()))?,
    );
    let cdot = provider
        .cdot_mapper()
        .ok_or_else(|| anyhow::anyhow!("manifest has no cdot; cannot project"))?
        .clone();

    let fixture_text = std::fs::read_to_string(CASES_JSON)
        .map_err(|e| anyhow::anyhow!("read {CASES_JSON}: {e}"))?;
    let cases: Fixture = serde_json::from_str(&fixture_text)
        .map_err(|e| anyhow::anyhow!("parse {CASES_JSON}: {e}"))?;

    let recorder = RecordingProvider::new(Arc::clone(&provider));
    let mut n_cases = 0usize;
    for case in &cases.cases {
        if !case.to_test || case.genomic.is_none() {
            continue;
        }
        n_cases += 1;
        run_case(&recorder, &cdot, &case.input);
    }

    let captured_from = format!("{}", manifest.display());
    let fixture = build_fixture(&recorder, captured_from);
    let json = fixture
        .to_json()
        .map_err(|e| anyhow::anyhow!("serialize fixture: {e}"))?;

    if cli.check {
        let committed = std::fs::read_to_string(WINDOWS_JSON).unwrap_or_default();
        if committed == json {
            println!("genomic-windows.json is up to date ({n_cases} genomic cases).");
            return Ok(ExitCode::SUCCESS);
        }
        eprintln!("genomic-windows.json is OUT OF DATE — regenerate with `cargo run --features dev --example extract_mutalyzer_genomic_windows`.");
        return Ok(ExitCode::FAILURE);
    }

    std::fs::write(WINDOWS_JSON, &json)
        .map_err(|e| anyhow::anyhow!("write {WINDOWS_JSON}: {e}"))?;
    println!(
        "Wrote {WINDOWS_JSON}: {} transcripts, {} windows, {} placements (from {n_cases} genomic cases).",
        fixture.transcripts.len(),
        fixture.genomic.len(),
        fixture.placements.len()
    );
    Ok(ExitCode::SUCCESS)
}
