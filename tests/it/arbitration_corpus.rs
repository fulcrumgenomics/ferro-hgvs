//! Corpus-based arbitration integration test (Task 14): proves the honesty
//! property of `arbitrate()` on real, documented ferro/Mutalyzer
//! discrepancies (`tests/fixtures/arbitration.json`), including cases where
//! Mutalyzer is itself spec non-compliant (e.g. a duplication of the wrong
//! base, or emitting `c.` on a non-coding `NR_` transcript). A tool authored
//! by ferro's maintainers claiming "ferro is right" is only trustworthy if
//! it provably never blames ferro on cases the corpus documents as
//! ferro-correct-or-equivalent — that is the invariant this test checks.
//!
//! FERRO_MANIFEST-gated exactly like the manifest-dependent tests added in
//! Task 10 (`tests/it/arbitrate_cli.rs`): the corpus uses real NM_/NR_/LRG_
//! accessions, so arbitrating against them needs a prepared reference. When
//! `FERRO_MANIFEST` is unset (and no `benchmark-output/manifest.json` is
//! present), the test skips rather than fails — same convention as the
//! conformance test suites.
//!
//! `arbitrate()` is called in-process (not via the CLI binary), feeding the
//! corpus's *recorded* `ferro_output` directly — this tests the arbitration
//! logic against the documented pair, decoupled from whatever the current
//! normalizer happens to produce today.

use ferro_hgvs::arbitrate::{arbitrate, ArbitrationCategory, FrameResolver, OtherResult};
use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::{FerroError, HgvsVariant, MultiFastaProvider, VariantProjector};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

/// One row of the documented discrepancy corpus. Only the fields this test
/// consumes are modeled; other fields (`spec_reference`, `timestamp`) are
/// tolerated as pass-through metadata (serde ignores unmodeled fields by
/// default).
#[derive(Debug, Deserialize)]
struct DecisionRow {
    input: String,
    ferro_output: String,
    mutalyzer_output: String,
    category: String,
    reason: String,
}

#[derive(Debug, Deserialize)]
struct Corpus {
    decisions: BTreeMap<String, DecisionRow>,
}

/// `FERRO_MANIFEST`, when set, is authoritative — no fallback to well-known
/// paths. Mirrors `tests/it/arbitrate_cli.rs::manifest_path` exactly, so
/// `FERRO_MANIFEST` can explicitly disable this test even on a host with a
/// reference mounted.
fn manifest_path() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let p = PathBuf::from(path);
        return if p.exists() { Some(p) } else { None };
    }
    let p = PathBuf::from("benchmark-output/manifest.json");
    if p.exists() {
        return Some(p);
    }
    None
}

/// The prepared reference *directory* (parent of `manifest.json`), derived
/// from `manifest_path()`.
fn reference_dir() -> Option<PathBuf> {
    manifest_path().and_then(|p| p.parent().map(|d| d.to_path_buf()))
}

/// Provider type shared by the arbitration's projector: `Arc<MultiFastaProvider>`
/// so it can be cheaply cloned into the projector while `arbitrate()` borrows
/// the same underlying reference (mirrors `run_arbitrate`'s provider idiom in
/// `src/bin/ferro.rs`).
type CorpusProvider = Arc<MultiFastaProvider>;

/// Thin wrapper around a `VariantProjector` so this test crate can implement
/// `ferro_hgvs::arbitrate::FrameResolver` for it. The orphan rule blocks
/// implementing that trait directly for `VariantProjector<CorpusProvider>`
/// here: both the trait and the type are defined in the `ferro_hgvs` library
/// crate, and this integration-test binary is a separate crate — neither is
/// local to it. Same newtype idiom as `CliFrameResolver` in
/// `src/bin/ferro.rs`.
struct CorpusFrameResolver(VariantProjector<CorpusProvider>);

impl FrameResolver for CorpusFrameResolver {
    fn resolve(&self, v: &HgvsVariant) -> Result<HgvsVariant, FerroError> {
        ferro_hgvs::arbitrate::oracle::resolve_frame(v, &self.0)
    }
}

/// Load and deserialize the documented discrepancy corpus.
fn load_corpus() -> Corpus {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/arbitration.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
    serde_json::from_str(&text)
        .unwrap_or_else(|e| panic!("failed to parse {}: {e}", path.display()))
}

#[test]
fn arbitration_corpus_never_falsely_blames_ferro() {
    let Some(reference) = reference_dir() else {
        println!(
            "arbitration_corpus: skipping — no manifest at FERRO_MANIFEST or benchmark-output/manifest.json"
        );
        return;
    };

    // Build the provider + projector (mirrors `run_arbitrate`'s idiom): a
    // `MultiFastaProvider` isn't `Clone`, and `VariantProjector` requires
    // `P: Clone`, so build the concrete provider, extract its owned cdot map,
    // then share it through an `Arc`.
    let manifest_path = reference.join("manifest.json");
    let provider = MultiFastaProvider::from_manifest(&manifest_path)
        .unwrap_or_else(|e| panic!("failed to load reference at {}: {e}", reference.display()));
    let cdot = provider
        .cdot_mapper()
        .cloned()
        .unwrap_or_else(CdotMapper::new);
    let provider: CorpusProvider = Arc::new(provider);
    let projector = VariantProjector::new(Projector::new(cdot), provider.clone());
    let resolver = CorpusFrameResolver(projector);

    let corpus = load_corpus();

    let mut asserted = 0usize;
    let mut skipped = 0usize;
    // Rows where arbitrate reached a *conclusive* verdict — i.e. the oracle
    // actually completed and produced a same/different-class answer rather
    // than `Inconclusive`. Tracked separately so a degraded reference that
    // turned every row into `Inconclusive` cannot vacuously satisfy the
    // one-sided honesty invariant below (the honesty assertion holds trivially
    // for an `Unknown`-category `Inconclusive` row). Note: these corpus rows
    // are genuine-difference cases (Mutalyzer produced a different/wrong
    // variant), so a conclusive verdict here is typically `Different` /
    // `OtherUnparseable`, not `Equivalent` — the deterministic engine leaves
    // the "who is right" call to spec interpretation, which is correct.
    let mut conclusive = 0usize;

    for (key, row) in &corpus.decisions {
        let other = OtherResult {
            tool: "mutalyzer".into(),
            status: "ok".into(),
            output: Some(row.mutalyzer_output.clone()),
        };
        let arb = match arbitrate(
            &row.input,
            &row.ferro_output,
            other,
            provider.as_ref(),
            Some(&resolver),
        ) {
            Ok(arb) => arb,
            Err(e) => {
                // A conversion/reference-coverage gap is inconclusive on this
                // particular reference, not a test failure: the corpus rows
                // reference real accessions that not every prepared
                // reference carries.
                eprintln!("arbitration_corpus: skipping {key} — arbitrate() error: {e}");
                skipped += 1;
                continue;
            }
        };
        asserted += 1;
        if !matches!(arb.verdict, Verdict::Inconclusive) {
            conclusive += 1;
        }

        // The honesty invariant: on a row the corpus documents as
        // ferro-correct-or-equivalent, arbitrate must never conclude that
        // Mutalyzer is right and ferro is wrong (`MutalyzerCorrect`), nor
        // that both are wrong (`BothIncorrect`). It may legitimately land on
        // `Equivalent`, `FerroCorrect`, `Unknown`, or a non-blaming verdict
        // (`Different`/`BasisMismatch`/`OtherUnparseable`).
        if row.category == "ferro_correct" || row.category == "equivalent" {
            assert_ne!(
                arb.category,
                ArbitrationCategory::MutalyzerCorrect,
                "{key}: documented category={:?} (reason: {}), but arbitrate falsely \
                 concluded MutalyzerCorrect — ferro was blamed on a case the corpus \
                 documents as ferro-correct-or-equivalent",
                row.category,
                row.reason
            );
            assert_ne!(
                arb.category,
                ArbitrationCategory::BothIncorrect,
                "{key}: documented category={:?} (reason: {}), but arbitrate falsely \
                 concluded BothIncorrect — ferro was blamed on a case the corpus \
                 documents as ferro-correct-or-equivalent",
                row.category,
                row.reason
            );
        }
    }

    eprintln!("arbitration_corpus: asserted={asserted} skipped={skipped} conclusive={conclusive}");
    assert!(
        asserted > 0,
        "arbitration_corpus: every one of {} corpus row(s) skipped (reference lacks all \
         corpus transcripts) — the honesty invariant could not be exercised on this \
         reference",
        corpus.decisions.len()
    );
    assert!(
        conclusive > 0,
        "arbitration_corpus: every row was Inconclusive ({asserted} asserted, {skipped} \
         skipped) — the reference is too degraded for the oracle to complete on any row, \
         so the honesty invariant would pass vacuously"
    );
}
