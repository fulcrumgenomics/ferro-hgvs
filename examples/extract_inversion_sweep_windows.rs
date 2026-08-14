//! Capture both committed artifacts behind the authored-inversion gate: the
//! pinned outputs (`cases.tsv`) and the hermetic reference slice
//! (`reference-windows.json`).
//!
//! # Why the gate cannot be manifest-backed
//!
//! Every row is sequence-dependent — the 3'rule reads the bases it would shift
//! over, and the retype-to-`=` outcome reads the span to discover it is its own
//! reverse complement. Running that under the prepared reference would make the
//! whole suite a `FERRO_MANIFEST`-or-skip test, and in CI that branch is always
//! "skip": a suite that reports PASS having evaluated nothing. So the usual
//! pattern applies — run the real pass against a real manifest through a
//! [`RecordingProvider`], and commit precisely the slice it touched.
//!
//! # `serve_transcripts_whole` is load-bearing, not a size tweak
//!
//! A window captured from a pass is exactly as wide as that pass happened to
//! read, which makes the fixture self-fulfilling: it serves the recorded pass
//! and fails any pass that reads further. `src/normalize/mod.rs` then turns that
//! failure into `Ok((input-unchanged, vec![]))` — success, empty warnings, input
//! returned. For this corpus that is catastrophic rather than merely wrong:
//! `Unchanged` is the *majority* outcome, so an under-serving slice reproduces
//! the pinned answer for three quarters of the rows by accident. Widening every
//! transcript's window to the whole stored sequence removes the failure mode by
//! construction, and `tests/it/inversion_sweep.rs` asserts the widening
//! hermetically plus zero failed reads through
//! [`AuditedProvider`](ferro_hgvs::conformance::audited_provider::AuditedProvider).
//!
//! Run:   `cargo run --features dev --example extract_inversion_sweep_windows -- --manifest PATH`
//! Check: `cargo run --features dev --example extract_inversion_sweep_windows -- --check --manifest PATH`

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use clap::Parser;
use ferro_hgvs::conformance::completeness::CaptureLedger;
use ferro_hgvs::conformance::inversion_sweep::{
    classify, sweep_inputs, sweep_normalize_config, PinnedCase, PinnedCases, ADJUDICATED_INPUTS,
    CASES_PATH, TRANSCRIPT, WINDOWS_PATH,
};
use ferro_hgvs::conformance::reference_window::{GenomicWindow, WindowFixture};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::{MultiFastaProvider, Normalizer};

#[path = "common/recording.rs"]
mod recording;

use recording::{build_window_fixture, RecordingProvider};

/// Named in every failure whose remedy is a regeneration, so a reader never has
/// to reconstruct the command.
const REGENERATE: &str = "cargo run --features dev --example extract_inversion_sweep_windows -- \
                          --manifest <manifest>";

/// Anchor a repo-relative artifact path on `CARGO_MANIFEST_DIR`.
///
/// Run from a subdirectory, a bare relative `CASES_PATH` resolves against the
/// working directory instead: `--check` then reports drift against a file that
/// does not exist, and a regeneration writes a stray fixture tree. The gate
/// side already anchors the same way (`tests/it/inversion_sweep.rs`).
fn repo_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

/// The corpus the generator pins, as the descriptions handed to the normalizer.
///
/// [`run`] iterates this rather than enumerating inline, so the generator and
/// the gate cannot drift onto different corpora — and an inlined replacement in
/// `run` would leave this unused, which `-D warnings` rejects rather than
/// silently accepts.
fn corpus() -> Vec<String> {
    sweep_inputs()
        .iter()
        .map(|input| input.description())
        .collect()
}

#[derive(Parser, Debug)]
#[command(about = "Capture the pinned outputs and hermetic slice for the authored-inversion gate")]
struct Cli {
    /// Prepared-reference manifest. Falls back to `$FERRO_MANIFEST`.
    #[arg(long)]
    manifest: Option<PathBuf>,
    /// Compare both committed artifacts to a fresh extraction; exit non-zero on
    /// drift and write nothing. Requires the manifest.
    #[arg(long)]
    check: bool,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    // Refuse a `FERRO_PARTITION` naming no arm before any work, so a
    // misspelling cannot be served by the shipped rule under a candidate's
    // name. See `tests/it/partition_switch_wiring.rs`.
    if let Some(message) = ferro_hgvs::normalize::partition_switch_startup_error() {
        eprintln!("error: {message}");
        return ExitCode::FAILURE;
    }
    match run(&cli) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> anyhow::Result<ExitCode> {
    let Some(manifest) = cli
        .manifest
        .clone()
        .or_else(|| std::env::var("FERRO_MANIFEST").ok().map(PathBuf::from))
    else {
        eprintln!("no manifest (pass --manifest or set FERRO_MANIFEST); cannot extract");
        return Ok(ExitCode::FAILURE);
    };

    let provider = Arc::new(
        MultiFastaProvider::from_manifest(&manifest)
            .map_err(|e| anyhow::anyhow!("load manifest {}: {e}", manifest.display()))?,
    );
    let recorder = RecordingProvider::new(Arc::clone(&provider));

    let mut rows = Vec::new();
    let mut census: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut non_conformant = Vec::new();
    let mut not_a_fixed_point = Vec::new();
    // Every row of the corpus is routed through the ledger (#1550/#1551), so a
    // pass that pinned fewer rows than it enumerated cannot write an artifact
    // that looks complete. `normalize` aborts rather than skipping, so the
    // shortfall this catches is the *silent* one — an empty or truncated
    // `sweep_inputs`, which would otherwise produce a well-formed `cases.tsv`
    // holding nothing and a census of zeros.
    let mut ledger = CaptureLedger::new("authored inversions");
    for description in corpus() {
        let Some((normalized, rendered)) =
            ledger.record(description.as_str(), normalize(&recorder, &description))
        else {
            continue;
        };
        // The gate asserts every answer is a fixed point, so the reads of
        // re-normalizing the answer are part of what the slice must serve — and
        // the *result* is checked here rather than discarded. Leaving it to
        // `every_pinned_output_is_a_fixed_point` means a careless regeneration
        // commits the unstable row first and the gate fails on it afterwards,
        // with the committed artifact — the thing a reader trusts — asserting
        // that the violation is expected. Same precedent as the conformance
        // refusal below: refuse to write it, do not report it later.
        let (_, again) = normalize(&recorder, &rendered)?;
        if again != rendered {
            not_a_fixed_point.push(format!("  {description} -> {rendered} -> {again}"));
        }
        let outcome = classify(&description, &normalized, &rendered);
        if !outcome.is_conformant() {
            non_conformant.push(format!("  {description} -> {rendered}  [{outcome:?}]"));
        }
        *census.entry(outcome.label()).or_default() += 1;
        rows.push(PinnedCase {
            input: description,
            output: rendered,
        });
    }
    let captured = ledger.finish()?;
    // Refuse to write a fixture that already encodes a violation. Otherwise a
    // careless regeneration commits the split or the retype, the invariant gate
    // fails on it, and the committed artifact — the thing a reader trusts —
    // says the violation is expected. The transcript-coverage check below sets
    // the same precedent: assert, do not print.
    if !non_conformant.is_empty() {
        anyhow::bail!(
            "{} row(s) normalized outside the inversion family; refusing to pin them. \
             `DNA/inversion.md:5` defines the replacement of the whole span as its reverse \
             complement, and a split into members or a retype asserts something else:\n{}",
            non_conformant.len(),
            non_conformant.join("\n")
        );
    }
    if !not_a_fixed_point.is_empty() {
        anyhow::bail!(
            "{} row(s) do not re-normalize to themselves; refusing to pin them. A pinned \
             answer that is not a fixed point makes `every_pinned_output_is_a_fixed_point` \
             fail on the committed artifact, which then reads as though the instability were \
             the expected behaviour:\n{}",
            not_a_fixed_point.len(),
            not_a_fixed_point.join("\n")
        );
    }

    // The three spec adjudications are pinned in the gate's source rather than
    // in `cases.tsv` — but their reads still have to be in the slice, or the
    // gate reports "unchanged" for them because the provider ran dry.
    for input in ADJUDICATED_INPUTS {
        let (_, rendered) = normalize(&recorder, input)?;
        normalize(&recorder, &rendered)?;
        println!("adjudicated: {input} -> {rendered}");
    }

    // `WindowFixture::captured_from` is documented as the manifest's
    // `prepared_at`, and that is the field worth recording: a bare file name is
    // `manifest.json` for every prepared reference ever built, so two different
    // references produce provenance a later `--check` cannot tell apart —
    // "the reference moved" reads identically to "ferro moved". The timestamp
    // identifies the reference without leaking the local path, which is what
    // ruled the absolute path out. Same treatment as
    // `extract_biocommons_windows.rs` and `build_conformance_snapshot.rs`.
    let captured_from = read_prepared_at(&manifest).unwrap_or_else(|| {
        manifest
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "manifest.json".to_string())
    });
    let description = format!(
        "GENERATED by `cargo run --features dev --example extract_inversion_sweep_windows`. \
         Hermetic reference slice for the authored-inversion gate: {TRANSCRIPT} served WHOLE, \
         plus any padded genomic window the normalize pass over the sweep touches. Do not \
         hand-edit; regenerate from a prepared-reference manifest. Each recorded range is \
         padded by {} bp on each side, and padded ranges whose gap is within {} bp are merged \
         into one window — so a window may span bases the pass never read, and its extent is \
         not a claim about what was touched.",
        recording::PAD,
        recording::MERGE_GAP
    );

    // mRNA-sized records only: an `NG_`/`LRG_` genomic record resolved as a
    // "transcript" carries hundreds of KB, and the pre-commit hook rejects an
    // added file over 500 KB.
    let mut fixture = build_window_fixture(&recorder, captured_from, description, |tx| {
        tx.id.starts_with("NM_") || tx.id.starts_with("NR_")
    });
    serve_transcripts_whole(&mut fixture);

    // Coverage, asserted rather than printed. A printed count is read by a
    // human who already believes the run worked, so a slice that lost the one
    // transcript the whole corpus is authored on would pass `--check` quietly —
    // certifying exactly the truncated fixture this generator exists to rule
    // out.
    if !fixture.transcripts.iter().any(|tx| tx.id == TRANSCRIPT) {
        anyhow::bail!(
            "the pass never resolved {TRANSCRIPT}, so the slice would under-serve every row of \
             the sweep. Check the manifest carries that transcript, then re-run."
        );
    }

    let cases = PinnedCases { rows };
    let rendered_cases = cases.to_tsv(&format!(
        "GENERATED by `{REGENERATE}` — do not hand-edit.\n\
         The authored-inversion sweep: {} rows of `{TRANSCRIPT}:c.<start>_<end>inv`, each pinned \
         to the exact string `ferro normalize --error-mode lenient` produces for it.\n\
         Outcome census: {}.\n\
         Columns: authored input <TAB> normalized output.",
        cases.rows.len(),
        census
            .iter()
            .map(|(label, count)| format!("{count} {label}"))
            .collect::<Vec<_>>()
            .join(", "),
    ));
    let rendered_windows = fixture
        .to_json()
        .map_err(|e| anyhow::anyhow!("serialize fixture: {e}"))?;

    println!(
        "captured {}/{} rows ({} dropped); census: {}; slice: {} transcript(s), {} window(s), {} bp",
        captured.succeeded,
        captured.attempted,
        captured.dropped,
        census
            .iter()
            .map(|(label, count)| format!("{count} {label}"))
            .collect::<Vec<_>>()
            .join(", "),
        fixture.transcripts.len(),
        fixture.genomic.len(),
        fixture.genomic.iter().map(|w| w.bases.len()).sum::<usize>(),
    );

    let cases_path = repo_path(CASES_PATH);
    let windows_path = repo_path(WINDOWS_PATH);

    if cli.check {
        let mut stale = Vec::new();
        for (relative, path, fresh) in [
            (CASES_PATH, &cases_path, &rendered_cases),
            (WINDOWS_PATH, &windows_path, &rendered_windows),
        ] {
            if std::fs::read_to_string(path).ok().as_ref() != Some(fresh) {
                stale.push(relative);
            }
        }
        if !stale.is_empty() {
            eprintln!("{stale:?} out of date — regenerate with `{REGENERATE}`.");
            return Ok(ExitCode::FAILURE);
        }
        println!("both committed artifacts are up to date.");
        return Ok(ExitCode::SUCCESS);
    }

    let dir = cases_path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("{} has no parent directory", cases_path.display()))?;
    std::fs::create_dir_all(dir).map_err(|e| anyhow::anyhow!("create {}: {e}", dir.display()))?;
    std::fs::write(&cases_path, &rendered_cases)
        .map_err(|e| anyhow::anyhow!("write {}: {e}", cases_path.display()))?;
    std::fs::write(&windows_path, &rendered_windows)
        .map_err(|e| anyhow::anyhow!("write {}: {e}", windows_path.display()))?;
    println!(
        "wrote {} and {}",
        cases_path.display(),
        windows_path.display()
    );
    Ok(ExitCode::SUCCESS)
}

/// Read the manifest's `prepared_at` field for provenance, if present.
fn read_prepared_at(manifest: &Path) -> Option<String> {
    let content = std::fs::read_to_string(manifest).ok()?;
    let value: serde_json::Value = serde_json::from_str(&content).ok()?;
    value
        .get("prepared_at")
        .and_then(|v| v.as_str())
        .map(str::to_string)
}

/// Normalize `input` through the shipped lenient path — the same entry point
/// `tests/it/inversion_sweep.rs` replays — returning the variant and its
/// rendering.
///
/// Every failure aborts. A row that does not parse is a generator bug, and a
/// row that does not normalize would have to be pinned as an error string,
/// which this corpus has no shape for: the sweep is 2,075 machine-generated
/// descriptions on one transcript, so a single failure means the reference is
/// wrong rather than the row.
fn normalize(
    recorder: &RecordingProvider,
    input: &str,
) -> anyhow::Result<(ferro_hgvs::HgvsVariant, String)> {
    let (config, normalize_config) = sweep_normalize_config();
    let parsed = parse_hgvs_with_config(input, config)
        .map_err(|e| anyhow::anyhow!("{input:?} does not parse: {e}"))?;
    let normalized = Normalizer::with_config(recorder.clone(), normalize_config)
        .normalize(&parsed.result)
        .map_err(|e| anyhow::anyhow!("{input:?} does not normalize: {e}"))?;
    let rendered = normalized.to_string();
    Ok((normalized, rendered))
}

/// Widen every window whose accession is also a stored transcript so it covers
/// that transcript whole.
///
/// This is a fidelity fix. The normalizer reaches a transcript through
/// `get_genomic_sequence` as well as `get_transcript`, so the recorder files
/// those reads under `contig_ranges` and the accession lands in the fixture's
/// `genomic` map. That makes `WindowProvider::is_known_contig` true for it,
/// which routes **`get_sequence` for the transcript** through the
/// genomic-window path instead of the whole stored sequence — so every read
/// outside the captured window errors, even though the bases are sitting right
/// there in `transcripts`.
///
/// Serving the accession whole removes the failure mode by construction.
fn serve_transcripts_whole(fixture: &mut WindowFixture) {
    for transcript in &fixture.transcripts {
        let Some(bases) = transcript.sequence.clone() else {
            continue;
        };
        fixture
            .genomic
            .retain(|window| window.contig != transcript.id);
        fixture.genomic.push(GenomicWindow {
            contig: transcript.id.clone(),
            start: 0,
            bases,
        });
        fixture
            .contig_lengths
            .insert(transcript.id.clone(), transcript.sequence_length());
    }
    // Deterministic order for a stable `--check`.
    fixture
        .genomic
        .sort_by(|a, b| (a.contig.as_str(), a.start).cmp(&(b.contig.as_str(), b.start)));
}

#[cfg(test)]
mod tests {
    use super::{corpus, serve_transcripts_whole};
    use ferro_hgvs::conformance::inversion_sweep::sweep_inputs;
    use ferro_hgvs::conformance::reference_window::{GenomicWindow, WindowFixture};
    use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};

    fn transcript(id: &str, bases: &str) -> Transcript {
        Transcript::new(
            id.to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            bases.to_string(),
            Some(1),
            Some(bases.len() as u64),
            vec![Exon::new(1, 1, bases.len() as u64)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        )
    }

    #[test]
    fn a_narrow_transcript_window_is_replaced_by_the_whole_sequence() {
        // The measured defect: the recorder captured a 4 bp window on the
        // transcript's own accession, which is what `WindowProvider` would then
        // serve `get_sequence` from — erroring on every read outside it, which
        // normalization converts into "input unchanged", which for this corpus
        // is the majority pinned answer.
        let mut fixture = WindowFixture {
            transcripts: vec![transcript("NM_TEST.1", "ACGTACGTACGT")],
            genomic: vec![GenomicWindow {
                contig: "NM_TEST.1".to_string(),
                start: 4,
                bases: "ACGT".to_string(),
            }],
            ..Default::default()
        };
        serve_transcripts_whole(&mut fixture);
        assert_eq!(fixture.genomic.len(), 1);
        assert_eq!(fixture.genomic[0].start, 0);
        assert_eq!(fixture.genomic[0].bases, "ACGTACGTACGT");
        assert_eq!(fixture.contig_lengths.get("NM_TEST.1"), Some(&12));
    }

    #[test]
    fn a_genomic_window_on_another_accession_is_left_alone() {
        let mut fixture = WindowFixture {
            transcripts: vec![transcript("NM_TEST.1", "ACGT")],
            genomic: vec![GenomicWindow {
                contig: "NC_000023.11".to_string(),
                start: 100,
                bases: "TTTT".to_string(),
            }],
            ..Default::default()
        };
        serve_transcripts_whole(&mut fixture);
        let contigs: Vec<&str> = fixture.genomic.iter().map(|w| w.contig.as_str()).collect();
        assert_eq!(contigs, ["NC_000023.11", "NM_TEST.1"]);
    }

    #[test]
    fn the_generator_and_the_gate_enumerate_one_corpus() {
        // `run` iterates `corpus()`, so this observes what the generator
        // actually writes rather than re-deriving it: the descriptions are
        // pinned as literal strings, and an inlined replacement inside `run`
        // would leave `corpus` unused, which `-D warnings` rejects.
        let corpus = corpus();
        assert_eq!(corpus.len(), 2075);
        assert_eq!(corpus.first().unwrap(), "NM_004006.2:c.101_104inv");
        assert_eq!(corpus.last().unwrap(), "NM_004006.2:c.2999_3014inv");
        let from_sweep: Vec<String> = sweep_inputs()
            .iter()
            .map(|input| input.description())
            .collect();
        assert_eq!(corpus, from_sweep, "the generator pins a different corpus");
    }
}
