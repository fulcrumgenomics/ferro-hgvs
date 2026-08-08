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
use ferro_hgvs::conformance::inversion_sweep::{
    classify, sweep_inputs, PinnedCase, PinnedCases, ADJUDICATED_INPUTS, CASES_PATH, TRANSCRIPT,
    WINDOWS_PATH,
};
use ferro_hgvs::conformance::reference_window::{GenomicWindow, WindowFixture};
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::{MultiFastaProvider, NormalizeConfig, Normalizer, ShuffleDirection};

#[path = "common/recording.rs"]
mod recording;

use recording::{build_window_fixture, RecordingProvider};

/// Named in every failure whose remedy is a regeneration, so a reader never has
/// to reconstruct the command.
const REGENERATE: &str = "cargo run --features dev --example extract_inversion_sweep_windows -- \
                          --manifest <manifest>";

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
    for input in sweep_inputs() {
        let description = input.description();
        let (normalized, rendered) = normalize(&recorder, &description)?;
        // The gate asserts every answer is a fixed point, so the reads of
        // re-normalizing the answer are part of what the slice must serve.
        normalize(&recorder, &rendered)?;
        *census
            .entry(classify(&description, &normalized, &rendered).label())
            .or_default() += 1;
        rows.push(PinnedCase {
            input: description,
            output: rendered,
        });
    }

    // The three spec adjudications are pinned in the gate's source rather than
    // in `cases.tsv` — but their reads still have to be in the slice, or the
    // gate reports "unchanged" for them because the provider ran dry.
    for input in ADJUDICATED_INPUTS {
        let (_, rendered) = normalize(&recorder, input)?;
        normalize(&recorder, &rendered)?;
        println!("adjudicated: {input} -> {rendered}");
    }

    // Record only the manifest's file name: an absolute path would leak the
    // local layout into the repo and churn `--check` per machine.
    let captured_from = manifest
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "manifest.json".to_string());
    let description = format!(
        "GENERATED by `cargo run --features dev --example extract_inversion_sweep_windows`. \
         Hermetic reference slice for the authored-inversion gate: {TRANSCRIPT} served WHOLE, \
         plus any padded genomic window the normalize pass over the sweep touches. Do not \
         hand-edit; regenerate from a prepared-reference manifest. Margin: {} bp.",
        recording::PAD
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
        "{} rows; census: {}; slice: {} transcript(s), {} window(s), {} bp",
        cases.rows.len(),
        census
            .iter()
            .map(|(label, count)| format!("{count} {label}"))
            .collect::<Vec<_>>()
            .join(", "),
        fixture.transcripts.len(),
        fixture.genomic.len(),
        fixture.genomic.iter().map(|w| w.bases.len()).sum::<usize>(),
    );

    if cli.check {
        let mut stale = Vec::new();
        for (path, fresh) in [
            (CASES_PATH, &rendered_cases),
            (WINDOWS_PATH, &rendered_windows),
        ] {
            if std::fs::read_to_string(path).ok().as_ref() != Some(fresh) {
                stale.push(path);
            }
        }
        if !stale.is_empty() {
            eprintln!("{stale:?} out of date — regenerate with `{REGENERATE}`.");
            return Ok(ExitCode::FAILURE);
        }
        println!("both committed artifacts are up to date.");
        return Ok(ExitCode::SUCCESS);
    }

    let dir = Path::new(CASES_PATH)
        .parent()
        .ok_or_else(|| anyhow::anyhow!("{CASES_PATH} has no parent directory"))?;
    std::fs::create_dir_all(dir).map_err(|e| anyhow::anyhow!("create {}: {e}", dir.display()))?;
    std::fs::write(CASES_PATH, &rendered_cases)
        .map_err(|e| anyhow::anyhow!("write {CASES_PATH}: {e}"))?;
    std::fs::write(WINDOWS_PATH, &rendered_windows)
        .map_err(|e| anyhow::anyhow!("write {WINDOWS_PATH}: {e}"))?;
    println!("wrote {CASES_PATH} and {WINDOWS_PATH}");
    Ok(ExitCode::SUCCESS)
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
    let config = ErrorConfig::lenient();
    let parsed = parse_hgvs_with_config(input, config.clone())
        .map_err(|e| anyhow::anyhow!("{input:?} does not parse: {e}"))?;
    let normalize_config = NormalizeConfig::for_entry_point(ShuffleDirection::ThreePrime, config);
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
    use super::serve_transcripts_whole;
    use ferro_hgvs::conformance::inversion_sweep::{sweep_inputs, SweepOutcome};
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
        // Both sides call `sweep_inputs`; this pins that the shared enumeration
        // is what the generator writes, so a future edit that inlines the loop
        // here shows up as a failure rather than as a silently different corpus.
        assert_eq!(sweep_inputs().len(), 2075);
        assert!(SweepOutcome::Unchanged.is_conformant());
    }
}
