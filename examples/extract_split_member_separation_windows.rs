//! Extract the hermetic reference windows for the issue #1539 separation guard.
//!
//! `tests/it/issue_1539_split_member_separation.rs` adjudicates a real
//! normalize pass over `tests/fixtures/split-member-separation/inputs.txt`, so
//! it needs real transcript bases — which live only in a multi-GB prepared
//! reference CI does not have. This generator runs *exactly* the pass the guard
//! runs (parse, normalize, then resolve each emitted member to its SPDI
//! reference/payload pair) against a real manifest through a
//! [`RecordingProvider`], and commits precisely the slice it touched as a
//! [`WindowFixture`](ferro_hgvs::conformance::reference_window::WindowFixture).
//!
//! Mirrors `extract_biocommons_windows.rs` — same recorder, same padded window
//! merging, same "committed and reviewed, never regenerated in CI" contract.
//! `--check` is a *local* guard: run against the same manifest it fails when the
//! committed fixture no longer matches a fresh extraction. The per-PR gate
//! consumes the committed file and needs no manifest.
//!
//! Run:   `cargo run --features dev --example extract_split_member_separation_windows -- --manifest PATH`
//! Check: `cargo run --features dev --example extract_split_member_separation_windows -- --check`

use std::any::Any;
use std::collections::BTreeSet;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::Arc;

use clap::Parser;
use ferro_hgvs::conformance::reference_window::{GenomicWindow, WindowFixture};
use ferro_hgvs::spdi::convert::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MultiFastaProvider, Normalizer};

#[path = "common/recording.rs"]
mod recording;

use recording::{build_window_fixture, RecordingProvider};

const INPUTS_TXT: &str = "tests/fixtures/split-member-separation/inputs.txt";
const WINDOWS_JSON: &str = "tests/fixtures/split-member-separation/reference-windows.json";
/// Named in **every** way `--check` can fail, because each of them is fixed by
/// running it and a reader should never have to reconstruct the command.
const REGENERATE: &str =
    "cargo run --features dev --example extract_split_member_separation_windows";

#[derive(Parser, Debug)]
#[command(about = "Extract hermetic reference windows for the issue #1539 separation guard")]
struct Cli {
    /// Prepared-reference manifest. Falls back to `$FERRO_MANIFEST` then
    /// `benchmark-output/manifest.json`.
    #[arg(long)]
    manifest: Option<PathBuf>,
    /// Compare the committed fixture to a fresh extraction; exit non-zero on
    /// drift and write nothing. Requires the manifest.
    #[arg(long)]
    check: bool,
}

/// Parse the committed input list: one HGVS string per line, `#` comments and
/// blank lines ignored. The guard applies the identical rule — keep the two in
/// step (`tests/it/issue_1539_split_member_separation.rs::inputs`).
fn read_inputs(path: &Path) -> std::io::Result<Vec<String>> {
    Ok(std::fs::read_to_string(path)?
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_string)
        .collect())
}

/// Flatten a normalized result into the members the guard adjudicates.
fn members_of(variant: &HgvsVariant) -> Vec<&HgvsVariant> {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.iter().flat_map(members_of).collect(),
        other => vec![other],
    }
}

/// Recover a printable message from a `catch_unwind` payload.
fn panic_message(payload: Box<dyn Any + Send>) -> String {
    if let Some(message) = payload.downcast_ref::<&'static str>() {
        return (*message).to_string();
    }
    if let Some(message) = payload.downcast_ref::<String>() {
        return message.clone();
    }
    "panic with a non-string payload".to_string()
}

/// Run one input through the recording provider, ignoring every *`Err`* outcome
/// — only the reference accesses matter, and they are recorded as they happen.
/// The member walk is what makes the capture sufficient for the guard:
/// adjudicating a member re-reads its reference bases through `hgvs_to_spdi`.
///
/// A **panic** is a different matter and is returned rather than swallowed. A
/// pass that unwinds partway through has made only some of its reads, so the
/// recording is short by exactly the reads it never reached — and a fixture
/// written from it under-serves the very pass the guard replays. That is not a
/// hypothetical: `tests/it/issue_1539_split_member_separation.rs` documents the
/// measured false negative an under-serving slice produces, where normalization
/// bails, returns its input unchanged, and scores as conformant. Catching per
/// input (rather than letting the first panic abort the process) is deliberate,
/// so one run reports *every* bad input; [`main`] then refuses to write.
fn run_input(provider: &RecordingProvider, input: &str) -> Result<(), String> {
    catch_unwind(AssertUnwindSafe(|| {
        let normalizer = Normalizer::new(provider.clone());
        let Ok(parsed) = parse_hgvs(input) else {
            return;
        };
        let Ok(normalized) = normalizer.normalize(&parsed) else {
            return;
        };
        for member in members_of(&normalized) {
            let _ = hgvs_to_spdi(member, provider);
        }
    }))
    .map_err(panic_message)
}

/// Widen every window whose accession is also a stored transcript to cover that
/// transcript whole.
///
/// This is a fidelity fix, not a size tweak, and it is load-bearing. The
/// normalizer reaches a transcript through `get_genomic_sequence` as well as
/// `get_transcript`, so the recorder files those reads under `contig_ranges` and
/// the accession lands in the fixture's `genomic` map. That makes
/// `WindowProvider::is_known_contig` true for it, which routes **`get_sequence`
/// for the transcript** through the genomic-window path instead of the whole
/// stored sequence — so every read outside the captured window errors, even
/// though the bases are sitting right there in `transcripts`.
///
/// A window captured from a pass is exactly as wide as that pass happened to
/// read, which makes the fixture self-fulfilling: it serves the recorded pass
/// and fails any pass that reads further. Measured against the manifest with
/// the padded windows in place, **every** 200 bp block of all five transcripts
/// failed under `WindowProvider` while succeeding under `MultiFastaProvider`.
/// A code change that makes normalization read more — precisely the change this
/// fixture exists to catch — therefore silently bails and returns the input
/// unchanged, which scores as conformant.
///
/// Serving the accession whole removes the failure mode by construction: there
/// is no "outside the window" left. `tests/it/issue_1539_split_member_separation.rs`
/// asserts it hermetically, so a regenerated fixture that lost the widening
/// fails rather than going quietly blind.
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

/// The accessions the corpus requires: the sequence each input names before its
/// `:`. Derived from `inputs.txt` rather than pinned to a count, so adding a row
/// raises the bar automatically instead of leaving a guessed floor behind.
fn required_accessions(inputs: &[String]) -> BTreeSet<String> {
    inputs
        .iter()
        .filter_map(|input| input.split(':').next())
        .filter(|accession| !accession.is_empty())
        .map(str::to_string)
        .collect()
}

/// Required accessions the fixture does not carry as a stored transcript.
fn missing_transcripts(fixture: &WindowFixture, required: &BTreeSet<String>) -> Vec<String> {
    let captured: BTreeSet<&str> = fixture
        .transcripts
        .iter()
        .map(|transcript| transcript.id.as_str())
        .collect();
    required
        .iter()
        .filter(|accession| !captured.contains(accession.as_str()))
        .cloned()
        .collect()
}

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

/// Read the manifest's `prepared_at` field for provenance, if present.
fn read_prepared_at(manifest_path: &Path) -> Option<String> {
    let content = std::fs::read_to_string(manifest_path).ok()?;
    let value: serde_json::Value = serde_json::from_str(&content).ok()?;
    value
        .get("prepared_at")
        .and_then(|v| v.as_str())
        .map(str::to_string)
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

    let Some(manifest_path) = resolve_manifest(cli.manifest) else {
        eprintln!(
            "extract_split_member_separation_windows: no manifest (pass --manifest, set \
             FERRO_MANIFEST, or place it at benchmark-output/manifest.json). This generator \
             requires the full prepared reference."
        );
        return ExitCode::FAILURE;
    };

    let inner = match MultiFastaProvider::from_manifest(&manifest_path) {
        Ok(p) => Arc::new(p),
        Err(e) => {
            eprintln!("from_manifest({}) failed: {e}", manifest_path.display());
            return ExitCode::FAILURE;
        }
    };

    let inputs = match read_inputs(Path::new(INPUTS_TXT)) {
        Ok(i) if !i.is_empty() => i,
        Ok(_) => {
            eprintln!("{INPUTS_TXT} holds no inputs; nothing to capture");
            return ExitCode::FAILURE;
        }
        Err(e) => {
            eprintln!("read {INPUTS_TXT}: {e}");
            return ExitCode::FAILURE;
        }
    };

    let recording = RecordingProvider::new(inner);
    let mut panicked = Vec::new();
    for input in &inputs {
        if let Err(message) = run_input(&recording, input) {
            panicked.push(format!("{input}: {message}"));
        }
    }
    // Report every panicking input, then refuse. A pass that unwound made only
    // part of its reads, so the recording is short by the rest — and a fixture
    // written from a short recording is exactly the under-serving slice this
    // generator's `serve_transcripts_whole` and the guard's
    // zero-failed-reads assertion exist to prevent.
    if !panicked.is_empty() {
        eprintln!(
            "{} of {} inputs panicked during capture, so the recording is missing every read \
             they had not yet made. Refusing to write a fixture built from a partial pass:\n  {}",
            panicked.len(),
            inputs.len(),
            panicked.join("\n  "),
        );
        return ExitCode::FAILURE;
    }

    let captured_from = read_prepared_at(&manifest_path).unwrap_or_default();
    let description = format!(
        "GENERATED by `cargo run --features dev --example \
         extract_split_member_separation_windows`. Hermetic reference windows for the issue \
         #1539 separation guard: every transcript and padded genomic window a normalize pass \
         over {INPUTS_TXT} — plus the per-member SPDI resolution the guard performs — touches. \
         Do not hand-edit; regenerate from the manifest."
    );
    // Transcripts are the whole point here (the inputs are all `c.` on `NM_`),
    // so keep every one that was resolved.
    let mut fixture = build_window_fixture(&recording, captured_from, description, |_| true);
    serve_transcripts_whole(&mut fixture);

    // The capture must carry every accession the corpus names, and this is an
    // assertion rather than the printed count it replaces: a printed count is
    // read by a human who already believes the run worked, so a fixture that
    // lost records passes `--check` quietly — which would certify precisely the
    // truncated slice this file exists to rule out. Coverage is checked instead
    // of a count because it cannot drift as rows are added, and it names what is
    // missing rather than reporting an arithmetic surprise.
    let required = required_accessions(&inputs);
    let missing = missing_transcripts(&fixture, &required);
    if !missing.is_empty() {
        eprintln!(
            "captured {} of the {} accession(s) {INPUTS_TXT} requires; missing {missing:?}. The \
             pass never resolved them, so the fixture would under-serve the guard — check the \
             manifest carries these transcripts, then re-run.",
            required.len() - missing.len(),
            required.len(),
        );
        return ExitCode::FAILURE;
    }

    let rendered = match fixture.to_json() {
        Ok(s) => s,
        Err(e) => {
            eprintln!("serialize window fixture: {e}");
            return ExitCode::FAILURE;
        }
    };

    println!(
        "extract_split_member_separation_windows: {} inputs, {} transcripts (covering all {} \
         required accessions), {} genomic windows ({} bp captured)",
        inputs.len(),
        fixture.transcripts.len(),
        required.len(),
        fixture.genomic.len(),
        fixture.genomic.iter().map(|w| w.bases.len()).sum::<usize>(),
    );

    if cli.check {
        // Every way of failing to *load* the committed fixture is one command
        // away from fixed, so every one of them names that command. Reading it
        // with `unwrap_or_default()` instead turned an absent file into an empty
        // string, which reached `serde_json` and surfaced as the bare
        // `EOF while parsing a value at line 1 column 0` — accurate, hintless,
        // and describing a parse problem in a file that was not there at all.
        let committed = match std::fs::read_to_string(WINDOWS_JSON) {
            Ok(content) => content,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                eprintln!(
                    "{WINDOWS_JSON} does not exist, so there is nothing to check against; \
                     generate it with `{REGENERATE}`"
                );
                return ExitCode::FAILURE;
            }
            Err(e) => {
                eprintln!("read {WINDOWS_JSON}: {e}; regenerate with `{REGENERATE}`");
                return ExitCode::FAILURE;
            }
        };
        // Coverage before equality: a truncated committed fixture is the failure
        // worth naming, and "out of date" would describe it as drift.
        match serde_json::from_str::<WindowFixture>(&committed) {
            Ok(on_disk) => {
                let lost = missing_transcripts(&on_disk, &required);
                if !lost.is_empty() {
                    eprintln!(
                        "{WINDOWS_JSON} is missing the transcript for {lost:?}, which \
                         {INPUTS_TXT} requires; regenerate with `{REGENERATE}`"
                    );
                    return ExitCode::FAILURE;
                }
            }
            // A corrupt or hand-edited fixture is the same class of problem as an
            // absent one — the file cannot be loaded, and regenerating is the
            // remedy — so it gets the same hint rather than a raw serde error.
            Err(e) => {
                eprintln!(
                    "{WINDOWS_JSON} is not readable as a window fixture ({e}); it is corrupt or \
                     hand-edited. Regenerate with `{REGENERATE}`"
                );
                return ExitCode::FAILURE;
            }
        }
        if committed == rendered {
            println!("{WINDOWS_JSON} is up to date");
            return ExitCode::SUCCESS;
        }
        eprintln!("{WINDOWS_JSON} is out of date; regenerate with `{REGENERATE}`");
        return ExitCode::FAILURE;
    }

    if let Err(e) = std::fs::write(WINDOWS_JSON, &rendered) {
        eprintln!("write {WINDOWS_JSON}: {e}");
        return ExitCode::FAILURE;
    }
    println!("wrote {WINDOWS_JSON}");
    ExitCode::SUCCESS
}
