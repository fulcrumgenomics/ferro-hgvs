//! Measure what the `FERRO_PARTITION` default flip moves over the **real**
//! corpora — ClinVar, CMRG and Paraphase — one arm per process, then diff.
//!
//! # Why this exists
//!
//! `rulings[delins-recommendation-reach-when-the-input-arrives-split]` says the
//! real-corpus disclosure is owed by "the change that makes that arm the
//! default". That change shipped with a **designed**-corpus measurement (cis
//! confluence classes) and with the manifest-backed conformance axes
//! dispositioned, but with no row count over stored expressions. This is the
//! instrument for that count, and the numbers it produced are in the PR that
//! added it.
//!
//! The sibling instrument is `measure_spec_conformance_per_arm`, which asks a
//! different question — agreement with the **spec's** stated form, per arm, over
//! the 934-row spec corpus. This one asks the consumer's question: *does the
//! string I have stored change?* Neither substitutes for the other, and the two
//! must not be added together.
//!
//! # What is measured, stated so nobody reads it as something else
//!
//! Two arms, two processes, one row set. For each row the pass records the
//! normalized string and its **SPDI key** (`equivalence::spdi_key`), and the
//! comparison then splits every difference in two:
//!
//! * **moved** — the string changed and both arms key to the *same* SPDI, so the
//!   description was respelled and denotes the same bases. This is the
//!   re-bucketing cost a consumer keying on the exact string pays.
//! * **denotes-different-bases** — the string changed and the two keys differ,
//!   or one arm has no key at all. Materially larger: the two descriptions are
//!   not two spellings of one variant.
//!
//! A row whose *status* changes (one arm normalizes, the other refuses) is
//! counted separately as `status-changed` rather than folded into either, since
//! it is neither a respelling nor a redenotation.
//!
//! # The effective denominator, and why a zero here would be a claim
//!
//! A difference count is only evidence if the rows reached the switch at all.
//! `normalize::partition_blocks_cut()` counts blocks cut under **every** arm
//! (`partition_decline_counts` cannot — it is structurally zero on `live`), and
//! this pass reads it either side of each row when it runs single-threaded. So
//! `rows_reaching_partitioner` is measured rather than predicted, and a run
//! reporting zero there has measured nothing about the arms whatever its diff
//! says.
//!
//! The counter is `#[cfg(debug_assertions)]`. Build with
//! `CARGO_PROFILE_DEV_OPT_LEVEL=3` to keep it while staying fast; a release run
//! reports `blocks_cut_available: false` and its denominator must not be quoted.
//!
//! # Bounded cost, and what is excluded
//!
//! #1846 records that `insertion_to_duplication` is Theta(payload^2), so a
//! cross-reference insert naming a megabase payload — ClinVar really does carry
//! `NC_000004.12:g.134850793_134850794ins[NC_000023.11:g.89555676_100352080]` —
//! does not terminate in any time worth waiting for. That is what stalled the
//! first attempt at this measurement.
//!
//! Rather than a timeout (which is not deterministic and cannot be reproduced
//! from the artifact), rows are excluded **before** normalizing, by a coarse
//! span bound read off the input string: see [`coordinate_span`]. The bound is
//! deliberately over-inclusive, every excluded row is written to the artifact
//! with status `excluded-oversize`, and the count is reported. An excluded row
//! is a row this measurement does not cover — never a row it found unchanged.
//!
//! # Sampling
//!
//! `--stride N` takes every Nth row of each corpus starting at `--offset`, which
//! is deterministic, reproducible from the artifact's own header, and spread
//! across each corpus's ordering rather than truncated to its head. `--stride 1`
//! is the full corpus. A run that hits `--time-budget-seconds` **removes its own
//! artifact** rather than leaving a short file that reads like a census, and the
//! comparison refuses any file holding fewer rows than its header plans.
//!
//! # Running it
//!
//! ```text
//! # arm A — the shipped default
//! FERRO_MANIFEST=<manifest> cargo run --features dev \
//!   --example measure_flip_corpus_impact -- --stride 100 --output a.tsv
//!
//! # arm B — the behaviour from before the flip
//! FERRO_PARTITION=live FERRO_MANIFEST=<manifest> cargo run --features dev \
//!   --example measure_flip_corpus_impact -- --stride 100 --output b.tsv
//!
//! # the diff
//! cargo run --features dev --example measure_flip_corpus_impact -- \
//!   --compare-a b.tsv --compare-b a.tsv --summary flip.json
//! ```

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};

use ferro_hgvs::conformance::completeness::CaptureLedger;
use ferro_hgvs::conformance::error_mode_stamp::ErrorModeStamp;
use ferro_hgvs::equivalence::spdi_key;
use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

/// The four corpora, in the order the disclosure quotes them.
///
/// `clinvar_hgvs_unique` is the broadest and is **not** part of the 5,761,302
/// total the v0.13.0/v0.13.1 disclosures used; it is listed here so a run can
/// include it deliberately rather than by accident.
const DEFAULT_CORPORA: [&str; 4] = [
    "tests/fixtures/bulk/clinvar_hgvs_500k.json.gz",
    "tests/fixtures/bulk/clinvar_hgvs_unique.json.gz",
    "tests/fixtures/validation/cmrg_genes_exhaustive.json.gz",
    "tests/fixtures/validation/paraphase_genes_exhaustive.json.gz",
];

/// Column header of the per-row artifact. Read back by the comparison, so the
/// two halves cannot drift.
const TSV_HEADER: &str = "corpus\tindex\tinput\tstatus\tnormalized\tspdi\tblocks_cut";

#[derive(Parser, Debug)]
#[command(
    about = "Measure the FERRO_PARTITION flip's effect on stored ClinVar/CMRG/Paraphase strings"
)]
struct Cli {
    /// A corpus to read. Repeatable. Defaults to all four bulk fixtures.
    #[arg(long)]
    corpus: Vec<PathBuf>,

    /// Prepared-reference manifest. Falls back to `$FERRO_MANIFEST`.
    #[arg(long)]
    manifest: Option<PathBuf>,

    /// Where to write the per-row artifact.
    #[arg(long)]
    output: Option<PathBuf>,

    /// Take every Nth row of each corpus. 1 is the full corpus.
    #[arg(long, default_value_t = 1)]
    stride: usize,

    /// First row of each corpus to take, before striding.
    #[arg(long, default_value_t = 0)]
    offset: usize,

    /// Widest coordinate span, in bases, this pass will normalize. Rows above it
    /// are excluded and reported; see [`coordinate_span`].
    #[arg(long, default_value_t = 1_000_000)]
    max_span: u64,

    /// Longest input, in bytes, this pass will normalize. The second half of the
    /// cost bound: a payload written out literally (`ins<200 kb of bases>`)
    /// carries no coordinates for [`coordinate_span`] to read, and it is the
    /// payload length #1846's Theta(payload^2) is quadratic in.
    #[arg(long, default_value_t = 10_000)]
    max_input_bytes: usize,

    /// Stop after this many wall seconds and mark the run partial. 0 disables.
    #[arg(long, default_value_t = 0)]
    time_budget_seconds: u64,

    /// Print a progress line every N rows. 0 disables.
    #[arg(long, default_value_t = 20_000)]
    progress_every: usize,

    /// Compare two artifacts instead of measuring. Both halves are required.
    #[arg(long, requires = "compare_b")]
    compare_a: Option<PathBuf>,

    /// The second artifact of a comparison — conventionally the shipped arm.
    #[arg(long, requires = "compare_a")]
    compare_b: Option<PathBuf>,

    /// Where a comparison writes its summary.
    #[arg(long)]
    summary: Option<PathBuf>,

    /// How many differing rows a comparison writes into its summary. The count
    /// itself is never capped — `differences_recorded` reports every one.
    #[arg(long, default_value_t = 5_000)]
    max_differences: usize,
}

/// The subset of a bulk-fixture row this pass reads. All four corpora share it.
#[derive(Debug, Deserialize)]
struct CorpusCase {
    input: String,
}

#[derive(Debug, Deserialize)]
struct CorpusFixture {
    test_cases: Vec<CorpusCase>,
}

/// Blocks this process has cut, or `None` in a release build, where the counter
/// does not exist.
#[cfg(debug_assertions)]
fn blocks_cut() -> Option<u64> {
    Some(ferro_hgvs::normalize::partition_blocks_cut())
}

#[cfg(not(debug_assertions))]
fn blocks_cut() -> Option<u64> {
    None
}

/// The widest span, in bases, named by any one coordinate region of `input`.
///
/// A coordinate region opens at an axis marker (`g.`, `c.`, `n.`, `r.`, `m.`,
/// `o.`, `p.`) and runs to the next one, so the integers of an accession
/// (`NC_000023.11`) never join those of a position, and a cross-reference
/// payload (`ins[NC_000023.11:g.89555676_100352080]`) is measured as its own
/// region rather than against the host's coordinates.
///
/// This is a **coarse** bound and is meant to be. It over-estimates a row whose
/// region mixes an interval with an offset, and it cannot see a payload whose
/// length is implied rather than written. What it buys is determinism: the same
/// input is excluded on every run, on every machine, and the exclusion is
/// legible from the input string alone — which a wall-clock timeout is not.
///
/// Brackets are read three ways, because three different things are written in
/// them and only two hold coordinates:
///
/// - **A group holding an axis marker** is a cross-reference payload
///   (`ins[NC_000023.11:g.89555676_100352080]`). The host region is closed at the
///   `[`, since the payload's own accession follows and its digits are not
///   coordinates, and the inner axis marker opens the payload's region.
/// - **A group opening immediately after the axis marker** is an allele or
///   compound (`g.[36_37insC;40del]`), whose members are numbered on the host
///   axis, so it is descended into with the region left open.
/// - **Anything else** is a repeat count or a literal payload list and is skipped
///   whole: `g.39665402TG[1]` writes a *count*, and reading that `1` as a
///   coordinate makes a two-base repeat look 39 Mb wide.
///
/// A `:` also closes the region, so an accession written after one (inside a
/// payload, or after a compound's parenthesised parent) cannot join it.
///
/// `None` means "unbounded": a `pter`/`qter`/`cen` endpoint names no coordinate,
/// so no numeric bound can be read off it and the row is excluded.
fn coordinate_span(input: &str) -> Option<u64> {
    for unbounded in ["pter", "qter", "cen"] {
        if input.contains(unbounded) {
            return None;
        }
    }
    let bytes = input.as_bytes();
    let mut widest = 0u64;
    // `in_region` and `bounds` are separate because a region can legitimately
    // contain no digits (`c.` followed immediately by another axis marker), and
    // folding an empty region's sentinel bounds underflows. That underflow
    // panicked a debug build 3,982,764 rows into the first full pass.
    let mut in_region = false;
    let mut bounds: Option<(u64, u64)> = None;
    fn fold(bounds: &mut Option<(u64, u64)>, widest: &mut u64) {
        if let Some((lo, hi)) = bounds.take() {
            *widest = (*widest).max(hi - lo);
        }
    }
    // Whether the byte just consumed was the axis marker, which is what tells an
    // allele group (`g.[…]`) from a repeat count (`…TG[1]`).
    let mut at_axis_marker = false;
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'[' {
            let close = matching_bracket(bytes, i);
            if holds_axis_marker(&bytes[i + 1..close.min(bytes.len())]) {
                fold(&mut bounds, &mut widest);
                in_region = false;
            } else if !at_axis_marker {
                i = close.saturating_add(1);
                at_axis_marker = false;
                continue;
            }
            at_axis_marker = false;
            i += 1;
            continue;
        }
        // A `:` separates an accession from its axis, so nothing before it on
        // this side belongs to the region in progress.
        if bytes[i] == b':' {
            fold(&mut bounds, &mut widest);
            in_region = false;
            at_axis_marker = false;
            i += 1;
            continue;
        }
        // An axis marker closes the region in progress and opens a new one.
        if i + 1 < bytes.len()
            && bytes[i + 1] == b'.'
            && matches!(bytes[i], b'g' | b'c' | b'n' | b'r' | b'm' | b'o' | b'p')
        {
            fold(&mut bounds, &mut widest);
            in_region = true;
            at_axis_marker = true;
            i += 2;
            continue;
        }
        at_axis_marker = false;
        if bytes[i].is_ascii_digit() {
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            // A number too long to be a coordinate is not one; treating it as
            // `u64::MAX` would exclude the row for the wrong reason.
            if in_region {
                if let Ok(value) = input[start..i].parse::<u64>() {
                    bounds = Some(match bounds {
                        None => (value, value),
                        Some((lo, hi)) => (lo.min(value), hi.max(value)),
                    });
                }
            }
            continue;
        }
        i += 1;
    }
    fold(&mut bounds, &mut widest);
    Some(widest)
}

/// The index of the `]` closing the `[` at `open`, or `bytes.len()` when the
/// string is unbalanced — an unbalanced group is treated as running to the end,
/// which is the conservative reading.
fn matching_bracket(bytes: &[u8], open: usize) -> usize {
    let mut depth = 0usize;
    for (offset, byte) in bytes[open..].iter().enumerate() {
        match byte {
            b'[' => depth += 1,
            b']' => {
                depth -= 1;
                if depth == 0 {
                    return open + offset;
                }
            }
            _ => {}
        }
    }
    bytes.len()
}

/// Whether a slice holds an axis marker (`g.`, `c.`, …), i.e. whether its digits
/// are coordinates. An accession's version dot never matches: the character
/// before it is a digit, not an axis letter.
fn holds_axis_marker(bytes: &[u8]) -> bool {
    bytes.windows(2).any(|pair| {
        pair[1] == b'.' && matches!(pair[0], b'g' | b'c' | b'n' | b'r' | b'm' | b'o' | b'p')
    })
}

/// One row of the per-row artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Row {
    corpus: String,
    index: usize,
    input: String,
    /// `ok`, `parse-error`, `normalize-error` or `excluded-oversize`.
    status: String,
    /// The normalized rendering, or empty on any non-`ok` status.
    normalized: String,
    /// The SPDI key of the normalized form, or empty when it has none.
    spdi: String,
    /// Blocks this row cut, or empty when the counter is unavailable.
    blocks_cut: String,
}

impl Row {
    fn to_tsv(&self) -> String {
        // No field can contain a tab or a newline: an HGVS string cannot, and
        // the statuses and keys are generated. Asserted rather than assumed.
        debug_assert!(!self.input.contains(['\t', '\n']));
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.corpus,
            self.index,
            self.input,
            self.status,
            self.normalized,
            self.spdi,
            self.blocks_cut
        )
    }

    fn from_tsv(line: &str) -> anyhow::Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 7 {
            anyhow::bail!("expected 7 tab-separated fields, found {}", fields.len());
        }
        Ok(Row {
            corpus: fields[0].to_string(),
            index: fields[1].parse()?,
            input: fields[2].to_string(),
            status: fields[3].to_string(),
            normalized: fields[4].to_string(),
            spdi: fields[5].to_string(),
            blocks_cut: fields[6].to_string(),
        })
    }
}

/// The artifact header, written as a `#`-prefixed JSON line so the file is both
/// self-describing and still a TSV a consumer can `grep -v '^#'`.
#[derive(Debug, Serialize, Deserialize)]
struct Header {
    /// The raw `FERRO_PARTITION` value this process ran under. `canonical-coalesced`
    /// when unset, which is the shipped default since the flip.
    arm: String,
    /// Whether `partition_blocks_cut` existed in this build.
    blocks_cut_available: bool,
    /// Manifest file name only — an absolute path would leak the local layout.
    manifest: String,
    /// The error-handling precondition these rows were normalized under. A
    /// census compared against one from another mode is as wrong as one built
    /// from a partial pass, and nothing else in the file would say so.
    normalized_under: ErrorModeStamp,
    stride: usize,
    offset: usize,
    max_span: u64,
    max_input_bytes: usize,
    corpora: Vec<String>,
    /// Rows each corpus holds, before striding.
    corpus_totals: BTreeMap<String, usize>,
    /// Rows this run set out to measure. The comparison refuses a file holding
    /// fewer, so a truncated or interrupted artifact cannot be read as a result.
    selected: usize,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    // The positive control the whole measurement rests on: a misspelled
    // `FERRO_PARTITION` would otherwise file one arm's numbers under another's
    // name.
    if let Some(message) = ferro_hgvs::normalize::partition_switch_startup_error() {
        eprintln!("error: {message}");
        return ExitCode::FAILURE;
    }
    let outcome = match (&cli.compare_a, &cli.compare_b) {
        (Some(a), Some(b)) => compare(a, b, cli.summary.as_deref(), cli.max_differences),
        _ => measure(&cli),
    };
    match outcome {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// Decompress a corpus and hand back the rows this run selects.
fn select_rows(path: &Path, stride: usize, offset: usize) -> anyhow::Result<(usize, Vec<String>)> {
    let file = File::open(path).map_err(|e| {
        anyhow::anyhow!(
            "open {}: {e}. Fetch the bulk fixtures first: `scripts/fetch-test-fixtures.sh`",
            path.display()
        )
    })?;
    let mut buf = Vec::new();
    GzDecoder::new(file)
        .read_to_end(&mut buf)
        .map_err(|e| anyhow::anyhow!("decompress {}: {e}", path.display()))?;
    let fixture: CorpusFixture = serde_json::from_slice(&buf)
        .map_err(|e| anyhow::anyhow!("parse {}: {e}", path.display()))?;
    let total = fixture.test_cases.len();
    let selected = fixture
        .test_cases
        .into_iter()
        .enumerate()
        .filter(|(i, _)| *i >= offset && (i - offset).is_multiple_of(stride))
        .map(|(_, case)| case.input)
        .collect();
    Ok((total, selected))
}

fn measure(cli: &Cli) -> anyhow::Result<()> {
    if cli.stride == 0 {
        anyhow::bail!("--stride 0 selects nothing");
    }
    let output = cli
        .output
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--output is required when measuring"))?;
    let manifest = cli
        .manifest
        .clone()
        .or_else(|| std::env::var("FERRO_MANIFEST").ok().map(PathBuf::from))
        .ok_or_else(|| {
            anyhow::anyhow!(
                "no manifest (pass --manifest or set FERRO_MANIFEST). Without a real reference \
                 nothing reaches the partitioner and this pass measures nothing."
            )
        })?;
    let corpora: Vec<PathBuf> = if cli.corpus.is_empty() {
        DEFAULT_CORPORA.iter().map(PathBuf::from).collect()
    } else {
        cli.corpus.clone()
    };

    let provider = MultiFastaProvider::from_manifest(&manifest)
        .map_err(|e| anyhow::anyhow!("load manifest {}: {e}", manifest.display()))?;
    // `NormalizeConfig::default()` is `ErrorConfig::lenient()`, which is what the
    // sibling per-arm measurement uses. Kept identical so the two censuses are
    // comparable; a census taken under another mode is not.
    let config = NormalizeConfig::default();
    let normalized_under = ErrorModeStamp::of(&config.error_config);
    let normalizer = Normalizer::with_config(provider, config);

    let arm = std::env::var("FERRO_PARTITION")
        .ok()
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| "canonical-coalesced".to_string());

    let mut corpus_totals = BTreeMap::new();
    let mut work: Vec<(String, Vec<String>)> = Vec::new();
    for path in &corpora {
        let name = corpus_name(path);
        let (total, rows) = select_rows(path, cli.stride, cli.offset)?;
        corpus_totals.insert(name.clone(), total);
        work.push((name, rows));
    }
    let selected: usize = work.iter().map(|(_, rows)| rows.len()).sum();
    if selected == 0 {
        anyhow::bail!("the selection is empty; there is nothing to measure");
    }
    eprintln!("selecting {selected} rows (stride {})", cli.stride);

    let header = Header {
        arm,
        blocks_cut_available: blocks_cut().is_some(),
        normalized_under,
        manifest: manifest
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default(),
        stride: cli.stride,
        offset: cli.offset,
        max_span: cli.max_span,
        max_input_bytes: cli.max_input_bytes,
        corpora: corpora.iter().map(|p| corpus_name(p)).collect(),
        corpus_totals,
        selected,
    };

    // Streamed rather than accumulated: at `--stride 1` the four corpora hold
    // 9,949,738 rows, and holding every rendered row in memory alongside the
    // provider is gigabytes for no reason. The consequence is that the header —
    // which carries the planned row count — is written before the run, and a run
    // that does not finish removes its own file rather than leaving a short one
    // that reads like a census.
    let started = std::time::Instant::now();
    let file =
        File::create(output).map_err(|e| anyhow::anyhow!("create {}: {e}", output.display()))?;
    let mut out = BufWriter::new(file);
    writeln!(out, "# {}", serde_json::to_string(&header)?)?;
    writeln!(out, "{TSV_HEADER}")?;

    // One row in, one row out. A refusal is a *result* here — `parse-error` and
    // `excluded-oversize` are buckets — so this pass has no legitimate drop, and
    // the ledger is what makes that visible rather than assumed.
    let mut ledger = CaptureLedger::new("corpus rows");
    let mut measured = 0usize;
    let mut stopped_early = false;
    'corpora: for (name, inputs) in &work {
        for (index, input) in inputs.iter().enumerate() {
            if cli.time_budget_seconds > 0 && started.elapsed().as_secs() >= cli.time_budget_seconds
            {
                stopped_early = true;
                break 'corpora;
            }
            if cli.progress_every > 0 && measured.is_multiple_of(cli.progress_every) {
                eprintln!(
                    "  [{measured}/{selected}] {:>5}s  {input}",
                    started.elapsed().as_secs()
                );
            }
            measured += 1;
            let row = measure_one(
                &normalizer,
                name,
                index,
                input,
                cli.max_span,
                cli.max_input_bytes,
            );
            writeln!(out, "{}", row.to_tsv())?;
            ledger.record_success();
        }
    }
    if stopped_early {
        drop(out);
        // A partial pass is not a result, and an artifact that could be mistaken
        // for one is worse than none. #1551 is the precedent.
        let _ = std::fs::remove_file(output);
        anyhow::bail!(
            "--time-budget-seconds reached after {measured} of {selected} rows; the partial \
             artifact has been removed. Re-run with a wider budget, or a wider --stride."
        );
    }
    // Stamped into the artifact as a trailing `#` line, since the header is
    // written before the pass and cannot carry a count the pass produces. A
    // shortfall removes the file rather than reporting itself in it.
    let counts = match ledger.finish() {
        Ok(counts) => counts,
        Err(shortfall) => {
            drop(out);
            let _ = std::fs::remove_file(output);
            anyhow::bail!("{shortfall}");
        }
    };
    writeln!(out, "# {}", serde_json::to_string(&counts)?)?;
    out.flush()?;
    drop(out);
    eprintln!(
        "wrote {measured} rows to {} in {}s",
        output.display(),
        started.elapsed().as_secs()
    );
    Ok(())
}

/// The corpus's file name, which is what the artifact and the disclosure name.
fn corpus_name(path: &Path) -> String {
    path.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string())
}

/// Normalize one input and record everything the comparison needs.
fn measure_one(
    normalizer: &Normalizer<MultiFastaProvider>,
    corpus: &str,
    index: usize,
    input: &str,
    max_span: u64,
    max_input_bytes: usize,
) -> Row {
    let row = |status: &str, normalized: String, spdi: String, cut: String| Row {
        corpus: corpus.to_string(),
        index,
        input: input.to_string(),
        status: status.to_string(),
        normalized,
        spdi,
        blocks_cut: cut,
    };
    let within_bounds = input.len() <= max_input_bytes
        && coordinate_span(input).is_some_and(|span| span <= max_span);
    if !within_bounds {
        // Excluded, and recorded as excluded. Never as unchanged.
        return row(
            "excluded-oversize",
            String::new(),
            String::new(),
            String::new(),
        );
    }
    let variant = match parse_hgvs(input) {
        Ok(variant) => variant,
        // Parsing precedes the switch, so this bucket is arm-independent by
        // construction and the comparison expects it to match row for row.
        Err(_) => return row("parse-error", String::new(), String::new(), String::new()),
    };
    let before = blocks_cut();
    // Contained: a real reference opens paths a hermetic run never takes, and one
    // panic would cost the whole pass. Recorded as its own status rather than
    // swallowed.
    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        normalizer.normalize_with_diagnostics(&variant)
    }));
    let cut = match (before, blocks_cut()) {
        (Some(before), Some(after)) => (after - before).to_string(),
        _ => String::new(),
    };
    match outcome {
        Ok(Ok(diagnosed)) => {
            let key = spdi_key(&diagnosed.result, normalizer.provider())
                .map(|key| key.as_spdi().to_string())
                .unwrap_or_default();
            row("ok", format!("{}", diagnosed.result), key, cut)
        }
        Ok(Err(_)) => row("normalize-error", String::new(), String::new(), cut),
        Err(_) => row("panicked", String::new(), String::new(), cut),
    }
}

/// One arm's artifact, opened for a streaming read.
///
/// Streaming rather than slurped: at `--stride 1` each artifact is ~1.7 GB, and
/// holding both as `Vec<Row>` costs more memory than the measurement itself did.
struct Artifact {
    header: Header,
    lines: std::io::Lines<std::io::BufReader<File>>,
}

fn open_artifact(path: &Path) -> anyhow::Result<Artifact> {
    let file = File::open(path).map_err(|e| anyhow::anyhow!("open {}: {e}", path.display()))?;
    let mut lines = std::io::BufRead::lines(std::io::BufReader::new(file));
    let header_line = lines
        .next()
        .transpose()?
        .ok_or_else(|| anyhow::anyhow!("{} is empty", path.display()))?;
    let header: Header = serde_json::from_str(
        header_line
            .strip_prefix("# ")
            .ok_or_else(|| anyhow::anyhow!("{}: no `# {{...}}` header line", path.display()))?,
    )
    .map_err(|e| anyhow::anyhow!("{}: parse header: {e}", path.display()))?;
    let column_line = lines
        .next()
        .transpose()?
        .ok_or_else(|| anyhow::anyhow!("{}: no column header", path.display()))?;
    if column_line != TSV_HEADER {
        anyhow::bail!("{}: unexpected columns `{column_line}`", path.display());
    }
    Ok(Artifact { header, lines })
}

/// A difference between the arms, kept for the artifact so every headline number
/// has rows behind it.
#[derive(Debug, Serialize)]
struct Difference {
    corpus: String,
    input: String,
    before: String,
    after: String,
    before_spdi: String,
    after_spdi: String,
    class: &'static str,
}

#[derive(Debug, Serialize)]
struct Summary {
    description: &'static str,
    before_arm: String,
    after_arm: String,
    stride: usize,
    offset: usize,
    max_span: u64,
    max_input_bytes: usize,
    corpora: Vec<String>,
    corpus_totals: BTreeMap<String, usize>,
    /// Rows in each artifact — the denominator of everything below.
    rows_compared: usize,
    /// Rows both arms normalized. Only these can move.
    rows_normalized_on_both_arms: usize,
    /// Rows one arm refused, per status, keyed `before-status/after-status`.
    by_status_pair: BTreeMap<String, usize>,
    /// Rows excluded by the cost bounds, which this measurement does not cover.
    excluded_oversize: usize,
    /// Rows whose normalization panicked on either arm. Recorded, never dropped.
    panicked: usize,
    /// Whether the counter existed. `rows_reaching_partitioner` is meaningless
    /// without it.
    blocks_cut_available: bool,
    /// Rows that cut at least one block on the shipped arm, i.e. actually
    /// reached the switch.
    rows_reaching_partitioner: usize,
    /// Blocks cut over the whole shipped-arm run.
    blocks_cut: u64,
    /// Strings that changed and denote the same bases.
    moved: usize,
    /// Strings that changed and do not key to the same SPDI.
    denotes_different_bases: usize,
    /// Changed strings neither arm could key. Reported apart rather than
    /// assumed into either: an absent key is not evidence of agreement.
    moved_unkeyable: usize,
    /// Rows whose status differs between the arms.
    status_changed: usize,
    /// Differences found, whether or not each was kept below.
    differences_recorded: usize,
    /// The kept differences, so the headline numbers have rows behind them.
    /// Capped by `--max-differences`; `differences_recorded` is the true count.
    differences: Vec<Difference>,
}

fn compare(
    a: &Path,
    b: &Path,
    summary: Option<&Path>,
    max_differences: usize,
) -> anyhow::Result<()> {
    let summary_path =
        summary.ok_or_else(|| anyhow::anyhow!("--summary is required when comparing"))?;
    let mut before = open_artifact(a)?;
    let mut after = open_artifact(b)?;
    if before.header.arm == after.header.arm {
        anyhow::bail!(
            "both artifacts were taken under `{}`; a diff of one arm against itself \
             reports zero for a reason that has nothing to do with the flip",
            before.header.arm
        );
    }
    if before.header.stride != after.header.stride || before.header.offset != after.header.offset {
        anyhow::bail!("the two artifacts sample different rows and cannot be compared");
    }
    if before.header.selected != after.header.selected {
        anyhow::bail!(
            "the artifacts plan different row counts ({} vs {})",
            before.header.selected,
            after.header.selected
        );
    }

    let mut summary_out = Summary {
        description: "FERRO_PARTITION default flip, measured over stored real-corpus expressions",
        before_arm: before.header.arm.clone(),
        after_arm: after.header.arm.clone(),
        stride: after.header.stride,
        offset: after.header.offset,
        max_span: after.header.max_span,
        max_input_bytes: after.header.max_input_bytes,
        corpora: after.header.corpora.clone(),
        corpus_totals: after.header.corpus_totals.clone(),
        rows_compared: 0,
        rows_normalized_on_both_arms: 0,
        by_status_pair: BTreeMap::new(),
        excluded_oversize: 0,
        panicked: 0,
        blocks_cut_available: after.header.blocks_cut_available,
        rows_reaching_partitioner: 0,
        blocks_cut: 0,
        moved: 0,
        denotes_different_bases: 0,
        moved_unkeyable: 0,
        status_changed: 0,
        differences_recorded: 0,
        differences: Vec::new(),
    };

    loop {
        let (before_line, after_line) = match (before.lines.next(), after.lines.next()) {
            (None, None) => break,
            (Some(x), Some(y)) => (x?, y?),
            _ => anyhow::bail!("the artifacts hold different numbers of rows"),
        };
        // The trailing completeness stamp, not a row.
        if before_line.starts_with('#') && after_line.starts_with('#') {
            continue;
        }
        let before_row =
            Row::from_tsv(&before_line).map_err(|e| anyhow::anyhow!("{}: {e}", a.display()))?;
        let after_row =
            Row::from_tsv(&after_line).map_err(|e| anyhow::anyhow!("{}: {e}", b.display()))?;
        summary_out.rows_compared += 1;
        if before_row.input != after_row.input || before_row.corpus != after_row.corpus {
            anyhow::bail!(
                "row {} disagrees between the artifacts ({} vs {}); they are not aligned",
                after_row.index,
                before_row.input,
                after_row.input
            );
        }
        if let Ok(cut) = after_row.blocks_cut.parse::<u64>() {
            summary_out.blocks_cut += cut;
            if cut > 0 {
                summary_out.rows_reaching_partitioner += 1;
            }
        }
        if after_row.status == "panicked" || before_row.status == "panicked" {
            summary_out.panicked += 1;
        }
        if after_row.status == "excluded-oversize" {
            summary_out.excluded_oversize += 1;
            continue;
        }
        let record = |class: &'static str, summary_out: &mut Summary| {
            summary_out.differences_recorded += 1;
            if summary_out.differences.len() < max_differences {
                summary_out.differences.push(Difference {
                    corpus: after_row.corpus.clone(),
                    input: after_row.input.clone(),
                    before: before_row.normalized.clone(),
                    after: after_row.normalized.clone(),
                    before_spdi: before_row.spdi.clone(),
                    after_spdi: after_row.spdi.clone(),
                    class,
                });
            }
        };
        if before_row.status != after_row.status {
            summary_out.status_changed += 1;
            let pair = format!("{}/{}", before_row.status, after_row.status);
            record("status-changed", &mut summary_out);
            *summary_out.by_status_pair.entry(pair).or_default() += 1;
            continue;
        }
        if after_row.status != "ok" {
            continue;
        }
        summary_out.rows_normalized_on_both_arms += 1;
        if before_row.normalized == after_row.normalized {
            continue;
        }
        let class = if before_row.spdi.is_empty() || after_row.spdi.is_empty() {
            summary_out.moved_unkeyable += 1;
            "moved-unkeyable"
        } else if before_row.spdi == after_row.spdi {
            summary_out.moved += 1;
            "moved"
        } else {
            summary_out.denotes_different_bases += 1;
            "denotes-different-bases"
        };
        record(class, &mut summary_out);
    }

    if summary_out.rows_compared != after.header.selected {
        anyhow::bail!(
            "compared {} rows against a planned {}. A partial pass is not a result.",
            summary_out.rows_compared,
            after.header.selected
        );
    }

    let text = serde_json::to_string_pretty(&summary_out)?;
    std::fs::write(summary_path, format!("{text}\n"))
        .map_err(|e| anyhow::anyhow!("write {}: {e}", summary_path.display()))?;
    eprintln!(
        "{} rows compared: {} moved, {} denote different bases, {} unkeyable, {} status changes",
        summary_out.rows_compared,
        summary_out.moved,
        summary_out.denotes_different_bases,
        summary_out.moved_unkeyable,
        summary_out.status_changed
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_plain_genomic_substitution_spans_nothing() {
        // The accession's own digits must not join the coordinate's, or every
        // high-coordinate row on chromosome 1 would be excluded.
        assert_eq!(coordinate_span("NC_000001.11:g.150000000A>G"), Some(0));
    }

    #[test]
    fn an_interval_spans_its_own_width() {
        assert_eq!(
            coordinate_span("NM_000083.3:c.2461_2464delinsCTCC"),
            Some(3)
        );
    }

    #[test]
    fn a_cross_reference_payload_is_measured_in_its_own_region() {
        // The #1846 shape. The host coordinates span 1; the payload spans ~10.8 Mb,
        // and it is the payload that does not terminate.
        let span = coordinate_span(
            "NC_000004.12:g.134850793_134850794ins[NC_000023.11:g.89555676_100352080]",
        );
        assert_eq!(span, Some(10_796_404));
    }

    #[test]
    fn a_repeat_count_is_not_a_coordinate() {
        // The bracket holds a repeat *count*, not a position. Reading it as one
        // made a two-base repeat look 39 Mb wide and excluded it — measured on
        // the ClinVar 500k corpus, where two rows of a 1,001-row probe were lost
        // this way.
        assert_eq!(coordinate_span("NC_000017.11:g.39665402TG[1]"), Some(0));
    }

    #[test]
    fn an_allele_keeps_the_span_of_its_members() {
        // The region opens before the bracket, so the members' own coordinates
        // are still read even though the group carries no axis marker.
        assert_eq!(coordinate_span("NC_000001.11:g.[36_37insC;40del]"), Some(4));
    }

    #[test]
    fn an_unbounded_endpoint_has_no_span() {
        assert_eq!(
            coordinate_span("NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825272]"),
            None
        );
    }

    #[test]
    fn a_row_round_trips_through_the_artifact_format() {
        let row = Row {
            corpus: "clinvar_hgvs_500k.json.gz".to_string(),
            index: 7,
            input: "NM_000083.3:c.2461_2464delinsCTCC".to_string(),
            status: "ok".to_string(),
            normalized: "NM_000083.3:c.2461_2464delinsCTCC".to_string(),
            spdi: "NM_000083.3:2460:TTCC:CTCC".to_string(),
            blocks_cut: "1".to_string(),
        };
        assert_eq!(Row::from_tsv(&row.to_tsv()).unwrap(), row);
    }

    #[test]
    fn a_short_row_is_refused_rather_than_silently_truncated() {
        assert!(Row::from_tsv("a\tb\tc").is_err());
    }
}
