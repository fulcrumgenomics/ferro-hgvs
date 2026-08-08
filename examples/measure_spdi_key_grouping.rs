//! Measure whether [`spdi_key`] buckets the confluence divergences.
//!
//! # Why
//!
//! `src/equivalence/key.rs` claims that descriptions asserting different
//! *partitions* of one edit — which this project's normalizer deliberately keeps
//! in distinct canonical forms — nonetheless share one grouping key. The unit
//! tests pin that on hand-written pairs. This asks it of the population that
//! actually exhibits the divergence: every class
//! `examples/dump_confluence_divergences.rs` reports, where two spellings of one
//! designed variant reached two different normalized strings.
//!
//! The distinction matters because a hand-written pair proves the mechanism and
//! says nothing about coverage. A claim of the form "the key groups the
//! divergent classes" is a claim about a corpus, and this is the only thing that
//! can settle it.
//!
//! # What is measured, and how it can be vacuous
//!
//! Per class: the key of every **distinct normalized output**. Those outputs are
//! what a downstream consumer stores, and the class is in the input precisely
//! because there is more than one of them. Four verdicts:
//!
//! - `grouped` — every output keyed, and to the same key. The divergence is
//!   invisible to a consumer bucketing by key.
//! - `split` — every output keyed, to two or more keys. The key does **not**
//!   reconcile this class.
//! - `undecided` — at least one output has no key, so the class cannot be
//!   answered either way. Reported separately rather than folded into `split`:
//!   those are different findings and one must not be read as the other.
//! - `unmeasured` — at least one output could not be *asked*: it does not parse,
//!   or keying it panicked. Split out from `undecided` for the same reason
//!   `undecided` is split out from `split`. `undecided` is `spdi_key` declining,
//!   which is evidence about the key; `unmeasured` is this harness failing,
//!   which is evidence about the corpus. A run with any of these exits non-zero,
//!   because every count in the report is then over a smaller population than
//!   the input and a partial run would otherwise read as a clean one.
//!
//! **The vacuity check is printed unconditionally.** A run whose input holds no
//! class with two or more distinct outputs measures nothing about grouping, and
//! a `0 split` from such a run reads exactly like success. `multi_output` is
//! that denominator; when it is zero the report says so instead of reporting a
//! ratio.
//!
//! # Usage
//!
//! ```text
//! cargo run --release --features dev --example generate_cis_confluence_corpus
//! cargo run --release --features dev --example dump_confluence_divergences -- --out /tmp/div.json
//! cargo run --release --features dev --example measure_spdi_key_grouping -- --divergences /tmp/div.json
//! ```

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use serde::Deserialize;

use ferro_hgvs::equivalence::spdi_key;
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::spdi::canonical_spdi;

/// Kept in step with `examples/dump_confluence_divergences.rs` and
/// `examples/generate_cis_confluence_corpus.rs`; all three must agree on the
/// padding, the accessions and the CDS or the coordinates mean different things.
const PAD_OFFSET: usize = 256;
const GENOMIC_CONTIG: &str = "NC_TEST.1";
const TX_ACCESSION: &str = "NM_TEST.1";
const TX_CONTIG: &str = "chr_synth";
const CDS_START: u64 = 1;
const CDS_END: u64 = 63;

#[derive(Parser, Debug)]
#[command(about = "Measure how the SPDI grouping key buckets the divergent confluence classes")]
struct Cli {
    /// The JSON written by `dump_confluence_divergences -- --out <path>`.
    #[arg(long)]
    divergences: PathBuf,
    /// How many non-grouping exemplars to print per verdict.
    #[arg(long, default_value_t = 5)]
    examples: usize,
}

// ---------------------------------------------------------------------------
// Input — a subset of `dump_confluence_divergences`'s record
// ---------------------------------------------------------------------------

/// Only the fields the measurement needs; serde ignores the rest, so this does
/// not have to track that generator's record shape.
#[derive(Deserialize)]
struct Record {
    id: String,
    axis: String,
    core: String,
    family: String,
    outputs: Vec<Output>,
}

#[derive(Deserialize)]
struct Output {
    text: String,
    /// The spellings that reached this output.
    from: Vec<String>,
}

// ---------------------------------------------------------------------------
// Providers — the same construction the generator and the dump use
// ---------------------------------------------------------------------------

fn padded(core: &str) -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{core}{pad}")
}

fn reference_for(axis: &str, core: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    if axis == "g" {
        provider.add_genomic_sequence(GENOMIC_CONTIG, padded(core));
        return provider;
    }
    let tx_len = core.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
    let transcript = Transcript::new(
        TX_ACCESSION.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        core.to_string(),
        Some(CDS_START),
        Some(CDS_END),
        vec![exon],
        Some(TX_CONTIG.to_string()),
        Some(g_start),
        Some(g_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    provider.add_genomic_sequence(TX_CONTIG, padded(core));
    provider.add_transcript(transcript);
    provider
}

// ---------------------------------------------------------------------------
// Measurement
// ---------------------------------------------------------------------------

/// What asking for a key produced, with the three cases kept apart.
///
/// Only [`KeyOutcome::Refused`] is a fact about `spdi_key`. The other two are
/// facts about the *corpus* or about a bug, and folding them into the refusal
/// count would report a measurement failure as an API limit — the shape
/// `CLAUDE.md` names under "a generator must account for what it dropped": a
/// fallible step whose failure is representable as a legitimate value, so a
/// partial run and a clean run write indistinguishable output.
enum KeyOutcome {
    /// `spdi_key` produced a key.
    Keyed(String),
    /// `spdi_key` declined, for the stated reason. A documented refusal class.
    Refused(String),
    /// The question could not be asked at all: the row does not parse, or
    /// keying it panicked. Not a refusal, and not evidence about the key.
    Failed(&'static str),
}

/// What the key concluded about one class.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Verdict {
    /// Every distinct output keyed, and to one key.
    Grouped,
    /// Every distinct output keyed, to two or more keys.
    Split,
    /// Some output has no key, so the class has no answer.
    Undecided,
    /// Some output could not be asked — unparseable, or it panicked. Ranked
    /// apart from `Undecided` because that one is the key declining and this one
    /// is the harness failing, and only the first says anything about `spdi_key`.
    Unmeasured,
}

impl Verdict {
    fn label(self) -> &'static str {
        match self {
            Verdict::Grouped => "grouped",
            Verdict::Split => "split",
            Verdict::Undecided => "undecided",
            Verdict::Unmeasured => "unmeasured",
        }
    }
}

/// One class's outcome, keyed and censused.
struct Measured {
    id: String,
    family: String,
    verdict: Verdict,
    /// Distinct normalized outputs the class reached.
    outputs: usize,
    /// Distinct keys those outputs produced, rendered, for an exemplar line.
    keys: BTreeSet<String>,
    /// Outputs that had no key, each paired with the reason.
    ///
    /// The reason is the whole value of an `undecided` row. Without it a reader
    /// has to guess which of [`spdi_key`]'s documented refusal classes fired, and
    /// a guess is what turns "declined for a stated reason" into "the key is
    /// incomplete in some unspecified way".
    declined: Vec<(String, String)>,
    /// Outputs and spellings the harness could not ask about, each paired with
    /// why — `unparseable` or `panicked`.
    ///
    /// Kept out of [`Measured::declined`] on purpose. A row the corpus cannot
    /// parse says nothing about `spdi_key`, and counting it as a refusal both
    /// overstates the refusal census and hides a broken corpus behind a number
    /// that looks like a finding.
    failed: Vec<(String, String)>,
    /// The same question asked of the class's raw *spellings* rather than its
    /// normalized outputs — `None` when any spelling declined.
    ///
    /// Tracked because the two are different claims and only one is about the
    /// normalizer. A consumer keying its call set directly, never normalizing,
    /// depends on the spellings agreeing; a consumer keying stored canonical
    /// forms depends on the outputs agreeing.
    spellings_grouped: Option<bool>,
}

fn measure(record: &Record) -> Measured {
    let provider = reference_for(&record.axis, &record.core);
    let key_of = |text: &str| {
        let Ok(variant) = parse_hgvs(text) else {
            return KeyOutcome::Failed("unparseable");
        };
        // The corpus is built to contain the shapes that break things, so one
        // panicking row must not take the whole measurement with it.
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            spdi_key(&variant, &provider)
        })) {
            Ok(Some(key)) => KeyOutcome::Keyed(key.to_string()),
            // The reason is read here rather than on a second pass, so the
            // refusal and its explanation come from one parse of one row.
            // `spdi_key` is deliberately reason-free — it is the bucketing API —
            // and `canonical_spdi` is the same decision carrying the described
            // error, which is what the module docs point a caller at.
            Ok(None) => KeyOutcome::Refused(
                match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    canonical_spdi(&variant, &provider)
                })) {
                    Ok(Err(e)) => e.to_string(),
                    Ok(Ok(_)) => "no reason — it keys on a second attempt".to_string(),
                    Err(_) => "panicked while explaining the refusal".to_string(),
                },
            ),
            Err(_) => KeyOutcome::Failed("panicked"),
        }
    };

    let mut keys = BTreeSet::new();
    let mut declined = Vec::new();
    let mut failed = Vec::new();
    for output in &record.outputs {
        match key_of(&output.text) {
            KeyOutcome::Keyed(key) => {
                keys.insert(key);
            }
            KeyOutcome::Refused(reason) => declined.push((output.text.clone(), reason)),
            KeyOutcome::Failed(why) => failed.push((output.text.clone(), why.to_string())),
        }
    }

    // The verdict is a statement about the *outputs*, so it is fixed here —
    // before the spelling pass, whose failures are reported but must not restate
    // themselves as a verdict about something they do not describe.
    let outputs_unmeasurable = !failed.is_empty();

    // The spelling pass gets the same three-way split, for the same reason: an
    // unparseable spelling must not be reported as "the spellings do not group".
    let mut spelling_keys = BTreeSet::new();
    let mut spellings_unanswerable = false;
    for spelling in record.outputs.iter().flat_map(|output| output.from.iter()) {
        match key_of(spelling) {
            KeyOutcome::Keyed(key) => {
                spelling_keys.insert(key);
            }
            KeyOutcome::Refused(_) => spellings_unanswerable = true,
            KeyOutcome::Failed(why) => {
                spellings_unanswerable = true;
                failed.push((spelling.clone(), why.to_string()));
            }
        }
    }
    let spellings_grouped = if spellings_unanswerable {
        None
    } else {
        Some(spelling_keys.len() <= 1)
    };

    // Ordered so a class the harness could not measure never masquerades as one
    // the key declined, and neither is counted as grouped or split.
    let verdict = if outputs_unmeasurable {
        Verdict::Unmeasured
    } else if !declined.is_empty() {
        Verdict::Undecided
    } else if keys.len() <= 1 {
        Verdict::Grouped
    } else {
        Verdict::Split
    };

    Measured {
        id: record.id.clone(),
        family: record.family.clone(),
        verdict,
        outputs: record.outputs.len(),
        keys,
        declined,
        failed,
        spellings_grouped,
    }
}

fn render(measured: &[Measured], examples: usize) -> String {
    let mut out = String::new();
    let total = measured.len();
    let multi_output = measured.iter().filter(|m| m.outputs > 1).count();

    let _ = writeln!(out, "divergent classes read: {total}");
    if multi_output == 0 {
        let _ = writeln!(
            out,
            "\nVACUOUS: no class in this input has two or more distinct normalized outputs, \
             so nothing here measures whether the key reconciles a divergence. Regenerate \
             the divergence dump before reading any figure below."
        );
    } else {
        let _ = writeln!(
            out,
            "classes with >1 distinct output (the population the claim is about): {multi_output}"
        );
    }

    let mut by_verdict: BTreeMap<Verdict, usize> = BTreeMap::new();
    for m in measured {
        *by_verdict.entry(m.verdict).or_default() += 1;
    }
    let _ = writeln!(out, "\nverdict over all {total} classes:");
    for verdict in [Verdict::Grouped, Verdict::Split, Verdict::Undecided] {
        let count = by_verdict.get(&verdict).copied().unwrap_or(0);
        let share = if total == 0 {
            0.0
        } else {
            100.0 * count as f64 / total as f64
        };
        let _ = writeln!(out, "  {:<10} {count:>7}  ({share:.1}%)", verdict.label());
    }

    // The same question asked of the raw spellings. Printed separately because a
    // consumer that never normalizes depends on this number and not the one
    // above, and conflating them would let one carry the other's evidence.
    let spelling_answered = measured
        .iter()
        .filter(|m| m.spellings_grouped.is_some())
        .count();
    let spelling_grouped = measured
        .iter()
        .filter(|m| m.spellings_grouped == Some(true))
        .count();
    let _ = writeln!(
        out,
        "\nraw spellings (not normalized) that key to one bucket: \
         {spelling_grouped} of {spelling_answered} answerable classes"
    );

    let mut families: BTreeMap<(&str, Verdict), usize> = BTreeMap::new();
    for m in measured {
        *families.entry((m.family.as_str(), m.verdict)).or_default() += 1;
    }
    let _ = writeln!(out, "\nby family:");
    let names: BTreeSet<&str> = measured.iter().map(|m| m.family.as_str()).collect();
    for name in names {
        let at = |v: Verdict| families.get(&(name, v)).copied().unwrap_or(0);
        let _ = writeln!(
            out,
            "  {name:<28} grouped {:>6}  split {:>6}  undecided {:>6}  unmeasured {:>6}",
            at(Verdict::Grouped),
            at(Verdict::Split),
            at(Verdict::Undecided),
            at(Verdict::Unmeasured)
        );
    }

    // Census the refusals by reason, so the report answers "which documented
    // refusal class fired" without the reader sampling exemplars and
    // extrapolating. Reasons are the error's own trailing clause, since the
    // leading part names the variant and would make every row unique.
    let mut reasons: BTreeMap<&str, usize> = BTreeMap::new();
    for m in measured {
        for (_, reason) in &m.declined {
            let tail = reason
                .rsplit_once(" — ")
                .map_or(reason.as_str(), |(_, t)| t);
            *reasons.entry(tail).or_default() += 1;
        }
    }
    if !reasons.is_empty() {
        let _ = writeln!(out, "\nrefusals by reason:");
        for (reason, count) in &reasons {
            let _ = writeln!(out, "  {count:>6}  {reason}");
        }
    }

    // Reported separately from the refusals above, and loudly: these rows were
    // never asked the question, so every count in this report is over a smaller
    // population than the header implies. Silence here is what would let a
    // partially-measured run read as a clean one.
    let mut failures: BTreeMap<&str, usize> = BTreeMap::new();
    for m in measured {
        for (_, why) in &m.failed {
            *failures.entry(why.as_str()).or_default() += 1;
        }
    }
    if !failures.is_empty() {
        let total: usize = failures.values().sum();
        let _ = writeln!(
            out,
            "\nUNMEASURED: {total} row(s) could not be asked, so they are excluded from every \
             count above — they are not refusals:"
        );
        for (why, count) in &failures {
            let _ = writeln!(out, "  {count:>6}  {why}");
        }
        for m in measured
            .iter()
            .filter(|m| !m.failed.is_empty())
            .take(examples)
        {
            for (text, why) in &m.failed {
                let _ = writeln!(out, "      {text}\n        {why}");
            }
        }
    }

    for verdict in [Verdict::Split, Verdict::Undecided] {
        let exemplars: Vec<&Measured> = measured
            .iter()
            .filter(|m| m.verdict == verdict)
            .take(examples)
            .collect();
        if exemplars.is_empty() {
            continue;
        }
        let _ = writeln!(out, "\n{} exemplars:", verdict.label());
        for m in exemplars {
            let _ = writeln!(out, "  {} [{}]", m.id, m.family);
            for key in &m.keys {
                let _ = writeln!(out, "      key {key}");
            }
            for (text, reason) in &m.declined {
                let _ = writeln!(out, "      no key for {text}\n        because {reason}");
            }
        }
    }
    out
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let text = match std::fs::read_to_string(&cli.divergences) {
        Ok(text) => text,
        Err(e) => {
            eprintln!(
                "error: reading {}: {e}\nGenerate it with:\n  cargo run --release --features dev \
                 --example dump_confluence_divergences -- --out {}",
                cli.divergences.display(),
                cli.divergences.display()
            );
            return ExitCode::FAILURE;
        }
    };
    let records: Vec<Record> = match serde_json::from_str(&text) {
        Ok(records) => records,
        Err(e) => {
            eprintln!("error: parsing {}: {e}", cli.divergences.display());
            return ExitCode::FAILURE;
        }
    };

    let measured: Vec<Measured> = records.iter().map(measure).collect();
    print!("{}", render(&measured, cli.examples));

    // Refuse, do not report. A row the harness could not ask about is a hole in
    // the measurement, and exiting 0 with the census printed is precisely the
    // "partial run and clean run are indistinguishable" failure this repository
    // treats as a defect in the generator rather than a caveat in its output.
    let unmeasured: usize = measured.iter().map(|m| m.failed.len()).sum();
    if unmeasured > 0 {
        eprintln!(
            "error: {unmeasured} row(s) could not be measured (unparseable, or keying panicked); \
             the counts above are over a smaller population than the input"
        );
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(axis: &str, core: &str, outputs: &[(&str, &[&str])]) -> Record {
        Record {
            id: "test".to_string(),
            axis: axis.to_string(),
            core: core.to_string(),
            family: "test-family".to_string(),
            outputs: outputs
                .iter()
                .map(|(text, from)| Output {
                    text: (*text).to_string(),
                    from: from.iter().map(|s| (*s).to_string()).collect(),
                })
                .collect(),
        }
    }

    /// The core of the measurement: two normalized outputs that partition one
    /// edit differently are `Grouped`, not `Split`.
    ///
    /// The core is written so its padded contig has the two spellings' bases
    /// where the coordinates say they are — `reference_for` prepends
    /// `PAD_OFFSET` bases, and the `g.` coordinates below are 1-based over the
    /// padded contig.
    #[test]
    fn two_partitions_of_one_edit_measure_as_grouped() {
        let core = "ATTACAG";
        let start = PAD_OFFSET + 1;
        let spanning = format!("{GENOMIC_CONTIG}:g.{start}_{}delinsGGCTA", start + 4);
        let decomposed = format!(
            "{GENOMIC_CONTIG}:g.[{start}A>G;{}T>G;{}T>C;{}A>T;{}C>A]",
            start + 1,
            start + 2,
            start + 3,
            start + 4
        );
        let record = record(
            "g",
            core,
            &[(&spanning, &["a"][..]), (&decomposed, &["b"][..])],
        );
        let measured = measure(&record);
        assert_eq!(
            measured.verdict,
            Verdict::Grouped,
            "keys {:?}",
            measured.keys
        );
        assert_eq!(measured.keys.len(), 1);
    }

    /// The control. Without it a `measure` that answered `Grouped` for
    /// everything would satisfy the test above, which is the exact failure the
    /// whole key is supposed to make impossible.
    #[test]
    fn two_genuinely_different_edits_measure_as_split() {
        let start = PAD_OFFSET + 1;
        let one = format!("{GENOMIC_CONTIG}:g.{start}A>G");
        let other = format!("{GENOMIC_CONTIG}:g.{start}A>C");
        let record = record("g", "ATTACAG", &[(&one, &["a"][..]), (&other, &["b"][..])]);
        assert_eq!(measure(&record).verdict, Verdict::Split);
    }

    /// An output with no key makes the class `Undecided`, never `Split` — the
    /// two are different findings and folding them together would report a
    /// missing answer as a wrong one.
    #[test]
    fn an_unkeyable_output_is_undecided_not_split() {
        let start = PAD_OFFSET + 1;
        let keyed = format!("{GENOMIC_CONTIG}:g.{start}A>G");
        let record = record(
            "g",
            "ATTACAG",
            &[(&keyed, &["a"][..]), ("NC_ABSENT.1:g.5A>G", &["b"][..])],
        );
        let measured = measure(&record);
        assert_eq!(measured.verdict, Verdict::Undecided);
        assert_eq!(measured.declined.len(), 1);
    }

    /// An output the harness cannot even parse is `Unmeasured`, never
    /// `Undecided`, and is kept out of the refusal census.
    ///
    /// The two look identical in a report that folds them together, and they are
    /// opposite findings: `Undecided` says `spdi_key` declined — a fact about the
    /// key — while this says the row never reached it. Counting an unparseable
    /// row as a refusal overstates what the key refuses and hides a broken corpus
    /// behind a number that reads like a result.
    #[test]
    fn an_unparseable_output_is_unmeasured_not_a_refusal() {
        let start = PAD_OFFSET + 1;
        let keyed = format!("{GENOMIC_CONTIG}:g.{start}A>G");
        // Real spellings, so the only unmeasurable row is the one under test.
        let record = record(
            "g",
            "ATTACAG",
            &[
                (&keyed, &[keyed.as_str()][..]),
                ("this is not HGVS at all", &[keyed.as_str()][..]),
            ],
        );
        let measured = measure(&record);
        assert_eq!(measured.verdict, Verdict::Unmeasured);
        assert_eq!(
            measured.declined.len(),
            0,
            "an unparseable row is not a refusal by spdi_key"
        );
        assert_eq!(measured.failed.len(), 1, "failed was {:?}", measured.failed);
        assert_eq!(measured.failed[0].1, "unparseable");
    }

    /// And the run says so on the way out, rather than printing a census over a
    /// population it silently shrank.
    #[test]
    fn the_report_names_the_rows_it_could_not_measure() {
        let start = PAD_OFFSET + 1;
        let keyed = format!("{GENOMIC_CONTIG}:g.{start}A>G");
        let measured = vec![measure(&record(
            "g",
            "ATTACAG",
            &[
                (&keyed, &["a"][..]),
                ("this is not HGVS at all", &["b"][..]),
            ],
        ))];
        let report = render(&measured, 3);
        assert!(report.contains("UNMEASURED"), "report was:\n{report}");
        assert!(
            report.contains("they are not refusals"),
            "the report must not let an unmeasured row read as a refusal:\n{report}"
        );
    }

    /// A run whose input has no multi-output class says so, because a `0 split`
    /// from such a run is indistinguishable from success.
    #[test]
    fn a_single_output_input_reports_itself_vacuous() {
        let start = PAD_OFFSET + 1;
        let only = format!("{GENOMIC_CONTIG}:g.{start}A>G");
        let measured = vec![measure(&record("g", "ATTACAG", &[(&only, &["a"][..])]))];
        assert!(render(&measured, 0).contains("VACUOUS"));
    }
}
