//! 2,075 authored DNA inversions on real reference bases, run hermetically.
//!
//! # The question
//!
//! `DNA/inversion.md` is one of the few recommendation pages whose whole
//! subject is sequence-dependent: the 3'rule (`:17`) reads the bases it would
//! shift over, and the definition itself (`:5`) is a statement about a span and
//! its reverse complement. Nothing in this repository was exercising that at
//! scale against real bases. The manifest-backed axes skip without
//! `FERRO_MANIFEST`, which in CI is always; the harvested spec fixture
//! normalizes through an empty `MockProvider`, so no sequence-dependent rule
//! can fire; and `examples/dump_normalized_corpus.rs` builds single-exon
//! transcripts (#1478), so no generated row can cross an exon junction.
//!
//! So: sweep `NM_004006.2:c.{start}_{start + len - 1}inv` — the DMD coding
//! sequence the spec uses for its own worked inversion examples — for every
//! `start` in `101, 108, …, 2999` and every `len` in `{4, 6, 8, 12, 16}`, and
//! run all 2,075 rows against a committed slice of the prepared reference.
//! **125 of them cross an exon junction**, which is the geometry #1478 records
//! as structurally unreachable elsewhere in the tree.
//!
//! # Four layers, none of them redundant
//!
//! 1. [`every_authored_inversion_produces_its_pinned_output`] — the exact
//!    string for all 2,075 rows, from `tests/fixtures/inversion-sweep/cases.tsv`.
//!    This is the only layer that notices the 3'-shift *arithmetic* moving:
//!    `c.101_108inv -> c.102_107inv` silently becoming `-> c.103_106inv` keeps
//!    the same outcome class, keeps the census, and denotes different bases
//!    only if the shift is wrong — so the other layers can stay green while a
//!    string a downstream consumer stores has moved, which is exactly the
//!    representation change this project requires to be declared.
//! 2. [`every_within_exon_output_reproduces_the_authored_inversion_bases`] — the
//!    sequence oracle, and the layer that carries the weight. It applies the
//!    output's own members to the committed bases and compares against the
//!    authored span replaced by its reverse complement, so an output that
//!    denotes different bases fails whatever shape it takes. Neither side
//!    calls the normalizer.
//! 3. [`no_authored_inversion_leaves_the_inversion_family`] — the shape
//!    invariant. It is what makes a *careless re-blessing* of layer 1 still
//!    fail: someone who regenerates `cases.tsv` to make a red build green
//!    cannot thereby accept a retype or a member that changes the sequence's
//!    length.
//! 4. [`the_outcome_census_is_unchanged`] — the four bucket counts as named
//!    constants, so a bulk re-bless that shifts the **distribution** is visible
//!    in a diff and has to be argued for.
//!
//! # A split was once a legal outcome here; #1703 removed it
//!
//! This corpus was authored against a tree where the `c.` axis did not run the
//! sequence-derived partitioner, and its invariant read "an authored `inv` is
//! never split into allele members". **#1484 widened that gate to
//! `c.`/`n.`/`r.`**, and 155 of these 2,075 rows repartitioned. **#1706's
//! post-hoc run scan put 62 of those back together, leaving 93**, and **#1835's
//! flip of the default partition rule to `canonical-coalesced` took one more,
//! leaving 92**.
//!
//! **That direction was then decided the other way, against every one of them,
//! and #1703 implemented it — so [`CENSUS_REPARTITIONED`] is now 0.** The
//! `decided` record `rulings[whole-span-reverse-complement-types-as-inv]`
//! (2026-08-13) types a whole-span reverse complement `inv` uniformly, keyed on
//! the span rather than on the type of the competing partition, superseding
//! `rulings[inversion-vs-two-delins-76-83]`'s note that "#1230's substitution
//! case is untouched and still splits" (`general.md:56` is declined on
//! `rulings[adjudication-precedence-order]`'s E1). `whole_span_reverse_complement`
//! in `merge.rs` types the exact whole-span shape before the competitor gate, so
//! no authored inversion in this corpus is repartitioned any more.
//!
//! **The population that moved was the all-substitutions case, and it is not the
//! one #1575 costed, in either direction.** At 155 the split was 87
//! substitutions-only, 43 an `inv` core with flanking substitutions, and 25
//! members including a `delins`. #1706 absorbed all 43 of the `inv`-core rows and
//! 19 of the 25, so **all 93 that then remained on `live` (92 on the default arm)
//! were the all-substitutions case** — the population `general.md:56` was read as
//! deciding, and the one this ruling overturns. Six of them used to *render* a
//! `delins` member (`c.367_374inv -> c.[367_369delinsTTC;372G>T;374A>C]`): the
//! coding codon carve-out merging two isolated substitutions after the partition,
//! not a `delins` in the partition itself — those are now one `inv` too. So the
//! ruling reached a **narrower** row set than the 155 the record was costed on
//! and a **wider** one than #1575's 25; the 92 moved to [`CENSUS_UNCHANGED`] (75)
//! and [`CENSUS_SHIFTED_STILL_INV`] (17). Every one still denotes the right bases
//! (layer 2 checks each), so this was a question about which legal spelling ships
//! rather than about correctness.
//!
//! # Which rows are adjudications, and which are characterization
//!
//! Say it plainly, because the distinction is the repo's own policy
//! (`CLAUDE.md`: "a test that merely pins today's output is not an adjudication
//! record"):
//!
//! - **Three rows are adjudications.** The spec publishes
//!   `NM_004006.2:c.5657_5660inv` and `NM_004006.2:c.4145_4160inv` outright, in
//!   `DNA/inversion.md:27-34`, quoted verbatim on
//!   [`the_spec_published_inversion_examples_produce_the_spec_form`]. The third
//!   is the `delins` spelling of the first converging on it, which
//!   `DNA/delins.md:5` settles by definitional exclusion.
//! - **The other 2,072 are characterization pins.** Their authority is the
//!   sequence oracle and the invariant, not the strings themselves. A
//!   future change may legitimately move any of them — it may not move them
//!   silently.
//!
//! # The pins are not merely self-consistent: they are grounded in the bases
//!
//! [`the_no_op_rows_are_exactly_the_self_reverse_complementary_spans`] recomputes,
//! from the committed transcript bases, which authored spans equal their own
//! reverse complement — and asserts that set is *exactly* the set of rows that
//! normalized to `=`. Both directions: no palindromic span may survive as an
//! `inv`, and no non-palindromic span may collapse to a no-op. That check
//! answers to `DNA/inversion.md:5` rather than to yesterday's output.
//!
//! # `c.a_b` spans introns, and four rows show it
//!
//! A `c.` DNA description is genomic with intron offsets, so `c.2789_2804`
//! denotes the genomic block from `c.2789` to `c.2804` **including intron 21**
//! — 12,625 bp, not 16. Ferro narrows such a block one base per *genomic* end
//! when the outermost pair is already complementary, which lands the endpoint
//! inside the intron:
//!
//! ```text
//! c.1319_1334inv -> c.1322_1332-1inv
//! c.1326_1333inv -> c.1328_1332-1inv
//! c.1473_1484inv -> c.1475_1483-1inv
//! c.2789_2804inv -> c.2790_2804-1inv
//! ```
//!
//! Each was checked by hand against the captured genomic window: for
//! `c.2789_2804inv` the block is `[32472309, 32484933]` on `NC_000023.11`, its
//! endpoints are `A`/`T` (complementary, so inverting leaves them untouched)
//! and the next pair in is `C`/`T` (not complementary, so the trim stops after
//! one base). The narrowing is symmetric on the genomic axis and identical in
//! kind to the within-exon narrowings. Recorded here because the output
//! *reads* like a defect and is not.
//!
//! # The mechanism this file must not fall for
//!
//! A window fixture captured from a normalize pass is exactly as wide as the
//! reads that pass made. A wider-reading pass runs off the end, the provider
//! errors, and `src/normalize/mod.rs:3884-3897` converts that into
//! `Ok((canonicalize_genome_variant(variant), vec![]))` — success, empty
//! warnings, **input returned unchanged**. For this corpus that is not merely
//! wrong, it is maximally deceptive: `Unchanged` is the majority outcome, so an
//! under-serving slice reproduces the pinned answer for 1,434 of 2,075 rows by
//! accident. Two mitigations, both mandatory and both asserted:
//!
//! 1. the generator's `serve_transcripts_whole` widens the transcript's window
//!    to the whole sequence, pinned hermetically by
//!    [`the_committed_slice_serves_the_transcript_whole`];
//! 2. every pass runs through
//!    [`AuditedProvider`](ferro_hgvs::conformance::audited_provider::AuditedProvider),
//!    and the zero-failed-reads assertion is evaluated **before** any verdict —
//!    "the fixture could not answer" must always outrank "ferro answered
//!    differently" — with a positive floor on served reads beside it, since
//!    reading nothing also yields zero failures.
//!
//! # Regenerating
//!
//! `cargo run --features dev --example extract_inversion_sweep_windows -- \
//!  --manifest <manifest>`, and never by hand. It rewrites both `cases.tsv` and
//! `reference-windows.json`; `--check` reports drift without writing.

use std::path::{Path, PathBuf};

use ferro_hgvs::conformance::audited_provider::AuditedProvider;
use ferro_hgvs::conformance::case_harvest::Citation;
use ferro_hgvs::conformance::inversion_sweep::{
    apply_member_edits, classify, member_edits, reverse_complement, sweep_inputs,
    sweep_normalize_config, MemberEdit, MemberReplacement, PinnedCases, SweepInput, SweepOutcome,
    CASES_PATH, SPEC_DIR, TRANSCRIPT, WINDOWS_PATH,
};
use ferro_hgvs::conformance::reference_window::{WindowFixture, WindowProvider};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::{Normalizer, ReferenceProvider};

/// Quoted in every failure whose remedy is a regeneration.
const REGENERATE: &str = "cargo run --features dev --example extract_inversion_sweep_windows -- \
                          --manifest <manifest>";

// ---------------------------------------------------------------------------
// The census
// ---------------------------------------------------------------------------

/// Rows returned character-for-character unchanged: already the 3'-most,
/// minimal spelling.
///
/// **1,434 before the post-hoc inversion run scan and 1,490 after it** (#1575):
/// 56 rows whose authored `inv` used to be shredded into an allele are now
/// returned exactly as authored.
///
/// **1,491 once `canonical-coalesced` became the default arm.** Exactly one
/// further row — `NM_004006.2:c.759_766inv` — joins them, moving out of
/// [`CENSUS_REPARTITIONED`]; see that constant for why the flip reaches it and
/// `live` does not.
///
/// **1,566 after #1703.** Implementing
/// `rulings[whole-span-reverse-complement-types-as-inv]` types every whole-span
/// reverse complement as one `inv`, so the 92 rows [`CENSUS_REPARTITIONED`]
/// counted stop being repartitioned; 75 of them return character-for-character
/// as authored and land here (the other 17 in [`CENSUS_SHIFTED_STILL_INV`]).
const CENSUS_UNCHANGED: usize = 1_566;
/// Rows that stayed a lone `inv` at a different span — narrowed to the minimal
/// block whose replacement is still the reverse complement, or 3'-shifted.
///
/// **460 before the run scan and 466 after it** (#1575): the other 6 of the 62
/// recovered rows come back as one `inv` at a narrowed span rather than at the
/// authored one.
///
/// **483 after #1703.** The 17 of [`CENSUS_REPARTITIONED`]'s 92 rows whose `inv`
/// falls at a narrowed or shifted span join here when the uniform whole-span
/// typing lands (the other 75 in [`CENSUS_UNCHANGED`]).
const CENSUS_SHIFTED_STILL_INV: usize = 483;
/// Rows that became a no-op `=`: the authored span is its own reverse
/// complement, so inverting it changes nothing.
const CENSUS_PALINDROMIC_NO_OP: usize = 26;
/// Rows repartitioned into allele members.
///
/// **Zero before #1484, 155 after it, and 93 after the post-hoc inversion run
/// scan.** #1484 widened the sequence-derived axis gate to `c.`/`n.`/`r.`, so an
/// authored `inv` whose block admits a smaller partition is re-derived into one;
/// the run scan puts 62 of those back together, which is what #1575 asked for.
///
/// **92 once `canonical-coalesced` became the default arm.** One row leaves:
/// `NM_004006.2:c.759_766inv`, which under `live` renders as the four lone
/// substitutions `c.[759G>T;761T>G;764C>A;766A>C]` and under the default arm is
/// returned as the authored `inv`.
///
/// **The `#1230` claim below is NOT re-verified for the flip arm, and must not
/// be read as if it were.** It was measured piece-by-piece on `live`, where the
/// two arms derive their partitions differently — which is the whole mechanism
/// of the flip — so "every piece is a one-base substitution" is a statement
/// about `live`'s pieces, not about the 92 the default arm now leaves here.
/// Re-deriving it needs piece-level instrumentation this corpus does not carry;
/// `dump_partitions` reads a corpus path, not one description. What DOES hold
/// unconditionally on both arms is the sequence oracle
/// [`every_within_exon_output_reproduces_the_authored_inversion_bases`] and the
/// shape invariant [`no_authored_inversion_leaves_the_inversion_family`], both
/// of which pass — so no row here has been retyped or moved off its bases.
///
/// **On `live`, all 93 were #1230's case**, which `general.md:56` ranks above
/// inversion and which `rulings[inversion-vs-two-delins-76-83]` records as
/// "untouched and still splits": every *piece* of the competing partition is a
/// one-base substitution. Measured, not assumed — every changed column in all 93
/// blocks is isolated, so there is no multi-column run for the run scan's gate to
/// admit on.
///
/// **That case used to be the endorsed one, and is now decided against.** The
/// paragraph directly above records that `general.md:56` ranks substitution
/// above inversion and that `rulings[inversion-vs-two-delins-76-83]` records
/// #1230's case as "untouched and still splits".
/// `rulings[whole-span-reverse-complement-types-as-inv]` (2026-08-13)
/// supersedes that: a whole-span reverse complement is typed `inv` uniformly,
/// keyed on the span rather than on what the competing partition looks like,
/// and `:56` is declined on `rulings[adjudication-precedence-order]`'s E1. So
/// every row this constant used to count is decided *against* — the 92 the
/// default arm left as much as the 93 `live` left.
///
/// **Now 0, implemented by #1703.** `whole_span_reverse_complement` in `merge.rs`
/// types every whole-span reverse complement as one `inv` before the competitor
/// gate is consulted, so no authored inversion in this corpus is repartitioned
/// any more: the 92 move to [`CENSUS_UNCHANGED`] (75) and
/// [`CENSUS_SHIFTED_STILL_INV`] (17). Re-blessed here with the implementing
/// change, per the instruction this doc used to carry; a move back above 0 is a
/// regression against the ruling, not drift to re-bless.
///
/// Before #1703, six of the 92 nevertheless *rendered* with a `delins` member
/// (`c.367_374inv -> c.[367_369delinsTTC;372G>T;374A>C]`); reading that as a
/// counter-example was the trap — the `delins` was the coding codon carve-out
/// (`general.md:35`, `DNA/delins.md:18`) merging two isolated substitutions one
/// base apart, applied *after* the partition, which itself held no `delins`.
/// Those rows are now one `inv`.
const CENSUS_REPARTITIONED: usize = 0;
/// Rows whose authored `c.` span crosses an exon junction. Pinned because it is
/// the geometry #1478 records as unreachable in the generated corpora — if it
/// ever reads 0, this corpus has stopped covering the thing it was built for.
const EXON_JUNCTION_CROSSING: usize = 125;

// ---------------------------------------------------------------------------
// Fixture plumbing
// ---------------------------------------------------------------------------

/// Anchored on `CARGO_MANIFEST_DIR` so the paths resolve regardless of the
/// test's working directory.
fn repo_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn window_fixture() -> WindowFixture {
    WindowFixture::from_json_path(&repo_path(WINDOWS_PATH))
        .unwrap_or_else(|e| panic!("load {WINDOWS_PATH}: {e}\nRegenerate with `{REGENERATE}`."))
}

fn pinned_cases() -> PinnedCases {
    PinnedCases::from_tsv_path(&repo_path(CASES_PATH))
        .unwrap_or_else(|e| panic!("{e}\nRegenerate with `{REGENERATE}`."))
}

/// The hermetic provider, wrapped so every failed reference read is reported
/// rather than silently becoming "input unchanged".
fn audited_provider() -> AuditedProvider<WindowProvider> {
    AuditedProvider::new(window_fixture().to_provider())
}

/// Normalize through the shipped `ferro normalize --error-mode lenient` path —
/// the same entry point the generator recorded — and render the answer.
fn normalize<P: ReferenceProvider + Clone>(provider: &P, input: &str) -> String {
    normalize_variant(provider, input).1
}

/// As [`normalize`], but keeps the parsed result so a caller can classify its
/// shape rather than re-parse the string.
///
/// The settings come from [`sweep_normalize_config`] rather than being built
/// here, so the gate cannot replay a configuration the generator never
/// recorded.
fn normalize_variant<P: ReferenceProvider + Clone>(
    provider: &P,
    input: &str,
) -> (ferro_hgvs::HgvsVariant, String) {
    let (config, normalize_config) = sweep_normalize_config();
    let parsed =
        parse_hgvs_with_config(input, config).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    let normalized = Normalizer::with_config(provider.clone(), normalize_config)
        .normalize(&parsed.result)
        .unwrap_or_else(|e| panic!("normalize `{input}`: {e}"));
    let rendered = normalized.to_string();
    (normalized, rendered)
}

/// One row of a completed pass.
struct Row {
    input: String,
    output: String,
    /// Kept alongside the rendering so the sequence oracle reads the members
    /// ferro produced rather than re-parsing its own output string.
    variant: ferro_hgvs::HgvsVariant,
    outcome: SweepOutcome,
}

/// Run the whole sweep against the hermetic slice.
///
/// Returns the rows *and* the audit, and every caller checks the audit first.
/// That ordering is the point of the type: a verdict computed from a pass that
/// ran short of reference data is not a verdict about ferro.
fn run_sweep() -> (Vec<Row>, AuditedProvider<WindowProvider>) {
    let provider = audited_provider();
    let rows = sweep_inputs()
        .into_iter()
        .map(|input| {
            let description = input.description();
            let (variant, output) = normalize_variant(&provider, &description);
            let outcome = classify(&description, &variant, &output);
            Row {
                input: description,
                output,
                variant,
                outcome,
            }
        })
        .collect();
    (rows, provider)
}

/// The committed transcript, reduced to what the two sequence-grounded tests
/// read from it.
///
/// Shared rather than derived twice: both tests must agree on which rows lie in
/// one exon, and two copies of that predicate is two ways for one of them to
/// quietly start covering a different population than its pinned count claims.
struct TranscriptBases {
    /// Stored bases, uppercased.
    bases: String,
    /// 1-based transcript position of `c.1`.
    cds_start: u64,
    /// 1-based inclusive transcript spans of the exons.
    exons: Vec<(u64, u64)>,
}

impl TranscriptBases {
    fn from_fixture(fixture: &WindowFixture) -> Self {
        let transcript = fixture
            .transcripts
            .iter()
            .find(|tx| tx.id == TRANSCRIPT)
            .unwrap_or_else(|| panic!("{WINDOWS_PATH} carries no {TRANSCRIPT} — `{REGENERATE}`"));
        Self {
            bases: transcript
                .sequence
                .as_deref()
                .unwrap_or_else(|| panic!("{TRANSCRIPT} carries no sequence — `{REGENERATE}`"))
                .to_ascii_uppercase(),
            cds_start: transcript
                .cds_start
                .unwrap_or_else(|| panic!("{TRANSCRIPT} carries no CDS start — `{REGENERATE}`")),
            exons: transcript
                .exons
                .iter()
                .map(|exon| (exon.start, exon.end))
                .collect(),
        }
    }

    /// 1-based transcript position of `c.n`.
    fn tx_pos(&self, c: u64) -> u64 {
        self.cds_start + c - 1
    }

    /// Whether the authored span lies inside a single exon.
    ///
    /// This is the gate on both sequence checks, and it is about the *authored*
    /// span rather than the output: a `c.` range spans the introns between its
    /// endpoints, so for a junction-crossing row the transcript-axis bases are
    /// not the variant's span and any comparison drawn from them is answering a
    /// different question.
    fn within_one_exon(&self, input: &SweepInput) -> bool {
        let (lo, hi) = (self.tx_pos(input.start), self.tx_pos(input.end()));
        self.exons
            .iter()
            .any(|&(start, end)| start <= lo && hi <= end)
    }
}

/// Fail unless the slice answered every read and answered a positive number of
/// them. Called *before* any verdict is trusted.
fn assert_the_slice_answered(provider: &AuditedProvider<WindowProvider>, floor: usize) {
    let failures = provider.failures();
    assert!(
        failures.is_empty(),
        "the hermetic slice failed {} reference read(s) during the pass, so every verdict below \
         is about a normalization that ran short of reference data, not about ferro's answer — \
         and `src/normalize/mod.rs` returns the input unchanged on such a failure, which for this \
         corpus is the majority pinned answer. Regenerate with `{REGENERATE}`:\n  {}",
        failures.len(),
        failures.join("\n  ")
    );
    assert!(
        provider.successful_reads() >= floor,
        "the pass served only {} reference read(s), fewer than the {floor} floor: zero failures is \
         trivially satisfiable by reading nothing, so the audit proves nothing here",
        provider.successful_reads()
    );
}

// ---------------------------------------------------------------------------
// Layer 1: the exact pins
// ---------------------------------------------------------------------------

/// Every one of the 2,075 rows produces the string committed for it.
///
/// Reported as one batch rather than row-by-row: when a normalizer change moves
/// these, *which* rows moved together is the diagnosis, and a per-row assertion
/// would report only the first.
#[test]
fn every_authored_inversion_produces_its_pinned_output() {
    let pinned = pinned_cases();
    let (rows, provider) = run_sweep();
    assert_the_slice_answered(&provider, rows.len());

    assert_eq!(
        pinned.rows.len(),
        rows.len(),
        "{CASES_PATH} holds {} rows but the sweep enumerates {} — regenerate with `{REGENERATE}`",
        pinned.rows.len(),
        rows.len()
    );

    let mut moved = Vec::new();
    for (pin, row) in pinned.rows.iter().zip(&rows) {
        // The input column is checked too: a fixture whose rows drifted out of
        // the generator's order would compare each output against the wrong
        // pin, which reads as a mass regression rather than as fixture drift.
        assert_eq!(
            pin.input, row.input,
            "{CASES_PATH} enumerates a different corpus than `sweep_inputs()` — regenerate with \
             `{REGENERATE}`"
        );
        if pin.output != row.output {
            moved.push(format!(
                "  {}\n    pinned: {}\n    actual: {}  [{}]",
                row.input,
                pin.output,
                row.output,
                row.outcome.label()
            ));
        }
    }
    assert!(
        moved.is_empty(),
        "{} of {} pinned inversion outputs moved. Each one is a normalized string a downstream \
         consumer may have stored, so this is a representation change that must be declared in \
         the PR description (`Representation-Change:`) and argued for — not re-blessed \
         silently:\n{}",
        moved.len(),
        rows.len(),
        moved.join("\n")
    );
}

// ---------------------------------------------------------------------------
// Layers 2 and 3: the sequence oracle, and the shape invariant that covers
// the rows it cannot reach
// ---------------------------------------------------------------------------

/// An authored `inv` normalizes to an `inv`, to a palindromic `=`, or to allele
/// members that are each on this transcript and each length-preserving — and to
/// nothing else.
///
/// > Inversion: a sequence change where, compared to a reference sequence,
/// > **more than one nucleotide** replacing the original sequence is the
/// > reverse complement of the original sequence.
/// >   — `assets/hgvs-nomenclature/docs/recommendations/DNA/inversion.md:5`
///
/// The authored description asserts exactly that relation between a span and
/// its replacement. Normalization may move where the span sits (the 3'rule,
/// `DNA/inversion.md:17`), may narrow it to the minimal block that still
/// carries the relation, may conclude nothing changed at all where the span is
/// its own reverse complement, and — since #1484 — may re-derive the block as
/// several members. What none of those may do is change how many bases the
/// span has: the reverse complement of a span is exactly as long as the span,
/// so a `del`, `dup`, `ins` or unbalanced `delins` asserts something `:5` does
/// not, whatever else is true of it.
///
/// **This is the weaker of the two invariants and it is deliberately kept.**
/// [`every_within_exon_output_reproduces_the_authored_inversion_bases`] is
/// strictly stronger where it applies — it compares bases, not shapes — but it
/// can only apply to the 1,950 rows that lie in one exon. This one covers all
/// 2,075, including the junction-crossing rows whose bases are not on the
/// transcript axis, and it survives a re-blessing of layer 1: regenerating
/// `cases.tsv` to make a red build green cannot make a retype conformant.
#[test]
fn no_authored_inversion_leaves_the_inversion_family() {
    let (rows, provider) = run_sweep();
    assert_the_slice_answered(&provider, rows.len());

    let violations: Vec<String> = rows
        .iter()
        .filter(|row| !row.outcome.is_conformant())
        .map(|row| {
            format!(
                "  {} -> {}  [{}] {:?}",
                row.input,
                row.output,
                row.outcome.label(),
                row.outcome
            )
        })
        .collect();
    assert!(
        violations.is_empty(),
        "{} of {} authored inversions normalized outside the inversion family. \
         `DNA/inversion.md:5` defines the replacement of the whole span as its reverse \
         complement; a split into members or a retype to another edit asserts something \
         else:\n{}",
        violations.len(),
        rows.len(),
        violations.join("\n")
    );
}

/// The grounding: the rows that became `=` are exactly the spans that equal
/// their own reverse complement, computed from the committed bases.
///
/// Without this, the census is a claim about ferro's output and nothing more.
/// With it, `=` has to be *right*: recomputing the palindrome directly from
/// `DNA/inversion.md:5`'s definition and comparing both directions catches a
/// palindromic span left as an `inv` just as it catches a non-palindromic span
/// collapsed to a no-op.
///
/// Restricted to rows whose span lies within a single exon, and deliberately so.
/// A `c.` range spans the introns between its endpoints, so for a
/// junction-crossing row the transcript-axis bases are *not* the variant's
/// span, and comparing them would be checking the wrong sequence. The count of
/// excluded rows is pinned rather than shrugged at, so the exclusion cannot
/// quietly grow to cover the whole corpus.
#[test]
fn the_no_op_rows_are_exactly_the_self_reverse_complementary_spans() {
    let transcript = TranscriptBases::from_fixture(&window_fixture());
    let bases = &transcript.bases;

    let (rows, provider) = run_sweep();
    assert_the_slice_answered(&provider, rows.len());

    let mut crossing = 0usize;
    let mut crossing_no_ops = Vec::new();
    let mut palindromic = 0usize;
    let mut disagreements = Vec::new();
    for (input, row) in sweep_inputs().into_iter().zip(&rows) {
        if !transcript.within_one_exon(&input) {
            crossing += 1;
            // A crossing row's transcript-axis bases are not its span, so it
            // cannot be judged here — but it still lands in the *census*, which
            // this test's `palindromic` count is compared against. Recording
            // them keeps that comparison honest instead of resting on an
            // unstated assumption; see the assertion below.
            if row.outcome == SweepOutcome::PalindromicNoOp {
                crossing_no_ops.push(format!("  {} -> {}", row.input, row.output));
            }
            continue;
        }
        let (lo, hi) = (
            transcript.tx_pos(input.start),
            transcript.tx_pos(input.end()),
        );
        assert!(
            lo >= 1 && (hi as usize) <= bases.len(),
            "{}: c.{}_{} maps to transcript [{lo}, {hi}], outside the {} committed bases of \
             {TRANSCRIPT}. The fixture and the sweep parameters disagree — regenerate with \
             `{REGENERATE}`",
            row.input,
            input.start,
            input.end(),
            bases.len()
        );
        let span = &bases[(lo - 1) as usize..hi as usize];
        let is_palindrome = span == reverse_complement(span);
        if is_palindrome {
            palindromic += 1;
        }
        let is_no_op = row.outcome == SweepOutcome::PalindromicNoOp;
        if is_palindrome != is_no_op {
            disagreements.push(format!(
                "  {} -> {} : span {span} {} its own reverse complement, but ferro {} it a no-op",
                row.input,
                row.output,
                if is_palindrome { "IS" } else { "is NOT" },
                if is_no_op { "DID call" } else { "did NOT call" }
            ));
        }
    }

    assert!(
        disagreements.is_empty(),
        "{} row(s) disagree with `DNA/inversion.md:5` recomputed from the committed bases. A span \
         equal to its own reverse complement is unchanged by inverting it, and only such a span \
         is:\n{}",
        disagreements.len(),
        disagreements.join("\n")
    );
    // `palindromic` below is a *within-exon* count, and it is compared against
    // `CENSUS_PALINDROMIC_NO_OP`, which `the_outcome_census_is_unchanged` takes
    // over all 2,075 rows. The two are the same number only while no crossing
    // row normalizes to `=`. Assert that rather than assume it: without this,
    // the first crossing no-op makes the census read 27 and this test read 26,
    // and the failure message below blames the committed bases — the wrong
    // diagnosis, pointing at a regeneration that would not fix anything.
    assert!(
        crossing_no_ops.is_empty(),
        "{} exon-junction-crossing row(s) normalized to a no-op `=`. Such a row is skipped by the \
         grounding loop above — its transcript-axis bases are not its span — so it cannot be \
         verified here, yet it still counts toward CENSUS_PALINDROMIC_NO_OP. Give the within-exon \
         population its own constant before accepting this:\n{}",
        crossing_no_ops.len(),
        crossing_no_ops.join("\n")
    );
    assert_eq!(
        crossing, EXON_JUNCTION_CROSSING,
        "the number of rows whose `c.` span crosses an exon junction changed. Those rows are the \
         geometry #1478 records as unreachable in the generated corpora, so a drop here is a loss \
         of the coverage this corpus was built for"
    );
    assert_eq!(
        palindromic, CENSUS_PALINDROMIC_NO_OP,
        "the number of self-reverse-complementary spans in the sweep changed, which can only mean \
         the committed transcript bases moved — regenerate with `{REGENERATE}` and re-read the \
         census"
    );
}

/// Every within-exon output must denote exactly the bases the authored inversion
/// denotes — whatever shape ferro gives it.
///
/// This is the sequence-preservation check that carries the weight. It once
/// mattered most for the repartitioned rows (a split into members that could drop
/// or double a base); #1703's uniform whole-span `inv` typing has since taken
/// [`CENSUS_REPARTITIONED`] to 0, so the rows it now guards are the single `inv`s
/// (unchanged, shifted, or narrowed). It is not a claim about ferro's
/// classification of anything: the expected sequence is built from the
/// **committed transcript bases** by replacing the authored span with its reverse
/// complement — `DNA/inversion.md:5` applied directly — and the observed sequence
/// by applying the output's own members to those same bases through
/// [`apply_member_edits`](ferro_hgvs::conformance::inversion_sweep::apply_member_edits),
/// which never calls the normalizer. A repartition that drops a base, doubles
/// one, or picks the wrong block fails here even though every member is
/// individually well-formed and the census is unmoved.
///
/// #1420's own words for what makes a partition right: "apply the edit set to
/// the reference, then re-derive the minimal canonical partition from the
/// resulting sequence". This asserts the first half of that — the edit set and
/// the inversion agree on the resulting sequence.
///
/// Restricted to rows whose authored span and whose every member lie on the
/// plain coding axis, and deliberately so: a `c.` range spans the introns
/// between its endpoints, so for a junction-crossing row the transcript-axis
/// bases are not the variant's span. Both exclusions are counted and the
/// junction-crossing count is pinned, so the restriction cannot quietly grow to
/// cover the corpus.
#[test]
fn every_within_exon_output_reproduces_the_authored_inversion_bases() {
    let transcript = TranscriptBases::from_fixture(&window_fixture());
    let (rows, provider) = run_sweep();
    assert_the_slice_answered(&provider, rows.len());

    let mut checked = 0usize;
    let mut repartitions_checked = 0usize;
    let mut crossing = 0usize;
    let mut divergences = Vec::new();
    for (input, row) in sweep_inputs().into_iter().zip(&rows) {
        if !transcript.within_one_exon(&input) {
            crossing += 1;
            continue;
        }
        // The authored inversion, applied to the committed bases. This is the
        // expected answer and it comes from `DNA/inversion.md:5`, not from
        // anything ferro returned.
        let authored = [MemberEdit {
            start: input.start,
            end: input.end(),
            replacement: MemberReplacement::ReverseComplement,
        }];
        let expected = apply_member_edits(&transcript.bases, transcript.cds_start, &authored)
            .unwrap_or_else(|e| panic!("{}: the authored span is unreadable: {e}", row.input));

        // A within-exon row whose output the coding axis cannot address is a
        // finding, not a row to skip: the block sits inside one exon, so an
        // intronic endpoint on the answer means normalization moved it
        // somewhere the authored variant never reached.
        let edits = match member_edits(&row.variant) {
            Ok(edits) => edits,
            Err(e) => {
                divergences.push(format!(
                    "  {} -> {}\n    lies in one exon, yet its output is not on the coding \
                     axis: {e}",
                    row.input, row.output
                ));
                continue;
            }
        };
        checked += 1;
        if matches!(row.outcome, SweepOutcome::Repartitioned(_)) {
            repartitions_checked += 1;
        }
        match apply_member_edits(&transcript.bases, transcript.cds_start, &edits) {
            Ok(observed) if observed == expected => {}
            Ok(observed) => {
                let at = observed
                    .bytes()
                    .zip(expected.bytes())
                    .position(|(a, b)| a != b);
                divergences.push(format!(
                    "  {} -> {}\n    first differing transcript base: {at:?}",
                    row.input, row.output
                ));
            }
            Err(e) => divergences.push(format!("  {} -> {}\n    {e}", row.input, row.output)),
        }
    }

    assert!(
        divergences.is_empty(),
        "{} of {checked} row(s) denote different bases than the inversion they came from. \
         `DNA/inversion.md:5` defines the replacement as the reverse complement of the whole \
         span, so an output that applies to different bases is not a spelling of that variant \
         whatever shape it takes:\n{}",
        divergences.len(),
        divergences.join("\n")
    );
    assert_eq!(
        crossing, EXON_JUNCTION_CROSSING,
        "the number of rows excluded for crossing an exon junction changed. Those rows are \
         excluded because a `c.` range spans its introns, so the transcript-axis bases are not \
         the variant's span — a rise here is coverage draining out of this check"
    );
    assert_eq!(
        checked + crossing,
        rows.len(),
        "every row must be either checked against the bases or counted as junction-crossing"
    );
    // Non-vacuity: the sequence check must actually have run over rows. Keyed on
    // `checked` rather than on `repartitions_checked`, because #1703's uniform
    // whole-span `inv` typing takes the repartitioned population to zero (see
    // [`CENSUS_REPARTITIONED`]) — a zero that used to read as "the check never met
    // its shape" is now the ruling's guarantee, asserted directly just below.
    assert!(
        checked > 0,
        "no within-exon row was checked against the bases, so this test proves nothing"
    );
    // The ruling's guarantee, made a live guard rather than a silent consequence:
    // `rulings[whole-span-reverse-complement-types-as-inv]` types every whole-span
    // reverse complement as one `inv`, so no authored inversion in this corpus is
    // repartitioned. A row reappearing in the `Repartitioned` bucket is a
    // regression against that ruling.
    assert_eq!(
        repartitions_checked, 0,
        "an authored inversion was repartitioned into members again — \
         `rulings[whole-span-reverse-complement-types-as-inv]` (#1703) types the whole span `inv`, \
         so this must stay 0"
    );
}

// ---------------------------------------------------------------------------
// Layer 4: the census
// ---------------------------------------------------------------------------

/// The three outcome buckets, pinned by count.
///
/// Layer 1 catches a single row moving; this catches a *bulk* re-bless that
/// keeps every row inside the invariant while shifting the distribution — 200
/// rows that used to be returned unchanged now shifting, say. That is a
/// representation change at scale and it should be impossible to land without
/// editing a number here.
#[test]
fn the_outcome_census_is_unchanged() {
    let (rows, provider) = run_sweep();
    assert_the_slice_answered(&provider, rows.len());

    let count = |wanted: &SweepOutcome| rows.iter().filter(|r| &r.outcome == wanted).count();
    let repartitioned = rows
        .iter()
        .filter(|r| matches!(r.outcome, SweepOutcome::Repartitioned(_)))
        .count();
    let observed = (
        count(&SweepOutcome::Unchanged),
        count(&SweepOutcome::ShiftedInversion),
        count(&SweepOutcome::PalindromicNoOp),
        repartitioned,
    );
    assert_eq!(
        observed,
        (
            CENSUS_UNCHANGED,
            CENSUS_SHIFTED_STILL_INV,
            CENSUS_PALINDROMIC_NO_OP,
            CENSUS_REPARTITIONED
        ),
        "the outcome census moved (unchanged, shifted-still-inv, palindromic-no-op, \
         repartitioned). Every row that changed bucket is a normalized string that moved, so this \
         needs a `Representation-Change:` declaration, not a new constant"
    );
    assert_eq!(
        observed.0 + observed.1 + observed.2 + observed.3,
        rows.len(),
        "the census buckets no longer total the corpus"
    );
}

// ---------------------------------------------------------------------------
// The adjudications: the spec publishes these answers itself
// ---------------------------------------------------------------------------

/// The two inversions `DNA/inversion.md` prints as its own worked examples are
/// returned unchanged.
///
/// > - **`NM_004006.2:c.5657_5660inv`**<br>
/// >   inversion of nucleotides `c.5657` to `c.5660` (coding DNA reference
/// >   sequence), changing `..AGG`<code class="sub">CTGA</code>`TG..` to
/// >   `..AGG`<code class="sub">TCAG</code>`TG..`.
/// >
/// > - **`NM_004006.2:c.4145_4160inv`**<br>
/// >   inversion of the 16 nucleotides from position `c.4145` to `c.4160`.
/// >   — `assets/hgvs-nomenclature/docs/recommendations/DNA/inversion.md:27-34`
///
/// Unlike the 2,072 characterization pins, these two are **adjudications**:
/// the spec names the input and prints the description to use, on the very
/// transcript the sweep is authored on, so there is nothing to interpret. A
/// divergence is a defect and a match is worth pinning. `TCAG` is indeed the
/// reverse complement of `CTGA` — the spec is checking the same relation
/// `:5` defines.
#[test]
fn the_spec_published_inversion_examples_produce_the_spec_form() {
    let provider = audited_provider();
    let published = ["NM_004006.2:c.5657_5660inv", "NM_004006.2:c.4145_4160inv"];
    // Collect first, audit second, compare third. The file's rule throughout:
    // "the fixture could not answer" outranks "ferro answered differently", and
    // asserting equality inside the loop would report the latter for what is
    // really the former.
    let actual: Vec<String> = published
        .iter()
        .map(|input| normalize(&provider, input))
        .collect();
    assert_the_slice_answered(&provider, published.len());
    for (input, actual) in published.iter().zip(&actual) {
        assert_eq!(
            actual, input,
            "the spec prints `{input}` as the description to use \
             (DNA/inversion.md:27-34); ferro rewrote it"
        );
    }
}

/// The `delins` spelling of the spec's own worked example converges on the
/// inversion.
///
/// `NM_004006.2:c.5657_5660delinsTCAG` and `NM_004006.2:c.5657_5660inv` are two
/// spellings of one variant — `TCAG` is the reverse complement of the reference
/// `CTGA` — and `DNA/delins.md:5` excludes the inversion case from `delins` by
/// definition ("**and which is not** a substitution or inversion"). So this is a
/// confluence claim with a clause behind it, not a preference: the pair must
/// agree, and it must agree *on the inversion*.
///
/// It is also the direction the repo has already decided: the
/// `inversion-vs-two-delins-76-83` ruling records `inversion.md:5` as governing
/// where the competing members are `delins`.
#[test]
fn the_delins_spelling_converges_on_the_spec_inversion() {
    let provider = audited_provider();
    let from_delins = normalize(&provider, "NM_004006.2:c.5657_5660delinsTCAG");
    let from_inv = normalize(&provider, "NM_004006.2:c.5657_5660inv");
    assert_the_slice_answered(&provider, 2);
    assert_eq!(
        from_delins, "NM_004006.2:c.5657_5660inv",
        "the delins spelling of the spec's worked example must converge on the inversion \
         (DNA/delins.md:5 excludes the inversion case from delins by definition)"
    );
    assert_eq!(
        from_delins, from_inv,
        "the two spellings of one variant did not converge"
    );
}

/// Both clauses this file rests on still say what it quotes.
///
/// A bare line number resolves against any file long enough to have that line,
/// so a spec-submodule bump could silently leave the citations above pointing
/// at unrelated prose while every other test here stayed green.
#[test]
fn the_cited_clauses_are_still_verbatim() {
    let spec_dir = repo_path(SPEC_DIR);
    assert!(
        spec_dir
            .join("docs/recommendations/DNA/inversion.md")
            .is_file(),
        "the vendored HGVS spec checkout at {SPEC_DIR} is empty. Initialise it:\n    \
         git submodule update --init {SPEC_DIR}"
    );
    for (clause, quote) in [
        (
            "docs/recommendations/DNA/inversion.md:5",
            "**more than one nucleotide** replacing the original sequence is the reverse \
             complement of the original sequence",
        ),
        (
            "docs/recommendations/DNA/inversion.md:27-34",
            "**`NM_004006.2:c.5657_5660inv`**",
        ),
        (
            "docs/recommendations/DNA/inversion.md:27-34",
            "**`NM_004006.2:c.4145_4160inv`**",
        ),
        (
            "docs/recommendations/DNA/inversion.md:17",
            "the **most 3' position** possible of the reference sequence is arbitrarily assigned \
             to have been changed",
        ),
        (
            "docs/recommendations/DNA/delins.md:5",
            "**and which is not** a substitution or inversion",
        ),
    ] {
        let citation = Citation {
            clause: clause.to_string(),
            quote: quote.to_string(),
        };
        if let Err(e) = citation.verify(&spec_dir) {
            panic!("{e}");
        }
    }
}

// ---------------------------------------------------------------------------
// Fidelity of the committed slice
// ---------------------------------------------------------------------------

/// Every pinned answer is a fixed point.
///
/// Cheap and worth having: each output is itself a legal description, so
/// re-normalizing it must be a no-op. A rule that reaches the right answer by
/// one shift too many and one back satisfies every test above and fails this.
#[test]
fn every_pinned_output_is_a_fixed_point() {
    let provider = audited_provider();
    let pinned = pinned_cases();
    let mut unstable = Vec::new();
    for row in &pinned.rows {
        let again = normalize(&provider, &row.output);
        if again != row.output {
            unstable.push(format!("  {} -> {} -> {again}", row.input, row.output));
        }
    }
    assert_the_slice_answered(&provider, pinned.rows.len());
    assert!(
        unstable.is_empty(),
        "{} of {} pinned answers are not fixed points:\n{}",
        unstable.len(),
        pinned.rows.len(),
        unstable.join("\n")
    );
}

/// The committed slice must serve `NM_004006.2` **whole**, through both entry
/// points, with bases identical to the stored record.
///
/// This encodes a real defect rather than a tidiness preference. The recorder
/// files the normalizer's `get_genomic_sequence` reads under the transcript's
/// own accession, which registers it in the fixture's `genomic` map, which makes
/// `WindowProvider::is_known_contig` true for it — and that routes `get_sequence`
/// for the transcript through the genomic-window path instead of the whole
/// stored sequence. Every read outside the captured window then errors, which
/// normalization converts into "input unchanged", which for 1,434 of these rows
/// is the pinned answer. `serve_transcripts_whole` in the generator removes the
/// failure mode; this test is what stops a future regeneration from dropping it.
#[test]
fn the_committed_slice_serves_the_transcript_whole() {
    let fixture = window_fixture();
    let provider = fixture.to_provider();
    let transcript = fixture
        .transcripts
        .iter()
        .find(|tx| tx.id == TRANSCRIPT)
        .unwrap_or_else(|| panic!("{WINDOWS_PATH} carries no {TRANSCRIPT} — `{REGENERATE}`"));
    let bases = transcript
        .sequence
        .as_ref()
        .unwrap_or_else(|| panic!("{TRANSCRIPT} carries no sequence — `{REGENERATE}`"));
    let length = bases.len() as u64;

    assert_eq!(
        provider.get_sequence_length(TRANSCRIPT).ok(),
        Some(length),
        "{TRANSCRIPT}: get_sequence_length disagrees with the stored sequence"
    );
    for (entry_point, served) in [
        ("get_sequence", provider.get_sequence(TRANSCRIPT, 0, length)),
        (
            "get_genomic_sequence",
            provider.get_genomic_sequence(TRANSCRIPT, 0, length),
        ),
    ] {
        let served = served.unwrap_or_else(|e| {
            panic!(
                "{TRANSCRIPT}: {entry_point}(0, {length}) failed against the committed slice: \
                 {e}\nThe slice does not serve this accession whole, so a pass that reads further \
                 than the recorded one bails and returns its input unchanged — which for 1,434 of \
                 these rows is the pinned answer. Regenerate with `{REGENERATE}`."
            )
        });
        assert_eq!(
            &served, bases,
            "{TRANSCRIPT}: {entry_point} served different bases than the stored transcript"
        );
    }
}

/// Ground truth: the hermetic slice must make ferro answer exactly what the
/// full prepared reference makes it answer.
///
/// Manifest-or-skip, in the Layer 2 / Layer 3 relationship this repo already
/// uses (`tests/fixtures/CORPUS_LAYOUT.md`): the hermetic tests above are what
/// block a PR, and this tier is what proves they are measuring the same thing.
/// It compares two providers running the *same* code, so it stays correct
/// through any legitimate representation change.
#[test]
fn hermetic_normalization_matches_the_prepared_reference() {
    // The per-PR gates are `every_authored_inversion_produces_its_pinned_output`
    // and `no_authored_inversion_leaves_the_inversion_family`, which never skip.
    let Some(manifest) = std::env::var_os("FERRO_MANIFEST").map(PathBuf::from) else {
        crate::common::manifest::absent(
            "inversion_sweep::hermetic_normalization_matches_the_prepared_reference",
        );
        return;
    };
    if !manifest.exists() {
        crate::common::manifest::absent(
            "inversion_sweep::hermetic_normalization_matches_the_prepared_reference \
             (FERRO_MANIFEST points at a path that does not exist)",
        );
        return;
    }
    let real = std::sync::Arc::new(
        ferro_hgvs::MultiFastaProvider::from_manifest(&manifest)
            .unwrap_or_else(|e| panic!("from_manifest({}): {e}", manifest.display())),
    );
    let hermetic = window_fixture().to_provider();

    let mut divergences = Vec::new();
    for input in sweep_inputs() {
        let description = input.description();
        let from_windows = normalize(&hermetic, &description);
        let from_manifest = normalize(&real, &description);
        if from_windows != from_manifest {
            divergences.push(format!(
                "  {description}\n    windows : {from_windows}\n    manifest: {from_manifest}"
            ));
        }
    }
    assert!(
        divergences.is_empty(),
        "{} of the sweep rows normalize differently under the committed slice than under the \
         prepared reference — the hermetic gate is not measuring ferro's real behaviour. \
         Regenerate with `{REGENERATE}`:\n{}",
        divergences.len(),
        divergences.join("\n")
    );
}
