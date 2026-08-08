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
//! `start` in `101, 108, …, 2996` and every `len` in `{4, 6, 8, 12, 16}`, and
//! run all 2,075 rows against a committed slice of the prepared reference.
//! **125 of them cross an exon junction**, which is the geometry #1478 records
//! as structurally unreachable elsewhere in the tree.
//!
//! # Three layers, none of them redundant
//!
//! 1. [`every_authored_inversion_produces_its_pinned_output`] — the exact
//!    string for all 2,075 rows, from `tests/fixtures/inversion-sweep/cases.tsv`.
//!    This is the only layer that notices the 3'-shift *arithmetic* moving:
//!    `c.101_108inv -> c.102_107inv` silently becoming `-> c.103_106inv` keeps
//!    the same outcome class, keeps the census, and never splits — so layers 2
//!    and 3 stay green while a string a downstream consumer stores has moved,
//!    which is exactly the representation change this project requires to be
//!    declared.
//! 2. [`no_authored_inversion_leaves_the_inversion_family`] — the invariant.
//!    It is what makes a *careless re-blessing* of layer 1 still fail: someone
//!    who regenerates `cases.tsv` to make a red build green cannot thereby
//!    accept a split or a retype.
//! 3. [`the_outcome_census_is_unchanged`] — the three bucket counts as named
//!    constants, so a bulk re-bless that shifts the **distribution** is visible
//!    in a diff and has to be argued for.
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
//!   invariant and the sequence grounding below, not the strings themselves. A
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
//! kind to the 481 within-exon narrowings. Recorded here because the output
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
//! under-serving slice reproduces the pinned answer for 1,566 of 2,075 rows by
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
    classify, reverse_complement, sweep_inputs, PinnedCases, SweepInput, SweepOutcome, CASES_PATH,
    SPEC_DIR, TRANSCRIPT, WINDOWS_PATH,
};
use ferro_hgvs::conformance::reference_window::{WindowFixture, WindowProvider};
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::{NormalizeConfig, Normalizer, ReferenceProvider, ShuffleDirection};

/// Quoted in every failure whose remedy is a regeneration.
const REGENERATE: &str = "cargo run --features dev --example extract_inversion_sweep_windows -- \
                          --manifest <manifest>";

// ---------------------------------------------------------------------------
// The census
// ---------------------------------------------------------------------------

/// Rows returned character-for-character unchanged: already the 3'-most,
/// minimal spelling.
const CENSUS_UNCHANGED: usize = 1_566;
/// Rows that stayed a lone `inv` at a different span — narrowed to the minimal
/// block whose replacement is still the reverse complement, or 3'-shifted.
const CENSUS_SHIFTED_STILL_INV: usize = 483;
/// Rows that became a no-op `=`: the authored span is its own reverse
/// complement, so inverting it changes nothing.
const CENSUS_PALINDROMIC_NO_OP: usize = 26;
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
fn normalize_variant<P: ReferenceProvider + Clone>(
    provider: &P,
    input: &str,
) -> (ferro_hgvs::HgvsVariant, String) {
    let config = ErrorConfig::lenient();
    let parsed = parse_hgvs_with_config(input, config.clone())
        .unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    let normalize_config = NormalizeConfig::for_entry_point(ShuffleDirection::ThreePrime, config);
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
                outcome,
            }
        })
        .collect();
    (rows, provider)
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
// Layer 2: the invariant
// ---------------------------------------------------------------------------

/// An authored `inv` normalizes to an `inv` or to a palindromic `=`, and to
/// nothing else. In particular it is **never** split into allele members.
///
/// > Inversion: a sequence change where, compared to a reference sequence,
/// > **more than one nucleotide** replacing the original sequence is the
/// > reverse complement of the original sequence.
/// >   — `assets/hgvs-nomenclature/docs/recommendations/DNA/inversion.md:5`
///
/// The authored description asserts exactly that relation between a span and
/// its replacement. Normalization may move where the span sits (the 3'rule,
/// `DNA/inversion.md:17`) and may narrow it to the minimal block that still
/// carries the relation, and where the whole span is its own reverse complement
/// it may conclude nothing changed at all. What it may not do is conclude the
/// replacement is something *other* than that span's reverse complement — which
/// is what a split into members, or a retype to a substitution / deletion /
/// duplication, asserts.
///
/// This survives a re-blessing of layer 1: regenerating `cases.tsv` to make a
/// red build green cannot make a split conformant.
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
    let fixture = window_fixture();
    let transcript = fixture
        .transcripts
        .iter()
        .find(|tx| tx.id == TRANSCRIPT)
        .unwrap_or_else(|| panic!("{WINDOWS_PATH} carries no {TRANSCRIPT} — `{REGENERATE}`"));
    let bases = transcript
        .sequence
        .as_deref()
        .unwrap_or_else(|| panic!("{TRANSCRIPT} carries no sequence — `{REGENERATE}`"))
        .to_ascii_uppercase();
    let cds_start = transcript
        .cds_start
        .unwrap_or_else(|| panic!("{TRANSCRIPT} carries no CDS start — `{REGENERATE}`"));

    // `c.n` is 1-based transcript position `cds_start + n - 1`.
    let tx_pos = |c: u64| cds_start + c - 1;
    let within_one_exon = |input: &SweepInput| {
        let (lo, hi) = (tx_pos(input.start), tx_pos(input.end()));
        transcript
            .exons
            .iter()
            .any(|exon| exon.start <= lo && hi <= exon.end)
    };

    let (rows, provider) = run_sweep();
    assert_the_slice_answered(&provider, rows.len());

    let mut crossing = 0usize;
    let mut palindromic = 0usize;
    let mut disagreements = Vec::new();
    for (input, row) in sweep_inputs().into_iter().zip(&rows) {
        if !within_one_exon(&input) {
            crossing += 1;
            continue;
        }
        let (lo, hi) = (tx_pos(input.start), tx_pos(input.end()));
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

// ---------------------------------------------------------------------------
// Layer 3: the census
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
    let observed = (
        count(&SweepOutcome::Unchanged),
        count(&SweepOutcome::ShiftedInversion),
        count(&SweepOutcome::PalindromicNoOp),
    );
    assert_eq!(
        observed,
        (
            CENSUS_UNCHANGED,
            CENSUS_SHIFTED_STILL_INV,
            CENSUS_PALINDROMIC_NO_OP
        ),
        "the outcome census moved (unchanged, shifted-still-inv, palindromic-no-op). Every row \
         that changed bucket is a normalized string that moved, so this needs a \
         `Representation-Change:` declaration, not a new constant"
    );
    assert_eq!(
        observed.0 + observed.1 + observed.2,
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
    for published in ["NM_004006.2:c.5657_5660inv", "NM_004006.2:c.4145_4160inv"] {
        let actual = normalize(&provider, published);
        assert_eq!(
            actual, published,
            "the spec prints `{published}` as the description to use \
             (DNA/inversion.md:27-34); ferro rewrote it"
        );
    }
    assert_the_slice_answered(&provider, 2);
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
/// normalization converts into "input unchanged", which for 1,566 of these rows
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
                 than the recorded one bails and returns its input unchanged — which for 1,566 of \
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
    let Some(manifest) = std::env::var_os("FERRO_MANIFEST").map(PathBuf::from) else {
        eprintln!(
            "hermetic_normalization_matches_the_prepared_reference: skipping — no FERRO_MANIFEST. \
             The per-PR gates are every_authored_inversion_produces_its_pinned_output and \
             no_authored_inversion_leaves_the_inversion_family, which never skip."
        );
        return;
    };
    if !manifest.exists() {
        eprintln!("hermetic_normalization_matches_the_prepared_reference: skipping — FERRO_MANIFEST does not exist");
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
