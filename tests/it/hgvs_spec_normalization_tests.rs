//! HGVS v21.0 spec normalization regression test.
//!
//! Companion to issue #84 / #83. The fixture pins ferro's current normalize()
//! output for every variant string in the v21.0 spec. Rows whose `current`
//! diverges from `spec_expected` carry a `todo` link to the #83 audit.
//!
//! Regenerate the fixture: `cargo run --features dev --bin generate_spec_fixture`.
//! That run is also the guard on the committed overrides, and is what CI and the
//! pre-push hook use. Adding `-- --check` instead asks only whether the local
//! (gitignored) artifact is current — useful for spotting whether a code change
//! moved the fixture, never a gate.
//!
//! **Which assertions here are real oracles (#1272).** The fixture is gitignored
//! and regenerated from the code under test before every CI run, so
//! [`pinned_v21_normalization_behavior`] compares ferro's `current` against
//! ferro — it detects a stale local artifact, not a regression. The guards that
//! judge behaviour are the committed ones:
//!
//! - [`ferro_produces_the_form_the_spec_states`] — compares against
//!   `spec_expected`, which the harvester takes from the spec text and the
//!   committed overrides, with [`KNOWN_DIVERGENT_INPUTS`] as the pinned
//!   exception list.
//! - [`status_census_is_unchanged`] — per-status row counts.
//! - [`spec_equivalence_classes_converge`] — the question the statuses never
//!   ask. `preserved` means only that ferro did not *change* a string; it says
//!   nothing about whether two rows are the same variant, so the spec's own
//!   non-confluent pair (`DNA/delins.md:44-47`) sat here as two `preserved`
//!   rows and passed. Equivalence classes assert one output per variant.
//! - [`ruling_records_are_intact`] — the curated record of how the project
//!   reads the spec where the spec conflicts with itself, and of what it has
//!   deliberately *not* decided.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Debug, Deserialize)]
struct Fixture {
    rows: Vec<Row>,
    /// Curated "these spellings denote one variant" declarations, resolved by
    /// the generator. Read by [`spec_equivalence_classes_converge`].
    #[serde(default)]
    equivalence_classes: Vec<EquivalenceClass>,
    /// Curated records of contested spec readings. Read by
    /// [`ruling_records_are_intact`].
    #[serde(default)]
    rulings: Vec<Ruling>,
}

#[derive(Debug, Deserialize)]
struct EquivalenceClass {
    id: String,
    members: Vec<ClassMember>,
}

#[derive(Debug, Deserialize)]
struct ClassMember {
    /// The harvested spelling. Read by the `by_input` lookup, and named in
    /// every failure message so a mismatch identifies its member.
    input: String,
    /// The string ferro is run against — the member's row target, re-anchored
    /// onto the class's declared reference when it has one. Taken from the
    /// fixture for the same reason `input_prefixed` is: it is an *input* to the
    /// comparison, not a claim about ferro's behaviour. What this test refuses
    /// to take from the fixture is the output.
    target: String,
}

#[derive(Debug, Deserialize)]
struct Ruling {
    id: String,
    status: String,
    clauses: Vec<Citation>,
    #[serde(default)]
    governing: Option<String>,
    #[serde(default)]
    deviates_from: Vec<String>,
    rationale: String,
}

#[derive(Debug, Deserialize)]
struct Citation {
    clause: String,
    #[allow(dead_code)] // verified in the generator, carried here for readers
    quote: String,
}

#[allow(dead_code)] // many fields are present for downstream tooling, not the assertion
#[derive(Debug, Deserialize)]
struct Row {
    input: String,
    /// Default-prefixed form for bare fragments (#4). When `Some`, this is
    /// what ferro is run against — `input` stays as the verbatim spec text.
    #[serde(default)]
    input_prefixed: Option<String>,
    current: String,
    /// `None` means the spec rejects this input (sentinel for #2).
    spec_expected: Option<String>,
    status: String,
    coordinate_system: String,
    source_kind: String,
    source_paths: Vec<String>,
    #[serde(default)]
    working_group: Option<String>,
    #[serde(default)]
    todo: Option<String>,
    /// `Some(true)` means the row needs a real reference sequence to
    /// evaluate (e.g. §2.1 3'-rule shifting). Skipped until #82 lands.
    #[serde(default)]
    requires_reference: Option<bool>,
    /// Warning codes ferro currently emits for this row, sorted alphabetically.
    /// Defaults to `vec![]` for rows with no warnings.
    #[serde(default)]
    expected_warnings: Vec<String>,
}

fn fixture_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests/fixtures/grammar/hgvs_spec_normalization.json");
    p
}

fn observe(normalizer: &Normalizer<MockProvider>, input: &str) -> (String, Vec<String>) {
    match parse_hgvs(input) {
        Err(e) => (format!("parse error: {e}"), Vec::new()),
        Ok(v) => match normalizer.normalize_with_diagnostics(&v) {
            Err(e) => (format!("normalize error: {e}"), Vec::new()),
            Ok(n) => {
                let mut codes: Vec<String> =
                    n.warnings.iter().map(|w| w.code().to_string()).collect();
                codes.sort();
                codes.dedup();
                (format!("{}", n.result), codes)
            }
        },
    }
}

/// The spellings *other than* its default rendering that the spec admits for
/// whatever `input` normalizes to — empty when it does not normalize at all.
///
/// Deliberately a sibling of [`observe`] rather than an extra return value from
/// it: only the spec-conformance oracle needs these, and the other callers of
/// `observe` compare ferro against ferro, where a second spelling has no meaning.
fn observe_equivalents(normalizer: &Normalizer<MockProvider>, input: &str) -> Vec<String> {
    match parse_hgvs(input) {
        Err(_) => Vec::new(),
        Ok(v) => match normalizer.normalize_with_diagnostics(&v) {
            Err(_) => Vec::new(),
            Ok(n) => n.result.spec_equivalent_renderings(),
        },
    }
}

/// True when `s` is a bare coordinate reference (coord prefix + only position
/// characters, no edit) — mirrors the harvester's drop predicate so this test
/// pins the *outcome* (a clean fixture) independent of the generator internals.
fn is_bare_coordinate(s: &str) -> bool {
    let core = s.rsplit(':').next().unwrap_or(s);
    match core.get(0..2) {
        Some("c." | "g." | "m." | "n." | "o." | "p." | "r.") => core[2..]
            .chars()
            .all(|ch| matches!(ch, '0'..='9' | '_' | '+' | '-' | '*' | '?')),
        _ => false,
    }
}

/// #955 harvester cleanup. The generated fixture must contain only real variant
/// strings: spec-prose position references (`c.5690`), `<code class="invalid">`
/// HTML remnants, and syntax-template fragments are dropped, while variants the
/// spec split across `<code class="spotN">` highlights are reassembled intact.
#[test]
fn harvester_yields_only_real_variants() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    // No HTML remnant leaks through — no real HGVS string contains `<` or `"`.
    for row in &fx.rows {
        assert!(
            !row.input.contains('<') && !row.input.contains('"'),
            "HTML remnant leaked into fixture: {:?}",
            row.input
        );
    }

    // Outcome check over *every* row: a bare coordinate is spec prose, not a
    // variant, whenever ferro either rejects it outright or "normalizes" it by
    // inventing an identity `=` it never had (`g.12345678` →
    // `NC_000023.11:g.12345678=`). Both shapes are dropped. This must scan all
    // rows, not fixed examples. Bare coordinates ferro parses and leaves alone
    // legitimately survive — the whole-variant markers (`p.?`, `p.0`) and the
    // spec's edit-less negatives such as `c.123-65_-50`. The predicate below
    // only pins the outcome; its char-set is unit-tested independently in the
    // generator (`is_edit_less_coordinate` tests). The generator's drop is also
    // gated on `ov.is_none()` — an explicit override deliberately keeps its row;
    // no current override is a bare coordinate, so this asserts over every row;
    // if one is ever added, exempt overridden inputs here.
    for row in &fx.rows {
        if !is_bare_coordinate(&row.input) {
            continue;
        }
        assert_ne!(
            row.status, "parse-error",
            "bare coordinate reference should have been dropped, not left as parse-error: {:?}",
            row.input
        );
        let target = row.input_prefixed.as_deref().unwrap_or(&row.input);
        assert_ne!(
            row.current,
            format!("{target}="),
            "bare coordinate reference stamped with an identity `=` should have been dropped: {:?}",
            row.input
        );
    }

    // Variants the spec split across a `<code class="spotN">` highlight are
    // reassembled and parse cleanly (previously lost as severed fragments).
    let by_input: std::collections::HashMap<&str, &Row> =
        fx.rows.iter().map(|r| (r.input.as_str(), r)).collect();
    for recovered in [
        "LRG_199t1:c.11T>G",
        "LRG_199p1:p.(Val25Gly)",
        "g.123_456dup",
    ] {
        let row = by_input
            .get(recovered)
            .unwrap_or_else(|| panic!("reassembled variant missing from fixture: {recovered}"));
        assert_eq!(
            row.status, "preserved",
            "reassembled {recovered} should round-trip (preserved)"
        );
    }
}

#[test]
fn pinned_v21_normalization_behavior() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    let normalizer = Normalizer::new(MockProvider::new());
    let mut diffs: Vec<String> = Vec::new();
    let mut skipped_needs_ref = 0usize;
    let mut skipped_requires_ref = 0usize;
    let mut tested = 0usize;

    for row in &fx.rows {
        if row.status == "needs-reference" {
            skipped_needs_ref += 1;
            continue;
        }
        if row.requires_reference == Some(true) {
            skipped_requires_ref += 1;
            continue;
        }
        let target = row.input_prefixed.as_deref().unwrap_or(&row.input);
        let (observed, mut observed_warnings) = observe(&normalizer, target);
        observed_warnings.sort();
        observed_warnings.dedup();
        let mut expected_warnings = row.expected_warnings.clone();
        expected_warnings.sort();
        expected_warnings.dedup();
        if observed != row.current || observed_warnings != expected_warnings {
            diffs.push(format!(
                "  input            : {}\n    target           : {}\n    expected         : {}\n    observed         : {}\n    expected_warnings: {:?}\n    observed_warnings: {:?}\n    status           : {}",
                row.input, target, row.current, observed, expected_warnings, observed_warnings, row.status,
            ));
        }
        tested += 1;
    }

    eprintln!(
        "hgvs_spec_normalization: tested {tested}, skipped(needs-reference) {skipped_needs_ref}, skipped(requires-reference) {skipped_requires_ref}"
    );

    assert!(
        tested > 0,
        "pinned_v21_normalization_behavior exercised no cases (fixture empty or all needs-reference)"
    );

    if !diffs.is_empty() {
        let preview = diffs
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n");
        panic!(
            "{} row(s) drifted from the pinned current behavior. First {} shown.\n\n{}\n\n\
             This means your *local* fixture is stale relative to your code — the fixture is \
             gitignored and regenerated by CI, so this can only fire on a workstation. It does \
             NOT mean the behaviour change is approved.\n\n\
             Before regenerating, decide which happened:\n  \
             (a) you changed behaviour deliberately — regenerate, and check \
             STATUS_CENSUS below, which is committed and will fail if the change \
             moved any row between statuses;\n  \
             (b) you changed behaviour accidentally — this is the regression, fix it.\n\n\
             Regenerating without deciding is how a real regression gets absorbed (#1272):\n  \
             cargo run --features dev --bin generate_spec_fixture",
            diffs.len(),
            diffs.len().min(20),
            preview,
        );
    }
}

/// Committed status census for the generated normalization fixture.
///
/// **Why this exists (#1272).** `hgvs_spec_normalization.json` is gitignored and
/// regenerated from the code under test, so [`pinned_v21_normalization_behavior`]
/// compares ferro against itself: on CI, where the fixture is regenerated before
/// the run, it can only ever pass. That test is a stale-artifact detector, not a
/// regression detector.
///
/// This census is the committed half. Unlike the fixture it lives in git, so a
/// change in ferro that moves a row from one status to another fails here and
/// has to be re-blessed in a reviewable diff.
///
/// The sibling enumeration driver has carried the equivalent guard
/// (`DIVERGENCE_BUDGET` in `spec_enumeration_tests.rs`) since #1085; this fixture
/// never got one, which is the gap #1272's Scope section asks about.
///
/// **What each status means.** `preserved`: ferro's output equals the form the
/// spec states. `correctly-rejected`: the spec marks it invalid and ferro
/// rejects it. `false-acceptance`: the spec marks it invalid and ferro accepts
/// it. `diverges`: ferro's output differs from the form the spec states.
/// `needs-reference`: parses, but normalization needs real reference bases.
///
/// Regenerate the numbers by reading `summary.by_status` from the generated
/// fixture after `cargo run --features dev --bin generate_spec_fixture`.
const STATUS_CENSUS: &[(&str, usize)] = &[
    ("correctly-rejected", 99),
    // 36 -> 19 in #1079: the generator now accepts either stop-codon spelling the
    // spec sanctions (`Ter` or `*`, per `checklist.md:63`) instead of only ferro's
    // default rendering, so 17 rows that differed by that glyph alone are
    // `preserved`. They moved to `preserved` below (697 -> 714); the 934-row total
    // is unchanged, and no ferro behaviour changed — the comparison did.
    ("diverges", 19),
    ("false-acceptance", 52),
    ("needs-reference", 50),
    ("preserved", 714),
];

/// The committed census must match the generated fixture exactly.
///
/// A count mismatch is not automatically a defect — ferro improving moves rows
/// out of `diverges` and `false-acceptance` too. It is a *decision point*: the
/// numbers must be updated deliberately, in a diff a reviewer can see.
#[test]
fn status_census_is_unchanged() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    let mut counts: std::collections::BTreeMap<&str, usize> = std::collections::BTreeMap::new();
    for row in &fx.rows {
        *counts.entry(row.status.as_str()).or_default() += 1;
    }

    let mut mismatches = Vec::new();
    for (status, budget) in STATUS_CENSUS {
        let actual = counts.get(status).copied().unwrap_or(0);
        if actual != *budget {
            mismatches.push(format!("  {status}: census {budget}, actual {actual}"));
        }
    }
    // A status the census does not name at all is the more dangerous direction:
    // it is invisible to the loop above, which only checks statuses it knows.
    for status in counts.keys() {
        if !STATUS_CENSUS.iter().any(|(s, _)| s == status) {
            mismatches.push(format!("  {status}: not in STATUS_CENSUS at all"));
        }
    }

    assert!(
        mismatches.is_empty(),
        "the spec-normalization status census changed. If ferro improved, update the \
         numbers; if it regressed, this is a real defect. Either way it is a deliberate, \
         reviewable edit — that is the point of the census (#1272).\n{}",
        mismatches.join("\n")
    );
}

/// Inputs the spec states a form for, where ferro currently produces a different
/// one — the `diverges` rows, pinned by **input string** rather than by count.
///
/// Keyed by input because that is what identifies the row across a regeneration;
/// the value that changes (`current`) is precisely what must not be trusted.
///
/// Shrinking this list is an improvement and always welcome. Growing it means a
/// spec-stated form that ferro used to produce, it no longer does — read that as
/// a regression until shown otherwise.
///
/// Every entry is one of six understood classes, grouped below. They are the
/// fixture's 19 `diverges` rows, which the status census counts but — before
/// #1272 — nothing checked the *values* of.
///
/// Class 1 used to hold 19 stop-glyph rows on its own. #1079 removed them by
/// fixing the *comparison* rather than the list: a row is no longer divergent for
/// picking one of two spellings the spec presents as equals. Prefer that shape of
/// fix — an entry here should record a difference someone could act on, not one
/// the comment next to it has to argue is fine.
const KNOWN_DIVERGENT_INPUTS: &[&str] = &[
    // 1. A stop-codon glyph difference **alone** is no longer listed here. `Ter`
    //    and `*` are co-equal in the spec (`checklist.md:63`), so #1079 taught the
    //    generator to accept either rather than pinning 17 rows as divergences
    //    that the comment here had to excuse as "deliberate render policy". Two
    //    rows survive, and neither is about the stop glyph:
    //
    //    the spec states a **one-letter** amino-acid code (`W24`), which ferro
    //    renders three-letter. `ProteinRenderStyle` can emit one-letter codes, but
    //    the generator deliberately does not offer that spelling: no other
    //    spec-stated form in the corpus uses one-letter codes, so admitting it
    //    would let rows pass on a form the spec never states here;
    "NP_003997.1:p.W24*",
    //    and a cis allele whose **members are reordered** into coordinate order
    //    (`p.[Lys79Asn;Lys79Ter]` for the spec's `p.[Lys79*;Lys79Asn]`) — the same
    //    render choice as class 6 below, which no stop spelling can account for.
    //    Its parenthesised sibling `p.[(Lys79*;Lys79Asn)]` *does* now pass, since
    //    there the glyph was the only difference.
    "p.[Lys79*;Lys79Asn]",
    // 2. Trans (`];[`) alleles: ferro renders the accession outside the first
    //    bracket, the spec writes it inside.
    "NM_004006.2:c.[2376G>C];[?]",
    "NM_004006.2:r.[76a>u];[?]",
    "NP_003997.1:p.[(Ser68Arg)];[?]",
    // 3. A single-member bracket is unwrapped to the bare member — the
    //    `A2-lone-single-member-bracket` advisory rule in the enumeration's
    //    invariant catalog, applied on render.
    "chrX:g.[89555676_100352080del]",
    "p.[Ser73Arg]",
    // 4. A bare residue with no edit renders with an explicit unknown-consequence
    //    marker rather than staying bare.
    "p.1",
    "p.2",
    "p.3",
    "p.Gln990",
    "p.Ser988",
    // 5. A repeat over a single residue is rendered as an explicit range.
    "NP_002102.4:p.(Gln18)[25]",
    "p.(Gln18)[(70_80)]",
    // 6. Assorted one-offs, each a deliberate render choice:
    //    the deleted bases are dropped from `del` (spec-conformant — the spec's
    //    own checklist forbids restating them);
    "c.35delG",
    //    uncertainty parentheses are dropped from a repeat count range;
    "NC_000003.12:g.63912687AGC[(60_?)]",
    //    an unknown-position identity renders `?=` rather than `?_?`;
    "g.?_?",
    //    cis members are sorted into coordinate order, moving the anchorless
    //    `insG` to the front.
    "NC_000002.12:g.[32310435_32310710del;32310711_171827243inv;insG]",
    "chr2:g.[32310435_32310710del;32310711_171827243inv;insG]",
];

/// The real conformance oracle: ferro must produce the form **the spec states**.
///
/// **Why this is not circular, where the other tests are.** `spec_expected` does
/// not come from ferro. The harvester takes it from the spec text and the
/// committed overrides, and `classify` (`examples/common/spec_harvest.rs`) then
/// derives `status` by comparing ferro's `current` against it. So the fixture
/// already contains a genuine spec-vs-ferro verdict on every row — and until
/// #1272 nothing re-checked it live. [`pinned_v21_normalization_behavior`]
/// compares `current` against ferro, which is ferro against itself;
/// this compares the spec against ferro.
///
/// Two directions, from `classify`'s own table:
/// - `spec_expected: Some(form)` — ferro must render exactly `form`.
/// - `spec_expected: None` — the spec marks the input invalid, so ferro must
///   reject it. A successful parse is a `false-acceptance`.
///
/// Rows needing real reference bases are skipped, as elsewhere in this file:
/// the hermetic `MockProvider` cannot decide them.
#[test]
fn ferro_produces_the_form_the_spec_states() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    let normalizer = Normalizer::new(MockProvider::new());
    let mut regressions = Vec::new();
    let mut unexpected_pass = Vec::new();
    let mut checked = 0usize;

    for row in &fx.rows {
        if row.status == "needs-reference" || row.requires_reference == Some(true) {
            continue;
        }
        // `false-acceptance` is a divergence this fixture already counts, and it
        // is the `spec_expected: None` direction; the census owns it.
        if row.status == "false-acceptance" {
            continue;
        }
        let target = row.input_prefixed.as_deref().unwrap_or(&row.input);
        let (observed, _) = observe(&normalizer, target);
        let equivalents = observe_equivalents(&normalizer, target);
        let known_divergent = KNOWN_DIVERGENT_INPUTS.contains(&row.input.as_str());

        let conformant = match row.spec_expected.as_deref() {
            // The spec states a form: ferro must render it — in any spelling the
            // spec presents as equal to it (#1079). Pinning this to ferro's
            // default rendering alone reported a display-convention choice as
            // non-conformance; see `HgvsVariant::spec_equivalent_renderings`.
            Some(expected) => observed == expected || equivalents.iter().any(|r| r == expected),
            // The spec forbids the input: ferro must refuse it.
            None => observed.starts_with("parse error") || observed.starts_with("normalize error"),
        };
        checked += 1;

        match (conformant, known_divergent) {
            (false, false) => regressions.push(format!(
                "  input    : {}\n    target   : {target}\n    spec says: {:?}\n    ferro    : {observed}",
                row.input, row.spec_expected,
            )),
            (true, true) => unexpected_pass.push(row.input.clone()),
            _ => {}
        }
    }

    eprintln!("spec conformance: {checked} rows re-checked against the spec's own stated forms");
    assert!(
        checked > 0,
        "ferro_produces_the_form_the_spec_states exercised no rows"
    );
    assert!(
        regressions.is_empty(),
        "{} row(s) no longer produce the form the spec states, and are not in \
         KNOWN_DIVERGENT_INPUTS. Unlike the pinned-behaviour test above, this one \
         compares against the spec rather than against ferro, so a regeneration \
         cannot absorb it (#1272).\n{}",
        regressions.len(),
        regressions
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert!(
        unexpected_pass.is_empty(),
        "{} input(s) in KNOWN_DIVERGENT_INPUTS now match the spec — that is an \
         improvement. Remove them from the list so the guard keeps its teeth: {:?}",
        unexpected_pass.len(),
        unexpected_pass
    );
}

/// Pinned verdict for every curated equivalence class, keyed by id.
///
/// **Why a verdict table and not a failure list.** Every transition matters in
/// both directions: a class that starts disagreeing is a new non-confluence, one
/// that stops disagreeing is an improvement whose pin must be released, and one
/// that quietly becomes `not-evaluable` — because ferro started refusing a
/// member — has disarmed itself while still looking present. A one-sided list
/// catches only the first.
///
/// `non-confluent` entries are **expected and not defects to fix here**.
/// Converging them is #1235's problem, and each is a representation change for a
/// downstream consumer that keys on the normalized string, so it belongs in a
/// deliberate migration rather than in a test-infrastructure PR. What this pin
/// buys is that the *set* cannot grow silently.
///
/// The verdicts:
///
/// - `confluent` — every evaluable member normalizes to one string.
/// - `non-confluent` — two or more distinct strings for one variant.
/// - `not-evaluable` — fewer than two members can be evaluated hermetically
///   (they need real reference bases, or ferro refuses them), so the class is
///   evidence of nothing either way.
const EQUIVALENCE_CLASS_VERDICTS: &[(&str, &str)] = &[
    // The spec's own worked non-confluent pair: it calls the bracketed split an
    // "alternative description" of the delins, and ferro keeps both. Which form
    // should win is settled: `delins-merge-vs-individual-gap-two-or-more` is
    // `decided` for `:47`, the spanning delins. Ferro has not converged on it
    // yet, so this class stays pinned `non-confluent` until it does.
    ("dna-delins-vs-aligned-split-850-901", "non-confluent"),
    // Codon carve-out, gap of one. Both spellings survive normalization here —
    // but READ THAT AS A HARNESS LIMIT, NOT AS A VIOLATED RULING. The rule is
    // settled (`delins-codon-carve-out-gap-one`, status `decided`) and ferro
    // does implement it, in `apply_coding_codon_exception`. The exception is
    // frame-dependent, and this harness runs every row against
    // `MockProvider::new()` — the *empty* constructor, no transcript and so no
    // `cds_start`. `merge::axis_frame` therefore returns `None` on the `c.`
    // axis and the codon path cannot execute at all.
    //
    // Measured rather than argued, on the identical shape (two substitutions
    // one base apart inside codon 1): against a CDS-bearing provider both
    // `c.[1G>T;3T>A]` and `c.1_3delinsTTA` normalize to `c.1_3delinsTTA` — the
    // class is confluent — while against the empty provider `c.[1G>T;3T>A]`
    // comes back byte-identical. So these three verdicts pin the hermetic
    // harness's reach, and converging them needs a provider, not a rule change.
    ("dna-delins-vs-two-substitutions-235-237", "non-confluent"),
    // Same carve-out, and here the spec marks the split spelling invalid — so
    // this class overlaps the `false-acceptance` census, which counts the
    // acceptance but cannot see that it is also a second fixed point.
    ("dna-delins-vs-two-substitutions-145-147", "non-confluent"),
    ("rna-delins-vs-two-substitutions-142-144", "non-confluent"),
    // Separation ZERO, and the only class here at that separation — the three
    // above are all separation one. `substitution.md:32` marks the split
    // spelling `class="invalid"`, and ferro merges it: both members reach
    // `LRG_199t1:c.79_80delinsTT`. So this is the positive control for ruling
    // `delins-adjacent-members-when-both-consume-reference`, which is `decided`
    // and which ferro conforms to as of #1537.
    //
    // Note it is `confluent` where its separation-one siblings are not, and the
    // reason is not that the ruling is stronger: this merge needs no reading
    // frame, so the empty-provider limit documented above does not bite. Do not
    // read the contrast as evidence about the carve-out.
    ("dna-delins-vs-two-substitutions-79-80", "confluent"),
    // Spelled-out run vs repeat count for an unspecified insertion. Nothing in
    // the spec ranks the two, so no rule is being broken — but a single variant
    // still has two stable outputs.
    (
        "dna-insertion-n-run-vs-repeat-count-761-762",
        "non-confluent",
    ),
    (
        "rna-insertion-n-run-vs-repeat-count-761-762",
        "non-confluent",
    ),
    // Repeat spelled by first-unit span vs by unit sequence; again co-equal in
    // the spec and non-confluent in ferro.
    ("rna-repeat-range-vs-unit-124", "non-confluent"),
    ("rna-repeat-range-vs-unit-124-two-alleles", "non-confluent"),
    (
        "protein-insertion-xaa-run-vs-repeat-count-582-583",
        "non-confluent",
    ),
    // Positive controls: co-equal stop glyphs, which ferro *does* converge.
    // They are what shows this assertion can pass, so a change that made every
    // class disagree could not hide in a wall of expected failures.
    ("protein-extension-stop-glyph-110", "confluent"),
    ("protein-extension-stop-glyph-327", "confluent"),
    // Both members are prose shorthand the spec does not admit as a variant;
    // ferro rejects both, so nothing is compared. Pinned so that a parser change
    // which starts accepting them turns this into a live class rather than
    // passing unnoticed.
    (
        "protein-position-less-repeat-shorthand-gln21",
        "not-evaluable",
    ),
];

/// Every curated equivalence class must reach **one** normalized string.
///
/// **The defect this closes.** Every other guard in this file asks whether ferro
/// *changed* a string. None asks whether two rows are the same variant, so the
/// spec's own non-confluent pair — `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` and
/// the `c.[850_869del;…]` split the spec names as its "alternative description"
/// — sat in the fixture as two `preserved` rows and passed. Two stable fixed
/// points for one variant is exactly #1235's defect, and the row-level status
/// bucket reported it as full conformance.
///
/// **Why the classes are curated.** Deciding that two descriptions denote one
/// variant means applying both to a reference sequence, which 50 `needs-reference`
/// rows cannot do. So the declarations live beside the other hand-curated spec
/// judgements, in `hgvs_spec_normalization_overrides.json`, and the generator
/// fails the build if a class names a row that no longer exists or cites spec
/// text that has moved.
///
/// **Why this is not circular.** The class membership comes from the spec text,
/// not from ferro, and the outputs are re-derived live here rather than read off
/// the fixture — the same shape as [`ferro_produces_the_form_the_spec_states`].
/// The fixture supplies only the *inputs* (which rows, under which accession).
#[test]
fn spec_equivalence_classes_converge() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    let by_input: std::collections::HashMap<&str, &Row> =
        fx.rows.iter().map(|r| (r.input.as_str(), r)).collect();
    let normalizer = Normalizer::new(MockProvider::new());

    // Verdict per class, plus the detail a failure needs to be actionable.
    let mut observed: std::collections::BTreeMap<&str, &str> = std::collections::BTreeMap::new();
    let mut detail: std::collections::BTreeMap<&str, String> = std::collections::BTreeMap::new();

    for class in &fx.equivalence_classes {
        let mut outputs: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
        let mut evaluable = 0usize;
        let mut lines = Vec::new();
        for member in &class.members {
            let row = by_input.get(member.input.as_str()).unwrap_or_else(|| {
                panic!(
                    "class {:?} names {:?}, which is not a fixture row",
                    class.id, member.input
                )
            });
            let (normalized, _) = observe(&normalizer, &member.target);
            // A member contributes only when it can be evaluated hermetically.
            // Reference-dependent rows and rows ferro refuses are neither
            // evidence of confluence nor of its absence, so they are dropped
            // rather than compared — otherwise an error string would count as a
            // distinct "representation" and every such class would fail for the
            // wrong reason.
            let skip = row.status == "needs-reference"
                || row.requires_reference == Some(true)
                || normalized.starts_with("parse error")
                || normalized.starts_with("normalize error");
            if skip {
                lines.push(format!("      {} -> (skipped) {normalized}", member.target));
                continue;
            }
            evaluable += 1;
            outputs.insert(normalized.clone());
            lines.push(format!("      {} -> {normalized}", member.target));
        }

        let verdict = if evaluable < 2 {
            "not-evaluable"
        } else if outputs.len() == 1 {
            "confluent"
        } else {
            "non-confluent"
        };
        observed.insert(class.id.as_str(), verdict);
        detail.insert(class.id.as_str(), lines.join("\n"));
    }

    let pinned: std::collections::BTreeMap<&str, &str> =
        EQUIVALENCE_CLASS_VERDICTS.iter().copied().collect();
    assert_eq!(
        pinned.len(),
        EQUIVALENCE_CLASS_VERDICTS.len(),
        "EQUIVALENCE_CLASS_VERDICTS has a duplicate id"
    );

    let mut mismatches = Vec::new();
    for (id, verdict) in &observed {
        match pinned.get(id) {
            None => mismatches.push(format!(
                "  {id}: new class, verdict {verdict} — add it to EQUIVALENCE_CLASS_VERDICTS\n{}",
                detail[id]
            )),
            Some(expected) if expected != verdict => mismatches.push(format!(
                "  {id}: pinned {expected}, observed {verdict}\n{}",
                detail[id]
            )),
            Some(_) => {}
        }
    }
    for id in pinned.keys() {
        if !observed.contains_key(id) {
            mismatches.push(format!(
                "  {id}: pinned but absent from the fixture — a curated class was deleted"
            ));
        }
    }

    eprintln!(
        "spec equivalence classes: {} checked, {} non-confluent",
        observed.len(),
        observed.values().filter(|v| **v == "non-confluent").count()
    );
    assert!(
        !observed.is_empty(),
        "spec_equivalence_classes_converge exercised no classes — the fixture carries none, \
         which means the curated declarations were lost, not that they all pass"
    );
    assert!(
        mismatches.is_empty(),
        "{} equivalence-class verdict(s) moved.\n\n{}\n\n\
         A class is a curated statement that its members are the SAME VARIANT, so more than one \
         distinct output is non-confluence — two stable fixed points for one variant (#1235).\n\n\
         `pinned confluent, observed non-confluent` or a new non-confluent class is a REGRESSION: \
         a variant that used to have one representation now has two.\n\
         `pinned non-confluent, observed confluent` is an IMPROVEMENT — but it is also a \
         representation change for anyone storing the losing form, so release the pin in the same \
         PR that says which form moved and how many corpus rows it touches.\n\
         `observed not-evaluable` means ferro stopped evaluating a member; the class has disarmed \
         itself and the reason is the thing to look at.",
        mismatches.len(),
        mismatches.join("\n"),
    );
}

/// Pinned id and status for every ruling record.
///
/// A ruling records how the project reads the spec where the spec conflicts with
/// itself. `undecided` means the project has **not** ruled — pinning the status
/// is what stops one being upgraded to `decided` without a reviewable diff, and
/// stops an inconvenient record being deleted.
const RULING_STATUSES: &[(&str, &str)] = &[
    // `delins.md:17` ("described individually … not as a delins") against
    // `delins.md:47` ("the delins format is recommended"), both reaching the
    // `:44-47` example. **Decided for `:47` by operator ruling (2026-08-07),
    // and SCOPED** — the scope is part of the ruling. It settles only the shape
    // `:44-47` describes: a *minimal* single `delins` split because payload
    // bases coincide with reference bases. It is not a licence to merge across
    // separation ≥2 generally; where the separation comes from anything else,
    // `general.md:34` still governs. Unscoped it would have reached ~2,112 of
    // 4,732 violations measured by a later corpus-wide audit — fifteen times
    // the 208 rows the argument was made on. Ratifies shipped behaviour
    // (v0.12.0 emits the unsplit form on 208/208), so it moves no row. Note
    // both clauses are *lowercase* prose — the record's former "equal RFC 2119
    // strength" argument was withdrawn, because a census of
    // `docs/recommendations/` finds an uppercase keyword in exactly one place
    // outside `style.md` (`RNA/adjoined_transcript.md:21`).
    ("delins-merge-vs-individual-gap-two-or-more", "decided"),
    // `general.md:34` against `general.md:56` on `c.76_83inv`, whose reverse
    // complement coincides with the reference at its 4 interior columns.
    // Renamed from `inversion-vs-two-substitutions-76-83`: the competing
    // members are two `delins` (`c.[76_77delinsTG;82_83delinsTT]`), and `:56`
    // ranks substitution above inversion while not listing `delins` at all, so
    // the old name carried #1230's answer onto a case `:56` does not reach.
    // **Decided for the inversion by operator ruling (2026-08-07)**, governing
    // `inversion.md:5` — an inversion is a whole-span property, and the spec's
    // own worked example (`inversion.md:33-34`, `NM_004006.2:c.4145_4160inv`)
    // has two three-base unchanged interior runs, so `general.md:34` does not
    // decompose one. `:56` cannot rank the alternative because `delins` is
    // absent from its list, which is what separates this from #1230 — whose
    // competing members really are substitutions, and whose guard stays green.
    ("inversion-vs-two-delins-76-83", "decided"),
    // Which authority ranks after the spec. Decided 2026-08-07: spec, then
    // confluence, then re-derivation from the sequence, then disclosure, with
    // stability only as a last-resort tiebreaker. Mutalyzer does not appear —
    // the filer himself calls it not spec-compliant, and its answers are
    // predicted better by a 2014 description-length weight model than by any
    // separation rule.
    ("adjudication-precedence-order", "decided"),
    // Which of two legal descriptions of one variant ships. Decided 2026-08-07:
    // re-derive from the resulting sequence, subject to the spec's explicit
    // tie-breaks. The spec's own model (`general.md:157-160`), and the method
    // the filer asked for — though the spec contradicts itself, conditioning
    // some preferences on provenance and allele frequency, which ferro cannot
    // see. See the record.
    ("canonical-form-choice-when-both-legal", "decided"),
    // `delins.md:18` names no edit type; ferro applied it only to
    // sub/unchanged/sub. **Decided by operator ruling (2026-08-10): WIDEN.**
    // The exception applies wherever its stated precondition holds — two
    // variants one nucleotide apart, together affecting one amino acid —
    // regardless of edit type, because edit type is a property of the SPELLING
    // and the precondition is a property of the RESULTING SEQUENCE (rule 3 of
    // the README ruleset; the same argument that decided
    // `separation-is-a-property-of-the-spelling-not-of-the-variant`). The two
    // limits survive: a frameshift pair fails the precondition, and
    // `duplication.md:18` still forbids merging a describable duplication
    // away. The 5.76M-row inertness measurement is kept in the record as a
    // measurement, not as support. #1599 is the implementation.
    ("codon-carve-out-shape-restriction", "decided"),
    // The one the spec does settle: `delins.md:18`'s explicit exception for two
    // variants one nucleotide apart affecting one amino acid.
    ("delins-codon-carve-out-gap-one", "decided"),
    // `DNA/duplication.md:26` names `c.3921dup` as *the* description "and not
    // `c.3922dup`", and `:60` and `:148` say the same thing independently —
    // `:148` calling the far-side position "the wrong nucleotide, in the wrong
    // exon", which refutes the reading that `c.3922dup` is a different
    // legitimate genomic event. **Decided by operator ruling (2026-08-10):
    // CONVERGE.** The canonical position is the most 3' position that does not
    // cross an exon/exon junction, reached from EITHER side — `general.md:44`
    // is directional, and ferro leaving both spellings alone is that exception
    // half-applied. Measured: both are fixed points, projecting 2,790 bp
    // apart, and cdot puts the junction exactly between `c.3921` (last base of
    // exon 27) and `c.3922` (first base of exon 28). This MOVES OUTPUT; the
    // implementation is #1621 and the guard for it is `#[ignore]`d in
    // `spec_worked_examples.rs`.
    ("exon-junction-dup-converge-from-the-far-side", "decided"),
    // `RNA/repeated.md:22` marks range-plus-unit invalid as redundant while
    // `:27` publishes exactly that shape as valid, five lines later. Ferro
    // answers both ways depending on input-hygiene mode. Upstream's conflict to
    // settle (#466), not ferro's.
    ("rna-repeat-range-plus-unit-redundancy", "undecided"),
    // A normalized allele must not leave two members flush against each other
    // when both consume reference bases. `substitution.md:32` marks exactly that
    // spelling `class="invalid"` for `LRG_199t1:c.79_80delinsTT`, and
    // `general.md:34` does not compete — it governs members separated by one or
    // more nucleotides, not by none. Ferro CONFORMS as of #1537, which closed
    // #1524; `cis_confluence_adjudication.rs` pins the spec's form as a
    // regression guard, not the deviation the record was first written against.
    (
        "delins-adjacent-members-when-both-consume-reference",
        "decided",
    ),
    // The separation `general.md:34` keys on is read off a decomposition, and
    // two spellings of one variant can present different ones — demonstrated,
    // not argued, in `cis_confluence_adjudication.rs`. #1537 settled the
    // separation-ZERO half and left this untouched. **Decided by operator
    // ruling (2026-08-10)**: the separation is read off the partition
    // RE-DERIVED FROM THE RESULTING SEQUENCE, never off the input's spelling —
    // rule 3 of the README ruleset — which on the record's genomic case makes
    // `g.[1001009_1001010del;1001013del]` the answer from both spellings. The
    // sequence-first arms (`FERRO_PARTITION=shadow`/`canonical`) already
    // produce it; the shipping path does not, and that gap is #1617, so the
    // pinned assertions still show today's divergent `live` output alongside an
    // `#[ignore]`d guard for the decided target. Which partition is canonical
    // in general remains `canonical-form-choice-when-both-legal`.
    (
        "separation-is-a-property-of-the-spelling-not-of-the-variant",
        "decided",
    ),
    // `general.md:58` (removing and replacing part of the same reference
    // sequence "are not allowed") against `DNA/complex.md:130` ("'::' is used to
    // indicate the join, instead of ';'") on whether the self-cancelling check
    // reaches a ring's `::`-joined segments. **Decided for `complex.md:130`**:
    // `:58` is a four-space sub-bullet of `general.md:56`'s *prioritisation*
    // bullet — sibling to `:57`, which is explicitly a prioritisation rule — so
    // it applies only "when a description is possible according to several
    // types" and works by redirecting to the preferred one. `complex.md:5`
    // defines complex as what *cannot* be described as a basic type, so no
    // competing single-type description exists for a ring and there is nothing
    // to redirect to. Ratifies shipped behaviour, so it moves no row.
    //
    // Scope: this decides only that `:58` does not cross `::`. It explicitly
    // does NOT bless the three malformed rings ferro currently accepts — those
    // are ring *well-formedness* defects (`complex.md:39`/`:51`/`:53`/`:55`)
    // needing their own rule. The negative is pinned by
    // `issue_1578_followup_self_cancelling_rings.rs` so nobody closes the escape
    // the cheap way by pointing `:58` at rings.
    ("self-cancelling-across-ring-junctions", "decided"),
    // `standards.md:36`/`:37` list `X` and `-` in the DNA symbol table while
    // `:39` daggers both as "used in alignment only". Decided for `:39` — the
    // dagger sits inside the symbol cell, so `:39` is the row's own annotation
    // rather than a competing clause. The strength grading is moot either way:
    // `general.md:48` admits only IUPAC-IUBMB symbols and `X` is not one, so the
    // shape fails the grammar half of rank-1 validity on its own. Ferro already
    // refuses `-` at parse and accepts `X` in both modes; ratifying the `-`
    // treatment moves no legal output.
    ("alignment-only-symbol-in-a-description", "decided"),
    // `checklist.md:20` conditions an intronic `c.` position on a genomic
    // wrapper; `checklist.md:45` then glosses a bare `c.12-14del` as an intronic
    // deletion, twenty-five lines later. Decided for `:20` AS A CONDITIONAL
    // clause: strict refuses (W4007), lenient accepts. Ratifies shipped
    // behaviour in both modes. Does not excuse the 371 bare-`NM_` intronic
    // descriptions ferro's own junction clamp EMITS.
    ("bare-transcript-intronic-position", "decided"),
    // Whether `general.md:58`'s member-vs-member prohibition reaches the class of
    // alleles whose members claim intersecting territory, or only the `del`+`dup`
    // pair it names. Decided for `DNA/alleles.md:5` — the definition, not `:58`,
    // is what reaches nested and coincident-insertion geometries; `:58`'s ground
    // ("replacing it with part of the same sequence") literally describes only its
    // own example. `general.md:56` is cited to record that it does NOT reach the
    // question: it ranks descriptions of ONE span, and citing it against a
    // multi-member allele is this repository's recorded cautionary error.
    ("conflicting-member-geometry-refusal-scope", "decided"),
    // WHERE an absolute prohibition is refused — at parse (unconditional) or at
    // strict-mode normalize (opt-outable). The spec addresses descriptions, not
    // stages, so there is nothing to defer to. Undecided, and recorded because
    // ferro answers it three ways for clauses of identical strength: within
    // `checklist.md:16` alone, `g.*10del` is refused at parse while `g.266+2del`
    // is refused nowhere.
    ("absolute-prohibition-enforcement-stage", "undecided"),
    // Must a ring's first `::` segment start at `pter` and its last end at
    // `qter`? **Undecided on purpose.** A ring loses both telomeres, and both
    // published ring shapes are anchored (`complex.md:127`, `:161`) — but that is
    // biology plus two examples, and no clause states it. `complex.md:28` only
    // defines what the markers mean; `:13` warns that a complex description can
    // be "literally correct" yet meaningless, which cuts against inventing a
    // structural requirement. Its sibling junction rule *was* decided, on
    // `:60-64`'s record that the committee withdrew `::` as a general join
    // operator — a materially stronger argument, which is exactly why this one is
    // not decided by analogy to it. Ferro accepts an unanchored ring today.
    ("ring-telomere-anchoring", "undecided"),
];

/// Every case where a preference the spec *states* was overridden, because the
/// spec supplied nothing that decides between the candidates satisfying it.
///
/// **This is a census, and its whole job is to make growth loud.** The canonical
/// ruleset is in `README.md`; `adjudication-precedence-order` points at it and
/// holds this register. Under that ruleset, producing the recommended form is
/// rule 2 — so an override of rule 2 is the rare, deliberate act this register
/// exists to count. Left uncounted, overrides would accumulate into a second,
/// undeclared ruleset, which is precisely the drift the pointer record guards
/// against.
///
/// **Do not file ordinary spec silence here.** Where the spec says *nothing at
/// all* the ruleset's rule 5 applies — file upstream first, cite it, ship a
/// provisional choice — and that zone is normal and not rare. This register is
/// only for a *stated* preference being overridden, of which there is one.
///
/// Adding an entry therefore fails [`precedence_escalations_are_censused`] until
/// the count below is updated deliberately. That friction is the point; it is the
/// same reason `KNOWN_DIVERGENT_INPUTS` is pinned rather than appended to freely.
///
/// Each entry must additionally carry, per the record: the clause that expresses
/// the preference, the reason no rule follows from it, and a unit test locking the
/// choice with an exact `file.md:line` citation. The second field names the ruling
/// record that governs the choice.
const PRECEDENCE_ESCALATIONS: &[(&str, &str)] = &[
    // E1. `general.md:56` ranks substitution, deletion, inversion, duplication
    // and insertion — and `delins` appears nowhere in it. So where the competing
    // forms are all `delins`, `:56` ranks nothing, while `general.md:34` gives a
    // separation test but no preference between two forms that both satisfy it.
    // The spec has stated a preference and supplied nothing that decides, which
    // is what puts this here rather than under ordinary spec silence.
    //
    // Recorded rather than decided. `canonical-form-choice-when-both-legal` is
    // `decided` and answers the general form of the question — derive from the
    // resulting sequence, then apply every explicit spec tie-break — so what E1
    // covers is the residue: two all-`delins` candidates that derivation itself
    // left tied and the spec ranks neither. That is a *typing* question over
    // candidates the derivation produced, which is exactly `:56`'s domain and
    // exactly where its `delins` gap bites. It is never a licence to keep the
    // partition the input arrived in.
    (
        "delins-is-unranked-by-general-56",
        "canonical-form-choice-when-both-legal",
    ),
];

/// The escalation register is exactly as long as the ruling says it is.
///
/// A bare length check on purpose. Asserting the *contents* would duplicate the
/// per-entry tests the record already requires; what is not otherwise guarded is
/// the count, because an entry can be appended without any existing test failing.
#[test]
fn precedence_escalations_are_censused() {
    const EXPECTED: usize = 1;
    assert_eq!(
        PRECEDENCE_ESCALATIONS.len(),
        EXPECTED,
        "the escalation register changed. Adding one is a decision for the operator, \
         not a contributor: `adjudication-precedence-order` requires every override of \
         a stated spec preference to be surfaced, documented in the register, and locked \
         by a unit test citing the clause. Update this count only alongside those three. \
         If what you are adding is ordinary spec silence, it does not belong here — that \
         is the ruleset's rule 5, and it is filed upstream instead."
    );
    let ids: std::collections::BTreeSet<&str> =
        PRECEDENCE_ESCALATIONS.iter().map(|(id, _)| *id).collect();
    assert_eq!(
        ids.len(),
        PRECEDENCE_ESCALATIONS.len(),
        "duplicate escalation id: a duplicate would let one entry shadow another \
         while the count still matched"
    );
}

/// Ruling records must stay well-formed, and `undecided` must stay undecided.
///
/// The generator enforces the same shape at build time (and additionally checks
/// that every citation still quotes the spec checkout). This is the committed
/// half: the generated fixture is gitignored and regenerated from the tree, so
/// without a pin here a record could be softened, re-statused or dropped and
/// nothing would notice.
#[test]
fn ruling_records_are_intact() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    let observed: std::collections::BTreeMap<&str, &str> = fx
        .rulings
        .iter()
        .map(|r| (r.id.as_str(), r.status.as_str()))
        .collect();
    let pinned: std::collections::BTreeMap<&str, &str> = RULING_STATUSES.iter().copied().collect();
    // Collecting into a map is last-wins, so a duplicated id silently drops one
    // pin and `observed == pinned` can still hold on a malformed table — the id
    // that survived is checked, the one it shadowed is not (#1530). Its sibling
    // `EQUIVALENCE_CLASS_VERDICTS` already carries this; the generator checks
    // duplicates on both of its own tables.
    //
    // Both sides need the guard, not just the pinned one. The hazard above is
    // symmetric, but only `RULING_STATUSES` was checked, so a duplicated id in
    // the *fixture* stayed invisible here: the shadowed record's status is never
    // compared against anything, which is the one property this assertion pair
    // exists to enforce. (The loop below iterates `fx.rulings` directly, so a
    // shadowed record still gets its clause/rationale/governing checks — it is
    // specifically the id-to-status pin that could be evaded.)
    assert_eq!(
        pinned.len(),
        RULING_STATUSES.len(),
        "RULING_STATUSES has a duplicate id"
    );
    assert_eq!(
        observed.len(),
        fx.rulings.len(),
        "the ruling fixture has a duplicate id; one record's status is shadowed and \
         therefore unpinned"
    );
    assert_eq!(
        observed, pinned,
        "the ruling records changed. Adding one is fine — pin it here. Changing an `undecided` \
         to `decided` means the project took a position, which must be a deliberate, reviewable \
         edit, not a side effect."
    );

    for ruling in &fx.rulings {
        assert!(
            !ruling.clauses.is_empty(),
            "ruling {:?} cites no clause",
            ruling.id
        );
        assert!(
            !ruling.rationale.trim().is_empty(),
            "ruling {:?} has no rationale",
            ruling.id
        );
        let cited: std::collections::BTreeSet<&str> =
            ruling.clauses.iter().map(|c| c.clause.as_str()).collect();
        match ruling.status.as_str() {
            "decided" => {
                let governing = ruling.governing.as_deref().unwrap_or_else(|| {
                    panic!(
                        "ruling {:?} is decided but names no governing clause",
                        ruling.id
                    )
                });
                assert!(
                    cited.contains(governing),
                    "ruling {:?} holds {governing:?} to govern, but does not cite it",
                    ruling.id
                );
                for deviated in &ruling.deviates_from {
                    assert!(
                        cited.contains(deviated.as_str()),
                        "ruling {:?} deviates from {deviated:?}, but does not cite it",
                        ruling.id
                    );
                }
            }
            // The load-bearing half. An undecided record states a conflict and
            // stops; if it names a governing clause it has smuggled in a ruling
            // nobody made, which is precisely what these records exist to
            // prevent being invented.
            "undecided" => {
                assert!(
                    ruling.governing.is_none() && ruling.deviates_from.is_empty(),
                    "ruling {:?} is undecided but names a governing or deviated-from clause — \
                     an undecided record must not imply a ruling",
                    ruling.id
                );
                assert!(
                    ruling.clauses.len() >= 2,
                    "ruling {:?} is undecided but names only one clause; an unsettled question \
                     needs both sides on the record",
                    ruling.id
                );
            }
            other => panic!("ruling {:?} has unknown status {other:?}", ruling.id),
        }
    }
}
