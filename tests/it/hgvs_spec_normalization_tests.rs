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

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Debug, Deserialize)]
struct Fixture {
    rows: Vec<Row>,
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
    ("diverges", 36),
    ("false-acceptance", 52),
    ("needs-reference", 50),
    ("preserved", 697),
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
/// fixture's 36 `diverges` rows, which the status census counts but — before
/// #1272 — nothing checked the *values* of.
const KNOWN_DIVERGENT_INPUTS: &[&str] = &[
    // 1. `*` is canonicalised to `Ter` in protein descriptions. Deliberate
    //    render policy: both spellings are legal on input (#1114 made the `*`
    //    form parse), and ferro emits one of them.
    "NP_003997.1:p.W24*",
    "NP_003997.2:p.(Asn444Lysfs*15)",
    "NP_003997.2:p.(His321Leufs*3)",
    "p.(*110Glnext*17)",
    "p.(Tyr4*)",
    "p.*110Glnext*17",
    "p.*327Argext*?",
    "p.Arg1459*",
    "p.Arg456Glyfs*17",
    "p.Arg97Profs*23",
    "p.Gln151Thrfs*9",
    "p.Ile327Argfs*?",
    "p.Trp24*",
    "p.Trp26*",
    "p.Trp41*",
    "p.Tyr4*",
    "p.[(Lys79*;Lys79Asn)]",
    "p.[(Ser868Argfs*2,Ser868=)]",
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
        let known_divergent = KNOWN_DIVERGENT_INPUTS.contains(&row.input.as_str());

        let conformant = match row.spec_expected.as_deref() {
            // The spec states a form: ferro must render it.
            Some(expected) => observed == expected,
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
