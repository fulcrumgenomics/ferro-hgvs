//! Harvested reproducers, run hermetically against real reference bases.
//!
//! # The question
//!
//! An audit of this repository's issues, PR descriptions and campaign notes
//! found roughly 150 concrete HGVS rows that carry a behavioural claim — "this
//! input produces that output, and here is why it is right or wrong" — and that
//! **nothing running in CI can fail on**. Each had been measured once, written
//! into prose, and left with no executable home. Can the verified subset be
//! given one that blocks a PR?
//!
//! Three existing homes each decline, for structural reasons:
//!
//! - the harvested spec fixture (`hgvs_spec_normalization_tests.rs`) normalizes
//!   through `MockProvider::new()`, so no sequence-dependent rule can fire, and
//!   its two committed guards skip the `requires-reference` and
//!   `false-acceptance` statuses — which is exactly where the spec's own
//!   published repair targets sit;
//! - the manifest-backed conformance axes skip without `FERRO_MANIFEST`, which
//!   in CI is always;
//! - the synthetic corpus behind `examples/dump_normalized_corpus.rs` cannot
//!   *build* several of these shapes at all (#1456 member geometry, #1460 scale,
//!   #1478 transcript geometry).
//!
//! So this file is the fourth: a curated corpus, a committed slice of the
//! prepared reference, and no skip path.
//!
//! # The rulings this file rests on
//!
//! Each row carries its own `file:line` citation into `assets/hgvs-nomenclature`
//! and the driver verifies the quote verbatim, so the authority is in the
//! fixture rather than restated here. The four that shape the file's structure:
//!
//! - `DNA/inversion.md:5` — "**more than one nucleotide** replacing the original
//!   sequence is the reverse complement of the original sequence" — and
//!   `DNA/delins.md:5` — "**and which is not** a substitution or inversion".
//!   A definitional-exclusion **pair**, not a ranking: it settles
//!   `NM_004006.2:c.76_83inv` without any appeal to prioritisation.
//! - `DNA/delins.md:16` — "changes involving two or more consecutive nucleotides
//!   are described as deletion/insertion (delins) variants" — the separation-0
//!   rule the corpus-wide structural check enforces.
//! - `general.md:41` — the 3'rule — and `DNA/repeated.md:21-24` /
//!   `RNA/repeated.md:24-25`, whose worked answers the spec prints outright.
//! - `background/basics.md:38` — "designed to be **stable**, **meaningful**,
//!   **memorable**, and **unequivocal**" — the authority for the *stability*
//!   pins, which are explicitly **not** spec-correctness claims. Stability is
//!   the value the spec leads with; it is not authority for any particular
//!   spelling, and the rows that cite it say so.
//!
//! # The mechanism this file must not fall for
//!
//! A window fixture captured from a normalize pass is exactly as wide as the
//! reads that pass made. A later, wider-reading pass runs off the end, the
//! provider errors, and `src/normalize/mod.rs` converts that into
//! `Ok((input-unchanged, vec![]))` — success, no warnings, input returned. For a
//! lone-member row that is byte-identical to a correct answer, so a guard built
//! naively on such a fixture is **green while the defect is fully present**.
//! That has happened here: a guard passed on three of five rows with the defect
//! in full.
//!
//! Two mitigations, both mandatory and both asserted below:
//!
//! 1. the generator's `serve_transcripts_whole` widens every transcript's window
//!    to the entire sequence, so there is no "outside the window" for a
//!    transcript; [`the_committed_windows_serve_every_transcript_whole`] pins
//!    that hermetically, so a regeneration that loses the widening fails loudly;
//! 2. every pass runs through
//!    [`AuditedProvider`](ferro_hgvs::conformance::audited_provider::AuditedProvider),
//!    and the zero-failed-reads assertion is evaluated **before** any verdict.
//!    Ordering is the point: "the fixture could not answer" must always outrank
//!    "ferro answered differently". A positive floor on *served* reads sits
//!    beside the zero ceiling, because reading nothing also yields zero
//!    failures.
//!
//! # Status
//!
//! Green: 29 of 34 rows produce their recorded answer, the two green confluence
//! classes converge, every answer is a fixed point, and no emitted output puts
//! two members on consecutive nucleotides.
//!
//! Red, `#[ignore]`d, and expected to stay red until fixed:
//!
//! - **#1541** — `NM_000500.9:c.[710T>A;713T>A]`, `c.710_713inv` and
//!   `c.710_713delinsACGA` do not converge. The target form is **decided, not
//!   open**: `tests/it/issue_1517_inv_priority_over_delins.rs` pins both
//!   directions of `general.md:56`, and because this block's competing members
//!   are *substitutions* (ranked 1, above inversion at 3) the preferred form is
//!   the split, `c.[710T>A;713T>A]` — the same answer as #1230. The blocker is
//!   the **axis gate**, not the typing rule: `is_splittable_single_member`
//!   (`src/normalize/mod.rs:2005`) matches only `HV::Genome`/`HV::Mt`, so a lone
//!   `c.` member never reaches the sequence-first pass, and #1230's control
//!   passes only because it is written on the `g.` axis. PR #1484 widens that
//!   gate; **un-ignore [`the_four_mer_inversion_pair_converges`] in the change
//!   that merges it.** The control for the other direction is
//!   [`the_inversion_preference_control_still_holds`], which is green.
//! - **#1542** — `NC_000017.11:g.80110044_80110047delinsGTTGG` emits members at
//!   separation 0 in a `dup`-shaped split, surviving #1537. The separation
//!   reading is itself part of the open question: the adjacency depends on
//!   treating a `dup` as occupying the position it duplicates, and under a
//!   zero-width (interbase) reading the pair is not adjacent at all.
//!
//! Neither red test may be weakened to green. A green re-classification is what
//! [`the_corpus_census_is_unchanged`] exists to catch.
//!
//! # One measurement in the brief did not reproduce
//!
//! `NM_000500.9:c.710_713delinsACGA` was reported as a fixed point on `main`
//! alongside its two siblings. It is not: on `main` it normalizes to
//! `NM_000500.9:c.710_713inv`. #1541's verdict is unchanged — the class still
//! yields two normal forms — but the split is **2 against 1**, not three
//! mutually disagreeing forms, and the row's `note` records the correction.
//!
//! # Regenerating the slice
//!
//! `cargo run --features dev --example extract_case_harvest_windows -- \
//!  --manifest <manifest>`, and never by hand.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use ferro_hgvs::conformance::audited_provider::AuditedProvider;
use ferro_hgvs::conformance::case_harvest::{Case, Fixture, CASES_PATH, SPEC_DIR, WINDOWS_PATH};
use ferro_hgvs::conformance::reference_window::{WindowFixture, WindowProvider};
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::{NormalizeConfig, Normalizer, ReferenceProvider, ShuffleDirection};

/// The command that rebuilds the committed slice, quoted in every failure whose
/// remedy is a regeneration.
const REGENERATE: &str = "cargo run --features dev --example extract_case_harvest_windows -- \
                          --manifest <manifest>";

// ---------------------------------------------------------------------------
// Fixture plumbing
// ---------------------------------------------------------------------------

/// Anchored on `CARGO_MANIFEST_DIR` so the paths resolve regardless of the
/// test's working directory.
fn repo_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

fn cases() -> Fixture {
    Fixture::from_json_path(&repo_path(CASES_PATH)).expect("load the harvested-case corpus")
}

fn window_fixture() -> WindowFixture {
    WindowFixture::from_json_path(&repo_path(WINDOWS_PATH))
        .expect("load the committed harvested-case reference slice")
}

/// The hermetic provider, wrapped so every failed reference read is recorded.
/// No manifest, no environment variable, no skip path.
fn audited_provider() -> AuditedProvider<WindowProvider> {
    AuditedProvider::new(window_fixture().to_provider())
}

/// Run one input through the shipped `--error-mode lenient` path, mirroring
/// `src/bin/ferro.rs`: the error configuration reaches both the preprocessor and
/// the normalizer (#1181/#1197). This is the mode every recorded answer in
/// `cases.json` was measured under.
fn run<P: ReferenceProvider + Clone>(provider: &P, input: &str) -> Result<String, String> {
    let config = ErrorConfig::lenient();
    let parsed = parse_hgvs_with_config(input, config.clone())
        .map_err(|e| format!("preprocess error: {e}"))?;
    let normalize_config = NormalizeConfig::for_entry_point(ShuffleDirection::ThreePrime, config);
    Normalizer::with_config(provider.clone(), normalize_config)
        .normalize(&parsed.result)
        .map(|n| n.to_string())
        .map_err(|e| format!("normalize error: {e}"))
}

/// The audit assertion, in the one order that is safe.
///
/// A read the committed slice cannot serve makes normalization bail and return
/// its input unchanged, which is indistinguishable at the output from a clean
/// answer — so this must be checked *before* any verdict, or a fixture defect
/// gets reported as a conformance result. The floor on served reads is the
/// other half: zero failures is trivially satisfiable by reading nothing.
fn assert_the_slice_answered_everything<P>(provider: &AuditedProvider<P>, floor: usize) {
    let failures = provider.failures();
    assert!(
        failures.is_empty(),
        "the hermetic slice failed {} reference read(s) during this pass, so any verdict below is \
         about a normalization that ran short of reference data, not about ferro's answer. \
         Regenerate with `{REGENERATE}`:\n  {}",
        failures.len(),
        failures.join("\n  ")
    );
    assert!(
        provider.successful_reads() >= floor,
        "the pass served only {} reference read(s), fewer than the {floor} expected — a zero \
         failure count means nothing when nothing was read",
        provider.successful_reads()
    );
}

// ---------------------------------------------------------------------------
// Member geometry (delins.md:16)
// ---------------------------------------------------------------------------

/// Every member boundary in `description`, as `(end of one, start of the next)`
/// on the HGVS axis. `None` when a member does not begin with a plain integer
/// coordinate (an intronic offset, a UTR position), which this crude reader
/// deliberately declines rather than mis-parses.
///
/// It reads the string rather than the parsed variant, on purpose: what
/// `delins.md:16` constrains is the *description*, and re-deriving the spans
/// from the code that produced them would check that code against itself.
fn member_bounds(description: &str) -> Option<Vec<(u64, u64)>> {
    let (_, body) = description.split_once(':')?;
    let body = body.get(2..)?; // strip the `g.` / `c.` / `n.` / `r.` / `m.` prefix
    let inner = body.strip_prefix('[')?.strip_suffix(']')?;
    let leading_number = |s: &str| -> Option<u64> {
        let digits: String = s.chars().take_while(char::is_ascii_digit).collect();
        digits.parse().ok()
    };
    let mut spans = Vec::new();
    for member in inner.split(';') {
        let start = leading_number(member)?;
        let end = match member.split_once('_') {
            Some((_, rest)) => leading_number(rest)?,
            None => start,
        };
        spans.push((start, end));
    }
    Some(spans.windows(2).map(|w| (w[0].1, w[1].0)).collect())
}

/// The members of `description` that sit on consecutive nucleotides.
fn adjacent_members(description: &str) -> Vec<(u64, u64)> {
    member_bounds(description)
        .unwrap_or_default()
        .into_iter()
        .filter(|(previous_end, next_start)| *next_start <= previous_end + 1)
        .collect()
}

// ---------------------------------------------------------------------------
// Reading the reference out of the committed slice
// ---------------------------------------------------------------------------

/// Bases `c.start` through `c.end` inclusive, read from the stored transcript.
///
/// Premise tests read the reference through this rather than restating it, so a
/// defect fix cannot be credited to a change of fixture.
fn coding_bases(fixture: &WindowFixture, transcript_id: &str, start: usize, end: usize) -> String {
    let transcript = fixture
        .transcripts
        .iter()
        .find(|t| t.id == transcript_id)
        .unwrap_or_else(|| panic!("{transcript_id} is not in the committed slice"));
    let sequence = transcript
        .sequence
        .as_deref()
        .unwrap_or_else(|| panic!("{transcript_id} carries no bases"));
    let cds_start = transcript
        .cds_start
        .unwrap_or_else(|| panic!("{transcript_id} carries no CDS start"))
        as usize;
    // `c.n` is transcript position `cds_start + n - 1`, 1-based; the index is
    // one less again.
    sequence[cds_start + start - 2..cds_start + end - 1].to_string()
}

/// Bases `g.start` through `g.end` inclusive, read from the committed windows.
fn genomic_bases(provider: &WindowProvider, contig: &str, start: u64, end: u64) -> String {
    provider
        .get_genomic_sequence(contig, start - 1, end)
        .unwrap_or_else(|e| panic!("{contig}:g.{start}_{end} is not in the committed slice: {e}"))
}

/// Watson-Crick reverse complement, uppercase.
fn reverse_complement(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// The gate
// ---------------------------------------------------------------------------

/// Every row produces its recorded answer.
///
/// Reported as one batch rather than row-by-row: when a normalizer change moves
/// these, *which* rows moved together is the diagnosis, and one assert per row
/// would surface only the first.
///
/// Rows whose `expected` is `null` are not compared here. That is not an
/// omission — see [`the_confluence_classes_converge`]: where the governing
/// ruling is open, asserting that two spellings agree is a correctness claim the
/// spec supports, while asserting *which* form they agree on would freeze a
/// representation nobody has adjudicated.
///
/// The two red rows' answers are pinned too, at today's **wrong** output. So
/// fixing #1542 turns this test red as well as turning its `#[ignore]`d guard
/// green — which is correct: a red row's answer moving is a representation
/// change, and it should require re-blessing the fixture rather than passing
/// unnoticed.
#[test]
fn every_row_produces_its_recorded_answer() {
    let cases = cases();
    let provider = audited_provider();

    let mut observed = Vec::new();
    for case in &cases.cases {
        observed.push(run(&provider, &case.input));
    }
    assert_the_slice_answered_everything(&provider, cases.cases.len());

    let mut failures = Vec::new();
    let mut compared = 0usize;
    for (case, actual) in cases.cases.iter().zip(&observed) {
        let Some(expected) = case.expected.as_deref() else {
            continue;
        };
        compared += 1;
        if actual.as_deref() != Ok(expected) {
            failures.push(format!(
                "  {}\n    expected: {expected}\n    actual:   {}\n    {}\n    {}",
                case.input,
                actual.as_deref().unwrap_or_else(|e| e),
                case.citation.clause,
                case.note
            ));
        }
    }
    assert_eq!(
        compared,
        cases.cases.len() - unpinned_inputs(&cases).len(),
        "the set of pinned rows changed size"
    );
    assert!(
        failures.is_empty(),
        "{} of {compared} pinned rows moved:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// What one confluence class produced: each member's input paired with its
/// normalized answer, or the error string that stood in for one.
type ClassOutcomes<'a> = Vec<(&'a str, Result<String, String>)>;

/// Two spellings of one variant normalize to one string.
///
/// This is the confluence half, and it is a genuine correctness claim: the rows
/// in a class were each verified sequence-equivalent from raw FASTA, so a
/// disagreement is two normal forms for one variant, full stop.
///
/// Only the green classes run here. The red class is #1541, asserted by
/// [`the_four_mer_inversion_pair_converges`], which is `#[ignore]`d and fails.
#[test]
fn the_confluence_classes_converge() {
    let cases = cases();
    let provider = audited_provider();

    let mut outcomes: Vec<(&str, ClassOutcomes<'_>)> = Vec::new();
    for class in cases.confluence_classes() {
        let members = cases.members_of(class);
        if members.iter().any(|case| case.is_red()) {
            continue;
        }
        outcomes.push((
            class,
            members
                .iter()
                .map(|case| (case.input.as_str(), run(&provider, &case.input)))
                .collect(),
        ));
    }
    assert_the_slice_answered_everything(&provider, outcomes.len());

    assert!(
        !outcomes.is_empty(),
        "no green confluence class ran — the corpus cannot report convergence it never measured"
    );
    for (class, members) in &outcomes {
        assert!(
            members.len() >= 2,
            "confluence class {class} has fewer than two members, so it asserts nothing"
        );
        let distinct: BTreeSet<String> = members
            .iter()
            .map(|(_, result)| match result {
                Ok(output) => output.clone(),
                Err(e) => format!("<{e}>"),
            })
            .collect();
        assert_eq!(
            distinct.len(),
            1,
            "confluence class {class} produced {} distinct normal forms for one variant:\n{}",
            distinct.len(),
            members
                .iter()
                .map(|(input, result)| format!("  {input}\n    -> {result:?}"))
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}

/// No emitted output puts two members on consecutive nucleotides
/// (`DNA/delins.md:16`).
///
/// A structural property, checked over the *whole* corpus rather than on the
/// three rows that motivated it, because that is what makes it a guard rather
/// than three more pins. Exactly one row is exempt — the #1542 reproducer — and
/// it is exempt because it *is* the violation this check looks for; the
/// `#[ignore]`d [`the_dup_shaped_split_does_not_touch`] holds it.
#[test]
fn no_emitted_output_puts_two_members_on_consecutive_nucleotides() {
    let cases = cases();
    let provider = audited_provider();

    let mut observed = Vec::new();
    for case in &cases.cases {
        observed.push(run(&provider, &case.input));
    }
    assert_the_slice_answered_everything(&provider, cases.cases.len());

    let mut violations = Vec::new();
    let mut multi_member_outputs = 0usize;
    for (case, actual) in cases.cases.iter().zip(&observed) {
        let Ok(output) = actual else { continue };
        let Some(bounds) = member_bounds(output) else {
            continue;
        };
        if bounds.is_empty() {
            continue;
        }
        if case.allows_adjacent_members() {
            continue;
        }
        multi_member_outputs += 1;
        for (previous_end, next_start) in adjacent_members(output) {
            violations.push(format!(
                "  {} -> {output}\n    members touch at {previous_end}/{next_start} — {}",
                case.input, case.citation.clause
            ));
        }
    }

    // A zero is only meaningful if multi-member outputs were actually examined.
    assert!(
        multi_member_outputs >= 5,
        "only {multi_member_outputs} multi-member output(s) were examined, so the separation check \
         cannot have measured what it claims to"
    );
    assert!(
        violations.is_empty(),
        "{} output(s) violate delins.md:16 (\"changes involving two or more consecutive \
         nucleotides are described as deletion/insertion (delins) variants\"):\n{}",
        violations.len(),
        violations.join("\n")
    );
}

/// Every answer the corpus produces is a fixed point.
///
/// Cheap, and worth having: each output is itself a legal description, so
/// re-normalizing it must be a no-op. A rule that reaches the right answer by
/// one shift too many and one back satisfies every pin above and fails this.
#[test]
fn every_answer_is_a_fixed_point() {
    let cases = cases();
    let provider = audited_provider();

    let mut drifting = Vec::new();
    for case in &cases.cases {
        let Ok(once) = run(&provider, &case.input) else {
            continue;
        };
        let twice = run(&provider, &once);
        if twice.as_deref() != Ok(once.as_str()) {
            drifting.push(format!("  {} -> {once} -> {twice:?}", case.input));
        }
    }
    assert_the_slice_answered_everything(&provider, cases.cases.len());
    assert!(
        drifting.is_empty(),
        "{} answer(s) are not fixed points:\n{}",
        drifting.len(),
        drifting.join("\n")
    );
}

/// Every row cites the spec, and the citation still quotes it.
///
/// A pinned output with no authority behind it is a change detector, not a
/// record: it tells a future reader that something moved, not whether the old
/// answer or the new one was right. Verifying the quote verbatim keeps the
/// citation honest across a spec-submodule bump — a bare line number resolves
/// against any file long enough to have that line.
#[test]
fn every_row_cites_the_spec_verbatim() {
    let spec_dir = repo_path(SPEC_DIR);
    assert!(
        spec_dir.join("docs/recommendations/general.md").is_file(),
        "the vendored HGVS spec checkout at {SPEC_DIR} is empty, so no citation can be verified. \
         Initialise it:\n    git submodule update --init {SPEC_DIR}"
    );
    let cases = cases();
    let mut broken = Vec::new();
    for case in &cases.cases {
        assert!(
            !case.note.trim().is_empty(),
            "{} records no reason for being here",
            case.input
        );
        if let Err(e) = case.citation.verify(&spec_dir) {
            broken.push(format!("  {}: {e}", case.input));
        }
    }
    assert!(
        broken.is_empty(),
        "{} citation(s) no longer resolve against {SPEC_DIR} \
         (submodule at {}):\n{}",
        broken.len(),
        cases.spec_commit,
        broken.join("\n")
    );
}

// ---------------------------------------------------------------------------
// The census
// ---------------------------------------------------------------------------

/// Inputs whose answer is deliberately not pinned.
fn unpinned_inputs(cases: &Fixture) -> Vec<&str> {
    cases
        .cases
        .iter()
        .filter(|case| case.expected.is_none())
        .map(|case| case.input.as_str())
        .collect()
}

/// The corpus is the size and shape it says it is.
///
/// Without this, a regression could be laundered into the fixture by
/// re-classifying a row: pinning a moved answer, marking a newly-failing row as
/// a known defect, or dropping a row that started failing. Every other test here
/// would stay green. So the identities — not merely the counts — of the red set,
/// the unpinned set and the confluence classes are asserted.
#[test]
fn the_corpus_census_is_unchanged() {
    let cases = cases();
    assert_eq!(cases.cases.len(), 34, "the corpus changed size");

    let inputs: BTreeSet<&str> = cases.cases.iter().map(|c| c.input.as_str()).collect();
    assert_eq!(
        inputs.len(),
        cases.cases.len(),
        "the corpus carries a duplicate input"
    );

    let red: Vec<&str> = cases
        .cases
        .iter()
        .filter(|case| case.is_red())
        .map(|case| case.input.as_str())
        .collect();
    assert_eq!(
        red,
        [
            "NM_000500.9:c.[710T>A;713T>A]",
            "NM_000500.9:c.710_713inv",
            "NM_000500.9:c.710_713delinsACGA",
            "NC_000017.11:g.80110044_80110047delinsGTTGG",
        ],
        "the red set changed. Adding a row to it converts a regression into a recorded defect, \
         which is exactly the laundering this assertion exists to block; removing one claims a \
         defect is fixed, which must come with the `#[ignore]` being lifted."
    );

    let issues: BTreeSet<&str> = cases
        .cases
        .iter()
        .filter_map(|case| case.defect.as_ref())
        .map(|defect| defect.issue.as_str())
        .collect();
    assert_eq!(
        issues,
        ["#1541", "#1542"].into_iter().collect::<BTreeSet<_>>(),
        "every red row must name a filed issue, and only these two are open"
    );

    assert_eq!(
        unpinned_inputs(&cases),
        [
            "NM_000143.3:c.[1A>T;3G>A]",
            "NM_000143.3:c.1_3delinsTTA",
            "NM_000500.9:c.[710T>A;713T>A]",
            "NM_000500.9:c.710_713inv",
            "NM_000500.9:c.710_713delinsACGA",
        ],
        "the set of rows that deliberately pin no output changed. Un-pinning a row silences a \
         movement the gate would otherwise catch, so it is an adjudication, not a fixture edit."
    );

    assert_eq!(
        cases.confluence_classes(),
        [
            "dmd-76-83-whole-block-inversion",
            "fh-start-codon-triplet",
            "cyp21a2-710-713-four-mer-inversion",
        ],
        "the confluence classes changed"
    );

    assert_eq!(
        cases.spec_commit, "6f85311989e76ead95d3547828f97ebaa3802e35",
        "the citations were verified against a different spec checkout; re-verify before bumping"
    );
}

// ---------------------------------------------------------------------------
// Fidelity of the committed slice
// ---------------------------------------------------------------------------

/// The committed slice carries every input's reference.
///
/// A missing one would make the gate fail on the audit assertion rather than
/// pass, but naming the regeneration command here turns that into a one-line
/// diagnosis.
#[test]
fn the_committed_windows_cover_every_input() {
    let fixture = window_fixture();
    let known: BTreeSet<&str> = fixture
        .transcripts
        .iter()
        .map(|t| t.id.as_str())
        .chain(fixture.genomic.iter().map(|w| w.contig.as_str()))
        .collect();
    let cases = cases();
    let missing: Vec<&str> = cases
        .cases
        .iter()
        .map(|case| case.input.as_str())
        .filter(|input| {
            let accession = input.split(':').next().unwrap_or_default();
            !known.contains(accession)
        })
        .collect();
    assert!(
        missing.is_empty(),
        "{WINDOWS_PATH} carries no reference for {missing:?} — regenerate with `{REGENERATE}`"
    );
}

/// Every captured transcript is served **whole**, through every entry point,
/// with bases identical to the stored record.
///
/// This encodes a measured defect, not a tidiness preference. The recorder files
/// the normalizer's `get_genomic_sequence` reads under the transcript's own
/// accession, which registers it in the fixture's `genomic` map, which makes
/// `WindowProvider::is_known_contig` true for it — and *that* routes
/// `get_sequence` for the transcript through the genomic-window path instead of
/// the whole stored sequence. With only the recorder's padded windows, every
/// 200 bp block of every transcript failed under `WindowProvider` while
/// succeeding under `MultiFastaProvider`.
///
/// The consequence was not a loud error but a silent pass: a normalization that
/// reads past the captured window bails and returns its input unchanged, and an
/// unchanged lone `delins` is exactly what a conformant answer looks like.
/// `serve_transcripts_whole` in the generator removes the failure mode; this
/// test is what stops a future regeneration from quietly dropping it.
#[test]
fn the_committed_windows_serve_every_transcript_whole() {
    let fixture = window_fixture();
    let provider = fixture.to_provider();
    assert!(
        !fixture.transcripts.is_empty(),
        "{WINDOWS_PATH} carries no transcripts at all — regenerate with `{REGENERATE}`"
    );
    for transcript in &fixture.transcripts {
        let id = &transcript.id;
        let bases = transcript
            .sequence
            .as_ref()
            .unwrap_or_else(|| panic!("{WINDOWS_PATH}: transcript {id} carries no sequence"));
        let length = bases.len() as u64;
        assert_eq!(
            provider.get_sequence_length(id).ok(),
            Some(length),
            "{id}: get_sequence_length disagrees with the stored sequence"
        );
        assert_eq!(
            provider
                .get_transcript(id)
                .unwrap_or_else(|e| panic!("{id}: get_transcript failed: {e}"))
                .sequence
                .as_ref(),
            Some(bases),
            "{id}: get_transcript served different bases than the stored transcript"
        );
        for (entry_point, served) in [
            ("get_sequence", provider.get_sequence(id, 0, length)),
            (
                "get_genomic_sequence",
                provider.get_genomic_sequence(id, 0, length),
            ),
        ] {
            let served = served.unwrap_or_else(|e| {
                panic!(
                    "{id}: {entry_point}(0, {length}) failed against the committed windows: {e}\n\
                     The slice does not serve this accession whole, so a pass that reads further \
                     than the recorded one bails and returns its input unchanged — which for a \
                     lone-member row is byte-identical to a correct answer. Regenerate with \
                     `{REGENERATE}`."
                )
            });
            assert_eq!(
                &served, bases,
                "{id}: {entry_point} served different bases than the stored transcript"
            );
        }
    }
}

/// The committed slice records no machine-local path.
///
/// The generator stores only the manifest's *basename*; an absolute path would
/// leak the local layout into a public-org repo and churn `--check` per machine.
#[test]
fn the_committed_slice_records_no_machine_local_path() {
    let fixture = window_fixture();
    for forbidden in ["/Users/", "/Volumes/", "/home/", "C:\\"] {
        assert!(
            !fixture.captured_from.contains(forbidden),
            "{WINDOWS_PATH} records a machine-local path in `captured_from`: {:?}",
            fixture.captured_from
        );
    }
    assert!(
        !fixture.captured_from.is_empty(),
        "{WINDOWS_PATH} records no provenance at all"
    );
}

// ---------------------------------------------------------------------------
// Premises: the geometry the records rest on, read from the bases
// ---------------------------------------------------------------------------

/// `NM_004006.2:c.76_83` is an inversion, and its two spellings coincide at the
/// interior.
///
/// The premise behind the flagship confluence class. Read from the committed
/// bases rather than restated, so a future fix cannot be credited to a change of
/// fixture: if the slice ever stops carrying `AATGCACA` here, this fails and
/// names the reason instead of silently re-basing the whole class.
#[test]
fn the_dmd_block_is_its_own_inversion_and_coincides_at_the_interior() {
    let fixture = window_fixture();
    let block = coding_bases(&fixture, "NM_004006.2", 76, 83);
    assert_eq!(block, "AATGCACA", "NM_004006.2:c.76_83");
    let inverted = reverse_complement(&block);
    assert_eq!(
        inverted, "TGTGCATT",
        "the delins payload the corpus pins must be revcomp(c.76_83) — DNA/inversion.md:5"
    );
    // Columns 78-81 are indices 2..6 of the 8-mer.
    assert_eq!(
        &block[2..6],
        &inverted[2..6],
        "the four interior columns 78-81 must coincide, which is what lets the block also be \
         spelled as two delins and is why this row exists at all"
    );
    assert_ne!(
        &block[..2],
        &inverted[..2],
        "the flanks must differ, or there is no variant here"
    );
}

/// `NM_000500.9:c.710_713` is an inversion whose interior is self-complementary.
///
/// The premise behind #1541, and the reason the pair is a genuine one-variant
/// case rather than two different edits: `TCGT` inverts to `ACGA`, the interior
/// `CG` is unchanged, and the flanks are exactly the two substitutions the
/// sibling spelling writes.
#[test]
fn the_cyp21a2_block_inverts_to_the_two_substitutions_the_sibling_row_writes() {
    let fixture = window_fixture();
    let block = coding_bases(&fixture, "NM_000500.9", 710, 713);
    assert_eq!(block, "TCGT", "NM_000500.9:c.710_713");
    let inverted = reverse_complement(&block);
    assert_eq!(inverted, "ACGA");
    assert_eq!(
        &block[1..3],
        &inverted[1..3],
        "the interior CG is self-complementary, so the inversion changes only the flanks"
    );
    // The sibling spelling is `c.[710T>A;713T>A]`: the reference base at both
    // flanks is T, and the inversion writes A at both.
    assert_eq!((&block[..1], &block[3..]), ("T", "T"));
    assert_eq!((&inverted[..1], &inverted[3..]), ("A", "A"));
}

/// `NM_000143.3:c.1_3` is the initiation codon, so `general.md:35` reaches it.
#[test]
fn the_fh_start_codon_triplet_reads_atg() {
    let fixture = window_fixture();
    assert_eq!(
        coding_bases(&fixture, "NM_000143.3", 1, 3),
        "ATG",
        "the codon exception of general.md:35 rests on both changes affecting one amino acid, and \
         the widely-quoted spelling `c.[1G>T;3T>A]` names reference bases this sequence does not \
         have (#1543)"
    );
}

/// #1542's row denotes what its input denotes.
///
/// This is what makes the row a *representation* defect rather than a denotation
/// one: applying ferro's three members to `CTGA` yields `GTTGG`, which is
/// exactly the input's payload.
#[test]
fn the_dup_shaped_split_denotes_the_same_bases_as_its_input() {
    let provider = window_fixture().to_provider();
    let span = genomic_bases(&provider, "NC_000017.11", 80110044, 80110047);
    assert_eq!(span, "CTGA");
    // g.80110044C>G, g.80110045dup, g.80110047A>G applied to CTGA.
    let applied = format!(
        "{}{}{}{}{}",
        "G",         // 80110044 C>G
        &span[1..2], // 80110045 T
        &span[1..2], // 80110045 duplicated
        &span[2..3], // 80110046 G
        "G",         // 80110047 A>G
    );
    assert_eq!(
        applied, "GTTGG",
        "the members must denote the input's own payload, or #1542 would be a denotation defect \
         rather than the representation question it is"
    );
}

/// The prior audit's claim about `NC_000002.11:g.47639670_47639673delinsTT` was
/// wrong, and this is why.
///
/// That audit recorded the row as returning itself. It splits — and the split is
/// sequence-equivalent: `AGTG` with `g.47639670_47639671del` and
/// `g.47639673G>T` applied yields `TT`, the input's payload.
#[test]
fn the_grch37_split_denotes_the_same_bases_as_its_input() {
    let provider = window_fixture().to_provider();
    let span = genomic_bases(&provider, "NC_000002.11", 47639670, 47639673);
    assert_eq!(span, "AGTG");
    let applied = format!("{}{}", &span[2..3], "T"); // AG deleted; 47639673 G>T
    assert_eq!(applied, "TT");
}

/// #1541's sibling question and #1542's are both representation questions, and
/// this pins the one that a prior audit misdiagnosed as a denotation bug.
///
/// `NM_006420.3:c.242..249` reads `CAAGGGTA`. The input
/// `c.[242C>A;247_249delinsTT]` and ferro's `c.[242del;245G>A;249A>T]` both
/// yield `AAAGGTT`, so the two spellings denote one sequence. The audit that
/// suspected a denotation bug (citing `delins.md:19`) was mistaken; what remains
/// is an unadjudicated representation question, which is why the row pins
/// stability only.
#[test]
fn the_arfgef2_pair_denotes_one_sequence() {
    let fixture = window_fixture();
    let span = coding_bases(&fixture, "NM_006420.3", 242, 249);
    assert_eq!(span, "CAAGGGTA");

    // `c.[242C>A;247_249delinsTT]` over CAAGGGTA: 242 C>A, 243..246 untouched,
    // 247..249 (GTA) replaced by TT.
    let from_input = format!("{}{}{}", "A", &span[1..5], "TT");
    // `c.[242del;245G>A;249A>T]` over CAAGGGTA: 242 dropped, 243..244 untouched,
    // 245 G>A, 246..248 untouched, 249 A>T.
    let from_output = format!("{}{}{}{}", &span[1..3], "A", &span[4..7], "T");
    assert_eq!(from_input, "AAAGGTT");
    assert_eq!(
        from_input, from_output,
        "the two spellings must denote one sequence, or this row is a denotation defect and not \
         the representation question the record says it is"
    );
}

// ---------------------------------------------------------------------------
// RED — live defects. These must stay red until fixed; never weaken them.
// ---------------------------------------------------------------------------

/// The decided target for #1541's class.
///
/// Not a guess and not a preference. `general.md:56` — "the preferred
/// description is: (1) substitution, (2) deletion, (3) inversion, (4)
/// duplication, (5) insertion" — ranks substitution above inversion, and
/// `tests/it/issue_1517_inv_priority_over_delins.rs` pins both directions of
/// that rule in one table. `c.710_713`'s competing members are two
/// substitutions, so the split wins, exactly as it does for #1230.
const CYP21A2_TARGET: &str = "NM_000500.9:c.[710T>A;713T>A]";

/// #1541 — `NM_000500.9:c.[710T>A;713T>A]`, `c.710_713inv` and
/// `c.710_713delinsACGA` are one variant and must produce
/// [`CYP21A2_TARGET`].
///
/// This asserts the **exact** string, not merely that the three agree, because
/// the target is decided (see [`CYP21A2_TARGET`]) and agreement on the wrong
/// form would otherwise satisfy the guard. On `main` the class splits 2 against
/// 1: the `inv` and `delins` spellings both settle on `c.710_713inv`, while the
/// substitution spelling stays put.
///
/// **The blocker is the axis gate, not the typing rule.**
/// `Normalizer::is_splittable_single_member` (`src/normalize/mod.rs:2005`)
/// matches only `HV::Genome`/`HV::Mt`, so a lone `c.` member never reaches the
/// sequence-first pass; #1230's control passes only because it is written on the
/// `g.` axis. PR #1484 widens the gate to `c.`/`n.`/`r.`, under which all three
/// spellings converge on the target. **Un-ignore this test in the change that
/// merges #1484** — do not weaken it.
///
/// Two green companions bound it. The premise is
/// [`the_cyp21a2_block_inverts_to_the_two_substitutions_the_sibling_row_writes`];
/// the control for the *opposite* direction of `general.md:56` is
/// [`the_inversion_preference_control_still_holds`].
#[test]
#[ignore = "#1541: c.[710T>A;713T>A], c.710_713inv and c.710_713delinsACGA do not converge. The \
            target is decided (general.md:56 ranks substitution above inversion); the blocker is \
            the c.-axis gate in is_splittable_single_member. Un-ignore when #1484 merges."]
fn the_four_mer_inversion_pair_converges() {
    let cases = cases();
    let provider = audited_provider();
    let members = cases.members_of("cyp21a2-710-713-four-mer-inversion");
    assert_eq!(members.len(), 3, "#1541 is a three-spelling class");

    let outputs: ClassOutcomes<'_> = members
        .iter()
        .map(|case| (case.input.as_str(), run(&provider, &case.input)))
        .collect();
    assert_the_slice_answered_everything(&provider, members.len());

    let rendered = outputs
        .iter()
        .map(|(input, result)| format!("  {input}\n    -> {result:?}"))
        .collect::<Vec<_>>()
        .join("\n");
    let distinct: BTreeSet<&str> = outputs
        .iter()
        .filter_map(|(_, result)| result.as_deref().ok())
        .collect();
    assert_eq!(
        distinct.len(),
        1,
        "#1541: one variant, {} normal forms. All three spellings were verified \
         sequence-equivalent (c.710_713 = TCGT, revcomp ACGA, self-complementary interior), so \
         this is not three edits:\n{rendered}",
        distinct.len(),
    );
    assert_eq!(
        distinct.iter().copied().next(),
        Some(CYP21A2_TARGET),
        "#1541: the class converged, but on the wrong form. general.md:56 ranks substitution (1) \
         above inversion (3), and the competing members here are two substitutions, so the \
         preferred description is {CYP21A2_TARGET} — the same answer as #1230, pinned in \
         tests/it/issue_1517_inv_priority_over_delins.rs:\n{rendered}"
    );
}

/// The control for the other direction of `general.md:56`: a whole-block reverse
/// complement competing with `delins` members must prefer the `inv`.
///
/// Green, and deliberately not folded into
/// [`the_confluence_classes_converge`]. Without it the rule #1541 asks for could
/// be implemented as "always prefer the split", which would satisfy the red test
/// and silently break the decided `inversion-vs-two-delins-76-83` ruling that
/// #1535 shipped. The two tests are a pair: the ranking discriminates on *what
/// the block competes with*, and only holding both sides pins that.
#[test]
fn the_inversion_preference_control_still_holds() {
    let cases = cases();
    let provider = audited_provider();
    let members = cases.members_of("dmd-76-83-whole-block-inversion");
    assert_eq!(members.len(), 3, "the control is a three-spelling class");

    let outputs: ClassOutcomes<'_> = members
        .iter()
        .map(|case| (case.input.as_str(), run(&provider, &case.input)))
        .collect();
    assert_the_slice_answered_everything(&provider, members.len());

    for (input, result) in &outputs {
        assert_eq!(
            result.as_deref(),
            Ok("NM_004006.2:c.76_83inv"),
            "`{input}` must stay typed as an inversion. Its competing members are `delins`, which \
             general.md:56 does not rank at all, so DNA/inversion.md:5 governs — the decided \
             ruling `inversion-vs-two-delins-76-83`, shipped by #1535. If this moved while \
             fixing #1541, the fix was written as \"always prefer the split\" rather than as the \
             ranking."
        );
    }
}

/// #1542 — a `dup`-shaped split still emits members at separation 0.
///
/// `NC_000017.11:g.80110044_80110047delinsGTTGG` becomes
/// `g.[80110044C>G;80110045dup;80110047A>G]`; under
/// `sep = next.start - prev.end - 1` the first two members are at separation 0,
/// which `DNA/delins.md:16` forbids. #1537 fixed three sibling rows of this
/// family and did not reach this one because of the `dup` shape.
///
/// The separation reading is itself part of the open question: the adjacency
/// depends on treating a `dup` as occupying the position it duplicates, and
/// under a zero-width (interbase) reading the pair is not adjacent at all. So
/// the fix may well be a ruling on how a `dup`'s extent is read rather than a
/// change to the splitter — but either way the row cannot stay as it is while
/// [`no_emitted_output_puts_two_members_on_consecutive_nucleotides`] enforces
/// the rule on every other row.
///
/// The premise is [`the_dup_shaped_split_denotes_the_same_bases_as_its_input`],
/// which passes: this is a representation defect, not a denotation one.
#[test]
#[ignore = "#1542: a dup-shaped split still emits members at separation 0 (delins.md:16), \
            surviving #1537. Red until the dup's extent is ruled on."]
fn the_dup_shaped_split_does_not_touch() {
    let cases = cases();
    let provider = audited_provider();
    let case: &Case = cases
        .cases
        .iter()
        .find(|case| case.input == "NC_000017.11:g.80110044_80110047delinsGTTGG")
        .expect("#1542's row is in the corpus");

    let output = run(&provider, &case.input);
    assert_the_slice_answered_everything(&provider, 1);

    let output = output.expect("#1542's row normalizes");
    let touching = adjacent_members(&output);
    assert!(
        touching.is_empty(),
        "#1542: `{}` -> `{output}` puts members on consecutive nucleotides at {touching:?}. \
         delins.md:16: \"changes involving two or more consecutive nucleotides are described as \
         deletion/insertion (delins) variants\".",
        case.input
    );
}

#[cfg(test)]
mod member_bounds_tests {
    use super::{adjacent_members, member_bounds, reverse_complement};

    #[test]
    fn a_single_member_description_has_no_boundaries() {
        assert_eq!(member_bounds("NM_004006.2:c.76_83inv"), None);
    }

    #[test]
    fn a_genomic_allele_yields_its_boundaries() {
        assert_eq!(
            member_bounds("NC_000017.11:g.[80110044C>G;80110045dup;80110047A>G]"),
            Some(vec![(80110044, 80110045), (80110045, 80110047)])
        );
    }

    #[test]
    fn a_range_member_uses_its_far_endpoint() {
        assert_eq!(
            member_bounds("NM_000038.6:c.[5265_5266delinsTG;5268T>G]"),
            Some(vec![(5266, 5268)])
        );
    }

    #[test]
    fn a_non_plain_coordinate_declines_rather_than_mis_parsing() {
        // A UTR member would otherwise parse as coordinate 6 rather than -6.
        assert_eq!(member_bounds("NM_024312.4:c.[-6_-3G[6];12A>T]"), None);
    }

    #[test]
    fn touching_members_are_reported_and_separated_ones_are_not() {
        assert_eq!(
            adjacent_members("NC_000017.11:g.[80110044C>G;80110045dup;80110047A>G]"),
            vec![(80110044, 80110045)]
        );
        assert!(adjacent_members("NM_000038.6:c.[5265_5266delinsTG;5268T>G]").is_empty());
    }

    #[test]
    fn reverse_complement_is_an_involution_on_acgt() {
        assert_eq!(reverse_complement("AATGCACA"), "TGTGCATT");
        assert_eq!(
            reverse_complement(&reverse_complement("AATGCACA")),
            "AATGCACA"
        );
    }
}
