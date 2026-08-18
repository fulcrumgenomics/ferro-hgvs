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
//! Green throughout, and **no row is red or `#[ignore]`d**: 29 of 34 rows
//! produce their recorded answer, all three confluence classes converge, every
//! answer is a fixed point, and no emitted output puts two members on
//! consecutive nucleotides — with no row exempt from that last check.
//!
//! A red row may never be weakened to green. A green re-classification is what
//! [`the_corpus_census_is_unchanged`] exists to catch, and it now asserts the
//! red set is **empty**, so re-introducing one is a deliberate edit.
//!
//! ## #1542 was ruled on, and its guard was replaced rather than deleted
//!
//! `NC_000017.11:g.80110044_80110047delinsGTTGG` was the last red row, pinned at
//! `g.[80110044C>G;80110045dup;80110047A>G]` and the one row exempted from the
//! corpus-wide separation check. Both of those are gone, for two different
//! reasons that are worth keeping apart:
//!
//! - **the instrument was wrong.** [`member_bounds`] read a `dup` as occupying
//!   the position it duplicates, so the pair reported separation 0. A `dup`
//!   inserts and consumes no reference position: on SPDI write footprints the
//!   members were at separation `[1, 1]`, and the exemption was buying nothing;
//! - **the real defect was elsewhere.** The partition depended on the *shuffle
//!   direction* — `live`/3' split where `live`/5' and all three sequence-first
//!   arms merged. Each answer was a fixed point, so no oracle could see it.
//!
//! So `the_dup_shaped_split_does_not_touch` is replaced by
//! [`the_dup_shaped_split_is_the_partition_of_its_own_resulting_sequence`],
//! which asserts the property the operator ruled on (2026-08-13) rather than
//! the one that never held.
//!
//! ## #1541 stopped reproducing, and its guard now runs
//!
//! The `cyp21a2-710-713-four-mer-inversion` class was red here: its three
//! spellings did not converge, because `is_splittable_single_member` matched
//! only `HV::Genome`/`HV::Mt`, so a lone `c.` member never reached the
//! sequence-first pass. **PR #1484 widened that gate to `c.`/`n.`/`r.` and
//! merged on 2026-08-08**, and all three spellings now converge. The rows
//! therefore no longer carry a `defect`, and
//! [`the_four_mer_inversion_pair_converges`] is no longer `#[ignore]`d — it had
//! been passing while gating nothing, since no CI job passes `--run-ignored`.
//!
//! **#1541 itself stays open, and this is not the claim that it is fixed.** The
//! issue asked for a *ruling on which form governs* and then convergence on it;
//! convergence arrived as a side effect and no ruling was made. So the class's
//! rows keep `expected: null` — the corpus's category for an open governing
//! ruling — and the exact form is held only as a **stability tripwire**, on
//! `background/basics.md:38`. See [`CYP21A2_TARGET`] for why the `general.md:56`
//! grounding this file used to give for that form does not stand.
//!
//! # One measurement in the brief did not reproduce
//!
//! `NM_000500.9:c.710_713delinsACGA` was reported as a fixed point when this
//! corpus was assembled. It was not: it normalized to `NM_000500.9:c.710_713inv`,
//! so the class of the day split **2 against 1** rather than three mutually
//! disagreeing ways. Kept because the figure was reported wrong once; the row's
//! `note` records the correction. Both readings are historical now — since PR
//! #1484 merged, all three spellings converge.
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
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{
    from_sequences_detailed, parse_hgvs, FromSequencesOptions, NormalizeConfig, Normalizer,
    ReferenceProvider, ShuffleDirection,
};

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
    run_in(provider, input, ShuffleDirection::ThreePrime)
}

/// [`run`] in an explicit shuffle direction.
///
/// Only [`the_dup_shaped_split_is_the_partition_of_its_own_resulting_sequence`]
/// asks for the 5' one, and it asks because #1542's defect was visible *only*
/// as a disagreement between the two: each direction's answer is a fixed point,
/// so a single-direction run cannot see it.
fn run_in<P: ReferenceProvider + Clone>(
    provider: &P,
    input: &str,
    direction: ShuffleDirection,
) -> Result<String, String> {
    let config = ErrorConfig::lenient();
    let parsed = parse_hgvs_with_config(input, config.clone())
        .map_err(|e| format!("preprocess error: {e}"))?;
    let normalize_config = NormalizeConfig::for_entry_point(direction, config);
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
/// on the reference's own axis, read from each member's **SPDI write
/// footprint**. Empty for a single-member description, which has a geometry but
/// no boundaries; `None` when `description` does not parse or when a member does
/// not convert — which declines rather than guesses.
///
/// # Why not the string (#1542)
///
/// This used to be a textual leading-number parse, on the stated ground that
/// `delins.md:16` constrains the *description* and re-deriving spans from the
/// producing code would check that code against itself. The second half of that
/// still holds and is why the footprint comes from
/// [`hgvs_to_spdi`](ferro_hgvs::spdi::hgvs_to_spdi) — which applies the
/// description to the reference and never asks the normalizer what a member
/// covers — rather than from `src/normalize/`. The first half was wrong: the
/// **positions a member writes** are not the integers it happens to spell.
///
/// Two shapes it got wrong, both load-bearing for the rule it feeds:
///
/// - a `dup` was measured as occupying the position it duplicates, so
///   `80110045dup` read as `[80110045, 80110045]`. A duplication inserts; it
///   consumes no reference position, and its footprint is empty at the junction
///   3' of the copied base. That single reading is what made #1542 present as a
///   `delins.md:16` violation at separation 0 when the members are in fact at
///   separation **1**;
/// - `101_102insT` read as `[101, 102]`, a two-base span. `insertion.md:15`
///   calls those two the **flanking** nucleotides — the insertion sits between
///   them and consumes neither — so the footprint is empty there too.
///
/// The interval is closed and 1-based so it reads against the HGVS coordinates
/// beside it, with an insertion returned empty (`hi == lo - 1`). That is what
/// makes `separation = next.lo - prev.hi - 1` uniform across edit kinds; the
/// same geometry `examples/dump_confluence_divergences.rs` documents, where
/// treating an insertion's `A_B` anchor as a consumed two-base span had already
/// invalidated one published distribution.
///
/// Members are sorted by footprint, because a separation is a property of the
/// geometry and not of the order the members are rendered in.
fn member_bounds<P: ReferenceProvider + ?Sized>(
    description: &str,
    provider: &P,
) -> Option<Vec<(i64, i64)>> {
    let spans = member_footprints(description, provider)?;
    Some(spans.windows(2).map(|w| (w[0].1, w[1].0)).collect())
}

/// Each member's own closed, 1-based SPDI write footprint, ascending — the
/// spans [`member_bounds`] reads its boundaries off. See that function for the
/// geometry and for why the footprint rather than the spelled endpoints.
///
/// A **single-member** description yields one footprint rather than `None`: it
/// has a geometry even though it has no boundaries, and
/// [`the_dup_shaped_split_is_the_partition_of_its_own_resulting_sequence`]
/// compares geometries between two descriptions that need not have the same
/// member count.
fn member_footprints<P: ReferenceProvider + ?Sized>(
    description: &str,
    provider: &P,
) -> Option<Vec<(i64, i64)>> {
    let members = match parse_hgvs(description).ok()? {
        HgvsVariant::Allele(allele) => allele.variants,
        single => vec![single],
    };
    let mut spans = Vec::with_capacity(members.len());
    for member in &members {
        let triple = hgvs_to_spdi(member, provider).ok()?;
        // `position` is a 0-based interbase offset, so `+ 1` is the first base
        // written; a pure insertion has an empty `deletion` and lands at
        // `hi == lo - 1`.
        let lo = i64::try_from(triple.position).ok()? + 1;
        let hi = lo + i64::try_from(triple.deletion.len()).ok()? - 1;
        spans.push((lo, hi));
    }
    spans.sort_unstable();
    Some(spans)
}

/// The members of `description` that sit on consecutive nucleotides.
fn adjacent_members<P: ReferenceProvider + ?Sized>(
    description: &str,
    provider: &P,
) -> Vec<(i64, i64)> {
    member_bounds(description, provider)
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
/// One sibling does assert such a form, and does it *as a declared exception*
/// rather than in spite of this rule: [`the_four_mer_inversion_pair_converges`]
/// pins [`CYP21A2_TARGET`] for #1541's class as a stability tripwire on
/// `background/basics.md:38`, not as an adjudication. Read that constant's doc
/// before adding a second one — a form pinned without that disclosure is the
/// freezing this paragraph forbids.
///
/// A red row's answer is pinned too, at today's **wrong** output, so fixing the
/// defect turns this test red as well as turning its `#[ignore]`d guard green.
/// That is correct and it is not a formality: it is what forced #1542's new
/// answer to be re-blessed here, with the ruling behind it written into the
/// row's `note`, rather than passing unnoticed. There is no red row today.
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
/// All three classes run here: none is red any more. #1541's class also carries
/// the dedicated [`the_four_mer_inversion_pair_converges`], which adds the form
/// pin this test deliberately withholds — where the governing ruling is open,
/// agreement is the claim and the form is not.
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

/// How many corpus rows produce a multi-member output, and so are actually
/// examined by [`no_emitted_output_puts_two_members_on_consecutive_nucleotides`].
///
/// The denominator of that guard's zero. It is pinned rather than floored
/// because the #1542 re-instrument changed *which* rows are readable, and a
/// floor cannot tell a widening from a narrowing.
///
/// **13 as of #2155, 15 before it, and 9 when this constant was first written.**
/// The 9 -> 15 widening predates #2155: it was measured on a base that predated
/// #1835's flip of the default partition rule to `canonical-coalesced`, under
/// which six further corpus rows were re-partitioned from the resulting
/// sequence into a multi-member output.
///
/// **15 -> 13 is #2155's own move, and it is a narrowing rather than the
/// widening the paragraph above describes.** `delins-payload-coincidence-carve-out-is-coding-dna-scoped`
/// (decided, scoped to the coding DNA axis) is superseded by #2155 to cover
/// every DNA axis, so `NC_000002.12:g.166052896_166052909delinsTGTG` and
/// `NC_000002.11:g.47639670_47639673delinsTT` — both previously split into
/// two-or-three-member outputs by the coding-only carve-out — now re-coalesce
/// on their `g.` axis into a single spanning `delins` member apiece. A
/// one-member output has no adjacent pair to examine at all, so both rows leave
/// the denominator: 15 - 2 = 13. The guard's own claim, that no emitted output
/// puts two members on consecutive nucleotides, is still **zero** over the
/// smaller set — the two rows that left were never violations, they merely
/// stopped being multi-member.
const MULTI_MEMBER_OUTPUTS_EXAMINED: usize = 13;

/// No emitted output puts two members on consecutive nucleotides
/// (`DNA/delins.md:16`).
///
/// A structural property, checked over the *whole* corpus rather than on the
/// three rows that motivated it, because that is what makes it a guard rather
/// than three more pins.
///
/// **No row is exempt any more (#1542).** The single exemption was the #1542
/// reproducer, and it existed because [`member_bounds`] measured a `dup` as
/// occupying the position it duplicates — so the row read as a violation that
/// its SPDI write footprints say it never was. With the instrument corrected
/// and the row's own re-derivation fixed, the corpus-wide guard covers every
/// row and there is no longer a hole in it. `Case::allows_adjacent_members` —
/// the predicate that granted the exemption — is deleted;
/// `DefectKind::AdjacentMembers` is **kept**, because classifying which guard a
/// future red row is expected to fail is a different job from waving that guard
/// off, and no row carries it today.
#[test]
fn no_emitted_output_puts_two_members_on_consecutive_nucleotides() {
    let cases = cases();
    let provider = audited_provider();

    let mut observed = Vec::new();
    for case in &cases.cases {
        observed.push(run(&provider, &case.input));
    }

    let mut violations = Vec::new();
    let mut declined = Vec::new();
    let mut multi_member_outputs = 0usize;
    for (case, actual) in cases.cases.iter().zip(&observed) {
        let Ok(output) = actual else { continue };
        // A `None` is a *decline*, and a decline reached by a bare `continue` is
        // a row dropped out of a corpus-wide guard with nothing counting it —
        // the shape `CLAUDE.md`'s "a generator must account for what it dropped"
        // is about. It is collected and refused below rather than skipped.
        let Some(bounds) = member_bounds(output, &provider) else {
            declined.push(format!("  {} -> {output}", case.input));
            continue;
        };
        if bounds.is_empty() {
            continue;
        }
        multi_member_outputs += 1;
        // Read off `bounds` rather than through `adjacent_members`, which would
        // re-derive the same geometry and fetch every member's footprint twice.
        for &(previous_end, next_start) in &bounds {
            if next_start > previous_end + 1 {
                continue;
            }
            violations.push(format!(
                "  {} -> {output}\n    members touch at {previous_end}/{next_start} — {}",
                case.input, case.citation.clause
            ));
        }
    }

    // After the geometry pass, not before it: reading each member's SPDI
    // footprint fetches reference bases too, so a window that cannot serve one
    // would otherwise turn into a silent `None` and a skipped row.
    assert_the_slice_answered_everything(&provider, cases.cases.len());

    assert!(
        declined.is_empty(),
        "{} output(s) yielded no member geometry, so this guard did not see them. A decline is \
         not a pass: `member_bounds` returns `None` only when the output fails to parse or a \
         member fails to convert to SPDI, and either is a finding rather than a row to skip. It \
         is also what evidences the #1542 widening — the textual reader declined a UTR position \
         or an intronic offset outright, and the footprint reader resolves them:\n{}",
        declined.len(),
        declined.join("\n")
    );

    // A zero violation count is only meaningful if multi-member outputs were
    // actually examined, and a *floor* cannot say whether the #1542 re-instrument
    // widened the guard or narrowed it — both readings satisfy `>= 5`. Pinned
    // exactly so a row that stops being examined is a failure rather than a
    // quietly smaller denominator.
    assert_eq!(
        multi_member_outputs, MULTI_MEMBER_OUTPUTS_EXAMINED,
        "the number of multi-member outputs this guard examines moved. Going up is usually a new \
         corpus row and needs re-blessing here; going down means a row stopped producing a \
         multi-member output or stopped being readable, and the zero above is then a smaller \
         claim than it was"
    );
    assert!(
        violations.is_empty(),
        "{} output(s) violate delins.md:16 (\"changes involving two or more consecutive \
         nucleotides are described as deletion/insertion (delins) variants\"):\n{}",
        violations.len(),
        violations.join("\n")
    );
}

/// How many members `description` has, reading no reference.
///
/// The count is the quantity #1542's rule is stated over, and it is a property
/// of the parsed description alone — so this deliberately does **not** go
/// through [`member_footprints`], whose `None` would conflate "not an allele"
/// with "a member did not convert".
fn member_count(description: &str) -> Option<usize> {
    Some(match parse_hgvs(description).ok()? {
        HgvsVariant::Allele(allele) => allele.variants.len(),
        _ => 1,
    })
}

/// The shuffle direction may move a member's placement; it may not change how
/// many members there are (#1542).
///
/// [`the_dup_shaped_split_is_the_partition_of_its_own_resulting_sequence`] pins
/// that rule on the one row that violated it. This states it as the **property**
/// it is, over every row in the corpus — which is the difference between a fix
/// and a regression pin, because `shift_and_coalesce_direction_symmetrically`
/// re-shifts an adopted mirror in the caller's direction and then coalesces
/// again, so the two directions' counts agreeing is an outcome of that pass
/// rather than something its shape guarantees.
///
/// Only the **count** is compared. The placement legitimately differs — that is
/// what the direction is *for* — so asserting string equality would pin 3' onto
/// 5' and forbid shuffling altogether.
///
/// Every row must answer in both directions: a row that errors on one and not
/// the other is a disagreement of the same kind, so declines are collected and
/// refused rather than skipped.
#[test]
fn the_member_count_does_not_depend_on_the_shuffle_direction() {
    let cases = cases();
    let provider = audited_provider();

    let mut compared = 0usize;
    let mut declined = Vec::new();
    let mut disagreements = Vec::new();
    for case in &cases.cases {
        let three_prime = run_in(&provider, &case.input, ShuffleDirection::ThreePrime);
        let five_prime = run_in(&provider, &case.input, ShuffleDirection::FivePrime);
        let (Ok(three_prime), Ok(five_prime)) = (&three_prime, &five_prime) else {
            declined.push(format!(
                "  {}\n    3': {three_prime:?}\n    5': {five_prime:?}",
                case.input
            ));
            continue;
        };
        let (Some(three_count), Some(five_count)) =
            (member_count(three_prime), member_count(five_prime))
        else {
            declined.push(format!(
                "  {}\n    3': {three_prime}\n    5': {five_prime}\n    (an output does not \
                 re-parse, so its member count is unreadable)",
                case.input
            ));
            continue;
        };
        compared += 1;
        if three_count != five_count {
            disagreements.push(format!(
                "  {}\n    3': {three_prime} ({three_count} member(s))\n    5': {five_prime} \
                 ({five_count} member(s))",
                case.input
            ));
        }
    }

    // Two normalizations per row, so the floor is twice the corpus.
    assert_the_slice_answered_everything(&provider, cases.cases.len() * 2);

    assert!(
        declined.is_empty(),
        "{} row(s) did not answer in both directions, so they say nothing either way:\n{}",
        declined.len(),
        declined.join("\n")
    );
    assert_eq!(
        compared,
        cases.cases.len(),
        "a zero below is only a claim about the rows that were compared"
    );
    assert!(
        disagreements.is_empty(),
        "{} row(s) emit a different number of members depending only on the shuffle direction. \
         Direction moves a member's placement; it may not change how many members there are, or \
         the partition is being read off a placement rather than off the sequence — which the \
         decided record `separation-is-a-property-of-the-spelling-not-of-the-variant` rules \
         against:\n{}",
        disagreements.len(),
        disagreements.join("\n")
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
    assert!(
        red.is_empty(),
        "the red set changed. Adding a row to it converts a regression into a recorded defect, \
         which is exactly the laundering this assertion exists to block; removing one claims a \
         defect is fixed, which must come with the `#[ignore]` being lifted. Two rows have left \
         this set, both that way. The three #1541 rows went when PR #1484 merged and \
         `the_four_mer_inversion_pair_converges` was un-ignored (they are still unpinned below, \
         because #1541's governing ruling is still open). #1542's went on the operator's \
         2026-08-13 ruling, with `the_dup_shaped_split_does_not_touch` replaced by the \
         un-ignored `the_dup_shaped_split_is_the_partition_of_its_own_resulting_sequence` — and \
         note the guard it left behind is a *different* one, because separation was never the \
         property that row violated. Found: {red:?}"
    );

    let issues: BTreeSet<&str> = cases
        .cases
        .iter()
        .filter_map(|case| case.defect.as_ref())
        .map(|defect| defect.issue.as_str())
        .collect();
    assert!(
        issues.is_empty(),
        "every red row must name a filed issue, and none still reproduces: {issues:?}"
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
// The adjudication-sensitive rows: #1541's class, its control, and #1542.
//
// **No live defect is left, and no guard here is `#[ignore]`d.** #1542 was ruled
// on 2026-08-13 and its guard was *replaced* rather than deleted — separation was
// never the property that row violated, so
// `the_dup_shaped_split_is_the_partition_of_its_own_resulting_sequence` asserts
// the ruled property and runs. #1541's class stopped reproducing when PR #1484
// merged, so its guard runs too — but the form it pins is unruled, which is why
// that pin is documented as a stability tripwire rather than a target.
// ---------------------------------------------------------------------------

/// The form #1541's class converges on today. **A stability pin, not an
/// adjudication** — read the second paragraph before citing this constant.
///
/// This file used to call it "the decided target", grounded on `general.md:56`
/// ("the preferred description is: (1) substitution, (2) deletion, (3)
/// inversion, (4) duplication, (5) insertion") ranking the competing
/// substitutions above the inversion. **That grounding is withdrawn.** The
/// `decided` ledger record `adjudication-precedence-order`, entry E1, says of
/// that clause: it "ranks single-variant TYPE LABELS FOR ONE SPAN — it never
/// ranks a multi-member allele against a spanning description, and an earlier
/// attempt to use it that way was refuted on exactly that ground — so `:56`
/// cannot settle a merge-versus-split question at all". `c.[710T>A;713T>A]`
/// against `c.710_713inv` is precisely a multi-member allele against a spanning
/// description, so `:56` does not reach it. The sibling `decided` record
/// `inversion-vs-two-delins-76-83` once read `:56` the other way ("a genuine
/// two-substitution shape splits"), and that tension **has since been resolved
/// by the operator, against the split**: the `decided` record
/// `whole-span-reverse-complement-types-as-inv` (2026-08-13) types a whole-span
/// reverse complement `inv` uniformly and supersedes the competitor-type
/// reasoning in both.
///
/// So the assertion below is held on `background/basics.md:38` — the corpus's
/// stated authority for its stability pins, which are explicitly **not**
/// spec-correctness claims. It fires if the form moves, which is what it is for:
/// #1541's own 2026-08-12 status records that this locus is a live member of the
/// #1703 family and that "whichever way that is ruled, this row moves with it".
/// **It has now been ruled, and this constant is on the losing side** — the
/// decided target is `NM_000500.9:c.710_713inv`. The pin is deliberately left
/// at today's value here, because the ruling ships no code and this file must
/// stay green until #1541/#1703 implement it; re-bless it **in that change**,
/// citing the ruling, rather than deleting the pin.
const CYP21A2_TARGET: &str = "NM_000500.9:c.[710T>A;713T>A]";

/// #1541 — `NM_000500.9:c.[710T>A;713T>A]`, `c.710_713inv` and
/// `c.710_713delinsACGA` are one variant and converge on [`CYP21A2_TARGET`].
///
/// **Two assertions with different standing, and the difference is the point.**
/// That the three agree is a correctness claim the project holds outright: the
/// `decided` record `confluence-gate-is-apply-equality-on-every-determined-axis`
/// makes "normalize is constant on each equivalence class" the release gate, and
/// these three spellings were verified sequence-equivalent from raw FASTA (see
/// [`the_cyp21a2_block_inverts_to_the_two_substitutions_the_sibling_row_writes`],
/// which reads the bases rather than restating them). That they agree on
/// *[`CYP21A2_TARGET`] specifically* is a **stability pin** on an unruled form —
/// read that constant's doc before treating it as an adjudication.
///
/// **Why it was `#[ignore]`d, and why it no longer is.** The class did not
/// converge while `Normalizer::is_splittable_single_member` matched only
/// `HV::Genome`/`HV::Mt`, so a lone `c.` member never reached the sequence-first
/// pass; #1230's control passed only because it is written on the `g.` axis. PR
/// #1484 widened that gate to `c.`/`n.`/`r.` and **merged on 2026-08-08**, and
/// the ignore reason — "un-ignore when #1484 merges" — was left in place. A
/// passing `#[ignore]`d test gates nothing: it ran on no PR for the five days
/// between, because no job in `ci.yml` passes `--run-ignored`. The one place CI
/// selects ignored tests at all is `external-validation.yml`'s weekly API
/// comparison, and it does not reach this test twice over — it filters to
/// `api_comparison_tests`, and `--ignored` selects *only* ignored tests, so
/// lifting the attribute here takes nothing away from it.
///
/// **This is not a claim that #1541 is fixed.** It asked for a ruling on which
/// form governs and then convergence on it; convergence arrived as a side
/// effect, no ruling was made, and the issue stays open.
///
/// Two green companions bound it. The premise is
/// [`the_cyp21a2_block_inverts_to_the_two_substitutions_the_sibling_row_writes`];
/// the control for the *opposite* direction — a whole-block reverse complement
/// competing with `delins` members, which must stay typed as an `inv` under the
/// `decided` `inversion-vs-two-delins-76-83` — is
/// [`the_inversion_preference_control_still_holds`].
#[test]
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
    // An error is a distinct form, not an absent one — the same reading
    // [`the_confluence_classes_converge`] uses. Dropping the `Err`s would let
    // two spellings stop normalizing while the third still emits the target,
    // and the class would report itself converged on one form. That mattered
    // little while this test was `#[ignore]`d; it gates CI now.
    let distinct: BTreeSet<String> = outputs
        .iter()
        .map(|(_, result)| match result {
            Ok(output) => output.clone(),
            Err(e) => format!("<{e}>"),
        })
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
        distinct.iter().map(String::as_str).next(),
        Some(CYP21A2_TARGET),
        "#1541: the class converged, but not on the form ferro has been shipping. This half is a \
         STABILITY pin on background/basics.md:38, not an adjudication — the form is unruled, so \
         do not simply re-point it at whatever came out. If a ruling on the #1703 family moved it \
         (see CYP21A2_TARGET), re-bless this constant citing that ruling; if nothing was ruled, \
         this is an undeclared representation change. Expected {CYP21A2_TARGET}:\n{rendered}"
    );
}

/// The control for the other direction of `general.md:56`: a whole-block reverse
/// complement competing with `delins` members must prefer the `inv`.
///
/// Green, and deliberately not folded into
/// [`the_confluence_classes_converge`]. Without it the rule #1541 asks for could
/// be implemented as "always prefer the split", which would silently break the
/// decided `inversion-vs-two-delins-76-83` ruling that #1535 shipped.
///
/// **What this control is now for has changed, and the change matters to
/// whoever implements #1541.** It used to pin the *discriminator* — that the
/// ranking turned on what the block competes with — so that only holding both
/// sides caught a fix written as "always prefer the split". The decided ruling
/// `whole-span-reverse-complement-types-as-inv` (2026-08-13) retires that
/// discriminator and types a whole-span reverse complement `inv` uniformly, so
/// the sibling row is now expected to move to the `inv` **and this row is
/// expected not to move at all**. The control therefore still earns its place,
/// against the opposite error: a fix written as "always prefer the split" is now
/// wrong in both directions, and a fix written as the uniform `inv` must leave
/// this row exactly where it is.
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
            "`{input}` must stay typed as an inversion. DNA/inversion.md:5 governs, because the \
             span is an exact reverse complement — the decided ruling \
             `whole-span-reverse-complement-types-as-inv`, which types this shape `inv` \
             uniformly and supersedes the competitor-type reasoning in \
             `inversion-vs-two-delins-76-83` (whose outcome for this row is unchanged). If this \
             moved while fixing #1541, the fix was written as \"always prefer the split\" — \
             which the ruling now rejects in both directions."
        );
    }
}

/// The window #1542's re-derivation is measured over.
///
/// Twenty bases of flank either side of the input's own
/// `g.80110044_80110047` span. The flank is not decoration: [`from_sequences_detailed`]
/// reads no reference, so a change resting on a window edge is placed at that
/// edge rather than where a wider read would put it — the condition its
/// `placement_bounded_by_window` flag reports. Twenty bases is far more than
/// this locus's ambiguous run, and the assertion below refuses a bounded
/// placement outright rather than trusting the margin.
const DUP_SHAPED_SPLIT_WINDOW: (u64, u64) = (80110024, 80110067);

/// #1542 — the emitted partition is the partition of its own resulting sequence.
///
/// # What this replaces, and why the old predicate was wrong
///
/// This test was `the_dup_shaped_split_does_not_touch`, `#[ignore]`d and red: it
/// asserted that `NC_000017.11:g.80110044_80110047delinsGTTGG` -> `g.[80110044C>G;
/// 80110045dup;80110047A>G]` puts no two members on consecutive nucleotides, on
/// `DNA/delins.md:16`.
///
/// **Separation was never the right predicate for this row.** The adjacency it
/// reported existed only in the instrument: [`member_bounds`] read a `dup` as
/// occupying the position it duplicates. A duplication inserts and consumes no
/// reference position, so the three members' SPDI write footprints are at
/// separation **[1, 1]** — the old output was not a `:16` violation, and the
/// corpus-wide guard's exemption for this row was buying nothing.
///
/// # What is asserted instead
///
/// The property the row does violate, and the one the operator ruled on
/// (2026-08-13): **the emitted partition must equal the partition re-derived
/// from the resulting sequence.** Both halves are checked, and they fail
/// independently:
///
/// 1. **Direction may not change the partition.** Seven of the eight
///    `FERRO_PARTITION` x direction configurations already emitted the merged
///    form; only `live`/3' produced the split, so one variant had two normal
///    forms differing by nothing but the shuffle direction. Each was a fixed
///    point, which is why no oracle saw it.
/// 2. **The member geometry is re-derivable.** [`from_sequences_detailed`] derives a
///    description from a `(reference, resulting)` pair alone — it holds no
///    provider, cuts with a different partitioner, and never sees this output —
///    so comparing member footprints against it is not the normalizer checking
///    itself. The pair handed to it is built from the *input's own stated
///    payload*, not from anything normalization produced.
///
/// The whole-string form is pinned separately, by
/// [`every_row_produces_its_recorded_answer`] against the fixture; what lives
/// here is the geometry, which is what the ruling is about.
///
/// The premise is [`the_dup_shaped_split_denotes_the_same_bases_as_its_input`],
/// which passes: this was a representation defect, not a denotation one.
#[test]
fn the_dup_shaped_split_is_the_partition_of_its_own_resulting_sequence() {
    let cases = cases();
    let provider = audited_provider();
    let case: &Case = cases
        .cases
        .iter()
        .find(|case| case.input == "NC_000017.11:g.80110044_80110047delinsGTTGG")
        .expect("#1542's row is in the corpus");

    let three_prime = run_in(&provider, &case.input, ShuffleDirection::ThreePrime);
    let five_prime = run_in(&provider, &case.input, ShuffleDirection::FivePrime);
    assert_the_slice_answered_everything(&provider, 2);

    assert_eq!(
        three_prime, five_prime,
        "#1542: `{}` normalizes to two different strings depending only on the shuffle \
         direction. This row is pinned on FULL STRING IDENTITY, which is stronger than the \
         corpus-wide rule — on this input the two directions must agree on placement as well \
         as on member count, because its merged form is placement-stable. The weaker \
         count-only property is \
         `the_member_count_does_not_depend_on_the_shuffle_direction`, and a failure there \
         means the partition is being read off a placement rather than off the sequence — \
         which the decided record \
         `separation-is-a-property-of-the-spelling-not-of-the-variant` rules against.",
        case.input
    );
    let output = three_prime.expect("#1542's row normalizes");

    // The `(reference, resulting)` pair, built without the normalizer: the
    // window's own bases, with the input's stated payload spliced over its
    // stated span. `the_dup_shaped_split_denotes_the_same_bases_as_its_input`
    // is the companion that pins those two strings against the fixture.
    //
    // Read through a plain provider rather than the audited one, because this
    // read cannot fail *quietly*: `genomic_bases` panics naming the span when
    // the committed slice cannot serve it, which is the fail-loud the audit
    // discipline exists to guarantee.
    let (window_start, window_end) = DUP_SHAPED_SPLIT_WINDOW;
    let (span_start, span_end) = (80110044u64, 80110047u64);
    let reference = genomic_bases(
        &window_fixture().to_provider(),
        "NC_000017.11",
        window_start,
        window_end,
    );
    // Derived from the two spans rather than written out, so widening the
    // window cannot silently splice over the wrong bases.
    let span_offset = (span_start - window_start) as usize;
    let span_len = (span_end - span_start + 1) as usize;
    assert_eq!(
        &reference[span_offset..span_offset + span_len],
        "CTGA",
        "the window is not aligned on the input's own span"
    );
    let resulting = format!(
        "{}GTTGG{}",
        &reference[..span_offset],
        &reference[span_offset + span_len..]
    );

    let derived = from_sequences_detailed(
        "NC_000017.11",
        window_start,
        &reference,
        &resulting,
        &FromSequencesOptions::default(),
    )
    .expect("the (reference, resulting) pair re-derives");
    assert!(
        !derived.placement_bounded_by_window(),
        "the re-derivation rested on a window edge, so it describes this window rather than \
         this locus — widen DUP_SHAPED_SPLIT_WINDOW rather than reading the comparison below"
    );
    let re_derived = derived.variant.to_string();

    // Unwrapped on each side before comparing, never compared as `Option`s.
    // `member_footprints` returns `None` when a member does not convert to SPDI,
    // and `assert_eq!(None, None)` is a **pass** — so an `Option` comparison
    // here would go green on the day both sides stopped converting, which is the
    // one outcome that must not read as agreement. Proven discriminating:
    // forcing either side to `None` fails this test rather than satisfying it.
    let output_footprints = member_footprints(&output, &provider).unwrap_or_else(|| {
        panic!(
            "#1542: the emitted `{output}` has no SPDI footprint, so this row's geometry cannot \
             be compared and the comparison below would pass vacuously"
        )
    });
    let re_derived_footprints = member_footprints(&re_derived, &provider).unwrap_or_else(|| {
        panic!(
            "#1542: the re-derived `{re_derived}` has no SPDI footprint, so this row's geometry \
             cannot be compared and the comparison below would pass vacuously"
        )
    });
    assert_eq!(
        output_footprints, re_derived_footprints,
        "#1542: `{}` -> `{output}`, whose members do not sit where a partition re-derived from \
         the resulting sequence puts them (`{re_derived}`). A member boundary interior to a run \
         of change is read off a placement, not off the sequence.",
        case.input
    );
}

#[cfg(test)]
mod member_bounds_tests {
    use super::{adjacent_members, member_bounds, reverse_complement, window_fixture};
    use ferro_hgvs::conformance::reference_window::WindowProvider;

    fn provider() -> WindowProvider {
        window_fixture().to_provider()
    }

    #[test]
    fn a_single_member_description_has_no_boundaries() {
        assert_eq!(
            member_bounds("NM_004006.2:c.76_83inv", &provider()),
            Some(Vec::new()),
            "one member has a footprint but no boundary between members"
        );
    }

    #[test]
    fn a_genomic_allele_yields_its_boundaries() {
        // **This pin was wrong until #1542**, and the correction is the whole
        // point of re-instrumenting `member_bounds`. It used to read
        // `[(80110044, 80110045), (80110045, 80110047)]`, from the leading
        // integer of each member's spelling — which puts `80110045dup` on the
        // base it duplicates and so reports the first pair at separation 0.
        //
        // A duplication *inserts*: it consumes no reference position, and its
        // SPDI write footprint is empty at the junction 3' of the copied base
        // (`lo = 80110046`, `hi = 80110045`). So the first boundary is
        // `(80110044, 80110046)` — separation 1 — and the second
        // `(80110045, 80110047)`, also 1. This output was never the
        // `delins.md:16` violation the old instrument reported.
        assert_eq!(
            member_bounds(
                "NC_000017.11:g.[80110044C>G;80110045dup;80110047A>G]",
                &provider()
            ),
            Some(vec![(80110044, 80110046), (80110045, 80110047)])
        );
    }

    #[test]
    fn an_insertion_occupies_no_position_between_its_flanking_nucleotides() {
        // `insertion.md:15` calls `a` and `b` in `a_b ins` the **flanking**
        // nucleotides — the insertion sits between them and consumes neither.
        // The old textual reader gave `101_102insT` the two-base span
        // `[101, 102]`, which double-counts and inverts the separation on both
        // sides. Sibling of the `dup` correction above, and the reason
        // `examples/dump_confluence_divergences.rs` documents this geometry:
        // reading the anchor as a consumed span had already invalidated one
        // published distribution.
        assert_eq!(
            member_bounds(
                "NC_000017.11:g.[80110040A>C;80110044_80110045insT]",
                &provider()
            ),
            Some(vec![(80110040, 80110045)]),
            "the insertion's footprint is empty at 80110044|80110045, so it starts at 80110045"
        );
    }

    #[test]
    fn a_range_member_uses_its_far_endpoint() {
        // The numbers are **transcript** offsets, not the `c.` numbers in the
        // string: an SPDI triple is stated on the reference sequence, so a `c.`
        // position resolves through `cds_start` (`NM_000038.6`'s is 60 in the
        // committed slice, so every member moves by the same +59). That is a frame
        // shift and not a geometry change — the separation this reads off is 1
        // either way — and it is uniform within one description, which is all a
        // boundary needs. The decided record
        // `c-and-n-positions-are-flat-transcript-offsets` is what makes that
        // resolution well defined.
        assert_eq!(
            member_bounds("NM_000038.6:c.[5265_5266delinsTG;5268T>G]", &provider()),
            Some(vec![(5325, 5327)])
        );
    }

    #[test]
    fn a_utr_member_is_resolved_rather_than_declined() {
        // The old textual reader declined here, because `-6` does not begin
        // with a digit and mis-parsing it as `6` would have been worse. Reading
        // the footprint resolves it properly against the transcript, so a whole
        // class of rows the corpus-wide guard used to skip is now covered.
        assert!(
            member_bounds("NM_024312.4:c.[-6_-3G[6];12A>T]", &provider()).is_some(),
            "a 5'UTR member resolves through the provider"
        );
    }

    #[test]
    fn touching_members_are_reported_and_separated_ones_are_not() {
        assert_eq!(
            adjacent_members("NC_000017.11:g.[80110044C>G;80110045G>T]", &provider()),
            vec![(80110044, 80110045)]
        );
        assert!(
            adjacent_members("NM_000038.6:c.[5265_5266delinsTG;5268T>G]", &provider()).is_empty()
        );
        assert!(
            adjacent_members(
                "NC_000017.11:g.[80110044C>G;80110045dup;80110047A>G]",
                &provider()
            )
            .is_empty(),
            "separation [1, 1] — see `a_genomic_allele_yields_its_boundaries`"
        );
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
