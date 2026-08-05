//! Property test for #1235's acceptance criterion: **confluence**.
//!
//! Idempotency is necessary but not sufficient. A normalizer can be perfectly
//! idempotent and still non-confluent — every spelling its own fixed point —
//! which is exactly the failure mode #1229-#1233 reported: several encodings of
//! one variant, each stable, none agreeing. A consumer keying on the normalized
//! string then silently over-splits one variant into several.
//!
//! `normalize_idempotency_proptest` already fuzzes confluence for *homopolymer
//! boundary* spellings of a single edit. This file covers the cis-allele case:
//! several encodings of one haplotype must all reach one form. Two models
//! generate those haplotypes — a substitution-only one and a length-changing
//! one — with different encoding sets and different reach; see **Two models**
//! below for what each can and cannot assert.
//!
//! Cases are **valid by construction** — the encodings are generated from a
//! chosen set of changed positions over a known core, so each one provably
//! describes the same resulting sequence. There is no skip path: an unparseable
//! render, a failed normalize, or an unexpected result shape is a bug, not a
//! case to discard.
//!
//! Convergence is checked against an **independent oracle**, not only against
//! the other encodings. Every encoding is normalized by the same implementation,
//! so agreement alone cannot distinguish "all correct" from "all wrong in the
//! same way"; each normalized result is therefore also re-applied to the
//! reference through `hgvs_to_spdi` — a different code path from the
//! sequence-first canonicalizer under test — and compared with the sequence the
//! haplotype was built to denote. Note `EquivalenceChecker` cannot serve as that
//! oracle: `check()` normalizes both sides, so comparing a normalized result
//! against anything returns `Identical` tautologically.
//!
//! **Two models.** The original `Haplotype` is substitution-only and so
//! length-preserving: no encoding and no normalized output can contain a bare
//! `del`, `ins` or `dup`. That bounded what this file could reach —
//! `shift_pieces_three_prime` only moves a piece for which `is_pure_indel()`
//! holds, so the sibling clamp of #1234 was unreachable from here, and
//! #1231/#1232/#1233's indel-bearing shapes were not covered.
//!
//! `IndelHaplotype` closes that gap: single-base deletions and one- or two-base
//! insertions, which do reach the shifting machinery. It is deliberately a
//! second model rather than a widening of the first — the substitution corpus
//! is tuned (see `core_strategy` and the non-vacuity guard) and runs at the
//! soak's 125k cases, and mixing the two would blur which shapes a failure came
//! from.
//!
//! What the indel model can assert today is narrower than what it should:
//!
//! - The **spanning-delins comparison is still held back in the indel model**,
//!   and the two issues that justified holding it back are now both fixed. A
//!   length-changing haplotype and its spanning `delins` reached two fixed points
//!   in general, not in two edge cases — the derivation reaches the *same* piece
//!   set from both spellings, and only the input-relative gates in
//!   `canonicalize_from_sequence` decided differently, because they were measured
//!   against how the input happened to be written. #1260 went with the two-gap
//!   alignment (which can express *insertion, retained base, insertion*, the
//!   separation threshold admitting the split across the retained base) and #1262
//!   with the removal of the input-separator veto, the last of those gates.
//!   Convergence is pinned by
//!   `adjacent_gap_insertions_and_the_delins_spelling_converge`.
//!
//!   So the restriction has outlived its reason. Restoring it is deliberately
//!   **not** done here: the test that holds it back,
//!   `an_indel_haplotype_normalizes_to_its_own_sequence`, is `#[ignore]`d because
//!   it already fails on #1394, so a comparison added to it could not be
//!   validated — a pass would be unobservable and a failure unattributable
//!   between #1394 and the new assertion. It waits on #1394, and is tracked
//!   rather than left to inattention. The substitution model already compares the
//!   delins encoding (`encodings()` includes `as_spanning_delins`), so the gap is
//!   specific to the indel model.
//!
//!   This is the **only** restriction those two issues justify. The generator
//!   used to carry a second one citing them — a filter excluding two insertions
//!   at adjacent gaps — which #1294 removed after re-measuring: neither #1260
//!   nor #1262 is about overlapping members, so the citation did not describe
//!   the shape it was excluding. What that shape hit was #1301, now fixed.
//! - `an_indel_haplotype_normalizes_to_its_own_sequence` is `#[ignore]`d
//!   because it finds #1394, a repeat whose growth exceeds its tract swallowing
//!   a sibling *deletion* that shifted into it. It found #1286, #1287, #1290,
//!   #1292, #1296, #1301, #1304, #1297, #1308, #1312, #1316, #1320, #1321,
//!   #1323 and #1325 first, all now fixed, each masked behind the one before
//!   it; see its doc comment for the history and the live reproduction.
//!
//! Neither restriction is silent: both are named, pinned, and fail loudly when
//! the underlying issue is fixed (#1268/#1283's complaint about coverage that
//! reads wider than it is).
//!
//! Axis scope: the genomic axis, matching where the sequence-first
//! canonicalization is enabled. Extending this to `c.`/`n.`/`r.` should follow
//! that work.

use crate::common::synthetic::{padded, SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};
use proptest::prelude::*;
use proptest::test_runner::Config as ProptestConfig;
use proptest::test_runner::TestCaseError;

/// HGVS position of 0-based index `i` into the core sequence.
fn hgvs_position(index: usize) -> u64 {
    PAD_OFFSET + 1 + index as u64
}

/// A haplotype: a core sequence plus the positions whose base changes, each
/// with its replacement. Substitution-only, so every encoding below is
/// length-preserving and provably describes the same resulting sequence.
#[derive(Debug, Clone)]
struct Haplotype {
    core: String,
    /// `(index into core, replacement base)`, strictly ascending by index.
    changes: Vec<(usize, char)>,
}

impl Haplotype {
    /// The core sequence with every change applied — what all encodings mean.
    ///
    /// This is the oracle's expected value. It is computed from the change list
    /// alone, without going near the normalizer, which is what lets it detect a
    /// canonicalization that silently altered the sequence.
    fn intended_sequence(&self) -> String {
        let mut bases: Vec<char> = self.core.chars().collect();
        for (index, alt) in &self.changes {
            bases[*index] = *alt;
        }
        bases.into_iter().collect()
    }

    /// One member per changed position: `g.[257A>G;259C>T]`.
    fn as_substitutions(&self) -> String {
        let members: Vec<String> = self
            .changes
            .iter()
            .map(|(index, alt)| {
                let reference = self.core.as_bytes()[*index] as char;
                format!("{}{reference}>{alt}", hgvs_position(*index))
            })
            .collect();
        render_members(&members)
    }

    /// One delins spanning the first to the last changed position, carrying the
    /// unchanged interior bases through unchanged.
    fn as_spanning_delins(&self) -> String {
        let first = self.changes.first().expect("non-empty").0;
        let last = self.changes.last().expect("non-empty").0;
        let mut replacement: Vec<char> = self.core[first..=last].chars().collect();
        for (index, alt) in &self.changes {
            replacement[index - first] = *alt;
        }
        let replacement: String = replacement.into_iter().collect();
        format!(
            "NC_TEST.1:g.{}_{}delins{replacement}",
            hgvs_position(first),
            hgvs_position(last)
        )
    }

    /// Maximal runs of consecutive changed positions, each rendered as one
    /// member: a lone position as a substitution, a run as a delins.
    fn as_grouped_runs(&self) -> String {
        let mut members: Vec<String> = Vec::new();
        let mut run: Vec<(usize, char)> = Vec::new();
        let flush = |run: &mut Vec<(usize, char)>, members: &mut Vec<String>, core: &str| {
            match run.len() {
                0 => {}
                1 => {
                    let (index, alt) = run[0];
                    let reference = core.as_bytes()[index] as char;
                    members.push(format!("{}{reference}>{alt}", hgvs_position(index)));
                }
                _ => {
                    let first = run[0].0;
                    let last = run[run.len() - 1].0;
                    let replacement: String = run.iter().map(|(_, alt)| *alt).collect();
                    members.push(format!(
                        "{}_{}delins{replacement}",
                        hgvs_position(first),
                        hgvs_position(last)
                    ));
                }
            }
            run.clear();
        };
        for change in &self.changes {
            if run
                .last()
                .is_some_and(|(previous, _)| change.0 != previous + 1)
            {
                flush(&mut run, &mut members, &self.core);
            }
            run.push(*change);
        }
        flush(&mut run, &mut members, &self.core);
        render_members(&members)
    }

    /// The distinct encodings of this haplotype, in a stable order.
    ///
    /// Deduplicated because two of the three renderings coincide for most
    /// haplotypes: `as_grouped_runs` is byte-identical to `as_substitutions`
    /// when no two changed positions touch (measured at 76% of generated cases),
    /// and to `as_spanning_delins` when they all do. Normalizing the same string
    /// twice and asserting it equals itself is a silent no-op that inflates the
    /// apparent coverage of the case budget; the non-vacuity guard below is what
    /// proves the three-distinct-encoding shape is actually reached.
    fn encodings(&self) -> Vec<String> {
        let mut out = Vec::with_capacity(3);
        for encoding in [
            self.as_substitutions(),
            self.as_spanning_delins(),
            self.as_grouped_runs(),
        ] {
            if !out.contains(&encoding) {
                out.push(encoding);
            }
        }
        out
    }
}

/// Render a cis-allele member list as one HGVS string.
///
/// A lone member is rendered **bare**, without the allele brackets. Bracketing
/// it would emit `g.[257_258delinsGT]`, which this crate's own W3026 rule
/// (`ErrorType::NonConformantBracketCardinality`) classifies as non-conformant
/// and `parse_hgvs_with_config` rejects under both the default and strict
/// configs — so a bracketed singleton would contradict this file's
/// valid-by-construction claim on the ~10% of cases whose changed positions form
/// a single run. It survives `parse_hgvs` only because the plain entry point
/// does not apply the cardinality rule.
fn render_members(members: &[String]) -> String {
    match members {
        [single] => format!("NC_TEST.1:g.{single}"),
        many => format!("NC_TEST.1:g.[{}]", many.join(";")),
    }
}

/// Parse and normalize, returning the variant rather than its rendering.
///
/// Deliberately not `common::synthetic::normalize_to_string`: these properties
/// assert on the normalized *variant* (its shape, and its SPDI projection), not
/// only on the string, and the shared helper discards it. It also takes the
/// `Normalizer` by reference so a case builds one and reuses it across
/// encodings — this file runs at the soak's 125k cases, where rebuilding the
/// provider and normalizer per encoding is the dominant cost.
fn normalize(normalizer: &Normalizer<MockProvider>, input: &str) -> HgvsVariant {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"))
}

/// The members of a normalized genomic result, as a flat list.
///
/// Fails rather than skips on any other shape: a genomic cis allele must
/// normalize either to a single genomic variant (every member coalesced) or to
/// an allele of genomic members. Anything else — a protein result, a trans
/// allele, a non-genomic member — is a defect in the code under test, and
/// `continue`-ing past it is what would let this property pass vacuously.
fn genomic_members(normalized: &HgvsVariant) -> Result<Vec<&HgvsVariant>, TestCaseError> {
    let members: Vec<&HgvsVariant> = match normalized {
        HgvsVariant::Allele(allele) => allele.variants.iter().collect(),
        single @ HgvsVariant::Genome(_) => vec![single],
        other => {
            return Err(TestCaseError::fail(format!(
                "normalizing a genomic cis allele produced `{other}`, which is neither a \
                 genomic variant nor an allele"
            )))
        }
    };
    for member in &members {
        if !matches!(member, HgvsVariant::Genome(_)) {
            return Err(TestCaseError::fail(format!(
                "`{normalized}`: member `{member}` is not a genomic variant"
            )));
        }
    }
    Ok(members)
}

/// The half-open interbase spans of a normalized result's members, ascending.
///
/// Spans come from `hgvs_to_spdi` rather than the printed HGVS endpoints
/// because the two disagree for insertions: `a_b ins` names two positions while
/// consuming neither, so two insertions at adjacent gaps share a printed
/// endpoint without overlapping. Interbase coordinates are the frame the
/// disjointness contract is actually stated in.
fn interbase_spans(
    normalized: &HgvsVariant,
    provider: &MockProvider,
) -> Result<Vec<(u64, u64)>, TestCaseError> {
    let mut spans = Vec::new();
    for member in genomic_members(normalized)? {
        let triple = hgvs_to_spdi(member, provider).map_err(|e| {
            TestCaseError::fail(format!(
                "`{normalized}`: member `{member}` has no SPDI: {e}"
            ))
        })?;
        let start = triple.position;
        spans.push((start, start + triple.deletion.len() as u64));
    }
    Ok(spans)
}

/// Re-apply a normalized variant to the reference through a code path that is
/// **not** the one under test, and return the sequence it denotes.
///
/// The sequence-first canonicalizer derives its answer with `normalize::merge`'s
/// own window applier, so comparing its output against itself proves nothing.
/// This goes through `spdi::convert::hgvs_to_spdi` instead — a separate
/// projection with a separate applier — so a canonicalizer that silently
/// changed the sequence shows up as a mismatch here even when every encoding
/// agrees with every other.
fn resulting_sequence(
    normalized: &HgvsVariant,
    provider: &MockProvider,
    reference: &str,
) -> Result<String, TestCaseError> {
    let mut triples = Vec::new();
    for member in genomic_members(normalized)? {
        let triple = hgvs_to_spdi(member, provider).map_err(|e| {
            TestCaseError::fail(format!(
                "`{normalized}`: member `{member}` has no SPDI: {e}"
            ))
        })?;
        triples.push(triple);
    }
    triples.sort_by_key(|triple| triple.position);

    let mut out = String::with_capacity(reference.len());
    let mut cursor = 0usize;
    for triple in &triples {
        let start = triple.position as usize;
        if start < cursor {
            return Err(TestCaseError::fail(format!(
                "`{normalized}`: members overlap — SPDI at {start} starts before {cursor}"
            )));
        }
        out.push_str(&reference[cursor..start]);
        out.push_str(&triple.insertion);
        cursor = start + triple.deletion.len();
    }
    out.push_str(&reference[cursor..]);
    Ok(out)
}

/// The building blocks a core is assembled from, biased toward repeats and
/// complementary pairs (see [`core_strategy`]).
///
/// [`UNION_CORE_PARTS`] holds the same parts for the indel model;
/// `the_two_core_part_lists_agree` fails if the two ever drift apart.
const CORE_PARTS: [&str; 8] = ["A", "C", "G", "T", "AA", "AT", "GC", "CAG"];

/// The parts [`union_core_strategy`] draws from, **frozen** at the values in
/// place when the seeds in
/// `tests/proptest-regressions/cis_allele_confluence_proptest.txt` were
/// recorded.
///
/// A separate constant rather than a reference to [`CORE_PARTS`], deliberately.
/// The union strategy has to keep drawing exactly what it drew when those seeds
/// were written (see its doc for why), so it must not follow a constant the
/// substitution model is free to edit: reordering [`CORE_PARTS`], or retitling
/// one entry, would otherwise silently retarget eight pinned regressions while
/// every test stayed green. `the_two_core_part_lists_agree` is what keeps the
/// duplication honest — it turns a silent divergence into a failing test.
const UNION_CORE_PARTS: [&str; 8] = ["A", "C", "G", "T", "AA", "AT", "GC", "CAG"];

/// The tripwire for the two constants above.
///
/// Both strategies' docs claim an identical distribution, and this is what
/// makes that claim checkable. If you changed [`CORE_PARTS`] on purpose you are
/// changing what the substitution model explores; decide separately whether the
/// indel model should follow, and if it should, re-record the pinned seeds in
/// the same commit rather than editing [`UNION_CORE_PARTS`] alone.
#[test]
fn the_two_core_part_lists_agree() {
    assert_eq!(
        CORE_PARTS, UNION_CORE_PARTS,
        "the core part lists diverged: `core_strategy` and `union_core_strategy` no longer draw \
         from the same set, so their docs' claim of an identical distribution is false. See the \
         note on `UNION_CORE_PARTS` before reconciling them."
    );
}

/// Cores biased toward repeats and complementary pairs, so inversion
/// recognition and repeat-aware partitioning are exercised rather than a random
/// string where nothing is a palindrome and no tract repeats.
///
/// This bias does **not** buy 3'-shifting: the corpus is substitution-only and
/// so length-preserving, and only a pure indel piece shifts (see the module
/// doc's scope note).
///
/// Drawn as **indices** into [`CORE_PARTS`] rather than as
/// `prop_oneof![Just("A".to_string()), ...]`. The two have the same
/// distribution — `prop_oneof!` weights its arms equally, so both are uniform
/// over the eight parts — but the union form dominated the soak's runtime.
/// Every case draws 6-13 parts, and each part cost a `TupleUnion` value tree
/// over eight `Just<String>` arms plus a `String` clone, none of it inlined in
/// the unoptimized test profile CI runs. Measured at 5,000 cases in that
/// profile: 2.52s of user CPU for this strategy alone against 0.14s for the
/// index form, an 18x difference that carried straight through to the
/// properties below (4.25s -> 1.9s each). At the soak's 125,000 cases it was
/// about 60% of the job's wall time.
fn core_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec(0..CORE_PARTS.len(), 6..14)
        .prop_map(|parts| parts.into_iter().map(|index| CORE_PARTS[index]).collect())
}

/// The original union-based draw, kept **only** for [`indel_haplotype_strategy`].
///
/// Identical in distribution to [`core_strategy`] — `prop_oneof!` weights its
/// arms equally, so both are uniform over the same eight parts, drawn here from
/// the frozen [`UNION_CORE_PARTS`] — and slower for the reason documented
/// there. It survives because
/// `tests/proptest-regressions/cis_allele_confluence_proptest.txt` pins eight
/// `IndelHaplotype` reproductions as RNG **seeds** — #1286, #1287, #1292, #1296,
/// #1297, #1308 and #1312, plus one unnumbered case the strategy's remaining
/// filter now rejects — and proptest replays a seed *through the strategy*:
/// change how the strategy consumes randomness and the same seed yields a
/// different value.
///
/// Verified rather than assumed, and the counterfactual is worse than a wrong
/// answer. With the draw preserved,
/// `an_indel_haplotype_normalizes_to_its_own_sequence --run-ignored all`
/// reproduces the recorded #1312 case (`core: "TAAAACCA"`) exactly and appends
/// nothing. Point the indel model at [`core_strategy`] instead and that same run
/// **passes**: #1312 is still open, but no seed reaches it any more, so the
/// property reads as fixed. Nothing goes red, because the property is
/// `#[ignore]`d and `the_ignored_indel_property_still_finds_its_defect` pins the
/// defect by hardcoded input rather than through the strategy.
///
/// The substitution model has no such pinned cases, so it is free to move.
/// Collapsing these two is safe once every pinned case is either fixed and
/// re-recorded or retired — all but #1312 are fixed and replay green — or once
/// the indel corpus runs at the soak's case count and needs the speed.
fn union_core_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec(
        prop_oneof![
            Just(UNION_CORE_PARTS[0].to_string()),
            Just(UNION_CORE_PARTS[1].to_string()),
            Just(UNION_CORE_PARTS[2].to_string()),
            Just(UNION_CORE_PARTS[3].to_string()),
            Just(UNION_CORE_PARTS[4].to_string()),
            Just(UNION_CORE_PARTS[5].to_string()),
            Just(UNION_CORE_PARTS[6].to_string()),
            Just(UNION_CORE_PARTS[7].to_string()),
        ],
        6..14,
    )
    .prop_map(|parts| parts.concat())
}

/// Replacement bases are chosen from the three that *differ* from the
/// reference, so every generated change is real and the strategy never has to
/// reject a case. Filtering after the fact instead (drop a change whose
/// replacement happened to equal the reference, then `prop_assume!` at least two
/// survive) rejects often enough to trip proptest's global-reject cap once the
/// case count is raised to the soak's 125k.
///
/// `choice` indexes the alternatives directly rather than modulo their count:
/// with an `ACGT` core there are always exactly three, and a silent `%` would
/// hide it if a core strategy ever introduced a base outside that alphabet.
fn differing_base(reference: u8, choice: usize) -> char {
    let alternatives: Vec<char> = "ACGT"
        .chars()
        .filter(|base| *base as u8 != reference)
        .collect();
    assert_eq!(
        alternatives.len(),
        3,
        "core base `{}` is not one of ACGT, so the replacement alphabet is not the assumed three",
        reference as char
    );
    alternatives[choice]
}

fn haplotype_strategy() -> impl Strategy<Value = Haplotype> {
    core_strategy().prop_flat_map(|core| {
        let length = core.len();
        // Two or three changed positions, which is where partitioning decisions
        // (adjacent vs separated) actually arise.
        prop::collection::btree_set(0..length, 2..=3).prop_flat_map(move |indices| {
            let core = core.clone();
            let indices: Vec<usize> = indices.into_iter().collect();
            let choices = prop::collection::vec(0..3usize, indices.len());
            choices.prop_map(move |choices| {
                let changes: Vec<(usize, char)> = indices
                    .iter()
                    .zip(choices.iter())
                    .map(|(index, choice)| {
                        (*index, differing_base(core.as_bytes()[*index], *choice))
                    })
                    .collect();
                Haplotype {
                    core: core.clone(),
                    changes,
                }
            })
        })
    })
}

/// A label for a normalized result's shape, for asserting that encodings agree
/// structurally and not merely in their rendering.
///
/// `AlleleVariant`'s `Display` short-circuits a single member to its bare
/// rendering, so a one-member `Allele` and the equivalent bare variant print
/// identically while remaining different types. A consumer that matches on
/// `HgvsVariant` — the Python bindings do — still sees the allele, which is the
/// over-splitting this file exists to rule out, one level up from the string.
fn variant_shape(variant: &HgvsVariant) -> String {
    match variant {
        HgvsVariant::Allele(allele) => format!("Allele({})", allele.variants.len()),
        HgvsVariant::Genome(_) => "Genome".to_string(),
        other => other.variant_type().to_string(),
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 512, ..ProptestConfig::default() })]

    /// Every encoding of one haplotype normalizes to one form — the same string
    /// *and* the same variant shape — and that form denotes the haplotype's
    /// intended sequence.
    ///
    /// The sequence check is what stops three encodings agreeing on a common
    /// wrong answer: all of them traverse the same implementation, so a defect
    /// that shifts the result uniformly leaves them agreeing and wrong.
    #[test]
    fn all_encodings_of_one_haplotype_converge(haplotype in haplotype_strategy()) {
        let provider = SyntheticBuilder::genomic(&haplotype.core).build();
        let normalizer = Normalizer::new(provider.clone());
        // The builder wraps the core in `PAD_OFFSET` bases on each side, so the
        // contig is `padded(core)` and SPDI's 0-based interbase positions index
        // straight into it.
        let reference = padded(&haplotype.core);
        let intended = padded(&haplotype.intended_sequence());

        let encodings = haplotype.encodings();
        let normalized: Vec<HgvsVariant> = encodings
            .iter()
            .map(|input| normalize(&normalizer, input))
            .collect();

        let first = &normalized[0];
        for (input, output) in encodings.iter().zip(normalized.iter()).skip(1) {
            prop_assert_eq!(
                output.to_string(),
                first.to_string(),
                "encodings of one haplotype diverged:\n  {} -> {}\n  {} -> {}",
                encodings[0], first, input, output
            );
            prop_assert_eq!(
                variant_shape(output),
                variant_shape(first),
                "encodings of one haplotype agree on the string `{}` but not on the \
                 variant shape:\n  {} -> {}\n  {} -> {}",
                first, encodings[0], variant_shape(first), input, variant_shape(output)
            );
        }

        for (input, output) in encodings.iter().zip(normalized.iter()) {
            let applied = resulting_sequence(output, &provider, &reference)?;
            prop_assert_eq!(
                &applied,
                &intended,
                "`{}` -> `{}` does not denote the haplotype's sequence",
                input, output
            );
        }
    }

    /// The canonical form is itself a fixed point, so convergence is to a
    /// genuine canonical form rather than to a shared way-point.
    ///
    /// Not redundant with `FERRO_ASSERT_IDEMPOTENT`: CI's oracle job runs
    /// `-E 'not test(proptest)'` (`ci.yml`), so the env-var oracle never sees
    /// this corpus. This property is its only fixed-point coverage.
    #[test]
    fn the_converged_form_is_a_fixed_point(haplotype in haplotype_strategy()) {
        let provider = SyntheticBuilder::genomic(&haplotype.core).build();
        let normalizer = Normalizer::new(provider);
        let once = normalize(&normalizer, &haplotype.as_substitutions()).to_string();
        let twice = normalize(&normalizer, &once).to_string();
        prop_assert_eq!(&twice, &once, "`{}` is not a fixed point", once);
    }

    /// No normalized cis allele contains overlapping or out-of-order members —
    /// the structural half of #1235's acceptance criteria.
    ///
    /// Ordering is compared on interbase SPDI spans, and every unexpected shape
    /// fails rather than being skipped, so the property cannot pass by asserting
    /// nothing. The sibling example test
    /// `issue_1235_cis_allele_confluence::normalized_cis_members_are_disjoint_and_ordered`
    /// pins the same contract on the specific #1229-#1234 reproducers; this is
    /// its fuzzed counterpart and must not be weaker.
    #[test]
    fn members_are_disjoint_and_ascending(haplotype in haplotype_strategy()) {
        let provider = SyntheticBuilder::genomic(&haplotype.core).build();
        let normalizer = Normalizer::new(provider.clone());
        for input in haplotype.encodings() {
            let normalized = normalize(&normalizer, &input);
            let spans = interbase_spans(&normalized, &provider)?;
            for pair in spans.windows(2) {
                let ((_, previous_end), (start, _)) = (pair[0], pair[1]);
                prop_assert!(
                    start >= previous_end,
                    "`{}` -> `{}`: member at interbase {} overlaps or does not follow the \
                     previous member ending at {}",
                    input, normalized, start, previous_end
                );
            }
        }
    }
}

/// Non-vacuity guard: the strategy must actually produce every shape the
/// properties above depend on — adjacent changed positions, separated ones, and
/// the mixed case that is the only one yielding three distinct encodings.
///
/// Without it the properties could pass while exercising a single partitioning
/// shape, and the deduplication in `encodings()` could silently reduce every
/// case to one encoding compared against itself.
#[test]
fn the_strategy_generates_every_shape_the_properties_depend_on() {
    use proptest::strategy::ValueTree;
    use proptest::test_runner::TestRunner;

    let mut runner = TestRunner::deterministic();
    let (mut adjacent, mut separated, mut three_encodings) = (false, false, false);
    for _ in 0..512 {
        let haplotype = haplotype_strategy()
            .new_tree(&mut runner)
            .expect("strategy produces a value")
            .current();
        for pair in haplotype.changes.windows(2) {
            if pair[1].0 == pair[0].0 + 1 {
                adjacent = true;
            } else {
                separated = true;
            }
        }
        if haplotype.encodings().len() == 3 {
            three_encodings = true;
        }
        if adjacent && separated && three_encodings {
            return;
        }
    }
    panic!(
        "strategy never produced every required shape \
         (adjacent={adjacent}, separated={separated}, three_encodings={three_encodings})"
    );
}

// ---------------------------------------------------------------------------
// Length-changing haplotypes (#1235 criterion 1, beyond substitutions)
// ---------------------------------------------------------------------------

/// One length-changing event on the core: a single-base deletion, or an
/// insertion of one or two bases at the gap *after* an index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Event {
    /// Delete `core[index]`.
    Delete { index: usize },
    /// Insert `bases` into the gap between `core[index]` and `core[index + 1]`.
    ///
    /// A `&'static str` payload rather than a `([char; 2], len)` pair: the
    /// latter made a `len` past the array's end constructible, and `inserted()`
    /// would then slice out of bounds.
    Insert { index: usize, bases: &'static str },
}

impl Event {
    fn index(&self) -> usize {
        match self {
            Event::Delete { index } | Event::Insert { index, .. } => *index,
        }
    }

    fn inserted(&self) -> &'static str {
        match self {
            Event::Delete { .. } => "",
            Event::Insert { bases, .. } => bases,
        }
    }
}

/// A haplotype built from length-changing events.
///
/// The substitution model above cannot reach a bare `del`, `ins` or `dup`, and
/// so cannot reach the shifting machinery at all: `shift_pieces_three_prime`
/// only moves a piece for which `is_pure_indel()` holds. This model reaches it.
#[derive(Debug, Clone)]
struct IndelHaplotype {
    core: String,
    /// Ascending by index, never two events on the same index.
    events: Vec<Event>,
}

impl IndelHaplotype {
    /// The core with every event applied — computed without the normalizer, so
    /// it can detect a canonicalization that silently altered the sequence.
    fn intended_sequence(&self) -> String {
        let mut out = String::new();
        for (index, base) in self.core.chars().enumerate() {
            let deleted = self
                .events
                .iter()
                .any(|e| matches!(e, Event::Delete { index: i } if *i == index));
            if !deleted {
                out.push(base);
            }
            for event in &self.events {
                if matches!(event, Event::Insert { index: i, .. } if *i == index) {
                    out.push_str(event.inserted());
                }
            }
        }
        out
    }

    /// The member string for one event.
    fn member(event: &Event) -> String {
        match event {
            Event::Delete { index } => format!("{}del", hgvs_position(*index)),
            Event::Insert { index, .. } => format!(
                "{}_{}ins{}",
                hgvs_position(*index),
                hgvs_position(*index) + 1,
                event.inserted()
            ),
        }
    }

    /// One member per event: `g.[261del;265_266insAC]`.
    fn as_separate_members(&self) -> String {
        let members: Vec<String> = self.events.iter().map(Self::member).collect();
        render_members(&members)
    }

    /// The same members authored in descending order — the input the
    /// order-independence property compares against. Shares `member` with
    /// `as_separate_members` so the two cannot drift.
    fn as_separate_members_reversed(&self) -> String {
        let members: Vec<String> = self.events.iter().rev().map(Self::member).collect();
        render_members(&members)
    }

    /// One delins spanning the first to the last affected base, carrying the
    /// unchanged interior through: the same haplotype, spelled as one member.
    fn as_spanning_delins(&self) -> String {
        let first = self.events.first().expect("non-empty").index();
        let last = self.events.last().expect("non-empty").index();
        let mut replacement = String::new();
        for (index, base) in self.core.chars().enumerate() {
            if index < first || index > last {
                continue;
            }
            let deleted = self
                .events
                .iter()
                .any(|e| matches!(e, Event::Delete { index: i } if *i == index));
            if !deleted {
                replacement.push(base);
            }
            for event in &self.events {
                if matches!(event, Event::Insert { index: i, .. } if *i == index) {
                    replacement.push_str(event.inserted());
                }
            }
        }
        if replacement.is_empty() {
            return format!(
                "NC_TEST.1:g.{}_{}del",
                hgvs_position(first),
                hgvs_position(last)
            );
        }
        format!(
            "NC_TEST.1:g.{}_{}delins{replacement}",
            hgvs_position(first),
            hgvs_position(last)
        )
    }

    fn encodings(&self) -> Vec<String> {
        let mut out = Vec::with_capacity(2);
        for encoding in [self.as_separate_members(), self.as_spanning_delins()] {
            if !out.contains(&encoding) {
                out.push(encoding);
            }
        }
        out
    }

    /// Whether two insertions sit at **adjacent gaps**.
    ///
    /// The strategy no longer filters this shape out (#1294). It is kept because
    /// `adjacent_gap_insertions_and_the_delins_spelling_converge` uses it to assert that the
    /// haplotype it pins really is one, so that test cannot drift onto a
    /// different shape than the one it claims to pin.
    fn has_adjacent_gap_insertions(&self) -> bool {
        self.events.windows(2).any(|pair| {
            matches!(pair[0], Event::Insert { .. })
                && matches!(pair[1], Event::Insert { .. })
                && pair[1].index() == pair[0].index() + 1
        })
    }
}

fn indel_haplotype_strategy() -> impl Strategy<Value = IndelHaplotype> {
    union_core_strategy().prop_flat_map(|core| {
        let length = core.len();
        // Leave the first and last base alone so an event never sits against
        // the core/padding seam, where a shift would leave the modelled region.
        prop::collection::btree_set(1..length.saturating_sub(1), 2..=3).prop_flat_map(
            move |indices| {
                let core = core.clone();
                let indices: Vec<usize> = indices.into_iter().collect();
                let kinds = prop::collection::vec(0..4usize, indices.len());
                let payloads = prop::collection::vec((0..4usize, 0..4usize), indices.len());
                (kinds, payloads)
                    .prop_map(move |(kinds, payloads)| {
                        // One- and two-base payloads over the full alphabet.
                        const ONE: [&str; 4] = ["A", "C", "G", "T"];
                        const TWO: [&str; 4] = ["AA", "AC", "GA", "TG"];
                        let events = indices
                            .iter()
                            .zip(kinds.iter())
                            .zip(payloads.iter())
                            .map(|((index, kind), (first, second))| match kind {
                                0 | 1 => Event::Delete { index: *index },
                                2 => Event::Insert {
                                    index: *index,
                                    bases: ONE[*first],
                                },
                                _ => Event::Insert {
                                    index: *index,
                                    bases: TWO[*second],
                                },
                            })
                            .collect();
                        IndelHaplotype {
                            core: core.clone(),
                            events,
                        }
                    })
                    // A haplotype whose events cancel denotes the reference itself.
                    // That is a real HGVS question (`=` and the self-cancelling
                    // merge, #1135) but not a *confluence* one: the spanning-delins
                    // encoding of a no-change region is degenerate, so the model
                    // would be comparing two spellings of "nothing happened".
                    .prop_filter("events cancel to no net change", |haplotype| {
                        haplotype.intended_sequence() != haplotype.core
                    })
            },
        )
    })
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 512, ..ProptestConfig::default() })]

    /// A length-changing haplotype normalizes to something that denotes it, is
    /// a fixed point, and renders disjoint members in ascending order.
    ///
    /// **`#[ignore]`d: this currently fails, and the failure is real.** Run it
    /// with `cargo test --features dev -- --ignored an_indel_haplotype`.
    ///
    /// It is committed rather than withheld because it keeps earning its place.
    /// Its first two runs found #1286 and #1287, both since fixed:
    ///
    /// ```text
    /// #1286  g.[258_259insA;259_260insA] -> g.[263dup;263dup]
    ///          two dup members at one junction              (fixed)
    /// #1287  g.[261_262insGA;263_264insAA] -> g.[262_263dup;263_266A[6]]
    ///          the dup's junction inside the repeat's span  (fixed)
    /// ```
    ///
    /// #1290 followed and is also fixed — a junction shifting past *another
    /// junction*, reordering two insertions. So are #1292, a duplication whose
    /// payload was rewritten when a clamp translated it; #1296, a repeat that
    /// did not report claiming the bases under its tract; #1301, two members
    /// sharing a span rendered out of order; #1304, a junction barrier reading a
    /// moved sibling's payload from the wrong snapshot; and #1297, a member that
    /// cancelled and was left as an identity member overlapping the repeat that
    /// absorbed it; #1308, a commuting payload sweeping past a sibling that
    /// stays put; and #1312, commuting tested against a payload the member no
    /// longer carries once it lands. The eleventh and last was **#1316**:
    ///
    /// ```text
    /// core "CAGCCAGTCAGCGCATCAG", insert "AA" after 4 and 5, delete index 12
    ///   g.[261_262insAA;262_263insAA;269del] -> g.[262A[3];262A[3];269del]
    ///     two identical repeats over one tract
    /// ```
    ///
    /// The seed's own shrink is smaller than that and carries no deletion at
    /// all — `core "ACAGCCAGTCAGCGCATCAG"` with `insAA` after 1 and 2 — so the
    /// deletion was incidental to the defect. The three-event form is kept here
    /// because it is the shape #1316 was filed from.
    ///
    /// The repeat-form analogue of #1286, which fixed the same shape spelled as
    /// duplications. (The output recorded when it was filed also lost the
    /// deletion and placed the repeats at 259; #1313 and #1317 moved both before
    /// #1316 itself was fixed.)
    ///
    /// The live failure is the seventeenth, **#1394**:
    ///
    /// ```text
    /// core "CAGGCAAACAGTGAAG", delete 5, insert "AA" after 6, "AA" after 7
    ///   g.[262del;263_264insAA;264_265insAA] -> g.[262_264A[7];264del]
    ///     the deletion claims 264, inside the tract 262-264
    /// ```
    ///
    /// #1325 was the sixteenth and is **fixed**, pinned by
    /// `issue_1325_repeat_growth_swallows_junction`: a repeat that grew its
    /// tract by more bases than the tract holds has no duplication form, so
    /// `demote_repeats_spanning_siblings` declined at its `start < a.start`
    /// guard and the tract kept a sibling's *junction*. It is now re-spelled as
    /// the insertion its growth stands for, which the clamp can bound.
    ///
    /// #1394 is the **base-claiming branch of that same guard**, and it was
    /// masked behind #1325: here the sibling is a deletion that 3'-shifted into
    /// the tract, so the decline stays. #1325's route is deliberately scoped
    /// away from it — an insertion blocks no sibling's shift, and widening the
    /// route costs `issue_1296_repeat_claims_its_bases`' two tests (measured).
    /// Verified pre-existing: `origin/main` at `dfea7d3` emits the identical
    /// output.
    ///
    /// Each of the seventeen was masked behind the one before it, which is the
    /// point of keeping this committed rather than deleting it until it passes.
    ///
    /// Enabling this test was #1292's acceptance criterion, then #1316's, then
    /// #1325's; all three are fixed and pinned, by
    /// `issue_1292_junction_payload_rotation`,
    /// `issue_1316_coincident_tract_repeats` and
    /// `issue_1325_repeat_growth_swallows_junction`, and by the `99d5d382`,
    /// `23c33174` and `4c6fab93` seeds below. None could carry the criterion,
    /// because the property's next failure was always behind it. Enabling now
    /// belongs to **#1394**, whose seed `6ec5686a` is committed below — so
    /// taking the ignore off replays it immediately and turns CI red, which is
    /// exactly what should gate taking it off.
    ///
    /// A caution about soak depth, and about reading a trend into it. With
    /// #1316 fixed this passed 200,000 cases and then failed a 1,000,000-case
    /// run at 27,633 successes (#1320); then at 114,694 (#1321), 143,034
    /// (#1323), 132,216 (#1325) and 231,801 (#1394). The early 4x jump did not
    /// continue — the middle four sit in one band — so depth has plateaued
    /// rather than trending toward an end, and it cannot be used to forecast
    /// one. A clean soak at one depth says nothing about the next, and a clean
    /// *seeded* soak says less than it looks: #1394 was found by an unseeded
    /// stream after all eight of `ci.yml`'s seeded shards passed clean at
    /// 1,000,000 cases with #1325 fixed.
    ///
    /// Do not weaken it to make it pass.
    #[test]
    #[ignore = "finds #1394, a real unfixed defect; see doc comment"]
    fn an_indel_haplotype_normalizes_to_its_own_sequence(
        haplotype in indel_haplotype_strategy()
    ) {
        let provider = SyntheticBuilder::genomic(&haplotype.core).build();
        let normalizer = Normalizer::new(provider.clone());
        let reference = padded(&haplotype.core);
        let intended = padded(&haplotype.intended_sequence());

        let input = haplotype.as_separate_members();
        let normalized = normalize(&normalizer, &input);

        let applied = resulting_sequence(&normalized, &provider, &reference)?;
        prop_assert_eq!(
            &applied,
            &intended,
            "`{}` -> `{}` does not denote the haplotype's sequence",
            input, normalized
        );

        let twice = normalize(&normalizer, &normalized.to_string()).to_string();
        prop_assert_eq!(
            &twice,
            &normalized.to_string(),
            "`{}` is not a fixed point",
            normalized
        );

        let spans = interbase_spans(&normalized, &provider)?;
        for pair in spans.windows(2) {
            let ((_, previous_end), (start, _)) = (pair[0], pair[1]);
            prop_assert!(
                start >= previous_end,
                "`{}` -> `{}`: member at interbase {} overlaps or does not follow the \
                 previous member ending at {}",
                input, normalized, start, previous_end
            );
        }
    }

    /// The rendered form does not depend on the order the members were authored
    /// in — #1098's contract, and the confluence this model *can* check without
    /// the delins spelling.
    #[test]
    fn an_indel_haplotype_is_authored_order_independent(
        haplotype in indel_haplotype_strategy()
    ) {
        let provider = SyntheticBuilder::genomic(&haplotype.core).build();
        let normalizer = Normalizer::new(provider);
        let forward = haplotype.as_separate_members();
        let reversed = haplotype.as_separate_members_reversed();
        prop_assume!(forward != reversed);
        prop_assert_eq!(
            normalize(&normalizer, &forward).to_string(),
            normalize(&normalizer, &reversed).to_string(),
            "authored order changed the result:\n  {}\n  {}",
            forward, reversed
        );
    }
}

/// Pins the shapes the indel model holds back, so the missing spanning-delins
/// comparison does not become permanent by inattention.
///
/// **#1260 and #1262 are both now fixed**, and this test pins the forms they
/// converged on. It used to assert that #1262 still diverged, with its failure
/// message naming the restoration of the delins encoding to
/// `an_indel_haplotype_normalizes_to_its_own_sequence`'s comparison set as the
/// follow-through. That restoration is still outstanding — see the module doc for
/// why it waits on #1394 rather than being done here.
///
/// The two rows are not two edge cases. They are the general behaviour of
/// "length-changing members versus the spanning delins": the derivation reaches
/// the *same* piece set from both spellings, and only the input-relative gates in
/// `canonicalize_from_sequence` — the `changed_columns_of_edits` bound and the
/// input-separator veto — decided differently, because both were measured against
/// how the input happened to be written. The veto is gone as of this change.
#[test]
fn adjacent_gap_insertions_and_the_delins_spelling_converge() {
    let provider = SyntheticBuilder::genomic("AAAAAA").build();
    let normalizer = Normalizer::new(provider);
    let normalize_str = |input: &str| normalize(&normalizer, input).to_string();

    // #1260, driven through the model itself so the held-back
    // `as_spanning_delins` encoding stays exercised rather than rotting: two
    // insertions of one base at adjacent gaps in a poly-A run.
    let haplotype = IndelHaplotype {
        core: "AAAAAA".to_string(),
        events: vec![
            Event::Insert {
                index: 1,
                bases: "C",
            },
            Event::Insert {
                index: 2,
                bases: "C",
            },
        ],
    };
    assert!(
        haplotype.has_adjacent_gap_insertions(),
        "the pinned haplotype must be an adjacent-gap-insertion one"
    );
    let encodings = haplotype.encodings();
    assert_eq!(
        encodings.len(),
        2,
        "the pinned haplotype must have two distinct encodings, got {encodings:?}"
    );
    let normalized: Vec<String> = encodings.iter().map(|e| normalize_str(e)).collect();
    // #1260 is FIXED. Both spellings converge, on the split form PR #1285 named:
    // the two-gap alignment `best_alignment` gained can express *insertion,
    // retained base, insertion*, and the separation threshold admits the split
    // across the one retained base.
    assert_eq!(
        normalized[0], normalized[1],
        "#1260's two spellings must converge"
    );
    assert_eq!(
        normalized[0], "NC_TEST.1:g.[258_259insC;259_260insC]",
        "and they converge on the split form, not the spanning delins"
    );

    // #1262: a substitution and a deletion against the spanning delins. Spelled
    // literally because the indel model has no substitution event.
    //
    // **Fixed** by removing the input-separator veto. The derivation always
    // reached one piece from both spellings; the veto refused the two-member one
    // because a derived piece covered the base that spelling had left between its
    // members, so the same variant was answered two ways depending on how it was
    // written.
    let split = normalize_str("NC_TEST.1:g.[258A>C;260del]");
    let spanning = normalize_str("NC_TEST.1:g.258_260delinsCA");
    assert_eq!(split, spanning, "#1262's two spellings must converge");
    assert_eq!(
        split, "NC_TEST.1:g.258_259delinsC",
        "and they converge on the form the block itself derives"
    );
}

/// The three clauses of `an_indel_haplotype_normalizes_to_its_own_sequence`,
/// evaluated on one hand-written case instead of over the strategy: the output
/// denotes the input's sequence, normalization is a fixed point, and the members
/// render disjoint and ascending.
///
/// `Ok(())` means the property holds for `input` against `core`; `Err` names the
/// clause that fails. [`the_ignored_indel_property_still_finds_its_defect`]
/// asserts `Err`, so fixing the pinned defect turns it red.
///
/// Stated over the property's own clauses rather than over a defect's *shape*
/// because the shapes differ: #1286 and #1287 rendered overlapping members,
/// while later defects change the sequence or break the fixed point instead. A
/// shape-specific assertion would have to be rewritten for each one, and — worse
/// — would report a defect as fixed while the property still failed on it some
/// other way.
fn indel_property_holds(core: &str, input: &str) -> Result<(), String> {
    let provider = SyntheticBuilder::genomic(core).build();
    let normalizer = Normalizer::new(provider.clone());
    let reference = padded(core);

    // The pinned inputs are well-formed by construction, so an input that does
    // not itself denote a sequence is a broken pin rather than a property
    // failure — panic rather than report it as the defect still standing.
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let intended = resulting_sequence(&parsed, &provider, &reference)
        .unwrap_or_else(|e| panic!("pinned input `{input}` has no resulting sequence: {e}"));

    let normalized = normalize(&normalizer, input);

    let applied = resulting_sequence(&normalized, &provider, &reference)
        .map_err(|e| format!("`{input}` -> `{normalized}` has no resulting sequence: {e}"))?;
    if applied != intended {
        return Err(format!(
            "`{input}` -> `{normalized}` does not denote the input's sequence"
        ));
    }

    let once = normalized.to_string();
    let twice = normalize(&normalizer, &once).to_string();
    if twice != once {
        return Err(format!(
            "`{once}` is not a fixed point: it normalizes to `{twice}`"
        ));
    }

    let spans = interbase_spans(&normalized, &provider)
        .map_err(|e| format!("`{input}` -> `{normalized}` has no interbase spans: {e}"))?;
    for pair in spans.windows(2) {
        let ((_, previous_end), (start, _)) = (pair[0], pair[1]);
        if start < previous_end {
            return Err(format!(
                "`{input}` -> `{normalized}`: member at interbase {start} overlaps or does not \
                 follow the previous member ending at {previous_end}"
            ));
        }
    }
    Ok(())
}

/// Pins the defect that currently keeps `an_indel_haplotype_normalizes_to_its_own_sequence`
/// `#[ignore]`d, so fixing it fails here and names the ignore to revisit.
///
/// The property's doc comment records the reproduction, but a doc comment does
/// not run: the live defect could be fixed and the property would simply stay
/// ignored, which is the coverage-that-reads-wider-than-it-is complaint in
/// #1268/#1283. This is the same guard
/// `adjacent_gap_insertions_and_the_delins_spelling_converge`
/// gives the two shapes the model holds back.
///
/// One row, not a list: the property stops at its first failure, so only the
/// *live* blocker is observable here. Each fix moves this pin to whatever the
/// property surfaces next, and the ignore comes off when there is nothing left
/// to pin.
#[test]
fn the_ignored_indel_property_still_finds_its_defect() {
    // #1394: the same `start < a.start` decline as #1325, on its base-claiming
    // branch — a sibling deletion 3'-shifts into the tract and the repeat is
    // left spanning it. Moved here from #1325, which this change fixes; #1325's
    // own junction case is pinned positively by
    // `issue_1325_repeat_growth_swallows_junction`.
    let (core, input, issue) = (
        "CAGGCAAACAGTGAAG",
        "NC_TEST.1:g.[262del;263_264insAA;264_265insAA]",
        "#1394",
    );
    assert!(
        indel_property_holds(core, input).is_err(),
        "{issue} appears fixed — the property now holds for `{input}`. Re-run \
         `an_indel_haplotype_normalizes_to_its_own_sequence` with `--ignored`, pin whatever \
         it finds next here, and drop its `#[ignore]` once it finds nothing."
    );
}

/// Non-vacuity guard for the indel model: it must actually produce deletions,
/// insertions, and haplotypes mixing the two, or the properties above are
/// asserting over a corpus narrower than they read as covering (#1268/#1283).
#[test]
fn the_indel_strategy_generates_every_shape_the_properties_depend_on() {
    use proptest::strategy::ValueTree;
    use proptest::test_runner::TestRunner;

    let mut runner = TestRunner::deterministic();
    let (mut deletion, mut insertion, mut mixed, mut multi_base) = (false, false, false, false);
    for _ in 0..512 {
        let haplotype = indel_haplotype_strategy()
            .new_tree(&mut runner)
            .expect("strategy produces a value")
            .current();
        let deletes = haplotype
            .events
            .iter()
            .filter(|e| matches!(e, Event::Delete { .. }))
            .count();
        let inserts = haplotype.events.len() - deletes;
        deletion |= deletes > 0;
        insertion |= inserts > 0;
        mixed |= deletes > 0 && inserts > 0;
        multi_base |= haplotype
            .events
            .iter()
            .any(|e| matches!(e, Event::Insert { bases, .. } if bases.len() == 2));
        if deletion && insertion && mixed && multi_base {
            return;
        }
    }
    panic!(
        "indel strategy never produced every required shape (deletion={deletion}, \
         insertion={insertion}, mixed={mixed}, multi_base={multi_base})"
    );
}
