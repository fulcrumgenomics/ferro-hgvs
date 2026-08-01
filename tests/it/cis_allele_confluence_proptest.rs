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
//! one haplotype spelled as separate substitutions, as one spanning delins, and
//! as runs grouped into delins members must all reach one form.
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
//! **Scope: substitution-only, and so length-preserving.** Every generated
//! haplotype replaces bases in place, which means no encoding and no normalized
//! output can contain a bare `del`, `ins` or `dup`. That bounds what this file
//! can reach: `shift_pieces_three_prime` only moves a piece for which
//! `is_pure_indel()` holds, so the sibling clamp of #1234 — where a member's
//! 3'-shift runs over a neighbour — is unreachable from here, and #1231/#1232/
//! #1233's indel-bearing shapes are not covered. Extending the model to
//! reference/alternate edit pairs reaches those; that work belongs with the
//! normalizer fixes it depends on and is tracked separately.
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
const CORE_PARTS: [&str; 8] = ["A", "C", "G", "T", "AA", "AT", "GC", "CAG"];

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
