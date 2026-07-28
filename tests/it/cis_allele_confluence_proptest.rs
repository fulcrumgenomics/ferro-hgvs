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
//! as runs grouped into delins members must all reach one string.
//!
//! Cases are **valid by construction** — the encodings are generated from a
//! chosen set of changed positions over a known core, so each one provably
//! describes the same resulting sequence. There is no skip path: an unparseable
//! render or a failed normalize is a bug, not a case to discard.
//!
//! Scope: the genomic axis, matching where the sequence-first canonicalization
//! is enabled. Extending this to `c.`/`n.`/`r.` should follow that work.

use crate::common::synthetic::{SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};
use proptest::prelude::*;
use proptest::test_runner::Config as ProptestConfig;

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
        format!("NC_TEST.1:g.[{}]", members.join(";"))
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
        format!("NC_TEST.1:g.[{}]", members.join(";"))
    }

    fn encodings(&self) -> Vec<String> {
        vec![
            self.as_substitutions(),
            self.as_spanning_delins(),
            self.as_grouped_runs(),
        ]
    }
}

fn normalize(provider: &MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider.clone());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"));
    format!("{normalized}")
}

/// Cores biased toward repeats and complementary pairs, so 3'-shifting and
/// inversion recognition are actually exercised rather than a random string
/// where nothing shifts and nothing is a palindrome.
fn core_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec(
        prop_oneof![
            Just("A".to_string()),
            Just("C".to_string()),
            Just("G".to_string()),
            Just("T".to_string()),
            Just("AA".to_string()),
            Just("AT".to_string()),
            Just("GC".to_string()),
            Just("CAG".to_string()),
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
fn differing_base(reference: u8, choice: usize) -> char {
    let alternatives: Vec<char> = "ACGT"
        .chars()
        .filter(|base| *base as u8 != reference)
        .collect();
    alternatives[choice % alternatives.len()]
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

proptest! {
    #![proptest_config(ProptestConfig { cases: 512, ..ProptestConfig::default() })]

    /// Every encoding of one haplotype normalizes to one string.
    #[test]
    fn all_encodings_of_one_haplotype_converge(haplotype in haplotype_strategy()) {
        let provider = SyntheticBuilder::genomic(&haplotype.core).build();
        let encodings = haplotype.encodings();
        let normalized: Vec<String> = encodings
            .iter()
            .map(|input| normalize(&provider, input))
            .collect();

        let first = &normalized[0];
        for (input, output) in encodings.iter().zip(normalized.iter()) {
            prop_assert_eq!(
                output,
                first,
                "encodings of one haplotype diverged:\n  {} -> {}\n  {} -> {}",
                encodings[0], first, input, output
            );
        }
    }

    /// The canonical form is itself a fixed point, so convergence is to a
    /// genuine canonical form rather than to a shared way-point.
    #[test]
    fn the_converged_form_is_a_fixed_point(haplotype in haplotype_strategy()) {
        let provider = SyntheticBuilder::genomic(&haplotype.core).build();
        let once = normalize(&provider, &haplotype.as_substitutions());
        let twice = normalize(&provider, &once);
        prop_assert_eq!(&twice, &once, "`{}` is not a fixed point", once);
    }

    /// No normalized cis allele contains overlapping or out-of-order members —
    /// the structural half of #1235's acceptance criteria.
    #[test]
    fn members_are_disjoint_and_ascending(haplotype in haplotype_strategy()) {
        use ferro_hgvs::hgvs::variant::HgvsVariant;
        let provider = SyntheticBuilder::genomic(&haplotype.core).build();
        let normalizer = Normalizer::new(provider);
        for input in haplotype.encodings() {
            let parsed = parse_hgvs(&input).expect("valid by construction");
            let normalized = normalizer.normalize(&parsed).expect("normalize");
            let HgvsVariant::Allele(allele) = &normalized else {
                continue;
            };
            let mut previous_end: Option<i64> = None;
            for member in &allele.variants {
                let HgvsVariant::Genome(genome) = member else {
                    continue;
                };
                let interval = &genome.loc_edit.location;
                let (Some(start), Some(end)) = (interval.start.inner(), interval.end.inner())
                else {
                    continue;
                };
                let (start, end) = (start.base as i64, end.base as i64);
                if let Some(previous) = previous_end {
                    prop_assert!(
                        start > previous,
                        "`{}` -> `{}`: member at {} overlaps the previous member ending at {}",
                        input, normalized, start, previous
                    );
                }
                previous_end = Some(end);
            }
        }
    }
}

/// Non-vacuity guard: the strategy must actually produce both adjacent and
/// separated changed positions, otherwise the properties above could pass while
/// exercising only one partitioning shape.
#[test]
fn the_strategy_generates_both_adjacent_and_separated_changes() {
    use proptest::strategy::ValueTree;
    use proptest::test_runner::TestRunner;

    let mut runner = TestRunner::deterministic();
    let (mut adjacent, mut separated) = (false, false);
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
        if adjacent && separated {
            return;
        }
    }
    panic!("strategy never produced both adjacent and separated changes (adjacent={adjacent}, separated={separated})");
}
