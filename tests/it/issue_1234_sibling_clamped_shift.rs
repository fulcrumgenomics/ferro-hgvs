//! #1234 — a cis member's 3'-shift must stop at a sibling, not cross it.
//!
//! `Normalizer::normalize` 3'-shifted each cis-allele member to its
//! *standalone* most-3' position, ignoring its siblings. A deletion with a
//! downstream substitution therefore shifted **into** the substituted position,
//! emitting an allele whose two members overlap. That is malformed — a base
//! cannot be both deleted and substituted — and ferro itself read the result as
//! a different sequence, silently swallowing the substitution.
//!
//! The spec says nothing about how 3'-shifting interacts with siblings (the
//! only text that ever addressed it, SVD-WG010, was rejected), but a
//! description must not change the sequence it denotes, and `general.md:58`
//! forbids a description that removes part of the reference and replaces it
//! with part of the same sequence.

use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use std::io::Write;

const SEQ: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

fn provider() -> JsonProvider {
    let n = SEQ.len() as u64;
    let json = serde_json::json!({
        "version": "1.0",
        "transcripts": [{
            "id": "TEMPLATE-gene.1",
            "chromosome": "TEMPLATE",
            "strand": "+",
            "sequence": SEQ,
            "cds_start": 1,
            "cds_end": n - (n % 3),
            "genomic_start": 1,
            "genomic_end": n,
            "protein_id": "TEMPLATE-gene.1",
            "exons": [{
                "number": 1, "start": 1, "end": n,
                "genomic_start": 1, "genomic_end": n
            }]
        }],
        "genomic_sequences": { "TEMPLATE": SEQ }
    });
    let mut f = tempfile::NamedTempFile::new().expect("tempfile");
    write!(f, "{json}").expect("write json");
    JsonProvider::from_json(f.path()).expect("load reference")
}

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(provider());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"));
    format!("{normalized}")
}

/// The reported case. Reference `CACCA` at 4-8: deleting 4-6 leaves the `C` at
/// 7 and the substituted `A` at 8. The deletion's standalone 3'-most position
/// is `6_8`, which collides with the substitution; clamped at `5_7` it is
/// adjacent to it instead, and the two coalesce into one delins.
#[test]
fn deletion_shift_stops_before_a_sibling_substitution() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.[4_6del;8A>T]"),
        "TEMPLATE:g.5_8delinsT"
    );
}

/// The same variant spelled with the deletion already at the clamp position
/// reaches the same string — it did so before this fix, which is what made the
/// unclamped output detectably wrong rather than merely non-canonical.
#[test]
fn the_clamped_spelling_reaches_the_same_string() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.[5_7del;8A>T]"),
        "TEMPLATE:g.5_8delinsT"
    );
}

/// A lone deletion has no sibling to clamp against and must still reach its
/// standalone 3'-most position. Guards against "fixing" the clamp by disabling
/// the shift.
#[test]
fn a_lone_deletion_still_shifts_to_its_three_prime_most_position() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.4_6del"),
        "TEMPLATE:g.6_8del"
    );
}

/// Members far apart must be untouched by the clamp: the deletion still shifts
/// fully because its sibling is nowhere near.
#[test]
fn a_distant_sibling_does_not_clamp_the_shift() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.[4_6del;60G>C]"),
        "TEMPLATE:g.[6_8del;60G>C]"
    );
}

/// The invariant behind the fix: no normalized cis allele may contain
/// overlapping or out-of-order members.
#[test]
fn normalized_cis_members_never_overlap() {
    use ferro_hgvs::hgvs::variant::HgvsVariant;

    let inputs = [
        "TEMPLATE:g.[4_6del;8A>T]",
        "TEMPLATE:g.[5_7del;8A>T]",
        "TEMPLATE:g.[4_6del;60G>C]",
        "TEMPLATE:g.[4_6del;9G>T]",
    ];
    let normalizer = Normalizer::new(provider());
    for input in inputs {
        let parsed = parse_hgvs(input).expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        let HgvsVariant::Allele(allele) = &normalized else {
            continue; // collapsed to one member: trivially disjoint
        };
        let mut previous_end: Option<i64> = None;
        for member in &allele.variants {
            let HgvsVariant::Genome(g) = member else {
                continue;
            };
            let interval = &g.loc_edit.location;
            let (Some(start), Some(end)) = (interval.start.inner(), interval.end.inner()) else {
                continue;
            };
            let (start, end) = (start.base as i64, end.base as i64);
            if let Some(previous) = previous_end {
                assert!(
                    start > previous,
                    "`{input}` -> `{normalized}`: member at {start} overlaps the previous \
                     member ending at {previous}"
                );
            }
            previous_end = Some(end);
        }
    }
}

/// Normalization must never change the resulting sequence. `EquivalenceChecker`
/// read the old overlapping output as a plain `6_8del` — a *different* variant
/// — so this is the check that actually catches the corruption.
#[test]
fn normalization_preserves_the_resulting_sequence() {
    use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};

    let inputs = [
        "TEMPLATE:g.[4_6del;8A>T]",
        "TEMPLATE:g.[5_7del;8A>T]",
        "TEMPLATE:g.[4_6del;60G>C]",
        "TEMPLATE:g.[4_6del;9G>T]",
    ];
    let checker = EquivalenceChecker::new(provider());
    let normalizer = Normalizer::new(provider());
    for input in inputs {
        let parsed = parse_hgvs(input).expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        let level = checker.check(&parsed, &normalized).expect("check").level;
        assert!(
            matches!(
                level,
                EquivalenceLevel::Identical
                    | EquivalenceLevel::NormalizedMatch
                    | EquivalenceLevel::SequenceMatch
            ),
            "`{input}` -> `{normalized}` changed the resulting sequence ({level:?})"
        );
    }
}

/// An uncertain allele — `g.[( … )]` — asserts only that its members are in
/// cis, not where they sit. Rewriting a member's position drops the predicted
/// marker and states something the input deliberately did not.
///
/// Every neighbouring transform in `normalize_allele` already gates on this
/// flag. The clamp did not, so `g.[(4_6del;8A>T)]` collapsed to a bare
/// `g.5_8delinsT`: the parentheses vanished, turning a prediction into an
/// assertion. The PR's own comment claimed the flag is always false here, which
/// holds for protein but not for a DNA `g.[( … )]`.
#[test]
fn an_uncertain_allele_keeps_its_predicted_marker() {
    let normalized = normalize_to_string("TEMPLATE:g.[(4_6del;8A>T)]");
    assert!(
        normalized.contains('('),
        "the predicted marker must survive normalization, got `{normalized}`"
    );
    assert_eq!(normalized, "TEMPLATE:g.[(6_8del;8A>T)]");
}

/// The certain spelling of the same allele is unaffected — the gate must not
/// disable the clamp generally.
#[test]
fn the_certain_spelling_of_that_allele_still_clamps() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.[4_6del;8A>T]"),
        "TEMPLATE:g.5_8delinsT"
    );
}
