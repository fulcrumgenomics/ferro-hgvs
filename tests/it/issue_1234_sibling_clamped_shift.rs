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
//!
//! ## Where the fix lives
//!
//! These cases are fixed by the sweep-window clamp added for #1254
//! (`merge::clamp_sibling_crossing_shifts`), not by a clamp of their own. #1234
//! and #1254 are one defect seen from two sides: #1234 is the shape where the
//! shifted member lands *on top of* a sibling, so the damage is visible as an
//! overlap; #1254 is the shape where it clears the sibling entirely, so the
//! damage is silent. Preventing the crossing upstream resolves both.
//!
//! That was measured rather than assumed. Every test in this file passed
//! against #1254's fix with the overlap-based clamp originally proposed for
//! #1234 entirely absent, and across 143,360 two-member cis alleles the two
//! together left exactly the same 526 overlapping outputs as #1254 alone —
//! against 1,444 for the overlap clamp alone and 2,578 unfixed. So the coverage
//! is kept here and the second clamp was dropped.
//!
//! The 526 residual is a *different* defect: a deletion canonicalised into a
//! whole-tract repeat (`g.[1_2del;4del]` -> `g.[1_9T[7];9del]`) that overlaps
//! its sibling. `NaEdit::Repeat` claims no reference bases as far as the clamp
//! is concerned, so it is out of scope here and not regressed by it.

use ferro_hgvs::hgvs::variant::HgvsVariant;
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
/// reaches the same string — it did so before #1254's fix too, which is what
/// made the unclamped output detectably wrong rather than merely
/// non-canonical.
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
///
/// The sibling was written `60G>C` until the independent oracle below was
/// introduced: `SEQ[60]` is `T`, so that was a ref-mismatched input the lenient
/// default silently accepted, and the old `EquivalenceChecker` assertion could
/// not see it either. Corrected to `60T>C`.
#[test]
fn a_distant_sibling_does_not_clamp_the_shift() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.[4_6del;60T>C]"),
        "TEMPLATE:g.[6_8del;60T>C]"
    );
}

/// The invariant behind the fix: no normalized cis allele may contain
/// overlapping or out-of-order members.
#[test]
fn normalized_cis_members_never_overlap() {
    let inputs = [
        "TEMPLATE:g.[4_6del;8A>T]",
        "TEMPLATE:g.[5_7del;8A>T]",
        "TEMPLATE:g.[4_6del;60T>C]",
        "TEMPLATE:g.[4_6del;9G>T]",
    ];
    let normalizer = Normalizer::new(provider());
    for input in inputs {
        let parsed = parse_hgvs(input).expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        let HgvsVariant::Allele(allele) = &normalized else {
            // Collapsed to one member: trivially disjoint, but assert the shape
            // rather than passing over anything that is not an allele.
            assert!(
                matches!(normalized, HgvsVariant::Genome(_)),
                "`{input}` -> `{normalized}`: expected a genomic variant or an allele of them"
            );
            continue;
        };
        let mut previous_end: Option<i64> = None;
        for member in &allele.variants {
            // Every member of a normalized genomic allele must itself be
            // genomic with resolvable bounds. Skipping instead would let this
            // invariant pass vacuously on exactly the malformed shapes it
            // exists to rule out.
            let HgvsVariant::Genome(g) = member else {
                panic!("`{input}` -> `{normalized}`: member `{member}` is not genomic");
            };
            let interval = &g.loc_edit.location;
            let (Some(start), Some(end)) = (interval.start.inner(), interval.end.inner()) else {
                panic!("`{input}` -> `{normalized}`: member `{member}` has unresolvable bounds");
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

/// Apply `descriptor` to [`SEQ`] **independently of the normalizer**, by
/// converting each member to its SPDI triple and splicing the reference.
///
/// Returns `None` when a member cannot be converted, when a stated deletion
/// disagrees with the reference, or when two members claim the same base — an
/// overlapping description has no single well-defined resulting sequence, and
/// double-splicing one would invent a sequence neither side denotes.
fn apply(descriptor: &str) -> Option<String> {
    crate::common::cis_apply_oracle::apply_with(&provider(), SEQ, descriptor)
}

/// Normalization must never change the resulting sequence.
///
/// The oracle is deliberately **independent of the normalizer**: both the input
/// and its normalized form are applied to the reference from their SPDI triples
/// and the resulting sequences compared.
///
/// An earlier version of this test used `EquivalenceChecker::check`, which is
/// circular — it normalizes both operands, so it passes whenever normalization
/// corrupts both sides consistently, which is exactly the failure mode #1234
/// and #1254 are about. It also silently accepted the old overlapping output by
/// reading `g.[6_8del;8A>T]` as a plain `g.6_8del`. `apply` instead declines an
/// overlapping descriptor, so that shape now fails here rather than passing.
#[test]
fn normalization_preserves_the_resulting_sequence() {
    let inputs = [
        "TEMPLATE:g.[4_6del;8A>T]",
        "TEMPLATE:g.[5_7del;8A>T]",
        "TEMPLATE:g.[4_6del;60T>C]",
        "TEMPLATE:g.[4_6del;9G>T]",
    ];
    for input in inputs {
        let normalized = normalize_to_string(input);
        let from_input = apply(input).unwrap_or_else(|| panic!("`{input}` does not apply"));
        let from_output = apply(&normalized).unwrap_or_else(|| {
            panic!("`{input}` -> `{normalized}` has no well-defined resulting sequence")
        });
        assert_eq!(
            from_output, from_input,
            "`{input}` -> `{normalized}` changed the resulting sequence"
        );
    }
}

/// An uncertain allele — `g.[( … )]` — asserts only that its members are in
/// cis, not where they sit. Rewriting a member's position drops the predicted
/// marker and states something the input deliberately did not.
///
/// Every transform in `normalize_allele` gates on this flag, and the clamp
/// must too: an ungated clamp collapsed `g.[(4_6del;8A>T)]` to a bare
/// `g.5_8delinsT`, so the parentheses vanished and a prediction became an
/// assertion. The flag is always false for protein but genuinely reachable for
/// a DNA `g.[( … )]`, so it has to be checked. Skipping the clamp leaves the
/// members overlapping, which is the correct trade here — the input asserts
/// nothing about where they sit, and the overlap stays visible to a consumer.
#[test]
fn an_uncertain_allele_keeps_its_predicted_marker() {
    let input = "TEMPLATE:g.[(4_6del;8A>T)]";
    let normalized = normalize_to_string(input);
    assert!(
        normalized.contains('('),
        "the predicted marker must survive normalization, got `{normalized}`"
    );
    assert_eq!(normalized, "TEMPLATE:g.[(6_8del;8A>T)]");

    // The pinned output has two members claiming base 8 — the overlap this
    // trade accepts, per the module header. Assert it *is* still overlapping,
    // so the pin cannot outlive the trade: teaching the clamp to run under
    // uncertainty makes `apply` start returning a sequence here and fails this
    // line, rather than leaving a stale expectation passing forever.
    assert!(
        apply(&normalized).is_none(),
        "`{input}` -> `{normalized}` no longer overlaps — the uncertainty trade \
         has been fixed, so pin the clamped form instead of this one"
    );
}

/// The certain spelling of the same allele is unaffected — the uncertainty
/// gate must not disable the clamp generally.
#[test]
fn the_certain_spelling_of_that_allele_still_clamps() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.[4_6del;8A>T]"),
        "TEMPLATE:g.5_8delinsT"
    );
}
