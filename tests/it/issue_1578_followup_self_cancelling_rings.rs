//! #1578 follow-up (item 3): `general.md:58`'s self-cancellation rule does
//! **not** reach a ring's `::`-joined segments. Adjudicated 2026-08-10; the
//! record is `rulings[self-cancelling-across-ring-junctions]` in
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`, whose
//! citations are quote-verified against the spec checkout at build time.
//!
//! The ruling in one line: `general.md:58` is a four-space-indented sub-bullet
//! of `general.md:56`, the **prioritisation** bullet — sibling to `:57`, which
//! is explicitly a prioritisation rule. So `:58` inherits `:56`'s antecedent,
//! "when a description is possible according to several types", and works by
//! *redirecting* to the preferred single-type description. `DNA/complex.md:5`
//! defines complex as a range of changes that **cannot** be described as a basic
//! type, so a ring has no competing single-type description to redirect to: a
//! `delins` over the union span would describe a *linear* product and discard
//! the join. And `DNA/complex.md:130` adjudicates the two operators directly,
//! on identical member content — `"::" is used to indicate the join, instead of
//! ";" to describe two not connected deletions` — so a `::` composite is not
//! reducible to its member set the way a `[a;b]` allele is.
//!
//! **This file exists to stop the escape being closed the cheap way.** Pointing
//! `detect_self_cancelling_pair` at ring segments looks free, because a
//! well-formed ring's segments are telomere-anchored and disjoint by
//! construction (`complex.md:127` and the `sup` form at `:161` are the spec's
//! only published ring shapes), so the check would never fire on a legitimate
//! ring. That is exactly why it needs a guard rather than a comment: the change
//! would pass every existing test while encoding the claim `complex.md:130`
//! forecloses, and E3006's repair hint ("rewrite as a single delins") is
//! actively wrong advice for a ring.
//!
//! **The gap this file used to pin is now closed.**
//! `no_ring_wellformedness_rule_yet_so_malformed_rings_are_still_accepted` lived
//! here and asserted that `g.100_200del::150_250dup`, `g.100_200del::300_400dup`
//! and `g.150_250dup::100_200del` were still accepted — the status quo the ruling
//! declined to change, recorded so it was visible rather than implied. All three
//! are now rejected by `validate_ring_segments_are_wellformed`, and
//! `ring_segment_wellformedness.rs` owns those rows. The test was deleted rather
//! than inverted, exactly as its own failure message instructed.
//!
//! The ruling is untouched: `general.md:58` still does not cross `::`, and the two
//! tests here pin that with overlapping `del::del`, which the new rules
//! deliberately leave alone. Both were `del::dup` originally; had they stayed so,
//! they would now pass by being rejected for the *wrong* reason.

use ferro_hgvs::parse_hgvs;

/// The spec's own two ring shapes must parse, and must not be mistaken for
/// self-cancelling alleles. These are the rows a `:58`-based check would have to
/// leave alone, and the reason it is tempting to believe extending it is free.
///
/// Both are quoted from `DNA/complex.md:127` and `:161` — the bare ring and the
/// supernumerary ring. Their segments are telomere-anchored (`pter_…del` and
/// `…_qterdel`) and therefore disjoint, so `detect_self_cancelling_pair` finds
/// nothing in them either way.
#[test]
fn the_specs_own_ring_shapes_are_accepted() {
    for input in [
        "NC_000022.11:g.pter_(12200001_14700000)del::(37600001_410000000)_qterdel",
        "NC_000022.11:g.[pter_(12200001_14700000)del::(37600001_410000000)_qterdel]sup",
    ] {
        let variant = parse_hgvs(input)
            .unwrap_or_else(|e| panic!("the spec publishes `{input}` as a ring: {e}"));
        assert_eq!(
            variant.to_string(),
            input,
            "`{input}` must round-trip byte-identically",
        );
    }
}

/// The adjudicated negative, and the load-bearing test in this file: a ring
/// whose segments overlap is **accepted**, because `general.md:58` does not
/// reach across `::`.
///
/// The same overlap spelled as a cis allele is rejected — that pairing is
/// asserted in `the_identical_overlap_is_still_rejected_as_a_cis_allele` below,
/// so the two halves of the ruling are pinned together rather than separately.
///
/// If someone extends `detect_self_cancelling_pair` to ring segments, this test
/// fails and points at the ruling record. That is its whole purpose; do not
/// "fix" it by relaxing the assertion.
///
/// **The inputs are all-`del` deliberately.** They were originally `del::dup`
/// pairs — the literal shape of `general.md:58`'s own example — but
/// `validate_ring_segments_are_wellformed` now rejects a non-deletion segment on
/// separate grounds (`DNA/complex.md:39`/`:60-64`), which would have made this
/// test pass for the wrong reason: rejected, but by the junction rule rather than
/// by `:58`. Overlapping `del::del` isolates the question this ruling answers,
/// since the junction rule has nothing to say about it.
#[test]
fn an_overlapping_ring_is_not_rejected_as_self_cancelling() {
    for input in [
        // overlap at 150–200, the spec example's geometry with both segments
        // spelled as deletions
        "NC_000022.11:g.100_200del::150_250del",
        // identical spans — the strongest "cancelling" shape available
        "NC_000022.11:g.100_200del::100_200del",
        // and under the `sup` marker, which is blind through the same arms
        "NC_000022.11:g.[100_200del::150_250del]sup",
    ] {
        let variant = parse_hgvs(input).unwrap_or_else(|e| {
            panic!(
                "`{input}` must NOT be rejected as self-cancelling: \
                 `general.md:58` is a sub-bullet of `general.md:56`'s \
                 prioritisation rule and has no single-type description to \
                 redirect a ring to (`DNA/complex.md:5`), and \
                 `DNA/complex.md:130` distinguishes `::` from `;` on exactly \
                 this member content. See \
                 rulings[self-cancelling-across-ring-junctions]. Got: {e}"
            )
        });
        assert_eq!(
            variant.to_string(),
            input,
            "`{input}` must round-trip byte-identically",
        );
    }
}

/// The other half of the ruling: the identical overlap **is** rejected inside a
/// `[a;b]` cis allele, where `:58`'s antecedent does hold — the members are
/// co-applied to one linear reference (`DNA/alleles.md:5`) and the union-span
/// `delins` prioritisation redirects to exists.
///
/// Pinned here rather than relying on the existing E3006 tests because the
/// ruling is precisely about the *asymmetry* between the two operators. A change
/// that accidentally disabled the allele check would otherwise make the test
/// above pass for the wrong reason.
#[test]
fn the_identical_overlap_is_still_rejected_as_a_cis_allele() {
    for input in [
        "NC_000022.11:g.[100_200del;150_250dup]",
        "NC_000022.11:g.[100_200del;100_200dup]",
    ] {
        let error = parse_hgvs(input).expect_err(&format!(
            "`{input}` violates `general.md:58` as a cis allele"
        ));
        let rendered = error.to_string();
        assert!(
            rendered.contains("Self-cancelling allele"),
            "`{input}` must be rejected by the self-cancelling check specifically, \
             not by some other validator — got: {rendered}"
        );
    }
}

/// Non-overlapping members are legal, and this is the row that shows `:58` is an
/// *overlap* predicate and therefore the wrong instrument for ring
/// well-formedness — it accepts a non-overlapping pair however it is wired.
///
/// The `::` case is spelled `del::del` for the reason given above: the
/// `del::dup` form this test originally used is now rejected by the junction rule
/// (`DNA/complex.md:39`/`:60-64`), which is a different question. The cis-allele
/// half keeps its `dup`, since `:58` is exactly what governs it there.
#[test]
fn a_non_overlapping_pair_is_legal_in_both_operators() {
    for input in [
        "NC_000022.11:g.[100_200del;300_400dup]",
        "NC_000022.11:g.100_200del::300_400del",
    ] {
        let variant = parse_hgvs(input)
            .unwrap_or_else(|e| panic!("`{input}` has no overlap and is legal: {e}"));
        assert_eq!(variant.to_string(), input, "`{input}` must round-trip");
    }
}
