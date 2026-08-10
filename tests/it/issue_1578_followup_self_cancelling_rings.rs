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
//! What the ruling does **not** settle is recorded in `no_ring_wellformedness_rule_yet`
//! below: three of the rings ferro accepts are malformed for reasons that have
//! nothing to do with self-cancellation.

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
#[test]
fn an_overlapping_ring_is_not_rejected_as_self_cancelling() {
    for input in [
        // overlap at 150–200, the shape of `general.md:58`'s own example
        "NC_000022.11:g.100_200del::150_250dup",
        // identical spans — the strongest "cancelling" shape available
        "NC_000022.11:g.100_200del::100_200dup",
        // and under the `sup` marker, which is blind through the same arms
        "NC_000022.11:g.[100_200del::150_250dup]sup",
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

/// Non-overlapping `del` + `dup` is legal both ways, and this is the row that
/// shows `:58` is the wrong instrument for the malformed rings noted below: a
/// `:58`-based check is an *overlap* predicate, so it would accept this one no
/// matter how it were wired.
#[test]
fn a_non_overlapping_del_dup_pair_is_legal_in_both_operators() {
    for input in [
        "NC_000022.11:g.[100_200del;300_400dup]",
        "NC_000022.11:g.100_200del::300_400dup",
    ] {
        let variant = parse_hgvs(input)
            .unwrap_or_else(|e| panic!("`{input}` has no overlap and is legal: {e}"));
        assert_eq!(variant.to_string(), input, "`{input}` must round-trip");
    }
}

/// **What the ruling does not settle, pinned as an open gap rather than as
/// correct behaviour.**
///
/// Each input below is accepted today and none is a well-formed ring, for
/// reasons independent of self-cancellation:
///
/// - a `dup` segment contributes no break junction (`DNA/complex.md:39`: "a
///   double colon is used to designate break point junctions creating a ring
///   chromosome"), and neither segment is telomere-anchored, so these describe
///   no ring geometry; and
/// - `g.150_250dup::100_200del` additionally lists its segments out of order,
///   which `DNA/complex.md:51` ("`pter` of the chromosome is to be listed
///   first"), `:53` and `:55` forbid.
///
/// This test asserts the **status quo**, deliberately, so the gap is visible and
/// counted rather than implied by the ruling's silence. When a ring
/// well-formedness rule lands, these inputs should start being rejected and this
/// test is expected to fail — that failure is the signal to delete it and move
/// the rows into the new rule's own tests, not to re-derive whether the ruling
/// above was wrong.
#[test]
fn no_ring_wellformedness_rule_yet_so_malformed_rings_are_still_accepted() {
    for input in [
        // no telomere anchor, and `dup` contributes no junction
        "NC_000022.11:g.100_200del::150_250dup",
        "NC_000022.11:g.100_200del::300_400dup",
        // segments out of pter->qter order
        "NC_000022.11:g.150_250dup::100_200del",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "`{input}` is accepted today. If this now fails, a ring \
             well-formedness rule has landed — move this row into that rule's \
             tests rather than reopening \
             rulings[self-cancelling-across-ring-junctions]",
        );
    }
}
