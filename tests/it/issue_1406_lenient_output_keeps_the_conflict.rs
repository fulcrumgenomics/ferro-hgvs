//! A conflict the input really had must survive into the output (#1406 row 3).
//!
//! W5002 says an overlap-conflicting allele has no canonical form, so ferro
//! preserves it verbatim and strict mode rejects it. But strict mode re-reads a
//! **description**, with no memory of what it was normalized from — so the
//! contract only holds if the emitted description still *looks* conflicting.
//!
//! It did not. The per-member pipeline repairs the members one at a time (an
//! inversion whose span is its own reverse complement cancels to an identity, an
//! interior insertion respells as a repeat) and the merge collapses what is
//! left, so lenient mode turned a description ferro refuses into one it admits:
//!
//! ```text
//! in      g.[11_12inv;11_12insAA]   strict: REJECT (W5002)
//! lenient g.10_11A[4]               strict: ACCEPT
//! ```
//!
//! ## Why preserving is the answer, and not repairing
//!
//! The independent applier declines an overlap-conflicting allele outright —
//! `ApplyFailure::Overlapping`, because "an overlapping description has no
//! single well-defined resulting sequence, and silently double-splicing one
//! would invent a sequence neither side denotes". So there are no "input's
//! bases" for a repair to preserve: any answer the repair produces comes from
//! *choosing* an apply order, which is precisely what W5002 refuses to do.
//!
//! The spec does not decide this either way — searched `alleles.md`,
//! `general.md` and every DNA/RNA recommendation, and no rule addresses two
//! members writing to the same reference bases. So the applier's refusal is the
//! strongest evidence available, and it points at preservation.
//!
//! ## The gate is on a PLURAL output, and that term is load-bearing
//!
//! Both detectors return nothing below two members, so for a single-member
//! output "the output does not conflict" is vacuously true — there is nothing
//! left for a second member to collide with. Two very different things reach
//! that state:
//!
//! * members **shifted apart** while still plural — the conflict was erased and
//!   the description launders; this is what the gate is for;
//! * members **composed into one edit** — the colliding writes were resolved
//!   into a single description denoting a definite sequence, which is what
//!   `delins.md:86-89` asks for and what the pure-deletion exemption in
//!   `overlap.rs` already relies on.
//!
//! Without the plural term the gate reverts #1423, whose whole point is that
//! `g.[11_12inv;11_12insAA]` collapses to the single member `g.10_11A[4]`.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::normalize::NormalizeConfig;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// 20 bases of `A`/`T`, the corpus sequence the `[inv;ins]` families settle on.
const AT_SEQUENCE: &str = "TTTTTTTTTAATATATTTTA";

/// 18 bases; `10` and `11` are `A`, `12` is `T`. #1416/#1423's reference.
const CG_SEQUENCE: &str = "CGCGCGCGCAATCGCGCG";

fn provider(sequence: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", sequence.to_string());
    provider
}

fn lenient(sequence: &str, input: &str) -> String {
    Normalizer::new(provider(sequence))
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("lenient normalization must not reject")
        .to_string()
}

/// Strict mode rejects `input` **as an overlap conflict**.
///
/// Keyed on the diagnostic, not on `is_err()`. A bare error check would count
/// any unrelated failure — a stated-reference mismatch, an out-of-bounds
/// coordinate — as a successful rejection, so the laundering assertion below
/// could pass while the conflict really had been laundered and the output was
/// refused for some other reason entirely.
fn strict_rejects_as_conflict(sequence: &str, input: &str) -> bool {
    match Normalizer::with_config(
        provider(sequence),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&parse_hgvs(input).expect("parse"))
    {
        Ok(_) => false,
        Err(error) => {
            let message = error.to_string();
            assert!(
                message.contains("W5002") || message.contains("OverlapConflictingEdits"),
                "`{input}` was rejected, but not as an overlap conflict; got: {message}"
            );
            true
        }
    }
}

/// The headline: lenient mode must not turn a rejected description into an
/// accepted one.
///
/// Asserted as the property rather than as a pinned string, because that is the
/// contract — whatever the output is, strict mode must judge it the same way it
/// judged the input.
#[test]
fn the_lenient_output_of_a_conflict_is_still_rejected_by_strict() {
    for input in [
        "TEMPLATE:g.[9_10insA;9_10inv]",
        "TEMPLATE:g.[11_12dup;12_13inv]",
        "TEMPLATE:g.[10_11inv;10dup]",
    ] {
        assert!(
            strict_rejects_as_conflict(AT_SEQUENCE, input),
            "`{input}` must be an overlap conflict for this test to mean anything"
        );
        let out = lenient(AT_SEQUENCE, input);
        assert!(
            strict_rejects_as_conflict(AT_SEQUENCE, &out),
            "lenient mode turned `{input}` into `{out}`, which strict mode accepts — \
             the conflict was laundered out of the description"
        );
    }
}

/// The preserved output is the input, so it is a fixed point in one pass.
#[test]
fn a_preserved_conflict_is_a_fixed_point() {
    for input in [
        "TEMPLATE:g.[9_10insA;9_10inv]",
        "TEMPLATE:g.[11_12dup;12_13inv]",
    ] {
        let once = lenient(AT_SEQUENCE, input);
        assert_eq!(
            once, input,
            "a conflicting allele must come back as authored"
        );
        assert_eq!(
            lenient(AT_SEQUENCE, &once),
            once,
            "and must settle in one pass"
        );
    }
}

/// An **uncertain** cis allele launders the same way, so the gate must reach it
/// too.
///
/// This is the one place the preserve gate deliberately differs from the sort
/// gate above it. That gate excludes `[(a;b)]` because member order inside an
/// uncertain group is authored presentation and not ours to reorder. Whether the
/// conflict survives into the output is a different question, and the answer
/// does not depend on the wrapper — measured before the gate was widened:
///
/// ```text
/// in  g.[(9_10insA;9_10inv)]   strict: REJECT (W5002)
/// out g.[(11dup;9_10=)]        strict: ACCEPT
/// ```
#[test]
fn an_uncertain_conflict_is_preserved_too() {
    for input in [
        "TEMPLATE:g.[(9_10insA;9_10inv)]",
        "TEMPLATE:g.[(11_12dup;12_13inv)]",
    ] {
        assert!(
            strict_rejects_as_conflict(AT_SEQUENCE, input),
            "`{input}` must be an overlap conflict for this test to mean anything"
        );
        let out = lenient(AT_SEQUENCE, input);
        assert!(
            strict_rejects_as_conflict(AT_SEQUENCE, &out),
            "the uncertain spelling laundered: `{input}` became `{out}`, which strict accepts"
        );
        assert!(
            out.contains("[("),
            "the uncertainty wrapper must survive preservation; got `{out}`"
        );
    }
}

/// The discriminating half: a conflict resolved by **composition into a single
/// member** is not preserved, because there is no laundering — the colliding
/// writes were merged into one description that denotes a definite sequence.
///
/// This is #1423's case, and it is the guard that stops the gate being written
/// as "the output no longer conflicts" alone. That predicate is vacuously true
/// for any single-member output, so without the plural term this input would be
/// handed back unchanged and #1423 would be undone.
#[test]
fn a_conflict_resolved_into_one_member_is_not_preserved() {
    let input = "TEMPLATE:g.[11_12inv;11_12insAA]";
    assert!(
        strict_rejects_as_conflict(CG_SEQUENCE, input),
        "the input must be a conflict for this test to discriminate"
    );
    assert_eq!(
        lenient(CG_SEQUENCE, input),
        "TEMPLATE:g.10_11A[4]",
        "a conflict composed into a single member must keep #1423's collapse, \
         not be handed back as authored"
    );
}

/// A junction interior to a **pure deletion** was never a conflict, so nothing
/// here is preserved and the merged form the spec asks for survives.
///
/// A deletion removes every base it spans, so an interior junction has nothing
/// left to be positioned against: the pair denotes the same bases in either
/// order and whichever interior junction the insertion named.
///
/// Deliberately on `CG_SEQUENCE`, where `2_3` is `GC`. On the `A`/`T` reference
/// the same input settles as `g.2_3inv` instead — `AA` is the reverse
/// complement of `TT`, so the merged form is a textbook inversion
/// (`inversion.md:5`) and the assertion would be testing the type dispatcher
/// rather than the exemption. `revcomp("GC")` is `GC`, so no payload spelled
/// `AA` can be read as an inversion there.
#[test]
fn a_junction_interior_to_a_deletion_still_merges() {
    let input = "TEMPLATE:g.[2_3del;2_3insAA]";
    assert!(
        !strict_rejects_as_conflict(CG_SEQUENCE, input),
        "a junction interior to a pure deletion composes uniquely, so it is not a conflict"
    );
    assert_eq!(
        lenient(CG_SEQUENCE, input),
        "TEMPLATE:g.2_3delinsAA",
        "`delins.md:86-89` asks for the merged form"
    );
}
