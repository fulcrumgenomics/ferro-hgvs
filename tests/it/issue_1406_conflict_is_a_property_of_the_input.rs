//! W5002 must describe the input, not whether the repair happened to succeed.
//!
//! `detect_overlap_conflicts` ran only *post-shift*, on the normalized members.
//! `detect_insertion_overlaps` already ran on the raw ones — deliberately, so a
//! merge could not hide an overlap — but its coincident-bounds sibling did not,
//! and that gap decided verdicts:
//!
//! ```text
//! reference  TTTTTTTTTAATATATTTTAATAC     24 bases, 23 = `A`, 24 = `C`
//!
//! g.[23dup;23A>G]   ->  g.22_23insG    no warning at all, strict ACCEPTS
//! g.[24dup;24C>G]   ->  unchanged      W5002, strict REJECTS
//! ```
//!
//! Both are a `dup` of one base beside a substitution of that same base — the
//! same shape, one position apart. At 23 the repair collapses the pair into a
//! single member, so by the time the post-shift detector looked there was
//! nothing left to conflict; at 24 the terminal decline (#1307) left the pair
//! intact and it was reported. Proximity to the end of the contig decided
//! whether a conflicting description was rejected, which is not a property of
//! the description.
//!
//! Running the coincident-bounds detector on the raw members too takes the
//! verdict from where the conflict lives, and a conflicting allele is preserved
//! as authored — ferro's standing answer since #395/#486/#1004: report the
//! conflict and leave the description alone rather than pick a winner among
//! orderings the spec does not rank.
//!
//! ## Deliberately not widened
//!
//! The *other* detector's conflicts — an insertion interior to a `del`/`delins`/
//! `inv` span — keep their existing treatment: each member is still normalized
//! on its own, so `c.[4_10inv;5_6insAA]` settles as `c.[5_9inv;6_7insAA]`.
//! Preserving those verbatim would drop the per-member 3' rule that #1276 and
//! #1235's transcript axes pin on purpose. Measured: widening the preserve to
//! them fails exactly those three tests.
//!
//! That leaves one of #1406's three rows open, and it is pinned below rather
//! than left implicit.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::normalize::NormalizeConfig;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// 24 bases; position 23 is `A` and 24 is the final `C`. The #1307 reference.
const SEQUENCE: &str = "TTTTTTTTTAATATATTTTAATAC";

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", SEQUENCE.to_string());
    provider
}

fn lenient(input: &str) -> String {
    Normalizer::new(provider())
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("lenient normalization must not reject")
        .to_string()
}

fn strict_rejects(input: &str) -> bool {
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&parse_hgvs(input).expect("parse"))
    .is_err()
}

#[test]
fn the_same_shape_gets_the_same_verdict_wherever_it_sits() {
    // The heart of #1406, and the verdict is ACCEPT at both positions.
    //
    // Neither position may be privileged: both are a `dup` of one base beside a
    // substitution of that base. The defect was that 23 merged and 24 was
    // rejected, decided only by proximity to the contig end.
    //
    // Which way the parity resolves is the other half, and it is the spec's
    // call rather than a preference. A `dup` writes at the junction 3' of the
    // span it names (`duplication.md:5`), not over it, so a `dup` and a
    // substitution of the same base have disjoint write footprints, compose
    // uniquely, and were never a conflict. `delins.md:86-89` asks for the merged
    // form outright. So the pair converges by merging, not by being refused.
    //
    // The merged forms are pinned literally, not just compared to each other:
    // two spellings settling on one *wrong* string would satisfy a parity check
    // silently. One base apart, the outputs are the same shape one base apart.
    for (input, expected) in [
        ("TEMPLATE:g.[23dup;23A>G]", "TEMPLATE:g.22_23insG"),
        ("TEMPLATE:g.[24dup;24C>G]", "TEMPLATE:g.23_24insG"),
    ] {
        assert!(
            !strict_rejects(input),
            "`{input}` is not an overlap conflict — its members write to \
             different places — so strict mode must accept it"
        );
        assert_eq!(
            lenient(input),
            expected,
            "`{input}` must reach the merged form the spec asks for"
        );
    }
}

#[test]
fn the_conflicting_allele_is_a_fixed_point() {
    // Preserving is only an answer if it is stable: re-normalizing the output
    // must not start the repair the first pass declined.
    //
    // The allele is two substitutions of ONE base. That is a genuine conflict —
    // both members write base 23 and the result depends on which wins — unlike
    // the `dup` pair above, whose members write to different places. This test
    // used that pair until #1406 established it was never a conflict; keeping it
    // here would have pinned the misclassification rather than the preserve.
    let input = "TEMPLATE:g.[23A>G;23A>T]";
    // Assert the premise first. Preservation is only the *right* answer for an
    // allele that genuinely conflicts, so a row that quietly stopped being one
    // would otherwise satisfy everything below while pinning nothing — the same
    // trap #1406 found in the `dup` rows this test used to use.
    assert!(
        strict_rejects(input),
        "`{input}` must still be a conflict strict mode rejects"
    );
    let once = lenient(input);
    // The form that must be stable is the *preserved* one. Asserting only
    // `lenient(once) == once` would pass on the repaired output too, which is
    // also a fixed point — the check would then guard nothing.
    assert_eq!(
        once, input,
        "the conflicting allele must be what is preserved"
    );
    assert_eq!(lenient(&once), once, "`{once}` must be a fixed point");
}

#[test]
fn a_disjoint_dup_and_substitution_still_normalize() {
    // The guard. The verdict must key on the bases actually shared, not on the
    // presence of a `dup` beside a substitution — #999's `g.[306dup;308C>A]` is
    // the shape that must keep flowing through.
    let input = "TEMPLATE:g.[23dup;24C>G]";
    assert!(
        !strict_rejects(input),
        "`{input}` shares no base and must not be reported as a conflict"
    );
    // The exact form, not merely "something other than the input". `assert_ne!`
    // passes on any change at all, including a wrong one, which is no guard for
    // a test whose whole point is that this allele takes the normal path.
    assert_eq!(
        lenient(input),
        "TEMPLATE:g.24delinsAG",
        "a non-conflicting allele must reach its canonical merged form"
    );
}

// The row of #1406 this change does **not** close.
//
// `g.[11_12inv;11_12insAA]` (on `CGCGCGCGCAATCGCGCG`) is strict-rejected, and
// its own lenient output `g.[11_12=;10_11A[4]]` is strict-*accepted* — lenient
// mode converts a description strict mode refuses into one it admits.
//
// It is not reachable from here. The laundering is produced by *per-member*
// normalization (the inversion cancels to `=`, the insertion respells as a
// repeat), not by the whole-allele canonicalization this change gates, and
// those per-member results are exactly what #1276 and #1235's transcript axes
// require for a conflicting allele. Closing it means stopping a member's own
// respelling from erasing the conflict, which is a different mechanism.
//
// **There is no executable pin for it here, deliberately.** Normalizing that
// input trips `FERRO_ASSERT_IDEMPOTENT`: the output `g.[11_12=;10_11A[4]]`
// re-normalizes to `g.[10_11A[4];11_12=]`, the same two members in genomic
// order rather than the order the first pass emitted. So the first output also
// carries out-of-order members, which is #1235's criterion 2. That is a
// *second*, separate defect on the same input, verified pre-existing by A/B
// against `origin/main`, and a test that normalized the input would redden the
// `test-oracle` job rather than pin anything. Both rows are recorded on #1406.
