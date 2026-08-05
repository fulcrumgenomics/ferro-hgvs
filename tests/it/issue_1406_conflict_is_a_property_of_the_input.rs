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
    // The heart of #1406. Neither position may be privileged: both are a `dup`
    // of one base beside a substitution of that base.
    for input in ["TEMPLATE:g.[23dup;23A>G]", "TEMPLATE:g.[24dup;24C>G]"] {
        assert!(
            strict_rejects(input),
            "`{input}` is an overlap conflict and strict mode must reject it"
        );
        assert_eq!(
            lenient(input),
            input,
            "a conflicting allele is preserved as authored, not canonicalized"
        );
    }
}

#[test]
fn the_conflicting_allele_is_a_fixed_point() {
    // Preserving is only an answer if it is stable: re-normalizing the output
    // must not start the repair the first pass declined.
    let input = "TEMPLATE:g.[23dup;23A>G]";
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
    assert_ne!(
        lenient(input),
        input,
        "a non-conflicting allele must still be canonicalized"
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
