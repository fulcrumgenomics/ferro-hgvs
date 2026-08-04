//! Lenient normalization must not launder away an overlap conflict.
//!
//! The property, stated once: **if strict mode rejects an input because the spec
//! defines no canonical form for it, strict mode must still reject whatever
//! lenient mode normalizes it into.** Lenient mode is allowed to re-spell a
//! conflicting allele; it is not allowed to resolve one. An output that strict
//! mode accepts is a claim that the description is well-formed, and for these
//! inputs that claim is false.
//!
//! This is a property of the normalizer, not of any one reproducer. It is pinned
//! here rather than only inside `issue_1307_terminal_dup_respell` because the
//! gate that enforces it — `Normalizer::sequence_first_pass`'s
//! `detect_overlap_conflicts` check — is general, and the next input to reach it
//! will not be #1307's. Individual shapes keep their own files; what lives here
//! is the invariant they are all instances of.
//!
//! ## Why the gate is needed at all
//!
//! Two refusals in the normalizer answer overlapping-looking questions in
//! different coordinate spaces, and they disagree:
//!
//! - `merge::apply_edits_to_window` refuses on **interbase** geometry. `24dup`
//!   occupies the zero-length junction `[24, 24]`; `24C>G` occupies the span
//!   `[23, 24]`. Flush, not overlapping — so it admits the pair, and so does the
//!   test-suite applier for the same reason.
//! - `normalize::overlap::detect_overlap_conflicts` works in **HGVS-coordinate**
//!   space, where both members name base 24. It reports the conflict, and it is
//!   this verdict that strict mode raises as `OverlapConflictingEdits / W5002`.
//!
//! With only the geometric refusal in place, the sequence-first derivation was
//! free to re-derive `g.[24dup;24C>G]` as `g.23_24insG` — in range, denoting the
//! same bases, and strict-**acceptable**. That is #1307, and it is the shape the
//! rows below start from.
//!
//! ## The control
//!
//! `a_non_conflicting_allele_of_the_same_shape_still_merges` is what makes the
//! rows above mean "no laundering" rather than "multi-member alleles stay
//! multi-member". It is the *same* `dup`-plus-substitution shape one base to the
//! left, where the duplication has a junction to be respelled at and so never
//! becomes a conflict — and it must still collapse to one member. A gate keyed
//! on the shape rather than on the conflict would refuse it too, and this test
//! is what would notice.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::normalize::NormalizeConfig;
use ferro_hgvs::{parse_hgvs, FerroError, MockProvider, Normalizer};

/// 24 bases ending in a unique `C`, so a member authored at position 24 stays
/// there and 24 is genuinely the last base. Shared with
/// `issue_1307_terminal_dup_respell`.
const SEQUENCE: &str = "TTTTTTTTTAATATATTTTAATAC";

fn provider(accession: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(accession, SEQUENCE.to_string());
    provider
}

/// Normalize in the default lenient mode.
fn normalize_lenient(accession: &str, input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::new(provider(accession))
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string()
}

/// `true` when strict mode rejects `description` **specifically as
/// `OverlapConflictingEdits / W5002`**.
///
/// Panics on any other error, rather than returning `true` for it. `is_err()`
/// alone cannot tell the contract under test from an unrelated rejection — a
/// past-end position or a reference mismatch would satisfy it just as well —
/// and this file's whole claim is about *that* gate surviving normalization.
///
/// Both halves of the tag are required because the promotion site emits them
/// together (`… (OverlapConflictingEdits / W5002)`); accepting either alone
/// would keep passing after the other half regressed. Same shape as
/// `strict_rejects_as_overlap` in `issue_1276_dup_junction_overlap.rs`.
fn strict_rejects(accession: &str, description: &str) -> bool {
    let variant = parse_hgvs(description).expect("description must parse");
    match Normalizer::with_config(
        provider(accession),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    {
        Ok(_) => false,
        Err(FerroError::InvalidCoordinates { msg }) => {
            assert!(
                msg.contains("W5002") && msg.contains("OverlapConflictingEdits"),
                "expected `OverlapConflictingEdits / W5002` for `{description}`; got: {msg}"
            );
            true
        }
        Err(other) => panic!("unexpected error variant for `{description}`: {other:?}"),
    }
}

/// Inputs strict mode rejects as `OverlapConflictingEdits / W5002`, one per
/// distinguishable route into the gate.
///
/// `(accession, input)`. The first three are the #1307 terminal `dup` collision
/// in both member orders and with a `delins` sibling; the fourth is its circular
/// twin, which reaches the same gate because `m.` has a last base like any other
/// sequence; the fifth is a different conflict class entirely — an insertion
/// interior to an inversion (#1276's shape) — so the property is not pinned on
/// one edit-type pair.
const CONFLICTING_INPUTS: &[(&str, &str)] = &[
    ("NC_TEST.1", "NC_TEST.1:g.[24dup;24C>G]"),
    ("NC_TEST.1", "NC_TEST.1:g.[24C>G;24dup]"),
    ("NC_TEST.1", "NC_TEST.1:g.[24dup;24delinsGG]"),
    ("NC_012920.1", "NC_012920.1:m.[24dup;24C>G]"),
    ("NC_TEST.1", "NC_TEST.1:g.[9_12inv;10_11insAA]"),
];

#[test]
fn a_conflict_strict_rejects_survives_lenient_normalization() {
    for (accession, input) in CONFLICTING_INPUTS {
        // Assert the premise, so a row that quietly stops being a conflict fails
        // loudly instead of satisfying the implication vacuously. This is the
        // half that would otherwise rot: if a future change makes one of these
        // well-formed, the row below becomes trivially true and goes on looking
        // like coverage.
        assert!(
            strict_rejects(accession, input),
            "`{input}` is listed here as an overlap conflict but strict mode now \
             accepts it, so this row pins nothing. Either it stopped being a \
             conflict — a result worth its own measurement — or the row is stale."
        );

        let lenient = normalize_lenient(accession, input);
        assert!(
            strict_rejects(accession, &lenient),
            "`{input}` -> `{lenient}`, which strict mode accepts. Lenient mode \
             may re-spell a conflicting allele but must not resolve one: strict \
             mode declares these have no canonical form, so an output it accepts \
             has laundered the conflict into a description that looks well-formed."
        );
    }
}

/// The gate must key on the conflict, not on the shape it was found in.
///
/// `g.[23dup;23A>G]` is the identical `dup`-plus-substitution pair one base to
/// the left of #1307's. There it is *not* a conflict — the duplication has a
/// junction 3' of itself to be respelled at, so the two members never end up
/// coincident — and the derivation is expected to merge it into a single
/// insertion exactly as it did before the gate existed.
#[test]
fn a_non_conflicting_allele_of_the_same_shape_still_merges() {
    let input = "NC_TEST.1:g.[23dup;23A>G]";
    assert!(
        !strict_rejects("NC_TEST.1", input),
        "`{input}` is the control precisely because it is *not* a conflict; if \
         strict mode has started rejecting it, it can no longer show that the \
         gate is scoped to conflicts rather than to the shape"
    );
    assert_eq!(
        normalize_lenient("NC_TEST.1", input),
        "NC_TEST.1:g.22_23insG",
        "`{input}` carries no conflict, so the derivation must still merge it \
         into one member — a gate that refused this would be keyed on the \
         `dup`-plus-substitution shape rather than on the conflict"
    );
}
