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
//! ## The companion row, and what it does *not* guard
//!
//! `a_non_conflicting_allele_of_the_same_shape_still_merges` pins that the
//! *same* `dup`-plus-substitution shape one base to the left — where the
//! duplication has a junction to be respelled at and so never becomes a
//! conflict — still collapses to one member. That is a real regression guard for
//! the merge outcome and worth keeping.
//!
//! **It is not a control on the gate's scoping**, though this header and that
//! test's own docstring both used to claim it was (#1435). The claim was that a
//! gate keyed on the `dup`-plus-substitution *shape* rather than on the conflict
//! would refuse it, and the test would notice. It would not: `sequence_first_pass`
//! runs *after* the legacy per-member collapse, so `g.[23dup;23A>G]` has already
//! become the single member `g.22_23insG` by the time the gate is reached, and a
//! shape-keyed gate would have nothing left to key on.
//!
//! Measured rather than argued, on this file's own fixture: instrumenting the
//! gate reports `sequence_first_pass sees 1 member` for that input, so its
//! `members.len() > 1` precondition is false and the gate never runs — and
//! disabling the gate outright leaves the test passing unchanged.
//!
//! Left corrected rather than deleted because the wrong version was load-bearing
//! in the wrong direction: it told a reader that gate scoping was already
//! covered, which is exactly what stops someone writing the control that would
//! cover it. **A genuine control needs an input that reaches
//! `sequence_first_pass` with two members still intact.**

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
/// `(accession, input)`, one per distinguishable route into the gate: two
/// substitutions of one base, a coincident deletion and inversion, and an
/// insertion interior to an inversion (#1276's shape). Three edit-type pairings,
/// so the property is not pinned on one.
///
/// **This list used to be four `dup` rows plus the inversion**, and #1406
/// removed all four — a `dup` writes at the junction 3' of the span it names
/// (`duplication.md:5`), so it never claims the base a coincident sibling
/// claims. The assertion below anticipated exactly this ("either it stopped
/// being a conflict — a result worth its own measurement — or the row is
/// stale"); it was the first, and the measurement is #1406's. Their merged
/// outputs are pinned in `issue_1307_terminal_dup_respell` and
/// `issue_1406_conflict_is_a_property_of_the_input` instead.
///
/// What replaced them are conflicts by write footprint rather than by read span,
/// in the two shapes the gate actually distinguishes. The first two rows are
/// coincident *bases*: both members write the same span, so the result depends
/// on which one wins. The third is the other route entirely — an insertion
/// junction lying inside an inversion span (#1276), which
/// `detect_insertion_overlaps` reports rather than `detect_overlap_conflicts`.
/// Keeping both means the property is not pinned on one detector.
const CONFLICTING_INPUTS: &[(&str, &str)] = &[
    // Two substitutions of base 23 — the least arguable conflict there is.
    ("NC_TEST.1", "NC_TEST.1:g.[23A>G;23A>T]"),
    // Coincident `del` and `inv` over one span. One of the rows #1406 measured
    // as silently laundered: it normalized to a single accepted `g.15_16del`,
    // dropping the inversion outright.
    ("NC_TEST.1", "NC_TEST.1:g.[15_16del;15_16inv]"),
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

/// A non-conflicting `dup`-plus-substitution pair still collapses to one member.
///
/// `g.[23dup;23A>G]` is the identical pair one base to the left of #1307's.
/// There it is *not* a conflict — the duplication has a junction 3' of itself to
/// be respelled at, so the two members never end up coincident — and the
/// **legacy per-member collapse** merges it into a single insertion.
///
/// **This does not guard `sequence_first_pass`'s gate (#1435.)** It used to say
/// it did — "a gate that refused this would be keyed on the shape rather than on
/// the conflict" — and that is not reachable here. The collapse runs first, so
/// the gate is handed one member, its `members.len() > 1` precondition is false,
/// and it never runs: instrumenting it reports `sees 1 member` for this input,
/// and disabling the gate outright leaves this test passing unchanged. A gate
/// erroneously keyed on the shape would have nothing to key on.
///
/// What it does guard is the merge outcome, which is worth keeping on its own
/// terms. See the module header for what a genuine control on the gate's
/// scoping would need.
#[test]
fn a_non_conflicting_allele_of_the_same_shape_still_merges() {
    let input = "NC_TEST.1:g.[23dup;23A>G]";
    assert!(
        !strict_rejects("NC_TEST.1", input),
        "`{input}` carries no conflict — the duplication writes at its own \
         junction, clear of the substitution — so strict mode must accept it"
    );
    assert_eq!(
        normalize_lenient("NC_TEST.1", input),
        "NC_TEST.1:g.22_23insG",
        "`{input}` must still collapse to one member through the legacy \
         per-member collapse"
    );
}
