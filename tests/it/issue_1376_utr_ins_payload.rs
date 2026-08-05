//! Issue #1376 — a UTR-marker `ins[...]` payload leaked a raw `nom` error.
//!
//! `c.249_250ins[*100_*200]` failed with
//!
//! ```text
//! Failed to parse variant: Error(Error { input: "ins[*100_*200]", code: Tag })
//! ```
//!
//! a combinator's internal state reaching the user, naming neither the construct
//! nor the reason. Its intronic sibling `ins[244-8_249]` parsed and was refused
//! downstream with a sentence that **already names this shape**: "CDS-offset
//! range (intronic *or UTR-marker*) is spec-undefined and not yet supported".
//!
//! So the message existed and the shape simply never reached it:
//! `parse_cds_position_range` required a leading digit, and stopped at the `*`.
//! Admitting the `*` routes the range to `InsertedPart::CdsPositionRange`, which
//! is the arm that carries the refusal.
//!
//! What this does **not** change: whether such a payload is expandable. #466's
//! candidate 2 asks SVD-WG to define the canonical expansion for these shapes,
//! and stays open either way. This is only about the diagnostic.

use ferro_hgvs::hgvs::edit::{InsertedPart, InsertedSequence, NaEdit};
use ferro_hgvs::hgvs::uncertainty::Mu;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;

/// The `InsertedPart`s a parsed `ins[...]` payload decomposed into.
///
/// Asserting on the parsed shape rather than on a normalization error message is
/// deliberate: the routing to `CdsPositionRange` *is* the fix, and it is what
/// determines which refusal a caller gets. Reaching the message itself needs a
/// provider that resolves the transcript's CDS, so it is verified end to end
/// against a real prepared reference (see the PR) rather than pinned against a
/// mock whose plumbing could change the answer for unrelated reasons.
fn inserted_parts(descriptor: &str) -> Vec<InsertedPart> {
    let variant = parse_hgvs(descriptor)
        .unwrap_or_else(|e| panic!("`{descriptor}` must now parse, but: {e}"));
    let HgvsVariant::Cds(cds) = variant else {
        panic!("`{descriptor}` must parse as a c. variant");
    };
    let Mu::Certain(NaEdit::Insertion { sequence }) = cds.loc_edit.edit else {
        panic!("`{descriptor}` must parse as a certain insertion");
    };
    match sequence {
        InsertedSequence::Complex(parts) => parts,
        other => panic!("`{descriptor}` must carry a complex payload, got {other:?}"),
    }
}

/// The UTR-marker payload parses, and routes to the arm that carries the curated
/// refusal — instead of dying in the grammar with a `nom` error.
#[test]
fn a_utr_marker_range_payload_routes_to_the_curated_refusal() {
    let parts = inserted_parts("NM_000088.3:c.249_250ins[*100_*200]");
    assert_eq!(parts.len(), 1, "one payload part, got {parts:?}");
    assert!(
        matches!(&parts[0], InsertedPart::CdsPositionRange(range) if range == "*100_*200"),
        "must route to CdsPositionRange, whose refusal already names UTR markers; \
         got {parts:?}"
    );
}

/// The intronic sibling is unchanged — it is the behaviour being matched, so a
/// change here would mean the two drifted apart again rather than converged.
#[test]
fn the_intronic_sibling_is_unchanged() {
    let parts = inserted_parts("NM_000088.3:c.249_250ins[244-8_249]");
    assert!(
        matches!(&parts[0], InsertedPart::CdsPositionRange(range) if range == "244-8_249"),
        "got {parts:?}"
    );
}

/// A plain positive-integer range must still be a plain same-reference range,
/// not routed to the refusal — `*` widened the grammar, and this pins that it
/// widened it only where intended.
#[test]
fn a_plain_range_payload_still_resolves() {
    let parts = inserted_parts("NM_000088.3:c.249_250ins[401_419]");
    assert!(
        matches!(
            &parts[0],
            InsertedPart::PositionRange {
                start: 401,
                end: 419
            }
        ),
        "a plain range must stay a plain same-reference range, got {parts:?}"
    );
    let variant = parse_hgvs("NM_000088.3:c.249_250ins[401_419]").expect("must parse");
    assert_eq!(
        variant.to_string(),
        "NM_000088.3:c.249_250ins401_419",
        "a plain range must keep its flat same-reference form"
    );
}

/// The rendered form must itself parse.
///
/// A second regression this change had to avoid, and did not at first: admitting
/// `ins[*100_*200]` in the *bracketed* grammar alone made it parse into something
/// whose own rendering — `ins*100_*200`, brackets dropped for a single part — was
/// still rejected. Parse-then-display then produced a string ferro could not read
/// back, which is the asymmetry `FERRO_ASSERT_REPARSE` exists to catch for
/// normalization and which nothing catches at the parse seam.
///
/// The intronic sibling never had it, because `ins244-8_249` starts with a digit
/// and so was already accepted unbracketed. Admitting `*` on that path too makes
/// the pair symmetric again.
#[test]
fn the_rendered_form_parses_back() {
    for descriptor in [
        "NM_000088.3:c.249_250ins[*100_*200]",
        // The unbracketed form directly, which is what the one above renders as.
        "NM_000088.3:c.249_250ins*100_*200",
    ] {
        let variant = parse_hgvs(descriptor).unwrap_or_else(|e| panic!("`{descriptor}`: {e}"));
        let rendered = variant.to_string();
        let reparsed = parse_hgvs(&rendered).unwrap_or_else(|e| {
            panic!("`{descriptor}` rendered `{rendered}`, which will not parse: {e}")
        });
        assert_eq!(
            reparsed.to_string(),
            rendered,
            "`{descriptor}` must reach a display fixed point"
        );
    }
}

/// Widening the unbracketed path must not disturb the other things that can
/// follow `ins`: a literal sequence, a bare count, a repeat, or an intronic
/// range. `*` was added to the dispatch arm those share.
#[test]
fn the_other_unbracketed_payload_shapes_are_unchanged() {
    for descriptor in [
        "NM_000088.3:c.249_250insATG",
        "NM_000088.3:c.249_250ins10",
        "NM_000088.3:c.249_250insN[15]",
        "NM_000088.3:c.249_250ins244-8_249",
        "NM_000088.3:c.249_250ins401_419",
    ] {
        let variant =
            parse_hgvs(descriptor).unwrap_or_else(|e| panic!("`{descriptor}` must parse: {e}"));
        assert_eq!(
            variant.to_string(),
            descriptor,
            "`{descriptor}` must round-trip unchanged"
        );
    }
}

/// The `*` is a marker on a position, not a position. A bare one, or an empty
/// token, must still be a parse error — otherwise admitting `*` would have made
/// the grammar accept nonsense.
#[test]
fn a_marker_without_a_position_is_still_rejected() {
    for descriptor in [
        "NM_000088.3:c.249_250ins[*_*200]",
        "NM_000088.3:c.249_250ins[*100_*]",
        "NM_000088.3:c.249_250ins[*100_]",
        "NM_000088.3:c.249_250ins[*]",
    ] {
        assert!(
            parse_hgvs(descriptor).is_err(),
            "`{descriptor}` names no second position and must not parse"
        );
    }
}

/// The mixed form the issue's real-world example uses — an intronic range beside
/// `N[k]` padding — routes the same way, and so does its UTR twin.
#[test]
fn a_padded_multi_part_payload_routes_the_same_way() {
    for (descriptor, expected) in [
        ("NM_000088.3:c.249_250ins[N[2800];244-8_249]", "244-8_249"),
        ("NM_000088.3:c.249_250ins[N[2800];*100_*200]", "*100_*200"),
    ] {
        let parts = inserted_parts(descriptor);
        assert_eq!(parts.len(), 2, "`{descriptor}` -> {parts:?}");
        assert!(
            matches!(&parts[0], InsertedPart::Repeat { .. }),
            "`{descriptor}`: the `N[2800]` padding must survive, got {parts:?}"
        );
        assert!(
            matches!(&parts[1], InsertedPart::CdsPositionRange(range) if range == expected),
            "`{descriptor}` -> {parts:?}"
        );
    }
}

/// Both positions of the range must admit exactly the same shapes.
///
/// A c. position carries at most ONE intronic offset, so `100+5-3` is not a
/// position at all. The refactor that admitted `*` extracted the position scan
/// into a helper that consumes the offset itself, and for a while the caller
/// still ran its own offset pass afterwards — but only for the FIRST position.
/// That made the two ends of one range disagree: `ins[100+5-3_200]` parsed while
/// `ins[100_200+5-3]` did not, so the grammar accepted a position the spec has no
/// reading for. Green suites do not catch this; nothing asserted the pair.
///
/// Asserted as a symmetry rather than as two independent cases, so a future
/// change that loosens one end has to loosen the other deliberately.
#[test]
fn both_ends_of_a_range_accept_the_same_position_shapes() {
    // Each row carries the expected verdict as well as the pair, because
    // symmetry alone is satisfied by BOTH ends being wrong together — a change
    // that rejected every offset would keep `lhs == rhs` and pass.
    for (first, second, accepted) in [
        // One offset per position is legal at either end.
        ("100+5_200", "100_200+5", true),
        ("100-5_200", "100_200-5", true),
        // Two offsets on one position is not a position, at either end.
        ("100+5-3_200", "100_200+5-3", false),
        ("100-5+3_200", "100_200-5+3", false),
        // An offset marker with no magnitude is not a position either.
        ("100+_200", "100_200+", false),
        ("100-_200", "100_200-", false),
        // UTR markers, the shape this file exists for.
        ("*100_*200", "*100_*200", true),
        ("*100-_*200", "*100_*200-", false),
    ] {
        let lhs = ferro_hgvs::parse_hgvs(&format!("NM_000088.3:c.249_250ins[{first}]")).is_ok();
        let rhs = ferro_hgvs::parse_hgvs(&format!("NM_000088.3:c.249_250ins[{second}]")).is_ok();
        assert_eq!(
            lhs,
            accepted,
            "`ins[{first}]` should {} parse",
            if accepted { "" } else { "not" }
        );
        assert_eq!(
            rhs,
            accepted,
            "`ins[{second}]` should {} parse",
            if accepted { "" } else { "not" }
        );
        assert_eq!(
            lhs, rhs,
            "`ins[{first}]` and `ins[{second}]` disagree ({lhs} vs {rhs}) — the two \
             positions of one range must admit the same shapes"
        );
    }

    // And the specific regression: a doubled offset is refused at BOTH ends.
    for descriptor in [
        "NM_000088.3:c.249_250ins[100+5-3_200]",
        "NM_000088.3:c.249_250ins[100_200+5-3]",
    ] {
        assert!(
            ferro_hgvs::parse_hgvs(descriptor).is_err(),
            "`{descriptor}` carries two offsets on one position and must not parse"
        );
    }
}
