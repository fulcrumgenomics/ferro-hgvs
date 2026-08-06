//! #1492 — a repeat decline must name coordinates the caller supplied.
//!
//! When a single-position repeat anchor names no tandem run, the tract search
//! falls back to the unit-wide span *at* the anchor, and the decline quoted
//! that span:
//!
//! ```text
//! g.260CAG[5]  ->  repeat span TEMPLATE:260-262 does not match repeat unit CAG
//! ```
//!
//! The caller wrote one position, `260`. The range `260-262` is the search's
//! own invention, so the message sent the reader to inspect a three-base window
//! they never named — and said nothing about the actual fault, which is that
//! the run begins at `259` and `260` is out of phase with it.
//!
//! Not a corner: **336 of the 560 rows** in #1452's census take this path — it
//! is every anchor that does not land on a unit boundary — and it behaves the
//! same on masked and unmasked references.
//!
//! The explicit-range spelling was already correct and must stay so, because
//! there the span *is* the caller's own. That asymmetry is the whole point, and
//! it is why the fix threads a `RepeatSpanOrigin` out of the search rather than
//! rewording the message unconditionally: a decline over a found tract, or over
//! a range the caller wrote, should still quote its span.
//!
//! Behaviour is unchanged — the same inputs convert and the same inputs are
//! refused. Only the text moves.

use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, MockProvider};

const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

fn padded(core: &str) -> String {
    format!("{PAD}{core}{PAD}")
}

fn spdi_err(sequence: &str, input: &str) -> String {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", sequence.to_string());
    let variant = parse_hgvs(input).expect("input must parse");
    match hgvs_to_spdi(&variant, &provider) {
        Ok(triple) => panic!("`{input}` must not convert, got {triple}"),
        Err(e) => e.to_string(),
    }
}

/// A 3-copy `CAG` tract at 259-267.
const TRACT: &str = "GGCAGCAGCAGGG";

/// The headline: an out-of-phase anchor is told about the position it wrote.
#[test]
fn an_out_of_phase_anchor_decline_names_the_anchor_not_a_fabricated_span() {
    let err = spdi_err(&padded(TRACT), "TEMPLATE:g.260CAG[5]");
    assert!(
        err.contains("no CAG repeat is anchored at TEMPLATE:260"),
        "decline must name the caller's anchor: {err}"
    );
    assert!(
        err.contains("g.<start>_<end>CAG[n]"),
        "decline must offer the explicit-range spelling as the way out: {err}"
    );
    // The specific regression: the fallback span must not appear at all.
    assert!(
        !err.contains("260-262"),
        "decline must not quote the unit-wide fallback span: {err}"
    );
}

/// The discriminating case, and the reason the fix is not an unconditional
/// reword: when the caller *did* write a range, quoting it is correct and
/// helpful, and that message must not change.
#[test]
fn an_explicit_range_decline_still_quotes_the_callers_own_span() {
    let err = spdi_err(&padded("GGATGCATGG"), "TEMPLATE:g.259_264AT[3]");
    assert!(
        err.contains("repeat span TEMPLATE:259-264 does not match repeat unit AT"),
        "an explicit range is the caller's own span and must still be quoted: {err}"
    );
    assert!(
        !err.contains("is anchored at"),
        "an explicit range has no anchor to report: {err}"
    );
}

/// Nothing that converted before stops converting: an in-phase anchor still
/// resolves to its whole tract.
///
/// Without this, a "fix" that simply refused every single-position anchor would
/// satisfy both tests above.
#[test]
fn an_in_phase_anchor_still_converts() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", padded(TRACT));
    let variant = parse_hgvs("TEMPLATE:g.259CAG[5]").expect("parse");
    assert_eq!(
        hgvs_to_spdi(&variant, &provider)
            .expect("an in-phase anchor must still convert")
            .to_string(),
        "TEMPLATE:258:CAGCAGCAG:CAGCAGCAGCAGCAG"
    );
}

/// The same decline, reached through the **divisibility** check instead of the
/// unit-match one — near the 3' end the unit-wide fallback is clamped inside
/// the contig, so a multi-base unit leaves a span whose length is not a
/// multiple of the unit and that check trips first.
///
/// Before the checks were unified, this path reported
/// `repeat span T2:9-10 length 2 is not a multiple of unit length 3` — quoting
/// `9-10`, which is the clamped fallback and not anything the caller wrote.
/// That is the same defect #1492 removes from the unit-match branch, so it must
/// produce the same anchor-naming message.
#[test]
fn a_failed_search_clamped_at_the_contig_end_still_names_the_anchor() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("T2", "GGGGGGGGTT".to_string());
    let variant = parse_hgvs("T2:g.9CAG[5]").expect("input must parse");
    let err = match hgvs_to_spdi(&variant, &provider) {
        Ok(triple) => panic!("must not convert, got {triple}"),
        Err(e) => e.to_string(),
    };
    assert!(
        err.contains("no CAG repeat is anchored at T2:9"),
        "must name the caller's anchor, got: {err}"
    );
    assert!(
        !err.contains("9-10"),
        "must not quote the clamped fallback span, got: {err}"
    );
}
