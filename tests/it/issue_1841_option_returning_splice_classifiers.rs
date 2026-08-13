//! Guards for #1841 — the two splice classifiers whose signatures promised an
//! answer they cannot give now return `Option`.
//!
//! #1806 closed #1767 by teaching every reachable classifier about the
//! unknown-offset sentinels (`+?` / `-?`, carried as [`OFFSET_UNKNOWN_POSITIVE`]
//! / [`OFFSET_UNKNOWN_NEGATIVE`]). Two were left with their signatures, and
//! recorded as needing an API decision:
//!
//! * [`IntronicRegion::from_offset`] returned `Self`, so a sentinel landed in
//!   `DeepIntronic` — every variant of that enum names a **distance band**, and
//!   a sentinel supports none of them.
//! * [`EffectPredictor::classify_splice_variant`] returned a bare
//!   `ProteinEffect`, so the best it could do was `IntronVariant`/`Modifier` —
//!   true, but still a positive classification to a caller matching on the
//!   consequence.
//!
//! Neither was a live defect at the point the break was taken: `from_offset`'s
//! two in-repo callers screened first, and `classify_splice_variant`'s answer
//! was true. The claim these guards make is therefore about the **contract**,
//! not about a wrong number — which is why each one below pins the discriminating
//! neighbour alongside the sentinel.
//!
//! **The discriminating case is `i64::MIN + 1` / `i64::MAX - 1`.** They are not
//! sentinels, so they must still classify. A decline keyed on "a very large
//! magnitude" rather than on the sentinel would pass every assertion that only
//! looked at `i64::MIN` and `i64::MAX`, and it is exactly the gap #1765 records
//! one layer up.
//!
//! [`OFFSET_UNKNOWN_POSITIVE`]: ferro_hgvs::hgvs::parser::position::OFFSET_UNKNOWN_POSITIVE
//! [`OFFSET_UNKNOWN_NEGATIVE`]: ferro_hgvs::hgvs::parser::position::OFFSET_UNKNOWN_NEGATIVE

use std::io::Write;
use std::process::Command;

use ferro_hgvs::convert::noncoding::IntronicRegion;
use ferro_hgvs::effect::{Consequence, EffectPredictor, Impact};
use ferro_hgvs::hgvs::location::{CdsPos, TxPos};
use ferro_hgvs::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;

/// Both sentinels, so each guard is exercised symmetrically — a fix that taught
/// the classifiers only about `i64::MIN` would leave `+?` answering by accident.
const SENTINELS: [i64; 2] = [OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE];

/// The neighbours of the sentinels, which are ordinary (if absurd) magnitudes.
/// Pinning these is what stops the decline from widening into "big offsets do
/// not classify".
const SENTINEL_NEIGHBOURS: [i64; 2] = [i64::MIN + 1, i64::MAX - 1];

fn cds_with_offset(offset: i64) -> CdsPos {
    CdsPos {
        base: 100,
        offset: Some(offset),
        utr3: false,
        special: None,
    }
}

fn tx_with_offset(offset: i64) -> TxPos {
    TxPos {
        base: 100,
        offset: Some(offset),
        downstream: false,
    }
}

/// `IntronicRegion::from_offset` declines a sentinel, and only a sentinel.
#[test]
fn intronic_region_from_offset_declines_exactly_the_sentinels() {
    for offset in SENTINELS {
        assert_eq!(
            IntronicRegion::from_offset(offset),
            None,
            "every variant of IntronicRegion names a distance band, and an \
             unknown offset ({offset}) supports none of them"
        );
    }

    for offset in SENTINEL_NEIGHBOURS {
        assert_eq!(
            IntronicRegion::from_offset(offset),
            Some(IntronicRegion::DeepIntronic),
            "{offset} is a measured distance, not a sentinel — the decline must \
             not widen into a magnitude test"
        );
    }
}

/// The two `Option`-returning constructors that already declined must keep
/// declining now that the rule lives in `from_offset` rather than in a guard of
/// their own. This is the discriminating test for the refactor: collapsing the
/// duplicated `has_unknown_offset()` screens is only safe if the behaviour they
/// provided survives, on both axes.
#[test]
fn the_position_constructors_still_decline_after_the_screens_were_collapsed() {
    for offset in SENTINELS {
        assert_eq!(IntronicRegion::from_cds_pos(&cds_with_offset(offset)), None);
        assert_eq!(IntronicRegion::from_tx_pos(&tx_with_offset(offset)), None);
    }

    for offset in SENTINEL_NEIGHBOURS {
        assert_eq!(
            IntronicRegion::from_cds_pos(&cds_with_offset(offset)),
            Some(IntronicRegion::DeepIntronic),
        );
        assert_eq!(
            IntronicRegion::from_tx_pos(&tx_with_offset(offset)),
            Some(IntronicRegion::DeepIntronic),
        );
    }

    // A non-intronic position has no offset at all, which is a different reason
    // to answer `None` and must not have been folded into the sentinel one.
    let exonic = CdsPos {
        base: 100,
        offset: None,
        utr3: false,
        special: None,
    };
    assert_eq!(IntronicRegion::from_cds_pos(&exonic), None);
}

/// `classify_splice_variant` declines a sentinel, and only a sentinel.
#[test]
fn classify_splice_variant_declines_exactly_the_sentinels() {
    let predictor = EffectPredictor::new();

    for offset in SENTINELS {
        assert!(
            predictor.classify_splice_variant(offset).is_none(),
            "an unknown offset ({offset}) states no distance, so there is \
             nothing to classify"
        );
    }

    for offset in SENTINEL_NEIGHBOURS {
        let effect = predictor
            .classify_splice_variant(offset)
            .unwrap_or_else(|| panic!("{offset} is a measured distance and must classify"));
        assert_eq!(
            effect.consequences,
            vec![Consequence::IntronVariant],
            "{offset} is far from any boundary"
        );
        assert_eq!(effect.impact, Impact::Modifier);
        assert_eq!(effect.intronic_offset, Some(offset));
    }
}

/// Every ordinary rung still answers. Without this the change could be read as
/// "the classifier declines more often than it used to" in general, when it
/// declines on exactly one input class.
#[test]
fn every_measured_rung_still_classifies() {
    let predictor = EffectPredictor::new();
    for offset in [-2, -1, 1, 2, 5, -8, 9, 50, -1_000] {
        assert!(
            predictor.classify_splice_variant(offset).is_some(),
            "offset {offset} is measured and must still classify"
        );
        assert!(
            IntronicRegion::from_offset(offset).is_some(),
            "offset {offset} is measured and must still classify"
        );
    }
}

/// The reachability half: the sentinel these functions decline is one a *parsed
/// description* actually carries. Asserted against the AST rather than a
/// rendered round-trip, so nothing here depends on how `-?` is printed.
///
/// Both spellings are exercised. `+?` is the half a fix keyed on `i64::MIN`
/// alone would leave answering by accident, which is #1767's own warning and
/// the reason `SENTINELS` exists.
#[test]
fn a_parsed_unknown_offset_reaches_the_classifiers_as_a_sentinel() {
    let predictor = EffectPredictor::new();

    for (description, expected) in [
        ("NM_000088.3:c.100-?del", OFFSET_UNKNOWN_NEGATIVE),
        ("NM_000088.3:c.100+?del", OFFSET_UNKNOWN_POSITIVE),
    ] {
        let parsed = parse_hgvs(description).expect("an unknown offset is legal HGVS");
        let HgvsVariant::Cds(cv) = &parsed else {
            panic!("expected a c. variant from {description}, got {parsed:?}");
        };
        let offset = cv
            .loc_edit
            .location
            .start
            .inner()
            .expect("a certain start position")
            .offset
            .expect("an unknown offset is still an intronic offset");

        assert_eq!(
            offset, expected,
            "the premise of both declines is that the parser hands {description}'s \
             sentinel on unchanged"
        );

        assert!(predictor.classify_splice_variant(offset).is_none());
        assert_eq!(IntronicRegion::from_offset(offset), None);
    }
}

// ---------------------------------------------------------------------------
// The CLI surface — the one shipped output this change actually moves
// ---------------------------------------------------------------------------
//
// The library half above pins a *contract*; `ferro effect` is the only place a
// user sees the change, and it was evidenced only by a before/after table in the
// PR description. A table is a measurement: it does not survive the merge and it
// cannot fail. These four turn each row of it into a guard, using the spawned-
// binary idiom the repo already runs in twenty `tests/it/cli_*.rs` files.
//
// The CLI path instruction treats stdout as an API, so each assertion names the
// stream it is about: the refusal is on **stderr**, stdout stays clean, and the
// batch loop keeps writing.

/// The `ferro` binary this test crate was built alongside.
fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

/// `ferro effect` refuses an unknown intronic offset rather than predicting one.
///
/// Both spellings, for the same reason `SENTINELS` has two entries: a change
/// keyed on `i64::MIN` alone leaves `+?` answering by accident.
#[test]
fn the_effect_cli_refuses_an_unknown_intronic_offset() {
    for description in ["NM_000088.3:c.100-?del", "NM_000088.3:c.100+?del"] {
        let out = ferro()
            .args(["effect", description])
            .output()
            .expect("the ferro binary runs");

        assert!(
            !out.status.success(),
            "{description} must exit non-zero — before #1841 it printed \
             `intron_variant (MODIFIER)` and exited 0"
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            stdout.trim().is_empty(),
            "a refusal must put nothing on stdout, which a pipeline parses: {stdout}"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("unknown intronic offset"),
            "the refusal must say why, on stderr: {stderr}"
        );
    }
}

/// The neighbour control. Without it every assertion above is satisfied by a
/// `ferro effect` that refuses *everything*, which is the failure mode a
/// decline-shaped change actually has.
#[test]
fn the_effect_cli_still_answers_a_measured_offset() {
    for (description, expected) in [
        ("NM_000088.3:c.100-2del", "splice_acceptor_variant"),
        ("NM_000088.3:c.100del", "coding_sequence_variant"),
    ] {
        let out = ferro()
            .args(["effect", description])
            .output()
            .expect("the ferro binary runs");

        assert!(
            out.status.success(),
            "{description} states a distance and must still be answered; stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            stdout.contains(expected),
            "{description} must still report {expected}: {stdout}"
        );
    }
}

/// Under `-f json` a declined row emits **no JSON object at all**, rather than
/// an object carrying a fabricated consequence.
///
/// Pinned because it is the consumer-visible half of the refusal and the one a
/// reader of the code would not predict: a pipeline reading JSON lines gets one
/// fewer line than it sent, and learns which one only from stderr.
#[test]
fn the_effect_cli_emits_no_json_object_for_a_declined_row() {
    let out = ferro()
        .args(["effect", "NM_000088.3:c.100-?del", "--format", "json"])
        .output()
        .expect("the ferro binary runs");

    assert!(!out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        !stdout.contains('{'),
        "a declined row must not emit a partial or fabricated JSON object: {stdout}"
    );
}

/// Batch mode does not regress into #1764: one declining line is reported on
/// stderr and the stream continues.
///
/// The PR that took this break verified it *by reading both loops*. Reading is
/// not a guard — this is, and it is cheap: the declining line is placed first so
/// a truncating implementation cannot pass by accident.
#[test]
fn the_effect_cli_batch_mode_does_not_truncate_on_a_declining_line() {
    let mut input = tempfile::Builder::new()
        .suffix(".txt")
        .tempfile()
        .expect("a temp input file");
    writeln!(input, "NM_000088.3:c.100-?del").unwrap();
    writeln!(input, "NM_000088.3:c.100-2del").unwrap();
    input.flush().unwrap();

    let out = ferro()
        .args([
            "effect",
            "--input",
            input.path().to_str().expect("a utf-8 temp path"),
        ])
        .output()
        .expect("the ferro binary runs");

    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stdout.contains("splice_acceptor_variant"),
        "the line AFTER the declining one must still be answered — a truncated \
         stream is the #1764 failure mode: stdout {stdout}, stderr {stderr}"
    );
    assert!(
        stderr.contains("ERROR:") && stderr.contains("unknown intronic offset"),
        "the declining line must be reported on stderr rather than dropped: {stderr}"
    );
}
