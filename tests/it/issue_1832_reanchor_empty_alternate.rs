//! #1832 — `reanchor` must be able to widen a window whose alternate is empty.
//!
//! `Normalizer::reanchor` narrows to the intersection of the requested window
//! and the pair's own before it widens, and delegates that narrowing to
//! `SequencePair::trim_to`. When the request *contains* the pair's window the
//! intersection is the pair's window, so nothing is trimmed — but `trim_to`'s
//! emptiness guard read `head + tail >= alt_len`, which is `0 >= 0` on an empty
//! alternate, and refused.
//!
//! An empty alternate is not exotic: it is what `Normalizer::to_sequences(v, 0)`
//! returns for any pure deletion, so the shape `reanchor` could not widen was
//! the minimal window for the commonest edit a caller holds.
//!
//! # What must NOT change
//!
//! The guard is still doing two jobs for `head + tail > 0`, and the second is
//! easy to miss: it stops a trim consuming a whole string, **and** it is what
//! keeps the "matching bases only" comparison from indexing `a[..head]` past the
//! end of an empty alternate. `a_narrowing_into_the_deleted_block_still_refuses`
//! is the discriminating case — a fix that dropped the guard outright passes
//! every widening assertion here and fails that one.

use crate::common::cis_apply_oracle::normalizer_for;
use ferro_hgvs::{
    from_sequences, FromSequencesOptions, MockProvider, Normalizer, SequencePair, ShuffleDirection,
};

/// 1-based:      1234567890
///
/// A 4-base `A` run at 12-15. The same contig the rest of the re-anchoring tests
/// use, deliberately: this module is about a pair shape they never build, not
/// about a sequence they never draw.
const SEQUENCE: &str = "GGATTACAGGCAAAAGCCTGAGGATTACAGGCATTAGCCT";

fn normalizer() -> Normalizer<MockProvider> {
    normalizer_for(SEQUENCE, ShuffleDirection::ThreePrime)
}

/// The whole `A` run deleted, with no flanking matched base — exactly the pair
/// `to_sequences(g.12_15del, 0)` produces.
fn whole_window_deleted() -> SequencePair {
    SequencePair::new("TEMPLATE", 12, "AAAA", "").expect("a well-formed pair")
}

/// **`reanchor` widens a pair whose alternate is empty**, on either edge or
/// both, and the identity request is not an error either.
///
/// Each result is pinned as bases rather than as `is_ok()`: the whole point of
/// widening is *which* flank arrives, and an assertion that only checked for
/// success would pass on a window that fetched the wrong bases.
#[test]
fn reanchor_widens_a_window_whose_alternate_is_empty() {
    let n = normalizer();
    let pair = whole_window_deleted();

    for (label, start, end, expect) in [
        ("identity", 12u64, 15u64, ("AAAA", "")),
        ("widen 5'", 9, 15, ("GGCAAAA", "GGC")),
        ("widen 3'", 12, 20, ("AAAAGCCTG", "GCCTG")),
        ("widen both", 9, 20, ("GGCAAAAGCCTG", "GGCGCCTG")),
    ] {
        let widened = n
            .reanchor(&pair, Some(start), Some(end))
            .unwrap_or_else(|e| panic!("{label}: reanchor({start}, {end}) refused: {e}"));
        assert_eq!(
            (
                widened.position,
                widened.end(),
                widened.reference.as_str(),
                widened.alternate.as_str()
            ),
            (start, end, expect.0, expect.1),
            "{label}: reanchor did not produce the requested window"
        );
    }
}

/// **A narrowing that would cut into the deleted block still refuses.**
///
/// The discriminating case for the guard. `head = 1` here, so there is genuinely
/// nothing in the alternate to trim, and refusing is both correct and what stops
/// the "matching bases only" comparison below it from indexing out of range.
#[test]
fn a_narrowing_into_the_deleted_block_still_refuses() {
    let n = normalizer();
    let error = n
        .reanchor(&whole_window_deleted(), Some(13), Some(15))
        .expect_err("there is no alternate base to trim")
        .to_string();
    assert!(
        error.contains("would leave nothing"),
        "the refusal is no longer the emptiness one: {error}"
    );
}

/// **A widened window derives, and derives the same description the unwidened
/// one does.**
///
/// The reason the refusal mattered rather than a restatement of it. `g.12_15del`
/// deletes the whole `A` run, so its placement is settled by the run's own
/// extent and does not move as flank arrives — which makes it the case where
/// widening must be a no-op on the *answer* while being the difference between
/// an answer and a refusal on the *call*.
#[test]
fn every_widened_window_derives_the_same_description() {
    let n = normalizer();
    let pair = whole_window_deleted();
    let options = FromSequencesOptions::default();

    let mut derived = std::collections::BTreeSet::new();
    for (start, end) in [(12u64, 15u64), (9, 15), (12, 20), (9, 20)] {
        let widened = n
            .reanchor(&pair, Some(start), Some(end))
            .expect("legal after #1832");
        derived.insert(
            from_sequences(
                &widened.accession,
                widened.position,
                &widened.reference,
                &widened.alternate,
                &options,
            )
            .unwrap_or_else(|e| panic!("window {start}-{end} did not derive: {e}"))
            .to_string(),
        );
    }

    assert_eq!(
        derived,
        ["TEMPLATE:g.12_15del".to_string()].into_iter().collect(),
        "widening the window moved the answer"
    );
}

/// **Widening carries `window_is_final` through, on a pair that had none.**
///
/// `reanchor` computes it as `end == length || (end == pair.end() &&
/// pair.window_is_final)`, and the source records that recomputing only the
/// first disjunct was a defect once — an identity re-anchor and a 5'-only widen
/// both downgraded a settled window. A caller-built pair starts `false`, so the
/// `end == length` arm is the only one reachable here, and it is the one worth
/// pinning: widening to the contig's own end settles a window that was not.
#[test]
fn widening_to_the_sequence_end_settles_the_window() {
    let n = normalizer();
    let pair = whole_window_deleted();
    assert!(
        !pair.window_is_final,
        "a caller-supplied pair carries no evidence about its 3' edge"
    );

    let length = SEQUENCE.len() as u64;
    assert!(
        n.reanchor(&pair, Some(9), Some(length))
            .expect("legal")
            .window_is_final,
        "widening 3' to the sequence end did not settle the window"
    );
    assert!(
        !n.reanchor(&pair, Some(9), Some(length - 1))
            .expect("legal")
            .window_is_final,
        "a window stopping one base short of the end reads as settled"
    );
}

/// **A refusal from the narrowing step names the caller's own request.**
///
/// #1832's second half. `reanchor` delegates narrowing to `trim_to`, whose
/// message names its own method and the *intersection* coordinates — so
/// `reanchor(13, 20)` reported `trim_to(13, 15)`, quoting an end the caller
/// never passed. That is the failure `reanchor`'s non-overlap branch already
/// carries a bespoke message to avoid; this pins the other exit.
#[test]
fn a_narrowing_refusal_names_the_callers_window() {
    let n = normalizer();
    let error = n
        .reanchor(&whole_window_deleted(), Some(13), Some(20))
        .expect_err("13 cuts into the deleted block")
        .to_string();

    assert!(
        error.contains("reanchor(13, 20)"),
        "the refusal does not name the caller's own request: {error}"
    );
    assert!(
        error.contains("[12, 15]"),
        "the refusal does not name the pair's own window, so the caller cannot \
         see why 13-20 was narrowed at all: {error}"
    );
    // The accurate half of `trim_to`'s message is kept rather than replaced —
    // it is what says *why* the narrowing was impossible.
    assert!(
        error.contains("would leave nothing"),
        "the underlying reason was dropped: {error}"
    );
}

/// **A pair carrying a single matched base was always widenable**, and still is.
///
/// The negative control: it is what made the defect look like a property of the
/// *request* rather than of the pair, since every existing re-anchoring test
/// builds its pairs through a helper that always includes flanking matched
/// bases, so none of them could construct an empty alternate.
#[test]
fn a_pair_with_one_matched_base_is_unchanged() {
    let n = normalizer();
    let pair = SequencePair::new("TEMPLATE", 11, "CAAAA", "C").expect("well-formed");
    let widened = n
        .reanchor(&pair, Some(9), Some(20))
        .expect("legal before and after");
    assert_eq!(
        (widened.reference.as_str(), widened.alternate.as_str()),
        ("GGCAAAAGCCTG", "GGCGCCTG"),
        "the control pair's widening moved"
    );
}
