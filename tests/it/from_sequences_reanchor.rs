//! Re-anchoring a window: [`SequencePair::trim_to`] and [`Normalizer::reanchor`].
//!
//! # What these are for
//!
//! Holding a derivation inside a region it must not leave — a target region, an
//! amplicon, a tiling window — starting from whatever raw window each caller
//! happens to have. `from_sequences` is a pure function of the window it is
//! given, so anchoring every input to one window makes them agree.
//!
//! # What they are not for, and this matters
//!
//! Making heterogeneous raw pairs agree *in general*. Two existing paths already
//! do that and both reach a **better** answer — the reference-anchored one,
//! which shifts as far as the sequence allows rather than as far as a chosen
//! window allows. [`the_reference_anchored_paths_already_converge`] measures
//! both, side by side with the raw derivations that do not converge, so the
//! comparison is on the record rather than in a doc comment.
//!
//! So a bound is worth reaching for when it is a *requirement*. Anchoring to a
//! window that cuts an ambiguous run makes every caller using that window agree
//! with each other and disagree with the reference — legitimate as a stated
//! contract, misleading as a default, and reported either way through
//! `placement_bounded_by_window`.
//!
//! # The split
//!
//! Trimming removes bases and needs no reference, so it is on `SequencePair`.
//! Padding needs bases the caller does not hold, so it is on `Normalizer`. That
//! is the same reference-free / reference-holding boundary that organises the
//! rest of this surface.

use crate::common::cis_apply_oracle::normalizer_for;
use ferro_hgvs::{
    from_sequences, FromSequencesOptions, MockProvider, Normalizer, SequencePair, ShuffleDirection,
};

/// 1-based:      1234567890
///
/// A 4-base `A` run at 12-15, so a deletion inside it can roll and a bound can
/// visibly stop it rolling.
const SEQUENCE: &str = "GGATTACAGGCAAAAGCCTGAGGATTACAGGCATTAGCCT";

/// The oracle's own provider, which serves `SEQUENCE` unpadded under the
/// accession `TEMPLATE` — reused rather than rebuilt, so this module cannot
/// drift from the fixture every other cis test is written against.
fn normalizer() -> Normalizer<MockProvider> {
    normalizer_for(SEQUENCE, ShuffleDirection::ThreePrime)
}

/// The window a read covering `[lo, hi]` reports when the sample is missing one
/// `A` from the run at 12-15.
fn read(lo: u64, hi: u64) -> SequencePair {
    let reference = &SEQUENCE[(lo - 1) as usize..hi as usize];
    let at = reference[(12 - lo) as usize..]
        .find('A')
        .expect("the read covers part of the run")
        + (12 - lo) as usize;
    let mut alternate = reference.to_string();
    alternate.remove(at);
    SequencePair::new("TEMPLATE", lo, reference, alternate).expect("a well-formed pair")
}

fn derive(pair: &SequencePair) -> String {
    from_sequences(
        &pair.accession,
        pair.position,
        &pair.reference,
        &pair.alternate,
        &FromSequencesOptions::default(),
    )
    .expect("derives")
    .to_string()
}

/// **A caller with bases and no description can build a pair.**
///
/// `SequencePair` is `#[non_exhaustive]` and was otherwise only ever *returned*,
/// which put both re-anchoring entry points out of reach of exactly the caller
/// they are for. Pinned as its own test because the gap is invisible from inside
/// the crate — a struct expression compiles here and nowhere else.
#[test]
fn a_pair_can_be_built_from_bases_alone() {
    let pair = SequencePair::new("TEMPLATE", 10, "GCAAAAG", "GCAAAG").expect("well-formed");
    assert_eq!((pair.position, pair.end()), (10, 16));
    assert!(
        !pair.window_is_final,
        "a caller-supplied window carries no evidence that its 3' edge is the sequence's"
    );
}

/// **`new` refuses exactly what `from_sequences` refuses**, through the same
/// check.
///
/// A constructor that admitted a pair the derivation would later reject would
/// move the error one call away from the argument that caused it.
#[test]
fn a_pair_that_constructs_is_a_pair_that_derives() {
    for (position, reference, alternate) in [(0, "AGCG", "AG"), (1, "", "AG"), (1, "AGXG", "AG")] {
        let built = SequencePair::new("TEMPLATE", position, reference, alternate);
        let derived = from_sequences(
            "TEMPLATE",
            position,
            reference,
            alternate,
            &FromSequencesOptions::default(),
        );
        assert_eq!(
            built.is_err(),
            derived.is_err(),
            "({position}, {reference:?}, {alternate:?}): the constructor and the derivation \
             disagree about whether this is usable"
        );
    }
}

/// **A bound holds the placement where it is put.**
///
/// The headline behaviour. Unbounded, the deletion rolls 3' to 15; bounded at
/// 14, it stays at 14 — and the derivation says so.
#[test]
fn a_trailing_bound_stops_the_roll_at_the_bound() {
    let pair = read(10, 16);
    assert_eq!(derive(&pair), "TEMPLATE:g.15del");

    let bounded = pair
        .trim_to(None, Some(14))
        .expect("14 is inside the window");
    assert_eq!(bounded.position, 10);
    assert_eq!(bounded.end(), 14);
    assert_eq!(derive(&bounded), "TEMPLATE:g.14del");

    // Not a wrong answer — the same variant, spelled at the bound.
    let n = normalizer();
    assert_eq!(
        n.canonical_spdi(&ferro_hgvs::parse_hgvs(&derive(&bounded)).unwrap())
            .unwrap(),
        n.canonical_spdi(&ferro_hgvs::parse_hgvs(&derive(&pair)).unwrap())
            .unwrap(),
        "a bounded placement must still denote the same bases"
    );
}

/// **Every raw window anchored to one region derives one description.**
///
/// The property the feature exists for, over reads that disagree without it.
#[test]
fn reads_anchored_to_one_region_agree() {
    let reads = [read(9, 16), read(10, 17), read(11, 18)];
    let raw: Vec<String> = reads.iter().map(derive).collect();

    let anchored: Vec<String> = reads
        .iter()
        .map(|pair| derive(&pair.trim_to(Some(11), Some(14)).expect("inside every read")))
        .collect();

    assert!(
        anchored.windows(2).all(|w| w[0] == w[1]),
        "anchored reads disagreed: {anchored:?} (raw was {raw:?})"
    );
    assert_eq!(anchored[0], "TEMPLATE:g.14del");
}

/// **Trimming refuses to cut through a base the two sequences disagree on.**
///
/// Trimming a mismatched base would change what the pair denotes, so it is a
/// refusal naming the coordinate rather than a silent narrowing.
///
/// A **substitution** is used rather than the homopolymer deletion the rest of
/// this module runs on, and the reason is worth recording because it caught a
/// wrong assumption here first. Trimming to 13 inside that deletion looks like
/// it should refuse and does not: the last three bases of `GCAAAAG` and of
/// `GCAAAG` are both `AAG`, so they are *matching* bases, and the narrowed pair
/// still denotes one deleted `A`. Legal, and correct. Only where the two
/// sequences genuinely disagree at the edge is there anything to refuse.
#[test]
fn trimming_refuses_to_cut_into_the_changed_block() {
    // 10-16 is `GCAAAAG`; the alternate substitutes the final base.
    let pair = SequencePair::new("TEMPLATE", 10, "GCAAAAG", "GCAAAAT").expect("well-formed");

    let error = pair
        .trim_to(None, Some(15))
        .expect_err("16 is the substituted base, so trimming it away changes the denotation")
        .to_string();
    assert!(
        error.contains("differ") && error.contains("16"),
        "the refusal must say why and where: {error}"
    );

    // The 5' direction is checked separately: the two ends are separate
    // comparisons in the implementation, and a bug in one would hide behind the
    // other if only one were tested.
    let leading = SequencePair::new("TEMPLATE", 10, "GCAAAAG", "TCAAAAG").expect("well-formed");
    let error = leading
        .trim_to(Some(11), None)
        .expect_err("10 is the substituted base")
        .to_string();
    assert!(
        error.contains("differ") && error.contains("10"),
        "the 5' refusal must say why and where: {error}"
    );
}

/// **Trimming to a bound over matching bases is legal even inside a run.**
///
/// The other half of the rule above, pinned so nobody "fixes" the refusal into
/// something stricter than it should be. The trimmed bases match, so the
/// narrowed pair denotes what the wide one denoted — and the derivation moves to
/// the bound, which is the whole feature.
#[test]
fn trimming_inside_a_homopolymer_is_legal_because_the_edges_match() {
    let pair = read(10, 16);
    let narrowed = pair
        .trim_to(None, Some(13))
        .expect("the trimmed bases are the same in both sequences");
    assert_eq!((narrowed.position, narrowed.end()), (10, 13));
    assert_eq!(derive(&narrowed), "TEMPLATE:g.13del");

    let n = normalizer();
    assert_eq!(
        n.canonical_spdi(&ferro_hgvs::parse_hgvs(&derive(&narrowed)).unwrap())
            .unwrap(),
        n.canonical_spdi(&ferro_hgvs::parse_hgvs(&derive(&pair)).unwrap())
            .unwrap(),
        "narrowing over matching bases must not change what is denoted"
    );
}

/// **Trimming refuses to widen, and says which method can.**
///
/// A dead-end refusal is worse than none; this one names its sibling.
#[test]
fn trimming_refuses_to_widen_and_names_reanchor() {
    let pair = read(10, 16);
    for (start, end) in [(Some(5), None), (None, Some(30))] {
        let error = pair
            .trim_to(start, end)
            .expect_err("outside the window")
            .to_string();
        assert!(
            error.contains("reanchor"),
            "the refusal does not name the method that can do this: {error}"
        );
    }
}

/// **`reanchor` pads from the reference and trims, in one call.**
///
/// Both rows are checked on the **bases denoted** as well as on the rendered
/// string, through `canonical_spdi` — which reaches the bases via
/// `apply_to_reference_padded` + `trim_common_flanks` and so shares nothing with
/// the alignment DAG the derivation runs on. Without it the widen-then-narrow
/// row asserts a string and nothing else, and that row is the one that moves the
/// answer (`g.14del` where the raw pair gives `g.15del`) — exactly where a
/// re-anchoring bug would show up as a *different variant* rather than as a
/// different spelling.
#[test]
fn reanchoring_pads_and_trims_to_reach_the_requested_window() {
    let n = normalizer();
    let pair = read(10, 16);
    let denoted = |p: &SequencePair| {
        n.canonical_spdi(&ferro_hgvs::parse_hgvs(&derive(p)).unwrap())
            .unwrap()
    };

    let widened = n
        .reanchor(&pair, Some(5), Some(25))
        .expect("inside the contig");
    assert_eq!((widened.position, widened.end()), (5, 25));
    assert_eq!(widened.reference, SEQUENCE[4..25]);
    assert_eq!(widened.reference.len(), widened.alternate.len() + 1);
    assert_eq!(derive(&widened), "TEMPLATE:g.15del");
    assert_eq!(
        denoted(&widened),
        denoted(&pair),
        "widening must not change what the pair denotes"
    );

    // One call that widens 5' and narrows 3' at once.
    let both = n
        .reanchor(&pair, Some(5), Some(14))
        .expect("inside the contig");
    assert_eq!((both.position, both.end()), (5, 14));
    assert_eq!(derive(&both), "TEMPLATE:g.14del");
    assert_eq!(
        denoted(&both),
        denoted(&pair),
        "the bound moved the spelling to 14; it must not have moved the variant"
    );
}

/// **A window disjoint from the pair's own is refused, and the message says so.**
///
/// `reanchor` moves a window's edges over the reference; it cannot relocate the
/// window, because the changed bases exist only in the pair. That constraint was
/// undocumented and enforced only as a side effect — the narrowing step ran over
/// the intersection, so a disjoint request surfaced as
/// `trim_to(1000, 16) on TEMPLATE: start is past end`: a sibling method's name
/// and two coordinates the caller never passed. The README taught the disjoint
/// call as *the* way to use the method, so this is the refusal a caller is most
/// likely to meet.
#[test]
fn reanchoring_refuses_a_window_disjoint_from_the_pair() {
    let n = normalizer();
    let pair = read(10, 16);

    for (start, end) in [(Some(20), Some(30)), (Some(1), Some(5))] {
        let error = n
            .reanchor(&pair, start, end)
            .expect_err("disjoint from [10, 16]")
            .to_string();
        assert!(
            error.contains("reanchor") && error.contains("does not overlap"),
            "the refusal must name this method and the reason: {error}"
        );
        assert!(
            !error.contains("trim_to"),
            "the refusal must not name a method the caller did not call: {error}"
        );
    }

    // A *partial* overlap is legal — the boundary is overlap, not containment,
    // and a test that only checked the disjoint side would admit an
    // implementation that demanded the pair contain the requested window.
    let partial = n
        .reanchor(&pair, Some(12), Some(30))
        .expect("[12, 30] overlaps [10, 16] and keeps the changed block");
    assert_eq!((partial.position, partial.end()), (12, 30));
    assert_eq!(derive(&partial), "TEMPLATE:g.15del");

    // Overlap is necessary and not sufficient: the overlap must still hold the
    // bases the two sequences disagree on. [16, 30] overlaps by one base and is
    // refused — by `trim_to`, whose rule that is, and whose message says the
    // reference or alternate would be left empty. Pinned so the refusal above
    // is not "fixed" into swallowing this one.
    let error = n
        .reanchor(&pair, Some(16), Some(30))
        .expect_err("the overlap keeps none of the changed block")
        .to_string();
    assert!(
        error.contains("would leave nothing"),
        "an overlap that keeps no changed bases is `trim_to`'s refusal, not the disjoint one: \
         {error}"
    );
}

/// **Case is not a disagreement, in either half of re-anchoring.**
///
/// A soft-masked reference against an upper-case alternate is an ordinary
/// pileup, and both halves got it wrong in opposite directions: `trim_to`
/// byte-compared, so it refused a legal trim claiming the two sequences "first
/// differ" at a coordinate where they do not; `reanchor` upper-cased only the
/// flank it fetched, so widening a masked window returned a mixed-case pair —
/// the very thing `to_sequences` folds its window to prevent.
#[test]
fn masked_bases_trim_and_reanchor_like_any_others() {
    // 10-16 is `gcaaaag` soft-masked; the alternate is the same window
    // upper-cased with one `A` of the run removed.
    let masked = SequencePair::new("TEMPLATE", 10, "gcaaaag", "GCAAAG").expect("well-formed");

    let bounded = masked
        .trim_to(None, Some(14))
        .expect("the trimmed bases match, case aside");
    assert_eq!((bounded.position, bounded.end()), (10, 14));
    assert_eq!(
        bounded.reference, "gcaaa",
        "`trim_to` fetches nothing, so it has no reason to rewrite the caller's bases"
    );
    assert_eq!(derive(&bounded), "TEMPLATE:g.14del");

    // A real disagreement is still refused, so the fold is not a blanket accept.
    let substituted = SequencePair::new("TEMPLATE", 10, "gcaaaag", "GCAAAAT").expect("well-formed");
    assert!(
        substituted.trim_to(None, Some(15)).is_err(),
        "16 differs in base, not in case"
    );

    let n = normalizer();
    let widened = n
        .reanchor(&masked, Some(5), Some(25))
        .expect("inside the contig");
    assert_eq!(
        widened.reference,
        SEQUENCE[4..25],
        "the whole window is folded, so provider bases and caller bases cannot be told apart"
    );
    assert!(
        !widened.alternate.chars().any(|c| c.is_ascii_lowercase()),
        "a mixed-case pair no caller wrote: {}",
        widened.alternate
    );
    assert_eq!(derive(&widened), "TEMPLATE:g.15del");
}

/// **`reanchor` refuses to run off the end of the sequence.**
///
/// An operator decision: a caller who asked for bases that do not exist has a
/// bug upstream, and a window silently clamped back to the contig would hide it.
#[test]
fn reanchoring_refuses_to_leave_the_sequence() {
    let n = normalizer();
    let pair = read(10, 16);
    let length = SEQUENCE.len() as u64;

    assert!(
        n.reanchor(&pair, Some(1), Some(length)).is_ok(),
        "the whole contig is a legal window"
    );
    assert!(
        n.reanchor(&pair, Some(1), Some(length + 1)).is_err(),
        "one base past the end must refuse, not clamp"
    );
    assert!(
        n.reanchor(&pair, Some(0), None).is_err(),
        "position is 1-based; 0 names no base"
    );
    assert!(
        n.reanchor(&pair, Some(14), Some(12)).is_err(),
        "start past end must refuse"
    );
}

/// **`window_is_final` survives a 3' edge that never moved, and is recomputed
/// when it did.**
///
/// A caller-chosen 3' bound short of the sequence's end is what stops the roll,
/// so it clears the flag. Leaving the 3' edge alone settles nothing new and
/// unsettles nothing either, so the input's answer stands — the same rule
/// `trim_to` states as `self.window_is_final && tail == 0`.
///
/// **This test could not fail before**, and that is worth recording rather than
/// counting as coverage. It built its pair with `read()`, i.e. through
/// `SequencePair::new`, which sets `window_is_final: false` unconditionally — so
/// its `false` assertion could not distinguish a correct recomputation from a
/// settled window silently downgraded, and no input it could construct was ever
/// `true`. The `true` side needs a pair from `to_sequences`, the one place a
/// settled window comes from, and nothing in either language fed one to
/// `reanchor` and then read the flag.
#[test]
fn reanchoring_reports_whether_the_window_reaches_the_sequence_end() {
    let n = normalizer();
    let pair = read(10, 16);
    let length = SEQUENCE.len() as u64;

    // The recompute, over an input that carries `false`.
    assert!(
        n.reanchor(&pair, None, Some(length))
            .expect("legal")
            .window_is_final,
        "a window reaching the sequence's own end is settled whatever it came from"
    );
    assert!(
        !n.reanchor(&pair, None, Some(length - 1))
            .expect("legal")
            .window_is_final,
        "one base short of the end, the bound is what stops the roll"
    );
    assert!(
        !n.reanchor(&pair, None, None)
            .expect("legal")
            .window_is_final,
        "an unsettled window is not settled by moving nothing"
    );

    // The carry-through, over an input that carries `true`. `to_sequences` is
    // the only producer of one; the two preconditions are asserted so this test
    // reports a fixture that stopped exercising the case rather than passing
    // vacuously.
    let settled = n
        .to_sequences(&ferro_hgvs::parse_hgvs(&derive(&pair)).unwrap(), 4)
        .expect("re-windows");
    assert!(
        settled.window_is_final,
        "the 3' pad was served in full, so `to_sequences` must report a settled window"
    );
    assert!(
        settled.end() < length,
        "the fixture must stop short of the sequence end, or `end == length` decides every row"
    );

    for (label, start, end) in [
        ("an identity call", None, None),
        ("a 5'-only widen", Some(1), None),
        ("a 5'-only narrow", Some(settled.position + 1), None),
    ] {
        assert!(
            n.reanchor(&settled, start, end)
                .expect("legal")
                .window_is_final,
            "{label} left the 3' edge where it was, so it cannot unsettle it"
        );
    }

    assert!(
        !n.reanchor(&settled, None, Some(settled.end() + 1))
            .expect("legal")
            .window_is_final,
        "widening 3' past a settled edge puts the window back against a chosen bound"
    );
    assert!(
        !n.reanchor(&settled, None, Some(settled.end() - 1))
            .expect("legal")
            .window_is_final,
        "narrowing 3' moves the window in off whatever settled it"
    );
    assert!(
        n.reanchor(&settled, None, Some(length))
            .expect("legal")
            .window_is_final,
        "widening 3' all the way to the sequence end settles it again"
    );
}

/// **The reference-anchored paths already converge, and reach further.**
///
/// Recorded as a measurement because it is the argument against reaching for a
/// bound by default. Three reads that disagree raw agree under either existing
/// path — and agree on `g.15del`, the 3'-most placement, which no fixed window
/// narrower than the run can produce.
#[test]
fn the_reference_anchored_paths_already_converge() {
    let n = normalizer();
    let reads = [read(9, 14), read(10, 15), read(11, 18)];

    let raw: Vec<String> = reads.iter().map(derive).collect();
    assert!(
        raw.windows(2).any(|w| w[0] != w[1]),
        "the reads no longer disagree raw, so this test compares nothing: {raw:?}"
    );

    let normalized: Vec<String> = reads
        .iter()
        .map(|p| {
            n.from_sequences(
                &p.accession,
                p.position,
                &p.reference,
                &p.alternate,
                &FromSequencesOptions::default(),
                true,
            )
            .expect("derives")
            .to_string()
        })
        .collect();

    let round_tripped: Vec<String> = reads
        .iter()
        .map(|p| {
            let variant = from_sequences(
                &p.accession,
                p.position,
                &p.reference,
                &p.alternate,
                &FromSequencesOptions::default(),
            )
            .expect("derives");
            derive(&n.to_sequences(&variant, 128).expect("re-windows"))
        })
        .collect();

    for (label, answers) in [
        ("normalize=true", &normalized),
        ("to_sequences", &round_tripped),
    ] {
        assert!(
            answers.windows(2).all(|w| w[0] == w[1]),
            "{label} did not converge: {answers:?}"
        );
        assert_eq!(
            answers[0], "TEMPLATE:g.15del",
            "{label} converged somewhere other than the 3'-most placement"
        );
    }
}
