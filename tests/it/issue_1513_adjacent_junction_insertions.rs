//! #1513 — two insertions at neighbouring junctions normalized to a description
//! of **different bases** than the input.
//!
//! ```text
//! core  TAAAATTATATTTATTATTT      NM_TEST.1, CDS 1..=15
//!
//! c.[1_2insAA;2_3insTT]   3' ->  c.2_3insTTAA
//!                         5' ->  c.1_2insATTA
//!   input  SPDI  NM_TEST.1:4::TTAA
//!   output SPDI  NM_TEST.1:2::TTAA        <- a different sequence
//! ```
//!
//! The members sit at *distinct, disjoint* junctions, so the allele is
//! well-formed: the applier accepts it and strict mode accepts it. This is not
//! the non-confluence of #1235 / #1419-#1421 — those are two spellings of one
//! variant disagreeing. This was one spelling normalizing to a *different
//! variant*.
//!
//! # The clamp was right and the write was silent
//!
//! Worth recording, because the pass that looked guilty was not.
//!
//! `c.1_2insAA` 3'-shifts through the `A` run to `c.4_5dup`, sweeping across the
//! `insTT` junction at 2. `AA` and `TT` do not commute, so
//! `clamp_sibling_crossing_junctions` bars the crossing and asks for junction 1
//! — instrumented, it computes `limit = Some(1)`, `destination = Some(1)`,
//! exactly right.
//!
//! Then nothing happened. Moving a **two-base `dup`** to junction 1 means the
//! span `c.0_1`, and the CDS axis has no zero, so `translate_member` reverted;
//! `translate_junction_member` read the unchanged variant as *"the translation
//! was refused; nothing to repair"* and returned, leaving the member at `4_5`.
//! The barrier was computed, and then dropped on the floor. The allele went on
//! to merge into a single insertion at the wrong junction.
//!
//! The fix is that a refusal is not "nothing to repair": the member falls
//! through to `respell_at_gap`, which writes it as a plain insertion carrying
//! the rotated payload. A zero-width junction has a spelling everywhere a span
//! may not.
//!
//! Note the sequence-first re-derivation is **not** implicated, though it looks
//! like the obvious suspect. Disabling `canonicalize_from_sequence` only changes
//! the spelling of the wrong answer: `c.[2_3insTT;4_5dup]`, which keys to the
//! same wrong `2::TTAA`. The error is upstream of it.
//!
//! # No oracle saw it
//!
//! `FERRO_ASSERT_REPARSE` passes (the output is well-formed),
//! `FERRO_ASSERT_IN_BOUNDS` passes (every coordinate exists), and
//! `FERRO_ASSERT_IDEMPOTENT` passes (the output is a fixed point). The applier is
//! not wired in as an oracle, so "denotes a different sequence" was invisible.
//! It was found by #1514's `--verify-spdi`, which asks the corpus that question.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// An `A` run at 2-5 for the insertion to shift through, then `TT`.
const CORE: &str = "TAAAATTATATTTATTATTT";

fn coding_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let length = CORE.len() as u64;
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        CORE.to_string(),
        Some(1),
        Some(15),
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn genomic_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_TEST.1", CORE.to_string());
    provider
}

fn normalizer(provider: MockProvider, direction: ShuffleDirection) -> Normalizer<MockProvider> {
    Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(direction),
    )
}

/// The property: the normalized output must denote the same bases as its input.
///
/// Asserted through `canonical_spdi` rather than against a pinned string.
/// That key is derived from the bases a description *produces* and is maximally
/// 3'-shifted, so two spellings of one edit key identically by construction —
/// which makes this survive a later change to which spelling is canonical, while
/// still failing the moment the meaning moves.
fn assert_denotes_the_same_bases(provider: MockProvider, input: &str, direction: ShuffleDirection) {
    let variant = parse_hgvs(input).expect("input must parse");
    let normalizer = normalizer(provider, direction);
    let output = normalizer
        .normalize(&variant)
        .expect("normalization must succeed");
    let key_in = normalizer
        .canonical_spdi(&variant)
        .expect("the input's members are disjoint, so it applies");
    let key_out = normalizer
        .canonical_spdi(&output)
        .expect("the output must apply too");
    assert_eq!(
        key_in.to_string(),
        key_out.to_string(),
        "{input} [{direction:?}] normalized to {output}, which denotes different bases"
    );
}

/// The reported allele, both directions and both input orders.
#[test]
fn adjacent_junction_insertions_keep_their_sequence() {
    for input in [
        "NM_TEST.1:c.[1_2insAA;2_3insTT]",
        "NM_TEST.1:c.[2_3insTT;1_2insAA]",
    ] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert_denotes_the_same_bases(coding_provider(), input, direction);
        }
    }
}

/// The same shape on the genomic axis.
///
/// Not redundant: the corpus reports **zero** `g.` rows for this defect, which
/// reads as "the genomic axis is fine" and is wrong — the corpus simply does not
/// build this geometry there. Measured directly, `g.[1_2insAA;2_3insTT]` gave
/// `g.2_3insTTAA` before the fix, exactly as the `c.` axis did. A corpus zero is
/// a claim about the corpus.
#[test]
fn the_genomic_axis_is_affected_too() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_denotes_the_same_bases(
            genomic_provider(),
            "NC_TEST.1:g.[1_2insAA;2_3insTT]",
            direction,
        );
    }
}

/// Two insertions far enough apart that neither sweeps the other, and two that
/// share a junction.
///
/// The controls for a change that makes a clamp *fire* where it used to no-op:
/// the risk is pulling members back that had no business moving. Both values
/// were measured on the pre-fix revision and are unchanged by it.
#[test]
fn insertions_that_do_not_sweep_a_sibling_are_untouched() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        // Junctions two apart: the shift never crosses the sibling.
        assert_denotes_the_same_bases(
            coding_provider(),
            "NM_TEST.1:c.[1_2insAA;3_4insTT]",
            direction,
        );
        // A lone insertion still shifts all the way, with no sibling to bound it.
        let normalizer = normalizer(coding_provider(), direction);
        let lone = parse_hgvs("NM_TEST.1:c.1_2insAA").unwrap();
        let expected = match direction {
            ShuffleDirection::ThreePrime => "NM_TEST.1:c.4_5dup",
            _ => "NM_TEST.1:c.2_3dup",
        };
        assert_eq!(
            normalizer.normalize(&lone).unwrap().to_string(),
            expected,
            "an unbounded insertion must keep shifting"
        );
    }
}

/// Two insertions on one junction are left as authored, and still are.
///
/// That pair genuinely has no defined order, so it is an overlap conflict rather
/// than something to merge — `detect_insertion_overlaps` reports it and the
/// allele is preserved. This change must not turn it into a merge.
#[test]
fn same_junction_insertions_are_still_preserved() {
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let normalizer = normalizer(coding_provider(), direction);
        let variant = parse_hgvs("NM_TEST.1:c.[2_3insAA;2_3insTT]").unwrap();
        assert_eq!(
            normalizer.normalize(&variant).unwrap().to_string(),
            "NM_TEST.1:c.[2_3insAA;2_3insTT]",
            "a same-junction pair has no defined order and is preserved verbatim"
        );
    }
}
