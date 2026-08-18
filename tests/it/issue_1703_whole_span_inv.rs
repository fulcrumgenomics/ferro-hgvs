//! #1703 — a whole-span reverse complement is typed as one `inv`, uniformly.
//!
//! # The adjudication this pins
//!
//! `rulings[whole-span-reverse-complement-types-as-inv]` (`DNA/inversion.md:5`,
//! 2026-08-13) types a span whose whole content is replaced by its exact reverse
//! complement as one `inv` — **uniformly**, however much of the interior
//! coincides with the reference, and whatever the competing partition is made of
//! (substitutions, `delins`, or anything else). It overturns #1230's
//! competitor-type distinction and supersedes the competitor-type reasoning in
//! `rulings[inversion-vs-two-delins-76-83]` (whose outcome is unchanged). The
//! load-bearing ground is that both sides of the old rule argued over
//! `general.md:56`, which `rulings[adjudication-precedence-order]`'s decided E1
//! holds "cannot settle a merge-versus-split question at all"; and the spec's own
//! published `NM_004006.2:c.4145_4160inv` has this shape (10 of 16 columns
//! changed, two unchanged interior runs) and is simply called `inv`.
//!
//! # Why a family, not one row
//!
//! The defect was competitor-gating: `coalesce_whole_block_inversion` and
//! `inversion_gate_admits` refused the `inv` when the competing partition was
//! lone substitutions whose changed columns did not dominate the span — the
//! all-substitutions case #1230/#1541 name. So the property under test is that
//! interior coincidence carries **no** weight: the cases below range from every
//! column changed to only the two end columns changed, and all must land on the
//! single `inv`. `merge.rs`'s `whole_span_reverse_complement` runs the exact
//! whole-span test ahead of and outside the competitor gate, which is what makes
//! this hold.
//!
//! Each genomic assertion goes through [`assert_padded_preserving`], which
//! projects both input and output through `hgvs_to_spdi` and refuses an output
//! that denotes different bases — so a re-typing that dropped or doubled a base
//! fails here even though every spelling agrees.

use crate::common::synthetic::{assert_padded_preserving, normalize_to_string, SyntheticBuilder};
use ferro_hgvs::reference::transcript::Strand;

const G: &str = "NC_TEST.1";

/// Assert every `spelling` of the whole-span reverse complement of `core`
/// (placed at genomic 257) normalizes to `expected`, denoting the same bases,
/// and that `expected` is itself a fixed point.
fn genomic_converges(core: &str, expected: &str, spellings: &[&str]) {
    for spelling in spellings {
        assert_eq!(
            assert_padded_preserving(core, spelling),
            expected,
            "`{spelling}` on core `{core}` should type as the whole-span inversion `{expected}`",
        );
    }
    assert_eq!(
        assert_padded_preserving(core, expected),
        expected,
        "`{expected}` must be a fixed point",
    );
}

#[test]
fn every_whole_span_reverse_complement_types_as_one_inv_regardless_of_interior() {
    // Two changed columns of four; interior `CG` self-complementary. #1541's
    // locus (`NM_000500.9:c.710_713`, `TCGT`/`ACGA`) in synthetic genomic form.
    genomic_converges(
        "TCGT",
        &format!("{G}:g.257_260inv"),
        &[
            &format!("{G}:g.257_260inv"),
            &format!("{G}:g.257_260delinsACGA"),
            &format!("{G}:g.[257T>A;260T>A]"),
        ],
    );

    // Two changed columns of five; interior columns 258/260 self-complementary,
    // so the substitution spelling is three separated subs. #1703's own locus.
    genomic_converges(
        "GTTAA",
        &format!("{G}:g.257_261inv"),
        &[
            &format!("{G}:g.257_261inv"),
            &format!("{G}:g.257_261delinsTTAAC"),
            &format!("{G}:g.[257G>T;259T>A;261A>C]"),
        ],
    );

    // Two changed columns of six, four unchanged interior bases between them —
    // the widest coincidence, and the case the old competitor gate refused hard.
    genomic_converges(
        "AAGCTA",
        &format!("{G}:g.257_262inv"),
        &[
            &format!("{G}:g.257_262inv"),
            &format!("{G}:g.257_262delinsTAGCTT"),
            &format!("{G}:g.[257A>T;262A>T]"),
        ],
    );

    // Two changed columns of eight, a six-base palindromic interior — a longer
    // span with the sparsest possible change, to show the rule carries no length
    // or density bound.
    genomic_converges(
        "AACCGGTA",
        &format!("{G}:g.257_264inv"),
        &[
            &format!("{G}:g.257_264inv"),
            &format!("{G}:g.257_264delinsTACCGGTT"),
            &format!("{G}:g.[257A>T;264A>T]"),
        ],
    );

    // Every column changed — the opposite extreme from the sparse cases above,
    // and the shape the old density routes already admitted. It must still land
    // on the `inv`, so the uniform rule agrees with the old one where they meet.
    genomic_converges(
        "AAGG",
        &format!("{G}:g.257_260inv"),
        &[
            &format!("{G}:g.257_260inv"),
            &format!("{G}:g.257_260delinsCCTT"),
        ],
    );
}

/// The #1541 locus on the axis it was filed on: `c.710_713` = `TCGT`,
/// `revcomp` `ACGA`, interior `CG` self-complementary. The synthetic CDS core is
/// the whole CDS, so the span is `c.1_4`. Confirms the ruling reaches the coding
/// axis, not only genomic.
#[test]
fn the_coding_axis_types_a_whole_span_reverse_complement_as_inv() {
    let provider = SyntheticBuilder::cds("TCGT", 1, 4, Strand::Plus).build();
    let expected = "NM_TEST.1:c.1_4inv";
    for spelling in [
        "NM_TEST.1:c.1_4inv",
        "NM_TEST.1:c.1_4delinsACGA",
        "NM_TEST.1:c.[1T>A;4T>A]",
    ] {
        assert_eq!(
            normalize_to_string(provider.clone(), spelling),
            expected,
            "`{spelling}` should type as `{expected}` on the coding axis",
        );
    }
}
