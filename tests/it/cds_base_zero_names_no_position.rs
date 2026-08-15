//! A CDS base of `0` names no nucleotide, and `CoordinateMapper::cds_to_tx`
//! must refuse it rather than compute on it — through the conversion entry
//! point, not only as a unit.
//!
//! # The collapse
//!
//! `CdsPos::is_unknown()` is `base == CDS_BASE_UNKNOWN && !utr3 &&
//! special.is_none()`, so a guard written on that predicate recognises exactly
//! one of the **three** shapes that carry the integer `0`:
//!
//! ```text
//! c.?              base 0, utr3=false, special=None    is_unknown
//! c.*0             base 0, utr3=TRUE                   NOT is_unknown
//! c.pter/qter/cen  base 0, special=Some                NOT is_unknown
//! ```
//!
//! The two it misses fell through to the arithmetic with a zero displacement.
//! Measured on `NM_TEST.1` (cds_start 6, cds_end 35) before the widening:
//! `c.*0` answered tx 35, which is `c.30`'s answer, and each terminus marker
//! answered tx 6, which is `c.1`'s.
//!
//! # Why this file exists beside the unit tests in `src/convert/mapper.rs`
//!
//! The unit tests pin the refusal at the function. This one pins that the
//! collapse was **reachable from a description a user can write**, and that the
//! refusal reaches that far too. Measured on this branch's parent, through
//! `hgvs_to_spdi` on `JsonProvider::with_test_data`:
//!
//! ```text
//! NM_001234.1:c.pterdel  -> NM_001234.1 : 4 : A :
//! NM_001234.1:c.qterdel  -> NM_001234.1 : 4 : A :
//! NM_001234.1:c.cendel   -> NM_001234.1 : 4 : A :
//! NM_001234.1:c.1del     -> NM_001234.1 : 4 : A :
//! ```
//!
//! Four descriptions, one bit-identical triple. `pter` and `qter` are opposite
//! ends of a chromosome, so this is not a near miss in either direction.
//!
//! # What is asserted, and why it is not `is_err()`
//!
//! The failure is a **collapse**, so each test pairs the no-base input with the
//! ordinary coordinate that used to receive its answer and requires the two to
//! be told apart. A test asserting only that some error occurs would stay green
//! against a guard that had widened far enough to refuse the ordinary
//! neighbour too — which would be a worse defect than the one being fixed.

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::spdi::hgvs_to_spdi;

/// The provider `NM_001234.1` comes from — the one every other coding-axis
/// conversion test on this path is written against. Its `cds_start` is 5, so
/// `c.1` is transcript base 5 and SPDI position 4.
fn provider() -> JsonProvider {
    JsonProvider::with_test_data()
}

/// The ordinary coding position the markers used to be answered with. Kept as
/// a named helper so each test below states the collapse in the same terms.
fn ordinary_first_base_triple() -> String {
    let parsed = parse_hgvs("NM_001234.1:c.1del").expect("c.1del parses");
    let triple = hgvs_to_spdi(&parsed, &provider()).expect("c.1 is an ordinary coding position");
    format!(
        "{}:{}:{}:{}",
        triple.sequence, triple.position, triple.deletion, triple.insertion
    )
}

/// `c.pter` — the first of the two arms `is_unknown()` does not see.
///
/// It carries `special = Some(Pter)`, so `is_unknown()` is false and it took
/// the `base < 1` 5'UTR arm with a zero displacement, answering `cds_start`.
#[test]
fn a_pter_marker_no_longer_denotes_the_first_coding_base() {
    let ordinary = ordinary_first_base_triple();
    assert_eq!(
        ordinary, "NM_001234.1:4:A:",
        "the ordinary neighbour must keep its own answer; if this line moves, \
         the comparison below is measuring something else",
    );

    let parsed = parse_hgvs("NM_001234.1:c.pterdel").expect("c.pterdel parses");
    let outcome = hgvs_to_spdi(&parsed, &provider());
    assert!(
        outcome.is_err(),
        "c.pter names a genomic landmark, not a position on NM_001234.1's CDS \
         axis, and must not be answered with c.1's triple ({ordinary}): {outcome:?}",
    );
}

/// `c.qter` — the opposite end of the chromosome from `c.pter`, and it received
/// the identical answer. Pinned separately from `pter` so a regression names
/// which marker re-opened.
#[test]
fn a_qter_marker_no_longer_denotes_the_first_coding_base() {
    let ordinary = ordinary_first_base_triple();
    let parsed = parse_hgvs("NM_001234.1:c.qterdel").expect("c.qterdel parses");
    let outcome = hgvs_to_spdi(&parsed, &provider());
    assert!(
        outcome.is_err(),
        "c.qter is the far end of the chromosome from c.pter and must not be \
         answered with c.1's triple ({ordinary}): {outcome:?}",
    );
}

/// `c.cen` — the third marker. `SpecialPosition` has exactly these three
/// variants, so together the three tests cover the arm exhaustively rather than
/// by sampling.
#[test]
fn a_cen_marker_no_longer_denotes_the_first_coding_base() {
    let ordinary = ordinary_first_base_triple();
    let parsed = parse_hgvs("NM_001234.1:c.cendel").expect("c.cendel parses");
    let outcome = hgvs_to_spdi(&parsed, &provider());
    assert!(
        outcome.is_err(),
        "c.cen is unresolvable on a transcript and must not be answered with \
         c.1's triple ({ordinary}): {outcome:?}",
    );
}

/// The three markers must not merely be refused — they must stop agreeing with
/// each other. Asserted directly, because "all three are errors" and "all three
/// are the same coordinate" are both uniform answers and only one of them is
/// correct.
#[test]
fn the_three_markers_no_longer_share_one_answer_with_each_other() {
    let answers: Vec<String> = ["pter", "qter", "cen"]
        .iter()
        .map(|marker| {
            let input = format!("NM_001234.1:c.{marker}del");
            let parsed = parse_hgvs(&input).expect("the marker description parses");
            match hgvs_to_spdi(&parsed, &provider()) {
                Ok(t) => format!(
                    "{}:{}:{}:{}",
                    t.sequence, t.position, t.deletion, t.insertion
                ),
                Err(_) => format!("refused:{marker}"),
            }
        })
        .collect();

    assert_eq!(
        answers,
        vec!["refused:pter", "refused:qter", "refused:cen"],
        "each marker must be refused on its own terms; before the fix all three \
         resolved to the single triple NM_001234.1:4:A:",
    );
}

/// The refusal is about the base, not about the entry point or the reference:
/// ordinary coding, 5'UTR and 3'UTR positions on the same transcript still
/// convert. Without this a guard that had widened to refuse the whole axis
/// would satisfy every assertion above.
#[test]
fn ordinary_positions_on_the_same_transcript_still_convert() {
    for input in [
        "NM_001234.1:c.1del",
        "NM_001234.1:c.4del",
        "NM_001234.1:c.-1del",
        "NM_001234.1:c.*1del",
    ] {
        let parsed = parse_hgvs(input).expect("the description parses");
        hgvs_to_spdi(&parsed, &provider())
            .unwrap_or_else(|e| panic!("{input} names a real nucleotide and must convert: {e}"));
    }
}

/// The `c.*0` arm's reach, stated as a measurement rather than as prose.
///
/// `parse_hgvs` refuses `c.*0del`, so unlike the terminus markers that arm is
/// **not** reachable from a description. It is reachable only by a **Rust**
/// caller constructing a `CdsPos` — the fields are `pub` — and specifically
/// **not** through the Python bindings: #1741 added `reject_zero_base`, which
/// `c_to_g`, `c_to_p`, `c_to_n` and `n_to_c` all call before converting and
/// which keys on the base alone, so `c_to_n(tx, 0, utr3=True)` is refused at
/// the boundary. That commit is on `main` and was not on this branch's original
/// base, so the reach narrowed at the rebase rather than at this guard. The
/// refusal itself is pinned at the function by
/// `cds_to_tx_refuses_a_zero_three_prime_utr_position` in
/// `src/convert/mapper.rs`.
///
/// Pinned here so the scoping claim fails loudly if the parser ever admits the
/// spelling, rather than sitting in a comment that has quietly gone false.
#[test]
fn a_zero_three_prime_utr_position_is_not_reachable_from_a_description() {
    assert!(
        parse_hgvs("NM_001234.1:c.*0del").is_err(),
        "c.*0 is expected to be unreachable from a string; if it now parses, the \
         3'UTR arm needs an end-to-end guard beside the terminus ones above",
    );
    // The neighbouring spelling the zone does number still parses, so the
    // assertion above is about the zero and not about `*` notation.
    parse_hgvs("NM_001234.1:c.*1del").expect("c.*1 is the first nucleotide the `*` zone numbers");
}
