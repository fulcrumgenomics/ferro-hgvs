//! Issue #1796 — `c.0` / `r.0` name no nucleotide, so the normalizer's
//! conversion sites must refuse them rather than invent a coordinate.
//!
//! `background/numbering.md:31` states, *inside* the definition of the `c.`
//! numbering axis and alongside `-1` and `*1`: *"there is no nucleotide
//! `c.0`."* That is an existence claim, not an input-hygiene rule — the
//! coordinate names no nucleotide, so there is nothing to convert it to.
//! biocommons/hgvs enforces the same in its position type
//! (`BaseOffsetPosition base may not be 0`), and mutalyzer-crossmapper opens by
//! naming this exact hazard: `-1` and `1` being adjacent *"gives rise to
//! various off by one errors when the conversion is not done properly."*
//!
//! Three sites in `src/normalize/mod.rs` answered it anyway, all with
//! `cds_start - 1` — which is the axis's own `-1`, a *different* position from
//! the one asked about. They now refuse. The mapper half of the same seam is
//! #1746/#1747.
//!
//! # Why the interesting assertions here are about what did NOT move
//!
//! These arms were **load-bearing**, and that is the trap this file exists to
//! keep pinned. `collapse_overlapping_cis_edits` computed a 5'-edge insertion
//! anchor as `c.1 - 1 = c.0` and built the intermediate member `c.?_1insG`;
//! `cds_start - 1` happens to be the correct transcript coordinate for `c.-1`,
//! so the arm silently repaired that producer's off-by-one. Refuse without
//! fixing the producer and `NM_TEST.1:c.[1A>G;1dup]` emits the malformed
//! `c.?_1insG` to the user — with **zero** tests going red, because the one
//! test that reached the arm asserted only that the output re-parses, and
//! `c.?_1insG` does.
//!
//! The producer is #1772, fixed by #1777, which this change is stacked on. So
//! `the_five_prime_edge_cis_collapse_is_unaffected` is not a change detector
//! for this diff — it is the guard that the *stack* is right, and it is the
//! single most important assertion in this file. If it ever reports `c.?_1insG`
//! the refusal is sitting on top of an unfixed producer.
//!
//! What this change does move is a position reached by construction rather than
//! by parsing — the shape the Python bindings and any other direct API consumer
//! can build. Reference-free (`MockProvider` via `SyntheticBuilder`), so these
//! hold with no manifest.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::hgvs::edit::NaEdit;
use ferro_hgvs::hgvs::interval::{CdsInterval, RnaInterval};
use ferro_hgvs::hgvs::location::{CdsPos, RnaPos};
use ferro_hgvs::hgvs::variant::{Accession, CdsVariant, LocEdit, RnaVariant};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

/// A 35-base transcript whose CDS is `13..=24`, so the 5'UTR runs `-1`..`-12`
/// and the 3'UTR `*1`..`*11` — real UTRs on both sides, which is what makes the
/// `cds_start - 1` step the arms used to take observable at all.
const CORE: &str = "GCAAAGCGCGCGATGAAACCCTAAGGCATTTTTAA";
const CDS: (u64, u64) = (13, 24);

fn provider() -> MockProvider {
    SyntheticBuilder::cds(CORE, CDS.0, CDS.1, Strand::Plus).build()
}

fn normalize(variant: &HgvsVariant) -> String {
    Normalizer::new(provider())
        .normalize(variant)
        .expect("lenient normalization must not reject")
        .to_string()
}

/// A `c.` variant built directly on base 0 — the shape a caller reaches through
/// the API rather than through the parser, which refuses `c.0` outright (#269).
fn built_cds_base_zero() -> HgvsVariant {
    HgvsVariant::Cds(CdsVariant {
        accession: Accession::new("NM", "TEST", Some(1)),
        gene_symbol: None,
        loc_edit: LocEdit::new(
            CdsInterval::new(CdsPos::new(0), CdsPos::new(0)),
            NaEdit::Deletion {
                sequence: None,
                length: None,
            },
        ),
    })
}

fn built_rna_base_zero() -> HgvsVariant {
    HgvsVariant::Rna(RnaVariant {
        accession: Accession::new("NM", "TEST", Some(1)),
        gene_symbol: None,
        loc_edit: LocEdit::new(
            RnaInterval::new(RnaPos::new(0), RnaPos::new(0)),
            NaEdit::Deletion {
                sequence: None,
                length: None,
            },
        ),
    })
}

// =====================================================================
// 1. The end-to-end behaviour that MOVES
// =====================================================================

/// A `c.` position on base 0 must come back as the caller wrote it, not
/// relabelled as `c.-1`.
///
/// Before this change `cds_to_tx_pos` mapped base 0 to `cds_start - 1` and the
/// whole shuffle proceeded as though the caller had written `c.-1`, so the
/// output was `NM_TEST.1:c.-1del` — a confident answer about a *different*
/// position, with no warning and nothing to tell the caller a substitution had
/// been made. Refusing routes the two call sites in `normalize_cds` to their
/// canonicalize-only fallback, which is the right outcome for a position naming
/// no nucleotide: ferro declines to move it rather than moving it somewhere
/// wrong.
///
/// Note what the input renders as. `CDS_BASE_UNKNOWN` **is** the integer 0
/// (`src/hgvs/location.rs:151`), so `CdsPos { base: 0 }` prints as `c.?` and a
/// hypothetical `c.0` is byte-identical to a genuine `c.?`. That collision is
/// why the old behaviour was doubly wrong here: it answered a position declared
/// *unknown* with a specific coordinate.
#[test]
fn a_cds_position_on_base_zero_is_not_relabelled_as_minus_one() {
    let built = built_cds_base_zero();
    assert_eq!(
        built.to_string(),
        "NM_TEST.1:c.?del",
        "precondition: base 0 and the `c.?` sentinel are the same value, and print alike"
    );
    let output = normalize(&built);
    assert_eq!(
        output, "NM_TEST.1:c.?del",
        "a position naming no nucleotide must be left as authored, not answered"
    );
    assert_ne!(
        output, "NM_TEST.1:c.-1del",
        "`cds_start - 1` is this axis's own `c.-1` — a different position from the one asked \
         about; answering with it is the defect"
    );
}

/// The **second** public exit must refuse base 0 the same way the first does.
///
/// `Normalizer` has exactly two public normalizing methods — `normalize` and
/// `normalize_with_diagnostics` — and a change can reach one without the other.
/// Every other test in this file goes through the `normalize` helper above, and
/// the guards in `src/normalize/mod.rs` exercise `cds_to_tx_pos` /
/// `rna_to_tx_pos` / `cds_pos_to_tx_boundary` *directly*, at the converter. So
/// before this test the diagnostics exit had no base-0 coverage at all.
///
/// The behaviour was already uniform across both exits — that was measured, not
/// assumed. What was missing was a test holding it there, which is the shape
/// that reddens CI when a later change reaches only one exit.
///
/// Asserted on the contract rather than on a convenient side effect: the
/// rendered result, the warnings vector and the infos vector, plus equality
/// with the plain exit's answer. (`NormalizeResult` carries `result` /
/// `warnings` / `infos`; it has no `status` or `changed` fields — those belong
/// to the Python and CLI surfaces, not to this type.)
#[test]
fn the_diagnostics_exit_refuses_base_zero_exactly_as_the_plain_exit_does() {
    let built = built_cds_base_zero();
    let normalizer = Normalizer::new(provider());

    let diagnosed = normalizer
        .normalize_with_diagnostics(&built)
        .expect("lenient normalization must not reject");

    assert_eq!(
        diagnosed.result.to_string(),
        "NM_TEST.1:c.?del",
        "the diagnostics exit must leave a position naming no nucleotide as authored"
    );
    assert_ne!(
        diagnosed.result.to_string(),
        "NM_TEST.1:c.-1del",
        "answering with `cds_start - 1` — this axis's own `c.-1` — is the defect, on either exit"
    );

    // The decline is currently silent, and this pins that as the *measured*
    // contract rather than asserting it is desirable. If the refusal is ever
    // surfaced to the caller (a warning, or the `Err` propagated), this is the
    // assertion that must be updated deliberately, which is the point of it.
    assert!(
        diagnosed.warnings.is_empty(),
        "base 0 currently declines without a warning; got: {:?}",
        diagnosed.warnings
    );
    assert!(
        diagnosed.infos.is_empty(),
        "nothing shuffled, so no shift info may be reported; got: {:?}",
        diagnosed.infos
    );

    // The property this test exists for: the two public exits agree.
    assert_eq!(
        diagnosed.result.to_string(),
        normalize(&built),
        "`normalize` and `normalize_with_diagnostics` must not diverge at base 0"
    );
}

/// The `r.` half of the same shape is pinned at the conversion itself, in
/// `normalize::tests::rna_to_tx_pos_refuses_a_zero_base`, and deliberately has
/// no end-to-end case here.
///
/// `RnaPos` has no unknown sentinel, so `RnaPos { base: 0 }` renders as the
/// literal `r.0` — a string `parse_hgvs` rejects. That makes the axis
/// unreachable end-to-end by any legitimate route: the parser cannot produce
/// base 0, and since #1777 neither can the cis collapse. The only way in is to
/// hand `Normalizer::normalize` a directly-built malformed AST, and doing that
/// trips `FERRO_ASSERT_REPARSE` — *"normalization was handed an input it cannot
/// re-parse, and that input is not one of the deliberate non-renderable
/// shapes"* — so such a test aborts the whole run under the armed CI job rather
/// than asserting anything.
///
/// Worth recording, because the naive reading of that abort is backwards. The
/// oracle is complaining about the **input**, not the output; before this
/// change the same input was silently answered with the parseable `r.-1del`,
/// which is exactly the invention being removed. A normalizer that satisfies
/// the oracle by fabricating a valid coordinate for an invalid one is the
/// defect, so the test was dropped rather than the refusal weakened.
///
/// The `c.` axis has no such problem, and that is not luck: `CDS_BASE_UNKNOWN`
/// is `0`, so `CdsPos { base: 0 }` renders as `c.?`, which round-trips.
#[test]
fn the_rna_axis_has_no_end_to_end_case_and_this_records_why() {
    let built = built_rna_base_zero();
    assert_eq!(
        built.to_string(),
        "NM_TEST.1:r.0del",
        "the premise of the note above: base 0 on `r.` renders as a string the parser rejects"
    );
    assert!(
        parse_hgvs(&built.to_string()).is_err(),
        "if `r.0` ever becomes parseable, this axis gains a legitimate end-to-end case and \
         one should be added alongside the `c.` test above"
    );
}

// =====================================================================
// 2. The stack guard — what must NOT move
// =====================================================================

/// **The load-bearing case.** A cis group netting to a pure insertion at the 5'
/// edge of its window must still normalize to the `-1` anchor, on every
/// transcript axis.
///
/// This is the assertion that says the refusal above is safe. The producer used
/// to hand `cds_to_tx_pos` a base-0 anchor and the arm quietly corrected it; now
/// #1777 names the anchor `-1` at the point the integer becomes a position, so
/// nothing reaches the arm and the output is unchanged. Were this stacked on a
/// base without #1777, `c.` and `r.` would report `?_1insG` / `0_1insg` here —
/// a malformed intermediate promoted to the user-visible answer.
///
/// `n.` is included as the control: `rna_to_tx_pos`/`cds_to_tx_pos` never
/// touched that axis, so it had no repairing arm and #1772 was *visible* there
/// (it emitted the literal `n.0_1insA`, which ferro's own parser rejects). It
/// must be unaffected by this change in either direction.
#[test]
fn the_five_prime_edge_cis_collapse_is_unaffected() {
    for (input, expected) in [
        ("NM_TEST.1:c.[1A>G;1dup]", "NM_TEST.1:c.-1dup"),
        ("NM_TEST.1:r.[1a>g;1dup]", "NM_TEST.1:r.-1dup"),
        ("NM_TEST.1:n.[1G>A;1dup]", "NM_TEST.1:n.-1_1insA"),
    ] {
        let variant = parse_hgvs(input).expect("input must parse");
        let output = normalize(&variant);
        assert_eq!(
            output, expected,
            "`{input}` must be unaffected by the base-0 refusal; a `?` or a `0` in this output \
             means the refusal is sitting on top of an unfixed producer (#1772/#1777)"
        );
        assert!(
            parse_hgvs(&output).is_ok(),
            "`{input}` -> `{output}` must re-parse"
        );
    }
}

/// The refusal is narrow: it removes exactly the one coordinate that names no
/// nucleotide, and every position either side of it still converts and still
/// shuffles.
///
/// Without this the previous two tests are satisfiable by a normalizer that has
/// simply stopped converting `c.`/`r.` positions at all.
#[test]
fn the_positions_either_side_of_the_absent_coordinate_still_normalize() {
    for (input, expected) in [
        // `c.-1` and `c.1` are the two the old arm collapsed onto one answer.
        ("NM_TEST.1:c.-1del", "NM_TEST.1:c.-1del"),
        ("NM_TEST.1:c.1del", "NM_TEST.1:c.1del"),
        ("NM_TEST.1:r.-1del", "NM_TEST.1:r.-1del"),
        ("NM_TEST.1:r.1del", "NM_TEST.1:r.1del"),
        // A 5'UTR deletion that actually shuffles, so the conversion is doing
        // real work rather than passing an already-canonical string through.
        // `c.-10` is transcript position 3, inside this core's `AAA` run
        // (transcript 3..=5), so the 3'-rule carries it to transcript 5 —
        // `c.-8`. Two bases of movement, entirely within the 5'UTR, which is
        // the region whose conversion goes through the arm this change touches.
        ("NM_TEST.1:c.-10del", "NM_TEST.1:c.-8del"),
        ("NM_TEST.1:c.[*1del;*2dup]", "NM_TEST.1:c.*2="),
    ] {
        let variant = parse_hgvs(input).expect("input must parse");
        assert_eq!(normalize(&variant), expected, "`{input}` must be unchanged");
    }
}

/// `c.-1` and `c.1` must still name **different** transcript positions.
///
/// The old arm's specific harm was collapsing the gap: it gave `c.0` `c.-1`'s
/// answer, so three `c.` spellings mapped onto two transcript bases. Deleting
/// one base at `c.-1` and one at `c.1` must therefore denote different edits —
/// asserted through the outputs rather than through the private conversion, so
/// it holds at the layer a consumer sees.
#[test]
fn the_gap_the_axis_skips_is_still_a_gap() {
    let minus_one = normalize(&parse_hgvs("NM_TEST.1:c.-1del").expect("must parse"));
    let plus_one = normalize(&parse_hgvs("NM_TEST.1:c.1del").expect("must parse"));
    assert_ne!(
        minus_one, plus_one,
        "`c.-1` and `c.1` are adjacent bases, not one base; the axis skips 0 between them"
    );
}
