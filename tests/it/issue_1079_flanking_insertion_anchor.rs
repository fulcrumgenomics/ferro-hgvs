//! Insertion anchors MUST name two *flanking* (adjacent) positions.
//!
//! Two sub-rules, both MUST-level, and both previously unenforced on the
//! paths exercised here.
//!
//! 1. **A single-position anchor is not allowed.** For protein,
//!    `protein/insertion.md:18`:
//!
//!    > an insertion can not be described using **one** amino acid
//!    > position, like `p.Lys23insAsp`.
//!
//!    The DNA half of this rule (`DNA/insertion.md:95-101`,
//!    `checklist.md:30-31` — "The format `c.52insT` is **ambiguous**, and
//!    not allowed") is already enforced by
//!    `reject_single_position_insertion.rs`; the protein axis was
//!    explicitly left as "handled elsewhere" and in fact was not handled
//!    at all. These tests close that gap.
//!
//! 2. **The two positions MUST be adjacent.** `DNA/insertion.md:15`:
//!
//!    > `positions flanking` should contain **two flanking nucleotides**,
//!    > e.g., `123_124`, not `123_125`.
//!
//!    The formal grammar states this as a hard requirement rather than a
//!    preference — `syntax.yaml` `dna.ins` notes "`range` must be
//!    specified using adjacent positions".
//!
//!    Sub-rule 2 is enforced on the **nucleotide axes only**. The protein
//!    twin (`protein/insertion.md:17`) is contradicted by the spec's own
//!    `protein/alleles.md:61` example; see the block comment below.
//!
//! Rejection (not repair) is the right disposition for both: given
//! `g.123_125insG` there is no way to recover which flanking pair the
//! author meant, exactly as `DNA/insertion.md:95-101` argues for the
//! single-position form. A wide anchor is a `delins`, not an `ins`.
//!
//! Enforcement is deliberately limited to anchors whose adjacency is
//! decidable from the description alone — plain positions with no
//! intronic offset, no `pter`/`qter`/`cen` marker, no `?`, and both
//! endpoints in the same coordinate region. Anything needing reference
//! data to settle (e.g. an exon/intron boundary such as
//! `c.100+1_101-1insG`, or a CDS-to-3'UTR span) is left alone; see the
//! guard tests at the bottom.
//!
//! Part of the #1079 spec-divergence burn-down.

use ferro_hgvs::error::ErrorCode;
use ferro_hgvs::parse_hgvs;

/// Assert `input` is refused with a structured `InvalidEdit` code and an
/// message that cites the spec, mirroring the existing single-position
/// DNA rejection so the refusal survives slash-form fallback paths.
fn assert_rejects(input: &str, expected_fragment: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "expected insertion anchor {input:?} to be rejected; got: {:?}",
        result.map(|v| v.to_string())
    );
    let err = result.unwrap_err();
    assert_eq!(
        err.code(),
        Some(ErrorCode::InvalidEdit),
        "insertion-anchor rejection must carry a structured InvalidEdit code for {input:?}",
    );
    let msg = err.to_string();
    assert!(
        msg.contains(expected_fragment),
        "error message must explain the problem; want {expected_fragment:?}, got: {msg}"
    );
    assert!(
        msg.contains("insertion.md"),
        "error message must cite the spec; got: {msg}"
    );
}

/// Assert `input` parses and round-trips to `expected`.
fn assert_accepts(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("{input:?} must parse: {e}"));
    assert_eq!(parsed.to_string(), expected, "round-trip mismatch");
}

// ---------------------------------------------------------------------
// Sub-rule 1 — single-position anchor, protein axis (protein/insertion.md:18)
// ---------------------------------------------------------------------

#[test]
fn rejects_protein_single_position_insertion() {
    // The spec's own counter-example, verbatim.
    assert_rejects("NP_003997.1:p.Lys23insAsp", "single-position insertion");
}

#[test]
fn rejects_protein_single_position_insertion_one_letter() {
    // One-letter input normalizes to three-letter internally; the refusal
    // must not depend on which spelling the author used.
    assert_rejects("NP_003997.1:p.K23insD", "single-position insertion");
}

#[test]
fn rejects_protein_single_position_insertion_predicted() {
    // Parenthesized (predicted) consequence is the same edit and is
    // equally ambiguous.
    assert_rejects("NP_003997.1:p.(Lys23insAsp)", "single-position insertion");
}

#[test]
fn rejects_protein_single_position_insertion_inside_allele() {
    // The rule applies per-variant inside allele brackets.
    assert_rejects(
        "NP_003997.1:p.[Trp24Cys;Lys23insAsp]",
        "single-position insertion",
    );
}

// ---------------------------------------------------------------------
// Sub-rule 2 on the PROTEIN axis — NOT enforced, and deliberately so.
//
// `protein/insertion.md:17` says the anchor should be `Lys23_Leu24`, not
// `Lys23_Asn25` — but the spec contradicts itself. `protein/alleles.md:61`
// publishes, unmarked and as a valid illustration of the allele-bracket
// format:
//
// > **NOTE**: for other variant types, the format is `p.[Ser68del];[Ser68=]`,
// > `p.[Ser68_Arg70dup];[Ser68_Arg70=]`,
// > `p.[Ser68_Ala74insSerGln];[Ser68_Ala74=]`, etc.
//
// whose insertion anchor `Ser68_Ala74` spans seven residues. Enforcing
// adjacency on `p.` would make ferro reject a description the spec itself
// prints, so the rule is left unenforced pending upstream clarification.
// The nucleotide axes carry no such contradicting example (this is the only
// conflict across the whole spec-derived corpus) and ARE enforced above.
//
// These tests pin the non-enforcement so the contradiction is visible and a
// future change to it is deliberate rather than accidental.
// ---------------------------------------------------------------------

#[test]
fn protein_non_flanking_insertion_is_not_rejected_pending_spec_conflict() {
    // `protein/insertion.md:17`'s own counter-example shape.
    assert_accepts(
        "NP_003997.1:p.Lys23_Asn25insAsp",
        "NP_003997.1:p.Lys23_Asn25insAsp",
    );
}

#[test]
fn protein_allele_spec_example_with_wide_insertion_anchor_parses() {
    // `protein/alleles.md:61`, verbatim — the description that forces the
    // non-enforcement above. This MUST keep parsing.
    let s = "NP_003997.1:p.[Ser68_Ala74insSerGln];[Ser68_Ala74=]";
    assert!(
        parse_hgvs(s).is_ok(),
        "spec-published example {s:?} must not be rejected"
    );
}

// ---------------------------------------------------------------------
// Sub-rule 2 — non-adjacent anchor, nucleotide axes (DNA/insertion.md:15)
// ---------------------------------------------------------------------

#[test]
fn rejects_genome_non_flanking_insertion() {
    // The spec's own counter-example shape: `123_125`, not `123_124`.
    assert_rejects("NC_000023.11:g.123_125insG", "not flanking");
}

#[test]
fn rejects_cds_non_flanking_insertion() {
    assert_rejects("NM_004006.2:c.100_102insT", "not flanking");
}

#[test]
fn rejects_cds_utr5_non_flanking_insertion() {
    // Both endpoints in the 5'UTR: `c.-10` and `c.-8` are two apart.
    assert_rejects("NM_004006.2:c.-10_-8insT", "not flanking");
}

#[test]
fn rejects_noncoding_non_flanking_insertion() {
    assert_rejects("NR_004430.2:n.789_791insC", "not flanking");
}

#[test]
fn rejects_mt_non_flanking_insertion() {
    assert_rejects("NC_012920.1:m.3243_3245insA", "not flanking");
}

#[test]
fn rejects_non_flanking_insertion_inside_allele() {
    assert_rejects("NM_004006.2:c.[100A>G;200_202insT]", "not flanking");
}

// ---------------------------------------------------------------------
// Guards — legal forms that MUST keep parsing
// ---------------------------------------------------------------------

#[test]
fn accepts_protein_adjacent_insertion() {
    assert_accepts(
        "NP_003997.1:p.Lys23_Leu24insAsp",
        "NP_003997.1:p.Lys23_Leu24insAsp",
    );
}

#[test]
fn accepts_protein_adjacent_insertion_spec_examples() {
    // Straight from `protein/insertion.md` — every one of these must stay
    // legal, including the length- and stop-bearing insert forms.
    for s in [
        "NP_004371.2:p.(Pro46_Asn47insSerSerTer)",
        "NP_003997.1:p.Arg78_Gly79insXaa[23]",
        "NP_003997.1:p.(Ser332_Ser333insXaa)",
        "NP_003997.1:p.His4_Gln5insAla",
    ] {
        assert_accepts(s, s);
    }
}

#[test]
fn accepts_nucleotide_adjacent_insertions() {
    // Straight from `DNA/insertion.md` examples.
    for s in [
        "NC_000023.10:g.32867861_32867862insT",
        "NM_004006.2:c.169_170insA",
        "NM_004006.2:c.849_850ins858_895",
    ] {
        assert_accepts(s, s);
    }
}

#[test]
fn accepts_utr5_to_cds_boundary_insertion() {
    // `c.-1` and `c.1` ARE adjacent — there is no `c.0`. The naive
    // "end == start + 1" test would wrongly reject this.
    assert_accepts("NM_004006.2:c.-1_1insT", "NM_004006.2:c.-1_1insT");
}

#[test]
fn accepts_intronic_adjacent_insertion() {
    assert_accepts(
        "NM_004006.2:c.100+1_100+2insT",
        "NM_004006.2:c.100+1_100+2insT",
    );
}

#[test]
fn accepts_exon_intron_boundary_insertion() {
    // Adjacency across an exon/intron junction cannot be decided from the
    // description alone, so it must be left alone rather than rejected.
    assert_accepts(
        "NM_004006.2:c.100+1_101-1insT",
        "NM_004006.2:c.100+1_101-1insT",
    );
}

#[test]
fn accepts_cds_to_utr3_boundary_insertion() {
    // Needs the CDS length to decide adjacency — not decidable here.
    assert_accepts("NM_004006.2:c.100_*1insT", "NM_004006.2:c.100_*1insT");
}

#[test]
fn accepts_uncertain_and_unknown_anchors() {
    // `?` and bounded-range boundaries carry no decidable adjacency; the
    // uncertain form is explicitly blessed by `DNA/insertion.md:23`.
    //
    // Acceptance is what matters here, not the exact rendering: ferro
    // independently collapses `g.?_?` to `g.?` on Display, which is out of
    // scope for this rule.
    for s in ["NC_000023.11:g.?_?insG", "NC_000023.11:g.(100_150)insN[25]"] {
        assert!(
            parse_hgvs(s).is_ok(),
            "{s:?} has no decidable adjacency and must not be rejected"
        );
    }
}

#[test]
fn unrelated_edits_are_unaffected() {
    // Non-insertion edits over a multi-position span stay legal — the new
    // rule must key on the edit being an insertion, not on span width.
    for s in [
        "NM_004006.2:c.79G>T",
        "NC_000023.11:g.123_125del",
        "NC_000023.11:g.123_125dup",
        "NC_000023.11:g.123_125inv",
        "NC_000023.11:g.123_125delinsTTT",
        "NP_003997.1:p.Arg97fsTer2",
        "NP_003997.1:p.Lys23_Asn25del",
        "NP_003997.1:p.Lys23_Asn25delinsTrp",
    ] {
        assert_accepts(s, s);
    }
}
