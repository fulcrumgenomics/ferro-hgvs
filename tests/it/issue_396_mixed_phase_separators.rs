//! Issue #396 item 2 — `parse_cds_allele_shorthand` and
//! `parse_rna_allele_shorthand` silently flatten genuinely mixed
//! `;` / `(;)` separators into a single `AllelePhase::Unknown` flat
//! list. There is no spec-defined HGVS shape that says "subgroup A;B
//! is cis to each other but unknown phase to C", so the parser must
//! reject the mixed form, not silently flatten it.
//!
//! `parse_protein_allele_shorthand` already does this (see
//! `parse_protein_allele_shorthand` line ~2755 in
//! `src/hgvs/parser/variant.rs` for the existing reject path); item 2
//! brings cds + rna in line.
//!
//! # Spec basis
//!
//! `assets/hgvs-nomenclature/docs/recommendations/DNA/alleles.md` and
//! `recommendations/RNA/alleles.md`: phase is either fully cis (`;`),
//! fully trans (`];[`), or fully unknown (`(;)`). The grammar does not
//! provide a way to express partial cis-within-unknown nesting in a
//! single bracket pair. A mix is malformed.
//!
//! # Negative controls
//!
//! Pure cis (`a;b`) and pure unknown-phase (`a(;)b`) must continue to
//! parse — only the mixed shape rejects.

use ferro_hgvs::parse_hgvs;

fn assert_parse_err(input: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "expected mixed-separator input to be rejected but it parsed: {input:?} → {:?}",
        result.ok().map(|v| v.to_string()),
    );
}

fn assert_parse_ok(input: &str) {
    parse_hgvs(input).unwrap_or_else(|e| panic!("parse({input:?}) must succeed; got {e:?}"));
}

// =============================================================================
// CDS axis (c.)
// =============================================================================

#[test]
fn cds_mixed_cis_then_unknown_rejected() {
    // `[A;B(;)C]` — A and B written as cis, then unknown phase before C.
    // No spec form for this nesting; must reject.
    assert_parse_err("NM_000088.3:c.[1A>G;2T>C(;)3G>A]");
}

#[test]
fn cds_mixed_unknown_then_cis_rejected() {
    // `[A(;)B;C]` — unknown phase between A and B, then cis B+C.
    assert_parse_err("NM_000088.3:c.[1A>G(;)2T>C;3G>A]");
}

#[test]
fn cds_mixed_three_groupings_rejected() {
    // `[A;B(;)C;D]`
    assert_parse_err("NM_000088.3:c.[1A>G;2T>C(;)3G>A;4C>T]");
}

#[test]
fn cds_pure_cis_still_accepted() {
    assert_parse_ok("NM_000088.3:c.[1A>G;2T>C]");
}

#[test]
fn cds_pure_unknown_phase_still_accepted() {
    assert_parse_ok("NM_000088.3:c.[1A>G(;)2T>C]");
}

// =============================================================================
// RNA axis (r.)
// =============================================================================

#[test]
fn rna_mixed_cis_then_unknown_rejected() {
    assert_parse_err("NM_004006.2:r.[100a>g;200u>c(;)300g>a]");
}

#[test]
fn rna_mixed_unknown_then_cis_rejected() {
    assert_parse_err("NM_004006.2:r.[100a>g(;)200u>c;300g>a]");
}

#[test]
fn rna_pure_cis_still_accepted() {
    assert_parse_ok("NM_004006.2:r.[100a>g;200u>c]");
}

#[test]
fn rna_pure_unknown_phase_still_accepted() {
    assert_parse_ok("NM_004006.2:r.[100a>g(;)200u>c]");
}

// =============================================================================
// Cross-axis: protein already rejects (regression — keep that working)
// =============================================================================

#[test]
fn protein_mixed_still_rejected() {
    assert_parse_err("NP_000079.2:p.[Ala1Val;Ser2Thr(;)Gly3Ala]");
}

// =============================================================================
// Nested-bracket members: the mixed-phase reject and the cis/unknown-phase
// splitter must scan separators at the OUTER bracket depth only. A member
// like `delins[A;T]` carries an inner `;` that is part of the edit, not
// an allele separator, so:
//   - `[X;<delins[A;T]>]`     → cis 2-member (inner `;` not a separator)
//   - `[X(;)<delins[A;T]>]`   → unknown-phase 2-member (no top-level `;`)
//   - `[X;<delins[A;T]>(;)Y]` → MUST reject as mixed-phase only when the
//      top-level depth-0 scan finds both `;` and `(;)`.
// See `parse_*_compound_allele` / `parse_*_allele_shorthand` and the
// `scan_allele_separators` helper for the depth rules.
// =============================================================================

#[test]
fn cds_cis_with_nested_delins_bracket_member() {
    // Inner `;` lives inside `delins[A;T]` and must NOT be treated as an
    // allele-level cis separator. The outer form is pure cis (one top-level `;`).
    assert_parse_ok("NM_000088.3:c.[100_200delins[A;T];300del]");
}

#[test]
fn cds_unknown_phase_with_nested_delins_bracket_member() {
    // Inner `;` inside `delins[A;T]`; outer separator is `(;)`. Pure
    // unknown-phase, NOT mixed-phase.
    assert_parse_ok("NM_000088.3:c.[100_200delins[A;T](;)300del]");
}

#[test]
fn cds_mixed_with_nested_delins_bracket_member_rejected() {
    // Top-level scan finds both `;` and `(;)`. Genuine mixed-phase even
    // though there's also a nested `;` inside `delins[A;T]`.
    assert_parse_err("NM_000088.3:c.[1A>G;100_200delins[A;T](;)300del]");
}

#[test]
fn rna_cis_with_nested_delins_bracket_member() {
    assert_parse_ok("NM_004006.2:r.[100_200delins[a;u];300del]");
}

#[test]
fn rna_unknown_phase_with_nested_delins_bracket_member() {
    assert_parse_ok("NM_004006.2:r.[100_200delins[a;u](;)300del]");
}

#[test]
fn rna_mixed_with_nested_delins_bracket_member_rejected() {
    assert_parse_err("NM_004006.2:r.[100a>g;100_200delins[a;u](;)300del]");
}

#[test]
fn genome_unknown_phase_with_nested_delins_bracket_member() {
    assert_parse_ok("NC_000001.11:g.[100_200delins[A;T](;)300del]");
}

#[test]
fn tx_unknown_phase_with_nested_delins_bracket_member() {
    assert_parse_ok("NM_002001.2:n.[100_200delins[A;T](;)300del]");
}
