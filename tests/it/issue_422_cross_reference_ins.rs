//! Issue #422 — expand cross-reference `ins[ACC:g.A_B]` payloads to
//! their literal sequence.
//!
//! `InsertedSequence::Reference` (bare cross-reference) and
//! `InsertedPart::ExternalRef` (cross-reference inside a Complex
//! bracket) are valid per the HGVS spec but ferro previously returned
//! `FerroError::UnsupportedVariant` instead of expanding them. This
//! file pins the new behavior:
//!
//!   - Same-accession bracketed range (the spec example
//!     `NC_000022.10:g.42522624_42522669delins[NC_000022.10:g.42536337_42536382]`).
//!   - Cross-chromosome translocation
//!     (`NC_000002.12:g.X_Y delins[NC_000011.10:g.A_B]`).
//!   - Complex brackets mixing literal and ExternalRef parts
//!     (`ins[A;NC_000022.11:g.100_200]`).
//!   - Graceful error when an inner accession isn't in the provider.
//!
//! Spec basis: `assets/hgvs-nomenclature/docs/recommendations/DNA/delins.md`
//! (same-accession bracketed range conversion in CYP2D6; cross-chromosome
//! translocation example pter_X delins[pter_Y]).

use ferro_hgvs::error::FerroError;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// Build a MockProvider with two genomic contigs preloaded. The
/// sequences are arbitrary IUPAC bases — only their lengths matter for
/// the position-range fetch.
fn two_genomic_provider() -> MockProvider {
    let mut p = MockProvider::new();
    // NC_000022.10 has bases "ACGTACGT..." padded to 50000+ positions
    // so the spec example's coordinates resolve. We don't need full
    // length — only enough for the range we'll query.
    p.add_genomic_sequence(
        "NC_000022.10",
        // 100 bp: positions 1..100 cyclic ACGT.
        "ACGT".repeat(25),
    );
    p.add_genomic_sequence("NC_000011.10", "GGGGTTTTAAAACCCC".repeat(10));
    p
}

fn normalize(input: &str, provider: MockProvider) -> Result<String, FerroError> {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input)?;
    normalizer.normalize(&variant).map(|v| format!("{}", v))
}

// =============================================================================
// Same-accession bare cross-reference (the spec CYP2D6 shape, simplified
// to coords that fit in our fixture).
// =============================================================================

/// `NC_000022.10:g.10_15delins[NC_000022.10:g.20_25]`: the inner range
/// at positions 20..25 spans 6 bases (1-based inclusive HGVS
/// semantics, `TACGTA` from cyclic `ACGTACGT...`); the outer delins
/// flattens to a literal `delins<6 chars>`. The exact bases depend on
/// the cyclic sequence at position 20.
#[test]
fn same_accession_bare_cross_reference_expands_to_literal() {
    let provider = two_genomic_provider();
    let out = normalize("NC_000022.10:g.10_15delins[NC_000022.10:g.20_25]", provider)
        .expect("must normalize cleanly with both accessions in provider");
    // Expanded form should contain `delins` followed by a literal IUPAC
    // run (the inner range), NOT a `[NC_` bracket reference. Check
    // exhaustively that the cross-ref bracket is gone.
    assert!(
        !out.contains("[NC_"),
        "cross-reference bracket must be flattened to a literal; got {out}",
    );
    assert!(
        out.contains("delins"),
        "expanded form should still have the `delins` keyword; got {out}",
    );
}

// =============================================================================
// Cross-chromosome translocation
// =============================================================================

/// `NC_000002.12:g.X_Y delins[NC_000011.10:g.A_B]` style. Both
/// chromosomes are in the provider, so the inner reference is
/// resolved against `NC_000011.10`.
#[test]
fn cross_chromosome_cross_reference_expands_to_literal() {
    let mut p = MockProvider::new();
    // Outer ref is a poly-C tract (not poly-A): the inner range at
    // NC_000011.10:g.20_25 is `TAAAAG`, so a poly-A outer ref would make the
    // flattened `CCCCCC`→`TAAAAG` delins share interior A-identities and
    // canonicalize to separate substitutions (spec-correct per #165, but not
    // what this test is checking). Poly-C differs from every inserted base, so
    // the flattened form stays a genuine `delins` and the test stays focused on
    // cross-reference flattening.
    p.add_genomic_sequence("NC_000002.12", "C".repeat(100));
    p.add_genomic_sequence("NC_000011.10", "GGGGTTTTAAAA".repeat(10));
    let out = normalize("NC_000002.12:g.5_10delins[NC_000011.10:g.20_25]", p)
        .expect("must normalize cleanly with both contigs in provider");
    assert!(
        !out.contains("[NC_"),
        "cross-chromosome reference must be flattened; got {out}",
    );
    // Inner range NC_000011.10:g.20_25 = `TAAAAG`; against the poly-C outer ref
    // it stays a literal delins over the full span. Pin the exact flattened
    // form, not just the `delins` keyword.
    assert_eq!(out, "NC_000002.12:g.5_10delinsTAAAAG");
}

// =============================================================================
// Complex bracket mixing literal and ExternalRef
// =============================================================================

/// `ins[A;NC_000022.11:g.100_200]` style — literal `A` followed by a
/// cross-reference range. Both parts flatten into a single
/// `Literal("A" + <bases>)`.
#[test]
fn complex_bracket_literal_then_cross_ref_expands_to_single_literal() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    let out = normalize("NC_000022.10:g.5_6ins[A;NC_000022.10:g.10_15]", p)
        .expect("Complex with literal + ExternalRef must expand");
    // The Complex bracket must be gone — either flattened into a
    // literal ins/delins, or further normalized into a `dup` if the
    // resolved bases happen to duplicate the preceding ref bases.
    // Either outcome means the cross-reference is no longer deferred.
    assert!(
        !out.contains("[A;"),
        "Complex bracket must be flattened to a single literal; got {out}",
    );
    assert!(
        !out.contains("NC_000022.10:g.10_15"),
        "the inner cross-reference range must be expanded away; got {out}",
    );
}

// =============================================================================
// Graceful error: inner accession not in provider
// =============================================================================

#[test]
fn cross_reference_unknown_inner_accession_surfaces_reference_not_found() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    // Inner `NC_999999.99` is not registered. The expansion path
    // calls `provider.get_sequence(inner_acc, ...)` which returns
    // `Err(FerroError::ReferenceNotFound | InvalidCoordinates)`. The
    // error must surface — not silently leave the cross-ref intact.
    let result = normalize("NC_000022.10:g.10_15delins[NC_999999.99:g.20_25]", p);
    assert!(
        result.is_err(),
        "unknown inner accession must surface as an error; got Ok({:?})",
        result.ok(),
    );
}

// =============================================================================
// Regression: existing same-reference bracketed shape (no cross-ref)
// =============================================================================

/// `ins[start_end]` with no foreign accession is the existing
/// `InsertedSequence::PositionRange` shape (issue #333). It must
/// continue to expand as before — the #422 changes are strictly
/// additive on the cross-reference paths.
#[test]
fn same_reference_position_range_still_expands() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    let out = normalize("NC_000022.10:g.5_6ins[10_15]", p)
        .expect("same-reference bracketed range must still expand");
    assert!(
        !out.contains("[10_15]"),
        "same-reference position range must flatten; got {out}",
    );
}

// =============================================================================
// Out-of-scope shapes continue to defer (spec-undefined / decoration)
// =============================================================================

/// `pter`/`qter` markers inside the cross-reference are out of scope
/// per the issue's "Out of scope" section. The resolver returns
/// `Unsupported` with a clear message.
#[test]
fn pter_marker_in_cross_reference_stays_deferred() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    p.add_genomic_sequence("NC_000011.10", "ACGT".repeat(50));
    // `pter` decoration → the inner reference isn't a simple
    // positive-integer range, so `parse_cross_reference` returns
    // `None` and the resolver surfaces `UnsupportedVariant`.
    let result = normalize("NC_000022.10:g.10_15delins[NC_000011.10:g.pter_25]", p);
    assert!(
        result.is_err(),
        "pter-marker cross-reference must continue to defer; got Ok({:?})",
        result.ok(),
    );
}

/// CDS-offset range inside a Complex bracket (`InsertedPart::CdsPositionRange`)
/// remains spec-undefined and continues to defer with the
/// `cross-reference is valid HGVS but not yet supported by ferro` /
/// `CDS-offset range` error message.
#[test]
fn cds_offset_range_part_still_deferred() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    // `244-8_249` is a CDS-offset range. Must continue to error.
    let result = normalize("NC_000022.10:g.5_6ins[N[2800];244-8_249]", p);
    assert!(
        result.is_err(),
        "CDS-offset range must continue to defer; got Ok({:?})",
        result.ok(),
    );
}

// =============================================================================
// Additional axis coverage: m. (mito) routed through `Direct`
// =============================================================================

#[test]
fn mito_cross_reference_expands_to_literal() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_012920.1", "ACGT".repeat(50));
    let out = normalize("NC_012920.1:m.5_6ins[NC_012920.1:m.10_15]", p)
        .expect("mito cross-reference must expand");
    assert!(
        !out.contains("[NC_"),
        "mito cross-reference must be flattened; got {out}",
    );
}

// =============================================================================
// r. (RNA) cross-reference now expands (coding-aware, transcript-relative for
// non-coding). #773.
// =============================================================================

#[test]
fn rna_cross_reference_expands_to_literal() {
    let mut p = MockProvider::with_test_data();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    // NR_000123.1 is non-coding ("ACGTACGTACGT"), so r. is transcript-relative:
    // r.1_4 == n.1_4 == "ACGT". The delins payload must flatten to that literal.
    let out = normalize("NC_000022.10:g.10_15delins[NR_000123.1:r.1_4]", p)
        .expect("r.-axis cross-reference must expand");
    assert!(
        !out.contains("[NR_"),
        "r. cross-reference must be flattened; got {out}",
    );
}

#[test]
fn rna_cross_reference_coding_expands_cds_relative() {
    let mut p = MockProvider::with_test_data();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    // NM_001234.1 is coding (cds_start=5 over "AAAAATGCCCAAG…"), so r. is
    // CDS-relative: r.1_3 == c.1_3 == "ATG" (NOT n.1_3 == "AAA"). After
    // expansion the normalizer suffix-trims the trailing "G" that matches the
    // reference base at g.15, collapsing "ATG" → "AT" and the outer range from
    // g.10_15 → g.10_14. If r. had been resolved transcript-relative the
    // inserted bases would have been "AAA", which shares no suffix with the ref
    // run "CGTACG" and would produce g.10_15delinsAAA — distinct from g.10_14.
    // Asserting the g.10_14 outcome therefore proves CDS-relative r. numbering
    // end to end. #773.
    let out = normalize("NC_000022.10:g.10_15delins[NM_001234.1:r.1_3]", p)
        .expect("coding r.-axis cross-reference must expand");
    assert!(
        !out.contains("[NM_"),
        "r. cross-reference must be flattened; got {out}"
    );
    assert!(
        !out.contains("AAA"),
        "transcript-relative AAA must not appear in output; got {out}",
    );
    // The normalizer suffix-trims "ATG" to "AT" (trailing G matches ref).
    // Presence of "g.10_14" in the output confirms CDS-relative expansion
    // succeeded (transcript-relative "AAA" would produce "g.10_15").
    assert!(
        out.contains("g.10_14"),
        "CDS-relative r.1_3 (ATG) must normalize to g.10_14 range; got {out}",
    );
}

// =============================================================================
// Out-of-scope axis: p. (protein) — structurally invalid as DNA payload
// =============================================================================

#[test]
fn protein_cross_reference_continues_to_defer() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    let result = normalize("NC_000022.10:g.10_15delins[NP_000079.2:p.10_15]", p);
    assert!(
        result.is_err(),
        "p.-axis cross-reference is structurally invalid as DNA-insertion \
         payload and must defer; got Ok({:?})",
        result.ok(),
    );
}

// =============================================================================
// Cross-reference with offset (`+N`/`-N`) continues to defer
// =============================================================================

#[test]
fn cross_reference_with_offset_continues_to_defer() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    p.add_genomic_sequence("NC_000011.10", "ACGT".repeat(50));
    // `+5` offset on the inner end position — out of scope per the
    // issue (offsets need exon/intron context the simple
    // position-range fetch doesn't provide).
    let result = normalize("NC_000022.10:g.10_15delins[NC_000011.10:g.20_25+5]", p);
    assert!(
        result.is_err(),
        "offset-bearing cross-reference must defer; got Ok({:?})",
        result.ok(),
    );
}

// =============================================================================
// Order-independence: an unsupported cross-reference must defer even when
// it follows a non-flattenable part (e.g. a Repeat). Regression for the
// order-dependent masking flagged on PR #437: a leading non-flatten part
// used to short-circuit `expand_complex_parts` to `Ok(None)` before the
// later `ExternalRef` arm was validated, silently passing the variant
// through instead of erroring on the unsupported cross-reference.
// =============================================================================

#[test]
fn unsupported_cross_reference_after_repeat_part_still_defers() {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000022.10", "ACGT".repeat(50));
    // `N[2800]` is a non-flattenable Repeat part; the trailing
    // `NP_000079.2:p.10_15` is an out-of-scope p.-axis cross-reference
    // (structurally invalid as a DNA payload). The p.-axis member must still
    // surface an error regardless of its position after the Repeat.
    let result = normalize("NC_000022.10:g.5_6ins[N[2800];NP_000079.2:p.10_15]", p);
    assert!(
        result.is_err(),
        "unsupported cross-reference after a Repeat part must defer; got Ok({:?})",
        result.ok(),
    );
}
