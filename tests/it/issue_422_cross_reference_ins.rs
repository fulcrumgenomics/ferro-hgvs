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
    // This used to assert `out.contains("delins")` as a proxy for "the payload
    // was expanded to literal bases". #1835 re-derived the block onto
    // `g.[10_11del;16_17dup]`, where the payload IS fully expanded and the
    // bracket IS gone but no member happens to be spelled `delins`, so the proxy
    // failed while the property it stood for held. It was replaced there with the
    // whole description, which is the stronger test and is kept.
    //
    // #1878 moves the value back. The 6 nt span `g.10_15` (`CGTACG`) is replaced
    // by the 6 nt `g.20_25` (`TACGTA`) — an equal-length block, so
    // `rulings[equal-length-block-column-correspondence-is-unique]` (decided)
    // says the correspondence is unique and there is no rotation alignment to
    // find. Every one of the six columns differs under it, and `delins.md:16`
    // types six consecutive changed nucleotides as one `delins`.
    assert_eq!(
        out, "NC_000022.10:g.10_15delinsTACGTA",
        "the cross-reference must expand to literal bases and then normalize from \
         the resulting sequence; got {out}",
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
    // The normalizer suffix-trims "ATG" to "AT" (trailing G matches ref), then
    // splits the result at an unchanged interior base, since changes separated by
    // unchanged nucleotides are described individually (`delins.md:17`,
    // `general.md:34`). The span still ends at 14, which is what confirms
    // CDS-relative expansion succeeded — transcript-relative "AAA" would reach
    // 15.
    //
    // #1835: was `NC_000022.10:g.[10_12del;14C>T]`; the partition default flip
    // moves it to `g.[10_11delinsA;13_14del]` — the same span, cut at a different
    // unchanged base, because the partition is now the minimal alignment
    // re-derived from the resulting sequence rather than the one
    // `partition_block` reached. Licensed by `canonical-form-choice-when-both-legal`:
    // both cuts are two members separated by an unchanged base, so `general.md:34`
    // is satisfied either way and no clause selects between them; the derivation
    // does.
    //
    // WHAT THIS ROW ACTUALLY PROVES IS UNCHANGED, and it is the reason to keep it
    // pinned exactly: the assertion is evidence that `r.1_3` resolved
    // CDS-relative. That rides on the SPAN ending at 14 rather than 15, which
    // both spellings do, and on "AAA" being absent, which the assertion above
    // checks separately. The cut point was never the evidence.
    assert_eq!(
        out, "NC_000022.10:g.[10_11delinsA;13_14del]",
        "CDS-relative r.1_3 (ATG) must normalize within the 10_14 range; got {out}",
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

// =============================================================================
// Issue #1184 — a trailing `inv` on a cross-reference range.
//
// `parse_external_ref_part` slurps the whole payload (including a trailing
// `inv`) into an opaque `ExternalRef`, so the suffix used to reach
// `parse_cross_reference` glued to the end position and fail its digits-only
// test — rejecting `[ACC:g.A_Binv]` as an unsupported shape even though the
// identical range without `inv` expanded fine. These pin the end-to-end
// behavior through BOTH entry points (bare `Reference` and bracketed
// `ExternalRef`), which matters because the `detect_deferred_part` pre-scan
// gates on `parse_cross_reference` independently of the expansion itself.
//
// `NC_000011.10` is `"GGGGTTTTAAAACCCC"` repeated, so `g.1_6` is `GGGGTT` —
// deliberately not a palindrome, so a bug that stripped the suffix without
// inverting would fail these rather than pass by coincidence.
// =============================================================================

/// Provider for the `inv` tests. The outer accession's replaced range
/// (`g.10_15` == `CTATAG`) is chosen so that **no** position matches either
/// payload — forward `GGGGTT` or inverted `AACCCC` — and both are the same
/// length as the range. Otherwise the normalizer's canonical split factors the
/// coincidental matches out into an allele (e.g.
/// `g.[10C>G;12_15delinsGGTT]`), which would hide the inserted literal this
/// test is about.
fn inv_payload_provider() -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence(
        "NC_000022.10",
        format!("{}CTATAG{}", "G".repeat(9), "A".repeat(10)),
    );
    p.add_genomic_sequence("NC_000011.10", "GGGGTTTTAAAACCCC".repeat(10));
    p
}

/// The forward control and the inverted form differ, and the inverted form is
/// exactly the reverse complement of the forward one.
#[test]
fn cross_reference_inv_suffix_inserts_the_reverse_complement() {
    let forward = normalize(
        "NC_000022.10:g.10_15delins[NC_000011.10:g.1_6]",
        inv_payload_provider(),
    )
    .expect("forward cross-reference must expand");
    assert_eq!(
        forward, "NC_000022.10:g.10_15delinsGGGGTT",
        "forward payload must be the literal GGGGTT",
    );

    let inverted = normalize(
        "NC_000022.10:g.10_15delins[NC_000011.10:g.1_6inv]",
        inv_payload_provider(),
    )
    .expect("inv-suffixed cross-reference must expand (#1184)");
    assert_eq!(
        inverted, "NC_000022.10:g.10_15delinsAACCCC",
        "inv payload must be the reverse complement of GGGGTT",
    );
    assert_ne!(
        forward, inverted,
        "inv must change the inserted sequence, not be silently ignored",
    );
}

/// The bracketed `Complex` path (`ins[<literal>;<cross-ref>inv]`) resolves
/// through `append_part_bases`' `ExternalRef` arm — a separate call site from
/// the bare-`Reference` path above, so it needs its own pin.
#[test]
fn cross_reference_inv_suffix_expands_inside_a_complex_bracket() {
    let out = normalize(
        "NC_000022.10:g.10_15delins[A;NC_000011.10:g.1_6inv]",
        inv_payload_provider(),
    )
    .expect("inv-suffixed cross-reference in a Complex bracket must expand (#1184)");
    // Pin the whole description, not just the `delins…` tail. The payload is the
    // literal `A` followed by the reverse complement `AACCCC`, and the
    // flattening leaves no bracket behind.
    //
    // # RE-PINNED BY THE PARTITION DEFAULT FLIP (#1835) — AND THE KNOWN
    // NON-CONFLUENCE BELOW IS CLOSED BY THE SAME MOVE
    //
    // This block used to pin the spanning `g.10_15delinsAAACCCC`, on the reading
    // that the lone match at `g.12` was a coincidental alignment and that
    // `delins.md:44-47` prefers the span when the payload merely "aligns". Two
    // decided rulings have since scoped `:47` out of this row, in two independent
    // ways — either alone would settle it:
    //
    // (1) THE AXIS. `delins-payload-coincidence-carve-out-is-coding-dna-scoped`
    // scopes `:47` to the coding DNA axis. This is `NC_000022.10:g.`, so `:47`
    // does not reach it and `general.md:34` governs unqualified: members
    // separated by an unchanged nucleotide are described individually.
    //
    // (2) THE DIRECTION. `delins-merge-vs-individual-gap-two-or-more` is scoped to
    // the NET-DELETION case and "does **not** reach net insertions, where the
    // split form stays canonical". A 6 nt span replaced by a 7 nt payload is a net
    // insertion, so even on a coding axis this row would not merge.
    //
    // So the "coincidental-alignment trap" the old comment warned about is a trap
    // only where `:47` reaches, and it reaches neither this axis nor this
    // direction. The #1235 revision that re-blessed this to two members reached
    // the right answer; what it lacked was the scope argument, which the ledger
    // now supplies.
    assert_eq!(
        out, "NC_000022.10:g.[10_11delinsAA;13_15delinsCCCC]",
        "literal `A` then the reverse complement `AACCCC`, described individually \
         across the unchanged base at g.12: `delins.md:47` is coding-axis-scoped \
         and does not reach a net insertion",
    );

    // THE NON-CONFLUENCE THIS TEST PINNED IS GONE. Both spellings now reach the
    // two-member form, so the two assertions below assert convergence where they
    // used to assert two stable strings for one variant.
    //
    // What the old comment said could not be done from here has been done
    // elsewhere: it recorded that the two-member input tripped the
    // `input_separator_positions` veto and was returned verbatim, that overriding
    // that veto would merge `g.[306dup;308C>A]` into `g.307_308delinsCGA` and
    // break #999, and that "the collapse cannot tell the two apart". The
    // canonical partitioner does not need to tell them apart, because it does not
    // start from either input's members — it re-derives the partition from the
    // resulting sequence, so both spellings present the same pieces and the veto
    // has nothing to arbitrate. #999 is unaffected and stays green; see the
    // ruling `separation-is-a-property-of-the-spelling-not-of-the-variant`, which
    // is the general form of this.
    assert_eq!(
        normalize("NC_000022.10:g.10_15delinsAAACCCC", inv_payload_provider()).expect("normalize"),
        "NC_000022.10:g.[10_11delinsAA;13_15delinsCCCC]",
        "the spanning spelling now converges on the split rather than being a \
         second fixed point",
    );
    assert_eq!(
        normalize(
            "NC_000022.10:g.[10_11delinsAA;13_15delinsCCCC]",
            inv_payload_provider()
        )
        .expect("normalize"),
        "NC_000022.10:g.[10_11delinsAA;13_15delinsCCCC]",
        "and the split spelling is the fixed point both now reach",
    );
}

/// A single position carries no orientation, and the same-reference form cannot
/// express one either (both parse paths build a `PositionRangeInv` only from a
/// two-part range). `Ainv` must therefore stay an out-of-scope shape rather
/// than quietly complement one base.
#[test]
fn cross_reference_inv_on_a_single_position_is_still_unsupported() {
    let result = normalize(
        "NC_000022.10:g.10_15delins[NC_000011.10:g.5inv]",
        inv_payload_provider(),
    );
    assert!(
        matches!(result, Err(FerroError::UnsupportedVariant { .. })),
        "single-position inv must stay unsupported; got {result:?}",
    );
}
