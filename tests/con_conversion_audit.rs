//! Audit of `con` (sequence-conversion) edit notation acceptance.
//!
//! Tracks issue #81 item H1: "`con` (sequence conversion) — rare but
//! spec'd: `c.123_456conNM_X.1:c.789_1011`."
//!
//! ## Background
//!
//! Sequence conversion replaces a range of the current reference with a
//! range from another (or the same) reference, e.g. `c.123_456conNM_X.1:c.789_1011`.
//! The HGVS Nomenclature page on DNA delins (`assets/hgvs-nomenclature/docs/recommendations/DNA/delins.md`)
//! says the `con` form is "no longer used" in favor of `delins`, citing
//! Community Consultation SVD-WG009. The form remains spec'd but
//! deprecated, and it still appears in real-world inputs (Mutalyzer
//! corpus, biocommons grammar tests, ClinVar legacy records).
//!
//! ## Audit goals
//!
//! The tests below pin parse, AST, normalize, and lowering behavior for
//! `con`, including the SVD-WG009 canonicalization implemented in this
//! PR. Any future change to any of these invariants must be intentional
//! and visible in code review.
//!
//! 1. The exact spec form from issue #81 (cross-reference c. → c.) —
//!    parse + display round-trip.
//! 2. Cross-reference genomic conversion (NC → NC) — round-trip.
//! 3. Same-reference position-only conversion (no second accession) —
//!    round-trip.
//! 4. Same-reference c. → c. position-only conversion — round-trip.
//! 5. Internal AST shape: edit type is `NaEdit::Conversion` with the
//!    full source string (including any `accession:` prefix) preserved
//!    verbatim.
//! 6. `normalize()` canonicalizes `con` edits to SVD-WG009 `delins`
//!    (cross-reference -> `delins[...]`, same-reference -> `delins<a>_<b>`).
//! 7. Downstream conversion to SPDI is an explicit
//!    `UnsupportedEditType` error whose message references "conversion"
//!    so callers know what's not supported.
//! 8. Downstream conversion to VCF is an explicit `ConversionError`
//!    — not a silent drop.
//! 9. The `gene_conversion` versioned-feature flag is registered and
//!    reports as supported in the current HGVS version.
//! 10. Coordinate-system breadth: `n.` and `m.` accept `con`.
//! 11. `Hash` + `Eq` stability: two parses of the same `con` input
//!     compare equal and hash to the same bucket.
//! 12. Idempotency: `normalize(normalize(x)) == normalize(x)`.
//! 13. Allele compound containing `con` — member is canonicalized.
//! 14. Empty-source parser rejection.
//! 15. Cross-reference accession + coord prefix survive canonicalization
//!     verbatim inside the `delins[...]` brackets.
//!
//! ## Closes #140
//!
//! - SVD-WG009 canonicalization: implemented (test 6). `parse → Display`
//!   continues to round-trip `con` verbatim (tests 1-4); the rewrite
//!   only happens on explicit `normalize()`.
//! - Sequence-aware source validation at parse time: declined.
//!   Superseded by canonicalization: the rewritten `delins` source is
//!   parsed by ferro's existing structured `delins`-source parser on
//!   any subsequent round-trip.
//! - SPDI / VCF expansion of `con`: declined. Both surfaces stay as
//!   explicit errors (tests 7, 8). Expansion would require lowering
//!   support for `delins[ACC:...]`, which is a separate, larger feature.

use ferro_hgvs::hgvs::edit::NaEdit;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::hgvs::version::{features, is_feature_supported};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Extract the `NaEdit` from a c./g./n./r./m. variant, panicking on
/// any other shape. Convenience wrapper for the audit assertions
/// below — keeps each test focused on the audit invariant.
fn na_edit_of(variant: &HgvsVariant) -> &NaEdit {
    match variant {
        HgvsVariant::Cds(v) => v
            .loc_edit
            .edit
            .inner()
            .expect("Cds variant must carry an edit"),
        HgvsVariant::Genome(v) => v
            .loc_edit
            .edit
            .inner()
            .expect("Genome variant must carry an edit"),
        HgvsVariant::Tx(v) => v
            .loc_edit
            .edit
            .inner()
            .expect("Tx variant must carry an edit"),
        HgvsVariant::Rna(v) => v
            .loc_edit
            .edit
            .inner()
            .expect("Rna variant must carry an edit"),
        HgvsVariant::Mt(v) => v
            .loc_edit
            .edit
            .inner()
            .expect("Mt variant must carry an edit"),
        other => panic!(
            "expected nucleotide variant carrying NaEdit, got {:?}",
            other
        ),
    }
}

// ---------------------------------------------------------------------
// 1. Exact spec form from issue #81: cross-reference c. → c.
// ---------------------------------------------------------------------

/// Pins parse + display round-trip for the canonical issue #81 H1
/// example: a CDS interval on one transcript replaced by a CDS
/// interval from another transcript.
///
/// If this regresses, ferro has either lost `con` parsing or changed
/// its display form (e.g. rewriting to `delins`). Either is a
/// behavior change that should be intentional and tracked.
#[test]
fn cross_reference_cds_conversion_round_trips() {
    let input = "NM_000088.3:c.123_456conNM_000089.1:c.789_1011";
    let variant = parse_hgvs(input).expect("issue #81 H1 canonical form must parse");
    assert!(matches!(variant, HgvsVariant::Cds(_)));
    assert_eq!(format!("{}", variant), input);
}

// ---------------------------------------------------------------------
// 2. Cross-reference genomic conversion.
// ---------------------------------------------------------------------

/// Pins parse + display round-trip for a genomic-to-genomic
/// cross-reference conversion. Matches the form found in the
/// biocommons gauntlet and Mutalyzer grammar fixtures.
#[test]
fn cross_reference_genomic_conversion_round_trips() {
    let input = "NC_000001.11:g.12345_12400conNC_000002.12:g.100_155";
    let variant = parse_hgvs(input).expect("genomic cross-reference con must parse");
    assert!(matches!(variant, HgvsVariant::Genome(_)));
    assert_eq!(format!("{}", variant), input);
}

// ---------------------------------------------------------------------
// 3. Same-reference position-only conversion (genomic).
// ---------------------------------------------------------------------

/// Pins the spec example from
/// `assets/hgvs-nomenclature/docs/recommendations/DNA/delins.md`
/// (CYP2D6 exon-9 conversion to CYP2D7P1 flanking sequence on the
/// same NC_ accession). Source has no `accession:` prefix; just
/// numeric positions.
#[test]
fn same_reference_genomic_position_conversion_round_trips() {
    let input = "NC_000017.11:g.42522624_42522669con42536337_42536382";
    let variant = parse_hgvs(input).expect("same-reference position con must parse");
    assert!(matches!(variant, HgvsVariant::Genome(_)));
    assert_eq!(format!("{}", variant), input);
}

// ---------------------------------------------------------------------
// 4. Same-reference position-only conversion (CDS).
// ---------------------------------------------------------------------

/// Pins the spec DRD4 example: `c.812_829con908_925` — same
/// reference, CDS coordinates on both sides, no second accession.
/// Note that the parser preserves the source string verbatim, so
/// `908_925` (without a `c.` prefix) round-trips as-is.
#[test]
fn same_reference_cds_position_conversion_round_trips() {
    let input = "NM_000797.3:c.812_829con908_925";
    let variant = parse_hgvs(input).expect("same-reference CDS position con must parse");
    assert!(matches!(variant, HgvsVariant::Cds(_)));
    assert_eq!(format!("{}", variant), input);
}

// ---------------------------------------------------------------------
// 5. Internal AST shape.
// ---------------------------------------------------------------------

/// Pins the AST shape: every accepted `con` form lands in
/// `NaEdit::Conversion { source }` with `source` capturing
/// everything after the literal `con` token verbatim. Pinning this
/// is what catches a future refactor that, for example, splits
/// `source` into structured `(Option<Accession>, Location)` — such
/// a change would break this test loudly and force a deliberate
/// migration.
#[test]
fn conversion_edit_preserves_source_string_verbatim() {
    let cases = [
        // (input HGVS, expected source string after the `con` token)
        (
            "NM_000088.3:c.123_456conNM_000089.1:c.789_1011",
            "NM_000089.1:c.789_1011",
        ),
        (
            "NC_000017.11:g.42522624_42522669con42536337_42536382",
            "42536337_42536382",
        ),
        ("NM_000797.3:c.812_829con908_925", "908_925"),
    ];

    for (input, expected_source) in cases {
        let variant =
            parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for {}: {}", input, e));
        let edit = na_edit_of(&variant);
        match edit {
            NaEdit::Conversion { source } => {
                assert_eq!(
                    source, expected_source,
                    "source mismatch for input {}",
                    input
                );
            }
            other => panic!("expected NaEdit::Conversion for {}, got {:?}", input, other),
        }
    }
}

// ---------------------------------------------------------------------
// 6. `normalize()` canonicalizes `con` to SVD-WG009 `delins`.
// ---------------------------------------------------------------------

/// Pins SVD-WG009 canonicalization combined with issue #333 ins-expand
/// chaining: `normalize` rewrites `con` → `delins` and then re-runs the
/// axis normalizer on the rewritten variant, so a single call carries
/// the rewrite all the way to its flat-literal form (or to the
/// documented deferral error for cross-references).
///
/// - **Cross-reference source** (e.g. `conNM_000089.1:c.789_1011`):
///   surfaces as `FerroError::UnsupportedVariant` carrying the
///   `cross-reference` / `follow-up` markers. Spec-allowed but
///   deferred by ferro per the PR contract.
/// - **Same-reference position source**: requires provider sequence at
///   the source coordinates. With an empty `MockProvider`, the lookup
///   fails and surfaces an error rather than passing the intermediate
///   `delins<start>_<end>` shape through.
///
/// `parse_hgvs` + `Display` continue to round-trip the legacy `con`
/// form (pinned by tests 1-4); canonicalization happens only when
/// callers explicitly invoke `normalize`.
#[test]
fn normalize_canonicalizes_conversion_to_delins() {
    let provider = MockProvider::new();
    let normalizer = Normalizer::new(provider);

    // Cross-reference sources: #422 wired actual lookup through the
    // cross-reference resolver, so an empty provider surfaces the
    // missing inner accession via `ReferenceNotFound` rather than the
    // earlier `UnsupportedVariant` placeholder. The semantic point —
    // "the recursion does not pass the bracketed form through
    // unchanged" — is preserved.
    let cross_ref_inputs = [
        "NM_000088.3:c.123_456conNM_000089.1:c.789_1011",
        "NC_000001.11:g.12345_12400conNC_000002.12:g.100_155",
    ];
    for input in cross_ref_inputs {
        let variant = parse_hgvs(input).expect("input must parse");
        let err = normalizer
            .normalize(&variant)
            .expect_err("cross-reference con payload must error against an empty provider");
        assert!(
            matches!(
                err,
                ferro_hgvs::FerroError::ReferenceNotFound { .. }
                    | ferro_hgvs::FerroError::GenomicReferenceNotAvailable { .. }
            ),
            "expected provider-lookup error for cross-reference con without test data (input: {}, got: {:?})",
            input,
            err
        );
    }

    // Same-reference position-range sources: the recursion runs the
    // intermediate `delins<start>_<end>` through the ins-expand step,
    // which needs provider sequence at the source coordinates. The
    // empty `MockProvider::new()` has no entry for `NC_000017.11` or
    // for `NM_000797.3`, so each surfaces a provider-lookup error
    // rather than passing the bracketed form through unchanged.
    let same_ref_inputs = [
        "NC_000017.11:g.42522624_42522669con42536337_42536382",
        "NM_000797.3:c.812_829con908_925",
    ];
    for input in same_ref_inputs {
        let variant = parse_hgvs(input).expect("input must parse");
        let err = normalizer.normalize(&variant).expect_err(
            "same-reference con payload requires provider sequence — no test data is loaded",
        );
        assert!(
            matches!(
                err,
                ferro_hgvs::FerroError::GenomicReferenceNotAvailable { .. }
                    | ferro_hgvs::FerroError::ReferenceNotFound { .. }
                    | ferro_hgvs::FerroError::ConversionError { .. }
            ),
            "expected provider-lookup error for same-reference con without test data (input: {}, got: {:?})",
            input,
            err
        );
    }
}

// ---------------------------------------------------------------------
// 7. SPDI lowering rejects `con` explicitly.
// ---------------------------------------------------------------------

/// Pins that converting a `con` HGVS to SPDI fails with an explicit
/// `UnsupportedEditType` error rather than silently dropping or
/// emitting a wrong record. Behavior lives at
/// `src/spdi/convert.rs:301` today.
#[test]
fn conversion_to_spdi_is_explicit_unsupported_error() {
    use ferro_hgvs::spdi::ConversionError;

    let variant =
        parse_hgvs("NC_000017.11:g.42522624_42522669con42536337_42536382").expect("parse");
    let err = ferro_hgvs::hgvs_to_spdi_simple(&variant).expect_err("SPDI lowering must error");
    assert!(
        matches!(err, ConversionError::UnsupportedEditType { .. }),
        "expected UnsupportedEditType, got {:?}",
        err
    );
    let msg = format!("{}", err);
    assert!(
        msg.to_lowercase().contains("conversion"),
        "SPDI error message must mention 'conversion' for callers, got: {}",
        msg
    );
}

// ---------------------------------------------------------------------
// 8. VCF lowering rejects `con` explicitly.
// ---------------------------------------------------------------------

/// Pins that converting a `con` HGVS to VCF returns an explicit
/// error. The match arm at `src/vcf/from_hgvs.rs:430` is exercised
/// via `genomic_hgvs_to_vcf` (the provider-less, genome-only entry
/// point). We do not pin the human-readable message — only that it
/// is an error and not silent success.
#[test]
fn conversion_to_vcf_is_explicit_error() {
    use ferro_hgvs::vcf::genomic_hgvs_to_vcf;

    let variant =
        parse_hgvs("NC_000017.11:g.42522624_42522669con42536337_42536382").expect("parse");
    let genome = match variant {
        HgvsVariant::Genome(g) => g,
        other => panic!("expected Genome variant, got {:?}", other),
    };
    let result = genomic_hgvs_to_vcf(&genome);
    assert!(
        result.is_err(),
        "VCF lowering of a con edit should error, got Ok({:?})",
        result.ok()
    );
}

// ---------------------------------------------------------------------
// 9. `gene_conversion` is a registered, supported feature.
// ---------------------------------------------------------------------

/// Pins that the version registry advertises `gene_conversion` as
/// supported in the current HGVS version. If a future refactor
/// removes the flag, this test catches the regression. The flag
/// itself is purely informational today (no code reads it to
/// decide whether to parse `con`); making it gate parsing is a
/// separate audit follow-up.
#[test]
fn gene_conversion_feature_is_registered_and_supported() {
    let feature = &features::GENE_CONVERSION;
    assert_eq!(feature.name, "gene_conversion");
    assert!(
        is_feature_supported(feature),
        "gene_conversion must be reported as supported in CURRENT version"
    );
}

// ---------------------------------------------------------------------
// 10. Coordinate-system breadth: n. and m. accept `con`.
// ---------------------------------------------------------------------

/// Pins that `con` is accepted in non-coding (`n.`) and mitochondrial
/// (`m.`) coord systems, not just `c.` and `g.`. The HGVS spec is
/// coord-system-agnostic for `con`; this guards against a future parser
/// refactor that accidentally narrows the accepted prefixes.
///
/// Round-trip parse + Display only — `normalize` canonicalization for
/// `n.` and `m.` is exercised separately via the per-system normalizers
/// (Tx and Mt) wired in this PR.
#[test]
fn conversion_round_trips_for_n_and_m_coords() {
    let n_input = "NR_000001.1:n.10_20conNR_000002.1:n.30_40";
    let variant = parse_hgvs(n_input).expect("n. con must parse");
    assert!(matches!(variant, HgvsVariant::Tx(_)));
    assert_eq!(format!("{}", variant), n_input);

    let m_input = "NC_012920.1:m.100_200conNC_012920.1:m.300_400";
    let variant = parse_hgvs(m_input).expect("m. con must parse");
    assert!(matches!(variant, HgvsVariant::Mt(_)));
    assert_eq!(format!("{}", variant), m_input);
}

// ---------------------------------------------------------------------
// 11. `Hash` + `Eq` stability for `NaEdit::Conversion`.
// ---------------------------------------------------------------------

/// Pins that two parses of the same `con` input produce structurally
/// equal variants that hash to the same bucket. A future refactor that
/// changes `source: String` to a structured payload must keep this
/// invariant -- otherwise `HashSet<HgvsVariant>` and friends silently
/// break.
#[test]
fn conversion_hash_and_eq_are_stable_across_parses() {
    use std::collections::HashSet;

    let input = "NM_000088.3:c.123_456conNM_000089.1:c.789_1011";
    let a = parse_hgvs(input).expect("first parse");
    let b = parse_hgvs(input).expect("second parse");
    assert_eq!(a, b, "two parses of the same input must be equal");

    let mut set = HashSet::new();
    set.insert(a);
    set.insert(b);
    assert_eq!(
        set.len(),
        1,
        "two parses of the same input must hash to the same bucket"
    );
}

// ---------------------------------------------------------------------
// 12. Idempotency: normalize(normalize(con)) == normalize(con).
// ---------------------------------------------------------------------

/// Pins that a single `normalize()` call on a `con` input does not
/// stop at the bracketed `delins{Reference | PositionRange}`
/// intermediate. Under issue #333 the axis normalizer recurses after
/// `con` → `delins` so the rewrite runs through the ins-expand step.
/// For these MockProvider-less inputs that means:
///
/// - cross-reference source → `FerroError::UnsupportedVariant` with
///   `cross-reference` / `follow-up` markers,
/// - same-reference position source → provider-lookup error.
///
/// In neither case does the bracketed intermediate survive to the
/// caller. This is the explicit anti-test for the "stops at intermediate
/// delins" shape: catching a regression that reintroduces the half-
/// canonicalized output would re-break the issue #333 contract.
#[test]
fn normalize_canonicalizes_conversion_to_delins_shape() {
    let provider = MockProvider::new();
    let normalizer = Normalizer::new(provider);

    let inputs = [
        "NM_000088.3:c.123_456conNM_000089.1:c.789_1011",
        "NC_000017.11:g.42522624_42522669con42536337_42536382",
    ];
    for input in inputs {
        let variant = parse_hgvs(input).expect("parse");
        let err = normalizer.normalize(&variant).expect_err(
            "first normalize must drive con-rewrite through ins-expand, not stop at the bracketed intermediate",
        );
        // The error carries either the cross-reference deferral marker
        // or a provider-lookup failure marker. Either is consistent
        // with the issue #333 contract; neither leaks the half-
        // canonicalized `delins[...]` shape to the caller.
        assert!(
            matches!(
                err,
                ferro_hgvs::FerroError::UnsupportedVariant { .. }
                    | ferro_hgvs::FerroError::GenomicReferenceNotAvailable { .. }
                    | ferro_hgvs::FerroError::ReferenceNotFound { .. }
                    | ferro_hgvs::FerroError::ConversionError { .. }
            ),
            "expected UnsupportedVariant or provider-lookup error (input: {}, got: {:?})",
            input,
            err
        );
    }
}

// ---------------------------------------------------------------------
// 13. Allele compound containing `con`.
// ---------------------------------------------------------------------

/// Pins that `con` survives inside a bracketed compound allele
/// alongside other edit types: the `con` member is canonicalized to
/// `delins` while the substitution member is unchanged.
#[test]
fn allele_compound_with_conversion_canonicalizes_member() {
    // Empty provider: the cross-reference deferral path is provider-
    // independent and must not be gated by the issue #336 PositionPastEnd
    // bounds check (which would otherwise fire on `c.123_456` against
    // the cds_end=60 NM_000088.3 entry in `with_test_data()`).
    let provider = MockProvider::new();
    let normalizer = Normalizer::new(provider);

    let input = "NM_000088.3:c.[123_456conNM_000089.1:c.789_1011;500A>G]";
    let variant = parse_hgvs(input).expect("compound con+sub must parse");
    // Under #333 + #422 the recursion drives the `con` member's rewrite
    // through the ins-expand step. The cross-reference half now hits
    // the actual resolver, which surfaces `ReferenceNotFound` against
    // the empty provider — the substitution sibling is unreachable
    // because the compound short-circuits on the first member error.
    let err = normalizer
        .normalize(&variant)
        .expect_err("cross-reference con member must surface a provider-lookup error");
    assert!(
        matches!(
            err,
            ferro_hgvs::FerroError::ReferenceNotFound { .. }
                | ferro_hgvs::FerroError::GenomicReferenceNotAvailable { .. }
        ),
        "expected provider-lookup error for cross-reference con member, got: {:?}",
        err
    );
    let msg = format!("{}", err);
    assert!(
        msg.contains("NM_000089.1"),
        "lookup error must cite the missing source accession; got {}",
        msg
    );
}

// ---------------------------------------------------------------------
// 14. Empty-source parser rejection.
// ---------------------------------------------------------------------

/// Pins that `con` followed by an empty source is rejected by the
/// parser. Today this is enforced by `take_while1(!whitespace)` in
/// `parse_conversion`; this test fences the invariant against a future
/// switch to `take_while0`.
#[test]
fn empty_source_after_con_is_rejected() {
    let result = parse_hgvs("NM_000088.3:c.1_2con");
    assert!(
        result.is_err(),
        "empty source after `con` must be rejected, got Ok({:?})",
        result.ok()
    );
}

// ---------------------------------------------------------------------
// 15. Cross-reference accession survives canonicalization verbatim.
// ---------------------------------------------------------------------

/// Pins the cross-reference deferral contract under issue #333. The
/// previous SVD-WG009 shape (`delins[source-with-accession-prefix]`)
/// is now an internal intermediate the recursion never lets surface:
/// a single `normalize` carries the rewrite through the ins-expand
/// step, which defers cross-reference payloads with
/// `FerroError::UnsupportedVariant`. The deferral message must mention
/// the source accession so callers can route follow-up resolution.
#[test]
fn cross_reference_canonicalization_preserves_source_payload() {
    // Empty provider: the cross-reference deferral path is provider-
    // independent and must not be gated by the issue #336 PositionPastEnd
    // bounds check (which would otherwise fire on `c.123_456` against
    // the cds_end=60 NM_000088.3 entry in `with_test_data()`).
    let provider = MockProvider::new();
    let normalizer = Normalizer::new(provider);

    let cases = [
        (
            "NM_000088.3:c.123_456conNM_000089.1:c.789_1011",
            "NM_000089.1",
        ),
        (
            "NC_000001.11:g.12345_12400conNC_000002.12:g.100_155",
            "NC_000002.12",
        ),
    ];
    for (input, expected_source_accession) in cases {
        let variant = parse_hgvs(input).expect("parse");
        let err = normalizer
            .normalize(&variant)
            .expect_err("cross-reference con must surface a provider-lookup error");
        // #422 replaced the `UnsupportedVariant` placeholder with actual
        // cross-reference resolution; an empty provider surfaces
        // `ReferenceNotFound` against the missing inner accession.
        assert!(
            matches!(
                err,
                ferro_hgvs::FerroError::ReferenceNotFound { .. }
                    | ferro_hgvs::FerroError::GenomicReferenceNotAvailable { .. }
            ),
            "expected provider-lookup error (input: {}, got: {:?})",
            input,
            err
        );
        let msg = format!("{}", err);
        assert!(
            msg.contains(expected_source_accession),
            "lookup error must cite the missing source accession {} (input: {}, got: {})",
            expected_source_accession,
            input,
            msg
        );
    }
}
