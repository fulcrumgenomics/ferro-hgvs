//! Audit tests for mitochondrial heteroplasmy notation acceptance.
//!
//! Tracks issue #81, item F2: heteroplasmy notation is common in
//! ClinVar but sparse in the HGVS spec; this file pins what ferro
//! currently accepts and emits, and which real-world inputs it
//! rejects.
//!
//! ## Background
//!
//! Heteroplasmy is the presence of multiple mitochondrial DNA
//! sequence states across mitochondria within the same cell. The
//! HGVS spec does **not** define a heteroplasmy-specific notation:
//! - The mtDNA `m.` prefix is the only mt-specific syntax (see
//!   `assets/hgvs-nomenclature/docs/background/refseq.md` and
//!   `consultation/SVD-WG006.md`).
//! - The closest spec analogue to "fraction of cells / fraction of
//!   molecules carrying the variant" is the **mosaic** (`var/=`) and
//!   **chimeric** (`var//=`) markers from
//!   `recommendations/DNA/substitution.md`. These are presence-only
//!   markers; the spec explicitly does **not** carry an allele
//!   fraction in the HGVS string itself.
//! - ClinVar reports heteroplasmy as a separate metadata field; the
//!   HGVS string is plain `NC_012920.1:m.<pos><ref>><alt>`.
//!
//! ## What this file pins
//!
//! 1. Plain ClinVar mt variants (the only well-supported shape) parse
//!    and round-trip exactly.
//! 2. Allele-bracket compounds combining several mt variants parse
//!    and emit in the spec-preferred compact form.
//! 3. Spec mosaic/chimeric markers (`var/var2`, `var//var2`) parse
//!    when each side is a fully qualified HGVS variant.
//! 4. The spec mosaic compact in-position form (`c.85=/T>C`,
//!    documented in `recommendations/DNA/substitution.md`) is **not
//!    supported** for any coordinate type, including `m.`. The
//!    same compact form is documented for deletion and duplication
//!    (`recommendations/DNA/deletion.md` lines 59-60,
//!    `recommendations/DNA/duplication.md` lines 56-57); those edit
//!    types are also rejected on `m.` and pinned here.
//! 5. The non-spec ClinVar prose multi-allelic shorthand
//!    `m.3243A>G/T` is rejected.
//! 6. Allele-fraction annotations (e.g. `[level=70%]`, `(80%)`) are
//!    rejected — and should be: HGVS does not encode allele fraction
//!    in the variant string. Heteroplasmy level belongs in
//!    accompanying metadata, not the HGVS expression.
//! 7. The bare `<acc>:m.<pos>=` whole-position identity edit
//!    (the LHS that the spec compact mosaic form would split off)
//!    parses standalone today.
//! 8. The non-spec workaround `m.<pos>=/m.<pos><ref>><alt>` (RHS
//!    has its own `m.` type prefix) is *accepted* by ferro's
//!    accession-inheritance branch. Pinned so #133's compact-form
//!    fix coexists with this existing acceptance.
//! 9. Slash mosaic/chimeric dispatchers are coord-agnostic — the
//!    same construction works on `c.` as on `m.`. Pinned so this
//!    audit reads as a slash-dispatcher contract surfaced through
//!    mt inputs, not as mt-only behavior.
//! 10. Heteroplasmy-prose copy-paste hazards (internal whitespace,
//!     non-ASCII A) are rejected.
//!
//! Items (4) and (5) are real-world gaps; (4) is the spec form, (5)
//! is the ClinVar prose form. Follow-up issue #133 tracks adding
//! support for the spec compact mosaic form.
//!
//! ## Boundary with sibling F1 (`tests/mito_circular_audit.rs`)
//!
//! Compound `m.[…]` allele forms (cis `[a;b]`, trans `[a];[b]`,
//! homozygous shorthand `[a](;)`) are all rejected today — pinned
//! by F1 §8, not duplicated here. Heteroplasmy may be expressed via
//! those forms in a future grammar surface, but ownership of the
//! decision lives in F1 / #129.

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

// ---------------------------------------------------------------------------
// 1. Plain ClinVar mt variants — must parse and round-trip exactly.
// ---------------------------------------------------------------------------

/// Every `m.` entry in `tests/fixtures/validation/clinvar.json` (as of
/// this audit) is a plain substitution against `NC_012920.1`. They
/// must all parse and round-trip without rewriting.
#[test]
fn plain_clinvar_mt_substitutions_round_trip() {
    // Verbatim from tests/fixtures/validation/clinvar.json (m. entries).
    let inputs = [
        "NC_012920.1:m.3243A>G",  // MELAS
        "NC_012920.1:m.8344A>G",  // MERRF
        "NC_012920.1:m.11778G>A", // LHON
        "NC_012920.1:m.14484T>C", // LHON
        "NC_012920.1:m.3460G>A",  // LHON
        "NC_012920.1:m.8993T>G",  // NARP/Leigh
        "NC_012920.1:m.1A>G",     // origin-adjacent edge
        "NC_012920.1:m.16569G>A", // last-base edge
    ];
    for input in inputs {
        let parsed =
            parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
        assert!(
            matches!(parsed, HgvsVariant::Mt(_)),
            "expected Mt variant for `{input}`, got {parsed:?}"
        );
        assert_eq!(format!("{parsed}"), input, "round-trip mismatch");
    }
}

/// Indels and the gene-symbol-suffixed form are also accepted on `m.`
/// (the latter is an explicit spec exception in `refseq.md`).
#[test]
fn mt_indels_and_gene_suffixed_forms_round_trip() {
    // (input, expected_emitted_form). Gene-symbol suffix is dropped on
    // emission today; this pins that current behavior.
    let cases = [
        ("NC_012920.1:m.8470_8482del", "NC_012920.1:m.8470_8482del"),
        ("NC_012920.1:m.302_303insC", "NC_012920.1:m.302_303insC"),
        ("NC_012920.1(MT-TL1):m.3243A>G", "NC_012920.1:m.3243A>G"),
        ("NC_012920.1(MT-ND4):m.11778G>A", "NC_012920.1:m.11778G>A"),
    ];
    for (input, expected) in cases {
        let parsed =
            parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
        assert_eq!(format!("{parsed}"), expected, "for input `{input}`");
    }
}

// ---------------------------------------------------------------------------
// 2. Allele compounds — the spec way to express co-occurrence on mtDNA.
// ---------------------------------------------------------------------------

/// Co-occurring mt variants are expressed via the bracket allele form
/// `[var1;var2]`. Per spec, mtDNA cells carry a single chromosome
/// type, so cis is the natural reading.
#[test]
fn mt_cis_allele_compound_round_trips_compact() {
    let input = "[NC_012920.1:m.3243A>G;NC_012920.1:m.8344A>G]";
    let parsed = parse_hgvs(input).expect("cis allele of mt variants must parse");
    if let HgvsVariant::Allele(allele) = &parsed {
        assert_eq!(allele.phase, AllelePhase::Cis);
        assert_eq!(allele.variants.len(), 2);
    } else {
        panic!("expected Allele variant, got {parsed:?}");
    }
    // Compact form is the spec-preferred output when both sides share
    // the accession and variant type.
    assert_eq!(
        format!("{parsed}"),
        "NC_012920.1:m.[3243A>G;8344A>G]",
        "expected compact bracketed form"
    );
}

// ---------------------------------------------------------------------------
// 3. Spec mosaic/chimeric markers on mt — fully qualified both sides.
// ---------------------------------------------------------------------------

/// Spec `var1/var2` mosaic syntax: works on `m.` when both sides are
/// fully qualified HGVS variants. This is the closest spec form to a
/// heteroplasmy "two states observed at one site" notation.
#[test]
fn mt_mosaic_two_fully_qualified_variants_round_trips() {
    let input = "NC_012920.1:m.3243A>G/NC_012920.1:m.3243A>T";
    let parsed = parse_hgvs(input).expect("mosaic mt allele must parse");
    if let HgvsVariant::Allele(allele) = &parsed {
        assert_eq!(allele.phase, AllelePhase::Mosaic);
        assert_eq!(allele.variants.len(), 2);
    } else {
        panic!("expected Allele variant, got {parsed:?}");
    }
    // No accession compaction is performed across the slash today;
    // both sides are emitted in full.
    assert_eq!(format!("{parsed}"), input);
}

/// The `var/=` shorthand (spec: "the second allele is the reference
/// sequence") parses on mt. The current emitted form expands `=` to a
/// concrete identity variant `NC_012920.1:m.1=` rather than the more
/// compact `m.=` — pinned here so a future cleanup is intentional.
#[test]
fn mt_mosaic_with_reference_shorthand_parses_with_known_render_quirk() {
    let input = "NC_012920.1:m.3243A>G/=";
    let parsed = parse_hgvs(input).expect("mosaic with `=` shorthand must parse");
    if let HgvsVariant::Allele(allele) = &parsed {
        assert_eq!(allele.phase, AllelePhase::Mosaic);
        assert_eq!(allele.variants.len(), 2);
    } else {
        panic!("expected Allele variant, got {parsed:?}");
    }
    // Pinned current rendering. A more spec-aligned form would emit
    // `NC_012920.1:m.3243A>G/=`; this test will fail (correctly) the
    // day that cleanup lands.
    assert_eq!(
        format!("{parsed}"),
        "NC_012920.1:m.3243A>G/NC_012920.1:m.1="
    );
}

// ---------------------------------------------------------------------------
// 4. Spec compact mosaic in-position form `c.85=/T>C` — REJECTED.
//
// Per `recommendations/DNA/substitution.md`:
//   `LRG_199t1:c.85=/T>C` — a mosaic case where at position 85
//   besides the normal sequence (a T, described as `=`) also
//   chromosomes are found containing a C.
// This is the spec's documented heteroplasmy/mosaic-friendly compact
// form. It is rejected today for every coordinate system, including
// `m.`, because the parser splits on `/` and tries to interpret the
// right-hand side as a full variant.
// ---------------------------------------------------------------------------

#[test]
fn spec_compact_mosaic_form_rejected_on_mt() {
    // Spec example shape, transposed to mt: at position 3243 the
    // reference (A) is observed alongside a G allele.
    let input = "NC_012920.1:m.3243=/A>G";
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "AUDIT: spec compact mosaic form is currently rejected; \
         update this test (and tracking issue) when support lands. \
         Got: {:?}",
        result
    );
}

#[test]
fn spec_compact_chimeric_form_rejected_on_mt() {
    // Same shape with the chimeric `//` separator from
    // substitution.md (`LRG_199t1:c.85=//T>C`).
    let input = "NC_012920.1:m.3243=//A>G";
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "AUDIT: spec compact chimeric form is currently rejected; \
         update this test (and tracking issue) when support lands. \
         Got: {:?}",
        result
    );
}

// ---------------------------------------------------------------------------
// 5. Non-spec ClinVar prose multi-allelic shorthand `m.3243A>G/T` — REJECTED.
//
// Some ClinVar release notes and submitter prose describe a mt site
// where two alternative alleles are observed as `m.3243A>G/T`. This
// is **not** a HGVS spec form: it elides the second-allele reference
// base and lacks an accession on the right. Today it is rejected by
// the parser (single-slash dispatcher tries to parse the right side
// as a fresh variant beginning with a base letter).
// ---------------------------------------------------------------------------

#[test]
fn clinvar_prose_multi_allelic_shorthand_rejected() {
    let inputs = ["NC_012920.1:m.3243A>G/T", "NC_012920.1:m.3243A>G/C"];
    for input in inputs {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "AUDIT: prose multi-allelic shorthand `{input}` is \
             currently rejected. If a follow-up issue elects to \
             accept this form (rewriting it to a spec-mosaic), \
             update this test. Got: {:?}",
            result
        );
    }
}

// ---------------------------------------------------------------------------
// 6. Allele-fraction annotations are out-of-scope and correctly rejected.
//
// HGVS does not encode allele fraction inside the variant string;
// heteroplasmy level (e.g. 70%) belongs in accompanying metadata
// (VCF FORMAT/AF, ClinVar's separate fields, etc.). These cases are
// pinned as **must-reject** so that a future change which
// accidentally accepts them — and tries to round-trip the suffix —
// trips the test.
// ---------------------------------------------------------------------------

#[test]
fn allele_fraction_annotations_must_be_rejected() {
    let inputs = [
        "NC_012920.1:m.3243A>G[level=70%]",
        "NC_012920.1:m.3243A>G(80%)",
    ];
    for input in inputs {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "HGVS does not carry allele fraction; `{input}` must be \
             rejected. Got: {:?}",
            result
        );
    }
}

// ---------------------------------------------------------------------------
// 7. Spec compact mosaic form on del/dup — REJECTED on mt.
//
// The spec documents the same `<pos>=/<edit>` shape for deletion and
// duplication, not just substitution:
//   - `recommendations/DNA/deletion.md` lines 59-60:
//       `NC_000023.11:g.33344590_33344592=/del`   (mosaic)
//       `NC_000023.11:g.33344590_33344592=//del`  (chimeric)
//   - `recommendations/DNA/duplication.md` lines 56-57:
//       `NC_000023.11:g.33344590_33344592=/dup`   (mosaic)
//       `NC_000023.11:g.33344590_33344592=//dup`  (chimeric)
// Pinned for `m.` so the future #133 fix has a complete contract
// (substitution, deletion, duplication — all four mosaic+chimeric
// pairs).
// ---------------------------------------------------------------------------

#[test]
fn spec_compact_mosaic_deletion_form_rejected_on_mt() {
    let inputs = [
        "NC_012920.1:m.8470_8482=/del",
        "NC_012920.1:m.8470_8482=//del",
    ];
    for input in inputs {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "AUDIT: spec compact mosaic/chimeric deletion form \
             `{input}` is currently rejected; update this test \
             when #133 lands. Got: {:?}",
            result
        );
    }
}

#[test]
fn spec_compact_mosaic_duplication_form_rejected_on_mt() {
    let inputs = ["NC_012920.1:m.302_303=/dup", "NC_012920.1:m.302_303=//dup"];
    for input in inputs {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "AUDIT: spec compact mosaic/chimeric duplication form \
             `{input}` is currently rejected; update this test \
             when #133 lands. Got: {:?}",
            result
        );
    }
}

// ---------------------------------------------------------------------------
// 8. Bare `m.<pos>=` whole-position identity edit — PARSES.
//
// The LHS of the spec compact mosaic form is a whole-position
// identity edit. Ferro already parses it standalone today; pinned
// here so any rework of the slash dispatcher (#133) does not
// regress this primitive.
// ---------------------------------------------------------------------------

#[test]
fn mt_position_identity_edit_round_trips() {
    let cases = [
        ("NC_012920.1:m.3243=", "NC_012920.1:m.3243="),
        ("NC_012920.1:m.3243A=", "NC_012920.1:m.3243A="),
        ("NC_012920.1:m.8470_8482=", "NC_012920.1:m.8470_8482="),
    ];
    for (input, expected) in cases {
        let parsed =
            parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
        assert!(
            matches!(parsed, HgvsVariant::Mt(_)),
            "expected Mt variant for `{input}`, got {parsed:?}"
        );
        assert_eq!(format!("{parsed}"), expected, "round-trip mismatch");
    }
}

// ---------------------------------------------------------------------------
// 9. Type-prefixed inheritance workaround `m.<pos>=/m.<pos><ref>><alt>`.
//
// Probe-confirmed: when the right-hand side carries its own `m.`
// type prefix, the slash dispatcher's accession-inheritance branch
// in `parse_mosaic_allele` succeeds. This is **not** the spec form
// (the spec writes the right side as a bare edit, no type prefix:
// `c.85=/T>C`), but it is what ferro accepts today.
//
// Pinned so any #133 implementation that adds the bare-edit form
// continues to accept this longer prefixed form too — losing it
// would silently break callers already relying on this quirk.
// ---------------------------------------------------------------------------

#[test]
fn mt_compact_mosaic_with_type_prefixed_rhs_parses_today() {
    let input = "NC_012920.1:m.3243=/m.3243A>G";
    let parsed = parse_hgvs(input)
        .unwrap_or_else(|e| panic!("AUDIT: type-prefixed RHS mosaic must parse today: {e}"));
    if let HgvsVariant::Allele(allele) = &parsed {
        assert_eq!(allele.phase, AllelePhase::Mosaic);
        assert_eq!(allele.variants.len(), 2);
    } else {
        panic!("expected Allele variant, got {parsed:?}");
    }
    // Inherited accession is materialized on the RHS during emission.
    assert_eq!(
        format!("{parsed}"),
        "NC_012920.1:m.3243=/NC_012920.1:m.3243A>G"
    );
}

// ---------------------------------------------------------------------------
// 10. Slash dispatcher is coord-agnostic.
//
// The mosaic/chimeric markers do not appear in any coord-system
// section of the spec other than via the compact forms; the slash
// dispatcher in `parse_mosaic_allele` / `parse_chimeric_allele`
// does not branch on coord type. Pinned with one positive `c.`
// case to make this contract explicit: the F2 audit is pinning
// slash-dispatcher behavior surfaced through mt inputs, not
// mt-specific parser logic.
// ---------------------------------------------------------------------------

#[test]
fn slash_dispatcher_is_coord_agnostic_between_c_and_mt() {
    // Two fully qualified variants on c. — must parse as a Mosaic
    // allele, the same shape as the mt case in section 3.
    let c_input = "NM_000088.3:c.100A>G/NM_000088.3:c.100A>T";
    let c_parsed = parse_hgvs(c_input).expect("c. mosaic must parse");
    if let HgvsVariant::Allele(allele) = &c_parsed {
        assert_eq!(allele.phase, AllelePhase::Mosaic);
        assert_eq!(allele.variants.len(), 2);
    } else {
        panic!("expected Allele variant on c., got {c_parsed:?}");
    }
    // The mt counterpart is already pinned in section 3; this test
    // pins that the c. form parses through the same dispatcher with
    // matching shape, so any change in #133 to the dispatcher must
    // affect both coord systems consistently.
}

// ---------------------------------------------------------------------------
// 11. Heteroplasmy-prose copy-paste hazards — must reject.
//
// Real-world heteroplasmy submissions sometimes carry these
// degraded shapes from manual transcription:
//   - internal whitespace inside the prefix (`m. 3243A>G`)
//   - non-ASCII Cyrillic А mistaken for Latin A (`m.3243А>G`)
// Both are rejected today and must stay rejected — silent
// acceptance would be a real bug because the parsed allele would
// not match the user's intent.
// ---------------------------------------------------------------------------

#[test]
fn heteroplasmy_prose_copy_paste_hazards_must_be_rejected() {
    // Non-exhaustive but covers the two hazards seen in submitter
    // prose. Each input intentionally carries one defect.
    let inputs = [
        "NC_012920.1:m. 3243A>G", // space after type prefix
        "NC_012920.1: m.3243A>G", // space after accession colon
        "NC_012920.1:m.3243А>G",  // Cyrillic А (U+0410), not Latin A
        "NC_012920.1:m.3243A>С",  // Cyrillic С (U+0421), not Latin C
    ];
    for input in inputs {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "AUDIT: copy-paste hazard `{input}` must reject. Got: {:?}",
            result
        );
    }
}
