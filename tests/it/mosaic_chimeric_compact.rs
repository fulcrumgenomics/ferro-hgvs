//! Tests for HGVS spec compact mosaic/chimeric form.
//!
//! Closes #133. The spec form is `<acc>:<type>.<pos>=/<edit>` (mosaic)
//! and `<acc>:<type>.<pos>=//<edit>` (chimeric). Documented for
//! substitution (`recommendations/DNA/substitution.md` lines 30-31),
//! deletion (`recommendations/DNA/deletion.md` lines 59-60), and
//! duplication (`recommendations/DNA/duplication.md` lines 56-57).

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

fn hash_of(variant: &HgvsVariant) -> u64 {
    let mut h = DefaultHasher::new();
    variant.hash(&mut h);
    h.finish()
}

/// Smoke: substitution compact mosaic on `m.` parses.
#[test]
fn compact_mosaic_substitution_mt_parses() {
    let parsed =
        parse_hgvs("NC_012920.1:m.3243=/A>G").expect("spec compact mosaic substitution must parse");
    let HgvsVariant::Allele(allele) = &parsed else {
        panic!("expected Allele, got {parsed:?}");
    };
    assert_eq!(allele.phase, AllelePhase::Mosaic);
    assert_eq!(allele.variants.len(), 2);
}

/// Smoke: substitution compact chimeric on `m.` parses.
#[test]
fn compact_chimeric_substitution_mt_parses() {
    let parsed = parse_hgvs("NC_012920.1:m.3243=//A>G")
        .expect("spec compact chimeric substitution must parse");
    let HgvsVariant::Allele(allele) = &parsed else {
        panic!("expected Allele, got {parsed:?}");
    };
    assert_eq!(allele.phase, AllelePhase::Chimeric);
    assert_eq!(allele.variants.len(), 2);
}

/// Spec compact form round-trips across all NA coord systems.
///
/// Each input is the canonical spec compact form; parse + Display
/// must yield the same string verbatim.
#[test]
fn compact_form_round_trips_across_coord_systems() {
    let inputs = [
        // Substitution — mosaic + chimeric, one per coord system.
        "NC_000023.11:g.33344590=/A>G",
        "NC_000023.11:g.33344590=//A>G",
        "NM_000088.3:c.85=/T>C",
        "NM_000088.3:c.85=//T>C",
        "NR_002196.2:n.100=/A>G",
        "NR_002196.2:n.100=//A>G",
        // RNA renders nucleotides in lowercase per the HGVS spec.
        "NR_002196.2:r.100=/a>g",
        "NR_002196.2:r.100=//a>g",
        "NC_012920.1:m.3243=/A>G",
        "NC_012920.1:m.3243=//A>G",
        // Deletion — mosaic + chimeric.
        "NC_000023.11:g.33344590_33344592=/del",
        "NC_000023.11:g.33344590_33344592=//del",
        "NM_000088.3:c.79_80=/del",
        "NC_012920.1:m.8470_8482=/del",
        // Duplication — mosaic + chimeric.
        "NC_000023.11:g.33344590_33344592=/dup",
        "NC_000023.11:g.33344590_33344592=//dup",
        "NM_000088.3:c.79_80=/dup",
        "NC_012920.1:m.302_303=/dup",
    ];
    for input in inputs {
        let parsed =
            parse_hgvs(input).unwrap_or_else(|e| panic!("compact form must parse `{input}`: {e}"));
        assert!(
            matches!(parsed, HgvsVariant::Allele(_)),
            "expected Allele for `{input}`, got {parsed:?}"
        );
        assert_eq!(
            format!("{parsed}"),
            input,
            "round-trip mismatch for `{input}`"
        );
    }
}

/// Protein compact mosaic `<acc>:p.<pos>=/<edit>` (the `=/` somatic
/// and/or form, HGVS general.md). The bare RHS is either an alternative
/// amino acid (a substitution inheriting the LHS reference + position) or
/// a bare protein edit (`del`, `dup`).
#[test]
fn protein_compact_mosaic_round_trips() {
    for input in [
        "LRG_199p1:p.Trp24=/Cys",
        "NP_003997.1:p.Val7=/del",
        "NP_003997.1:p.Val7=/dup",
    ] {
        let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("must parse `{input}`: {e}"));
        assert!(
            matches!(parsed, HgvsVariant::Allele(_)),
            "expected Allele for `{input}`, got {parsed:?}"
        );
        assert_eq!(format!("{parsed}"), input, "round-trip for `{input}`");
    }
}

/// A bare alternative amino acid RHS only forms a substitution when the LHS
/// is a single point. A range identity LHS (`p.Trp24_Cys26=`) must not be
/// coerced into a substitution carrying the whole interval — the bare-AA
/// compact form is rejected so the caller falls through to its error path.
#[test]
fn protein_compact_mosaic_bare_aa_requires_point_lhs() {
    for input in [
        "NP_003997.1:p.Trp24_Cys26=/Gly",
        "NP_003997.1:p.Val7_Leu9=/Cys",
    ] {
        assert!(
            parse_hgvs(input).is_err(),
            "range-LHS bare-AA compact form must be rejected, but `{input}` parsed"
        );
    }
}

/// A single-base substitution RHS (`A>G`) only forms a compact mosaic when
/// the inherited LHS identity is a point. A range identity LHS
/// (`c.79_80=`) must not be coerced into a substitution spanning the whole
/// interval — the compact form is rejected. del / dup span ranges fine.
#[test]
fn na_compact_mosaic_substitution_requires_point_lhs() {
    for input in [
        "NM_000088.3:c.79_80=/A>G",
        "NM_000088.3:c.79_80=//A>G",
        "NC_000023.11:g.33344590_33344592=/A>G",
        "NC_012920.1:m.8470_8482=/A>G",
    ] {
        assert!(
            parse_hgvs(input).is_err(),
            "range-LHS substitution compact form must be rejected, but `{input}` parsed"
        );
    }
    // Range del / dup remain valid (they legitimately span the interval).
    for input in ["NM_000088.3:c.79_80=/del", "NM_000088.3:c.79_80=/dup"] {
        assert!(
            parse_hgvs(input).is_ok(),
            "range-LHS del/dup compact form must still parse: `{input}`"
        );
    }
}

/// RNA compact mosaic/chimeric forms must preserve lowercase bases on the
/// bare-edit RHS — HGVS renders RNA in lowercase. Previously the bare-edit
/// RHS routed through `NaEdit`'s uppercase `Display`, losing the RNA case.
#[test]
fn rna_compact_mosaic_chimeric_preserves_lowercase() {
    for input in [
        "NM_004006.3:r.85=/u>c",
        "NM_004006.3:r.85=//u>c",
        "NR_002196.2:r.100=/a>g",
        "NR_002196.2:r.100=//a>g",
    ] {
        let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("must parse `{input}`: {e}"));
        assert_eq!(
            format!("{parsed}"),
            input,
            "RNA case round-trip for `{input}`"
        );
    }
}

/// parse → Display → parse cycle is stable across two iterations.
#[test]
fn compact_form_is_idempotent_under_repeated_display_parse() {
    let inputs = [
        "NC_012920.1:m.3243=/A>G",
        "NC_012920.1:m.3243=//A>G",
        "NC_012920.1:m.8470_8482=/del",
        "NC_012920.1:m.302_303=/dup",
    ];
    for input in inputs {
        let p1 = parse_hgvs(input).expect("parse 1");
        let s1 = format!("{p1}");
        let p2 = parse_hgvs(&s1).expect("parse 2");
        let s2 = format!("{p2}");
        assert_eq!(s1, s2, "idempotency 1->2 for `{input}`");
        let p3 = parse_hgvs(&s2).expect("parse 3");
        let s3 = format!("{p3}");
        assert_eq!(s2, s3, "idempotency 2->3 for `{input}`");
        assert_eq!(p1, p2, "equality 1<->2 for `{input}`");
        assert_eq!(p2, p3, "equality 2<->3 for `{input}`");
    }
}

/// Three input forms producing the same spec-compact Allele compare
/// equal and hash equal. This pins that `parse_hgvs` does not embed
/// input-shape state into the parsed variant.
#[test]
fn equivalent_compact_inputs_compare_and_hash_equal() {
    // (1) spec compact, (2) inherited-accession with type prefix on RHS,
    // (3) two fully qualified variants — all denote the same mosaic.
    let inputs = [
        "NC_012920.1:m.3243=/A>G",
        "NC_012920.1:m.3243=/m.3243A>G",
        "NC_012920.1:m.3243=/NC_012920.1:m.3243A>G",
    ];
    let parsed: Vec<HgvsVariant> = inputs
        .iter()
        .map(|i| parse_hgvs(i).expect("parse"))
        .collect();
    for (i, p) in parsed.iter().enumerate() {
        assert_eq!(
            parsed[0], *p,
            "Allele equality across input forms (input {i})"
        );
        assert_eq!(
            hash_of(&parsed[0]),
            hash_of(p),
            "Allele hash stability across input forms (input {i})"
        );
    }
    // All three Display in the spec compact form.
    let canonical = "NC_012920.1:m.3243=/A>G";
    for (i, p) in parsed.iter().enumerate() {
        assert_eq!(
            format!("{p}"),
            canonical,
            "canonical Display for input {i} (`{}`)",
            inputs[i]
        );
    }
}

/// Malformed compact forms must reject.
#[test]
fn malformed_compact_forms_reject() {
    let inputs = [
        "NC_012920.1:m.3243=/",    // empty RHS
        "NC_012920.1:m.3243=//",   // empty RHS, chimeric
        "NC_012920.1:m.3243=/X>Y", // RHS edit with non-IUPAC chars (X is not IUPAC)
    ];
    for input in inputs {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "malformed compact form `{input}` must reject. Got: {:?}",
            result
        );
    }
}

/// LHS that is not a position-bound identity edit must NOT trigger
/// the compact-form parser branch (and must NOT trigger the compact
/// Display branch either). The pre-#133 long-form behavior is preserved.
#[test]
fn non_identity_lhs_falls_through_to_long_form() {
    // Two fully qualified, neither side is `=`.
    let input = "NC_012920.1:m.3243A>G/NC_012920.1:m.3243A>T";
    let parsed = parse_hgvs(input).expect("long-form mosaic must parse");
    assert_eq!(
        format!("{parsed}"),
        input,
        "long form preserved when LHS is not position-identity"
    );
}

/// When LHS is position-identity but RHS interval differs, compact
/// Display does NOT apply — the Allele is still valid (parsed via the
/// fully-qualified branch), but emits in long form because the spec
/// compact form requires same position on both sides.
#[test]
fn intervals_differ_falls_through_to_long_form_on_display() {
    // LHS `=` at 3243, RHS substitution at 3244 — same accession but
    // different position. Parses via the fully-qualified branch.
    let input = "NC_012920.1:m.3243=/NC_012920.1:m.3244A>G";
    let parsed = parse_hgvs(input).expect("long-form must parse");
    // Long form is preserved on emission because intervals differ.
    assert_eq!(format!("{parsed}"), input);
}

/// Three input forms producing the same spec-compact chimeric Allele
/// compare equal and Display in spec compact form. Mirrors the mosaic
/// equivalence test above; pins that the chimeric path also accepts
/// the inherited-accession+type-prefix RHS (`<acc>:m.<pos>=//<type>.<edit>`)
/// in addition to spec-compact and fully-qualified inputs.
#[test]
fn equivalent_compact_chimeric_inputs_compare_and_hash_equal() {
    let inputs = [
        "NC_012920.1:m.3243=//A>G",
        "NC_012920.1:m.3243=//m.3243A>G",
        "NC_012920.1:m.3243=//NC_012920.1:m.3243A>G",
    ];
    let parsed: Vec<HgvsVariant> = inputs
        .iter()
        .map(|i| parse_hgvs(i).unwrap_or_else(|e| panic!("parse failed for `{i}`: {e}")))
        .collect();
    for (i, p) in parsed.iter().enumerate() {
        assert_eq!(
            parsed[0], *p,
            "Allele equality across chimeric input forms (input {i})"
        );
        assert_eq!(
            hash_of(&parsed[0]),
            hash_of(p),
            "Allele hash stability across chimeric input forms (input {i})"
        );
    }
    let canonical = "NC_012920.1:m.3243=//A>G";
    for (i, p) in parsed.iter().enumerate() {
        assert_eq!(
            format!("{p}"),
            canonical,
            "canonical Display for input {i} (`{}`)",
            inputs[i]
        );
    }
}

/// `var/=` shorthand: emits bare `=` on RHS post-#133 (no synthetic
/// `<acc>:<type>.1=` expansion).
#[test]
fn var_slash_eq_renders_bare_equals() {
    let inputs = [
        ("NC_012920.1:m.3243A>G/=", AllelePhase::Mosaic),
        ("NC_012920.1:m.3243A>G//=", AllelePhase::Chimeric),
        ("NM_000088.3:c.85T>C/=", AllelePhase::Mosaic),
        ("NC_000023.11:g.33344590A>G/=", AllelePhase::Mosaic),
    ];
    for (input, phase) in inputs {
        let parsed = parse_hgvs(input).expect("var/= must parse");
        let HgvsVariant::Allele(allele) = &parsed else {
            panic!("expected Allele for `{input}`, got {parsed:?}");
        };
        assert_eq!(allele.phase, phase, "phase for `{input}`");
        assert_eq!(format!("{parsed}"), input, "round-trip for `{input}`");
        // Idempotency.
        let s2 = format!("{}", parse_hgvs(&format!("{parsed}")).expect("reparse"));
        assert_eq!(s2, input, "idempotency for `{input}`");
    }
}
