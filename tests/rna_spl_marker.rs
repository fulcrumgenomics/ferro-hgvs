//! E2: HGVS v21.0 RNA splicing markers — `r.spl`, `r.spl?`, `r.(spl)`, `r.(spl?)`.
//!
//! Spec: `assets/hgvs-nomenclature/docs/recommendations/RNA/splicing.md` and
//! `docs/recommendations/uncertain.md` (submodule sha `a32f970d`):
//!
//! - `r.spl` — RNA not analysed; splicing is **very likely** affected
//!   (canonical donor/acceptor: intron +1, +2, -2, -1, excl. GT↔GC).
//! - `r.spl?` — RNA not analysed; splicing **might** be affected
//!   (first/last exon nucleotide, intron +3..+6, new acceptor-adjacent AG).
//! - `r.(spl)` — predicted/uncertain wrapper around `spl`.
//! - `r.(spl?)` — predicted/uncertain wrapper around `spl?` (canonical
//!   predicted-splicing notation per HGVS v21.0).
//!
//! All four forms are now supported. This file pins:
//! - parse + Display round-trip (single and double pass) with an accession,
//! - `Normalizer::normalize` idempotency,
//! - hash + Eq distinctness across the four forms,
//! - non-equivalence with `r.?` (whole-RNA unknown),
//! - rejection of malformed inputs (`r.spl??`, `r.spline`, `r.SPL`, whitespace),
//! - round-trip of whole-entity edits inside allele brackets
//!   (`r.[spl?];[spl]`) as of #396 item 3.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

// ---------------------------------------------------------------------------
// `r.spl` — round-trip + normalize idempotency.
// ---------------------------------------------------------------------------

/// `r.spl` round-trips through `parse_hgvs` -> `Display` with an accession.
#[test]
fn rna_spl_round_trips_with_accession() {
    let input = "NM_004006.2:r.spl";
    let variant = parse_hgvs(input).expect("r.spl with accession should parse");
    assert_eq!(format!("{}", variant), input);
}

/// `r.spl` round-trips through a second parse — idempotent re-emission.
#[test]
fn rna_spl_double_round_trip_with_accession() {
    let input = "NM_004006.2:r.spl";
    let once = parse_hgvs(input).expect("first parse should succeed");
    let displayed = format!("{}", once);
    let twice = parse_hgvs(&displayed).expect("re-parse should succeed");
    assert_eq!(format!("{}", twice), input);
}

/// `Normalizer::normalize` is a no-op (idempotent) on `r.spl` since the marker
/// has no positional content to shift. Pin this so any future change to
/// normalization of whole-entity RNA edits is caught.
#[test]
fn rna_spl_normalize_idempotent() {
    let input = "NM_004006.2:r.spl";
    let variant = parse_hgvs(input).expect("r.spl should parse");
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize should succeed for whole-entity r.spl");
    assert_eq!(format!("{}", normalized), input);

    // Normalize a second time to confirm idempotency.
    let normalized_twice = normalizer
        .normalize(&normalized)
        .expect("second normalize should succeed");
    assert_eq!(format!("{}", normalized_twice), input);
}

// ---------------------------------------------------------------------------
// `r.spl?` — round-trip + normalize idempotency.
// ---------------------------------------------------------------------------

/// `r.spl?` ("splicing might be affected") parses and round-trips with an accession.
#[test]
fn rna_spl_question_mark_round_trips_with_accession() {
    let input = "NM_004006.2:r.spl?";
    let variant = parse_hgvs(input).expect("r.spl? with accession should parse");
    assert_eq!(format!("{}", variant), input);
}

/// `r.spl?` round-trips through a second parse — idempotent re-emission.
#[test]
fn rna_spl_question_mark_double_round_trip_with_accession() {
    let input = "NM_004006.2:r.spl?";
    let once = parse_hgvs(input).expect("first parse should succeed");
    let displayed = format!("{}", once);
    let twice = parse_hgvs(&displayed).expect("re-parse should succeed");
    assert_eq!(format!("{}", twice), input);
}

/// `Normalizer::normalize` is a no-op (idempotent) on `r.spl?`.
#[test]
fn rna_spl_question_mark_normalize_idempotent() {
    let input = "NM_004006.2:r.spl?";
    let variant = parse_hgvs(input).expect("r.spl? should parse");
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize should succeed for whole-entity r.spl?");
    assert_eq!(format!("{}", normalized), input);

    let normalized_twice = normalizer
        .normalize(&normalized)
        .expect("second normalize should succeed");
    assert_eq!(format!("{}", normalized_twice), input);
}

// ---------------------------------------------------------------------------
// `r.(spl)` and `r.(spl?)` — predicted wrappers.
// ---------------------------------------------------------------------------

/// `r.(spl)` (predicted wrapper around `spl`) round-trips with an accession.
#[test]
fn rna_predicted_spl_round_trips_with_accession() {
    let input = "NM_004006.2:r.(spl)";
    let variant = parse_hgvs(input).expect("r.(spl) with accession should parse");
    assert_eq!(format!("{}", variant), input);
}

/// `r.(spl?)` (predicted wrapper around `spl?`, the canonical predicted-splicing form per
/// HGVS v21.0) round-trips with an accession.
#[test]
fn rna_predicted_spl_question_mark_round_trips_with_accession() {
    let input = "NM_004006.2:r.(spl?)";
    let variant = parse_hgvs(input).expect("r.(spl?) with accession should parse");
    assert_eq!(format!("{}", variant), input);
}

/// `r.(spl)` round-trips through a second parse — idempotent re-emission.
#[test]
fn rna_predicted_spl_double_round_trip_with_accession() {
    let input = "NM_004006.2:r.(spl)";
    let once = parse_hgvs(input).expect("first parse should succeed");
    let displayed = format!("{}", once);
    let twice = parse_hgvs(&displayed).expect("re-parse should succeed");
    assert_eq!(format!("{}", twice), input);
}

/// `r.(spl?)` round-trips through a second parse — idempotent re-emission.
#[test]
fn rna_predicted_spl_question_mark_double_round_trip_with_accession() {
    let input = "NM_004006.2:r.(spl?)";
    let once = parse_hgvs(input).expect("first parse should succeed");
    let displayed = format!("{}", once);
    let twice = parse_hgvs(&displayed).expect("re-parse should succeed");
    assert_eq!(format!("{}", twice), input);
}

/// `Normalizer::normalize` is idempotent on `r.(spl)` and `r.(spl?)`.
#[test]
fn rna_predicted_spl_normalize_idempotent() {
    let normalizer = Normalizer::new(MockProvider::new());
    for input in ["NM_004006.2:r.(spl)", "NM_004006.2:r.(spl?)"] {
        let variant =
            parse_hgvs(input).unwrap_or_else(|e| panic!("{} should parse: {:?}", input, e));
        let normalized = normalizer
            .normalize(&variant)
            .unwrap_or_else(|e| panic!("normalize {} should succeed: {:?}", input, e));
        assert_eq!(
            format!("{}", normalized),
            input,
            "first normalize round-trip"
        );
        let normalized_twice = normalizer
            .normalize(&normalized)
            .unwrap_or_else(|e| panic!("second normalize {} should succeed: {:?}", input, e));
        assert_eq!(
            format!("{}", normalized_twice),
            input,
            "second normalize round-trip"
        );
    }
}

// ---------------------------------------------------------------------------
// Hash + Eq stability: the four forms are pairwise distinct.
// ---------------------------------------------------------------------------

/// `r.spl`, `r.spl?`, `r.(spl)`, `r.(spl?)` all hash and compare distinctly.
#[test]
fn rna_spl_four_forms_are_pairwise_distinct() {
    use std::collections::HashSet;

    let inputs = [
        "NM_004006.2:r.spl",
        "NM_004006.2:r.spl?",
        "NM_004006.2:r.(spl)",
        "NM_004006.2:r.(spl?)",
    ];
    let variants: Vec<_> = inputs
        .iter()
        .map(|s| parse_hgvs(s).unwrap_or_else(|e| panic!("{} should parse: {:?}", s, e)))
        .collect();

    // All four Display strings are pairwise distinct.
    let displays: HashSet<_> = variants.iter().map(|v| v.to_string()).collect();
    assert_eq!(
        displays.len(),
        4,
        "expected 4 distinct Display strings, got {:?}",
        displays
    );

    // All four hashed values are pairwise distinct.
    let by_hash: HashSet<_> = variants.iter().collect();
    assert_eq!(by_hash.len(), 4, "expected 4 distinct variant hashes");

    // Pairwise inequality.
    for (i, a) in variants.iter().enumerate() {
        for (j, b) in variants.iter().enumerate() {
            if i != j {
                assert_ne!(
                    a, b,
                    "variant {} should not equal variant {}",
                    inputs[i], inputs[j]
                );
            }
        }
    }
}

/// `r.spl?` (splicing-marker uncertain) is NOT equal to `r.?` (whole-RNA unknown).
/// They are different `NaEdit` variants — pin the non-equivalence.
#[test]
fn rna_spl_question_mark_not_equal_to_rna_unknown() {
    let spl_q = parse_hgvs("NM_004006.2:r.spl?").expect("r.spl? should parse");
    let r_q = parse_hgvs("NM_004006.2:r.?").expect("r.? should parse");
    assert_ne!(spl_q, r_q, "r.spl? and r.? must remain distinct");
    assert_ne!(spl_q.to_string(), r_q.to_string());
}

// ---------------------------------------------------------------------------
// Allele compounds — supported as of #396 item 3.
//
// `parse_rna_bracket_member` (the per-bracket dispatcher used by both
// `parse_rna_allele_shorthand` cis/unknown paths and
// `parse_rna_trans_allele_shorthand`) now recognises the whole-entity RNA
// edits `=` / `?` / `0` / `spl` / `spl?` / `(spl)` / `(spl?)` in addition to
// the previously-supported position-based edits and predicted-wrapper form.
// Whole-entity members attach a dummy interval at position 1 (mirroring
// `parse_rna_variant`'s top-level handling).
// ---------------------------------------------------------------------------

/// `r.[spl?];[spl]` must round-trip after #396 item 3. Replaces the older
/// "currently rejected" pin.
#[test]
fn rna_spl_allele_compound_roundtrip() {
    let variant = parse_hgvs("NM_004006.2:r.[spl?];[spl]")
        .expect("r.[spl?];[spl] must parse after #396 item 3");
    // Exact-string equality so a regression that drops the `?` (yielding
    // `[spl];[spl]`) or that loses the trans-form split (`[spl?];[spl]` →
    // `[spl?;spl]`) is caught — a `contains("spl")` check would silently
    // match `spl?` and miss either regression.
    assert_eq!(variant.to_string(), "NM_004006.2:r.[spl?];[spl]");
}

// ---------------------------------------------------------------------------
// Negative cases — these MUST fail to parse.
// ---------------------------------------------------------------------------

/// `r.spl??` (double `?`) is not a valid HGVS form.
#[test]
fn rna_spl_double_question_mark_rejected() {
    let result = parse_hgvs("NM_004006.2:r.spl??");
    assert!(
        result.is_err(),
        "r.spl?? must be rejected; got Ok({:?})",
        result.ok().map(|v| v.to_string()),
    );
}

/// `r.spline` (extra letters after `spl`) — the parser must not accept the prefix and
/// silently drop the trailing characters.
#[test]
fn rna_spl_trailing_letters_rejected() {
    let result = parse_hgvs("NM_004006.2:r.spline");
    assert!(
        result.is_err(),
        "r.spline must be rejected; got Ok({:?})",
        result.ok().map(|v| v.to_string()),
    );
}

/// Uppercase `SPL` violates the HGVS lowercase-RNA rule.
#[test]
fn rna_spl_uppercase_rejected() {
    let result = parse_hgvs("NM_004006.2:r.SPL");
    assert!(
        result.is_err(),
        "r.SPL must be rejected (RNA must be lowercase); got Ok({:?})",
        result.ok().map(|v| v.to_string()),
    );
}

/// Whitespace between `spl` and `?` is not valid.
#[test]
fn rna_spl_whitespace_before_question_rejected() {
    let result = parse_hgvs("NM_004006.2:r.spl ?");
    assert!(
        result.is_err(),
        "r.spl <space>? must be rejected; got Ok({:?})",
        result.ok().map(|v| v.to_string()),
    );
}

/// Whitespace inside the predicted parens is not valid.
#[test]
fn rna_predicted_spl_inner_whitespace_rejected() {
    let result = parse_hgvs("NM_004006.2:r.( spl?)");
    assert!(
        result.is_err(),
        "r.( spl?) must be rejected; got Ok({:?})",
        result.ok().map(|v| v.to_string()),
    );
}
