//! Issue #1097 — `RefSeqMismatch → Reject` must fire for the asserted reference
//! bases of a multi-base range substitution (`g.pos_pos REF>ALT`), exactly as it
//! already does for the single-base (`g.posA>C`) and per-base allele
//! (`g.[posA>C;posA>G]`) forms.
//!
//! # The defect (reported against v0.9.1)
//!
//! With reference validation explicitly turned on
//! (`ErrorConfig::silent().with_override(RefSeqMismatch, Reject)`), ferro
//! silently accepted a *wrong* asserted reference base in the multi-base range
//! form. The preprocessor rewrites `pos_pos REF>ALT` to `delins` before parsing,
//! and the asserted `REF` run was dropped — so the reject the caller opted into
//! never saw those bases to compare. The identical mismatch *was* rejected in the
//! single-base and per-base allele forms, so validation was inconsistent purely
//! by input syntax.
//!
//! #1092 introduced the machinery that keeps the asserted run alive on the
//! rewritten `Delins` (`substitution_reference`) and routes it through
//! `validate_reference`. This test pins the specific guarantee #1097 asks for:
//! under `RefSeqMismatch → Reject`, all four syntactic forms of the *same*
//! wrong assertion must reject — and a correct-ref range substitution must still
//! normalize cleanly.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorOverride, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::{FerroError, JsonProvider, NormalizeConfig, Normalizer};
use std::io::Write;

/// Genome-capable provider: one `contig` of `len` cyclic ACGT bytes with
/// `payload` written 1-based at `pos1`. (Mirrors tests/it/issue_1052_substitution_refseq.rs.)
fn g_provider(contig: &str, len: usize, pos1: usize, payload: &str) -> JsonProvider {
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(len).collect();
    for (i, b) in payload.bytes().enumerate() {
        bases[pos1 - 1 + i] = b;
    }
    let seq = String::from_utf8(bases).unwrap();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// The reproduction reference from the issue: 1-based positions 22..24 are `CGC`.
fn ref_cgc_at_22() -> JsonProvider {
    g_provider("c", 300, 22, "CGC")
}

/// `silent()` base — so the deprecated multi-base `>` form is rewritten to
/// `delins` rather than rejected as `OldSubstitutionSyntax` (W3003) — with
/// `RefSeqMismatch` escalated to `Reject`. This is the exact opt-in the issue
/// describes: reference validation on, everything else quiet.
fn silent_but_reject_ref_mismatch() -> ErrorConfig {
    ErrorConfig::silent().with_override(ErrorType::RefSeqMismatch, ErrorOverride::Reject)
}

/// Parse under the opt-in config, then normalize under it, returning the
/// normalize result. Mirrors the single-`Normalizer(error_config=...)` call the
/// issue's Python reproduction makes.
fn normalize_under_reject(input: &str) -> Result<ferro_hgvs::HgvsVariant, FerroError> {
    let cfg = silent_but_reject_ref_mismatch();
    let parsed = parse_hgvs_with_config(input, cfg.clone())
        .unwrap_or_else(|e| panic!("{input:?} must parse/preprocess: {e}"));
    Normalizer::with_config(
        ref_cgc_at_22(),
        NormalizeConfig::silent().with_error_config(cfg),
    )
    .normalize(&parsed.result)
}

/// Assert `err` is a `ReferenceMismatch` naming the asserted-vs-real bases and
/// span the case is testing. The `ReferenceMismatch` fields carry `expected` =
/// the *asserted* reference run the input stated, `found` = the *actual*
/// reference bases; pinning them keeps the test from passing on an incidental
/// mismatch from some other span (the whole point of #1097 is that the asserted
/// run — not the syntax — drives the reject).
fn assert_ref_mismatch(err: FerroError, location: &str, asserted: &str, real: &str) {
    match err {
        FerroError::ReferenceMismatch {
            location: loc,
            expected,
            found,
        } => {
            assert_eq!(loc, location, "mismatch location");
            assert_eq!(expected, asserted, "asserted reference run");
            assert_eq!(found, real, "actual reference bases");
        }
        other => panic!("expected ReferenceMismatch, got {other:?}"),
    }
}

// =====================================================================
// The wrong-assertion matrix — all four forms must reject identically.
// Real bases: C@22, G@23, C@24.
// =====================================================================

#[test]
fn single_base_sub_wrong_ref_rejects() {
    // Baseline that already worked: asserts A@22, real is C.
    let err = normalize_under_reject("c:g.22A>C")
        .expect_err("single-base wrong-ref sub must reject under RefSeqMismatch=Reject");
    assert_ref_mismatch(err, "22-22", "A", "C");
}

#[test]
fn per_base_allele_subs_wrong_ref_rejects() {
    // Baseline that already worked: asserts A@22 (real C) and A@24 (real C).
    // The first mismatch (position 22) is the one surfaced.
    let err = normalize_under_reject("c:g.[22A>C;24A>G]")
        .expect_err("per-base allele wrong-ref subs must reject under RefSeqMismatch=Reject");
    assert_ref_mismatch(err, "22-22", "A", "C");
}

#[test]
fn three_base_range_sub_wrong_ref_rejects() {
    // The bug: asserts AAA at 22_24, real bases are CGC.
    let err = normalize_under_reject("c:g.22_24AAA>CAG")
        .expect_err("3-base range wrong-ref sub must reject under RefSeqMismatch=Reject");
    assert_ref_mismatch(err, "22-24", "AAA", "CGC");
}

#[test]
fn two_base_range_sub_wrong_ref_rejects() {
    // The bug: asserts AA at 22_23, real bases are CG.
    let err = normalize_under_reject("c:g.22_23AA>CT")
        .expect_err("2-base range wrong-ref sub must reject under RefSeqMismatch=Reject");
    assert_ref_mismatch(err, "22-23", "AA", "CG");
}

// =====================================================================
// A correct-ref range substitution must still normalize cleanly — the
// reject must key off the asserted bases, not the syntax.
// =====================================================================

#[test]
fn correct_ref_range_sub_normalizes_cleanly() {
    // Asserts the TRUE bases CGC at 22_24 → must not reject.
    let out = normalize_under_reject("c:g.22_24CGC>TAG")
        .expect("a correct-ref range substitution must normalize, not reject");
    // The canonical delins form is emitted (the `>` provenance never renders);
    // CGC>TAG shares no prefix/suffix, so the whole span is a delins.
    assert_eq!(out.to_string(), "c:g.22_24delinsTAG");
}
