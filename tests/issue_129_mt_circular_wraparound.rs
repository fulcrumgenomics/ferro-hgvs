//! Issue #129 — Mitochondrial circular reference handling (Path 1:
//! strict spec compliance).
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/129>.
//!
//! HGVS spec analysis:
//! - `deletion.md:17` — for `o.`/`m.` circular references, positions
//!   MAY be listed 3'→5' "when the deletion includes both the last and
//!   first nucleotides of the reference sequence."
//! - `docs/consultation/SVD-WG006.md` (accepted into HGVS v19.01,
//!   founding proposal for circular DNA) — line 15 grants the
//!   exception for any "rearrangement" that includes both ends, with
//!   explicit examples for both `del`
//!   (`NC_012920.1:m.16563_13del`) and `dup`
//!   (`J01749.1:o.4344_197dup`). Line 33 names
//!   "deletions/duplications" as the spec-authorised scope.
//! - `delins` inherits the exception via the `del + ins` composition.
//! - `insertion.md` and `inversion.md` are spec-silent on wraparound
//!   forms; no examples and no rationale. The general "5'→3'" rule
//!   applies.
//!
//! Comparison: both mutalyzer and biocommons hgvs reject ALL reversed
//! ranges including the spec-authorized `del`/`dup` exception (see
//! `mutalyzer/description.py:889`, `hgvs/location.py:437`). ferro
//! goes more spec-compliant by accepting the explicitly-authorized
//! `del`/`dup`/`delins` wraparound forms while rejecting the
//! unauthorized `ins`/`inv` reversed forms — closing a silent-accept
//! gap pinned by `tests/mito_circular_audit.rs`.
//!
//! What this PR delivers:
//! - Parser accepts `m.<high>_<low>del`, `m.<high>_<low>delins`, and
//!   `m.<high>_<low>dup` (and the same on `o.`).
//! - Parser **rejects** `m.<high>_<low>ins`/`inv` (were silently
//!   accepted; mis-computed indel length).
//! - 3'-rule wraparound shift remains **disabled** — matches
//!   mutalyzer + biocommons + the strict spec reading. A homopolymer
//!   straddling the origin won't shift across it.

use ferro_hgvs::{parse_hgvs, FerroError, MockProvider, Normalizer};

const MT: &str = "NC_012920.1";

// ============================================================================
// Wraparound `del` and `delins` — spec-authorized, must parse.
// ============================================================================

/// Helper: parse + round-trip-display. Tight `assert_eq!` catches any
/// stray leading/trailing characters in the Display path that a
/// `contains` check would miss.
fn assert_roundtrip(input: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse({input:?}) failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, input, "round-trip mismatch for {input:?}: got {out:?}");
}

#[test]
fn mt_wraparound_del_origin_pair_parses() {
    // `m.16569_1del` — the 2-nt deletion of "both the last and first
    // nucleotides of the reference sequence". Spec-authorized.
    assert_roundtrip(&format!("{MT}:m.16569_1del"));
}

#[test]
fn mt_wraparound_del_longer_parses() {
    // A longer wraparound del. Per the spec exception "includes both
    // the last and first nucleotides" means the deletion spans the
    // origin; any `<high>_<low>` form on a circular ref qualifies.
    assert_roundtrip(&format!("{MT}:m.16500_100del"));
}

#[test]
fn mt_wraparound_delins_parses() {
    // `delins` inherits the deletion-side exception via its `del + ins`
    // composition.
    assert_roundtrip(&format!("{MT}:m.16569_1delinsAA"));
}

// ============================================================================
// Wraparound `dup` — spec-authorized per SVD-WG006 (line 15 + line 23
// example `J01749.1:o.4344_197dup`). Must parse.
// ============================================================================

#[test]
fn mt_wraparound_dup_parses() {
    assert_roundtrip(&format!("{MT}:m.16569_1dup"));
}

#[test]
fn circular_wraparound_dup_parses_pbr322() {
    // Verbatim from SVD-WG006:23.
    assert_roundtrip("J01749.1:o.4344_197dup");
}

// ============================================================================
// Wraparound `ins` / `inv` — spec-silent on circular references (no
// examples or rationale in SVD-WG006). ferro previously silently
// accepted these with mis-computed length; this PR rejects them per
// the strict reading of the spec's "5'→3'" general rule.
// ============================================================================

#[test]
fn mt_wraparound_ins_is_rejected() {
    // No spec exception for `ins` on circular refs. (Note: ins is the
    // standard `<X>_<X+1>ins<seq>` adjacent-position form; `<high>_<low>`
    // is a separate spec-unauthorized shape.)
    //
    // Pin the rejection class to `FerroError::Parse` — not just `is_err()` —
    // so the test stays meaningful if a downstream stage starts emitting a
    // different error variant. `check_circular_reversed_range` lifts the
    // rejection via nom's `ErrorKind::Verify`, which `parse_variant` wraps
    // as `FerroError::Parse`.
    let err = parse_hgvs(&format!("{MT}:m.16569_1insA"))
        .expect_err("wraparound `ins` has no spec exception; must be rejected");
    assert!(
        matches!(err, FerroError::Parse { .. }),
        "wraparound `ins` must reject at parse stage; got {err:?}",
    );
}

#[test]
fn mt_wraparound_inv_is_rejected() {
    // No spec exception for `inv` on circular refs.
    let err = parse_hgvs(&format!("{MT}:m.16569_1inv"))
        .expect_err("wraparound `inv` has no spec exception; must be rejected");
    assert!(
        matches!(err, FerroError::Parse { .. }),
        "wraparound `inv` must reject at parse stage; got {err:?}",
    );
}

// ============================================================================
// `o.` circular references — same rules per `deletion.md:17` (which
// names both "o" and "m" prefixes).
// ============================================================================

#[test]
fn circular_wraparound_del_parses() {
    assert_roundtrip("NC_000000.0:o.500_100del");
}

#[test]
fn circular_wraparound_inv_is_rejected() {
    // `inv` is spec-silent on circular refs; reject.
    let err = parse_hgvs("NC_000000.0:o.500_100inv")
        .expect_err("o. wraparound `inv` has no spec exception; must be rejected");
    assert!(
        matches!(err, FerroError::Parse { .. }),
        "o. wraparound `inv` must reject at parse stage; got {err:?}",
    );
}

// ============================================================================
// Non-wraparound `m.` variants (start < end) must continue to parse
// for all edit kinds — this is the common case and shouldn't change.
// ============================================================================

#[test]
fn mt_non_wraparound_dup_still_parses() {
    assert_roundtrip(&format!("{MT}:m.100_103dup"));
}

#[test]
fn mt_non_wraparound_inv_still_parses() {
    assert_roundtrip(&format!("{MT}:m.100_103inv"));
}

// ============================================================================
// `g.` axis must NOT get the circular exception — genomic references
// are linear and the spec authorizes the wraparound form only for
// `o.`/`m.` prefixes. This test pins that the `g.` rejection stays in
// place.
// ============================================================================

// ============================================================================
// Bracketed compound alleles — wraparound must propagate through the
// shared bracket-member parser too (CodeRabbit follow-up).
// ============================================================================

#[test]
fn mt_compound_cis_allele_with_wraparound_del_parses() {
    // SVD-WG006 wraparound del as a member of a cis-allele.
    assert_roundtrip(&format!("{MT}:m.[16569_1del;100A>G]"));
}

#[test]
fn circular_compound_cis_allele_with_wraparound_dup_parses() {
    assert_roundtrip("J01749.1:o.[4344_197dup;500A>G]");
}

#[test]
fn mt_compound_cis_allele_with_wraparound_ins_is_rejected() {
    // `ins` is spec-silent on circular refs; reject inside brackets
    // just like outside.
    let err = parse_hgvs(&format!("{MT}:m.[16569_1insA;100A>G]"))
        .expect_err("compound-bracket wraparound ins must reject");
    assert!(
        matches!(err, FerroError::Parse { .. }),
        "compound-bracket wraparound ins must reject at parse stage; got {err:?}",
    );
}

#[test]
fn mt_trans_allele_with_wraparound_del_parses() {
    // `m.[edit1];[edit2]` shorthand — both members must accept the
    // circular exception per the same code path.
    assert_roundtrip(&format!("{MT}:m.[16569_1del];[100A>G]"));
}

#[test]
fn genome_wraparound_del_is_still_rejected() {
    // Linear g. has no circular exception. Rejected by the standard
    // (non-circular) `parse_genome_interval` inverted-range check rather
    // than `check_circular_reversed_range`; both surface as
    // `FerroError::Parse`.
    let err = parse_hgvs("NC_000001.11:g.16569_1del")
        .expect_err("linear g. has no circular exception; reversed range must stay rejected");
    assert!(
        matches!(err, FerroError::Parse { .. }),
        "linear g. reversed range must reject at parse stage; got {err:?}",
    );
}

// ============================================================================
// Normalization safety on reversed `m.`/`o.` ranges (CodeRabbit follow-up).
//
// `normalize_mt`'s reversed-range branch routes to `mt_fallback`, which
// historically still invoked `apply_canonical_split`. The canonical-split
// path computes `expected_span = hgvs_end - hgvs_start + 1` as `u64`, so a
// reversed range (`end < start`) would underflow into a giant span and
// trigger a `CanonicalSplitSkipped` warning (or worse, depending on the
// reference provider's tolerance for `start > end`). Today's shipped
// providers reject `start > end` at the `get_sequence` boundary, masking
// the issue, but the latent dependency on every provider doing that guard
// correctly is fragile. The fix in `normalize/mod.rs::normalize_mt`
// short-circuits the reversed-range fallback before `apply_canonical_split`
// runs. This test pins that contract.
// ============================================================================

#[test]
fn mt_wraparound_delins_normalize_does_not_emit_canonical_split_skipped() {
    let variant =
        parse_hgvs(&format!("{MT}:m.16569_1delinsT")).expect("wraparound delins must parse");
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize must succeed on no-reference fallback");
    assert_eq!(
        format!("{}", normalized.result),
        format!("{MT}:m.16569_1delinsT"),
        "wraparound delins must round-trip unchanged through the no-reference fallback",
    );
    assert!(
        !normalized
            .warnings
            .iter()
            .any(|w| w.code() == "CANONICAL_SPLIT_SKIPPED"),
        "wraparound mt delins must skip `apply_canonical_split` entirely; \
         got warnings {:?}",
        normalized
            .warnings
            .iter()
            .map(|w| w.code())
            .collect::<Vec<_>>(),
    );
}

#[test]
fn mt_wraparound_del_normalize_does_not_emit_canonical_split_skipped() {
    let variant =
        parse_hgvs(&format!("{MT}:m.16563_13del")).expect("wraparound del must parse (SVD-WG006)");
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize must succeed on no-reference fallback");
    assert!(
        !normalized
            .warnings
            .iter()
            .any(|w| w.code() == "CANONICAL_SPLIT_SKIPPED"),
        "wraparound mt del must skip `apply_canonical_split` entirely; \
         got warnings {:?}",
        normalized
            .warnings
            .iter()
            .map(|w| w.code())
            .collect::<Vec<_>>(),
    );
}
