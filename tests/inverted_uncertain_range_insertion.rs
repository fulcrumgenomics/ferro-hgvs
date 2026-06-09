//! Inverted insertion of an uncertain-boundary range — the spec-sanctioned way
//! (HGVS DNA/complex.md) to describe an orientation-reversed duplication. HGVS
//! writes this as `ins(a_b)_(c_d)inv`, NOT as the invalid `dupinv`
//! (DNA/inversion.md:19 marks `dupinv` `class="invalid"`). Issue #546.

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

/// The canonical orientation-reversed intra-chromosomal duplication from
/// `DNA/complex.md` (ISCN `dup(8)(q24.22q24.21)`, breakpoint not sequenced).
#[test]
fn inverted_uncertain_range_insertion_round_trips() {
    let s = "NC_000008.11:g.(131500001_136400000)ins(127300001_131500000)_(131500001_136400000)inv";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
}

/// Guard: the bare (non-`inv`) compound paren range keeps its established
/// canonical bracket form (`[a_b;c_d]`) — the `inv` support must not disturb it.
#[test]
fn bare_compound_paren_insertion_keeps_bracket_canonical_form() {
    let s = "NG_007406.1:g.?_?ins(23632682_23625413)_(23625324_23619334)";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(
        format!("{v}"),
        "NG_007406.1:g.?ins[23632682_23625413;23625324_23619334]",
        "bare compound paren insertion canonicalizes to brackets"
    );
}

/// The normalizer must pass the `UncertainRangeInv` form through unchanged:
/// `expand_inserted_sequence` returns `Ok(None)`, and `normalize_na_edit`'s
/// `Insertion` arm bails early when `sequence.bases()` returns `None`.
#[test]
fn inverted_uncertain_range_insertion_normalizes_idempotently() {
    let s = "NC_000008.11:g.(131500001_136400000)ins(127300001_131500000)_(131500001_136400000)inv";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    let provider = MockProvider::new();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::default());
    let result = normalizer.normalize(&v).expect("normalize");
    assert_eq!(
        format!("{result}"),
        s,
        "normalize must pass through unchanged for `{s}`"
    );
}

/// `delins` with an `UncertainRangeInv` inserted payload — both `ins` and
/// `delins` share `parse_parenthesized_count` as the entry point, so the new
/// variant parses and round-trips for `delins` too.
#[test]
fn delins_uncertain_range_inv_round_trips() {
    let s =
        "NC_000008.11:g.(131500001_136400000)delins(127300001_131500000)_(131500001_136400000)inv";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
}
