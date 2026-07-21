//! MUST-level rejection of the multi-nucleotide substitution form
//! (`c.79GC>TT`, `c.79_80GC>TT`, `g.4GC>TG`) — issue #1079.
//!
//! # Spec
//!
//! `recommendations/DNA/delins.md:73`:
//!
//! > Can I describe a `GC` to `TG` variant as a di-nucleotide substitution
//! > (`g.4GC>TG`)? **No, this is not allowed.** By definition, a substitution
//! > changes **one** nucleotide into **one** other nucleotide. The change
//! > `TGT`GC`CA` to `TGT`TG`CA` should be described as `g.4_5delinsTG`.
//!
//! `recommendations/DNA/substitution.md:30` names both spellings explicitly:
//!
//! > based on the definition of a substitution, i.e. **one** nucleotide
//! > replaced by **one** other nucleotide, this change can not be described as
//! > a substitution like `c.79_80GC>TT` or `c.79GC>TT`.
//!
//! "not allowed" / "can not be described" is MUST-level under the spec's
//! RFC 2119 reading (`recommendations/style.md:9`).
//!
//! # The defect this file pins
//!
//! The preprocessor (lenient/silent) already rewrote both spellings to
//! `c.79_80delinsTT`, but the raw grammar path used by [`parse_hgvs`] accepted
//! the single-position spelling and **silently dropped a reference base**:
//! `c.79GC>TT` parsed as `c.79delinsTT`, which deletes one base and inserts
//! two — a different variant from the `c.79_80delinsTT` the CLI produced for
//! the same input. The two entry points must agree.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// The spec's own invalid spellings, paired with the canonical delins form.
const MULTIBASE_SUBS: &[(&str, &str)] = &[
    ("NM_004006.2:c.79GC>TT", "NM_004006.2:c.79_80delinsTT"),
    ("NM_004006.2:c.79_80GC>TT", "NM_004006.2:c.79_80delinsTT"),
    ("NC_000023.11:g.4GC>TG", "NC_000023.11:g.4_5delinsTG"),
    ("NG_012232.1:g.12GC>TG", "NG_012232.1:g.12_13delinsTG"),
];

// =====================================================================
// Default (strict) parse rejects the multi-base substitution form
// =====================================================================

#[test]
fn default_parse_rejects_multibase_substitution() {
    for (input, _) in MULTIBASE_SUBS {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "DNA/substitution.md:30 forbids {input:?}; got {:?}",
            result.map(|v| v.to_string())
        );
    }
}

#[test]
fn rejection_names_the_canonical_delins_form() {
    for (input, canonical) in MULTIBASE_SUBS {
        let msg = parse_hgvs(input).unwrap_err().to_string();
        let edit = canonical.split(':').nth(1).unwrap();
        assert!(
            msg.contains(edit),
            "the diagnostic for {input:?} should offer {edit:?}; got: {msg}"
        );
    }
}

// =====================================================================
// Lenient / silent canonicalise, and both spellings converge
// =====================================================================

#[test]
fn lenient_canonicalises_multibase_substitution() {
    for (input, canonical) in MULTIBASE_SUBS {
        let out = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .unwrap_or_else(|e| panic!("lenient must repair {input:?}: {e}"));
        assert_eq!(out.result.to_string(), *canonical, "input {input:?}");
    }
}

#[test]
fn silent_canonicalises_multibase_substitution() {
    for (input, canonical) in MULTIBASE_SUBS {
        let out = parse_hgvs_with_config(input, ErrorConfig::silent())
            .unwrap_or_else(|e| panic!("silent must repair {input:?}: {e}"));
        assert_eq!(out.result.to_string(), *canonical, "input {input:?}");
    }
}

#[test]
fn point_and_range_spellings_converge() {
    // The whole point of the fix: `c.79GC>TT` and `c.79_80GC>TT` describe the
    // same two-base change and must not yield different variants.
    let point = parse_hgvs_with_config("NM_004006.2:c.79GC>TT", ErrorConfig::lenient())
        .expect("lenient must repair the point spelling");
    let range = parse_hgvs_with_config("NM_004006.2:c.79_80GC>TT", ErrorConfig::lenient())
        .expect("lenient must repair the range spelling");
    assert_eq!(point.result.to_string(), range.result.to_string());
}

// =====================================================================
// The normalizer path (`parse_hgvs` + `Normalizer`) agrees with the CLI
// =====================================================================

#[test]
fn normalizer_path_never_sees_a_truncated_delins() {
    // `examples/generate_spec_fixture.rs` drives `parse_hgvs` then
    // `Normalizer<MockProvider>`. Before the fix this path produced
    // `c.79delinsTT` — one reference base instead of two. It must now either
    // reject (strict) or, on the repaired input, span both bases.
    let normalizer = Normalizer::new(MockProvider::new());
    for (input, canonical) in MULTIBASE_SUBS {
        assert!(
            parse_hgvs(input).is_err(),
            "{input:?} must not reach the normalizer as a truncated delins"
        );
        let repaired = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .unwrap_or_else(|e| panic!("lenient must repair {input:?}: {e}"));
        let normalized = normalizer
            .normalize_with_diagnostics(&repaired.result)
            .unwrap_or_else(|e| panic!("normalize {input:?}: {e}"));
        assert_eq!(normalized.result.to_string(), *canonical, "input {input:?}");
    }
}

// =====================================================================
// Canonical forms are untouched
// =====================================================================

#[test]
fn single_base_substitution_still_parses() {
    for input in [
        "NM_004006.2:c.79G>T",
        "NC_000023.10:g.33038255C>A",
        "NG_012232.1(NM_004006.2):c.93+1G>T",
        "NM_004006.3:r.79g>u",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "{input:?} is a canonical one-to-one substitution and must parse"
        );
    }
}

#[test]
fn canonical_delins_forms_still_parse() {
    for input in [
        "LRG_199t1:c.79_80delinsTT",
        "NM_004006.2:c.145_147delinsTGG",
        "NM_007294.3:c.2077delinsATA",
    ] {
        assert!(parse_hgvs(input).is_ok(), "{input:?} must parse");
    }
}

#[test]
fn single_base_reference_is_not_affected() {
    // Only a *multi-base reference* loses information at a point anchor. A
    // one-base reference with a longer insert (`c.79G>TT`) already spans
    // exactly the anchor base, so this rule does not touch it.
    assert!(parse_hgvs("NM_004006.2:c.79G>TT").is_ok());
}
