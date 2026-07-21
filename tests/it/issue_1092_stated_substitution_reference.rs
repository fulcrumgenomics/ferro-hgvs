//! The stated reference run of a `NN>MM` substitution survives parsing — issue #1092.
//!
//! # The defect
//!
//! `parse_multibase_substitution` folded `GC>TT` into
//! `NaEdit::Delins { deleted: None, deleted_length: None, inserted: "TT" }`,
//! discarding the stated reference run `GC` at parse time. Three consequences:
//!
//! 1. `docs/syntax.yaml:209` (`dna.sub`) has **no** production omitting
//!    `reference_sequence` — the bases are a required element of that input
//!    form, so ferro dropped a mandatory field.
//! 2. ferro already validates the neighbouring spelling: `validate_reference`
//!    checks a stated `delGCinsTT` deleted run against the reference (#486).
//!    So the *deprecated* spelling was verified while the *forbidden* one,
//!    asserting the identical thing, was not.
//! 3. The MUST-level rejection added for #1079 was **source-keyed** — it
//!    re-read the input string, so it could not fire for an AST composed
//!    programmatically, by projection, or by round-trip.
//!
//! # The fix under test
//!
//! `NaEdit::Delins` gains a `substitution_reference: Option<Sequence>`
//! provenance field, set only by the `>` spelling. `Display` ignores it (the
//! canonical short form `delinsTT` is preserved — rendering `delGCinsTT`
//! would swap one non-conformant output for another, see `DNA/delins.md:29`),
//! the rejection rule keys off it instead of the source string, and
//! `validate_reference` feeds it to the same check `deleted` already gets.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::edit::{InsertedSequence, NaEdit, Sequence};
use ferro_hgvs::hgvs::interval::CdsInterval;
use ferro_hgvs::hgvs::location::CdsPos;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::hgvs::parser::variant::validate_no_multibase_substitution;
use ferro_hgvs::hgvs::variant::{Accession, CdsVariant, LocEdit};
use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::{parse_hgvs, HgvsVariant, Normalizer};
use std::str::FromStr;

use crate::common::synthetic::SyntheticBuilder;

/// The spec's own invalid spellings, paired with the canonical delins form.
const MULTIBASE_SUBS: &[(&str, &str)] = &[
    ("NM_004006.2:c.79GC>TT", "NM_004006.2:c.79_80delinsTT"),
    ("NM_004006.2:c.79_80GC>TT", "NM_004006.2:c.79_80delinsTT"),
    ("NC_000023.11:g.4GC>TG", "NC_000023.11:g.4_5delinsTG"),
];

/// `NM_TEST.1:c.79_80GC>TT` as an AST, built by hand — never parsed.
fn programmatic_multibase_sub() -> HgvsVariant {
    HgvsVariant::Cds(CdsVariant {
        accession: Accession::new("NM", "TEST", Some(1)),
        gene_symbol: None,
        loc_edit: LocEdit::new(
            CdsInterval::new(CdsPos::new(79), CdsPos::new(80)),
            NaEdit::Delins {
                sequence: InsertedSequence::Literal(Sequence::from_str("TT").unwrap()),
                deleted: None,
                deleted_length: None,
                substitution_reference: Some(Sequence::from_str("GC").unwrap()),
            },
        ),
    })
}

// =====================================================================
// 1. Display is unchanged — the stated run is preserved, never rendered
// =====================================================================

#[test]
fn display_never_renders_the_deprecated_del_ins_spelling() {
    // The whole reason for a separate provenance field: reusing `deleted`
    // would render `c.79_80delGCinsTT`, the spelling `DNA/delins.md:29-30`
    // tells you not to use.
    let rendered = programmatic_multibase_sub().to_string();
    assert_eq!(rendered, "NM_TEST.1:c.79_80delinsTT");
    assert!(
        !rendered.contains("delGCins"),
        "the stated run must be preserved in the AST but never rendered; got {rendered}"
    );
}

#[test]
fn lenient_and_silent_display_is_unchanged() {
    for (input, canonical) in MULTIBASE_SUBS {
        for config in [ErrorConfig::lenient(), ErrorConfig::silent()] {
            let out = parse_hgvs_with_config(input, config)
                .unwrap_or_else(|e| panic!("must repair {input:?}: {e}"));
            let rendered = out.result.to_string();
            assert_eq!(rendered, *canonical, "input {input:?}");
            assert!(
                !rendered.contains("delGCinsTT") && !rendered.contains("delGCinsTG"),
                "input {input:?} rendered the deprecated explicit form: {rendered}"
            );
        }
    }
}

// =====================================================================
// 2. Strict mode still rejects, naming the canonical repair (#1079)
// =====================================================================

#[test]
fn strict_still_rejects_and_names_the_canonical_repair() {
    for (input, canonical) in MULTIBASE_SUBS {
        let err = parse_hgvs(input)
            .expect_err(&format!("DNA/substitution.md:30 forbids {input:?}"))
            .to_string();
        let edit = canonical.split(':').nth(1).unwrap();
        assert!(
            err.contains(edit),
            "the diagnostic for {input:?} should offer {edit:?}; got: {err}"
        );
    }
}

// =====================================================================
// 3. The rule is AST-keyed: it fires on a variant that was never parsed
// =====================================================================

#[test]
fn programmatically_built_multibase_substitution_is_rejected() {
    let variant = programmatic_multibase_sub();
    let err = validate_no_multibase_substitution(&variant)
        .expect_err("an AST carrying a multi-base stated reference run must be rejected")
        .to_string();
    assert!(
        err.contains("c.79_80delinsTT"),
        "the diagnostic must name the canonical repair; got: {err}"
    );
}

#[test]
fn ordinary_delins_ast_is_accepted() {
    // Same edit without the `>` provenance — a plain `c.79_80delinsTT` — is
    // conformant and must pass the same AST-keyed rule.
    let variant = HgvsVariant::Cds(CdsVariant {
        accession: Accession::new("NM", "TEST", Some(1)),
        gene_symbol: None,
        loc_edit: LocEdit::new(
            CdsInterval::new(CdsPos::new(79), CdsPos::new(80)),
            NaEdit::Delins {
                sequence: InsertedSequence::Literal(Sequence::from_str("TT").unwrap()),
                deleted: None,
                deleted_length: None,
                substitution_reference: None,
            },
        ),
    });
    assert!(validate_no_multibase_substitution(&variant).is_ok());
}

// =====================================================================
// 4. The legal single-base substitution is untouched
// =====================================================================

#[test]
fn single_base_substitution_is_unaffected() {
    for input in [
        "NM_004006.2:c.79G>T",
        "NC_000023.10:g.33038255C>A",
        "NM_004006.3:r.79g>u",
    ] {
        let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input:?} must parse: {e}"));
        assert_eq!(variant.to_string(), *input);
        assert!(validate_no_multibase_substitution(&variant).is_ok());
    }
}

// =====================================================================
// 5. Lenient/silent now report a false stated reference
// =====================================================================

/// A synthetic genomic contig whose bases at 257_258 are `AC`, so a stated
/// `GC` run is a *false* claim about the reference.
fn synthetic_genome() -> ferro_hgvs::MockProvider {
    // Core starts at 1-based position 257 (see `common::synthetic`).
    SyntheticBuilder::genomic("ACGTACGTACGTACGT").build()
}

#[test]
fn lenient_reports_a_false_stated_reference_run() {
    let normalizer = Normalizer::new(synthetic_genome());
    // Reference at 257_258 is "AC"; the input claims "GC".
    let parsed = parse_hgvs_with_config("NC_TEST.1:g.257_258GC>TT", ErrorConfig::lenient())
        .expect("lenient must repair the multi-base substitution");
    let out = normalizer
        .normalize_with_diagnostics(&parsed.result)
        .expect("normalize");
    let mismatch = out.warnings.iter().find_map(|w| match w {
        NormalizationWarning::RefSeqMismatch {
            stated_ref,
            actual_ref,
            ..
        } => Some((stated_ref.clone(), actual_ref.clone())),
        _ => None,
    });
    let (stated, actual) = mismatch.unwrap_or_else(|| {
        panic!(
            "a false stated reference run must be diagnosed, not silently dropped; \
             warnings were {:?}",
            out.warnings
        )
    });
    assert_eq!(stated, "GC");
    assert_eq!(actual, "AC");
}

#[test]
fn lenient_stays_quiet_when_the_stated_run_is_true() {
    let normalizer = Normalizer::new(synthetic_genome());
    // Reference at 257_258 really is "AC".
    let parsed = parse_hgvs_with_config("NC_TEST.1:g.257_258AC>TT", ErrorConfig::lenient())
        .expect("lenient must repair the multi-base substitution");
    let out = normalizer
        .normalize_with_diagnostics(&parsed.result)
        .expect("normalize");
    assert!(
        !out.warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. })),
        "a true stated reference run must not warn; got {:?}",
        out.warnings
    );
}
