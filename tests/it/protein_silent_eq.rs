//! Round-trip pinning for protein silent / synonymous notation (`p.=` family) — issue #81 D5.
//!
//! HGVS spec context (assets/hgvs-nomenclature/docs/recommendations/protein/):
//!
//! - `p.Leu54=`   — observed, position-specific synonymous change.
//! - `p.(Leu54=)` — predicted (uncertain) position-specific synonymous.
//! - `p.=`        — observed whole-protein identity (no change).
//! - `p.(=)`      — predicted whole-protein identity.
//!
//! Distinct from `p.0` / `p.0?` (no protein, D6) and `p.?` / `p.(?)` / `p.Met1?` /
//! `p.(Met1?)` (unknown effect, D7). Each row of Section D in #81 is independently
//! enforceable; this file pins every D5 invariant the existing parse / Display /
//! normalize stack provides today, so a future refactor that collapses any pair
//! (predicted↔certain, whole-protein↔position-specific, silent↔no-protein,
//! silent↔unknown) fails loudly here.
//!
//! Spec corpus cross-link: `tests/fixtures/grammar/hgvs_spec_normalization.json`
//! row 7457 covers `NP_003997.1:p.(Lys874=)` with `status: "preserved"`. We do
//! not duplicate that here.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use ferro_hgvs::hgvs::edit::ProteinEdit;
use ferro_hgvs::hgvs::interval::ProtInterval;
use ferro_hgvs::hgvs::location::{AminoAcid, ProtPos};
use ferro_hgvs::hgvs::uncertainty::Mu;
use ferro_hgvs::hgvs::variant::ProteinVariant;
use ferro_hgvs::{parse_hgvs, parse_hgvs_fast, HgvsVariant, MockProvider, Normalizer};

/// One canonical silent-= input plus the structural shape we expect after parse.
///
/// `outer_uncertain`: true iff the parser wraps the inner `ProteinEdit` in
/// `Mu::Uncertain` (currently only `p.(Leu54=)` — the position-specific
/// predicted form — takes that path; the whole-protein `p.(=)` keeps the
/// `Mu::Certain` outer with `predicted: true` on the inner edit).
struct SilentCase {
    input: &'static str,
    /// True for whole-protein (`p.=`, `p.(=)`); false for position-specific.
    whole_protein: bool,
    /// True for predicted (`p.(=)`, `p.(Leu54=)`); false for certain.
    predicted: bool,
    /// True iff the outer `Mu` wrapper is `Mu::Uncertain`. See struct comment.
    outer_uncertain: bool,
}

/// All canonical D5 inputs, covering the 2x2 cross product (whole/position) x
/// (certain/predicted) plus a second position-specific accession (`Val600=`,
/// the BRAF-style canonical).
const CASES: &[SilentCase] = &[
    // Whole-protein, certain.  parser path: parse_whole_protein_identity → Mu::Certain
    SilentCase {
        input: "NP_000079.2:p.=",
        whole_protein: true,
        predicted: false,
        outer_uncertain: false,
    },
    // Whole-protein, predicted.  parser path: parse_whole_protein_identity → Mu::Certain,
    // predicted bit on the inner edit (asymmetric vs the position-specific predicted form below).
    SilentCase {
        input: "NP_000079.2:p.(=)",
        whole_protein: true,
        predicted: true,
        outer_uncertain: false,
    },
    // Position-specific, certain.
    SilentCase {
        input: "NP_000079.2:p.Leu54=",
        whole_protein: false,
        predicted: false,
        outer_uncertain: false,
    },
    // Position-specific, predicted.  parser path: predicted-protein-change branch wraps
    // the inner edit in Mu::Uncertain; the inner ProteinEdit::Identity carries
    // predicted: false.  This asymmetry vs p.(=) is intentional and pinned in the
    // shape test below.
    SilentCase {
        input: "NP_000079.2:p.(Leu54=)",
        whole_protein: false,
        predicted: true,
        outer_uncertain: true,
    },
    // BRAF-style canonical position-specific synonymous example.  Uses a
    // different accession to ensure the rule is accession-independent.
    SilentCase {
        input: "NP_004324.2:p.Val600=",
        whole_protein: false,
        predicted: false,
        outer_uncertain: false,
    },
];

/// Convenience: every input string from CASES, for tests that don't need the shape flags.
const SILENT_EQ_FORMS: &[&str] = &[
    "NP_000079.2:p.=",
    "NP_000079.2:p.(=)",
    "NP_000079.2:p.Leu54=",
    "NP_000079.2:p.(Leu54=)",
    "NP_004324.2:p.Val600=",
];

/// Borrow the inner `ProteinVariant` from an `HgvsVariant`, panicking with a
/// useful message for any other arm.
fn protein_variant(variant: &HgvsVariant) -> &ProteinVariant {
    match variant {
        HgvsVariant::Protein(p) => p,
        other => panic!("expected Protein variant, got {:?}", other),
    }
}

/// Borrow the inner `ProteinEdit`, panicking if the loc_edit somehow lacks one.
fn protein_edit(variant: &HgvsVariant) -> &ProteinEdit {
    protein_variant(variant)
        .loc_edit
        .edit
        .inner()
        .expect("silent-= variants always carry a concrete inner ProteinEdit")
}

/// Convenience: hash a value to a `u64` so we can compare hashes across passes.
fn hash_of<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

#[test]
fn silent_eq_round_trips_through_display() {
    for input in SILENT_EQ_FORMS {
        let parsed = parse_hgvs(input)
            .unwrap_or_else(|e| panic!("failed to parse silent-= form {input:?}: {e}"));
        let rendered = format!("{}", parsed);
        assert_eq!(
            rendered, *input,
            "Display did not round-trip silent = form {input:?}",
        );
    }
}

#[test]
fn silent_eq_parse_is_idempotent() {
    // parse(format(parse(x))) == parse(x). This catches subtle parser/Display
    // drift where a re-parse would produce a structurally different variant
    // (e.g. predicted vs. certain, position-specific vs. whole-protein).
    for input in SILENT_EQ_FORMS {
        let first = parse_hgvs(input).expect("first parse");
        let rendered = format!("{}", first);
        let second = parse_hgvs(&rendered).expect("re-parse");
        assert_eq!(
            first, second,
            "parse(format(parse({input:?}))) drifted from parse({input:?})",
        );
    }
}

#[test]
fn silent_eq_normalize_is_identity() {
    // Protein-level identity edits have nothing to shift or rewrite, so
    // Normalizer::normalize must be a no-op for each canonical form.
    let normalizer = Normalizer::new(MockProvider::new());
    for input in SILENT_EQ_FORMS {
        let variant = parse_hgvs(input).expect("parse");
        let normalized = normalizer
            .normalize(&variant)
            .expect("normalize should succeed for protein-identity edits");
        assert_eq!(
            normalized, variant,
            "normalize altered structure for silent = form {input:?}",
        );
        let rendered = format!("{}", normalized);
        assert_eq!(
            rendered, *input,
            "normalize altered Display for silent = form {input:?}",
        );

        // And one more round of normalize is also a fixed point.
        let renormalized = normalizer.normalize(&normalized).expect("renormalize");
        assert_eq!(
            renormalized, normalized,
            "normalize was not idempotent for silent = form {input:?}",
        );
    }
}

#[test]
fn silent_eq_classifies_each_case_correctly() {
    for case in CASES {
        let parsed = parse_hgvs(case.input)
            .unwrap_or_else(|e| panic!("parse({:?}) failed: {}", case.input, e));

        // (1) Must be a Protein variant.
        assert!(
            matches!(parsed, HgvsVariant::Protein(_)),
            "{:?} did not classify as HgvsVariant::Protein: got {:?}",
            case.input,
            parsed
        );

        // (2) Inner edit must be ProteinEdit::Identity with the expected flags.
        // Note: the position-specific predicted form (`p.(Leu54=)`) carries the
        // "predicted" qualifier on the OUTER Mu wrapper, NOT the inner Identity.
        // The inner Identity therefore has predicted=false even though the case
        // overall is predicted. This asymmetry is current behavior and is pinned
        // separately by silent_eq_internal_shape_is_pinned below.
        let edit = protein_edit(&parsed);
        match edit {
            ProteinEdit::Identity {
                predicted,
                whole_protein,
            } => {
                assert_eq!(
                    *whole_protein, case.whole_protein,
                    "{:?}: whole_protein mismatch on inner edit (got {}, want {})",
                    case.input, *whole_protein, case.whole_protein
                );
                // Inner `predicted` flag is set only when the parser took the
                // whole-protein-identity path (which encodes the parens on the
                // inner edit).  The position-specific predicted form takes the
                // outer-Mu::Uncertain path and leaves the inner predicted=false.
                let expected_inner_predicted = case.predicted && !case.outer_uncertain;
                assert_eq!(
                    *predicted, expected_inner_predicted,
                    "{:?}: inner predicted mismatch (got {}, want {})",
                    case.input, *predicted, expected_inner_predicted
                );
            }
            other => panic!(
                "{:?}: expected ProteinEdit::Identity, got {:?}",
                case.input, other
            ),
        }
    }
}

#[test]
fn silent_eq_internal_shape_is_pinned() {
    // The internal shape under silent-= is asymmetric:
    //
    // - `p.=`         → Mu::Certain(Identity { whole_protein: true,  predicted: false })
    //                   with the canonical Met1 dummy point interval.
    // - `p.(=)`       → Mu::Certain(Identity { whole_protein: true,  predicted: true })
    //                   with the canonical Met1 dummy point interval.
    //                   (Predicted bit lives on the inner edit; outer Mu is Certain.)
    // - `p.Leu54=`    → Mu::Certain(Identity { whole_protein: false, predicted: false })
    //                   with a real position interval (Leu54 point).
    // - `p.(Leu54=)`  → Mu::Uncertain(Identity { whole_protein: false, predicted: false })
    //                   with a real position interval. (Predicted bit lives on the
    //                   OUTER Mu wrapper here, not the inner edit. This asymmetry vs
    //                   `p.(=)` is current behavior; pinning it stops a refactor from
    //                   silently aligning the two paths in either direction without
    //                   updating both the parser and the Display logic.)
    //
    // A regression that flipped any of these (e.g. promoting `p.(=)` to
    // `Mu::Uncertain`, or demoting `p.(Leu54=)` to `Mu::Certain` with the
    // predicted bit on the inner edit) would Display the same string but break
    // structural-equality consumers (caches, sets, JSON serialisers).
    let met1 = ProtInterval::point(ProtPos::new(AminoAcid::Met, 1));

    for case in CASES {
        let parsed = parse_hgvs(case.input)
            .unwrap_or_else(|e| panic!("parse({:?}) failed: {}", case.input, e));
        let pv = protein_variant(&parsed);

        // (a) outer Mu wrapper matches expected.
        match (&pv.loc_edit.edit, case.outer_uncertain) {
            (Mu::Certain(_), false) | (Mu::Uncertain(_), true) => {}
            (Mu::Certain(_), true) => panic!(
                "{:?}: expected Mu::Uncertain outer wrapper, got Mu::Certain",
                case.input
            ),
            (Mu::Uncertain(_), false) => panic!(
                "{:?}: expected Mu::Certain outer wrapper, got Mu::Uncertain",
                case.input
            ),
            (other, _) => panic!("{:?}: unexpected Mu variant {:?}", case.input, other),
        }

        // (b) whole-protein cases use the canonical Met1 dummy point interval;
        //     position-specific cases must not collapse to it.
        if case.whole_protein {
            assert_eq!(
                pv.loc_edit.location, met1,
                "{:?}: expected Met1 dummy point interval, got {:?}",
                case.input, pv.loc_edit.location
            );
        } else {
            assert_ne!(
                pv.loc_edit.location, met1,
                "{:?}: position-specific form must NOT use the Met1 dummy interval",
                case.input
            );
        }
    }
}

#[test]
fn silent_eq_normalize_three_passes_with_hash_agreement() {
    // Three passes of Normalizer::normalize must each be a no-op AND each pass's
    // result must Hash equal to the previous pass.  `Hash` is derived but the
    // derive-blocker check isn't redundant: a future change to ProteinVariant or
    // ProteinEdit that broke `Hash` consistency (without breaking `Eq`) would
    // corrupt the normalizer's caches and any downstream HashSet/HashMap.
    let normalizer = Normalizer::new(MockProvider::new());

    for input in SILENT_EQ_FORMS {
        let parsed = parse_hgvs(input).expect("parse");
        let h0 = hash_of(&parsed);

        let pass1 = normalizer.normalize(&parsed).expect("pass1");
        assert_eq!(parsed, pass1, "{input:?}: pass1 altered structure");
        assert_eq!(h0, hash_of(&pass1), "{input:?}: pass1 hash drift");

        let pass2 = normalizer.normalize(&pass1).expect("pass2");
        assert_eq!(pass1, pass2, "{input:?}: pass2 altered structure");
        assert_eq!(h0, hash_of(&pass2), "{input:?}: pass2 hash drift");

        let pass3 = normalizer.normalize(&pass2).expect("pass3");
        assert_eq!(pass2, pass3, "{input:?}: pass3 altered structure");
        assert_eq!(h0, hash_of(&pass3), "{input:?}: pass3 hash drift");
    }
}

#[test]
fn silent_eq_hash_distinguishes_each_canonical_form() {
    // Four shapes — `p.=`, `p.(=)`, `p.Leu54=`, `p.(Leu54=)` — must hash distinctly
    // so HashSet/HashMap dedup keeps them apart.  This is a stronger guard than
    // Eq: a derive change that produced colliding hashes (without breaking Eq)
    // would silently corrupt downstream caches.
    let inputs = [
        "NP_000079.2:p.=",
        "NP_000079.2:p.(=)",
        "NP_000079.2:p.Leu54=",
        "NP_000079.2:p.(Leu54=)",
    ];
    let parsed: Vec<HgvsVariant> = inputs
        .iter()
        .map(|s| parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?}: {e}")))
        .collect();
    let hashes: Vec<u64> = parsed.iter().map(hash_of).collect();

    for i in 0..hashes.len() {
        for j in (i + 1)..hashes.len() {
            assert_ne!(
                hashes[i], hashes[j],
                "hash collision between {:?} and {:?}",
                inputs[i], inputs[j]
            );
        }
    }
}

#[test]
fn silent_eq_is_distinct_from_p0_and_q_forms() {
    // Section D's whole point is that each row is independently enforceable.  A
    // refactor that collapsed `Identity` into `NoProtein` or `Unknown` would
    // pass D5's existing round-trip tests (Display still works) but corrupt
    // semantics.  These cross-row guards stop that.
    let silent_whole_certain = parse_hgvs("NP_000079.2:p.=").unwrap();
    let silent_whole_predicted = parse_hgvs("NP_000079.2:p.(=)").unwrap();
    let silent_pos_certain = parse_hgvs("NP_000079.2:p.Leu54=").unwrap();
    let silent_pos_predicted = parse_hgvs("NP_000079.2:p.(Leu54=)").unwrap();

    // D6 territory.
    let no_protein_certain = parse_hgvs("NP_000079.2:p.0").unwrap();
    let no_protein_predicted = parse_hgvs("NP_000079.2:p.0?").unwrap();

    // D7 territory.
    let unknown_whole_certain = parse_hgvs("NP_000079.2:p.?").unwrap();
    let unknown_whole_predicted = parse_hgvs("NP_000079.2:p.(?)").unwrap();
    let unknown_pos_certain = parse_hgvs("NP_000079.2:p.Met1?").unwrap();
    let unknown_pos_predicted = parse_hgvs("NP_000079.2:p.(Met1?)").unwrap();

    let silents = [
        ("p.=", &silent_whole_certain),
        ("p.(=)", &silent_whole_predicted),
        ("p.Leu54=", &silent_pos_certain),
        ("p.(Leu54=)", &silent_pos_predicted),
    ];
    let others = [
        ("p.0", &no_protein_certain),
        ("p.0?", &no_protein_predicted),
        ("p.?", &unknown_whole_certain),
        ("p.(?)", &unknown_whole_predicted),
        ("p.Met1?", &unknown_pos_certain),
        ("p.(Met1?)", &unknown_pos_predicted),
    ];

    for (silent_label, silent) in &silents {
        for (other_label, other) in &others {
            assert_ne!(
                silent, other,
                "{} must remain structurally distinct from {}",
                silent_label, other_label
            );
            assert_ne!(
                format!("{}", silent),
                format!("{}", other),
                "{} must render distinctly from {}",
                silent_label,
                other_label
            );
        }
    }
}

#[test]
fn silent_eq_distinguishes_predicted_from_certain() {
    // The "predicted" parens are load-bearing; dropping them silently drops
    // information the spec sanctions for distinguishing observed vs predicted
    // consequences (recommendations/general/uncertainties.html).
    let whole_certain = parse_hgvs("NP_000079.2:p.=").unwrap();
    let whole_predicted = parse_hgvs("NP_000079.2:p.(=)").unwrap();
    assert_ne!(
        whole_certain, whole_predicted,
        "p.= and p.(=) must parse to distinct ASTs"
    );

    let pos_certain = parse_hgvs("NP_000079.2:p.Leu54=").unwrap();
    let pos_predicted = parse_hgvs("NP_000079.2:p.(Leu54=)").unwrap();
    assert_ne!(
        pos_certain, pos_predicted,
        "p.Leu54= and p.(Leu54=) must parse to distinct ASTs"
    );
}

#[test]
fn silent_eq_distinguishes_whole_from_position() {
    // Whole-protein and position-specific must stay distinct.  Collapsing
    // `p.Leu54=` onto `p.=` (or vice versa) loses the position information.
    let whole = parse_hgvs("NP_000079.2:p.=").unwrap();
    let position = parse_hgvs("NP_000079.2:p.Leu54=").unwrap();
    assert_ne!(
        whole, position,
        "p.= (whole-protein) and p.Leu54= (position-specific) must parse to distinct ASTs"
    );

    let whole_predicted = parse_hgvs("NP_000079.2:p.(=)").unwrap();
    let position_predicted = parse_hgvs("NP_000079.2:p.(Leu54=)").unwrap();
    assert_ne!(
        whole_predicted, position_predicted,
        "p.(=) and p.(Leu54=) must parse to distinct ASTs"
    );
}

#[test]
fn silent_eq_rejects_malformed_adjacent_inputs() {
    // The strings below are NOT in the HGVS spec for silent / synonymous
    // notation.  `parse_hgvs` requires consuming the whole input and rejects
    // each.  Pinning these stops a permissive future parser from accepting a
    // form that round-trips to a Display the spec doesn't sanction.
    //
    // Note: `parse_hgvs` trims OUTER whitespace but rejects INTERIOR whitespace
    // (matching D6's pinned behavior).  Bare inputs without an accession also
    // reject (the prefix `<acc>:p.` is required).
    let rejects = [
        // Trailing junk after a valid silent prefix.
        "NP_000079.2:p.Leu54==",
        "NP_000079.2:p.==",
        "NP_000079.2:p.=Leu54",
        "NP_000079.2:p.=X",
        // Reversed / un-spec parenthesised content.
        "NP_000079.2:p.(=Leu54)",
        "NP_000079.2:p.(Leu54=", // unclosed paren
        // Interior whitespace before the edit body.
        "NP_000079.2:p. =",
        "NP_000079.2:p. Leu54=",
        // Bare (no accession).
        "p.=",
        "p.(=)",
        "p.Leu54=",
        "p.(Leu54=)",
    ];

    for input in &rejects {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "{input:?} unexpectedly parsed: {:?}",
            result.ok()
        );
    }
}

/// Since #468, whole-entity protein edits (`=`, `(=)`, `(?)`) are
/// admitted as bracket members so compound forms like `p.[Ser68Arg];[=]`
/// parse (HGVS v21 `recommendations/protein/alleles.md`). Bare `?` is
/// NOT admitted here — it remains the trans whole-allele unknown marker
/// (`UnknownAllele`), handled outside `parse_protein_bracket_member`. As a
/// consequence, a *singleton* bracketed whole-protein identity `p.[=]`
/// now parses too — unwrapping to the spec-sanctioned `p.=`, exactly as
/// every other singleton bracket unwraps its sole member
/// (`p.[Ser68Arg]` → `p.Ser68Arg`). This replaces the prior pin that
/// rejected `p.[=]` outright (the rejection was an artifact of
/// whole-entity members not being admitted in bracket position at all).
#[test]
fn silent_eq_singleton_bracket_unwraps_to_whole_protein_identity() {
    let parsed = parse_hgvs("NP_000079.2:p.[=]").expect("p.[=] must parse since #468");
    assert_eq!(parsed.to_string(), "NP_000079.2:p.=");
    // Re-parse fixed point on the canonical (unwrapped) form.
    let reparsed = parse_hgvs("NP_000079.2:p.=").expect("p.= must parse");
    assert_eq!(reparsed.to_string(), "NP_000079.2:p.=");
}

#[test]
fn silent_eq_outer_whitespace_is_trimmed() {
    // Outer whitespace is tolerated by `parse_hgvs` (the variant-level parser
    // checks `final_remaining.trim().is_empty()` rather than strict
    // `is_empty()`, so trailing whitespace after a fully-formed variant is
    // silently swallowed).  This pins that behavior so a future tightening of
    // the trim doesn't silently regress.
    //
    // Note: leading whitespace inside `p.<edit>` (e.g. `p. =` or `p. Leu54=`)
    // is NOT tolerated — see `silent_eq_rejects_malformed_adjacent_inputs`.
    for input in SILENT_EQ_FORMS {
        let padded = format!("  {}  ", input);
        let parsed =
            parse_hgvs(&padded).unwrap_or_else(|e| panic!("parse({padded:?}) failed: {e}"));
        let direct = parse_hgvs(input).unwrap();
        assert_eq!(parsed, direct, "outer whitespace trim drift for {input:?}");
        assert_eq!(
            format!("{}", parsed),
            *input,
            "Display dropped accession or shape after outer trim of {input:?}"
        );

        // Trailing whitespace alone (no leading) is also tolerated.
        let trailing_only = format!("{} ", input);
        let parsed_trailing = parse_hgvs(&trailing_only)
            .unwrap_or_else(|e| panic!("parse({trailing_only:?}) failed: {e}"));
        assert_eq!(
            parsed_trailing, direct,
            "trailing whitespace trim drift for {input:?}"
        );
    }
}

#[test]
fn silent_eq_fast_path_matches_full_path() {
    // `parse_hgvs_fast` is a hot-path entry that should produce a structurally
    // identical AST to `parse_hgvs` for every D5 input.  Pinning this stops a
    // fast-path optimisation from silently diverging on the silent-= family
    // (e.g. shortcutting `p.=` to a different Mu wrapper).
    for input in SILENT_EQ_FORMS {
        let full =
            parse_hgvs(input).unwrap_or_else(|e| panic!("parse_hgvs({input:?}) failed: {e}"));
        let fast = parse_hgvs_fast(input)
            .unwrap_or_else(|e| panic!("parse_hgvs_fast({input:?}) failed: {e}"));
        assert_eq!(
            full, fast,
            "parse_hgvs vs parse_hgvs_fast diverged for {input:?}"
        );
        assert_eq!(
            format!("{}", full),
            format!("{}", fast),
            "Display drift between parse_hgvs and parse_hgvs_fast for {input:?}"
        );
        assert_eq!(
            hash_of(&full),
            hash_of(&fast),
            "Hash drift between parse_hgvs and parse_hgvs_fast for {input:?}"
        );
    }
}

#[test]
fn silent_eq_trans_allele_compact_form_round_trips() {
    // PR #146 (sibling C1) added the compact-prefix protein trans-allele
    // shorthand `ACC:p.[a];[b]`, fixing the parse-error drift previously
    // tracked under issue #83.  Both expanded and compact forms now parse
    // and the AST round-trips structurally:
    //
    //   parse("[ACC:p.Leu54=];[ACC:p.Arg97Trp]")
    //     == parse("ACC:p.[Leu54=];[Arg97Trp]")
    //
    // The cis form `ACC:p.[Leu54=;Arg97Trp]` (single bracket pair,
    // semicolons inside) was already working and is pinned alongside.
    let expanded = "[NP_000079.2:p.Leu54=];[NP_000079.2:p.Arg97Trp]";
    let parsed = parse_hgvs(expanded).unwrap_or_else(|e| panic!("parse({expanded:?}) failed: {e}"));

    let displayed = format!("{}", parsed);
    let compact = "NP_000079.2:p.[Leu54=];[Arg97Trp]";
    assert_eq!(
        displayed, compact,
        "expanded trans allele did not Display as the compact form"
    );

    let reparsed = parse_hgvs(compact)
        .unwrap_or_else(|e| panic!("parse({compact:?}) failed after PR #146: {e}"));
    assert_eq!(
        parsed, reparsed,
        "compact-form trans allele AST drifted on re-parse"
    );
    assert_eq!(
        format!("{}", reparsed),
        compact,
        "compact-form trans allele Display did not round-trip"
    );

    // Cis compact form does round-trip cleanly today; pin that contract too.
    let cis = "NP_000079.2:p.[Leu54=;Arg97Trp]";
    let cis_parsed = parse_hgvs(cis).unwrap_or_else(|e| panic!("parse({cis:?}) failed: {e}"));
    assert_eq!(
        format!("{}", cis_parsed),
        cis,
        "cis-form silent + missense allele Display did not round-trip"
    );
    let cis_reparsed = parse_hgvs(cis).unwrap();
    assert_eq!(
        cis_parsed, cis_reparsed,
        "cis-form silent + missense allele AST drifted on re-parse"
    );
}
