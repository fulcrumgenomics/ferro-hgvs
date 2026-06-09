//! Round-trip pinning for protein unknown-effect notation (`p.?` family) — issue #81 D7.
//!
//! HGVS spec context (`assets/hgvs-nomenclature/docs/recommendations/protein/`):
//!
//! - `p.?` — translation product is unknown (whole-protein, certain).
//! - `p.(?)` — predicted whole-protein unknown.
//! - `p.Met1?` — initiation Met affected, effect unknown (position-specific,
//!   certain).
//! - `p.(Met1?)` — predicted unknown effect on initiation Met. This is the
//!   form the spec uses in `substitution.md:43`: "the consequence, at the
//!   protein level, of a variant affecting the translation initiation codon
//!   can not be predicted (i.e. is unknown)".
//!
//! Distinct from `p.0` / `p.0?` (no protein produced — D6) and from
//! `p.=` / `p.(=)` / `p.Leu54=` / `p.(Leu54=)` (silent — D5). Each row of
//! Section D is independently enforceable; this file pins every D7 invariant
//! the existing parse / Display / normalize stack provides today, so a future
//! refactor that collapses any pair (predicted↔certain, whole↔position,
//! unknown↔no-protein, unknown↔silent) fails loudly here.
//!
//! Spec corpus cross-link: `tests/fixtures/grammar/hgvs_spec_normalization.json`
//! already covers `p.?`, `p.(?)`, `LRG_199p1:p.(Met1?)`, `NP_003997.1:p.?`,
//! and `NP_003997.2:p.?` with `status: "preserved"`. This file does not
//! duplicate the v21.0 corpus rows — it pins behavior the corpus harness
//! cannot (Mu shape, Hash agreement, cross-form non-equivalence, the
//! `[?]`-vs-`p.(?)` distinction, the trans compact-form drift gate).
//!
//! Sibling pinners: `tests/protein_silent_eq.rs` (D5) and
//! `tests/protein_no_protein_roundtrip.rs` (D6).

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use ferro_hgvs::hgvs::edit::ProteinEdit;
use ferro_hgvs::hgvs::interval::ProtInterval;
use ferro_hgvs::hgvs::location::{AminoAcid, ProtPos};
use ferro_hgvs::hgvs::uncertainty::Mu;
use ferro_hgvs::hgvs::variant::{AllelePhase, ProteinVariant};
use ferro_hgvs::{parse_hgvs, parse_hgvs_fast, HgvsVariant, MockProvider, Normalizer};

/// One canonical D7 input plus the structural shape we expect after parse.
///
/// `outer_uncertain`: true iff the parser wraps the inner `ProteinEdit` in
/// `Mu::Uncertain`. Currently only `p.(Met1?)` — the position-specific
/// predicted form — takes that path; the whole-protein `p.(?)` keeps the
/// `Mu::Certain` outer with `predicted: true` on the inner edit. This
/// asymmetry mirrors the D5 silent-`=` family and is intentional.
struct UnknownCase {
    /// HGVS input (always with an accession prefix).
    input: &'static str,
    /// True for whole-protein (`p.?`, `p.(?)`); false for position-specific.
    whole_protein: bool,
    /// True for predicted (`p.(?)`, `p.(Met1?)`); false for certain.
    predicted: bool,
    /// True iff the outer `Mu` wrapper is `Mu::Uncertain`. See struct comment.
    outer_uncertain: bool,
}

/// All canonical D7 inputs, covering the 2x2 cross product
/// (whole/position) x (certain/predicted).
///
/// We use `NP_000079.2` (the same NP accession exercised by the existing
/// parser unit tests in `src/hgvs/parser/variant.rs`) so this fixture is
/// self-contained and does not depend on any reference data.
const CASES: &[UnknownCase] = &[
    // Whole-protein, certain. parser path: parse_whole_protein_unknown -> Mu::Certain
    UnknownCase {
        input: "NP_000079.2:p.?",
        whole_protein: true,
        predicted: false,
        outer_uncertain: false,
    },
    // Whole-protein, predicted. parser path: parse_whole_protein_unknown matches
    // "(?)" via tag("(?)") and returns `Unknown { predicted: true, whole_protein: true }`,
    // which is wrapped in Mu::Certain (NOT Mu::Uncertain) — the predicted bit
    // lives on the inner edit. Asymmetric vs the position-specific predicted
    // form below.
    UnknownCase {
        input: "NP_000079.2:p.(?)",
        whole_protein: true,
        predicted: true,
        outer_uncertain: false,
    },
    // Position-specific, certain. parser path: position+edit ->
    // parse_protein_unknown -> position_unknown() -> Mu::Certain.
    UnknownCase {
        input: "NP_000079.2:p.Met1?",
        whole_protein: false,
        predicted: false,
        outer_uncertain: false,
    },
    // Position-specific, predicted. parser path: predicted-wrapping `(...)`
    // branch (variant.rs:1771) wraps in Mu::Uncertain. The inner edit's
    // `predicted` field is FALSE — the predicted qualifier lives on the OUTER
    // Mu::Uncertain wrapper, not the inner edit.
    UnknownCase {
        input: "NP_000079.2:p.(Met1?)",
        whole_protein: false,
        predicted: true,
        outer_uncertain: true,
    },
];

/// Helper: return the `(ProteinVariant, &Mu<ProteinEdit>, &ProteinEdit)`
/// triple for an input we expect to parse to a single (non-allele) protein
/// variant. Panics with a clear message otherwise.
fn parse_as_protein_variant(input: &str) -> ProteinVariant {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse({input:?}) failed: {e}"));
    match parsed {
        HgvsVariant::Protein(v) => v,
        other => panic!("parse({input:?}) -> {other:?}, expected HgvsVariant::Protein"),
    }
}

/// Hash of any `Hash` value via `DefaultHasher` (used for `Hash`-agreement
/// assertions independent of `Eq`).
fn hash_of<T: Hash>(value: &T) -> u64 {
    let mut h = DefaultHasher::new();
    value.hash(&mut h);
    h.finish()
}

/// Run the full round-trip battery on a single input string.
///
/// All four `p.?` shapes are ref-independent (no sequence is consulted), so
/// an empty `MockProvider` is sufficient for the normalize step.
fn assert_roundtrips(input: &str) {
    // (1) parse succeeds
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse({input:?}) failed: {e}"));

    // The variant must classify as a protein variant — the empty-MockProvider
    // shortcut below relies on this, and a regression that demoted `p.?` to
    // (say) `UnknownAllele` would also be a real bug.
    assert!(
        matches!(parsed, HgvsVariant::Protein(_)),
        "parse({input:?}) -> {parsed:?}: expected HgvsVariant::Protein"
    );

    // (2) Display is idempotent: format(parse(s)) == s
    let displayed = format!("{}", parsed);
    assert_eq!(
        displayed, input,
        "Display round-trip lost information: parse({input:?}) displays as {displayed:?}"
    );

    // (3) re-parse: parse(format(parse(s))) == parse(s)
    let reparsed =
        parse_hgvs(&displayed).unwrap_or_else(|e| panic!("re-parse({displayed:?}) failed: {e}"));
    assert_eq!(
        reparsed, parsed,
        "AST drifted on re-parse: original={parsed:?}, reparsed={reparsed:?}"
    );

    // (4) normalize is a no-op for all four shapes (no ref-dependent rewrite).
    //     `MockProvider::new()` has no protein data, so the protein-reference
    //     validator short-circuits and `normalize_protein` falls through to
    //     the pass-through arm.
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize({input:?}) failed: {e}"));
    assert_eq!(
        normalized, parsed,
        "normalize({input:?}) was not a no-op: before={parsed:?}, after={normalized:?}"
    );

    // (5) the normalized variant still displays as the original input, and
    //     re-parsing it yields the same AST. This guards against a future
    //     change that rewrites the AST during normalize but happens to land
    //     on a Display-equivalent shape.
    let normalized_displayed = format!("{}", normalized);
    assert_eq!(
        normalized_displayed, input,
        "Display after normalize differs from input: got {normalized_displayed:?}, want {input:?}"
    );
    let normalized_reparsed = parse_hgvs(&normalized_displayed).unwrap_or_else(|e| {
        panic!("re-parse-after-normalize({normalized_displayed:?}) failed: {e}")
    });
    assert_eq!(
        normalized_reparsed, parsed,
        "AST drifted after normalize+display+reparse: before={parsed:?}, after={normalized_reparsed:?}"
    );
}

// =====================================================================
// Round-trip core: each canonical form parses, Displays, re-parses, and
// normalizes as a no-op. Carried over from the prior PR revision.
// =====================================================================

#[test]
fn protein_unknown_whole_certain_roundtrips() {
    // p.? — whole-protein, certain ("unknown effect").
    assert_roundtrips("NP_000079.2:p.?");
}

#[test]
fn protein_unknown_whole_predicted_roundtrips() {
    // p.(?) — whole-protein, predicted unknown.
    // Distinct from p.? and must NOT be flattened to it (the parens carry the
    // "predicted" qualifier and dropping them silently loses information).
    assert_roundtrips("NP_000079.2:p.(?)");
}

#[test]
fn protein_unknown_position_certain_roundtrips() {
    // p.Met1? — position-specific, certain. Initiation Met affected, effect
    // unknown.
    assert_roundtrips("NP_000079.2:p.Met1?");
}

#[test]
fn protein_unknown_position_predicted_roundtrips() {
    // p.(Met1?) — position-specific, predicted. This is the form the HGVS
    // spec gives in `protein/substitution.md:43`.
    assert_roundtrips("NP_000079.2:p.(Met1?)");
}

#[test]
fn protein_unknown_all_canonical_forms_roundtrip() {
    // Belt-and-braces: drive every canonical form through the same battery
    // in one place. If a future contributor adds a fifth canonical shape (or
    // collapses one), this loop is the single edit point.
    for case in CASES {
        assert_roundtrips(case.input);
    }
}

// =====================================================================
// Edit-shape pin — each canonical form's `ProteinEdit::Unknown` carries
// the expected `(predicted, whole_protein)` flags. A regression that
// flattens any pair would fail here even before round-trip noticed.
// =====================================================================

#[test]
fn protein_unknown_edit_shape_is_pinned() {
    for case in CASES {
        let pv = parse_as_protein_variant(case.input);
        let inner = pv.loc_edit.edit.inner().unwrap_or_else(|| {
            panic!(
                "edit on {:?} should not be Mu::Unknown; got {:?}",
                case.input, pv.loc_edit.edit
            )
        });
        match inner {
            ProteinEdit::Unknown {
                predicted,
                whole_protein,
            } => {
                assert_eq!(
                    *whole_protein, case.whole_protein,
                    "{:?}: whole_protein mismatch",
                    case.input
                );
                // The `predicted` field on the inner edit ONLY carries the
                // predicted qualifier when the outer Mu is Certain (i.e. for
                // the whole-protein predicted form `p.(?)`). For the position-
                // specific predicted form `p.(Met1?)` the predicted bit lives
                // on the OUTER `Mu::Uncertain` wrapper instead. So the inner
                // `predicted` field equals `case.predicted XOR case.outer_uncertain`.
                let expected_inner_predicted = case.predicted ^ case.outer_uncertain;
                assert_eq!(
                    *predicted, expected_inner_predicted,
                    "{:?}: inner predicted bit mismatch (expected {} = predicted({}) XOR outer_uncertain({}))",
                    case.input, expected_inner_predicted, case.predicted, case.outer_uncertain
                );
            }
            other => panic!(
                "{:?}: expected ProteinEdit::Unknown, got {other:?}",
                case.input
            ),
        }
    }
}

// =====================================================================
// Internal `Mu` shape pin — asymmetric, intentional. Locks down the
// outer-Mu wrapper variant for each canonical form so a refactor that
// aligns the two predicted paths in either direction fails loudly.
// =====================================================================

#[test]
fn protein_unknown_internal_mu_shape_is_pinned() {
    for case in CASES {
        let pv = parse_as_protein_variant(case.input);
        match (case.outer_uncertain, &pv.loc_edit.edit) {
            (false, Mu::Certain(_)) | (true, Mu::Uncertain(_)) => {} // expected
            (false, other) => panic!(
                "{:?}: expected outer Mu::Certain, got {other:?}",
                case.input
            ),
            (true, other) => panic!(
                "{:?}: expected outer Mu::Uncertain, got {other:?}",
                case.input
            ),
        }
        // The outer Mu must NEVER be `Mu::Unknown` — `?` is a property of the
        // protein edit, not of the outer Mu wrapper. `Mu::Unknown` is reserved
        // for nucleotide-edit `c.?` / `r.?` paths (variant.rs:701).
        assert!(
            !matches!(pv.loc_edit.edit, Mu::Unknown),
            "{:?}: outer Mu must not be Mu::Unknown; that variant is for NaEdit, not ProteinEdit",
            case.input
        );
    }
}

// =====================================================================
// Dummy `Met1` interval pin — every canonical D7 form (whole AND
// position-specific) carries the canonical `ProtPos::new(Met, 1)` point
// interval. For whole-protein forms this is the parser's dummy fill
// (variant.rs:1749-1750); for position-specific forms it happens because
// the input literally says `Met1`. These come from different parser
// paths but produce equal intervals.
// =====================================================================

#[test]
fn protein_unknown_met1_dummy_interval_is_pinned() {
    let canonical_met1_interval = ProtInterval::point(ProtPos::new(AminoAcid::Met, 1));
    for case in CASES {
        let pv = parse_as_protein_variant(case.input);
        assert_eq!(
            pv.loc_edit.location, canonical_met1_interval,
            "{:?}: expected canonical Met1 point interval, got {:?}",
            case.input, pv.loc_edit.location
        );
    }
}

// =====================================================================
// Multi-pass normalize idempotency + Hash agreement.
//
// Three passes guard against a future change that converges only after the
// second iteration. Hash agreement guards `Hash`/`Eq` drift independent of
// `Eq` (relevant for normalizer caches and downstream `HashMap`/`HashSet`
// dedup; `Hash` is derived but a derive-blocker change wouldn't be caught
// by `Eq`-only assertions).
// =====================================================================

#[test]
fn protein_unknown_normalize_is_multi_pass_idempotent() {
    let normalizer = Normalizer::new(MockProvider::new());
    for case in CASES {
        let parsed = parse_hgvs(case.input).unwrap();
        let pass1 = normalizer.normalize(&parsed).unwrap();
        let pass2 = normalizer.normalize(&pass1).unwrap();
        let pass3 = normalizer.normalize(&pass2).unwrap();

        assert_eq!(parsed, pass1, "{:?}: pass1 differed", case.input);
        assert_eq!(pass1, pass2, "{:?}: pass2 differed", case.input);
        assert_eq!(pass2, pass3, "{:?}: pass3 differed", case.input);

        let h0 = hash_of(&parsed);
        assert_eq!(hash_of(&pass1), h0, "{:?}: hash drift on pass1", case.input);
        assert_eq!(hash_of(&pass2), h0, "{:?}: hash drift on pass2", case.input);
        assert_eq!(hash_of(&pass3), h0, "{:?}: hash drift on pass3", case.input);
    }
}

// =====================================================================
// Hash distinguishes the four canonical shapes — guards downstream
// HashSet / HashMap dedup against a regression that hashes equal but
// compares non-equal (or vice versa).
// =====================================================================

#[test]
fn protein_unknown_hash_distinguishes_canonical_shapes() {
    let parsed: Vec<HgvsVariant> = CASES.iter().map(|c| parse_hgvs(c.input).unwrap()).collect();
    for (i, ai) in parsed.iter().enumerate() {
        for (j, aj) in parsed.iter().enumerate() {
            if i == j {
                continue;
            }
            assert_ne!(
                hash_of(ai),
                hash_of(aj),
                "{:?} and {:?} hash-collide",
                CASES[i].input,
                CASES[j].input
            );
        }
    }
}

// =====================================================================
// Cross-form non-equivalence (D5 / D6 guards).
//
// Every D7 form is structurally and Display-distinct from every D5 silent
// form (`p.=`, `p.(=)`, `p.Leu54=`, `p.(Leu54=)`) and every D6 no-protein
// form (`p.0`, `p.0?`). Without these, a regression that collapsed
// `Unknown` onto `Identity` or `NoProtein` would slip past D7's existing
// round-trip tests.
// =====================================================================

#[test]
fn protein_unknown_distinct_from_d5_silent_and_d6_no_protein() {
    let d5_inputs = [
        "NP_000079.2:p.=",
        "NP_000079.2:p.(=)",
        "NP_000079.2:p.Leu54=",
        "NP_000079.2:p.(Leu54=)",
    ];
    let d6_inputs = ["NP_000079.2:p.0", "NP_000079.2:p.0?"];

    for case in CASES {
        let d7_parsed = parse_hgvs(case.input).unwrap();
        let d7_display = format!("{}", d7_parsed);

        for other in d5_inputs.iter().chain(d6_inputs.iter()) {
            let other_parsed = parse_hgvs(other).unwrap();
            assert_ne!(
                d7_parsed, other_parsed,
                "D7 form {:?} structurally equals D5/D6 form {other:?}",
                case.input
            );
            assert_ne!(
                d7_display,
                format!("{}", other_parsed),
                "D7 form {:?} Display-equals D5/D6 form {other:?}",
                case.input
            );
            assert_ne!(
                hash_of(&d7_parsed),
                hash_of(&other_parsed),
                "D7 form {:?} hashes equal to D5/D6 form {other:?}",
                case.input
            );
        }
    }
}

// =====================================================================
// Predicted ↔ certain non-equivalence — predicted parens are load-bearing.
// =====================================================================

#[test]
fn protein_unknown_predicted_distinct_from_certain() {
    // Whole-protein
    let certain = parse_hgvs("NP_000079.2:p.?").unwrap();
    let predicted = parse_hgvs("NP_000079.2:p.(?)").unwrap();
    assert_ne!(certain, predicted, "p.? must differ from p.(?)");
    assert_ne!(
        format!("{certain}"),
        format!("{predicted}"),
        "p.? must Display differently from p.(?)"
    );

    // Position-specific
    let pos_certain = parse_hgvs("NP_000079.2:p.Met1?").unwrap();
    let pos_predicted = parse_hgvs("NP_000079.2:p.(Met1?)").unwrap();
    assert_ne!(
        pos_certain, pos_predicted,
        "p.Met1? must differ from p.(Met1?)"
    );
    assert_ne!(
        format!("{pos_certain}"),
        format!("{pos_predicted}"),
        "p.Met1? must Display differently from p.(Met1?)"
    );
}

// =====================================================================
// Whole-protein ↔ position-specific non-equivalence.
// =====================================================================

#[test]
fn protein_unknown_whole_distinct_from_position() {
    let pairs = [
        ("NP_000079.2:p.?", "NP_000079.2:p.Met1?"),
        ("NP_000079.2:p.(?)", "NP_000079.2:p.(Met1?)"),
    ];
    for (whole_input, pos_input) in pairs {
        let whole = parse_hgvs(whole_input).unwrap();
        let pos = parse_hgvs(pos_input).unwrap();
        assert_ne!(
            whole, pos,
            "{whole_input:?} (whole) must differ from {pos_input:?} (position-specific)"
        );
        assert_ne!(
            format!("{whole}"),
            format!("{pos}"),
            "{whole_input:?} must Display differently from {pos_input:?}"
        );
    }
}

// =====================================================================
// `parse_hgvs_fast` parity for every D7 input.
// =====================================================================

#[test]
fn protein_unknown_parse_hgvs_fast_matches_normal() {
    for case in CASES {
        let normal = parse_hgvs(case.input).unwrap();
        let fast = parse_hgvs_fast(case.input)
            .unwrap_or_else(|e| panic!("parse_hgvs_fast({:?}) failed: {e}", case.input));
        assert_eq!(
            fast, normal,
            "{:?}: parse_hgvs_fast and parse_hgvs disagree",
            case.input
        );
        assert_eq!(
            hash_of(&fast),
            hash_of(&normal),
            "{:?}: hash disagreement between parse_hgvs_fast and parse_hgvs",
            case.input
        );
    }
}

// =====================================================================
// Outer whitespace tolerance — `parse_hgvs` swallows leading and trailing
// whitespace. Pinned so a future tightening doesn't silently regress
// (matches D5 / D6).
// =====================================================================

#[test]
fn protein_unknown_outer_whitespace_is_trimmed() {
    let canonical = parse_hgvs("NP_000079.2:p.?").unwrap();

    let leading = parse_hgvs("  NP_000079.2:p.?").unwrap();
    assert_eq!(leading, canonical, "leading whitespace was not trimmed");

    let trailing = parse_hgvs("NP_000079.2:p.?  ").unwrap();
    assert_eq!(trailing, canonical, "trailing whitespace was not trimmed");

    let both = parse_hgvs("  NP_000079.2:p.?  ").unwrap();
    assert_eq!(
        both, canonical,
        "leading+trailing whitespace was not trimmed"
    );
}

// =====================================================================
// Malformed / adjacent-input rejection. Each input must be `is_err()`.
// =====================================================================

#[test]
fn protein_unknown_malformed_inputs_reject() {
    // (label, input) — the label is for the panic message.
    let bad = [
        // Trailing junk after the canonical form.
        ("p.?? (double ?)", "NP_000079.2:p.??"),
        ("p.?X (junk)", "NP_000079.2:p.?X"),
        ("p.?Met1 (reversed)", "NP_000079.2:p.?Met1"),
        ("p.(?)? (junk after predicted)", "NP_000079.2:p.(?)?"),
        // Interior whitespace.
        ("p. ? (interior space, whole)", "NP_000079.2:p. ?"),
        ("p. (?) (interior space, predicted)", "NP_000079.2:p. (?)"),
        (
            "p.Met1 ? (interior space, position)",
            "NP_000079.2:p.Met1 ?",
        ),
        // Unclosed parens.
        ("p.(? (unclosed)", "NP_000079.2:p.(?"),
        ("p.(Met1? (unclosed)", "NP_000079.2:p.(Met1?"),
        // Position 0 — invalid (a `Met0?` is meaningless; positions are 1-based).
        ("p.Met0? (zero position)", "NP_000079.2:p.Met0?"),
        // No accession.
        ("bare p.? (no accession)", "p.?"),
        ("bare p.Met1? (no accession)", "p.Met1?"),
        // Bracket forms — the spec only sanctions `[?]` as the trans whole-
        // allele marker (variant.rs:2553-2554), never inside the protein
        // bracket-shorthand path. `parse_protein_allele_shorthand` requires
        // position+edit on each element and rejects bare `?`.
        ("p.[?] (singleton bracket)", "NP_000079.2:p.[?]"),
        (
            "p.[?;Arg97Trp] (cis with bare ?)",
            "NP_000079.2:p.[?;Arg97Trp]",
        ),
        (
            "p.[Arg97Trp;?] (cis with bare ? trailing)",
            "NP_000079.2:p.[Arg97Trp;?]",
        ),
        // Empty.
        ("p. (empty)", "NP_000079.2:p."),
    ];
    for (label, input) in bad {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "expected parse error for {label} ({input:?}), got Ok({:?})",
            result.ok()
        );
    }
}

// =====================================================================
// Display-canonicalization: single-letter and lowercase amino acid codes
// in `p.M1?` / `p.m1?` are accepted by the parser and Display-canonicalize
// to the spec's 3-letter form. Pinned so a future change either upholds
// the canonicalization (preferred) or fails loudly if it stops accepting
// the lenient input forms.
// =====================================================================

#[test]
fn protein_unknown_single_letter_residue_canonicalizes_to_three_letter() {
    // p.M1? -> p.Met1?
    let ml = parse_hgvs("NP_000079.2:p.M1?").unwrap();
    let canonical = parse_hgvs("NP_000079.2:p.Met1?").unwrap();
    assert_eq!(ml, canonical, "p.M1? must parse to the same AST as p.Met1?");
    assert_eq!(
        format!("{ml}"),
        "NP_000079.2:p.Met1?",
        "p.M1? must Display in the spec's 3-letter form"
    );

    // p.(M1?) -> p.(Met1?)
    let mlp = parse_hgvs("NP_000079.2:p.(M1?)").unwrap();
    let canonical_p = parse_hgvs("NP_000079.2:p.(Met1?)").unwrap();
    assert_eq!(mlp, canonical_p);
    assert_eq!(format!("{mlp}"), "NP_000079.2:p.(Met1?)");

    // Lowercase is also accepted (the underlying amino-acid lookup is
    // case-insensitive).
    let lc = parse_hgvs("NP_000079.2:p.met1?").unwrap();
    assert_eq!(lc, canonical, "p.met1? must canonicalize to p.Met1?");
    assert_eq!(format!("{lc}"), "NP_000079.2:p.Met1?");
}

// =====================================================================
// `[?]` is the spec's whole-allele unknown marker (alleles.md:43) and
// parses to `HgvsVariant::UnknownAllele`, NOT a `ProteinEdit::Unknown`
// arm. They are semantically distinct: a refactor that folded them
// would corrupt downstream interpretation (e.g. recessive-disease
// reporting where one allele is "expected but not yet identified" gets
// miscoded as "this allele has the protein-level unknown variant `p.(?)`").
// Mirrors D6's `[0]` vs `p.0` finding.
// =====================================================================

#[test]
fn protein_unknown_bracket_marker_is_unknown_allele_not_protein_edit() {
    // The spec's canonical heterozygous example (alleles.md:43):
    //   NP_003997.1:p.[(Ser68Arg)];[?]
    // Arm 0 is a predicted protein substitution, arm 1 is the bare
    // `[?]` whole-allele marker. Compact form: the shared accession is
    // written once (the `[?]` marker contributes none).
    let input = "NP_003997.1:p.[(Ser68Arg)];[?]";
    let parsed = parse_hgvs(input).unwrap();

    let allele = match &parsed {
        HgvsVariant::Allele(a) => a,
        other => panic!("expected HgvsVariant::Allele, got {other:?}"),
    };
    assert_eq!(allele.phase, AllelePhase::Trans, "phase must be Trans");
    assert_eq!(allele.variants.len(), 2, "expected two arms");

    // Arm 0: a protein variant.
    assert!(
        matches!(allele.variants[0], HgvsVariant::Protein(_)),
        "arm 0 must be a Protein variant; got {:?}",
        allele.variants[0]
    );

    // Arm 1: bare `[?]` is the whole-allele unknown marker, NOT a
    // ProteinEdit::Unknown variant.
    assert!(
        matches!(allele.variants[1], HgvsVariant::UnknownAllele),
        "arm 1 must be HgvsVariant::UnknownAllele (the spec's whole-allele marker), got {:?}",
        allele.variants[1]
    );

    // Distinct from a protein variant carrying `ProteinEdit::Unknown`.
    let protein_unknown_arm = parse_hgvs("NP_003997.1:p.(?)").unwrap();
    assert_ne!(
        allele.variants[1], protein_unknown_arm,
        "[?] (UnknownAllele) must NOT equal p.(?) (ProteinEdit::Unknown)"
    );

    // Display + reparse is the identity.
    let displayed = format!("{}", parsed);
    assert_eq!(displayed, input, "Display drift on [?]-marker compound");
    let reparsed = parse_hgvs(&displayed).unwrap();
    assert_eq!(reparsed, parsed, "AST drift on [?]-marker compound reparse");
}

// =====================================================================
// Trans-allele compound: when any arm is a WHOLE-protein unknown
// (`p.?` or `p.(?)`), the compact form is suppressed by
// `is_loc_edit_unknown` (variant.rs:691-708) — Display emits the
// expanded form `[ACC:p.X];[ACC:p.?]`, which round-trips cleanly. This
// is the intended contract and is the inverse of the position-
// specific drift pinned below.
// =====================================================================

#[test]
fn protein_unknown_trans_compound_with_whole_protein_unknown_arm_round_trips() {
    let inputs = [
        "[NP_000079.2:p.Arg97Trp];[NP_000079.2:p.?]",
        "[NP_000079.2:p.Arg97Trp];[NP_000079.2:p.(?)]",
    ];
    for input in inputs {
        let parsed = parse_hgvs(input).unwrap();

        // Phase = Trans.
        let allele = match &parsed {
            HgvsVariant::Allele(a) => a,
            other => panic!("{input:?}: expected HgvsVariant::Allele, got {other:?}"),
        };
        assert_eq!(
            allele.phase,
            AllelePhase::Trans,
            "{input:?}: phase mismatch"
        );
        assert_eq!(allele.variants.len(), 2, "{input:?}: arm count mismatch");

        // Arm 1 must carry a `ProteinEdit::Unknown` with `whole_protein: true`.
        let arm1 = match &allele.variants[1] {
            HgvsVariant::Protein(v) => v,
            other => panic!("{input:?}: arm 1 expected Protein, got {other:?}"),
        };
        match arm1.loc_edit.edit.inner() {
            Some(ProteinEdit::Unknown {
                whole_protein: true,
                ..
            }) => {} // expected
            other => panic!(
                "{input:?}: arm 1 inner edit expected Unknown(whole_protein=true), got {other:?}"
            ),
        }

        // Compact-form suppression: Display preserves the expanded form.
        let displayed = format!("{}", parsed);
        assert_eq!(
            displayed, input,
            "{input:?}: Display did NOT preserve expanded form (compact suppression broke)"
        );

        // Reparse equals original.
        let reparsed = parse_hgvs(&displayed).unwrap();
        assert_eq!(reparsed, parsed, "{input:?}: AST drifted on reparse");
    }
}

// =====================================================================
// Trans-allele compound: when an arm is a POSITION-SPECIFIC unknown
// (`p.Met1?`), `is_whole_protein_unknown()` returns false, so
// compact-form suppression does NOT trigger — Display produces the
// compact form `ACC:p.[Arg97Trp];[Met1?]`.
//
// PR #146 (sibling C1) added the compact-prefix protein trans-allele
// shorthand parser, so the compact form now round-trips structurally
// (previously this was tracked under issue #83 as a parse-error drift).
// =====================================================================

#[test]
fn protein_unknown_trans_compound_with_position_unknown_arm_round_trips() {
    let input = "[NP_000079.2:p.Arg97Trp];[NP_000079.2:p.Met1?]";
    let parsed = parse_hgvs(input).unwrap();

    // Phase = Trans, arm 1 is position-specific Unknown.
    let allele = match &parsed {
        HgvsVariant::Allele(a) => a,
        other => panic!("{input:?}: expected HgvsVariant::Allele, got {other:?}"),
    };
    assert_eq!(
        allele.phase,
        AllelePhase::Trans,
        "{input:?}: phase mismatch"
    );

    let arm1 = match &allele.variants[1] {
        HgvsVariant::Protein(v) => v,
        other => panic!("{input:?}: arm 1 expected Protein, got {other:?}"),
    };
    match arm1.loc_edit.edit.inner() {
        Some(ProteinEdit::Unknown {
            whole_protein: false,
            ..
        }) => {} // expected
        other => panic!(
            "{input:?}: arm 1 inner edit expected Unknown(whole_protein=false), got {other:?}"
        ),
    }

    // Display compacts to ACC:p.[X];[...] and round-trips structurally.
    let displayed = format!("{}", parsed);
    assert!(
        displayed.starts_with("NP_000079.2:p.["),
        "{input:?}: expected compact same-accession Display, got {displayed:?}"
    );

    let reparsed = parse_hgvs(&displayed)
        .unwrap_or_else(|e| panic!("{input:?}: parse({displayed:?}) failed after PR #146: {e}"));
    assert_eq!(
        reparsed, parsed,
        "{input:?}: compact-form Display {displayed:?} did not round-trip"
    );
}

// =====================================================================
// Trans-allele compound, predicted-uncertainty arm `(Met1?)`: the
// predicted-uncertainty wrapper `(...)` is NOT recognized inside any
// allele-shorthand bracket today (cis `p.[(Met1?);...]` and trans
// `p.[X];[(Met1?)]` both fail to re-parse). PR #146 explicitly scoped
// out predicted-form arms; this gap is independent of the same-accession
// trans compact-form fix and remains tracked under issue #83 alongside
// the cis bracket path.
//
// We pin the drift here so a future fix (extending the bracket parsers
// to dispatch on `(`) lights this up and forces a flip to round-trip.
// =====================================================================

#[test]
fn protein_unknown_trans_compound_with_predicted_arm_round_trips() {
    let input = "[NP_000079.2:p.Arg97Trp];[NP_000079.2:p.(Met1?)]";
    let parsed = parse_hgvs(input).unwrap();

    // Phase = Trans, arm 1 is position-specific Unknown (uncertainty
    // wrapper is preserved on the loc_edit, not the inner edit).
    let allele = match &parsed {
        HgvsVariant::Allele(a) => a,
        other => panic!("{input:?}: expected HgvsVariant::Allele, got {other:?}"),
    };
    assert_eq!(
        allele.phase,
        AllelePhase::Trans,
        "{input:?}: phase mismatch"
    );

    let arm1 = match &allele.variants[1] {
        HgvsVariant::Protein(v) => v,
        other => panic!("{input:?}: arm 1 expected Protein, got {other:?}"),
    };
    match arm1.loc_edit.edit.inner() {
        Some(ProteinEdit::Unknown {
            whole_protein: false,
            ..
        }) => {} // expected
        other => panic!(
            "{input:?}: arm 1 inner edit expected Unknown(whole_protein=false), got {other:?}"
        ),
    }

    // Display compacts to ACC:p.[Arg97Trp];[(Met1?)]. The bracket parsers
    // now dispatch on `(`, so the compact form re-parses cleanly (#544).
    let displayed = format!("{}", parsed);
    assert!(
        displayed.starts_with("NP_000079.2:p.["),
        "{input:?}: expected compact same-accession Display, got {displayed:?}"
    );
    assert!(
        displayed.contains("[(Met1?)]"),
        "{input:?}: expected predicted-form arm in compact Display, got {displayed:?}"
    );

    let reparsed = parse_hgvs(&displayed).unwrap_or_else(|e| {
        panic!("{input:?}: compact-form Display {displayed:?} must re-parse: {e}")
    });
    assert_eq!(
        reparsed, parsed,
        "{input:?}: re-parsing the compact Display must yield the same allele"
    );
}
