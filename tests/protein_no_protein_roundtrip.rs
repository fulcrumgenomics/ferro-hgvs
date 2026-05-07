//! Round-trip pinning for protein "no product" notation: `p.0` and `p.0?` (#81 D6).
//!
//! HGVS spec context (assets/hgvs-nomenclature/docs/recommendations/protein/):
//!
//! - `p.0` — no protein is produced (definite), used for null alleles
//!   (substitution.md: "as a consequence of a variant in the translation
//!   initiation codon no protein is produced").
//! - `p.0?` — predicted no-protein (deletion.md: "Possibly no transcript is
//!   generated and no protein (`p.0?`)").
//!
//! These are distinct from `p.?` (unknown effect — D7) and `p.=` (silent — D5).
//! Section D of issue #81 already has parse-only coverage for these forms; this
//! file pins the full parse / Display / re-parse / normalize round-trip so that
//! any future change touching the protein edit pipeline must keep the certain
//! vs. predicted distinction intact.
//!
//! For each canonical input, we assert:
//!
//! 1. `parse_hgvs(s)` succeeds and produces a `Protein` variant.
//! 2. `format!("{}", parsed) == s` (Display is the spec form).
//! 3. The parsed edit is `ProteinEdit::NoProtein { predicted }` with the
//!    expected `predicted` flag (so `p.0` and `p.0?` are not collapsed).
//! 4. `parse_hgvs(format!("{}", parsed)) == parsed` (re-parse stability,
//!    structural equality).
//! 5. `normalize(parsed)` is a no-op on these whole-protein null forms and
//!    is idempotent.

use ferro_hgvs::hgvs::edit::ProteinEdit;
use ferro_hgvs::hgvs::variant::ProteinVariant;
use ferro_hgvs::{parse_hgvs, parse_hgvs_fast, HgvsVariant, MockProvider, Normalizer};

/// One canonical no-protein input plus the `predicted` flag we expect after parsing.
struct NoProteinCase {
    input: &'static str,
    predicted: bool,
}

/// Canonical inputs covering the cross product of {certain, predicted} and
/// {`NP_` accession, `LRG_` accession}. Both accession styles appear verbatim
/// in the HGVS spec for `p.0` (substitution.md uses `LRG_199p1:p.0`,
/// deletion.md uses `NP_003997.2:p.0?`).
const CASES: &[NoProteinCase] = &[
    NoProteinCase {
        input: "NP_003997.2:p.0",
        predicted: false,
    },
    NoProteinCase {
        input: "NP_003997.2:p.0?",
        predicted: true,
    },
    NoProteinCase {
        input: "LRG_199p1:p.0",
        predicted: false,
    },
    NoProteinCase {
        input: "LRG_199p1:p.0?",
        predicted: true,
    },
];

/// Borrow the inner `ProteinVariant` from an `HgvsVariant`, panicking with a
/// useful message for any other arm. We only call this on values we have
/// already round-tripped through the protein parse path.
fn protein_variant(variant: &HgvsVariant) -> &ProteinVariant {
    match variant {
        HgvsVariant::Protein(p) => p,
        other => panic!("expected Protein variant, got {:?}", other),
    }
}

/// Borrow the inner `ProteinEdit`, panicking if the loc_edit somehow lacks an
/// inner edit (which would itself be a regression for `p.0` / `p.0?`).
fn protein_edit(variant: &HgvsVariant) -> &ProteinEdit {
    protein_variant(variant)
        .loc_edit
        .edit
        .inner()
        .expect("p.0 / p.0? should always carry a concrete (non-uncertain) edit")
}

#[test]
fn parses_displays_and_classifies_no_protein() {
    for case in CASES {
        let parsed = parse_hgvs(case.input)
            .unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", case.input, e));

        // (1) The result is a Protein variant.
        assert!(
            matches!(parsed, HgvsVariant::Protein(_)),
            "{:?} did not parse as a Protein variant: {:?}",
            case.input,
            parsed
        );

        // (2) Display is the spec form.
        let displayed = format!("{}", parsed);
        assert_eq!(
            displayed, case.input,
            "Display round-trip drift for {:?}",
            case.input
        );

        // (3) Edit is NoProtein with the expected predicted flag — pins the
        // certain (`p.0`) vs predicted (`p.0?`) distinction.
        match protein_edit(&parsed) {
            ProteinEdit::NoProtein { predicted } => assert_eq!(
                *predicted, case.predicted,
                "predicted flag mismatch for {:?}: expected {}, got {}",
                case.input, case.predicted, *predicted
            ),
            other => panic!(
                "{:?} did not parse as ProteinEdit::NoProtein, got {:?}",
                case.input, other
            ),
        }
    }
}

#[test]
fn reparse_is_structurally_stable() {
    for case in CASES {
        let parsed = parse_hgvs(case.input)
            .unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", case.input, e));

        // (4) parse(Display(parsed)) == parsed structurally — and Display is
        // stable across the round-trip.
        let displayed = format!("{}", parsed);
        let reparsed = parse_hgvs(&displayed)
            .unwrap_or_else(|e| panic!("re-parse of {:?} failed: {}", displayed, e));

        assert_eq!(
            parsed, reparsed,
            "structural equality lost on re-parse of {:?}",
            case.input
        );
        assert_eq!(
            format!("{}", reparsed),
            displayed,
            "Display drifted on re-parse of {:?}",
            case.input
        );
    }
}

#[test]
fn normalize_is_noop_and_idempotent() {
    // No-protein variants describe the absence of a product and have no
    // genomic coordinates to shift, so normalize() must be a structural no-op.
    // We use a bare MockProvider (no transcript data) — the normalizer's
    // protein path validates against a provider only when one has protein data,
    // so this exercises the no-op pass-through path.
    let provider = MockProvider::new();
    let normalizer = Normalizer::new(provider);

    for case in CASES {
        let parsed = parse_hgvs(case.input)
            .unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", case.input, e));

        // (5a) normalize is a no-op on p.0 / p.0?.
        let normalized = normalizer
            .normalize(&parsed)
            .unwrap_or_else(|e| panic!("normalize({:?}) failed: {}", case.input, e));
        assert_eq!(
            parsed, normalized,
            "normalize unexpectedly altered {:?}",
            case.input
        );
        assert_eq!(
            format!("{}", normalized),
            case.input,
            "normalize Display drift for {:?}",
            case.input
        );

        // (5b) normalize(normalize(x)) == normalize(x).
        let renormalized = normalizer
            .normalize(&normalized)
            .unwrap_or_else(|e| panic!("re-normalize({:?}) failed: {}", case.input, e));
        assert_eq!(
            normalized, renormalized,
            "normalize is not idempotent on {:?}",
            case.input
        );

        // And the post-normalize re-parse path stays stable too.
        let reparsed_normalized = parse_hgvs(&format!("{}", normalized))
            .unwrap_or_else(|e| panic!("re-parse of normalize({:?}) failed: {}", case.input, e));
        assert_eq!(
            normalized, reparsed_normalized,
            "structural equality lost on re-parse of normalize({:?})",
            case.input
        );
    }
}

#[test]
fn predicted_and_certain_are_distinct_variants() {
    // Guards against any future "simplification" that drops the predicted bit
    // and silently collapses p.0? into p.0.
    let certain = parse_hgvs("NP_003997.2:p.0").expect("parse NP:p.0");
    let predicted = parse_hgvs("NP_003997.2:p.0?").expect("parse NP:p.0?");

    assert_ne!(
        certain, predicted,
        "p.0 and p.0? must remain structurally distinct"
    );
    assert_ne!(
        format!("{}", certain),
        format!("{}", predicted),
        "p.0 and p.0? must render distinctly"
    );

    assert!(matches!(
        protein_edit(&certain),
        ProteinEdit::NoProtein { predicted: false }
    ));
    assert!(matches!(
        protein_edit(&predicted),
        ProteinEdit::NoProtein { predicted: true }
    ));
}

/// Pin the internal shape of `p.0` / `p.0?`:
///
/// - The outer `Mu` wrapper is `Mu::Certain` for *both* the certain (`p.0`) and
///   the predicted (`p.0?`) form. The `?` is part of the inner `ProteinEdit`,
///   not an outer `Mu::Uncertain`. A refactor that promoted `p.0?` to
///   `Mu::Uncertain` while clearing the inner `predicted` flag would Display
///   the same string but break consumers that pattern-match on `Mu::Uncertain`
///   to detect predicted variants.
/// - The dummy location is the canonical `ProtPos::new(AminoAcid::Met, 1)`
///   point interval, matching what the parser builds at
///   src/hgvs/parser/variant.rs:1735-1736 (and the duplicate path at :2909-2910).
///   This is an internal detail today, but pinning it stops a refactor from
///   silently changing the dummy in a way that breaks structural-equality
///   consumers (caches, sets, JSON serializers).
#[test]
fn internal_shape_is_mu_certain_with_dummy_met1() {
    use ferro_hgvs::hgvs::interval::ProtInterval;
    use ferro_hgvs::hgvs::location::{AminoAcid, ProtPos};
    use ferro_hgvs::hgvs::uncertainty::Mu;

    for case in CASES {
        let parsed = parse_hgvs(case.input)
            .unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", case.input, e));
        let pv = protein_variant(&parsed);

        // Outer Mu wrapper must be Certain for both forms.
        match &pv.loc_edit.edit {
            Mu::Certain(_) => {}
            Mu::Uncertain(_) => panic!(
                "{:?}: edit Mu wrapper must be Certain (the `?` is part of the inner edit, \
                 not an outer Mu::Uncertain)",
                case.input
            ),
            other => panic!(
                "{:?}: edit Mu wrapper must be Certain, got {:?}",
                case.input, other
            ),
        }

        // Dummy location is the canonical Met1 point interval, exactly as the
        // parser builds it.
        let expected_pos = ProtPos::new(AminoAcid::Met, 1);
        let expected_interval = ProtInterval::point(expected_pos);
        assert_eq!(
            pv.loc_edit.location, expected_interval,
            "{:?}: expected dummy Met1 point interval, got {:?}",
            case.input, pv.loc_edit.location
        );
    }
}

/// Multi-pass normalize idempotency + Hash/Eq agreement.
///
/// Three passes pin "stable fixed point under repeated normalization", which
/// is what the normalizer's caches and any downstream dedup tables rely on.
/// We also assert hash equality across passes — `Hash` is derived on
/// `HgvsVariant` and `ProteinEdit`, but a future change to the derive (e.g.
/// adding a non-hash field) would not be caught by the `Eq` assertions alone.
#[test]
fn normalize_is_idempotent_across_multiple_passes_with_hash_agreement() {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    fn hash_of(v: &HgvsVariant) -> u64 {
        let mut h = DefaultHasher::new();
        v.hash(&mut h);
        h.finish()
    }

    let provider = MockProvider::new();
    let normalizer = Normalizer::new(provider);

    for case in CASES {
        let parsed = parse_hgvs(case.input)
            .unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", case.input, e));

        let n1 = normalizer.normalize(&parsed).expect("pass 1");
        let n2 = normalizer.normalize(&n1).expect("pass 2");
        let n3 = normalizer.normalize(&n2).expect("pass 3");

        assert_eq!(parsed, n1, "{:?}: pass 1 must be a no-op", case.input);
        assert_eq!(n1, n2, "{:?}: pass 2 must be idempotent", case.input);
        assert_eq!(n2, n3, "{:?}: pass 3 must be idempotent", case.input);

        // Hash agrees with Eq across every pair.
        let h0 = hash_of(&parsed);
        assert_eq!(
            h0,
            hash_of(&n1),
            "{:?}: hash drift after pass 1",
            case.input
        );
        assert_eq!(
            h0,
            hash_of(&n2),
            "{:?}: hash drift after pass 2",
            case.input
        );
        assert_eq!(
            h0,
            hash_of(&n3),
            "{:?}: hash drift after pass 3",
            case.input
        );
    }
}

/// `Hash` must distinguish `p.0` from `p.0?`. PR #130's
/// `predicted_and_certain_are_distinct_variants` already pins `Eq`; pin `Hash`
/// independently because a derive change that hashed only the enum
/// discriminant could silently collide them.
#[test]
fn certain_and_predicted_have_distinct_hashes() {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    fn hash_of(v: &HgvsVariant) -> u64 {
        let mut h = DefaultHasher::new();
        v.hash(&mut h);
        h.finish()
    }

    let certain = parse_hgvs("NP_003997.2:p.0").expect("parse NP:p.0");
    let predicted = parse_hgvs("NP_003997.2:p.0?").expect("parse NP:p.0?");

    assert_ne!(
        hash_of(&certain),
        hash_of(&predicted),
        "p.0 and p.0? must hash distinctly so caches/dedup tables don't collide them"
    );

    // And two parses of the same input hash identically (Hash agrees with Eq).
    let certain2 = parse_hgvs("NP_003997.2:p.0").expect("parse NP:p.0");
    assert_eq!(certain, certain2);
    assert_eq!(hash_of(&certain), hash_of(&certain2));
}

/// `p.0` and `p.0?` (D6) must be structurally and Display-distinct from the
/// adjacent whole-protein forms in spec Section D:
///
/// - `p.=` (D5 — silent / no change), and the predicted form `p.(=)`
/// - `p.?` (D7 — unknown effect),     and the predicted form `p.(?)`
///
/// Without this guard, a future "simplification" that folds `NoProtein` into
/// either neighbour would Display-roundtrip its own form correctly but pass
/// PR #130's "p.0 vs p.0? are distinct" assertion. We pin all six pairwise
/// inequalities here so D6 stands on its own.
#[test]
fn p0_is_distinct_from_neighbouring_whole_protein_forms() {
    let p0 = parse_hgvs("NP_003997.2:p.0").expect("parse p.0");
    let p0q = parse_hgvs("NP_003997.2:p.0?").expect("parse p.0?");
    let p_eq = parse_hgvs("NP_003997.2:p.=").expect("parse p.=");
    let p_eq_pred = parse_hgvs("NP_003997.2:p.(=)").expect("parse p.(=)");
    let p_q = parse_hgvs("NP_003997.2:p.?").expect("parse p.?");
    let p_q_pred = parse_hgvs("NP_003997.2:p.(?)").expect("parse p.(?)");

    let neighbours: &[(&str, &HgvsVariant)] = &[
        ("p.=", &p_eq),
        ("p.(=)", &p_eq_pred),
        ("p.?", &p_q),
        ("p.(?)", &p_q_pred),
    ];

    for (label, n) in neighbours {
        assert_ne!(&p0, *n, "p.0 must be structurally distinct from {}", label);
        assert_ne!(
            &p0q, *n,
            "p.0? must be structurally distinct from {}",
            label
        );
        assert_ne!(
            format!("{}", p0),
            format!("{}", *n),
            "p.0 must Display distinctly from {}",
            label
        );
        assert_ne!(
            format!("{}", p0q),
            format!("{}", *n),
            "p.0? must Display distinctly from {}",
            label
        );
    }
}

/// Malformed inputs adjacent to `p.0` / `p.0?` must reject. The HGVS spec
/// (assets/hgvs-nomenclature/docs/recommendations/protein/) defines exactly
/// two whole-protein no-product forms: `p.0` (certain) and `p.0?` (predicted).
/// Anything else in this neighbourhood is malformed:
///
/// - Trailing junk: `p.0X`, `p.0a`, `p.0!`, `p.00`.
/// - Parenthesised: `p.(0)`, `p.(0?)` — spec uses `?` suffix, never parens.
/// - Internal whitespace: `p. 0`.
/// - No accession: `p.0`, `p.0?` — `parse_hgvs` always requires accession.
/// - Cis-bracket form: `p.[0]`, `p.[Arg97Trp;0]`, `p.[0;Arg97Trp]` — `alleles.md`
///   only allows `[0]` in trans, never cis.
#[test]
fn malformed_adjacent_inputs_reject() {
    let bad: &[&str] = &[
        // Trailing junk after the canonical 0 / 0?.
        "NP_003997.2:p.0X",
        "NP_003997.2:p.0a",
        "NP_003997.2:p.0!",
        "NP_003997.2:p.00",
        // Parenthesised forms — not in spec.
        "NP_003997.2:p.(0)",
        "NP_003997.2:p.(0?)",
        // Interior whitespace — outer trim is fine, but the parser must reject
        // `p. 0`.
        "NP_003997.2:p. 0",
        "NP_003997.2:p. 0?",
        // Bare (no accession): parse_hgvs always requires accession+coordinate.
        "p.0",
        "p.0?",
        // Cis-bracket: spec allows [0] only in TRANS allele context.
        "NP_003997.2:p.[0]",
        "NP_003997.2:p.[Arg97Trp;0]",
        "NP_003997.2:p.[0;Arg97Trp]",
    ];

    for input in bad {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "expected parse_hgvs({:?}) to reject as malformed, but it accepted: {:?}",
            input,
            result.ok()
        );
    }
}

/// `parse_hgvs` trims leading/trailing whitespace on the input. Pin this
/// for D6 inputs so the outer-trim contract can't silently regress for
/// `p.0` / `p.0?`. (Interior whitespace, e.g. `p. 0`, is covered by the
/// malformed-rejection test.)
#[test]
fn outer_whitespace_is_trimmed() {
    for case in CASES {
        let padded = format!("   {}   ", case.input);
        let parsed = parse_hgvs(&padded)
            .unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", padded, e));
        // Display matches the un-padded canonical form.
        assert_eq!(
            format!("{}", parsed),
            case.input,
            "outer-trim Display drift for {:?}",
            padded
        );
        // ...and matches what the un-padded parse produces, structurally.
        let unpadded = parse_hgvs(case.input).expect("unpadded parse");
        assert_eq!(
            parsed, unpadded,
            "outer-trim structural drift for {:?}",
            padded
        );
    }
}

/// `parse_hgvs_fast` falls back to the standard parser for protein inputs and
/// must produce a structurally identical result on D6 inputs. Pinning this
/// keeps the fast-path entry honest as it grows.
#[test]
fn fast_path_parity_for_no_protein_inputs() {
    for case in CASES {
        let std = parse_hgvs(case.input)
            .unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", case.input, e));
        let fast = parse_hgvs_fast(case.input)
            .unwrap_or_else(|e| panic!("parse_hgvs_fast({:?}) failed: {}", case.input, e));
        assert_eq!(
            std, fast,
            "{:?}: fast-path produced structurally different result",
            case.input
        );
        assert_eq!(
            format!("{}", std),
            format!("{}", fast),
            "{:?}: fast-path Display drift",
            case.input
        );
    }
}

/// `ProteinEdit::NoProtein` must survive being placed inside a trans-allele
/// compound. The expanded form `[ACC:p.X];[ACC:p.0(?)]` parses each arm
/// through the protein dispatch, so the second arm carries
/// `ProteinEdit::NoProtein { predicted }` (with the right flag) and the outer
/// allele is `phase: Trans`. This is the spec-blessed "p.0 inside an allele"
/// pattern from `alleles.md`.
///
/// We use the expanded form (full prefix on each arm) because that's the
/// shape PR #130's D6 contract owns: each `p.0(?)` element goes through the
/// no-protein parse path. The compact-form Display→reparse round-trip
/// (`ACC:p.[X];[0]`) is governed by the bare-bracket trans dispatch, which
/// is out of D6 scope.
#[test]
fn trans_allele_carries_no_protein_edit_in_expanded_form() {
    use ferro_hgvs::hgvs::variant::AllelePhase;

    let cases: &[(&str, bool)] = &[
        ("[NP_003997.2:p.Arg97Trp];[NP_003997.2:p.0]", false),
        ("[NP_003997.2:p.Arg97Trp];[NP_003997.2:p.0?]", true),
    ];

    for (input, expected_predicted) in cases {
        let parsed =
            parse_hgvs(input).unwrap_or_else(|e| panic!("parse_hgvs({:?}) failed: {}", input, e));

        let allele = match &parsed {
            HgvsVariant::Allele(a) => a,
            other => panic!("{:?}: expected Allele variant, got {:?}", input, other),
        };

        assert_eq!(
            allele.phase,
            AllelePhase::Trans,
            "{:?}: expected Trans phase, got {:?}",
            input,
            allele.phase
        );
        assert_eq!(
            allele.variants.len(),
            2,
            "{:?}: expected exactly 2 sub-variants, got {}",
            input,
            allele.variants.len()
        );

        // Second arm carries ProteinEdit::NoProtein with the right `predicted`.
        let second = &allele.variants[1];
        match protein_edit(second) {
            ProteinEdit::NoProtein { predicted } => assert_eq!(
                *predicted, *expected_predicted,
                "{:?}: predicted flag mismatch on second arm: expected {}, got {}",
                input, expected_predicted, *predicted
            ),
            other => panic!(
                "{:?}: second arm must be ProteinEdit::NoProtein, got {:?}",
                input, other
            ),
        }
    }
}

/// Spec-critical D6 distinction: bare `[0]` in a trans allele is the
/// cross-coordinate "absent allele" marker (`HgvsVariant::NullAllele`),
/// **not** a protein-coordinate `p.0` (no-protein-produced) edit. They are
/// semantically distinct:
///
/// - `[ACC:p.Ser86Arg];[0]` — second X-chromosome is *absent*
///   (alleles.md line 78).
/// - `[ACC:p.Ser86Arg];[ACC:p.0]` — second allele *exists* but produces no
///   protein due to a start-codon variant (substitution.md line 42).
///
/// PR #130 doesn't pin this distinction. A refactor that folds `[0]` into
/// `ProteinEdit::NoProtein` would Display-roundtrip identically but corrupt
/// downstream interpretation (a hemizygous male becomes "no protein from this
/// allele due to a Met1 variant", which is biologically wrong).
#[test]
fn bare_bracketed_zero_is_null_allele_not_no_protein() {
    use ferro_hgvs::hgvs::variant::AllelePhase;

    // (a) Bare [0] -> NullAllele.
    let bare = parse_hgvs("[LRG_199p1:p.Ser86Arg];[0]").expect("parse bare-[0] trans");
    let bare_allele = match &bare {
        HgvsVariant::Allele(a) => a,
        other => panic!("expected Allele, got {:?}", other),
    };
    assert_eq!(bare_allele.phase, AllelePhase::Trans);
    assert_eq!(bare_allele.variants.len(), 2);
    assert!(
        matches!(bare_allele.variants[1], HgvsVariant::NullAllele),
        "bare [0] arm must be HgvsVariant::NullAllele, got {:?}",
        bare_allele.variants[1]
    );

    // (b) Prefixed [ACC:p.0] -> Protein(... NoProtein).
    let prefixed =
        parse_hgvs("[LRG_199p1:p.Ser86Arg];[LRG_199p1:p.0]").expect("parse prefixed trans");
    let prefixed_allele = match &prefixed {
        HgvsVariant::Allele(a) => a,
        other => panic!("expected Allele, got {:?}", other),
    };
    assert_eq!(prefixed_allele.phase, AllelePhase::Trans);
    assert_eq!(prefixed_allele.variants.len(), 2);
    assert!(
        matches!(&prefixed_allele.variants[1], HgvsVariant::Protein(_)),
        "prefixed [ACC:p.0] arm must be HgvsVariant::Protein, got {:?}",
        prefixed_allele.variants[1]
    );
    match protein_edit(&prefixed_allele.variants[1]) {
        ProteinEdit::NoProtein { predicted: false } => {}
        other => panic!(
            "prefixed [ACC:p.0] arm must carry ProteinEdit::NoProtein {{ predicted: false }}, got {:?}",
            other
        ),
    }

    // (c) The two are structurally NOT equal.
    assert_ne!(
        bare, prefixed,
        "bare [0] and prefixed [ACC:p.0] must remain structurally distinct"
    );
}
