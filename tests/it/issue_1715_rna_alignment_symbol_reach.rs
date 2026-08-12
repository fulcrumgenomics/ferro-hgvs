//! Issue #1715 — the alignment-only-symbol refusal on the `r.` axis, measured.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1715>. Sibling
//! of `issue_1627_named_element_alphabet_reach.rs`, which does the same job for
//! the DNA axes and should be read first.
//!
//! # This file is the MEASUREMENT half, and it is settled
//!
//! Every test here measures **which spellings reach the AST at all** on `r.`.
//! That is a property of `parse_inserted_sequence`'s arms, so it is true
//! whichever way the open question below is answered, and it is the half worth
//! keeping either way. Whether ferro should *refuse* the shapes it finds is an
//! adjudication nobody has made; it is proposed separately, so that backing the
//! proposal out does not take the measurement with it.
//!
//! # What the question is, and what it is not
//!
//! It is **not** whether an uppercase `X` may appear in a description — the
//! decided `rulings[alignment-only-symbol-in-a-description]` settles that, and
//! it already reaches `r.10delinsX`. It is **not** at which stage a refusal
//! lives — the decided `rulings[absolute-prohibition-enforcement-stage]` settles
//! that: strict fails at PARSE, lenient does not validate input conformance and
//! fails only when it cannot NORMALIZE, silent is lenient without messages.
//!
//! It is whether the **lower-case** spelling on the `r.` axis is the same
//! finding. `general.md:48` is the DNA-level bullet ("nucleotides in CAPITALS");
//! `general.md:50` is its RNA sibling ("nucleotides in lower case"). So the
//! premise #1684 declined on — "a lower-case letter is indistinguishable from an
//! ordinary name character" — inverts on `r.`, where lower case *is* the
//! sequence alphabet.
//!
//! # The reach on `r.` is narrower than on the DNA axes, and that matters twice
//!
//! `parse_inserted_sequence` dispatches on the **first** byte, and no arm
//! accepts a leading lower-case non-IUPAC letter. So a **lone** or **leading**
//! `x` never reaches the AST: the grammar refuses it in every mode, exactly as
//! it refuses `-`, and the two-stage mode schedule does **not** describe it.
//! Reading the schedule as covering every `x` would be the same misattribution
//! `issue_1627_named_element_alphabet_reach.rs` exists to prevent for `-`.
//!
//! It matters a second time for the *population*: on the DNA axes all four
//! offsets reach `Named`, here only the non-leading ones do, so a guard copied
//! across from that file without re-measuring would assert a refusal the grammar
//! makes and call it a rule.
//!
//! # The corpus figure is a STRUCTURAL zero and must not be quoted as evidence
//!
//! `RefShape::all()` enumerates `Genomic`, `CodingSingleExon`, `CodingMultiExon`
//! and `NonCodingMultiExon` — `g.`, `c.` and `n.`. There is no `r.` shape, so of
//! the 58,552 spellings the corpus builds, **0** are on this axis and **0**
//! contain a lower-case `x`. No census counter can move, in either direction,
//! and a `0 moved` from that harness is a claim about the generator rather than
//! about this change. Pinned by [`the_corpus_cannot_reach_the_rna_axis_at_all`]
//! so the zero is asserted as structural rather than read as reassurance.

use ferro_hgvs::conformance::spec_corpus::{corpus, CorpusBounds};
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::edit::{InsertedSequence, NaEdit};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

/// The `InsertedSequence` a single `r.` `delins`/`ins` carries, or `None` when
/// the input does not parse at all.
///
/// Reads through the bare `parse_hgvs`, which applies no `ErrorConfig` (#1632),
/// so it observes the **grammar** alone and is unaffected by whichever way the
/// ruling goes.
fn inserted(input: &str) -> Option<InsertedSequence> {
    let variant = parse_hgvs(input).ok()?;
    let HgvsVariant::Rna(rna) = variant else {
        panic!("{input}: expected an r. variant");
    };
    match rna.loc_edit.edit.inner() {
        Some(NaEdit::Insertion { sequence }) | Some(NaEdit::Delins { sequence, .. }) => {
            Some(sequence.clone())
        }
        other => panic!("{input}: expected an ins or delins, got {other:?}"),
    }
}

/// True when the insert was classified as a named mobile element.
fn is_named(input: &str) -> bool {
    matches!(inserted(input), Some(InsertedSequence::Named(_)))
}

/// The shapes on `r.` that carry a non-leading lower-case `x` into the AST.
/// Every offset the grammar admits is represented, including the `ins` spelling
/// and an uppercase-led name, because "the trailing case works" does not entail
/// the others — that asymmetry is exactly what `X` turned out to have.
const REACHING_SHAPES: &[&str] = &[
    "NM_U.1:r.10delinsacgux",  // stray at the end
    "NM_U.1:r.10delinsacxgu",  // stray in the middle
    "NM_U.1:r.10delinsax",     // shortest reaching form
    "NM_U.1:r.10delinsaxa",    // interior, short
    "NM_U.1:r.10delinsacguxx", // two strays
    "NM_U.1:r.10_11insacgux",  // the `ins` spelling
    "NM_U.1:r.10delinsAcgux",  // via the uppercase arm rather than the IUPAC one
];

/// The shapes the **grammar** refuses outright, in every mode. Named separately
/// from [`REACHING_SHAPES`] because the reason is different in kind.
const GRAMMAR_REFUSED_SHAPES: &[&str] = &[
    "NM_U.1:r.10delinsx",     // lone
    "NM_U.1:r.10delinsxa",    // leading, short
    "NM_U.1:r.10delinsxacgu", // leading
    "NM_U.1:r.10_11insx",     // lone, the `ins` spelling
    "NM_U.1:r.10delins-",     // `standards.md:37`, for contrast
    "NM_U.1:r.10delinsac-gu", // …and embedded
];

// ---------------------------------------------------------------------------
// The measurement — true whichever way the ruling goes
// ---------------------------------------------------------------------------

/// **The defect, at AST level.** One stray lower-case `x` reclassifies an
/// otherwise-literal RNA run *in its entirety* as a mobile-element name, through
/// the `is_iupac_base` arm of `parse_inserted_sequence` — the same wholesale
/// reclassification #1627 found for uppercase `X`, one axis over.
///
/// This is a grammar property and is deliberately asserted through the bare
/// `parse_hgvs`, so it neither depends on nor asserts the ruling.
#[test]
fn the_grammar_reach_on_the_rna_axis() {
    // The control: a pure RNA IUPAC run is `Literal`, so the reclassification
    // really is triggered by the stray byte and not by the arm always winning.
    assert!(
        matches!(
            inserted("NM_U.1:r.10delinsacgu"),
            Some(InsertedSequence::Literal(_))
        ),
        "a pure `acgu` run must be Literal, or every assertion below is vacuous"
    );

    for input in REACHING_SHAPES {
        assert!(
            is_named(input),
            "{input}: a non-leading `x` reclassifies the whole run as Named"
        );
    }
}

/// **A lone or leading `x` never reaches the AST**, so no conformance rule is
/// consulted for it at either stage — the grammar refuses it, in every mode, for
/// grammar reasons. That is the standing fact `-` has, and describing either as
/// mode-dependent is the misattribution this test exists to prevent.
///
/// Measured in every position because "refused when leading" would not have
/// entailed "refused when trailing", which is precisely the asymmetry `X` had.
#[test]
fn a_lone_or_leading_lowercase_x_is_refused_by_the_grammar_rather_than_by_a_rule() {
    for input in GRAMMAR_REFUSED_SHAPES {
        let err = parse_hgvs(input)
            .map(|v| panic!("{input} must be refused by the grammar, but parsed to {v:?}"))
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("trailing") || err.contains("Parse error"),
            "{input}: the refusal is a grammar one, not an alphabet finding; got: {err}"
        );
        for config in [
            ErrorConfig::strict(),
            ErrorConfig::lenient(),
            ErrorConfig::silent(),
        ] {
            assert!(
                parse_hgvs_with_config(input, config).is_err(),
                "{input}: a grammar refusal is not mode-gated"
            );
        }
    }
}

/// **The corpus cannot reach this axis, so its zero is structural.**
///
/// `RefShape::all()` builds `g.`, `c.` and `n.` references and nothing else, so
/// no row it can produce is on the `r.` axis or carries a lower-case `x`. Any
/// "0 rows move" from that harness is therefore a claim about the generator, and
/// the population for this change has to be established by construction — which
/// is what [`REACHING_SHAPES`] is.
///
/// Asserted rather than asserted-about: the denominator is checked non-zero too,
/// so a corpus that silently stopped building anything could not pass this as a
/// clean bill of health.
#[test]
fn the_corpus_cannot_reach_the_rna_axis_at_all() {
    let built = corpus(&CorpusBounds::default());
    let spellings: Vec<&String> = built.rows.iter().flat_map(|row| &row.spellings).collect();

    assert!(
        spellings.len() > 10_000,
        "the corpus built {} spellings — the zeroes below would be vacuous",
        spellings.len()
    );
    assert_eq!(
        spellings.iter().filter(|s| s.contains(":r.")).count(),
        0,
        "`RefShape::all()` has no `r.` shape; if that changes, re-measure this change's \
         movement instead of quoting the structural zero"
    );
    assert_eq!(
        spellings.iter().filter(|s| s.contains('x')).count(),
        0,
        "no corpus spelling states a lower-case `x` on any axis"
    );
}
