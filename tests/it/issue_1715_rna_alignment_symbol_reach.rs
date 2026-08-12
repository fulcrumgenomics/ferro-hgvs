//! Issue #1715 — the alignment-only-symbol refusal on the `r.` axis, measured.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/1715>. Sibling
//! of `issue_1627_named_element_alphabet_reach.rs`, which does the same job for
//! the DNA axes and should be read first.
//!
//! # The two halves of this file rest on different things — read the split
//!
//! [`the_grammar_reach_on_the_rna_axis`] and its two neighbours below it
//! measure **which spellings reach the AST at all** on `r.`. That is a property
//! of `parse_inserted_sequence`'s arms, and it is true independently of any
//! adjudication — it would stand even had the ruling gone the other way.
//!
//! Everything under **the refusal** asserts behaviour that rests on
//! `rulings[rna-axis-alignment-only-symbol-reach]`, ruled 2026-08-12: refuse,
//! with `standards.md:47`–`:61` governing. That record quotes both symbol
//! tables and both `general.md` bullets, and states the reading it was weighed
//! against. Read it before arguing the point from the spec text.
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
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

/// A single-exon transcript with `cds_start = 1`, so `r.N == c.N == n.N` and the
/// axis is the only thing that differs between the probes below.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_U.1".to_string(),
        Some("UG".to_string()),
        Strand::Plus,
        "AATTTGCCAATTTGCCAATTTGCC".to_string(),
        Some(1u64),
        Some(24u64),
        vec![Exon::new(1, 1, 24)],
        None,
        None,
        None,
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

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

/// Normalize `input` in `config`'s mode, from parse through normalize, as a
/// caller would. `Err` carries whichever stage refused.
fn round_trip(input: &str, config: ErrorConfig) -> Result<String, String> {
    let parsed = parse_hgvs_with_config(input, config).map_err(|e| e.to_string())?;
    Normalizer::new(provider())
        .normalize(&parsed.result)
        .map(|v| v.to_string())
        .map_err(|e| e.to_string())
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

// ---------------------------------------------------------------------------
// The refusal — asserts the ruled behaviour
// ---------------------------------------------------------------------------

/// **Rests on `rulings[rna-axis-alignment-only-symbol-reach]`**, ruled
/// 2026-08-12: `standards.md`'s RNA table (`:47`–`:61`) plus `general.md:50`
/// do reach a lower-case `x`, so an `r.` description stating one is refused.
///
/// The stage schedule is **not** this record's to choose — it is the decided
/// `rulings[absolute-prohibition-enforcement-stage]`, and it is asserted here in
/// the same shape `corpus_prohibited_inputs::`
/// `an_alignment_only_symbol_is_refused_in_every_mode_for_both_x_and_dash`
/// asserts it for the DNA axes:
///
/// | mode | parse | normalize |
/// |---|---|---|
/// | strict | **fails**, `W3028` | (unreached) |
/// | lenient | accepts, warns `W3028` | **fails** |
/// | silent | accepts, quiet | **fails** |
#[test]
fn a_lowercase_masked_nucleotide_on_the_rna_axis_is_refused_at_the_ruled_stage() {
    for input in REACHING_SHAPES {
        // STRICT — refused at PARSE. Strict validates input conformance, not
        // merely parseability.
        let err = parse_hgvs_with_config(input, ErrorConfig::strict())
            .map(|p| {
                panic!(
                    "{input}: strict must refuse at parse, got {p:?}",
                    p = p.result
                )
            })
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("W3028"),
            "{input}: strict must refuse with W3028; got: {err}"
        );
        // The diagnostic must send an `r.` author to the RNA clause and the RNA
        // spelling of the unknown base. Telling them to use `general.md:48`'s
        // CAPITALS, or `N`, would prescribe the one thing their axis forbids.
        assert!(
            err.contains("general.md:50") && err.contains("`n` is the IUPAC symbol"),
            "{input}: the refusal must cite the RNA ground; got: {err}"
        );

        // LENIENT and SILENT — accepted at parse, because neither validates
        // input conformance.
        for config in [ErrorConfig::lenient(), ErrorConfig::silent()] {
            assert!(
                parse_hgvs_with_config(input, config).is_ok(),
                "{input}: lenient/silent do not validate input conformance"
            );
        }

        // The parse-stage difference between the two is the message, not the
        // verdict.
        let lenient = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .expect("lenient accepts at parse");
        assert!(
            lenient
                .warnings
                .iter()
                .any(|w| w.error_type.code() == "W3028"),
            "{input}: lenient must say why it will fail later"
        );
        let silent =
            parse_hgvs_with_config(input, ErrorConfig::silent()).expect("silent accepts at parse");
        assert!(
            !silent
                .warnings
                .iter()
                .any(|w| w.error_type.code() == "W3028"),
            "{input}: silent is lenient without messages"
        );

        // All three modes refuse to NORMALIZE. Output conformance is rule 1 of
        // the README ruleset and has no mode escape, so this rung is not
        // mode-gated — lenient fails on the ruling's own ground, because a
        // masked base names no nucleotide to normalize.
        for (label, config) in [
            ("strict", ErrorConfig::strict()),
            ("lenient", ErrorConfig::lenient()),
            ("silent", ErrorConfig::silent()),
        ] {
            let err = round_trip(input, config)
                .expect_err("every mode must refuse an alignment-only symbol");
            assert!(
                err.contains("W3028"),
                "{input}: the {label} refusal must name W3028; got: {err}"
            );
        }
    }
}

/// A member of an `r.` composite is reached on **its own** axis, so the symbol
/// cannot hide behind a bracketed spelling. All three shapes, because "the cis
/// case works" does not entail the other two: a trans pair nests one `Allele`
/// inside another, and an uncertain group sets `uncertain` on it.
///
/// This is also what fails if the axis is ever resolved from the enclosing group
/// rather than from the member.
#[test]
fn a_member_of_an_rna_composite_is_refused_too() {
    for input in [
        "NM_U.1:r.[10delinsacgu;20delinsacgux]",
        "NM_U.1:r.[10delinsacgux];[20del]",
        "NM_U.1:r.(10delinsacgux)",
    ] {
        let err = parse_hgvs_with_config(input, ErrorConfig::strict())
            .map(|_| panic!("{input}: strict must refuse"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("W3028"), "{input}: got {err}");
        let err = round_trip(input, ErrorConfig::lenient())
            .expect_err("lenient cannot normalize a masked base");
        assert!(err.contains("W3028"), "{input}: got {err}");
    }
}

/// **The scope line, and the constraint on any widening.** The check is keyed on
/// the *symbol* `x`, never on the RNA alphabet as a closed set — a
/// closed-alphabet rule would refuse every name below, and whether a named
/// mobile element is conformant at all is the separate question
/// `DNA/complex.md:169` raises and no record settles.
///
/// `alu` is the load-bearing row: it is lower case **and** reaches `Named`, so
/// it is exactly what an alphabet rule catches and a symbol rule must not.
#[test]
fn a_genuine_named_element_on_the_rna_axis_is_untouched() {
    for input in [
        "NM_U.1:r.10delinsAluYb8",
        "NM_U.1:r.10delinsLINE1",
        "NM_U.1:r.10delinsalu",
        "NM_U.1:r.10delinsAlu",
    ] {
        assert!(
            is_named(input),
            "{input}: this is the shape the named-element arm exists for"
        );
        for config in [
            ErrorConfig::strict(),
            ErrorConfig::lenient(),
            ErrorConfig::silent(),
        ] {
            assert_eq!(
                round_trip(input, config).as_deref(),
                Ok(input),
                "{input}: a name carrying no daggered symbol is re-emitted unchanged"
            );
        }
    }
}

/// **The DNA axes do not move**, which is what makes this an `r.`-axis change
/// rather than a global widening to lower case. `n.` is on this side of the line
/// — it addresses a non-coding *DNA* reference sequence, so `general.md:48`'s
/// CAPITALS governs it and a lower-case letter is an ordinary name character.
///
/// This is the test that fails if the `Rna` arm is ever implemented as an
/// unconditional lower-case widening.
#[test]
fn a_lowercase_x_on_every_dna_axis_is_still_accepted() {
    for input in [
        "NM_U.1:c.10delinsACGTx",
        "NM_U.1:n.10delinsACGTx",
        "NC_TEST.1:g.10delinsACGTx",
    ] {
        assert!(
            parse_hgvs_with_config(input, ErrorConfig::strict()).is_ok(),
            "{input}: `general.md:48` puts a DNA description in CAPITALS, so a lower-case \
             letter is an ordinary name character — #1684 declined to rule on it and this \
             change does not either"
        );
    }
}

/// The uppercase spelling on the `r.` axis was already refused before #1715 and
/// **its diagnostic must not move**: it is refused on the DNA table's own
/// clause, which holds whichever way the RNA question is ruled.
#[test]
fn an_uppercase_masked_nucleotide_on_the_rna_axis_keeps_its_dna_diagnostic() {
    let err = parse_hgvs_with_config("NM_U.1:r.10delinsX", ErrorConfig::strict())
        .map(|_| panic!("strict already refused this before #1715"))
        .unwrap_err()
        .to_string();
    assert!(err.contains("W3028"), "{err}");
    assert!(
        err.contains("general.md:48") && err.contains("`N` is the IUPAC symbol"),
        "the uppercase spelling is refused on the DNA clause and must keep citing it; got: {err}"
    );
}
