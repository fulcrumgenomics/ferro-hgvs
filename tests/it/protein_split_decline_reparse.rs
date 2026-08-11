//! Every decline of the protein-axis split move leaves a description
//! `parse_hgvs` accepts.
//!
//! # Why this guard exists
//!
//! `Normalizer::try_protein_split_delins` (added by #1606) turns an
//! equal-length `p.` delins whose interior is unchanged into its members, and
//! declines on seven distinct conditions. One of those declines — **residue 1
//! is the changed one** — exists *only* to stop normalization emitting
//! `p.Met1Xxx`, which `protein/substitution.md:49` forbids and `parse_hgvs`
//! refuses. So that guard's whole value is a re-parse property, and nothing was
//! asserting the property directly: #1606's own suite pins the exact output
//! string for each decline, which goes green just as happily on a string that
//! cannot be read back.
//!
//! The other six declines are covered here for the same reason, cheaply: a
//! future decline that starts emitting a rewritten form rather than the input
//! as authored would otherwise have nothing checking that the rewrite is
//! legal.
//!
//! # What this file deliberately does NOT cover
//!
//! Ferro *does* emit an unparseable `p.Met1Gly` for a protein delins whose
//! common-affix trimming leaves residue 1 as the residual — `p.Met1_Ser3delins`
//! `GlyAlaSer` is the smallest case. That is **#1607**, it lives in
//! `try_protein_delins_canonicalize`, and it is not the split: measured over a
//! 4,336-input sweep, every such output enters `try_protein_split_delins`
//! already collapsed to a `Substitution` and is turned away on the
//! "not a delins" gate, and disabling the split entirely reproduces the same 16
//! unparseable outputs byte for byte. Pinning it belongs to whatever closes
//! #1607, not here — a guard that failed today would make this file a
//! duplicate bug report rather than a regression guard for the split.
//!
//! Hermetic: no `FERRO_MANIFEST`, no network, no environment variable.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// The accession every case is written against.
const NP: &str = "NP_003997.1";

/// Residues 1..=20. Residue 1 is `Met`; 2, 3, 4 are `Ala`, `Ser`, `Leu`, which
/// is the span the splittable and declined cases below both name.
const PROTEIN: &str = "MASLWEKTRNDQHIPFYCAA";

/// The same protein with a stop **inserted** at position 4, for the
/// reference-side `Ter` decline. It is 21 residues to [`PROTEIN`]'s 20 —
/// nothing was replaced, so `Leu` and everything after it shift one position 3'
/// (`Leu` is residue 5 here, not 4).
const PROTEIN_WITH_STOP_AT_4: &str = "MAS*LWEKTRNDQHIPFYCAA";

/// The same protein with an unknown residue **inserted** at position 4, for the
/// reference-side `Xaa` decline. Same shape and same caveat as
/// [`PROTEIN_WITH_STOP_AT_4`]: 21 residues, and `Leu` is residue 5.
const PROTEIN_WITH_UNKNOWN_AT_4: &str = "MASXLWEKTRNDQHIPFYCAA";

fn normalizer_over(sequences: &[(&str, &str)]) -> Normalizer<MockProvider> {
    let mut provider = MockProvider::new();
    for (accession, sequence) in sequences {
        provider.add_protein(*accession, *sequence);
    }
    Normalizer::new(provider)
}

fn normalizer() -> Normalizer<MockProvider> {
    normalizer_over(&[(NP, PROTEIN)])
}

/// Normalize `input`, assert the result is exactly `expected`, and assert the
/// result **round-trips** through `parse_hgvs`.
///
/// The two assertions answer different questions and the second is the one this
/// file is about: an exact-string assertion alone cannot tell a legal output
/// from an illegal one.
///
/// The round-trip is a re-render comparison, not an `is_ok()`. `is_ok()` only
/// says the parser did not refuse the string — it would pass on an output the
/// parser silently read as some *other* description, which for the `p.` axis is
/// a live shape (a member list flattened into the wrong bracket parses fine and
/// means something else). Comparing `parse_hgvs(&output)?.to_string()` against
/// `output` asks whether the parser and the renderer agree on what was written.
fn assert_normalizes_and_reparses(
    normalizer: &Normalizer<MockProvider>,
    input: &str,
    expected: &str,
) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("input {input} must parse: {e}"));
    let output = normalizer
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalizing {input}: {e}"))
        .to_string();
    assert_eq!(output, expected, "normalized form of {input}");
    let reparsed = parse_hgvs(&output).unwrap_or_else(|e| {
        panic!("normalization of {input} emitted {output}, which parse_hgvs refuses: {e}")
    });
    assert_eq!(
        reparsed.to_string(),
        output,
        "normalization of {input} emitted {output}, which parse_hgvs reads back as a \
         different description"
    );
}

/// The positive control, and it is load-bearing: without it every assertion
/// below could pass on a build where the split move never runs at all, and this
/// file would report a structural zero as a clean bill of health.
///
/// Residues 2..=4 are `Ala-Ser-Leu`; the payload `Gly-Ser-Gly` leaves `Ser`
/// unchanged between two changed residues, which is the separation
/// `protein/delins.md:21` is about.
#[test]
fn the_split_move_is_reached_and_fires() {
    assert_normalizes_and_reparses(
        &normalizer(),
        &format!("{NP}:p.Ala2_Leu4delinsGlySerGly"),
        &format!("{NP}:p.[Ala2Gly;Leu4Gly]"),
    );
}

/// The predicted form of the same input, so the `( )` marker's placement on the
/// split's output is checked for legality too — it moves from the edit onto the
/// whole allele, which is a different string than any decline can produce.
#[test]
fn the_split_move_fires_on_a_predicted_description() {
    assert_normalizes_and_reparses(
        &normalizer(),
        &format!("{NP}:p.(Ala2_Leu4delinsGlySerGly)"),
        &format!("{NP}:p.[(Ala2Gly;Leu4Gly)]"),
    );
}

/// **Residue 1 changed.** The decline whose entire purpose is this property:
/// splitting would emit `p.[Met1Gly;Ser3Gly]`, and `p.Met1Gly` is a description
/// `parse_hgvs` refuses (`protein/substitution.md:49`).
#[test]
fn a_changed_residue_one_declines_to_a_description_that_reparses() {
    let normalizer = normalizer();
    for (input, expected) in [
        (
            format!("{NP}:p.Met1_Ser3delinsGlyAlaGly"),
            format!("{NP}:p.Met1_Ser3delinsGlyAlaGly"),
        ),
        (
            format!("{NP}:p.Met1_Trp5delinsGlyAlaGlyLeuGly"),
            format!("{NP}:p.Met1_Trp5delinsGlyAlaGlyLeuGly"),
        ),
    ] {
        assert_normalizes_and_reparses(&normalizer, &input, &expected);
    }
}

/// **A `Ter` or an `Xaa` in the payload** (`delins.md:45`) and **a `Ter` or an
/// `Xaa` in the reference span** (`extension.md:18`, the stop-loss ranking).
///
/// All four arms are exercised because the production guard is one predicate
/// per side — `seq.0.contains(Ter) || seq.0.contains(Xaa)` and the matching
/// `ref_aas` test in `try_protein_split_delins` — so covering only `Ter` would
/// leave half of each disjunction unasserted, and dropping `Xaa` from either
/// would go unnoticed.
///
/// The two reference-side arms use their own providers, because the residue
/// they turn on is *inserted* into the protein rather than substituted in: the
/// span named is `Ala2_Ter4` / `Ala2_Xaa4`, whose third residue is the one the
/// guard rejects.
#[test]
fn a_ter_or_xaa_on_either_side_declines_to_a_description_that_reparses() {
    // Payload side: the last residue of an otherwise splittable payload.
    for payload in ["GlySerTer", "GlySerXaa"] {
        let description = format!("{NP}:p.Ala2_Leu4delins{payload}");
        assert_normalizes_and_reparses(&normalizer(), &description, &description);
    }

    // Reference side: the span's own residues carry the `Ter`/`Xaa`.
    for (protein, span) in [
        (PROTEIN_WITH_STOP_AT_4, "Ala2_Ter4"),
        (PROTEIN_WITH_UNKNOWN_AT_4, "Ala2_Xaa4"),
    ] {
        let description = format!("{NP}:p.{span}delinsGlySerGly");
        assert_normalizes_and_reparses(
            &normalizer_over(&[(NP, protein)]),
            &description,
            &description,
        );
    }
}

/// **Unequal length** (no residue-wise correspondence) and **a span under
/// three** (no interior to be unchanged).
#[test]
fn an_ineligible_geometry_declines_to_a_description_that_reparses() {
    let normalizer = normalizer();
    for (input, expected) in [
        (
            format!("{NP}:p.Ala2_Leu4delinsGlySer"),
            format!("{NP}:p.Ala2_Leu4delinsGlySer"),
        ),
        (
            format!("{NP}:p.Ala2_Ser3delinsGlyGly"),
            format!("{NP}:p.Ala2_Ser3delinsGlyGly"),
        ),
    ] {
        assert_normalizes_and_reparses(&normalizer, &input, &expected);
    }
}

/// **Fewer than two changed runs** — the change really is one contiguous
/// `delins`, and here the affix trimming that precedes the split has already
/// reduced it to a substitution away from residue 1.
#[test]
fn one_changed_run_declines_to_a_description_that_reparses() {
    assert_normalizes_and_reparses(
        &normalizer(),
        &format!("{NP}:p.Ala2_Leu4delinsGlySerLeu"),
        &format!("{NP}:p.Ala2Gly"),
    );
}

/// **No protein reference.** Both shapes of it: a provider carrying no protein
/// data at all (`has_protein_data()` false) and one carrying a *different*
/// accession, where the fetch is what declines (#1131).
#[test]
fn an_unavailable_protein_reference_declines_to_a_description_that_reparses() {
    let input = format!("{NP}:p.Ala2_Leu4delinsGlySerGly");
    for normalizer in [
        normalizer_over(&[]),
        normalizer_over(&[("NP_999999.9", PROTEIN)]),
    ] {
        assert_normalizes_and_reparses(&normalizer, &input, &input);
    }
}

/// **Endpoints that are not certain single positions**, and a span whose window
/// runs off the end of the protein.
#[test]
fn an_unreadable_span_declines_to_a_description_that_reparses() {
    let normalizer = normalizer();
    for input in [
        format!("{NP}:p.Ala2_?delinsGlySerGly"),
        format!("{NP}:p.(Ala2_Leu4)delinsGlySerGly"),
        format!("{NP}:p.Ala19_Ala25delinsGlyAlaGlyAlaGlyAlaGly"),
    ] {
        assert_normalizes_and_reparses(&normalizer, &input, &input);
    }
}

/// A splittable delins **inside** an allele is left alone — #1606 declines
/// there rather than emit `p.[A;B];[C]`, which ferro's protein allele grammar
/// cannot read back. Asserting the re-parse is the point: this decline is the
/// one that is *about* an unparseable string, so a future change that flattened
/// the members into the outer bracket would fail here rather than only on an
/// exact-string comparison.
#[test]
fn a_delins_inside_an_allele_declines_to_a_description_that_reparses() {
    let normalizer = normalizer();
    for input in [
        format!("{NP}:p.[Ala2_Leu4delinsGlySerGly;Ala19Gly]"),
        format!("{NP}:p.[Ala2_Leu4delinsGlySerGly];[Ala19Gly]"),
    ] {
        assert_normalizes_and_reparses(&normalizer, &input, &input);
    }
}
