//! The protein render style must reach an allele's members.
//!
//! Two rules govern how a `p.` description is written, and both were originally
//! implemented only on the single-variant path — an allele of protein members
//! reaches a different renderer and silently diverged from each:
//!
//! 1. **No gene-symbol selector** (`refseq.md:41`). Fixed for alleles by #1142
//!    and single-sourced by #1141's `gene_selector_visible`, which composes the
//!    accession-level and axis-level halves. Pinned here across the *styled*
//!    path too, which reaches alleles for the first time now that they are
//!    styled at all.
//! 2. **The amino-acid render style** — what this file's remaining tests fix.
//!    `to_styled_string` applied a caller's one-letter/three-letter and
//!    `Ter`/`*` preferences to a single protein variant but sent an allele to
//!    unstyled `Display`, so a bracketed consequence ignored the request. The
//!    justification was that "the projector's `.protein` field is always a
//!    single `HgvsVariant::Protein`, never an allele" — untrue since separated
//!    residues began rendering as a cis allele (#1095 for changed residues,
//!    #1099 for silent ones), which the projector emits routinely.

use ferro_hgvs::{parse_hgvs, AaCode, ProteinRenderStyle, TerStyle};

fn one_letter() -> ProteinRenderStyle {
    ProteinRenderStyle {
        aa_code: AaCode::One,
        stop: TerStyle::Ter,
    }
}

fn styled(input: &str) -> String {
    parse_hgvs(input)
        .unwrap_or_else(|e| panic!("{input} should parse: {e}"))
        .to_styled_string(one_letter())
}

/// The baseline that already worked: a single protein variant honours the style.
#[test]
fn a_single_protein_variant_honours_the_render_style() {
    assert_eq!(styled("NP_003997.2:p.(Phe41Cys)"), "NP_003997.2:p.(F41C)");
}

/// An allele must honour it too — the members are the things carrying the
/// amino acids, so falling through to unstyled `Display` ignores the request.
#[test]
fn a_protein_cis_allele_honours_the_render_style() {
    assert_eq!(
        styled("NP_003997.2:p.[(Phe41Cys);(Gly47Arg)]"),
        "NP_003997.2:p.[(F41C);(G47R)]"
    );
}

#[test]
fn a_protein_trans_allele_honours_the_render_style() {
    assert_eq!(
        styled("NP_003997.2:p.[(Phe41Cys)];[(Gly47Arg)]"),
        "NP_003997.2:p.[(F41C)];[(G47R)]"
    );
}

/// The remaining parseable phases, so the style reaches every renderer arm and
/// not just the two common ones. `(;)` exercises the unknown-phase compact
/// path; `^` renders *expanded*, so it exercises the full-member arm (accession
/// included) rather than the loc-edit-only one.
#[test]
fn the_remaining_allele_phases_honour_the_render_style() {
    assert_eq!(
        styled("NP_003997.2:p.(Phe41Cys)(;)(Gly47Arg)"),
        "NP_003997.2:p.(F41C)(;)(G47R)"
    );
    assert_eq!(
        styled("NP_003997.2:p.(Phe41Cys)^(Gly47Arg)"),
        "NP_003997.2:p.(F41C)^NP_003997.2:p.(G47R)"
    );
}

/// The identity shapes #1099 emits are protein alleles too.
#[test]
fn a_silent_protein_allele_honours_the_render_style() {
    assert_eq!(
        styled("NP_003997.2:p.[(Phe41=);(Gly47=)]"),
        "NP_003997.2:p.[(F41=);(G47=)]"
    );
}

/// `Ter` styling reaches allele members as well.
#[test]
fn ter_style_reaches_allele_members() {
    let star = ProteinRenderStyle {
        aa_code: AaCode::Three,
        stop: TerStyle::Star,
    };
    let v = parse_hgvs("NP_003997.2:p.[(Phe41Ter);(Gly47Arg)]").expect("parses");
    assert_eq!(
        v.to_styled_string(star),
        "NP_003997.2:p.[(Phe41*);(Gly47Arg)]"
    );
}

/// A non-protein allele is unaffected: the style has no meaning off the `p.`
/// axis, so rendering must be byte-identical to `Display`.
#[test]
fn a_non_protein_allele_is_unchanged_by_the_style() {
    for input in [
        "NC_000013.11(FLT3):g.[12345A>G;12400C>T]",
        "NM_000088.3:c.[100A>G;200C>T]",
    ] {
        let v = parse_hgvs(input).expect("parses");
        assert_eq!(v.to_styled_string(one_letter()), v.to_string(), "{input}");
    }
}

/// The gene-symbol rule holds on every protein rendering path, styled included
/// — one predicate decides it, so the styled renderer cannot reintroduce the
/// selector that #1142 removed from the plain one.
#[test]
fn no_protein_rendering_path_emits_a_gene_symbol_selector() {
    for input in [
        "NP_003997.2(DMD):p.(Phe41Cys)",
        "NP_003997.2(DMD):p.[(Phe41Cys);(Gly47Arg)]",
        "NP_003997.2(DMD):p.[(Phe41Cys)];[(Gly47Arg)]",
        // The unknown-phase and product phases reach the styled renderer via
        // different arms (compact and expanded); both must still suppress.
        "NP_003997.2(DMD):p.(Phe41Cys)(;)(Gly47Arg)",
        "NP_003997.2(DMD):p.(Phe41Cys)^(Gly47Arg)",
        "YP_003024026.1(MT-ND1):p.[(Ala52Thr);(Gly300Arg)]",
    ] {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} should parse: {e}"));
        for rendered in [v.to_string(), v.to_styled_string(one_letter())] {
            assert!(
                !rendered.contains("(DMD)") && !rendered.contains("(MT-ND1)"),
                "{input} rendered a gene-symbol selector: {rendered}"
            );
        }
    }
}
