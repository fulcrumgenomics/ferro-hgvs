//! Pin protein stop-codon glyph canonicalization (#453).
//!
//! Per HGVS v21 `general.md:53-54`:
//!
//! > "three-letter amino acid code is **preferred** ... the `*` can be
//! > used to indicate the translation stop codon in **both** one- and
//! > three-letter amino acid code descriptions."
//!
//! Both `Ter` and `*` are explicitly valid; spec prefers the
//! three-letter form. ferro Display canonicalizes `*` to `Ter` on parse,
//! honoring the spec's preference.
//!
//! Mutalyzer preserves the input glyph verbatim. ferro's canonicalize-
//! to-`Ter` policy is more spec-directive but loses input-glyph
//! information; this divergence is documented (not a bug).
//!
//! Source: phase4 spec-coverage triage divergence #9.

use ferro_hgvs::parse_hgvs;

fn round_trip(input: &str) -> String {
    parse_hgvs(input)
        .unwrap_or_else(|e| panic!("{input:?} must parse: {e}"))
        .to_string()
}

// ---- Substitution: * canonicalizes to Ter ----

/// Asterisk glyph on substitution canonicalizes to three-letter `Ter`.
#[test]
fn sub_asterisk_canonicalizes_to_ter() {
    assert_eq!(round_trip("NP_003997.1:p.Trp41*"), "NP_003997.1:p.Trp41Ter");
}

/// Three-letter `Ter` round-trips verbatim.
#[test]
fn sub_ter_round_trips_verbatim() {
    assert_eq!(
        round_trip("NP_003997.1:p.Trp41Ter"),
        "NP_003997.1:p.Trp41Ter"
    );
}

/// Asterisk glyph on a 1-letter-mode substitution still canonicalizes.
#[test]
fn sub_asterisk_one_letter_mode_canonicalizes() {
    assert_eq!(round_trip("NP_003997.1:p.W41*"), "NP_003997.1:p.Trp41Ter");
}

// ---- Extension: same rule applies per general.md:54 ----
//
// extension.md:30 explicitly shows both forms as equivalent:
// `p.Ter110GlnextTer17` ≡ `p.*110Glnext*17`.

#[test]
fn ext_asterisk_canonicalizes_to_ter() {
    // Pin the equivalence; ferro picks the three-letter form. Use the
    // spec example `p.*110Glnext*17` (both glyphs as `*`) so the test
    // matches the cited extension.md:30 form above.
    let star_out = round_trip("NP_003997.1:p.*110Glnext*17");
    let ter_out = round_trip("NP_003997.1:p.Ter110GlnextTer17");
    assert_eq!(star_out, ter_out, "* and Ter equivalent on extension");
    assert!(
        ter_out.contains("Ter110") && ter_out.contains("extTer17"),
        "canonical form uses three-letter Ter; got: {ter_out}"
    );
    // Also cover the mixed-glyph form (`*` first stop, `Ter` second):
    // canonicalization must normalize both glyphs uniformly.
    assert_eq!(
        round_trip("NP_003997.1:p.*110GlnextTer17"),
        ter_out,
        "mixed * / Ter extension form canonicalizes identically",
    );
}

// ---- Frameshift: frameshift.md:33-34 also shows both forms equivalent ----
//
// `p.Arg97ProfsTer23` ≡ `p.Arg97Profs*23`.

#[test]
fn fs_asterisk_canonicalizes_to_ter() {
    let star_out = round_trip("NP_003997.1:p.Arg97Profs*23");
    let ter_out = round_trip("NP_003997.1:p.Arg97ProfsTer23");
    assert_eq!(star_out, ter_out, "* and Ter equivalent on frameshift");
    assert!(
        ter_out.contains("ProfsTer23"),
        "canonical form uses three-letter Ter; got: {ter_out}"
    );
}

// ---- Inside allele brackets the rule still applies ----

#[test]
fn allele_bracket_asterisk_canonicalizes_to_ter() {
    let out = round_trip("NP_003997.1:p.[Trp24Cys;Lys41*]");
    assert!(
        out.contains("Lys41Ter"),
        "* inside allele bracket should canonicalize to Ter; got: {out}"
    );
}

// ---- Mixed-glyph round-trip (Ter outside, * inside the same allele) ----

#[test]
fn allele_mixed_glyphs_normalize_uniformly() {
    let mixed = round_trip("NP_003997.1:p.[Trp24Ter;Lys41*]");
    // Both arms must come out three-letter on canonical Display.
    assert!(
        mixed.contains("Trp24Ter") && mixed.contains("Lys41Ter"),
        "both stops canonicalize to Ter; got: {mixed}"
    );
}
