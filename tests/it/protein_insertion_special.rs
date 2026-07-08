//! Protein insertion special inserted-sequences (issue #955).
//!
//! The spec (`protein/insertion.md:21,45-61`) allows an inserted protein
//! sequence to be described by `insXaa[n]` (n unknown residues, open reading
//! frame insertion) or `ins*<n>` / `insTer<n>` (an inserted sequence ending in a
//! translation stop at 1-based position n within the insertion). Stop is
//! canonicalized to `Ter<n>` on Display, matching frameshift/extension.

use ferro_hgvs::parse_hgvs;

/// Parse `input`, render it back, and assert the canonical Display equals
/// `expected` (equal to `input` for forms already in canonical shape).
fn assert_canonical(input: &str, expected: &str) {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("failed to parse {input:?}: {e}"));
    assert_eq!(
        variant.to_string(),
        expected,
        "canonical Display mismatch for {input:?}"
    );
}

#[test]
fn ins_xaa_repeat_count_round_trips() {
    assert_canonical(
        "NP_003997.1:p.(Val582_Asn583insXaa[5])",
        "NP_003997.1:p.(Val582_Asn583insXaa[5])",
    );
    assert_canonical(
        "NP_003997.1:p.Arg78_Gly79insXaa[23]",
        "NP_003997.1:p.Arg78_Gly79insXaa[23]",
    );
    assert_canonical(
        "NP_003997.1:p.Lys2_Leu3insXaa[34]",
        "NP_003997.1:p.Lys2_Leu3insXaa[34]",
    );
}

#[test]
fn ins_xaa_with_uncertain_interval_parses() {
    // Parens wrap the *interval* here (edit outside): `(Ala123_Pro131)insXaa[4]`.
    // Confirm the enriched inserted-sequence parses and round-trips byte-exactly
    // under an uncertain interval (not just a substring match).
    assert_canonical(
        "NP_003997.1:p.(Ala123_Pro131)insXaa[4]",
        "NP_003997.1:p.(Ala123_Pro131)insXaa[4]",
    );
}

#[test]
fn ins_stop_length_round_trips_and_canonicalizes_to_ter() {
    // `insTer<n>` is already canonical.
    assert_canonical(
        "NP_003997.1:p.Lys2_Leu3insTer12",
        "NP_003997.1:p.Lys2_Leu3insTer12",
    );
    // `ins*<n>` is an accepted input alternative; Display canonicalizes `*` to
    // `Ter`, matching frameshift/extension stop rendering.
    assert_canonical(
        "NP_060250.2:p.Gln746_Lys747ins*63",
        "NP_060250.2:p.Gln746_Lys747insTer63",
    );
}

#[test]
fn ins_repeat_rejects_non_spec_forms() {
    // The spec sanctions only `insXaa[n]` with an unknown residue and an exact
    // positive count. Everything else must be rejected, not silently accepted as
    // a `Repeat` (which would, e.g., let `insTer[5]` evade stop detection).
    for bad in [
        "NP_003997.1:p.Lys2_Leu3insXaa[0]",   // zero count
        "NP_003997.1:p.Lys2_Leu3insXaa[?]",   // uncertain count
        "NP_003997.1:p.Lys2_Leu3insXaa[2_5]", // range count
        "NP_003997.1:p.Lys2_Leu3insTer[5]",   // stop is insTer<n>, not a repeat
        "NP_003997.1:p.Lys2_Leu3insGln[5]",   // known residue is not the Xaa repeat form
    ] {
        assert!(
            parse_hgvs(bad).is_err(),
            "non-spec protein insertion repeat form must be rejected: {bad:?}"
        );
    }
}

#[test]
fn plain_amino_acid_insertion_still_round_trips() {
    // Regression: ordinary literal insertions (incl. single Xaa / Ter) unchanged.
    assert_canonical(
        "NP_003997.1:p.Lys2_Leu3insGlnSerLys",
        "NP_003997.1:p.Lys2_Leu3insGlnSerLys",
    );
    assert_canonical(
        "NP_003997.1:p.Ser332_Ser333insXaa",
        "NP_003997.1:p.Ser332_Ser333insXaa",
    );
}
