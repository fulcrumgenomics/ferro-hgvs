//! Issue #92 — `p.delins` → canonical-form rewrites via affix-trim.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/92>.
//!
//! HGVS spec mandates that a `delins` whose inserted residues share a
//! prefix or suffix with the deleted residues must be re-described
//! with the trimmed range and trimmed insert (`delins.md:53-56`). After
//! affix-trim, the residual is canonicalized per Prioritization
//! (`general.md:56-57`): `sub > del > dup > ins > delins`.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// Synthetic protein `MKVAAALELE` (10 AAs): poly-A at p.4..p.6,
/// LELE motif at p.7..p.10. Matches the fixture in
/// `tests/issue_92_protein_ins_to_dup.rs`.
fn provider_with_polya_protein() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein("NP_TESTPROT.1", "MKVAAALELE");
    provider
}

/// Affix-trim collapses a delins whose inserted sequence equals the
/// deleted sequence to a silent identity. `p.Leu7_Glu8delinsLeuGlu`
/// over ref `LE` should normalize to `p.Leu7_Glu8=`.
#[test]
fn protein_delins_collapses_to_identity_when_insert_equals_delete() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8delinsLeuGlu").expect("parse p.Leu7_Glu8delinsLeuGlu");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8delinsLeuGlu");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu7_Glu8=") || out.ends_with(":p.="),
        "delins where insert == delete must collapse to identity; got {out:?}",
    );
    assert!(
        !out.contains("delins"),
        "delins must NOT survive the affix-trim collapse; got {out:?}",
    );
}

/// Affix-trim with shared *trailing* residue collapses a delins to a
/// pure deletion. `p.Leu7_Glu8delinsGlu` over ref `LE` (trailing `E`
/// matches both sides) should normalize to `p.Leu7del`.
#[test]
fn protein_delins_collapses_to_deletion_with_shared_suffix() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8delinsGlu").expect("parse p.Leu7_Glu8delinsGlu");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8delinsGlu");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu7del"),
        "delins with shared suffix Glu must trim to p.Leu7del; got {out:?}",
    );
    assert!(
        !out.contains("delins"),
        "delins must be rewritten as pure del; got {out:?}",
    );
}

/// Affix-trim leaving a 1-residue del + 1-residue ins must collapse to
/// a substitution per Prioritization `sub > delins`. Construct on the
/// `LELE` motif: `p.Leu7_Glu8delinsValGlu` (ref `LE`, trim trailing
/// `Glu`) leaves del `L` + ins `Val` → `p.Leu7Val`.
#[test]
fn protein_delins_collapses_to_substitution_when_residual_is_single() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8delinsValGlu").expect("parse p.Leu7_Glu8delinsValGlu");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8delinsValGlu");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu7Val"),
        "delins with shared suffix Glu and 1-residue residual must \
         collapse to sub p.Leu7Val; got {out:?}",
    );
}

/// Affix-trim leaving a zero-width del + non-empty ins must route the
/// residual through the existing `try_protein_ins_to_dup`. Construct
/// on the `LELE` motif: `p.Leu7_Glu8delinsLeuGluLeuGlu` (ref `LE`,
/// trim leading `LeuGlu`) leaves an ins of `LeuGlu` between p.8 and
/// p.9. The upstream window ref[7..8] = `LE` matches → dup over
/// p.7..p.8. Then 3'-shift through the LELE motif lands at
/// `p.Leu9_Glu10dup`.
#[test]
fn protein_delins_routes_residual_ins_to_dup() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8delinsLeuGluLeuGlu")
        .expect("parse p.Leu7_Glu8delinsLeuGluLeuGlu");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8delinsLeuGluLeuGlu");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu9_Glu10dup"),
        "delins whose residual is a tandem ins must canonicalize to \
         dup at the 3'-shifted anchor; got {out:?}",
    );
}

/// A genuine delins with no shared affixes and no tandem-match must
/// stay as delins. `p.Leu7_Glu8delinsAlaGly` over ref `LE` shares no
/// prefix/suffix with `AlaGly`, and `AlaGly` does not match any
/// upstream or downstream window. Per Prioritization, delins is the
/// last-resort form, so it survives.
#[test]
fn protein_delins_without_shared_affix_stays_as_delins() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8delinsAlaGly").expect("parse p.Leu7_Glu8delinsAlaGly");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8delinsAlaGly");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu7_Glu8delinsAlaGly"),
        "genuine delins (no shared affix, no tandem match) must \
         survive normalization; got {out:?}",
    );
}

/// Predicted-edit wrapper: `p.(Leu7_Glu8delinsLeuGluLeuGlu)` must
/// canonicalize the inner form (delins → trim → ins → dup → 3'-shift)
/// while preserving the outer `(...)` predicted-edit marker. Pins
/// that `Mu::Uncertain` wrapping is not silently dropped by the
/// delins helper, mirroring the existing
/// `protein_predicted_ins_canonicalizes_inside_parens` invariant.
#[test]
fn protein_predicted_delins_canonicalizes_inside_parens() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.(Leu7_Glu8delinsLeuGluLeuGlu)")
        .expect("parse p.(Leu7_Glu8delinsLeuGluLeuGlu)");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.(Leu7_Glu8delinsLeuGluLeuGlu)");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.(Leu9_Glu10dup)"),
        "predicted delins must canonicalize to predicted dup with \
         parens preserved; got {out:?}",
    );
}

/// Uncertain position endpoints carry semantics that the affix-trim
/// rewrite would silently drop (e.g. `(Leu7)_Glu8` parenthesizing
/// only the start would not survive a rewrite to a smaller range).
/// The helper must reject these and let the input pass through.
#[test]
fn protein_delins_with_uncertain_start_passes_through() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.(Leu7)_Glu8delinsLeuGlu")
        .expect("parse p.(Leu7)_Glu8delinsLeuGlu");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.(Leu7)_Glu8delinsLeuGlu");
    let out = format!("{}", normalized);

    // The exact Display form may render the partial-uncertainty
    // differently across parser/Display revisions, but the helper
    // must not collapse this to a smaller form.
    assert!(
        out.contains("delins"),
        "delins with partial-uncertain endpoint must NOT be \
         canonicalized; got {out:?}",
    );
}
