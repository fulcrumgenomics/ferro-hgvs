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

/// Shared assertion for the "must stay as delins" negative
/// controls below: the normalized form must keep the `delins`
/// token and must NOT have been rewritten to a `dup`. Mirrors
/// `assert_stays_as_ins` in `tests/issue_92_protein_ins_to_dup.rs`.
fn assert_stays_as_delins(out: &str) {
    assert!(
        out.contains("delins"),
        "variant must remain delins; got {out:?}",
    );
    assert!(
        !out.contains("dup"),
        "variant must NOT canonicalize to dup; got {out:?}",
    );
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

    assert_stays_as_delins(&out);
}

/// Regression: a delins with **no shared affix** whose inserted
/// residues happen to match an upstream window adjacent to the
/// deleted range must NOT be rewritten to a bare `dup` — that
/// would silently drop the deletion side of the edit and return
/// a non-equivalent variant.
///
/// On the fixture `MKVAAALELE`, `p.Ala4_Ala5delinsLysVal` deletes
/// `AA` (at p.4-5) and inserts `KV`. The inserted `KV` matches
/// the upstream window p.2-3 (also `KV`), so an unguarded
/// "treat the residual ins as a free-floating insertion and look
/// for a tandem dup" path would emit `p.Lys2_Val3dup`, dropping
/// the `p.4_5del` half of the edit. The canonical form must keep
/// the `delins` token.
#[test]
fn protein_delins_without_shared_affix_does_not_collapse_to_upstream_dup() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Ala4_Ala5delinsLysVal").expect("parse p.Ala4_Ala5delinsLysVal");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Ala4_Ala5delinsLysVal");
    let out = format!("{}", normalized);

    assert_stays_as_delins(&out);
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

/// Regression: when affix-trim leaves a zero-width del + non-empty
/// ins, and `try_protein_ins_to_dup` finds no tandem-match, the helper
/// falls through to emitting a plain insertion anchored at the
/// trimmed-range flanks. The AA letters attached to those flank
/// positions MUST match the reference at the labelled positions —
/// otherwise the helper produces a self-contradictory HGVS string
/// (positions saying one thing, AA labels saying another).
///
/// Construct on the fixture `MKVAAALELE`: `p.Leu7_Glu8delinsLeuMetGluGlu`
/// trims leading `Leu` and trailing `Glu` (lcp = lcs = 1), leaving a
/// zero-width del and a `MetGlu` residual ins anchored between p.Leu7
/// and p.Glu8 (new_start = 8). The residual `MetGlu` matches no upstream
/// (`AlaLeu`) or downstream (`GluLeu`) window, so the ins→dup helper
/// returns None and execution falls through to the plain insertion emit
/// path. Expected: `p.Leu7_Glu8insMetGlu` — Leu at p.7 and Glu at p.8
/// both match the reference. A buggy fallback that fetches the wrong
/// window would emit `p.Glu7_Leu8insMetGlu` (E and L swapped) or fail.
#[test]
fn protein_delins_fallback_ins_uses_correct_flank_residues() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8delinsLeuMetGluGlu")
        .expect("parse p.Leu7_Glu8delinsLeuMetGluGlu");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8delinsLeuMetGluGlu");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu7_Glu8insMetGlu"),
        "delins fallback must emit insertion with flank AAs matching \
         the reference at the labelled positions; got {out:?}",
    );
}

/// Regression: the fallback insertion emit path runs unsigned
/// arithmetic on `new_start`. When `new_start <= 1` the fallback
/// cannot label a left flank (position 0 is invalid in HGVS protein
/// coordinates) and the helper must bail out cleanly — without the
/// guard, the fallback emits an HGVS string anchored at position 0
/// (e.g. `p.Met0_Lys1ins...`) which downstream tools cannot parse.
///
/// Construct on the fixture `MKVAAALELE`: `p.Met1delinsValMet` has
/// ref = `Met` and ins = `ValMet`. Affix-trim finds lcp = 0 (V vs M)
/// and lcs = 1 (trailing Met), leaving residual_del = empty and
/// residual_seq = `Val`. The trimmed anchor is new_start = 1, so
/// flank_start = 0 and the upstream ins→dup block is skipped
/// (`flank_start >= 1` is false). Execution reaches the fallback code
/// path with new_start == 1, where the helper must not emit a
/// position-0 label.
#[test]
fn protein_delins_fallback_ins_rejects_position_zero_label() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Met1delinsValMet").expect("parse p.Met1delinsValMet");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Met1delinsValMet");
    let out = format!("{}", normalized);

    assert!(
        !out.contains("Met0") && !out.contains("_0"),
        "delins fallback must not emit a position-0 protein label \
         (positions are 1-based in HGVS); got {out:?}",
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
    assert_stays_as_delins(&out);
}
