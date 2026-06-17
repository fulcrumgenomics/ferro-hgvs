//! Issue #92 — `p.ins` → `p.dup` canonicalization.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/92>.
//!
//! HGVS spec (Prioritization rule): when both `ins` and `dup` can
//! describe the same change, the canonical form is `dup`. Concretely:
//! if the inserted residues of `p.<X>_<Y>ins<seq>` equal the residues
//! immediately upstream (after 3' shifting), the variant must be
//! expressed as `p.<a>_<b>dup` for the corresponding 3' anchor.
//!
//! Depends on #91 (protein 3'-shifting) so the "most 3' position" is
//! well-defined. This PR is stacked on `feat/issue-91-protein-three-prime-shift`.
//!
//! Examples (against the synthetic protein `MKVAAALELE` used in #91's
//! tests):
//! - `p.Ala4_Ala5insAla` (insert one A between p.4 and p.5) duplicates
//!   the preceding A — canonical form is `p.Ala4dup`. After 3'-shift
//!   on the AAA homopolymer at p.4..p.6, the canonical anchor is
//!   `p.Ala6dup`.
//! - `p.Glu8_Leu9insLeuGlu` (insert LE between p.8 and p.9) duplicates
//!   the preceding LE — canonical form is `p.Leu7_Glu8dup`. After
//!   3'-shift on the LELE motif, the canonical anchor is
//!   `p.Leu9_Glu10dup`.
//! - `p.Ala4_Ala5insGly` (insert one G between p.4 and p.5) does NOT
//!   duplicate any preceding residue — canonical form stays as `ins`.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// Same fixture as #91: `MKVAAALELE` (10 AAs) with a poly-A tract at
/// p.4..p.6 and a repeated `LE` motif at p.7..p.8 / p.9..p.10.
fn provider_with_polya_protein() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein("NP_TESTPROT.1", "MKVAAALELE");
    provider
}

/// Shared assertion for the three "must stay as ins" negative
/// controls below: the normalized form must keep the `ins` token
/// and must NOT have been rewritten to a `dup`. Extracted to one
/// helper so the two-half pattern (`contains("ins")` plus
/// `!contains("dup")`) and its failure messages stay consistent.
fn assert_stays_as_ins(out: &str) {
    assert!(out.contains("ins"), "variant must remain ins; got {out:?}",);
    assert!(
        !out.contains("dup"),
        "variant must NOT canonicalize to dup; got {out:?}",
    );
}

/// Single-residue ins → dup canonicalization with 3' shift.
///
/// `p.Ala4_Ala5insAla` inserts one A between p.4 and p.5. The
/// preceding residue at p.4 is A, so the insertion duplicates it. The
/// canonical `dup` form gets 3'-shifted through the AAA homopolymer to
/// `p.Ala6dup`.
#[test]
fn protein_single_ins_becomes_dup_at_3prime_anchor() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Ala4_Ala5insAla").expect("parse p.Ala4_Ala5insAla");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Ala4_Ala5insAla");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Ala6dup"),
        "ins of preceding residue must canonicalize to dup at 3' anchor; got {out:?}",
    );
}

/// Multi-residue ins → dup with range 3' shift.
///
/// `p.Glu8_Leu9insLeuGlu` inserts LE between p.8 and p.9. The
/// preceding two residues at p.7..p.8 are LE, so the insertion
/// duplicates them. The canonical `dup` form gets 3'-shifted through
/// the LELE motif to `p.Leu9_Glu10dup`.
#[test]
fn protein_range_ins_becomes_range_dup_at_3prime_anchor() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Glu8_Leu9insLeuGlu").expect("parse p.Glu8_Leu9insLeuGlu");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Glu8_Leu9insLeuGlu");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu9_Glu10dup"),
        "ins of preceding range must canonicalize to range dup at 3' anchor; got {out:?}",
    );
}

/// Negative control: an insertion that does NOT duplicate any
/// preceding residue must stay as `ins`. `p.Ala4_Ala5insGly` inserts
/// one G between p.4 and p.5 — neither neighbour is G, so the
/// canonical form remains `ins`.
#[test]
fn protein_ins_without_duplicate_stays_as_ins() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Ala4_Ala5insGly").expect("parse p.Ala4_Ala5insGly");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Ala4_Ala5insGly");
    let out = format!("{}", normalized);

    assert_stays_as_ins(&out);
}

/// Downstream-window dup canonicalization with 3' shift.
///
/// `p.Val3_Ala4insAla` inserts A between p.3 (V) and p.4 (A); the
/// preceding residue is V, not A — but the *following* residue at
/// p.4 is A, so the change still describes the same protein as
/// duplicating an A in the AAA tract. Per the HGVS Prioritization
/// rule "If both 'dup' and 'ins' can be used to describe the same
/// change, the 'dup' nomenclature is recommended", the canonical
/// form is `p.Ala4dup` → 3'-shift → `p.Ala6dup`.
#[test]
fn protein_ins_matching_following_residue_canonicalizes_to_dup() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Val3_Ala4insAla").expect("parse p.Val3_Ala4insAla");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Val3_Ala4insAla");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Ala6dup"),
        "ins of A between V and A duplicates the following A; canonical form is p.Ala6dup; got {out:?}",
    );
}

/// Partial-match negative control: only the **full** inserted
/// sequence must match a contiguous window. `p.Ala4_Ala5insAlaGly`
/// — the leading A matches the upstream A at p.4, but the full 2-AA
/// window `[p.3, p.4] = "VA"` does not match `AG`, so the rewrite
/// does NOT fire and the variant stays as `ins`.
#[test]
fn protein_ins_partial_match_stays_as_ins() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Ala4_Ala5insAlaGly").expect("parse p.Ala4_Ala5insAlaGly");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Ala4_Ala5insAlaGly");
    let out = format!("{}", normalized);

    assert_stays_as_ins(&out);
}

/// Predicted-edit wrapper: `p.(Ala4_Ala5insAla)` must canonicalize
/// the inner form (ins → dup → 3'-shift to `p.Ala6dup`) while
/// preserving the outer `(...)` predicted-edit marker. Pins that
/// the `Mu::Uncertain` wrapping is not silently dropped by the
/// canonicalization pipeline.
#[test]
fn protein_predicted_ins_canonicalizes_inside_parens() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.(Ala4_Ala5insAla)").expect("parse p.(Ala4_Ala5insAla)");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.(Ala4_Ala5insAla)");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.(Ala6dup)"),
        "predicted ins must canonicalize to predicted dup with parens preserved; got {out:?}",
    );
}

/// Insertion longer than the variant's upstream window must not
/// crash the helper. `p.Lys2_Val3insMetLysVal` inserts the whole
/// upstream window from p.1..p.3 — the upstream window check
/// requires `start_pos.number >= len`, so for len=3 and X=2 the
/// upstream window doesn't exist; the downstream window
/// `[p.3, p.5] = "VAA"` does not match `MKV`, so the variant stays
/// as `ins`.
#[test]
fn protein_ins_longer_than_upstream_window_stays_as_ins() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Lys2_Val3insMetLysVal").expect("parse p.Lys2_Val3insMetLysVal");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Lys2_Val3insMetLysVal");
    let out = format!("{}", normalized);

    assert_stays_as_ins(&out);
}

/// Met1-touching dup canonicalization must emit the
/// `InitiatorMetCanonicalization` soft-warning so downstream
/// consumers can detect the start-codon interaction. The protein
/// fixture starts with `M K V A A A L E L E`; inserting M between p.1
/// and p.2 canonicalizes to `p.Met1dup` (upstream window matches),
/// which includes position 1 → warning fires.
#[test]
fn protein_met1_dup_emits_initiator_warning() {
    use ferro_hgvs::normalize::NormalizationWarning;

    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Met1_Lys2insMet").expect("parse p.Met1_Lys2insMet");

    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize p.Met1_Lys2insMet");
    let out = format!("{}", result.result);
    assert!(
        out.ends_with(":p.Met1dup"),
        "ins→dup must canonicalize to p.Met1dup; got {out:?}",
    );
    assert!(
        result
            .warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::InitiatorMetCanonicalization { .. })),
        "Met1-touching dup must emit InitiatorMetCanonicalization; \
         got warnings: {:?}",
        result.warnings,
    );
}

/// Negative control: a non-Met1 dup must NOT emit the warning.
/// `p.Glu8_Leu9insLeuGlu` canonicalizes to `p.Leu9_Glu10dup` (from
/// the existing test `protein_range_ins_becomes_range_dup_at_3prime_anchor`),
/// which is entirely past p.1 → no warning.
#[test]
fn protein_non_met1_dup_does_not_emit_initiator_warning() {
    use ferro_hgvs::normalize::NormalizationWarning;

    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Glu8_Leu9insLeuGlu").expect("parse p.Glu8_Leu9insLeuGlu");

    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize p.Glu8_Leu9insLeuGlu");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| matches!(w, NormalizationWarning::InitiatorMetCanonicalization { .. })),
        "non-Met1 dup must NOT emit InitiatorMetCanonicalization; \
         got warnings: {:?}",
        result.warnings,
    );
}
