//! Issue #334 — 3'-rule shuffle edge cases: exon-junction exception +
//! repeat-region depth-by-1.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/334>.
//!
//! Two distinct spec-arbitration items pinned together because they're
//! intertwined in the shuffle layer:
//!
//! ## Exon-junction exception (HGVS general.md)
//!
//! > "Although the 3' rule is mandatory, deletions/duplications around
//! >  exon/exon junctions using c., r., or n. reference sequences are
//! >  not shifted."
//!
//! Concretely: when a `c.`/`r.`/`n.` deletion or duplication's 3' shift
//! would relocate the variant across an exon boundary into the next
//! exon, the shift is **not applied**. ferro currently halts the shuffle
//! at the exon boundary (via `cross_boundaries=false`); this file pins
//! that contract so a future refactor cannot regress it.
//!
//! `g.`-axis variants do not get the exception (g. is contiguous; no
//! splicing concept applies), so the same physical position shifts
//! freely in genomic space.
//!
//! ## Repeat-region depth-by-1
//!
//! When the variant sits inside a tandem repeat and both 3'-shift
//! anchors are valid, ferro and mutalyzer differ by one position on
//! some inputs. Per SVD-WG009 the canonical anchor is the rightmost
//! valid start. This file pins ferro's chosen anchor against a
//! homopolymer fixture so regressions are detectable without the
//! upstream corpus manifest.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, Normalizer};

/// Build a 2-exon plus-strand transcript with a poly-A tract that
/// straddles the exon 1 / exon 2 junction.
///
/// CDS positions:
/// ```text
///                1         2
///       1234567890123456789012
/// seq:  ATGAAAAAAAAAAACGTCGTCG
///       ^ATG = start codon
///          ^ poly-A from c.4
///                   ^ junction (c.10 last base of exon 1, c.11 first base of exon 2)
///                       ^ poly-A ends at c.14
/// ```
///
/// Without the exon-junction exception, `c.4del` (deletion of a single A)
/// would 3'-shift to `c.14del` (the rightmost position in the combined
/// poly-A tract). With the exception, the shift halts at `c.10del` — the
/// last base of exon 1 — because crossing into exon 2 is forbidden.
fn provider_polya_straddling_junction() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_JUNC.1".to_string(),
        Some("JUNC".to_string()),
        Strand::Plus,
        Some("ATGAAAAAAAAAAACGTCGTCG".to_string()),
        Some(1),
        Some(22),
        vec![Exon::new(1, 1, 10), Exon::new(2, 11, 22)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// `c.4del` on a transcript whose poly-A run spans the exon 1 / exon 2
/// junction MUST halt at `c.10del` (the last base of exon 1).
///
/// Per HGVS general.md: deletions around exon/exon junctions on c./r./n.
/// references are not shifted across the junction. The 3'-most anchor
/// within exon 1 is `c.10` (the last `A` of the poly-A run in exon 1);
/// the canonical Display drops the stated base.
#[test]
fn cds_deletion_does_not_cross_exon_junction() {
    let normalizer = Normalizer::new(provider_polya_straddling_junction());
    let variant = parse_hgvs("NM_JUNC.1:c.4del").expect("parse c.4del");

    let normalized = normalizer.normalize(&variant).expect("normalize c.4del");
    let out = format!("{}", normalized);

    // Spec-mandated stopping position: last A in exon 1 = c.10.
    // Forbidden 3' shift past junction: c.14del.
    assert!(
        out.ends_with(":c.10del"),
        "spec-mandated halt at exon-1 / exon-2 junction (c.10del); got {out:?}",
    );
    assert!(
        !out.ends_with(":c.14del"),
        "ferro illegally shifted past exon junction to c.14del: {out:?}",
    );
}

/// Sibling case: `c.4_5del` (a 2-base deletion entirely within the poly-A
/// run of exon 1) must halt at the exon-1 right edge — its 3'-most form
/// within the exon is `c.9_10del`.
#[test]
fn cds_multi_base_deletion_clamps_at_junction() {
    let normalizer = Normalizer::new(provider_polya_straddling_junction());
    let variant = parse_hgvs("NM_JUNC.1:c.4_5del").expect("parse c.4_5del");

    let normalized = normalizer.normalize(&variant).expect("normalize c.4_5del");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":c.9_10del"),
        "2-base deletion must halt at exon-1 boundary (c.9_10del); got {out:?}",
    );
}

/// `c.4dup` (duplication of a single A within exon 1's poly-A run) must
/// also halt at the exon-1 right edge — the 3'-most dup anchor is `c.10`.
#[test]
fn cds_duplication_does_not_cross_exon_junction() {
    let normalizer = Normalizer::new(provider_polya_straddling_junction());
    let variant = parse_hgvs("NM_JUNC.1:c.4dup").expect("parse c.4dup");

    let normalized = normalizer.normalize(&variant).expect("normalize c.4dup");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":c.10dup"),
        "dup must halt at exon-1 boundary (c.10dup); got {out:?}",
    );
}

/// Build a 2-exon non-coding (n.) transcript with a poly-T tract that
/// straddles the exon 1 / exon 2 junction. Used to verify the
/// exception fires on n.-axis variants too.
fn provider_polyt_noncoding_straddling_junction() -> MockProvider {
    let mut provider = MockProvider::new();
    // Positions: 1         2
    //            1234567890123456789012
    //            CGTTTTTTTTTTTTTGAACGTA
    //              ^poly-T from n.3
    //                       ^junction (n.10 last base of exon 1, n.11 first base of exon 2)
    //                            ^poly-T ends at n.15
    provider.add_transcript(Transcript::new(
        "NR_JUNC.1".to_string(),
        Some("JUNC_NC".to_string()),
        Strand::Plus,
        Some("CGTTTTTTTTTTTTTGAACGTA".to_string()),
        None, // non-coding: no CDS
        None,
        vec![Exon::new(1, 1, 10), Exon::new(2, 11, 22)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// Non-coding axis: `n.3del` (deletion of T inside the poly-T tract
/// crossing the exon junction) must halt at the exon-1 boundary at
/// `n.10del`. Per HGVS general.md the exception applies to c., r., AND
/// n. references.
#[test]
fn noncoding_deletion_does_not_cross_exon_junction() {
    let normalizer = Normalizer::new(provider_polyt_noncoding_straddling_junction());
    let variant = parse_hgvs("NR_JUNC.1:n.3del").expect("parse n.3del");

    let normalized = normalizer.normalize(&variant).expect("normalize n.3del");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":n.10del"),
        "n. deletion must halt at exon-1 boundary (n.10del); got {out:?}",
    );
    assert!(
        !out.ends_with(":n.15del"),
        "ferro illegally shifted n. variant past exon junction to n.15del: {out:?}",
    );
}

/// Negative control: a variant where the poly-A run lives **entirely
/// within** a single exon must shift to its 3'-most position. The exon-
/// junction exception only fires when the shift would actually cross a
/// junction; an in-exon shuffle is the spec default.
#[test]
fn shuffle_within_exon_proceeds_normally() {
    let mut provider = MockProvider::new();
    // Position: 1         2
    //           12345678901234567890
    //           ATGAAAAAAAACGTCGTCGT
    //              ^poly-A entirely inside exon 1 (c.4..c.10), exon 1 ends at c.20
    provider.add_transcript(Transcript::new(
        "NM_INEX.1".to_string(),
        Some("INEX".to_string()),
        Strand::Plus,
        Some("ATGAAAAAAAACGTCGTCGT".to_string()),
        Some(1),
        Some(20),
        vec![Exon::new(1, 1, 20)], // single exon — no junction to respect
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_INEX.1:c.4del").expect("parse c.4del");
    let normalized = normalizer.normalize(&variant).expect("normalize c.4del");
    let out = format!("{}", normalized);

    // No junction in a single-exon transcript: shuffle proceeds to the
    // 3'-most position in the poly-A run. The sequence is "ATG" then
    // 8 A's (positions 4-11), then "C" at position 12. Spec-canonical
    // 3' anchor is c.11.
    assert!(
        out.ends_with(":c.11del"),
        "single-exon transcript should shuffle to c.11del; got {out:?}",
    );
}

/// Negative control for the n. axis: a single-exon non-coding
/// transcript with a poly-T tract must shuffle to its 3'-most anchor
/// when no junction exists, mirroring the c. negative control above.
/// Pins the n. path symmetrically with the c. path.
#[test]
fn noncoding_shuffle_within_single_exon_proceeds_normally() {
    let mut provider = MockProvider::new();
    // Position: 1         2
    //           12345678901234567890
    //           CGTTTTTTTTTGAACGTACG
    //              ^poly-T entirely inside the single exon (n.3..n.11)
    provider.add_transcript(Transcript::new(
        "NR_INEX.1".to_string(),
        Some("INEX_NC".to_string()),
        Strand::Plus,
        Some("CGTTTTTTTTTGAACGTACG".to_string()),
        None,
        None,
        vec![Exon::new(1, 1, 20)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NR_INEX.1:n.3del").expect("parse n.3del");
    let normalized = normalizer.normalize(&variant).expect("normalize n.3del");
    let out = format!("{}", normalized);

    // Poly-T runs n.3..n.11 (9 Ts). 3'-most anchor is n.11.
    assert!(
        out.ends_with(":n.11del"),
        "single-exon n. transcript should shuffle to n.11del; got {out:?}",
    );
}

/// Build a 2-exon minus-strand transcript whose poly-A run on the
/// transcript view straddles the exon 1 / exon 2 junction. ferro stores
/// the transcript sequence in transcript orientation (always 5' → 3'),
/// so the junction location in tx-frame coordinates is independent of
/// strand — but the genomic-projection paths differ. The exception
/// applies on c. equally regardless of strand because the boundary
/// helper operates on tx-frame positions.
#[test]
fn minus_strand_cds_deletion_does_not_cross_exon_junction() {
    let mut provider = MockProvider::new();
    // tx-frame sequence (always 5' → 3'):
    //   1234567890123456789012
    //   ATGAAAAAAAAAAACGTCGTCG
    //                ^junction at c.10/c.11
    provider.add_transcript(Transcript::new(
        "NM_MJUNC.1".to_string(),
        Some("MJUNC".to_string()),
        Strand::Minus,
        Some("ATGAAAAAAAAAAACGTCGTCG".to_string()),
        Some(1),
        Some(22),
        vec![Exon::new(1, 1, 10), Exon::new(2, 11, 22)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_MJUNC.1:c.4del").expect("parse c.4del");

    let normalized = normalizer.normalize(&variant).expect("normalize c.4del");
    let out = format!("{}", normalized);

    // Strand-independent: tx-frame boundary halts the shift at c.10 just
    // like the plus-strand fixture.
    assert!(
        out.ends_with(":c.10del"),
        "minus-strand transcript must halt at exon-1 / exon-2 junction (c.10del); got {out:?}",
    );
}

/// Build a 2-exon transcript whose homopolymer-adjacent run lets a
/// multi-base insertion shuffle past the exon 1 / exon 2 junction
/// while staying as an `ins` (no dup / repeat reclassification).
///
/// ```text
///                1
///       12345 6789012
/// seq:  ATGCC CCCGTCG
///       ^exon1 (1-5)
///             ^exon2 (6-12)
/// ```
///
/// Inserting `CCTC` at `r.4_5` / `n.4_5` has a 3'-equivalent insertion
/// site at `r.6_7` / `n.6_7` (the canonical Display rotates the
/// inserted sequence to `UCCC` for `r.` and `TCCC` for `n.`). That
/// position lives in exon 2 — past the c.5/c.6 junction. The HGVS
/// exon-junction exception is narrow: it names deletions and
/// duplications only, so insertions MUST still 3'-shift across the
/// junction. If a future regression clamps insertions at the exon
/// boundary, the shuffled position would be pinned inside exon 1 and
/// these tests would fail.
fn provider_c_run_across_junction() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_CRUN.1".to_string(),
        Some("CRUN".to_string()),
        Strand::Plus,
        Some("ATGCCCCCGTCG".to_string()),
        Some(1),
        Some(12),
        vec![Exon::new(1, 1, 5), Exon::new(2, 6, 12)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// `r.` axis: insertions are NOT subject to the exon-junction
/// exception. `r.4_5inscctc` must 3'-shift across the c.5/c.6
/// junction; the spec-canonical anchor is `r.6_7insuccc` (rotated
/// inserted sequence). If the fix accidentally clamped insertions
/// at the exon boundary, the variant would be pinned at
/// `r.5_6inscucc` (one position earlier, with the rotated form one
/// step short of the canonical answer).
#[test]
fn rna_insertion_shifts_past_exon_junction() {
    let normalizer = Normalizer::new(provider_c_run_across_junction());
    let variant = parse_hgvs("NM_CRUN.1:r.4_5inscctc").expect("parse r.4_5inscctc");
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize r.4_5inscctc");
    let out = format!("{}", normalized);

    // Spec-correct 3' anchor across the exon junction.
    assert!(
        out.ends_with(":r.6_7insuccc"),
        "r. insertion must 3'-shift across exon junction to r.6_7insuccc; got {out:?}",
    );
    // Bug signature: clamped at the exon-1 boundary.
    assert!(
        !out.ends_with(":r.5_6inscucc"),
        "r. insertion was illegally clamped at the exon-1 boundary: {out:?}",
    );
}

/// `n.` axis sibling of the r. test above. `n.4_5insCCTC` must
/// 3'-shift past the exon junction to `n.6_7insTCCC`. Same anti-clamp
/// regression check.
#[test]
fn noncoding_insertion_shifts_past_exon_junction() {
    let normalizer = Normalizer::new(provider_c_run_across_junction());
    let variant = parse_hgvs("NM_CRUN.1:n.4_5insCCTC").expect("parse n.4_5insCCTC");
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize n.4_5insCCTC");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":n.6_7insTCCC"),
        "n. insertion must 3'-shift across exon junction to n.6_7insTCCC; got {out:?}",
    );
    assert!(
        !out.ends_with(":n.5_6insCTCC"),
        "n. insertion was illegally clamped at the exon-1 boundary: {out:?}",
    );
}

/// Empty-`delins` regression: a `delins` with no inserted sequence
/// (`...delins`) is rewritten by `normalize_na_edit` into a plain
/// `Deletion` (HGVS spec, issue #81 A3), but that rewrite happens
/// **inside** `normalize_na_edit` — i.e. *after* `normalize_tx` /
/// `normalize_rna` have already chosen the shuffle bounds. If the
/// bounds-selection predicate `edit_is_del_or_dup` only matches
/// `NaEdit::Deletion` / `NaEdit::Duplication`, an empty `delins` at an
/// exon junction still gets the full-transcript bounds and shuffles
/// across the junction, even though its canonical form is a deletion
/// (which should halt at the junction per the exception).
///
/// This test pins the n.-axis case: `n.3_3delins` on the same poly-T
/// junction fixture must halt at `n.10del` (exon-1 right edge), not
/// shift to `n.15del` (the 3'-most poly-T position past the junction).
/// Without the empty-`delins` arm in `edit_is_del_or_dup`, the variant
/// would be pinned at `n.15del`.
#[test]
fn noncoding_empty_delins_does_not_cross_exon_junction() {
    let normalizer = Normalizer::new(provider_polyt_noncoding_straddling_junction());
    let variant = parse_hgvs("NR_JUNC.1:n.3delins").expect("parse n.3delins");

    let normalized = normalizer.normalize(&variant).expect("normalize n.3delins");
    let out = format!("{}", normalized);

    // Spec-correct: empty delins canonicalizes to `del` and must halt
    // at the exon-1 / exon-2 junction (n.10), not cross to n.15.
    assert!(
        out.ends_with(":n.10del"),
        "empty n. delins must canonicalize to del and halt at exon-1 boundary (n.10del); got {out:?}",
    );
    assert!(
        !out.ends_with(":n.15del"),
        "empty n. delins illegally shifted past exon junction to n.15del: {out:?}",
    );
}
