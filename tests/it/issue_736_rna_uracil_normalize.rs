//! Regression test for issue #736: `normalize_rna` compares an `r.` edit's
//! literal bases against the DNA-stored reference byte-for-byte, but the parser
//! keeps the RNA base `u` as a distinct [`Base::U`] (`b'U'`), separate from
//! [`Base::T`] (`b'T'`). So any `r.` edit carrying a literal `u` never matches
//! the reference `T`, and insertions / delins fail to canonicalize or 3'-shift.
//!
//! The spec-correct `u` form produced a *wrong* result while the
//! technically-invalid `T` form produced the right one — e.g. on a `TTT` run,
//! `r.3_4insu` stayed `r.3_4insu` instead of collapsing to `r.5dup` (which is
//! exactly what `r.3_4insT` already does).
//!
//! The fix maps the `r.` edit's literal `U` bases to `T` before the DNA-based
//! normalization, then `RnaVariant`'s Display (`to_rna_string`) renders `T`→`u`
//! on output. These tests pin the `u`-form results against the DNA-equivalent
//! `n.`/`r.`-with-`T` outputs on the same transcript.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Single-exon coding transcript `AATTTGCC` with a `TTT` run at tx 3-5 and
/// `cds_start = 1` (so `r.N == c.N == n.N`). An insertion of one more thymine
/// into the run canonicalizes to a duplication at the run's 3' end.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_U.1".to_string(),
        Some("UG".to_string()),
        Strand::Plus,
        "AATTTGCC".to_string(),
        Some(1u64),
        Some(8u64),
        vec![Exon::new(1, 1, 8)],
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

fn normalize(input: &str) -> String {
    let normalizer = Normalizer::new(provider());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize {input:?} failed: {e}"));
    format!("{}", normalized)
}

#[test]
fn rna_insertion_with_u_collapses_to_dup() {
    // Inserting a single `u` into the `TTT` run must collapse to the 3'-most
    // duplication, exactly as `r.3_4insT` / `n.3_4insT` -> `*.5dup`.
    assert_eq!(normalize("NM_U.1:r.3_4insu"), "NM_U.1:r.5dup");
}

#[test]
fn rna_two_base_u_insertion_is_gated_like_the_c_axis() {
    // Re-blessed by #1192. The original expectation here was `r.3_5u[5]`,
    // justified as "the `r.` analog of `n.3_4insTT` -> `n.3_5T[5]`". That
    // compared `r.` against the wrong axis: `NM_U.1` is coding
    // (`cds_start = 1`), so `r.` is CDS-relative — the `c.` axis, not the `n.`
    // axis (#469). `RNA/repeated.md` L24-27 forbids repeat notation on a coding
    // RNA reference for units whose length is not a multiple of 3, and a
    // homopolymer unit is length 1, so the codon-frame gate must fire. On the
    // same transcript and bases the three axes now read:
    //   c.3_4insTT  -> c.4_5dup     (gated: coding, unit not codon-aligned)
    //   r.3_4insuu  -> r.4_5dup     (must match c.; this assertion)
    //   n.3_4insTT  -> n.3_5T[5]    (no reading frame, so repeat notation is fine)
    // Re-blessed again by #1204: the gate refuses repeat notation, and the
    // fallback is the `dup` the same spec sentence prescribes rather than an `ins`
    // literal — here over `r.4_5`, the 3'-most pair of the `TTT` run at r.3_5.
    // Note the canonical `dup` carries no bases at all, so the `u`-spelling
    // round trip that #736 is about is pinned by the single-`u` case above (which
    // has always collapsed to `r.5dup`) and by the `insgg` case below, not here.
    assert_eq!(normalize("NM_U.1:r.3_4insuu"), "NM_U.1:r.4_5dup");
}

/// The two sibling axes the re-blessing above is justified against. They were
/// asserted only in a comment, which is what let the original `r.3_5u[5]`
/// expectation look defensible: nothing pinned the `c.` value it was supposed to
/// match, or the `n.` value it was actually copied from. Pinning both makes the
/// three-way relationship checkable — `r.` must track `c.` (CDS-relative on a
/// coding transcript, so gated) and must *not* track `n.` (no reading frame, so
/// repeat notation is legitimate there).
#[test]
fn c_and_n_axes_bracket_the_gated_r_axis_expectation() {
    assert_eq!(
        normalize("NM_U.1:c.3_4insTT"),
        "NM_U.1:c.4_5dup",
        "c. is gated: coding footprint, unit length 1 is not codon-aligned, so \
         repeat notation is refused and the fallback is the prescribed dup (#1204)",
    );
    assert_eq!(
        normalize("NM_U.1:n.3_4insTT"),
        "NM_U.1:n.3_5T[5]",
        "n. has no reading frame, so repeat notation stays available",
    );
}

#[test]
fn rna_delins_with_u_simplifies() {
    // `n.3_5delinsTT` -> `n.5del`; the `u` form must simplify identically.
    assert_eq!(normalize("NM_U.1:r.3_5delinsuu"), "NM_U.1:r.5del");
}

#[test]
fn rna_insertion_with_t_already_works() {
    // Control: the technically-invalid `T` spelling already normalized
    // correctly (the comparison matched the DNA reference); it must keep doing so.
    assert_eq!(normalize("NM_U.1:r.3_4insT"), "NM_U.1:r.5dup");
}

#[test]
fn rna_insertion_without_u_is_unaffected() {
    // Control: a non-`u` insertion that does not touch the run must be untouched
    // by the U->T mapping.
    assert_eq!(normalize("NM_U.1:r.3_4insgg"), "NM_U.1:r.3_4insgg");
}

#[test]
fn rna_delins_uracil_reverse_complement_recognized_as_inversion() {
    // The same root cause also blocked RNA inversion recognition when the
    // reverse complement contains `u`. `r.1_2` ref is `aa`, whose reverse
    // complement is `uu`, so the canonical form is `inv`. Pre-#736 the inserted
    // `uu` (`Base::U`) was compared against the DNA reverse complement `TT` and
    // never matched, leaving a `delinsuu`; with the U->T normalization the RNA
    // path now matches the DNA axes (`n.1_2delinsTT` -> `n.1_2inv`).
    assert_eq!(normalize("NM_U.1:r.1_2delinsuu"), "NM_U.1:r.1_2inv");
}
