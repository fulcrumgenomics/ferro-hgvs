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
fn rna_two_base_u_insertion_canonicalizes_to_repeat() {
    // `n.3_4insTT` -> `n.3_5T[5]`; the `r.` analog renders the repeat unit as `u`.
    assert_eq!(normalize("NM_U.1:r.3_4insuu"), "NM_U.1:r.3_5u[5]");
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
