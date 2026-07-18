//! Issue #1063 — the predicted/uncertain `(...)` wrapper (`Mu::Uncertain`)
//! was silently dropped by the shift/canonicalization pipeline.
//!
//! Root cause: the final-result reconstructions in the `normalize_*`
//! functions (one per coordinate system) rebuild the output variant with
//! `LocEdit::new(...)`, which hardcodes `Mu::Certain` — discarding whatever
//! certainty the input variant's edit carried. So `g.(8del)` normalized to
//! `g.12del` instead of `g.(12del)`.
//!
//! The policy (per the design doc): uncertainty is the epistemic status of
//! the *described event*.
//! - For single-edit outputs it attaches to the edit: `g.(12del)`.
//! - For canonical-split outputs (a delins that decomposes into a
//!   multi-member allele) it attaches to the *whole allele*, not to the
//!   individual members: `g.[(200C>T;202C>G)]` (the predicted cis-allele
//!   notation from `recommendations/uncertain.md`, already exercised by
//!   `tests/it/predicted_cis_allele.rs`), never `g.[(200C>T);(202C>G)]`.
//!
//! This suite exercises all three shapes (3' shift, delins->inv
//! canonicalization, delins->split) plus idempotency, using the genomic
//! `JsonProvider` helper pattern from `issue_1041_repro.rs`.

use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

/// Genome-capable provider: one contig named `c` of `len` cyclic ACGT bytes
/// with `payload` written 1-based at `pos1`. Mirrors the helper in
/// `issue_1041_repro.rs`.
fn provider_len(len: usize, pos1: usize, payload: &str) -> JsonProvider {
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(len).collect();
    for (i, b) in payload.bytes().enumerate() {
        bases[pos1 - 1 + i] = b;
    }
    let seq = String::from_utf8(bases).unwrap();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { "c": seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// A contig named `contig` with a homopolymer run of `A`s at 1-based [195, 205]
/// in an otherwise cyclic-ACGT 250bp contig. Used to exercise 3' shifting.
fn provider_homopolymer_named(contig: &str) -> JsonProvider {
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(250).collect();
    for b in bases.iter_mut().take(205).skip(194) {
        *b = b'A';
    }
    let seq = String::from_utf8(bases).unwrap();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// The genome-axis homopolymer provider (contig `c`).
fn provider_homopolymer() -> JsonProvider {
    provider_homopolymer_named("c")
}

fn norm(p: &JsonProvider, input: &str) -> String {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e}"));
    Normalizer::new(p.clone())
        .normalize(&v)
        .unwrap_or_else(|e| panic!("normalize {input:?} failed: {e}"))
        .to_string()
}

#[track_caller]
fn assert_idempotent(p: &JsonProvider, input: &str) {
    let once = norm(p, input);
    let twice = norm(p, &once);
    assert_eq!(
        twice, once,
        "not idempotent for {input:?}: {once:?} -> {twice:?}"
    );
}

// =============================================================================
// Shift: uncertain del/dup 3'-shifts and keeps its wrapper.
// =============================================================================

#[test]
fn certain_del_shifts_unwrapped() {
    let p = provider_homopolymer();
    assert_eq!(norm(&p, "c:g.196del"), "c:g.205del");
}

#[test]
fn uncertain_del_shifts_and_stays_wrapped() {
    let p = provider_homopolymer();
    assert_eq!(norm(&p, "c:g.(196del)"), "c:g.(205del)");
}

#[test]
fn uncertain_del_shift_is_idempotent() {
    let p = provider_homopolymer();
    assert_idempotent(&p, "c:g.(196del)");
}

#[test]
fn certain_dup_shifts_unwrapped() {
    // A single-base dup inside the [195,205] `A` run 3'-shifts to the 3'-most
    // position (g.205dup). Assert the exact form directly so a wrong shift is
    // caught — not just "ends in dup, no parens".
    let p = provider_homopolymer();
    assert_eq!(norm(&p, "c:g.196dup"), "c:g.205dup");
}

#[test]
fn uncertain_dup_shifts_and_stays_wrapped() {
    // Assert the exact 3'-shifted, wrapped form directly (not derived from the
    // certain output) so both tests fail independently if the dup shift breaks.
    let p = provider_homopolymer();
    assert_eq!(norm(&p, "c:g.(196dup)"), "c:g.(205dup)");
}

#[test]
fn uncertain_dup_shift_is_idempotent() {
    let p = provider_homopolymer();
    assert_idempotent(&p, "c:g.(196dup)");
}

// =============================================================================
// Canonicalization: uncertain whole-run-revcomp delins -> inv, still wrapped.
// =============================================================================

#[test]
fn certain_revcomp_delins_canonicalizes_to_inv_unwrapped() {
    // ref[200..202] = GCA, revcomp = TGC.
    let p = provider_len(600, 200, "GCA");
    assert_eq!(norm(&p, "c:g.200_202delinsTGC"), "c:g.200_202inv");
}

#[test]
fn uncertain_revcomp_delins_canonicalizes_to_inv_and_stays_wrapped() {
    let p = provider_len(600, 200, "GCA");
    assert_eq!(norm(&p, "c:g.(200_202delinsTGC)"), "c:g.(200_202inv)");
}

#[test]
fn uncertain_revcomp_delins_to_inv_is_idempotent() {
    let p = provider_len(600, 200, "GCA");
    assert_idempotent(&p, "c:g.(200_202delinsTGC)");
}

// =============================================================================
// Canonical split: uncertain delins that decomposes into a 2-member allele.
// The wrapper moves to the WHOLE ALLELE, not to individual members.
// =============================================================================

#[test]
fn certain_delins_splits_into_unwrapped_allele() {
    // ref[200..202] = CAC, alt = TAG: pos0 C->T (diff), pos1 A==A (identity,
    // gap), pos2 C->G (diff) -> decomposes into two independent subs
    // (general.md:56 edit-priority / #165 sub-only decomposition).
    let p = provider_len(600, 200, "CAC");
    assert_eq!(norm(&p, "c:g.200_202delinsTAG"), "c:g.[200C>T;202C>G]");
}

#[test]
fn uncertain_delins_splits_into_allele_level_wrapped_form() {
    let p = provider_len(600, 200, "CAC");
    // Predicted cis-allele notation (`recommendations/uncertain.md`):
    // `[(<members>)]`, matching `tests/it/predicted_cis_allele.rs`.
    assert_eq!(norm(&p, "c:g.(200_202delinsTAG)"), "c:g.[(200C>T;202C>G)]");
}

#[test]
fn uncertain_delins_split_members_are_not_individually_wrapped() {
    let p = provider_len(600, 200, "CAC");
    let out = norm(&p, "c:g.(200_202delinsTAG)");
    // The member bodies (`200C>T`, `202C>G`) must not carry their own
    // parens; only the allele-level `[(...)]` wrapper is present.
    assert!(
        !out.contains("(200C>T)") && !out.contains("(202C>G)"),
        "members must not be individually wrapped: {out}"
    );
}

#[test]
fn uncertain_delins_split_is_idempotent() {
    let p = provider_len(600, 200, "CAC");
    assert_idempotent(&p, "c:g.(200_202delinsTAG)");
}

// =============================================================================
// Axis coverage: the wrapper reconstruction is per-coordinate-system (one
// `LocEdit::with_uncertainty(... with_same_certainty ...)` site per axis, each
// with its own position type). The genome cases above exercise `normalize_genome`
// (GenomePos); these exercise the `normalize_mt` (GenomePos, distinct function)
// and `normalize_cds` (CdsPos, distinct position type) reconstruction sites so a
// copy-paste slip in one axis can't hide behind the genome tests.
// =============================================================================

#[test]
fn uncertain_mito_del_shifts_and_stays_wrapped() {
    // `m.` axis (`normalize_mt`): a del in the [195,205] `A` run 3'-shifts to
    // g.205 and must keep its `(...)` wrapper.
    let p = provider_homopolymer_named("NC_012920.1");
    assert_eq!(norm(&p, "NC_012920.1:m.(196del)"), "NC_012920.1:m.(205del)");
}

#[test]
fn uncertain_mito_del_shift_is_idempotent() {
    let p = provider_homopolymer_named("NC_012920.1");
    assert_idempotent(&p, "NC_012920.1:m.(196del)");
}

mod cds_axis {
    use super::{assert_idempotent, norm};
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::MockProvider;

    /// Single-exon plus-strand coding transcript `ATGGCATAA` (ATG·GCA·TAA):
    /// c.4_6 = `GCA`, whose reverse-complement is `TGC`, so
    /// `c.4_6delinsTGC` canonicalizes to `c.4_6inv` — a deterministic
    /// (shift-free) exercise of the `normalize_cds` wrapper reconstruction.
    fn tx_provider() -> MockProvider {
        let mut p = MockProvider::new();
        let seq = "ATGGCATAA".to_string();
        let len = seq.len() as u64;
        let tx = Transcript::new(
            "NM_TEST.1".to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            seq,
            Some(1),
            Some(len),
            vec![Exon::new(1, 1, len)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::None,
            None,
            None,
        );
        p.add_transcript(tx);
        p
    }

    #[test]
    fn certain_cds_revcomp_delins_canonicalizes_to_inv_unwrapped() {
        let p = tx_provider();
        assert_eq!(norm(&p, "NM_TEST.1:c.4_6delinsTGC"), "NM_TEST.1:c.4_6inv");
    }

    #[test]
    fn uncertain_cds_revcomp_delins_canonicalizes_to_inv_and_stays_wrapped() {
        let p = tx_provider();
        assert_eq!(
            norm(&p, "NM_TEST.1:c.(4_6delinsTGC)"),
            "NM_TEST.1:c.(4_6inv)"
        );
    }

    #[test]
    fn uncertain_cds_revcomp_delins_to_inv_is_idempotent() {
        let p = tx_provider();
        assert_idempotent(&p, "NM_TEST.1:c.(4_6delinsTGC)");
    }
}
