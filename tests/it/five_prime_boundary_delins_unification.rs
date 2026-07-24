//! Regression tests for the 5'/3' homopolymer-boundary insertion/dup → delins
//! unification (idempotency-campaign follow-up to PR #1169; resolves the #387
//! CDS-end saturation latent and the dup-path non-confluence).
//!
//! At a CDS-interior homopolymer that saturates against a coding boundary, every
//! spelling of the same haplotype — whether written as a dup or a direct
//! insertion, at any start position — must normalize to ONE valid, idempotent
//! `c.<cds_boundary>delins…` form. ferro's design forbids re-axing a CDS-interior
//! edit onto the 3'UTR (`c.*N`) axis (see the CDS-start/-end clamps in
//! `src/normalize/mod.rs`), so the canonical boundary form is the delins, not a
//! `c.<cds_end>_*1ins…` spanning insertion. Spanning *duplications*
//! (`c.7_*1dup`) remain the spec-canonical form and are preserved.
//!
//! Before this fix the dup path returned the 3'-anchored `c.7_*1insAA` verbatim
//! (bypassing both boundary clamps), the direct-insertion path was non-confluent
//! (`c.1delinsAAA` / `c.1_3delinsAAAAA` / `c.1_4delinsAAAAAA`), and every 5' form
//! re-normalized to the INVALID single-position `c.1insAA` (which ferro's own
//! parser rejects).

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// Normalize `input` against a plus-strand CDS whose transcript is `core`, with
/// the CDS spanning `cds_start..=cds_end` (1-based), in shuffle `dir`.
fn norm(core: &str, cds_start: u64, cds_end: u64, input: &str, dir: ShuffleDirection) -> String {
    let normalizer = Normalizer::with_config(
        SyntheticBuilder::cds(core, cds_start, cds_end, Strand::Plus).build(),
        NormalizeConfig::default().with_direction(dir),
    );
    normalizer
        .normalize(&parse_hgvs(input).expect("parse input"))
        .expect("normalize")
        .to_string()
}

/// Assert `s` is valid HGVS that ferro can re-parse (guards against emitting the
/// invalid single-position `c.1insAA`).
fn assert_reparses(s: &str) {
    parse_hgvs(s).unwrap_or_else(|e| panic!("normalizer emitted un-parseable HGVS {s:?}: {e}"));
}

const A7: &str = "AAAAAAA"; // 7-A CDS, cds 1..=7, no UTR either side.

// ---- 5' direction: all spellings converge to c.1delinsAAA -------------------

#[test]
fn five_prime_boundary_dup_becomes_start_delins() {
    for input in ["c.1_2dup", "c.3_4dup", "c.6_7dup"] {
        let out = norm(
            A7,
            1,
            7,
            &format!("NM_TEST.1:{input}"),
            ShuffleDirection::FivePrime,
        );
        assert_eq!(out, "NM_TEST.1:c.1delinsAAA", "5' {input}");
    }
}

#[test]
fn five_prime_boundary_insertion_becomes_start_delins() {
    for input in ["c.1_2insAA", "c.3_4insAA", "c.5_6insAA", "c.6_7insAA"] {
        let out = norm(
            A7,
            1,
            7,
            &format!("NM_TEST.1:{input}"),
            ShuffleDirection::FivePrime,
        );
        assert_eq!(out, "NM_TEST.1:c.1delinsAAA", "5' {input}");
    }
}

#[test]
fn five_prime_boundary_delins_is_idempotent() {
    let once = norm(
        A7,
        1,
        7,
        "NM_TEST.1:c.1delinsAAA",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(
        once, "NM_TEST.1:c.1delinsAAA",
        "c.1delinsAAA must be a 5' fixed point"
    );
}

#[test]
fn five_prime_boundary_output_always_reparses() {
    for input in ["c.1_2dup", "c.6_7dup", "c.1_2insAA", "c.6_7insAA"] {
        let out = norm(
            A7,
            1,
            7,
            &format!("NM_TEST.1:{input}"),
            ShuffleDirection::FivePrime,
        );
        assert_reparses(&out);
    }
}

// ---- 3' direction: all spellings converge to c.7delinsAAA -------------------

#[test]
fn three_prime_boundary_dup_becomes_end_delins() {
    for input in ["c.1_2dup", "c.3_4dup", "c.6_7dup"] {
        let out = norm(
            A7,
            1,
            7,
            &format!("NM_TEST.1:{input}"),
            ShuffleDirection::ThreePrime,
        );
        assert_eq!(out, "NM_TEST.1:c.7delinsAAA", "3' {input}");
    }
}

#[test]
fn three_prime_boundary_insertion_becomes_end_delins() {
    for input in ["c.1_2insAA", "c.3_4insAA", "c.5_6insAA", "c.6_7insAA"] {
        let out = norm(
            A7,
            1,
            7,
            &format!("NM_TEST.1:{input}"),
            ShuffleDirection::ThreePrime,
        );
        assert_eq!(out, "NM_TEST.1:c.7delinsAAA", "3' {input}");
    }
}

#[test]
fn three_prime_boundary_delins_is_idempotent() {
    let once = norm(
        A7,
        1,
        7,
        "NM_TEST.1:c.7delinsAAA",
        ShuffleDirection::ThreePrime,
    );
    assert_eq!(
        once, "NM_TEST.1:c.7delinsAAA",
        "c.7delinsAAA must be a 3' fixed point"
    );
}

// ---- WITH a real 3'UTR: spanning dups preserved, not over-clamped -----------

/// `AAAAAAAGCG`: A-run CDS c.1..=7, 3'UTR = `GCG` (so c.*1 = 'G' EXISTS).
/// A homopolymer insertion at the boundary still clamps to the CDS-end delins
/// (CDS-interior edit stays on the CDS axis)…
#[test]
fn with_3utr_homopolymer_insertion_clamps_to_end_delins() {
    let out = norm(
        "AAAAAAAGCG",
        1,
        7,
        "NM_TEST.1:c.6_7insAA",
        ShuffleDirection::ThreePrime,
    );
    assert_eq!(out, "NM_TEST.1:c.7delinsAAA");
}

/// …but a genuine boundary-spanning *duplication* (`AG` = dup of c.7_*1) is the
/// spec-canonical form and must be PRESERVED, never clamped to a delins.
#[test]
fn with_3utr_spanning_duplication_is_preserved() {
    let out = norm(
        "AAAAAAAGCG",
        1,
        7,
        "NM_TEST.1:c.6_7insAG",
        ShuffleDirection::ThreePrime,
    );
    assert_eq!(out, "NM_TEST.1:c.7_*1dup");
}
