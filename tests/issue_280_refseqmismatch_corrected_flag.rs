//! Audit for issue #280 — `RefSeqMismatch.corrected` flag must reflect
//! whether the normalizer actually rewrote the stated bases.
//!
//! PR #215 introduced a `RefSeqMismatch` warning that hardcoded
//! `corrected: true` to match the existing sub/del/dup convention.
//! But the `Repeat` consistency mismatches that #215 added pass the
//! description through verbatim — the normalizer cannot and does not
//! rewrite them — so the flag was lying for that case.
//!
//! After the fix the flag is honest: `true` for the sub/del/dup/inv
//! paths whose canonical Display drops the stated bases, and `false`
//! for the `Repeat` / `MultiRepeat` consistency paths that pass the
//! input through verbatim.
//!
//! Coverage — each section drives `validate::validate_reference` to an
//! invalid result, then asserts the `corrected` flag on the emitted
//! `RefSeqMismatch` warning:
//! - SECTION 1: `del` with mismatched stated-ref → `corrected: true`
//!   (canonical Display drops the explicit `sequence`).
//! - SECTION 2: `dup` with mismatched stated-ref → `corrected: true`
//!   (canonical Display drops the explicit `sequence`).
//! - SECTION 3: `inv` with mismatched stated-ref → `corrected: true`
//!   (the Inversion arm in `normalize_na_edit` emits `sequence: None`).
//! - SECTION 4: `Repeat` (single-repeat) span/content mismatch from
//!   PR #215 → `corrected: false` (passes through verbatim). This is
//!   the case that motivated this issue.
//! - SECTION 5: `MultiRepeat` consistency mismatch → `corrected: false`
//!   (also passes through verbatim).

use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

// =============================================================================
// Test infrastructure
// =============================================================================

/// 256 bp of `N` flanking the variant tract on each side. Inert against
/// byte-equality scans (any non-IUPAC base stops `count_tandem_repeats`)
/// and against the repeat consistency check (`N` ≠ any unit byte, so
/// the boundary never spuriously extends the apparent tract).
const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

fn padded(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}

fn g_provider(accession: &str, padded_seq: &str) -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence(accession, padded_seq);
    p
}

/// Build a transcript provider holding a single 60bp transcript whose
/// sequence is `ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT`.
///
/// Positions inside the `CCCCC` homopolymer at c.10..=c.14 exercise the
/// dup mismatched-stated-ref path. Position c.5 is `A`, so `c.5delG`
/// exercises the del mismatched-stated-ref path. Position c.20..=c.22
/// covers `TTT`, so `c.20_22invAAA` exercises the inv mismatched-
/// stated-ref path.
fn tx_provider() -> MockProvider {
    let mut p = MockProvider::new();
    let seq = "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
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

/// Return the `corrected` flag of the single `RefSeqMismatch` warning,
/// panicking if zero or multiple are present.
fn corrected_flag(warnings: &[NormalizationWarning]) -> bool {
    let mismatches: Vec<bool> = warnings
        .iter()
        .filter_map(|w| match w {
            NormalizationWarning::RefSeqMismatch { corrected, .. } => Some(*corrected),
            _ => None,
        })
        .collect();
    assert_eq!(
        mismatches.len(),
        1,
        "expected exactly one RefSeqMismatch, got warnings={:?}",
        warnings
    );
    mismatches[0]
}

// =============================================================================
// SECTION 1 — Deletion with mismatched stated-ref → corrected: true
// =============================================================================
//
// `c.5delG` against a transcript where c.5 = `A` triggers
// `validate_sequence` → mismatch. The canonical Display drops the
// explicit `sequence` (via `canonicalize_edit`), so the stated bases
// never reach the output — the flag is honest at `true`.

mod del_mismatched_stated_ref {
    use super::*;

    #[test]
    fn del_mismatched_stated_ref_corrected_true() {
        let normalizer = Normalizer::with_config(tx_provider(), NormalizeConfig::lenient());
        let v = parse_hgvs("NM_TEST.1:c.5delG").expect("parse");
        let r = normalizer
            .normalize_with_warnings(&v)
            .expect("normalize must not reject in lenient");
        assert!(
            corrected_flag(&r.warnings),
            "del with mismatched stated-ref must report corrected=true; output={}, warnings={:?}",
            r.result,
            r.warnings
        );
    }
}

// =============================================================================
// SECTION 2 — Duplication with mismatched stated-ref → corrected: true
// =============================================================================
//
// `c.10dupA` against a `CCCCC` homopolymer at c.10 triggers
// `validate_sequence` → mismatch. Per issue #219 the dup arm re-reads
// duplicated bases from reference; the canonical Display also drops the
// explicit `sequence` regardless. The flag is honest at `true`.

mod dup_mismatched_stated_ref {
    use super::*;

    #[test]
    fn dup_mismatched_stated_ref_corrected_true() {
        let normalizer = Normalizer::with_config(tx_provider(), NormalizeConfig::lenient());
        let v = parse_hgvs("NM_TEST.1:c.10dupA").expect("parse");
        let r = normalizer
            .normalize_with_warnings(&v)
            .expect("normalize must not reject in lenient");
        assert!(
            corrected_flag(&r.warnings),
            "dup with mismatched stated-ref must report corrected=true; output={}, warnings={:?}",
            r.result,
            r.warnings
        );
    }
}

// =============================================================================
// SECTION 3 — Inversion with mismatched stated-ref → corrected: true
// =============================================================================
//
// `c.20_22invAAA` against a `TTT` tract triggers `validate_sequence` →
// mismatch. The Inversion arm in `normalize_na_edit` emits
// `NaEdit::Inversion { sequence: None, length: None }`, so the stated
// bases are dropped. The flag is honest at `true`.

mod inv_mismatched_stated_ref {
    use super::*;

    #[test]
    fn inv_mismatched_stated_ref_corrected_true() {
        let normalizer = Normalizer::with_config(tx_provider(), NormalizeConfig::lenient());
        let v = parse_hgvs("NM_TEST.1:c.20_22invAAA").expect("parse");
        let r = normalizer
            .normalize_with_warnings(&v)
            .expect("normalize must not reject in lenient");
        assert!(
            corrected_flag(&r.warnings),
            "inv with mismatched stated-ref must report corrected=true; output={}, warnings={:?}",
            r.result,
            r.warnings
        );
    }
}

// =============================================================================
// SECTION 4 — Repeat consistency mismatch → corrected: false
// =============================================================================
//
// This is the case that motivated the issue. PR #215 added a
// span/content consistency check on `NaEdit::Repeat` that emits a
// `RefSeqMismatch` warning. In lenient mode the description is passed
// through verbatim (the normalizer cannot rewrite it without violating
// the user's intent), so the flag must be `false`.

mod repeat_single_consistency_mismatch {
    use super::*;

    /// `g.257_263AC[3]` against `ACACAC` (span 7, unit_len 2 → not
    /// divisible). The validator emits a mismatch; the normalizer
    /// declines to rewrite. The flag must be `false`.
    #[test]
    fn span_indivisible_corrected_false() {
        let core = "ACACAC";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
        let v = parse_hgvs("NC_000001.11:g.257_263AC[3]").expect("parse");
        let r = normalizer
            .normalize_with_warnings(&v)
            .expect("normalize must not reject in lenient");
        assert!(
            !corrected_flag(&r.warnings),
            "Repeat span-indivisible passes through verbatim; corrected must be false; \
             output={}, warnings={:?}",
            r.result,
            r.warnings
        );
    }

    /// `g.257_262AC[3]` against `ACACAT` (span 6 ✓ but ref content
    /// mismatches the unit×k expansion). The validator emits a
    /// mismatch; the normalizer may still rewrite the form (e.g.
    /// emit `dup`), but the original *stated repeat description* is
    /// the input; the warning is about that input being inconsistent
    /// with the reference span. The fix is "the normalizer did not
    /// correct the stated repeat description" → corrected: false.
    #[test]
    fn content_mismatch_corrected_false() {
        let core = "ACACAT";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
        let v = parse_hgvs("NC_000001.11:g.257_262AC[3]").expect("parse");
        let r = normalizer
            .normalize_with_warnings(&v)
            .expect("normalize must not reject in lenient");
        assert!(
            !corrected_flag(&r.warnings),
            "Repeat ref-content mismatch passes through verbatim; corrected must be false; \
             output={}, warnings={:?}",
            r.result,
            r.warnings
        );
    }

    /// Single-base unit homopolymer ref-content mismatch — exercises
    /// the `unit_len == 1` path of the new check.
    #[test]
    fn homopolymer_content_mismatch_corrected_false() {
        let core = "AAAT";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
        let v = parse_hgvs("NC_000001.11:g.257_260A[4]").expect("parse");
        let r = normalizer
            .normalize_with_warnings(&v)
            .expect("normalize must not reject in lenient");
        assert!(
            !corrected_flag(&r.warnings),
            "Repeat homopolymer content-mismatch passes through verbatim; corrected must be false; \
             output={}, warnings={:?}",
            r.result,
            r.warnings
        );
    }
}

// =============================================================================
// SECTION 5 — MultiRepeat consistency mismatch → corrected: false
// =============================================================================
//
// `validate::validate_reference`'s `MultiRepeat` arm emits a mismatch
// when the declared units disagree with the reference span (sum ≠ span,
// or order disagrees). In lenient mode the normalizer passes the
// verbatim multi-unit description through, so the flag must be `false`.

mod multirepeat_consistency_mismatch {
    use super::*;

    /// Sum (42) ≠ span (24).
    #[test]
    fn multirepeat_span_mismatch_corrected_false() {
        // 2 CTG + 1 TTG + 11 CTG = 14 trimers = 42 bp.
        let core = "CTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
        assert_eq!(core.len(), 42);
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
        let v = parse_hgvs("NC_000001.11:g.257_280CTG[2]TTG[1]CTG[11]").expect("parse");
        let r = normalizer
            .normalize_with_warnings(&v)
            .expect("normalize must not reject in lenient");
        assert!(
            !corrected_flag(&r.warnings),
            "MultiRepeat span-mismatch passes through verbatim; corrected must be false; \
             output={}, warnings={:?}",
            r.result,
            r.warnings
        );
    }

    /// Sum = span but the per-unit ordering disagrees with the reference.
    #[test]
    fn multirepeat_content_mismatch_corrected_false() {
        // Reference: 1 CTG + 1 TTG + 12 CTG = 42 bp. User declares
        // 2 CTG + 1 TTG + 11 CTG: same length, same bases, but
        // *positions* of the TTG don't line up.
        let core = "CTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
        assert_eq!(core.len(), 42);
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
        let v = parse_hgvs("NC_000001.11:g.257_298CTG[2]TTG[1]CTG[11]").expect("parse");
        let r = normalizer
            .normalize_with_warnings(&v)
            .expect("normalize must not reject in lenient");
        assert!(
            !corrected_flag(&r.warnings),
            "MultiRepeat content-mismatch passes through verbatim; corrected must be false; \
             output={}, warnings={:?}",
            r.result,
            r.warnings
        );
    }
}
