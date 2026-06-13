//! Audit for issue #214 — normalize-time consistency check for
//! `unit_len × span` and reference-base content on `NaEdit::Repeat` and
//! `NaEdit::MultiRepeat`.
//!
//! Spec: `assets/hgvs-nomenclature/docs/recommendations/DNA/repeated.md`.
//! For `unit[N]` over `[start, end]`:
//!   1. `(end - start + 1) % unit_len == 0` — the span must divide cleanly
//!      into whole units.
//!   2. The reference bases at `[start, end]` must equal `unit` repeated
//!      `(end - start + 1) / unit_len` times. The variant count `N` is the
//!      alt-allele count and is independent of the reference repeat count
//!      (which `normalize_repeat` discovers by scanning); the consistency
//!      check is between the *span* and the *reference bases*.
//!
//! Before this fix the parser was the only gatekeeper and it has no
//! reference handle, so structurally-valid-but-reference-inconsistent
//! descriptions (e.g. `g.100_106AC[3]` — span 7, unit_len 2; or
//! `g.100_105AC[3]` against reference `ACACAT` — span ok, last unit
//! mismatches) sailed through and either silently misnormalized or
//! emitted positions the spec considers undefined.
//!
//! After this fix the check runs inside `validate::validate_reference`
//! and threads through the existing strict/lenient/silent
//! `RefSeqMismatch` flow: strict rejects with
//! `FerroError::ReferenceMismatch`; lenient emits a warning and keeps
//! the input verbatim; silent does not reject (and matches the existing
//! sub/del/dup behavior of still recording the warning in
//! `result.warnings` — `silent` only suppresses the rejection and the
//! stderr-print, not the in-memory record).

use ferro_hgvs::error::FerroError;
use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

// =============================================================================
// Test infrastructure
// =============================================================================

/// 256 bp of `N` flanking the variant tract on each side. The `coverage_gap_tests`
/// audit uses `ACGT`-based padding, but `ACGT` starts with `AC`, which would
/// form a spurious extra `AC` copy at the tract boundary for `AC[N]`-style
/// tests. `N` is inert against byte-equality scans (`count_tandem_repeats`
/// stops at any non-IUPAC-base boundary) and against the new repeat
/// consistency check (which compares the unit*k expansion against the
/// reference span; N ≠ any unit byte, so the boundary never spuriously
/// extends the apparent tract).
const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
const PAD_OFFSET: u64 = 256;

fn padded(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}

/// Build a MockProvider holding a single genomic accession.
fn g_provider(accession: &str, padded_seq: &str) -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence(accession, padded_seq);
    p
}

fn normalize_strict(provider: MockProvider, input: &str) -> Result<String, FerroError> {
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs(input).expect("parse failed");
    normalizer.normalize(&variant).map(|v| format!("{}", v))
}

fn normalize_lenient_warnings(
    provider: MockProvider,
    input: &str,
) -> (String, Vec<NormalizationWarning>) {
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    let variant = parse_hgvs(input).expect("parse failed");
    let r = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient normalize failed");
    (format!("{}", r.result), r.warnings)
}

fn has_refseq_mismatch(warnings: &[NormalizationWarning]) -> bool {
    warnings
        .iter()
        .any(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. }))
}

fn is_ref_mismatch(err: &FerroError) -> bool {
    matches!(err, FerroError::ReferenceMismatch { .. })
}

fn is_eintronic(err: &FerroError) -> bool {
    matches!(err, FerroError::IntronicVariant { .. })
}

// =============================================================================
// SECTION 1 — g. unit_len does not divide span
// =============================================================================

mod g_unit_does_not_divide_span {
    use super::*;

    /// Build a g. provider with reference tract `ACACAC` at HGVS
    /// positions 257..=262 (the PAD is 256 bp of N, so HGVS pos 257 is
    /// the first base of the core).
    fn provider_with_acacac() -> MockProvider {
        let core = "ACACAC";
        let padded_seq = padded(core);
        g_provider("NC_000001.11", &padded_seq)
    }

    /// `g.257_263AC[3]` — span 7, unit_len 2: 7 % 2 != 0. The user typed
    /// one position too many, so the description is internally
    /// inconsistent and must be rejected in strict mode.
    #[test]
    fn span_indivisible_by_unit_len_rejected_strict() {
        let err = normalize_strict(provider_with_acacac(), "NC_000001.11:g.257_263AC[3]")
            .expect_err("strict mode must reject span-not-divisible-by-unit_len");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch error, got {:?}",
            err
        );
    }

    /// Same inconsistency in lenient mode emits a `RefSeqMismatch`
    /// warning and continues without rewriting the description (keeping
    /// the input verbatim is the defensive fallback shared with sub/del
    /// stated-ref mismatches).
    #[test]
    fn span_indivisible_by_unit_len_warns_lenient() {
        let (_output, warnings) =
            normalize_lenient_warnings(provider_with_acacac(), "NC_000001.11:g.257_263AC[3]");
        assert!(
            has_refseq_mismatch(&warnings),
            "expected RefSeqMismatch warning, got warnings={:?}",
            warnings
        );
    }

    /// Silent mode does NOT reject the inconsistent description — the
    /// caller has opted in to "treat the input as authoritative". (Note:
    /// the warning is still recorded in `result.warnings` for
    /// observability, matching the existing sub/del/dup behavior;
    /// `silent` only suppresses the rejection / stderr-print, not the
    /// in-memory warning record.)
    #[test]
    fn span_indivisible_by_unit_len_silent_does_not_reject() {
        let normalizer = Normalizer::with_config(provider_with_acacac(), NormalizeConfig::silent());
        let variant = parse_hgvs("NC_000001.11:g.257_263AC[3]").expect("parse");
        let _ = normalizer
            .normalize(&variant)
            .expect("silent mode must not reject");
    }
}

// =============================================================================
// SECTION 2 — g. span divides but reference bases mismatch unit*k
// =============================================================================

mod g_reference_content_mismatch {
    use super::*;

    /// Reference tract at HGVS positions 257..=262 is `ACACAT` — the
    /// last base is `T` rather than the `C` that `AC[3]` requires.
    fn provider_with_acacat() -> MockProvider {
        let core = "ACACAT";
        let padded_seq = padded(core);
        g_provider("NC_000001.11", &padded_seq)
    }

    /// `g.257_262AC[3]` — span 6, unit_len 2, span / unit_len = 3 ✓
    /// (the divisibility gate passes), but the reference bases at that
    /// span are `ACACAT`, which is not `AC` repeated three times.
    #[test]
    fn content_mismatch_rejected_strict() {
        let err = normalize_strict(provider_with_acacat(), "NC_000001.11:g.257_262AC[3]")
            .expect_err("strict mode must reject ref-content mismatch");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }

    #[test]
    fn content_mismatch_warns_lenient() {
        let (_output, warnings) =
            normalize_lenient_warnings(provider_with_acacat(), "NC_000001.11:g.257_262AC[3]");
        assert!(
            has_refseq_mismatch(&warnings),
            "expected RefSeqMismatch warning, got {:?}",
            warnings
        );
    }

    /// Single-base unit, ref-content mismatch (`A[4]` against `AAAT`):
    /// homopolymer test — exercises the unit_len = 1 path.
    #[test]
    fn homopolymer_content_mismatch_rejected_strict() {
        let core = "AAAT";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let err = normalize_strict(provider, "NC_000001.11:g.257_260A[4]")
            .expect_err("strict mode must reject homopolymer ref-content mismatch");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }
}

// =============================================================================
// SECTION 3 — consistent descriptions still round-trip (regression guard)
// =============================================================================

mod g_consistent_descriptions {
    use super::*;

    /// Reference tract at HGVS positions 257..=262 is `ACACAC`. The
    /// description `g.257_262AC[3]` is internally consistent (span 6 =
    /// unit_len 2 × 3, ref bases = `AC` × 3) — must not be rejected by
    /// the new check.
    /// The B5 contract is: a consistent description (span = unit_len × k,
    /// reference bases = unit × k) must not be rejected by the new
    /// consistency check. The specific output shape (preserve repeat /
    /// rewrite to del / dup / identity) is determined by the existing
    /// `normalize_repeat` logic and is out of scope for this test —
    /// what matters is that strict mode accepts the input rather than
    /// surfacing a `ReferenceMismatch`.
    #[test]
    fn ac3_consistent_with_acacac_strict_accepted() {
        let core = "ACACAC";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        normalize_strict(provider, "NC_000001.11:g.257_262AC[3]")
            .expect("strict mode must accept a consistent AC[3] over ACACAC");
    }

    /// Expansion: variant count (4) > ref count (3). Consistency check
    /// looks at the span vs the reference, not at the variant count, so
    /// this must pass.
    #[test]
    fn expansion_consistent_strict_accepted() {
        let core = "ACACAC";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        normalize_strict(provider, "NC_000001.11:g.257_262AC[4]")
            .expect("strict mode must accept consistent AC[4] expansion over ACACAC");
    }

    /// Contraction: variant count (1) < ref count (3). Same logic.
    #[test]
    fn contraction_consistent_strict_accepted() {
        let core = "ACACAC";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        normalize_strict(provider, "NC_000001.11:g.257_262AC[1]")
            .expect("strict mode must accept consistent AC[1] contraction over ACACAC");
    }

    /// Tri-nucleotide unit (the canonical c. coding repeat shape).
    #[test]
    fn trinucleotide_consistent_strict_accepted() {
        let core = "CAGCAGCAGCAGCAG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        normalize_strict(provider, "NC_000001.11:g.257_271CAG[5]")
            .expect("strict mode must accept consistent CAG[5]");
    }
}

// =============================================================================
// SECTION 4 — MultiRepeat consistency
// =============================================================================

mod g_multirepeat {
    use super::*;

    /// Reference tract `CTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG` is
    /// the spec's _FMR1_-style mixed-repeat example (here adapted to
    /// CTG[2]TTG[1]CTG[11] for a clean 42-base tract on the genome).
    fn provider_with_mixed_ctg_ttg() -> MockProvider {
        // 2 CTG + 1 TTG + 11 CTG = 14 trimers = 42 bp.
        let core = "CTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
        assert_eq!(core.len(), 42, "test setup: mixed tract must be 42 bp");
        let padded_seq = padded(core);
        g_provider("NC_000001.11", &padded_seq)
    }

    /// `g.257_298CTG[2]TTG[1]CTG[11]` — span = 298-257+1 = 42, sum of
    /// unit_len × N = 6 + 3 + 33 = 42 ✓, ref bases match the
    /// concatenation. Must round-trip in strict mode.
    #[test]
    fn consistent_multirepeat_round_trips_strict() {
        let provider = provider_with_mixed_ctg_ttg();
        let out = normalize_strict(provider, "NC_000001.11:g.257_298CTG[2]TTG[1]CTG[11]")
            .expect("must accept consistent multirepeat");
        // No rejection is the assertion. Output shape is otherwise
        // unconstrained — normalize may emit the same multirepeat or
        // collapse adjacent CTG runs.
        assert!(!out.is_empty());
    }

    /// `g.257_280CTG[2]TTG[1]CTG[11]` — span = 24, sum = 42, span
    /// doesn't equal sum: reject in strict.
    #[test]
    fn multirepeat_sum_does_not_match_span_rejected_strict() {
        let provider = provider_with_mixed_ctg_ttg();
        let err = normalize_strict(provider, "NC_000001.11:g.257_280CTG[2]TTG[1]CTG[11]")
            .expect_err("strict mode must reject multirepeat with span≠sum");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }

    /// Span = sum but reference content doesn't match the concatenation
    /// — e.g., declare `CTG[2]TTG[1]CTG[11]` but the reference is
    /// `CTGTTGCTG...` with the order swapped.
    #[test]
    fn multirepeat_content_mismatch_rejected_strict() {
        // Reference: 1 CTG + 1 TTG + 12 CTG = 42 bp. User declares
        // 2 CTG + 1 TTG + 11 CTG: same length, same bases, but the
        // *positions* of the TTG don't line up.
        let core = "CTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
        assert_eq!(core.len(), 42);
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let err = normalize_strict(provider, "NC_000001.11:g.257_298CTG[2]TTG[1]CTG[11]")
            .expect_err("strict mode must reject multirepeat with content mismatch");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }
}

// =============================================================================
// SECTION 5 — no reference data → check skipped (graceful degradation)
// =============================================================================

mod no_reference_data {
    use super::*;

    /// When the provider has no entry for the accession, the validation
    /// gate has nothing to compare against. The shared convention with
    /// the existing stated-ref check is to skip silently rather than
    /// emit a spurious mismatch — exercised here so we lock that the
    /// new check follows the same convention.
    #[test]
    fn missing_accession_does_not_emit_mismatch() {
        let provider = MockProvider::new(); // no genomic sequence added
        let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
        let variant = parse_hgvs("NC_000001.11:g.257_263AC[3]").expect("parse must succeed");
        let result = normalizer.normalize(&variant);
        // The normalize call may surface a "Reference not found" /
        // sequence-fetch error, but it must NOT surface a
        // ReferenceMismatch on an unfetched window. Either Ok or a
        // non-mismatch error is acceptable.
        if let Err(err) = result {
            assert!(
                !is_ref_mismatch(&err),
                "must not emit ReferenceMismatch when reference data is unavailable, got {:?}",
                err
            );
        }
    }
}

// =============================================================================
// SECTION 6 — c. UTR context (codon-frame exempt, content check still applies)
// =============================================================================

mod c_utr {
    use super::*;

    /// Build a single-exon plus-strand transcript with a 5'-UTR
    /// containing a `G[6]` homopolymer. CDS starts at tx position 31.
    ///
    /// Transcript layout (in transcript / coding order, plus strand):
    ///   tx 1..=30   5' UTR — contains `GGGGGG` at tx positions 6..=11
    ///                         (i.e. c.-25..=c.-20)
    ///   tx 31..=60  CDS proper
    ///   tx 61..=90  3' UTR
    ///
    /// Genomic mapping: 1:1 onto the padded genomic sequence starting at
    /// `PAD_OFFSET`. So tx pos 6 ↔ genomic pos (PAD_OFFSET + 5).
    fn provider_with_utr_g6() -> MockProvider {
        let utr5 = "AAAAA".to_string() + "GGGGGG" + "AAAAAAAAAAAAAAAAAAA"; // 5+6+19 = 30
        let cds = "ATGAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string(); // 30 bp, in-frame
        let utr3 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string(); // 30 bp
        let tx_seq = format!("{}{}{}", utr5, cds, utr3);
        assert_eq!(tx_seq.len(), 90);

        // Genomic is identical to transcript (plus strand, single exon).
        let padded_seq = padded(&tx_seq);

        let transcript = Transcript::new(
            "NM_UTR.1".to_string(),
            Some("UTRGENE".to_string()),
            Strand::Plus,
            tx_seq,
            Some(31),
            Some(60),
            vec![Exon::with_genomic(1, 1, 90, PAD_OFFSET, PAD_OFFSET + 89)],
            Some("chr_utr".to_string()),
            Some(PAD_OFFSET),
            Some(PAD_OFFSET + 89),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut p = MockProvider::new();
        p.add_genomic_sequence("chr_utr", padded_seq);
        p.add_transcript(transcript);
        p
    }

    /// 5'UTR is exempt from the codon-frame `unit_len % 3 == 0`
    /// restriction (issue #209 / repeated.md). A consistent `G[6]`
    /// description over the UTR tract must NOT be rejected on
    /// codon-frame grounds; the content check is the only gate.
    #[test]
    fn consistent_utr_g6_round_trips_strict() {
        let provider = provider_with_utr_g6();
        // tx pos 6..=11 ↔ c.-25..=c.-20 (CDS starts at tx 31, so c.1 is
        // tx 31 and c.-1 is tx 30; tx 6 = c.-25).
        let out =
            normalize_strict(provider, "NM_UTR.1:c.-25_-20G[6]").expect("UTR G[6] must round-trip");
        assert!(
            out.contains("G[6]") || out.contains("="),
            "unexpected UTR repeat normalization: {}",
            out
        );
    }

    /// Reference at the same span has a non-`G` base, so the content
    /// check must reject in strict mode.
    #[test]
    fn utr_content_mismatch_rejected_strict() {
        // Replace one of the G's with a T in the UTR tract.
        let utr5 = "AAAAA".to_string() + "GGGTGG" + "AAAAAAAAAAAAAAAAAAA"; // bad
        let cds = "ATGAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string();
        let utr3 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_string();
        let tx_seq = format!("{}{}{}", utr5, cds, utr3);
        let padded_seq = padded(&tx_seq);

        let transcript = Transcript::new(
            "NM_UTR.1".to_string(),
            Some("UTRGENE".to_string()),
            Strand::Plus,
            tx_seq,
            Some(31),
            Some(60),
            vec![Exon::with_genomic(1, 1, 90, PAD_OFFSET, PAD_OFFSET + 89)],
            Some("chr_utr".to_string()),
            Some(PAD_OFFSET),
            Some(PAD_OFFSET + 89),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );

        let mut p = MockProvider::new();
        p.add_genomic_sequence("chr_utr", padded_seq);
        p.add_transcript(transcript);

        let err = normalize_strict(p, "NM_UTR.1:c.-25_-20G[6]")
            .expect_err("UTR G[6] over GGGTGG must be rejected");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }
}

// =============================================================================
// SECTION 7 — c. intronic plus-strand
// =============================================================================
//
// Per repeated.md, the codon-frame `unit_len % 3 == 0` restriction
// applies only to the coding sequence — introns are exempt (issue #209).
// The content check still applies: an inconsistent intronic repeat must
// be rejected.

mod c_intronic_plus_strand {
    use super::*;

    /// Build a plus-strand transcript with an intron containing an
    /// `AT[4]`-shaped tract — `ATATATAT` at intronic positions
    /// c.30+1..=c.30+8.
    ///
    /// Layout (plus strand, two exons + one intron):
    ///   Exon 1: tx 1..=30  ↔ genomic (PAD..=PAD+29)
    ///   Intron 1:           genomic (PAD+30..=PAD+39)   — 10 bp,
    ///                       seeded with `ATATATATGG`
    ///   Exon 2: tx 31..=60 ↔ genomic (PAD+40..=PAD+69)
    fn provider_with_intronic_at4() -> MockProvider {
        let exon1 = "ATGCATGCATGCATGCATGCATGCATGCAT"; // 30 bp
        let intron = "ATATATATGG"; // 10 bp, `AT` tract followed by `GG`
        let exon2 = "AAACAACATGGAAAAAAAAAAAAAAAAAAA"; // 30 bp
        let core = format!("{}{}{}", exon1, intron, exon2);
        let padded_seq = padded(&core);

        let tx_seq = format!("{}{}", exon1, exon2);
        assert_eq!(tx_seq.len(), 60);

        let p_off = PAD_OFFSET;
        let transcript = Transcript::new(
            "NM_INTRON.1".to_string(),
            Some("INTRONGENE".to_string()),
            Strand::Plus,
            tx_seq,
            Some(1),
            Some(60),
            vec![
                Exon::with_genomic(1, 1, 30, p_off, p_off + 29),
                Exon::with_genomic(2, 31, 60, p_off + 40, p_off + 69),
            ],
            Some("chr_intron".to_string()),
            Some(p_off),
            Some(p_off + 69),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );
        let mut p = MockProvider::new();
        p.add_genomic_sequence("chr_intron", padded_seq);
        p.add_transcript(transcript);
        p
    }

    /// Consistent intronic `AT[4]` on the spec-valid genomic-context form:
    /// strict mode must accept (no rejection). Output shape is out of scope.
    #[test]
    fn consistent_intronic_at4_strict_accepted() {
        let provider = provider_with_intronic_at4();
        normalize_strict(provider, "NG_012337.1(NM_INTRON.1):c.30+1_30+8AT[4]")
            .expect("strict mode must accept consistent intronic AT[4] on NG_(NM_)");
    }

    /// Bare-transcript form of the same input is spec-invalid: strict mode
    /// rejects it as EINTRONIC (#486) before any repeat validation runs.
    #[test]
    fn consistent_intronic_at4_bare_rejects_eintronic_strict() {
        let provider = provider_with_intronic_at4();
        let err = normalize_strict(provider, "NM_INTRON.1:c.30+1_30+8AT[4]")
            .expect_err("bare-NM intronic must reject as EINTRONIC in strict mode");
        assert!(
            is_eintronic(&err),
            "expected IntronicVariant, got {:?}",
            err
        );
    }

    /// `c.30+1_30+9AT[4]` — span 9, unit_len 2: doesn't divide. On the
    /// genomic-context form, strict mode rejects with ReferenceMismatch.
    #[test]
    fn intronic_span_indivisible_rejected_strict() {
        let provider = provider_with_intronic_at4();
        let err = normalize_strict(provider, "NG_012337.1(NM_INTRON.1):c.30+1_30+9AT[4]")
            .expect_err("intronic span-not-divisible-by-unit_len must be rejected");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }

    /// `c.30+3_30+10AT[4]` — span 8 ✓ but reference is `ATATATGG`, not
    /// `ATATATAT`. Content mismatch on the genomic-context form.
    #[test]
    fn intronic_content_mismatch_rejected_strict() {
        let provider = provider_with_intronic_at4();
        let err = normalize_strict(provider, "NG_012337.1(NM_INTRON.1):c.30+3_30+10AT[4]")
            .expect_err("intronic content mismatch must be rejected");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }
}

// =============================================================================
// SECTION 8 — c. intronic minus-strand
// =============================================================================

mod c_intronic_minus_strand {
    use super::*;

    /// Plus-strand intron contains `CCCCATATAT` at genomic (PAD+30..=PAD+39).
    /// Reverse-complemented to transcript view, the intron looks like
    /// `ATATATGGGG` — so transcript-view AT[3] sits at intronic positions
    /// c.30+1..=c.30+6 (note: minus-strand intronic numbering puts +1 at
    /// the 3' end of the genomic intron when read in transcript order).
    fn provider_with_minus_intronic_at3() -> MockProvider {
        // Plus-strand genomic (read 5' to 3'):
        //   exon1 RC (30) | intron (10) | exon2 RC (30)
        // Transcript reads the reverse complement of all of this.
        let intron = "CCCCATATAT"; // plus-strand; tx view = RC = ATATATGGGG
        let exon1_rc = "ATGCATGCATGCATGCATGCATGCATGCAT"; // arbitrary 30 bp
        let exon2_rc = "GCATGCATGCATGCATGCATGCATGCATGC"; // arbitrary 30 bp
        let core = format!("{}{}{}", exon2_rc, intron, exon1_rc);
        let padded_seq = padded(&core);

        // Transcript view: RC of each exon. We don't care about exact tx
        // content for these tests beyond the intronic tract — these
        // tests check intronic positions which are derived from the
        // intronic genomic bases anyway.
        let tx_e1 = "ATGCATGCATGCATGCATGCATGCATGCAT"; // 30
        let tx_e2 = "GCATGCATGCATGCATGCATGCATGCATGC"; // 30
        let tx_seq = format!("{}{}", tx_e1, tx_e2);
        assert_eq!(tx_seq.len(), 60);

        let p_off = PAD_OFFSET;
        let transcript = Transcript::new(
            "NM_MINTRON.1".to_string(),
            Some("MINTRONGENE".to_string()),
            Strand::Minus,
            tx_seq,
            Some(1),
            Some(60),
            vec![
                Exon::with_genomic(1, 1, 30, p_off + 40, p_off + 69),
                Exon::with_genomic(2, 31, 60, p_off, p_off + 29),
            ],
            Some("chr_mintron".to_string()),
            Some(p_off),
            Some(p_off + 69),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );
        let mut p = MockProvider::new();
        p.add_genomic_sequence("chr_mintron", padded_seq);
        p.add_transcript(transcript);
        p
    }

    /// `c.30+1_30+9AT[4]` — span 9 ≠ 4×2, inconsistent regardless of strand.
    /// On the genomic-context form, strict mode rejects with ReferenceMismatch.
    #[test]
    fn minus_strand_intronic_span_indivisible_rejected_strict() {
        let provider = provider_with_minus_intronic_at3();
        let err = normalize_strict(provider, "NG_012337.1(NM_MINTRON.1):c.30+1_30+9AT[4]")
            .expect_err("minus-strand intronic span-indivisible must be rejected");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }

    /// Bare-transcript form rejects as EINTRONIC (#486) before strand-aware
    /// repeat validation runs.
    #[test]
    fn minus_strand_intronic_bare_rejects_eintronic_strict() {
        let provider = provider_with_minus_intronic_at3();
        let err = normalize_strict(provider, "NM_MINTRON.1:c.30+1_30+9AT[4]")
            .expect_err("bare-NM minus-strand intronic must reject as EINTRONIC");
        assert!(
            is_eintronic(&err),
            "expected IntronicVariant, got {:?}",
            err
        );
    }
}

// =============================================================================
// SECTION 8b — VEP `trailing` + genotype `additional_counts` skip-list
// =============================================================================
//
// `normalize_na_edit`'s repeat arm bails out when `trailing.is_some()` or
// `additional_counts` is non-empty (VEP-style `c.212-18[3]T` and
// genotype notation `A[6][1]`). The new validation gate must mirror
// that skip-list — running the strict `span × unit_len` check on a
// span that includes a trailing base or a genotype-only second count
// would emit a spurious `ReferenceMismatch` while the normalizer
// declines to rewrite the shape.

mod vep_and_genotype_skip {
    use super::*;

    /// VEP-style trailing base: the variant's span includes the trailing
    /// `T` past the `[3]` count. Reference at HGVS 257..=263 is
    /// `ACACACT` (3 AC copies + trailing T). The naive divisibility
    /// check would fail (span 7 vs unit_len 2), but the gate must skip.
    #[test]
    fn vep_style_trailing_does_not_emit_mismatch_strict() {
        let core = "ACACACT";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        // `g.257_263AC[3]T` — VEP trailing-base form. Must not be
        // rejected by the new validation gate.
        let _ = normalize_strict(provider, "NC_000001.11:g.257_263AC[3]T")
            .expect("VEP-style trailing-base repeat must not be rejected");
    }

    /// Genotype notation: `A[6][1]` declares two allele counts. The
    /// span/unit/ref check is undefined when additional counts are
    /// present (the spec doesn't specify which count drives the span)
    /// — validation must skip per the shared skip-list.
    #[test]
    fn genotype_additional_counts_does_not_emit_mismatch_strict() {
        let core = "AAAAAA";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let _ = normalize_strict(provider, "NC_000001.11:g.257_262A[6][1]")
            .expect("genotype additional_counts must not be rejected by the consistency check");
    }
}

// =============================================================================
// SECTION 9 — n. (non-coding transcript) consistency
// =============================================================================

mod n_noncoding {
    use super::*;

    /// Plain non-coding transcript (no CDS). All positions are "n.X";
    /// codon-frame gate never applies. Content check must still gate.
    fn provider_noncoding() -> MockProvider {
        // 60 bp transcript with a 6-bp AC tract starting at n.10.
        // bases 1..=9   = AAAAAAAA + A
        // bases 10..=15 = ACACAC
        // bases 16..=60 = filler
        let tx_seq = "AAAAAAAAA\
                      ACACAC\
                      AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            .to_string();
        assert_eq!(tx_seq.len(), 60);

        let p_off = PAD_OFFSET;
        let padded_seq = padded(&tx_seq);

        let transcript = Transcript::new(
            "NR_NONCODING.1".to_string(),
            Some("NRGENE".to_string()),
            Strand::Plus,
            tx_seq,
            None,
            None,
            vec![Exon::with_genomic(1, 1, 60, p_off, p_off + 59)],
            Some("chr_nr".to_string()),
            Some(p_off),
            Some(p_off + 59),
            GenomeBuild::GRCh38,
            ManeStatus::None,
            None,
            None,
        );
        let mut p = MockProvider::new();
        p.add_genomic_sequence("chr_nr", padded_seq);
        p.add_transcript(transcript);
        p
    }

    /// `n.10_15AC[3]` consistent with reference — round-trip.
    #[test]
    fn n_consistent_round_trips_strict() {
        let provider = provider_noncoding();
        let out = normalize_strict(provider, "NR_NONCODING.1:n.10_15AC[3]").expect("must accept");
        assert!(
            out.contains("AC[3]") || out.contains("="),
            "unexpected n. output: {}",
            out
        );
    }

    /// `n.10_16AC[3]` — span 7, unit_len 2: must reject.
    #[test]
    fn n_span_indivisible_rejected_strict() {
        let provider = provider_noncoding();
        let err = normalize_strict(provider, "NR_NONCODING.1:n.10_16AC[3]")
            .expect_err("n. span-indivisible must be rejected");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }
}
