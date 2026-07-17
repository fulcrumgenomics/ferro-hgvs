//! Issue #1052 — a substitution asserting a reference base that disagrees with
//! the loaded reference (`g.8A>C` where g.8 is `T`) must surface a
//! `RefSeqMismatch` warning in lenient mode and reject in strict mode. Before
//! this fix, real substitutions skipped `normalize_na_edit` entirely
//! (`needs_normalization` returned `false` for them), so the reference-base
//! assertion was never validated on any axis.
//!
//! Regression guard (adversarial review): intronic substitutions must stay on
//! their pre-#1052 silent pass-through — no warning, no `ConversionError`, no
//! `ReducedCapabilityNoGenome` — because routing them through the genomic
//! intronic projection would error on variants that used to normalize cleanly.

use ferro_hgvs::{parse_hgvs, FerroError, JsonProvider, NormalizeConfig, Normalizer};
use std::io::Write;

/// Genome-capable provider: one `contig` of `len` cyclic ACGT bytes with
/// `payload` written 1-based at `pos1`. (Mirrors tests/it/issue_1041_repro.rs.)
fn g_provider(contig: &str, len: usize, pos1: usize, payload: &str) -> JsonProvider {
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(len).collect();
    for (i, b) in payload.bytes().enumerate() {
        bases[pos1 - 1 + i] = b;
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

#[test]
fn wrong_ref_genomic_sub_warns_in_lenient() {
    // 300 bp cyclic ACGT: base at 1-based pos 8 is `T` (ACGTACG|T). Assert `A>C`.
    let p = g_provider("c", 300, 8, "T");
    let v = parse_hgvs("c:g.8A>C").unwrap();
    let r = Normalizer::with_config(p, NormalizeConfig::lenient())
        .normalize_with_diagnostics(&v)
        .expect("lenient must not reject");
    assert_eq!(
        r.result.to_string(),
        "c:g.8A>C",
        "output unchanged in lenient"
    );
    assert!(
        r.warnings.iter().any(|w| matches!(
            w,
            ferro_hgvs::normalize::NormalizationWarning::RefSeqMismatch { .. }
        )),
        "expected RefSeqMismatch; got {:?}",
        r.warnings
    );
}

#[test]
fn correct_ref_genomic_sub_no_warning() {
    let p = g_provider("c", 300, 8, "T");
    let v = parse_hgvs("c:g.8T>C").unwrap();
    let r = Normalizer::with_config(p, NormalizeConfig::lenient())
        .normalize_with_diagnostics(&v)
        .expect("lenient");
    assert_eq!(r.result.to_string(), "c:g.8T>C");
    assert!(
        r.warnings.is_empty(),
        "correct ref → no warning; got {:?}",
        r.warnings
    );
}

#[test]
fn wrong_ref_genomic_sub_rejects_in_strict() {
    let p = g_provider("c", 300, 8, "T");
    let v = parse_hgvs("c:g.8A>C").unwrap();
    let err = Normalizer::with_config(p, NormalizeConfig::strict())
        .normalize(&v)
        .expect_err("strict must reject wrong-ref sub");
    assert!(
        matches!(err, FerroError::ReferenceMismatch { .. }),
        "expected ReferenceMismatch, got {err:?}"
    );
}

#[test]
fn uncertain_wrong_ref_genomic_sub_stays_silent_in_lenient_and_strict() {
    // Regression guard (whole-branch review): the uncertain/predicted-
    // wrapped-substitution carve-out was initially added only to
    // normalize_cds/tx/rna, missing normalize_genome — the primary axis
    // this issue targets. `g.(8A>C)` (wrong ref, uncertain) must stay a
    // silent pass-through in EVERY mode, exactly like the c./n./r. cases in
    // `transcript_axis`: no RefSeqMismatch in lenient, no ReferenceMismatch
    // rejection in strict.
    let p = g_provider("c", 300, 8, "T");
    let v = parse_hgvs("c:g.(8A>C)").unwrap();

    let r = Normalizer::with_config(p.clone(), NormalizeConfig::lenient())
        .normalize_with_diagnostics(&v)
        .expect("lenient must not reject an uncertain sub");
    assert_eq!(r.result.to_string(), "c:g.(8A>C)", "unchanged in lenient");
    assert!(
        r.warnings.is_empty(),
        "uncertain sub must stay silent regardless of ref mismatch; got {:?}",
        r.warnings
    );

    let out = Normalizer::with_config(p, NormalizeConfig::strict())
        .normalize(&v)
        .expect("strict mode must not reject an uncertain sub (no ReferenceMismatch)");
    assert_eq!(out.to_string(), "c:g.(8A>C)", "unchanged in strict");
}

#[test]
fn uncertain_correct_ref_genomic_sub_preserves_wrapper() {
    // Silent-data-loss guard: `g.(8T>C)` has a CORRECT stated ref (g.8 is
    // `T`). Before the genome-axis carve-out, this uncertain substitution
    // still reached `normalize_na_edit`'s validation-then-pass-through arm,
    // which returns the edit unchanged — but the CALLING function
    // (`normalize_genome`) then reconstructed the result via
    // `LocEdit::new(...)`, which unconditionally emits `Mu::Certain`,
    // silently dropping the `(...)` predicted-change wrapper even though
    // the ref matched and zero warnings fired. That's silent semantic data
    // loss on entirely valid input. The carve-out returns the ORIGINAL
    // variant, so the wrapper survives.
    let p = g_provider("c", 300, 8, "T");
    let v = parse_hgvs("c:g.(8T>C)").unwrap();
    let r = Normalizer::with_config(p, NormalizeConfig::lenient())
        .normalize_with_diagnostics(&v)
        .expect("lenient must not reject");
    assert_eq!(
        r.result.to_string(),
        "c:g.(8T>C)",
        "the uncertain wrapper must survive normalization"
    );
    assert!(
        r.warnings.is_empty(),
        "correct-ref uncertain sub → no warnings; got {:?}",
        r.warnings
    );
}

#[test]
fn uncertain_wrong_ref_mito_sub_stays_silent() {
    // `m.` axis counterpart (`normalize_mt`), mirroring the `g.` guards
    // above. Mirrors the provider pattern from
    // `tests/it/issue_1044_mito_window_clamp.rs` (a plain `JsonProvider`
    // genomic-sequence entry keyed by the mitochondrial accession).
    let contig = "NC_012920.1";
    let p = g_provider(contig, 300, 8, "T");
    let v = parse_hgvs(&format!("{contig}:m.(8A>C)")).unwrap();

    let r = Normalizer::with_config(p.clone(), NormalizeConfig::lenient())
        .normalize_with_diagnostics(&v)
        .expect("lenient must not reject an uncertain mito sub");
    assert_eq!(
        r.result.to_string(),
        format!("{contig}:m.(8A>C)"),
        "unchanged in lenient"
    );
    assert!(
        r.warnings.is_empty(),
        "uncertain mito sub must stay silent regardless of ref mismatch; got {:?}",
        r.warnings
    );

    let out = Normalizer::with_config(p, NormalizeConfig::strict())
        .normalize(&v)
        .expect("strict mode must not reject an uncertain mito sub (no ReferenceMismatch)");
    assert_eq!(
        out.to_string(),
        format!("{contig}:m.(8A>C)"),
        "unchanged in strict"
    );
}

#[test]
fn wrong_ref_mito_sub_warns_in_lenient() {
    // Certain `m.` counterpart of `wrong_ref_genomic_sub_warns_in_lenient`.
    // A plain (non-uncertain, non-intronic) mitochondrial substitution flows
    // through `normalize_mt`'s window fetch into `normalize_na_edit`, exactly
    // like `g.`, so a mismatched stated ref must surface `RefSeqMismatch`.
    // m.8 is `T` (cyclic ACGT); assert `A>C`.
    let contig = "NC_012920.1";
    let p = g_provider(contig, 300, 8, "T");
    let v = parse_hgvs(&format!("{contig}:m.8A>C")).unwrap();
    let r = Normalizer::with_config(p, NormalizeConfig::lenient())
        .normalize_with_diagnostics(&v)
        .expect("lenient must not reject");
    assert_eq!(
        r.result.to_string(),
        format!("{contig}:m.8A>C"),
        "output unchanged in lenient"
    );
    assert!(
        r.warnings.iter().any(|w| matches!(
            w,
            ferro_hgvs::normalize::NormalizationWarning::RefSeqMismatch { .. }
        )),
        "expected RefSeqMismatch on m. sub; got {:?}",
        r.warnings
    );
}

#[test]
fn wrong_ref_mito_sub_rejects_in_strict() {
    let contig = "NC_012920.1";
    let p = g_provider(contig, 300, 8, "T");
    let v = parse_hgvs(&format!("{contig}:m.8A>C")).unwrap();
    let err = Normalizer::with_config(p, NormalizeConfig::strict())
        .normalize(&v)
        .expect_err("strict must reject a wrong-ref mito sub");
    assert!(
        matches!(err, FerroError::ReferenceMismatch { .. }),
        "expected ReferenceMismatch, got {err:?}"
    );
}

#[test]
fn correct_ref_mito_sub_no_warning() {
    let contig = "NC_012920.1";
    let p = g_provider(contig, 300, 8, "T");
    let v = parse_hgvs(&format!("{contig}:m.8T>C")).unwrap();
    let r = Normalizer::with_config(p, NormalizeConfig::lenient())
        .normalize_with_diagnostics(&v)
        .expect("lenient");
    assert_eq!(r.result.to_string(), format!("{contig}:m.8T>C"));
    assert!(
        r.warnings.is_empty(),
        "correct ref → no warning; got {:?}",
        r.warnings
    );
}

#[test]
fn genomic_sub_without_reference_passes_through() {
    // Provider with a DIFFERENT contig → fetch fails for `nc` → pass-through.
    let p = g_provider("other", 300, 8, "T");
    let v = parse_hgvs("nc:g.8A>C").unwrap();
    let r = Normalizer::with_config(p, NormalizeConfig::lenient())
        .normalize_with_diagnostics(&v)
        .expect("lenient");
    assert_eq!(r.result.to_string(), "nc:g.8A>C");
    assert!(
        r.warnings.is_empty(),
        "no reference → no warning; got {:?}",
        r.warnings
    );
}

mod transcript_axis {
    use super::*;
    use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
    use ferro_hgvs::MockProvider;

    // 60 bp transcript; c.5 = 'A' (ATGCA…). Build a genomic backing so intronic
    // guard tests can run WITH and WITHOUT genomic data.
    fn tx_only() -> MockProvider {
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

    /// 2-exon plus-strand transcript with genomic exon coordinates + a
    /// chromosome, but NO genomic sequence loaded (mirrors
    /// `tests/it/issue_1012_reduced_capability_no_genome.rs`'s
    /// `provider_no_genomic_sequence`). `has_genomic_data()` is false, so
    /// pre-#1052 this shape is exactly what made a former-error intronic path
    /// emit `ReducedCapabilityNoGenome` (#1012) once it was reachable.
    fn intronic_tx_no_genomic_sequence() -> MockProvider {
        let mut p = MockProvider::new();
        let tx = Transcript::new(
            "NM_INTR.1".to_string(),
            Some("INTR".to_string()),
            Strand::Plus,
            Some("AAAATGCCCCGGGGTAGAATAA".to_string()),
            Some(1),
            Some(22),
            vec![
                Exon::with_genomic(1, 1, 10, 100, 109),
                Exon::with_genomic(2, 11, 22, 130, 141),
            ],
            Some("chr_intr".to_string()),
            Some(100),
            Some(141),
            Default::default(),
            ManeStatus::None,
            None,
            None,
        );
        p.add_transcript(tx);
        // Deliberately no add_genomic_sequence — has_genomic_data() == false.
        p
    }

    #[test]
    fn wrong_ref_cds_sub_warns() {
        // c.5 is 'A'; assert 'G>C' (wrong ref).
        let v = parse_hgvs("NM_TEST.1:c.5G>C").unwrap();
        let r = Normalizer::with_config(tx_only(), NormalizeConfig::lenient())
            .normalize_with_diagnostics(&v)
            .expect("lenient");
        assert!(
            r.warnings.iter().any(|w| matches!(
                w,
                ferro_hgvs::normalize::NormalizationWarning::RefSeqMismatch { .. }
            )),
            "expected RefSeqMismatch on c. sub; got {:?}",
            r.warnings
        );
    }

    #[test]
    fn uncertain_wrong_ref_cds_sub_stays_silent_in_lenient() {
        // Regression guard: an uncertain/predicted-wrapped substitution
        // (`c.(5G>C)`) with a MISMATCHED stated ref must stay a silent
        // pass-through — no `RefSeqMismatch` warning — even though the
        // CERTAIN counterpart (`wrong_ref_cds_sub_warns`, same mismatch)
        // does warn. `Mu::inner()` unwraps both `Certain` and `Uncertain`,
        // so without the certainty-gated carve-out this would validate and
        // emit a NEW warning on input that previously normalized cleanly.
        // c.5 is 'A'; stated ref 'G' mismatches.
        let v = parse_hgvs("NM_TEST.1:c.(5G>C)").unwrap();
        let r = Normalizer::with_config(tx_only(), NormalizeConfig::lenient())
            .normalize_with_diagnostics(&v)
            .expect("lenient must not reject an uncertain sub");
        assert_eq!(r.result.to_string(), "NM_TEST.1:c.(5G>C)", "unchanged");
        assert!(
            r.warnings.is_empty(),
            "uncertain sub must stay silent regardless of ref mismatch; got {:?}",
            r.warnings
        );
    }

    #[test]
    fn uncertain_wrong_ref_cds_sub_stays_silent_in_strict() {
        // Strict-mode counterpart: before the certainty gate, this would
        // have escalated to a hard `FerroError::ReferenceMismatch` — a new
        // rejection of input that previously normalized cleanly.
        let v = parse_hgvs("NM_TEST.1:c.(5G>C)").unwrap();
        let r = Normalizer::with_config(tx_only(), NormalizeConfig::strict())
            .normalize(&v)
            .expect("strict mode must not reject an uncertain sub (no ReferenceMismatch)");
        assert_eq!(r.to_string(), "NM_TEST.1:c.(5G>C)", "unchanged");
    }

    #[test]
    fn correct_ref_cds_sub_no_warning() {
        // c.5 is 'A'.
        let v = parse_hgvs("NM_TEST.1:c.5A>C").unwrap();
        let r = Normalizer::with_config(tx_only(), NormalizeConfig::lenient())
            .normalize_with_diagnostics(&v)
            .expect("lenient");
        assert!(
            r.warnings.is_empty(),
            "correct c. ref → no warning; got {:?}",
            r.warnings
        );
    }

    #[test]
    fn correct_ref_exonic_sub_at_transcript_boundary_no_capability_warning() {
        // Regression guard (#1052 follow-up finding 1): a correct-ref
        // substitution at the LAST base of a single-exon transcript
        // (`tx_only()`'s exon spans the whole 60bp transcript, so c.60 sits
        // exactly at `exon_only.right`) used to spuriously trip the #670
        // exon/intron junction 3'-shuffle-continuation check inside
        // `normalize_cds` — a check that only makes sense for shuffle-capable
        // edit kinds (del/dup/etc.), never for a point substitution. Before
        // the post-`normalize_na_edit` short-circuit, this fixture (no
        // genomic sequence loaded) emitted a spurious `ReducedCapabilityNoGenome`
        // warning in lenient mode even though the substitution is entirely
        // valid and needs no shuffle. c.60 is `T` in
        // "ATGC...GGGGGT" — assert the matching ref.
        let v = parse_hgvs("NM_TEST.1:c.60T>C").unwrap();
        let r = Normalizer::with_config(tx_only(), NormalizeConfig::lenient())
            .normalize_with_diagnostics(&v)
            .expect("lenient must not reject a valid boundary substitution");
        assert_eq!(r.result.to_string(), "NM_TEST.1:c.60T>C", "unchanged");
        assert!(
            r.warnings.is_empty(),
            "correct-ref boundary sub must have zero warnings \
             (specifically no ReducedCapabilityNoGenome); got {:?}",
            r.warnings
        );
    }

    #[test]
    fn correct_ref_exonic_sub_at_transcript_boundary_passes_strict() {
        // Strict-mode counterpart: before the fix, the spurious
        // `ReducedCapabilityNoGenome` warning above would have escalated to a
        // hard `FerroError::ReducedReferenceCapability`, rejecting a
        // perfectly valid, in-bounds, correct-ref substitution.
        let v = parse_hgvs("NM_TEST.1:c.60T>C").unwrap();
        Normalizer::with_config(tx_only(), NormalizeConfig::strict())
            .normalize(&v)
            .expect(
                "strict mode must accept a valid boundary substitution \
                 (no ReducedReferenceCapability)",
            );
    }

    #[test]
    fn intronic_sub_bare_nm_stays_on_pre_existing_eintronic_warning_only() {
        // `NM_TEST.1:c.10+2A>G` is a bare (no genomic-context) NM_ accession
        // at an intronic offset. The pre-existing #486 EINTRONIC/W4007 check
        // (`intronic_on_bare_transcript_warning`) runs in `normalize_core`
        // BEFORE dispatch and fires for ANY edit type on this shape — see
        // `tests/it/issue_486_eintronic.rs::
        // lenient_warns_and_keeps_value_for_bare_nm_intronic_substitution`.
        // That warning is orthogonal to #1052. The regression this guards
        // against is specifically the sub ALSO tripping `ReducedCapabilityNoGenome`
        // or a hard error by falling into the genomic intronic dispatch — so
        // assert the warning set is *exactly* the one pre-existing warning,
        // not "any warning at all".
        let v = parse_hgvs("NM_TEST.1:c.10+2A>G").unwrap();
        let r = Normalizer::with_config(tx_only(), NormalizeConfig::lenient())
            .normalize_with_diagnostics(&v)
            .expect("intronic sub must not error in lenient");
        assert_eq!(r.result.to_string(), "NM_TEST.1:c.10+2A>G", "unchanged");
        assert_eq!(
            r.warnings.len(),
            1,
            "expected exactly the pre-existing EINTRONIC warning, got {:?}",
            r.warnings
        );
        assert!(
            matches!(
                &r.warnings[0],
                ferro_hgvs::normalize::NormalizationWarning::IntronicOnBareTranscript { .. }
            ),
            "expected IntronicOnBareTranscript, got {:?}",
            r.warnings
        );
    }

    #[test]
    fn intronic_sub_without_genomic_sequence_is_silent_passthrough() {
        // True regression guard (adversarial review, HIGH): before the #1052
        // intronic guard, a real substitution reaching this genomic-context
        // (`NG_(NM_)`) intronic form — spec-valid, so it does NOT also raise
        // the orthogonal EINTRONIC warning (mirrors
        // `tests/it/issue_1012_reduced_capability_no_genome.rs`'s
        // `intronic_variant_warns_and_degrades_without_genome`, but for a
        // deletion) — would have been routed through
        // `normalize_intronic_cds`/`normalize_boundary_spanning_cds`, which
        // requires genomic sequence data. With none loaded, that would emit
        // `ReducedCapabilityNoGenome`: a NEW warning on a variant that
        // previously normalized cleanly. The guard must keep this silent.
        let v = parse_hgvs("NG_012337.1(NM_INTR.1):c.10+5A>G").unwrap();
        let r = Normalizer::with_config(
            intronic_tx_no_genomic_sequence(),
            NormalizeConfig::lenient(),
        )
        .normalize_with_diagnostics(&v)
        .expect("intronic sub must not error in lenient");
        assert_eq!(
            r.result.to_string(),
            "NG_012337.1(NM_INTR.1):c.10+5A>G",
            "unchanged"
        );
        assert!(
            r.warnings.is_empty(),
            "intronic sub without genomic sequence → zero warnings; got {:?}",
            r.warnings
        );
    }
}
