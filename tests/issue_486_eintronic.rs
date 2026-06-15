//! #486 EINTRONIC — reject intronic offsets on a bare transcript reference.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
const PAD_OFFSET: u64 = 256;

fn padded(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}

/// Plus-strand 2-exon coding transcript with a 10-bp intron seeded
/// `ATATATATGG` at intronic positions c.30+1..=c.30+10, with full genomic
/// sequence so indels resolve (mirrors issue_214's provider_with_intronic_at4).
fn provider_intronic_nm() -> MockProvider {
    let exon1 = "ATGCATGCATGCATGCATGCATGCATGCAT"; // 30 bp
    let intron = "ATATATATGG"; // 10 bp
    let exon2 = "AAACAACATGGAAAAAAAAAAAAAAAAAAA"; // 30 bp
    let core = format!("{}{}{}", exon1, intron, exon2);
    let tx_seq = format!("{}{}", exon1, exon2);
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
    p.add_genomic_sequence("chr_intron", padded(&core));
    p.add_transcript(transcript);
    p
}

fn normalize(
    provider: MockProvider,
    cfg: NormalizeConfig,
    input: &str,
) -> Result<String, FerroError> {
    let normalizer = Normalizer::with_config(provider, cfg);
    let variant = parse_hgvs(input).expect("parse failed");
    normalizer.normalize(&variant).map(|v| format!("{}", v))
}

#[test]
fn strict_rejects_bare_nm_intronic_substitution() {
    let err = normalize(
        provider_intronic_nm(),
        NormalizeConfig::strict(),
        "NM_INTRON.1:c.30+1A>G",
    )
    .expect_err("strict mode must reject a bare-NM intronic substitution as EINTRONIC");
    assert!(
        matches!(err, FerroError::IntronicVariant { .. }),
        "expected IntronicVariant (EINTRONIC), got {:?}",
        err
    );
}

#[test]
fn intronic_variant_field_stays_raw_hgvs_clarifier_in_detail() {
    // The spec-form bare-transcript rejection must keep `variant` the raw HGVS
    // string (machine-readable — reparseable, comparable) and carry the
    // human-facing clarifier in `detail` / the rendered message instead (#681).
    let input = "NM_INTRON.1:c.30+1A>G";
    let err = normalize(provider_intronic_nm(), NormalizeConfig::strict(), input)
        .expect_err("strict mode must reject a bare-NM intronic substitution as EINTRONIC");

    let FerroError::IntronicVariant { variant, detail } = &err else {
        panic!("expected IntronicVariant, got {:?}", err);
    };

    // `variant` must be machine-readable HGVS: it reparses, and it round-trips
    // to the raw input (so it carries no prose clarifier).
    assert!(
        parse_hgvs(variant).is_ok(),
        "`variant` must be a reparseable HGVS string, got {variant:?}"
    );
    assert_eq!(
        variant, input,
        "`variant` must be the raw input HGVS, got {variant:?}"
    );

    // The clarifier lives in `detail` and still reaches users via Display, in
    // parentheses after the raw variant.
    let detail = detail
        .as_deref()
        .expect("spec-form rejection should carry a `detail` clarifier");
    assert!(
        detail.contains("a genomic reference is required"),
        "`detail` should explain why, got {detail:?}"
    );
    let rendered = err.to_string();
    assert!(
        rendered.contains("a genomic reference is required") && rendered.contains('('),
        "rendered error should append the clarifier in parentheses, got {rendered:?}"
    );
}

/// 2-exon plus-strand non-coding transcript (NR_, no CDS) with a 20-bp
/// intron-1; genomic coords only (substitutions never reach the
/// genomic-sequence-dependent intronic pass, so no genomic sequence needed).
fn provider_intronic_nr() -> MockProvider {
    let mut provider = MockProvider::new();
    let transcript = Transcript::new(
        "NR_INT.1".to_string(),
        Some("NRINT".to_string()),
        Strand::Plus,
        Some("AAAATGCCCCGGGGTAGAATAA".to_string()),
        None,
        None,
        vec![
            Exon::with_genomic(1, 1, 10, 100, 109),
            Exon::with_genomic(2, 11, 22, 130, 141),
        ],
        Some("chr_int".to_string()),
        Some(100),
        Some(141),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

#[test]
fn strict_rejects_bare_nr_intronic_substitution() {
    let err = normalize(
        provider_intronic_nr(),
        NormalizeConfig::strict(),
        "NR_INT.1:n.10+5C>A",
    )
    .expect_err("strict mode must reject a bare-NR intronic substitution as EINTRONIC");
    assert!(
        matches!(err, FerroError::IntronicVariant { .. }),
        "expected IntronicVariant (EINTRONIC), got {:?}",
        err
    );
}

#[test]
fn lenient_warns_and_keeps_value_for_bare_nm_intronic_substitution() {
    let normalizer = Normalizer::with_config(provider_intronic_nm(), NormalizeConfig::lenient());
    let variant = parse_hgvs("NM_INTRON.1:c.30+1A>G").expect("parse");
    // Lenient must NOT reject.
    let out = normalizer
        .normalize(&variant)
        .expect("lenient must not reject");
    // Output is unchanged from the lenient pre-change behavior (a substitution
    // is echoed; the EINTRONIC change only adds a warning).
    assert_eq!(out.to_string(), "NM_INTRON.1:c.30+1A>G");
    // The W4007 warning IS surfaced via normalize_with_diagnostics.
    let diag = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient diagnostics");
    assert!(
        diag.warnings
            .iter()
            .any(|w| w.code() == "INTRONIC_ON_BARE_TRANSCRIPT"),
        "expected INTRONIC_ON_BARE_TRANSCRIPT warning, got {:?}",
        diag.warnings
    );
}

#[test]
fn ng_nm_intronic_substitution_not_rejected_strict() {
    // Genomic-context form is the spec-valid description: EINTRONIC must NOT
    // fire even in strict mode.
    let out = normalize(
        provider_intronic_nm(),
        NormalizeConfig::strict(),
        "NG_012337.1(NM_INTRON.1):c.30+1A>G",
    )
    .expect("NG_(NM_) intronic must be accepted in strict mode");
    assert!(out.contains("NM_INTRON.1"));
}

#[test]
fn bare_nm_exonic_substitution_unaffected_strict() {
    // An exonic position on a bare NM_ is fine — no EINTRONIC.
    normalize(
        provider_intronic_nm(),
        NormalizeConfig::strict(),
        "NM_INTRON.1:c.5C>A",
    )
    .expect("exonic substitution on a bare NM must not be rejected as EINTRONIC");
}

#[test]
fn strict_rejects_bare_nm_intronic_deletion() {
    // Indel path: del reaches the intronic resolution pass; with genomic
    // sequence present it resolves, then the W4007 warning escalates.
    let err = normalize(
        provider_intronic_nm(),
        NormalizeConfig::strict(),
        "NM_INTRON.1:c.30+2_30+3del",
    )
    .expect_err("strict mode must reject a bare-NM intronic deletion as EINTRONIC");
    assert!(
        matches!(err, FerroError::IntronicVariant { .. }),
        "expected IntronicVariant (EINTRONIC), got {:?}",
        err
    );
}

#[test]
fn strict_rejects_cis_allele_with_intronic_member() {
    // The intronic member inside a cis allele must trigger the reject
    // (normalize_allele recurses through normalize_with_diagnostics per member).
    let err = normalize(
        provider_intronic_nm(),
        NormalizeConfig::strict(),
        "NM_INTRON.1:c.[5C>A;30+1A>G]",
    )
    .expect_err("strict mode must reject a cis allele containing a bare-NM intronic member");
    assert!(
        matches!(err, FerroError::IntronicVariant { .. }),
        "expected IntronicVariant (EINTRONIC), got {:?}",
        err
    );
}

#[test]
fn lenient_warns_and_keeps_value_when_intronic_indel_cannot_resolve() {
    // #682 regression: a substitution never reaches the genomic-sequence-dependent
    // intronic resolution pass, so warn-only mode has always echoed it with W4007.
    // An *indel* does reach that pass; with no genomic sequence to anchor it,
    // `normalize_tx` errored and — because warn-only mode (unlike reject mode) did
    // NOT short-circuit before normalization — that capability error propagated via
    // `?`, dropping the W4007 warning entirely. The spec-invalidity of an intronic
    // offset on a bare transcript is independent of provider capability, so warn-only
    // mode must still surface W4007 and echo the input, not hard-fail.
    //
    // `provider_intronic_nr` has the transcript but no genomic sequence, so the
    // intronic deletion cannot resolve.
    let normalizer = Normalizer::with_config(provider_intronic_nr(), NormalizeConfig::lenient());
    let variant = parse_hgvs("NR_INT.1:n.10+2_10+3del").expect("parse");

    // Warn-only must NOT propagate the capability error.
    let out = normalizer
        .normalize(&variant)
        .expect("warn-only must not fail when an intronic indel cannot resolve (#682)");
    // The spec-invalid form is echoed unchanged (not normalized/shifted).
    assert_eq!(out.to_string(), "NR_INT.1:n.10+2_10+3del");

    // And W4007 is surfaced rather than silently dropped.
    let diag = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("warn-only diagnostics");
    assert!(
        diag.warnings
            .iter()
            .any(|w| w.code() == "INTRONIC_ON_BARE_TRANSCRIPT"),
        "expected INTRONIC_ON_BARE_TRANSCRIPT warning, got {:?}",
        diag.warnings
    );
}

#[test]
fn silent_accepts_bare_nm_intronic_substitution_without_warning() {
    // Silent mode (`warn_accept()` → Accept): the intronic offset on a bare
    // NM_ is accepted without rejection and without emitting W4007.
    let normalizer = Normalizer::with_config(provider_intronic_nm(), NormalizeConfig::silent());
    let variant = parse_hgvs("NM_INTRON.1:c.30+1A>G").expect("parse");
    // Must not reject.
    let _ = normalizer
        .normalize(&variant)
        .expect("silent mode must accept a bare-NM intronic substitution without error");
    // Must not emit the INTRONIC_ON_BARE_TRANSCRIPT warning.
    let diag = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("silent diagnostics");
    assert!(
        !diag
            .warnings
            .iter()
            .any(|w| w.code() == "INTRONIC_ON_BARE_TRANSCRIPT"),
        "silent mode must not emit INTRONIC_ON_BARE_TRANSCRIPT, got {:?}",
        diag.warnings
    );
}

#[test]
fn silent_propagates_error_when_intronic_indel_cannot_resolve() {
    // Silent-mode counterpart to `lenient_warns_and_keeps_value_when_intronic_indel_cannot_resolve`.
    // Silent mode produces no EINTRONIC warning, so the `None => return Err(e)` arm in
    // `normalize_core` must still propagate the capability error for an intronic *indel* that
    // reaches the genomic-sequence-dependent resolution pass — no echo, no W4007 (#682).
    //
    // `provider_intronic_nr` has the transcript but no genomic sequence, so the intronic
    // deletion cannot resolve.
    let normalizer = Normalizer::with_config(provider_intronic_nr(), NormalizeConfig::silent());
    let variant = parse_hgvs("NR_INT.1:n.10+2_10+3del").expect("parse");

    // Silent mode has no EINTRONIC warning to recover with, so the capability error propagates.
    assert!(
        normalizer.normalize(&variant).is_err(),
        "silent mode must propagate the resolve-failure error for an intronic indel (#682)"
    );

    // And it must not emit the INTRONIC_ON_BARE_TRANSCRIPT warning (no echo path taken).
    let diag = normalizer.normalize_with_diagnostics(&variant);
    if let Ok(diag) = diag {
        assert!(
            !diag
                .warnings
                .iter()
                .any(|w| w.code() == "INTRONIC_ON_BARE_TRANSCRIPT"),
            "silent mode must not emit INTRONIC_ON_BARE_TRANSCRIPT, got {:?}",
            diag.warnings
        );
    }
}
