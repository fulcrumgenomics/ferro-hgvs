//! Issue #392 — extend W4004 PositionPastEnd to intronic offsets.
//!
//! PR #342 wired W4004 for plain `c.<N>` (cds-end), `c.*<N>`
//! (transcript-end), and `c.-<N>` (5'utr-start), and PR #347 added the
//! `n.<N>` (transcript-end) case. All four short-circuit when the
//! position is intronic (offset present) because the bound depends on
//! intron size, which requires genomic alignment data. This file pins
//! the extension: when intron length is computable from the exons'
//! genomic coords, `|offset| > intron_length` must fire `W4004`.
//!
//! Spec basis (`assets/hgvs-nomenclature/docs/background/numbering.md`):
//! - "nucleotides at the 5' end of an intron are numbered relative to
//!   the last nucleotide of the directly upstream exon, followed by a
//!   `+`"; "in the middle of the intron, nucleotide numbering changes
//!   from `+` to `-`". So `+M` and `-M` are both bounded by the intron
//!   length (any `|M|` exceeding the intron length is past-end).
//! - "a coding DNA reference sequence does not contain intron or 5'
//!   and 3' gene flanking sequences, and can therefore not be used as
//!   a reference to describe variants in these regions"; "Correct
//!   descriptions refer to a genomic reference sequence" — so
//!   intronic bound checks are conditional on the genomic alignment
//!   being available (`Exon::genomic_start`/`genomic_end` populated).
//!   When unavailable, ferro conservatively skips per the existing
//!   W4004 missing-metadata convention.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

/// 2-exon plus-strand transcript with explicit genomic coords on each
/// exon, so the intron between them has a known 20-bp length.
///
/// Transcript layout:
/// ```text
///   tx        1         2
///         1234567890 ⇢⇢ 1234567890 12
///   seq:  AAAATGCCCC    GGGGTAGAAT AA
///   genomic
///   exon 1:    100-109                  (10 bases)
///   intron 1:           110-129          (20 bases — the test bound)
///   exon 2:                       130-141 (12 bases)
/// ```
///
/// Because cds_start=1 and cds_end=22, `c.<N>` maps 1:1 to tx `<N>`.
/// The exon-1/exon-2 junction is between c.10 and c.11; intron 1 is
/// the 20-bp gap between them.
///
/// Intron-1 length = 129 - 110 + 1 = 20 bp. So:
/// - `c.10+M` valid for `M` in `1..=20`
/// - `c.10+21` past intron-end → W4004 fires
/// - `c.11-M` valid for `M` in `1..=20`
/// - `c.11-21` past intron-end → W4004 fires
fn provider_with_intronic_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let transcript = Transcript::new(
        "NM_INT.1".to_string(),
        Some("INT".to_string()),
        Strand::Plus,
        Some("AAAATGCCCCGGGGTAGAATAA".to_string()),
        Some(1),
        Some(22),
        vec![
            Exon::with_genomic(1, 1, 10, 100, 109),
            Exon::with_genomic(2, 11, 22, 130, 141),
        ],
        Some("chr_int".to_string()),
        Some(100),
        Some(141),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

/// Same shape as `provider_with_intronic_transcript` but with
/// `Exon::new` (no genomic coords) so the intron-length computation
/// returns None and the check must conservatively skip.
fn provider_with_intronic_transcript_no_genomic_coords() -> MockProvider {
    let mut provider = MockProvider::new();
    let transcript = Transcript::new(
        "NM_NOGEN.1".to_string(),
        Some("NOGEN".to_string()),
        Strand::Plus,
        Some("AAAATGCCCCGGGGTAGAATAA".to_string()),
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
    );
    provider.add_transcript(transcript);
    provider
}

/// 2-exon plus-strand non-coding transcript fixture (no CDS) — same
/// 20-bp intron-1 shape as `provider_with_intronic_transcript`, used
/// for the `n.` axis cases.
fn provider_with_intronic_noncoding_transcript() -> MockProvider {
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
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

/// 3-exon plus-strand transcript with introns inside both the 5'UTR
/// and the 3'UTR, used to exercise `cds_pos_to_tx_boundary`'s
/// `pos.base < 0` (5'UTR-intronic) and `pos.utr3` (3'UTR-intronic)
/// branches that the `provider_with_intronic_transcript` fixture (CDS
/// = entire transcript) does not reach.
///
/// Transcript layout (tx coords on a 33-bp transcript):
/// ```text
///   tx pos: 1..6  | 7..18    | 19..33
///   region: 5'UTR | CDS      | 3'UTR
///   exons : exon1 | exon2    | exon3
///   genomic: 100..105 (6 bp) | 130..141 (12 bp) | 170..184 (15 bp)
///   intron 1 (between exons 1 and 2): genomic 106..129 = 24 bp
///   intron 2 (between exons 2 and 3): genomic 142..169 = 28 bp
/// ```
///
/// cds_start = 7 (first base of exon 2 = c.1), cds_end = 18 (last
/// base of exon 2 = c.12). Intron 1 sits entirely in the 5'UTR
/// (between tx 6 / c.-1 and tx 7 / c.1); intron 2 sits entirely in
/// the 3'UTR (between tx 18 / c.12 and tx 19 / c.*1). Both UTR
/// introns therefore live at the *exon-end / next-exon-start*
/// boundary, addressed via `c.-1+M` (UTR-5 intron 5' side),
/// `c.1-M` (UTR-5 intron 3' side), `c.12+M` (UTR-3 intron 5' side),
/// and `c.*1-M` (UTR-3 intron 3' side).
fn provider_with_utr_intronic_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let transcript = Transcript::new(
        "NM_UTRINT.1".to_string(),
        Some("UTRINT".to_string()),
        Strand::Plus,
        Some("CCCCCCATGAAATAACCCCCCCCCCCCCCCAAA".to_string()),
        Some(7),
        Some(18),
        vec![
            Exon::with_genomic(1, 1, 6, 100, 105),
            Exon::with_genomic(2, 7, 18, 130, 141),
            Exon::with_genomic(3, 19, 33, 170, 184),
        ],
        Some("chr_int".to_string()),
        Some(100),
        Some(184),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

/// 2-exon **minus-strand** transcript with explicit genomic coords —
/// mirror of `provider_with_intronic_transcript` to pin
/// strand-invariance of the intron-length math.
///
/// Per `compute_introns` (`src/reference/transcript.rs:649-655`),
/// minus-strand intron genomic coords are derived from
/// `downstream.genomic_end + 1 .. upstream.genomic_start - 1`. So with:
/// ```text
///   exon 1 (tx 1-10, coding-start side):  genomic 200..209  (10 bp)
///   exon 2 (tx 11-22, coding-end side):   genomic 130..141  (12 bp)
///   intron 1: genomic 142..199 = 58 bp
/// ```
/// (exon 1 maps to the higher genomic range because the transcript
/// reads 3'→5' along the genome.)
fn provider_with_minus_strand_intronic_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let transcript = Transcript::new(
        "NM_MINT.1".to_string(),
        Some("MINT".to_string()),
        Strand::Minus,
        Some("AAAATGCCCCGGGGTAGAATAA".to_string()),
        Some(1),
        Some(22),
        vec![
            Exon::with_genomic(1, 1, 10, 200, 209),
            Exon::with_genomic(2, 11, 22, 130, 141),
        ],
        Some("chr_mint".to_string()),
        Some(130),
        Some(209),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

/// Returns `true` iff the lenient-mode normalize call emits a
/// `POSITION_PAST_END` (W4004) warning for `input`.
///
/// When W4004 fires the bounds gate at `normalize_cds` / `normalize_tx`
/// short-circuits with `Ok(canonical, vec![warning])`, so this helper
/// reports `true`. When it does not fire the normalize call may
/// proceed to the downstream intronic-projection pass (which requires
/// genomic sequence data and surfaces `FerroError::IntronicVariant` on
/// fixtures that omit it) — that `Err` is treated as "W4004 did not
/// fire" because if it had, the bounds gate would have intercepted
/// before the intronic pass ran.
fn position_past_end_fires(provider: MockProvider, input: &str) -> bool {
    let config = NormalizeConfig::lenient();
    let normalizer = Normalizer::with_config(provider, config);
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for {:?}: {:?}", input, e));
    match normalizer.normalize_with_diagnostics(&variant) {
        Ok(result) => result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        // The bounds gate runs BEFORE the downstream intronic pass, so any
        // downstream `IntronicVariant` (or similar) error implies W4004 did
        // not fire — if it had, the gate would have returned Ok with the
        // warning intercepted at the top of normalize_cds / normalize_tx.
        Err(_) => false,
    }
}

// --- c. axis: + offset past intron end --------------------------------------

#[test]
fn c_dot_plus_offset_past_intron_end_emits_w4004() {
    // `c.10+21del`: intron-1 is 20 bp, so offset 21 is past-end.
    assert!(
        position_past_end_fires(provider_with_intronic_transcript(), "NM_INT.1:c.10+21del"),
        "expected POSITION_PAST_END (W4004) for c.10+21 on a 20-bp intron",
    );
}

#[test]
fn c_dot_plus_offset_within_intron_does_not_emit_w4004() {
    // `c.10+5del`: well within the 20-bp intron — the bounds gate must
    // not fire. Downstream normalize may still error on missing genomic
    // sequence data (intronic projection requires it); that's fine and
    // distinct from a W4004 false-positive (the helper interprets the
    // downstream error as "no W4004 fired").
    assert!(
        !position_past_end_fires(provider_with_intronic_transcript(), "NM_INT.1:c.10+5del"),
        "unexpected POSITION_PAST_END for c.10+5 (within 20-bp intron)",
    );
}

#[test]
fn c_dot_plus_offset_at_intron_end_boundary_does_not_emit_w4004() {
    // `c.10+20del`: offset exactly equal to intron length — the last
    // base of the intron, addressed via the upstream-exon "+" form.
    // Within bound, must NOT fire.
    assert!(
        !position_past_end_fires(provider_with_intronic_transcript(), "NM_INT.1:c.10+20del"),
        "unexpected POSITION_PAST_END for c.10+20 (at intron-end boundary)",
    );
}

// --- c. axis: - offset past intron start ------------------------------------

#[test]
fn c_dot_minus_offset_past_intron_end_emits_w4004() {
    // `c.11-21del`: 20-bp intron, offset 21 → past-end via the
    // downstream-exon "-" form.
    assert!(
        position_past_end_fires(provider_with_intronic_transcript(), "NM_INT.1:c.11-21del"),
        "expected POSITION_PAST_END (W4004) for c.11-21 on a 20-bp intron",
    );
}

#[test]
fn c_dot_minus_offset_within_intron_does_not_emit_w4004() {
    // `c.11-5del`: well within the 20-bp intron.
    assert!(
        !position_past_end_fires(provider_with_intronic_transcript(), "NM_INT.1:c.11-5del"),
        "unexpected POSITION_PAST_END for c.11-5 (within 20-bp intron)",
    );
}

// --- n. axis: intronic offset symmetry --------------------------------------

#[test]
fn n_dot_plus_offset_past_intron_end_emits_w4004() {
    // `n.10+21del`: n. axis intronic offset on a non-coding transcript.
    // Same 20-bp intron as the c. fixture.
    assert!(
        position_past_end_fires(
            provider_with_intronic_noncoding_transcript(),
            "NR_INT.1:n.10+21del"
        ),
        "expected POSITION_PAST_END (W4004) for n.10+21 on a 20-bp intron",
    );
}

#[test]
fn n_dot_minus_offset_past_intron_end_emits_w4004() {
    assert!(
        position_past_end_fires(
            provider_with_intronic_noncoding_transcript(),
            "NR_INT.1:n.11-21del"
        ),
        "expected POSITION_PAST_END (W4004) for n.11-21 on a 20-bp intron",
    );
}

#[test]
fn n_dot_plus_offset_within_intron_does_not_emit_w4004() {
    assert!(
        !position_past_end_fires(
            provider_with_intronic_noncoding_transcript(),
            "NR_INT.1:n.10+5del"
        ),
        "unexpected POSITION_PAST_END for n.10+5 (within 20-bp intron)",
    );
}

// --- Conservative skip when genomic coords missing --------------------------

#[test]
fn intronic_offset_check_skips_when_genomic_coords_missing() {
    // Same shape as the c.10+21 test but the transcript's exons have no
    // genomic coords. Intron length is not computable → conservative
    // skip per the existing W4004 missing-metadata convention.
    assert!(
        !position_past_end_fires(
            provider_with_intronic_transcript_no_genomic_coords(),
            "NM_NOGEN.1:c.10+21del"
        ),
        "expected conservative skip when genomic coords are missing",
    );
}

// --- Strict mode promotes the intronic past-end warning to an error --------

#[test]
fn strict_mode_promotes_intronic_past_end_to_error() {
    // On the spec-valid genomic-context form, a past-intron-end intronic
    // offset still promotes the W4004 PositionPastEnd warning to
    // FerroError::InvalidCoordinates in strict mode. (The bare form now
    // rejects earlier as EINTRONIC — see the companion test below.)
    let provider = provider_with_intronic_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NG_012337.1(NM_INT.1):c.10+21del").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("strict mode must reject c.10+21del as past-intron-end");
    match err {
        FerroError::InvalidCoordinates { msg } => {
            assert!(
                msg.contains("W4004") || msg.contains("intron"),
                "strict-mode error message must reference W4004 or intron; got: {msg}",
            );
        }
        other => panic!("unexpected error variant: {other:?}"),
    }
}

#[test]
fn strict_mode_rejects_bare_intronic_as_eintronic() {
    // The bare-transcript form is spec-invalid (#486 EINTRONIC) and rejects
    // before the past-intron-end check runs.
    let provider = provider_with_intronic_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_INT.1:c.10+21del").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("strict mode must reject bare-NM intronic as EINTRONIC");
    assert!(
        matches!(err, FerroError::IntronicVariant { .. }),
        "expected IntronicVariant (EINTRONIC), got {err:?}",
    );
}

// --- No-intron edge case: offset off the last exon -------------------------

#[test]
fn intronic_offset_past_last_exon_skips_no_intron_exists() {
    // `c.22+5`: c.22 is the last base of the last exon, so there is no
    // intron 3' of it. `find_intron_at_tx_boundary` returns None →
    // conservative skip (no false positive). This pins the "no intron
    // exists" branch of the new logic.
    assert!(
        !position_past_end_fires(provider_with_intronic_transcript(), "NM_INT.1:c.22+5del"),
        "unexpected POSITION_PAST_END for c.22+5 when no intron exists 3' of the last exon",
    );
}

// --- 5'UTR-intronic (c.-N+M / c.-N-M) --------------------------------------

#[test]
fn c_dot_5utr_intronic_plus_offset_past_intron_end_emits_w4004() {
    // Intron 1 of `NM_UTRINT.1` is the 24-bp gap between exon 1 (tx 6 =
    // c.-1, last 5'UTR base) and exon 2 (tx 7 = c.1, first CDS base).
    // `c.-1+25` (25 > 24) is past-end via the 5'UTR-intron's "+" form.
    assert!(
        position_past_end_fires(
            provider_with_utr_intronic_transcript(),
            "NM_UTRINT.1:c.-1+25del"
        ),
        "expected POSITION_PAST_END (W4004) for c.-1+25 on a 24-bp 5'UTR intron",
    );
}

#[test]
fn c_dot_5utr_intronic_plus_offset_within_intron_does_not_emit_w4004() {
    // `c.-1+5del` is well within the 24-bp 5'UTR intron.
    assert!(
        !position_past_end_fires(
            provider_with_utr_intronic_transcript(),
            "NM_UTRINT.1:c.-1+5del"
        ),
        "unexpected POSITION_PAST_END for c.-1+5 (within 24-bp 5'UTR intron)",
    );
}

#[test]
fn c_dot_cds_start_minus_offset_past_intron_5utr_emits_w4004() {
    // `c.1-25del` — `c.1` is the first CDS base; the intron immediately
    // 5' of it is the same 24-bp 5'UTR intron. Past-end via the "-" form.
    assert!(
        position_past_end_fires(
            provider_with_utr_intronic_transcript(),
            "NM_UTRINT.1:c.1-25del"
        ),
        "expected POSITION_PAST_END (W4004) for c.1-25 on a 24-bp 5'UTR intron",
    );
}

// --- 3'UTR-intronic (c.*N+M / c.NNN+M off last CDS base) -------------------

#[test]
fn c_dot_3utr_intronic_plus_offset_past_intron_end_emits_w4004() {
    // Intron 2 of `NM_UTRINT.1` is the 28-bp gap between exon 2 (tx 18
    // = c.12, last CDS base) and exon 3 (tx 19 = c.*1, first 3'UTR
    // base). `c.12+29` (29 > 28) is past-end on the 3'UTR-intron's "+"
    // form.
    assert!(
        position_past_end_fires(
            provider_with_utr_intronic_transcript(),
            "NM_UTRINT.1:c.12+29del"
        ),
        "expected POSITION_PAST_END (W4004) for c.12+29 on a 28-bp 3'UTR intron",
    );
}

#[test]
fn c_dot_3utr_intronic_minus_offset_past_intron_end_emits_w4004() {
    // `c.*1-29del` — `c.*1` is the first 3'UTR base; the intron 5' of
    // it is the same 28-bp 3'UTR intron. Past-end via the "-" form.
    assert!(
        position_past_end_fires(
            provider_with_utr_intronic_transcript(),
            "NM_UTRINT.1:c.*1-29del"
        ),
        "expected POSITION_PAST_END (W4004) for c.*1-29 on a 28-bp 3'UTR intron",
    );
}

#[test]
fn c_dot_3utr_intronic_within_intron_does_not_emit_w4004() {
    // `c.*1-5del` is well within the 28-bp 3'UTR intron.
    assert!(
        !position_past_end_fires(
            provider_with_utr_intronic_transcript(),
            "NM_UTRINT.1:c.*1-5del"
        ),
        "unexpected POSITION_PAST_END for c.*1-5 (within 28-bp 3'UTR intron)",
    );
}

// --- Minus-strand symmetry --------------------------------------------------

#[test]
fn minus_strand_c_dot_plus_offset_past_intron_end_emits_w4004() {
    // Minus-strand fixture: intron 1 spans genomic 142..199 = 58 bp.
    // `c.10+59` (59 > 58) is past-end. Pins that the intron-length
    // computation works through the minus-strand `compute_introns`
    // branch (`downstream.genomic_end + 1 .. upstream.genomic_start - 1`).
    assert!(
        position_past_end_fires(
            provider_with_minus_strand_intronic_transcript(),
            "NM_MINT.1:c.10+59del"
        ),
        "expected POSITION_PAST_END (W4004) for c.10+59 on a minus-strand 58-bp intron",
    );
}

#[test]
fn minus_strand_c_dot_plus_offset_within_intron_does_not_emit_w4004() {
    // `c.10+5del` is well within the 58-bp minus-strand intron.
    assert!(
        !position_past_end_fires(
            provider_with_minus_strand_intronic_transcript(),
            "NM_MINT.1:c.10+5del"
        ),
        "unexpected POSITION_PAST_END for minus-strand c.10+5 (within 58-bp intron)",
    );
}
