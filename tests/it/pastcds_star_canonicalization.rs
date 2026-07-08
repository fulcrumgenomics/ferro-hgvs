//! A plain past-CDS coordinate that maps into the 3'UTR must render in the
//! canonical `c.*N` form in lenient/silent mode (not echoed as an out-of-scheme
//! plain integer such as `c.600`), while the W4004 `PositionPastEnd` warning
//! still fires. Strict mode still rejects (see `issue_336_position_past_end`).
//! A position past the whole transcript has no valid `c.*N` form and is left
//! as-is.
//!
//! Motivated by `NG_012337.1(NM_003002.2):c.[274;600]`, where mutalyzer emits
//! `c.*120` (CDS length 480) and ferro previously echoed the non-canonical,
//! LOVD-invalid `c.600` (#920/#985).

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Single-exon synthetic transcript (mirrors `issue_336_position_past_end`):
///   CDS: tx 4..12 (`cds_len = 9`); 3'UTR: tx 13..20 (`c.*1..c.*8`).
fn provider_short_cds() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "AAAATGAAATAGCCCCCCCC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(4),
        Some(12),
        vec![Exon::new(1, 1, len)],
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

/// Normalize `input` under the default (lenient) config; return the rendered
/// variant string and whether a `POSITION_PAST_END` (W4004) warning fired.
fn norm(input: &str) -> (String, bool) {
    let provider = provider_short_cds();
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse");
    let out = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must not error on a past-end position");
    let warned = out.warnings.iter().any(|w| w.code() == "POSITION_PAST_END");
    (out.result.to_string(), warned)
}

#[test]
fn plain_past_cds_in_3utr_renders_star_and_still_warns() {
    // c.10 = *1 (cds_len 9; 3'UTR c.*1..c.*8). Canonicalize the plain past-CDS
    // coordinate to c.*1, but keep the W4004 flag (the input was out-of-scheme).
    let (s, warned) = norm("NM_TEST.1:c.10=");
    assert_eq!(s, "NM_TEST.1:c.*1=");
    assert!(warned, "the W4004 past-CDS warning must still fire");
}

#[test]
fn plain_past_cds_last_3utr_base_renders_star() {
    // c.17 = *8, the last base of the 3'UTR (still in-bounds).
    let (s, _) = norm("NM_TEST.1:c.17=");
    assert_eq!(s, "NM_TEST.1:c.*8=");
}

#[test]
fn position_past_whole_transcript_is_not_canonicalized() {
    // c.18 = *9, but the 3'UTR is only 8 long -> past the transcript end. There
    // is no valid c.*N for it, so leave it echoed (and still warned) rather than
    // inventing an out-of-range 3'UTR coordinate.
    let (s, warned) = norm("NM_TEST.1:c.18=");
    assert_eq!(s, "NM_TEST.1:c.18=");
    assert!(warned);
}

#[test]
fn plain_past_cds_del_in_3utr_homopolymer_is_shuffled_and_idempotent() {
    // The 3'UTR is CCCCCCCC (c.*1..c.*8) — a homopolymer, so a deletion 3'-shifts
    // to its 3'-most position (c.*8del). A plain past-CDS `c.10del` (= *1del) must
    // therefore normalize straight to the fully-shuffled canonical `c.*8del` on
    // the FIRST pass — not the merely position-rewritten `c.*1del` — otherwise a
    // second normalize pass would move it again, breaking idempotency.
    let (s1, warned) = norm("NM_TEST.1:c.10del");
    assert_eq!(
        s1, "NM_TEST.1:c.*8del",
        "past-CDS del must render fully 3'-shifted"
    );
    assert!(warned, "the W4004 past-CDS warning must still fire");
    // Idempotent: normalize(normalize(x)) == normalize(x).
    let (s2, _) = norm(&s1);
    assert_eq!(s2, s1, "second normalize pass must not move it again");
}
