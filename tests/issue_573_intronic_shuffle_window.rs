//! Regression tests for issue #573: intronic shuffle window must be sized to
//! cover the enclosing intron.
//!
//! When the huge-intron cap in `intronic_window_bounds` fires (intron span
//! exceeds `MAX_INTRONIC_SHUFFLE_WINDOW = 64 KiB`), the fallback window is the
//! variant-only window, which may not include both intron edges. Before the fix
//! the caller silently used `saturating_sub` to compute `intron_rel_start`,
//! clamping to 1 when the near intron edge lay outside the fetched window.
//! That 1-valued boundary passed the `flip_intronic_for_strand` in-window guard
//! on minus-strand transcripts, producing incorrect normalization. After the
//! fix both callers (`normalize_intronic_cds` and `normalize_intronic_tx`)
//! explicitly check that both intron edges lie within the resolved window and
//! return a `ConversionError` when they do not.
//!
//! Fixture geometry (minus strand, `c.`/`n.` notation):
//!
//! ```text
//! genomic:  1001  1021                           69980  70000
//!                 |<——————— intron (~68960 bp) ————————>|
//!           [exon2: 1001..1020]                         [exon1: 69981..70000]
//!           tx 21-40 (3'-end, small genomic)            tx 1-20 (5'-end, large genomic)
//! ```
//!
//! The intron (genomic 1021..69980) is 68960 bases — above the 64 KiB cap.
//! Two variants are tested, each hitting one of the two out-of-window cases:
//!
//! - `c.20+3`: near the 5'-donor end (large genomic coords); after the cap
//!   fires, `intron_g_start` (1021) < `seq_start` (~69878). Before the fix
//!   this silently clamped `intron_rel_start` to 1.
//!
//! - `c.21-3`: near the 3'-acceptor end (small genomic coords); after the cap
//!   fires, `intron_g_end + 1` (69981) > `seq_end` (~1123). Before the fix
//!   this silently clamped `intron_rel_end` to `seq_len`.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, FerroError, MockProvider, Normalizer};

/// Build a mock provider with:
/// - A minus-strand transcript `NM_573TEST.1` whose single intron is
///   ~68960 bp (above the 64 KiB `MAX_INTRONIC_SHUFFLE_WINDOW` cap).
/// - A genomic sequence long enough to cover any variant-sized fetch window
///   used by the normalizer.
fn make_huge_intron_fixture() -> MockProvider {
    // Genomic sequence: 70 100 arbitrary bases (enough to cover up to position
    // 70099 inclusive; the largest fetch window in these tests ends at ~70078).
    let genomic_seq = "ACGT".repeat(17525); // 17525 * 4 = 70100 bytes
    assert_eq!(genomic_seq.len(), 70100);

    // Transcript sequence: 40 bases (two 20-bp exons).
    let tx_seq = "ACGT".repeat(10); // 40 bytes

    // Exons in transcript order.  On a minus-strand transcript the first exon
    // has the largest genomic coordinates.
    //   Exon 1: tx 1-20, genomic 69981-70000
    //   Exon 2: tx 21-40, genomic 1001-1020
    //
    // Computed intron (by `compute_introns`):
    //   intron.genomic_start = downstream.genomic_end + 1 = 1020 + 1 = 1021
    //   intron.genomic_end   = upstream.genomic_start - 1 = 69981 - 1 = 69980
    //   size = 69980 - 1021 + 1 = 68960  (> 65536 cap)
    let transcript = Transcript::new(
        "NM_573TEST.1".to_string(),
        Some("HUGEGENE".to_string()),
        Strand::Minus,
        tx_seq,
        Some(1),
        Some(40),
        vec![
            Exon::with_genomic(1, 1, 20, 69981, 70000),
            Exon::with_genomic(2, 21, 40, 1001, 1020),
        ],
        Some("chr_huge".to_string()),
        Some(1001),
        Some(70000),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );

    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr_huge", genomic_seq);
    provider.add_transcript(transcript);
    provider
}

/// Helper: normalize and return the `Err` variant, asserting it is a
/// `ConversionError`.
fn expect_conversion_error(provider: MockProvider, input: &str) -> FerroError {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse should succeed");
    let result = normalizer.normalize(&variant);
    match result {
        Err(e @ FerroError::ConversionError { .. }) => e,
        other => panic!("expected ConversionError for {input:?}, got {other:?}"),
    }
}

/// **Regression — near 5'-donor (`c.20+3`): `intron_g_start` outside the cap
/// fallback window.**
///
/// Variant: `c.20+3del` — 3 bases into the intron from the upstream exon.
/// On minus strand, this is at genomic position 69981 - 3 = 69978.
///
/// After the cap fires: `seq_start ≈ 69878`, `seq_end ≈ 70078`.
/// `intron_g_start = 1021 < seq_start` → guard fires → `ConversionError`.
///
/// Before the fix: `saturating_sub(1021, 69878)` = 0, `intron_rel_start` = 1,
/// the in-window guard in `flip_intronic_for_strand` silently passed,
/// and minus-strand normalization proceeded with a wrong `boundaries.left`.
#[test]
fn near_donor_huge_intron_cds_returns_error() {
    let provider = make_huge_intron_fixture();
    let err = expect_conversion_error(provider, "NM_573TEST.1:c.20+3del");
    assert!(
        format!("{err}").contains("intronic shuffle window"),
        "error should mention 'intronic shuffle window', got: {err}"
    );
}

/// **Regression — near 3'-acceptor (`c.21-3`): `intron_g_end + 1` outside the
/// cap fallback window.**
///
/// Variant: `c.21-3del` — 3 bases into the intron from the downstream exon.
/// On minus strand, this is at genomic position 1020 + 3 = 1023.
///
/// After the cap fires: `seq_start ≈ 923`, `seq_end ≈ 1123`.
/// `intron_g_end + 1 = 69981 > seq_end` → guard fires → `ConversionError`.
///
/// Before the fix: `69980.saturating_sub(923)` = 69057, which produces a
/// relative boundary far past the fetched window length, but the *left*
/// boundary (`intron_g_start = 1021`) would also be clamped via `saturating_sub`
/// only if `< seq_start` (here it isn't). The right boundary (69057) is >> seq_len,
/// so `flip_intronic_for_strand` would have caught it — but the explicit guard
/// makes the rejection happen earlier with a clearer diagnostic.
#[test]
fn near_acceptor_huge_intron_cds_returns_error() {
    let provider = make_huge_intron_fixture();
    let err = expect_conversion_error(provider, "NM_573TEST.1:c.21-3del");
    assert!(
        format!("{err}").contains("intronic shuffle window"),
        "error should mention 'intronic shuffle window', got: {err}"
    );
}

/// **`n.` parity — near 5'-donor via `normalize_intronic_tx`.**
///
/// Mirrors `near_donor_huge_intron_cds_returns_error` but uses `n.` notation,
/// which routes through `normalize_intronic_tx` instead of
/// `normalize_intronic_cds`. Both callers now carry the same guard.
#[test]
fn near_donor_huge_intron_tx_returns_error() {
    let provider = make_huge_intron_fixture();
    let err = expect_conversion_error(provider, "NM_573TEST.1:n.20+3del");
    assert!(
        format!("{err}").contains("intronic shuffle window"),
        "error should mention 'intronic shuffle window', got: {err}"
    );
}
