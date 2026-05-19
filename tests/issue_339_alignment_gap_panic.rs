//! Regression test for issue #339: `Normalizer::normalize` must not
//! panic when the reference window returned by the provider is shorter
//! than the HGVS interval span (e.g. because a cdot alignment includes
//! gaps that collapse the byte length). Pre-fix, ferro fires a
//! `debug_assert_eq!` at `src/normalize/mod.rs:2779`.
//!
//! Biocommons case (from #324 baseline-failures):
//!   `NG_032871.1:g.32476_53457delinsAATTAAGGTATA`
//!     bio:   rejects with HGVSInvalidVariantError
//!     ferro: panics with "ref_bytes length must match HGVS interval span"
//!
//! Reproducing the panic without the manifest requires a custom
//! `ReferenceProvider` that intentionally returns fewer bytes than the
//! requested half-open span. `TruncatingProvider` below does exactly
//! that — for any `get_sequence` request it returns a fixed-length
//! sequence regardless of the requested span. Once the fix is in
//! place, normalize() returns the input unchanged (no canonical
//! split could be applied without an authoritative ref window).

use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::{parse_hgvs, FerroError, NormalizeConfig, Normalizer, ReferenceProvider};

/// Provider that returns a short fixed sequence regardless of the
/// requested span — models the alignment-gap behavior of real cdot data
/// where the FASTA window is shorter than the HGVS-coordinate span.
struct TruncatingProvider;

impl ReferenceProvider for TruncatingProvider {
    fn get_transcript(&self, _id: &str) -> Result<Transcript, FerroError> {
        Err(FerroError::ReferenceNotFound {
            id: "<truncating-provider>".to_string(),
        })
    }

    fn get_sequence(&self, _id: &str, _start: u64, _end: u64) -> Result<String, FerroError> {
        // Return a fixed 12-byte slice regardless of the requested span.
        // A delins requesting a much larger window (e.g. 20kb) would have
        // its `ref_bytes.len()` here at 12 ≪ the HGVS span — exactly the
        // mismatch that fires the `debug_assert_eq!` panic pre-fix.
        Ok("AATTAAGGTATA".to_string())
    }

    fn get_genomic_sequence(
        &self,
        _contig: &str,
        _start: u64,
        _end: u64,
    ) -> Result<String, FerroError> {
        Ok("AATTAAGGTATA".to_string())
    }
}

#[test]
fn normalize_does_not_panic_when_ref_window_shorter_than_hgvs_span() {
    // Reproducer for the biocommons NG_032871.1 case. The variant span
    // (32476..53457 = 20,982 bp) is far larger than the 12-byte window
    // the provider returns. Pre-fix this triggers
    // `debug_assert_eq!(n, (hgvs_end - hgvs_start + 1) as usize)` at
    // src/normalize/mod.rs:2779. Post-fix the canonical-split helper
    // returns the input variant unchanged.
    let normalizer = Normalizer::with_config(TruncatingProvider, NormalizeConfig::default());
    let variant = parse_hgvs("NG_032871.1:g.32476_53457delinsAATTAAGGTATA").expect("parse");

    // Must not panic in debug builds. We don't pin the exact normalized
    // output — the test's purpose is to lock in graceful degradation.
    let _ = normalizer
        .normalize(&variant)
        .expect("normalize must not panic on alignment-gap-spanning delins");
}

#[test]
fn normalize_does_not_panic_on_smaller_alignment_gap_mismatch() {
    // Same shape, smaller span: 100_200delins (101 bp) with the
    // provider still returning only 12 bytes. The fix must handle any
    // mismatch, not just the extreme 20kb case.
    let normalizer = Normalizer::with_config(TruncatingProvider, NormalizeConfig::default());
    let variant = parse_hgvs("NC_000001.11:g.100_200delinsAATTAAGGTATA").expect("parse");
    let _ = normalizer
        .normalize(&variant)
        .expect("normalize must not panic on smaller ref-window mismatch");
}
