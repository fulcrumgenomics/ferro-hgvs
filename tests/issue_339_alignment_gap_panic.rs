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

/// Provider whose `get_sequence` returns a fixed-length sequence
/// regardless of the requested span. Models the alignment-gap behavior
/// of real cdot data where the FASTA window is shorter than the
/// HGVS-coordinate span. The fixed-length payload is configurable so
/// individual tests can drive both the "ref shorter than span" and
/// "ref longer than span" branches of the bail-out.
struct FixedLengthProvider {
    payload: &'static str,
}

impl ReferenceProvider for FixedLengthProvider {
    fn get_transcript(&self, _id: &str) -> Result<Transcript, FerroError> {
        Err(FerroError::ReferenceNotFound {
            id: "<fixed-length-provider>".to_string(),
        })
    }

    fn get_sequence(&self, _id: &str, _start: u64, _end: u64) -> Result<String, FerroError> {
        Ok(self.payload.to_string())
    }

    fn get_genomic_sequence(
        &self,
        _contig: &str,
        _start: u64,
        _end: u64,
    ) -> Result<String, FerroError> {
        Ok(self.payload.to_string())
    }
}

#[test]
fn normalize_returns_input_when_ref_window_shorter_than_hgvs_span() {
    // Reproducer for the biocommons NG_032871.1 case. The variant span
    // (32476..53457 = 20,982 bp) is far larger than the 12-byte window
    // the provider returns. Pre-fix this triggered
    // `debug_assert_eq!(n, (hgvs_end - hgvs_start + 1) as usize)` at
    // src/normalize/mod.rs:2779 (debug panic) and silently skipped the
    // split via decompose_delins' own bounds check in release. Post-fix
    // the canonical-split helper returns the input variant unchanged in
    // both build profiles.
    let provider = FixedLengthProvider {
        payload: "AATTAAGGTATA",
    };
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::default());
    let variant = parse_hgvs("NG_032871.1:g.32476_53457delinsAATTAAGGTATA").expect("parse");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize must not panic on alignment-gap-spanning delins");
    // Bail-out contract: input round-trips unchanged when the helper
    // cannot trust the ref window. We assert the canonical Display
    // round-trip rather than relying on a parsed-tree equality, since
    // upstream canonicalization may strip redundant fields without
    // shifting the variant.
    assert_eq!(
        format!("{}", normalized),
        "NG_032871.1:g.32476_53457delinsAATTAAGGTATA",
    );
}

#[test]
fn normalize_returns_input_when_ref_window_shorter_small_span() {
    // Same mismatch shape (ref short of span) at a smaller scale —
    // 101-bp span vs 12-byte provider response. Locks in that the fix
    // is independent of span magnitude.
    let provider = FixedLengthProvider {
        payload: "AATTAAGGTATA",
    };
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::default());
    let variant = parse_hgvs("NC_000001.11:g.100_200delinsAATTAAGGTATA").expect("parse");
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize must not panic on small-span ref-window mismatch");
    assert_eq!(
        format!("{}", normalized),
        "NC_000001.11:g.100_200delinsAATTAAGGTATA",
    );
}

#[test]
fn normalize_returns_input_when_ref_window_longer_than_hgvs_span() {
    // Inverse mismatch direction: provider returns MORE bytes than the
    // HGVS span requests. This shouldn't happen under the trait's
    // half-open `[start, end)` contract, but a misbehaving provider
    // would also have tripped the pre-fix debug_assert. The post-fix
    // `n != expected_span` predicate handles both directions; this
    // test locks that property in so a future refactor changing the
    // comparator to `n < expected_span` is caught.
    let provider = FixedLengthProvider {
        // 50 bytes, but the HGVS span below is only 11.
        payload: "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTACGTACGTAC",
    };
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::default());
    let variant = parse_hgvs("NC_000001.11:g.100_110delinsAAGCTT").expect("parse");
    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize must not panic on over-long ref window");
    assert_eq!(
        format!("{}", normalized),
        "NC_000001.11:g.100_110delinsAAGCTT",
    );
}
