//! #1472 — a zero-length `get_sequence` read is legal on a contig and rejected
//! on a transcript, so `apply_to_reference` declines every insertion and `dup`
//! on the transcript axes.
//!
//! ```text
//! genomic     get_sequence("NC_TEST.1", 9, 9)  ->  Ok("")
//! transcript  get_sequence("NR_TEST.1", 9, 9)  ->  Err("Position 9-9 out of range")
//! ```
//!
//! A zero-length window is what an **insertion** or a **duplication** has for
//! its reference span, so the transcript path declined both while accepting
//! substitutions and deletions on the same record.
//!
//! ## Why it matters
//!
//! `apply_to_reference` is documented as the encoding-**invariant** ground truth
//! — "two descriptions denote the same edit exactly when they produce the same
//! resulting bases". It was not, on the transcript axes: two spellings of one
//! variant got opposite answers, and the one it rejected is the one ferro's own
//! normalizer produces.
//!
//! ```text
//! NR_TEST.1:n.[9dup;10dup;11dup]   canonical_spdi = Ok("NR_TEST.1:9::TAA")
//! NR_TEST.1:n.9_10insTAA           canonical_spdi = Err("... window could not be read")
//! ```
//!
//! ## Four paths, not two
//!
//! The issue framed this as contig-versus-transcript. Surveying every
//! `get_sequence` implementation, the split is **2–2 and cuts across that
//! framing** — `FastaProvider` disagrees with *itself* depending on whether the
//! read is served from the index or from an mmap:
//!
//! | path | zero-length read at an in-range start |
//! |---|---|
//! | `MockProvider` genomic (`start > end` rejects) | `Ok("")` |
//! | `FastaProvider` indexed (`start >= actual_end` returns empty) | `Ok("")` |
//! | `Transcript::get_sequence` (`start < end` required) | **rejected** |
//! | `FastaProvider::get_sequence_from_mmap` (`start >= end` rejects) | **rejected** |
//!
//! So this fixes the two that reject, adopting the convention the other two
//! already implement rather than inventing a third.
//!
//! ## What stays an error
//!
//! Zero *length* and out of *range* are different questions, and only the first
//! is being relaxed. A zero-length read whose start does not name a base the
//! sequence has is still an error, matching the genomic path's `start >= len`
//! guard — including at exactly `len`, the junction past the last base.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::ReferenceProvider;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// 63 bases: 1-9 are `T`, 10 and 11 are `A`, 12 is `T`.
const CORE: &str = "TTTTTTTTTAATATATTTTAATATAATTAAAAAAATAATTTTTATAAATATATTATTTTAAAAA";

/// One provider carrying `CORE` twice: once as a contig, once as a non-coding
/// transcript. Serving identical bases under both is what makes the parity
/// assertions below a controlled comparison rather than two unrelated reads.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_TEST.1", CORE);
    let length = CORE.len() as u64;
    provider.add_transcript(Transcript::new(
        "NR_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        CORE.to_string(),
        None,
        None,
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// The reproducer, asserted as **agreement** between the two paths rather than
/// as two independent expectations: a change that made the contig path start
/// rejecting zero-length reads would satisfy a pair of separate `is_ok` checks
/// read the wrong way, but cannot satisfy this.
#[test]
fn a_zero_length_read_behaves_the_same_on_a_contig_and_a_transcript() {
    let provider = provider();
    let contig = provider.get_sequence("NC_TEST.1", 9, 9);
    let transcript = provider.get_sequence("NR_TEST.1", 9, 9);
    assert_eq!(
        contig.as_deref().ok(),
        Some(""),
        "a zero-length read on a contig is the empty string — this is the \
         established behaviour the transcript path is being aligned to"
    );
    assert_eq!(
        transcript.as_deref().ok(),
        Some(""),
        "and a transcript must answer identically: `[9, 9)` is an empty \
         interval, not an invalid one (#1472)"
    );
}

/// Non-empty reads must be untouched by the relaxation — the bases, not merely
/// the `Ok`-ness.
#[test]
fn ordinary_reads_are_unchanged() {
    let provider = provider();
    assert_eq!(provider.get_sequence("NR_TEST.1", 9, 11).unwrap(), "AA");
    assert_eq!(provider.get_sequence("NC_TEST.1", 9, 11).unwrap(), "AA");
    assert_eq!(
        provider.get_sequence("NR_TEST.1", 0, 9).unwrap(),
        "TTTTTTTTT"
    );
}

/// The discriminating half: zero *length* is now legal, out of *range* is not.
///
/// Without these a fix could simply drop the bounds check and pass the
/// reproducer. `999` is past the end; `63` is exactly the length, the junction
/// past the last base, which the contig path already refuses via `start >= len`
/// — so the transcript path must refuse it too or the two are still not aligned.
#[test]
fn an_out_of_range_zero_length_read_is_still_an_error() {
    let provider = provider();
    for (id, start) in [
        ("NR_TEST.1", 999u64),
        ("NC_TEST.1", 999),
        ("NR_TEST.1", CORE.len() as u64),
        ("NC_TEST.1", CORE.len() as u64),
    ] {
        assert!(
            provider.get_sequence(id, start, start).is_err(),
            "{id}: a zero-length read at {start} names no position the sequence \
             has, and must stay an error — only zero *length* is being relaxed"
        );
    }
}

/// An inverted range must stay rejected. `Transcript::get_sequence` guarded this
/// with the same `start < end` test that rejected the zero-length case, so
/// relaxing one must not relax the other.
#[test]
fn an_inverted_range_is_still_an_error() {
    let provider = provider();
    assert!(provider.get_sequence("NR_TEST.1", 11, 9).is_err());
    assert!(provider.get_sequence("NC_TEST.1", 11, 9).is_err());
}

/// The consequence the issue is really about: `apply_to_reference` must resolve
/// a pure insertion on a transcript axis, and must agree with the allele
/// spelling of the same variant.
///
/// Pinned as agreement between the two spellings, which is the property
/// `canonical_spdi` actually claims — asserting only that the insertion is `Ok`
/// would pass even if it resolved to the wrong bases.
#[test]
fn a_transcript_axis_insertion_denotes_the_same_variant_as_its_allele_spelling() {
    let normalizer = Normalizer::new(provider());
    let allele = parse_hgvs("NR_TEST.1:n.[9dup;10dup;11dup]").expect("fixture parses");
    let insertion = parse_hgvs("NR_TEST.1:n.9_10insTAA").expect("fixture parses");

    let allele_key = normalizer
        .canonical_spdi(&allele)
        .expect("the allele spelling already resolved before #1472");
    let insertion_key = normalizer
        .canonical_spdi(&insertion)
        .expect("the insertion spelling must resolve too — it is what ferro's own normalizer returns for the allele (#1472)");

    assert_eq!(
        insertion_key.to_string(),
        allele_key.to_string(),
        "`canonical_spdi` is the encoding-invariant key, so two spellings of one \
         variant must produce one triple"
    );
    // ...and that one triple is pinned, not merely shared. Agreement alone is
    // satisfied by two identically *wrong* keys, which is reachable here
    // because both spellings resolve through the same zero-length window this
    // fix relaxes: a window read at the wrong offset would move both sides
    // together and the equality above would still hold.
    assert_eq!(
        insertion_key.to_string(),
        "NR_TEST.1:9::TAA",
        "the deleted span is empty at interbase 9 and the inserted bases are \
         the three the allele duplicates"
    );
}

/// A `dup` on a transcript axis has the same zero-length reference window as an
/// insertion, and was declined for the same reason.
#[test]
fn a_transcript_axis_duplication_applies() {
    let normalizer = Normalizer::new(provider());
    let dup = parse_hgvs("NR_TEST.1:n.9dup").expect("fixture parses");
    let applied = normalizer
        .apply_to_reference(&dup)
        .expect("a lone `dup` on a transcript axis must apply (#1472)");
    assert_eq!(
        (applied.reference.as_str(), applied.resulting.as_str()),
        ("", "T"),
        "the reference window is the empty one this fix admits, and the result \
         is the single duplicated base — pinned as bases rather than as a \
         length relation, which `(\"AC\", \"ACG\")` would satisfy just as well"
    );
}

// ----------------------------------------------------------------------------
// The mmap path (`MmapFastaProvider`)
// ----------------------------------------------------------------------------

/// `MmapFastaProvider` is `#[cfg(feature = "mmap")]`, so these are gated to
/// match — as `issue_314_mmap_prefix_collision` is.
///
/// Gated as a **module** rather than with a file-level `#![cfg]` like that one,
/// because the rest of this file tests `Transcript::get_sequence` and the
/// `MockProvider` contig path, neither of which needs the feature. A file-level
/// gate would compile those away too, and the `soak` archive builds `--test it`
/// with **no** `--features` flag (`default = []`), so #1472's transcript-side
/// coverage would silently vanish from that job while it stayed green.
#[cfg(feature = "mmap")]
mod mmap_path {
    use super::CORE;
    use ferro_hgvs::reference::fasta::MmapFastaProvider;
    use ferro_hgvs::reference::ReferenceProvider;

    /// [`CORE`] as a single-contig FASTA, memory-mapped.
    ///
    /// The tests above reach only `Transcript::get_sequence` and the `MockProvider`
    /// contig path. `MmapFastaProvider::get_sequence_from_mmap` is the **second**
    /// of the two paths this fix relaxes — the module table above names it — and it
    /// needs a real file to map, so it cannot be reached through `MockProvider` at
    /// all. Written unwrapped so a coordinate is a byte offset into the one
    /// sequence line.
    fn mmap_provider() -> (tempfile::TempDir, MmapFastaProvider) {
        use std::io::Write;
        let dir = tempfile::tempdir().expect("create tempdir");
        let fa_path = dir.path().join("core.fa");
        let mut f = std::fs::File::create(&fa_path).expect("create fasta");
        writeln!(f, ">chrCORE").expect("write fasta header");
        writeln!(f, "{CORE}").expect("write fasta seq");
        f.sync_all().expect("sync fasta");
        drop(f);
        let provider = MmapFastaProvider::new(&fa_path).expect("mmap fasta");
        (dir, provider)
    }

    /// The mmap path serves a zero-length read as empty, like the other three.
    ///
    /// Without this the fix's `fasta.rs` half is untested: every other assertion in
    /// this file is served by `MockProvider` or by `Transcript::get_sequence`, so
    /// reverting `get_sequence_from_mmap` to `start >= end` leaves them all green.
    #[test]
    fn a_zero_length_read_is_empty_on_the_mmap_path_too() {
        let (_dir, provider) = mmap_provider();

        assert_eq!(
            provider.get_sequence("chrCORE", 9, 9).as_deref(),
            Ok(""),
            "a zero-length read at an in-range start reads no bases, and the \
             indexed read on the same provider already answers `\"\"` — the two \
             disagreeing meant an insertion resolved or not depending on whether \
             the record happened to be mmap-backed (#1472)"
        );
        assert_eq!(provider.get_sequence("chrCORE", 0, 0).as_deref(), Ok(""));

        // Ordinary reads are unaffected, bases included.
        assert_eq!(provider.get_sequence("chrCORE", 0, 3).as_deref(), Ok("TTT"));
        assert_eq!(
            provider.get_sequence("chrCORE", 9, 12).as_deref(),
            Ok("AAT"),
            "positions 10 and 11 are `A` and 12 is `T` on this core"
        );
    }

    /// The discriminating half on the mmap path: only zero *length* was relaxed.
    ///
    /// `start >= entry.length` still rejects, so the junction past the last base is
    /// an error here exactly as it is on the transcript path — the two conventions
    /// have to agree about where a sequence ends, not just about zero length.
    #[test]
    fn the_mmap_path_still_rejects_out_of_range_and_inverted_ranges() {
        let (_dir, provider) = mmap_provider();
        let len = CORE.len() as u64;

        assert!(
            provider.get_sequence("chrCORE", len, len).is_err(),
            "a zero-length read at the junction past the last base names no base \
             the contig has"
        );
        assert!(provider
            .get_sequence("chrCORE", len + 100, len + 100)
            .is_err());
        assert!(
            provider.get_sequence("chrCORE", 11, 9).is_err(),
            "an inverted range shared the old `start >= end` test, so relaxing it \
             had to keep rejecting this explicitly"
        );
        assert!(provider.get_sequence("chrCORE", 60, len + 1).is_err());
    }
}
