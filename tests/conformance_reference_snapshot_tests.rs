//! Integrity tests for the committed mutalyzer-normalize transcript snapshot
//! (#719, increment I2).
//!
//! These run in CI with **no** prepared reference — they read only the two
//! committed snapshot files. They guard that the snapshot the builder
//! (`examples/build_conformance_snapshot.rs`) produced stays internally
//! consistent and that the headline `NM_003002.2` carries its OWN, version-exact
//! CDS frame (the #714 fix the whole snapshot arc exists to make reproducible).

use std::path::PathBuf;

use ferro_hgvs::conformance::reference_snapshot::{load_sequences, TranscriptSnapshot};
use sha2::{Digest, Sha256};

fn snapshot_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/mutalyzer-normalize/reference-snapshot")
}

fn sha256_hex(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    let mut s = String::with_capacity(digest.len() * 2);
    for byte in digest {
        s.push_str(&format!("{byte:02x}"));
    }
    s
}

/// The metadata and the FASTA describe exactly the same set of accessions, every
/// entry's recorded length and `sha256` match the actual bases, and bases are
/// uppercase ACGTN. This is the audit the committed `provenance.sha256` exists
/// for: if anyone edits the FASTA by hand without rebuilding, this fails.
#[test]
fn snapshot_is_internally_consistent() {
    let dir = snapshot_dir();
    let metadata = TranscriptSnapshot::from_json_path(dir.join("transcripts.metadata.json"))
        .expect("snapshot metadata loads");
    let sequences = load_sequences(&dir).expect("snapshot FASTA loads");

    assert!(
        !metadata.transcripts.is_empty(),
        "snapshot should carry at least one transcript"
    );
    assert_eq!(
        metadata.transcripts.len(),
        sequences.len(),
        "metadata and FASTA must describe the same number of accessions"
    );

    for (accession, entry) in &metadata.transcripts {
        let bases = sequences
            .get(accession)
            .unwrap_or_else(|| panic!("{accession}: in metadata but missing from FASTA"));
        assert_eq!(
            bases.len() as u64,
            entry.length,
            "{accession}: FASTA length disagrees with metadata length"
        );
        assert_eq!(
            sha256_hex(bases.as_bytes()),
            entry.provenance.sha256,
            "{accession}: provenance sha256 does not match the committed bases"
        );
        assert!(
            bases
                .bytes()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')),
            "{accession}: bases must be uppercase ACGTN"
        );
        // CDS is either fully absent (non-coding) or a sane in-bounds 1-based
        // range. A half-populated CDS (exactly one bound) is a malformed frame
        // the builder must never emit — assert the invariant, don't skip it.
        match (entry.cds_start, entry.cds_end) {
            (Some(start), Some(end)) => assert!(
                1 <= start && start <= end && end <= entry.length,
                "{accession}: CDS {start}..={end} out of 1..={} bounds",
                entry.length
            ),
            (None, None) => {}
            _ => panic!("{accession}: CDS must have both bounds or neither"),
        }
    }
}

/// The headline #714 case: `NM_003002.2` must carry its OWN frame
/// (`cds_start = 62`, 1-based), harvested cross-build from the GRCh37 cdot —
/// NOT `NM_003002.3`/`.4`'s frame, which the old fuzzy fallback would have
/// substituted. With `cds_start = 62`, `c.273del` rolls to `c.274del`.
#[test]
fn headline_nm_003002_2_has_own_cds_frame() {
    let dir = snapshot_dir();
    let metadata = TranscriptSnapshot::from_json_path(dir.join("transcripts.metadata.json"))
        .expect("snapshot metadata loads");

    let entry = metadata
        .transcripts
        .get("NM_003002.2")
        .expect("NM_003002.2 is in the snapshot");
    assert_eq!(entry.cds_start, Some(62), "NM_003002.2 CDS start (1-based)");
    assert_eq!(entry.cds_end, Some(541), "NM_003002.2 CDS end (1-based)");
    assert_eq!(entry.length, 1382);
    assert_eq!(entry.gene_symbol.as_deref(), Some("SDHD"));
    assert_eq!(entry.protein.as_deref(), Some("NP_002993.1"));
}
