//! Issue #944: CDS<->tx is transcript-native across intra-exon CIGAR indels.
//!
//! Reference-backed (skips unless `FERRO_MANIFEST` points at a prepared
//! reference). Ground truth confirmed against the RefSeq mRNA and
//! VariantValidator: NM_015120.4 exon 0 CIGAR `M185 I3 M250`, so c./n.
//! numbering counts the 3 inserted bases and positions 3' of the insertion
//! must NOT shift by +3.

use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

use ferro_hgvs::convert::CoordinateMapper;
use ferro_hgvs::hgvs::location::{CdsPos, TxPos};
use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
use ferro_hgvs::reference::ReferenceProvider;

fn manifest_path() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let p = PathBuf::from(path);
        return if p.exists() { Some(p) } else { None };
    }
    let p = Path::new("benchmark-output/manifest.json");
    if p.exists() {
        return Some(p.to_path_buf());
    }
    None
}

fn provider() -> Option<Arc<MultiFastaProvider>> {
    static PROVIDER: OnceLock<Option<Arc<MultiFastaProvider>>> = OnceLock::new();
    PROVIDER
        .get_or_init(|| {
            let path = manifest_path()?;
            Some(Arc::new(
                MultiFastaProvider::from_manifest(&path)
                    .unwrap_or_else(|e| panic!("from_manifest({}) failed: {e}", path.display())),
            ))
        })
        .clone()
}

#[test]
fn nm_015120_cds_to_tx_transcript_native_across_insertion() {
    let Some(provider) = provider() else {
        eprintln!("issue_944: skipping — no manifest at FERRO_MANIFEST");
        return;
    };
    let tx = provider
        .get_transcript("NM_015120.4")
        .expect("NM_015120.4 present in reference");
    let mapper = CoordinateMapper::new(&tx);

    // c.1 -> n.112 (pre-insertion control); c.76 -> n.187 (unaffected).
    assert_eq!(mapper.cds_to_tx(&CdsPos::new(1)).unwrap().base, 112);
    assert_eq!(mapper.cds_to_tx(&CdsPos::new(76)).unwrap().base, 187);
    // c.77 -> n.188, c.87 -> n.198 (would be n.191/n.201 under the +3 bug).
    assert_eq!(mapper.cds_to_tx(&CdsPos::new(77)).unwrap().base, 188);
    assert_eq!(mapper.cds_to_tx(&CdsPos::new(87)).unwrap().base, 198);

    // Round-trip: tx 198 -> c.87.
    assert_eq!(mapper.tx_to_cds(&TxPos::new(198)).unwrap().base, 87);

    // Base-exactness: c.87 reference base is G (transcript 0-based index 197).
    let base = provider
        .get_sequence("NM_015120.4", 197, 198)
        .expect("read c.87 base");
    assert_eq!(base.to_ascii_uppercase(), "G");
}

#[test]
fn nm_000277_cds_to_tx_naive_across_deletion() {
    // PAH — carries a CIGAR deletion; the transcript axis must stay naive
    // (deletions add no transcript bases), i.e. contiguous with no shift.
    let Some(provider) = provider() else {
        eprintln!("issue_944: skipping — no manifest at FERRO_MANIFEST");
        return;
    };
    let tx = provider
        .get_transcript("NM_000277.3")
        .expect("NM_000277.3 present in reference");
    let mapper = CoordinateMapper::new(&tx);
    // Independent oracle: hardcoded expected transcript positions (not derived
    // from the method under test). c.1 anchors at n.115; a naive,
    // deletion-transparent axis places c.N at exactly n.(114 + N) with no
    // CIGAR-deletion shift. Under the pre-#944 bug, positions past the deletion
    // were pulled 5' by the deletion length and would miss these anchors.
    for (n, expected_tx) in [
        (1_i64, 115_i64),
        (50, 164),
        (100, 214),
        (300, 414),
        (500, 614),
    ] {
        assert_eq!(
            mapper.cds_to_tx(&CdsPos::new(n)).unwrap().base,
            expected_tx,
            "c.{n} must map to n.{expected_tx} (naive across the deletion)",
        );
        // Round-trip through the inverse re-derives c.N.
        assert_eq!(mapper.tx_to_cds(&TxPos::new(expected_tx)).unwrap().base, n);
    }
}
