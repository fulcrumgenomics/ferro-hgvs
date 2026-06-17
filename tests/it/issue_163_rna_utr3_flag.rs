//! Regression test for issue #163: `normalize_rna` drops the `*` UTR flag
//! on 3'UTR `r.*N` variants, and reads the wrong slice of the transcript
//! sequence because it treats `*N` as a transcript-1 position rather
//! than an offset from `cds_end`.
//!
//! Before the fix, normalize on `r.*1_*2del` over a transcript whose
//! 5'UTR begins with an `A` homopolymer would shuffle the deletion
//! against the 5'UTR slice and emit a non-UTR `r.M_Ndel`, silently
//! corrupting both the position and the region.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Build a single-exon transcript with:
/// - 5'UTR (tx 1-10):  `AAAAAAAAAC`
/// - CDS    (tx 11-50): `ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC`
/// - 3'UTR  (tx 51-60): `ACTTTTTTTT`
///
/// `r.*1` → tx 51 = `A`; `r.*2` → tx 52 = `C`. The neighbouring byte
/// `r.*3` → tx 53 = `T`, so a deletion at `r.*1_*2` cannot 3'-shift
/// (the rotation predicate `ref[del_start] == ref[del_end]` fails).
fn provider_with_utr_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence: String =
        "AAAAAAAAACATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCACTTTTTTTT".to_string();
    assert_eq!(sequence.len(), 60, "fixture length must be 60");
    let len = sequence.len() as u64;
    let exons = vec![Exon::new(1, 1, len)];
    let transcript = Transcript::new(
        "NM_TESTUTR.1".to_string(),
        Some("UTR_TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(11),
        Some(50),
        exons,
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

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(provider_with_utr_transcript());
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

#[test]
fn rna_utr3_del_preserves_star_flag_and_position() {
    // ref at *1='A', *2='C', *3='T'; rotation predicate fails so the
    // canonical 3'-most form is the input itself.
    assert_eq!(
        normalize_to_string("NM_TESTUTR.1:r.*1_*2del"),
        "NM_TESTUTR.1:r.*1_*2del",
    );
}

#[test]
fn rna_utr3_del_in_homopolymer_shifts_within_utr() {
    // 3'UTR T-tract at tx 53..60 (eight T's). Deleting r.*3 must shift
    // 3' to the last T (r.*10), staying inside the 3'UTR — both `*`
    // markers must be preserved.
    assert_eq!(
        normalize_to_string("NM_TESTUTR.1:r.*3del"),
        "NM_TESTUTR.1:r.*10del",
    );
}

#[test]
fn rna_utr3_delins_preserves_star_flag() {
    // delins(ac→ug) at *1_*2: rev-comp of "ac" is "gu" (RNA), so the
    // canonical form is `inv` (covered by the existing
    // delins-canonicalization rule). Either way the `*` flag must
    // survive normalization.
    let result = normalize_to_string("NM_TESTUTR.1:r.*1_*2delinsug");
    assert!(
        result.contains("*1") && result.contains("*2"),
        "expected `*1` and `*2` in normalized output, got {}",
        result,
    );
}

#[test]
fn rna_utr3_ins_preserves_star_flag() {
    // Insertion `gg` between r.*1 (A) and r.*2 (C). The inserted bases
    // match neither neighbor, so no ins→dup canonicalization fires and
    // the result stays as `ins`. The `*` markers must survive — covers
    // the insertion path through `normalize_rna`, which is a separate
    // canonicalization path from `del`/`delins`.
    assert_eq!(
        normalize_to_string("NM_TESTUTR.1:r.*1_*2insgg"),
        "NM_TESTUTR.1:r.*1_*2insgg",
    );
}

#[test]
fn rna_utr3_dup_preserves_star_flag() {
    // 3'UTR T-tract at tx 53..60 (eight T's). Duplicating the first T
    // (r.*3) must 3'-shift to the last T (r.*10). The `*` flag must
    // survive — covers the duplication path through `normalize_rna`,
    // which is a separate canonicalization path from `del`/`delins`/`ins`.
    assert_eq!(
        normalize_to_string("NM_TESTUTR.1:r.*3dup"),
        "NM_TESTUTR.1:r.*10dup",
    );
}

#[test]
fn rna_mixed_cds_utr3_del_shifts_into_utr() {
    // Mixed-interval: start in CDS, end in 3'UTR. `r.` is CDS-relative
    // (== `c.`, #469) with cds_start = 11, so the last CDS base (tx 50) is
    // `r.40` (= c.40), and `r.*1` = tx 51 = 'A'. Deletion of "CA" (tx 50_51)
    // can 3'-shift one position because tx 52 ('C') equals the first deleted
    // byte; the new range [tx 51, tx 52] = "AC" lies entirely in the 3'UTR,
    // so the canonical output must use `*1_*2` notation.
    assert_eq!(
        normalize_to_string("NM_TESTUTR.1:r.40_*1del"),
        "NM_TESTUTR.1:r.*1_*2del",
    );
}
