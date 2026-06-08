//! Reject a coding/genomic/mito coordinate system on a non-coding RNA
//! (`NR_`/`XR_`) reference (#486, `ECOORDINATESYSTEMMISMATCH`).
//!
//! An `NR_`/`XR_` accession is a non-coding RNA transcript: it has no CDS (so
//! `c.` is meaningless), it is not a genomic reference (`g.`), and it is not
//! the mitochondrion (`m.`). The only valid coordinate systems are `n.`
//! (non-coding transcript) and `r.` (RNA).

use ferro_hgvs::error::ErrorCode;
use ferro_hgvs::parse_hgvs;

fn assert_rejects(input: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "expected {input:?} to be rejected; got: {:?}",
        result.map(|v| v.to_string())
    );
    assert_eq!(
        result.unwrap_err().code(),
        Some(ErrorCode::CoordinateSystemMismatch),
        "rejection must carry CoordinateSystemMismatch for {input:?}",
    );
}

fn assert_accepts(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("{input:?} must parse: {e}"));
    assert_eq!(parsed.to_string(), expected, "round-trip mismatch");
}

// ---- Rejected: c./g./m. on a non-coding RNA reference ----

#[test]
fn rejects_coding_on_noncoding() {
    assert_rejects("NR_038420.1:c.10del");
}

#[test]
fn rejects_genomic_on_noncoding() {
    assert_rejects("NR_038420.1:g.10del");
}

#[test]
fn rejects_mito_on_noncoding() {
    assert_rejects("NR_038420.1:m.10del");
}

#[test]
fn rejects_coding_on_predicted_noncoding_xr() {
    assert_rejects("XR_001234.1:c.10del");
}

#[test]
fn rejects_circular_on_noncoding() {
    assert_rejects("NR_038420.1:o.10del");
}

#[test]
fn rejects_coord_mismatch_inside_allele() {
    assert_rejects("NR_038420.1:c.[10del;20del]");
}

// ---- Accepted: n./r. on NR_; and c./g./m. on their proper references ----

#[test]
fn accepts_noncoding_n_on_nr() {
    assert_accepts("NR_038420.1:n.10del", "NR_038420.1:n.10del");
}

#[test]
fn accepts_noprefix_inference_on_nr_as_noncoding() {
    // A bare (no coordinate prefix) description on a non-coding RNA accession
    // must infer `n.` (the only valid system), not `c.` — otherwise the new
    // coordinate-system validator would hard-fail legitimate inferred input.
    assert_accepts("NR_038420.1:10del", "NR_038420.1:n.10del");
}

#[test]
fn accepts_noprefix_inference_on_xr_as_noncoding() {
    assert_accepts("XR_001234.1:10del", "XR_001234.1:n.10del");
}

#[test]
fn accepts_noprefix_inference_on_nm_stays_coding() {
    // Coding transcripts (`NM_`/`XM_`) keep inferring `c.` for bare input.
    assert_accepts("NM_000088.3:10del", "NM_000088.3:c.10del");
}

#[test]
fn accepts_rna_r_on_nr() {
    assert_accepts("NR_038420.1:r.10del", "NR_038420.1:r.10del");
}

#[test]
fn accepts_coding_on_nm() {
    assert_accepts("NM_000088.3:c.10del", "NM_000088.3:c.10del");
}

#[test]
fn accepts_noncoding_form_on_coding_transcript() {
    // A coding transcript may still be described in n. — must NOT be rejected
    // (the rule keys on the NR_/XR_ accession, not the coordinate system).
    assert_accepts("NM_000088.3:n.10del", "NM_000088.3:n.10del");
}

#[test]
fn accepts_genomic_on_nc() {
    assert_accepts("NC_000001.11:g.10del", "NC_000001.11:g.10del");
}

#[test]
fn accepts_mito_on_mito_accession() {
    assert_accepts("NC_012920.1:m.10del", "NC_012920.1:m.10del");
}
