//! Reject `dupins<seq>` per HGVS DNA/duplication.md:92.
//!
//! The spec says:
//!
//! > "the variant is not described using `dupins`, a format not used in
//! > HGVS nomenclature."
//!
//! Marked `class="invalid"` in the spec — the strongest prohibition
//! register. Canonical alternatives:
//! - describe `del` and `ins` as separate edits,
//! - use `delins` for a contiguous deletion-insertion.
//!
//! Closes #445.

use ferro_hgvs::error::ErrorCode;
use ferro_hgvs::parse_hgvs;

fn assert_rejects(input: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "expected `dupins` form {input:?} to be rejected; got: {:?}",
        result.map(|v| v.to_string())
    );
    let err = result.unwrap_err();
    // The rejection must carry a structured `InvalidEdit` code so the
    // semantic refusal survives slash-form fallback paths rather than
    // being masked as an unstructured parse error.
    assert_eq!(
        err.code(),
        Some(ErrorCode::InvalidEdit),
        "dupins rejection must carry a structured InvalidEdit code for {input:?}",
    );
    let msg = err.to_string();
    assert!(
        msg.contains("dupins") && msg.contains("DNA/duplication.md"),
        "diagnostic must cite the spec; got: {msg}"
    );
}

fn assert_accepts(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("{input:?} must parse: {e}"));
    assert_eq!(parsed.to_string(), expected, "round-trip mismatch");
}

// ---- Rejected: `dupins<seq>` across coord systems ----

#[test]
fn rejects_genome_dupins() {
    assert_rejects("NC_000001.11:g.5207_5208dupinsATC");
}

#[test]
fn rejects_genome_single_position_dupins() {
    // The combination of `dupins` + single-position is doubly invalid;
    // reject for the dupins reason (#445), not for the single-position
    // reason (#446).
    let result = parse_hgvs("NC_000001.11:g.1000dupinsTAG");
    assert!(result.is_err());
    let msg = result.unwrap_err().to_string();
    assert!(
        msg.contains("dupins"),
        "diagnostic must mention dupins first; got: {msg}"
    );
}

#[test]
fn rejects_cds_dupins() {
    assert_rejects("NM_004006.2:c.20_21dupinsATG");
}

#[test]
fn rejects_rna_dupins() {
    assert_rejects("NM_004006.2:r.20_21dupinsacg");
}

#[test]
fn rejects_dupins_with_bracketed_payload() {
    assert_rejects("NC_000001.11:g.5207_5208dupins[ATC]");
}

// ---- Inside allele brackets ----

#[test]
fn rejects_dupins_inside_cis_allele() {
    assert_rejects("NM_004006.2:c.[100A>G;20_21dupinsATG]");
}

// ---- Canonical alternatives stay accepted ----

#[test]
fn accepts_canonical_delins_alternative() {
    assert_accepts(
        "NC_000001.11:g.5207_5208delinsATC",
        "NC_000001.11:g.5207_5208delinsATC",
    );
}

#[test]
fn accepts_canonical_dup_alone() {
    assert_accepts("NC_000001.11:g.5207_5208dup", "NC_000001.11:g.5207_5208dup");
}
