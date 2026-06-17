//! Reject single-position insertions per HGVS DNA/insertion.md:95-101.
//!
//! Q&A in the spec:
//!
//! > "Can I describe a variant as `g.123insG`?
//! > **No**, since the description is not unequivocal, it is not allowed.
//! > What does the description mean, the insertion of a `G` **at**
//! > position `g.123` or the insertion of a `G` **after** position
//! > `g.123`?"
//!
//! The canonical form requires a 2-position adjacent range:
//! `g.123_124insG`. Same rule applies to c., n., r., m., o. coordinates
//! by construct-symmetry.
//!
//! Closes #446.

use ferro_hgvs::error::ErrorCode;
use ferro_hgvs::parse_hgvs;

fn assert_rejects(input: &str) {
    let result = parse_hgvs(input);
    assert!(
        result.is_err(),
        "expected single-position insertion {input:?} to be rejected; got: {:?}",
        result.map(|v| v.to_string())
    );
    let err = result.unwrap_err();
    // The rejection must carry a structured `InvalidEdit` code (mirroring
    // `validate_no_dupins`, #445) so the semantic refusal survives
    // slash-form fallback paths rather than being masked as an
    // unstructured parse error.
    assert_eq!(
        err.code(),
        Some(ErrorCode::InvalidEdit),
        "single-position insertion rejection must carry a structured InvalidEdit code for {input:?}",
    );
    let msg = err.to_string();
    assert!(
        msg.contains("single-position insertion") && msg.contains("DNA/insertion.md"),
        "error message must cite the spec and explain the problem; got: {msg}"
    );
}

fn assert_accepts(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("{input:?} must parse: {e}"));
    assert_eq!(parsed.to_string(), expected, "round-trip mismatch");
}

// ---- Rejected: single-position insertion across every coord system ----

#[test]
fn rejects_genome_single_position_insertion() {
    assert_rejects("NC_000023.11:g.123insG");
}

#[test]
fn rejects_cds_single_position_insertion() {
    assert_rejects("NM_004006.2:c.456insT");
}

#[test]
fn rejects_cds_intronic_single_position_insertion() {
    // c.8010+30 is a single intronic position; with `ins<seq>` it's
    // still a single-position anchor and rejected (matches the
    // LRG_135t1:c.8010+30insALU case from ClinVar fixtures).
    assert_rejects("NM_000051.3:c.8010+30insALU");
}

#[test]
fn rejects_noncoding_single_position_insertion() {
    assert_rejects("NR_004430.2:n.789insC");
}

#[test]
fn rejects_rna_single_position_insertion() {
    assert_rejects("NM_004006.2:r.456insu");
}

#[test]
fn rejects_mt_single_position_insertion() {
    assert_rejects("NC_012920.1:m.3243insA");
}

#[test]
fn rejects_organelle_single_position_insertion() {
    // o. (organelle/circular) is named in the module-doc rule list;
    // pin that the single-position rejection actually covers it.
    assert_rejects("NC_012920.1:o.3243insA");
}

// ---- Accepted: spec-canonical 2-position range form ----
// Coverage spans every coord system the rejection tests above touch
// (g./c./n./r./m./o.) so the accept/reject pairing is symmetric.

#[test]
fn accepts_genome_two_position_insertion() {
    assert_accepts("NC_000023.11:g.123_124insG", "NC_000023.11:g.123_124insG");
}

#[test]
fn accepts_cds_two_position_insertion() {
    assert_accepts("NM_004006.2:c.456_457insT", "NM_004006.2:c.456_457insT");
}

#[test]
fn accepts_rna_two_position_insertion() {
    assert_accepts("NM_004006.2:r.456_457insu", "NM_004006.2:r.456_457insu");
}

#[test]
fn accepts_noncoding_two_position_insertion() {
    assert_accepts("NR_004430.2:n.789_790insC", "NR_004430.2:n.789_790insC");
}

#[test]
fn accepts_mt_two_position_insertion() {
    assert_accepts("NC_012920.1:m.3243_3244insA", "NC_012920.1:m.3243_3244insA");
}

#[test]
fn accepts_organelle_two_position_insertion() {
    assert_accepts("NC_012920.1:o.3243_3244insA", "NC_012920.1:o.3243_3244insA");
}

// ---- Inside allele brackets the rule applies per-variant ----

#[test]
fn rejects_single_position_insertion_inside_cis_allele() {
    assert_rejects("NM_004006.2:c.[100A>G;456insT]");
}

#[test]
fn accepts_two_position_insertion_inside_cis_allele() {
    assert_accepts(
        "NM_004006.2:c.[100A>G;456_457insT]",
        "NM_004006.2:c.[100A>G;456_457insT]",
    );
}
