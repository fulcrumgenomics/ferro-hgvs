//! Tests for HGVS unknown-phase compound `(;)` separator across every
//! coordinate system. Closes #123 / #81 C3.
//!
//! Per `assets/hgvs-nomenclature/docs/recommendations/DNA/alleles.md` and
//! `recommendations/protein/alleles.md`, the unknown-phase compound
//! `[a(;)b]` (parens around `;`) is symmetric across coordinate systems.
//! Prior to this PR ferro accepted it only for `c.` and `r.`; this file
//! pins parity across `g.`, `n.`, `m.`, `o.`, `p.`.
//!
//! See: docs/superpowers/plans/2026-05-06-issue-123.md

use ferro_hgvs::hgvs::variant::{AllelePhase, HgvsVariant};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Parse and assert the result is an unknown-phase compound allele with the
/// expected number of sub-variants.
fn assert_unknown_phase(input: &str, expected_n_variants: usize) -> HgvsVariant {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("failed to parse {input:?}: {e}"));
    match &v {
        HgvsVariant::Allele(allele) => {
            assert_eq!(
                allele.phase,
                AllelePhase::Unknown,
                "expected unknown phase for {input:?}, got {:?}",
                allele.phase
            );
            assert_eq!(
                allele.variants.len(),
                expected_n_variants,
                "wrong sub-variant count for {input:?}",
            );
        }
        other => panic!("expected Allele variant for {input:?}, got {other:?}"),
    }
    v
}

/// Round-trip parse → Display → parse and assert the second parse equals
/// the first.
fn assert_round_trip(input: &str) {
    let v1 = parse_hgvs(input).unwrap_or_else(|e| panic!("first parse of {input:?} failed: {e}"));
    let display = format!("{v1}");
    let v2 = parse_hgvs(&display)
        .unwrap_or_else(|e| panic!("re-parse of display {display:?} failed: {e}"));
    assert_eq!(
        v1, v2,
        "round-trip mismatch: input={input:?} display={display:?}",
    );
}

fn hash_of<T: Hash>(t: &T) -> u64 {
    let mut h = DefaultHasher::new();
    t.hash(&mut h);
    h.finish()
}

// ---------------------------------------------------------------------------
// Bracketless compact form across all coord systems
// ---------------------------------------------------------------------------

#[test]
fn cds_bracketless_unknown_phase() {
    // c. is already supported; pin it as a parity baseline.
    assert_unknown_phase("NM_004006.2:c.2376G>C(;)3103del", 2);
    assert_round_trip("NM_004006.2:c.2376G>C(;)3103del");
}

#[test]
fn rna_bracketless_unknown_phase() {
    // r. is already supported; baseline.
    assert_unknown_phase("NM_004006.3:r.123c>a(;)345del", 2);
    assert_round_trip("NM_004006.3:r.123c>a(;)345del");
}

#[test]
fn genome_bracketless_unknown_phase() {
    // Spec example (syntax.yaml line 145).
    assert_unknown_phase("NC_000001.11:g.123G>A(;)345del", 2);
    assert_round_trip("NC_000001.11:g.123G>A(;)345del");
}

#[test]
fn tx_bracketless_unknown_phase() {
    // Non-coding transcript.
    assert_unknown_phase("NM_000088.3:n.100A>G(;)200T>C", 2);
    assert_round_trip("NM_000088.3:n.100A>G(;)200T>C");
}

#[test]
fn mt_bracketless_unknown_phase() {
    // Mitochondrial.
    assert_unknown_phase("NC_012920.1:m.100A>G(;)200T>C", 2);
    assert_round_trip("NC_012920.1:m.100A>G(;)200T>C");
}

#[test]
fn circular_bracketless_unknown_phase() {
    // Circular DNA (SVD-WG006).
    assert_unknown_phase("NC_001416.1:o.100A>G(;)200T>C", 2);
    assert_round_trip("NC_001416.1:o.100A>G(;)200T>C");
}

#[test]
fn protein_bracketless_unknown_phase() {
    // Spec example (recommendations/protein/alleles.md line 48).
    // NOTE: per spec, brackets are NOT used for unknown-phase protein.
    assert_unknown_phase("NP_003997.1:p.(Ser68Arg)(;)(Asn594del)", 2);
    assert_round_trip("NP_003997.1:p.(Ser68Arg)(;)(Asn594del)");
}

#[test]
fn protein_bracketless_unknown_phase_homozygous_predicted() {
    // Homozygous predicted form with same variant on both sides.
    assert_unknown_phase("NP_003997.1:p.(Ser68Arg)(;)(Ser68Arg)", 2);
    assert_round_trip("NP_003997.1:p.(Ser68Arg)(;)(Ser68Arg)");
}

// ---------------------------------------------------------------------------
// Bracketed form `[a(;)b]` per coord system (DNA/RNA only — protein uses
// the bare bracketless shape per spec).
// ---------------------------------------------------------------------------

#[test]
fn cds_bracketed_unknown_phase() {
    // Baseline; already supported.
    assert_unknown_phase("NM_004006.2:c.[2376G>C(;)3103del]", 2);
    assert_round_trip("NM_004006.2:c.[2376G>C(;)3103del]");
}

#[test]
fn rna_bracketed_unknown_phase() {
    assert_unknown_phase("NM_004006.3:r.[123c>a(;)345del]", 2);
    assert_round_trip("NM_004006.3:r.[123c>a(;)345del]");
}

#[test]
fn genome_bracketed_unknown_phase() {
    assert_unknown_phase("NC_000001.11:g.[123G>A(;)345del]", 2);
    assert_round_trip("NC_000001.11:g.[123G>A(;)345del]");
}

#[test]
fn tx_bracketed_unknown_phase() {
    assert_unknown_phase("NM_000088.3:n.[100A>G(;)200T>C]", 2);
    assert_round_trip("NM_000088.3:n.[100A>G(;)200T>C]");
}

#[test]
fn mt_bracketed_unknown_phase() {
    assert_unknown_phase("NC_012920.1:m.[100A>G(;)200T>C]", 2);
    assert_round_trip("NC_012920.1:m.[100A>G(;)200T>C]");
}

#[test]
fn circular_bracketed_unknown_phase() {
    assert_unknown_phase("NC_001416.1:o.[100A>G(;)200T>C]", 2);
    assert_round_trip("NC_001416.1:o.[100A>G(;)200T>C]");
}

// ---------------------------------------------------------------------------
// 3+ variant compounds
// ---------------------------------------------------------------------------

#[test]
fn three_variant_unknown_phase_genome() {
    assert_unknown_phase("NC_000001.11:g.100A>G(;)200T>C(;)300del", 3);
    assert_round_trip("NC_000001.11:g.100A>G(;)200T>C(;)300del");
}

#[test]
fn three_variant_unknown_phase_protein() {
    assert_unknown_phase("NP_003997.1:p.(Ser68Arg)(;)(Asn594del)(;)(Pro100Leu)", 3);
    assert_round_trip("NP_003997.1:p.(Ser68Arg)(;)(Asn594del)(;)(Pro100Leu)");
}

#[test]
fn three_variant_unknown_phase_bracketed_genome() {
    assert_unknown_phase("NC_000001.11:g.[100A>G(;)200T>C(;)300del]", 3);
    assert_round_trip("NC_000001.11:g.[100A>G(;)200T>C(;)300del]");
}

// ---------------------------------------------------------------------------
// Phase distinctness: Cis vs Trans vs Unknown all hash and compare distinctly.
// ---------------------------------------------------------------------------

#[test]
fn unknown_phase_distinct_from_cis_genome() {
    let cis = parse_hgvs("NC_000001.11:g.[100A>G;200T>C]").unwrap();
    let unknown = parse_hgvs("NC_000001.11:g.[100A>G(;)200T>C]").unwrap();
    assert_ne!(cis, unknown, "cis and unknown phase must be distinct");
    assert_ne!(
        hash_of(&cis),
        hash_of(&unknown),
        "cis and unknown phase must hash distinctly"
    );
}

#[test]
fn unknown_phase_distinct_from_trans_cds() {
    // Trans-form (`[a];[b]`) parser support across all coord systems is the
    // sister effort tracked in #81 C1 / PR #146 — landed for `c.` and `r.`,
    // pending for `g./n./m./o./p.`. Pin distinctness against `c.` to keep
    // this test orthogonal to the other PR.
    let trans = parse_hgvs("NM_004006.2:c.[2376G>C];[3103del]").unwrap();
    let unknown = parse_hgvs("NM_004006.2:c.[2376G>C(;)3103del]").unwrap();
    assert_ne!(trans, unknown, "trans and unknown phase must be distinct");
    assert_ne!(
        hash_of(&trans),
        hash_of(&unknown),
        "trans and unknown phase must hash distinctly"
    );
}

#[test]
fn unknown_phase_distinct_from_cis_protein() {
    let cis = parse_hgvs("NP_003997.1:p.[Ser68Arg;Asn594del]").unwrap();
    let unknown = parse_hgvs("NP_003997.1:p.(Ser68Arg)(;)(Asn594del)").unwrap();
    assert_ne!(cis, unknown);
    assert_ne!(hash_of(&cis), hash_of(&unknown));
}

#[test]
fn unknown_phase_eq_self_genome() {
    // Eq + Hash stability: two parses of the same string compare and hash equal.
    let a = parse_hgvs("NC_000001.11:g.[100A>G(;)200T>C]").unwrap();
    let b = parse_hgvs("NC_000001.11:g.[100A>G(;)200T>C]").unwrap();
    assert_eq!(a, b);
    assert_eq!(hash_of(&a), hash_of(&b));
}

// ---------------------------------------------------------------------------
// Idempotency under multi-pass normalize.
// ---------------------------------------------------------------------------

fn idempotent_normalize(input: &str) {
    let normalizer = Normalizer::new(MockProvider::new());
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for {input:?}: {e}"));
    let n1 = normalizer.normalize(&v).expect("first normalize");
    let n2 = normalizer.normalize(&n1).expect("second normalize");
    assert_eq!(
        format!("{n1}"),
        format!("{n2}"),
        "normalize is not idempotent for {input:?}",
    );
}

#[test]
fn normalize_idempotent_unknown_phase_genome() {
    idempotent_normalize("NC_000001.11:g.123G>A(;)345del");
    idempotent_normalize("NC_000001.11:g.[123G>A(;)345del]");
}

#[test]
fn normalize_idempotent_unknown_phase_tx() {
    idempotent_normalize("NM_000088.3:n.100A>G(;)200T>C");
    idempotent_normalize("NM_000088.3:n.[100A>G(;)200T>C]");
}

#[test]
fn normalize_idempotent_unknown_phase_mt() {
    idempotent_normalize("NC_012920.1:m.100A>G(;)200T>C");
}

#[test]
fn normalize_idempotent_unknown_phase_circular() {
    idempotent_normalize("NC_001416.1:o.100A>G(;)200T>C");
}

#[test]
fn normalize_idempotent_unknown_phase_protein() {
    idempotent_normalize("NP_003997.1:p.(Ser68Arg)(;)(Asn594del)");
}

// ---------------------------------------------------------------------------
// No-merge contract: an unknown-phase compound of *adjacent* edits MUST NOT
// be merged into a single delins (cis-only contract from #80).
// ---------------------------------------------------------------------------

#[test]
fn unknown_phase_does_not_merge_adjacent_subs_genome() {
    let normalizer = Normalizer::new(MockProvider::new());
    // Cis baseline merges to delins — proves the merge path is wired.
    let cis = parse_hgvs("NC_000001.11:g.[1000G>A;1001A>C]").unwrap();
    let cis_norm = normalizer.normalize(&cis).expect("normalize cis");
    assert_eq!(
        format!("{cis_norm}"),
        "NC_000001.11:g.1000_1001delinsAC",
        "cis-pair sanity (must merge so the unknown-phase assertion is not vacuous)",
    );

    // Unknown-phase parallel input MUST NOT merge.
    let unk = parse_hgvs("NC_000001.11:g.[1000G>A(;)1001A>C]").unwrap();
    let unk_norm = normalizer.normalize(&unk).expect("normalize unknown");
    let s = format!("{unk_norm}");
    assert!(
        s.contains("(;)"),
        "unknown-phase compound must preserve `(;)` (got {s:?})",
    );
    assert!(
        !s.contains("delinsAC"),
        "unknown-phase compound must NOT merge into delins (got {s:?})",
    );
}

// ---------------------------------------------------------------------------
// Mixed-accession compound (ClinVar-style) accepts unknown phase.
// ---------------------------------------------------------------------------

#[test]
fn mixed_accession_unknown_phase() {
    // Two coding variants on different transcripts under unknown-phase. The
    // top-level bracketed form `[A:c.x(;)B:c.y]` is dispatched through
    // `parse_unknown_phase_allele` (no leading shared accession).
    let v = parse_hgvs("[NM_000088.3:c.100A>G(;)NM_000089.4:c.200T>C]").unwrap_or_else(|e| {
        panic!("mixed-accession unknown-phase parse failed: {e}");
    });
    match v {
        HgvsVariant::Allele(allele) => {
            assert_eq!(allele.phase, AllelePhase::Unknown);
            assert_eq!(allele.variants.len(), 2);
        }
        other => panic!("expected Allele, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Malformed rejection
// ---------------------------------------------------------------------------

#[test]
fn malformed_unbalanced_bracket_genome_rejected() {
    // missing closing bracket
    assert!(parse_hgvs("NC_000001.11:g.[100A>G(;)200T>C").is_err());
}

#[test]
fn malformed_trailing_garbage_protein_rejected() {
    assert!(parse_hgvs("NP_003997.1:p.(Ser68Arg)(;)(Asn594del)garbage").is_err());
}

#[test]
fn malformed_empty_segment_unknown_phase_rejected() {
    // Empty `(;)` group inside brackets must hard-fail. Previously the
    // parser silently dropped empty parts and accepted this as a valid
    // 2-variant allele.
    assert!(parse_hgvs("NC_000001.11:g.[100A>G(;)(;)200T>C]").is_err());
}

// Empty `(;)` segments must hard-fail across every coord system, both for
// the bracketed (`[a(;)(;)b]`) and bracketless (`a(;)(;)b`) forms. CodeRabbit
// (PR #148) noted that previously empty parts were silently dropped, so a
// malformed input like `(;)(;)` between two real entries still parsed as a
// valid 2-variant allele.

#[test]
fn malformed_empty_segment_bracketless_rejected_genome() {
    assert!(parse_hgvs("NC_000001.11:g.100A>G(;)(;)200T>C").is_err());
}

#[test]
fn malformed_empty_segment_bracketless_rejected_tx() {
    assert!(parse_hgvs("NR_004430.2:n.100A>G(;)(;)200T>C").is_err());
}

#[test]
fn malformed_empty_segment_bracketless_rejected_mt() {
    assert!(parse_hgvs("NC_012920.1:m.100A>G(;)(;)200T>C").is_err());
}

#[test]
fn malformed_empty_segment_bracketless_rejected_circular() {
    assert!(parse_hgvs("NC_000913.3:o.100A>G(;)(;)200T>C").is_err());
}

#[test]
fn malformed_empty_segment_bracketless_rejected_cds() {
    assert!(parse_hgvs("NM_004006.2:c.100A>G(;)(;)200T>C").is_err());
}

#[test]
fn malformed_empty_segment_bracketless_rejected_rna() {
    assert!(parse_hgvs("NM_004006.3:r.100a>g(;)(;)200u>c").is_err());
}

#[test]
fn malformed_empty_segment_bracketless_rejected_protein() {
    // Bracketless protein form: each sub-variant is `(...)`; empty `(;)`
    // group between them must error.
    assert!(parse_hgvs("NP_003997.1:p.(Ser68Arg)(;)(;)(Asn594del)").is_err());
}

#[test]
fn malformed_empty_segment_bracketed_rejected_tx() {
    assert!(parse_hgvs("NR_004430.2:n.[100A>G(;)(;)200T>C]").is_err());
}

#[test]
fn malformed_empty_segment_bracketed_rejected_protein() {
    assert!(parse_hgvs("NP_003997.1:p.[Ser68Arg(;)(;)Asn594del]").is_err());
}

// Trailing/leading empty `(;)` groups inside brackets (e.g. `[a(;)]`) must
// error rather than silently produce a phantom 1-variant unknown-phase allele.
// Note: bare singletons like `[a]` are intentionally accepted by ferro and
// expanded upstream — see test_singleton_cis_allele_preserves_wrapper.

#[test]
fn malformed_trailing_empty_unknown_group_rejected_genome() {
    assert!(parse_hgvs("NC_000001.11:g.[100A>G(;)]").is_err());
}

#[test]
fn malformed_trailing_empty_unknown_group_rejected_tx() {
    assert!(parse_hgvs("NR_004430.2:n.[100A>G(;)]").is_err());
}

#[test]
fn malformed_trailing_empty_unknown_group_rejected_protein() {
    assert!(parse_hgvs("NP_003997.1:p.[Ser68Arg(;)]").is_err());
}

// ---------------------------------------------------------------------------
// Empty `(;)` segments must be rejected on the mixed-accession path too
// (CodeRabbit PR #148): `parse_unknown_phase_allele` previously skipped empty
// splits, so `[A:c.x(;)(;)B:c.y]` parsed as a valid 2-variant allele.
// ---------------------------------------------------------------------------

#[test]
fn malformed_empty_segment_mixed_accession_rejected() {
    assert!(
        parse_hgvs("[NM_000088.3:c.100A>G(;)(;)NM_000089.4:c.200T>C]").is_err(),
        "double `(;)` in a mixed-accession unknown-phase compound must be rejected",
    );
}

// ---------------------------------------------------------------------------
// Mixed `;`/`(;)` inside one bracket pair (CodeRabbit PR #148): per HGVS spec
// (`recommendations/{DNA,protein}/alleles.md`), a single `[...]` carries either
// all-cis (`;`) or all-unknown-phase (`(;)`) members. Top-level mixing
// (`[a;b];[c](;)d`) uses the trans + bracketless paths separately. Mixing
// inside one bracket is not spec-valid — must be rejected, not flattened.
// ---------------------------------------------------------------------------

#[test]
fn malformed_mixed_separators_rejected_genome() {
    assert!(parse_hgvs("NC_000001.11:g.[100A>G;200T>C(;)300del]").is_err());
}

#[test]
fn malformed_mixed_separators_rejected_tx() {
    assert!(parse_hgvs("NM_000088.3:n.[100A>G;200T>C(;)300del]").is_err());
}

#[test]
fn malformed_mixed_separators_rejected_protein() {
    assert!(parse_hgvs("NP_003997.1:p.[Ser68Arg;Asn594del(;)Pro100Leu]").is_err());
}

// ---------------------------------------------------------------------------
// Compound allele members containing nested bracketed edits (CodeRabbit
// PR #148): the compound parser must locate the FIRST TOP-LEVEL `]`, not the
// first byte match — otherwise valid members like `delins[ACC:g.x_y]`,
// `delins[A;T]`, etc. are truncated mid-edit and the whole compound rejects.
// ---------------------------------------------------------------------------

#[test]
fn compound_member_with_xref_delins_genome() {
    let v = parse_hgvs("NC_000008.11:g.[86587460_86650711delins[KY923049.1:g.1_466];86660000A>G]")
        .expect("compound with delins[xref] member must parse");
    match v {
        HgvsVariant::Allele(allele) => {
            assert_eq!(allele.variants.len(), 2);
        }
        other => panic!("expected Allele, got {other:?}"),
    }
}

#[test]
fn compound_member_with_inner_bracket_delins_tx() {
    let v = parse_hgvs("NM_002001.2:c.[12_17delins[28_39];100A>G]")
        .expect("tx compound with delins[range] member must parse");
    match v {
        HgvsVariant::Allele(allele) => {
            assert_eq!(allele.variants.len(), 2);
        }
        other => panic!("expected Allele, got {other:?}"),
    }
}

#[test]
fn compound_member_with_inner_bracket_unknown_phase_genome() {
    // Same nested-bracket hazard but with `(;)` separator.
    let v =
        parse_hgvs("NC_000008.11:g.[86587460_86650711delins[KY923049.1:g.1_466](;)86660000A>G]")
            .expect("unknown-phase compound with delins[xref] member must parse");
    match v {
        HgvsVariant::Allele(allele) => {
            assert_eq!(allele.phase, AllelePhase::Unknown);
            assert_eq!(allele.variants.len(), 2);
        }
        other => panic!("expected Allele, got {other:?}"),
    }
}
