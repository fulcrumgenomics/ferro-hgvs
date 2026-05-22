//! Regression locks for canonical HGVS v21 spec examples that aren't yet
//! pinned by ferro's test suite.
//!
//! Source: phase4_final_triage.md families F3 + F4 + F5 (13 net-new corners,
//! all POS-style — the spec gives an explicit canonical form and we want a
//! test to catch any future drift). Each test names the construct,
//! coordinate system, and the originating axes-doc corner; the spec citation
//! is in the per-test comment.
//!
//! Every test in this file pins parse + Display round-trip against the
//! spec-canonical form. If a future change regresses any input to a parse
//! error or to a non-canonical Display, the corresponding test fails — that
//! is the entire point of the file.

use ferro_hgvs::parse_hgvs;

// =====================================================================
// F3 — CNV / large-range / uncertain endpoint pairs (7 corners)
// =====================================================================
//
// `DNA/duplication.md` gives several canonical-shape examples that combine
// intronic offsets, uncertain endpoints, and large ranges. SVD-WG003
// proposed the nested-paren uncertain endpoint form
// `(A+1_B-1)_(C+1_D-1)dup`; ferro parses + round-trips it.

/// Exon/intron border duplication: spec uses `c.1704+1dup` (the intronic
/// `+1` position), NOT `c.1704dup`. Source: dna-duplication:C-C3.
#[test]
fn dup_exon_intron_border_round_trips() {
    let input = "NM_000088.3:c.1704+1dup";
    match parse_hgvs(input) {
        Ok(v) => assert_eq!(v.to_string(), input, "round-trip mismatch"),
        Err(e) => panic!("ferro should parse spec-canonical {input:?}: {e}"),
    }
}

/// Intron/exon border duplication: spec uses `c.1813dup`, NOT
/// `c.1813-1dup`. Source: dna-duplication:C-C4.
#[test]
fn dup_intron_exon_border_round_trips() {
    let input = "NM_000088.3:c.1813dup";
    match parse_hgvs(input) {
        Ok(v) => assert_eq!(v.to_string(), input, "round-trip mismatch"),
        Err(e) => panic!("ferro should parse spec-canonical {input:?}: {e}"),
    }
}

/// Two-exon-spanning duplication with intronic endpoints. Source:
/// dna-duplication:C-C5.
#[test]
fn dup_two_exon_spanning_intronic_endpoints_round_trips() {
    let input = "NM_000088.3:c.4072-1234_5155-246dup";
    let parsed = parse_hgvs(input).expect("must parse");
    assert_eq!(parsed.to_string(), input, "round-trip mismatch");
}

/// SVD-WG003 uncertain endpoint pair form. Source: dna-duplication:C-C6.
#[test]
fn dup_uncertain_endpoint_pairs_svd_wg003_round_trips() {
    let input = "NM_000088.3:c.(4071+1_4072-1)_(5154+1_5155-1)dup";
    let parsed = parse_hgvs(input).expect("must parse");
    assert_eq!(
        parsed.to_string(),
        input,
        "round-trip mismatch for SVD-WG003 uncertain endpoint pair"
    );
}

/// Triplication CNV with uncertain endpoint pairs. Source:
/// dna-duplication:C-C7.
#[test]
fn cnv_triplication_uncertain_endpoints_round_trips() {
    let input = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)[3]";
    let parsed = parse_hgvs(input).expect("must parse");
    assert_eq!(parsed.to_string(), input, "round-trip mismatch");
}

/// Whole-gene CNV duplication form. Source: dna-duplication:C-C8.
#[test]
fn cnv_whole_gene_dup_round_trips() {
    let input = "NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup";
    let parsed = parse_hgvs(input).expect("must parse");
    assert_eq!(parsed.to_string(), input, "round-trip mismatch");
}

/// Very large inversion (megabase scale). Source: dna-inversion:C-C6.
#[test]
fn large_inversion_round_trips() {
    let input = "NC_000023.11:g.111754331_111966764inv";
    match parse_hgvs(input) {
        Ok(v) => assert_eq!(v.to_string(), input),
        Err(e) => panic!("ferro should parse large inversion {input:?}: {e}"),
    }
}

// =====================================================================
// F4 — Combined-insertion bracketed forms (4 corners)
// =====================================================================
//
// `DNA/inversion.md` shows insertions whose payload is a bracketed
// composition of inv / sub / del. These are spec-canonical for describing
// the inverted-copy and complex-rearrangement cases.

/// Inverted-copy at 5' boundary of original. Source: dna-inversion:C-C1.
/// Spec canonical: `ins123_234inv` without brackets (DNA/insertion.md:18,
/// DNA/inversion.md:39). Locks the fix for the single-payload Display.
#[test]
fn ins_inverted_copy_5prime_boundary_round_trips() {
    let input = "NC_000023.11:g.122_123ins123_234inv";
    let parsed = parse_hgvs(input).expect("must parse");
    assert_eq!(
        parsed.to_string(),
        input,
        "single-payload ins<range>inv must not be bracketed (spec form)"
    );
}

/// Inverted-copy at 3' boundary of original. Source: dna-inversion:C-C2.
#[test]
fn ins_inverted_copy_3prime_boundary_round_trips() {
    let input = "NC_000023.11:g.234_235ins123_234inv";
    let parsed = parse_hgvs(input).expect("must parse");
    assert_eq!(
        parsed.to_string(),
        input,
        "single-payload ins<range>inv must not be bracketed (spec form)"
    );
}

/// Combination inserted form: two inversions + an inserted A.
/// Source: dna-inversion:C-C3.
#[test]
fn ins_combination_inv_sub_inv_round_trips() {
    let input = "NM_004006.2:c.940_941ins[885_940inv;A;851_883inv]";
    let parsed = parse_hgvs(input).expect("must parse");
    assert_eq!(parsed.to_string(), input, "round-trip mismatch");
}

/// Combination inserted form: two inversions composed.
/// Source: dna-inversion:C-C4.
#[test]
fn ins_combination_two_inversions_round_trips() {
    let input = "NM_004006.2:c.940_941ins[903_940inv;851_885inv]";
    let parsed = parse_hgvs(input).expect("must parse");
    assert_eq!(parsed.to_string(), input, "round-trip mismatch");
}

// =====================================================================
// F5 — Boundary-stacking subs (UTR ∩ intronic) (2 corners)
// =====================================================================
//
// Substitutions at positions that span two boundary classes simultaneously
// (e.g. 5' UTR + intronic donor). `numbering.md` defines the position
// kinds; the substitution doc doesn't cross them with a worked example.

/// Substitution at 5' UTR intronic donor. Source: dna-substitution:C-C3.
#[test]
fn sub_5utr_intronic_donor_round_trips() {
    let input = "NM_004006.2:c.-85+1A>G";
    match parse_hgvs(input) {
        Ok(v) => assert_eq!(v.to_string(), input, "round-trip mismatch"),
        Err(e) => panic!("ferro should parse 5'UTR-intronic sub {input:?}: {e}"),
    }
}

/// Substitution at 5' UTR intronic donor (variant position).
/// Source: dna-substitution:C-C3.
#[test]
fn sub_5utr_intronic_donor_variant_round_trips() {
    let input = "NM_004006.2:c.-15+1A>G";
    match parse_hgvs(input) {
        Ok(v) => assert_eq!(v.to_string(), input, "round-trip mismatch"),
        Err(e) => panic!("ferro should parse {input:?}: {e}"),
    }
}

/// Substitution at first intronic nt of intron immediately after stop.
/// Source: dna-substitution:C-C8.
#[test]
fn sub_stop_plus_1_intronic_acceptor_round_trips() {
    let input = "NM_004006.2:c.*1-1A>G";
    match parse_hgvs(input) {
        Ok(v) => assert_eq!(v.to_string(), input, "round-trip mismatch"),
        Err(e) => panic!("ferro should parse stop+1 intronic sub {input:?}: {e}"),
    }
}

/// Substitution at deep 3'UTR intronic. Source: dna-substitution:C-C8.
#[test]
fn sub_3utr_deep_intronic_round_trips() {
    let input = "NM_004006.2:c.*639-1G>A";
    match parse_hgvs(input) {
        Ok(v) => assert_eq!(v.to_string(), input, "round-trip mismatch"),
        Err(e) => panic!("ferro should parse 3'UTR deep intronic {input:?}: {e}"),
    }
}
