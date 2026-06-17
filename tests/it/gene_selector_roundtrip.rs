//! Round-trip tests for the optional gene-symbol selector.
//!
//! HGVS Nomenclature allows an informational gene-symbol selector in
//! parentheses after the accession to disambiguate when the same accession
//! spans multiple genes (e.g. `NM_000088.3(COL1A1):c.459A>G`).
//!
//! PR #70 (`b34ce35`) extended the parser so the selector is accepted on
//! non-RefSeq accessions and is captured into the variant struct's
//! `gene_symbol: Option<String>`. PR #144 (closes #121) extended Display
//! to preserve the selector across every concrete variant kind and every
//! compound-allele form. This file pins the *end-to-end* round-trip
//! identity (parse → normalize → Display → reparse) across the full
//! audit matrix: every accession family, every coordinate system, every
//! allele compound form, every LRG flavor, and the `Hash` + `Eq`
//! contract.
//!
//! Spec rationale
//! --------------
//! Per <https://hgvs-nomenclature.org>, the gene-symbol selector is
//! informational disambiguation. ferro's round-trip policy is
//! **preserve when present in input; do not synthesize when absent** —
//! aligning ferro with Mutalyzer and VariantValidator.
//!
//! Scope split with `tests/gene_selector_display_preserve.rs` (#144):
//! - That file targets the *formatter*: assert each per-kind `Display`
//!   impl emits `accession(gene)` when `gene_symbol = Some(g)` and the
//!   bare form when `None`.
//! - This file targets *end-to-end identity*: parse + normalize +
//!   Display + reparse, including the full coordinate-system and
//!   allele-form matrix and the negative half of the policy
//!   (no synthesis, multi-pass-normalize idempotency, Hash/Eq
//!   distinctness on round-trip).

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, HgvsVariant, Normalizer};

/// Return the parsed variant's `gene_symbol`, regardless of variant kind.
fn gene_symbol_of(variant: &HgvsVariant) -> Option<&str> {
    match variant {
        HgvsVariant::Genome(v) => v.gene_symbol.as_deref(),
        HgvsVariant::Cds(v) => v.gene_symbol.as_deref(),
        HgvsVariant::Tx(v) => v.gene_symbol.as_deref(),
        HgvsVariant::Rna(v) => v.gene_symbol.as_deref(),
        HgvsVariant::Protein(v) => v.gene_symbol.as_deref(),
        HgvsVariant::Mt(v) => v.gene_symbol.as_deref(),
        HgvsVariant::Circular(v) => v.gene_symbol.as_deref(),
        HgvsVariant::GenomeRing(g) => g.gene_symbol.as_deref(),
        HgvsVariant::Supernumerary(inner) => inner.gene_symbol(),
        // RnaFusion, Allele: selectors live on sub-variants/breakpoints.
        HgvsVariant::RnaFusion(_)
        | HgvsVariant::Allele(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// Assert a parse → display → reparse cycle is byte-stable and the selector
/// survives both the first parse and the reparse.
///
/// Use this for any input that should yield a non-`Allele` variant carrying
/// `gene_symbol = Some(g)`. For nested-kind variants (`Allele`, `RnaFusion`)
/// the top-level `gene_symbol` is `None` by design — assert the round-trip
/// directly without this helper.
fn assert_display_preserves_gene_selector(input: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    let expected_gene = gene_symbol_of(&v)
        .unwrap_or_else(|| panic!("test bug: input {input:?} must have a gene_symbol after parse"))
        .to_string();

    let displayed = v.to_string();
    assert_eq!(
        displayed, input,
        "Display must preserve the gene-symbol selector verbatim (#121)",
    );

    let reparsed =
        parse_hgvs(&displayed).unwrap_or_else(|e| panic!("reparse {displayed:?} failed: {e:?}"));
    assert_eq!(reparsed.to_string(), displayed, "Display is not idempotent");
    assert_eq!(
        gene_symbol_of(&reparsed).map(str::to_string),
        Some(expected_gene),
        "reparse must observe the preserved gene_symbol",
    );
}

// =============================================================================
// Parse: gene_symbol is captured on every supported accession family.
// =============================================================================

#[test]
fn parse_captures_gene_symbol_refseq_nm() {
    // Canonical example from the HGVS Nomenclature page.
    let v = parse_hgvs("NM_000088.3(COL1A1):c.459A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("COL1A1"));
    assert_eq!(v.accession().unwrap().full(), "NM_000088.3");
}

#[test]
fn parse_captures_gene_symbol_refseq_nr() {
    // NR_ is a non-coding RNA RefSeq accession; coordinate type is `n.`.
    let v = parse_hgvs("NR_046018.2(DDX11L1):n.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("DDX11L1"));
    assert_eq!(v.accession().unwrap().full(), "NR_046018.2");
}

#[test]
fn parse_captures_gene_symbol_refseq_np() {
    // Protein accession with predicted protein change uses `p.(...)` form.
    let v = parse_hgvs("NP_000079.2(COL1A1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("COL1A1"));
    assert_eq!(v.accession().unwrap().full(), "NP_000079.2");
}

#[test]
fn parse_captures_gene_symbol_ensembl() {
    let v = parse_hgvs("ENST00000380152.7(BRCA2):c.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("BRCA2"));
    assert_eq!(v.accession().unwrap().full(), "ENST00000380152.7");
}

#[test]
fn parse_captures_gene_symbol_lrg() {
    // LRG_<n>t<m> denotes an LRG transcript reference.
    let v = parse_hgvs("LRG_199t1(BRCA1):c.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("BRCA1"));
}

#[test]
fn parse_captures_gene_symbol_simple_accession() {
    // PR #70's motivating case: a custom (non-RefSeq, non-Ensembl, non-LRG)
    // accession with a gene-symbol selector must round-trip the selector.
    let v = parse_hgvs("MYREF_SEQ(GENE1):c.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("GENE1"));
}

#[test]
fn parse_without_selector_leaves_gene_symbol_none() {
    // The selector must not be synthesized when absent from input.
    let v = parse_hgvs("NM_000088.3:c.459A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), None);
}

// =============================================================================
// Display: every supported accession family preserves the selector verbatim
// (#121). Replaces the pre-#144 strip-pinned cases.
// =============================================================================

#[test]
fn display_preserves_gene_selector_refseq_nm() {
    assert_display_preserves_gene_selector("NM_000088.3(COL1A1):c.459A>G");
}

#[test]
fn display_preserves_gene_selector_refseq_nr() {
    assert_display_preserves_gene_selector("NR_046018.2(DDX11L1):n.100A>G");
}

#[test]
fn display_drops_gene_selector_refseq_np() {
    // Per HGVS syntax.yaml 119–128 the `accession(selector):p.` form is not
    // part of the protein-variant grammar; the parenthesized selector is only
    // defined for genomic-to-transcript projection on c./n. (syntax.yaml 213).
    // Parsing remains permissive — `gene_symbol` round-trips on the struct —
    // but Display emits the spec-compliant `accession:p.` form. See #310.
    let v = parse_hgvs("NP_000079.2(COL1A1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("COL1A1"));
    assert_eq!(v.to_string(), "NP_000079.2:p.(Arg8Gln)");
}

#[test]
fn display_preserves_gene_selector_ensembl() {
    assert_display_preserves_gene_selector("ENST00000380152.7(BRCA2):c.100A>G");
}

#[test]
fn display_preserves_gene_selector_simple_accession() {
    assert_display_preserves_gene_selector("MYREF_SEQ(GENE1):c.100A>G");
}

#[test]
fn display_no_selector_unchanged() {
    // Sanity check: Display is idempotent when the selector is absent and
    // is not synthesized from somewhere else (e.g. a transcript provider
    // lookup). This is the negative half of the #121 round-trip policy.
    let input = "NM_000088.3:c.459A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), input);
    assert_eq!(gene_symbol_of(&v), None);
}

// =============================================================================
// Coordinate-system matrix: every supported `:type.` value round-trips with a
// selector. `c`, `n`, `p` are exercised in the section above; this section
// adds `g`, `r`, `m`, `o`.
// =============================================================================

#[test]
fn display_preserves_gene_selector_genome() {
    // `(FLT3)` here is a gene-symbol selector, not a compound-ref wrapper —
    // a single-segment accession with a selector at the same position.
    assert_display_preserves_gene_selector("NC_000013.11(FLT3):g.12345A>G");
}

#[test]
fn display_preserves_gene_selector_rna() {
    // `r.` form: lowercase nucleotides, otherwise selector handling is uniform.
    assert_display_preserves_gene_selector("NM_000088.3(COL1A1):r.100a>g");
}

#[test]
fn display_preserves_gene_selector_mitochondrial() {
    // `m.` form: spec explicitly endorses `NC_012920.1(MT-…):m.…`
    // (see assets/hgvs-nomenclature/docs/background/refseq.md lines 197-203).
    // This case flipped from `diverges` to `preserved` in PR #144's spec
    // fixture; pinning it here gives the audit a direct end-to-end anchor.
    assert_display_preserves_gene_selector("NC_012920.1(MT-ND1):m.3460G>A");
}

#[test]
fn display_preserves_gene_selector_circular() {
    // `o.` form (circular DNA, SVD-WG006). ferro accepts the selector
    // verbatim across this coordinate type too.
    assert_display_preserves_gene_selector("NC_001416.1(GENE):o.1A>G");
}

// =============================================================================
// LRG flavor matrix: bare (`LRG_<n>` -> g.), transcript (`LRG_<n>t<m>` -> c.,
// also valid with n. and r.), protein (`LRG_<n>p<m>` -> p.). Each form must
// round-trip the selector.
// =============================================================================

#[test]
fn display_preserves_gene_selector_lrg_bare_genome() {
    assert_display_preserves_gene_selector("LRG_199(BRCA1):g.100A>G");
}

#[test]
fn display_preserves_gene_selector_lrg_transcript_cds() {
    // Already covered above; included here for completeness within the LRG
    // matrix so a future LRG regression surfaces in this section directly.
    assert_display_preserves_gene_selector("LRG_199t1(BRCA1):c.100A>G");
}

#[test]
fn display_preserves_gene_selector_lrg_transcript_noncoding() {
    assert_display_preserves_gene_selector("LRG_199t1(BRCA1):n.100A>G");
}

#[test]
fn display_preserves_gene_selector_lrg_transcript_rna() {
    assert_display_preserves_gene_selector("LRG_199t1(BRCA1):r.100a>g");
}

#[test]
fn display_drops_gene_selector_lrg_protein() {
    // Per HGVS syntax.yaml 119–128 the protein-variant grammar has no
    // parenthesized selector position; Display emits the spec-compliant
    // `LRG_*p*:p.` form even when the parser accepted a `(GENE)` selector.
    // See #310.
    let v = parse_hgvs("LRG_199p1(BRCA1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("BRCA1"));
    assert_eq!(v.to_string(), "LRG_199p1:p.(Arg8Gln)");
}

// =============================================================================
// Compound-ref + selector: `NC_X(NM_Y)(GENE):c.…` — a genomic-context
// accession (`NC(NM)`) followed by a gene selector. Both must round-trip
// independently. Bare-compound-ref negative case pins that no selector is
// synthesized for the inner accession.
// =============================================================================

#[test]
fn display_preserves_gene_selector_with_compound_ref() {
    assert_display_preserves_gene_selector("NC_000013.11(NM_004119.3)(FLT3):c.100A>G");
}

#[test]
fn display_compound_ref_without_selector_unchanged() {
    let input = "NC_000013.11(NM_004119.3):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), None);
    assert_eq!(v.to_string(), input);
}

// =============================================================================
// Allele compound-form matrix: the compact form (`ACC(GENE):c.[a;b]`,
// `ACC(GENE):c.[a];[b]`, `ACC(GENE):c.a(;)b`) carries the selector once on
// the prefix; the expanded form (`[ACC(GENE):c.a];[ACC(GENE):c.b]`) carries
// it per-sub-variant. PR #144 tightened `all_share_accession_and_type` so
// the compact form is only used when every sub-variant agrees on
// `gene_symbol`.
// =============================================================================

#[test]
fn display_preserves_gene_selector_in_compact_cis_allele() {
    // Compact cis: `ACC(GENE):c.[a;b]`. Selector emitted once on the prefix.
    assert_display_preserves_gene_selector_in_allele(
        "NM_000088.3(COL1A1):c.[100A>G;200C>T]",
        &[Some("COL1A1"), Some("COL1A1")],
    );
}

#[test]
fn display_preserves_gene_selector_in_compact_trans_allele() {
    // Compact trans: `ACC(GENE):c.[a];[b]`. Selector emitted once on the prefix.
    assert_display_preserves_gene_selector_in_allele(
        "NM_000088.3(COL1A1):c.[100A>G];[200C>T]",
        &[Some("COL1A1"), Some("COL1A1")],
    );
}

#[test]
fn display_preserves_gene_selector_in_unknown_phase_compact_allele() {
    // Unknown-phase compact: `ACC(GENE):c.a(;)b`. No surrounding brackets;
    // selector lives on the prefix once.
    assert_display_preserves_gene_selector_in_allele(
        "NM_000088.3(COL1A1):c.100A>G(;)200C>T",
        &[Some("COL1A1"), Some("COL1A1")],
    );
}

#[test]
fn display_canonicalizes_expanded_same_acc_trans_to_compact_with_selector() {
    // Same accession + same gene + trans phase: ferro canonicalizes the
    // expanded `[v1];[v2]` input to the compact-trans form
    // `ACC(GENE):c.[a];[b]` on Display, since `all_share_accession_and_type`
    // approves the consensus. The selector survives on every sub-variant
    // and on the compact prefix; reparse is byte-stable thereafter.
    let input = "[NM_000088.3(COL1A1):c.100A>G];[NM_000088.3(COL1A1):c.200C>T]";
    let v = parse_hgvs(input).unwrap();
    let displayed = v.to_string();
    assert_eq!(
        displayed, "NM_000088.3(COL1A1):c.[100A>G];[200C>T]",
        "matching selectors must canonicalize to compact-trans form",
    );

    let reparsed = parse_hgvs(&displayed).unwrap();
    assert_eq!(
        reparsed.to_string(),
        displayed,
        "reparse of compact-trans must be byte-stable"
    );

    if let HgvsVariant::Allele(a) = &v {
        assert_eq!(a.variants.len(), 2);
        for sub in &a.variants {
            assert_eq!(gene_symbol_of(sub), Some("COL1A1"));
        }
    } else {
        panic!("expected Allele; got {:?}", v);
    }
}

#[test]
fn display_preserves_independent_gene_selectors_in_mixed_acc_trans() {
    // Different accessions with different gene symbols — the canonical
    // compound-heterozygosity report shape. Per #121 the expanded form
    // preserves both selectors, one per sub-variant.
    let input = "[NM_000088.3(COL1A1):c.100A>G];[NM_000089.3(COL1A2):c.200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), input);

    if let HgvsVariant::Allele(a) = &v {
        assert_eq!(a.variants.len(), 2);
        assert_eq!(gene_symbol_of(&a.variants[0]), Some("COL1A1"));
        assert_eq!(gene_symbol_of(&a.variants[1]), Some("COL1A2"));
    } else {
        panic!("expected Allele; got {:?}", v);
    }

    let reparsed = parse_hgvs(&v.to_string()).unwrap();
    assert_eq!(reparsed, v, "structural identity across reparse");
}

#[test]
fn display_falls_back_to_expanded_when_gene_symbols_mismatch() {
    // End-to-end pin of the `all_share_accession_and_type` consensus rule
    // tightened in #144: same accession but mismatched gene_symbol →
    // compact form would lose information → emit expanded form so each
    // sub-variant's selector survives. Parsing the expanded form and
    // round-tripping must produce the same expanded form.
    let input = "[NM_000088.3(COL1A1):c.100A>G];[NM_000088.3(COL1A2):c.200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(
        v.to_string(),
        input,
        "mismatched selectors must keep expanded form, not collapse to compact"
    );

    if let HgvsVariant::Allele(a) = &v {
        assert_eq!(a.variants.len(), 2);
        assert_eq!(gene_symbol_of(&a.variants[0]), Some("COL1A1"));
        assert_eq!(gene_symbol_of(&a.variants[1]), Some("COL1A2"));
    } else {
        panic!("expected Allele; got {:?}", v);
    }
}

/// Helper for allele-form tests: assert byte-stable round-trip plus per-sub
/// gene_symbol expectations.
fn assert_display_preserves_gene_selector_in_allele(input: &str, expected: &[Option<&str>]) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    let displayed = v.to_string();
    assert_eq!(displayed, input, "Display must preserve verbatim");

    if let HgvsVariant::Allele(a) = &v {
        assert_eq!(
            a.variants.len(),
            expected.len(),
            "expected {} sub-variants for {input}",
            expected.len(),
        );
        for (i, (sub, want)) in a.variants.iter().zip(expected).enumerate() {
            assert_eq!(
                gene_symbol_of(sub),
                *want,
                "sub-variant {i} gene_symbol mismatch for {input}",
            );
        }
    } else {
        panic!("expected Allele for {input}; got {:?}", v);
    }

    let reparsed = parse_hgvs(&displayed).unwrap();
    assert_eq!(
        reparsed, v,
        "structural identity across reparse for {input}"
    );
}

// =============================================================================
// Edge inputs: empty `()`, whitespace padding, casing, hyphen. ferro's
// policy is "preserve verbatim" — the parser does not sanitize gene_symbol,
// and Display emits whatever was captured, except `()` collapses to None
// per the parser's "treat empty strings as None" rule.
// =============================================================================

#[test]
fn empty_selector_round_trips_to_bare_form() {
    // `NM_X():c.…` parses with gene_symbol = None; Display emits the bare
    // form (does not synthesize an empty `()`).
    let input = "NM_000088.3():c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), None);

    let displayed = v.to_string();
    assert_eq!(
        displayed, "NM_000088.3:c.100A>G",
        "empty selector must drop to bare form, not synthesize `()`"
    );

    let reparsed = parse_hgvs(&displayed).unwrap();
    assert_eq!(reparsed.to_string(), displayed);
    assert_eq!(gene_symbol_of(&reparsed), None);
}

#[test]
fn whitespace_padded_selector_round_trips_verbatim() {
    // ferro is bytes-in / bytes-out: the spec uses HGNC symbols (no
    // whitespace), but ferro does not sanitize so callers can detect
    // malformed input.
    let input = "NM_000088.3( COL1A1 ):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), Some(" COL1A1 "));
    assert_eq!(v.to_string(), input);

    let reparsed = parse_hgvs(&v.to_string()).unwrap();
    assert_eq!(reparsed, v);
}

#[test]
fn lowercase_selector_round_trips_verbatim() {
    // ferro does not normalize case on gene_symbol.
    let input = "NM_000088.3(col1a1):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), Some("col1a1"));
    assert_eq!(v.to_string(), input);

    let reparsed = parse_hgvs(&v.to_string()).unwrap();
    assert_eq!(reparsed, v);
}

#[test]
fn hyphenated_selector_round_trips() {
    // Mitochondrial genes carry hyphens (MT-ND1, MT-TL1). The selector
    // parser uses `take_while(c != ')')`, so any non-paren bytes are
    // accepted verbatim. The spec endorses this exact form (see
    // assets/hgvs-nomenclature/docs/background/refseq.md lines 197-203).
    assert_display_preserves_gene_selector("NC_012920.1(MT-TL1):m.3243A>G");
    assert_display_preserves_gene_selector("NC_012920.1(MT-TL1):n.14A>G");
}

// =============================================================================
// Normalize: the gene_symbol field is preserved through normalization, even
// while Display now also preserves it on output. Multi-pass idempotency
// pins that gene_symbol does not drift on repeated normalize calls.
// =============================================================================

#[test]
fn normalize_preserves_gene_symbol_field_cds() {
    // MockProvider::with_test_data has NM_000088.3 (COL1A1), so this
    // exercises the real CDS normalize path.
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3(COL1A1):c.4A>G").unwrap();
    assert_eq!(gene_symbol_of(&variant), Some("COL1A1"));

    let normalized = normalizer.normalize(&variant).unwrap();

    // The gene_symbol field must survive normalization, regardless of any
    // shifting the normalizer might do to the position.
    assert_eq!(
        gene_symbol_of(&normalized),
        Some("COL1A1"),
        "normalize must not drop gene_symbol on CdsVariant",
    );

    // And Display must continue to emit it after normalize.
    assert!(
        normalized.to_string().contains("(COL1A1)"),
        "normalized Display must preserve the gene selector: {}",
        normalized
    );
}

#[test]
fn normalize_preserves_gene_symbol_field_genome() {
    // Genomic substitutions normalize without a transcript lookup.
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NC_000017.11(BRCA1):g.43094466G>A").unwrap();
    assert_eq!(gene_symbol_of(&variant), Some("BRCA1"));

    let normalized = normalizer.normalize(&variant).unwrap();
    assert_eq!(
        gene_symbol_of(&normalized),
        Some("BRCA1"),
        "normalize must not drop gene_symbol on GenomeVariant",
    );
}

#[test]
fn normalize_does_not_synthesize_gene_symbol() {
    // When the input lacked a selector, normalization must not invent one
    // (even though the underlying transcript record knows the gene symbol),
    // and Display must continue to emit the bare form.
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3:c.4A>G").unwrap();
    assert_eq!(gene_symbol_of(&variant), None);

    let normalized = normalizer.normalize(&variant).unwrap();
    assert_eq!(
        gene_symbol_of(&normalized),
        None,
        "normalize must not synthesize gene_symbol from transcript metadata",
    );
    assert_eq!(
        normalized.to_string(),
        "NM_000088.3:c.4A>G",
        "Display must not synthesize a gene selector after normalize",
    );
}

#[test]
fn normalize_is_idempotent_with_gene_symbol() {
    // Multi-pass idempotency: normalize-of-normalize must be a fixed point
    // under the derived `Eq`. This pins the contract so a future change
    // that re-routes normalization through a transcript provider (which
    // knows the gene symbol and could synthesize one) does not silently
    // alter `gene_symbol`.
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let cases = ["NM_000088.3(COL1A1):c.4A>G", "NM_000088.3:c.4A>G"];
    for input in cases {
        let v = parse_hgvs(input).unwrap();
        let n1 = normalizer.normalize(&v).unwrap();
        let n2 = normalizer.normalize(&n1).unwrap();
        assert_eq!(
            n1, n2,
            "normalize must be idempotent under derived Eq for {input}"
        );
        assert_eq!(
            gene_symbol_of(&n1),
            gene_symbol_of(&n2),
            "normalize must not change gene_symbol on a second pass for {input}"
        );
    }
}

// =============================================================================
// Re-parse identity: parse → display → parse a second time produces a
// structurally equal variant — the second parse observes Some(<gene>),
// not None. Pins that the round-trip is identity-preserving for the
// selector.
// =============================================================================

#[test]
fn reparse_preserves_gene_symbol() {
    // Only non-protein accession families round-trip byte-stably — `p.`
    // Display intentionally drops the gene-symbol selector (#310), so a
    // protein variant is not identity-preserving across Display in this
    // sense and is covered separately in `display_drops_gene_selector_refseq_np`.
    let cases: &[(&str, &str)] = &[
        ("NM_000088.3(COL1A1):c.459A>G", "COL1A1"),
        ("NR_046018.2(DDX11L1):n.100A>G", "DDX11L1"),
        ("ENST00000380152.7(BRCA2):c.100A>G", "BRCA2"),
        ("LRG_199t1(BRCA1):c.100A>G", "BRCA1"),
        ("MYREF_SEQ(GENE1):c.100A>G", "GENE1"),
    ];

    for (input, gene) in cases {
        let first = parse_hgvs(input).unwrap();
        let displayed = first.to_string();
        assert_eq!(displayed, *input, "Display must be byte-stable for {input}");

        let second = parse_hgvs(&displayed).unwrap();
        assert_eq!(
            second.to_string(),
            displayed,
            "reparse not idempotent: {input}"
        );
        assert_eq!(
            gene_symbol_of(&second),
            Some(*gene),
            "reparse must observe the preserved gene_symbol for {input}"
        );

        // Identity: parse(input) == parse(parse(input).to_string()).
        assert_eq!(
            first, second,
            "structural equality must survive parse → display → parse for {input}"
        );
    }
}

// =============================================================================
// Hash + Eq round-trip identity: `Some(g)`, `Some(g')`, and `None` are all
// distinct under the derived Hash/Eq. The round-trip parse → display →
// parse preserves identity (covered above as part of the per-input
// asserts via `assert_eq!(first, second, ...)`); this test pins the
// negative half — that the round-trip does *not* collapse `Some(g)` and
// `None` to equal values.
// =============================================================================

#[test]
fn round_trip_identity_distinguishes_selector_presence() {
    use std::collections::HashSet;

    let with_gene_a = parse_hgvs("NM_000088.3(COL1A1):c.100A>G").unwrap();
    let with_gene_b = parse_hgvs("NM_000088.3(COL1A2):c.100A>G").unwrap();
    let bare = parse_hgvs("NM_000088.3:c.100A>G").unwrap();

    // Round-trip each through Display + reparse and verify identity holds.
    let with_gene_a_rt = parse_hgvs(&with_gene_a.to_string()).unwrap();
    let with_gene_b_rt = parse_hgvs(&with_gene_b.to_string()).unwrap();
    let bare_rt = parse_hgvs(&bare.to_string()).unwrap();

    assert_eq!(with_gene_a, with_gene_a_rt);
    assert_eq!(with_gene_b, with_gene_b_rt);
    assert_eq!(bare, bare_rt);

    // After round-trip the three remain distinct.
    assert_ne!(with_gene_a_rt, bare_rt);
    assert_ne!(with_gene_a_rt, with_gene_b_rt);
    assert_ne!(with_gene_b_rt, bare_rt);

    let mut set: HashSet<HgvsVariant> = HashSet::new();
    set.insert(with_gene_a_rt);
    set.insert(with_gene_b_rt);
    set.insert(bare_rt);
    assert_eq!(
        set.len(),
        3,
        "round-tripped variants with distinct selector presence must hash distinctly"
    );
}
