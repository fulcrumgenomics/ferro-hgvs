//! Display preserves the gene-symbol selector on round-trip (#121).
//!
//! Per HGVS Nomenclature, the gene-symbol selector (e.g. `(COL1A1)` in
//! `NM_000088.3(COL1A1):c.459A>G`) is informational disambiguation. ferro's
//! parsing layer captures it into `variant.<X>.gene_symbol: Option<String>`
//! and normalization propagates it. Until #121, every concrete `Display`
//! impl except `RnaFusionBreakpoint` dropped the selector.
//!
//! This file pins the principled round-trip policy:
//!
//! 1. **Preserve when present in input** — Display emits `accession(gene):`
//!    whenever `gene_symbol` is `Some(g)`, across all coordinate types
//!    (`g`/`c`/`n`/`r`/`p`/`m`/`o`).
//! 2. **Do not synthesize when absent** — Display emits the bare
//!    `accession:` form when `gene_symbol` is `None`, even if the underlying
//!    transcript provider could supply one.
//! 3. **Idempotent re-parse** — `parse → display → parse` is byte-stable
//!    and the second parse sees `gene_symbol = Some(<gene>)`.
//!
//! These cases are the inverse of the strip-pinned cases in PR #135's
//! `tests/gene_selector_roundtrip.rs`; the `// TODO #121` markers there
//! become preserve-assertions as a small follow-up after this PR lands.
//!
//! Spec rationale: <https://hgvs-nomenclature.org>. Mutalyzer and
//! VariantValidator both preserve the selector on round-trip; ferro
//! aligns with that policy here.

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
        HgvsVariant::RnaFusion(_)
        | HgvsVariant::Allele(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// Assert that `input` round-trips through `parse → display` byte-stably and
/// that the selector survives a second parse.
fn assert_display_preserves_gene_selector(input: &str, expected_gene: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    assert_eq!(
        gene_symbol_of(&v),
        Some(expected_gene),
        "parse must capture gene_symbol={expected_gene:?} for {input:?}",
    );

    let displayed = v.to_string();
    assert_eq!(
        displayed, input,
        "Display must preserve the gene-symbol selector verbatim",
    );

    // Re-parsing the displayed form must succeed, remain stable, and the
    // selector must still be present (not synthesized — actually preserved).
    let reparsed =
        parse_hgvs(&displayed).unwrap_or_else(|e| panic!("reparse {displayed:?} failed: {e:?}"));
    assert_eq!(reparsed.to_string(), displayed, "Display is not idempotent");
    assert_eq!(
        gene_symbol_of(&reparsed),
        Some(expected_gene),
        "reparse must observe the preserved gene_symbol",
    );
}

// =============================================================================
// Preserve: every supported accession family emits `(gene)` after the accession.
// =============================================================================

#[test]
fn display_preserves_gene_selector_refseq_nm_cds() {
    // Canonical example from the HGVS Nomenclature page.
    assert_display_preserves_gene_selector("NM_000088.3(COL1A1):c.459A>G", "COL1A1");
}

#[test]
fn display_preserves_gene_selector_refseq_nr_noncoding() {
    // NR_ is a non-coding RNA RefSeq accession; coordinate type is `n.`.
    assert_display_preserves_gene_selector("NR_046018.2(DDX11L1):n.100A>G", "DDX11L1");
}

#[test]
fn display_drops_gene_selector_refseq_np_protein() {
    // Per HGVS syntax.yaml 119–128 the protein-variant grammar has no
    // parenthesized selector position; the `accession(selector):c.` form
    // (syntax.yaml 213) is reserved for genomic-to-transcript projection on
    // c./n., not p. The parser remains permissive and `gene_symbol` round-trips
    // on the struct, but Display emits the spec-compliant `NP_*:p.` form.
    // See #310.
    use ferro_hgvs::hgvs::parser::parse_hgvs;
    let v = parse_hgvs("NP_000079.2(COL1A1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("COL1A1"));
    assert_eq!(v.to_string(), "NP_000079.2:p.(Arg8Gln)");
}

#[test]
fn display_preserves_gene_selector_genome() {
    // `(FLT3)` here is a gene-symbol selector, not a compound-ref wrapper.
    assert_display_preserves_gene_selector("NC_000013.11(FLT3):g.12345A>G", "FLT3");
}

#[test]
fn display_preserves_gene_selector_ensembl() {
    assert_display_preserves_gene_selector("ENST00000380152.7(BRCA2):c.100A>G", "BRCA2");
}

#[test]
fn display_preserves_gene_selector_lrg() {
    // LRG references come in three flavors: bare (`LRG_<n>`), transcript
    // (`LRG_<n>t<m>`), and protein (`LRG_<n>p<m>`). The bare and transcript
    // forms preserve the gene selector verbatim on `c.`/`n.`. The protein
    // form does NOT round-trip the selector through Display: the HGVS
    // protein-variant grammar (syntax.yaml 119–128) has no parenthesized
    // selector position, so Display emits the spec-compliant `LRG_*p*:p.`
    // form. See #310.
    assert_display_preserves_gene_selector("LRG_199(BRCA1):c.100A>G", "BRCA1");
    assert_display_preserves_gene_selector("LRG_199t1(BRCA1):c.100A>G", "BRCA1");

    use ferro_hgvs::hgvs::parser::parse_hgvs;
    let v = parse_hgvs("LRG_199p1(BRCA1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("BRCA1"));
    assert_eq!(v.to_string(), "LRG_199p1:p.(Arg8Gln)");
}

#[test]
fn display_preserves_gene_selector_simple_accession() {
    // Custom (non-RefSeq, non-Ensembl, non-LRG) accession with a selector.
    assert_display_preserves_gene_selector("MYREF_SEQ(GENE1):c.100A>G", "GENE1");
}

// =============================================================================
// Do-not-synthesize: bare input stays bare on Display.
// =============================================================================

/// Assert Display emits the bare form when the input had no selector and
/// re-parsing produces a variant with `gene_symbol == None`.
fn assert_display_does_not_synthesize(input: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    assert_eq!(
        gene_symbol_of(&v),
        None,
        "parse must not synthesize gene_symbol for bare input {input:?}",
    );
    assert_eq!(
        v.to_string(),
        input,
        "Display must not synthesize a gene selector when input had none",
    );
}

#[test]
fn display_does_not_synthesize_refseq_nm_cds() {
    assert_display_does_not_synthesize("NM_000088.3:c.459A>G");
}

#[test]
fn display_does_not_synthesize_refseq_nr_noncoding() {
    assert_display_does_not_synthesize("NR_046018.2:n.100A>G");
}

#[test]
fn display_does_not_synthesize_refseq_np_protein() {
    assert_display_does_not_synthesize("NP_000079.2:p.(Arg8Gln)");
}

#[test]
fn display_does_not_synthesize_genome() {
    assert_display_does_not_synthesize("NC_000013.11:g.12345A>G");
}

#[test]
fn display_does_not_synthesize_ensembl() {
    assert_display_does_not_synthesize("ENST00000380152.7:c.100A>G");
}

// =============================================================================
// Compound-ref + selector: a bare gene selector composes with the existing
// `outer(inner)` compound-ref wrapper.
// =============================================================================

#[test]
fn display_preserves_gene_selector_with_compound_ref() {
    // `NC_…(NM_…)` is a compound-ref accession (genomic_context); the
    // trailing `(GENE)` is the gene-symbol selector. Both must round-trip.
    assert_display_preserves_gene_selector("NC_000013.11(NM_004119.3)(FLT3):c.100A>G", "FLT3");
}

#[test]
fn display_compound_ref_without_selector_unchanged() {
    // Compound-ref alone (no gene selector) must not be confused with a
    // selector. The display is unchanged from PR #70's behavior.
    let input = "NC_000013.11(NM_004119.3):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), None);
    assert_eq!(v.to_string(), input);
}

// =============================================================================
// Allele compact form: `ACC(GENE):c.[edit1;edit2]` carries the selector once.
// =============================================================================

#[test]
fn display_preserves_gene_selector_in_compact_allele() {
    // Compact cis allele form: when both sub-variants share an accession
    // and gene selector, the prefix emits `acc(gene):type.[…]` once.
    let input = "NM_000088.3(COL1A1):c.[100A>G;200C>T]";
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    assert_eq!(v.to_string(), input);

    // Idempotent re-parse.
    let reparsed = parse_hgvs(&v.to_string()).unwrap();
    assert_eq!(reparsed.to_string(), input);
}

// =============================================================================
// Normalize → Display: the selector survives both stages.
// =============================================================================

#[test]
fn normalize_then_display_preserves_gene_selector() {
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let input = "NM_000088.3(COL1A1):c.4A>G";
    let variant = parse_hgvs(input).unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    // gene_symbol survives normalization, and Display preserves it on output.
    assert_eq!(gene_symbol_of(&normalized), Some("COL1A1"));
    assert!(
        normalized.to_string().contains("(COL1A1)"),
        "normalized Display must preserve the gene selector: {}",
        normalized
    );
}

#[test]
fn normalize_does_not_synthesize_gene_selector_in_display() {
    // When the input lacked a selector, normalization must not invent one,
    // even though the underlying transcript record knows the gene symbol.
    // The byte-equality assertion pins the exact canonical bare form so any
    // future change that adds parens for an unrelated reason cannot
    // silently introduce a synthesized selector.
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let input = "NM_000088.3:c.4A>G";
    let variant = parse_hgvs(input).unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    assert_eq!(gene_symbol_of(&normalized), None);
    assert_eq!(
        normalized.to_string(),
        "NM_000088.3:c.4A>G",
        "Display must not synthesize a gene selector"
    );
}

// =============================================================================
// Per-coordinate-type Preserve coverage for the kinds not exercised in the top
// section (the top section already covers `c`, `n`, `p`, `g`).
// =============================================================================

#[test]
fn display_preserves_gene_selector_rna() {
    // `r.` form (RNA): selector follows the same rule as DNA forms.
    assert_display_preserves_gene_selector("NM_000088.3(COL1A1):r.100a>g", "COL1A1");
}

#[test]
fn display_preserves_gene_selector_mitochondrial() {
    // `m.` form: spec explicitly endorses `NC_012920.1(MT-…):m.…` — see
    // assets/hgvs-nomenclature/docs/background/refseq.md lines 197–203.
    assert_display_preserves_gene_selector("NC_012920.1(MT-ND1):m.3460G>A", "MT-ND1");
}

#[test]
fn display_preserves_gene_selector_circular() {
    // `o.` form (circular DNA, SVD-WG006). ferro accepts the selector
    // verbatim; Display preserves it.
    assert_display_preserves_gene_selector("NC_001416.1(genE):o.1A>G", "genE");
}

// =============================================================================
// Trans-allele expanded form: each sub-variant carries its own selector through
// the `[v1];[v2]` wrapping (no shared prefix; no information loss).
// =============================================================================

#[test]
fn display_preserves_gene_selector_in_trans_allele_same_acc_same_gene() {
    // Same accession + same gene + trans phase. The compact-trans form
    // `ACC(GENE):c.[a];[b]` is legal here, and ferro emits it because both
    // sub-variants agree on accession and gene_symbol. The `Allele` wrapper
    // returns `None` from `gene_symbol()` (selectors live on sub-variants),
    // so we assert via byte-stable round-trip + per-sub-variant inspection
    // rather than the generic helper.
    let input = "NM_000088.3(COL1A1):c.[100A>G];[200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(
        v.to_string(),
        input,
        "Display must preserve compact-trans selector"
    );

    if let HgvsVariant::Allele(a) = &v {
        for sub in &a.variants {
            assert_eq!(gene_symbol_of(sub), Some("COL1A1"));
        }
    } else {
        panic!("expected Allele; got {:?}", v);
    }

    let reparsed = parse_hgvs(&v.to_string()).unwrap();
    assert_eq!(reparsed.to_string(), input, "reparse must be byte-stable");
}

#[test]
fn display_preserves_gene_selector_in_trans_allele_expanded_form() {
    // Parse the expanded `[v1];[v2]` notation; each sub-variant captures the
    // gene_symbol independently. Because both sub-variants share accession and
    // gene_symbol, `all_share_accession_and_type` permits the compact-trans
    // output `ACC(GENE):c.[a];[b]`, so this test does NOT exercise the
    // per-kind expanded-form Display path (that case lives in
    // `display_preserves_independent_gene_selectors_in_mixed_acc_trans`, where
    // mismatched accessions force the expanded path). What this test pins is
    // sub-variant `gene_symbol` capture under the expanded-input → compact-
    // trans-output normalization, plus byte-stable re-parse of the output.
    let input = "[NM_000088.3(COL1A1):c.100A>G];[NM_000088.3(COL1A1):c.200C>T]";
    let v = parse_hgvs(input).unwrap();
    let s = v.to_string();
    let reparsed = parse_hgvs(&s).unwrap();
    assert_eq!(reparsed.to_string(), s, "reparse must be byte-stable");

    // Both inner sub-variants must observe the gene_symbol.
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
    // Different accessions, each with its own gene symbol. Per #121 the
    // expanded trans form must preserve both selectors, one per sub-variant.
    // This is the round-trip-stability case for cross-gene compound
    // heterozygosity reports.
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
    assert_eq!(reparsed.to_string(), input, "reparse must be byte-stable");
}

// =============================================================================
// Allele compact forms (cis, trans, unknown-phase) all carry the selector once
// on the prefix via `write_compact_prefix`.
// =============================================================================

#[test]
fn display_preserves_gene_selector_in_compact_trans_allele() {
    // Compact-trans form: `ACC(GENE):c.[edit1];[edit2]`. The selector lives
    // on the prefix once (via `write_compact_prefix`); each bracketed group
    // is a single trans variant.
    let input = "NM_000088.3(COL1A1):c.[100A>G];[200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), input);
    let reparsed = parse_hgvs(&v.to_string()).unwrap();
    assert_eq!(reparsed.to_string(), input);
}

#[test]
fn display_preserves_gene_selector_in_unknown_phase_compact_allele() {
    // Unknown-phase compact form: `ACC(GENE):c.edit1(;)edit2`. No surrounding
    // brackets; the selector still lives on the prefix once.
    let input = "NM_000088.3(COL1A1):c.100A>G(;)200C>T";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), input);
    let reparsed = parse_hgvs(&v.to_string()).unwrap();
    assert_eq!(reparsed.to_string(), input);
}

// =============================================================================
// Hash + Eq: gene_symbol participates in the derived Hash/Eq on every variant
// struct. Two variants with the same accession + edit but different selectors
// must compare unequal and hash distinctly so set/map deduplication respects
// the user-supplied disambiguation.
// =============================================================================

#[test]
fn variants_with_different_gene_symbols_are_distinct() {
    use std::collections::HashSet;

    let with_gene = parse_hgvs("NM_000088.3(COL1A1):c.100A>G").unwrap();
    let bare = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
    let other_gene = parse_hgvs("NM_000088.3(COL1A2):c.100A>G").unwrap();
    let with_gene_dup = parse_hgvs("NM_000088.3(COL1A1):c.100A>G").unwrap();

    assert_ne!(
        with_gene, bare,
        "Some(gene) != None on otherwise-equal variants"
    );
    assert_ne!(
        with_gene, other_gene,
        "different gene symbols must compare unequal"
    );
    assert_eq!(
        with_gene, with_gene_dup,
        "byte-equal inputs must compare equal"
    );

    let mut set = HashSet::new();
    set.insert(with_gene.clone());
    set.insert(bare);
    set.insert(other_gene);
    set.insert(with_gene_dup);
    assert_eq!(
        set.len(),
        3,
        "set should contain three distinct variants (the COL1A1 dup collapses)"
    );
}

// =============================================================================
// `HgvsVariant::gene_symbol()` accessor contract.
//
// The accessor is documented as returning `None` for variant kinds that don't
// carry a single top-level selector (`Allele`, `RnaFusion`, `NullAllele`,
// `UnknownAllele`). Selectors on those nested kinds are emitted by their own
// `Display` impls (e.g. `RnaFusionBreakpoint::Display`). This test pins that
// contract so a future "surface nested gene_symbol" change is a deliberate
// API decision rather than an accident.
// =============================================================================

#[test]
fn gene_symbol_accessor_returns_none_for_nested_kinds() {
    // Allele: selector lives on each sub-variant, not on the AlleleVariant.
    let allele = parse_hgvs("NM_000088.3(COL1A1):c.[100A>G;200C>T]").unwrap();
    assert!(matches!(allele, HgvsVariant::Allele(_)));
    assert_eq!(
        allele.gene_symbol(),
        None,
        "Allele's top-level gene_symbol() is None; selectors live on sub-variants"
    );

    // RnaFusion: selector lives on each breakpoint, not on the fusion.
    let fusion =
        parse_hgvs("NM_152263.2(ACTN1):r.-115_775::NM_002609.3(BCR):r.1580_*1924").unwrap();
    assert!(matches!(fusion, HgvsVariant::RnaFusion(_)));
    assert_eq!(
        fusion.gene_symbol(),
        None,
        "RnaFusion's top-level gene_symbol() is None; selectors live on breakpoints"
    );
    // The breakpoint selectors must still appear in the Display output.
    let s = fusion.to_string();
    assert!(
        s.contains("(ACTN1)"),
        "5' breakpoint gene must round-trip: {}",
        s
    );
    assert!(
        s.contains("(BCR)"),
        "3' breakpoint gene must round-trip: {}",
        s
    );
}

// =============================================================================
// Edge cases: empty selector, whitespace padding, casing, special characters.
// ferro's policy is "preserve verbatim" — the parser does not sanitize
// gene_symbol, and Display emits whatever was captured. The HGVS spec's
// canonical form uses approved HGNC symbols (no whitespace, no quotes), but
// ferro stays bytes-in / bytes-out so callers can detect malformed input.
// =============================================================================

#[test]
fn empty_selector_is_treated_as_no_selector() {
    // `NM_X():c.…` parses with gene_symbol = None (per the parser's
    // "Treat empty strings as None" rule). Display must emit the bare form.
    let v = parse_hgvs("NM_000088.3():c.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), None);
    assert_eq!(
        v.to_string(),
        "NM_000088.3:c.100A>G",
        "empty selector must drop to bare form, not synthesize an empty `()`"
    );
}

#[test]
fn whitespace_padded_selector_round_trips_verbatim() {
    // Padding is preserved in `gene_symbol` and round-trips on Display.
    // Spec uses HGNC symbols which contain no whitespace; this pin documents
    // ferro's bytes-in/bytes-out policy.
    let input = "NM_000088.3( COL1A1 ):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), Some(" COL1A1 "));
    assert_eq!(v.to_string(), input);
}

#[test]
fn lowercase_selector_preserved_verbatim() {
    // Spec uses uppercase HGNC symbols; ferro does not normalize case.
    let input = "NM_000088.3(col1a1):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), Some("col1a1"));
    assert_eq!(v.to_string(), input);
}

#[test]
fn hyphenated_selector_round_trips() {
    // Mitochondrial genes carry hyphens (MT-ND1, MT-TL1). The selector
    // parser uses `take_while(c != ')')`, so any non-paren bytes are
    // accepted verbatim. The spec endorses this exact form (see
    // assets/hgvs-nomenclature/docs/background/refseq.md lines 197–203).
    assert_display_preserves_gene_selector("NC_012920.1(MT-TL1):m.3243A>G", "MT-TL1");
    assert_display_preserves_gene_selector("NC_012920.1(MT-ND1):m.3460G>A", "MT-ND1");
    assert_display_preserves_gene_selector("NC_012920.1(MT-TL1):n.14A>G", "MT-TL1");
}
