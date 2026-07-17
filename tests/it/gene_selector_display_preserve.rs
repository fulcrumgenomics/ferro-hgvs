//! Display's reference-type policy for the gene-symbol selector (#1051).
//!
//! Per HGVS Nomenclature (`refseq.md`), the parenthetical after a reference is
//! a *specification* whose accepted values are transcript/protein accessions,
//! **not** gene symbols; and a transcript reference needs no specification at
//! all (bare `NM_…:c.…` is canonical). ferro's `Display` therefore splits by
//! reference type — the parser still captures `gene_symbol` on every family,
//! but the formatter:
//!
//! 1. **Drops** the selector on a **transcript** reference
//!    (`NM`/`NR`/`XM`/`XR`/`ENST`/`LRG_<n>t<m>`), emitting the canonical bare
//!    form — the same treatment `p.` has always had (#310).
//! 2. **Preserves** it on a **non-transcript** reference: genomic
//!    (`NC`/`NG`/bare `LRG_<n>`) where the gene stands in for a transcript, the
//!    spec-sanctioned mitochondrial exception, and unclassifiable custom
//!    references.
//! 3. **Never synthesizes** it — a bare input stays bare, and normalization
//!    does not invent a selector from transcript metadata.
//!
//! The value always remains on the struct for programmatic access regardless
//! of whether Display emits it (`variant.<X>.gene_symbol`). This file targets
//! the *formatter* per-kind; end-to-end round-trip identity lives in
//! `tests/gene_selector_roundtrip.rs`.

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
        HgvsVariant::RnaFusion(_)
        | HgvsVariant::Allele(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// Transcript reference: the parser captures `gene_symbol = Some(gene)`, but
/// Display drops it, emitting `bare` (input with the `(GENE)` selector removed).
fn assert_display_drops_gene_selector(input: &str, gene: &str, bare: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    assert_eq!(
        gene_symbol_of(&v),
        Some(gene),
        "parse must still capture gene_symbol={gene:?} for {input:?}",
    );
    assert_eq!(
        v.to_string(),
        bare,
        "Display must drop the selector on a transcript reference (#1051)",
    );
    // The bare form is idempotent and carries no selector on reparse.
    let reparsed = parse_hgvs(bare).unwrap_or_else(|e| panic!("reparse {bare:?} failed: {e:?}"));
    assert_eq!(reparsed.to_string(), bare, "bare form must be idempotent");
    assert_eq!(gene_symbol_of(&reparsed), None);
}

/// Non-transcript reference: Display preserves the selector verbatim and the
/// round-trip is byte-stable, with the reparse re-observing it.
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
        "Display must preserve the selector on a non-transcript reference",
    );

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
// Drop: a transcript reference emits the canonical bare form.
// =============================================================================

#[test]
fn display_drops_gene_selector_refseq_nm_cds() {
    // Canonical example from the HGVS Nomenclature page.
    assert_display_drops_gene_selector(
        "NM_000088.3(COL1A1):c.459A>G",
        "COL1A1",
        "NM_000088.3:c.459A>G",
    );
}

#[test]
fn display_drops_gene_selector_refseq_nr_noncoding() {
    // NR_ is a non-coding RNA RefSeq (transcript) accession; coordinate `n.`.
    assert_display_drops_gene_selector(
        "NR_046018.2(DDX11L1):n.100A>G",
        "DDX11L1",
        "NR_046018.2:n.100A>G",
    );
}

#[test]
fn display_drops_gene_selector_refseq_np_protein() {
    // Protein `p.` has always dropped the selector (#310); #1051 makes the
    // transcript DNA/RNA forms consistent with it.
    let v = parse_hgvs("NP_000079.2(COL1A1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("COL1A1"));
    assert_eq!(v.to_string(), "NP_000079.2:p.(Arg8Gln)");
}

#[test]
fn display_drops_gene_selector_ensembl() {
    assert_display_drops_gene_selector(
        "ENST00000380152.7(BRCA2):c.100A>G",
        "BRCA2",
        "ENST00000380152.7:c.100A>G",
    );
}

#[test]
fn display_drops_gene_selector_rna() {
    // `r.` form on a transcript reference: dropped like `c.`/`n.`.
    assert_display_drops_gene_selector(
        "NM_000088.3(COL1A1):r.100a>g",
        "COL1A1",
        "NM_000088.3:r.100a>g",
    );
}

#[test]
fn display_drops_gene_selector_lrg_transcript() {
    // An LRG transcript (`LRG_<n>t<m>`) is a transcript reference → drop.
    assert_display_drops_gene_selector("LRG_199t1(BRCA1):c.100A>G", "BRCA1", "LRG_199t1:c.100A>G");

    // The LRG protein form drops it via the protein grammar (#310).
    let v = parse_hgvs("LRG_199p1(BRCA1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("BRCA1"));
    assert_eq!(v.to_string(), "LRG_199p1:p.(Arg8Gln)");
}

// =============================================================================
// Preserve: a non-transcript reference keeps the selector.
// =============================================================================

#[test]
fn display_preserves_gene_selector_genome() {
    // `(FLT3)` is a gene-symbol selector on a single-segment genomic reference
    // (no transcript named), so ferro preserves it.
    assert_display_preserves_gene_selector("NC_000013.11(FLT3):g.12345A>G", "FLT3");
}

#[test]
fn display_preserves_gene_selector_mitochondrial() {
    // `m.` form: the spec explicitly endorses `NC_012920.1(MT-…):m.…` — the
    // gene has no transcript reference (refseq.md lines 197–203), so the
    // selector is the only way to specify it.
    assert_display_preserves_gene_selector("NC_012920.1(MT-ND1):m.3460G>A", "MT-ND1");
}

#[test]
fn display_preserves_gene_selector_circular() {
    // `o.` form (circular DNA, SVD-WG006) on a genomic reference.
    assert_display_preserves_gene_selector("NC_001416.1(genE):o.1A>G", "genE");
}

#[test]
fn display_preserves_gene_selector_lrg_bare_genome() {
    // Bare LRG (`LRG_<n>`) is the genomic record — a non-transcript reference,
    // so the selector is preserved.
    assert_display_preserves_gene_selector("LRG_199(BRCA1):g.100A>G", "BRCA1");
}

#[test]
fn display_preserves_gene_selector_simple_accession() {
    // Custom (non-RefSeq, non-Ensembl, non-LRG) accession: unclassifiable, so
    // ferro conservatively preserves the selector rather than drop information.
    assert_display_preserves_gene_selector("MYREF_SEQ(GENE1):c.100A>G", "GENE1");
}

// =============================================================================
// Do-not-synthesize: bare input stays bare on Display.
// =============================================================================

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
// Compound-ref + selector: `NC_…(NM_…)(GENE):c.…`. The compound ref already
// names the transcript (the spec-correct specification), so `self` is the
// transcript `NM_` and the redundant `(GENE)` is dropped, leaving the
// canonical `NC_(NM_):c.…` (#1051).
// =============================================================================

#[test]
fn display_drops_redundant_gene_selector_with_compound_ref() {
    let v = parse_hgvs("NC_000013.11(NM_004119.3)(FLT3):c.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("FLT3"));
    assert_eq!(v.to_string(), "NC_000013.11(NM_004119.3):c.100A>G");
}

#[test]
fn display_compound_ref_without_selector_unchanged() {
    // Compound-ref alone (no gene selector) is unchanged.
    let input = "NC_000013.11(NM_004119.3):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), None);
    assert_eq!(v.to_string(), input);
}

// =============================================================================
// Allele forms on a transcript reference: the selector is dropped from the
// compact prefix / expanded members, but the compaction decision still keys
// off `gene_symbol` equality.
// =============================================================================

#[test]
fn display_drops_gene_selector_in_compact_allele() {
    let input = "NM_000088.3(COL1A1):c.[100A>G;200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), "NM_000088.3:c.[100A>G;200C>T]");
    let reparsed = parse_hgvs(&v.to_string()).unwrap();
    assert_eq!(reparsed.to_string(), "NM_000088.3:c.[100A>G;200C>T]");
}

#[test]
fn display_drops_gene_selector_in_compact_trans_allele() {
    let input = "NM_000088.3(COL1A1):c.[100A>G];[200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), "NM_000088.3:c.[100A>G];[200C>T]");
}

#[test]
fn display_drops_gene_selector_in_unknown_phase_compact_allele() {
    let input = "NM_000088.3(COL1A1):c.100A>G(;)200C>T";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), "NM_000088.3:c.100A>G(;)200C>T");
}

#[test]
fn display_drops_gene_selector_in_trans_allele_same_acc_same_gene() {
    // Same accession + same gene + trans phase → compact-trans, selector
    // dropped from the transcript-reference prefix. The field survives on each
    // sub-variant.
    let input = "NM_000088.3(COL1A1):c.[100A>G];[200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), "NM_000088.3:c.[100A>G];[200C>T]");

    if let HgvsVariant::Allele(a) = &v {
        for sub in &a.variants {
            assert_eq!(gene_symbol_of(sub), Some("COL1A1"));
        }
    } else {
        panic!("expected Allele; got {:?}", v);
    }
}

#[test]
fn display_drops_independent_gene_selectors_in_mixed_acc_trans() {
    // Different transcript accessions with different gene symbols. Both
    // selectors are dropped from Display; the fields survive on the members.
    let input = "[NM_000088.3(COL1A1):c.100A>G];[NM_000089.3(COL1A2):c.200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(
        v.to_string(),
        "[NM_000088.3:c.100A>G];[NM_000089.3:c.200C>T]"
    );

    if let HgvsVariant::Allele(a) = &v {
        assert_eq!(a.variants.len(), 2);
        assert_eq!(gene_symbol_of(&a.variants[0]), Some("COL1A1"));
        assert_eq!(gene_symbol_of(&a.variants[1]), Some("COL1A2"));
    } else {
        panic!("expected Allele; got {:?}", v);
    }
}

// =============================================================================
// Normalize → Display: the selector survives normalization on the struct, and
// Display applies the reference-type policy afterward.
// =============================================================================

#[test]
fn normalize_then_display_drops_gene_selector_on_transcript_ref() {
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let input = "NM_000088.3(COL1A1):c.4A>G";
    let variant = parse_hgvs(input).unwrap();
    let normalized = normalizer.normalize(&variant).unwrap();

    // gene_symbol survives normalization on the struct...
    assert_eq!(gene_symbol_of(&normalized), Some("COL1A1"));
    // ...but Display drops it on the transcript reference (#1051). Pin the exact
    // bare form, not just the absence of the selector, so malformed output can't
    // slip past (a substitution does not shift under normalization).
    assert_eq!(
        normalized.to_string(),
        "NM_000088.3:c.4A>G",
        "Display on a transcript reference must emit the canonical bare form"
    );
}

#[test]
fn normalize_does_not_synthesize_gene_selector_in_display() {
    // When the input lacked a selector, normalization must not invent one,
    // even though the underlying transcript record knows the gene symbol.
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
// Hash + Eq: gene_symbol participates in the derived Hash/Eq on every variant
// struct (parser-level; independent of the Display drop). Two variants with the
// same accession + edit but different selectors must compare unequal.
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
// `HgvsVariant::gene_symbol()` accessor contract: `None` for nested kinds
// (`Allele`, `RnaFusion`, …); selectors on those live on the sub-parts.
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
    // Both breakpoints are transcript references (`NM_`), so Display drops the
    // selectors (#1051) — the breakpoint genes must NOT appear in the output.
    let s = fusion.to_string();
    assert_eq!(
        s, "NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924",
        "transcript-reference fusion breakpoints drop their gene selectors: {s}",
    );
    // The captured gene_symbol still lives on each breakpoint struct.
    if let HgvsVariant::RnaFusion(fx) = &fusion {
        assert_eq!(fx.five_prime.gene_symbol.as_deref(), Some("ACTN1"));
        assert_eq!(fx.three_prime.gene_symbol.as_deref(), Some("BCR"));
    } else {
        unreachable!()
    }
}

// =============================================================================
// Edge cases: empty selector, whitespace/case (captured verbatim on parse,
// then subject to the reference-type Display policy).
// =============================================================================

#[test]
fn empty_selector_is_treated_as_no_selector() {
    // `NM_X():c.…` parses with gene_symbol = None; Display emits the bare form.
    let v = parse_hgvs("NM_000088.3():c.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), None);
    assert_eq!(
        v.to_string(),
        "NM_000088.3:c.100A>G",
        "empty selector must drop to bare form, not synthesize an empty `()`"
    );
}

#[test]
fn whitespace_padded_selector_captured_then_dropped() {
    // Padding is preserved in `gene_symbol` on parse; a transcript reference
    // then drops it on Display (#1051).
    let input = "NM_000088.3( COL1A1 ):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), Some(" COL1A1 "));
    assert_eq!(v.to_string(), "NM_000088.3:c.100A>G");
}

#[test]
fn lowercase_selector_captured_then_dropped() {
    // ferro does not normalize case; the transcript reference drops it anyway.
    let input = "NM_000088.3(col1a1):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), Some("col1a1"));
    assert_eq!(v.to_string(), "NM_000088.3:c.100A>G");
}

#[test]
fn hyphenated_selector_preserved_on_mitochondrial_reference() {
    // Mitochondrial genes carry hyphens (MT-ND1, MT-TL1) and have no transcript
    // reference, so the selector is the spec-sanctioned exception and preserved
    // (refseq.md lines 197–203).
    assert_display_preserves_gene_selector("NC_012920.1(MT-TL1):m.3243A>G", "MT-TL1");
    assert_display_preserves_gene_selector("NC_012920.1(MT-ND1):m.3460G>A", "MT-ND1");
    assert_display_preserves_gene_selector("NC_012920.1(MT-TL1):n.14A>G", "MT-TL1");
}
