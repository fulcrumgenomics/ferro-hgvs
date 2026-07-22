//! Round-trip tests for the optional gene-symbol selector.
//!
//! HGVS Nomenclature allows an informational gene-symbol selector in
//! parentheses after the accession to disambiguate when the same accession
//! spans multiple genes (e.g. `NM_000088.3(COL1A1):c.459A>G`).
//!
//! PR #70 (`b34ce35`) extended the parser so the selector is accepted on
//! non-RefSeq accessions and captured into the variant struct's
//! `gene_symbol: Option<String>`. The parser is unchanged by #1051 and still
//! captures the selector on every accession family.
//!
//! Display policy (#1051)
//! ----------------------
//! Per HGVS Nomenclature (`refseq.md`), the parenthetical after a reference is
//! a *specification* whose accepted values are transcript/protein accessions,
//! **not** gene symbols; and a transcript reference is fully specified on its
//! own (bare `NM_…:c.…` is canonical) so it takes no selector at all. ferro
//! therefore splits by **reference type** on Display:
//!
//! - **Transcript reference** (`NM`/`NR`/`XM`/`XR`/`ENST`/`LRG_<n>t<m>`):
//!   Display *drops* the selector, emitting the canonical bare form — the
//!   same treatment `p.` has always had (#310). The value stays on the struct
//!   for programmatic access; parse → Display → reparse is therefore lossy for
//!   the selector (by design, mirroring a normalization).
//! - **Non-transcript reference** (genomic `NC`/`NG`/bare `LRG_<n>`,
//!   mitochondrial, custom): Display *preserves* the selector — it is either
//!   the spec-sanctioned exception (a gene with no transcript, e.g.
//!   mitochondrial) or a value ferro cannot safely drop. Round-trip is
//!   byte-stable and reparse re-observes the selector.
//!
//! The negative half of the policy — "do not synthesize when absent" — is
//! unchanged: a bare input stays bare and normalization never invents a
//! selector from transcript metadata.

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

/// Assert the #1051 policy for a **transcript reference**: the parser captures
/// `gene_symbol = Some(gene)`, but Display drops it, emitting `bare`. The bare
/// form is idempotent and its reparse observes no selector (the Display leg is
/// intentionally lossy for the selector, like `p.`).
fn assert_transcript_ref_drops_selector(input: &str, gene: &str, bare: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    assert_eq!(
        gene_symbol_of(&v),
        Some(gene),
        "parse must still capture gene_symbol for {input:?}",
    );

    let displayed = v.to_string();
    assert_eq!(
        displayed, bare,
        "Display must drop the selector on a transcript reference (#1051)",
    );

    let reparsed =
        parse_hgvs(&displayed).unwrap_or_else(|e| panic!("reparse {displayed:?} failed: {e:?}"));
    assert_eq!(reparsed.to_string(), bare, "bare form must be idempotent");
    assert_eq!(
        gene_symbol_of(&reparsed),
        None,
        "reparse of the dropped form must observe no selector",
    );
}

/// Assert the #1051 policy for a **non-transcript reference** (genomic, mito,
/// custom): Display preserves the selector verbatim and the round-trip is
/// byte-stable, with the reparse re-observing it.
fn assert_nontranscript_ref_preserves_selector(input: &str) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    let expected_gene = gene_symbol_of(&v)
        .unwrap_or_else(|| panic!("test bug: input {input:?} must have a gene_symbol after parse"))
        .to_string();

    let displayed = v.to_string();
    assert_eq!(
        displayed, input,
        "Display must preserve the selector on a non-transcript reference (#1051)",
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
// Parse: gene_symbol is captured on every supported accession family. The
// parser is unchanged by #1051 — only Display differs by reference type.
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
// Display on a TRANSCRIPT reference: the selector is dropped to the canonical
// bare form (#1051), mirroring the long-standing `p.` behavior (#310).
// =============================================================================

#[test]
fn display_drops_gene_selector_refseq_nm() {
    assert_transcript_ref_drops_selector(
        "NM_000088.3(COL1A1):c.459A>G",
        "COL1A1",
        "NM_000088.3:c.459A>G",
    );
}

#[test]
fn display_drops_gene_selector_refseq_nr() {
    assert_transcript_ref_drops_selector(
        "NR_046018.2(DDX11L1):n.100A>G",
        "DDX11L1",
        "NR_046018.2:n.100A>G",
    );
}

#[test]
fn display_drops_gene_selector_refseq_np() {
    // Protein `p.` has always dropped the selector (#310); #1051 makes the
    // transcript DNA/RNA forms consistent with it.
    let v = parse_hgvs("NP_000079.2(COL1A1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("COL1A1"));
    assert_eq!(v.to_string(), "NP_000079.2:p.(Arg8Gln)");
}

#[test]
fn display_drops_gene_selector_ensembl() {
    assert_transcript_ref_drops_selector(
        "ENST00000380152.7(BRCA2):c.100A>G",
        "BRCA2",
        "ENST00000380152.7:c.100A>G",
    );
}

#[test]
fn display_drops_gene_selector_ensembl_species_qualified() {
    // #1057: a species-qualified Ensembl transcript (mouse ENSMUST, rat ENSRNOT,
    // …) is just as much a transcript reference as the human ENST, so the
    // gene-symbol selector must be dropped identically — not left dangling.
    assert_transcript_ref_drops_selector(
        "ENSMUST00000012345.1(Trp53):c.100A>G",
        "Trp53",
        "ENSMUST00000012345.1:c.100A>G",
    );
    assert_transcript_ref_drops_selector(
        "ENSRNOT00000012345.1(Tp53):c.100A>G",
        "Tp53",
        "ENSRNOT00000012345.1:c.100A>G",
    );
}

#[test]
fn display_drops_gene_selector_rna() {
    // `r.` form on a transcript reference: selector dropped like `c.`/`n.`.
    assert_transcript_ref_drops_selector(
        "NM_000088.3(COL1A1):r.100a>g",
        "COL1A1",
        "NM_000088.3:r.100a>g",
    );
}

// =============================================================================
// Display on a NON-transcript reference: the selector is preserved. Genomic
// references (`NC`/`NG`/bare `LRG`) carry the gene as an exception/legacy
// selector; the mitochondrial form is the spec-sanctioned exception; a custom
// accession cannot be classified so ferro keeps the value.
// =============================================================================

#[test]
fn display_preserves_gene_selector_genome() {
    // `(FLT3)` here is a gene-symbol selector on a single-segment genomic
    // reference (no transcript named), so ferro preserves it.
    assert_nontranscript_ref_preserves_selector("NC_000013.11(FLT3):g.12345A>G");
}

#[test]
fn display_preserves_gene_selector_mitochondrial() {
    // `m.` form: spec explicitly endorses `NC_012920.1(MT-…):m.…` — the gene
    // has no transcript reference (refseq.md lines 197-203), so the selector
    // is the *only* way to specify it and must be preserved.
    assert_nontranscript_ref_preserves_selector("NC_012920.1(MT-ND1):m.3460G>A");
}

#[test]
fn display_preserves_gene_selector_circular() {
    // `o.` form (circular DNA, SVD-WG006) on a genomic reference.
    assert_nontranscript_ref_preserves_selector("NC_001416.1(GENE):o.1A>G");
}

#[test]
fn display_preserves_gene_selector_simple_accession() {
    // Custom (unclassifiable) accession: ferro cannot prove it is a transcript,
    // so it conservatively preserves the selector rather than drop information.
    assert_nontranscript_ref_preserves_selector("MYREF_SEQ(GENE1):c.100A>G");
}

#[test]
fn display_no_selector_unchanged() {
    // Sanity check: Display is idempotent when the selector is absent and
    // is not synthesized from somewhere else (e.g. a transcript provider
    // lookup). This is the negative half of the policy.
    let input = "NM_000088.3:c.459A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), input);
    assert_eq!(gene_symbol_of(&v), None);
}

// =============================================================================
// LRG flavor matrix: bare (`LRG_<n>` -> g., NON-transcript, preserves),
// transcript (`LRG_<n>t<m>` -> c./n./r., transcript, DROPS), protein
// (`LRG_<n>p<m>` -> p., drops via the protein grammar).
// =============================================================================

#[test]
fn display_preserves_gene_selector_lrg_bare_genome() {
    // Bare LRG is the genomic record itself — a non-transcript reference, so
    // the selector is preserved (it stands in for a transcript that is not
    // named here).
    assert_nontranscript_ref_preserves_selector("LRG_199(BRCA1):g.100A>G");
}

#[test]
fn display_drops_gene_selector_lrg_transcript_cds() {
    assert_transcript_ref_drops_selector(
        "LRG_199t1(BRCA1):c.100A>G",
        "BRCA1",
        "LRG_199t1:c.100A>G",
    );
}

#[test]
fn display_drops_gene_selector_lrg_transcript_noncoding() {
    assert_transcript_ref_drops_selector(
        "LRG_199t1(BRCA1):n.100A>G",
        "BRCA1",
        "LRG_199t1:n.100A>G",
    );
}

#[test]
fn display_drops_gene_selector_lrg_transcript_rna() {
    assert_transcript_ref_drops_selector(
        "LRG_199t1(BRCA1):r.100a>g",
        "BRCA1",
        "LRG_199t1:r.100a>g",
    );
}

#[test]
fn display_drops_gene_selector_lrg_protein() {
    // Protein grammar has no selector slot (#310); the LRG protein form drops
    // it like any other `p.`.
    let v = parse_hgvs("LRG_199p1(BRCA1):p.(Arg8Gln)").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("BRCA1"));
    assert_eq!(v.to_string(), "LRG_199p1:p.(Arg8Gln)");
}

// =============================================================================
// Compound-ref + selector: `NC_X(NM_Y)(GENE):c.…`. The `NC_(NM_)` compound
// reference already names the transcript (the spec-correct specification), so
// `self` is the transcript `NM_` and the redundant `(GENE)` selector is
// dropped, leaving the canonical `NC_(NM_):c.…` form (#1051).
// =============================================================================

#[test]
fn display_drops_redundant_gene_selector_with_compound_ref() {
    let v = parse_hgvs("NC_000013.11(NM_004119.3)(FLT3):c.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("FLT3"));
    assert_eq!(
        v.to_string(),
        "NC_000013.11(NM_004119.3):c.100A>G",
        "the compound ref already names the transcript; the gene selector is redundant and dropped",
    );
}

#[test]
fn display_compound_ref_without_selector_unchanged() {
    let input = "NC_000013.11(NM_004119.3):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(gene_symbol_of(&v), None);
    assert_eq!(v.to_string(), input);
}

#[test]
fn display_drops_gene_selector_with_compound_ref_whose_inner_is_not_a_transcript() {
    // #1139: the drop is driven by the compound reference having already spent
    // its one permitted specification (`refseq.md:40-41`), NOT by the inner
    // accession being a transcript. `NG_012232.1` is a genomic reference, so the
    // #1051 `is_transcript_reference` gate alone would keep `(BRCA1)` and render
    // a doubled `(...)(...)` specification. The symbol stays on the struct.
    let v = parse_hgvs("NC_000013.11(NG_012232.1)(BRCA1):c.100A>G").unwrap();
    assert_eq!(gene_symbol_of(&v), Some("BRCA1"));
    assert_eq!(
        v.to_string(),
        "NC_000013.11(NG_012232.1):c.100A>G",
        "a compound reference already carries its one specification; the gene selector is dropped",
    );
}

#[test]
fn compound_ref_with_nontranscript_inner_display_is_idempotent() {
    // The dropped form must be a fixed point — re-parsing the rendered string
    // has to reproduce it byte-identically (this repo is output-stability
    // critical, and #1058 turns a non-idempotent Display into a real bug).
    let once = parse_hgvs("NC_000013.11(NG_012232.1)(BRCA1):c.100A>G")
        .unwrap()
        .to_string();
    let twice = parse_hgvs(&once).unwrap().to_string();
    assert_eq!(once, twice);
}

// =============================================================================
// Allele compound-form matrix. On a transcript reference the per-sub-variant
// selector is dropped (#1051), so the compact/expanded prefix loses `(GENE)`.
// Because that selector is never displayed on a transcript reference, the
// compaction decision is gene-agnostic there — a gene mismatch does not block
// compaction, keeping Display idempotent (#1058). On a non-transcript
// reference the selector is shown, so a mismatch still forces expanded form.
// =============================================================================

#[test]
fn display_drops_gene_selector_in_compact_cis_allele() {
    // Compact cis on a transcript ref: `ACC:c.[a;b]` (selector dropped).
    let input = "NM_000088.3(COL1A1):c.[100A>G;200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(v.to_string(), "NM_000088.3:c.[100A>G;200C>T]");
    if let HgvsVariant::Allele(a) = &v {
        for sub in &a.variants {
            assert_eq!(
                gene_symbol_of(sub),
                Some("COL1A1"),
                "field preserved on struct"
            );
        }
    } else {
        panic!("expected Allele; got {:?}", v);
    }
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
fn display_canonicalizes_expanded_same_acc_trans_to_compact_dropping_selector() {
    // Same accession + same gene + trans phase: ferro canonicalizes the
    // expanded `[v1];[v2]` input to the compact-trans form, and the selector
    // is dropped from the transcript-reference prefix (#1051). The field
    // survives on every sub-variant; reparse is byte-stable thereafter.
    let input = "[NM_000088.3(COL1A1):c.100A>G];[NM_000088.3(COL1A1):c.200C>T]";
    let v = parse_hgvs(input).unwrap();
    let displayed = v.to_string();
    assert_eq!(
        displayed, "NM_000088.3:c.[100A>G];[200C>T]",
        "matching selectors canonicalize to compact-trans; selector dropped on transcript ref",
    );

    let reparsed = parse_hgvs(&displayed).unwrap();
    assert_eq!(
        reparsed.to_string(),
        displayed,
        "reparse must be byte-stable"
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
fn display_drops_independent_gene_selectors_in_mixed_acc_trans() {
    // Different transcript accessions with different gene symbols — the
    // canonical compound-heterozygosity report shape. Both selectors are
    // dropped from Display (both accessions are transcripts), but the fields
    // survive on the sub-variants.
    let input = "[NM_000088.3(COL1A1):c.100A>G];[NM_000089.3(COL1A2):c.200C>T]";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(
        v.to_string(),
        "[NM_000088.3:c.100A>G];[NM_000089.3:c.200C>T]",
    );

    if let HgvsVariant::Allele(a) = &v {
        assert_eq!(a.variants.len(), 2);
        assert_eq!(gene_symbol_of(&a.variants[0]), Some("COL1A1"));
        assert_eq!(gene_symbol_of(&a.variants[1]), Some("COL1A2"));
    } else {
        panic!("expected Allele; got {:?}", v);
    }
}

#[test]
fn mismatched_gene_symbols_compact_on_transcript_ref() {
    // Same transcript accession, mismatched gene_symbol. On a transcript
    // reference the selector is never displayed (#1051), so the mismatch is
    // invisible in both the compact and the expanded form. Declining
    // compaction on that invisible mismatch made Display non-idempotent (the
    // expanded form reparses to all-`None` genes and then compacts), so the
    // compaction decision is gene-agnostic on a transcript reference (#1058).
    let input = "[NM_000088.3(COL1A1):c.100A>G];[NM_000088.3(COL1A2):c.200C>T]";
    let v = parse_hgvs(input).unwrap();
    let displayed = v.to_string();
    assert_eq!(
        displayed, "NM_000088.3:c.[100A>G];[200C>T]",
        "mismatched selectors on a transcript ref compact (selector not shown)",
    );

    // The parsed struct still carries both selectors for programmatic access,
    // even though Display suppresses them.
    if let HgvsVariant::Allele(a) = &v {
        assert_eq!(a.variants.len(), 2);
        assert_eq!(gene_symbol_of(&a.variants[0]), Some("COL1A1"));
        assert_eq!(gene_symbol_of(&a.variants[1]), Some("COL1A2"));
    } else {
        panic!("expected Allele; got {:?}", v);
    }

    // Display -> parse -> Display is now stable (#1058).
    assert_eq!(parse_hgvs(&displayed).unwrap().to_string(), displayed);
}

// =============================================================================
// Edge inputs: empty `()`, whitespace padding, casing, hyphen. The parser
// preserves the captured bytes verbatim on the struct; Display then applies
// the reference-type policy.
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
fn whitespace_padded_selector_captured_verbatim_then_dropped() {
    // ferro is bytes-in on parse: the captured gene_symbol keeps the padding.
    // On a transcript reference Display then drops it (#1051), so the output
    // is the canonical bare form.
    let input = "NM_000088.3( COL1A1 ):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(
        gene_symbol_of(&v),
        Some(" COL1A1 "),
        "parse captures verbatim"
    );
    assert_eq!(
        v.to_string(),
        "NM_000088.3:c.100A>G",
        "transcript ref drops selector"
    );
}

#[test]
fn lowercase_selector_captured_verbatim_then_dropped() {
    // ferro does not normalize case on the captured gene_symbol; a transcript
    // reference then drops it on Display.
    let input = "NM_000088.3(col1a1):c.100A>G";
    let v = parse_hgvs(input).unwrap();
    assert_eq!(
        gene_symbol_of(&v),
        Some("col1a1"),
        "parse captures verbatim"
    );
    assert_eq!(
        v.to_string(),
        "NM_000088.3:c.100A>G",
        "transcript ref drops selector"
    );
}

#[test]
fn hyphenated_selector_preserved_on_mitochondrial_reference() {
    // Mitochondrial genes carry hyphens (MT-ND1, MT-TL1) and have no
    // transcript reference, so the selector is the spec-sanctioned exception
    // and is preserved verbatim (refseq.md lines 197-203).
    assert_nontranscript_ref_preserves_selector("NC_012920.1(MT-TL1):m.3243A>G");
    assert_nontranscript_ref_preserves_selector("NC_012920.1(MT-TL1):n.14A>G");
}

// =============================================================================
// Normalize: the gene_symbol FIELD is preserved through normalization on every
// reference type (only Display differs by reference type). Multi-pass
// idempotency pins that gene_symbol does not drift on repeated normalize.
// =============================================================================

#[test]
fn normalize_preserves_gene_symbol_field_cds() {
    // MockProvider::with_test_data has NM_000088.3 (COL1A1), so this exercises
    // the real CDS normalize path.
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NM_000088.3(COL1A1):c.4A>G").unwrap();
    assert_eq!(gene_symbol_of(&variant), Some("COL1A1"));

    let normalized = normalizer.normalize(&variant).unwrap();

    // The gene_symbol field must survive normalization on the struct.
    assert_eq!(
        gene_symbol_of(&normalized),
        Some("COL1A1"),
        "normalize must not drop the gene_symbol field on CdsVariant",
    );

    // But Display drops it on the transcript reference (#1051).
    assert!(
        !normalized.to_string().contains("(COL1A1)"),
        "Display on a transcript reference must not emit the selector: {}",
        normalized
    );
}

#[test]
fn normalize_preserves_gene_symbol_field_genome() {
    // Genomic substitutions normalize without a transcript lookup; the field
    // is preserved AND Display emits it (non-transcript reference).
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let variant = parse_hgvs("NC_000017.11(BRCA1):g.43094466G>A").unwrap();
    assert_eq!(gene_symbol_of(&variant), Some("BRCA1"));

    let normalized = normalizer.normalize(&variant).unwrap();
    assert_eq!(
        gene_symbol_of(&normalized),
        Some("BRCA1"),
        "normalize must not drop the gene_symbol field on GenomeVariant",
    );
    assert!(
        normalized.to_string().contains("(BRCA1)"),
        "Display on a genomic reference must preserve the selector: {}",
        normalized
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
    // under the derived `Eq`, and must not alter `gene_symbol`.
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
// Re-parse identity. For a NON-transcript reference the selector survives
// Display, so parse → display → parse is identity-preserving. For a TRANSCRIPT
// reference Display intentionally drops the selector (#1051), so the reparse
// yields a variant whose `gene_symbol` is None — identity holds only on the
// struct's non-selector content.
// =============================================================================

#[test]
fn reparse_preserves_gene_symbol_on_nontranscript_reference() {
    // Only non-transcript references round-trip byte-stably with the selector.
    let cases: &[(&str, &str)] = &[
        ("NC_000013.11(FLT3):g.12345A>G", "FLT3"),
        ("NC_012920.1(MT-ND1):m.3460G>A", "MT-ND1"),
        ("LRG_199(BRCA1):g.100A>G", "BRCA1"),
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
        assert_eq!(
            first, second,
            "structural equality must survive parse → display → parse for {input}"
        );
    }
}

#[test]
fn reparse_drops_gene_symbol_on_transcript_reference() {
    // On a transcript reference the Display leg is lossy for the selector, so
    // the reparse observes None and the round-tripped variant differs from the
    // original only in `gene_symbol` (like re-reading a `p.` string).
    let cases: &[(&str, &str)] = &[
        ("NM_000088.3(COL1A1):c.459A>G", "COL1A1"),
        ("NR_046018.2(DDX11L1):n.100A>G", "DDX11L1"),
        ("ENST00000380152.7(BRCA2):c.100A>G", "BRCA2"),
        ("LRG_199t1(BRCA1):c.100A>G", "BRCA1"),
    ];

    for (input, gene) in cases {
        let first = parse_hgvs(input).unwrap();
        assert_eq!(gene_symbol_of(&first), Some(*gene));

        let second = parse_hgvs(&first.to_string()).unwrap();
        assert_eq!(
            gene_symbol_of(&second),
            None,
            "reparse of the dropped form must observe no selector for {input}"
        );
        assert_ne!(
            first, second,
            "the dropped selector makes the round-trip non-identity for {input}"
        );
    }
}

// =============================================================================
// Hash + Eq distinctness. The derived Hash/Eq still distinguishes `Some(g)`,
// `Some(g')`, and `None` at the struct level. On a NON-transcript reference the
// round-trip preserves that distinctness; this pins that the parser itself
// keeps the three variants distinct before any Display drop.
// =============================================================================

#[test]
fn parsed_variants_distinguish_selector_presence() {
    use std::collections::HashSet;

    let with_gene_a = parse_hgvs("NM_000088.3(COL1A1):c.100A>G").unwrap();
    let with_gene_b = parse_hgvs("NM_000088.3(COL1A2):c.100A>G").unwrap();
    let bare = parse_hgvs("NM_000088.3:c.100A>G").unwrap();

    // The three parsed structs are distinct under derived Hash/Eq.
    assert_ne!(with_gene_a, bare);
    assert_ne!(with_gene_a, with_gene_b);
    assert_ne!(with_gene_b, bare);

    let mut set: HashSet<HgvsVariant> = HashSet::new();
    set.insert(with_gene_a);
    set.insert(with_gene_b);
    set.insert(bare);
    assert_eq!(
        set.len(),
        3,
        "variants with distinct selector presence must hash distinctly"
    );
}

#[test]
fn round_trip_identity_holds_on_nontranscript_reference() {
    // On a genomic reference the selector survives Display, so the round-trip
    // preserves full structural identity for all three selector states.
    let with_gene_a = parse_hgvs("NC_000013.11(FLT3):g.100A>G").unwrap();
    let with_gene_b = parse_hgvs("NC_000013.11(KIT):g.100A>G").unwrap();
    let bare = parse_hgvs("NC_000013.11:g.100A>G").unwrap();

    let a_rt = parse_hgvs(&with_gene_a.to_string()).unwrap();
    let b_rt = parse_hgvs(&with_gene_b.to_string()).unwrap();
    let bare_rt = parse_hgvs(&bare.to_string()).unwrap();

    assert_eq!(with_gene_a, a_rt);
    assert_eq!(with_gene_b, b_rt);
    assert_eq!(bare, bare_rt);

    assert_ne!(a_rt, bare_rt);
    assert_ne!(a_rt, b_rt);
}
