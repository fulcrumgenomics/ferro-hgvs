//! Issue #1146 — compound reference syntax `outer(inner)` was recognised only
//! for RefSeq-shaped, Ensembl, and LRG **outer** accessions, and only for those
//! same families as the **inner** token. A custom (SAM-refname) accession, or an
//! assembly/chromosome reference, could not carry a transcript specification.
//!
//! Two consequences, both fixed here:
//!
//! 1. **Asymmetric acceptance.** `NC_000013.11(NM_004119.3)(FLT3):c.…` parsed,
//!    but `MYSEQ.1(MYTX.1)(GENE1):c.…` and `Template(Template-gene.1)(GENE1):c.…`
//!    failed with `Expected ':' after accession`, purely because of the prefix
//!    whitelist in `looks_like_accession_start`.
//!
//! 2. **Display output was string-stable but not struct-stable.** The
//!    spec-conformant string the projector emits for a custom reference (#1139),
//!    `Template(Template-gene.1):c.5A>G`, re-parsed with the *transcript id* in
//!    `gene_symbol` and `genomic_context: None` — so a consumer that re-parsed
//!    `projection.c_name` and read `.gene_symbol` got a transcript accession. The
//!    same held for `GRCh38(chr1)(NM_004006.2):c.…` (#1145).
//!
//! The disambiguator between "inner is an accession" and "inner is a gene-symbol
//! selector" is the **version suffix**: HGNC gene symbols do not carry a
//! `.<digits>` version, so `MYREF_SEQ(GENE1)` keeps its long-standing
//! gene-symbol reading while `MYREF_SEQ(MYTX.1)` is read as a compound
//! reference.

use ferro_hgvs::hgvs::parser::variant::parse_variant;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

/// The parsed variant's `gene_symbol`, regardless of variant kind.
fn gene_symbol_of(variant: &HgvsVariant) -> Option<&str> {
    match variant {
        HgvsVariant::Genome(v) => v.gene_symbol.as_deref(),
        HgvsVariant::Cds(v) => v.gene_symbol.as_deref(),
        HgvsVariant::Tx(v) => v.gene_symbol.as_deref(),
        _ => None,
    }
}

/// `(genomic_context, transcript_accession, gene_symbol)` for an input.
fn parts(input: &str) -> (Option<String>, String, Option<String>) {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?} failed: {e:?}"));
    let acc = v.accession().expect("accession").clone();
    (
        acc.genomic_context.as_ref().map(|c| c.full()),
        acc.transcript_accession(),
        gene_symbol_of(&v).map(str::to_string),
    )
}

// =============================================================================
// (1) A custom outer accession accepts a versioned inner as a specification.
// =============================================================================

#[test]
fn custom_outer_with_versioned_inner_is_a_compound_reference() {
    let (context, transcript, gene) = parts("MYSEQ.1(MYTX.1):c.100A>G");
    assert_eq!(context.as_deref(), Some("MYSEQ.1"));
    assert_eq!(transcript, "MYTX.1");
    assert_eq!(
        gene, None,
        "the inner token is an accession, not a selector"
    );
}

#[test]
fn custom_outer_with_versioned_inner_and_selector_parses() {
    // Previously `Expected ':' after accession`.
    let (context, transcript, gene) = parts("MYSEQ.1(MYTX.1)(GENE1):c.100A>G");
    assert_eq!(context.as_deref(), Some("MYSEQ.1"));
    assert_eq!(transcript, "MYTX.1");
    assert_eq!(gene.as_deref(), Some("GENE1"));
}

#[test]
fn projector_style_custom_reference_parses_as_compound() {
    // The exact shape `VariantProjector` emits for a custom reference (#1139).
    let (context, transcript, gene) = parts("Template(Template-gene.1):c.5A>G");
    assert_eq!(context.as_deref(), Some("Template"));
    assert_eq!(
        transcript, "Template-gene.1",
        "the transcript id must land in the accession, not in gene_symbol",
    );
    assert_eq!(gene, None);
}

#[test]
fn projector_style_custom_reference_with_stacked_selector_parses() {
    // The pre-#1139 malformed output must now parse (and, per #1139, Display
    // drops the stacked selector rather than reproducing it).
    let (context, transcript, gene) = parts("Template(Template-gene.1)(GENE1):c.5A>G");
    assert_eq!(context.as_deref(), Some("Template"));
    assert_eq!(transcript, "Template-gene.1");
    assert_eq!(gene.as_deref(), Some("GENE1"));
}

// =============================================================================
// (2) An assembly/chromosome outer reference accepts a specification too (the
// #1145 probe: the transcript used to land in `gene_symbol`).
// =============================================================================

#[test]
fn assembly_outer_with_transcript_inner_is_a_compound_reference() {
    let (context, transcript, gene) = parts("GRCh38(chr1)(NM_004006.2):c.100A>G");
    assert_eq!(context.as_deref(), Some("GRCh38(chr1)"));
    assert_eq!(transcript, "NM_004006.2");
    assert_eq!(
        gene, None,
        "a transcript specification must not be modelled as a gene symbol",
    );
}

// =============================================================================
// (3) Every outer dispatch path admits a custom versioned inner, not just the
// custom/assembly ones — the shared `try_compound_suffix` is reached from the
// standard and Ensembl parsers too.
// =============================================================================

#[test]
fn refseq_outer_with_custom_versioned_inner_is_a_compound_reference() {
    let (context, transcript, gene) = parts("NC_000013.11(MYTX.1):c.100A>G");
    assert_eq!(context.as_deref(), Some("NC_000013.11"));
    assert_eq!(transcript, "MYTX.1");
    assert_eq!(gene, None);
}

#[test]
fn ensembl_outer_with_custom_versioned_inner_is_a_compound_reference() {
    let (context, transcript, gene) = parts("ENSG00000139618.15(MYTX.1):c.100A>G");
    assert_eq!(context.as_deref(), Some("ENSG00000139618.15"));
    assert_eq!(transcript, "MYTX.1");
    assert_eq!(gene, None);
}

#[test]
fn alternate_assembly_and_wgs_outers_admit_a_specification() {
    // `AC_` / `NZ_` are accepted by `is_valid_compound_outer`'s explicit arm,
    // which this PR restructured into an early return — pin that the arm still
    // fires.
    let (context, transcript, _) = parts("AC_000001.1(NM_004006.2):c.100A>G");
    assert_eq!(context.as_deref(), Some("AC_000001.1"));
    assert_eq!(transcript, "NM_004006.2");
    let (context, transcript, _) = parts("NZ_CP007265.1(NM_004006.2):c.100A>G");
    assert_eq!(context.as_deref(), Some("NZ_CP007265.1"));
    assert_eq!(transcript, "NM_004006.2");
}

#[test]
fn lrg_outer_admits_a_specification_on_every_inner_family() {
    // A bare `LRG_<N>` is a genomic reference, so it takes a specification like
    // any other — including the custom versioned inner this PR adds. (The `ferro
    // parse` CLI rejects these for a *downstream* reason unrelated to parsing,
    // so assert at the parser level, not through the CLI.)
    for (input, inner) in [
        ("LRG_199(NM_004006.2):c.100A>G", "NM_004006.2"),
        ("LRG_199(MYTX.1):c.100A>G", "MYTX.1"),
        ("LRG_199(LRG_199t1):c.100A>G", "LRG_199t1"),
    ] {
        let (context, transcript, gene) = parts(input);
        assert_eq!(context.as_deref(), Some("LRG_199"), "context for {input:?}");
        assert_eq!(transcript, inner, "inner for {input:?}");
        assert_eq!(gene, None, "gene for {input:?}");
    }
}

// =============================================================================
// Non-regression: an UNVERSIONED inner token is still a gene-symbol selector.
// This is the ambiguity the version suffix resolves.
// =============================================================================

#[test]
fn custom_outer_with_unversioned_inner_is_still_a_gene_selector() {
    // PR #70's motivating case must not change meaning.
    let (context, transcript, gene) = parts("MYREF_SEQ(GENE1):c.100A>G");
    assert_eq!(context, None);
    assert_eq!(transcript, "MYREF_SEQ");
    assert_eq!(gene.as_deref(), Some("GENE1"));
}

#[test]
fn assembly_ref_with_unversioned_inner_is_still_a_gene_selector() {
    // #1145 depends on this reading.
    let (context, transcript, gene) = parts("GRCh38(chr1)(BRCA1):g.100A>G");
    assert_eq!(context, None);
    assert_eq!(gene.as_deref(), Some("BRCA1"));
    assert_eq!(transcript, "GRCh38(chr1)");
}

#[test]
fn refseq_compound_reference_is_unchanged_by_the_refactor() {
    // The pre-existing `NC_(NM_)` path must be untouched by extracting the rule
    // into the shared `try_compound_suffix` and by the nesting cap.
    let (context, transcript, gene) = parts("NC_000013.11(NM_004119.3):c.100A>G");
    assert_eq!(context.as_deref(), Some("NC_000013.11"));
    assert_eq!(transcript, "NM_004119.3");
    assert_eq!(gene, None);
}

#[test]
fn nested_compound_reference_is_still_not_a_compound() {
    // #1151 nesting cap: a nested form never parsed as a compound reference
    // (the inner-is-compound rejection), and still does not. Pinning it proves
    // the cap changed no behavior.
    assert!(
        parse_hgvs("NC_000013.11(NC_000013.11(NM_004119.3)):c.100A>G").is_err(),
        "a nested compound reference must not parse",
    );
}

#[test]
fn refseq_outer_with_unversioned_inner_is_still_a_gene_selector() {
    let (context, transcript, gene) = parts("NC_000023.10(DMD):c.100A>G");
    assert_eq!(context, None);
    assert_eq!(transcript, "NC_000023.10");
    assert_eq!(gene.as_deref(), Some("DMD"));
}

// =============================================================================
// Display round-trips on the newly-parseable forms.
// =============================================================================

#[test]
fn custom_compound_reference_round_trips() {
    // Cross-check against the generic `parse_variant` — a *different* internal
    // code path than the fast-path-fronted `parse_hgvs` under test — rather than
    // deriving the expectation from the same function (repo test policy).
    let input = "Template(Template-gene.1):c.5A>G";
    let via_fast_path = parse_hgvs(input).unwrap();
    let via_generic = parse_variant(input).expect("generic parser accepts the form");
    assert_eq!(
        via_fast_path, via_generic,
        "the two parse paths must agree on the compound reference",
    );
    assert_eq!(via_fast_path.to_string(), input);
    assert_eq!(via_generic.to_string(), input);
}

#[test]
fn custom_compound_display_is_idempotent() {
    const EXPECTED: &str = "MYSEQ.1(MYTX.1):c.100A>G";
    let once = parse_hgvs("MYSEQ.1(MYTX.1)(GENE1):c.100A>G")
        .unwrap()
        .to_string();
    assert_eq!(
        once, EXPECTED,
        "#1139: the stacked selector is dropped, leaving one specification",
    );

    // Assert the reparsed *structure*, not just the string: a render that
    // round-trips byte-identically can still re-parse into a different shape —
    // that is precisely the bug #1146 fixes, so string equality alone would not
    // pin it.
    let reparsed = parse_hgvs(&once).unwrap();
    let acc = reparsed.accession().expect("accession");
    assert_eq!(
        acc.genomic_context.as_ref().map(|c| c.full()).as_deref(),
        Some("MYSEQ.1")
    );
    assert_eq!(acc.transcript_accession(), "MYTX.1");
    assert_eq!(gene_symbol_of(&reparsed), None);
    assert_eq!(reparsed.to_string(), EXPECTED, "second render must match");
}

/// An over-wide version suffix must not be silently dropped.
///
/// `Accession::version` is a `u32`, and the parse used to be
/// `digits.parse::<u32>().ok()` — so a digit run that did not fit became
/// `None`, splitting the token and rendering `MYTX.4294967296` back as `MYTX`.
/// The string parsed "successfully" and named a *different reference sequence*,
/// the same class of silent rewrite this file's compound handling exists to
/// prevent. An unrepresentable run is now simply not a version.
///
/// The two positions resolve it differently, and both are lossless:
///
/// * **outer** — the whole token stays the accession name, so `version` is
///   `None` and the digits survive in the name;
/// * **inner** — not being version-shaped, it no longer satisfies
///   `versioned_inner_span`, so the parens fall through to the gene-symbol
///   reading. Structurally that is the shape #1146 complains about, but it is
///   the better of the two available answers while `version` is a `u32`:
///   the alternative renamed the reference. Pinned so the choice is deliberate.
#[test]
fn an_unrepresentable_version_is_kept_in_the_accession_name() {
    // Outer: kept in the name, no version parsed.
    let v = parse_hgvs("MYTX.4294967296:c.100A>G").expect("parses");
    assert_eq!(v.to_string(), "MYTX.4294967296:c.100A>G");
    let accession = v.accession().expect("accession");
    assert_eq!(accession.version, None, "over-wide run stored as a version");
    assert!(
        accession.to_string().contains("4294967296"),
        "digits lost from the accession name: {accession}"
    );

    // Inner: falls back to the gene-symbol reading, string intact either way.
    for input in [
        "MYSEQ.1(MYTX.4294967296):c.100A>G",
        "MYSEQ.1(MYTX.99999999999999999999):c.100A>G",
    ] {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} should parse: {e}"));
        assert_eq!(v.to_string(), input, "version dropped from {input}");
        assert!(
            v.accession().expect("accession").genomic_context.is_none(),
            "{input} should not be read as a compound reference"
        );
    }
}

/// A representable version still parses as one, so the guard above did not
/// simply stop recognising versions.
///
/// Asserted on the `version` **field**, not only the rendering: keeping the
/// digits in the accession *name* renders byte-identically, so a string
/// round-trip cannot tell "parsed as a version" from "kept as text" — which is
/// the whole distinction this fix turns on.
#[test]
fn a_representable_version_is_still_a_version() {
    for (input, want) in [
        ("MYTX.5:c.100A>G", 5u32),
        ("MYSEQ.1(MYTX.5):c.100A>G", 5),
        // `u32::MAX` itself is representable and must still split off.
        ("MYSEQ.1(MYTX.4294967295):c.100A>G", u32::MAX),
    ] {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} should parse: {e}"));
        let accession = v.accession().expect("accession");
        assert_eq!(
            accession.version,
            Some(want),
            "version not parsed into the field for {input}"
        );
        assert_eq!(v.to_string(), input);
    }
}

/// A zero-padded version is not a version **on a custom accession**.
///
/// `Accession::version` is a `u32`, which cannot represent the padding, so
/// splitting `MYTX.007` into name + `Some(7)` re-rendered it as `MYTX.7` — a
/// different reference sequence, the same silent rewrite as the overflow case
/// above. Keeping the token whole is lossless.
///
/// A single `0` is canonical and still parses as a version.
#[test]
fn a_zero_padded_version_is_kept_in_a_custom_accession_name() {
    for input in [
        "MYTX.01:c.100A>G",
        "MYTX.007:c.100A>G",
        "MYSEQ.1(MYTX.007):c.100A>G",
    ] {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} should parse: {e}"));
        assert_eq!(v.to_string(), input, "zero-padded version rewritten");
    }

    let v = parse_hgvs("MYTX.0:c.100A>G").expect("parses");
    assert_eq!(v.accession().expect("accession").version, Some(0));
    assert_eq!(v.to_string(), "MYTX.0:c.100A>G");
}

/// A RefSeq accession is deliberately NOT covered by that rule.
///
/// A RefSeq version is a number whose canonical rendering is unpadded, so
/// `NM_000088.03` → `NM_000088.3` is a correction rather than a rewrite: the
/// accession still names the same sequence. A custom SAM refname is an opaque
/// identifier, where the same edit changes *which* sequence is named. Pinned so
/// the asymmetry stays deliberate.
#[test]
fn a_zero_padded_refseq_version_is_still_normalized() {
    let v = parse_hgvs("NM_000088.03:c.100A>G").expect("parses");
    assert_eq!(v.to_string(), "NM_000088.3:c.100A>G");
    assert_eq!(v.accession().expect("accession").version, Some(3));
}
