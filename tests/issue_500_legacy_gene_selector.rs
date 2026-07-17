//! Issue #500/#637: resolve a legacy LOVD gene-model selector on a genomic
//! reference to the spec-preferred transcript accession.
//!
//! HGVS prefers a transcript accession in the selector slot of a
//! genomic-reference coding variant (`NG_012337.1(NM_003002.2):c.…`), not a gene
//! symbol (`refseq.md:38-42`). `normalize` rewrites a legacy gene-model selector
//! — `NG_/NC_/LRG(GENE[_v001]):c.…` — to the gene's reference-standard
//! transcript, resolved Symbol-keyed (so a secondary gene cited under another
//! gene's RefSeqGene record resolves via its own reference standard), while
//! leaving `NM_(GENE)` and unresolvable selectors untouched (#121).

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{ManeStatus, Strand};
use ferro_hgvs::reference::Transcript;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// A `MockProvider` carrying the SDHD (`NM_003002.4`) and TIMM8B (`NM_012459.4`)
/// transcripts and their legacy gene-model → reference-standard mappings. The
/// sequences are all-`A` filler long enough to validate the `c.` coordinate
/// used below.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    for (id, gene) in [("NM_003002.4", "SDHD"), ("NM_012459.4", "TIMM8B")] {
        let tx = Transcript::new(
            id.to_string(),
            Some(gene.to_string()),
            Strand::Plus,
            Some("A".repeat(300)),
            Some(1),
            Some(300),
            Vec::new(),
            None,
            None,
            None,
            Default::default(),
            ManeStatus::Select,
            None,
            None,
        );
        provider.add_transcript(tx);
        provider.add_legacy_gene_model(gene, id);
    }
    provider
}

fn normalize_str(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"));
    format!("{normalized}")
}

#[test]
fn ng_gene_v001_selector_rewritten_to_nm() {
    assert_eq!(
        normalize_str(provider(), "NG_012337.3(SDHD_v001):c.92A>T"),
        "NG_012337.3(NM_003002.4):c.92A>T"
    );
}

#[test]
fn ng_bare_gene_selector_rewritten_to_nm() {
    // A bare gene name (no `_vNNN`) resolves to the same reference-standard transcript.
    assert_eq!(
        normalize_str(provider(), "NG_012337.3(SDHD):c.92A>T"),
        "NG_012337.3(NM_003002.4):c.92A>T"
    );
}

#[test]
fn nc_gene_selector_rewritten_to_nm() {
    // `is_genomic_ref` also accepts an `NC_` genomic reference in the selector
    // slot; the gene resolves via the same Symbol-keyed reference standard and
    // the `NC_` wrapper is preserved as the genomic context.
    assert_eq!(
        normalize_str(provider(), "NC_000011.10(SDHD):c.92A>T"),
        "NC_000011.10(NM_003002.4):c.92A>T"
    );
}

#[test]
fn lrg_gene_selector_rewritten_to_nm() {
    // An `LRG_` genomic reference is likewise in scope; the gene resolves to the
    // reference-standard transcript and the `LRG_` wrapper is preserved.
    assert_eq!(
        normalize_str(provider(), "LRG_9(SDHD):c.92A>T"),
        "LRG_9(NM_003002.4):c.92A>T"
    );
}

#[test]
fn secondary_gene_resolves_via_its_own_reference_standard() {
    // TIMM8B is a *secondary* gene annotated on the SDHD RefSeqGene record
    // (NG_012337); the Symbol-keyed lookup resolves it via TIMM8B's own
    // reference-standard transcript, not anything scoped to NG_012337.
    assert_eq!(
        normalize_str(provider(), "NG_012337.3(TIMM8B_v001):c.92A>T"),
        "NG_012337.3(NM_012459.4):c.92A>T"
    );
}

#[test]
fn rewrite_is_idempotent() {
    // Re-normalizing the rewritten output is a no-op: the selector already names
    // a transcript (it carries a genomic context), so no resolution is attempted.
    let once = normalize_str(provider(), "NG_012337.3(SDHD_v001):c.92A>T");
    let twice = normalize_str(provider(), &once);
    assert_eq!(once, twice);
    assert_eq!(twice, "NG_012337.3(NM_003002.4):c.92A>T");
}

#[test]
fn nm_gene_selector_is_not_resolved_away() {
    // `NM_(GENE)` is out of scope for legacy-selector resolution: the accession
    // is already a named transcript, so it is NOT rewritten to a different
    // accession. The gene symbol is dropped from Display on the transcript
    // reference (#1051), so the canonical output is the bare `NM_…:c.…` form.
    let out = normalize_str(provider(), "NM_003002.4(SDHD):c.92A>T");
    assert_eq!(
        out, "NM_003002.4:c.92A>T",
        "NM_ accession kept (not resolved away); gene selector dropped on the transcript ref: {out}"
    );
}

#[test]
fn higher_locus_version_is_declined_and_preserved() {
    // `_v002` is not expressible from the reference-standard map → no rewrite;
    // the input selector is preserved (never an error on the unresolvable path).
    let out = normalize_str(provider(), "NG_012337.3(SDHD_v002):c.92A>T");
    assert_eq!(out, "NG_012337.3(SDHD_v002):c.92A>T");
}

#[test]
fn unknown_gene_without_mapping_is_preserved() {
    // An empty provider (no legacy gene-model summary) resolves nothing; the
    // input selector is preserved verbatim, never an error.
    let out = normalize_str(MockProvider::new(), "NG_012337.3(SDHD_v001):c.92A>T");
    assert_eq!(out, "NG_012337.3(SDHD_v001):c.92A>T");
}
