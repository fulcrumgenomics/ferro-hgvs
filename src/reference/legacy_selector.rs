//! Resolving legacy LOVD gene-model selectors (`GENE`, `GENE_v001`) on a
//! genomic reference (`NG_`/`NC_`/`LRG_`) to a transcript accession (#500/#637).
//!
//! HGVS strongly prefers a transcript accession in the selector slot of a
//! genomic-reference coding variant (`NG_012337.1(NM_003002.2):c.…`), not a
//! gene symbol (`refseq.md:38-42`). Inputs that use the legacy LOVD
//! gene-model selector — `NG_012337.1(TIMM8B_v001):c.…` — should resolve to the
//! gene's transcript.
//!
//! The `_v00N` numbering is positional on the RefSeqGene record. NCBI's
//! `LRG_RefSeqGene` association table flags each gene's *reference standard*
//! transcript (`Category == "reference standard"`); for the ~96% of genes with
//! a single reference-standard `NM_`, that transcript is the LOVD `_v001`, so a
//! gene-symbol → reference-standard `NM_` map resolves `_v001`/bare selectors
//! authoritatively.
//!
//! It does **not** resolve everything: ~4.3% of genes (e.g. APC, ABL1, BCL2)
//! carry *multiple* reference-standard transcripts on one record, and the table
//! does not encode the LOVD `_v` ordering (its row order is neither the `_v`
//! order nor accession-numeric). So a gene with more than one reference-standard
//! `NM_` is **excluded from the map** — it declines (preserves the input)
//! rather than guess a transcript that could shift the `c.` numbering. Higher
//! locus versions (`_v002`…), which name the 2nd/3rd transcript on exactly those
//! multi-transcript records, are likewise declined. Reliable `_vNNN` resolution
//! for multi-transcript genes needs the `NG_` GenBank feature order (the
//! definitional `_v` source) — a future enhancement.

use rustc_hash::FxHashMap;

/// Parse NCBI's `LRG_RefSeqGene` association table into a
/// `gene Symbol (upper-case) → reference-standard NM_ accession` map.
///
/// The tab-separated file's columns are
/// `tax_id  GeneID  Symbol  RSG  LRG  RNA  t  Protein  p  Category`.
/// Only rows whose `Category` marks the gene's reference standard and whose
/// `RNA` is an `NM_` transcript are collected. A gene is kept **only when it has
/// exactly one** such transcript: that is its unambiguous LOVD `_v001`. Genes
/// with multiple reference-standard transcripts are dropped (the table does not
/// encode which is `_v001`), so they decline rather than resolve to an
/// arbitrary file-order pick.
pub fn parse_refseqgene_summary(tsv: &str) -> FxHashMap<String, String> {
    // Collect all distinct reference-standard NM_ transcripts per gene first;
    // only genes with a single candidate are resolvable (see module docs).
    let mut by_gene: FxHashMap<String, Vec<String>> = FxHashMap::default();
    for line in tsv.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 10 {
            continue;
        }
        let symbol = cols[2];
        let rna = cols[5];
        let category = cols[9];
        if !category.contains("reference standard") || !rna.starts_with("NM_") {
            continue;
        }
        let entry = by_gene.entry(symbol.to_ascii_uppercase()).or_default();
        if !entry.iter().any(|existing| existing == rna) {
            entry.push(rna.to_string());
        }
    }
    by_gene
        .into_iter()
        .filter_map(|(gene, mut nms)| (nms.len() == 1).then(|| (gene, nms.pop().unwrap())))
        .collect()
}

/// Split a legacy gene-model selector into its gene symbol and optional locus
/// version: `"TIMM8B_v001" → ("TIMM8B", Some(1))`, `"TIMM8B" → ("TIMM8B", None)`,
/// `"SDHD_v002" → ("SDHD", Some(2))`. The `_v<n>` (transcript) and `_i<n>`
/// (protein isoform) suffixes are recognized; anything else is treated as part
/// of the gene symbol.
pub fn split_legacy_gene_selector(selector: &str) -> (&str, Option<u32>) {
    if let Some(idx) = selector.rfind('_') {
        let suffix = &selector[idx + 1..];
        let mut chars = suffix.chars();
        if let Some(first) = chars.next() {
            let digits = &suffix[1..];
            if (first == 'v' || first == 'i')
                && !digits.is_empty()
                && digits.bytes().all(|b| b.is_ascii_digit())
            {
                return (&selector[..idx], digits.parse::<u32>().ok());
            }
        }
    }
    (selector, None)
}

/// Resolve a legacy gene-model selector to a transcript accession via `lookup`
/// (a gene Symbol, upper-cased, → transcript accession map).
///
/// Returns the gene's reference-standard transcript for a bare gene or `_v001`;
/// declines (`None`) for higher locus versions (`_v002`…), which the
/// reference-standard map cannot disambiguate, and for an unknown gene.
pub fn resolve_legacy_selector_in<F>(selector: &str, lookup: F) -> Option<String>
where
    F: Fn(&str) -> Option<String>,
{
    let (gene, version) = split_legacy_gene_selector(selector);
    if version.is_some_and(|n| n != 1) {
        return None;
    }
    lookup(&gene.to_ascii_uppercase())
}

/// Parent-aware legacy selector resolution (#792). For a bare/`_v001` selector
/// on an `NG_` parent, prefer the transcript that exact parent uniquely hosts;
/// otherwise (no parent, gene not uniquely hosted there, or `_v002+`) fall back
/// to the global reference-standard `lookup`. `_v002+` still declines via the
/// existing rule.
pub fn resolve_legacy_selector_with_parent<H, G>(
    selector: &str,
    ng_parent: Option<&str>,
    hosted_unique: H,
    lookup: G,
) -> Option<String>
where
    H: Fn(&str, &str) -> Option<String>,
    G: Fn(&str) -> Option<String>,
{
    let (gene, version) = split_legacy_gene_selector(selector);
    if version.is_some_and(|n| n != 1) {
        return None;
    }
    let gene_upper = gene.to_ascii_uppercase();
    if let Some(ng) = ng_parent {
        if let Some(tx) = hosted_unique(ng, &gene_upper) {
            return Some(tx);
        }
    }
    lookup(&gene_upper)
}

#[cfg(test)]
mod tests {
    use super::*;

    // A representative slice of the real `LRG_RefSeqGene` table: a primary gene
    // (SDHD on its own record), a secondary gene on a different record (TIMM8B),
    // non-reference-standard rows, and a non-NM_ (NR_) reference row to skip.
    const SAMPLE: &str = "\
#tax_id\tGeneID\tSymbol\tRSG\tLRG\tRNA\tt\tProtein\tp\tCategory
9606\t6392\tSDHD\tNG_012337.3\tLRG_9\tNR_077060.2\t\t\t\taligned: Selected
9606\t6392\tSDHD\tNG_012337.3\tLRG_9\tNM_003002.4\tt1\tNP_002993.1\t\treference standard
9606\t1788\tTIMM8B\tNG_033145.1\t\tNM_012459.4\t\tNP_036591.3\t\treference standard
9606\t1788\tTIMM8B\tNG_033145.1\t\tNR_028383.2\t\t\t\taligned: Selected
9606\t675\tBRCA2\tNG_012772.3\t\tNM_000059.3\tt1\tNP_000050.2\t\treference standard
9606\t675\tBRCA2\tNG_012772.3\t\tNM_000059.4\t\tNP_000050.3\t\taligned: Selected";

    #[test]
    fn parses_reference_standard_nm_per_gene() {
        let map = parse_refseqgene_summary(SAMPLE);
        assert_eq!(map.get("SDHD").map(String::as_str), Some("NM_003002.4"));
        assert_eq!(map.get("TIMM8B").map(String::as_str), Some("NM_012459.4"));
        assert_eq!(map.get("BRCA2").map(String::as_str), Some("NM_000059.3"));
        // The NR_ reference row and the "aligned: Selected" rows are not kept.
        assert_eq!(map.len(), 3);
    }

    #[test]
    fn multi_reference_standard_gene_is_excluded() {
        // APC carries three reference-standard transcripts on one record; the
        // table doesn't encode which is `_v001`, so APC must be DROPPED (decline)
        // rather than resolve to an arbitrary file-order pick.
        let sample = "\
#tax_id\tGeneID\tSymbol\tRSG\tLRG\tRNA\tt\tProtein\tp\tCategory
9606\t324\tAPC\tNG_008481.4\tLRG_130\tNM_000038.6\t\tNP_000029.2\t\treference standard
9606\t324\tAPC\tNG_008481.4\tLRG_130\tNM_001127510.3\t\tNP_001120982.2\t\treference standard
9606\t324\tAPC\tNG_008481.4\tLRG_130\tNM_001127511.3\t\tNP_001120983.2\t\treference standard
9606\t6392\tSDHD\tNG_012337.3\tLRG_9\tNM_003002.4\tt1\tNP_002993.1\t\treference standard";
        let map = parse_refseqgene_summary(sample);
        assert!(
            !map.contains_key("APC"),
            "ambiguous multi-transcript gene must be dropped"
        );
        assert_eq!(map.get("SDHD").map(String::as_str), Some("NM_003002.4"));
        let lookup = |g: &str| map.get(g).cloned();
        assert_eq!(resolve_legacy_selector_in("APC_v001", lookup), None);
        assert_eq!(resolve_legacy_selector_in("APC", lookup), None);
    }

    #[test]
    fn splits_locus_versions() {
        assert_eq!(
            split_legacy_gene_selector("TIMM8B_v001"),
            ("TIMM8B", Some(1))
        );
        assert_eq!(split_legacy_gene_selector("SDHD_v002"), ("SDHD", Some(2)));
        assert_eq!(split_legacy_gene_selector("BRCA2_i001"), ("BRCA2", Some(1)));
        assert_eq!(split_legacy_gene_selector("TIMM8B"), ("TIMM8B", None));
        // A trailing `_` group that isn't `_v<n>`/`_i<n>` stays part of the gene.
        assert_eq!(split_legacy_gene_selector("C9orf72"), ("C9orf72", None));
        assert_eq!(split_legacy_gene_selector("GENE_x1"), ("GENE_x1", None));
    }

    #[test]
    fn resolves_v001_and_bare_declines_higher_and_unknown() {
        let map = parse_refseqgene_summary(SAMPLE);
        let lookup = |g: &str| map.get(g).cloned();
        // _v001 and bare gene resolve to the reference standard (case-insensitive).
        assert_eq!(
            resolve_legacy_selector_in("TIMM8B_v001", lookup),
            Some("NM_012459.4".to_string())
        );
        assert_eq!(
            resolve_legacy_selector_in("timm8b", lookup),
            Some("NM_012459.4".to_string())
        );
        assert_eq!(
            resolve_legacy_selector_in("SDHD", lookup),
            Some("NM_003002.4".to_string())
        );
        // _v002 is not expressible from the reference-standard map → decline.
        assert_eq!(resolve_legacy_selector_in("SDHD_v002", lookup), None);
        // Unknown gene → decline.
        assert_eq!(resolve_legacy_selector_in("NOTAGENE_v001", lookup), None);
    }

    #[test]
    fn ng_parent_hosted_wins_over_global_for_v001() {
        // hosted: NG_012337.1 hosts TIMM8B → NM_012459.2 (the parent-relative answer)
        let hosted = |ng: &str, g: &str| -> Option<String> {
            (ng == "NG_012337.1" && g == "TIMM8B").then(|| "NM_012459.2".to_string())
        };
        // global reference-standard map would say NM_012459.4 (the bug)
        let global =
            |g: &str| -> Option<String> { (g == "TIMM8B").then(|| "NM_012459.4".to_string()) };
        // _v001 on the NG parent → hosted wins
        assert_eq!(
            resolve_legacy_selector_with_parent("TIMM8B_v001", Some("NG_012337.1"), hosted, global),
            Some("NM_012459.2".to_string())
        );
        // _v002 → still declines (no positional)
        assert_eq!(
            resolve_legacy_selector_with_parent("TIMM8B_v002", Some("NG_012337.1"), hosted, global),
            None
        );
        // no NG parent → falls back to global
        assert_eq!(
            resolve_legacy_selector_with_parent("TIMM8B", None, hosted, global),
            Some("NM_012459.4".to_string())
        );
        // NG parent present but gene not hosted there → falls back to global
        let nohost = |_: &str, _: &str| None;
        assert_eq!(
            resolve_legacy_selector_with_parent("TIMM8B", Some("NG_999999.9"), nohost, global),
            Some("NM_012459.4".to_string())
        );
    }
}
