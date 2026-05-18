//! Audit for #81 § I1 — LRG references round-trip across edit types.
//!
//! LRG accessions come in three flavors:
//! - `LRG_<n>` — genomic; pairs with `g.`
//! - `LRG_<n>t<m>` — transcript; pairs with `c.`, `n.`, `r.`
//! - `LRG_<n>p<m>` — protein; pairs with `p.`
//!
//! Pin the round-trip of each LRG flavor across the edit types the
//! parser already supports. Cases known to be broken for non-LRG
//! reasons (predicted-wrapper `g.(=)` / `c.(=)` / `r.(0)` and protein
//! bracketed insertion `p.…ins[AA;AA]`) are tracked elsewhere and are
//! deliberately omitted from this audit.

use ferro_hgvs::parse_hgvs;

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

// =============================================================================
// SECTION 1 — `LRG_<n>` genomic (g.) across edit types
// =============================================================================

mod lrg_genomic {
    use super::*;

    #[test]
    fn sub_round_trips() {
        assert_round_trips("LRG_199:g.100A>G");
    }

    #[test]
    fn del_round_trips() {
        assert_round_trips("LRG_199:g.100_105del");
    }

    #[test]
    fn dup_round_trips() {
        assert_round_trips("LRG_199:g.100_105dup");
    }

    #[test]
    fn ins_round_trips() {
        assert_round_trips("LRG_199:g.100_101insAT");
    }

    #[test]
    fn inv_round_trips() {
        assert_round_trips("LRG_199:g.100_105inv");
    }

    #[test]
    fn delins_round_trips() {
        assert_round_trips("LRG_199:g.100_105delinsACGT");
    }

    #[test]
    fn repeat_round_trips() {
        assert_round_trips("LRG_199:g.100_102[5]");
    }

    #[test]
    fn identity_round_trips() {
        assert_round_trips("LRG_199:g.=");
    }

    #[test]
    fn unknown_round_trips() {
        assert_round_trips("LRG_199:g.?");
    }
}

// =============================================================================
// SECTION 2 — `LRG_<n>t<m>` transcript (c.) across edit types
// =============================================================================

mod lrg_transcript_cds {
    use super::*;

    #[test]
    fn sub_round_trips() {
        assert_round_trips("LRG_199t1:c.100A>G");
    }

    #[test]
    fn del_round_trips() {
        assert_round_trips("LRG_199t1:c.100_105del");
    }

    #[test]
    fn dup_round_trips() {
        assert_round_trips("LRG_199t1:c.100_105dup");
    }

    #[test]
    fn ins_round_trips() {
        assert_round_trips("LRG_199t1:c.100_101insAT");
    }

    #[test]
    fn inv_round_trips() {
        assert_round_trips("LRG_199t1:c.100_105inv");
    }

    #[test]
    fn delins_round_trips() {
        assert_round_trips("LRG_199t1:c.100_105delinsACGT");
    }

    #[test]
    fn identity_round_trips() {
        assert_round_trips("LRG_199t1:c.=");
    }

    #[test]
    fn unknown_round_trips() {
        assert_round_trips("LRG_199t1:c.?");
    }

    #[test]
    fn intronic_position_round_trips() {
        assert_round_trips("LRG_741t1:c.10078+1_10078+5delinsAAAAACCCTTGCAGAATGAAAA");
        assert_round_trips("LRG_741t1:c.1433-12_1433-10dup");
        assert_round_trips("LRG_741t1:c.7007+30del");
    }
}

// =============================================================================
// SECTION 3 — `LRG_<n>t<m>` transcript (n.) — non-coding RNA
// =============================================================================

mod lrg_transcript_noncoding {
    use super::*;

    #[test]
    fn sub_round_trips() {
        assert_round_trips("LRG_199t1:n.100A>G");
    }

    #[test]
    fn del_round_trips() {
        assert_round_trips("LRG_199t1:n.100_105del");
    }
}

// =============================================================================
// SECTION 4 — `LRG_<n>t<m>` transcript (r.) — RNA
// =============================================================================
//
// RNA uses lowercase bases.

mod lrg_transcript_rna {
    use super::*;

    #[test]
    fn sub_round_trips() {
        assert_round_trips("LRG_199t1:r.100a>g");
    }

    #[test]
    fn del_round_trips() {
        assert_round_trips("LRG_199t1:r.100_105del");
    }

    /// `r.0` is the spec's RNA-no-product form (no transcript produced).
    /// Distinct from the predicted form `r.(0)` (tracked elsewhere).
    #[test]
    fn rna_no_product_round_trips() {
        assert_round_trips("LRG_199t1:r.0");
    }
}

// =============================================================================
// SECTION 5 — `LRG_<n>p<m>` protein (p.) across edit types
// =============================================================================

mod lrg_protein {
    use super::*;

    #[test]
    fn sub_round_trips() {
        assert_round_trips("LRG_199p1:p.Arg8Gln");
    }

    /// Protein-level predicted wrapper has always worked.
    #[test]
    fn predicted_sub_round_trips() {
        assert_round_trips("LRG_199p1:p.(Arg8Gln)");
    }

    #[test]
    fn single_residue_del_round_trips() {
        assert_round_trips("LRG_199p1:p.Arg8del");
    }

    #[test]
    fn range_del_round_trips() {
        assert_round_trips("LRG_199p1:p.Arg8_Lys10del");
    }

    #[test]
    fn delins_round_trips() {
        assert_round_trips("LRG_199p1:p.Arg8_Lys10delinsAlaThr");
    }

    #[test]
    fn dup_round_trips() {
        assert_round_trips("LRG_199p1:p.Arg8_Lys10dup");
    }

    #[test]
    fn extension_round_trips() {
        assert_round_trips("LRG_199p1:p.Met1ext-5");
    }

    #[test]
    fn met1_unknown_round_trips() {
        assert_round_trips("LRG_199p1:p.Met1?");
    }

    #[test]
    fn frameshift_round_trips() {
        assert_round_trips("LRG_199p1:p.Arg8fs");
    }

    #[test]
    fn identity_round_trips() {
        assert_round_trips("LRG_199p1:p.=");
    }

    #[test]
    fn unknown_round_trips() {
        assert_round_trips("LRG_199p1:p.?");
    }

    /// `p.0` — no protein produced; explicitly endorsed for LRG by spec.
    #[test]
    fn no_protein_round_trips() {
        assert_round_trips("LRG_199p1:p.0");
    }

    /// `p.0?` — predicted no protein.
    #[test]
    fn no_protein_predicted_round_trips() {
        assert_round_trips("LRG_199p1:p.0?");
    }
}

// =============================================================================
// SECTION 6 — Gene-symbol selector composition with each LRG flavor
// =============================================================================

mod lrg_with_gene_symbol {
    use super::*;

    #[test]
    fn genomic_with_gene_round_trips() {
        assert_round_trips("LRG_199(BRCA1):g.100A>G");
    }

    #[test]
    fn transcript_cds_with_gene_round_trips() {
        assert_round_trips("LRG_199t1(BRCA1):c.100A>G");
    }

    #[test]
    fn transcript_rna_with_gene_round_trips() {
        assert_round_trips("LRG_199t1(BRCA1):r.100a>g");
    }

    #[test]
    fn protein_with_gene_drops_selector_on_display() {
        // Per HGVS syntax.yaml 119–128 the protein-variant grammar has no
        // parenthesized selector position. The parser still accepts the
        // `LRG_*p*(GENE):p.` form and stores the selector text in
        // `gene_symbol`, but Display emits the spec-compliant form without
        // the selector. See #310.
        let v = parse_hgvs("LRG_199p1(BRCA1):p.Arg8Gln")
            .unwrap_or_else(|e| panic!("parse failed: {e}"));
        // Parser remains permissive: gene selector text is retained on the
        // struct even though Display drops it.
        assert_eq!(v.gene_symbol(), Some("BRCA1"));
        assert_eq!(format!("{}", v), "LRG_199p1:p.Arg8Gln");
    }
}

// =============================================================================
// SECTION 7 — Compound alleles with LRG
// =============================================================================

mod lrg_compound {
    use super::*;

    #[test]
    fn cis_two_subs_round_trips() {
        assert_round_trips("LRG_199t1:c.[100A>G;105C>T]");
    }

    #[test]
    fn trans_two_subs_round_trips() {
        assert_round_trips("LRG_199t1:c.[100A>G];[105C>T]");
    }

    /// Two different LRG transcripts in trans — pins that the parser
    /// treats `t1` and `t2` as distinct accessions.
    #[test]
    fn trans_two_lrg_transcripts_round_trips() {
        assert_round_trips("[LRG_199t1:c.100A>G];[LRG_199t2:c.50C>T]");
    }

    /// Mixed-accession trans with one LRG and one RefSeq. Pins that the
    /// allele path accepts cross-flavor accessions on either side.
    #[test]
    fn trans_lrg_and_refseq_round_trips() {
        assert_round_trips("[LRG_199t1:c.1A>G];[NM_000088.3:c.2A>G]");
    }
}
