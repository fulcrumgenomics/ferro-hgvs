//! Issue #182: post-canonicalization adjacency grouping of substitution runs
//! produced by `decompose_delins` (as corrected by issue #1034).
//!
//! `decompose_delins` splits a delins only where a *maximal contiguous run* is
//! entirely an inversion (`inv`) or where an interior identity separates
//! independent variants (`general.md:34`). A reverse-complement **sub-run** of
//! a longer contiguous change is never carved out — that change stays a single
//! `delins` (issue #1034). So a run of two or more strictly-adjacent
//! `Substitution` sub-edits only ever arises between an identity boundary and
//! another edit (an `inv`, an identity, or the span edge).
//!
//! When it does arise, those adjacent substitutions represent a
//! multi-nucleotide change and must render as a single `delins` per HGVS
//! (`substitution.md`: "changes involving two or more consecutive nucleotides
//! are described as deletion/insertion (delins)"). The grouping is local to
//! `build_split_variants` (in `src/normalize/mod.rs`). `Inversion` sub-edits
//! remain a hard barrier — sub flanks on either side of an inv stay split,
//! preserving the sub > del > **inv** > dup > ins priority (`general.md:56`).
//! `IdentityAt` also breaks a run and drops cleanly (no `=` emission).

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};

const SEQID: &str = "NC_TEST.1";

// ---------------------------------------------------------------------------
// Two adjacent sub flanks AFTER a full-run inv (identity between) → group
// ---------------------------------------------------------------------------

#[test]
fn two_subs_trailing_an_inv_merge_to_delins() {
    // Core "GATCC" at 257-261:
    //   257 G  258 A : full-run inv (revcomp(GA)=TC)
    //   259 T        : identity (bounds the run; dropped, no `259=`)
    //   260 C  261 C : two strictly-adjacent subs (CC→AA; revcomp(CC)=GG, so
    //                  NOT an inv) → group into one delins.
    // User-typed delins carries the interior identity so the decomposition
    // reaches build_split_variants.
    let p = SyntheticBuilder::genomic("GATCC").build();
    let result = normalize_to_string(p, &format!("{}:g.257_261delinsTCTAA", SEQID));
    assert_eq!(result, format!("{}:g.[257_258inv;260_261delinsAA]", SEQID));
}

// ---------------------------------------------------------------------------
// Two adjacent sub flanks BEFORE a full-run inv (identity between) → group
// ---------------------------------------------------------------------------

#[test]
fn two_subs_leading_an_inv_merge_to_delins() {
    // Core "CCTGA" at 257-261:
    //   257 C  258 C : two adjacent subs (CC→AA) → group into one delins
    //   259 T        : identity (dropped)
    //   260 G  261 A : full-run inv (revcomp(GA)=TC)
    let p = SyntheticBuilder::genomic("CCTGA").build();
    let result = normalize_to_string(p, &format!("{}:g.257_261delinsAATTC", SEQID));
    assert_eq!(result, format!("{}:g.[257_258delinsAA;260_261inv]", SEQID));
}

// ---------------------------------------------------------------------------
// Inv flanked on BOTH sides by two-sub runs (identities between) → both group
// ---------------------------------------------------------------------------

#[test]
fn inv_flanked_by_two_sub_runs_each_side_groups_both_flanks() {
    // Core "CCTGATCC" at 257-264:
    //   257 C  258 C : leading two-sub run (CC→AA) → delins
    //   259 T        : identity (dropped)
    //   260 G  261 A : full-run inv (revcomp(GA)=TC)
    //   262 T        : identity (dropped)
    //   263 C  264 C : trailing two-sub run (CC→AA) → delins
    // The inv stays between the two grouped flanks (not re-merged across it).
    let p = SyntheticBuilder::genomic("CCTGATCC").build();
    let result = normalize_to_string(p, &format!("{}:g.257_264delinsAATTCTAA", SEQID));
    assert_eq!(
        result,
        format!("{}:g.[257_258delinsAA;260_261inv;263_264delinsAA]", SEQID)
    );
}

// ---------------------------------------------------------------------------
// Single-sub flank (length-1 run) stays as Substitution
// ---------------------------------------------------------------------------

#[test]
fn single_sub_flank_stays_as_substitution() {
    // Core "CTGA" at 257-260:
    //   257 C        : length-1 sub run (C→A) → stays a substitution, NOT a
    //                  delins (the spec rule requires "two or more consecutive
    //                  nucleotides")
    //   258 T        : identity (dropped)
    //   259 G  260 A : full-run inv (revcomp(GA)=TC)
    let p = SyntheticBuilder::genomic("CTGA").build();
    let result = normalize_to_string(p, &format!("{}:g.257_260delinsATTC", SEQID));
    assert_eq!(result, format!("{}:g.[257C>A;259_260inv]", SEQID));
}

// ---------------------------------------------------------------------------
// IdentityAt breaks a substitution run (subs across an identity do not group)
// ---------------------------------------------------------------------------

#[test]
fn identity_at_between_subs_prevents_grouping() {
    // Core "GATCTC" at 257-262:
    //   257 G  258 A : full-run inv (revcomp(GA)=TC)
    //   259 T        : identity (dropped; bounds the inv run)
    //   260 C        : sub (C→A), length-1 run
    //   261 T        : identity (dropped; breaks the run)
    //   262 C        : sub (C→A), length-1 run
    // The subs at 260 and 262 are separated by the identity at 261, so they
    // must NOT group into a delins; each stays a singleton substitution and
    // no spurious `259=` / `261=` is emitted.
    let p = SyntheticBuilder::genomic("GATCTC").build();
    let result = normalize_to_string(p, &format!("{}:g.257_262delinsTCTATA", SEQID));
    assert_eq!(result, format!("{}:g.[257_258inv;260C>A;262C>A]", SEQID));
}
