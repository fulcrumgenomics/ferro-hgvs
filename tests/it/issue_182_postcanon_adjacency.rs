//! Issue #182: post-canonicalization adjacency of substitution flanks
//! produced by the inv-split pass (#160 / #166).
//!
//! When `decompose_delins` splits a merged delins into `[…; inv; …]`,
//! the non-inv positions are emitted as per-position `Substitution`
//! sub-edits. Two or more strictly-adjacent (no gap) Substitutions in that
//! emission stream represent a multi-nucleotide substitution and must
//! render as a single `delins` per HGVS spec — `substitution.md`: "changes
//! involving two or more consecutive nucleotides are described as
//! deletion/insertion (delins)" / `delins.md`: same rule.
//!
//! The fix is local to `build_split_variants` (in `src/normalize/mod.rs`):
//! it now groups runs of strictly-adjacent `Substitution` sub-edits into a
//! single `Delins` variant. `Inversion` sub-edits remain a hard barrier —
//! sub flanks on either side of an inv stay split, preserving the
//! sub > del > **inv** > dup > ins priority rule (`general.md:45`) and the
//! mutalyzer-aligned position adopted in #166. `IdentityAt` also breaks a
//! run (a position-gap means the flanking subs are no longer "consecutive"
//! per the spec's wording).
//!
//! Scope (deliberately narrow): this fix does NOT re-merge across an `inv`
//! to re-form a flat delins. The broader "should adjacency to inv force
//! delins?" question is tracked separately, gated on what SVD-WG010
//! settles — see the pushback comment on issue #182 for rationale.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};

const SEQID: &str = "NC_TEST.1";

// ---------------------------------------------------------------------------
// Two adjacent sub flanks AFTER an inv → group into one delins
// ---------------------------------------------------------------------------

#[test]
fn two_subs_trailing_an_inv_merge_to_delins() {
    // Core "GACTTTTGAAG" places:
    //   pos 257 G  258 A  259 C  (first sub run, no inv)
    //   pos 260-263 TTTT (gap of 4 — barrier between the two runs)
    //   pos 264 G  265 A  266 A  267 G  (second run, with 264_265 inv-eligible)
    //
    // Input merges to two delins after #80:
    //   257_259delinsCCT, 264_267delinsTCCC
    //
    // 264_267delins: ref GAAG, alt TCCC.
    //   - 264_265 ref GA / alt TC: revcomp(GA)=TC ✓ INV.
    //   - 266 A → C : SUB.
    //   - 267 G → C : SUB.
    // The two trailing subs are strictly adjacent (266→267, gap=0) and
    // must render as a single delinsCC per spec. The inv stays between
    // the leading delins and the new trailing delins.
    let p = SyntheticBuilder::genomic("GACTTTTGAAG").build();
    let result = normalize_to_string(
        p,
        &format!(
            "{}:g.[257G>C;258A>C;259C>T;264G>T;265A>C;266A>C;267G>C]",
            SEQID
        ),
    );
    assert_eq!(
        result,
        format!("{}:g.[257_259delinsCCT;264_265inv;266_267delinsCC]", SEQID)
    );
}

// ---------------------------------------------------------------------------
// Two adjacent sub flanks BEFORE an inv → group into one delins
// ---------------------------------------------------------------------------

#[test]
fn two_subs_leading_an_inv_merge_to_delins() {
    // Core "TCGACG" places ref T C G A C G at 257-262.
    //   257 T>A  258 C>G  : two-base sub flank to the left of the inv
    //   259 G>T  260 A>C  : 259_260 inv (revcomp(GA)=TC)
    //   261 C>X  262 G>X  : (omitted to keep this test focused on the leading flank)
    //
    // Inputs cover only 257-260 so the trailing positions are unedited.
    // Step #80 merges to 257_260delinsAGTC. Decompose:
    //   - 257_258 TC→AG: revcomp(TC)=GA, no.
    //   - 259_260 GA→TC: revcomp(GA)=TC ✓ INV.
    //   Sub(257 T→A); Sub(258 C→G); Inv(259_260).
    // The two leading subs (257, 258) are strictly adjacent → delinsAG.
    let p = SyntheticBuilder::genomic("TCGACG").build();
    let result = normalize_to_string(p, &format!("{}:g.[257T>A;258C>G;259G>T;260A>C]", SEQID));
    assert_eq!(result, format!("{}:g.[257_258delinsAG;259_260inv]", SEQID));
}

// ---------------------------------------------------------------------------
// Inv flanked on BOTH sides by two-sub runs → both flanks group, inv stays
// ---------------------------------------------------------------------------

#[test]
fn inv_flanked_by_two_sub_runs_each_side_groups_both_flanks() {
    // Core "TCGACG" places ref T C G A C G at 257-262.
    //   257 T>A  258 C>G  : two-base leading flank
    //   259 G>T  260 A>C  : 259_260 inv (revcomp(GA)=TC)
    //   261 C>A  262 G>T  : two-base trailing flank
    //
    // Step #80 merges to 257_262delinsAGTCAT. Decompose:
    //   - 259_260 GA→TC ✓ INV (longest, picked greedily at i=259).
    //   - i=257: no inv → Sub(257 T→A).
    //   - i=258: no inv → Sub(258 C→G).
    //   - i=261: no inv → Sub(261 C→A).
    //   - i=262: out of range for inv (i+2 > n) → Sub(262 G→T).
    // Both flanks are length-2 strictly-adjacent runs → each becomes a delins.
    // The inv stays between them (not re-merged across — see #166 / general.md:45).
    let p = SyntheticBuilder::genomic("TCGACG").build();
    let result = normalize_to_string(
        p,
        &format!("{}:g.[257T>A;258C>G;259G>T;260A>C;261C>A;262G>T]", SEQID),
    );
    assert_eq!(
        result,
        format!("{}:g.[257_258delinsAG;259_260inv;261_262delinsAT]", SEQID)
    );
}

// ---------------------------------------------------------------------------
// Single-sub flank (length-1 run) stays as Substitution
// ---------------------------------------------------------------------------

#[test]
fn single_sub_flank_stays_as_substitution() {
    // Core "TGAGC" at 257-261:
    //   257 T  258 G  259 A  260 G  261 C
    //
    // Inputs cover 257-260 only. Step #80 merges to 257_260delinsATCC.
    // Decompose: Sub(257 T→A); Inv(258_259); Sub(260 G→C). Each flank is a
    // length-1 run → emits as Substitution, not delins (the spec rule
    // requires "two or more consecutive nucleotides").
    //
    // 258_259 is non-palindromic (GA, revcomp TC). Length-2 palindromic
    // candidates (CG, AT, GC, TA) collapse under `shorten_inversion` and
    // are not emitted as Inversion sub-edits.
    let p = SyntheticBuilder::genomic("TGAGC").build();
    let result = normalize_to_string(p, &format!("{}:g.[257T>A;258G>T;259A>C;260G>C]", SEQID));
    assert_eq!(result, format!("{}:g.[257T>A;258_259inv;260G>C]", SEQID));
}

// ---------------------------------------------------------------------------
// IdentityAt breaks a substitution run
// ---------------------------------------------------------------------------

#[test]
fn identity_at_in_decomposition_breaks_sub_run() {
    // Core "TGAGCC" at 257-262 produces a delins where the decomposer
    // emits an IdentityAt between two Substitution sub-edits. The IdentityAt
    // both (a) drops cleanly (no `g.<pos>=` emission) and (b) breaks the
    // sub run so the bracketing subs do NOT group into a delins.
    //
    //   pos 257 T  258 G  259 A  260 G  261 C  262 C
    //
    // Inputs:
    //   257 T>A     : leading sub flank (length 1)
    //   258 G>T     : start of inv 258_259 (alt T)
    //   259 A>C     : end   of inv 258_259 (alt C; revcomp(GA)=TC ✓)
    //   260 G>G     : alt = ref → IdentityAt at 260
    //   261 C>A     : trailing sub at 261
    //   262 C>T     : (omitted — keeping the trailing run length 1 isolates
    //                  the identity-breaks-run case from the run-grouping case)
    //
    // Expected merged delins: 257_261delinsATCGA.
    // Decompose: Sub(257 T→A); Inv(258_259); IdentityAt(260); Sub(261 C→A).
    // Output: [257T>A; 258_259inv; 261C>A] — IdentityAt drops, no spurious
    // `260=` emission, and 261 stays a singleton Substitution (no run to
    // group with — it's bracketed by the inv on the left and end-of-allele
    // on the right).
    let p = SyntheticBuilder::genomic("TGAGCC").build();
    let result = normalize_to_string(p, &format!("{}:g.[257T>A;258G>T;259A>C;261C>A]", SEQID));
    assert_eq!(result, format!("{}:g.[257T>A;258_259inv;261C>A]", SEQID));
}
