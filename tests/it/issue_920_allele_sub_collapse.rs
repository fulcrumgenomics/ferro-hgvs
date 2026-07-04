//! Issue #920: collapse an overlapping cis del+ins on the transcript axis
//! into the equivalent substitution.
//!
//! `collapse_overlapping_cis_edits` already folds an overlapping del+ins cis
//! group into a single `delins` on the genomic (`g.`/`m.`) axis, and a
//! single-base `delins` reduces to a substitution downstream. Before #920 the
//! pass bailed for any non-`Genome`/`Mt` variant, so a transcript-axis group
//! like `NG_012337.1(NM_012459.2):c.[8del;7_8insC]` was left as a two-member
//! allele instead of collapsing to `c.8A>C`.
//!
//! #920 generalizes the collapse to the transcript axes (`c.`/`n.`/`r.`),
//! scoped to the **positive body** only (`Cds` for `c.`, `Tx` for `n.`,
//! `Rna` for `r.`). Members in the 5'UTR (`c.-N`), 3'UTR (`c.*N`), an
//! intronic offset, upstream/downstream, or a mixed-region group refuse the
//! whole collapse (all-or-nothing), leaving the allele unchanged.
//!
//! The collapse operates purely in transcript coordinates (it fetches the
//! transcript sequence directly), so it is strand-agnostic; the minus-strand
//! case below exercises that.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
use ferro_hgvs::reference::transcript::Strand;

/// Plus-strand coding transcript with a non-trivial `cds_start` (= 4), so the
/// collapse's CDS→transcript coordinate translation (delta = `cds_start - 1`)
/// is genuinely exercised rather than an identity map.
///
/// Layout of `core = "TTTGGGGGGGACTTTTT"` (1-based transcript pos → base):
///   pos 1..3   = 5'UTR (`c.-3`,`c.-2`,`c.-1`)
///   pos 4..12  = CDS (`c.1`..`c.9`) = `GGGGGGGAC`, so `c.7 = G`, `c.8 = A`
///   pos 13..17 = 3'UTR
///
/// `c.[8del;7_8insC]`: delete `c.8` (A) and insert `C` between `c.7`/`c.8`.
/// Over the window `c.7..c.8` the variant is `G C` vs reference `G A`, so the
/// net edit is `c.8delinsC`, which reduces to the substitution `c.8A>C`.
#[test]
fn cds_body_overlapping_del_ins_collapses_to_sub_plus_strand() {
    let p = SyntheticBuilder::cds("TTTGGGGGGGACTTTTT", 4, 12, Strand::Plus).build();
    let result = normalize_to_string(p, "NM_TEST.1:c.[8del;7_8insC]");
    assert_eq!(result, "NM_TEST.1:c.8A>C");
}

/// Same collapse on a minus-strand transcript. The collapse works in
/// transcript (CDS-relative) coordinates, so strand does not change the ref
/// base read at `c.8`. `cds_start = 1` here (delta = 0).
///
/// `core = "GGGGGGGACGGGGGGG"`: `c.7 = G`, `c.8 = A`. `c.[8del;7_8insC]`
/// collapses to `c.8A>C`.
#[test]
fn cds_body_overlapping_del_ins_collapses_to_sub_minus_strand() {
    let p = SyntheticBuilder::cds("GGGGGGGACGGGGGGG", 1, 16, Strand::Minus).build();
    let result = normalize_to_string(p, "NM_TEST.1:c.[8del;7_8insC]");
    assert_eq!(result, "NM_TEST.1:c.8A>C");
}

/// The same collapse on the non-coding transcript (`n.`) axis, whose body is
/// `Region::Tx` and whose coordinate is the transcript sequence directly
/// (delta = 0). `core = "GGGGGGGACGGGGGGG"`: `n.7 = G`, `n.8 = A`.
#[test]
fn tx_body_overlapping_del_ins_collapses_to_sub() {
    let p = SyntheticBuilder::noncoding("GGGGGGGACGGGGGGG", Strand::Plus).build();
    let result = normalize_to_string(p, "NR_TEST.1:n.[8del;7_8insC]");
    assert_eq!(result, "NR_TEST.1:n.8A>C");
}

/// Refusal: a group containing a 3'UTR (`c.*N`) member is outside the positive
/// CDS body, so the whole collapse is refused (all-or-nothing) and the allele
/// is left as its separate members rather than collapsing to a substitution.
///
/// `core = "GGGGGGGACGCC"`, `cds_start = 1`, `cds_end = 9`: `c.1..c.9` are the
/// CDS body and `c.*1` is the first 3'UTR base. `c.[7_8insC;*1del]` mixes a
/// CDS-body insertion with a 3'UTR deletion — no collapse.
#[test]
fn utr_member_refuses_collapse() {
    let p = SyntheticBuilder::cds("GGGGGGGACGCC", 1, 9, Strand::Plus).build();
    let result = normalize_to_string(p, "NM_TEST.1:c.[7_8insC;*1del]");
    // The allele must be returned unchanged — the 3'UTR member forces the
    // whole collapse to be refused (all-or-nothing), so both members survive
    // verbatim rather than collapsing to a single substitution. Pin the exact
    // string so a regression that reorders, drops, or partially collapses a
    // member is caught (not just the presence of a `;`).
    assert_eq!(result, "NM_TEST.1:c.[7_8insC;*1del]");
}
