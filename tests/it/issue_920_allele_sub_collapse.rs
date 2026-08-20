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

/// A group that straddles the CDS/3'UTR seam now derives its sequence-canonical
/// form instead of being refused (#1816, combined with #2174's contiguous-run
/// collapse).
///
/// `core = "GGGGGGGACGCC"`, `cds_start = 1`, `cds_end = 9`: `c.1..c.9` are the
/// CDS body and `c.*1` is the first 3'UTR base. `c.[7_8insC;*1del]` denotes
/// `GGGGGGGCACCC`, whose diff against the reference is **three consecutive
/// changed bases** at `c.8`/`c.9`/`c.*1` (`A→C`, `C→A`, `G→C`) — one spanning
/// `delins` by `delins.md:16`, not two members separated by an unchanged base.
///
/// This used to be refused ("all-or-nothing, because the 3'UTR member is outside
/// the positive CDS body"), which was the straddle-anchor limitation, not a
/// principle: `Anchor` carried one `Region`, so a member spanning `c.9_*1` had
/// no anchor. #1816 gives `Anchor` a region per endpoint, so the derivation
/// renders the spanning `delins` across the seam — the same sequence-first
/// canonical form (`canonical-form-choice-when-both-legal`) the CDS-interior
/// twins above collapse to. Meaning-preserving: both spellings denote
/// `GGGGGGGCACCC`. The result is a run of consecutive changes, so the
/// separation-two-or-more census is untouched.
#[test]
fn utr_straddle_collapses_to_its_spanning_delins() {
    let p = SyntheticBuilder::cds("GGGGGGGACGCC", 1, 9, Strand::Plus).build();
    let result = normalize_to_string(p, "NM_TEST.1:c.[7_8insC;*1del]");
    assert_eq!(result, "NM_TEST.1:c.8_*1delinsCAC");
}
