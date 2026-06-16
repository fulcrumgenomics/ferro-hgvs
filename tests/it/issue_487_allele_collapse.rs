//! Issue #487: collapse overlapping cis-allele edits into a single canonical edit.
//!
//! HGVS `general.md` requires consecutive edits in cis to render as a single
//! edit. ferro's `merge_consecutive_edits` handles strictly-consecutive
//! same-region edits, but *overlapping* mixed edits at one locus (insertions
//! flanking a deletion) were left as separate sub-variants. Applying the cis
//! edits to the reference and re-deriving the minimal edit collapses them.
//!
//! Mutalyzer-verified shape (corpus): `NG_012337.1:g.[104_105insA;
//! 105_106insC;105del]` -> `NG_012337.1:g.105delinsAC`.
//!
//! ## Implementation
//!
//! Implemented as `collapse_overlapping_cis_edits` in src/normalize/merge.rs,
//! invoked from `normalize_allele` before `merge_consecutive_edits`.
//!
//! `merge_consecutive_edits` (src/normalize/merge.rs) handles strictly-
//! consecutive *non-overlapping* same-region edits via an anchor chain
//! (`prev.end + 1 == next.start`). The remaining gap is *overlapping* cis
//! edits — insertions flanking a deletion/sub at one locus — handled by
//! `collapse_overlapping_cis_edits`: it builds a `GEdit` window over the
//! changed interval plus insertion flanks, applies each edit into the
//! window's `cell`/`after`/`before` buffers, then minimal-trims the result
//! against the reference window to emit a single `delins`. The pass is
//! all-or-nothing and refuses any group it cannot collapse unambiguously
//! (unsupported edit, non-contiguous changed interval, or two insertions
//! sharing a gap). It covers both the genomic (`g.`) and mitochondrial
//! (`m.`) single-axis coordinate systems.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};

const SEQID: &str = "NC_TEST.1";

/// Core `AACGT` -> HGVS 257=A, 258=A, 259=C, 260=G, 261=T.
/// At locus 259..261: insert `A` between C(259)/G(260), insert `C` between
/// G(260)/T(261), delete G(260). The variant sequence over 259..261 becomes
/// `C A C T` vs reference `C G T`, so the net edit is `g.260delinsAC`.
#[test]
fn overlapping_ins_ins_del_collapses_to_delins() {
    let p = SyntheticBuilder::genomic("AACGT").build();
    let result = normalize_to_string(p, &format!("{}:g.[259_260insA;260_261insC;260del]", SEQID));
    assert_eq!(result, format!("{}:g.260delinsAC", SEQID));
}

/// Reordered input (deletion first) collapses to the same canonical edit —
/// cis-allele member order must not affect the result.
#[test]
fn overlapping_edits_collapse_is_order_independent() {
    let p = SyntheticBuilder::genomic("AACGT").build();
    let result = normalize_to_string(p, &format!("{}:g.[260del;259_260insA;260_261insC]", SEQID));
    assert_eq!(result, format!("{}:g.260delinsAC", SEQID));
}

/// Mitochondrial (`m.`) cis edits share the same single-axis coordinate
/// system as `g.` (`Region::Genome`), so the overlap-collapse pass handles
/// them too. Same `AACGT` core as the genomic case: insert `A` between
/// C(259)/G(260), insert `C` between G(260)/T(261), delete G(260) — the net
/// edit over 259..261 is `m.260delinsAC`. (`NC_TEST.1:m.` parses to an
/// `MtVariant`; the synthetic genomic provider serves its sequence by contig.)
#[test]
fn overlapping_mito_ins_ins_del_collapses_to_delins() {
    let p = SyntheticBuilder::genomic("AACGT").build();
    let result = normalize_to_string(p, &format!("{}:m.[259_260insA;260_261insC;260del]", SEQID));
    assert_eq!(result, format!("{}:m.260delinsAC", SEQID));
}

/// Two insertions at the *same* gap would concatenate in member order, making
/// the collapsed result order-dependent. The pass must refuse to collapse such
/// a group rather than emit an arbitrary ordering — it passes through as the
/// original separate sub-variants. (`g.260del` 5'-shifts/renormalizes to
/// `g.260del`; the two same-gap inserts stay as authored.)
#[test]
fn same_gap_insertions_are_not_collapsed() {
    let p = SyntheticBuilder::genomic("AACGT").build();
    let result = normalize_to_string(p, &format!("{}:g.[260_261insA;260_261insC;260del]", SEQID));
    // No single-delins collapse: the result keeps the inserts as distinct
    // members, so it is not the collapsed `g.260delins...` form.
    assert!(
        !result.contains("delins"),
        "same-gap insertions must not collapse to a single delins, got: {result}"
    );
}
