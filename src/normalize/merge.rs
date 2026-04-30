//! Merge consecutive sub-variants in a cis allele into a single edit per HGVS spec.
//!
//! See `docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md`.

use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{
    AllelePhase, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, MtVariant, RnaVariant, TxVariant,
};

/// Anchor for a single sub-variant.
///
/// `start` and `end` are 1-based inclusive position bounds. For an `Insertion`
/// between positions `p` and `p+1`, `start = p+1` and `end = p` (empty range
/// at a boundary).
#[derive(Debug, Clone)]
struct Anchor {
    start: u64,
    end: u64,
    alt: Vec<Base>,
}

/// Merge consecutive sub-variants in an allele.
///
/// Returns the input unchanged unless `phase == Cis` and at least one merge
/// is possible. Sub-variants that aren't merge-eligible (different accession,
/// non-NaEdit, uncertain edit, disqualifying position, etc.) act as merge
/// barriers — they pass through unchanged in their original input order.
pub(crate) fn merge_consecutive_edits(
    variants: Vec<HgvsVariant>,
    phase: AllelePhase,
) -> Vec<HgvsVariant> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return variants;
    }

    let mut output: Vec<HgvsVariant> = Vec::with_capacity(variants.len());
    for next in variants {
        if let Some(prev) = output.last() {
            if let Some(merged) = try_merge_pair(prev, &next) {
                let last = output.last_mut().unwrap();
                *last = merged;
                continue;
            }
        }
        output.push(next);
    }
    output
}

fn try_merge_pair(a: &HgvsVariant, b: &HgvsVariant) -> Option<HgvsVariant> {
    match (a, b) {
        (HgvsVariant::Genome(av), HgvsVariant::Genome(bv)) => {
            try_merge_genome(av, bv).map(HgvsVariant::Genome)
        }
        (HgvsVariant::Cds(av), HgvsVariant::Cds(bv)) => try_merge_cds(av, bv).map(HgvsVariant::Cds),
        (HgvsVariant::Tx(av), HgvsVariant::Tx(bv)) => try_merge_tx(av, bv).map(HgvsVariant::Tx),
        (HgvsVariant::Rna(av), HgvsVariant::Rna(bv)) => try_merge_rna(av, bv).map(HgvsVariant::Rna),
        (HgvsVariant::Mt(av), HgvsVariant::Mt(bv)) => try_merge_mt(av, bv).map(HgvsVariant::Mt),
        _ => None,
    }
}

fn try_merge_genome(a: &GenomeVariant, b: &GenomeVariant) -> Option<GenomeVariant> {
    if a.accession != b.accession {
        return None;
    }
    let aa = anchor_from_loc_edit(&a.loc_edit, simple_genome_range)?;
    let ba = anchor_from_loc_edit(&b.loc_edit, simple_genome_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_genome_merged(a, merged))
}

fn try_merge_cds(a: &CdsVariant, b: &CdsVariant) -> Option<CdsVariant> {
    if a.accession != b.accession {
        return None;
    }
    let aa = anchor_from_loc_edit(&a.loc_edit, simple_cds_range)?;
    let ba = anchor_from_loc_edit(&b.loc_edit, simple_cds_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_cds_merged(a, merged))
}

fn try_merge_tx(a: &TxVariant, b: &TxVariant) -> Option<TxVariant> {
    if a.accession != b.accession {
        return None;
    }
    let aa = anchor_from_loc_edit(&a.loc_edit, simple_tx_range)?;
    let ba = anchor_from_loc_edit(&b.loc_edit, simple_tx_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_tx_merged(a, merged))
}

fn try_merge_rna(a: &RnaVariant, b: &RnaVariant) -> Option<RnaVariant> {
    if a.accession != b.accession {
        return None;
    }
    let aa = anchor_from_loc_edit(&a.loc_edit, simple_rna_range)?;
    let ba = anchor_from_loc_edit(&b.loc_edit, simple_rna_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_rna_merged(a, merged))
}

fn try_merge_mt(a: &MtVariant, b: &MtVariant) -> Option<MtVariant> {
    if a.accession != b.accession {
        return None;
    }
    let aa = anchor_from_loc_edit(&a.loc_edit, simple_genome_range)?;
    let ba = anchor_from_loc_edit(&b.loc_edit, simple_genome_range)?;
    let merged = merge_anchors(aa, ba)?;
    Some(build_mt_merged(a, merged))
}

/// Combine two adjacency-checked anchors into one merged anchor.
/// Returns None if the pair is not strictly adjacent (`a.end + 1 == b.start`).
/// For two `Insertion` anchors at the same boundary `p|p+1` (both `start=p+1, end=p`),
/// `a.end + 1 = p + 1 == b.start = p + 1`, so adjacency holds.
fn merge_anchors(a: Anchor, b: Anchor) -> Option<Anchor> {
    if a.end.checked_add(1)? != b.start {
        return None;
    }
    let mut alt = a.alt;
    alt.extend(b.alt);
    Some(Anchor {
        start: a.start.min(b.start),
        end: a.end.max(b.end),
        alt,
    })
}

/// Extract an anchor from a sub-variant's location+edit. The `range_fn`
/// callback returns `(start, end)` for the location only when it is "simple"
/// (no offsets, no UTR markers, certain). Returns None when the edit is
/// uncertain, unknown, or not a merge-eligible NaEdit variant.
fn anchor_from_loc_edit<L>(
    loc_edit: &LocEdit<Interval<L>, NaEdit>,
    range_fn: impl Fn(&Interval<L>) -> Option<(u64, u64)>,
) -> Option<Anchor> {
    if !loc_edit.edit.is_certain() {
        return None;
    }
    let edit = loc_edit.edit.inner()?;
    let (start, end) = range_fn(&loc_edit.location)?;
    match edit {
        NaEdit::Substitution { alternative, .. } | NaEdit::SubstitutionNoRef { alternative } => {
            Some(Anchor {
                start,
                end,
                alt: vec![*alternative],
            })
        }
        NaEdit::Deletion { .. } => Some(Anchor {
            start,
            end,
            alt: Vec::new(),
        }),
        NaEdit::Delins { sequence } => {
            let bases = sequence.bases()?.to_vec();
            Some(Anchor {
                start,
                end,
                alt: bases,
            })
        }
        NaEdit::Insertion { sequence } => {
            if end != start.checked_add(1)? {
                return None;
            }
            let bases = sequence.bases()?.to_vec();
            Some(Anchor {
                start: end,
                end: start,
                alt: bases,
            })
        }
        _ => None,
    }
}

fn simple_genome_range(interval: &Interval<GenomePos>) -> Option<(u64, u64)> {
    Some((
        simple_genome_pos(&interval.start)?,
        simple_genome_pos(&interval.end)?,
    ))
}

fn simple_genome_pos(boundary: &UncertainBoundary<GenomePos>) -> Option<u64> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_special() || pos.offset.is_some() {
        return None;
    }
    Some(pos.base)
}

fn simple_cds_range(interval: &Interval<CdsPos>) -> Option<(u64, u64)> {
    Some((
        simple_cds_pos(&interval.start)?,
        simple_cds_pos(&interval.end)?,
    ))
}

fn simple_cds_pos(boundary: &UncertainBoundary<CdsPos>) -> Option<u64> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_unknown() || pos.is_intronic() || pos.is_5utr() || pos.is_3utr() || pos.base < 1 {
        return None;
    }
    Some(pos.base as u64)
}

fn simple_tx_range(interval: &Interval<TxPos>) -> Option<(u64, u64)> {
    Some((
        simple_tx_pos(&interval.start)?,
        simple_tx_pos(&interval.end)?,
    ))
}

fn simple_tx_pos(boundary: &UncertainBoundary<TxPos>) -> Option<u64> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_intronic() || pos.is_upstream() || pos.is_downstream() || pos.base < 1 {
        return None;
    }
    Some(pos.base as u64)
}

fn simple_rna_range(interval: &Interval<RnaPos>) -> Option<(u64, u64)> {
    Some((
        simple_rna_pos(&interval.start)?,
        simple_rna_pos(&interval.end)?,
    ))
}

fn simple_rna_pos(boundary: &UncertainBoundary<RnaPos>) -> Option<u64> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_intronic() || pos.is_3utr() || pos.is_5utr() || pos.base < 1 {
        return None;
    }
    Some(pos.base as u64)
}

fn build_genome_merged(template: &GenomeVariant, merged: Anchor) -> GenomeVariant {
    let (location, edit) = build_naedit(merged, GenomePos::new);
    GenomeVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_cds_merged(template: &CdsVariant, merged: Anchor) -> CdsVariant {
    let (location, edit) = build_naedit(merged, |b| CdsPos::new(b as i64));
    CdsVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_tx_merged(template: &TxVariant, merged: Anchor) -> TxVariant {
    let (location, edit) = build_naedit(merged, |b| TxPos::new(b as i64));
    TxVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_rna_merged(template: &RnaVariant, merged: Anchor) -> RnaVariant {
    let (location, edit) = build_naedit(merged, |b| RnaPos::new(b as i64));
    RnaVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_mt_merged(template: &MtVariant, merged: Anchor) -> MtVariant {
    let (location, edit) = build_naedit(merged, GenomePos::new);
    MtVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_naedit<P>(merged: Anchor, mut to_pos: impl FnMut(u64) -> P) -> (Interval<P>, NaEdit) {
    let edit = if merged.start > merged.end {
        debug_assert_eq!(
            merged.start,
            merged.end + 1,
            "invariant: insertion anchor span = 1 nt"
        );
        NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::new(merged.alt)),
        }
    } else if merged.alt.is_empty() {
        NaEdit::Deletion {
            sequence: None,
            length: None,
        }
    } else {
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(merged.alt)),
        }
    };
    let (lo, hi) = if merged.start > merged.end {
        (merged.end, merged.start)
    } else {
        (merged.start, merged.end)
    };
    (Interval::new(to_pos(lo), to_pos(hi)), edit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::variant::Accession;

    #[test]
    fn empty_input_returns_empty() {
        assert!(merge_consecutive_edits(vec![], AllelePhase::Cis).is_empty());
    }

    #[test]
    fn single_input_passes_through() {
        let acc = Accession::new("NC", "000001", Some(11));
        let v = HgvsVariant::Genome(GenomeVariant {
            accession: acc,
            gene_symbol: None,
            loc_edit: LocEdit::new(
                Interval::new(GenomePos::new(100), GenomePos::new(100)),
                NaEdit::Substitution {
                    reference: Base::G,
                    alternative: Base::A,
                },
            ),
        });
        let out = merge_consecutive_edits(vec![v], AllelePhase::Cis);
        assert_eq!(out.len(), 1);
    }
}
