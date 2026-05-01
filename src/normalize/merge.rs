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
///
/// Variants are pushed to `output` immediately on arrival. While a chain is
/// growing, only the running anchor is updated (alt extended in place). At
/// chain end the head of `output` is rebuilt once from the merged anchor.
/// This keeps non-merging iterations as cheap as the no-merge baseline and
/// makes a chain of N consecutive merges O(N) total instead of O(N^2).
pub(crate) fn merge_consecutive_edits(
    variants: Vec<HgvsVariant>,
    phase: AllelePhase,
) -> Vec<HgvsVariant> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return variants;
    }

    let mut output: Vec<HgvsVariant> = Vec::with_capacity(variants.len());
    // Anchor of `output.last()` while the head is merge-eligible. None once
    // the head is non-eligible (a Duplication, uncertain edit, etc.) or
    // output is empty.
    let mut head_anchor: Option<Anchor> = None;
    // Whether at least one merge has folded into `head_anchor`. While true,
    // `output.last()` is stale relative to `head_anchor` and must be
    // rebuilt before any external observation.
    let mut head_merged = false;

    for next in variants {
        let merged = 'try_merge: {
            let Some(prev_a) = head_anchor.as_mut() else {
                break 'try_merge false;
            };
            let Some((next_start, _)) = simple_range_for_variant(&next) else {
                break 'try_merge false;
            };
            if prev_a.end.checked_add(1) != Some(next_start) {
                break 'try_merge false;
            }
            if !same_accession_and_kind(output.last().unwrap(), &next) {
                break 'try_merge false;
            }
            let Some(next_anchor) = anchor_for_variant(&next) else {
                break 'try_merge false;
            };
            // Extend in place — amortized O(1) per base added.
            prev_a.alt.extend(next_anchor.alt);
            prev_a.start = prev_a.start.min(next_anchor.start);
            prev_a.end = prev_a.end.max(next_anchor.end);
            true
        };
        if merged {
            head_merged = true;
            continue;
        }
        // Chain ends here. If we'd been growing one, rebuild the head once
        // from the final merged anchor before moving on.
        if head_merged {
            reconcile_head(&mut output, head_anchor.take().unwrap());
            head_merged = false;
        }
        head_anchor = anchor_for_variant(&next);
        output.push(next);
    }
    if head_merged {
        reconcile_head(&mut output, head_anchor.take().unwrap());
    }
    output
}

/// Replace `output.last()` with a freshly-built variant from `anchor`.
/// Caller has established that the head is merge-eligible (so kind dispatch
/// in `build_merged` is safe).
fn reconcile_head(output: &mut [HgvsVariant], anchor: Anchor) {
    let last = output.last_mut().expect("head must exist when head_merged");
    *last = build_merged(last, anchor);
}

/// Same accession and same `HgvsVariant` discriminant.
fn same_accession_and_kind(a: &HgvsVariant, b: &HgvsVariant) -> bool {
    match (a, b) {
        (HgvsVariant::Genome(av), HgvsVariant::Genome(bv)) => av.accession == bv.accession,
        (HgvsVariant::Cds(av), HgvsVariant::Cds(bv)) => av.accession == bv.accession,
        (HgvsVariant::Tx(av), HgvsVariant::Tx(bv)) => av.accession == bv.accession,
        (HgvsVariant::Rna(av), HgvsVariant::Rna(bv)) => av.accession == bv.accession,
        (HgvsVariant::Mt(av), HgvsVariant::Mt(bv)) => av.accession == bv.accession,
        _ => false,
    }
}

/// Position-only extraction for the cheap adjacency precheck. Returns the
/// anchor's `(start, end)` tuple if the variant is merge-eligible by type
/// and position, without allocating alt bases.
fn simple_range_for_variant(v: &HgvsVariant) -> Option<(u64, u64)> {
    match v {
        HgvsVariant::Genome(g) => simple_range_for_loc_edit(&g.loc_edit, simple_genome_range),
        HgvsVariant::Cds(c) => simple_range_for_loc_edit(&c.loc_edit, simple_cds_range),
        HgvsVariant::Tx(t) => simple_range_for_loc_edit(&t.loc_edit, simple_tx_range),
        HgvsVariant::Rna(r) => simple_range_for_loc_edit(&r.loc_edit, simple_rna_range),
        HgvsVariant::Mt(m) => simple_range_for_loc_edit(&m.loc_edit, simple_genome_range),
        _ => None,
    }
}

/// Per-coordinate-system dispatch for full anchor extraction.
fn anchor_for_variant(v: &HgvsVariant) -> Option<Anchor> {
    match v {
        HgvsVariant::Genome(g) => anchor_from_loc_edit(&g.loc_edit, simple_genome_range),
        HgvsVariant::Cds(c) => anchor_from_loc_edit(&c.loc_edit, simple_cds_range),
        HgvsVariant::Tx(t) => anchor_from_loc_edit(&t.loc_edit, simple_tx_range),
        HgvsVariant::Rna(r) => anchor_from_loc_edit(&r.loc_edit, simple_rna_range),
        HgvsVariant::Mt(m) => anchor_from_loc_edit(&m.loc_edit, simple_genome_range),
        _ => None,
    }
}

/// Position-only counterpart of `anchor_from_loc_edit`. Mirrors its
/// edit-kind and edit-certainty filters but does not touch alt bases.
fn simple_range_for_loc_edit<L>(
    loc_edit: &LocEdit<Interval<L>, NaEdit>,
    range_fn: impl Fn(&Interval<L>) -> Option<(u64, u64)>,
) -> Option<(u64, u64)> {
    if !loc_edit.edit.is_certain() {
        return None;
    }
    let edit = loc_edit.edit.inner()?;
    let (start, end) = range_fn(&loc_edit.location)?;
    match edit {
        NaEdit::Substitution { .. }
        | NaEdit::SubstitutionNoRef { .. }
        | NaEdit::Deletion { .. } => Some((start, end)),
        NaEdit::Delins { sequence } => {
            sequence.bases()?;
            Some((start, end))
        }
        NaEdit::Insertion { sequence } => {
            sequence.bases()?;
            // Anchor for an insertion is [end, start] — empty range at boundary.
            (end == start.checked_add(1)?).then_some((end, start))
        }
        _ => None,
    }
}

/// Build the merged variant from the head of output and the merged anchor.
/// Caller has already established that `head` and the partner are the same
/// kind via `same_accession_and_kind`. Takes `merged` by value so the alt
/// vec moves into the new variant rather than being cloned.
fn build_merged(head: &HgvsVariant, merged: Anchor) -> HgvsVariant {
    match head {
        HgvsVariant::Genome(g) => HgvsVariant::Genome(build_genome_merged(g, merged)),
        HgvsVariant::Cds(c) => HgvsVariant::Cds(build_cds_merged(c, merged)),
        HgvsVariant::Tx(t) => HgvsVariant::Tx(build_tx_merged(t, merged)),
        HgvsVariant::Rna(r) => HgvsVariant::Rna(build_rna_merged(r, merged)),
        HgvsVariant::Mt(m) => HgvsVariant::Mt(build_mt_merged(m, merged)),
        _ => unreachable!("build_merged called with non-NaEdit variant kind"),
    }
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
