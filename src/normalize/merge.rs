//! Merge consecutive sub-variants in a cis allele into a single edit per HGVS spec.
//!
//! See `docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md`.

use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::GenomePos;
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{AllelePhase, GenomeVariant, HgvsVariant, LocEdit};

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
    variants: &[HgvsVariant],
    phase: AllelePhase,
) -> Vec<HgvsVariant> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return variants.to_vec();
    }

    let mut output: Vec<HgvsVariant> = Vec::with_capacity(variants.len());
    for next in variants {
        if let Some(prev) = output.last() {
            if let Some(merged) = try_merge_pair(prev, next) {
                let last = output.last_mut().unwrap();
                *last = merged;
                continue;
            }
        }
        output.push(next.clone());
    }
    output
}

fn try_merge_pair(a: &HgvsVariant, b: &HgvsVariant) -> Option<HgvsVariant> {
    // Phase 1: only Genome-Genome same-accession merging.
    let (av, bv) = match (a, b) {
        (HgvsVariant::Genome(av), HgvsVariant::Genome(bv)) => (av, bv),
        _ => return None,
    };
    if av.accession != bv.accession {
        return None;
    }
    let a_anchor = genome_anchor(av)?;
    let b_anchor = genome_anchor(bv)?;
    if a_anchor.end + 1 != b_anchor.start {
        return None;
    }
    let mut alt = a_anchor.alt;
    alt.extend(b_anchor.alt);
    Some(HgvsVariant::Genome(build_genome_merged(
        av,
        a_anchor.start,
        b_anchor.end,
        alt,
    )))
}

fn genome_anchor(v: &GenomeVariant) -> Option<Anchor> {
    let edit = v.loc_edit.edit.inner()?;
    if !v.loc_edit.edit.is_certain() {
        return None;
    }
    let (start, end) = simple_genome_range(&v.loc_edit.location)?;
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
        _ => None,
    }
}

fn simple_genome_range(interval: &Interval<GenomePos>) -> Option<(u64, u64)> {
    let start = simple_genome_pos(&interval.start)?;
    let end = simple_genome_pos(&interval.end)?;
    Some((start, end))
}

fn simple_genome_pos(boundary: &UncertainBoundary<GenomePos>) -> Option<u64> {
    let mu = boundary.as_single()?;
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_special() || pos.offset.is_some() {
        return None;
    }
    Some(pos.base)
}

fn build_genome_merged(
    template: &GenomeVariant,
    start: u64,
    end: u64,
    alt: Vec<Base>,
) -> GenomeVariant {
    let location = Interval::new(GenomePos::new(start), GenomePos::new(end));
    let edit = if alt.is_empty() {
        // Spec-recommended no-sequence/no-length form.
        NaEdit::Deletion {
            sequence: None,
            length: None,
        }
    } else {
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::new(alt)),
        }
    };
    GenomeVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::variant::Accession;

    #[test]
    fn empty_input_returns_empty() {
        assert!(merge_consecutive_edits(&[], AllelePhase::Cis).is_empty());
    }

    #[test]
    fn single_input_passes_through() {
        // Single-variant alleles are not Phase-1 mergeable.
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
        let out = merge_consecutive_edits(std::slice::from_ref(&v), AllelePhase::Cis);
        assert_eq!(out.len(), 1);
    }
}
