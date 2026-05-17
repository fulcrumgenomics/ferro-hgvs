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
use crate::reference::ReferenceProvider;

/// Coordinate-system region used as the merge-eligibility key.
///
/// Adjacency in the merge pass requires both ends to share a region.
/// HGVS forbids ranges that span the 5'UTR↔CDS, CDS↔3'UTR (for `c.`/`r.`)
/// or upstream↔transcript / transcript↔downstream (for `n.`) boundaries
/// — `c.-1_1`, `c.40_*1`, etc. do not exist as valid range syntax — so
/// we tag each position with its region and refuse to merge across.
/// The integer `start`/`end` axis is region-local: positive for `Cds` /
/// `ThreePrimeUtr` / `Tx` / `TxDownstream`, negative for `FivePrimeUtr`
/// / `TxUpstream`. Adjacency `prev.end + 1 == next.start` works
/// naturally for all six because the axis is monotonic 5'→3' within
/// each region (`c.-3 → c.-2 → c.-1` maps to `-3 → -2 → -1`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum Region {
    /// `g.` (genomic) and `m.` (mitochondrial) — single axis, no
    /// sub-regions.
    Genome,
    /// `c.` CDS proper (positive non-UTR base).
    Cds,
    /// `r.` CDS proper (positive non-UTR base on a transcript).
    Rna,
    /// `c.` / `r.` 5'UTR (negative base, e.g. `c.-3`).
    FivePrimeUtr,
    /// `c.` / `r.` 3'UTR (`utr3` flag, e.g. `c.*1`).
    ThreePrimeUtr,
    /// `n.` transcript body (positive non-downstream base).
    Tx,
    /// `n.` upstream of transcript (negative base, e.g. `n.-3`).
    TxUpstream,
    /// `n.` downstream of transcript (`downstream` flag, e.g. `n.*1`).
    TxDownstream,
}

/// Anchor for a single sub-variant.
///
/// `start` and `end` are 1-based inclusive position bounds on the
/// `region` axis. For an `Insertion` between positions `p` and `p+1`,
/// `start = p+1` and `end = p` (empty range at a boundary). The
/// integer axis is signed because 5'UTR / upstream regions use
/// negative values.
#[derive(Debug, Clone)]
struct Anchor {
    region: Region,
    start: i64,
    end: i64,
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
pub(crate) fn merge_consecutive_edits<P: ReferenceProvider>(
    variants: Vec<HgvsVariant>,
    phase: AllelePhase,
    provider: &P,
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
            let Some((next_region, next_start, _)) = simple_range_for_variant(&next) else {
                break 'try_merge false;
            };
            // Region match is required for both adjacency rules below.
            if prev_a.region != next_region {
                break 'try_merge false;
            }

            // Strictly-consecutive adjacency: `prev.end + 1 == next.start`
            // on the region's signed integer axis. Works uniformly for
            // every region (5'UTR / upstream use negative axis values, so
            // `c.-3 + 1 == c.-2` is `-3 + 1 == -2`).
            let strict_adjacent = prev_a.end.checked_add(1) == Some(next_start);

            // Codon-frame exception (issue #79, extended by issue #275):
            // two exonic edits on the CDS axis (`c.` or `r.`), separated
            // by exactly one nucleotide, that fall within the same codon
            // merge into a delins with the unchanged middle reference base
            // preserved verbatim. Eligibility:
            //   * same accession + same kind (checked below).
            //   * region is `Cds` (for `c.`) or `Rna` (for `r.`) — both
            //     share a codon-relative axis. `g.` / `n.` / UTR / intron
            //     are excluded.
            //   * `prev_a.end + 2 == next_start` (gap of exactly 1 nt).
            //   * `same_codon(prev_a.end, next_start)` — both endpoints
            //     fall in the same 1-indexed codon. The right edge of
            //     `prev_a` is what matters, so a chain that has already
            //     grown via earlier strict merges can still extend via
            //     codon-frame (issue #275 item 2).
            //   * `prev_a.start <= prev_a.end` — `prev_a` is not an
            //     insertion anchor (those use `start == end + 1`).
            //   * `next_anchor` is a 1-position substitution OR deletion
            //     (checked below): `next_anchor.start == next_anchor.end`
            //     and `next_anchor.alt.len() <= 1`. Issue #275 item 3
            //     relaxes the original `alt.len() == 1` requirement to
            //     allow `sub`+`del` (and `del`+`sub` mirror) pairs.
            let codon_frame_eligible = !strict_adjacent
                && (prev_a.region == Region::Cds || prev_a.region == Region::Rna)
                && prev_a.end.checked_add(2) == Some(next_start)
                && prev_a.start <= prev_a.end
                && same_codon(prev_a.end, next_start);

            if !strict_adjacent && !codon_frame_eligible {
                break 'try_merge false;
            }
            if !same_accession_and_kind(output.last().unwrap(), &next) {
                break 'try_merge false;
            }
            let Some(next_anchor) = anchor_for_variant(&next) else {
                break 'try_merge false;
            };
            if codon_frame_eligible {
                // Next must be a 1-position substitution OR deletion
                // (issue #275 item 3 — `del` partners are eligible too).
                // Insertions and multi-position delins are excluded:
                // their position semantics ("between p and p+1" for ins,
                // ranged for delins) don't satisfy "one variant
                // separated by one nucleotide". A 1-position del has
                // `start == end` and `alt.len() == 0`; a 1-position sub
                // has `start == end` and `alt.len() == 1`.
                if next_anchor.start != next_anchor.end || next_anchor.alt.len() > 1 {
                    break 'try_merge false;
                }
                // Look up the unchanged middle reference base. If the
                // provider has no transcript or no sequence covering the
                // gap position, gracefully decline the codon-frame merge
                // rather than erroring.
                let Some(middle_ref) =
                    lookup_codon_middle_ref(provider, output.last().unwrap(), prev_a.end + 1)
                else {
                    break 'try_merge false;
                };
                prev_a.alt.push(middle_ref);
            }
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
/// anchor's `(region, start, end)` tuple if the variant is merge-eligible by
/// type and position, without allocating alt bases.
fn simple_range_for_variant(v: &HgvsVariant) -> Option<(Region, i64, i64)> {
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
    range_fn: impl Fn(&Interval<L>) -> Option<(Region, i64, i64)>,
) -> Option<(Region, i64, i64)> {
    if !loc_edit.edit.is_certain() {
        return None;
    }
    let edit = loc_edit.edit.inner()?;
    let (region, start, end) = range_fn(&loc_edit.location)?;
    match edit {
        NaEdit::Substitution { .. }
        | NaEdit::SubstitutionNoRef { .. }
        | NaEdit::Deletion { .. } => Some((region, start, end)),
        NaEdit::Delins { sequence, .. } => {
            sequence.bases()?;
            Some((region, start, end))
        }
        NaEdit::Insertion { sequence } => {
            sequence.bases()?;
            // Anchor for an insertion is [end, start] — empty range at boundary.
            (end == start.checked_add(1)?).then_some((region, end, start))
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
/// callback returns `(region, start, end)` for the location only when it
/// is "simple" (no intronic offsets, certain edit). Returns None when
/// the edit is uncertain, unknown, or not a merge-eligible NaEdit
/// variant. The `region` is propagated through to the merged anchor so
/// `build_*_merged` can reconstruct the right `CdsPos` / `TxPos` /
/// `RnaPos` shape (negative base for 5'UTR / upstream, `utr3` /
/// `downstream` flag for 3'UTR / downstream).
fn anchor_from_loc_edit<L>(
    loc_edit: &LocEdit<Interval<L>, NaEdit>,
    range_fn: impl Fn(&Interval<L>) -> Option<(Region, i64, i64)>,
) -> Option<Anchor> {
    if !loc_edit.edit.is_certain() {
        return None;
    }
    let edit = loc_edit.edit.inner()?;
    let (region, start, end) = range_fn(&loc_edit.location)?;
    match edit {
        NaEdit::Substitution { alternative, .. } | NaEdit::SubstitutionNoRef { alternative } => {
            Some(Anchor {
                region,
                start,
                end,
                alt: vec![*alternative],
            })
        }
        NaEdit::Deletion { .. } => Some(Anchor {
            region,
            start,
            end,
            alt: Vec::new(),
        }),
        NaEdit::Delins { sequence, .. } => {
            let bases = sequence.bases()?.to_vec();
            Some(Anchor {
                region,
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
                region,
                start: end,
                end: start,
                alt: bases,
            })
        }
        _ => None,
    }
}

/// Both endpoints of an interval must share the same `Region` for the
/// interval to be merge-eligible. Cross-region ranges (`c.-1_1`, etc.)
/// have no valid HGVS syntax, so failing this check on a parsed
/// `Interval` indicates upstream malformedness rather than a normal
/// merge barrier; we treat it as ineligible just in case.
fn join_pos(
    start: Option<(Region, i64)>,
    end: Option<(Region, i64)>,
) -> Option<(Region, i64, i64)> {
    let (rs, s) = start?;
    let (re, e) = end?;
    if rs != re {
        return None;
    }
    Some((rs, s, e))
}

fn simple_genome_range(interval: &Interval<GenomePos>) -> Option<(Region, i64, i64)> {
    join_pos(
        simple_genome_pos(&interval.start),
        simple_genome_pos(&interval.end),
    )
}

fn simple_genome_pos(boundary: &UncertainBoundary<GenomePos>) -> Option<(Region, i64)> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_special() || pos.offset.is_some() {
        return None;
    }
    // Genomic coordinates fit comfortably within i64 (positive only).
    let v = i64::try_from(pos.base).ok()?;
    Some((Region::Genome, v))
}

fn simple_cds_range(interval: &Interval<CdsPos>) -> Option<(Region, i64, i64)> {
    join_pos(
        simple_cds_pos(&interval.start),
        simple_cds_pos(&interval.end),
    )
}

fn simple_cds_pos(boundary: &UncertainBoundary<CdsPos>) -> Option<(Region, i64)> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_unknown() || pos.is_intronic() {
        return None;
    }
    if pos.is_3utr() {
        // 3'UTR axis is the `*N` count (positive).
        if pos.base < 1 {
            return None;
        }
        return Some((Region::ThreePrimeUtr, pos.base));
    }
    if pos.base < 0 {
        // 5'UTR axis is the signed CDS coord (negative).
        return Some((Region::FivePrimeUtr, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Cds, pos.base));
    }
    // pos.base == 0 — invalid c. position (HGVS skips c.0).
    None
}

fn simple_tx_range(interval: &Interval<TxPos>) -> Option<(Region, i64, i64)> {
    join_pos(simple_tx_pos(&interval.start), simple_tx_pos(&interval.end))
}

fn simple_tx_pos(boundary: &UncertainBoundary<TxPos>) -> Option<(Region, i64)> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_intronic() {
        return None;
    }
    if pos.is_downstream() {
        if pos.base < 1 {
            return None;
        }
        return Some((Region::TxDownstream, pos.base));
    }
    if pos.base < 0 {
        return Some((Region::TxUpstream, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Tx, pos.base));
    }
    None
}

fn simple_rna_range(interval: &Interval<RnaPos>) -> Option<(Region, i64, i64)> {
    join_pos(
        simple_rna_pos(&interval.start),
        simple_rna_pos(&interval.end),
    )
}

fn simple_rna_pos(boundary: &UncertainBoundary<RnaPos>) -> Option<(Region, i64)> {
    let pos = boundary.as_single().and_then(|mu| match mu {
        Mu::Certain(p) => Some(p),
        _ => None,
    })?;
    if pos.is_intronic() {
        return None;
    }
    if pos.is_3utr() {
        if pos.base < 1 {
            return None;
        }
        return Some((Region::ThreePrimeUtr, pos.base));
    }
    if pos.base < 0 {
        return Some((Region::FivePrimeUtr, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Rna, pos.base));
    }
    None
}

fn build_genome_merged(template: &GenomeVariant, merged: Anchor) -> GenomeVariant {
    debug_assert_eq!(merged.region, Region::Genome);
    let (location, edit) = build_naedit(merged, |_, b| {
        GenomePos::new(u64::try_from(b).expect("genome anchor base must be non-negative"))
    });
    GenomeVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_cds_merged(template: &CdsVariant, merged: Anchor) -> CdsVariant {
    let (location, edit) = build_naedit(merged, |region, b| match region {
        Region::Cds | Region::FivePrimeUtr => CdsPos::new(b),
        Region::ThreePrimeUtr => CdsPos {
            base: b,
            offset: None,
            utr3: true,
        },
        _ => unreachable!("non-c. region {:?} on CdsVariant", region),
    });
    CdsVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_tx_merged(template: &TxVariant, merged: Anchor) -> TxVariant {
    let (location, edit) = build_naedit(merged, |region, b| match region {
        Region::Tx | Region::TxUpstream => TxPos::new(b),
        Region::TxDownstream => TxPos::downstream(b),
        _ => unreachable!("non-n. region {:?} on TxVariant", region),
    });
    TxVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_rna_merged(template: &RnaVariant, merged: Anchor) -> RnaVariant {
    let (location, edit) = build_naedit(merged, |region, b| match region {
        Region::Rna | Region::FivePrimeUtr => RnaPos::new(b),
        Region::ThreePrimeUtr => RnaPos::utr3(b),
        _ => unreachable!("non-r. region {:?} on RnaVariant", region),
    });
    RnaVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

fn build_mt_merged(template: &MtVariant, merged: Anchor) -> MtVariant {
    debug_assert_eq!(merged.region, Region::Genome);
    let (location, edit) = build_naedit(merged, |_, b| {
        GenomePos::new(u64::try_from(b).expect("mitochondrial anchor base must be non-negative"))
    });
    MtVariant {
        accession: template.accession.clone(),
        gene_symbol: template.gene_symbol.clone(),
        loc_edit: LocEdit::new(location, edit),
    }
}

/// Two CDS positions share a codon if they fall in the same 1-indexed
/// triplet: `c.1..3` is codon 1, `c.4..6` is codon 2, etc. The standard
/// codon-number formula is `(base - 1) / 3` (integer division).
///
/// Exposed at `pub(crate)` so the post-merge `decompose_delins`
/// pass can re-apply the spec's codon-frame exception (issue #165 /
/// tracking issue #81 item A10).
pub(crate) fn same_codon(a: i64, b: i64) -> bool {
    if a < 1 || b < 1 {
        return false;
    }
    (a - 1) / 3 == (b - 1) / 3
}

/// Look up the reference base at a CDS-axis position on the head variant's
/// transcript. Used by the codon-frame merge (issue #79, extended for
/// `r.` and chains in #275) to insert the unchanged middle nucleotide
/// between two edits separated by one base on the codon axis.
///
/// Returns `None` if:
///   * the head is not a coding-axis variant (`c.` or `r.`);
///   * the provider has no transcript for the accession;
///   * the transcript has no `cds_start`;
///   * `cds_axis < 1` (the codon-frame eligibility check already gates
///     on `same_codon`, which requires positive values, but we still
///     guard defensively);
///   * the resolved transcript index falls outside the sequence.
///
/// Both `c.` and `r.` share the CDS-relative axis: in this codebase,
/// `RnaPos` (see `hgvs/location.rs`) represents positions on the
/// CDS-relative axis, and `Normalizer::rna_to_tx_pos` (in
/// `normalize/mod.rs`) maps `r.N` to `cds_start + N - 1` in transcript
/// coordinates — identical to the `c.N` mapping. This is a property
/// of the in-repo coordinate model, not a generic HGVS-spec claim:
/// other implementations sometimes use a transcript-relative `r.`
/// axis instead. A single lookup therefore handles both `c.` and
/// `r.`. Failure to look up the byte makes the codon-frame merge
/// silently decline — the input pair passes through unmerged with
/// no error.
fn lookup_codon_middle_ref<P: ReferenceProvider>(
    provider: &P,
    head: &HgvsVariant,
    cds_axis: i64,
) -> Option<Base> {
    let accession = match head {
        HgvsVariant::Cds(v) => &v.accession,
        HgvsVariant::Rna(v) => &v.accession,
        _ => return None,
    };
    if cds_axis < 1 {
        return None;
    }
    let tx = provider
        .get_transcript(&accession.transcript_accession())
        .ok()?;
    let cds_start = tx.cds_start?;
    let tx_idx_1based = cds_start.checked_add(cds_axis as u64)?.checked_sub(1)?;
    let byte = *tx
        .sequence
        .as_deref()?
        .as_bytes()
        .get(tx_idx_1based.checked_sub(1)? as usize)?;
    Base::from_char(byte as char)
}

fn build_naedit<P>(
    merged: Anchor,
    mut to_pos: impl FnMut(Region, i64) -> P,
) -> (Interval<P>, NaEdit) {
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
            deleted: None,
            deleted_length: None,
        }
    };
    let (lo, hi) = if merged.start > merged.end {
        (merged.end, merged.start)
    } else {
        (merged.start, merged.end)
    };
    (
        Interval::new(to_pos(merged.region, lo), to_pos(merged.region, hi)),
        edit,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::variant::Accession;
    use crate::reference::MockProvider;

    #[test]
    fn empty_input_returns_empty() {
        assert!(merge_consecutive_edits(vec![], AllelePhase::Cis, &MockProvider::new()).is_empty());
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
        let out = merge_consecutive_edits(vec![v], AllelePhase::Cis, &MockProvider::new());
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn same_codon_classifies_correctly() {
        assert!(same_codon(1, 3));
        assert!(same_codon(4, 6));
        assert!(same_codon(145, 147));
        assert!(!same_codon(3, 5)); // codon 1 vs codon 2
        assert!(!same_codon(0, 2)); // 0 invalid
        assert!(!same_codon(-3, -1));
        // Large positions exercise the `(base - 1) / 3` formula far from
        // its near-zero edge. c.99997..c.99999 falls in codon 33333;
        // c.99998..c.100000 straddles codon 33333 and 33334.
        assert!(same_codon(99997, 99999));
        assert!(!same_codon(99998, 100000));
    }
}
