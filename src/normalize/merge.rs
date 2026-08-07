//! Merge consecutive sub-variants in a cis allele into a single edit per HGVS spec.
//!
//! See `docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md`.

use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, RepeatCount, Sequence};
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{
    Accession, AllelePhase, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, MtVariant, RnaVariant,
    TxVariant,
};
use crate::normalize::boundary::Boundaries;
use crate::normalize::config::ShuffleDirection;
use crate::normalize::shuffle::shuffle;
use crate::normalize::{boundary_delins_bases, BoundarySide, SequenceEnds};
use crate::reference::ReferenceProvider;
use crate::sequence::complement_base;

/// Coordinate-system region used as the merge-eligibility key.
///
/// Adjacency in the merge pass requires both ends to share a region, so we
/// tag each position with its region and refuse to merge across.
///
/// The integer `start`/`end` axis is region-local: positive for `Cds` /
/// `ThreePrimeUtr` / `Tx` / `TxDownstream`, negative for `FivePrimeUtr`
/// / `TxUpstream`. Adjacency `prev.end + 1 == next.start` works
/// naturally for all six because the axis is monotonic 5'→3' within
/// each region (`c.-3 → c.-2 → c.-1` maps to `-3 → -2 → -1`).
///
/// That refusal is therefore a **property of this pass's arithmetic, not of
/// HGVS**: `prev.end + 1 == next.start` is only meaningful within one region,
/// since `c.-1` and `c.1` are adjacent nucleotides whose axis values differ by
/// 2. A region-spanning range is perfectly valid nomenclature — the spec writes
/// `NM_001849.3:c.-1_*1=` (`consultation/SVD-WG001.md:37`) and
/// `NM_004006.2:c.-244_*2691del` (`consultation/open-issues.md:53`), and
/// ferro's own parser accepts both. (This comment previously asserted the
/// opposite, that `c.-1_1` "does not exist as valid range syntax"; corrected
/// 2026-08-02, and [`respell_at_gap`] now deliberately *emits* such a range
/// when a repair's junction lands on the CDS start.)
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

/// How an [`Anchor`]'s span and `alt` relate — i.e. which HGVS form
/// [`build_naedit`] must emit for it.
///
/// Everything the *merge* pass produces is a [`AnchorForm::Replacement`]: it
/// coalesces members into one span-plus-replacement and lets the two lengths
/// pick substitution / deletion / insertion / delins. The sequence-first
/// derivation additionally recognises a shape those lengths cannot express — a
/// tandem duplication, which `duplication.md:18` states as a MUST — so the form
/// travels with the anchor rather than being re-inferred downstream, where the
/// duplicated source span is no longer in hand.
///
/// An **inversion** is deliberately *not* a form here, even though it is the
/// other shape the two lengths cannot see. It does not need to be: a `dup` is
/// unrecoverable downstream because its source span is gone by then, whereas an
/// `inv` is a span plus its own bases, so [`crate::normalize::rules`]'s
/// single-span typing re-derives it from the rendered member. A variant was
/// tried and measured inert — 849,752 normalizations across the whole-span
/// inversion and random-haplotype corpora were byte-identical with and without
/// it — so it was dropped rather than kept as a second, silent authority for a
/// typing another module already owns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AnchorForm {
    /// The span is replaced by `alt`; typing follows from the two lengths.
    Replacement,
    /// The span is duplicated in tandem. `alt` is empty and unused: a `dup`
    /// names only the source bases, which the span already gives.
    Duplication,
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
    form: AnchorForm,
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
            reconcile_head(&mut output, head_anchor.take().unwrap(), provider);
            head_merged = false;
        }
        head_anchor = anchor_for_variant(&next);
        output.push(next);
    }
    if head_merged {
        reconcile_head(&mut output, head_anchor.take().unwrap(), provider);
    }
    output
}

/// One extracted genomic cis edit for the overlap-collapse pass.
enum GEdit {
    /// Insertion between reference positions `gap` and `gap + 1`.
    Ins { gap: i64, alt: Vec<Base> },
    /// Deletion of the inclusive 1-based span `[s, e]`.
    Del { s: i64, e: i64 },
    /// Single-base substitution at `pos`.
    Sub { pos: i64, alt: Base },
    /// `delins` of the inclusive span `[s, e]` with replacement `alt`.
    Delins { s: i64, e: i64, alt: Vec<Base> },
    /// Tandem duplication of the inclusive span `[s, e]`. Semantically an
    /// insertion of `ref[s..=e]` at gap `e` (the 3' copy); its alt bases are
    /// read from the reference window once fetched, so it participates in the
    /// collapse exactly like an `Ins { gap: e, .. }`.
    Dup { s: i64, e: i64 },
    /// Inversion of the inclusive span `[s, e]`: replaced by its own reverse
    /// complement, read from the reference window once fetched. Only the
    /// sequence-first canonicalizer builds this; `collapse_overlapping_cis_edits`
    /// refuses inversion members.
    Inv { s: i64, e: i64 },
    /// A tandem repeat member (`274A[3]`) that has **not been lowered yet**:
    /// `unit` repeated `count` times, anchored on the inclusive span `[s, e]`.
    ///
    /// The only `GEdit` that does not yet denote a sequence. A repeat names its
    /// reference tract implicitly — "the run of `unit` covering `[s, e]`" — and
    /// the tract's extent is in the reference, not in the description, so this
    /// cannot be resolved until the window is fetched. [`lower_repeat_edits`]
    /// rewrites every one of these into the `Ins`/`Del` it denotes; nothing
    /// downstream of that call can see one, and both the applier and
    /// `collapse_overlapping_cis_edits` refuse it outright rather than guess.
    Repeat {
        s: i64,
        e: i64,
        unit: Vec<u8>,
        count: u64,
    },
}

/// Report whether `variants` contains two or more insertion-like cis members
/// (`ins` or tandem `dup`) that attach at the **same** reference gap.
///
/// Such a pair is order-ambiguous (`[g insA; g insC]` describes `...AC...` or
/// `...CA...`), so the collapse deliberately refuses it (see the
/// `seen_insertion_gaps` guard below, issue #487). That refusal keys off the
/// members' *current* positions, but an independent per-member 3'-shift can
/// move two originally-same-gap insertions apart — after which a re-collapse
/// would no longer see the collision and would silently pick an order. Callers
/// that iterate collapse across shifts (see `normalize_allele`) use this to
/// detect the ambiguity on the **pre-shift** input and decline the extra
/// passes, preserving both the #487 refusal and the #180 merge-first invariant.
///
/// Conservative: only the head member's axis kind is considered, and any member
/// that isn't a simple certain `ins`/`dup` on that axis is ignored (it cannot be
/// half of a same-gap insertion pair).
pub(crate) fn has_same_gap_insertions(variants: &[HgvsVariant], phase: AllelePhase) -> bool {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return false;
    }
    let kind = match &variants[0] {
        HgvsVariant::Genome(_) => CisKind::Genome,
        HgvsVariant::Mt(_) => CisKind::Mt,
        HgvsVariant::Cds(_) => CisKind::Cds,
        HgvsVariant::Tx(_) => CisKind::Tx,
        HgvsVariant::Rna(_) => CisKind::Rna,
        _ => return false,
    };
    let mut seen = std::collections::HashSet::new();
    for v in variants {
        let Some((_, _, s, e, edit)) = cis_axis_parts(v, kind) else {
            continue;
        };
        let gap = match edit {
            NaEdit::Insertion { .. } if e == s + 1 => s,
            NaEdit::Duplication {
                uncertain_extent: None,
                ..
            } => e,
            _ => continue,
        };
        if !seen.insert(gap) {
            return true;
        }
    }
    false
}

/// Collapse a contiguous run of *overlapping* cis genomic edits — insertions
/// flanking/within a deletion or substitution at one locus — into a single
/// `delins`, by applying them to the reference and re-deriving the minimal
/// edit. This is the case the strict-adjacency chain in
/// `merge_consecutive_edits` cannot express (it only merges non-overlapping
/// consecutive edits). Example (mutalyzer-verified): `g.[104_105insA;
/// 105_106insC;105del]` -> `g.105delinsAC`.
///
/// Deliberately conservative and all-or-nothing — it returns the input
/// unchanged unless **every** sub-variant participates in one collapsible
/// group:
///   * all sub-variants are same-accession genomic `NaEdit`s with simple
///     (certain, non-special, non-offset) positions, and one of
///     ins/**dup**/del/sub/delins, and where the edit carries a payload at all
///     that payload is literal (a `dup` carries none — it copies the range it
///     names from the reference);
///   * the del/sub/delins spans form one contiguous interval with no
///     unchanged interior base (no holes);
///   * every insertion-like member attaches to that interval (its gap lies in
///     `[c_lo - 1, c_hi]`);
///   * the group has BOTH an insertion-like member and a del/sub/delins.
///
/// **`dup` is insertion-like here, not a span (#1436.)** The list above used to
/// omit it, which read as an exhaustive precondition and told anyone tracing
/// whether a `dup` can reach this collapse that it cannot. It can: the body
/// treats it as one of the members that attaches to the del/sub/delins interval
/// —
///
/// ```ignore
/// let is_insertion_like = |e: &GEdit| matches!(e, GEdit::Ins { .. } | GEdit::Dup { .. });
/// ```
///
/// — plus four further `GEdit::Dup` arms in the span and gap computations. That
/// is the same read/write distinction #1411 drew: a `dup` writes at the junction
/// 3' of the run it reads, so it attaches to an interval rather than forming
/// one.
///
/// Anything else passes through untouched, so pure-insertion groups (including
/// pure-`dup` ones, which have no del/sub/delins to attach to), pure-deletion
/// groups (owned by `merge_consecutive_edits`), and non-overlapping inputs are
/// unaffected.
pub(crate) fn collapse_overlapping_cis_edits<P: ReferenceProvider>(
    variants: Vec<HgvsVariant>,
    phase: AllelePhase,
    provider: &P,
) -> Vec<HgvsVariant> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return variants;
    }
    // The single-axis genomic (`g.`) / mitochondrial (`m.`) systems share
    // `Region::Genome`; the transcript axes (`c.`, `n.`, `r.`) each have a
    // positive body region (`Cds` / `Tx` / `Rna`). Track which kind the group
    // is so we rebuild the correct variant at the end and translate window
    // coordinates correctly; a mixed-kind group is refused (mirrors the
    // accession check). Transcript members are collapsed only within the
    // positive body (issue #920) — `collect_canonical_edits` enforces that with
    // its `region != body` refusal.
    let Some(kind) = cis_kind_of(&variants[0]) else {
        return variants;
    };
    // The positive body region every member must lie in for this axis.
    let body = body_region(kind);
    // Borrow the head's accession as the template once; the helper rejects any
    // variant whose kind doesn't match `kind`.
    let Some((template_accession, _, _, _, _)) = cis_axis_parts(&variants[0], kind) else {
        return variants;
    };
    let template_accession = template_accession.clone();

    // Lowering is `collect_canonical_edits`': same `cis_axis_parts` call, same
    // accession / body-region refusals, same per-`NaEdit` guards. It admits two
    // shapes this path does not — `NaEdit::Inversion` and `NaEdit::Repeat`,
    // both of which the sequence-first canonicalizer serves — so reject those
    // explicitly and the two passes stay one implementation instead of two
    // copies ten screens apart.
    //
    // A `GEdit::Repeat` does not yet denote a sequence at all (its tract is in
    // the reference, and only `lower_repeat_edits` reads it), so this pass —
    // which reasons purely about spans — could not use one even if it wanted
    // to.
    let Some(edits) = collect_canonical_edits(&variants, kind, body, &template_accession) else {
        return variants;
    };
    if edits
        .iter()
        .any(|e| matches!(e, GEdit::Inv { .. } | GEdit::Repeat { .. }))
    {
        return variants;
    }

    // Require a genuinely mixed group: at least one insertion AND at least one
    // replacement (del/sub/delins). This excludes pure-insertion groups (which
    // must stay separate when separated by unchanged bases) and pure-deletion /
    // pure-substitution groups (owned by `merge_consecutive_edits`).
    let is_insertion_like = |e: &GEdit| matches!(e, GEdit::Ins { .. } | GEdit::Dup { .. });
    let has_ins = edits.iter().any(&is_insertion_like);
    let has_repl = edits.iter().any(|e| !is_insertion_like(e));
    if !has_ins || !has_repl {
        return variants;
    }

    // Changed interval = union of replacement spans; must be contiguous.
    let covered: Vec<(i64, i64)> = edits
        .iter()
        .filter_map(|e| match e {
            // `Inv` and `Repeat` are unreachable here — the rejection above
            // drops any group carrying either — but both hold reference bases
            // over `[s, e]`, so classify them as replacements rather than leave
            // a misleading `None` that would let the contiguity check pass over
            // territory a member owns.
            GEdit::Del { s, e }
            | GEdit::Delins { s, e, .. }
            | GEdit::Inv { s, e }
            | GEdit::Repeat { s, e, .. } => Some((*s, *e)),
            GEdit::Sub { pos, .. } => Some((*pos, *pos)),
            GEdit::Ins { .. } | GEdit::Dup { .. } => None,
        })
        .collect();
    let c_lo = covered.iter().map(|(s, _)| *s).min().unwrap();
    let c_hi = covered.iter().map(|(_, e)| *e).max().unwrap();
    // No unchanged interior base: every position in [c_lo, c_hi] is covered.
    for p in c_lo..=c_hi {
        if !covered.iter().any(|(s, e)| *s <= p && p <= *e) {
            return variants;
        }
    }
    // Every insertion must attach to the changed interval. Two insertions at
    // the *same* gap would concatenate into `after[idx(gap)]` in member order,
    // making the collapsed result order-dependent (e.g. `[gap insA; gap insB]`
    // vs the reverse yields `...AB...` vs `...BA...`). Refuse such a group —
    // matches the conservative all-or-nothing philosophy and preserves the
    // member-order invariance the collapse otherwise guarantees. A `Dup { s, e }`
    // attaches at gap `e` (the 3' copy) and reads its source span `[s, e]`.
    let mut seen_insertion_gaps = std::collections::HashSet::new();
    for e in &edits {
        let gap = match e {
            GEdit::Ins { gap, .. } => *gap,
            GEdit::Dup { e, .. } => *e,
            _ => continue,
        };
        if gap < c_lo - 1 || gap > c_hi {
            return variants;
        }
        if !seen_insertion_gaps.insert(gap) {
            return variants;
        }
    }

    // Window covers the changed interval plus all insertion flanks (and, for a
    // dup, its full duplicated source span so its alt bases can be read).
    let mut w_lo = c_lo;
    let mut w_hi = c_hi;
    for e in &edits {
        match e {
            GEdit::Ins { gap, .. } => {
                w_lo = w_lo.min(*gap);
                w_hi = w_hi.max(*gap + 1);
            }
            GEdit::Dup { s, e } => {
                w_lo = w_lo.min(*s);
                w_hi = w_hi.max(*e + 1);
            }
            _ => {}
        }
    }
    if w_lo < 1 {
        return variants;
    }

    // Fetch the reference window [w_lo, w_hi] (1-based inclusive on the axis)
    // as a 0-based half-open slice of the underlying sequence. Bail on any
    // provider miss or short read.
    //
    // The window is expressed in the axis's own coordinates; the underlying
    // sequence offset differs per axis:
    //   * `g.`/`m.` and `n.` (`Tx`): the axis IS the fetched sequence, so a
    //     1-based position `N` is sequence offset `N` (delta = 0).
    //   * `c.` (`Cds`) and `r.` (`Rna`): the axis is CDS-relative, so a
    //     1-based position `N` maps to transcript position `cds_start + N - 1`
    //     (delta = `cds_start - 1`) — the same mapping `lookup_codon_middle_ref`
    //     uses. Both regions have been restricted to the positive body above.
    let accession = template_accession.transcript_accession();
    // Only the offset matters here — this collapse has no codon exception to
    // gate, so the frame's `reading_frame` half is not consulted.
    let Some(delta) = axis_frame(kind, &template_accession, provider).map(|frame| frame.delta)
    else {
        return variants;
    };
    let start0 = w_lo + delta - 1;
    let end0 = w_hi + delta;
    if start0 < 0 {
        return variants;
    }
    let Ok(ref_seq) = provider.get_sequence(&accession, start0 as u64, end0 as u64) else {
        return variants;
    };
    let ref_bytes = ref_seq.as_bytes();
    let n = (w_hi - w_lo + 1) as usize;
    if ref_bytes.len() != n {
        return variants;
    }

    // Apply edits over the window: `cell[i]` is the kept ref byte (or None if
    // deleted) for position `w_lo + i`; `after[i]` is bases inserted after that
    // position; `before` is bases inserted before the window's first base.
    let idx = |p: i64| -> usize { (p - w_lo) as usize };
    let mut cell: Vec<Option<u8>> = ref_bytes.iter().map(|b| Some(*b)).collect();
    let mut after: Vec<Vec<u8>> = vec![Vec::new(); n];
    let mut before: Vec<u8> = Vec::new();
    for e in &edits {
        match e {
            GEdit::Del { s, e } => {
                for p in *s..=*e {
                    cell[idx(p)] = None;
                }
            }
            GEdit::Sub { pos, alt } => cell[idx(*pos)] = Some(alt.to_u8()),
            // Unreachable — the group was rejected above — and unapplicable:
            // a repeat's tract is resolved only by `lower_repeat_edits`, which
            // this pass never runs. Decline rather than approximate.
            GEdit::Repeat { .. } => return variants,
            GEdit::Delins { s, e, alt } => {
                for p in *s..=*e {
                    cell[idx(p)] = None;
                }
                let bytes: Vec<u8> = alt.iter().map(|b| b.to_u8()).collect();
                if *s > w_lo {
                    after[idx(*s - 1)].extend(bytes);
                } else {
                    before.extend(bytes);
                }
            }
            GEdit::Ins { gap, alt } => {
                let bytes: Vec<u8> = alt.iter().map(|b| b.to_u8()).collect();
                if *gap >= w_lo {
                    after[idx(*gap)].extend(bytes);
                } else {
                    before.extend(bytes);
                }
            }
            GEdit::Dup { s, e } => {
                // Tandem copy of ref[s..=e] inserted immediately 3' of `e`.
                let bytes: Vec<u8> = (*s..=*e).map(|p| ref_bytes[idx(p)]).collect();
                after[idx(*e)].extend(bytes);
            }
            // Unreachable: the explicit rejection above drops any group with an
            // inversion member. Refusing keeps the match exhaustive without
            // carrying an inversion implementation no caller can reach — and
            // refusing, not panicking, is this module's convention.
            GEdit::Inv { .. } => return variants,
        }
    }
    let mut variant: Vec<u8> = before;
    for i in 0..n {
        if let Some(b) = cell[i] {
            variant.push(b);
        }
        variant.extend(&after[i]);
    }

    // Minimal-trim diff vs the reference window.
    let (lo, hi_ref, hi_var) = trim_common_flanks(ref_bytes, &variant);
    let alt_bases: Vec<Base> = variant[lo..hi_var]
        .iter()
        .filter_map(|b| Base::from_char(*b as char))
        .collect();
    if alt_bases.len() != hi_var - lo {
        return variants; // non-IUPAC byte from the reference; refuse.
    }
    let del_start = w_lo + lo as i64;
    let (a_start, a_end) = if hi_ref == lo {
        // Net pure insertion at the boundary: anchor span = 1 nt (start = end+1).
        (del_start, del_start - 1)
    } else {
        (del_start, w_lo + hi_ref as i64 - 1)
    };
    // The collapse has its own way of naming a position past the end (#1327):
    // when the group nets to a pure insertion at the window's 3' edge,
    // `del_start` is one past the last base, and the insertion-shaped anchor
    // `(del_start, del_start - 1)` renders as `m.16569_16570=` on rCRS.
    //
    // Distinct from the `respell_at_gap` overrun this issue names — a different
    // pass, reached without it — but the same defect, so it is fixed here too
    // rather than left for the symptom to survive the fix. Refusing the
    // collapse is the right answer *here* (unlike in `respell_at_gap`, where it
    // reintroduced #1286): the members are left as they are, which is in range,
    // and refusing is this module's convention for a group it cannot place.
    if let Ok(length) = provider.get_sequence_length(&accession) {
        let past_end = |p: i64| u64::try_from(p + delta).is_ok_and(|v| v > length);
        if past_end(a_start) || past_end(a_end) {
            return variants;
        }
    }
    let anchor = Anchor {
        region: body,
        start: a_start,
        end: a_end,
        alt: alt_bases,
        form: AnchorForm::Replacement,
    };
    // Rebuild the variant kind the group came in as. `g.`/`m.` use the
    // single-axis `Region::Genome` anchor; `c.`/`n.`/`r.` use the positive
    // body region threaded into `anchor` so the builder reconstructs the
    // right position shape. The head carries the accession / gene symbol to
    // seed the merged variant. `kind` was derived from `variants[0]` and
    // every member passed `cis_axis_parts(_, kind)`, so the head necessarily
    // matches `kind` here.
    match (kind, &variants[0]) {
        (CisKind::Genome, HgvsVariant::Genome(g)) => {
            vec![HgvsVariant::Genome(build_genome_merged(g, anchor))]
        }
        (CisKind::Mt, HgvsVariant::Mt(m)) => vec![HgvsVariant::Mt(build_mt_merged(m, anchor))],
        (CisKind::Cds, HgvsVariant::Cds(c)) => vec![HgvsVariant::Cds(build_cds_merged(c, anchor))],
        (CisKind::Tx, HgvsVariant::Tx(t)) => vec![HgvsVariant::Tx(build_tx_merged(t, anchor))],
        (CisKind::Rna, HgvsVariant::Rna(r)) => vec![HgvsVariant::Rna(build_rna_merged(r, anchor))],
        _ => variants,
    }
}

/// Which variant kind a cis-collapse group is. `g.` and `m.` share
/// `Region::Genome` and the same `GenomeInterval`/`NaEdit` machinery; the
/// transcript kinds (`c.`/`n.`/`r.`) each collapse within their positive body
/// region (issue #920). Every kind rebuilds to a distinct variant kind, and a
/// mixed-kind group is refused.
#[derive(Clone, Copy, PartialEq, Eq)]
enum CisKind {
    Genome,
    Mt,
    Cds,
    Tx,
    Rna,
}

/// The positive body region a cis-collapse group of this kind must lie in.
/// The collapse only ever operates within this single region — 5'UTR / 3'UTR
/// / intronic-offset / upstream / downstream members are refused (issue #920).
fn body_region(kind: CisKind) -> Region {
    match kind {
        CisKind::Genome | CisKind::Mt => Region::Genome,
        CisKind::Cds => Region::Cds,
        CisKind::Tx => Region::Tx,
        CisKind::Rna => Region::Rna,
    }
}

/// Extract the `(accession, region, start, end, edit)` of a variant for the
/// cis-collapse pass, but only when its kind matches `kind` and the edit is
/// certain. The `(region, start, end)` are the location's *raw* endpoints on
/// the axis (unlike `simple_range_for_variant`, which swaps insertion
/// endpoints into anchor form); the collapse reads insertion gaps from the
/// raw endpoints. Returns `None` for any other variant kind, a kind mismatch,
/// a disqualifying position (offset / special / mixed-region), or an
/// uncertain / non-`NaEdit` edit — the caller treats that as "decline to
/// collapse".
///
/// Written as [`member_axis_endpoints`] plus [`join_pos`], because that *is*
/// what it was: the per-axis `simple_*_range` helpers each fold their two ends
/// through `join_pos`, and folding them here instead keeps the two readers from
/// drifting apart. The single-region requirement is the whole of the difference
/// between them (#1482) — which is worth being able to see in one place, since
/// it is also the whole of the defect that motivated the other reader.
fn cis_axis_parts(
    v: &HgvsVariant,
    kind: CisKind,
) -> Option<(&Accession, Region, i64, i64, &NaEdit)> {
    let (accession, start, end, edit) = member_axis_endpoints(v, kind)?;
    let (region, s, e) = join_pos(Some(start), Some(end))?;
    Some((accession, region, s, e, edit))
}

/// Whether a merged anchor restates the reference bases under its own span —
/// that is, whether the members it merged cancelled each other out.
///
/// `g.[262_263insA;263del]` over a reference whose base 263 is `A` inserts an
/// `A` and deletes the `A` beside it. The merge is faithful: it yields
/// `g.263delinsA`, which denotes exactly the reference. Per-member
/// normalization then renders that as `g.263=`, and inside a cis allele with
/// real members beside it, that `=` is a residue of the merge rather than
/// anything the caller said (#1321).
///
/// Only the axes whose positions index the fetched sequence directly are
/// answered; anything else returns `false` so the member is built as before.
fn merged_anchor_restates_reference<P: ReferenceProvider>(
    head: &HgvsVariant,
    anchor: &Anchor,
    provider: &P,
) -> bool {
    // An insertion anchor with nothing left to insert is the empty-payload
    // cancellation `build_naedit` renders as an identity; it has no span to
    // compare against, and `start > end` is its shape.
    if anchor.start > anchor.end {
        return anchor.alt.is_empty();
    }
    let Some(kind) = cis_kind_of(head) else {
        return false;
    };
    let Some((accession, region, _start, _end, _edit)) = cis_axis_parts(head, kind) else {
        return false;
    };
    // `region_sequence_delta`, not `axis_frame(..).delta`: the anchor's positions
    // are on the member's own region axis, and on `c.` that axis is not one
    // affine shift. A `c.-n` 5'UTR or `c.*n` 3'UTR position needs the
    // region-aware conversion — the same one `first_out_of_bounds_coordinate`
    // runs before any comparison — or the bases fetched below are simply the
    // wrong ones, and a member that merely *looks* like a cancellation against
    // them would be dropped.
    let provider_key = accession.transcript_accession();
    let Some(delta) = region_sequence_delta(region, &provider_key, provider) else {
        return false;
    };
    // Checked for the same reason the conversion below is: `anchor.start`/
    // `anchor.end` come off a parsed description, so an adversarial span can
    // overflow the subtraction or the inclusive `+ 1` and panic in a debug build
    // where declining is the answer every other unrepresentable coordinate gets.
    let Some(span) = anchor
        .end
        .checked_sub(anchor.start)
        .and_then(|d| d.checked_add(1))
        .and_then(|d| usize::try_from(d).ok())
    else {
        return false;
    };
    // A cancellation replaces each base it claims with that same base, so a
    // payload of a different length cannot be one.
    if anchor.alt.len() != span {
        return false;
    }
    // Checked, matching the sibling conversions in this file: `anchor.start`/
    // `anchor.end` come off a parsed description and `delta` off the record's CDS
    // bounds, so neither is bounded by anything here — an unchecked `i64` add
    // panics in a debug build where declining is the answer every other
    // unrepresentable coordinate gets.
    let Some(start0) = anchor
        .start
        .checked_add(delta)
        .and_then(|s| s.checked_sub(1))
    else {
        return false;
    };
    let Some(end0) = anchor.end.checked_add(delta) else {
        return false;
    };
    if start0 < 0 || end0 < start0 {
        return false;
    }
    let Ok(ref_seq) = provider.get_sequence(&provider_key, start0 as u64, end0 as u64) else {
        return false;
    };
    let ref_bytes = ref_seq.as_bytes();
    if ref_bytes.len() != span {
        return false;
    }
    anchor
        .alt
        .iter()
        .zip(ref_bytes)
        .all(|(base, byte)| (base.to_char() as u8).eq_ignore_ascii_case(byte))
}

/// Replace `output.last()` with a freshly-built variant from `anchor`, or drop
/// it when the members it merged cancelled each other out.
///
/// Caller has established that the head is merge-eligible (so kind dispatch
/// in `build_merged` is safe).
///
/// A chain that cancels leaves nothing to describe, so keeping it costs
/// confluence: the same variant spelled without the cancelling members
/// normalizes without the `=` they collapse to, and the two spellings then
/// disagree (#1321). Dropping it here — where the merge is what produced it — is
/// what separates it from a **user-authored** identity. `g.[258=;260del]` and
/// `g.[259=;259_260insG]` never pass through a merge and keep their `=`, which is
/// the guard #1297 installed; `drop_identity_members_covered_by_siblings`
/// approximates the same distinction afterwards by asking whether a sibling
/// covers the member, and a sibling that shifted away from it (as the
/// duplication does under the 5' rule) leaves the residue unrecognised.
///
/// The last member is never dropped: an allele that cancels away entirely still
/// has to render as something, and that something is the identity.
fn reconcile_head<P: ReferenceProvider>(
    output: &mut Vec<HgvsVariant>,
    anchor: Anchor,
    provider: &P,
) {
    let head = output.last().expect("head must exist when head_merged");
    if output.len() > 1 && merged_anchor_restates_reference(head, &anchor, provider) {
        output.pop();
        return;
    }
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

/// Raw location `(region, start, end)` for a sub-variant, regardless of edit
/// kind. Unlike [`simple_range_for_variant`] it does not filter by edit kind or
/// swap insertion endpoints — it is the positional fallback used by
/// [`cis_merge_order_key`] to place non-merge-eligible members (`dup` / `inv` /
/// `repeat` / …) into genomic order too, so a merge barrier lands at its own
/// coordinate instead of wherever it was authored.
fn location_range_for_variant(v: &HgvsVariant) -> Option<(Region, i64, i64)> {
    match v {
        HgvsVariant::Genome(g) => simple_genome_range(&g.loc_edit.location),
        HgvsVariant::Cds(c) => simple_cds_range(&c.loc_edit.location),
        HgvsVariant::Tx(t) => simple_tx_range(&t.loc_edit.location),
        HgvsVariant::Rna(r) => simple_rna_range(&r.loc_edit.location),
        HgvsVariant::Mt(m) => simple_genome_range(&m.loc_edit.location),
        _ => None,
    }
}

/// Total-order key that puts cis members into the canonical **merge** order
/// (#1103), so `merge_consecutive_edits`' greedy left-to-right adjacency walk
/// fires every valid merge regardless of input order.
///
/// The primary axis is the merge **anchor**'s `(end, start)` (not the display
/// start point [`crate::hgvs::variant::cis_member_order_key`] uses): an edit
/// ending at `X` must precede one starting at `X + 1`, and — crucially for a
/// co-located ins/del — a span edit *at* a locus `p` must precede an insertion
/// in the gap 3' of it. Insertions anchor as `[p + 1, p]`, so their `end = p`
/// ties with a span at `p` and their `start = p + 1` breaks the tie *after* it;
/// an insertion in the gap 5' of `p` (`(p-1)_p ins`) anchors as `[p, p - 1]`,
/// so its `end = p - 1` sorts it *before* the span. That is exactly what turns
/// `[104_105insA; 105del; 105_106insC]` into the chainable order
/// `104_105insA → 105del → 105_106insC` (→ `105delinsAC`) for every input
/// permutation. The rendered descriptor is the final tie-break, making the
/// order total so it never falls back to input order. Members with no simple
/// location (special/uncertain/offset-only positions this pass cannot place)
/// get a `no_position` sentinel and sort last, still deterministically.
fn cis_merge_order_key(v: &HgvsVariant) -> (bool, Region, i64, i64, String) {
    let descriptor = format!("{v}");
    // Prefer the merge anchor (insertion endpoints swapped to `[p + 1, p]`);
    // fall back to the raw location for non-merge-eligible kinds.
    match simple_range_for_variant(v).or_else(|| location_range_for_variant(v)) {
        Some((region, anchor_start, anchor_end)) => {
            (false, region, anchor_end, anchor_start, descriptor)
        }
        // `Region::Genome` is an unused placeholder here — `no_position == true`
        // already segregates these, and the `i64::MAX` ties defer to the
        // descriptor tie-break.
        None => (true, Region::Genome, i64::MAX, i64::MAX, descriptor),
    }
}

/// Sort cis members into the canonical merge order (see [`cis_merge_order_key`])
/// so `merge_consecutive_edits` is input-order-independent (#1103).
///
/// A no-op unless every member shares a single accession: `merge_consecutive_edits`
/// only combines same-accession members, and a mixed-accession bracketed allele
/// (`[NC_…;NM_…]`, #218/#219) is left in authored order by the post-merge display
/// sort too, so reordering it here would leak an arbitrary order into the output.
/// Callers additionally gate on cis / certain / not-overlap-conflicting (see
/// `Normalizer::normalize_allele`).
pub(crate) fn sort_cis_members_for_merge(variants: &mut [HgvsVariant]) {
    let first_accession = variants
        .first()
        .and_then(|v| v.accession().map(|a| a.full()));
    let single_accession = first_accession.is_some()
        && variants
            .iter()
            .all(|v| v.accession().map(|a| a.full()) == first_accession);
    if single_accession {
        variants.sort_by_key(cis_merge_order_key);
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
/// Returns `true` when a `Region::Genome` interval has `start > end`,
/// which identifies a wraparound mitochondrial range (e.g. `m.16569_1del`)
/// whose raw endpoints must not participate in linear-adjacency merging.
///
/// Defensive: span/indel math is covered separately by `tests/mito_circular_audit.rs`
/// and `tests/issue_399_mt_circular_followup.rs`; this guard exists so a future
/// refactor doesn't accidentally re-introduce silent-wrong merges across the
/// origin. `HgvsVariant::Circular` (`o.`) variants are excluded from merging
/// at the dispatch level above and never reach this check.
#[inline]
fn is_wraparound_genome(region: Region, start: i64, end: i64) -> bool {
    start > end && matches!(region, Region::Genome)
}

fn simple_range_for_loc_edit<L>(
    loc_edit: &LocEdit<Interval<L>, NaEdit>,
    range_fn: impl Fn(&Interval<L>) -> Option<(Region, i64, i64)>,
) -> Option<(Region, i64, i64)> {
    if !loc_edit.edit.is_certain() {
        return None;
    }
    let edit = loc_edit.edit.inner()?;
    let (region, start, end) = range_fn(&loc_edit.location)?;
    if is_wraparound_genome(region, start, end) {
        return None;
    }
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
    if is_wraparound_genome(region, start, end) {
        return None;
    }
    match edit {
        NaEdit::Substitution { alternative, .. } | NaEdit::SubstitutionNoRef { alternative } => {
            Some(Anchor {
                region,
                start,
                end,
                alt: vec![*alternative],
                form: AnchorForm::Replacement,
            })
        }
        NaEdit::Deletion { .. } => Some(Anchor {
            region,
            start,
            end,
            alt: Vec::new(),
            form: AnchorForm::Replacement,
        }),
        NaEdit::Delins { sequence, .. } => {
            let bases = sequence.bases()?.to_vec();
            Some(Anchor {
                region,
                start,
                end,
                alt: bases,
                form: AnchorForm::Replacement,
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
                form: AnchorForm::Replacement,
            })
        }
        _ => None,
    }
}

/// Fold two endpoints into one interval, **refusing a pair that straddles a
/// region boundary**.
///
/// The refusal is a limit of the *merge* path, not a property of HGVS.
///
/// This comment used to say cross-region ranges "have no valid HGVS syntax, so
/// failing this check on a parsed `Interval` indicates upstream malformedness".
/// **That is false**, and #1482 is what it cost. `c.15_*1del` deletes across a
/// stop codon — an ordinary and common real variant — and `c.-1_1`, which the
/// old comment offered as its example of the impossible thing, is written out in
/// `consultation/SVD-WG001.md:37`.
///
/// Refusing here is still right: the merge coalesces members onto one axis and
/// has no representation for a span numbered on two. What was wrong was
/// inheriting the refusal *and its justification* elsewhere. [`member_span`]
/// did, and because every sibling-awareness pass drops a `None` through
/// `filter_map` rather than declining, a boundary-crossing member became
/// invisible **as a sibling** module-wide — which is how
/// `c.[15_*1del;15_*1insCC]` came to normalize into an allele ferro's own parser
/// rejects. That reader now converts each endpoint onto the sequence axis
/// instead; see [`member_span`] for the whole argument.
///
/// So: use this when the caller genuinely needs one region for one axis, and
/// read the `None` as "this pass cannot express that span", never as "that span
/// is malformed".
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
            special: None,
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
    // A duplication is decided by the *builder*, from the reference, not by the
    // two span lengths below — `dup` and the `ins` of the same bases have the
    // same lengths and `duplication.md:18` says only one of them is allowed.
    let edit = if merged.form == AnchorForm::Duplication {
        NaEdit::Duplication {
            // Canonical form states neither the bases nor the length: the span
            // already names them, and every other builder here likewise leaves
            // the optional reference-bases channel empty.
            sequence: None,
            length: None,
            uncertain_extent: None,
        }
    } else if merged.start > merged.end {
        debug_assert_eq!(
            merged.start,
            merged.end + 1,
            "invariant: insertion anchor span = 1 nt"
        );
        if merged.alt.is_empty() {
            // The members cancelled out: an insertion anchor with nothing left
            // to insert (e.g. `g.[100del;100_101insA]` over a homopolymer,
            // where the deletion 3'-shifts onto the insertion point). Emitting
            // `NaEdit::Insertion` here would build an edit with an empty
            // sequence — not a valid HGVS description (it renders as `ins`
            // with no bases), and it later divides by zero when
            // `normalize_na_edit` rotates the inserted sequence through a
            // repeat. The sequence is unchanged across the anchor, so describe
            // it as an identity (`DNA/other.md`, `c.123_145=`).
            //
            // Mirrors the protein coalescer below, which likewise refuses to
            // emit an empty `delins` when every member of a run is a deletion.
            NaEdit::Identity {
                sequence: None,
                whole_entity: false,
            }
        } else {
            NaEdit::Insertion {
                sequence: InsertedSequence::Literal(Sequence::new(merged.alt)),
            }
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
            substitution_reference: None,
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

/// Coalesce runs of consecutive-residue changes in a **cis protein allele**
/// into a single delins (or a pure range deletion), per
/// `protein/substitution.md:23` / `protein/delins.md:18`:
///
/// > changes involving two or more consecutive amino acids are described as a
/// > deletion/insertion variant (delins) […] the description
/// > `p.Arg76_Cys77delinsSerTrp` is correct, the description
/// > `p.[Arg76Ser;Cys77Trp]` is not correct.
///
/// "changes" is general: each member may be a single-residue **substitution**
/// (contributing its alternative residue to the merged insert) or a
/// single-residue **deletion** (contributing nothing) — e.g.
/// `p.[Arg76Ser;Cys77del]` → `p.Arg76_Cys77delinsSer`, and `delins.md:41`'s
/// `p.Cys28_Lys29delinsTrp` (del+sub). A run whose members are *all* deletions
/// has an empty insert, so it collapses to a pure range deletion
/// (`p.Arg76_Cys77del`) rather than an invalid empty `delins`.
///
/// This is the protein-axis counterpart of the nucleotide adjacency merge that
/// [`merge_consecutive_edits`] performs — the nucleotide path is positional on
/// a `Base` sequence and never fires on protein members, so the protein axis is
/// handled here instead.
///
/// Only **strictly adjacent** residues (position `n` and `n+1`) coalesce.
/// A gap of one or more unchanged residues keeps the members separate — the
/// spec pins that exact non-example (`protein/delins.md:63`: `p.[Ser44Arg;Trp46Arg]`
/// is *not* described as `p.Ser44_Trp46delinsArgLeuArg`).
///
/// Member **input order** does not affect the result (#1116): the members are
/// sorted into ascending residue order before runs are detected, so every
/// permutation of one member set coalesces identically. Reordering is
/// meaning-preserving here because every member has already been established to
/// be a plain single-residue substitution or deletion on one shared accession —
/// and the post-normalize display sort (#1098/#1101) puts cis members into
/// exactly that residue order anyway, so an allele that does *not* coalesce
/// still renders ascending. This is the protein-axis counterpart of #1103's
/// sort-before-merge on the nucleotide axis; without it a descending-order
/// allele fell through to the bracket form `protein/substitution.md:23` calls
/// "not correct".
///
/// A to-`Ter` (nonsense) member **bounds** coalescing rather than vetoing the
/// whole allele (#1125): a run merges only when it **ends at or before the
/// earliest `Ter`**. An unrelated upstream run still collapses, e.g.
/// `p.[Cys100Ser;Asp101Gly;Arg200Ter]` → `p.[Cys100_Asp101delinsSerGly;Arg200Ter]`.
///
/// A run that ends *exactly* at that `Ter` merges too (#1129): its insert then
/// carries the stop in final position, which the spec endorses
/// (`protein/delins.md:47`: `p.(Pro578_Lys579delinsLeuTer)`), so
/// `p.[Cys100Ser;Asp101Ter]` → `p.Cys100_Asp101delinsSerTer`. Only the two
/// shapes the spec actually forbids stay declined, and one comparison covers
/// both since `first_ter_pos` is the earliest stop:
///   - a run whose **first** member is the `Ter` — an *immediate* stop, which
///     must be a nonsense substitution (`protein/substitution.md:20`: not
///     `p.Cys5_Ser6delinsTerGluAsp`, but `p.Tyr4Ter`);
///   - a run carrying the `Ter` **interior** — residues after the stop
///     (`protein/delins.md:45`: not `delinsSerSerTerAlaAsp`).
///
/// A run sitting wholly downstream of an earlier stop is likewise left as
/// authored. Collapsing a nonsense change toward `p.<ref><pos>Ter` is a
/// separate rule.
///
/// Every other restriction is likewise scoped to the member or run it concerns,
/// so one out-of-scope member never leaves a well-formed run elsewhere in the
/// allele as the bracket `protein/substitution.md:23` calls "not correct"
/// (#1130). A member is barred from merging — and re-emitted verbatim — when:
///   - its edit is not a single-residue substitution or deletion (`dup` / `ins`
///     / `fs` / `ext` / multi-residue are out of scope for this narrow rule);
///     such a member is **opaque**, and it also bars any member it overlaps;
///   - it overlaps another member — two edits at one residue are the
///     contradictory same-residue case, not a delins.
///
/// A run additionally **splits** at any predicted `( )` ↔ certain boundary, so a
/// mixed stretch never fabricates a certainty its members do not jointly assert,
/// while each same-certainty stretch of ≥2 still merges.
///
/// Returns `Some` only when at least one run of ≥2 adjacent members is merged;
/// a fully-merged allele collapses to a bare [`HgvsVariant::Protein`], a
/// partial merge to a smaller [`HgvsVariant::Allele`]. Returns `None` — leave
/// the allele untouched — when no run merges, or when the allele as a whole is
/// out of scope:
///   - the phase is not cis (trans / unknown / mosaic / … are independent),
///   - a member is not a protein variant, or does not share the first member's
///     accession and gene symbol, or
///   - a member has an uncertain or range endpoint, leaving it unplaceable —
///     without a total order over the members, run detection is unsound, so
///     this one really does decline the whole allele.
pub(crate) fn coalesce_protein_adjacent_changes(
    allele: &crate::hgvs::variant::AlleleVariant,
) -> Option<HgvsVariant> {
    use crate::hgvs::edit::{AminoAcidSeq, ProteinEdit};
    use crate::hgvs::interval::ProtInterval;
    use crate::hgvs::location::{AminoAcid, ProtPos};
    use crate::hgvs::variant::ProteinVariant;

    if allele.phase != AllelePhase::Cis || allele.variants.len() < 2 {
        return None;
    }

    /// What a member contributes to a merged run.
    enum Change {
        /// A single-residue substitution; contributes `aa` to the insert.
        Substitution(AminoAcid),
        /// A single-residue deletion; contributes nothing to the insert.
        Deletion,
        /// Any other member — a `dup` / `ins` / `fs` / `ext` / multi-residue
        /// edit, i.e. an edit kind this narrow rule does not merge. It is
        /// ordered with the rest and re-emitted verbatim, but never joins a run
        /// (#1130).
        Opaque,
    }

    /// One member of the allele, reduced to what run detection needs.
    ///
    /// `lo`/`hi` are the member's occupied residue span (1-based, inclusive);
    /// a single-residue member has `lo == hi`. `source` is the member's index
    /// in `allele.variants`, kept so a member that ends up in no run is
    /// re-emitted verbatim even after the residue-order sort (#1116).
    struct Member {
        lo: u64,
        hi: u64,
        reference: AminoAcid,
        change: Change,
        predicted: bool,
        /// Set when this member overlaps another: a second edit at one residue
        /// (contradictory), or a residue covered by an opaque member's span.
        /// Blocked members are emitted verbatim and break any run through them
        /// (#1130).
        blocked: bool,
        source: usize,
    }

    let first = match &allele.variants[0] {
        HgvsVariant::Protein(pv) => pv,
        _ => return None,
    };
    let accession = first.accession.clone();
    let gene_symbol = first.gene_symbol.clone();

    let mut members: Vec<Member> = Vec::with_capacity(allele.variants.len());
    for (source, v) in allele.variants.iter().enumerate() {
        let HgvsVariant::Protein(pv) = v else {
            return None;
        };
        // Every member must share the first member's reference. A cis allele
        // states its accession once, so this normally holds by construction, but
        // the merged delins is stamped with `accession`/`gene_symbol` below —
        // guard so a mixed-accession member can never be silently re-labelled.
        if pv.accession != accession || pv.gene_symbol != gene_symbol {
            return None;
        }
        // Both endpoints must be certain single points so the member can be
        // ordered and its span known. An uncertain or range endpoint leaves the
        // member unplaceable, and run detection would be unsound without a
        // total order — so that one really does decline the whole allele.
        let (Some(Mu::Certain(s)), Some(Mu::Certain(e))) = (
            pv.loc_edit.location.start.as_single(),
            pv.loc_edit.location.end.as_single(),
        ) else {
            return None;
        };
        if e.number < s.number {
            return None;
        }
        let single_residue = s == e;
        let (change, predicted) = match &pv.loc_edit.edit {
            Mu::Certain(ProteinEdit::Substitution { alternative, .. }) if single_residue => {
                (Change::Substitution(*alternative), false)
            }
            Mu::Uncertain(ProteinEdit::Substitution { alternative, .. }) if single_residue => {
                (Change::Substitution(*alternative), true)
            }
            // A single-residue deletion contributes no residue to the merged
            // insert. A deletion carrying an explicit deleted-sequence or count
            // annotation still describes a single residue here, so it merges
            // the same way; a multi-residue `del` range falls through to
            // `Opaque` below.
            Mu::Certain(ProteinEdit::Deletion { .. }) if single_residue => {
                (Change::Deletion, false)
            }
            Mu::Uncertain(ProteinEdit::Deletion { .. }) if single_residue => {
                (Change::Deletion, true)
            }
            // Every other edit kind is out of scope for this narrow rule, but
            // only for itself — it no longer vetoes the allele (#1130).
            other => (Change::Opaque, matches!(other, Mu::Uncertain(_))),
        };
        members.push(Member {
            lo: s.number,
            hi: e.number,
            reference: s.aa,
            change,
            predicted,
            blocked: false,
            source,
        });
    }

    // Sort into ascending residue order so run detection — and therefore the
    // coalesced payload — is independent of the author's member order (#1116).
    // The sort is stable, so members sharing a start stay in input order.
    members.sort_by_key(|m| (m.lo, m.hi));

    // Mark every member that overlaps another. Two edits at one residue are a
    // contradiction rather than a delins, and a residue covered by an opaque
    // member's span cannot also be merged into a run — so both members of any
    // overlapping pair are emitted verbatim, and neither can join a run
    // (#1130). Quadratic in the member count, which is tiny for an allele.
    for i in 0..members.len() {
        for j in (i + 1)..members.len() {
            if members[i].hi >= members[j].lo {
                members[i].blocked = true;
                members[j].blocked = true;
            }
        }
    }

    // The earliest residue changed to `Ter` (a nonsense substitution), if any.
    // Only a run ending at or before this position may coalesce; every other
    // run is left as authored (see the to-`Ter` paragraph on this function).
    // Members are already in ascending residue order, so the first to-`Ter`
    // member is the earliest.
    let first_ter_pos = members
        .iter()
        .find(|m| matches!(m.change, Change::Substitution(AminoAcid::Ter)))
        .map(|m| m.lo);

    /// Whether `m` may take part in a merged run at all: it must be one of the
    /// two mergeable single-residue edit kinds and must not overlap another
    /// member (#1130).
    fn is_runnable(m: &Member) -> bool {
        !m.blocked && !matches!(m.change, Change::Opaque)
    }

    let mut merged: Vec<HgvsVariant> = Vec::new();
    let mut changed = false;
    let mut i = 0;
    while i < members.len() {
        // Grow the longest run of strictly-adjacent runnable members that also
        // agree on certainty. Splitting at a predicted/certain boundary keeps a
        // mixed stretch from fabricating a certainty its members do not jointly
        // assert, while still letting each same-certainty stretch of ≥2 merge
        // (#1130); before, the mixed case declined the whole allele.
        let mut j = i;
        if is_runnable(&members[i]) {
            while j + 1 < members.len()
                && is_runnable(&members[j + 1])
                && members[j + 1].lo == members[j].lo + 1
                && members[j + 1].predicted == members[j].predicted
            {
                j += 1;
            }
        }
        // A run may merge when it ends at or before the earliest stop.
        //
        // `<` is the #1125 case: the run lies wholly 5′ of the stop, which does
        // not constrain it. `==` is the #1129 case: the earliest `Ter` is the
        // run's own **last** member, so the merged insert carries the `Ter` in
        // final position — the spec's `p.(Pro578_Lys579delinsLeuTer)` form
        // (`protein/delins.md:47`).
        //
        // `>` covers both forbidden shapes at once, because `first_ter_pos` is
        // the *earliest* stop and run positions strictly ascend: a run reaching
        // past it either starts at that `Ter` (an immediate stop, which must be
        // a nonsense substitution — `protein/substitution.md:20`) or carries it
        // interior (residues after the stop — `protein/delins.md:45`), or else
        // sits wholly downstream of an earlier stop, which stays out of scope.
        let run_ends_at_or_before_any_ter = first_ter_pos.is_none_or(|ter| members[j].lo <= ter);
        if j > i && run_ends_at_or_before_any_ter {
            // A run of ≥2 strictly-adjacent members → one edit spanning
            // residues `i..=j`. The inserted sequence is each member's
            // alternative in order, with deletions contributing nothing.
            // Certainty is uniform across the run by construction (the walk
            // above splits on it), so the whole run is predicted iff its first
            // member is.
            let all_predicted = members[i].predicted;
            let loc = ProtInterval::new(
                ProtPos::new(members[i].reference, members[i].lo),
                ProtPos::new(members[j].reference, members[j].lo),
            );
            let inserted: Vec<AminoAcid> = members[i..=j]
                .iter()
                .filter_map(|m| match m.change {
                    Change::Substitution(aa) => Some(aa),
                    Change::Deletion | Change::Opaque => None,
                })
                .collect();
            let edit = if inserted.is_empty() {
                // Every member in the run is a deletion, so there is nothing to
                // insert: emit a pure range deletion, not an (invalid) empty
                // `delins`.
                ProteinEdit::Deletion {
                    sequence: None,
                    count: None,
                }
            } else {
                ProteinEdit::Delins {
                    sequence: AminoAcidSeq::new(inserted),
                }
            };
            let loc_edit = if all_predicted {
                LocEdit::new_predicted(loc, edit)
            } else {
                LocEdit::new(loc, edit)
            };
            merged.push(HgvsVariant::Protein(ProteinVariant {
                accession: accession.clone(),
                gene_symbol: gene_symbol.clone(),
                loc_edit,
            }));
            changed = true;
        } else {
            // A lone member (substitution or deletion) — or every member of a
            // declined run — keeps its original variant verbatim, indexed
            // through `source` since `members` is in residue order, not input
            // order.
            for m in &members[i..=j] {
                merged.push(allele.variants[m.source].clone());
            }
        }
        i = j + 1;
    }

    if !changed {
        return None;
    }
    // A full merge to a single member drops the allele wrapper — but only when
    // the allele carried no whole-bracket uncertainty, so the outer predicted
    // marker is not silently lost. In practice `allele.uncertain` is always
    // false here (the `p.(…;…)` whole-allele form does not parse into a live
    // `AlleleVariant`), so the bare return is what actually fires; the guard is
    // defensive future-proofing against that parser gap ever closing.
    if merged.len() == 1 && !allele.uncertain {
        merged.into_iter().next()
    } else {
        let mut coalesced = crate::hgvs::variant::AlleleVariant::new(merged, allele.phase);
        coalesced.uncertain = allele.uncertain;
        Some(HgvsVariant::Allele(coalesced))
    }
}

// ----------------------------------------------------------------------------
// Sequence-first canonicalization (#1229-#1235)
// ----------------------------------------------------------------------------

/// Longest reference window the sequence-first canonicalizer will fetch.
///
/// Bounds the cost of a group whose members are far apart; a wider group is
/// refused and falls back to the per-member pipeline.
const MAX_CANONICAL_WINDOW: i64 = 4096;

/// Padding either side of the changed interval, giving the 3'-shift room.
const CANONICAL_PAD: i64 = 128;

/// Longest **length-changing** block the canonicalizer will attempt to
/// partition.
///
/// Equal-length blocks are exempt: [`best_alignment`] returns their
/// position-wise pairing immediately, so there is no gap placement to search for
/// and nothing here to protect. They are bounded instead by
/// [`MAX_CANONICAL_WINDOW`], the widest window the canonicalizer will fetch at
/// all.
///
/// This is a cost bound, not a policy, and a block this long is a structural
/// event rather than a few nucleotide changes. It is deliberately far above the
/// size at which the separation rule (`delins.md:17`) is the interesting
/// question.
///
/// It does **not** bound a quadratic cost, despite what an earlier revision of
/// this comment claimed. A placement scores in O(1) because the score is
/// separable, so [`best_alignment`]'s search is linear; the one quadratic step —
/// ranking a tie by what it separates — carries its own and tighter bound in
/// [`MAX_TIE_BREAK_SWEEP`], which binds first at every length this constant
/// reaches.
///
/// It is *not* the `delins.md:44-47` "coincidental alignment" guard. That
/// concern — a large replacement whose interior only accidentally aligns, which
/// the spec keeps as one delins — is handled by [`separations_are_meaningful`]
/// and the [`changed_columns_of_pieces`] bound, not here.
///
/// Restricting the search to single-gap alignments is **not** that guard, and an
/// earlier revision of this comment claimed it was. The restriction stops the
/// aligner shredding a contiguous replacement by opening compensating gaps, which
/// is a real and necessary bound — the #1034/#1040/#182 "stays delins" cases
/// depend on it. But within one gap [`best_alignment`] still scores *every*
/// placement and keeps the highest, so it actively hunts for the position with
/// the most matches and will happily seize on a base that survives by
/// coincidence. Structure alone therefore cannot tell a real separation from an
/// accidental one; [`separations_are_meaningful`] is what draws that line.
///
/// Every length-changing regime now has that guard, so this is the single bound
/// for all of them (#1271). It used to be joined by a second, smaller bound for
/// net deletions, which `separations_are_meaningful` did not reach; extending
/// that rule retired it. See its doc comment for the measurement.
const MAX_SPLIT_BLOCK: usize = 1024;

/// Unchanged reference bases two pieces of a net insertion must be separated by
/// before the split between them is believed. See `separations_are_meaningful`.
///
/// `1`, per `general.md:34`: "two variants separated by one or more nucleotides
/// should be described individually and **not** as a 'delins'". One base is one
/// or more, so a split across one unchanged base is what the spec asks for, not
/// something to disbelieve.
///
/// It was `2`, which granted every axis `general.md:35`'s exception — "separated
/// by one nucleotide, **together affecting one amino acid**" — although that
/// exception is doubly conditioned and cannot reach a genomic axis at all. The
/// codon half is applied where it belongs, triplet-precisely, by
/// `apply_coding_codon_exception` *after* partitioning; see
/// [`MIN_SEPARATION_NO_FRAME`], which reached the same value from the same
/// reading for the sequence-first splitter.
///
/// Lowering it is only safe alongside [`split_buys_no_higher_priority_type`]:
/// #422 and #999's negative control are the *same shape* at this threshold —
/// net insertion, one coincidentally matched interior base — and the type of
/// the resulting members is the only thing that separates them.
///
/// It applies only where the block's net length change is at most
/// [`MAX_SINGLE_BASE_SEPARATION_CHANGE`]; beyond that the threshold rises to
/// [`RAISED_PIECE_SEPARATION`].
///
/// # Mutalyzer does not have a separation rule, so do not read one out of it
///
/// Recorded here because the question — "maybe the reference implementation
/// knows something the recommendations do not" — recurs on every divergence
/// at this threshold, and answering it costs a source dig each time.
///
/// Mutalyzer's `/api/normalize` re-derives a description from the *mutated
/// sequence*: `Description.normalize()` -> `mutate()` -> `extract()` ->
/// `describe_dna(reference, observed)`. `extract()` minimizes a **weighted
/// description length**; the weights are constants in the extractor, dated
/// 2014 — seven years before SVD-WG010, the (rejected) proposal that would
/// have restated this rule positionally. The strings "separated by", "amino
/// acid" and "SVD-WG" occur nowhere in its source, its tests or its commit
/// history. A stratified 72-row probe found the weight model predicts its
/// answers on **91.7%** (66/72), against **79.2%** for the best achievable
/// separation threshold — i.e. cost, not distance, is what it is computing.
///
/// Two consequences. A Mutalyzer/ferro disagreement about separation is two
/// objectives meeting, not evidence about `general.md:34`. And where the spec
/// speaks plainly Mutalyzer is measurably wrong: it splits
/// `LRG_199t1:c.145_147delinsTGG` into `c.[145C>T;147C>G]`, the form
/// `DNA/delins.md:42` says in as many words "is not correct".
///
/// Ferro's split at separation 2 is therefore a **deliberate** divergence:
/// `general.md:34` splits there (and the rejected SVD-WG010's own example at
/// `:45` splits `c.[235A>T;238G>T]` at exactly two nucleotides) while Mutalyzer
/// merges on 23 of 25 measured rows. Which authority governs that choice is
/// itself unsettled — see `rulings[adjudication-precedence-order]` in
/// `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`.
const MIN_PIECE_SEPARATION: usize = 1;

/// The separation [`separations_are_meaningful`] requires once the block's net
/// length change exceeds [`MAX_SINGLE_BASE_SEPARATION_CHANGE`].
///
/// Named rather than written inline because it is a policy choice doing
/// semantic work, not an arithmetic constant: beyond that much change a single
/// matched base is not believable as structure, so two are required before a
/// split is believed. `2` is the value the threshold had for every block before
/// [`MIN_PIECE_SEPARATION`] was lowered to 1, so this is the old behaviour
/// retained where the evidence for the new one runs out — not a second
/// calibration.
const RAISED_PIECE_SEPARATION: usize = 2;

/// Largest block, in aligned positions, over which a tied gap placement is
/// broken by what it separates rather than by being 5'-most.
///
/// The score for one placement is O(1) — that is what makes `best_alignment`'s
/// search linear — but ranking a *tie* needs the placement's member count, which
/// costs a pass over the block. In a repeat tract every placement ties, so the
/// tie-break is O(n²) there, which at [`MAX_SPLIT_BLOCK`] would be ~1e6 column
/// operations for one call.
///
/// **This is a bound on coverage, not a silent cap.** Above it the 5'-most
/// placement wins as before, so a tied block longer than this keeps whatever
/// description the old rule gave. Nothing in the corpus is near it — #1262's
/// block is 3 aligned positions — and the bound exists so a pathological
/// homopolymer cannot make one call quadratic in the window size.
const MAX_TIE_BREAK_SWEEP: usize = 256;

/// Placements [`two_gap_insertion_alignment`] will examine before giving up.
///
/// The search sweeps `a` over the two blocks' common prefix and `b` over their
/// common suffix. Both are **zero on a trimmed block**, which is what every
/// caller passes today, so the sweep costs a handful of placements — but the
/// point of generalising that search was to stop depending on the caller's trim,
/// and leaving the *cost* dependent on it would reintroduce the same coupling.
///
/// Each placement costs an `O(block)` slice comparison. Measured worst case (no
/// solution exists, so the sweep runs to exhaustion) at [`MAX_SPLIT_BLOCK`] with
/// a 6-base insertion:
///
/// ```text
/// common affixes 0, 0        5 placements     (trimmed — today)
/// common affixes 1, 1       20 placements     (one flank base restored per side)
/// common affixes 512, 512    1 315 840        (~1.3e9 byte comparisons)
/// ```
///
/// # Why this declines rather than narrowing the sweep
///
/// An earlier revision clamped the affix *lengths* to a small constant instead.
/// That is not a cost bound, it is a **silent change of answer**, and it was
/// measured doing exactly that: `AAAAAAA -> AACACAAAA` (#1260's untrimmed
/// window) decomposes only at `a = 2, b = 3`, which needs a common suffix of 4.
/// Clamping to 2 excluded the sole solution, the function returned `None`, and
/// the block fell back to a spanning `delins`. It went unnoticed because on that
/// row the single-gap search happens to reach the same member count by another
/// route — so the clamp was wrong in a way no test could see until something
/// else moved.
///
/// Exceeding this budget therefore abandons the whole two-gap search and lets
/// [`best_alignment`] keep its single-gap winner — the answer it gave before this
/// function existed — rather than returning a decomposition found under a
/// narrower search than the one documented. The budget is sized far above the
/// worked cases above; nothing in the corpus approaches it.
const MAX_TWO_GAP_PLACEMENTS: usize = 65_536;

/// Largest alignment grid the sequence-first splitter will build, in cells.
///
/// `AlignmentDag::build` is `Θ(n·m)` in **space** as well as time, and it
/// allocates four grids of that size: two `u32` edit grids, the on-path flags,
/// and the adjacency mask. Only the reference side is bounded upstream
/// ([`MAX_CANONICAL_WINDOW`]); the alternate side is whatever the members
/// produce, and a duplication over the whole window doubles it. A 4096 x 8192
/// grid is 33.5 M cells, or roughly 340 MB across the four — enough to matter,
/// and reachable from a single description rather than from a pathological
/// corpus.
///
/// Derived from [`MAX_CANONICAL_WINDOW`] rather than picked: the budget is a
/// *square* grid at that bound, so it admits any block whose alternate side is
/// within the same limit already imposed on the reference side, and refuses one
/// where the alternate side runs past it. The 4096 x 8192 duplication above is
/// 33.5 M cells against this 16.8 M, so it is declined.
///
/// The honest cost of stating it this way: at the very top of the reference range
/// an alternate side even slightly longer than the reference is declined too.
/// That is a cost bound, not a policy — above it
/// [`partition_block_sequence_first`] returns `None` and the caller keeps the
/// per-member result, the same fallback every other bound here takes.
const MAX_SEQFIRST_GRID_CELLS: usize =
    (MAX_CANONICAL_WINDOW as usize + 1) * (MAX_CANONICAL_WINDOW as usize + 1);

/// Largest **net length change**, in bases, for which one unchanged base counts
/// as separation.
///
/// A single matched base is evidence of structure when the surrounding change is
/// small and evidence of nothing when it is large, where a one- or two-base run
/// recurs by chance. The Mutalyzer conformance corpus is what forces this: a 6 nt
/// block replaced by a 21 nt payload
/// (`NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC` and nine siblings) has
/// interior bases that happen to match, and the oracle keeps it as one `delins`.
/// Without a bound, a threshold of 1 shreds it into
/// `g.[5207_5209delinsGTC;5211C>T;5213_5214insGCTCATTATCTGGCT]` — mixed member
/// types, so [`split_buys_no_higher_priority_type`] does not reach it either.
///
/// # Why the change and not the block length
///
/// Keying on block length looks equivalent and is not, because the block is
/// about to be handed the *supremal* extent rather than the trimmed one. A
/// variant rolled out over a homopolymer has a long block and a tiny change —
/// #1260's window is 7 nt of reference for 2 inserted bases — while the
/// Mutalyzer block is long *and* changes a lot. Length cannot tell those apart;
/// the net change can, and it is the quantity the belief is actually about.
///
/// **The value is under-determined and the bounds are what is measured.** The
/// corpus constrains it only to `[2, 15)`: #1260's window changes 2 bases and
/// must split at one unchanged base, the Mutalyzer block changes 15 and must
/// not. Every value in between scores identically, so `4` is a choice inside a
/// measured window, not a calibrated constant.
///
/// **No threshold at all satisfies #1157 and #1232 simultaneously**, and that
/// is a measured non-existence rather than a value nobody has found yet: the
/// two demand that a 7-nt block stay whole while a 5-nt block splits, which no
/// monotone bound on block length or net change can deliver. That pair is
/// resolved instead by adopting the split as the default and carrying
/// `LRG_199` and #422 as *named* exceptions — see
/// [`split_buys_no_higher_priority_type`] and the `issue_422_*` tests. Do not
/// go looking for the value that reconciles them.
const MAX_SINGLE_BASE_SEPARATION_CHANGE: usize = 4;

/// [`partition_block_sequence_first`]'s separation threshold on an axis with no
/// reading frame (`g.`, `m.`, `n.`, and `r.` on a non-coding transcript).
///
/// `1`, not `seqfirst::MIN_SEPARATION`'s `2`. `general.md:34`'s general rule is
/// "two variants separated by one or more nucleotides should be described
/// individually" — one or more unchanged bases between two runs of change means
/// they split, so runs merge only when separated by *fewer than one*, i.e. zero,
/// unchanged bases. `general.md:35`'s exception — "separated by one nucleotide,
/// together affecting one amino acid" — is what licenses merging across exactly
/// one unchanged base, and it is scoped to "one amino acid": a reading frame,
/// which `AxisFrame::reading_frame` already isolates for
/// `apply_coding_codon_exception`. Applying `seqfirst::MIN_SEPARATION` (`2`)
/// unconditionally therefore grants every axis the coding exception, which is
/// what the `FERRO_SEQFIRST_SHADOW` audit's class `B2` measured: 3 088 blocks
/// where the sequence-first splitter merged across an unchanged base the live
/// splitter, correctly, kept split. Pinned by
/// `sequence_first_split_axis_separation_matches_general_rule`.
const MIN_SEPARATION_NO_FRAME: u32 = 1;

/// One derived edit: a maximal run of change over the reference window, as
/// 0-based half-open offsets into that window plus its replacement bases.
///
/// A pure insertion has `ref_start == ref_end` (a zero-width span at the gap).
#[derive(Debug, Clone, PartialEq, Eq)]
struct Piece {
    ref_start: usize,
    ref_end: usize,
    alt: Vec<u8>,
}

impl Piece {
    /// Whether this piece is a pure deletion or a pure insertion.
    ///
    /// Only these may be 3'-shifted. `shuffle`'s rotation invariant — advance
    /// one position when the base leaving the span equals the base entering it
    /// — holds only when the span is *either* removed or inserted, never both.
    /// Rotating a `delins` changes the resulting sequence, and the spec scopes
    /// the rule the same way: "for deletions, duplications, and insertions, the
    /// most 3' position possible is arbitrarily assigned to have been changed"
    /// (`checklist.md:37`).
    fn is_pure_indel(&self) -> bool {
        let ref_len = self.ref_end - self.ref_start;
        (self.alt.is_empty() && ref_len > 0) || (ref_len == 0 && !self.alt.is_empty())
    }

    /// Whether this piece is a pure insertion sitting immediately outside one
    /// end of the fetched window.
    ///
    /// A cheap necessary condition for the terminal-insertion clamp: only such
    /// a piece can be resting on a *sequence* bound, so only for such a piece is
    /// it worth asking the provider whether the window edge is the sequence's
    /// edge (see [`window_sequence_ends`]).
    fn rests_on_a_window_edge(&self, ref_bytes: &[u8]) -> bool {
        self.ref_start == self.ref_end
            && !self.alt.is_empty()
            && (self.ref_start == 0 || self.ref_start == ref_bytes.len())
    }
}

/// Whether to run both block splitters and report their disagreements.
///
/// A **development-only audit switch**, enabled by `FERRO_SEQFIRST_SHADOW=1`.
/// Read once and cached, like the idempotency gate in [`crate::normalize`], so a
/// disabled run pays one relaxed atomic load per canonicalization. When on, the
/// sequence-first splitter runs *in addition to* the live one and every
/// comparison is reported; it never changes what is returned.
///
/// Reports go through `log::debug!` rather than `eprintln!`, so ordinary log
/// filtering applies — the switch fires once per canonicalized block, which over
/// a whole suite run is tens of thousands of lines that no level or target filter
/// could reach when they were written straight to stderr. Capture them with
/// `RUST_LOG=ferro_hgvs::normalize::merge=debug`.
///
/// # Promoting it is a churn event, not a drop-in — measured
///
/// The sequence-first splitter is more principled (confluence is structural
/// rather than repaired), which makes it tempting to flip on. Measured against
/// the live `partition_block` on the 14 blocks drawn from the #1419/#1420/#1421
/// examples, it **disagrees on 13**, and on several it yields a *third*
/// representation matching neither side of the current non-confluence. So
/// promoting it does not converge the disputed pairs onto one of their existing
/// forms; it re-buckets them onto new ones, which for the downstream consumer
/// is a re-normalization of the whole stored library.
///
/// Note also that "more principled" is not a spec argument here: no minimality
/// or confluence principle appears in the recommendations at all — see
/// `changed_columns_of_pieces` and `background/basics.md:38`. Measure the blast
/// radius on shipped forms before proposing the flip.
fn seqfirst_shadow_enabled() -> bool {
    static ENABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var("FERRO_SEQFIRST_SHADOW").as_deref() == Ok("1"))
}

/// Which block partitioner [`canonicalize_from_sequence`] cuts with.
///
/// Three rules exist and they are not variations on one — see
/// [`crate::normalize::seqfirst::partition::CanonicalAlignment`] for what
/// separates the latter two.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PartitionRule {
    /// [`partition_block`]: a single-gap alignment search plus its two narrow
    /// escapes. The shipped rule, and the default.
    Live,
    /// [`partition_block_sequence_first`]: cut at the steps common to every
    /// minimal alignment (`mutalyzer-algebra`'s `local_supremal`).
    Shadow,
    /// [`partition_block_canonical`]: the member-count-minimal minimal alignment
    /// (`mutalyzer-algebra`'s `canonical`).
    Canonical,
    /// `Canonical`, plus the terminal `delins.md:44-47` re-spelling applied
    /// after the downstream passes — see
    /// [`coalesce_payload_alignment_split`].
    CanonicalCoalesced,
}

/// Read a [`PartitionRule`] from what `FERRO_PARTITION` was set to.
///
/// Unset, empty, or unrecognised all yield [`PartitionRule::Live`], so no
/// environment can make ferro emit something other than its shipped
/// representation by accident. An unrecognised value is warned about rather than
/// rejected: this is a measurement knob, and a typo that silently selected a
/// different partitioner would poison a bake-off far more quietly than one that
/// silently ran the default and said so.
///
/// Split out from [`partition_rule`] so the mapping is testable without an
/// environment variable — the cached read below can only ever be exercised once
/// per process.
fn partition_rule_from_env(value: Option<&str>) -> PartitionRule {
    match value {
        None | Some("") | Some("live") => PartitionRule::Live,
        Some("shadow") => PartitionRule::Shadow,
        Some("canonical") => PartitionRule::Canonical,
        Some("canonical-coalesced") => PartitionRule::CanonicalCoalesced,
        Some(other) => {
            log::warn!(
                "FERRO_PARTITION={other} is not one of \
                 live|shadow|canonical|canonical-coalesced; using live"
            );
            PartitionRule::Live
        }
    }
}

/// Which partitioner this process cuts blocks with, from `FERRO_PARTITION`.
///
/// A **development-only bake-off switch**. Read once and cached, like
/// [`seqfirst_shadow_enabled`], so the default path pays one relaxed atomic load
/// per canonicalized block and nothing else; with the variable unset the result
/// is byte-identical to having no switch at all.
fn partition_rule() -> PartitionRule {
    static RULE: std::sync::OnceLock<PartitionRule> = std::sync::OnceLock::new();
    *RULE.get_or_init(|| partition_rule_from_env(std::env::var("FERRO_PARTITION").ok().as_deref()))
}

/// The unchanged-base separation threshold an axis merges runs below.
///
/// A reading-frame axis gets `general.md:35`'s coding one-amino-acid exception
/// ([`crate::normalize::seqfirst::MIN_SEPARATION`], 2); every other axis gets the
/// general rule ([`MIN_SEPARATION_NO_FRAME`], 1). See the latter for why applying
/// the exception everywhere is wrong.
fn axis_min_separation(reading_frame: bool) -> u32 {
    if reading_frame {
        crate::normalize::seqfirst::MIN_SEPARATION
    } else {
        MIN_SEPARATION_NO_FRAME
    }
}

/// A piece list rendered on one line with `alt` as text rather than bytes.
///
/// Only used by the `FERRO_SEQFIRST_SHADOW` report, whose whole value is being
/// greppable: `Piece`'s derived `Debug` prints `alt` as `[65, 67, …]`, which
/// makes a report of a 40-base block unreadable and long.
struct DebugPieces<'a>(&'a [Piece]);

impl std::fmt::Debug for DebugPieces<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("[")?;
        for (i, piece) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_str(",")?;
            }
            write!(
                f,
                "{}..{}>{}",
                piece.ref_start,
                piece.ref_end,
                String::from_utf8_lossy(&piece.alt)
            )?;
        }
        f.write_str("]")
    }
}

/// One column of an alignment between the reference and result blocks.
///
/// `None` on either side is a gap. Both-`Some` columns are matches when the
/// bases agree and substitutions when they do not.
type Column = (Option<usize>, Option<usize>);

/// Re-derive a cis allele (or a lone member) from the sequence it produces.
///
/// This is the sequence-first canonicalization of #1235. Rather than
/// normalizing each member in isolation — which makes the result depend on how
/// the variant was spelled and lets one member's 3'-shift run over a sibling —
/// it applies the whole group to the reference, re-derives the minimal edit set
/// from the (reference, result) pair, and partitions and types that set once,
/// globally. Two encodings of one variant therefore converge, and members
/// cannot overlap by construction.
///
/// "Cannot overlap" is a statement about the *output*, and it is narrower than
/// it sounds about the input. The geometry this path applies is **interbase**,
/// via [`apply_edits_to_window`]: a zero-length junction and a span that merely
/// touches it are flush rather than overlapping, so `g.[24dup;24C>G]` is
/// admitted here even though both members name base 24.
/// `overlap::detect_overlap_conflicts` asks the other question — coincident
/// spans in **HGVS-coordinate** space — and reports that pair as
/// `OverlapConflictingEdits` / W5002. The two are not redundant and the
/// coincident-span check is not implied by this one, which is why the caller
/// runs `detect_overlap_conflicts` *before* reaching here rather than relying on
/// this function to refuse (#1307).
///
/// The spec is explicit that the description follows from the sequence and not
/// from the input spelling: `delins.md:86-89` *deleted* its former carve-out
/// permitting a two-member spelling for a variant "likely a combination of two
/// other variants".
///
/// Returns `None` whenever the path does not apply — no reference window, a
/// repeat/uncertain/protein member, a mixed axis, and so on — and `Some` only
/// when the re-derivation both succeeded and differs from the input. Refusal is
/// the established convention here (see `collapse_overlapping_cis_edits`).
///
/// `variants` is borrowed and the result optional so that the far more common
/// refusal costs the caller nothing: by value it would have to clone the whole
/// member list — a thousand members for a large allele — just to have it handed
/// straight back.
pub(crate) fn canonicalize_from_sequence<P: ReferenceProvider>(
    variants: &[HgvsVariant],
    phase: AllelePhase,
    provider: &P,
    direction: ShuffleDirection,
) -> Option<Vec<HgvsVariant>> {
    if variants.is_empty() || (variants.len() > 1 && phase != AllelePhase::Cis) {
        return None;
    }
    let kind = cis_kind_of(&variants[0])?;
    // Every axis, genomic and transcript. #1237 gated this to `g.`/`m.` because
    // the transcript axes layer extra rules on a raw re-derivation; each is now
    // either implemented here or shown to be already covered:
    //
    // * the coding one-amino-acid exception (`general.md:35`, still current law
    //   — SVD-WG010 was rejected) — `apply_coding_codon_exception` below;
    // * the `r.` `U`/`T` alphabet — `canonical_base_byte` folds the window into
    //   one alphabet before any comparison;
    // * the CDS/UTR and transcript-end insertion clamps (#1202, #1205, #1207,
    //   #1209, #383, #387) and the exon-junction exception to the 3' rule
    //   (`general.md:44`) — the body region check in `collect_canonical_edits`
    //   keeps every member inside the positive body, and the terminal-insertion
    //   clamp is applied here, by `boundary_delins_anchor`. (It used to be
    //   *deferred* here, by refusing every derivation that collapsed to one pure
    //   insertion; that refusal is gone — see where it stood, below.)
    let body = body_region(kind);
    let (template_accession, _, _, _, _) = cis_axis_parts(&variants[0], kind)?;
    let template_accession = template_accession.clone();

    let mut edits = collect_canonical_edits(variants, kind, body, &template_accession)?;

    // A **lone** repeat is not what the lockout is about, and re-deriving one
    // can only cost. Its whole change is a single contiguous tract, so there is
    // no partition decision for this pass to make — the one thing it could do
    // is hand back the `ins`/`del` that denotes the same bases, trading a
    // spelling the spec explicitly leaves open (`open-issues.md:88`) for
    // another. Sometimes the per-member renderer promotes that straight back
    // (`274_275insAA` → `274A[3]`) and the alternation settles; sometimes it
    // does not — a `B2` contraction re-derives as a plain deletion and
    // `g.258CAG[2]` came back as `g.258_263del`. The value is in a *group*,
    // where the tract's bases can merge with a sibling's, so that is where
    // lowering is allowed to run.
    if variants.len() < 2 && edits.iter().any(|e| matches!(e, GEdit::Repeat { .. })) {
        return None;
    }

    // Window: the union of every member's footprint, padded for the 3'-shift.
    let (c_lo, c_hi) = edit_span_union(&edits)?;
    let w_lo = (c_lo - CANONICAL_PAD).max(1);
    let w_hi = c_hi + CANONICAL_PAD;
    if w_hi - w_lo + 1 > MAX_CANONICAL_WINDOW {
        return None;
    }

    let frame = axis_frame(kind, &template_accession, provider)?;
    let accession = template_accession.transcript_accession();

    // `general.md:44` exempts deletions and duplications around an exon/exon
    // junction from the 3' rule on the transcript axes, and the per-member
    // pipeline honours it — so the members reaching this function are *already*
    // clamped. Re-deriving across a junction throws that away: `c.[10del;13del]`
    // on a transcript whose G-run spans c.9-c.33 arrives here as the correctly
    // clamped `[c.11del;c.26del]` and comes back as `c.32_33del`, both members
    // carried past their junctions and merged into one (#1450).
    //
    // The constraint belongs on the **input span**, here, and nowhere later. By
    // the time `trim_common_flanks` has run, the two deletions have collapsed
    // into a single run and which exon each came from is no longer recoverable
    // from the resulting sequence; the merged piece then sits wholly inside one
    // exon, so no check on the block, the partition or the shift can see that a
    // junction was crossed. (`shift_pieces` in particular moves nothing for this
    // shape, so a clamp there is unreachable — measured on #1450.)
    //
    // Declining rather than deriving per exon: the per-member result that
    // survives a refusal *is* the expected answer, so refusing costs this shape
    // nothing. Segmenting the block at exon edges would additionally buy
    // confluence for cross-junction alleles; that is materially larger and is
    // left open.
    if crosses_exon_junction(kind, &accession, provider, &frame, c_lo, c_hi) {
        return None;
    }
    let start0 = w_lo + frame.delta - 1;
    let end0 = w_hi + frame.delta;
    if start0 < 0 {
        return None;
    }
    // A short read at the 3' end is normal near the sequence end; retry with the
    // window the provider can actually serve so terminal variants still resolve.
    let ref_bytes = fetch_canonical_window(provider, &accession, start0, end0)?;
    // Compare in one alphabet. On the `r.` axis the submitted bases are RNA
    // (`u`) while the transcript reference is DNA (`T`), so a position holding
    // the same nucleotide would otherwise read as changed and corrupt the whole
    // partition. Folding `U` to `T` is safe to do unconditionally, and costs
    // nothing on the DNA axes. Emitting `T` is likewise correct for `r.`: its
    // formatter lower-cases and maps `T` to `u` when it renders an alt
    // sequence, the same assumption `apply_canonical_split` already relies on.
    // (`fetch_canonical_window` has already upper-cased, so this only folds.)
    let ref_bytes: Vec<u8> = ref_bytes.iter().map(|b| canonical_base_byte(*b)).collect();

    // A repeat member names its tract implicitly, so it is only now — with the
    // window in hand — that it denotes anything. Every `GEdit::Repeat` becomes
    // the `ins`/`del` it means, and nothing past this line can see one (#1296).
    lower_repeat_edits(&mut edits, &ref_bytes, w_lo)?;
    let edits = edits;

    // A member that misstates its reference base is a reference mismatch, not a
    // canonicalization problem: strict mode must still reject it and lenient
    // mode must still warn. Both live in the per-member pipeline, so refuse and
    // let it run (#1052, #1097).
    if !stated_reference_bases_match(variants, kind, &ref_bytes, w_lo) {
        return None;
    }
    let result = apply_edits_to_window(&edits, &ref_bytes, w_lo)?;

    // Trim to the minimal changed block. Trimming first is what makes the
    // window choice canonical: `[4_6del;8A>T]` and `[5_7del;8A>T]` span
    // different windows but trim to the same block, which is precisely why they
    // must converge (#1234).
    let (lo, hi_ref, hi_alt) = trim_common_flanks(&ref_bytes, &result);
    if lo == hi_ref && lo == hi_alt {
        return None; // no net change; leave it to the existing pipeline.
    }

    // Which rule cuts the block. `FERRO_PARTITION` unset — the only configuration
    // that ships — takes the `Live` arm, so this is `partition_block` and nothing
    // else. The other two arms fall back to it when their splitter declines,
    // rather than abandoning the canonicalization, so a bake-off run measures the
    // partitioners and not their decline rates.
    let mut pieces = match partition_rule() {
        PartitionRule::Live => partition_block(&ref_bytes[lo..hi_ref], &result[lo..hi_alt]),
        PartitionRule::Shadow => partition_block_sequence_first(
            &ref_bytes[lo..hi_ref],
            &result[lo..hi_alt],
            axis_min_separation(frame.reading_frame),
        )
        .unwrap_or_else(|| partition_block(&ref_bytes[lo..hi_ref], &result[lo..hi_alt])),
        PartitionRule::Canonical => {
            partition_block_canonical(&ref_bytes[lo..hi_ref], &result[lo..hi_alt])
                .unwrap_or_else(|| partition_block(&ref_bytes[lo..hi_ref], &result[lo..hi_alt]))
        }
        PartitionRule::CanonicalCoalesced => {
            // Identical to `Canonical` here on purpose. The `delins.md:44-47`
            // merge is applied *after* the downstream passes below, not as part
            // of partitioning — see `coalesce_payload_alignment_split`.
            partition_block_canonical(&ref_bytes[lo..hi_ref], &result[lo..hi_alt])
                .unwrap_or_else(|| partition_block(&ref_bytes[lo..hi_ref], &result[lo..hi_alt]))
        }
    };
    // Shadow the sequence-first splitter on the very same trimmed block, before
    // any 3'-shift or coalescing, so a reported disagreement is a disagreement
    // about the *split* and not about a downstream pass. Never affects `pieces`.
    if seqfirst_shadow_enabled() {
        // The audit's live baseline, computed independently of `partition_rule()`.
        //
        // It cannot be `pieces`: since `FERRO_PARTITION` was added, `pieces` is
        // whatever rule the environment selected, so
        // `FERRO_SEQFIRST_SHADOW=1 FERRO_PARTITION=shadow` would compare the
        // sequence-first splitter with *itself* and report `SEQFIRST_AGREE` on
        // every block — a full denominator and zero disagreements, which reads
        // as the strongest possible result and is produced by comparing a thing
        // with itself. That is the precise failure the denominator below exists
        // to rule out. `canonical` and `canonical-coalesced` would not be
        // vacuous but would label their output `live=`, which is worse than
        // useless in a bake-off log.
        //
        // The extra `partition_block` call is paid only when the shadow audit is
        // on, which is already a measurement-only path.
        let live = partition_block(&ref_bytes[lo..hi_ref], &result[lo..hi_alt]);
        // Per-axis separation threshold (see `axis_min_separation`) — the same
        // distinction `apply_coding_codon_exception` uses below, via the same
        // `AxisFrame`.
        let min_separation = axis_min_separation(frame.reading_frame);
        let block = &ref_bytes[lo..hi_ref];
        let mut shadow = partition_block_sequence_first(block, &result[lo..hi_alt], min_separation);
        // `min_separation` alone captures only the distance half of
        // `general.md:35`'s exception; this captures the codon half, splitting
        // back apart any merge the exception does not actually license. `lo`
        // is added because `block`'s offsets (unlike `pieces`' at this point)
        // have not yet been shifted into the untrimmed window `w_lo` is
        // relative to.
        if let Some(shadow_pieces) = shadow.as_mut() {
            split_codon_incompatible_triplets(
                shadow_pieces,
                frame.reading_frame,
                w_lo + lo as i64,
                block,
            );
        }
        // Both sides are also compared after `shrink_pieces_to_differences`,
        // because a split that only differs in how wide its members are is a
        // difference the rendering stage removes. Reported as a third outcome
        // rather than folded into either of the other two, so one run yields
        // both the raw disagreement count and the count that survives
        // minimisation. The shrink runs here before the 3'-shift, whereas the
        // returned `pieces` are shrunk after it — see
        // `shrinking_before_and_after_the_three_prime_shift_agree`.
        let mut min_live = live.clone();
        shrink_pieces_to_differences(&mut min_live, block);
        let mut min_shadow = shadow.clone();
        if let Some(min_shadow) = min_shadow.as_mut() {
            shrink_pieces_to_differences(min_shadow, block);
        }
        if shadow.as_ref() == Some(&live) {
            // The denominator. Without it, `grep -c SEQFIRST_SHADOW` returning 0
            // cannot be told apart from a shadow that never ran — the failure
            // mode that would make an audit of the disagreements vacuous.
            log::debug!("SEQFIRST_AGREE");
        } else if min_shadow.as_ref() == Some(&min_live) {
            log::debug!(
                "SEQFIRST_MINAGREE ref={} alt={} live={:?} shadow={:?}",
                String::from_utf8_lossy(block),
                String::from_utf8_lossy(&result[lo..hi_alt]),
                DebugPieces(&live),
                shadow.as_deref().map(DebugPieces),
            );
        } else {
            log::debug!(
                // Space-separated, and no space inside a piece list, on purpose:
                // `cargo nextest` strips tab characters when it re-emits captured
                // test output (`--success-output immediate`), so a tab-delimited
                // report is unparseable in the one run mode that covers the whole
                // suite. Verified by `od -c` on both run modes.
                "SEQFIRST_SHADOW ref={} alt={} live={:?} shadow={:?} min_live={:?} min_shadow={:?}",
                String::from_utf8_lossy(block),
                String::from_utf8_lossy(&result[lo..hi_alt]),
                DebugPieces(&live),
                shadow.as_deref().map(DebugPieces),
                DebugPieces(&min_live),
                min_shadow.as_deref().map(DebugPieces),
            );
        }
    }
    for piece in &mut pieces {
        piece.ref_start += lo;
        piece.ref_end += lo;
    }
    shift_pieces(&mut pieces, &ref_bytes, direction);
    coalesce_adjacent_pieces(&mut pieces);
    // Partitioning decides where the members are; this decides how wide each
    // one is spelled, and the two are not the same question — see
    // `shrink_pieces_to_differences`. Placed before the weight bound below so
    // that the bound judges the pieces that will actually be rendered.
    shrink_pieces_to_differences(&mut pieces, &ref_bytes);

    // A canonicalization may re-partition and re-type the change. It may not
    // describe *more* change than the input already did.
    //
    // `best_alignment` considers only single-gap alignments — one contiguous
    // indel plus substitutions — because letting it place compensating gaps lets
    // it manufacture matches from coincidence and shred a genuinely contiguous
    // replacement (#1034/#1040/#182). That restriction is sound only while the
    // block *has* a single-gap explanation as economical as the input's own.
    // When two members' length changes do not cancel within one gap — a `+3`
    // insertion and a `-3` deletion 11 nt apart, or simply two separate
    // deletions — the lone gap parks at one end and the remaining columns pair
    // up position-wise, so the offset between them reads as a run of
    // substitutions. The derived form then marks more columns changed than the
    // input did (12 against 7 for
    // `NG_012337.1(NM_012459.2):c.[10del;23_25del;36_37insAAT]`) and merges
    // across ten bases the input left unchanged (`general.md:34`). That is a
    // worse description, not a canonical one, so refuse and let the per-member
    // pipeline — which never re-aligns across members — answer.
    //
    // The bound is one-sided and self-limiting: it can never refuse a
    // single-member input, because the pieces of a single-gap alignment
    // partition that alignment's columns, so `sum(max(ref_i, alt_i))` is at most
    // `max(ref, alt)`; the 3'-shift preserves each piece's weight and coalescing
    // can only lower it. Every "stays delins" case (#1034, #1040, #182, #422) is
    // single-member and therefore untouched.
    //
    // Strict on purpose: #1229, #1231, #1233 and #1234 all sit exactly at
    // equality, and equality is where prioritisation (`general.md:56`) may
    // legitimately prefer the re-derived shape.
    if changed_columns_of_pieces(&pieces) > changed_columns_of_edits(&edits) {
        return None;
    }

    // `delins.md:44-47`, applied AFTER the weight bound, for exactly the reason
    // the note on `coalesce_whole_block_inversion` below gives: this pass widens,
    // and "a widening judged against the *input's* weight is accepted for one
    // spelling of a variant and refused for another."
    //
    // Measured, because it was first placed before the bound: judged there it
    // costs 427 converged classes per direction over the 11,272-class corpus
    // (5': 8,387 -> 7,960), with 427 regressions and 0 gains. Every one is a
    // pair whose merged form is accepted for the spelling that was already a
    // single `delins` -- whose input weight is the wide one -- and refused for
    // the spelling that arrived as separate members. The re-spelling itself is
    // spelling-independent; the bound it was being judged against is not.
    // That input-relative comparand is #1440.
    //
    // Placement relative to the 3'-shift and the width minimisation is NOT the
    // mechanism: moving the call across those produced byte-identical censuses.
    // Only the bound matters.
    if partition_rule() == PartitionRule::CanonicalCoalesced {
        coalesce_payload_alignment_split(&mut pieces, &ref_bytes);
    }

    // Applied *after* the weight bound for the same reason the codon exception
    // below is: it widens, and a widening judged against the *input's* weight is
    // accepted for one spelling of a variant and refused for another. See
    // `coalesce_whole_block_inversion`'s own doc for the worked case.
    //
    // It sits *before* `apply_coding_codon_exception`, and on a `c.`/`n.` axis
    // the two do reach for the same pieces: a 5 nt whole-block inversion whose
    // changed columns fall at offsets 0/2/4 inside one codon is exactly the
    // `[Sub@p; Identity@p+1; Sub@p+2]` triplet the codon exception merges
    // (`general.md:35`). Stub this call and that is what happens —
    // `c.[1G>T;3T>A;5A>C]` comes out as `c.[1_3delinsTTA;5A>C]`.
    //
    // **The order between them is nevertheless not observable, and that is
    // measured, not assumed.** Moving this call to *after* the codon exception
    // leaves the whole suite green, the new `c.`-axis test included. The reason
    // is that this rule reads only three things — the hull the pieces span, the
    // sequence they denote, and whether every separation is a single base — and
    // the codon exception preserves all three: it splices the unchanged middle
    // base into the payload explicitly (denoted sequence unchanged), it grows
    // the left piece by exactly what it takes from the right (hull unchanged),
    // and it only ever *closes* a one-base gap, never opens one, so the gaps
    // that remain are a subset of the gaps that were there
    // (`every_separation_is_a_single_base`'s verdict unchanged). A codon-merged
    // partition therefore reconstructs to the same span, and this rule reaches
    // the same `inv` from either side.
    //
    // Two reasons to keep the order written down anyway. It stops being neutral
    // the moment the codon exception can merge across more than one unchanged
    // base, since that would break the third invariant above and let a widened
    // partition fail the gate. And where both fire, the `inv` is the answer
    // regardless of who reaches it — the codon exception's product is a
    // `delins`, which `delins.md:5` defines as a replacement "which is not a
    // substitution or inversion", so over a reverse-complement span it is not a
    // description the spec admits.
    //
    // Pinned by `issue_1040_inv_overrecognition_probes::
    // a_derived_whole_block_inversion_outranks_the_codon_exception`.
    coalesce_whole_block_inversion(&mut pieces, &ref_bytes);

    // Applied *after* the weight bound, deliberately. The bound is a statement
    // about the re-derived partition — it may not describe more change than the
    // input did — and the codon exception (`general.md:35`) is a licensed
    // widening on top of an already-accepted partition, not part of what the
    // bound judges. Running it first would let a legitimate codon merge inflate
    // the weight and trip a refusal.
    apply_coding_codon_exception(&mut pieces, frame.reading_frame, w_lo, &ref_bytes);
    // A veto stood here that refused the whole group whenever any derived piece
    // was an inversion. It was not an input-relative gate — it read only the
    // pieces and the reference, so two spellings of one variant always got the
    // same verdict. It was a capability gap: the derivation could shred a whole
    // reverse complement into a run of substitutions (`ACGTGC -> GCACGT` splits
    // into `[A>G, GT>AC, C>T]`), and declining was better than emitting that.
    //
    // `coalesce_whole_block_inversion` above closes the gap at its source, by
    // putting such a block back into one piece before it is rendered. That piece
    // is then a single span whose replacement is its own reverse complement,
    // which `crate::normalize::rules`' single-span typing renders as `inv` —
    // so there is nothing left for a veto to protect.

    // There was a veto here, and it was the last gate in this pass that read the
    // *input spelling* rather than the sequence: when the derivation produced
    // fewer pieces than there were members, and any piece covered a base the
    // input had left between two of them, the pass refused and returned the
    // input verbatim. It cited `general.md:34` — two variants separated by one
    // or more unchanged nucleotides are described individually.
    //
    // The citation is sound and the mechanism was not. `general.md:34` is about
    // what the *sequence* separates; the veto measured what the *input's
    // spelling* separated, which is a different quantity whenever a member sits
    // in an ambiguous run. So two spellings of one variant were refused
    // differently, which is precisely the non-confluence #1235 is about — a gate
    // installed to enforce a spec rule was breaking the property the whole pass
    // exists to provide.
    //
    // Where the sequence really does separate two runs of change, the partition
    // already keeps them apart, and it does so from the block alone: that is
    // `MIN_PIECE_SEPARATION` plus `separations_are_meaningful` on this — the
    // *live* — path through `partition_block`. The veto added nothing there; it
    // only ever fired where the block said "one member" and the input spelling
    // said "two".
    //
    // Not `MIN_SEPARATION_NO_FRAME`, which an earlier revision of this comment
    // cited. That constant has the same value (`1`) and the same derivation, but
    // its only non-test use is inside the `seqfirst_shadow_enabled()` block, so
    // it holds no live line and could not have been the net this argument leans
    // on. The two are twins by design (see its own doc comment); the live one is
    // the one worth naming here.

    // A derivation collapsing to a single pure insertion used to be refused
    // here, on the grounds that the per-member pipeline owned two capabilities
    // this one lacked: the terminal-insertion clamps (#1205, #1217, and #129's
    // ban on a bare `ins` at the mito origin) and duplication recovery. Both now
    // live in the piece builder — `boundary_delins_anchor` and
    // `duplication_anchor` — so the refusal has nothing left to protect, and it
    // was the single largest source of #1235's non-confluence: every allele
    // whose members denote one insertion was handed back to a pipeline that
    // normalizes members in isolation and therefore cannot merge them.
    //
    // Do not restore it. If a terminal or tandem case regresses, the fix belongs
    // in the builder; refusing here also refuses every merge the derivation gets
    // right.

    // Only a pure insertion resting on a window edge can be clamped, and only
    // then is it worth asking the provider how long the whole sequence is — the
    // one question `ref_bytes` alone cannot answer, since a window that stops
    // short of the end and one that is flush with it look identical from here.
    let ends = if pieces
        .iter()
        .any(|piece| piece.rests_on_a_window_edge(&ref_bytes))
    {
        window_sequence_ends(provider, &accession, start0, ref_bytes.len())
    } else {
        SequenceEnds::INTERIOR
    };

    // The CDS end on the member axis, for `cds_end_delins_anchor`. `None`
    // wherever there is no CDS to change axis at — every genomic kind, and a
    // non-coding transcript — which makes that clamp inert there.
    //
    // `CisKind::Cds` only: `r.` reaches the same junction, but `normalize_cds`'s
    // #387 clamp is gated on `AxisRegion::Cds`, so claiming the `r.` case here
    // would make the derivation *disagree* with the per-member pipeline in the
    // opposite direction — the very fault this fixes. Left for its own evidence.
    let cds_end_axis = if matches!(kind, CisKind::Cds) {
        provider
            .get_transcript(&accession)
            .ok()
            .and_then(|tx| tx.cds_end)
            .and_then(|end| i64::try_from(end).ok())
            .map(|end| end - frame.delta)
    } else {
        None
    };

    let rebuilt = rebuild_members(
        &pieces,
        &variants[0],
        body,
        w_lo,
        &ref_bytes,
        ends,
        cds_end_axis,
    )?;
    if rebuilt == variants {
        return None;
    }
    // Never let canonicalization change what the variant means. This is a
    // *runtime* refusal, not a debug assertion: the whole point of #1234 is that
    // a normalizer silently emitting a different variant is the worst failure
    // this crate can have, and a `debug_assert` is compiled out of exactly the
    // release builds that process real data. Re-deriving the result from the
    // rebuilt members and comparing costs one window-sized apply on a path that
    // has already fetched and applied that window once.
    //
    // If it does not round-trip, fall back to the per-member pipeline's output
    // rather than emitting the re-derivation. A missed canonicalization is a
    // cosmetic defect; a changed sequence is data corruption.
    //
    // The three ways this can decline are kept apart on purpose. A rebuilt form
    // that will not lower back to edits, and one whose edits will not re-apply,
    // are both "cannot tell" — the pass declines and the per-member output
    // stands, which is unremarkable. A round trip that *succeeds* and disagrees
    // is the interesting one: it means the canonicalizer built a form denoting
    // different bases, which is a defect in this module rather than an input it
    // cannot serve. Collapsing all three into one boolean would hide that case
    // among the other two, so it gets a `debug_assert` of its own — loud in
    // tests and development, with the runtime refusal still carrying release.
    let rebuilt_edits = collect_canonical_edits(&rebuilt, kind, body, &template_accession)?;
    let reapplied = apply_edits_to_window(&rebuilt_edits, &ref_bytes, w_lo)?;
    debug_assert_eq!(
        reapplied, result,
        "sequence-first canonicalization changed the resulting sequence"
    );
    if reapplied != result {
        return None;
    }
    Some(rebuilt)
}
/// The `CisKind` a variant belongs to, or `None` for an axis the sequence-first
/// path does not serve (protein has no apply-to-reference; `o.` is circular).
fn cis_kind_of(variant: &HgvsVariant) -> Option<CisKind> {
    match variant {
        HgvsVariant::Genome(_) => Some(CisKind::Genome),
        HgvsVariant::Mt(_) => Some(CisKind::Mt),
        HgvsVariant::Cds(_) => Some(CisKind::Cds),
        HgvsVariant::Tx(_) => Some(CisKind::Tx),
        HgvsVariant::Rna(_) => Some(CisKind::Rna),
        _ => None,
    }
}

/// Lower every member to a [`GEdit`], or refuse the group.
///
/// Refuses anything whose description carries information the resulting
/// sequence does not: conversions, methylation, copy number, identity/`=`,
/// mixed repeats (`MultiRepeat`), and uncertain edits. Also refuses a mixed
/// axis, a second accession, and any member outside the positive body region (a
/// contiguous window would happily produce `c.-1_1`, which is not valid HGVS).
///
/// # Repeats are lowered, not refused (#1296)
///
/// A single-unit tandem repeat *does* denote a sequence, so it becomes a
/// [`GEdit::Repeat`] here and [`lower_repeat_edits`] resolves it against the
/// window into the `ins`/`del` it means. That is a change of what the
/// derivation may **read**, not of what it emits: repeat *output* stays with
/// the per-member pipeline's renderer, which promotes a derived `ins` straight
/// back to `unit[N]`.
///
/// Refusing was a permanent lockout rather than a single declined case. The
/// per-member pipeline *promotes* members into repeat notation, so
/// `[272_273del;274_275insAA]` normalizes to `[272_273del;274A[3]]` — and on
/// the next pass the derivation refused the group it had just produced, so it
/// could never run on that variant again.
///
/// This takes no side in the repeat-vs-del/dup question — `open-issues.md:88`
/// says both spellings are correct, and that is why the old refusal cited it.
/// Reading a repeat does not choose a spelling for it.
///
/// `MultiRepeat` stays refused, and the reason is **not** that a mixed repeat
/// has no single unit to lower. That describes the shape of [`lowered_repeat`],
/// not an impossibility: a mixed repeat could be expanded unit by unit and the
/// concatenation diffed against its tract.
///
/// The reason is that the lockout above **cannot arise for one**. Nothing in
/// ferro ever *produces* a `MultiRepeat`:
/// [`normalize_na_edit`](crate::normalize) gives it an explicit pass-through arm
/// — deliberately not 3'-shifted, routed only for tract validation — and the
/// only other construction site is a case-folding rebuild in
/// [`rules`](crate::normalize::rules). A `MultiRepeat` appears in ferro's output
/// only if the submitter wrote one, so the "refused the group it had just
/// produced" failure has no mixed-repeat analogue. Refusing is not free, though:
/// a **submitter-written** `MultiRepeat` beside a sibling still cannot be
/// merged. That single case is the whole of what it costs.
///
/// Two further grounds, each covering only part of the space:
///
/// * a `MultiRepeat` carrying any non-`Exact` count has no determined resulting
///   sequence at all — the same rule this function already applies to uncertain
///   edits;
/// * **for the all-`Exact` case only**,
///   [`validate_multirepeat_tract`](crate::normalize::validate) requires the
///   declared units to expand to exactly the reference bases over the stated
///   span, so under ferro's semantic such a member denotes no change — an
///   identity, which this pass refuses on principle because the resulting
///   sequence cannot express it. This argument does not generalise past that
///   branch, and is not what settles the question.
fn collect_canonical_edits(
    variants: &[HgvsVariant],
    kind: CisKind,
    body: Region,
    template_accession: &Accession,
) -> Option<Vec<GEdit>> {
    let mut edits = Vec::with_capacity(variants.len());
    for v in variants {
        let (accession, region, s, e, edit) = cis_axis_parts(v, kind)?;
        if accession != template_accession || region != body {
            return None;
        }
        match edit {
            NaEdit::Insertion { sequence } => {
                let bases = sequence.bases()?;
                if e != s + 1 {
                    return None;
                }
                edits.push(GEdit::Ins {
                    gap: s,
                    alt: bases.to_vec(),
                });
            }
            NaEdit::Deletion { .. } => edits.push(GEdit::Del { s, e }),
            NaEdit::Substitution { alternative, .. } => {
                if s != e {
                    return None;
                }
                edits.push(GEdit::Sub {
                    pos: s,
                    alt: *alternative,
                });
            }
            NaEdit::Delins { sequence, .. } => {
                let bases = sequence.bases()?;
                edits.push(GEdit::Delins {
                    s,
                    e,
                    alt: bases.to_vec(),
                });
            }
            NaEdit::Inversion { .. } => {
                if e <= s {
                    return None; // a 1-nt "inversion" is not valid (inversion.md:16)
                }
                edits.push(GEdit::Inv { s, e });
            }
            NaEdit::Duplication {
                uncertain_extent: None,
                ..
            } => edits.push(GEdit::Dup { s, e }),
            // Deferred until the window exists — see [`GEdit::Repeat`]. The
            // guards here are exactly the per-member pipeline's own admission
            // test for a shufflable repeat (`normalize_na_edit`'s `Repeat`
            // arm): a stated unit, one exact count, and no genotype or VEP
            // trailing baggage. Anything else has no single resulting sequence
            // to derive from.
            // `e >= s` joins the admission test rather than being left to
            // `lowered_repeat`. That function skips its stated-end branch when
            // `e < s`, so a reversed anchor does not crash — it is silently
            // re-read as the single-position `g.10AC[3]`, which is a different
            // variant from the `g.10_5AC[3]` that was written. SVD-WG006's
            // reversed-range provision covers circular deletions and
            // duplications, not repeats, so there is nothing here to honour and
            // declining is the answer.
            NaEdit::Repeat {
                sequence: Some(sequence),
                count: RepeatCount::Exact(count),
                additional_counts,
                trailing: None,
            } if additional_counts.is_empty() && e >= s => {
                let unit: Vec<u8> = sequence
                    .bases()
                    .iter()
                    .map(|b| canonical_base_byte(b.to_u8()))
                    .collect();
                if unit.is_empty() {
                    return None;
                }
                edits.push(GEdit::Repeat {
                    s,
                    e,
                    unit,
                    count: *count,
                });
            }
            _ => return None,
        }
    }
    Some(edits)
}

/// Resolve every [`GEdit::Repeat`] into the `ins`/`del` it denotes, or refuse
/// the group.
///
/// A repeat says "the run of `unit` covering my anchor holds `count` copies".
/// The run itself is in the reference, so this is the first point at which the
/// member denotes anything: `ref_bytes`/`w_lo` are the window
/// [`canonicalize_from_sequence`] has just fetched.
///
/// Tract-finding is [`count_tandem_repeats`]', the same one
/// [`normalize_repeat`](crate::normalize::rules::normalize_repeat) uses, over
/// the same [`smallest_repeat_unit`] canonicalization and the same flank
/// absorption for an under-specified explicit range — **minus that function's
/// off-phase rotation search**, which is the one place the two part company.
///
/// For a single-anchor member whose unit is two or more bases, `normalize_repeat`
/// (#852) tries every rotation of the unit at every anchor offset and keeps the
/// longest tract that spans the anchor; this takes the literal lookup only. So
/// the two agree on which tract a member names whenever it is written in phase
/// with its tract, and an off-phase single-anchor member is **refused** here
/// rather than lowered onto some other tract — it fails the spanning check or
/// the room check below. That is a capability gap, not a disagreement: the group
/// falls back to the per-member pipeline, which does run the rotation search.
///
/// What is deliberately *not* shared is that function's choice between `dup`,
/// `del`, `ins` and `unit[N]`: that is a rendering decision, and the derivation
/// only needs the bases.
///
/// Where inside the tract the difference is placed does not matter **as long as
/// no sibling sits inside it**. The tract is a whole number of copies of one
/// unit, so adding or removing copies at either end yields the same sequence —
/// which is why this can otherwise ignore the shuffle direction that
/// `normalize_repeat` must honour. The 3' end is chosen so the footprint stays
/// as small and as far from the group's other members as possible.
///
/// That premise fails when another member's junction falls strictly inside the
/// tract, and it fails silently. `g.[263_265A[7];264_265insC]` puts the
/// insertion's `C` between the tract's own bases, so the four added `A`s are
/// *before* it; lowering the repeat to an insertion at the tract's 3' junction
/// puts them after, and the applied window comes back as `TAACAAAAAT` instead of
/// `TAAAAAACAT`. Nothing downstream can notice, because by then the repeat is an
/// ordinary `Ins` and the tract it came from is gone — the derivation simply
/// canonicalizes the wrong sequence, and the caller gets a well-formed
/// description of bases the input never denoted.
///
/// Such a pair overlaps and has no single resulting sequence to derive from, so
/// this refuses rather than picks an order — the same answer
/// `apply_edits_to_window` already gives an insertion interior to a `del`,
/// `delins` or `inv` span (#486), applied to the one span-bearing member that
/// reaches here still carrying its tract.
///
/// Refuses (leaving the group to the per-member pipeline) when:
///
/// * the anchor is outside the window, or the tract runs into either window
///   edge — a truncated tract would under-count the copies and denote the wrong
///   sequence. [`count_tandem_repeats`] stops scanning when it has less than a
///   whole unit left in that direction, so "a whole unit of room on each side"
///   is exactly the condition under which it stopped for a real reason;
/// * the tract does not span the anchor, which `count_tandem_repeats` can
///   return (it back-scans first, so a run lying entirely 5' of the anchor is
///   reachable);
/// * another member's junction falls strictly inside a repeat's stated span,
///   per the paragraph above — the pair overlaps and denotes no single
///   sequence, so there is nothing to derive from;
/// * `count` equals the reference count, i.e. the member is an identity. The
///   resulting sequence cannot express one, so a group carrying one would come
///   back with that member silently dropped — the same reason `NaEdit::Identity`
///   is refused above.
fn lower_repeat_edits(edits: &mut [GEdit], ref_bytes: &[u8], w_lo: i64) -> Option<()> {
    // The junction each member adds sequence at, for the interior test below.
    // An insertion's is its gap; a duplication's is the 3' end of what it
    // copies. Every other kind claims bases rather than a junction and is
    // already covered by the applier's own interior test.
    let junction_of = |edit: &GEdit| match edit {
        GEdit::Ins { gap, .. } => Some(*gap),
        GEdit::Dup { e, .. } => Some(*e),
        _ => None,
    };
    for i in 0..edits.len() {
        let GEdit::Repeat { s, e, .. } = edits[i] else {
            continue;
        };
        // Strictly inside, matching `demote_repeats_spanning_siblings`'
        // `spans_sibling_junction`: a junction at the tract's 3' end is flush
        // against it rather than within it, and that adjacency is legitimate.
        if (0..edits.len())
            .filter(|&j| j != i)
            .filter_map(|j| junction_of(&edits[j]))
            .any(|junction| junction >= s && junction < e)
        {
            return None;
        }
    }
    for edit in edits.iter_mut() {
        let GEdit::Repeat { s, e, unit, count } = edit else {
            continue;
        };
        *edit = lowered_repeat(*s, *e, unit, *count, ref_bytes, w_lo)?;
    }
    Some(())
}

/// The `ins`/`del` one repeat member denotes, or `None`. See
/// [`lower_repeat_edits`] for the refusal conditions.
fn lowered_repeat(
    s: i64,
    e: i64,
    unit: &[u8],
    count: u64,
    ref_bytes: &[u8],
    w_lo: i64,
) -> Option<GEdit> {
    // Canonicalize to the smallest period, exactly as `normalize_repeat` does,
    // so a member spelling a non-minimal unit (`ATAT[1]` over an `AT[4]` tract)
    // finds the same tract the per-member pipeline would.
    let canonical = crate::normalize::rules::smallest_repeat_unit(unit);
    let copies_per_stated_unit = (unit.len() / canonical.len()) as u64;
    let mut count = count.checked_mul(copies_per_stated_unit)?;
    let unit_len = canonical.len();

    let anchor = usize::try_from(s.checked_sub(w_lo)?).ok()?;
    if anchor >= ref_bytes.len() {
        return None;
    }
    let (ref_count, tract_start, tract_end) =
        crate::normalize::rules::count_tandem_repeats(ref_bytes, anchor, canonical)?;
    if tract_start > anchor || anchor >= tract_end {
        return None;
    }
    if tract_start < unit_len || tract_end + unit_len > ref_bytes.len() {
        return None;
    }

    // Flank absorption, mirroring `normalize_repeat`: an explicit range states
    // which reference copies the count applies to, and copies inside the tract
    // but outside that range are untouched by the edit, so they survive into the
    // result and belong in the count. A single-position anchor (`e == s`) names
    // only the tract's start and means "the whole run becomes N", so it absorbs
    // nothing.
    if e > s {
        let stated_end = usize::try_from(e.checked_sub(w_lo)?).ok()?;
        let three_prime = tract_end.saturating_sub(stated_end + 1);
        let five_prime = anchor.saturating_sub(tract_start);
        // `checked_add`, not `+=`: `count` is the submitter's copy number times
        // `copies_per_stated_unit`, so a description like `g.5_6A[18446744073709551615]`
        // reaches here already near the top of the range and absorbing a flank
        // would wrap. Every other arithmetic step in this function is checked;
        // this was the one that was not.
        count = count.checked_add(((three_prime + five_prime) / unit_len) as u64)?;
    }

    match count.cmp(&ref_count) {
        // Identity: the member denotes no change at all.
        std::cmp::Ordering::Equal => None,
        std::cmp::Ordering::Greater => {
            let added = usize::try_from(count - ref_count).ok()?;
            // A repeat count is the one place a *short* description names an
            // arbitrarily *long* payload: `g.274A[4000000000]` is thirteen
            // characters and would size a four-gigabyte `Vec`, which
            // `with_capacity` answers by aborting the process — in release as
            // well as debug, and reachable from `Normalizer::normalize`. Every
            // other `alt` here is bounded by the length of the string the
            // submitter actually wrote. Bound it by the window instead: an
            // insertion wider than the whole window has no partition to
            // contribute to, so refusing costs nothing.
            let payload = added.checked_mul(unit_len)?;
            if payload > MAX_CANONICAL_WINDOW as usize {
                return None;
            }
            let mut alt = Vec::with_capacity(payload);
            for _ in 0..added {
                for byte in canonical {
                    alt.push(Base::from_char(*byte as char)?);
                }
            }
            Some(GEdit::Ins {
                gap: w_lo + tract_end as i64 - 1,
                alt,
            })
        }
        std::cmp::Ordering::Less => {
            let removed = usize::try_from(ref_count - count)
                .ok()?
                .checked_mul(unit_len)?;
            let del_start = tract_end.checked_sub(removed)?;
            Some(GEdit::Del {
                s: w_lo + del_start as i64,
                e: w_lo + tract_end as i64 - 1,
            })
        }
    }
}

/// Whether every member's *stated* reference base agrees with the window.
///
/// `NaEdit` optionally carries the reference bases the submitter asserted
/// (`c.5A>T`, `c.5_7delACG`, `c.263A=`). Where present they must match, otherwise
/// the variant is a reference mismatch and belongs to the per-member pipeline's
/// strict-reject / warn path rather than to canonicalization.
fn stated_reference_bases_match(
    variants: &[HgvsVariant],
    kind: CisKind,
    ref_bytes: &[u8],
    w_lo: i64,
) -> bool {
    for v in variants {
        let Some((_, _, s, _, edit)) = cis_axis_parts(v, kind) else {
            return false;
        };
        // Every `NaEdit` channel that can carry submitter-asserted reference
        // bases starting at `s` is listed here. `canonicalize_edit` strips most
        // of them before this pass runs, so several arms are inert today — but a
        // validator that silently ignores half the channels it claims to cover is
        // one pipeline reordering away from letting a reference mismatch be
        // canonicalized out from under strict mode.
        let stated = match edit {
            NaEdit::Substitution { reference, .. } => Some(vec![*reference]),
            // `c.263A=` asserts that position 263 is `A`, so an identity that
            // spells its sequence is a stated-reference channel like any other
            // (#1352). It was missing while the doc comment above claimed the
            // list was complete, and the safety therefore rested on a refusal in
            // a *different* function — `collect_canonical_edits` declines
            // `Identity` at its catch-all — which is exactly the coupling that
            // comment warns against. Inert until something upstream stops
            // refusing identity members, which is the point: the validator
            // should be correct before the refusal is relaxed, not as part of
            // relaxing it.
            //
            // `..` rather than `whole_entity: false`: a whole-entity identity
            // (`c.=`) carries no sequence in anything the parser produces, but if
            // one ever arrives with bases, comparing them is the fail-safe
            // answer — a mismatch refuses canonicalization, whereas skipping the
            // arm is the very hole being closed.
            NaEdit::Identity {
                sequence: Some(seq),
                ..
            }
            | NaEdit::Deletion {
                sequence: Some(seq),
                ..
            }
            | NaEdit::Duplication {
                sequence: Some(seq),
                ..
            }
            | NaEdit::Inversion {
                sequence: Some(seq),
                ..
            }
            | NaEdit::Delins {
                deleted: Some(seq), ..
            }
            | NaEdit::Delins {
                substitution_reference: Some(seq),
                ..
            } => Some(seq.bases().to_vec()),
            _ => None,
        };
        let Some(stated) = stated else { continue };
        for (offset, base) in stated.iter().enumerate() {
            let index = s - w_lo + offset as i64;
            if index < 0 || index as usize >= ref_bytes.len() {
                return false;
            }
            if !ref_bytes[index as usize].eq_ignore_ascii_case(&base.to_u8()) {
                return false;
            }
        }
    }
    true
}

/// Every reference mismatch the **authored** members of a cis allele carry.
///
/// # Why this exists at all
///
/// The per-member pipeline validates a member's stated reference bases in
/// `normalize_na_edit`, and that is where `RefSeqMismatch` (W5001) — and so
/// strict mode's rejection and lenient mode's warning — comes from. But
/// `normalize_allele` runs [`collapse_overlapping_cis_edits`] and
/// [`merge_consecutive_edits`] over the raw members **before** any member
/// reaches that validator, and a merge keeps only each member's *alt* bases:
/// [`Anchor`] has no channel for the reference the submitter asserted. So a
/// member consumed by a merge has its assertions dropped, unexamined, and the
/// merged anchor is then rebuilt from the **real** bases — turning a
/// description whose reference claims are false into a well-formed one, with
/// `status=ok`, in the mode whose purpose is to reject exactly that (#1543).
///
/// Measured on `main` at `5616cdb9`, reference bases `c.1 = A`, `c.3 = G`:
///
/// ```text
/// NM_000143.3:c.1G>T          -> rejected  (no sibling to merge with)
/// NM_000143.3:c.[1G>T;500G>C] -> rejected  (sibling too far to merge)
/// NM_000143.3:c.[1G>T;3T>A]   -> c.1_3delinsTTA, status=ok   <- the defect
/// ```
///
/// The discriminator was merge distance and nothing else.
///
/// # Why here, and not at the merge sites
///
/// This is the sixth defect in one family — #1052 (substitutions), #1068 (the
/// `m.` axis), #1092 (the parser discarded the bases), #1097 (multi-base range
/// substitutions), #1352 ([`NaEdit::Identity`] omitted) — and each of the first
/// five was fixed one shape at a time. A per-shape guard inside each merge pass
/// would be the sixth such patch and would leave the same ordering hazard for
/// the seventh: any future pass that consumes a member before the validator
/// runs re-opens the hole.
///
/// So the question is asked **on the authored input**, once, before any pass
/// has had the chance to strip anything. That is ordering-immune by
/// construction: no rearrangement of the pipeline below can make this check
/// vacuous, which is precisely what happened to [`stated_reference_bases_match`]
/// — whose own doc comment predicted it.
///
/// Nothing about the merge itself changes, deliberately: the warnings this
/// returns are *added* to the member-local ones, so no correctly-spelled input
/// moves and the only behavioural change is that a false reference assertion is
/// now reported wherever it was authored.
///
/// # Contract
///
/// Cis only — merging is cis-only, so that is the whole exposure, and a
/// trans/unknown-phase allele already has every member normalized in isolation.
/// The verdict itself comes from [`crate::normalize::validate::validate_reference`],
/// the same function the per-member pipeline uses, so the channel list cannot
/// drift between the two sites the way #1352's could.
///
/// Members this cannot answer for — an offset (intronic) position, a region
/// with no served bases, an unreadable or short provider window — are skipped
/// rather than guessed at, exactly as every other refusal in this module does.
pub(crate) fn authored_member_reference_mismatches<P: ReferenceProvider>(
    variants: &[HgvsVariant],
    phase: AllelePhase,
    provider: &P,
) -> Vec<crate::normalize::NormalizationWarning> {
    if phase != AllelePhase::Cis {
        return Vec::new();
    }
    variants
        .iter()
        .filter_map(|v| authored_reference_mismatch(v, provider))
        .collect()
}

/// One member's authored reference mismatch, if it has one.
///
/// The returned warning is built to be **indistinguishable** from the one the
/// per-member pipeline raises for the same member: same `position` (the
/// sequence-axis span, via [`region_sequence_delta`]), same `corrected` flag,
/// same `details` text. That is what lets the caller drop it when the member
/// also reached the per-member validator, rather than reporting one finding
/// twice.
fn authored_reference_mismatch<P: ReferenceProvider>(
    v: &HgvsVariant,
    provider: &P,
) -> Option<crate::normalize::NormalizationWarning> {
    let kind = cis_kind_of(v)?;
    let (accession, region, s, e, edit) = cis_axis_parts(v, kind)?;
    // `r.` spells its bases in the RNA alphabet while the transcript is served
    // as DNA, so a truthful `r.2u>a` over a `T` would otherwise read as a
    // mismatch. `normalize_rna` makes the same rewrite before it validates
    // (#736); making it here too is what keeps the two verdicts identical.
    let rewritten = matches!(kind, CisKind::Rna)
        .then(|| crate::normalize::rules::rna_uracil_to_thymine(edit))
        .flatten();
    let edit = rewritten.as_ref().unwrap_or(edit);

    let provider_key = accession.transcript_accession();
    let delta = region_sequence_delta(region, &provider_key, provider)?;
    // Checked, like every other coordinate conversion in this file: `s`/`e`
    // come off a parsed description and `delta` off the record's CDS bounds, so
    // an adversarial span could otherwise overflow and panic in a debug build
    // where declining is the answer.
    let start = s.checked_add(delta)?;
    let end = e.checked_add(delta)?;
    if start < 1 || end < start {
        return None;
    }
    let span = usize::try_from(end.checked_sub(start)?.checked_add(1)?).ok()?;
    let bases = provider
        .get_sequence(
            &provider_key,
            u64::try_from(start - 1).ok()?,
            u64::try_from(end).ok()?,
        )
        .ok()?;
    // A short read means the member runs off the end of the sequence. That is
    // W4004 `PositionPastEnd`'s finding, not this one, so decline rather than
    // compare against a truncated window.
    if bases.len() != span {
        return None;
    }

    // `1`/`span` rather than the member's own coordinates: `validate_reference`
    // indexes `ref_seq` from `start`, and the window fetched above starts *at*
    // the member. Every arm it dispatches to reads only `[start - 1, end)` or
    // the span width, so re-basing the pair is equivalent to handing it the
    // whole sequence and costs one small read instead of the entire transcript.
    let result =
        crate::normalize::validate::validate_reference(edit, bases.as_bytes(), 1, span as u64);
    if result.valid {
        return None;
    }
    // Mirrors the `corrected` reasoning in `normalize_na_edit`: the canonical
    // `Display` keeps a substitution's stated base and a repeat's unit, so
    // claiming those were corrected would be dishonest.
    let corrected = !matches!(
        edit,
        NaEdit::Repeat { .. } | NaEdit::MultiRepeat { .. } | NaEdit::Substitution { .. }
    );
    Some(crate::normalize::NormalizationWarning::RefSeqMismatch {
        stated_ref: result.stated_ref.unwrap_or_default(),
        actual_ref: result.actual_ref.unwrap_or_default(),
        position: format!("{start}-{end}"),
        corrected,
        details: result.warning,
    })
}

/// Inclusive 1-based span covering every edit's reference footprint.
fn edit_span_union(edits: &[GEdit]) -> Option<(i64, i64)> {
    let mut lo = i64::MAX;
    let mut hi = i64::MIN;
    for edit in edits {
        let (s, e) = match edit {
            GEdit::Ins { gap, .. } => (*gap, *gap + 1),
            GEdit::Del { s, e }
            | GEdit::Delins { s, e, .. }
            | GEdit::Inv { s, e }
            | GEdit::Dup { s, e }
            // The *stated* anchor, which is all an unlowered repeat has. Its
            // tract can reach past it in either direction, so this is the one
            // shape whose real footprint is wider than its span — that is what
            // `CANONICAL_PAD` gives it room for, and `lower_repeat_edits`
            // refuses any tract that reaches a window edge anyway.
            | GEdit::Repeat { s, e, .. } => (*s, *e),
            GEdit::Sub { pos, .. } => (*pos, *pos),
        };
        lo = lo.min(s);
        hi = hi.max(e);
    }
    (lo <= hi).then_some((lo, hi))
}

/// How a cis group's own axis sits on the sequence the provider serves.
///
/// The two facts travel together because both are answered by the same
/// `cds_start` lookup, and a second provider round-trip to re-answer the other
/// half would be pure waste on a path that already fetches a window.
/// Answers per *axis*, which is enough for its one caller because
/// [`canonicalize_from_sequence`] restricts itself to the positive body region.
/// The repair passes are not so restricted and must use
/// [`region_sequence_delta`], which answers per *region* — see its doc comment
/// for why a per-axis offset cannot be right in the UTRs.
struct AxisFrame {
    /// Offset between the member axis and the fetched sequence: 0 for
    /// `g.`/`m.`/`n.`, `cds_start - 1` for the CDS-relative `c.`/`r.` axes.
    ///
    /// Correct only inside the positive body region; see the note above.
    delta: i64,
    /// Whether positions on this axis carry a **reading frame** — i.e. whether
    /// [`same_codon`] means anything here. True only when the axis is
    /// CDS-relative *and* the transcript actually has a CDS.
    ///
    /// Not the same question as `delta`: a non-coding `r.` falls back to
    /// `delta == 0` and is still an `r.` description, so keying the codon
    /// exception off the axis alone stamps a reading frame onto a transcript
    /// that has none (#1241).
    reading_frame: bool,
}

/// Whether `[lo, hi]` — a member-span union in axis coordinates — crosses an
/// exon/exon junction of its own transcript (#1450).
///
/// The genomic axes have no exons and so never cross one. A transcript the
/// provider cannot serve is reported as not crossing: the derivation is then no
/// worse off than it was before this check existed, which is the same fallback
/// `axis_frame` takes for a missing transcript.
fn crosses_exon_junction<P: ReferenceProvider>(
    kind: CisKind,
    accession: &str,
    provider: &P,
    frame: &AxisFrame,
    lo: i64,
    hi: i64,
) -> bool {
    if matches!(kind, CisKind::Genome | CisKind::Mt) {
        return false;
    }
    let Ok(transcript) = provider.get_transcript(accession) else {
        return false;
    };
    // `axis + delta` is the 1-based transcript coordinate — the same conversion
    // the window fetch performs immediately below the call site.
    //
    // Checked, and `try_from` rather than `as`, so neither an extreme authored
    // coordinate nor an `exon.end` past `i64::MAX` can wrap this predicate into
    // a wrong answer — a wrapped `junction` compares as negative and silently
    // reports "no junction crossed", which is the direction that reinstates
    // #1450. Both failures take the same conservative exit as an unservable
    // transcript above.
    //
    // Note this hardens the predicate, not the whole path: an extreme
    // coordinate still panics upstream at `w_hi = c_hi + CANONICAL_PAD`, which
    // runs before this is reached and is untouched by this change. Tracked
    // separately.
    let (Some(tx_lo), Some(tx_hi)) = (lo.checked_add(frame.delta), hi.checked_add(frame.delta))
    else {
        return false;
    };
    transcript.exons.iter().any(|exon| {
        i64::try_from(exon.end).is_ok_and(|junction| tx_lo <= junction && junction < tx_hi)
    })
}

/// Resolve the group's [`AxisFrame`], or refuse the group.
fn axis_frame<P: ReferenceProvider>(
    kind: CisKind,
    accession: &Accession,
    provider: &P,
) -> Option<AxisFrame> {
    match kind {
        CisKind::Genome | CisKind::Mt | CisKind::Tx => Some(AxisFrame {
            delta: 0,
            reading_frame: false,
        }),
        CisKind::Cds => {
            let tx = provider
                .get_transcript(&accession.transcript_accession())
                .ok()?;
            Some(AxisFrame {
                delta: i64::try_from(tx.cds_start?).ok()? - 1,
                reading_frame: true,
            })
        }
        // `r.` is CDS-relative on a coding transcript but transcript-relative on
        // a non-coding one, which has no `cds_start` at all. Refusing there
        // would leave every non-coding `r.` description outside the
        // canonicalization; fall back to transcript-1, matching what
        // `fetch_ref_for_canonical_split` does for the same case.
        CisKind::Rna => {
            let cds_start = provider
                .get_transcript(&accession.transcript_accession())
                .ok()
                .and_then(|tx| tx.cds_start);
            match cds_start {
                Some(cds_start) => Some(AxisFrame {
                    delta: i64::try_from(cds_start).ok()? - 1,
                    reading_frame: true,
                }),
                // Transcript-relative: either a non-coding transcript, or a
                // provider that cannot resolve it. Both mean "no CDS offset",
                // not "refuse" — `Tx` reaches the same conclusion without ever
                // consulting the provider. Neither has a reading frame: the
                // non-coding transcript has no CDS, and an unresolvable one
                // gives no grounds to claim it has.
                None => Some(AxisFrame {
                    delta: 0,
                    reading_frame: false,
                }),
            }
        }
    }
}

/// Offset from a member's own **region** axis to the 1-based sequence the
/// provider serves, or `None` for a region that is not part of that sequence.
///
/// `position_on_the_sequence = axis_position + delta`.
///
/// # Why this is keyed on `Region` and not on `CisKind`
///
/// [`AxisFrame`] answers the same question one level coarser — per *axis* —
/// and that is all its caller needs, because [`canonicalize_from_sequence`]
/// restricts itself to the positive body region. The repair passes have no
/// such restriction: they operate wherever a collision settles, including
/// both UTRs.
///
/// Per-axis is the wrong granularity there, and #1284 records the resulting
/// dead end. A single `cds_start - 1` for the whole `c.` axis is right at
/// `c.1` and off by one below it, because the CDS axis has **no zero**:
/// `c.-1` is `cds_start - 1`, not `cds_start - 2`. Patching that with a sign
/// test (`axis + delta + 1` when negative) gets the endpoints right and still
/// leaves the 3'UTR — a third offset entirely, `cds_end` — wrong.
///
/// Keyed on the region the discontinuity disappears rather than being patched
/// across, because each region *is* individually affine onto the transcript:
///
/// | region | axis → transcript | delta |
/// |---|---|---|
/// | `Genome`, `Tx` | the axis is the sequence | `0` |
/// | `Cds`, `Rna` (`c.N`) | `cds_start + N - 1` | `cds_start - 1` |
/// | `FivePrimeUtr` (`c.-N`) | `cds_start - N` | `cds_start` |
/// | `ThreePrimeUtr` (`c.*N`) | `cds_end + N` | `cds_end` |
///
/// and every pass that uses it already groups members within one region, so a
/// span never straddles two.
///
/// `TxUpstream` / `TxDownstream` (`n.-N` / `n.*N`) are refused: those positions
/// lie outside the transcript the provider serves, so there are no bases to
/// read at all. That is a genuine refusal, not a missing conversion.
fn region_sequence_delta<P: ReferenceProvider>(
    region: Region,
    provider_key: &str,
    provider: &P,
) -> Option<i64> {
    // `ordered_cds_bounds` for every CDS-relative region, not just
    // `ThreePrimeUtr`. Reading `cds_start` on its own accepts a record whose
    // bounds are inverted, and on such a record the `c.` axis is incoherent —
    // its 5'UTR and 3'UTR halves overlap, so one transcript position carries two
    // names. Refusing in one region and not the others is the worst of both: the
    // 3'UTR declined while the 5'UTR and the CDS body silently picked a
    // spelling, which is precisely what `ordered_cds_bounds`' own doc comment
    // says this conversion must not do.
    //
    // Written as "refuse if inverted" rather than as `ordered_cds_bounds`
    // directly, because that would also require `cds_end` to be *present*: a
    // record carrying `cds_start` but no `cds_end` has an unknown 3'UTR yet a
    // perfectly well-defined `c.N` body, and refusing it would be a new
    // restriction rather than the guard-parity fix.
    // Read lazily: `Genome`/`Tx` need no transcript at all, and asking for one
    // under a genomic accession is a guaranteed miss.
    let bounds = || -> Option<(Option<i64>, Option<i64>)> {
        let transcript = provider.get_transcript(provider_key).ok()?;
        Some((
            transcript.cds_start.and_then(|v| i64::try_from(v).ok()),
            transcript.cds_end.and_then(|v| i64::try_from(v).ok()),
        ))
    };
    // Bounds that exist and are the wrong way round. Distinguished from "no CDS
    // on this record", because the two want opposite answers: absent bounds are
    // a non-coding transcript, which `Rna` legitimately falls back for; inverted
    // bounds are a malformed one, which every region must refuse.
    let inverted = || matches!(bounds(), Some((Some(start), Some(end))) if end < start);
    let cds_start = || match bounds()? {
        (Some(start), Some(end)) if end < start => None,
        (start, _) => start,
    };
    match region {
        Region::Genome | Region::Tx => Some(0),
        Region::Cds => Some(cds_start()? - 1),
        // Falls back to transcript-relative when there is no `cds_start` —
        // either a non-coding transcript, whose `r.` numbers it directly, or a
        // provider that cannot resolve it. That is the same fallback
        // `axis_frame` makes and for the same reason: refusing would put every
        // non-coding `r.` description outside these repairs. `Cds` above
        // refuses instead, also matching `axis_frame`, because a `c.`
        // description asserts a CDS this record does not have.
        //
        // The fallback must not swallow the *inverted* case, though: `map_or(0)`
        // on its own turns a refusal into a silent delta of 0, which is the
        // guess this conversion exists not to make. Absent bounds fall back;
        // inverted bounds refuse, like every other region here.
        Region::Rna if inverted() => None,
        Region::Rna => Some(cds_start().map_or(0, |start| start - 1)),
        Region::FivePrimeUtr => cds_start(),
        Region::ThreePrimeUtr => Some(ordered_cds_bounds(provider_key, provider)?.1),
        Region::TxUpstream | Region::TxDownstream => None,
    }
}

/// [`region_sequence_delta`] widened to the two `n.` regions that lie *outside*
/// the served sequence, for **comparison only**.
///
/// `region_sequence_delta` refuses `TxUpstream` (`n.-N`) and `TxDownstream`
/// (`n.*N`) because there are no bases to read there, and that is the right
/// answer for every caller that reads bases. It is the wrong answer for
/// [`member_span`], whose callers only ever *compare* members: refusing there
/// does not make a pass conservative, it makes the member invisible, which is
/// the whole of #1482.
///
/// Measured, before this existed: `n.[-5=;-6_-4del]` normalized to
/// `n.-6_-4del` on `main` and to `n.[-6_-4del;-5=]` with `member_span`
/// converting through `region_sequence_delta` alone —
/// `drop_identity_members_covered_by_siblings` stopped seeing that the deletion
/// covers the identity. `n.[-5=;-5del]` was worse: it kept both members on one
/// position, which is the overlap that pass exists to remove.
///
/// The offsets are the ones the axes actually mean, so a cross-region
/// comparison is meaningful rather than merely consistent: `n.-1` is the base
/// immediately 5' of `n.1`, so `n.-N` is transcript position `1 - N`; `n.*N` is
/// `N` bases past the end, so it is `length + N`.
///
/// **The positions it yields for those two regions are virtual** — outside
/// `[1, length]` by construction. Nothing may read bases through them, and
/// nothing does: every base read in this module goes through
/// `provider.get_sequence`, whose `u64::try_from` refuses a non-positive
/// coordinate, and [`respell_at_gap`] refuses such a member outright rather
/// than relying on that.
fn region_span_delta<P: ReferenceProvider>(
    region: Region,
    provider_key: &str,
    provider: &P,
) -> Option<i64> {
    match region {
        Region::TxUpstream => Some(1),
        Region::TxDownstream => {
            i64::try_from(provider.get_sequence_length(provider_key).ok()?).ok()
        }
        other => region_sequence_delta(other, provider_key, provider),
    }
}

/// Both CDS bounds, **refused unless they are ordered**.
///
/// `cds_end < cds_start` is a malformed record whose `c.` axis is incoherent:
/// the 5'UTR and 3'UTR halves overlap, so one transcript position has two
/// names. On such a record `c.*1` and `c.-10` can denote the same base, and a
/// conversion that does not refuse simply picks one of the two spellings —
/// silently, since the result re-parses and so is invisible to the
/// `FERRO_ASSERT_REPARSE` oracle.
///
/// Refusing is not a fix for that ambiguity, which predates this conversion and
/// is reachable without it; it keeps *this* conversion from being one more
/// place that resolves it arbitrarily.
/// `boundary::get_cds_boundaries_with_axis_info` guards the same condition the
/// same way (`(Some(s), Some(e)) if e >= s`); this applies that precedent here.
fn ordered_cds_bounds<P: ReferenceProvider>(
    provider_key: &str,
    provider: &P,
) -> Option<(i64, i64)> {
    let transcript = provider.get_transcript(provider_key).ok()?;
    let start = i64::try_from(transcript.cds_start?).ok()?;
    let end = i64::try_from(transcript.cds_end?).ok()?;
    (end >= start).then_some((start, end))
}

/// The inverse of [`region_sequence_delta`] for the CDS-relative axes: which
/// `(region, axis position)` names transcript position `tx`.
///
/// Needed because a repair's *junction* can land on a region boundary even
/// though its span does not. A duplication resting at `c.-1` copies its bases
/// to the junction between `c.-1` and `c.1` — one interbase point, two
/// regions — and naming it requires converting through the transcript rather
/// than adding 1 to an axis value that has no zero to add it to.
///
/// `body` is the region name this axis gives the CDS proper: `Region::Cds` on
/// `c.`, `Region::Rna` on `r.`. Passed in rather than fixed, because the two
/// axes share every other region name and `member_endpoints` reads the body
/// back under the axis's own label.
fn cds_axis_position(tx: i64, cds_start: i64, cds_end: i64, body: Region) -> Option<(Region, i64)> {
    if tx < 1 {
        return None;
    }
    if tx < cds_start {
        // `c.-N`, counting back from the base before the CDS start.
        Some((Region::FivePrimeUtr, tx - cds_start))
    } else if tx <= cds_end {
        Some((body, tx - cds_start + 1))
    } else {
        Some((Region::ThreePrimeUtr, tx - cds_end))
    }
}

/// Fetch the window `[w_lo, w_hi]` on the member axis, **upper-cased**.
///
/// The padding routinely runs past the end of the sequence, and providers
/// reject an out-of-range request rather than truncating, so clamp to the known
/// length and retry. The window only needs to cover the edits plus whatever
/// shift room exists; a shorter window simply means less room. An empty read
/// means the window is unusable.
///
/// # Why upper-case here
///
/// Reference FASTAs are routinely soft-masked, so a provider may hand back
/// lower-case bases for a repeat-rich region. Case must not reach the rest of
/// this pass, because the pass mixes provider bytes with bytes it writes itself:
/// `apply_edits_to_window` emits `Base::to_u8()`, which is always upper-case,
/// while copying unchanged and duplicated bases through verbatim. The two then
/// meet in comparisons that are exact rather than case-folded —
/// `best_alignment` scores a column as matched only on `reference[i] ==
/// result[j]`, and `pieces_from_columns` calls a column changed on the same
/// test. On a soft-masked window every upper-case replacement base therefore
/// mismatches a lower-case reference base it is in fact identical to, so
/// `best_alignment` finds no interior matches and a block that splits into
/// members anywhere else stays one spanning `delins` here: the same variant
/// gets a different canonical string purely because its region happens to be
/// soft-masked.
///
/// Folding case at the single point where provider bytes enter is both the
/// smallest fix and the only one that cannot be forgotten by a later comparison
/// added downstream. It loses nothing: the two places that already had to cope
/// with case (`stated_reference_bases_match`, `is_tandem_duplication`) compare
/// with `eq_ignore_ascii_case` and are unaffected, and upper-case is what the
/// rebuilt members must carry anyway — `Base::from_char` does not accept a
/// lower-case byte, so a lower-case base reaching `rebuild_members` refused the
/// canonicalization outright.
fn fetch_canonical_window<P: ReferenceProvider>(
    provider: &P,
    accession: &str,
    start0: i64,
    end0: i64,
) -> Option<Vec<u8>> {
    debug_assert!(end0 > start0, "canonical window must be non-empty");
    let read = |end: i64| -> Option<String> {
        let seq = provider
            .get_sequence(accession, start0 as u64, end as u64)
            .ok()?;
        (!seq.is_empty()).then_some(seq)
    };
    let seq = match read(end0) {
        Some(seq) => seq,
        None => {
            let length = i64::try_from(provider.get_sequence_length(accession).ok()?).ok()?;
            let clamped = end0.min(length);
            if clamped <= start0 {
                return None;
            }
            read(clamped)?
        }
    };
    // Uppercased in place rather than via `to_ascii_uppercase()`, which would
    // allocate and copy a second window-sized buffer on top of the one the
    // provider just built.
    let mut bytes = seq.into_bytes();
    bytes.make_ascii_uppercase();
    Some(bytes)
}

/// Which edges of the fetched window are edges of the sequence itself.
///
/// [`fetch_canonical_window`] serves 0-based `[start0, start0 + len)` of
/// `accession`, and a window is padded on both sides, so an edge is the
/// sequence's own only when the padding was clipped there. The 5' side is
/// decidable from `start0` alone; the 3' side needs the sequence's length, which
/// is the single provider round-trip this costs — and the caller only pays it
/// when a derived piece is actually resting on an edge.
///
/// A length the provider cannot answer yields [`SequenceEnds::INTERIOR`]: not
/// knowing where the sequence ends is a reason to leave a derived insertion
/// alone, never a reason to rewrite one.
fn window_sequence_ends<P: ReferenceProvider>(
    provider: &P,
    accession: &str,
    start0: i64,
    window_len: usize,
) -> SequenceEnds {
    let three_prime = provider
        .get_sequence_length(accession)
        .ok()
        .and_then(|len| i64::try_from(len).ok())
        .is_some_and(|len| start0 + window_len as i64 == len);
    SequenceEnds {
        five_prime: start0 == 0,
        three_prime,
    }
}

/// Apply every edit to the window, returning the resulting bytes.
///
/// Returns `None` if an edit falls outside the window, a stated reference base
/// cannot be read, or **two edits claim the same reference position**. That last
/// case is the #1234 corruption: an allele whose members overlap has no
/// well-defined resulting sequence — applying it depends on member order — so
/// refusing is the only honest answer. It is also what stops this pass from
/// laundering an already-corrupt allele into a plausible-looking one.
fn apply_edits_to_window(edits: &[GEdit], ref_bytes: &[u8], w_lo: i64) -> Option<Vec<u8>> {
    // Refuse a descending span before any arithmetic touches it.
    //
    // `e < s` means the member wraps the origin of a circular genome — the `m.`
    // axis admits `m.16569_1del`. This window is linear, so such a member cannot
    // be applied, and every way of trying is silently wrong:
    //
    // * `Dup` sizes a `Vec` from `(e - s + 1) as usize`, and `(1 - 16569 + 1)`
    //   reinterprets as ~1.8e19 — `Vec::with_capacity` then **aborts the
    //   process**, in release as well as debug, on the valid input
    //   `NC_012920.1:m.[16569_1dup;5A>G]`. That is reachable from
    //   `Normalizer::normalize`, so it is a denial-of-service surface for the
    //   `web-service` build, not merely a test-time panic.
    // * `Del`/`Inv`/`Delins` iterate `s..=e` and call `claim(.., s, e)`, both of
    //   which are *empty* ranges when `e < s`. The member claims nothing,
    //   changes nothing, and `rebuild_members` then emits the allele with it
    //   **deleted** — `m.[16569_1del;5A>G]` came back as plain `m.5A>G`.
    //
    // `edit_span_union` only catches this when the wraparound member is alone
    // (its `lo <= hi` check fails); paired with any sibling the bogus union
    // passes. Refusing here sends the allele to the per-member pipeline, which
    // applies the `is_wraparound_genome` guard that `cis_axis_parts` does not.
    for edit in edits {
        let descending = match edit {
            GEdit::Del { s, e }
            | GEdit::Dup { s, e }
            | GEdit::Inv { s, e }
            | GEdit::Delins { s, e, .. } => e < s,
            // A repeat does not denote a sequence until `lower_repeat_edits`
            // has resolved its tract, so there is nothing here to apply. Both
            // callers lower first; refusing rather than approximating keeps
            // that a precondition instead of a silent wrong answer.
            GEdit::Repeat { .. } => return None,
            GEdit::Ins { .. } | GEdit::Sub { .. } => false,
        };
        if descending {
            return None;
        }
    }

    // An insertion — or a duplication's copy — landing *strictly inside*
    // territory a sibling deletes, replaces, or inverts is the #486 overlap
    // conflict: `c.[4_10inv;5_6insAA]` says nothing about whether the inserted
    // bases are themselves inverted, so the combination has no defined
    // resulting sequence. Ferro's standing answer since #395/#486/#1004 is to
    // report `OverlapConflict` (an error in strict mode) and leave the
    // description alone rather than pick a winner. This pass would instead
    // apply both edits in member order and hand back a tidy single edit that no
    // longer looks conflicted at all (`c.5_9delinsCAAAATT`) — laundering the
    // conflict out of the string while the warning still says there is one.
    //
    // The interior test is `overlap::detect_insertion_overlaps`' own — the
    // junction after `gap` is interior to `[s, e]` when `s <= gap < e` — so
    // this refuses exactly the set that detector reports and no more. An
    // insertion merely *abutting* a span's edge is not interior and stays
    // accepted, which is what keeps the spec-valid shape that detector
    // documents (an insertion on each side of a substitution) flowing through.
    //
    // Keyed off the edit geometry rather than off the `OverlapConflict` warning
    // the earlier layers raise, deliberately: the geometry is a property of the
    // edits themselves and so survives any respelling of them, whereas the
    // diagnostic has to be re-derived by a detector on each pass. That
    // distinction is not hypothetical — the per-member pipeline rewrites
    // `[4_10inv;5_6insAA]` into `[5_9inv;5_6dup]`, and a detector that keys on
    // `NaEdit::Insertion` alone stops seeing the conflict once the interior
    // `ins` has become a `dup`. A warning-driven refusal would then decline on
    // the first pass and accept on the second, costing idempotency
    // (`FERRO_ASSERT_IDEMPOTENT` catches it).
    //
    // `detect_insertion_overlaps` now registers `dup` and `repeat` as junction
    // occupants too, exactly as the `GEdit::Dup { e, .. } => *e` arm below does,
    // so the two agree on those shapes. Keying off geometry keeps them from
    // drifting apart again.
    let n = ref_bytes.len();
    let idx = |p: i64| -> Option<usize> {
        let i = p - w_lo;
        (i >= 0 && (i as usize) < n).then_some(i as usize)
    };
    // `interior_gap[i]` marks the junction *after* window position `i` as
    // enclosed by some span edit. Marked as a bitmap in one pass rather than
    // re-scanned per insertion, mirroring `claimed` below: a 4096-wide window
    // may carry as many members, and the pairwise form is quadratic in that.
    let mut interior_gap = vec![false; n];
    for edit in edits {
        let (s, e) = match edit {
            GEdit::Del { s, e } | GEdit::Delins { s, e, .. } | GEdit::Inv { s, e } => (*s, *e),
            // Unreachable: the descending-span loop above returns `None` on
            // any unlowered repeat before this one runs.
            GEdit::Ins { .. } | GEdit::Dup { .. } | GEdit::Sub { .. } | GEdit::Repeat { .. } => {
                continue
            }
        };
        for p in s..e {
            if let Some(i) = idx(p) {
                interior_gap[i] = true;
            }
        }
    }
    for edit in edits {
        let gap = match edit {
            GEdit::Ins { gap, .. } => *gap,
            // A `Dup` writes its copy into the slot after its own last base.
            GEdit::Dup { e, .. } => *e,
            GEdit::Del { .. }
            | GEdit::Delins { .. }
            | GEdit::Inv { .. }
            | GEdit::Sub { .. }
            // Unreachable for the same reason as the loop above.
            | GEdit::Repeat { .. } => continue,
        };
        if idx(gap).is_some_and(|i| interior_gap[i]) {
            return None;
        }
    }

    let mut cell: Vec<Option<u8>> = ref_bytes.iter().map(|b| Some(*b)).collect();
    let mut after: Vec<Vec<u8>> = vec![Vec::new(); n];
    let mut claimed = vec![false; n];
    // Claim `[s, e]` for one edit, refusing if another already holds any of it.
    let claim = |claimed: &mut Vec<bool>, s: i64, e: i64| -> Option<()> {
        for p in s..=e {
            let i = idx(p)?;
            if std::mem::replace(&mut claimed[i], true) {
                return None;
            }
        }
        Some(())
    };

    for edit in edits {
        match edit {
            GEdit::Del { s, e } => {
                claim(&mut claimed, *s, *e)?;
                for p in *s..=*e {
                    cell[idx(p)?] = None;
                }
            }
            GEdit::Sub { pos, alt } => {
                claim(&mut claimed, *pos, *pos)?;
                cell[idx(*pos)?] = Some(canonical_base_byte(alt.to_u8()));
            }
            GEdit::Delins { s, e, alt } => {
                claim(&mut claimed, *s, *e)?;
                for p in *s..=*e {
                    cell[idx(p)?] = None;
                }
                // Same collision guard as `Ins`/`Dup` below. A `Delins` writes
                // its replacement into the `after` slot of its first base, so an
                // insertion at that same gap (or a duplication landing there)
                // competes for one slot — and without this check the outcome
                // depended on member order: `g.[5_6delinsTT;5insAAA]` was
                // accepted as `AAATT` when the insertion was applied first and
                // refused when the delins was. One input, two answers, which is
                // precisely the non-confluence this pass exists to remove.
                let slot = idx(*s)?;
                if !after[slot].is_empty() {
                    return None;
                }
                after[slot].extend(alt.iter().map(|b| canonical_base_byte(b.to_u8())));
            }
            GEdit::Ins { gap, alt } => {
                // An insertion occupies no reference position, so it claims
                // none; two insertions at one gap are caught below.
                let slot = idx(*gap)?;
                if !after[slot].is_empty() {
                    return None;
                }
                after[slot].extend(alt.iter().map(|b| canonical_base_byte(b.to_u8())));
            }
            GEdit::Dup { s, e } => {
                let mut copied = Vec::with_capacity((*e - *s + 1) as usize);
                for p in *s..=*e {
                    copied.push(ref_bytes[idx(p)?]);
                }
                let slot = idx(*e)?;
                if !after[slot].is_empty() {
                    return None;
                }
                after[slot].extend(copied);
            }
            // Unreachable — refused by the descending-span loop at the top —
            // and restated here so a future edit that removes that loop cannot
            // silently start applying a repeat as if it denoted nothing.
            GEdit::Repeat { .. } => return None,
            GEdit::Inv { s, e } => {
                claim(&mut claimed, *s, *e)?;
                let span: Vec<u8> = (*s..=*e)
                    .map(|p| idx(p).map(|i| ref_bytes[i]))
                    .collect::<Option<_>>()?;
                for (offset, base) in span.iter().rev().enumerate() {
                    cell[idx(*s + offset as i64)?] = Some(complement_base(*base)?);
                }
            }
        }
    }

    // A `Delins` writes its replacement into `after[start]` while blanking the
    // whole span, so the replacement lands 3' of the (now deleted) first base —
    // the same position it occupies in the reference.
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        if let Some(b) = cell[i] {
            out.push(b);
        }
        out.extend(&after[i]);
    }
    Some(out)
}

/// Fold a base byte into the single alphabet used for comparison: upper case,
/// with `U` mapped to `T`.
///
/// The `r.` axis submits RNA bases while its reference is the DNA transcript
/// sequence, so uracil and thymine denote the same nucleotide and must compare
/// equal — otherwise every `u` in an `r.` description reads as a change.
fn canonical_base_byte(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    }
}

/// Trim the common prefix and suffix of the reference and result blocks.
///
/// Returns `(prefix_len, ref_end, alt_end)` as offsets into the two slices.
fn trim_common_flanks(reference: &[u8], result: &[u8]) -> (usize, usize, usize) {
    let mut lo = 0;
    while lo < reference.len() && lo < result.len() && reference[lo] == result[lo] {
        lo += 1;
    }
    let mut hi_ref = reference.len();
    let mut hi_alt = result.len();
    while hi_ref > lo && hi_alt > lo && reference[hi_ref - 1] == result[hi_alt - 1] {
        hi_ref -= 1;
        hi_alt -= 1;
    }
    (lo, hi_ref, hi_alt)
}

/// Partition a changed block into maximal runs of change separated by at least
/// one unchanged base.
///
/// Equal-length blocks compare position-wise — there is no alignment choice to
/// make, which is why #1230/#1231/#1233 have a single answer. Unequal-length
/// blocks do have a choice: the minimal edit set is *not* unique (#1232 has four
/// equally minimal, equally stable encodings). We pick the alignment that
/// maximizes matched bases — the separation rule (`delins.md:17`) wants the
/// pieces — and break ties by placing the indel 5'-most. The tie-break is an
/// implementer's choice; the spec does not reach it. Each piece is 3'-shifted
/// afterwards, so the 3' rule is still honoured per member.
fn partition_block(reference: &[u8], result: &[u8]) -> Vec<Piece> {
    let whole = || {
        vec![Piece {
            ref_start: 0,
            ref_end: reference.len(),
            alt: result.to_vec(),
        }]
    };
    // Scoped to length-changing blocks, because that is the only regime the cap
    // reaches at all. An equal-length block has no gap to place, so
    // `best_alignment` returns the position-wise pairing on its first statement
    // and the whole partition is linear — there is no search here to bound.
    //
    // Note what the cap does *not* do, even where it applies: it does not bound
    // a quadratic cost. Scoring one placement is O(1) because the score is
    // separable, which makes `best_alignment`'s search linear, and the one
    // quadratic step — ranking a tie by what it separates — carries its own and
    // tighter bound in `MAX_TIE_BREAK_SWEEP` (256), which binds first.
    //
    // Capping it anyway cost confluence rather than time: the un-partitioned
    // whole block is refused by the weight bound whenever the input was spelled
    // as its individual changes, so a 1100 nt near-palindrome left
    // `g.257_1356inv` and `g.[257C>A;267A>C;1346G>T;1356T>G]` both stable — one
    // variant, two normal forms, with an exact boundary at the cap (1024
    // confluent, 1026 not) that gave the length short-circuit away.
    if reference.len() != result.len()
        && (reference.len() > MAX_SPLIT_BLOCK || result.len() > MAX_SPLIT_BLOCK)
    {
        return whole();
    }
    let Some(columns) = best_alignment(reference, result) else {
        return whole();
    };
    let pieces = pieces_from_columns(&columns, reference, result);
    if pieces.is_empty() {
        return whole();
    }
    // Every length-changing block, not just net insertions (#1271). Equal-length
    // blocks are exempt because `best_alignment` compares position-wise there —
    // there is no gap to place, so no search to seize on a coincidental match.
    if reference.len() != result.len()
        && !separations_are_meaningful(&pieces, result.len().abs_diff(reference.len()))
    {
        return whole();
    }
    // Scoped to length-changing blocks: an equal-length block has no gap to
    // place, so every matched base is a genuine coordinate-wise identity rather
    // than an artefact of where the gap landed, and a split across one is real.
    // `separations_are_meaningful` above draws the same line.
    if reference.len() != result.len()
        && pieces.len() > 1
        && every_separation_is_a_single_base(&pieces)
        && split_buys_no_higher_priority_type(&pieces, reference)
    {
        return whole();
    }
    pieces
}

/// Merge a run of pieces that together span one inversion back into that single
/// inversion.
///
/// An inversion is defined by `inversion.md:5` as a span replaced by its reverse
/// complement, and that is a property of the *whole* span. Nothing forbids some
/// of its columns coinciding: `ACGTGC -> GCACGT` is a textbook 6 nt inversion
/// whose 2nd and 5th columns happen to match, so a position-wise partition finds
/// `[A>G, GT>AC, C>T]` in it and the inversion disappears — which is what the
/// since-removed `needs_unsupported_form` used to prevent, by refusing the
/// derivation outright.
///
/// This is the same reading [`crate::normalize::rules`]'s single-span typing
/// already applies (a shared-affix-trimmed span that is a whole reverse
/// complement is an `inv`), so putting the block back together here is what
/// makes the derivation agree with the per-member pipeline instead of shredding
/// what that pipeline types correctly.
///
/// Gated on [`every_separation_is_a_single_base`], and that gate is load-bearing
/// rather than defensive. `AAGCTA -> TAGCTT` is also a whole reverse complement,
/// but only its first and last bases change and four unchanged bases lie between
/// them. `general.md:34` describes two variants separated by one or more
/// nucleotides individually; its letter names only `delins`, but its rationale
/// (`delins.md:81-84` — the variants may have been reported, or may occur,
/// individually) reaches any single spanning description, and at four unchanged
/// bases these are two independent changes rather than interior columns of one
/// reverse-complement relation. At exactly one unchanged base the match is far
/// more likely to be the alignment's coincidence than structure — the same line
/// [`every_separation_is_a_single_base`] draws inside [`partition_block`].
///
/// Expressed through [`is_inversion`] rather than re-deriving the test, so the
/// block-level and piece-level answers cannot drift apart.
///
/// # Why this runs *after* the changed-columns weight bound
///
/// Load-bearing placement, not an accident of ordering. This rule **widens**:
/// the merged piece claims every column between the first and the last, including
/// the ones the partition found unchanged. The weight bound compares the derived
/// weight against *the input's*, so a widening applied before it is judged
/// against a quantity that differs between two spellings of one variant — the
/// 5 nt inversion `GTTAA -> TTAAC` weighs 3 spelled as `g.[257G>T;259T>A;261A>C]`
/// and 5 spelled as `g.257_261delinsTTAAC`, so a pre-bound widening to 5 columns
/// is accepted for the second and refused for the first. That is exactly the
/// non-confluence #1235 exists to remove, and it would be *introduced* here.
/// (Offsets are the synthetic reference's, which puts the core at 257, so the
/// case is quotable verbatim from
/// `issue_1040_inv_overrecognition_probes::every_spelling_of_a_derived_whole_block_inversion_converges_on_inv`.)
///
/// Applied afterwards, the bound judges the un-widened partition — whose weight
/// is what it was before this rule existed — and the widening is a licensed
/// re-typing on top of an already-accepted partition. That is the same argument
/// [`apply_coding_codon_exception`] is placed by.
///
/// Two citations carry the licence, and they answer different questions.
///
/// The **widening** — collapsing changed columns back into one spanning
/// description when the only thing between them is a base that happens to match
/// — is what `delins.md:46-47` recommends. Against
/// `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` it notes that "parts of the
/// inserted sequence 'align' with the reference sequence, giving an alternative
/// description like `c.[850_869del;874_881del;887_897del;901_902insG]`", and
/// answers: "**The 'delins' format is recommended**: it is simpler and prevents
/// software tools making incorrect predictions for the consequences on protein
/// level." That is this shape exactly — interior columns coinciding by alignment
/// rather than by structure — and it endorses the spanning form over the split.
/// [`every_separation_is_a_single_base`] cites the same passage for the regime
/// it admits, so both ends of this rule rest on one line of the spec.
///
/// The **type** of that spanning description is then not a preference but a
/// definition: `delins.md:5` defines a delins as a replacement "which is not a
/// substitution or inversion", so the very form `delins.md:46-47` recommends may
/// not be spelled `delins` over a reverse-complement span, and `inversion.md:5`
/// admits the `inv` on sequence alone. Between them the widening is recommended
/// and the resulting type is forced; reaching the `inv` describes the same change
/// under the only spanning type the spec allows for it — not more change.
///
/// Deliberately **not** `general.md:56`. That list ranks single-variant type
/// labels for one span; it never ranks a multi-member allele against a spanning
/// description, it omits `delins` altogether, and it places substitution *above*
/// inversion — so if it reached this comparison at all it would argue the other
/// way.
///
/// The general rule, worth stating because it caught this once: **any partition
/// rule that widens a piece must run after the weight bound, never inside
/// [`partition_block`].**
fn coalesce_whole_block_inversion(pieces: &mut Vec<Piece>, ref_bytes: &[u8]) {
    if pieces.len() < 2
        || !(every_separation_is_a_single_base(pieces)
            || changed_columns_dominate_the_span(pieces)
            || no_piece_is_a_lone_substitution(pieces))
    {
        return;
    }
    let start = pieces[0].ref_start;
    let end = pieces[pieces.len() - 1].ref_end;
    // Defence in depth, not a live case, and checked **before** anything is
    // allocated or sliced. The pieces reaching here are disjoint and ascending:
    // `coalesce_adjacent_pieces` asserts it, and the only step between that call
    // and this one is `shrink_pieces_to_differences`, which can narrow a piece
    // but never widen one. But that assertion is a `debug_assert!` and so is
    // compiled out of release, and the gate above cannot substitute for it: it
    // measures separations with `saturating_sub`, so a strict overlap
    // (`ref_start < ref_end` of its predecessor) reads as `0 <= 1` and *passes*.
    //
    // Two things downstream would then panic inside a library rather than
    // decline: `end - start` underflows if the pieces are out of order at all,
    // and the reconstruction slice indexes with `start > end`. Both are why this
    // runs first. Declining instead leaves the un-coalesced pieces, which the
    // downstream round-trip guard (`reapplied == result`) is already able to
    // refuse.
    if end < start
        || pieces
            .windows(2)
            .any(|pair| pair[1].ref_start < pair[0].ref_end)
    {
        return;
    }
    // Reconstruct what the span denotes: each piece's payload, with the
    // unchanged reference bases separating them spliced back in.
    let mut alt = Vec::with_capacity(end - start);
    let mut cursor = start;
    for piece in pieces.iter() {
        alt.extend_from_slice(&ref_bytes[cursor..piece.ref_start]);
        alt.extend_from_slice(&piece.alt);
        cursor = piece.ref_end;
    }
    let whole = Piece {
        ref_start: start,
        ref_end: end,
        alt,
    };
    if is_inversion(&whole, ref_bytes) {
        *pieces = vec![whole];
    }
}

/// Whether more of the span is changed than is left unchanged.
///
/// The second admission route into [`coalesce_whole_block_inversion`], and the
/// one that recognises a long inversion (#1461). It is deliberately **not** a
/// wider separation threshold, because separation cannot decide this case:
///
/// | | span | changed | widest gap |
/// |---|---|---|---|
/// | `NC_000013.10:g.100809575_100810031inv` (#1461) | 457 nt | 307 (67.2%) | **5** |
/// | `AAGCTA -> TAGCTT` (the #1040 control) | 6 nt | 2 (33.3%) | **4** |
///
/// The inversion that must be recognised has *wider* gaps than the pair that
/// must stay individual, so every threshold on separation either refuses #1461
/// or admits the #1040 control. That is why widening
/// [`every_separation_is_a_single_base`] breaks
/// `a_whole_span_reverse_complement_is_not_merged_across_a_multi_base_separation`.
///
/// Density does separate the two rows of that table, and this rule is what
/// admits #1461. What it is *not* is a test that the reverse-complement relation
/// is real, and an earlier revision of this comment claimed it was. That claim
/// is corrected here rather than deleted, because the mistake in it is an easy
/// one to make again.
///
/// # Why an unchanged column count is weak evidence at small `n`
///
/// Unchanged columns are not independent, and they do not arrive one at a time.
/// Column `i` of a whole-span reverse complement coincides exactly when column
/// `n-1-i` does — both say `b[i] == complement(b[n-1-i])` — so the coincidences
/// come in **mirror pairs**, and an odd `n` has a centre column that can never
/// coincide. Both halves of that hold over a uniform A/C/G/T alphabet and not
/// beyond it: [`crate::sequence::complement_base`] maps the self-complementary
/// IUPAC codes to themselves (`N -> N`, `S -> S`, `W -> W`), so a centre column
/// holding one of those *does* coincide, and real reference sequences carry `N`
/// at assembly gaps. For a uniform random block the count of unchanged columns
/// is therefore
///
/// ```text
/// m = floor(n/2)                         the number of mirror pairs
/// U = 2 · Binomial(m, 1/4)               E[U] = m/2      sd(U) = sqrt(3m/4)
/// ```
///
/// Stated in `m`, not `n`, because the centre column of an odd span is not part
/// of a pair: for even `n` the mean is `n/4` — a quarter of the span — but for
/// odd `n` it is `(n−1)/4`, and the same `−1` runs through the standard
/// deviation. An earlier revision of this comment gave only the even-`n` forms
/// and asserted them for "every `n`"; #1461's span is 457, so the mean it
/// implies is 114.25 rather than 114 and its `z` comes out at `+2.73` rather
/// than the `+2.75` tabulated below. The table was computed correctly and only
/// the prose over-generalised, which is exactly the failure #1461 is about — a
/// statistic stated more broadly than it holds, then cited as evidence.
/// [`the_unchanged_column_moments_reproduce_the_tabulated_z_scores`] now pins
/// the two against each other so the generalisation cannot recur silently.
///
/// The other misreading, which the earlier revision also made, is to treat the
/// mean as if it were the whole distribution and conclude that an observation
/// above it is "far above chance". Measured against the actual distribution:
///
/// | block | span | unchanged | z | P(at least this many) |
/// |---|---|---|---|---|
/// | `AAGCTA -> TAGCTT` (the #1040 control) | 6 | 4 (66.7%) | +1.67 | **15.6%** |
/// | `AATGCACA -> TGTGCATT` (#1517) | 8 | 4 (50.0%) | +1.15 | **26.2%** |
/// | `NC_000013.10:g.100809575_100810031inv` (#1461) | 457 | 150 (32.8%) | +2.75 | **0.45%** |
///
/// So roughly **one in six** genuine 6 nt inversions leaves two thirds of its
/// columns unchanged, and roughly **one in four** genuine 8 nt inversions leaves
/// half of them unchanged. Neither is a near-palindrome signature; both are
/// ordinary draws. Only #1461's is rare, and it is rare because the span is 457
/// nt, not because 32.8% is a low fraction.
///
/// The discriminating quantity is `n`. The distribution's relative spread is
/// `sd/mean = sqrt(3/m)`, so it is 100% of the mean at `n = 6` and 11.5% at
/// `n = 457`: a density reading carries almost no information about a short
/// block and a great deal about a long one. Any future rule that wants to argue
/// from coincidence rates must scale with `n`; this one does not.
///
/// What this predicate honestly is, then, is a **length-correlated proxy** that
/// separates the two cases in the table and is calibrated on them. It is not
/// evidence about whether a given reverse-complement relation is structural, and
/// no comment here should say that it is.
///
/// Additive by construction: this is an `||` alternative to the single-base gate,
/// so every block that coalesces today still coalesces. Only a block that is
/// currently *refused* can change, which is what confines the blast radius to
/// the class #1461 reports.
fn changed_columns_dominate_the_span(pieces: &[Piece]) -> bool {
    let (start, end) = (pieces[0].ref_start, pieces[pieces.len() - 1].ref_end);
    let Some(span) = end.checked_sub(start).filter(|span| *span > 0) else {
        return false;
    };
    let changed: usize = pieces
        .iter()
        .map(|piece| piece.ref_end.saturating_sub(piece.ref_start))
        .sum();
    changed * 2 > span
}

/// Whether every gap between consecutive pieces is a single unchanged base —
/// the regime where a match is most likely to be coincidence rather than
/// structure.
///
/// This is what keeps [`split_buys_no_higher_priority_type`] from overriding
/// `general.md:34` outright. That rule — "two variants separated by one or more
/// nucleotides should be described individually and **not** as a 'delins'" — is
/// unambiguous once there is real distance between the members, and an
/// all-`delins` split at, say, two unchanged bases is one the rule plainly
/// wants: `AGAGT -> GAAGAG` divides into `AG>GA` and `T>AG` for 4 changed
/// columns against the spanning form's 6, so collapsing it would be both less
/// minimal and against the rule.
///
/// At exactly one base the two readings genuinely compete, which is why #422 is
/// a *deliberately pinned* non-confluence rather than a plain bug, and why
/// `delins.md:44-47` prefers the spanning form where a payload merely *aligns*.
/// Restricting the collapse to that regime is what lets both rules stand.
fn every_separation_is_a_single_base(pieces: &[Piece]) -> bool {
    pieces
        .windows(2)
        .all(|pair| pair[1].ref_start.saturating_sub(pair[0].ref_end) <= 1)
}

/// No piece spells a lone substitution — every one covers two or more columns.
///
/// The third alternative admitting [`coalesce_whole_block_inversion`], and the
/// one that closes #1517. It exists because the other two are geometric — they
/// ask how far apart the pieces are, or what fraction of the span they cover —
/// and geometry cannot separate #1517 from the case that must stay split. The
/// two blocks are the same shape:
///
/// ```text
/// #1230   GATG     -> CATC       changes at 2 of 4 columns, interior `AT` coincides
/// #1517   AATGCACA -> TGTGCATT   changes at 4 of 8 columns, interior `TGCA` coincides
/// ```
///
/// Both are whole-span reverse complements whose interior columns coincide, so
/// `delins.md:44-47`'s carve-out — prefer the spanning description when the
/// alternative exists only because parts of the payload align with the reference
/// — reaches both. What differs is the **type of the members the spanning form
/// competes with**, and that is the one axis `general.md:56` speaks to:
///
/// | | competing members | rank vs inversion (3) |
/// |---|---|---|
/// | #1230 | two substitutions | substitution is (1) — **above** |
/// | #1517 | two `delins` | `delins` is absent from the list — **not above** |
///
/// So `:44-47` licenses the widening for both, and `:56` withdraws it for #1230
/// only. A run of one column is a substitution by `delins.md:15`; a run of two or
/// more is a `delins` by `:16`. This predicate is that distinction, expressed on
/// the pieces.
///
/// # This reading is contested, and the contest is recorded rather than hidden
///
/// [`coalesce_whole_block_inversion`]'s own comment argues `general.md:56` must
/// **not** reach this comparison — that it "ranks single-variant type labels for
/// one span" and "never ranks a multi-member allele against a spanning
/// description". That objection is real and is not answered here. Two things
/// weigh against it: the 208-row net-insertion arbitration used `:56` in exactly
/// this way to conclude that a split into `delins`/insertion members is
/// unsupported, and the alternative — applying `:44-47` with no type test at all
/// — merges #1230 too, contradicting a guard filed against that very behaviour.
///
/// The honest statement is that ranking descriptions of differing **arity** by
/// member type is an implementer's reading of `:56`, not something its text
/// compels. `tests/it/issue_1517_inv_priority_over_delins.rs` records the
/// question, both readings, and the choice.
fn no_piece_is_a_lone_substitution(pieces: &[Piece]) -> bool {
    pieces
        .iter()
        .all(|piece| piece.ref_end.saturating_sub(piece.ref_start) >= 2)
}

/// Partition a changed block by deriving member boundaries from the **denoted
/// sequence**, via [`crate::normalize::seqfirst`].
///
/// Same contract as [`partition_block`] — pieces are disjoint, ascending, and
/// reconstruct `result` when substituted into `reference` — but derived from
/// the alignment steps common to *every* minimal alignment rather than from a
/// single-gap search. Two spellings of one variant produce the same
/// `(reference, result)` pair and therefore the same pieces.
///
/// `reference` and `result` must already have their common flanks trimmed;
/// `canonicalize_from_sequence` does that with [`trim_common_flanks`] before
/// calling in. An untrimmed flank changes the partition, so this is a
/// precondition and not a convenience.
///
/// Returns `None` when the block cannot be partitioned this way, which today
/// means only that [`member_alt_spans`] declined. The caller falls back.
///
/// `min_separation` is the unchanged-base threshold `partition_members_with`
/// merges below (see that function's doc). It is **not** hardcoded to
/// `seqfirst::MIN_SEPARATION` — the caller must choose it per axis. See
/// `MIN_SEPARATION_NO_FRAME` for why: the value that is correct on a
/// reading-frame axis is wrong on a genomic one.
///
/// Called only from the `FERRO_SEQFIRST_SHADOW` comparison in
/// [`canonicalize_from_sequence`]; it does not yet decide any result.
fn partition_block_sequence_first(
    reference: &[u8],
    result: &[u8],
    min_separation: u32,
) -> Option<Vec<Piece>> {
    use crate::normalize::seqfirst::align::AlignmentDag;
    use crate::normalize::seqfirst::partition::{member_alt_spans, partition_members_with};

    // Refuse an oversized grid before building it, not after: the allocation is
    // the cost. Checked, because `(n + 1) * (m + 1)` on two `usize` lengths is
    // itself an overflow site.
    let cells = reference
        .len()
        .checked_add(1)?
        .checked_mul(result.len().checked_add(1)?)?;
    if cells > MAX_SEQFIRST_GRID_CELLS {
        return None;
    }

    let dag = AlignmentDag::build(reference, result);
    let dominators = dag.dominators();
    let members = partition_members_with(&dag, &dominators, min_separation);
    let spans = member_alt_spans(&dag, &dominators, &members)?;

    Some(
        members
            .iter()
            .zip(spans)
            .map(|(member, (alt_start, alt_end))| Piece {
                ref_start: member.ref_start as usize,
                ref_end: member.ref_end as usize,
                alt: result[alt_start as usize..alt_end as usize].to_vec(),
            })
            .collect(),
    )
}

/// Partition a changed block by the **canonical** rule: the member-count-minimal
/// minimal alignment, via
/// [`crate::normalize::seqfirst::partition::CanonicalAlignment`].
///
/// Same contract as [`partition_block`] — pieces are disjoint, ascending, and
/// reconstruct `result` when substituted into `reference` — and the same
/// precondition as [`partition_block_sequence_first`]: `reference` and `result`
/// must already have their common flanks trimmed.
///
/// Takes no `min_separation`. The canonical rule describes exactly what one
/// alignment changes, so merging runs across a short unchanged gap would make it
/// claim more changed columns than the block's edit distance —
/// `partition_members_canonical`'s doc works the case through.
///
/// Returns `None` only when the alignment grid would exceed
/// [`MAX_SEQFIRST_GRID_CELLS`]; the caller falls back to [`partition_block`].
///
/// Reached only under `FERRO_PARTITION=canonical` and from the `dump_partitions`
/// example; it does not decide any shipped result.
fn partition_block_canonical(reference: &[u8], result: &[u8]) -> Option<Vec<Piece>> {
    use crate::normalize::seqfirst::align::AlignmentDag;
    use crate::normalize::seqfirst::partition::CanonicalAlignment;

    // Refuse an oversized grid before building it, for the same reason
    // `partition_block_sequence_first` does: the allocation is the cost, and
    // `(n + 1) * (m + 1)` is itself an overflow site.
    let cells = reference
        .len()
        .checked_add(1)?
        .checked_mul(result.len().checked_add(1)?)?;
    if cells > MAX_SEQFIRST_GRID_CELLS {
        return None;
    }

    let dag = AlignmentDag::build(reference, result);
    let canonical = CanonicalAlignment::of(&dag);
    Some(
        canonical
            .members()
            .iter()
            .zip(canonical.alt_spans())
            .map(|(member, (alt_start, alt_end))| Piece {
                ref_start: member.ref_start as usize,
                ref_end: member.ref_end as usize,
                alt: result[alt_start as usize..alt_end as usize].to_vec(),
            })
            .collect(),
    )
}

/// Whether `result` embeds in `reference` as a **subsequence**, i.e. whether the
/// whole block can be described by deletion alone.
///
/// This is `delins.md:44-47`'s "parts of the inserted sequence *align* with the
/// reference sequence" stated exactly: when it holds, every base of the payload
/// can be read off the reference in order, so a partitioner that cuts at those
/// coincidences produces extra members that exist **only** because payload bases
/// happen to match reference bases between them. The spec's answer to that shape
/// is unambiguous — "The 'delins' format is recommended".
///
/// # Why a two-pointer walk and not the counting DP
///
/// The adjudication's predicate is `count_minimal_embeddings(ref, payload).1 > 0`,
/// where the DP returns `(minimum deletion runs, number of embeddings achieving
/// it)`. The count is zero exactly when no embedding exists, so **for the
/// predicate the DP is equivalent to a subsequence test** and this runs in `O(n)`
/// with no allocation instead of `O(n·m)`. The DP is still the right tool for the
/// *non-uniqueness* count (how many equally-compliant splits exist), which this
/// pass does not need — it merges regardless of whether the split it replaces was
/// unique.
///
/// Vacuously true for an empty `result`: a wholly deleted block is one member
/// already, so the caller's member-count gate is what declines it.
fn payload_embeds_as_subsequence(reference: &[u8], result: &[u8]) -> bool {
    let mut reference = reference.iter();
    result
        .iter()
        .all(|base| reference.any(|candidate| candidate == base))
}

/// How many substituted positions the coalesce pass tolerates inside an
/// otherwise coincidental split.
///
/// **One**, and the value is measured rather than chosen. The spec's own worked
/// example, `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`, needs exactly one: its
/// minimal split is four coincidental deletions plus a single substitution, so a
/// pure-subsequence test cannot see it. The dominant churn class needs zero. Two
/// substitutions six bases apart — the case `delins.md:17` requires to stay
/// individual — needs two, and also fails the net-deletion condition, so it is
/// excluded twice over.
///
/// Raising this is a representation change: it merges blocks that are currently
/// spelled as separate members. Measure before moving it.
const COALESCE_MISMATCH_BUDGET: usize = 1;

/// Widest run of unchanged reference bases the coalesce pass will merge across.
///
/// **Eight**, and the value is measured, not chosen.
///
/// `delins.md:44-47`'s own worked example is the lower bound: the canonical
/// arm's split of it leaves a **6**-base gap, so any cap below that rejects the
/// one case the spec states outright. (The 4-member alternative the spec quotes
/// has gaps of 4/5/4; canonical cuts it differently.)
///
/// The upper bound comes from the corpora. Over 500k ClinVar rows the 78 splits
/// this pass correctly rejoins have a maximum separation of **8** (median 2,
/// p90 3). Over the 592 real multi-member cis alleles, the members it must NOT
/// merge sit at `[1, 4, 9, 21, 26, 58, 166, 234, 604, 973, 1562]`. Eight is the
/// widest value that admits every correct rejoin, and it excludes 9 of those 11.
///
/// **The two at 1 and 4 are not separable by separation** and are a known gap:
/// they are genuinely distinct variants that happen to sit as close together as
/// a coincidental split. Distinguishing them needs something this pass does not
/// look at.
///
/// Without this cap the pass merges genuinely separate members of a multi-member
/// cis allele, which `delins.md:17` and `general.md:34` require to stay
/// individual. Measured over the 592 real multi-member cis alleles, the members
/// it wrongly merged sat a median of **58** unchanged bases apart (p90 604, max
/// 1,562), while the 78 coincidental splits it correctly rejoined on ClinVar sat
/// a median of **2** apart (p90 3, max 8). The two populations barely overlap,
/// and separation is what separates them.
const COALESCE_MAX_SEPARATION: usize = 8;

/// Whether `payload` can be read out of `span` by deleting bases and
/// substituting at most `budget` of them.
///
/// `budget == 0` is exactly [`payload_embeds_as_subsequence`], and is dispatched
/// to it: that case is the common one and answers in `O(n)` with no allocation.
///
/// Above zero this is an `O(n·m)` rolling DP over the minimum number of
/// substituted positions across every embedding — the same shape as the
/// adjudication's counting DP, but minimising mismatches instead of counting
/// deletion runs. An oversized grid falls back to the exact subsequence test
/// rather than declining outright: that is conservative in the safe direction,
/// since it can only make the pass fire *less*.
fn payload_embeds_within_budget(span: &[u8], payload: &[u8], budget: usize) -> bool {
    if budget == 0 {
        return payload_embeds_as_subsequence(span, payload);
    }
    // Each payload base must consume a span base, so a longer payload cannot
    // embed however many substitutions are allowed.
    if payload.len() > span.len() {
        return false;
    }
    let cells = match span
        .len()
        .checked_add(1)
        .and_then(|n| n.checked_mul(payload.len().checked_add(1)?))
    {
        Some(cells) if cells <= MAX_SEQFIRST_GRID_CELLS => cells,
        _ => return payload_embeds_as_subsequence(span, payload),
    };
    let _ = cells;

    // `best[j]` is the fewest substitutions needed to embed `payload[j..]` into
    // the span suffix under consideration. Consuming the whole span with payload
    // left over is impossible, so those cells stay saturated.
    const UNREACHABLE: u32 = u32::MAX / 2;
    let m = payload.len();
    let mut best = vec![UNREACHABLE; m + 1];
    best[m] = 0;
    for &span_base in span.iter().rev() {
        let mut next = vec![UNREACHABLE; m + 1];
        // Deleting every remaining span base is always available and free.
        next[m] = 0;
        for j in (0..m).rev() {
            let deleted = best[j];
            let consumed = best[j + 1].saturating_add(u32::from(span_base != payload[j]));
            next[j] = deleted.min(consumed);
        }
        best = next;
    }
    best[0] as usize <= budget
}

/// Re-spell a multi-member allele as the single `delins` the spec recommends,
/// when the split exists only because payload bases coincide with reference
/// bases between the members.
///
/// This is `DNA/delins.md:44-47` — "parts of the inserted sequence *align* with
/// the reference sequence ... **The 'delins' format is recommended**".
///
/// # It runs last, and that placement is the whole design
///
/// Called after the 3'-shift, the adjacency coalesce and
/// `shrink_pieces_to_differences`, so every earlier pass sees exactly what it
/// would have seen under [`PartitionRule::Canonical`] and this is a terminal
/// re-spelling of an answer already derived from the sequence.
///
/// Running it inside the partitioner instead was measured, and it costs **427
/// converged classes per direction** over the 11,272-class corpus (5':
/// 8,387 -> 7,960, with `split 3` going 0 -> 59 and `split 4+` 0 -> 9). The
/// reason is that the passes above then operate on a merged block for one
/// spelling and an unmerged block for another, so a merge that is itself
/// deterministic still moves two spellings apart.
///
/// That is the separation issue #1430 asks for: step 2 derives a representative
/// variant from the sequence, step 3 makes *that variant* spec-compliant. The
/// coalesce belongs in step 3. (The exon-junction clamp, `general.md:44`, really
/// is a step-2 constraint — an illegal merge destroys a member boundary — but
/// this one destroys nothing and re-spells only.)
///
/// # The two conditions
///
/// **1. The payload embeds as an ordered subsequence of the span.** Then every
/// base of the merged payload can be read off the reference in order, so the
/// members exist only through coincidence. See
/// [`payload_embeds_as_subsequence`].
///
/// **2. The span must lose length** (`payload.len() < span.len()`). Two
/// independent lines of evidence require this and neither is optional:
///
/// * Two plain substitutions six bases apart embed trivially with no deletion,
///   and `delins.md:17` requires those to stay individual.
/// * Measured per spelling, #1421's family runs *opposite* to #1419/#1420's: for
///   a **net insertion** the split form is the canonical one, not the span form.
///   A blanket "prefer fewer members" gets that third of the family backwards.
///
/// Condition 1 already implies `payload.len() <= span.len()`, and equality means
/// the payload *is* the span (an unchanged block), so the strict `<` both
/// enforces net deletion and excludes the degenerate case.
///
/// # Confluence
///
/// Its input is the fully-derived piece list, which under the canonical rule is a
/// function of the denoted sequence rather than of the input's spelling. A
/// deterministic function of a spelling-independent object is spelling
/// independent, so two spellings that converged still converge. The block-level
/// version of this pass did *not* have that property, which is what the 427-class
/// measurement above records.
fn coalesce_payload_alignment_split(pieces: &mut Vec<Piece>, reference: &[u8]) {
    if pieces.len() < 2 {
        return;
    }
    let (start, end) = (pieces[0].ref_start, pieces[pieces.len() - 1].ref_end);
    // Defensive: a malformed piece list is a bug upstream, not something to
    // re-spell. Declining leaves the derived answer untouched.
    if start >= end || end > reference.len() {
        return;
    }

    // What the pieces denote over `[start, end)`: each piece's payload, with the
    // untouched reference between them spliced back in.
    let mut payload = Vec::new();
    let mut cursor = start;
    for piece in pieces.iter() {
        if piece.ref_start < cursor || piece.ref_end < piece.ref_start || piece.ref_end > end {
            return;
        }
        payload.extend_from_slice(&reference[cursor..piece.ref_start]);
        payload.extend_from_slice(&piece.alt);
        cursor = piece.ref_end;
    }
    payload.extend_from_slice(&reference[cursor..end]);

    // Refuse to merge across a wide run of unchanged bases: that is a genuine
    // multi-member allele, not one variant split by coincidence. See
    // `COALESCE_MAX_SEPARATION`.
    if pieces
        .windows(2)
        .any(|pair| pair[1].ref_start.saturating_sub(pair[0].ref_end) > COALESCE_MAX_SEPARATION)
    {
        return;
    }

    let span = &reference[start..end];
    if payload.len() < span.len()
        && payload_embeds_within_budget(span, &payload, COALESCE_MISMATCH_BUDGET)
    {
        *pieces = vec![Piece {
            ref_start: start,
            ref_end: end,
            alt: payload,
        }];
    }
}

/// The four block partitioners, callable by name for measurement.
///
/// **Dev-only measurement surface, not API.** It exists so the `dump_partitions`
/// example can run all three rules over the same blocks and print them side by
/// side; nothing here is part of ferro's supported interface, none of it is
/// covered by semantic versioning, and it may change or vanish with the
/// bake-off it serves. Production code must go through
/// `canonicalize_from_sequence`, which is where the axis, the window and the
/// downstream passes are applied.
///
/// Gated on the `dev` feature, and re-exported as
/// `ferro_hgvs::normalize::dev_partitioners` only there, because the underlying
/// functions are private to this module and should stay that way.
#[cfg(feature = "dev")]
pub mod dev_partitioners {
    use super::Piece;

    /// The default unchanged-base separation threshold for the `shadow` rule:
    /// the coding one-amino-acid exception (`general.md:35`).
    ///
    /// A raw block carries no axis, so a caller measuring blocks in isolation has
    /// to choose. See `MIN_SEPARATION_NO_FRAME` for why the other value is the
    /// right one on an axis with no reading frame.
    pub const DEFAULT_MIN_SEPARATION: u32 = crate::normalize::seqfirst::MIN_SEPARATION;

    /// One derived member of a partitioned block, as 0-based half-open offsets
    /// into the reference block plus its replacement bases.
    ///
    /// A pure insertion has `ref_start == ref_end`; a pure deletion has an empty
    /// `alt`.
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct DevPiece {
        pub ref_start: usize,
        pub ref_end: usize,
        pub alt: Vec<u8>,
    }

    /// Re-shape the module-private `Piece` into the public one.
    ///
    /// A free function rather than a `From` impl on purpose: an
    /// `impl From<&Piece> for DevPiece` would put a module-private type into a
    /// public trait implementation, which is exactly the leak this whole module
    /// exists to avoid.
    fn convert(pieces: &[Piece]) -> Vec<DevPiece> {
        pieces
            .iter()
            .map(|piece| DevPiece {
                ref_start: piece.ref_start,
                ref_end: piece.ref_end,
                alt: piece.alt.clone(),
            })
            .collect()
    }

    /// The shipped rule: `partition_block`'s single-gap alignment search.
    ///
    /// Never declines — an unsplittable block comes back as one spanning member.
    pub fn live(reference: &[u8], result: &[u8]) -> Vec<DevPiece> {
        convert(&super::partition_block(reference, result))
    }

    /// The dominator rule: cut at the steps common to every minimal alignment.
    ///
    /// `None` when the sequence-first splitter declines — an oversized alignment
    /// grid, or an alternate-span derivation that refused.
    pub fn shadow(reference: &[u8], result: &[u8], min_separation: u32) -> Option<Vec<DevPiece>> {
        super::partition_block_sequence_first(reference, result, min_separation)
            .as_deref()
            .map(convert)
    }

    /// The canonical rule: the member-count-minimal minimal alignment.
    ///
    /// `None` only when the alignment grid would be oversized.
    pub fn canonical(reference: &[u8], result: &[u8]) -> Option<Vec<DevPiece>> {
        super::partition_block_canonical(reference, result)
            .as_deref()
            .map(convert)
    }

    /// The block's unit-cost Levenshtein distance — the floor every arm's
    /// changed-column count is measured against.
    pub fn edit_distance(reference: &[u8], result: &[u8]) -> u32 {
        crate::normalize::seqfirst::align::AlignmentDag::build(reference, result).edit_distance()
    }

    /// Changed columns a member set claims, `sum(max(ref_span, alt_len))` — the
    /// same measure `changed_columns_of_pieces` applies inside the normalizer.
    pub fn changed_columns(pieces: &[DevPiece]) -> usize {
        pieces
            .iter()
            .map(|piece| (piece.ref_end - piece.ref_start).max(piece.alt.len()))
            .sum()
    }
}

/// Alignment columns a description marks as changed.
///
/// One member occupies `max(|span|, |replacement|)` columns: a substitution costs
/// 1, an `n`-base deletion or insertion costs `n`, an `n`->`m` `delins` costs
/// `max(n, m)`. Summed over the members this is the count of non-matching columns
/// the description implies — the quantity HGVS asks a description to minimise,
/// and so a search-free bound on how much change a re-derivation may claim.
///
/// Spans are clamped at zero rather than trusted, for the same reason
/// `apply_edits_to_window` refuses them: a wraparound member has `e < s`, and
/// `(e - s + 1) as usize` on those wraps to ~1.8e19.
fn changed_columns_of_edits(edits: &[GEdit]) -> usize {
    fn span(s: i64, e: i64) -> usize {
        (e - s + 1).max(0) as usize
    }
    edits
        .iter()
        .map(|edit| match edit {
            GEdit::Sub { .. } => 1,
            GEdit::Ins { alt, .. } => alt.len(),
            GEdit::Del { s, e } | GEdit::Dup { s, e } | GEdit::Inv { s, e } => span(*s, *e),
            // Unreachable (this measure runs on lowered edit sets only); its
            // stated anchor is the widest thing that can be said about an
            // unresolved tract without reading the reference.
            GEdit::Repeat { s, e, .. } => span(*s, *e),
            GEdit::Delins { s, e, alt } => span(*s, *e).max(alt.len()),
        })
        .sum()
}

/// The same measure for a derived piece set (see `changed_columns_of_edits`).
///
/// # Minimizing this is ferro policy, not compliance
///
/// The recommendations state their design values at `background/basics.md:38`
/// — "designed to be **stable**, **meaningful**, **memorable**, and
/// **unequivocal**". Minimality is **not** among them, and stability is first.
/// So a column count is a tie-break this project chose; it is never the answer
/// to "what does the spec require", and no comment here should imply otherwise.
///
/// The spec in fact prefers a non-minimal description in its own worked
/// example: `DNA/delins.md:44-47`'s spanning
/// `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` covers 66 columns where the split
/// alternative at `:46` needs 40, and `:47` recommends the spanning form
/// anyway. Whether that carve-out governs the general instruction at `:17` is
/// open — `rulings[delins-merge-vs-individual-gap-two-or-more]`.
fn changed_columns_of_pieces(pieces: &[Piece]) -> usize {
    pieces
        .iter()
        .map(|piece| (piece.ref_end - piece.ref_start).max(piece.alt.len()))
        .sum()
}

/// Whether the unchanged runs separating these pieces are long enough to be real
/// separation rather than a byproduct of the alignment search.
///
/// `best_alignment` scores every single-gap placement and keeps the highest, so
/// it actively hunts for the placement with the most matches. A lone unchanged
/// base inside a heavily rewritten block is very often that hunt's own artefact
/// rather than a base the variant left alone.
///
/// **This threshold is an implementer's calibration, not a spec mandate**, and it
/// is worth being exact about that. `delins.md:44-47` does rule that a
/// replacement whose payload merely "aligns" is better described as one `delins`,
/// and the harm it names — tools drawing wrong protein consequences — is exactly
/// what the conformance corpus shows, where splitting
/// `g.5207_5212delinsGTCCTGTGCTCATTATCTGGC` turns a single
/// `p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)` into three bogus ones. But that
/// worked example is a net *deletion* (52 nt -> 14 nt), so the spec does not
/// itself draw the line where this does.
///
/// One property is measured rather than assumed:
///
/// * **Measured before the 3'-shift.** Re-checking afterwards lets every one of
///   the corpus's coincidental splits back through, because a shift widens the
///   gap to a piece's left neighbour.
///
/// Alternatives that look more principled do not survive measurement. Gating on
/// `best_alignment`'s score margin cannot work at all: in these cases the winning
/// placement *is* the 5'-most one, so the margin is zero. Match-density variants
/// trade corpus rows against unit tests without a threshold that satisfies both.
///
/// # Net deletions are covered too, since #1271
///
/// They were not always. This rule applied only to net insertions, and a long
/// net *deletion* was held instead by a second, smaller block bound
/// (`MAX_UNGUARDED_SPLIT_BLOCK`, 32) that guarded by accident — it was simply
/// below the sizes at risk. The spec's own worked example is what that accident
/// was protecting: on `delins.md:44-47`'s `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`
/// (52 nt -> 14 nt) the aligner, left to itself, returns **four** pieces —
///
/// ```text
/// [0, 39) [40, 41) [42, 49) [50, 52)
/// ```
///
/// — which is exactly the harm that passage names, and which the corpus shows as
/// one correct `p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)` becoming three bogus
/// consequences.
///
/// Extending this rule to every length-changing regime replaces that accident
/// with a rule about the derived pieces: on that block `net_change` is 38, the
/// required separation is therefore [`RAISED_PIECE_SEPARATION`], and the widest
/// gap between consecutive pieces is one base — so the split is refused on its
/// merits and the block stays whole. `partition_block_refuses_a_coincidental_net_deletion_split`
/// pins both halves of that.
///
/// The extension was expected to be costly and measured not to be. #1271 records
/// that a naive attempt "breaks #1232 and #1157 … and breaks confluence
/// outright"; re-measured on this tree, all 13 tests of those two issues pass
/// (including `sequence_identical_delins_and_allele_normalize_equal`, the
/// confluence one), all 168 manifest-backed conformance tests pass with the
/// FAIL-set ledgers unmoved, and the spec enumeration *improves* by one row —
/// `projection-splits-single-member` 10 -> 9 divergences, `projection-pinned`
/// 1167 -> 1168 passing. Whatever made the naive attempt fail was fixed
/// elsewhere in the meantime, plausibly by #1237's regime-aware bound and the
/// sequence-first splitter that followed it.
///
/// # Measuring separation over a corpus: get the insertion arithmetic right
///
/// This function works on `Piece`s, where a pure insertion is already a
/// zero-width span. A corpus census works on *HGVS members*, where it is not,
/// and the difference has produced two published-then-retracted numbers.
///
/// Model each member as a closed interval and an insertion `A_B` as the
/// **empty** interval `lo = A + 1, hi = A` — an insertion consumes no reference
/// base, it names the junction it sits in. Then
///
/// ```text
/// separation = next.lo - prev.hi - 1
/// ```
///
/// Treating `A_B` as a consumed two-base span instead reported "121 pairs at
/// separation 0 — a spec violation", where the true count over 5.76M rows is
/// **zero**. Correcting it also moved the whole gap distribution materially
/// (gap 2 had been quoted at 42.7%, gap 3 at 14.8%), so any figure taken from a
/// pre-correction distribution must be **re-derived rather than adjusted**.
fn separations_are_meaningful(pieces: &[Piece], net_change: usize) -> bool {
    // `ref_end` is exclusive and a pure insertion has `ref_start == ref_end`, so
    // this already counts unchanged *bases* rather than event indices: a
    // junction contributes no width because it occupies none. The sequence-first
    // splitter had to be corrected for exactly this, because it works in event
    // space where a junction and a changed position are both one offset.
    let required = if net_change > MAX_SINGLE_BASE_SEPARATION_CHANGE {
        RAISED_PIECE_SEPARATION
    } else {
        MIN_PIECE_SEPARATION
    };
    pieces
        .windows(2)
        .all(|pair| pair[1].ref_start.saturating_sub(pair[0].ref_end) >= required)
}

/// Whether a proposed split buys nothing over the spanning `delins` it came
/// from — that is, whether **every** member would render as a `delins` too.
///
/// # This is ferro policy, not a spec rule
///
/// It is tempting to cite `general.md:56`'s prioritisation here, and an earlier
/// revision of this comment did — claiming `delins` was "last" in it. That
/// citation does not reach. The list is "(1) substitution, (2) deletion, (3)
/// inversion, (4) duplication, (5) insertion", and `delins` is not in it at
/// all. What it ranks is the type of *one description of one change*, which is
/// how msto's #1231/#1233 use it — comparing two two-member descriptions
/// member-wise. It never ranks a multi-member allele against a single spanning
/// description.
///
/// What licenses the collapse is an analogy plus the corpus. `delins.md:44-47`
/// keeps the spanning form where the payload merely *"aligns"* against a
/// coincidentally matched interior base, and that is the shape of every corpus
/// row required to stay whole. A split whose members are all `delins` buys no
/// descriptive gain over the spanning `delins` it came from, so in the one
/// regime where the matched base is least believable — a single unchanged base,
/// see [`every_separation_is_a_single_base`] — this prefers the form the
/// analogous spec case prefers. Stated as a policy choice inside a spec gap so
/// it is not mistaken for a mandate.
///
/// This is what separates #422 from #999's negative control. The two are the
/// same shape at [`MIN_PIECE_SEPARATION`] — a net insertion with one
/// coincidentally matched interior base — and their split *member types* are the
/// only difference: `delins` + `delins` for #422, `ins` + `sub` for #999.
/// `issue_422_cross_reference_ins.rs` calls that conjunction impossible; it is
/// not, and #1235's own comment asks for precisely this fix.
///
/// # Typing is block-local
///
/// `ins`, `del` and `sub` are structural. So are `dup` and repeat notations —
/// they add or remove copies, so they arrive as insertions or deletions and are
/// already excluded. `inv` is the one renderer type that hides inside a
/// structural `delins`, and it is detectable from the member's own span:
/// [`is_inversion`] does it, on [`reverse_complement_bytes`], which refuses a
/// byte it cannot complement rather than passing it through. Nothing here needs
/// the renderer, so the test runs on the trimmed block before rendering.
///
/// Measured over an exhaustive reverse-complement-closed sweep, omitting the
/// `inv` and self-repeat cases would collapse 694 of 1 094 firings (`AT`, ≤7 nt)
/// and 2 592 of 26 244 (`ACGT`, ≤5 nt) — so both are load-bearing, not defensive.
fn split_buys_no_higher_priority_type(pieces: &[Piece], reference: &[u8]) -> bool {
    pieces.iter().all(|piece| {
        let span = piece.ref_end - piece.ref_start;
        span > 0
            && !piece.alt.is_empty()
            && !(span == 1 && piece.alt.len() == 1)
            && !is_inversion(piece, reference)
            && !is_tandem_duplication(piece, reference)
            && !alt_repeats_its_own_span(piece, reference)
    })
}

/// A member whose payload is its own reference span repeated, which renders as
/// a `dup` (two copies) or a repeat notation rather than a `delins`.
///
/// [`is_tandem_duplication`] covers the *insertion* spelling of a duplication,
/// where the copy sits outside the member's span; this covers the `delins`
/// spelling, where the member's span carries the copies itself. Both are
/// block-local. Case-folded, because reference FASTAs arrive soft-masked.
fn alt_repeats_its_own_span(piece: &Piece, reference: &[u8]) -> bool {
    let span = &reference[piece.ref_start..piece.ref_end];
    !span.is_empty()
        && piece.alt.len() > span.len()
        && piece.alt.len().is_multiple_of(span.len())
        && piece
            .alt
            .chunks(span.len())
            .all(|chunk| chunk.eq_ignore_ascii_case(span))
}

/// The alignment maximizing matched bases, ties broken by the 5'-most gap.
///
/// The minimal edit set is **not** unique when the blocks differ in length —
/// #1232 has four equally minimal, equally stable encodings of one variant — so
/// a deterministic tie-break is required, not optional. Maximizing matches is
/// what the separation rule (`delins.md:17`) asks for; the 5'-most tie-break is
/// arbitrary, since the spec does not reach it, but each piece is 3'-shifted
/// afterwards so the 3' rule is still honoured per member.
///
/// Only single-gap alignments are considered: one contiguous indel plus
/// substitutions. Allowing compensating gaps (a deletion *and* an insertion
/// inside one block) lets the aligner manufacture matches from coincidence and
/// shred a genuinely contiguous replacement — the regime `delins.md:44-47`
/// warns about, seen concretely in the #1034/#1040/#182 "stays delins" cases.
///
/// # Cost
///
/// Scoring is arithmetic rather than materialized. A single-gap alignment is
/// fully described by its gap offset `k`, so a placement's matched-column count
/// can be read straight off the two slices and only the winner is ever built.
/// Building every candidate instead cost one `Vec<Column>` per placement — at
/// [`MAX_SPLIT_BLOCK`] that is 1025 allocations of 1024 columns each, roughly
/// 33 MB of churn for one call, all but one dropped immediately.
///
/// The score is also *separable*: the columns before the gap pair `i` with `i`
/// and those after it pair `i` with `i + gap_len`, and neither set depends on
/// `k` beyond where it is cut. Two running sums therefore answer every placement
/// in O(1), making the search linear rather than quadratic — at
/// [`MAX_SPLIT_BLOCK`] that is ~1024 comparisons against ~1e6.
fn best_alignment(reference: &[u8], result: &[u8]) -> Option<Vec<Column>> {
    if reference.len() == result.len() {
        return Some((0..reference.len()).map(|i| (Some(i), Some(i))).collect());
    }
    let (shorter_is_ref, gap_len) = if reference.len() < result.len() {
        (true, result.len() - reference.len())
    } else {
        (false, reference.len() - result.len())
    };
    let anchored = if shorter_is_ref {
        reference.len()
    } else {
        result.len()
    };

    // Matched columns for the placement whose gap follows `k` aligned positions.
    // Columns before the gap pair `i` with `i`; those after it pair across the
    // gap. Exactly what `single_gap_alignment` lays out, counted instead of
    // built.
    //
    // `leading(k)` grows with `k` and `trailing(k)` shrinks with it, so both are
    // swept once: `leading` accumulates the column `k` just crossed, and
    // `trailing` starts at its `k == 0` value and sheds that same column.
    let ungapped = |i: usize| -> bool {
        if shorter_is_ref {
            reference[i] == result[i + gap_len]
        } else {
            reference[i + gap_len] == result[i]
        }
    };
    let mut leading = 0usize;
    let mut trailing = (0..anchored).filter(|&i| ungapped(i)).count();

    // How many members a placement separates the change into. Computed only to
    // break a tie, because it costs a pass over the block where the score costs
    // O(1) — see the tie-break note below.
    let members_at = |k: usize| -> usize {
        let columns =
            single_gap_alignment(k, gap_len, shorter_is_ref, reference.len(), result.len());
        pieces_from_columns(&columns, reference, result).len()
    };

    // `(matched columns, members separated, k)`.
    let mut best: Option<(usize, usize, usize)> = None;
    for k in 0..=anchored {
        if k > 0 {
            leading += usize::from(reference[k - 1] == result[k - 1]);
            trailing -= usize::from(ungapped(k - 1));
        }
        let matches = leading + trailing;
        match best {
            Some((score, _, _)) if matches < score => continue,
            // A TIE, and the score cannot choose. In a repeat tract every
            // placement matches the same number of columns, so this is the
            // common case there rather than a corner: `AAA -> CA` ties three
            // ways. `delins.md:17` decides it — "two variants separated by one
            // or more nucleotides should be described individually and not as a
            // 'delins'" — so prefer the placement that leaves a matched base
            // *between* the changes, which is the one separating them into more
            // members. `general.md:56` is *not* a second authority for this:
            // `delins` does not appear in its list at all, and what it ranks is
            // the type of one description of one change — see
            // `split_buys_no_higher_priority_type`'s doc.
            //
            // This is not the 5'-most/3'-most flip that was tried and rejected:
            // it does not rank placements by position at all. `AAA -> CA`'s
            // winner happens to be 3'-most, but the criterion is what the
            // placement separates, and `split_buys_no_higher_priority_type`
            // still collapses a split that buys nothing.
            Some((score, separated, _)) if matches == score => {
                if anchored <= MAX_TIE_BREAK_SWEEP {
                    let candidate = members_at(k);
                    if candidate > separated {
                        best = Some((matches, candidate, k));
                    }
                }
            }
            _ => {
                let separated = if anchored <= MAX_TIE_BREAK_SWEEP {
                    members_at(k)
                } else {
                    0
                };
                best = Some((matches, separated, k));
            }
        }
    }
    // A single gap cannot express *insertion, retained reference, insertion*,
    // so consider the two-gap form when the best single gap leaves a reference
    // base unmatched. Scoped to net insertions by `shorter_is_ref`, which is
    // what keeps it clear of the corpus the single-gap restriction protects —
    // see [`two_gap_insertion_alignment`].
    if shorter_is_ref && best.is_some_and(|(matches, _, _)| matches < reference.len()) {
        if let Some(columns) = two_gap_insertion_alignment(reference, result) {
            return Some(columns);
        }
    }

    // Only the winner is materialized.
    best.map(|(_, _, k)| {
        single_gap_alignment(k, gap_len, shorter_is_ref, reference.len(), result.len())
    })
}

/// The alignment that keeps the **whole** reference intact across two
/// insertions, or `None` when the pair does not have that shape (#1260).
///
/// PR #1285 named the gap this fills: for #1260 the answer is "a two-gap
/// alignment — insertion, retained base, insertion — which `best_alignment`'s
/// single-gap restriction cannot express". One contiguous indel must either
/// swallow the retained base into a spanning `delins` or leave one insertion
/// unaccounted for.
///
/// # The shape, and why it is not just `I1 + reference + I2`
///
/// The general all-matched two-gap alignment is
///
/// ```text
/// result = ref[..a] + I1 + ref[a..b] + I2 + ref[b..]
/// ```
///
/// with `0 <= a < b <= ref.len()` and both insertions non-empty. An earlier
/// revision searched only `a == 0, b == ref.len()`, deriving that from the
/// caller having trimmed common flanks: on a trimmed block `a > 0` would mean
/// the two blocks share a first base and `b < ref.len()` a last one, both
/// contradicting the trim.
///
/// That derivation is sound *only while the precondition holds*, and it makes
/// the function silently wrong the moment a caller widens a block — which is
/// exactly what any attempt to give back an ambiguously-trimmed flank base
/// must do. #1260 is the worked example: trimmed its block is `A -> CAC` with
/// `a = 0`, and with one flank base restored per side it is `AAA -> ACACA`,
/// whose only decomposition is `a = 1, b = 2`. A search fixed at `a == 0`
/// finds nothing there, `best_alignment` falls through to the single-gap
/// winner, and #1260 reverts to a spanning `delins`.
///
/// No caller widens a block today, so this generalisation changes no current
/// result — it removes a hidden coupling between this function and the
/// caller's trim, and it is verified against blocks that exercise `a > 0` and
/// `b < ref.len()` directly rather than through a caller.
///
/// `b > a` is required rather than incidental: with `b == a` the two insertions
/// are adjacent in `result`, which is a single gap the existing search already
/// places.
///
/// # Cost
///
/// Every reference base survives verbatim in this form, so `a` cannot exceed the
/// two blocks' common prefix and `ref.len() - b` cannot exceed their common
/// suffix. Both are 0 on a trimmed block — collapsing the search to the single
/// case the earlier revision hardcoded — and at most one after a restored flank.
/// The sweep is therefore bounded by the *residual* affix lengths, not by the
/// block length. A block whose affixes would make that sweep large is declined
/// outright — see [`MAX_TWO_GAP_PLACEMENTS`] — rather than searched under a
/// narrower bound than this doc describes.
///
/// # The tie-break is an implementer's choice
///
/// The decomposition is **not unique**: `AA -> AAAA` admits `(a, b)` of `(0, 1)`,
/// `(0, 2)` and `(1, 2)`, all with `I1 = I2 = "A"`, all matching every reference
/// base and all yielding two members. Nothing in the spec reaches this. The
/// smallest `a`, then the smallest `|I1|`, is taken — the 5'-most placement,
/// preserving the earlier revision's behaviour. Deliberately *not* the 3'-most
/// rule [`best_alignment`] applies to a member-count tie: adopting it here would
/// move every such tie in a homopolymer.
///
/// # Why this does not reopen the corpus the single-gap rule protects
///
/// `merge.rs`'s single-gap restriction exists to stop *compensating* gaps — a
/// deletion **and** an insertion in one block — from manufacturing matches out
/// of coincidence and shredding a contiguous replacement (#1034/#1040/#182).
/// Both gaps here are insertions in a net-insertion block, so no deletion is
/// introduced and no match is manufactured: every matched column is a reference
/// base surviving verbatim. That argument never depended on `a` being 0 — all
/// three chunks are verbatim reference either way. Cluster A is structurally out
/// of reach — #1040 and the `inv` cases are equal-length, #1157A is a net
/// deletion, and an equal-length or net-deletion block never reaches this
/// function.
///
/// Comparison is byte-exact, matching [`best_alignment`]'s own `==` throughout.
/// Making only this path case-insensitive would let the two-gap form win on a
/// soft-masked block where the single-gap search sees a mismatch.
fn two_gap_insertion_alignment(reference: &[u8], result: &[u8]) -> Option<Vec<Column>> {
    // `reference.len() + 2` is the shortest result that leaves room for two
    // insertions of one base each; an empty reference has nothing to retain.
    if reference.is_empty() || result.len() < reference.len() + 2 {
        return None;
    }
    let ref_len = reference.len();
    let inserted = result.len() - ref_len;

    // Bounds, not shortcuts: a matched `ref[..a]` is a common prefix of the two
    // blocks, and a matched `ref[b..]` a common suffix. Not clamped — narrowing
    // them changes the answer rather than bounding the cost, which is what
    // [`MAX_TWO_GAP_PLACEMENTS`] is for.
    let common_prefix = reference
        .iter()
        .zip(result)
        .take_while(|(r, a)| r == a)
        .count();
    let common_suffix = reference
        .iter()
        .rev()
        .zip(result.iter().rev())
        .take_while(|(r, a)| r == a)
        .count();

    // Refuse the whole search rather than truncate it, so a block over budget
    // gets `best_alignment`'s single-gap answer instead of one found under a
    // narrower sweep than this function documents.
    let lowest_b = (ref_len - common_suffix.min(ref_len)).max(1);
    let placements = (common_prefix + 1)
        .saturating_mul(inserted.saturating_sub(1))
        .saturating_mul(ref_len + 1 - lowest_b.min(ref_len + 1));
    if placements > MAX_TWO_GAP_PLACEMENTS {
        return None;
    }

    // 5'-most: smallest `a`, then smallest `|I1|`, then smallest `b`.
    for a in 0..=common_prefix {
        for first_insert in 1..inserted {
            for b in (a + 1).max(lowest_b)..=ref_len {
                let retained_start = a + first_insert;
                let retained_end = retained_start + (b - a);
                if result[retained_start..retained_end] != reference[a..b] {
                    continue;
                }
                // `result[..a] == reference[..a]` already holds by `a <=
                // common_prefix`; the suffix is checked rather than inferred so
                // the three chunk equalities read together.
                let tail_start = retained_end + (inserted - first_insert);
                if result[tail_start..] != reference[b..] {
                    continue;
                }
                let mut columns = Vec::with_capacity(result.len());
                columns.extend((0..a).map(|i| (Some(i), Some(i))));
                columns.extend((a..retained_start).map(|j| (None, Some(j))));
                columns.extend((a..b).map(|i| (Some(i), Some(retained_start + i - a))));
                columns.extend((retained_end..tail_start).map(|j| (None, Some(j))));
                columns.extend((b..ref_len).map(|i| (Some(i), Some(tail_start + i - b))));
                return Some(columns);
            }
        }
    }
    None
}

/// Columns for an alignment with a single gap of `gap_len` after `k` aligned
/// positions. The gap is on the reference side when `shorter_is_ref` (a net
/// insertion) and on the result side otherwise (a net deletion).
fn single_gap_alignment(
    k: usize,
    gap_len: usize,
    shorter_is_ref: bool,
    ref_len: usize,
    alt_len: usize,
) -> Vec<Column> {
    let mut columns = Vec::with_capacity(ref_len.max(alt_len));
    for i in 0..k {
        columns.push((Some(i), Some(i)));
    }
    if shorter_is_ref {
        for j in 0..gap_len {
            columns.push((None, Some(k + j)));
        }
        for i in k..ref_len {
            columns.push((Some(i), Some(i + gap_len)));
        }
    } else {
        for j in 0..gap_len {
            columns.push((Some(k + j), None));
        }
        for i in k..alt_len {
            columns.push((Some(i + gap_len), Some(i)));
        }
    }
    columns
}

/// Group consecutive changed columns into pieces.
fn pieces_from_columns(columns: &[Column], reference: &[u8], result: &[u8]) -> Vec<Piece> {
    let changed = |(r, a): &Column| match (r, a) {
        (Some(ri), Some(ai)) => reference[*ri] != result[*ai],
        _ => true,
    };
    let mut pieces: Vec<Piece> = Vec::new();
    let mut run: Vec<Column> = Vec::new();
    let mut next_ref = 0usize;

    let flush = |run: &mut Vec<Column>, boundary: usize, pieces: &mut Vec<Piece>| {
        if run.is_empty() {
            return;
        }
        let ref_indices: Vec<usize> = run.iter().filter_map(|(r, _)| *r).collect();
        let (ref_start, ref_end) = match (ref_indices.first(), ref_indices.last()) {
            (Some(first), Some(last)) => (*first, *last + 1),
            // Pure insertion: a zero-width span at the following reference base.
            _ => (boundary, boundary),
        };
        let alt = run
            .iter()
            .filter_map(|(_, a)| *a)
            .map(|ai| result[ai])
            .collect();
        pieces.push(Piece {
            ref_start,
            ref_end,
            alt,
        });
        run.clear();
    };

    for column in columns {
        if let Some(ri) = column.0 {
            next_ref = ri;
        }
        if changed(column) {
            run.push(*column);
        } else {
            flush(&mut run, next_ref, &mut pieces);
        }
    }
    flush(&mut run, reference.len(), &mut pieces);
    pieces
}

/// Shift each length-changing piece in `direction`, clamped so it cannot reach
/// a sibling.
///
/// The clamp is the #1234 fix. The 3' rule (`general.md:41`) is stated per
/// description with no allele-awareness, and the only text that ever addressed
/// cross-member shifting was the *rejected* SVD-WG010. Left unbounded, a
/// deletion shifts over a downstream substitution and the two members overlap —
/// malformed, and denoting a different sequence. Substitutions are anchored and
/// never shift.
///
/// **`direction` is a parameter, not a constant (#1429).** This was hardcoded to
/// `ThreePrime`, and `canonicalize_from_sequence` never received the caller's
/// direction at all — so under a 5' shuffle the derivation handed back a
/// *3'-shifted* allele, which the per-member pipeline then moved on the next
/// pass:
///
/// ```text
/// c.[16dup;18_19insC]  ->  c.*1_*2insCT  ->  c.18_*1insTC
/// ```
///
/// Two passes, two directions. The 5' answer was also simply wrong on the first
/// pass, independently of the idempotency: a caller asking for a 5' shuffle got
/// a 3'-shifted merged allele.
fn shift_pieces(pieces: &mut [Piece], ref_bytes: &[u8], direction: ShuffleDirection) {
    for i in 0..pieces.len() {
        if !pieces[i].is_pure_indel() {
            continue;
        }
        let left = if i == 0 { 0 } else { pieces[i - 1].ref_end };
        let right = pieces
            .get(i + 1)
            .map_or(ref_bytes.len(), |next| next.ref_start);
        let bounds = Boundaries::new(left as u64, right as u64);
        let shuffled = shuffle(
            ref_bytes,
            &pieces[i].alt,
            pieces[i].ref_start as u64,
            pieces[i].ref_end as u64,
            &bounds,
            direction,
        );
        let old_start = pieces[i].ref_start;
        let new_start = shuffled.start as usize;
        pieces[i].ref_start = new_start;
        pieces[i].ref_end = shuffled.end as usize;

        // `shuffle` returns coordinates only. Moving a *pure insertion*'s point
        // 3' by `k` rotates the inserted sequence left by `k mod len`: inserting
        // `TG` before position 16 is not the variant that inserts `TG` before
        // 17. Leaving `alt` alone emits a description denoting a *different*
        // sequence, and because the wrong form is a stable fixed point no
        // idempotency check can see it — only the pass's own round-trip assert
        // caught it, on 56 conformance rows.
        //
        // `k mod len` is not a guess: `shuffle`'s insertion predicate compares
        // `ref[start + k]` against `alt[(new_end - start) % alt.len()]`, so
        // `shuffle` already *defines* the post-shift payload as
        // `rotate_left(k mod len)`. The same modulus emits exactly the variant
        // `shuffle` proved legal. A pure deletion has an empty `alt`.
        //
        // This is the hazard the per-member pipeline documents inside
        // `Normalizer::normalize_na_edit` (the #1157 follow-up).
        if !pieces[i].alt.is_empty() {
            // Each direction moves the point only its own way — `shuffle`
            // initialises its cursor to `start` and mutates it solely via
            // `+= 1` or `-= 1` — so the two are asserted separately rather
            // than with one `abs`, which would mask a `shuffle` that moved a
            // piece the wrong way.
            match direction {
                ShuffleDirection::ThreePrime => debug_assert!(
                    new_start >= old_start,
                    "a 3'-shuffle cannot move a piece 5' (old={old_start}, new={new_start})",
                ),
                ShuffleDirection::FivePrime => debug_assert!(
                    new_start <= old_start,
                    "a 5'-shuffle cannot move a piece 3' (old={old_start}, new={new_start})",
                ),
            }
            let distance = new_start.abs_diff(old_start);
            let rotation = distance % pieces[i].alt.len();
            if rotation > 0 {
                // Mirror images. Moving the point 3' by `k` rotates the payload
                // left by `k`; moving it 5' by `k` rotates it right by the same
                // `k`, because the payload's phase against the reference is what
                // the rotation tracks and the point moved the other way.
                match direction {
                    ShuffleDirection::ThreePrime => pieces[i].alt.rotate_left(rotation),
                    ShuffleDirection::FivePrime => pieces[i].alt.rotate_right(rotation),
                }
            }
        }
    }
}

/// Re-merge pieces the coding one-amino-acid exception keeps together.
///
/// The separation rule says changes separated by one or more unchanged
/// nucleotides are described individually — but `general.md:35` carves out the
/// case where they are separated by *exactly one* nucleotide and together affect
/// **one amino acid**, which must be a delins:
///
/// > **exception**: two variants separated by one nucleotide, together affecting
/// > one amino acid, should be described as a "delins".
///
/// `NM_004006.2:c.145_147delinsTGG` is the spec's worked example, and it
/// explicitly rejects `c.[145C>T;147C>G]` even though position 146 is unchanged.
/// This is still current law: SVD-WG010, which would have replaced it with a
/// context-free "separated by less than two nucleotides" rule at every level,
/// was **rejected** (`SVD-WG010.md:5-8`).
///
/// The exception needs reading-frame context, so it applies only where one
/// exists — never on `g.`/`m.`/`n.`, and on `r.` only when the transcript
/// actually has a CDS (see [`AxisFrame::reading_frame`]). Keying it off the
/// axis alone was #1241: `CisKind::Rna` says nothing about whether the
/// transcript is coding, so a non-coding `r.` got a codon frame stamped onto it
/// and stopped converging with the `n.` spelling of the very same change.
///
/// `merge.rs`'s codon-frame merge already implements the exception for separate
/// input members; without this the re-derivation would split such a delins
/// straight back apart, leaving the two spellings converging in opposite
/// directions.
fn apply_coding_codon_exception(
    pieces: &mut Vec<Piece>,
    reading_frame: bool,
    w_lo: i64,
    ref_bytes: &[u8],
) {
    if !reading_frame {
        return;
    }
    // A single changed reference base replaced by a single alt base — the
    // `Substitution` the spec's "two variants" are.
    let is_substitution =
        |piece: &Piece| piece.ref_end == piece.ref_start + 1 && piece.alt.len() == 1;

    let mut index = 1;
    while index < pieces.len() {
        let previous_end = pieces[index - 1].ref_end;
        let next_start = pieces[index].ref_start;
        // Exactly one unchanged reference base between the two pieces.
        let separated_by_one = next_start == previous_end + 1;
        // The exception is a *triplet*, not a span: `build_split_variants`
        // matches exactly `[Sub@p; Identity@p+1; Sub@p+2]` with `p` and `p+2`
        // in one codon, and `general.md:35` says "two variants separated by one
        // nucleotide", not "any two edits whose hull fits a codon". So the
        // left piece must be a lone substitution, the right piece must *begin*
        // with one, and only `p`/`p+2` are codon-tested.
        //
        // Testing the hull instead is what made this rule disagree with the
        // per-member pipeline. On `c.10_13delinsTCAG` (ref `CCCC`, codon 4 =
        // 10-12) the changed columns are 10, 12 and 13, so the pieces arrive as
        // `[10]` and `[12,13]`; the hull 10..13 crosses into codon 5, the merge
        // declined, and the result was `[10C>T;12_13delinsAG]` against the
        // pipeline's `[10_12delinsTCA;13C>G]`.
        let starts_with_substitution =
            pieces[index].ref_end - pieces[index].ref_start == pieces[index].alt.len();
        let triplet_start = w_lo + pieces[index - 1].ref_start as i64;
        if separated_by_one
            && is_substitution(&pieces[index - 1])
            && starts_with_substitution
            && same_codon(triplet_start, triplet_start + 2)
        {
            let middle = ref_bytes[previous_end];
            // Take only the triplet's third base. Anything after it belongs to
            // the next codon and stays a member of its own — that split is the
            // whole point, and the remainder is re-examined on the next
            // iteration so a second triplet can still form from it.
            let next = &mut pieces[index];
            let tail_alt = next.alt.split_off(1);
            let third = next.alt[0];
            next.ref_start += 1;
            next.alt = tail_alt;
            let exhausted = next.ref_start == next.ref_end;
            if exhausted {
                pieces.remove(index);
            }
            let previous = &mut pieces[index - 1];
            // The unchanged base is now interior to the merged replacement, so
            // it has to be carried through explicitly.
            previous.alt.push(middle);
            previous.alt.push(third);
            previous.ref_end = previous.ref_start + 3;
        }
        index += 1;
    }
}

/// Split a sequence-first merge back apart when it crosses a codon boundary.
///
/// `min_separation` (`seqfirst::MIN_SEPARATION`, `2`, on a reading-frame axis)
/// captures only the *distance* half of `general.md:35`'s exception — merging
/// any two runs separated by fewer than 2 unchanged bases, whichever codons
/// they fall in. The exception itself is narrower: "two variants separated by
/// one nucleotide, **together affecting one amino acid**" — both conditions,
/// not one. A pair in different codons affects *two* amino acids and
/// `general.md:34`'s plain rule governs it instead: described individually.
///
/// Scoped to the exact shape the exception describes — two single-nucleotide
/// substitutions with their shared unchanged base carried in the middle
/// (`ref_end - ref_start == 3`, `alt.len() == 3`, the middle alt base equal to
/// the middle reference base, both flanking bases actually changed) — using
/// the same [`same_codon`] arithmetic [`apply_coding_codon_exception`] already
/// uses for the live splitter's own triplet, on the same axis distinction
/// (`reading_frame`, from [`AxisFrame`]). A wider merged piece is not "two
/// variants" in `general.md:35`'s sense at all — `LRG_199t1:c.850_901`
/// (`delins.md:44`) merges a 52-base member for an unrelated reason (its
/// `partition_block` sibling reaches the aligner and is then refused by
/// [`separations_are_meaningful`], whose gaps there are all one base wide —
/// before #1271 it was refused earlier still, by a block bound), not this
/// exception, so it must not be touched here and is not: 52 ≠ 3.
///
/// Pinned by `sequence_first_split_declines_the_codon_exception_across_a_boundary`
/// (`GCA -> TCC`, `general.md:34`'s split rule wins) and
/// `sequence_first_split_keeps_the_codon_exception_within_one_codon`
/// (`CCC -> GCA`, `general.md:35`'s exception wins) — the same two blocks
/// `merge_consecutive_edits_tests::test_no_codon_frame_pair_straddles_codon_boundary`
/// and `test_codon_frame_two_subs_one_codon` pin for the live splitter, so
/// both sides are checked against the identical spec-derived expectation.
fn split_codon_incompatible_triplets(
    pieces: &mut Vec<Piece>,
    reading_frame: bool,
    w_lo: i64,
    ref_bytes: &[u8],
) {
    if !reading_frame {
        return;
    }
    let mut index = 0;
    while index < pieces.len() {
        let piece = &pieces[index];
        // `ref_end <= ref_bytes.len()` is tested before any indexing. Every
        // caller today passes pieces carved from this same window, so the span is
        // in range — but this reads three bytes by raw index off a piece it does
        // not own, and a piece whose span outran the window would panic rather
        // than decline. An out-of-range span is simply not a merged triplet.
        let is_merged_triplet = piece.ref_end == piece.ref_start + 3
            && piece.ref_end <= ref_bytes.len()
            && piece.alt.len() == 3
            && piece.alt[1] == ref_bytes[piece.ref_start + 1]
            && piece.alt[0] != ref_bytes[piece.ref_start]
            && piece.alt[2] != ref_bytes[piece.ref_start + 2];
        if !is_merged_triplet {
            index += 1;
            continue;
        }
        let triplet_start = w_lo + piece.ref_start as i64;
        if same_codon(triplet_start, triplet_start + 2) {
            index += 1;
            continue;
        }
        let ref_start = piece.ref_start;
        let first_alt = piece.alt[0];
        let third_alt = piece.alt[2];
        pieces.splice(
            index..=index,
            [
                Piece {
                    ref_start,
                    ref_end: ref_start + 1,
                    alt: vec![first_alt],
                },
                Piece {
                    ref_start: ref_start + 2,
                    ref_end: ref_start + 3,
                    alt: vec![third_alt],
                },
            ],
        );
        index += 2;
    }
}

/// Merge pieces the 3'-shift left touching.
///
/// Two or more consecutive changed nucleotides are one delins
/// (`delins.md:16`), so once a shift closes the gap the pieces must combine —
/// this is what turns #1234's clamped `5_7del` plus `8A>T` into `5_8delinsT`.
fn coalesce_adjacent_pieces(pieces: &mut Vec<Piece>) {
    let mut i = 1;
    while i < pieces.len() {
        debug_assert!(
            pieces[i - 1].ref_end <= pieces[i].ref_start,
            "strict overlap between adjacent pieces: the 3'-shuffle cannot produce this",
        );
        if pieces[i - 1].ref_end == pieces[i].ref_start {
            let next = pieces.remove(i);
            let prev = &mut pieces[i - 1];
            prev.ref_end = next.ref_end;
            prev.alt.extend(next.alt);
        } else {
            i += 1;
        }
    }
}

/// Shrink each piece to the hull of its own differences.
///
/// A piece is a *region*, and a region is not obliged to be the narrowest span
/// that spells the change it holds. Where an edit sits inside an ambiguous run
/// the two come apart: on `ACAA -> CAACAG` the inserted `CA` can sit at more
/// than one offset, so a partition derived from the alignment steps common to
/// every minimal alignment confines it to `reference[0..1]` — which spells a
/// `delins` of `CAA` over one base for what those same bases denote as an
/// insertion of `CA`.
///
/// This is [`trim_common_flanks`] applied *within* a piece rather than to the
/// whole block, and like that call it reads only the reference span and the
/// payload — never the input spelling — so it cannot reintroduce a dependence
/// on how the variant was written.
///
/// Only a `delins` can move. A deletion has no payload to share a flank with,
/// an insertion spans no reference base, and a substitution's single base
/// differs by construction, so each is returned untouched; the widened forms
/// this exists to narrow are the only ones narrowed. Growing a run back to a
/// `dup` or a repeat spelling is a rendering decision made downstream, on the
/// narrowed member.
fn shrink_pieces_to_differences(pieces: &mut [Piece], ref_bytes: &[u8]) {
    for piece in pieces.iter_mut() {
        let (lo, hi_ref, hi_alt) =
            trim_common_flanks(&ref_bytes[piece.ref_start..piece.ref_end], &piece.alt);
        // A piece that trims away entirely spells no change; leaving it whole
        // keeps `anchor_for_piece` from building a member with no extent.
        if lo == hi_ref && lo == hi_alt {
            continue;
        }
        piece.alt.truncate(hi_alt);
        piece.alt.drain(..lo);
        piece.ref_end = piece.ref_start + hi_ref;
        piece.ref_start += lo;
    }
}

/// Turn pieces back into HGVS members on the group's own axis.
///
/// Dispatch is [`build_merged`]'s, not a second copy of the kind→builder table:
/// `template` is `variants[0]`, which [`cis_kind_of`] already accepted, so its
/// `unreachable!` arm cannot fire here.
fn rebuild_members(
    pieces: &[Piece],
    template: &HgvsVariant,
    body: Region,
    w_lo: i64,
    ref_bytes: &[u8],
    ends: SequenceEnds,
    cds_end_axis: Option<i64>,
) -> Option<Vec<HgvsVariant>> {
    let mut members = Vec::with_capacity(pieces.len());
    // Window offset one past the last reference base the *previous* piece
    // claimed, for `duplication_anchor`'s disjointness check. `0` before the
    // first piece. Correct because `pieces` is in ascending offset order —
    // `partition_block` builds it that way and neither `shift_pieces`
    // nor `coalesce_adjacent_pieces` reorders — which the assertion below pins
    // rather than assumes, since a silent reordering would turn the check into a
    // no-op instead of a failure.
    let mut previous_ref_end = 0usize;
    for piece in pieces {
        debug_assert!(
            piece.ref_start >= previous_ref_end,
            "pieces must be in ascending offset order for the dup disjointness check"
        );
        let anchor = anchor_for_piece(
            piece,
            body,
            w_lo,
            ref_bytes,
            ends,
            previous_ref_end,
            cds_end_axis,
        )?;
        previous_ref_end = piece.ref_end;
        members.push(build_merged(template, anchor));
    }
    Some(members)
}

/// Build the [`Anchor`] for one piece, typing it under the spec's rules.
///
/// [`build_naedit`] already picks insertion / deletion / delins from the span
/// and replacement, which covers the mandates that 1 nt → 1 nt is a
/// substitution (`delins.md:15`) and that consecutive changes are one delins
/// (`delins.md:16`). A pure insertion needs one rule that span alone cannot
/// express, applied here before falling back to that typing: a payload resting
/// outside the sequence's own first or last base has no second anchor to name,
/// and is re-spelled as a single-position `delins` ([`boundary_delins_anchor`],
/// #1205/#1217).
///
/// One further shape cannot be read off the two lengths at all, because it
/// leaves them unchanged, so it is recognised from the reference before that
/// fallback: a tandem duplication ([`duplication_anchor`], `duplication.md:18`).
///
/// An inversion is the other such shape and is deliberately *not* typed here —
/// see [`AnchorForm`]. Both were once refused wholesale, sending the group back
/// to the per-member pipeline; the `dup` is now built, and the `inv` is left to
/// [`crate::normalize::rules`]'s single-span typing, which can still see it in
/// the rendered member.
///
/// `previous_ref_end` is the window offset one past the last reference base the
/// preceding piece claimed, which only [`duplication_anchor`] needs — see there.
fn anchor_for_piece(
    piece: &Piece,
    body: Region,
    w_lo: i64,
    ref_bytes: &[u8],
    ends: SequenceEnds,
    previous_ref_end: usize,
    cds_end_axis: Option<i64>,
) -> Option<Anchor> {
    let alt: Vec<Base> = piece
        .alt
        .iter()
        .filter_map(|b| Base::from_char(*b as char))
        .collect();
    if alt.len() != piece.alt.len() {
        return None; // non-IUPAC byte; refuse rather than mangle.
    }
    if piece.ref_start == piece.ref_end && !alt.is_empty() {
        if let Some(anchor) = duplication_anchor(piece, body, w_lo, ref_bytes, previous_ref_end) {
            return Some(anchor);
        }
        if let Some(anchor) = boundary_delins_anchor(piece, &alt, body, w_lo, ref_bytes, ends) {
            return Some(anchor);
        }
        if let Some(anchor) =
            cds_end_delins_anchor(piece, &alt, body, w_lo, ref_bytes, cds_end_axis)
        {
            return Some(anchor);
        }
    }
    let start = w_lo + piece.ref_start as i64;
    let end = w_lo + piece.ref_end as i64 - 1;
    Some(Anchor {
        region: body,
        start,
        end,
        alt,
        form: AnchorForm::Replacement,
    })
}

/// The `dup` anchor for a piece that is a tandem duplication, or `None`.
///
/// `duplication.md:18` — "when a variant can be described as a duplication, it
/// must be described as a duplication" — is a MUST, not a preference, so this is
/// tried before every other typing a pure insertion can take.
///
/// [`is_tandem_duplication`] already recognised the shape; until #1235's
/// sequence-first ladder it was used **only** to refuse, via the since-removed
/// `needs_unsupported_form`, and the `dup` had to come back from the
/// per-member pipeline — which types each member in isolation and so cannot
/// produce one for a duplication the derivation assembled out of several
/// members. It now builds as well as recognises.
///
/// The duplicated span is the reference run the payload repeats: the
/// `alt.len()` bases immediately 5' of the insertion point, which is exactly
/// what [`is_tandem_duplication`] compared the payload against.
///
/// # Why the span has to be checked against the preceding piece
///
/// Unlike every other form this module builds, a `dup` **names reference bases
/// its own piece does not claim** — its span reaches `alt.len()` bases back from
/// an insertion point that consumes nothing. In a multi-piece group that reach
/// can extend into the piece before it, and the result would be an allele whose
/// members contradict each other: `g.[10_12delinsX;9_12dup]` duplicates three
/// bases the sibling has just replaced.
///
/// Nothing downstream catches that. `GEdit::Dup` neither claims its span nor
/// reads the window through the applied cells — it copies from the untouched
/// `ref_bytes`, which is where the piece's own `alt` bytes came from — so the
/// round-trip guard's `reapplied == result` holds by construction. The test
/// helpers' disjointness checks go through `hgvs_to_spdi`, where a `dup` is a
/// zero-width interbase insert, so they see no overlap either. Before the `dup`
/// builder existed this whole class was refused wholesale by the since-removed
/// `needs_unsupported_form`; dropping that refusal is what makes it reachable.
///
/// No input reaching the partition has been shown to produce it — candidates
/// collapse to a single piece under `trim_common_flanks` — but "unreachable, I
/// think" is not a guard, and the fallback costs nothing: the piece is spelled as
/// the plain insertion it already is, which denotes the same bases.
fn duplication_anchor(
    piece: &Piece,
    body: Region,
    w_lo: i64,
    ref_bytes: &[u8],
    previous_ref_end: usize,
) -> Option<Anchor> {
    if !is_tandem_duplication(piece, ref_bytes) {
        return None;
    }
    let source_start = piece.ref_start.checked_sub(piece.alt.len())?;
    if source_start < previous_ref_end {
        return None;
    }
    Some(Anchor {
        region: body,
        start: w_lo + source_start as i64,
        end: w_lo + piece.ref_start as i64 - 1,
        // A `dup` names its source bases by span; there is nothing to insert.
        alt: Vec::new(),
        form: AnchorForm::Duplication,
    })
}

/// The single-position `delins` anchor for a derived insertion that has come to
/// rest outside the sequence's own bounds, or `None` when it has not.
///
/// The per-member pipeline owns this rewrite through
/// [`clamp_insertion_at_sequence_bounds`](crate::normalize); it is the reason
/// the sequence-first pass used to refuse every derivation that collapsed to one
/// pure insertion. An insertion resting immediately 5' of base 1, or immediately
/// 3' of the last base, has only one flanking anchor — `g.1ins…` is the form the
/// spec itself rejects (`DNA/insertion.md:95-101`) and `g.<len>_<len+1>ins…`
/// names a position past the end. On `m.` the tempting wraparound spelling is not
/// available either (#129/#1217).
///
/// **Which boundary this is, and which it is not.** `ends` is derived from
/// [`window_sequence_ends`], i.e. from the *entity's* own length, so this covers
/// exactly the **sequence** bounds — the contig ends on `g.`/`m.` (#1205, #1217)
/// and the transcript ends on `c.`/`n.`/`r.` (the transcript-bound half of
/// #1202). It does **not** cover the CDS/UTR *axis-change* boundaries at
/// `cds_start` / `cds_end` (#1170, #387, #1209): those have a representable far
/// side (`c.-1`, `c.*1`), so clamping there is a policy with carve-outs rather
/// than a representability rule, and `normalize_cds` keeps it — see
/// [`clamp_insertion_at_sequence_bounds`](crate::normalize)'s table.
///
/// On `c.`/`r.` the transcript's last base may sit in the **3'UTR**, and the
/// anchor this returns is still placed on `body` at `w_lo + offset`. That is
/// sound: the `AxisFrame` identity `sequence = axis + delta` is discontinuous
/// only at the *5'* end of the CDS (`c.-1` is `cds_start - 1`, since the axis has
/// no zero — the reason [`region_sequence_delta`] exists), and running *past*
/// `cds_end` is continuous, `c.<n>` denoting exactly the base `c.*(n - cds_len)`
/// does. The formatter then renders the canonical `*` form. Pinned by
/// `issue_1235_transcript_axes::cds_axis_clamps_a_derived_insertion_at_the_transcript_end`,
/// on a transcript with a real UTR at both ends.
///
/// `ends` says which of the fetched window's edges are the *entity's* edges; a
/// window that merely stopped where the fetch stopped must not be clamped at, or
/// a perfectly valid interior insertion is rewritten at the wrong place.
///
/// The rewrite itself is
/// [`boundary_delins_bases`](crate::normalize::boundary_delins_bases) — the same
/// coordinate identity the per-member clamp applies, shared rather than copied.
fn boundary_delins_anchor(
    piece: &Piece,
    alt: &[Base],
    body: Region,
    w_lo: i64,
    ref_bytes: &[u8],
    ends: SequenceEnds,
) -> Option<Anchor> {
    let (offset, side) = if ends.five_prime && piece.ref_start == 0 {
        (0usize, BoundarySide::FivePrime)
    } else if ends.three_prime && piece.ref_start == ref_bytes.len() {
        (ref_bytes.len().checked_sub(1)?, BoundarySide::ThreePrime)
    } else {
        return None;
    };
    // `boundary_delins_bases` numbers its anchor 1-based into the slice.
    let bases = boundary_delins_bases(ref_bytes, alt, offset as u64 + 1, side)?;
    let position = w_lo + offset as i64;
    Some(Anchor {
        region: body,
        start: position,
        end: position,
        alt: bases,
        form: AnchorForm::Replacement,
    })
}

/// The **axis-change** counterpart to [`boundary_delins_anchor`]: a derived pure
/// insertion resting exactly on the CDS/3'UTR junction is re-spelled as
/// `c.<cds_end>delins<ref[cds_end] ++ payload>` (#387, #1398).
///
/// [`boundary_delins_anchor`] deliberately covers only the *representability*
/// bounds, on the grounds that `normalize_cds` keeps the axis-change policy. For
/// a lone member that reasoning holds — the per-member pipeline re-clamps it. It
/// does **not** hold for a derivation, which is authoritative for the group and
/// whose output nothing re-clamps: `c.[11_12dup;16_17insC]` came back as
/// `c.18_*1insCAT`, and normalizing *that* produced `c.18delinsTCAT`, so the
/// first answer was not a fixed point and the two shuffle directions disagreed
/// about the second (#1398). The clamp therefore has to exist on both paths.
///
/// The identity is the same one `normalize_cds` applies and is shared, not
/// copied: inserting `A'` between `cds_end` and `cds_end + 1` is exactly
/// "delete `ref[cds_end]`, insert `ref[cds_end] ++ A'`", true for any `A'`.
///
/// **Spanning duplications are not clamped.** Per the spec a `dup` is the
/// priority form even when its span bridges the boundary (`c.9171_*1dup`), which
/// #401 pins. That carve-out needs no gate here: [`anchor_for_piece`] tries
/// [`duplication_anchor`] first, so a piece that a `dup` describes never reaches
/// this function — the same positive-list-`Insertion` effect `normalize_cds`
/// spells out explicitly, obtained from the call order instead.
///
/// `cds_end_axis` is the CDS end expressed on the member axis (`cds_end - delta`,
/// the `AxisFrame` identity), and is `None` for every axis without a CDS, which
/// makes this inert there. Positions run continuously past it — `c.<n>` beyond
/// the CDS length denotes the base `c.*(n - cds_len)` does — so the junction is
/// simply `cds_end_axis` and `cds_end_axis + 1`, with no discontinuity to correct
/// for. That is only true at the **3'** end; the 5' end has one (`c.-1` is
/// `cds_start - 1`, the reason `region_sequence_delta` exists), which is why this
/// has no `cds_start` mirror here.
fn cds_end_delins_anchor(
    piece: &Piece,
    alt: &[Base],
    body: Region,
    w_lo: i64,
    ref_bytes: &[u8],
    cds_end_axis: Option<i64>,
) -> Option<Anchor> {
    let cds_end_axis = cds_end_axis?;
    // A pure insertion's anchor spans `(start - 1, start)`, so it rests on the
    // junction exactly when its insertion point is the first 3'UTR base.
    let insertion_point = w_lo + i64::try_from(piece.ref_start).ok()?;
    if insertion_point != cds_end_axis + 1 {
        return None;
    }
    // Same 1-based-into-the-slice convention `boundary_delins_anchor` uses; the
    // anchor base is `ref[cds_end]`, one before the insertion point.
    let anchor_1b = u64::try_from(piece.ref_start).ok()?;
    let bases = boundary_delins_bases(ref_bytes, alt, anchor_1b, BoundarySide::ThreePrime)?;
    Some(Anchor {
        region: body,
        start: cds_end_axis,
        end: cds_end_axis,
        alt: bases,
        form: AnchorForm::Replacement,
    })
}

/// A piece is an inversion when its replacement is the reverse complement of
/// its own span and spans more than one nucleotide (`inversion.md:5,16` — a
/// 1-nt "inversion" is a substitution).
fn is_inversion(piece: &Piece, ref_bytes: &[u8]) -> bool {
    let span = &ref_bytes[piece.ref_start..piece.ref_end];
    span.len() > 1
        && span.len() == piece.alt.len()
        && reverse_complement_bytes(span).is_some_and(|rc| rc == piece.alt)
}

/// A piece is a tandem duplication when it is a pure insertion whose bases
/// repeat the reference immediately 5' of the insertion point.
fn is_tandem_duplication(piece: &Piece, ref_bytes: &[u8]) -> bool {
    if piece.ref_start != piece.ref_end || piece.alt.is_empty() {
        return false;
    }
    let Some(source_start) = piece.ref_start.checked_sub(piece.alt.len()) else {
        return false;
    };
    ref_bytes[source_start..piece.ref_start].eq_ignore_ascii_case(&piece.alt)
}

/// Reverse complement of a byte slice, uppercased, or `None` on a non-IUPAC
/// byte.
///
/// Built on [`complement_base`] rather than [`crate::sequence::reverse_complement`]
/// for the reason that function documents: the display helper's fallback arm is
/// `other => other`, so routing through it made this function total. It could
/// not return `None`, its own doc comment was therefore false, and
/// [`is_inversion`]'s `is_some_and` guard was unreachable — a non-base byte
/// would have been "complemented" to itself and folded into an inversion.
fn reverse_complement_bytes(bases: &[u8]) -> Option<Vec<u8>> {
    bases.iter().rev().map(|b| complement_base(*b)).collect()
}

// ----------------------------------------------------------------------------
// Sibling-crossing shift clamp (#1254)
// ----------------------------------------------------------------------------

/// A cis member's span on the shared coordinate axis, copied out of the
/// variant so the caller can mutate the member afterwards.
struct MemberSpan {
    /// The member's accession as written, used to decide whether two members
    /// sit on the same reference. Keeps any genomic-context wrapper, so
    /// members differing only in that wrapper are not treated as siblings.
    accession: String,
    /// The same accession reduced to the bare form a provider is keyed by, per
    /// `Accession::transcript_accession`. `accession` itself may carry a
    /// wrapper (`NC_000013.11(BRCA2)`) that no provider has an entry for, so a
    /// lookup through it would silently fail and skip the rewrite; every other
    /// provider lookup in this module already reduces first.
    provider_key: String,
    /// The region the span **starts** in. For a member wholly inside one region
    /// this is its region; for one crossing a boundary it is where the span
    /// begins, which is the axis its written coordinates are numbered on.
    region: Region,
    /// **Sequence** coordinates, 1-based inclusive — the representation in which
    /// any two members of one molecule are comparable and every length is right.
    ///
    /// The `c.`/`n.` axis is piecewise: three regions numbered by three
    /// independent integer sequences (`-n`, `n`, `*n`), so an axis coordinate
    /// does not order against one in another region — `c.*1` is 3' of `c.15`
    /// while `1 < 15` — and a span crossing a boundary has no meaningful axis
    /// length: `c.15_*1` is two bases, but `1 - 15` is `-14`. Converting each
    /// endpoint through `region_sequence_delta` removes both problems at the
    /// source, which is why this is the default representation and the axis one
    /// below has to be asked for by name (#1482).
    start: i64,
    end: i64,
    /// The coordinates as **written**, on the member's own axis.
    ///
    /// Only for rendering a position back onto that axis; `respell_as` is the
    /// sole consumer. Meaningless as a length or an ordering whenever
    /// `crosses_regions` is set, which is precisely why comparisons use the
    /// sequence pair above.
    axis_start: i64,
    axis_end: i64,
    /// Whether the two endpoints lie in different regions.
    crosses_regions: bool,
    /// Whether the edit consumes the reference bases under its span.
    claims_bases: bool,
    /// Whether an insertion landing flush against this edit is the adjacency the
    /// collapse pass merges — see [`merges_with_a_flush_insertion`].
    absorbs_a_flush_insertion: bool,
    /// Whether the edit must stop a sibling's shift at its boundary. Superset of
    /// `claims_bases` — see `blocks_sibling_shift`.
    blocks_shift: bool,
    /// For an insertion-like edit, the position its added sequence sits
    /// immediately 3' of; `None` for anything that consumes bases instead.
    junction: Option<i64>,
}

/// Read a member's span for the sibling-awareness passes, or `None` for a member
/// they cannot place (wrong kind, uncertain edit, offset / special position, or
/// a region the provider cannot convert onto the sequence axis).
///
/// # Why this takes a provider (#1482)
///
/// It used to route through `cis_axis_parts` -> `join_pos`, which refuses an
/// interval whose two endpoints are in different regions. `join_pos`' comment
/// justified that by claiming cross-region ranges "have no valid HGVS syntax".
/// **That is false.** `c.15_*1del` spans the CDS/3'UTR boundary and `c.-1_1del`
/// the 5'UTR/CDS boundary; both are ordinary descriptions, and deleting across a
/// stop codon is a common real variant. The comment even names `c.-1_1` as its
/// example of the impossible thing.
///
/// The refusal is defensible for the *merge* path, which coalesces members onto
/// one axis and has no representation for such a span. It was never defensible
/// here, and here it did not read as a refusal at all: this function returned
/// `None`, and every pass that consumes it drops a `None` through `filter_map`
/// rather than declining. A boundary-crossing member was therefore invisible
/// **as a sibling** across the whole module — it could not collide with a
/// duplication, block a sibling shift, bound a repeat tract, or take part in a
/// junction merge.
///
/// What that produced was an allele ferro's own parser rejects:
///
/// ```text
/// NM_TEST.1:c.[15_*1del;15_*1insCC]  ->  NM_TEST.1:c.[14_15dup;15_*1del]
/// ```
///
/// `c.15_*1insCC` alone correctly becomes `c.14_15dup` (`duplication.md:18`
/// makes `dup` mandatory when a variant can be described as one). Beside a
/// sibling deleting `c.15` the two claim the same base, and
/// `respell_colliding_duplications` — which exists to turn exactly that
/// duplication back into an insertion, a form that "claims nothing and so cannot
/// collide" — never saw the deletion.
///
/// # Why the whole span converts, rather than the collision test being patched
///
/// Patching one predicate fixes one symptom. The defect is that the *span
/// representation* cannot express these members, so every pass keyed on it is
/// blind in the same way. Converting here fixes them together and makes the
/// class unable to recur: a future pass gets a comparable span for free.
///
/// # Why the conversion is [`region_span_delta`] and not [`region_sequence_delta`]
///
/// Requiring a conversion here is a way to *lose* visibility as well as gain it,
/// and the two `n.` regions outside the transcript are where that bites:
/// `region_sequence_delta` refuses them because they hold no bases, which is
/// right for reading and wrong for comparing. Refusing here would not make a
/// pass conservative, it would make an `n.-5del` invisible to its siblings in
/// exactly the way this function exists to stop — measured, it un-fixed
/// `n.[-5=;-5del]` into an allele with two members on one position. So the
/// conversion used here widens to those regions with virtual positions; see
/// `region_span_delta` for why nothing can read bases through one.
///
/// What it does still decline — `Cds` on a record with no CDS, an inverted CDS —
/// every pass that compares members already declined for itself, so moving the
/// conversion here only stops each pass repeating it.
fn member_span<P: ReferenceProvider>(
    v: &HgvsVariant,
    kind: CisKind,
    provider: &P,
) -> Option<MemberSpan> {
    let (accession, (start_region, axis_start), (end_region, axis_end), edit) =
        member_axis_endpoints(v, kind)?;
    let provider_key = accession.transcript_accession();
    let to_sequence = |region: Region, coord: i64| -> Option<i64> {
        region_span_delta(region, &provider_key, provider)?.checked_add(coord)
    };
    let start = to_sequence(start_region, axis_start)?;
    let end = to_sequence(end_region, axis_end)?;
    Some(MemberSpan {
        accession: accession.full(),
        provider_key,
        region: start_region,
        start,
        end,
        axis_start,
        axis_end,
        crosses_regions: start_region != end_region,
        claims_bases: claims_reference_bases(edit),
        absorbs_a_flush_insertion: merges_with_a_flush_insertion(edit),
        blocks_shift: blocks_sibling_shift(edit),
        // Sequence coordinates in, sequence coordinate out: `junction_of` is
        // arithmetic on the two ends, so it carries whichever representation it
        // is given, and every consumer of `junction` compares it with `start` /
        // `end` above.
        junction: junction_of(edit, start, end),
    })
}

/// What [`member_axis_endpoints`] reads off one member: its accession, each
/// endpoint as `(region, written position)`, and the edit.
type MemberAxisEndpoints<'v> = (&'v Accession, (Region, i64), (Region, i64), &'v NaEdit);

/// Both endpoints of a member's location with **each keeping its own region**.
///
/// The cross-region counterpart of `cis_axis_parts`, which routes through
/// `join_pos` and so refuses exactly the members [`member_span`] must see.
fn member_axis_endpoints(v: &HgvsVariant, kind: CisKind) -> Option<MemberAxisEndpoints<'_>> {
    let (accession, start, end, edit) = match (kind, v) {
        (CisKind::Genome, HgvsVariant::Genome(g)) => (
            &g.accession,
            simple_genome_pos(&g.loc_edit.location.start)?,
            simple_genome_pos(&g.loc_edit.location.end)?,
            &g.loc_edit.edit,
        ),
        (CisKind::Mt, HgvsVariant::Mt(m)) => (
            &m.accession,
            simple_genome_pos(&m.loc_edit.location.start)?,
            simple_genome_pos(&m.loc_edit.location.end)?,
            &m.loc_edit.edit,
        ),
        (CisKind::Cds, HgvsVariant::Cds(c)) => (
            &c.accession,
            simple_cds_pos(&c.loc_edit.location.start)?,
            simple_cds_pos(&c.loc_edit.location.end)?,
            &c.loc_edit.edit,
        ),
        (CisKind::Tx, HgvsVariant::Tx(t)) => (
            &t.accession,
            simple_tx_pos(&t.loc_edit.location.start)?,
            simple_tx_pos(&t.loc_edit.location.end)?,
            &t.loc_edit.edit,
        ),
        (CisKind::Rna, HgvsVariant::Rna(rv)) => (
            &rv.accession,
            simple_rna_pos(&rv.loc_edit.location.start)?,
            simple_rna_pos(&rv.loc_edit.location.end)?,
            &rv.loc_edit.edit,
        ),
        _ => return None,
    };
    if !edit.is_certain() {
        return None;
    }
    Some((accession, start, end, edit.inner()?))
}

/// Whether an edit consumes the reference bases under its span.
///
/// Substitutions, deletions, delins and inversions replace or remove the bases
/// they cover, so a shift that carries one of them over another's bases changes
/// what the allele denotes. `NPaddedDeletion` (`delN[15]`) is a deletion that
/// states its length as a count rather than spelling the bases, and removes them
/// just the same; `SubstitutionNoRef` (`9>G`) replaces the base under its span
/// exactly as `9A>G` does, and omitting it made the clamp *spelling-sensitive* —
/// the same variant clamped or not depending on whether the reference base
/// happened to be written out. Insertions and duplications add sequence at a junction without
/// consuming a base, so a member landing flush against one is the *adjacency*
/// the collapse pass is meant to catch (#999, #1135) — they are not clamped
/// here. A duplication is nonetheless a 5' barrier, because the bases under its
/// span are what it copies; `blocks_sibling_shift` carries that, keeping this
/// predicate to the narrower "consumes bases" question its name asks.
///
/// `Repeat`/`MultiRepeat` claim the bases under their tract too: `g.1_9T[7]`
/// covers positions 1–9 whatever copy count it asserts. They were excluded until
/// the change could be measured (#1269), which left a repeat invisible as a
/// barrier and let a sibling's 3' shift land inside its tract (#1296).
fn claims_reference_bases(edit: &NaEdit) -> bool {
    matches!(
        edit,
        NaEdit::Substitution { .. }
            | NaEdit::SubstitutionNoRef { .. }
            | NaEdit::Deletion { .. }
            | NaEdit::NPaddedDeletion { .. }
            | NaEdit::Delins { .. }
            | NaEdit::Inversion { .. }
            | NaEdit::Repeat { .. }
            | NaEdit::MultiRepeat { .. }
    )
}

/// Whether an insertion placed flush against `edit` is the #999 adjacency the
/// collapse pass merges.
///
/// Strictly narrower than [`claims_reference_bases`], and narrower again than
/// "removes bases": only a **plain deletion** qualifies. `respell_at_gap`'s
/// terminal-overrun gate is the only caller, and the question it is really
/// asking is not what the sibling does to the reference but whether something
/// downstream is guaranteed to consume an out-of-range junction coordinate
/// before it renders (#1344 vs #1307). Nothing here may answer "yes" on a guess:
/// a wrong "yes" emits a position the sequence does not have, and that output
/// re-parses, so `FERRO_ASSERT_REPARSE` does not see it either. (The idempotency
/// oracle does, as it happens — see `tests/it/issue_1307_terminal_dup_respell.rs`
/// — but only for an input some test actually normalizes.)
///
/// Deliberately excluded, each measured on a 24-base contig with the member at
/// the last base:
///
/// - **`Delins`** — removes its bases, yet nothing merged the pair:
///   `g.[24dup;24delinsGG]` gave `g.[24delinsGG;24_25insA]`, the #1307 defect
///   with a different sibling. Removing bases is therefore not sufficient; the
///   collapse pass wants a deletion to absorb the payload into.
/// - **`NPaddedDeletion`** (`delN[15]`) — a deletion that states its length as a
///   count. Plausibly mergeable, but untested here, and the conservative branch
///   costs nothing.
/// - **`Repeat`/`MultiRepeat`** — whether they remove bases depends on the copy
///   count, so it is not a property of the edit kind at all.
///
/// Answering "no" for any of these yields `LeaveMemberUnchanged`, which is
/// always in range; the cost is only that a repairable collision goes
/// unrepaired.
fn merges_with_a_flush_insertion(edit: &NaEdit) -> bool {
    matches!(edit, NaEdit::Deletion { .. })
}

/// Whether `edit` must stop a *sibling's* shift at its boundary, even though it
/// may not consume the bases under its span.
///
/// EXPERIMENT (#1279/#1286): `claims_reference_bases` plus `Duplication`/
/// `DupIns`. A duplication reads its payload from the reference under its own
/// span, so a sibling deletion landing on those bases changes what the
/// duplication copies — it is a barrier even though it claims nothing.
///
/// Used for the *5'* half of `clamp_sibling_crossing_shifts` only; the 3' half
/// deliberately stays on `claims_reference_bases`. See the comment on its
/// `onto_bases` filter for why.
fn blocks_sibling_shift(edit: &NaEdit) -> bool {
    claims_reference_bases(edit)
        || matches!(edit, NaEdit::Duplication { .. } | NaEdit::DupIns { .. })
}

/// Pull back any cis member whose shift carried it over a sibling's bases.
///
/// The 3' rule (`general.md:41`) is stated per description, with no
/// allele-awareness, so a member is shifted to its *standalone* most-3'
/// position. When a sibling edits a base inside the tract the member rotates
/// through, that is not meaning-preserving: the member leapfrogs the sibling
/// and the pair now describes a different sequence. `g.[3_4del;9del]` over
/// `TAATATATATAATATATATT` shifted `3_4del` to `10_11del`, straight past the
/// `9del`, and the two then merged to `g.12_14del` — a well-formed,
/// warning-free description of *different* bases (#1254).
///
/// The invariant this restores: a member's **sweep window** — every position
/// between its pre-shift and post-shift span — must contain no base claimed by
/// another member. A shift is a rotation, so within a window no sibling edits,
/// every position in it describes the same sequence; the members then sit on
/// disjoint windows and compose exactly as the input did. Stopping at the last
/// position before the nearest sibling base is the most 3' choice that keeps
/// that true, which leaves the member *adjacent* to its sibling rather than
/// past it — and the caller's next merge pass coalesces the two
/// (`delins.md:16`). For the case above the clamp yields `[7_8del;9del]`, which
/// merges to `g.7_9del`: three contiguous deleted bases, the same sequence the
/// input denotes.
///
/// Both shift directions are handled — [`ShuffleDirection::FivePrime`] crosses
/// siblings the same way, in the other direction.
///
/// `before` is the member list as it entered per-member normalization and
/// supplies both the pull-back floor and the barrier positions; `after` is the
/// shifted list, mutated in place. The two must be index-aligned. Barriers are
/// read from a snapshot of `after` as well, so a sibling that shifted *into*
/// this member's window also blocks; taking them from snapshots rather than the
/// live list keeps the pass independent of member order.
///
/// [`ShuffleDirection::FivePrime`]: crate::normalize::ShuffleDirection
///
/// Siblings are compared against the moving member on the **sequence** axis,
/// which is what [`MemberSpan`] already carries (#1482). The two need not share
/// a region — a member shifting out of the CDS is bounded by siblings still
/// inside it (#1418) — and region axes are not comparable position-by-position
/// (`c.12` and `c.*1` are adjacent bases with unrelated numbers), so every
/// comparison below is on the converted coordinates rather than the written
/// ones. This pass used to re-convert each span itself through
/// [`region_sequence_delta`]; that conversion now happens once, in
/// [`member_span`].
pub(crate) fn clamp_sibling_crossing_shifts<P: ReferenceProvider>(
    before: &[HgvsVariant],
    after: &mut [HgvsVariant],
    phase: AllelePhase,
    uncertain: bool,
    provider: &P,
) {
    // An uncertain allele — `g.[(…)]` — asserts only that the members are in
    // cis, not where they sit, so repositioning one states something the input
    // did not. Every neighbouring transform in `normalize_allele` gates on this.
    if phase != AllelePhase::Cis || uncertain || after.len() < 2 || before.len() != after.len() {
        return;
    }
    let Some(kind) = cis_kind_of(&after[0]) else {
        return;
    };
    let pre: Vec<Option<MemberSpan>> = before
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();
    let post: Vec<Option<MemberSpan>> = after
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();

    for i in 0..after.len() {
        let (Some(b), Some(a)) = (pre[i].as_ref(), post[i].as_ref()) else {
            continue;
        };
        if !b.blocks_shift || !a.blocks_shift {
            continue;
        }
        // A shift that changes region is still a shift, and still has to respect
        // siblings (#1418). It used to be skipped here, which let a member cross
        // the CDS/3'UTR boundary — and every sibling on the way — unchecked:
        // `c.[11_12del;14del]` emitted `c.[14del;*1_*2del]`, straight past the
        // `14del`, denoting different bases.
        //
        // Regions cannot be compared position-by-position: `c.12` and `c.*1` are
        // adjacent bases with unrelated numbers. So the whole pass runs on the
        // **sequence** axis, which is the representation `MemberSpan` now carries
        // (#1482). Within one region that is a constant offset from the written
        // numbers and changes nothing, which is why the arithmetic below is
        // untouched; across regions it is what makes the comparison mean anything
        // at all.
        //
        // A region with no conversion — an `n.` position outside the transcript —
        // has no span at all, so it was already skipped by the `pre`/`post` reads
        // above.
        let ((b_start, b_end), (a_start, a_end)) = ((b.start, b.end), (a.start, a.end));
        // A member that did not move cannot have crossed anything.
        let delta = a_start - b_start;
        if delta == 0 {
            continue;
        }
        // A pure translation is the clear case: the whole span rotated.
        //
        // A member that *shrank* while its start moved is admitted too (#1281).
        // Merging an adjacent substitution and deletion yields a `delins` that
        // reduces to a plain deletion and shifts in the same pass — `g.
        // 9_10delinsA` becomes `g.1del` on a nine-`T` tract, where `delta` is -8
        // at the start but -9 at the end. Refusing that let the member cross
        // eight positions unchecked, straight onto a `1del` sibling, and
        // `g.[1del;9_10delinsA]` emitted `g.[1del;1del]` — two members claiming
        // one base, which the apply oracle declines. It still swept a window
        // this pass can bound, and the cap below keeps the pull inside it.
        //
        // A member that *grew* is still refused. An insertion canonicalising to
        // a duplication moves its span and its junction together, so bounding
        // the span mis-places the copy — measured on the 5' sweep, that turns
        // correct outputs into silently wrong ones (#1266, #1279).
        let translated = a_end - b_end == delta;
        let shrank = a_end - a_start < b_end - b_start;
        if !translated && !shrank {
            continue;
        }
        // Every sibling, from either snapshot, that shares this member's
        // coordinate line. A sibling blocks in one of two ways: it *claims*
        // bases, so this member may not land on them; or it is insertion-like
        // and sits at a *junction*, which this member may not sweep across —
        // sliding a deletion past the point where a sibling adds sequence
        // reorders the two edits and changes what the allele denotes (`g.
        // [258del;259_260insC]` emitted `g.[259_260insC;263del]`).
        //
        // Matched on accession alone, not on region: a sibling in the 3'UTR is
        // on this member's coordinate line just as much as one in the CDS, and
        // filtering by region was the other half of #1418 — it hid exactly the
        // siblings a boundary-crossing member sweeps over. Their spans are
        // already on the sequence axis, which is what makes them comparable; a
        // sibling whose region has no conversion has no span and so is absent.
        let siblings: Vec<&MemberSpan> = (0..pre.len())
            .filter(|&j| j != i)
            .flat_map(|j| [pre[j].as_ref(), post[j].as_ref()])
            .flatten()
            .filter(|s| s.accession == a.accession)
            .collect();

        let pull = if delta > 0 {
            // 3'-shift: the window is `[b.start, a.end]`. A sibling starting
            // beyond this member's original end but at or before its new end
            // was rotated over; stop one position short of the nearest.
            //
            // `claims_bases`, not `blocks_shift`: a duplication is a barrier in
            // the 5' direction only. Widening it here stops a member *reaching*
            // a downstream duplication, which is the adjacency the collapse
            // pass needs — measured, `g.[1004del;1008_1009insA]` then emits
            // `g.[1007del;1008dup]` instead of collapsing to `g.1009_1010=`,
            // because the insertion canonicalises to `1008dup` and the barrier
            // holds `1004del` off it (#1135). A duplication's end is already a
            // junction barrier in this direction via `across_junctions` below.
            let onto_bases = siblings
                .iter()
                .filter(|s| s.claims_bases)
                .map(|s| s.start)
                .filter(|&start| start > b_end && start <= a_end)
                .map(|start| start - 1);
            // A junction between `j` and `j+1` is crossed once this member's
            // end reaches `j + 1`; stopping at `j` leaves it flush against the
            // junction, which is the adjacency the collapse pass wants.
            let across_junctions = siblings
                .iter()
                .filter_map(|s| s.junction)
                .filter(|&j| j >= b_end && j < a_end);
            onto_bases
                .chain(across_junctions)
                .min()
                .map(|limit| a_end - limit)
        } else {
            // 5'-shift: mirror image, window `[a.start, b.end]`.
            let onto_bases = siblings
                .iter()
                .filter(|s| s.blocks_shift)
                .map(|s| s.end)
                .filter(|&end| end < b_start && end >= a_start)
                .map(|end| end + 1);
            // Mirror of the 3' junction barrier. A junction between `j` and
            // `j + 1` is crossed once this member's start reaches `j`; stopping
            // at `j + 1` leaves it flush against the junction on the 3' side,
            // which is the adjacency the collapse pass wants.
            //
            // Base-claiming members only. The two branches bound opposite edges
            // — 3' limits `a.end`, 5' limits `a.start` — and for a member that
            // carries its own junction those are not interchangeable. A `dup`'s
            // junction is its *end* (`junction_of`), so bounding its start
            // holds it back by `len - 1` more than the invariant needs: on
            // `TTAATATATAATAATATTAT`, `g.[4_5insC;5_6dup]` stops at its input
            // position instead of reaching the 5'-most legal `4_5dup`.
            // Sequence-preserving, so the exhaustive sweeps cannot see it —
            // measured at 905 de-canonicalised outputs across a 64,512-case
            // dup corpus, every one with `dup_len > 1`.
            //
            // Junction-carrying members are already bounded correctly, on the
            // right edge, by `clamp_sibling_crossing_junctions` — which runs
            // immediately after this pass against the same snapshot and carries
            // a `commutes` escape hatch this bound would otherwise override.
            let across_junctions = siblings
                .iter()
                .filter(|_| a.junction.is_none())
                .filter_map(|s| s.junction)
                .filter(|&j| j < b_start && j >= a_start)
                .map(|j| j + 1);
            onto_bases
                .chain(across_junctions)
                .max()
                .map(|limit| a_start - limit)
        };
        // Never move past where the member started: those positions were not
        // reachable by this shift and are not equivalent to it.
        let Some(pull) = pull.map(|p| {
            if delta > 0 {
                p.min(delta)
            } else {
                p.max(delta)
            }
        }) else {
            continue;
        };
        if pull != 0 {
            // `crosses_regions` joins the region test for the same reason the
            // test is there (#1482): both translates add `-pull` to the written
            // axis number, and a member whose two ends sit in different regions
            // has no single axis to add to — `c.15_*1` shifted 5' by one is
            // `c.14_15`, not `c.14_*0`. Such a member is *visible* here now,
            // which is the point, but it is repaired by the restore below rather
            // than by arithmetic on a number that means two things.
            if a.region == b.region && !a.crosses_regions && !b.crosses_regions {
                // A duplication blocks a sibling's shift without claiming bases,
                // so it reaches this pass too — and its payload is read from
                // under its own span, so it cannot simply be slid (#1280, #1292).
                match a.junction {
                    Some(_) => translate_junction_member(&mut after[i], -pull, kind, a, provider),
                    None => translate_member(&mut after[i], -pull, kind, a, provider),
                }
            } else {
                // The member crossed a region boundary — either by shifting into
                // a new one (#1418) or by spanning one to begin with (#1482) —
                // and the pull would carry it back across. Neither translate can
                // express that:
                // both add `-pull` to the axis number and then verify the result
                // landed in `from.region`, so a `c.*1` pulled seven bases 5'
                // would be written as `c.*-6` and reverted — silently leaving the
                // uncorrected output in place.
                //
                // Restore the pre-shift member instead. That is a *stronger* pull
                // than the bound requires, so it is always inside the swept
                // window and cannot cross anything itself; and the input spelling
                // is sequence-correct by construction, since it is what the
                // allele denoted before this pass moved it. The cost is
                // canonicalisation, not correctness: the member may sit 5' of the
                // most-3' position a boundary-aware translate would reach.
                //
                // Re-spelling it properly means a translate that converts between
                // region axes rather than adding to one. That is worth doing, and
                // is deliberately not done here — it is a change to shared
                // machinery every clamp uses, and this pass's job is to stop the
                // allele denoting the wrong bases.
                after[i] = before[i].clone();
            }
        }
    }
}

/// Re-spell as the plain edit it grew from any cis member whose repeat notation
/// spans a sibling's bases.
///
/// Covers all three sources, not deletions alone: a repeat grown from a
/// `Deletion` re-spells as a deletion, and one grown from a `Duplication` or an
/// `Insertion` as a duplication over the tract's 3'-most bases. Demoting an
/// insertion to a *duplication* rather than back to an insertion is deliberate
/// and load-bearing: a duplication's payload is read from the reference under
/// its span, so it is in phase at whatever position it lands on, whereas an
/// insertion carries a fixed literal that would need rotating to stay correct in
/// a non-homopolymer tract.
///
/// `deletion_to_repeat` (`repeated.md` B2) re-expresses a deletion inside a
/// tandem tract as a copy count over the **whole tract**, not over the deleted
/// bases: on a nine-`T` run, `g.1_2del` becomes `g.1_9T[7]`. That widened span
/// is correct for a lone variant, but in a cis allele it can swallow a sibling
/// — `g.[1_2del;4del]` emitted `g.[1_9T[7];9del]`, where `9del` sits inside
/// `1_9`. Overlapping members are malformed, and the pair denotes a different
/// sequence than the input.
///
/// The conversion happens per member, deep inside `normalize_na_edit`, with no
/// sibling context reaching it, so the decision is undone here instead: the
/// repeat is re-spelled as the equivalent 3'-most deletion of the same number
/// of bases, ending where the tract ends. That is the form the member would
/// have taken had B2 declined, and it is exactly equivalent — B2 only re-writes
/// notation, never the bases removed.
///
/// Restoring the deletion is what lets the rest of the pipeline finish the job:
/// the member is back to a span the caller's next merge pass can coalesce with
/// its sibling. `g.[1_2del;4del]` ends up as `g.1_9T[6]` — one member, three
/// bases removed from the tract, repeat notation restored now that there is no
/// sibling to span.
///
/// Visibility to [`clamp_sibling_crossing_shifts`] is no longer part of that:
/// since #1296 a repeat reports that it claims the bases under its tract, so
/// the clamp sees the widened span whether or not this pass has run. What this
/// pass still supplies is the *narrower, coalescible* spelling.
///
/// Runs before the clamp, so the clamp sees the deletion rather than the
/// repeat. Members whose pre-normalization form was not a plain deletion are
/// left alone: without a deleted-base count there is nothing to re-spell to.
///
/// A repeat that grew its tract by more bases than the tract holds has no
/// duplication form at all — the copy would have to start 5' of the tract — and
/// for a swallowed *junction* it is re-spelled as the insertion its growth
/// stands for instead (#1325). Running before the clamp is what makes that
/// worth doing: the clamp bounds the restored junction against the sibling.
pub(crate) fn demote_repeats_spanning_siblings<P: ReferenceProvider>(
    before: &[HgvsVariant],
    after: &mut [HgvsVariant],
    phase: AllelePhase,
    uncertain: bool,
    provider: &P,
) {
    if phase != AllelePhase::Cis || uncertain || after.len() < 2 || before.len() != after.len() {
        return;
    }
    let Some(kind) = cis_kind_of(&after[0]) else {
        return;
    };
    let pre: Vec<Option<MemberSpan>> = before
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();
    let post: Vec<Option<MemberSpan>> = after
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();

    for i in 0..after.len() {
        let (Some(b), Some(a)) = (pre[i].as_ref(), post[i].as_ref()) else {
            continue;
        };
        // Only a repeat is re-spellable, and only back into the kind of edit it
        // grew out of.
        if !matches!(
            cis_axis_parts(&after[i], kind).map(|(_, _, _, _, e)| e),
            Some(NaEdit::Repeat { .. })
        ) {
            continue;
        }
        let Some(source) = cis_axis_parts(&before[i], kind).and_then(|(_, _, _, _, e)| match e {
            // A deletion removes the bases under its span.
            NaEdit::Deletion { .. } => Some((RepeatSource::Removed, b.end - b.start + 1)),
            // A duplication adds a copy of the bases under its span.
            NaEdit::Duplication { .. } => Some((RepeatSource::Added, b.end - b.start + 1)),
            // An insertion adds its own sequence at a junction; its span is the
            // gap, so the length has to come from the sequence itself.
            NaEdit::Insertion { sequence } => sequence
                .len()
                .and_then(|n| i64::try_from(n).ok())
                .map(|n| (RepeatSource::Added, n)),
            _ => None,
        }) else {
            continue;
        };
        let (source, length) = source;
        if length <= 0 || b.region != a.region {
            continue;
        }
        // Does the tract actually swallow a sibling? Both snapshots, so a
        // sibling that moved into the tract counts too.
        //
        // Two ways to be swallowed. A sibling that *claims bases* is swallowed
        // when its span meets the tract. A sibling that claims none is still
        // swallowed when its **junction** falls strictly inside the tract
        // (#1287): a copy count over `1_9` cannot express another member adding
        // sequence in the middle of those nine bases, so the pair overlaps. The
        // 3' end is exclusive on purpose — a junction at `a.end` is flush against
        // the tract rather than inside it, which is the adjacency the collapse
        // pass exists to catch (#999, #1135), and the same strictness
        // `detect_insertion_overlaps` uses for an interior junction (#1276).
        let siblings = || {
            (0..pre.len())
                .filter(|&j| j != i)
                .flat_map(|j| [pre[j].as_ref(), post[j].as_ref()])
                .flatten()
                .filter(|s| s.region == a.region && s.accession == a.accession)
        };
        let spans_sibling_bases = siblings()
            .filter(|s| s.claims_bases)
            .any(|s| s.start <= a.end && s.end >= a.start);
        let spans_sibling_junction = siblings()
            .filter(|s| !s.claims_bases)
            .filter_map(|s| s.junction)
            .any(|junction| junction >= a.start && junction < a.end);
        let spans_sibling = spans_sibling_bases || spans_sibling_junction;
        if !spans_sibling {
            continue;
        }
        // Either way the equivalent edit sits on the 3'-most `length` bases of
        // the tract: removing them, or duplicating them.
        let start = a.end - length + 1;
        if start < a.start {
            // The member grew the tract by more bases than the tract holds, so
            // neither target form fits — a duplication over the tract can add
            // at most `tract` bases. Declining here is what left #1325's repeat
            // spanning its sibling, and the loss is larger than one spelling:
            // the first merge pass had already bounded this member correctly
            // (the clamp held both duplications short of the sibling's
            // junction, and the coalesce combined them), and the next pass
            // re-spelled the result as a tract-wide copy count. A repeat
            // carries no junction, so the clamp can no longer see it — this
            // demotion is the only route back to it.
            //
            // The growth is expressible as an insertion at the tract's 3'
            // junction, which is the payload `demote_coincident_tract_repeats`
            // already builds. Emitting it restores the junction and
            // `clamp_sibling_crossing_junctions`, immediately after, pulls it
            // back to the last position before the sibling.
            //
            // An insertion claims no bases and so blocks no sibling's shift, so
            // this route may only fire where that barrier has nothing left to
            // do. For a base-claiming sibling that means the tract holds it in
            // *both* snapshots: where it started and where this pass's
            // per-member normalization has put it.
            //
            // That is what separates #1394 from #1296. In #1296 the deletion
            // reaches the tract only after its own 3' shift (272_273 ->
            // 273_274), and `clamp_sibling_crossing_shifts` — which runs
            // immediately after this pass — pulls it back out to 272_273. The
            // overlap is the ordinary pre-clamp state the clamp exists to
            // repair, so demoting on it destroys the barrier one step before it
            // does its job, which is the defect #1296 fixed. In #1394 the
            // deletion is inside the tract at both 263 and 264: there is no
            // position outside it for the clamp to reach, so the barrier buys
            // nothing and the tract is simply swallowing a sibling.
            //
            // A swallowed junction needs no such care and keeps the plain test:
            // a junction member claims no bases, so there is nothing to release.
            let traps_sibling_bases = (0..pre.len()).filter(|&j| j != i).any(|j| {
                let overlaps = |s: &MemberSpan| {
                    s.claims_bases
                        && s.region == a.region
                        && s.accession == a.accession
                        && s.start <= a.end
                        && s.end >= a.start
                };
                pre[j].as_ref().is_some_and(overlaps) && post[j].as_ref().is_some_and(overlaps)
            });
            if !spans_sibling_junction && !traps_sibling_bases {
                continue;
            }
            let Some((payload, _)) = repeat_growth_as_insertion(&after[i], kind, a) else {
                continue;
            };
            respell_at_gap(
                &mut after[i],
                a,
                a.end,
                NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(Sequence::new(payload)),
                },
                provider,
                TerminalOverrun::RespellAtBoundary,
            );
            continue;
        }
        // A member spanning two regions is *seen* by the sibling tests above —
        // that is #1482 — but cannot be re-spelled here: `respell_as` writes one
        // pair of numbers onto one axis, and a tract running from the CDS into
        // the 3'UTR has no such pair. The `respell_at_gap` route above resolves
        // both ends through the transcript and does handle it.
        if a.crosses_regions {
            continue;
        }
        // The axis pair, not the sequence one: this writes positions back onto
        // the member's own axis, while `start`/`end` are sequence coordinates
        // (#1482). Within one region the two differ by a constant, so the
        // `start < a.start` test above and this write agree about the span.
        respell_as(
            &mut after[i],
            kind,
            a.region,
            a.axis_end - length + 1,
            a.axis_end,
            source,
            provider,
        );
    }
}

/// The bases a repeat added to its own tract, spelled as an insertion payload,
/// together with the unit it repeats.
///
/// A tract is tandem, so appending the added bases at its 3' end is just the
/// unit repeated — and the 3' end is where the tract's own shift rule places
/// them anyway. `g.262A[3]` over the one-base `A` at 262 added two copies, so
/// its payload is `AA` and its junction is 262.
///
/// `None` unless the repeat is a plain exact count over a whole number of units
/// that actually grew: a genotype count (`A[6][1]`), an inexact count, a VEP
/// trailing sequence, a partial unit, or a tract that lost copies all describe
/// something this equivalence does not cover, and fabricating a payload for
/// them would state bases nothing asserted.
fn repeat_growth_as_insertion(
    variant: &HgvsVariant,
    kind: CisKind,
    span: &MemberSpan,
) -> Option<(Vec<Base>, Vec<Base>)> {
    let NaEdit::Repeat {
        sequence: Some(unit),
        count: RepeatCount::Exact(count),
        additional_counts,
        trailing: None,
    } = cis_axis_parts(variant, kind)?.4
    else {
        return None;
    };
    if !additional_counts.is_empty() {
        return None;
    }
    let unit = unit.bases();
    let unit_len = i64::try_from(unit.len()).ok()?;
    let tract = span.end - span.start + 1;
    if unit_len <= 0 || tract <= 0 || tract % unit_len != 0 {
        return None;
    }
    let added = i64::try_from(*count).ok()? * unit_len - tract;
    if added <= 0 {
        return None;
    }
    let payload: Vec<Base> = unit
        .iter()
        .copied()
        .cycle()
        .take(usize::try_from(added).ok()?)
        .collect();
    Some((payload, unit.to_vec()))
}

/// Re-spell two or more repeats that grew the *same* tract as the insertions
/// they stand for, so the junction merge can combine their copies.
///
/// Each member is re-spelled per member, with no sibling context, and a member
/// that grows a tract is spelled as a copy count over the whole tract. Two
/// members that grew one tract therefore both name it, and the pair claims it
/// twice — an overlap with no resulting sequence at all:
///
/// ```text
/// reference  ("ACGT" x 64) + "CAGCCAGTCAGCGCATCAG" + ("ACGT" x 64)
///            the `A` at 262 is a one-base tract, between `C` at 261 and `G` at 263
///
/// g.[261_262insAA;262_263insAA]  ->  g.[262A[3];262A[3]]
///   each member grew the tract to three, and says so over the same one base
/// ```
///
/// This is the repeat-form analogue of #1286, and nothing else catches it. A
/// repeat has no junction (`junction_of` returns `None` for `NaEdit::Repeat`),
/// so [`coalesce_members_at_one_junction`] never considers the pair; and
/// [`demote_repeats_spanning_siblings`] re-spells a repeat back into a deletion
/// or a duplication *over the tract*, neither of which can express two bases
/// added to a one-base tract — it declines exactly this shape.
///
/// A repeat's growth is expressible, as an insertion at the tract's 3' junction.
/// Emitting that gives both members a junction and the existing merge finishes
/// the job: `g.262_263insAA` twice, coalesced to `insAAAA`, re-spelled `g.262A[5]`
/// on the next pass (#1316).
///
/// Only a whole *group* sharing one tract is re-spelled, and only when every
/// member of it yields a payload — a lone repeat is the correct spelling for a
/// member with no sibling on its tract, and half a group demoted would trade one
/// overlap for another. Requiring a shared unit as well as a shared span makes
/// the payloads powers of one word, so they commute and the coalesce below
/// cannot decline on order-ambiguity after this pass has already rewritten them.
///
/// Runs late, immediately before that coalesce, and **not** inside
/// [`demote_repeats_spanning_siblings`] despite the kinship: the demotion there
/// preserves the barrier a repeat presents to a sibling's shift, because a
/// deletion and a duplication both still block one. An insertion claims no
/// bases and blocks nothing, so demoting early released a sibling deletion that
/// #1296 had deliberately clamped out of the tract. Two repeats on one tract
/// have no such barrier to preserve — they occupy the same bases — so the
/// rewrite is safe here and only here.
pub(crate) fn demote_coincident_tract_repeats<P: ReferenceProvider>(
    members: &mut [HgvsVariant],
    phase: AllelePhase,
    uncertain: bool,
    provider: &P,
) {
    if phase != AllelePhase::Cis || uncertain || members.len() < 2 {
        return;
    }
    let Some(kind) = cis_kind_of(&members[0]) else {
        return;
    };
    let spans: Vec<Option<MemberSpan>> = members
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();
    // A member participates only if it is a repeat that grew its tract by a
    // whole number of units; `growth` carries the payload and the unit.
    let growth: Vec<Option<(Vec<Base>, Vec<Base>)>> = (0..members.len())
        .map(|i| {
            let span = spans[i].as_ref()?;
            repeat_growth_as_insertion(&members[i], kind, span)
        })
        .collect();

    // Group by the tract itself — accession, region, span — plus the unit.
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for i in 0..members.len() {
        let (Some(span), Some((_, unit))) = (spans[i].as_ref(), growth[i].as_ref()) else {
            continue;
        };
        match groups.iter_mut().find(|at| {
            let first = at[0];
            spans[first].as_ref().is_some_and(|s| {
                s.accession == span.accession
                    && s.region == span.region
                    && s.start == span.start
                    && s.end == span.end
            }) && growth[first].as_ref().is_some_and(|(_, u)| u == unit)
        }) {
            Some(at) => at.push(i),
            None => groups.push(vec![i]),
        }
    }

    for at in &groups {
        if at.len() < 2 {
            continue;
        }
        // All or nothing. `respell_at_gap` reverts a member it cannot place, and
        // a group left half re-spelled is one insertion overlapping one repeat —
        // the same defect in a new spelling, and one no later pass looks for.
        let originals: Vec<HgvsVariant> = at.iter().map(|&i| members[i].clone()).collect();
        let rewritten = at.iter().enumerate().all(|(slot, &i)| {
            let (Some(span), Some((payload, _))) = (spans[i].as_ref(), growth[i].as_ref()) else {
                return false;
            };
            let edit = NaEdit::Insertion {
                sequence: InsertedSequence::Literal(Sequence::new(payload.clone())),
            };
            respell_at_gap(
                &mut members[i],
                span,
                span.end,
                edit,
                provider,
                TerminalOverrun::RespellAtBoundary,
            );
            members[i] != originals[slot]
        });
        if !rewritten {
            for (slot, &i) in at.iter().enumerate() {
                members[i] = originals[slot].clone();
            }
        }
    }
}

/// Pull back any cis member whose *junction* was carried over a sibling's bases.
///
/// [`clamp_sibling_crossing_shifts`] governs edits that consume the bases under
/// their span. An insertion or duplication consumes none — it adds sequence at
/// the junction between two bases — so it is excluded there, and deliberately:
/// a member landing flush against one is the adjacency the collapse pass exists
/// to catch (#999, #1135).
///
/// It can still cross. The junction itself 3'-shifts through a tract, and when
/// a sibling substitutes or deletes a base inside that tract, moving the
/// junction past it changes what the allele denotes:
///
/// ```text
/// reference    ACAAAAAAAACGTACGTACG        A-run at 3-10
/// g.[2_3insA;5A>G]  ->  g.[5A>G;10dup]     junction moved from 2|3 to 10|11
/// input applied     ACAAAGAAAAACGTACGTACG
/// output applied    ACAAGAAAAAACGTACGTACG  ← the substituted base moved
/// ```
///
/// No overlap, no warning, well-formed — and wrong. The constraint is that the
/// junction may reach a sibling's 5' edge but not pass it: for a sibling
/// claiming `[s, e]`, a 3'-shifted junction must stay at or below `s - 1`, and
/// a 5'-shifted one at or above `e`. That is exactly the adjacency #999 needs —
/// its insertion lands at `306|307` against a substitution at `307`, which is
/// the last position the rule permits — so the collapse still fires.
///
/// A sibling's own junction bounds it as well, not only a sibling's bases
/// (#1290): two junction-occupying members add their payloads at two interbase
/// points, and their relative order is the only thing fixing the order of those
/// payloads, so sweeping past one reorders them. That bound stops one position
/// short of the sibling's junction unless the two payloads commute.
///
/// Both shift directions are governed. The 5' half was missing until #1267 —
/// which is the shape where the sibling is *upstream* and the junction travels
/// toward it — and the two bounds needed it independently: base-claiming
/// siblings accounted for half the measured defect, junction-occupying siblings
/// the other half.
///
/// A junction sits immediately 3' of one position: `end` for a duplication
/// (the copy follows the duplicated bases), `start` for an insertion (whose
/// span is the gap `start_end`).
pub(crate) fn clamp_sibling_crossing_junctions<P: ReferenceProvider>(
    before: &[HgvsVariant],
    after: &mut [HgvsVariant],
    phase: AllelePhase,
    uncertain: bool,
    provider: &P,
) {
    if phase != AllelePhase::Cis || uncertain || after.len() < 2 || before.len() != after.len() {
        return;
    }
    let Some(kind) = cis_kind_of(&after[0]) else {
        return;
    };
    let pre: Vec<Option<MemberSpan>> = before
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();
    let post: Vec<Option<MemberSpan>> = after
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();
    // The members `post` was measured from. The loop below translates `after`
    // in place, so by the time member `i` inspects a sibling `j < i` that has
    // already moved, `after[j]` no longer matches `post[j]` — pairing them
    // reads a duplication's payload from bases under a span it has left, which
    // is the #1304 mismatch one iteration later rather than one snapshot over.
    // Keeping the pair together also keeps the barrier evaluated against the
    // state the whole pass started from, so the result does not depend on the
    // order the members happen to be visited in.
    let post_snapshot: Vec<HgvsVariant> = after.to_vec();

    for i in 0..after.len() {
        let (Some(b), Some(a)) = (pre[i].as_ref(), post[i].as_ref()) else {
            continue;
        };
        // Taken off the spans rather than re-derived from the parsed locations,
        // so that these two junctions and every sibling bound below are on one
        // axis. `MemberSpan::junction` is `junction_of` applied to the same edit,
        // on **sequence** coordinates (#1482); reading it from `cis_axis_parts`
        // instead yields the written axis number, which no longer compares
        // against `s.start`/`s.end`.
        let (Some(before_junction), Some(after_junction)) = (b.junction, a.junction) else {
            continue;
        };
        if before_junction == after_junction || b.region != a.region {
            continue;
        }
        let siblings: Vec<&MemberSpan> = (0..pre.len())
            .filter(|&j| j != i)
            .flat_map(|j| [pre[j].as_ref(), post[j].as_ref()])
            .flatten()
            .filter(|s| s.claims_bases && s.region == a.region && s.accession == a.accession)
            .collect();

        // Both directions, and the rule is the same in each: a junction may
        // reach a sibling's edge but not pass it (#1267).
        //
        // The 5' half was absent, and its absence read as deliberate — #1259
        // measured a 5' mirror turning 80 previously-correct outputs silently
        // wrong, on the reasoning that a duplication's span and junction move
        // together so bounding the junction mis-places the copy. But that mirror
        // bounded the junction against a sibling the member was moving *away*
        // from, not one it was moving *toward*; it was a different rule, and its
        // failure does not bear on this one.
        //
        // Restricting the 5' half by input spelling — admitting only a member the
        // input wrote as an insertion — was tried and refuted by measurement.
        // Over 73,376 duplication-mover cases with an upstream sibling, bounding
        // every junction mover leaves **zero** sequence changes while the
        // restricted form leaves **12**, all of them a `dup` whose junction swept
        // past an upstream sibling's junction. Deciding this from the input's
        // encoding was both more complex and less correct, so the direction is
        // simply mirrored.
        let five_prime = after_junction < before_junction;
        // The junction crossed a base-claiming sibling if it started 5' of the
        // sibling's first base and ended at or past it. Mirrored under a 5'
        // shift: the junction started 3' of the sibling's *last* base and ended
        // at or before it, and must stay at or above that base.
        let onto_bases: Vec<i64> = if five_prime {
            siblings
                .iter()
                .filter(|s| s.end > after_junction && s.end <= before_junction)
                .map(|s| s.end)
                .collect()
        } else {
            siblings
                .iter()
                .filter(|s| s.start > before_junction && s.start <= after_junction)
                .map(|s| s.start - 1)
                .collect()
        };
        // It can equally cross another *junction* (#1290). Two junction-
        // occupying members add their payloads at two interbase points, and
        // their relative order is the only thing fixing the order of those
        // payloads — so a junction sweeping past another reorders them and
        // changes what the allele denotes, while leaving both members
        // well-formed and disjoint. Nothing else catches it: the overlap
        // detector sees two distinct junctions, and the clamp above sees a
        // sibling with no bases to land on.
        //
        // Stop one position **short** of the sibling's junction, not at it.
        // Landing on it makes the two share a junction, which is the case with
        // no defined payload order at all — `coalesce_members_at_one_junction`
        // merges that only when the payloads commute, and declines otherwise,
        // which would leave an overlap here rather than the correct ordered
        // pair. One short keeps them distinct and in their original order.
        //
        // Only a sibling whose payload does not *commute* with this member's.
        // Two payloads that commute have no observable order, so crossing is
        // harmless, and letting them settle on one junction is better: they then
        // merge into a single member (#1286), the canonical form the
        // sequence-first derivation would also reach. Barring the crossing
        // unconditionally would split `g.[258_259insA;259_260insA]` into
        // `g.[262dup;263dup]` — correct, but two members where one will do — and
        // would leave the three-insertion case at two members for the same
        // reason, since a merged `AA` is not *equal* to a remaining `A` even
        // though it commutes with it.
        // Each span is carried with the variant it was read from. A duplication
        // reads its payload from the reference *under its own span*, so pairing
        // the post-normalization edit with the pre-normalization span reads
        // bases belonging to neither: `g.[260_261insGA;261_262insA;264del]` had
        // the sibling's `insA` span measured against its own normalized
        // `265dup`, which reported the payload `GA` rather than `A`. `GA`
        // commutes with this member's `GA`, so the barrier vanished and the
        // member swept past it (#1304).
        // Commuting payloads are exempt from the barrier, but only as far as the
        // sibling's **settled** junction. Commuting says `p ++ q == q ++ p` — the
        // two payloads have no observable order *once they meet*. It says
        // nothing about a member sweeping past a sibling that stays put and
        // landing beyond it, where reference bases now separate them differently:
        //
        //     core "CAGAAGATGAATAA"
        //     g.[263_264insTG;264_265insTG]
        //       the second member does not move; the first 3'-shifts 263 -> 265
        //       intended  A T G T T G G
        //       emitted   A T T G T G G          (#1308)
        //
        // So a commuting sibling bounds this member at the junction it settles
        // on — landing *on* it is what lets the pair merge (#1286) — and only
        // its post-normalization position counts, since where it started is
        // exactly the order that commuting makes unobservable. A sibling whose
        // payload does not commute keeps the stricter rule above: one short, and
        // from either snapshot.
        let payload = junction_payload(&after[i], kind, a, provider);
        // Every junction compared below is already on the sequence axis (#1482),
        // so no per-region offset is applied here. A region with no conversion —
        // an `n.` position outside the transcript — has no span at all and so
        // never reaches this point, which is the same conservative answer the
        // per-pass conversion used to give.
        let across_junctions = (0..pre.len())
            .filter(|&j| j != i)
            .flat_map(|j| {
                [
                    (false, &before[j], pre[j].as_ref()),
                    (true, &post_snapshot[j], post[j].as_ref()),
                ]
            })
            .filter_map(|(settled, variant, span)| span.map(|s| (settled, variant, s)))
            .filter(|(_, _, s)| {
                !s.claims_bases && s.region == a.region && s.accession == a.accession
            })
            .filter_map(|(settled, variant, s)| {
                let junction = s.junction?;
                // Was it swept over? Under a 3' shift the sibling's junction
                // must lie above where this one started and at or below where it
                // ended; under a 5' shift, mirrored.
                let swept = if five_prime {
                    junction >= after_junction && junction < before_junction
                } else {
                    junction > before_junction && junction <= after_junction
                };
                if !swept {
                    return None;
                }
                // Commuting is tested against the payload this member would
                // carry **at that junction**, not the one it carries now. The
                // two differ: landing there rotates the payload into phase, and
                // a rotation is exactly what can destroy the identity —
                // `260_261insAC` moved onto junction 261 becomes `insCA`, and
                // `CA` does not commute with the `AC` already sitting there,
                // though `AC` and `AC` did. Sharing a junction with a payload it
                // does not commute with leaves the pair's order undefined, and
                // the merge downstream concatenates them in rendered order:
                // `g.[260_261insAC;261_262insAC]` became `g.261_262insACCA`,
                // denoting different bases (#1312).
                let commutes = match (
                    payload
                        .as_ref()
                        .and_then(|mine| {
                            payload_at_junction(
                                mine,
                                after_junction,
                                junction,
                                &a.provider_key,
                                provider,
                            )
                        })
                        .as_deref(),
                    junction_payload(variant, kind, s, provider),
                ) {
                    (Some(landed), Some(theirs)) => payloads_commute(landed, &theirs),
                    // A payload that cannot be read, or cannot legally reach the
                    // junction, blocks: refusing to cross is the conservative
                    // answer when the reordering cannot be shown to be harmless.
                    _ => false,
                };
                match commutes {
                    // Where it settles, and it may be met rather than only
                    // approached. Its pre-normalization position is exactly the
                    // order commuting makes unobservable, so it is not a bound.
                    true => settled.then_some(junction),
                    // One short of it, on the side this member came from.
                    false if five_prime => Some(junction + 1),
                    false => Some(junction - 1),
                }
            });
        // The most restrictive bound: the lowest ceiling under a 3' shift, the
        // highest floor under a 5' one.
        let bounds = onto_bases.into_iter().chain(across_junctions);
        let limit = if five_prime {
            bounds.max()
        } else {
            bounds.min()
        };
        let Some(limit) = limit else {
            continue;
        };
        let Some(payload) = payload.as_ref() else {
            continue;
        };
        // Never move past where the junction started, and never onto a junction
        // this member's payload is out of phase with: the walk back finds the
        // most extreme position it can legally occupy in the direction it was
        // travelling. `before_junction` always can, since that is where the
        // shift began.
        let legal = |j: i64| {
            payload_at_junction(payload, after_junction, j, &a.provider_key, provider).is_some()
        };
        let destination = if five_prime {
            let floor = limit.min(before_junction).max(after_junction);
            (floor..=before_junction).find(|&j| legal(j))
        } else {
            let ceiling = limit.max(before_junction).min(after_junction);
            (before_junction..=ceiling).rev().find(|&j| legal(j))
        };
        let Some(destination) = destination else {
            continue;
        };
        let delta = destination - after_junction;
        if delta != 0 {
            translate_junction_member(&mut after[i], delta, kind, a, provider);
        }
    }
}

/// The position a junction-inserting edit sits immediately 3' of, or `None` for
/// an edit that is not junction-inserting.
fn junction_of(edit: &NaEdit, start: i64, end: i64) -> Option<i64> {
    match edit {
        // The copy follows the duplicated bases. `DupIns` is a duplication that
        // also carries inserted sequence, so its junction sits in the same
        // place: omitting it would leave a `DupIns` sibling invisible as a
        // barrier and let a 3'-shifting member sweep across it.
        NaEdit::Duplication { .. } | NaEdit::DupIns { .. } => Some(end),
        // The span is the gap `start_end`, so the sequence lands after `start`.
        NaEdit::Insertion { .. } => Some(start),
        _ => None,
    }
}

/// Which plain edit a tract-wide repeat was standing in for.
#[derive(Clone, Copy)]
enum RepeatSource {
    /// The tract lost copies — re-spell as a deletion.
    Removed,
    /// The tract gained copies — re-spell as a duplication.
    Added,
}

impl RepeatSource {
    fn edit(self) -> NaEdit {
        match self {
            RepeatSource::Removed => NaEdit::Deletion {
                sequence: None,
                length: None,
            },
            RepeatSource::Added => NaEdit::Duplication {
                sequence: None,
                length: None,
                uncertain_extent: None,
            },
        }
    }
}

/// Replace `variant`'s edit with `source`'s plain form over `[start, end]`.
///
/// `start` and `end` are positions on the member's **own axis** — the numbers
/// that get written — not the sequence coordinates [`MemberSpan`] carries
/// (#1482). The caller supplies them from `axis_start`/`axis_end`.
///
/// Reverts on any axis this pass cannot rewrite, or if the result does not land
/// exactly where intended on the same region.
fn respell_as<P: ReferenceProvider>(
    variant: &mut HgvsVariant,
    kind: CisKind,
    region: Region,
    start: i64,
    end: i64,
    source: RepeatSource,
    provider: &P,
) {
    let original = variant.clone();
    let placed = match variant {
        HgvsVariant::Genome(g) => place_genome(&mut g.loc_edit, start, end, source),
        HgvsVariant::Mt(m) => place_genome(&mut m.loc_edit, start, end, source),
        HgvsVariant::Cds(c) => {
            place_signed(&mut c.loc_edit, start, end, source, |p: &mut CdsPos, v| {
                p.base = v;
            })
        }
        HgvsVariant::Tx(t) => {
            place_signed(&mut t.loc_edit, start, end, source, |p: &mut TxPos, v| {
                p.base = v;
            })
        }
        HgvsVariant::Rna(r) => {
            place_signed(&mut r.loc_edit, start, end, source, |p: &mut RnaPos, v| {
                p.base = v;
            })
        }
        _ => false,
    };
    // Read back on the axis it was written on, and required to sit in one
    // region: a write that landed with its two ends in different regions is not
    // the span this function was asked for, and before #1482 `member_span`
    // reported such a result as `None` and so reverted it. Keeping that refusal
    // explicit is what stops the wider span reader from quietly widening what
    // this pass is willing to write.
    let landed = placed
        && member_span(variant, kind, provider).is_some_and(|s| {
            !s.crosses_regions && s.region == region && s.axis_start == start && s.axis_end == end
        });
    if !landed {
        *variant = original;
    }
}

/// Set a genomic `LocEdit` to `source`'s plain form over `[start, end]`.
fn place_genome(
    loc_edit: &mut LocEdit<Interval<GenomePos>, NaEdit>,
    start: i64,
    end: i64,
    source: RepeatSource,
) -> bool {
    let (Ok(start), Ok(end)) = (u64::try_from(start), u64::try_from(end)) else {
        return false;
    };
    if !set_endpoints(
        &mut loc_edit.location,
        |p: &mut GenomePos, v: u64| p.base = v,
        start,
        end,
    ) {
        return false;
    }
    loc_edit.edit = Mu::Certain(source.edit());
    true
}

/// Set a signed-axis `LocEdit` to `source`'s plain form over `[start, end]`.
fn place_signed<P, L>(
    loc_edit: &mut LocEdit<Interval<P>, NaEdit>,
    start: i64,
    end: i64,
    source: RepeatSource,
    assign: L,
) -> bool
where
    L: Fn(&mut P, i64) + Copy,
{
    if start == 0 || end == 0 {
        return false; // no `c.0` / `n.0` position exists
    }
    if !set_endpoints(&mut loc_edit.location, assign, start, end) {
        return false;
    }
    loc_edit.edit = Mu::Certain(source.edit());
    true
}

/// Assign both certain endpoints of an interval, refusing anything else.
fn set_endpoints<P, V, L>(interval: &mut Interval<P>, assign: L, start: V, end: V) -> bool
where
    L: Fn(&mut P, V),
    V: Copy,
{
    for (boundary, value) in [(&mut interval.start, start), (&mut interval.end, end)] {
        match boundary {
            UncertainBoundary::Single(Mu::Certain(pos)) => assign(pos, value),
            _ => return false,
        }
    }
    true
}

/// Move `variant`'s span `delta` positions along its axis (negative is 5').
///
/// Reverts and leaves the variant untouched unless the result lands exactly
/// where intended on the *same* region — the guard against an axis whose
/// integer coordinate is not simply translatable (a `c.` span that would cross
/// into the 5'UTR or `*`-region, a genomic base that would underflow).
fn translate_member<P: ReferenceProvider>(
    variant: &mut HgvsVariant,
    delta: i64,
    kind: CisKind,
    from: &MemberSpan,
    provider: &P,
) {
    // A member whose ends sit in different regions has no single axis number to
    // add `delta` to (#1482). Adding it to both writes a span that is not a
    // translation of this one: `c.15_*1` moved 3' by one becomes `c.16_*2`, and
    // `c.16` is a coding position on a CDS that ends at 15. The sequence pair
    // *does* translate, so the `landed` check below would accept it — the
    // refusal has to be explicit rather than left to the verification.
    if from.crosses_regions {
        return;
    }
    let original = variant.clone();
    let moved = match variant {
        HgvsVariant::Genome(g) => {
            translate_interval(&mut g.loc_edit.location, |p: &mut GenomePos| {
                shift_unsigned(&mut p.base, delta)
            })
        }
        HgvsVariant::Mt(m) => translate_interval(&mut m.loc_edit.location, |p: &mut GenomePos| {
            shift_unsigned(&mut p.base, delta)
        }),
        HgvsVariant::Cds(c) => translate_interval(&mut c.loc_edit.location, |p: &mut CdsPos| {
            shift_signed(&mut p.base, delta)
        }),
        HgvsVariant::Tx(t) => translate_interval(&mut t.loc_edit.location, |p: &mut TxPos| {
            shift_signed(&mut p.base, delta)
        }),
        HgvsVariant::Rna(r) => translate_interval(&mut r.loc_edit.location, |p: &mut RnaPos| {
            shift_signed(&mut p.base, delta)
        }),
        _ => false,
    };
    let landed = moved
        && member_span(variant, kind, provider).is_some_and(|s| {
            s.region == from.region && s.start == from.start + delta && s.end == from.end + delta
        });
    if !landed {
        *variant = original;
    }
}

/// Apply `step` to both certain endpoints of an interval, refusing any
/// boundary this pass cannot move (uncertain, unknown, or a range endpoint).
fn translate_interval<T>(interval: &mut Interval<T>, mut step: impl FnMut(&mut T) -> bool) -> bool {
    let mut moved = true;
    for boundary in [&mut interval.start, &mut interval.end] {
        match boundary {
            UncertainBoundary::Single(Mu::Certain(pos)) => moved &= step(pos),
            _ => return false,
        }
    }
    moved
}

/// Shift a 1-based unsigned coordinate, refusing a move off the 5' end.
fn shift_unsigned(base: &mut u64, delta: i64) -> bool {
    match base.checked_add_signed(delta) {
        Some(shifted) if shifted >= 1 => {
            *base = shifted;
            true
        }
        _ => false,
    }
}

/// Shift a signed axis coordinate, refusing a move onto the non-existent `0`.
fn shift_signed(base: &mut i64, delta: i64) -> bool {
    match base.checked_add(delta) {
        Some(shifted) if shifted != 0 => {
            *base = shifted;
            true
        }
        _ => false,
    }
}

/// Re-spell a duplication that collides with a base-claiming sibling as the
/// equivalent insertion.
///
/// A duplication's span *is* the base it copies, so `g.262dup` claims position
/// 262 — and a sibling deletion covering 262 makes the pair contradictory.
/// ferro's own parser rejects the result ("Self-cancelling allele: variants at
/// index 0 (del) and 1 (dup) describe overlapping reference positions"), so
/// normalization emits a description it cannot read back.
///
/// The collision is a *spelling* problem, not a positional one. `Xdup` and
/// `X_X+1ins<ref[X]>` denote exactly the same edit — a copy of the reference
/// bases inserted at the junction 3' of them — but the insertion form claims no
/// base, so it does not collide. Re-spelling therefore fixes the output without
/// moving anything.
///
/// This is why the positional clamps cannot solve it. Pulling the deletion one
/// base short of the junction also works, and breaks #1135, where a deletion
/// must reach the junction for the self-cancelling del+ins merge to fire.
/// Telling those apart needs the reference; so does writing the inserted bases,
/// which is why this pass takes a provider where the clamps do not.
///
/// Runs **inside** the fixed-point loop, last of the four sibling passes. Run
/// once after the loop instead, the result is not a fixed point: the `del` and
/// `ins` this re-spelling exposes cancel further on the next pass
/// (`g.[261_262del;262_263insA;…]` reduces to `g.[261del;…]`), so the loop has
/// to see the re-spelled form to settle on that reduction.
///
/// Runs on **every** axis. The pass reads the duplicated bases over the
/// member's own coordinates, so all it needs is the offset from those
/// coordinates onto the sequence the provider serves —
/// [`region_sequence_delta`], which answers per *region* and so has an answer
/// for the CDS-relative axes too.
///
/// This was gated to the genomic axes, then to those plus `n.` (#1315), on the
/// grounds that only they index the served sequence directly. That was true of
/// the arithmetic as it stood and is no longer a reason to refuse: `c.` and
/// `r.` are CDS-relative, and keying the offset on the region rather than on
/// the axis makes each of `FivePrimeUtr` / `Cds` / `ThreePrimeUtr`
/// individually affine — see [`region_sequence_delta`] for why the earlier
/// single-delta attempt recorded in #1284 could not work.
///
/// There is no axis gate left, because the conversion itself is the gate: a
/// region it cannot place (`n.` outside the transcript, a `c.` on a record
/// with no CDS) yields `None` and the member is skipped.
pub(crate) fn respell_colliding_duplications<P: ReferenceProvider>(
    members: &mut [HgvsVariant],
    phase: AllelePhase,
    uncertain: bool,
    provider: &P,
) {
    if phase != AllelePhase::Cis || uncertain || members.len() < 2 {
        return;
    }
    let Some(kind) = cis_kind_of(&members[0]) else {
        return;
    };
    let spans: Vec<Option<MemberSpan>> = members
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();

    for i in 0..members.len() {
        let Some(dup) = spans[i].as_ref() else {
            continue;
        };
        if !matches!(
            cis_axis_parts(&members[i], kind).map(|(_, _, _, _, e)| e),
            Some(NaEdit::Duplication { .. })
        ) {
            continue;
        }
        let same_molecule =
            |s: &&MemberSpan| s.region == dup.region && s.accession == dup.accession;
        let siblings = || {
            (0..spans.len())
                .filter(|&j| j != i)
                .filter_map(|j| spans[j].as_ref())
                .filter(same_molecule)
        };
        // A sibling that claims bases collides when its span meets the
        // duplication's — the two then name the same positions.
        let claims_the_bases = siblings()
            .filter(|s| s.claims_bases)
            .any(|s| s.start <= dup.end && s.end >= dup.start);
        // A sibling that claims none still collides when its **junction** falls
        // strictly between the duplication's two ends (#1320).
        //
        // A duplication *sorts* by its span's start and *acts* at its span's
        // end — it copies the bases under the span and lands them at the
        // junction 3' of them. Those two positions are the same for a one-base
        // dup and diverge as the span widens, and a sibling junction strictly
        // between them falls on the wrong side of both: the duplication is
        // rendered first, because 263 precedes 264, while the bases it adds
        // arrive at 266, after the sibling's. The pair then reads as out of
        // order however it is spelled. Re-spelling settles it, because the
        // insertion form sorts where it acts.
        //
        // **Both** ends are exclusive, and each for its own reason:
        //
        // - at `dup.end` the junction is the duplication's own landing junction,
        //   the shared-junction case `coalesce_members_at_one_junction` merges
        //   (#1286) — re-spelling would pull it out from under that merge;
        // - at `dup.start` the duplication sorts to the same position as the
        //   sibling, so `junction_rank` already orders the pair correctly
        //   (#1301) and there is no discrepancy to repair. `g.[264_265insCA;
        //   264_265dup]` and `g.[258dup;258_259dup]` are both well-ordered and
        //   must survive untouched.
        //
        // This is narrower than the equivalent test in
        // `demote_repeats_spanning_siblings` (#1287), which takes `>=` at the 5'
        // end. A repeat has no comparable sort/act split — it spans its tract
        // and is rendered there — so an interior junction is a problem for it
        // wherever it sits.
        let swallows_a_junction = siblings()
            .filter(|s| !s.claims_bases)
            .filter_map(|s| s.junction)
            .any(|junction| junction > dup.start && junction < dup.end);
        let collides = claims_the_bases || swallows_a_junction;
        if !collides {
            continue;
        }
        // The duplicated bases, read from the reference: `[start, end]` 1-based
        // inclusive is `[start - 1, end)` half-open 0-based. The span is already
        // on the served sequence's axis (#1284, #1482), so there is no per-region
        // conversion left to apply here.
        //
        // Checked because `dup.start` comes off a parsed description by way of a
        // conversion that adds the record's CDS bounds, so it is not bounded by
        // anything this function controls; an unchecked subtraction panics in a
        // debug build where the refusal one line below is the answer every other
        // unrepresentable coordinate here already gets.
        let Some(from_sequence) = dup.start.checked_sub(1) else {
            continue;
        };
        let (Ok(from), Ok(to)) = (u64::try_from(from_sequence), u64::try_from(dup.end)) else {
            continue;
        };
        let Ok(copied) = provider.get_sequence(&dup.provider_key, from, to) else {
            continue;
        };
        if copied.len() as i64 != dup.end - dup.start + 1 {
            continue;
        }
        let bases: Option<Vec<Base>> = copied.chars().map(Base::from_char).collect();
        let Some(bases) = bases else {
            continue;
        };
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::new(bases)),
        };
        // The copy lands at the junction 3' of the duplicated span, which is
        // the gap `[end, end + 1]`.
        //
        // The terminal boundary re-spelling (#1327) is refused when a
        // base-claiming sibling already occupies the base that gap runs past
        // (#1344). Re-spelling turns this member into a `delins` on that base,
        // so the pair goes from an adjacency the collapse pass merges to an
        // overlap it cannot; the zero-width junction form is what lets the
        // cancellation happen. The out-of-range coordinate is transient there —
        // the merge consumes it before anything is rendered.
        // `dup.end`, not `dup.end + 1`: the boundary identity deletes the *last*
        // base and re-inserts it carrying the payload, so the base the re-spelt
        // member claims is the one at the 5' side of the gap.
        //
        // Whether that sibling would *absorb* the junction form decides between
        // the remaining two answers (#1307) — a stricter test than "removes the
        // base", since a `delins` removes it and still does not absorb. Only an
        // absorbing sibling makes the out-of-range coordinate transient; against
        // one that does not — a substitution, a `delins`, an
        // inversion — the repair has to be abandoned instead.
        //
        // Quantified over **every** claimant, not the first one found: it takes
        // only one non-absorbing claimant to leave the junction form unconsumed,
        // so a single absorbing sibling is not licence to write it. Taking the
        // first match instead made the answer depend on member order, which is
        // authoring order and carries no meaning.
        //
        // No input reaches a two-claimant set today — two claimants on one base
        // is the overlap the parser rejects up front (`SelfCancellingAllele` for
        // a `dup`/`del` pair), and where the pair arrives by *shifting* instead
        // the cancellation runs before this pass. Measured by instrumenting this
        // site and running the suite: 1068 hits, 832 with one claimant and 236
        // with none, never more. So this is order-independence for its own sake
        // rather than a fix for an observable defect, and it is the safe
        // direction — it can only turn `WriteTransientJunction` into
        // `LeaveMemberUnchanged`, which is always in range.
        let mut landing_claimants = siblings()
            .filter(|s| s.claims_bases)
            .filter(|s| s.start <= dup.end && s.end >= dup.end)
            .peekable();
        let on_overrun = if landing_claimants.peek().is_none() {
            TerminalOverrun::RespellAtBoundary
        } else if landing_claimants.all(|s| s.absorbs_a_flush_insertion) {
            TerminalOverrun::WriteTransientJunction
        } else {
            TerminalOverrun::LeaveMemberUnchanged
        };
        respell_at_gap(&mut members[i], dup, dup.end, edit, provider, on_overrun);
    }
}

/// Coalesce cis members that settled onto **one junction** into a single
/// insertion carrying their combined payload.
///
/// A member that consumes no base occupies the interbase junction where its
/// added sequence lands. Two of them can shift onto the same junction from
/// distinct starting points — two one-base insertions at adjacent gaps inside a
/// homopolymer both travel to the tract's 3'-most junction:
///
/// ```text
/// reference  ("ACGT" x 64) + "AAAAAA" + ("ACGT" x 64)   poly-A run at 257-263
/// g.[258_259insA;259_260insA]  ->  g.[263dup;263dup]
/// ```
///
/// That output claims one position twice, so it has no well-defined resulting
/// sequence — `parse_hgvs` accepts it, and the SPDI apply oracle declines it
/// (#1286). Nothing upstream prevents it: the sibling clamps bound a member
/// against a sibling's *bases*, and neither of these has any.
///
/// It cannot be repaired downstream either. The sequence-first canonicalization
/// would derive the correct single member, but it runs *after* the per-member
/// pipeline and derives from the allele's resulting sequence — which an
/// overlapping allele does not have. It declines, so the corruption survives the
/// one pass that would have merged it.
///
/// So merge here, before that pass sees it. Both members add their payload at
/// the same interbase point, so one insertion carrying the concatenation denotes
/// exactly what the pair denoted. Re-spelling as an `Insertion` rather than a
/// wider `Duplication` is the same choice `respell_colliding_duplications`
/// makes: an insertion's payload is a literal and is correct wherever it lands,
/// whereas `g.262_263dup` would only be equivalent when `ref[262] == ref[263]`.
/// The next pass re-canonicalises it to a `dup` or a repeat where the reference
/// permits.
///
/// **Only when the colliding payloads commute** (`p ++ q == q ++ p`, see
/// [`payloads_commute`]). Concatenation is otherwise order-dependent, and the
/// order is exactly what the shift destroyed: two insertions at distinct
/// junctions had a defined order, and once they share a junction they do not.
/// Commuting payloads make the concatenation order-independent, so the repair is
/// unambiguous.
///
/// A non-commuting pair never reaches here: `clamp_sibling_crossing_junctions`
/// bars precisely that crossing (#1290), so the two stay at distinct junctions
/// in their original order. The two predicates are deliberately the same one.
///
/// Every axis, matching `respell_colliding_duplications`: the coordinate
/// conversion the transcript axes need is [`region_sequence_delta`], applied
/// in `junction_payload` where this pass reads bases (#1284).
pub(crate) fn coalesce_members_at_one_junction<P: ReferenceProvider>(
    before: &[HgvsVariant],
    members: &mut Vec<HgvsVariant>,
    phase: AllelePhase,
    uncertain: bool,
    provider: &P,
) {
    if phase != AllelePhase::Cis || uncertain || members.len() < 2 {
        return;
    }
    let Some(kind) = cis_kind_of(&members[0]) else {
        return;
    };
    let spans: Vec<Option<MemberSpan>> = members
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();

    // Group by (accession, region, junction). Only members that occupy a
    // junction and claim no bases participate.
    let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
    for i in 0..members.len() {
        let Some(span) = spans[i].as_ref() else {
            continue;
        };
        if span.claims_bases {
            continue;
        }
        let Some(junction) = span.junction else {
            continue;
        };
        match groups.iter_mut().find(|(first, _)| {
            spans[*first].as_ref().is_some_and(|s| {
                s.junction == Some(junction)
                    && s.region == span.region
                    && s.accession == span.accession
            })
        }) {
            Some((_, members_at)) => members_at.push(i),
            None => groups.push((i, vec![i])),
        }
    }

    // Where each member sat *before* this pass's per-member normalization, used
    // to order a group whose payloads do not commute. Aligned by index with
    // `members`: the caller builds both from one list and every pass between
    // mutates in place, so `before[i]` is member `i`'s origin. An empty or
    // mismatched snapshot disables the ordering rather than mis-attributing it.
    let origins: Vec<Option<(i64, i64)>> = if before.len() == members.len() {
        before
            .iter()
            .map(|v| member_span(v, kind, provider).map(|s| (s.start, s.end)))
            .collect()
    } else {
        vec![None; members.len()]
    };

    let mut remove: Vec<usize> = Vec::new();
    for (_, at) in &groups {
        if at.len() < 2 {
            continue;
        }
        // Concatenate in the order the members *came from*, not the order they
        // happen to sit in now. Two members can reach one junction without
        // either crossing the other — a collapse can create one there while a
        // sibling re-spells into it in place — and then the rendered order is
        // an artefact of their spans rather than of what the allele means
        // (#1323). Their origins are what the input actually said.
        let mut at: Vec<usize> = at.clone();
        let ordered_by_origin = at.iter().all(|&i| origins[i].is_some());
        if ordered_by_origin {
            at.sort_by_key(|&i| (origins[i], i));
        }
        let payloads: Option<Vec<Vec<Base>>> = at
            .iter()
            .map(|&i| junction_payload(&members[i], kind, spans[i].as_ref()?, provider))
            .collect();
        let Some(payloads) = payloads else {
            continue;
        };
        // Order-independence guard, needed only when the origins could not
        // supply an order. Commuting payloads (`p ++ q == q ++ p`) concatenate
        // the same either way, so this is exactly the case where "whatever
        // order they happen to be in" is safe — and it is the same predicate
        // `clamp_sibling_crossing_junctions` uses to decide whether the two were
        // allowed to meet here at all.
        //
        // With origins available the guard is unnecessary rather than merely
        // relaxed: a commuting group reaches the same concatenation under
        // either ordering, so ordering by origin subsumes the old behaviour
        // instead of overriding it.
        if !ordered_by_origin
            && payloads.iter().enumerate().any(|(x, left)| {
                payloads[x + 1..]
                    .iter()
                    .any(|right| !payloads_commute(left, right))
            })
        {
            continue;
        }
        let at = &at;
        let combined: Vec<Base> = payloads.concat();
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::new(combined)),
        };
        let target = at[0];
        let Some(span) = spans[target].as_ref() else {
            continue;
        };
        // The merged payload belongs at the junction the group shares, which is
        // not `span.end` when the target member is spelled as an insertion.
        let Some(junction) = span.junction else {
            continue;
        };
        let before = members[target].clone();
        respell_at_gap(
            &mut members[target],
            span,
            junction,
            edit,
            provider,
            TerminalOverrun::RespellAtBoundary,
        );
        if members[target] == before {
            continue; // the re-spell did not take; leave the group alone
        }
        remove.extend(at.iter().skip(1).copied());
    }
    if remove.is_empty() {
        return;
    }
    remove.sort_unstable();
    let mut index = 0;
    members.retain(|_| {
        let keep = !remove.contains(&index);
        index += 1;
        keep
    });
}

/// The payload an equivalent junction-inserting member carries at `to`, given
/// that it carries `payload` at `from`.
///
/// A member that adds sequence at a junction denotes those bases *at that
/// junction*, so the same payload spelled one position over denotes a different
/// sequence. The two agree only under a rotation, and only when the base being
/// stepped over is the one the rotation brings around:
///
/// ```text
/// insert P at junction j  ==  insert (ref[j] ++ P[..n-1]) at junction j - 1
///                             only when P[n-1] == ref[j]
/// ```
///
/// with the mirror image 3'. A multi-position move composes single steps, which
/// is why a payload can rotate through more than its own length.
///
/// `None` when some step is not payload-preserving, or the reference cannot be
/// read — the caller's signal to leave the member where it is rather than emit
/// a description of different bases.
///
/// `from` and `to` are **sequence** junctions, which is what
/// [`MemberSpan::junction`] carries (#1482): the bases stepped over are read
/// straight off the served sequence, with no per-region conversion left to
/// apply. Every caller passes junctions taken from a span, so they arrive
/// already converted.
fn payload_at_junction<P: ReferenceProvider>(
    payload: &[Base],
    from: i64,
    to: i64,
    key: &str,
    provider: &P,
) -> Option<Vec<Base>> {
    let base_at = |position: i64| -> Option<Base> {
        let (start, end) = (
            u64::try_from(position.checked_sub(1)?).ok()?,
            u64::try_from(position).ok()?,
        );
        let read = provider.get_sequence(key, start, end).ok()?;
        let mut chars = read.chars();
        let base = Base::from_char(chars.next()?)?;
        chars.next().is_none().then_some(base)
    };
    let mut rotated = payload.to_vec();
    let mut junction = from;
    while junction > to {
        // 5' by one: the payload's last base is the one it steps back over.
        if *rotated.last()? != base_at(junction)? {
            return None;
        }
        rotated.rotate_right(1);
        junction -= 1;
    }
    while junction < to {
        // 3' by one: the payload's first base is the one it steps over.
        if *rotated.first()? != base_at(junction + 1)? {
            return None;
        }
        rotated.rotate_left(1);
        junction += 1;
    }
    Some(rotated)
}

/// Move a junction-occupying member `delta` positions along its axis, keeping
/// the bases it adds unchanged.
///
/// [`translate_member`] alone is not enough for these members, because a
/// spelling can be position-dependent: a `Duplication` reads its payload from
/// the reference **under its own span**, so translating the span silently
/// rewrites what the member copies. On `g.[263_264insAC;265del]` that turned a
/// clamped `g.265_266dup` (payload `CA`) into `g.263_264dup` (payload `AA`) and
/// dropped the inserted `C` from the allele (#1292); the same rewrite corrupts
/// a repeat demoted to a duplication that the clamp then pulls back (#1280).
///
/// So the payload is rotated into phase at the destination first, and the
/// member is re-spelled as a plain insertion carrying it whenever the
/// translated spelling would denote something else. A move whose payload cannot
/// be rotated, read, or re-spelled leaves the member untouched: an unclamped
/// crossing is the defect this pass is repairing, but a corrupted payload is a
/// worse one.
fn translate_junction_member<P: ReferenceProvider>(
    variant: &mut HgvsVariant,
    delta: i64,
    kind: CisKind,
    from: &MemberSpan,
    provider: &P,
) {
    let (Some(junction), Some(payload)) = (
        from.junction,
        junction_payload(variant, kind, from, provider),
    ) else {
        return;
    };
    // Both junctions are sequence coordinates (#1482), so `delta` is a plain
    // count of bases and the rotation reads the reference directly.
    let destination = junction + delta;
    let Some(moved) = payload_at_junction(
        &payload,
        junction,
        destination,
        &from.provider_key,
        provider,
    ) else {
        return;
    };

    let original = variant.clone();
    translate_member(variant, delta, kind, from, provider);
    // A refused translation is **not** "nothing to repair" (#1513). It used to
    // return here, and that turned a correctly-computed barrier into a silent
    // no-op: the caller had already decided this member must not sit where it
    // does, and returning leaves it exactly there.
    //
    // `translate_member` refuses whenever the destination has no spelling *as
    // this edit kind* — and for a duplication near the 5' end of an axis with no
    // zero, that is routine rather than exotic. On `TAAAATTATATTTATTATTT`,
    // `c.[1_2insAA;2_3insTT]` shifts its insertion 3' to `c.4_5dup`, sweeping
    // across the `insTT` junction it does not commute with;
    // `clamp_sibling_crossing_junctions` computes the barrier correctly and asks
    // for junction 1, which for a two-base `dup` means the span `c.0_1`. There
    // is no `c.0`, so the write reverted, the member stayed at `4_5`, and the
    // allele went on to merge into `c.2_3insTTAA` — a description of **different
    // bases** than the input.
    //
    // The repair the caller wants exists: the member re-spelled as a plain
    // insertion carrying `moved` at `destination`, which is what the rest of this
    // function already does for a translation that landed but changed meaning.
    // A zero-width junction has a spelling everywhere a span may not, so falling
    // through reaches it. `respell_at_gap` reverts anything it cannot place, and
    // the guard at the end restores `original` when it does, so a member that
    // truly cannot be moved is left as it was — the old behaviour, now reached by
    // failing to place rather than by declining to try.
    let translation_refused = *variant == original;
    let landed = member_span(variant, kind, provider)
        .filter(|s| s.junction == Some(destination))
        .and_then(|s| junction_payload(variant, kind, &s, provider));
    // Only meaningful when the translation actually happened: on a refusal the
    // member is still `original`, so this reads the payload it carries *where it
    // was*, and if the barrier asked for the position it already occupies there
    // was nothing to do in the first place.
    if !translation_refused && landed.as_deref() == Some(moved.as_slice()) {
        return; // the spelling still denotes the same bases
    }
    let edit = NaEdit::Insertion {
        sequence: InsertedSequence::Literal(Sequence::new(moved)),
    };
    let translated = variant.clone();
    respell_at_gap(
        variant,
        from,
        destination,
        edit,
        provider,
        TerminalOverrun::RespellAtBoundary,
    );
    if *variant == translated {
        *variant = original;
    }
}

/// Drop a cis member asserting **no change** at a position another member
/// changes.
///
/// An identity member (`265=`) says a position was examined and found
/// unchanged. That is real information when it stands alone, so `g.[1002=;
/// 1005del]` keeps both. It is a contradiction when a sibling claims the same
/// bases, and the pair also overlaps — two members on one position, which
/// #1235's second criterion forbids and the apply oracle declines.
///
/// The shape is not authored, it is derived: a `del` and an `ins` that merge
/// into a `delins` restating the reference cancel to `=`, and the member that
/// absorbed the change then grows over it (#1297).
///
/// ```text
/// g.[261_262insAA;263del;264_265insA]  ->  g.[262_265A[6];265=]
///   the repeat alone denotes the input; the `265=` residue overlaps it
/// ```
///
/// Only the identity member is dropped, and only when a sibling's span *overlaps*
/// it under [`blocks_sibling_shift`] — a sibling that merely adds sequence at a
/// junction inside it contradicts nothing.
///
/// **Overlap, not containment** (#1416). The rule this enforces is "two members
/// must not claim one base", and a partial overlap breaks it exactly as a
/// containment does — the applier declines both alike:
///
/// ```text
/// g.[11_12inv;11_12insAA]  ->  g.[11_12=;10_11A[4]]  over `CGCGCGCGCAATCGCGCG`
///   the inversion of `AT` cancels to `=`; the repeat grows the `A` tract at
///   10_11 to four copies and so claims base 11, which the identity also names.
///   The repeat's span (10_11) does not *contain* the identity's (11_12), so a
///   containment test left the pair overlapping and the allele denoted nothing.
/// ```
///
/// Dropping is sequence-preserving even where the identity reaches past the
/// sibling: an identity asserts only that its bases are unchanged, so removing
/// it can never change what the allele denotes — the bases outside the sibling
/// are unchanged whether or not anything says so. What is lost is an assertion,
/// and keeping it costs a description no consumer can splice.
///
/// `blocks_sibling_shift` rather than `claims_reference_bases`, because a
/// duplication covers the positions it copies without consuming them (#1321):
///
/// ```text
/// g.[261_262insGA;262_263insA;263del]  ->  g.[262_263dup;263=]
///   `262_263dup` alone denotes the input; the `263=` overlaps it
/// ```
///
/// That is #1297's shape one spelling over, and the two predicates differ by
/// exactly the case that matters — `blocks_sibling_shift` is
/// `claims_reference_bases` plus `Duplication`/`DupIns`, which are barriers
/// precisely because they read their payload from the reference under their own
/// span. Naming those positions is what makes an identity member inside them
/// redundant.
///
/// Insertions stay out of it under either predicate, and must: an insertion's
/// span *is* the gap it occupies, so `g.264_265insA` nominally covers position
/// 264 while changing nothing there.
pub(crate) fn drop_identity_members_covered_by_siblings<P: ReferenceProvider>(
    members: &mut Vec<HgvsVariant>,
    phase: AllelePhase,
    uncertain: bool,
    provider: &P,
) {
    if phase != AllelePhase::Cis || uncertain || members.len() < 2 {
        return;
    }
    let Some(kind) = cis_kind_of(&members[0]) else {
        return;
    };
    let spans: Vec<Option<MemberSpan>> = members
        .iter()
        .map(|v| member_span(v, kind, provider))
        .collect();
    let is_identity = |v: &HgvsVariant| {
        matches!(
            cis_axis_parts(v, kind).map(|(_, _, _, _, edit)| edit),
            Some(NaEdit::Identity { .. })
        )
    };

    let covered: Vec<bool> = (0..members.len())
        .map(|i| {
            let (Some(span), true) = (spans[i].as_ref(), is_identity(&members[i])) else {
                return false;
            };
            (0..members.len()).filter(|&j| j != i).any(|j| {
                spans[j].as_ref().is_some_and(|s| {
                    s.blocks_shift
                        && s.region == span.region
                        && s.accession == span.accession
                        // Both spans forward. A reversed `<high>_<low>` range is
                        // SVD-WG006's circular deletion/duplication form, and
                        // `cis_axis_parts` hands the endpoints over raw, so one
                        // can reach here with `start > end`. This pass has no
                        // wraparound semantics: the interval test below reads
                        // such a span as an ordinary one, and widening
                        // containment to overlap made that newly *matter* — a
                        // reversed sibling against a forward identity now
                        // satisfies the test where containment refused it.
                        // Declining is the honest answer, and it is also the
                        // conservative one: the identity is kept, which is the
                        // behaviour that predates this rule.
                        && s.start <= s.end
                        && span.start <= span.end
                        && s.start <= span.end
                        && s.end >= span.start
                })
            })
        })
        .collect();
    if !covered.iter().any(|&drop| drop) {
        return;
    }
    let mut index = 0;
    members.retain(|_| {
        let keep = !covered[index];
        index += 1;
        keep
    });
}

/// Whether two junction payloads may be reordered without changing what the
/// allele denotes.
///
/// Exactly when they commute: `p ++ q == q ++ p`. That identity holds if and
/// only if both are powers of one common word, which is precisely the case
/// where the order two insertions apply in is unobservable — `A` and `AA` in a
/// poly-A tract commute, `A` and `C` do not.
///
/// Used in two places that must agree. `clamp_sibling_crossing_junctions` bars a
/// junction from crossing a sibling's only when the payloads do *not* commute,
/// and `coalesce_members_at_one_junction` merges members sharing a junction only
/// when they *do* — so a pair is either kept apart and ordered, or allowed
/// together and merged, never left overlapping.
fn payloads_commute(left: &[Base], right: &[Base]) -> bool {
    left.iter().chain(right).eq(right.iter().chain(left))
}

/// The bases a junction-occupying member adds, or `None` when they cannot be
/// determined without guessing.
///
/// A `Duplication` copies the reference under its own span, so its payload is
/// read from the provider. An `Insertion` carries a literal. Every other shape —
/// an unspecified count, a bracketed range, a `DupIns` — has no single literal
/// payload to concatenate, so the caller declines to merge it.
fn junction_payload<P: ReferenceProvider>(
    variant: &HgvsVariant,
    kind: CisKind,
    span: &MemberSpan,
    provider: &P,
) -> Option<Vec<Base>> {
    match cis_axis_parts(variant, kind).map(|(_, _, _, _, edit)| edit)? {
        NaEdit::Insertion {
            sequence: InsertedSequence::Literal(literal),
        } => Some(literal.bases().to_vec()),
        NaEdit::Duplication { .. } => {
            // The span is already on the served sequence's axis (#1284, #1482),
            // so it is read directly. Checked for the same reason as the sibling
            // conversion in `respell_colliding_duplications`: the span reaches
            // here through an addition of the record's CDS bounds, so an
            // unchecked subtraction would panic in a debug build where every
            // other unrepresentable coordinate here declines.
            let (from, to) = (
                u64::try_from(span.start.checked_sub(1)?).ok()?,
                u64::try_from(span.end).ok()?,
            );
            let copied = provider.get_sequence(&span.provider_key, from, to).ok()?;
            if copied.len() as i64 != span.end - span.start + 1 {
                return None;
            }
            copied.chars().map(Base::from_char).collect()
        }
        _ => None,
    }
}

/// Re-spell a member whose landing junction is one past the sequence end as the
/// equivalent boundary `delins` on the last base (#1327).
///
/// `edit` must be the literal `Insertion` the repair passes build; anything else
/// has no payload to fold into the `delins` and is left alone. Reverts unless
/// the rewritten member reads back as a single position at the last base — the
/// same all-or-nothing discipline [`respell_at_gap`] applies to its own writes.
fn respell_at_sequence_end<P: ReferenceProvider>(
    variant: &mut HgvsVariant,
    span: &MemberSpan,
    edit: &NaEdit,
    delta: i64,
    length: u64,
    provider: &P,
) {
    let NaEdit::Insertion {
        sequence: InsertedSequence::Literal(payload),
    } = edit
    else {
        return;
    };
    // The last base, read as the one-byte slice `insertion_to_boundary_delins`
    // anchors at position 1.
    let Ok(last) = provider.get_sequence(&span.provider_key, length.saturating_sub(1), length)
    else {
        return;
    };
    let Some(delins) = crate::normalize::insertion_to_boundary_delins(
        last.as_bytes(),
        payload,
        1,
        crate::normalize::BoundarySide::ThreePrime,
    ) else {
        return;
    };
    // Back onto the member's own axis: the last base is `length - delta` there.
    // `raw_length` keeps the unsigned form for the CDS-relative conversion below.
    let raw_length = length;
    let Ok(length) = i64::try_from(length) else {
        return;
    };
    let axis = length - delta;
    let original = variant.clone();
    let placed = match variant {
        HgvsVariant::Genome(g) => {
            u64::try_from(axis).is_ok_and(|v| {
                set_endpoints(
                    &mut g.loc_edit.location,
                    |p: &mut GenomePos, v: u64| p.base = v,
                    v,
                    v,
                )
            }) && {
                g.loc_edit.edit = Mu::Certain(delins.clone());
                true
            }
        }
        HgvsVariant::Mt(m) => {
            u64::try_from(axis).is_ok_and(|v| {
                set_endpoints(
                    &mut m.loc_edit.location,
                    |p: &mut GenomePos, v: u64| p.base = v,
                    v,
                    v,
                )
            }) && {
                m.loc_edit.edit = Mu::Certain(delins.clone());
                true
            }
        }
        // The transcript axes reach this too, and must be handled rather than
        // refused. `c.*11` on a 35-base transcript *is* the last base, so a
        // duplication resting there lands its copy past the end exactly as the
        // mitochondrial case does — and refusing left #1284's collision
        // unrepaired and the output unparseable, which the re-parse oracle
        // caught immediately. (An earlier draft of this arm asserted no such
        // case had been measured; the suite falsified that.)
        HgvsVariant::Tx(t) => {
            axis != 0
                && set_endpoints(
                    &mut t.loc_edit.location,
                    |p: &mut TxPos, v: i64| p.base = v,
                    axis,
                    axis,
                )
                && {
                    t.loc_edit.edit = Mu::Certain(delins.clone());
                    true
                }
        }
        HgvsVariant::Cds(c) => {
            let Some((region, axis)) = cds_axis_end(span, provider, raw_length, Region::Cds) else {
                return;
            };
            let utr3 = region == Region::ThreePrimeUtr;
            axis != 0
                && set_endpoints(
                    &mut c.loc_edit.location,
                    |p: &mut CdsPos, v: i64| {
                        p.base = v;
                        p.utr3 = utr3;
                        p.offset = None;
                    },
                    axis,
                    axis,
                )
                && {
                    c.loc_edit.edit = Mu::Certain(delins.clone());
                    true
                }
        }
        HgvsVariant::Rna(r) => {
            let Some((region, axis)) = rna_axis_end(span, provider, raw_length) else {
                return;
            };
            let utr3 = region == Region::ThreePrimeUtr;
            axis != 0
                && set_endpoints(
                    &mut r.loc_edit.location,
                    |p: &mut RnaPos, v: i64| {
                        p.base = v;
                        p.utr3 = utr3;
                        p.offset = None;
                    },
                    axis,
                    axis,
                )
                && {
                    r.loc_edit.edit = Mu::Certain(delins.clone());
                    true
                }
        }
        _ => false,
    };
    // Compare against the region the arm actually wrote, not against
    // `span.region`. For `c.`/`r.` those differ exactly when the last base lies
    // outside the member's own region — which is the case this guard exists for,
    // so keying the check on `span.region` would make it vacuously true for the
    // one input it has to judge.
    //
    // `rna_axis_end` rather than `cds_axis_end` for the same reason the arm
    // above uses it: on a non-coding record `cds_axis_end` declines, which made
    // this verification reject a write the arm had just made correctly (#1453).
    // The two must be the same resolution or the guard judges the arm against a
    // convention the arm does not use.
    let written = match variant {
        HgvsVariant::Cds(_) => cds_axis_end(span, provider, raw_length, Region::Cds),
        HgvsVariant::Rna(_) => rna_axis_end(span, provider, raw_length),
        _ => Some((span.region, axis)),
    };
    let landed = placed && written.is_some_and(|w| member_endpoints(variant) == Some((w, w)));
    if !landed {
        *variant = original;
    }
}

/// Which `(region, axis position)` names the sequence's **last base** on a
/// CDS-relative axis.
///
/// `respell_at_sequence_end` otherwise reaches the last base as
/// `length - delta`, using the *member's own* region delta. That is right only
/// while the last base lies in the member's region. It does not have to: a
/// positive coordinate outside the CDS is classified as a body region, so on a
/// 35-base record with `cds_start = 13`, `cds_end = 24` the last base came out
/// as `c.23` — which converts back the same way, so `member_endpoints` accepted
/// it — where the record's own axis calls it `c.*11`. `c.23` names a coding
/// position in a CDS that has only 12 bases.
///
/// Resolving through the transcript instead is what `respell_at_gap`'s
/// `cds_relative_gap` already does for the junction case; this is the same
/// conversion for the terminal case.
fn cds_axis_end<P: ReferenceProvider>(
    span: &MemberSpan,
    provider: &P,
    length: u64,
    body: Region,
) -> Option<(Region, i64)> {
    let (cds_start, cds_end) = ordered_cds_bounds(&span.provider_key, provider)?;
    cds_axis_position(i64::try_from(length).ok()?, cds_start, cds_end, body)
}

/// [`cds_axis_end`] for the `r.` axis, with the non-coding fallback (#1453).
///
/// A record with no CDS has nothing to resolve through, and needs none: its
/// `r.` axis is the transcript's own, so the last base is `length` in the
/// member's own region — the same answer [`respell_at_sequence_end`]'s `Tx` arm
/// reaches directly. `cds_axis_end` declines that record, and it is called
/// twice in that function — once to place the write and once to verify it — so
/// the fallback has to be in one shared place or the two disagree.
fn rna_axis_end<P: ReferenceProvider>(
    span: &MemberSpan,
    provider: &P,
    length: u64,
) -> Option<(Region, i64)> {
    if rna_axis_is_transcript_relative(span, provider) {
        return Some((span.region, i64::try_from(length).ok()?));
    }
    cds_axis_end(span, provider, length, Region::Rna)
}

/// Whether `span` sits on an `r.` axis that numbers its **transcript** rather
/// than its CDS — i.e. a `Region::Rna` member on a record with no CDS at all
/// (#1453).
///
/// The two `r.` repair arms — [`respell_at_gap`] and [`respell_at_sequence_end`]
/// — resolve their coordinates through the CDS, via [`cds_relative_gap`] and
/// [`cds_axis_end`] respectively. That is right on a coding record and is what
/// names a junction on the CDS start `r.-1_1` rather than the non-existent
/// `r.0` (#1284). Both route through [`ordered_cds_bounds`], which requires
/// **both** bounds, so on a non-coding record they returned `None`, the write
/// was reverted, and *every* repair became a silent no-op on a non-coding `r.`
/// member.
///
/// What that cost, measured on `coalesce_members_at_one_junction`: it grouped
/// the two colliding members and read both payloads before the revert threw its
/// write away, so `NR_TEST.1:r.[9dup;10dup;11dup]` normalized to
/// `r.[9dup;11dup;11dup]`. Two members claiming one interbase point denote no
/// sequence, and no seam oracle sees it — the output is well-formed, so the
/// re-parse oracle accepts it, and every coordinate is in range, so the
/// in-bounds oracle does too. Over a 3,200-row corpus of nearby junction-
/// occupying cis members on a non-coding transcript, 849 outputs repeated a
/// member and 1,046 disagreed with the `n.` axis of the same record; both go to
/// 0 with this fallback in place.
///
/// # Why this predicate and not "`ordered_cds_bounds` declined"
///
/// It is deliberately the same fallback [`region_sequence_delta`] and
/// [`axis_frame`] already make, and deliberately no wider — the three have to
/// agree about which records and regions are transcript-relative, or a repair
/// reads its bases through one convention and names them under another. So:
///
/// * **`Region::Rna` only.** `region_sequence_delta` falls back for that region
///   alone; `FivePrimeUtr` (`r.-N`) and `ThreePrimeUtr` (`r.*N`) still require
///   CDS bounds, because without a CDS those regions do not exist and there is
///   no position for a fallback to name.
/// * **`cds_start` alone**, rather than any `ordered_cds_bounds` failure — and
///   rather than "both bounds absent", which is the narrower test this
///   predicate first carried. `region_sequence_delta` and `axis_frame` both key
///   on `cds_start` only, so requiring `cds_end` to be absent too made the
///   three *disagree* on a record annotated with a stop codon but no start
///   codon. `CdotTranscript::from_*` reaches that shape from real data: its CDS
///   guard is an **or** (`start_codon.is_some() || stop_codon.is_some()`) and
///   then assigns both fields straight through, so a 5'-incomplete record
///   (`cds_start_NF`, #972) arrives as `(None, Some(end))`. On those records the
///   repair declined and #1453's repeated member survived — verbatim
///   `r.[9dup;11dup;11dup]`, against `n.9_10insTAA` on the same transcript.
///   The other two failures must keep refusing, and still do, because both
///   carry a `cds_start` for [`ordered_cds_bounds`] to reject: **inverted**
///   bounds are a malformed record whose CDS-relative axis is incoherent (see
///   [`ordered_cds_bounds`]), and a record carrying `cds_start` without
///   `cds_end` has a genuinely CDS-relative axis whose 3'UTR simply cannot be
///   named.
/// * **The provider must resolve the record.** `is_ok_and` is false for one it
///   cannot, which keeps refusing: an unresolvable transcript gives no grounds
///   to claim its axis is transcript-relative.
fn rna_axis_is_transcript_relative<P: ReferenceProvider>(span: &MemberSpan, provider: &P) -> bool {
    span.region == Region::Rna
        && provider
            .get_transcript(&span.provider_key)
            .is_ok_and(|tx| tx.cds_start.is_none())
}

// ----------------------------------------------------------------------------
// In-bounds seam oracle (#1353)
// ----------------------------------------------------------------------------

/// A coordinate an output names that its own sequence does not have.
///
/// Carries what the panic message needs and nothing more; every field is read
/// by `Normalizer::assert_in_bounds`.
///
/// Gated like the rest of this chain: `assert_in_bounds` and its
/// `in_bounds_self_check_enabled` flag are `#[cfg(debug_assertions)]`, and so is
/// the `use merge::OutOfBoundsCoordinate` that imports this, so in release these
/// items had no caller at all and only added compile time and dead code to the
/// binary. `test` rides along so the unit tests below still build under
/// `cargo test --release`, where `cfg(test)` holds but `debug_assertions` does not.
#[cfg(any(debug_assertions, test))]
#[derive(Debug)]
pub(crate) struct OutOfBoundsCoordinate {
    /// The provider key the length was read under, not the accession as
    /// written — that is the thing the number was actually compared against.
    pub(crate) accession: String,
    /// The offending position, converted onto the served sequence's own 1-based
    /// axis, so `c.*12` on an 11-base 3'UTR reports the transcript position
    /// rather than `12`.
    pub(crate) position: i64,
    /// The sequence's length, i.e. the largest position that does exist.
    pub(crate) length: u64,
    /// The member as rendered, so a failure inside an allele names which member.
    pub(crate) member: String,
}

/// The first coordinate in `variant` naming a position past the end of its own
/// sequence, or `None` if every coordinate exists.
///
/// This is the predicate behind `FERRO_ASSERT_IN_BOUNDS` (#1353). Three
/// instances of the class had been found by hand — #1274, #1343 and #1307 — each
/// filed and fixed separately, which is the argument for asking the question
/// once at the seam.
///
/// # Why it needs the axis conversion rather than a raw comparison
///
/// A `c.` position may legitimately exceed the CDS length: `*n` counts into the
/// 3'UTR and `-n` into the 5'UTR, so `c.*11` on a 35-base transcript with a
/// 12-base CDS is in range while `12` compared against `24` says nothing. Every
/// endpoint is therefore put onto the **served sequence's** axis with
/// [`region_sequence_delta`] before it is compared, which is the same conversion
/// the repair passes use.
///
/// # What it deliberately does not check
///
/// - **A reversed range is not an error.** SVD-WG006 admits `<high>_<low>` for a
///   circular deletion or duplication (`NC_012920.1:m.16563_13del`), so the two
///   endpoints are each checked against `[1, length]` independently and their
///   order is never compared. That also makes the check identical on `m.`/`o.`
///   and on the linear axes, which is the honest answer: a circular sequence has
///   a last base like any other, and #1327 established that the wrapped
///   *insertion* spelling that would blur it is not valid HGVS anyway.
/// - **A junction insertion naming `end + 1` is fine when that base exists.**
///   The question is whether a position exists, not whether an edit touches the
///   last one — `g.23_24insC` on a 24-base contig is correct and must stay
///   silent, while `g.24_25insC` is the #1307 defect.
/// - **Intronic offsets, unknown, uncertain and special positions** (`pter`,
///   `qter`) yield no plain axis position, so they are skipped — see
///   [`readable_endpoints`], which reads each end independently for a reason
///   worth knowing.  An intronic offset is outside the transcript by
///   construction and has no length to be compared with.
/// - **Protein axes**, for the same reason: `member_endpoints` covers the
///   nucleotide axes only. The three known instances of this class are all
///   nucleotide, so widening it is a separate change with its own evidence.
/// - **An inserted range payload** (`g.10_11ins[20_30]`), which names positions
///   in the location of no member. Not covered; if a defect of that shape turns
///   up, this is where to extend.
/// - **A provider that cannot report a length** leaves the member unchecked,
///   rather than failing a run for a reference that simply does not know. The
///   oracle only ever fires on a *known* overrun.
#[cfg(any(debug_assertions, test))]
pub(crate) fn first_out_of_bounds_coordinate<P: ReferenceProvider>(
    variant: &HgvsVariant,
    provider: &P,
) -> Option<OutOfBoundsCoordinate> {
    let members: &[HgvsVariant] = match variant {
        HgvsVariant::Allele(allele) => &allele.variants,
        single => std::slice::from_ref(single),
    };
    members
        .iter()
        .find_map(|member| member_out_of_bounds_coordinate(member, provider))
}

/// [`first_out_of_bounds_coordinate`] for one member.
#[cfg(any(debug_assertions, test))]
fn member_out_of_bounds_coordinate<P: ReferenceProvider>(
    member: &HgvsVariant,
    provider: &P,
) -> Option<OutOfBoundsCoordinate> {
    let accession = member.accession()?.transcript_accession();
    let length = provider.get_sequence_length(&accession).ok()?;
    let ceiling = i64::try_from(length).ok()?;
    readable_endpoints(member)
        .into_iter()
        .flatten()
        .find_map(|(region, axis)| {
            let delta = region_sequence_delta(region, &accession, provider)?;
            let position = axis.checked_add(delta)?;
            // The floor is part of the same invariant: a position at or below
            // zero exists on no sequence either, and the repairs that produce
            // `g.0_1` (#1282) are the same family as the ones that produce
            // `g.24_25`.
            (position < 1 || position > ceiling).then(|| OutOfBoundsCoordinate {
                accession: accession.clone(),
                position,
                length,
                member: member.to_string(),
            })
        })
}

/// Each end of a member's interval that has a plain `(region, axis position)`,
/// independently of whether the other end does.
///
/// [`member_endpoints`] deliberately demands **both** ends, because its callers
/// place a two-ended gap and cannot act on half of one. This oracle wants the
/// opposite: an endpoint it can read is one it can check, and an endpoint it
/// cannot read is simply not its business.
///
/// The difference is load-bearing rather than stylistic, and a measured false
/// positive is what showed it. `NC_PASTEND.1:g.pter_5000del` on a 200-base contig
/// normalizes to `g.1_5000del` — `pter` resolves to `1`, and `5000` is carried
/// through from the input. Asking `member_endpoints` about that *input* yields
/// `None`, because `pter` is a special position: the authored overrun becomes
/// invisible, the "did normalization introduce this?" exemption never fires, and
/// the oracle blames normalization for a coordinate the caller wrote. Reading the
/// ends independently sees the input's `5000`, and W4004 `PositionPastEnd`
/// remains the mechanism that reports it.
#[cfg(any(debug_assertions, test))]
fn readable_endpoints(v: &HgvsVariant) -> [Option<(Region, i64)>; 2] {
    fn ends<T>(
        interval: &Interval<T>,
        read: impl Fn(&UncertainBoundary<T>) -> Option<(Region, i64)>,
    ) -> [Option<(Region, i64)>; 2] {
        [read(&interval.start), read(&interval.end)]
    }
    match v {
        HgvsVariant::Genome(g) => ends(&g.loc_edit.location, simple_genome_pos),
        HgvsVariant::Mt(m) => ends(&m.loc_edit.location, simple_genome_pos),
        // `o.` reads exactly like `m.`: both are circular genomic axes carrying a
        // plain `GenomePos` interval. Omitting it left every `o.` member with no
        // readable endpoint, so the oracle could not fire on one at all — while
        // the doc above claimed the check was "identical on `m.`/`o.`". This is
        // the branch that makes that true. Note the sibling `member_endpoints`
        // deliberately still omits `o.`: its callers are the repair passes, which
        // exclude circular members by design (see `respell_colliding_duplications`
        // and the note at the top of this module), and widening it would change
        // what those passes rewrite rather than what this oracle inspects.
        HgvsVariant::Circular(o) => ends(&o.loc_edit.location, simple_genome_pos),
        HgvsVariant::Cds(c) => ends(&c.loc_edit.location, simple_cds_pos),
        HgvsVariant::Tx(t) => ends(&t.loc_edit.location, simple_tx_pos),
        HgvsVariant::Rna(r) => ends(&r.loc_edit.location, simple_rna_pos),
        _ => [None, None],
    }
}

fn member_endpoints(v: &HgvsVariant) -> Option<((Region, i64), (Region, i64))> {
    fn ends<T>(
        interval: &Interval<T>,
        read: impl Fn(&UncertainBoundary<T>) -> Option<(Region, i64)>,
    ) -> Option<((Region, i64), (Region, i64))> {
        Some((read(&interval.start)?, read(&interval.end)?))
    }
    match v {
        HgvsVariant::Genome(g) => ends(&g.loc_edit.location, simple_genome_pos),
        HgvsVariant::Mt(m) => ends(&m.loc_edit.location, simple_genome_pos),
        HgvsVariant::Cds(c) => ends(&c.loc_edit.location, simple_cds_pos),
        HgvsVariant::Tx(t) => ends(&t.loc_edit.location, simple_tx_pos),
        HgvsVariant::Rna(r) => ends(&r.loc_edit.location, simple_rna_pos),
        _ => None,
    }
}

/// Replace a member with `edit` over the gap `[junction, junction + 1]`.
///
/// The junction is passed rather than read off `span`, because where it sits
/// depends on the spelling: `span.end` is the junction for a duplication, but
/// one position past it for an insertion, whose span is the gap itself. It is a
/// **sequence** coordinate, like every other position [`MemberSpan`] carries
/// (#1482); on the genomic and `n.` axes that is the axis number itself, and on
/// the CDS-relative ones `cds_relative_gap` below names it back.
///
/// Reverts unless the result reads back exactly as intended.
///
/// # The CDS-relative axes, and the boundary junction (#1284)
///
/// On `g.`/`m.`/`n.` the gap is `[junction, junction + 1]` on the axis itself.
/// On `c.`/`r.` that arithmetic is wrong at exactly one place, and it is a
/// place these repairs reach: a duplication resting at `c.-1` lands its copy at
/// the junction between `c.-1` and `c.1`, and `-1 + 1` is `c.0`, which does not
/// exist. Adding 1 on an axis with no zero is the off-by-one #1284 records.
///
/// So the two ends are resolved *through the transcript* —
/// [`region_sequence_delta`] out, [`cds_axis_position`] back — which names each
/// end in whichever region it actually falls in and yields `c.-1_1` here. That
/// is valid nomenclature (`consultation/SVD-WG001.md:37` writes
/// `NM_001849.3:c.-1_*1=`), and the same machinery names the `cds_end`
/// boundary `c.<last>_*1` without a second case.
///
/// The verification is [`member_endpoints`] rather than [`member_span`] for
/// that reason: the correct answer here is a range whose ends sit in two
/// regions, and what this has to confirm is that each end was *written* in the
/// region it belongs to. `member_span` answers in sequence coordinates (#1482),
/// which is the right question for comparing members and the wrong one for
/// checking a spelling.
/// What [`respell_at_gap`] does when the gap it is told to write runs one base
/// past the end of the sequence.
///
/// Three outcomes, not two, because the right answer depends on what else the
/// allele has on that last base — and each of the three was arrived at by a
/// separate measurement. The gate in `respell_at_gap` documents each in full;
/// in brief:
///
/// | variant | when | why |
/// | --- | --- | --- |
/// | `RespellAtBoundary` | nothing else claims the last base | the boundary `delins` denotes the same bases and is in range (#1327) |
/// | `WriteTransientJunction` | **every** sibling claiming the last base *absorbs a flush insertion* | the pair is the #999 adjacency the collapse pass merges, so the out-of-range coordinate is consumed before anything renders (#1344) |
/// | `LeaveMemberUnchanged` | **any** sibling claiming the last base does not absorb one | nothing merges that pair, so both other outcomes are wrong: the boundary `delins` would overlap the sibling and the junction form would reach the output (#1307) |
///
/// "Absorbs a flush insertion" ([`merges_with_a_flush_insertion`]) rather than
/// "removes the last base": a `delins` removes its bases and still does not
/// absorb, so it takes the `LeaveMemberUnchanged` branch. The looser wording
/// named a superset and was what the earlier cut of #1307 got wrong.
///
/// Only the overrun is affected. An in-range gap is placed identically under
/// all three.
#[derive(Clone, Copy, PartialEq, Eq)]
enum TerminalOverrun {
    RespellAtBoundary,
    WriteTransientJunction,
    LeaveMemberUnchanged,
}

fn respell_at_gap<P: ReferenceProvider>(
    variant: &mut HgvsVariant,
    span: &MemberSpan,
    junction: i64,
    edit: NaEdit,
    provider: &P,
    on_overrun: TerminalOverrun,
) {
    // A member in a region with no bases of its own — `n.-N` / `n.*N` — carries
    // a **virtual** sequence position (see `region_span_delta`), which is enough
    // to compare it against a sibling and not enough to write it back. The
    // genomic and `n.` arms below take `gap_start` as an axis number directly,
    // which is only true where the two coincide, so refuse here rather than let
    // a virtual coordinate reach them and re-spell the member several bases
    // over. `region_sequence_delta`, deliberately, not `region_span_delta`:
    // "has real bases" is exactly the question being asked.
    if region_sequence_delta(span.region, &span.provider_key, provider).is_none() {
        return;
    }
    let original = variant.clone();
    let (gap_start, gap_end) = (junction, junction + 1);

    // Where each end of the gap lands, in `(region, axis)` form. Genomic and
    // `n.` gaps stay in the member's own region by construction; a
    // CDS-relative one is resolved through the transcript.
    let cds_relative_gap = |body: Region| -> Option<((Region, i64), (Region, i64))> {
        let (cds_start, cds_end) = ordered_cds_bounds(&span.provider_key, provider)?;
        // `junction` is already a sequence coordinate (#1482); only the way back
        // onto the axis is left to do.
        let right = junction.checked_add(1)?;
        Some((
            cds_axis_position(junction, cds_start, cds_end, body)?,
            cds_axis_position(right, cds_start, cds_end, body)?,
        ))
    };

    // The gap's 3' end must name a base the sequence actually has (#1327).
    //
    // Nothing else checks it. `landed` below only confirms the member reads
    // back at the gap it was *told* to use, which an out-of-range gap does; and
    // `parse_hgvs` does not know the sequence's length, so `m.16569_16570insAA`
    // re-parses cleanly and the `FERRO_ASSERT_REPARSE` oracle stays silent. A
    // member that 3'-shifts onto the final base is exactly how the repair
    // reaches this: its landing junction is one past the end.
    //
    // The wrapped spelling is not available. SVD-WG006 authorises the reversed
    // `<high>_<low>` range for **deletions and duplications only** (line 33,
    // examples `NC_012920.1:m.16563_13del` and `J01749.1:o.4344_197dup`);
    // `insertion.md` is silent, so the general 5'→3' rule stands and ferro's
    // parser rejects `m.16569_1insA` (#129). Every re-spelling here produces an
    // insertion, so there is no legal wrapped form to emit.
    //
    // **Declining is not good enough either**, which measurement is what
    // showed: refusing here left `coalesce_members_at_one_junction`'s two
    // members unmerged, and they had each already been canonicalised to
    // `m.16569dup` — an allele claiming one position twice, which is the very
    // #1286 defect that pass exists to remove, and which re-parses so no oracle
    // sees it. Refusing traded an out-of-range spelling for an overlapping one.
    //
    // So re-spell at the boundary instead, via the identity every other clamp
    // in this crate already uses: inserting `A'` immediately 3' of the last
    // base *is* deleting that base and inserting `ref[last] ++ A'`
    // (`insertion_to_boundary_delins`, shared with #1170 / #387 / #1202 /
    // #1205 / #1217). That is in range, valid on a circular or linear
    // reference alike, and denotes exactly what the out-of-range insertion did.
    //
    // Only acts when the overrun is *known*: a provider that cannot report a
    // length leaves the check unevaluated rather than rewriting a repair it
    // cannot judge.
    //
    // The trigger is the last base *exactly* (`junction + delta == length`), not
    // "at or past it". A junction genuinely beyond the end is a different case:
    // re-spelling it at the last base would silently relocate an edit the author
    // placed elsewhere, which is a worse answer than declining. Only the
    // one-past-the-end overrun has a boundary identity that denotes the same
    // bases.
    //
    // Skipped when a sibling that **removes** the last base already claims it
    // (#1344). The boundary identity trades a zero-width junction insertion for
    // a `delins` that *claims* that base, and those two are not interchangeable
    // when something else is already there: an insertion flush against a
    // sibling's deleted base is the #999 adjacency the collapse pass merges,
    // while two members claiming one base is an overlap it cannot. Re-spelling
    // therefore converted a mergeable pair into an unmergeable one, and the
    // out-of-range coordinate it was avoiding never reached the output in that
    // case anyway — the merge consumed it. Measured: `c.[*10del;*11dup]` gave
    // `c.[*11del;*11delinsAA]`, which ferro's own strict mode rejects as
    // `OverlapConflictingEdits / W5002`, where the junction form settles on the
    // correct `c.*11=`.
    //
    // **A sibling that claims the last base without absorbing a flush insertion
    // is a third case** (#1307), and both of the other two answers are wrong for
    // it. "Without absorbing" rather than "without removing": a `delins` removes
    // its bases and still does not absorb, so it lands here too — that shape is
    // what showed "removes the base" to be the wrong predicate, and
    // `merges_with_a_flush_insertion` is the one that holds. There
    // is no deletion for the junction form to abut, so nothing merges the pair
    // and the out-of-range coordinate is not transient — it reaches the output,
    // where it re-parses and is a fixed point, so neither existing oracle sees
    // it. Measured on a 24-base contig: `g.[24dup;24C>G]` gave
    // `g.[24C>G;24_25insC]`, naming a position the contig does not have. The
    // boundary `delins` is no better, for the #1344 reason one paragraph up —
    // the substitution already claims that base, so it would overlap. What is
    // left is to leave the member as the duplication it was: the collision the
    // caller wanted repaired stays unrepaired, which is the pre-existing
    // spelling of the input rather than a new defect, and every coordinate it
    // names exists.
    //
    // `junction` is a sequence coordinate (#1482), so the overrun test is a
    // direct comparison with the length. `delta` is still read, because
    // `respell_at_sequence_end` names the last base back on the *member's* axis
    // and needs the offset to get there.
    if let Some(delta) = region_sequence_delta(span.region, &span.provider_key, provider) {
        if let Ok(length) = provider.get_sequence_length(&span.provider_key) {
            if u64::try_from(junction) == Ok(length) {
                match on_overrun {
                    TerminalOverrun::RespellAtBoundary => {
                        respell_at_sequence_end(variant, span, &edit, delta, length, provider);
                        return;
                    }
                    TerminalOverrun::LeaveMemberUnchanged => return,
                    // Falls through to the placement below, which writes the
                    // out-of-range junction form deliberately.
                    TerminalOverrun::WriteTransientJunction => {}
                }
            }
        }
    }

    // The unsigned genomic pair, `None` on a negative axis. Computed per-arm
    // rather than as an early return for the whole function: hoisting it was
    // the shape this function had while it was genomic-only, and it silently
    // refused *every* CDS-relative repair in the 5'UTR, where the axis value is
    // negative by definition — the arms below never ran.
    // Rejects a zero endpoint as well as a negative one. `u64::try_from(0)`
    // succeeds, so without the explicit test a junction at 0 would name `g.0_1` /
    // `m.0_1`, and neither axis has a position 0 — the same refusal the `Tx` arm
    // below makes for `n.0`, which this closure has to match rather than leave to
    // one axis.
    let unsigned_gap = || -> Option<(u64, u64)> {
        (gap_start != 0 && gap_end != 0)
            .then(|| Some((u64::try_from(gap_start).ok()?, u64::try_from(gap_end).ok()?)))?
    };
    let same_region_gap = ((span.region, gap_start), (span.region, gap_end));

    // Where an `r.` gap lands, which depends on whether the record has a CDS
    // (#1453).
    //
    // On a **coding** transcript the `r.` axis is CDS-relative, so the gap is
    // resolved through the transcript exactly as `c.`'s is — that is what names
    // the junction on the CDS start `r.-1_1` rather than the non-existent `r.0`
    // (#1284). A **non-coding** record has no CDS to resolve through and its
    // `r.` axis numbers the transcript directly, so the gap is
    // `[junction, junction + 1]` on the axis itself: the `n.` arm's arithmetic
    // below, including its refusal of a gap that would name position 0, since
    // this axis has no zero either.
    //
    // `cds_relative_gap` alone *refused* that second case, which reverted the
    // write and made every repair routed through here a silent no-op on a
    // non-coding `r.` member. See [`rna_axis_is_transcript_relative`] for what
    // that cost, and for why the fallback is keyed on that predicate rather
    // than on `cds_relative_gap` returning `None`.
    let rna_gap = || -> Option<((Region, i64), (Region, i64))> {
        if rna_axis_is_transcript_relative(span, provider) {
            return (gap_start != 0 && gap_end != 0).then_some(same_region_gap);
        }
        cds_relative_gap(Region::Rna)
    };

    let (placed, intended) = match variant {
        HgvsVariant::Genome(g) => (
            unsigned_gap().is_some_and(|(start, end)| {
                set_endpoints(
                    &mut g.loc_edit.location,
                    |p: &mut GenomePos, v: u64| p.base = v,
                    start,
                    end,
                )
            }) && {
                g.loc_edit.edit = Mu::Certain(edit.clone());
                true
            },
            same_region_gap,
        ),
        HgvsVariant::Mt(m) => (
            unsigned_gap().is_some_and(|(start, end)| {
                set_endpoints(
                    &mut m.loc_edit.location,
                    |p: &mut GenomePos, v: u64| p.base = v,
                    start,
                    end,
                )
            }) && {
                m.loc_edit.edit = Mu::Certain(edit.clone());
                true
            },
            same_region_gap,
        ),
        // `n.` is a signed axis with no zero, so a gap that would name `n.0` is
        // refused — the same guard `place_signed` applies for `respell_as`.
        // Unlike `c.`/`r.` there is nothing on the far side to name it as: the
        // `n.` regions either side of `0` are the transcript and the sequence
        // *outside* it, so a conversion would have no bases to offer.
        HgvsVariant::Tx(t) => (
            gap_start != 0
                && gap_end != 0
                && set_endpoints(
                    &mut t.loc_edit.location,
                    |p: &mut TxPos, v: i64| p.base = v,
                    gap_start,
                    gap_end,
                )
                && {
                    t.loc_edit.edit = Mu::Certain(edit.clone());
                    true
                },
            same_region_gap,
        ),
        HgvsVariant::Cds(c) => match cds_relative_gap(Region::Cds) {
            Some(intended) => (
                set_endpoints(
                    &mut c.loc_edit.location,
                    |p: &mut CdsPos, (region, base): (Region, i64)| {
                        p.base = base;
                        p.utr3 = region == Region::ThreePrimeUtr;
                        p.offset = None;
                    },
                    intended.0,
                    intended.1,
                ) && {
                    c.loc_edit.edit = Mu::Certain(edit.clone());
                    true
                },
                intended,
            ),
            None => (false, same_region_gap),
        },
        HgvsVariant::Rna(r) => match rna_gap() {
            Some(intended) => (
                set_endpoints(
                    &mut r.loc_edit.location,
                    |p: &mut RnaPos, (region, base): (Region, i64)| {
                        p.base = base;
                        p.utr3 = region == Region::ThreePrimeUtr;
                        p.offset = None;
                    },
                    intended.0,
                    intended.1,
                ) && {
                    r.loc_edit.edit = Mu::Certain(edit.clone());
                    true
                },
                intended,
            ),
            None => (false, same_region_gap),
        },
        _ => (false, same_region_gap),
    };
    let landed = placed && member_endpoints(variant) == Some(intended);
    if !landed {
        *variant = original;
    }
}

#[cfg(test)]
mod splitter_reproducer_corpus;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::hgvs::variant::Accession;
    use crate::reference::MockProvider;

    /// The two geometries #1461 turns on, pinned as numbers so the reasoning
    /// survives without a 457-base literal in the test.
    ///
    /// Real, measured from the shipped output of
    /// `NC_000013.10:g.100809575_100810031inv` and from the `AAGCTA -> TAGCTT`
    /// control in `issue_1040_inv_overrecognition_probes`.
    ///
    /// The reported geometry is rebuilt from its measured gap histogram rather
    /// than approximated, because the margin is what the test is for: #1461 is
    /// 67.2% changed, so a fixture that is denser than the real case would keep
    /// passing under a threshold that refuses the case the fix exists for.
    #[test]
    fn density_separates_the_inversion_from_the_over_recognition_control() {
        // #1461: 109 members over 457 nt, 307 columns changed. The measured gap
        // histogram is {1: 78, 2: 24, 3: 2, 4: 2, 5: 2} — 108 gaps holding 150
        // unchanged columns, so the span is 307 + 150 = 457 exactly as reported.
        // Only the two totals reach the rule, so the 307 changed columns are
        // spread as evenly over the members as the total allows.
        let gaps: Vec<usize> = [(1, 78), (2, 24), (3, 2), (4, 2), (5, 2)]
            .into_iter()
            .flat_map(|(width, count)| std::iter::repeat_n(width, count))
            .collect();
        let member_count = gaps.len() + 1;
        let mut cursor = 0usize;
        let mut inversion: Vec<Piece> = Vec::with_capacity(member_count);
        for member in 0..member_count {
            let len = 307 / member_count + usize::from(member < 307 % member_count);
            inversion.push(Piece {
                ref_start: cursor,
                ref_end: cursor + len,
                alt: vec![b'N'; len],
            });
            cursor += len + gaps.get(member).copied().unwrap_or(0);
        }

        let span = inversion[member_count - 1].ref_end - inversion[0].ref_start;
        let changed: usize = inversion.iter().map(|p| p.ref_end - p.ref_start).sum();
        assert_eq!(inversion.len(), 109, "the reported member count");
        assert_eq!(span, 457, "the reported span");
        assert_eq!(changed, 307, "the reported changed-column total");
        assert_eq!(
            gaps.iter().copied().max(),
            Some(5),
            "the reported widest separation"
        );
        assert!(
            changed_columns_dominate_the_span(&inversion),
            "a dense inversion must be admitted even though its widest gap is 5"
        );
        assert!(
            !every_separation_is_a_single_base(&inversion),
            "and it is admitted by density alone -- the single-base gate refuses it"
        );

        // The #1040 control, `AAGCTA -> TAGCTT`: 6 nt, 2 changed, one gap of 4.
        let control = vec![
            Piece {
                ref_start: 0,
                ref_end: 1,
                alt: vec![b'T'],
            },
            Piece {
                ref_start: 5,
                ref_end: 6,
                alt: vec![b'T'],
            },
        ];
        assert!(
            !changed_columns_dominate_the_span(&control),
            "two isolated edits inside a near-palindrome must stay individual"
        );

        // The ordering that makes a separation threshold impossible is asserted
        // structurally above: the case to admit fails `every_separation_is_a_
        // single_base` (its widest real gap is 5) while the case to refuse has a
        // gap of 4, so no threshold on separation admits one without the other.
    }

    /// The moments quoted in [`changed_columns_dominate_the_span`]'s doc must
    /// reproduce the z column of the table beneath them (#1530).
    ///
    /// That doc is the one place in this file that argues from a distribution
    /// rather than from a pinned string, and its own thesis is that a statistic
    /// stated more broadly than it holds gets cited later as evidence. It then
    /// did exactly that: `E[U] = n/4` and `sd(U) = sqrt(3n/8)` are the **even**-`n`
    /// case of `U = 2·Binomial(floor(n/2), 1/4)`, asserted for "**every** `n`" two
    /// lines after the comment notes that an odd `n` has a centre column that can
    /// never coincide. #1461's span is odd (457), so the displayed formulas give
    /// `z = +2.73` where the table says `+2.75`.
    ///
    /// Pinned executably rather than corrected in prose alone, because prose is
    /// what was wrong: the table was right the whole time, which is precisely
    /// why nothing caught the generalisation.
    #[test]
    fn the_unchanged_column_moments_reproduce_the_tabulated_z_scores() {
        /// `U = 2 · Binomial(m, 1/4)` where `m = floor(n/2)` is the number of
        /// mirror pairs — the centre column of an odd span is not one, and over
        /// an A/C/G/T alphabet can never coincide (a self-complementary IUPAC
        /// code such as `N` does, which is why the model is stated over A/C/G/T).
        fn moments(n: usize) -> (f64, f64) {
            let m = (n / 2) as f64;
            (m / 2.0, (3.0 * m / 4.0).sqrt())
        }

        // (span, unchanged columns, the z the doc's table states)
        for (n, unchanged, tabulated) in [
            (6usize, 4.0, 1.67), // `AAGCTA -> TAGCTT`, the #1040 control
            (8, 4.0, 1.15),      // `AATGCACA -> TGTGCATT` (#1517)
            (457, 150.0, 2.75),  // #1461, and the only odd span of the three
        ] {
            let (mean, sd) = moments(n);
            let z = (unchanged - mean) / sd;
            assert!(
                (z - tabulated).abs() < 0.005,
                "n = {n}: the stated moments give z = {z:.4}, but the table says {tabulated}"
            );
        }

        // The discriminating half. Restating the moments in `n` alone is the
        // mistake being fixed, so pin that the two forms genuinely disagree —
        // otherwise the simpler-looking `n/4` reads as an equivalent rewrite.
        assert_eq!(moments(8), (2.0, 3f64.sqrt()), "even n: m/2 == n/4");
        let (odd_mean, _) = moments(457);
        assert_eq!(odd_mean, 114.0);
        assert_ne!(
            odd_mean,
            457.0 / 4.0,
            "an odd span's centre column is unpaired, so `n/4` overstates the mean"
        );
    }

    // ------------------------------------------------------------------
    // Exon-junction guard on the derivation input (#1450)
    // ------------------------------------------------------------------

    /// The guard is asserted here, next to the code, rather than only through
    /// `normalize()`. An integration assertion cannot isolate it: the derivation
    /// returns `Some` only when its result *differs* from what the per-member
    /// pipeline already produced, so a same-exon pair — the case that must stay
    /// admitted — declines for an unrelated reason and would make an
    /// integration-level "still derived" test pass vacuously.
    ///
    /// `MockProvider::with_test_data`'s `NM_001234.1` has exons at transcript
    /// 1-15 / 16-30 / 31-44 and `cds_start = 5`, so the axis-to-transcript delta
    /// is 4 and the junctions sit at `c.11|c.12` and `c.26|c.27`.
    #[test]
    fn the_exon_junction_guard_fires_only_across_a_junction() {
        let provider = MockProvider::with_test_data();
        let frame = AxisFrame {
            delta: 4,
            reading_frame: true,
        };
        let crosses =
            |lo, hi| crosses_exon_junction(CisKind::Cds, "NM_001234.1", &provider, &frame, lo, hi);

        // Spans that cross a junction — the defect.
        assert!(
            crosses(10, 13),
            "c.10..c.13 crosses the exon-1/exon-2 junction"
        );
        assert!(
            crosses(11, 12),
            "the junction itself is crossed by its two flanking bases"
        );
        assert!(
            crosses(13, 30),
            "c.13..c.30 crosses the exon-2/exon-3 junction"
        );

        // Spans wholly inside one exon must stay admitted, or the guard would
        // remove the sequence-first path from the transcript axes altogether.
        assert!(!crosses(6, 11), "c.6..c.11 is inside exon 1");
        assert!(
            !crosses(12, 26),
            "c.12..c.26 is exactly exon 2, edge to edge"
        );
        assert!(!crosses(27, 34), "c.27..c.34 is inside exon 3");
        assert!(!crosses(20, 20), "a single position cannot cross anything");
    }

    /// The genomic axes have no exons, so the guard must never fire for them —
    /// the whole `g.`/`m.` corpus depends on the derivation continuing to run.
    #[test]
    fn the_exon_junction_guard_never_fires_on_a_genomic_axis() {
        let provider = MockProvider::with_test_data();
        let frame = AxisFrame {
            delta: 0,
            reading_frame: false,
        };
        for kind in [CisKind::Genome, CisKind::Mt] {
            assert!(
                !crosses_exon_junction(kind, "NM_001234.1", &provider, &frame, 1, 10_000),
                "a genomic axis has no exon junctions to cross"
            );
        }
    }

    // ------------------------------------------------------------------
    // In-bounds seam oracle (#1353)
    // ------------------------------------------------------------------

    /// A 24-base genomic contig, matching #1307's fixture length.
    fn bounds_genomic_provider() -> MockProvider {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_TEST.1", "TTTTTTTTTAATATATTTTAATAC");
        provider
    }

    /// A 35-base transcript whose CDS is `13..=24`, so the 3'UTR runs `*1`..`*11`
    /// and the 5'UTR `-1`..`-12`. The axis conversion is the whole point of the
    /// oracle, and only a transcript with real UTRs on both sides exercises it.
    fn bounds_transcript_provider() -> MockProvider {
        use crate::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("GCAAAGCGCGCGATGAAACCCTAAGGCATTTTTAA".to_string()),
            cds_start: Some(13),
            cds_end: Some(24),
            exons: vec![Exon::new(1, 1, 35)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            cds_start_incomplete: false,
            exon_cigars: Vec::new(),
            cached_introns: std::sync::OnceLock::new(),
        });
        provider
    }

    fn overrun(descriptor: &str, provider: &MockProvider) -> Option<OutOfBoundsCoordinate> {
        let variant = parse_hgvs(descriptor).expect("fixture must parse");
        first_out_of_bounds_coordinate(&variant, provider)
    }

    /// The three shapes the class was found in, spelled as descriptions rather
    /// than reached through a defect, so they keep testing the predicate after
    /// each producer is fixed.
    #[test]
    fn the_known_out_of_range_shapes_are_flagged() {
        let genomic = bounds_genomic_provider();
        // #1307 on a 24-base contig, and the same overrun inside an allele.
        for descriptor in [
            "NC_TEST.1:g.24_25insC",
            "NC_TEST.1:g.[24C>G;24_25insC]",
            // #1274's shape: a two-base identity running one past the end.
            "NC_TEST.1:g.24_25=",
            "NC_TEST.1:g.25del",
            "NC_TEST.1:g.30_40del",
        ] {
            let found = overrun(descriptor, &genomic)
                .unwrap_or_else(|| panic!("`{descriptor}` must be flagged as out of bounds"));
            assert_eq!(
                found.length, 24,
                "`{descriptor}` compared against {found:?}"
            );
            assert!(
                found.position > 24,
                "`{descriptor}` must report the offending position, got {found:?}"
            );
        }
        // #1343's shape on the transcript axis: `*12` does not exist when the
        // 3'UTR ends at `*11`.
        let transcript = bounds_transcript_provider();
        let found = overrun("NM_TEST.1:c.*11_*12insAA", &transcript)
            .expect("`c.*11_*12insAA` must be flagged");
        assert_eq!(
            (found.position, found.length),
            (36, 35),
            "the report must be on the transcript's own axis, not the `*` count: {found:?}"
        );
    }

    /// Everything legitimate must stay silent. These are the cases that make the
    /// check non-trivial, and each would be a false positive under a naive
    /// comparison.
    #[test]
    fn legitimate_coordinates_are_not_flagged() {
        let genomic = bounds_genomic_provider();
        for descriptor in [
            // A junction insertion whose right anchor IS the last base. The
            // question is whether a position exists, not whether an edit reaches
            // the end.
            "NC_TEST.1:g.23_24insC",
            "NC_TEST.1:g.24dup",
            "NC_TEST.1:g.24C>G",
            "NC_TEST.1:g.1_24del",
            "NC_TEST.1:g.[24dup;24C>G]",
            // An accession the provider does not serve: no length, no verdict.
            "NC_ABSENT.1:g.9999del",
        ] {
            assert!(
                overrun(descriptor, &genomic).is_none(),
                "`{descriptor}` names only positions that exist and must not be flagged"
            );
        }
        let transcript = bounds_transcript_provider();
        for descriptor in [
            // `*n` and `-n` legitimately exceed the CDS length — the reason the
            // conversion exists. `*11` is the transcript's last base.
            "NM_TEST.1:c.*11del",
            "NM_TEST.1:c.*10_*11insAA",
            "NM_TEST.1:c.-12del",
            "NM_TEST.1:c.13del",
            // An intronic offset has no plain axis position and is skipped.
            "NM_TEST.1:c.13+5del",
            // So is an uncertain boundary.
            "NM_TEST.1:c.(13_15)del",
        ] {
            assert!(
                overrun(descriptor, &transcript).is_none(),
                "`{descriptor}` must not be flagged"
            );
        }
    }

    // Stated-reference validator (#1352)
    // ------------------------------------------------------------------

    /// `stated_reference_bases_match` against a window that starts at position 1,
    /// so a member's stated bases are compared with `reference[s - 1..]`.
    fn stated_bases_agree(descriptor: &str, reference: &str) -> bool {
        let variant = parse_hgvs(descriptor).expect("fixture must parse");
        let members: Vec<HgvsVariant> = match variant {
            HgvsVariant::Allele(allele) => allele.variants.clone(),
            single => vec![single],
        };
        stated_reference_bases_match(&members, CisKind::Genome, reference.as_bytes(), 1)
    }

    /// An identity spelling its base asserts that base, so a wrong assertion must
    /// be caught (#1352).
    ///
    /// The channel was missing from the alternation while the function's own doc
    /// comment claimed the list was complete, so a wrongly-stated `g.3T=` on a
    /// reference `A` fell through to `_ => None` and was accepted. The validator's
    /// correctness rested instead on `collect_canonical_edits` refusing identity
    /// members — a guard in a different function, which is the coupling the doc
    /// comment warns about.
    #[test]
    fn a_wrongly_stated_identity_base_is_rejected() {
        //                     1234567
        let reference = "GGATTAC";
        // Position 3 is `A`.
        assert!(
            stated_bases_agree("NC_TEST.1:g.3A=", reference),
            "a correctly stated identity base must be accepted"
        );
        assert!(
            !stated_bases_agree("NC_TEST.1:g.3T=", reference),
            "a wrongly stated identity base must be rejected — before #1352 this \
             channel fell through to the catch-all and was accepted"
        );
    }

    /// An identity that states no base asserts nothing, so it must stay neutral.
    /// The arm keys off `sequence: Some(..)`, and a bare `g.3=` has `None`.
    #[test]
    fn a_bare_identity_asserts_nothing() {
        assert!(
            stated_bases_agree("NC_TEST.1:g.3=", "GGATTAC"),
            "a bare identity carries no assertion and must not be judged"
        );
    }

    /// The channels that were already covered must keep working — the new arm
    /// shares its binding with four of them, so a mis-edit there would be silent.
    #[test]
    fn the_pre_existing_stated_reference_channels_still_match() {
        //                     1234567
        let reference = "GGATTAC";
        for (descriptor, expected) in [
            ("NC_TEST.1:g.3A>G", true),
            ("NC_TEST.1:g.3T>G", false),
            ("NC_TEST.1:g.3_5delATT", true),
            ("NC_TEST.1:g.3_5delAAA", false),
            ("NC_TEST.1:g.3_5dupATT", true),
            ("NC_TEST.1:g.3_5dupAAA", false),
            ("NC_TEST.1:g.3_5delATTinsGG", true),
            ("NC_TEST.1:g.3_5delAAAinsGG", false),
        ] {
            assert_eq!(
                stated_bases_agree(descriptor, reference),
                expected,
                "`{descriptor}` against `{reference}`"
            );
        }
    }

    /// The 5' end of the same invariant. A position at or below zero exists on no
    /// sequence, and the repairs that reach `g.0_1` (#1282) are the same family as
    /// the ones that reach `g.24_25`.
    #[test]
    fn a_five_prime_overrun_is_flagged_too() {
        let transcript = bounds_transcript_provider();
        // `c.-13` is one before the transcript's first base (the 5'UTR is 12 long).
        let found = overrun("NM_TEST.1:c.-13del", &transcript).expect("`c.-13del` must be flagged");
        assert_eq!(
            found.position, 0,
            "expected the transcript-axis 0: {found:?}"
        );
    }

    /// SVD-WG006 admits the reversed `<high>_<low>` range for a circular deletion
    /// or duplication, so the endpoints are checked independently and their order
    /// is never compared. Both must still exist.
    #[test]
    fn a_reversed_circular_range_is_judged_only_on_its_endpoints() {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_012920.1", "TTTTTTTTTAATATATTTTAATAC");
        // In range, wrapping: legal, and the check must not read it as 24 > 3.
        assert!(
            overrun("NC_012920.1:m.24_3del", &provider).is_none(),
            "a wrapped range whose ends both exist must not be flagged"
        );
        // Out of range on the high end, wrapping or not.
        assert!(
            overrun("NC_012920.1:m.25_3del", &provider).is_some(),
            "a wrapped range naming a position past the end must still be flagged"
        );
    }

    /// The whole diagnostic payload, not just the two numeric fields the tests
    /// above read.
    ///
    /// `accession` and `member` exist so a failure inside an allele says *which*
    /// member overran and *what length it was compared against*; they are read
    /// only by `Normalizer::assert_in_bounds`, which is `#[cfg(debug_assertions)]`.
    /// That left them provably dead in a release build — `cargo clippy --release
    /// --features dev --all-targets -- -D warnings` failed with "fields
    /// `accession` and `member` are never read" — and, more to the point, unasserted
    /// anywhere. Reading them here fixes both.
    #[test]
    fn the_report_names_the_member_and_the_accession_it_compared() {
        let genomic = bounds_genomic_provider();
        let found = overrun("NC_TEST.1:g.[24C>G;24_25insC]", &genomic)
            .expect("the allele's second member overruns");
        assert_eq!(
            found.accession, "NC_TEST.1",
            "the report must name the provider key the length was read under"
        );
        assert_eq!(
            found.member, "NC_TEST.1:g.24_25insC",
            "the report must name the offending member, not the whole allele: {found:?}"
        );
        assert_eq!((found.position, found.length), (25, 24), "{found:?}");
    }

    /// The `o.` half of the case above. `o.` is a circular genomic axis like
    /// `m.`, and the predicate's doc claims the check is "identical on `m.`/`o.`"
    /// — but `readable_endpoints` had no `Circular` branch, so every `o.` member
    /// yielded no readable endpoint and the oracle could not fire on one at all.
    /// Measured before the fix: `o.25_3del` on a 24-base sequence returned
    /// `None`, i.e. an out-of-bounds `o.` output was invisible to
    /// `FERRO_ASSERT_IN_BOUNDS`.
    ///
    /// Both directions are pinned, because a branch that returned the endpoints
    /// in the wrong order would pass a one-sided test.
    #[test]
    fn a_reversed_circular_o_range_is_judged_only_on_its_endpoints() {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("J01749.1", "TTTTTTTTTAATATATTTTAATAC");
        // In range, wrapping: legal under SVD-WG006, and the check must not read
        // it as 24 > 3.
        assert!(
            overrun("J01749.1:o.24_3del", &provider).is_none(),
            "a wrapped `o.` range whose ends both exist must not be flagged"
        );
        // Past the end on the high side — the endpoint the wrap makes `start`.
        let found = overrun("J01749.1:o.25_3del", &provider)
            .expect("a wrapped `o.` range naming a position past the end must be flagged");
        assert_eq!((found.position, found.length), (25, 24), "{found:?}");
        // And on the *other* endpoint, so the branch is not reading only one end.
        let found = overrun("J01749.1:o.20_25del", &provider)
            .expect("an `o.` range whose second end is past the end must be flagged");
        assert_eq!((found.position, found.length), (25, 24), "{found:?}");
    }

    /// An authored overrun must be reported by the same predicate on the *input*,
    /// which is what lets the assertion exempt it instead of blaming
    /// normalization. Reading each end independently is load-bearing here: `pter`
    /// is unreadable, and demanding both ends made the authored `5000` invisible.
    #[test]
    fn an_authored_overrun_is_visible_through_a_special_position() {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_PASTEND.1", "A".repeat(200));
        let found = overrun("NC_PASTEND.1:g.pter_5000del", &provider)
            .expect("the readable end must be checked even though `pter` is not");
        assert_eq!((found.position, found.length), (5000, 200), "{found:?}");
    }

    /// A stated base running past the window is a refusal, not a panic. Shared
    /// with every channel, but the new arm is a fresh way to reach it.
    #[test]
    fn a_stated_identity_past_the_window_is_refused() {
        assert!(
            !stated_bases_agree("NC_TEST.1:g.99A=", "GGATTAC"),
            "an assertion outside the window must be refused rather than indexed"
        );
    }

    /// #1135: an insertion-shaped anchor with nothing left to insert (a
    /// cancelling del+ins after 3'-shifting) must not become an `Insertion`
    /// carrying an empty sequence — that is not valid HGVS and made
    /// `normalize_na_edit` divide by zero rotating the empty sequence.
    ///
    /// Pinned on `build_naedit` with the anchor state directly: reaching it
    /// through `merge_consecutive_edits` requires the normalizer to have
    /// shifted the pair together first, so a merge-level test exercises the
    /// `Delins` branch instead and would not cover this.
    /// A repeat count is the one description that names an unbounded payload
    /// in a bounded number of characters, so the lowering has to bound it
    /// itself. `A[4000000000]` sizes a four-gigabyte `Vec` through
    /// `with_capacity`, which aborts the process rather than unwinding — in
    /// release as well as debug — and `Normalizer::normalize` reaches here.
    ///
    /// Asserted on `lowered_repeat` directly. Driving it through
    /// `canonicalize_from_sequence` would prove the same thing by *surviving*,
    /// and a test whose evidence is "the process did not die" cannot tell a
    /// working bound from a lowering that never ran.
    ///
    /// The just-over-the-window count is asserted alongside the absurd one on
    /// purpose: it is the case that can be checked *against* an unbounded
    /// build without allocating anything dangerous, so it is what makes this
    /// test demonstrably sensitive to the bound rather than merely consistent
    /// with it.
    #[test]
    fn an_absurd_repeat_count_is_refused_rather_than_allocated() {
        let reference = b"CCCCAAAAAGGGG";
        // Sanity: the same tract with a sane count does lower, so the refusals
        // below are the bound and not a mis-specified fixture.
        assert!(matches!(
            lowered_repeat(5, 5, b"A", 7, reference, 1),
            Some(GEdit::Ins { .. })
        ));
        // The tract's own five copies are kept, so this asks for 95 bases
        // more than the widest window the derivation will ever fetch.
        let just_past_the_window = MAX_CANONICAL_WINDOW as u64 + 100;
        assert!(lowered_repeat(5, 5, b"A", just_past_the_window, reference, 1).is_none());
        assert!(lowered_repeat(5, 5, b"A", 4_000_000_000, reference, 1).is_none());
    }

    #[test]
    fn empty_insertion_anchor_builds_an_identity() {
        // Insertion-shaped anchor (`start == end + 1`) with an empty payload.
        let anchor = Anchor {
            region: Region::Genome,
            start: 1010,
            end: 1009,
            alt: Vec::new(),
            form: AnchorForm::Replacement,
        };
        let (interval, edit) = build_naedit(anchor, |_, p| GenomePos::new(p as u64));

        match edit {
            NaEdit::Identity { .. } => {}
            other => panic!("expected NaEdit::Identity, got {other:?}"),
        }
        assert_eq!(interval.to_string(), "1009_1010");
    }

    /// The non-empty sibling still builds a real insertion — the guard must not
    /// swallow a payload-carrying anchor.
    #[test]
    fn non_empty_insertion_anchor_still_builds_an_insertion() {
        let anchor = Anchor {
            region: Region::Genome,
            start: 1010,
            end: 1009,
            alt: vec![Base::A],
            form: AnchorForm::Replacement,
        };
        let (_, edit) = build_naedit(anchor, |_, p| GenomePos::new(p as u64));
        match edit {
            NaEdit::Insertion { .. } => {}
            other => panic!("expected NaEdit::Insertion, got {other:?}"),
        }
    }

    /// A `dup`'s span reaches back from an insertion point that consumes
    /// nothing, so it is the one form here that can name reference bases a
    /// *preceding* piece already claims. It must decline rather than emit
    /// `[10_12delinsX;9_12dup]`.
    ///
    /// Pinned on `duplication_anchor` directly. No input has been shown to reach
    /// this through `partition_block` — candidates collapse to a single piece
    /// under `trim_common_flanks` — and every rail that would otherwise catch it
    /// is blind to it: `GEdit::Dup` copies from the untouched window, so the
    /// round-trip guard passes by construction, and `hgvs_to_spdi` renders a
    /// `dup` as a zero-width insert, so the disjointness helpers in `tests/it`
    /// see no overlap. An end-to-end test is therefore not available to write;
    /// this is what stands in for one.
    #[test]
    fn a_duplication_reaching_into_a_preceding_piece_declines() {
        // Window `AAACAG`. The piece is a pure insertion of `CAG` at offset 6,
        // whose source span is [3, 6) — a genuine tandem duplication.
        let ref_bytes = b"AAACAG";
        let piece = Piece {
            ref_start: 6,
            ref_end: 6,
            alt: b"CAG".to_vec(),
        };
        assert!(
            is_tandem_duplication(&piece, ref_bytes),
            "fixture must be a tandem duplication, or the test proves nothing"
        );

        // Preceding piece ends at 3: source span [3, 6) is clear of it.
        let clear = duplication_anchor(&piece, Region::Genome, 1, ref_bytes, 3)
            .expect("a disjoint source span still builds a dup");
        assert_eq!(clear.form, AnchorForm::Duplication);
        assert_eq!((clear.start, clear.end), (4, 6));

        // Preceding piece ends at 4: source span [3, 6) reaches one base into it.
        assert!(
            duplication_anchor(&piece, Region::Genome, 1, ref_bytes, 4).is_none(),
            "a dup whose source span overlaps the preceding piece must decline",
        );
    }

    /// The fallback the decline above lands on is not a refusal of the whole
    /// group — `anchor_for_piece` still returns the piece, spelled as the plain
    /// insertion it already is, which denotes the same bases.
    #[test]
    fn a_declined_duplication_falls_back_to_the_insertion_spelling() {
        let ref_bytes = b"AAACAG";
        let piece = Piece {
            ref_start: 6,
            ref_end: 6,
            alt: b"CAG".to_vec(),
        };
        let anchor = anchor_for_piece(
            &piece,
            Region::Genome,
            1,
            ref_bytes,
            SequenceEnds::INTERIOR,
            4,
            None,
        )
        .expect("the piece is still built");
        assert_eq!(anchor.form, AnchorForm::Replacement);
        // Insertion shape: `start == end + 1`, payload intact.
        assert_eq!((anchor.start, anchor.end), (7, 6));
        assert_eq!(anchor.alt, vec![Base::C, Base::A, Base::G]);
    }

    /// The terminal base of a CDS-relative record is named by the record's own
    /// axis, not by the member's region delta.
    ///
    /// `respell_at_sequence_end` used to reach the last base as
    /// `length - delta`, taking `delta` from the *member's* region. On a 35-base
    /// record with `cds_start = 13`, `cds_end = 24` a member classified as CDS
    /// body has `delta = cds_start - 1 = 12`, so the last base came out as
    /// `c.23` — and `member_endpoints` converts `c.23` back the same way, so the
    /// self-consistency check accepted it. `c.23` names a coding position in a
    /// CDS that has only 12 bases; the record's own axis calls transcript 35
    /// `c.*11`.
    ///
    /// Asserted on the conversion rather than end-to-end, deliberately: no
    /// *reachable* input distinguishes the two today, because positive
    /// out-of-CDS coordinates are re-spelled to `*N` upstream (`c.23dup`
    /// normalizes to `c.*11dup`), so `span.region` is already `ThreePrimeUtr`
    /// by the time the helper runs. Every candidate input was measured
    /// byte-identical with and without the fix. The guard is kept as hardening —
    /// it matches what `respell_at_gap`'s `cds_relative_gap` already does for
    /// the junction case — and pinned here rather than dressed up as an
    /// end-to-end regression it is not.
    #[test]
    fn the_terminal_base_is_named_by_the_records_own_axis() {
        // Transcript 35 on a record whose CDS is transcript 13..24.
        assert_eq!(
            cds_axis_position(35, 13, 24, Region::Cds),
            Some((Region::ThreePrimeUtr, 11)),
            "transcript 35 is `c.*11`, not the `c.23` that `length - delta` gave"
        );
        // The `r.` axis names the body differently and the UTRs identically.
        assert_eq!(
            cds_axis_position(35, 13, 24, Region::Rna),
            Some((Region::ThreePrimeUtr, 11))
        );
        // Inside the CDS the two agree, which is why the defect only showed at
        // the ends: transcript 20 is `c.8`, and `20 - 12` is also 8.
        assert_eq!(
            cds_axis_position(20, 13, 24, Region::Cds),
            Some((Region::Cds, 8))
        );
        // And the 5' side, where the axis has no zero.
        assert_eq!(
            cds_axis_position(12, 13, 24, Region::Cds),
            Some((Region::FivePrimeUtr, -1))
        );
        // No transcript position 0 or below.
        assert_eq!(cds_axis_position(0, 13, 24, Region::Cds), None);
    }

    /// A member whose two ends lie in different regions has a span (#1482).
    ///
    /// This is the root cause, pinned at the reader rather than only through an
    /// end-to-end description. `member_span` used to route through
    /// `cis_axis_parts` -> `join_pos`, which refuses a cross-region interval on
    /// the (false) grounds that such a range "has no valid HGVS syntax". It
    /// returned `None`, and every sibling-awareness pass drops a `None` through
    /// `filter_map` rather than declining — so `c.12_*1del` was invisible **as a
    /// sibling** module-wide, and `respell_colliding_duplications` could not see
    /// the collision it exists to repair.
    ///
    /// Pinned as sequence coordinates because that is what makes the member
    /// comparable with siblings in other regions at all: `c.12_*1` is two bases,
    /// while its written endpoints differ by `1 - 12 = -11`.
    #[test]
    fn a_member_spanning_a_region_boundary_has_a_span() {
        // CDS is transcript 13..=24, so `c.12` is transcript 24 and `c.*1` is 25;
        // `c.-1` is 12 and `c.1` is 13.
        let provider = bounds_transcript_provider();
        let span = |descriptor: &str| {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            member_span(&variant, CisKind::Cds, &provider)
        };

        let across_the_stop = span("NM_TEST.1:c.12_*1del").expect("CDS/3'UTR span must be read");
        assert_eq!((across_the_stop.start, across_the_stop.end), (24, 25));
        assert_eq!(
            (across_the_stop.axis_start, across_the_stop.axis_end),
            (12, 1),
            "the written pair is kept for whoever has to render it back"
        );
        assert!(across_the_stop.crosses_regions);
        assert_eq!(across_the_stop.region, Region::Cds, "where the span begins");
        assert_eq!(
            across_the_stop.end - across_the_stop.start + 1,
            2,
            "two bases; the written endpoints would have said -11"
        );

        let across_the_start = span("NM_TEST.1:c.-1_1del").expect("5'UTR/CDS span must be read");
        assert_eq!((across_the_start.start, across_the_start.end), (12, 13));
        assert!(across_the_start.crosses_regions);
        assert_eq!(across_the_start.region, Region::FivePrimeUtr);

        // A member wholly inside one region is unaffected, and is the control
        // that says the conversion did not simply shift everything.
        let inside = span("NM_TEST.1:c.11_12del").expect("in-CDS span must be read");
        assert_eq!((inside.start, inside.end), (23, 24));
        assert_eq!((inside.axis_start, inside.axis_end), (11, 12));
        assert!(!inside.crosses_regions);
    }

    /// The collision the boundary-crossing member was invisible to.
    ///
    /// `respell_colliding_duplications` turns a duplication that claims a
    /// sibling's bases back into an insertion, "a form that claims nothing and so
    /// cannot collide" (#1320/#1323). With `c.12_*1del` reading as `None` it saw
    /// no sibling at all and left `c.[11_12dup;12_*1del]` standing — two members
    /// claiming transcript 24, which ferro's own `parse_hgvs` rejects as a
    /// self-cancelling allele.
    ///
    /// Asserted on the *predicate*, not on a normalized string, so it pins the
    /// step that was broken rather than any downstream spelling of the repair.
    #[test]
    fn a_duplication_collides_with_a_boundary_crossing_sibling() {
        let provider = bounds_transcript_provider();
        let span = |descriptor: &str| {
            let variant = parse_hgvs(descriptor).expect("fixture must parse");
            member_span(&variant, CisKind::Cds, &provider).expect("span must be read")
        };
        let dup = span("NM_TEST.1:c.11_12dup");
        let del = span("NM_TEST.1:c.12_*1del");
        assert!(del.claims_bases);
        assert!(
            del.start <= dup.end && del.end >= dup.start,
            "the deletion claims transcript {}..={} and the duplication {}..={}; \
             they meet at 24, which is the collision the repair pass looks for",
            del.start,
            del.end,
            dup.start,
            dup.end
        );
    }

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

    #[test]
    fn merge_skips_wraparound_mt_anchor() {
        // Wraparound m. del + an adjacent sub at position 2 — adjacency
        // arithmetic on the raw (16569, 1) endpoints would suggest a
        // strict merge (1+1 == 2). The defensive guard returns None for
        // wraparound endpoints so merge declines to attempt it.
        let v1 = parse_hgvs("NC_012920.1:m.16569_1del").unwrap();
        let v2 = parse_hgvs("NC_012920.1:m.2A>G").unwrap();
        let out = merge_consecutive_edits(
            vec![v1.clone(), v2.clone()],
            AllelePhase::Cis,
            &MockProvider::new(),
        );
        assert_eq!(out.len(), 2, "expected no merge across origin wraparound");
        assert_eq!(out[0], v1);
        assert_eq!(out[1], v2);
    }

    #[test]
    fn merge_skips_wraparound_mt_when_wraparound_arrives_second() {
        // Reverse ordering: a linear sub at position 1 followed by a wraparound
        // del whose start (16569) is the natural "next" position after end (1).
        // Exercises the `simple_range_for_loc_edit` guard explicitly (the prior
        // test exercises `anchor_from_loc_edit`).
        let v1 = parse_hgvs("NC_012920.1:m.1A>G").unwrap();
        let v2 = parse_hgvs("NC_012920.1:m.16569_1del").unwrap();
        let out = merge_consecutive_edits(
            vec![v1.clone(), v2.clone()],
            AllelePhase::Cis,
            &MockProvider::new(),
        );
        assert_eq!(
            out.len(),
            2,
            "expected no merge when wraparound arrives as next"
        );
        assert_eq!(out[0], v1);
        assert_eq!(out[1], v2);
    }

    /// `apply_edits_to_window` must not depend on the order its edits arrive in.
    ///
    /// A `Delins` writes its replacement into the `after` slot of its first
    /// base, which is the same slot an insertion at that gap uses. Without the
    /// emptiness guard the two silently concatenated when the insertion came
    /// first and were refused when the delins did — one input with two answers,
    /// the exact non-confluence the sequence-first pass exists to remove. A
    /// slot collision has no well-defined resulting sequence, so the only
    /// order-independent answer is to refuse both ways.
    #[test]
    fn applying_a_delins_and_a_colliding_insertion_refuses_either_order() {
        let reference = b"ACGTACGTAC";
        // Rebuilt per call: `GEdit` is deliberately not `Clone`, and adding a
        // derive for a test's convenience is the wrong trade.
        let delins = || GEdit::Delins {
            s: 5,
            e: 6,
            alt: vec![Base::T, Base::T],
        };
        let insertion = || GEdit::Ins {
            gap: 5,
            alt: vec![Base::A, Base::A, Base::A],
        };

        let insertion_first = apply_edits_to_window(&[insertion(), delins()], reference, 1);
        let delins_first = apply_edits_to_window(&[delins(), insertion()], reference, 1);

        assert_eq!(
            insertion_first, delins_first,
            "member order changed the resulting sequence"
        );
        assert_eq!(
            insertion_first, None,
            "two edits writing one `after` slot must be refused"
        );
    }

    /// The guard must not refuse an insertion that merely sits *next to* a
    /// `delins` rather than in its slot, which would turn a bounded correctness
    /// fix into a silent loss of canonicalization.
    #[test]
    fn a_delins_and_a_non_colliding_insertion_still_apply() {
        let reference = b"ACGTACGTAC";
        // `delins` over [5, 6] writes into `after[5]`; the insertion at gap 8
        // writes into `after[8]`. Disjoint slots, so both apply.
        let delins = || GEdit::Delins {
            s: 5,
            e: 6,
            alt: vec![Base::T, Base::T],
        };
        let insertion = || GEdit::Ins {
            gap: 8,
            alt: vec![Base::G],
        };

        let forward = apply_edits_to_window(&[delins(), insertion()], reference, 1);
        let reversed = apply_edits_to_window(&[insertion(), delins()], reference, 1);

        assert_eq!(forward, reversed, "member order changed the result");
        // `ACGT` | `TT` replacing [5, 6] | `GT` | `G` inserted after 8 | `AC`
        assert_eq!(forward.as_deref(), Some(&b"ACGTTTGTGAC"[..]));
    }

    /// A long net deletion's split is refused **on its merits**, not by a length
    /// bound (#1271).
    ///
    /// This test used to pin the opposite mechanism: net deletions were the one
    /// unguarded regime, so a second bound of 32 nt kept them whole while its
    /// siblings were free to split at 1024. Extending `separations_are_meaningful`
    /// to every length-changing regime retired that bound, and this now pins that
    /// the rule reaches the same answer for the right reason.
    ///
    /// The distinction is worth asserting rather than assuming, because both
    /// mechanisms return the identical single piece: `whole()` is the refusal path
    /// for either. So the first case below checks the *aligner's* unrefused output
    /// too — four pieces, which is precisely `delins.md:44-47`'s harm — and checks
    /// that `separations_are_meaningful` is what rejects it.
    ///
    /// Tested here rather than through `normalize` because the end-to-end path
    /// routes through `best_alignment`'s max-match search, which would make the
    /// *fixture* rather than the rule decide the outcome.
    #[test]
    fn partition_block_refuses_a_coincidental_net_deletion_split() {
        let reference: Vec<u8> = b"ACGT".repeat(13); // 52 nt, past the 32 nt bound
        assert_eq!(reference.len(), 52);

        // Net deletion, the shape of `delins.md:44-47`'s worked example.
        let deletion_result = b"TTCCTCGATGCCTG".to_vec(); // 14 nt
        assert!(deletion_result.len() < reference.len());

        // What the aligner proposes when nothing refuses it: a four-way split of
        // a single delins, the corruption that passage warns about.
        let columns = best_alignment(&reference, &deletion_result)
            .expect("the aligner must reach this block now that no bound stops it");
        let proposed = pieces_from_columns(&columns, &reference, &deletion_result);
        assert!(
            proposed.len() > 1,
            "fixture must actually tempt the aligner into a split, else the \
             refusal below proves nothing; got {proposed:?}"
        );
        // And the rule — not a length bound — is what refuses it.
        let net_change = deletion_result.len().abs_diff(reference.len());
        assert!(
            !separations_are_meaningful(&proposed, net_change),
            "the {} proposed pieces are separated by coincidence and must be \
             refused on their merits; got {proposed:?}",
            proposed.len()
        );

        let pieces = partition_block(&reference, &deletion_result);
        assert_eq!(
            pieces.len(),
            1,
            "a 52 -> 14 net deletion must stay one piece; got {pieces:?}"
        );
        assert_eq!(pieces[0].ref_start, 0);
        assert_eq!(pieces[0].ref_end, reference.len());

        // Equal length, same span: no gap to place, so the raise applies and an
        // unchanged interior base still separates two pieces.
        let mut equal_result: Vec<u8> = reference.iter().map(|b| flip_base(*b)).collect();
        equal_result[25] = reference[25];
        let pieces = partition_block(&reference, &equal_result);
        assert_eq!(
            pieces.len(),
            2,
            "an equal-length block past 32 nt must still split; got {pieces:?}"
        );

        // Net insertion of the same span: guarded by `separations_are_meaningful`,
        // so the raise applies here too — the bound is not what stops it.
        //
        // The fixture is a lone substitution near the 5' end plus one contiguous
        // 2 nt insertion near the 3' end, far enough apart that the separation
        // *is* meaningful. That matters: `whole()` returns exactly one piece on
        // every refusal path, so `!pieces.is_empty()` cannot tell "reached the
        // aligner" from "the bound refused it" — only a genuine multi-piece
        // partition can.
        let mut insertion_result = reference.clone();
        insertion_result[2] = flip_base(reference[2]);
        insertion_result.splice(45..45, b"GG".iter().copied());
        assert!(insertion_result.len() > reference.len());
        let pieces = partition_block(&reference, &insertion_result);
        assert!(
            pieces.len() >= 2,
            "a net insertion past 32 nt must reach the aligner and keep its \
             separated pieces; a single piece is the bound-refusal `whole()` \
             fallback. got {pieces:?}"
        );
        assert_eq!(
            pieces[0].ref_start, 2,
            "the first piece must be the lone 5' substitution; got {pieces:?}"
        );

        // And the case that proves the retired bound is really gone: a net
        // *deletion* past 32 nt whose separations ARE meaningful must still
        // split.
        //
        // Without this the suite cannot tell "the rule refused a coincidental
        // split" from "a 32 nt bound refused every net deletion" — both return
        // one piece via `whole()`, so the refusal case above is satisfied by the
        // very mechanism this change retires. Only a net deletion that splits
        // distinguishes them.
        //
        // Same shape as the insertion fixture, so the two regimes are compared
        // like for like: a lone substitution near the 5' end plus one contiguous
        // 2 nt deletion near the 3' end. `net_change` is 2, which is at or below
        // `MAX_SINGLE_BASE_SEPARATION_CHANGE`, so `MIN_PIECE_SEPARATION` applies
        // and the wide unchanged run between them clears it easily.
        let mut deletion_split_result = reference.clone();
        deletion_split_result[2] = flip_base(reference[2]);
        deletion_split_result.drain(45..47);
        assert!(deletion_split_result.len() < reference.len());
        assert_eq!(reference.len() - deletion_split_result.len(), 2);
        let pieces = partition_block(&reference, &deletion_split_result);
        assert!(
            pieces.len() >= 2,
            "a net deletion past 32 nt with meaningful separations must reach \
             the aligner and keep its pieces; a single piece here would mean a \
             length bound is still refusing net deletions wholesale. got {pieces:?}"
        );
        assert_eq!(
            pieces[0].ref_start, 2,
            "the first piece must be the lone 5' substitution; got {pieces:?}"
        );
    }

    /// An insertion occupies no reference base, so it contributes no width to
    /// the separation between its neighbours.
    ///
    /// This is the arithmetic every corpus census of "how far apart are two
    /// members" has to get right, and getting it wrong is not visible in the
    /// output — it is visible only as a wrong *number*. Modelling an
    /// insertion's `A_B` anchor as a consumed two-base span once produced a
    /// reported "121 pairs at separation 0 — a spec violation", where the true
    /// count over 5.76M rows is zero, and it shifted the whole gap distribution
    /// with it. See `separations_are_meaningful`'s doc comment.
    ///
    /// Asserted against the function rather than through `normalize` on
    /// purpose: end-to-end the answer would be decided by `best_alignment`'s
    /// choice of gap placement, which is what the fixture would really be
    /// pinning.
    #[test]
    fn an_insertion_contributes_no_width_to_the_separation_it_sits_in() {
        let ins = |at: usize, alt: &[u8]| Piece {
            ref_start: at,
            ref_end: at,
            alt: alt.to_vec(),
        };
        let sub = |at: usize, alt: &[u8]| Piece {
            ref_start: at,
            ref_end: at + 1,
            alt: alt.to_vec(),
        };

        // Two insertions at reference offsets 10 and 11: exactly one unchanged
        // base (offset 10) lies between them, so `next.ref_start - prev.ref_end`
        // is 1 and `MIN_PIECE_SEPARATION` is met. Were the anchors read as
        // two-base spans the same pair would compute as separation 0 and be
        // refused.
        assert!(
            separations_are_meaningful(&[ins(10, b"A"), ins(11, b"T")], 2),
            "one unchanged base separates insertions at 10 and 11"
        );
        // Two insertions at the *same* junction genuinely have nothing between
        // them, and that is the only way a pair of insertions reaches 0.
        assert!(
            !separations_are_meaningful(&[ins(10, b"A"), ins(10, b"T")], 2),
            "insertions at one junction are not separated at all"
        );
        // Mixed pair: a substitution consuming offset 10 leaves offsets 11..12
        // unchanged before the insertion at 12.
        assert!(
            separations_are_meaningful(&[sub(10, b"A"), ins(13, b"T")], 1),
            "two unchanged bases separate a substitution at 10 from an \
             insertion at 13"
        );
        assert!(
            !separations_are_meaningful(&[sub(10, b"A"), ins(11, b"T")], 1),
            "an insertion at the junction immediately 3' of a substitution is \
             adjacent to it, not separated"
        );

        // The two cases above both clear `MIN_PIECE_SEPARATION` (1), so they
        // pin the junction-has-no-width rule without ever distinguishing the
        // two thresholds. Raise `net_change` past
        // `MAX_SINGLE_BASE_SEPARATION_CHANGE` so `RAISED_PIECE_SEPARATION` (2)
        // applies, and the one-base gap that passed above now fails — which is
        // what makes the insertion's zero width load-bearing rather than
        // incidental.
        assert!(
            !separations_are_meaningful(
                &[sub(10, b"A"), ins(12, b"T")],
                MAX_SINGLE_BASE_SEPARATION_CHANGE + 1
            ),
            "one unchanged base does not satisfy the raised separation"
        );
        assert!(
            separations_are_meaningful(
                &[sub(10, b"A"), ins(13, b"T")],
                MAX_SINGLE_BASE_SEPARATION_CHANGE + 1
            ),
            "two unchanged bases do satisfy the raised separation"
        );
    }

    /// Flank absorption must not wrap the submitter's copy number.
    ///
    /// `count` arrives as the stated count times `copies_per_stated_unit`, so a
    /// description written at the top of the range reaches the absorption step
    /// already saturated and `count += flank` wraps. Every other arithmetic step
    /// in [`lowered_repeat`] is checked; this pins the one that was not.
    ///
    /// The anchor spans an explicit two-position range (`e > s`) inside a
    /// five-`A` tract, which is what makes the 3' flank non-zero and the add
    /// reachable at all.
    #[test]
    fn lowered_repeat_refuses_a_count_that_would_overflow_flank_absorption() {
        assert!(
            lowered_repeat(5, 6, b"A", u64::MAX, b"CCCCAAAAAGGGG", 1).is_none(),
            "a count at the top of the range must be refused, not wrapped"
        );
        // Control: the same tract and anchor with a sane count still lowers, so
        // the guard above is rejecting the overflow rather than the shape.
        assert!(
            lowered_repeat(5, 6, b"A", 7, b"CCCCAAAAAGGGG", 1).is_some(),
            "the same shape with a small count must still lower"
        );
    }

    /// The two splitters must agree wherever the old one is already right.
    /// A single changed base has one unambiguous answer, so any disagreement
    /// here is a bug in the new path rather than a behavioural change.
    #[test]
    fn sequence_first_split_agrees_on_an_unambiguous_single_change() {
        let reference = b"GACTGACTGA";
        let result = b"GACTAACTGA";
        let old = partition_block(reference, result);
        let new = partition_block_sequence_first(
            reference,
            result,
            crate::normalize::seqfirst::MIN_SEPARATION,
        )
        .expect("an unambiguous single substitution must partition");
        assert_eq!(new.len(), 1, "expected one member, got {new:?}");
        assert_eq!(
            (new[0].ref_start, new[0].ref_end, new[0].alt.clone()),
            (old[0].ref_start, old[0].ref_end, old[0].alt.clone())
        );
    }

    /// The spec's own worked example (`delins.md:44`, `LRG_199t1:c.850_901`)
    /// must stay one member. This is the case that discriminates counting
    /// unchanged bases from comparing event offsets.
    #[test]
    fn sequence_first_split_keeps_the_spec_delins_example_whole() {
        let reference = b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT";
        let result = b"TTCCTCGATGCCTG";
        let pieces = partition_block_sequence_first(
            reference,
            result,
            crate::normalize::seqfirst::MIN_SEPARATION,
        )
        .expect("the spec example must partition");
        assert_eq!(pieces.len(), 1, "spec requires one delins, got {pieces:?}");
        assert_eq!((pieces[0].ref_start, pieces[0].ref_end), (0, 52));
    }

    /// Task 3, class `B2`: the separation threshold is axis-dependent, not a
    /// uniform constant. `general.md:34`'s general rule — two variants
    /// separated by one or more nucleotides are described individually — governs
    /// a genomic axis, so `MIN_SEPARATION_NO_FRAME` (`1`) must keep two runs one
    /// unchanged base apart split. `general.md:35`'s exception — separated by one
    /// nucleotide, together affecting one amino acid — governs a reading-frame
    /// axis, so `seqfirst::MIN_SEPARATION` (`2`) must merge the very same block.
    ///
    /// `CAA -> TAG` is the block `issue_1235_transcript_axes::` `n.`
    /// (non-coding, hence `AxisFrame::reading_frame == false`)
    /// `noncoding_axis_converges_on_separated_changes` produces, taken from the
    /// `FERRO_SEQFIRST_SHADOW` audit's class `B2` — a test named for exactly the
    /// genomic-axis split this threshold must preserve. Its equal-length blocks
    /// have a single dominator match at the middle base (only one minimal
    /// alignment achieves the edit distance), so the split is a pure threshold
    /// question with no alignment ambiguity involved.
    ///
    /// Not chosen: `AAT -> CAA` and `CAATT -> TA`, two of `B2`'s other curated
    /// blocks. Both stay one member even at `MIN_SEPARATION_NO_FRAME` — their
    /// apparent single unchanged base is not a dominator match at all (a second,
    /// equally minimal alignment matches the same reference position to a
    /// *different* alternate offset, so `Dominators::matched` is empty; see
    /// `align.rs`'s `AT -> AAT` example). That is `B1`'s mechanism, not `B2`'s,
    /// even though the audit's live-gap-width heuristic filed them under `B2`.
    /// So this fix does not close every `B2` block — only the ones where a
    /// dominator genuinely exists and the threshold alone was suppressing it.
    #[test]
    fn sequence_first_split_axis_separation_matches_general_rule() {
        let reference = b"CAA";
        let result = b"TAG";
        let genomic = partition_block_sequence_first(reference, result, MIN_SEPARATION_NO_FRAME)
            .expect("must partition");
        assert_eq!(
            genomic.len(),
            2,
            "a genomic axis must split at the one unchanged base: {genomic:?}"
        );

        let coding = partition_block_sequence_first(
            reference,
            result,
            crate::normalize::seqfirst::MIN_SEPARATION,
        )
        .expect("must partition");
        assert_eq!(
            coding.len(),
            1,
            "a reading-frame axis must merge across the one unchanged base: {coding:?}"
        );
        assert_eq!(
            (coding[0].ref_start, coding[0].ref_end),
            (0, reference.len())
        );
    }

    /// Codon-precision follow-up to Task 3: the coding-axis merge must also
    /// check the codon, not just the distance. `ref=GCA alt=TCC` is the block
    /// `issue_275_codon_frame_extensions::no_codon_frame_rna_pair_straddles_codon_boundary`
    /// produces, from the `FERRO_SEQFIRST_SHADOW` audit's residual `B2` list —
    /// and it is the exact block
    /// `merge_consecutive_edits_tests::test_no_codon_frame_pair_straddles_codon_boundary`
    /// already pins for the live splitter, at the same real position (`c.3G>T`,
    /// `c.5A>C`: `(3-1)/3 = 0`, `(5-1)/3 = 1`, different codons).
    /// `min_separation = 2` alone merges it into one triplet;
    /// `split_codon_incompatible_triplets` must split it straight back apart.
    #[test]
    fn sequence_first_split_declines_the_codon_exception_across_a_boundary() {
        let reference = b"GCA";
        let result = b"TCC";
        let mut merged = partition_block_sequence_first(
            reference,
            result,
            crate::normalize::seqfirst::MIN_SEPARATION,
        )
        .expect("must partition");
        assert_eq!(
            merged.len(),
            1,
            "the distance rule alone merges this into one triplet: {merged:?}"
        );

        split_codon_incompatible_triplets(&mut merged, true, 3, reference);
        assert_eq!(
            merged.len(),
            2,
            "different codons must not merge, even one nucleotide apart: {merged:?}"
        );
        assert_eq!(
            (
                merged[0].ref_start,
                merged[0].ref_end,
                merged[0].alt.as_slice()
            ),
            (0, 1, b"T".as_slice())
        );
        assert_eq!(
            (
                merged[1].ref_start,
                merged[1].ref_end,
                merged[1].alt.as_slice()
            ),
            (2, 3, b"C".as_slice())
        );

        // No reading frame: never touched, regardless of codon arithmetic.
        let mut no_frame = partition_block_sequence_first(
            reference,
            result,
            crate::normalize::seqfirst::MIN_SEPARATION,
        )
        .expect("must partition");
        split_codon_incompatible_triplets(&mut no_frame, false, 3, reference);
        assert_eq!(
            no_frame.len(),
            1,
            "no reading frame means no split: {no_frame:?}"
        );
    }

    /// The other side of the same pin: two substitutions genuinely within one
    /// codon must still merge. `ref=CCC alt=GCA` at `w_lo=10` is
    /// `merge_consecutive_edits_tests::test_codon_frame_two_subs_one_codon`'s
    /// exact block (`c.10C>G`, `c.12C>A`: `(10-1)/3 = 3`, `(12-1)/3 = 3`, same
    /// codon) — the spec's `c.[145C>T;147C>G]` → `c.145_147delinsTGG` example,
    /// re-based. `split_codon_incompatible_triplets` must leave it merged.
    #[test]
    fn sequence_first_split_keeps_the_codon_exception_within_one_codon() {
        let reference = b"CCC";
        let result = b"GCA";
        let mut merged = partition_block_sequence_first(
            reference,
            result,
            crate::normalize::seqfirst::MIN_SEPARATION,
        )
        .expect("must partition");
        assert_eq!(merged.len(), 1, "must merge into one triplet: {merged:?}");

        split_codon_incompatible_triplets(&mut merged, true, 10, reference);
        assert_eq!(merged.len(), 1, "same codon must stay merged: {merged:?}");
        assert_eq!((merged[0].ref_start, merged[0].ref_end), (0, 3));
        assert_eq!(merged[0].alt, b"GCA");
    }

    /// The calibration this task must not break: `LRG_199t1:c.850_901`
    /// (`delins.md:44`) merges a 52-base member for a reason unrelated to the
    /// codon exception (its `partition_block` sibling reaches the aligner and is
    /// then refused by `separations_are_meaningful` — see #1271), so
    /// `split_codon_incompatible_triplets` must not touch it: the piece is 52
    /// nt wide, not the 3 nt the codon exception's shape requires.
    #[test]
    fn split_codon_incompatible_triplets_leaves_the_spec_delins_example_whole() {
        let reference = b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT";
        let result = b"TTCCTCGATGCCTG";
        let mut pieces = partition_block_sequence_first(
            reference,
            result,
            crate::normalize::seqfirst::MIN_SEPARATION,
        )
        .expect("the spec example must partition");
        assert_eq!(pieces.len(), 1);

        // The real accession's CDS offset (`c.850` = this block's reference
        // offset 0; `cds_start = 245`, `LRG_199t1` transcript offset 1094,
        // `spec_enumeration_transcripts.fna`), so this is the exact codon
        // arithmetic production would run, not a stand-in value.
        split_codon_incompatible_triplets(&mut pieces, true, 850, reference);
        assert_eq!(
            pieces.len(),
            1,
            "spec still requires one delins, got {pieces:?}"
        );
        assert_eq!((pieces[0].ref_start, pieces[0].ref_end), (0, 52));
    }

    /// Every piece the new splitter emits must reconstruct the result block
    /// exactly. A split that does not round-trip denotes a different variant.
    #[test]
    fn sequence_first_split_round_trips_every_block_it_accepts() {
        fn words(len: u32) -> Vec<Vec<u8>> {
            if len == 0 {
                return vec![Vec::new()];
            }
            let mut out = Vec::new();
            for shorter in words(len - 1) {
                for base in *b"AC" {
                    let mut w = shorter.clone();
                    w.push(base);
                    out.push(w);
                }
            }
            out
        }
        let all: Vec<Vec<u8>> = (0..=5u32).flat_map(words).collect();
        // Counted, because the `continue` below is silent: if a regression made
        // `member_alt_spans` decline every pair, the round-trip assertion would
        // never run and this test would still pass. The sibling sweep
        // `shrinking_before_and_after_the_three_prime_shift_agree` guards itself
        // the same way.
        let mut accepted = 0usize;
        for reference in &all {
            for result in &all {
                let Some(pieces) = partition_block_sequence_first(
                    reference,
                    result,
                    crate::normalize::seqfirst::MIN_SEPARATION,
                ) else {
                    continue;
                };
                accepted += 1;
                let mut rebuilt = Vec::new();
                let mut cursor = 0usize;
                for piece in &pieces {
                    rebuilt.extend_from_slice(&reference[cursor..piece.ref_start]);
                    rebuilt.extend_from_slice(&piece.alt);
                    cursor = piece.ref_end;
                }
                rebuilt.extend_from_slice(&reference[cursor..]);
                assert_eq!(
                    rebuilt,
                    *result,
                    "{} -> {} rebuilt as {}",
                    String::from_utf8_lossy(reference),
                    String::from_utf8_lossy(result),
                    String::from_utf8_lossy(&rebuilt)
                );
            }
        }
        assert!(
            accepted > 3_000,
            "the sweep accepted only {accepted} blocks; a mass decline would make \
             this round-trip assertion vacuous"
        );
    }

    /// An oversized grid is declined rather than allocated.
    ///
    /// The reference side is bounded upstream by `MAX_CANONICAL_WINDOW` (4096),
    /// but the alternate side is not: a duplication over the window doubles it,
    /// and `AlignmentDag::build` allocates four `Θ(n·m)` grids. The shape
    /// `4096 -> 8192` is 33.5 M cells, past the bound, and must come back `None`
    /// so the caller keeps the per-member result.
    ///
    /// Asserted on the boundary from both sides, so the test fails if the bound
    /// is removed *or* set low enough to refuse ordinary blocks.
    #[test]
    fn sequence_first_split_declines_an_oversized_grid() {
        // The budget is a square grid at the reference bound, so the largest
        // equal-length block is exactly admissible.
        assert_eq!(MAX_SEQFIRST_GRID_CELLS, 4097 * 4097);

        // Just inside: a 4096 x 4096 grid is the boundary square.
        let reference = b"ACGT".repeat(1024);
        let result = b"ACGA".repeat(1024);
        assert_eq!(reference.len(), 4096);
        assert!(
            partition_block_sequence_first(
                &reference,
                &result,
                crate::normalize::seqfirst::MIN_SEPARATION
            )
            .is_some(),
            "a 4096 x 4096 grid is within budget and must still be partitioned"
        );

        // Past it: the same reference against a duplicated alternate, the shape a
        // whole-window `dup` produces.
        let doubled = b"ACGT".repeat(2048);
        assert_eq!(doubled.len(), 8192);
        assert_eq!(
            partition_block_sequence_first(
                &reference,
                &doubled,
                crate::normalize::seqfirst::MIN_SEPARATION
            ),
            None,
            "a 4096 x 8192 grid is 33.5 M cells and must be declined, not allocated"
        );
    }

    /// Pieces must be disjoint and ascending — the property that makes the
    /// overlap and ordering repair passes unnecessary.
    #[test]
    fn sequence_first_split_emits_disjoint_ascending_pieces() {
        let reference = b"ACGTACGTACGTACGT";
        let result = b"ATGTACGTACGTACTT";
        let pieces = partition_block_sequence_first(
            reference,
            result,
            crate::normalize::seqfirst::MIN_SEPARATION,
        )
        .expect("must partition");
        for pair in pieces.windows(2) {
            assert!(
                pair[0].ref_end <= pair[1].ref_start,
                "pieces overlap or are out of order: {pieces:?}"
            );
        }
    }

    /// `coalesce_adjacent_pieces` must reject a strictly overlapping pair
    /// rather than merge it.
    ///
    /// Unreachable in production (the 3'-shuffle cannot produce it — see the
    /// function's doc comment), but the old `>=` condition's `max()`/`extend()`
    /// body would silently double-count the shared columns, changing the
    /// denoted sequence instead of merely mis-describing it.
    ///
    /// Gated on `debug_assertions` because what it asserts *is* a
    /// `debug_assert!`, which compiles out of a release build — the function
    /// then falls through to the `else` arm and leaves the pair un-coalesced
    /// (safe, just silent), so `#[should_panic]` would fail there for the right
    /// reason and the wrong outcome. The release behaviour is deliberately not
    /// asserted here: it is the fail-safe, not the contract.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "strict overlap")]
    fn coalesce_rejects_a_strictly_overlapping_pair() {
        let mut pieces = vec![
            Piece {
                ref_start: 0,
                ref_end: 4,
                alt: b"AAAA".to_vec(),
            },
            Piece {
                ref_start: 2,
                ref_end: 6,
                alt: b"CCCC".to_vec(),
            },
        ];
        coalesce_adjacent_pieces(&mut pieces);
    }

    /// The enumerated classes of disagreement between the two block splitters.
    ///
    /// Harvested with `FERRO_SEQFIRST_SHADOW=1` over the whole suite plus the
    /// exhaustive two-member sweeps: **30 261 distinct `(reference, alternate)`
    /// blocks** disagree. Every one of the 60 522 piece lists — both sides, every
    /// block — reproduces the alternate block when its payloads are substituted
    /// into the reference. So **neither splitter is ever wrong** in the only sense
    /// that admits a definitive answer; each disagreement is a difference in how
    /// the change is *divided*, and the tie-breaker below is minimality: the
    /// changed-column count `changed_columns_of_pieces` measures, which is the
    /// quantity HGVS asks a description to minimise.
    ///
    /// Counts are distinct blocks, and sum to 30 261.
    ///
    /// | class | mechanism | n | more minimal |
    /// |-------|-----------|---|--------------|
    /// | `A1` | shadow splits a block live kept whole | 7 035 | shadow |
    /// | `A2` | shadow replaces live's coincidental substitution run with one indel pair | 6 915 | shadow |
    /// | `A3` | shadow narrows a piece live over-widened | 4 567 | shadow |
    /// | `A4` | shadow widens one piece but lowers the total | 379 | shadow |
    /// | `B1` | shadow widens a member over an indel it will not localise | 8 044 | live |
    /// | `B2` | shadow merges pieces live kept separate | 3 088 | live |
    /// | `B3` | shadow narrows one piece and pays for it in the next | 1 | live |
    /// | `C1` | equal weight, different division | 232 | neither |
    ///
    /// `A1`–`A4` (18 896, 62.4%) are the migration's payoff. The mechanism, in the
    /// blocks pinned below, is `partition_block`'s single-gap restriction: two
    /// independent length changes have no one-gap explanation, so the lone gap
    /// parks at one end and the offset columns pair up as coincidental
    /// substitutions or as one whole-block `delins`.
    ///
    /// `B1`–`B3` (11 133, 36.8%) are the risk, and they have two distinct causes.
    /// `B1` is inherent to deriving the split from the alignment steps common to
    /// every minimal alignment: when the inserted bases can sit at more than one
    /// offset inside a run, no single alignment node is common to all of them, so
    /// the member absorbs the whole ambiguous span instead of committing to a
    /// tie-broken offset.
    ///
    /// `B2` splits in two by the narrowest unchanged run between live's pieces:
    /// 1 453 blocks are separated by exactly one unchanged base (1 106 of those
    /// equal-length, where `partition_block` compares position-wise and there is
    /// no alignment choice to make), and 1 635 by two or more. The one-base half
    /// is a threshold mismatch rather than an ambiguity —
    /// `seqfirst::MIN_SEPARATION` is 2 on every axis, whereas this module's own
    /// `MIN_PIECE_SEPARATION` of 2 gates only net insertions.
    ///
    /// `B2` is the class that bears on flipping the default, because it is what
    /// the curated issue corpora overwhelmingly reach: of the 49 distinct blocks
    /// the non-sweep, non-proptest tests produce, 38 are a `B2` merge. See
    /// `shadow_merges_across_one_unchanged_base_where_live_splits`, which pins two
    /// blocks taken straight from the #1232 and #1235 transcript-axis tests.
    ///
    /// This module's tests call `partition_block_sequence_first` directly at the
    /// uniform `seqfirst::MIN_SEPARATION` (`2`), on purpose: they document the
    /// raw algorithm-level disagreement this audit found, unchanged. The fix —
    /// choosing the separation per axis at the `canonicalize_from_sequence` call
    /// site rather than changing this default — is pinned separately, by
    /// `sequence_first_split_axis_separation_matches_general_rule`.
    /// The `FERRO_PARTITION` bake-off knob.
    ///
    /// The knob exists to make four partitioners measurable side by side. Its
    /// most important property is therefore the one nobody looks at: that leaving
    /// it alone changes nothing at all.
    mod partition_rule_knob {
        use super::*;

        /// The audit's `live=` baseline must be `partition_block`, not whatever
        /// `FERRO_PARTITION` selected.
        ///
        /// Deleted once already: this branch's rebase resolved a conflict in
        /// favour of the pre-fix side, dropping both the independent `live`
        /// binding and this guard, so the vacuity returned silently. Under
        /// `FERRO_SEQFIRST_SHADOW=1 FERRO_PARTITION=shadow` the audit then
        /// compared the sequence-first splitter with itself and logged
        /// `SEQFIRST_AGREE` on every block.
        ///
        /// What this test can and cannot do: `partition_rule()` caches in a
        /// `OnceLock`, so one process cannot exercise two rules and the audit
        /// itself is unreachable from here. What it *can* pin is the premise the
        /// fix rests on — that live and canonical genuinely disagree on some
        /// block, so which one the audit names is observable rather than moot.
        /// `A -> ACA` is the smallest such block: both spellings denote `ACA`,
        /// but live anchors the inserted `AC` at reference 0 while canonical
        /// takes the 3'-most optimal step and anchors `CA` at reference 1.
        /// The end-to-end guard needs an isolated process per rule (#1444
        /// follow-up).
        #[test]
        fn live_and_canonical_disagree_so_the_audit_baseline_is_observable() {
            let live = partition_block(b"A", b"ACA");
            let canonical =
                partition_block_canonical(b"A", b"ACA").expect("grid is far below the bound");
            assert_ne!(
                live, canonical,
                "if these ever agree, pick another block — this test is vacuous otherwise"
            );
        }

        /// Unset means `live`, and so does anything the knob does not recognise.
        ///
        /// This is the property the default path depends on. A knob that fell
        /// through to a different rule on an empty string, or that treated a
        /// typo (`FERRO_PARTITION=cannonical`) as an instruction rather than as
        /// noise, would silently re-partition every block for whoever set it —
        /// which for this repo is a representation change, not a configuration.
        #[test]
        fn unset_and_unrecognised_values_select_the_live_rule() {
            for value in [None, Some(""), Some("live")] {
                assert_eq!(
                    partition_rule_from_env(value),
                    PartitionRule::Live,
                    "{value:?}"
                );
            }
            for typo in ["cannonical", "LIVE", "shadow ", "1", "true"] {
                assert_eq!(
                    partition_rule_from_env(Some(typo)),
                    PartitionRule::Live,
                    "{typo} must not select a rule"
                );
            }
        }

        /// The two non-default rules are reachable, and only by their own names.
        #[test]
        fn the_alternative_rules_are_selected_by_name() {
            assert_eq!(
                partition_rule_from_env(Some("shadow")),
                PartitionRule::Shadow
            );
            assert_eq!(
                partition_rule_from_env(Some("canonical")),
                PartitionRule::Canonical
            );
            // The rule this change adds. Without it the knob's newest value was
            // the one value nothing exercised, so a typo in the match arm would
            // have surfaced only as a silent fallback to `Live` during a
            // bake-off — exactly the failure the knob's own doc says it is
            // designed to make loud.
            assert_eq!(
                partition_rule_from_env(Some("canonical-coalesced")),
                PartitionRule::CanonicalCoalesced
            );
        }

        /// This suite's expectations are the live rule's, and the cached read
        /// agrees with the pure mapping above.
        ///
        /// Every other test in this file asserts a *shipped* representation, and
        /// switching the rule moves dozens of them at once (measured: 36 under
        /// `shadow`, 33 under `canonical`). This is the one that names the cause.
        /// It therefore **fails by design** during a deliberate
        /// `FERRO_PARTITION=…` bake-off run — that failure is the message
        /// "you are not measuring the shipped rule", not a defect.
        #[test]
        fn the_suite_runs_on_the_live_rule() {
            assert!(
                std::env::var_os("FERRO_PARTITION").is_none(),
                "this suite's expectations are the live rule's; FERRO_PARTITION must be unset"
            );
            assert_eq!(partition_rule(), PartitionRule::Live);
        }

        /// The canonical arm produces pieces that denote the block, so the arm is
        /// wired to something that works rather than merely something that
        /// compiles.
        ///
        /// Uses the spec's own worked example (`delins.md:44`,
        /// `LRG_199t1:c.850_901`) because it is the block with the most
        /// alternative minimal alignments in the calibration corpus — the regime
        /// where a mis-selected path would still round-trip but claim more
        /// columns than the block's distance.
        #[test]
        fn the_canonical_arm_denotes_the_block_it_partitions() {
            let reference = b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT";
            let result = b"TTCCTCGATGCCTG";
            let pieces =
                partition_block_canonical(reference, result).expect("grid is far below the bound");
            let mut rebuilt = Vec::new();
            let mut cursor = 0usize;
            for piece in &pieces {
                rebuilt.extend_from_slice(&reference[cursor..piece.ref_start]);
                rebuilt.extend_from_slice(&piece.alt);
                cursor = piece.ref_end;
            }
            rebuilt.extend_from_slice(&reference[cursor..]);
            assert_eq!(rebuilt, result.to_vec(), "canonical pieces: {pieces:?}");
            assert_eq!(
                changed_columns_of_pieces(&pieces),
                crate::normalize::seqfirst::align::AlignmentDag::build(reference, result)
                    .edit_distance() as usize,
                "the canonical arm must claim exactly the block's edit distance"
            );
        }

        /// Substitute each piece's payload back into `reference`.
        fn denotation(reference: &[u8], pieces: &[Piece]) -> Vec<u8> {
            let mut rebuilt = Vec::new();
            let mut cursor = 0usize;
            for piece in pieces {
                rebuilt.extend_from_slice(&reference[cursor..piece.ref_start]);
                rebuilt.extend_from_slice(&piece.alt);
                cursor = piece.ref_end;
            }
            rebuilt.extend_from_slice(&reference[cursor..]);
            rebuilt
        }

        /// Apply the terminal coalesce to a piece list, as
        /// `canonicalize_from_sequence` does after its downstream passes.
        fn coalesced(reference: &[u8], pieces: &[Piece]) -> Vec<Piece> {
            let mut pieces = pieces.to_vec();
            coalesce_payload_alignment_split(&mut pieces, reference);
            pieces
        }

        /// The real block behind `NC_000001.10:g.240370952_240370985delinsT`
        /// (GRCh37), the largest single class of representation churn the
        /// canonical arm introduces: 121 of the moved rows are this shape.
        ///
        /// The block holds four `T`s, so the single payload base can align to
        /// any of them — the split is coincidence, and `delins.md:44-47`
        /// recommends the `delins`.
        #[test]
        fn the_coalesce_pass_merges_a_payload_alignment_split() {
            let reference = b"CTCTACCCGGAGCGGCAATACCCCCTCCGCCCCC";
            let result = b"T";

            let canonical =
                partition_block_canonical(reference, result).expect("grid is below the bound");
            assert_eq!(
                canonical.len(),
                2,
                "precondition: canonical splits this block ({canonical:?})"
            );

            let merged = coalesced(reference, &canonical);
            assert_eq!(
                merged.len(),
                1,
                "must restore the single delins ({merged:?})"
            );
            assert_eq!(merged[0].ref_start, 0);
            assert_eq!(merged[0].ref_end, reference.len());
            assert_eq!(
                denotation(reference, &merged),
                denotation(reference, &canonical),
                "re-spelling must not change what the pieces denote"
            );
        }

        /// Two genuinely separate substitutions six unchanged bases apart.
        /// `delins.md:17` requires these to stay individual, and the
        /// net-deletion condition is what enforces it: the payload embeds
        /// trivially, but the span loses no length.
        #[test]
        fn the_net_deletion_condition_keeps_two_substitutions_split() {
            let reference = b"ACCCCCCA";
            let result = b"TCCCCCCT";
            let canonical =
                partition_block_canonical(reference, result).expect("grid is below the bound");
            assert_eq!(canonical.len(), 2, "precondition: two members");
            assert_eq!(
                coalesced(reference, &canonical),
                canonical,
                "a span that loses no length must not be merged"
            );
        }

        /// A **net insertion** must not be merged either. Measured per spelling,
        /// #1421's family runs opposite to #1419/#1420's: for a net insertion the
        /// split form is the canonical one, so a blanket preference for fewer
        /// members gets that family backwards. The net-deletion condition
        /// excludes it structurally — a subsequence embedding cannot lengthen.
        #[test]
        fn a_net_insertion_is_never_merged() {
            let reference = b"ACGT";
            let pieces = vec![
                Piece {
                    ref_start: 0,
                    ref_end: 1,
                    alt: b"AAAA".to_vec(),
                },
                Piece {
                    ref_start: 3,
                    ref_end: 4,
                    alt: b"TTTT".to_vec(),
                },
            ];
            assert_eq!(
                coalesced(reference, &pieces),
                pieces,
                "a net insertion must survive as separate members"
            );
        }

        /// The pass declines a single member rather than re-spelling it, so a
        /// block that was never split cannot be widened.
        #[test]
        fn a_single_member_is_left_alone() {
            let reference = b"ACGTACGT";
            let pieces = vec![Piece {
                ref_start: 2,
                ref_end: 6,
                alt: b"T".to_vec(),
            }];
            assert_eq!(coalesced(reference, &pieces), pieces);
        }

        /// The spec's own worked example, `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`
        /// (real transcript bases). `delins.md:44-47` pins the single `delins`
        /// and explicitly names the four-member split as the description to
        /// avoid; the canonical arm emits *five* members.
        ///
        /// The pure-subsequence predicate cannot see this case — its minimal
        /// split carries one substitution — which is what
        /// `COALESCE_MISMATCH_BUDGET` exists for.
        #[test]
        fn the_coalesce_pass_reaches_the_spec_worked_example() {
            let reference = b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT";
            let result = b"TTCCTCGATGCCTG";

            let canonical =
                partition_block_canonical(reference, result).expect("grid is below the bound");
            assert_eq!(
                canonical.len(),
                5,
                "precondition: canonical fragments the spec's example ({canonical:?})"
            );
            assert!(
                !payload_embeds_as_subsequence(reference, result),
                "precondition: one substituted position puts it out of the exact test's reach"
            );

            let merged = coalesced(reference, &canonical);
            assert_eq!(
                merged.len(),
                1,
                "the spec's recommended single delins must be restored ({merged:?})"
            );
            assert_eq!(merged[0].alt, result.to_vec());
            assert_eq!(
                denotation(reference, &merged),
                denotation(reference, &canonical),
                "re-spelling must not change what the pieces denote"
            );
        }

        /// The budget is one substituted position, and the value is what keeps
        /// two genuinely separate substitutions apart even before the
        /// net-deletion condition is consulted.
        #[test]
        fn the_mismatch_budget_is_one_substituted_position() {
            assert_eq!(COALESCE_MISMATCH_BUDGET, 1);
            // the spec's example: four coincidental deletions plus one substitution
            assert!(payload_embeds_within_budget(
                b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT",
                b"TTCCTCGATGCCTG",
                1
            ));
            // two substitutions six bases apart need two, and are refused at one
            assert!(!payload_embeds_within_budget(b"ACCCCCCA", b"TCCCCCCT", 1));
            assert!(payload_embeds_within_budget(b"ACCCCCCA", b"TCCCCCCT", 2));
            // a budget cannot rescue a payload longer than its span
            assert!(!payload_embeds_within_budget(b"AC", b"ACGT", 9));
            // budget 0 is exactly the subsequence test
            assert!(payload_embeds_within_budget(b"ACGT", b"AT", 0));
            assert!(!payload_embeds_within_budget(b"ACGT", b"TA", 0));
        }

        /// A deterministic, seeded byte drawn from a small alphabet.
        ///
        /// Not `rand`: the table must enumerate the *same* cases on every run, or
        /// a failure cannot be reproduced from the printed case description.
        fn draw(seed: &mut u64, alphabet: &[u8]) -> u8 {
            // xorshift64*, inline so the table has no dependency and no RNG state
            // shared with anything else.
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;
            alphabet[(*seed % alphabet.len() as u64) as usize]
        }

        /// What the payload does to the span's length, which is the axis
        /// `delins.md:17` and #1421 turn on.
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        enum Shape {
            /// Pure deletion members: the payload embeds by construction.
            Deletion,
            /// Deletion members plus one substituted position.
            OneMismatch,
            /// Members that replace as many bases as they consume.
            LengthNeutral,
            /// Members that add bases.
            NetInsertion,
        }

        /// One generated case, carrying enough context to reproduce a failure.
        struct Case {
            span: Vec<u8>,
            pieces: Vec<Piece>,
            gap: usize,
            shape: Shape,
        }

        impl Case {
            fn describe(&self) -> String {
                format!(
                    "gap={} shape={:?} members={} span={}",
                    self.gap,
                    self.shape,
                    self.pieces.len(),
                    String::from_utf8_lossy(&self.span)
                )
            }
        }

        /// Build one case: `members` members of `width` reference bases each,
        /// separated by `gap` unchanged bases, over a span drawn from `alphabet`.
        fn build(
            members: usize,
            width: usize,
            gap: usize,
            shape: Shape,
            alphabet: &[u8],
            seed: &mut u64,
        ) -> Case {
            let mut span = Vec::new();
            let mut pieces = Vec::new();
            for i in 0..members {
                if i > 0 {
                    for _ in 0..gap {
                        span.push(draw(seed, alphabet));
                    }
                }
                let start = span.len();
                for _ in 0..width {
                    span.push(draw(seed, alphabet));
                }
                let alt = match shape {
                    Shape::Deletion => Vec::new(),
                    // One substituted position, in the first member only.
                    Shape::OneMismatch if i == 0 => vec![b'N'],
                    Shape::OneMismatch => Vec::new(),
                    Shape::LengthNeutral => vec![b'N'; width],
                    Shape::NetInsertion => vec![b'N'; width + 2],
                };
                pieces.push(Piece {
                    ref_start: start,
                    ref_end: span.len(),
                    alt,
                });
            }
            Case {
                span,
                pieces,
                gap,
                shape,
            }
        }

        /// The generated corpus. Geometry axes are the ones the corpora showed
        /// matter: separations spanning the coincidental range (0-8) and the real
        /// multi-member cis range (9-1562), several member counts and widths, and
        /// every length shape.
        fn cases() -> Vec<Case> {
            const GAPS: &[usize] = &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 21, 26, 32, 58, 100, 166, 234, 400, 604,
                973, 1562,
            ];
            const WIDTHS: &[usize] = &[1, 2, 3, 5, 8];
            const MEMBERS: &[usize] = &[2, 3, 4, 5];
            const SHAPES: &[Shape] = &[
                Shape::Deletion,
                Shape::OneMismatch,
                Shape::LengthNeutral,
                Shape::NetInsertion,
            ];
            // A homopolymer maximises accidental embedding; four bases is the
            // ordinary case. Both matter — coincidence is the whole subject.
            const ALPHABETS: &[&[u8]] = &[b"A", b"AC", b"ACGT"];

            let mut seed = 0x9E37_79B9_7F4A_7C15u64;
            let mut out = Vec::new();
            for &gap in GAPS {
                for &width in WIDTHS {
                    for &members in MEMBERS {
                        for &shape in SHAPES {
                            for alphabet in ALPHABETS {
                                out.push(build(members, width, gap, shape, alphabet, &mut seed));
                            }
                        }
                    }
                }
            }
            out
        }

        /// Drive the generated corpus against the pass's *invariants*.
        ///
        /// Deliberately not a lookup table of expected merge/no-merge: computing
        /// the expectation from the same conditions the pass applies would only
        /// restate the implementation. These are properties that hold for
        /// reasons outside it — sequence preservation, the spec's rules on when
        /// members stay individual, and internal consistency.
        #[test]
        fn the_coalesce_pass_holds_its_invariants_over_the_generated_corpus() {
            let cases = cases();
            assert!(cases.len() >= 5000, "corpus is {} cases", cases.len());

            let mut merged_count = 0usize;
            for case in &cases {
                let mut pieces = case.pieces.clone();
                coalesce_payload_alignment_split(&mut pieces, &case.span);
                let merged = pieces.len() == 1 && case.pieces.len() > 1;
                if merged {
                    merged_count += 1;
                }
                let ctx = case.describe();

                // 1. It is a re-spelling: the denoted sequence never moves.
                assert_eq!(
                    denotation(&case.span, &pieces),
                    denotation(&case.span, &case.pieces),
                    "denotation changed [{ctx}]"
                );

                // 2. It never merges across a wide separation -- those are
                //    separate variants (`delins.md:17`, `general.md:34`).
                if case.gap > COALESCE_MAX_SEPARATION {
                    assert!(!merged, "merged across a {}-base gap [{ctx}]", case.gap);
                }

                // 3. It never merges a shape that does not lose length. A
                //    length-neutral block is substitutions; a longer one is
                //    #1421's net-insertion family, where the split form is the
                //    canonical one.
                if matches!(case.shape, Shape::LengthNeutral | Shape::NetInsertion) {
                    assert!(!merged, "merged a non-shortening block [{ctx}]");
                }

                // 4. It only ever reduces the member count.
                assert!(
                    pieces.len() <= case.pieces.len(),
                    "member count grew [{ctx}]"
                );

                // 5. Idempotent: a second application changes nothing, so the
                //    output is a fixed point of the pass.
                let mut again = pieces.clone();
                coalesce_payload_alignment_split(&mut again, &case.span);
                assert_eq!(again, pieces, "not idempotent [{ctx}]");
            }

            // The corpus must actually exercise the merging path, or every
            // assertion above is vacuous.
            assert!(
                merged_count > 100,
                "only {merged_count} of {} cases merged -- the corpus is not exercising the pass",
                cases.len()
            );
        }

        /// Widening the separation can only ever *withdraw* a merge, never create
        /// one. Checked by holding everything else fixed and walking the gap up.
        #[test]
        fn widening_the_separation_is_monotone() {
            let mut seed = 0x0DDB_1A5E_5BAD_5EEDu64;
            for &width in &[1usize, 3, 8] {
                for &members in &[2usize, 3] {
                    let mut still_merging = true;
                    for gap in 0..=(COALESCE_MAX_SEPARATION + 4) {
                        let case = build(members, width, gap, Shape::Deletion, b"ACGT", &mut seed);
                        let mut pieces = case.pieces.clone();
                        coalesce_payload_alignment_split(&mut pieces, &case.span);
                        let merged = pieces.len() == 1;
                        // Asserted **before** `still_merging` is updated, and on
                        // `!merged || still_merging`. Both halves matter: the
                        // earlier form cleared the flag first and tested
                        // `merged || !still_merging`, which is a tautology — a
                        // merge satisfies the left side, and a non-merge has just
                        // set `still_merging = false` and so satisfies the right.
                        // It could not fail for any input. The condition was also
                        // inverted: a *non*-merge while still merging is the
                        // normal transition this walk is looking for, and a merge
                        // *after* one is the violation.
                        assert!(
                            !merged || still_merging,
                            "merge reappeared after it stopped [{}]",
                            case.describe()
                        );
                        if !merged {
                            still_merging = false;
                        }
                    }
                }
            }
        }

        /// **Known limit, pinned rather than fixed.**
        ///
        /// Three of the 592 real multi-member cis alleles are merged by this pass
        /// when they should stay individual. They are genuine separate variants
        /// that happen to sit as close together as a coincidental split — their
        /// separations are 1 and 4 unchanged bases, inside the range where every
        /// correct rejoin also lives (ClinVar's 78 correct rejoins run 1..8).
        ///
        /// **Separation cannot distinguish them, and neither can any other
        /// property of the sequence.** `delins.md:79-84` gives the spec's own
        /// discriminator: *"the two variants may have been reported (or might
        /// occur) individually"* — provenance, which lives in the input's
        /// spelling and nowhere in the reference or the observed bases.
        ///
        /// Reading the input's spelling to recover it is exactly what forfeits
        /// confluence, and the cost is measured rather than assumed: the weight
        /// bound above compares against the input's own edits, and that single
        /// input-relative comparand cost **427 converged classes** per direction
        /// until this pass was moved past it. Confluence and provenance are
        /// incompatible for this class; #1235 chooses confluence.
        ///
        /// If provenance is ever needed, carry it beside the string as metadata —
        /// not in it. This assertion exists so that anyone reaching for the input
        /// form has to move it deliberately and read this first.
        #[test]
        fn close_separate_variants_are_a_known_limit_not_a_bug() {
            let mut seed = 0x5EED_1235_5EED_1235u64;
            for &gap in &[1usize, 4] {
                let case = build(2, 3, gap, Shape::Deletion, b"ACGT", &mut seed);
                let mut pieces = case.pieces.clone();
                coalesce_payload_alignment_split(&mut pieces, &case.span);
                assert_eq!(
                    pieces.len(),
                    1,
                    "pinned: the pass merges at separation {gap}, which is right for a \
                     coincidental split and wrong for two genuinely separate variants — \
                     the sequence cannot tell them apart [{}]",
                    case.describe()
                );
            }
        }

        #[test]
        fn payload_embedding_is_an_ordered_subsequence_test() {
            assert!(payload_embeds_as_subsequence(b"ACGT", b"AT"));
            assert!(payload_embeds_as_subsequence(b"ACGT", b"ACGT"));
            // order matters: the bases are present but not in this order
            assert!(!payload_embeds_as_subsequence(b"ACGT", b"TA"));
            // a base the block does not hold at all
            assert!(!payload_embeds_as_subsequence(b"AAAA", b"C"));
            // each payload base must consume its own reference base
            assert!(!payload_embeds_as_subsequence(b"A", b"AA"));
            // vacuous for an empty payload; the caller's member gate declines it
            assert!(payload_embeds_as_subsequence(b"ACGT", b""));
        }
    }

    mod seqfirst_shadow_audit {
        use super::*;

        /// Substitute each piece's payload into `reference`.
        ///
        /// Returns `None` if the pieces are not disjoint and ascending, so a
        /// caller cannot mistake an ill-formed piece list for a denotation.
        fn denotes(reference: &[u8], pieces: &[Piece]) -> Option<Vec<u8>> {
            let mut cursor = 0;
            let mut out = Vec::new();
            for piece in pieces {
                if piece.ref_start < cursor
                    || piece.ref_end < piece.ref_start
                    || piece.ref_end > reference.len()
                {
                    return None;
                }
                out.extend_from_slice(&reference[cursor..piece.ref_start]);
                out.extend_from_slice(&piece.alt);
                cursor = piece.ref_end;
            }
            out.extend_from_slice(&reference[cursor..]);
            Some(out)
        }

        /// `(ref_start, ref_end, alt)` triples, for terse expectations.
        fn shape(pieces: &[Piece]) -> Vec<(usize, usize, &[u8])> {
            pieces
                .iter()
                .map(|p| (p.ref_start, p.ref_end, p.alt.as_slice()))
                .collect()
        }

        /// Run both splitters on one block and check the two invariants that hold
        /// for every disagreement found: both sides denote the alternate block,
        /// and the sequence-first side never declines.
        ///
        /// Returns `(live, shadow)` for the caller to pin exactly.
        fn both(reference: &[u8], result: &[u8]) -> (Vec<Piece>, Vec<Piece>) {
            let live = partition_block(reference, result);
            let shadow = partition_block_sequence_first(
                reference,
                result,
                crate::normalize::seqfirst::MIN_SEPARATION,
            )
            .expect("no harvested disagreement had the sequence-first side decline");
            assert_eq!(
                denotes(reference, &live).as_deref(),
                Some(result),
                "live pieces must denote the alternate block: {live:?}"
            );
            assert_eq!(
                denotes(reference, &shadow).as_deref(),
                Some(result),
                "shadow pieces must denote the alternate block: {shadow:?}"
            );
            assert_ne!(shape(&live), shape(&shadow), "this block must disagree");
            (live, shadow)
        }

        /// `A1` (7 035 blocks). `partition_block` has no single-gap explanation
        /// for two separate deletions, so it emits the whole block as one
        /// `delins`; the sequence-first split finds both deletions.
        #[test]
        fn shadow_splits_a_block_live_kept_whole() {
            let (live, shadow) = both(b"ACAC", b"CA");
            assert_eq!(shape(&live), [(0, 4, b"CA".as_slice())]);
            assert_eq!(
                shape(&shadow),
                [(0, 1, b"".as_slice()), (3, 4, b"".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&live), 4);
            assert_eq!(changed_columns_of_pieces(&shadow), 2);
        }

        /// `A2` (6 915 blocks). The block is a 5' insertion plus a 3' deletion.
        /// `partition_block`'s lone gap parks at one end, so the offset columns
        /// pair up position-wise and read as a run of substitutions — three
        /// pieces claiming three changed columns for what is one base in and one
        /// base out.
        #[test]
        fn shadow_replaces_a_coincidental_substitution_run_with_one_indel_pair() {
            let (live, shadow) = both(b"AACCA", b"TAACC");
            assert_eq!(
                shape(&live),
                [
                    (0, 1, b"T".as_slice()),
                    (2, 3, b"A".as_slice()),
                    (4, 5, b"C".as_slice()),
                ]
            );
            assert_eq!(
                shape(&shadow),
                [(0, 0, b"T".as_slice()), (4, 5, b"".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&live), 3);
            assert_eq!(changed_columns_of_pieces(&shadow), 2);
        }

        /// `A3` (4 567 blocks). Same piece count; the live first piece spans one
        /// base more than it needs to, and pays for it with a substitution the
        /// sequence-first split renders as a deletion.
        #[test]
        fn shadow_narrows_a_piece_live_over_widened() {
            let (live, shadow) = both(b"ACCA", b"CC");
            assert_eq!(
                shape(&live),
                [(0, 2, b"".as_slice()), (3, 4, b"C".as_slice())]
            );
            assert_eq!(
                shape(&shadow),
                [(0, 1, b"".as_slice()), (3, 4, b"".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&live), 3);
            assert_eq!(changed_columns_of_pieces(&shadow), 2);
        }

        /// `A4` (379 blocks). The rarest of the four: the sequence-first first
        /// piece is *wider* than live's, and the total is still lower.
        #[test]
        fn shadow_widens_one_piece_and_still_lowers_the_total() {
            let (live, shadow) = both(b"AACAC", b"CA");
            assert_eq!(
                shape(&live),
                [(0, 1, b"C".as_slice()), (2, 5, b"".as_slice())]
            );
            assert_eq!(
                shape(&shadow),
                [(0, 2, b"".as_slice()), (4, 5, b"".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&live), 4);
            assert_eq!(changed_columns_of_pieces(&shadow), 3);
        }

        /// `B1` (8 044 blocks) — **the sequence-first answer is less minimal.**
        ///
        /// `CA` is inserted into a run where it can sit at more than one offset,
        /// so no alignment node is common to every minimal alignment there and
        /// the member absorbs `reference[0..1]` rather than committing to the
        /// 5'-most tie-break live picks. Both denote `CAACAG`; live spends 3
        /// changed columns, the sequence-first split 4.
        #[test]
        fn shadow_widens_a_member_over_an_indel_it_will_not_localise() {
            let (live, shadow) = both(b"ACAA", b"CAACAG");
            assert_eq!(
                shape(&live),
                [(0, 0, b"CA".as_slice()), (3, 4, b"G".as_slice())]
            );
            assert_eq!(
                shape(&shadow),
                [(0, 1, b"CAA".as_slice()), (3, 4, b"G".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&live), 3);
            assert_eq!(changed_columns_of_pieces(&shadow), 4);
        }

        /// `B2` (3 088 blocks, 1 265 of them equal-length) — **the sequence-first
        /// answer is less minimal, and this is the class that reaches the curated
        /// corpora.**
        ///
        /// `seqfirst::MIN_SEPARATION` is 2, so one unchanged reference base
        /// between two runs of change merges them. `partition_block` splits
        /// there: it compares equal-length blocks position-wise, and its own
        /// `MIN_PIECE_SEPARATION` of 2 gates only net insertions.
        ///
        /// The three blocks below are, in order, the smallest instance, the block
        /// `issue_1235_transcript_axes::noncoding_axis_converges_on_separated_changes`
        /// produces, and the block
        /// `issue_1235_cis_allele_confluence::issue_1232_spanning_delins_splits_at_unchanged_base`
        /// produces. Both of those tests are named for the split the
        /// sequence-first splitter declines to make.
        #[test]
        fn shadow_merges_across_one_unchanged_base_where_live_splits() {
            let (live, shadow) = both(b"ACA", b"CC");
            assert_eq!(
                shape(&live),
                [(0, 1, b"".as_slice()), (2, 3, b"C".as_slice())]
            );
            assert_eq!(shape(&shadow), [(0, 3, b"CC".as_slice())]);
            assert_eq!(changed_columns_of_pieces(&live), 2);
            assert_eq!(changed_columns_of_pieces(&shadow), 3);

            // #1235's `noncoding_axis_converges_on_separated_changes`: an
            // equal-length block, two substitutions one unchanged base apart.
            let (live, shadow) = both(b"CAA", b"TAG");
            assert_eq!(
                shape(&live),
                [(0, 1, b"T".as_slice()), (2, 3, b"G".as_slice())]
            );
            assert_eq!(shape(&shadow), [(0, 3, b"TAG".as_slice())]);
            assert_eq!(changed_columns_of_pieces(&live), 2);
            assert_eq!(changed_columns_of_pieces(&shadow), 3);

            // #1232's `spanning_delins_splits_at_unchanged_base`: the split the
            // issue is named for is the one that goes away.
            let (live, shadow) = both(b"CAATT", b"TA");
            assert_eq!(
                shape(&live),
                [(0, 3, b"".as_slice()), (4, 5, b"A".as_slice())]
            );
            assert_eq!(shape(&shadow), [(0, 5, b"TA".as_slice())]);
            assert_eq!(changed_columns_of_pieces(&live), 4);
            assert_eq!(changed_columns_of_pieces(&shadow), 5);
        }

        /// `B3` (1 block) — **the sequence-first answer is less minimal.** The
        /// only harvested block where the two splits have the same piece count,
        /// the sequence-first first piece is narrower, and the total is still
        /// higher because the second piece grows by more than the first shrank.
        #[test]
        fn shadow_narrows_one_piece_and_pays_for_it_in_the_next() {
            let (live, shadow) = both(b"GCCATA", b"CATAAC");
            assert_eq!(
                shape(&live),
                [(0, 3, b"CAT".as_slice()), (4, 6, b"AC".as_slice())]
            );
            assert_eq!(
                shape(&shadow),
                [(0, 3, b"C".as_slice()), (5, 6, b"AAC".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&live), 5);
            assert_eq!(changed_columns_of_pieces(&shadow), 6);
        }

        /// `C1` (232 blocks). Both denote the alternate block and both spend the
        /// same number of changed columns; they only divide the change
        /// differently. Which is preferable is not decided here — `delins.md:17`
        /// and ferro's own `split_buys_no_higher_priority_type` policy are what
        /// reach these, and the two blocks below fall to them in opposite
        /// directions.
        #[test]
        fn equal_weight_different_division() {
            let (live, shadow) = both(b"ACAGCG", b"CGCTGT");
            assert_eq!(shape(&live), [(0, 6, b"CGCTGT".as_slice())]);
            assert_eq!(
                shape(&shadow),
                [(0, 3, b"C".as_slice()), (5, 6, b"TGT".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&live), 6);
            assert_eq!(changed_columns_of_pieces(&shadow), 6);

            let (live, shadow) = both(b"AGAGT", b"GAAGAG");
            assert_eq!(
                shape(&live),
                [(0, 2, b"GA".as_slice()), (4, 5, b"AG".as_slice())]
            );
            assert_eq!(
                shape(&shadow),
                [(0, 1, b"GAA".as_slice()), (4, 5, b"".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&live), 4);
            assert_eq!(changed_columns_of_pieces(&shadow), 4);
        }

        /// The one property that held over all 30 261 harvested disagreements:
        /// both piece lists are disjoint, ascending, and denote the alternate
        /// block. `denotes` returns `None` on a piece list that is not, so this
        /// is what makes "neither splitter is ever wrong" a checked claim rather
        /// than a reading of the log.
        #[test]
        fn every_pinned_class_denotes_the_alternate_block_on_both_sides() {
            for (reference, result) in [
                (b"ACAC".as_slice(), b"CA".as_slice()),
                (b"AACCA", b"TAACC"),
                (b"ACCA", b"CC"),
                (b"AACAC", b"CA"),
                (b"ACAA", b"CAACAG"),
                (b"ACA", b"CC"),
                (b"CAA", b"TAG"),
                (b"CAATT", b"TA"),
                (b"GCCATA", b"CATAAC"),
                (b"ACAGCG", b"CGCTGT"),
                (b"AGAGT", b"GAAGAG"),
            ] {
                // `both` asserts the denotation on each side and that the block
                // does in fact disagree.
                both(reference, result);
            }
        }
    }

    /// Minimal member rendering — [`shrink_pieces_to_differences`].
    ///
    /// The decision this pins: keep the sequence-derived partition, and spell
    /// each member as the hull of its own differences. A whole ambiguous run is
    /// described as one unit only where that unit has an edit type of its own
    /// (`dup`, `NNN[k]`), never as a `delins` spanning the run.
    mod minimal_member_rendering {
        use super::*;

        /// `(ref_start, ref_end, alt)` triples, for terse expectations.
        fn shape(pieces: &[Piece]) -> Vec<(usize, usize, &[u8])> {
            pieces
                .iter()
                .map(|p| (p.ref_start, p.ref_end, p.alt.as_slice()))
                .collect()
        }

        /// The sequence-first split of `reference -> result`, shrunk.
        ///
        /// At [`MIN_SEPARATION_NO_FRAME`], the threshold
        /// `canonicalize_from_sequence` chooses on an axis with no reading
        /// frame — which is the axis every block below was harvested from.
        fn shrunk_sequence_first(reference: &[u8], result: &[u8]) -> Vec<Piece> {
            let mut pieces =
                partition_block_sequence_first(reference, result, MIN_SEPARATION_NO_FRAME)
                    .expect("the sequence-first split must not decline on these blocks");
            shrink_pieces_to_differences(&mut pieces, reference);
            pieces
        }

        /// Every trimmed `(block, alternate)` pair drawn from `{A,C,G}` with
        /// both sides of length 1..=4, excluding the pairs that trim away to no
        /// change at all (which `canonicalize_from_sequence` returns early on).
        ///
        /// Small on purpose: the blocks the shrink acts on are already trimmed
        /// and rarely long, and an exhaustive sweep of short blocks covers the
        /// shapes — homopolymer, tandem 2-mer, net insertion, net deletion,
        /// equal length — more evenly than a longer random sample would.
        fn enumerated_blocks() -> Vec<(Vec<u8>, Vec<u8>)> {
            fn words(alphabet: &[u8], max_len: usize) -> Vec<Vec<u8>> {
                let mut out: Vec<Vec<u8>> = Vec::new();
                let mut level = vec![Vec::new()];
                for _ in 0..max_len {
                    let mut next = Vec::new();
                    for word in &level {
                        for base in alphabet {
                            let mut grown = word.clone();
                            grown.push(*base);
                            next.push(grown);
                        }
                    }
                    out.extend(next.iter().cloned());
                    level = next;
                }
                out
            }

            let words = words(b"ACG", 4);
            let mut out = Vec::new();
            for reference in &words {
                for result in &words {
                    let (lo, hi_ref, hi_alt) = trim_common_flanks(reference, result);
                    if lo == hi_ref && lo == hi_alt {
                        continue;
                    }
                    out.push((reference[lo..hi_ref].to_vec(), result[lo..hi_alt].to_vec()));
                }
            }
            out
        }

        /// Substitute each piece's payload into `reference` (see the audit
        /// module's `denotes`, duplicated here so this module stands alone).
        fn denotes(reference: &[u8], pieces: &[Piece]) -> Vec<u8> {
            let mut cursor = 0;
            let mut out = Vec::new();
            for piece in pieces {
                assert!(piece.ref_start >= cursor, "pieces overlap: {pieces:?}");
                out.extend_from_slice(&reference[cursor..piece.ref_start]);
                out.extend_from_slice(&piece.alt);
                cursor = piece.ref_end;
            }
            out.extend_from_slice(&reference[cursor..]);
            out
        }

        /// Class `B1`: the sequence-first split widens a member over an indel it
        /// will not localise. Shrinking each member to its own differences
        /// recovers the minimal spelling — on all three blocks, the very piece
        /// list `partition_block` produces.
        ///
        /// The blocks are the audit's pinned `B1` example plus two more taken
        /// from the same harvest, chosen to cover both directions and the widest
        /// widening seen: a net insertion widened by one base, a net deletion
        /// widened by one base, and an 8-base homopolymer where the derived
        /// member is 8 columns wider than the change it holds.
        #[test]
        fn shrinking_collapses_the_b1_blocks_onto_the_minimal_split() {
            for (reference, result) in [
                (b"ACAA".as_slice(), b"CAACAG".as_slice()),
                (b"ACCAA", b"CAG"),
                (b"AAAAAAAATA", b"TAAAAAAAAATG"),
            ] {
                let live = partition_block(reference, result);
                let shrunk = shrunk_sequence_first(reference, result);
                assert_eq!(
                    shape(&shrunk),
                    shape(&live),
                    "block {} -> {}",
                    String::from_utf8_lossy(reference),
                    String::from_utf8_lossy(result)
                );
                assert_eq!(denotes(reference, &shrunk), result);
                assert_eq!(
                    changed_columns_of_pieces(&shrunk),
                    changed_columns_of_pieces(&live)
                );
            }
        }

        /// The audit's `B1` block in full, with the weights it recorded: the
        /// sequence-first split spends 4 changed columns, shrinking brings it
        /// back to live's 3, and the widened `delins` of `CAA` over one base
        /// becomes the insertion of `CA` those bases denote.
        #[test]
        fn the_pinned_b1_block_loses_its_spanning_delins() {
            let reference = b"ACAA";
            let raw = partition_block_sequence_first(reference, b"CAACAG", MIN_SEPARATION_NO_FRAME)
                .expect("must partition");
            assert_eq!(
                shape(&raw),
                [(0, 1, b"CAA".as_slice()), (3, 4, b"G".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&raw), 4);

            let shrunk = shrunk_sequence_first(reference, b"CAACAG");
            assert_eq!(
                shape(&shrunk),
                [(0, 0, b"CA".as_slice()), (3, 4, b"G".as_slice())]
            );
            assert_eq!(changed_columns_of_pieces(&shrunk), 3);
        }

        /// The guard against over-shrinking, at the piece level: a deletion, an
        /// insertion, a duplicating insertion and a substitution are returned
        /// exactly as they arrived. None of them can share a flank between span
        /// and payload, so the shrink reaches only `delins`-shaped pieces.
        #[test]
        fn a_deletion_insertion_dup_and_substitution_are_left_alone() {
            let ref_bytes = b"AAAAAG";
            for piece in [
                // A 2-base deletion inside the A-run: no payload.
                Piece {
                    ref_start: 1,
                    ref_end: 3,
                    alt: Vec::new(),
                },
                // An insertion of two A into the same run: no reference span.
                Piece {
                    ref_start: 4,
                    ref_end: 4,
                    alt: b"AA".to_vec(),
                },
                // A duplicating insertion — `is_tandem_duplication` holds.
                Piece {
                    ref_start: 3,
                    ref_end: 3,
                    alt: b"A".to_vec(),
                },
                // A substitution: one base, and it differs.
                Piece {
                    ref_start: 5,
                    ref_end: 6,
                    alt: b"T".to_vec(),
                },
            ] {
                let mut pieces = [piece.clone()];
                shrink_pieces_to_differences(&mut pieces, ref_bytes);
                assert_eq!(pieces[0], piece);
            }
        }

        /// A piece whose span and payload are identical would trim to nothing.
        /// It is left whole instead, so `anchor_for_piece` never sees a member
        /// with no extent.
        #[test]
        fn a_piece_that_trims_away_entirely_is_left_alone() {
            let piece = Piece {
                ref_start: 0,
                ref_end: 1,
                alt: b"A".to_vec(),
            };
            let mut pieces = [piece.clone()];
            shrink_pieces_to_differences(&mut pieces, b"A");
            assert_eq!(pieces[0], piece);
        }

        /// The live splitter's own pieces are already minimal, so adding the
        /// shrink to `canonicalize_from_sequence` cannot move them.
        ///
        /// `best_alignment` maximises matched columns, and a piece that shared a
        /// flank between its span and its payload would mean a better placement
        /// of the single gap was available — so the property is a consequence of
        /// the search, not a coincidence. Enumerated over the same 14 280 blocks
        /// as the ordering test: `partition_block`'s answer is a fixed point of
        /// the shrink on every one.
        #[test]
        fn shrinking_never_moves_a_piece_the_live_splitter_produced() {
            let mut compared = 0usize;
            for (block, alt) in enumerated_blocks() {
                let pieces = partition_block(&block, &alt);
                let mut shrunk = pieces.clone();
                shrink_pieces_to_differences(&mut shrunk, &block);
                assert_eq!(
                    shape(&shrunk),
                    shape(&pieces),
                    "the shrink moved a live piece on {} -> {}",
                    String::from_utf8_lossy(&block),
                    String::from_utf8_lossy(&alt)
                );
                compared += 1;
            }
            assert_eq!(compared, 14_280);
        }

        /// Shrinking commutes with the 3'-shift.
        ///
        /// `canonicalize_from_sequence` shifts, coalesces, then shrinks; the
        /// opposite order — shrink, then shift and coalesce — is not obviously
        /// the same pipeline, because shrinking can turn a `delins` into a pure
        /// indel and only pure indels shift. Enumerated over every block pair
        /// drawn from `{A,C,G}` of length 1..=4 on each side, at both separation
        /// thresholds (14 280 partitioned blocks per threshold), the two orders
        /// agree on every block.
        #[test]
        fn shrinking_before_and_after_the_three_prime_shift_agree() {
            let blocks = enumerated_blocks();
            for min_separation in [
                MIN_SEPARATION_NO_FRAME,
                crate::normalize::seqfirst::MIN_SEPARATION,
            ] {
                let mut compared = 0usize;
                for (block, alt) in &blocks {
                    let Some(pieces) = partition_block_sequence_first(block, alt, min_separation)
                    else {
                        continue;
                    };

                    let mut shift_first = pieces.clone();
                    shift_pieces(&mut shift_first, block, ShuffleDirection::ThreePrime);
                    coalesce_adjacent_pieces(&mut shift_first);
                    shrink_pieces_to_differences(&mut shift_first, block);

                    let mut shrink_first = pieces.clone();
                    shrink_pieces_to_differences(&mut shrink_first, block);
                    shift_pieces(&mut shrink_first, block, ShuffleDirection::ThreePrime);
                    coalesce_adjacent_pieces(&mut shrink_first);

                    assert_eq!(
                        shape(&shift_first),
                        shape(&shrink_first),
                        "orders disagree at min_separation {min_separation} on {} -> {}",
                        String::from_utf8_lossy(block),
                        String::from_utf8_lossy(alt)
                    );
                    compared += 1;
                }
                // Every enumerated block partitions; none declines.
                assert_eq!(compared, blocks.len());
            }
        }
    }

    fn flip_base(base: u8) -> u8 {
        match base {
            b'A' => b'C',
            b'C' => b'A',
            b'G' => b'T',
            _ => b'G',
        }
    }

    /// `complement_base` must refuse a byte that is not a base.
    ///
    /// Delegating to `sequence::reverse_complement` made this function total —
    /// that helper's fallback arm is `other => other`, so an unrecognised byte
    /// came back unchanged and every `?` refusal downstream was unreachable. An
    /// inversion over such a byte would then have been accepted with the byte
    /// "complemented" to itself.
    #[test]
    fn complement_base_refuses_a_non_base_byte() {
        // Recognised, both cases, folded to uppercase.
        assert_eq!(complement_base(b'A'), Some(b'T'));
        assert_eq!(complement_base(b'a'), Some(b'T'));
        assert_eq!(complement_base(b'N'), Some(b'N'));
        // `U` complements as `A`, matching the RNA alphabet.
        assert_eq!(complement_base(b'U'), Some(b'A'));

        // Not bases — the refusal path this test exists for.
        for byte in *b"XZ-*0 " {
            assert_eq!(
                complement_base(byte),
                None,
                "`{}` is not a base and must be refused",
                byte as char
            );
        }
    }

    /// The refusal reaches `apply_edits_to_window`, which must decline rather
    /// than fabricate an inversion over bytes it cannot complement.
    #[test]
    fn inverting_a_span_with_a_non_base_byte_is_refused() {
        let reference = b"ACGTXCGTAC";
        let edits = [GEdit::Inv { s: 3, e: 6 }];
        assert_eq!(
            apply_edits_to_window(&edits, reference, 1),
            None,
            "an inversion spanning a non-base byte must be refused"
        );
        // The same inversion clear of the bad byte still applies.
        let clean = [GEdit::Inv { s: 6, e: 9 }];
        assert!(apply_edits_to_window(&clean, reference, 1).is_some());
    }

    /// Ties between gap placements, broken by what they separate (#1262).
    ///
    /// In a repeat tract every placement of the gap matches the same number of
    /// columns, so the score cannot choose and the tie-break decides the
    /// description. msto's `high` corpus states the rule the spec gives for
    /// that: two variants separated by one or more unchanged nucleotides are
    /// described individually and **not** as a `delins` (`delins.md:17`, cited
    /// by #1232). `general.md:56` does not rank `delins` last — `delins` is
    /// absent from that list, which ranks the type of one description of one
    /// change, the way #1231/#1233 apply it member-wise.
    mod tie_break_by_separation {
        use super::*;

        #[test]
        fn a_tied_placement_that_separates_the_change_wins() {
            // #1262's authored block. Deleting `ref[0]`, `ref[1]` or `ref[2]`
            // each leaves exactly one matched column, so all three placements
            // tie. Only the third leaves the match *between* the two changes,
            // which is what makes them separate members.
            assert_eq!(
                partition_block(b"AAA", b"CA"),
                vec![
                    Piece {
                        ref_start: 0,
                        ref_end: 1,
                        alt: b"C".to_vec()
                    },
                    Piece {
                        ref_start: 2,
                        ref_end: 3,
                        alt: Vec::new()
                    },
                ],
                "a substitution and a deletion, separated by the unchanged base"
            );
        }

        #[test]
        fn a_three_way_tie_still_yields_the_separating_placement() {
            // #1232's block at an untrimmed extent. All three placements match
            // exactly one column, so they DO tie and the tie-break runs; only
            // the placement that leaves the match between the two changes gives
            // two members. An earlier revision of this comment claimed the
            // placements did not tie and that the score decided — measured,
            // that is false, and the test's name is the part that is wrong
            // rather than its assertion.
            assert_eq!(partition_block(b"TTTTT", b"TA").len(), 2);
        }

        #[test]
        fn the_tie_break_leaves_the_contiguous_corpus_alone() {
            // Cluster A stays whole: these have no tie to break, because their
            // best placement is unique.
            assert_eq!(
                partition_block(b"CAGTGACTAG", b"TGTCACGACT").len(),
                1,
                "#1040"
            );
            assert_eq!(partition_block(b"GCT", b"AGC").len(), 1, "#1034");
            assert_eq!(partition_block(b"CTATAG", b"AAACCCC").len(), 1, "#422");
        }
    }

    /// The two-gap insertion form (#1260).
    ///
    /// A single gap cannot express *insertion, retained reference, insertion*:
    /// it has one contiguous indel to place, so it must either swallow the
    /// retained base into a spanning `delins` or leave one of the insertions
    /// unaccounted for. PR #1285 named this directly — for #1260 the answer is
    /// "a two-gap alignment — insertion, retained base, insertion — which
    /// `best_alignment`'s single-gap restriction cannot express".
    mod two_gap_insertion_alignment {
        use super::*;

        #[test]
        fn a_retained_reference_between_two_insertions_is_expressible() {
            // #1260's trimmed block. The single-gap search scores 0 matched
            // columns either way it places the 2 nt gap; keeping the reference
            // `A` intact between two 1 nt insertions matches it, which is the
            // most any alignment of this pair can match.
            assert_eq!(
                best_alignment(b"A", b"CAC"),
                Some(vec![(None, Some(0)), (Some(0), Some(1)), (None, Some(2))]),
                "the whole reference should survive between the two insertions"
            );
        }

        #[test]
        fn the_two_gap_form_partitions_into_two_pure_insertions() {
            // The point of expressing it: two members, each a pure insertion,
            // separated by the retained base. That preference is `delins.md:17`
            // plus ferro policy, not `general.md:56`, whose list does not
            // contain `delins`.
            let columns = best_alignment(b"A", b"CAC").expect("aligned");
            assert_eq!(
                pieces_from_columns(&columns, b"A", b"CAC"),
                vec![
                    Piece {
                        ref_start: 0,
                        ref_end: 0,
                        alt: b"C".to_vec()
                    },
                    Piece {
                        ref_start: 1,
                        ref_end: 1,
                        alt: b"C".to_vec()
                    },
                ]
            );
        }

        #[test]
        fn a_retained_chunk_need_not_be_the_whole_reference() {
            // #1260's block with one flank base restored per side — the shape
            // any widening of a trimmed block produces inside a poly-A run.
            // The only all-matched decomposition is `a = 1, b = 2`:
            // `A` + `C` + `A` + `C` + `A`. A search fixed at `a == 0` finds
            // nothing here, and #1260 reverts to a spanning `delins`.
            assert_eq!(
                best_alignment(b"AAA", b"ACACA"),
                Some(vec![
                    (Some(0), Some(0)),
                    (None, Some(1)),
                    (Some(1), Some(2)),
                    (None, Some(3)),
                    (Some(2), Some(4)),
                ]),
                "every reference base must survive, with `a = 1`"
            );
            assert_eq!(
                partition_block(b"AAA", b"ACACA").len(),
                2,
                "and it must still partition into the two insertions"
            );
        }

        #[test]
        fn an_untrimmed_flank_on_one_side_only_is_also_expressible() {
            // The asymmetric case: a restored base on the 5' side alone, so
            // `a = 1` and `b == ref.len()`.
            assert_eq!(
                partition_block(b"AA", b"ACAC").len(),
                2,
                "a > 0 with the retained chunk running to the block's end"
            );
            // ...and on the 3' side alone, so `a == 0` and `b < ref.len()`.
            assert_eq!(
                partition_block(b"AA", b"CACA").len(),
                2,
                "b < ref.len() with the retained chunk starting at the block's start"
            );
        }

        #[test]
        fn two_adjacent_insertions_are_left_to_the_single_gap_search() {
            // `b == a` would put the two insertions next to each other in the
            // result, which is one contiguous gap. The single-gap search already
            // places it, so this must decline rather than re-express it as two
            // members that a coalesce would immediately merge back.
            assert_eq!(
                partition_block(b"AA", b"AACC").len(),
                1,
                "an insertion at one junction is one member, not two"
            );
        }

        #[test]
        fn the_decomposition_tie_is_broken_five_prime() {
            // `AA -> AAAA` admits (a,b) of (0,1), (0,2) and (1,2), all matching
            // both reference bases. Nothing in the spec reaches this, so the
            // smallest `a` wins — the 5'-most placement, as before.
            assert_eq!(
                two_gap_insertion_alignment(b"AA", b"AAAA"),
                Some(vec![
                    (None, Some(0)),
                    (Some(0), Some(1)),
                    (None, Some(2)),
                    (Some(1), Some(3)),
                ]),
                "a = 0, b = 1, |I1| = 1"
            );
        }

        #[test]
        fn a_long_shared_affix_does_not_hide_the_decomposition() {
            // #1260's UNTRIMMED window. Its only all-matched decomposition is
            // `a = 2, b = 3` — `AA` + `C` + `A` + `C` + `AAAA` — which needs a
            // common prefix of 2 and a common suffix of 4. An earlier revision
            // bounded cost by clamping those affix lengths to 2, which excluded
            // the sole solution and silently returned `None`.
            //
            // It was invisible: on this block the single-gap search reaches the
            // same member count by another route, so no corpus row moved. Only
            // a change elsewhere in the tie-break exposed it. Hence this asserts
            // the exact columns rather than a piece count.
            assert_eq!(
                two_gap_insertion_alignment(b"AAAAAAA", b"AACACAAAA"),
                Some(vec![
                    (Some(0), Some(0)),
                    (Some(1), Some(1)),
                    (None, Some(2)),
                    (Some(2), Some(3)),
                    (None, Some(4)),
                    (Some(3), Some(5)),
                    (Some(4), Some(6)),
                    (Some(5), Some(7)),
                    (Some(6), Some(8)),
                ]),
                "every reference base must survive, with a = 2 and b = 3"
            );
        }

        #[test]
        fn a_single_gap_that_already_keeps_the_whole_reference_is_left_alone() {
            // `AT -> ACT`: the 1 nt gap placed after `A` keeps both reference
            // bases, so there is nothing for a second gap to buy. One gap, and
            // the 5'-most placement, per `best_alignment`'s existing tie-break.
            assert_eq!(
                best_alignment(b"AT", b"ACT"),
                Some(vec![
                    (Some(0), Some(0)),
                    (None, Some(1)),
                    (Some(1), Some(2))
                ])
            );
        }

        #[test]
        fn a_split_at_one_unchanged_base_is_admitted() {
            // `general.md:34`: "two variants separated by one or more
            // nucleotides should be described individually and not as a
            // delins". One base is one or more. #1260's trimmed block reaches
            // this only because the two-gap alignment can express it.
            assert_eq!(
                partition_block(b"A", b"CAC").len(),
                2,
                "#1260a: two insertions either side of the retained A"
            );
            assert_eq!(
                partition_block(b"GC", b"CGA").len(),
                2,
                "#999-neg: the two-member count the corpus records as required"
            );
        }

        #[test]
        fn a_split_that_buys_no_higher_priority_type_collapses() {
            // #422. The same shape as #999-neg — net insertion, one
            // coincidentally matched interior base — but its split is two
            // `delins` members, so it buys no higher-priority type than the
            // spanning `delins` it came from — ferro policy, per
            // `split_buys_no_higher_priority_type`'s doc, not `general.md:56`.
            // #1235's own comment asks for
            // exactly this: a fix that resolves #422 and keeps #999 green.
            assert_eq!(partition_block(b"CTATAG", b"AAACCCC").len(), 1);
        }

        #[test]
        fn the_collapse_spares_an_equal_length_block() {
            // `long-delins-40nt`: both members are `delins`, but the block is
            // equal-length, so there is no gap to place and every matched base
            // is a genuine coordinate-wise identity rather than a placement
            // artefact. Collapsing here would lose a real two-member answer.
            assert_eq!(
                partition_block(
                    b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                    b"CATGCATGCATGCATGCATTCATGCATGCATGCATGCATG"
                )
                .len(),
                2
            );
        }

        #[test]
        fn the_collapse_spares_a_member_that_renders_as_inv_or_dup() {
            // A structural `delins` can still render as an `inv` — the one
            // renderer type that hides inside one — and `inv` outranks `delins`,
            // so such a split is worth keeping. Measured: a naive all-`delins`
            // test collapses 694 of 1094 firings over an exhaustive `AT` sweep.
            let pieces = partition_block(b"AAAAA", b"TATT");
            assert!(
                pieces.len() > 1,
                "a member that is the reverse complement of its own span must \
                 not be collapsed away: {pieces:?}"
            );
        }

        #[test]
        fn the_second_gap_is_scoped_out_of_cluster_a() {
            // Equal-length and net-deletion blocks never reach the two-gap
            // search, which is what keeps it clear of the "stays delins" corpus
            // (#1034/#1040/#182) that `merge.rs`'s single-gap restriction
            // exists to protect. Pinned as behaviour, not as a code path.
            assert_eq!(
                partition_block(b"CAGTGACTAG", b"TGTCACGACT").len(),
                1,
                "#1040"
            );
            assert_eq!(
                partition_block(b"GCT", b"AGC").len(),
                1,
                "#1034 whole-run revcomp"
            );
            assert_eq!(
                partition_block(b"AGTCAGT", b"GATTA").len(),
                3,
                "#1157A, net deletion"
            );
        }
    }
}
