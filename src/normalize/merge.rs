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
///     ins/del/sub/delins (literal payloads only);
///   * the del/sub/delins spans form one contiguous interval with no
///     unchanged interior base (no holes);
///   * every insertion attaches to that interval (its gap lies in
///     `[c_lo - 1, c_hi]`);
///   * the group has BOTH an insertion and a del/sub/delins.
///
/// Anything else passes through untouched, so pure-insertion groups,
/// pure-deletion groups (owned by `merge_consecutive_edits`), and
/// non-overlapping inputs are unaffected.
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
    // accession / body-region refusals, same per-`NaEdit` guards. It admits one
    // shape this path does not — `NaEdit::Inversion`, which the sequence-first
    // canonicalizer serves — so reject that explicitly and the two passes stay
    // one implementation instead of two copies ten screens apart.
    let Some(edits) = collect_canonical_edits(&variants, kind, body, &template_accession) else {
        return variants;
    };
    if edits.iter().any(|e| matches!(e, GEdit::Inv { .. })) {
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
            // `Inv` is unreachable here — the rejection above drops any group
            // carrying one — but it is a replacement over `[s, e]`, so classify
            // it as one rather than leaving a misleading `None`.
            GEdit::Del { s, e } | GEdit::Delins { s, e, .. } | GEdit::Inv { s, e } => {
                Some((*s, *e))
            }
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
fn cis_axis_parts(
    v: &HgvsVariant,
    kind: CisKind,
) -> Option<(&Accession, Region, i64, i64, &NaEdit)> {
    let (accession, region, s, e, edit) = match (kind, v) {
        (CisKind::Genome, HgvsVariant::Genome(g)) => {
            let (r, s, e) = simple_genome_range(&g.loc_edit.location)?;
            (&g.accession, r, s, e, &g.loc_edit.edit)
        }
        (CisKind::Mt, HgvsVariant::Mt(m)) => {
            let (r, s, e) = simple_genome_range(&m.loc_edit.location)?;
            (&m.accession, r, s, e, &m.loc_edit.edit)
        }
        (CisKind::Cds, HgvsVariant::Cds(c)) => {
            let (r, s, e) = simple_cds_range(&c.loc_edit.location)?;
            (&c.accession, r, s, e, &c.loc_edit.edit)
        }
        (CisKind::Tx, HgvsVariant::Tx(t)) => {
            let (r, s, e) = simple_tx_range(&t.loc_edit.location)?;
            (&t.accession, r, s, e, &t.loc_edit.edit)
        }
        (CisKind::Rna, HgvsVariant::Rna(rv)) => {
            let (r, s, e) = simple_rna_range(&rv.loc_edit.location)?;
            (&rv.accession, r, s, e, &rv.loc_edit.edit)
        }
        _ => return None,
    };
    if !edit.is_certain() {
        return None;
    }
    Some((accession, region, s, e, edit.inner()?))
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
    let edit = if merged.start > merged.end {
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

/// Longest changed block the canonicalizer will attempt to partition.
///
/// This is a cost bound, not a policy: the alignment is quadratic in the block
/// length, and a block this long is a structural event rather than a few
/// nucleotide changes. It is deliberately far above the size at which the
/// separation rule (`delins.md:17`) is the interesting question.
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
/// This bound therefore applies only where a real guard exists — see
/// [`MAX_UNGUARDED_SPLIT_BLOCK`] for the regime where one does not.
const MAX_SPLIT_BLOCK: usize = 1024;

/// The same bound for a **net deletion**, where nothing else guards the split.
///
/// [`separations_are_meaningful`] is restricted to net insertions, so the three
/// block regimes are not equally protected:
///
/// | regime | guard against a coincidental split |
/// | --- | --- |
/// | equal length | none needed — [`best_alignment`] compares position-wise, so there is no gap to place and no search to go wrong |
/// | net insertion | [`separations_are_meaningful`] |
/// | **net deletion** | **nothing** |
///
/// In the third regime this constant *is* the guard, which is why it stays at
/// the original 32. That was never the intent — it guarded by accident, being
/// smaller than the blocks at risk — but raising it to [`MAX_SPLIT_BLOCK`]
/// across the board removed the cover and split the spec's own
/// `delins.md:44-47` worked example:
///
/// ```text
/// LRG_199t1:c.850_901delinsTTCCTCGATGCCTG        52 nt -> 14 nt
///   spanning  LRG_199:g.646630_646681delins…
///   split     LRG_199:g.[646630_646636delinsTTCCTCG;646640_646678delinsC;
///                        646680_646681delinsTG]
/// ```
///
/// Which is precisely the harm that passage names, on that passage's own
/// example: the corpus shows one correct
/// `p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)` becoming three bogus
/// consequences.
///
/// Splitting the bound by regime keeps the raise where it is safe — a long
/// equal-length block has no alignment choice to get wrong, and a long net
/// insertion is guarded — while leaving net deletions exactly as they were.
///
/// The principled fix is to extend [`separations_are_meaningful`] to net
/// deletions and then retire this constant. That is not a local change: the
/// notes there record that a naive extension breaks #1232 and #1157 and breaks
/// confluence outright, so it needs its own measurement pass.
const MAX_UNGUARDED_SPLIT_BLOCK: usize = 32;

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
/// It applies only to blocks no larger than
/// [`MAX_SINGLE_BASE_SEPARATION_BLOCK`]; above that a single unchanged base is
/// not believable and the threshold stays at 2.
const MIN_PIECE_SEPARATION: usize = 1;

/// Largest block, in bases, for which one unchanged base counts as separation.
///
/// A single matched base is evidence of structure in a short block and evidence
/// of nothing in a long one, where a two- or three-base run recurs by chance.
/// The Mutalyzer conformance corpus is what forces this: a 6 nt block replaced
/// by a 21 nt payload (`NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC` and
/// nine siblings) has interior bases that happen to match, and the oracle keeps
/// it as one `delins`. Without a bound, a threshold of 1 shreds it into
/// `g.[5207_5209delinsGTC;5211C>T;5213_5214insGCTCATTATCTGGCT]` — mixed member
/// types, so [`split_buys_no_higher_priority_type`] does not reach it either.
///
/// **The value is under-determined and the bounds are what is measured.** The
/// corpus constrains it only to `(4, 21]`: `#1260b`'s block is 4 nt and must
/// split at one base, the Mutalyzer block is 21 nt and must not. Every value in
/// between scores identically, so `8` is a choice inside a measured window, not
/// a calibrated constant. Narrowing it needs a case in that gap, not a
/// re-derivation.
const MAX_SINGLE_BASE_SEPARATION_BLOCK: usize = 8;

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
    /// Whether this piece claims any offset inside one of `ranges` as changed.
    ///
    /// A pure insertion spans no reference base, so it overlaps nothing — which
    /// is exactly right for the separator veto: an insertion never swallows a
    /// base the input left between two members.
    fn overlaps_any(&self, ranges: &[(usize, usize)]) -> bool {
        ranges
            .iter()
            .any(|&(lo, hi)| self.ref_start < hi && lo < self.ref_end)
    }

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
}

/// Whether to run both block splitters and report their disagreements.
///
/// Enabled by `FERRO_SEQFIRST_SHADOW=1`. Read once and cached, like the
/// idempotency gate in [`crate::normalize`], so a disabled run pays one relaxed
/// atomic load per canonicalization. When on, the sequence-first splitter runs
/// *in addition to* the live one and any disagreement goes to stderr; it never
/// changes what is returned.
fn seqfirst_shadow_enabled() -> bool {
    static ENABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var("FERRO_SEQFIRST_SHADOW").as_deref() == Ok("1"))
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
    //   (`general.md:44`) — covered by the pass's existing refusals: the body
    //   region check in `collect_canonical_edits` keeps every member inside the
    //   positive body, and the lone-pure-insertion bail hands terminal
    //   insertions back to the per-member pipeline that owns those clamps.
    let body = body_region(kind);
    let (template_accession, _, _, _, _) = cis_axis_parts(&variants[0], kind)?;
    let template_accession = template_accession.clone();

    let edits = collect_canonical_edits(variants, kind, body, &template_accession)?;

    // Window: the union of every member's footprint, padded for the 3'-shift.
    let (c_lo, c_hi) = edit_span_union(&edits)?;
    let w_lo = (c_lo - CANONICAL_PAD).max(1);
    let w_hi = c_hi + CANONICAL_PAD;
    if w_hi - w_lo + 1 > MAX_CANONICAL_WINDOW {
        return None;
    }

    let frame = axis_frame(kind, &template_accession, provider)?;
    let accession = template_accession.transcript_accession();
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

    let mut pieces = partition_block(&ref_bytes[lo..hi_ref], &result[lo..hi_alt]);
    // Shadow the sequence-first splitter on the very same trimmed block, before
    // any 3'-shift or coalescing, so a reported disagreement is a disagreement
    // about the *split* and not about a downstream pass. Never affects `pieces`.
    if seqfirst_shadow_enabled() {
        // Reading-frame axes get the coding one-amino-acid exception
        // (`seqfirst::MIN_SEPARATION`, 2); every other axis gets the general
        // rule (`MIN_SEPARATION_NO_FRAME`, 1). Same distinction
        // `apply_coding_codon_exception` uses below, via the same `AxisFrame`.
        let min_separation = if frame.reading_frame {
            crate::normalize::seqfirst::MIN_SEPARATION
        } else {
            MIN_SEPARATION_NO_FRAME
        };
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
        let mut min_live = pieces.clone();
        shrink_pieces_to_differences(&mut min_live, block);
        let mut min_shadow = shadow.clone();
        if let Some(min_shadow) = min_shadow.as_mut() {
            shrink_pieces_to_differences(min_shadow, block);
        }
        if shadow.as_ref() == Some(&pieces) {
            // The denominator. Without it, `grep -c SEQFIRST_SHADOW` returning 0
            // cannot be told apart from a shadow that never ran — the failure
            // mode that would make an audit of the disagreements vacuous.
            eprintln!("SEQFIRST_AGREE");
        } else if min_shadow.as_ref() == Some(&min_live) {
            eprintln!(
                "SEQFIRST_MINAGREE ref={} alt={} live={:?} shadow={:?}",
                String::from_utf8_lossy(block),
                String::from_utf8_lossy(&result[lo..hi_alt]),
                DebugPieces(&pieces),
                shadow.as_deref().map(DebugPieces),
            );
        } else {
            eprintln!(
                // Space-separated, and no space inside a piece list, on purpose:
                // `cargo nextest` strips tab characters when it re-emits captured
                // test output (`--success-output immediate`), so a tab-delimited
                // report is unparseable in the one run mode that covers the whole
                // suite. Verified by `od -c` on both run modes.
                "SEQFIRST_SHADOW ref={} alt={} live={:?} shadow={:?} min_live={:?} min_shadow={:?}",
                String::from_utf8_lossy(block),
                String::from_utf8_lossy(&result[lo..hi_alt]),
                DebugPieces(&pieces),
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
    shift_pieces_three_prime(&mut pieces, &ref_bytes);
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

    // Applied *after* the weight bound, deliberately. The bound is a statement
    // about the re-derived partition — it may not describe more change than the
    // input did — and the codon exception (`general.md:35`) is a licensed
    // widening on top of an already-accepted partition, not part of what the
    // bound judges. Running it first would let a legitimate codon merge inflate
    // the weight and trip a refusal.
    apply_coding_codon_exception(&mut pieces, frame.reading_frame, w_lo, &ref_bytes);
    if needs_unsupported_form(&pieces, &ref_bytes) {
        return None;
    }
    // Re-partitioning is the point of this pass, but *merging* two members
    // across a base the input left between them is not: two variants separated
    // by one or more unchanged nucleotides are described individually
    // (`general.md:34`). Only a shift that closes the gap may merge them, and
    // that shows up as pieces still outnumbering nothing — so the check applies
    // exactly when the derivation has produced fewer pieces than there were
    // members.
    if pieces.len() < edits.len() {
        let gaps = input_separator_gaps(&edits, w_lo);
        if pieces.iter().any(|piece| piece.overlaps_any(&gaps)) {
            return None;
        }
    }
    // A derivation that collapses to a single pure insertion belongs to the
    // per-member pipeline, which owns the terminal-insertion clamps (#1205,
    // #1217) and duplication recovery. Re-deriving one here would undo a clamp
    // the pipeline had already applied.
    if pieces.len() == 1 && pieces[0].is_pure_indel() && !pieces[0].alt.is_empty() {
        return None;
    }

    let rebuilt = rebuild_members(&pieces, &variants[0], body, w_lo)?;
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
/// sequence does not: repeats (the spec says a repeat and a del/dup spelling of
/// one change are *both* correct and declines to choose — `open-issues.md:88`),
/// conversions, methylation, copy number, identity/`=`, and uncertain edits.
/// Also refuses a mixed axis, a second accession, and any member outside the
/// positive body region (a contiguous window would happily produce `c.-1_1`,
/// which is not valid HGVS).
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
            _ => return None,
        }
    }
    Some(edits)
}

/// Whether every member's *stated* reference base agrees with the window.
///
/// `NaEdit` optionally carries the reference bases the submitter asserted
/// (`c.5A>T`, `c.5_7delACG`). Where present they must match, otherwise the
/// variant is a reference mismatch and belongs to the per-member pipeline's
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
            NaEdit::Deletion {
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
            | GEdit::Dup { s, e } => (*s, *e),
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
            GEdit::Ins { .. } | GEdit::Dup { .. } | GEdit::Sub { .. } => continue,
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
            GEdit::Del { .. } | GEdit::Delins { .. } | GEdit::Inv { .. } | GEdit::Sub { .. } => {
                continue
            }
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
    // The cap depends on whether anything else guards this block against a
    // coincidental split. See `MAX_UNGUARDED_SPLIT_BLOCK`.
    let cap = if result.len() < reference.len() {
        MAX_UNGUARDED_SPLIT_BLOCK
    } else {
        MAX_SPLIT_BLOCK
    };
    if reference.len() > cap || result.len() > cap {
        return whole();
    }
    let Some(columns) = best_alignment(reference, result) else {
        return whole();
    };
    let pieces = pieces_from_columns(&columns, reference, result);
    if pieces.is_empty() {
        return whole();
    }
    if result.len() > reference.len()
        && !separations_are_meaningful(&pieces, reference.len().max(result.len()))
    {
        return whole();
    }
    // Scoped to length-changing blocks: an equal-length block has no gap to
    // place, so every matched base is a genuine coordinate-wise identity rather
    // than an artefact of where the gap landed, and a split across one is real.
    // `MAX_UNGUARDED_SPLIT_BLOCK`'s own doc draws the same line.
    if reference.len() != result.len()
        && pieces.len() > 1
        && every_separation_is_a_single_base(&pieces)
        && split_buys_no_higher_priority_type(&pieces, reference)
    {
        return whole();
    }
    pieces
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
            GEdit::Delins { s, e, alt } => span(*s, *e).max(alt.len()),
        })
        .sum()
}

/// The same measure for a derived piece set (see `changed_columns_of_edits`).
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
/// Two properties were measured rather than assumed, and both pin the shape:
///
/// * **Restricted to net insertions.** Extending it to net deletions breaks
///   #1232 and #1157, which split a shorter payload at a genuinely unchanged
///   base, and breaks confluence outright: a two-member input collapsing to one
///   piece re-arms the `input_separator_positions` veto and returns verbatim,
///   while the spanning spelling stays spanning — two stable strings for one
///   variant.
/// * **Measured before the 3'-shift.** Re-checking afterwards lets every one of
///   the corpus's coincidental splits back through, because a shift widens the
///   gap to a piece's left neighbour.
///
/// Alternatives that look more principled do not survive measurement. Gating on
/// `best_alignment`'s score margin cannot work at all: in these cases the winning
/// placement *is* the 5'-most one, so the margin is zero. Match-density variants
/// trade corpus rows against unit tests without a threshold that satisfies both.
///
/// # Net deletions are not covered here
///
/// Deliberately, per the first bullet above. A long net deletion is instead held
/// by [`MAX_UNGUARDED_SPLIT_BLOCK`], which stays at the original 32 precisely
/// because this rule does not reach it — see that constant for the spec example
/// that fails when the bound is raised across the board. Extending this rule to
/// cover them, and retiring that constant, is the principled fix and needs its
/// own measurement pass.
fn separations_are_meaningful(pieces: &[Piece], block_len: usize) -> bool {
    // `ref_end` is exclusive and a pure insertion has `ref_start == ref_end`, so
    // this already counts unchanged *bases* rather than event indices: a
    // junction contributes no width because it occupies none. The sequence-first
    // splitter had to be corrected for exactly this, because it works in event
    // space where a junction and a changed position are both one offset.
    let required = if block_len > MAX_SINGLE_BASE_SEPARATION_BLOCK {
        2
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

    let mut best: Option<(usize, usize)> = None;
    for k in 0..=anchored {
        if k > 0 {
            leading += usize::from(reference[k - 1] == result[k - 1]);
            trailing -= usize::from(ungapped(k - 1));
        }
        let matches = leading + trailing;
        // Strictly greater keeps the first (5'-most) alignment on a tie.
        if best.as_ref().is_none_or(|(score, _)| matches > *score) {
            best = Some((matches, k));
        }
    }
    // A single gap cannot express *insertion, retained reference, insertion*,
    // so consider the two-gap form when the best single gap leaves a reference
    // base unmatched. Scoped to net insertions by `shorter_is_ref`, which is
    // what keeps it clear of the corpus the single-gap restriction protects —
    // see [`two_gap_insertion_alignment`].
    if shorter_is_ref && best.is_some_and(|(matches, _)| matches < reference.len()) {
        if let Some(columns) = two_gap_insertion_alignment(reference, result) {
            return Some(columns);
        }
    }

    // Only the winner is materialized.
    best.map(|(_, k)| {
        single_gap_alignment(k, gap_len, shorter_is_ref, reference.len(), result.len())
    })
}

/// The alignment that keeps the **whole** reference intact between two
/// insertions, or `None` when the pair does not have that shape (#1260).
///
/// PR #1285 named the gap this fills: for #1260 the answer is "a two-gap
/// alignment — insertion, retained base, insertion — which `best_alignment`'s
/// single-gap restriction cannot express". One contiguous indel must either
/// swallow the retained base into a spanning `delins` or leave one insertion
/// unaccounted for.
///
/// # Why one shape is the whole search
///
/// The caller trims common flanks before partitioning
/// (`canonicalize_from_sequence` via [`trim_common_flanks`]), and on a trimmed
/// block the general three-chunk form collapses to this one. Writing an
/// all-matched two-gap alignment as `result = ref[..a] + I1 + ref[a..b] + I2 +
/// ref[b..]`: if `a > 0` the two blocks share a first base, and if `b <
/// ref.len()` they share a last one — both contradict trimming. So `a == 0` and
/// `b == ref.len()`, leaving `result = I1 + reference + I2`. That is a substring
/// search, not a quadratic sweep over gap placements.
///
/// # Why this does not reopen the corpus the single-gap rule protects
///
/// `merge.rs`'s single-gap restriction exists to stop *compensating* gaps — a
/// deletion **and** an insertion in one block — from manufacturing matches out
/// of coincidence and shredding a contiguous replacement (#1034/#1040/#182).
/// Both gaps here are insertions in a net-insertion block, so no deletion is
/// introduced and no match is manufactured: every matched column is a reference
/// base surviving verbatim. Cluster A is structurally out of reach — #1040 and
/// the `inv` cases are equal-length, #1157A is a net deletion, and an
/// equal-length or net-deletion block never reaches this function.
///
/// Both insertions must be non-empty; otherwise this is a single gap and the
/// existing search already places it, 5'-most.
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
    // 5'-most placement, matching `best_alignment`'s tie-break. The bounds put
    // at least one base on each side of the retained reference.
    let last_start = result.len() - reference.len() - 1;
    let start = (1..=last_start).find(|&p| result[p..p + reference.len()] == *reference)?;

    let mut columns = Vec::with_capacity(result.len());
    columns.extend((0..start).map(|a| (None, Some(a))));
    columns.extend((0..reference.len()).map(|i| (Some(i), Some(start + i))));
    columns.extend((start + reference.len()..result.len()).map(|a| (None, Some(a))));
    Some(columns)
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

/// 3'-shift each length-changing piece, clamped so it cannot reach a sibling.
///
/// The clamp is the #1234 fix. The 3' rule (`general.md:41`) is stated per
/// description with no allele-awareness, and the only text that ever addressed
/// cross-member shifting was the *rejected* SVD-WG010. Left unbounded, a
/// deletion shifts over a downstream substitution and the two members overlap —
/// malformed, and denoting a different sequence. Substitutions are anchored and
/// never shift.
fn shift_pieces_three_prime(pieces: &mut [Piece], ref_bytes: &[u8]) {
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
            ShuffleDirection::ThreePrime,
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
            // A 3'-shuffle only ever advances: it initialises its cursor to
            // `start` and mutates it solely via `new_start += 1`. Unlike the
            // `mod.rs` precedent, which serves both directions and once
            // regressed 5' insertions by clamping leftward moves with
            // `saturating_sub`, this call site is `ThreePrime`-only, so a
            // leftward move would be a bug in `shuffle`, not a case to handle.
            debug_assert!(
                new_start >= old_start,
                "a 3'-shuffle cannot move a piece 5' (old={old_start}, new={new_start})",
            );
            let rotation = new_start.saturating_sub(old_start) % pieces[i].alt.len();
            if rotation > 0 {
                pieces[i].alt.rotate_left(rotation);
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
/// `partition_block` sibling never reaches `best_alignment` in the first place;
/// see [`MAX_UNGUARDED_SPLIT_BLOCK`]), not this exception, so it must not be
/// touched here and is not: 52 ≠ 3.
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
        let is_merged_triplet = piece.ref_end == piece.ref_start + 3
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
        if pieces[i - 1].ref_end >= pieces[i].ref_start {
            let next = pieces.remove(i);
            let prev = &mut pieces[i - 1];
            prev.ref_end = prev.ref_end.max(next.ref_end);
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
) -> Option<Vec<HgvsVariant>> {
    let mut members = Vec::with_capacity(pieces.len());
    for piece in pieces {
        members.push(build_merged(template, anchor_for_piece(piece, body, w_lo)?));
    }
    Some(members)
}

/// Build the [`Anchor`] for one piece, typing it under the spec's rules.
///
/// [`build_naedit`] already picks insertion / deletion / delins from the span
/// and replacement, which covers the mandates that 1 nt → 1 nt is a
/// substitution (`delins.md:15`) and that consecutive changes are one delins
/// (`delins.md:16`). Inversion and tandem duplication are *not* built here —
/// [`needs_unsupported_form`] refuses the whole group when a piece would need
/// either, leaving them to the per-member pipeline.
fn anchor_for_piece(piece: &Piece, body: Region, w_lo: i64) -> Option<Anchor> {
    let alt: Vec<Base> = piece
        .alt
        .iter()
        .filter_map(|b| Base::from_char(*b as char))
        .collect();
    if alt.len() != piece.alt.len() {
        return None; // non-IUPAC byte; refuse rather than mangle.
    }
    let start = w_lo + piece.ref_start as i64;
    let end = w_lo + piece.ref_end as i64 - 1;
    Some(Anchor {
        region: body,
        start,
        end,
        alt,
    })
}

/// Reference runs the input left unchanged *between* two of its members, as
/// half-open 0-based offset ranges into the window.
///
/// Returned as ranges rather than as the individual positions they contain: two
/// members at opposite ends of a 4096 nt window are separated by thousands of
/// unchanged bases, and materializing every one of them — then rescanning the
/// list once per derived piece — is quadratic work for what is an interval
/// overlap test.
///
/// Two variants separated by one or more unchanged nucleotides are described
/// individually and **not** as one delins (`general.md:34`, restated in every
/// DNA page). So a re-derivation may split a member and may merge members the
/// 3'-shift made adjacent, but it must never swallow a base the input itself
/// held between two members: `[1006_1008del;1009_1010insCCC]` leaves 1009
/// untouched, and collapsing it to `1006_1009delinsTCCC` merges across it
/// (#180).
///
/// An insertion consumes no reference base, so its footprint is the junction it
/// sits in rather than a span.
fn input_separator_gaps(edits: &[GEdit], w_lo: i64) -> Vec<(usize, usize)> {
    let mut footprints: Vec<(i64, i64)> = edits
        .iter()
        .map(|edit| match edit {
            // An insertion between `gap` and `gap + 1` touches neither base;
            // represent it as the empty span just past `gap`.
            GEdit::Ins { gap, .. } => (*gap + 1, *gap),
            GEdit::Dup { s: _, e } => (*e + 1, *e),
            GEdit::Del { s, e } | GEdit::Delins { s, e, .. } | GEdit::Inv { s, e } => (*s, *e),
            GEdit::Sub { pos, .. } => (*pos, *pos),
        })
        .collect();
    if footprints.len() < 2 {
        return Vec::new();
    }
    footprints.sort_unstable();
    let mut gaps = Vec::new();
    for pair in footprints.windows(2) {
        let (_, previous_end) = pair[0];
        let (next_start, _) = pair[1];
        // Half-open `[previous_end + 1, next_start)` in axis coordinates,
        // clamped to the window's own 0-based offsets.
        let lo = (previous_end + 1 - w_lo).max(0);
        let hi = next_start - w_lo;
        if hi > lo {
            gaps.push((lo as usize, hi as usize));
        }
    }
    gaps
}

/// Whether any piece would have to be rendered as an inversion or a tandem
/// duplication.
///
/// [`build_naedit`] only produces insertion / deletion / delins / identity, and
/// the spec is emphatic that these two forms are not interchangeable with those:
/// an inversion is the reverse complement of its span (`inversion.md:5`), and a
/// tandem duplication **must** be described as a `dup`, never an insertion
/// (`duplication.md:18`). Rather than teach the builder two more shapes, refuse
/// the group and let the existing per-member pipeline — which already gets both
/// right — handle it. Refusing costs nothing here: these are exactly the cases
/// where a single member is already canonical.
fn needs_unsupported_form(pieces: &[Piece], ref_bytes: &[u8]) -> bool {
    pieces
        .iter()
        .any(|piece| is_inversion(piece, ref_bytes) || is_tandem_duplication(piece, ref_bytes))
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
    region: Region,
    start: i64,
    end: i64,
    /// Whether the edit consumes the reference bases under its span.
    claims_bases: bool,
    /// Whether the edit must stop a sibling's shift at its boundary. Superset of
    /// `claims_bases` — see `blocks_sibling_shift`.
    blocks_shift: bool,
    /// For an insertion-like edit, the position its added sequence sits
    /// immediately 3' of; `None` for anything that consumes bases instead.
    junction: Option<i64>,
}

/// Read a member's span for the clamp pass, or `None` for a member this pass
/// cannot place (wrong kind, uncertain edit, offset / special position).
fn member_span(v: &HgvsVariant, kind: CisKind) -> Option<MemberSpan> {
    let (accession, region, start, end, edit) = cis_axis_parts(v, kind)?;
    Some(MemberSpan {
        accession: accession.full(),
        provider_key: accession.transcript_accession(),
        region,
        start,
        end,
        claims_bases: claims_reference_bases(edit),
        blocks_shift: blocks_sibling_shift(edit),
        junction: junction_of(edit, start, end),
    })
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
    let pre: Vec<Option<MemberSpan>> = before.iter().map(|v| member_span(v, kind)).collect();
    let post: Vec<Option<MemberSpan>> = after.iter().map(|v| member_span(v, kind)).collect();

    for i in 0..after.len() {
        let (Some(b), Some(a)) = (pre[i].as_ref(), post[i].as_ref()) else {
            continue;
        };
        if !b.blocks_shift || !a.blocks_shift || b.region != a.region {
            continue;
        }
        // A member that did not move cannot have crossed anything.
        let delta = a.start - b.start;
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
        let translated = a.end - b.end == delta;
        let shrank = a.end - a.start < b.end - b.start;
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
        let siblings: Vec<&MemberSpan> = (0..pre.len())
            .filter(|&j| j != i)
            .flat_map(|j| [pre[j].as_ref(), post[j].as_ref()])
            .flatten()
            .filter(|s| s.region == a.region && s.accession == a.accession)
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
                .filter(|&start| start > b.end && start <= a.end)
                .map(|start| start - 1);
            // A junction between `j` and `j+1` is crossed once this member's
            // end reaches `j + 1`; stopping at `j` leaves it flush against the
            // junction, which is the adjacency the collapse pass wants.
            let across_junctions = siblings
                .iter()
                .filter_map(|s| s.junction)
                .filter(|&j| j >= b.end && j < a.end);
            onto_bases
                .chain(across_junctions)
                .min()
                .map(|limit| a.end - limit)
        } else {
            // 5'-shift: mirror image, window `[a.start, b.end]`.
            let onto_bases = siblings
                .iter()
                .filter(|s| s.blocks_shift)
                .map(|s| s.end)
                .filter(|&end| end < b.start && end >= a.start)
                .map(|end| end + 1);
            // No junction barrier in this direction: measured, it introduced
            // fresh 5'-shuffle failures rather than removing them (an
            // insertion-like sibling and a 5'-shifting span interact through
            // the merge's cancellation path, not through the crossing this
            // rule models). 3'-only, deliberately — see the module tests.
            onto_bases.max().map(|limit| a.start - limit)
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
            // A duplication blocks a sibling's shift without claiming bases, so
            // it reaches this pass too — and its payload is read from under its
            // own span, so it cannot simply be slid (#1280, #1292).
            match a.junction {
                Some(_) => translate_junction_member(&mut after[i], -pull, kind, a, provider),
                None => translate_member(&mut after[i], -pull, kind, a),
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
pub(crate) fn demote_repeats_spanning_siblings(
    before: &[HgvsVariant],
    after: &mut [HgvsVariant],
    phase: AllelePhase,
    uncertain: bool,
) {
    if phase != AllelePhase::Cis || uncertain || after.len() < 2 || before.len() != after.len() {
        return;
    }
    let Some(kind) = cis_kind_of(&after[0]) else {
        return;
    };
    let pre: Vec<Option<MemberSpan>> = before.iter().map(|v| member_span(v, kind)).collect();
    let post: Vec<Option<MemberSpan>> = after.iter().map(|v| member_span(v, kind)).collect();

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
        // Either way the equivalent edit sits on the 3'-most `length` bases of
        // the tract: removing them, or duplicating them.
        let start = a.end - length + 1;
        if length <= 0 || start < a.start || b.region != a.region {
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
        respell_as(&mut after[i], kind, a.region, start, a.end, source);
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
    let spans: Vec<Option<MemberSpan>> = members.iter().map(|v| member_span(v, kind)).collect();
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
                TerminalRespell::Allowed,
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
    let pre: Vec<Option<MemberSpan>> = before.iter().map(|v| member_span(v, kind)).collect();
    let post: Vec<Option<MemberSpan>> = after.iter().map(|v| member_span(v, kind)).collect();
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
        let (Some(before_junction), Some(after_junction)) = (
            cis_axis_parts(&before[i], kind).and_then(|(_, _, s, e, edit)| junction_of(edit, s, e)),
            cis_axis_parts(&after[i], kind).and_then(|(_, _, s, e, edit)| junction_of(edit, s, e)),
        ) else {
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

        // 3' only, and the 5' case is genuinely unsolved rather than merely
        // unwritten. A duplication's span and its junction move together under a
        // 5' shift, so bounding the junction there mis-places the copy: a mirror
        // of the rule below was measured to turn 80 previously-correct outputs
        // into silently wrong ones, and was removed as defective rather than
        // deferred.
        //
        // A 5' junction *does* cross, in a shape this rule would not have caught
        // anyway — when the sibling is **upstream** and the junction travels
        // toward it. The duplication half of that (`g.[4A>G;9dup]`, which used
        // to emit `g.4_5insG`) is closed by `blocks_sibling_shift`, which bounds
        // the duplication's *span* rather than its junction. The insertion half
        // remains open, narrowed and characterised in
        // `a_five_prime_insertion_junction_with_an_upstream_sibling_is_a_known_gap`.
        if after_junction < before_junction {
            continue;
        }
        // The junction crossed a base-claiming sibling if it started 5' of the
        // sibling's first base and ended at or past it.
        let onto_bases = siblings
            .iter()
            .filter(|s| s.start > before_junction && s.start <= after_junction)
            .map(|s| s.start - 1);
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
        // Every junction compared below sits on `a`'s own region axis, so one
        // offset serves them all (#1284). A region with no conversion — an
        // `n.` position outside the transcript — has no bases to rotate
        // against, and refusing here is the same conservative answer the
        // payload reads below give.
        let Some(axis_delta) = region_sequence_delta(a.region, &a.provider_key, provider) else {
            continue;
        };
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
                if junction <= before_junction || junction > after_junction {
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
                                axis_delta,
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
                    false => Some(junction - 1),
                }
            });
        let limit = onto_bases.chain(across_junctions).min();
        let Some(limit) = limit else {
            continue;
        };
        // Never move past where the junction started, and never onto a junction
        // this member's payload is out of phase with: the walk back finds the
        // most 3' position it can legally occupy. `before_junction` always can,
        // since that is where the shift began.
        let ceiling = limit.max(before_junction).min(after_junction);
        let Some(payload) = payload.as_ref() else {
            continue;
        };
        let Some(destination) = (before_junction..=ceiling).rev().find(|&j| {
            payload_at_junction(
                payload,
                after_junction,
                j,
                &a.provider_key,
                axis_delta,
                provider,
            )
            .is_some()
        }) else {
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
/// Reverts on any axis this pass cannot rewrite, or if the result does not land
/// exactly where intended on the same region.
fn respell_as(
    variant: &mut HgvsVariant,
    kind: CisKind,
    region: Region,
    start: i64,
    end: i64,
    source: RepeatSource,
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
    let landed = placed
        && member_span(variant, kind)
            .is_some_and(|s| s.region == region && s.start == start && s.end == end);
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
fn translate_member(variant: &mut HgvsVariant, delta: i64, kind: CisKind, from: &MemberSpan) {
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
        && member_span(variant, kind).is_some_and(|s| {
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
    let spans: Vec<Option<MemberSpan>> = members.iter().map(|v| member_span(v, kind)).collect();

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
        // inclusive is `[start - 1, end)` half-open 0-based, after `delta` puts
        // the member's region axis onto the served sequence (#1284).
        let Some(delta) = region_sequence_delta(dup.region, &dup.provider_key, provider) else {
            continue;
        };
        // Checked, matching the two sibling conversions (`cds_relative_gap` and
        // `payload_at_junction`). `dup.start`/`dup.end` come off a parsed
        // description and `delta` off the record's CDS bounds, so neither is
        // bounded by anything this function controls; an unchecked `i64` add
        // panics in a debug build where the refusal one line below is the answer
        // every other unrepresentable coordinate here already gets.
        let (Some(from_axis), Some(to_axis)) = (
            dup.start.checked_sub(1).and_then(|s| s.checked_add(delta)),
            dup.end.checked_add(delta),
        ) else {
            continue;
        };
        let (Ok(from), Ok(to)) = (u64::try_from(from_axis), u64::try_from(to_axis)) else {
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
        let sibling_claims_the_landing_base = siblings()
            .filter(|s| s.claims_bases)
            .any(|s| s.start <= dup.end && s.end >= dup.end);
        let terminal = if sibling_claims_the_landing_base {
            TerminalRespell::Refused
        } else {
            TerminalRespell::Allowed
        };
        respell_at_gap(&mut members[i], dup, dup.end, edit, provider, terminal);
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
    let spans: Vec<Option<MemberSpan>> = members.iter().map(|v| member_span(v, kind)).collect();

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
            .map(|v| member_span(v, kind).map(|s| (s.start, s.end)))
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
            TerminalRespell::Allowed,
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
fn payload_at_junction<P: ReferenceProvider>(
    payload: &[Base],
    from: i64,
    to: i64,
    key: &str,
    delta: i64,
    provider: &P,
) -> Option<Vec<Base>> {
    // `from`/`to` are junctions on the member's own region axis; `delta` puts
    // the base they step over onto the sequence the provider serves.
    let base_at = |position: i64| -> Option<Base> {
        let position = position.checked_add(delta)?;
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
    let destination = junction + delta;
    let Some(axis_delta) = region_sequence_delta(from.region, &from.provider_key, provider) else {
        return;
    };
    let Some(moved) = payload_at_junction(
        &payload,
        junction,
        destination,
        &from.provider_key,
        axis_delta,
        provider,
    ) else {
        return;
    };

    let original = variant.clone();
    translate_member(variant, delta, kind, from);
    if *variant == original {
        return; // the translation was refused; nothing to repair
    }
    let landed = member_span(variant, kind)
        .filter(|s| s.junction == Some(destination))
        .and_then(|s| junction_payload(variant, kind, &s, provider));
    if landed.as_deref() == Some(moved.as_slice()) {
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
        TerminalRespell::Allowed,
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
/// Only the identity member is dropped, and only when a sibling's span *covers*
/// it under [`blocks_sibling_shift`] — a sibling that merely adds sequence at a
/// junction inside it contradicts nothing.
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
pub(crate) fn drop_identity_members_covered_by_siblings(
    members: &mut Vec<HgvsVariant>,
    phase: AllelePhase,
    uncertain: bool,
) {
    if phase != AllelePhase::Cis || uncertain || members.len() < 2 {
        return;
    }
    let Some(kind) = cis_kind_of(&members[0]) else {
        return;
    };
    let spans: Vec<Option<MemberSpan>> = members.iter().map(|v| member_span(v, kind)).collect();
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
                        && s.start <= span.start
                        && s.end >= span.end
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
            // The span is on the member's own region axis; shift it onto the
            // served sequence before reading (#1284).
            let delta = region_sequence_delta(span.region, &span.provider_key, provider)?;
            // Checked, for the same reason as the sibling conversions in
            // `respell_colliding_duplications` and `cds_relative_gap`: the span
            // comes off a parsed description and `delta` off the record's CDS
            // bounds, so an unchecked `i64` add would panic in a debug build
            // where every other unrepresentable coordinate here declines.
            let (from, to) = (
                u64::try_from(span.start.checked_sub(1)?.checked_add(delta)?).ok()?,
                u64::try_from(span.end.checked_add(delta)?).ok()?,
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

/// Both endpoints of a member's location as `(region, axis position)`, **not**
/// requiring the two to share a region.
///
/// [`member_span`] folds them through `join_pos`, which refuses a
/// region-spanning range because the pass's own adjacency arithmetic is
/// region-local. That refusal is right for a *span* and wrong for verifying a
/// *gap*: an interbase point on the CDS start is legitimately named `c.-1_1`,
/// and reading it back through `member_span` would report `None` and force a
/// correct repair to revert. Used only by [`respell_at_gap`].
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
            let Some((region, axis)) = cds_axis_end(span, provider, raw_length, Region::Rna) else {
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
    let written = match variant {
        HgvsVariant::Cds(_) => cds_axis_end(span, provider, raw_length, Region::Cds),
        HgvsVariant::Rna(_) => cds_axis_end(span, provider, raw_length, Region::Rna),
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
/// one position past it for an insertion, whose span is the gap itself.
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
/// regions, and `member_span` refuses those by construction.
/// Whether [`respell_at_gap`] may apply the terminal boundary re-spelling
/// (#1327) when the gap runs one past the last base.
///
/// `Refused` when a sibling already claims that base: see the comment at the
/// gate for why the two spellings are not interchangeable there (#1344).
#[derive(Clone, Copy, PartialEq, Eq)]
enum TerminalRespell {
    Allowed,
    Refused,
}

fn respell_at_gap<P: ReferenceProvider>(
    variant: &mut HgvsVariant,
    span: &MemberSpan,
    junction: i64,
    edit: NaEdit,
    provider: &P,
    terminal_respell: TerminalRespell,
) {
    let original = variant.clone();
    let (gap_start, gap_end) = (junction, junction + 1);

    // Where each end of the gap lands, in `(region, axis)` form. Genomic and
    // `n.` gaps stay in the member's own region by construction; a
    // CDS-relative one is resolved through the transcript.
    let cds_relative_gap = |body: Region| -> Option<((Region, i64), (Region, i64))> {
        let (cds_start, cds_end) = ordered_cds_bounds(&span.provider_key, provider)?;
        let delta = region_sequence_delta(span.region, &span.provider_key, provider)?;
        let left = junction.checked_add(delta)?;
        let right = left.checked_add(1)?;
        Some((
            cds_axis_position(left, cds_start, cds_end, body)?,
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
    // Skipped entirely when a sibling already claims the last base (#1344). The
    // boundary identity trades a zero-width junction insertion for a `delins`
    // that *claims* that base, and those two are not interchangeable when
    // something else is already there: an insertion flush against a sibling's
    // deleted base is the #999 adjacency the collapse pass merges, while two
    // members claiming one base is an overlap it cannot. Re-spelling therefore
    // converted a mergeable pair into an unmergeable one, and the out-of-range
    // coordinate it was avoiding never reached the output in that case anyway —
    // the merge consumed it. Measured: `c.[*10del;*11dup]` gave
    // `c.[*11del;*11delinsAA]`, which ferro's own strict mode rejects as
    // `OverlapConflictingEdits / W5002`, where the junction form settles on the
    // correct `c.*11=`.
    if terminal_respell == TerminalRespell::Allowed {
        if let Some(delta) = region_sequence_delta(span.region, &span.provider_key, provider) {
            let Some(last) = junction.checked_add(delta) else {
                return;
            };
            if let Ok(length) = provider.get_sequence_length(&span.provider_key) {
                if u64::try_from(last) == Ok(length) {
                    respell_at_sequence_end(variant, span, &edit, delta, length, provider);
                    return;
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
        HgvsVariant::Rna(r) => match cds_relative_gap(Region::Rna) {
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

    /// #1135: an insertion-shaped anchor with nothing left to insert (a
    /// cancelling del+ins after 3'-shifting) must not become an `Insertion`
    /// carrying an empty sequence — that is not valid HGVS and made
    /// `normalize_na_edit` divide by zero rotating the empty sequence.
    ///
    /// Pinned on `build_naedit` with the anchor state directly: reaching it
    /// through `merge_consecutive_edits` requires the normalizer to have
    /// shifted the pair together first, so a merge-level test exercises the
    /// `Delins` branch instead and would not cover this.
    #[test]
    fn empty_insertion_anchor_builds_an_identity() {
        // Insertion-shaped anchor (`start == end + 1`) with an empty payload.
        let anchor = Anchor {
            region: Region::Genome,
            start: 1010,
            end: 1009,
            alt: Vec::new(),
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
        };
        let (_, edit) = build_naedit(anchor, |_, p| GenomePos::new(p as u64));
        match edit {
            NaEdit::Insertion { .. } => {}
            other => panic!("expected NaEdit::Insertion, got {other:?}"),
        }
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

    /// `partition_block` bounds by *regime*, not by length alone.
    ///
    /// Only net deletions are unguarded (`separations_are_meaningful` covers net
    /// insertions, and an equal-length block has no alignment choice to make), so
    /// only they keep the original 32 nt bound. Tested here rather than through
    /// `normalize` because the end-to-end path routes through `best_alignment`'s
    /// max-match search, which makes the *fixture* rather than the bound decide
    /// the outcome.
    #[test]
    fn partition_block_bounds_a_long_net_deletion_but_not_its_siblings() {
        let reference: Vec<u8> = b"ACGT".repeat(13); // 52 nt, past the 32 nt bound
        assert_eq!(reference.len(), 52);

        // Net deletion: unguarded, so the bound must keep it whole.
        let deletion_result = b"TTCCTCGATGCCTG".to_vec(); // 14 nt
        assert!(deletion_result.len() < reference.len());
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
        assert!(
            reference.len() > MAX_UNGUARDED_SPLIT_BLOCK,
            "fixture must exceed the unguarded bound or it proves nothing"
        );
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
    /// codon exception (its `partition_block` sibling never reaches
    /// `best_alignment` at all — see `MAX_UNGUARDED_SPLIT_BLOCK`), so
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
        for reference in &all {
            for result in &all {
                let Some(pieces) = partition_block_sequence_first(
                    reference,
                    result,
                    crate::normalize::seqfirst::MIN_SEPARATION,
                ) else {
                    continue;
                };
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
        /// differently. Which is preferable is not decided here — the spec's
        /// prioritisation rule (`general.md:56`) is what reaches these, and the
        /// two blocks below fall to it in opposite directions.
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
                    shift_pieces_three_prime(&mut shift_first, block);
                    coalesce_adjacent_pieces(&mut shift_first);
                    shrink_pieces_to_differences(&mut shift_first, block);

                    let mut shrink_first = pieces.clone();
                    shrink_pieces_to_differences(&mut shrink_first, block);
                    shift_pieces_three_prime(&mut shrink_first, block);
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
            // separated by the retained base. `general.md:56` prefers those to
            // one spanning `delins`.
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
            // `delins` members, so it buys nothing `general.md:56` ranks above
            // the spanning `delins` it came from. #1235's own comment asks for
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
