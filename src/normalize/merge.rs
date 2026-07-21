//! Merge consecutive sub-variants in a cis allele into a single edit per HGVS spec.
//!
//! See `docs/superpowers/specs/2026-04-30-merge-consecutive-allele-edits-design.md`.

use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{
    Accession, AllelePhase, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, MtVariant, RnaVariant,
    TxVariant,
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
    // positive body (issue #920) — see the `region != body` refusal below.
    let kind = match &variants[0] {
        HgvsVariant::Genome(_) => CisKind::Genome,
        HgvsVariant::Mt(_) => CisKind::Mt,
        HgvsVariant::Cds(_) => CisKind::Cds,
        HgvsVariant::Tx(_) => CisKind::Tx,
        HgvsVariant::Rna(_) => CisKind::Rna,
        _ => return variants,
    };
    // The positive body region every member must lie in for this axis.
    let body = body_region(kind);
    // Borrow the head's accession as the template once; the helper rejects any
    // variant whose kind doesn't match `kind`.
    let Some((template_accession, _, _, _, _)) = cis_axis_parts(&variants[0], kind) else {
        return variants;
    };
    let template_accession = template_accession.clone();

    let mut edits: Vec<GEdit> = Vec::with_capacity(variants.len());
    for v in &variants {
        let Some((accession, region, s, e, edit)) = cis_axis_parts(v, kind) else {
            return variants;
        };
        if *accession != template_accession {
            return variants;
        }
        // CDS-body-only scope (issue #920): every member must lie in the
        // positive body region for the axis — `Cds` for `c.`, `Tx` for `n.`,
        // `Rna` for `r.`, `Genome` for `g.`/`m.`. Any 5'UTR (`c.-N`), 3'UTR
        // (`c.*N`), intronic-offset, upstream/downstream, or mixed-region
        // member refuses the whole group (all-or-nothing, mirroring the other
        // conservative refusals). For `g.`/`m.` this is always satisfied.
        if region != body {
            return variants;
        }
        match edit {
            NaEdit::Insertion { sequence } => {
                let Some(bases) = sequence.bases() else {
                    return variants;
                };
                // Insertion location is the two flanking positions `[gap, gap+1]`.
                if e != s + 1 {
                    return variants;
                }
                edits.push(GEdit::Ins {
                    gap: s,
                    alt: bases.to_vec(),
                });
            }
            NaEdit::Deletion { .. } => edits.push(GEdit::Del { s, e }),
            NaEdit::Substitution { alternative, .. } => {
                if s != e {
                    return variants;
                }
                edits.push(GEdit::Sub {
                    pos: s,
                    alt: *alternative,
                });
            }
            NaEdit::Delins { sequence, .. } => {
                let Some(bases) = sequence.bases() else {
                    return variants;
                };
                edits.push(GEdit::Delins {
                    s,
                    e,
                    alt: bases.to_vec(),
                });
            }
            // A certain, simple tandem duplication is an insertion of its own
            // reference span at gap `e` (#999/#1000: a per-member 3'-shift can
            // canonicalise an insertion to a `dup` that lands flush against an
            // adjacent substitution/deletion). Refuse the uncertain-extent form
            // (`dup?`, `dup(731_741)`) — its span is not pinned down.
            NaEdit::Duplication {
                uncertain_extent: None,
                ..
            } => edits.push(GEdit::Dup { s, e }),
            // repeat / inversion / identity / uncertain dup etc. — not handled.
            _ => return variants,
        }
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
            GEdit::Del { s, e } | GEdit::Delins { s, e, .. } => Some((*s, *e)),
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
    let delta: i64 = match kind {
        CisKind::Genome | CisKind::Mt | CisKind::Tx => 0,
        CisKind::Cds | CisKind::Rna => {
            let Ok(tx) = provider.get_transcript(&accession) else {
                return variants;
            };
            let Some(cds_start) = tx.cds_start else {
                return variants;
            };
            let Ok(cds_start) = i64::try_from(cds_start) else {
                return variants;
            };
            cds_start - 1
        }
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
    let mut lo = 0usize;
    while lo < ref_bytes.len() && lo < variant.len() && ref_bytes[lo] == variant[lo] {
        lo += 1;
    }
    let mut hi_ref = ref_bytes.len();
    let mut hi_var = variant.len();
    while hi_ref > lo && hi_var > lo && ref_bytes[hi_ref - 1] == variant[hi_var - 1] {
        hi_ref -= 1;
        hi_var -= 1;
    }
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

/// Coalesce runs of consecutive-residue substitutions in a **cis protein
/// allele** into a single delins, per `protein/substitution.md:23`:
///
/// > changes involving two or more consecutive amino acids are described as a
/// > deletion/insertion variant (delins) […] the description
/// > `p.Arg76_Cys77delinsSerTrp` is correct, the description
/// > `p.[Arg76Ser;Cys77Trp]` is not correct.
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
/// Returns `Some` only when at least one run of ≥2 adjacent substitutions is
/// merged; a fully-merged allele collapses to a bare [`HgvsVariant::Protein`],
/// a partial merge to a smaller [`HgvsVariant::Allele`]. Returns `None` — leave
/// the allele untouched — when it is conservatively out of scope:
///   - the phase is not cis (trans / unknown / mosaic / … are independent),
///   - any member is not a single-residue protein substitution (mixed edit
///     kinds are out of scope for this narrow rule),
///   - members are not in ascending, duplicate-free residue order (a `n,n`
///     pair is the contradictory same-residue case, not a delins; an
///     out-of-order allele is left as the author wrote it), or
///   - a run would mix predicted `( )` and certain members (ambiguous).
pub(crate) fn coalesce_protein_adjacent_substitutions(
    allele: &crate::hgvs::variant::AlleleVariant,
) -> Option<HgvsVariant> {
    use crate::hgvs::edit::{AminoAcidSeq, ProteinEdit};
    use crate::hgvs::interval::ProtInterval;
    use crate::hgvs::location::{AminoAcid, ProtPos};
    use crate::hgvs::variant::ProteinVariant;

    if allele.phase != AllelePhase::Cis || allele.variants.len() < 2 {
        return None;
    }

    /// One member reduced to a single-residue substitution.
    struct Sub {
        pos: u64,
        reference: AminoAcid,
        alternative: AminoAcid,
        predicted: bool,
    }

    let first = match &allele.variants[0] {
        HgvsVariant::Protein(pv) => pv,
        _ => return None,
    };
    let accession = first.accession.clone();
    let gene_symbol = first.gene_symbol.clone();

    let mut subs: Vec<Sub> = Vec::with_capacity(allele.variants.len());
    for v in &allele.variants {
        let HgvsVariant::Protein(pv) = v else {
            return None;
        };
        // Exactly one certain residue (a point, start == end, both certain).
        let (Some(Mu::Certain(s)), Some(Mu::Certain(e))) = (
            pv.loc_edit.location.start.as_single(),
            pv.loc_edit.location.end.as_single(),
        ) else {
            return None;
        };
        if s != e {
            return None;
        }
        let (alternative, predicted) = match &pv.loc_edit.edit {
            Mu::Certain(ProteinEdit::Substitution { alternative, .. }) => (*alternative, false),
            Mu::Uncertain(ProteinEdit::Substitution { alternative, .. }) => (*alternative, true),
            _ => return None,
        };
        subs.push(Sub {
            pos: s.number,
            reference: s.aa,
            alternative,
            predicted,
        });
    }

    // Require ascending, duplicate-free residue order. A descending or
    // duplicate (`n,n`) allele is out of scope — do not reorder the author's
    // description, and the same-residue pair is a contradiction, not a delins.
    if subs.windows(2).any(|w| w[1].pos <= w[0].pos) {
        return None;
    }

    let mut merged: Vec<HgvsVariant> = Vec::new();
    let mut changed = false;
    let mut i = 0;
    while i < subs.len() {
        let mut j = i;
        while j + 1 < subs.len() && subs[j + 1].pos == subs[j].pos + 1 {
            j += 1;
        }
        if j > i {
            // A run of ≥2 strictly-adjacent substitutions → one delins spanning
            // residues `i..=j`, inserting each member's alternative in order.
            let all_predicted = subs[i..=j].iter().all(|s| s.predicted);
            let any_predicted = subs[i..=j].iter().any(|s| s.predicted);
            if any_predicted && !all_predicted {
                return None;
            }
            let loc = ProtInterval::new(
                ProtPos::new(subs[i].reference, subs[i].pos),
                ProtPos::new(subs[j].reference, subs[j].pos),
            );
            let edit = ProteinEdit::Delins {
                sequence: AminoAcidSeq::new(subs[i..=j].iter().map(|s| s.alternative).collect()),
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
            // A lone substitution keeps its original member verbatim.
            merged.push(allele.variants[i].clone());
        }
        i = j + 1;
    }

    if !changed {
        return None;
    }
    if merged.len() == 1 {
        merged.into_iter().next()
    } else {
        let mut coalesced = crate::hgvs::variant::AlleleVariant::new(merged, allele.phase);
        coalesced.uncertain = allele.uncertain;
        Some(HgvsVariant::Allele(coalesced))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
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
}
