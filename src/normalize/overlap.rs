//! Detect overlapping cis-allele sub-variants.
//!
//! See `docs/superpowers/specs/2026-05-05-allele-overlap-same-position-design.md`.
//!
//! Two cis-phase, same-accession sub-variants overlap when they describe edits
//! on the same reference territory — a case with no canonical HGVS form.
//! Rather than silently picking a winner, ferro preserves the input verbatim
//! and emits one [`NormalizationWarning::OverlapConflict`] per overlap (strict
//! mode promotes it to an error).
//!
//! Two detectors run at different points in the allele pipeline:
//!
//! - [`detect_overlap_conflicts`] — *post-shift* coincident bounds: span edits
//!   (`sub`/`del`/`delins`/`dup`/`inv`/`repeat`) whose ranked start/end keys are
//!   identical. Insertions are excluded here (they anchor at a boundary, not a
//!   single-base span).
//! - [`detect_insertion_overlaps`] — *pre-merge* junction overlaps: two
//!   junction-anchored edits at one junction *whose order nothing determines*,
//!   or one interior to a span edit
//!   (mutalyzer `EOVERLAP`, #486). Must run before the normalizer's merge step
//!   collapses overlapping cis edits into one combined edit. "Junction-
//!   anchored" covers a true `ins` plus the `dup` and `repeat` spellings (which
//!   land their extra copies at the junction after their own span), so the
//!   detector sees the same conflict whichever way the pipeline has spelled it.

use std::collections::BTreeMap;

use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{AllelePhase, HgvsVariant, LocEdit};
use crate::normalize::footprint::WriteFootprint;
use crate::normalize::merge::Region;
use crate::normalize::NormalizationWarning;
use smol_str::SmolStr;

/// One position on a transcript-relative axis, ordered the way the axis reads.
///
/// # Why a rank rather than a bare `(Region, i64)`
///
/// The `c.`/`r.`/`n.` axes are piecewise: three regions numbered by three
/// independent integer sequences, so a bare coordinate does not order against
/// one in another region — `c.*1` is the base immediately 3' of the last CDS
/// base, and `1 < 15`. Both detectors in this module previously handled that by
/// refusing to read a range whose two endpoints were in different regions
/// (`if rs != re { return None }` in `cds_range` / `tx_range` / `rna_range`) and
/// by comparing only members that shared a region.
///
/// That was not a decline. `None` reaches `let Some(..) = .. else { continue }`,
/// so the member was dropped from the analysis rather than treated
/// conservatively, and **strict mode silently accepted overlaps it rejects when
/// the same pair sits inside one region** (#1508). Measured on a 20-base
/// transcript with CDS `1..=15`:
///
/// ```text
/// c.[11_12dup;12_13inv]      REJECTED  W5002
/// c.[14_15dup;15_*1inv]      ACCEPTED            <-- same shape, across the boundary
/// c.[11_12del;11_12inv]      REJECTED  W5002
/// c.[15_*1del;15_*1inv]      ACCEPTED
/// c.[11_12insAA;11_12insCC]  REJECTED  W5002
/// c.[15_*1insAA;15_*1insCC]  ACCEPTED  -> c.[15delinsCAA;15_*1insCC]
/// ```
///
/// Ranking the regions gives the whole axis one order, so the boundary stops
/// being a place where members become invisible to each other.
///
/// # Why this is sound here and was rejected for `merge.rs`
///
/// #1482 considered exactly this key for `merge::MemberSpan` and rejected it:
/// 42 sites there do *arithmetic* on a span, many computing `end - start + 1`
/// as a length, and for `c.15_*1` that is `1 - 15`. A rank-ordered key would
/// have made all 42 compile while silently producing wrong lengths, so that
/// module converts onto the sequence axis with a provider instead.
///
/// This module has no provider and needs none, because it never takes a length
/// or an offset — every use is a **comparison**: do two spans coincide, do they
/// intersect, does a junction fall interior to a span. The one place that looked
/// like arithmetic, `gap + 1 <= end`, is already written `gap < end`.
///
/// The ordering itself is the one [`crate::hgvs::variant`]'s `CanonicalPos`
/// already gives the parser's E3006 self-cancelling check — which is why that
/// check gets a cross-region `c.[14_15dup;15_*1del]` right while these
/// detectors do not. This adopts the precedent rather than inventing an order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct AxisPos {
    /// `0` = 5' of the body, `1` = the body, `2` = 3' of it. First in
    /// declaration order so the derived `Ord` is lexicographic.
    rank: u8,
    /// The written coordinate, whose numeric order is already correct *within*
    /// a region — including the 5' regions, where `-3 < -2 < -1` runs toward
    /// the body exactly as the axis does.
    coord: i64,
}

impl AxisPos {
    fn new(region: Region, coord: i64) -> Self {
        let rank = match region {
            Region::FivePrimeUtr | Region::TxUpstream => 0,
            Region::Genome | Region::Cds | Region::Rna | Region::Tx => 1,
            Region::ThreePrimeUtr | Region::TxDownstream => 2,
        };
        Self { rank, coord }
    }

    /// Whether `next` is the position immediately 3' of `self`.
    ///
    /// Needed because an insertion's anchor is recognised by its two endpoints
    /// being adjacent, and at a region boundary they are adjacent without being
    /// consecutive integers: `c.15_*1insCC` anchors at the junction between the
    /// last CDS base and the first 3'UTR base.
    ///
    /// Crossing *into* the body is fully checkable — the 5' regions end at `-1`.
    /// Crossing *out* of it is not: whether `c.15` is the last CDS base depends
    /// on the record, and this module holds no provider. HGVS insertion syntax
    /// requires the two flanking positions to be adjacent, so the author's
    /// spelling is taken at its word there — the same latitude the parser
    /// already extends, since it accepts `c.15_*1insCC` without consulting a
    /// reference either.
    fn is_immediately_followed_by(self, next: Self) -> bool {
        if self.rank == next.rank {
            return next.coord == self.coord + 1;
        }
        next.rank == self.rank + 1 && next.coord == 1 && (self.rank != 0 || self.coord == -1)
    }
}

/// Detect coincident-bounds groups in a cis allele.
///
/// Returns one warning per group of `>= 2` same-accession sub-variants whose
/// ranked start/end positions ([`AxisPos`]) are identical.
/// Insertions are excluded — they anchor at boundaries (`[end, start]`) and
/// have no single-base location to coincide on. Trans / mosaic / chimeric /
/// unknown phases short-circuit (the warning only applies to a same-haplotype
/// allele, where the conflict is real).
///
/// Input `variants` are the post-normalization sub-variants of an
/// [`AlleleVariant`]. The pass is purely observational: it does not mutate
/// or reorder its input. Warnings are emitted in deterministic key order
/// (BTreeMap iteration), so two equivalent inputs yield identical warning
/// sequences regardless of source-line ordering of the conflicting edits.
pub(crate) fn detect_overlap_conflicts(
    variants: &[HgvsVariant],
    phase: AllelePhase,
) -> Vec<NormalizationWarning> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return Vec::new();
    }

    // Group sub-variants by (accession, coord_system, start, end), the endpoints
    // ranked so a member spanning a region boundary groups with the rest rather
    // than vanishing (#1508). The `Vec<usize>` collects 0-based indices into
    // `variants` in input order so each group's `edit_kinds` reflect the source
    // ordering of the conflict.
    let mut groups: BTreeMap<GroupKey, Vec<usize>> = BTreeMap::new();
    for (idx, variant) in variants.iter().enumerate() {
        let Some(key) = group_key(variant) else {
            continue;
        };
        groups.entry(key).or_default().push(idx);
    }

    let mut warnings = Vec::new();
    for (key, indices) in &groups {
        if indices.len() < 2 {
            continue;
        }
        let edit_kinds: Vec<String> = indices
            .iter()
            .filter_map(|&i| edit_kind(&variants[i]).map(|s| s.to_string()))
            .collect();
        // group_key already filtered to variants with a known edit kind, so
        // edit_kinds.len() == indices.len() in practice. Preserved as a
        // debug assert against future drift.
        debug_assert_eq!(edit_kinds.len(), indices.len());
        let location_str = location_for_variant(&variants[indices[0]])
            .expect("group_key established a renderable location");

        warnings.push(NormalizationWarning::OverlapConflict {
            accession: key.accession.to_string(),
            coordinate_system: key.coord_system.to_string(),
            location: location_str,
            edit_kinds,
        });
    }

    // Footprints that *collide* without coinciding (#1448).
    //
    // The grouping above keys on exact `(accession, coord_system, start, end)`
    // equality, so two members that overlap only partially — or that nest —
    // land in different groups and were reported by nothing. Strict mode
    // therefore accepted them, and they denote no sequence at all: the
    // independent applier declines them in either member order.
    //
    //     g.[10_14del;12_16del]        ACCEPTED, applies to nothing
    //     g.[10_14inv;12_16inv]        ACCEPTED, applies to nothing
    //     g.[10_14del;10_16del]        ACCEPTED, applies to nothing (nested)
    //
    // Measured over two-member span-edit alleles on four sequences: 10,499 of
    // the 25,848 strict mode accepted — 41% — denoted nothing.
    //
    // **The geometry is [`WriteFootprint`]'s, not this loop's** (#1749). This
    // pass used to carry its own `span_writer_footprint`, hard-restricted to
    // `sub`/`del`/`delins`/`inv`, which silently dropped `repeat` — so a repeat
    // whose tract intersected a sibling deletion was reported by nothing at
    // all. That was one of three answers this file gave for `repeat`; there is
    // now one, and it is stated in exactly one place.
    //
    // Pairwise rather than grouped, and skipping pairs whose keys are equal, so
    // an exactly-coincident pair is reported once — by the loop above, which
    // names every member of its group — rather than twice.
    for (i, first) in variants.iter().enumerate() {
        let Some((coord_system, first_accession, first_footprint)) = member_footprint(first) else {
            continue;
        };
        // Junction geometry belongs to [`detect_insertion_overlaps`], which runs
        // *pre-merge*, where a junction-anchored sibling is still observable as
        // its own member. By the time this pass runs the merge step has
        // collapsed them, so a junction test here would report a collision
        // whose other half has left the list.
        let Some(first_bases) = first_footprint.bases else {
            continue;
        };
        for (j, second) in variants.iter().enumerate().skip(i + 1) {
            let Some((second_system, second_accession, second_footprint)) =
                member_footprint(second)
            else {
                continue;
            };
            let Some(second_bases) = second_footprint.bases else {
                continue;
            };
            // No region-equality gate (#1508). It was never a molecule test —
            // `coord_system` is that — and gating on it hid every pair one of
            // whose members crosses a region boundary. Ranked positions compare
            // across the boundary, so the geometry test below decides.
            if first_accession != second_accession || coord_system != second_system {
                continue;
            }
            // Coincident pairs belong to the grouped loop above.
            if first_bases == second_bases {
                continue;
            }
            // Base ranges only, for the reason given above — so the comparison
            // is made on span-only views of the two footprints. A reversed
            // range is declined rather than judged, inside `conflicts_with`.
            let bases_only = |bases: (AxisPos, AxisPos), preserves: bool| {
                WriteFootprint::spanning(bases.0, bases.1, preserves)
            };
            if !bases_only(first_bases, first_footprint.preserves_bases)
                .conflicts_with(&bases_only(second_bases, second_footprint.preserves_bases))
            {
                continue;
            }
            let kinds: Vec<String> = [i, j]
                .iter()
                .filter_map(|&idx| member_kind(&variants[idx]).map(|s| s.to_string()))
                .collect();
            warnings.push(NormalizationWarning::OverlapConflict {
                accession: first_accession.to_string(),
                coordinate_system: coord_system.to_string(),
                location: location_for_variant(first)
                    .expect("member_footprint established a renderable location"),
                edit_kinds: kinds,
            });
        }
    }
    warnings
}

/// Where a member writes, on this module's ranked axis — the single definition
/// every pass in this file shares (#1749).
///
/// Returns the accession and coordinate system alongside it, because a
/// footprint means nothing except against another footprint on the same
/// molecule and the same axis.
///
/// The per-edit mapping, and the argument for each entry, is documented on
/// [`crate::normalize::footprint`]. `None` for anything with no definite
/// footprint: a non-`NaEdit` variant, an uncertain edit, a position this module
/// declines to read (intronic offset, `?` sentinel, `pter`/`qter`), or an edit
/// kind with no reference geometry (`conversion`, `methylation`, `copy number`,
/// an N-padded deletion over an uncertain range, …).
fn member_footprint(
    variant: &HgvsVariant,
) -> Option<(&'static str, SmolStr, WriteFootprint<AxisPos>)> {
    let (coord_system, start, end) = simple_span(variant)?;
    // Rendered into a `SmolStr`, which is inline for every real accession, and
    // identical bytes: `full_smol` drives the same writer `Display` does.
    //
    // This is the one place the accession is rendered for this whole module —
    // #1749 collapsed four call sites into this one — and the *inner* loop of
    // `detect_overlap_conflicts` calls it, so a `String` here is an allocation
    // per **pair** of members. It is only ever compared; the owned form is built
    // where a warning is actually emitted.
    let accession = variant.accession()?.full_smol();
    let edit = inner_edit(variant)?;
    let footprint = match edit {
        // Rewrite the span they name, and leave bases standing there.
        NaEdit::Substitution { .. }
        | NaEdit::SubstitutionNoRef { .. }
        | NaEdit::Delins { .. }
        | NaEdit::Inversion { .. } => WriteFootprint::spanning(start, end, true),
        // Rewrites its span and leaves nothing, which is what makes an interior
        // junction meaningless against it (#1406).
        NaEdit::Deletion { .. } => WriteFootprint::spanning(start, end, false),
        // Reads its span and writes the copy directly 3' of it
        // (`duplication.md:5`) — but only when the extent is certain. An
        // uncertain-extent `dup` is not known to write at that junction, so it
        // is not entitled to the narrow footprint and keeps its span.
        NaEdit::Duplication {
            uncertain_extent: None,
            ..
        } => WriteFootprint::at_junction(end),
        NaEdit::Duplication { .. } => WriteFootprint::spanning(start, end, true),
        // A canonical insertion anchors at two adjacent positions; anything
        // else (e.g. a malformed single-position insertion) has no junction.
        //
        // Adjacency, not `end == start + 1` (#1508): at a region boundary the
        // two flanking positions are adjacent without being consecutive
        // integers, so the integer test silently declined to register
        // `c.15_*1insCC` as occupying a junction at all — and two of those at
        // one junction were accepted where the in-region pair is rejected.
        NaEdit::Insertion { .. } => {
            if !start.is_immediately_followed_by(end) {
                return None;
            }
            WriteFootprint::at_junction(start)
        }
        // A repeat states both the tract it replaces and what it becomes, so
        // which of the two it does is derivable from the description alone.
        NaEdit::Repeat {
            sequence,
            count,
            additional_counts,
            ..
        } => repeat_footprint(sequence.as_ref(), count, additional_counts, start, end)?,
        _ => return None,
    };
    Some((coord_system, accession, footprint))
}

/// Where a `repeat` writes, from the description alone.
///
/// A repeat states the tract it replaces (its span) *and* what that tract
/// becomes (`unit x count`), so whether it grows or shrinks is arithmetic on
/// the description — no provider needed, which is what makes this better than
/// either answer the three passes used to give:
///
/// - **grows, or stays the same length** (`unit x count >= tract`): the tract's
///   own bases are asserted unchanged and the extra copies land at the junction
///   3' of it. That is exactly a `dup`'s geometry, which is what
///   `issue_1437_dup_read_span_interior_sibling` requires: the per-member
///   pipeline picks between the `dup` and `repeat` spellings **by axis** (the
///   CDS axis rewrites `c.5_6insAA` to `c.5_6dup`, the genomic axis rewrites the
///   same shape to `g.1005_1006A[4]`), so a detector that answered differently
///   for the two would be sensitive to a choice the author never made.
/// - **shrinks** (`unit x count < tract`): the trailing `tract - unit x count`
///   bases are removed and the leading ones are untouched reference. The
///   footprint is that trailing range, and its bases do **not** survive — so an
///   interior junction is meaningless against it for the #1406 reason, while a
///   junction among the *kept* prefix bases has a well-defined position and is
///   likewise no conflict.
///
/// This supersedes the note that used to sit on `writes_only_at_a_junction`,
/// which held that `A[8]`, `G[6]` and `TA[3]` are "indistinguishable at this
/// layer" and so a repeat must conservatively keep its whole span. They are
/// indistinguishable in *content*, which is what a `REFSEQ_MISMATCH` check is
/// for — `g.4_9G[6]` against a `TTTTTT` reference is a mismatched description,
/// not an overlap — but they are perfectly distinguishable in *length*, and
/// length is all this geometry needs.
///
/// `None` when the arithmetic is not available: no stated unit
/// (`g.262_263[6]`), an inexact count (`[10_15]`, `[(10_15)]`, `[10_?]`),
/// genotype notation carrying more than one count, or endpoints in different
/// regions, where a tract length is not a subtraction. Declining leaves the
/// member out of the analysis entirely, which is the same conservative choice
/// every other undecidable shape in this module makes.
fn repeat_footprint(
    sequence: Option<&crate::hgvs::edit::Sequence>,
    count: &crate::hgvs::edit::RepeatCount,
    additional_counts: &[crate::hgvs::edit::RepeatCount],
    start: AxisPos,
    end: AxisPos,
) -> Option<WriteFootprint<AxisPos>> {
    use crate::hgvs::edit::RepeatCount;

    if !additional_counts.is_empty() {
        return None;
    }
    let RepeatCount::Exact(copies) = count else {
        return None;
    };
    let unit_len = u64::try_from(sequence?.0.len()).ok()?;
    // A tract length is a subtraction, so both endpoints must be numbered by
    // the same sequence — the one thing a ranked `AxisPos` deliberately does
    // *not* let you do across a region boundary (#1482).
    if start.rank != end.rank || start.coord > end.coord {
        return None;
    }
    let tract = end.coord.checked_sub(start.coord)?.checked_add(1)?;
    let resulting = i64::try_from(unit_len.checked_mul(*copies)?).ok()?;

    if resulting >= tract {
        // Grows (or is length-neutral): the write is the junction 3' of the
        // tract, exactly as for a `dup`.
        return Some(WriteFootprint::at_junction(end));
    }
    // Shrinks: the trailing `tract - resulting` bases go, the leading
    // `resulting` are untouched reference.
    let removed_from = AxisPos {
        rank: start.rank,
        coord: start.coord.checked_add(resulting)?,
    };
    Some(WriteFootprint::spanning(removed_from, end, false))
}

/// The label the emitted warning reports for a member (`"ins"` / `"dup"` /
/// `"repeat"` / `"del"` / …).
///
/// [`edit_kind`] predates insertions participating in these detectors and still
/// declines them, so the insertion label is supplied here rather than widening
/// a function several other call sites use as an "is this a span edit" filter.
fn member_kind(variant: &HgvsVariant) -> Option<&'static str> {
    if matches!(inner_edit(variant), Some(NaEdit::Insertion { .. })) {
        return Some("ins");
    }
    edit_kind(variant)
}

/// Detect overlaps that involve at least one junction-anchored edit within a
/// cis allele.
///
/// An insertion `a_(a+1)ins…` occupies the zero-width junction between
/// reference positions `a` and `a+1`. It overlaps:
///
/// - **another insertion** at the *same* junction (`[4_5insT;4_5insA]`) — two
///   payloads into one slot with nothing to order them, for which
///   `general.md:78` supplies the designated spelling `ins[A;B]`; and
/// - **a span edit** (`del`/`delins`/`dup`/`inv`/`sub`/`repeat`) whose
///   reference range *strictly encloses* the junction, i.e.
///   `range.start <= a` and `a + 1 <= range.end`
///   (`[274_275delinsT;274_275insA]`).
///
/// It does **not** overlap a `dup`/`repeat` that writes into the same junction,
/// when that is the only other occupant. `DNA/duplication.md:90` publishes such
/// a pair as a correct description and glosses its order — the duplication
/// "**followed by** the insertion" — so the composition is determined rather
/// than ambiguous. Two `dup`s at one junction (`[3_6dup;5_6dup]`) carry no such
/// gloss and remain a conflict.
///
/// An insertion abutting a span's edge — e.g. `100_101ins` next to a single-
/// base `100` substitution — is *not* interior, so it does not overlap. This
/// keeps the spec-valid `[273_274insT;274G>T;274_275insA]` accepted.
///
/// A **`dup` and a `repeat` also occupy a junction** (see
/// [`junction_writing_kind`]): each writes its extra copies directly 3' of its
/// own span, i.e. at the junction after `end`. Both are therefore registered as
/// junction occupants *in addition to* being spans, so `[5_9inv;5_6dup]` and
/// `[1005_1009inv;1005_1006A[4]]` are each flagged just as the `ins`-spelled
/// input they normalize from is.
///
/// Without this the detector is spelling-sensitive and normalization is not
/// idempotent on a conflicting allele. The first pass respells the interior
/// `insAA` — as a `dup` on the CDS axis, as a `repeat` on the genomic axis —
/// the second pass no longer recognises either as a conflict, and so it
/// reorders members that the first pass deliberately left in authored order
/// (`#395`). Covering only one of the two respellings fixes only one axis.
///
/// An occupant's own span never encloses its own junction (`gap == end`, and
/// the interior test is `gap < end`), so a lone `dup`/`repeat` cannot conflict
/// with itself.
///
/// One warning is emitted per same-junction occupant group and one per span
/// edit that encloses ≥1 occupant. Iteration is in deterministic order
/// (BTreeMap junction key, then input index) so equivalent inputs yield
/// identical warning sequences.
///
/// This must run on the *pre-merge* allele members: the normalizer collapses
/// overlapping cis edits into a single combined edit before
/// [`detect_overlap_conflicts`] sees the post-shift list, so by then the
/// overlap is no longer observable as two sub-variants.
pub(crate) fn detect_insertion_overlaps(
    variants: &[HgvsVariant],
    phase: AllelePhase,
) -> Vec<NormalizationWarning> {
    junction_overlaps(variants, phase, true)
}

/// Only the *interior-junction* half of [`detect_insertion_overlaps`]: a
/// junction strictly inside a span edit that keeps the bases it spans.
///
/// This is the subset that is genuinely ambiguous about what the allele
/// denotes, and it is the only subset the #1406 conflict-preservation gates may
/// act on. The same-junction half must stay out of them: two insertions at one
/// junction are a two-member **spelling** of a variant that also has a
/// single-member spelling, and merging them is what makes the two converge.
/// Treating that as unresolvable made `g.[263_264insAC;264_265insAA]` settle
/// apart from `g.264_265insCAAA` — measured, and caught by
/// `cis_spelling_confluence_gap::converged_pairs_stay_converged`, which reports
/// it as "#1301 regressed". Losing confluence to gain a strict-mode property is
/// the wrong trade: confluence is what downstream consumers key on.
pub(crate) fn detect_interior_junction_conflicts(
    variants: &[HgvsVariant],
    phase: AllelePhase,
) -> Vec<NormalizationWarning> {
    junction_overlaps(variants, phase, false)
}

/// Shared body. `include_same_junction` selects branch (a) below.
fn junction_overlaps(
    variants: &[HgvsVariant],
    phase: AllelePhase,
    include_same_junction: bool,
) -> Vec<NormalizationWarning> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return Vec::new();
    }
    /// A member writing at a zero-width junction: a true `ins`, the copy a
    /// certain-extent `dup` lands 3' of its own span, or the expansion a
    /// `repeat` lands there. `kind` is the label the emitted warning reports.
    struct Insertion {
        idx: usize,
        accession: SmolStr,
        coord_system: &'static str,
        gap: AxisPos,
        kind: &'static str,
    }
    struct Span {
        idx: usize,
        accession: SmolStr,
        coord_system: &'static str,
        start: AxisPos,
        end: AxisPos,
        /// Whether the bases this member spans are still there afterwards.
        /// See branch (b).
        preserves_bases: bool,
    }

    // One pass, one definition (#1749). A member is registered wherever its
    // footprint says it writes, and the two registrations are independent: the
    // footprint carries `bases` and `junction` as separate `Option`s, and each
    // is consulted on its own rather than through a ladder that picks one.
    //
    // **No constructor produces both today**, and that is worth stating so the
    // `ins.idx != span.idx` self-exclusion in branch (b) is not read as
    // protection against a shape that exists. `repeat_footprint` answers
    // `at_junction(end)` when the repeat grows or is length-neutral and
    // `spanning(removed_from, end, false)` when it shrinks — a `repeat` writes
    // at a junction **or** over a span, never both — and `WriteFootprint` has
    // only those two constructors, with no call site building one by struct
    // literal. The independence is the model's, kept because it is the honest
    // geometry and because a future edit type could need it; it is not a live
    // capability, so the self-exclusion is defence-in-depth.
    //
    // That is the correction. This loop used to decide with its own
    // `junction_writing_kind` ladder and `continue` after a junction
    // registration, which made a `repeat` junction-*only* — so an insertion
    // strictly interior to a repeat tract was reported by nothing, while the
    // very same shape against a `delins` was rejected. The `continue` was right
    // for `dup`, whose write footprint genuinely *is* only that junction
    // (#1437), and was silently inherited by `repeat`, whose is not.
    let mut insertions: Vec<Insertion> = Vec::new();
    let mut spans: Vec<Span> = Vec::new();
    for (idx, variant) in variants.iter().enumerate() {
        let Some((coord_system, accession, footprint)) = member_footprint(variant) else {
            continue;
        };
        let Some(kind) = member_kind(variant) else {
            continue;
        };
        if let Some(gap) = footprint.junction {
            insertions.push(Insertion {
                idx,
                accession: accession.clone(),
                coord_system,
                gap,
                kind,
            });
        }
        if let Some((start, end)) = footprint.bases {
            spans.push(Span {
                idx,
                accession,
                coord_system,
                start,
                end,
                preserves_bases: footprint.preserves_bases,
            });
        }
    }

    let mut warnings = Vec::new();

    // (a) Two or more junction-anchored edits sharing a junction. Dups and
    // repeats count alongside true insertions: `[3_6dup;5_6dup]` writes both
    // copies into the slot after position 6 with no defined order between them,
    // exactly as `[5_6insA;5_6insT]` does. Branch (b) does not cover that pair —
    // the inner dup's junction sits at the outer dup's *edge*, not interior to
    // it — so omitting them here would leave the shape unreported even though
    // the sequence-first guard in `merge::…` refuses it (one `after[slot]`, two
    // writers). The same holds for the `repeat` spelling.
    // Each group carries `(index, kind)` so the warning can name every member by
    // its own spelling rather than labelling the whole group `ins`.
    /// `(accession, coordinate system, gap)` — one zero-width junction on one
    /// molecule. The ranked gap carries its own region (#1508).
    type JunctionKey = (SmolStr, &'static str, AxisPos);
    /// The members writing at a junction, as `(index into `variants`, kind)`.
    type Occupants = Vec<(usize, &'static str)>;

    let mut by_junction: BTreeMap<JunctionKey, Occupants> = BTreeMap::new();
    for ins in &insertions {
        by_junction
            .entry((ins.accession.clone(), ins.coord_system, ins.gap))
            .or_default()
            .push((ins.idx, ins.kind));
    }
    for ((accession, coord_system, _gap), occupants) in &by_junction {
        if !include_same_junction || occupants.len() < 2 {
            continue;
        }
        // …but a `dup`/`repeat` sharing a junction with ONE true `ins` is
        // ordered by the spec, not undetermined.
        //
        // `DNA/duplication.md:90` publishes exactly this pair as a correct
        // description — `NC_000001.11(NM_206933.2):c.[675-542_1211-703dup;
        // 1211-703_1211-702insGTAAA]`, whose insertion sits at the junction the
        // duplication writes its copy into — and glosses it "a duplication of
        // the sequence from … **followed by** the insertion of the sequence
        // `GTAAA`". That gloss supplies the very thing this branch's refusal
        // says is missing, and the only thing the attached NOTE rejects is the
        // `dupins` spelling.
        //
        // Keeping it split is also the reading that satisfies the *other*
        // clause in play. `duplication.md:18` requires a variant describable as
        // a duplication to be described as one, which is why the ledger's
        // `delins-adjacent-members-when-both-consume-reference` left its `dup`
        // half open rather than merging such a pair into one `delins` — merging
        // destroys the `dup`. Accepting the split pair is the only answer that
        // neither destroys the duplication nor refuses a shape the spec
        // publishes.
        //
        // **Narrow deliberately.** Two true insertions at one junction stay a
        // conflict: their order really is undetermined, and `general.md:78`
        // supplies a different designated spelling for that content
        // (`ins[A;B]`). Two junction-writing span edits stay a conflict too —
        // `[3_6dup;5_6dup]` writes both copies into one slot with nothing to
        // order them, and no clause glosses that pair. Only the mixed
        // one-of-each case is ordered.
        //
        // Branch (b) is untouched: `dup`/`repeat` remain registered as junction
        // occupants there, which is what #1446 is actually protecting — an
        // interior `ins` that the per-member pipeline respells as a `dup` (CDS
        // axis) or a `repeat` (genomic axis) must stay visible to the interior
        // test, or the detector becomes spelling-sensitive and normalization
        // stops being idempotent on a conflicting allele (#395).
        let insertions_here = occupants.iter().filter(|(_, k)| *k == "ins").count();
        let junction_writers_here = occupants.len() - insertions_here;
        if insertions_here <= 1 && junction_writers_here <= 1 {
            continue;
        }
        // Render the junction via the occupant's HGVS Display (like branch (b)
        // and `detect_overlap_conflicts`) so region prefixes (`*`/`-`) survive;
        // the raw signed `gap` drops them (e.g. 3'UTR `*1_*2` → `1_2`).
        let location = location_for_variant(&variants[occupants[0].0])
            .expect("same-junction occupant has a renderable location");
        warnings.push(NormalizationWarning::OverlapConflict {
            accession: accession.to_string(),
            coordinate_system: coord_system.to_string(),
            location,
            edit_kinds: occupants.iter().map(|(_, k)| k.to_string()).collect(),
        });
    }

    // (b) An insertion junction strictly interior to a span edit's range.
    //
    // **Except a pure deletion (#1406).** The conflict this branch reports is
    // that the junction's position within the span is meaningful and the
    // combination therefore has no single answer. A deletion removes every base
    // it spans, so an interior junction has nothing left to be positioned
    // against: `g.[2_3del;2_3insAA]` denotes `AA` in place of `2_3` whichever
    // order the two members are applied in, and whichever interior junction the
    // insertion had named. It composes uniquely, so it was never a conflict —
    // the same argument #1411 used to stop rejecting `g.[24dup;24C>G]`, whose
    // members likewise write to disjoint places.
    //
    // Reporting it anyway had a concrete cost beyond the false verdict: strict
    // mode rejected the input while accepting its own lenient output
    // (`g.2_3delinsAA`), which is the laundering #1406 row 3 is about, and the
    // merged form is what `delins.md:86-89` asks for outright.
    //
    // Deletion only, deliberately. A `delins` does **not** qualify: its payload
    // survives, so an interior insertion has a position relative to it —
    // `g.[2_3delinsGG;2_3insA]` genuinely does not say whether the `A` lands
    // before, inside or after the `GG`. Nor does `inv`, `dup` or `repeat`,
    // each of which keeps the spanned bases and so keeps the interior junction
    // meaningful. Those remain conflicts.
    for span in &spans {
        if !span.preserves_bases {
            continue;
        }
        let interior = insertions.iter().filter(|ins| {
            // Defence-in-depth, not protection against a shape the model can
            // build: no member is registered in both lists today, because
            // every `WriteFootprint` comes from `at_junction` or `spanning`
            // and neither sets both fields (see the registration loop above).
            // If one ever were, the interior test would already exclude it —
            // a `repeat`'s junction is its own `end` and `gap < end` is strict
            // — so this conjunct costs nothing and states the pairing rule
            // where it is relied on.
            ins.idx != span.idx
                && ins.accession == span.accession
                && ins.coord_system == span.coord_system
                && span.start <= ins.gap
                // `gap + 1 <= end` (junction interior); `gap < end` for ints.
                && ins.gap < span.end
        });
        let mut edit_kinds = vec![member_kind(&variants[span.idx])
            .expect("span edit has a known kind")
            .to_string()];
        edit_kinds.extend(interior.map(|ins| ins.kind.to_string()));
        if edit_kinds.len() < 2 {
            continue;
        }
        warnings.push(NormalizationWarning::OverlapConflict {
            accession: span.accession.to_string(),
            coordinate_system: span.coord_system.to_string(),
            location: location_for_variant(&variants[span.idx])
                .expect("span edit has a renderable location"),
            edit_kinds,
        });
    }

    warnings
}

/// Group key for the coincident-bounds detector.
///
/// `accession` is rendered as the canonical `Accession::Display` string so
/// equality is value-based (two variants pointing at distinct `Arc<str>`
/// instances of the same accession compare equal). `coord_system` is the
/// HGVS prefix character (`g`/`c`/`n`/`r`/`m`) — note `m.` differs from
/// `g.` even though both share `Region::Genome`, so we key on it explicitly
/// rather than collapsing.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct GroupKey {
    /// A [`SmolStr`], not a `String`: this is a map key and a comparand, never
    /// mutated, and every real accession fits the inline capacity — so grouping
    /// an allele's members no longer allocates one string per member. `SmolStr`
    /// orders as the `str` does, so the `BTreeMap` iteration order this pass's
    /// determinism rests on is unchanged.
    accession: SmolStr,
    coord_system: &'static str,
    /// Ranked endpoints, which subsume the old separate `region` field: two
    /// members coincide only if both ends agree, and a ranked position already
    /// carries which region it is in (#1508).
    start: AxisPos,
    end: AxisPos,
}

/// Build the group key for a sub-variant, or `None` if it can't participate
/// in coincident-bounds detection.
///
/// Returns `None` for:
/// - non-NaEdit variants (`Protein`, `RnaFusion`, `Allele`, `NullAllele`,
///   `UnknownAllele`, `Circular`)
/// - uncertain edits (`Mu::Uncertain` / `Mu::Unknown`)
/// - `NaEdit` kinds with no fixed reference span (`Insertion`, `DupIns`,
///   `Identity`, `Conversion`, `Unknown`, `Methylation`, `CopyNumber`,
///   `MultiRepeat`)
/// - positions with intronic offsets, `?` sentinels, or special anchors
///   (`pter`/`qter`/`cen`)
fn group_key(variant: &HgvsVariant) -> Option<GroupKey> {
    // The span is the shared definition's `bases` (#1749). A member whose write
    // is a junction rather than a span — a certain-extent `dup` — therefore has
    // no coincident-bounds key at all, and cannot be grouped with a sibling
    // that merely *reads* the same range. That is the #1411 reading, now taken
    // from one place instead of restated as a local `writes_only_at_a_junction`
    // ladder that had already drifted from its two neighbours.
    let (coord_system, accession, footprint) = member_footprint(variant)?;
    let (start, end) = footprint.bases?;
    edit_kind(variant)?;
    Some(GroupKey {
        accession,
        coord_system,
        start,
        end,
    })
}

/// The certain inner [`NaEdit`] of an NaEdit-bearing variant, or `None` for
/// non-NaEdit variants and uncertain edits.
fn inner_edit(variant: &HgvsVariant) -> Option<&NaEdit> {
    match variant {
        HgvsVariant::Genome(g) => g.loc_edit.edit.inner(),
        HgvsVariant::Cds(c) => c.loc_edit.edit.inner(),
        HgvsVariant::Tx(t) => t.loc_edit.edit.inner(),
        HgvsVariant::Rna(r) => r.loc_edit.edit.inner(),
        HgvsVariant::Mt(m) => m.loc_edit.edit.inner(),
        _ => None,
    }
}

/// Extract `(coord_system, ranked start, ranked end)` for a certain NaEdit-bearing
/// variant *regardless of edit kind*.
///
/// This is [`simple_range`] without the `is_overlap_edit` gate: insertion-
/// overlap detection needs the flanking-position span of an `Insertion`,
/// which `simple_range` deliberately drops. Position-validity gates (special
/// anchors, intronic offsets, `?` sentinels, region splits) still apply via
/// the per-coordinate-system range helpers.
fn simple_span(variant: &HgvsVariant) -> Option<(&'static str, AxisPos, AxisPos)> {
    fn na_span<L>(
        loc_edit: &LocEdit<Interval<L>, NaEdit>,
        range_fn: impl Fn(&Interval<L>) -> Option<(AxisPos, AxisPos)>,
    ) -> Option<(AxisPos, AxisPos)> {
        if !loc_edit.edit.is_certain() {
            return None;
        }
        range_fn(&loc_edit.location)
    }
    match variant {
        HgvsVariant::Genome(g) => na_span(&g.loc_edit, genome_range).map(|(s, e)| ("g", s, e)),
        HgvsVariant::Cds(c) => na_span(&c.loc_edit, cds_range).map(|(s, e)| ("c", s, e)),
        HgvsVariant::Tx(t) => na_span(&t.loc_edit, tx_range).map(|(s, e)| ("n", s, e)),
        HgvsVariant::Rna(r) => na_span(&r.loc_edit, rna_range).map(|(s, e)| ("r", s, e)),
        HgvsVariant::Mt(m) => na_span(&m.loc_edit, genome_range).map(|(s, e)| ("m", s, e)),
        _ => None,
    }
}

/// Short tag for the edit kind reported in the warning.
fn edit_kind(variant: &HgvsVariant) -> Option<&'static str> {
    let inner: Option<&NaEdit> = match variant {
        HgvsVariant::Genome(g) => g.loc_edit.edit.inner(),
        HgvsVariant::Cds(c) => c.loc_edit.edit.inner(),
        HgvsVariant::Tx(t) => t.loc_edit.edit.inner(),
        HgvsVariant::Rna(r) => r.loc_edit.edit.inner(),
        HgvsVariant::Mt(m) => m.loc_edit.edit.inner(),
        _ => None,
    };
    let edit = inner?;
    Some(match edit {
        NaEdit::Substitution { .. } | NaEdit::SubstitutionNoRef { .. } => "sub",
        NaEdit::Deletion { .. } => "del",
        NaEdit::Delins { .. } => "delins",
        NaEdit::Duplication { .. } => "dup",
        NaEdit::Inversion { .. } => "inv",
        NaEdit::Repeat { .. } => "repeat",
        // Excluded by `is_overlap_edit`; should never reach here for variants
        // that survived `group_key`. Returning None keeps `edit_kind` honest
        // when called from a future caller without the same prefilter.
        _ => return None,
    })
}

/// Render the canonical span text for the warning's `location` field.
///
/// Matches HGVS Display: a point uses bare `100` / `*1` / `-3`, a range
/// uses `start_end`. Returns `None` for non-NaEdit variants.
fn location_for_variant(variant: &HgvsVariant) -> Option<String> {
    match variant {
        HgvsVariant::Genome(g) => Some(format_interval(&g.loc_edit.location)),
        HgvsVariant::Cds(c) => Some(format_interval(&c.loc_edit.location)),
        HgvsVariant::Tx(t) => Some(format_interval(&t.loc_edit.location)),
        HgvsVariant::Rna(r) => Some(format_interval(&r.loc_edit.location)),
        HgvsVariant::Mt(m) => Some(format_interval(&m.loc_edit.location)),
        _ => None,
    }
}

/// Render an interval using the underlying position type's Display impl,
/// collapsing `start == end` to a single point.
fn format_interval<P: std::fmt::Display + PartialEq>(interval: &Interval<P>) -> String {
    let start = render_boundary(&interval.start);
    let end = render_boundary(&interval.end);
    if start == end {
        start
    } else {
        format!("{}_{}", start, end)
    }
}

/// Render a single side of an interval. `member_footprint` rejects any
/// variant whose boundary isn't `Single(Certain(_))`, so render_boundary
/// is only ever called on certain boundaries — the other arms exist to
/// document the invariant and panic loudly if it's ever violated.
fn render_boundary<P: std::fmt::Display>(boundary: &UncertainBoundary<P>) -> String {
    match boundary {
        UncertainBoundary::Single(Mu::Certain(p)) => p.to_string(),
        UncertainBoundary::Single(Mu::Uncertain(_))
        | UncertainBoundary::Single(Mu::Unknown)
        | UncertainBoundary::Range { .. } => {
            unreachable!("simple_range gates these out")
        }
    }
}

// --- Per-coordinate-system range extraction ----------------------------------
//
// These helpers mirror the equivalents in `merge.rs`. They are duplicated
// (rather than shared) because `merge.rs`'s versions are intentionally
// scoped to that module — re-exporting them all just to share four
// six-line helpers would widen merge.rs's surface area for no real
// benefit. The shared piece is `Region`, which is `pub(crate)`.

fn genome_range(interval: &Interval<GenomePos>) -> Option<(AxisPos, AxisPos)> {
    let s = simple_genome(interval.start.as_single()?)?;
    let e = simple_genome(interval.end.as_single()?)?;
    Some((
        AxisPos::new(Region::Genome, s),
        AxisPos::new(Region::Genome, e),
    ))
}

fn simple_genome(mu: &Mu<GenomePos>) -> Option<i64> {
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_special() || pos.offset.is_some() {
        return None;
    }
    i64::try_from(pos.base).ok()
}

fn cds_range(interval: &Interval<CdsPos>) -> Option<(AxisPos, AxisPos)> {
    // No cross-region refusal (#1508): each endpoint keeps its own region and
    // is ranked, so a span running from the body into a UTR is readable and
    // therefore visible to the detectors above.
    let (rs, s) = simple_cds(interval.start.as_single()?)?;
    let (re, e) = simple_cds(interval.end.as_single()?)?;
    Some((AxisPos::new(rs, s), AxisPos::new(re, e)))
}

fn simple_cds(mu: &Mu<CdsPos>) -> Option<(Region, i64)> {
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_unknown() || pos.is_intronic() {
        return None;
    }
    if pos.is_3utr() {
        return (pos.base >= 1).then_some((Region::ThreePrimeUtr, pos.base));
    }
    if pos.base < 0 {
        return Some((Region::FivePrimeUtr, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Cds, pos.base));
    }
    None
}

fn tx_range(interval: &Interval<TxPos>) -> Option<(AxisPos, AxisPos)> {
    // No cross-region refusal (#1508): each endpoint keeps its own region and
    // is ranked, so a span running from the body into a UTR is readable and
    // therefore visible to the detectors above.
    let (rs, s) = simple_tx(interval.start.as_single()?)?;
    let (re, e) = simple_tx(interval.end.as_single()?)?;
    Some((AxisPos::new(rs, s), AxisPos::new(re, e)))
}

fn simple_tx(mu: &Mu<TxPos>) -> Option<(Region, i64)> {
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_intronic() {
        return None;
    }
    if pos.is_downstream() {
        return (pos.base >= 1).then_some((Region::TxDownstream, pos.base));
    }
    if pos.base < 0 {
        return Some((Region::TxUpstream, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Tx, pos.base));
    }
    None
}

fn rna_range(interval: &Interval<RnaPos>) -> Option<(AxisPos, AxisPos)> {
    // No cross-region refusal (#1508): each endpoint keeps its own region and
    // is ranked, so a span running from the body into a UTR is readable and
    // therefore visible to the detectors above.
    let (rs, s) = simple_rna(interval.start.as_single()?)?;
    let (re, e) = simple_rna(interval.end.as_single()?)?;
    Some((AxisPos::new(rs, s), AxisPos::new(re, e)))
}

fn simple_rna(mu: &Mu<RnaPos>) -> Option<(Region, i64)> {
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_intronic() {
        return None;
    }
    if pos.is_3utr() {
        return (pos.base >= 1).then_some((Region::ThreePrimeUtr, pos.base));
    }
    if pos.base < 0 {
        return Some((Region::FivePrimeUtr, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Cds, pos.base));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;

    fn parse_allele(s: &str) -> (Vec<HgvsVariant>, AllelePhase) {
        let v = parse_hgvs(s).expect("parse");
        match v {
            HgvsVariant::Allele(a) => (a.variants, a.phase),
            other => panic!("expected allele in test, got {:?}", other),
        }
    }

    #[test]
    fn same_position_two_subs_emits_one_warning() {
        let (variants, phase) = parse_allele("NC_000001.11:g.[100G>A;100A>C]");
        let warnings = detect_overlap_conflicts(&variants, phase);
        assert_eq!(
            warnings.len(),
            1,
            "expected one warning, got {:?}",
            warnings
        );
        let NormalizationWarning::OverlapConflict {
            accession,
            coordinate_system,
            location,
            edit_kinds,
            ..
        } = &warnings[0]
        else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(accession, "NC_000001.11");
        assert_eq!(coordinate_system, "g");
        assert_eq!(location, "100");
        assert_eq!(edit_kinds, &vec!["sub".to_string(), "sub".to_string()]);
    }

    #[test]
    fn same_position_sub_plus_del_emits_one_warning() {
        let (variants, phase) = parse_allele("NC_000001.11:g.[100del;100A>C]");
        let warnings = detect_overlap_conflicts(&variants, phase);
        assert_eq!(warnings.len(), 1);
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!();
        };
        assert!(edit_kinds.contains(&"del".to_string()));
        assert!(edit_kinds.contains(&"sub".to_string()));
    }

    #[test]
    fn coincident_range_del_inv_emits_one_warning() {
        let (variants, phase) = parse_allele("NM_TEST.1:c.[100_103del;100_103inv]");
        let warnings = detect_overlap_conflicts(&variants, phase);
        assert_eq!(warnings.len(), 1);
        let NormalizationWarning::OverlapConflict {
            location,
            edit_kinds,
            ..
        } = &warnings[0]
        else {
            panic!();
        };
        assert_eq!(location, "100_103");
        assert!(edit_kinds.contains(&"del".to_string()));
        assert!(edit_kinds.contains(&"inv".to_string()));
    }

    #[test]
    fn three_subs_at_one_base_emit_one_group_warning() {
        let (variants, phase) = parse_allele("NC_000001.11:g.[100A>C;100A>G;100A>T]");
        let warnings = detect_overlap_conflicts(&variants, phase);
        assert_eq!(warnings.len(), 1, "groups, not pairs");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!();
        };
        assert_eq!(edit_kinds.len(), 3);
    }

    #[test]
    fn adjacent_subs_no_warning() {
        let (variants, phase) = parse_allele("NC_000001.11:g.[100A>C;101A>G]");
        assert!(detect_overlap_conflicts(&variants, phase).is_empty());
    }

    #[test]
    fn multi_accession_no_warning() {
        let (variants, phase) = parse_allele("[NC_000001.11:g.100A>C;NC_000002.11:g.100A>G]");
        assert!(detect_overlap_conflicts(&variants, phase).is_empty());
    }

    #[test]
    fn insertion_at_boundary_no_warning() {
        // Insertions anchor at p_p+1 — no single-base location to coincide.
        let (variants, phase) = parse_allele("NC_000001.11:g.[100A>C;100_101insT]");
        assert!(detect_overlap_conflicts(&variants, phase).is_empty());
    }

    // --- Insertion overlap detection (#486 EOVERLAP) -------------------------

    #[test]
    fn two_insertions_at_same_junction_emit_one_warning() {
        // Two insertions at the same interspace `4_5` overlap (mutalyzer
        // EOVERLAP). The inserted sequence is irrelevant — the junction is
        // shared.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_5insT;4_5insA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict {
            location,
            edit_kinds,
            ..
        } = &warnings[0]
        else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(location, "4_5");
        assert_eq!(edit_kinds, &vec!["ins".to_string(), "ins".to_string()]);
    }

    #[test]
    fn same_junction_insertions_render_utr_location_with_star_prefix() {
        // Regression: the same-junction branch must render the location via
        // HGVS Display, not the raw signed base. A 3'UTR junction `*1_*2`
        // would otherwise drop the `*` prefix and print `1_2`.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[*1_*2insT;*1_*2insA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { location, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(location, "*1_*2");
    }

    /// Two **overlapping `dup` spans** are not a coincident-bounds conflict (#1448).
    ///
    /// The intersection test added for #1448 is restricted to edits whose write
    /// footprint *is* their span. A `dup` reads its span and writes at the
    /// junction 3' of it (`duplication.md:5`), so overlapping read spans do not
    /// mean colliding writes — and keying the test on the read span instead,
    /// the reading #1411 corrected, would report these.
    ///
    /// Asserted here rather than end to end so it pins what *this* detector
    /// decides rather than the pipeline's net verdict.
    ///
    /// **The reason given here used to be wrong** (#1749). It read "because
    /// `detect_insertion_overlaps` reports an overlapping `dup` pair for its own
    /// separate reasons (#1437)". It does not, and #1437 is why: that change
    /// made `junction_overlaps` register a `dup` as a junction occupant and then
    /// stop, so it never becomes a span, and two `dup`s at different junctions
    /// match neither branch. Measured end to end on real hg38 in both strict and
    /// lenient, the only warning `g.[43045721dup;43045721_43045722dup]` raises is
    /// `MEMBERS_COALESCED_FROM_REPORTED_FORM`. The sibling test
    /// `overlapping_duplication_spans_are_silent_in_the_junction_detector_too`
    /// pins that, so the correction cannot drift back.
    #[test]
    fn overlapping_duplication_spans_are_not_a_coincident_conflict() {
        for input in [
            "NG_012337.1:g.[10_14dup;12_16dup]",
            "NG_012337.1:g.[10_14dup;10_16dup]",
        ] {
            let (variants, phase) = parse_allele(input);
            let warnings = detect_overlap_conflicts(&variants, phase);
            assert!(
                warnings.is_empty(),
                "`{input}`: two duplications write at different junctions, so \
                 their overlapping read spans are not a coincident-bounds \
                 conflict; got {warnings:?}"
            );
        }
    }

    #[test]
    fn insertion_interior_to_delins_emits_one_warning() {
        // An insertion `274_275ins` whose junction sits strictly inside a
        // `274_275delins` range overlaps it (mutalyzer EOVERLAP).
        let (variants, phase) = parse_allele("NM_TEST.1:c.[274_275delinsT;274_275insA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert!(edit_kinds.contains(&"delins".to_string()));
        assert!(edit_kinds.contains(&"ins".to_string()));
    }

    #[test]
    fn insertion_interior_to_deletion_is_not_a_conflict() {
        // Insertion junction `5_6` is strictly interior to the deleted range
        // `4_7`, but a deletion removes every base it spans, so the insertion
        // has nothing left to be positioned against: the pair denotes the same
        // bases whichever order it is applied in, and whichever interior
        // junction the insertion had named. Not a conflict (#1406).
        //
        // `insertion_interior_to_delins_emits_one_warning` is the
        // discriminating sibling: a `delins` keeps a payload, so an
        // interior insertion *does* have a position relative to it. Only a
        // pure deletion is exempt.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_7del;5_6insAA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert!(
            warnings.is_empty(),
            "an insertion interior to a pure deletion composes uniquely, so it \
             must not be reported; got {warnings:?}"
        );
    }

    #[test]
    fn insertions_sharing_junction_and_interior_to_span_emit_two_warnings() {
        // `5_6insA` and `5_6insT` share junction `5` (branch (a)) *and* both
        // sit strictly interior to the `4_7delinsGG` range positions 4..=7
        // (branch (b)), so the pass emits two `OverlapConflict` warnings for
        // the one overlap cluster. The combined-shape outcome is still a
        // correct rejection; this pins the per-branch warning count so a
        // future dedup refactor doesn't silently change it.
        //
        // The span is a `delins` rather than the `4_7del` this used to use:
        // a pure deletion no longer contributes a branch-(b) warning (#1406),
        // so with `del` here only branch (a) would fire and the test would no
        // longer exercise the two-branch shape it exists for.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_7delinsGG;5_6insA;5_6insT]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 2, "expected two warnings, got {warnings:?}");
        assert!(
            warnings
                .iter()
                .all(|w| matches!(w, NormalizationWarning::OverlapConflict { .. })),
            "all warnings must be OverlapConflict: {warnings:?}"
        );
    }

    #[test]
    fn duplication_interior_to_inversion_emits_one_warning() {
        // A `dup` lands its copy at the junction after its own last base, so
        // `5_6dup` occupies junction 6, strictly interior to `5_9inv`
        // (positions 5..=9). This is the `ins`-spelled conflict
        // `[4_10inv;5_6insAA]` after the per-member pipeline has respelled it,
        // and it must be recognised as the same conflict — otherwise
        // normalization is not idempotent on the conflicting allele.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[5_9inv;5_6dup]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict {
            location,
            edit_kinds,
            ..
        } = &warnings[0]
        else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(location, "5_9", "the warning names the enclosing span");
        assert_eq!(edit_kinds, &vec!["inv".to_string(), "dup".to_string()]);
    }

    #[test]
    fn duplication_conflict_is_detected_in_either_member_order() {
        // Detection must not depend on authored order: the sort that renders
        // cis members in genomic order is gated *off* by this very warning, so
        // both spellings have to be reported or the two orders normalize to
        // each other and idempotency breaks.
        for allele in ["NM_TEST.1:c.[5_9inv;5_6dup]", "NM_TEST.1:c.[5_6dup;5_9inv]"] {
            let (variants, phase) = parse_allele(allele);
            assert_eq!(
                detect_insertion_overlaps(&variants, phase).len(),
                1,
                "expected one warning for {allele}"
            );
        }
    }

    #[test]
    fn lone_duplication_does_not_conflict_with_itself() {
        // A `dup` is registered as both a span and a junction occupant. Its own
        // junction (`gap == end`) is not strictly interior to its own span
        // (`gap < end` fails), and the self-pairing is skipped explicitly, so a
        // duplication alongside an unrelated edit stays clean.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[5_6dup;20A>G]");
        assert!(detect_insertion_overlaps(&variants, phase).is_empty());
    }

    #[test]
    fn duplication_abutting_a_span_edge_no_warning() {
        // `1_4dup` writes at junction 4, which is the *edge* of `5_9inv`
        // (interior junctions are 5..=8), not interior to it. Mirrors
        // `insertions_flanking_a_sub_no_warning` for the dup spelling.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[1_4dup;5_9inv]");
        assert!(detect_insertion_overlaps(&variants, phase).is_empty());
    }

    #[test]
    fn two_duplications_at_one_junction_emit_one_warning() {
        // Both dups write their copy into the slot after position 6, with no
        // defined order between the two payloads — the same ambiguity branch
        // (a) reports for `[5_6insA;5_6insT]`. Branch (b) does *not* cover this
        // pair: the inner dup's junction (6) is the outer dup's end, i.e. its
        // edge, and the interior test is strict (`gap < end`).
        let (variants, phase) = parse_allele("NM_TEST.1:c.[3_6dup;5_6dup]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(edit_kinds, &vec!["dup".to_string(), "dup".to_string()]);
    }

    #[test]
    fn repeat_interior_to_inversion_emits_one_warning() {
        // The genomic axis respells an interior insertion as a repeat
        // expansion rather than a dup, so the `repeat` spelling has to be
        // recognised as a junction occupant for the same reason: `1005_1006A[4]`
        // writes its extra copies at junction 1006, interior to `1005_1009inv`.
        let (variants, phase) = parse_allele("NC_000001.11:g.[1005_1009inv;1005_1006A[4]]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(edit_kinds, &vec!["inv".to_string(), "repeat".to_string()]);
    }

    #[test]
    fn duplication_and_insertion_at_one_junction_emit_one_warning() {
        // Mixed spellings compete for the same slot too, and the warning must
        // name each member by its own kind rather than labelling both `ins`.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[5_6dup;6_7insA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(edit_kinds, &vec!["dup".to_string(), "ins".to_string()]);
    }

    #[test]
    fn two_insertions_at_different_junctions_no_warning() {
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_5insT;8_9insA]");
        assert!(detect_insertion_overlaps(&variants, phase).is_empty());
    }

    #[test]
    fn insertions_flanking_a_sub_no_warning() {
        // `273_274ins` and `274_275ins` are at distinct junctions either side
        // of the substitution at `274`; none overlaps the single-base sub.
        // Mutalyzer normalizes this to `c.274delinsTTA` (no EOVERLAP).
        let (variants, phase) = parse_allele("NM_TEST.1:c.[273_274insT;274G>T;274_275insA]");
        assert!(
            detect_insertion_overlaps(&variants, phase).is_empty(),
            "non-overlapping flanking insertions must not warn: {:?}",
            detect_insertion_overlaps(&variants, phase)
        );
    }

    #[test]
    fn insertion_overlap_only_in_cis() {
        // Trans phase: the edits are on different haplotypes, so coincident
        // junctions are not a conflict.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_5insT];[4_5insA]");
        assert!(detect_insertion_overlaps(&variants, phase).is_empty());
    }

    #[test]
    fn end_to_end_normalize_emits_warning() {
        use crate::normalize::Normalizer;
        use crate::reference::mock::MockProvider;

        let normalizer = Normalizer::new(MockProvider::new());
        let v = parse_hgvs("NC_000001.11:g.[100G>A;100A>C]").expect("parse");
        let result = normalizer
            .normalize_with_diagnostics(&v)
            .expect("normalize");
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.code() == "OVERLAP_CONFLICTING_EDITS"),
            "expected OVERLAP_CONFLICTING_EDITS in warnings, got {:?}",
            result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
        );
        let out = result.result.to_string();
        assert!(
            out.contains("100G>A") && out.contains("100A>C"),
            "expected pass-through, got {out}"
        );
    }

    // --- #1749: one overlap definition, on write footprints ---------------
    //
    // The matrix below is the guard that stops a third opinion coming back.
    // Its `conflict` column is NOT a recording of what ferro does today: it is
    // derived independently, from the write-footprint principle #1411/#1437/#1448
    // established but never centralised —
    //
    //   * `sub`/`delins`/`inv`  write over the span they name, and keep bases there;
    //   * `del`                 writes over its span and keeps nothing;
    //   * `dup` (certain)       reads its span and writes at the junction 3' of it
    //                           (`duplication.md:5`);
    //   * `ins`                 writes at the junction between its two flanks;
    //   * `repeat`              may rewrite its whole span AND lands any expansion
    //                           at the junction 3' of it — its footprint is
    //                           reference-dependent, so it must be read as both.
    //
    // Two members conflict when those footprints collide: base ranges that
    // intersect, one junction claimed twice, or a junction strictly interior to a
    // span whose bases SURVIVE (a pure `del` leaves nothing for an interior
    // junction to be positioned against — #1406).
    //
    // Every row is a pair, so exactly one verdict is well defined for it, and the
    // union of the two detectors must equal that verdict. Where a detector is
    // silent on a shape by design, the OTHER must speak: that is what makes this
    // a test of the definition rather than of either pass's internal division of
    // labour.
    const FOOTPRINT_MATRIX: &[(&str, bool, &str)] = &[
        // --- span x span: footprints are the spans, so plain intersection -----
        ("NC_TEST.1:g.[10_14del;12_16del]", true, "partial overlap"),
        ("NC_TEST.1:g.[10_14del;10_16del]", true, "nested"),
        ("NC_TEST.1:g.[10_14del;10_14del]", true, "coincident"),
        ("NC_TEST.1:g.[10_14del;15_19del]", false, "flush, disjoint"),
        ("NC_TEST.1:g.[10_14inv;12_16inv]", true, "partial overlap"),
        (
            "NC_TEST.1:g.[10_14delinsAA;12_16del]",
            true,
            "partial overlap",
        ),
        ("NC_TEST.1:g.[12G>A;12G>C]", true, "coincident single base"),
        ("NC_TEST.1:g.[12G>A;13G>C]", false, "adjacent single bases"),
        // --- dup: writes at a junction, so its READ span is not its footprint --
        (
            "NC_TEST.1:g.[10_14dup;12_16dup]",
            false,
            "junctions 14 and 16 differ (#1448)",
        ),
        (
            "NC_TEST.1:g.[10_14dup;10_16dup]",
            false,
            "nested reads, junctions 14 and 16 differ",
        ),
        (
            "NC_TEST.1:g.[10_14dup;12_14dup]",
            true,
            "one junction, two writers",
        ),
        (
            "NC_TEST.1:g.[10_14dup;12G>A]",
            false,
            "sub inside the READ span, not the write (#1411)",
        ),
        (
            "NC_TEST.1:g.[10_14dup;14G>A]",
            false,
            "sub at the dup's last base; junction is 3' of it",
        ),
        (
            "NC_TEST.1:g.[10_20inv;14_15dup]",
            true,
            "dup junction interior to a span that keeps its bases",
        ),
        // --- ins: writes at the junction between its flanks --------------------
        (
            "NC_TEST.1:g.[14_15insA;14_15insT]",
            true,
            "one junction, two writers",
        ),
        (
            "NC_TEST.1:g.[14_15insA;15_16insT]",
            false,
            "different junctions",
        ),
        (
            "NC_TEST.1:g.[10_20delinsAA;14_15insT]",
            true,
            "interior junction, bases survive",
        ),
        (
            "NC_TEST.1:g.[10_20del;14_15insT]",
            false,
            "interior junction, del keeps nothing (#1406)",
        ),
        (
            "NC_TEST.1:g.[10_14del;14_15insT]",
            false,
            "junction flush against the 3' edge",
        ),
        (
            "NC_TEST.1:g.[10_14del;9_10insT]",
            false,
            "junction flush against the 5' edge",
        ),
        // --- repeat: whether it writes a span or a junction is arithmetic on the
        // description (`unit x count` against the tract), not a coin toss. Both
        // branches are exercised, because a rule that only ever took one of them
        // would pass a matrix that only ever asked for one.
        //
        // SHRINKING — `A[3]` over a 5-base tract keeps 10-12 and removes 13-14,
        // so the footprint is that trailing range and its bases do not survive.
        (
            "NC_TEST.1:g.[10_14A[3];12_16del]",
            true,
            "removed range 13-14 intersects the deletion",
        ),
        (
            "NC_TEST.1:g.[10_14A[3];13_17del]",
            true,
            "removed range 13-14 intersects the deletion",
        ),
        (
            "NC_TEST.1:g.[10_14A[3];10_11del]",
            false,
            "the deletion sits in the KEPT prefix, untouched reference",
        ),
        // `A[3]` over an 11-base tract removes 13-20. A junction inside a range
        // that is going away has nothing left to be positioned against (#1406) —
        // the same reason a pure `del` exempts one.
        (
            "NC_TEST.1:g.[10_20A[3];14_15insT]",
            false,
            "interior to the REMOVED range, which keeps nothing",
        ),
        // Its junction is not interior to the removed range 13-14 (`gap < end`
        // is strict), so the two writes are disjoint.
        (
            "NC_TEST.1:g.[10_14A[3];12_14dup]",
            false,
            "dup junction 14 sits at the 3' edge, not inside",
        ),
        // GROWING — `A[3]` over a 2-base tract asserts 10-11 unchanged and lands
        // one extra copy at the junction 3' of 11. That is a `dup`'s geometry,
        // which is what makes the two spellings of one shape agree (#1437).
        (
            "NC_TEST.1:g.[10_11A[3];10_11insT]",
            false,
            "repeat writes at junction 11, insertion at 10",
        ),
        (
            "NC_TEST.1:g.[10_11A[3];11_12insT]",
            true,
            "both write at junction 11",
        ),
        (
            "NC_TEST.1:g.[10_11A[3];10_11dup]",
            true,
            "the two spellings of one shape collide identically",
        ),
    ];

    /// The union of the two detectors must equal the write-footprint verdict, on
    /// every row of [`FOOTPRINT_MATRIX`].
    ///
    /// Phrased on observable detector output rather than on the shared predicate
    /// so that it keeps judging behaviour across any later refactor of how the
    /// definition is spelled.
    #[test]
    fn the_detectors_agree_with_the_write_footprint_definition() {
        let mut wrong = Vec::new();
        for (input, expected, why) in FOOTPRINT_MATRIX {
            let (variants, phase) = parse_allele(input);
            let coincident = detect_overlap_conflicts(&variants, phase);
            let junction = detect_insertion_overlaps(&variants, phase);
            let got = !coincident.is_empty() || !junction.is_empty();
            if got != *expected {
                wrong.push(format!(
                    "  {input}\n      expected conflict={expected} ({why}), got {got} \
                     [coincident={}, junction={}]",
                    coincident.len(),
                    junction.len()
                ));
            }
        }
        assert!(
            wrong.is_empty(),
            "{} of {} rows disagree with the write-footprint definition:\n{}",
            wrong.len(),
            FOOTPRINT_MATRIX.len(),
            wrong.join("\n")
        );
    }

    /// The #1448 unit test above says it is asserted at unit level rather than
    /// end to end "because `detect_insertion_overlaps` reports an overlapping
    /// `dup` pair for its own separate reasons (#1437)".
    ///
    /// It does not, and #1437 is the reason it does not: that change made
    /// [`junction_overlaps`] register a `dup`/`repeat` as a junction occupant
    /// and then `continue`, so it is never pushed into `spans`. Two `dup`s at
    /// different junctions therefore match neither branch (a) (the junctions
    /// differ) nor branch (b) (there are no spans). The comment describes
    /// pre-#1437 behaviour and cites #1437 as its cause.
    ///
    /// Pinned so the corrected comment cannot drift back.
    #[test]
    fn overlapping_duplication_spans_are_silent_in_the_junction_detector_too() {
        for input in [
            "NG_012337.1:g.[10_14dup;12_16dup]",
            "NG_012337.1:g.[10_14dup;10_16dup]",
        ] {
            let (variants, phase) = parse_allele(input);
            let warnings = detect_insertion_overlaps(&variants, phase);
            assert!(
                warnings.is_empty(),
                "`{input}`: since #1437 a `dup` is a junction occupant and not a \
                 span, so two `dup`s at different junctions collide in neither \
                 branch; got {warnings:?}"
            );
        }
    }
}
