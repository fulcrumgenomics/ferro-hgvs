//! Dial B: the L2 merge-back carve-outs (design §5.5).
//!
//! Three carve-outs, each a pure function of `(reference, resulting, frame)`,
//! applied as a post-pass over the member list a Dial-A typer produced. Each
//! decides whether a run of members separated by single unchanged bases should be
//! re-merged into one spanning `delins`.
//!
//! **Composition is SEQUENTIAL, C1 -> C2/C3**, mirroring production's pipeline
//! order (`coalesce_coding_frame_separation` at `merge.rs:4104` precedes
//! `coalesce_payload_alignment_split`/`_compensating_gap_split` at `:4245`/`:4330`).
//! C1 runs first per run, and C2/C3 act on its output — not as a mutually-exclusive
//! `C2/C3 else C1`. This is the principled model (a C1 codon merge is a real
//! re-typing later passes must see) and it keeps the sweep's dial attribution
//! near-additive. Each is a SEPARATE pass that re-groups the member list its own way
//! (C1/C3 single-base, C2 a bounded gap post-M4); the full rationale is in
//! [`apply_dial_b`].
//!
//! **Round-trip soundness is preserved by construction.** Every merge replaces a
//! run of members with one member over the *same* reference span whose payload is
//! exactly the bases the run (and the unchanged reference between its members)
//! contributed — so applying the merged partition rebuilds the identical resulting
//! sequence. Nothing here can break the round-trip gate; it can only change how
//! the same sequence is spelled.
//!
//! The three gates, and the ledger record each mirrors:
//!
//! * **C1 — codon gap-1 merge** (`delins-codon-carve-out-gap-one`, widened by
//!   `codon-carve-out-shape-restriction`). On a coding axis, two or more members
//!   separated by one unchanged base that together fall within a single codon —
//!   i.e. together affect one amino acid, a property of the resulting sequence,
//!   *any* edit type — merge. `delins.md:18` / `general.md:34` govern; reaches
//!   only a frame-declaring (coding) axis per
//!   `projection-codon-exception-is-decided-by-the-rendered-axis`.
//! * **C2 — net-deletion payload-coincidence merge**
//!   (`delins-merge-vs-individual-gap-two-or-more`,
//!   `delins-recommendation-reach-when-the-input-arrives-split`). A net-deletion
//!   run (merged payload shorter than its span) that carries a gap-bearing insert
//!   — some member supplies bases while consuming a *different* number of
//!   reference bases — merges. Mirrors ferro's `split_carries_a_gap_bearing_insert`
//!   (`merge.rs`, #1698). Net-deletion is the direction scope from #2155; C2/C3
//!   are explicitly gated on `ctx.molecule == Molecule::Dna` — no RNA counterpart
//!   (`RNA/delins.md:17` governs unqualified; M1).
//! * **C3 — lone unequal-length placed-gap merge**
//!   (`unequal-length-block-a-placed-gap-is-not-a-separation`, #1610). A *lone*,
//!   minimal, unequal-length net-deletion block whose members are single-base
//!   separated, where some member supplies bases (excludes W58's all-deletion
//!   split) and every member renders as a `del`/`delins`/`sub` (the split buys no
//!   higher-priority label) — keep the block whole. Mirrors ferro's
//!   `split_is_a_placed_gap_coincidence` collapse (`merge.rs:5924`). Deliberately
//!   scoped to the case C2 *misses*: a block with no gap-bearing member (its
//!   canonical reproduction `n.2_5delinsAAC`, `CGCG`->`AAC`, splits into
//!   `delins(CG->AA)` + `del(G)`, neither gap-bearing), so `{C2}` and `{C3}`
//!   isolate cleanly in the sweep.

use crate::partition::arm::{DialBConfig, MergeContext};
use crate::partition::block_ctx::{FrameContext, IndividuationPolicy, Molecule};
use crate::partition::coincidence::{self, CoincidencePiece};
use crate::partition::output::{EditKind, Member};

/// Bases a member contributes to the resulting sequence — the exact inverse of
/// `RefApplier` (see `metrics.rs`), so a merge built from these is round-trip
/// exact.
fn member_content(m: &Member, reference: &[u8]) -> Vec<u8> {
    match m.kind {
        EditKind::Identity => reference[m.ref_start..m.ref_end].to_vec(),
        EditKind::Del => Vec::new(),
        EditKind::Ins | EditKind::Delins | EditKind::Dup | EditKind::Sub => m.inserted.clone(),
        EditKind::Inv => revcomp(&reference[m.ref_start..m.ref_end]),
    }
}

fn revcomp(b: &[u8]) -> Vec<u8> {
    b.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => *other,
        })
        .collect()
}

/// Merge a run of (reference-ordered, non-overlapping) members into one member
/// over `[first.ref_start, last.ref_end)`. The payload is each member's content
/// interleaved with the unchanged reference between them, so the merged member
/// denotes exactly the same bases over exactly the same span. Empty payload
/// renders as `Del`; otherwise `Delins` — the spanning form the carve-outs
/// recommend by name (`delins.md:47`).
fn merge_run(run: &[Member], reference: &[u8]) -> Member {
    debug_assert!(!run.is_empty(), "merge_run needs at least one member");
    let first = &run[0];
    let last = &run[run.len() - 1];
    let mut payload = Vec::new();
    let mut cursor = first.ref_start;
    for m in run {
        payload.extend_from_slice(&reference[cursor..m.ref_start]);
        payload.extend(member_content(m, reference));
        cursor = m.ref_end;
    }
    let kind = if payload.is_empty() {
        EditKind::Del
    } else {
        EditKind::Delins
    };
    Member {
        kind,
        ref_start: first.ref_start,
        ref_end: last.ref_end,
        inserted: payload,
    }
}

/// Reference bases the merged run consumes.
fn span_len(run: &[Member]) -> usize {
    run[run.len() - 1].ref_end - run[0].ref_start
}

/// Bases the run supplies to the resulting sequence (sum of member contents plus
/// the unchanged reference between members).
fn payload_len(run: &[Member], reference: &[u8]) -> usize {
    let mut n = 0usize;
    let mut cursor = run[0].ref_start;
    for m in run {
        n += m.ref_start - cursor; // unchanged reference between members
        n += member_content(m, reference).len();
        cursor = m.ref_end;
    }
    n
}

/// A run is a net deletion iff it supplies fewer bases than it consumes.
#[cfg(test)]
fn is_net_deletion(run: &[Member], reference: &[u8]) -> bool {
    payload_len(run, reference) < span_len(run)
}

/// Ferro's `split_carries_a_gap_bearing_insert` (`merge.rs:9032`): some member
/// supplies bases while consuming a *different* number of reference bases. Such a
/// member owes its content to where the alignment placed a gap, which is the
/// coincidence `delins.md:46` describes and `:47` overrides.
#[cfg(test)]
fn carries_gap_bearing_insert(run: &[Member], reference: &[u8]) -> bool {
    run.iter().any(|m| {
        let content = member_content(m, reference);
        let consumed = m.ref_end - m.ref_start;
        !content.is_empty() && content.len() != consumed
    })
}

/// Codon (0-indexed) a reference offset falls in, given the CDS phase at the
/// block's reference start. `same_codon` in `merge.rs` uses `(base - 1) / 3` on
/// the 1-indexed CDS coordinate; here `cds_phase + offset` is the 0-indexed CDS
/// coordinate, so the codon is `(cds_phase + offset) / 3`.
fn codon_of(offset: usize, cds_phase: u8) -> usize {
    (cds_phase as usize + offset) / 3
}

/// The reference offset a member's amino-acid effect anchors on. For a span it is
/// the first consumed base; for a pure insertion (zero-width) it is the insertion
/// point.
fn anchor_offset(m: &Member) -> usize {
    m.ref_start
}

/// The last reference offset a member touches (for codon cohesion). A zero-width
/// insertion touches only its anchor.
fn last_offset(m: &Member) -> usize {
    if m.ref_end > m.ref_start {
        m.ref_end - 1
    } else {
        m.ref_start
    }
}

/// Split members into maximal runs where every consecutive pair is separated by
/// exactly one unchanged reference base. Identity members are ignored (they are
/// unchanged content); a single-base separation may be an unchanged reference base
/// or a length-1 identity member. Returns runs of length >= 2 as `(start, end)`
/// index ranges into `reals`, the non-identity members in reference order.
fn single_base_runs(reals: &[Member]) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    let mut i = 0usize;
    while i < reals.len() {
        let mut j = i;
        while j + 1 < reals.len() && reals[j + 1].ref_start.saturating_sub(reals[j].ref_end) == 1 {
            j += 1;
        }
        if j > i {
            runs.push((i, j + 1));
        }
        i = j + 1;
    }
    runs
}

/// C1 gate helper: is the whole run within a single codon (does it together affect
/// one amino acid)?
fn within_one_codon(run: &[Member], cds_phase: u8) -> bool {
    let lo = codon_of(anchor_offset(&run[0]), cds_phase);
    let hi = codon_of(last_offset(&run[run.len() - 1]), cds_phase);
    lo == hi
}

/// C3 kind guard: every member is either the placed-gap `del` or a member that renders
/// as a `delins` — the split exposes no higher-priority label (`sub`/`dup`/`ins`/`inv`).
/// Mirrors production's `piece_renders_as_delins` (`src/normalize/merge.rs:9740`), whose
/// `!(span == 1 && alt.len() == 1)` clause EXCLUDES a substitution: a `sub` is "a real
/// type the split buys" (`unequal-length-block-a-placed-gap-is-not-a-separation`,
/// condition 4), so a sub-bearing placed-gap run must stay split, exactly as production's
/// `split_is_a_placed_gap_coincidence` keeps it. (Fixed 2026-08-25: the `Sub` arm was a
/// bakeoff-only over-merge — finding #3 `c13-g-shape2` `GACG->TA`, where op-extract's
/// correct `Sub+Del` split was being collapsed to one spanning `delins`.) Kept in the
/// Dial-B WRAPPER (not the
/// molecule-blind shared `placed_gap_extension`, which sees only ref-span + bases
/// supplied and so cannot tell an `inv` from a same-length `delins`): collapsing an
/// inv-bearing run into a `delins` would destroy the `inv` label. L1's `#1610` branch
/// admits an `inv` here (span-based `all(ref_len>0)`); that the two wrappers differ on
/// `inv` is a real, pre-existing divergence, preserved by this refactor rather than
/// silently reconciled.
fn all_delins_shaped(run: &[Member]) -> bool {
    run.iter()
        .all(|m| matches!(m.kind, EditKind::Del | EditKind::Delins))
}

/// Does some member of the run supply at least one base? (Excludes W58's
/// all-deletion split, which must stay split.)
#[cfg(test)]
fn some_member_supplies(run: &[Member], reference: &[u8]) -> bool {
    run.iter().any(|m| !member_content(m, reference).is_empty())
}

/// Apply the enabled Dial-B carve-outs to a Dial-A member list.
///
/// `is_lone` for C3 is "this run is the entire changed content of the case",
/// which the caller establishes: C3 only fires when the whole partition is one
/// single-base-separated run.
pub fn apply_dial_b(
    members: Vec<Member>,
    reference: &[u8],
    _resulting: &[u8],
    ctx: &MergeContext,
    cfg: &DialBConfig,
) -> Vec<Member> {
    if !cfg.c1 && !cfg.c2 && !cfg.c3 {
        return members;
    }

    // Partition into identity members (passed through) and real members (analysed),
    // preserving reference order.
    let mut ordered = members;
    ordered.sort_by_key(|m| m.ref_start);
    let reals: Vec<Member> = ordered
        .iter()
        .filter(|m| m.kind != EditKind::Identity)
        .cloned()
        .collect();
    let identities: Vec<Member> = ordered
        .iter()
        .filter(|m| m.kind == EditKind::Identity)
        .cloned()
        .collect();

    let dna = ctx.molecule == Molecule::Dna;
    // The C2/C3 payload-coincidence merges reach DNA always, and RNA iff R5 is armed
    // (`cfg.coincidence_reaches_rna`). Protein is never reached. Pre-R5 this is just `dna`.
    let coincidence_axis = dna || (cfg.coincidence_reaches_rna && ctx.molecule == Molecule::Rna);
    let pol = ctx.provenance.individuation;

    // Spans that a merge collapsed into one member, `[first.ref_start, last.ref_end)`.
    // Used to drop identities the merge subsumed (M3).
    let mut merged_spans: Vec<(usize, usize)> = Vec::new();

    // THREE SEQUENTIAL PASSES, in production's pipeline order C1 -> C2 -> C3
    // (`coalesce_coding_frame_separation` at `merge.rs:4104`, then
    // `coalesce_payload_alignment_split` at `:4245`, then the single-base placed-gap
    // rule), each RE-GROUPING the member list its own way. One shared
    // `single_base_runs` grouper cannot serve all three once C2 widens to a bounded
    // gap (M4 commit 2b): C2's payload-coincidence merge reaches a separation of up to
    // `COALESCE_MAX_SEPARATION`, while C1's codon exception and C3's placed-gap rule
    // are single-base by their OWN ratified antecedents (`delins-codon-carve-out-gap-
    // one`; `unequal-length-block-a-placed-gap-is-not-a-separation`), not as a leftover
    // of a shared grouper. So each pass owns its grouping. The passes run over tiny
    // member lists (cheap), and toggling a dial maps 1:1 to running a pass, which keeps
    // the sweep's dial attribution clean. Sequential composition is also the principled
    // model — a C1 codon merge is a real re-typing a later pass must see — and the
    // shift-free reduction of production's codon BRACKET (its second codon pass at
    // `:4357` exists only because the 3'-shift moves pieces; the bakeoff has no shift).
    // C1 is frame-only (`projection-codon-exception-is-decided-by-the-rendered-axis`);
    // C2/C3 are DNA carve-outs (#2155's all-DNA scope, no RNA counterpart).
    let mut reals = reals;
    if cfg.c1 {
        reals = pass_c1_codon(&reals, reference, ctx.frame, &mut merged_spans);
    }
    if cfg.c2 && coincidence_axis {
        reals = pass_c2_coincidence(
            &reals,
            reference,
            cfg.coincidence.max_separation,
            cfg.coincidence.mismatch_budget,
            pol,
            &mut merged_spans,
        );
    }
    if cfg.c3 && coincidence_axis {
        reals = pass_c3_placed_gap(&reals, reference, pol, &mut merged_spans);
    }
    let mut out = reals;

    // M3: re-insert identity members, but DROP any the merge subsumed. A merge over
    // `[first.ref_start, last.ref_end)` swallows every unchanged base inside it; a
    // zero-width Identity a typer emitted there would otherwise be re-appended
    // OVERLAPPING the merged member, and `RefApplier` rejects the overlap (the hard
    // round-trip gate fails on a partition that was valid before Dial-B). An identity
    // that merely touches a merged span's boundary, or lies between unmerged members,
    // does not intersect any span's open interval and is kept.
    out.extend(identities.into_iter().filter(|id| {
        !merged_spans
            .iter()
            .any(|&(s, e)| s < id.ref_end && id.ref_start < e)
    }));
    // Sort by (ref_start, ref_end) so a zero-width member (e.g. an Identity spdi
    // emits for an unchanged-but-shifted block, or an Ins) precedes a span that
    // starts at the same reference position. Keying on ref_start alone would leave
    // a zero-width member after a co-starting span, which the applier reads as a
    // spurious overlap — the clcs-inv × spdi-canonical unsoundness.
    out.sort_by_key(|m| (m.ref_start, m.ref_end));
    out
}

/// Project members into the shared `CoincidencePiece` view (ref span + bases supplied).
fn project_pieces(run: &[Member], reference: &[u8]) -> Vec<CoincidencePiece> {
    run.iter()
        .map(|m| CoincidencePiece {
            ref_start: m.ref_start,
            ref_end: m.ref_end,
            alt: member_content(m, reference),
        })
        .collect()
}

/// Like [`single_base_runs`] but adjacency is 1..=`max_sep` unchanged reference bases
/// — C2's bounded-gap grouper (M4), mirroring production's per-pair separation gate
/// (`coalesce_payload_alignment_split` refuses any pair separated by more than
/// `COALESCE_MAX_SEPARATION`). A separation of zero (touching members) is not a
/// coincidence wall and does not join a run, exactly as the single-base grouper's
/// `== 1` excludes it.
fn bounded_gap_runs(reals: &[Member], max_sep: usize) -> Vec<(usize, usize)> {
    let spans: Vec<(usize, usize)> = reals.iter().map(|m| (m.ref_start, m.ref_end)).collect();
    coincidence::bounded_gap_runs_over(&spans, max_sep)
}

/// C1 pass: within each single-base-separated run on a coding frame, merge each
/// maximal codon-cohesive sub-run (`merge_codon_cohesive`). Off a coding frame this
/// pass is the identity. Single-base grouping is C1's own antecedent
/// (`delins-codon-carve-out-gap-one`), never widened.
fn pass_c1_codon(
    reals: &[Member],
    reference: &[u8],
    frame: &FrameContext,
    merged_spans: &mut Vec<(usize, usize)>,
) -> Vec<Member> {
    let FrameContext::Coding {
        cds_start, cds_end, ..
    } = frame
    else {
        return reals.to_vec();
    };
    // The zone model (design §D2): the CDS occupies `[cds_lo, cds_hi)`, flanked by
    // the 5'UTR `[0, cds_lo)` and 3'UTR `[cds_hi, len)`. C1 respects the zones — no
    // codon merge in a UTR, none across the `c.-1`/`c.1` or `c.72`/`c.*1` seams
    // (`general.md:36` "in a coding sequence"). Two facts make this fall out:
    // (1) `grid_phase` pins a codon boundary at `cds_lo` (and at `cds_hi`), so a
    //     UTR base and its adjacent CDS base land in DIFFERENT grid-codons —
    //     `within_one_codon` already refuses to merge across a seam.
    // (2) `merge_codon_cohesive`'s CDS-membership guard blocks a merge whose members
    //     are not all inside `[cds_lo, cds_hi)`, which is what stops a UTR-INTERNAL
    //     pair (same grid-codon, but not a real CDS codon) from merging.
    // `None` bounds mean the boundary is at/beyond the window, so the zone extends
    // to the window edge — on the all-`None` corpus this is `[0, len)`, the whole
    // window, and the guard is inert (every existing coding row is unmoved).
    let grid = frame.grid_phase().expect("a coding frame has a grid phase");
    let cds_lo = cds_start.unwrap_or(0);
    let cds_hi = cds_end.unwrap_or(reference.len());
    let mut out = Vec::new();
    let mut next = 0usize;
    for &(start, end) in &single_base_runs(reals) {
        out.extend(reals[next..start].iter().cloned());
        next = end;
        out.extend(merge_codon_cohesive(
            &reals[start..end],
            reference,
            grid,
            cds_lo,
            cds_hi,
            merged_spans,
        ));
    }
    out.extend(reals[next..].iter().cloned());
    out
}

/// C2 pass: merge each bounded-gap run (separation 1..=`COALESCE_MAX_SEPARATION`) that
/// is a net-deletion payload coincidence whose payload embeds in the span within the
/// mismatch budget (`base_carve_out`, decided by the shared core). Coalesces a
/// MULTI-member run only. DNA-gated by the caller.
fn pass_c2_coincidence(
    reals: &[Member],
    reference: &[u8],
    max_separation: usize,
    mismatch_budget: usize,
    pol: IndividuationPolicy,
    merged_spans: &mut Vec<(usize, usize)>,
) -> Vec<Member> {
    let mut out = Vec::new();
    let mut next = 0usize;
    for &(start, end) in &bounded_gap_runs(reals, max_separation) {
        out.extend(reals[next..start].iter().cloned());
        next = end;
        let run = &reals[start..end];
        let pieces = project_pieces(run, reference);
        if run.len() >= 2 && coincidence::base_carve_out(&pieces, reference, mismatch_budget, pol) {
            let merged = merge_run(run, reference);
            merged_spans.push((merged.ref_start, merged.ref_end));
            out.push(merged);
        } else {
            out.extend(run.iter().cloned());
        }
    }
    out.extend(reals[next..].iter().cloned());
    out
}

/// C3 pass: the #1610 lone placed-gap extension. Fires only when the ENTIRE changed
/// content is one single-base-separated run (`whole_is_one_run`, computed here from the
/// members this pass receives), every member is `del`/`delins`/`sub` (no
/// higher-priority label to destroy), and the shared `placed_gap_extension` holds.
/// Single-base by its own antecedent
/// (`unequal-length-block-a-placed-gap-is-not-a-separation`) — never widened. Its
/// production analog is the split-time `split_is_a_placed_gap_coincidence` rule
/// (`merge.rs`, ~:5980), NOT `coalesce_compensating_gap_split` (a distinct,
/// density-licensed mechanism with no bakeoff dial today).
fn pass_c3_placed_gap(
    reals: &[Member],
    reference: &[u8],
    pol: IndividuationPolicy,
    merged_spans: &mut Vec<(usize, usize)>,
) -> Vec<Member> {
    let runs = single_base_runs(reals);
    let whole_is_one_run = runs.len() == 1 && runs[0] == (0, reals.len());
    let mut out = Vec::new();
    let mut next = 0usize;
    for &(start, end) in &runs {
        out.extend(reals[next..start].iter().cloned());
        next = end;
        let run = &reals[start..end];
        let pieces = project_pieces(run, reference);
        if run.len() >= 2
            && whole_is_one_run
            && all_delins_shaped(run)
            && coincidence::placed_gap_extension(&pieces, pol)
        {
            let merged = merge_run(run, reference);
            merged_spans.push((merged.ref_start, merged.ref_end));
            out.push(merged);
        } else {
            out.extend(run.iter().cloned());
        }
    }
    out.extend(reals[next..].iter().cloned());
    out
}

/// C1: within one single-base-separated run, merge each maximal sub-run whose
/// members all fall in one codon. Sub-runs that do not cohere (span two codons)
/// pass through unmerged.
fn merge_codon_cohesive(
    run: &[Member],
    reference: &[u8],
    cds_phase: u8,
    cds_lo: usize,
    cds_hi: usize,
    merged_spans: &mut Vec<(usize, usize)>,
) -> Vec<Member> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < run.len() {
        let mut j = i;
        // Explicit gap-1 adjacency guard on each hop, in addition to the codon-hull
        // test. For base-consuming members `within_one_codon` already caps the gap at
        // 1 by pigeonhole (a 3-base codon holds no two separated changed offsets more
        // than one base apart), but a zero-width `Ins` is a POINT to `anchor_offset`/
        // `last_offset`: `Ins(0,0) + Del(2,3)` is two bases apart yet shares a codon
        // and passes the net-0 gate. Post-M4 widening (gap up to 8) that pair would
        // merge, past `delins-codon-carve-out-gap-one`'s gap-one antecedent — so bound
        // the hop here. Whether that gap-2 compensated pair SHOULD merge is an open
        // ledger question, deliberately not decided by a silent widening. At gap-1
        // grouping every intra-run separation is already exactly 1, so this is
        // zero-effect until the grouper widens.
        while j + 1 < run.len()
            && run[j + 1].ref_start.saturating_sub(run[j].ref_end) <= 1
            && within_one_codon(&run[i..=j + 1], cds_phase)
        {
            j += 1;
        }
        let sub = &run[i..=j];
        // Zone guard (design §D2): every member of the sub-run must lie inside the
        // CDS `[cds_lo, cds_hi)`. A span member `[s, e)` is inside iff `cds_lo <= s`
        // and `last_offset < cds_hi`; a zero-width insertion at `p` iff
        // `cds_lo <= p < cds_hi`. This blocks a UTR-internal pair (which shares a
        // grid-codon but not a real CDS codon) from merging; seam-crossing pairs are
        // already blocked by `within_one_codon` (the grid boundary sits at each
        // seam). On the all-`None` corpus `[cds_lo, cds_hi) == [0, len)`, so the
        // guard is vacuously true and C1 is unchanged.
        let sub_in_cds = sub
            .iter()
            .all(|m| cds_lo <= m.ref_start && last_offset(m) < cds_hi);
        // M2: the codon carve-out (`codon-carve-out-shape-restriction`) requires the
        // members to TOGETHER AFFECT ONE AMINO ACID — a property of the resulting
        // sequence, not merely of their reference positions. A net-length-changing
        // sub-run inside one codon FRAMESHIFTS every downstream amino acid, so it does
        // not qualify however cohesive its positions are. Gate the merge on net-0
        // (`payload_len == span_len`) in addition to codon cohesion; a frameshifting
        // sub-run passes through as individual members.
        if j > i && sub_in_cds && payload_len(sub, reference) == span_len(sub) {
            let merged = merge_run(sub, reference);
            merged_spans.push((merged.ref_start, merged.ref_end));
            out.push(merged);
        } else {
            out.extend(sub.iter().cloned());
        }
        i = j + 1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::block_ctx::{IndividuationPolicy, Molecule, Provenance};

    /// The pre-threading reading: DNA molecule, ∅ provenance.
    static NONE: Provenance = Provenance {
        individuation: IndividuationPolicy::Unspecified,
        examined: Vec::new(),
    };
    fn dna(frame: &FrameContext) -> MergeContext<'_> {
        MergeContext {
            frame,
            molecule: Molecule::Dna,
            provenance: &NONE,
        }
    }

    /// ∅ provenance, RNA molecule — the M1 negative-control context.
    fn rna(frame: &FrameContext) -> MergeContext<'_> {
        MergeContext {
            frame,
            molecule: Molecule::Rna,
            provenance: &NONE,
        }
    }

    static KEEP: Provenance = Provenance {
        individuation: IndividuationPolicy::KeepSeparate,
        examined: Vec::new(),
    };
    /// DNA molecule, `KeepSeparate` provenance — the individuation-suppression
    /// context.
    fn dna_keep(frame: &FrameContext) -> MergeContext<'_> {
        MergeContext {
            frame,
            molecule: Molecule::Dna,
            provenance: &KEEP,
        }
    }

    fn m(kind: EditKind, rs: usize, re: usize, ins: &str) -> Member {
        Member {
            kind,
            ref_start: rs,
            ref_end: re,
            inserted: ins.bytes().collect(),
        }
    }

    // Applies a member list the way RefApplier does, so tests can assert
    // round-trip equivalence of a merge locally.
    fn apply(reference: &[u8], members: &[Member]) -> Vec<u8> {
        let mut out = Vec::new();
        let mut cursor = 0usize;
        let mut ms = members.to_vec();
        ms.sort_by_key(|m| m.ref_start);
        for mm in &ms {
            out.extend_from_slice(&reference[cursor..mm.ref_start]);
            out.extend(member_content(mm, reference));
            cursor = mm.ref_end;
        }
        out.extend_from_slice(&reference[cursor..]);
        out
    }

    #[test]
    fn all_off_is_identity() {
        let reference = b"ACGTACGT";
        let members = vec![m(EditKind::Sub, 1, 2, "T"), m(EditKind::Del, 3, 4, "")];
        let out = apply_dial_b(
            members.clone(),
            reference,
            b"",
            &dna(&FrameContext::NonCoding),
            &DialBConfig::all_off(),
        );
        assert_eq!(out, members);
    }

    /// W1/C2 parity pin (M4 2c): the L1 wall seam (`arms::collapse_bounded_net_del`)
    /// and the Dial-B merge-back seam (C2 via `apply_dial_b`) make the SAME merge
    /// decision on the SAME net-deletion run at the SAME params — they share one
    /// grouper (`bounded_gap_runs_over`) and one predicate (`base_carve_out`). This
    /// guards the drift a future edit could introduce by changing one seam's
    /// predicate/grouper call but not the other's. Two runs: one merges (2-base
    /// net-deletion, payload embeds) and one declines (embed budget exceeded); both
    /// seams must agree in both directions, and on the merged span.
    #[test]
    fn w1_and_c2_agree_on_the_same_net_deletion_run() {
        use crate::partition::arms::arms::collapse_bounded_net_del;
        use crate::partition::coincidence::CoincidenceParams;
        use crate::partition::output::Segment;

        let params = CoincidenceParams::default();
        let c2_cfg = DialBConfig {
            c1: false,
            c2: true,
            c3: false,
            ..Default::default()
        };

        // Case 1 — CGGAG -> CGA: [delins CG->C; 2-base wall "GA"; del G], net-deletion,
        // "CGA" embeds in "CGGAG" with 0 substitutions. BOTH seams merge over ref [0,5).
        {
            let reference = b"CGGAG";
            let resulting = b"CGA";
            let members = vec![m(EditKind::Delins, 0, 2, "C"), m(EditKind::Del, 4, 5, "")];
            assert_eq!(apply(reference, &members), resulting);
            let segments = vec![
                Segment {
                    ref_start: 0,
                    ref_end: 2,
                    res_start: 0,
                    res_end: 1,
                },
                Segment {
                    ref_start: 4,
                    ref_end: 5,
                    res_start: 3,
                    res_end: 3,
                },
            ];
            let c2 = apply_dial_b(
                members,
                reference,
                resulting,
                &dna(&FrameContext::NonCoding),
                &c2_cfg,
            );
            let w1 = collapse_bounded_net_del(
                segments,
                reference,
                resulting,
                params,
                IndividuationPolicy::Unspecified,
            );
            assert_eq!(c2.len(), 1, "C2 merges the 2-base net-deletion coincidence");
            assert_eq!(w1.len(), 1, "W1 merges it too (parity)");
            assert_eq!((c2[0].ref_start, c2[0].ref_end), (0, 5));
            assert_eq!((w1[0].ref_start, w1[0].ref_end), (0, 5));
        }

        // Case 2 — AAAAA -> XYA: [delins AAA->XY; gap1; del A], net-deletion but "XYA"
        // needs 2 substitutions to embed in "AAAAA" (> budget 1). BOTH seams DECLINE.
        {
            let reference = b"AAAAA";
            let resulting = b"XYA";
            let members = vec![m(EditKind::Delins, 0, 3, "XY"), m(EditKind::Del, 4, 5, "")];
            assert_eq!(apply(reference, &members), resulting);
            let segments = vec![
                Segment {
                    ref_start: 0,
                    ref_end: 3,
                    res_start: 0,
                    res_end: 2,
                },
                Segment {
                    ref_start: 4,
                    ref_end: 5,
                    res_start: 3,
                    res_end: 3,
                },
            ];
            let c2 = apply_dial_b(
                members,
                reference,
                resulting,
                &dna(&FrameContext::NonCoding),
                &c2_cfg,
            );
            let w1 = collapse_bounded_net_del(
                segments,
                reference,
                resulting,
                params,
                IndividuationPolicy::Unspecified,
            );
            assert_eq!(c2.len(), 2, "C2 declines: embed budget exceeded");
            assert_eq!(w1.len(), 2, "W1 declines too (parity)");
        }
    }

    #[test]
    fn c3_merges_the_1610_lone_unequal_length_block() {
        // reference CGCG (offsets 0..4) -> AAC. A cutter splits it into
        // delins(CG->AA) at 0..2 and del(G) at 3..4, with the unchanged C at
        // offset 2 as the single-base separation. C3 keeps the block whole.
        let reference = b"CGCG";
        let split = vec![m(EditKind::Delins, 0, 2, "AA"), m(EditKind::Del, 3, 4, "")];
        assert_eq!(apply(reference, &split), b"AAC");
        // Neither member is gap-bearing (2->2, 1->0), so C2 does NOT fire.
        assert!(!carries_gap_bearing_insert(&split, reference));

        let c3 = DialBConfig {
            c1: false,
            c2: false,
            c3: true,
            ..Default::default()
        };
        let out = apply_dial_b(
            split.clone(),
            reference,
            b"AAC",
            &dna(&FrameContext::NonCoding),
            &c3,
        );
        assert_eq!(out.len(), 1, "block kept whole");
        assert_eq!(out[0].kind, EditKind::Delins);
        assert_eq!(out[0].ref_start, 0);
        assert_eq!(out[0].ref_end, 4);
        assert_eq!(out[0].inserted, b"AAC");
        // Round-trip preserved.
        assert_eq!(apply(reference, &out), b"AAC");

        // C2 alone leaves it split (no gap-bearing member).
        let c2 = DialBConfig {
            c1: false,
            c2: true,
            c3: false,
            ..Default::default()
        };
        let out2 = apply_dial_b(
            split.clone(),
            reference,
            b"AAC",
            &dna(&FrameContext::NonCoding),
            &c2,
        );
        assert_eq!(out2, split);
    }

    #[test]
    fn c2_merges_a_gap_bearing_net_deletion_split() {
        // reference ACGTA (0..5) -> ATA. delins(ACGT->AT) is gap-bearing
        // (4 consumed, 2 supplied), a net deletion. Split it as delins(CG->"") is
        // awkward; use a two-member single-base-separated split:
        //   delins(AC->A) at 0..2 (gap-bearing: 2->1), unchanged G at 2,
        //   del(T) at 3..4.
        let reference = b"ACGTA";
        let split = vec![m(EditKind::Delins, 0, 2, "A"), m(EditKind::Del, 3, 4, "")];
        assert_eq!(apply(reference, &split), b"AGA");
        assert!(carries_gap_bearing_insert(&split, reference));
        assert!(is_net_deletion(&split, reference));

        let c2 = DialBConfig {
            c1: false,
            c2: true,
            c3: false,
            ..Default::default()
        };
        let out = apply_dial_b(
            split.clone(),
            reference,
            b"AGA",
            &dna(&FrameContext::NonCoding),
            &c2,
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].kind, EditKind::Delins);
        assert_eq!(out[0].ref_start, 0);
        assert_eq!(out[0].ref_end, 4);
        assert_eq!(apply(reference, &out), b"AGA");
    }

    #[test]
    fn c2_does_not_fire_on_an_rna_molecule() {
        // Same gap-bearing net-deletion split DNA merges (c2_merges_a_gap_bearing_
        // net_deletion_split): on RNA the separation is genuine — RNA/delins.md:17
        // governs, there is no :47 counterpart (M1).
        let reference = b"ACGTA";
        let split = vec![m(EditKind::Delins, 0, 2, "A"), m(EditKind::Del, 3, 4, "")];
        let c2 = DialBConfig {
            c1: false,
            c2: true,
            c3: false,
            ..Default::default()
        };
        let out = apply_dial_b(
            split.clone(),
            reference,
            b"AGA",
            &rna(&FrameContext::NonCoding),
            &c2,
        );
        assert_eq!(out, split, "RNA keeps the split C2 would merge on DNA");
    }

    #[test]
    fn c3_does_not_fire_on_an_rna_molecule() {
        let reference = b"CGCG";
        let split = vec![m(EditKind::Delins, 0, 2, "AA"), m(EditKind::Del, 3, 4, "")];
        let c3 = DialBConfig {
            c1: false,
            c2: false,
            c3: true,
            ..Default::default()
        };
        let out = apply_dial_b(
            split.clone(),
            reference,
            b"AAC",
            &rna(&FrameContext::NonCoding),
            &c3,
        );
        assert_eq!(out, split, "RNA keeps the #1610 block split");
    }

    #[test]
    fn keep_separate_provenance_suppresses_c2_on_dna() {
        let reference = b"ACGTA";
        let split = vec![m(EditKind::Delins, 0, 2, "A"), m(EditKind::Del, 3, 4, "")];
        let c2 = DialBConfig {
            c1: false,
            c2: true,
            c3: false,
            ..Default::default()
        };
        let out = apply_dial_b(
            split.clone(),
            reference,
            b"AGA",
            &dna_keep(&FrameContext::NonCoding),
            &c2,
        );
        assert_eq!(
            out, split,
            "KeepSeparate suppresses the carve-out even on DNA"
        );
    }

    #[test]
    fn c1_merges_within_a_codon_only_on_a_coding_frame() {
        // reference AAAAAA. cds_phase 0 => codon 0 = offsets 0..3, codon 1 = 3..6.
        // Two subs at offsets 0 and 2 (one unchanged base at 1 between them), both
        // in codon 0 => C1 merges. NonCoding => untouched.
        let reference = b"AAAAAA";
        let split = vec![m(EditKind::Sub, 0, 1, "C"), m(EditKind::Sub, 2, 3, "G")];
        assert_eq!(apply(reference, &split), b"CAGAAA");

        let c1 = DialBConfig {
            c1: true,
            c2: false,
            c3: false,
            ..Default::default()
        };
        let coding = FrameContext::Coding {
            cds_phase: 0,
            cds_start: None,
            cds_end: None,
        };
        let merged = apply_dial_b(split.clone(), reference, b"CAGAAA", &dna(&coding), &c1);
        assert_eq!(merged.len(), 1, "within one codon => merged");
        assert_eq!(merged[0].ref_start, 0);
        assert_eq!(merged[0].ref_end, 3);
        assert_eq!(apply(reference, &merged), b"CAGAAA");

        // NonCoding: C1 cannot fire.
        let noncoding = apply_dial_b(
            split.clone(),
            reference,
            b"CAGAAA",
            &dna(&FrameContext::NonCoding),
            &c1,
        );
        assert_eq!(noncoding, split);
    }

    #[test]
    fn c1_does_not_merge_across_a_codon_boundary() {
        // subs at offsets 2 and 4 (unchanged at 3), phase 0: offset 2 in codon 0,
        // offset 4 in codon 1 => different codons => NOT merged.
        let reference = b"AAAAAA";
        let split = vec![m(EditKind::Sub, 2, 3, "C"), m(EditKind::Sub, 4, 5, "G")];
        let c1 = DialBConfig {
            c1: true,
            c2: false,
            c3: false,
            ..Default::default()
        };
        let out = apply_dial_b(
            split.clone(),
            reference,
            b"",
            &dna(&FrameContext::Coding {
                cds_phase: 0,
                cds_start: None,
                cds_end: None,
            }),
            &c1,
        );
        assert_eq!(out, split);
    }

    // ---- T4: C1 zone-awareness (design §D2) ----

    /// A gap-1 codon-cohesive net-0 pair sitting entirely in the 5'UTR does NOT
    /// merge: it shares a grid-codon, but not a real CDS codon.
    #[test]
    fn c1_does_not_merge_a_utr_internal_pair() {
        // cds_start = 3 => 5'UTR is [0,3); grid phase (3 - 3%3) % 3 = 0. The pair at
        // offsets 0 and 2 shares grid-codon 0, but both are in the 5'UTR.
        let reference = b"AAAAAAAAA";
        let split = vec![m(EditKind::Sub, 0, 1, "C"), m(EditKind::Sub, 2, 3, "G")];
        let c1 = DialBConfig {
            c1: true,
            c2: false,
            c3: false,
            ..Default::default()
        };
        let f = FrameContext::coding(0, Some(3), Some(9));
        let out = apply_dial_b(split.clone(), reference, b"CAGAAAAAA", &dna(&f), &c1);
        assert_eq!(out, split, "a 5'UTR-internal pair is not a CDS codon merge");
    }

    /// A gap-1 pair straddling the c.-1/c.1 seam (cds_start) does not merge: the
    /// grid boundary sits at cds_start, so the two members fall in different codons.
    #[test]
    fn c1_does_not_merge_across_the_cds_start_seam() {
        // cds_start = 3, grid 0. Sub at 1 (5'UTR, codon 0) and Sub at 3 (CDS, codon 1),
        // separated by the unchanged base at 2 (gap 1).
        let reference = b"AAAAAAAAA";
        let split = vec![m(EditKind::Sub, 1, 2, "C"), m(EditKind::Sub, 3, 4, "G")];
        let c1 = DialBConfig {
            c1: true,
            c2: false,
            c3: false,
            ..Default::default()
        };
        let f = FrameContext::coding(0, Some(3), Some(9));
        let out = apply_dial_b(split.clone(), reference, b"", &dna(&f), &c1);
        assert_eq!(out, split, "no merge across the c.-1/c.1 seam");
    }

    /// A gap-1 pair straddling the c.72/c.*1 seam (cds_end) does not merge.
    #[test]
    fn c1_does_not_merge_across_the_cds_end_seam() {
        // cds_end = 6 (cds_start None, so the CDS runs from the window start), grid 0.
        // Sub at 4 (CDS, codon 1) and Sub at 6 (3'UTR, codon 2), gap-1 (unchanged 5).
        let reference = b"AAAAAAAAA";
        let split = vec![m(EditKind::Sub, 4, 5, "C"), m(EditKind::Sub, 6, 7, "G")];
        let c1 = DialBConfig {
            c1: true,
            c2: false,
            c3: false,
            ..Default::default()
        };
        let f = FrameContext::coding(0, None, Some(6));
        let out = apply_dial_b(split.clone(), reference, b"", &dna(&f), &c1);
        assert_eq!(out, split, "no merge across the c.72/c.*1 seam");
    }

    /// F1 (decided): the terminal/stop codon is inside `[cds_start, cds_end)`, so a
    /// net-0 pair within it merges with no translation — `cds_end` is placed AFTER it.
    #[test]
    fn c1_merges_a_pair_within_the_stop_codon() {
        // CDS [0,9); codon 2 (offsets 6..9) is the terminal codon — TAA here. Two
        // subs within it, gap-1 (unchanged base at 7).
        let reference = b"AAAAAATAA";
        let split = vec![m(EditKind::Sub, 6, 7, "C"), m(EditKind::Sub, 8, 9, "G")];
        assert_eq!(apply(reference, &split), b"AAAAAACAG");
        let c1 = DialBConfig {
            c1: true,
            c2: false,
            c3: false,
            ..Default::default()
        };
        let f = FrameContext::coding(0, Some(0), Some(9));
        let out = apply_dial_b(split.clone(), reference, b"AAAAAACAG", &dna(&f), &c1);
        assert_eq!(out.len(), 1, "stop-codon pair merges (F1)");
        assert_eq!((out[0].ref_start, out[0].ref_end), (6, 9));
    }

    /// The codon grid is derived from `cds_start` (`grid_phase`): an in-CDS gap-1
    /// net-0 pair within one codon merges at every phase.
    #[test]
    fn c1_merges_in_cds_at_all_three_grid_phases() {
        let reference = b"AAAAAAAAA";
        let c1 = DialBConfig {
            c1: true,
            c2: false,
            c3: false,
            ..Default::default()
        };
        // (cds_start, phase, off_a, off_b): the pair lies inside one grid-codon 1.
        //   cds_start 3 -> phase 0, codon 1 = [3,6): offsets 3,5
        //   cds_start 1 -> phase 2, codon 1 = [1,4): offsets 1,3
        //   cds_start 2 -> phase 1, codon 1 = [2,5): offsets 2,4
        for (cds_start, phase, a, b) in [(3usize, 0u8, 3usize, 5usize), (1, 2, 1, 3), (2, 1, 2, 4)]
        {
            let split = vec![
                m(EditKind::Sub, a, a + 1, "C"),
                m(EditKind::Sub, b, b + 1, "G"),
            ];
            let f = FrameContext::coding(phase, Some(cds_start), Some(9));
            let out = apply_dial_b(split.clone(), reference, b"", &dna(&f), &c1);
            assert_eq!(out.len(), 1, "phase {phase}: in-CDS codon pair merges");
            assert_eq!((out[0].ref_start, out[0].ref_end), (a, b + 1));
        }
    }

    /// C2 is zone-BLIND (design §D2): the payload-coincidence carve-out is
    /// axis-scoped, not CDS-scoped, so a net-deletion coincidence in the 5'UTR of a
    /// coding window still merges.
    #[test]
    fn c2_fires_in_the_utr_of_a_coding_window() {
        // The c2-merges shape (reference ACGTA), in a longer window whose CDS starts
        // at offset 5 — so the whole run sits in the 5'UTR. C2 merges it regardless.
        let reference = b"ACGTAAAAA";
        let split = vec![m(EditKind::Delins, 0, 2, "A"), m(EditKind::Del, 3, 4, "")];
        let c2 = DialBConfig {
            c1: false,
            c2: true,
            c3: false,
            ..Default::default()
        };
        let f = FrameContext::coding(1, Some(5), Some(9));
        let out = apply_dial_b(split.clone(), reference, b"", &dna(&f), &c2);
        assert_eq!(out.len(), 1, "C2 ignores the zone");
        assert_eq!((out[0].ref_start, out[0].ref_end), (0, 4));
    }

    #[test]
    fn w58_all_deletion_split_stays_split_under_c3() {
        // A net-deletion block whose members are all pure deletions: no member
        // supplies a base, so C3 must not merge it (W58's shape).
        let reference = b"ACACAC";
        let split = vec![m(EditKind::Del, 0, 1, ""), m(EditKind::Del, 2, 3, "")];
        assert!(!some_member_supplies(&split, reference));
        let c3 = DialBConfig {
            c1: false,
            c2: false,
            c3: true,
            ..Default::default()
        };
        let out = apply_dial_b(
            split.clone(),
            reference,
            b"",
            &dna(&FrameContext::NonCoding),
            &c3,
        );
        assert_eq!(out, split, "all-deletion split is not a coincidence merge");
    }

    #[test]
    fn c1_does_not_merge_a_frameshifting_run_in_one_codon() {
        // M2: reference AAAAAA, cds_phase 0. Sub(0,1,"C") at codon 0 and Del(2,3) at
        // codon 0 are gap-1 and codon-cohesive, but the run is a NET DELETION
        // (payload "CA" = 2, span 3), so it frameshifts every downstream amino acid.
        // The codon carve-out must NOT merge it — it stays split.
        let reference = b"AAAAAA";
        let split = vec![m(EditKind::Sub, 0, 1, "C"), m(EditKind::Del, 2, 3, "")];
        // Precondition: the run is a net deletion (frameshift), so M2's gate applies.
        assert_ne!(payload_len(&split, reference), span_len(&split));
        let c1 = DialBConfig {
            c1: true,
            c2: false,
            c3: false,
            ..Default::default()
        };
        let out = apply_dial_b(
            split.clone(),
            reference,
            b"CAAAA",
            &dna(&FrameContext::Coding {
                cds_phase: 0,
                cds_start: None,
                cds_end: None,
            }),
            &c1,
        );
        assert_eq!(
            out, split,
            "a frameshifting codon-cohesive run is not merged"
        );
        assert!(
            !out.iter().any(|mm| mm.kind == EditKind::Delins),
            "no delins merge: {out:?}",
        );
    }

    #[test]
    fn c3_drops_an_identity_subsumed_by_the_merged_span() {
        // M3: reference CGCG -> AAC, split as delins(CG->AA) + Del(G) with a zero-width
        // Identity(3,3) sitting INSIDE the block. C3 merges the two reals into
        // Delins(0,4,"AAC"); the identity would otherwise be re-appended overlapping
        // that span, and RefApplier would reject the overlap (round-trip gate fails).
        // The subsumed identity must be dropped.
        let reference = b"CGCG";
        let members = vec![
            m(EditKind::Delins, 0, 2, "AA"),
            m(EditKind::Identity, 3, 3, ""),
            m(EditKind::Del, 3, 4, ""),
        ];
        let c3 = DialBConfig {
            c1: false,
            c2: false,
            c3: true,
            ..Default::default()
        };
        let out = apply_dial_b(
            members,
            reference,
            b"AAC",
            &dna(&FrameContext::NonCoding),
            &c3,
        );
        assert_eq!(out.len(), 1, "identity dropped, one merged member: {out:?}");
        assert_eq!(out[0].kind, EditKind::Delins);
        assert_eq!((out[0].ref_start, out[0].ref_end), (0, 4));
        assert!(
            !out.iter().any(|mm| mm.kind == EditKind::Identity),
            "the subsumed identity is gone: {out:?}",
        );
        // No overlap remains, so the merged partition reconstructs.
        assert_eq!(apply(reference, &out), b"AAC");
    }

    #[test]
    fn an_identity_between_unmerged_members_is_kept() {
        // Guard the other direction: with no merge, a separating identity survives.
        let reference = b"ACGT";
        let members = vec![
            m(EditKind::Sub, 0, 1, "T"),
            m(EditKind::Identity, 1, 2, ""),
            m(EditKind::Sub, 3, 4, "A"),
        ];
        let out = apply_dial_b(
            members.clone(),
            reference,
            b"TCGA",
            &dna(&FrameContext::NonCoding),
            &DialBConfig::all_off(),
        );
        assert!(
            out.iter().any(|mm| mm.kind == EditKind::Identity),
            "a non-subsumed identity is kept: {out:?}",
        );
    }
}
