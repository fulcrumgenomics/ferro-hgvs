//! The one place the DNA payload-coincidence carve-out and its #1610 placed-gap
//! extension are decided, shared by the L1 wall filter (`arms::disbelieve_dna_
//! coincidences`) and the Dial-B merge-back (`dial_b`'s C2/C3). Decision #1: one
//! shared core, two thin wrappers — so the two realizations cannot drift silently
//! the way they did (the `>=`-vs-`<` equal-length drift that #2155/#1a corrected),
//! while each wrapper can still be tuned independently if a measurement ever calls
//! for it.
//!
//! A caller projects its native type — an L1 `Segment` or a Dial-B `Member` — into
//! [`CoincidencePiece`], a molecule-blind `(ref span, bases supplied)` view, and the
//! predicates below decide from that alone. Ledger scope: net-deletion ONLY
//! (`delins-merge-vs-individual-gap-two-or-more`'s direction scope); equal-length
//! collapse is a different rule (#2174) owned by the typer, never here.
//!
//! The individuation channel (`IndividuationPolicy`) is decided HERE, once, for
//! both wrappers — every predicate below takes the caller's policy as its last
//! parameter and refuses to fire when it suppresses the carve-out
//! (`IndividuationPolicy::suppresses_coincidence_collapse`). A per-wrapper guard
//! is exactly the drift the shared core exists to prevent (decision #1).

use crate::partition::block_ctx::IndividuationPolicy;

/// Maximum unchanged-base separation C2's payload-coincidence merge reaches — the
/// bounded gap that M4 widens single-base to. Production's `COALESCE_MAX_SEPARATION`
/// (`merge.rs:8829`), which production concedes "does not merely lack a source"
/// (`merge.rs:9024-9026`) — the width is a house choice, not a clause. Exposed as
/// [`CoincidenceParams::max_separation`] so a sweep can vary it; this const is only its
/// default. Applies to C2 ONLY: C1's codon exception and C3's placed-gap rule are
/// single-base by their own antecedents.
pub(crate) const COALESCE_MAX_SEPARATION: usize = 8;

/// Substitution budget for C2's embed test — production's `COALESCE_MISMATCH_BUDGET`
/// (`merge.rs:8753`). A widened coincidence merge is licensed only if the merged
/// payload reads out of the span by deletion plus at most this many substitutions.
/// Production concedes this "has no counterpart in the spec at all" (`merge.rs:9028`):
/// it is a house choice, exposed as [`CoincidenceParams::mismatch_budget`] so a sweep can
/// vary it; this const is only its default.
pub(crate) const COALESCE_MISMATCH_BUDGET: usize = 1;

/// The two house-choice numbers of the DNA payload-coincidence merge, bundled so a
/// single sweep cell drives BOTH seams at once. `WallPolicy` (L1's W1) and
/// `DialBConfig` (Dial-B's C2) each embed one, so "one `(max_separation,
/// mismatch_budget)` pair per cell" is a type-level fact rather than a convention two
/// construction sites could silently violate — the same drift the shared predicate
/// core exists to prevent (decision #1). Both fields default to the production-derived
/// consts above.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoincidenceParams {
    /// Maximum unchanged-base separation the bounded-gap grouper joins (`1..=this`).
    pub max_separation: usize,
    /// Substitution budget for the embed test.
    pub mismatch_budget: usize,
}

impl Default for CoincidenceParams {
    fn default() -> Self {
        Self {
            max_separation: COALESCE_MAX_SEPARATION,
            mismatch_budget: COALESCE_MISMATCH_BUDGET,
        }
    }
}

impl CoincidenceParams {
    /// The RATIFIED coincidence parameters (operator decision, 2026-08-25, closing the
    /// Step-7+ deferred-plan item 2 confluence sweep). `max_separation = 5`,
    /// `mismatch_budget = 1` — derived from first principles, owing nothing to
    /// production's `(8, 1)` default:
    ///
    /// * **`mismatch_budget = 1`** is forced by the F2 codon-first ceiling
    ///   (`codon-exception-vs-coincidence-carve-out-precedence`, `budget < 2`) plus the
    ///   fact that `budget = 0` is a degenerate (inert) C2 column.
    /// * **`max_separation = 5`** is the spec-grounded FLOOR: `delins.md:44-47`'s own
    ///   worked example is a chain of gaps 4/5/3, and reproducing its single spanning
    ///   `delins` requires bridging the widest gap, 5 (`bounded_gap_runs_over` needs
    ///   *every* interior gap `<= max_separation`). Measured on the production-WINNING arm
    ///   `trim-only/op-extract/ledger`, the gap-5 gap-bearing coincidence witnesses flip
    ///   from incorrect to correct exactly at 5, on BOTH a coding and a genomic axis (the
    ///   carve-out is DNA-general since #2155). Production's `8` over-reaches with no
    ///   adjudicated-correct benefit.
    ///
    /// **Disclosed open item:** the value *above* the floor has no first-principles
    /// source — the spec and ledger are silent on any upper bound, and no corpus (real or
    /// synthetic) yet contains a gap-6/7/8 witness of the genuine gap-bearing carve-out.
    /// A targeted real-corpus filter (net-deletion + ≥1 gap-bearing member + interior gap
    /// 6-13) is the follow-up. The pure-deletion stratum does NOT bear on this: it is
    /// excluded by the `is_gap_bearing` gate at every `max_separation`, so it argues for
    /// neither raising nor lowering the parameter.
    pub fn ratified() -> Self {
        Self {
            max_separation: 5,
            mismatch_budget: 1,
        }
    }
}

/// Maximal runs of consecutive pieces whose successive gap (unchanged reference bases
/// between them) is in `1..=max_sep`. Keyed on `(ref_start, ref_end)` spans so the
/// two seams that group — Dial-B's C2 (from `Member`s) and L1's W1 pass 1 (from
/// `Segment`s) — share ONE grouper and cannot drift on the separation bound (decision
/// #1). A separation of zero (touching pieces) is not a coincidence wall and does not
/// join a run, exactly as the single-base grouper's `== 1` excludes it. Returns
/// half-open index ranges `[start, end)` into `spans`, only for runs of length ≥ 2.
pub(crate) fn bounded_gap_runs_over(
    spans: &[(usize, usize)],
    max_sep: usize,
) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    let mut i = 0usize;
    while i < spans.len() {
        let mut j = i;
        while j + 1 < spans.len() {
            let gap = spans[j + 1].0.saturating_sub(spans[j].1);
            if (1..=max_sep).contains(&gap) {
                j += 1;
            } else {
                break;
            }
        }
        if j > i {
            runs.push((i, j + 1));
        }
        i = j + 1;
    }
    runs
}

/// One piece of a single-base-separated (or, post-M4, bounded-gap) run, projected
/// from a `Segment` or a `Member`. `ref_start..ref_end` is the reference span it
/// consumes; `alt` is the bases it supplies to the resulting sequence. Owned (not
/// borrowed) because a projected member's content is computed, not a slice — an
/// `Inv` reverse-complements, a `Del` supplies nothing; the allocation is negligible
/// on these tiny runs. `alt_len` is derived so it cannot disagree with `alt`.
pub(crate) struct CoincidencePiece {
    pub ref_start: usize,
    pub ref_end: usize,
    pub alt: Vec<u8>,
}

impl CoincidencePiece {
    fn ref_len(&self) -> usize {
        self.ref_end - self.ref_start
    }
    fn alt_len(&self) -> usize {
        self.alt.len()
    }
    /// Supplies at least one base (not a pure deletion).
    fn supplies_bases(&self) -> bool {
        !self.alt.is_empty()
    }
    /// Supplies bases while consuming a *different* number of reference bases — the
    /// member that owes its content to where the alignment placed a gap
    /// (`delins.md:46`). Mirrors L1's `segment_reaches_carve_out` and Dial-B's
    /// `carries_gap_bearing_insert`, which are the same test.
    fn is_gap_bearing(&self) -> bool {
        !self.alt.is_empty() && self.alt_len() != self.ref_len()
    }
    fn is_pure_deletion(&self) -> bool {
        self.alt.is_empty() && self.ref_len() > 0
    }
}

/// Whether `payload` reads out of `span` by deleting bases and substituting at most
/// `budget` of them — production's `payload_embeds_within_budget` (`merge.rs:8843`),
/// a rolling DP over the minimum substitutions across every embedding. Ported (the
/// bakeoff is isolated from production `merge.rs`), and simplified: the bakeoff's
/// blocks are small, so the production grid-size fallback is unnecessary here.
fn payload_embeds_within_budget(span: &[u8], payload: &[u8], budget: usize) -> bool {
    // Each payload base consumes a span base, so a longer payload cannot embed.
    if payload.len() > span.len() {
        return false;
    }
    const UNREACHABLE: u32 = u32::MAX / 2;
    let m = payload.len();
    // `best[j]` = fewest substitutions to embed `payload[j..]` into the span suffix
    // seen so far; consuming the whole span with payload left over is impossible.
    let mut best = vec![UNREACHABLE; m + 1];
    best[m] = 0;
    for &span_base in span.iter().rev() {
        let mut next = vec![UNREACHABLE; m + 1];
        next[m] = 0; // deleting every remaining span base is free
        for j in (0..m).rev() {
            let deleted = best[j];
            let consumed = best[j + 1].saturating_add(u32::from(span_base != payload[j]));
            next[j] = deleted.min(consumed);
        }
        best = next;
    }
    best[0] as usize <= budget
}

/// The bases the merged run supplies, spelled out: each piece's `alt` plus the
/// unchanged reference between consecutive pieces (the coincidence walls the merge
/// writes literally). Its length is exactly [`merged_alt_len`].
fn merged_payload(run: &[CoincidencePiece], reference: &[u8]) -> Vec<u8> {
    let mut payload = Vec::new();
    let mut cursor = run[0].ref_start;
    for p in run {
        payload.extend_from_slice(&reference[cursor..p.ref_start]);
        payload.extend_from_slice(&p.alt);
        cursor = p.ref_end;
    }
    payload
}

/// Reference bases the merged run consumes: last piece's end minus first's start.
fn merged_ref_len(run: &[CoincidencePiece]) -> usize {
    run[run.len() - 1].ref_end - run[0].ref_start
}

/// Bases the merged run supplies: each piece's `alt_len` plus the unchanged
/// reference between consecutive pieces (the coincidence walls, which the merge
/// spells literally).
fn merged_alt_len(run: &[CoincidencePiece]) -> usize {
    let mut n = 0usize;
    let mut cursor = run[0].ref_start;
    for p in run {
        n += p.ref_start - cursor; // unchanged reference between members
        n += p.alt_len();
        cursor = p.ref_end;
    }
    n
}

/// A run is a net deletion iff it supplies fewer bases than it consumes — the
/// direction scope of `delins-merge-vs-individual-gap-two-or-more`.
pub(crate) fn is_net_deletion(run: &[CoincidencePiece]) -> bool {
    merged_alt_len(run) < merged_ref_len(run)
}

/// A run is a net insertion iff it supplies MORE bases than it consumes
/// (`DNA/duplication.md:90-92` keeps those split, on every wrapper).
fn is_net_insertion(run: &[CoincidencePiece]) -> bool {
    merged_alt_len(run) > merged_ref_len(run)
}

/// Base rule (`delins-recommendation-reach-when-the-input-arrives-split` +
/// `delins-merge-vs-individual-gap-two-or-more`): the run contains a gap-bearing
/// member — so its interior walls are alignment coincidences, not position-anchored —
/// the merged span is a net deletion, AND the merged payload reads out of the span
/// within [`COALESCE_MISMATCH_BUDGET`] substitutions (M4). This is exactly Dial-B's
/// C2 gate, mirroring production's `coalesce_payload_alignment_split` (`merge.rs:9353`):
/// gap-bearing + net-deletion + `payload_embeds_within_budget`. The embed test is what
/// makes a BOUNDED-gap merge principled — with walls up to `COALESCE_MAX_SEPARATION`
/// wide, "gap-bearing + net-deletion" alone would over-merge; the budget requires the
/// span to actually spell the payload. The run's separation bound is enforced upstream
/// by the caller's grouper (Dial-B's `bounded_gap_runs`), matching production's
/// per-pair separation gate.
pub(crate) fn base_carve_out(
    run: &[CoincidencePiece],
    reference: &[u8],
    budget: usize,
    policy: IndividuationPolicy,
) -> bool {
    if policy.suppresses_coincidence_collapse() {
        return false;
    }
    if !(run.iter().any(CoincidencePiece::is_gap_bearing) && is_net_deletion(run)) {
        return false;
    }
    let span = &reference[run[0].ref_start..run[run.len() - 1].ref_end];
    payload_embeds_within_budget(span, &merged_payload(run, reference), budget)
}

/// L1-ONLY variant of the base carve-out, widened from net-deletion to "not a net
/// insertion" — i.e. it also collapses the EQUAL-LENGTH case. This is load-bearing
/// and is why L1's wrapper and Dial-B's C2 legitimately differ (decision #1 keeps
/// two wrappers precisely so they can): a splitting L1 (`AllAlignmentSplitL1`) cuts
/// an equal-length block at its forced-unchanged (minimal-alignment) columns, whose
/// interior matches are shift artifacts, not zero-shift separators (#2174). L1 must
/// re-collapse it here, because it split BEFORE the typer's `op_extract` step-6
/// (#2174) could see the whole block. Dial-B runs post-typer on already-typed
/// members and so uses the strict net-deletion [`base_carve_out`] instead — the
/// equal-length case is resolved by the typer there, never re-merged.
///
/// Two further reasons Dial-B must NOT adopt this widening, beyond the typer
/// pre-emption above: (1) grid attribution — a widened C2 would collapse the
/// equal-length #2155 block even in W1-off / axis-blind cells, contaminating the
/// baseline that isolates the wall dial (see `wall_policy_w1_off_does_not_
/// collapse_the_octamer`); and (2) rule ownership — the equal-length collapse is
/// #2174's rule (typer step 6, or L1 standing in for it when a splitting L1
/// pre-empted the typer), while Dial-B's C2 mirrors ferro's production
/// `split_carries_a_gap_bearing_insert` (#1698), whose ledger direction scope is
/// net-deletion ONLY. Converging the two predicates would blur that boundary in
/// both directions.
pub(crate) fn base_carve_out_or_equal_length(
    run: &[CoincidencePiece],
    policy: IndividuationPolicy,
) -> bool {
    if policy.suppresses_coincidence_collapse() {
        return false;
    }
    run.iter().any(CoincidencePiece::is_gap_bearing) && !is_net_insertion(run)
}

/// #1610 placed-gap extension (`unequal-length-block-a-placed-gap-is-not-a-
/// separation`): a run with no gap-bearing member still collapses when it is a net
/// deletion, some member supplies bases (excludes W58's all-deletion split), every
/// member consumes at least one reference base (excludes a pure insertion), and at
/// most one member is a pure deletion (the "one member wider" gap the aligner
/// placed). This is Dial-B's C3 gate minus its caller-side "the whole partition is
/// one run" condition, which stays with the caller. When a gap-bearing member is
/// present the base rule already fires, so this branch matters only in its absence.
pub(crate) fn placed_gap_extension(run: &[CoincidencePiece], policy: IndividuationPolicy) -> bool {
    if policy.suppresses_coincidence_collapse() {
        return false;
    }
    is_net_deletion(run)
        && run.iter().any(CoincidencePiece::supplies_bases)
        && run.iter().all(|p| p.ref_len() > 0)
        && run.iter().filter(|p| p.is_pure_deletion()).count() <= 1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn piece(rs: usize, re: usize, alt: &[u8]) -> CoincidencePiece {
        CoincidencePiece {
            ref_start: rs,
            ref_end: re,
            alt: alt.to_vec(),
        }
    }

    const NONE: IndividuationPolicy = IndividuationPolicy::Unspecified;

    #[test]
    fn base_rule_collapses_a_gap_bearing_net_deletion() {
        // ref CGAG: [delins CG->C (2->1, gap-bearing)] wall A [del G (1->0)]. Payload
        // "CA" (C + unchanged A), span "CGAG" — net deletion, and "CA" is a subsequence
        // of the span (0 substitutions), so the embed test passes.
        let run = [piece(0, 2, b"C"), piece(3, 4, b"")];
        assert!(base_carve_out(
            &run,
            b"CGAG",
            COALESCE_MISMATCH_BUDGET,
            NONE
        ));
    }

    #[test]
    fn base_rule_declines_a_net_insertion() {
        // [delins A->TTT (1->3, gap-bearing)] wall [sub] — net INSERTION, kept split
        // (declined before the embed test even runs).
        let run = [piece(0, 1, b"TTT"), piece(2, 3, b"X")];
        assert!(!base_carve_out(
            &run,
            b"AGC",
            COALESCE_MISMATCH_BUDGET,
            NONE
        ));
    }

    #[test]
    fn base_rule_declines_an_equal_length_run() {
        // net-0 (equal length): NOT this carve-out's business (that is #2174, in the
        // typer). [delins AC->TTT (2->3)] wall [del X (1->0)] over span 4 -> alt 4.
        let run = [piece(0, 2, b"TTT"), piece(3, 4, b"")];
        assert!(!is_net_deletion(&run));
        assert!(!base_carve_out(
            &run,
            b"AACG",
            COALESCE_MISMATCH_BUDGET,
            NONE
        ));
    }

    #[test]
    fn base_rule_declines_when_the_payload_does_not_embed_within_budget() {
        // M4 embed test: gap-bearing AND net-deletion, but the span cannot spell the
        // payload within one substitution. ref AAAAA, [delins AAA->XY (3->2, gap-
        // bearing)] wall A [del A]. Payload "XYA" needs X->A and Y->A = 2 subs > budget
        // 1, so a naive "gap-bearing + net-deletion" would over-merge and the embed
        // test correctly declines.
        let run = [piece(0, 3, b"XY"), piece(4, 5, b"")];
        assert!(run.iter().any(CoincidencePiece::is_gap_bearing));
        assert!(is_net_deletion(&run));
        assert!(!base_carve_out(
            &run,
            b"AAAAA",
            COALESCE_MISMATCH_BUDGET,
            NONE
        ));
    }

    #[test]
    fn placed_gap_collapses_the_1610_shape() {
        // CGCG->AAC = [delins CG->AA (2->2)] wall [del G (1->0)]: no gap-bearing
        // member, net deletion, one pure deletion -> #1610 collapses.
        let run = [piece(0, 2, b"AA"), piece(3, 4, b"")];
        assert!(!run.iter().any(CoincidencePiece::is_gap_bearing));
        assert!(placed_gap_extension(&run, NONE));
    }

    #[test]
    fn placed_gap_declines_two_pure_deletions() {
        // W58-adjacent: two pure deletions and nothing supplies -> not a coincidence.
        let run = [piece(0, 1, b""), piece(2, 3, b"")];
        assert!(!placed_gap_extension(&run, NONE));
    }

    #[test]
    fn an_individuation_claim_suppresses_every_carve_out() {
        // Runs on which the ∅ policy fires (copied from the passing cases above).
        let net_del = [piece(0, 2, b"C"), piece(3, 4, b"")]; // base rule fires at ∅
        let placed = [piece(0, 2, b"AA"), piece(3, 4, b"")]; // #1610 extension fires at ∅
        for policy in [
            IndividuationPolicy::KeepSeparate,
            IndividuationPolicy::KnownPolymorphism,
        ] {
            assert!(!base_carve_out(
                &net_del,
                b"CGAG",
                COALESCE_MISMATCH_BUDGET,
                policy
            ));
            assert!(!base_carve_out_or_equal_length(&net_del, policy));
            assert!(!placed_gap_extension(&placed, policy));
        }
    }
}
