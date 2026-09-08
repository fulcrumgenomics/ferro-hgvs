//! Arm traits: L1 segmentation, L2 typing (Dial A × Dial B), and monolithic.

use crate::partition::block_ctx::{BlockCtx, FrameContext, Molecule, Provenance};
use crate::partition::output::{Member, Partition, Segment};

/// L1: find segment walls. Pure fn of `(reference, resulting, molecule, provenance)`
/// — no reading frame, no input spelling. `molecule` is one axis fact L1 needs: the
/// payload-coincidence carve-out (`DNA/delins.md:44-47`) is in reach on a DNA axis
/// and never on RNA, so where the walls fall is molecule-dependent (`DNA` collapses
/// a coincidence separation, `RNA` keeps it). Frame stays out — the wall decision
/// keys on molecule, not on coding-ness (see [`Molecule`]). `provenance` is the
/// second axis fact: an individuation claim suppresses the DNA coincidence carve-out
/// (wired in the shared `coincidence` core; every other segmenter ignores it).
pub trait L1Segmenter {
    fn name(&self) -> &str;
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        molecule: Molecule,
        provenance: &Provenance,
    ) -> Vec<Segment>;
}

/// L2: type each segment (Dial A) then merge-back across walls (Dial B).
pub trait L2Typer {
    fn name(&self) -> &str;
    /// Dial A: type one segment into member(s).
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        frame: &FrameContext,
    ) -> Vec<Member>;
    /// Dial B: apply the merge-back carve-outs enabled in `cfg` over the full
    /// member list produced across all segments.
    fn merge_back(
        &self,
        members: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member>;
}

/// Monolithic: decides L1 + L2 jointly (design §5.6).
///
/// **Residual (Steps 2+3 scope note):** monolithic arms remain molecule- and
/// provenance-blind — `partition` takes only `(reference, resulting, frame)`.
/// `MaxSplitMerge` pins `Molecule::Dna` / `Provenance::none()` internally at its
/// one Dial-B call site (`l2_typers.rs`), preserving today's molecule-blind
/// reading exactly. Threading the full `Case` facts through `Monolithic` is a
/// future step, out of Steps 2+3's scope (design §3, Appendix B).
pub trait Monolithic {
    fn name(&self) -> &str;
    fn partition(&self, reference: &[u8], resulting: &[u8], frame: &FrameContext) -> Partition;
    /// Whether this arm takes the canonical coalesce passes (the `delins.md:44-47`
    /// payload-coincidence merge-back). No default: every arm must declare it
    /// (design §6, the `Partitioner::cuts_with_canonical` tripwire). The four
    /// arms implementing this trait today are bakeoff-experimental or
    /// retired-in-Task-5 legacy wrappers, so their values here are not asserted
    /// by any test — Task 3 introduces the legacy-rule wrapper types that carry
    /// the *asserted* canonical-coalesce values.
    fn cuts_with_canonical(&self) -> bool;
}

/// Which merge-back carve-outs are active (design §5.5). C1 = codon gap-1,
/// C2 = net-deletion payload-coincidence, C3 = unequal-length placed-gap.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DialBConfig {
    pub c1: bool,
    pub c2: bool,
    pub c3: bool,
    /// C2's bounded-gap grouper reach and embed tolerance (M4), bundled in the shared
    /// [`CoincidenceParams`] so a sweep cell drives this seam AND L1's W1
    /// (`WallPolicy::coincidence`) from one pair. These two are **house-choice numbers
    /// with no spec grounding** — production says so of its own
    /// `COALESCE_MAX_SEPARATION`/`COALESCE_MISMATCH_BUDGET` (`merge.rs:9024-9028`): the
    /// width and the substitution tolerance are not stated by any clause. They default
    /// to production's values but are config fields precisely so a later sweep can vary
    /// them (the settlement criterion is the confluence objective, not agreement with
    /// production — see the ledger disclosure). Only C2 reads them; C1 and C3 are
    /// single-base by their own antecedents.
    pub coincidence: crate::partition::coincidence::CoincidenceParams,
    /// Extend the C2/C3 payload-coincidence merges to the **RNA** axis (`r.`), not
    /// only DNA. `false` (the default) is the pre-R5 policy the ledger shipped: C2/C3
    /// are DNA-only (`#2155` all-DNA scope, no RNA counterpart), so an RNA
    /// payload-coincidence block stays split. `true` implements the operator-blessed
    /// ruling **R5** (`codon-carve-out-…` sibling; 2026-08-25): the payload-coincidence
    /// carve-out reaches `r.` too, on the anti-provenance ground that the split asserts
    /// conservation of bases whose identity to the payload is coincidental. Protein is
    /// never reached either way. Kept a config flag (not a mutation of
    /// `ledger_current`) so the pre-R5 arms stay in the tournament for comparison.
    pub coincidence_reaches_rna: bool,
}

impl DialBConfig {
    pub fn all_off() -> Self {
        Self::default()
    }
    pub fn ledger_current() -> Self {
        Self {
            c1: true,
            c2: true,
            c3: true,
            ..Self::default()
        }
    }
    /// `ledger_current` with R5 applied: the C2/C3 payload-coincidence merges also
    /// fire on the RNA axis. The op-extract v-next candidate for bakeoff finding #1.
    pub fn ledger_r5() -> Self {
        Self {
            coincidence_reaches_rna: true,
            ..Self::ledger_current()
        }
    }
    /// A stable short label for the JSONL record, e.g. "all-off" / "ledger-current" / "C2".
    pub fn name(&self) -> String {
        match (self.c1, self.c2, self.c3) {
            (false, false, false) => "all-off".into(),
            (true, true, true) => "ledger-current".into(),
            (c1, c2, c3) => [("C1", c1), ("C2", c2), ("C3", c3)]
                .iter()
                .filter(|(_, on)| *on)
                .map(|(n, _)| *n)
                .collect::<Vec<_>>()
                .join("+"),
        }
    }
}

/// The per-case axis facts Dial-B needs, bundled so the `merge_back` signature
/// stays within the lint budget: reading frame (C1), molecule (C2/C3 are
/// DNA-scoped — #2155), and provenance (individuation suppresses C2/C3).
pub struct MergeContext<'a> {
    pub frame: &'a FrameContext,
    pub molecule: Molecule,
    pub provenance: &'a Provenance,
}

impl<'a> MergeContext<'a> {
    /// The context a [`BlockCtx`] implies — what `Strategy::run` threads.
    pub fn of(ctx: &BlockCtx<'a>) -> Self {
        Self {
            frame: ctx.frame,
            molecule: ctx.molecule,
            provenance: ctx.provenance,
        }
    }
}
