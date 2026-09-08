//! The neutral per-block context an arm consumes, plus the axis facts it carries.
//!
//! A [`BlockCtx`] is `(reference, resulting, frame, molecule, provenance)` — the
//! sole input to [`Strategy::run`](crate::partition::strategy::Strategy::run). It
//! carries NO input spelling and NO corpus bookkeeping (that is `bakeoff::Case`);
//! every arm is a pure function of it, which is what makes confluence provable by
//! construction.

use serde::{Deserialize, Serialize};

/// Reading-frame / axis context threaded into L2 (design §5.4). Derived from the
/// projection (reference + axis), NEVER from input spelling.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum FrameContext {
    /// Genomic, noncoding, or unknown — no coding frame available.
    NonCoding,
    /// Coding axis, with the CDS zone model (design §D2). Three window-relative
    /// facts, each derived from the projection:
    /// * `cds_phase` — the codon-grid phase (0, 1, or 2) at the window start
    ///   (reference offset 0): the number of bases before offset 0 that the
    ///   enclosing codon began. `codon_of(offset)` is `(cds_phase + offset) / 3`.
    /// * `cds_start` — window offset of the first CDS base (`c.1`), or `None` when
    ///   the CDS begins at or before the window (no 5' boundary in view). When
    ///   `Some(s)`, the codon grid is pinned by it and `cds_phase` is redundant:
    ///   the consistency invariant is `cds_phase == (3 - s % 3) % 3` (the grid
    ///   phase at offset 0 implied by a codon boundary at `s`). Use
    ///   [`FrameContext::grid_phase`], which prefers the `cds_start`-derived value.
    /// * `cds_end` — window offset one past the last CDS base (the `c.*1` seam), or
    ///   `None` when the CDS ends at or after the window (no 3' boundary in view).
    ///
    /// `cds_start`/`cds_end` carve the window into 5'UTR `[0, cds_start)`, CDS
    /// `[cds_start, cds_end)`, and 3'UTR `[cds_end, len)` — the zones C1 respects
    /// (no codon merge in a UTR or across either seam) and C2/C3 ignore (they are
    /// axis-scoped, not CDS-scoped).
    Coding {
        cds_phase: u8,
        #[serde(default)]
        cds_start: Option<usize>,
        #[serde(default)]
        cds_end: Option<usize>,
    },
}

impl FrameContext {
    /// A coding frame, asserting the codon-grid consistency invariant. When
    /// `cds_start` is in-window (`Some(s)`), the grid phase at offset 0 is fixed at
    /// `(3 - s % 3) % 3`; a `cds_phase` that disagrees is a construction bug (a
    /// generator that varied one without the other). The assert is
    /// `debug_assert!` — release corpora that hand-build literals are not slowed,
    /// but every test and dev run catches an inconsistent frame.
    pub fn coding(cds_phase: u8, cds_start: Option<usize>, cds_end: Option<usize>) -> Self {
        if let Some(s) = cds_start {
            debug_assert_eq!(
                cds_phase as usize,
                (3 - s % 3) % 3,
                "cds_phase {cds_phase} disagrees with the grid phase implied by cds_start {s}"
            );
        }
        Self::Coding {
            cds_phase,
            cds_start,
            cds_end,
        }
    }

    /// The effective codon-grid phase at the window start: the single source of
    /// truth for the codon grid. Prefers the `cds_start`-derived value when the CDS
    /// start is in-window (`(3 - s % 3) % 3`), because a pinned boundary determines
    /// the grid; falls back to the stored `cds_phase` when `cds_start` is `None`
    /// (the CDS start is off-window, so only the phase is known). `None` on a
    /// [`FrameContext::NonCoding`] frame — there is no codon grid to speak of.
    pub fn grid_phase(&self) -> Option<u8> {
        match self {
            Self::NonCoding => None,
            Self::Coding {
                cds_phase,
                cds_start,
                ..
            } => Some(
                cds_start
                    .map(|s| ((3 - s % 3) % 3) as u8)
                    .unwrap_or(*cds_phase),
            ),
        }
    }
}

/// The molecule an axis addresses. **Orthogonal to the reading frame**
/// (`FrameContext`): a molecule is DNA or RNA regardless of whether the axis is
/// coding. This is the fact L1 needs — the payload-coincidence carve-out
/// (`DNA/delins.md:44-47`) reaches every DNA axis and never the RNA axis, keyed on
/// molecule and *not* on the reading frame (`merge.rs::CoincidenceCarveOut::for_axis`
/// keys on `is_dna`, so a coding `r.` is still out of reach — jurisdiction, not
/// frame). `Protein` is reserved for the residue-core interface, wired once DNA and
/// RNA work; no bakeoff corpus emits it today.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Molecule {
    Dna,
    Rna,
    /// Reserved. No corpus produces protein cases yet; L1 treats it like a
    /// non-DNA molecule (no coincidence carve-out) until the residue model lands.
    Protein,
}

impl Molecule {
    /// The molecule an HGVS axis-letter addresses. `c/g/m/n/o` are all DNA-axis
    /// (ACGT alphabet, coincidence carve-out in reach — the #2155 all-DNA scope);
    /// `r` is RNA (out of reach → separations stand); `p` is protein (reserved).
    /// An unrecognised letter defaults to `Dna`, matching the DNA-only corpus; the
    /// loaders own the axis strings, so this is a safe fallback, not a live branch.
    pub fn from_axis(axis: &str) -> Self {
        match axis {
            "r" => Self::Rna,
            "p" => Self::Protein,
            _ => Self::Dna,
        }
    }
}

/// How the changes at this locus were OBSERVED to arrive — split-vs-merge
/// individuation (`DNA/delins.md:79-84`'s "reported (or might occur)
/// individually" discriminator; `:86-89` / `RNA/delins.md:41`'s known-polymorphism
/// reading). This is the #2155/#1419 lever: the ONLY provenance channel that may
/// affect the member partition (design §3 ①). Phase/zygosity are upstream of the
/// harness; annotation-only channels (`=` examined, predicted `p.(…)`) are out of
/// scope here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum IndividuationPolicy {
    /// ∅ — no individuation claim was made. The sequence-derived carve-outs run.
    #[default]
    Unspecified,
    /// The changes were reported/observed individually — the payload-coincidence
    /// carve-out is suppressed even on a DNA molecule.
    KeepSeparate,
    /// A member is a known frequently-occurring polymorphism — likewise suppressed.
    KnownPolymorphism,
}

impl IndividuationPolicy {
    /// Whether this policy suppresses the DNA payload-coincidence carve-out
    /// (`base_carve_out*` / `placed_gap_extension` in `coincidence.rs`). Any
    /// non-∅ individuation claim suppresses: the carve-out exists to disbelieve a
    /// separation the ALIGNMENT manufactured, and an explicit observation claim is
    /// exactly the evidence that the separation was not manufactured.
    pub fn suppresses_coincidence_collapse(self) -> bool {
        !matches!(self, Self::Unspecified)
    }
}

/// A caller-asserted examined-`=` region: a half-open span of the reference
/// window the observer states it looked at and found unchanged (`background/`
/// analysis-scope `=`; design §3 ③). Offsets are into the block's reference
/// window, the coordinate space every arm already works in — the bakeoff analog
/// of the production channel's authored `c.`-axis regions. Annotation-only: it
/// re-attaches an identity member after the sequence-first derivation, and never
/// steers the partition (contrast [`IndividuationPolicy`], the one channel that
/// does). Validated against the sequence before it is honoured — a region an edit
/// touches, or one out of range, is dropped, not re-attached.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExaminedRegion {
    /// Inclusive half-open start offset into the reference window.
    pub start: usize,
    /// Exclusive half-open end offset into the reference window.
    pub end: usize,
}

/// Epistemic provenance of the observation — a REQUIRED, first-class input to the
/// harness, ∅ by default (`Provenance::none()`). It is typed and caller-supplied:
/// no arm may ever derive it from an input spelling (confluence-by-construction,
/// design §2 invariant 1). Two channels exist: `individuation`
/// (partition-affecting — the #2155/#1419 lever) and `examined` (annotation-only
/// `=` regions, re-attached by
/// [`reattach_examined`](crate::partition::examined::reattach_examined)); the phase
/// channel is added here later, each ∅-default.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Provenance {
    pub individuation: IndividuationPolicy,
    /// Examined-`=` regions the caller asserts unchanged. ∅ (empty) by default,
    /// so a corpus that supplies none is untouched by the re-attach pass.
    #[serde(default)]
    pub examined: Vec<ExaminedRegion>,
}

impl Provenance {
    /// The ∅ provenance: no claims. `canonical(sequence, axis, ∅)` is strict
    /// sequence-confluence by construction.
    pub fn none() -> Self {
        Self::default()
    }
    /// A provenance carrying only examined-`=` regions (no individuation claim) —
    /// the annotation-only channel, for the re-attach tests and callers.
    pub fn with_examined(examined: Vec<ExaminedRegion>) -> Self {
        Self {
            examined,
            ..Self::default()
        }
    }
    /// Whether this provenance makes no claim of any kind — neither an
    /// individuation policy nor an examined region.
    pub fn is_none(&self) -> bool {
        self.individuation == IndividuationPolicy::Unspecified && self.examined.is_empty()
    }
}

/// The neutral per-block context every arm consumes: the reference/resulting byte
/// windows plus the three axis facts (frame, molecule, provenance). This is what a
/// `bakeoff::Case` carries minus the corpus bookkeeping (`id`/`origin`), and it is
/// the sole input to [`Strategy::run`](crate::partition::strategy::Strategy::run) —
/// an arm is a pure function of it. All fields borrow, so `BlockCtx` is `Copy` and
/// threads through the pipeline without allocation.
#[derive(Clone, Copy, Debug)]
pub struct BlockCtx<'a> {
    /// ACGT bytes of the reference window.
    pub reference: &'a [u8],
    /// ACGT bytes of the resulting sequence.
    pub resulting: &'a [u8],
    /// Reading-frame / axis context (design §5.4).
    pub frame: &'a FrameContext,
    /// The molecule the axis addresses (design §5.4).
    pub molecule: Molecule,
    /// Epistemic provenance (∅-default via `Provenance::none()`).
    pub provenance: &'a Provenance,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn provenance_default_is_none_and_suppresses_nothing() {
        let p = Provenance::none();
        assert_eq!(p, Provenance::default());
        assert!(p.is_none());
        assert_eq!(p.individuation, IndividuationPolicy::Unspecified);
        assert!(!p.individuation.suppresses_coincidence_collapse());
    }

    #[test]
    fn an_examined_region_is_a_claim_so_provenance_is_not_none() {
        let p = Provenance::with_examined(vec![ExaminedRegion { start: 1, end: 3 }]);
        assert!(!p.is_none(), "an examined region is a claim");
        assert_eq!(p.individuation, IndividuationPolicy::Unspecified);
        // …and it must NOT suppress the coincidence carve-out: examined is
        // annotation-only, only individuation steers the partition.
        assert!(!p.individuation.suppresses_coincidence_collapse());
    }

    #[test]
    fn both_non_empty_individuation_policies_suppress_the_carve_out() {
        assert!(IndividuationPolicy::KeepSeparate.suppresses_coincidence_collapse());
        assert!(IndividuationPolicy::KnownPolymorphism.suppresses_coincidence_collapse());
    }

    #[test]
    fn grid_phase_prefers_cds_start_derivation_when_in_window() {
        // A codon boundary at window offset s implies grid phase (3 - s%3)%3 at
        // offset 0. Check s ∈ {0,1,2,3}: phases 0, 2, 1, 0.
        for (s, want) in [(0usize, 0u8), (1, 2), (2, 1), (3, 0)] {
            let f = FrameContext::coding(want, Some(s), None);
            assert_eq!(
                f.grid_phase(),
                Some(want),
                "cds_start {s} => grid phase {want}"
            );
        }
    }

    #[test]
    fn grid_phase_falls_back_to_cds_phase_when_start_off_window() {
        for phase in 0u8..3 {
            let f = FrameContext::Coding {
                cds_phase: phase,
                cds_start: None,
                cds_end: None,
            };
            assert_eq!(f.grid_phase(), Some(phase));
        }
    }

    #[test]
    fn grid_phase_is_none_off_a_coding_frame() {
        assert_eq!(FrameContext::NonCoding.grid_phase(), None);
    }

    #[test]
    #[should_panic(expected = "disagrees with the grid phase")]
    fn coding_constructor_asserts_the_consistency_invariant() {
        // cds_start = 1 implies grid phase 2, so cds_phase = 0 is inconsistent.
        let _ = FrameContext::coding(0, Some(1), None);
    }
}
