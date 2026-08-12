//! The reference/alternate pair a description denotes.

use crate::error::FerroError;

/// A description's denoted bases, in the shape
/// [`crate::from_sequences`] consumes.
///
/// Positions here are **1-based**, matching `VcfRecord`, HGVS and
/// [`crate::from_sequences`]. [`crate::spdi::AppliedVariant`] reports the same
/// window **0-based**, and this type exists so that the two halves of the round
/// trip cannot disagree about which: an off-by-one between them would shift
/// every derived coordinate by one base and would be silent, since both
/// conventions produce a description that parses.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct SequencePair {
    /// The accession every member acts on.
    pub accession: String,
    /// 1-based position of the first base of [`Self::reference`].
    pub position: u64,
    /// The reference bases over the window.
    pub reference: String,
    /// The same window after every member has been applied.
    pub alternate: String,
    /// Whether the window is as wide 3' as it can usefully be — it ends where
    /// the pad asked, **or** where the sequence itself ends.
    ///
    /// So `true` is the ordinary answer, including for a variant near a sequence
    /// end where the pad was clipped: there was nothing further to read, so the
    /// roll is settled either way. `false` means the provider stopped short of
    /// both, and the window may have cut an ambiguous run in half — which is the
    /// one case where a caller should not trust a 3' placement against this
    /// window's edge.
    ///
    /// **This is not "the full pad was served"**, and it was called
    /// `full_pad_served` until a Python test asserted the literal reading and
    /// failed. It is also 3'-only: `Normalizer::to_sequences` prepends its 5'
    /// flank separately and that half is not reflected here. The flag is
    /// `apply_to_reference_padded`'s `window_is_final`, carried through under
    /// its own name so the two cannot drift.
    pub window_is_final: bool,
}

impl SequencePair {
    /// Build a pair from bases the caller already holds — a window out of a BAM,
    /// a VCF row, an aligner's output.
    ///
    /// This type is otherwise only ever *returned*, by
    /// [`crate::Normalizer::to_sequences`]. Without a constructor it is
    /// `#[non_exhaustive]` and so unbuildable outside this crate, which would put
    /// [`Self::derive`] out of reach of exactly the caller it is for: one who has
    /// bases and no description.
    ///
    /// `position` is **1-based** and names the first base of `reference`.
    ///
    /// [`Self::window_is_final`] is set `false`, which is the only honest answer
    /// here: a caller-supplied window carries no evidence about whether its 3'
    /// edge is where the sequence stops or where the read did.
    ///
    /// # Errors
    ///
    /// Exactly what [`crate::from_sequences`] refuses, through the same check —
    /// a zero `position`, an empty `reference`, or a symbol outside the
    /// IUPAC-IUBMB nucleotide set (`general.md:48`), which `U` is excluded from
    /// here because this surface's axis is DNA. Sharing it means a pair that
    /// constructs is a pair that derives, so the refusal arrives at the argument
    /// that caused it rather than one call later.
    ///
    /// This doc said `ACGTN` until 2026-08-11, which was the alphabet before
    /// this very change widened it — a stale claim in the one place a caller
    /// looks to find out what is accepted. The ambiguity codes (`Y`, `R`, `S`, …)
    /// are admitted; real ClinVar rows carry them.
    ///
    /// Note what is *not* checked: whether `reference` really is the reference
    /// there. That needs a provider, and the whole point of the derivation is
    /// that it has no hidden inputs. Pass bases that are not the reference and
    /// you get a faithful description of the pair you passed.
    pub fn new(
        accession: impl Into<String>,
        position: u64,
        reference: impl Into<String>,
        alternate: impl Into<String>,
    ) -> Result<Self, FerroError> {
        let (accession, reference, alternate) =
            (accession.into(), reference.into(), alternate.into());
        crate::normalize::from_sequences::validate(
            position,
            reference.as_bytes(),
            alternate.as_bytes(),
        )?;
        Ok(Self {
            accession,
            position,
            reference,
            alternate,
            window_is_final: false,
        })
    }

    /// [`crate::from_sequences_detailed`] over this window.
    ///
    /// The pair already carries all four arguments the free function takes, so
    /// spreading them back out at the call site re-passes an accession and a
    /// position this type is holding — and invites pairing one window's position
    /// with another window's bases.
    pub fn derive(
        &self,
        options: &crate::normalize::from_sequences::FromSequencesOptions,
    ) -> Result<crate::normalize::from_sequences::DerivedDescription, FerroError> {
        crate::from_sequences_detailed(
            &self.accession,
            self.position,
            &self.reference,
            &self.alternate,
            options,
        )
    }

    /// 1-based position of the last base of [`Self::reference`].
    ///
    /// Saturating rather than wrapping. [`Self::new`] refuses a zero `position`
    /// and an empty `reference`, but every field here is `pub`, so a caller can
    /// clear `reference` or zero `position` on a pair the constructor already
    /// blessed — and `0 + 0 - 1` underflows, which panics in debug and wraps to
    /// `u64::MAX` in release, since `[profile.release]` sets no
    /// `overflow-checks`. A position one before the window is a wrong answer;
    /// `u64::MAX` is a wrong answer that then propagates into a coordinate.
    #[must_use]
    pub fn end(&self) -> u64 {
        self.position
            .saturating_add(self.reference.len() as u64)
            .saturating_sub(1)
    }
}
