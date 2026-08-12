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
    /// [`Self::trim_to`] and [`crate::Normalizer::reanchor`] out of reach of
    /// exactly the caller they are for: one who has bases and no description.
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
    /// position this type is holding — and, after a [`Self::trim_to`], invites
    /// pairing the *original* position with the *trimmed* bases.
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

    /// Narrow this window to `[start, end]`, **trimming matching bases only**.
    ///
    /// The reference-free half of re-anchoring. Use it to hold a derivation
    /// inside a region it must not leave — a target region, an amplicon, a
    /// tiling window — when every caller's raw window is at least as wide as
    /// that region.
    ///
    /// `None` leaves that edge where it is. Both bounds are 1-based inclusive,
    /// matching [`Self::position`].
    ///
    /// # This can only narrow, and narrowing costs the 3' rule
    ///
    /// Widening needs bases this type does not hold, so it lives on
    /// [`crate::Normalizer::reanchor`], which has a reference to read them from.
    ///
    /// And narrowing is not free: `from_sequences` may not shift outside the
    /// window it is given, so a bound that cuts an ambiguous run pulls the
    /// placement back to the bound. That is the *point* when the bound is a
    /// region the variant must not leave, and a **footgun** when it is an
    /// arbitrary window — every caller anchored to it will agree with each other
    /// and disagree with the reference-anchored answer. The derivation reports
    /// this as `placement_bounded_by_window`; it is not silent.
    ///
    /// # Errors
    ///
    /// Refuses rather than clamping, in every case — a silent clamp would hand
    /// back a window the caller did not ask for and cannot detect:
    ///
    /// * a bound that would *widen* the window (there are no bases to widen
    ///   with — the message names `reanchor`);
    /// * a bound that would cut into a base the two sequences disagree on,
    ///   naming the coordinate where the disagreement starts. **Case is not a
    ///   disagreement** — a soft-masked reference against an upper-case
    ///   alternate is ordinary input, and the comparison folds;
    /// * a bound that would leave the reference empty, or `start` past `end`.
    ///
    /// Case is otherwise preserved: this method fetches nothing, so it cannot
    /// manufacture a mixed-case pair and has no reason to rewrite the caller's
    /// bases. [`crate::Normalizer::reanchor`], which does fetch, folds instead.
    pub fn trim_to(&self, start: Option<u64>, end: Option<u64>) -> Result<Self, FerroError> {
        let (current_start, current_end) = (self.position, self.end());
        let start = start.unwrap_or(current_start);
        let end = end.unwrap_or(current_end);

        let invalid = |msg: String| FerroError::InvalidCoordinates { msg };
        if start > end {
            return Err(invalid(format!(
                "trim_to({start}, {end}) on {}: start is past end",
                self.accession
            )));
        }
        for (label, wanted, widens) in [
            ("start", start, start < current_start),
            ("end", end, end > current_end),
        ] {
            if widens {
                return Err(invalid(format!(
                    "trim_to would widen the window: {label} {wanted} is outside the window \
                     [{current_start}, {current_end}] on {}. `trim_to` only removes bases; use \
                     `Normalizer::reanchor`, which holds a reference and can supply them",
                    self.accession
                )));
            }
        }

        let head = (start - current_start) as usize;
        let tail = (current_end - end) as usize;
        // The two sequences differ in length, so a 3' trim is counted from each
        // one's own end rather than from a shared index.
        let (ref_len, alt_len) = (self.reference.len(), self.alternate.len());
        if head + tail >= ref_len || head + tail >= alt_len {
            return Err(invalid(format!(
                "trim_to({start}, {end}) on {} would leave nothing of the {ref_len}-base \
                 reference or the {alt_len}-base alternate",
                self.accession
            )));
        }

        // "Matching bases only": a trimmed base must be one both sequences
        // agree on, or the trim would silently change what the pair denotes.
        //
        // Compared **case-insensitively**, because case is not a difference
        // anywhere else on this surface: `Base::from_char` folds, `validate`
        // says so in as many words, and `Normalizer::to_sequences` upper-cases
        // its whole window. A byte comparison here made a soft-masked reference
        // against an upper-case alternate — an ordinary pileup — refuse a legal
        // trim with "the reference and alternate first differ at position N",
        // naming a coordinate where they do not differ.
        let (r, a) = (self.reference.as_bytes(), self.alternate.as_bytes());
        let differs = |x: u8, y: u8| !x.eq_ignore_ascii_case(&y);
        if !r[..head].eq_ignore_ascii_case(&a[..head]) {
            let at = (0..head).find(|i| differs(r[*i], a[*i])).unwrap_or(0);
            return Err(invalid(format!(
                "trim_to cannot trim to start {start} on {}: the reference and alternate first \
                 differ at position {}, so trimming there would change what the pair denotes",
                self.accession,
                current_start + at as u64
            )));
        }
        if !r[ref_len - tail..].eq_ignore_ascii_case(&a[alt_len - tail..]) {
            let at = (0..tail)
                .find(|i| differs(r[ref_len - 1 - *i], a[alt_len - 1 - *i]))
                .unwrap_or(0);
            return Err(invalid(format!(
                "trim_to cannot trim to end {end} on {}: the reference and alternate differ at \
                 position {}, so trimming there would change what the pair denotes",
                self.accession,
                current_end - at as u64
            )));
        }

        Ok(Self {
            accession: self.accession.clone(),
            position: start,
            reference: self.reference[head..ref_len - tail].to_string(),
            alternate: self.alternate[head..alt_len - tail].to_string(),
            // Trimming the 3' edge moves the window in off wherever it used to
            // stop, so whatever settled the old edge no longer settles this one.
            window_is_final: self.window_is_final && tail == 0,
        })
    }
}
