//! `background/numbering.md:50`–`:54`'s numbering zones for the `n.` axis, and
//! the two markers a description on that axis may therefore not state.
//!
//! # The rule
//!
//! `numbering.md:50`–`:54` enumerates the non-coding DNA axis in full:
//!
//! ```text
//! 52: - nucleotide numbering is `n.1`, `n.2`, `n.3`, ..., etc., from the first
//!      to the last nucleotide of the reference sequence.
//! 53: - nucleotides in introns are numbered as for coding DNA reference
//!      sequences (see above), although proceeded by `n.` (not `c.`).
//! 54: - it is **not** allowed to describe variants in nucleotides beyond the
//!      boundaries of a transcript reference sequence, using that transcript
//!      reference sequence.
//! ```
//!
//! `:52` runs from the first nucleotide to the last and names no `*` zone and no
//! negative one. Read alone that could be terseness rather than exclusion. `:53`
//! is what settles it: the spec knew how to add a zone to this axis and added
//! exactly one, so the omission of the other two is a choice. `:54` then
//! forbids naming a nucleotide past the boundary outright, and `numbering.md:45`
//! records that the sibling proposal — extending the numbering "to specifically
//! mark non-transcribed nucleotides" — "[has] been made but [was] rejected".
//!
//! # The counter-reading, and why `:53` does not carry it
//!
//! `:53` incorporates the coding numbering by reference, and the coding
//! numbering is where `-N` and `*N` are defined — so it can look as though `:53`
//! pulls those two zones onto `n.` as well. The structure of the section it
//! points at rules that out. `:18`–`:44` is three named bullets: **protein
//! coding region** (`:20`), **untranslated region (UTR)** (`:28`, which is where
//! `c.-1`/`c.*1` are defined) and **introns** (`:33`). `:53`'s subject is
//! "nucleotides **in introns**", so its "see above" is scoped by its own
//! sentence to the third of those; the `-`/`*` zones sit in the second, which it
//! never reaches.
//!
//! The coding section's closing bullet confirms the reading, being `:54`'s
//! sentence verbatim: `:43`–`:44` **transcript flanking** — "it is **not**
//! allowed to describe variants in nucleotides beyond the boundaries of a
//! transcript reference sequence". `c.*1` and that prohibition coexist because
//! the UTR is *inside* the transcript. `n.` has no UTR bullet to put anything
//! inside, so `:54` stands unopposed.
//!
//! # Why the coding axis keeps both, and this is not an inconsistency
//!
//! On `c.` the two zones are anchored to the **CDS**, not to the sequence:
//! `c.-1` is the base before the ATG and `c.*1` the base after the stop codon,
//! and both are still *inside* the transcript, which is why `:54` does not
//! reach them. `n.` has no CDS to anchor them to, so the same spelling can only
//! name something outside the reference sequence. Nothing here touches `c.` —
//! including ferro's deliberate past-the-end `c.*` in the #797 poly-A carve-out.
//!
//! # Why `r.` is out of scope, on the spec's own words
//!
//! `numbering.md:58`: "nucleotide numbering for an RNA reference sequence
//! follows that of the associated **coding or non-coding** DNA reference
//! sequence; nucleotide `r.123` relates to `c.123` **or** `n.123`." `:61` gives
//! a coding RNA reference the coding numbering — so `r.*5` is conformant there —
//! and `:60` gives a non-coding one this axis's numbering, where it is not.
//!
//! The zone set is therefore a property of the **reference**, and the parser
//! holds no provider, so on `r.` it cannot tell the conformant case from the
//! prohibited one. Refusing there would refuse conformant descriptions.
//!
//! `n.` is decidable at parse for a reason that does not transfer:
//! `general.md:26` makes the prefix a claim about the reference *type* — "`n`
//! for a non-coding DNA reference sequence" — so `:52` binds on the prefix
//! alone, with nothing to look up.
//!
//! # THE TWO ZONES ARE ENFORCED AT DIFFERENT STAGES, AND THE SPLIT IS MEASURED
//!
//! This is the part to read before changing anything here. One clause reading
//! covers both markers; the **stages differ**, and they differ on evidence
//! rather than on argument.
//!
//! | marker | stage | code |
//! |---|---|---|
//! | `n.*N` | refused at **parse, in every mode** — including the bare [`parse_hgvs`](crate::parse_hgvs) entry | `E1003 InvalidPosition` |
//! | `n.-N` | refused at **strict** parse only; lenient and silent accept | `W4008` |
//!
//! **`n.*N` costs nothing to refuse.** Measured over ferro's four committed
//! corpora — `clinvar_hgvs_500k`, `clinvar_hgvs_unique`, `cmrg_genes_exhaustive`
//! and `paraphase_genes_exhaustive` — **0 of 103,762 `n.`-axis rows** state a
//! `*` marker. That zero is a real one and not a structural artifact: the same
//! scan over the same rows finds the `-` marker **5** times, so the counter can
//! demonstrably see this shape when it is there.
//!
//! **`n.-N` costs five real ClinVar rows**, which is why it is not refused
//! unconditionally: `NR_003051.3:n.-57T>C`, `NR_003051.3:n.-30_-7dup` and
//! `LRG_163t1:n.-5delins17` (all RMRP, whose upstream promoter variants are the
//! clinically conventional case for this spelling), `NR_029595.1:n.-4771G>T`
//! (MIR208A) and `NR_033294.1:n.-6G>A` (SNORD118). These are descriptions NCBI
//! publishes today. Lenient and silent therefore keep accepting them, and so
//! does the bare `parse_hgvs` entry every `ferro` subcommand and all four corpus
//! consumers call.
//!
//! So the asymmetry is not a hedge about the clause. Both spellings are read the
//! same way; the deviation from the reading is **tolerated in the permissive
//! modes exactly where real clinical data relies on it**, and nowhere else.
//!
//! # The `n.*N` stage departs from a decided ruling. Say so; do not restate it.
//!
//! `rulings[absolute-prohibition-enforcement-stage]` (decided 2026-08-10) rules
//! that enforcement is **mode-dependent, uniformly**: strict fails at parse,
//! lenient and silent do not validate input conformance. Refusing `n.*N` in
//! **every** mode does not fit that schedule, and it is worth being exact about
//! why, because two tempting reconciliations are weaker than they look:
//!
//! - *"`n.*5` is a grammar matter, not a conformance one — it names no position
//!   at all, like `n.0`."* Arguable, and **not conclusive**. `checklist.md:16`
//!   says a genomic reference "can not have nucleotides with additions like a
//!   `+`, `-`, or `*`", making `g.*10` exactly parallel — a marker on an axis
//!   that has no such zone — and that record's own census places `g.*10` under
//!   mode-gating and states that lenient should newly *accept* it. A scoping
//!   argument that would also carve out `g.*10` proves too much.
//! - *"The record only governs `checklist.md` items."* True of its question
//!   sentence, but its rationale generalises ("uniformly across every absolute
//!   prohibition"), and `:54` is an absolute prohibition in the same register.
//!
//! What does hold is narrower and empirical. The record rejected unconditional
//! refusal for **one stated reason** — it "would newly refuse inputs ferro
//! accepts today, with no escape for a caller round-tripping a real-world
//! corpus". For `n.*N` that objection is measurably empty: there is no such
//! caller in 103,762 rows. So this is a **maintainer's decision** under rule 6
//! of the `README.md` ruleset, disclosed under rule 7, and it is **revisitable
//! on user demand** — if someone reports a real `n.*N` corpus, the honest
//! response is to move it onto the `n.-N` schedule, not to defend the refusal.
//!
//! The carve-out is **recorded in the ledger**, as the `AMENDMENT, 2026-08-12`
//! section of `rulings[absolute-prohibition-enforcement-stage]`'s rationale, so
//! the departure and the record it departs from sit in one place.
//! `issue_1748_noncoding_axis_zones` pins its *shape* —
//! `the_unconditional_arm_is_a_disclosed_departure` — so neither arm can be
//! quietly moved onto the other's schedule without the amendment moving too.
//!
//! # What it is keyed on
//!
//! [`TxPos`] carries the two zones as `downstream: bool` (the `*`) and a
//! negative `base`, so the predicate reads them off the AST rather than
//! scanning the rendered text — a description's accession may contain `-`, and
//! an intronic offset is spelled with one too.
//!
//! **`base == 0` is deliberately not reported here.** `n.0` is refused at the
//! grammar in every mode and already surfaces as `E1003 InvalidPosition`;
//! classifying it here as well would give one input two diagnostics.
//!
//! **An offset is never a marker either**, whatever its sign: `n.6-3` is
//! `:53`'s explicitly granted intronic numbering, and `-3` there is a distance
//! from position 6 rather than a zone. The predicate therefore tests `base`,
//! never `offset`.

use crate::hgvs::interval::{Interval, TxInterval};
use crate::hgvs::location::TxPos;
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::HgvsVariant;

/// Which of the two zones `numbering.md:52` does not define a description
/// states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NonCodingZone {
    /// A `*N` marker — a nucleotide after the last one `:52` numbers.
    Downstream,
    /// A `-N` marker — a nucleotide before the first one `:52` numbers.
    Upstream,
}

impl NonCodingZone {
    /// The marker as it is spelled, for a diagnostic.
    #[must_use]
    pub fn marker(self) -> &'static str {
        match self {
            Self::Downstream => "*",
            Self::Upstream => "-",
        }
    }

    /// Where the prohibited nucleotide sits relative to the reference sequence.
    #[must_use]
    pub fn relation(self) -> &'static str {
        match self {
            Self::Downstream => "after the last nucleotide of",
            Self::Upstream => "before the first nucleotide of",
        }
    }
}

/// One out-of-zone position found on an `n.` axis, with the text that states it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NonCodingZoneMarker {
    /// Which zone the position claims.
    pub zone: NonCodingZone,
    /// The position as spelled, e.g. `*5` or `-100`.
    pub stated: String,
    /// Whether the transcript being numbered is also the **reference sequence**.
    ///
    /// `numbering.md:54` prohibits describing a variant in a nucleotide beyond a
    /// transcript's boundaries **"using that transcript reference sequence"** —
    /// a conditional whose subject is the reference. On a bare transcript
    /// (`NR_003051.3:n.-57T>C`, `LRG_163t1:n.-5delins17`) the condition holds.
    /// In the selector form `NG_007485.1(NR_003529.3):n.-40000del` it does
    /// **not**: the reference sequence is the genomic `NG_007485.1`, which does
    /// contain the nucleotide, and the transcript is only the numbering
    /// selector. The refusal there stands on `:52` alone — the axis has no such
    /// zone to number from — and the message must not assert `:54`'s clause for
    /// a shape it does not reach.
    pub numbers_a_bare_transcript: bool,
}

impl NonCodingZoneMarker {
    /// `:54`'s sentence, included only where that clause's own condition holds.
    ///
    /// See [`Self::numbers_a_bare_transcript`]. Returns a trailing space so it
    /// can be spliced into a message or omitted without leaving a double space.
    fn boundary_clause(&self) -> &'static str {
        if self.numbers_a_bare_transcript {
            ":54 separately forbids describing a variant in a nucleotide beyond that \
             transcript's boundaries, using that transcript reference sequence. "
        } else {
            ""
        }
    }

    /// The refusal message, less any error-code prefix.
    ///
    /// Shared by the unconditional `*` refusal in
    /// [`crate::hgvs::parser::variant::parse_variant`] and the strict-mode `-`
    /// refusal in [`crate::hgvs::parser::parse_hgvs_with_config`], so the `:54`
    /// scoping above cannot be right in one message and wrong in the other.
    #[must_use]
    pub fn refusal(&self) -> String {
        format!(
            "`n.{stated}` uses the `{marker}` marker, which background/numbering.md:52 does not \
             put on the non-coding DNA axis: that axis is numbered `n.1`, `n.2`, `n.3`, ..., from \
             the first to the last nucleotide of the reference sequence, and :53 grants intronic \
             offsets as its only other zone. There is no `{marker}` zone here for `{stated}` to \
             be numbered from, so it names a nucleotide {relation} the reference sequence. \
             {boundary}The `-`/`*` zones exist on the CODING axis (`c.`), where they are anchored \
             to the CDS and are still inside the transcript; name the nucleotide on a reference \
             that contains it instead.",
            stated = self.stated,
            marker = self.zone.marker(),
            relation = self.zone.relation(),
            boundary = self.boundary_clause(),
        )
    }

    /// The lenient-mode warning message. Same clause scoping as [`Self::refusal`].
    #[must_use]
    pub fn warning(&self) -> String {
        format!(
            "`n.{stated}` names a nucleotide {relation} the reference sequence; \
             background/numbering.md:52 numbers the non-coding axis from the first nucleotide to \
             the last, and :53 grants intronic offsets as its only other zone. {boundary}Accepted \
             as authored: this mode does not validate input conformance.",
            stated = self.stated,
            relation = self.zone.relation(),
            boundary = self.boundary_clause(),
        )
    }
}

/// Classify one `n.`-axis position, if it states `want`.
fn zone_of(pos: &TxPos, want: NonCodingZone, bare: bool) -> Option<NonCodingZoneMarker> {
    // `base == 0` is refused at the grammar and reported as `E1003
    // InvalidPosition`; an offset of either sign is `:53`'s intronic numbering
    // and legal. Neither is a zone.
    let states = match want {
        NonCodingZone::Downstream => pos.downstream,
        // A `*N` position also carries a positive `base`, so the upstream test
        // must exclude it rather than rely on the sign alone.
        NonCodingZone::Upstream => !pos.downstream && pos.base < 0,
    };
    states.then(|| NonCodingZoneMarker {
        zone: want,
        stated: pos.to_string(),
        numbers_a_bare_transcript: bare,
    })
}

/// The first endpoint of `interval` that states `want`.
///
/// Both endpoints are read, and each boundary's full shape is walked — a
/// `Range` boundary (`(100_200)_300`) holds two positions of its own, and an
/// uncertain position (`(100)`) is still a position. `Mu::Unknown` (`?`) states
/// no coordinate at all and so cannot state a zone.
fn zone_in_interval(
    interval: &TxInterval,
    want: NonCodingZone,
    bare: bool,
) -> Option<NonCodingZoneMarker> {
    let in_mu = |mu: &Mu<TxPos>| mu.inner().and_then(|p| zone_of(p, want, bare));
    let in_boundary = |boundary: &crate::hgvs::interval::UncertainBoundary<TxPos>| match boundary {
        crate::hgvs::interval::UncertainBoundary::Single(mu) => in_mu(mu),
        crate::hgvs::interval::UncertainBoundary::Range { start, end } => {
            in_mu(start).or_else(|| in_mu(end))
        }
    };
    let Interval { start, end } = interval;
    in_boundary(start).or_else(|| in_boundary(end))
}

/// The first `numbering.md:52`-undefined position of zone `want` that `variant`
/// states on an `n.` axis, if any.
///
/// Walks allele members and the inner description of a `sup` marker, so a
/// prohibited position cannot be laundered behind a composite spelling — a cis
/// allele (`n.[*5A>G;100del]`) and a trans pair (`n.[*5A>G];[100del]`) are both
/// reached. A **ring** is deliberately *not* walked: its segments are genomic,
/// and `GenomeRing` yields `None` for the reason the arm list below gives.
///
/// The zone is a **parameter** rather than a scan for either, because the two
/// are enforced at different stages: `Downstream` is refused unconditionally at
/// parse and `Upstream` only in strict mode. See the module docs for the
/// measurement that split them.
///
/// The match is exhaustive rather than wildcarded so a new axis has to be
/// classified here rather than defaulted to "states no marker". Every arm other
/// than [`HgvsVariant::Tx`] yields `None`, and the reasons differ:
///
/// - `Cds` — `c.` genuinely has both zones, anchored to the CDS
///   (`numbering.md:20`–`:34`). Refusing them would be a serious regression.
/// - `Rna` — `numbering.md:58` makes the zone set a property of the underlying
///   reference, which the parser cannot resolve. See the module docs.
/// - `Genome`, `Mt`, `Circular`, `GenomeRing` — a genomic reference is not a
///   transcript, so neither `:52` nor `:54` addresses it. `checklist.md:16`
///   prohibits `g.*10`, and that is a *different* clause with its own row in
///   `rulings[absolute-prohibition-enforcement-stage]`'s census (#1628); folding
///   it in here would attribute a genomic refusal to a non-coding clause.
/// - `Protein` — no nucleotide numbering at all.
/// - `RnaFusion` — its breakpoints hold an accession, an optional gene symbol
///   and an interval on the `r.` axis, which is the `Rna` case above.
/// - `NullAllele` / `UnknownAllele` — whole-allele markers with no position.
#[must_use]
pub fn noncoding_zone_marker(
    variant: &HgvsVariant,
    want: NonCodingZone,
) -> Option<NonCodingZoneMarker> {
    match variant {
        HgvsVariant::Tx(v) => zone_in_interval(
            &v.loc_edit.location,
            want,
            v.accession.genomic_context.is_none(),
        ),
        HgvsVariant::Allele(allele) => allele
            .variants
            .iter()
            .find_map(|inner| noncoding_zone_marker(inner, want)),
        HgvsVariant::Supernumerary(inner) => noncoding_zone_marker(inner, want),
        HgvsVariant::Genome(_)
        | HgvsVariant::Cds(_)
        | HgvsVariant::Rna(_)
        | HgvsVariant::Mt(_)
        | HgvsVariant::Circular(_)
        | HgvsVariant::GenomeRing(_)
        | HgvsVariant::Protein(_)
        | HgvsVariant::RnaFusion(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;

    /// `n.*N` no longer parses, so a `Downstream` fixture has to be built rather
    /// than parsed. `TxPos::downstream` is public API and directly
    /// constructible, which is exactly why the walk still has to classify it.
    fn tx_with(accession: &str, start: TxPos, end: TxPos) -> HgvsVariant {
        let base = parse_hgvs(&format!("{accession}:n.10_20del")).expect("template parses");
        let HgvsVariant::Tx(mut v) = base else {
            unreachable!("`n.` parses as Tx")
        };
        v.loc_edit.location = TxInterval::new(start, end);
        HgvsVariant::Tx(v)
    }

    fn upstream_marker_of(input: &str) -> Option<NonCodingZoneMarker> {
        noncoding_zone_marker(
            &parse_hgvs(input).expect("input parses"),
            NonCodingZone::Upstream,
        )
    }

    #[test]
    fn a_downstream_marker_is_found() {
        let v = tx_with("NR_037639.1", TxPos::downstream(5), TxPos::downstream(5));
        let found = noncoding_zone_marker(&v, NonCodingZone::Downstream).expect("`*5` is stated");
        assert_eq!(found.zone, NonCodingZone::Downstream);
        assert_eq!(found.stated, "*5");
    }

    /// The two zones are separate queries, so neither may answer for the other.
    #[test]
    fn the_two_zones_do_not_answer_for_each_other() {
        let star = tx_with("NR_037639.1", TxPos::downstream(5), TxPos::downstream(5));
        assert!(noncoding_zone_marker(&star, NonCodingZone::Upstream).is_none());
        assert!(upstream_marker_of("NR_037639.1:n.-100_-50del").is_some());
        assert!(noncoding_zone_marker(
            &parse_hgvs("NR_037639.1:n.-100_-50del").unwrap(),
            NonCodingZone::Downstream
        )
        .is_none());
    }

    #[test]
    fn an_upstream_marker_is_found() {
        let found = upstream_marker_of("NR_037639.1:n.-100_-50del").expect("`-100` is stated");
        assert_eq!(found.zone, NonCodingZone::Upstream);
        assert_eq!(found.stated, "-100");
    }

    /// One endpoint outside is enough, from either end — these are #255's
    /// boundary-spanning shapes.
    #[test]
    fn either_endpoint_alone_is_enough() {
        assert_eq!(
            upstream_marker_of("NR_037639.1:n.-3_5del")
                .expect("start is outside")
                .zone,
            NonCodingZone::Upstream
        );
        let end_outside = tx_with("NR_037639.1", TxPos::new(40), TxPos::downstream(3));
        assert_eq!(
            noncoding_zone_marker(&end_outside, NonCodingZone::Downstream)
                .expect("end is outside")
                .zone,
            NonCodingZone::Downstream
        );
    }

    /// `:53` grants introns on this axis explicitly. A negative *offset* is a
    /// distance from the base, not a zone, and must not be confused for one.
    #[test]
    fn intronic_offsets_are_not_markers() {
        assert!(upstream_marker_of("NG_012337.1(NR_037639.1):n.5+3A>G").is_none());
        assert!(upstream_marker_of("NG_012337.1(NR_037639.1):n.6-3A>G").is_none());
        assert!(upstream_marker_of("NG_012337.1(NR_037639.1):n.100+10del").is_none());
    }

    #[test]
    fn in_transcript_positions_are_not_markers() {
        assert!(upstream_marker_of("NR_024540.1:n.1A>G").is_none());
        assert!(upstream_marker_of("NR_024540.1:n.1786A>G").is_none());
        assert!(upstream_marker_of("NR_037639.1:n.100_200del").is_none());
    }

    /// The coding axis has both zones for real; nothing here may reach them.
    #[test]
    fn the_coding_axis_states_no_marker() {
        for input in ["NM_000492.4:c.*1A>G", "NM_000492.4:c.-1A>G"] {
            let v = parse_hgvs(input).expect("parses");
            assert!(noncoding_zone_marker(&v, NonCodingZone::Upstream).is_none());
            assert!(noncoding_zone_marker(&v, NonCodingZone::Downstream).is_none());
        }
    }

    /// `numbering.md:58` — an `r.` description's zone set depends on a reference
    /// the parser does not hold, so it is left alone in both directions.
    #[test]
    fn the_rna_axis_states_no_marker() {
        for input in ["NM_003002.4:r.*5a>g", "NM_003002.4:r.-5a>g"] {
            let v = parse_hgvs(input).expect("parses");
            assert!(noncoding_zone_marker(&v, NonCodingZone::Upstream).is_none());
            assert!(noncoding_zone_marker(&v, NonCodingZone::Downstream).is_none());
        }
    }

    #[test]
    fn an_allele_member_is_reached() {
        assert_eq!(
            upstream_marker_of("NR_037639.1:n.[100del;-5A>G]")
                .expect("the second member states it")
                .zone,
            NonCodingZone::Upstream
        );
    }

    /// `:54`'s clause is conditioned on the transcript being the reference
    /// sequence. In the selector form it is not, and the message must not claim
    /// it is — that is the scoping this flag exists for.
    #[test]
    fn the_boundary_clause_is_scoped_to_a_bare_transcript() {
        let bare = upstream_marker_of("NR_003529.3:n.-40000del").expect("bare transcript");
        assert!(bare.numbers_a_bare_transcript);
        assert!(bare.refusal().contains(":54"));
        assert!(bare.warning().contains(":54"));

        let selector =
            upstream_marker_of("NG_007485.1(NR_003529.3):n.-40000del").expect("selector form");
        assert!(!selector.numbers_a_bare_transcript);
        assert!(
            !selector.refusal().contains(":54"),
            "the genomic reference DOES contain the nucleotide, so :54 is not \
             the operative clause: {}",
            selector.refusal()
        );
        assert!(!selector.warning().contains(":54"));
        // Both still stand on `:52`, which is what actually refuses them.
        assert!(selector.refusal().contains("numbering.md:52"));
    }
}
