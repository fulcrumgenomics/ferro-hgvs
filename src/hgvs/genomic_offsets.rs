//! `background/numbering.md:6`, `:8` and `:11` — a `+`/`-` offset on a
//! genomic-family position, and where a description can still state one.
//!
//! # The rule, stated three times
//!
//! `background/numbering.md` opens by numbering the three genomic-family axes,
//! one bullet each, and ends every bullet with the same sentence:
//!
//! | axis | clause | wording |
//! |---|---|---|
//! | `g.` | `:6` | "Nucleotide numbers based on a **genomic** reference sequence **do not include** `+`, `-`, `*`, or other prefixes." |
//! | `o.` | `:8` | "Nucleotide numbers based on a **circular** reference sequence **do not include** `+`, `-`, `*` or other prefixes." |
//! | `m.` | `:11` | "Nucleotide numbers based on a **mitochondrial** reference sequence **do not include** `+`, `-`, `*`, or other prefixes." |
//!
//! So the `m.`/`o.` question is not a widening of a `g.`-only clause. The spec
//! states the prohibition per axis, three times, in three consecutive bullets,
//! and `:9` additionally makes a mitochondrial reference "a special **circular**
//! genomic reference sequence" — the two would inherit `:8` even if `:11` were
//! absent, which it is not.
//!
//! `docs/recommendations/checklist.md:16` says the same thing about `g.` in the
//! checklist's own register — "genomic (`g.`) reference sequences start with
//! nucleotide 1 and can not have nucleotides with additions like a `+`, `-`, or
//! `*`" — and `checklist.md:45` supplies the reason a *hyphen* is the shape
//! that actually turns up: it is a range written with `-` where the spec's
//! range separator is `_`, which the checklist marks `Not correct`. On a
//! coding axis `c.12-14` at least *denotes* something (nucleotide `-14` in the
//! intron 5' of `c.12`); on these three axes there is no intron to denote it
//! in, so the same spelling denotes nothing at all.
//!
//! **Why there is no offset to be measured from.** An intronic `c.100+2` is
//! anchored to an exon boundary supplied by the transcript's exon table. A
//! genomic, circular or mitochondrial accession has no exon table — it *is* the
//! sequence — so `g.266+2` names a base by an offset from nothing. That is why
//! neither half of ferro's own conversion layer could ever honour it, and why
//! both were fixed to refuse rather than to guess: see
//! `spdi::convert::reject_unresolvable_genomic_position` (#1628, #1641) and
//! `vcf::from_hgvs::reject_unresolvable_genomic_position` (#1729, #1734), which
//! cite `checklist.md:16` for the identical verdict on all three axes.
//!
//! # Enforcement stage
//!
//! Per the decided `rulings[absolute-prohibition-enforcement-stage]`, the stage
//! is **mode-dependent** at the input and **unconditional** at the output:
//!
//! - **strict** refuses at parse, because strict validates input conformance
//!   rather than merely parseability
//!   ([`apply_genomic_offset_rule`](crate::hgvs::parser::parse_hgvs_with_config));
//! - **lenient** accepts at parse with a `W4009` warning, **silent** accepts
//!   without one, and both then fail at **normalize**, because there is no
//!   base for the offset to name and so nothing to normalize
//!   (`Normalizer::normalize_core_checked`).
//!
//! There is no repair arm, and the absence is deliberate. Two repairs look
//! available and both invent a variant the caller did not describe: dropping
//! the offset answers for `g.266del`, which is a *different* nucleotide and is
//! precisely the silent flattening #1641 and #1734 were filed to stop; and
//! rewriting `g.266-268del` as the range `g.266_268del` turns a one-base
//! deletion into a three-base one on nothing more than a guess about intent.
//! Rule 1 of the README ruleset governs the output half and has no mode escape,
//! so the honest answer in the permissive modes is to fail rather than to
//! re-emit or to invent.
//!
//! # What this deliberately does not reach
//!
//! **The `c.`, `n.` and `r.` axes.** An offset is meaningful there — `:53`
//! grants intronic offsets on the non-coding axis explicitly, and the coding
//! axis is built on them — so nothing here touches them. The check is
//! **AST-keyed on [`GenomePos::offset`]**, never a scan of the rendered text: a
//! description's accession carries `-` (`NC_000023.10`), a repeat count can be
//! written `CAG[(54-68)]`, and an uncertain-boundary form spells `?`. None of
//! those reaches `GenomePos::offset`, and a textual gate would refuse all three.
//!
//! **`pter`/`qter`/`cen`.** A special position is a different defect with its
//! own history (#1643, #1662) and its own refusal; `GenomePos::Display`
//! short-circuits the offset when `special` is set, so the two cannot both be
//! rendered anyway.
//!
//! **`+?` / `-?`, the unknown-offset sentinels, are IN scope and not special.**
//! `g.100+?` renders a `+` on a genomic position exactly as `g.100+2` does, and
//! `:6` prohibits the symbol rather than a magnitude. The sentinel pair is only
//! an encoding of *which* offset is unknown; there is still no anchor for it to
//! be unknown relative to.

use crate::hgvs::interval::{GenomeInterval, UncertainBoundary};
use crate::hgvs::location::GenomePos;
use crate::hgvs::variant::{CoordinateAxis, HgvsVariant};

/// A `+`/`-` offset stated on a position of a genomic-family axis.
///
/// Carries the axis so the diagnostic can cite the clause that governs *that*
/// axis rather than a `g.`-only one, and the rendered position so the message
/// can quote what the caller actually wrote.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicAxisOffset {
    /// The axis the offending position sits on (`g.`, `m.` or `o.`).
    pub axis: CoordinateAxis,
    /// The position as it renders, e.g. `266+2` or `266-268`.
    pub stated: String,
}

impl GenomicAxisOffset {
    /// The `background/numbering.md` line prohibiting an offset on this axis.
    ///
    /// One line per axis, because the spec states the rule per axis. Citing
    /// `:6` for an `m.` description would be attributing a genomic clause to a
    /// mitochondrial one when `:11` says it directly.
    pub fn clause(&self) -> &'static str {
        match self.axis {
            CoordinateAxis::Circular => "background/numbering.md:8",
            CoordinateAxis::Mitochondrial => "background/numbering.md:11",
            // `Genomic` and — defensively — anything the walker could not have
            // produced. Only the three genomic-family axes reach this type.
            _ => "background/numbering.md:6",
        }
    }

    /// How that clause names this axis's reference sequence.
    pub fn reference_kind(&self) -> &'static str {
        match self.axis {
            CoordinateAxis::Circular => "circular",
            CoordinateAxis::Mitochondrial => "mitochondrial",
            _ => "genomic",
        }
    }

    /// The single-letter coordinate prefix, for quoting back at the caller.
    pub fn prefix(&self) -> &'static str {
        self.axis.code()
    }

    /// The refusal text, owned here so the parse-stage and normalize-stage
    /// messages cannot drift apart on which clause governs which axis.
    pub fn refusal(&self) -> String {
        format!(
            "`{}.{}` states a `+`/`-` offset on a {} position. {} says nucleotide numbers based \
             on a {} reference sequence \"do not include\" `+`, `-`, `*`, or other prefixes: a {} \
             accession carries no exon table, so there is no boundary for the offset to be \
             measured from and the position names no nucleotide. Give the position a plain \
             number, or describe the variant on a transcript accession, where an offset is \
             meaningful",
            self.prefix(),
            self.stated,
            self.reference_kind(),
            self.clause(),
            self.reference_kind(),
            self.reference_kind(),
        )
    }
}

/// Find the first `+`/`-` offset stated on a genomic-family position, if any.
///
/// Walks every position a description can carry on the `g.`, `m.` and `o.`
/// axes, including both endpoints of an interval, **both endpoints of a
/// complex `(a_b)` boundary**, every segment of a `::`-joined ring, and every
/// member of an allele.
///
/// The complex-boundary arm matters and is easy to miss:
/// [`UncertainBoundary::inner`] returns `None` for a `Range`, so a walker built
/// on `inner()` alone silently skips `g.(100+1_101-1)_(200_300)del`. Both
/// existing genomic-offset guards — SPDI's and VCF's — are written that way and
/// so have that hole; this one is not.
pub fn genomic_axis_offset(variant: &HgvsVariant) -> Option<GenomicAxisOffset> {
    match variant {
        HgvsVariant::Genome(v) => offset_in_interval(&v.loc_edit.location, CoordinateAxis::Genomic),
        HgvsVariant::Mt(v) => {
            offset_in_interval(&v.loc_edit.location, CoordinateAxis::Mitochondrial)
        }
        HgvsVariant::Circular(v) => {
            offset_in_interval(&v.loc_edit.location, CoordinateAxis::Circular)
        }
        // A ring is genomic: `GenomeRing` renders with a `g.` prefix and its
        // segments are `GenomeInterval`s on one genomic accession.
        HgvsVariant::GenomeRing(ring) => ring
            .segments
            .iter()
            .find_map(|segment| offset_in_interval(&segment.location, CoordinateAxis::Genomic)),
        HgvsVariant::Allele(allele) => allele.variants.iter().find_map(genomic_axis_offset),
        HgvsVariant::Supernumerary(inner) => genomic_axis_offset(inner),
        // Offsets are legitimate on the transcript axes and meaningless on
        // protein; neither is this clause's business.
        HgvsVariant::Cds(_)
        | HgvsVariant::Tx(_)
        | HgvsVariant::Rna(_)
        | HgvsVariant::Protein(_)
        | HgvsVariant::RnaFusion(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// Both endpoints of one interval, each of which may itself be a range.
fn offset_in_interval(
    interval: &GenomeInterval,
    axis: CoordinateAxis,
) -> Option<GenomicAxisOffset> {
    offset_in_boundary(&interval.start, axis).or_else(|| offset_in_boundary(&interval.end, axis))
}

/// Every [`GenomePos`] one boundary can hold — one for `Single`, two for a
/// complex `(a_b)` `Range`.
fn offset_in_boundary(
    boundary: &UncertainBoundary<GenomePos>,
    axis: CoordinateAxis,
) -> Option<GenomicAxisOffset> {
    match boundary {
        UncertainBoundary::Single(mu) => mu.inner().and_then(|pos| offset_at(pos, axis)),
        UncertainBoundary::Range { start, end } => start
            .inner()
            .and_then(|pos| offset_at(pos, axis))
            .or_else(|| end.inner().and_then(|pos| offset_at(pos, axis))),
    }
}

/// One position.
///
/// `offset: Some(0)` is not matched: [`GenomePos`]'s `Display` prints nothing
/// for it, so such a position renders as a bare base and states no prohibited
/// symbol. Refusing it would refuse a description that is textually conformant.
fn offset_at(pos: &GenomePos, axis: CoordinateAxis) -> Option<GenomicAxisOffset> {
    match pos.offset {
        Some(0) | None => None,
        Some(_) => Some(GenomicAxisOffset {
            axis,
            stated: pos.to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;

    fn found(input: &str) -> Option<GenomicAxisOffset> {
        genomic_axis_offset(&parse_hgvs(input).expect("input parses"))
    }

    #[test]
    fn a_plus_offset_on_the_genomic_axis_is_found() {
        let hit = found("NC_TEST.1:g.266+2del").expect("`+2` is stated");
        assert_eq!(hit.axis, CoordinateAxis::Genomic);
        assert_eq!(hit.stated, "266+2");
        assert_eq!(hit.clause(), "background/numbering.md:6");
    }

    #[test]
    fn a_hyphen_range_on_the_genomic_axis_is_found_as_an_offset() {
        // `checklist.md:45`'s shape: a range written with `-`. The grammar
        // reads it as base 266 offset -268, which is why it is this check's
        // business rather than a separate range rule.
        let hit = found("NC_TEST.1:g.266-268del").expect("`-268` is stated");
        assert_eq!(hit.stated, "266-268");
    }

    #[test]
    fn each_genomic_family_axis_cites_its_own_clause() {
        for (input, axis, clause) in [
            (
                "NC_000001.11:g.100+2del",
                CoordinateAxis::Genomic,
                "background/numbering.md:6",
            ),
            (
                "NC_012920.1:m.100+2del",
                CoordinateAxis::Mitochondrial,
                "background/numbering.md:11",
            ),
            (
                "NC_001416.1:o.100+2del",
                CoordinateAxis::Circular,
                "background/numbering.md:8",
            ),
        ] {
            let hit = found(input).unwrap_or_else(|| panic!("{input} states an offset"));
            assert_eq!(hit.axis, axis, "{input}");
            assert_eq!(hit.clause(), clause, "{input}");
            assert!(
                hit.refusal().contains(clause),
                "the refusal must cite the axis's own clause: {}",
                hit.refusal()
            );
        }
    }

    #[test]
    fn an_unknown_offset_sentinel_is_still_an_offset() {
        // `+?`/`-?` render a prohibited symbol just as `+2` does; `:6`
        // prohibits the symbol, not a magnitude.
        assert_eq!(
            found("NC_000001.11:g.100+?del").expect("`+?`").stated,
            "100+?"
        );
        assert_eq!(
            found("NC_000001.11:g.100-?del").expect("`-?`").stated,
            "100-?"
        );
    }

    #[test]
    fn a_plain_genomic_position_is_not_matched() {
        assert!(found("NC_TEST.1:g.266del").is_none());
        assert!(found("NC_000001.11:g.100_200del").is_none());
        // The accession itself carries no offset even though it carries digits
        // and a `.`; a textual gate is what this must not become.
        assert!(found("NC_000023.10:g.100A>G").is_none());
    }

    #[test]
    fn a_transcript_axis_offset_is_left_alone() {
        // The whole point of the axis scoping: `:53` grants these.
        assert!(found("NM_000088.3:c.100+5A>G").is_none());
        assert!(found("NR_003051.3:n.100+5A>G").is_none());
    }

    #[test]
    fn an_offending_member_is_found_inside_an_allele() {
        let hit = found("NC_TEST.1:g.[265del;266+2del]").expect("the second member states one");
        assert_eq!(hit.stated, "266+2");
    }
}
