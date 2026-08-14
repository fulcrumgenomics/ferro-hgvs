//! `docs/recommendations/checklist.md:20` — an intronic position stated on a
//! **bare** transcript reference, with no genomic reference to annotate it on.
//!
//! # The clause
//!
//! > an `NM_` reference sequence can only be used to describe variants in
//! > introns using a `c.` prefix when a genomic reference sequence is given on
//! > which the coding DNA reference sequence is annotated
//!
//! The reason is definitional rather than stylistic: a coding or non-coding DNA
//! reference sequence *is* the spliced transcript, so it contains no introns for
//! an offset to name. `NM_000088.3:c.589-1G>T` states a position one nucleotide
//! 5' of an exon boundary on a sequence that has no such nucleotide. The
//! wrapper forms the clause asks for — `NG_(NM_)`, `NC_(NM_)` — supply the
//! genomic sequence the offset is measured on, which is why they are accepted
//! here and the bare form is not.
//!
//! # This is a CONDITIONAL clause, and the mode split follows from that
//!
//! `rulings[bare-transcript-intronic-position]` decided this, and nothing in
//! this module reopens it. Two things make the clause conditional in form
//! rather than an absolute prohibition on a shape:
//!
//! - It says "can only be used … **when**", which states a condition on the
//!   reference sequence rather than banning a spelling.
//! - The spec reads a bare-`c.` intronic description itself, four items later:
//!   `checklist.md:45` glosses `c.12-14del` as "a deletion of nucleotide -14 in
//!   the intron directly 5' of nucleotide `c.12`", with no genomic wrapper in
//!   sight. A clause the spec does not apply to its own worked gloss cannot be
//!   read as an absolute bar.
//!
//! So strict input hygiene refuses the bare form and lenient accepts it. That
//! is the ruling; the corpus classifies these rows `Strength::Conditional` for
//! the same reason.
//!
//! # Enforcement stage
//!
//! Per the decided `rulings[absolute-prohibition-enforcement-stage]`, the stage
//! is **mode-dependent** and the check belongs at **parse**:
//!
//! - **strict** refuses at parse, because strict validates input conformance
//!   rather than merely parseability
//!   ([`apply_bare_transcript_intron_rule`](crate::hgvs::parser::parse_hgvs_with_config));
//! - **lenient** accepts at parse with a `W4007` warning, and **silent**
//!   accepts without one.
//!
//! That record's census named this clause's stage as the third of the three
//! things it had to fix — "`checklist.md:20`'s refusal is at normalize, not at
//! parse … moving the stage is cosmetic there, but it should move for
//! uniformity" — and #1630 is where it moved. The ground it gives is the
//! general one: "whether the INPUT conforms is answered before the input is
//! accepted, not part-way through normalizing it."
//!
//! **The normalize-stage rung is kept, and is not redundant.** It answers for a
//! caller who reaches the normalizer through some *other* door — the
//! config-less `parse_hgvs`, or a lenient parse followed by a strict
//! `Normalizer` — and it is the rung that carries the `EINTRONIC` tag the
//! Mutalyzer conformance map keys off. The two rungs answer for two different
//! callers and neither subsumes the other.
//!
//! # There is no repair arm
//!
//! Re-expressing a bare `NM_…:c.20+2del` on a genomic reference needs a
//! genomic accession and an exon table, neither of which the parser holds — and
//! choosing *which* genomic parent is a question with more than one answer for
//! any transcript with several placements. `rulings[bare-transcript-intronic-
//! position]` also settles the direction: "An input that already names one is
//! still left as authored", so a repair here would re-spell a description the
//! ruling says to leave alone. Ferro does re-parent an intronic offset it
//! *manufactured* itself (#1704), which is a different question decided on
//! provenance the normalizer has and the parser does not.

use crate::hgvs::interval::UncertainBoundary;
use crate::hgvs::location::{CdsPos, TxPos};
use crate::hgvs::variant::HgvsVariant;

/// A bare-transcript intronic position found in a description.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BareTranscriptIntron {
    /// The coordinate system it was found on — `"c"` or `"n"`.
    pub coordinate_system: &'static str,
}

impl BareTranscriptIntron {
    /// The clause this finding cites.
    pub fn clause(&self) -> &'static str {
        "checklist.md:20"
    }

    /// A genomic-wrapper example on the same axis as the finding, so the
    /// message names a form the caller can actually write.
    pub fn wrapper_example(&self) -> &'static str {
        if self.coordinate_system == "n" {
            "NG_(NR_)/NC_(NR_)"
        } else {
            "NG_(NM_)/NC_(NM_)"
        }
    }

    /// The refusal sentence, shared by the parse and lenient-warning arms so
    /// the two cannot state the clause differently.
    pub fn refusal(&self) -> String {
        format!(
            "intronic offset on a bare {}. transcript reference; {} says a transcript \
             reference may describe an intronic variant only when a genomic reference \
             sequence is given on which it is annotated, e.g. {}",
            self.coordinate_system,
            self.clause(),
            self.wrapper_example(),
        )
    }
}

/// Find a bare-transcript intronic position anywhere in `variant`, walking
/// allele members.
///
/// This is the **description-wide** entry, used by the parse-stage rule: a
/// hygiene check that runs once per description rather than once per member
/// stops firing exactly on an allele, which is why the walk is here rather than
/// left to the caller. It delegates to [`bare_transcript_intronic_leaf`] so the
/// clause has one reading; see that function's note on why the split matters.
pub fn bare_transcript_intron(variant: &HgvsVariant) -> Option<BareTranscriptIntron> {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.iter().find_map(bare_transcript_intron),
        HgvsVariant::Supernumerary(inner) => bare_transcript_intron(inner),
        leaf => bare_transcript_intronic_leaf(leaf),
    }
}

/// The **per-leaf** predicate: the coordinate system (`"c"` or `"n"`) whose bare
/// transcript reference `variant` names an intronic position on, or `None`.
///
/// In scope: a bare coding transcript (`NM_`/`XM_` or LRG `LRG_<N>t<k>`) used
/// with `c.`, and a bare non-coding transcript (`NR_`/`XR_` via
/// `is_noncoding_rna()`, or LRG) used with `n.`, in both cases with
/// `Accession.genomic_context == None`. The spec's "a (non-)coding DNA
/// reference sequence does not contain introns" rule applies equally to curated
/// (`NM_`/`NR_`), predicted (`XM_`/`XR_`) and LRG transcript references — an LRG
/// transcript is itself a bare reference with no `NG_`/`NC_` genomic context —
/// so all are covered on each axis (#834).
///
/// Out of scope: genomic-context forms (`genomic_context: Some`), `NG_`/`NC_`
/// references (which never reach the `c.`/`n.` transcript path), Ensembl `ENST`,
/// and the `r.` axis.
///
/// Both `Single` and `Range` (uncertain-breakpoint) position boundaries are
/// inspected, so `c.(100+1_101-1)_(200+1_201-1)del` is covered; an unknown (`?`)
/// offset still counts as intronic (`CdsPos::is_intronic` treats the
/// unknown-offset sentinel as intronic), which is correct — it is an intronic
/// position whose exact offset is unspecified.
///
/// # One rule, three callers, and the split is deliberate
///
/// The scope prose above documents this function as much as its callers, and
/// that is load-bearing in three directions at once: the same question decides
/// whether **strict parse refuses an input** (#1630), whether **strict normalize
/// refuses one** (#486), and whether **#1704 must re-parent an output**. Two
/// readings of one clause is how ferro came to refuse a description in strict
/// mode while manufacturing the identical description in lenient. So every
/// caller derives from this predicate rather than restating it, which is what
/// makes them unable to drift.
///
/// The leaf/description split exists because the normalizer wants the leaf form
/// — `normalize_allele` already recurses per member, so a walking predicate
/// would ask the same question twice — while the parser holds the whole
/// description and has no other pass that would reach the members.
pub fn bare_transcript_intronic_leaf(variant: &HgvsVariant) -> Option<BareTranscriptIntron> {
    bare_transcript_intronic_axis(variant)
        .map(|coordinate_system| BareTranscriptIntron { coordinate_system })
}

/// The bare axis label, for callers that want the answer without building a
/// [`BareTranscriptIntron`].
///
/// The reason this exists separately is cost, not taste. Every `normalize()`
/// asks the question at least twice — once on the strict ladder and once at
/// `Normalizer::reparent_junction_exit` — on a path that wants a `bool`.
pub fn bare_transcript_intronic_axis(variant: &HgvsVariant) -> Option<&'static str> {
    match variant {
        HgvsVariant::Cds(v) => {
            // A bare coding-DNA reference does not contain introns, so an
            // intronic offset on it is a spec-invalid form regardless of which
            // transcript namespace addresses it. Curated/predicted RefSeq
            // (`NM_`/`XM_`) and LRG transcript references (`LRG_<N>t<k>`, which
            // carry no `NG_`/`NC_` genomic context — the LRG *is* the reference)
            // are all bare coding transcripts; treat them uniformly (#834).
            let is_bare_coding_transcript =
                matches!(&*v.accession.prefix, "NM" | "XM") || v.accession.is_lrg();
            if !is_bare_coding_transcript || v.accession.genomic_context.is_some() {
                return None;
            }
            let intronic = boundary_has_intronic(&v.loc_edit.location.start, CdsPos::is_intronic)
                || boundary_has_intronic(&v.loc_edit.location.end, CdsPos::is_intronic);
            intronic.then_some("c")
        }
        HgvsVariant::Tx(v) => {
            // Same rule on the non-coding axis: a bare `NR_`/`XR_` or LRG
            // non-coding transcript reference used with `n.` has no introns
            // (#834 extends the LRG coverage to match the `c.` arm).
            let is_bare_noncoding_transcript =
                v.accession.is_noncoding_rna() || v.accession.is_lrg();
            if !is_bare_noncoding_transcript || v.accession.genomic_context.is_some() {
                return None;
            }
            let intronic = boundary_has_intronic(&v.loc_edit.location.start, TxPos::is_intronic)
                || boundary_has_intronic(&v.loc_edit.location.end, TxPos::is_intronic);
            intronic.then_some("n")
        }
        // Enumerated rather than caught by `_`, to match the two rules that
        // share this seam (`genomic_offsets.rs`, `noncoding_zones.rs`) and for
        // the reason those give: a catch-all silently swallows a future
        // `HgvsVariant` variant, and the failure mode is this clause quietly
        // ceasing to reach a shape it should. The `r.` axis and the genomic
        // family carry no bare-transcript `c.`/`n.` position, and the wrapper
        // walk above has already taken `Allele`/`Supernumerary`.
        HgvsVariant::Genome(_)
        | HgvsVariant::Rna(_)
        | HgvsVariant::Protein(_)
        | HgvsVariant::Mt(_)
        | HgvsVariant::Circular(_)
        | HgvsVariant::GenomeRing(_)
        | HgvsVariant::RnaFusion(_)
        | HgvsVariant::Supernumerary(_)
        | HgvsVariant::Allele(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// True iff any concrete position reachable from `boundary` is intronic —
/// covering both a `Single` position and either endpoint of a `Range` boundary
/// (uncertain breakpoints like `c.(100+1_101-1)_(200+1_201-1)del`). `Unknown`
/// (`?`) and otherwise-absent inner positions contribute `false`.
fn boundary_has_intronic<T>(
    boundary: &UncertainBoundary<T>,
    is_intronic: impl Fn(&T) -> bool,
) -> bool {
    let mu_intronic = |mu: &crate::hgvs::uncertainty::Mu<T>| mu.inner().is_some_and(&is_intronic);
    match boundary {
        UncertainBoundary::Single(mu) => mu_intronic(mu),
        UncertainBoundary::Range { start, end } => mu_intronic(start) || mu_intronic(end),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;

    fn found(input: &str) -> Option<BareTranscriptIntron> {
        bare_transcript_intron(&parse_hgvs(input).expect("input parses"))
    }

    #[test]
    fn a_bare_coding_transcript_intronic_position_is_found() {
        let hit =
            found("NM_000088.3:c.589-1G>T").expect("a bare NM_ intronic position is in scope");
        assert_eq!(hit.coordinate_system, "c");
        assert_eq!(hit.clause(), "checklist.md:20");
        assert!(hit.refusal().contains("NG_(NM_)"));
    }

    #[test]
    fn a_bare_noncoding_transcript_intronic_position_is_found_on_the_n_axis() {
        let hit = found("NR_000019.1:n.20+2del").expect("a bare NR_ intronic position is in scope");
        assert_eq!(hit.coordinate_system, "n");
        assert!(hit.refusal().contains("NG_(NR_)"));
    }

    #[test]
    fn the_genomic_wrapper_form_the_clause_names_is_not_found() {
        assert_eq!(found("NC_000017.11(NM_000088.3):c.589-1G>T"), None);
    }

    #[test]
    fn an_exonic_position_on_a_bare_transcript_is_not_found() {
        // The control that keeps this a rule about INTRONIC offsets rather than
        // about bare accessions.
        assert_eq!(found("NM_000088.3:c.589G>T"), None);
    }

    #[test]
    fn an_allele_member_is_reached() {
        // The whole reason the description-wide walk exists: a check that ran
        // once per description without walking members would miss this.
        let hit = found("NM_000088.3:c.[589-1G>T;600G>A]").expect("an allele member is in scope");
        assert_eq!(hit.coordinate_system, "c");
    }

    #[test]
    fn an_allele_with_no_intronic_member_is_not_found() {
        assert_eq!(found("NM_000088.3:c.[589G>T;600G>A]"), None);
    }

    #[test]
    fn a_genomic_axis_description_is_not_this_clauses_business() {
        assert_eq!(found("NC_000001.11:g.12345A>G"), None);
    }
}
