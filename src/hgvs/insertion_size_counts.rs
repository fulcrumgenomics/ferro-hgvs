//! `checklist.md:33`'s size-count insertions, and where a description can still state one.
//!
//! # The rule
//!
//! `docs/recommendations/checklist.md` item 3 ("Insertions") asks, at `:32`,
//! "**do you provide the inserted sequence?**", and answers at `:33`:
//!
//! > Describing a variant as `c.5439_5430ins6` is not allowed, the inserted
//! > sequence (for `ins6`, e.g., `TGCCAT`) should be specified.
//!
//! `:33` prints its example inside `class="invalid"`, and "not allowed" is
//! MUST-level under the spec's RFC 2119 reading (`recommendations/style.md:9`).
//! Its sibling `:49` is the deletion case ferro has enforced since #1079
//! (`validate_no_point_size_suffix` in the parser), so before #1789 one numbered
//! checklist item was enforced for `del`/`dup` and not for `ins`.
//!
//! # Two axis clauses say it constructively, which is what scopes the rule
//!
//! `checklist.md` is not molecule-specific, and a `DNA/` clause cannot scope
//! `r.` (see `CLAUDE.md`, "Cite the clause exactly"). Both axes state the rule
//! in their own jurisdiction, and state it the same way — by enumerating what an
//! inserted sequence may be, with no bare count among the alternatives:
//!
//! - `DNA/insertion.md:22` — the "inserted_sequence" can be given as *the
//!   nucleotides inserted* (`insAGC`), *a range of the same reference*
//!   (`c.849_850ins858_895`), or *a range of another reference*.
//! - `RNA/insertion.md:20` — the same sentence, with `insagc` and
//!   `r.849_850ins858_895`.
//!
//! And both name the conformant spelling for the case a size-count is reaching
//! for — a known *number* of nucleotides whose identity is unknown:
//! `DNA/insertion.md:77` gives `g.32717298_32717299insN[100]` ("the insertion of
//! 100 nucleotides (not specified)") and `RNA/insertion.md:41` gives
//! `r.1149_1150insn[100]`. So refusing `ins6` costs no expressive power: it has
//! an exact conformant translation, `insN[6]`, which ferro already parses,
//! normalizes and re-emits.
//!
//! `DNA/insertion.md:119` reaches the same verdict about a real reported
//! description: "`c.23ins24` is not correct since the position of the insertion
//! is not described properly and because \"ins24\" does not define the sequence
//! inserted."
//!
//! # Where the shape actually survives, which is narrower than the text reads
//!
//! Two surface spellings reach [`InsertedSequence::Count`], through two
//! different arms of `parse_inserted_sequence`:
//!
//! | input | arm | result |
//! |---|---|---|
//! | `c.10_11ins6` | `parse_simple_count` | `Count(6)` |
//! | `c.10_11ins(6)` | the parenthesized arm | `Count(6)` |
//!
//! That is why this predicate is keyed on the **AST** and not on the rendered
//! description. A textual `ins<digits>` scan would have caught one spelling and
//! missed the other, and — worse — would have false-positived on
//! `c.849_850ins858_895`, the range form `DNA/insertion.md:22` explicitly
//! sanctions, which is also `ins` followed by digits.
//!
//! # What it deliberately does not reach
//!
//! **[`InsertedSequence::Range`] (`ins(10_20)`).** In ferro's AST that node is
//! ambiguous between "a count range" and "an uncertain *position* range", and
//! the second is spec-sanctioned — `LRG_308:g.?_?ins(23632682_23625413)_…` is a
//! position range under `DNA/complex.md`, and its single-bracket sibling parses
//! to the same `Range`. Refusing the node would risk refusing a conformant
//! description, so widening to it needs its own measurement and its own record.
//! (The conformant spelling for a count range is `insN[(80_120)]`,
//! `DNA/insertion.md:80`.) Pinned by
//! `tests/it/issue_1789_insertion_size_count.rs::`
//! `a_count_range_insert_is_deliberately_untouched`.
//!
//! **`delins<number>`.** `checklist.md:33` sits under item 3, "Insertions", and
//! `DNA/delins.md` states no equivalent prohibition — its `class="invalid"`
//! notes at `:29` and `:34` are about restating the *deleted* sequence, not
//! about the insert's form. Extending there would be an adjudication no clause
//! supports, so `Delins` and `DupIns` are classified `None` below with that
//! reason attached.
//!
//! **`inv<number>`.** [`NaEdit::Inversion`] carries an explicit `length` field
//! and a range anchor already names both endpoints, so `:49`'s "name the first
//! and last residue" reasoning does not reach it and no clause names an
//! `inv<number>` form at all. It is also not a live output defect: normalization
//! drops the count. Measured in
//! `tests/it/issue_1789_insertion_size_count.rs::`
//! `an_inversion_sized_by_a_number_is_a_different_shape`.
//!
//! # The stage
//!
//! Set by the decided `rulings[absolute-prohibition-enforcement-stage]`, which
//! names `checklist.md:33` in its own clause list: strict fails at PARSE,
//! lenient does not validate input conformance and fails only when it cannot
//! NORMALIZE, silent is lenient without messages. The normalize rung is not
//! mode-gated, because rule 1 of the README ruleset is about *output* and no
//! mode may trade it. See `parser::apply_insertion_size_count_rule` and
//! `Normalizer::normalize_core_checked`.

use crate::hgvs::edit::{InsertedSequence, NaEdit};
use crate::hgvs::variant::HgvsVariant;

/// An insertion whose payload is a bare count of unspecified nucleotides.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InsertionSizeCount {
    /// The count the description states, e.g. `6` for `ins6`.
    pub count: u64,
}

impl InsertionSizeCount {
    /// The clause carrying the prohibition, for a diagnostic.
    #[must_use]
    pub fn clause(&self) -> &'static str {
        "checklist.md:33"
    }

    /// The conformant translation of this description's payload — the spelling
    /// `DNA/insertion.md:77` and `RNA/insertion.md:41` use for a known number of
    /// unspecified nucleotides.
    ///
    /// Rendered on the DNA axis (`N`); an `r.` description spells the same
    /// symbol in lower case (`general.md:50`), which
    /// [`Self::conformant_spelling_for_rna`] supplies.
    #[must_use]
    pub fn conformant_spelling(&self) -> String {
        format!("insN[{}]", self.count)
    }

    /// The `r.` spelling of [`Self::conformant_spelling`]. `general.md:50` says
    /// an RNA description states its nucleotides in lower case, so telling an
    /// `r.` author to write `insN[6]` would prescribe the one thing that axis
    /// forbids.
    #[must_use]
    pub fn conformant_spelling_for_rna(&self) -> String {
        format!("insn[{}]", self.count)
    }
}

/// The clause enumerating what an inserted sequence may be, in the axis's own
/// jurisdiction.
///
/// A `DNA/` clause cannot scope `r.` (module docs), so a diagnostic aimed at an
/// `r.` author must cite `RNA/insertion.md:20` rather than its DNA twin. Same
/// discipline the `W3028` sibling applies — its lower-case RNA arm swaps the
/// citation as well as the suggested symbol.
#[must_use]
pub fn payload_forms_clause(rna_axis: bool) -> &'static str {
    if rna_axis {
        "RNA/insertion.md:20"
    } else {
        "DNA/insertion.md:22"
    }
}

/// The clause naming the spelling for a known number of *unspecified*
/// nucleotides — the conformant translation of the refused payload — in the
/// axis's own jurisdiction.
#[must_use]
pub fn unspecified_run_clause(rna_axis: bool) -> &'static str {
    if rna_axis {
        "RNA/insertion.md:41"
    } else {
        "DNA/insertion.md:77"
    }
}

/// The worked "range of the same reference" example each axis clause publishes,
/// so the diagnostic offers one spelled on the axis the author is writing.
#[must_use]
pub fn reference_range_example(rna_axis: bool) -> &'static str {
    if rna_axis {
        "r.849_850ins858_895"
    } else {
        "c.849_850ins858_895"
    }
}

/// The count a nucleic-acid edit's insertion payload states, if any.
///
/// Exhaustive rather than wildcarded, for the same reason
/// [`crate::hgvs::alignment_symbols`]'s sibling is: the four arms carrying an
/// [`InsertedSequence`] are enumerated, and a new [`NaEdit`] variant has to be
/// classified here rather than silently defaulting to "states a sequence".
///
/// `Delins` and `DupIns` also carry an `InsertedSequence` and are deliberately
/// **not** matched — see the module docs. `Inversion` and `Duplication` carry a
/// `length` of their own, which is a different clause (`checklist.md:49`,
/// enforced by `validate_no_point_size_suffix`) and a different geometry.
fn count_in_edit(edit: &NaEdit) -> Option<InsertionSizeCount> {
    match edit {
        NaEdit::Insertion { sequence } | NaEdit::BreakpointInsertion { sequence } => {
            match sequence {
                InsertedSequence::Count(n) => Some(InsertionSizeCount { count: *n }),
                // `Range` is the ambiguous node (module docs); everything else
                // either states bases, names a position range, or states
                // nothing at all — and `Empty` is W3027's business.
                InsertedSequence::Literal(_)
                | InsertedSequence::Range(_, _)
                | InsertedSequence::Repeat { .. }
                | InsertedSequence::SequenceRepeat { .. }
                | InsertedSequence::Complex(_)
                | InsertedSequence::Named(_)
                | InsertedSequence::Reference(_)
                | InsertedSequence::PositionRange { .. }
                | InsertedSequence::PositionRangeInv { .. }
                | InsertedSequence::UncertainRangeInv { .. }
                | InsertedSequence::SpecialPositionRange { .. }
                | InsertedSequence::Uncertain
                | InsertedSequence::Empty => None,
            }
        }
        NaEdit::Delins { .. }
        | NaEdit::DupIns { .. }
        | NaEdit::Substitution { .. }
        | NaEdit::SubstitutionNoRef { .. }
        | NaEdit::Deletion { .. }
        | NaEdit::NPaddedDeletion { .. }
        | NaEdit::Duplication { .. }
        | NaEdit::Inversion { .. }
        | NaEdit::Repeat { .. }
        | NaEdit::MultiRepeat { .. }
        | NaEdit::Identity { .. }
        | NaEdit::Conversion { .. }
        | NaEdit::Unknown { .. }
        | NaEdit::Methylation { .. }
        | NaEdit::CopyNumber { .. }
        | NaEdit::Splice { .. }
        | NaEdit::NoProduct
        | NaEdit::PositionOnly => None,
    }
}

/// The first size-count insertion `variant` states, if any.
///
/// Walks every nucleic-acid member, including allele members, ring segments and
/// the inner description of a `sup` marker, so the spelling cannot hide behind a
/// composite one — a cis allele (`c.[10_11ins6;20del]`), a trans pair
/// (`c.[10_11ins6];[20del]`) and an uncertain group (`c.(10_11ins6)`) are all
/// reached.
///
/// The match is exhaustive rather than wildcarded so a new axis must be
/// classified here. The four that yield `None` do so for stated reasons:
///
/// - `Protein` carries a [`ProteinInsSeq`](crate::hgvs::edit::ProteinInsSeq),
///   not an [`InsertedSequence`], and that type has no count variant — a protein
///   insertion states 3-letter codes or nothing. There is no shape here to
///   classify.
/// - `RnaFusion` holds two breakpoints — an accession, an optional gene symbol
///   and an interval — and **no edit at all**.
/// - `NullAllele` (`0`) and `UnknownAllele` (`?`) are whole-allele markers with
///   neither position nor edit.
///
/// Unlike [`crate::hgvs::alignment_symbols::alignment_only_symbol`] this is not
/// axis-keyed: the prohibition is on the payload's *form*, not on an alphabet,
/// so `DNA/insertion.md:22` and `RNA/insertion.md:20` yield the same verdict and
/// only the suggested spelling differs (`insN[6]` vs `insn[6]`).
#[must_use]
pub fn insertion_size_count(variant: &HgvsVariant) -> Option<InsertionSizeCount> {
    match variant {
        HgvsVariant::Genome(v) => v.loc_edit.edit.inner().and_then(count_in_edit),
        HgvsVariant::Cds(v) => v.loc_edit.edit.inner().and_then(count_in_edit),
        HgvsVariant::Tx(v) => v.loc_edit.edit.inner().and_then(count_in_edit),
        HgvsVariant::Rna(v) => v.loc_edit.edit.inner().and_then(count_in_edit),
        HgvsVariant::Mt(v) => v.loc_edit.edit.inner().and_then(count_in_edit),
        HgvsVariant::Circular(v) => v.loc_edit.edit.inner().and_then(count_in_edit),
        HgvsVariant::Allele(allele) => allele.variants.iter().find_map(insertion_size_count),
        HgvsVariant::GenomeRing(ring) => ring
            .segments
            .iter()
            .find_map(|segment| segment.edit.inner().and_then(count_in_edit)),
        HgvsVariant::Supernumerary(inner) => insertion_size_count(inner),
        HgvsVariant::Protein(_)
        | HgvsVariant::RnaFusion(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// True when `variant` is on the `r.` axis, so a diagnostic can offer the
/// lower-case spelling `general.md:50` requires there — and cite that axis's own
/// clauses (see [`payload_forms_clause`]).
///
/// Resolved on the **member**, not on an enclosing group, so a mixed composite
/// cannot inherit the wrong case. An allele's members are full variants, so the
/// first member carrying the offending payload is the one asked — the same
/// member [`insertion_size_count`] reports, since both walk in order and stop at
/// the first hit.
///
/// Exhaustive rather than wildcarded, for the same reason [`count_in_edit`] and
/// [`insertion_size_count`] are: a new axis must be classified here rather than
/// silently defaulting to the DNA spelling and the DNA citations. `Tx` (`n.`)
/// and `Mt`/`Circular` are DNA axes and take the upper-case form; `Protein`,
/// `RnaFusion` and the two whole-allele markers carry no
/// [`InsertedSequence`](crate::hgvs::edit::InsertedSequence) to reach here at
/// all.
#[must_use]
pub fn states_rna_axis(variant: &HgvsVariant) -> bool {
    match variant {
        HgvsVariant::Rna(_) => true,
        HgvsVariant::Allele(allele) => allele
            .variants
            .iter()
            .find(|member| insertion_size_count(member).is_some())
            .is_some_and(states_rna_axis),
        HgvsVariant::Supernumerary(inner) => states_rna_axis(inner),
        HgvsVariant::Genome(_)
        | HgvsVariant::Cds(_)
        | HgvsVariant::Tx(_)
        | HgvsVariant::Mt(_)
        | HgvsVariant::Circular(_)
        | HgvsVariant::GenomeRing(_)
        | HgvsVariant::Protein(_)
        | HgvsVariant::RnaFusion(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;

    fn count_of(input: &str) -> Option<InsertionSizeCount> {
        insertion_size_count(&parse_hgvs(input).expect("input parses"))
    }

    #[test]
    fn a_bare_count_insert_is_found() {
        assert_eq!(
            count_of("NC_TEST.1:g.10_11ins6"),
            Some(InsertionSizeCount { count: 6 })
        );
    }

    #[test]
    fn the_parenthesized_spelling_reaches_the_same_node() {
        assert_eq!(
            count_of("NC_TEST.1:g.10_11ins(6)"),
            Some(InsertionSizeCount { count: 6 })
        );
    }

    /// The range form `DNA/insertion.md:22` sanctions is also `ins` followed by
    /// digits, so this is the false positive a textual scan would produce.
    #[test]
    fn a_position_range_insert_is_not_a_count() {
        assert_eq!(count_of("NC_TEST.1:g.10_11ins858_895"), None);
    }

    /// The conformant translation of the very shape this rule refuses.
    #[test]
    fn the_conformant_unspecified_run_is_not_a_count() {
        assert_eq!(count_of("NC_TEST.1:g.10_11insN[6]"), None);
    }

    #[test]
    fn a_literal_insert_is_not_a_count() {
        assert_eq!(count_of("NC_TEST.1:g.10_11insTGCCAT"), None);
    }

    /// `Range` is the deliberate exclusion; see the module docs.
    #[test]
    fn a_count_range_is_deliberately_not_matched() {
        assert_eq!(count_of("NC_TEST.1:g.10_11ins(6_9)"), None);
    }

    /// `delins` carries the same node and is out of jurisdiction.
    #[test]
    fn a_delins_count_is_out_of_scope() {
        assert_eq!(count_of("NC_TEST.1:g.10_11delins6"), None);
    }

    #[test]
    fn a_member_of_a_composite_is_reached() {
        assert_eq!(
            count_of("NC_TEST.1:g.[10_11ins6;20del]"),
            Some(InsertionSizeCount { count: 6 })
        );
        assert_eq!(
            count_of("NC_TEST.1:g.[10_11ins6];[20del]"),
            Some(InsertionSizeCount { count: 6 })
        );
    }

    #[test]
    fn the_rna_axis_takes_the_lower_case_spelling() {
        let found = count_of("NM_TEST.1:r.10_11ins6").expect("stated");
        assert_eq!(found.conformant_spelling(), "insN[6]");
        assert_eq!(found.conformant_spelling_for_rna(), "insn[6]");
        assert!(states_rna_axis(
            &parse_hgvs("NM_TEST.1:r.10_11ins6").expect("parses")
        ));
        assert!(!states_rna_axis(
            &parse_hgvs("NC_TEST.1:g.10_11ins6").expect("parses")
        ));
    }
}
