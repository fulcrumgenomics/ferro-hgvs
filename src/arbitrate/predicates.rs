//! Deterministic, normalizer-independent compliance predicates (§5.1.4).
//! Only rules that are a direct reference check live here; the 3'-shift rule
//! is intentionally NOT here (it reduces to `shuffle`; see §5.1.4/§A2).

use crate::arbitrate::category::ArbitrationCategory;
use crate::arbitrate::oracle::{edit_tuple, EditTuple};
use crate::hgvs::edit::NaEdit;
use crate::hgvs::variant::HgvsVariant;
use crate::reference::provider::ReferenceProvider;

/// Spec (`duplication.md`): "when a variant can be described as a duplication
/// it must be described as a duplication and not as, e.g., an insertion."
///
/// Detects the case by direct reference comparison — no [`crate::normalize`]
/// involved. Both variants are first reduced to reference-anchored SPDI
/// tuples ([`edit_tuple`], the same normalizer-independent reduction the
/// same/different oracle uses); the rule applies only when both reduce to an
/// identical pure insertion (same accession, same inserted bases, no deleted
/// bases) whose inserted bases equal the reference bases immediately 5' of
/// the insertion point — i.e. the edit is a tandem duplication regardless of
/// how either tool spelled it.
///
/// Returns `Some(FerroCorrect)` when ferro wrote the edit as a duplication
/// (`NaEdit::Duplication`/`NaEdit::DupIns`) and the other tool wrote the
/// equivalent insertion, `Some(MutalyzerCorrect)` for the mirror image, and
/// `None` when the rule does not apply: the two edits are not the same
/// insertion, the inserted bases are not a tandem copy of the directly-5'
/// reference context, or neither/both tools wrote it as a `dup`.
pub fn dup_vs_ins_compliance(
    ferro: &HgvsVariant,
    other: &HgvsVariant,
    provider: &(impl ReferenceProvider + ?Sized),
) -> Option<ArbitrationCategory> {
    let tf = edit_tuple(ferro, provider).ok()?;
    let to = edit_tuple(other, provider).ok()?;

    // Rule applies only to two spellings of the SAME pure insertion: same
    // accession, same interbase position, identical (non-empty) inserted
    // bases, and no deleted bases on either side (a `dup`/`ins` never
    // deletes reference bases; ruling out `deleted` also excludes delins).
    if tf.accession != to.accession
        || tf.position != to.position
        || tf.inserted != to.inserted
        || tf.inserted.is_empty()
        || !tf.deleted.is_empty()
        || !to.deleted.is_empty()
    {
        return None;
    }

    if !is_duplication_of_5prime_copy(&tf, provider) {
        return None; // inserted bases are not a directly-5' copy -> rule silent
    }

    match (is_written_as_dup(ferro), is_written_as_dup(other)) {
        (true, false) => Some(ArbitrationCategory::FerroCorrect),
        (false, true) => Some(ArbitrationCategory::MutalyzerCorrect),
        _ => None,
    }
}

/// True if the inserted bases equal the reference bases immediately 5' of the
/// insertion point (i.e. the edit is a tandem duplication). Direct reference
/// comparison via [`ReferenceProvider::get_sequence`] — no `Normalizer`.
fn is_duplication_of_5prime_copy(
    t: &EditTuple,
    provider: &(impl ReferenceProvider + ?Sized),
) -> bool {
    let n = t.inserted.len() as u64;
    if n == 0 || t.position < n {
        return false; // guard against underflow when fetching the flank
    }
    let preceding = provider
        .get_sequence(&t.accession, t.position - n, t.position)
        .unwrap_or_default();
    preceding.eq_ignore_ascii_case(&t.inserted)
}

/// The certain inner [`NaEdit`] of an NaEdit-bearing variant, or `None` for
/// non-NaEdit variants (e.g. protein, fusion) and uncertain/unknown edits.
fn na_edit(variant: &HgvsVariant) -> Option<&NaEdit> {
    match variant {
        HgvsVariant::Genome(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Cds(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Tx(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Rna(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Mt(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Circular(v) => v.loc_edit.edit.inner(),
        _ => None,
    }
}

/// True when `v`'s edit is typed as a duplication (`dup`) or the
/// dup-then-insert compound (`dupins`) — a typed check on the parsed
/// [`NaEdit`], not a string match on the rendered HGVS text.
fn is_written_as_dup(v: &HgvsVariant) -> bool {
    matches!(
        na_edit(v),
        Some(NaEdit::Duplication { .. } | NaEdit::DupIns { .. })
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::mock::MockProvider;

    /// A 10 bp synthetic genomic contig: `ATG|CG|CG|TAA` (1-based 1..=10).
    /// The tandem `CG|CG` at 1-based 4..=7 makes a duplication of the 2-base
    /// unit unambiguous: `g.4_5dup` (duplicating the first `CG`, 1-based
    /// 4..=5) and `g.5_6insCG` (inserting `CG` between 1-based 5 and 6) both
    /// reduce to the identical SPDI insertion `NC_000001.11:5::CG`, and the
    /// inserted `CG` equals the reference bases immediately 5' of position 5
    /// (1-based 4..=5) — the tandem-duplication condition the rule keys on.
    fn provider() -> MockProvider {
        let mut p = MockProvider::new();
        let seq = "ATGCGCGTAA";
        assert_eq!(&seq[3..5], "CG"); // 0-based [3,5) == 1-based 4..=5
        p.add_genomic_sequence("NC_000001.11", seq);
        p
    }

    #[test]
    fn ins_of_directly_3prime_copy_must_be_dup_ferro_correct() {
        let p = provider();
        let ferro = parse_hgvs("NC_000001.11:g.4_5dup").unwrap(); // dup of CG
        let other = parse_hgvs("NC_000001.11:g.5_6insCG").unwrap(); // same edit, as ins
        assert_eq!(
            dup_vs_ins_compliance(&ferro, &other, &p),
            Some(ArbitrationCategory::FerroCorrect)
        );
    }

    #[test]
    fn mutalyzer_correct_when_other_writes_the_dup() {
        let p = provider();
        // Mirror image of the case above: ferro spells the same edit as an
        // insertion while the other tool correctly spells it as a dup.
        let ferro = parse_hgvs("NC_000001.11:g.5_6insCG").unwrap();
        let other = parse_hgvs("NC_000001.11:g.4_5dup").unwrap();
        assert_eq!(
            dup_vs_ins_compliance(&ferro, &other, &p),
            Some(ArbitrationCategory::MutalyzerCorrect)
        );
    }

    #[test]
    fn non_flanking_insertion_rule_does_not_apply() {
        let p = provider();
        let ferro = parse_hgvs("NC_000001.11:g.4_5dup").unwrap();
        let other = parse_hgvs("NC_000001.11:g.5_6insTT").unwrap(); // not a dup; different edit
        assert_eq!(dup_vs_ins_compliance(&ferro, &other, &p), None);
    }

    #[test]
    fn flank_mismatch_guard_returns_none_even_for_same_insertion() {
        // Same accession/position/inserted-bases/empty-deletions as the
        // matching-pair tests above, so this reaches
        // `is_duplication_of_5prime_copy` rather than short-circuiting on
        // one of the earlier equality checks (the coverage gap
        // `non_flanking_insertion_rule_does_not_apply` above leaves: that
        // test differs in *inserted bases*, so it never reaches the guard).
        //
        // `g.7_8dupTT` (explicit-sequence dup, trusted as given — not
        // reference-verified) and `g.8_9insTT` both reduce to the identical
        // edit tuple `NC_000001.11:8::TT` (dup's SPDI position is
        // `end_one_based` = 8; ins's SPDI position is `start_one_based` = 8
        // for `g.8_9insTT`). But the reference bases directly 5' of
        // position 8 (1-based 7..=8) are `GT`, not `TT` — the inserted
        // bases are NOT a copy of the 5' flank, so the guard must return
        // `false` and the rule stays silent (`None`), even though `ferro`
        // is typed as a genuine `NaEdit::Duplication`: a broken guard would
        // let the `is_written_as_dup` match below decide instead and
        // return `Some(FerroCorrect)`.
        let p = provider();
        let ferro = parse_hgvs("NC_000001.11:g.7_8dupTT").unwrap();
        let other = parse_hgvs("NC_000001.11:g.8_9insTT").unwrap();
        assert_eq!(dup_vs_ins_compliance(&ferro, &other, &p), None);
    }
}
