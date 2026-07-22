//! Describing an all-silent change by the codons it rewrote (#1099).
//!
//! When a change leaves the protein identical to the reference, the description
//! must still say *what was interrogated*: `=` means a sequence "was tested but
//! found unchanged" (`general.md:91`), and a description "should include the
//! position(s) screened" (`SVD-WG001:17`). The spec's preferred form for a
//! predicted silent change is residue-scoped — `p.(Leu54=)`, not `p.Leu54Leu`
//! (`checklist.md:64`), and `SVD-WG001:54-55` localises `p.(Phe41=)` from the
//! DNA variant `c.123C>T`.
//!
//! # The rule
//!
//! Describe **the codons whose reference triplet was rewritten**, grouped the way
//! changed residues already are (`DNA/other.md:29-37`, which states this exact
//! three-tier structure one axis over):
//!
//! | rewritten codons | output |
//! |---|---|
//! | zero | `p.(=)` |
//! | one | `p.(Phe41=)` |
//! | one consecutive run | `p.(Phe41_Ser42=)` |
//! | separated runs | `p.[(Phe41=);(Gly47=)]` |
//!
//! The whole-molecule `p.(=)` is legal (`uncertain.md:212-213`) but least
//! informative, so it is kept for the case where it is also *accurate*: nothing
//! in the CDS was rewritten at all. One site emits it for a second reason — the
//! indel identity path has no fallback (declining there makes the projection
//! emit no `p.` at all), so it also renders `p.(=)` when the affected set cannot
//! be derived. Every other site declines instead.
//!
//! # Why triplets and never member spans
//!
//! The affected set is derived by comparing the mutated CDS against the reference
//! CDS triplet-by-triplet. Member spans are **not** usable, for two reasons:
//!
//! * **3'-shift.** `g.[4_6del;15_16insGCA]` normalizes to `c.[6_8del;16_18dup]`;
//!   the deletion of codon 2 becomes a span straddling codons 2 *and* 3, which
//!   would invent an `Ala3=`. Triplet comparison is shift-invariant — the two
//!   spellings splice byte-identical mutated CDSs.
//! * **Insertions have no codon.** A cis member's insertion span is collapsed to
//!   its anchor, so the inserted codon lands *between* residues.
//!
//! # Accepted limitation
//!
//! A codon interior to a rewritten span whose triplet is coincidentally unchanged
//! is not named: codons 41–43 with 42 unchanged give `p.[(Phe41=);(Cys43=)]`,
//! which is true but less complete than `p.(Phe41_Cys43=)`. Closing it needs
//! span-derived reasoning, which is exactly what the two bullets above forbid.

use crate::hgvs::edit::ProteinEdit;
use crate::hgvs::interval::ProtInterval;
use crate::hgvs::location::{AminoAcid, ProtPos};
use crate::hgvs::variant::{
    Accession, AllelePhase, AlleleVariant, HgvsVariant, LocEdit, ProteinVariant,
};

/// The 1-based codon positions whose reference triplet `mut_cds` rewrote.
///
/// Returns them ascending and distinct. Returns empty when nothing was rewritten
/// (the caller then describes whole-molecule identity) and, defensively, when the
/// two CDSs differ in length: this is the *silent* path, where a length change
/// cannot leave the protein identical, and a per-triplet comparison across a
/// length change would compare shifted frames.
///
/// A trailing partial codon (a CDS whose length is not a multiple of three) is
/// ignored — it translates to nothing, so it has no residue to name.
pub(crate) fn rewritten_codons(ref_cds: &str, mut_cds: &str) -> Vec<u64> {
    if ref_cds.len() != mut_cds.len() {
        return Vec::new();
    }
    ref_cds
        .as_bytes()
        .chunks_exact(3)
        .zip(mut_cds.as_bytes().chunks_exact(3))
        .enumerate()
        .filter(|(_, (ref_codon, mut_codon))| ref_codon != mut_codon)
        .map(|(i, _)| (i + 1) as u64)
        .collect()
}

/// Pair each rewritten codon with its reference amino acid.
///
/// `ref_residues` is indexed by codon position (1-based), so pass the translation
/// that **includes the terminator** where one is available: a silent change in the
/// stop codon rewrites a real triplet and is named `p.(Ter66=)`, and a translation
/// without the stop has no residue for it.
///
/// Returns `None` if any position has no reference residue — a codon past the end
/// of the translated protein, which happens when the annotated CDS runs past its
/// first stop. The caller decides what an underivable set means: sites that can
/// fall back decline, and the whole-molecule identity site keeps `p.(=)`.
pub(crate) fn codon_residues(
    codons: &[u64],
    ref_residues: &[AminoAcid],
) -> Option<Vec<(u64, AminoAcid)>> {
    codons
        .iter()
        .map(|&pos| {
            let aa = *ref_residues.get(usize::try_from(pos).ok()?.checked_sub(1)?)?;
            Some((pos, aa))
        })
        .collect()
}

/// Group items carrying ascending, distinct 1-based positions into runs of
/// consecutive positions. Items separated by an unused position land in
/// different runs.
pub(crate) fn group_consecutive_by<T>(items: Vec<T>, position: impl Fn(&T) -> u64) -> Vec<Vec<T>> {
    let mut groups: Vec<Vec<T>> = Vec::new();
    for item in items {
        match groups.last_mut() {
            Some(g) if position(&item) == position(g.last().expect("non-empty group")) + 1 => {
                g.push(item)
            }
            _ => groups.push(vec![item]),
        }
    }
    groups
}

/// Render the rewritten codons of an all-silent change as a protein consequence.
///
/// `residues` are the `(1-based position, reference amino acid)` pairs of the
/// rewritten codons, ascending and distinct — see [`rewritten_codons`] for how
/// the set must be derived. Empty renders the whole-molecule `p.(=)`; one run
/// renders a single residue or a range identity; separated runs are bracketed
/// into a cis allele, mirroring how `render_cis_substitution_groups` renders
/// separated *changed* residues.
///
/// `initiator` is the reference's residue 1, used only to anchor the
/// whole-molecule form. `Display` drops that position, but the variant is
/// serialized structurally through the Python API, so it must carry the
/// reference's own residue rather than a presumed `Met`: `RefProteinBundle`
/// forces `Met` only on a recognized initiator, leaving a drifted `cds_start`
/// untouched (#625). Pass `None` where no translated protein is in scope.
///
/// This is the one renderer for every site that can produce a silent
/// consequence, so the codon-level and residue-level combiners cannot drift.
pub(crate) fn render_silent_identity(
    residues: Vec<(u64, AminoAcid)>,
    initiator: Option<AminoAcid>,
    accession: &Accession,
    gene_symbol: &Option<String>,
) -> HgvsVariant {
    if residues.is_empty() {
        return HgvsVariant::Protein(ProteinVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::new_predicted(
                ProtInterval::point(ProtPos::new(initiator.unwrap_or(AminoAcid::Met), 1)),
                ProteinEdit::Identity {
                    predicted: false,
                    whole_protein: true,
                },
            ),
        });
    }

    let render_group = |g: &[(u64, AminoAcid)]| -> ProteinVariant {
        let (first_pos, first_aa) = g[0];
        let (last_pos, last_aa) = g[g.len() - 1];
        let location = if g.len() == 1 {
            ProtInterval::point(ProtPos::new(first_aa, first_pos))
        } else {
            ProtInterval::new(
                ProtPos::new(first_aa, first_pos),
                ProtPos::new(last_aa, last_pos),
            )
        };
        ProteinVariant {
            accession: accession.clone(),
            gene_symbol: gene_symbol.clone(),
            loc_edit: LocEdit::new_predicted(
                location,
                ProteinEdit::Identity {
                    predicted: false,
                    whole_protein: false,
                },
            ),
        }
    };

    let groups = group_consecutive_by(residues, |(pos, _)| *pos);
    if groups.len() == 1 {
        HgvsVariant::Protein(render_group(&groups[0]))
    } else {
        let variants: Vec<HgvsVariant> = groups
            .iter()
            .map(|g| HgvsVariant::Protein(render_group(g)))
            .collect();
        HgvsVariant::Allele(AlleleVariant::new(variants, AllelePhase::Cis))
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::project::accession::parse_accession;

    fn render(codons: &[u64], residues: &[AminoAcid]) -> String {
        let pairs = codon_residues(codons, residues).expect("residues cover every codon");
        render_silent_identity(pairs, None, &parse_accession("NP_000001.1"), &None).to_string()
    }

    /// Met-Leu-Arg-Leu-Ala-Ter, the shape the integration fixtures use.
    fn protein() -> Vec<AminoAcid> {
        vec![
            AminoAcid::Met,
            AminoAcid::Leu,
            AminoAcid::Arg,
            AminoAcid::Leu,
            AminoAcid::Ala,
            AminoAcid::Ter,
        ]
    }

    #[test]
    fn a_rewritten_triplet_is_reported_by_codon_position() {
        // Codon 2 CTG -> CTA and codon 4 CTC -> CTT; codon 3 untouched.
        assert_eq!(
            rewritten_codons("ATGCTGCGTCTCGCCTAA", "ATGCTACGTCTTGCCTAA"),
            vec![2, 4]
        );
    }

    #[test]
    fn an_unchanged_cds_rewrites_nothing() {
        assert!(rewritten_codons("ATGCTGCGTTAA", "ATGCTGCGTTAA").is_empty());
    }

    /// A span that straddles a codon boundary names only the codons it actually
    /// rewrote — the discriminator against deriving the set from member spans.
    ///
    /// Both spellings below splice the same mutated CDS, one within codon 2 and
    /// one straddling codons 2 and 3. A span-derived implementation would name
    /// `{2, 3}` for the straddling spelling and `{2}` for the other; the triplet
    /// rule names `{2}` for both, because codon 3 is `GCC` either way.
    #[test]
    fn a_straddling_span_names_only_the_codons_it_rewrote() {
        let reference = "ATGGCAGCCTAA"; // Met-Ala-Ala-Ter
        let within_codon_two = format!("{}GCC{}", &reference[..3], &reference[6..]); // c.4_6delinsGCC
        let straddling = format!("{}CG{}", &reference[..5], &reference[7..]); // c.6_7delinsCG
        assert_eq!(within_codon_two, straddling, "same mutated CDS either way");
        assert_eq!(rewritten_codons(reference, &within_codon_two), vec![2]);
        assert_eq!(rewritten_codons(reference, &straddling), vec![2]);
    }

    #[test]
    fn a_length_change_derives_no_codons() {
        assert!(rewritten_codons("ATGCTGCGTTAA", "ATGCTGTAA").is_empty());
    }

    #[test]
    fn a_trailing_partial_codon_is_ignored() {
        // Both CDSs carry a dangling 2-base tail that translates to nothing.
        assert_eq!(rewritten_codons("ATGCTGCG", "ATGCTACG"), vec![2]);
    }

    #[test]
    fn no_rewritten_codons_render_whole_molecule_identity() {
        assert_eq!(render(&[], &protein()), "NP_000001.1:p.(=)");
    }

    #[test]
    fn one_rewritten_codon_renders_a_single_residue() {
        assert_eq!(render(&[2], &protein()), "NP_000001.1:p.(Leu2=)");
    }

    #[test]
    fn a_consecutive_run_renders_as_a_range() {
        assert_eq!(render(&[2, 3], &protein()), "NP_000001.1:p.(Leu2_Arg3=)");
        assert_eq!(render(&[2, 3, 4], &protein()), "NP_000001.1:p.(Leu2_Leu4=)");
    }

    #[test]
    fn separated_runs_render_as_a_cis_bracket() {
        assert_eq!(
            render(&[2, 5], &protein()),
            "NP_000001.1:p.[(Leu2=);(Ala5=)]"
        );
        assert_eq!(
            render(&[2, 3, 5], &protein()),
            "NP_000001.1:p.[(Leu2_Arg3=);(Ala5=)]"
        );
    }

    #[test]
    fn a_silent_change_in_the_stop_codon_names_the_terminator() {
        assert_eq!(render(&[6], &protein()), "NP_000001.1:p.(Ter6=)");
    }

    #[test]
    fn a_codon_past_the_translated_protein_has_no_residue() {
        assert!(codon_residues(&[2, 9], &protein()).is_none());
    }

    /// The gene symbol is carried onto every rendered variant. Whether it is
    /// *emitted* is the renderer's caller's business: a single protein variant
    /// suppresses it (#310 — the p. grammar has no selector production) while an
    /// allele emits it, exactly as the changed-residue sibling renderer behaves.
    /// The whole-molecule form anchors on the REFERENCE's residue 1, not a
    /// presumed `Met`. `Display` drops the position, but the variant is
    /// serialized structurally through the Python API, and a transcript whose
    /// `cds_start` drifted off a recognized initiator keeps its own residue 1
    /// (#625) — so a fabricated `Met` would be a silent output divergence.
    #[test]
    fn whole_molecule_identity_anchors_on_the_reference_initiator() {
        let anchored = render_silent_identity(
            Vec::new(),
            Some(AminoAcid::Ser),
            &parse_accession("NP_000001.1"),
            &None,
        );
        assert_eq!(anchored.to_string(), "NP_000001.1:p.(=)");
        let HgvsVariant::Protein(pv) = &anchored else {
            panic!("expected a protein variant");
        };
        let start = pv
            .loc_edit
            .location
            .start
            .inner()
            .expect("certain position");
        assert_eq!(start.aa, AminoAcid::Ser, "anchored on a presumed Met");

        // With no protein in scope the presumed initiator is the fallback.
        let unanchored =
            render_silent_identity(Vec::new(), None, &parse_accession("NP_000001.1"), &None);
        let HgvsVariant::Protein(pv) = &unanchored else {
            panic!("expected a protein variant");
        };
        assert_eq!(
            pv.loc_edit
                .location
                .start
                .inner()
                .expect("certain position")
                .aa,
            AminoAcid::Met
        );
    }

    #[test]
    fn the_gene_symbol_is_carried_onto_every_rendered_variant() {
        let single = render_silent_identity(
            codon_residues(&[2, 3], &protein()).expect("covered"),
            None,
            &parse_accession("NP_000001.1"),
            &Some("GENE1".to_string()),
        );
        assert_eq!(single.gene_symbol(), Some("GENE1"));
        assert_eq!(single.to_string(), "NP_000001.1:p.(Leu2_Arg3=)");

        let bracketed = render_silent_identity(
            codon_residues(&[2, 5], &protein()).expect("covered"),
            None,
            &parse_accession("NP_000001.1"),
            &Some("GENE1".to_string()),
        );
        assert_eq!(bracketed.to_string(), "NP_000001.1:p.[(Leu2=);(Ala5=)]");
    }
}
