//! Issue #1099: an all-silent change must name the codons it interrogated.
//!
//! When a change leaves the protein identical to the reference, ferro collapsed
//! the consequence to a single arbitrary residue's `=` — dropping every other
//! silent codon — and a silent `delins` escaped to the whole-molecule `p.(=)`.
//!
//! The rule (see the issue's verdict comment): describe the codons whose
//! reference triplet was rewritten, grouped the way changed residues already
//! are — one codon → `p.(Xxx=)`, a consecutive run → a range identity, and
//! separated runs → a cis bracket. `p.(=)` is reserved for the case where
//! nothing in the CDS was rewritten at all, because only then does it not
//! over-claim the scope of what was analysed (`general.md:91`,
//! `SVD-WG001:17`, `checklist.md:64`).

use ferro_hgvs::data::{CdotMapper, CdotTranscript, Projector};
use ferro_hgvs::project::VariantProjector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// Build a single-exon plus-strand coding transcript on contig `REF` whose CDS
/// is `cds`, placed so that `g.N == c.N`.
fn projector_for(cds: &str) -> VariantProjector<MockProvider> {
    let len = cds.len() as u64;
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "tx".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("GENE1099".to_string()),
            contig: "REF".to_string(),
            strand: Strand::Plus,
            exons: vec![[1, len + 1, 0, len]],
            cds_start: Some(0),
            cds_end: Some(len),
            gene_id: None,
            protein: Some("txp".to_string()),
            exon_cigars: Vec::new(),
        },
    );

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "tx".to_string(),
        Some("GENE1099".to_string()),
        TxStrand::Plus,
        cds.to_string(),
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        Some("REF".to_string()),
        Some(1),
        Some(len),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    provider.add_genomic_sequence("REF", format!("{cds}NN"));

    VariantProjector::new(Projector::new(cdot), provider)
}

/// The issue's own transcript: Met-Leu-Arg-Trp-Ala-Stop.
/// Codon 2 `CTG`, codon 3 `CGT`, codon 5 `GCC` all have synonymous third-base
/// alternatives; codon 4 `TGG` (Trp) has none, which is why the separated-run
/// case below skips to codon 5.
const ISSUE_CDS: &str = "ATGCTGCGTTGGGCCTAA";

/// Met-Leu-Arg-Leu-Ala-Stop — four consecutive degenerate codons, so a run of
/// three silent codons is constructible (the issue's CDS cannot express one,
/// because `TGG` blocks codon 4).
const RUN_CDS: &str = "ATGCTGCGTCTCGCCTAA";

fn protein_of(projector: &VariantProjector<MockProvider>, input: &str) -> String {
    projector
        .project(input, "tx")
        .unwrap_or_else(|e| panic!("{input} failed to project: {e}"))
        .protein
        .unwrap_or_else(|| panic!("{input} produced no protein consequence"))
        .to_string()
}

/// Every all-silent cis allele names each codon it rewrote, grouped into runs.
/// The first row is the input the issue was filed on.
#[test]
fn all_silent_cis_allele_names_every_rewritten_codon() {
    let vp = projector_for(ISSUE_CDS);
    for (input, expected) in [
        // Issue headline: three members, two codons (2 and 3), consecutive.
        ("REF:g.[4C>T;6G>A;9T>C]", "txp:p.(Leu2_Arg3=)"),
        // Two members, one per codon — the case the issue called "fine" and
        // which regressed to a single residue under #1083.
        ("REF:g.[6G>A;9T>C]", "txp:p.(Leu2_Arg3=)"),
        ("REF:g.[4C>T;9T>C]", "txp:p.(Leu2_Arg3=)"),
    ] {
        assert_eq!(protein_of(&vp, input), expected, "input: {input}");
    }
}

/// A run of three consecutive silent codons renders as one range identity —
/// not the two endpoints, and not a three-member bracket.
#[test]
fn three_consecutive_silent_codons_render_as_one_range() {
    let vp = projector_for(RUN_CDS);
    assert_eq!(
        protein_of(&vp, "REF:g.[6G>A;9T>C;12C>A]"),
        "txp:p.(Leu2_Leu4=)"
    );
}

/// Silent codons separated by an untouched codon stay individual, bracketed —
/// the same rule `delins.md` applies to separated *changed* residues.
///
/// The accession renders bare: the `p.` axis takes no gene-symbol selector in
/// any shape (#1142), so an allele drops it exactly as a single protein variant
/// always has.
#[test]
fn separated_silent_codons_render_as_a_bracket() {
    let vp = projector_for(ISSUE_CDS);
    assert_eq!(
        protein_of(&vp, "REF:g.[6G>A;15C>A]"),
        "txp:p.[(Leu2=);(Ala5=)]"
    );
}

/// The **codon-level** combiner groups separated codons too.
///
/// Two members sharing a codon route the allele to
/// `combine_cis_substitutions_by_codon` rather than
/// `combine_cis_members_by_residue` — a separate implementation that derives its
/// codon set one reference codon at a time instead of from a whole mutated CDS.
/// Every other bracket assertion here goes through the residue-level site, so
/// without this case a grouping regression in the codon-level path would ship
/// green. `c.4C>T` and `c.6G>A` share codon 2; `c.15C>A` rewrites codon 5.
#[test]
fn the_codon_level_combiner_also_brackets_separated_codons() {
    let vp = projector_for(ISSUE_CDS);
    assert_eq!(
        protein_of(&vp, "REF:g.[4C>T;6G>A;15C>A]"),
        "txp:p.[(Leu2=);(Ala5=)]"
    );
}

/// Silent and changed residues render their accession the same way, so the
/// identity renderer cannot drift from the substitution renderer it mirrors.
///
/// Both suppress the gene-symbol selector, in both the bracketed-allele and the
/// single-variant shape — the `p.` axis takes none (#310/#1142). The transcript
/// fixture does set a gene symbol (`GENE1099`), so these assertions would catch
/// it leaking back into either renderer.
#[test]
fn gene_symbol_rendering_matches_the_changed_residue_sibling() {
    let vp = projector_for(ISSUE_CDS);
    // Separated: bracket, selector suppressed on both.
    assert_eq!(
        protein_of(&vp, "REF:g.[6G>A;15C>A]"),
        "txp:p.[(Leu2=);(Ala5=)]"
    );
    assert_eq!(
        protein_of(&vp, "REF:g.[5T>A;14C>A]"),
        "txp:p.[(Leu2Gln);(Ala5Asp)]"
    );
    // Consecutive: single variant, selector suppressed on both.
    assert_eq!(protein_of(&vp, "REF:g.[6G>A;9T>C]"), "txp:p.(Leu2_Arg3=)");
    assert_eq!(
        protein_of(&vp, "REF:g.[5T>A;8G>A]"),
        "txp:p.(Leu2_Arg3delinsGlnHis)"
    );
}

/// Members inside one codon still collapse to that single residue — the fix
/// must not turn a one-codon allele into a bracket or a range.
#[test]
fn all_silent_members_in_one_codon_stay_a_single_residue() {
    let vp = projector_for(ISSUE_CDS);
    assert_eq!(protein_of(&vp, "REF:g.[4C>T;6G>A]"), "txp:p.(Leu2=)");
}

/// A single silent substitution is unchanged by this fix.
#[test]
fn single_silent_substitution_is_unchanged() {
    let vp = projector_for(ISSUE_CDS);
    assert_eq!(protein_of(&vp, "REF:g.9T>C"), "txp:p.(Arg3=)");
    assert_eq!(protein_of(&vp, "REF:g.6G>A"), "txp:p.(Leu2=)");
}

/// A silent `delins` spanning two codons must name them, not emit the
/// whole-molecule `p.(=)`. This is reachable from a **single variant** — no
/// allele involved — and is the half of the defect the issue does not report.
#[test]
fn silent_delins_names_its_codons_instead_of_whole_protein_identity() {
    let vp = projector_for(ISSUE_CDS);
    // c.6_9 -> ACGC: codon 2 CTG->CTA (Leu), codon 3 CGT->CGC (Arg). Both
    // triplets rewritten, both silent.
    let rendered = protein_of(&vp, "REF:g.6_9delinsACGC");
    assert_eq!(rendered, "txp:p.(Leu2_Arg3=)");
    assert!(
        rendered != "txp:p.(=)",
        "a targeted delins must not claim whole-protein identity"
    );
}

/// ACCEPTED LIMITATION, pinned so it stays deliberate: a codon *interior* to a
/// rewritten span whose triplet is coincidentally unchanged is not named, so the
/// run is reported as two separated codons rather than one range.
///
/// `c.4_12delinsCTACGTCTT` rewrites codon 2 (`CTG`→`CTA`, Leu) and codon 4
/// (`CTC`→`CTT`, Leu) but leaves codon 3's `CGT` byte-identical. The output
/// `p.[(Leu2=);(Leu4=)]` is true but less complete than `p.(Leu2_Leu4=)`.
/// Closing the gap would require deriving the set from the edit's *span*, which
/// is exactly what a 3'-shifted spelling makes unsound — so the limitation is
/// the price of the shift-invariance the rest of these tests depend on.
#[test]
fn an_unchanged_interior_codon_is_not_named() {
    let vp = projector_for(RUN_CDS);
    assert_eq!(
        protein_of(&vp, "REF:g.4_12delinsCTACGTCTT"),
        "txp:p.[(Leu2=);(Leu4=)]"
    );
}

/// Order and grouping of members must not affect the result — the defect the
/// issue describes is precisely that they did.
#[test]
fn output_is_independent_of_member_order() {
    let vp = projector_for(ISSUE_CDS);
    let forward = protein_of(&vp, "REF:g.[6G>A;9T>C]");
    let reverse = protein_of(&vp, "REF:g.[9T>C;6G>A]");
    assert_eq!(forward, reverse);
}

/// Every emitted shape survives `ferro normalize` byte-for-byte.
///
/// The two new shapes — a range identity and a bracketed cis allele of
/// identities — are protein cis runs flowing through the same merge machinery
/// that coalesces `p.[(Arg76Ser);(Cys77Trp)]` into a delins (#1127/#1128/#1133/
/// #1134). An `=` is not a change, so nothing may be coalesced here; this pins
/// that the coalescer leaves them alone rather than inventing a delins whose
/// inserted residues equal the reference (which `protein/delins.md:55` marks
/// invalid, "effectively is no change").
#[test]
fn emitted_identities_survive_normalization_unchanged() {
    let normalizer = Normalizer::new(MockProvider::new());
    for rendered in [
        "txp:p.(Leu2=)",
        "txp:p.(Leu2_Arg3=)",
        "txp:p.(Leu2_Leu4=)",
        "txp:p.[(Leu2=);(Ala5=)]",
        "txp:p.(=)",
    ] {
        let parsed =
            parse_hgvs(rendered).unwrap_or_else(|e| panic!("{rendered} does not parse: {e}"));
        let normalized = normalizer
            .normalize(&parsed)
            .unwrap_or_else(|e| panic!("{rendered} does not normalize: {e}"));
        assert_eq!(normalized.to_string(), rendered, "normalization moved it");
    }
}

/// Emitted consequences have to be real HGVS.
#[test]
fn emitted_identities_reparse() {
    let vp = projector_for(ISSUE_CDS);
    for input in ["REF:g.[4C>T;6G>A;9T>C]", "REF:g.[6G>A;15C>A]", "REF:g.9T>C"] {
        let rendered = protein_of(&vp, input);
        ferro_hgvs::parse_hgvs(&rendered)
            .unwrap_or_else(|e| panic!("{input} -> {rendered} does not re-parse: {e}"));
    }
}
