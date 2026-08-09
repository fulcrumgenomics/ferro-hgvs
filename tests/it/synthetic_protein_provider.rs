//! The synthetic protein reference provider and the amino-acid denotation
//! oracle (`src/conformance/synthetic_protein.rs`).
//!
//! # What is being judged, and against what
//!
//! Two different things, deliberately kept apart:
//!
//! - **The provider** is judged against itself, structurally: what it serves as
//!   `NP_TEST.1` must be what `NM_TEST.1`'s CDS translates to. That is a
//!   closed property with no external authority, so a round trip is the right
//!   test.
//! - **The oracle** is judged against the **spec**, not against the applier
//!   that implements it. Every expectation below that names a protein sequence
//!   is a string the HGVS recommendations publish, quoted verbatim in the
//!   test's doc comment with an exact `file.md:line`. Re-deriving an
//!   expectation by running the applier would make these change detectors
//!   rather than conformance tests — the distinction this repository's
//!   `CLAUDE.md` draws when it says a test that merely pins today's output is
//!   not an adjudication record.
//!
//! # The one substitution made to the spec's examples
//!
//! The **accession** is rewritten to `NP_TEST.1`. The spec's examples name
//! `LRG_199p1`, `NP_003997.1`, `NP_004371.2` and others, whose sequences are
//! not distributed with the recommendations; positions, residues and edit
//! spellings are preserved exactly. Where an example names a position the
//! reference must carry a specific residue at, the synthetic peptide is built
//! to carry it — see [`bench`].
//!
//! Line numbers are against the pinned `assets/hgvs-nomenclature` submodule
//! (`6f85311`). Paths are relative to `assets/hgvs-nomenclature/docs/`.

use ferro_hgvs::backtranslate::codon::{Codon, CodonTable};
use ferro_hgvs::conformance::synthetic_protein::{
    protein_denotation_of, Indeterminacy, Indeterminate, NoSingleSequence, ProteinDenotation,
    ProteinFrame, ProteinFrameError, ProteinRefShape, PROTEIN_ACCESSION,
};
use ferro_hgvs::hgvs::location::AminoAcid;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::reference::ReferenceProvider;

// ---------------------------------------------------------------------------
// Peptide construction
// ---------------------------------------------------------------------------

/// A peptide of `len` residues with `fixed` residues placed at the 1-based
/// positions given, `Met` at position 1, and a non-repeating filler elsewhere.
///
/// The filler cycles through eight residues rather than repeating one, so no
/// single-residue run and no short tandem repeat exists by accident. That
/// matters twice over: a deletion in a homopolymer run is ambiguous by design
/// (`protein/deletion.md:39`), and a repeat description keys on a tract the
/// reference actually holds (`protein/repeated.md:20-23`). A filler that
/// created either would make a test pass or fail for a reason it does not name.
fn bench(len: usize, fixed: &[(usize, AminoAcid)]) -> Vec<AminoAcid> {
    const CYCLE: [AminoAcid; 8] = [
        AminoAcid::Ala,
        AminoAcid::Gly,
        AminoAcid::Ser,
        AminoAcid::Thr,
        AminoAcid::Val,
        AminoAcid::Leu,
        AminoAcid::Ile,
        AminoAcid::Pro,
    ];
    let mut peptide: Vec<AminoAcid> = (0..len).map(|index| CYCLE[index % CYCLE.len()]).collect();
    peptide[0] = AminoAcid::Met;
    for (position, aa) in fixed {
        peptide[*position - 1] = *aa;
    }
    peptide
}

/// Parse a peptide written the way the spec writes one, e.g. `MetLeuTrpTrpGlu`.
///
/// # Panics
///
/// If `spelled` is not a whole number of codes from the table at
/// `background/standards.md:216`-`:242`.
fn peptide(spelled: &str) -> Vec<AminoAcid> {
    ProteinFrame::peptide_from_three_letter(spelled)
        .unwrap_or_else(|| panic!("{spelled} is not a three-letter peptide"))
}

/// A single-exon frame over `peptide`.
///
/// # Panics
///
/// If the peptide cannot be back-translated and round-tripped.
fn frame(peptide: &[AminoAcid]) -> ProteinFrame {
    ProteinFrame::from_peptide(peptide).expect("a Met-initiated peptide builds a frame")
}

/// The denotation of `NP_TEST.1:p.<suffix>` on `frame`.
fn denote(frame: &ProteinFrame, suffix: &str) -> ProteinDenotation {
    protein_denotation_of(frame.provider(), &ProteinFrame::protein_descriptor(suffix))
}

/// The denoted sequence, or a panic naming what came back instead.
///
/// # Panics
///
/// If the description does not denote a single sequence.
fn sequence(frame: &ProteinFrame, suffix: &str) -> String {
    match denote(frame, suffix) {
        ProteinDenotation::Sequence(spelled) => spelled,
        other => panic!("p.{suffix} denoted {other:?}, expected a single sequence"),
    }
}

// ---------------------------------------------------------------------------
// The provider: the protein is the transcript's own translation
// ---------------------------------------------------------------------------

/// The served protein is the translation of the served transcript's CDS.
///
/// This is the whole point of the frame. The check reads both sides back out of
/// the provider — the protein through `get_protein_sequence`, the CDS through
/// the transcript record — so it cannot pass by comparing the builder's inputs
/// to themselves.
#[test]
fn the_served_protein_is_the_served_transcripts_translation() {
    for shape in ProteinRefShape::all() {
        let requested = bench(40, &[(24, AminoAcid::Trp)]);
        let built = ProteinFrame::build(shape, &requested).expect("frame builds");

        let transcript = built
            .provider()
            .get_transcript("NM_TEST.1")
            .expect("the frame serves its transcript");
        let bases = transcript
            .sequence
            .as_deref()
            .expect("the frame's transcript carries bases");
        let (cds_start, cds_end) = built.cds_bounds();
        assert_eq!(
            &bases[cds_start - 1..cds_end],
            built.cds(),
            "{}: the CDS the provider serves is the CDS the frame reports",
            shape.label()
        );

        let served = built
            .provider()
            .get_protein_sequence(PROTEIN_ACCESSION, 0, requested.len() as u64)
            .expect("the frame serves its protein");
        let expected: String = requested.iter().map(AminoAcid::to_one_letter).collect();
        assert_eq!(
            served,
            expected,
            "{}: the served protein is the requested peptide",
            shape.label()
        );

        // And the CDS really does encode it. The length check alone does not say
        // that — it says the CDS is the right SIZE for "one codon per residue
        // plus a terminator", which a CDS whose last codon coded an amino acid
        // would satisfy just as well. This comment claimed the terminator was
        // checked while only the arithmetic was, so the claim is now the
        // assertion.
        //
        // The per-residue half of the claim is genuinely elsewhere
        // (`a_residues_codon_is_where_the_cds_says_it_is`), and is cross-linked
        // rather than duplicated.
        assert_eq!(
            built.cds().len(),
            (requested.len() + 1) * 3,
            "{}: the CDS is one codon per residue plus a terminator",
            shape.label()
        );
        let cds = built.cds();
        let table = CodonTable::standard();
        let last = Codon::parse(&cds[cds.len() - 3..]).expect("the final CDS codon parses");
        assert!(
            table.is_stop(&last),
            "{}: the CDS's final codon must be a terminator, not merely leave room for one",
            shape.label()
        );
        // `is_stop`, not `amino_acid_for(..) == Some(Ter)`: the table models a
        // terminator as the ABSENCE of an amino acid, so the latter is `None`
        // here and reads as "this codon is unknown". Pinned alongside, because
        // that asymmetry is what makes the natural spelling of this assertion
        // the wrong one.
        assert_eq!(
            table.amino_acid_for(&last),
            None,
            "{}: a terminator codes no amino acid",
            shape.label()
        );
    }
}

/// The transcript record names the protein, so a consumer resolving the protein
/// *through* the transcript reaches the same accession.
#[test]
fn the_transcript_records_its_protein_accession() {
    let built = frame(&bench(30, &[]));
    let transcript = built
        .provider()
        .get_transcript("NM_TEST.1")
        .expect("the frame serves its transcript");
    assert_eq!(
        transcript.protein_id.as_deref(),
        Some(PROTEIN_ACCESSION),
        "the transcript points at the protein the frame serves"
    );
}

/// Every geometry serves the same protein: exon structure and strand change the
/// contig, never the translation.
///
/// The multi-exon shapes are the ones #1478 says a corpus must be able to vary;
/// this pins that varying them does not silently vary the molecule as well.
#[test]
fn every_transcript_geometry_serves_one_protein() {
    let requested = bench(45, &[(24, AminoAcid::Trp), (30, AminoAcid::Cys)]);
    let proteins: Vec<String> = ProteinRefShape::all()
        .into_iter()
        .map(|shape| {
            ProteinFrame::build(shape, &requested)
                .expect("frame builds")
                .spelled()
        })
        .collect();
    assert_eq!(proteins.len(), 3, "single-exon plus both strands");
    assert!(
        proteins.windows(2).all(|pair| pair[0] == pair[1]),
        "geometries disagreed about the protein: {proteins:?}"
    );
    assert_eq!(proteins[0], ProteinFrame::spell(&requested));
}

/// `cds_position_of_residue` is the seam a paired `c.`/`p.` stratum needs, and
/// it must agree with the CDS the frame actually serves.
#[test]
fn a_residues_codon_is_where_the_cds_says_it_is() {
    let requested = peptide("MetLeuTrpTrpGlu");
    let built = frame(&requested);
    let cds = built.cds();
    let table = CodonTable::standard();
    for (index, aa) in requested.iter().enumerate() {
        let residue = index as u64 + 1;
        let position = built
            .cds_position_of_residue(residue)
            .expect("an in-range residue has a codon");
        let start = position as usize - 1;
        let codon = Codon::parse(&cds[start..start + 3]).expect("a CDS codon parses");
        assert_eq!(
            table.amino_acid_for(&codon),
            Some(*aa),
            "the codon `c.{position}` names encodes residue {residue}"
        );
    }
    // The terminator has a codon too, and nothing past it does.
    assert!(built
        .cds_position_of_residue(requested.len() as u64 + 1)
        .is_some());
    assert!(built
        .cds_position_of_residue(requested.len() as u64 + 2)
        .is_none());
}

/// A peptide the builder cannot honestly realise is refused, not approximated.
#[test]
fn an_unrealisable_peptide_is_refused() {
    assert_eq!(
        ProteinFrame::from_peptide(&[]).err(),
        Some(ProteinFrameError::Empty)
    );

    // `background/standards.md:229` lists `ATG` as "(translation initiation)",
    // and an authoritative `NP_` therefore always begins with `Met`.
    assert_eq!(
        ProteinFrame::from_peptide(&[AminoAcid::Ala, AminoAcid::Gly]).err(),
        Some(ProteinFrameError::NotMetInitiated(AminoAcid::Ala))
    );

    // A `Ter` inside the peptide would truncate the translation, so the round
    // trip could not hold.
    assert_eq!(
        ProteinFrame::from_peptide(&[AminoAcid::Met, AminoAcid::Ter, AminoAcid::Gly]).err(),
        Some(ProteinFrameError::InternalStop(2))
    );

    // `background/standards.md:239` gives `Xaa` the codons `NNN` — it names no
    // codon the standard table can spell.
    assert_eq!(
        ProteinFrame::from_peptide(&[AminoAcid::Met, AminoAcid::Xaa]).err(),
        Some(ProteinFrameError::NotBackTranslatable(AminoAcid::Xaa))
    );
}

/// Three-letter spelling round-trips, and rejects a string that is not one.
#[test]
fn three_letter_spelling_round_trips() {
    let residues = peptide("MetLeuTrpTrpGlu");
    assert_eq!(ProteinFrame::spell(&residues), "MetLeuTrpTrpGlu");
    assert_eq!(ProteinFrame::peptide_from_three_letter("MetLeuX"), None);
    assert_eq!(ProteinFrame::peptide_from_three_letter("MetZzz"), None);
    // `Ter` and `Xaa` both spell, which is why the denoted form is three-letter:
    // neither has an unambiguous one-letter code in a sequence that may also
    // hold a stop (`background/standards.md:213`, `:244`).
    assert_eq!(
        ProteinFrame::spell(&[AminoAcid::Ter, AminoAcid::Xaa]),
        "TerXaa"
    );
}

/// No denoted sequence ever spells a stop as `X`.
///
/// `background/standards.md:213`: "HGVS nomenclature uses `Ter` (three-letter
/// amino acid code) and `*` (three- and one-letter amino acid code) to indicate
/// a translation termination (stop) codon (in older versions, `X` was used
/// instead)."
#[test]
fn a_denoted_sequence_never_spells_a_stop_as_x() {
    let built = frame(&bench(30, &[(24, AminoAcid::Trp)]));
    let nonsense = sequence(&built, "Trp24Ter");
    assert!(nonsense.ends_with("Ter"), "{nonsense}");
    assert!(
        !nonsense.contains('X'),
        "a denoted sequence must not contain X: {nonsense}"
    );
    assert!(
        !nonsense.contains('*'),
        "the three-letter form spells the stop `Ter`: {nonsense}"
    );
}

// ---------------------------------------------------------------------------
// The oracle, against the spec's own worked examples
// ---------------------------------------------------------------------------

/// `protein/deletion.md:37-38`:
///
/// > **`NP_003997.2:p.Trp4del`**
/// > a deletion of amino acid `Trp4` in the sequence `MetLeuTrp`Trp`Glu` to
/// > `MetLeuTrpGlu`.
///
/// The spec gives both sides verbatim, which is what makes this an oracle test
/// and not a change detector. `:39` adds the ambiguity the reference itself
/// carries — "for deletions in single amino acid stretches or tandem repeats,
/// the most C-terminal residue is arbitrarily assigned to have been deleted" —
/// so `p.Trp3del` denotes the same protein and is checked alongside.
#[test]
fn spec_deletion_of_one_residue_in_a_two_residue_stretch() {
    let built = frame(&peptide("MetLeuTrpTrpGlu"));
    assert_eq!(sequence(&built, "Trp4del"), "MetLeuTrpGluTer");
    assert_eq!(
        sequence(&built, "Trp3del"),
        sequence(&built, "Trp4del"),
        "the two placements of one deletion in a `Trp` stretch denote one protein"
    );
}

/// `protein/deletion.md:106`:
///
/// > E.g., when `MetTrpSerSer`Ser`HisAsp..` changes to `MetTrpSerSerHisAsp..`,
/// > this is described as `p.Ser5del`.
#[test]
fn spec_deletion_derived_from_the_protein_not_the_dna() {
    let built = frame(&peptide("MetTrpSerSerSerHisAsp"));
    assert_eq!(sequence(&built, "Ser5del"), "MetTrpSerSerHisAspTer");
    // `:107` says the DNA-derived `p.Ser4del` "is not correct" — but it denotes
    // the same protein, which is exactly why the spec has to rule on spelling.
    assert_eq!(sequence(&built, "Ser4del"), sequence(&built, "Ser5del"));
}

/// `protein/duplication.md:38-39`:
///
/// > **`NP_003997.2:p.Trp4dup`**
/// > a duplication of amino acid `Trp4` in the sequence `MetLeuTrpTrpGlu` to
/// > `MetLeuTrpTrp`Trp`Glu`.
#[test]
fn spec_duplication_of_one_residue() {
    let built = frame(&peptide("MetLeuTrpTrpGlu"));
    assert_eq!(sequence(&built, "Trp4dup"), "MetLeuTrpTrpTrpGluTer");
}

/// `protein/duplication.md:97`:
///
/// > E.g., when `MetTrpSerSerSerHisAsp..` changes to `MetTrpSerSerSer`Ser`HisAsp..`,
/// > this is described as `p.Ser5dup`.
#[test]
fn spec_duplication_derived_from_the_protein_not_the_dna() {
    let built = frame(&peptide("MetTrpSerSerSerHisAsp"));
    assert_eq!(sequence(&built, "Ser5dup"), "MetTrpSerSerSerSerHisAspTer");
}

/// `protein/insertion.md:31-32`:
///
/// > **`p.His4_Gln5insAla`**
/// > the insertion of amino acid Ala between amino acids `His4` and `Gln5`,
/// > changing `MetLysGlyHisGlnGlnCys` to `MetLysGlyHis`Ala`GlnGlnCys`.
///
/// and `:34-35`:
///
/// > **`p.Lys2_Gly3insGlnSerLys`**
/// > the insertion of amino acids GlnSerLys between amino acids `Lys2` and
/// > `Gly3`, changing `MetLysGlyHisGlnGlnCys` to `MetLys`GlnSerLys`GlyHisGlnGlnCys`.
#[test]
fn spec_insertions_between_two_flanking_residues() {
    let built = frame(&peptide("MetLysGlyHisGlnGlnCys"));
    assert_eq!(
        sequence(&built, "His4_Gln5insAla"),
        "MetLysGlyHisAlaGlnGlnCysTer"
    );
    assert_eq!(
        sequence(&built, "Lys2_Gly3insGlnSerLys"),
        "MetLysGlnSerLysGlyHisGlnGlnCysTer"
    );
}

/// `protein/insertion.md:37-39`:
///
/// > **`p.(Met3_His4insGlyTer)`**
/// > the predicted consequence on the protein level of an insertion on the DNA
/// > level (`c.9_10insGGGTAG`), is the insertion of `GlyTer` (alternatively
/// > `Gly*`).
/// > **NOTE**: this is not described as `p.(Met3_Ile3418delinsGly)`, a
/// > deletion-insertion replacing the entire C-terminal protein coding sequence
/// > downstream of `Met3` with a `Gly`.
///
/// The note is the interesting half: the two descriptions denote the **same**
/// protein, so the spec is ruling on spelling. Both are checked, because an
/// oracle that disagreed with the note would be measuring the wrong thing.
#[test]
fn spec_insertion_encoding_a_stop_truncates_the_protein() {
    let built = frame(&peptide("MetLysMetHisGlnGlnCys"));
    assert_eq!(sequence(&built, "(Met3_His4insGlyTer)"), "MetLysMetGlyTer");
    assert_eq!(
        sequence(&built, "Met3_Cys7delinsMetGly"),
        sequence(&built, "(Met3_His4insGlyTer)"),
        "the delins the note rejects denotes the same protein, which is why the \
         spec has to rule on the spelling rather than on the denotation"
    );
}

/// `protein/insertion.md:17-18` — an insertion's positions must flank one
/// junction:
///
/// > the "positions_flanking" should contain **two flanking residues**, e.g.,
/// > `Lys23_Leu24`, not two non-flanking residues (`Lys23_Asn25`).
#[test]
fn spec_non_flanking_insertion_denotes_nothing() {
    let built = frame(&bench(30, &[(23, AminoAcid::Lys), (25, AminoAcid::Asn)]));
    assert_eq!(
        denote(&built, "Lys23_Asn25insAsp"),
        ProteinDenotation::NoSingleSequence(NoSingleSequence::NonFlankingInsertion)
    );
}

/// `protein/delins.md:37-38` and `:40-41`:
///
/// > **`p.Cys28delinsTrpVal`**
/// > a deletion of amino acid `Cys28`, replaced by `TrpVal`.
///
/// > **`p.Cys28_Lys29delinsTrp`**
/// > a deletion of amino acids `Cys28` and `Lys29`, replaced by `Trp`.
#[test]
fn spec_delins_replaces_the_named_residues() {
    let built = frame(&bench(32, &[(28, AminoAcid::Cys), (29, AminoAcid::Lys)]));
    let reference = built.spelled();

    let widened = sequence(&built, "Cys28delinsTrpVal");
    assert_eq!(
        widened.len(),
        reference.len() + 3 /* one extra code */ + 3 /* Ter */
    );
    assert_eq!(&widened[..27 * 3], &reference[..27 * 3]);
    assert_eq!(&widened[27 * 3..29 * 3], "TrpVal");
    assert_eq!(&widened[29 * 3..], &format!("{}Ter", &reference[28 * 3..]));

    let narrowed = sequence(&built, "Cys28_Lys29delinsTrp");
    assert_eq!(&narrowed[27 * 3..28 * 3], "Trp");
    assert_eq!(&narrowed[..27 * 3], &reference[..27 * 3]);
    assert_eq!(&narrowed[28 * 3..], &format!("{}Ter", &reference[29 * 3..]));
}

/// `protein/delins.md:62-64`:
///
/// > **`p.[Ser44Arg;Trp46Arg]`**
/// > the change of two variants affecting amino acids separated by another
/// > amino acid.
/// > **NOTE**: the variant is not described as `p.Ser44_Trp46delinsArgLeuArg`.
///
/// The note names a second spelling of the **same** variant, so the two must
/// denote one protein. That is the cis-allele case the oracle's multi-member
/// walk exists for, and the only place in the protein recommendations where a
/// two-member allele and its spanning delins are published side by side.
#[test]
fn spec_cis_allele_and_its_rejected_spanning_delins_denote_one_protein() {
    let built = frame(&bench(
        50,
        &[
            (44, AminoAcid::Ser),
            (45, AminoAcid::Leu),
            (46, AminoAcid::Trp),
        ],
    ));
    let members = sequence(&built, "[Ser44Arg;Trp46Arg]");
    let spanning = sequence(&built, "Ser44_Trp46delinsArgLeuArg");
    assert_eq!(members, spanning);
    assert_eq!(&members[43 * 3..46 * 3], "ArgLeuArg");
}

/// `protein/delins.md:53-56`:
///
/// > **`NP_000213.1:p.(Val559_Glu561del)`**
/// > … The variant is **not** described as `p.(Val559_Glu562delinsGlu)`, where
/// > `Glu562` would be replaced by a `Glu`, which effectively is no change.
///
/// Again a spelling ruling over two descriptions of one protein, so the oracle
/// must agree they denote the same thing.
#[test]
fn spec_delins_with_a_no_change_tail_denotes_the_deletion() {
    let peptide = bench(
        570,
        &[
            (559, AminoAcid::Val),
            (561, AminoAcid::Glu),
            (562, AminoAcid::Glu),
        ],
    );
    let built = frame(&peptide);

    // The absolute anchor first. Both sides of the equivalence below are oracle
    // output, so on their own they stay equal under a defect that moved BOTH
    // spans by the same amount — an off-by-one in `range_len`, or a wrong
    // `start` out of `resolve`. Pin what the protein actually is, built here
    // from the reference peptide rather than read back out of the oracle.
    let mut expected = peptide.clone();
    expected.drain(558..561); // 1-based 559..=561, the three residues deleted
    let deletion = sequence(&built, "(Val559_Glu561del)");
    assert_eq!(
        deletion,
        format!("{}Ter", ProteinFrame::spell(&expected)),
        "the deletion must remove exactly residues 559..=561, leaving the reference's \
         `Glu562` as the new 559"
    );

    assert_eq!(deletion, sequence(&built, "(Val559_Glu562delinsGlu)"));
}

/// `protein/substitution.md:29-30` and `:32-33`:
///
/// > **`LRG_199p1:p.Trp24Cys`**
/// > amino acid `Trp24` is changed to a `Cys`.
///
/// > **`NP_003997.1:p.(Trp24Cys)`**
/// > amino acid `Trp24` is **predicted** to change to a `Cys` (no experimental
/// > proof, e.g., based on DNA level data).
///
/// The parenthesised form is the same change with weaker evidence
/// (`:16`: "predicted consequences, i.e. without experimental evidence … should
/// be given in parentheses"), so it denotes the same single protein. Treating
/// the wrapper as ambiguity would drop most of the spec's own examples from any
/// future census, on a misreading of what the parentheses record.
#[test]
fn spec_missense_substitution_and_its_predicted_spelling() {
    let built = frame(&bench(30, &[(24, AminoAcid::Trp)]));
    let observed = sequence(&built, "Trp24Cys");
    assert_eq!(&observed[23 * 3..24 * 3], "Cys");
    assert_eq!(sequence(&built, "(Trp24Cys)"), observed);
}

/// `protein/substitution.md:36-38`:
///
/// > **`LRG_199p1:p.Trp24Ter` (`p.Trp24*`)**
/// > amino acid `Trp24` is changed to a stop codon (`Ter`, `*`).
/// > **NOTE**: this change is **not** described as a deletion of the C-terminal
/// > end of the protein (i.e. `p.Trp24_Met36853del`).
///
/// Both spellings of the stop must denote one protein — that is the whole
/// content of "(`Ter`, `*`)".
///
/// The note is the same kind of ruling as `protein/delins.md:27` and
/// `protein/insertion.md:39`: the C-terminal deletion it rejects denotes
/// **exactly the same protein**, because the deletion stops short of the
/// terminator and the terminator survives. That is why the recommendation has
/// to be stated at all — nothing about the denoted molecule distinguishes the
/// two spellings, so the choice cannot be made by an oracle and has to be made
/// by the nomenclature. The last assertion pins that, rather than asserting a
/// difference the reference does not have.
#[test]
fn spec_nonsense_substitution_terminates_translation() {
    let built = frame(&bench(30, &[(24, AminoAcid::Trp), (30, AminoAcid::Pro)]));
    let truncated = sequence(&built, "Trp24Ter");
    assert_eq!(truncated, sequence(&built, "Trp24*"));
    assert_eq!(
        truncated.len(),
        24 * 3,
        "23 residues plus the new terminator"
    );
    assert!(truncated.ends_with("Ter"));

    assert_eq!(
        sequence(&built, "Trp24_Pro30del"),
        truncated,
        "the C-terminal deletion the note rejects denotes the same protein"
    );
}

/// `protein/substitution.md:41-43`:
///
/// > **`NP_003997.1:p.Cys188=`**
/// > amino acid `Cys188` is not changed …
/// > **NOTE**: the description `p.=` means the **entire** protein coding region
/// > was analysed and no variant was found that changes (or is predicted to
/// > change) the protein sequence.
#[test]
fn spec_silent_descriptions_denote_the_reference() {
    let built = frame(&bench(200, &[(188, AminoAcid::Cys)]));
    let reference = format!("{}Ter", built.spelled());
    assert_eq!(sequence(&built, "Cys188="), reference);
    assert_eq!(sequence(&built, "="), reference);
    assert_eq!(sequence(&built, "(=)"), reference);
}

/// `protein/substitution.md:46-48`:
///
/// > **no protein: `LRG_199p1:p.0`**
/// > as a consequence of a variant in the translation initiation codon, no
/// > protein is produced.
/// > **NOTE**: `LRG_199p1:p.0?` can be used when you predict that no protein is
/// > produced.
///
/// "No protein" is not "a protein of length zero", so it gets its own answer
/// rather than an empty sequence.
#[test]
fn spec_no_protein_is_not_an_empty_sequence() {
    let built = frame(&bench(30, &[]));
    assert_eq!(
        denote(&built, "0"),
        ProteinDenotation::NoProtein { predicted: false }
    );
    assert_eq!(
        denote(&built, "0?"),
        ProteinDenotation::NoProtein { predicted: true }
    );
}

/// `protein/substitution.md:51-52` and `:71-73`:
///
/// > **unknown: `LRG_199p1:p.(Met1?)`**
/// > the consequence, on the protein level, of a variant affecting the
/// > translation initiation codon can not be predicted (i.e. is unknown).
///
/// > **`NP_003997.1:p.?`**
/// > … it can not be excluded that the variant affects splicing, having unknown
/// > consequences.
#[test]
fn spec_unknown_consequences_denote_no_single_protein() {
    let built = frame(&bench(30, &[]));
    for suffix in ["?", "(?)", "(Met1?)", "Met1?"] {
        assert_eq!(
            denote(&built, suffix),
            ProteinDenotation::NoSingleSequence(NoSingleSequence::UnknownConsequence),
            "p.{suffix}"
        );
    }
}

/// `protein/substitution.md:77-78` and `:81-82`:
///
/// > **`NP_003997.1:p.(Gly56Ala^Ser^Cys)`**
/// > amino acid `Gly56` is changed to an `Ala`, `Ser`, or `Cys`.
///
/// > **`LRG_199p1:p.Trp24=/Cys`**
/// > a mosaic case where at amino acid position `24`, besides the normal amino
/// > acid (a `Trp`, described as `=`), also protein is found containing a `Cys`.
///
/// Both name a set of proteins. Returning one of them would make a future
/// census wrong in the flattering direction, so both are refused as set-valued.
#[test]
fn spec_set_valued_descriptions_are_refused() {
    let built = frame(&bench(60, &[(24, AminoAcid::Trp), (56, AminoAcid::Gly)]));
    for suffix in ["(Gly56Ala^Ser^Cys)", "Trp24=/Cys", "(Trp24=/Cys)"] {
        assert_eq!(
            denote(&built, suffix),
            ProteinDenotation::NoSingleSequence(NoSingleSequence::SetValued),
            "p.{suffix}"
        );
    }
}

/// `protein/frameshift.md:33-34`:
///
/// > **`p.Arg97ProfsTer23` (short `p.Arg97fs`) / `p.Arg97Profs*23`**
/// > a variant with `Arg97` as the first amino acid changed, shifting the
/// > reading frame, replacing it for a `Pro` and terminating at position
/// > `Ter23`.
///
/// and `:20`:
///
/// > the position of the translation termination (stop) codon in the new
/// > reading frame is calculated **starting** at the first amino acid changed by
/// > the frameshift (codon 1), and **ending** at the first stop codon
/// > encountered (`Ter#` or `*#`).
///
/// The residues between the first new one and the new stop are read in a
/// different frame, which a protein reference does not carry. So the oracle
/// reports what the description fixes — the unchanged prefix, the first new
/// residue, and the total length — and nothing else. Inventing the tail would
/// be a guess.
#[test]
fn spec_frameshift_determines_a_prefix_and_a_length_only() {
    let built = frame(&bench(140, &[(97, AminoAcid::Arg)]));
    let prefix = format!("{}Pro", &built.spelled()[..96 * 3]);

    let expected = ProteinDenotation::Indeterminate(Indeterminate {
        prefix: prefix.clone(),
        length: Some(96 + 23),
        reason: Indeterminacy::FrameshiftTail,
    });
    assert_eq!(denote(&built, "Arg97ProfsTer23"), expected);
    assert_eq!(denote(&built, "Arg97Profs*23"), expected);

    // The short format states no termination detail (`:23`), so no length.
    assert_eq!(
        denote(&built, "Arg97fs"),
        ProteinDenotation::Indeterminate(Indeterminate {
            prefix: built.spelled()[..96 * 3].to_string(),
            length: None,
            reason: Indeterminacy::FrameshiftTail,
        })
    );
}

/// `protein/frameshift.md:44-45`:
///
/// > **`p.Ile327Argfs*?` (short `p.Ile327fs`)**
/// > the predicted consequence of a frameshifting variant changes `Ile327` to an
/// > `Arg`, but the new reading frame does not encounter a new translation
/// > termination (stop) codon.
#[test]
fn spec_frameshift_reaching_no_stop_has_no_length() {
    let built = frame(&bench(400, &[(327, AminoAcid::Ile)]));
    assert_eq!(
        denote(&built, "Ile327Argfs*?"),
        ProteinDenotation::Indeterminate(Indeterminate {
            prefix: format!("{}Arg", &built.spelled()[..326 * 3]),
            length: None,
            reason: Indeterminacy::FrameshiftTail,
        })
    );
}

/// `protein/extension.md:30-31`:
///
/// > **`p.Ter110GlnextTer17` (alternatively `p.*110Glnext*17`)**
/// > a variant in the stop codon (`Ter` / `*`) at position 110, changing it to a
/// > `Gln`-codon (a no-stop variant) and adding a tail of new amino acids to the
/// > protein's C-terminus, ending at a new stop codon (`Ter` / `*`) at position
/// > 17 of the added sequence.
///
/// The added tail comes from the 3'UTR, which a protein reference does not
/// carry, so — as for a frameshift — the prefix and the length are reported and
/// the tail is not invented. Note the description addresses position 110 on a
/// **109-residue** protein: the terminator is addressable, which is why the
/// oracle's reference carries it.
#[test]
fn spec_c_terminal_extension_determines_a_prefix_and_a_length() {
    let built = frame(&bench(109, &[]));
    let expected = ProteinDenotation::Indeterminate(Indeterminate {
        prefix: format!("{}Gln", built.spelled()),
        length: Some(110 + 17),
        reason: Indeterminacy::CTerminalExtension,
    });
    assert_eq!(denote(&built, "Ter110GlnextTer17"), expected);
    assert_eq!(denote(&built, "*110Glnext*17"), expected);
}

/// `protein/extension.md:33-34`:
///
/// > **`p.Ter327ArgextTer?` (alternatively `p.*327Argext*?`)**
/// > a variant in the stop codon (`Ter` / `*`) at position 327, changing it to
/// > an `Arg`-codon and adding a tail of new amino acids of unknown length
/// > (position `Ter?`) since the shifted frame does not contain a new stop
/// > codon.
#[test]
fn spec_c_terminal_extension_reaching_no_stop_has_no_length() {
    let built = frame(&bench(326, &[]));
    assert_eq!(
        denote(&built, "Ter327ArgextTer?"),
        ProteinDenotation::Indeterminate(Indeterminate {
            prefix: format!("{}Arg", built.spelled()),
            length: None,
            reason: Indeterminacy::CTerminalExtension,
        })
    );
}

/// `protein/extension.md:22-23`:
///
/// > **`p.Met1ext-5`**
/// > a variant in the 5' UTR activates a new upstream translation initiation
/// > site starting with amino acid `Met-5`.
///
/// The added residues are upstream of residue 1, so **no** N-terminal prefix is
/// determined even though the reference protein survives intact as a suffix.
/// Reporting an empty prefix is the honest answer; reporting the reference
/// would claim the description had fixed the new N-terminus.
#[test]
fn spec_n_terminal_extension_determines_no_prefix() {
    let built = frame(&bench(40, &[]));
    assert_eq!(
        denote(&built, "Met1ext-5"),
        ProteinDenotation::Indeterminate(Indeterminate {
            prefix: String::new(),
            // 40 residues + the terminator, plus the 5 added upstream.
            length: Some(41 + 5),
            reason: Indeterminacy::NTerminalExtension,
        })
    );
}

/// `protein/insertion.md:45-47`:
///
/// > **`p.Arg78_Gly79insXaa[23]`**
/// > the in-frame insertion of a 23 amino acid sequence between amino acids
/// > `Arg78` and `Gly79`.
/// > **NOTE**: it must be possible to deduce the 23 inserted amino acids from
/// > the description given on the DNA or RNA level.
///
/// and `:49-51`:
///
/// > **`NP_060250.2:p.Gln746_Lys747ins*63`**
/// > the in-frame insertion of a 62 amino acid sequence ending at a stop codon
/// > at position `*63` between amino acids `Gln746` and `Lys747`.
///
/// The note is explicit that the residues are recoverable only from the
/// nucleotide description, so the oracle reports the length and refuses to name
/// them. The two forms differ in one way that matters: `ins*63` terminates
/// **inside** the insertion, so the reference tail does not survive.
#[test]
fn spec_insertions_stated_by_count_determine_only_a_length() {
    let short = frame(&bench(90, &[(78, AminoAcid::Arg), (79, AminoAcid::Gly)]));
    assert_eq!(
        denote(&short, "Arg78_Gly79insXaa[23]"),
        ProteinDenotation::Indeterminate(Indeterminate {
            prefix: short.spelled()[..78 * 3].to_string(),
            // 90 residues + terminator + 23 inserted; the tail survives.
            length: Some(91 + 23),
            reason: Indeterminacy::ResiduesStatedByCount,
        })
    );

    let long = frame(&bench(760, &[(746, AminoAcid::Gln), (747, AminoAcid::Lys)]));
    assert_eq!(
        denote(&long, "Gln746_Lys747ins*63"),
        ProteinDenotation::Indeterminate(Indeterminate {
            prefix: long.spelled()[..746 * 3].to_string(),
            // 746 residues N-terminal of the junction, then a 63-position run
            // ending in the new stop.
            length: Some(746 + 63),
            reason: Indeterminacy::ResiduesStatedByCount,
        })
    );
}

/// `protein/repeated.md:20-23`:
///
/// > **`p.Ala2[10]`**
/// > a repeated amino acid sequence, with the first `Ala`-residue located at
/// > position 2, is present in 10 copies.
/// > **NOTE**: when the repeat is variable in the population and the reference
/// > sequence has 10 units, the description `p.Ala2[9]` is preferred over
/// > `p.Ala11del`.
/// > **NOTE**: … the description `p.Ala2[12]` is preferred over
/// > `p.Ala10_Ala11dup`.
///
/// Both notes name two descriptions of one variant, so both pairs must denote
/// one protein. That makes this a genuine oracle test rather than a check that
/// the repeat applier agrees with itself.
#[test]
fn spec_repeat_counts_agree_with_the_deletion_and_duplication_they_replace() {
    // Ten `Ala` at positions 2..11 on a reference that holds no other run.
    let mut residues = bench(30, &[]);
    for position in 2..=11 {
        residues[position - 1] = AminoAcid::Ala;
    }
    let built = frame(&residues);

    assert_eq!(
        sequence(&built, "Ala2[10]"),
        format!("{}Ter", built.spelled()),
        "the reference's own count denotes the reference"
    );
    // ABSOLUTE ANCHORS FIRST, for the reason
    // `spec_delins_with_a_no_change_tail_denotes_the_deletion` already applies
    // above: both sides of each equivalence below are oracle output, so on their
    // own they stay equal under a defect that moves BOTH by the same amount. If
    // `repeat_span` miscounted the reference tract by one, `Ala2[9]` and
    // `Ala11del` would shift together and the pair would still match. Build the
    // expected peptide from `residues` and compare against that.
    let mut one_fewer = residues.clone();
    one_fewer.remove(10); // 1-based 11, the tract's last `Ala`
    let shortened = sequence(&built, "Ala2[9]");
    assert_eq!(
        shortened,
        format!("{}Ter", ProteinFrame::spell(&one_fewer)),
        "`Ala2[9]` must leave a nine-residue tract at 2..=10, not merely agree with `Ala11del`"
    );
    assert_eq!(shortened, sequence(&built, "Ala11del"));

    let mut one_more = residues.clone();
    one_more.insert(10, AminoAcid::Ala); // an eleventh `Ala`, extending the tract
    let mut two_more = one_more.clone();
    two_more.insert(10, AminoAcid::Ala);
    let lengthened = sequence(&built, "Ala2[12]");
    assert_eq!(
        lengthened,
        format!("{}Ter", ProteinFrame::spell(&two_more)),
        "`Ala2[12]` must leave a twelve-residue tract at 2..=13"
    );
    assert_eq!(lengthened, sequence(&built, "Ala10_Ala11dup"));
}

// ---------------------------------------------------------------------------
// Refusals
// ---------------------------------------------------------------------------

/// A description whose stated residue is not the reference's residue denotes
/// nothing on this reference, and the oracle says which position and which two
/// residues rather than just declining.
#[test]
fn a_stated_residue_that_the_reference_contradicts_is_named() {
    let built = frame(&bench(30, &[(24, AminoAcid::Trp)]));
    assert_eq!(
        denote(&built, "Cys24Trp"),
        ProteinDenotation::NoSingleSequence(NoSingleSequence::ReferenceMismatch {
            position: 24,
            stated: "Cys".to_string(),
            actual: "Trp".to_string(),
        })
    );
}

/// A position past the terminator is out of range, and the terminator itself is
/// not — `p.Ter31=` is addressable on a 30-residue protein.
#[test]
fn a_position_past_the_terminator_is_out_of_range() {
    let built = frame(&bench(30, &[]));
    assert_eq!(
        denote(&built, "Ala32Gly"),
        ProteinDenotation::NoSingleSequence(NoSingleSequence::PositionOutOfRange {
            position: 32,
            length: 31,
        })
    );
    assert!(matches!(
        denote(&built, "Ter31="),
        ProteinDenotation::Sequence(_)
    ));
}

/// `protein/deletion.md:18`:
///
/// > the "positions_deleted" should be listed from **5' to 3'**, i.e.
/// > `Cys76_Glu79`, not `Glu79_Cys76`.
#[test]
fn spec_a_range_written_three_prime_to_five_prime_denotes_nothing() {
    let built = frame(&bench(90, &[(76, AminoAcid::Cys), (79, AminoAcid::Glu)]));
    assert_eq!(
        sequence(&built, "Cys76_Glu79del").len(),
        (90 - 4 + 1) * 3,
        "the 5'→3' form deletes four residues"
    );
    assert_eq!(
        denote(&built, "Glu79_Cys76del"),
        ProteinDenotation::NoSingleSequence(NoSingleSequence::ReversedRange)
    );
}

/// Two members of one cis allele that claim the same residue denote no single
/// protein — the amino-acid twin of the corpus's overlapping-member refusal.
#[test]
fn two_members_claiming_one_residue_denote_nothing() {
    let built = frame(&bench(50, &[(44, AminoAcid::Ser)]));
    assert_eq!(
        denote(&built, "[Ser44Arg;Ser44Cys]"),
        ProteinDenotation::NoSingleSequence(NoSingleSequence::OverlappingMembers)
    );
}

/// The oracle refuses a reference it is not holding, rather than scoring the
/// description against whatever protein it does hold.
#[test]
fn an_accession_the_frame_does_not_serve_is_refused() {
    let built = frame(&bench(30, &[(24, AminoAcid::Trp)]));
    assert_eq!(
        protein_denotation_of(built.provider(), "NP_999999.9:p.Trp24Cys"),
        ProteinDenotation::NoSingleSequence(NoSingleSequence::ReferenceUnavailable(
            "NP_999999.9".to_string()
        ))
    );
}

/// A description that is not on the protein axis, and one that does not parse
/// at all, are different answers.
#[test]
fn a_non_protein_description_is_not_a_parse_failure() {
    let built = frame(&bench(30, &[]));
    assert_eq!(
        protein_denotation_of(built.provider(), "NM_TEST.1:c.5A>T"),
        ProteinDenotation::NotProteinAxis
    );
    assert_eq!(
        protein_denotation_of(built.provider(), "this is not HGVS"),
        ProteinDenotation::Unparseable
    );
}

/// An uncertain **position** is refused even though an uncertain **edit** is
/// not. The two parentheses mean different things: `p.(Trp24Cys)` records
/// absent experimental evidence (`protein/substitution.md:16`), while a `?`
/// endpoint leaves the edit's territory open.
#[test]
fn an_uncertain_position_is_refused_where_a_predicted_edit_is_not() {
    let built = frame(&bench(60, &[(24, AminoAcid::Trp)]));
    assert!(matches!(
        denote(&built, "(Trp24Cys)"),
        ProteinDenotation::Sequence(_)
    ));
    assert_eq!(
        denote(&built, "Trp24_?del"),
        ProteinDenotation::NoSingleSequence(NoSingleSequence::UncertainPosition)
    );
}

// ---------------------------------------------------------------------------
// Round trip
// ---------------------------------------------------------------------------

/// Apply, then compare: a denoted protein is a real protein, and a description
/// that undoes the first one on that protein returns the original.
///
/// The second half goes back through the **provider** — the denoted protein is
/// built into a fresh frame — so it exercises the provider and the oracle
/// together rather than the applier alone. The pair is the spec's own
/// `p.Trp4del` / `p.Trp4dup` on `MetLeuTrpTrpGlu` (`protein/deletion.md:38`,
/// `protein/duplication.md:39`).
#[test]
fn a_deletion_and_the_duplication_that_undoes_it_round_trip() {
    let original = peptide("MetLeuTrpTrpGlu");
    let built = frame(&original);

    let deleted = sequence(&built, "Trp4del");
    assert_eq!(deleted, "MetLeuTrpGluTer");

    let shorter = frame(&peptide("MetLeuTrpGlu"));
    assert_eq!(
        sequence(&shorter, "Trp3dup"),
        format!("{}Ter", ProteinFrame::spell(&original)),
        "duplicating the surviving Trp restores the original protein"
    );
}

/// Every edit kind the oracle claims to handle denotes a protein that is itself
/// buildable — i.e. the oracle never emits a sequence the provider could not
/// serve. Cheap, and it catches an applier that produces a stray terminator or
/// drops one.
#[test]
fn every_denoted_protein_is_itself_a_buildable_protein() {
    let built = frame(&bench(
        60,
        &[
            (24, AminoAcid::Trp),
            (28, AminoAcid::Cys),
            (29, AminoAcid::Lys),
            (44, AminoAcid::Ser),
            (46, AminoAcid::Trp),
        ],
    ));
    let suffixes = [
        "Trp24Cys",
        "Trp24del",
        "Trp24dup",
        "Cys28_Lys29delinsTrp",
        "Cys28delinsTrpVal",
        "Trp24_Cys28del",
        "[Ser44Arg;Trp46Arg]",
        "=",
        "Trp24=",
    ];
    for suffix in suffixes {
        let denoted = sequence(&built, suffix);
        assert!(
            denoted.ends_with("Ter"),
            "p.{suffix} denoted {denoted}, which does not terminate"
        );
        let residues = ProteinFrame::peptide_from_three_letter(&denoted)
            .unwrap_or_else(|| panic!("p.{suffix} denoted an unparseable peptide: {denoted}"));
        let without_stop = &residues[..residues.len() - 1];
        assert!(
            !without_stop.contains(&AminoAcid::Ter),
            "p.{suffix} denoted a protein with an internal stop: {denoted}"
        );
        ProteinFrame::from_peptide(without_stop)
            .unwrap_or_else(|error| panic!("p.{suffix} denoted an unbuildable protein: {error}"));
    }
}

/// The three-exon minus-strand frame answers exactly as the single-exon one
/// does, so a future paired `c.`/`p.` stratum can vary the geometry without the
/// protein answers moving under it.
#[test]
fn the_oracle_answers_the_same_on_every_transcript_geometry() {
    let residues = bench(60, &[(24, AminoAcid::Trp)]);
    let answers: Vec<ProteinDenotation> = ProteinRefShape::all()
        .into_iter()
        .map(|shape| {
            let built = ProteinFrame::build(shape, &residues).expect("frame builds");
            denote(&built, "Trp24Cys")
        })
        .collect();
    assert!(answers.windows(2).all(|pair| pair[0] == pair[1]));
    assert!(matches!(answers[0], ProteinDenotation::Sequence(_)));
    assert_eq!(
        ProteinRefShape::ThreeExon(Strand::Minus).label(),
        "p3m",
        "labels stay distinguishable in a census"
    );
}

/// **Every** spelling of an implausible count is refused, and refused the same
/// way.
///
/// Four spellings state a length against this reference: `insXaa[n]`, `ins*n`,
/// `ext*n` and `ext-n`. Only `insXaa[n]` routed its count through
/// `checked_count`, so one implausible count was a refusal and the other three
/// were determinate lengths on the same reference — the oracle answering one
/// thing four different ways.
///
/// The two insertion arms were reconciled first and this test covered only
/// those, which left the class half-swept: the extension arms have the same
/// shape (a stated count with no residues) and were still unbounded. Widened to
/// the whole class rather than paired with a second test, so a fifth
/// count-stating arm has one obvious place to be added.
///
/// **Nothing overflowed, in any arm.** `ByCount` carries only a number and
/// `stated_insert_length` saturates; `count` is an `i64`, so the largest length
/// the extension arms could build was `i64::MAX + start + 1`, measured at
/// `Some(9223372036854775812)` rather than a panic. What was wrong is that
/// `RepeatCountOutOfRange`'s own doc claims counts are validated against the
/// reference before use, and for three of the four arms that was false. (The
/// C-terminal arm had a second, quieter failure: `usize::try_from(count)
/// .unwrap_or(0)` degraded an out-of-range count to a length of **zero** on a
/// 32-bit target, which is worse than refusing it.)
///
/// The bound is the reference length, which keeps `protein/insertion.md:49-51`'s
/// own `ins*63` against a 746-residue protein valid, and
/// `protein/extension.md:30-31`'s `p.Ter110Glnext*17` valid against its own —
/// both pinned below so the bound cannot be tightened into rejecting a
/// published example.
#[test]
fn every_spelling_of_an_implausible_count_is_refused() {
    let built = frame(&bench(30, &[]));
    let absurd = 4_000_000_000u64;

    for suffix in [
        format!("Gly2_Ser3insXaa[{absurd}]"),
        format!("Gly2_Ser3ins*{absurd}"),
        format!("Ter31Glnext*{absurd}"),
        format!("Met1ext-{absurd}"),
    ] {
        assert_eq!(
            denote(&built, &suffix),
            ProteinDenotation::NoSingleSequence(NoSingleSequence::RepeatCountOutOfRange {
                count: absurd,
                limit: 31,
            }),
            "p.{suffix} must be refused, not resolved to a determinate length"
        );
    }

    // The spec's own example shape stays REACHABLE: a stop-terminated insertion
    // whose count is within the reference length is still answered on its
    // merits. It is `Indeterminate` rather than a `Sequence` because `ins*n`
    // states a length and no residues — which is the distinction the refusals
    // above must not be confused with. Asserted so the new bound cannot be
    // tightened into rejecting `protein/insertion.md:49-51`'s published
    // `ins*63`.
    for suffix in ["Gly2_Ser3ins*20", "Ter31Glnext*17", "Met1ext-5"] {
        assert!(
            matches!(denote(&built, suffix), ProteinDenotation::Indeterminate(_)),
            "a count within the reference length must still be answered, not refused: \
             p.{suffix} gave {:?}",
            denote(&built, suffix)
        );
    }
}
