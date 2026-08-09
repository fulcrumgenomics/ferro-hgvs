//! A synthetic **protein** reference derived from a synthetic transcript, plus
//! an amino-acid-level denotation oracle over it.
//!
//! # The measurement gap this closes
//!
//! [`spec_corpus`](super::spec_corpus) makes a nucleotide description's
//! denotation checkable with no published answer: build a synthetic reference,
//! resolve every member through `hgvs_to_spdi`, apply the triples, compare. The
//! protein slice of the spec gets none of that. Its 149 clauses (102 of them
//! validity requirements) are blanket-exempted from the corpus with an honest
//! reason — the corpus's oracle is nucleotide-only, and a `p.` description's
//! denotation depends on **translation** rather than on splicing bases. The
//! blessed prepared reference cannot close it either: its `protein_fastas` is
//! empty, so it serves no protein sequence at all.
//!
//! This module is the enabling half. It supplies the two things a protein
//! stratum needs and nothing else:
//!
//! 1. [`ProteinFrame`] — a `MockProvider` serving an `NP_`-style protein that is
//!    **translated from the transcript it ships with**, so a `c.` row and its
//!    `p.` consequence are about one molecule.
//! 2. [`protein_denotation_of`] — the amino-acid analogue of
//!    [`spec_corpus::denotation_of`](super::spec_corpus::denotation_of).
//!
//! It deliberately adds **no corpus rows**. Wiring a protein stratum into
//! [`SpecCorpus`](super::spec_corpus::SpecCorpus) would move the censuses that
//! `tests/it/spec_conformance_axis.rs` pins, and is a separate change.
//!
//! # Why the coupling is enforced rather than documented
//!
//! A provider that serves an unrelated protein string is worse than no provider:
//! every comparison against it passes or fails for reasons that have nothing to
//! do with the transcript. So [`ProteinFrame::from_peptide`] back-translates the
//! requested peptide into a CDS, builds the transcript around it, and then
//! **re-translates the CDS out of the built transcript** and asserts the result
//! is the peptide it was asked for. The round trip is the guarantee; there is no
//! path that produces a frame whose protein is not its transcript's translation.
//!
//! Both directions reuse ferro's own codon table
//! ([`CodonTable`](crate::backtranslate::codon::CodonTable) and
//! [`crate::project::protein::translate`]) rather than a private copy, so a
//! frame cannot disagree with the projector about what a codon means.
//!
//! # How a denoted protein is spelled
//!
//! As a concatenation of **three-letter** codes — `MetLeuTrpGlu` — which is the
//! form the spec's own worked examples use (`protein/insertion.md:32`,
//! `protein/deletion.md:38`). Two reasons, both from
//! `background/standards.md`:
//!
//! - `:213` — "HGVS nomenclature uses `Ter` (three-letter amino acid code) and
//!   `*` (three- and one-letter amino acid code) to indicate a translation
//!   termination (stop) codon (in older versions, `X` was used instead)". A
//!   one-letter rendering would put `X` (the code for `Xaa`, `:239`) into the
//!   same alphabet as a stop, which is exactly the confusion the footnote at
//!   `:244` tells readers to avoid: "since 'X' has been used to indicate a
//!   translation stop codon, use 'Xaa' only".
//! - Every code in the table at `:216`–`:242` is exactly three characters, so
//!   the concatenation parses unambiguously.
//!
//! The **provider** stores the one-letter form, because that is what
//! `MockProvider::add_protein` and every existing ferro protein path use.

use std::fmt;

use crate::backtranslate::codon::CodonTable;
use crate::conformance::spec_corpus::{three_exon_layout, transcript_provider, CODING_ACCESSION};
use crate::hgvs::edit::{
    AminoAcidSeq, ExtDirection, FrameshiftTer, ProteinEdit, ProteinInsSeq, RepeatCount,
};
use crate::hgvs::interval::UncertainBoundary;
use crate::hgvs::location::{AminoAcid, ProtPos};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{AllelePhase, ProteinVariant};
use crate::hgvs::HgvsVariant;
use crate::parse_hgvs;
use crate::project::protein::translate;
use crate::reference::transcript::Strand;
use crate::reference::{MockProvider, ReferenceProvider};

// ---------------------------------------------------------------------------
// Synthetic protein reference conventions
// ---------------------------------------------------------------------------

/// The protein a synthetic `p.` row is drawn against.
///
/// Deliberately an `NP_` accession: `background/refseq.md` makes the protein
/// sequence accession the reference of a `p.` description, and every example in
/// `docs/recommendations/protein/*` uses that bare form.
pub const PROTEIN_ACCESSION: &str = "NP_TEST.1";

/// Untranslated bases either side of the CDS, in transcript coordinates.
///
/// Twelve is four codons' worth, so `c.-1`..`c.-12` and `c.*1`..`c.*12` all
/// exist and an N-terminal extension has somewhere to start.
const UTR_LEN: usize = 12;

/// The 5'UTR/3'UTR filler.
///
/// Period-4 `ACGT`, matching [`spec_corpus::padded`](super::spec_corpus::padded)
/// so a UTR cannot be mistaken for a designed tract. It carries no `ATG`, so it
/// introduces no upstream open reading frame that would make the frame's own
/// start codon ambiguous.
const UTR_FILLER: &str = "ACGT";

/// Which transcript geometry a protein frame's CDS is laid out on.
///
/// [`Self::SingleExon`] is the control and [`Self::ThreeExon`] is the shape
/// #1478 says a corpus must be able to vary: a CDS on one exon can never place
/// a residue across a splice junction, so a projection defect that only appears
/// at a junction is structurally invisible to it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProteinRefShape {
    /// One exon spanning the whole transcript, plus strand.
    SingleExon,
    /// Three exons on the given strand, so the CDS crosses two junctions.
    ThreeExon(Strand),
}

impl ProteinRefShape {
    /// Every shape, in a deterministic order.
    #[must_use]
    pub fn all() -> Vec<Self> {
        vec![
            Self::SingleExon,
            Self::ThreeExon(Strand::Plus),
            Self::ThreeExon(Strand::Minus),
        ]
    }

    /// Stable label, for row ids and test names.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::SingleExon => "p1",
            Self::ThreeExon(Strand::Minus) => "p3m",
            Self::ThreeExon(_) => "p3p",
        }
    }

    fn strand(self) -> Strand {
        match self {
            Self::SingleExon => Strand::Plus,
            Self::ThreeExon(strand) => strand,
        }
    }
}

/// Why a peptide cannot be made into a [`ProteinFrame`].
///
/// Every variant names a property of the *request*, never a partial frame: a
/// frame either round-trips or is not built.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinFrameError {
    /// The peptide is empty.
    Empty,
    /// Residue 1 is not `Met`. Translation initiation installs `Met` whatever
    /// the start codon is (`background/standards.md:229`, which lists `ATG` as
    /// "translation initiation"), so an authoritative `NP_` always begins with
    /// one and a frame that did not would be modelling something that does not
    /// occur.
    NotMetInitiated(AminoAcid),
    /// The peptide names a residue the standard codon table cannot spell —
    /// `Xaa` (`:239`, "unknown or 'other'"), or a residue whose only codon is
    /// read as a stop.
    NotBackTranslatable(AminoAcid),
    /// The peptide contains a `Ter` other than as an implicit terminator. The
    /// terminator is supplied by the builder; a `Ter` inside the peptide would
    /// truncate the translation and break the round trip.
    InternalStop(usize),
    /// The re-translation of the built CDS is not the requested peptide. Not
    /// reachable through the checks above; it exists so the round trip is a
    /// returned error rather than a comment.
    RoundTripMismatch {
        requested: String,
        translated: String,
    },
    /// The provider did not serve back the transcript it was just built with,
    /// so the transcript→protein link could not be recorded. Not reachable
    /// through a well-formed provider; it exists so a broken one is a returned
    /// error rather than a frame that silently carries no `protein_id`.
    TranscriptNotServed(String),
}

impl fmt::Display for ProteinFrameError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "a synthetic peptide must have at least one residue"),
            Self::NotMetInitiated(aa) => write!(
                f,
                "residue 1 must be Met (translation initiation), got {}",
                aa.to_three_letter()
            ),
            Self::NotBackTranslatable(aa) => write!(
                f,
                "{} has no codon in the standard table",
                aa.to_three_letter()
            ),
            Self::InternalStop(at) => {
                write!(f, "the peptide has an internal Ter at residue {at}")
            }
            Self::RoundTripMismatch {
                requested,
                translated,
            } => write!(
                f,
                "the built CDS translates to {translated}, not the requested {requested}"
            ),
            Self::TranscriptNotServed(accession) => write!(
                f,
                "the provider does not serve {accession} back after it was registered, so the \
                 transcript-to-protein link cannot be recorded"
            ),
        }
    }
}

impl std::error::Error for ProteinFrameError {}

/// A materialized synthetic protein reference: the provider, the protein it
/// serves, and the transcript that protein is the translation of.
///
/// The nucleotide twin is [`spec_corpus::Frame`](super::spec_corpus::Frame), and
/// the shape follows it deliberately — a provider plus the sequence the axis's
/// coordinates address plus the arithmetic that maps an index to an HGVS
/// position label.
#[derive(Clone)]
pub struct ProteinFrame {
    shape: ProteinRefShape,
    /// The reference protein, **without** a terminating `Ter` — the form an
    /// `NP_` record holds and the form the provider serves.
    residues: Vec<AminoAcid>,
    /// The whole transcript: 5'UTR, CDS (including its stop codon), 3'UTR.
    transcript: String,
    /// 1-based inclusive CDS bounds within [`Self::transcript`].
    cds: (usize, usize),
    provider: MockProvider,
}

impl ProteinFrame {
    /// Build a single-exon frame whose protein is `peptide`.
    ///
    /// # Errors
    ///
    /// See [`ProteinFrameError`].
    pub fn from_peptide(peptide: &[AminoAcid]) -> Result<Self, ProteinFrameError> {
        Self::build(ProteinRefShape::SingleExon, peptide)
    }

    /// Build a frame whose protein is the peptide spelled in concatenated
    /// three-letter codes, e.g. `"MetLeuTrpTrpGlu"`.
    ///
    /// Returns `None` when `spelled` is not a whole number of three-letter
    /// codes, or names one that is not in `background/standards.md:216`–`:242`.
    /// Kept separate from the fallible-for-other-reasons
    /// [`Self::from_peptide`] so a typo in a test's peptide is not reported as a
    /// modelling error.
    #[must_use]
    pub fn peptide_from_three_letter(spelled: &str) -> Option<Vec<AminoAcid>> {
        if !spelled.is_ascii() || !spelled.len().is_multiple_of(3) {
            return None;
        }
        (0..spelled.len() / 3)
            .map(|index| AminoAcid::from_three_letter(&spelled[index * 3..index * 3 + 3]))
            .collect()
    }

    /// Spell `residues` as concatenated three-letter codes.
    ///
    /// The inverse of [`Self::peptide_from_three_letter`], and the single place
    /// a denoted protein is turned into a string — see the module docs for why
    /// the one-letter form is not used.
    #[must_use]
    pub fn spell(residues: &[AminoAcid]) -> String {
        residues.iter().map(AminoAcid::to_three_letter).collect()
    }

    /// Build a frame with an explicit transcript geometry.
    ///
    /// # Errors
    ///
    /// See [`ProteinFrameError`].
    pub fn build(shape: ProteinRefShape, peptide: &[AminoAcid]) -> Result<Self, ProteinFrameError> {
        let cds_bases = back_translate_cds(peptide)?;

        let filler = |len: usize| UTR_FILLER.chars().cycle().take(len).collect::<String>();
        let utr5 = filler(UTR_LEN);
        let utr3 = filler(UTR_LEN);
        let transcript = format!("{utr5}{cds_bases}{utr3}");
        let cds = (UTR_LEN + 1, UTR_LEN + cds_bases.len());

        let exons = match shape {
            ProteinRefShape::SingleExon => vec![(1, transcript.len())],
            ProteinRefShape::ThreeExon(_) => three_exon_layout(transcript.len()),
        };
        let mut provider = transcript_provider(
            CODING_ACCESSION,
            shape.strand(),
            &transcript,
            Some(cds),
            &exons,
        );

        // Re-read the CDS out of the built transcript and translate it. This is
        // the coupling: what the provider serves as `NP_TEST.1` is what
        // `NM_TEST.1`'s CDS translates to, not a string handed in alongside it.
        let served_cds = &transcript[cds.0 - 1..cds.1];
        let residues = translate_cds(served_cds);
        if residues != peptide {
            return Err(ProteinFrameError::RoundTripMismatch {
                requested: Self::spell(peptide),
                translated: Self::spell(&residues),
            });
        }

        // Record the transcript→protein link in the data, not only in the
        // accession spelling, so a consumer that resolves the protein through
        // the transcript reaches the same record.
        //
        // The failure is NOT skipped. `transcript_provider` registered this
        // accession immediately above, so a lookup failure here means the
        // provider is broken — exactly the case that must be loud. Swallowing it
        // would return `Ok` with a frame whose transcript carries no
        // `protein_id`, which is the inconsistent frame this module's docs claim
        // cannot be built.
        let built = provider
            .get_transcript(CODING_ACCESSION)
            .map_err(|_| ProteinFrameError::TranscriptNotServed(CODING_ACCESSION.to_string()))?;
        let linked = (*built)
            .clone()
            .with_protein_id(Some(PROTEIN_ACCESSION.to_string()));
        provider.add_transcript(linked);
        provider.add_protein(
            PROTEIN_ACCESSION,
            residues
                .iter()
                .map(AminoAcid::to_one_letter)
                .collect::<String>(),
        );

        Ok(Self {
            shape,
            residues,
            transcript,
            cds,
            provider,
        })
    }

    /// The provider, serving `NP_TEST.1`, `NM_TEST.1` and the contig it sits on.
    #[must_use]
    pub fn provider(&self) -> &MockProvider {
        &self.provider
    }

    /// The transcript geometry the CDS is laid out on.
    #[must_use]
    pub fn shape(&self) -> ProteinRefShape {
        self.shape
    }

    /// The reference protein, without a terminating `Ter`.
    #[must_use]
    pub fn residues(&self) -> &[AminoAcid] {
        &self.residues
    }

    /// The reference protein **with** its terminating `Ter`, which is what a
    /// `p.` description addresses: `p.Ter110Glnext*17` names position 110 on a
    /// 109-residue protein (`protein/extension.md:30-31`).
    #[must_use]
    pub fn residues_with_stop(&self) -> Vec<AminoAcid> {
        let mut with_stop = self.residues.clone();
        with_stop.push(AminoAcid::Ter);
        with_stop
    }

    /// The reference protein spelled in three-letter codes, without the stop.
    #[must_use]
    pub fn spelled(&self) -> String {
        Self::spell(&self.residues)
    }

    /// The whole transcript sequence: 5'UTR, CDS, 3'UTR.
    #[must_use]
    pub fn transcript(&self) -> &str {
        &self.transcript
    }

    /// The CDS, including its stop codon.
    #[must_use]
    pub fn cds(&self) -> &str {
        &self.transcript[self.cds.0 - 1..self.cds.1]
    }

    /// 1-based inclusive CDS bounds within [`Self::transcript`].
    #[must_use]
    pub fn cds_bounds(&self) -> (usize, usize) {
        self.cds
    }

    /// The `c.` position of the first base of the codon for 1-based residue
    /// `residue`, or `None` when the residue is past the protein's terminator.
    ///
    /// This is the seam a future paired `c.`/`p.` stratum needs: it is what lets
    /// a nucleotide row and a protein row be stated about the same residue
    /// without either side re-deriving the CDS offset.
    #[must_use]
    pub fn cds_position_of_residue(&self, residue: u64) -> Option<i64> {
        let index = usize::try_from(residue).ok()?.checked_sub(1)?;
        // `residues.len() + 1` is the terminator, which has a codon too.
        if index > self.residues.len() {
            return None;
        }
        i64::try_from(index * 3 + 1).ok()
    }

    /// A `p.` description against this frame's protein accession.
    #[must_use]
    pub fn protein_descriptor(suffix: &str) -> String {
        format!("{PROTEIN_ACCESSION}:p.{suffix}")
    }

    /// A `c.` description against this frame's transcript accession.
    #[must_use]
    pub fn coding_descriptor(suffix: &str) -> String {
        format!("{CODING_ACCESSION}:c.{suffix}")
    }
}

/// Translate a CDS with ferro's own codon table, stopping **after** the first
/// terminator and dropping any trailing partial codon.
///
/// Unlike `project::protein::helpers::translate_full_cds` this keeps the `Ter`
/// out of the result: the reference protein an `NP_` record holds is the
/// residues, and the terminator is re-attached by
/// [`ProteinFrame::residues_with_stop`] where a description needs to address it.
fn translate_cds(cds: &str) -> Vec<AminoAcid> {
    let mut residues = Vec::with_capacity(cds.len() / 3);
    for index in 0..cds.len() / 3 {
        match translate(&cds[index * 3..index * 3 + 3]) {
            Some(AminoAcid::Ter) | None => break,
            Some(aa) => residues.push(aa),
        }
    }
    residues
}

/// Back-translate `peptide` into a CDS, terminator included.
///
/// Deterministic: the first codon the standard table lists for each residue.
/// The table is ferro's (`CodonTable::standard`), so a frame's CDS cannot encode
/// a residue the projector would translate differently.
fn back_translate_cds(peptide: &[AminoAcid]) -> Result<String, ProteinFrameError> {
    if peptide.is_empty() {
        return Err(ProteinFrameError::Empty);
    }
    if peptide[0] != AminoAcid::Met {
        return Err(ProteinFrameError::NotMetInitiated(peptide[0]));
    }
    if let Some(at) = peptide.iter().position(|aa| *aa == AminoAcid::Ter) {
        return Err(ProteinFrameError::InternalStop(at + 1));
    }

    let table = CodonTable::standard();
    let mut cds = String::with_capacity((peptide.len() + 1) * 3);
    for aa in peptide {
        let codon = table
            .codons_for(aa)
            .iter()
            // A residue whose only codon the table also reads as a stop (`Sec`,
            // `background/standards.md:236`) would truncate the translation, so
            // it is refused rather than silently mis-encoded.
            .find(|codon| !table.is_stop(codon))
            .ok_or(ProteinFrameError::NotBackTranslatable(*aa))?;
        cds.push_str(&codon.to_string());
    }
    let stop = table
        .stop_codons()
        .first()
        .ok_or(ProteinFrameError::NotBackTranslatable(AminoAcid::Ter))?;
    cds.push_str(&stop.to_string());
    Ok(cds)
}

// ---------------------------------------------------------------------------
// The denotation oracle
// ---------------------------------------------------------------------------

/// What the protein oracle can say about a `p.` description.
///
/// The four non-[`Sequence`](Self::Sequence) answers are kept apart for the same
/// reason [`spec_corpus::Denotation`](super::spec_corpus::Denotation) keeps
/// `NoSequence` and `Inexpressible` apart: collapsing them turns "the instrument
/// cannot express this" into "this description is wrong", and a census built on
/// that is wrong in the direction that flatters the implementation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProteinDenotation {
    /// The description denotes exactly this protein, spelled in concatenated
    /// three-letter codes and **including** its terminating `Ter` when
    /// translation reaches one.
    Sequence(String),
    /// The description fixes some of the protein and leaves the rest open. See
    /// [`Indeterminate`]. This is **not** a refusal: the determined part is
    /// still comparable, and the reason names exactly what is missing.
    Indeterminate(Indeterminate),
    /// `p.0` / `p.0?` — no protein is produced
    /// (`protein/substitution.md:46-48`). Distinct from an empty
    /// [`Sequence`](Self::Sequence), which would be a protein of length zero.
    NoProtein {
        /// `p.0?` rather than `p.0`.
        predicted: bool,
    },
    /// `parse_hgvs` rejected it.
    Unparseable,
    /// It parses, but it is not a `p.` description this oracle can be asked
    /// about — a nucleotide axis, or an allele mixing axes or accessions.
    NotProteinAxis,
    /// It is a `p.` description on this reference, and it denotes no single
    /// protein.
    NoSingleSequence(NoSingleSequence),
}

/// A description that determines part of a protein and not the rest.
///
/// The `prefix`/`length` split is what keeps this from being a guess: both
/// fields are read straight out of the description, and the residues between
/// them are simply absent rather than invented. A frameshift is the motivating
/// case — `p.Arg97ProfsTer23` fixes residues 1–96 (unchanged), residue 97
/// (`Pro`) and the total length (119), and says nothing at all about residues
/// 98–118, which are only recoverable from the *nucleotide* sequence in the
/// shifted frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Indeterminate {
    /// The longest determinate N-terminal prefix, in three-letter codes.
    pub prefix: String,
    /// Total residues including the terminating `Ter`, when the description
    /// states enough to fix it (`fsTer23`, `ext*17`, `insXaa[23]`). `None` for
    /// the open forms (`fs`, `fsTer?`, `ext*?`).
    pub length: Option<usize>,
    /// What is missing.
    pub reason: Indeterminacy,
}

/// Why part of a protein is not determined by its description.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Indeterminacy {
    /// A frameshift: the residues after the first new one are read in a
    /// different frame, which the protein reference does not carry
    /// (`protein/frameshift.md:18-22`).
    FrameshiftTail,
    /// A C-terminal extension past the terminator
    /// (`protein/extension.md:30-34`).
    CTerminalExtension,
    /// An N-terminal extension from an upstream initiation site
    /// (`protein/extension.md:22-24`); the added residues are upstream of the
    /// reference protein entirely, so no prefix is determined.
    NTerminalExtension,
    /// The description states residues by count rather than identity —
    /// `insXaa[23]` (`protein/insertion.md:45-47`) or `ins*63`
    /// (`:49-51`), whose note is explicit that the residues must be
    /// "deduce[d] … from the description given on the DNA or RNA level".
    ResiduesStatedByCount,
    /// An edit that reads past the terminator, so the residues that follow come
    /// from the 3'UTR (`protein/frameshift.md:25`).
    ReadsPastStop,
}

/// Why a `p.` description denotes no single protein on this reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NoSingleSequence {
    /// The provider serves no protein for the description's accession.
    ReferenceUnavailable(String),
    /// `p.?`, `p.(?)`, `p.Met1?` — the consequence is stated to be unknown
    /// (`protein/substitution.md:51-52`, `:71-73`).
    UnknownConsequence,
    /// A position is `?` or a bounded range, so the edit's territory is not
    /// fixed.
    UncertainPosition,
    /// The description names a set of proteins rather than one: `^`
    /// alternatives (`protein/substitution.md:77-78`), a mosaic `=/`
    /// (`:81-84`), a repeat count that is a range
    /// (`protein/repeated.md:31-33`), or an allele that is not in cis.
    SetValued,
    /// The residue the description states at a position is not the residue the
    /// reference has there.
    ReferenceMismatch {
        /// 1-based protein position.
        position: u64,
        /// The three-letter code the description states.
        stated: String,
        /// The three-letter code the reference has there.
        actual: String,
    },
    /// An insertion whose two positions do not flank one junction
    /// (`protein/insertion.md:17-18`: the positions "should contain **two
    /// flanking residues**, e.g. `Lys23_Leu24`, not two non-flanking residues
    /// (`Lys23_Asn25`)", and "an insertion can not be described using **one**
    /// amino acid position").
    NonFlankingInsertion,
    /// A position past the reference protein's terminator.
    PositionOutOfRange {
        /// 1-based protein position.
        position: u64,
        /// Residues in the reference, terminator included.
        length: usize,
    },
    /// The stated range runs 3'→5', which names no interval
    /// (`protein/deletion.md:18`: positions "should be listed from **5' to
    /// 3'**, i.e. `Cys76_Glu79`, not `Glu79_Cys76`").
    ReversedRange,
    /// A stated deleted sequence or length that the reference contradicts.
    StatedExtentMismatch,
    /// Two members of one allele claim the same residue, or two insertions
    /// claim one junction — the protein twin of
    /// [`spec_corpus::apply_triples`](super::spec_corpus::apply_triples)'s
    /// refusal.
    OverlappingMembers,
    /// A repeat or by-count insertion whose stated count is larger than the
    /// reference protein makes meaningful.
    ///
    /// `protein_denotation_of` is `pub` and takes an arbitrary descriptor, so
    /// the count is attacker-shaped input rather than only corpus-shaped input:
    /// `p.Ala2[4000000000]` would otherwise be materialised before it was
    /// judged. The count is validated against the reference **before** any
    /// allocation or multiplication, so an implausible one is a refusal rather
    /// than four billion residues or an overflowed multiply.
    RepeatCountOutOfRange {
        /// The count the description states.
        count: u64,
        /// The largest count this reference makes meaningful — its own length
        /// in residues, terminator included.
        limit: usize,
    },
    /// A shape this oracle does not model. Named rather than folded into
    /// another variant so a census can report it as instrument coverage rather
    /// than as a verdict on the description.
    Unsupported(&'static str),
}

/// Refuse a stated repeat count the reference cannot make meaningful, before it
/// is turned into a length or an allocation.
fn checked_count(count: u64, reference_len: usize) -> Result<usize, NoSingleSequence> {
    let limit = reference_len;
    match usize::try_from(count) {
        Ok(count) if count <= limit => Ok(count),
        _ => Err(NoSingleSequence::RepeatCountOutOfRange { count, limit }),
    }
}

/// What `descriptor` denotes on the protein `provider` serves for its accession.
///
/// The protein analogue of
/// [`spec_corpus::denotation_of`](super::spec_corpus::denotation_of). The
/// reference is looked up through the provider rather than passed in, so a
/// description naming an accession the frame does not serve is reported as
/// [`NoSingleSequence::ReferenceUnavailable`] instead of being silently scored
/// against the wrong molecule.
#[must_use]
pub fn protein_denotation_of(provider: &MockProvider, descriptor: &str) -> ProteinDenotation {
    let Ok(parsed) = parse_hgvs(descriptor) else {
        return ProteinDenotation::Unparseable;
    };
    let members = match &parsed {
        HgvsVariant::Protein(variant) => vec![variant.clone()],
        HgvsVariant::Allele(allele) => {
            // Only a **cis** allele describes one molecule. `[a;b]` is the
            // shape `protein/delins.md:62` publishes; trans, mosaic, `^` and
            // `,` each name more than one product.
            if allele.phase != AllelePhase::Cis {
                return ProteinDenotation::NoSingleSequence(NoSingleSequence::SetValued);
            }
            let mut protein_members = Vec::with_capacity(allele.variants.len());
            for member in &allele.variants {
                match member {
                    HgvsVariant::Protein(variant) => protein_members.push(variant.clone()),
                    _ => return ProteinDenotation::NotProteinAxis,
                }
            }
            if protein_members.is_empty() {
                return ProteinDenotation::NotProteinAxis;
            }
            protein_members
        }
        _ => return ProteinDenotation::NotProteinAxis,
    };

    let accession = members[0].accession.full();
    if members
        .iter()
        .any(|member| member.accession.full() != accession)
    {
        return ProteinDenotation::NotProteinAxis;
    }
    let Some(reference) = reference_protein(provider, &accession) else {
        return ProteinDenotation::NoSingleSequence(NoSingleSequence::ReferenceUnavailable(
            accession,
        ));
    };

    let mut spans = Vec::with_capacity(members.len());
    for member in &members {
        match member_span(&reference, member) {
            Ok(MemberOutcome::Span(span)) => spans.push(span),
            Ok(MemberOutcome::Terminal(denotation)) => {
                // `p.0`, an extension, a frameshift and the unknown forms each
                // describe the whole molecule, so they do not compose with a
                // sibling member. A cis allele containing one is refused rather
                // than resolved by dropping the sibling.
                return if members.len() == 1 {
                    denotation
                } else {
                    ProteinDenotation::NoSingleSequence(NoSingleSequence::Unsupported(
                        "a whole-protein member combined with another member",
                    ))
                };
            }
            Err(reason) => return ProteinDenotation::NoSingleSequence(reason),
        }
    }

    match apply_spans(&reference, &spans) {
        Some(applied) => {
            ProteinDenotation::Sequence(ProteinFrame::spell(&truncate_at_stop(&applied)))
        }
        None => ProteinDenotation::NoSingleSequence(NoSingleSequence::OverlappingMembers),
    }
}

/// The reference protein for `accession`, terminator included.
fn reference_protein(provider: &MockProvider, accession: &str) -> Option<Vec<AminoAcid>> {
    let length = usize::try_from(provider.get_protein_length(accession).ok()?).ok()?;
    if length == 0 {
        return None;
    }
    let sequence = provider
        .get_protein_sequence(accession, 0, u64::try_from(length).ok()?)
        .ok()?;
    let mut residues: Vec<AminoAcid> = sequence
        .chars()
        .map(AminoAcid::from_one_letter)
        .collect::<Option<_>>()?;
    residues.push(AminoAcid::Ter);
    Some(residues)
}

/// One member's replacement of `[start, start + removed)` with `inserted`,
/// in 0-based indices into the reference protein.
#[derive(Debug, Clone)]
struct Span {
    start: usize,
    removed: usize,
    inserted: Vec<AminoAcid>,
}

/// What reading one member produced.
enum MemberOutcome {
    /// An edit that composes with siblings.
    Span(Span),
    /// A description of the whole molecule, which does not.
    Terminal(ProteinDenotation),
}

/// Apply `spans` to `reference`, or decline.
///
/// The walk is [`spec_corpus::apply_triples`](super::spec_corpus::apply_triples)
/// transposed onto residues, and for the same three reasons: a 3'→5' walk with a
/// `claimed` cursor so an overlapping description is declined rather than
/// double-applied; a longer-removal-first tie-break so a zero-width member flush
/// against a deletion does not read as an overlap; and rejection of two pure
/// insertions at one junction, whose joint denotation is undefined.
fn apply_spans(reference: &[AminoAcid], spans: &[Span]) -> Option<Vec<AminoAcid>> {
    let mut ordered: Vec<&Span> = spans.iter().collect();
    ordered.sort_by_key(|span| {
        (
            std::cmp::Reverse(span.start),
            std::cmp::Reverse(span.removed),
        )
    });
    let mut edited = reference.to_vec();
    let mut claimed = reference.len();
    let mut insertion_at: Option<usize> = None;
    for span in ordered {
        let end = span.start.checked_add(span.removed)?;
        if end > reference.len() || end > claimed {
            return None;
        }
        if span.removed == 0 && insertion_at == Some(span.start) {
            return None;
        }
        edited.splice(span.start..end, span.inserted.iter().copied());
        if span.removed == 0 {
            insertion_at = Some(span.start);
        }
        claimed = span.start;
    }
    Some(edited)
}

/// Read one member into either a composable [`Span`] or a whole-molecule answer.
///
/// # The uncertainty wrapper is deliberately not a refusal
///
/// `Mu::Uncertain` on the edit is the `p.(Trp24Cys)` form, and
/// `protein/substitution.md:16` says what it means: "predicted consequences,
/// i.e. without experimental evidence (no RNA or protein sequence analysed),
/// should be given in parentheses". That is a statement about **evidence**, not
/// about which protein — `p.(Trp24Cys)` and `p.Trp24Cys` denote the same single
/// sequence. Refusing the parenthesised form would drop most of the spec's own
/// worked examples from any future census, on a misreading. `Mu::Unknown` (`?`
/// in the edit slot) is a different thing and is refused.
fn member_span(
    reference: &[AminoAcid],
    member: &ProteinVariant,
) -> Result<MemberOutcome, NoSingleSequence> {
    let edit = match &member.loc_edit.edit {
        Mu::Certain(edit) | Mu::Uncertain(edit) => edit,
        Mu::Unknown => return Err(NoSingleSequence::UnknownConsequence),
    };

    // Whole-molecule forms first: their location slot carries no meaning, so it
    // must not be resolved.
    match edit {
        ProteinEdit::NoProtein { predicted } => {
            return Ok(MemberOutcome::Terminal(ProteinDenotation::NoProtein {
                predicted: *predicted,
            }))
        }
        ProteinEdit::Unknown { .. } => return Err(NoSingleSequence::UnknownConsequence),
        ProteinEdit::Identity {
            whole_protein: true,
            ..
        } => {
            return Ok(MemberOutcome::Terminal(ProteinDenotation::Sequence(
                ProteinFrame::spell(&truncate_at_stop(reference)),
            )))
        }
        _ => {}
    }

    let start = resolve(reference, endpoint(&member.loc_edit.location.start)?)?;
    let end = resolve(reference, endpoint(&member.loc_edit.location.end)?)?;

    let span = match edit {
        ProteinEdit::Substitution {
            reference: stated,
            alternative,
        } => {
            if start != end {
                return Err(NoSingleSequence::Unsupported(
                    "a substitution names exactly one position",
                ));
            }
            check_residue(reference, start, *stated)?;
            Span {
                start,
                removed: 1,
                inserted: vec![*alternative],
            }
        }

        // `^` alternatives and a mosaic `=/` name a set of proteins, not one:
        // `protein/substitution.md:77-78` and `:81-84`.
        ProteinEdit::SubstitutionAlternatives { .. }
        | ProteinEdit::FrameshiftAlternatives { .. } => return Err(NoSingleSequence::SetValued),

        ProteinEdit::Deletion { sequence, count } => {
            let removed = range_len(start, end)?;
            check_stated_sequence(reference, start, removed, sequence.as_ref())?;
            check_stated_count(removed, *count)?;
            Span {
                start,
                removed,
                inserted: Vec::new(),
            }
        }

        ProteinEdit::Duplication => {
            let removed = range_len(start, end)?;
            let copy = &reference[start..start + removed];
            let mut inserted = copy.to_vec();
            inserted.extend_from_slice(copy);
            Span {
                start,
                removed,
                inserted,
            }
        }

        ProteinEdit::Insertion { sequence } => {
            if end != start + 1 {
                return Err(NoSingleSequence::NonFlankingInsertion);
            }
            match insertion_payload(sequence, reference.len())? {
                Payload::Residues(inserted) => Span {
                    start: end,
                    removed: 0,
                    inserted,
                },
                Payload::ByCount(added) => {
                    return Ok(MemberOutcome::Terminal(ProteinDenotation::Indeterminate(
                        Indeterminate {
                            prefix: ProteinFrame::spell(&reference[..end]),
                            length: added.map(|count| stated_insert_length(reference, end, count)),
                            reason: Indeterminacy::ResiduesStatedByCount,
                        },
                    )))
                }
            }
        }

        // One residue replaced by one is a substitution rather than a delins
        // (`protein/delins.md:17`), but that is a rule about *spelling*: the
        // denotation is well defined either way, and this is a denotation
        // oracle. Validity is a separate question and is not answered here.
        ProteinEdit::Delins { sequence } => {
            let removed = range_len(start, end)?;
            Span {
                start,
                removed,
                inserted: sequence.0.clone(),
            }
        }

        ProteinEdit::Identity { .. } => {
            let removed = range_len(start, end)?;
            Span {
                start,
                removed,
                inserted: reference[start..start + removed].to_vec(),
            }
        }

        ProteinEdit::Frameshift { new_aa, ter } => {
            let mut prefix = reference[..start].to_vec();
            prefix.extend(new_aa.iter().copied());
            // The stop's position is counted "starting at the first amino acid
            // changed by the frameshift (codon 1)" (`protein/frameshift.md:20`),
            // so a `Ter`n lands on protein position `start + n` with `start`
            // 0-based.
            let length = match ter {
                FrameshiftTer::At(count) => {
                    Some(start.saturating_add(usize::try_from(*count).unwrap_or(0)))
                }
                FrameshiftTer::Unknown | FrameshiftTer::Unspecified => None,
            };
            return Ok(MemberOutcome::Terminal(ProteinDenotation::Indeterminate(
                Indeterminate {
                    prefix: ProteinFrame::spell(&prefix),
                    length,
                    reason: Indeterminacy::FrameshiftTail,
                },
            )));
        }

        ProteinEdit::Extension {
            new_aa,
            direction,
            count,
        } => {
            return Ok(MemberOutcome::Terminal(extension(
                reference, start, *new_aa, *direction, *count,
            )?))
        }

        ProteinEdit::Repeat { sequence, count } => repeat_span(reference, start, sequence, count)?,

        ProteinEdit::MultiRepeat { .. } => {
            return Err(NoSingleSequence::Unsupported("a multi-unit protein repeat"))
        }

        // Handled above, before the location was resolved.
        ProteinEdit::NoProtein { .. } | ProteinEdit::Unknown { .. } => unreachable!(),
    };

    // An edit that removes the terminator without supplying a new one reads on
    // into the 3'UTR (`protein/frameshift.md:25`), which the protein reference
    // does not carry. The determinate part is still reportable.
    let terminator = reference.len() - 1;
    let removes_terminator = span.start <= terminator && terminator < span.start + span.removed;
    if removes_terminator && !span.inserted.contains(&AminoAcid::Ter) {
        let mut prefix = reference[..span.start].to_vec();
        prefix.extend(span.inserted.iter().copied());
        return Ok(MemberOutcome::Terminal(ProteinDenotation::Indeterminate(
            Indeterminate {
                prefix: ProteinFrame::spell(&prefix),
                length: None,
                reason: Indeterminacy::ReadsPastStop,
            },
        )));
    }
    Ok(MemberOutcome::Span(span))
}

/// The determinate endpoint of one interval boundary.
fn endpoint(boundary: &UncertainBoundary<ProtPos>) -> Result<ProtPos, NoSingleSequence> {
    match boundary {
        UncertainBoundary::Single(Mu::Certain(pos) | Mu::Uncertain(pos)) => Ok(*pos),
        UncertainBoundary::Single(Mu::Unknown) | UncertainBoundary::Range { .. } => {
            Err(NoSingleSequence::UncertainPosition)
        }
    }
}

/// The 0-based index `pos` names, checking both the bound and the residue.
fn resolve(reference: &[AminoAcid], pos: ProtPos) -> Result<usize, NoSingleSequence> {
    let number = usize::try_from(pos.number).unwrap_or(usize::MAX);
    if number == 0 || number > reference.len() {
        return Err(NoSingleSequence::PositionOutOfRange {
            position: pos.number,
            length: reference.len(),
        });
    }
    let index = number - 1;
    check_residue(reference, index, pos.aa)?;
    Ok(index)
}

/// Check the residue a description states at `index` against the reference.
///
/// `Xaa` is not tested: it is the code for "unknown or 'other'"
/// (`background/standards.md:239`), so a description using it is declining to
/// state the residue rather than stating a wrong one.
fn check_residue(
    reference: &[AminoAcid],
    index: usize,
    stated: AminoAcid,
) -> Result<(), NoSingleSequence> {
    let actual = reference[index];
    if stated == AminoAcid::Xaa || stated == actual {
        return Ok(());
    }
    Err(NoSingleSequence::ReferenceMismatch {
        position: index as u64 + 1,
        stated: stated.to_three_letter().to_string(),
        actual: actual.to_three_letter().to_string(),
    })
}

/// Residues in `[start, end]`, refusing a range written 3'→5'.
fn range_len(start: usize, end: usize) -> Result<usize, NoSingleSequence> {
    if end < start {
        return Err(NoSingleSequence::ReversedRange);
    }
    Ok(end - start + 1)
}

/// Check a deletion's explicitly-spelled residues against the reference.
fn check_stated_sequence(
    reference: &[AminoAcid],
    start: usize,
    removed: usize,
    stated: Option<&AminoAcidSeq>,
) -> Result<(), NoSingleSequence> {
    let Some(stated) = stated else { return Ok(()) };
    if stated.0.len() != removed || stated.0 != reference[start..start + removed] {
        return Err(NoSingleSequence::StatedExtentMismatch);
    }
    Ok(())
}

/// Check a deletion's explicit length suffix (`p.Lys228_Met259del32`).
fn check_stated_count(removed: usize, stated: Option<u64>) -> Result<(), NoSingleSequence> {
    match stated {
        Some(count) if usize::try_from(count).unwrap_or(usize::MAX) != removed => {
            Err(NoSingleSequence::StatedExtentMismatch)
        }
        _ => Ok(()),
    }
}

/// What an insertion's payload resolves to.
enum Payload {
    /// The residues are named.
    Residues(Vec<AminoAcid>),
    /// Only a count is named. `Some(n)` when the count is stated, `None` when
    /// even that is open.
    ByCount(Option<usize>),
}

/// Resolve an insertion payload.
///
/// `reference_len` bounds a stated repeat count; see [`checked_count`].
fn insertion_payload(
    sequence: &ProteinInsSeq,
    reference_len: usize,
) -> Result<Payload, NoSingleSequence> {
    match sequence {
        ProteinInsSeq::Literal(seq) => Ok(Payload::Residues(seq.0.clone())),
        ProteinInsSeq::Repeat { aa, count } => {
            let RepeatCount::Exact(count) = count else {
                return Err(NoSingleSequence::SetValued);
            };
            let count = checked_count(*count, reference_len)?;
            if *aa == AminoAcid::Xaa {
                // `p.Arg78_Gly79insXaa[23]`: the count is stated and the
                // residues are not (`protein/insertion.md:45-47`).
                Ok(Payload::ByCount(Some(count)))
            } else {
                Ok(Payload::Residues(vec![*aa; count]))
            }
        }
        // `p.Gln746_Lys747ins*63`: an inserted run ending in a stop at position
        // 63 of the insertion (`protein/insertion.md:49-51`). The length is
        // fixed, the residues are not.
        ProteinInsSeq::Stop { position } => Ok(Payload::ByCount(Some(
            usize::try_from(*position).unwrap_or(usize::MAX),
        ))),
    }
}

/// Total residues, terminator included, after inserting `count` unnamed
/// residues at the junction after 0-based `at`.
///
/// A `ins*n` payload terminates inside itself, so the result ends there; an
/// `insXaa[n]` payload does not, so the reference tail survives.
fn stated_insert_length(reference: &[AminoAcid], at: usize, count: usize) -> usize {
    reference.len().saturating_add(count).max(at + count)
}

/// The whole-molecule answer for an extension.
///
/// Both arms route their stated count through [`checked_count`], for the same
/// reason the `ProteinInsSeq::Stop` arm of [`insertion_payload`] does and for
/// the same reason it originally did not: without it the oracle answers one
/// implausible count three different ways. `insXaa[4000000000]` was refused as
/// [`NoSingleSequence::RepeatCountOutOfRange`], `ins*4000000000` was a
/// determinate length until that arm was bounded, and `ext*4000000000` was an
/// `Indeterminate` of length 4000000005 — three spellings of "a run this
/// reference cannot make meaningful", adjudicated three ways.
///
/// **Nothing overflowed, and the fix is not about overflow.** `count` is an
/// `i64` from the parser, so the largest `start + 1 + count` this can build is
/// `i64::MAX + start + 1`, which is comfortably inside `usize` on a 64-bit
/// target; measured, `ext*9223372036854775807` returned
/// `Some(9223372036854775812)` rather than panicking. On a 32-bit target the
/// old `usize::try_from(count).unwrap_or(0)` degraded to a silent **zero**,
/// which is worse than a refusal and is the second reason to bound it here.
///
/// The bound is `reference.len()`, which keeps the spec's own published example
/// valid: `p.Ter110Glnext*17` (`protein/extension.md:30-31`) states 17 against
/// a reference of at least 110 residues.
fn extension(
    reference: &[AminoAcid],
    start: usize,
    new_aa: Option<AminoAcid>,
    direction: ExtDirection,
    count: Option<i64>,
) -> Result<ProteinDenotation, NoSingleSequence> {
    match direction {
        ExtDirection::CTerminal => {
            let mut prefix = reference[..start].to_vec();
            prefix.extend(new_aa.iter().copied());
            // `p.Ter110Glnext*17` puts the new stop "at position 17 of the added
            // sequence" (`protein/extension.md:30-31`), i.e. protein position
            // 110 + 17 with `start` 0-based at 109.
            let length = match count.filter(|count| *count > 0) {
                Some(count) => {
                    Some(start + 1 + checked_count(count.unsigned_abs(), reference.len())?)
                }
                None => None,
            };
            Ok(ProteinDenotation::Indeterminate(Indeterminate {
                prefix: ProteinFrame::spell(&prefix),
                length,
                reason: Indeterminacy::CTerminalExtension,
            }))
        }
        // `p.Met1ext-5` adds residues *upstream* of residue 1
        // (`protein/extension.md:22-24`), so no N-terminal prefix is determined
        // even though the reference survives as a suffix.
        ExtDirection::NTerminal => {
            let length = match count {
                Some(count) => {
                    Some(reference.len() + checked_count(count.unsigned_abs(), reference.len())?)
                }
                None => None,
            };
            Ok(ProteinDenotation::Indeterminate(Indeterminate {
                prefix: String::new(),
                length,
                reason: Indeterminacy::NTerminalExtension,
            }))
        }
    }
}

/// The span a `p.Ala2[10]` repeat denotes.
///
/// The reference tract is measured rather than assumed: the description names
/// the unit and where it starts, and how many copies the *reference* holds is a
/// property of the reference (`protein/repeated.md:20-23`).
fn repeat_span(
    reference: &[AminoAcid],
    start: usize,
    unit: &AminoAcidSeq,
    count: &RepeatCount,
) -> Result<Span, NoSingleSequence> {
    let RepeatCount::Exact(count) = count else {
        return Err(NoSingleSequence::SetValued);
    };
    // Validate before anything is allocated or multiplied. `count *
    // unit.0.len()` below panics in debug and wraps in release once `count` is
    // large, and the `take` would materialise it regardless.
    let count = checked_count(*count, reference.len())?;
    if unit.0.is_empty() {
        return Err(NoSingleSequence::Unsupported("an empty repeat unit"));
    }
    let mut copies = 0usize;
    while reference
        .get(start + copies * unit.0.len()..start + (copies + 1) * unit.0.len())
        .is_some_and(|window| window == unit.0)
    {
        copies += 1;
    }
    if copies == 0 {
        return Err(NoSingleSequence::StatedExtentMismatch);
    }
    Ok(Span {
        start,
        removed: copies * unit.0.len(),
        inserted: unit
            .0
            .iter()
            .copied()
            .cycle()
            .take(count * unit.0.len())
            .collect(),
    })
}

/// Cut the protein at its first terminator, keeping it.
///
/// Translation stops there, so a nonsense substitution really does denote a
/// shorter molecule. `protein/substitution.md:17-18` rules only on how that is
/// *described* — "A nonsense variant is not described as a Deletion of the
/// C-terminal end of the protein" — not on what it denotes.
fn truncate_at_stop(residues: &[AminoAcid]) -> Vec<AminoAcid> {
    match residues.iter().position(|aa| *aa == AminoAcid::Ter) {
        Some(at) => residues[..=at].to_vec(),
        None => residues.to_vec(),
    }
}
