//! A codon-designed **protein** conformance corpus, enumerated deterministically.
//!
//! # The measurement gap this closes
//!
//! [`synthetic_protein`](super::synthetic_protein) supplies a protein reference
//! coupled to its transcript, and an amino-acid-level denotation oracle over it.
//! It deliberately adds no corpus rows, so the protein axis has been
//! *measurable* and not *measured*. This module is the rows.
//!
//! # Why the cores had to be redesigned rather than reused
//!
//! [`spec_corpus`](super::spec_corpus) draws its cores from a xorshift stream
//! over an `AT` or `ACGT` alphabet and lays a CDS over them at
//! `(tx_len/8 + 1, tx_len - tx_len/8)` — a *ratio*, not a codon boundary. Such a
//! core is an open reading frame only by accident: it need not start with a
//! start codon, its CDS length need not be a multiple of three, and any of its
//! interior codons may be a terminator, which truncates the translation. Almost
//! none of them translate to anything a `p.` description can be written against,
//! so a protein stratum could not have been drawn from them.
//!
//! [`codon_core`] designs the CDS instead: an `ATG`, a body of sense codons
//! drawn from ferro's own [`CodonTable`], and a terminator. It is handed to
//! [`ProteinFrame::from_cds`], which **translates** it — so the peptide is
//! derived from the codons rather than the codons back-derived from a peptide,
//! and the frame's round-trip check still guarantees that what the provider
//! serves as `NP_TEST.1` is what `NM_TEST.1`'s CDS translates to.
//!
//! Back-translation would not have done. [`ProteinFrame::build`] takes the
//! table's *first* codon for each residue, so every CDS it can build holds at
//! most 20 distinct codons and no synonymous pair at all — and the single most
//! interesting property a protein corpus has over a nucleotide one is
//! unreachable in it: **a tandem repeat of one residue spelled in four
//! different synonymous codons**, which is an ambiguous run on the protein axis
//! and no run whatsoever on the nucleotide axis. [`RUN_RESIDUE`]'s tract is
//! exactly that shape, and `tests::the_ambiguous_run_is_invisible_to_the_nucleotide_axis`
//! is the demonstration.
//!
//! # #1810's generator was read and NOT reused
//!
//! `tests/it/codon_frame_deep_window_offset.rs` builds `"ATG" + "GCT".repeat(199)`
//! for the *coding* axis's codon exception. Its two design choices are both
//! deliberately wrong for this corpus: a **uniform** codon (so every residue is
//! `Ala` and no residue identity varies) and **no UTR** with `cds_start = 1` (so
//! there is no 5'UTR for an N-terminal extension to start in and no 3'UTR for a
//! C-terminal one to read into). What is taken from it is its finding rather
//! than its code — that a site near `c.1` has its window clamped and therefore
//! its phase forced to zero — which is why [`CoreDesign::body_codons`] is large
//! enough to place the designed sites deep in the CDS.
//!
//! # The structural-blindness audit
//!
//! Five blindnesses have been found in this repository's generators, each
//! invisible until the one before it was fixed. A protein stratum is exactly the
//! kind of generator that acquires a sixth, so every property it varies is named
//! here **and demonstrated by a test**, and every property it provably cannot
//! reach is stated as a limit rather than left to read as coverage.
//!
//! | property | varied by | demonstrated by |
//! |---|---|---|
//! | residue position class | [`Site`] — initiator, second residue, deep interior, last residue, the terminator, and past the end | `every_site_class_is_reached` |
//! | reading-frame phase at a splice junction | [`CoreDesign::body_codons`] mod 3; see [`CoreDesign::junction_phase`] | `both_junction_phases_are_reached_and_so_is_neither` |
//! | transcript geometry | [`ProteinRefShape`] — one exon (the #1478 control) and three exons | `every_shape_is_reached` |
//! | strand | [`ProteinRefShape::ThreeExon`] on both strands | `every_shape_is_reached` |
//! | synonymous codon identity | [`codon_core`] draws filler codons from the whole table, not one per residue | `the_generator_reaches_most_of_the_codon_table` |
//! | protein-axis ambiguity with no nucleotide-axis ambiguity | [`RUN_RESIDUE`]'s synonymous tract | `the_ambiguous_run_is_invisible_to_the_nucleotide_axis` |
//! | terminator codon identity | [`CoreDesign::stop_index`] | `every_stop_codon_is_reached` |
//! | edit type | [`Kind`] | `every_kind_is_reached` |
//! | determinacy class | [`ProteinDenotation`] — a sequence, an [`Indeterminate`](super::synthetic_protein::Indeterminate), `p.0`, or none | `every_determinacy_class_is_reached` |
//! | member count | 1 and 2, in cis | `every_kind_is_reached` |
//! | epistemic claim, at one denotation | [`RowKind::Distinguished`] — a description beside its parenthesised twin | `every_kind_is_reached` |
//!
//! ## What it provably cannot reach — declared, not omitted
//!
//! - **`Sec` and `Pyl`.** `Sec`'s only codon is read as a terminator by the
//!   standard table and `Pyl` has none at all, so neither can appear in a
//!   designed CDS. Any defect specific to them is outside this corpus.
//! - **`Xaa` in the reference.** A reference residue is whatever a codon
//!   translates to, and no codon translates to "unknown". `Xaa` is reachable
//!   only in a *description* (`insXaa[n]`), which the corpus does generate.
//! - **A non-standard genetic code.** [`CodonTable::standard`] is the only table
//!   used, so the mitochondrial reassignments are unreachable and an `m.`-axis
//!   protein defect is invisible here.
//! - **An incomplete or non-triplet CDS.** [`ProteinFrame::from_cds`] refuses
//!   one, by construction — the frame's round trip is what makes the oracle's
//!   answers mean anything, and it cannot hold for a CDS that does not
//!   translate. A transcript annotated with a partial CDS is therefore a real
//!   shape this corpus cannot express.
//! - **Codon-level ambiguity *within* one residue.** Two spellings of a `p.`
//!   description cannot differ in codon, because a `p.` description names no
//!   codons — see the next section.
//! - **An intronic or UTR position.** `p.` numbering has no zone outside the
//!   protein, so [`Region`](super::spec_corpus::Region)'s nine placements
//!   collapse to the seven in [`Site`] here, and the two axes' position
//!   taxonomies overlap only at "somewhere in the middle".
//!
//! ## And one limit that is a property of the AXIS, not of the generator
//!
//! A `p.` description has no codons of its own — `tests/it/protein_axis_split_move.rs`
//! states it in as many words, arguing from `general.md:35-38` being a
//! *reading-frame* exception. So the phase and strand dimensions above are
//! varied in the **frame** and may well be inert in a `p.`-only measurement.
//! That is measured rather than assumed: `shape_is_inert_for_the_protein_axis`
//! in `tests/it/protein_conformance_axis.rs` reports whether any row's outcome
//! differs across the three geometries. A measured inertness is a fact about the
//! axis worth recording; an *assumed* one would be the sixth blindness.
//!
//! # What a green run does NOT say
//!
//! The oracle judges **denotation** and not validity, deliberately: a
//! one-for-one `delins`, a `Cys76_Cys76` range and `p.[Arg76Ser;Cys77Trp]` all
//! denote a sequence perfectly well while the recommendations rule against each
//! spelling. The protein slice is **205 clause units** across nine files —
//! re-derived here from the checkout with the inventory's own definition, "a
//! line opening a bullet (`- `/`* `), a numbered item (`1. `), or an admonition
//! block (`!!! `)". **All 205 are classified `not-generatable`** by the
//! inventory, and this change does not reclassify one of them: a denotation
//! oracle cannot test a requirement on how a description is *written*, and that
//! is what a validity clause is.
//!
//! **A figure inherited from `synthetic_protein.rs` does not reconcile and is
//! not restated here.** That module's docs say "149 clauses (102 of them
//! validity requirements)". The committed rule inventory
//! (`tests/fixtures/spec-corpus/spec_rule_inventory.json`) counts **205** clause
//! units over `docs/recommendations/protein/*` at spec checkout `6f85311`, which
//! an independent re-derivation reproduces file for file; and nothing in the
//! repository classifies a clause as a "validity requirement", so the 102 has no
//! recoverable denominator either. Quote 205, or re-derive; do not quote 149.
//!
//! The rule inventory's `not-generatable` classification for
//! `docs/recommendations/protein/*` is therefore **left in place** by this
//! module, and must not be lifted on the strength of a green census.

use std::collections::BTreeMap;
use std::fmt;

use crate::backtranslate::codon::{Codon, CodonTable};
use crate::conformance::synthetic_protein::{
    protein_denotation_of, ProteinDenotation, ProteinFrame, ProteinRefShape,
};
use crate::hgvs::location::AminoAcid;

// ---------------------------------------------------------------------------
// The designed core
// ---------------------------------------------------------------------------

/// The residue the ambiguous tract is built from.
///
/// `Ala` because the standard table gives it four codons (`GCA/GCC/GCG/GCT`),
/// which is enough for every residue of a four-copy tract to be spelled
/// differently — the property that makes the tract ambiguous on the protein axis
/// and not on the nucleotide one.
pub const RUN_RESIDUE: AminoAcid = AminoAcid::Ala;

/// Copies of [`RUN_RESIDUE`] in the ambiguous tract.
///
/// Four, because that is the largest tract every four-codon residue can spell
/// without repeating a codon, and because a tract of `n` copies yields `n`
/// equivalent single-residue deletions and `n` equivalent duplications — the
/// corpus's whole confluence signal comes from here.
pub const RUN_COPIES: usize = 4;

/// Residues the filler is drawn from, in a fixed order.
///
/// Chosen so no two adjacent filler residues are equal (which would create an
/// undesigned ambiguous run beside the designed one) and so that none of them is
/// [`RUN_RESIDUE`] (which would extend the designed tract past its intended
/// bounds). `Met` is included so an **interior** `Met` exists and start-codon
/// reasoning has something other than residue 1 to look at.
const FILLER: &[AminoAcid] = &[
    AminoAcid::Arg,
    AminoAcid::Ser,
    AminoAcid::Trp,
    AminoAcid::Gly,
    AminoAcid::Cys,
    AminoAcid::Met,
    AminoAcid::Leu,
    AminoAcid::Thr,
    AminoAcid::Tyr,
    AminoAcid::Pro,
    AminoAcid::Val,
    AminoAcid::His,
    AminoAcid::Phe,
    AminoAcid::Gln,
    AminoAcid::Ile,
    AminoAcid::Asn,
    AminoAcid::Glu,
    AminoAcid::Lys,
    AminoAcid::Asp,
];

/// 1-based residue the ambiguous tract starts at.
///
/// Deep enough that the tract sits past the regime
/// `tests/it/codon_frame_deep_window_offset.rs` shows is structurally blind: a
/// site within `CANONICAL_PAD` (128) of the transcript start has its canonical
/// window clamped there, which forces the trim offset to a multiple of three and
/// therefore the codon phase to zero. Residue 50's codon starts at `c.148`,
/// transcript position 160, which clears it.
///
/// **It buys nothing for the `p.`-only census, and it is here on purpose.** A
/// protein description names no codons, so the nucleotide canonical window is
/// not consulted at all for the rows this corpus measures today. What the depth
/// removes is a caveat: any later stratum pairing a `c.` spelling with its `p.`
/// consequence inherits a frame that is already out of the clamped regime,
/// rather than one whose zeros would have to be re-derived.
///
/// See `tests::the_tract_clears_the_canonical_window_clamp`, which checks the
/// arithmetic rather than trusting this comment.
pub const RUN_START: u64 = 50;

/// One designed reading frame.
///
/// Three knobs, and each one exists to vary a property that would otherwise be
/// structurally unreachable — see [`Self::junction_phase`] for the third, which
/// is the subtle one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CoreDesign {
    /// Which deterministic draw the filler residues and their codons come from.
    pub seed: u32,
    /// Codons in the CDS, terminator included.
    ///
    /// **This is the phase knob.** See [`Self::junction_phase`].
    pub body_codons: usize,
    /// Which of the three terminators ends the CDS.
    pub stop_index: usize,
}

impl CoreDesign {
    /// Where the exon-1/exon-2 junction falls relative to a codon, on a
    /// [`ProteinRefShape::ThreeExon`] frame: `0` when it falls on a codon
    /// boundary, otherwise the number of the straddling codon's bases that lie
    /// in exon 1.
    ///
    /// # Why this is derivable at all, which is the point
    ///
    /// [`ProteinFrame`] wraps the CDS in a fixed 12-base UTR each side and lays
    /// the result out with `three_exon_layout`, whose first exon ends at
    /// `tx_len / 3`. With `c` codons, `tx_len = 24 + 3c`, so the junction sits
    /// after transcript position `8 + c` — and residue `r`'s codon starts at
    /// transcript position `3r + 10`. Solving `8 + c = 3r + 9 + k` gives
    /// `c = 3r + 1 + k`, so the phase is `(c - 1) mod 3` and **nothing else**.
    ///
    /// A generator that fixed the CDS length would therefore have fixed the
    /// phase, and every row it built would have agreed with every other on the
    /// one dimension a splice-junction defect lives in. That is the #1478 shape
    /// exactly, one level down. Varying `body_codons` over three consecutive
    /// values reaches all three phases; [`all_designs`] does.
    #[must_use]
    pub fn junction_phase(self) -> usize {
        (self.body_codons + 2) % 3
    }

    /// Stable label, used in row ids.
    #[must_use]
    pub fn label(self) -> String {
        format!("s{}c{}t{}", self.seed, self.body_codons, self.stop_index)
    }

    /// Residues in the peptide the design translates to — the terminator is a
    /// codon of the CDS and not a residue of the protein.
    #[must_use]
    pub fn residues(self) -> usize {
        self.body_codons - 1
    }
}

/// Every design the corpus enumerates.
///
/// Three consecutive `body_codons` so all three junction phases are reached
/// (see [`CoreDesign::junction_phase`]), each with its own terminator so the
/// three stop codons are reached too, and two seeds so no measured figure rests
/// on one draw of the filler.
#[must_use]
pub fn all_designs() -> Vec<CoreDesign> {
    let mut designs = Vec::new();
    for seed in 0..2u32 {
        for (offset, stop_index) in (0..3usize).enumerate() {
            designs.push(CoreDesign {
                seed,
                // 100, 101, 102 — phases 0, 1, 2 respectively, and all of them
                // long enough to put `RUN_START`'s tract past the canonical
                // window clamp and to leave a determinate interior, a last
                // residue and a terminator beyond it.
                body_codons: 100 + offset,
                stop_index,
            });
        }
    }
    designs
}

/// Build the CDS for `design`: a start codon, a body of sense codons holding the
/// designed ambiguous tract, and a terminator.
///
/// The tract's `RUN_COPIES` residues are all [`RUN_RESIDUE`] and each is spelled
/// with a **different** synonymous codon, which is the corpus's one shape that a
/// nucleotide-axis corpus structurally cannot hold. Everything else is filler,
/// drawn so that no residue adjacent to the tract equals [`RUN_RESIDUE`] — an
/// accidental fifth copy would silently change every ambiguity the tract is
/// there to create.
///
/// # Panics
///
/// If `design.body_codons` is too small to hold the tract at [`RUN_START`], or
/// if ferro's own codon table cannot spell a residue this module names — both of
/// which are programming errors in this module rather than inputs.
#[must_use]
pub fn codon_core(design: CoreDesign) -> String {
    let table = CodonTable::standard();
    let run_codons = table.codons_for(&RUN_RESIDUE);
    assert!(
        run_codons.len() >= RUN_COPIES,
        "{RUN_RESIDUE:?} has {} codons, too few for a {RUN_COPIES}-copy tract",
        run_codons.len()
    );
    let residues = design.residues();
    let run_start = usize::try_from(RUN_START).expect("RUN_START fits a usize");
    assert!(
        run_start + RUN_COPIES <= residues,
        "design {} has {residues} residues, too few for a tract at {run_start}",
        design.label()
    );

    let mut state = u64::from(design.seed)
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(0x0123_4567_89AB_CDEF)
        | 1;
    let mut draw = move |modulus: usize| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % modulus as u64) as usize
    };

    let mut cds = String::with_capacity(design.body_codons * 3);
    // Residue 1 is the initiator. `ATG` is the only `Met` codon, so this is not
    // a choice.
    cds.push_str(&codon_for(&table, AminoAcid::Met, 0).to_string());
    let mut previous = AminoAcid::Met;
    for residue in 2..=residues {
        let (aa, codon) = if (run_start..run_start + RUN_COPIES).contains(&residue) {
            let copy = residue - run_start;
            (RUN_RESIDUE, run_codons[copy].clone())
        } else {
            // Never `RUN_RESIDUE` and never the residue before it, so the only
            // tandem run in the protein is the designed one.
            let mut aa = FILLER[draw(FILLER.len())];
            let mut guard = 0usize;
            while aa == previous && guard < FILLER.len() {
                aa = FILLER[draw(FILLER.len())];
                guard += 1;
            }
            let options = table.codons_for(&aa);
            let pick = draw(options.len().max(1));
            (aa, codon_for(&table, aa, pick))
        };
        previous = aa;
        cds.push_str(&codon.to_string());
    }
    let stops = table.stop_codons();
    cds.push_str(&stops[design.stop_index % stops.len()].to_string());
    cds
}

/// The `pick`-th codon the standard table lists for `aa`.
///
/// # Panics
///
/// If the table lists none, or lists only terminators — neither of which is true
/// for any residue [`FILLER`] or [`RUN_RESIDUE`] names.
fn codon_for(table: &CodonTable, aa: AminoAcid, pick: usize) -> Codon {
    let options: Vec<&Codon> = table
        .codons_for(&aa)
        .iter()
        .filter(|codon| !table.is_stop(codon))
        .collect();
    assert!(
        !options.is_empty(),
        "the standard table cannot spell {aa:?} without a terminator"
    );
    options[pick % options.len()].clone()
}

/// Build the frame for `design` on `shape`.
///
/// # Panics
///
/// If the designed CDS does not round-trip through [`ProteinFrame::from_cds`],
/// which would mean [`codon_core`] built something that is not an open reading
/// frame — a defect in this module, and one that must be loud rather than
/// silently dropped from the corpus.
#[must_use]
pub fn frame_for(shape: ProteinRefShape, design: CoreDesign) -> ProteinFrame {
    let cds = codon_core(design);
    ProteinFrame::from_cds(shape, &cds)
        .unwrap_or_else(|why| panic!("designed core {} is not a CDS: {why}", design.label()))
}

// ---------------------------------------------------------------------------
// The row taxonomy
// ---------------------------------------------------------------------------

/// Where in the protein a row's edit sits.
///
/// The protein analogue of [`Region`](super::spec_corpus::Region), and shorter
/// for a reason worth stating: `p.` numbering has exactly one zone, so the
/// 5'UTR / intronic / junction placements a `c.` row can take do not exist
/// here. What replaces them is the terminator, which a `c.` row has no analogue
/// of at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Site {
    /// Residue 1, the initiator `Met`.
    ///
    /// **A substitution here is structurally undescribable, and that is
    /// correct.** `protein/substitution.md:49` forbids `p.Met1Xxx` and ferro's
    /// parser enforces it, naming the legal alternatives (`p.0`, `p.0?`,
    /// `p.(Met1?)`, or the insertion form). So the substitution and uncertainty
    /// strata both drop their initiator row, 36 drops in all, and the drop is
    /// evidence the clause is enforced rather than a hole in the corpus. The
    /// site is still reached, by `p.0`, `p.0?` and the N-terminal extension.
    Initiator,
    /// Residue 2 — adjacent to the initiator, so an edit here is the first that
    /// does not engage `protein/substitution.md:49`'s `p.Met1` rule.
    SecondResidue,
    /// Inside the designed ambiguous tract.
    AmbiguousRun,
    /// A determinate interior residue, deep in the CDS.
    DeepInterior,
    /// The last residue before the terminator.
    LastResidue,
    /// The terminator itself, at protein position `len + 1`.
    Terminator,
    /// A position past the terminator, which the reference does not have.
    PastEnd,
}

impl Site {
    /// Every site, in a deterministic order.
    #[must_use]
    pub fn all() -> Vec<Self> {
        vec![
            Self::Initiator,
            Self::SecondResidue,
            Self::AmbiguousRun,
            Self::DeepInterior,
            Self::LastResidue,
            Self::Terminator,
            Self::PastEnd,
        ]
    }

    /// Stable label, used in row ids and censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Initiator => "initiator",
            Self::SecondResidue => "second",
            Self::AmbiguousRun => "run",
            Self::DeepInterior => "interior",
            Self::LastResidue => "last",
            Self::Terminator => "terminator",
            Self::PastEnd => "past-end",
        }
    }

    /// The 1-based protein position this site names on `frame`.
    #[must_use]
    pub fn position(self, frame: &ProteinFrame) -> u64 {
        let residues = frame.residues().len() as u64;
        match self {
            Self::Initiator => 1,
            Self::SecondResidue => 2,
            Self::AmbiguousRun => RUN_START,
            // Between the tract and the last residue, so it is neither.
            Self::DeepInterior => RUN_START + RUN_COPIES as u64 + 2,
            Self::LastResidue => residues,
            Self::Terminator => residues + 1,
            Self::PastEnd => residues + 2,
        }
    }
}

/// The edit types a row is built from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Kind {
    /// A one-residue substitution, and its `delins` and parenthesised spellings.
    Substitution,
    /// A deletion inside the ambiguous tract, spelled at every equivalent
    /// position.
    Deletion,
    /// A duplication inside the ambiguous tract, likewise.
    Duplication,
    /// An insertion of the tract's own residue at a junction inside it.
    Insertion,
    /// A repeat count over the tract, against the explicit del/dup spellings of
    /// the same protein.
    Repeat,
    /// An edit that denotes the reference.
    Identity,
    /// A frameshift.
    Frameshift,
    /// An extension, N- or C-terminal.
    Extension,
    /// A payload stated by count rather than by identity.
    ByCount,
    /// `p.0` / `p.0?`.
    NoProtein,
    /// Two members in cis.
    CisAllele,
    /// A description beside its parenthesised, evidence-declining twin.
    Uncertainty,
}

impl Kind {
    /// Stable label, used in row ids and censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Substitution => "sub",
            Self::Deletion => "del",
            Self::Duplication => "dup",
            Self::Insertion => "ins",
            Self::Repeat => "repeat",
            Self::Identity => "identity",
            Self::Frameshift => "fs",
            Self::Extension => "ext",
            Self::ByCount => "by-count",
            Self::NoProtein => "no-protein",
            Self::CisAllele => "cis",
            Self::Uncertainty => "uncertainty",
        }
    }
}

/// Which properties a row can be asked about.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum RowKind {
    /// Two or more spellings the oracle agrees denote one protein **and** which
    /// make the same claim about it. Confluence, idempotency and protein
    /// preservation all apply.
    Family,
    /// One spelling. Idempotency and protein preservation apply; confluence is
    /// vacuous and is reported as such rather than counted as a pass.
    Single,
    /// Two or more spellings that denote one protein and **must not** converge,
    /// because they make different *epistemic* claims about it.
    ///
    /// # Why this is a kind of its own rather than a family
    ///
    /// `rulings[confluence-gate-is-apply-equality-on-every-determined-axis]`
    /// states the hazard in its own words: "two descriptions can be apply-equal
    /// on every axis and still make different epistemic claims, which is a
    /// canonical-form question and must not be encoded as a rung." The
    /// motivating pair is here — `protein/substitution.md:16` says "predicted
    /// consequences, i.e. without experimental evidence (no RNA or protein
    /// sequence analysed), should be given in parentheses", so `p.(Arg727Ser)`
    /// asserts less than `p.Arg727Ser` while denoting the same molecule.
    ///
    /// The oracle deliberately reads through the wrapper — refusing the
    /// parenthesised form would drop most of the spec's own worked examples
    /// from any census — so denotation equality **cannot** be the confluence
    /// relation for this pair. Filed as a `Family` it produced 90 of the first
    /// run's 90 `split_two` rows, every one of them correct behaviour: a
    /// confluence figure whose bulk is the normalizer doing the right thing
    /// hides the rows where it does not.
    ///
    /// The property measured over these is **preservation**, not convergence.
    Distinguished,
    /// A description that denotes no single protein. The property is that the
    /// implementation refuses it.
    Conflict,
}

impl fmt::Display for RowKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Family => "family",
            Self::Single => "single",
            Self::Distinguished => "distinguished",
            Self::Conflict => "conflict",
        })
    }
}

/// How determinate a row's denotation is.
///
/// **[`Self::Indeterminate`] is a first-class verdict and gets its own census
/// column.** That is not a stylistic choice: `rulings[confluence-gate-is-apply-equality-on-every-determined-axis]`
/// rules that "confluence is asserted only over DECIDED pairs; the
/// `Indeterminate` count is reported alongside and never folded into either
/// side". Folding it into `converged` would flatter the implementation, and
/// folding it into the divergent side would invent failures — the description
/// genuinely does not fix the protein, and no normalizer can be blamed for that.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Determinacy {
    /// The description fixes the whole protein.
    Determined,
    /// It fixes part of it and leaves the rest open.
    Indeterminate,
    /// It says no protein is produced.
    NoProtein,
    /// It denotes no single protein at all.
    Undenoted,
}

impl Determinacy {
    /// Read the class off an oracle answer.
    #[must_use]
    pub fn of(denotation: &ProteinDenotation) -> Self {
        match denotation {
            ProteinDenotation::Sequence(_) => Self::Determined,
            ProteinDenotation::Indeterminate(_) => Self::Indeterminate,
            ProteinDenotation::NoProtein { .. } => Self::NoProtein,
            ProteinDenotation::Unparseable
            | ProteinDenotation::NotProteinAxis
            | ProteinDenotation::NoSingleSequence(_) => Self::Undenoted,
        }
    }

    /// Stable label, used in censuses.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Determined => "determined",
            Self::Indeterminate => "indeterminate",
            Self::NoProtein => "no-protein",
            Self::Undenoted => "undenoted",
        }
    }
}

/// One corpus row.
#[derive(Debug, Clone)]
pub struct Row {
    /// Deterministic identifier, derived from the design parameters. Unique
    /// within a corpus and independent of iteration order, so a failing row can
    /// be named in a committed regression test.
    pub id: String,
    /// Which stratum enumerated it.
    pub stratum: &'static str,
    /// Which properties apply.
    pub kind: RowKind,
    /// The transcript geometry the protein's CDS is laid on.
    pub shape: ProteinRefShape,
    /// The designed reading frame.
    pub design: CoreDesign,
    /// Where in the protein the edit sits.
    pub site: Site,
    /// The edit type.
    pub edit: Kind,
    /// Members in the authored description.
    pub members: usize,
    /// What the authored spelling denotes, per the oracle.
    pub denoted: ProteinDenotation,
    /// Distinct spellings the oracle agrees denote [`Self::denoted`]. At least
    /// two for a [`RowKind::Family`], exactly one otherwise.
    pub spellings: Vec<String>,
    /// Clauses this row exercises, for a coverage join.
    pub rules: Vec<&'static str>,
}

impl Row {
    /// The design as authored — the first spelling.
    #[must_use]
    pub fn authored_spelling(&self) -> &str {
        self.spellings
            .first()
            .map_or("", std::string::String::as_str)
    }

    /// How determinate the row's denotation is.
    #[must_use]
    pub fn determinacy(&self) -> Determinacy {
        Determinacy::of(&self.denoted)
    }

    /// The frame the row's coordinates resolve against.
    ///
    /// Rebuilt from [`Self::shape`] and [`Self::design`] rather than stored, so
    /// the corpus does not carry one synthetic transcript per row. Rows sharing
    /// a `(shape, design)` share a frame, which is what keeps a census's cost
    /// linear in rows rather than in bases.
    #[must_use]
    pub fn frame(&self) -> ProteinFrame {
        frame_for(self.shape, self.design)
    }

    /// The key rows sharing one frame agree on.
    #[must_use]
    pub fn frame_key(&self) -> (ProteinRefShape, CoreDesign) {
        (self.shape, self.design)
    }
}

/// Why a design produced no row.
///
/// Counted rather than skipped, for the reason
/// [`spec_corpus::DropReason`](super::spec_corpus::DropReason) is: a generator
/// whose designs all collapsed into one of these would produce an empty corpus
/// that reads as "nothing to find".
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum DropReason {
    /// The site does not exist on this frame.
    SiteUnavailable,
    /// A spelling the generator built does not parse, so it cannot be measured.
    Unspellable,
    /// The oracle could not name what the authored spelling denotes, in a
    /// stratum that expected an answer.
    NoDenotation,
    /// Fewer than two spellings survived the oracle's agreement filter, so
    /// confluence is not decidable for the row.
    Singleton,
    /// The design denotes the reference protein, in a stratum that expected a
    /// change.
    DenotesTheReference,
}

impl DropReason {
    /// Stable label, grouped in the census.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::SiteUnavailable => "site unavailable",
            Self::Unspellable => "spelling does not parse",
            Self::NoDenotation => "no denotation",
            Self::Singleton => "fewer than two spellings",
            Self::DenotesTheReference => "denotes the reference",
        }
    }
}

/// One enumerated design: either a row, or the reason it produced none.
pub type Attempt = Result<Row, (String, DropReason)>;

/// A corpus, folded from [`enumerate`].
#[derive(Debug, Clone)]
pub struct ProteinCorpus {
    /// Designs the enumeration proposed, before any filtering.
    pub designs_considered: usize,
    /// Drops, by reason.
    pub drops: BTreeMap<&'static str, usize>,
    /// Every drop, by design id and reason.
    ///
    /// Kept rather than folded away because a drop is where a generator goes
    /// silently blind: `spelling does not parse` counted 45 on the first run
    /// and the *reason* — a `p.Met1Xxx` the parser refuses, and a short-format
    /// frameshift whose residue makes it ambiguous — was only recoverable by
    /// naming the designs. A census that reports a drop count without the ids
    /// cannot answer "which shape stopped being generated".
    pub dropped: Vec<(String, DropReason)>,
    /// The rows.
    pub rows: Vec<Row>,
}

impl ProteinCorpus {
    /// Fold [`enumerate`]'s attempts.
    #[must_use]
    pub fn from_attempts(attempts: Vec<Attempt>) -> Self {
        let designs_considered = attempts.len();
        let mut drops: BTreeMap<&'static str, usize> = BTreeMap::new();
        let mut dropped = Vec::new();
        let mut rows = Vec::new();
        for attempt in attempts {
            match attempt {
                Ok(row) => rows.push(row),
                Err((id, reason)) => {
                    *drops.entry(reason.label()).or_default() += 1;
                    dropped.push((id, reason));
                }
            }
        }
        rows.sort_by(|a, b| a.id.cmp(&b.id));
        dropped.sort();
        Self {
            designs_considered,
            drops,
            dropped,
            rows,
        }
    }

    /// Spellings across every row — the number of normalizations a consumer
    /// will perform.
    #[must_use]
    pub fn spellings(&self) -> usize {
        self.rows.iter().map(|row| row.spellings.len()).sum()
    }

    /// Rows per [`RowKind`].
    #[must_use]
    pub fn by_kind(&self) -> BTreeMap<RowKind, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.kind).or_default() += 1;
        }
        counts
    }

    /// Rows per [`Determinacy`] — the census column
    /// `rulings[confluence-gate-is-apply-equality-on-every-determined-axis]`
    /// requires be kept apart.
    #[must_use]
    pub fn by_determinacy(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.determinacy().label()).or_default() += 1;
        }
        counts
    }

    /// Rows per stratum.
    #[must_use]
    pub fn by_stratum(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.stratum).or_default() += 1;
        }
        counts
    }

    /// Rows per [`Site`].
    #[must_use]
    pub fn by_site(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.site.label()).or_default() += 1;
        }
        counts
    }

    /// Rows per [`Kind`].
    #[must_use]
    pub fn by_edit(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.edit.label()).or_default() += 1;
        }
        counts
    }

    /// Rows per transcript geometry.
    #[must_use]
    pub fn by_shape(&self) -> BTreeMap<&'static str, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            *counts.entry(row.shape.label()).or_default() += 1;
        }
        counts
    }

    /// Rows per junction phase, which is the dimension #1478's protein analogue
    /// lives in.
    #[must_use]
    pub fn by_junction_phase(&self) -> BTreeMap<usize, usize> {
        let mut counts = BTreeMap::new();
        for row in &self.rows {
            // A single-exon frame has no junction, so it has no phase.
            if matches!(row.shape, ProteinRefShape::ThreeExon(_)) {
                *counts.entry(row.design.junction_phase()).or_default() += 1;
            }
        }
        counts
    }
}

/// Enumerate every design, in a fixed order.
///
/// The order is set by the loop nesting and does not depend on hashing, so two
/// runs produce byte-identical output.
#[must_use]
pub fn enumerate() -> Vec<Attempt> {
    let mut attempts = Vec::new();
    for shape in ProteinRefShape::all() {
        for design in all_designs() {
            let frame = frame_for(shape, design);
            enumerate_substitutions(&frame, shape, design, &mut attempts);
            enumerate_run_shifts(&frame, shape, design, &mut attempts);
            enumerate_identities(&frame, shape, design, &mut attempts);
            enumerate_indeterminates(&frame, shape, design, &mut attempts);
            enumerate_alleles(&frame, shape, design, &mut attempts);
            enumerate_conflicts(&frame, shape, design, &mut attempts);
        }
    }
    attempts
}

/// Build the corpus.
#[must_use]
pub fn corpus() -> ProteinCorpus {
    ProteinCorpus::from_attempts(enumerate())
}

// ---------------------------------------------------------------------------
// The strata
// ---------------------------------------------------------------------------

/// The residue at 1-based protein `position`, terminator included, or `None`
/// past the end.
fn residue_at(frame: &ProteinFrame, position: u64) -> Option<AminoAcid> {
    let with_stop = frame.residues_with_stop();
    usize::try_from(position)
        .ok()
        .and_then(|index| index.checked_sub(1))
        .and_then(|index| with_stop.get(index).copied())
}

/// Three-letter code for the residue at `position`.
fn code_at(frame: &ProteinFrame, position: u64) -> Option<&'static str> {
    residue_at(frame, position).map(|aa| aa.to_three_letter())
}

/// A residue that is not the one at `position` and is spellable, so a
/// substitution is never an identity.
fn alternative_to(aa: AminoAcid) -> AminoAcid {
    if aa == AminoAcid::Gly {
        AminoAcid::Ala
    } else {
        AminoAcid::Gly
    }
}

/// Assemble a row, keeping only the spellings the oracle agrees denote the same
/// thing as the authored one.
///
/// The filter is what makes a family's confluence claim mean something: a
/// spelling that does *not* denote the authored protein would make every output
/// difference a property of the corpus rather than of the normalizer. It is also
/// where a generator bug surfaces, because a dropped spelling is counted rather
/// than silently discarded.
#[allow(clippy::too_many_arguments)]
fn build_row(
    frame: &ProteinFrame,
    shape: ProteinRefShape,
    design: CoreDesign,
    stratum: &'static str,
    site: Site,
    edit: Kind,
    tag: &str,
    members: usize,
    rules: Vec<&'static str>,
    spellings: Vec<String>,
) -> Attempt {
    let id = format!(
        "{}-{}-{stratum}-{}-{}-{tag}",
        shape.label(),
        design.label(),
        site.label(),
        edit.label()
    );
    let Some(authored) = spellings.first() else {
        return Err((id, DropReason::Unspellable));
    };
    let denoted = protein_denotation_of(frame.provider(), authored);
    if matches!(
        denoted,
        ProteinDenotation::Unparseable | ProteinDenotation::NotProteinAxis
    ) {
        return Err((id, DropReason::Unspellable));
    }
    let kept: Vec<String> = spellings
        .iter()
        .filter(|spelling| protein_denotation_of(frame.provider(), spelling) == denoted)
        .cloned()
        .collect();
    let kind = match &denoted {
        ProteinDenotation::NoSingleSequence(_) => RowKind::Conflict,
        // The uncertainty stratum's rows are the only ones whose spellings
        // denote one protein while claiming different things about it; see
        // `RowKind::Distinguished`.
        _ if stratum == "uncertainty" && kept.len() > 1 => RowKind::Distinguished,
        _ if kept.len() > 1 => RowKind::Family,
        _ => RowKind::Single,
    };
    Ok(Row {
        id,
        stratum,
        kind,
        shape,
        design,
        site,
        edit,
        members,
        denoted,
        spellings: kept,
        rules,
    })
}

/// A single-member `p.` descriptor.
fn p(suffix: &str) -> String {
    ProteinFrame::protein_descriptor(suffix)
}

/// The substitution stratum: one residue replaced, spelled three ways.
///
/// `protein/delins.md:17` rules that one residue replaced by one is a
/// substitution rather than a `delins`, and `protein/substitution.md:16` that a
/// predicted consequence is parenthesised. Both are rules about **spelling**, so
/// all three spellings denote one protein and the family is a legitimate
/// confluence question rather than a validity one.
fn enumerate_substitutions(
    frame: &ProteinFrame,
    shape: ProteinRefShape,
    design: CoreDesign,
    attempts: &mut Vec<Attempt>,
) {
    for site in Site::all() {
        let position = site.position(frame);
        let Some(from) = residue_at(frame, position) else {
            attempts.push(Err((
                format!(
                    "{}-{}-substitution-{}",
                    shape.label(),
                    design.label(),
                    site.label()
                ),
                DropReason::SiteUnavailable,
            )));
            continue;
        };
        let to = alternative_to(from);
        let (from, to) = (from.to_three_letter(), to.to_three_letter());
        attempts.push(build_row(
            frame,
            shape,
            design,
            "substitution",
            site,
            Kind::Substitution,
            "three-spellings",
            1,
            vec![
                "protein/substitution.md:5-sub-format",
                "protein/delins.md:17-one-for-one-is-a-substitution",
            ],
            vec![
                p(&format!("{from}{position}{to}")),
                p(&format!("{from}{position}_{from}{position}delins{to}")),
            ],
        ));
        // The parenthesised twin is deliberately NOT a member of the family
        // above — see `RowKind::Distinguished`.
        attempts.push(build_row(
            frame,
            shape,
            design,
            "uncertainty",
            site,
            Kind::Uncertainty,
            "parenthesised-twin",
            1,
            vec!["protein/substitution.md:16-predicted-in-parentheses"],
            vec![
                p(&format!("{from}{position}{to}")),
                p(&format!("({from}{position}{to})")),
            ],
        ));
    }
}

/// The ambiguous-tract stratum — the corpus's confluence workhorse.
///
/// Every spelling below denotes one protein **by construction**, because the
/// tract is `RUN_COPIES` copies of one residue: deleting any one of them leaves
/// the same protein, duplicating any one of them leaves the same protein, and
/// inserting the residue at any junction inside the tract leaves the same
/// protein. `general.md:43` makes the 3' rule reach "ALL descriptions (genome,
/// gene, transcript, and protein)", so exactly one of each family is the
/// recommended form and a divergence is a real finding.
fn enumerate_run_shifts(
    frame: &ProteinFrame,
    shape: ProteinRefShape,
    design: CoreDesign,
    attempts: &mut Vec<Attempt>,
) {
    let code = RUN_RESIDUE.to_three_letter();
    let first = RUN_START;
    let last = RUN_START + RUN_COPIES as u64 - 1;

    // Deleting any one copy.
    attempts.push(build_row(
        frame,
        shape,
        design,
        "ambiguous-run",
        Site::AmbiguousRun,
        Kind::Deletion,
        "one-copy",
        1,
        vec![
            "protein/deletion.md:5-del-format",
            "general.md:43-three-prime-rule-reaches-protein",
        ],
        (first..=last)
            .map(|position| p(&format!("{code}{position}del")))
            .collect(),
    ));

    // Duplicating any one copy.
    attempts.push(build_row(
        frame,
        shape,
        design,
        "ambiguous-run",
        Site::AmbiguousRun,
        Kind::Duplication,
        "one-copy",
        1,
        vec![
            "protein/duplication.md:5-dup-format",
            "general.md:43-three-prime-rule-reaches-protein",
        ],
        (first..=last)
            .map(|position| p(&format!("{code}{position}dup")))
            .collect(),
    ));

    // Inserting the tract's own residue at any junction inside it — the same
    // protein a duplication of any copy denotes, spelled as an insertion.
    // `general.md:57`'s "dup outranks ins" analogue is what a normalizer is
    // expected to apply here.
    attempts.push(build_row(
        frame,
        shape,
        design,
        "ambiguous-run",
        Site::AmbiguousRun,
        Kind::Insertion,
        "inside-the-tract",
        1,
        vec![
            "protein/insertion.md:5-ins-flanking-format",
            "protein/duplication.md:18-dup-not-ins",
        ],
        (first..last)
            .map(|position| p(&format!("{code}{position}_{code}{}ins{code}", position + 1)))
            .collect(),
    ));

    // A two-residue deletion, anchored at each equivalent start.
    attempts.push(build_row(
        frame,
        shape,
        design,
        "ambiguous-run",
        Site::AmbiguousRun,
        Kind::Repeat,
        "two-copies",
        1,
        vec![
            "protein/repeated.md:5-repeat-format",
            "protein/deletion.md:5-del-format",
        ],
        (first..last)
            .map(|position| p(&format!("{code}{position}_{code}{}del", position + 1)))
            .chain(std::iter::once(p(&format!(
                "{code}{first}[{}]",
                RUN_COPIES - 2
            ))))
            .collect(),
    ));
}

/// The identity stratum: descriptions that denote the reference protein.
///
/// Kept apart from the others because an identity is where a normalizer's
/// collapse rules meet a description that must not disappear:
/// `protein/substitution.md:36-38` publishes `p.Cys188=` for exactly this.
fn enumerate_identities(
    frame: &ProteinFrame,
    shape: ProteinRefShape,
    design: CoreDesign,
    attempts: &mut Vec<Attempt>,
) {
    for site in [Site::SecondResidue, Site::DeepInterior, Site::LastResidue] {
        let position = site.position(frame);
        let Some(code) = code_at(frame, position) else {
            continue;
        };
        attempts.push(build_row(
            frame,
            shape,
            design,
            "identity",
            site,
            Kind::Identity,
            "three-spellings",
            1,
            vec!["protein/substitution.md:36-silent-equals"],
            vec![
                p(&format!("{code}{position}=")),
                p(&format!("{code}{position}delins{code}")),
                p(&format!("{code}{position}_{code}{position}delins{code}")),
            ],
        ));
    }
}

/// The indeterminate stratum — the rows whose census column exists because of
/// `rulings[confluence-gate-is-apply-equality-on-every-determined-axis]`.
///
/// Each shape here fixes part of the protein and leaves the rest open, and the
/// oracle says exactly which part. They are generated as their own stratum so
/// the count can be reported beside the confluence figures and never folded into
/// them.
fn enumerate_indeterminates(
    frame: &ProteinFrame,
    shape: ProteinRefShape,
    design: CoreDesign,
    attempts: &mut Vec<Attempt>,
) {
    let residues = frame.residues().len() as u64;
    let interior = Site::DeepInterior.position(frame);
    let Some(interior_code) = code_at(frame, interior) else {
        return;
    };
    // `protein/frameshift.md:19`: "the description of a frameshift starts with
    // the **first new** amino acid". A generator that hard-coded `Gly` there
    // emitted `p.Gly36GlyfsTer12` on every design whose residue 36 is already
    // `Gly` — a frameshift whose first new residue is the old one, which ferro
    // refuses at parse. Nine rows vanished into the drop ledger before the ids
    // were read back out of it.
    let new_residue =
        alternative_to(residue_at(frame, interior).unwrap_or(AminoAcid::Gly)).to_three_letter();
    let terminator = residues + 1;

    // A frameshift, long and short form. `protein/frameshift.md:18-23`.
    for (tag, suffix) in [
        (
            "stated-stop",
            format!("{interior_code}{interior}{new_residue}fsTer12"),
        ),
        // `protein/frameshift.md:23`'s short format names the first amino acid
        // changed, its position and `fs` — and no new residue.
        ("short-format", format!("{interior_code}{interior}fs")),
    ] {
        attempts.push(build_row(
            frame,
            shape,
            design,
            "indeterminate",
            Site::DeepInterior,
            Kind::Frameshift,
            tag,
            1,
            vec!["protein/frameshift.md:18-frameshift-format"],
            vec![p(&suffix)],
        ));
    }

    // A C-terminal extension past the terminator (`protein/extension.md:30-34`)
    // and an N-terminal one from upstream (`:22-24`).
    attempts.push(build_row(
        frame,
        shape,
        design,
        "indeterminate",
        Site::Terminator,
        Kind::Extension,
        "c-terminal",
        1,
        vec!["protein/extension.md:30-c-terminal-extension"],
        vec![p(&format!("Ter{terminator}GlnextTer17"))],
    ));
    attempts.push(build_row(
        frame,
        shape,
        design,
        "indeterminate",
        Site::Initiator,
        Kind::Extension,
        "n-terminal",
        1,
        vec!["protein/extension.md:22-n-terminal-extension"],
        vec![p("Met1ext-5")],
    ));

    // A payload stated by count rather than by identity
    // (`protein/insertion.md:45-51`).
    let next = interior + 1;
    if let Some(next_code) = code_at(frame, next) {
        attempts.push(build_row(
            frame,
            shape,
            design,
            "indeterminate",
            Site::DeepInterior,
            Kind::ByCount,
            "unknown-residues",
            1,
            vec!["protein/insertion.md:45-insert-by-count"],
            vec![p(&format!(
                "{interior_code}{interior}_{next_code}{next}insXaa[7]"
            ))],
        ));
        attempts.push(build_row(
            frame,
            shape,
            design,
            "indeterminate",
            Site::DeepInterior,
            Kind::ByCount,
            "to-a-stop",
            1,
            vec!["protein/insertion.md:49-insert-to-a-stop"],
            vec![p(&format!(
                "{interior_code}{interior}_{next_code}{next}ins*9"
            ))],
        ));
    }

    // `p.0` and `p.0?` — no protein at all (`protein/substitution.md:46-48`).
    attempts.push(build_row(
        frame,
        shape,
        design,
        "indeterminate",
        Site::Initiator,
        Kind::NoProtein,
        "stated",
        1,
        vec!["protein/substitution.md:46-no-protein"],
        vec![p("0")],
    ));
    attempts.push(build_row(
        frame,
        shape,
        design,
        "indeterminate",
        Site::Initiator,
        Kind::NoProtein,
        "predicted",
        1,
        vec!["protein/substitution.md:47-no-protein-predicted"],
        vec![p("0?")],
    ));
}

/// The cis-allele stratum: two members on one molecule.
///
/// `protein/delins.md:62` publishes the `p.[Ser44Arg;Trp46Arg]` form. The two
/// members are placed far apart so neither the merge nor the split move is in
/// play, and the row's two spellings differ **only** in member order — which
/// denotes the same protein and must reach one output.
fn enumerate_alleles(
    frame: &ProteinFrame,
    shape: ProteinRefShape,
    design: CoreDesign,
    attempts: &mut Vec<Attempt>,
) {
    let left = Site::SecondResidue.position(frame);
    let right = Site::DeepInterior.position(frame);
    let (Some(left_aa), Some(right_aa)) = (residue_at(frame, left), residue_at(frame, right))
    else {
        return;
    };
    let member = |position: u64, aa: AminoAcid| {
        format!(
            "{}{position}{}",
            aa.to_three_letter(),
            alternative_to(aa).to_three_letter()
        )
    };
    let first = member(left, left_aa);
    let second = member(right, right_aa);
    attempts.push(build_row(
        frame,
        shape,
        design,
        "cis-allele",
        Site::DeepInterior,
        Kind::CisAllele,
        "member-order",
        2,
        vec!["protein/alleles.md:5-cis-allele-format"],
        vec![
            p(&format!("[{first};{second}]")),
            p(&format!("[{second};{first}]")),
        ],
    ));
}

/// The conflict stratum: descriptions that denote no single protein.
///
/// The property is **refusal**, so these rows carry one spelling and the census
/// counts how many the implementation accepts. Each shape cites the clause it
/// violates, and every one of them is a shape the oracle independently reports
/// as [`ProteinDenotation::NoSingleSequence`] — so the corpus and the clause
/// agree before the normalizer is consulted.
fn enumerate_conflicts(
    frame: &ProteinFrame,
    shape: ProteinRefShape,
    design: CoreDesign,
    attempts: &mut Vec<Attempt>,
) {
    let residues = frame.residues().len() as u64;
    let low = Site::AmbiguousRun.position(frame);
    let high = Site::DeepInterior.position(frame);
    let (Some(low_code), Some(high_code)) = (code_at(frame, low), code_at(frame, high)) else {
        return;
    };

    // A range written 3'->5' (`protein/deletion.md:18`).
    attempts.push(build_row(
        frame,
        shape,
        design,
        "conflict",
        Site::DeepInterior,
        Kind::Deletion,
        "reversed-range",
        1,
        vec!["protein/deletion.md:18-list-five-to-three"],
        vec![p(&format!("{high_code}{high}_{low_code}{low}del"))],
    ));

    // An insertion whose two positions do not flank one junction
    // (`protein/insertion.md:17-18`).
    let far = low + 2;
    if let Some(far_code) = code_at(frame, far) {
        attempts.push(build_row(
            frame,
            shape,
            design,
            "conflict",
            Site::AmbiguousRun,
            Kind::Insertion,
            "non-flanking",
            1,
            vec!["protein/insertion.md:17-two-flanking-residues"],
            vec![p(&format!("{low_code}{low}_{far_code}{far}insGly"))],
        ));
    }

    // Two members claiming one residue — the protein twin of the nucleotide
    // corpus's overlapping geometry.
    attempts.push(build_row(
        frame,
        shape,
        design,
        "conflict",
        Site::DeepInterior,
        Kind::CisAllele,
        "overlapping-members",
        2,
        vec!["protein/alleles.md:5-cis-allele-format"],
        // BOTH members must differ from the reference residue, or one of them
        // is a silent identity and the row is a weaker conflict than it reads
        // as. `Gly`/`Ala` are the two `alternative_to` ever returns, so the
        // reference residue is excluded by construction.
        vec![p(&format!(
            "[{high_code}{high}{first};{high_code}{high}{second}]",
            first =
                alternative_to(residue_at(frame, high).unwrap_or(AminoAcid::Gly)).to_three_letter(),
            second = if residue_at(frame, high) == Some(AminoAcid::Ser) {
                AminoAcid::Trp.to_three_letter()
            } else {
                AminoAcid::Ser.to_three_letter()
            },
        ))],
    ));

    // A position past the terminator.
    attempts.push(build_row(
        frame,
        shape,
        design,
        "conflict",
        Site::PastEnd,
        Kind::Substitution,
        "past-the-terminator",
        1,
        vec!["protein/substitution.md:5-sub-format"],
        vec![p(&format!("Gly{}Ala", residues + 5))],
    ));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::spec_corpus::corpus_cores;
    use crate::reference::transcript::Strand;

    #[test]
    fn a_designed_core_is_an_open_reading_frame_that_round_trips() {
        for design in all_designs() {
            let cds = codon_core(design);
            assert_eq!(cds.len(), design.body_codons * 3, "{}", design.label());
            assert!(cds.starts_with("ATG"), "{}", design.label());
            let frame = frame_for(ProteinRefShape::SingleExon, design);
            assert_eq!(
                frame.residues().len(),
                design.residues(),
                "{} translated to the wrong length",
                design.label()
            );
            assert_eq!(frame.residues()[0], AminoAcid::Met, "{}", design.label());
        }
    }

    /// The blocker this module exists to remove, demonstrated rather than
    /// asserted: `spec_corpus`'s cores do not translate, so a protein stratum
    /// could not have been drawn from them.
    #[test]
    fn the_nucleotide_corpuss_cores_are_not_reading_frames() {
        let table = CodonTable::standard();
        let mut translatable = 0usize;
        let cores = corpus_cores(16, 96);
        for core in &cores {
            let starts = core.starts_with("ATG");
            let whole_codons = core.len().is_multiple_of(3);
            let no_internal_stop = (0..core.len() / 3).all(|index| {
                Codon::parse(&core[index * 3..index * 3 + 3])
                    .is_some_and(|codon| !table.is_stop(&codon))
            });
            if starts && whole_codons && no_internal_stop {
                translatable += 1;
            }
        }
        assert_eq!(
            translatable,
            0,
            "{translatable} of {} nucleotide-corpus cores are open reading frames; if this is \
             no longer zero the protein stratum's premise has changed",
            cores.len()
        );
    }

    /// The designed tract is ambiguous on the protein axis and holds **no**
    /// repeat at all on the nucleotide axis, because each copy is spelled with a
    /// different synonymous codon.
    ///
    /// This is the one property a nucleotide corpus structurally cannot carry,
    /// and it is checked rather than described because "the codons differ" is
    /// exactly the kind of claim that survives a refactor as a stale comment.
    #[test]
    fn the_ambiguous_run_is_invisible_to_the_nucleotide_axis() {
        for design in all_designs() {
            let frame = frame_for(ProteinRefShape::SingleExon, design);
            let start = usize::try_from(RUN_START).expect("fits");
            let residues = frame.residues();
            for offset in 0..RUN_COPIES {
                assert_eq!(
                    residues[start - 1 + offset],
                    RUN_RESIDUE,
                    "{} residue {} is not the tract residue",
                    design.label(),
                    start + offset
                );
            }
            // The tract is bounded: neither flank is the tract residue, so the
            // ambiguity has exactly the extent the corpus believes it has.
            assert_ne!(residues[start - 2], RUN_RESIDUE, "{}", design.label());
            assert_ne!(
                residues[start - 1 + RUN_COPIES],
                RUN_RESIDUE,
                "{}",
                design.label()
            );

            // And the codons are all different, so the CDS holds no tandem
            // repeat over the tract.
            let cds = frame.cds();
            let codons: Vec<&str> = (0..RUN_COPIES)
                .map(|offset| {
                    let at = (start - 1 + offset) * 3;
                    &cds[at..at + 3]
                })
                .collect();
            let mut distinct = codons.clone();
            distinct.sort_unstable();
            distinct.dedup();
            assert_eq!(
                distinct.len(),
                RUN_COPIES,
                "{}: the tract's codons are not all distinct: {codons:?}",
                design.label()
            );
        }
    }

    /// Both ways a codon can straddle a splice junction are reached, **and** so
    /// is the codon-aligned case that is the control.
    ///
    /// The phase is derived from the frame rather than from
    /// [`CoreDesign::junction_phase`]'s arithmetic, so the two are independent
    /// and a mistake in the arithmetic fails here rather than being confirmed by
    /// itself.
    #[test]
    fn both_junction_phases_are_reached_and_so_is_neither() {
        let mut seen: BTreeMap<usize, usize> = BTreeMap::new();
        for design in all_designs() {
            let frame = frame_for(ProteinRefShape::ThreeExon(Strand::Plus), design);
            let observed: Vec<usize> = (1..=frame.residues().len() as u64)
                .filter_map(|residue| frame.junction_phase_of_residue(residue))
                .collect();
            let phase = observed.first().copied().unwrap_or(0);
            assert_eq!(
                phase,
                design.junction_phase(),
                "{}: measured phase {phase}, `junction_phase()` says {}",
                design.label(),
                design.junction_phase()
            );
            *seen.entry(design.junction_phase()).or_default() += 1;
        }
        assert_eq!(
            seen.keys().copied().collect::<Vec<_>>(),
            vec![0, 1, 2],
            "the designs do not reach all three junction phases: {seen:?}"
        );
    }

    /// The designed tract sits past `CANONICAL_PAD`, so a frame built here is
    /// outside the regime `codon_frame_deep_window_offset.rs` shows is blind.
    ///
    /// Checked rather than commented, because "deep enough" is a claim about
    /// three constants that can each move independently.
    #[test]
    fn the_tract_clears_the_canonical_window_clamp() {
        /// `merge.rs`'s private window pad, mirrored for the assertion only —
        /// exactly as `codon_frame_deep_window_offset.rs` mirrors it.
        const CANONICAL_PAD: usize = 128;
        for design in all_designs() {
            let frame = frame_for(ProteinRefShape::SingleExon, design);
            let coding = frame
                .cds_position_of_residue(RUN_START)
                .expect("the tract's first residue has a codon");
            let (cds_start, _) = frame.cds_bounds();
            let transcript_position = cds_start + usize::try_from(coding).expect("fits") - 1;
            assert!(
                transcript_position > CANONICAL_PAD,
                "{}: the tract starts at transcript position {transcript_position}, inside \
                 the clamped window (CANONICAL_PAD = {CANONICAL_PAD})",
                design.label()
            );
        }
    }

    #[test]
    fn every_stop_codon_is_reached() {
        let table = CodonTable::standard();
        let stops: Vec<String> = table
            .stop_codons()
            .iter()
            .map(ToString::to_string)
            .collect();
        let mut seen: Vec<String> = all_designs()
            .into_iter()
            .map(|design| {
                let cds = codon_core(design);
                cds[cds.len() - 3..].to_string()
            })
            .collect();
        seen.sort();
        seen.dedup();
        assert_eq!(
            seen.len(),
            stops.len(),
            "stops reached: {seen:?} of {stops:?}"
        );
    }

    /// The filler draws widely from the codon table, so the corpus is not a
    /// measurement of twenty codons wearing a disguise.
    #[test]
    fn the_generator_reaches_most_of_the_codon_table() {
        let mut seen: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
        for design in all_designs() {
            let cds = codon_core(design);
            for index in 0..cds.len() / 3 {
                seen.insert(cds[index * 3..index * 3 + 3].to_string());
            }
        }
        assert!(
            seen.len() >= 40,
            "only {} distinct codons are reachable, which is close enough to \
             back-translation's 20 that the codon design buys nothing: {seen:?}",
            seen.len()
        );
    }

    #[test]
    fn every_site_class_is_reached() {
        let built = corpus();
        let by_site = built.by_site();
        for site in Site::all() {
            assert!(
                by_site.get(site.label()).copied().unwrap_or(0) > 0,
                "no row at site {}: {by_site:?}",
                site.label()
            );
        }
    }

    #[test]
    fn every_shape_is_reached() {
        let built = corpus();
        let by_shape = built.by_shape();
        for shape in ProteinRefShape::all() {
            assert!(
                by_shape.get(shape.label()).copied().unwrap_or(0) > 0,
                "no row on shape {}: {by_shape:?}",
                shape.label()
            );
        }
        let by_phase = built.by_junction_phase();
        assert_eq!(
            by_phase.keys().copied().collect::<Vec<_>>(),
            vec![0, 1, 2],
            "the corpus's multi-exon rows do not span every junction phase: {by_phase:?}"
        );
    }

    #[test]
    fn every_determinacy_class_is_reached() {
        let built = corpus();
        let by_determinacy = built.by_determinacy();
        for class in [
            Determinacy::Determined,
            Determinacy::Indeterminate,
            Determinacy::NoProtein,
            Determinacy::Undenoted,
        ] {
            assert!(
                by_determinacy.get(class.label()).copied().unwrap_or(0) > 0,
                "no row in determinacy class {}: {by_determinacy:?}",
                class.label()
            );
        }
    }

    #[test]
    fn every_kind_is_reached() {
        let built = corpus();
        let by_edit = built.by_edit();
        for edit in [
            Kind::Substitution,
            Kind::Deletion,
            Kind::Duplication,
            Kind::Insertion,
            Kind::Repeat,
            Kind::Identity,
            Kind::Frameshift,
            Kind::Extension,
            Kind::ByCount,
            Kind::NoProtein,
            Kind::CisAllele,
            Kind::Uncertainty,
        ] {
            assert!(
                by_edit.get(edit.label()).copied().unwrap_or(0) > 0,
                "no row of kind {}: {by_edit:?}",
                edit.label()
            );
        }
        assert!(
            built.rows.iter().any(|row| row.members == 2),
            "no two-member row"
        );
    }

    /// Every family's spellings really do denote one protein, checked through
    /// the oracle rather than trusted from the generator.
    #[test]
    fn every_family_denotes_one_protein() {
        let built = corpus();
        let mut families = 0usize;
        for row in &built.rows {
            if row.kind != RowKind::Family {
                continue;
            }
            families += 1;
            let frame = row.frame();
            for spelling in &row.spellings {
                assert_eq!(
                    protein_denotation_of(frame.provider(), spelling),
                    row.denoted,
                    "{}: {spelling} does not denote the row's protein",
                    row.id
                );
            }
        }
        assert!(families > 0, "the corpus holds no families at all");
    }

    /// A row id names its design uniquely, so a failure can be pinned.
    #[test]
    fn row_ids_are_unique() {
        let built = corpus();
        let mut ids: Vec<&str> = built.rows.iter().map(|row| row.id.as_str()).collect();
        let total = ids.len();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), total, "duplicate row ids");
    }

    /// Enumeration is deterministic.
    #[test]
    fn the_corpus_is_deterministic() {
        let first = corpus();
        let second = corpus();
        assert_eq!(first.rows.len(), second.rows.len());
        for (a, b) in first.rows.iter().zip(second.rows.iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.spellings, b.spellings);
        }
    }
}
