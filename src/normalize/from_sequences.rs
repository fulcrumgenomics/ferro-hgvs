//! Derive an HGVS description from a reference/alternate sequence pair.
//!
//! The caller supplies the bases; this module supplies the description. It reads
//! **no reference sequence**, so its output is a pure function of its inputs —
//! which is what makes it deterministic in the sense a BAM post-processor needs:
//! the same bases give the same description on any machine, against any
//! reference build, with no hidden input.
//!
//! "Its inputs" is five values, not four, and the fifth is not inert: the
//! accession, the position, the two sequences **and** the
//! [`FromSequencesOptions`] — whose `direction` moves a placement within the
//! window and whose `max_grid_cells` decides whether an answer is produced at
//! all. Purity is the claim; a four-argument count was never accurate and is
//! withdrawn wherever it was written.
//!
//! # Case is folded, and the axis is DNA
//!
//! Both sequences are upper-cased before anything reads them. A soft-masked
//! reference arrives lower-case (`fetch_window` does not case-fold), while a
//! rendered payload is always upper-case, so without folding a masked window
//! aligns against its own payload as if the two disagreed — and the round trip
//! below then refuses a derivation that was correct. Folding once, at the
//! entry, means every later stage sees one alphabet.
//!
//! `U` is refused rather than folded to `T`. This surface emits `g.`/`m.`
//! descriptions, which are DNA; admitting `U` produced `NC_TEST.1:g.8U>T`,
//! a well-formed string naming a base its own axis does not have.
//!
//! # What it owns, and what it does not
//!
//! `README.md`'s four normalization rules split cleanly here. This module
//! delivers rules **1 (conformant)** and **4 (deterministic)** — the two the
//! README calls always achievable. Rules **2 (recommended form)** and
//! **3 (confluent)** stay with [`crate::Normalizer::normalize`], because both
//! need the reference: rule 2's scope names the 3' rule explicitly, and a
//! reference-anchored shift is precisely what a window-local function cannot do.
//!
//! So an output here may be 3'-shiftable further than this module shifted it,
//! and that is not a defect. Run `normalize` afterwards if you want it.
//!
//! # Why it partitions at all
//!
//! Strictly it need not: a single spanning `delins` is a legal description, so
//! this module could emit one and leave any split to `normalize`. It cannot
//! today — handed a spanning `delins`, the shipped normalizer leaves it alone,
//! because [`super::merge::partition_block`] searches single-gap alignments only
//! and the weight bound then refuses the re-derivation. Measured over two 1 bp
//! deletions swept across separations 0-20 on `NC_000001.11:g.1000050`, 17 of
//! the 21 separations come back as one spanning `delins`.
//!
//! **That is a workaround and not the design.** When #1419 / #1420 / #1421 /
//! #1440 land, `normalize` will re-derive over a spanning `delins` and this
//! reason expires. The partitioning stays, because a partition derived from the
//! sequence is what this module is for; the justification above should not be
//! mistaken for a permanent one.
//!
//! # Why the alignment DAG, and not an aligner
//!
//! [`crate::normalize::seqfirst::align::AlignmentDag`] holds **every**
//! minimal-cost alignment at once, and `CanonicalAlignment` takes the
//! member-count-minimal one, 3'-most among ties. The tie-break is
//! `general.md:41`, an explicit spec tie-break, which the
//! `canonical-form-choice-when-both-legal` ruling requires be applied.
//!
//! Affine-gap Needleman-Wunsch was implemented and measured against it over
//! 28,639 synthetic two-event blocks at separations 0-12. bwa-mem's defaults
//! (`A=1 B=4 O=6 E=1`) return an alignment above the block's edit distance on
//! **6.7%** of them, at up to double the necessary changed columns, because a
//! mismatch costs 4 where a gap-extend costs 1 — so it prefers deleting and
//! reinserting a run to substituting it. minimap2's dual-affine `lr:hq` scheme
//! reaches 2.2%. `CanonicalAlignment` is 0 of 28,639, not because it is better
//! tuned but because every path it walks is distance-minimal by construction.
//! Both the scoring matrix and the traceback order are free parameters and both
//! are visible in the output; this has neither.
//!
//! **Distance-minimality is ferro's policy, not the spec's.** `basics.md:38`
//! lists stability, meaning, memorability and unequivocality and does not
//! mention minimality at all; `DNA/delins.md:44-47` recommends a *non-minimal*
//! description in its own worked example. Cite this module's own reasoning for
//! it, never the recommendations.

use crate::error::FerroError;
use crate::hgvs::edit::{InsertedSequence, NaEdit};
use crate::hgvs::interval::UncertainBoundary;
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{AllelePhase, AlleleVariant, HgvsVariant};
use crate::normalize::merge::{
    denoted_bases, derive_block_members, BlockDecline, MAX_SEQFIRST_GRID_CELLS,
};
use crate::normalize::ShuffleDirection;

/// Cost knobs for [`from_sequences`].
///
/// Nothing here changes which *forms* the function is willing to emit — that is
/// fixed by the rules, not by the caller (`README.md` rule 6). `max_grid_cells`
/// is a memory bound. `direction` is **not** a caller-facing knob: it mirrors
/// `NormalizeConfig`'s internal test instrument, is `#[doc(hidden)]` for the
/// same reason, and is always `ThreePrime` on every shipped path.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct FromSequencesOptions {
    /// Largest alignment grid, in cells, the partitioner will build.
    ///
    /// A cell costs roughly **18 bytes** on this arm, so the default —
    /// `(4096 + 1)^2`, about 16.8 M cells — admits roughly **310 MB** at the
    /// limit. Lower it when reads are short and the budget matters; raise it for
    /// a long-read haplotype, having read that figure. Exceeding it **refuses**:
    /// this is a cost bound, and refusing rather than silently answering with a
    /// weaker rule is the policy.
    pub max_grid_cells: usize,
    /// Which end of an ambiguous run a pure indel is placed at, within the
    /// caller's window. Always `ThreePrime`, which matches `general.md:41`.
    ///
    /// **Internal test instrument, not a supported knob** — see
    /// [`crate::normalize::ShuffleDirection`].
    #[doc(hidden)]
    pub direction: ShuffleDirection,
}

impl Default for FromSequencesOptions {
    fn default() -> Self {
        Self {
            max_grid_cells: MAX_SEQFIRST_GRID_CELLS,
            direction: ShuffleDirection::ThreePrime,
        }
    }
}

/// The two setters exist because the struct is `#[non_exhaustive]`, which
/// forbids a struct expression **outside this crate** — so without them a
/// downstream caller can reach `Default::default()` and nothing else, and the
/// knobs are documented but unreachable. The bug was found by moving one test
/// into `tests/it`, which is a separate crate and so sees the same surface a
/// user does.
///
/// Named after [`crate::NormalizeConfig::with_direction`], which is the same
/// pattern for the same reason.
impl FromSequencesOptions {
    /// Set which end of an ambiguous run a pure indel is placed at.
    ///
    /// **Internal test instrument, not a supported knob** — see
    /// [`crate::normalize::ShuffleDirection`].
    #[doc(hidden)]
    #[must_use]
    pub fn with_direction(mut self, direction: ShuffleDirection) -> Self {
        self.direction = direction;
        self
    }

    /// Set the alignment-grid budget, in cells. See
    /// [`FromSequencesOptions::max_grid_cells`] for what a cell costs.
    #[must_use]
    pub fn with_max_grid_cells(mut self, max_grid_cells: usize) -> Self {
        self.max_grid_cells = max_grid_cells;
        self
    }
}

/// A derived description, plus the one caveat a window-local derivation owes its
/// caller.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct DerivedDescription {
    /// The description.
    pub variant: HgvsVariant,
    /// Whether a member rests on the window's **5' edge**.
    ///
    /// Ask [`Self::placement_bounded_by_window`] for the plain "could this move
    /// at all" answer; this and [`Self::bounded_at_end`] say *which* side, which
    /// is what a caller widening one side at a time needs.
    pub bounded_at_start: bool,
    /// Whether a member rests on the window's **3' edge**.
    ///
    /// The 3' counterpart of [`Self::bounded_at_start`]; see it.
    pub bounded_at_end: bool,
}

impl DerivedDescription {
    /// Whether a member's placement was bounded by the window rather than
    /// settled by the sequence — i.e. it rests on *either* edge of the bases
    /// supplied. The OR of [`Self::bounded_at_start`] and
    /// [`Self::bounded_at_end`].
    ///
    /// **This is a "could move" flag, not a "is wrong" flag**, and it is
    /// deliberately conservative in that direction. Distinguishing a placement
    /// that merely *reaches* the edge from one that would have gone *past* it
    /// requires knowing what lies outside the window — the reference — which
    /// this function does not read. So it reports the uncertainty rather than
    /// resolving it.
    ///
    /// And note what a flagged answer is *not*: it is never wrong. `g.14del`
    /// and `g.15del` over that run denote the same bases and share a canonical
    /// SPDI. What a clipped placement costs is the **recommended form** (rule 2)
    /// and **agreement with a wider read** (rule 3) — the two rules this
    /// function never claimed, because both need the reference. Rules 1 and 4
    /// hold regardless.
    ///
    /// The condition under which the flag stops mattering is exact: two windows
    /// that both contain the whole interval over which the change can be placed
    /// derive the same description, and a window that cuts that interval places
    /// the change at its own edge instead.
    #[must_use]
    pub fn placement_bounded_by_window(&self) -> bool {
        self.bounded_at_start || self.bounded_at_end
    }
}

/// Derive an HGVS description from a reference/alternate sequence pair.
///
/// `position` is **1-based** and names the first base of `reference`, matching
/// `VcfRecord` and HGVS.
///
/// `reference` is taken on trust to be the reference bases over
/// `[position, position + reference.len())`. It is **not** verified, because
/// verifying it would need the reference and would make the provider a hidden
/// fifth input — which would cost exactly the determinism this function exists
/// to provide. Pass bases that are not the reference and you get a faithful
/// description of the pair you passed.
///
/// See [`crate::Normalizer::from_sequences`] for the sibling that does hold a
/// provider, and so can additionally refuse an unknown accession or an interval
/// running past the end of the sequence.
///
/// # Errors
///
/// Refuses a zero `position`, an empty `reference`, a non-nucleotide symbol in
/// either sequence, an alignment grid over `options.max_grid_cells`, a partition
/// that will not render, and — as a runtime check, never a debug assertion — a
/// derivation that does not re-apply to `alternate`.
pub fn from_sequences(
    accession: &str,
    position: u64,
    reference: &str,
    alternate: &str,
    options: &FromSequencesOptions,
) -> Result<HgvsVariant, FerroError> {
    from_sequences_detailed(accession, position, reference, alternate, options)
        .map(|derived| derived.variant)
}

/// [`from_sequences`], reporting also whether the derivation reached a window
/// edge. See [`DerivedDescription::placement_bounded_by_window`].
pub fn from_sequences_detailed(
    accession: &str,
    position: u64,
    reference: &str,
    alternate: &str,
    options: &FromSequencesOptions,
) -> Result<DerivedDescription, FerroError> {
    validate(position, reference.as_bytes(), alternate.as_bytes())?;
    // Folded **after** validation, so a refusal quotes the symbol the caller
    // actually passed rather than an upper-cased rewrite of it. See the module
    // docs for why folding happens at all.
    let (reference, alternate) = (
        reference.to_ascii_uppercase().into_bytes(),
        alternate.to_ascii_uppercase().into_bytes(),
    );
    let (reference, alternate) = (reference.as_slice(), alternate.as_slice());

    let template = reference_template(accession)?;
    let w_lo = i64::try_from(position).map_err(|_| FerroError::InvalidCoordinates {
        msg: format!("position {position} does not fit a signed coordinate"),
    })?;

    // This surface's block partitioner is a **pin**, not `FERRO_PARTITION`'s
    // reading: `merge::DERIVED_BLOCK_PARTITION_RULE` names it, and
    // `merge::partition_block_for_derivation` is where it is consulted. So a
    // change to the rule `Normalizer::normalize` cuts with does not reach here,
    // which is #1834 — the tripwire that fails when either moves alone is
    // `merge`'s `the_two_surfaces_cut_with_a_pinned_pair_of_rules`.
    let block = derive_block_members(
        reference,
        alternate,
        &template,
        w_lo,
        options.direction,
        options.max_grid_cells,
    )
    .map_err(|decline| match decline {
        BlockDecline::GridTooLarge { ref_len, alt_len } => {
            grid_refusal(ref_len, alt_len, options.max_grid_cells)
        }
        BlockDecline::WouldNotRender => FerroError::ConversionError {
            msg: format!(
                "could not render the derived partition of {accession}:{position} as HGVS \
                 members"
            ),
        },
    })?;

    if block.anchors_before_window {
        return Err(five_prime_anchor_refusal(accession, position));
    }

    let mut members = block.members;
    retype_inversions(&mut members, reference, w_lo);

    let variant = match members.len() {
        0 => template,
        1 => members.into_iter().next().expect("length checked"),
        _ => HgvsVariant::Allele(AlleleVariant::new(members, AllelePhase::Cis)),
    };

    verify_round_trip(&variant, w_lo, reference, alternate)?;

    Ok(DerivedDescription {
        variant,
        bounded_at_start: block.bounded_at_start,
        bounded_at_end: block.bounded_at_end,
    })
}

/// The refusal for a pure insertion resting on the window's 5' edge.
///
/// See [`crate::normalize::merge::DerivedBlock::anchors_before_window`] for why
/// this cannot be answered rather than refused: the only HGVS spelling of such a
/// payload names `position - 1`, which the caller supplied no bases for and
/// which does not exist at all when `position` is 1.
///
/// Named as its own refusal rather than left to `verify_round_trip`, which
/// catches it as a side effect and reports it as an internal invariant failure
/// ("could not be re-applied to its window") — true, unactionable, and a
/// misattribution of a caller-fixable problem to ferro. Found by running the
/// #1419/#1420/#1421 reported pairs at a zero pad, where `g.[37dup;41del]`
/// derives to `g.[37_38insA;41del]` against a window starting at 38.
///
/// **How much flank is enough depends on the direction**, which is why the
/// message says "more" rather than naming a number. One base suffices for the
/// row above under `ThreePrime`. Under `FivePrime` the same payload is placed
/// 5'-most, so it walks to the start of whatever ambiguous run it sits in and
/// needs enough flank to clear it — the same row still refuses at a pad of one
/// and derives at six. Any fixed figure here would be wrong for one of the two
/// directions.
fn five_prime_anchor_refusal(accession: &str, position: u64) -> FerroError {
    FerroError::InvalidCoordinates {
        msg: format!(
            "{accession}:{position} — the derivation places an inserted payload immediately \
             5' of the window's first base, so its only HGVS anchor is position {} (an \
             insertion is written between the two positions it falls between), which is \
             outside the window supplied. Supply more 5' flank — \
             `Normalizer::to_sequences` pads both sides, and how much is enough depends on \
             the shuffle direction, since a 5'-most placement walks to the start of an \
             ambiguous run",
            position.saturating_sub(1)
        ),
    }
}

/// Reject the inputs no description can be derived from.
///
/// Deliberately short. Range-checking `position` against the sequence needs the
/// reference, which this path does not hold; that check belongs to
/// [`crate::Normalizer::from_sequences`] and is documented as living there.
///
/// Shared with [`crate::SequencePair::new`] rather than copied into it. The two
/// entry points accept the same four values and must accept exactly the same
/// ones — a constructor that admitted a pair `from_sequences` would later refuse
/// would just move the error somewhere less useful.
///
/// Runs on the caller's own bytes, before the upper-casing
/// [`from_sequences_detailed`] applies, so a message quotes what was passed.
/// Case itself is not a reason to refuse: `Base::from_char` folds, and a
/// soft-masked window is ordinary input.
pub(crate) fn validate(
    position: u64,
    reference: &[u8],
    alternate: &[u8],
) -> Result<(), FerroError> {
    if position == 0 {
        return Err(FerroError::InvalidCoordinates {
            msg: "position is 1-based; 0 does not name a base".to_string(),
        });
    }
    if reference.is_empty() {
        return Err(FerroError::InvalidCoordinates {
            msg: "reference is empty, so there is no interval to describe".to_string(),
        });
    }
    // The alphabet is `Base::from_char`'s — the codebase's own IUPAC-IUBMB
    // table — rather than a list local to this function.
    //
    // It was `ACGTN` at first, which is stricter than the spec and stricter
    // than real data. `general.md:48` admits the IUPAC-IUBMB symbol set, and
    // the `alignment-only-symbol-in-a-description` ruling cites `:48` for that
    // scope while excluding only `X` and `-`; ambiguity codes are in. Measured
    // against the harvested ClinVar/CMRG/Paraphase multi-member alleles, the
    // narrow alphabet refused real submitted rows outright —
    // `NM_000518.4:c.[20A>T;249G>Y]` among them.
    //
    // `X` and `-` are still refused, because `from_char` does not admit them.
    //
    // `U` is refused too, and that exclusion is this surface's own rather than
    // `from_char`'s. `from_char` admits `U` because it also serves the `r.`
    // axis; this surface emits `g.` and `m.` descriptions, which are DNA. It was
    // admitted asymmetrically before: `ACGU` -> `ACGT` derived
    // `NC_TEST.1:g.8U>T`, a well-formed string naming a base the axis does not
    // have, while `ACGT` -> `AUGT` refused with the round trip's internal
    // "denotes different bases" message — a caller-fixable input reported as an
    // invariant failure. One refusal, stated at the argument, replaces both.
    for (label, bases) in [("reference", reference), ("alternate", alternate)] {
        if let Some(bad) = bases
            .iter()
            .find(|b| crate::hgvs::edit::Base::from_char(**b as char).is_none())
        {
            return Err(FerroError::ConversionError {
                msg: format!(
                    "{label} contains '{}', which is not a nucleotide — standards.md:39 admits \
                     no alignment-only symbol in a description; pass IUPAC-IUBMB nucleotide \
                     codes (general.md:48)",
                    *bad as char
                ),
            });
        }
        if bases.iter().any(|b| b.eq_ignore_ascii_case(&b'U')) {
            return Err(FerroError::ConversionError {
                msg: format!(
                    "{label} contains 'U', which is RNA; this surface derives g./m. \
                     descriptions, whose axis is DNA. Pass 'T' instead, or project onto an r. \
                     axis afterwards"
                ),
            });
        }
    }
    Ok(())
}

/// The refusal for a block whose alignment grid exceeds the budget.
///
/// Names the knob and its per-cell cost, so the caller can decide rather than
/// guess. A refusal that does not say how to proceed is a dead end.
fn grid_refusal(ref_len: usize, alt_len: usize, budget: usize) -> FerroError {
    FerroError::ConversionError {
        msg: format!(
            "alignment grid for a {ref_len} x {alt_len} window exceeds max_grid_cells \
             ({budget}); raise FromSequencesOptions::max_grid_cells — about 18 bytes per cell — \
             or narrow the window"
        ),
    }
}

/// A minimal `g.` (or `m.`) variant carrying only the accession, for
/// [`super::merge::rebuild_members`]' template argument and as the identity
/// description.
///
/// Built by parsing rather than by constructing the types directly: `g.=` is the
/// identity description this function must return for an unchanged pair anyway,
/// and `build_merged` reads only the accession and gene symbol off a template,
/// so one parse serves both jobs and neither depends on internal constructors.
///
/// # Why the accession is classified before the template is built
///
/// The classification used to run *after* `parse_hgvs("{accession}:g.=")`, which
/// made it unreachable for exactly the accessions that most need it: the parser
/// refuses a `g.` axis on an `NR_`/`XR_` outright (`Accession::is_noncoding_rna`
/// — #486), so `NR_000001.1` fell out of the parse arm with the generic
/// "'NR_000001.1' is not a usable genomic accession", which reads as "ferro does
/// not know this accession" rather than "this is a non-coding transcript, use
/// `n.`". Classifying first gives every non-genomic class the same informative
/// refusal.
///
/// # The axis is `m.` on a mitochondrial accession
///
/// `Accession::is_mitochondrial` is the codebase's own predicate, and its doc
/// states the rule: "HGVS requires the `m.` coordinate system for these
/// accessions, so `normalize()` coerces a `g.` variant on one of them to `m.`."
/// This surface emitted `NC_012920.1:g.6C>T` for the rCRS — a rule-1 conformance
/// defect on the most-used mitochondrial accession, in a module that claims rule
/// 1 as one of the two it delivers. `MtVariant` carries a `GenomeInterval`, so
/// the template is the only thing that has to change; `build_merged` already
/// dispatches on `HgvsVariant::Mt`.
fn reference_template(accession: &str) -> Result<HgvsVariant, FerroError> {
    // A `g.` description on a transcript or protein accession is not a
    // description the spec admits — `checklist.md:20` is explicit that an `NM_`
    // needs a genomic reference to carry genomic coordinates at all — and the
    // parser will happily build one, because the prefix and the coordinate axis
    // are independent to it.
    //
    // Found by running the cis confluence corpus through this function: its
    // `c.`-axis classes are drawn against `NM_TEST.1`, and every one of them
    // came back as `NM_TEST.1:g.<n>del`. That is a well-formed string denoting
    // nothing, which is the worst shape an output can take.
    //
    // Classified through `Accession::inferred_variant_type`, the codebase's own
    // accession-class table, rather than a prefix list local to this function.
    // A local list was the first attempt and it leaked exactly the shapes a
    // second table always leaks: it knew `ENST` but not `ENSMUST`/`ENSRNOT`
    // (Ensembl feature letters are species-independent, #1057), and knew
    // nothing of LRG's `t<M>`/`p<M>` discriminator or of UniProt — so
    // `ENSMUST00000123.1:g.12_13del` and `LRG_1p1:g.12_13del` were both emitted.
    //
    // Parse first, then ask the parsed accession: `parse_hgvs` is needed anyway
    // for the template, and `HgvsVariant::accession` hands the classified
    // `Accession` straight back.
    // Parsed on its own, not read back off a `g.=` template, so the
    // classification below is reachable for the accessions the parser refuses a
    // `g.` axis on. `parse_accession` is nom-shaped, so the leftover must be
    // checked: a partial parse would classify a prefix the caller did not write.
    // Checked before the parse, not after: `parse_accession` does not consume a
    // bare `P12345` cleanly, so a check placed downstream of it never ran and
    // the accession fell out with the generic "not a usable genomic accession" —
    // technically a refusal, but one that reads as "ferro does not know this
    // accession" rather than "this is a protein".
    if is_uniprot_shaped(accession) {
        return Err(non_genomic_refusal(accession));
    }
    // Parsed with the `:g.=` suffix attached, which is how the parser is
    // designed to see an accession: `parse_simple_accession` locates the HGVS
    // separator to know where the accession ends, so a **bare** SAM refname
    // (`chr1`, `scaffold_123`, `my-contig.v2`) is not fully consumed on its own.
    // Classifying off the bare string refused every one of those, which the
    // shipped surface accepts — `inferred_variant_type` returns `None` for an
    // unclassifiable accession and `None` is not a refusal.
    //
    // The suffix is required back verbatim rather than merely allowed, so a
    // partial parse cannot classify a prefix the caller did not write.
    let probe = format!("{accession}:g.=");
    let parsed = match crate::hgvs::parser::accession::parse_accession(&probe) {
        Ok((":g.=", parsed)) => parsed,
        _ => {
            return Err(FerroError::InvalidCoordinates {
                msg: format!("'{accession}' is not a usable genomic accession"),
            })
        }
    };

    // `inferred_variant_type` covers the classes it knows; the model-prediction
    // (`XM`/`XR`/`XP`) and RefSeq-protein (`YP`/`AP`) prefixes are not in that
    // table, so they are named here explicitly rather than assumed absent. That
    // asymmetry is the argument for widening the shared table one day, not for
    // keeping a second one.
    const NON_GENOMIC_PREFIXES: [&str; 5] = ["XM", "XR", "XP", "YP", "AP"];
    let prefix = parsed.prefix.to_string();
    let inferred = parsed.inferred_variant_type();
    // `is_uniprot` is checked explicitly because `inferred_variant_type`'s
    // UniProt arm keys on a **one-character** prefix, and this parser hands back
    // `P12345` as a single prefix with no number — so the table's `p` verdict
    // never fired and `P12345:g.6C>G` was emitted. `is_uniprot` reads the same
    // shape off the accession as a whole.
    //
    // Deliberately NOT extended to `ENSG`: `inferred_variant_type` classifies an
    // Ensembl gene as `g`, and that is the codebase's own table. Overriding it
    // here would put a second, disagreeing classifier in exactly the place the
    // comment above argues against one.
    if inferred.is_some_and(|axis| axis != "g")
        || NON_GENOMIC_PREFIXES.contains(&&*prefix)
        || parsed.is_uniprot()
    {
        return Err(non_genomic_refusal(accession));
    }

    // `m.` on the two rCRS accessions; see this function's doc.
    let axis = if parsed.is_mitochondrial() { "m" } else { "g" };
    crate::parse_hgvs(&format!("{accession}:{axis}.=")).map_err(|_| {
        FerroError::InvalidCoordinates {
            msg: format!("'{accession}' is not a usable genomic accession"),
        }
    })
}

/// The refusal for an accession that names a transcript or a protein.
fn non_genomic_refusal(accession: &str) -> FerroError {
    FerroError::InvalidCoordinates {
        msg: format!(
            "'{accession}' names a transcript or protein rather than a genomic sequence; \
             from_sequences derives genomic (g.) and mitochondrial (m.) descriptions only, and a \
             g. description on such a reference is not one the recommendations admit \
             (checklist.md:20). Project the result onto a transcript axis instead."
        ),
    }
}

/// Whether `accession` has UniProt's shape — one upper-case letter followed by
/// five alphanumerics — read off the raw string.
///
/// `Accession::is_uniprot` asks the same question of a *split* accession
/// (`prefix` one character, `number` five). This parser does not always split a
/// UniProt accession that way, so the predicate answered `false` for `P12345`
/// and the protein sailed through the genomic gate. Kept next to the gate rather
/// than pushed into `Accession`, because widening that predicate is a change to
/// a type every axis reads and is not this change's to make.
fn is_uniprot_shaped(accession: &str) -> bool {
    let bytes = accession.as_bytes();
    bytes.len() == 6
        && bytes[0].is_ascii_uppercase()
        && bytes[1..].iter().all(u8::is_ascii_alphanumeric)
}

/// Re-type any member whose payload is the reverse complement of the bases it
/// replaces as an `inv`.
///
/// # Why this pass exists at all
///
/// [`super::merge::anchor_for_piece`] deliberately does **not** type inversions.
/// Its doc says so and gives the reason: on the shipped path the `inv` "is left
/// to `crate::normalize::rules`'s single-span typing, which can still see it in
/// the rendered member" — that is, `canonicalize_from_sequence`'s output is an
/// *intermediate*, handed back to `normalize_core` for exactly this kind of
/// typing. This module does not run `normalize_core`, so without this pass it
/// ships the intermediate.
///
/// Measured before the fix: `AACC` -> `GGTT` came back as
/// `g.11_14delinsGGTT` where `normalize` gives `g.11_14inv`, and both
/// Mutalyzer 3 and VariantValidator rewrite the `delins` spelling to `inv` on
/// `NC_000001.11:g.1000000_1000003delinsCACC`.
///
/// **Whether the un-typed form is non-conformant (rule 1) or merely
/// non-preferred (rule 2) is NOT settled here, and the fix does not depend on
/// it.** `DNA/delins.md:5` defines a delins as a replacement "**and which is
/// not** a substitution or inversion", which reads as a definitional exclusion;
/// but `general.md:56`'s preference list does not rank `delins` at all, and no
/// `class="invalid"` marks this shape. Both tools *rewrote* rather than
/// rejected, which is what one does with a valid-but-non-preferred form. The
/// question belongs in the ruling ledger; emitting `inv` is right either way.
///
/// # Scope, stated so the gap is not mistaken for completeness
///
/// This closes the **inversion** half of the typing gap. A sweep of 6 000
/// random shapes against `normalize` also found **repeat notation**
/// (`g.27_28insAAA` -> `g.27A[4]`) still un-typed here, because a tandem tract
/// can extend past the caller's window and so is not always decidable from it.
/// Everything else the sweep surfaced is reference-anchored member
/// re-derivation — rules 2 and 3, which this module never claimed.
fn retype_inversions(members: &mut [HgvsVariant], reference: &[u8], w_lo: i64) {
    for member in members.iter_mut() {
        // `Mt` alongside `Genome`: a mitochondrial member carries the same
        // `GenomeInterval`/`NaEdit` shape, and leaving it out would have made
        // `inv` typing silently axis-dependent the moment `m.` started being
        // emitted.
        let loc_edit = match member {
            HgvsVariant::Genome(variant) => &mut variant.loc_edit,
            HgvsVariant::Mt(variant) => &mut variant.loc_edit,
            _ => continue,
        };
        // `Mu::Certain` only: an edit the input marked uncertain must not be
        // sharpened into a definite `inv`.
        let Mu::Certain(NaEdit::Delins { sequence, .. }) = &loc_edit.edit else {
            continue;
        };
        let InsertedSequence::Literal(payload) = sequence else {
            continue;
        };
        let payload: Vec<u8> = payload.to_string().into_bytes();
        // `inversion.md:14` — an inversion covers more than one nucleotide; a
        // one-base complement is a substitution, which `build_naedit` already
        // renders from the single reference base.
        if payload.len() < 2 {
            continue;
        }
        // Only a plain, certain two-endpoint span is re-typed. An uncertain or
        // ranged boundary states less than an `inv` would claim, so it is left
        // exactly as it was rather than sharpened.
        let (
            UncertainBoundary::Single(Mu::Certain(start)),
            UncertainBoundary::Single(Mu::Certain(end)),
        ) = (&loc_edit.location.start, &loc_edit.location.end)
        else {
            continue;
        };
        if start.special.is_some() || end.special.is_some() {
            continue;
        }
        // Window offsets. Both endpoints are 1-based inclusive on the same axis
        // as `w_lo`, so the span is `[start - w_lo, end - w_lo]`.
        let (Ok(lo), Ok(hi)) = (
            usize::try_from(start.base as i64 - w_lo),
            usize::try_from(end.base as i64 - w_lo),
        ) else {
            continue;
        };
        if hi < lo || hi >= reference.len() {
            continue;
        }
        // Delegated to `rules::canonicalize_delins`, the same classifier the
        // shipped path uses, rather than a local reverse-complement test. A
        // local test was the first attempt and it missed 85 of 6 000 sampled
        // shapes, because the classifier trims shared affixes *before* testing
        // and can therefore shorten the span — `g.12_16delinsCTTTT` is an
        // inversion of a sub-range, which a whole-span comparison cannot see.
        let outcome = crate::normalize::rules::canonicalize_delins(reference, lo, hi + 1, &payload);
        let crate::normalize::rules::DelinsCanonical::Inversion {
            start: inv_lo,
            end: inv_hi,
        } = outcome
        else {
            continue;
        };
        // `canonicalize_delins` reports 0-based half-open offsets into the
        // window; the description carries 1-based inclusive positions on the
        // caller's axis.
        let (Ok(new_start), Ok(new_end)) = (
            u64::try_from(inv_lo as i64 + w_lo),
            u64::try_from(inv_hi as i64 + w_lo - 1),
        ) else {
            continue;
        };
        loc_edit.location.start =
            UncertainBoundary::Single(Mu::Certain(crate::hgvs::location::GenomePos {
                base: new_start,
                special: None,
                offset: None,
            }));
        loc_edit.location.end =
            UncertainBoundary::Single(Mu::Certain(crate::hgvs::location::GenomePos {
                base: new_end,
                special: None,
                offset: None,
            }));
        loc_edit.edit = Mu::Certain(NaEdit::Inversion {
            sequence: None,
            length: None,
        });
    }
}

/// Re-apply the derived description to the supplied reference and require it to
/// reproduce `alternate` byte for byte.
///
/// A **runtime** check, not a `debug_assert`. Emitting a description that
/// denotes different bases is the worst failure available here, and a debug
/// assertion is compiled out of exactly the builds that process real data — the
/// same reasoning `canonicalize_from_sequence` gives for its own round trip.
///
/// **It is not one of the four seam oracles, and it is not a substitute for
/// them.** Those run from `Normalizer::assert_seam_oracles`, at the single exit
/// of `normalize_core_checked`, which only
/// [`crate::Normalizer::from_sequences`] with `normalize = true` reaches — the
/// free functions and [`crate::SequencePair::derive`] hold no provider and so
/// cannot reach it at all. In particular this shares
/// [`super::merge::apply_edits_to_window`] with the derivation it is checking,
/// where `FERRO_ASSERT_SEQUENCE`'s
/// [`crate::spdi::compare_denoted_sequences`] deliberately routes through
/// `hgvs_to_spdi` so that it does not.
///
/// The gap is stated rather than papered over because it cannot be closed here:
/// the three reference-free oracles have nothing to compare against (there is no
/// *input* description — the input is bases), and the denoted-sequence oracle
/// needs a reference this path does not read. The corpus and multi-member axes
/// are where the independent comparison is made instead, and both make it
/// through `hgvs_to_spdi`.
fn verify_round_trip(
    variant: &HgvsVariant,
    w_lo: i64,
    reference: &[u8],
    alternate: &[u8],
) -> Result<(), FerroError> {
    // Borrowed, not cloned. This runs on **every** derivation — it is a runtime
    // check by design, not a `debug_assert` — and `denoted_bases` takes a
    // slice, so cloning every member here bought nothing but a `Vec` and one
    // deep copy per member on the hot path.
    let members: &[HgvsVariant] = match variant {
        HgvsVariant::Allele(allele) => &allele.variants,
        // The identity description denotes the reference; there is nothing to
        // apply, and `collect_canonical_edits` has no edit to collect.
        _ if reference == alternate => return Ok(()),
        other => std::slice::from_ref(other),
    };
    let rebuilt =
        denoted_bases(members, reference, w_lo).ok_or_else(|| FerroError::ConversionError {
            msg: format!("derived description {variant} could not be re-applied to its window"),
        })?;
    if rebuilt != alternate {
        return Err(FerroError::ConversionError {
            msg: format!(
                "derived description {variant} denotes different bases: expected {}, got {}",
                String::from_utf8_lossy(alternate),
                String::from_utf8_lossy(&rebuilt)
            ),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The block from `fulcrumgenomics/ferro-hgvs#1420`'s comment `5253702249`:
    /// one read, four aligner spellings, one variant. `from_sequences` sees only
    /// the bases, so there is no spelling left for an aligner to perturb.
    #[test]
    fn the_1420_read_block_derives_its_alignment_form() {
        let derived = from_sequences(
            "NC_TEST.1",
            130,
            "GGCGAC",
            "TACCGAGCTT",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert_eq!(
            derived.to_string(),
            "NC_TEST.1:g.[130_131delinsTAC;134_135insG;135_136insTT]"
        );
    }

    /// Two 1 bp deletions separated by two unchanged bases. `partition_block`
    /// cannot cut this — it searches single-gap alignments only — which is why
    /// the shipped normalizer leaves the equivalent spanning `delins` alone.
    #[test]
    fn two_deletions_two_bases_apart_are_two_members() {
        let derived = from_sequences(
            "NC_TEST.1",
            5,
            "AGCG",
            "GC",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert_eq!(derived.to_string(), "NC_TEST.1:g.[5del;8del]");
    }

    /// An unchanged pair is not an error; it denotes the reference.
    #[test]
    fn an_unchanged_pair_is_the_identity_description() {
        let derived = from_sequences(
            "NC_TEST.1",
            5,
            "ACGT",
            "ACGT",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert_eq!(derived.to_string(), "NC_TEST.1:g.=");
    }

    #[test]
    fn a_zero_position_is_refused() {
        let err = from_sequences("NC_TEST.1", 0, "A", "G", &FromSequencesOptions::default())
            .expect_err("refuses");
        assert!(err.to_string().contains("1-based"), "{err}");
    }

    #[test]
    fn an_empty_reference_is_refused() {
        let err = from_sequences("NC_TEST.1", 5, "", "G", &FromSequencesOptions::default())
            .expect_err("refuses");
        assert!(err.to_string().contains("no interval"), "{err}");
    }

    /// **The rCRS gets `m.`, not `g.`** — `Accession::is_mitochondrial`'s own
    /// doc says HGVS requires the `m.` coordinate system for these accessions,
    /// so emitting `g.` here was a rule-1 defect on the most-used mitochondrial
    /// reference. `NC_001807` is the older rCRS draft and is in the same table.
    #[test]
    fn a_mitochondrial_accession_derives_on_the_m_axis() {
        for accession in ["NC_012920.1", "NC_001807.4"] {
            let derived = from_sequences(accession, 6, "C", "T", &FromSequencesOptions::default())
                .expect("derives");
            assert_eq!(derived.to_string(), format!("{accession}:m.6C>T"));
        }
        // The neighbouring non-mitochondrial `NC_` is untouched.
        let genomic = from_sequences(
            "NC_000001.11",
            6,
            "C",
            "T",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert_eq!(genomic.to_string(), "NC_000001.11:g.6C>T");
    }

    /// A multi-base mitochondrial inversion still types as `inv`, which pins
    /// that `retype_inversions` sees an `Mt` member and not only a `Genome` one.
    #[test]
    fn a_mitochondrial_inversion_is_typed() {
        let derived = from_sequences(
            "NC_012920.1",
            11,
            "AACC",
            "GGTT",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert_eq!(derived.to_string(), "NC_012920.1:m.11_14inv");
    }

    /// **A soft-masked window derives, and its round trip passes.**
    ///
    /// `validate` always admitted lower case (`Base::from_char` folds), but
    /// `apply_edits_to_window` copies the reference verbatim while splicing an
    /// upper-case payload, so `verify_round_trip`'s byte comparison rejected a
    /// derivation that was correct: `acgt` -> `atgt` refused with "expected
    /// atgt, got aTgt". `del` and `dup` passed, since neither carries a payload
    /// — which is why this pins a substitution, a delins and an inversion.
    #[test]
    fn a_soft_masked_window_derives() {
        for (reference, alternate, expected) in [
            ("acgt", "atgt", "NC_TEST.1:g.6C>T"),
            ("acgt", "act", "NC_TEST.1:g.7del"),
            ("acgt", "acggt", "NC_TEST.1:g.7dup"),
            ("acgtacgt", "acgttgca", "NC_TEST.1:g.9_12delinsTGCA"),
            ("ttaacct", "ttggttt", "NC_TEST.1:g.7_10inv"),
        ] {
            let derived = from_sequences(
                "NC_TEST.1",
                5,
                reference,
                alternate,
                &FromSequencesOptions::default(),
            )
            .unwrap_or_else(|e| panic!("{reference}/{alternate} must derive: {e}"));
            assert_eq!(derived.to_string(), expected);
        }
    }

    /// Case is folded, so a masked window and its upper-case twin are the same
    /// input — including the mixed-case pair `to_sequences` used to manufacture
    /// by upper-casing only its 5' flank.
    #[test]
    fn case_does_not_change_the_answer() {
        let options = FromSequencesOptions::default();
        let upper = from_sequences("NC_TEST.1", 5, "ACGTACGT", "ACGTTGCA", &options)
            .expect("derives")
            .to_string();
        for (reference, alternate) in [
            ("acgtacgt", "acgttgca"),
            ("ACGTacgt", "ACGTtgca"),
            ("acgtACGT", "acgtTGCA"),
        ] {
            let derived = from_sequences("NC_TEST.1", 5, reference, alternate, &options)
                .expect("derives")
                .to_string();
            assert_eq!(derived, upper, "{reference}/{alternate}");
        }
    }

    /// `U` is RNA and this surface's axis is DNA. It used to be admitted
    /// asymmetrically — `ACGU` -> `ACGT` emitted `g.8U>T` while `ACGT` ->
    /// `AUGT` refused with the round trip's internal message.
    #[test]
    fn uracil_is_refused_on_both_sides() {
        for (reference, alternate) in [("ACGT", "AUGT"), ("ACGU", "ACGT"), ("acgu", "acgt")] {
            let err = from_sequences(
                "NC_TEST.1",
                5,
                reference,
                alternate,
                &FromSequencesOptions::default(),
            )
            .expect_err("refuses");
            assert!(err.to_string().contains("which is RNA"), "{err}");
        }
    }

    /// Every non-genomic class reaches the informative refusal, including the
    /// two that used to fall out of the parse arm with "not a usable genomic
    /// accession" (`NR_`/`XR_`, which the parser refuses a `g.` axis on) and the
    /// UniProt shape, which passed the gate outright and emitted `P12345:g.6C>G`.
    #[test]
    fn every_non_genomic_class_gets_the_named_refusal() {
        for accession in [
            "NM_000518.4",
            "NR_000001.1",
            "XR_000001.1",
            "XM_000001.1",
            "NP_000509.1",
            "ENST00000012345.1",
            "ENSMUST00000012345.1",
            "LRG_1t1",
            "LRG_1p1",
            "P12345",
        ] {
            let err = from_sequences(
                accession,
                5,
                "ACGT",
                "AGGT",
                &FromSequencesOptions::default(),
            )
            .expect_err("refuses");
            assert!(
                err.to_string().contains("names a transcript or protein"),
                "{accession}: {err}"
            );
        }
    }

    /// **An accession the classifier cannot place is accepted, not refused.**
    ///
    /// `inferred_variant_type` returns `None` for a SAM refname, a custom contig
    /// and an assembly reference, and `None` is not a verdict of "non-genomic" —
    /// ferro cannot call `chr1:g.6C>G` wrong. Pinned because classifying off the
    /// **bare** accession string broke every one of these at once: the parser
    /// finds the end of an accession by locating the HGVS separator, so a bare
    /// `chr1` is not fully consumed and fell into the "not a usable genomic
    /// accession" arm. The gate now probes with the `:g.=` suffix attached.
    #[test]
    fn an_unclassifiable_accession_is_still_derived() {
        for accession in [
            "chr1",
            "chrM",
            "1",
            "scaffold_123",
            "my-contig.v2",
            "GRCh38(chr1)",
            "hg38(chr1)",
            "NZ_CP012345.1",
            "AC_000001.1",
            "NW_000001.1",
            "LRG_1",
        ] {
            let derived = from_sequences(
                accession,
                5,
                "ACGT",
                "AGGT",
                &FromSequencesOptions::default(),
            )
            .unwrap_or_else(|e| panic!("{accession} must derive: {e}"));
            assert_eq!(derived.to_string(), format!("{accession}:g.6C>G"));
        }
    }

    #[test]
    fn an_alignment_only_symbol_is_refused_by_clause() {
        let err = from_sequences(
            "NC_TEST.1",
            5,
            "ACGT",
            "ACXT",
            &FromSequencesOptions::default(),
        )
        .expect_err("refuses");
        assert!(err.to_string().contains("standards.md:39"), "{err}");
    }

    /// The budget refuses rather than degrading, and the message names the knob
    /// and its cost — a refusal that does not say how to proceed is a dead end.
    #[test]
    fn a_window_over_the_grid_budget_is_refused_with_the_knob_named() {
        let options = FromSequencesOptions {
            max_grid_cells: 16,
            ..Default::default()
        };
        let err =
            from_sequences("NC_TEST.1", 5, "ACGTACGT", "TGCATGCA", &options).expect_err("refuses");
        assert!(err.to_string().contains("max_grid_cells"), "{err}");
    }

    /// A member on the window edge is where a window-local answer and a
    /// reference-anchored one can differ, so it is reported rather than absorbed.
    #[test]
    fn a_member_on_the_window_edge_is_reported() {
        let interior = from_sequences_detailed(
            "NC_TEST.1",
            5,
            "TAAAAG",
            "TAAAG",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert!(
            !interior.placement_bounded_by_window(),
            "an interior deletion must not flag: {}",
            interior.variant
        );

        let edge = from_sequences_detailed(
            "NC_TEST.1",
            5,
            "AAAA",
            "AAA",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert!(
            edge.placement_bounded_by_window(),
            "a deletion flush with the window end must flag: {}",
            edge.variant
        );
    }

    /// The combined flag says *that* a member is on an edge; the per-side flags
    /// say *which*. A deletion in a homopolymer that fills the window rolls, by
    /// default, to the 3' edge — so `bounded_at_end` alone must fire, and a
    /// caller widening one side at a time must not be told to widen 5'.
    #[test]
    fn a_deletion_at_the_three_prime_edge_flags_only_that_side() {
        let d = from_sequences_detailed(
            "NC_TEST.1",
            5,
            "AAAA",
            "AAA",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert!(d.bounded_at_end, "3' edge must flag: {}", d.variant);
        assert!(!d.bounded_at_start, "5' edge must not flag: {}", d.variant);
        assert!(d.placement_bounded_by_window(), "OR must hold");
    }

    /// The 5' mirror: under 5'-shuffle the same deletion rolls to the window's
    /// 5' edge instead, so `bounded_at_start` alone fires. This is the side a
    /// window pinned to base 1 has permanently, and the reason
    /// `sequence_normalize` reads the two flags apart.
    #[test]
    fn a_deletion_at_the_five_prime_edge_flags_only_that_side() {
        let d = from_sequences_detailed(
            "NC_TEST.1",
            5,
            "AAAA",
            "AAA",
            &FromSequencesOptions::default().with_direction(ShuffleDirection::FivePrime),
        )
        .expect("derives");
        assert!(d.bounded_at_start, "5' edge must flag: {}", d.variant);
        assert!(!d.bounded_at_end, "3' edge must not flag: {}", d.variant);
        assert!(d.placement_bounded_by_window(), "OR must hold");
    }

    /// **Contig-start escape.** An insertion in an ambiguous run reaching the
    /// accession's first base has no 5'-most HGVS anchor — the 5'-shuffle rolls
    /// it to interbase 0, which is spelled `0_1ins…` and names a position that
    /// does not exist. Because `w_lo == 1` here, that interbase *is* the contig
    /// start, so widening cannot rescue it. It is instead re-presented at the
    /// leftmost nameable interbase (offset 1), the terminal analogue of the
    /// 3'-most rule: inserting `AT` before base 1 equals inserting `TA` before
    /// base 2, and only the latter can be written. This is the shape the real
    /// run's `1A>C`/`1A>T`-anchored multi-substitution alleles hit.
    #[test]
    fn a_contig_start_insertion_escapes_to_the_leftmost_nameable_interbase() {
        let derived = from_sequences(
            "NC_TEST.1",
            1,
            "ATATATG",
            "ATATATATG",
            &FromSequencesOptions::default().with_direction(ShuffleDirection::FivePrime),
        )
        .expect("the ambiguous-run insertion escapes to interbase 1 rather than refusing");
        assert_eq!(derived.to_string(), "NC_TEST.1:g.1_2insTA");
    }

    /// A pure insertion before base 1 with **no piece to fold into** is a
    /// genuine insertion at a non-existent anchor and still refuses. Here the
    /// whole block is a single inserted `G` (`CATATG` -> `GCATATG` share the
    /// suffix `CATATG`), so there is no following delins to absorb it, and the
    /// leading base is not itself rewritten — the span-form escape does not
    /// apply. The base-1-rewriting span form is covered in the
    /// `sequence_normalize` integration suite, which drives this same path.
    #[test]
    fn a_lone_contig_start_insertion_still_refuses() {
        let err = from_sequences(
            "NC_TEST.1",
            1,
            "CATATG",
            "GCATATG",
            &FromSequencesOptions::default().with_direction(ShuffleDirection::FivePrime),
        )
        .expect_err("an insertion genuinely before base 1 has no HGVS anchor");
        assert!(
            err.to_string().contains("5' of the window's first base"),
            "{err}"
        );
    }

    /// A change spanning the whole window rests on both edges at once, so both
    /// per-side flags fire. Neither side can be dismissed as settled, which is
    /// what makes a two-sided widen the right response.
    #[test]
    fn a_change_filling_the_window_flags_both_sides() {
        let d = from_sequences_detailed(
            "NC_TEST.1",
            5,
            "ACGT",
            "TGCA",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert!(d.bounded_at_start, "5' edge must flag: {}", d.variant);
        assert!(d.bounded_at_end, "3' edge must flag: {}", d.variant);
    }

    /// An interior change touches neither edge, so both per-side flags — and
    /// therefore the combined flag — are clear.
    #[test]
    fn an_interior_change_flags_neither_side() {
        let d = from_sequences_detailed(
            "NC_TEST.1",
            5,
            "TAAAAG",
            "TAAAG",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert!(!d.bounded_at_start, "5' edge must not flag: {}", d.variant);
        assert!(!d.bounded_at_end, "3' edge must not flag: {}", d.variant);
        assert!(!d.placement_bounded_by_window(), "OR must be clear");
    }

    /// The same four arguments give the same string. Asserted rather than
    /// assumed, because "deterministic" is the one property this whole surface
    /// is sold on.
    #[test]
    fn the_same_arguments_give_the_same_description() {
        let once = from_sequences(
            "NC_TEST.1",
            130,
            "GGCGAC",
            "TACCGAGCTT",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        let twice = from_sequences(
            "NC_TEST.1",
            130,
            "GGCGAC",
            "TACCGAGCTT",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert_eq!(once.to_string(), twice.to_string());
    }

    /// Padding the window with matching flank must not move the answer, until
    /// the edge stops binding. This is the property that makes a caller's choice
    /// of window size safe when the change is interior.
    #[test]
    fn matching_flank_does_not_move_an_interior_answer() {
        let tight = from_sequences(
            "NC_TEST.1",
            10,
            "TAGCGA",
            "TGCA",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        let padded = from_sequences(
            "NC_TEST.1",
            8,
            "CCTAGCGAGG",
            "CCTGCAGG",
            &FromSequencesOptions::default(),
        )
        .expect("derives");
        assert_eq!(tight.to_string(), padded.to_string());
    }
}
