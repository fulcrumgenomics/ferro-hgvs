//! `standards.md:39`'s alignment-only symbols, and where a description can still state one.
//!
//! # The rule
//!
//! `background/standards.md`'s DNA symbol table daggers two of its rows — `X†`
//! ("masked nucleotide", `:36`) and `-†` ("gap of indeterminate length", `:37`)
//! — and footnotes the dagger at `:39`: *"used in alignment only."* The dagger
//! is printed **inside** the symbol cell, so `:39` annotates those two rows
//! rather than competing with them, and `general.md:48` independently admits
//! only "nucleotides in CAPITALS using [IUPAC-IUBMB assigned nucleotide
//! symbols]", of which neither is one.
//!
//! The decided `rulings[alignment-only-symbol-in-a-description]` reads those
//! together: **neither symbol may appear in a description, and ferro must
//! refuse both.** The decided `rulings[absolute-prohibition-enforcement-stage]`
//! says *where*: strict fails at PARSE, lenient does not validate input
//! conformance and fails only when it cannot NORMALIZE, silent is lenient
//! without messages.
//!
//! # Where the symbol actually survives, which is narrower than it reads
//!
//! Measured on `78c43230`, `-` never reaches the AST at all: it matches no arm
//! of `parse_inserted_sequence`, so it is left unconsumed and the
//! variant-level trailing-character check refuses it — in every mode, as a
//! *grammar* accident rather than as a conformance rule. `X` reaches the AST
//! through exactly one carrier, [`InsertedSequence::Named`], and through **two**
//! arms of that function:
//!
//! | input | arm | result |
//! |---|---|---|
//! | `g.10delinsX` | `c if c.is_ascii_uppercase()` | `Named("X")` |
//! | `g.10delinsACGTX` | `c if is_iupac_base(c)` | `Named("ACGTX")` |
//! | `g.10delinsXACGT` | `c if c.is_ascii_uppercase()` | `Named("XACGT")` |
//! | `g.10delinsACXGT` | `c if is_iupac_base(c)` | `Named("ACXGT")` |
//!
//! The second arm walks the whole run and sets `has_non_iupac` on any uppercase
//! byte that is not an IUPAC base, so an otherwise-literal run carrying one
//! stray `X` is reclassified **in its entirety** as a mobile-element name. That
//! is why the predicate below tests for the *presence* of the symbol rather
//! than for a name that *is* the symbol: a check keyed on a lone letter closes
//! `delinsX` and leaves `delinsACGTX` open, which is how every "24 rows" figure
//! quoted against `standards.md:39` came to be a lone-`X` count.
//!
//! # What it deliberately does not reach
//!
//! **Genuine mobile-element names.** `AluYb8`, `LINE1`, `L1` and `Alu` are the
//! shape the named-element arm exists for, and none contains `X` or `-`, so
//! none is touched. Whether a named element is conformant at all is a separate
//! question — `DNA/complex.md:169` answers "No, not really, it is not exact"
//! about `insL1.603bp` — and **no ruling settles it**, so it is left alone here
//! rather than folded into this one.
//!
//! **Lowercase `x`.** `g.10delinsACGTx` also reaches `Named`, via the
//! `is_ascii_lowercase` branch of the same arm. It is not matched, because the
//! symbol table is uppercase and a lowercase letter inside a named element is
//! indistinguishable from an ordinary name character (`Alu` has three). Ruling
//! on it would be adjudicating `general.md:48`'s CAPITALS requirement against
//! the named-element arm, which nothing has decided.
//!
//! **That argument is DNA-axis reasoning, and it does not transfer to `r.` —
//! measured, and left open deliberately.** `general.md:48` is the **DNA-level**
//! bullet; its RNA-level sibling `general.md:50` says "nucleotides in **lower
//! case** using [IUPAC-IUBMB assigned nucleotide symbols]". So on the RNA axis
//! the ordinary sequence alphabet *is* lowercase, and the "indistinguishable
//! from a name character" premise inverts. Measured directly:
//!
//! | input | parses to | this predicate |
//! |---|---|---|
//! | `r.10delinsacgu` | `Literal([A, C, G, U])` | n/a — not a name |
//! | `r.10delinsacgux` | `Named("acgux")` | **misses** |
//! | `r.10delinsX` | `Named("X")` | matches |
//!
//! So on `r.` the predicate catches the spelling the axis does not use and
//! misses the one it does, and `r.10delinsacgux` is accepted and re-emitted in
//! all three modes exactly as `g.10delinsACGTX` was before this change — the
//! same wholesale reclassification, one axis over.
//!
//! It is **not** fixed here, and that is a scope decision rather than an
//! oversight — but note which ground it would have to rest on, because the two
//! are not the same and only one of them is ruled on. `standards.md:39` daggers
//! rows of the **DNA** symbol table, so it does not reach a lowercase RNA `x` at
//! all; a refusal there would rest on `general.md:50`'s *alphabet* requirement
//! alone, exactly as `general.md:48` is the alphabet half of the DNA case. That
//! is a different clause from the one the decided
//! `rulings[alignment-only-symbol-in-a-description]` names, and **no ruling
//! covers it** — neither that record, which is written against the DNA table,
//! nor anything under `RNA/`. Widening the predicate to lowercase without one
//! would also reach a lowercase letter in a genuine element name, which is the
//! false positive the uppercase-only test exists to avoid. Pinned as an open
//! question by `a_lowercase_rna_masked_nucleotide_is_still_accepted` below, so
//! it cannot rot into a silent gap.
//!
//! The converse over-reach is deliberate and harmless: an **uppercase** `X` on
//! the `r.` axis (`r.10delinsX`) is matched here on a DNA-table clause. It is
//! non-conformant on `general.md:50`'s ground too — the RNA alphabet is
//! lowercase — so both readings refuse it, and the refusal does not depend on
//! which clause is doing the work.
//!
//! **Accessions.** The predicate is AST-keyed, never a scan of the rendered
//! description, because `XM_`/`XR_`/`XP_` RefSeq accessions carry an `X` that
//! has nothing to do with the symbol table.

use crate::hgvs::edit::{InsertedSequence, NaEdit};
use crate::hgvs::variant::HgvsVariant;

/// `standards.md:36` — masked nucleotide, daggered "used in alignment only".
const MASKED_NUCLEOTIDE: char = 'X';

/// `standards.md:37` — gap of indeterminate length, daggered likewise.
const INDETERMINATE_GAP: char = '-';

/// One alignment-only symbol found in a description, with the text that states it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentOnlySymbol {
    /// The daggered symbol itself — `X` or `-`.
    pub symbol: char,
    /// The stated text carrying it, e.g. `ACGTX` for `g.10delinsACGTX`.
    pub stated: String,
}

impl AlignmentOnlySymbol {
    /// The symbol's meaning in `standards.md`'s table, for a diagnostic.
    #[must_use]
    pub fn meaning(&self) -> &'static str {
        if self.symbol == INDETERMINATE_GAP {
            "gap of indeterminate length"
        } else {
            "masked nucleotide"
        }
    }

    /// The clause this violates, quoted where it is cited.
    #[must_use]
    pub fn clause(&self) -> &'static str {
        "background/standards.md:39"
    }
}

/// The first alignment-only symbol stated by `text`, if any.
fn symbol_in(text: &str) -> Option<char> {
    text.chars()
        .find(|c| *c == MASKED_NUCLEOTIDE || *c == INDETERMINATE_GAP)
}

/// The first alignment-only symbol an inserted sequence states, if any.
///
/// **Only [`InsertedSequence::Named`] can carry one**, and the exhaustive match
/// is deliberate so a new variant has to be classified rather than defaulted:
///
/// - `Literal`, `Repeat`, `SequenceRepeat` are built from
///   [`Base`](crate::hgvs::edit::Base), whose discriminants are the fifteen
///   IUPAC symbols plus `U`. Neither `X` nor `-` is representable.
/// - `Count`, `Range`, `PositionRange`, `PositionRangeInv`,
///   `SpecialPositionRange` are numeric.
/// - `Reference` and `Complex`'s `ExternalRef`/`CdsPositionRange` are
///   **positions and accessions, not stated sequences**, and matching them
///   would be a false positive twice over: `ins[244-8_249]` states an intronic
///   offset whose `-` is position syntax, and `XM_`/`XR_`/`XP_` accessions
///   carry an `X` that has nothing to do with the symbol table. The first of
///   those is not hypothetical — an earlier revision of this file rendered each
///   `InsertedPart` and scanned the text, and
///   `ins_intronic_offset_range_returns_unsupported_error` caught it.
/// - `Complex`'s remaining parts (`Literal`, `Repeat`) are `Base`-built, as above.
fn symbol_in_inserted(sequence: &InsertedSequence) -> Option<AlignmentOnlySymbol> {
    match sequence {
        InsertedSequence::Named(name) => symbol_in(name).map(|symbol| AlignmentOnlySymbol {
            symbol,
            stated: name.clone(),
        }),
        InsertedSequence::Empty
        | InsertedSequence::Literal(_)
        | InsertedSequence::Count(_)
        | InsertedSequence::Range(_, _)
        | InsertedSequence::Repeat { .. }
        | InsertedSequence::SequenceRepeat { .. }
        | InsertedSequence::Complex(_)
        | InsertedSequence::Reference(_)
        | InsertedSequence::PositionRange { .. }
        | InsertedSequence::PositionRangeInv { .. }
        | InsertedSequence::SpecialPositionRange { .. }
        | InsertedSequence::UncertainRangeInv { .. }
        | InsertedSequence::Uncertain => None,
    }
}

/// The first alignment-only symbol a nucleic-acid edit states, if any.
///
/// Exhaustive for the same reason [`symbol_in_inserted`] is: the four arms below
/// are the *only* [`NaEdit`] variants carrying an [`InsertedSequence`], and a new
/// one has to be classified here rather than silently defaulted to "states no
/// symbol". Every remaining variant either states no sequence at all (numeric,
/// positional or marker-only) or states one built from
/// [`Sequence`](crate::hgvs::edit::Sequence)/[`Base`](crate::hgvs::edit::Base),
/// whose discriminants are the fifteen IUPAC symbols plus `U` — neither `X` nor
/// `-` is representable. `Conversion`'s `source` is an accession-and-range
/// string, so matching it would be the `XM_`/`XR_`/`XP_` false positive the
/// module docs describe; lenient mode rewrites `con` into a `Delins` carrying an
/// `ExternalRef`, which is classified by `symbol_in_inserted` on that path.
fn symbol_in_edit(edit: &NaEdit) -> Option<AlignmentOnlySymbol> {
    match edit {
        NaEdit::Insertion { sequence }
        | NaEdit::BreakpointInsertion { sequence }
        | NaEdit::Delins { sequence, .. }
        | NaEdit::DupIns { sequence } => symbol_in_inserted(sequence),
        NaEdit::Substitution { .. }
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

/// The first alignment-only symbol `variant` states, if any.
///
/// Walks every nucleic-acid member, including allele members, ring segments and
/// the inner description of a `sup` marker, so a symbol cannot hide behind a
/// composite spelling — a cis allele (`g.[10delinsACGTX;20del]`), a trans pair
/// (`g.[10delinsACGTX];[20del]`) and an uncertain group (`g.(10delinsACGTX)`)
/// are all reached, measured.
///
/// The match is exhaustive rather than wildcarded, on the same reasoning as
/// [`symbol_in_edit`]: a new axis must be classified here rather than defaulted
/// to "carries no symbol". The four that yield `None` do so for stated reasons:
///
/// - `Protein` is out of scope. `standards.md:39` annotates the **DNA** symbol
///   table, and `X` on a protein axis is the parser-native `Xaa` — the
///   deprecated stop-codon spelling is `W3008`/`W3010`'s business.
/// - `RnaFusion` carries two `RnaFusionBreakpoint`s, which hold an accession,
///   an optional gene symbol and an interval — **no edit at all**, so there is
///   no stated sequence for a symbol to sit in.
/// - `NullAllele` (`0`) and `UnknownAllele` (`?`) are whole-allele markers with
///   neither position nor edit.
#[must_use]
pub fn alignment_only_symbol(variant: &HgvsVariant) -> Option<AlignmentOnlySymbol> {
    match variant {
        HgvsVariant::Genome(v) => v.loc_edit.edit.inner().and_then(symbol_in_edit),
        HgvsVariant::Cds(v) => v.loc_edit.edit.inner().and_then(symbol_in_edit),
        HgvsVariant::Tx(v) => v.loc_edit.edit.inner().and_then(symbol_in_edit),
        HgvsVariant::Rna(v) => v.loc_edit.edit.inner().and_then(symbol_in_edit),
        HgvsVariant::Mt(v) => v.loc_edit.edit.inner().and_then(symbol_in_edit),
        HgvsVariant::Circular(v) => v.loc_edit.edit.inner().and_then(symbol_in_edit),
        HgvsVariant::Allele(allele) => allele.variants.iter().find_map(alignment_only_symbol),
        HgvsVariant::GenomeRing(ring) => ring
            .segments
            .iter()
            .find_map(|segment| segment.edit.inner().and_then(symbol_in_edit)),
        HgvsVariant::Supernumerary(inner) => alignment_only_symbol(inner),
        HgvsVariant::Protein(_)
        | HgvsVariant::RnaFusion(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;

    fn symbol_of(input: &str) -> Option<AlignmentOnlySymbol> {
        alignment_only_symbol(&parse_hgvs(input).expect("input parses"))
    }

    #[test]
    fn a_lone_masked_nucleotide_is_found() {
        let found = symbol_of("NC_TEST.1:g.10delinsX").expect("X is stated");
        assert_eq!(found.symbol, 'X');
        assert_eq!(found.stated, "X");
        assert_eq!(found.meaning(), "masked nucleotide");
    }

    /// The half no committed count covered: one stray symbol reclassifies the
    /// whole run, so the check must test for presence, not for a lone letter.
    #[test]
    fn an_embedded_masked_nucleotide_is_found_at_every_offset() {
        for (input, stated) in [
            ("NC_TEST.1:g.10delinsACGTX", "ACGTX"),
            ("NC_TEST.1:g.10delinsXACGT", "XACGT"),
            ("NC_TEST.1:g.10delinsACXGT", "ACXGT"),
            ("NC_TEST.1:g.10_11insACGTX", "ACGTX"),
        ] {
            let found = symbol_of(input).unwrap_or_else(|| panic!("{input}: X is stated"));
            assert_eq!(found.symbol, 'X', "{input}");
            assert_eq!(found.stated, stated, "{input}");
        }
    }

    /// The constraint on any narrowing of the named-element arm: these are what
    /// it exists for, they contain no daggered symbol, and they must stay clear.
    #[test]
    fn genuine_mobile_element_names_state_no_alignment_symbol() {
        for input in [
            "NC_TEST.1:g.10delinsAluYb8",
            "NC_TEST.1:g.10delinsLINE1",
            "NC_TEST.1:g.10delinsL1",
            "NC_TEST.1:g.10delinsAlu",
            "NC_TEST.1:g.10_11insAluYb8",
        ] {
            assert_eq!(
                symbol_of(input),
                None,
                "{input} is a legitimate named element"
            );
        }
    }

    #[test]
    fn a_literal_run_states_no_alignment_symbol() {
        assert_eq!(symbol_of("NC_TEST.1:g.10delinsACGT"), None);
        assert_eq!(symbol_of("NC_TEST.1:g.10delinsN"), None);
    }

    /// An `XM_` accession carries an `X` that is not the symbol table's. The
    /// predicate is AST-keyed precisely so this cannot be a false positive.
    #[test]
    fn an_x_prefixed_accession_is_not_a_masked_nucleotide() {
        assert_eq!(symbol_of("XM_005260378.1:c.10delinsACGT"), None);
    }

    /// A member inside an allele is reached, so the symbol cannot hide behind a
    /// composite spelling. All three composite shapes are covered, because "the
    /// cis case works" does not entail the other two: a trans pair nests one
    /// `Allele` inside another, and an uncertain group sets `uncertain` on it.
    #[test]
    fn a_member_of_any_composite_spelling_is_reached() {
        for input in [
            "NC_TEST.1:g.[10delinsACGT;20delinsACGTX]",
            "NC_TEST.1:g.[10delinsACGTX];[20del]",
            "NC_TEST.1:g.(10delinsACGTX)",
        ] {
            let found = symbol_of(input).unwrap_or_else(|| panic!("{input}: X is stated"));
            assert_eq!(found.symbol, 'X', "{input}");
            assert_eq!(found.stated, "ACGTX", "{input}");
        }
    }

    /// **An open question, pinned as one.** On the `r.` axis the ordinary
    /// sequence alphabet is lowercase, so `r.10delinsacgux` is the same
    /// wholesale reclassification `g.10delinsACGTX` was — and this predicate
    /// misses it, while matching the uppercase spelling that axis does not use.
    ///
    /// Asserted here as today's behaviour, **not** as correct behaviour: no
    /// ruling extends `standards.md:39`'s DNA symbol table to a lowercase RNA
    /// `x`, and the clause that *would* reach it — `general.md:50`'s lowercase
    /// IUPAC-IUBMB alphabet, the RNA sibling of the `general.md:48` this file
    /// argues from — is a different ground that nothing has ruled on. Inventing
    /// the ruling in a fix PR would be exactly the unadjudicated widening the
    /// module docs decline. If a ruling later closes it, this test fails and
    /// names the reason, which is the point of pinning it.
    #[test]
    fn a_lowercase_rna_masked_nucleotide_is_still_accepted() {
        assert_eq!(
            symbol_of("NM_TEST.1:r.10delinsacgux"),
            None,
            "unadjudicated: `standards.md:39` daggers the DNA table, and no ruling reaches \
             a lowercase `x` on the RNA axis"
        );
        // The uppercase spelling on the same axis *is* matched, which is what
        // makes the gap a lowercase one rather than an axis one.
        assert_eq!(
            symbol_of("NM_TEST.1:r.10delinsX").map(|s| s.symbol),
            Some('X')
        );
    }

    /// The DNA-axis half of the same scope line, and the constraint on widening
    /// it: a lowercase letter is an ordinary name character here.
    #[test]
    fn a_lowercase_x_on_a_dna_axis_is_still_accepted() {
        assert_eq!(symbol_of("NC_TEST.1:g.10delinsACGTx"), None);
    }
}
