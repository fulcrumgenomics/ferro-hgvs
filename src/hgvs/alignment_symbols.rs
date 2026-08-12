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
//! **That argument is DNA-axis reasoning, and it does not transfer to `r.`** —
//! which is the whole of `#1715` and the reason this predicate is keyed on the
//! axis rather than on one symbol set. `general.md:48` is the **DNA-level**
//! bullet; its RNA-level sibling `general.md:50` says "nucleotides in **lower
//! case** using [IUPAC-IUBMB assigned nucleotide symbols]". So on the RNA axis
//! the ordinary sequence alphabet *is* lowercase, and the "indistinguishable
//! from a name character" premise inverts. Measured on `2f4e3bb9`, before this
//! change:
//!
//! | input | parses to | the DNA-only predicate |
//! |---|---|---|
//! | `r.10delinsacgu` | `Literal([A, C, G, U])` | n/a — not a name |
//! | `r.10delinsacgux` | `Named("acgux")` | **missed** |
//! | `r.10delinsX` | `Named("X")` | matched |
//!
//! So on `r.` it caught the spelling the axis does not use and missed the one it
//! does, and `r.10delinsacgux` was accepted and re-emitted in all three modes
//! with an empty warning vector, exactly as `g.10delinsACGTX` was before #1627 —
//! the same wholesale reclassification, one axis over.
//!
//! # The RNA table is the DNA table minus its two daggered rows
//!
//! Which is what lets the `r.` half rest on the RNA axis's *own* jurisdiction
//! rather than on an extension of a DNA clause. `standards.md`'s RNA table
//! (`:45`–`:61`) lists exactly `a c g u b d h k m n r s v w y` — the DNA table's
//! fifteen **non**-daggered rows, lowercased, with `u` for `t` — and **stops**.
//! `X†` and `-†` have no RNA counterpart at all. So the RNA alphabet's exclusion
//! of `x` is stated by the spec's own tabulation, and `general.md:50` supplies
//! the case: an `r.` description states those symbols in lower case.
//!
//! That matters because a `DNA/` clause cannot scope `r.` (see `CLAUDE.md`,
//! "Cite the clause exactly"), and `standards.md` is a `background/` document
//! carrying **both** tables rather than a DNA-only one — so `:47`–`:61`, the
//! table's rows, is an RNA-jurisdiction citation and not a borrowed DNA one.
//! That is the range the ledger record cites; `:45` is the header row.
//!
//! **This is an ADJUDICATION and it was the operator's**, recorded as
//! `rulings[rna-axis-alignment-only-symbol-reach]` and ruled on 2026-08-12:
//! REFUSE, with `standards.md:47`–`:61` governing. That record states the
//! question, quotes the clauses, and gives the reading it was weighed against.
//! Read it before arguing the point from the spec text.
//!
//! The converse over-reach is deliberate and harmless: an **uppercase** `X` on
//! the `r.` axis (`r.10delinsX`) is matched on a DNA-table clause. It is
//! non-conformant on `general.md:50`'s ground too — the RNA alphabet is
//! lowercase — so both readings refuse it, and the refusal does not depend on
//! which clause is doing the work. Its message is left citing `standards.md:39`
//! and `general.md:48` for exactly that reason: nothing about the uppercase
//! spelling's verdict or diagnostic moves here.
//!
//! # What the `r.` half deliberately does NOT do
//!
//! **It is not an alphabet rule.** The predicate looks for the *symbol* `x`, not
//! for "a character outside the fifteen RNA symbols". A closed-alphabet check
//! would reach `Named("alu")` — measured, `r.10delinsalu` parses to exactly that
//! — and whether a named mobile element is conformant at all is the separate
//! question `DNA/complex.md:169` raises and **no ruling settles**. Widening to
//! the alphabet is a strictly larger adjudication than the one being asked for.
//!
//! **Lowercase `x` on a DNA axis is still accepted**, unchanged: `g.10delinsACGTx`
//! and `n.10delinsACGTx` reach `Named` and are left alone, because there
//! `general.md:48`'s CAPITALS makes a lowercase letter an ordinary name
//! character (`Alu` has three). `n.` is a non-coding **DNA** reference sequence,
//! so it takes the DNA table; only `r.` does not.
//!
//! # The reachable surface on `r.` is NARROWER than on the DNA axes, measured
//!
//! `parse_inserted_sequence` dispatches on the *first* byte, and no arm accepts a
//! leading lowercase non-IUPAC letter. So on `r.`, unlike the uppercase case:
//!
//! | input | reaches the AST as | why |
//! |---|---|---|
//! | `r.10delinsx` | — | grammar refuses it, like `-` |
//! | `r.10delinsxacgu` | — | grammar refuses it, like `-` |
//! | `r.10delinsacgux` | `Named("acgux")` | `is_iupac_base` arm, stray at the end |
//! | `r.10delinsacxgu` | `Named("acxgu")` | `is_iupac_base` arm, stray inside |
//! | `r.10delinsAcgux` | `Named("Acgux")` | uppercase arm, stray inside |
//!
//! **A lone or leading `x` is therefore refused for grammar reasons in every
//! mode, and the two-stage mode schedule does not describe it** — the same
//! standing fact `-` has, and the same misattribution to avoid. Only a
//! *non-leading* `x` is mode-dependent. Pinned in
//! `tests/it/issue_1715_rna_alignment_symbol_reach.rs`.
//!
//! **Accessions.** The predicate is AST-keyed, never a scan of the rendered
//! description, because `XM_`/`XR_`/`XP_` RefSeq accessions carry an `X` that
//! has nothing to do with the symbol table.

use crate::hgvs::edit::{InsertedSequence, NaEdit};
use crate::hgvs::variant::HgvsVariant;

/// `standards.md:36` — masked nucleotide, daggered "used in alignment only".
const MASKED_NUCLEOTIDE: char = 'X';

/// The same row's RNA spelling (`#1715`).
///
/// `standards.md`'s RNA table (`:45`–`:61`) is the DNA table's fifteen
/// **non**-daggered rows lowercased with `u` for `t`, and omits both daggered
/// rows; `general.md:50` puts an `r.` description's nucleotides in lower case.
/// Reached only on the `r.` axis — see [`SymbolTable`].
const MASKED_NUCLEOTIDE_RNA: char = 'x';

/// `standards.md:37` — gap of indeterminate length, daggered likewise. Caseless,
/// so it is one character on both axes.
const INDETERMINATE_GAP: char = '-';

/// Which of `standards.md`'s two symbol tables an axis's descriptions are
/// written in.
///
/// `general.md:46` opens "descriptions on DNA, RNA, and protein level are
/// clearly different" and the two sub-bullets say how: `:48` CAPITALS for DNA,
/// `:50` lower case for RNA. The split is therefore a property of the *axis*,
/// which is why it is resolved once in [`alignment_only_symbol`] and threaded
/// down rather than re-derived per edit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SymbolTable {
    /// `standards.md:19`–`:37`. Every nucleic-acid axis except `r.` — including
    /// `n.`, which is a non-coding **DNA** reference sequence.
    Dna,
    /// `standards.md:45`–`:61`. The `r.` axis alone.
    Rna,
}

impl SymbolTable {
    /// The symbols a description on this axis may not state.
    ///
    /// `r.` carries the DNA pair as well as its own lower-case spelling. That is
    /// the deliberate over-reach the module docs describe: an uppercase `X` on
    /// `r.` is non-conformant under either reading, so nothing about its verdict
    /// depends on which clause is doing the work.
    fn daggered(self) -> &'static [char] {
        match self {
            Self::Dna => &[MASKED_NUCLEOTIDE, INDETERMINATE_GAP],
            Self::Rna => &[MASKED_NUCLEOTIDE, MASKED_NUCLEOTIDE_RNA, INDETERMINATE_GAP],
        }
    }
}

/// One alignment-only symbol found in a description, with the text that states it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentOnlySymbol {
    /// The daggered symbol itself — `X`, `x` or `-`.
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
    ///
    /// Keyed on the **symbol**, not on the axis, so every diagnostic the DNA
    /// half already emitted is byte-identical: only the lower-case spelling —
    /// which no axis could state before `#1715` — cites the RNA ground.
    #[must_use]
    pub fn clause(&self) -> &'static str {
        if self.symbol == MASKED_NUCLEOTIDE_RNA {
            "background/standards.md:39 read through recommendations/general.md:50"
        } else {
            "background/standards.md:39"
        }
    }

    /// The alphabet clause that independently excludes the symbol — the DNA
    /// CAPITALS bullet or its RNA lower-case sibling.
    #[must_use]
    pub fn alphabet_clause(&self) -> &'static str {
        if self.symbol == MASKED_NUCLEOTIDE_RNA {
            "general.md:50"
        } else {
            "general.md:48"
        }
    }

    /// The IUPAC symbol for an unknown base, spelled for the axis that stated
    /// this one — so the repair a diagnostic suggests is itself conformant.
    #[must_use]
    pub fn unknown_base(&self) -> char {
        if self.symbol == MASKED_NUCLEOTIDE_RNA {
            'n'
        } else {
            'N'
        }
    }
}

/// The first alignment-only symbol stated by `text`, if any.
fn symbol_in(text: &str, table: SymbolTable) -> Option<char> {
    let daggered = table.daggered();
    text.chars().find(|c| daggered.contains(c))
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
fn symbol_in_inserted(
    sequence: &InsertedSequence,
    table: SymbolTable,
) -> Option<AlignmentOnlySymbol> {
    match sequence {
        InsertedSequence::Named(name) => symbol_in(name, table).map(|symbol| AlignmentOnlySymbol {
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
fn symbol_in_edit(edit: &NaEdit, table: SymbolTable) -> Option<AlignmentOnlySymbol> {
    match edit {
        NaEdit::Insertion { sequence }
        | NaEdit::BreakpointInsertion { sequence }
        | NaEdit::Delins { sequence, .. }
        | NaEdit::DupIns { sequence } => symbol_in_inserted(sequence, table),
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
/// **This is where the axis is resolved**, once, into the [`SymbolTable`] the
/// axis's descriptions are written in. `r.` takes the RNA table (`general.md:50`,
/// lower case) and every other nucleic-acid axis takes the DNA one
/// (`general.md:48`, CAPITALS) — including `n.`, which addresses a non-coding
/// **DNA** reference sequence. An allele's members are full variants, so each
/// carries its own axis through this same match rather than inheriting the
/// group's.
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
        HgvsVariant::Genome(v) => v
            .loc_edit
            .edit
            .inner()
            .and_then(|edit| symbol_in_edit(edit, SymbolTable::Dna)),
        HgvsVariant::Cds(v) => v
            .loc_edit
            .edit
            .inner()
            .and_then(|edit| symbol_in_edit(edit, SymbolTable::Dna)),
        HgvsVariant::Tx(v) => v
            .loc_edit
            .edit
            .inner()
            .and_then(|edit| symbol_in_edit(edit, SymbolTable::Dna)),
        HgvsVariant::Rna(v) => v
            .loc_edit
            .edit
            .inner()
            .and_then(|edit| symbol_in_edit(edit, SymbolTable::Rna)),
        HgvsVariant::Mt(v) => v
            .loc_edit
            .edit
            .inner()
            .and_then(|edit| symbol_in_edit(edit, SymbolTable::Dna)),
        HgvsVariant::Circular(v) => v
            .loc_edit
            .edit
            .inner()
            .and_then(|edit| symbol_in_edit(edit, SymbolTable::Dna)),
        HgvsVariant::Allele(allele) => allele.variants.iter().find_map(alignment_only_symbol),
        HgvsVariant::GenomeRing(ring) => ring.segments.iter().find_map(|segment| {
            segment
                .edit
                .inner()
                .and_then(|edit| symbol_in_edit(edit, SymbolTable::Dna))
        }),
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

    /// The `r.` half — `#1715`, and the reason [`SymbolTable`] exists.
    ///
    /// On the `r.` axis the ordinary sequence alphabet is lower case
    /// (`general.md:50`), so `r.10delinsacgux` is the same wholesale
    /// reclassification `g.10delinsACGTX` was. Its predecessor
    /// `a_lowercase_rna_masked_nucleotide_is_still_accepted` pinned the miss
    /// while the question was still live; this is the guard that replaces it.
    ///
    /// **Rests on `rulings[rna-axis-alignment-only-symbol-reach]`**, ruled
    /// 2026-08-12: refuse, `standards.md:47`–`:61` governing.
    #[test]
    fn a_lowercase_rna_masked_nucleotide_is_found() {
        let found = symbol_of("NM_TEST.1:r.10delinsacgux").expect("x is stated");
        assert_eq!(found.symbol, 'x');
        assert_eq!(found.stated, "acgux");
        assert_eq!(found.meaning(), "masked nucleotide");
        // The citation is the RNA ground, not the DNA one — a diagnostic that
        // sent an `r.` author to `general.md:48`'s CAPITALS bullet would be
        // telling them to do the one thing their axis forbids.
        assert_eq!(found.alphabet_clause(), "general.md:50");
        assert_eq!(found.unknown_base(), 'n');
        assert!(
            found.clause().contains("general.md:50"),
            "{}",
            found.clause()
        );
    }

    /// Every reachable offset, by construction. One stray `x` reclassifies the
    /// whole run, so — exactly as for uppercase `X` — the check must test for
    /// presence rather than for a lone letter.
    ///
    /// A **leading** `x` is absent from this list on purpose: it does not reach
    /// the AST at all (the grammar refuses it, like `-`), which is measured in
    /// `tests/it/issue_1715_rna_alignment_symbol_reach.rs`. Asserting it here
    /// would be asserting a grammar accident through the wrong door.
    #[test]
    fn an_embedded_lowercase_masked_nucleotide_is_found_at_every_reachable_offset() {
        for (input, stated) in [
            ("NM_TEST.1:r.10delinsacgux", "acgux"),
            ("NM_TEST.1:r.10delinsacxgu", "acxgu"),
            ("NM_TEST.1:r.10delinsax", "ax"),
            ("NM_TEST.1:r.10delinsaxa", "axa"),
            ("NM_TEST.1:r.10_11insacgux", "acgux"),
            ("NM_TEST.1:r.10delinsAcgux", "Acgux"),
        ] {
            let found = symbol_of(input).unwrap_or_else(|| panic!("{input}: x is stated"));
            assert_eq!(found.symbol, 'x', "{input}");
            assert_eq!(found.stated, stated, "{input}");
        }
    }

    /// A member inside an `r.` composite is reached, on the member's own axis.
    /// An allele's members are full variants, so this is what fails if the axis
    /// is ever resolved from the group instead of from the member.
    #[test]
    fn a_member_of_any_rna_composite_spelling_is_reached() {
        for input in [
            "NM_TEST.1:r.[10delinsacgu;20delinsacgux]",
            "NM_TEST.1:r.[10delinsacgux];[20del]",
            "NM_TEST.1:r.(10delinsacgux)",
        ] {
            let found = symbol_of(input).unwrap_or_else(|| panic!("{input}: x is stated"));
            assert_eq!(found.symbol, 'x', "{input}");
            assert_eq!(found.stated, "acgux", "{input}");
        }
    }

    /// The uppercase spelling on the `r.` axis was already matched before
    /// `#1715` and its diagnostic must not move: it is refused on the DNA
    /// table's own clause, which is true of it whichever way the RNA question is
    /// ruled. Pinned so the RNA arm cannot quietly re-cite it.
    #[test]
    fn an_uppercase_masked_nucleotide_on_the_rna_axis_keeps_its_dna_citation() {
        let found = symbol_of("NM_TEST.1:r.10delinsX").expect("X is stated");
        assert_eq!(found.symbol, 'X');
        assert_eq!(found.clause(), "background/standards.md:39");
        assert_eq!(found.alphabet_clause(), "general.md:48");
        assert_eq!(found.unknown_base(), 'N');
    }

    /// **The constraint on the `r.` arm, and the reason it is keyed on the
    /// symbol rather than on the alphabet.** A closed-alphabet check — "every
    /// character must be one of the fifteen RNA symbols" — would refuse all of
    /// these, and whether a named mobile element is conformant at all is a
    /// separate question `DNA/complex.md:169` raises and no ruling settles.
    ///
    /// `alu` is the sharp one: it is lower case *and* reaches `Named`, so it is
    /// exactly what an alphabet rule would catch and a symbol rule must not.
    #[test]
    fn a_named_element_without_the_symbol_is_untouched_on_the_rna_axis() {
        for input in [
            "NM_TEST.1:r.10delinsAluYb8",
            "NM_TEST.1:r.10delinsLINE1",
            "NM_TEST.1:r.10delinsalu",
            "NM_TEST.1:r.10delinsAlu",
        ] {
            assert_eq!(symbol_of(input), None, "{input} states no daggered symbol");
        }
    }

    /// A pure RNA IUPAC run is `Literal`, never `Named`, so the reclassification
    /// above really is triggered by the stray byte. Without this the RNA arm
    /// could be vacuous — matching nothing because nothing reaches it.
    #[test]
    fn a_pure_rna_iupac_run_states_no_alignment_symbol() {
        assert_eq!(symbol_of("NM_TEST.1:r.10delinsacgu"), None);
        assert_eq!(symbol_of("NM_TEST.1:r.10delinsn"), None);
    }

    /// The DNA-axis half of the same scope line, and the constraint on widening
    /// it: a lowercase letter is an ordinary name character here, because
    /// `general.md:48` puts a DNA description in CAPITALS. **`n.` is on this
    /// side of the line**, not the RNA one — it addresses a non-coding *DNA*
    /// reference sequence — which is the axis-assignment error this pins.
    ///
    /// This is also what fails if the `r.` arm is ever implemented as a global
    /// widening to lower case rather than as an axis-keyed one.
    #[test]
    fn a_lowercase_x_on_a_dna_axis_is_still_accepted() {
        for input in [
            "NC_TEST.1:g.10delinsACGTx",
            "NM_TEST.1:c.10delinsACGTx",
            "NM_TEST.1:n.10delinsACGTx",
            "NC_TEST.1:m.10delinsACGTx",
        ] {
            assert_eq!(
                symbol_of(input),
                None,
                "{input}: `general.md:48` is the CAPITALS bullet, so a lower-case letter is an \
                 ordinary name character on every DNA axis"
            );
        }
    }
}
