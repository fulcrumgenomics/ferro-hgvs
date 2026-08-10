//! Why two `g.` delins that look identical get opposite verdicts (`delins.md:47`).
//!
//! # Read this before re-litigating either row
//!
//! This module exists to stop one question being asked a fourth time. Two real
//! genomic `delins` inputs present the same way — a short span, a short payload,
//! one interior reference base that survives into the result because the payload
//! happens to carry the same nucleotide — and ferro splits both. For one of them
//! that is right and for the other it is wrong, and the discriminator is not
//! visible in the shape of the two strings. It is **net length**.
//!
//! Both rows below were measured on this base against the prepared reference,
//! and the reference bases each argument rests on are asserted here rather than
//! quoted, so a reference that does not carry them fails the test instead of
//! quietly changing what the rows mean.
//!
//! ## Row 1 — BRCA2, EQUAL length. Ferro follows `delins.md:17`.
//!
//! ```text
//! NC_000013.11:g.32340866_32340868delinsATC
//!   -> NC_000013.11:g.[32340866G>A;32340868G>C]
//! ```
//!
//! The span is 3 nt and the payload is 3 nt, so the description denotes a
//! column-for-column reading: reference `G,T,G` at 32340866..32340868 becomes
//! `A,T,C`. Denotation therefore *forces* `ref[32340867] == T` — and it is a `T`
//! (asserted below). The changed columns are {32340866} and {32340868}; the
//! column between them is unchanged, so the separation is one nucleotide.
//!
//! Three clauses are in play and all three point the same way:
//!
//! * `delins.md:16` cannot apply — it needs the changed positions to be
//!   consecutive, and these are not.
//! * `delins.md:17` — "two variants separated by one or more nucleotides should
//!   be described individually and **not** as a `delins`" — recommends the
//!   individual description. It is lowercase prose, so read strictly it
//!   *requires* nothing (see the repository `CLAUDE.md` on RFC 2119 keywords in
//!   this spec). What makes it decisive here is not its strength but that it is
//!   the only one of the three whose preconditions this shape actually meets.
//! * `delins.md:18` / `general.md:35`'s codon exception ("two variants separated
//!   by one nucleotide, together affecting one amino acid") cannot reach a `g.`
//!   description at all: a genomic description has no reading frame, so there is
//!   no amino acid for the exception to be about.
//!
//! Nor does the decided ruling below reach this row: it is scoped to a split
//! that exists only because payload bases *coincide* with reference bases, and
//! at equal length the interior column is retained by denotation, not by
//! coincidence.
//!
//! So ferro's output is the one the spec's own recommendation points at, and it
//! is the *input* that departs from it. What this row pins is therefore not that
//! the spanning form is forbidden — no clause here forbids anything — but that
//! merging this row back into one would contradict the only clause that reaches
//! it. That needs a ruling; it is not a confluence fix.
//!
//! ## Row 2 — MSH2, UNEQUAL length. Ferro VIOLATES a decided ruling.
//!
//! ```text
//! NC_000002.11:g.47639670_47639673delinsTT
//!   -> NC_000002.11:g.[47639670_47639671del;47639673G>T]
//! ```
//!
//! The span is 4 nt and the payload is 2 nt — net **-2**. There is no
//! column-for-column reading available, so nothing about the input forces any
//! interior base to be retained. Reference `A,G,T,G` at 47639670..47639673
//! becomes `T,T`, and the `T` at 47639672 survives into ferro's output *only*
//! because the payload happens to carry a `T` that an aligner can pair it with.
//!
//! That coincidence is precisely the construction `delins.md:46` builds — "parts
//! of the inserted sequence *align* with the reference sequence, giving an
//! alternative description" — and `delins.md:47` answers it: **"The `delins`
//! format is recommended"**. Ferro emits the alignment-driven split that `:47`
//! advises against.
//!
//! **The ruling is `delins-merge-vs-individual-gap-two-or-more`**, in
//! `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`, status
//! **`decided`** (operator ruling, 2026-08-07, governing `delins.md:47`). Its
//! scope paragraph is what this row falls under, and it is worth quoting rather
//! than paraphrasing, because the ruling is deliberately narrow:
//!
//! > This record settles ONLY the case `:44-47` describes: a MINIMAL single
//! > `delins` that would be split because payload bases coincide with reference
//! > bases — the alignment-driven split `:46` constructs and `:47` advises
//! > against. It is NOT a general licence to merge changes separated by two or
//! > more nucleotides.
//!
//! `g.47639670_47639673delinsTT` is minimal (reference `A,G,T,G` against payload
//! `T,T` shares neither a prefix nor a suffix, so it cannot be trimmed), and the
//! split exists only because of the payload/reference coincidence. That is the
//! ruling's case.
//!
//! ### One boundary caveat, recorded so the next reader does not trip on it
//!
//! The ruling's **id and question** are framed at a separation of "two or more"
//! nucleotides. This row's separation is **one** — ferro's two members are
//! separated by the single unchanged column 47639672. So the row matches the
//! ruling's *scope paragraph* (an alignment-driven split off a minimal delins)
//! while sitting outside its *title*. The scope paragraph is the operative text
//! — it is what the 2026-08-07 operator ruling wrote to bound the decision —
//! but anyone acting on this row should know the two readings are not identical
//! and say which one they are using.
//!
//! ### Why "it is undecided" is not an available answer any more
//!
//! The ruling's own rationale records that **three different answers were live**
//! in the project's working record when it was rewritten on 2026-08-07: (1) that
//! record itself, holding no ruling had been made; (2) a campaign document
//! calling it implementor's choice, and therefore not a defect in either
//! direction; (3) a measurement-driven session note holding that `:44-47`
//! governs `:17` empirically because Mutalyzer converges on the merge. The
//! operator ruling closed that, for `:47`, scoped. Citing any of the three
//! superseded positions is citing a record that has been overtaken.
//!
//! # Assert-then-flip
//!
//! Row 2's expectation below is ferro's **current, wrong** answer. When the
//! ruling is implemented, it becomes:
//!
//! ```text
//! NC_000002.11:g.47639670_47639673delinsTT   (unchanged — the spanning delins)
//! ```
//!
//! and this test flips to assert that instead. It is asserted rather than
//! `#[ignore]`d so that the day the behaviour moves, this file is one of the
//! things that has to be read. Note the flip is an output-moving change and owes
//! the release a `Representation-Change:` trailer.
//!
//! # Gating
//!
//! Both rows need real reference bases, so they need `FERRO_MANIFEST`. The gate
//! is the strict form (the one `issue_806_effect_real_residues.rs` uses): an
//! **unset** `FERRO_MANIFEST` is a legitimate skip, but a manifest that is set
//! and unusable — or that does not carry these two loci with these bases — is a
//! **failure**. A test that skips green on a half-present reference is worse
//! than no test, because it reports coverage it does not have.

use std::path::PathBuf;

use ferro_hgvs::reference::ReferenceProvider;
use ferro_hgvs::{parse_hgvs, MultiFastaProvider, Normalizer};

/// The manifest, when the operator opted in.
///
/// Returns `Some` for any non-empty value — even a path that does not exist — so
/// that an explicit opt-in pointing at a broken manifest fails loudly instead of
/// silently skipping. Only an *unset* `FERRO_MANIFEST` is a legitimate skip.
fn manifest_path() -> Option<PathBuf> {
    std::env::var_os("FERRO_MANIFEST").map(PathBuf::from)
}

/// One pinned row: the input, the output measured on this base, and the
/// reference bases the row's spec argument depends on.
struct Row {
    input: &'static str,
    expected: &'static str,
    accession: &'static str,
    /// 1-based inclusive, matching the HGVS span in `input`.
    span: (u64, u64),
    bases: &'static str,
}

/// `delins.md:17` is the only clause reaching this shape: equal length, so the
/// description is column-for-column and the unchanged interior column is
/// described individually — which is what ferro emits.
const BRCA2_EQUAL_LENGTH: Row = Row {
    input: "NC_000013.11:g.32340866_32340868delinsATC",
    expected: "NC_000013.11:g.[32340866G>A;32340868G>C]",
    accession: "NC_000013.11",
    span: (32_340_866, 32_340_868),
    bases: "GTG",
};

/// The decided ruling on `delins.md:47` governs: unequal length, so the retained
/// interior base is an alignment coincidence rather than a denoted identity.
///
/// **Assert-then-flip.** `expected` is the current, wrong output; the right one
/// is the input string unchanged.
const MSH2_UNEQUAL_LENGTH: Row = Row {
    input: "NC_000002.11:g.47639670_47639673delinsTT",
    expected: "NC_000002.11:g.[47639670_47639671del;47639673G>T]",
    accession: "NC_000002.11",
    span: (47_639_670, 47_639_673),
    bases: "AGTG",
};

/// The right answer for [`MSH2_UNEQUAL_LENGTH`] once the ruling is implemented.
///
/// Named rather than left in prose so the flip is a two-line edit and so this
/// string is greppable from the ruling record.
const MSH2_RULING_CONFORMANT: &str = "NC_000002.11:g.47639670_47639673delinsTT";

/// Both rows, side by side. See the module docs for the whole argument.
///
/// They are one test rather than two because the point is the *contrast*: read
/// apart, each looks like an arbitrary call, and the pair has already been
/// re-argued three times on the strength of that.
#[test]
fn equal_and_unequal_length_delins_get_opposite_verdicts() {
    let Some(path) = manifest_path() else {
        eprintln!(
            "delins_equal_vs_unequal_length_discriminator: skipping — FERRO_MANIFEST unset. \
             Set it to the prepared reference's manifest.json to run these two rows."
        );
        return;
    };
    // FERRO_MANIFEST is an explicit opt-in: once set, a missing or unloadable
    // manifest is a failure, not a skip.
    assert!(
        path.is_file(),
        "FERRO_MANIFEST points at a missing manifest: {}",
        path.display()
    );
    let provider = MultiFastaProvider::from_manifest(&path)
        .unwrap_or_else(|e| panic!("from_manifest({}) failed: {e}", path.display()));

    // Assert the bases before normalizing. Every clause argument in the module
    // docs is an argument about these bases; if the reference does not carry
    // them, the expectations below are pinning something else and the right
    // outcome is a red test, not a green one.
    for row in [BRCA2_EQUAL_LENGTH, MSH2_UNEQUAL_LENGTH] {
        let (start, end) = row.span;
        let observed = provider
            .get_sequence(row.accession, start - 1, end)
            .unwrap_or_else(|e| panic!("get_sequence({} {start}..{end}): {e}", row.accession));
        assert_eq!(
            observed.to_ascii_uppercase(),
            row.bases,
            "{}:{start}_{end} must read {} for {} to mean what this module says it means",
            row.accession,
            row.bases,
            row.input
        );
    }

    let normalizer = Normalizer::new(provider);
    let normalize = |input: &str| -> String {
        let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        normalizer
            .normalize(&parsed)
            .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
            .to_string()
    };

    // Row 1 — ferro follows `delins.md:17`, the only clause whose preconditions
    // this shape meets: the equal-length span denotes an unchanged interior
    // column.
    assert_eq!(
        normalize(BRCA2_EQUAL_LENGTH.input),
        BRCA2_EQUAL_LENGTH.expected,
        "equal-length {} is described individually, which is what `delins.md:17` recommends and \
         the only one of the three clauses that reaches this shape — `delins.md:16` needs \
         consecutive columns, and `delins.md:18`'s codon exception cannot reach a `g.` description",
        BRCA2_EQUAL_LENGTH.input
    );

    // Row 2 — WRONG, pinned deliberately. See "Assert-then-flip" above.
    assert_eq!(
        normalize(MSH2_UNEQUAL_LENGTH.input),
        MSH2_UNEQUAL_LENGTH.expected,
        "unequal-length {} is currently split on a payload/reference coincidence, which is the \
         alignment-driven split `delins.md:46` constructs and `delins.md:47` advises against. \
         The decided ruling `delins-merge-vs-individual-gap-two-or-more` says the spanning form \
         wins, i.e. {MSH2_RULING_CONFORMANT}. If this assertion just failed with that string, \
         the defect is FIXED — flip the expectation and declare the representation change",
        MSH2_UNEQUAL_LENGTH.input
    );

    // Guard the flip against being confused with a different move: state
    // explicitly that today's output is not yet the conformant one.
    assert_ne!(
        MSH2_UNEQUAL_LENGTH.expected, MSH2_RULING_CONFORMANT,
        "the pinned output equals the ruling-conformant form — the assert-then-flip note above \
         is stale and this test now pins correct behaviour"
    );
}
