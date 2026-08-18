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
//! ## Row 2 — MSH2, UNEQUAL length. Since #2155, ferro FOLLOWS `delins.md:47`.
//!
//! ```text
//! NC_000002.11:g.47639670_47639673delinsTT
//!   -> NC_000002.11:g.47639670_47639673delinsTT
//! ```
//!
//! (Unchanged: the spanning form is already what the input names. See the
//! 2026-08-17 correction section below for why this row moved and what it
//! used to read.)
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
//! format is recommended"**. Since #2155, ferro emits exactly that spanning
//! form — see the 2026-08-17 correction section below; before that date it
//! emitted the alignment-driven split `:47` advises against.
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
//! ruling's case — and, since #2155 widened the axis the ruling's carve-out
//! reaches, it is now honoured on this `g.` row too.
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
//! # The assert-then-flip is WITHDRAWN (#1835) — row 2 was never going to merge
//!
//! This section used to say that row 2's expectation was ferro's "current,
//! **wrong**" answer, and that when `delins-merge-vs-individual-gap-two-or-more`
//! was implemented the row would become the spanning
//! `NC_000002.11:g.47639670_47639673delinsTT`. **That instruction is wrong and
//! must not be followed.** It was written before the axis scope was decided.
//!
//! `delins-payload-coincidence-carve-out-is-coding-dna-scoped` (decided) scopes
//! `DNA/delins.md:47` — the clause that recommends the spanning form over a
//! payload-coincidence split — to the **coding DNA axis**. Row 2 is
//! `NC_000002.11:g.`, a genomic reference. So `:47` does not reach it, and
//! `general.md:34` governs unqualified: two variants separated by an unchanged
//! nucleotide are described individually. Row 2's split is therefore
//! **conformant**, and it always was; what was wrong was the note, not the
//! output.
//!
//! Note the direction scope points the same way about *which* record applies but
//! not about the outcome: this row is a net deletion (4 nt span, 2 nt payload),
//! so `delins-merge-vs-individual-gap-two-or-more` would reach it by direction —
//! the AXIS scope is what excludes it, and the axis scope is the later and
//! narrower of the two. Read them in that order.
//!
//! **The row's SPLIT nevertheless moved**, which is a separate and much smaller
//! fact — as of the history through #1835, before the correction below: the
//! partition default flip re-derives the block from the resulting sequence, so
//! the cut lands at a different place inside the same span. Both the old and
//! the new spelling were two members separated by an unchanged base, so
//! `general.md:34` was satisfied by each and `canonical-form-choice-when-both-legal`
//! selected between them by derivation. That was an output-moving change and
//! owed the release a `Representation-Change:` trailer, which #1835 carried.
//!
//! **This section's own claim that "row 2 was never going to merge" is now
//! superseded — see the 2026-08-17 correction below.** It was true of the
//! ruling as scoped through #1835; it stopped being true the day the axis scope
//! widened.
//!
//! # CORRECTED 2026-08-17 (#2155): row 2 DOES merge now — the axis scope widened
//!
//! `delins-payload-coincidence-carve-out-is-coding-dna-scoped` was **superseded
//! to all DNA axes** on 2026-08-17 (#2155): `:47`'s carve-out now reaches
//! `c./g./m./n.`, not `c.` alone. Row 2 is `NC_000002.11:g.`, one of the axes
//! newly in reach. So both paragraphs above — "So `:47` does not reach it" and
//! "Row 2's split is therefore conformant" — describe the ruling as it stood
//! before that date, not as it stands now.
//!
//! # #1610 predicted this outcome BEFORE it happened, and the axis is why it now does
//!
//! `unequal-length-block-a-placed-gap-is-not-a-separation` (decided) keeps a
//! lone unequal-length net-deletion `delins` whole where the derived split is a
//! placed gap plus `delins`-rendering members — and row 2 meets every one of its
//! four *shape* conditions: `AGTG` against `TT` is a net deletion, the two
//! derived members sit one unchanged base apart, and that base is a payload `T`
//! landing on a reference `T`.
//!
//! This section used to say row 2 was "nevertheless unmoved", because #1610's
//! rule is gated on `CoincidenceCarveOut::may_disbelieve_a_separation` — the
//! same gate this module already argues from — and that gate answered `false`
//! for `g.`. **It now answers `true`.** #2155 widened the same
//! `CoincidenceCarveOut` the axis scope above names, so #1610's rule reaches
//! row 2's shape and merges it: exactly the shape-match this section already
//! established, now actually taken.
//!
//! An earlier revision of #1610 shipped the rule axis-blind — reaching every
//! axis rather than none off `c.` — and that measured a **rank-1 regression**:
//! `spec_conformance_axis`'s `guard_violations` moved 0 to 5, the rejected
//! SVD-WG010 shape on the frameless axes generally, not scoped to this row's own
//! shape. #2155 is a different change: a **deliberate, adjudicated** widening of
//! the carve-out itself (see `delins-payload-coincidence-carve-out-is-coding-
//! dna-scoped`'s supersession note), not an accidental axis-blind gate. Its own
//! measured cost on the same guard is `guard_violations` 0 -> 12 in both
//! directions (`spec_conformance_axis.rs`), disclosed and accepted as the
//! ruling's own rule-2 deviation — not reverted. So the guard below now asserts
//! the **positive** — that row 2 reaches the spanning form — rather than the
//! negative it asserted through #1835.
//!
//! # Gating, and the hermetic companion that covers what this file cannot
//!
//! Both rows need real reference bases, so they need `FERRO_MANIFEST`, and PR
//! CI provisions none (`.github/workflows/ci.yml` says so and emits a `::notice`
//! about it). So the test below runs in the **nightly** and skips on every pull
//! request. Its PR-time counterpart is
//! `tests/it/delins_equal_vs_unequal_length_hermetic.rs`
//! (`it::delins_equal_vs_unequal_length_hermetic::equal_and_unequal_length_delins_get_opposite_verdicts_hermetically`),
//! which reproduces the same discriminator on a synthetic contig — one span, two
//! payloads differing by a single dropped base, so net length is provably the
//! only thing separating the two verdicts. Read the two together: the hermetic
//! file is the gate, and this file is the reality check that the synthetic shape
//! still matches biology.
//!
//! **The two files' unequal rows no longer land the same way, and that is
//! expected, not a divergence to reconcile.** Row 2 here (`AGTG` -> `TT`)
//! derives a placed gap plus a `delins`-rendering member, meets every one of
//! #1610's shape conditions, and merges under the widened axis. The hermetic
//! file's unequal arm (`ATG` -> `TC`) derives a placed gap plus a plain
//! **substitution**, which #1610's own fourth condition excludes regardless of
//! axis, so it stays split. The hermetic file's own positive control
//! (`the_widened_carve_out_now_merges_a_placed_gap_shape_on_g`) reproduces this
//! row's shape hermetically and merges, which is what actually verifies the
//! widening on every PR — this file only verifies it against real biology, in
//! the nightly.
//!
//! The gate
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

/// Unequal length, and — since #2155 (2026-08-17) — MERGED: `delins.md:47`'s
/// carve-out now reaches this `g.` row, so `general.md:34`'s default no longer
/// governs it unqualified. `expected` is the spanning form, identical to
/// `input`.
///
/// **History, kept because it explains why `expected` changed twice.** Through
/// #1835 the split moved *within the span* when the partition default flipped
/// (from `g.[47639670_47639671del;47639673G>T]` to
/// `g.[47639670_47639671delinsT;47639673del]`) while staying a split — both
/// spellings were two members separated by the unchanged base at 47639672, and
/// `canonical-form-choice-when-both-legal` selected between them by
/// derivation. #2155 is a second, independent move: the axis scope widened, so
/// the split — whichever spelling it was — no longer survives at all. See the
/// module docs' 2026-08-17 correction section for the full argument.
const MSH2_UNEQUAL_LENGTH: Row = Row {
    input: "NC_000002.11:g.47639670_47639673delinsTT",
    expected: "NC_000002.11:g.47639670_47639673delinsTT",
    accession: "NC_000002.11",
    span: (47_639_670, 47_639_673),
    bases: "AGTG",
};

/// The split form `general.md:34` used to require and `delins.md:47` always
/// recommended against, kept as a negative guard now that `expected` **is**
/// the spanning form: the assertion below requires ferro's answer to differ
/// from this, so a regression back to the pre-#2155 split is caught rather
/// than silently re-blessed.
///
/// Named rather than left in prose so it stays greppable from the ruling
/// record. Previously named `MSH2_SPANNING_FORM_OFF_AXIS` and used the other
/// way around (asserted equal to nothing, not-equal to `expected`); the rename
/// and the pinned string both changed to the split spelling #1835 last shipped,
/// since that is now the form to guard *against* rather than the form to guard
/// *for*.
const MSH2_SPLIT_FORM_PRE_2155: &str = "NC_000002.11:g.[47639670_47639671delinsT;47639673del]";

/// Both rows, side by side. See the module docs for the whole argument.
///
/// They are one test rather than two because the point is the *contrast*: read
/// apart, each looks like an arbitrary call, and the pair has already been
/// re-argued three times on the strength of that.
#[test]
fn equal_and_unequal_length_delins_get_opposite_verdicts() {
    let Some(path) = manifest_path() else {
        crate::common::manifest::absent("delins_equal_vs_unequal_length_discriminator");
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

    // Row 2 — since #2155, MERGED. See the module docs' 2026-08-17 correction.
    assert_eq!(
        normalize(MSH2_UNEQUAL_LENGTH.input),
        MSH2_UNEQUAL_LENGTH.expected,
        "unequal-length {} must collapse to the spanning form: the derived split is a placed \
         gap (`47639673del`) plus a `delins`-rendering member (`47639670_47639671delinsT`), \
         which is exactly `unequal-length-block-a-placed-gap-is-not-a-separation`'s (#1610) \
         shape, and since #2155 widened `delins-payload-coincidence-carve-out-is-coding-dna-\
         scoped` to every DNA axis, `:47`'s carve-out now reaches this `g.` row and recommends \
         the spanning form it always would have off that axis restriction. If this just failed \
         with the split form {MSH2_SPLIT_FORM_PRE_2155}, that is a REGRESSION against the \
         widened axis scope, not the pre-#2155 conformant answer an older revision of this \
         message asserted — read the module docs before re-pinning",
        MSH2_UNEQUAL_LENGTH.input
    );

    // The positive guard: ferro must not regress back to the pre-#2155 split.
    assert_ne!(
        MSH2_UNEQUAL_LENGTH.expected, MSH2_SPLIT_FORM_PRE_2155,
        "the pinned output equals the pre-#2155 split form — `delins.md:47`'s carve-out reaches \
         every DNA axis now, so pinning the split here would encode a regression against the \
         widened `delins-payload-coincidence-carve-out-is-coding-dna-scoped`"
    );
}
