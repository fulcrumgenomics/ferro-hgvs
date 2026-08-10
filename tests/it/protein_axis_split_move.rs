//! The protein-axis split move: an **equal-length** `delins` whose interior
//! leaves a residue unchanged is two or more separate changes, not one `delins`.
//!
//! # The clause
//!
//! `recommendations/protein/delins.md:21` states it directly — "two variants
//! separated by one or more amino acids should be described individually and not
//! as a `delins`" — and `:62-64` publishes the worked answer,
//! `p.[Ser44Arg;Trp46Arg]`, with the note that "the variant is **not** described
//! as `p.Ser44_Trp46delinsArgLeuArg`". Unlike the DNA/RNA axes there is no codon
//! carve-out to compete with it: `general.md:35-38` is a *reading-frame*
//! exception, and a protein description has no codons of its own.
//!
//! The spec's own two rows (W48, W59) are executed against the hermetic spec
//! fixtures in `spec_worked_example_rules.rs`. This file covers the parts of the
//! move that the spec does not print an example for — the boundary at which it
//! must **decline** — because that is where a split move goes wrong.
//!
//! # The negative half is the point
//!
//! Every case below that must not split names the exact string ferro must not
//! emit. A positive-only suite would go green on an implementation that split
//! everything in reach: the unequal-length case, the stop-codon case and the
//! no-reference case all *look* like the W48 shape and all have a different
//! right answer.
//!
//! Four declines carry an argument rather than a convention:
//!
//! - **Unequal length.** A `delins` whose span and payload are the same length
//!   has exactly one residue-wise correspondence, so an interior unchanged
//!   residue is a fact about the two sequences. An unequal-length one has none:
//!   locating a run inside it means first choosing an alignment, and the
//!   reference and payload do not determine which. The nucleotide axis settles
//!   that by applying the edits and re-deriving from the resulting sequence;
//!   the protein axis has no apply-to-reference, so there is nothing to derive
//!   the answer from and ferro declines rather than invent one.
//! - **No protein reference.** The unchanged middle residue is a claim about the
//!   *reference*. The payload cannot testify to it, so an accession the provider
//!   cannot serve is declined outright rather than split on the assumption that
//!   `payload[i]` is what the reference holds.
//! - **A `Ter` in the reference span.** Replacing the stop is a stop-loss, which
//!   `protein/extension.md:18` ranks first — "prioritisation: (1) extension, (2)
//!   frameshift or deletion-insertion" — and `:30` spells `p.Ter110GlnextTer17`.
//!   Splitting would render it as a plain substitution and lose the extension.
//! - **A `Ter` in the payload.** `protein/delins.md:47` publishes
//!   `p.(Pro578_Lys579delinsLeuTer)`, so a stop-introducing change stays a
//!   `delins` carrying its `Ter`, and the merge pass converges *toward* that
//!   form. An interior `Ter` is not even authorable — `:45` ("amino acids after
//!   the translation termination codon are **not** listed") is enforced by the
//!   parser — so the whole family is left alone.
//! - **A changed residue 1.** The member would render as `p.Met1Xxx`, which
//!   `protein/substitution.md:49` forbids and the parser refuses, so the split
//!   would emit a description ferro cannot read back. Unlike the three above,
//!   this illegality is *created* by the split — the range input parses fine.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// The accession every case is written against.
const NP: &str = "NP_003997.1";

/// An unrelated accession the provider also carries, so `has_protein_data()` is
/// true even in the cases where `NP` itself is absent (#1131's shape).
const OTHER: &str = "NP_999999.9";

/// A protein whose residues 44, 45, 46 are `Ser`, `Leu`, `Trp` and whose 47 is
/// `Glu` — the spec's W48 context, extended by one residue so a two-residue
/// changed run can be exercised beside a single-residue one.
///
/// Residues 48..=60 are all `Ala`, which makes the sequence **60** residues long
/// and residue 60 an `Ala`. That is load-bearing rather than padding: the allele
/// cases below use `Ala60Gly` as the companion member sitting away from the
/// delins, and the reference residue an allele member names is only checked when
/// it is **in range** — a position past the end is passed through with no
/// mismatch error and no warning. At 57 residues the `Ala60Gly` member was
/// therefore never resolved against the reference at all, so `input == output`
/// held without the companion member ever having been verified. Extending to 60
/// makes the same assertion measure what it reads as: `Ala60Gly` now resolves,
/// its `Ala` matches, and a wrong reference residue there would fail loudly
/// (`Amino acid mismatch at position 60`) instead of passing silently.
fn spec_w48_protein() -> String {
    let mut seq = String::from("M");
    seq.push_str(&"A".repeat(42)); // residues 2..=43
    seq.push_str("SLWE"); // residues 44, 45, 46, 47
    seq.push_str(&"A".repeat(13)); // residues 48..=60 — 60 is the `Ala` the allele cases name
    seq
}

/// A provider carrying `NP`'s protein sequence.
fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein(NP, spec_w48_protein());
    provider
}

/// A protein whose residue 46 is the stop (`Ter`), residues 44 and 45 being
/// `Ser` and `Leu` as above. Replacing that stop is a stop-loss, so this is the
/// reference-side counterpart of the payload-`Ter` fixture.
fn stop_at_46_protein() -> String {
    let mut seq = String::from("M");
    seq.push_str(&"A".repeat(42)); // residues 2..=43
    seq.push_str("SL*"); // residues 44, 45, 46 — 46 is the stop
    seq
}

/// A provider carrying a protein whose span ends at the stop codon.
fn provider_with_a_reference_stop() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein(NP, stop_at_46_protein());
    provider
}

/// A provider with protein data for a *different* accession, so
/// `has_protein_data()` is true while `NP` cannot be fetched.
fn provider_without_the_accession() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein(OTHER, "MAAAAKGSTVE");
    provider
}

/// Normalize `input` through `provider`, rendered.
fn normalize_with(provider: MockProvider, input: &str) -> String {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    Normalizer::new(provider)
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input:?}: {e}"))
        .to_string()
}

/// Assert `input` normalizes to `expected` **and** that `expected` is a fixed
/// point under the same provider. A split that re-merges on the next pass
/// satisfies the first half and fails the second.
fn splits_to(provider: MockProvider, input: &str, expected: &str) {
    assert_eq!(
        normalize_with(provider.clone(), input),
        expected,
        "for {input:?}"
    );
    assert_eq!(
        normalize_with(provider, expected),
        expected,
        "{expected:?} is not a normalization fixed point"
    );
}

/// Assert `input` is preserved verbatim, naming the split ferro must not make.
fn must_not_split(provider: MockProvider, input: &str, forbidden: &str, clause: &str) {
    let actual = normalize_with(provider, input);
    assert_ne!(
        actual, forbidden,
        "{clause}: {input:?} must not be split into {forbidden:?}"
    );
    assert_eq!(actual, input, "{input:?} must be preserved as authored");
}

// ---------------------------------------------------------------------------
// The positive half — the move fires
// ---------------------------------------------------------------------------

/// `protein/delins.md:64` — the spec's own worked example, on a hand-built
/// reference rather than the spec-fixture provider.
#[test]
fn one_unchanged_residue_splits_into_two_substitutions() {
    splits_to(
        provider(),
        &format!("{NP}:p.Ser44_Trp46delinsArgLeuArg"),
        &format!("{NP}:p.[Ser44Arg;Trp46Arg]"),
    );
}

/// `protein/delins.md:21` says "one **or more** amino acids", so a separation of
/// two is the same rule, not a wider one. Residues 45 and 46 (`Leu`, `Trp`) are
/// both carried through unchanged; only 44 and 47 change.
#[test]
fn two_unchanged_residues_split_the_same_way() {
    splits_to(
        provider(),
        &format!("{NP}:p.Ser44_Glu47delinsArgLeuTrpLys"),
        &format!("{NP}:p.[Ser44Arg;Glu47Lys]"),
    );
}

/// A changed run longer than one residue stays a `delins` — a *smaller* one,
/// over the residues it really claims. `general.md:56` ranks substitution above
/// delins, so a one-residue run becomes a substitution and a two-residue run
/// does not.
#[test]
fn a_multi_residue_run_becomes_a_smaller_delins_beside_a_substitution() {
    splits_to(
        provider(),
        &format!("{NP}:p.Ser44_Glu47delinsArgLeuAlaLys"),
        &format!("{NP}:p.[Ser44Arg;Trp46_Glu47delinsAlaLys]"),
    );
}

/// A predicted input keeps its `( )` — carried onto the whole allele, never onto
/// the members (`protein/alleles.md:34`, and #1063 for the nucleotide analogue).
#[test]
fn a_predicted_delins_splits_inside_its_parentheses() {
    splits_to(
        provider(),
        &format!("{NP}:p.(Ser44_Trp46delinsArgLeuArg)"),
        &format!("{NP}:p.[(Ser44Arg;Trp46Arg)]"),
    );
}

// ---------------------------------------------------------------------------
// The negative half — the move declines, and what it must not emit
// ---------------------------------------------------------------------------

/// An **unequal-length** `delins` has no residue-wise correspondence to read, so
/// there is no unchanged interior residue to find — only an alignment somebody
/// would have to choose. Left exactly as authored.
#[test]
fn an_unequal_length_delins_is_not_split() {
    must_not_split(
        provider(),
        &format!("{NP}:p.Ser44_Trp46delinsArgLeuArgLys"),
        &format!("{NP}:p.[Ser44Arg;Trp46_Trp46delinsArgLys]"),
        "protein/delins.md:21 does not license re-aligning a length-changing delins",
    );
}

/// The reference is what says the middle residue is unchanged. With no protein
/// data at all the question cannot be answered, and guessing from the payload is
/// exactly the error the rule must not make.
#[test]
fn a_provider_with_no_protein_data_declines_rather_than_guessing() {
    must_not_split(
        MockProvider::new(),
        &format!("{NP}:p.Ser44_Trp46delinsArgLeuArg"),
        &format!("{NP}:p.[Ser44Arg;Trp46Arg]"),
        "the unchanged middle residue is a claim about the reference",
    );
}

/// `has_protein_data()` is a **provider-wide** flag (#1131), so a provider that
/// carries some other protein still cannot answer for this accession. The fetch,
/// not the flag, is what decides.
#[test]
fn a_provider_missing_this_accession_declines_rather_than_guessing() {
    must_not_split(
        provider_without_the_accession(),
        &format!("{NP}:p.Ser44_Trp46delinsArgLeuArg"),
        &format!("{NP}:p.[Ser44Arg;Trp46Arg]"),
        "a provider-wide capability flag is not a reference for this accession",
    );
}

/// A `Ter` in the payload is declined. `protein/delins.md:47` publishes
/// `p.(Pro578_Lys579delinsLeuTer)` — a stop-introducing change *stays* a `delins`
/// carrying its `Ter` — and `merge::coalesce_protein_adjacent_changes` merges
/// toward exactly that form, admitting a run whose last member is the earliest
/// `Ter` (#1129). Splitting the same shape apart would have one axis pulling both
/// ways at once, so the split move stays out of it.
#[test]
fn a_payload_carrying_a_stop_is_declined() {
    must_not_split(
        provider(),
        &format!("{NP}:p.Ser44_Trp46delinsArgLeuTer"),
        &format!("{NP}:p.[Ser44Arg;Trp46Ter]"),
        "protein/delins.md:47",
    );
}

/// And an *interior* `Ter` is not reachable at all: `protein/delins.md:45` — "amino
/// acids after the translation termination codon are **not** listed" — is enforced
/// by the parser, so the shape the split move would have had to reason about
/// cannot be authored. Asserted rather than assumed, because the decline above is
/// only the whole story while this holds.
#[test]
fn an_interior_stop_is_refused_by_the_parser() {
    let input = format!("{NP}:p.Ser44_Trp46delinsTerLeuArg");
    assert!(
        parse_hgvs(&input).is_err(),
        "protein/delins.md:45: {input} lists residues after the stop and must be refused"
    );
    // The control that makes the refusal attributable to the stop's *position*
    // rather than to anything else in the string: the same shape with the `Ter`
    // last carries no residue after it and parses.
    let trailing = format!("{NP}:p.Ser44_Trp46delinsArgLeuTer");
    assert!(
        parse_hgvs(&trailing).is_ok(),
        "{trailing} lists nothing after the stop and must parse — without this the \
         assertion above cannot tell an interior `Ter` from any other refusal"
    );
}

/// A `Ter` in the **reference span** declines too, and for a different reason
/// than the payload case above: replacing the stop is a **stop-loss**.
///
/// `protein/extension.md:18` ranks the descriptions — "prioritisation: (1)
/// extension, (2) frameshift or deletion-insertion" — and `:30` publishes the
/// form, `p.Ter110GlnextTer17`. So the residue that replaces a stop is described
/// with an extension, never a substitution.
///
/// Without this guard the split fires on the reference `Ser-Leu-Ter` against a
/// payload of `Ala-Leu-His` and emits `p.[Ser44Ala;Ter46His]`, which states the
/// wrong consequence: it spells a stop-loss as an ordinary substitution and
/// silently drops the extension. The payload guard cannot catch it — the payload
/// here carries no `Ter` at all — which is why the reference span needs its own.
#[test]
fn a_reference_stop_is_not_split_into_a_substitution() {
    must_not_split(
        provider_with_a_reference_stop(),
        &format!("{NP}:p.Ser44_Ter46delinsAlaLeuHis"),
        &format!("{NP}:p.[Ser44Ala;Ter46His]"),
        "protein/extension.md:18 — a stop-loss is an extension, not a substitution",
    );
}

/// One contiguous changed run is a genuine `delins`. Nothing separates anything
/// here, so the move must not fire — this is the case the whole rule is *not*
/// about, and the one a too-eager splitter would break.
#[test]
fn a_fully_changed_run_stays_one_delins() {
    must_not_split(
        provider(),
        &format!("{NP}:p.Ser44_Trp46delinsArgAlaArg"),
        &format!("{NP}:p.[Ser44Arg;Leu45Ala;Trp46Arg]"),
        "protein/delins.md:19 — consecutive changed residues are one delins",
    );
}

/// Two residues is the smallest `delins` there is and has no interior at all.
/// Guarded so the run scan can never be relaxed into splitting an adjacent pair,
/// which is the merge direction `protein/substitution.md:23` mandates.
#[test]
fn an_adjacent_pair_has_no_interior_to_split() {
    must_not_split(
        provider(),
        &format!("{NP}:p.Ser44_Leu45delinsArgMet"),
        &format!("{NP}:p.[Ser44Arg;Leu45Met]"),
        "protein/substitution.md:23 — adjacent changes are one delins, not two members",
    );
}

// ---------------------------------------------------------------------------
// Where the move runs, and where it deliberately does not
// ---------------------------------------------------------------------------

/// A separated `delins` **inside a cis allele** keeps its delins.
///
/// The move runs at the whole-description seam
/// (`Normalizer::split_protein_separation`), not per allele member: splitting a
/// member would nest one bracket inside another, and flattening the pieces into
/// the outer bracket first needs a decision about where the member's `( )`
/// marker goes (`protein/alleles.md:34` puts it per member, `:90` round the
/// group). Recorded as a limitation rather than left to emit a nested bracket —
/// see that helper for the two measured outputs.
#[test]
fn a_delins_inside_a_cis_allele_is_left_alone() {
    let input = format!("{NP}:p.[Ser44_Trp46delinsArgLeuArg;Ala60Gly]");
    let actual = normalize_with(provider(), &input);
    assert!(
        !actual.contains("];"),
        "a member's split must not nest a bracket inside the allele, got {actual}"
    );
    assert_eq!(
        actual,
        format!("{NP}:p.[Ser44_Trp46delinsArgLeuArg;Ala60Gly]"),
        "the separated delins stays a delins inside an allele"
    );
}

/// The same for a **trans** allele, and here the reason is harder: the flat form
/// `p.[Ser44Arg;Trp46Arg];[Ala60Gly]` is legal HGVS and the DNA axis both emits
/// and reads `g.[A;B];[C]`, but ferro's protein allele grammar accepts only
/// single-member arms — so emitting it would break the re-parse invariant
/// (`FERRO_ASSERT_REPARSE`). The gap is asserted below so this decline stops
/// being necessary the moment the parser closes it.
#[test]
fn a_delins_inside_a_trans_allele_is_left_alone() {
    let input = format!("{NP}:p.[Ser44_Trp46delinsArgLeuArg];[Ala60Gly]");
    assert_eq!(
        normalize_with(provider(), &input),
        input,
        "the arm's split would produce a description ferro cannot read back"
    );
}

/// The control the two allele declines above rest on: `Ala60Gly` is a member the
/// reference really answers for, not a position past the end of the fixture.
///
/// This is not pedantry about a fixture length. An out-of-range protein position
/// is passed through with **no** mismatch error and **no** warning, so an
/// `input == output` assertion holds for a member the reference was never
/// consulted about — which is the same "an assertion that cannot distinguish the
/// failure from the success" shape the three `is_err()` sites in this file were
/// already corrected for. Both halves are asserted, because either alone is
/// satisfiable by the wrong thing: an *in-range* position is what makes the
/// mismatch check run at all, and a matching residue is what makes the declines
/// above measure the delins member rather than an error path.
#[test]
fn residue_sixty_is_a_residue_the_reference_answers_for() {
    assert_eq!(
        spec_w48_protein().len(),
        60,
        "the allele cases name residue 60; a shorter fixture makes them vacuous"
    );
    assert_eq!(
        normalize_with(provider(), &format!("{NP}:p.Ala60Gly")),
        format!("{NP}:p.Ala60Gly"),
        "residue 60 is `Ala`, so the companion member the allele cases carry is sound"
    );
    // The half that proves the reference is *read* there: name the wrong residue
    // and normalization must fail. At 57 residues this returned `Ok` unchanged.
    let wrong = parse_hgvs(&format!("{NP}:p.Gly60Ala")).expect("parses");
    let err = Normalizer::new(provider())
        .normalize(&wrong)
        .expect_err("a wrong reference residue at an in-range position must not pass");
    assert!(
        err.to_string().contains("position 60"),
        "expected a reference mismatch naming position 60, got {err}"
    );
}

/// The parser gap the trans decline rests on, asserted rather than assumed: a
/// protein allele arm with two members is refused, while its DNA counterpart is
/// accepted. If this ever starts parsing, revisit
/// `a_delins_inside_a_trans_allele_is_left_alone`.
#[test]
fn the_protein_allele_grammar_reads_only_single_member_arms() {
    assert!(
        parse_hgvs(&format!("{NP}:p.[Ser44Arg;Trp46Arg];[Ala60Gly]")).is_err(),
        "if `p.[A;B];[C]` now parses, the trans decline above is no longer needed"
    );
    // The control that pins the *arm's member count* as the cause: the same
    // trans description with a single-member first arm parses, so the refusal
    // above is not about the trans form, the accession or these residues.
    assert!(
        parse_hgvs(&format!("{NP}:p.[Ser44Arg];[Ala60Gly]")).is_ok(),
        "a single-member protein arm must parse — without this the assertion above \
         cannot tell a multi-member arm from any other refusal"
    );
    assert!(
        parse_hgvs("NC_000001.11:g.[100A>T;102A>T];[200A>T]").is_ok(),
        "the DNA axis reads a multi-member arm, which is why this is a protein-only gap"
    );
}

/// And the whole-description case still splits, which is what makes the two
/// declines above a *scope* rather than a regression.
#[test]
fn the_whole_description_case_still_splits() {
    assert_eq!(
        normalize_with(provider(), &format!("{NP}:p.Ser44_Trp46delinsArgLeuArg")),
        format!("{NP}:p.[Ser44Arg;Trp46Arg]")
    );
}

/// A changed residue 1 is declined, because splitting it emits a description
/// the parser refuses.
///
/// `p.Met1_Ala3delinsValAlaTrp` has the W48 shape exactly — equal length, an
/// unchanged interior residue, two changed runs — so every other gate passes it
/// through. But residue 1 is the translation initiation codon, and the member it
/// would produce is `p.Met1Val`: `protein/substitution.md:49` says a variant
/// changing the initiator is "not described as a substitution or as an
/// extension", and `parse_hgvs` enforces that.
///
/// The input itself parses, because `validate_no_start_loss_substitution` keys
/// on a single certain position and this is a range. So the illegality is
/// **created by the split**, which is why the guard has to sit in the split
/// rather than being inherited from the parser — and why the assertion below is
/// on the *re-parse*, not just on the string.
///
/// Ferro declines rather than choosing among `p.0` / `p.0?` / `p.(Met1?)` /
/// `p.Met1_Leu2ins…`: which of those is right is a claim about the consequence,
/// and two protein sequences do not settle it.
#[test]
fn a_changed_initiation_codon_is_left_alone() {
    let input = format!("{NP}:p.Met1_Ala3delinsValAlaTrp");
    let out = normalize_with(provider(), &input);
    assert_eq!(
        out, input,
        "splitting here would emit `p.[Met1Val;Ala3Trp]`, which is not a legal \
         description of a start-loss (`protein/substitution.md:49`)"
    );
    assert!(
        parse_hgvs(&out).is_ok(),
        "the output must re-parse; {out:?} did not"
    );
    assert!(
        parse_hgvs(&format!("{NP}:p.[Met1Val;Ala3Trp]")).is_err(),
        "if the split form now parses, this decline needs revisiting rather than keeping"
    );
    // The control that pins `Met1` as the cause rather than the `p.[A;B]`
    // bracket form: the same two-member cis allele off the initiator parses.
    assert!(
        parse_hgvs(&format!("{NP}:p.[Ala2Val;Ala3Trp]")).is_ok(),
        "a two-member cis protein allele away from residue 1 must parse — without this \
         the assertion above cannot tell a start-loss substitution from a refusal of \
         the bracket form itself"
    );
}

/// The guard is scoped to residue 1 *changing*, not to any span that reaches it.
///
/// Here residue 1 is unchanged and the changed runs are at 2 and 4, so no member
/// lands on the initiator and the split proceeds. Without this, a guard written
/// as "decline any span starting at 1" would look correct and silently stop
/// splitting a family it has no business declining.
#[test]
fn a_span_reaching_residue_one_still_splits_when_residue_one_is_unchanged() {
    splits_to(
        provider(),
        &format!("{NP}:p.Met1_Ala4delinsMetGlyAlaTrp"),
        &format!("{NP}:p.[Ala2Gly;Ala4Trp]"),
    );
}
