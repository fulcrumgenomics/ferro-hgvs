//! Where the spec **does** rank two legal descriptions, ferro must emit the
//! ranked one.
//!
//! # Why this module exists separately from the corpus modules
//!
//! `insertion_adjacency_corpus` and `copy_range_payload_corpus` pin what ferro
//! currently emits, having checked independently that it denotes the right
//! sequence. That is a regression guard, and it is indifferent to *form*: it
//! would happily freeze a legal-but-non-recommended spelling forever.
//!
//! The expectations here are written from the spec clause instead, and are
//! **not** harvested from ferro. A test in this file going red means ferro emits
//! a form the spec ranks below the one it should — which is a finding, not a
//! fixture that needs updating. Read the clause before changing an expectation.
//!
//! # The line between this module and policy pins
//!
//! Only clauses that actually rank forms belong here. The spec is silent on far
//! more than it settles — it does not rank literal against copy-range payloads
//! (`open-issues.md:77`), does not define member overlap beyond one `del`+`dup`
//! example (`general.md:58`), and says nothing about normalizers rewriting one
//! legal description into another. Those live in the corpus modules labelled
//! POLICY, and must not be relabelled as recommended just because ferro is
//! consistent about them.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, NormalizeConfig, Normalizer};

use crate::common::hg38_window::{base_at, provider, LOCAL_CONTIG};

fn normalize_g(body: &str) -> String {
    let input = format!("{LOCAL_CONTIG}:g.{body}");
    let variant: HgvsVariant = parse_hgvs(&input).expect("parse");
    let normalizer = Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    );
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("{body}: {e}"))
        .to_string()
        .strip_prefix(&format!("{LOCAL_CONTIG}:g."))
        .expect("accession and axis survive rendering")
        .to_string()
}

/// Assert the emitted form, quoting the clause that ranks it.
#[track_caller]
fn assert_recommended(body: &str, expected: &str, clause: &str) {
    let got = normalize_g(body);
    assert_eq!(
        got, expected,
        "\n  input:    {body}\n  emitted:  {got}\n  expected: {expected}\n  clause:   {clause}\n\
         \n  This is a recommended-form pin. Before changing the expectation, read the clause."
    );
}

// ---------------------------------------------------------------------------
// general.md:56 — the prioritisation ladder.
//
//   "when a description is possible according to several types, the preferred
//    description is: (1) substitution, (2) deletion, (3) inversion,
//    (4) duplication, (5) insertion"
// ---------------------------------------------------------------------------

/// A one-for-one replacement is a **substitution**, not a one-base `delins`.
///
/// `DNA/delins.md:12` states it directly and `general.md:56` ranks substitution
/// first, so this is the least ambiguous rung on the ladder.
#[test]
fn a_one_base_replacement_renders_as_a_substitution() {
    assert_eq!(base_at(302), b'T');
    assert_recommended(
        "302delinsA",
        "302T>A",
        "general.md:56 (1) substitution; DNA/delins.md:12",
    );
}

/// A tandem copy is a **duplication**, not an insertion.
///
/// `DNA/insertion.md:17`: "tandem duplications are described as a duplication
/// (`g.123_456dup`), not an insertion (`g.456_457ins123_456`)" — the insertion
/// spelling is marked invalid there. `general.md:56` ranks duplication (4) above
/// insertion (5).
#[test]
fn a_tandem_copy_renders_as_a_duplication_not_an_insertion() {
    assert_eq!(
        (base_at(302), base_at(303), base_at(304)),
        (b'T', b'C', b'T'),
        "case assumes g.302_304 is TCT"
    );
    assert_recommended(
        "304_305insTCT",
        "302_304dup",
        "DNA/insertion.md:17 (tandem copy is a dup, not an ins); general.md:56 (4) before (5)",
    );
}

// ---------------------------------------------------------------------------
// inversion.md:69 — an inverted duplication is `ins<range>inv`.
// ---------------------------------------------------------------------------

/// A short match is a **coincidence**, and is left as a literal.
///
/// `g.306_308` is `ATC`, so an insertion of `GAT` at the junction after 308 does
/// reverse-complement-match its neighbour - but a three-base match happens by
/// chance far too often to be evidence of an inverted duplication, and
/// `inversion.md`'s own definition of an inversion requires "**more than one
/// nucleotide**" without licensing every longer match.
///
/// This guard is what stops the rule firing on noise. It is a **ceiling** on
/// the derivation, and the derivation now exists — so it is a live assertion
/// rather than the trivially-passing placeholder it was written as.
/// `the_coincidence_floor_sits_between_eight_and_nine_bases` pins the other
/// side, where the rule must fire.
#[test]
fn a_short_reverse_complement_match_is_not_treated_as_an_inverted_duplication() {
    assert_eq!(
        (base_at(306), base_at(307), base_at(308)),
        (b'A', b'T', b'C'),
        "case assumes g.306_308 is ATC, whose reverse complement is GAT"
    );
    assert_eq!(normalize_g("308_309insGAT"), "308_309insGAT");
}

/// A single-base match is a coincidence roughly one time in four.
#[test]
fn a_single_base_reverse_complement_match_is_not_an_inverted_duplication() {
    let complement_of_neighbour = match base_at(302) {
        b'A' => "T",
        b'C' => "G",
        b'G' => "C",
        _ => "A",
    };
    let body = format!("302_303ins{complement_of_neighbour}");
    assert_eq!(normalize_g(&body), body);
}

// ---------------------------------------------------------------------------
// delins.md:16 / substitution.md:32 — consecutive changes are one delins.
// ---------------------------------------------------------------------------

/// Two changes at **consecutive** positions become one `delins`.
///
/// `DNA/substitution.md:32`: "changes involving two or more consecutive
/// nucleotides are described as deletion-insertion (delins) so the description
/// `c.[79G>T;80C>T]` is not correct."
#[test]
fn two_consecutive_substitutions_become_one_delins() {
    assert_eq!((base_at(302), base_at(303)), (b'T', b'C'));
    assert_recommended(
        "[302T>A;303C>G]",
        "302_303delinsAG",
        "DNA/substitution.md:32; DNA/delins.md:16",
    );
}

// ---------------------------------------------------------------------------
// general.md:34 — separation of one or more nucleotides stays split.
// ---------------------------------------------------------------------------

/// Two changes **separated by a base** stay individual on the genomic axis.
///
/// `general.md:34`: "two variants separated by one or more nucleotides should be
/// described individually and **not** as a 'delins'." The stated exception —
/// "two variants separated by one nucleotide, together affecting one amino acid"
/// — has no purchase on `g.`, which has no reading frame. The coding-axis
/// carve-out is a disclosed deviation and is scoped to `c.`
/// (ledger: `coding-axis-merges-are-a-disclosed-general-34-deviation`), so the
/// genomic axis must not inherit it.
#[test]
fn two_substitutions_separated_by_one_base_stay_split_on_the_genomic_axis() {
    assert_recommended(
        "[302T>A;304T>G]",
        "[302T>A;304T>G]",
        "general.md:34 (separation >= 1 is described individually; the amino-acid \
         exception is coding-axis only)",
    );
}

// ---------------------------------------------------------------------------
// delins.md:86-89 — a substitution and the insertion abutting it are ONE delins.
//
// The single clause in the corpus that rules directly on an insertion abutting a
// neighbouring member, and the strongest register the spec uses: the split form
// is marked invalid, and the note records that the sentence which had permitted
// it was REMOVED.
// ---------------------------------------------------------------------------

/// A coding transcript whose CDS reproduces the `delins.md:86` local context.
///
/// The spec's example sits at `c.2074..c.2080` of a real transcript; the
/// geometry, not the coordinate, is what the ruling turns on, so it is
/// reproduced at `c.10..c.16` over a short synthetic CDS. Building it by hand
/// (rather than via `SyntheticBuilder`) keeps the CDS bounds and the reading
/// frame visible in this file, since the ruling's exception clause elsewhere is
/// frame-dependent.
///
/// | c. range | bases      | note                        |
/// |----------|------------|-----------------------------|
/// | 1..=3    | `ATG`      | initiator                   |
/// | 4..=9    | `CCCGGG`   | filler                      |
/// | 10..=16  | `CATGACA`  | the `..CAT`G`ACA..` context |
/// | 17..=18  | `CC`       | filler                      |
/// | 19..=21  | `TAA`      | terminator                  |
fn spec_context_transcript() -> MockProvider {
    const UTR5: &str = "CCCCCCCCCC";
    const CDS: &str = "ATGCCCGGGCATGACACCTAA";
    const UTR3: &str = "GTGTGTGTGTGTGTGTGTGT";
    assert_eq!(CDS.len() % 3, 0, "CDS must be whole codons");
    assert_eq!(&CDS[9..16], "CATGACA", "c.10_16 must be the spec's context");

    let sequence = format!("{UTR5}{CDS}{UTR3}");
    let cds_start = (UTR5.len() + 1) as u64;
    let cds_end = (UTR5.len() + CDS.len()) as u64;
    let tx_len = sequence.len() as u64;

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_CTX.1".to_string(),
        Some("CTX_TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(cds_start),
        Some(cds_end),
        vec![Exon::new(1, 1, tx_len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// `[G>A; insTA]` at one position is `delinsATA`.
///
/// `DNA/delins.md:86-89`, verbatim: "The correct description of this variant is
/// `NM_007294.3:c.2077delinsATA`", with `c.[2077G>A;2077_2078insTA]` marked
/// invalid and this note attached — "**NOTE**: the answer was modified, i.e. the
/// addition 'However, since the variant is likely a combination of two other
/// variants, it is acceptable to describe it as `c.[2077G>A;2077_2078insTA]`.'
/// was removed."
///
/// Mirrored at `DNA/substitution.md:93-96` and `RNA/delins.md:68-71`.
#[test]
fn substitution_and_abutting_insertion_merge_per_delins_md_86() {
    let input = "NM_CTX.1:c.[13G>A;13_14insTA]";
    let variant: HgvsVariant = parse_hgvs(input).expect("parse");
    let normalizer = Normalizer::with_config(
        spec_context_transcript(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    );
    let got = normalizer
        .normalize(&variant)
        .expect("delins.md:86 shape must normalize, not be refused")
        .to_string();
    assert_eq!(
        got, "NM_CTX.1:c.13delinsATA",
        "\n  DNA/delins.md:86-89 gives one delins as THE correct description and marks the \
         split form invalid;\n  the permission to keep it split was explicitly removed from \
         the spec."
    );
}

// ---------------------------------------------------------------------------
// complex.md:113 — a zero-separation pair the spec publishes SPLIT.
//
// The counterweight to the clause above, and the reason "merge everything
// adjacent" is not the rule. Kept here so the tension is visible in one place.
// ---------------------------------------------------------------------------

/// A large inversion and the substitution abutting its 3' end stay split —
/// **at ISCN scale**.
///
/// `DNA/complex.md:113` publishes `NC_000006.12:g.[776788_93191545inv;
/// 93191546T>C]` as the recommended description of an ISCN inversion "with
/// substitution at break point": separation zero, two members, correct.
///
/// Asserted at that scale rather than in the local window, and the distinction
/// is not cosmetic. In the window — where the members are within the
/// normalizer's reach — ferro does **not** keep the published shape: it
/// re-derives from the resulting sequence and returns a different partition
/// entirely. That behaviour is pinned, with its sequence checked, as
/// `insertion_adjacency_defects::policy_an_abutting_inv_and_sub_are_repartitioned`.
///
/// Splitting the two is what keeps this test honest: it asserts the thing
/// `complex.md:113` actually witnesses (a published pair at breakpoint scale
/// survives), and does not smuggle in a claim about short-range partitioning
/// that the clause does not make.
#[test]
fn an_iscn_scale_inversion_and_its_abutting_substitution_stay_split() {
    let input = "NC_000006.12:g.[776788_93191545inv;93191546T>C]";
    let parsed: HgvsVariant = parse_hgvs(input).expect("the published form is valid HGVS");
    assert_eq!(
        parsed.to_string(),
        input,
        "DNA/complex.md:113's published description must round-trip unchanged"
    );
}

// ---------------------------------------------------------------------------
// duplication.md:17 + :90 — a payload whose PREFIX is a tandem copy carries a
// duplication, and the duplication must be labelled.

/// A duplication at its own write junction, beside an unrelated insertion,
/// keeps its `dup` label.
///
/// `DNA/duplication.md:17` is the operative clause: "duplication may only be
/// used when the additional copy is **directly 3'-flanking** of the original
/// copy". That test passes here — in `…CCAGGGAGTCGGC | CCAGGGAGTCGGC GGCC |
/// T…` the extra copy abuts its original and the unrelated bases fall 3' of the
/// copy, not between it and the original. `:18` then requires the label.
///
/// `DNA/duplication.md:90-92` publishes exactly this geometry, split:
/// `NC_000001.11(NM_206933.2):c.[675-542_1211-703dup;1211-703_1211-702insGTAAA]`,
/// glossed "a duplication of the sequence from … **followed by** the
/// insertion" — which also fixes the member order pinned here.
///
/// The corpus discriminates by that flanking test, not by whether the *whole*
/// payload is a copy: every published complex insertion whose payload contains
/// a nearby copy is one where the copy is **not** directly 3'-flanking
/// (`DNA/insertion.md:42-43`, `:45-46`, `:86-87`, `DNA/complex.md:169-175`),
/// each glossed with the flanking rule as the stated reason.
///
/// Local `209_221` is `CCAGGGAGTCGGC` and `222` is `T`, so the payload opens
/// with `G` and cannot 3'-shift: this pins form, not shuffling.
#[test]
#[ignore = "#1946: needs the render stage — the canonical spelling must be chosen after the member set is final, which is not possible in the partitioner."]
fn a_duplication_beside_an_unrelated_insertion_keeps_its_dup_label() {
    let span: String = (209..=221).map(|p| base_at(p) as char).collect();
    assert_eq!(span, "CCAGGGAGTCGGC", "precondition: the duplicated span");
    assert_eq!(base_at(222), b'T', "payload opens with G, so no 3'-shift");
    assert_recommended(
        "[209_221dup;221_222insGGCC]",
        "[209_221dup;221_222insGGCC]",
        "DNA/duplication.md:17 the extra copy is directly 3'-flanking, so :18 requires \
         the dup label; :90 publishes this geometry split",
    );
}

/// The split and the merged literal are two spellings of one variant and must
/// settle together — on the one that carries the label.
///
/// This is what makes the pin above safe to hold. Refusing to merge without
/// also re-deriving the literal would leave `221_222insCCAGGGAGTCGGCGGCC`
/// settling apart from the split, which
/// `confluence-gate-is-apply-equality-on-every-determined-axis` forbids.
#[test]
#[ignore = "#1946: needs the render stage — see the sibling pin above."]
fn both_spellings_of_a_duplication_plus_insertion_converge_on_the_split() {
    let from_split = normalize_g("[209_221dup;221_222insGGCC]");
    let from_literal = normalize_g("221_222insCCAGGGAGTCGGCGGCC");
    assert_eq!(
        from_split, from_literal,
        "the two spellings of one duplication-plus-insertion settled apart"
    );
    assert_eq!(from_split, "[209_221dup;221_222insGGCC]");
}

/// A short tandem prefix is **not** promoted to a duplication.
///
/// The recommendations supply no floor — `DNA/duplication.md:5` admits "one or
/// more nucleotides" — so a maximal-prefix rule with no coincidence test would
/// mint a `dup` from a one-base match, which occurs in roughly one payload in
/// four. The floor is ours, not the spec's, and it is the same
/// composition-aware threshold the inverted-copy derivation uses.
#[test]
fn a_short_tandem_prefix_is_not_promoted_to_a_duplication() {
    let out = normalize_g("221_222insCGGCC");
    assert!(
        !out.contains("dup"),
        "a short prefix match is coincidence, not evidence of a duplication; got {out}"
    );
}

/// An insertion whose payload is the reverse complement of the **adjacent**
/// span is spelled `ins<range>inv`, naming where the copy came from.
///
/// `DNA/inversion.md:69`: "the recommendation is to describe inverted
/// duplication using the format `g.122_123ins123_234inv` or
/// `g.234_235ins123_234inv` depending on whether the inverted copy is 5' or 3'
/// of the original copy (reference sequence)." Restated as a Note in five more
/// places, each rejecting the `dupinv` shorthand.
///
/// This is the mirror of a rule ferro already applies: a payload equal to the
/// adjacent span *uninverted* is re-spelled `dup` (`insertion.md:17`, pinned
/// above).
///
/// # Why the span must abut the junction
///
/// `:69` answers a question about an inverted **duplication**, and both of its
/// forms place the source immediately beside the insertion point. Restricting
/// the derivation to an abutting span is what the clause describes, and it
/// avoids a tie-break the spec does not supply.
///
/// # Why this payload is ten bases and not three
///
/// A match between a payload and the reverse complement of its neighbour is
/// only evidence of an inverted duplication when it is unlikely to have happened
/// by chance. `g.212_221` is `GGGAGTCGGC`, so the payload is its reverse
/// complement `GCCGACTCCC`; at that length the match is far below the
/// coincidence floor. `g.222` is `T` and the payload opens with `G`, so the
/// insertion cannot 3'-shift and the assertion is about spelling rather than
/// about where the edit settled.
#[test]
fn an_inverted_copy_of_the_adjacent_span_renders_as_ins_range_inv() {
    let span: String = (212..=221).map(|p| base_at(p) as char).collect();
    assert_eq!(span, "GGGAGTCGGC", "case assumes the reference span");
    assert_eq!(base_at(222), b'T', "payload opens with G, so no 3'-shift");
    assert_recommended(
        "221_222insGCCGACTCCC",
        "221_222ins212_221inv",
        "DNA/inversion.md:69 - an inverted duplication is described ins<range>inv",
    );
}

/// The literal and copy-range spellings of one inverted duplication converge.
///
/// The derivation is from the **reference**, not from the input's spelling, so
/// both authors reach the same string - which is what
/// `canonical-form-choice-when-both-legal` requires and what a
/// preserve-the-input implementation would not give.
#[test]
fn both_spellings_of_an_inverted_duplication_converge_on_the_range_form() {
    let from_literal = normalize_g("221_222insGCCGACTCCC");
    let from_range = normalize_g("221_222ins212_221inv");
    assert_eq!(
        from_literal, from_range,
        "the two spellings of one inverted duplication settled apart"
    );
    assert_eq!(from_literal, "221_222ins212_221inv");
}

/// The **other** form in `inversion.md:69` — the source span starting at the
/// junction, when the inverted copy lands 5' of its original.
///
/// `:69` gives two spellings, "depending on whether the inverted copy is 5' or
/// 3' of the original copy": `g.234_235ins123_234inv`, where the source *ends*
/// at the junction (pinned above), and `g.122_123ins123_234inv`, where the
/// source *starts* at it. Implementing one and not the other would make the
/// spelling depend on which side of the original the copy happened to land, for
/// no reason the clause gives.
///
/// `g.201_210` is `TTTCTGAGCC`, so the payload is `GGCTCAGAAA`; `g.201` is `T`
/// and the payload opens with `G`, so it cannot 3'-shift. The 5'-abutting span
/// does not also match, so this case is unambiguous.
#[test]
fn an_inverted_copy_of_the_span_starting_at_the_junction_renders_as_ins_range_inv() {
    let span: String = (201..=210).map(|p| base_at(p) as char).collect();
    assert_eq!(span, "TTTCTGAGCC", "case assumes the reference span");
    assert_eq!(base_at(201), b'T', "payload opens with G, so no 3'-shift");
    assert_recommended(
        "200_201insGGCTCAGAAA",
        "200_201ins201_210inv",
        "DNA/inversion.md:69 - the inverted copy 5' of its original",
    );
}

/// The derived range form and the literal it replaces resolve to the **same
/// SPDI**, so the re-spelling is a change of notation and not of denotation.
///
/// This is the guard behind the ruling record's downstream-consequence
/// paragraph. That paragraph originally named `hgvs_to_spdi` as a blocker, on
/// the strength of the function not yet understanding a range payload; #1966
/// ("allow ref span insertions when converting to SPDI") taught it the shape,
/// and `src/spdi/convert.rs` now decomposes
/// `InsertedSequence::PositionRangeInv` and reads the span through
/// `read_reference_span`. A record that states a downstream fact with nothing
/// checking it is how that paragraph came to outlive the thing it described, so
/// the fact is pinned here rather than only asserted in prose.
///
/// The provider-bearing entry point is the one pinned: `hgvs_to_spdi_simple`
/// and `hgvs_to_spdi(_, None)` decline any payload that must be read from the
/// reference, exactly as they decline the `del`/`dup`/`delins` short forms, and
/// the normalizer that mints this shape holds a provider by construction.
#[test]
fn the_derived_range_form_and_its_literal_resolve_to_the_same_spdi() {
    use ferro_hgvs::spdi::hgvs_to_spdi;

    let spdi_of = |body: &str| -> String {
        let input = format!("{LOCAL_CONTIG}:g.{body}");
        let variant: HgvsVariant = parse_hgvs(&input).expect("parse");
        hgvs_to_spdi(&variant, &provider())
            .unwrap_or_else(|e| panic!("{body}: {e}"))
            .to_string()
    };

    for (literal, range) in [
        ("221_222insGCCGACTCCC", "221_222ins212_221inv"),
        ("200_201insGGCTCAGAAA", "200_201ins201_210inv"),
    ] {
        // The range form is what the normalizer now emits for the literal.
        assert_eq!(normalize_g(literal), range, "derivation moved");
        assert_eq!(
            spdi_of(range),
            spdi_of(literal),
            "the derived form `{range}` must denote what `{literal}` denotes; \
             a divergence here means the re-spelling changed the sequence"
        );
    }
}

/// **POLICY, not recommended form.** Where `MAX_CHANCE_MATCH_PROBABILITY`
/// actually sits: an eight-base inverted copy is left literal, a nine-base one
/// is derived.
///
/// The spec ranks the *forms* — that is what the rest of this module is for —
/// but it supplies no minimum length for an inverted duplication, so the
/// coincidence floor is a house choice. The record says so in capitals, and
/// warns against re-tuning the constant "to whatever makes a test pass". Nothing
/// enforced that: **weakening the threshold 85x, from `1/85_642` to `1/1_000`,
/// left all 41 tests in the three dedicated corpora green.** The existing
/// negatives are one and three bases, and at three bases in this window the
/// chance-match probability is about `1.5e-2` — so anything up to roughly `1/70`
/// still passed them. The suite pinned the threshold only from *above*, via the
/// ten-base positives, and not at all from *below*, which is the side that mints
/// spurious inverted duplications and the whole reason the floor exists.
///
/// So this pins the **boundary itself**, which is the constant's only observable
/// content. Both halves are load-bearing and neither may be relaxed alone:
/// loosening the constant moves `len = 8` into the range form, tightening it
/// moves `len = 9` out.
///
/// The spans are contiguous 5' of the junction at `g.221`, so this walks one
/// length step across the floor with everything else held fixed:
///
/// ```text
/// len=8  span=GAGTCGGC   payload=GCCGACTC   -> 221_222insGCCGACTC   (literal)
/// len=9  span=GGAGTCGGC  payload=GCCGACTCC  -> 221_222ins213_221inv (derived)
/// ```
#[test]
fn the_coincidence_floor_sits_between_eight_and_nine_bases() {
    let span_of =
        |len: u64| -> String { ((222 - len)..=221).map(|p| base_at(p) as char).collect() };
    let payload_of =
        |len: u64| -> String { ferro_hgvs::sequence::reverse_complement(&span_of(len)) };

    assert_eq!(span_of(8), "GAGTCGGC", "case assumes the reference span");
    assert_eq!(span_of(9), "GGAGTCGGC", "case assumes the reference span");

    // Below the floor: the match is coincidence, and the literal stands.
    let eight = payload_of(8);
    assert_eq!(
        normalize_g(&format!("221_222ins{eight}")),
        format!("221_222ins{eight}"),
        "an eight-base reverse-complement match is below the coincidence floor \
         and must stay literal; if this now derives, the floor was loosened"
    );

    // At the floor: the match is evidence, and the range form is derived.
    let nine = payload_of(9);
    assert_eq!(
        normalize_g(&format!("221_222ins{nine}")),
        "221_222ins213_221inv",
        "a nine-base reverse-complement match must clear the coincidence floor; \
         if this now stays literal, the floor was tightened"
    );
}
