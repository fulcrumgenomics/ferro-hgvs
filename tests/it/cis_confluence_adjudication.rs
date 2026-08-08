//! Adjudication records for the divergent cis confluence classes.
//!
//! `tests/it/cis_confluence_axis.rs` measures the pile — 11 272 designed classes,
//! **4 643** of which reach more than one normalized output in the 3' direction,
//! with `declined` and `sequence_changed` both zero. That census says how large
//! the disagreement is and nothing about what the right answer would be. This
//! module carries the answers that have been reached, and names the ones that
//! have not, on concrete inputs rather than on aggregate counts.
//!
//! The shape census behind these records is produced by
//! `examples/dump_confluence_divergences.rs`, which clusters the 4 643 by the
//! *disagreement* — how the competing outputs' arities relate, and what interior
//! gaps the more-split one leaves — rather than by the design parameters. Eight
//! families result; the top four cover **92.6 %**:
//!
//! ```text
//!   3066 (66.0%)  spanning-vs-split/min-2+
//!    696 (15.0%)  spanning-vs-split/min-1
//!    271 ( 5.8%)  split-vs-split/min-2+
//!    267 ( 5.8%)  split-vs-split/min-1
//!    181 ( 3.9%)  same-arity/min-2+
//!    131 ( 2.8%)  same-arity/min-1
//!     28 ( 0.6%)  spanning-vs-split/has-0
//!      3 ( 0.1%)  split-vs-split/has-0
//! ```
//!
//! **These are the post-#1537 figures, and the movement is the point.** Before
//! that PR the census read 4 639 classes over *nine* families with a
//! `same-arity/has-0` arm, and the three `has-0` arms held 104 classes between
//! them. #1537 ("never split a delins into members on consecutive nucleotides")
//! merged the flush members, so the `has-0` arms collapsed to 31 and the
//! `same-arity/has-0` family disappeared entirely. What remains under `has-0` is
//! the carve-out — see `a_dup_flush_against_a_del_is_left_alone`.
//!
//! **One mechanism accounts for 4 568 of the 4 643 (98.4 %)**: the single
//! spanning spelling reaches an output disjoint from every output the
//! multi-member spellings reach. Ferro re-partitions a lone `delins` from the
//! sequence and preserves the member boundaries it is handed otherwise, so two
//! spellings of one variant keep two partitions. Only 205 classes (4.4 %) have
//! two *multi-member* spellings disagreeing with each other, and in every one of
//! those the disagreeing spellings have the same member count — a placement
//! disagreement, not a partitioning one.
//!
//! Four records live here, in **three kinds**. They are not interchangeable
//! (see the repository `CLAUDE.md`), so each test states its own:
//!
//!  - **One adjudicated-correct.**
//!    `two_adjacent_members_that_both_consume_reference_are_one_delins` pins the
//!    exact form a cited clause requires, so it fails if behaviour regresses
//!    away from a decided answer.
//!  - **One scope carve-out.** `a_dup_flush_against_a_del_is_left_alone` pins
//!    the limit of that same ruling. Ferro's output there is not adjudicated
//!    right; it is adjudicated *out of scope*, which is a third thing.
//!  - **Two pinned open disagreements.**
//!    `the_separation_two_members_present_is_not_a_property_of_the_variant`
//!    demonstrates a question its record leaves `undecided`, and
//!    `a_spanning_delins_and_its_aligned_split_are_two_fixed_points` pins the
//!    largest family with the ruling it waits on named, so a future session does
//!    not re-derive it.
//!
//! Nothing in this module asserts a behaviour change: every expectation below
//! is ferro's output on `main` today.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

/// The synthetic transcript every case here is drawn against, taken verbatim
/// from `examples/generate_cis_confluence_corpus.rs`'s first core so the cases
/// are the corpus's own rather than newly invented ones. Its first bases are
///
/// ```text
///   c.  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
///       T T T T T T T T T  A  A  T  A  T  A  T  T
/// ```
///
/// which is what makes the `c.9`/`c.10`/`c.11` cases below hand-checkable: the
/// `A` at `c.10` and the `A` at `c.11` are one two-base run, so deleting either
/// denotes the same sequence.
const CORE: &str = "TTTTTTTTTAATATATTTTAATATAATTAAAAAAATAATTTTTATAAATATATTATTTTAAAAA";

const TX_ACCESSION: &str = "NM_TEST.1";
const TX_CONTIG: &str = "chr_synth";
const PAD_OFFSET: usize = 256;
const CDS_START: u64 = 1;
const CDS_END: u64 = 63;

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let tx_len = CORE.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
    let transcript = Transcript::new(
        TX_ACCESSION.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        CORE.to_string(),
        Some(CDS_START),
        Some(CDS_END),
        vec![exon],
        Some(TX_CONTIG.to_string()),
        Some(g_start),
        Some(g_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    provider.add_genomic_sequence(TX_CONTIG, format!("{pad}{CORE}{pad}"));
    provider.add_transcript(transcript);
    provider
}

fn normalized(provider: &MockProvider, description: &str) -> String {
    let variant = parse_hgvs(description).unwrap_or_else(|e| panic!("{description}: {e}"));
    let normalizer = Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("{description}: {e}"))
        .to_string()
}

/// The sequence a description denotes, established by applying it to the
/// reference through `hgvs_to_spdi` — the corpus's own ground-truth oracle, and
/// the only one in the repository that is independent of the normalizer by
/// construction. Every "these two spellings are one variant" claim below is
/// checked with this rather than asserted.
fn denotes(provider: &MockProvider, description: &str) -> String {
    let members: Vec<HgvsVariant> = match parse_hgvs(description).expect("parses") {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let mut triples: Vec<(usize, String, String)> = members
        .iter()
        .map(|member| {
            let triple = hgvs_to_spdi(member, provider).expect("applies");
            (
                usize::try_from(triple.position).expect("non-negative"),
                triple.deletion.clone(),
                triple.insertion.clone(),
            )
        })
        .collect();
    // 3' to 5', longer deletion first — the `cis_apply_oracle` walk, so a member
    // is never spliced over one that has already been applied.
    triples.sort_by_key(|t| (std::cmp::Reverse(t.0), std::cmp::Reverse(t.1.len())));
    let mut edited = CORE.as_bytes().to_vec();
    for (position, deletion, insertion) in &triples {
        let end = position + deletion.len();
        assert_eq!(
            &CORE.as_bytes()[*position..end],
            deletion.as_bytes(),
            "{description}: stated bases disagree with the reference"
        );
        edited.splice(*position..end, insertion.bytes());
    }
    String::from_utf8(edited).expect("ascii")
}

// ---------------------------------------------------------------------------
// Record 1 — settled: separation is a property of the spelling
// ---------------------------------------------------------------------------

/// **Adjudication.** `general.md:34` and `DNA/delins.md:17` — "two variants
/// separated by one or more nucleotides should be described individually and
/// **not** as a 'delins'" — key on a quantity the *variant* does not carry.
/// Two spellings of one variant can present different separations, so the rule
/// returns two answers for one variant and cannot be evaluated confluently.
///
/// **Authority.** Ruling record
/// `separation-is-a-property-of-the-spelling-not-of-the-variant`, status
/// `undecided` as to which partition should be canonical — this test pins the
/// *demonstration*, which is settled, not the choice, which is not.
///
/// The case: `c.10` and `c.11` are both `A`, so deleting one or the other is the
/// same edit. `c.[9del;10del;13del]` therefore denotes exactly what
/// `c.[9del;11del;13del]` denotes — checked here through the apply oracle, not
/// assumed. But the first spells its members at separations 0 and 2 and the
/// second at 1 and 1, and `general.md:34` merges a separation-0 pair while
/// splitting a separation-1 pair. Ferro follows it, and produces two outputs.
#[test]
fn the_separation_two_members_present_is_not_a_property_of_the_variant() {
    let provider = provider();
    let adjacent = "NM_TEST.1:c.[9del;10del;13del]";
    let apart = "NM_TEST.1:c.[9del;11del;13del]";
    let spanning = "NM_TEST.1:c.9_13delinsAT";

    // One variant: verified by application, independently of the normalizer.
    let sequence = denotes(&provider, adjacent);
    assert_eq!(denotes(&provider, apart), sequence);
    assert_eq!(denotes(&provider, spanning), sequence);

    // Three outputs. The separation each spelling presents is what selects them,
    // and no two spellings present the same one.
    assert_eq!(
        normalized(&provider, adjacent),
        "NM_TEST.1:c.[9_10del;13del]",
        "separations 0 and 2: `delins.md:16` merges the adjacent pair"
    );
    assert_eq!(
        normalized(&provider, apart),
        "NM_TEST.1:c.[9del;11del;13del]",
        "separations 1 and 1: `general.md:34` keeps all three individual"
    );
    assert_eq!(
        normalized(&provider, spanning),
        "NM_TEST.1:c.9_13delinsAT",
        "no members at all, so no separation to read: the spanning form survives"
    );
}

// ---------------------------------------------------------------------------
// Record 2 — settled, and ferro conforms: adjacent members
// ---------------------------------------------------------------------------

/// **Adjudication.** A normalized allele must not leave two members with no
/// unchanged nucleotide between them when both consume reference bases.
/// `DNA/delins.md:16` covers consecutive changed nucleotides — they are one
/// delins — and `DNA/substitution.md:32` marks the split spelling of exactly
/// that shape `class="invalid"`, saying of `LRG_199t1:c.79_80delinsTT` that
/// "the description `c.[79G>T;80C>T]` is not correct". `general.md:34` does not
/// compete: it governs members separated by *one or more* nucleotides and says
/// nothing about members separated by none.
///
/// **Authority.** Ruling record
/// `delins-adjacent-members-when-both-consume-reference`, status `decided`,
/// governing `DNA/substitution.md:32`, deviating from `general.md:56`.
///
/// **Adjudicated-correct**, and it did not start that way. When this record was
/// first written ferro cut `c.10_13delinsTAAT` at the codon boundary between
/// codon 4 (`c.10_12`) and codon 5 (`c.13_15`) and emitted
/// `c.[10_12delinsTAA;13A>T]` — two flush members — so the test asserted the
/// wrong form and labelled it a deviation. **PR #1537 fixed it and closed
/// #1524**, so this now asserts the form the spec states and fails if that
/// regresses. The corpus figure that came with the old framing — all 65 classes
/// with this shape on the `c.` axis, none on `g.` — is gone with it: after
/// #1537 the count in this ruling's scope is **zero**.
///
/// **Scope.** Both adjacent members must consume reference bases; see
/// `a_dup_flush_against_a_del_is_left_alone` for the other side of that line.
#[test]
fn two_adjacent_members_that_both_consume_reference_are_one_delins() {
    let provider = provider();
    let spanning = "NM_TEST.1:c.10_13delinsTAAT";
    let authored = "NM_TEST.1:c.[9dup;13del]";

    assert_eq!(denotes(&provider, spanning), denotes(&provider, authored));

    // `c.10_12` and `c.13` would be adjacent — no unchanged nucleotide between
    // them — so `delins.md:16` asks for the single delins this input already
    // was, and `substitution.md:32` marks the split invalid by name.
    assert_eq!(
        normalized(&provider, spanning),
        "NM_TEST.1:c.10_13delinsTAAT",
        "adjudicated CORRECT against `DNA/substitution.md:32`; fixed by #1537, closing #1524"
    );
    // The other spelling of the same variant still keeps its own partition, so
    // the class remains divergent. That residue is not this record's — it is
    // `separation-is-a-property-of-the-spelling-not-of-the-variant`, undecided.
    assert_eq!(normalized(&provider, authored), "NM_TEST.1:c.[9dup;13del]");
}

// ---------------------------------------------------------------------------
// Record 3 — the scope limit of record 2, pinned on a real corpus class
// ---------------------------------------------------------------------------

/// **Adjudication.** The merge above stops where `DNA/duplication.md:18` starts:
/// "when a variant can be described as a duplication, it **must** be described
/// as a duplication and not as, e.g., an insertion". Merging a `dup` into a
/// flush neighbour destroys the duplication, so the ruling does not reach an
/// adjacency with a `dup` or an `ins` on either side. That conflict is real and
/// is deliberately left open.
///
/// **Authority.** The `SCOPE` paragraph of ruling record
/// `delins-adjacent-members-when-both-consume-reference` — a carve-out from a
/// `decided` record, not a deviation from it. Ferro's output here is not
/// adjudicated right; it is adjudicated *out of scope*.
///
/// The case is corpus class `s00-c-m4-sep3-p8-rot3`, taken verbatim rather than
/// invented, and it is one of the 33 classes that still leave a zero interior
/// gap after #1537. `c.16_23dup` and `c.24_31del` are flush and stay flush.
#[test]
fn a_dup_flush_against_a_del_is_left_alone() {
    let provider = provider();
    let authored = "NM_TEST.1:c.[9T>A;13_20dup;24_31del;34_35insCACCAAAA]";

    assert_eq!(
        normalized(&provider, authored),
        "NM_TEST.1:c.[9T>A;16_23dup;24_31del;34_35insCACCAAAA]",
        "separation 0 between `16_23dup` and `24_31del`, and `DNA/duplication.md:18` \
         competes — outside the scope of the adjacency ruling, not a deviation from it"
    );
}

// ---------------------------------------------------------------------------
// Record 4 — open on the split-input side: the spanning delins against its
// aligned split
// ---------------------------------------------------------------------------

/// **Not settled on the side this test pins, and it names the record that turns
/// on.** This is the largest family in the pile — 3 066 of 4 643 classes,
/// 66.0 % — and it is
/// the spec's own `DNA/delins.md:44-47` shape reproduced synthetically: parts of
/// the inserted sequence align with the reference, giving an alternative
/// description as separate members. `:17` read literally demands the split;
/// `:47` recommends the delins three lines later, in the spec's own worked
/// example.
///
/// **Authority.** Ruling record `delins-merge-vs-individual-gap-two-or-more`,
/// whose committed status is **`decided`** — for `:47`, and *scoped*. It
/// settles only "a MINIMAL single `delins` that would be split because payload
/// bases coincide with reference bases", and says in as many words that it is
/// "NOT a general licence to merge changes separated by two or more
/// nucleotides". Ferro emitting `:47`'s form from the spanning spelling is that
/// ruling working. What is **not** settled is the other side: whether the scope
/// reaches a spelling that arrives already split, so that
/// `c.[9_12del;14_17del]` ought to converge on the delins too. Nothing here
/// decides that — both outputs are pinned as observations, and converging them
/// is a representation change for a stored library. The equivalence class
/// `dna-delins-vs-aligned-split-850-901` is the spec's own instance of the pair.
///
/// (An earlier revision of this comment described that record as `undecided`.
/// It was already `decided` when this test was written; only the status word
/// is corrected here, not the disposition of the case.)
///
/// The premise of `:44-47` — that parts of the payload align with the
/// reference — holds for **every** class in this corpus, not just this family:
/// two spellings can only denote one sequence if the coarser one's payload
/// re-states the finer one's interior gap bases. So the whole 4 643 is one
/// question wearing eight shapes, and this record is where it is named.
#[test]
fn a_spanning_delins_and_its_aligned_split_are_two_fixed_points() {
    let provider = provider();
    let spanning = "NM_TEST.1:c.9_17delinsA";
    let split = "NM_TEST.1:c.[9_12del;14_17del]";

    assert_eq!(denotes(&provider, spanning), denotes(&provider, split));

    // `:47` recommends this one.
    assert_eq!(normalized(&provider, spanning), "NM_TEST.1:c.9_17delinsA");
    // `:17` and `general.md:34` demand this one — the members are two unchanged
    // nucleotides apart, and `general.md:35`'s codon exception needs a gap of
    // exactly one, so it cannot rescue either side here.
    assert_eq!(
        normalized(&provider, split),
        "NM_TEST.1:c.[9_12del;15_18del]",
        "pinned as an observation: whether `delins-merge-vs-individual-gap-two-or-more`'s \
         scope reaches a spelling that arrives already split is not settled"
    );
}
