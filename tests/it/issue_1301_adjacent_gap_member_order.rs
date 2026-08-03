//! Two members that share a span but add sequence at different junctions must
//! render in junction order.
//!
//! `cis_member_order_key` broke a start tie with the member's **end** (#1261).
//! That is not enough for two members whose spans are identical:
//!
//! ```text
//! reference  ("ACGT" x 64) + "GCATGAAAAT" + ("ACGT" x 64)
//! g.[263_264insAC;264_265insAA]  ->  g.[264_265dup;264_265insCA]
//! ```
//!
//! `264_265insCA` fills the gap `264|265`, so it adds at interbase 264;
//! `264_265dup` places its copy after its last base, so it adds at 265. Both
//! span `264_265`, so start and end tie and the order fell to the descriptor
//! tie-break — which sorts `264_265dup` first, the reverse of the order the two
//! apply in. The pair then rendered out of order and their interbase spans read
//! as overlapping, which is #1235's second acceptance criterion.
//!
//! The sequence was never wrong here: both members applied denote the input's
//! bases. What was wrong is the rendering. `junction_rank` adds the missing
//! discriminator — a member adding at its span's *start* sorts before one acting
//! across the whole span.

use crate::common::synthetic::{padded, SyntheticBuilder};
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, Normalizer};

const CORE: &str = "GCATGAAAAT";

fn normalize(input: &str) -> String {
    let provider = SyntheticBuilder::genomic(CORE).build();
    Normalizer::new(provider)
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize")
        .to_string()
}

/// The interbase start of each member, in rendered order, via `hgvs_to_spdi` —
/// a different code path from the ordering key under test.
fn member_starts(descriptor: &str) -> Vec<u64> {
    let provider = SyntheticBuilder::genomic(CORE).build();
    let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).expect("parse") {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    members
        .iter()
        .map(|member| {
            hgvs_to_spdi(member, &provider)
                .expect("member converts to SPDI")
                .position
        })
        .collect()
}

/// Assert the output denotes the input's bases, so an ordering test cannot pass
/// by rendering a *wrong* pair tidily.
fn assert_preserving(input: &str, output: &str) {
    let provider = SyntheticBuilder::genomic(CORE).build();
    let reference = padded(CORE);
    let denote = |descriptor: &str| -> Option<String> {
        let members: Vec<HgvsVariant> = match parse_hgvs(descriptor).ok()? {
            HgvsVariant::Allele(allele) => allele.variants.clone(),
            single => vec![single],
        };
        let mut triples = Vec::new();
        for member in &members {
            triples.push(hgvs_to_spdi(member, &provider).ok()?);
        }
        triples.sort_by_key(|t| std::cmp::Reverse(t.position));
        let mut edited = reference.as_bytes().to_vec();
        let mut claimed = reference.len();
        for triple in &triples {
            let start = usize::try_from(triple.position).ok()?;
            let end = start.checked_add(triple.deletion.len())?;
            if end > claimed {
                return None;
            }
            edited.splice(start..end, triple.insertion.bytes());
            claimed = start;
        }
        String::from_utf8(edited).ok()
    };
    assert_eq!(
        denote(output).expect("output applies"),
        denote(input).expect("input applies"),
        "`{input}` -> `{output}` changed the sequence"
    );
}

#[test]
fn an_insertion_sorts_before_a_duplication_sharing_its_span() {
    let input = "NC_TEST.1:g.[263_264insAC;264_265insAA]";
    let output = normalize(input);
    // Since #1235 the pair reaches one member: the merged form is derived from
    // the sequence rather than assembled per member, so the tie this file is
    // about never has to be broken on *this* input. It is still broken, and
    // still asserted, on the protein cases below — where both members tie on
    // start and end — and the interbase-order check below is unchanged and
    // trivially satisfied by a single member.
    assert_eq!(output, "NC_TEST.1:g.264_265insCAAA");
    assert_preserving(input, &output);

    let starts = member_starts(&output);
    assert!(
        starts.windows(2).all(|pair| pair[0] <= pair[1]),
        "`{output}` renders its members out of interbase order: {starts:?}"
    );
}

#[test]
fn the_order_does_not_depend_on_how_it_was_authored() {
    assert_eq!(
        normalize("NC_TEST.1:g.[263_264insAC;264_265insAA]"),
        normalize("NC_TEST.1:g.[264_265insAA;263_264insAC]"),
    );
}

// ---------------------------------------------------------------------------
// The protein axis reaches the same tie (#1328).
//
// These three assert on the rendered descriptor alone, where the nucleotide
// cases above also cross-check against SPDI (`assert_preserving`,
// `member_starts`). That oracle is unavailable here rather than skipped:
// `hgvs_to_spdi` rejects `HgvsVariant::Protein` outright — "protein variants
// cannot be represented in SPDI" (`src/spdi/convert.rs`) — so there is no
// second code path on this axis to derive a member's junction from. What is
// under test is also precisely the discriminator positions cannot supply: both
// members tie on start *and* end, which is what drops the order onto
// `junction_rank`.
// ---------------------------------------------------------------------------

/// Normalize a protein descriptor. Ordering is positional, so no reference
/// sequence is needed — the same approach
/// `issue_1261_cis_member_order::protein_members_sharing_a_start_render_in_end_order`
/// takes.
fn normalize_protein(input: &str) -> String {
    use ferro_hgvs::{MockProvider, Normalizer};
    Normalizer::new(MockProvider::new())
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize")
        .to_string()
}

/// The protein axis hits this tie for exactly the reasons the nucleotide axes
/// do, and was excluded only because #1301 did not want to move protein output.
///
/// Both preconditions hold: `cis_member_range` returns a real range for
/// `HgvsVariant::Protein`, so two protein members sharing a span tie on
/// `(start, end)`; and `sort_cis_members_by_genomic_order` is axis-agnostic.
/// With `junction_rank` matching only the six nucleotide-bearing arms, every
/// protein member ranked `1` and the tie fell to the rendered-descriptor
/// tie-break — which sorts `dup` before `ins` alphabetically, reversing the
/// order the two apply in.
///
/// `Gly4_Ala5insSer` fills the gap between residues 4 and 5, so it adds at the
/// start of its span; `Gly4_Ala5dup` copies both residues and places the copy
/// after residue 5, so it adds at the end. The insertion must render first.
#[test]
fn a_protein_insertion_sorts_before_a_duplication_sharing_its_span() {
    assert_eq!(
        normalize_protein("NP_000001.1:p.[Gly4_Ala5insSer;Gly4_Ala5dup]"),
        "NP_000001.1:p.[Gly4_Ala5insSer;Gly4_Ala5dup]",
    );
}

/// The same pair authored the other way round must render identically — the
/// property that makes the order a *total* one rather than input-dependent.
///
/// The *property* already held before the fix — the descriptor tie-break is
/// itself total, so both authorings agreed on `[dup;ins]`. The **test** did
/// not: it pins the post-fix rendering, so it fails on the pre-fix code like
/// its sibling above (verified by removing the protein arm and re-running).
/// What it adds is that the two failure modes cannot be traded for one another
/// — a rank that left protein order depending on how the allele was written
/// would be a worse bug than the one being fixed, and agreeing-but-wrong is
/// caught by the sibling test rather than by this one.
#[test]
fn the_protein_order_does_not_depend_on_how_it_was_authored() {
    let expected = "NP_000001.1:p.[Gly4_Ala5insSer;Gly4_Ala5dup]";
    for input in [
        "NP_000001.1:p.[Gly4_Ala5insSer;Gly4_Ala5dup]",
        "NP_000001.1:p.[Gly4_Ala5dup;Gly4_Ala5insSer]",
    ] {
        assert_eq!(normalize_protein(input), expected, "authored as `{input}`");
    }
}

/// Protein members that are *not* insertions still rank with the 3' end, so a
/// pair of them sharing a span keeps the descriptor tie-break that #1261
/// established. This is the arm's negative control: it passes both before and
/// after the fix, which is the point — the new rank must move `ins` and nothing
/// else.
///
/// The pair has to share **both** endpoints to reach the contract. `Ala3dup`
/// and `Ala3_Ala4dup` do not: their ends differ, so `cis_member_order_key`
/// settles them on range alone and neither `junction_rank` nor the descriptor
/// tie-break is consulted. `Ala3_Ala4del` and `Ala3_Ala4dup` tie on
/// `(start, end)`, so the rank runs, both members must come out `1`, and the
/// descriptors break it — `del` before `dup`, alphabetically. Authoring them
/// the other way round makes the assertion sensitive rather than a restatement
/// of the input: a rank of `0` on either member would put `dup` first and fail
/// here.
///
/// Both authorings are checked for the same reason the insertion pair's are —
/// the order has to be a property of the members, not of how they were written.
///
/// Scope note: despite what this comment used to claim, nothing here exercises
/// the **nucleotide** ranking. That parity is guarded by the nucleotide tests
/// earlier in this file, which the protein arm cannot reach — it returns before
/// the `NaEdit` match rather than altering it.
#[test]
fn protein_members_that_are_not_insertions_keep_their_existing_order() {
    let expected = "NP_000001.1:p.[Ala3_Ala4del;Ala3_Ala4dup]";
    for input in [
        "NP_000001.1:p.[Ala3_Ala4dup;Ala3_Ala4del]",
        "NP_000001.1:p.[Ala3_Ala4del;Ala3_Ala4dup]",
    ] {
        assert_eq!(normalize_protein(input), expected, "authored as `{input}`");
    }
}
