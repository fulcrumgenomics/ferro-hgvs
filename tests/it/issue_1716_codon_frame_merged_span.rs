//! Issue #1716: the codon-frame merge must test the span it authorises.
//!
//! `general.md:34` — "two variants separated by one or more nucleotides should
//! be described individually and **not** as a 'delins'" (`DNA/delins.md:17`
//! restates it verbatim) — is the governing clause at separation one. The only
//! thing that lifts it is `general.md:35`'s exception, and that clause has
//! **two** conjuncts:
//!
//! > **exception**: two variants separated by one nucleotide, together
//! > affecting one amino acid, should be described as a "delins".
//!
//! `merge_consecutive_edits`' codon-frame predicate implemented the second
//! conjunct as `same_codon(prev_a.end, next_start)` — the left member's *right*
//! edge — while the `delins` it authorises spans `prev_a.start ..= next.end`.
//! For a one-position left member the two are the same test, which is why the
//! term was harmless when #104 wrote it under `prev_a.start == prev_a.end`.
//! #292 relaxed that to `prev_a.start <= prev_a.end` to admit chain extension,
//! and from then on a **wider** left member could satisfy the endpoint test
//! while the merged span demonstrably covered four or more reference positions
//! — which cannot lie in one codon, so the exception's second conjunct fails
//! and nothing licenses the merge.
//!
//! # The fixture is the real reproducer, hermetically
//!
//! The transcript below is the first 210 nucleotides of the **real** DMD CDS
//! (`NM_004006.2` / `LRG_199t1`, the transcript the spec itself uses for
//! `c.145_147delinsTGG`), with the 5'UTR trimmed so `cds_start = 1` and every
//! `c.` coordinate here is the coordinate on the real transcript. Each string
//! below was first observed against the prepared reference through
//! `ferro normalize --reference <dir>`; running them through this fixture
//! reproduces those outputs byte for byte, which is what makes a `MockProvider`
//! legitimate here — the guard needs no `FERRO_MANIFEST` and therefore cannot
//! skip green.
//!
//! # What each case is for
//!
//! Two cases are the defect, and three are the surrounding population that a
//! narrowed predicate must **not** disturb: the spec's own worked example, the
//! chain extension #292 was written for, and a frameshift pair inside one
//! codon. The ledger record `codon-carve-out-shape-restriction` is `decided`
//! (WIDEN — the exception reaches every edit type), and it separately records
//! that it does not implement the frameshift half of its own reasoning, leaving
//! that question unanswered. This change does not answer it either.
//!
//! **All three of those are end-to-end pins on the emitted string, and none of
//! them guards this predicate.** The PR that added them said so for the chain
//! extension and not for the other two; it is true of all three, and measuring
//! it is what found the gap. Forcing `codon_frame_eligible` to `false` — the
//! codon-frame arm disabled outright — leaves every one of them **passing**,
//! because the sequence-first canonicalizer re-derives those strings from the
//! resulting sequence rather than from the members `merge_consecutive_edits` was
//! handed. They are worth having as composite-output pins; they are not evidence
//! about the arm.
//!
//! **The guards on the arm's two call-site conjuncts are therefore not in this
//! file at all.** `prev_a.start <= prev_a.end` (which keeps an insertion anchor
//! out) and the 1-position `next_anchor` width check are both invisible from
//! here for the same measured reason — the canonicalizer re-derives the same
//! string either way — so they are pinned as unit tests beside the predicate,
//! in `normalize::merge`'s own `mod tests`:
//! `an_insertion_anchor_is_refused_by_the_codon_frame_arm`,
//! `a_multi_position_next_anchor_is_refused_by_the_codon_frame_arm`, and the
//! positive control `the_codon_frame_arm_still_merges_the_shape_it_is_for`
//! that stops the other two passing on a dead arm. Each states its own
//! measurement.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// First 210 nt of the DMD coding sequence (`NM_004006.2` c.1..c.210).
///
/// Bases that matter below, 1-indexed on the `c.` axis: `c.18=A`, `c.19=G`,
/// `c.20=T`, `c.21=A`; `c.22=G`, `c.23=A`, `c.24=G`; `c.30=T`, `c.31=T`,
/// `c.32=A`, `c.33=T`; `c.145=C`, `c.146=G`, `c.147=C`.
const DMD_CDS_PREFIX: &str = concat!(
    "ATGCTTTGGTGGGAAGAAGTAGAGGACTGTTATGAAAGAGAAGATGTTCAAAAGAAAACA",
    "TTCACAAAATGGGTAAATGCACAATTTTCTAAGTTTGGGAAGCAGCATATTGAGAACCTC",
    "TTCAGTGACCTACAGGATGGGAGGCGCCTCCTAGACCTCCTCGAAGGCCTGACAGGGCAA",
    "AAACTGCCAAAAGAAAAAGGATCCACAAGA",
);

/// The DMD CDS prefix as a single-exon coding transcript with `cds_start = 1`,
/// so `c.N` indexes byte `N` of [`DMD_CDS_PREFIX`].
fn provider_dmd_cds() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = DMD_CDS_PREFIX.to_string();
    let len = sequence.len() as u64;
    let exons = vec![Exon::new(1, 1, len)];
    let transcript = Transcript::new(
        "NM_004006.2".to_string(),
        Some("DMD".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(len),
        exons,
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

fn normalize(input: &str) -> String {
    let normalizer = Normalizer::new(provider_dmd_cds());
    let variant = parse_hgvs(input).expect("parse failed");
    format!(
        "{}",
        normalizer.normalize(&variant).expect("normalize failed")
    )
}

// ---------------------------------------------------------------------------
// The defect: a net-length-changing pair whose merged span crosses a codon
// ---------------------------------------------------------------------------

/// `c.[18_19del;21A>C]` merged into `c.18_21delinsTC` on the shipping arm.
///
/// The merged span is `c.18..c.21`: `c.18` is in **codon 6** (positions 16-18)
/// and `c.19_21` is **codon 7**. Two codons, so the pair cannot "together
/// affect one amino acid" and `general.md:35` does not reach it under any
/// reading — the merge is unlicensed whatever force `general.md:34` carries.
/// With the exception out, `general.md:34` / `DNA/delins.md:17` govern and the
/// two variants are described individually.
#[test]
fn del_then_sub_across_a_codon_boundary_stays_individual() {
    assert_eq!(
        normalize("NM_004006.2:c.[18_19del;21A>C]"),
        "NM_004006.2:c.[18_19del;21A>C]",
    );
}

/// Second real reproducer, same shape one codon along: the merged span
/// `c.30..c.33` covers **codon 10** (positions 28-30) and **codon 11**
/// (31-33). Merged into `c.30_33delinsAG` on the shipping arm.
#[test]
fn del_then_sub_across_a_codon_boundary_stays_individual_second_row() {
    assert_eq!(
        normalize("NM_004006.2:c.[30_31del;33T>G]"),
        "NM_004006.2:c.[30_31del;33T>G]",
    );
}

// ---------------------------------------------------------------------------
// The population that must not move
// ---------------------------------------------------------------------------

/// The spec's own worked example (`DNA/delins.md:37` — `LRG_199t1:c.145_147delinsTGG`),
/// on the transcript the spec quotes it against. Both changed positions and the
/// unchanged one lie in **codon 49** (`c.145_147`), so the exception applies and
/// it is the *split* that `DNA/delins.md:42` marks `class="invalid"`. Narrowing
/// the predicate to the merged span must leave this merging.
///
/// **This is an end-to-end pin on the emitted string, not a guard on the
/// predicate** — the same disclosure [`equal_length_chain_extension_is_unchanged`]
/// carries, and for the same reason. Measured by mutation: with the codon-frame
/// arm forced off entirely, this test still **passes**, because the
/// sequence-first canonicalizer re-derives `c.145_147delinsTGG` from the
/// resulting sequence without ever consulting `merge_consecutive_edits`. What it
/// pins is that the composite output does not move; what it cannot detect is a
/// change to the arm this PR narrows. The two guards that *are* sensitive to it
/// are the two reproducers above (they fail when the pre-#1716 endpoint test is
/// restored) and `merge_consecutive_edits_tests::test_codon_frame_then_strict_chain`.
#[test]
fn spec_worked_example_still_merges_within_one_codon() {
    assert_eq!(
        normalize("NM_004006.2:c.[145C>T;147C>G]"),
        "NM_004006.2:c.145_147delinsTGG",
    );
}

/// Chain extension, which is what #292 relaxed `prev_a.start == prev_a.end` to
/// `prev_a.start <= prev_a.end` in order to admit. `c.18` and `c.19` merge by
/// strict adjacency first, so the left anchor is already two positions wide
/// when `c.21` is offered — the exact configuration the narrowed predicate
/// refuses. This output must survive the narrowing, and it does, because the
/// sequence-first canonicalizer re-derives the same grouping from the resulting
/// sequence rather than from the members it was handed.
#[test]
fn equal_length_chain_extension_is_unchanged() {
    assert_eq!(
        normalize("NM_004006.2:c.[18A>C;19G>T;21A>C]"),
        "NM_004006.2:c.[18_19inv;21A>C]",
    );
}

/// A net-length-changing pair whose merged span *does* lie in one codon
/// (`c.22_24` is **codon 8**). Whether `general.md:35`'s second conjunct can be
/// satisfied by a pair that frameshifts every downstream residue is a question
/// the ledger's `codon-carve-out-shape-restriction` states and then declines to
/// implement — see the module doc. This fix answers only the span question, so
/// the row must not move.
///
/// Same disclosure as [`spec_worked_example_still_merges_within_one_codon`]: an
/// end-to-end pin on the string, measured to survive the codon-frame arm being
/// forced off, so it is not a guard on the predicate either.
#[test]
fn frameshift_pair_inside_one_codon_is_left_open() {
    assert_eq!(
        normalize("NM_004006.2:c.[22del;24G>C]"),
        "NM_004006.2:c.22_24delinsAC",
    );
}
