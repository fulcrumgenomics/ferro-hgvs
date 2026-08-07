//! #1524 — the canonical split must never emit two members on *consecutive*
//! nucleotides.
//!
//! # The question
//!
//! Ferro split a lone `delins` into members with nothing unchanged between
//! them:
//!
//! ```text
//! NM_000038.6:c.5265_5268delinsTGCG -> c.[5265G>T;5266_5268delinsGCG]
//! NM_000083.3:c.2461_2464delinsCTCC -> c.[2461_2463delinsCTC;2464G>C]
//! ```
//!
//! In the first, member one ends at 5265 and member two begins at 5266; in the
//! second, at 2463 and 2464. May a description put two changed nucleotides that
//! touch into two members?
//!
//! # The ruling: no
//!
//! `assets/hgvs-nomenclature/docs/recommendations/DNA/delins.md:16`:
//!
//! > changes involving two or more consecutive nucleotides are described as
//! > deletion/insertion (delins) variants.
//!
//! Nothing competes with it. `general.md:34` ("two variants **separated by one
//! or more** nucleotides should be described individually and **not** as a
//! `delins`") governs members that something separates, so it does not reach
//! separation 0; and `delins.md:47`'s "the delins format is recommended"
//! escape hatch is about a *wider* merge, not about splitting a consecutive
//! run. So this is a plain rule with a single clause on point, not a
//! conflict needing a `rulings` record.
//!
//! **The issue and the campaign note cite this clause as `delins.md:15`. That
//! is off by one** — `:15` is "by definition, when **one** nucleotide is
//! replaced by **one** other nucleotide, the change is a substitution". The
//! clause quoted above is `:16`, which is also what
//! `merge::coalesce_adjacent_pieces` has cited all along.
//!
//! # Where it came from, and the second ruling it forced
//!
//! Both violations came from the codon-frame exception (`general.md:35`) in
//! `Normalizer::build_split_variants`, not from `merge::coalesce_adjacent_pieces`
//! (which implements this rule correctly and which a lone `c.` delins never
//! reaches — `is_splittable_single_member` admits only `g.`/`m.`). The branch
//! matched `[Sub@p; Identity@p+1; Sub@p+2]` with `p` and `p+2` in one codon and
//! emitted it as its own member, flushing whatever run was pending. That left
//! an adjacent member on the **left** when `p - 1` was changed, and on the
//! **right** when `p + 3` was.
//!
//! Fixing the right edge alone is wrong, and that was measured rather than
//! reasoned. Simply joining the triplet to the pending run merges across the
//! triplet's unchanged centre even when the left endpoint belongs to a longer
//! run — `c.9_13delinsACAAC` collapsing to one member whose two real variants
//! (`9_10delinsAC`, `12_13delinsAC`) span three codons — which cost **44
//! classes** of `cis_confluence_axis`'s 3' census (converged 6633 -> 6589) and
//! 44 of its 5' census.
//!
//! So the left edge is a second, separate ruling. `general.md:35` licenses the
//! merge for "two **variants** separated by one nucleotide, together affecting
//! one amino acid". When `p - 1` is changed, `p` is not a variant —
//! `delins.md:16` puts it inside the `delins` spanning its run — and that run
//! reaches beyond the codon the exception is about. The exception therefore
//! does not apply, which is exactly the precondition
//! `merge::apply_coding_codon_exception` already enforced on the other side of
//! the seam (`is_substitution(&pieces[index - 1])`). The two paths implement
//! one rule and only one of them was testing for it.
//!
//! # Consequence for representation stability
//!
//! This moves shipped output. Over the 5,761,894-row corpus audit that found
//! the defect (ClinVar 500k, multi-member cis, CMRG exhaustive, Paraphase
//! exhaustive), **58 of the 59 distinct violating outputs move**; of those, 12
//! become the input's own spelling and 46 become a differently-placed split.
//! The 59th (`NC_000017.11:g.[80110044C>G;80110045dup;80110047A>G]`) is not
//! this defect: its "adjacency" is an artefact of reading a `dup`'s named span
//! as its changed extent, and under the interbase reading the separation is 1.
//!
//! Note that the issue's framing — "every move is split -> merged, toward the
//! input form" — does not survive the second ruling. `c.5265_5268delinsTGCG`
//! does **not** come back as itself: its two variants are `5265_5266delinsTG`
//! and `5268T>G`, which affect two amino acids, so `general.md:34` splits them.
//! The input was merging across an unchanged nucleotide without the exception
//! applying, so it was not correct either.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Test transcript, shared with `issue_165_delins_sub_only_decompose` so the
/// two files argue about the same geometry.
///
/// ```text
///   c. axis: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 ...
///   base:    A T G C A A A A A  C  C  C  C  C  G  G  G  G  G  T ...
/// ```
///
/// CDS = `c.1..=c.60`, so `c.N`, `n.N` and `r.N` all name transcript position
/// `N` and the cross-axis rows below differ only in whether a reading frame
/// exists. Codons are `(base - 1) / 3 + 1`: codon 3 = `c.7..9`, codon 4 =
/// `c.10..12`, codon 5 = `c.13..15`.
const CORE: &str = "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT";

fn provider(accession: &str, coding: bool) -> MockProvider {
    let mut provider = MockProvider::new();
    let length = CORE.len() as u64;
    let (cds_start, cds_end) = if coding {
        (Some(1), Some(60))
    } else {
        (None, None)
    };
    provider.add_transcript(Transcript::new(
        accession.to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        CORE.to_string(),
        cds_start,
        cds_end,
        vec![Exon::new(1, 1, length)],
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

/// Normalize, and assert the result is its own normal form.
///
/// The fixed-point half is asserted here rather than left to
/// `FERRO_ASSERT_IDEMPOTENT`, because CI's gating `Test` job runs without the
/// oracles — an oracle-only guard would be green in the one job that blocks a
/// merge.
fn normalized_fixed_point(provider: &MockProvider, descriptor: &str) -> String {
    let variant = parse_hgvs(descriptor).expect("fixture must parse");
    let normalizer = Normalizer::new(provider.clone());
    let once = normalizer
        .normalize(&variant)
        .expect("fixture must normalize")
        .to_string();
    let twice = normalizer
        .normalize(&parse_hgvs(&once).expect("output must re-parse"))
        .expect("output must normalize")
        .to_string();
    assert_eq!(
        twice, once,
        "`{descriptor}` normalized to a form that is not a fixed point",
    );
    once
}

/// Every member boundary in `description`, as `(end of one, start of the next)`
/// on the HGVS axis, for the adjacency assertions below.
///
/// Deliberately crude — it reads the string rather than the parsed variant —
/// because what `delins.md:16` constrains is the *description*, and a test that
/// re-derived the spans from the same code that produced them would be checking
/// that code against itself.
fn member_bounds(description: &str) -> Vec<(u64, u64)> {
    let body = description
        .split_once(":c.")
        .or_else(|| description.split_once(":n."))
        .or_else(|| description.split_once(":r."))
        .expect("fixture is a c./n./r. description")
        .1;
    let Some(inner) = body.strip_prefix('[').and_then(|s| s.strip_suffix(']')) else {
        return Vec::new();
    };
    let spans: Vec<(u64, u64)> = inner
        .split(';')
        .map(|member| {
            let digits = |s: &str| -> u64 {
                s.chars()
                    .take_while(char::is_ascii_digit)
                    .collect::<String>()
                    .parse()
                    .expect("member starts with a plain coordinate")
            };
            let start = digits(member);
            let end = match member.split_once('_') {
                Some((_, rest)) => digits(rest),
                None => start,
            };
            (start, end)
        })
        .collect();
    spans.windows(2).map(|w| (w[0].1, w[1].0)).collect()
}

/// No two members of `description` may sit on consecutive nucleotides.
fn assert_no_adjacent_members(description: &str) {
    for (previous_end, next_start) in member_bounds(description) {
        assert!(
            next_start > previous_end + 1,
            "`{description}` puts members on consecutive nucleotides \
             ({previous_end} then {next_start}) — delins.md:16",
        );
    }
}

// ---------------------------------------------------------------------------
// Adjudicated-correct: the two shapes the corpus audit found.
// ---------------------------------------------------------------------------

/// Right-edge adjacency — the `NM_000083.3:c.2461_2464delinsCTCC` shape.
///
/// `c.10..13` is `CCCC`; the alt makes 10, 12 and 13 mismatch and leaves 11
/// unchanged. 10 and 12 are both in codon 4, so `general.md:35` merges them —
/// and 13 touches 12, so `delins.md:16` puts it in the same member. The
/// output is the input's own spelling.
#[test]
fn a_run_touching_the_triplet_on_the_right_is_one_member() {
    let provider = provider("NM_TEST.1", true);
    let output = normalized_fixed_point(&provider, "NM_TEST.1:c.10_13delinsTCAG");
    assert_no_adjacent_members(&output);
    assert_eq!(output, "NM_TEST.1:c.10_13delinsTCAG");
}

/// Left-edge adjacency — the `NM_000038.6:c.5265_5268delinsTGCG` shape.
///
/// `c.9..12` is `ACCC`; the alt makes 9, 10 and 12 mismatch and leaves 11
/// unchanged. 10 and 12 share codon 4, but 9 touches 10, so the left endpoint
/// of the would-be triplet is not a variant of its own (`delins.md:16`) and the
/// pair the exception is offered — `9_10delins` and `12C>A` — affects codons 3,
/// 4 and 4, i.e. two amino acids. `general.md:34` therefore governs and they
/// are described individually.
///
/// This is the ruling that costs the most downstream: the corpus rows of this
/// shape do **not** come back as their inputs.
#[test]
fn a_run_touching_the_triplet_on_the_left_declines_the_codon_exception() {
    let provider = provider("NM_TEST.1", true);
    let output = normalized_fixed_point(&provider, "NM_TEST.1:c.9_12delinsTACA");
    assert_no_adjacent_members(&output);
    assert_eq!(output, "NM_TEST.1:c.[9_10delinsTA;12C>A]");
}

/// Two codon triplets back to back — the `NM_000392.5:c.4465_4473delinsGGCCCACAG`
/// shape, where the two merged members touched at 4470/4471.
#[test]
fn back_to_back_codon_triplets_do_not_touch() {
    let provider = provider("NM_TEST.1", true);
    let output = normalized_fixed_point(&provider, "NM_TEST.1:c.10_15delinsTCAGCT");
    assert_no_adjacent_members(&output);
    assert_eq!(output, "NM_TEST.1:c.[10_13delinsTCAG;15G>T]");
}

// ---------------------------------------------------------------------------
// The discriminating cases. Each fails if the fix is over-generalised into a
// blanket merge.
// ---------------------------------------------------------------------------

/// Separation 1 across a codon boundary still splits (`general.md:34`).
///
/// `c.12` is in codon 4 and `c.14` in codon 5, so the pair does not "together
/// affect one amino acid" and the exception does not reach it. A blanket
/// "merge anything close" rule would fail here.
#[test]
fn a_pair_separated_by_one_across_a_codon_boundary_still_splits() {
    let provider = provider("NM_TEST.1", true);
    let output = normalized_fixed_point(&provider, "NM_TEST.1:c.12_14delinsACA");
    assert_no_adjacent_members(&output);
    assert_eq!(output, "NM_TEST.1:c.[12C>A;14C>A]");
}

/// Separation 2 still splits, in one codon or not.
///
/// The exception is for a gap of *exactly* one; `c.10` and `c.13` are two apart
/// and `general.md:34`'s plain rule applies.
#[test]
fn a_pair_separated_by_two_still_splits() {
    let provider = provider("NM_TEST.1", true);
    let output = normalized_fixed_point(&provider, "NM_TEST.1:c.10_13delinsTCCA");
    assert_no_adjacent_members(&output);
    assert_eq!(output, "NM_TEST.1:c.[10C>T;13C>A]");
}

/// A run **separated** from the triplet is still its own member.
///
/// This is the boundary the left-edge ruling turns on: the exception is
/// declined only when the pending run *touches* `p`. Here `c.9` is unchanged,
/// so the run ends two before the triplet, `general.md:34` keeps it apart, and
/// the triplet at `c.10..12` still merges under `general.md:35`.
#[test]
fn a_run_separated_from_the_triplet_still_splits() {
    let provider = provider("NM_TEST.1", true);
    let output = normalized_fixed_point(&provider, "NM_TEST.1:c.8_12delinsGATCA");
    assert_no_adjacent_members(&output);
    assert_eq!(output, "NM_TEST.1:c.[8A>G;10_12delinsTCA]");
}

/// The codon exception is untouched where its left endpoint *is* a lone
/// variant — the `general.md:35` case the spec's own `c.9002_9009delinsTTT`
/// note is about.
#[test]
fn a_lone_codon_pair_still_merges() {
    let provider = provider("NM_TEST.1", true);
    assert_eq!(
        normalized_fixed_point(&provider, "NM_TEST.1:c.[10C>T;12C>A]"),
        "NM_TEST.1:c.10_12delinsTCA",
    );
}

/// An axis with no reading frame has no codon exception, so it has no branch to
/// fix — and must not acquire one.
///
/// `n.` on the same coding transcript and on a non-coding one both split at the
/// unchanged base per `general.md:34`, exactly as they did before.
#[test]
fn an_axis_without_a_reading_frame_is_unchanged() {
    let coding = provider("NM_TEST.1", true);
    let noncoding = provider("NR_TEST.1", false);
    let output = normalized_fixed_point(&coding, "NM_TEST.1:n.10_13delinsTCAG");
    assert_no_adjacent_members(&output);
    assert_eq!(output, "NM_TEST.1:n.[10C>T;12_13delinsAG]");
    assert_eq!(
        normalized_fixed_point(&noncoding, "NR_TEST.1:n.10_13delinsTCAG"),
        "NR_TEST.1:n.[10C>T;12_13delinsAG]",
    );
}

/// The `r.` axis is codon-frame-aware on a coding transcript (#275), so it takes
/// the same branch and must reach the same answer as `c.`.
#[test]
fn the_rna_axis_reaches_the_same_form() {
    let provider = provider("NM_TEST.1", true);
    let output = normalized_fixed_point(&provider, "NM_TEST.1:r.10_13delinsucag");
    assert_no_adjacent_members(&output);
    assert_eq!(output, "NM_TEST.1:r.10_13delinsucag");
}
