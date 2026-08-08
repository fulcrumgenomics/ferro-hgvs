//! The del+dup collision repair reaches the non-coding transcript axis.
//!
//! `respell_colliding_duplications` (#1259) resolves the #1263 collision: a
//! duplication whose junction lands inside a sibling's span, which ferro's own
//! parser rejects as a self-cancelling allele. It was gated to
//! `CisKind::Genome | CisKind::Mt`, so on the transcript axes `normalize` still
//! emitted a description `parse_hgvs` refuses:
//!
//! ```text
//! reference   ATGGCAACAGCGTAAACGGCTAGCTAGCTA        (the g. motif at n.4)
//! in    n.[7_8del;9_10insA;11delinsAC]
//! out   n.[8_9del;9dup;10_11insA]
//!       Parse error: Self-cancelling allele: variants at index 0 (del) and
//!       1 (dup) describe overlapping reference positions
//! ```
//!
//! #1284 records why lifting that gate looked hard: the repair reads the
//! duplicated bases straight out of the provider over the member's own
//! coordinates, and the `c.` axis does not index the transcript sequence —
//! `c.-1` is `cds_start - 1`, not `cds_start - 2`, so a naive `delta` is off by
//! one below the CDS start, and wiring it up that way was tried and reverted.
//!
//! That blocker is specific to the **CDS-relative** axes. `n.` is
//! transcript-relative: `axis_frame` already returns `delta: 0` for
//! `CisKind::Tx`, exactly as it does for the genomic axes, so the member's
//! coordinates index the fetched sequence directly and the repair needs no
//! conversion at all. Admitting `Tx` therefore closes the whole non-coding
//! transcript axis without touching the discontinuity.
//!
//! `c.` and `r.` are still refused, and #1284 stays open for them.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// A transcript carrying the `GCAACAGCGTAAAC` motif of the genomic #1263
/// reproducer at position 4, so the `n.` and `g.` cases describe the same bases.
const TX: &str = "ATGGCAACAGCGTAAACGGCTAGCTAGCTA";

fn normalize(provider: &MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider.clone());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize `{input}`: {e}"))
        .to_string()
}

/// Assert `input` normalizes to `expected`, that the result is one the parser
/// accepts, and that it is a fixed point.
///
/// The re-parse is the load-bearing half: the defect is not a wrong spelling
/// but an **unparseable** one, so a test that only compared strings would still
/// pass against the broken output if someone pinned it.
fn assert_repairs_to(provider: &MockProvider, input: &str, expected: &str) {
    let once = normalize(provider, input);
    assert_eq!(once, expected, "normalizing `{input}`");
    let reparsed = parse_hgvs(&once)
        .unwrap_or_else(|e| panic!("`{input}` -> `{once}`, which ferro cannot re-parse: {e}"));
    let twice = Normalizer::new(provider.clone())
        .normalize(&reparsed)
        .unwrap_or_else(|e| panic!("re-normalizing `{once}`: {e}"))
        .to_string();
    assert_eq!(twice, once, "`{input}` -> `{once}` is not a fixed point");
}

#[test]
fn a_noncoding_del_dup_collision_is_repaired() {
    // The #1284 reproduction on the axis that needs no coordinate conversion.
    let provider = SyntheticBuilder::noncoding(TX, Strand::Plus).build();
    assert_repairs_to(
        &provider,
        "NR_TEST.1:n.[7_8del;9_10insA;11delinsAC]",
        "NR_TEST.1:n.[8del;10_11insA]",
    );
}

#[test]
fn the_repair_preserves_the_sequence_the_allele_denotes() {
    // The repair rewrites the members, so the binding check is not the spelling
    // but the bases: applied to the transcript, input and output must agree.
    // `apply_with` converts each member through `hgvs_to_spdi` — a different
    // code path from the repair — and declines an overlapping description, so
    // the pre-fix output has no resulting sequence at all rather than a wrong
    // one.
    let provider = SyntheticBuilder::noncoding(TX, Strand::Plus).build();
    for input in [
        "NR_TEST.1:n.[7_8del;9_10insA;11delinsAC]",
        "NR_TEST.1:n.[7_8del;11delinsAC]",
    ] {
        let output = normalize(&provider, input);
        let from_input = crate::common::cis_apply_oracle::apply_with(&provider, TX, input)
            .unwrap_or_else(|| panic!("`{input}` does not apply"));
        let from_output = crate::common::cis_apply_oracle::apply_with(&provider, TX, &output)
            .unwrap_or_else(|| {
                panic!("`{input}` -> `{output}`, which has no well-defined resulting sequence")
            });
        assert_eq!(
            from_output, from_input,
            "`{input}` -> `{output}` changed the sequence"
        );
    }
}

/// The shorter route to the same collision, without the middle insertion.
///
/// **Re-blessed from `n.[8C>G;10del]` under the partition model.** The old form
/// was *derived*: the two sequences differ over `n.8_10` (`CAG` -> `GA`), and
/// minimising edit distance over that block reads it position-wise as a
/// substitution at 8 and a deletion at 10, separated by the unchanged `A` at 9.
/// Nothing in the input said that. The input asserts a **two-member partition** —
/// a deletion at `n.7_8` and a delins at `n.11` — separated by the two unchanged
/// nucleotides `n.9`/`n.10`, so `general.md:34` ("two variants separated by one or
/// more nucleotides should be described individually and **not** as a "delins"")
/// licenses no merge and both members survive.
///
/// Each is then placed and re-spelled, which is where the new string comes from
/// and neither step is new:
///
/// - `n.7_8del` deletes `AC`; `n.7` and `n.9` are both `A`, so the 3' rule
///   (`general.md:41`, "the most 3' position possible of the reference sequence is
///   arbitrarily assigned to have been changed") moves it one base to
///   **`n.8_9del`**, which deletes `CA` for the same resulting sequence.
/// - `n.11delinsAC` replaces `C` with `AC`; the payload's `C` is the reference
///   base itself, so trimming the shared suffix respells the member as the
///   insertion it is — **`n.10_11insA`** — and the `A` cannot shift 3' because
///   `n.11` is `C`.
///
/// The two then sit one unchanged nucleotide apart (`n.10`), which `general.md:34`
/// again keeps individual — the `n.` axis is non-coding, so `general.md:35`'s
/// one-amino-acid exception cannot reach it.
///
/// **What this row no longer exercises is the del+dup collision itself, and nor
/// does its sibling above.** The #1263 shape came out of a *derived* `dup` — the
/// module header's `out n.[8_9del;9dup;10_11insA]` — and no partition is derived
/// any more, so both `n.` rows now reach a disjoint, re-parseable output without
/// `respell_colliding_duplications` having to fire. What they still assert is the
/// invariant #1284 was filed for and the one `assert_repairs_to` checks: that the
/// `n.` axis emits a description ferro's own parser accepts, that it is a fixed
/// point, and (in [`the_repair_preserves_the_sequence_the_allele_denotes`]) that
/// it denotes the bases the input did. The repair's own coverage is the genomic
/// #1263 suite, not this file.
#[test]
fn the_simple_two_member_collision_is_repaired_too() {
    let provider = SyntheticBuilder::noncoding(TX, Strand::Plus).build();
    assert_repairs_to(
        &provider,
        "NR_TEST.1:n.[7_8del;11delinsAC]",
        "NR_TEST.1:n.[8_9del;10_11insA]",
    );
}
