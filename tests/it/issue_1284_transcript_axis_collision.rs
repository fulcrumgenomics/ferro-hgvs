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
        // The collision is repaired to one spanning `delins`. The block is
        // `CAG` -> `AGA`, equal length 3/3:
        //
        //     position-wise:            C->A, A->G, G->A            cost 3
        //     del C, A=A, G=G, ins A                                cost 2  <- minimal
        //
        // `unchanged-is-read-over-every-minimal-alignment` originally read `9`
        // and `10` as unchanged (matched in every cost-2 alignment) and kept the
        // members individual. #2174 scopes that record: this is an **equal-length
        // run with no interior zero-shift fixed point** — no base survives at its
        // own coordinate, the position-wise C/A, A/G, G/A all differ — so the
        // minimal-alignment interior matches are shift artifacts, not separators.
        // The run is one contiguous change and `delins.md:16` types it
        // `n.8_10delinsAGA`. The record now carries the superseding scope.
        //
        // Same block and same reasoning as
        // `normalize_reparse_invariant::a_del_beside_a_dup_re_spells_instead_of_colliding`,
        // reached from an unrelated input; see that row for the full write-up.
        "NR_TEST.1:n.8_10delinsAGA",
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

#[test]
fn the_simple_two_member_collision_is_repaired_too() {
    // The shorter route to the same collision, without the middle insertion.
    let provider = SyntheticBuilder::noncoding(TX, Strand::Plus).build();
    assert_repairs_to(
        &provider,
        "NR_TEST.1:n.[7_8del;11delinsAC]",
        "NR_TEST.1:n.[8C>G;10del]",
    );
}
