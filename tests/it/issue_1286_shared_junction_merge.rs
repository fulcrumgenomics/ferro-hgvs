//! Two cis members that settle on one junction must merge, not claim it twice.
//!
//! A member that consumes no base occupies the interbase junction where its
//! payload lands. Two of them can reach the same junction from distinct starting
//! points — two one-base insertions at adjacent gaps inside a homopolymer both
//! travel to the tract's 3'-most junction:
//!
//! ```text
//! reference  ("ACGT" x 64) + "AAAAAA" + ("ACGT" x 64)   poly-A run at 257-263
//! g.[258_259insA;259_260insA]  ->  g.[263dup;263dup]
//! ```
//!
//! `parse_hgvs` **accepts** that, so `FERRO_ASSERT_REPARSE` does not catch it,
//! while the SPDI apply oracle correctly declines it — two members claiming one
//! position have no well-defined resulting sequence.
//!
//! Nothing upstream prevented it: the sibling clamps bound a member against a
//! sibling's *bases*, and neither of these has any. And nothing downstream could
//! repair it: the sequence-first canonicalization would derive the correct single
//! member, but it runs after the per-member pipeline and derives from the
//! allele's resulting sequence — which an overlapping allele does not have. It
//! declines, so the corruption survives the one pass that would have merged it.
//! Traced: the pass is entered with `[263dup;263dup]` and exits without
//! deriving.
//!
//! `coalesce_members_at_one_junction` merges them before that pass runs. Both add
//! their payload at the same interbase point, so one insertion carrying the
//! concatenation denotes exactly what the pair denoted; the next pass
//! re-canonicalises it to a `dup` or a repeat where the reference permits.

use crate::common::cis_apply_oracle::{apply, assert_normalizes_preserving, normalize};
use crate::common::synthetic::{padded, SyntheticBuilder};
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, Normalizer};

/// Normalize against the padded synthetic contig and assert the result denotes
/// the same bases the input did — through `hgvs_to_spdi`, a different code path
/// from the normalizer.
fn assert_padded_preserving(core: &str, input: &str) -> String {
    let provider = SyntheticBuilder::genomic(core).build();
    let normalizer = Normalizer::new(provider.clone());
    let reference = padded(core);
    let variant = parse_hgvs(input).expect("parse");
    let output = normalizer
        .normalize(&variant)
        .expect("normalize")
        .to_string();

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
                return None; // overlapping members
            }
            edited.splice(start..end, triple.insertion.bytes());
            claimed = start;
        }
        String::from_utf8(edited).ok()
    };

    let from_input = denote(input).expect("input applies");
    let from_output = denote(&output)
        .unwrap_or_else(|| panic!("`{input}` -> `{output}` has no resulting sequence"));
    assert_eq!(
        from_output, from_input,
        "`{input}` -> `{output}` changed the sequence"
    );
    output
}

#[test]
fn two_insertions_reaching_one_junction_merge() {
    // #1286. Seven A at 257-263 plus two inserted A is a nine-A run, which
    // canonicalises to the repeat spelling.
    let output = assert_padded_preserving("AAAAAA", "NC_TEST.1:g.[258_259insA;259_260insA]");
    assert_eq!(output, "NC_TEST.1:g.257_263A[9]");
}

#[test]
fn the_merge_is_independent_of_the_authored_order() {
    let forward = assert_padded_preserving("AAAAAA", "NC_TEST.1:g.[258_259insA;259_260insA]");
    let reverse = assert_padded_preserving("AAAAAA", "NC_TEST.1:g.[259_260insA;258_259insA]");
    assert_eq!(forward, reverse);
}

#[test]
fn three_insertions_reaching_one_junction_merge() {
    // The group is not limited to pairs.
    let output = assert_padded_preserving(
        "AAAAAA",
        "NC_TEST.1:g.[258_259insA;259_260insA;260_261insA]",
    );
    assert_eq!(output, "NC_TEST.1:g.257_263A[10]");
}

#[test]
fn insertions_at_distinct_junctions_are_left_alone() {
    // The pass must fire on a *shared* junction only. These two land in
    // different tracts and stay two members.
    let seq = "ACGTTTTTTACGTAAAAAAACGT";
    let output = normalize(seq, "TEMPLATE:g.[5_6insT;15_16insA]");
    assert!(
        output.starts_with("TEMPLATE:g.["),
        "two disjoint insertions must stay two members, got `{output}`"
    );
    assert_eq!(
        apply(seq, "TEMPLATE:g.[5_6insT;15_16insA]"),
        apply(seq, &output),
        "sequence changed"
    );
}

#[test]
fn a_single_insertion_is_untouched() {
    // No group, nothing to merge — the ordinary path must be unaffected.
    assert_normalizes_preserving("ACGTTTTTTACGT", "TEMPLATE:g.5_6insT", "TEMPLATE:g.9dup");
}
