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
    assert_eq!(output, "NC_TEST.1:g.[264_265insCA;264_265dup]");
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
