//! Commuting is tested against the payload a member would carry **at** the
//! junction it lands on, not the one it carries now.
//!
//! `clamp_sibling_crossing_junctions` lets a member meet a sibling whose payload
//! commutes with its own, because sharing a junction is well-defined there and
//! the pair then merges (#1286, #1308). Landing rotates the payload into phase —
//! and a rotation is exactly what can destroy the commuting identity:
//!
//! ```text
//! reference  ("ACGT" x 64) + "TAAAACCA" + ("ACGT" x 64)
//! g.[260_261insAC;261_262insAC]
//!   `AC` and `AC` commute, so the first member was allowed onto junction 261
//!   at 261 its payload is `CA`, and `CA ++ AC != AC ++ CA`
//!
//!   intended  T A A A A C A A C C C A
//!   emitted   T A A A A A C C A C C A
//! ```
//!
//! Two members then share a junction with no defined order between them, and the
//! merge downstream concatenates them in *rendered* order — producing `ACCA`,
//! which is neither payload nor any rotation of the pair.
//!
//! Testing the landed payload keeps them apart: the bound falls back to one
//! position short, and the allele renders as its input does.

use crate::common::synthetic::{padded, SyntheticBuilder};
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, Normalizer};

/// Normalize and assert the output denotes the input's bases, projected through
/// `hgvs_to_spdi` rather than the normalizer's own applier.
fn assert_padded_preserving(core: &str, input: &str) -> String {
    let provider = SyntheticBuilder::genomic(core).build();
    let normalizer = Normalizer::new(provider.clone());
    let reference = padded(core);
    let output = normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
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
                return None;
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
fn a_rotation_that_breaks_commuting_keeps_the_pair_apart() {
    // #1312. `AC` at 260 becomes `CA` at 261, which does not commute with the
    // `AC` already there, so the member may not land on that junction.
    let output = assert_padded_preserving("TAAAACCA", "NC_TEST.1:g.[260_261insAC;261_262insAC]");
    assert_eq!(output, "NC_TEST.1:g.[260_261insAC;261_262insAC]");
}

#[test]
fn a_rotation_that_preserves_commuting_still_lets_them_merge() {
    // The guard. In a poly-A tract every rotation of `A` is `A`, so the landed
    // payload still commutes and the pair merges into one member.
    let output = assert_padded_preserving("AAAAAA", "NC_TEST.1:g.[258_259insA;259_260insA]");
    assert_eq!(output, "NC_TEST.1:g.257_263A[9]");
}
