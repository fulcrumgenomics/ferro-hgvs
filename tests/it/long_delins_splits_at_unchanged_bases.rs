//! A delins must split at its unchanged bases regardless of how long it is.
//!
//! The canonicalizer originally refused to partition any block longer than
//! 32 nt, on the reading that `delins.md:44-47` keeps a long replacement as one
//! spanning delins. That conflated two different things. The spec's concern
//! there is *coincidental* alignment — "parts of the inserted sequence align" —
//! not length as such, and the separation rule (`delins.md:17`) carries no
//! length qualifier at all. The coincidence concern is now handled
//! structurally: only single-gap alignments are considered, so a block is split
//! only where a base genuinely survives unchanged.
//!
//! The remaining cap is a cost bound on a quadratic alignment, far above the
//! sizes at which this question is interesting.

use crate::common::synthetic::{hgvs, SyntheticBuilder};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    format!(
        "{}",
        normalizer
            .normalize(&variant)
            .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"))
    )
}

/// A 40 nt delins — comfortably past the old 32 nt cap — whose interior leaves
/// one base untouched. It must be described as two members, exactly as the same
/// shape is at 5 nt.
#[test]
fn a_delins_longer_than_the_old_cap_still_splits() {
    // Core: 40 bases. The replacement changes everything except position 20.
    let core = "ACGT".repeat(10);
    let mut replacement: Vec<char> = core
        .chars()
        .map(|base| match base {
            'A' => 'C',
            'C' => 'A',
            'G' => 'T',
            _ => 'G',
        })
        .collect();
    // Keep index 19 (HGVS position 20 within the core) identical to the
    // reference, so the block has a genuine unchanged interior base.
    replacement[19] = core.as_bytes()[19] as char;
    let replacement: String = replacement.into_iter().collect();

    let provider = SyntheticBuilder::genomic(&core).build();
    let input = hgvs(
        &format!("NC_TEST.1:g.{{0}}_{{1}}delins{replacement}"),
        &[1, core.len() as u64],
    );
    let out = normalize(provider, &input);
    assert!(
        out.contains(';'),
        "a 40 nt delins with an unchanged interior base must split into members, got `{out}`"
    );
    assert!(
        !out.contains(&replacement),
        "the spanning replacement must not survive intact, got `{out}`"
    );
}
