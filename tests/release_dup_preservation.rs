//! Featureless regression guard for the ruled dup-preservation consult.
//!
//! # What it pins
//!
//! At separation zero a tandem duplication abutting a substitution
//! (`g.[4_5dup;6C>A]`) must be KEPT, not folded into a spanning `delins` — the
//! decided ruling `separation-zero-dup-member-is-preserved-not-merged`
//! (`duplication.md:18` governs). `collapse_overlapping_cis_edits` in
//! `src/normalize/merge.rs` consults the ruled partitioner to preserve it.
//!
//! # Why this test exists SEPARATELY from the dev-gated one
//!
//! `merge.rs::the_default_preserves_a_dup_abutting_a_change` pins the same
//! behaviour, but it is `#[cfg(feature = "dev")]`-gated — and the bug it now
//! guards was ITSELF a `#[cfg(feature = "dev")]` gate on the consult. With
//! `DEFAULT_PARTITION_RULE = Ruled`, that gate made the shipped RELEASE binary
//! skip the consult and fold the dup (8 of 500,004 ClinVar rows; 62 of the
//! 9.95M-row corpus), while the dev build preserved it. A test that only runs
//! under `--features dev` shares the bug's own predicate, so it could never
//! observe the divergence: both the test and the bug were on the same side of
//! the gate.
//!
//! This test is DELIBERATELY not feature-gated and uses only the public API
//! (`MockProvider`, `Normalizer`, `parse_hgvs` — none dev-gated), so it can be
//! run in a build with the default feature set. Run featureless it fails the
//! moment the consult is put back behind any feature gate the shipped binary
//! does not carry. CI runs it that way in the `build-lint` job; locally:
//!
//! ```sh
//! cargo test --test release_dup_preservation          # default features (the guard)
//! cargo test --features dev --test release_dup_preservation   # also passes
//! ```
//!
//! Running it WITH `--features dev` still passes but proves nothing about the
//! release binary — the point is the featureless run.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// `NC_TEST:g.[4_5dup;6C>A]` over `GGGTACGGG` — `4_5dup`'s source `TA` has no
/// adjacent copy, so the pre-fix member-merge collapse folded it to
/// `g.6delinsTAA`. The ruled consult keeps the dup.
#[test]
fn the_release_default_preserves_a_dup_abutting_a_change() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_TEST", "GGGTACGGG".to_string());
    let normalizer = Normalizer::new(provider);
    let exposed = parse_hgvs("NC_TEST:g.[4_5dup;6C>A]").expect("parse");

    assert_eq!(
        normalizer
            .normalize(&exposed)
            .expect("normalize")
            .to_string(),
        "NC_TEST:g.[4_5dup;6C>A]",
        "the shipped default must preserve a dup abutting a change \
         (separation-zero-dup-member-is-preserved-not-merged). A fold to a \
         spanning delins here means the ruled consult in \
         collapse_overlapping_cis_edits is skipped in this build — most likely \
         re-gated behind a feature the release binary does not carry.",
    );
}
