//! Issue #1135: merging a self-cancelling del+ins cis allele must not panic.
//!
//! `build_naedit` chose the insertion branch on the anchor *shape*
//! (`start == end + 1`) without checking whether anything was left to insert.
//! When a deletion 3'-shifts onto an adjacent insertion of the same base — the
//! two cancel — `merged.alt` is empty while the anchor keeps insertion shape,
//! so an `NaEdit::Insertion` was built with an **empty** sequence. That is not
//! a valid HGVS edit, and `normalize_na_edit` then divided by zero rotating the
//! (empty) inserted sequence through the repeat:
//! `attempt to calculate the remainder with a divisor of zero`.
//!
//! These run without a prepared reference so they gate every PR, not just the
//! nightly manifest job. Expected strings are pinned exactly rather than
//! pattern-matched: the identity's span comes from the merged anchor
//! arithmetic, so a change there has to fail loudly instead of sliding past a
//! suffix check.

use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// A contig carrying an `AAAAA` homopolymer at `g.1004..1008`, long enough for
/// a deletion to 3'-shift across it onto an insertion anchor.
fn homopolymer_normalizer() -> Normalizer<MockProvider> {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(
        "NC_000001.11",
        format!("{}ATGAAAAAGCCTAA{}", "N".repeat(1000), "N".repeat(100)),
    );
    Normalizer::new(provider)
}

/// The reported input and two neighbouring spellings: a deletion and an
/// insertion of the same base inside one homopolymer cancel out, leaving the
/// sequence unchanged across the merged anchor.
#[test]
fn self_cancelling_del_ins_allele_normalizes_to_an_identity() {
    let normalizer = homopolymer_normalizer();
    for (input, expected) in [
        (
            "NC_000001.11:g.[1004del;1008_1009insA]",
            "NC_000001.11:g.1009_1010=",
        ),
        (
            "NC_000001.11:g.[1004del;1004_1005insA]",
            "NC_000001.11:g.1005_1006=",
        ),
        (
            "NC_000001.11:g.[1005del;1008_1009insA]",
            "NC_000001.11:g.1009_1010=",
        ),
    ] {
        let variant = parse_hgvs(input).unwrap();
        let normalized = normalizer
            .normalize(&variant)
            .unwrap_or_else(|e| panic!("{input} failed to normalize: {e}"));
        let rendered = normalized.to_string();

        assert_eq!(rendered, expected, "{input} normalized unexpectedly");

        // Whatever we emit has to be real HGVS, and stable under a second pass.
        parse_hgvs(&rendered)
            .unwrap_or_else(|e| panic!("{input} -> {rendered} does not re-parse: {e}"));
        assert_eq!(
            normalizer.normalize(&normalized).unwrap().to_string(),
            rendered,
            "{input} is not idempotent"
        );
    }
}

/// Guard against the fix over-triggering: a del+ins pair that does **not**
/// cancel must still merge into a real edit, not collapse to an identity.
#[test]
fn non_cancelling_del_ins_allele_still_merges() {
    let normalizer = homopolymer_normalizer();
    let input = "NC_000001.11:g.[1004del;1008_1009insGG]";
    let variant = parse_hgvs(input).unwrap();
    let rendered = normalizer.normalize(&variant).unwrap().to_string();

    // RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
    // 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`),
    // back to the two-member form this row held before #1296.
    //
    //   was  NC_000001.11:g.1008delinsGG           one member, re-derived
    //   now  NC_000001.11:g.[1008del;1009G[3]]     two members, as authored
    //
    // The input asserts a deletion inside the `A` run and an insertion of `GG`
    // at the junction `1008|1009`, which is what `general.md:34` describes
    // individually. The single-member `delins` was reached by re-deriving the
    // partition from the resulting sequence — the note this replaces says so in
    // as many words ("the sequence-first pass no longer refuses to read the
    // promoted repeat, and what the pair denotes is one reference base replaced
    // by two") — and that is the move the ruling removes.
    //
    // The sequence is unchanged. Measured with `spec_corpus::denotation_of`,
    // whose applier goes through `hgvs_to_spdi` rather than the normalizer: both
    // the input and this output denote `…G AAAA GGG CC…` across 1003-1011, and
    // the output is a fixed point. The guard this test exists for is unchanged
    // in substance and in strength — the answer is still a real edit and pinned
    // exactly, not the identity the #1135 fix could wrongly collapse it to, and
    // the identity cases above still pin `g.…=` for the pairs that do cancel.
    assert_eq!(
        rendered, "NC_000001.11:g.[1008del;1009G[3]]",
        "{input} changes the sequence and must not collapse to an identity"
    );
    parse_hgvs(&rendered)
        .unwrap_or_else(|e| panic!("{input} -> {rendered} does not re-parse: {e}"));
}
