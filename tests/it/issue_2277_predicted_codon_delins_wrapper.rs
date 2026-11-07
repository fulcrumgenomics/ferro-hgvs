//! #2277: normalization keeps a predicted `(...)` wrapper when it rebuilds a
//! predicted change as a single member.
//!
//! Two paths rebuilt one member certain:
//!
//! - the canonical split on the coding axis, where a codon triplet
//!   `[Sub; Identity; Sub]` re-merges into one 3-base `delins`
//!   (`general.md:34`), so `c.(4_6delinsAGA)` came back `c.4_6delinsAGA`; and
//! - a predicted cis allele `[(a;b)]` whose members merge to one, which kept a
//!   one-member allele that renders bare, so `g.[(1003C>A;1004G>T)]` came back
//!   `g.1003_1004delinsAT` on every axis.
//!
//! The `(...)` marks a predicted change (`recommendations/general.md:12`, and
//! `:85` for the parentheses); dropping it states an observed one.
//!
//! One 9-base transcript, `ATGCGCTAA`, on `NC_000001.11:1000-1008` with the
//! whole transcript as CDS, so `c.N` is `n.N` is `g.(999+N)`. Codon 2 is
//! `CGC` at `c.4_6`.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::Normalizer;

fn normalizer() -> Normalizer<MockProvider> {
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_TX1.1".to_string(),
        Some("GENE1".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    provider.add_genomic_sequence(
        "NC_000001.11",
        format!("{}ATGCGCTAA{}", "N".repeat(999), "N".repeat(100)),
    );
    Normalizer::new(provider)
}

fn normalize(input: &str) -> String {
    let v = ferro_hgvs::parse_hgvs(input).unwrap_or_else(|e| panic!("{input} parses: {e}"));
    normalizer()
        .normalize(&v)
        .unwrap_or_else(|e| panic!("{input} normalizes: {e}"))
        .to_string()
}

/// Every pinned output is its own fixed point, so the wrapper also survives a
/// second pass.
fn assert_normalizes_to(input: &str, expected: &str) {
    assert_eq!(normalize(input), expected, "normalize({input})");
    assert_eq!(normalize(expected), expected, "normalize({expected})");
}

/// The issue's reproducer: the codon exception keeps `4_6` as one `delins`,
/// and that `delins` stays predicted.
#[test]
fn a_predicted_codon_delins_keeps_its_wrapper_on_the_coding_axis() {
    assert_normalizes_to("NM_TX1.1:c.(4_6delinsAGA)", "NM_TX1.1:c.(4_6delinsAGA)");
    // A different payload through the same triplet: `CGC` -> `AGT`.
    assert_normalizes_to("NM_TX1.1:c.(4_6delinsAGT)", "NM_TX1.1:c.(4_6delinsAGT)");
}

/// `r.` on a coding transcript takes the same codon-frame branch and dropped
/// the wrapper the same way.
#[test]
fn a_predicted_codon_delins_keeps_its_wrapper_on_the_rna_axis() {
    assert_normalizes_to("NM_TX1.1:r.(4_6delinsaga)", "NM_TX1.1:r.(4_6delinsaga)");
}

/// The observed control: a certain input does not gain a wrapper, in either
/// spelling.
#[test]
fn an_observed_codon_delins_stays_observed() {
    assert_normalizes_to("NM_TX1.1:c.4_6delinsAGA", "NM_TX1.1:c.4_6delinsAGA");
    assert_normalizes_to("NM_TX1.1:c.[4C>A;6C>A]", "NM_TX1.1:c.4_6delinsAGA");
}

/// Unchanged on purpose: with no codon frame the change splits into two
/// members, and the wrapper goes round the allele.
#[test]
fn the_noncoding_and_genomic_spellings_are_unchanged() {
    assert_normalizes_to("NM_TX1.1:n.(4_6delinsAGA)", "NM_TX1.1:n.[(4C>A;6C>A)]");
    assert_normalizes_to(
        "NC_000001.11:g.(1003_1005delinsAGA)",
        "NC_000001.11:g.[(1003C>A;1005C>A)]",
    );
}

/// The allele spelling of the same predicted codon change settles on the
/// same string as the single-edit spelling.
#[test]
fn a_predicted_allele_merged_by_the_codon_exception_stays_predicted() {
    assert_normalizes_to("NM_TX1.1:c.[(4C>A;6C>A)]", "NM_TX1.1:c.(4_6delinsAGA)");
}

/// The allele half is not codon-specific: a predicted allele whose adjacent
/// members merge to one `delins` lost the wrapper on every axis.
#[test]
fn a_predicted_allele_merged_to_one_member_stays_predicted_on_every_axis() {
    assert_normalizes_to(
        "NC_000001.11:g.[(1003C>A;1004G>T)]",
        "NC_000001.11:g.(1003_1004delinsAT)",
    );
    assert_normalizes_to("NM_TX1.1:n.[(4C>A;5G>T)]", "NM_TX1.1:n.(4_5delinsAT)");
    assert_normalizes_to("NM_TX1.1:c.[(4C>A;5G>T)]", "NM_TX1.1:c.(4_5delinsAT)");
}
