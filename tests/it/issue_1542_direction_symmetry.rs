//! The shuffle direction may not change how many members a variant normalizes
//! to — asserted **end to end, through `Normalizer`** (#1542).
//!
//! # Why this exists beside the unit guard
//!
//! `merge.rs`'s `direction_symmetry` module states the same property over
//! *every* `PartitionRule`, which is the durable half: a future flip of
//! `DEFAULT_PARTITION_RULE` cannot neuter it. It cannot be the whole of the
//! answer, though, because it lives in the same file as the code it guards —
//! reverting `src/normalize/merge.rs` wholesale deletes the guard along with the
//! fix, and a deletion is not a failure.
//!
//! This module is the half that survives that. It drives the **public**
//! normalization API on whatever arm the process was started with (the shipped
//! default, in CI), so a revert of `merge.rs` reddens it rather than removing
//! it. Measured against `origin/main` at `2d8b490b`, the corpus below carries
//! **36** direction-dependent rows on the shipped default; the two guards
//! therefore fail for different reasons and neither can be satisfied by
//! deleting the other.
//!
//! # What the defect looked like
//!
//! ```text
//! reference   ACGTACGTACGTACGTACGT AGCA ACGTACGTACGTACGTACGT
//! input       NC_TEST.1:g.21_24delinsCATGC
//!
//! 3'  ->  g.[20_21insC;21_22insT;25del]     three members
//! 5'  ->  g.[20_21insC;22_24inv]            two members
//! ```
//!
//! Both directions cut the same block into the same three pieces and still hold
//! three after the shift; it is a *later* merging pass that sees an inversion at
//! one placement and not at the other. Each answer is a fixed point and both
//! denote the same bases, so no idempotency, re-parse, in-bounds or
//! denoted-sequence oracle can see it — which is why the property has to be
//! asserted directly.
//!
//! # Authority
//!
//! `canonical-form-choice-when-both-legal` (decided 2026-08-07): a form is
//! derived from the **resulting sequence**, not preserved from the input's
//! spelling. A partition that differs by shuffle direction is derived from the
//! traversal that produced it — a fact about how the normalizer got there, not
//! about the sequence it arrived at. `separation-is-a-property-of-the-spelling-
//! not-of-the-variant` says the same thing of the separation `general.md:34`
//! keys on: it is read off the partition re-derived from the resulting sequence.

use rayon::prelude::*;

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

/// The flanking bases every core is embedded in.
///
/// Long enough that a shift in either direction stays inside the served
/// sequence, and periodic so a pure indel at the block edge has somewhere to
/// roll — a core in unique flanks would pin every piece and leave the corpus
/// unable to express the defect at all.
const PAD: &str = "ACGTACGTACGTACGTACGT";

/// The core length. With [`PAYLOAD_LEN`] this is the shape that reaches the
/// re-derivation: one spanning `delins` whose payload differs in length from
/// its span arrives as a single member and is re-partitioned from the resulting
/// sequence.
const CORE_LEN: usize = 4;

/// The payload length. Differs from [`CORE_LEN`] on purpose — an equal-length
/// payload is a pure substitution block with no indel to place, so the shift
/// has nothing to move and the corpus could not vary the property under test.
const PAYLOAD_LEN: usize = 5;

/// Every word of length `n` over the DNA alphabet, in a fixed order.
fn words(n: usize) -> Vec<String> {
    let mut out = vec![String::new()];
    for _ in 0..n {
        out = out
            .iter()
            .flat_map(|word| {
                ["A", "C", "G", "T"]
                    .into_iter()
                    .map(move |base| format!("{word}{base}"))
            })
            .collect();
    }
    out
}

/// Normalize `input` against `provider` in `direction`, through the public API.
fn normalize_in(
    provider: &MockProvider,
    input: &str,
    direction: ShuffleDirection,
) -> Option<String> {
    let variant = parse_hgvs(input).ok()?;
    let config = NormalizeConfig::for_entry_point(direction, ErrorConfig::lenient());
    Normalizer::with_config(provider.clone(), config)
        .normalize(&variant)
        .ok()
        .map(|normalized| normalized.to_string())
}

/// How many members a normalized description renders as.
///
/// Read off the re-parsed output rather than off any internal structure, so the
/// count this guard compares is the one a consumer would see.
fn member_count(description: &str) -> Option<usize> {
    Some(match parse_hgvs(description).ok()? {
        HgvsVariant::Allele(allele) => allele.variants.len(),
        _ => 1,
    })
}

/// One row's verdict: how many comparisons it contributed, whether it reached a
/// multi-member re-derivation, and any disagreement.
struct Row {
    compared: usize,
    declined: usize,
    /// Whether the 3' answer has more than one member. Summed into the
    /// non-vacuity floor below: a corpus on which every row is a single member
    /// cannot vary the quantity this module compares, so its zero would be
    /// structural rather than a result.
    multi_member: bool,
    disagreement: Option<String>,
}

/// The member count of a normalized variant does not depend on the shuffle
/// direction (#1542).
///
/// Only the **count** is compared. The placement legitimately differs — that is
/// what the direction is *for* — so asserting string equality would pin 3' onto
/// 5' and forbid shuffling altogether.
///
/// Every row must answer in both directions: a row that errors on one and not
/// the other is a disagreement of the same kind, so declines are counted and
/// refused rather than skipped.
#[test]
fn the_member_count_does_not_depend_on_the_shuffle_direction() {
    let lo = PAD.len() + 1;
    let hi = PAD.len() + CORE_LEN;
    let payloads = words(PAYLOAD_LEN);

    let rows: Vec<Row> = words(CORE_LEN)
        .par_iter()
        .flat_map_iter(|core| {
            let reference = format!("{PAD}{core}{PAD}");
            let mut provider = MockProvider::new();
            provider.add_genomic_sequence("NC_TEST.1", reference);
            payloads.iter().map(move |payload| {
                let input = format!("NC_TEST.1:g.{lo}_{hi}delins{payload}");
                let three = normalize_in(&provider, &input, ShuffleDirection::ThreePrime);
                let five = normalize_in(&provider, &input, ShuffleDirection::FivePrime);
                let (Some(three), Some(five)) = (three, five) else {
                    return Row {
                        compared: 0,
                        declined: 1,
                        multi_member: false,
                        disagreement: None,
                    };
                };
                let (Some(three_count), Some(five_count)) =
                    (member_count(&three), member_count(&five))
                else {
                    return Row {
                        compared: 0,
                        declined: 1,
                        multi_member: false,
                        disagreement: None,
                    };
                };
                Row {
                    compared: 1,
                    declined: 0,
                    multi_member: three_count > 1,
                    disagreement: (three_count != five_count).then(|| {
                        format!(
                            "  core={core} {input}\n    3': {three} ({three_count} member(s))\n    \
                             5': {five} ({five_count} member(s))"
                        )
                    }),
                }
            })
        })
        .collect();

    let compared: usize = rows.iter().map(|row| row.compared).sum();
    let declined: usize = rows.iter().map(|row| row.declined).sum();
    let multi_member = rows.iter().filter(|row| row.multi_member).count();
    let disagreements: Vec<&String> = rows
        .iter()
        .filter_map(|row| row.disagreement.as_ref())
        .collect();

    // A zero below is only a claim about the rows that were compared, so the
    // denominator is asserted first. `declined` is held at zero rather than
    // tolerated: a change that made this corpus stop normalizing would empty the
    // comparison set while leaving the assertion below trivially satisfiable.
    assert_eq!(
        declined, 0,
        "{declined} row(s) failed to normalize in both directions, so they say nothing either way"
    );
    assert_eq!(
        compared,
        4usize.pow(CORE_LEN as u32) * 4usize.pow(PAYLOAD_LEN as u32),
        "the corpus did not enumerate every core x payload pair"
    );
    // The other half of "a zero is a claim about the corpus": a corpus that
    // never re-derives more than one member could not disagree about the count
    // however broken the partitioner was, so its zero would be structural. This
    // asserts the generator can vary the very thing the guard keys on.
    assert!(
        multi_member > 1_000,
        "only {multi_member} of {compared} row(s) re-derive into more than one member, so the \
         corpus cannot vary the member count compared below and a zero there would say nothing"
    );

    assert!(
        disagreements.is_empty(),
        "{} of {compared} row(s) normalize to a different number of members depending only on the \
         shuffle direction. Direction moves a member's placement; it may not change how many \
         members there are, or the partition is being read off a placement rather than off the \
         sequence — which `canonical-form-choice-when-both-legal` (derive from the resulting \
         sequence) and `separation-is-a-property-of-the-spelling-not-of-the-variant` both rule \
         against. Showing the first 20:\n{}",
        disagreements.len(),
        disagreements
            .iter()
            .take(20)
            .map(|row| row.as_str())
            .collect::<Vec<_>>()
            .join("\n")
    );
}
