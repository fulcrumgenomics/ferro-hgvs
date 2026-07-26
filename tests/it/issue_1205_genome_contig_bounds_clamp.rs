//! Issue #1205 — a `g.`-axis insertion that saturates at a contig bound must be
//! clamped to a `delins`, exactly as the transcript axes already do (#1202).
//!
//! `normalize_genome` had no boundary clamp. #1202 introduced the rewrite and
//! wired it into `normalize_cds`, `normalize_tx` and `normalize_rna` — three call
//! sites, all numbered against a transcript — leaving the most heavily used axis
//! uncovered. Both ends escape, in the two shapes #1202 describes:
//!
//! | direction | input | before | why it is wrong |
//! |---|---|---|---|
//! | 5' | `g.1_2insGC` | `g.1insCG` | single-position insertion; ferro's own parser rejects it |
//! | 3' | `g.13_14insGGA` | `g.15_16insAGG` | `g.16` is past the end of a 15-base contig |
//!
//! The 5' shape is the one the spec names, on this very axis
//! (`DNA/insertion.md:95-101`):
//!
//! > **Can I describe a variant as `g.123insG`?**
//! > No, since the description is not unequivocal, it is not allowed.
//! > What does the description mean, the insertion of a `G` **at** position
//! > `g.123` or the insertion of a `G` **after** position `g.123`?
//!
//! An insertion needs two adjacent flanking anchors; `g.1ins…` has one.
//!
//! **Reachability.** GRCh38 primary chromosomes open with telomeric `N` runs, so
//! this is unlikely on a primary assembly. It is reachable on small or custom
//! references whose first bases are a homopolymer — plasmids, viral genomes,
//! scaffolds, synthetic contigs — and on any `NC_`/`NG_` slice a user supplies.

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// The issue's probe contig: `C`×4, `A`×7, `G`×4. 15 bases, so both ends are
/// within one fetch window and both bounds are exercisable.
const CONTIG: &str = "CCCCAAAAAAAGGGG";

fn normalize(input: &str, dir: ShuffleDirection) -> String {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_PROBE.1", CONTIG);
    let normalizer =
        Normalizer::with_config(provider, NormalizeConfig::default().with_direction(dir));
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
        .to_string()
}

/// Assert the exact canonical form, that it re-parses, and that it is a fixed
/// point.
///
/// Re-parsing is not boilerplate here: the 5' defect's output was rejected by
/// ferro's own parser, so "can ferro read back what it just wrote" is the
/// property under test as much as the string itself.
fn assert_canonical(input: &str, dir: ShuffleDirection, expected: &str) {
    let once = normalize(input, dir);
    assert_eq!(once, expected, "{input} ({dir:?}) canonical form");
    parse_hgvs(&once).unwrap_or_else(|e| panic!("{once} must re-parse: {e}"));
    assert_eq!(normalize(&once, dir), once, "{once} must be a fixed point");
}

/// Apply "replace 1-based inclusive `start..=end` of [`CONTIG`] with `insert`",
/// returning the edited sequence. `start > end` denotes a pure insertion between
/// the two flanking positions.
///
/// Deliberately dumb string splicing, independent of the normalizer: it is what
/// lets the equivalence assertions below prove the rewrite preserves the
/// haplotype rather than merely restating what the normalizer did.
fn splice(start: usize, end: usize, insert: &str) -> String {
    let (before, rest) = CONTIG.split_at(start - 1);
    let after = &rest[(end + 1).saturating_sub(start)..];
    format!("{before}{insert}{after}")
}

// ---------------------------------------------------------------------------
// The defect at each contig bound
// ---------------------------------------------------------------------------

/// 5' saturation. The shuffle drives the insertion to rest immediately 5' of
/// `g.1`, where the coordinate identity gives `delete ref[1]` (= `C`),
/// `insert CG ++ C`.
#[test]
fn five_prime_saturated_insertion_clamps_at_contig_start() {
    assert_canonical(
        "NC_PROBE.1:g.1_2insGC",
        ShuffleDirection::FivePrime,
        "NC_PROBE.1:g.1delinsCGC",
    );
}

/// 3' saturation, the mirror: `delete ref[15]` (= `G`), `insert G ++ AGG`.
///
/// The alt is `GGA` rather than a clean homopolymer because a tandem alt never
/// reaches this shape on `g.` — with no codon-frame gate on the genomic axis, two
/// or more copies become repeat notation and a single copy becomes a `dup`, both
/// before any saturation. `GGA` shifts twice through the `G`×4 tract and then
/// mismatches, which is what leaves it an insertion resting on the bound.
#[test]
fn three_prime_saturated_insertion_clamps_at_contig_end() {
    assert_canonical(
        "NC_PROBE.1:g.13_14insGGA",
        ShuffleDirection::ThreePrime,
        "NC_PROBE.1:g.15delinsGAGG",
    );
}

/// The rewrite must describe the same haplotype, not merely a parseable one.
///
/// Both sides are spliced into the reference by hand here, so this does not go
/// through the normalizer and cannot be satisfied by an internally-consistent but
/// wrong rewrite.
#[test]
fn the_clamped_delins_describes_the_same_sequence() {
    // 5': insert `GC` between g.1 and g.2  ==  replace g.1 with `CGC`.
    assert_eq!(splice(2, 1, "GC"), splice(1, 1, "CGC"));
    // 3': insert `GGA` between g.13 and g.14  ==  replace g.15 with `GAGG`.
    //
    // Note the input is stated at g.13_14 but the shuffled payload is `AGG`; the
    // 3'-shifted insertion between g.15 and g.16 is the same haplotype, which the
    // first equality below pins before the delins comparison.
    assert_eq!(splice(14, 13, "GGA"), splice(16, 15, "AGG"));
    assert_eq!(splice(16, 15, "AGG"), splice(15, 15, "GAGG"));
}

// ---------------------------------------------------------------------------
// Controls: the clamp must not fire away from the bounds
// ---------------------------------------------------------------------------

/// The same alt one position 5' of the saturating case comes to rest inside the
/// contig, so it keeps two valid anchors and must stay an insertion.
#[test]
fn interior_insertion_is_left_alone() {
    assert_canonical(
        "NC_PROBE.1:g.12_13insGGA",
        ShuffleDirection::ThreePrime,
        "NC_PROBE.1:g.14_15insAGG",
    );
}

/// …and on the 5' side.
#[test]
fn five_prime_interior_insertion_is_left_alone() {
    assert_canonical(
        "NC_PROBE.1:g.3_4insGC",
        ShuffleDirection::FivePrime,
        "NC_PROBE.1:g.2_3insCG",
    );
}

/// A tandem expansion whose tract runs to the contig start becomes repeat
/// notation, which is not an insertion and so is not the clamp's business. Guards
/// against the clamp firing on position alone rather than on shape.
#[test]
fn repeat_notation_reaching_the_contig_start_is_untouched() {
    assert_canonical(
        "NC_PROBE.1:g.2_3insCC",
        ShuffleDirection::FivePrime,
        "NC_PROBE.1:g.1_4C[6]",
    );
}

/// An insertion that cannot shuffle at all is unaffected, even written against
/// the first bases of the contig.
#[test]
fn non_shuffling_insertion_at_the_contig_start_is_untouched() {
    assert_canonical(
        "NC_PROBE.1:g.1_2insCG",
        ShuffleDirection::FivePrime,
        "NC_PROBE.1:g.1_2insCG",
    );
}

// ---------------------------------------------------------------------------
// Scope boundary: circular references are NOT clamped
// ---------------------------------------------------------------------------

/// `m.` (and `o.`) are deliberately left out of this fix, and this pins that
/// boundary so a later change does not extend the clamp there by reflex.
///
/// A circular reference wraps: position 1 HAS a valid 5' neighbour — the last
/// base — so the correct output for a 5'-saturated insertion there is a
/// wraparound description, not this `delins` rewrite. Applying the genomic clamp
/// would produce a well-formed but semantically wrong answer, which is worse than
/// the status quo.
///
/// To be explicit about what that status quo is: `m.` currently emits the same
/// invalid single-position form this PR fixes on `g.` (`m.1insCG`). That is a real
/// defect and is **not** fixed here — it needs the circular answer, tracked
/// separately. This test therefore asserts only that the clamp did not fire; it
/// does not bless the output.
#[test]
fn mito_axis_is_not_clamped() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_012920.1", CONTIG);
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NC_012920.1:m.1_2insGC").expect("parse");
    let rendered = normalizer
        .normalize(&variant)
        .expect("normalize")
        .to_string();
    assert!(
        !rendered.contains("delins"),
        "the genomic contig clamp must not reach the circular axis, got {rendered}",
    );
}
