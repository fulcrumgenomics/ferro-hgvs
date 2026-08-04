//! Issue #1217 — an `m.`-axis insertion that saturates at a mitochondrial
//! terminus must be clamped to a `delins`, exactly as the transcript axes
//! (#1202) and the `g.` axis (#1205) already do.
//!
//! `normalize_mt` was the one caller of `clamp_insertion_at_sequence_bounds`
//! never wired up, so both ends escaped in the two shapes #1202 describes:
//!
//! | direction | input | before | why it is wrong |
//! |---|---|---|---|
//! | 5' | `m.1_2insGC` | `m.1insCG` | single-position insertion; ferro's own parser rejects it |
//! | 3' | `m.13_14insGGA` | `m.15_16insAGG` | `m.16` is past the end of a 15-base contig |
//!
//! The 5' shape is the one the spec names by example
//! (`DNA/insertion.md:95-101`): an insertion needs two adjacent flanking
//! anchors, and `m.1ins…` has one.
//!
//! ## Why not a wraparound `ins`
//!
//! #1205 left `m.` out on the grounds that a circular reference wraps, so
//! position 1 has a valid 5' neighbour and the answer ought to be a wraparound
//! description. That rationale does not survive contact with #129, which
//! established — and `issue_129_mt_circular_wraparound.rs` pins — that ferro
//! **rejects** `m.<high>_<low>ins` by design: the spec's reversed-range
//! exception (`deletion.md:17`, SVD-WG006) is granted to `del`/`dup`, and
//! `delins` inherits it by composition, but `insertion.md` is silent so the
//! general 5'→3' rule applies. Both mutalyzer and biocommons reject all
//! reversed ranges. So `m.16569_1insCG` is not an available answer.
//!
//! The single-position `delins` is, and it needs no reversed range at all:
//! it is the same coordinate identity `insertion_to_boundary_delins` already
//! applies on the other four axes. #129 also disabled 3'-rule shifting across
//! the origin, so the shuffle coming to rest at the terminus is already the
//! intended behaviour — what was broken is only how that resting place was
//! rendered.
//!
//! ## `o.` is out of scope
//!
//! Circular `o.` variants are returned unchanged by `normalize_core` — no
//! shuffle runs at all — so they cannot reach this clamp. Pinned below.

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};
use std::path::{Path, PathBuf};

/// Probe contig: `C`×4, `A`×7, `G`×4 — 15 bases, so both ends sit inside one
/// fetch window and both bounds are exercisable. Deliberately the same shape
/// as #1205's contig, so the `m.` expectations below are directly comparable
/// to the `g.` ones and any divergence between the axes is visible.
const CONTIG: &str = "CCCCAAAAAAAGGGG";

/// The mitochondrial accession, which is what routes these through
/// `normalize_mt` — including the `g.`-spelled inputs, which `normalize_core`
/// coerces to the mito path when the accession `is_mitochondrial()`.
const MT: &str = "NC_012920.1";

fn normalize(input: &str, dir: ShuffleDirection) -> String {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(MT, CONTIG);
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
/// Re-parsing is the property under test as much as the string is: the 5'
/// defect emitted a form ferro's own parser rejects, and the wraparound
/// spelling this fix deliberately avoids would be rejected too (#129).
fn assert_canonical(input: &str, dir: ShuffleDirection, expected: &str) {
    let once = normalize(input, dir);
    assert_eq!(once, expected, "{input} ({dir:?}) canonical form");
    parse_hgvs(&once).unwrap_or_else(|e| panic!("{once} must re-parse: {e}"));
    assert_eq!(normalize(&once, dir), once, "{once} must be a fixed point");
}

/// Apply "replace 1-based inclusive `start..=end` of [`CONTIG`] with `insert`",
/// returning the edited sequence. `start > end` denotes a pure insertion
/// between the two flanking positions.
///
/// Deliberately dumb string splicing, independent of the normalizer: it is
/// what lets the equivalence assertion prove the rewrite preserves the
/// haplotype rather than merely restating what the normalizer did.
fn splice(start: usize, end: usize, insert: &str) -> String {
    let (before, rest) = CONTIG.split_at(start - 1);
    let after = &rest[(end + 1).saturating_sub(start)..];
    format!("{before}{insert}{after}")
}

// ---------------------------------------------------------------------------
// The defect at each mitochondrial terminus
// ---------------------------------------------------------------------------

/// 5' saturation. The shuffle drives the insertion to rest immediately 5' of
/// `m.1`, where the coordinate identity gives `delete ref[1]` (= `C`),
/// `insert CG ++ C`.
#[test]
fn five_prime_saturated_insertion_clamps_at_mito_start() {
    assert_canonical(
        &format!("{MT}:m.1_2insGC"),
        ShuffleDirection::FivePrime,
        &format!("{MT}:m.1delinsCGC"),
    );
}

/// 3' saturation, the mirror: `delete ref[15]` (= `G`), `insert G ++ AGG`.
///
/// The alt is `GGA` rather than a clean homopolymer because a tandem alt never
/// reaches this shape: with `CodonGate::NotApplicable` on the mito axis, two or
/// more copies become repeat notation and a single copy becomes a `dup`, both
/// before any saturation. `GGA` shifts twice through the `G`×4 tract and then
/// mismatches, which is what leaves it an insertion resting on the bound.
#[test]
fn three_prime_saturated_insertion_clamps_at_mito_end() {
    assert_canonical(
        &format!("{MT}:m.13_14insGGA"),
        ShuffleDirection::ThreePrime,
        &format!("{MT}:m.15delinsGAGG"),
    );
}

/// A `g.`-spelled variant on a mitochondrial accession is affected too:
/// `normalize_core` coerces `HV::Genome` to the mito path when the accession
/// `is_mitochondrial()` (#487), so #1205's contig clamp never saw it and the
/// fix has to land in `normalize_mt` to cover this spelling.
///
/// Note the output is `m.`, not `g.` — the coercion is what #487 requires.
#[test]
fn genome_spelled_mito_insertion_is_clamped_too() {
    assert_canonical(
        &format!("{MT}:g.1_2insGC"),
        ShuffleDirection::FivePrime,
        &format!("{MT}:m.1delinsCGC"),
    );
}

/// The rewrite must describe the same haplotype, not merely a parseable one.
///
/// Both sides are spliced into the reference by hand here, so this does not go
/// through the normalizer and cannot be satisfied by an internally-consistent
/// but wrong rewrite.
#[test]
fn the_clamped_delins_describes_the_same_sequence() {
    // 5': insert `GC` between m.1 and m.2  ==  replace m.1 with `CGC`.
    assert_eq!(splice(2, 1, "GC"), splice(1, 1, "CGC"));
    // 3': insert `GGA` between m.13 and m.14  ==  replace m.15 with `GAGG`.
    //
    // The input is stated at m.13_14 but the shuffled payload is `AGG`; the
    // 3'-shifted insertion between m.15 and m.16 is the same haplotype, which
    // the first equality pins before the delins comparison.
    assert_eq!(splice(14, 13, "GGA"), splice(16, 15, "AGG"));
    assert_eq!(splice(16, 15, "AGG"), splice(15, 15, "GAGG"));
}

/// The clamp must not reach for the wraparound spelling #129 rejects.
///
/// Stated as its own test because it is the judgement call this issue records:
/// the tempting "circular references wrap, so use the last base as the 5'
/// neighbour" answer would emit `m.15_1ins…`, which ferro's parser refuses.
/// Asserting on the absence of a reversed range is weaker than the exact-form
/// assertions above, but it fails loudly if a later change swaps strategies.
#[test]
fn the_clamp_does_not_emit_a_wraparound_insertion() {
    for (input, dir) in [
        (format!("{MT}:m.1_2insGC"), ShuffleDirection::FivePrime),
        (format!("{MT}:m.13_14insGGA"), ShuffleDirection::ThreePrime),
    ] {
        let out = normalize(&input, dir);
        assert!(
            !out.contains("ins") || out.contains("delins"),
            "{input} must not come back as a bare ins at the origin, got {out}",
        );
        // A reversed range is the specific shape #129 rejects; `15_1` and
        // `1_15`-with-a-wrap are both spelled with an underscore, so the
        // exact-form tests above carry the precision. This checks the parser
        // agrees, which is the property #129 actually pins.
        parse_hgvs(&out).unwrap_or_else(|e| panic!("{out} must re-parse: {e}"));
    }
}

// ---------------------------------------------------------------------------
// Controls: the clamp must not fire away from the termini
// ---------------------------------------------------------------------------

/// The same alt one position 5' of the saturating case comes to rest inside the
/// contig, so it keeps two valid anchors and must stay an insertion.
#[test]
fn interior_insertion_is_left_alone() {
    assert_canonical(
        &format!("{MT}:m.12_13insGGA"),
        ShuffleDirection::ThreePrime,
        &format!("{MT}:m.14_15insAGG"),
    );
}

/// …and on the 5' side.
#[test]
fn five_prime_interior_insertion_is_left_alone() {
    assert_canonical(
        &format!("{MT}:m.3_4insGC"),
        ShuffleDirection::FivePrime,
        &format!("{MT}:m.2_3insCG"),
    );
}

/// A tandem expansion whose tract runs to the contig start becomes repeat
/// notation, which is not an insertion and so is not the clamp's business.
/// Guards against the clamp firing on position alone rather than on shape.
#[test]
fn repeat_notation_reaching_the_mito_start_is_untouched() {
    assert_canonical(
        &format!("{MT}:m.2_3insCC"),
        ShuffleDirection::FivePrime,
        &format!("{MT}:m.1_4C[6]"),
    );
}

/// An insertion that cannot shuffle at all is unaffected, even written against
/// the first bases of the contig.
#[test]
fn non_shuffling_insertion_at_the_mito_start_is_untouched() {
    assert_canonical(
        &format!("{MT}:m.1_2insCG"),
        ShuffleDirection::FivePrime,
        &format!("{MT}:m.1_2insCG"),
    );
}

/// `o.` is returned unchanged by `normalize_core` — no shuffle runs, so no
/// saturation and nothing for the clamp to fix. Pins that this fix does not
/// extend there by reflex; a real circular normalizer is #466's circular
/// candidate (the closed #951 was the previous tracker).
#[test]
fn circular_o_axis_is_untouched() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("J01749.1", CONTIG);
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("J01749.1:o.1_2insGC").expect("parse");
    let rendered = normalizer
        .normalize(&variant)
        .expect("normalize")
        .to_string();
    assert_eq!(
        rendered, "J01749.1:o.1_2insGC",
        "o. variants pass through unchanged",
    );
}

// ---------------------------------------------------------------------------
// The real rCRS, reference-gated
// ---------------------------------------------------------------------------
//
// The mock contig above proves the clamp's mechanics; this proves the defect
// was reachable on the actual mitochondrial reference rather than only on a
// contrived one. The rCRS terminal bases are what make it reachable without a
// contrived contig: it ends `…ACATCACGATG` and begins `GATCACAGGT…`, so a `G`
// sits on both sides of the origin and a payload can shuffle onto the terminus.
//
// Skips unless `FERRO_MANIFEST` points at a prepared reference, the same gate
// the other reference-backed axis tests use.

fn manifest_path() -> Option<PathBuf> {
    // FERRO_MANIFEST, when set, is authoritative — no fallback to well-known
    // paths, so CI can disable the runner with `FERRO_MANIFEST=/nonexistent`.
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let p = PathBuf::from(path);
        return if p.exists() { Some(p) } else { None };
    }
    let p = Path::new("benchmark-output/manifest.json");
    if p.exists() {
        return Some(p.to_path_buf());
    }
    None
}

/// Check a group of `(input, expected)` cases that share one shuffle direction
/// against the real rCRS, or skip when no manifest is available.
///
/// Cases are grouped by direction because the direction is fixed when the
/// `Normalizer` is constructed and `MultiFastaProvider` is not `Clone` — so one
/// group is one manifest load. (The other reference-backed tests dodge this with
/// a hand-rolled `Arc<MultiFastaProvider>` newtype forwarding all of
/// `ReferenceProvider`; two copies of that boilerplate already exist and a third
/// is not worth two cases.)
fn check_rcrs(dir: ShuffleDirection, cases: &[(&str, &str)]) {
    let Some(path) = manifest_path() else {
        eprintln!("issue_1217: skipping — no manifest at FERRO_MANIFEST");
        return;
    };
    let provider = ferro_hgvs::MultiFastaProvider::from_manifest(&path)
        .unwrap_or_else(|e| panic!("from_manifest({}) failed: {e}", path.display()));
    let normalizer =
        Normalizer::with_config(provider, NormalizeConfig::default().with_direction(dir));

    for (input, expected) in cases {
        let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        let out = normalizer
            .normalize(&variant)
            .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
            .to_string();
        assert_eq!(&out, expected, "{input} ({dir:?}) on the real rCRS");
        parse_hgvs(&out).unwrap_or_else(|e| panic!("{out} must re-parse: {e}"));
    }
}

/// 3' (the default direction) on the real rCRS: delete `m.16569` (= `G`),
/// insert `G ++ CG`. The `g.`-spelled input is the same case routed through
/// #487's coercion.
#[test]
fn real_rcrs_three_prime_terminal_insertion_is_clamped() {
    check_rcrs(
        ShuffleDirection::ThreePrime,
        &[
            (
                "NC_012920.1:m.16568_16569insGC",
                "NC_012920.1:m.16569delinsGCG",
            ),
            (
                "NC_012920.1:g.16568_16569insGC",
                "NC_012920.1:m.16569delinsGCG",
            ),
        ],
    );
}

/// 5' on the real rCRS: delete `m.1` (= `G`), insert `GC ++ G`.
#[test]
fn real_rcrs_five_prime_terminal_insertion_is_clamped() {
    check_rcrs(
        ShuffleDirection::FivePrime,
        &[("NC_012920.1:m.1_2insCG", "NC_012920.1:m.1delinsGCG")],
    );
}
