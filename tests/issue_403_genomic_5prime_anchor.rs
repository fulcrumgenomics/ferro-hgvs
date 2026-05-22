//! Regression for #403.
//!
//! Surfaced by `tests/fixtures/biocommons-normalize/baseline-failures/normalized.txt`:
//! `NC_000006.11:g.49917121_49917122insGA` (5'+cross) should canonicalise
//! to `NC_000006.11:g.49917121_49917122dup` per biocommons; ferro
//! currently emits `NC_000006.11:g.49917122_49917123dup` (the 3'-anchored
//! form).
//!
//! Reproducer uses a synthetic genomic contig (`NC_TEST.1`, 250 bp) whose
//! bases around the insertion site mirror the real reference at
//! `NC_000006.11:49917119-49917124`: `AAAGAA` (a single isolated `AG`
//! pair with `A`-only flanks). On that layout, both rotations of the
//! inserted `GA` find a `ref_count=1` tandem at different anchor
//! positions — the rotation-selection tie-break in
//! `insertion_to_duplication` must be direction-aware so the
//! `5'-direction` path picks the more-5' anchor and the `3'-direction`
//! path picks the more-3' one.

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

fn provider() -> MockProvider {
    let mut p = MockProvider::new();
    // Ref is 250 bp so the genomic normalize-path ±100 window around
    // 1-based positions 101-102 fits inside the contig. Bases at
    // 1-based positions 99-104 = `AAAGAA`, with the insertion at
    // 101_102 (between the middle `A` and the only `G`). The flanking
    // sequence is `C`-only so no other AG / GA tandems are introduced.
    let mut seq = "C".repeat(98);
    seq.push_str("AAAGAA"); // 1-based positions 99-104
    seq.push_str(&"C".repeat(146)); // 1-based positions 105-250
    assert_eq!(seq.len(), 250);
    p.add_genomic_sequence("NC_TEST.1", seq);
    p
}

fn normalize(input: &str, direction: ShuffleDirection, cross: bool) -> String {
    let mut cfg = NormalizeConfig::default().with_direction(direction);
    if cross {
        cfg = cfg.allow_crossing_boundaries();
    }
    let normalizer = Normalizer::with_config(provider(), cfg);
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{normalized}")
}

#[test]
fn five_prime_cross_emits_five_prime_anchored_dup() {
    // Most-5' anchor of the AG-tandem dup. Same shape as the failing
    // biocommons row `NC_000006.11:g.49917121_49917122insGA`.
    let out = normalize(
        "NC_TEST.1:g.101_102insGA",
        ShuffleDirection::FivePrime,
        true,
    );
    assert_eq!(out, "NC_TEST.1:g.101_102dup");
}

#[test]
fn three_prime_cross_emits_three_prime_anchored_dup() {
    // Sanity: under 3' direction the equivalent dup is at the 3'-most
    // tandem anchor — the `GA` rotation at 1-based 102-103. Same input
    // as the 5'-direction case, opposite direction → opposite tie-break.
    let out = normalize(
        "NC_TEST.1:g.101_102insGA",
        ShuffleDirection::ThreePrime,
        true,
    );
    assert_eq!(out, "NC_TEST.1:g.102_103dup");
}

#[test]
fn five_prime_no_cross_emits_five_prime_anchored_dup() {
    // No crossing boundaries: this input is entirely within the contig,
    // no axis boundary involved, so cross=false must agree with cross=true.
    let out = normalize(
        "NC_TEST.1:g.101_102insGA",
        ShuffleDirection::FivePrime,
        false,
    );
    assert_eq!(out, "NC_TEST.1:g.101_102dup");
}
