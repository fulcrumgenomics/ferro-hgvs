//! Audit for issue #221 — 3+ variants in one allele bracket where the
//! merge pass does NOT apply (#81 item C5).
//!
//! Background.
//! The cis-allele merge pass (#72 / #80) collapses adjacent same-coord
//! sub-variants: `g.[1000G>A;1001A>C]` → `g.1000_1001delinsAC`. Three
//! consecutive adjacent subs also merge: `g.[1000G>A;1001A>C;1002C>T]`
//! → `g.1000_1002delinsACT`. C5 asks for symmetric coverage — what
//! happens when the merge pass cannot fire because the sub-variants
//! aren't adjacent, or aren't combinable, or aren't all the same edit
//! kind, and the cis compound has 3+ members.
//!
//! Invariants pinned here:
//!   * Non-adjacent triples / quadruples / 5-way compounds round-trip
//!     unchanged through parse + normalize — the merge pass does not
//!     bridge gaps.
//!   * Mixed-edit-kind compounds (sub + del + dup + inv + delins)
//!     round-trip unchanged when the kinds are not merge-eligible AND the
//!     members are already in genomic order (all inputs here are ascending).
//!   * Single-accession cis members are rendered in genomic (coordinate)
//!     order regardless of input order (#1098) — see the `ordering` module.
//!     (This supersedes the original audit's "input order preserved verbatim"
//!     invariant, which pinned the pre-#1098 order-dependent behavior.)
//!   * An adjacent pair surrounded by non-adjacent siblings merges
//!     correctly without touching the siblings (selective merge).
//!   * Each sub-variant 3'-shifts against its own neighborhood,
//!     even when sandwiched between non-mergeable siblings.
//!   * Separators (`;` for cis, `(;)` for unknown phase) are not
//!     interchangeable inside one bracket — mixing them rejects.
//!     Double-semicolons and trailing semicolons reject.
//!
//! The original audit (#221) landed no production-code changes and pinned the
//! then-current behavior. The `ordering` module was later updated for #1098,
//! which added genomic-order sorting of single-accession cis members; the
//! other sections here were already in ascending order, so they still
//! round-trip unchanged and continue to guard against regressions.

use ferro_hgvs::hgvs::variant::{AllelePhase, HgvsVariant};
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{MockProvider, Normalizer};

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let n = normalizer
        .normalize(&v)
        .unwrap_or_else(|e| panic!("normalize failed for `{}`: {}", input, e));
    format!("{}", n)
}

fn assert_parse_and_normalize_round_trip(input: &str) -> HgvsVariant {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let parsed_out = format!("{}", v);
    assert_eq!(
        parsed_out, input,
        "parse round-trip mismatch:\n  input:  `{}`\n  output: `{}`",
        input, parsed_out
    );
    let normalizer = Normalizer::new(MockProvider::new());
    let normalized = normalizer
        .normalize(&v)
        .unwrap_or_else(|e| panic!("normalize failed for `{}`: {}", input, e));
    let norm_out = format!("{}", normalized);
    assert_eq!(
        norm_out, input,
        "normalize round-trip mismatch:\n  input:  `{}`\n  output: `{}`",
        input, norm_out
    );
    normalized
}

fn expect_phase(variant: &HgvsVariant, phase: AllelePhase, count: usize) {
    let allele = match variant {
        HgvsVariant::Allele(a) => a,
        other => panic!("expected AlleleVariant, got {:?}", other),
    };
    assert_eq!(allele.phase, phase, "wrong phase");
    assert_eq!(allele.variants.len(), count, "wrong sub-variant count");
}

// =============================================================================
// SECTION 1 — non-adjacent 3-way (sub / del / dup / inv / delins)
// =============================================================================
//
// The merge pass requires adjacency (`prev.end + 1 == next.start`). All
// of these have gaps and must round-trip verbatim.

mod nonadjacent_three_way {
    use super::*;

    #[test]
    fn three_substitutions_round_trip() {
        let v = assert_parse_and_normalize_round_trip("NC_000001.11:g.[1000G>A;2000A>C;3000C>T]");
        expect_phase(&v, AllelePhase::Cis, 3);
    }

    #[test]
    fn mixed_kinds_round_trip() {
        let v = assert_parse_and_normalize_round_trip(
            "NC_000001.11:g.[1000G>A;2000_2001del;3000_3001dup]",
        );
        expect_phase(&v, AllelePhase::Cis, 3);
    }

    #[test]
    fn inversion_in_the_middle_round_trips() {
        let v =
            assert_parse_and_normalize_round_trip("NC_000001.11:g.[1000G>A;2000_2005inv;3000C>T]");
        expect_phase(&v, AllelePhase::Cis, 3);
    }

    #[test]
    fn three_delins_round_trip() {
        let v = assert_parse_and_normalize_round_trip(
            "NC_000001.11:g.[1000delinsAT;2000delinsGC;3000G>A]",
        );
        expect_phase(&v, AllelePhase::Cis, 3);
    }
}

// =============================================================================
// SECTION 2 — 4-way and 5-way non-adjacent
// =============================================================================
//
// Real-world haplotype reports can list many variants on a single
// chromosome. Pins that the bracket parser, merge pass, and normalize
// pipeline all scale beyond three sub-variants.

mod multi_way {
    use super::*;

    #[test]
    fn four_way_nonadjacent_round_trips() {
        let v =
            assert_parse_and_normalize_round_trip("NC_000001.11:g.[100G>A;200A>C;300C>T;400T>G]");
        expect_phase(&v, AllelePhase::Cis, 4);
    }

    #[test]
    fn five_way_nonadjacent_round_trips() {
        let v = assert_parse_and_normalize_round_trip(
            "NC_000001.11:g.[100G>A;200A>C;300C>T;400T>G;500A>G]",
        );
        expect_phase(&v, AllelePhase::Cis, 5);
    }
}

// =============================================================================
// SECTION 3 — selective merging in mixed compounds
// =============================================================================
//
// The merge pass must collapse the adjacent pair without touching the
// non-adjacent siblings. The output is therefore SHORTER than the input
// (one variant fewer) but the non-merged variants stay verbatim.

mod selective_merge {
    use super::*;

    #[test]
    fn adjacent_pair_then_nonadjacent_third() {
        // `[1000G>A;1001A>C;5000G>T]`: positions 1000+1001 are adjacent
        // and merge to `1000_1001delinsAC`; 5000 stands alone.
        assert_eq!(
            normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C;5000G>T]"),
            "NC_000001.11:g.[1000_1001delinsAC;5000G>T]",
        );
    }

    #[test]
    fn nonadjacent_then_adjacent_pair() {
        assert_eq!(
            normalize_to_string("NC_000001.11:g.[100G>A;5000G>A;5001A>C]"),
            "NC_000001.11:g.[100G>A;5000_5001delinsAC]",
        );
    }

    #[test]
    fn nonadjacent_then_adjacent_then_nonadjacent() {
        // Middle pair (`5000G>A;5001A>C`) merges; siblings at 100 and
        // 9000 stand alone. Pin both the rendered shape and the final
        // sub-variant count: the merge must shrink 4-way to 3-way, and
        // the wrapper must stay a Cis AlleleVariant.
        let v = parse_hgvs("NC_000001.11:g.[100G>A;5000G>A;5001A>C;9000C>T]").expect("parse");
        let normalizer = Normalizer::new(MockProvider::new());
        let normalized = normalizer.normalize(&v).expect("normalize");
        assert_eq!(
            format!("{}", normalized),
            "NC_000001.11:g.[100G>A;5000_5001delinsAC;9000C>T]",
        );
        expect_phase(&normalized, AllelePhase::Cis, 3);
    }
}

// =============================================================================
// SECTION 4 — order is preserved (no implicit sort)
// =============================================================================
//
// The spec does not define a canonical order for cis sub-variants.
// ferro preserves input order; descending and shuffled inputs must
// round-trip verbatim through normalize.

mod ordering {
    use super::*;

    /// #1098: single-accession cis members are reordered into ascending
    /// genomic order. These non-adjacent subs do not merge, so member count
    /// is preserved — only the order changes.
    #[test]
    fn descending_position_order_sorts_to_ascending() {
        let out = normalize_to_string("NC_000001.11:g.[300C>T;200A>C;100G>A]");
        assert_eq!(out, "NC_000001.11:g.[100G>A;200A>C;300C>T]");
        expect_phase(&parse_hgvs(&out).unwrap(), AllelePhase::Cis, 3);
    }

    #[test]
    fn shuffled_position_order_sorts_to_ascending() {
        let out = normalize_to_string("NC_000001.11:g.[200A>C;100G>A;300C>T]");
        assert_eq!(out, "NC_000001.11:g.[100G>A;200A>C;300C>T]");
        expect_phase(&parse_hgvs(&out).unwrap(), AllelePhase::Cis, 3);
    }

    /// Five-element scramble — neither ascending nor descending, neither
    /// alphabetical by edit type. Order-independence only holds under all
    /// permutations, so a 3-way shuffle can accidentally match the sort key;
    /// 5 elements in a non-monotonic order forecloses that loophole.
    #[test]
    fn five_way_scrambled_order_sorts_to_ascending() {
        let out = normalize_to_string("NC_000001.11:g.[300C>T;100G>A;500A>G;200A>C;400T>G]");
        assert_eq!(out, "NC_000001.11:g.[100G>A;200A>C;300C>T;400T>G;500A>G]");
        expect_phase(&parse_hgvs(&out).unwrap(), AllelePhase::Cis, 5);
    }

    /// Order-independence: *every* permutation of the same members normalizes to
    /// one canonical string — the core guarantee #1098 restores. All 3! = 6
    /// orderings are generated programmatically (Heap's algorithm) so no
    /// permutation is silently omitted.
    #[test]
    fn all_permutations_normalize_to_same_canonical_string() {
        let members = ["100G>A", "200A>C", "300C>T"];
        let canonical = "NC_000001.11:g.[100G>A;200A>C;300C>T]";
        for perm in permutations(&members) {
            let input = format!("NC_000001.11:g.[{}]", perm.join(";"));
            assert_eq!(
                normalize_to_string(&input),
                canonical,
                "input `{input}` must normalize to the canonical genomic order"
            );
        }
    }

    /// All permutations of `items` (Heap's algorithm), each as a `Vec` of the
    /// original elements. Kept local to the test so no dependency is added.
    fn permutations<'a>(items: &[&'a str]) -> Vec<Vec<&'a str>> {
        fn generate<'a>(k: usize, arr: &mut Vec<&'a str>, out: &mut Vec<Vec<&'a str>>) {
            if k <= 1 {
                out.push(arr.clone());
                return;
            }
            for i in 0..k {
                generate(k - 1, arr, out);
                if k.is_multiple_of(2) {
                    arr.swap(i, k - 1);
                } else {
                    arr.swap(0, k - 1);
                }
            }
        }
        let mut arr = items.to_vec();
        let mut out = Vec::new();
        let n = arr.len();
        generate(n, &mut arr, &mut out);
        out
    }
}

// =============================================================================
// SECTION 5 — normalize-each-independently with a transcript provider
// =============================================================================

fn provider_with_test_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    // Sequence layout: ATGC | AAAAA | CCCCC | GGGGG | TTTTT | AAAAA | CCCCC | GGGGG | TTTTT | AAAAA | CCCCC | GGGGG | T
    // Per-position (1-based): pos 1=A, 2=T, 3=G, 4=C, then 5 As at
    // c.5..=c.9, 5 Cs at c.10..=c.14, 5 Gs at c.15..=c.19, etc. A
    // `c.6dup` (duplicate the A at position 6) 3'-shifts to the
    // rightmost A in the run, emitted as `c.9dup`.
    let sequence = "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

mod normalize_each_independently {
    use super::*;

    /// `c.[6dup;30A>G;55C>G]` — three sub-variants in cis. The first
    /// shifts 3' within its homopolymer; the substitutions are fixed
    /// points. The shift on `c.6dup` must not perturb the siblings,
    /// and the siblings must not block the shift.
    #[test]
    fn middle_subs_do_not_block_outer_dup_3prime_shift() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("NM_TEST.1:c.[6dup;30A>G;55C>G]").expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        assert_eq!(format!("{}", normalized), "NM_TEST.1:c.[9dup;30A>G;55C>G]");
    }

    /// Two sub-variants that BOTH require independent 3'-shifts in
    /// different homopolymer regions, plus a fixed-point sub. The
    /// transcript has 5-A runs at c.5..=c.9 and c.25..=c.29; `c.6dup`
    /// shifts to `c.9dup` (last A in first run); `c.26dup` shifts to
    /// `c.29dup` (last A in second run). Pins that multiple shifting
    /// sub-variants do not interfere with each other.
    #[test]
    fn two_shifting_dups_plus_fixed_point_sub_shift_independently() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("NM_TEST.1:c.[6dup;26dup;55C>G]").expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        assert_eq!(format!("{}", normalized), "NM_TEST.1:c.[9dup;29dup;55C>G]");
    }

    /// Four-way cis with three independent shifts. Pins that scaling
    /// beyond three sub-variants does not break the per-variant
    /// shuffle isolation.
    #[test]
    fn four_way_three_shifts_plus_fixed_point_sub() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("NM_TEST.1:c.[6dup;26dup;46dup;55C>G]").expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        assert_eq!(
            format!("{}", normalized),
            "NM_TEST.1:c.[9dup;29dup;49dup;55C>G]"
        );
    }
}

// =============================================================================
// SECTION 6 — unknown-phase 3+ canonicalization
// =============================================================================
//
// Per `recommendations/DNA/alleles.md` and #148, unknown-phase compounds
// canonicalize to the bare form (no outer brackets). 3+ variant inputs
// must canonicalize the same way; the `(;)` separator must survive at
// every position.

mod unknown_phase_three_plus {
    use super::*;

    #[test]
    fn three_way_bracketed_input_canonicalizes_to_bare() {
        // Bracketed input drifts to bare (the canonical unknown form);
        // bare input round-trips. Both shapes are accepted, the bare
        // shape is canonical. Pin both parse/display and normalize so
        // regressions in either path are caught.
        let parsed = parse_hgvs("NC_000001.11:g.[100G>A(;)200A>C(;)300C>T]").expect("parse");
        assert_eq!(
            format!("{}", parsed),
            "NC_000001.11:g.100G>A(;)200A>C(;)300C>T"
        );
        let normalizer = Normalizer::new(MockProvider::new());
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        assert_eq!(
            format!("{}", normalized),
            "NC_000001.11:g.100G>A(;)200A>C(;)300C>T"
        );
    }

    #[test]
    fn three_way_bare_form_round_trips() {
        let normalizer = Normalizer::new(MockProvider::new());
        let parsed = parse_hgvs("NC_000001.11:g.100G>A(;)200A>C(;)300C>T").expect("parse");
        assert_eq!(
            format!("{}", parsed),
            "NC_000001.11:g.100G>A(;)200A>C(;)300C>T"
        );
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        assert_eq!(
            format!("{}", normalized),
            "NC_000001.11:g.100G>A(;)200A>C(;)300C>T",
        );
    }

    #[test]
    fn four_way_unknown_phase_round_trips() {
        let normalizer = Normalizer::new(MockProvider::new());
        let input = "NC_000001.11:g.100G>A(;)200A>C(;)300C>T(;)400T>G";
        let parsed = parse_hgvs(input).expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        assert_eq!(format!("{}", normalized), input);
    }
}

// =============================================================================
// SECTION 7 — separator integrity
// =============================================================================

mod separator_integrity {
    use super::*;

    /// Mixing `;` (cis) and `(;)` (unknown-phase) inside one bracket
    /// must reject — the spec treats them as distinct phase markers
    /// and gives no semantics for a mixed compound.
    #[test]
    fn mixed_semicolon_and_unknown_separator_rejected() {
        assert!(parse_hgvs("NC_000001.11:g.[100G>A;200A>C(;)300C>T]").is_err());
    }

    /// Empty segment between two semicolons rejects.
    #[test]
    fn double_semicolon_rejected() {
        assert!(parse_hgvs("NC_000001.11:g.[100G>A;;200A>C]").is_err());
    }

    /// Trailing semicolon (empty trailing segment) rejects.
    #[test]
    fn trailing_semicolon_rejected() {
        assert!(parse_hgvs("NC_000001.11:g.[100G>A;200A>C;]").is_err());
    }
}

// =============================================================================
// SECTION 8 — merge must not cross bracket / phase-marker boundaries
// =============================================================================
//
// Two adjacent positions sitting in DIFFERENT bracket groups (trans
// phase) must NOT merge — they're on different alleles. The
// 2-variant case is pinned by
// `merge_consecutive_edits_tests::test_no_merge_trans_phase`. The 3+
// trans-phase / mosaic / chimeric cases below extend the negative
// invariant to the lattice this audit owns.

mod merge_does_not_cross_phase_boundaries {
    use super::*;

    /// `[100G>A];[101A>C];[102C>T]` — three trans groups whose
    /// positions are pairwise adjacent. Merge must not bridge any pair.
    #[test]
    fn trans_three_way_adjacent_positions_no_merge() {
        let v = assert_parse_and_normalize_round_trip("NC_000001.11:g.[100G>A];[101A>C];[102C>T]");
        expect_phase(&v, AllelePhase::Trans, 3);
    }

    /// `100G>A/101A>C/102C>T` — three mosaic cell populations with
    /// adjacent positions. Merge must not bridge.
    #[test]
    fn mosaic_three_way_adjacent_positions_no_merge() {
        let v = assert_parse_and_normalize_round_trip(
            "NC_000001.11:g.100G>A/NC_000001.11:g.101A>C/NC_000001.11:g.102C>T",
        );
        expect_phase(&v, AllelePhase::Mosaic, 3);
    }

    /// `100G>A//101A>C//102C>T` — chimeric counterpart.
    #[test]
    fn chimeric_three_way_adjacent_positions_no_merge() {
        let v = assert_parse_and_normalize_round_trip(
            "NC_000001.11:g.100G>A//NC_000001.11:g.101A>C//NC_000001.11:g.102C>T",
        );
        expect_phase(&v, AllelePhase::Chimeric, 3);
    }
}

// =============================================================================
// SECTION 9 — non-mergeable edit kinds in 3+ compounds
// =============================================================================
//
// Inversions are the canonical non-mergeable kind for adjacent
// compounds (see `merge_consecutive_edits_tests::test_no_merge_inversion_adjacent_to_sub`).
// When an inversion sits adjacent to two mergeable subs, the merge
// pass collapses only the subs and leaves the inversion intact.

mod nonmergeable_kinds {
    use super::*;

    /// `[1000_1005inv;1006G>A;1007A>C]` — inversion at 1000-1005,
    /// adjacent subs at 1006+1007. The inversion stays put (inv does
    /// not merge with subs); the two adjacent subs merge to delins.
    #[test]
    fn inversion_then_adjacent_subs_partial_merge() {
        let out = normalize_to_string("NC_000001.11:g.[1000_1005inv;1006G>A;1007A>C]");
        assert!(
            out.contains("1000_1005inv"),
            "inv was wrongly absorbed; got `{}`",
            out
        );
        assert!(
            out.contains("1006_1007delinsAC"),
            "adjacent subs did not merge; got `{}`",
            out
        );
        // No spanning delins absorbing the inversion.
        assert!(
            !out.contains("1000_1007delins"),
            "inv + subs wrongly merged into spanning delins; got `{}`",
            out
        );
    }

    /// `[1000G>A;1001A>C;1005_1010inv]` — adjacent subs at start +
    /// inversion off the end. Subs merge; inversion untouched.
    #[test]
    fn adjacent_subs_then_inversion_partial_merge() {
        // Pin the full output: a loose `contains` check could pass even
        // if an extra unintended merge crept in.
        let out = normalize_to_string("NC_000001.11:g.[1000G>A;1001A>C;1005_1010inv]");
        assert_eq!(
            out, "NC_000001.11:g.[1000_1001delinsAC;1005_1010inv]",
            "got `{}`",
            out
        );
    }

    /// Sandwich: adjacent subs around a non-adjacent inversion. The
    /// inversion is NOT positionally between the subs (which are far
    /// apart), so neither pair merges. All three sub-variants
    /// round-trip verbatim.
    #[test]
    fn sub_then_distant_inv_then_sub_no_merge() {
        let v = assert_parse_and_normalize_round_trip("NC_000001.11:g.[100G>A;500_510inv;9000C>T]");
        expect_phase(&v, AllelePhase::Cis, 3);
    }
}
