//! #1282 — a cis allele whose payload lands before the first base.
//!
//! A 5'-shuffled cis allele can denote "the reference with bases added at the
//! very front". That is a real variant, but it has no `ins` spelling: the
//! anchors would be `g.0_1ins`, and no HGVS axis has a position 0. Before this
//! fix the derived member reached `hgvs_pos_to_index` with `start == 0` and the
//! subtraction underflowed — a debug build panicked from inside
//! `Normalizer::normalize` (a `Result`-returning library entry point, so the
//! caller could not catch it without `catch_unwind`), and a release build,
//! where `[profile.release]` sets no `overflow-checks`, wrapped to an index
//! near `usize::MAX` instead.
//!
//! The answer is the boundary-delins identity the rest of `src/normalize/mod.rs`
//! already uses at every other bound: "insert `A'` immediately 5' of `anchor`"
//! is "delete `anchor`, insert `A' ++ ref[anchor]`".

use crate::common::cis_apply_oracle::{apply, normalize_in};
use ferro_hgvs::{parse_hgvs, ShuffleDirection};

/// Leading `T` run, so a 5' shuffle drives the payload to the contig start.
const SEQ: &str = "TTTTTTTTTAATATATTTTA";

/// The reported reproducer and two siblings, all of which panicked.
///
/// Each denotes the reference with an `A` added before base 1: the substitution
/// rewrites base 1 to `A` and the insertion lengthens the `T` run, so the
/// applied sequence is `A` ++ reference. The insertion's own start position is
/// irrelevant to the outcome, which is why all three collapse to one answer —
/// and why fixing only the reported spelling would have missed the class.
const REPRODUCERS: &[&str] = &[
    "TEMPLATE:g.[3_4insT;1T>A]",
    "TEMPLATE:g.[2_3insT;1T>A]",
    "TEMPLATE:g.[5_6insT;1T>A]",
];

#[test]
fn a_payload_before_the_first_base_becomes_a_boundary_delins() {
    for input in REPRODUCERS {
        let actual = normalize_in(SEQ, input, ShuffleDirection::FivePrime);
        assert_eq!(
            actual, "TEMPLATE:g.1delinsAT",
            "{input} should clamp to the 5' boundary delins"
        );
    }
}

/// The rewrite must denote what the input denoted — the point of the whole
/// exercise, and the check that would catch a clamp that merely stops the panic.
///
/// Uses the SPDI-backed apply oracle rather than re-normalizing, so the
/// comparison does not run through the code under test.
#[test]
fn the_boundary_delins_preserves_the_sequence() {
    for input in REPRODUCERS {
        let want = apply(SEQ, input).expect("input applies");
        let actual = normalize_in(SEQ, input, ShuffleDirection::FivePrime);
        let got = apply(SEQ, &actual).expect("output applies");
        assert_eq!(
            got, want,
            "{input} -> {actual} changed the denoted sequence"
        );
        // And it is what it claims: the reference with an `A` in front.
        assert_eq!(got, format!("A{SEQ}"), "{input}");
    }
}

/// The output must be re-readable and stable, so the clamp is a fixed point
/// rather than a one-pass patch.
#[test]
fn the_boundary_delins_re_parses_and_is_idempotent() {
    for input in REPRODUCERS {
        let once = normalize_in(SEQ, input, ShuffleDirection::FivePrime);
        parse_hgvs(&once).unwrap_or_else(|e| panic!("{input} -> {once} does not re-parse: {e}"));
        let twice = normalize_in(SEQ, &once, ShuffleDirection::FivePrime);
        assert_eq!(twice, once, "{input} is not a fixed point");
    }
}

/// A payload that shares an affix with `ref[1]` must still come out canonical.
///
/// `insertion_to_boundary_delins` spells `1delins<A' ++ ref[1]>` without asking
/// whether that reduces. When `A'` shares an affix with `ref[1]` it does:
/// `canonicalize_delins` trims `1delinsTAT` on a leading `T` down to an
/// insertion after base 1, and an `A'` *equal* to `ref[1]` would make it a
/// `dup`. Every other call site of that helper runs after a completed shuffle
/// walk, whose completion property already rules a reducible result out. This
/// one fires on a cis-derived member *before* any shuffle, so it does not
/// inherit that protection — which is why the clamp routes its result back
/// through `normalize_na_edit` instead of returning it directly.
///
/// The reproducers above cannot catch this: their payload is a bare `A` against
/// a leading `T`, which shares nothing and so is canonical either way. These
/// inputs put the substitution on base 2 and the payload's first base on `T`,
/// which is the trimmable shape.
///
/// A sweep of two-member cis alleles over 83 references produced 3569 distinct
/// arrivals at this clamp, 239 of them with an affix-sharing payload, and no
/// `A' == ref[1]` case at all — so the `dup` reduction looks unreachable today.
/// "Looks unreachable" is not a guarantee, and the recursion costs one call, so
/// the clamp does not rely on it.
#[test]
fn a_reducible_boundary_payload_still_lands_canonical() {
    for (input, expected) in [
        ("TEMPLATE:g.[2_3insTT;2T>A]", "TEMPLATE:g.1delinsTAT"),
        ("TEMPLATE:g.[2_3insTT;2T>C]", "TEMPLATE:g.1delinsTCT"),
    ] {
        let once = normalize_in(SEQ, input, ShuffleDirection::FivePrime);
        assert_eq!(once, expected, "{input}");

        // Canonical, so a second pass cannot trim it to an insertion.
        let twice = normalize_in(SEQ, &once, ShuffleDirection::FivePrime);
        assert_eq!(twice, once, "{input} -> {once} is not a fixed point");

        // And the rewrite still denotes what the input denoted.
        let want = apply(SEQ, input).expect("input applies");
        let got = apply(SEQ, &once).expect("output applies");
        assert_eq!(got, want, "{input} -> {once} changed the denoted sequence");
    }
}

/// Neighbouring shapes must be untouched — the clamp fires only at interbase 0.
///
/// Without this, a clamp that fired too eagerly would rewrite valid insertions
/// into `delins` and still pass every test above.
#[test]
fn an_interior_payload_keeps_its_insertion_spelling() {
    for (input, expected) in [
        ("TEMPLATE:g.3_4insT", "TEMPLATE:g.1dup"),
        ("TEMPLATE:g.1T>A", "TEMPLATE:g.1T>A"),
    ] {
        let actual = normalize_in(SEQ, input, ShuffleDirection::FivePrime);
        assert_eq!(actual, expected, "{input} must not be clamped");
        let want = apply(SEQ, input).expect("input applies");
        let got = apply(SEQ, &actual).expect("output applies");
        assert_eq!(
            got, want,
            "{input} -> {actual} changed the denoted sequence"
        );
    }
}

/// Moving the substitution off base 1 takes the allele out of this clamp's
/// scope — the negative control for a two-member input.
///
/// `[3_4insT;2T>A]` differs from the first reproducer only in which base the
/// substitution hits, and that is enough: the derived payload no longer lands
/// before base 1, so interbase 0 is never reached and the boundary delins must
/// not appear. `an_interior_payload_keeps_its_insertion_spelling` above makes
/// the same point for *single*-member inputs; this is the two-member case,
/// which is the shape that produced the panic.
///
/// What this deliberately does **not** assert is the value it normalizes to.
/// That output is separately wrong — it is #1267, where under a 5' shuffle the
/// insertion's junction travels toward the *upstream* substitution and crosses
/// it, silently changing what the allele denotes:
///
/// ```text
/// in   TEMPLATE:g.[3_4insT;2T>A]   applies to  TATTTTTTTTAATATATTTTA
/// out  TEMPLATE:g.2_3insA          applies to  TTATTTTTTTAATATATTTTA
/// ```
///
/// #1259 measured a 5' mirror of the junction clamp turning 80
/// previously-correct outputs into silently wrong ones, so it is not a local
/// fix, and it is being addressed by the sequence-first rewrite rather than
/// here — clamping it in this PR would collide with that work.
///
/// The tripwire for that defect already exists, on its own reference, in
/// `cis_junction_crossing_shift.rs`'s
/// `a_five_prime_insertion_junction_with_an_upstream_sibling_is_a_known_gap`
/// (`g.[4A>G;9_10insA]` -> `g.4_5insG`). Pinning `g.2_3insA` here as well would
/// put a second guard on one defect in a file named for a different issue, so
/// the fix has two tests to reconcile instead of one. This asserts only what is
/// true both before and after #1267 lands.
#[test]
fn a_substitution_off_the_first_base_is_out_of_scope() {
    let input = "TEMPLATE:g.[3_4insT;2T>A]";
    let actual = normalize_in(SEQ, input, ShuffleDirection::FivePrime);
    assert_ne!(
        actual, "TEMPLATE:g.1delinsAT",
        "{input} must not reach the interbase-0 clamp"
    );
    parse_hgvs(&actual).unwrap_or_else(|e| panic!("{input} -> {actual} does not re-parse: {e}"));
}

/// The circular axis takes the same answer, deliberately.
///
/// `m.` is the one axis where "before base 1" has a candidate alternative
/// spelling — a wraparound to the last base. It is not available: #129 (pinned
/// by `issue_129_mt_circular_wraparound.rs`) established that ferro rejects
/// `m.<high>_<low>ins`, because the spec's reversed-range exception
/// (`deletion.md:17`, SVD-WG006) is granted to `del`/`dup` and `insertion.md`
/// is silent. #1217 settled the same question for the 3' terminus and chose the
/// single-position `delins`, which needs no reversed range at all.
///
/// So the clamp firing on `m.` is the established behaviour, not an oversight —
/// this pins it so a future circular normalizer (#466's circular candidate,
/// which supersedes the closed #951) has to decide
/// deliberately rather than change it by accident.
#[test]
fn the_mitochondrial_axis_takes_the_same_boundary_delins() {
    use ferro_hgvs::{MockProvider, NormalizeConfig, Normalizer};

    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_012920.1", SEQ.to_string());
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let variant = parse_hgvs("NC_012920.1:m.[3_4insT;1T>A]").expect("parses");
    let normalized = normalizer
        .normalize(&variant)
        .expect("must not panic or error")
        .to_string();
    assert_eq!(normalized, "NC_012920.1:m.1delinsAT");
    parse_hgvs(&normalized).expect("re-parses");
}
