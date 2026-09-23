//! A number too large for the type that holds it must be rejected, never wrapped,
//! truncated, or replaced by a default: any of those turns the input into a
//! different, valid-looking variant.
//!
//! Found by the `fuzz_parse_hgvs` target (#2239).

use ferro_hgvs::hgvs::parser::fast_path::{try_fast_path, FastPathResult};
use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
use ferro_hgvs::parse_hgvs;
use rstest::rstest;

const U64_MAX: &str = "18446744073709551615";
const U64_MAX_PLUS_1: &str = "18446744073709551616";
const I64_MAX: &str = "9223372036854775807";
const I64_MAX_PLUS_1: &str = "9223372036854775808";

/// Every edit-grammar slot that holds a `u64`, with `{}` for the number, paired
/// with how the parser renders the input when the number fits.
const U64_SLOTS: &[(&str, &str)] = &[
    (
        "NC_000016.10:g.3243405CAT[{}]",
        "NC_000016.10:g.3243405CAT[{}]",
    ),
    (
        "NC_000016.10:g.3243405CAT[{}_20]",
        "NC_000016.10:g.3243405CAT[{}_20]",
    ),
    (
        "NC_000016.10:g.3243405CAT[10_{}]",
        "NC_000016.10:g.3243405CAT[10_{}]",
    ),
    (
        "NC_000016.10:g.3243405CAT[{}_?]",
        "NC_000016.10:g.3243405CAT[{}_?]",
    ),
    (
        "NC_000016.10:g.3243405CAT[?_{}]",
        "NC_000016.10:g.3243405CAT[?_{}]",
    ),
    (
        "NC_000016.10:g.3243405CAT[({}_20)]",
        "NC_000016.10:g.3243405CAT[({}_20)]",
    ),
    (
        "NC_000016.10:g.3243405CAT[(10_{})]",
        "NC_000016.10:g.3243405CAT[(10_{})]",
    ),
    (
        "NC_000016.10:g.3243405CAT[({}_?)]",
        "NC_000016.10:g.3243405CAT[{}_?]",
    ),
    (
        "NC_000016.10:g.3243405CAT[(?_{})]",
        "NC_000016.10:g.3243405CAT[?_{}]",
    ),
    ("NC_000001.11:g.100_200({})", "NC_000001.11:g.100_200[{}]"),
    (
        "NC_000001.11:g.100_200(10_{})",
        "NC_000001.11:g.100_200[10_{}]",
    ),
    (
        "NC_000001.11:g.100_200({}_?)",
        "NC_000001.11:g.100_200[{}_?]",
    ),
    ("NC_000001.11:g.100_200del{}", "NC_000001.11:g.100_200del{}"),
    (
        "NC_000001.11:g.100_200del{}insA",
        "NC_000001.11:g.100_200del{}insA",
    ),
    ("NC_000001.11:g.100_200dup{}", "NC_000001.11:g.100_200dup{}"),
    (
        "NC_000001.11:g.100_200dup(150_{})",
        "NC_000001.11:g.100_200dup(150_{})",
    ),
    (
        "NC_000001.11:g.100_200copy{}",
        "NC_000001.11:g.100_200copy{}",
    ),
    ("NC_000001.11:g.100_101ins{}", "NC_000001.11:g.100_101ins{}"),
    (
        "NC_000001.11:g.100_101ins10_{}",
        "NC_000001.11:g.100_101ins10_{}",
    ),
    (
        "NC_000001.11:g.100_101ins10_{}inv",
        "NC_000001.11:g.100_101ins10_{}inv",
    ),
    (
        "NC_000001.11:g.100_101ins({})",
        "NC_000001.11:g.100_101ins{}",
    ),
    (
        "NC_000001.11:g.100_101ins({}_?)",
        "NC_000001.11:g.100_101ins({}_18446744073709551615)",
    ),
    (
        "NC_000001.11:g.100_101ins(10_{})",
        "NC_000001.11:g.100_101ins(10_{})",
    ),
    (
        "NC_000001.11:g.100_101ins(10_20)_({}_30)",
        "NC_000001.11:g.100_101ins[10_20;{}_30]",
    ),
    (
        "NC_000001.11:g.100_101ins[(10_{})]",
        "NC_000001.11:g.100_101ins(10_{})",
    ),
    (
        "NC_000001.11:g.100_101ins[10_{}]",
        "NC_000001.11:g.100_101ins10_{}",
    ),
    (
        "NC_000001.11:g.100_101ins[10_{};A]",
        "NC_000001.11:g.100_101ins[10_{};A]",
    ),
    (
        "NP_001120979.1:p.Met(?_{})_Gln40del",
        "NP_001120979.1:p.(?_Met{})_Gln40del",
    ),
    (
        "NP_001120979.1:p.Gln({}_?)_Gln40del",
        "NP_001120979.1:p.(Gln{}_?)_Gln40del",
    ),
    (
        "NP_001120979.1:p.Gln{}(?)_Gln40del",
        "NP_001120979.1:p.(Gln{})_Gln40del",
    ),
];

/// Every protein-extension slot, which holds an `i64`.
const I64_SLOTS: &[(&str, &str)] = &[
    (
        "NP_001120979.1:p.Ter110GlnextTer{}",
        "NP_001120979.1:p.Ter110GlnextTer{}",
    ),
    (
        "NP_001120979.1:p.Ter110Glnext*{}",
        "NP_001120979.1:p.Ter110GlnextTer{}",
    ),
    ("NP_001120979.1:p.Met1ext-{}", "NP_001120979.1:p.Met1ext-{}"),
    (
        "NP_001120979.1:p.Ter110delextTer{}",
        "NP_001120979.1:p.Ter110extTer{}",
    ),
    (
        "NP_001120979.1:p.Ter110delext*{}",
        "NP_001120979.1:p.Ter110extTer{}",
    ),
    (
        "NP_001120979.1:p.Gln40({}_?)",
        "NP_001120979.1:p.Gln40extTer{}",
    ),
    (
        "NP_001120979.1:p.Gln40(?_{})",
        "NP_001120979.1:p.Gln40extTer{}",
    ),
];

#[rstest]
fn a_u64_slot_rejects_a_number_past_u64_max(
    #[values(U64_MAX_PLUS_1, "99999999999999999999999999")] number: &str,
) {
    for (shape, _) in U64_SLOTS {
        let input = shape.replace("{}", number);
        assert!(parse_hgvs(&input).is_err(), "{input} must be rejected");
    }
}

#[test]
fn a_u64_slot_keeps_u64_max_exactly() {
    for (shape, rendered) in U64_SLOTS {
        let input = shape.replace("{}", U64_MAX);
        let parsed = parse_hgvs(&input).unwrap_or_else(|e| panic!("{input} should parse: {e:?}"));
        assert_eq!(
            parsed.to_string(),
            rendered.replace("{}", U64_MAX),
            "{input}"
        );
    }
}

#[rstest]
fn an_i64_slot_rejects_a_number_past_i64_max(
    #[values(I64_MAX_PLUS_1, U64_MAX_PLUS_1)] number: &str,
) {
    for (shape, _) in I64_SLOTS {
        let input = shape.replace("{}", number);
        assert!(parse_hgvs(&input).is_err(), "{input} must be rejected");
    }
}

#[test]
fn an_i64_slot_keeps_i64_max_exactly() {
    for (shape, rendered) in I64_SLOTS {
        let input = shape.replace("{}", I64_MAX);
        let parsed = parse_hgvs(&input).unwrap_or_else(|e| panic!("{input} should parse: {e:?}"));
        assert_eq!(
            parsed.to_string(),
            rendered.replace("{}", I64_MAX),
            "{input}"
        );
    }
}

/// A version that fits in `u32` still takes the fast path; only the one past it
/// defers (the parity of the two is pinned in `fast_path_differential`).
#[test]
fn a_u32_max_version_still_takes_the_fast_path() {
    assert!(matches!(
        try_fast_path("NM_000088.4294967295:c.459A>G"),
        FastPathResult::Success(_)
    ));
    assert!(matches!(
        try_fast_path("NM_000088.4294967296:c.459A>G"),
        FastPathResult::Fallback
    ));
}

/// A multi-base substitution at the top of the coordinate range is refused with
/// a suggested `delins` repair. Widening the anchor into that repair must not
/// overflow: strict parsing refuses it, and lenient parsing keeps the anchor.
#[rstest]
#[case(
    format!("NC_000001.11:g.{U64_MAX}AC>T"),
    format!("NC_000001.11:g.{U64_MAX}delinsT")
)]
#[case(
    format!("NM_000088.3:c.{I64_MAX}AC>T"),
    format!("NM_000088.3:c.{I64_MAX}delinsT")
)]
#[case(
    format!("NR_000001.1:n.{I64_MAX}AC>T"),
    format!("NR_000001.1:n.{I64_MAX}delinsT")
)]
#[case(
    format!("NM_000088.3:r.{I64_MAX}ac>u"),
    format!("NM_000088.3:r.{I64_MAX}delinsu")
)]
fn a_multibase_substitution_at_the_coordinate_limit_keeps_its_anchor(
    #[case] input: String,
    #[case] expected: String,
) {
    assert!(parse_hgvs(&input).is_err(), "{input}");
    let lenient = parse_hgvs_lenient(&input)
        .unwrap_or_else(|e| panic!("{input} should parse leniently: {e:?}"));
    assert_eq!(lenient.result.to_string(), expected, "{input}");
}
