//! #1631 — the W3013 repeat-label repair must not change the variant.
//!
//! `RNA/repeated.md:22` calls a repeat carrying BOTH a position range and the
//! repeat unit redundant, because the range already fixes the unit. The repair
//! used to act on that by dropping the **unit** and keeping the **range** — but
//! `:20` fixes the reading of a bare range ("`r.123_125[23]` describes a
//! tri-nucleotide repeat of 23 units"), so a bare range denotes one unit whose
//! length is the range's length. Dropping the unit therefore preserves the
//! variant only while the two lengths agree.
//!
//! Measured on `main` before this change, which is what makes the class worth a
//! file of its own:
//!
//! ```text
//! NM_024312.4:r.-125_-123cug[4]  --error-mode lenient -> r.-125_-123[4]  CORRECT (3 nt == 3 nt)
//! NM_024312.4:r.-6_-3g[6]        --error-mode lenient -> r.-6_-3[6]      WRONG   (4 nt vs 1 nt)
//! ```
//!
//! `r.-6_-3g[6]` is a single `g` six times, 6 nucleotides; `r.-6_-3[6]` is a
//! four-nucleotide unit six times, 24 nucleotides. The repaired description
//! denoted four times the input's bases, and `--error-mode silent` produced it
//! with an empty warning vector, so nothing on that path said the variant had
//! moved.
//!
//! The repair now keeps whichever half carries the variant: the unit when the
//! two disagree (reducing the range to its anchor), the range when they agree.
//! Both arms stay W3013 — the input defect is the same, only the surviving half
//! differs.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorType};
use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_silent, parse_hgvs_with_config};

/// The reported case: range 4 nt, unit 1 nt. The unit wins and the range is
/// reduced to its anchor, which is the issue's suggested acceptance.
#[test]
fn a_range_longer_than_its_unit_keeps_the_unit() {
    let result = parse_hgvs_lenient("NM_024312.4:r.-6_-3g[6]").unwrap();
    assert_eq!(result.preprocessed_input, "NM_024312.4:r.-6g[6]");
    assert_eq!(
        result
            .warnings
            .iter()
            .filter(|w| w.error_type == ErrorType::RedundantRepeatLabel)
            .count(),
        1,
        "the repair is still reported as W3013"
    );
}

/// The same on the silent path, which is where the old behaviour was worst: it
/// produced the wrong variant with an empty warning vector.
#[test]
fn the_silent_path_repairs_the_same_way() {
    let result = parse_hgvs_silent("NM_024312.4:r.-6_-3g[6]").unwrap();
    assert!(!result.has_warnings());
    assert_eq!(result.preprocessed_input, "NM_024312.4:r.-6g[6]");
}

/// The sibling case, and the control: `:22`'s own exemplar has a 3-nt range and
/// a 3-nt unit, so the range really does restate the unit and dropping the unit
/// is variant-preserving. This arm must not move.
#[test]
fn a_range_that_matches_its_unit_still_drops_the_unit() {
    let result = parse_hgvs_lenient("NM_024312.4:r.-125_-123cug[4]").unwrap();
    assert_eq!(result.preprocessed_input, "NM_024312.4:r.-125_-123[4]");
}

/// A positive-coordinate range, so the fix is not keyed on 5'UTR numbering.
#[test]
fn the_length_test_holds_on_cds_coordinates() {
    for (input, expected) in [
        ("NM_000088.3:r.100_102cug[4]", "NM_000088.3:r.100_102[4]"),
        ("NM_000088.3:r.100_101cug[4]", "NM_000088.3:r.100cug[4]"),
        ("NM_000088.3:r.100_104cug[4]", "NM_000088.3:r.100cug[4]"),
    ] {
        let result = parse_hgvs_lenient(input).unwrap();
        assert_eq!(result.preprocessed_input, expected, "input {input}");
    }
}

/// A range spanning the `-1`/`1` boundary spans `|start| + end` positions,
/// because `c.`/`r.` numbering has no position zero. `r.-2_1` is three
/// positions (`-2`, `-1`, `1`), so a 3-nt unit beside it IS redundant.
#[test]
fn a_range_crossing_the_utr_boundary_counts_positions_not_arithmetic() {
    let redundant = parse_hgvs_lenient("NM_000088.3:r.-2_1cug[4]").unwrap();
    assert_eq!(redundant.preprocessed_input, "NM_000088.3:r.-2_1[4]");

    // Naive `end - start + 1` would read `-2_2` as three positions too and
    // wrongly call this one redundant; it is four (`-2`, `-1`, `1`, `2`).
    let mismatched = parse_hgvs_lenient("NM_000088.3:r.-2_2cug[4]").unwrap();
    assert_eq!(mismatched.preprocessed_input, "NM_000088.3:r.-2cug[4]");
}

/// Strict still refuses the input — the ruling that `:22` governs settles
/// *whether* to repair, and that is unchanged here. Only the repaired string
/// moves.
#[test]
fn strict_still_refuses_both_shapes() {
    for input in ["NM_024312.4:r.-6_-3g[6]", "NM_024312.4:r.-125_-123cug[4]"] {
        assert!(
            parse_hgvs_with_config(input, ErrorConfig::strict()).is_err(),
            "strict accepted {input}"
        );
    }
}

/// The repaired string is a fixed point of the repair: re-running the lenient
/// path over it neither changes it nor warns again.
#[test]
fn the_repair_is_idempotent_on_both_arms() {
    for input in ["NM_024312.4:r.-6_-3g[6]", "NM_024312.4:r.-125_-123cug[4]"] {
        let once = parse_hgvs_lenient(input).unwrap();
        let twice = parse_hgvs_lenient(&once.preprocessed_input).unwrap();
        assert_eq!(twice.preprocessed_input, once.preprocessed_input);
        assert!(
            !twice
                .warnings
                .iter()
                .any(|w| w.error_type == ErrorType::RedundantRepeatLabel),
            "{input} repaired to a string that still trips W3013"
        );
    }
}
