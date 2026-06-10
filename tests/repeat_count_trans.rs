//! Trans alleles written as shared-position repeat counts (#544).
//!
//! `p.Ala2[10];[11]` is a trans allele where both alleles repeat the same
//! unit at the same position, differing only in count — the position+unit is
//! written once and each allele's count is bracketed (`[10];[11]`). HGVS
//! repeated.md. Distinct from `p.[a];[b]` (which brackets whole members).

use ferro_hgvs::hgvs::AllelePhase;
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn repeat_count_trans_round_trips() {
    for s in [
        "NP_003997.1:p.Ala2[10];[11]",
        "NM_004006.3:r.-124_-123[14];[18]",
        "NM_004006.3:r.-124ug[14];[18]",
    ] {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
        let HgvsVariant::Allele(a) = &v else {
            panic!("expected Allele for `{s}`, got {v:?}");
        };
        assert_eq!(a.phase, AllelePhase::Trans, "phase for `{s}`");
        assert_eq!(a.variants.len(), 2, "member count for `{s}`");
        assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
    }
}

/// Three or more alleles exercise the N-member parse loop and the `[1..]`
/// Display loop, not just the two-member base case.
#[test]
fn repeat_count_trans_three_members() {
    let s = "NP_003997.1:p.Ala2[10];[11];[12]";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele for `{s}`, got {v:?}");
    };
    assert_eq!(a.phase, AllelePhase::Trans, "phase for `{s}`");
    assert_eq!(a.variants.len(), 3, "member count for `{s}`");
    assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
}

/// The shared-position repeat-trans path must still run `parse_variant`'s
/// post-parse validators rather than short-circuiting with an early `return`.
/// Both validators the reviewer flagged as bypassed (#544) are exercised:
///
/// - `validate_coordinate_system` (#486): a `c.` (coding) coordinate on an
///   `NR_` non-coding RNA reference is a coordinate-system mismatch.
/// - `validate_no_u_in_dna`: a `U` (RNA base) in a DNA-context repeat unit.
///
/// Each must be rejected even when written in the compact repeat-trans
/// shorthand, not silently accepted by a path that skips the validators.
#[test]
fn repeat_count_trans_runs_post_parse_validators() {
    let coord_mismatch = "NR_003051.3:c.100_102CAG[10];[11]";
    let err = parse_hgvs(coord_mismatch).expect_err(
        "c. coordinate on an NR_ non-coding RNA must be rejected even in repeat-trans shorthand",
    );
    let msg = format!("{err}");
    assert!(
        msg.contains("non-coding RNA"),
        "expected coordinate-system-mismatch error, got: {msg}"
    );

    let u_in_dna = "NC_000001.11:g.123_124AU[10];[11]";
    let err = parse_hgvs(u_in_dna).expect_err(
        "U in a DNA-context repeat unit must be rejected even in repeat-trans shorthand",
    );
    let msg = format!("{err}");
    assert!(
        msg.contains("RNA base") && msg.contains("not valid in a DNA"),
        "expected U-in-DNA rejection, got: {msg}"
    );
}

/// When trailing text follows the `;[count]` continuations, the diagnostic
/// must point at the *post-loop* leftover, not the pre-loop tail (#544 review).
/// After consuming `;[11]`, the real leftover of `…[10];[11](80%)` is `(80%)`,
/// an allele-fraction annotation (W3017) — surfacing that requires passing the
/// loop's `rest`, not the original `remaining` (which still starts with `;[11]`
/// and would only yield a generic trailing-characters error).
#[test]
fn repeat_count_trans_trailing_diagnostic_uses_post_loop_tail() {
    let input = "NP_003997.1:p.Ala2[10];[11](80%)";
    let err =
        parse_hgvs(input).expect_err("trailing (80%) after a repeat-trans member must be rejected");
    let msg = format!("{err}");
    assert!(
        msg.contains("Allele-fraction") && msg.contains("(80%)"),
        "expected the W3017 allele-fraction diagnostic on the post-loop tail, got: {msg}"
    );
    assert!(
        !msg.contains(";[11]"),
        "diagnostic should point at the leftover `(80%)`, not the consumed `;[11]`: {msg}"
    );
}

/// The compact `<pos><unit>[c1];[c2]` Display is gated to RNA and protein
/// (`recommendations/{RNA,protein}/repeated.md`); the DNA axes spell each
/// allele out in full (`DNA/repeated.md:33`). A `g.` shared-repeat trans is
/// still *accepted* by the parser, but it must normalize to the expanded
/// `[member];[member]` spec form rather than the compact shape — pinning that
/// the RNA/protein Display gate keeps DNA on its explicit form.
#[test]
fn repeat_count_trans_dna_axis_renders_expanded() {
    let input = "NC_000001.11:g.123_124AC[10];[11]";
    let expanded = "NC_000001.11:g.[123_124AC[10]];[123_124AC[11]]";
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("must parse `{input}`: {e}"));
    let HgvsVariant::Allele(a) = &v else {
        panic!("expected Allele for `{input}`, got {v:?}");
    };
    assert_eq!(a.phase, AllelePhase::Trans, "phase for `{input}`");
    assert_eq!(a.variants.len(), 2, "member count for `{input}`");
    assert_eq!(
        format!("{v}"),
        expanded,
        "DNA shared-repeat trans must render expanded, not compact"
    );
}
