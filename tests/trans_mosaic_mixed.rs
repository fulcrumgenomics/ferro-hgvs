//! Trans alleles with a mosaic `(;)` tail on one allele (#544).
//!
//! `c.[296T>G;476T>C];[476T>C](;)1083A>C` — a trans allele where the second
//! allele is itself an unknown-phase group: a bracketed member followed by
//! `(;)`-separated siblings. HGVS DNA/alleles.md + mosaicism.

use ferro_hgvs::hgvs::{AllelePhase, AlleleVariant};
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn trans_with_mosaic_tail_round_trips() {
    for s in [
        "NM_004006.2:c.[296T>G;476T>C];[476T>C](;)1083A>C",
        "NM_004006.2:c.[296T>G];[476T>C](;)1083G>C(;)1406del",
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

/// First trans group carries a `(;)` mosaic tail — symmetric to the
/// already-supported second-member tail. When the leading group ends with a
/// `(;)` tail the literal `];[` no longer appears in the input, so the old
/// `starts_with('[') && contains("];[")` gate misrouted these to compound
/// parsing and they failed. #559 / CodeRabbit.
#[test]
fn first_member_mosaic_tail_round_trips() {
    for s in [
        "NC_000002.12:g.[100A>G](;)150del;[200del]",
        "NM_004006.2:c.[296T>G](;)1083A>C;[476T>C]",
        "NR_004006.2:n.[100A>G](;)300del;[200T>C]",
        "NC_012920.1:m.[100A>G](;)300del;[200T>C]",
        "AC_000001.1:o.[100A>G](;)300del;[200T>C]",
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

// ---------------------------------------------------------------------------
// g. / n. / m. / o. axis round-trips (#544 nitpick)
// ---------------------------------------------------------------------------

/// Genomic trans+mosaic round-trips. Exercises the same `[a;b];[c](;)d`
/// shape as the c. test above but for `g.` coordinates.
#[test]
fn trans_with_mosaic_tail_g_round_trips() {
    for s in [
        "NC_000001.11:g.[100A>G;200T>C];[300del](;)400dup",
        "NC_000001.11:g.[100A>G];[200T>C](;)300del(;)400dup",
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

/// Non-coding transcript trans+mosaic round-trips (`n.` axis).
#[test]
fn trans_with_mosaic_tail_n_round_trips() {
    for s in [
        "NR_004006.2:n.[100A>G;200T>C];[300del](;)400dup",
        "NR_004006.2:n.[100A>G];[200T>C](;)300del",
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

/// Mitochondrial trans+mosaic round-trips (`m.` axis, heteroplasmy).
#[test]
fn trans_with_mosaic_tail_m_round_trips() {
    for s in [
        "NC_012920.1:m.[100A>G;200T>C];[300del](;)400dup",
        "NC_012920.1:m.[100A>G];[200T>C](;)300del",
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

/// Circular DNA trans+mosaic round-trips (`o.` axis, SVD-WG006).
#[test]
fn trans_with_mosaic_tail_o_round_trips() {
    for s in [
        "AC_000001.1:o.[100A>G;200T>C];[300del](;)400dup",
        "AC_000001.1:o.[100A>G];[200T>C](;)300del",
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

// ---------------------------------------------------------------------------
// Mixed-reference trans+mosaic Display fix (#544 — Finding 1)
// ---------------------------------------------------------------------------

/// A programmatically-constructed trans allele where one member is an
/// `Unknown`-phase sub-allele and the members have different reference
/// accessions (so `trans_compact_anchor` returns `None`, forcing the
/// mixed-reference expanded path).
///
/// Before the fix, the else-branch used `write!(f, "[{}]", v)?`, which
/// called the full `Display` on the Unknown-phase sub-allele (producing the
/// compact bracketless form `ACC:c.B(;)C`) and then wrapped it in `[…]`,
/// yielding the syntactically-invalid `[ACC:c.B(;)C]`. The fix replaces
/// that with `write_trans_member_bracketed`, which recognises the
/// Unknown-phase inner allele and emits `[B](;)C` with the bracket only
/// around the first member — exactly the trans+mosaic spec shape.
#[test]
fn mixed_accession_trans_mosaic_display_not_double_wrapped() {
    // Build the inner Unknown-phase allele from two separately-parsed
    // leaf variants that share accession NM_001005922.2 (different from
    // the outer leaf below).
    let inner_a = parse_hgvs("NM_001005922.2:c.476T>C").expect("inner leaf a must parse");
    let inner_b = parse_hgvs("NM_001005922.2:c.1083A>C").expect("inner leaf b must parse");
    let unknown_member = HgvsVariant::Allele(AlleleVariant::new(
        vec![inner_a, inner_b],
        AllelePhase::Unknown,
    ));

    // The outer trans allele has two members with different accessions:
    // NM_004006.2 (leaf) vs NM_001005922.2 (Unknown-phase group above).
    // trans_compact_anchor will return None, so the mixed-reference else
    // branch is taken.
    let leaf = parse_hgvs("NM_004006.2:c.296T>G").expect("leaf must parse");
    let trans = HgvsVariant::Allele(AlleleVariant::new(
        vec![leaf, unknown_member],
        AllelePhase::Trans,
    ));
    let display = format!("{trans}");

    // The Unknown-phase member must render as [B](;)C, not [ACC:type.B(;)C].
    // In the mixed-reference path the accession+type prefix is omitted by
    // write_trans_member_bracketed (it calls fmt_loc_edit for each sub),
    // so the expected form is the bare edit without a coord-system prefix.
    assert!(
        !display.contains("[NM_001005922.2:"),
        "display must not wrap the accession-prefixed form in brackets; got: {display}"
    );
    // The first sub of the Unknown-phase member should be bracket-wrapped
    // and the second appended with (;) — the trans+mosaic bracket shape.
    assert!(
        display.contains("[476T>C](;)1083A>C"),
        "Unknown-phase member must render as [476T>C](;)1083A>C; got: {display}"
    );
    // No double brackets anywhere.
    assert!(
        !display.contains("[["),
        "display must not contain double brackets; got: {display}"
    );
}
