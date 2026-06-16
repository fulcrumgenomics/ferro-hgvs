#![cfg(feature = "web-service")]
//! Issue #394 item 1 — effect classifier classifies delins frameshift
//! correctly.
//!
//! `src/service/handlers/effect.rs:216-219` previously hardcoded
//! `is_frameshift = true` for every `NaEdit::Delins`, with the comment
//! "Can't determine ref length easily". The information IS available
//! via the position interval (`end - start + 1`). This file pins the
//! spec-correct behavior.
//!
//! # Spec basis
//!
//! `assets/hgvs-nomenclature/docs/recommendations/protein/frameshift.md`:
//! frameshift is a CDS-level concept tied to `net_delta % 3 != 0`,
//! where `net_delta = alt_len - ref_len`.
//!
//! # Span-length computability
//!
//! The new `span_len_from_*_interval` helpers compute span length only
//! when both endpoints are integer-only (no offset / no special marker
//! / non-mixed axis sides). Anything else returns `None` and the
//! delins frameshift signal falls back to `false` — conservative, safer
//! than the previous always-true default.
//!
//! The `analyze_na_edit` helper is `pub` so this integration test can
//! exercise it directly without invoking the async `predict_effect`
//! HTTP handler.

use ferro_hgvs::hgvs::edit::{InsertedSequence, NaEdit};
use ferro_hgvs::service::handlers::effect::{
    analyze_na_edit, span_len_from_cds_interval, span_len_from_genome_interval,
    span_len_from_tx_interval,
};
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn delins_in_frame_genomic_same_length() {
    // g.10_12delinsATG: ref_len = 3, alt_len = 3 → in-frame.
    let v = parse_hgvs("NC_000001.11:g.10_12delinsATG").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "delins");
    assert!(!is_fs, "3->3 delins is in-frame");
}

#[test]
fn delins_frameshift_genomic_3_to_1() {
    // g.10_12delinsA: ref_len = 3, alt_len = 1 → net_delta = -2 → frameshift.
    let v = parse_hgvs("NC_000001.11:g.10_12delinsA").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "delins");
    assert!(is_fs, "3->1 delins (-2) is a frameshift");
}

#[test]
fn delins_in_frame_genomic_3_to_6() {
    // g.10_12delinsATGCAT: ref_len = 3, alt_len = 6 → net_delta = +3 → in-frame.
    let v = parse_hgvs("NC_000001.11:g.10_12delinsATGCAT").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert!(!is_fs, "3->6 delins (+3) is in-frame");
}

#[test]
fn delins_frameshift_cds_2_to_3() {
    // c.10_11delinsATG: ref_len = 2, alt_len = 3 → net_delta = +1 → frameshift.
    let v = parse_hgvs("NM_000001.1:c.10_11delinsATG").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert!(is_fs, "2->3 delins (+1) is a frameshift");
}

#[test]
fn delins_intronic_span_unknown_frameshift_false() {
    // c.10+5_10+7delinsA: both endpoints intronic. Span can't be
    // computed from base alone (needs intron length); span_len = None →
    // conservative `is_frameshift = false`. Pins the new behavior
    // against the prior always-true default.
    let v = parse_hgvs("NM_000001.1:c.10+5_10+7delinsA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(
        span.is_none(),
        "intronic-endpoint interval must produce None span_len",
    );
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(
        !is_fs,
        "intronic delins with unknown span_len must conservatively report \
         is_frameshift = false (NOT the prior hardcoded true)",
    );
}

#[test]
fn delins_single_position_alt1_is_substitution_shape_but_classifier_says_delins() {
    // c.10delinsA: ref_len = 1, alt_len = 1. This is the substitution
    // canonical form per spec, but the parser may preserve it as
    // delins. Whichever form survives, the frameshift signal must be
    // `false`.
    let v = parse_hgvs("NM_000001.1:c.10delinsA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let edit = cv.loc_edit.edit.inner().unwrap();
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(
        !is_fs,
        "1->1 delins is in-frame regardless of canonical form"
    );
}

// --- UTR-only and Tx-axis coverage (review fixes) --------------------------

#[test]
fn delins_5utr_only_in_frame_minus3_to_minus1() {
    // `c.-3_-1delinsAAAAAA`: ref_len = 3 (5'UTR positions -3, -2, -1
    // contiguous), alt_len = 6 → net_delta = +3 → in-frame. Pins the
    // 5'UTR-only branch of `span_len_from_cds_interval` (both negative
    // bases, same utr3=false).
    let v = parse_hgvs("NM_000001.1:c.-3_-1delinsAAAAAA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("5'UTR span computable");
    assert_eq!(
        span, 3,
        "5'UTR span_len = end - start + 1 = -1 - -3 + 1 = 3"
    );
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 3);
    assert!(!is_fs, "5'UTR 3->6 delins (+3) is in-frame");
}

#[test]
fn delins_3utr_only_frameshift_star1_star2_to_1() {
    // `c.*1_*2delinsA`: ref_len = 2 (3'UTR positions *1, *2), alt_len = 1
    // → net_delta = -1 → frameshift. Pins the 3'UTR-only branch
    // (both utr3=true).
    let v = parse_hgvs("NM_000001.1:c.*1_*2delinsA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("3'UTR span computable");
    assert_eq!(span, 2, "3'UTR span_len = end - start + 1 = 2 - 1 + 1 = 2");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 2);
    assert!(is_fs, "3'UTR 2->1 delins (-1) is a frameshift");
}

#[test]
fn delins_mixed_5utr_to_cds_returns_none_span() {
    // `c.-1_1delinsAAA`: 5'UTR-to-CDS straddle. HGVS skips `c.0`, so
    // the tx-frame span is 2 — but the helper cannot derive that from
    // the c.-1 vs c.1 sign asymmetry alone. Conservative-skip per
    // documented contract; pin to `None`.
    let v = parse_hgvs("NM_000001.1:c.-1_1delinsAAA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(
        span.is_none(),
        "mixed 5'UTR-to-CDS span_len must be None (HGVS skips c.0)",
    );
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(
        !is_fs,
        "mixed-axis delins with unknown span_len must conservatively report \
         is_frameshift = false",
    );
}

#[test]
fn delins_tx_axis_in_frame_10_to_12() {
    // `n.10_12delinsAAA`: ref_len = 3, alt_len = 3 → in-frame. Pins
    // `span_len_from_tx_interval`'s happy-path branch which had no
    // direct test coverage in the initial commit (called out as a
    // MAJOR coverage gap in review).
    let v = parse_hgvs("NR_000001.1:n.10_12delinsAAA").unwrap();
    let HgvsVariant::Tx(tv) = v else {
        panic!("expected Tx")
    };
    let span = span_len_from_tx_interval(&tv.loc_edit.location).expect("tx span computable");
    assert_eq!(span, 3);
    let edit = tv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert!(!is_fs, "n. 3->3 delins is in-frame");
}

#[test]
fn delins_tx_axis_frameshift_10_to_12_alt5() {
    // `n.10_12delinsAAAAA`: ref_len = 3, alt_len = 5 → net_delta = +2
    // → frameshift on the n. axis. The frameshift concept is a
    // CDS-only spec construct, but the classifier's `is_frameshift`
    // signal still flips here because the math is shared.
    let v = parse_hgvs("NR_000001.1:n.10_12delinsAAAAA").unwrap();
    let HgvsVariant::Tx(tv) = v else {
        panic!("expected Tx")
    };
    let span = span_len_from_tx_interval(&tv.loc_edit.location).expect("tx span computable");
    let edit = tv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert!(is_fs, "n. 3->5 delins (+2) flips the frameshift signal");
}

#[test]
fn span_len_from_tx_interval_rejects_pure_downstream() {
    // The helper's contract states: "Returns `Some(...)` only when both
    // endpoints are present, integer-only (no offset), and non-downstream
    // (no `n.*N`)." Programmatic construction of a pure-downstream span
    // such as `n.*5_*10` must return `None`, not a base-difference span
    // (the `*` axis is a separate numbering frame from the n. body and
    // `end.base - start.base + 1` is not a meaningful nt count).
    //
    // Pinned at the helper level because the parser does not currently
    // construct `TxInterval` with both endpoints downstream from this
    // path, but downstream callers (e.g. coordinate-system conversions)
    // could.
    use ferro_hgvs::hgvs::interval::TxInterval;
    use ferro_hgvs::hgvs::location::TxPos;
    let interval = TxInterval::new(TxPos::downstream(5), TxPos::downstream(10));
    assert!(
        span_len_from_tx_interval(&interval).is_none(),
        "pure downstream `n.*5_*10` violates the non-downstream contract → None",
    );
}

#[test]
fn span_len_from_tx_interval_rejects_mixed_downstream() {
    // Symmetric guard: one downstream endpoint and one body endpoint
    // straddles the transcript-end boundary. `end.base - start.base + 1`
    // would silently span across that boundary, which is meaningless.
    use ferro_hgvs::hgvs::interval::TxInterval;
    use ferro_hgvs::hgvs::location::TxPos;
    let interval = TxInterval::new(TxPos::new(95), TxPos::downstream(5));
    assert!(
        span_len_from_tx_interval(&interval).is_none(),
        "mixed body↔downstream `n.95_*5` straddles tx-end → None",
    );
}

// --- Non-literal delins inserted sequences (CodeRabbit thread 1) ----------
//
// Previously the Delins arm used `sequence.to_string().len()` to derive
// `alt_len`. For non-literal `InsertedSequence` variants
// (`Reference("...")`, `PositionRange`, `Count`, `Repeat`, etc.), the
// rendered notation length is NOT the nucleotide count. Counting it
// silently lies to the frameshift classifier — e.g. `delins[NC_...:g.500_510]`
// renders to ~25 chars, not 11 nt.
//
// Spec basis: `assets/hgvs-nomenclature/docs/recommendations/DNA/delins.md`
// allows `delinsN[12]`, `delins<cross_ref>`, and `delins<pos_range>`.
// `InsertedSequence::len()` already encodes the spec-correct nt-count
// semantics (Some(n) only when computable). Use it.

#[test]
fn delins_reference_insert_alt_len_unknown_falls_back_to_false() {
    // `g.100_105delins[NC_000002.12:g.500_510]` — Reference cross-ref
    // insert. `InsertedSequence::len()` returns None (requires lookup).
    // Span_len is known (= 6), but alt_len is unknown → conservative
    // `is_frameshift = false`. The previous `to_string().len()` would
    // count ~25 characters and report `(25 - 6).rem_euclid(3) = 1 != 0`
    // → frameshift, which is wrong.
    let v = parse_hgvs("NC_000001.11:g.100_105delins[NC_000002.12:g.500_510]").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    assert_eq!(span, 6, "g.100_105 span_len = 6");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "delins");
    assert_eq!(ref_len, 6, "ref_len comes from span_len");
    assert_eq!(
        alt_len, 0,
        "non-literal alt must be reported as unknown (0), NOT counted from rendered notation",
    );
    assert!(
        !is_fs,
        "Reference-insert delins has unknown alt nt count → conservative is_frameshift = false",
    );
}

#[test]
fn delins_position_range_insert_alt_len_unknown_falls_back_to_false() {
    // `g.100_105delins200_210` — the parser wraps a bare position
    // range into `InsertedSequence::Complex(vec![PositionRange{...}])`
    // (see `parse_inserted_sequence` in `src/hgvs/parser/edit.rs`).
    // `InsertedSequence::len()` returns None for Complex (would need
    // to sum parts), so alt_len is reported as unknown and
    // `is_frameshift` falls back to false. The previous
    // `to_string().len()` would count "200_210" = 7 characters
    // (rendered notation length) — that math is wrong by construction
    // for non-literal inserts, regardless of whether Complex could in
    // principle yield a known length. Pin the conservative-skip
    // contract.
    let v = parse_hgvs("NC_000001.11:g.100_105delins200_210").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 6, "ref_len comes from span_len");
    assert_eq!(
        alt_len, 0,
        "Complex/PositionRange alt is unknown length, NOT len(\"200_210\") = 7",
    );
    assert!(
        !is_fs,
        "PositionRange-insert delins has unknown alt nt count → conservative is_frameshift = false",
    );
}

#[test]
fn delins_count_insert_uses_count_value_not_digit_string_length() {
    // `g.100_102delins12` — Count(12) insert. nt count = 12.
    // `InsertedSequence::len()` returns Some(12). Net delta = +9 →
    // in-frame. The previous `to_string().len()` would count "12" = 2
    // characters → net delta = -1 → frameshift (wrong).
    let v = parse_hgvs("NC_000001.11:g.100_102delins12").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(
        alt_len, 12,
        "Count alt_len = count value, NOT the digit-string render length",
    );
    assert!(!is_fs, "3 -> 12 delins (+9) is in-frame");
}

#[test]
fn delins_repeat_insert_uses_total_repeat_length() {
    // `g.100_102delinsN[12]` — Repeat { count: Exact(12) } insert.
    // nt count = 12. `InsertedSequence::len()` returns Some(12). Net
    // delta = +9 → in-frame. The previous `to_string().len()` would
    // count "N[12]" = 5 characters → net delta = +2 → frameshift
    // (wrong).
    let v = parse_hgvs("NC_000001.11:g.100_102delinsN[12]").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(alt_len, 12, "Repeat alt_len uses count, not display chars");
    assert!(!is_fs, "3 -> 12 delins (+9) is in-frame");
}

#[test]
fn delins_uncertain_range_insert_alt_len_unknown_falls_back_to_false() {
    // Direct struct construction: `Delins { sequence: Range(10, 20),
    // ... }`. `InsertedSequence::len()` returns None for Range
    // (unknown exact length). With span_len = Some(3), the result must
    // still be `is_frameshift = false` — we can't compute net_delta
    // without alt_len.
    let edit = NaEdit::Delins {
        sequence: InsertedSequence::Range(10, 20),
        deleted: None,
        deleted_length: None,
    };
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(&edit, Some(3));
    assert_eq!(kind, "delins");
    assert_eq!(ref_len, 3, "ref_len comes from span_len");
    assert_eq!(alt_len, 0, "Range alt nt count is unknown");
    assert!(
        !is_fs,
        "Range-insert delins with span_len known but alt_len unknown → conservative false",
    );
}

#[test]
fn delins_empty_insert_with_known_span_is_pure_deletion() {
    // `Delins { sequence: Empty, ... }` is shape-equivalent to a
    // deletion. `InsertedSequence::len()` returns Some(0). With
    // span_len = 3, net_delta = -3 → in-frame (deletion of full codon).
    // Pins that `Empty` does flow into the frameshift math (alt_len = 0
    // is *known*, not unknown).
    let edit = NaEdit::Delins {
        sequence: InsertedSequence::Empty,
        deleted: None,
        deleted_length: None,
    };
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(&edit, Some(3));
    assert_eq!(ref_len, 3);
    assert_eq!(alt_len, 0);
    assert!(!is_fs, "3 -> 0 delins (-3) is in-frame");

    // And span_len = 4 with Empty alt → net_delta = -4 → frameshift.
    let (_, is_fs, _, _) = analyze_na_edit(&edit, Some(4));
    assert!(is_fs, "4 -> 0 delins (-4) is a frameshift");
}

// --- Insertion arm sibling-audit: same to_string() bug shape ---------------
//
// `NaEdit::Insertion { sequence: InsertedSequence }` had the same
// rendered-notation-length bug as Delins. Fix it the same way.

#[test]
fn insertion_count_uses_count_value_not_digit_string_length() {
    // `g.100_101ins12` — Count(12) insert. The insertion arm's
    // frameshift heuristic is `len % 3 != 0`. nt count = 12 → in-frame.
    // The previous `to_string().len()` would count "12" = 2 characters
    // → frameshift (wrong).
    let v = parse_hgvs("NC_000001.11:g.100_101ins12").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, _, alt_len) = analyze_na_edit(edit, None);
    assert_eq!(kind, "insertion");
    assert_eq!(alt_len, 12, "Count alt_len = count value");
    assert!(!is_fs, "ins of 12 nt is in-frame");
}

#[test]
fn insertion_unknown_length_falls_back_to_false() {
    // `Insertion { sequence: Range(10, 20) }` — length unknown.
    // Conservative `is_frameshift = false`. The previous code would
    // return `("?_?")".len() % 3 != 0` which is meaningless.
    let edit = NaEdit::Insertion {
        sequence: InsertedSequence::Range(10, 20),
    };
    let (_, is_fs, _, alt_len) = analyze_na_edit(&edit, None);
    assert_eq!(alt_len, 0, "Range alt nt count is unknown");
    assert!(
        !is_fs,
        "ins with unknown alt nt count → conservative is_frameshift = false",
    );
}
