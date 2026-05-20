//! Regression / convention pin for issue #291: the r.-axis used by
//! `fetch_ref_for_canonical_split`, `normalize_na_edit` (via
//! `normalize_rna::map_in`/`map_out`), and `simple_rna_pos` is
//! **transcript-1-relative** for positive non-UTR bases in this codebase,
//! NOT CDS-relative.
//!
//! Issue #291 was filed based on a claim that "RNA shares the CDS-relative
//! coordinate axis"; investigation showed the claim is incorrect. The r.
//! arm of `fetch_ref_for_canonical_split` is internally consistent with the
//! rest of the r. normalization path: `r.N` for `N > 0` and non-UTR maps
//! directly to tx index `N`, and only `r.*N`/`r.-N` (UTR) translate through
//! `cds_start` / `cds_end` via `rna_to_tx_pos`.
//!
//! These tests pin the convention with a fixture that has a non-trivial
//! 5'UTR (`cds_start = 100`) so the tx-1 axis and a hypothetical
//! CDS-relative axis yield distinct, easily-distinguishable outputs.
//!
//! Closes #291. Follow-up audit to PR #275.
//!
//! Related pins:
//! - `tests/issue_163_rna_utr3_flag.rs` (UTR `*N` translation via cds_end)
//! - `tests/issue_233_rna_cds_consistency.rs` (r./c. Display parity across
//!   edit shapes — relies on the tx-1 axis being internally consistent
//!   for matched-position pairs)
//! - `src/normalize/mod.rs::test_normalize_rna_deletion_shifts_3prime`
//!   (the `r.14del` → `r.15del` shift uses tx-1 indices; the result
//!   halts at the exon-1 right edge per the HGVS exon-junction
//!   exception — see PR #374 / issue #334)

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Build a 200-base, single-exon transcript with `cds_start = 100` and
/// `cds_end = 180`. The transcript layout:
///
/// | tx range  | bases                       | r. region          |
/// |-----------|-----------------------------|--------------------|
/// | 1..=9     | `CCCCCCCCC` (C * 9)         | 5'UTR              |
/// | 10        | `G`                         | 5'UTR              |
/// | 11..=98   | `C * 88`                    | 5'UTR              |
/// | 99        | `T`                         | 5'UTR (r.-1)       |
/// | 100..=108 | `C * 9`                     | CDS (r.1..r.9)     |
/// | 109       | `G`                         | CDS                |
/// | 110..=130 | `G * 21`                    | CDS                |
/// | 131..=180 | `C * 50`                    | CDS                |
/// | 181       | `A`                         | 3'UTR (r.*1)       |
/// | 182..=200 | `T * 19`                    | 3'UTR              |
///
/// Distinguishing features under each candidate axis convention:
///
/// - Tx-1-relative (the actual codebase convention): `r.10` → tx 10 = `G`
///   surrounded by `C`'s; `r.110` → tx 110 inside the G-run at tx 109..=130;
///   `r.1` → tx 1 in the C-run at tx 1..=9.
/// - CDS-relative (the hypothetical "fix" #291 would have applied): `r.10`
///   → tx 109 = `G` adjacent to the G-run; `r.110` → tx 209 (out of range);
///   `r.1` → tx 100 in the C-run at tx 100..=108.
///
/// Output positions differ measurably between the two conventions, so a
/// passing test pins the tx-1 convention unambiguously.
fn provider_with_cds_start_100() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence: String = "CCCCCCCCC".to_string()        // tx 1..=9
        + "G"                                              // tx 10
        + &"C".repeat(88)                                  // tx 11..=98
        + "T"                                              // tx 99
        + &"C".repeat(9)                                   // tx 100..=108
        + "G"                                              // tx 109
        + &"G".repeat(21)                                  // tx 110..=130
        + &"C".repeat(50)                                  // tx 131..=180
        + "A"                                              // tx 181
        + &"T".repeat(19); // tx 182..=200
    assert_eq!(sequence.len(), 200, "fixture length must be 200");
    // Spot-check the byte layout so future edits to the constants above can
    // be caught by the fixture itself rather than by a downstream test.
    let bytes = sequence.as_bytes();
    assert_eq!(bytes[0], b'C', "tx 1 must be 'C'");
    assert_eq!(bytes[9], b'G', "tx 10 must be 'G'");
    assert_eq!(bytes[10], b'C', "tx 11 must be 'C'");
    assert_eq!(bytes[98], b'T', "tx 99 must be 'T' (r.-1 reference)");
    assert_eq!(bytes[99], b'C', "tx 100 must be 'C' (r.1)");
    assert_eq!(bytes[108], b'G', "tx 109 must be 'G' (start of G-run)");
    assert_eq!(bytes[109], b'G', "tx 110 must be 'G'");
    assert_eq!(bytes[129], b'G', "tx 130 must be 'G' (end of G-run)");
    assert_eq!(bytes[130], b'C', "tx 131 must be 'C' (after G-run)");
    assert_eq!(bytes[180], b'A', "tx 181 must be 'A' (r.*1 reference)");
    assert_eq!(bytes[181], b'T', "tx 182 must be 'T'");

    let len = sequence.len() as u64;
    let exons = vec![Exon::new(1, 1, len)];
    let transcript = Transcript::new(
        "NM_TESTCDS100.1".to_string(),
        Some("CDS100_TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(100),
        Some(180),
        exons,
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

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(provider_with_cds_start_100());
    let variant = parse_hgvs(input).expect("parse failed");
    let normalized = normalizer.normalize(&variant).expect("normalize failed");
    format!("{}", normalized)
}

/// Pins that `normalize_na_edit` (via `normalize_rna::map_in`) reads byte
/// `seq[hgvs_pos - 1]` for positive non-UTR `r.` bases — the tx-1
/// convention. `r.10` → tx 10 = `G`, surrounded by `C`'s, so no 3'-shift.
/// On the way out, `map_out` sees tx 10 < cds_start (100) and canonicalizes
/// the position into 5'UTR notation: tx 10 → `r.{10 - 100}` = `r.-90`.
///
/// Two facets of the tx-1 convention are pinned by this one assertion:
/// 1. `map_in` for positive non-UTR `r.N` maps directly to tx `N` (no
///    `cds_start` offset) — otherwise the read byte would be tx 109 = `G`
///    inside the 22-base G-run at tx 109..=130 and the deletion would shift
///    3' to tx 130, producing `r.31del` or similar.
/// 2. `map_out` honors the actual transcript region: tx 10 is in 5'UTR
///    (since `cds_start = 100`), so it must come back out as `r.-90`, not
///    as `r.10`. This is the asymmetry that makes the convention robust —
///    the input notation `r.10` is interpreted as tx 10 regardless of
///    whether tx 10 happens to be inside the CDS, but the canonical output
///    uses the appropriate region marker.
#[test]
fn rna_positive_nonutr_uses_tx1_axis_no_shift() {
    assert_eq!(
        normalize_to_string("NM_TESTCDS100.1:r.10del"),
        "NM_TESTCDS100.1:r.-90del",
        "r.10 must read tx 10 (tx-1 axis), not tx 109 (CDS-rel); since tx \
         10 lies in the 5'UTR (cds_start = 100), map_out re-emits as r.-90"
    );
}

/// Pins that `normalize_na_edit`'s 3'-shift loop uses tx-1 indices on the
/// way in AND `map_out` re-emits the raw tx index (NOT a CDS-relative
/// translation) for positions inside the CDS-proper window. `r.110` → tx
/// 110 inside the G-run at tx 109..=130, shifts 3' to tx 130. tx 130 lies
/// in CDS (`100 <= 130 <= 180`), so `map_out` returns `r.130` verbatim, not
/// `r.{130 - cds_start + 1}` = `r.31`.
///
/// Under a hypothetical CDS-relative axis, `r.110` would map to tx 209,
/// which is past the 200-base transcript end — `map_in` would return `None`
/// and normalization would fall back to the input, yielding `r.110del`
/// unchanged. The output here (`r.130del`, a distinct shifted position) is
/// only reachable under the tx-1 convention.
#[test]
fn rna_positive_nonutr_shift_uses_tx1_axis() {
    assert_eq!(
        normalize_to_string("NM_TESTCDS100.1:r.110del"),
        "NM_TESTCDS100.1:r.130del",
        "r.110del must shift 3' through the tx 109..=130 G-run to r.130del \
         (proves both map_in and map_out skip cds_start translation for \
         positive non-UTR bases)"
    );
}

/// Pins that `r.1` (the lowest positive non-UTR base) is treated as tx 1,
/// NOT as `cds_start = 100`. tx 1 = `C` at the start of a 9-base C-run;
/// `r.1del` shifts 3' through the run to tx 9. tx 9 lies in 5'UTR (since
/// `cds_start = 100`), so `map_out` re-emits as `r.{9 - 100}` = `r.-91`.
///
/// Under a hypothetical CDS-relative axis, `r.1` → tx 100 = `C` at the
/// start of a different 9-base C-run (tx 100..=108); the shift would land
/// at tx 108 and `map_out` would emit `r.108del` (since 100 <= 108 <= 180
/// is CDS-proper). The tx-1 convention's `r.-91del` is unambiguously
/// distinguishable from the hypothetical CDS-rel `r.108del`.
#[test]
fn rna_first_positive_base_is_tx1_not_cds_start() {
    assert_eq!(
        normalize_to_string("NM_TESTCDS100.1:r.1del"),
        "NM_TESTCDS100.1:r.-91del",
        "r.1del must shift through the tx 1..=9 C-run; tx 9 lies in 5'UTR \
         so map_out re-emits as r.-91 (proves r.1 maps to tx 1, not to \
         cds_start)"
    );
}

/// Pins that `fetch_ref_for_canonical_split` reads the r. ref window from
/// the tx-1 axis (NOT through `cds_start`). The variant below is a
/// 3-base `delins` whose ref slice lands inside the CDS-interior G-run
/// (tx 109..=130). Under the tx-1 convention, the ref bytes are `GGG`
/// and the alt `ccu` (normalized to `CCT`) decomposes as `[Inv(0,2),
/// Sub@2 G>U]` via `decompose_delins`'s revcomp scan — the resulting
/// canonical split renders as a `[..inv;..g>u]` cis-allele.
///
/// Under the hypothetical CDS-relative axis, `r.109_111` would index
/// tx `(100 + 109 - 2)..(100 + 111 - 1)` = tx 207..210, which is past
/// the 200-base transcript end. `fetch_ref_for_canonical_split` would
/// return `None`, `apply_canonical_split` would short-circuit to the
/// un-split variant, and Display would emit the input shape unchanged
/// (`r.109_111delinsccu`). The two outputs are mechanically
/// distinguishable, so the assertion below pins the tx-1 axis at the
/// canonical-split entry point even though the parse, normalize, and
/// Display layers are all unchanged from the simpler `del` cases above.
#[test]
fn rna_canonical_split_fetch_uses_tx1_axis() {
    let out = normalize_to_string("NM_TESTCDS100.1:r.109_111delinsccu");
    assert_eq!(
        out, "NM_TESTCDS100.1:r.[109_110inv;111g>u]",
        "fetch_ref_for_canonical_split must slice tx 108..111 (= GGG) for \
         r.109_111 — not tx 207..210 (out-of-range under a hypothetical \
         CDS-relative axis); the inv + sub split is only reachable when \
         the slice actually contains the CDS-interior G-run"
    );
}

/// Pins that UTR positions (`r.*N`) DO translate through `cds_end` via
/// `rna_to_tx_pos` — this is the one place where r.-axis math is NOT
/// transcript-1-relative. `r.*1` → tx `cds_end + 1` = 181 = `A`; next byte
/// tx 182 = `T`; no shift. Result stays `r.*1del`.
///
/// This is the asymmetry the codebase commits to: positive non-UTR bases
/// are tx-1, but `*N` / `-N` UTR bases go through `cds_end` / `cds_start`.
/// Together with the three tests above, this pins the full r.-axis
/// convention end-to-end.
#[test]
fn rna_utr3_position_translates_through_cds_end() {
    assert_eq!(
        normalize_to_string("NM_TESTCDS100.1:r.*1del"),
        "NM_TESTCDS100.1:r.*1del",
        "r.*1 must translate to tx 181 (cds_end + 1); the 'A' there cannot \
         3'-shift against the following T-run because A != T"
    );
}
