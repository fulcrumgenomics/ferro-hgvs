//! Convention pin for the r.-axis used by `fetch_ref_for_canonical_split`,
//! `normalize_na_edit` (via `normalize_rna::map_in`/`map_out`), and
//! `simple_rna_pos`: on a coding transcript `r.` is **CDS-relative** — the
//! same axis as `c.` — so `r.N` maps to tx `cds_start + N - 1` and `r.123`
//! relates to `c.123` (HGVS `background/numbering.md` L58/L61).
//!
//! History: #291 was filed to make r. CDS-relative; PR #304 instead pinned
//! the then-current **transcript-1-relative** behavior as a "convention",
//! deciding on internal-consistency grounds without checking the spec.
//! **#469 corrected this** — r. is CDS-relative — and these tests were
//! rewritten to assert the spec-correct axis, superseding PR #304's pin.
//!
//! The fixture has a non-trivial 5'UTR (`cds_start = 100`) so the
//! CDS-relative axis and the old tx-1 axis yield distinct, easily
//! distinguishable outputs.
//!
//! Related pins:
//! - `tests/issue_163_rna_utr3_flag.rs` (UTR `*N` translation via cds_end)
//! - `tests/issue_233_rna_cds_consistency.rs` (r./c. Display parity across
//!   edit shapes)
//! - `tests/rna_coding_consistency.rs` (r./c. CDS-relative parity, #469)
//! - `src/normalize/mod.rs::test_normalize_rna_deletion_shifts_3prime`
//!   (the `r.10del` → `r.11del` shift, CDS-relative; the result halts at the
//!   exon-1 right edge per the HGVS exon-junction exception — PR #374 / #334)

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

/// `r.` is CDS-relative on a coding transcript (#469): `r.10` = `c.10` =
/// tx `cds_start + 10 - 1` = tx 109, the first base of the 22-base G-run at
/// tx 109..=130. The deletion 3'-shifts through the run to tx 130, which
/// `map_out` emits CDS-relative as `r.{130 - cds_start + 1}` = `r.31`.
/// (Under the superseded tx-1 pin from PR #304 this read tx 10 — an isolated
/// G in the 5'UTR — and emitted `r.-90del`.)
#[test]
fn rna_positive_cds_base_is_cds_relative() {
    assert_eq!(
        normalize_to_string("NM_TESTCDS100.1:r.10del"),
        "NM_TESTCDS100.1:r.31del",
        "r.10 must read tx 109 (= c.10, CDS-relative) and shift 3' through \
         the G-run to tx 130 = r.31"
    );
}

/// The CDS-relative axis means `r.N` and `c.N` are the same position, so a
/// shifting edit normalizes identically across the two coordinate systems
/// (modulo alphabet). `c.10del` on this fixture shifts the same G-run to
/// tx 130 = `c.31`, the c. twin of `r.10del` → `r.31del` above (#469).
#[test]
fn rna_and_cds_numbering_agree() {
    assert_eq!(
        normalize_to_string("NM_TESTCDS100.1:c.10del"),
        "NM_TESTCDS100.1:c.31del",
        "c.10 == r.10 (both CDS-relative); shifts through the G-run to c.31, \
         the c. twin of the r.31del pinned above"
    );
}

/// `r.1` is the first coding base = `cds_start` = tx 100 (#469), at the
/// start of the 9-base C-run tx 100..=108. `r.1del` 3'-shifts through the
/// run to tx 108, emitted CDS-relative as `r.{108 - cds_start + 1}` = `r.9`.
/// (Under the superseded tx-1 pin, `r.1` was tx 1 and this emitted
/// `r.-91del`.)
#[test]
fn rna_first_coding_base_maps_to_cds_start() {
    assert_eq!(
        normalize_to_string("NM_TESTCDS100.1:r.1del"),
        "NM_TESTCDS100.1:r.9del",
        "r.1 must map to cds_start (tx 100); the C-run shift lands at tx 108 \
         = r.9 (proves r.1 == c.1 == cds_start, not tx 1)"
    );
}

/// Pins that `fetch_ref_for_canonical_split` reads the r. ref window
/// CDS-relative (through `cds_start`, exactly like the `c.` arm — #469).
/// The 3-base `delins` below targets `r.10_12` = tx 109..=111 = the `GGG`
/// at the start of the CDS-interior G-run; alt `ccu` (→ `CCT`) decomposes
/// via the revcomp scan into `[Inv; Sub g>u]`, rendering as a cis-allele.
#[test]
fn rna_canonical_split_fetch_is_cds_relative() {
    let out = normalize_to_string("NM_TESTCDS100.1:r.10_12delinsccu");
    assert_eq!(
        out, "NM_TESTCDS100.1:r.[10_11inv;12g>u]",
        "fetch_ref_for_canonical_split must slice tx 109..=111 (= GGG) for \
         r.10_12 via cds_start, mirroring the c. arm"
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
