//! Issue #2018 — a past-`cds_end` `c.` coordinate is non-conformant, and both
//! the lone-position path and the cis-allele path apply one rule, mode-dependent.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/2018>.
//!
//! # The defect
//!
//! `c.22` on a transcript whose CDS ends at `c.21` names no base on the coding
//! axis — `background/numbering.md:21` ends `c.` numbering at "the last
//! nucleotide of the translation termination (stop) codon", and the next base
//! is `*1`. Ferro used to give that one coordinate two answers:
//!
//! ```text
//! NM_SEAM.1:c.22del              -> REFUSED  (W4004), in strict mode
//! NM_SEAM.1:c.[22_23insG;*1G>A]  -> ACCEPTED -> c.21dup, in every mode
//! ```
//!
//! The lone path ran the `normalize_cds` bounds gate; the allele path merged the
//! members first, which converts `c.22` to the same sequence coordinate as
//! `c.*1` with no bounds check, so the coordinate was silently remapped across
//! the CDS/3'UTR seam. The disagreement is the defect, independent of which
//! answer is right (issue #2018).
//!
//! # The ruling
//!
//! `rulings[past-cds-end-coordinate-is-non-conformant]` (decided) with
//! `rulings[absolute-prohibition-enforcement-stage]`: a past-`cds_end` `c.`
//! coordinate is non-conformant, and the enforcement stage is mode-dependent and
//! applied IDENTICALLY on both paths:
//!
//! - **STRICT** — refuse (`W4004`) on both the lone-position and the cis-allele
//!   path.
//! - **LENIENT** — repair `c.N` to `c.*(N - cds_end)` (deterministic) on both
//!   paths, but as a RECORDED correction (`W4004`), never a silent remap. The
//!   output is the repaired `*N` form; the warning discloses it.
//! - **SILENT** — lenient with the message suppressed (silent mode's contract),
//!   so the repair is applied silently on both paths — consistent, since the
//!   lone path behaves the same way in silent mode.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

/// The CDS/3'UTR seam fixture. `cds_end` is `c.21`, so `c.22` is one base past
/// the coding axis and its canonical spelling is `c.*1`; `c.23` is `c.*2`. The
/// 3'UTR opens with `G`, so `c.*1G>A` is well-formed against it. The layout is
/// identical to `insertion_adjacency_defects`' `seam_transcript` on purpose, so
/// the guard flipped there and the tests here speak of the same coordinates.
fn seam_provider() -> MockProvider {
    const UTR5: &str = "CCCCCCCCCC";
    const CDS: &str = "ATGCCCGGGCATGACACCTAA";
    const UTR3: &str = "GTGTGTGTGTGTGTGTGTGT";
    assert_eq!(CDS.len(), 21);
    assert_eq!(UTR3.as_bytes()[0], b'G', "c.*1 must be G");

    let sequence = format!("{UTR5}{CDS}{UTR3}");
    let cds_start = (UTR5.len() + 1) as u64;
    let cds_end = (UTR5.len() + CDS.len()) as u64;
    let tx_len = sequence.len() as u64;

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_SEAM.1".to_string(),
        Some("SEAM_TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(cds_start),
        Some(cds_end),
        vec![Exon::new(1, 1, tx_len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalizer(mode: ErrorMode) -> Normalizer<MockProvider> {
    Normalizer::with_config(
        seam_provider(),
        NormalizeConfig::default().with_error_mode(mode),
    )
}

/// `normalize` (strict) — returns the error string, or panics if it accepted.
#[track_caller]
fn refuse(input: &str) -> String {
    let variant = parse_hgvs(input).expect("parse");
    match normalizer(ErrorMode::Strict).normalize(&variant) {
        Ok(v) => panic!("strict accepted {input:?}, emitting {v}"),
        Err(e) => e.to_string(),
    }
}

/// `normalize_with_diagnostics` in the given mode — returns `(output, warnings)`.
#[track_caller]
fn normalize_diag(mode: ErrorMode, input: &str) -> (String, Vec<String>) {
    let variant = parse_hgvs(input).expect("parse");
    let result = normalizer(mode)
        .normalize_with_diagnostics(&variant)
        .unwrap_or_else(|e| panic!("{mode:?} refused {input:?}: {e}"));
    let warnings = result.warnings.iter().map(|w| format!("{w:?}")).collect();
    (result.result.to_string(), warnings)
}

/// A refusal that names a past-`cds_end` `c.` coordinate as `W4004`.
#[track_caller]
fn is_past_cds_end_refusal(err: &str) {
    assert!(
        err.contains("W4004") && err.contains("PositionPastEnd") && err.contains("cds-end"),
        "expected a W4004 past-cds-end refusal, got: {err}"
    );
}

/// A `W4004 PositionPastEnd` disclosure naming `position` on the `cds-end` bound.
#[track_caller]
fn discloses_past_cds_end(warnings: &[String], position: &str) {
    assert!(
        warnings.iter().any(|w| w.contains("PositionPastEnd")
            && w.contains(&format!("position: {position:?}"))
            && w.contains("bound_kind: \"cds-end\"")),
        "expected a PositionPastEnd disclosure for c.{position}, got: {warnings:?}"
    );
}

#[track_caller]
fn discloses_no_past_cds_end(warnings: &[String]) {
    assert!(
        !warnings.iter().any(|w| w.contains("PositionPastEnd")),
        "expected no PositionPastEnd disclosure, got: {warnings:?}"
    );
}

// ---------------------------------------------------------------------------
// STRICT — both paths refuse with W4004.
// ---------------------------------------------------------------------------

/// The lone-position path: `c.22del` refuses. This is the pre-existing,
/// correct behaviour; it is the yardstick the allele path must match.
#[test]
fn strict_refuses_a_past_cds_end_coordinate_alone() {
    is_past_cds_end_refusal(&refuse("NM_SEAM.1:c.22del"));
}

/// The cis-allele path, collapsing case: `c.[22_23insG;*1G>A]` used to be
/// accepted and remapped to `c.21dup`. It now refuses with the same W4004 the
/// lone position gets.
#[test]
fn strict_refuses_a_past_cds_end_coordinate_in_a_collapsing_cis_allele() {
    is_past_cds_end_refusal(&refuse("NM_SEAM.1:c.[22_23insG;*1G>A]"));
}

/// The cis-allele path, non-collapsing case: `c.[22del;5C>A]` keeps two
/// members, and the past-`cds_end` one still refuses. This proves the fix is
/// not specific to the merge that collapses the members — any past-`cds_end`
/// member is caught.
#[test]
fn strict_refuses_a_past_cds_end_member_in_a_non_collapsing_cis_allele() {
    is_past_cds_end_refusal(&refuse("NM_SEAM.1:c.[22del;5C>A]"));
}

// ---------------------------------------------------------------------------
// LENIENT — both paths repair to the `*N` spelling AND record W4004.
// ---------------------------------------------------------------------------

/// The lone-position path repairs `c.22del` to `c.*1del` and discloses W4004.
/// Pre-existing behaviour, asserted here as the yardstick for the allele path.
#[test]
fn lenient_repairs_and_discloses_a_past_cds_end_coordinate_alone() {
    let (out, warnings) = normalize_diag(ErrorMode::Lenient, "NM_SEAM.1:c.22del");
    assert_eq!(out, "NM_SEAM.1:c.*1del");
    discloses_past_cds_end(&warnings, "22");
}

/// The cis-allele path repairs the past-`cds_end` member to its `c.*N` spelling
/// and records W4004 — the correction is disclosed rather than silent.
#[test]
fn lenient_repairs_and_discloses_a_past_cds_end_member_in_a_cis_allele() {
    let (out, warnings) = normalize_diag(ErrorMode::Lenient, "NM_SEAM.1:c.[22del;5C>A]");
    assert_eq!(out, "NM_SEAM.1:c.[5C>A;*1del]");
    discloses_past_cds_end(&warnings, "22");
}

/// The collapsing allele also discloses W4004 in lenient mode while still
/// emitting its repaired `c.21dup`. `c.22_23ins` names two past-`cds_end`
/// endpoints, so both are disclosed.
#[test]
fn lenient_discloses_past_cds_end_in_a_collapsing_cis_allele() {
    let (out, warnings) = normalize_diag(ErrorMode::Lenient, "NM_SEAM.1:c.[22_23insG;*1G>A]");
    assert_eq!(out, "NM_SEAM.1:c.21dup");
    discloses_past_cds_end(&warnings, "22");
    discloses_past_cds_end(&warnings, "23");
}

// ---------------------------------------------------------------------------
// The two paths AGREE — the property #2018 is about.
// ---------------------------------------------------------------------------

/// The whole point of the issue: the same coordinate gets the same verdict
/// whether it stands alone or has a sibling. Strict refuses both; lenient
/// repairs-and-discloses both.
#[test]
fn the_lone_and_allele_paths_agree_on_a_past_cds_end_coordinate() {
    // STRICT: both refuse.
    is_past_cds_end_refusal(&refuse("NM_SEAM.1:c.22del"));
    is_past_cds_end_refusal(&refuse("NM_SEAM.1:c.[22del;5C>A]"));

    // LENIENT: both repair and disclose.
    let (lone_out, lone_warnings) = normalize_diag(ErrorMode::Lenient, "NM_SEAM.1:c.22del");
    assert_eq!(lone_out, "NM_SEAM.1:c.*1del");
    discloses_past_cds_end(&lone_warnings, "22");

    let (allele_out, allele_warnings) =
        normalize_diag(ErrorMode::Lenient, "NM_SEAM.1:c.[22del;5C>A]");
    assert_eq!(allele_out, "NM_SEAM.1:c.[5C>A;*1del]");
    discloses_past_cds_end(&allele_warnings, "22");
}

// ---------------------------------------------------------------------------
// Controls — an in-bounds 3'UTR member is untouched (no false positive).
// ---------------------------------------------------------------------------

/// A genuine `c.*1` member (in bounds) in a cis allele is not flagged and not
/// refused, in any mode. The gate keys on past-`cds_end`, not on "touches the
/// 3'UTR".
#[test]
fn an_in_bounds_utr3_member_is_never_flagged() {
    for mode in [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent] {
        let (_, warnings) = normalize_diag(mode, "NM_SEAM.1:c.[*1del;5C>A]");
        discloses_no_past_cds_end(&warnings);
    }
}

/// A lone in-bounds `c.*1del` is accepted unchanged in strict mode — the
/// negative control for the refusal helper.
#[test]
fn a_lone_in_bounds_utr3_position_is_accepted_in_strict() {
    let (out, warnings) = normalize_diag(ErrorMode::Strict, "NM_SEAM.1:c.*1del");
    assert_eq!(out, "NM_SEAM.1:c.*1del");
    discloses_no_past_cds_end(&warnings);
}

// ---------------------------------------------------------------------------
// SILENT — the repair is applied without a message, identically on both paths.
// ---------------------------------------------------------------------------

/// Silent mode suppresses the W4004 message (its contract) and applies the
/// repair on both paths, so the two remain consistent. This is the one mode in
/// which the pre-fix allele behaviour was already consistent with the lone
/// path; the fix must not change it.
#[test]
fn silent_repairs_without_a_message_on_both_paths() {
    let (lone_out, lone_warnings) = normalize_diag(ErrorMode::Silent, "NM_SEAM.1:c.22del");
    assert_eq!(lone_out, "NM_SEAM.1:c.*1del");
    discloses_no_past_cds_end(&lone_warnings);

    let (allele_out, allele_warnings) =
        normalize_diag(ErrorMode::Silent, "NM_SEAM.1:c.[22del;5C>A]");
    assert_eq!(allele_out, "NM_SEAM.1:c.[5C>A;*1del]");
    discloses_no_past_cds_end(&allele_warnings);
}
