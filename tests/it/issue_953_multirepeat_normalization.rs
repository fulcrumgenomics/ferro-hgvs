//! Issue #953 — `NaEdit::MultiRepeat` is flagged `needs_normalization == true`,
//! but a compound multi-unit tandem repeat has no canonical 3'-shift target.
//! The flag is `true` for the *validation* half (`validate_multirepeat_tract`),
//! and the shuffle dispatch now has an explicit `MultiRepeat` arm that passes
//! the (validated) edit through unchanged — reconciling the two so the
//! pass-through is intentional and documented, not a silent contradiction.
//!
//! These tests pin that intended behavior: a MultiRepeat normalizes to itself
//! and is idempotent. To make that assertion non-vacuous, each test registers
//! **real reference bases** whose span exactly matches the declared tract, so
//! the normalizer reaches the `MultiRepeat` arm and runs
//! `validate_multirepeat_tract` (a validated pass-through) rather than bailing
//! early on a missing-reference error the way a bare, empty `MockProvider`
//! would. The `GT[2]GC[2]` tract occupies exactly 8 reference bases
//! (`GTGTGCGC`), so the described span is 8 bp wide and the reference at that
//! span literally spells `GTGTGCGC`; a coordinate/base mismatch would surface
//! as a `RefSeqMismatch` (asserted absent here — see
//! `issue_279_multirepeat_partial_validation` for the mismatch direction), so
//! these tests genuinely exercise the validation path.
//!
//! (The parse assertion guards that the input really is a `MultiRepeat`, so the
//! test keeps exercising the reconciled code path.)

use ferro_hgvs::hgvs::edit::NaEdit;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// 256 bp of `N` flanking the variant tract on each side, mirroring
/// `issue_279_multirepeat_partial_validation`. `N` is inert against the
/// byte-equality scans the normalizer / validator uses, so the boundaries
/// never spuriously extend the apparent tract.
const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

fn padded(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}

/// HGVS position of the first base of `core` in `padded(core)`.
/// `PAD` is 256 bp of `N`, so the first core base is the 257th overall
/// (1-based HGVS coordinates).
const CORE_START: u64 = PAD.len() as u64 + 1;

/// Extract the concrete `NaEdit` from a genome/cds variant for a type assertion.
fn edit_of(variant: &HgvsVariant) -> Option<&NaEdit> {
    match variant {
        HgvsVariant::Genome(v) => v.loc_edit.edit.inner(),
        HgvsVariant::Cds(v) => v.loc_edit.edit.inner(),
        _ => None,
    }
}

fn has_refseq_mismatch(warnings: &[NormalizationWarning]) -> bool {
    warnings
        .iter()
        .any(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. }))
}

/// A genomic provider whose contig `NC_000001.11` carries `core` (the exact
/// reference tract) flanked by inert `N` padding, so a description anchored at
/// `CORE_START` spanning `core.len()` bases finds `core` verbatim under it.
fn g_provider(core: &str) -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000001.11", padded(core));
    p
}

/// A provider carrying `NM_000088.3` with `cds_start == 1` (so `c.N` maps to
/// transcript position `N`) and the reference tract `GTGTGCGC` at `c.100_107`
/// (transcript bases 100..=107 → 0-based indices 99..107). The tract is
/// surrounded by `A` padding within the CDS so the described 8 bp span lands
/// exactly on `GTGTGCGC`.
fn cds_gt2gc2_provider() -> MockProvider {
    const LEN: usize = 120;
    let mut bases = vec![b'A'; LEN];
    // c.100 with cds_start == 1 is transcript position 100 == 0-based index 99.
    for (i, b) in b"GTGTGCGC".iter().enumerate() {
        bases[99 + i] = *b;
    }
    let seq = String::from_utf8(bases).expect("ASCII bases are valid UTF-8");
    let tx = Transcript::new(
        "NM_000088.3".to_string(),
        Some("COL1A1".to_string()),
        Strand::Plus,
        seq,
        Some(1),          // cds_start: c.1 == transcript position 1
        Some(LEN as u64), // cds_end
        vec![Exon::new(1, 1, LEN as u64)],
        None,
        None,
        None,
        GenomeBuild::GRCh38,
        ManeStatus::Select,
        None,
        None,
    );
    let mut p = MockProvider::new();
    p.add_transcript(tx);
    p
}

#[test]
fn genome_multirepeat_is_a_normalization_passthrough() {
    // GT[2]GC[2] expands to GTGTGCGC (8 bp), so the described span is exactly
    // 8 bp wide and the reference under it literally spells GTGTGCGC.
    let end = CORE_START + 8 - 1;
    let input = format!("NC_000001.11:g.{}_{}GT[2]GC[2]", CORE_START, end);
    let normalizer = Normalizer::new(g_provider("GTGTGCGC"));
    let variant = parse_hgvs(&input).expect("parse genome MultiRepeat");
    assert!(
        matches!(edit_of(&variant), Some(NaEdit::MultiRepeat { .. })),
        "input must parse as a MultiRepeat so the test exercises the #953 arm",
    );
    let outcome = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize MultiRepeat");
    assert_eq!(
        format!("{}", outcome.result),
        input,
        "a compound MultiRepeat has no canonical 3'-shift; normalization must \
         return it unchanged (#953)",
    );
    // The reference bases match the declared tract, so the validated
    // pass-through must NOT emit a RefSeqMismatch. Absence here (over a real,
    // matching reference) confirms the arm ran and validated rather than
    // bailing early on a missing-reference error.
    assert!(
        !has_refseq_mismatch(&outcome.warnings),
        "matching reference tract must validate cleanly, got {:?}",
        outcome.warnings,
    );
}

#[test]
fn cds_multirepeat_is_a_normalization_passthrough() {
    // c.100_107 is an 8 bp span matching GT[2]GC[2] = GTGTGCGC.
    let input = "NM_000088.3:c.100_107GT[2]GC[2]";
    let normalizer = Normalizer::new(cds_gt2gc2_provider());
    let variant = parse_hgvs(input).expect("parse cds MultiRepeat");
    assert!(
        matches!(edit_of(&variant), Some(NaEdit::MultiRepeat { .. })),
        "input must parse as a MultiRepeat",
    );
    let outcome = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("normalize");
    assert_eq!(
        format!("{}", outcome.result),
        input,
        "cds MultiRepeat pass-through (#953)",
    );
    assert!(
        !has_refseq_mismatch(&outcome.warnings),
        "matching transcript tract must validate cleanly, got {:?}",
        outcome.warnings,
    );
}

#[test]
fn multirepeat_normalization_is_idempotent() {
    // AC[3]TG[2] expands to ACACACTGTG (10 bp); the span is 10 bp wide and the
    // reference under it spells ACACACTGTG.
    let end = CORE_START + 10 - 1;
    let input = format!("NC_000001.11:g.{}_{}AC[3]TG[2]", CORE_START, end);
    let normalizer = Normalizer::new(g_provider("ACACACTGTG"));
    let variant = parse_hgvs(&input).expect("parse");
    assert!(
        matches!(edit_of(&variant), Some(NaEdit::MultiRepeat { .. })),
        "input must parse as a MultiRepeat",
    );
    let once = normalizer.normalize(&variant).expect("first normalize");
    let twice = normalizer.normalize(&once).expect("second normalize");
    assert_eq!(
        format!("{once}"),
        input,
        "a validated MultiRepeat over a matching reference must return unchanged",
    );
    assert_eq!(
        format!("{once}"),
        format!("{twice}"),
        "MultiRepeat idempotent",
    );
}
