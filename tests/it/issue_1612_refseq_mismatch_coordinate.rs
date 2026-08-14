//! Issue #1612 — `REFSEQ_MISMATCH` must not report a window-relative offset as
//! if it were a coordinate.
//!
//! `Normalizer::normalize_na_edit` indexes the slice it is handed, so its
//! `start`/`end` are in the frame of that slice. Most callers hand it the whole
//! sequence the description addresses, and for those the two coincide. The
//! genomic and mitochondrial callers hand it a **fetched window**, and nothing
//! converted the pair back before it was formatted into the warning's
//! `position` field. Measured on `main` before the fix, against a 400-base
//! `ACGT`-cycling contig whose base 201 is `A`:
//!
//! ```text
//! NC_000001.11:g.201C>T   ->  reference sequence mismatch at 100-100
//! ```
//!
//! `100` is `201 - window_start`, with `window_start = 201 - window_size`. The
//! mismatch is real; only the coordinate was wrong.
//!
//! # What the field promises
//!
//! **The 1-based inclusive span on the reference sequence that was actually
//! read**, in that sequence's own numbering. That is what the cis producer
//! (`merge::authored_reference_mismatch`) already documented and implemented —
//! "the sequence-axis span, via `region_sequence_delta`" — and the two must
//! agree, because the duplicate-suppression in `normalize_core` is keyed on
//! this string and nothing else.
//!
//! So the promise is deliberately **not** "the description's own axis". A
//! `c.5C>T` on a transcript whose CDS starts at 21 reports the transcript
//! offset `25-25` on both producers; moving it to `5-5` would have been a
//! second wrong answer and would have broken the dedupe. Those axes were
//! already correct and are pinned below so a future change cannot quietly
//! "fix" them apart.
//!
//! # Three things this actually repairs
//!
//! 1. The genomic/mitochondrial coordinate itself.
//! 2. The `location` of strict mode's `FerroError::ReferenceMismatch`, which is
//!    the same string.
//! 3. The dedupe against the cis producer, which had silently stopped matching
//!    on the genomic axis — a two-member allele with two distinct mismatches
//!    reported **four** warnings, two of them at the *same* fabricated
//!    `100-100`, because each member is windowed about itself.
//!
//! The intronic and boundary-spanning transcript paths are a fourth case and a
//! different answer: they shuffle against the genomic **contig**, which is not
//! the sequence the accession names, so a bare number there would be a
//! coordinate in a space the reader cannot guess. Those report the genomic span
//! qualified by the contig accession. Nothing dedupes against them —
//! `merge::simple_cds_pos` declines every offset position, so no cis-side
//! duplicate of an intronic finding can exist.

use std::io::Write;

use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, FerroError, JsonProvider, MockProvider, NormalizeConfig, Normalizer};

/// The fixture alphabet, shared by every provider below: base `n` (1-based) is
/// `"ACGT"[(n - 1) % 4]`.
fn base_at(pos1: usize) -> char {
    ['A', 'C', 'G', 'T'][(pos1 - 1) % 4]
}

/// A genome-capable provider serving one `len`-base `ACGT`-cycling contig.
fn genomic_provider(contig: &str, len: usize) -> JsonProvider {
    let seq: String = (1..=len).map(base_at).collect();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// Every `RefSeqMismatch` position in `warnings`, in order.
fn mismatch_positions(warnings: &[NormalizationWarning]) -> Vec<String> {
    warnings
        .iter()
        .filter_map(|w| match w {
            NormalizationWarning::RefSeqMismatch { position, .. } => Some(position.clone()),
            _ => None,
        })
        .collect()
}

/// The single `RefSeqMismatch` in `warnings`, or a panic naming what was there.
fn sole_mismatch(warnings: &[NormalizationWarning]) -> &NormalizationWarning {
    let mut found = warnings
        .iter()
        .filter(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. }));
    let first = found
        .next()
        .unwrap_or_else(|| panic!("expected a RefSeqMismatch; got {warnings:?}"));
    assert!(
        found.next().is_none(),
        "expected exactly one RefSeqMismatch; got {warnings:?}"
    );
    first
}

/// Normalize leniently and hand back the warnings.
fn lenient_warnings<P>(provider: P, input: &str) -> Vec<NormalizationWarning>
where
    P: ferro_hgvs::reference::ReferenceProvider,
{
    let variant = parse_hgvs(input).unwrap();
    Normalizer::with_config(provider, NormalizeConfig::lenient())
        .normalize_with_diagnostics(&variant)
        .unwrap_or_else(|e| panic!("lenient must accept `{input}`: {e:?}"))
        .warnings
}

/// The bases the inputs below assume. Asserted rather than trusted, so a change
/// to `base_at` fails here — where the cause is obvious — rather than
/// downstream as a mystery about coordinates.
#[test]
fn the_fixture_carries_the_bases_the_inputs_assume() {
    assert_eq!(base_at(201), 'A', "the g. reproducer states C over an A");
    assert_eq!(base_at(300), 'T', "the second cis member states G over a T");
    assert_eq!(base_at(314), 'C', "the plus-strand intronic base");
    assert_eq!(
        base_at(395),
        'G',
        "the minus-strand intronic base, genomic view"
    );
}

// ---------------------------------------------------------------------------
// The reported issue
// ---------------------------------------------------------------------------

/// The issue's own reproducer. `g.201` is `A`; stating `C` is a real mismatch,
/// and it must be reported at `201-201` rather than at the window offset
/// `100-100`.
#[test]
fn a_genomic_mismatch_reports_the_genomic_coordinate() {
    let warnings = lenient_warnings(
        genomic_provider("NC_000001.11", 400),
        "NC_000001.11:g.201C>T",
    );
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        stated_ref,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!("sole_mismatch only returns this variant")
    };

    assert_eq!(position, "201-201", "warning: {warning}");
    assert_eq!((stated_ref.as_str(), actual_ref.as_str()), ("C", "A"));
    // The rendered message is the user-visible surface, and it is what the CLI
    // prints, so pin it too rather than only the field behind it.
    assert!(
        warning
            .to_string()
            .starts_with(r#"reference sequence mismatch at 201-201: stated "C", actual "A""#),
        "rendered message: {warning}"
    );
}

/// The same string is strict mode's `FerroError::ReferenceMismatch` `location`,
/// so a caller parsing it to point at the offending base is pointed at `g.201`.
#[test]
fn strict_mode_locates_the_error_at_the_same_coordinate() {
    let variant = parse_hgvs("NC_000001.11:g.201C>T").unwrap();
    let err = Normalizer::with_config(
        genomic_provider("NC_000001.11", 400),
        NormalizeConfig::strict(),
    )
    .normalize(&variant)
    .expect_err("strict must reject a wrong-ref substitution");

    match err {
        FerroError::ReferenceMismatch {
            location,
            expected,
            found,
        } => {
            assert_eq!(location, "201-201");
            assert_eq!((expected.as_str(), found.as_str()), ("C", "A"));
        }
        other => panic!("expected ReferenceMismatch, got {other:?}"),
    }
}

/// The mitochondrial axis reaches the same windowed helper
/// (`normalize_in_grown_window`) and carried the same defect.
#[test]
fn a_mitochondrial_mismatch_reports_the_mitochondrial_coordinate() {
    let warnings = lenient_warnings(genomic_provider("NC_012920.1", 400), "NC_012920.1:m.201C>T");
    assert_eq!(mismatch_positions(&warnings), vec!["201-201".to_string()]);
}

// ---------------------------------------------------------------------------
// The dedupe the wrong coordinate had broken
// ---------------------------------------------------------------------------

/// Two cis members, each misstating its own reference base, must produce
/// **two** warnings — one per member, at the member's own coordinate.
///
/// Before the fix this produced **four**: the cis producer's `201-201` and
/// `300-300`, plus two per-member warnings *both* reading `100-100`, because
/// each member is windowed about itself and so lands at the same relative
/// offset. `normalize_core`'s suppression is keyed on `position` alone, so
/// neither of the real findings matched the fabricated ones and every mismatch
/// was reported twice — while the two fabricated coordinates were, on their
/// face, indistinguishable from each other.
#[test]
fn each_genomic_cis_member_is_reported_exactly_once() {
    let warnings = lenient_warnings(
        genomic_provider("NC_000001.11", 400),
        "NC_000001.11:g.[201C>T;300G>A]",
    );
    let mut positions = mismatch_positions(&warnings);
    positions.sort();
    assert_eq!(
        positions,
        vec!["201-201".to_string(), "300-300".to_string()],
        "one warning per misstating member, at its own coordinate; got {warnings:?}"
    );
}

// ---------------------------------------------------------------------------
// The axes that were already right, and must stay where they are
// ---------------------------------------------------------------------------

/// A coding transcript whose CDS starts at 21, so the `c.` axis and the
/// transcript-sequence axis differ by 20 and a confusion between them is
/// visible.
fn coding_transcript_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let seq: String = (1..=300).map(base_at).collect();
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        Some(seq),
        Some(21),
        Some(200),
        vec![Exon::new(1, 1, 300)],
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

/// `c.`, `n.` and `r.` all address the transcript sequence directly, so their
/// positions were never window-relative and must not move.
///
/// `c.5` is transcript offset 25 (`cds_start = 21`), whose base is `A`; stating
/// `C` mismatches. The reported span is `25-25` — the **sequence-axis** span,
/// which is what the cis producer reports for the same member and therefore
/// what the dedupe requires. Re-spelling it as `5-5` would look like a fix and
/// would break that agreement, so it is pinned here as a deliberate answer
/// rather than left to be re-derived.
#[test]
fn transcript_axes_report_the_sequence_axis_span_unchanged() {
    for (input, expected) in [
        ("NM_TEST.1:c.5C>T", "25-25"),
        ("NM_TEST.1:n.25C>T", "25-25"),
        ("NM_TEST.1:r.5c>u", "25-25"),
    ] {
        let warnings = lenient_warnings(coding_transcript_provider(), input);
        assert_eq!(
            mismatch_positions(&warnings),
            vec![expected.to_string()],
            "`{input}` must report the transcript-sequence span; got {warnings:?}"
        );
    }
}

/// The agreement itself, rather than the two halves separately: a `c.` cis
/// allele whose members both misstate is reported once per member, which can
/// only happen if the per-member and cis producers spell the position
/// identically.
#[test]
fn each_coding_cis_member_is_reported_exactly_once() {
    let warnings = lenient_warnings(coding_transcript_provider(), "NM_TEST.1:c.[5C>T;7A>G]");
    let mut positions = mismatch_positions(&warnings);
    positions.sort();
    assert_eq!(
        positions,
        vec!["25-25".to_string(), "27-27".to_string()],
        "got {warnings:?}"
    );
}

// ---------------------------------------------------------------------------
// The case where the sequence read is not the one the accession names
// ---------------------------------------------------------------------------

/// A two-exon transcript with genomic backing, on `strand`. Exon 1 covers
/// transcript 1..10 and exon 2 covers 11..22; the genomic placement is chosen
/// so `c.10+5` lands on a base whose coordinate is checkable by hand.
fn intronic_provider(
    accession: &str,
    contig: &str,
    strand: Strand,
    exon1: (u64, u64),
    exon2: (u64, u64),
) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(contig, (1..=1000).map(base_at).collect::<String>());
    let lo = exon1.0.min(exon2.0);
    let hi = exon1.1.max(exon2.1);
    provider.add_transcript(Transcript::new(
        accession.to_string(),
        Some("INTR".to_string()),
        strand,
        Some("AAAATGCCCCGGGGTAGAATAA".to_string()),
        Some(1),
        Some(22),
        vec![
            Exon::with_genomic(1, 1, 10, exon1.0, exon1.1),
            Exon::with_genomic(2, 11, 22, exon2.0, exon2.1),
        ],
        Some(contig.to_string()),
        Some(lo),
        Some(hi),
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// An intronic description is shuffled against the genomic contig, so the span
/// it reports is genomic — and says so, because the accession the description
/// names is a transcript and a bare number would be unreadable.
///
/// Plus strand: exon 1 is transcript 1..10 at genomic 300..309, so transcript
/// position 10 is genomic 309 and `c.10+5` is genomic **314**. Stating `A`
/// there misstates a `C`. Before the fix this read `101-101`, the offset into a
/// window opening at genomic 214.
#[test]
fn a_plus_strand_intronic_mismatch_names_the_contig_span() {
    let provider = intronic_provider("NM_INTR.1", "chr_i", Strand::Plus, (300, 309), (400, 411));
    let warnings = lenient_warnings(provider.clone(), "NM_INTR.1:c.10+5delA");
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };
    assert_eq!(position, "chr_i:314-314", "warning: {warning}");
    // Independent of the coordinate arithmetic above: whatever base the warning
    // reports as actual must be the base standing at the coordinate it reports.
    assert_eq!(actual_ref.as_str(), "C");
    assert_eq!(base_at(314), 'C');

    // A two-base span keeps both endpoints, and they stay in genomic order.
    let warnings = lenient_warnings(provider, "NM_INTR.1:c.10+5_10+6delAA");
    assert_eq!(
        mismatch_positions(&warnings),
        vec!["chr_i:314-315".to_string()]
    );
}

/// Minus strand, where the window is additionally reverse-complemented before
/// the shuffle sees it — so the relative offset is not even in the same
/// direction as the contig.
///
/// Exon 1 is transcript 1..10 at genomic 400..409 read backwards, so transcript
/// position 10 is genomic 400 and `c.10+5` — five bases further along the
/// transcript — is genomic **395**. The contig's base there is `G`; the
/// transcript reads its reverse complement, `C`, which is what the warning must
/// report as actual. That pairing is the independent check on the coordinate:
/// only the right genomic position revcomps to the reported base.
#[test]
fn a_minus_strand_intronic_mismatch_names_the_contig_span() {
    let provider = intronic_provider("NM_MINUS.1", "chr_m", Strand::Minus, (400, 409), (300, 311));
    let warnings = lenient_warnings(provider.clone(), "NM_MINUS.1:c.10+5delA");
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };
    assert_eq!(position, "chr_m:395-395", "warning: {warning}");
    assert_eq!(base_at(395), 'G', "the contig's own base");
    assert_eq!(
        actual_ref.as_str(),
        "C",
        "the transcript reads the reverse complement of the contig base at the \
         coordinate the warning names"
    );

    // Still low-to-high in genomic order, which is the opposite end of the
    // transcript's `10+5_10+6`.
    let warnings = lenient_warnings(provider, "NM_MINUS.1:c.10+5_10+6delAA");
    assert_eq!(
        mismatch_positions(&warnings),
        vec!["chr_m:394-395".to_string()]
    );
}

/// The `n.` mirror of the plus/minus intronic pair above.
///
/// `normalize_intronic_tx` builds its own `MismatchFrame::GenomicWindow` —
/// the `c.` path's does not stand in for it, and a wrong `seq_start` in either
/// one is invisible to the other. `n.10+5` is the same locus as `c.10+5` here
/// (`cds_start = 1`), so the coordinates below are the ones the `c.` tests
/// already state, arrived at through a different function.
#[test]
fn an_intronic_noncoding_mismatch_names_the_contig_span() {
    let plus = intronic_provider("NM_INTR.1", "chr_i", Strand::Plus, (300, 309), (400, 411));
    let warnings = lenient_warnings(plus.clone(), "NM_INTR.1:n.10+5delA");
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };
    assert_eq!(position, "chr_i:314-314", "warning: {warning}");
    // The same independent pairing the `c.` tests make: the base reported as
    // actual must be the base standing at the coordinate reported.
    assert_eq!(actual_ref.as_str(), "C");
    assert_eq!(base_at(314), 'C');

    let warnings = lenient_warnings(plus, "NM_INTR.1:n.10+5_10+6delAA");
    assert_eq!(
        mismatch_positions(&warnings),
        vec!["chr_i:314-315".to_string()]
    );

    // Minus strand, where the window is reverse-complemented before the shuffle
    // sees it, so a `seq_start` error and a flip error are separable.
    let minus = intronic_provider("NM_MINUS.1", "chr_m", Strand::Minus, (400, 409), (300, 311));
    let warnings = lenient_warnings(minus, "NM_MINUS.1:n.10+5delA");
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };
    assert_eq!(position, "chr_m:395-395", "warning: {warning}");
    assert_eq!(base_at(395), 'G', "the contig's own base");
    assert_eq!(
        actual_ref.as_str(),
        "C",
        "the transcript reads the reverse complement of the contig base at the \
         coordinate the warning names"
    );
}

/// The reverse complement of `bases`, so a minus-strand expectation is derived
/// from the contig rather than transcribed by hand.
fn revcomp(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// The contig's own bases over the 1-based inclusive span `[lo, hi]`.
fn contig_bases(lo: usize, hi: usize) -> String {
    (lo..=hi).map(base_at).collect()
}

/// A span with one **exonic** and one intronic endpoint does not reach the
/// intronic path at all — it is `normalize_boundary_spanning_cds`, a third
/// `GenomicWindow` construction site with its own `seq_start`.
///
/// `c.9` is genomic 308 and `c.10+2` is genomic 311 on the plus-strand fixture,
/// so the four bases actually read are `g.308_311`. The transcript is CDS-only
/// (`cds_start = 1`), which is what lets the `n.` mirror below name the same
/// locus and so compare like with like.
#[test]
fn a_boundary_spanning_coding_mismatch_names_the_contig_span() {
    let plus = intronic_provider("NM_INTR.1", "chr_i", Strand::Plus, (300, 309), (400, 411));
    let warnings = lenient_warnings(plus, "NM_INTR.1:c.9_10+2delAAAA");
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };
    assert_eq!(position, "chr_i:308-311", "warning: {warning}");
    assert_eq!(
        actual_ref.as_str(),
        contig_bases(308, 311),
        "the bases reported as actual must be the contig's own bases over the \
         span the warning names",
    );

    // Minus strand: exon 1 is transcript 1..10 at genomic 400..409 read
    // backwards, so `c.9` is genomic 401 and `c.10+2` is genomic 398.
    let minus = intronic_provider("NM_MINUS.1", "chr_m", Strand::Minus, (400, 409), (300, 311));
    let warnings = lenient_warnings(minus, "NM_MINUS.1:c.9_10+2delAAAA");
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };
    assert_eq!(position, "chr_m:398-401", "warning: {warning}");
    assert_eq!(
        actual_ref.as_str(),
        revcomp(&contig_bases(398, 401)),
        "the transcript reads the reverse complement of the contig bases over \
         the span the warning names",
    );
}

/// The `n.` mirror of the test above — `normalize_boundary_spanning_tx`, the
/// fourth and last `GenomicWindow` site. Same locus, same expected span,
/// different function.
#[test]
fn a_boundary_spanning_noncoding_mismatch_names_the_contig_span() {
    let plus = intronic_provider("NM_INTR.1", "chr_i", Strand::Plus, (300, 309), (400, 411));
    let warnings = lenient_warnings(plus, "NM_INTR.1:n.9_10+2delAAAA");
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };
    assert_eq!(position, "chr_i:308-311", "warning: {warning}");
    assert_eq!(actual_ref.as_str(), contig_bases(308, 311));

    let minus = intronic_provider("NM_MINUS.1", "chr_m", Strand::Minus, (400, 409), (300, 311));
    let warnings = lenient_warnings(minus, "NM_MINUS.1:n.9_10+2delAAAA");
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };
    assert_eq!(position, "chr_m:398-401", "warning: {warning}");
    assert_eq!(actual_ref.as_str(), revcomp(&contig_bases(398, 401)));
}

// ---------------------------------------------------------------------------
// The second, independent window offset: the growth cap
// ---------------------------------------------------------------------------

/// Flanking bases either side of the tract in [`capped_growth_provider`]. A
/// multiple of 4 so the cyclic `ACGT` lead ends on `T` rather than on the
/// tract's own base, which would extend the run.
const CAPPED_FLANK: usize = 64;

/// Long enough that the geometric growth in `normalize_in_grown_window`
/// exhausts `MAX_SHUFFLE_FETCH_WINDOW` (64 KiB) with contig left to read on
/// both sides of the positions below — which is the only way to reach
/// `canonicalize_without_shifting`. Modelled on
/// `tests/it/issue_1691_homopolymer_convergence.rs`, which uses 70_000 for the
/// same purpose; this is larger so that a position can sit a full window past
/// the tract's 5' end and still be a window short of its 3' end.
const CAPPED_TRACT: usize = 200_000;

/// One contig: `CAPPED_FLANK` bases of cyclic `ACGT`, then `CAPPED_TRACT` `A`s,
/// then `CAPPED_FLANK` bases of cyclic `CGTA`. Base `CAPPED_FLANK + 1` is the
/// first `A` of the tract, so stating any other base there is a mismatch.
fn capped_growth_provider(contig: &str) -> JsonProvider {
    let lead: String = "ACGT".chars().cycle().take(CAPPED_FLANK).collect();
    let tract: String = std::iter::repeat_n('A', CAPPED_TRACT).collect();
    let trail: String = "CGTA".chars().cycle().take(CAPPED_FLANK).collect();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: format!("{lead}{tract}{trail}") },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// Past the growth cap the warning is built by
/// `canonicalize_without_shifting`, which normalizes against a **slice** of the
/// window it was handed and so carries a *second* frame offset
/// (`slice_start = window_start + lo`). It is not a refactor of the growth
/// loop's own `MismatchFrame::Window`: it adds `lo`, so a correct growth loop
/// says nothing about it, and a wrong `lo` reproduces #1612 exactly — a
/// plausible number in a user-facing coordinate with nothing red.
///
/// Two positions, because the two terms of that sum are separable:
///
/// - `g.65`, one base into the tract, is closer to the contig start than the
///   64 KiB cap, so the window opens at base 1 and `window_start` is **0**.
///   Only `lo` can move the answer.
/// - `g.100000` is a full cap-width into the contig, so `window_start` is
///   non-zero **and** `lo` is added on top of it. A frame that carried either
///   term alone lands somewhere else.
///
/// Both must come back unshifted — that is the capped answer's signature, and
/// asserting it is what pins these two inputs to this path rather than to the
/// growth loop's own site. A position the growth loop *can* settle (`g.150000`
/// below) shifts to the tract's 3' end instead.
#[test]
fn a_capped_growth_mismatch_reports_the_contig_coordinate() {
    let provider = capped_growth_provider("c");
    let first = CAPPED_FLANK + 1;

    for authored in [first, 100_000] {
        let input = format!("c:g.{authored}delC");
        let variant = parse_hgvs(&input).unwrap();
        let outcome = Normalizer::with_config(provider.clone(), NormalizeConfig::lenient())
            .normalize_with_diagnostics(&variant)
            .unwrap_or_else(|e| panic!("lenient must accept `{input}`: {e:?}"));

        assert_eq!(
            outcome.result.to_string(),
            format!("c:g.{authored}del"),
            "`{input}` must come back unshifted; a shifted answer means the \
             growth cap was never reached and this test is not on the path it \
             claims to cover",
        );
        assert_eq!(
            mismatch_positions(&outcome.warnings),
            vec![format!("{authored}-{authored}")],
            "`{input}` must report the contig coordinate; got {:?}",
            outcome.warnings,
        );
    }

    // The control: far enough from the tract's 3' end that the grown window
    // contains the shuffle's resting place, so this one never reaches the cap
    // and its warning comes from the growth loop's own frame. It shifts.
    let input = "c:g.150000delG";
    let warnings = lenient_warnings(provider, input);
    assert_eq!(
        mismatch_positions(&warnings),
        vec!["150000-150000".to_string()],
        "{input}",
    );
}

// ---------------------------------------------------------------------------
// The same window offsets, reached through the validator's own message text
// ---------------------------------------------------------------------------

/// `validate_stated_length` built the *same* window offsets into the warning's
/// `actual_ref`, in HGVS range notation. With `position` repaired the two
/// halves of one warning contradicted each other — measured on
/// `g.201_202del5`:
///
/// ```text
/// at 201-202: stated "length 5", actual "span 2 (100_101)"
/// ```
///
/// The endpoints are `position`'s job and it now states them correctly, so the
/// restatement is dropped rather than duplicated. The width, which is what this
/// check is actually about, stays.
#[test]
fn a_numeric_length_mismatch_states_the_width_and_no_window_offset() {
    let warnings = lenient_warnings(
        genomic_provider("NC_000001.11", 400),
        "NC_000001.11:g.201_202del5",
    );
    let warning = sole_mismatch(&warnings);
    let NormalizationWarning::RefSeqMismatch {
        position,
        stated_ref,
        actual_ref,
        ..
    } = warning
    else {
        unreachable!()
    };

    assert_eq!(position, "201-202");
    assert_eq!(
        (stated_ref.as_str(), actual_ref.as_str()),
        ("length 5", "span 2")
    );
    let rendered = warning.to_string();
    assert!(
        !rendered.contains("100"),
        "no window offset may survive anywhere in the message: {rendered}"
    );
}
