//! Issue #1209: an insertion whose 3'-shift can carry it past `cds_end` must
//! complete that shift on the **first** pass, instead of stopping on the
//! boundary and being rewritten to a `delins` by the issue #387 CDS-end clamp.
//!
//! ## The defect
//!
//! The #918 bound relaxation in `normalize_cds` lets the 3'-rule shuffle run
//! past `cds_end` into the 3'UTR, but originally applied only to
//! `del`/`dup`/`delins`. An `Insertion` therefore kept the `cds_end` bound,
//! stopped on the boundary, and hit the #387 clamp — whose output is a
//! `Delins`, which *does* get the relaxation. A second pass then shifted it the
//! rest of the way, so the two passes took different branches:
//!
//! ```text
//! c.25_26insGAT -> c.26delinsGATG -> c.*1_*2insTGA
//! ```
//!
//! ## Why the clamp is still needed
//!
//! Relaxing the bound does not disarm the #387 clamp. An insertion that
//! genuinely *saturates* `cds_end` — one with nowhere left to shift — is left
//! resting on the boundary by the shuffle, so the clamp fires exactly as
//! before. [`saturating_insertion_at_cds_end_is_still_clamped`] pins that, and
//! it is the case the clamp exists for. Only insertions that had somewhere to
//! go stop reaching it.
//!
//! ## Cross-axis oracle
//!
//! On a coding transcript `r.` *is* the CDS-relative axis — the same axis as
//! `c.` (HGVS `background/numbering.md` L58/L61, issue #469) — so the `r.`
//! spelling of an edit names the same bases over the byte-identical transcript
//! and must shift identically. `r.` already completed this shift in one pass,
//! which is what showed the bound (not the shuffle) was at fault, exactly as it
//! did for #1185. [`c_and_r_spellings_agree_after_the_fix`] pins the agreement.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

/// Build a coding transcript as `C * 40` (5'UTR) ++ `cds` ++ `utr3`, by hand
/// rather than via `SyntheticBuilder` so no property of the synthetic padding
/// can be blamed for the outcome.
///
/// **These fixtures are deliberately coordinate-only.** The `cds` here is not
/// an in-frame, stop-terminated ORF, and it is not meant to be: this module
/// exercises the *coordinate* behaviour of the `cds_end` axis transition
/// (which bound the shuffle respects, and when the #387 clamp fires), and
/// nothing on that path reads the reading frame or the stop codon — only
/// `cds_start` / `cds_end`. Nor could it be made biologically valid without
/// destroying the case: `cds_end` is by definition the last base of the stop
/// codon, so a stop-terminated CDS ends in `TAA`/`TAG`/`TGA`, and the defect
/// requires a CDS ending in a run that an insertion can rotate *through* into
/// the 3'UTR. `issue_387_canon_cds_end_clamp`'s own anchor fixture
/// (`NM_TEST387.1`, CDS `ATGAA`) is coordinate-only for the same reason.
///
/// Translation-level conformance is therefore explicitly out of scope here.
/// What keeps these assertions honest is not the fixture's realism but
/// [`c_and_r_spellings_agree_after_the_fix`], which re-derives the same answer
/// through a genuinely different code path — `normalize_rna` rather than
/// `normalize_cds` — over the byte-identical transcript.
fn coding_transcript(cds: &str, utr3: &str, strand: Strand) -> MockProvider {
    const UTR5_LEN: usize = 40;
    let sequence = format!("{}{cds}{utr3}", "C".repeat(UTR5_LEN));
    let cds_start = (UTR5_LEN + 1) as u64;
    let cds_end = (UTR5_LEN + cds.len()) as u64;
    let tx_len = sequence.len() as u64;

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_UTR3.1".to_string(),
        Some("UTR3_TEST".to_string()),
        strand,
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

/// A CDS ending in a `G` run, followed by an `ACGT…` 3'UTR.
///
/// The `GGGG` at `c.23_26` abuts `cds_end` and the UTR continues `A C G T`, so
/// an inserted `GAT` can rotate 3' *through* the boundary rather than
/// saturating at it:
///
/// ```text
/// …G G G |GAT| G  A  C…      c.25_26insGAT
/// …G G G  G  |ATG| A  C…     c.26_*1insATG
/// …G G G  G   A  |TGA| C…    c.*1_*2insTGA   <- stops here (*2 is C)
/// ```
fn transcript_with_g_run_at_cds_end(strand: Strand) -> MockProvider {
    const CDS: &str = "CCCCATATTTTTCACACACACAGGGG";
    // The run must be the CDS *suffix*, not merely present at index 22 — the
    // whole point is that it abuts `cds_end`.
    assert_eq!(CDS.len(), 26, "c.26 must be the last CDS base");
    assert_eq!(&CDS[22..], "GGGG", "c.23_26 must be the G run at cds_end");
    coding_transcript(CDS, &"ACGT".repeat(10), strand)
}

/// Normalize `input` twice and return `(once, twice)` as rendered strings.
fn normalize_twice(provider: MockProvider, input: &str) -> (String, String) {
    let normalizer = Normalizer::new(provider);
    let variant: HgvsVariant = parse_hgvs(input).expect("parse failed");
    let once = normalizer
        .normalize(&variant)
        .expect("first normalize failed");
    let twice = normalizer.normalize(&once).expect("re-normalize failed");
    (once.to_string(), twice.to_string())
}

/// The reported case: the shift completes on the first pass and is a fixed
/// point, rather than stopping at `cds_end` as `c.26delinsGATG`.
#[test]
fn insertion_shifts_past_cds_end_in_one_pass() {
    let (once, twice) = normalize_twice(
        transcript_with_g_run_at_cds_end(Strand::Plus),
        "NM_UTR3.1:c.25_26insGAT",
    );
    assert_eq!(
        once, "NM_UTR3.1:c.*1_*2insTGA",
        "pass 1 must complete the shift into the 3'UTR"
    );
    assert_eq!(twice, once, "and it is a fixed point");
}

/// The intermediate the old code stopped at normalizes to the same answer, so
/// the two spellings are confluent rather than each its own fixed point.
#[test]
fn the_old_clamped_output_converges_on_the_shifted_form() {
    let (once, twice) = normalize_twice(
        transcript_with_g_run_at_cds_end(Strand::Plus),
        "NM_UTR3.1:c.26delinsGATG",
    );
    assert_eq!(once, "NM_UTR3.1:c.*1_*2insTGA");
    assert_eq!(twice, once, "and is stable from there");
}

/// Cross-axis oracle: `r.` on a coding transcript is the same axis as `c.`
/// (#469), so the `r.` spelling over the byte-identical transcript must give
/// the same answer. It already did before the fix; now `c.` agrees with it.
#[test]
fn c_and_r_spellings_agree_after_the_fix() {
    let (r_once, r_twice) = normalize_twice(
        transcript_with_g_run_at_cds_end(Strand::Plus),
        "NM_UTR3.1:r.25_26insgau",
    );
    assert_eq!(r_once, "NM_UTR3.1:r.*1_*2insuga");
    assert_eq!(r_twice, r_once, "the r. axis is idempotent");

    let (c_once, _) = normalize_twice(
        transcript_with_g_run_at_cds_end(Strand::Plus),
        "NM_UTR3.1:c.25_26insGAT",
    );
    // Same positions, same rotated payload, RNA alphabet vs DNA.
    assert_eq!(
        c_once.replace("c.", "").replace("insTGA", "insuga"),
        r_once.replace("r.", ""),
        "the two spellings of one axis must name the same edit"
    );
}

/// Present on the minus strand too, so the fix is not strand-specific.
#[test]
fn insertion_shifts_past_cds_end_in_one_pass_on_minus_strand() {
    let (once, twice) = normalize_twice(
        transcript_with_g_run_at_cds_end(Strand::Minus),
        "NM_UTR3.1:c.25_26insGAT",
    );
    assert_eq!(once, "NM_UTR3.1:c.*1_*2insTGA");
    assert_eq!(twice, once, "and is a fixed point");
}

// ---------------------------------------------------------------------------
// Guards: the relaxation must not disturb what the #387 clamp is actually for.
// ---------------------------------------------------------------------------

/// **The clamp is still armed.** This mirrors #387's anchor case
/// (`issue_387_canon_cds_end_clamp::three_prime_insertion_at_cds_end_clamps_to_delins_at_cds_end`,
/// itself modelled on `NM_212556.2:c.1400_1401insAC`) on a hand-built
/// transcript: CDS `ATGAA`, 3'UTR all `G`.
///
/// Inserting `AC` between `c.4` and `c.5` rotates once past the `A` at `c.5`
/// to `c.5_*1insCA`, and then stops — `*1` is `G` and the payload starts `C`.
/// The insertion therefore *saturates* `cds_end`, the shuffle leaves it resting
/// there, and the #387 clamp rewrites it to `c.5delinsACA`.
///
/// This is the case the clamp exists for. Relaxing the bound for insertions
/// (#1209) does not disturb it, because an insertion with nowhere left to go
/// stops on the boundary whether or not the bound allows it to pass.
#[test]
fn saturating_insertion_at_cds_end_is_still_clamped() {
    let provider = coding_transcript("ATGAA", &"G".repeat(40), Strand::Plus);
    let (once, twice) = normalize_twice(provider, "NM_UTR3.1:c.4_5insAC");
    assert_eq!(
        once, "NM_UTR3.1:c.5delinsACA",
        "an insertion with nowhere to shift still reaches the #387 clamp"
    );
    assert_eq!(twice, once, "and the clamped form is a fixed point");
}

/// An insertion whose shift stays **inside** the CDS is untouched by the
/// relaxation, so the change is scoped to boundary-crossing shifts.
#[test]
fn insertion_shifting_within_the_cds_is_unaffected() {
    let provider = coding_transcript("GGAAAAAACCCCGGGG", &"ACGT".repeat(10), Strand::Plus);
    let (once, twice) = normalize_twice(provider, "NM_UTR3.1:c.3_4insA");
    assert_eq!(
        once, "NM_UTR3.1:c.8dup",
        "the A-run shift lands inside the CDS and is spelled as a duplication"
    );
    assert_eq!(twice, once, "and is a fixed point");
}
