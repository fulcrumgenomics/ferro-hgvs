//! Regression pins for `c.`-axis 3'-shifts that cross `cds_end` into the 3'UTR
//! (`*N`) — a family that needed **two** normalization passes to settle and now
//! settles in one.
//!
//! Found by the extended `normalize_idempotency_proptest` generator (the
//! `n.`/`r.` axis + `inv`/`con` edit-kind extension). The `dup`/`delins`/`con`
//! finding was filed as #1185 and fixed by PR #1189 (complete the shift in one
//! pass) together with the #1192 codon-frame gate. Every case below is now
//! idempotent on the first pass and is pinned here as a regression over a
//! **hand-built** transcript, complementing the random synthetic fuzz that
//! found them.
//!
//! The `ins` sibling — an insertion resting at `cds_end` that the #387 clamp
//! rewrote instead of letting it shift on — was #1209, fixed by PR #1211. It is
//! pinned in its own module, `issue_1209_cds_end_insertion_shift`, rather than
//! duplicated here. With both fixed, the `c.`-with-3'UTR systems are fuzzed
//! directly as `Sys::CdsUtr3Plus` / `Sys::CdsUtr3Minus` in
//! `tests/it/normalize_idempotency_proptest.rs`, so this module is now a
//! deterministic backstop rather than a stand-in for missing fuzz coverage.
//!
//! ## Why these cases were real bugs and not fixture artifacts
//!
//! [`rna_axis_shifts_across_cds_end_in_one_pass`] is the control. It runs the
//! **same edit over the byte-identical transcript**, differing only in the axis
//! it is spelled on (`r.` instead of `c.`). On a coding transcript `r.` *is*
//! the CDS-relative axis — the same axis as `c.` (HGVS
//! `background/numbering.md` L58/L61, issue #469, pinned by
//! `issue_291_rna_axis_convention`) — so both spellings name the same bases and
//! must shift identically. `r.` settled in one pass while `c.` did not, which
//! localised the defect to the `c.` code path rather than to the reference
//! bases or to 3'-shifting across `cds_end` in general. Both spellings now
//! agree, and that agreement is what these tests assert.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

/// A plain coding transcript with a 5'UTR, a short CDS, and a 3'UTR — built by
/// hand rather than via `SyntheticBuilder` so no property of the synthetic
/// padding can be blamed for the outcome.
///
/// | tx range   | bases                                      | region        |
/// |------------|--------------------------------------------|---------------|
/// | 1..=40     | `C` * 40                                   | 5'UTR         |
/// | 41..=49    | `GGAAAAAAC`                                | CDS c.1..c.9  |
/// | 50..=89    | `AC` * 20                                  | 3'UTR `*1..`  |
///
/// The CDS *ends* with `AC` (`c.8_9`) and the 3'UTR *begins* with `AC`, so a
/// deletion of that `AC` has an `AC` tandem to 3'-shift along — a shift that
/// necessarily crosses `cds_end` and must be rendered with `*N` positions.
fn transcript_with_utr3_ac_tandem(strand: Strand) -> MockProvider {
    let provider = coding_transcript("GGAAAAAAC", &"AC".repeat(20), strand);
    // c.8_9 is the `AC` at the end of the CDS, and the 3'UTR continues it.
    assert_eq!(&"GGAAAAAAC"[7..9], "AC", "c.8_9 must be 'AC'");
    provider
}

/// Same shape, but the 3'UTR is a period-4 `ACGT` tandem, so an `AC` repeat
/// running out of the CDS *terminates* two bases into the UTR (`*3` is `G`)
/// instead of continuing to the transcript end. That bounded-tract layout is
/// what exposes the `dup` symptom.
fn transcript_with_utr3_acgt_padding(strand: Strand) -> MockProvider {
    coding_transcript("GAAAAAACAC", &"ACGT".repeat(10), strand)
}

/// Build a coding transcript as `C * 40` (5'UTR) ++ `cds` ++ `utr3`.
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

/// A `delins` that reduces to a deletion whose 3'-shift crosses `cds_end`
/// completes that shift on the **first** pass.
///
/// `c.7_9delinsA`: ref `c.7_9` = `AAC`, alt `A`; the shared `A` prefix trims it
/// to a deletion of `AC` at `c.8_9`, which 3'-shifts along the `AC` tandem
/// running from `c.8` into the 3'UTR. The one-pass answer is the fully shifted
/// `c.*39_*40del`. Before #1189 pass 1 stopped at the unshifted `c.8_9del`.
#[test]
fn delins_reducing_to_deletion_shifts_across_cds_end_in_one_pass() {
    let (once, twice) = normalize_twice(
        transcript_with_utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:c.7_9delinsA",
    );
    assert_eq!(
        once, "NM_UTR3.1:c.*39_*40del",
        "pass 1 shifts all the way past cds_end"
    );
    assert_eq!(twice, once, "and is a fixed point");
}

/// The intermediate spelling pass 1 used to stop at is itself normalized to the
/// same fully shifted form — so the two spellings are confluent, not merely
/// each idempotent.
#[test]
fn unshifted_deletion_at_cds_end_shifts_to_the_same_fixed_point() {
    let (once, twice) = normalize_twice(
        transcript_with_utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:c.8_9del",
    );
    assert_eq!(once, "NM_UTR3.1:c.*39_*40del");
    assert_eq!(
        twice, "NM_UTR3.1:c.*39_*40del",
        "and it is stable from there"
    );
}

/// A `dup` whose 3'-shift crosses `cds_end` used to diverge across passes in
/// *edit kind*, not merely position: `dup` -> `ins` -> repeat. It now reaches
/// the repeat spelling directly.
///
/// CDS `GAAAAAACAC` (`c.1..c.10`) followed by an `ACGT…` 3'UTR puts an `AC[3]`
/// tandem at `c.7..*2`, bounded by the `G` at `*3`. Duplicating `c.7_10`
/// (`ACAC`) 3'-shifts within that tandem, and the shift crosses `cds_end`.
#[test]
fn duplication_crossing_cds_end_settles_in_one_pass() {
    let (once, twice) = normalize_twice(
        transcript_with_utr3_acgt_padding(Strand::Plus),
        "NM_UTR3.1:c.7_10dup",
    );
    assert_eq!(
        once, "NM_UTR3.1:c.7_*2AC[5]",
        "pass 1 spells the shifted duplication as the repeat directly, rather \
         than parking an insertion in the 3'UTR for a second pass to re-spell"
    );
    assert_eq!(twice, once, "and is a fixed point");

    // Cross-axis oracle: `r.` reaches the same answer through `normalize_rna`
    // rather than `normalize_cds`, so this is corroboration by a different code
    // path over the byte-identical transcript — not another self-consistency
    // check on the same one.
    let (r_once, r_twice) = normalize_twice(
        transcript_with_utr3_acgt_padding(Strand::Plus),
        "NM_UTR3.1:r.7_10dup",
    );
    assert_eq!(r_once, "NM_UTR3.1:r.7_*2ac[5]");
    assert_eq!(r_twice, r_once, "the r. spelling is a fixed point too");
    assert_eq!(
        once.split_once("c.").expect("c. prefix").1.to_lowercase(),
        r_once.split_once("r.").expect("r. prefix").1,
        "the c. and r. spellings of one axis must name the same edit"
    );
}

/// `con` canonicalizes to a `delins`, so it inherits the fix along with it.
#[test]
fn conversion_crossing_cds_end_settles_in_one_pass() {
    let (once, twice) = normalize_twice(
        transcript_with_utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:c.7_9con7_7",
    );
    assert_eq!(
        once, "NM_UTR3.1:c.*39_*40del",
        "con reduces to the same fully shifted deletion as the delins spelling"
    );
    assert_eq!(twice, once, "and is a fixed point");
}

/// Fixed on the minus strand too, so the repair is not strand-specific.
#[test]
fn delins_reducing_to_deletion_settles_in_one_pass_on_minus_strand() {
    let (once, twice) = normalize_twice(
        transcript_with_utr3_ac_tandem(Strand::Minus),
        "NM_UTR3.1:c.7_9delinsA",
    );
    assert_eq!(once, "NM_UTR3.1:c.*39_*40del");
    assert_eq!(twice, once, "and is a fixed point");
}

// ---------------------------------------------------------------------------
// Controls. These passed even while the cases above were broken, and MUST keep
// passing: they are what established that the defect was a `c.`-axis code-path
// bug rather than a property of the reference sequence or of 3'-shifting across
// `cds_end` in general. They now also pin the agreement between axes.
// ---------------------------------------------------------------------------

/// **The control that ruled out a fixture artifact.** `r.` on a coding
/// transcript is the CDS-relative axis — the same axis as `c.` (#469) — so
/// `r.7_9delinsa` names exactly the bases `c.7_9delinsA` names, over the
/// byte-identical transcript. It shifts across `cds_end` and settles in ONE
/// pass, which the `c.` spelling of the same edit now does too (see
/// [`delins_reducing_to_deletion_shifts_across_cds_end_in_one_pass`], whose
/// answer is this one re-spelled on the `c.` axis).
#[test]
fn rna_axis_shifts_across_cds_end_in_one_pass() {
    let (once, twice) = normalize_twice(
        transcript_with_utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:r.7_9delinsa",
    );
    assert_eq!(
        once, "NM_UTR3.1:r.*39_*40del",
        "the r. axis completes the cross-cds_end shift on the first pass"
    );
    assert_eq!(once, twice, "and is idempotent");
}

/// A plain `del` on the `c.` axis crosses `cds_end` in one pass, so the bug is
/// not "the `c.` axis cannot shift into the 3'UTR".
#[test]
fn plain_deletion_crosses_cds_end_in_one_pass() {
    let (once, twice) = normalize_twice(
        transcript_with_utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:c.*1_*2del",
    );
    assert_eq!(once, "NM_UTR3.1:c.*39_*40del");
    assert_eq!(once, twice);
}

/// The same `delins` reduction is idempotent when the shift target stays
/// **inside** the CDS, so the trigger is specifically crossing `cds_end`.
/// `c.3_5delinsA` reduces to a deletion of `AA` that shifts along the interior
/// `A`-run (`c.3..c.8`) and lands at `c.7_8`, still inside the CDS.
#[test]
fn delins_reduction_is_idempotent_when_shift_stays_inside_cds() {
    let (once, twice) = normalize_twice(
        transcript_with_utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:c.3_5delinsA",
    );
    assert_eq!(once, "NM_UTR3.1:c.7_8del");
    assert_eq!(once, twice);

    // Cross-axis oracle, as above: the `r.` spelling of the same edit resolves
    // through `normalize_rna` and must land on the same two bases.
    let (r_once, r_twice) = normalize_twice(
        transcript_with_utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:r.3_5delinsa",
    );
    assert_eq!(r_once, "NM_UTR3.1:r.7_8del");
    assert_eq!(r_twice, r_once, "the r. spelling is a fixed point too");
    assert_eq!(
        once.split_once("c.").expect("c. prefix").1.to_lowercase(),
        r_once.split_once("r.").expect("r. prefix").1,
        "an interior shift must agree across the two spellings too"
    );
}

/// The genomic axis over the same sequence is idempotent, so the bug is not in
/// the shuffle itself.
#[test]
fn genomic_axis_shifts_in_one_pass() {
    let mut provider = MockProvider::new();
    // Same layout as the transcript, addressed as a contig: the "CDS" bases sit
    // at g.41..=49 and the `AC` tandem continues from g.48.
    provider.add_genomic_sequence(
        "NC_UTR3.1",
        format!("{}{}{}", "C".repeat(40), "GGAAAAAAC", "AC".repeat(20)),
    );
    let (once, twice) = normalize_twice(provider, "NC_UTR3.1:g.47_49delinsA");
    // The `AC` tandem runs g.48..g.89 (the contig ends at g.89), so the
    // reduced `AC` deletion 3'-shifts to the final copy.
    assert_eq!(
        once, "NC_UTR3.1:g.88_89del",
        "the genomic axis shifts to the end of the tandem in one pass"
    );
    assert_eq!(twice, once, "and is idempotent over the same bases");
}
