//! Issue #1210: the codon-frame gate must be answered for the span the emitted
//! **repeat** occupies, not for the span the input edit occupied.
//!
//! ## The defect
//!
//! `normalize_na_edit` receives `is_coding` — a verdict precomputed for the
//! *input* edit's span — and the insertion→repeat site passed it straight to
//! `rules::insertion_to_repeat`. That is the wrong question whenever the edit
//! 3'-shifts out of the CDS: the input is CDS-resident, so the gate fires and
//! suppresses repeat notation, but the tract the repeat would occupy is in the
//! 3'UTR, where `RNA/repeated.md` L24-27's codon constraint does not apply.
//!
//! A second pass, whose input is by then UTR-resident, computes `is_coding =
//! false` and emits the repeat — so the two passes disagreed:
//!
//! ```text
//! r.15delinsgaa -> r.*1_*2insaa -> r.*1a[3]
//! ```
//!
//! `cds_span` already existed for exactly this — its contract is that it "lets
//! a site that decides about a DIFFERENT span ask about that span instead"
//! (#1185). This is the second site to need it.
//!
//! ## Which answer is correct
//!
//! `r.*1a[3]`. The tract is entirely 3'UTR-resident, so the CDS-only
//! prohibition on a 1-nt repeat unit does not reach it, and it is the fixed
//! point both passes now agree on. Pass 1's suppression was the error, not pass
//! 2's spelling.
//!
//! Found by the widened idempotency fuzz in #1180, at ~150k cases — well past
//! the 12k that suite runs by default, which is why `FERRO_ASSERT_IDEMPOTENT=1`
//! plus a soak are the gates that actually cover this class.

use crate::common::synthetic::{SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider, Normalizer};

/// The proptest-reduced fixture: a coding transcript whose CDS is exactly
/// `core`, wrapped in a real `ACGT…` 5'UTR and 3'UTR.
///
/// The padding is built locally rather than via `common::synthetic::padded`,
/// which is private on this branch. The `ACGT` period-4 pad is what bounds the
/// tract: the base immediately 3' of the CDS is `A` and the next is `C`, so an
/// `a` tract running out of the CDS terminates one base into the UTR.
fn coding_transcript_with_utrs(core: &str, strand: Strand) -> MockProvider {
    let pad = "ACGT".repeat(64);
    assert_eq!(pad.len() as u64, PAD_OFFSET, "pad must be PAD_OFFSET bases");
    let sequence = format!("{pad}{core}{pad}");
    let len = core.len() as u64;
    SyntheticBuilder::cds(&sequence, PAD_OFFSET + 1, PAD_OFFSET + len, strand).build()
}

/// `CCCCAACACACGGGG` — 15 bases, so `cds_end` is `r.15`/`c.15` (a `G`), and the
/// 3'UTR begins `A C G T`. A `delins` at the last CDS base that reduces to an
/// insertion of `aa` therefore shifts into the UTR and meets the lone `A` at
/// `*1`, forming an `a[3]` tract wholly outside the CDS.
const CORE: &str = "CCCCAACACACGGGG";

fn normalize_twice(provider: MockProvider, input: &str) -> (String, String) {
    let normalizer = Normalizer::new(provider);
    let variant: HgvsVariant = parse_hgvs(input).expect("parse failed");
    let once = normalizer
        .normalize(&variant)
        .expect("first normalize failed");
    let twice = normalizer.normalize(&once).expect("re-normalize failed");
    (once.to_string(), twice.to_string())
}

/// The reported case: repeat notation is emitted on the **first** pass, because
/// the tract it describes is 3'UTR-resident even though the input was not.
#[test]
fn rna_repeat_gate_uses_the_tract_span_not_the_input_span() {
    assert_eq!(CORE.len(), 15, "cds_end must be r.15");
    let (once, twice) = normalize_twice(
        coding_transcript_with_utrs(CORE, Strand::Plus),
        "NM_TEST.1:r.15delinsgaa",
    );
    assert_eq!(
        once, "NM_TEST.1:r.*1a[3]",
        "pass 1 must emit the repeat: the tract is in the 3'UTR, not the CDS"
    );
    assert_eq!(twice, once, "and it is a fixed point");
}

/// The `c.` spelling goes through the same site, so it must move with it. On a
/// coding transcript `r.` *is* the `c.` axis (#469), so this is the same edit
/// over the byte-identical transcript, and the two must agree.
#[test]
fn cds_spelling_agrees_with_the_rna_spelling() {
    let (c_once, c_twice) = normalize_twice(
        coding_transcript_with_utrs(CORE, Strand::Plus),
        "NM_TEST.1:c.15delinsGAA",
    );
    assert_eq!(c_once, "NM_TEST.1:c.*1A[3]");
    assert_eq!(c_twice, c_once, "and is a fixed point");

    let (r_once, _) = normalize_twice(
        coding_transcript_with_utrs(CORE, Strand::Plus),
        "NM_TEST.1:r.15delinsgaa",
    );
    // Compare only the position+edit tail, so the accession is not folded: the
    // axes differ exactly by the prefix and the DNA/RNA alphabet.
    let c_tail = c_once.split_once("c.").expect("c. prefix").1.to_lowercase();
    let r_tail = r_once.split_once("r.").expect("r. prefix").1;
    assert_eq!(
        c_tail, r_tail,
        "the two spellings of one axis must name the same edit"
    );
}

/// Present on the minus strand too, so the fix is not strand-specific.
#[test]
fn rna_repeat_gate_tract_span_on_minus_strand() {
    let (once, twice) = normalize_twice(
        coding_transcript_with_utrs(CORE, Strand::Minus),
        "NM_TEST.1:r.15delinsgaa",
    );
    assert_eq!(once, "NM_TEST.1:r.*1a[3]");
    assert_eq!(twice, once, "and is a fixed point");
}

// ---------------------------------------------------------------------------
// Guards: the gate must still fire where it should. Re-answering the question
// for the tract is only correct if a CDS-resident tract still gets gated.
// ---------------------------------------------------------------------------

/// **The gate is still armed.** The same shape with the tract kept *inside* the
/// CDS must still be refused repeat notation, because a 1-nt repeat unit in the
/// CDS is spec-forbidden (`RNA/repeated.md` L24-27, the #1192 gate).
///
/// `CCCCAAAAAAAAGGGG` puts a long `A` run at `c.5..c.12`, comfortably interior,
/// so an inserted `aa` forms a tract that never approaches `cds_end`.
#[test]
fn interior_tract_is_still_codon_gated() {
    let (once, twice) = normalize_twice(
        coding_transcript_with_utrs("CCCCAAAAAAAAGGGG", Strand::Plus),
        "NM_TEST.1:r.11_12insaa",
    );
    // Gated: the tract is CDS-resident and the unit is 1 nt, so repeat notation
    // is refused — which is what this guard is for, and is unchanged.
    assert!(
        !once.contains('['),
        "a CDS-interior 1-nt tract must not become a repeat, got {once}",
    );
    // The fallback form is the `dup` `repeated.md` L22 prescribes for this shape
    // ("use `NM_024312.4:c.2692_2693dup`"): the added `aa` copies `r.11_12`, the
    // 3'-most pair of the `A` run at `r.5_12`. Re-blessed from `r.12_13insaa` by
    // #1204, which is about the fallback and not about whether the gate fires.
    assert_eq!(
        once, "NM_TEST.1:r.11_12dup",
        "a CDS-interior 1-nt tract must stay gated, falling back to the dup"
    );
    assert_eq!(twice, once, "and is a fixed point");
}

/// The non-coding axis has no CDS at all, so `cds_span` is `None` and the
/// helper keeps the caller's verdict verbatim — the `None` branch of the fix.
///
/// This doubles as an **independent oracle for the headline case**. `NR_TEST.1`
/// carries the byte-identical sequence, and `n.271` is the same base as `r.15`
/// (both are transcript position `PAD_OFFSET + 15`), so `n.271delinsGAA` is the
/// same edit over the same bases — but reached with the gate never consulted at
/// all, rather than consulted and answered `false`. It lands on the repeat, at
/// the transcript position `r.*1` denotes:
///
/// ```text
/// n.271delinsGAA -> n.272A[3]      (no CDS, gate never consulted)
/// r.15delinsgaa  -> r.*1a[3]       (CDS present, gate answered for the tract)
/// ```
///
/// Before the fix these two disagreed, which is the clearest statement that the
/// `r.` answer was wrong rather than merely differently spelled.
#[test]
fn noncoding_axis_is_unaffected() {
    let pad = "ACGT".repeat(64);
    let provider = SyntheticBuilder::noncoding(&format!("{pad}{CORE}{pad}"), Strand::Plus).build();
    let (once, twice) = normalize_twice(provider, "NR_TEST.1:n.271delinsGAA");
    assert_eq!(
        once, "NR_TEST.1:n.272A[3]",
        "with no CDS the gate never fires, so the repeat is emitted directly"
    );
    assert_eq!(twice, once, "n. normalization is idempotent here");

    // The `n.` position and the `r.` `*N` position name the same base:
    // `n.272` == transcript `PAD_OFFSET + 16` == `r.*1` (cds_end is at +15).
    assert_eq!(PAD_OFFSET + 16, 272, "n.272 must be the base r.*1 denotes");
}
