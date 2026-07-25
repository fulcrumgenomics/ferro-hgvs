//! Issue #1185 — a `delins` that reduces to a pure deletion must complete its
//! 3'-shift across `cds_end` in a **single** normalization pass.
//!
//! The HGVS 3' rule is unconditional — "for all descriptions, the most 3'
//! position possible **of the reference sequence** is arbitrarily assigned to
//! have been changed" (`general.md`; `DNA/deletion.md` L20) — and the coding DNA
//! reference is one contiguous string (5'UTR+CDS+3'UTR, `background/refseq.md`),
//! so a shift from `c.<N>` into `c.*<M>` stays inside that single sequence.
//! \#918 already relaxed the CDS↔3'UTR shuffle bound for a CDS-resident
//! `del`/`dup` for exactly this reason.
//!
//! That relaxation was gated on the **input's** edit-kind spelling, so a
//! `delins` whose alt/ref share a prefix — which reduces to a pure deletion
//! during normalization — never received it. Its shift was clamped at
//! `cds_end`, and the reduced deletion only shifted on a *second* pass, once it
//! re-entered spelled as a `Deletion`. Normalization was therefore not
//! idempotent for these inputs.
//!
//! The controls below are what make this a bug rather than a fixture artifact:
//! the same edit on the `r.` axis, the already-reduced `del` spelling, the
//! `n.` axis over the same bases, and a `delins` whose shift stays inside the
//! CDS all settle in one pass — before and after the fix.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Length of the 5'UTR, so `c.1` is transcript base `UTR5_LEN + 1`.
const UTR5_LEN: usize = 40;

/// Build `NM_UTR3.1` as `C * 40` (5'UTR) ++ `cds` ++ `utr3`, single exon.
///
/// Hand-rolled rather than built with `SyntheticBuilder` so that no property of
/// the synthetic `ACGT` padding can be blamed for the outcome — a padding
/// artifact of exactly this shape produced a false finding once before.
fn coding_transcript(cds: &str, utr3: &str, strand: Strand) -> MockProvider {
    let sequence = format!("{}{cds}{utr3}", "C".repeat(UTR5_LEN));
    let tx_len = sequence.len() as u64;
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_UTR3.1".to_string(),
        Some("UTR3_TEST".to_string()),
        strand,
        sequence,
        Some((UTR5_LEN + 1) as u64),
        Some((UTR5_LEN + cds.len()) as u64),
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

/// CDS `GGAAAAAAC` (`c.1..c.9`) ending in `AC`, followed by an `AC` 3'UTR — so
/// a deletion of that `AC` has an unbroken `AC` tandem to 3'-shift along, and
/// the shift necessarily crosses `cds_end` and must render with `*N`.
fn utr3_ac_tandem(strand: Strand) -> MockProvider {
    coding_transcript("GGAAAAAAC", &"AC".repeat(20), strand)
}

/// Normalize once and twice, returning both renderings.
fn normalize_twice(provider: MockProvider, input: &str) -> (String, String) {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse");
    let once = normalizer.normalize(&variant).expect("first normalize");
    let twice = normalizer.normalize(&once).expect("re-normalize");
    (once.to_string(), twice.to_string())
}

/// Assert `input` reaches `expected` on the FIRST pass and stays there.
fn assert_one_pass(provider: MockProvider, input: &str, expected: &str) {
    let (once, twice) = normalize_twice(provider, input);
    assert_eq!(
        once, expected,
        "{input} must reach its most-3' form in one pass",
    );
    assert_eq!(
        twice, once,
        "{input} normalized to {once}, which is not a fixed point",
    );
}

/// `c.7_9delinsA`: ref `c.7_9` = `AAC`, alt `A`. The shared `A` prefix reduces
/// it to a deletion of `AC` at `c.8_9`, which 3'-shifts along the `AC` tandem
/// running out of the CDS into the 3'UTR.
///
/// Before the fix this stopped at `c.8_9del` and only shifted on a second pass.
#[test]
fn delins_reducing_to_deletion_shifts_across_cds_end_in_one_pass() {
    assert_one_pass(
        utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:c.7_9delinsA",
        "NM_UTR3.1:c.*39_*40del",
    );
}

/// Strand must not matter: the coding reference is read in transcript
/// orientation either way.
#[test]
fn delins_reducing_to_deletion_shifts_in_one_pass_on_minus_strand() {
    assert_one_pass(
        utr3_ac_tandem(Strand::Minus),
        "NM_UTR3.1:c.7_9delinsA",
        "NM_UTR3.1:c.*39_*40del",
    );
}

// A `con` inherits this fix for free rather than needing its own test: the
// SVD-WG009 `con` → `delins` rewrite (`src/normalize/mod.rs`, the
// `canonicalize_conversion_to_delins` branch) recurses into `normalize_cds`, so
// every downstream pass — including the relaxation fixed here — sees the
// `delins` form and nothing else.
//
// There is deliberately no `con` case below. Reaching one would need a donor
// range whose copied bases reduce the target to a pure deletion crossing
// `cds_end`, which on this fixture means a `*N`-marked payload
// (`c.7_9con*39_*41`) — and `ins[...]` expansion declines `*N` markers by
// design ("Out-of-scope: … UTR markers (*N)"). So the `con` spelling errors on
// the payload long before the shuffle, independently of this issue.

// ---------------------------------------------------------------------------
// Controls — these all passed BEFORE the fix. They are what localise the
// defect to the `delins`-on-`c.` path rather than to the shuffle, the fixture,
// or the 3'UTR rendering.
// ---------------------------------------------------------------------------

/// The already-reduced spelling. The shuffle machinery itself was never broken:
/// fed a `Deletion` directly, it crosses `cds_end` in one pass.
#[test]
fn control_plain_deletion_crosses_cds_end_in_one_pass() {
    assert_one_pass(
        utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:c.8_9del",
        "NM_UTR3.1:c.*39_*40del",
    );
}

/// The decisive control: the SAME edit over the BYTE-IDENTICAL transcript,
/// differing only in the axis it is spelled on. On a coding transcript `r.` is
/// the CDS-relative axis — the same axis as `c.` (#469) — so both spellings
/// name the same bases and must shift identically.
#[test]
fn control_same_edit_on_r_axis_shifts_in_one_pass() {
    assert_one_pass(
        utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:r.7_9delinsa",
        "NM_UTR3.1:r.*39_*40del",
    );
}

/// The transcript-relative axis over the same bases (`n.48_49` is `c.8_9`),
/// where there is no CDS↔3'UTR boundary to clamp against at all.
#[test]
fn control_n_axis_over_the_same_bases_shifts_in_one_pass() {
    assert_one_pass(
        utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:n.48_49del",
        "NM_UTR3.1:n.88_89del",
    );
}

/// A `delins` that reduces to a deletion whose 3'-shift stays INSIDE the CDS
/// was always one-pass — so the defect was specifically about crossing
/// `cds_end`, not about `delins` reduction in general.
#[test]
fn control_delins_reduction_staying_inside_cds_is_one_pass() {
    assert_one_pass(
        utr3_ac_tandem(Strand::Plus),
        "NM_UTR3.1:c.3_5delinsA",
        "NM_UTR3.1:c.7_8del",
    );
}
