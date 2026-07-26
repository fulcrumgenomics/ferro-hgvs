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

/// The SAME edit over the BYTE-IDENTICAL transcript, differing only in the axis
/// it is spelled on. On a coding transcript `r.` is the CDS-relative axis — the
/// same axis as `c.` (#469) — so both spellings name the same bases.
///
/// **Scope, deliberately narrow.** This is a valid control for the *shuffle*
/// (this test), because no repeat notation is involved. It is **not** a valid
/// control for the codon-frame gate exercised by the `dup` tests below:
/// `normalize_rna` passes no CDS span at all, so the gate never applies on the
/// `r.` axis — even though `RNA/repeated.md:24-27` imposes the same
/// multiple-of-3 restriction as the DNA page. That gap is a separate defect and
/// is deliberately not addressed here; treating `r.` as an oracle for the gate
/// would have "confirmed" whatever `r.` happened to do.
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

// ---------------------------------------------------------------------------
// Site 2 — the codon-frame gate must be evaluated against the tract the repeat
// description will NAME (the shuffled tract), not the input span.
//
// `DNA/repeated.md`: "using a coding DNA reference sequence ("c." description),
// a repeated sequence variant description can be used only for repeat units
// with a length which is a multiple of 3, i.e. which can not affect the reading
// frame. … **This restriction only applies to the coding sequence, which does
// not include the introns or the UTR sequence.** As such,
// `NM_024312.4:c.-6_-3G[6]` is valid as the reading frame is not affected."
//
// So a tract that 3'-shifts out of the CDS into the 3'UTR is exempt: the added
// copies land past `cds_end` and cannot move the reading frame.
// ---------------------------------------------------------------------------

/// CDS `GAAAAAACAC` (`c.1..c.10`) + an `ACGT` 3'UTR puts an `AC[3]` tract at
/// `c.7_*2`, bounded by the `G` at `*3`.
fn utr3_acgt_bounded() -> MockProvider {
    coding_transcript("GAAAAAACAC", &"ACGT".repeat(10), Strand::Plus)
}

/// `c.7_10dup` duplicates `ACAC`, expanding the `AC[3]` tract at `c.7_*2` to
/// `AC[5]`. The unit is 2 nt — not a multiple of 3 — but the tract runs out of
/// the CDS into the 3'UTR, so the codon-frame restriction does not apply and
/// repeat notation is correct.
///
/// Before the fix the gate was keyed on the input span (`c.7_10`, entirely
/// CDS-resident), so pass 1 emitted the `ins` literal `c.*2_*3insACAC` and only
/// pass 2 — by then seeing a straddling span — collapsed it to `AC[5]`.
#[test]
fn dup_whose_shifted_tract_leaves_the_cds_gets_repeat_notation_in_one_pass() {
    assert_one_pass(
        utr3_acgt_bounded(),
        "NM_UTR3.1:c.7_10dup",
        "NM_UTR3.1:c.7_*2AC[5]",
    );
}

/// The control that proves the gate was narrowed, not disabled. A non-triplet
/// unit whose tract stays **entirely inside** the CDS must still be refused
/// repeat notation — no `AC[6]` here, unlike the straddling case above.
///
/// Re-blessed by #1204 from the `ins` literal `c.9_10insACAC`. The gate forbids
/// repeat notation, and of the spec's two replacements this input takes the first:
/// the added `ACAC` copies `c.6_9`, the 3'-most four bases of the `c.2_9` tract, so
/// it is a duplication ("use `NM_024312.4:c.2692_2693dup`"). The `ins` remedy
/// ("use `c.1741_1742insTATATATA` and **not** `c.1738TA[6]`") is for added bases no
/// duplication covers, which are longer than the tract they land in; that shape is
/// pinned in `issue_1204_gated_dup_fallback` and `ins_shift_matrix`'s case 2.
#[test]
fn control_non_triplet_dup_inside_the_cds_is_still_refused_repeat_notation() {
    assert_one_pass(
        coding_transcript("GACACACACGGGTTT", &"GGGTTT".repeat(4), Strand::Plus),
        "NM_UTR3.1:c.2_5dup",
        "NM_UTR3.1:c.6_9dup",
    );
}

/// And a codon-aligned unit inside the CDS still gets repeat notation, so the
/// gate is discriminating on unit length rather than blanket-refusing.
#[test]
fn control_codon_aligned_dup_inside_the_cds_gets_repeat_notation() {
    assert_one_pass(
        coding_transcript("GCAGCAGCAGTTT", &"GGGTTT".repeat(4), Strand::Plus),
        "NM_UTR3.1:c.2_7dup",
        "NM_UTR3.1:c.2_10CAG[5]",
    );
}
