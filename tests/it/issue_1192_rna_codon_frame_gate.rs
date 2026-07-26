//! Issue #1192 — the `r.` axis must apply the codon-frame gate to repeat
//! notation, exactly as the `c.` axis does.
//!
//! `RNA/repeated.md` L24-27 states the restriction for an RNA reference in the
//! same terms the DNA page uses, including the UTR carve-out:
//!
//! > **exception**: using a coding RNA reference sequence, a repeated sequence
//! > variant description can be used only for repeat units with a length which
//! > is a multiple of 3, i.e. which can not affect the reading frame.
//! > Consequently, use `NM_024312.4:r.2692_2693dup` and **not**
//! > `NM_024312.4:r.2686a[10]` … This restriction only applies to the coding
//! > sequence, which does not include the UTR sequence. As such,
//! > `NM_024312.4:r.-6_-3g[6]` is valid as the reading frame is not affected.
//!
//! `r.` on a coding transcript is CDS-relative — the same axis as `c.`
//! (`background/numbering.md` L58/L61, pinned by #469) — so a `c.`/`r.` pair
//! over the same bases must reach the same decision. `normalize_rna` used to
//! pass `is_coding = false` unconditionally, so the gate never fired and ferro
//! emitted, on `r.`, the very form the spec marks invalid: with the real
//! reference, `NM_024312.4:r.2692_2693dup` rendered as `r.2686_2693a[10]`.
//!
//! Each case is asserted on **both** axes. That pairing is what makes this a
//! spec defect rather than a fixture artifact: the two axes address identical
//! reference bases, so any divergence between them is the bug itself.
//!
//! The repro is a homopolymer duplication because that is the shape the spec's
//! own counter-example takes — a unit length of 1 is never a multiple of 3, so
//! it can never be codon-aligned.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Length of the 5'UTR, so `c.1` / `r.1` is transcript base `UTR5_LEN + 1`.
const UTR5_LEN: usize = 30;

/// Build a single-exon coding transcript `C * 30` (5'UTR) ++ `cds` ++ `utr3`.
///
/// Hand-rolled rather than built with `SyntheticBuilder`: that helper pads with
/// `ACGT` repeats, which a tandem unit can cycle straight into, and a padding
/// artifact of exactly that shape produced a retracted finding once before.
/// Every tract below is flanked by a base that cannot continue it.
fn coding_transcript(cds: &str, utr3: &str) -> MockProvider {
    let sequence = format!("{}{cds}{utr3}", "C".repeat(UTR5_LEN));
    let tx_len = sequence.len() as u64;
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_RFRAME.1".to_string(),
        Some("RFRAME_TEST".to_string()),
        Strand::Plus,
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

/// Build a single-exon **non-coding** transcript, with no CDS annotation at all.
fn noncoding_transcript(sequence: &str) -> MockProvider {
    let tx_len = sequence.len() as u64;
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NR_RFRAME.1".to_string(),
        Some("RFRAME_NC_TEST".to_string()),
        Strand::Plus,
        sequence.to_string(),
        None,
        None,
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

/// CDS `ATG` ++ `A`×6 ++ `GGGTAA`: a homopolymer tract at `c.4_9`, flanked by
/// `G` on both sides and well clear of both UTRs. A repeat unit of length 1 is
/// never a multiple of 3, so the codon-frame gate must block repeat notation.
fn homopolymer_tract_in_cds() -> MockProvider {
    coding_transcript("ATGAAAAAAGGGTAA", &"G".repeat(30))
}

/// CDS `ATG` ++ `AGC`×4 ++ `TAA`: a codon-aligned tri-nucleotide tract at
/// `c.4_15` that the gate must leave alone.
fn codon_aligned_tract_in_cds() -> MockProvider {
    coding_transcript("ATGAGCAGCAGCAGCTAA", &"G".repeat(30))
}

/// A CDS with no repeat structure, and the homopolymer tract moved into the
/// 3'UTR at `r.*4_*9` — where the spec's carve-out exempts it from the gate.
fn homopolymer_tract_in_utr3() -> MockProvider {
    coding_transcript("ATGGGGTTTTAA", &format!("GGG{}GGGGGGGGGG", "A".repeat(6)))
}

/// Normalize `input` and render the result.
fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).expect("parse");
    normalizer
        .normalize(&variant)
        .expect("normalize")
        .to_string()
}

// ---------------------------------------------------------------------------
// The defect: a non-triplet repeat unit inside the CDS proper
// ---------------------------------------------------------------------------

/// The `c.` spelling, pinning the behaviour the `r.` axis has to match: the
/// gate fires, so `A[8]` is refused and the description stays a `dup`.
///
/// Re-blessed by #1204 from `c.9_10insAA`. What this test is for is unchanged —
/// the absence of repeat notation is the assertion — but the *fallback* is now the
/// form the same spec sentence prescribes: "use `NM_024312.4:c.2692_2693dup` and
/// **not** `c.2686A[10]`", i.e. a `dup` over the tract's 3'-most pair, which here
/// is `c.8_9` and makes the input its own canonical form.
#[test]
fn cds_homopolymer_expansion_inside_cds_is_gated() {
    assert_eq!(
        normalize(homopolymer_tract_in_cds(), "NM_RFRAME.1:c.8_9dup"),
        "NM_RFRAME.1:c.8_9dup",
    );
}

/// The same bases on `r.`. Before #1194 this rendered as `r.4_9a[8]` — the shape
/// the spec explicitly forbids on a coding RNA reference — and before #1204 as
/// `r.9_10insaa`. Both axes moved together at each step, which is this file's
/// contract.
#[test]
fn rna_homopolymer_expansion_inside_cds_is_gated() {
    assert_eq!(
        normalize(homopolymer_tract_in_cds(), "NM_RFRAME.1:r.8_9dup"),
        "NM_RFRAME.1:r.8_9dup",
    );
}

// ---------------------------------------------------------------------------
// Controls — each of these must hold both before and after the fix
// ---------------------------------------------------------------------------

/// A codon-aligned (multiple-of-3) unit inside the CDS is untouched by the
/// gate and still earns repeat notation on `c.`.
#[test]
fn cds_codon_aligned_repeat_inside_cds_is_allowed() {
    let rendered = normalize(codon_aligned_tract_in_cds(), "NM_RFRAME.1:c.10_15dup");
    assert!(
        rendered.contains('['),
        "a 3-base repeat unit inside the CDS keeps repeat notation on c., got {rendered}",
    );
}

/// …and likewise on `r.`. Wiring the gate must not cost the `r.` axis repeat
/// notation it is entitled to.
#[test]
fn rna_codon_aligned_repeat_inside_cds_is_allowed() {
    let rendered = normalize(codon_aligned_tract_in_cds(), "NM_RFRAME.1:r.10_15dup");
    assert!(
        rendered.contains('['),
        "a 3-base repeat unit inside the CDS keeps repeat notation on r., got {rendered}",
    );
}

/// The spec's carve-out on `c.`: the restriction "only applies to the coding
/// sequence, which does not include the UTR sequence".
#[test]
fn cds_homopolymer_expansion_in_utr3_is_allowed() {
    let rendered = normalize(homopolymer_tract_in_utr3(), "NM_RFRAME.1:c.*8_*9dup");
    assert!(
        rendered.contains('['),
        "a homopolymer in the 3'UTR keeps repeat notation on c., got {rendered}",
    );
}

/// …and the same carve-out on `r.`: the gate must key on the footprint's
/// position, not merely on the transcript being coding.
#[test]
fn rna_homopolymer_expansion_in_utr3_is_allowed() {
    let rendered = normalize(homopolymer_tract_in_utr3(), "NM_RFRAME.1:r.*8_*9dup");
    assert!(
        rendered.contains('['),
        "a homopolymer in the 3'UTR keeps repeat notation on r., got {rendered}",
    );
}

/// The **5'UTR** half of the carve-out, and the case the spec actually cites
/// (`NM_024312.4:r.-6_-3g[6]` "is valid as the reading frame is not affected").
/// This exercises the *other* arm of the `is_coding` condition from the 3'UTR
/// tests above — `tx_start < cds_start` rather than `tx_end > cds_end` — so
/// neither test alone covers both sides of the footprint check.
///
/// Derived, not observed: the 5'UTR is `C`×30 at tx 1..30 == `r.-30_-1`, and
/// `r.-2_-1` is tx 29_30, so `dup` lands wholly in the UTR (29 < cds_start 31)
/// and the gate must not fire. The CDS opens with `A`, which cannot continue the
/// `C` run, so the tract is exactly 30 copies of unit `c` and the dup takes it
/// to 32.
#[test]
fn cds_homopolymer_expansion_in_utr5_is_allowed() {
    assert_eq!(
        normalize(homopolymer_tract_in_cds(), "NM_RFRAME.1:c.-2_-1dup"),
        "NM_RFRAME.1:c.-30_-1C[32]",
    );
}

/// …and the same on `r.`, where the spec states the rule.
#[test]
fn rna_homopolymer_expansion_in_utr5_is_allowed() {
    assert_eq!(
        normalize(homopolymer_tract_in_cds(), "NM_RFRAME.1:r.-2_-1dup"),
        "NM_RFRAME.1:r.-30_-1c[32]",
    );
}

/// A non-coding transcript has no reading frame to preserve, so `r.` — which is
/// `n.`-equivalent there — must never gate.
#[test]
fn rna_homopolymer_expansion_on_noncoding_transcript_is_allowed() {
    let provider = noncoding_transcript("GGGAAAAAAGGGGGGGGGG");
    let rendered = normalize(provider, "NR_RFRAME.1:r.8_9dup");
    assert!(
        rendered.contains('['),
        "a non-coding transcript has no reading frame and must not gate, got {rendered}",
    );
}
