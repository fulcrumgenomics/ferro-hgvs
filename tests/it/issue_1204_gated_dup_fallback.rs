//! Issue #1204 — where the codon-frame gate refuses repeat notation inside the
//! CDS, the fallback must be the spec's `dup`, not a flat `ins` literal.
//!
//! `DNA/repeated.md` L21-22 states the exception and, crucially, gives **two**
//! replacements of **different shapes**:
//!
//! > **exception**: using a coding DNA reference sequence ("c." description), a
//! > repeated sequence variant description can be used only for repeat units
//! > with a length which is a multiple of 3 […]
//! > Consequently, use `NM_024312.4:c.2692_2693dup` and **not**
//! > `NM_024312.4:c.2686A[10]`; use `NM_024312.4:c.1741_1742insTATATATA` and
//! > **not** `NM_024312.4:c.1738TA[6]`.
//!
//! `RNA/repeated.md` L25 states it identically for a coding RNA reference
//! (`r.2692_2693dup`, `r.1741_1742insuauauaua`).
//!
//! The two replacements are not interchangeable, and which one applies is
//! decided by the reference, not by the gate:
//!
//! - `c.2692_2693dup` — the added bases duplicate the immediately adjacent
//!   reference bases (`A[8]` + `AA`: the two added `A`s are a copy of the tract's
//!   two 3'-most bases), so the change *is* a duplication and `dup` describes it.
//! - `c.1741_1742insTATATATA` — four `TA` copies added to a two-copy `TA` tract.
//!   The added 8 bases are longer than the whole adjacent tract, so no `dup`
//!   describes them and the flat literal is the only form left.
//!
//! ferro emitted the `ins` form for both: the gate that (correctly) refuses
//! repeat notation also suppressed the `ins` → `dup` promotion, so the
//! duplication half of the guidance was unreachable. The gate belongs only in
//! the repeat path — `dup` is not repeat notation and the exception says nothing
//! against it.
//!
//! Every case is asserted on **both** the `c.` and `r.` axes. `r.` on a coding
//! transcript IS the `c.` axis (`background/numbering.md` L58/L61, #469), so the
//! pairing is what makes a divergence a defect rather than a fixture artifact —
//! the same reasoning #1192 and #1194 used for the gate itself.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Length of the 5'UTR, so `c.1` / `r.1` is transcript base `UTR5_LEN + 1`.
const UTR5_LEN: usize = 30;

/// Build a single-exon coding transcript `C * 30` (5'UTR) ++ `cds` ++ `utr3`.
///
/// Padded with `C` rather than via `SyntheticBuilder`'s `ACGT` repeats: an
/// `ACGT` pad contains both an `A` and a `TA`, either of which a tract under
/// test could cycle straight into, and a padding artifact of exactly that shape
/// produced a retracted finding once before (see #1192's fixture note). Every
/// tract below is flanked by a base that cannot continue it.
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

/// The spec's **first** replacement, transposed onto a short fixture.
///
/// CDS `ATG` ++ `A`×6 ++ `GGGTAA`, i.e. an `A[6]` tract at `c.4_9` flanked by
/// `G` on both sides (`c.3` = `G`, `c.10` = `G`) and clear of both UTRs. Adding
/// two `A`s takes it to `A[8]`, which the gate must refuse (unit length 1 is
/// never a multiple of 3) — and the two added bases duplicate `c.8_9`, so `dup`
/// is the form the spec prescribes. Same tract as #1192's fixture, so the two
/// files pin the same bases.
fn homopolymer_tract_in_cds() -> MockProvider {
    coding_transcript("ATGAAAAAAGGGTAA", &"G".repeat(30))
}

/// The spec's **second** replacement, transposed onto a short fixture.
///
/// CDS `ATG` ++ `TATA` ++ `GGGTAA`, i.e. a `TA[2]` tract at `c.4_7` flanked by
/// `G` at `c.3` and `c.8`. Adding four `TA` copies takes it to `TA[6]`, which the
/// gate must refuse (unit length 2). The eight added bases are twice the length
/// of the whole adjacent tract, so no duplication describes them and the flat
/// `ins` literal must survive. This is the half of the guidance that was already
/// correct, and the guard that the fix does not over-reach into it.
fn dinucleotide_tract_in_cds() -> MockProvider {
    coding_transcript("ATGTATAGGGTAA", &"G".repeat(30))
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

/// Normalize `input`, assert it renders as `expected`, and assert `expected` is
/// a re-parseable fixed point.
///
/// The fixed-point half matters here beyond the usual idempotency hygiene: the
/// defect's `ins` output and the correct `dup` output are both stable, so an
/// idempotency-only assertion would pass on the bug. The exact-string assertion
/// is what distinguishes them; the fixed point rules out trading this defect for
/// a normalize-twice divergence.
fn assert_normalizes_to(provider: fn() -> MockProvider, input: &str, expected: &str) {
    assert_eq!(
        normalize(provider(), input),
        expected,
        "normalizing {input}"
    );
    assert_eq!(
        normalize(provider(), expected),
        expected,
        "{expected} must be a fixed point",
    );
}

// ---------------------------------------------------------------------------
// The defect: the gated fallback must be `dup` where the change is one
// ---------------------------------------------------------------------------

/// A `dup` spelling of the spec's first case stays a `dup`. Before the fix this
/// rendered as `c.9_10insAA` — the second replacement's shape applied to the
/// first replacement's input.
#[test]
fn cds_gated_homopolymer_dup_stays_a_dup() {
    assert_normalizes_to(
        homopolymer_tract_in_cds,
        "NM_RFRAME.1:c.8_9dup",
        "NM_RFRAME.1:c.8_9dup",
    );
}

/// …and the same bases on `r.`, where #1194 wired the gate. Before the fix:
/// `r.9_10insaa`.
#[test]
fn rna_gated_homopolymer_dup_stays_a_dup() {
    assert_normalizes_to(
        homopolymer_tract_in_cds,
        "NM_RFRAME.1:r.8_9dup",
        "NM_RFRAME.1:r.8_9dup",
    );
}

/// The `ins`-spelled twin of the same change canonicalizes to that same `dup`.
///
/// This is the assertion that matches the spec's own wording most literally —
/// "use `c.2692_2693dup`" is an instruction about the *output* form, whatever
/// the input spelling — and it pins spelling-independence: a `dup` and an `ins`
/// describing identical bases must not normalize to different descriptions.
#[test]
fn cds_gated_homopolymer_insertion_becomes_a_dup() {
    assert_normalizes_to(
        homopolymer_tract_in_cds,
        "NM_RFRAME.1:c.9_10insAA",
        "NM_RFRAME.1:c.8_9dup",
    );
}

/// …and on `r.`, with the RNA `u`/lowercase spelling on input.
#[test]
fn rna_gated_homopolymer_insertion_becomes_a_dup() {
    assert_normalizes_to(
        homopolymer_tract_in_cds,
        "NM_RFRAME.1:r.9_10insaa",
        "NM_RFRAME.1:r.8_9dup",
    );
}

/// The gate itself is untouched: repeat notation is still refused. Redundant
/// with the exact-string assertions above by construction, asserted separately
/// so a future re-blessing cannot quietly hand `A[8]` back — the form the spec
/// marks invalid, and the whole point of #1192/#1194.
#[test]
fn the_gated_dup_is_not_repeat_notation() {
    for input in [
        "NM_RFRAME.1:c.8_9dup",
        "NM_RFRAME.1:r.8_9dup",
        "NM_RFRAME.1:c.9_10insAA",
        "NM_RFRAME.1:r.9_10insaa",
    ] {
        let rendered = normalize(homopolymer_tract_in_cds(), input);
        assert!(
            !rendered.contains('['),
            "repeat notation is forbidden inside the CDS for a unit of length 1, \
             got {rendered} for {input}",
        );
    }
}

// ---------------------------------------------------------------------------
// The other half of the spec's guidance, which must not regress
// ---------------------------------------------------------------------------

/// Four `TA` copies added to a `TA[2]` tract: no `dup` describes eight added
/// bases against a four-base tract, so the flat `ins` literal is correct and
/// must survive. This is the guard that the `dup` fallback is driven by the
/// reference rather than applied to every gated expansion.
#[test]
fn cds_gated_dinucleotide_expansion_stays_an_insertion() {
    assert_normalizes_to(
        dinucleotide_tract_in_cds,
        "NM_RFRAME.1:c.7_8insTATATATA",
        "NM_RFRAME.1:c.7_8insTATATATA",
    );
}

/// …and its `r.` sibling, the form `RNA/repeated.md` L25 spells out.
#[test]
fn rna_gated_dinucleotide_expansion_stays_an_insertion() {
    assert_normalizes_to(
        dinucleotide_tract_in_cds,
        "NM_RFRAME.1:r.7_8insuauauaua",
        "NM_RFRAME.1:r.7_8insuauauaua",
    );
}

/// A *single* added `TA` copy against the same tract is a duplication, and the
/// gate must not suppress it either — the `TA` unit is no more codon-aligned
/// than the homopolymer above, so this pins that the fix keys on "the added
/// bases duplicate the adjacent reference", not on the unit being a homopolymer.
///
/// Derived from the fixture: the tract is `c.4_7` = `TATA`, `c.8` = `G` cannot
/// continue it, so the 3'-most duplicated unit is `c.6_7`.
#[test]
fn cds_gated_dinucleotide_duplication_stays_a_dup() {
    assert_normalizes_to(
        dinucleotide_tract_in_cds,
        "NM_RFRAME.1:c.6_7dup",
        "NM_RFRAME.1:c.6_7dup",
    );
}

/// …and on `r.`.
#[test]
fn rna_gated_dinucleotide_duplication_stays_a_dup() {
    assert_normalizes_to(
        dinucleotide_tract_in_cds,
        "NM_RFRAME.1:r.6_7dup",
        "NM_RFRAME.1:r.6_7dup",
    );
}

// ---------------------------------------------------------------------------
// Controls — the gate's carve-outs are unaffected
// ---------------------------------------------------------------------------

/// Outside the CDS the gate never fires, so repeat notation — not `dup` and not
/// `ins` — remains the canonical form. `DNA/repeated.md` L23-24: "This
/// restriction only applies to the coding sequence, which does not include the
/// introns or the UTR sequence."
#[test]
fn utr5_homopolymer_expansion_still_earns_repeat_notation() {
    assert_normalizes_to(
        homopolymer_tract_in_cds,
        "NM_RFRAME.1:c.-2_-1dup",
        "NM_RFRAME.1:c.-30_-1C[32]",
    );
    assert_normalizes_to(
        homopolymer_tract_in_cds,
        "NM_RFRAME.1:r.-2_-1dup",
        "NM_RFRAME.1:r.-30_-1c[32]",
    );
}
