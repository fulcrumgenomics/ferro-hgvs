//! Issue #1207 — an `r.`-axis insertion that saturates at a CDS boundary must be
//! clamped to a `delins`, exactly as the `c.` axis already does. (See the #1204
//! note below: the `dup`-spelled repros this was filed with now resolve to a `dup`
//! before reaching the clamp, so the clamp itself is pinned on the insertion
//! shapes no duplication describes.)
//!
//! `normalize_cds` has rewritten a saturated insertion at `cds_start` since
//! #1170 and at `cds_end` since #387. `normalize_rna` never got either, which
//! was survivable while the `r.` axis passed `is_coding = false` — the
//! codon-frame gate never fired, so pass 1 went straight to repeat notation and
//! stayed there. #1192 correctly enabled the gate on `r.` (`RNA/repeated.md`
//! L24-27 forbids non-triplet repeat notation on a coding RNA reference), and
//! that turned the missing clamp into a **non-idempotency**:
//!
//! | direction | input | pass 1 | pass 2 |
//! |---|---|---|---|
//! | 5' | `r.1_2dup` | `r.-1_1inscc` | `r.1_4c[6]` |
//! | 3' | `r.13_14dup` | `r.14_*1insgg` | `r.11_14g[6]` |
//!
//! Pass 1 refuses repeat notation because the *input* span is CDS-resident and
//! emits an unclamped insertion across the axis boundary; pass 2 sees a span
//! touching a UTR, does not gate, and produces repeat notation instead. The `c.`
//! spelling of the byte-identical edit was already stable, because the clamp
//! exists there.
//!
//! Every expectation below is derived from the fixture by hand and asserted
//! exactly, and each `r.` case is pinned against its `c.` sibling — the axis
//! parity is the actual contract (#469: on a coding transcript `r.` IS `c.`).
//!
//! **Updated by #1204.** The `dup` inputs in the table above no longer reach the
//! clamp at all: their added bases copy an adjacent same-length tract, so the
//! gated fallback is the `dup` that `repeated.md` L22 prescribes, and a `dup` has
//! two in-CDS anchors — there is no boundary to repair. The non-idempotency this
//! file was filed for is fixed either way, and both regimes are now pinned side by
//! side:
//!
//! | input | added vs. tract | result |
//! |---|---|---|
//! | `r.1_2dup` (5') | 2 added, `C[4]` tract | `r.1_2dup` — no clamp needed |
//! | `r.2_3insccccc` (5') | 5 added, `C[4]` tract | `r.1delinscccccc` — clamped |
//! | `r.13_14dup` (3') | 2 added, `G[4]` tract | `r.13_14dup` — no clamp needed |
//! | `r.12_13insggggg` (3') | 5 added, `G[4]` tract | `r.14delinsgggggg` — clamped |
//!
//! The clamp rows are what keep this file's subject covered: an alt longer than
//! the tract it lands in cannot be spelled as a duplication, so it stays an
//! insertion, saturates the boundary, and must be rewritten.

use crate::common::synthetic::{SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// The core from #1180's minimal failing case: `C`×4, `A`×6, `G`×4.
const CORE: &str = "CCCCAAAAAAGGGG";

/// Coding transcript with a [`PAD_OFFSET`]-base 5'UTR and 3'UTR, so
/// `cds_start == PAD_OFFSET + 1` and both CDS boundaries have a representable
/// far side (`r.-1` / `r.*1`) — which is exactly why the transcript-bound clamp
/// from #1202 does not fire on these cases.
///
/// The padding is built here rather than through `common::synthetic::padded`
/// (private on this branch) so this test needs no visibility change to shared
/// test infrastructure — #1180 widens it for its own use, and this PR is a
/// prerequisite of that one. `"ACGT" × 64` is byte-identical to the helper's
/// `PAD`, so `cds_start`/`cds_end` land exactly where `PAD_OFFSET` says.
/// The core's flanks are safe against that padding: the pad ends `…ACGT` before
/// `CCCC`, and `GGGG` is followed by `ACGT…`, so neither the `C` nor the `G`
/// tract can cycle into the pad (the padding-artifact trap noted in #1192).
fn coding() -> MockProvider {
    let pad = "ACGT".repeat(64);
    assert_eq!(pad.len() as u64, PAD_OFFSET, "pad must match PAD_OFFSET");
    let sequence = format!("{pad}{CORE}{pad}");
    SyntheticBuilder::cds(
        &sequence,
        PAD_OFFSET + 1,
        PAD_OFFSET + CORE.len() as u64,
        Strand::Plus,
    )
    .build()
}

fn normalize(input: &str, dir: ShuffleDirection) -> String {
    let normalizer =
        Normalizer::with_config(coding(), NormalizeConfig::default().with_direction(dir));
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
        .to_string()
}

/// Assert the exact canonical form, that it re-parses, and that a second pass is
/// a fixed point — the last of which is the property #1207 is about.
fn assert_canonical(input: &str, dir: ShuffleDirection, expected: &str) {
    let once = normalize(input, dir);
    assert_eq!(once, expected, "{input} ({dir:?}) canonical form");
    parse_hgvs(&once).unwrap_or_else(|e| panic!("{once} must re-parse: {e}"));
    let twice = normalize(&once, dir);
    assert_eq!(twice, once, "{input} ({dir:?}) must be a fixed point");
}

// ---------------------------------------------------------------------------
// The 5' boundary at cds_start
// ---------------------------------------------------------------------------

/// `r.1_2dup` duplicates `cc` at the CDS start.
///
/// Re-blessed by #1204 from `r.1delinsccc`. The 5' shuffle still drives the
/// derived insertion to rest immediately 5' of `cds_start`, but the two added
/// `c`s copy `r.1_2` — the 5'-most pair of the `C[4]` tract — so the change is a
/// duplication and `repeated.md` L22 prescribes `dup` for it. A `dup` needs no
/// clamp: both its anchors are inside the CDS, which is exactly the boundary
/// problem the clamp exists to repair. #1207's non-idempotency is fixed either
/// way; the clamp regime it was filed against is pinned below, on the shape that
/// no duplication describes.
#[test]
fn rna_five_prime_boundary_dup_resolves_without_clamping() {
    assert_canonical(
        "NM_TEST.1:r.1_2dup",
        ShuffleDirection::FivePrime,
        "NM_TEST.1:r.1_2dup",
    );
}

/// The `c.` sibling. Parity with `r.` is the contract; that both axes moved
/// together under #1204 is what shows the re-blessing is not `r.`-only drift.
#[test]
fn cds_five_prime_boundary_dup_resolves_without_clamping() {
    assert_canonical(
        "NM_TEST.1:c.1_2dup",
        ShuffleDirection::FivePrime,
        "NM_TEST.1:c.1_2dup",
    );
}

/// Five added `c`s against the four-base `C[4]` tract: there is no five-base `C`
/// tract to copy, so no duplication describes this. It stays an insertion, the 5'
/// shuffle saturates it at `cds_start`, and the clamp applies the identity
/// `delete ref[cds_start]` (= `c`), `insert ccccc ++ c`.
///
/// This is #1207's actual contract — an `r.`-axis insertion saturated at a CDS
/// boundary must be clamped, as the `c.` axis already did — restated on an input
/// the dup promotion cannot absorb. Without a row of this shape #1204 would leave
/// the clamp this file exists for with no direct coverage.
#[test]
fn rna_five_prime_saturated_insertion_clamps_at_cds_start() {
    assert_canonical(
        "NM_TEST.1:r.2_3insccccc",
        ShuffleDirection::FivePrime,
        "NM_TEST.1:r.1delinscccccc",
    );
}

/// The `c.` sibling, which was already correct — this is the parity #1207
/// restores, and the row that shows the `r.` expectation above is not invented.
#[test]
fn cds_five_prime_saturated_insertion_clamps_at_cds_start() {
    assert_canonical(
        "NM_TEST.1:c.2_3insCCCCC",
        ShuffleDirection::FivePrime,
        "NM_TEST.1:c.1delinsCCCCCC",
    );
}

// ---------------------------------------------------------------------------
// The mirror: the 3' boundary at cds_end
// ---------------------------------------------------------------------------

/// `r.13_14dup` duplicates `gg` at the CDS end; the added pair copies the
/// 3'-most two bases of the `G[4]` tract, so this is the `dup` regime.
/// Re-blessed by #1204 from `r.14delinsggg`.
#[test]
fn rna_three_prime_boundary_dup_resolves_without_clamping() {
    assert_canonical(
        "NM_TEST.1:r.13_14dup",
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:r.13_14dup",
    );
}

/// The `c.` sibling for the 3' boundary.
#[test]
fn cds_three_prime_boundary_dup_resolves_without_clamping() {
    assert_canonical(
        "NM_TEST.1:c.13_14dup",
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:c.13_14dup",
    );
}

/// The clamp regime at `cds_end`: five added `g`s against the `G[4]` tract give
/// `delete ref[cds_end]` (= `g`), `insert g ++ ggggg`.
#[test]
fn rna_three_prime_saturated_insertion_clamps_at_cds_end() {
    assert_canonical(
        "NM_TEST.1:r.12_13insggggg",
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:r.14delinsgggggg",
    );
}

/// The `c.` sibling for the 3' clamp.
#[test]
fn cds_three_prime_saturated_insertion_clamps_at_cds_end() {
    assert_canonical(
        "NM_TEST.1:c.12_13insGGGGG",
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:c.14delinsGGGGGG",
    );
}

// ---------------------------------------------------------------------------
// Axis parity, stated as a property rather than per-value
// ---------------------------------------------------------------------------

/// The `r.` result must be the `c.` result re-rendered as RNA (lowercase, `u`
/// for `T`) for both boundaries, both directions, and both regimes. Asserting the
/// relationship rather than only the values keeps this honest when the canonical
/// spelling of a boundary result changes — as #1204 changed the `dup` rows here:
/// the two axes must move together, whatever they move to.
///
/// Note the earlier version of this test required the `c.` answer to contain
/// `delins` and panicked otherwise, which hard-coded the clamp regime into what is
/// supposed to be a shape-agnostic parity property. It now finds whichever edit
/// keyword the answer carries, so a `dup` row is compared the same way a `delins`
/// row is.
#[test]
fn rna_and_cds_axes_agree_at_both_cds_boundaries() {
    /// Re-render a `c.` description as its `r.` sibling: swap the axis, then
    /// lowercase from the edit keyword on and map `t` to `u`. Positions are left
    /// alone, so this is safe for keywords that carry no bases (`dup`, `del`).
    fn as_rna(cds: &str) -> String {
        let swapped = cds.replace(":c.", ":r.");
        let split = ["delins", "dup", "ins", "del", "inv"]
            .iter()
            .filter_map(|kw| swapped.find(kw))
            .min()
            .unwrap_or_else(|| panic!("no edit keyword in {cds}"));
        let (prefix, edit) = swapped.split_at(split);
        format!("{prefix}{}", edit.to_lowercase().replace('t', "u"))
    }

    for (dir, input) in [
        // the dup regime
        (ShuffleDirection::FivePrime, "1_2dup"),
        (ShuffleDirection::ThreePrime, "13_14dup"),
        // the clamp regime
        (ShuffleDirection::FivePrime, "2_3insCCCCC"),
        (ShuffleDirection::ThreePrime, "12_13insGGGGG"),
    ] {
        let rna = normalize(&format!("NM_TEST.1:r.{}", input.to_lowercase()), dir);
        let cds = normalize(&format!("NM_TEST.1:c.{input}"), dir);
        assert_eq!(
            rna,
            as_rna(&cds),
            "r. and c. must agree at the CDS boundary ({input}, {dir:?})",
        );
    }
}

// ---------------------------------------------------------------------------
// Non-regression: the clamp must not fire away from the boundaries
// ---------------------------------------------------------------------------

/// An interior dup must be untouched by the clamp. `r.5_6dup` sits in the `A`×6
/// tract (`r.5`..`r.10`), so the 3' shuffle rests the derived insertion just past
/// the tract; the codon-frame gate refuses repeat notation there (1-nt unit,
/// wholly inside the CDS, #1192), and the resting place is nowhere near `cds_end`
/// (`r.14`), so the boundary clamp must not fire.
///
/// Re-blessed by #1204 from `r.10_11insaa` / `c.10_11insAA`: this row never had
/// anything to do with the clamp, and its `ins` expectation was itself an instance
/// of the #1204 defect — the added `aa` copies `r.9_10`, the 3'-most pair of the
/// tract, so the gated fallback is a `dup`. The clamp still does not fire, which is
/// what this test is for.
///
/// Both axes are pinned, so this also shows the interior case keeps `r.`/`c.`
/// parity — the property, not just "no delins appeared".
#[test]
fn interior_dup_is_untouched_by_the_clamp() {
    assert_canonical(
        "NM_TEST.1:r.5_6dup",
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:r.9_10dup",
    );
    assert_canonical(
        "NM_TEST.1:c.5_6dup",
        ShuffleDirection::ThreePrime,
        "NM_TEST.1:c.9_10dup",
    );
}

/// The #383 / #387 carve-out must survive on `r.` too: a `Delins` **input**
/// whose shared-affix trim pushes the residual past a CDS boundary keeps its own
/// form instead of being clamped.
///
/// Fixture mirrors `issue_387_canon_cds_end_clamp`'s exactly (`c.1..c.5` =
/// `ATGAA`, `cds_end` = `c.5`, first 3'UTR base `G`), where
/// `three_prime_delins_at_cds_end_suppresses_rewrite` pins
/// `c.5delinsAC -> c.5delinsAC`. Without the restore arm the new clamp would
/// rewrite the reduced insertion into a boundary `delins` on `r.` only — trading
/// the divergence this PR fixes for a different one.
#[test]
fn rna_delins_input_at_cds_end_keeps_its_own_form_like_the_c_axis() {
    use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Transcript};
    let mut seq = String::from("GGG"); // 5'UTR
    seq.push_str("ATGAA"); // c.1..c.5
    while seq.len() < 203 {
        seq.push('G'); // 3'UTR
    }
    let len = seq.len() as u64;
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_TEST387.1".to_string(),
        Some("TEST387".to_string()),
        Strand::Plus,
        seq,
        Some(4),
        Some(8),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    ));
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let run = |input: &str| {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        normalizer
            .normalize(&v)
            .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
            .to_string()
    };
    assert_eq!(
        run("NM_TEST387.1:c.5delinsAC"),
        "NM_TEST387.1:c.5delinsAC",
        "premise: the c. axis suppresses the rewrite (#387)",
    );
    assert_eq!(
        run("NM_TEST387.1:r.5delinsac"),
        "NM_TEST387.1:r.5delinsac",
        "r. must suppress it too, not clamp to a boundary delins",
    );
}

/// A non-coding transcript has no CDS, so the new clamp has nothing to key on.
/// This guards against the new gate accidentally firing on a `None` CDS.
///
/// The mechanism is worth stating, because it is not "the clamp declined" — it is
/// that the clamp is never even consulted. With no CDS, `is_coding` is `false`, so
/// the codon-frame gate does not fire (#1192: nothing to keep in frame) and the
/// dup converts straight to repeat notation. No saturated insertion is ever
/// produced, so neither this clamp nor #1202's transcript-bound one sees an
/// `Insertion` to rewrite.
///
/// Derived: the fixture is the bare `CORE` (14 bases, unpadded), whose leading
/// `C` tract is `r.1_4` at 4 copies; duplicating two of them gives 6.
#[test]
fn noncoding_transcript_is_unaffected() {
    let provider = SyntheticBuilder::noncoding(CORE, Strand::Plus).build();
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
    );
    let run = |input: &str| {
        let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        normalizer
            .normalize(&v)
            .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
            .to_string()
    };
    let once = run("NR_TEST.1:r.1_2dup");
    assert_eq!(
        once, "NR_TEST.1:r.1_4c[6]",
        "a non-coding r. transcript keeps repeat notation: no reading frame to gate on",
    );
    assert_eq!(run(&once), once, "and must be a fixed point");
}
