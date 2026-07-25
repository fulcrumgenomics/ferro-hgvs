//! Idempotency of **repeat-notation (`unit[N]`) inputs**.
//!
//! The idempotency proptest generates only `del`/`dup`/`ins`/`delins` inputs, so
//! `unit[N]` spelled *as input* was never fuzzed. `normalize_repeat` must be a
//! fixed point for these too: `norm(norm(x)) == norm(x)`.
//!
//! Core `TCAGCAGCAGCAGTT` places exactly 4 clean copies of `CAG` at g.258..269,
//! bounded by `T` on both sides so neither the 5' nor the 3' phase walk can leak
//! into the `ACGT…` padding (a `G`-prefixed or `C`-suffixed core would continue
//! the CAG rotation into the pad and measure a harness artifact instead).
//!
//! Three count regimes exercise three of the `normalize_repeat` result arms:
//!   * `[5]` = ref_count + 1  -> Duplication arm
//!   * `[3]` = ref_count - 1  -> Deletion arm
//!   * `[2]` = ref_count - 2  -> Repeat (contraction-with-survivors, B2) arm
//!
//! The fourth — the codon-gated `Insertion` arm, reachable only in `c.` context
//! with a non-codon-aligned unit — is covered separately at the bottom of the
//! file, on its own CDS fixture.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// 4 clean copies of `CAG` at g.258..269, `T`-bounded on both flanks.
const CORE: &str = "TCAGCAGCAGCAGTT";

fn norm(input: &str, dir: ShuffleDirection) -> String {
    Normalizer::with_config(
        SyntheticBuilder::genomic(CORE).build(),
        NormalizeConfig::default().with_direction(dir),
    )
    .normalize(&parse_hgvs(input).expect("parse"))
    .expect("normalize")
    .to_string()
}

fn norm_cds(core: &str, input: &str, dir: ShuffleDirection) -> String {
    Normalizer::with_config(
        SyntheticBuilder::cds(core, 1, core.len() as u64, Strand::Plus).build(),
        NormalizeConfig::default().with_direction(dir),
    )
    .normalize(&parse_hgvs(input).expect("parse"))
    .expect("normalize")
    .to_string()
}

/// `norm(x) == expected` AND `norm(norm(x)) == norm(x)`.
///
/// Pinning the exact canonical form matters as much as the fixed point: a
/// normalizer that returned its input untouched would satisfy idempotency alone,
/// so a fixed-point-only assertion cannot tell a working direction rule from a
/// no-op. `expected` also encodes *which* copy each direction names, which is the
/// behavior this PR changes.
fn assert_normalizes_to(input: &str, dir: ShuffleDirection, expected: &str) {
    let once = norm(input, dir);
    assert_eq!(
        once, expected,
        "WRONG CANONICAL FORM ({dir:?})\n  input   ={input}\n  got     ={once}\n  expected={expected}"
    );
    let twice = norm(&once, dir);
    assert_eq!(
        once, twice,
        "NOT IDEMPOTENT ({dir:?})\n  input={input}\n  once ={once}\n  twice={twice}"
    );
}

// The tract is 4 copies of `CAG` at g.258..269 — copies at 258_260, 261_263,
// 264_266, 267_269. Every expectation below follows from that layout plus the
// direction rule: `ThreePrime` names the 3'-most copy/copies, `FivePrime` the
// 5'-most. The `[2]` (contraction-with-survivors) case is the exception — the
// B2 `Repeat` arm emits a window spanning the WHOLE tract, so it is
// direction-invariant by construction and both directions agree.

#[test]
fn three_prime_repeat_expansion_names_the_last_copy() {
    // +1 copy -> Duplication arm; 3' names the 3'-most copy.
    assert_normalizes_to(
        "NC_TEST.1:g.258CAG[5]",
        ShuffleDirection::ThreePrime,
        "NC_TEST.1:g.267_269dup",
    );
}

#[test]
fn three_prime_repeat_single_contraction_names_the_last_copy() {
    // -1 copy -> Deletion arm; 3' names the 3'-most copy.
    assert_normalizes_to(
        "NC_TEST.1:g.258CAG[3]",
        ShuffleDirection::ThreePrime,
        "NC_TEST.1:g.267_269del",
    );
}

#[test]
fn three_prime_repeat_double_contraction_spans_the_whole_tract() {
    // -2 copies with survivors -> B2 Repeat arm; window covers the full tract.
    assert_normalizes_to(
        "NC_TEST.1:g.258CAG[2]",
        ShuffleDirection::ThreePrime,
        "NC_TEST.1:g.258_269CAG[2]",
    );
}

#[test]
fn five_prime_repeat_expansion_names_the_first_copy() {
    // Same haplotype as the 3' case, but 5' names the 5'-most copy.
    assert_normalizes_to(
        "NC_TEST.1:g.258CAG[5]",
        ShuffleDirection::FivePrime,
        "NC_TEST.1:g.258_260dup",
    );
}

#[test]
fn five_prime_repeat_single_contraction_names_the_first_copy() {
    assert_normalizes_to(
        "NC_TEST.1:g.258CAG[3]",
        ShuffleDirection::FivePrime,
        "NC_TEST.1:g.258_260del",
    );
}

#[test]
fn five_prime_repeat_double_contraction_spans_the_whole_tract() {
    // Direction-invariant: identical to the 3' expectation above.
    assert_normalizes_to(
        "NC_TEST.1:g.258CAG[2]",
        ShuffleDirection::FivePrime,
        "NC_TEST.1:g.258_269CAG[2]",
    );
}

/// The codon-frame gate (`repeated.md`: `c.` repeat notation requires
/// `unit_len % 3 == 0`) routes a >=2-copy expansion of a non-codon-aligned unit
/// to a literal `ins`. That insertion point was named against the tract's **3'**
/// flank regardless of direction, so under 5' it re-shuffled to the 5' flank
/// (`c.6_7insAA` -> `c.3_4insAA`). Found by the repeat-input fuzz.
///
/// Core `CCCAAACCC`: an `AAA` tract at c.4..6 (unit `A`, len 1, gate blocks),
/// expanded to 5 copies = +2 copies.
#[test]
fn five_prime_codon_gated_repeat_expansion_is_idempotent() {
    let once = norm_cds(
        "CCCAAACCC",
        "NM_TEST.1:c.4A[5]",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(once, "NM_TEST.1:c.3_4insAA", "must name the 5' flank");
    assert_eq!(
        norm_cds("CCCAAACCC", &once, ShuffleDirection::FivePrime),
        once,
        "codon-gated 5' repeat expansion must be a fixed point"
    );
}

/// 3' is unchanged: it names the tract's 3' flank and is idempotent.
#[test]
fn three_prime_codon_gated_repeat_expansion_unchanged() {
    let once = norm_cds(
        "CCCAAACCC",
        "NM_TEST.1:c.4A[5]",
        ShuffleDirection::ThreePrime,
    );
    assert_eq!(once, "NM_TEST.1:c.6_7insAA");
    assert_eq!(
        norm_cds("CCCAAACCC", &once, ShuffleDirection::ThreePrime),
        once,
        "3' must remain idempotent"
    );
}
