//! The sequence-first derivation must shift its pieces in the **configured**
//! direction (#1429).
//!
//! `canonicalize_from_sequence` never received the caller's shuffle direction,
//! and its piece shift was hardcoded to `ThreePrime`. A caller asking for a 5'
//! shuffle therefore got a *3'-shifted* merged allele, which the per-member
//! pipeline then moved on the next pass:
//!
//! ```text
//! c.[16dup;18_19insC]  ->  c.*1_*2insCT  ->  c.18_*1insTC
//! ```
//!
//! Two things were wrong there and the idempotency failure is the lesser one:
//! the first-pass answer was simply not the 5' answer. So these tests assert the
//! *value*, not merely that it is a fixed point — a derivation that stably
//! returned the 3' form would satisfy idempotency while still ignoring the
//! direction it was asked for.
//!
//! Found by the transcript-axis sweep (#1400), the fifth distinct defect behind
//! its single `#[ignore]`.
//!
//! The case below is deliberately **not** the one from the sweep. That input
//! lands on the CDS/3'UTR junction, so its idempotency also depends on #1427 and
//! #1428 (the two halves of #1426) and a guard here would be red for reasons
//! that are not this change. This one sits nine bases inside the CDS, where the
//! only thing deciding the answer is which way the derivation shifts. Verified
//! end-to-end separately: with all five fixes applied the full 829,440-case
//! transcript sweep passes under all three oracles.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::normalize::{NormalizeConfig, ShuffleDirection};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// 20 bases, CDS `1..=18`. The nine-`T` run at the 5' end is the tract the
/// merged duplication travels, and it is long enough that the two directions
/// land seven bases apart — far too much to read as a rounding difference.
const CORE: &str = "TTTTTTTTTAATATATTTTA";

fn normalize_once(input: &str, dir: ShuffleDirection) -> String {
    let normalizer = Normalizer::with_config(
        SyntheticBuilder::cds(CORE, 1, 18, Strand::Plus).build(),
        NormalizeConfig::default()
            .with_direction(dir)
            .allow_crossing_boundaries(),
    );
    normalizer
        .normalize(&parse_hgvs(input).expect("parse input"))
        .expect("normalize")
        .to_string()
}

/// A 5' request produces the 5'-most merged form.
///
/// Before the fix this returned `c.8_9dup` — the 3' answer — because the
/// derivation shifted its pieces 3' whatever it was asked for.
#[test]
fn a_five_prime_request_gets_the_five_prime_merged_form() {
    let once = normalize_once("NM_TEST.1:c.[1dup;2dup]", ShuffleDirection::FivePrime);
    assert_eq!(once, "NM_TEST.1:c.1_2dup");
    assert_eq!(
        normalize_once(&once, ShuffleDirection::FivePrime),
        once,
        "and is stable there"
    );
}

/// The 3' direction is unchanged: it passed `ThreePrime` before and still does.
///
/// Pinned so the new parameter cannot be threaded through backwards — an
/// inversion would swap these two answers and nothing else in the suite would
/// notice.
#[test]
fn a_three_prime_request_is_unaffected() {
    let once = normalize_once("NM_TEST.1:c.[1dup;2dup]", ShuffleDirection::ThreePrime);
    assert_eq!(once, "NM_TEST.1:c.8_9dup");
    assert_eq!(
        normalize_once(&once, ShuffleDirection::ThreePrime),
        once,
        "and is stable there"
    );
}
