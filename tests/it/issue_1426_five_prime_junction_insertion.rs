//! The 5' mirror of #1426: a 3'UTR insertion shifting toward the CDS must
//! finish in one pass.
//!
//! ```text
//! c.*1_*2insCTA  ->  c.18_*1insACT  ->  c.17_18insTAC
//! ```
//!
//! The cause is not the 3' half's — that was the #918 relaxation asking an
//! insertion for both its flanks. Here it is the axis clamp itself: its left
//! bound is derived from the region the variant currently sits in, and for a
//! zero-width insertion that region flips as the insertion moves, so the bound
//! that stops one pass is one the next no longer reads.
//!
//! ```text
//! c.*1_*2insCTA  tx=(19,20)  start_axis=ThreeUtr  clamped.left=18  -> stops at c.18_*1
//! c.18_*1insACT  tx=(18,19)  start_axis=Cds       clamped.left=0   -> reaches c.17_18
//! ```
//!
//! So the 5' shuffle already crossed `cds_end` before this fix — one base per
//! pass. The fix makes it do so in one, which is why this is a defect rather
//! than a policy change about whether it may cross at all.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::normalize::{NormalizeConfig, ShuffleDirection};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// 20 bases, CDS `1..=18`; `c.18` is the last coding base, `c.*1` the first
/// 3'UTR base. From the transcript-axis sweep corpus, where this was found.
const CORE: &str = "TAATTTATTAATATATATAT";

fn normalizer() -> Normalizer<impl ferro_hgvs::reference::ReferenceProvider> {
    Normalizer::with_config(
        SyntheticBuilder::cds(CORE, 1, 18, Strand::Plus).build(),
        NormalizeConfig::default()
            .with_direction(ShuffleDirection::FivePrime)
            .allow_crossing_boundaries(),
    )
}

fn normalize_twice(input: &str) -> (String, String) {
    let once = normalize_once(input);
    (once.clone(), normalize_once(&once))
}

fn normalize_once(input: &str) -> String {
    normalizer()
        .normalize(&parse_hgvs(input).expect("parse input"))
        .expect("normalize")
        .to_string()
}

/// The same input through `normalize_with_diagnostics`.
///
/// ferro has **two** public normalization exits and they are not the same code
/// path: `normalize()` reaches `normalize_core_checked`, while
/// `normalize_with_diagnostics` historically reached `normalize_core_canonical`
/// directly and so returned results no oracle had inspected (#1382). A bound
/// relaxation like this one has to hold on both, or the fix is real for one
/// caller and absent for the other.
///
/// Mirrors `normalize_once_with_diagnostics` in
/// `issue_1426_junction_insertion_settles_in_one_pass` — the 3' half of this
/// same fix covers both exits, and the two arms should not differ in what they
/// pin any more than they should differ in their predicate.
fn normalize_once_with_diagnostics(input: &str) -> String {
    normalizer()
        .normalize_with_diagnostics(&parse_hgvs(input).expect("parse input"))
        .expect("normalize")
        .result
        .to_string()
}

/// Both exits must agree, on the pinned answer and with each other.
#[test]
fn both_public_exits_settle_on_the_same_form() {
    for input in ["NM_TEST.1:c.*1_*2insCTA", "NM_TEST.1:c.18_*1insACT"] {
        let plain = normalize_once(input);
        let with_diagnostics = normalize_once_with_diagnostics(input);
        assert_eq!(
            plain, with_diagnostics,
            "`{input}` settles differently through the two public exits"
        );
        assert_eq!(
            with_diagnostics, "NM_TEST.1:c.17_18insTAC",
            "`{input}` must reach the 5'-most form through the diagnostics exit too"
        );
    }
}

/// The 3'UTR-resident insertion reaches its 5'-most position immediately.
#[test]
fn a_three_prime_utr_insertion_settles_in_one_pass() {
    let (once, twice) = normalize_twice("NM_TEST.1:c.*1_*2insCTA");
    assert_eq!(once, "NM_TEST.1:c.17_18insTAC");
    assert_eq!(twice, once, "and is stable from there");
}

/// The intermediate the old code stopped on now settles to the same place.
///
/// Stated separately: a fix that merely made `c.18_*1insACT` stable would
/// satisfy idempotency for that input while leaving the two entry points
/// disagreeing about where the same edit belongs.
#[test]
fn the_junction_intermediate_reaches_the_same_answer() {
    let (once, twice) = normalize_twice("NM_TEST.1:c.18_*1insACT");
    assert_eq!(once, "NM_TEST.1:c.17_18insTAC");
    assert_eq!(twice, once);
}
