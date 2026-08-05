//! An insertion resting on the CDS/3'UTR junction must complete its 3'-shift on
//! the **first** pass (#1426).
//!
//! This is #1209's defect one shape over, and its file records the identical
//! three-step alternation:
//!
//! ```text
//! #1209:  c.25_26insGAT  ->  c.26delinsGATG  ->  c.*1_*2insTGA
//! #1426:  c.18_*1insACT  ->  c.18delinsTACT  ->  c.*1_*2insCTA
//! ```
//!
//! #1209 fixed it by admitting `Insertion` to the #918 bound relaxation. That
//! relaxation also asks both endpoints to be CDS-resident, which an insertion
//! *on* the junction can never satisfy: its endpoints are the two flanks of a
//! gap, so `end_axis` reads `ThreeUtr` while it covers no 3'UTR base at all. So
//! the first pass was barred, stopped on the boundary and was rewritten by the
//! #387 clamp; the resulting `Delins` has both flanks on `c.18`, so the second
//! pass was admitted and shifted. Two passes, two branches.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::normalize::{NormalizeConfig, ShuffleDirection};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// 20 bases, CDS `1..=18`, so `c.18` is the last coding base and `c.*1` the
/// first 3'UTR base. Drawn from the transcript-axis sweep corpus, which is where
/// this was found.
const CORE: &str = "TAATTTATTAATATATATAT";

fn normalize_twice(input: &str) -> (String, String) {
    let once = normalize_once(input);
    (once.clone(), normalize_once(&once))
}

fn normalizer() -> Normalizer<impl ferro_hgvs::reference::ReferenceProvider> {
    Normalizer::with_config(
        SyntheticBuilder::cds(CORE, 1, 18, Strand::Plus).build(),
        NormalizeConfig::default()
            .with_direction(ShuffleDirection::ThreePrime)
            .allow_crossing_boundaries(),
    )
}

fn normalize_once(input: &str) -> String {
    normalizer()
        .normalize(&parse_hgvs(input).expect("parse input"))
        .expect("normalize")
        .to_string()
}

/// The same input through `normalize_with_diagnostics`.
///
/// ferro has **two** public normalization exits, and they are not the same code
/// path: `normalize()` reaches `normalize_core_checked`, while
/// `normalize_with_diagnostics` historically reached `normalize_core_canonical`
/// directly and so returned results no oracle had inspected (#1382). A bound
/// relaxation like this one has to hold on both, or the fix is real for one
/// caller and absent for the other.
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
    for input in ["NM_TEST.1:c.18_*1insACT", "NM_TEST.1:c.*1_*2insCTA"] {
        let plain = normalize_once(input);
        let with_diagnostics = normalize_once_with_diagnostics(input);
        assert_eq!(
            plain, with_diagnostics,
            "`{input}` settles differently through the two public exits"
        );
    }
}

/// The junction insertion shifts the whole way immediately.
#[test]
fn a_junction_insertion_settles_in_one_pass() {
    let (once, twice) = normalize_twice("NM_TEST.1:c.18_*1insACT");
    assert_eq!(once, "NM_TEST.1:c.*1_*2insCTA");
    assert_eq!(twice, once, "and is stable from there");
}

/// The clamped form is no longer produced, and converges on the shifted one.
///
/// The mirror of `issue_1209_cds_end_insertion_shift::
/// the_old_clamped_output_converges_on_the_shifted_form`, and stated separately
/// because it is the assertion that would fail if someone "fixed" this by
/// widening the #387 clamp instead — that makes `c.18delins…` a fixed point,
/// which is the behaviour #1209 deliberately removed.
#[test]
fn the_clamped_form_converges_on_the_shifted_form() {
    let (once, twice) = normalize_twice("NM_TEST.1:c.18delinsTACT");
    assert_eq!(once, "NM_TEST.1:c.*1_*2insCTA");
    assert_eq!(twice, once, "and is stable from there");
}

/// A 3'UTR-resident insertion is untouched — it was already stable.
#[test]
fn a_three_prime_utr_insertion_is_unchanged() {
    let (once, twice) = normalize_twice("NM_TEST.1:c.*1_*2insCTA");
    assert_eq!(once, "NM_TEST.1:c.*1_*2insCTA");
    assert_eq!(twice, once);
}
