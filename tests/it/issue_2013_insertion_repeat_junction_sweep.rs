//! #2013: an insertion, an adjacent `dup` and a nearby `del` normalized to a
//! description denoting a *different* sequence.
//!
//! Against a 600-base GRCh38 window served as `NC_TESTWIN.1`, the cis allele
//! `g.[213_214insGAGA;214dup;218_220del]` normalized to
//! `g.[214dup;217_219delinsAGA;220_221insT]`, which `parse_hgvs` accepts — so
//! the re-parse oracle never saw it — but which denotes a different sequence
//! than the input.
//!
//! # The mechanism, so the guard is understood rather than merely pinned
//!
//! The three members are normalized per member, then a chain of sibling-clamp
//! passes repairs any interaction the standalone shifts created. Standalone,
//! `213_214insGAGA` 3'-shifts into the *repeat* `215_216AG[3]`, carrying its
//! inserted `GAGA` clear over the `dup` at 214. That crossing changes what the
//! allele denotes, and it is exactly what `clamp_sibling_crossing_junctions`
//! exists to bound — but a repeat states no junction, so that clamp skipped it;
//! `clamp_sibling_crossing_shifts` declines a member that grew; and the tract
//! (215-216) overlaps no sibling, so `demote_repeats_spanning_siblings`'s tract
//! test did not fire either. It was the un-owned `before = ins`, `after =
//! repeat` cell of the ownership table on `reduced_before_junction`.
//!
//! The fix demotes such a repeat back to the insertion at its 3' end when the
//! insertion's payload junction swept across a shift-blocking sibling, so the
//! junction clamp — immediately after — pulls it back.
//!
//! # Why the oracle is not the normalizer
//!
//! The denoted sequence is computed with [`apply_with`], which converts each
//! member to SPDI and splices the reference itself. Checking the normalizer's
//! output against the normalizer would pass for any self-consistent defect;
//! this asks whether the output *means* the same thing as the input.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer};

use crate::common::cis_apply_oracle::apply_with;
use crate::common::hg38_window::{local_desc, provider, HG38_WINDOW};

/// Normalize a local-coordinate `g.` allele body against the window in strict
/// mode, returning the rendered output.
fn normalize(body: &str) -> String {
    let input = local_desc(body);
    let variant: HgvsVariant = parse_hgvs(&input).expect("input parses");
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    .expect("input normalizes")
    .to_string()
}

/// The reported allele normalizes to a description denoting the same sequence,
/// and that output re-parses.
#[test]
fn the_reported_allele_preserves_its_denoted_sequence() {
    let body = "[213_214insGAGA;214dup;218_220del]";
    let out = normalize(body);

    assert!(parse_hgvs(&out).is_ok(), "output {out} does not re-parse");

    let p = provider();
    let before = apply_with(&p, HG38_WINDOW, &local_desc(body)).expect("input denotes a sequence");
    let after = apply_with(&p, HG38_WINDOW, &out).expect("output denotes a sequence");
    assert_eq!(
        before, after,
        "normalizing {body} to {out} changed the denoted sequence",
    );
}

/// Member order is not part of what the allele denotes: the reverse spelling
/// normalizes to the same sequence-preserving output. Order-independence is the
/// property #2013's fall-through most easily masks, since the buggy pass keyed
/// off the order the standalone shifts happened to leave the members in.
#[test]
fn the_reported_allele_is_order_independent() {
    let forward = normalize("[213_214insGAGA;214dup;218_220del]");
    let reversed = normalize("[218_220del;214dup;213_214insGAGA]");
    assert_eq!(
        forward, reversed,
        "member order leaked into the normalized form",
    );
}
