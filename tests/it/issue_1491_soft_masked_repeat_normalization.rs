//! #1491 — the canonical form of a repeat must not depend on whether the
//! reference happens to be soft-masked.
//!
//! `normalize::rules::count_tandem_repeats` compared reference bytes to the
//! description's unit with raw slice equality. A reference FASTA is routinely
//! soft-masked and a description is conventionally written upper-case, so on a
//! masked tract the search found no run — and because the function returns
//! *coordinates* rather than a verdict, that did not surface as an error. It
//! surfaced as a different canonical form:
//!
//! | input | uppercase reference | soft-masked reference |
//! |---|---|---|
//! | `g.259_267CAG[2]` | `g.265_267del` | `g.259_267CAG[2]` |
//! | `g.259CAG[5]` | `g.259_267CAG[5]` | `g.259CAG[5]` |
//! | `g.259A[7]` | `g.259_263A[7]` | `g.259A[7]` |
//! | `g.259AT[6]` | `g.259_266AT[6]` | `g.259AT[6]` |
//!
//! Two shapes. A **contraction fails to lower**: `CAG[2]` over a 3-copy tract
//! is a deletion of one copy and becomes `g.265_267del` on an unmasked
//! reference, but stayed an unlowered repeat on the masked twin. And a
//! **single-position anchor is not widened to its tract**, which leaves the
//! bare anchor standing where the whole-range form is the canonical answer.
//!
//! Measured 4 of 15 repeat inputs moved with case alone; the other 11 — plain
//! `del`/`dup`/`ins`, and explicit-range repeats whose count is not a
//! contraction — agreed, which is the control showing this is the tract search
//! and not a general masking problem.
//!
//! No oracle could see it. Every masked output above is well-formed, in bounds
//! and a fixed point, so `FERRO_ASSERT_REPARSE`, `FERRO_ASSERT_IN_BOUNDS` and
//! `FERRO_ASSERT_IDEMPOTENT` all pass on it. The violated property is *parity
//! between two references differing only in case*, which nothing asserted —
//! the same blind spot #1318 and #1250 were each filed against for a different
//! comparison site.
//!
//! Reference FASTAs are routinely soft-masked and masked bytes reach
//! normalization from ordinary input (`issue_1318_soft_masked_delins.rs`,
//! `issue_1235_cis_allele_confluence.rs::soft_masked_reference_yields_the_same_canonical_form`).

use std::io::Write;

use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};

/// 256 `N`s, so core position 1 is `g.257` and each tract below starts at 259.
/// `N` matches no unit used here, so the flanks cannot look like part of a run.
const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

fn padded(core: &str) -> String {
    format!("{PAD}{core}{PAD}")
}

fn provider_for(sequence: &str) -> JsonProvider {
    let n = sequence.len() as u64;
    let json = serde_json::json!({
        "version": "1.0",
        "transcripts": [{
            "id": "TEMPLATE-gene.1",
            "chromosome": "TEMPLATE",
            "strand": "+",
            "sequence": sequence,
            "cds_start": 1,
            "cds_end": n - (n % 3),
            "genomic_start": 1,
            "genomic_end": n,
            "protein_id": "TEMPLATE-gene.1",
            "exons": [{"number": 1, "start": 1, "end": n,
                       "genomic_start": 1, "genomic_end": n}],
        }],
        "genomic_sequences": {"TEMPLATE": sequence},
    });
    let mut f = tempfile::NamedTempFile::new().expect("tempfile");
    write!(f, "{json}").expect("write json");
    JsonProvider::from_json(f.path()).expect("load reference")
}

fn normalize_with(sequence: &str, input: &str) -> String {
    let normalizer = Normalizer::new(provider_for(sequence));
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse `{input}`: {e}"));
    normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize `{input}`: {e}"))
        .to_string()
}

/// The four inputs whose canonical form moved with case, each with the answer
/// pinned as a literal.
///
/// Both assertions per row are needed, and the pinned literal is the primary
/// one: a parity-only check is satisfied by two *equally wrong* answers, and
/// the whole defect class here is one comparison folding case while its
/// siblings do. The parity assertion is the supplement — it says the reason the
/// answer is right is not "the uppercase path happens to agree today".
#[test]
fn a_repeats_canonical_form_does_not_depend_on_soft_masking() {
    for (core, input, expected) in [
        // A contraction: `CAG[2]` over a 3-copy tract deletes one copy.
        (
            "GGCAGCAGCAGGG",
            "TEMPLATE:g.259_267CAG[2]",
            "TEMPLATE:g.265_267del",
        ),
        // Single-position anchors, widened to the run they name.
        (
            "GGCAGCAGCAGGG",
            "TEMPLATE:g.259CAG[5]",
            "TEMPLATE:g.259_267CAG[5]",
        ),
        ("GGAAAAAGG", "TEMPLATE:g.259A[7]", "TEMPLATE:g.259_263A[7]"),
        (
            "GGATATATATGG",
            "TEMPLATE:g.259AT[6]",
            "TEMPLATE:g.259_266AT[6]",
        ),
    ] {
        let masked = normalize_with(&padded(&core.to_ascii_lowercase()), input);
        assert_eq!(
            masked, expected,
            "canonical form of `{input}` when soft-masked"
        );
        let upper = normalize_with(&padded(core), input);
        assert_eq!(
            masked, upper,
            "soft-masking must not change the canonical form of `{input}`"
        );
    }
}

/// The control: the eleven inputs that already agreed must keep agreeing, and
/// keep their current answers.
///
/// This is what says the fix is the tract search rather than a blanket
/// case-fold that happened to move everything. `del`/`dup`/`ins` reach their
/// own (already folded) comparison sites, and an explicit-range repeat whose
/// count is not a contraction never needs the tract's extent.
#[test]
fn the_shapes_that_already_agreed_are_untouched() {
    // The expected form is pinned alongside the parity check on purpose. Parity
    // alone (`masked == upper`) is satisfied by *both* sides moving together,
    // which is exactly the regression this control exists to rule out — the
    // "blanket case-fold that happened to move everything" in the doc above
    // would keep every row of this table equal to itself while changing all of
    // them. Only the literal says the answers did not move.
    for (core, input, expected) in [
        (
            "GGCAGCAGCAGGG",
            "TEMPLATE:g.259_267CAG[5]",
            "TEMPLATE:g.259_267CAG[5]",
        ),
        (
            "GGCAGCAGCAGGG",
            "TEMPLATE:g.259_267del",
            "TEMPLATE:g.259_267del",
        ),
        (
            "GGCAGCAGCAGGG",
            "TEMPLATE:g.259_261dup",
            "TEMPLATE:g.265_267dup",
        ),
        // Converges on the same `dup` as the row above, from an `ins` spelling.
        (
            "GGCAGCAGCAGGG",
            "TEMPLATE:g.267_268insCAG",
            "TEMPLATE:g.265_267dup",
        ),
        (
            "GGAAAAAGG",
            "TEMPLATE:g.259_263A[7]",
            "TEMPLATE:g.259_263A[7]",
        ),
        (
            "GGAAAAAGG",
            "TEMPLATE:g.259_263A[2]",
            "TEMPLATE:g.259_263A[2]",
        ),
        ("GGAAAAAGG", "TEMPLATE:g.259dup", "TEMPLATE:g.263dup"),
        ("GGAAAAAGG", "TEMPLATE:g.263_264insA", "TEMPLATE:g.263dup"),
        (
            "GGATATATATGG",
            "TEMPLATE:g.259_266AT[6]",
            "TEMPLATE:g.259_266AT[6]",
        ),
        (
            "GGATATATATGG",
            "TEMPLATE:g.259_266AT[2]",
            "TEMPLATE:g.259_266AT[2]",
        ),
        (
            "GGATATATATGG",
            "TEMPLATE:g.259_260dup",
            "TEMPLATE:g.265_266dup",
        ),
    ] {
        let masked = normalize_with(&padded(&core.to_ascii_lowercase()), input);
        assert_eq!(masked, expected, "canonical form of `{input}` when masked");
        let upper = normalize_with(&padded(core), input);
        assert_eq!(
            masked, upper,
            "`{input}` agreed across case before the fix and must still"
        );
    }
}

/// The discriminating case: folding must not invent a tract out of bases that
/// are not the declared unit.
///
/// The over-general version of this fix — comparing only lengths, or treating
/// any base as matching — would widen `g.259A[5]` over `ACGTA` into a run that
/// is not there. `N` is the sharpest probe because it is neither the unit nor
/// its lower-case form, and the pad is built from it.
#[test]
fn folding_does_not_invent_a_tract() {
    // `A` appears at 259 and 263 but the run is one base long in both cases, so
    // the one-copy tract at 259 is the whole run and the repeat stands as
    // written. (Measured, not inferred: an earlier draft of this test guessed
    // `g.259_260insAA` and was wrong.)
    let core = "GGACGTAGG";
    for seq in [padded(core), padded(&core.to_ascii_lowercase())] {
        assert_eq!(
            normalize_with(&seq, "TEMPLATE:g.259A[3]"),
            "TEMPLATE:g.259A[3]",
            "a single `A` must not become a longer tract"
        );
    }
}
