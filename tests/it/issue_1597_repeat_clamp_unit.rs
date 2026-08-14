//! #1597: a `Repeat` member must never be clamped onto bases that are not its unit.
//!
//! `clamp_sibling_crossing_shifts` pulls a member back when its standalone
//! 3'-shift carried it over a sibling's bases. The pull is a plain translation
//! of the member's coordinates, and `translate_member`'s `landed` check verifies
//! only that the span arrived where it was aimed — never what the reference
//! holds there. For every other edit type that is enough, because the span is
//! the only claim the spelling makes.
//!
//! A repeat is the exception: `g.269A[3]` asserts that the bases under its span
//! are a tandem tract of `A`. Slide it one base 5' onto `g.268` — a `C` on this
//! contig — and the description is no longer about the sequence it sits on. It
//! still parses, still occupies a legal coordinate and is still a fixed point,
//! so `FERRO_ASSERT_REPARSE`, `FERRO_ASSERT_IN_BOUNDS` and
//! `FERRO_ASSERT_IDEMPOTENT` all pass on the result; only the bases disagree.
//!
//! Measured on `origin/main` at `439617c2`:
//!
//! ```text
//! in                    NC_TEST.1:g.[267_268insCA;268C>A;269_270insAA]
//! before the clamp      [NC_TEST.1:g.269A[3];NC_TEST.1:g.269A[3]]
//! after  the clamp      [NC_TEST.1:g.268A[3];NC_TEST.1:g.269A[3]]   <-- g.268 is C
//! emitted               NC_TEST.1:g.[268A[3];269A[3]]  (members overlap; denotes nothing)
//! ```
//!
//! `junction_of` reports no junction for a `Repeat`, so such a member never
//! reaches `translate_junction_member` — the arm that already carries this
//! hazard for `dup` (#1280, #1292). This is the same class for an edit type that
//! arm does not cover.
//!
//! **Not an adjudication.** There is no competing-legal-forms question here, so
//! `canonical-form-choice-when-both-legal` does not reach it: the output denotes
//! a different sequence than the input, which `README.md`'s ruleset makes a
//! rule-7 **bug** rather than a house-style choice ("best effort is bounded by
//! the spec's determinacy … not by ferro's implementation quality"). Hence an
//! ordinary regression test and no `rulings` record.

use crate::common::cis_apply_oracle::{apply_parsed_reason, apply_reason};
use crate::common::synthetic::{padded, SyntheticBuilder};
use ferro_hgvs::hgvs::edit::NaEdit;
use ferro_hgvs::hgvs::uncertainty::Mu;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

/// The contig the sweep in #1597 drew these rows from. Core base 1 sits at
/// `g.257`, so `g.268` is core index 11 — a `C`, and the reason this core is the
/// one pinned here rather than a homopolymer where every slide lands on an `A`
/// and the defect is invisible.
const CORE: &str = "TCAACGTATAGCAGCACTAC";

/// Every `(start, end, unit)` triple a normalized description states for the
/// repeat members it contains, on the genomic axis it was written on.
fn repeat_members(variant: &HgvsVariant) -> Vec<(u64, u64, String)> {
    let members: Vec<&HgvsVariant> = match variant {
        HgvsVariant::Allele(allele) => allele.variants.iter().collect(),
        single => vec![single],
    };
    members
        .into_iter()
        .filter_map(|member| {
            let HgvsVariant::Genome(g) = member else {
                return None;
            };
            let NaEdit::Repeat {
                sequence: Some(unit),
                ..
            } = g.loc_edit.edit.inner()?
            else {
                return None;
            };
            let (Some(Mu::Certain(start)), Some(Mu::Certain(end))) = (
                g.loc_edit.location.start.as_single(),
                g.loc_edit.location.end.as_single(),
            ) else {
                return None;
            };
            Some((start.base, end.base, unit.to_string()))
        })
        .collect()
}

/// Assert that every repeat member of `input`'s normalized form sits on a tract
/// of its own unit, and that the allele still denotes the input's sequence.
fn assert_repeats_sit_on_their_unit(input: &str, direction: ShuffleDirection) -> String {
    let provider = SyntheticBuilder::genomic(CORE).build();
    let reference = padded(CORE);
    let normalizer = Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(direction),
    );
    let variant = parse_hgvs(input).expect("input parses");
    let want =
        apply_parsed_reason(&provider, &reference, &variant).expect("input denotes a sequence");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    let output = normalized.to_string();

    for (start, end, unit) in repeat_members(&normalized) {
        // `padded` is 1-based genomic, so `g.N` is `reference[N - 1]`.
        let bases = &reference[(start - 1) as usize..end as usize];
        // Cyclic, matching what the fix enforces: the span must read as the unit
        // repeated, allowing a partial trailing copy. Asserting a whole number of
        // copies instead would pin a stricter contract than the code makes, and
        // would go red on a legitimate `g.1_9AT[4]`-shaped member.
        let expected: String = unit.chars().cycle().take(bases.len()).collect();
        assert_eq!(
            bases, expected,
            "{input} ({direction:?}) -> {output}: repeat g.{start}_{end}{unit}[..] claims a \
             tract of {unit}, but the reference reads {bases} there"
        );
    }

    // A repeat on the wrong bases is not merely mis-spelled — it changes what the
    // allele denotes, which is the rule-1 failure #1597 is about. Assert that too,
    // so the guard cannot be satisfied by a differently-wrong output.
    let got = apply_reason(&provider, &reference, &output);
    assert_eq!(
        got.as_deref(),
        Ok(want.as_str()),
        "{input} ({direction:?}) -> {output}: output does not denote the input's sequence"
    );
    output
}

#[test]
fn a_repeat_is_not_clamped_onto_a_base_that_is_not_its_unit() {
    // The row the clamp trace was taken from. Both directions emit the same
    // wrong allele on `main`, so both are pinned.
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        let output = assert_repeats_sit_on_their_unit(
            "NC_TEST.1:g.[267_268insCA;268C>A;269_270insAA]",
            direction,
        );
        assert_eq!(output, "NC_TEST.1:g.269A[5]");
    }
}

#[test]
fn the_clamped_repeat_survives_beside_a_sibling_it_cannot_merge_with() {
    // Same slide, but the third member's payload keeps the two members distinct,
    // so the repeat stays a member rather than being coalesced away. This is the
    // discriminating case: a fix that repaired the row above only by letting the
    // two repeats collapse would leave this one wrong.
    let output = assert_repeats_sit_on_their_unit(
        "NC_TEST.1:g.[267_268insCA;268C>A;269_270insTA]",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(output, "NC_TEST.1:g.[269A[3];269_270insTA]");
}
