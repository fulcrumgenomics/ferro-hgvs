//! Cis-allele confluence: #1229-#1235.
//!
//! `Normalizer::normalize` used to normalize each cis-allele member in
//! isolation rather than re-deriving the allele from the resulting sequence.
//! Two consequences, both covered here:
//!
//! * **non-confluence** — several encodings of one variant each normalized to a
//!   different stable string (#1229, #1230, #1231, #1232, #1233), so any
//!   consumer keying on the normalized string silently over-split the variant;
//! * **silent corruption** — a deletion 3'-shifted across a sibling
//!   substitution, emitting overlapping members that denote a *different*
//!   sequence (#1234).
//!
//! Only the first is covered here. #1234 needs a *sibling clamp* in the
//! per-member pipeline, which is #1259's work: this pass runs after the damage
//! is done, and `apply_edits_to_window` refuses an overlapping member list by
//! design, since such an allele has no well-defined resulting sequence to
//! re-derive from.
//!
//! The synthetic `TEMPLATE` reference is the one from the issue reproductions,
//! kept verbatim so these tests and the reports stay comparable.

use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use std::io::Write;

const SEQ: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

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
            "exons": [{
                "number": 1, "start": 1, "end": n,
                "genomic_start": 1, "genomic_end": n
            }]
        }],
        "genomic_sequences": { "TEMPLATE": sequence }
    });
    let mut f = tempfile::NamedTempFile::new().expect("tempfile");
    write!(f, "{json}").expect("write json");
    JsonProvider::from_json(f.path()).expect("load reference")
}

fn provider() -> JsonProvider {
    provider_for(SEQ)
}

fn normalize_with(provider: JsonProvider, input: &str) -> String {
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"));
    format!("{normalized}")
}

fn normalize_to_string(input: &str) -> String {
    normalize_with(provider(), input)
}

/// Assert every encoding normalizes to `canonical`, and that `canonical` is
/// itself a fixed point. The second half is what makes this a confluence
/// check rather than a table of expected strings.
fn converges_to(canonical: &str, encodings: &[&str]) {
    for input in encodings {
        assert_eq!(
            normalize_to_string(input),
            canonical,
            "`{input}` should normalize to `{canonical}`"
        );
    }
    assert_eq!(
        normalize_to_string(canonical),
        canonical,
        "`{canonical}` should be a normalization fixed point"
    );
}

// ----------------------------------------------------------------------------
// Per-issue regressions
// ----------------------------------------------------------------------------

/// #1229 — an `inv` member must not block coalescing of an adjacent run.
///
/// Ref `CAC` at 4-6 becomes `AGT`. Three consecutive changed nucleotides are
/// one delins (`delins.md:16`); `AGT` is not the reverse complement of `CAC`
/// (that is `GTG`), so this is not an inversion.
#[test]
fn issue_1229_adjacent_inv_member_coalesces() {
    converges_to(
        "TEMPLATE:g.4_6delinsAGT",
        &[
            "TEMPLATE:g.[4C>A;5A>G;6C>T]",
            "TEMPLATE:g.[4C>A;5_6inv]",
            "TEMPLATE:g.4_6delinsAGT",
        ],
    );
}

/// #1230 — an `inv` whose interior bases are unchanged is really the end
/// substitutions.
///
/// Ref `GATG` at 20-23 becomes `CATC`. Only positions 20 and 23 change; the
/// interior `AT` is a complementary pair and is untouched. Changes separated by
/// unchanged nucleotides are described individually (`inversion.md:21`), and
/// prioritisation ranks substitution above inversion (`general.md:56`).
///
/// The spec has no worked example of an inversion with an unchanged interior,
/// so this direction is a documented implementer's choice, not a spec citation.
#[test]
fn issue_1230_inv_with_unchanged_interior_becomes_substitutions() {
    converges_to(
        "TEMPLATE:g.[20G>C;23G>C]",
        &[
            "TEMPLATE:g.20_23inv",
            "TEMPLATE:g.20_23delinsCATC",
            "TEMPLATE:g.[20G>C;23G>C]",
        ],
    );
}

/// #1231 — `[dup;del]` reduces to the substitutions it actually describes.
///
/// Positions 36 and 38 change, separated by unchanged 37. Prioritisation puts
/// substitution above duplication (`general.md:56`).
#[test]
fn issue_1231_dup_del_reduces_to_substitutions() {
    converges_to(
        "TEMPLATE:g.[36A>C;38T>A]",
        &["TEMPLATE:g.[35dup;38del]", "TEMPLATE:g.[36A>C;38T>A]"],
    );
}

/// #1232 — a spanning delins must not be retained across an unchanged interior
/// base (`delins.md:17`).
///
/// All four encodings below are stable fixed points on v0.11.0 — the issue
/// reported two of them. The tie-break between `[35_37del;39T>A]` and
/// `[35C>T;37_39del]` (equally minimal, equally stable) is an implementer's
/// choice: the partition is taken 5'-most, then each member is 3'-shifted.
#[test]
fn issue_1232_spanning_delins_splits_at_unchanged_base() {
    converges_to(
        "TEMPLATE:g.[35_37del;39T>A]",
        &[
            "TEMPLATE:g.35_39delinsTA",
            "TEMPLATE:g.[35_37del;39T>A]",
            "TEMPLATE:g.[35C>T;37_39del]",
            "TEMPLATE:g.[35_36delinsT;38_39del]",
        ],
    );
}

/// Soft-masking must not change the canonical form.
///
/// Reference FASTAs are routinely soft-masked, so the same variant can be
/// fetched as `acgt…` in a repeat-rich region and `ACGT…` everywhere else. The
/// sequence-first pass mixes provider bytes with bytes it writes itself
/// (`Base::to_u8()` is always upper-case) and then compares them exactly:
/// `best_alignment` scores a matched column on `reference[i] == result[j]`. On a
/// lower-case window every upper-case replacement base therefore mismatched a
/// reference base it was in fact identical to, so the aligner found no interior
/// match and #1232's split silently collapsed back to one spanning `delins` —
/// one variant with two canonical strings, decided by whether its region
/// happened to be masked.
///
/// Pinned against the #1232 shape specifically, because that is the case whose
/// answer depends on an interior base surviving unchanged.
#[test]
fn soft_masked_reference_yields_the_same_canonical_form() {
    let input = "TEMPLATE:g.35_39delinsTA";
    let upper = normalize_with(provider_for(SEQ), input);
    let masked = normalize_with(provider_for(&SEQ.to_lowercase()), input);

    assert_eq!(
        masked, upper,
        "a soft-masked reference must give the same canonical form as an \
         upper-case one for `{input}`"
    );
    // Pin the shared answer too, so this cannot pass by both sides regressing
    // to the same spanning delins.
    assert_eq!(
        upper, "TEMPLATE:g.[35_37del;39T>A]",
        "the #1232 split is the expected canonical form"
    );
}

/// #1233 — `[ins;del]` reduces to substitutions. Highest-frequency shape in the
/// reporter's variant-counting run.
#[test]
fn issue_1233_ins_del_reduces_to_substitutions() {
    converges_to(
        "TEMPLATE:g.[36A>T;38T>A]",
        &["TEMPLATE:g.[35_36insT;38del]", "TEMPLATE:g.[36A>T;38T>A]"],
    );
}

/// A derived piece that is a tandem duplication must be **typed** as one.
///
/// `duplication.md:18` — "when a variant can be described as a duplication, it
/// must be described as a duplication" — is one of the spec's genuine MUSTs, and
/// the sequence-first derivation used to be unable to honour it: it recognised
/// the shape (`is_tandem_duplication`) only in order to **refuse** the whole
/// group, handing the allele back to the per-member pipeline, which types each
/// member in isolation and so never merges the two insertions below.
///
/// `g.[74_75insC;75_76insG]` and `g.73_74dup` denote one variant: reference
/// `…C(73) G(74) G(75) C(76)…` becomes `…C G C G G C…` either way, because the
/// interleaved `C`/`G` of the split spelling and the contiguous `CG` of the dup
/// commute across the `G` at 75. The substitution at 90 is not decoration — it
/// is what keeps the derivation at two pieces, so this exercises the `dup`
/// typing rather than the lone-insertion path.
#[test]
fn a_derived_tandem_insertion_renders_as_a_duplication() {
    converges_to(
        "TEMPLATE:g.[73_74dup;90A>T]",
        &[
            "TEMPLATE:g.[74_75insC;75_76insG;90A>T]",
            "TEMPLATE:g.[73_74dup;90A>T]",
        ],
    );
}

// #1234's clamp-dependent assertions are **not** here.
//
// `[4_6del;8A>T]` -> `[6_8del;8A>T]` — a deletion shifting over its sibling and
// emitting overlapping members — needs a sibling clamp in the per-member
// pipeline. This PR's sequence-first pass cannot repair it: `apply_edits_to_window`
// deliberately refuses an overlapping member list, because an allele whose
// members overlap has no well-defined resulting sequence to re-derive from. The
// corruption has already happened by the time this pass sees it.
//
// That clamp is #1259's, so its regressions belong with it. They were carried
// here while this branch was stacked on the (now closed) #1236 and are dropped
// on the rebase onto `main` rather than left failing or weakened to pass.

/// #1234 — a lone deletion with no sibling still shifts to its standalone
/// 3'-most position. Guards against fixing the clamp by disabling the shift.
#[test]
fn issue_1234_lone_deletion_still_shifts_fully() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.4_6del"),
        "TEMPLATE:g.6_8del"
    );
}

// ----------------------------------------------------------------------------
// #1235 acceptance criteria
// ----------------------------------------------------------------------------

/// #1235 — no normalized cis allele may contain overlapping or out-of-order
/// members. `[4_6del;8A>T]` violated this before the fix.
#[test]
fn normalized_cis_members_are_disjoint_and_ordered() {
    use ferro_hgvs::hgvs::variant::HgvsVariant;

    // `[4_6del;8A>T]` and `[5_7del;8A>T]` are deliberately absent: both need the
    // sibling clamp (#1259) to be disjoint at all, and this pass cannot supply
    // it — see the note above. Every input below is a partition decision this
    // PR owns.
    let inputs = [
        "TEMPLATE:g.[35dup;38del]",
        "TEMPLATE:g.[35_36insT;38del]",
        "TEMPLATE:g.[4C>A;5_6inv]",
        "TEMPLATE:g.[35_37del;39T>A]",
        "TEMPLATE:g.[20G>C;23G>C]",
    ];
    let normalizer = Normalizer::new(provider());
    for input in inputs {
        let parsed = parse_hgvs(input).expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        let HgvsVariant::Allele(allele) = &normalized else {
            continue; // collapsed to a single member: trivially disjoint
        };
        let mut prev_end: Option<i64> = None;
        for member in &allele.variants {
            let (start, end) = member_span(member)
                .unwrap_or_else(|| panic!("no span for member of `{normalized}`"));
            assert!(
                start <= end,
                "`{normalized}`: member span {start}..{end} is inverted"
            );
            if let Some(prev) = prev_end {
                assert!(
                    start > prev,
                    "`{input}` -> `{normalized}`: member starting at {start} overlaps or \
                     precedes the previous member ending at {prev}"
                );
            }
            prev_end = Some(end);
        }
    }
}

/// Span of a simple genomic member, as `(start, end)` inclusive.
fn member_span(variant: &ferro_hgvs::hgvs::variant::HgvsVariant) -> Option<(i64, i64)> {
    use ferro_hgvs::hgvs::variant::HgvsVariant;
    let HgvsVariant::Genome(g) = variant else {
        return None;
    };
    let interval = &g.loc_edit.location;
    let start = interval.start.inner()?.base as i64;
    let end = interval.end.inner()?.base as i64;
    Some((start, end))
}

/// #1235 — `normalize` is idempotent for every case in this file.
#[test]
fn normalization_is_idempotent() {
    let inputs = [
        "TEMPLATE:g.[4C>A;5A>G;6C>T]",
        "TEMPLATE:g.[4C>A;5_6inv]",
        "TEMPLATE:g.20_23inv",
        "TEMPLATE:g.20_23delinsCATC",
        "TEMPLATE:g.[35dup;38del]",
        "TEMPLATE:g.35_39delinsTA",
        "TEMPLATE:g.[35C>T;37_39del]",
        "TEMPLATE:g.[35_36insT;38del]",
        "TEMPLATE:g.[4_6del;8A>T]",
        "TEMPLATE:g.4_8delinsCT",
        "TEMPLATE:g.4_6del",
    ];
    for input in inputs {
        let once = normalize_to_string(input);
        let twice = normalize_to_string(&once);
        assert_eq!(once, twice, "`{input}` is not idempotent");
    }
}

/// #1234 — `normalize` must never change the resulting sequence. Checked
/// through ferro's own `EquivalenceChecker`, which reported the old corrupt
/// output as equal to a *different* variant.
#[test]
fn normalization_preserves_the_resulting_sequence() {
    use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};

    // `[4_6del;8A>T]` and `[5_7del;8A>T]` are deliberately absent: both need the
    // sibling clamp (#1259) to be disjoint at all, and this pass cannot supply
    // it — see the note above. Every input below is a partition decision this
    // PR owns.
    let inputs = [
        "TEMPLATE:g.[35dup;38del]",
        "TEMPLATE:g.[35_36insT;38del]",
        "TEMPLATE:g.[4C>A;5_6inv]",
        "TEMPLATE:g.35_39delinsTA",
        "TEMPLATE:g.20_23inv",
        // Derives a pure-insertion piece that actually shifts, so its payload
        // must rotate. This is the only input here that covers that, and this
        // assertion is the only oracle for it that survives `--release` (the
        // round-trip check inside the pass is `debug_assert`-only).
        "TEMPLATE:g.[4C>A;6_7insCGA]",
    ];
    let checker = EquivalenceChecker::new(provider());
    let normalizer = Normalizer::new(provider());
    for input in inputs {
        let parsed = parse_hgvs(input).expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        let level = checker.check(&parsed, &normalized).expect("check").level;
        // Any rung at or above `SequenceMatch` means the resulting sequence
        // survived; `NormalizedMatch` is the stronger claim that the two also
        // normalize to one string.
        assert!(
            matches!(
                level,
                EquivalenceLevel::Identical
                    | EquivalenceLevel::NormalizedMatch
                    | EquivalenceLevel::SequenceMatch
            ),
            "`{input}` -> `{normalized}` changed the resulting sequence ({level:?})"
        );
    }
}

/// Both public normalization exits must emit the same canonical variant.
///
/// `normalize()` routes through `normalize_core_checked`; `normalize_with_diagnostics()`
/// is the lenient exit the Python bindings and the web service call. Before the
/// canonicalization pass they differed only in error and diagnostic behavior,
/// never in the emitted variant. Wiring the pass into just one of them made them
/// disagree on the canonical string for exactly the shapes #1229-#1235 exist to
/// canonicalize — a Python caller got the old non-confluent answer, and the
/// shift `infos` were computed against a variant the library no longer
/// considered canonical. The idempotency oracle cannot see this: it lives inside
/// `normalize_core_checked`, on the side that was wired.
#[test]
fn both_public_exits_emit_the_same_canonical_variant() {
    let inputs = [
        "TEMPLATE:g.35_39delinsTA",
        "TEMPLATE:g.20_23inv",
        "TEMPLATE:g.[35dup;38del]",
        "TEMPLATE:g.[4C>A;5_6inv]",
        "TEMPLATE:g.[35_36insT;38del]",
    ];
    let normalizer = Normalizer::new(provider());
    for input in inputs {
        let parsed = parse_hgvs(input).expect("parse");
        let plain = normalizer
            .normalize(&parsed)
            .expect("normalize")
            .to_string();
        let with_diagnostics = normalizer
            .normalize_with_diagnostics(&parsed)
            .expect("normalize_with_diagnostics")
            .result
            .to_string();
        assert_eq!(
            with_diagnostics, plain,
            "`{input}`: normalize_with_diagnostics must agree with normalize",
        );
    }
}

/// A derived pure-insertion piece must rotate its payload when it 3'-shifts.
///
/// `shuffle` returns coordinates only. Moving an insertion point 3' by `k`
/// rotates the inserted sequence left by `k mod len` — `shuffle`'s own predicate
/// defines it that way, comparing `ref[start + k]` against
/// `alt[(new_end - start) % alt.len()]`. Writing back the shifted coordinates
/// while leaving `alt` alone therefore emits a description denoting a
/// *different* sequence, and because the wrong form is a stable fixed point no
/// idempotency check can see it.
///
/// Here the derived insertion of `CGA` before `g.7` shifts one position, because
/// `ref[7] == 'C' == payload[0]`; a second step is blocked (`payload[1] == 'G'`
/// against `ref[8] == 'A'`). So the payload rotates to `GAC` at `7_8`. The
/// sibling `4C>A` keeps the piece count above the lone-insertion bail.
///
/// A three-base payload is deliberate: with one base every rotation is the
/// identity, so a shorter payload would pass even unrotated and prove nothing.
#[test]
fn shifted_insertion_piece_rotates_its_payload() {
    assert_eq!(
        normalize_to_string("TEMPLATE:g.[4C>A;6_7insCGA]"),
        "TEMPLATE:g.[4C>A;7_8insGAC]",
        "the shifted insertion must be spelled from the ROTATED payload",
    );
}
