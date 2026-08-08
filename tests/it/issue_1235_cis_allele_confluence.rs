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

/// #1232's minimal spanning `delins` is RETAINED. Adjudicated 2026-08-08.
///
/// `g.35_39delinsTA` asserts one changed block: reference `CAATT` at 35–39
/// replaced by `TA`. Splitting it into `[35_37del;39T>A]` requires aligning 5
/// reference nt against a 2-nt payload and then reading the one column that
/// happens to agree (the `T` at 38) as an unchanged run.
///
/// That is the construction `DNA/delins.md:46` performs on its own worked
/// example — "parts of the inserted sequence "align" with the reference
/// sequence, giving an alternative description like
/// `c.[850_869del;874_881del;887_897del;901_902insG]`" — and `DNA/delins.md:47`
/// rejects it one line later: "**The "delins" format is recommended**: it is
/// simpler and prevents software tools making incorrect predictions for the
/// consequences on protein level." The ruling record
/// `delins-merge-vs-individual-gap-two-or-more` (DECIDED 2026-08-07) settles it
/// for `:47` and scopes the ruling to exactly this class: "a MINIMAL single
/// `delins` that would be split because payload bases coincide with reference
/// bases".
///
/// The strongest evidence is in the old doc comment for this very case, which
/// conceded that the tie-break between `[35_37del;39T>A]` and
/// `[35C>T;37_39del]` was "an implementer's choice". Non-uniqueness of the
/// split is ground (4) of the ruling: adopting one trades a unique canonical
/// form for an arbitrary pick out of a family.
///
/// DISCLOSURE: this moves the output for `g.35_39delinsTA` from
/// `[35_37del;39T>A]` back to itself.
#[test]
fn issue_1232_minimal_spanning_delins_is_retained() {
    converges_to("TEMPLATE:g.35_39delinsTA", &["TEMPLATE:g.35_39delinsTA"]);
}

/// NOTE ON THE NAME: it says "splits" and the assertions say nothing splits.
/// The name is pinned by `msto_regression_corpus::MSTO_ISSUES` (#1232, #1235)
/// and `every_cataloged_test_exists_in_its_source_file` fails if it moves, so it
/// is left alone; renaming it and the catalog entry together is a follow-up.
///
/// #1232's THREE SPLIT spellings each keep their own partition. RE-BLESSED
/// 2026-08-08.
///
/// ```text
/// TEMPLATE:g.[35_37del;39T>A]        -> TEMPLATE:g.[35_37del;39T>A]
/// TEMPLATE:g.[35C>T;37_39del]        -> TEMPLATE:g.[35C>T;37_39del]
/// TEMPLATE:g.[35_36delinsT;38_39del] -> TEMPLATE:g.[35_36delinsT;38_39del]
/// ```
///
/// All three previously converged on `[35_37del;39T>A]`, picked by an
/// edit-distance minimiser out of a family the old doc comment for this case
/// already admitted was an "implementer's choice". They now each stay put.
///
/// `general.md:34` governs each independently — "two variants separated by one
/// or more nucleotides should be described individually and **not** as a
/// "delins"" — and each spelling asserts exactly two blocks one unchanged
/// nucleotide apart, so each is already in the form `:34` asks for.
/// `general.md:35`'s codon exception cannot reach a `g.` description. What no
/// clause does is rank one *partition* against another: `general.md:56`'s
/// ladder ("(1) substitution, (2) deletion, (3) inversion, (4) duplication, (5)
/// insertion") ranks edit types within a member, not whole partitions, and
/// `delins` is not even on it.
///
/// The ruling record `canonical-form-choice-when-both-legal` (OPERATOR RULING
/// 2026-08-08, which REVERSES the 2026-08-07 ruling on that record) answers the
/// remaining question directly: "Two descriptors asserting different partitions
/// are different variants under this ruling and legitimately reach different
/// canonical forms." These three assert different partitions — they disagree
/// about which nucleotides changed, not merely how to spell one agreed change.
///
/// DISCLOSURE: this is a representation change and it is the costly direction.
/// Sequence-level confluence across these three spellings is gone. The ruling
/// names that cost and assigns it elsewhere: "Sequence-level equivalence still
/// needs an answer for consumers who dedupe; that is
/// `EquivalenceLevel::SequenceMatch` and a groupable SPDI key, not the
/// canonical string."
///
/// Each is pinned as an exact string AND as a fixed point, and their mutual
/// distinctness is pinned too — so a future change that re-converges them fails
/// on a named string rather than silently restoring the minimiser.
#[test]
fn issue_1232_spanning_delins_splits_at_unchanged_base() {
    let spellings = [
        "TEMPLATE:g.[35_37del;39T>A]",
        "TEMPLATE:g.[35C>T;37_39del]",
        "TEMPLATE:g.[35_36delinsT;38_39del]",
    ];
    for spelling in spellings {
        converges_to(spelling, &[spelling]);
    }
    // Not merely "each is a fixed point" — they must remain three, so this
    // cannot pass by two of them collapsing together.
    for (i, a) in spellings.iter().enumerate() {
        for b in &spellings[i + 1..] {
            assert_ne!(
                normalize_to_string(a),
                normalize_to_string(b),
                "`{a}` and `{b}` assert different partitions and must stay \
                 distinct canonical forms"
            );
        }
    }
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
///
/// RE-BLESSED 2026-08-08. The masking-invariance property is untouched and is
/// still the point of this test; only the shared answer moved, from the split
/// `[35_37del;39T>A]` to the retained `g.35_39delinsTA` — see
/// [`issue_1232_minimal_spanning_delins_is_retained`] for the adjudication
/// (`DNA/delins.md:47`, ruling `delins-merge-vs-individual-gap-two-or-more`).
///
/// The pin is *strengthened* rather than merely updated: both the masked and
/// the upper-case answer are now pinned to the exact string, so neither side
/// can drift while the equality still holds. The failure mode the pin exists
/// for is the reverse of the original one — case comparison could now make the
/// masked window split where the upper-case one does not, which the equality
/// catches and the exact pins name.
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
    // Pin the shared answer on BOTH sides, so this cannot pass by both
    // regressing together.
    assert_eq!(
        upper, "TEMPLATE:g.35_39delinsTA",
        "the retained minimal delins is the expected canonical form \
         (DNA/delins.md:47)"
    );
    assert_eq!(masked, "TEMPLATE:g.35_39delinsTA");
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
///
/// # RE-BLESSED 2026-08-08 — the `dup` typing survives, the merge does not.
///
/// ```text
/// TEMPLATE:g.[74_75insC;75_76insG;90A>T]  ->  TEMPLATE:g.[74_75insC;75dup;90A>T]
/// TEMPLATE:g.[73_74dup;90A>T]             ->  TEMPLATE:g.[73_74dup;90A>T]
/// ```
///
/// The test's own claim — that a derived tandem insertion must be **typed** as a
/// duplication — still holds, and is what the first row shows: the member
/// `g.75_76insG` inserts a `G` directly 3' of the `G` at 75 and is emitted as
/// `g.75dup`, at the 3'-most position of the `GG` tract (`general.md:42`, "the
/// 3'rule also applies for changes in single residue stretches"). So
/// `DNA/duplication.md:18` — "when a variant can be described as a duplication,
/// it **must** be described as a duplication and not as, e.g., an insertion" —
/// is honoured on the member it reaches. What no longer happens is the *merge*
/// of the two asserted insertion members into one `g.73_74dup`.
///
/// WHY THE MERGE IS NOT REQUIRED, since `:18` is one of only two
/// lowercase-**must** clauses in the recommendations and this is the reading
/// that changed. `:18`'s antecedent is "when a variant **can** be described as
/// a duplication", and the very next line tells you that is not a question
/// about sequence identity: `DNA/duplication.md:19` — "when there is **no
/// evidence** that the extra copy of a sequence detected is in tandem (directly
/// 3'-flanking the original copy), the change can not be described as a
/// duplication; it should be described as **an insertion**". Tandemness is
/// evidence, and it is not recoverable from a reference plus a resulting
/// sequence. The evidence available here is the input, which asserts two
/// separate insertions straddling an unchanged `G` — not one tandem copy at
/// 73_74. `general.md:34` then keeps those two blocks individual, since one
/// unchanged nucleotide separates them and `general.md:35`'s codon exception
/// cannot reach a `g.` description.
///
/// The governing record is `canonical-form-choice-when-both-legal` (OPERATOR
/// RULING 2026-08-08, REVERSING the 2026-08-07 ruling on that record), which
/// cites `duplication.md:19` for exactly this point and concludes: "Two
/// descriptors asserting different partitions are different variants under this
/// ruling and legitimately reach different canonical forms."
///
/// DISCLOSURE: `g.[74_75insC;75_76insG;90A>T]` moves from
/// `g.[73_74dup;90A>T]` to `g.[74_75insC;75dup;90A>T]`, and the two spellings no
/// longer converge. Both are pinned exactly, and their distinctness is pinned,
/// so a change that re-merges them names the string it moved.
#[test]
fn a_derived_tandem_insertion_renders_as_a_duplication() {
    // The derived piece is still TYPED as a duplication: `75_76insG` -> `75dup`.
    converges_to(
        "TEMPLATE:g.[74_75insC;75dup;90A>T]",
        &[
            "TEMPLATE:g.[74_75insC;75_76insG;90A>T]",
            "TEMPLATE:g.[74_75insC;75dup;90A>T]",
        ],
    );
    // The one-member spelling of the same resulting sequence is its own
    // canonical form, and stays a `dup` rather than being split.
    converges_to(
        "TEMPLATE:g.[73_74dup;90A>T]",
        &["TEMPLATE:g.[73_74dup;90A>T]"],
    );
    assert_ne!(
        normalize_to_string("TEMPLATE:g.[74_75insC;75_76insG;90A>T]"),
        normalize_to_string("TEMPLATE:g.[73_74dup;90A>T]"),
        "the two-block and one-block spellings assert different partitions"
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
        // The insertion-beside-duplication row, restored. It normalizes to
        // `g.[74_75insC;75dup;90A>T]`, where the insertion's 3' anchor and the
        // duplication's base are both 75. That is flush, not overlapping — see
        // `member_span`, which now returns an empty range for an insertion. The
        // row was excluded while the helper still read the location interval
        // verbatim; excluding it cost the coverage this test exists for, so the
        // helper was fixed instead.
        "TEMPLATE:g.[74_75insC;75_76insG;90A>T]",
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
                // Half-open, so `start == prev` is adjacency and is allowed;
                // `start < prev` is a genuine overlap of claimed bases.
                assert!(
                    start >= prev,
                    "`{input}` -> `{normalized}`: member claiming [{start}, {end}) overlaps the \
                     previous member, which claimed up to {prev}"
                );
            }
            prev_end = Some(end);
        }
    }
}

/// Span of a simple genomic member, as `(start, end)` inclusive.
/// The reference territory a member CLAIMS, as a HALF-OPEN range `[start, end)`.
///
/// Reading the location interval verbatim is wrong for an insertion, and that is
/// not a detail: `A_B ins` names its two flanking bases as **anchors it reads**,
/// not as bases it changes. It occupies the *junction* between them, so the next
/// member may legitimately begin at `B` — which is exactly what
/// `g.[74_75insC;75dup;…]` does. Reading the interval verbatim reported that pair
/// as overlapping and made this test fail on a representation question rather
/// than on the disjointness property it guards.
///
/// So a pure insertion returns the EMPTY range `[B, B)`, positioned at its 3'
/// anchor, and every base-claiming edit returns `[A, B + 1)`. This is the same
/// interbase-versus-HGVS-coordinate distinction #1307 records: an insertion is
/// flush with a neighbour that starts where it ends, and only a genuine
/// base-claiming overlap is a violation.
fn member_span(variant: &ferro_hgvs::hgvs::variant::HgvsVariant) -> Option<(i64, i64)> {
    use ferro_hgvs::hgvs::edit::NaEdit;
    use ferro_hgvs::hgvs::variant::HgvsVariant;
    let HgvsVariant::Genome(g) = variant else {
        return None;
    };
    let interval = &g.loc_edit.location;
    let start = interval.start.inner()?.base as i64;
    let end = interval.end.inner()?.base as i64;
    // A pure insertion consumes no reference base.
    if matches!(g.loc_edit.edit.inner(), Some(NaEdit::Insertion { .. })) {
        return Some((end, end));
    }
    Some((start, end + 1))
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
