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
//! The synthetic `TEMPLATE` reference is the one from the issue reproductions,
//! kept verbatim so these tests and the reports stay comparable.

use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use std::io::Write;

const SEQ: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

fn provider() -> JsonProvider {
    let n = SEQ.len() as u64;
    let json = serde_json::json!({
        "version": "1.0",
        "transcripts": [{
            "id": "TEMPLATE-gene.1",
            "chromosome": "TEMPLATE",
            "strand": "+",
            "sequence": SEQ,
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
        "genomic_sequences": { "TEMPLATE": SEQ }
    });
    let mut f = tempfile::NamedTempFile::new().expect("tempfile");
    write!(f, "{json}").expect("write json");
    JsonProvider::from_json(f.path()).expect("load reference")
}

fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(provider());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"));
    format!("{normalized}")
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

/// #1233 — `[ins;del]` reduces to substitutions. Highest-frequency shape in the
/// reporter's saturation-mutagenesis counting run.
#[test]
fn issue_1233_ins_del_reduces_to_substitutions() {
    converges_to(
        "TEMPLATE:g.[36A>T;38T>A]",
        &["TEMPLATE:g.[35_36insT;38del]", "TEMPLATE:g.[36A>T;38T>A]"],
    );
}

/// #1234 — a deletion's 3'-shift must stop at a sibling, not cross it.
///
/// `[4_6del;8A>T]` used to emit `[6_8del;8A>T]`: both members claim position 8,
/// which is malformed *and* denotes a different sequence (ferro read it as a
/// bare `6_8del`, swallowing the substitution). The deletion may shift only to
/// `5_7`, after which it is adjacent to the substitution and the two coalesce.
#[test]
fn issue_1234_deletion_shift_clamps_at_sibling() {
    converges_to(
        "TEMPLATE:g.5_8delinsT",
        &[
            "TEMPLATE:g.[4_6del;8A>T]",
            "TEMPLATE:g.[5_7del;8A>T]",
            "TEMPLATE:g.4_8delinsCT",
            "TEMPLATE:g.[5A>T;6_8del]",
            "TEMPLATE:g.[4_5del;7_8delinsT]",
            "TEMPLATE:g.5_8delinsT",
        ],
    );
}

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

    let inputs = [
        "TEMPLATE:g.[4_6del;8A>T]",
        "TEMPLATE:g.[5_7del;8A>T]",
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

    let inputs = [
        "TEMPLATE:g.[4_6del;8A>T]",
        "TEMPLATE:g.[5_7del;8A>T]",
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
