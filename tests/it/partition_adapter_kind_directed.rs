//! Task 4: the Direction-2 kind-directed adapter (`anchor_for_member`) honors the
//! `EditKind` an arm chose, and the bake-off arms hard-error through the seam.
//!
//! The whole point of the adapter (the finding this task exists to fix): a member
//! an arm typed `Delins` must render as `delins`, NOT be promoted by the geometry
//! ladder into `inv`/`dup`; a `Dup` member whose bytes are not a duplication must be
//! REFUSED, never demoted to an insertion; and an arm whose partition does not
//! round-trip on designed net-insertion geometry must surface `Unsound` through the
//! seam rather than shipping a wrong string.

use std::fs;
use std::sync::OnceLock;

use serde::Deserialize;

use ferro_hgvs::normalize::{
    normalize_block_via_arm, render_member_via_adapter, MemberRenderError,
};
use ferro_hgvs::partition::output::{EditKind, Member};
use ferro_hgvs::partition::partitioner::PartitionError;

use crate::common::fixture_gen;

/// An L2 control arm that types the whole trimmed block as one `delins`, registered
/// as `some-l1/always-one-delins/off`.
const ALWAYS_DELINS_ARM: &str = "some-l1/always-one-delins/off";

/// A graded R2-winner arm. `op-extract-v2` is 0-invalid on the real corpus but
/// UNSOUND on designed net-insertion geometry — the case exercised below.
const TRIM_ONLY_ARM: &str = "trim-only/op-extract-v2/ledger-r5";

/// A member the arm explicitly typed `Delins` must render as `delins` and must NOT
/// be re-typed to `inv` by the geometry ladder.
///
/// The block `TT ATCG TT` -> `TT CGAT TT` replaces `ATCG` with its reverse
/// complement `CGAT`: the recommended-form ladder (`LadderL2`) would type that whole
/// span as `inv`. The `always-one-delins` control instead declares `Delins`, and the
/// adapter must preserve that choice. `contains("delins")` alone is sufficient
/// (ferro's rendered kinds are keyword-exclusive); the `!inv` check is a redundant
/// belt-and-braces.
#[test]
fn an_arm_that_types_delins_is_not_re_typed_to_inv_or_dup() {
    let reference = b"TTATCGTT";
    let resulting = b"TTCGATTT";
    let out = normalize_block_via_arm(ALWAYS_DELINS_ARM, reference, resulting)
        .expect("the always-one-delins control is sound on this block");
    assert!(
        out.contains("delins"),
        "the arm's Delins choice must survive the adapter, got {out:?}"
    );
    assert!(
        !out.contains("inv"),
        "the geometry ladder must not override the arm, got {out:?}"
    );
}

/// A `Dup` member whose payload does not copy the 5' reference flank must FAIL the
/// render gate (`NotADuplication`), never be silently demoted to an `Insertion`.
#[test]
fn a_dup_member_whose_bytes_do_not_match_the_5prime_flank_is_refused_not_demoted() {
    // Insertion point at offset 3; the 5' flank is reference[1..3] == "GG", which is
    // not the payload "AC", so this is not a duplication.
    let reference = b"GGGGG";
    let not_a_dup = Member {
        kind: EditKind::Dup,
        ref_start: 3,
        ref_end: 3,
        inserted: b"AC".to_vec(),
    };
    let refused = render_member_via_adapter(&not_a_dup, reference);
    assert!(
        matches!(refused, Err(MemberRenderError::NotADuplication)),
        "a Dup whose bytes are not a duplication must be refused, got {refused:?}"
    );

    // Positive control: a genuine tandem duplication (payload == 5' flank) renders
    // as `dup`, so the refusal above is discriminating, not a blanket decline.
    let real_dup = Member {
        kind: EditKind::Dup,
        ref_start: 3,
        ref_end: 3,
        inserted: b"AC".to_vec(),
    };
    let rendered =
        render_member_via_adapter(&real_dup, b"TAC").expect("a genuine tandem duplication renders");
    assert!(
        rendered.contains("dup") && !rendered.contains("ins"),
        "a real duplication renders as dup, got {rendered:?}"
    );
}

#[test]
fn a_sub_whose_reference_base_is_unreadable_is_refused_not_delins() {
    // The Sub arm renders `X>Y` only if it can read the replaced reference base. If
    // it cannot — `ref_start` out of window here — it must REFUSE, never let
    // `build_naedit` silently fall through to a one-base `delins`. That silent kind
    // change is exactly the override this whole adapter exists to prevent, so the Sub
    // arm must refuse it as squarely as the Dup arm refuses a non-duplication
    // (whole-branch review m9). The geometry is a valid Sub (span 1, one payload
    // base); only the reference read fails.
    let reference = b"A"; // one base — offset 3 is past the end
    let unreadable_sub = Member {
        kind: EditKind::Sub,
        ref_start: 3,
        ref_end: 4,
        inserted: b"G".to_vec(),
    };
    let refused = render_member_via_adapter(&unreadable_sub, reference);
    assert!(
        matches!(refused, Err(MemberRenderError::KindGeometryMismatch)),
        "a Sub whose reference base is unreadable must be refused, got {refused:?}"
    );

    // Positive control: the same Sub geometry with a readable reference base renders
    // as a substitution, so the refusal above is discriminating, not a blanket
    // decline of every Sub.
    let readable_sub = Member {
        kind: EditKind::Sub,
        ref_start: 0,
        ref_end: 1,
        inserted: b"G".to_vec(),
    };
    let rendered = render_member_via_adapter(&readable_sub, b"A")
        .expect("a Sub with a readable reference base renders");
    assert!(
        rendered.contains('>') && !rendered.contains("delins"),
        "a readable Sub renders as a substitution, got {rendered:?}"
    );
}

#[test]
fn an_identity_member_renders_as_an_equals_not_a_del() {
    // The `Identity` arm (an examined re-attachment) renders `=`, not `del`. It is
    // dead on the ∅-provenance seam (no examined re-attachment is produced there), so
    // it is covered here directly rather than left as an unverified assumption
    // (whole-branch review m10).
    let identity = Member {
        kind: EditKind::Identity,
        ref_start: 2,
        ref_end: 2,
        inserted: Vec::new(),
    };
    let rendered =
        render_member_via_adapter(&identity, b"ACGT").expect("a zero-width identity renders");
    assert!(
        rendered.contains('=') && !rendered.contains("del"),
        "an Identity member renders as `=`, got {rendered:?}"
    );
}

/// An arm's UNSOUND output on designed net-insertion geometry must hard-error through
/// the seam, not ship garbage.
///
/// GROUNDED: `op-extract-v2` measured `unsound=56` on `cis_confluence_corpus`, and the
/// diagnostic pinned `s00-c-m2-sep1-p8-all-ins` — two adjacent insertions that
/// reconstruct to the wrong bases. Its `(core, denoted)` is the block; the seam must
/// return `PartitionError::Unsound`.
#[test]
fn known_unsound_arms_hard_error_through_the_seam() {
    let class = unsound_case();
    let result = normalize_block_via_arm(
        TRIM_ONLY_ARM,
        class.core.as_bytes(),
        class.denoted.as_bytes(),
    );
    assert!(
        matches!(result, Err(PartitionError::Unsound)),
        "an arm's unsound output on designed geometry must hard-error, got {result:?}"
    );
}

// ---------------------------------------------------------------------------
// Corpus loading — a minimal mirror of the generator's `Class`, so this test can
// read `denoted` (private to `cis_confluence_axis`). The corpus is regenerated on
// demand exactly as `cis_confluence_axis` does, off the same `--axes g,c` default.
// ---------------------------------------------------------------------------

const CORPUS_RELATIVE_PATH: &str = "tests/fixtures/cis/cis_confluence_corpus.json";
const UNSOUND_CASE_ID: &str = "s00-c-m2-sep1-p8-all-ins";

#[derive(Deserialize)]
struct Class {
    id: String,
    core: String,
    denoted: String,
}

#[derive(Deserialize)]
struct Corpus {
    classes: Vec<Class>,
}

fn corpus() -> &'static Corpus {
    static CORPUS: OnceLock<Corpus> = OnceLock::new();
    CORPUS.get_or_init(|| {
        let path = fixture_gen::fixture_path(CORPUS_RELATIVE_PATH);
        fixture_gen::ensure_generated_example_fixture(
            &path,
            "generate_cis_confluence_corpus",
            &["--axes", "g,c"],
            ".cis_confluence_corpus",
            "cis confluence corpus",
            || {},
        );
        let text = fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
        serde_json::from_str(&text)
            .unwrap_or_else(|e| panic!("failed to parse {}: {e}", path.display()))
    })
}

fn unsound_case() -> &'static Class {
    corpus()
        .classes
        .iter()
        .find(|c| c.id == UNSOUND_CASE_ID)
        .unwrap_or_else(|| panic!("corpus is missing the pinned case {UNSOUND_CASE_ID:?}"))
}
