//! Issue #1103: cis-allele member order must not leak into the merge step.
//!
//! #1101 (closes #1098) renders single-accession cis members in genomic order,
//! restoring order-independence for **disjoint** members by sorting *after*
//! `merge_consecutive_edits`. But `merge_consecutive_edits` runs *before* that
//! sort and merges only pairs that happen to be adjacent *in the input list*, so
//! which merges fire depends on input order — and the post-merge sort cannot
//! undo a merge that already collapsed two members into one (or re-fire a merge
//! that a stray sibling blocked). Adjacency-mergeable member sets therefore
//! still normalized order-dependently: several permutations of the *same* set
//! produced different strings, with different member **counts**.
//!
//! The fix sorts cis members into genomic order *before* the merge, so the
//! merge always sees a canonical member order and fires the same merges
//! regardless of input order.
//!
//! These cases run under `MockProvider::new()` (no reference bases), so no
//! 3'-shift or reference-window collapse is possible — the order-dependence
//! being pinned is purely `merge_consecutive_edits`' list-adjacency behavior.
//! The trio is the mutalyzer-corpus `NG_012337.1` example from the issue.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Normalize `input` against a base-free `MockProvider` and return the rendered
/// string. Panics on parse or normalize failure.
fn normalize_to_string(input: &str) -> String {
    let normalizer = Normalizer::new(MockProvider::new());
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let normalized = normalizer
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"));
    format!("{normalized}")
}

/// Generate every permutation of `items` via Heap's algorithm.
fn permutations<T: Clone>(items: &[T]) -> Vec<Vec<T>> {
    let mut result = Vec::new();
    let mut arr = items.to_vec();
    let n = arr.len();
    let mut c = vec![0usize; n];
    result.push(arr.clone());
    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 {
                arr.swap(0, i);
            } else {
                arr.swap(c[i], i);
            }
            result.push(arr.clone());
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
    result
}

/// The mergeable analogue of #1101's `all_permutations_normalize_to_same_canonical_string`
/// (which uses disjoint members): every permutation of the adjacency-mergeable
/// trio `{104_105insA, 105del, 105_106insC}` must normalize to one canonical
/// string with a stable member count.
///
/// Sorted into genomic order the members are `104_105insA`, `105del`,
/// `105_106insC`; the merge chains all three (insertion flush against the
/// deletion, deletion flush against the trailing insertion) into a single
/// `105delinsAC`. The two flanking insertions land on distinct junctions
/// (one 5' of the deleted base, one 3' of it), so their genomic order — not
/// their input order — fixes the `AC` ordering, making the result unambiguous.
#[test]
fn all_permutations_of_mergeable_trio_normalize_to_one_canonical_string() {
    let members = ["104_105insA", "105del", "105_106insC"];
    let canonical = "NG_012337.1:g.105delinsAC";

    for perm in permutations(&members) {
        let input = format!("NG_012337.1:g.[{}]", perm.join(";"));
        let output = normalize_to_string(&input);
        assert_eq!(
            output, canonical,
            "permutation `{input}` must normalize to the canonical `{canonical}` \
             (order must not leak into the merge)"
        );
        // Idempotent: re-normalizing the canonical form is a no-op.
        assert_eq!(
            normalize_to_string(canonical),
            canonical,
            "canonical form must be a normalization fixed point"
        );
    }
}
