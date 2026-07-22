//! Broaden protein cis-allele coalescing to mixed edit-types: a consecutive
//! run mixing single-residue substitutions with single-residue deletions
//! collapses into one `delins` (or a pure range `del` when every member is a
//! deletion) — #1120, extending #1095 (substitution-only).
//!
//! # Spec
//!
//! `protein/delins.md:18` / `protein/substitution.md:23`:
//!
//! > changes involving **two or more consecutive amino acids** are described
//! > as a deletion/insertion variant (delins).
//!
//! "changes" is general — a substitution and a deletion are both changes.
//! `delins.md:41` gives the mixed example directly: `p.Cys28_Lys29delinsTrp`
//! is "a deletion of amino acids `Cys28` and `Lys29`, replaced by `Trp`" — a
//! del+sub run collapsed to one delins.
//!
//! These rewrites are reference-independent (all residues are named in the
//! AST), so an empty `MockProvider` exercises them.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Pin the canonical form ferro's normalizer produces for `input`, and require
/// that form to be a normalization fixed point (idempotent). The idempotency
/// leg is load-bearing: because the coalesce now sorts members before merging,
/// re-normalizing a canonical form must not shift it again (#1120/#1116).
fn canonicalizes_to(input: &str, expected: &str) {
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input:?}: {e}"));
    let normalized = Normalizer::new(MockProvider::new())
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input:?}: {e}"));
    assert_eq!(normalized.to_string(), expected, "for {input:?}");

    let reparsed =
        parse_hgvs(expected).unwrap_or_else(|e| panic!("parse expected {expected:?}: {e}"));
    let renormalized = Normalizer::new(MockProvider::new())
        .normalize(&reparsed)
        .unwrap_or_else(|e| panic!("re-normalize {expected:?}: {e}"));
    assert_eq!(
        renormalized.to_string(),
        expected,
        "canonical form {expected:?} is not a normalization fixed point",
    );
}

/// Every permutation of `members` via Heap's algorithm.
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

/// Assert EVERY permutation of `members` (joined into `NP_003997.1:p.[…]`)
/// canonicalizes to `expected` — and, via `canonicalizes_to`, that the result
/// is an idempotent fixed point. This pins member-order independence
/// exhaustively rather than for a single hand-picked reversed order (#1116).
fn all_orders_canonicalize_to(members: &[&str], expected: &str) {
    for perm in permutations(members) {
        let input = format!("NP_003997.1:p.[{}]", perm.join(";"));
        canonicalizes_to(&input, expected);
    }
}

/// sub + del over adjacent residues → one delins; the deleted residue
/// contributes nothing to the inserted sequence.
#[test]
fn sub_then_del_coalesces_to_delins() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76Ser;Cys77del]",
        "NP_003997.1:p.Arg76_Cys77delinsSer",
    );
}

/// del + sub over adjacent residues → one delins.
#[test]
fn del_then_sub_coalesces_to_delins() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76del;Cys77Trp]",
        "NP_003997.1:p.Arg76_Cys77delinsTrp",
    );
}

/// A run whose members are ALL deletions has an empty inserted sequence, so it
/// collapses to a pure range deletion, not an (invalid) empty `delins`.
#[test]
fn adjacent_dels_coalesce_to_range_deletion() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76del;Cys77del]",
        "NP_003997.1:p.Arg76_Cys77del",
    );
}

/// Three-member mixed run (sub, del, sub) → one delins spanning the whole run.
#[test]
fn three_member_mixed_run_coalesces() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76Ser;Cys77del;Gly78Asp]",
        "NP_003997.1:p.Arg76_Gly78delinsSerAsp",
    );
}

/// Predicted `( )` members with a del sibling stay predicted through the merge.
#[test]
fn predicted_mixed_run_coalesces_to_predicted_delins() {
    canonicalizes_to(
        "NP_003997.1:p.[(Arg76Ser);(Cys77del)]",
        "NP_003997.1:p.(Arg76_Cys77delinsSer)",
    );
}

/// Non-adjacent members are never coalesced (a gap of unchanged residues keeps
/// them separate), even when one is a deletion.
#[test]
fn non_adjacent_sub_and_del_are_not_coalesced() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76Ser;Cys90del]",
        "NP_003997.1:p.[Arg76Ser;Cys90del]",
    );
}

/// A to-`Ter` (nonsense) substitution still bails the whole run even with a
/// del sibling — a `delins…Ter…` run would list residues after the stop
/// (substitution.md:20, delins.md:45).
#[test]
fn to_ter_run_with_del_sibling_is_not_coalesced() {
    canonicalizes_to(
        "NP_003997.1:p.[Arg76Ter;Cys77del]",
        "NP_003997.1:p.[Arg76Ter;Cys77del]",
    );
}

/// Member **input order** must not affect a mixed sub+del coalesce, and the
/// result must be idempotent. Before #1116's sort-before-coalesce reached this
/// mixed path, a descending-order allele (e.g. `del` authored before its
/// preceding `sub`) bailed the ascending-order gate and only the display sort
/// reordered it — so `normalize` was NOT idempotent: a second pass would then
/// coalesce it. Sorting inside the coalesce makes every permutation land on the
/// same delins in one pass.
#[test]
fn mixed_sub_del_coalesces_identically_regardless_of_member_order() {
    // sub+del (both orders) → the same delins; the deleted residue contributes
    // nothing to the insert.
    all_orders_canonicalize_to(
        &["Arg76Ser", "Cys77del"],
        "NP_003997.1:p.Arg76_Cys77delinsSer",
    );
    // A pure-del run collapses to a range deletion, order-independent.
    all_orders_canonicalize_to(&["Arg76del", "Cys77del"], "NP_003997.1:p.Arg76_Cys77del");
    // All 6 permutations of a three-member mixed run land on one delins.
    all_orders_canonicalize_to(
        &["Arg76Ser", "Cys77del", "Gly78Asp"],
        "NP_003997.1:p.Arg76_Gly78delinsSerAsp",
    );
    // Partial coalesce whose surviving lone member is a DELETION: the run
    // [Arg76,Cys77] collapses while the separated `Gly90del` must re-emit
    // verbatim as a `del` — indexed through `source`, not the post-sort
    // position — identically from EVERY input permutation. This is the exact
    // #1116-sort × #1120-mixed-edit seam the hand-merge introduced.
    all_orders_canonicalize_to(
        &["Arg76Ser", "Cys77del", "Gly90del"],
        "NP_003997.1:p.[Arg76_Cys77delinsSer;Gly90del]",
    );
}
