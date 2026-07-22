//! A protein cis-allele run whose to-`Ter` (nonsense) member is **last**
//! coalesces to a trailing-`Ter` `delins` — #1129, narrowing the conservative
//! decline #1095 introduced and #1125 scoped to the run.
//!
//! # Spec
//!
//! `protein/substitution.md:23` / `protein/delins.md:18` require two or more
//! consecutive changed amino acids to be described as a `delins`, and pin the
//! bracket form as wrong:
//!
//! > the description `p.Arg76_Cys77delinsSerTrp` is correct, the description
//! > `p.[Arg76Ser;Cys77Trp]` is not correct.
//!
//! A trailing `Ter` in the inserted peptide is explicitly endorsed —
//! `protein/delins.md:47`:
//!
//! > **`p.(Pro578_Lys579delinsLeuTer)`** — a deletion of amino acids `Pro578`
//! > and `Lys579` replaced by `LeuTer` (alternatively `Leu*`).
//!
//! and `protein/delins.md:43` gives the single-residue form
//! (`NP_004371.2:p.(Asn47delinsSerSerTer)`).
//!
//! Only two `Ter` shapes stay forbidden, and both still decline here:
//!
//! - residues listed **after** the stop — `protein/delins.md:45`: "the
//!   deletion-insertion is not described as `delinsSerSerTerAlaAsp`, amino acids
//!   after the translation termination codon are **not** listed". So a run whose
//!   `Ter` member is interior must not merge.
//! - an **immediate** stop — `protein/substitution.md:20`: "not
//!   `p.Cys5_Ser6delinsTerGluAsp` but `p.Tyr4Ter`". So a run whose *first*
//!   member is the `Ter` must not merge either.
//!
//! These rewrites are reference-independent (all residues are named in the
//! AST), so an empty `MockProvider` exercises them.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Pin the canonical form of `input` and require it to be a normalization fixed
/// point. Idempotency is load-bearing here: the merged `delins…Ter` must not be
/// re-trimmed or re-expanded on a second pass.
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

/// Every permutation of `items` via Heap's algorithm.
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
/// canonicalizes to `expected` — the `Ter`'s position in the *run* decides, not
/// its position in the author's input (#1116's order-independence contract).
fn all_orders_canonicalize_to(members: &[&str], expected: &str) {
    for perm in permutations(members) {
        let input = format!("NP_003997.1:p.[{}]", perm.join(";"));
        canonicalizes_to(&input, expected);
    }
}

/// The #1129 repro: a two-member run whose second member is the nonsense
/// substitution merges to the spec's `delins…Ter` form.
#[test]
fn a_run_ending_in_ter_coalesces_to_a_trailing_ter_delins() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Ter"],
        "NP_003997.1:p.Cys100_Asp101delinsSerTer",
    );
}

/// A longer run ending in `Ter` merges the same way — every preceding residue
/// contributes to the insert, which terminates at the stop.
#[test]
fn a_longer_run_ending_in_ter_coalesces() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Gly", "Arg102Ter"],
        "NP_003997.1:p.Cys100_Arg102delinsSerGlyTer",
    );
}

/// A deletion member inside a run ending in `Ter` contributes nothing to the
/// insert, exactly as in #1120's mixed sub/del runs.
#[test]
fn a_del_member_inside_a_run_ending_in_ter_contributes_nothing() {
    all_orders_canonicalize_to(
        &["Cys100del", "Asp101Ter"],
        "NP_003997.1:p.Cys100_Asp101delinsTer",
    );
}

/// Predicted `( )` members stay predicted through the merge.
#[test]
fn a_predicted_run_ending_in_ter_coalesces_predicted() {
    all_orders_canonicalize_to(
        &["(Cys100Ser)", "(Asp101Ter)"],
        "NP_003997.1:p.(Cys100_Asp101delinsSerTer)",
    );
}

/// An **immediate** stop must stay a nonsense substitution: a run whose FIRST
/// member is the `Ter` would merge to `delinsTer…`, which
/// `protein/substitution.md:20` forbids.
#[test]
fn a_run_starting_with_ter_still_declines() {
    all_orders_canonicalize_to(
        &["Arg76Ter", "Cys77Trp"],
        "NP_003997.1:p.[Arg76Ter;Cys77Trp]",
    );
}

/// An **interior** `Ter` would list residues after the stop
/// (`protein/delins.md:45`), so the run still declines.
#[test]
fn a_run_with_an_interior_ter_still_declines() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Ter", "Arg102Trp"],
        "NP_003997.1:p.[Cys100Ser;Asp101Ter;Arg102Trp]",
    );
}

/// #1125's scoping is preserved: a run 5′ of an unrelated downstream `Ter`
/// still coalesces, and that lone `Ter` re-emits verbatim.
#[test]
fn an_unrelated_downstream_ter_still_does_not_block_an_upstream_run() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Gly", "Arg200Ter"],
        "NP_003997.1:p.[Cys100_Asp101delinsSerGly;Arg200Ter]",
    );
}

/// A run 3′ of an earlier stop still declines even though it ends in its own
/// `Ter` — the residues it describes are already downstream of a translation
/// termination, which stays out of scope for this narrow rule (#1125).
#[test]
fn a_run_ending_in_ter_but_downstream_of_an_earlier_ter_declines() {
    all_orders_canonicalize_to(
        &["Arg50Ter", "Cys100Ser", "Asp101Ter"],
        "NP_003997.1:p.[Arg50Ter;Cys100Ser;Asp101Ter]",
    );
}
