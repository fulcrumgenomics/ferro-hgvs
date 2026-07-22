//! Scope the protein cis-allele coalesce's to-`Ter` decline to the affected
//! run(s) instead of the whole allele — #1125.
//!
//! # Spec
//!
//! `protein/substitution.md:23` / `protein/delins.md:18`:
//!
//! > changes involving **two or more consecutive amino acids** are described
//! > as a deletion/insertion variant (delins).
//!
//! The countervailing rules are about a translation stop: residues may not be
//! listed *after* one (`protein/delins.md:45` — not `delinsSerSerTerAlaAsp`),
//! and an *immediate* stop is a nonsense substitution, not a delins
//! (`protein/substitution.md:20` — not `p.Cys5_Ser6delinsTerGluAsp`). A run
//! carrying a to-`Ter` member is declined conservatively on that basis.
//!
//! Either way the constraint is about the **run**, not the allele: a run lying
//! entirely 5′ of the earliest `Ter` is unaffected by the stop and still owes
//! the spec its `delins`. That is what this suite pins.
//!
//! Known gap (pre-existing, tracked separately): a run whose to-`Ter` member is
//! **last** would merge to a *trailing*-`Ter` delins, a form the spec endorses
//! (`protein/delins.md:47`: `p.(Pro578_Lys579delinsLeuTer)`). It is still
//! declined here; the tests below pin the current conservative behavior.
//!
//! Every rewrite here is reference-independent (all residues are named in the
//! AST), so an empty `MockProvider` exercises them.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// Pin the canonical form ferro's normalizer produces for `input`, and require
/// that form to be a normalization fixed point (idempotent). The idempotency
/// leg matters here because a partially-coalesced allele is re-parsed as a
/// bracket whose surviving `Ter` member must not then block a second pass.
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
/// canonicalizes to `expected`, so the `Ter` scoping is independent of where
/// the nonsense member was authored (#1116's order-independence contract).
fn all_orders_canonicalize_to(members: &[&str], expected: &str) {
    for perm in permutations(members) {
        let input = format!("NP_003997.1:p.[{}]", perm.join(";"));
        canonicalizes_to(&input, expected);
    }
}

/// Control: with no `Ter` anywhere, the adjacent run coalesces (the behavior
/// #1125's upstream run is being denied).
#[test]
fn adjacent_run_without_ter_coalesces() {
    canonicalizes_to(
        "NP_003997.1:p.[Cys100Ser;Asp101Gly]",
        "NP_003997.1:p.Cys100_Asp101delinsSerGly",
    );
}

/// The #1125 repro: a lone downstream nonsense member must not block an
/// adjacent run that lies entirely 5′ of it. The run coalesces; the `Ter`
/// re-emits verbatim alongside it.
#[test]
fn downstream_ter_does_not_block_an_upstream_run() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Gly", "Arg200Ter"],
        "NP_003997.1:p.[Cys100_Asp101delinsSerGly;Arg200Ter]",
    );
}

/// A run that *contains* the nonsense member still declines — the conservative
/// pre-existing behavior this PR does not change (see the module header's known
/// gap note for the trailing-`Ter` form the spec would allow).
#[test]
fn a_run_containing_the_ter_still_declines() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Ter"],
        "NP_003997.1:p.[Cys100Ser;Asp101Ter]",
    );
}

/// A run whose members sit at/after the earliest `Ter` also declines — the
/// residues it would describe are downstream of a stop, so coalescing them
/// into a new `delins` is out of scope for this narrow rule.
#[test]
fn a_run_downstream_of_the_ter_declines() {
    all_orders_canonicalize_to(
        &["Arg100Ter", "Cys200Ser", "Asp201Gly"],
        "NP_003997.1:p.[Arg100Ter;Cys200Ser;Asp201Gly]",
    );
}

/// With two nonsense members, the EARLIEST one governs: a run between them is
/// downstream of a stop and declines, while a run 5′ of both coalesces.
#[test]
fn the_earliest_ter_governs_when_several_are_present() {
    all_orders_canonicalize_to(
        &[
            "Cys50Ser",
            "Asp51Gly",
            "Arg100Ter",
            "Gly150Ala",
            "Ala151Val",
            "Trp200Ter",
        ],
        "NP_003997.1:p.[Cys50_Asp51delinsSerGly;Arg100Ter;Gly150Ala;Ala151Val;Trp200Ter]",
    );
}

/// A `Ter` strictly adjacent to a run is *part of* that run (adjacency is what
/// defines a run), so the run contains the stop and takes the same conservative
/// decline — the scoping does not let an adjacent stop be merged in.
#[test]
fn a_ter_adjacent_to_a_run_joins_it_and_declines() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Gly", "Arg102Ter"],
        "NP_003997.1:p.[Cys100Ser;Asp101Gly;Arg102Ter]",
    );
}

/// Tightest coalescible boundary: one unchanged residue separates the run from
/// the `Ter`, so the run is a separate run AND ends strictly before the stop.
#[test]
fn a_run_ending_just_before_the_ter_coalesces() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Gly", "Arg103Ter"],
        "NP_003997.1:p.[Cys100_Asp101delinsSerGly;Arg103Ter]",
    );
}

/// A deletion sibling downstream of a coalescible run is unaffected by the
/// scoping: the run merges, the `Ter` and the separated `del` re-emit verbatim.
#[test]
fn scoping_preserves_unrelated_lone_members() {
    all_orders_canonicalize_to(
        &["Cys100Ser", "Asp101Gly", "Arg200Ter", "Gly300del"],
        "NP_003997.1:p.[Cys100_Asp101delinsSerGly;Arg200Ter;Gly300del]",
    );
}
