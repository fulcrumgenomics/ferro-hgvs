//! Issue #1116: protein cis-allele member order must not leak into the
//! consecutive-residue delins coalescing.
//!
//! #1103/#1106 make the **nucleotide** merge input-order-independent by sorting
//! cis members into a canonical merge order before `merge_consecutive_edits`.
//! The **protein** axis is canonicalized by a separate pass,
//! `coalesce_protein_adjacent_substitutions`, which runs *before* that sort and
//! declined outright on non-ascending members. A descending-order protein allele
//! therefore only got the #1098/#1101 post-normalize display sort applied and
//! came out as the bracket form `protein/substitution.md:23` explicitly calls
//! "not correct":
//!
//! > […] the description `p.Arg76_Cys77delinsSerTrp` is correct, the
//! > description `p.[Arg76Ser;Cys77Trp]` is not correct.
//!
//! The fix sorts the members into ascending residue order inside the coalescing
//! pass, so every permutation of one adjacent-residue substitution set yields
//! the same canonical delins. Duplicate residue positions still decline — a
//! `n,n` pair is a contradiction, not a delins.
//!
//! These cases run under `MockProvider::new()` (no reference sequence), so the
//! only behavior being pinned is the coalescing pass' member-order handling.

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

/// Assert every permutation of `members` normalizes to `canonical`.
fn all_permutations_normalize_to(accession: &str, members: &[&str], canonical: &str) {
    for perm in permutations(members) {
        let input = format!("{accession}:p.[{}]", perm.join(";"));
        assert_eq!(
            normalize_to_string(&input),
            canonical,
            "permutation `{input}` must normalize to the canonical `{canonical}` \
             (member order must not leak into the protein coalescing pass)"
        );
    }
    // The canonical form is a normalization fixed point (idempotent).
    assert_eq!(
        normalize_to_string(canonical),
        canonical,
        "canonical `{canonical}` must be a normalization fixed point"
    );
}

/// The protein analogue of #1103's mergeable-trio test: every permutation of an
/// adjacent-residue substitution pair must coalesce to the one delins
/// `protein/substitution.md:23` requires. `p.[Arg76Ser;Cys77Trp]` and
/// `p.Arg76_Cys77delinsSerTrp` are the spec's own correct/not-correct pair.
#[test]
fn both_orders_of_an_adjacent_pair_coalesce_to_one_delins() {
    all_permutations_normalize_to(
        "NP_003997.1",
        &["Arg76Ser", "Cys77Trp"],
        "NP_003997.1:p.Arg76_Cys77delinsSerTrp",
    );
}

/// A three-residue run coalesces identically from all six input orders, and the
/// inserted payload follows the residues' ascending order — not the input's.
#[test]
fn all_permutations_of_an_adjacent_run_coalesce_to_one_delins() {
    all_permutations_normalize_to(
        "NP_003997.1",
        &["Arg76Ser", "Cys77Trp", "Gly78Ala"],
        "NP_003997.1:p.Arg76_Gly78delinsSerTrpAla",
    );
}

/// A partial allele — one adjacent run plus a separated member — coalesces only
/// the run, and the separated member is re-emitted verbatim in residue order
/// from every input permutation. This pins that a member which ends up in no run
/// is looked up through its *source* index, not its post-sort position.
#[test]
fn a_partial_run_coalesces_identically_from_every_order() {
    all_permutations_normalize_to(
        "NP_003997.1",
        &["Arg76Ser", "Cys77Trp", "Gly90Asp"],
        "NP_003997.1:p.[Arg76_Cys77delinsSerTrp;Gly90Asp]",
    );
}

/// The leading member of the allele is not privileged: a run that starts at the
/// *last*-authored member coalesces the same way. (The coalesced delins takes
/// its accession from the allele's first member; every member is required to
/// share it, so the stamp is order-independent.)
#[test]
fn a_trailing_authored_run_coalesces_identically_from_every_order() {
    all_permutations_normalize_to(
        "NP_003997.1",
        &["Gly90Asp", "Arg76Ser", "Cys77Trp"],
        "NP_003997.1:p.[Arg76_Cys77delinsSerTrp;Gly90Asp]",
    );
}

/// A gapped run is not adjacent, so it never coalesces (`protein/delins.md:63`
/// pins `p.[Ser44Arg;Trp46Arg]` as *not* describable as a single delins). Both
/// orders must still land on the same canonical bracket — the #1098/#1101
/// display sort's job, unchanged by this fix.
#[test]
fn a_gapped_run_stays_a_bracket_in_both_orders() {
    all_permutations_normalize_to(
        "NP_003997.1",
        &["Ser44Arg", "Trp46Arg"],
        "NP_003997.1:p.[Ser44Arg;Trp46Arg]",
    );
}

/// Two substitutions at the *same* residue are a contradiction, not a delins:
/// the coalescing pass must still decline regardless of input order. Sorting is
/// stable, so the duplicate pair reaches the decline check adjacent, exactly as
/// before.
#[test]
fn duplicate_residue_positions_never_coalesce() {
    // Both input orders decline coalescing and land on the SAME canonical
    // bracket: both edits are retained at residue 76 and the #1098/#1101 display
    // sort orders them deterministically (Cys before Ser), so neither the
    // coalescing decline nor the rendering leaks input order.
    for input in [
        "NP_003997.1:p.[Arg76Ser;Arg76Cys]",
        "NP_003997.1:p.[Arg76Cys;Arg76Ser]",
    ] {
        let output = normalize_to_string(input);
        assert_eq!(
            output, "NP_003997.1:p.[Arg76Cys;Arg76Ser]",
            "`{input}` must keep both same-residue edits as a canonical bracket, \
             not coalesce them into a delins"
        );
        assert!(
            !output.contains("delins"),
            "`{input}` must not coalesce two edits at one residue; got `{output}`"
        );
    }
}

/// A run containing a to-`Ter` member never coalesces (the spec forbids
/// listing residues after the stop, `protein/substitution.md:20` /
/// `protein/delins.md:45`), in either input order. The decline is scoped to
/// the run, not the allele — see `issue_1125_ter_scoped_coalesce` for a run
/// that is 5′ of an unrelated `Ter` and does coalesce.
#[test]
fn a_to_ter_run_never_coalesces_in_either_order() {
    all_permutations_normalize_to(
        "NP_003997.1",
        &["Arg76Ter", "Cys77Trp"],
        "NP_003997.1:p.[Arg76Ter;Cys77Trp]",
    );
}

/// An adjacent run that mixes a predicted `( )` member with a certain one is
/// ambiguous and must decline coalescing (merging would fabricate a
/// certainty the members do not jointly assert). With run membership now keyed
/// off sorted residue position, that decline must hold in either input order —
/// this is the one decline branch the permutation suite otherwise leaves
/// untouched.
#[test]
fn a_mixed_predicted_and_certain_run_never_coalesces_in_either_order() {
    all_permutations_normalize_to(
        "NP_003997.1",
        &["(Arg76Ser)", "Cys77Trp"],
        "NP_003997.1:p.[(Arg76Ser);Cys77Trp]",
    );
}
