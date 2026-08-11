//! A protein cis allele whose adjacency is *created* by the 3' shift must
//! coalesce in the same pass that created it.
//!
//! Found by #1614's protein corpus axis under `FERRO_ASSERT_IDEMPOTENT=1`, which
//! is the point of adding that axis: the shape is unreachable from any
//! hand-written protein test in the suite, because it needs a reference peptide
//! with a residue run long enough for a member to shift *into* its sibling.
//!
//! The defect. `coalesce_protein_adjacent_changes` runs once, on the **authored**
//! members, before the collapse/merge/shift loop — so it only ever saw adjacency
//! the input already had. The per-member 3' shift then created more of it, and
//! the nucleotide `merge_consecutive_edits` inside the loop never fires on
//! protein members, so nothing closed the gap:
//!
//! ```text
//! input   NP_TEST.1:p.[Gly13del;Gly16Ala]
//! pass 1  NP_TEST.1:p.[Gly16Ala;Gly17del]     <- the deletion shifted 13 -> 17
//! pass 2  NP_TEST.1:p.Gly16_Gly17delinsAla    <- and only now merged
//! ```
//!
//! Two passes disagreeing is exactly what `FERRO_ASSERT_IDEMPOTENT` exists to
//! catch, and pass 2 is the correct answer rather than pass 1:
//! `protein/substitution.md:32` marks the split spelling of two adjacent members
//! `class="invalid"`, and the decided
//! `rulings[delins-adjacent-members-when-both-consume-reference]` record scopes
//! that ruling to members which both consume reference bases — a substitution and
//! a single-residue deletion both do. So the fix is to coalesce inside the loop,
//! not to accept the bracket form.
//!
//! The *shift* half — what makes `Gly13del` -> `Gly17del` obligatory rather than
//! merely something ferro happens to do — is `general.md:157-160`. The argument
//! and its verbatim quote live at the coalesce site in `src/normalize/mod.rs`,
//! beside the code they justify, deliberately not restated here: one rule written
//! in two places and then drifting apart is this project's named recurring
//! failure mode (`rulings[adjudication-precedence-order]`). The control test
//! below exercises that clause directly.
//!
//! **The `Gly13del` control below is what keeps these assertions honest.** The
//! merge can only be observed if the deletion really travels from 13 to 17; were
//! the shift to stop happening, the two members would never become adjacent and a
//! bare "the output is the merged form" assertion would be satisfied by the
//! absence of the shift rather than by the presence of the merge. Assert the
//! movement separately, so a regression in either half is named.

use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

const ACCESSION: &str = "NP_TEST.1";

/// Positions 13..=17 are a five-residue glycine run, so a deletion authored at
/// 13 shifts to 17 and lands flush against a substitution at 16. Residue 18 is
/// not glycine, which is what bounds the run.
///
/// ```text
/// 1  2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
/// M  A C D E F H I K L  M  N  G  G  G  G  G  A  C  D  E  F
/// ```
const PEPTIDE: &str = "MACDEFHIKLMNGGGGGACDEF";

fn normalizer() -> Normalizer<MockProvider> {
    let mut provider = MockProvider::new();
    provider.add_protein(ACCESSION, PEPTIDE);
    Normalizer::new(provider)
}

fn normalize_to_string(input: &str) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{input}`: {e}"));
    let normalized = normalizer()
        .normalize(&variant)
        .unwrap_or_else(|e| panic!("normalize failed for `{input}`: {e}"));
    format!("{normalized}")
}

/// The control: the deletion really does travel the glycine run.
///
/// Without this, the merge assertions below could pass because the shift stopped
/// happening rather than because the coalesce started.
#[test]
fn the_deletion_shifts_to_the_end_of_the_glycine_run() {
    assert_eq!(
        normalize_to_string(&format!("{ACCESSION}:p.Gly13del")),
        format!("{ACCESSION}:p.Gly17del"),
        "the 3' rule must carry a deletion to the last position of its run — if it \
         does not, `p.[Gly13del;Gly16Ala]` never becomes adjacent and the merge \
         assertions in this file are vacuous"
    );
}

/// The defect itself: one pass, not two.
#[test]
fn shift_created_adjacency_coalesces_in_the_pass_that_created_it() {
    let once = normalize_to_string(&format!("{ACCESSION}:p.[Gly13del;Gly16Ala]"));
    assert_eq!(
        once,
        format!("{ACCESSION}:p.Gly16_Gly17delinsAla"),
        "the shifted deletion is flush against the substitution, so \
         `protein/substitution.md:32` makes the merged delins the correct form; \
         emitting the bracket spelling here is the #1614 defect"
    );
    let twice = normalize_to_string(&once);
    assert_eq!(
        once, twice,
        "and the merged form must be its own fixed point, or the pass merely moved \
         the non-idempotency one step later"
    );
}

/// The bracket form the defect used to emit must itself normalize to the merged
/// form, so a consumer holding the old string converges rather than sitting on a
/// second fixed point.
#[test]
fn the_previously_emitted_bracket_form_converges_on_the_merged_form() {
    assert_eq!(
        normalize_to_string(&format!("{ACCESSION}:p.[Gly16Ala;Gly17del]")),
        format!("{ACCESSION}:p.Gly16_Gly17delinsAla"),
        "the string the defect shipped must converge, not persist as a rival form"
    );
}

/// Authored adjacency still reaches the merged form — by whichever pass reaches
/// it first.
///
/// This deliberately does **not** attribute the merge to the pre-loop coalesce.
/// The in-loop coalesce this change adds runs on every cis pass and would merge
/// authored adjacency on its own, so no assertion here can tell the two routes
/// apart; what is pinned is the outcome. (An earlier version of this comment
/// claimed the test guarded the pre-loop coalesce specifically. It never could.)
///
/// **Pinned against the merged literal, not against the sibling spelling's
/// output.** Comparing the two orders to each other proves only
/// order-independence, and order-independence is satisfied by the coalesce *not
/// firing at all*: cis members are sorted by `cis_member_order_key` before
/// rendering, so both spellings converge on one bracket string either way and the
/// assertion could not fail for the reason its name gives. Pinning the literal
/// here, with the forward spelling pinned in
/// `the_previously_emitted_bracket_form_converges_on_the_merged_form`, keeps
/// #1116's order-independence — two literals that agree — while making each
/// direction falsifiable on its own.
#[test]
fn authored_adjacency_still_coalesces() {
    assert_eq!(
        normalize_to_string(&format!("{ACCESSION}:p.[Gly17del;Gly16Ala]")),
        format!("{ACCESSION}:p.Gly16_Gly17delinsAla"),
        "authored adjacency must reach the merged form from either member order \
         (#1116) — and reaching it is what this asserts, not merely agreeing with \
         the other order, which the coordinate-order render guarantees anyway"
    );
}
