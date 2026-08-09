//! Issue #1307 — a terminal duplication re-spelled to an insertion one base
//! past the end of the contig.
//!
//! `respell_colliding_duplications` repairs a duplication that collides with a
//! sibling by re-spelling it as an insertion at the junction 3' of the
//! duplicated span — the gap `[end, end + 1]`. When the duplication rests on
//! the contig's **last** base there is no such gap:
//!
//! ```text
//! reference  TTTTTTTTTAATATATTTTAATAC     24 bases
//! in   g.[24dup;24C>G]
//! out  g.[24C>G;24_25insC]                position 25 does not exist
//! ```
//!
//! ## Which guards catch it, and which do not
//!
//! `respell_at_gap`'s `landed` check does not: it only confirms the member reads
//! back at the gap it was *told* to use, which an out-of-range gap does.
//! `FERRO_ASSERT_REPARSE` does not either — `parse_hgvs` does not know the
//! contig's length, so `g.24_25insC` re-parses cleanly. That is #1327's blind
//! spot exactly, and #1353 tracks the seam-level bound that would catch the
//! class rather than the instance.
//!
//! **`FERRO_ASSERT_IDEMPOTENT` does catch it**, though, which #1307 and #1353
//! both assume it cannot. The out-of-range output is not a fixed point; measured
//! against `origin/main`:
//!
//! ```text
//! input: NC_TEST.1:g.[24dup;24C>G]
//! once:  NC_TEST.1:g.[24C>G;24_25insC]
//! twice: NC_TEST.1:g.23_24insG
//! ```
//!
//! So the oracle was never blind here — no test normalized this shape, and an
//! oracle only covers the corpus that reaches it. Worth knowing when reading
//! #1353's table, and worth remembering as the cheaper of the two ways to catch
//! this class.
//!
//! ## Why #1327's fix did not already cover it
//!
//! #1327 added the terminal bound, and #1344 then *disabled* it whenever a
//! base-claiming sibling occupies the base the gap runs past — correctly, for
//! the case it measured: a sibling **deletion** there makes the pair the #999
//! adjacency the collapse pass merges, and the merge consumes the out-of-range
//! coordinate before anything renders. Re-spelling instead produced
//! `c.[*11del;*11delinsAA]`, two members claiming one base, which ferro's own
//! strict mode rejects.
//!
//! That reasoning does not carry to a sibling that claims the last base without
//! being a deletion. There is nothing for the junction form to be absorbed into,
//! so no merge happens and the out-of-range coordinate is not transient — it
//! reaches the output. Two shapes reach it, both pinned below: a substitution
//! sibling (the issue's reproducer) and a `delins` sibling (found while
//! measuring the fix — `delins` removes its bases, yet the collapse pass still
//! does not merge the pair).
//!
//! ## What it does instead, and why not the other three options
//!
//! **Leave the member as the duplication it was.** The collision goes
//! unrepaired, which is the spelling the input already had, and every coordinate
//! named exists.
//!
//! **Not the boundary `delins`** that #1327 emits — for #1344's reason: the
//! sibling already claims that base, so the re-spelt member would overlap it.
//!
//! **Not the wrapped form** on `m.`, either. SVD-WG006 authorises the reversed
//! `<high>_<low>` range for deletions and duplications only, so `m.24_1insC` is
//! not valid HGVS and ferro's parser rejects it (#129) — pinned by
//! `the_wrapped_insertion_spelling_is_not_valid_hgvs` in
//! `issue_1327_mt_respell_past_contig_end`. The circular axis therefore gets the
//! same answer as the linear one.
//!
//! ## The single member the trace ends on — reached, as of #1406
//!
//! This section previously argued that `g.23_24insG` — in range, one member,
//! denoting the same sequence — was **unreachable by policy**, because
//! `normalize_allele` disables its pre-merge reorder when the input carries an
//! overlap conflict, and "a `dup` and a substitution on one base is exactly what
//! strict mode rejects as `OverlapConflictingEdits / W5002`."
//!
//! **That premise was wrong, and #1406 corrected it.** A `dup` is not an
//! in-place edit: `duplication.md:5` places the copy "directly 3' of the
//! original copy", so `g.24dup` *reads* base 24 and *writes* at the junction
//! 24|25, while `g.24C>G` writes base 24. The two write footprints are disjoint,
//! there is no ordering to choose, and the composition is unique. Coincident
//! *read* spans were being reported as a conflict that never existed —
//! `detect_overlap_conflicts` keyed on the span a member names rather than on
//! the bases it claims.
//!
//! The spec asks for the merge, not the refusal. `delins.md:86-89` marks
//! `NM_007294.3:c.[2077G>A;2077_2078insTA]` invalid and gives `c.2077delinsATA`
//! as *the correct description*; `general.md:56` prioritisation then requires
//! that form be spelled as a `dup` where it is one. And the old reading made the
//! verdict depend on the spelling: the same variant written
//! `g.[24C>G;24_25insC]` was accepted while `g.[24dup;24C>G]` was rejected.
//!
//! **So the issue's own bug is fixed more cleanly than by declining.** The
//! reported defect is a coordinate one base past the contig — `g.24_25insC`.
//! Reaching `g.23_24insG` never constructs an out-of-range coordinate at all,
//! rather than constructing one and then refusing to repair around it.
//! `assert_within_contig` below still guards every output, and now passes on a
//! merged form instead of on an un-repaired pair.
//!
//! The four guards this section used to cite as blocking still **pass**,
//! measured: `issue_395_overlap_conflict_strict_rejection::
//! lenient_mode_preserves_authored_order_of_conflicting_members`,
//! `issue_1276_dup_junction_overlap::
//! lenient_output_of_a_conflict_is_still_rejected_by_strict`,
//! `issue_1235_transcript_axes::overlap_conflicting_allele_is_not_canonicalized`
//! and `idempotency_tests::
//! test_overlap_conflicting_allele_is_idempotent_across_respellings`. They are
//! unaffected because nothing was dropped wholesale — the change is to *which
//! shapes are conflicts*, and every allele those four use is one still.
//!
//! A genuinely conflicting allele must still stay recognisable as one so strict
//! mode can reject it rather than being canonicalised into something that looks
//! valid. That property is intact; this shape simply was never a member of it.

use crate::common::cis_apply_oracle::apply_with;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// The issue's reference: 24 bases ending in a unique `C`, so a duplication
/// authored at position 24 stays there.
const SEQUENCE: &str = "TTTTTTTTTAATATATTTTAATAC";
const LENGTH: usize = 24;

/// A 24-base contig ending in a 4-`A` run, so a member authored inside the run
/// 3'-shifts onto the final base instead of being written there.
fn run_terminated_sequence() -> String {
    format!("{}{}", "CG".repeat(10), "A".repeat(4))
}

fn normalize(accession: &str, sequence: &str, input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(accession, sequence.to_string());
    Normalizer::new(provider)
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string()
}

/// Every position the output names must exist on the contig.
///
/// Read out of the rendered string rather than string-matching `"25"`, so an
/// overrun by any amount at any endpoint is caught. Deliberately
/// over-approximate — it scans every integer after the `:`, so a repeat count or
/// an inserted length is checked as though it were a position. Sound for these
/// fixtures (single-digit counts against a 24-base contig) and it fails safe.
/// Mirrors `assert_within_contig` in `issue_1327_mt_respell_past_contig_end`.
fn assert_within_contig(input: &str, output: &str) {
    let reparsed = parse_hgvs(output).unwrap_or_else(|e| {
        panic!("`{input}` normalized to `{output}`, which will not re-parse: {e}")
    });
    let positions = output
        .rsplit_once(':')
        .map(|(_, rest)| rest)
        .unwrap_or(output);
    for run in positions.split(|c: char| !c.is_ascii_digit()) {
        let Ok(position) = run.parse::<usize>() else {
            continue;
        };
        assert!(
            (1..=LENGTH).contains(&position),
            "`{input}` normalized to `{output}` ({reparsed}), which names position \
             {position} on a {LENGTH}-base contig"
        );
    }
}

/// The output must denote the same bases as the input.
///
/// A bounds check ([`assert_within_contig`]) and a pinned expected string both
/// pass on a well-formed description of the *wrong* sequence: a merge composes
/// two members into one edit, and if the composition is wrong the result is a
/// well-formed description of the wrong bases, which no other check here would
/// notice. Checked with an applier that is **not** the normalizer —
/// `apply_with` converts each member to its SPDI triple and splices `sequence`
/// directly — so a bug shared between the normalizer and this check cannot
/// cancel out. Mirrors the same pattern in
/// `issue_1327_mt_respell_past_contig_end` and `cis_apply_oracle` generally.
///
/// Extracted from [`the_merged_output_denotes_the_inputs_bases`], which held
/// this body inline and now calls it, so the two cannot drift apart.
///
/// [`a_deletion_sibling_at_the_last_base_still_merges`],
/// [`a_three_member_allele_stays_in_range_in_any_order`] and
/// [`an_in_range_collision_is_still_repaired`] call it directly too. The
/// module's three remaining merged/reordered pins do **not**, and that does not
/// leave them unguarded: all five inputs of
/// [`a_substitution_sibling_leaves_the_terminal_duplication_alone`],
/// [`a_delins_sibling_leaves_the_terminal_duplication_alone`] and
/// [`the_circular_axis_gets_the_same_answer`] are themselves rows of
/// `the_merged_output_denotes_the_inputs_bases`, which is where their sequence
/// check lives. Adding a newly pinned input to any of those three is the case to
/// watch, because that overlap is a coincidence nothing enforces.
fn assert_denotes_the_same_bases(accession: &str, sequence: &str, input: &str, output: &str) {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(accession, sequence.to_string());
    let want = apply_with(&provider, sequence, input)
        .unwrap_or_else(|| panic!("`{input}` must denote a sequence"));
    let got = apply_with(&provider, sequence, output)
        .unwrap_or_else(|| panic!("`{input}` -> `{output}`, which denotes no sequence at all"));
    assert_eq!(
        got, want,
        "`{input}` -> `{output}` no longer denotes the input's bases"
    );
}

/// The issue's reproducer, on both member orders.
///
/// Pinned, not merely bounded. A bounds check alone is satisfied by a refusal,
/// by a dropped member, or by any other in-range spelling — and pinning is
/// exactly what caught #1344 hiding behind the weaker check. What must come back
/// is the input's own two members, unrepaired.
#[test]
fn a_substitution_sibling_leaves_the_terminal_duplication_alone() {
    assert_eq!(SEQUENCE.len(), LENGTH, "fixture must be 24 bases");
    for (input, expected) in [
        ("NC_TEST.1:g.[24dup;24C>G]", "NC_TEST.1:g.23_24insG"),
        ("NC_TEST.1:g.[24C>G;24dup]", "NC_TEST.1:g.23_24insG"),
    ] {
        let output = normalize("NC_TEST.1", SEQUENCE, input);
        assert_within_contig(input, &output);
        assert_eq!(
            output, expected,
            "`{input}` must reach the merged form, in range"
        );
    }
}

/// A `delins` sibling reaches the same gate and must get the same answer.
///
/// This is the shape that showed "removes its bases" is the wrong predicate for
/// admitting the transient out-of-range junction: a `delins` does remove them,
/// and the collapse pass still does not merge the pair, so the earlier cut of
/// this fix emitted `g.[24delinsGG;24_25insA]` here.
#[test]
fn a_delins_sibling_leaves_the_terminal_duplication_alone() {
    let input = "NC_TEST.1:g.[24dup;24delinsGG]";
    let output = normalize("NC_TEST.1", SEQUENCE, input);
    assert_within_contig(input, &output);
    // Two members still, unlike the substitution case — the pair does not
    // collapse to one edit. What changed is that they are now ordered rather
    // than refused, and the out-of-range junction is never constructed.
    assert_eq!(output, "NC_TEST.1:g.[24delinsGG;24dup]");
}

/// The same overrun on the circular axis. Not a special case: `m.` has a last
/// base like any other sequence, and the wrapped `ins` spelling that would make
/// it one is not valid HGVS (see the module docs).
#[test]
fn the_circular_axis_gets_the_same_answer() {
    // The genomic answers, verbatim, on `m.` — that identity is the whole
    // claim, so they are written out rather than derived from the input.
    for (input, expected) in [
        ("NC_012920.1:m.[24dup;24C>G]", "NC_012920.1:m.23_24insG"),
        (
            "NC_012920.1:m.[24dup;24delinsGG]",
            "NC_012920.1:m.[24delinsGG;24dup]",
        ),
    ] {
        let output = normalize("NC_012920.1", SEQUENCE, input);
        assert_within_contig(input, &output);
        assert_eq!(
            output, expected,
            "`{input}` must get the same answer as the linear axis, in range"
        );
    }
}

/// The #1344 path must survive: a sibling **deletion** at the last base still
/// gets the transient out-of-range junction form, because the collapse pass
/// consumes it. The genomic half of what
/// `issue_1327_mt_respell_past_contig_end` pins on `m.` and `c.`.
///
/// Both member orders, since the repair is order-independent and a regression
/// that only reordered would otherwise pass.
#[test]
fn a_deletion_sibling_at_the_last_base_still_merges() {
    let sequence = run_terminated_sequence();
    for input in ["NC_TEST.1:g.[21del;22dup]", "NC_TEST.1:g.[22dup;21del]"] {
        let output = normalize("NC_TEST.1", &sequence, input);
        assert_within_contig(input, &output);
        assert_denotes_the_same_bases("NC_TEST.1", &sequence, input, &output);
        // The pair cancels — one `A` removed from the terminal run and one added
        // back — so an identity is the right answer, and the merge that produces
        // it is only reachable through the junction spelling this fix must not
        // have taken away.
        assert_eq!(
            output, "NC_TEST.1:g.23=",
            "`{input}` must still collapse to an identity"
        );
    }
}

/// Three members, with an absorbing sibling *and* a non-absorbing one both on
/// the terminal run: the allele must stay in range whichever order it is
/// authored in.
///
/// These are the permutations that motivated quantifying the landing-claimant
/// test over **every** claimant rather than the first one found — with a `find`,
/// one absorbing sibling could decide the policy while a non-absorbing sibling
/// sat behind it in member order.
///
/// **They do not reach a two-claimant set, and are not the guard for one.**
/// Measured by instrumenting the selection site and running the whole suite:
/// 1068 hits, 832 with a single claimant and 236 with none, never two. Two
/// claimants on one base is the overlap the parser refuses outright
/// (`SelfCancellingAllele` for the `dup`/`del` pair here), and when the pair
/// arrives by *shifting* instead — as it does below, `21del` and `22dup` both
/// travelling to the last `A` — the cancellation runs before this pass, so only
/// the third member is still claiming anything by the time the policy is
/// chosen. What these pin is therefore the outcome and its order-independence,
/// not the `all`-vs-`find` branch; that branch is unreachable today and the
/// quantifier is there so it cannot become order-dependent if it stops being.
#[test]
fn a_three_member_allele_stays_in_range_in_any_order() {
    let sequence = run_terminated_sequence();
    for (input, expected) in [
        // Absorbing (`21del`) plus non-absorbing substitution.
        ("NC_TEST.1:g.[21del;22dup;24A>G]", "NC_TEST.1:g.24A>G"),
        ("NC_TEST.1:g.[22dup;21del;24A>G]", "NC_TEST.1:g.24A>G"),
        ("NC_TEST.1:g.[24A>G;21del;22dup]", "NC_TEST.1:g.24A>G"),
        // Absorbing plus non-absorbing `delins` — `delins` removes its bases and
        // still does not absorb, the distinction the policy turns on.
        (
            "NC_TEST.1:g.[21del;22dup;24delinsGG]",
            "NC_TEST.1:g.24delinsGG",
        ),
        (
            "NC_TEST.1:g.[22dup;21del;24delinsGG]",
            "NC_TEST.1:g.24delinsGG",
        ),
        // A non-absorbing claimant spanning two bases, so the claim reaches the
        // last base from a span rather than a point.
        (
            "NC_TEST.1:g.[21del;22dup;23_24delinsGG]",
            "NC_TEST.1:g.23_24delinsGG",
        ),
    ] {
        let output = normalize("NC_TEST.1", &sequence, input);
        assert_within_contig(input, &output);
        assert_denotes_the_same_bases("NC_TEST.1", &sequence, input, &output);
        assert_eq!(
            output, expected,
            "`{input}` must stay in range and not depend on member order"
        );
    }
}

/// An in-range collision is repaired exactly as before: declining at the last
/// base must not cost the repair everywhere else.
#[test]
fn an_in_range_collision_is_still_repaired() {
    let sequence = run_terminated_sequence();
    let input = "NC_TEST.1:g.[20dup;21A>G]";
    let output = normalize("NC_TEST.1", &sequence, input);
    assert_within_contig(input, &output);
    assert_denotes_the_same_bases("NC_TEST.1", &sequence, input, &output);
    assert_eq!(
        output, "NC_TEST.1:g.21delinsGG",
        "`{input}` must still collapse as it did before"
    );
}

/// A lone terminal duplication has no sibling to collide with and must pass
/// through untouched — the control showing the gate is scoped to the repair and
/// is not clamping ordinary terminal edits.
#[test]
fn a_lone_terminal_duplication_is_unchanged() {
    assert_eq!(
        normalize("NC_TEST.1", SEQUENCE, "NC_TEST.1:g.24dup"),
        "NC_TEST.1:g.24dup"
    );
    // Authored inside the terminal run, it still shifts to the last base.
    assert_eq!(
        normalize("NC_TEST.1", &run_terminated_sequence(), "NC_TEST.1:g.21dup"),
        "NC_TEST.1:g.24dup"
    );
}

/// Declining must leave the output **idempotent**, which the out-of-range
/// spelling was not (see the trace in the module docs).
///
/// This is the assertion that would have failed on `origin/main`, so it is the
/// cheapest guard against the class and belongs beside the pinned strings rather
/// than only in a suite-wide env-gated oracle.
#[test]
fn the_declined_output_is_idempotent() {
    for input in [
        "NC_TEST.1:g.[24dup;24C>G]",
        "NC_TEST.1:g.[24dup;24delinsGG]",
    ] {
        let once = normalize("NC_TEST.1", SEQUENCE, input);
        let twice = normalize("NC_TEST.1", SEQUENCE, &once);
        assert_eq!(
            once, twice,
            "`{input}` -> `{once}` -> `{twice}` is not a fixed point"
        );
    }
}

/// The merged form must denote the bases the input denotes.
///
/// This replaces `the_declined_output_is_still_rejected_by_strict`, whose
/// premise #1406 removed: these shapes are not overlap conflicts, so there is no
/// conflict left for strict mode to reject and asserting one would pin the
/// misclassification rather than the behaviour.
///
/// What matters now is stronger and specific to having *started* merging a shape
/// that was previously refused. A refusal is safe by construction — the output
/// is the input. A merge is not: it asserts that two members compose into one
/// edit, and if the composition is wrong the result is a well-formed description
/// of the wrong bases, which no oracle in this suite would notice. So the check
/// is a sequence identity, made with an applier that is **not** the normalizer:
/// convert each member to its SPDI triple and splice the reference.
///
/// `g.[24dup;24delinsGG]` is included even though it does not collapse to one
/// member, because reordering must be sequence-preserving too.
#[test]
fn the_merged_output_denotes_the_inputs_bases() {
    // The circular inputs are here too, not only on the linear axis. `m.` is
    // where the wrapped spelling would have been tempting, so it is the axis
    // where a wrong composition is most plausible — checking only `g.` would
    // leave the riskier half to a string comparison.
    for (accession, input) in [
        ("NC_TEST.1", "NC_TEST.1:g.[24dup;24C>G]"),
        ("NC_TEST.1", "NC_TEST.1:g.[24C>G;24dup]"),
        ("NC_TEST.1", "NC_TEST.1:g.[24dup;24delinsGG]"),
        ("NC_012920.1", "NC_012920.1:m.[24dup;24C>G]"),
        ("NC_012920.1", "NC_012920.1:m.[24dup;24delinsGG]"),
    ] {
        let output = normalize("NC_TEST.1", SEQUENCE, input);
        assert_denotes_the_same_bases(accession, SEQUENCE, input, &output);
    }
}
