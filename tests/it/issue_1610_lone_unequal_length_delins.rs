//! #1610 — a lone, minimal, **unequal-length** `delins` must not be split into
//! members by an alignment coincidence.
//!
//! # The defect
//!
//! ```text
//! NM_TEST.1:n.2_5delinsAAC   ->   NM_TEST.1:n.[2_3delinsAA;9del]
//! ```
//!
//! Four reference bases (`CGCG`) against a three-base payload (`AAC`). The only
//! interior "unchanged" column is the payload's own `C` landing on the reference
//! `C` at `n.4`; ferro cut there and 3'-shifted the residual deletion out along
//! the `G` run to `n.9`. The output asserts two variants seven bases apart where
//! the input asserted one contiguous change.
//!
//! It reproduced on **every** partition arm — `live`, `canonical` and the shipped
//! `canonical-coalesced` — so it is not a consequence of the #1835 default flip.
//!
//! # Why the spanning form is the one to emit
//!
//! `DNA/delins.md:46` constructs this description itself — "parts of the
//! inserted sequence 'align' with the reference sequence, giving an alternative
//! description like `c.[850_869del;874_881del;887_897del;901_902insG]`" — and
//! `:47` answers it: "**The 'delins' format is recommended**: it is simpler and
//! prevents software tools making incorrect predictions for the consequences on
//! protein level." So the spec names the decomposition as the *alternative* and
//! recommends the span; ferro emitted the alternative.
//!
//! Both forms are conformant, so this is a `README.md` **rule 6** choice — the
//! maintainers choose among conformant forms — disclosed under rule 7. It is not
//! a conformance finding and must not be cited as one. `README.md`'s own ruleset
//! supplies the ground:
//!
//! > The variant's decomposition is not recoverable. Recovering one means
//! > *choosing* an alignment, and the spec does not say which, so there is no
//! > derivable form to converge on.
//!
//! An earlier revision of this module quoted a second sentence from that bullet,
//! contrasting equal-length blocks with unequal-length ones, as the ground the
//! rule keys on. **#1937 deleted that sentence from `README.md` as false**, on
//! the authority of `rulings[unchanged-is-read-over-every-minimal-alignment]`:
//! `CAG -> AGA` is equal length with edit distance 2, not the position-wise 3.
//! Do not reintroduce it in any form. The length gate's ground is the decided
//! direction scope of `delins-merge-vs-individual-gap-two-or-more`.
//!
//! # THE SHIPPED RULE IS SCOPED TO `c.`, SO THE ISSUE'S OWN REPRODUCTION IS NOT CLOSED
//!
//! Read this before citing #1610 as fixed. The rule is gated on
//! `CoincidenceCarveOut::may_disbelieve_a_separation`, so it runs on the coding
//! DNA axis and nowhere else, per
//! `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`.
//!
//! The reproduction above is spelled `n.`, and `n.` is outside that scope. So
//! `NM_TEST.1:n.2_5delinsAAC` **still** returns `n.[2_3delinsAA;9del]`, and
//! `a_lone_unequal_length_delins_is_not_split` pins that rather than the fix.
//! What this change does close is the same block on `c.`
//! (`the_rule_is_scoped_to_the_coding_dna_axis`).
//!
//! The scope is not a hedge. Ungated, this rule sat directly after the sibling
//! collapse in `partition_block` — which carries the same gate — and re-admitted
//! on the frameless axes exactly what that gate excludes, taking
//! `spec_conformance_axis`'s `guard_violations` from **0 to 5**, all of them the
//! **rejected** SVD-WG010 shape. That is a rank-1 conformance regression, which
//! rule 2 outranks a rule 6 preference. Closing the `n.` and `g.` halves means
//! answering that first, and it is not answered here.
//!
//! # Where the rule lives
//!
//! `merge::split_is_a_placed_gap_coincidence`, consulted by `partition_block`
//! (the `live` arm's own escape ladder) and by `partition_block_for_rule` (the
//! sequence-first arms). Its four conditions and the negative controls each one
//! answers are on that function; the three controls below are the ones that
//! would go wrong first.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// One single-exon transcript carrying `core`, coding or not.
///
/// The coding arm sets the CDS to the whole transcript, so `c.N` and `n.N` name
/// the same transcript position and the `c.`/`n.` comparisons below differ only
/// in whether a reading frame exists.
fn provider(accession: &str, core: &str, coding: bool) -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = core.to_string();
    let length = sequence.len() as u64;
    let (cds_start, cds_end) = if coding {
        (Some(1), Some(length))
    } else {
        (None, None)
    };
    provider.add_transcript(Transcript::new(
        accession.to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        sequence,
        cds_start,
        cds_end,
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalize(provider: &MockProvider, descriptor: &str) -> String {
    let variant = parse_hgvs(descriptor).expect("fixture must parse");
    Normalizer::new(provider.clone())
        .normalize(&variant)
        .expect("fixture must normalize")
        .to_string()
}

/// Normalize, and assert the result is its own normal form.
///
/// Asserted here rather than left to `FERRO_ASSERT_IDEMPOTENT`, so these tests
/// fail in a plain `cargo nextest run` too — CI's gating `Test` job runs without
/// the oracles.
fn normalized_fixed_point(provider: &MockProvider, descriptor: &str) -> String {
    let once = normalize(provider, descriptor);
    let twice = normalize(provider, &once);
    assert_eq!(
        twice, once,
        "`{descriptor}` normalized to a form that is not a fixed point",
    );
    once
}

/// `n.2_5` is `CGCG`; `n.5`-`n.9` is a five-base `G` run, which is what let the
/// residual deletion shift out to `n.9`.
const CORE: &str = "ACGCGGGGGTTTTTTTTTTTTTTTTTTTT";

/// The issue's reproduction, on the axis it was filed against — **still split**.
///
/// This test pins the defect, not the fix. `n.` is outside the shipped rule's
/// `c.` scope (see the module docs), so the exact string #1610 reports is
/// unchanged by this change. It is kept, and kept under the issue's own name, so
/// that the open half of the issue has a committed guard rather than living only
/// in a PR description — and so that closing it later is a visible flip here.
///
/// The counterpart that this change *does* close is
/// `the_rule_is_scoped_to_the_coding_dna_axis` below, on the same block.
#[test]
fn a_lone_unequal_length_delins_is_not_split() {
    let coding = provider("NM_TEST.1", CORE, true);
    assert_eq!(
        normalized_fixed_point(&coding, "NM_TEST.1:n.2_5delinsAAC"),
        "NM_TEST.1:n.[2_3delinsAA;9del]",
        "the `n.` axis is outside `delins.md:47`'s carve-out, so `general.md:34` \
         governs; #1610's own reproduction is NOT closed by this change",
    );
}

/// The same block on `c.`, on `n.`, and on a transcript with no CDS at all.
///
/// This is the axis scope pinned end to end. All three use the same core, the
/// same span and the same payload, so the axis is the only thing that varies —
/// which is what makes the `c.`/`n.` disagreement below evidence about the
/// carve-out rather than about the block.
///
/// The scope is `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`:
/// `delins.md:47` reaches `c.` and nothing else, its stated reason (preventing
/// "incorrect predictions for the consequences on protein level") having nothing
/// to bite on where no protein is coded.
#[test]
fn the_rule_is_scoped_to_the_coding_dna_axis() {
    // In reach: the coding DNA axis, where `delins.md:47` governs.
    let coding = provider("NM_TEST.1", CORE, true);
    assert_eq!(
        normalized_fixed_point(&coding, "NM_TEST.1:c.2_5delinsAAC"),
        "NM_TEST.1:c.2_5delinsAAC",
        "on `c.` the block is kept whole (`DNA/delins.md:47`)",
    );
    // Out of reach: `n.` on the very same coding transcript, so the reading
    // frame's presence is not what differs — only the axis the caller spelled.
    assert_eq!(
        normalized_fixed_point(&coding, "NM_TEST.1:n.2_5delinsAAC"),
        "NM_TEST.1:n.[2_3delinsAA;9del]",
        "`n.` on a coding transcript is still outside the carve-out",
    );
    // Out of reach: a transcript with no CDS at all.
    let noncoding = provider("NR_TEST.1", CORE, false);
    assert_eq!(
        normalized_fixed_point(&noncoding, "NR_TEST.1:n.2_5delinsAAC"),
        "NR_TEST.1:n.[2_3delinsAA;9del]",
        "a CDS-less transcript is outside the carve-out",
    );
}

/// The other side of the boundary: an **equal**-length block still splits.
///
/// `c.2_5` is `CGCG` and the payload `AGGT` is four bases. Ferro derives `c.2`
/// changed, `c.3` unchanged, `c.4_5` changed consecutively, and
/// `DNA/delins.md:17` requires those members individually. Nothing about #1610
/// may reach this: the two cases read identically as strings and only the net
/// length tells them apart.
///
/// Do **not** justify this by saying an equal-length block's columns pair
/// uniquely — `rulings[unchanged-is-read-over-every-minimal-alignment]` holds
/// otherwise and #1937 deleted that claim from `README.md`. This block simply is
/// not a net deletion, which is condition 1.
///
/// # WHAT THIS DOES AND DOES NOT GUARD — measured by sabotage, not assumed
///
/// Relaxing the length gate alone (`result.len() >= reference.len()` to `>`, so
/// equal-length blocks are admitted) leaves this test **green**. Its split is
/// held apart by a different condition: `c.2C>A` is a lone substitution, which
/// is a rank-1 type the split buys, so the rule declines on condition 4 whatever
/// the lengths are. That is not a defect in this test — an equal-length block a
/// full `normalize` emits as a split is always held by a substitution or a
/// compensating insertion (condition 4 or the pure-insertion refusal), never by
/// the length gate. Measured closing #1914: the `>=`->`>` sabotage moves zero
/// rows over ~1.6M swept spelling inputs, so no spelling this module can write
/// isolates the gate end to end.
///
/// **Two guards do pin the gate itself, and neither is reachable by a spelling:**
/// `merge::tests::the_placed_gap_predicate_discriminates` asks the predicate
/// directly with hand-built pieces and an artificial equal-length `result`, and
/// `merge::tests::the_length_gate_holds_a_partitioner_derived_equal_length_block`
/// (#1914) drives the **canonical partitioner** on `AACAAC -> CGAGGA`, whose
/// derived `[delins;delins;placed-gap]` split makes the length gate the only
/// condition that declines the collapse. Both go red under the `>=`->`>`
/// sabotage. Keep all three: this one pins the *behaviour* at the boundary, and
/// those two pin the *rule* — the second on a block the partitioner actually
/// produces.
///
/// # ASSERTED ON `c.`, DELIBERATELY
///
/// This and the three controls below are spelled `c.` because the rule is gated
/// to the coding DNA axis. On `n.` they would still pass — and would pass
/// **vacuously**, declined by the axis before their own condition was ever
/// reached, which is a control that cannot fail for the reason it names.
#[test]
fn an_equal_length_block_keeps_its_split() {
    let coding = provider("NM_TEST.1", CORE, true);
    assert_eq!(
        normalized_fixed_point(&coding, "NM_TEST.1:c.2_5delinsAGGT"),
        "NM_TEST.1:c.[2C>A;4_5delinsGT]",
    );
}

/// A split whose members are **all pure deletions** stays split.
///
/// `:46`'s mechanism is an *inserted sequence* re-aligning, so a split that
/// supplies nothing cannot be its artefact — the reasoning
/// `delins-recommendation-reach-when-the-input-arrives-split` (decided) applies
/// to the sibling coalesce pass, and it holds here for the same reason. This is
/// W58's shape (`c.[992_993del;995_997del;999_1004del]`) at its smallest.
///
/// On `c.`, so the axis gate is satisfied and condition 3 is what declines —
/// see `an_equal_length_block_keeps_its_split`.
#[test]
fn an_all_deletion_split_stays_split() {
    let two = provider("NM_DEL.1", "CAGTCTCTCTCTCTCTCTCT", true);
    assert_eq!(
        normalized_fixed_point(&two, "NM_DEL.1:c.2_4delinsG"),
        "NM_DEL.1:c.[2del;4del]",
    );
    let three = provider("NM_DEL.2", "CAGTACTCTCTCTCTCTCTCT", true);
    assert_eq!(
        normalized_fixed_point(&three, "NM_DEL.2:c.2_5delinsGA"),
        "NM_DEL.2:c.[2del;4del]",
    );
}

/// A **net insertion** stays split — the direction scope of
/// `delins-merge-vs-individual-gap-two-or-more` (decided 2026-08-11), which
/// reaches the net-deletion case and does not reach net insertions, where the
/// split form remains canonical. It is also what keeps #422 and #999 outside
/// this rule entirely.
///
/// # WHAT THIS DOES AND DOES NOT GUARD — measured by sabotage, not assumed
///
/// Removing the length gate entirely leaves this test **green**. Its split is
/// held by condition 4, not condition 1: the second derived member is a `dup`
/// (`c.[2_3delinsAA;9dup]`), which `piece_renders_as_delins` excludes via
/// `is_tandem_duplication`. The guard for the direction scope is
/// `merge::tests::the_placed_gap_predicate_discriminates`, which does redden —
/// the same honest disclosure `an_equal_length_block_keeps_its_split` carries.
///
/// On `c.`, so the axis gate does not decline it first.
#[test]
fn a_net_insertion_stays_split() {
    let coding = provider("NM_TEST.1", CORE, true);
    assert_eq!(
        normalized_fixed_point(&coding, "NM_TEST.1:c.2_4delinsAACG"),
        "NM_TEST.1:c.[2_3delinsAA;9dup]",
    );
}

/// A member that is a lone **substitution** is a rank-1 type the split genuinely
/// buys, so the block stays split even though the other member is a placed gap.
///
/// The name reads as "an allele containing a substitution stays split", and that
/// is **not** what condition 4 says. It is read off the **derived** partition,
/// per `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]`,
/// so an input whose *authored* spelling carries a substitution can still merge
/// when its derivation does not — the conformance corpus has such rows. What
/// this pins is a block whose derivation yields a lone substitution.
///
/// On `c.`, so the axis gate does not decline it first.
#[test]
fn a_substitution_beside_a_placed_gap_stays_split() {
    let provider = provider("NM_SUB.1", "CAGTCTCTCTCTCTCTCTCT", true);
    assert_eq!(
        normalized_fixed_point(&provider, "NM_SUB.1:c.2_4delinsCG"),
        "NM_SUB.1:c.[2A>C;4del]",
    );
}
