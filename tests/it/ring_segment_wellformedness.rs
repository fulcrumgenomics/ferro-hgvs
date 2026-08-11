//! A ring's `::`-joined segments must actually describe a ring: listed in
//! `pter`→`qter` order, and each contributing a break junction.
//!
//! Follow-up to #1578 item 3. Adjudicating that `general.md:58`'s
//! self-cancellation rule does **not** cross `::`
//! (`rulings[self-cancelling-across-ring-junctions]`) left a hole it explicitly
//! declined to fill: ferro accepted `g.100_200del::150_250dup`,
//! `g.100_200del::300_400dup` and `g.150_250dup::100_200del`, none of which is a
//! well-formed ring. `:58` is the wrong instrument for them — it is an *overlap*
//! predicate, so it accepts the non-overlapping `del::300_400dup` however it is
//! wired. These are the right instruments.
//!
//! **Two rules, with the citations they rest on:**
//!
//! - **Ordering.** `DNA/complex.md:51` — "Break point location is determined by
//!   the first break point encountered, i.e. `pter` of the chromosome is to be
//!   listed first"; `:53` — "multiple breakpoints in one chromosome are listed in
//!   order of occurrence from `pter` to `qter`"; `:55` — "variant descriptions
//!   are always in the forward orientation (from `pter` to `qter`, the end of the
//!   chromosome)".
//! - **Junction-contributing segments.** `DNA/complex.md:39` — "a double colon
//!   is used to designate break point junctions creating a ring chromosome" —
//!   read with `:60-64`, where the committee withdrew `::` as a general join
//!   operator because it gave "one identical derivative chromosome … two
//!   different descriptions". A `dup`, `inv`, substitution or identity segment
//!   designates no junction, so it is not a thing `::` can join; the spec spells
//!   rearrangement-plus-local-edit composites as `;` cis alleles instead
//!   (`:113`, `:117`), where `general.md:58` governs them normally.
//!
//!   **Both deletion spellings count**, and this is the easy half to get wrong:
//!   `delN[15]` is a deletion whose removed size is a count of unknown bases, and
//!   an unsequenced breakpoint is precisely what the spec's ring examples describe
//!   ("breakpoint not sequenced", `:128`/`:163`). A first cut matched only
//!   `NaEdit::Deletion` and refused it.
//!
//! **What is deliberately NOT enforced: telomere anchoring.** A ring loses both
//! telomeres, so a well-formed ring's first segment starts at `pter` and its last
//! ends at `qter` — which is true of both ring shapes the spec publishes
//! (`:127`, `:161`). But no clause *says* so, and enforcing it would be ferro
//! legislating from two worked examples plus biology. Recorded as
//! `rulings[ring-telomere-anchoring]`, status `undecided`, and
//! `a_ring_with_no_telomere_anchor_is_still_accepted` pins that it stays
//! accepted so the open question does not get quietly closed by a later change.
//!
//! **Both rules are conservative: they reject only what is *definitely* wrong.**
//! Ordering compares an upper bound of the later segment's start against a lower
//! bound of the earlier one's, so an uncertain boundary is never enough to
//! reject — which is what keeps `:127`'s own `pter_(12200001_14700000)del::…`
//! accepted. Junction-contributing skips a segment carrying no edit to judge —
//! though note that is unreachable through the grammar (a `?` edit fails to parse
//! inside a `::`-join, measured), so it guards programmatic construction only and
//! is not exercised by a test here. A segment with no edit *at all* is a different
//! thing and **is** rejected: that is the withdrawn ISCN2016 join, pinned by
//! `a_position_only_segment_is_the_withdrawn_iscn2016_join_and_is_rejected`.

use ferro_hgvs::parse_hgvs;

/// The three shapes `rulings[self-cancelling-across-ring-junctions]` named as
/// malformed-but-accepted. This test is the reason this change exists, and it
/// replaces `no_ring_wellformedness_rule_yet_so_malformed_rings_are_still_accepted`
/// in `issue_1578_followup_self_cancelling_rings.rs`, which pinned the gap.
#[test]
fn the_three_malformed_rings_the_ruling_named_are_now_rejected() {
    for input in [
        // a `dup` segment designates no break junction
        "NC_000022.11:g.100_200del::150_250dup",
        // …and non-overlapping, so `general.md:58` could never have caught it
        "NC_000022.11:g.100_200del::300_400dup",
        // segments listed out of `pter`->`qter` order
        "NC_000022.11:g.150_250dup::100_200del",
    ] {
        parse_hgvs(input).expect_err(&format!(
            "`{input}` is not a well-formed ring and must be rejected"
        ));
    }
}

/// Ordering, on its own. `DNA/complex.md:53` requires breakpoints "listed in
/// order of occurrence from `pter` to `qter`", so a descending pair is wrong
/// regardless of edit type — which is why this is checked with `del` segments
/// that the junction rule accepts.
#[test]
fn ring_segments_listed_out_of_order_are_rejected() {
    for input in [
        "NC_000022.11:g.300_400del::100_200del",
        "NC_000022.11:g.37600001_qterdel::pter_14700000del",
        // three segments, only the last pair inverted
        "NC_000022.11:g.100_200del::500_600del::300_400del",
    ] {
        let error =
            parse_hgvs(input).expect_err(&format!("`{input}` lists its segments out of order"));
        let rendered = error.to_string();
        assert!(
            rendered.contains("pter") && rendered.contains("qter"),
            "the diagnostic must name the ordering rule it enforces, got: {rendered}"
        );
    }
}

/// Junction-contributing, on its own: every non-`del` edit kind ferro's ring
/// grammar accepts today, each in ascending order so the ordering rule cannot be
/// what rejects it.
#[test]
fn a_ring_segment_that_designates_no_break_junction_is_rejected() {
    for input in [
        "NC_000022.11:g.100_200del::300_400dup",
        "NC_000022.11:g.100_200del::300_400inv",
        "NC_000022.11:g.100_200del::300_400delinsAT",
        "NC_000022.11:g.100_200del::300A>T",
        "NC_000022.11:g.100_200del::300_400=",
        "NC_000022.11:g.100_200del::300_301insAT",
        // and when the offending segment comes first
        "NC_000022.11:g.100_200dup::300_400del",
        "NC_000022.11:g.100A>T::300C>G",
    ] {
        let error = parse_hgvs(input).expect_err(&format!(
            "`{input}` has a segment that designates no break junction"
        ));
        let rendered = error.to_string();
        assert!(
            rendered.contains("break point junction") || rendered.contains("junction"),
            "the diagnostic must name the junction rule it enforces, got: {rendered}"
        );
    }
}

/// The spec's own ring shapes must survive both rules, byte-identically. These
/// are the rows the change is measured against — `:127` in particular carries
/// **uncertain boundaries** on both segments, which is what forces the ordering
/// comparison to be conservative rather than a naive integer compare.
#[test]
fn the_specs_own_ring_shapes_still_round_trip() {
    for input in [
        "NC_000022.11:g.pter_(12200001_14700000)del::(37600001_410000000)_qterdel",
        "NC_000022.11:g.[pter_(12200001_14700000)del::(37600001_410000000)_qterdel]sup",
        // the same shape with certain bounds, and the `chr`-prefixed corpus row
        "NC_000022.11:g.pter_14700000del::37600001_qterdel",
        "chr22:g.pter_10000000del::40000000_qterdel",
    ] {
        let variant = parse_hgvs(input)
            .unwrap_or_else(|e| panic!("the spec publishes `{input}` as a ring: {e}"));
        assert_eq!(
            variant.to_string(),
            input,
            "`{input}` must round-trip byte-identically",
        );
    }
}

/// An ascending, all-`del` ring with no telomere anchor stays **accepted**: the
/// telomere-anchoring rule is `undecided`, and this pins that it was left open
/// rather than forgotten. See `rulings[ring-telomere-anchoring]`.
///
/// If a later change decides anchoring, this test is expected to fail — that is
/// the signal to move it into that rule's tests with the ruling flipped to
/// `decided`, not to relax it in place.
#[test]
fn a_ring_with_no_telomere_anchor_is_still_accepted() {
    for input in [
        "NC_000022.11:g.100_200del::300_400del",
        "NC_000022.11:g.pter_14700000del::37600001_40000000del",
        "NC_000022.11:g.100_14700000del::37600001_qterdel",
    ] {
        let variant = parse_hgvs(input).unwrap_or_else(|e| {
            panic!(
                "`{input}` is ascending and all-`del`; telomere anchoring is \
                 `undecided` (rulings[ring-telomere-anchoring]) so it must still \
                 parse: {e}"
            )
        });
        assert_eq!(variant.to_string(), input, "`{input}` must round-trip");
    }
}

/// Overlap between `del` segments is **not** what these rules judge, and stays
/// accepted. It is the question `rulings[self-cancelling-across-ring-junctions]`
/// answered in the negative, so rejecting it here would re-open a decided
/// ruling through the back door — the exact failure mode that ruling's guard was
/// written to prevent.
#[test]
fn overlapping_del_segments_are_not_what_these_rules_judge() {
    for input in [
        "NC_000022.11:g.100_200del::150_250del",
        "NC_000022.11:g.100_200del::100_200del",
    ] {
        let variant = parse_hgvs(input).unwrap_or_else(|e| {
            panic!(
                "`{input}` is ascending and all-`del`; overlap is governed by \
                 rulings[self-cancelling-across-ring-junctions], which ruled it \
                 out of scope for `::`: {e}"
            )
        });
        assert_eq!(variant.to_string(), input, "`{input}` must round-trip");
    }
}

/// Ring rules must reach a `sup`-wrapped ring too — the same delivery #1578's
/// walker established for the edit-form validators.
#[test]
fn the_rules_reach_a_sup_wrapped_ring() {
    for input in [
        "NC_000022.11:g.[100_200del::300_400dup]sup",
        "NC_000022.11:g.[300_400del::100_200del]sup",
    ] {
        parse_hgvs(input).expect_err(&format!(
            "`{input}`: a `sup`-wrapped ring is still a ring and must be checked"
        ));
    }
}

/// `delN[15]` is a **deletion** — one whose removed size is given as a count of
/// unknown bases — so it designates a break point junction and must be accepted.
///
/// This was a live defect in the first cut of the junction rule, which matched
/// only `NaEdit::Deletion` and so refused `NaEdit::NPaddedDeletion`. It is the
/// worst possible case to get wrong here: both ring shapes the spec publishes are
/// annotated *"breakpoint not sequenced"* (`DNA/complex.md:128`, `:163`), and an
/// unsequenced breakpoint is precisely what `delN[]` spells. The rule was
/// rejecting the form the spec's own examples describe in prose.
///
/// It also produced a diagnostic that misdescribed the input, calling it "a
/// duplication, inversion, substitution or identity segment" when it was none of
/// those.
#[test]
fn an_n_padded_deletion_segment_designates_a_junction() {
    for input in [
        "NC_000022.11:g.100_200delN[15]::300_400del",
        "NC_000022.11:g.100_200del::300_400delN[15]",
        "NC_000022.11:g.100_200delN[15]::300_400delN[20]",
        "NC_000022.11:g.[100_200delN[15]::300_400del]sup",
    ] {
        let variant = parse_hgvs(input).unwrap_or_else(|e| {
            panic!(
                "`{input}`: `delN[]` is a deletion with an unsequenced size, which \
                 is exactly what the spec's ring examples describe \
                 (\"breakpoint not sequenced\"): {e}"
            )
        });
        assert_eq!(variant.to_string(), input, "`{input}` must round-trip");
    }
}

/// A **position-only** segment — an interval with no edit at all — is rejected,
/// and this is the shape that most directly vindicates the junction rule.
///
/// It is the withdrawn ISCN2016 join form. `DNA/complex.md:72` marks
/// `NC_000002.12:g.pter_8247756::NC_000011.10:g.15825273_cen_qter` as
/// `<code class="invalid">` and says it was "corrected" to the `delins` format;
/// `:84` explains what such a coupling used to mean — "coupling
/// `chr11:108111981` to `108111987` implies nucleotides `108111982_108111986` are
/// deleted", i.e. the deletion was *implied* by the join rather than stated. That
/// implicit-deletion reading is what `:60-64` withdrew, on the ground that one
/// derivative chromosome was getting two descriptions.
///
/// So a position-only `::` segment is not an under-specified ring — it is the
/// notation the committee removed.
#[test]
fn a_position_only_segment_is_the_withdrawn_iscn2016_join_and_is_rejected() {
    for input in [
        "NC_000022.11:g.pter_100::200_qter",
        "NC_000022.11:g.100_200::300_400",
    ] {
        let error = parse_hgvs(input).expect_err(&format!(
            "`{input}` is the withdrawn ISCN2016 join form (DNA/complex.md:72, :84)"
        ));
        assert!(
            error.to_string().contains("junction"),
            "the diagnostic must name the junction rule, got: {error}"
        );
    }
}
