//! Issue #1284 — the del+dup collision repair was genomic-only, so the
//! CDS-relative axes emitted descriptions ferro's own parser rejects.
//!
//! `respell_colliding_duplications` (and the two sibling repairs added behind
//! the same gate) read the duplicated bases *over the member's own
//! coordinates*, which is only sound on an axis whose 1-based position is a
//! direct offset into the sequence the provider serves. That is true of `g.`,
//! `m.` and — since #1315 — `n.`, and false of `c.`/`r.`, which are
//! CDS-relative. The gate therefore refused those axes, and a del+dup pair that
//! 3'-shifted onto one position stayed collided:
//!
//! ```text
//! NM_TEST.1:c.[-8del;-6dup]  ->  NM_TEST.1:c.[-1del;-1dup]
//! parse_hgvs(out) -> Err: Self-cancelling allele … overlapping reference positions
//! ```
//!
//! ## Why this is an abort rather than a diff
//!
//! `FERRO_ASSERT_REPARSE=1` is armed in `ci.yml` and `nightly-mutalyzer.yml`.
//! That oracle asserts normalization emits something `parse_hgvs` accepts, so
//! every case below **panics the run** under it rather than failing with a
//! comparison. That is why #1283's transcript-axis sweep gap is blocked on this
//! issue, and why the invariant — not any single expected string — is the
//! contract these tests pin.
//!
//! ## The conversion, and why the earlier attempt did not work
//!
//! The issue records a naive `axis_sequence_delta` (a single `cds_start - 1`)
//! being tried and reverted: it is right at `c.1` and off by one below, because
//! the CDS axis has no zero. The fix here keys the offset on the member's
//! **`Region`** instead of on the axis kind, which removes the discontinuity
//! rather than patching across it — `FivePrimeUtr`, `Cds` and `ThreePrimeUtr`
//! are each individually affine onto the transcript, and every pass already
//! groups members within one region.
//!
//! Reference-free (`MockProvider` via `SyntheticBuilder`), so these hold with
//! no manifest.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// The issue's own transcript: 5'UTR is a T-run reaching the CDS start, so both
/// members shuffle to `c.-1` and collide at the UTR/CDS boundary — the hardest
/// placement, because the junction 3' of `c.-1` is `c.-1_1`, spanning two
/// regions.
const ISSUE_CORE: &str = "TTTTTTTTTAATATATTTTAATAT";
const ISSUE_CDS: (u64, u64) = (9, 24);

/// A 5'UTR with interior structure (an A-run at tx 3..5, well inside the UTR)
/// plus a real 3'UTR, so a collision can settle away from either boundary.
const STRUCTURED_CORE: &str = "GCAAAGCGCGCGATGAAACCCTAAGGCATTTTTAA";
const STRUCTURED_CDS: (u64, u64) = (13, 24);

fn provider(core: &str, cds: (u64, u64)) -> MockProvider {
    SyntheticBuilder::cds(core, cds.0, cds.1, Strand::Plus).build()
}

/// Normalize, and return the rendered output.
fn normalize(input: &str, core: &str, cds: (u64, u64)) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::new(provider(core, cds))
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string()
}

/// The contract: whatever normalization emits, `parse_hgvs` must accept it.
///
/// Asserted through the real parser rather than by string-matching the known
/// `Self-cancelling` message, so a future collision that renders some *other*
/// unparseable shape is caught by the same test.
fn assert_reparses(input: &str, core: &str, cds: (u64, u64)) -> String {
    let output = normalize(input, core, cds);
    if let Err(err) = parse_hgvs(&output) {
        panic!("`{input}` normalized to `{output}`, which ferro cannot re-parse: {err}");
    }
    output
}

/// The issue's verbatim reproducer.
#[test]
fn the_issues_own_reproducer_reparses() {
    assert_reparses("NM_TEST.1:c.[-8del;-6dup]", ISSUE_CORE, ISSUE_CDS);
}

/// The same shape on `r.`, which shares the CDS-relative axis and was gated
/// identically.
///
/// **Deliberately not run on `ISSUE_CORE`.** It passes there — and proves
/// nothing. That transcript's 5'UTR is a T-run reaching the CDS start, and
/// `normalize_rna` does not apply the CDS↔UTR axis clamp that `normalize_cds`
/// does (an explicit #334 decision, on the grounds that `r.` natively spans the
/// sub-axes; see the comment on `normalize_rna`'s `boundaries` and the test it
/// names, `issue_163_rna_utr3_flag::rna_mixed_cds_utr3_del_shifts_into_utr`).
/// So on `ISSUE_CORE` the `r.` members shuffle *out* of the 5'UTR into the CDS
/// body, where the sequence-first canonicalization resolves them — the repair
/// under test never runs, and the assertion is vacuous. `STRUCTURED_CORE`'s
/// A-run sits well inside the 5'UTR, so the members stay there and the
/// CDS-relative conversion is what has to be right.
#[test]
fn the_reproducer_reparses_on_the_rna_axis() {
    for input in [
        "NM_TEST.1:r.[-10del;-9dup]",
        "NM_TEST.1:r.[-9del;-8dup]",
        "NM_TEST.1:r.[-10del;-8dup]",
    ] {
        assert_reparses(input, STRUCTURED_CORE, STRUCTURED_CDS);
    }
}

/// Characterizes the `c.`/`r.` divergence the test above works around, so it is
/// recorded in code rather than only in a PR body.
///
/// The same molecular event — delete one base from a run that reaches the CDS
/// start — settles at `c.-1` on the `c.` axis and at `r.1` on the `r.` axis,
/// because only `normalize_cds` intersects the shuffle window with the CDS↔UTR
/// axis bound. Deliberate (#334) and out of scope for #1284, which is about the
/// repair passes' coordinate conversion, not about where the shuffle stops.
/// Pinned as *current* behavior: if #334 is ever revisited this test should
/// move with it, and until then it stops the divergence being rediscovered as a
/// surprise.
#[test]
fn the_cds_and_rna_axes_stop_a_five_prime_shuffle_in_different_places() {
    assert_eq!(
        normalize("NM_TEST.1:c.-8del", ISSUE_CORE, ISSUE_CDS),
        "NM_TEST.1:c.-1del",
        "`c.` clamps the shuffle at the CDS start"
    );
    assert_eq!(
        normalize("NM_TEST.1:r.-8del", ISSUE_CORE, ISSUE_CDS),
        "NM_TEST.1:r.1del",
        "`r.` does not (#334), so the same event crosses into the CDS body"
    );
}

/// A collision settling in the *interior* of the 5'UTR, where the repair's
/// junction is wholly inside one region. Distinct from the reproducer above,
/// whose junction straddles the UTR/CDS boundary: an implementation that
/// handled only the interior case would still abort on the issue's own input,
/// and one that special-cased only the boundary would miss these.
#[test]
fn an_interior_five_prime_utr_collision_reparses() {
    for input in [
        "NM_TEST.1:c.[-10del;-9dup]",
        "NM_TEST.1:c.[-10del;-8dup]",
        "NM_TEST.1:c.[-9dup;-8del]",
    ] {
        assert_reparses(input, STRUCTURED_CORE, STRUCTURED_CDS);
    }
}

/// The 3'UTR half of the same axis. Its offset is `cds_end`, not `cds_start`,
/// so it is a genuinely different conversion and not covered by the 5' cases.
#[test]
fn a_three_prime_utr_collision_reparses() {
    for input in ["NM_TEST.1:c.[*1del;*2dup]", "NM_TEST.1:r.[*1del;*2dup]"] {
        assert_reparses(input, STRUCTURED_CORE, STRUCTURED_CDS);
    }
}

/// The census. Every two-member combination over both UTRs **and** the CDS
/// body, on both CDS-relative axes, must re-parse.
///
/// Measured on this grid: **120 of 12,110 outputs failed to re-parse before the
/// fix, 0 after**, every one of them the `Self-cancelling` del+dup collision.
/// (A first grid, whose 5'UTR was a single homopolymer reaching the CDS start,
/// put 296 of 3080 in the same one shape.) The CDS body contributed zero
/// failures either way — there the sequence-first canonicalization resolves the
/// collision before it can be rendered, which is exactly why the defect was
/// confined to the UTRs.
///
/// The grid pairs coordinates **across** regions as well as within them. Same
/// region only was 3,950 cases and the same 120 failures, so cross-region pairs
/// contribute none — which is the expected answer rather than a wasted 8,000
/// cases: a collision needs both members on one position, and members in
/// different regions cannot reach one. Covering them is what lets the census
/// claim that, instead of leaving it assumed.
///
/// Two things are asserted besides re-parsing. A sweep that silently generated
/// nothing would pass vacuously, so the attempted count is checked against the
/// measured population; and each output is normalized a **second** time, because
/// these repairs rewrite members in place and a repair that is not a fixed point
/// is what the idempotency invariant exists to catch. That second pass is
/// deliberately not left to `FERRO_ASSERT_IDEMPOTENT`: CI is the only place that
/// oracle is armed, so a plain local `cargo nextest run --features dev` got no
/// idempotency signal from this grid at all.
#[test]
fn no_transcript_axis_two_member_allele_fails_to_reparse() {
    let utr5: Vec<String> = (1..=12).map(|i| format!("-{i}")).collect();
    let utr3: Vec<String> = (1..=11).map(|i| format!("*{i}")).collect();
    let body: Vec<String> = (1..=12).map(|i| i.to_string()).collect();

    let mut attempted = 0usize;
    let mut failures: Vec<String> = Vec::new();

    // Cross-region pairs as well as same-region ones. A pair drawn from one
    // region can only collide inside that region, so restricting to those left
    // the whole cross-region class — `c.[-8del;1dup]`, a 5'UTR member beside a
    // CDS-body one — outside a census that reads as exhaustive. Those are the
    // pairs whose two members take *different* `region_sequence_delta` values,
    // which is the arithmetic this PR adds.
    let regions = [&utr5, &utr3, &body];
    for first_coords in regions {
        for second_coords in regions {
            for axis in ["c", "r"] {
                for a in first_coords {
                    for b in second_coords {
                        for (first, second) in [
                            (format!("{a}del"), format!("{b}dup")),
                            (format!("{a}dup"), format!("{b}del")),
                            (format!("{a}del"), format!("{b}del")),
                            (format!("{a}dup"), format!("{b}dup")),
                            (format!("{a}A>G"), format!("{b}dup")),
                        ] {
                            let input = format!("NM_TEST.1:{axis}.[{first};{second}]");
                            // Not every generated string is a legal input (a stated
                            // reference base that does not match, say); those are
                            // simply not part of the population under test.
                            let Ok(variant) = parse_hgvs(&input) else {
                                continue;
                            };
                            let Ok(normalized) =
                                Normalizer::new(provider(STRUCTURED_CORE, STRUCTURED_CDS))
                                    .normalize(&variant)
                            else {
                                continue;
                            };
                            attempted += 1;
                            let output = normalized.to_string();
                            let Ok(reparsed) = parse_hgvs(&output) else {
                                failures.push(format!("{input} -> {output} (does not re-parse)"));
                                continue;
                            };
                            // Idempotency, asserted here rather than left to
                            // `FERRO_ASSERT_IDEMPOTENT`. These repairs rewrite
                            // members in place, so a repair that is not a fixed
                            // point is the failure mode the invariant exists to
                            // catch — and CI is the only place that oracle is armed,
                            // so a plain local `cargo nextest run --features dev`
                            // got no idempotency signal at all from this 12,110-case
                            // grid.
                            //
                            // An `Err` on the second pass is recorded, not skipped.
                            // The first pass succeeded and its output re-parsed, so a
                            // re-normalization that now refuses is itself the defect —
                            // skipping it would hide exactly the case this check adds.
                            match Normalizer::new(provider(STRUCTURED_CORE, STRUCTURED_CDS))
                                .normalize(&reparsed)
                            {
                                Ok(twice) if twice.to_string() == output => {}
                                Ok(twice) => failures.push(format!(
                                    "{input} -> {output} -> {twice} (not a fixed point)"
                                )),
                                Err(err) => failures.push(format!(
                                    "{input} -> {output} (re-normalizing the output failed: {err})"
                                )),
                            }
                        }
                    }
                }
            }
        }
    }

    // Close to the measured 12,110 rather than a token floor. At `> 1000` roughly
    // nine tenths of the grid could start erroring out of `normalize` before this
    // test noticed, and the `continue` on `Err` above makes that silent.
    assert!(
        attempted > 11_000,
        "the sweep must actually normalize a population; only {attempted} inputs survived"
    );
    assert!(
        failures.is_empty(),
        "{} of {attempted} normalized outputs do not re-parse; first 10:\n{}",
        failures.len(),
        failures
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// Anti-vacuity control: the grid above really does contain colliding pairs, so
/// the census is measuring something. Pinned on a pair whose two members
/// 3'-shift onto one position.
#[test]
fn the_swept_grid_really_does_produce_collisions() {
    let collided = normalize(
        "NM_TEST.1:c.[-10del;-9dup]",
        STRUCTURED_CORE,
        STRUCTURED_CDS,
    );
    assert!(
        !collided.contains("-10"),
        "precondition: the members must 3'-shift together before colliding; got `{collided}`"
    );
    // Movement is not collision. The above holds whenever *either* member
    // shifts at all, so on its own it does not support the test's name. The
    // collision itself is directly observable: the `del` and the `dup` land on
    // one position, cancel, and the two-member allele collapses to a single
    // identity member. A pair that merely shifted would still be rendered as
    // `[...;...]`.
    assert_eq!(
        collided, "NM_TEST.1:c.-9=",
        "precondition: the pair must actually collide and cancel to one identity          member, which is what the repair produces; got `{collided}`"
    );
}

/// The conversion must not depend on strand.
///
/// `cds_start`/`cds_end` are transcript coordinates and the provider serves
/// transcript-oriented bases, so a minus-strand record should convert
/// identically — but "wrong shift direction or strand handling" is exactly the
/// class this file's arithmetic could get wrong, and every other test here runs
/// on the plus strand. Asserted as *parity* rather than against pinned strings,
/// so it keeps holding if the settled forms move.
#[test]
fn the_conversion_is_strand_independent() {
    for input in [
        "NM_TEST.1:c.[-10del;-9dup]",
        "NM_TEST.1:c.[-10del;-8dup]",
        "NM_TEST.1:c.[*1del;*2dup]",
        "NM_TEST.1:r.[-9del;-8dup]",
    ] {
        let plus = normalize(input, STRUCTURED_CORE, STRUCTURED_CDS);
        let variant = parse_hgvs(input).expect("input must parse");
        let minus = Normalizer::new(
            SyntheticBuilder::cds(
                STRUCTURED_CORE,
                STRUCTURED_CDS.0,
                STRUCTURED_CDS.1,
                Strand::Minus,
            )
            .build(),
        )
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string();
        assert_eq!(
            plus, minus,
            "`{input}` must normalize the same on either strand"
        );
    }
}

/// A record whose CDS bounds are inverted (`cds_end < cds_start`) has an
/// incoherent `c.` axis — its 5'UTR and 3'UTR halves overlap, so one transcript
/// position carries two names. The repair must refuse rather than pick one.
///
/// This does not claim the ambiguity itself is fixed: it predates this change
/// and is reachable without it (`c.[*1del;*2dup]` comes back spelled `c.[-10del;
/// -9dup]` on `origin/main` too, both naming transcript position 14). What is
/// pinned is that the new conversion declines instead of adding another site
/// that resolves it arbitrarily — and that declining still leaves a
/// re-parseable output.
#[test]
fn inverted_cds_bounds_are_refused_not_guessed() {
    let inverted = || SyntheticBuilder::cds(STRUCTURED_CORE, 24, 13, Strand::Plus).build();
    let variant = parse_hgvs("NM_TEST.1:c.[*1del;*2dup]").expect("input must parse");
    let output = Normalizer::new(inverted())
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string();
    assert!(
        parse_hgvs(&output).is_ok(),
        "even on a malformed record the output must re-parse; got `{output}`"
    );
    // Refusal, asserted directly. `is_ok()` alone passes against the exact
    // behaviour this test excludes: a repair that *guesses* one of the two names
    // an inverted-bounds record admits still emits a parseable string.
    //
    // `ordered_cds_bounds` returning `None` is observable as "no member is
    // respelled": both keep their authored `del`/`dup` shapes and the allele
    // stays two members. Contrast the well-formed record, where the same input
    // collapses to the single identity member `c.*2=`.
    assert_eq!(
        output, "NM_TEST.1:c.[-10del;-9dup]",
        "inverted CDS bounds must make the conversion decline and leave both members \
         authored, not pick a spelling; got `{output}`"
    );

    // Every CDS-relative region must refuse alike, `r.` included. Its
    // `region_sequence_delta` arm falls back to a transcript-relative delta of 0
    // when the record has no CDS at all — right for a non-coding transcript —
    // and that fallback must not also swallow an *inverted* CDS, which would turn
    // the refusal into exactly the silent guess this test excludes. The 5'UTR and
    // CDS-body arms had the mirror-image gap: they read `cds_start` alone, so
    // only the 3'UTR declined while those two picked a spelling.
    let rna = parse_hgvs("NM_TEST.1:r.[*1del;*2dup]").expect("input must parse");
    let rna_output = Normalizer::new(inverted())
        .normalize(&rna)
        .expect("lenient normalization must not reject")
        .to_string();
    assert_eq!(
        rna_output, "NM_TEST.1:r.[-10del;-9dup]",
        "the r. axis must decline on an inverted record too; got `{rna_output}`"
    );
}

/// A collision that settles on the **final transcript base** is spelled in
/// range — the #1343 gap, closed by #1327's terminal re-spelling.
///
/// `c.[*10dup;*11dup]` 3'-shifts both members onto `*11`, the last base of this
/// transcript. Anchoring the repair's insertion at the junction 3' of it would
/// need `*12`, which the record does not have: `cds_axis_position` has a floor
/// (`tx < 1`) and no ceiling, so it answered `ThreePrimeUtr(12)` rather than
/// refusing, and the result rendered as `c.*11_*12insAA` — readable, so
/// `FERRO_ASSERT_REPARSE` could not see it, and wrong.
///
/// `respell_at_sequence_end` now spells the terminal case as the duplication it
/// is, which needs no position past the end. Three answers for one input, in
/// order:
///
/// ```text
/// origin/main   c.[*11dup;*11dup]   two members on one position (#1234 class)
/// #1284 alone   c.*11_*12insAA      names a base the transcript does not have
/// with #1327    c.*10_*11A[4]       in range, one member, same bases
/// ```
#[test]
fn a_collision_on_the_final_transcript_base_is_spelled_in_range() {
    // 3'UTR is `*1..*11` (transcript 25..35 of a 35-base record), so `*12` is
    // one past the end and must not appear.
    for input in ["NM_TEST.1:c.[*10dup;*11dup]", "NM_TEST.1:c.[*11dup;*10dup]"] {
        let output = normalize(input, STRUCTURED_CORE, STRUCTURED_CDS);
        assert_eq!(
            output, "NM_TEST.1:c.*10_*11A[4]",
            "`{input}` must be spelled without a position past the transcript end"
        );
        assert!(
            !output.contains("*12"),
            "`{input}` -> `{output}` names a position past the transcript end"
        );
        assert!(parse_hgvs(&output).is_ok(), "`{output}` must re-parse");
    }
}

/// The other terminal shape: a collision at the last base whose two members
/// **cancel**.
///
/// `*10` and `*11` are both `A`, so deleting the first and duplicating the second
/// cancel exactly and the allele denotes the reference. On `origin/main` this
/// settled as the unparseable `c.[*11del;*11dup]` — one of the 120 — and #1284
/// resolves it to the single identity member asserted here.
///
/// **This is the guard that caught #1344.** The terminal re-spelling #1327 adds
/// is right for the `dup`+`dup` case above, but applied here it turned the
/// cancellation into `c.[*11del;*11delinsAA]` — two coincident edits on one
/// position, which ferro's own strict mode rejects (`OverlapConflictingEdits /
/// W5002`) and which was a stable fixed point, so a permanent answer rather than
/// a transient one.
///
/// The cause was that the boundary identity trades a zero-width junction
/// insertion for a `delins` that *claims* the last base, and those are not
/// interchangeable when a sibling already claims it: an insertion flush against
/// a deleted base is the #999 adjacency the collapse pass merges, two members on
/// one base is an overlap it cannot. `respell_colliding_duplications` now refuses
/// the re-spelling in that case, which leaves the junction form for the merge to
/// consume — so the out-of-range coordinate #1327 was avoiding still never
/// reaches the output.
#[test]
fn a_cancelling_collision_on_the_final_transcript_base_is_repaired() {
    for input in ["NM_TEST.1:c.[*10del;*11dup]", "NM_TEST.1:r.[*10del;*11dup]"] {
        let output = assert_reparses(input, STRUCTURED_CORE, STRUCTURED_CDS);
        let axis = if input.contains(":c.") { "c" } else { "r" };
        assert_eq!(
            output,
            format!("NM_TEST.1:{axis}.*11="),
            "`{input}` must collapse to one identity member"
        );
    }
}

/// Shapes that already worked must not move. The CDS body is resolved by the
/// sequence-first canonicalization rather than by this repair, and a
/// non-colliding UTR pair has nothing to repair — both are regressions if they
/// change.
#[test]
fn already_correct_shapes_are_unchanged() {
    for (input, expected) in [
        // Non-colliding 5'UTR pair: distinct positions, no shuffle onto each other.
        ("NM_TEST.1:c.[-6del;-5dup]", "NM_TEST.1:c.[-6del;-5dup]"),
        ("NM_TEST.1:c.[-2del;-1dup]", "NM_TEST.1:c.[-2del;-1dup]"),
    ] {
        assert_eq!(
            normalize(input, STRUCTURED_CORE, STRUCTURED_CDS),
            expected,
            "`{input}` must be left alone"
        );
    }
}
