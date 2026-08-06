//! #1431 — a single-position repeat anchor was read two ways, and the two
//! readings denoted different sequences.
//!
//! `hgvs_to_spdi` read `g.263A[7]` as *"the one base named at 263 becomes 7
//! copies"* — a deletion span of one — so on a 3-base `A` tract the two
//! untouched tract bases survived and the description denoted **nine** `A`s.
//! `merge::lowered_repeat` documents the opposite reading, explicitly:
//!
//! > A single-position anchor (`e == s`) names only the tract's start and means
//! > "the whole run becomes N", so it absorbs nothing.
//!
//! Under that reading the answer is seven, matching the range spelling
//! `g.263_265A[7]`. Two sites, two answers, one variant.
//!
//! **The spec settles it.** `DNA/repeated.md` presents the two spellings as two
//! *formats of one variant*:
//!
//! > a Community Consultation proposal is being prepared which will suggest to
//! > allow only the format where the **entire range** of the repeated sequence
//! > is indicated; so `g.123_191CAG[23]`, **not** `g.123CAG[23]`.
//!
//! `123_191` is 69 bases — the whole 23-copy tract — so the start-only form
//! addresses the run, not one base of it.
//!
//! Two things were wrong, not one. Besides the arithmetic above, the one-base
//! reading did not survive a multi-base unit **at all**: `g.263CAG[5]` fetched a
//! single base and died on the divisibility check ("repeat span length 1 is not
//! a multiple of unit length 3"), so the spec's own start-only format was
//! unusable for every unit the spec illustrates it with.
//!
//! Why it matters beyond the arithmetic: it made the apply oracle unreliable for
//! single-position repeats, and that oracle is the independent check several
//! sweeps and regression tests grade the normalizer against. It produced a
//! well-formed output denoting different bases, caught by none of
//! `FERRO_ASSERT_REPARSE` (both parse), `FERRO_ASSERT_IN_BOUNDS` (both in range)
//! or `FERRO_ASSERT_IDEMPOTENT` (the output is a fixed point).

use crate::common::cis_apply_oracle::apply;
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, MockProvider};

/// 256 `N`s, so core position 1 is `g.257` and the `N` flanks can never look
/// like part of a tract — `N` matches no repeat unit used here.
const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

fn padded(core: &str) -> String {
    format!("{PAD}{core}{PAD}")
}

fn provider(sequence: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", sequence.to_string());
    provider
}

/// The SPDI triple, rendered.
fn spdi(sequence: &str, input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    hgvs_to_spdi(&variant, &provider(sequence))
        .unwrap_or_else(|e| panic!("`{input}` must convert to SPDI: {e}"))
        .to_string()
}

/// The core bases a description denotes, with the `N` padding trimmed off.
fn denotes(sequence: &str, input: &str) -> String {
    apply(sequence, input)
        .unwrap_or_else(|| panic!("`{input}` must apply"))
        .trim_matches('N')
        .to_string()
}

/// The headline: the two spellings of one repeat denote the same sequence.
///
/// A 3-base `A` tract at 263-265. Both descriptions say "this `A` tract has
/// seven copies", and both must mean it.
#[test]
fn the_two_spellings_of_one_repeat_denote_one_sequence() {
    let sequence = padded("GATCATAAATTCAGC");
    let anchored = denotes(&sequence, "TEMPLATE:g.263A[7]");
    let ranged = denotes(&sequence, "TEMPLATE:g.263_265A[7]");
    assert_eq!(
        anchored, ranged,
        "the start-only and full-range spellings of one repeat must denote one sequence"
    );
    // Pinned, not merely equal: equality alone would be satisfied if both
    // regressed to nine `A`s together.
    assert_eq!(
        anchored, "GATCATAAAAAAATTCAGC",
        "seven `A`s, not the nine the one-base reading produced"
    );
    // The *tract*, not every `A` in the core — `GATCAT…TTCAGC` has three more.
    assert_eq!(
        longest_run(&anchored, 'A'),
        7,
        "the tract must hold exactly the seven copies the description states"
    );
}

/// …and the triples themselves are identical, not merely equivalent.
///
/// The width of the SPDI deletion is what decides overlap
/// (`triples_are_disjoint`), so two triples that denote the same bases at
/// different widths still disagree about whether an allele is well-formed. That
/// is the failure in `an_anchored_repeat_and_its_range_agree_on_overlap` below,
/// and this pins the property that removes it.
#[test]
fn the_two_spellings_produce_the_same_triple() {
    let sequence = padded("GATCATAAATTCAGC");
    assert_eq!(
        spdi(&sequence, "TEMPLATE:g.263A[7]"),
        spdi(&sequence, "TEMPLATE:g.263_265A[7]"),
    );
    assert_eq!(
        spdi(&sequence, "TEMPLATE:g.263A[7]"),
        "TEMPLATE:262:AAA:AAAAAAA",
        "the deletion must span the whole 3-base tract, not the anchor base alone"
    );
}

/// A multi-base unit — the shape the spec actually illustrates the start-only
/// format with (`g.123CAG[23]`).
///
/// This did not merely give a wrong answer before; it was a hard error, because
/// a one-base deletion span is not divisible by a 3-base unit. So the format the
/// spec documents was unusable for the unit the spec documents it with.
#[test]
fn a_multi_base_unit_works_from_the_anchor_alone() {
    // A 3-copy `CAG` tract at 263-271.
    let sequence = padded("GATCATCAGCAGCAGTTC");
    let anchored = denotes(&sequence, "TEMPLATE:g.263CAG[5]");
    let ranged = denotes(&sequence, "TEMPLATE:g.263_271CAG[5]");
    assert_eq!(
        anchored, ranged,
        "the start-only spelling of a multi-base repeat must reach the range spelling's answer"
    );
    assert_eq!(anchored, "GATCATCAGCAGCAGCAGCAGTTC");
    assert_eq!(
        anchored.matches("CAG").count(),
        5,
        "five copies, as the description states"
    );
}

/// A range that names only *part* of a tract is unchanged: it still states which
/// copies the count applies to, and the copies outside it survive.
///
/// This is the discriminating half. Widening every repeat to its physical tract
/// would break this case, and it would pass every assertion above.
#[test]
fn a_partial_range_still_absorbs_only_its_own_copies() {
    let sequence = padded("GATCATCAGCAGCAGTTC");
    // `263_265` names the first of three `CAG` copies. Five copies there plus
    // the two the range does not name is seven, not five.
    assert_eq!(
        denotes(&sequence, "TEMPLATE:g.263_265CAG[5]"),
        "GATCATCAGCAGCAGCAGCAGCAGCAGTTC"
    );
    assert_eq!(
        spdi(&sequence, "TEMPLATE:g.263_265CAG[5]"),
        "TEMPLATE:262:CAG:CAGCAGCAGCAGCAG",
        "an explicit range must keep its own width"
    );
}

/// The overlap verdict now agrees between the two spellings.
///
/// The deletion width decides `triples_are_disjoint`, so before this fix
/// `g.[263A[7];264dup]` was *admitted* — its repeat triple was one base wide and
/// so missed the duplication at 264 — while `g.[263_265A[7];264dup]`, the same
/// variant, was declined as overlapping. One variant, two spellings, opposite
/// answers to "is this even well-formed".
///
/// Both decline now, which is the right verdict: the repeat spans 263-265 and
/// the duplication claims 264, so they genuinely collide.
#[test]
fn an_anchored_repeat_and_its_range_agree_on_overlap() {
    let sequence = padded("GATCATAAATTCAGC");
    assert_eq!(
        apply(&sequence, "TEMPLATE:g.[263A[7];264dup]"),
        None,
        "a repeat whose tract contains its sibling must be declined as overlapping"
    );
    assert_eq!(
        apply(&sequence, "TEMPLATE:g.[263_265A[7];264dup]"),
        None,
        "and the range spelling of that same allele must get the same verdict"
    );
}

/// An anchor with no tract under it is untouched: the unit does not match the
/// reference there, and the existing "does not match repeat unit" diagnostic is
/// what must report it — not a second error introduced by the tract search.
#[test]
fn an_anchor_with_no_matching_tract_still_reports_the_mismatch() {
    let sequence = padded("GATCATAAATTCAGC");
    // Position 257 is `G`; a `C` unit does not match it.
    let variant = parse_hgvs("TEMPLATE:g.257C[4]").expect("parse");
    let error = hgvs_to_spdi(&variant, &provider(&sequence))
        .expect_err("a unit that does not match the reference must be refused");
    let message = error.to_string();
    assert!(
        message.contains("does not match repeat unit"),
        "the mismatch must be reported by the existing diagnostic; got: {message}"
    );
}

/// A single isolated copy is a one-unit tract, so the anchored spelling keeps
/// the width it always had. Guards the boundary between "no tract" and "a tract
/// of one".
#[test]
fn a_lone_copy_is_a_one_unit_tract() {
    // Core `GATCATAAATTCAGC` starts at `g.257`, so its 14th base — the lone `G`
    // between `A` and `C` — is `g.270`.
    let sequence = padded("GATCATAAATTCAGC");
    assert_eq!(
        spdi(&sequence, "TEMPLATE:g.270G[3]"),
        "TEMPLATE:269:G:GGG",
        "a lone copy must keep its one-unit width"
    );
}

/// A tract longer than the initial search window is found whole, for a unit
/// length that does not divide the window evenly.
///
/// The window doubles while the run still reaches an edge, and both edge tests
/// have to be stated in **units** rather than bytes: `count_tandem_repeats`
/// steps by `unit_len`, so a clipped run stops up to `unit_len - 1` bytes short
/// of the edge instead of on it. A byte-equality test therefore only fires when
/// the bytes remaining happen to be a multiple of the unit — which at the
/// initial half-width is true for a 1- and 3-base unit and false for a 4- and
/// 5-base one.
///
/// That asymmetry is why this guard uses `CAGT`: with a byte-equality test the
/// 3-base units the spec illustrates (`g.123_191CAG[23]`) pass while a 4-base
/// unit silently returns the clipped span — 128 bases of a 240-base tract,
/// which is itself a multiple of 4 and so clears both the divisibility and
/// unit-match checks downstream. A well-formed triple denoting the wrong bases,
/// which is the very defect this whole change removes.
#[test]
fn a_tract_longer_than_the_search_window_is_found_whole() {
    const COPIES: usize = 60;
    let core = "CAGT".repeat(COPIES);
    let sequence = padded(&core);
    let tract_end = 257 + core.len() - 1;

    // The anchored spelling and the explicit-range spelling must agree, and
    // both must denote five copies rather than five copies plus the tail of an
    // under-counted tract.
    let anchored = denotes(&sequence, "TEMPLATE:g.257CAGT[5]");
    let ranged = denotes(&sequence, &format!("TEMPLATE:g.257_{tract_end}CAGT[5]"));
    assert_eq!(
        anchored, ranged,
        "the two spellings of one repeat must denote one sequence at any tract length"
    );
    assert_eq!(
        anchored,
        "CAGT".repeat(5),
        "the whole {}-base tract must be replaced, not the window-clipped prefix",
        core.len()
    );
}

/// The 5' mirror: an anchor *inside* a tract that extends past the window's
/// left edge must still find the run's true start.
///
/// The backward scan stops at `anchor_offset % unit_len`, which is 0 only by
/// coincidence, so a `tract_start == 0` test misses the same way the 3' one did.
#[test]
fn a_tract_extending_past_the_left_edge_is_found_whole() {
    const COPIES: usize = 60;
    let core = "CAGT".repeat(COPIES);
    let sequence = padded(&core);
    let tract_end = 257 + core.len() - 1;
    // Anchor on the last copy, so the run extends 236 bases 5' of it — well
    // past the 128-base initial half-window.
    let anchor = tract_end - 3;

    let anchored = denotes(&sequence, &format!("TEMPLATE:g.{anchor}CAGT[5]"));
    let ranged = denotes(&sequence, &format!("TEMPLATE:g.257_{tract_end}CAGT[5]"));
    assert_eq!(
        anchored, ranged,
        "an anchor anywhere in the run must address the whole run"
    );
}

/// The length of the longest run of `base`.
fn longest_run(sequence: &str, base: char) -> usize {
    sequence
        .split(|c| c != base)
        .map(str::len)
        .max()
        .unwrap_or(0)
}
