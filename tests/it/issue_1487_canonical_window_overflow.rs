//! Issue #1487: a parseable description carrying a coordinate near `i64::MAX`
//! panicked with `attempt to add with overflow` instead of being refused.
//!
//! # The bug
//!
//! `canonicalize_from_sequence` opens its derivation window by padding the
//! union of every member's footprint by `CANONICAL_PAD`:
//!
//! ```text
//! let (c_lo, c_hi) = edit_span_union(&edits)?;
//! let mut w_lo = (c_lo - CANONICAL_PAD).max(1);
//! let mut w_hi = c_hi + CANONICAL_PAD;      // <- overflows
//! ```
//!
//! Nothing bounds `c_hi` beforehand — the width test that would refuse an
//! over-wide window runs on the *padded* values, one line too late — so a
//! member at `c.9223372036854775807` wraps the addition.
//!
//! # Why it is reached through the ordinary public API
//!
//! The description parses. `parse_hgvs` holds no provider and imposes no
//! magnitude limit on a coordinate, so the panic is reachable from
//! `Normalizer::normalize` with nothing unusual in the caller. Through the
//! Python bindings a Rust panic surfaces as `pyo3_runtime.PanicException`,
//! which subclasses `BaseException` and so escapes `except Exception` — the
//! same contract break `issue_1244_equivalence_overlap_panic` records.
//!
//! # Why the release build is worse than the panic
//!
//! Debug builds panic. Release builds wrap instead: `w_hi` becomes negative,
//! the window is nonsense, and the derivation proceeds against it silently.
//! A crash is at least observable.
//!
//! # The three sites, and why no one of them covers the others
//!
//! The padding is one of **three** unchecked additions on this path, not one of
//! two. The second converts the padded window onto the served sequence's axis:
//!
//! ```text
//! let start0 = w_lo + frame.delta - 1;
//! let end0 = w_hi + frame.delta;
//! ```
//!
//! The width test between them refuses a *wide* window, so it looks like it
//! guards this. It does not: a span that sits near `i64::MAX` while remaining
//! narrower than `MAX_CANONICAL_WINDOW` passes the width test and then wraps
//! against a positive `frame.delta`.
//!
//! The third is `enclosing_exon`'s own `hi + frame.delta`, reached from the
//! exon clamp that sits *between* the two above. It runs on the **raw** member
//! union rather than on the padded window, so the pad's saturation cannot reach
//! it, and it is the sibling of `crosses_exon_junction` — which had already
//! been hardened with `checked_add` while this one was missed.
//!
//! No test overflows at all three by accident, so the attribution is recorded
//! here rather than left to be re-derived. It was **measured**, by instrumenting
//! the three refusals and running this module, not inferred from the arithmetic.
//!
//! Read the column as *"whose addition overflows there"*, which is not the same
//! question as *"who executes that line"* — only the first is coverage, and
//! conflating the two is how the third row came to be written too narrowly:
//!
//! | site | whose addition overflows there |
//! |---|---|
//! | the pad, `c_hi + CANONICAL_PAD` | every test below that carries an extreme coordinate. It saturates rather than refusing, and for `an_extreme_coordinate_is_refused_instead_of_panicking` and `the_genomic_axis_is_covered_too` it is also the only one of the three sites reached at all, both being refused by the width test immediately after it |
//! | the axis conversion, `w_hi + frame.delta` | `a_narrow_span_at_an_extreme_position_is_refused_instead_of_panicking`, and the two exon-clamp tests below it |
//! | `enclosing_exon`'s `hi + frame.delta` | `two_adjacent_members_at_the_top_of_the_range_are_refused_instead_of_panicking` and `a_lone_member_at_an_extreme_coordinate_is_refused_instead_of_panicking`, and nothing else in the suite. The narrow-span test *executes* both of that site's `checked_add`s and overflows at neither |
//!
//! Deliberately no `merge.rs` line numbers: they drift on the next insertion,
//! and each site is named by an expression that does not.
//!
//! # The fix
//!
//! Refuse, do not wrap. The pad saturates — `w_lo` is clamped to `1`
//! immediately and a saturated `w_hi` can only make the width test larger, so
//! saturation lands on the refusal that test already gives an over-wide
//! window — and both conversions are checked, taking the same conservative
//! exit the neighbouring `start0 < 0` guard already takes.
//!
//! Declining is the correct outcome, not a workaround: refusing the
//! sequence-first derivation hands the description back to the per-member
//! pipeline, which is what decides whether such a coordinate is an error.
//! `crosses_exon_junction` already reasons this way — it takes the same exit
//! for an unservable transcript — and it is *that* function, not
//! `enclosing_exon`, whose comment named this very panic as still live upstream
//! and "tracked separately". This change is that tracking, so the comment is
//! rewritten with it; `enclosing_exon` is the sibling it did not mention.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::JsonProvider;
use ferro_hgvs::Normalizer;

/// The issue's reproducer: a two-member cis allele whose first member sits at
/// the top of the `i64` range, so `c_hi + CANONICAL_PAD` overflows.
const EXTREME_ALLELE: &str = "NM_001234.1:c.[9223372036854775806_9223372036854775807del;10del]";

/// Two **adjacent** members at the top of the range — `i64::MAX` and the base
/// below it — so the raw member union is two wide while sitting as high as a
/// coordinate can go.
///
/// The issue calls this "the single-member spelling" and reports that it
/// "normalizes fine". Both halves are wrong: the string is a two-member bracket
/// allele, and it panicked identically. Recorded under a name that matches its
/// value, because the mislabel is what made the panic look specific to the
/// span-union path — the discriminator is how *high* the union sits, not
/// whether there is more than one member.
///
/// It is also what keeps the third site covered. A two-wide union survives the
/// width test, and `crosses_exon_junction` — hardened first, so already
/// checked — answers "no junction crossed" on its own overflow rather than
/// refusing, so control reaches `enclosing_exon` and its addition.
const EXTREME_ADJACENT_MEMBERS: &str =
    "NM_001234.1:c.[9223372036854775807del;9223372036854775806del]";

/// The genuine single-member spelling the issue meant to name. Added rather
/// than substituted for the constant above, which owns the third site's
/// coverage and cannot be swapped out without losing it.
///
/// A lone member reaches the same pass — `sequence_first_pass` admits any
/// single member whose edit has an inner `NaEdit` — so this is not a second
/// spelling of an allele test: it pins that the refusal is not conditional on
/// bracket syntax, and it reaches the exon clamp too.
const EXTREME_LONE_MEMBER: &str = "NM_001234.1:c.9223372036854775807del";

/// `NM_001234.1` from `JsonProvider::with_test_data` — the provider the issue's
/// own reproducer uses, and the one every other coding-axis derivation test on
/// this path is written against.
fn normalizer() -> Normalizer<JsonProvider> {
    Normalizer::new(JsonProvider::with_test_data())
}

// ---------------------------------------------------------------------------
// The reported crash.
// ---------------------------------------------------------------------------

#[test]
fn an_extreme_coordinate_is_refused_instead_of_panicking() {
    let parsed = parse_hgvs(EXTREME_ALLELE).expect("the description parses");

    // Whether the normalizer accepts or refuses is its own choice; what it may
    // not do is panic. `normalize` returns `Result`, so a refusal is a verdict
    // the caller can handle and a panic is not.
    let outcome = normalizer().normalize(&parsed);

    // Deliberately not asserting `is_err()`. The coordinate is past the end of
    // any real sequence, so a refusal is the expected answer today — but if a
    // later change makes the per-member pipeline answer it some other way,
    // that is a behaviour question, not a regression of this fix. The contract
    // this test owns is only that control returns.
    let _ = outcome;
}

// ---------------------------------------------------------------------------
// The second site: narrow span, extreme position.
//
// This is the case the width test does not catch. Both members sit within a
// few bases of each other — so the padded window is ~260 wide, far under
// `MAX_CANONICAL_WINDOW` — but near `i64::MAX`, so the axis conversion
// `w_hi + frame.delta` is what wraps. A fix that only saturates the pad leaves
// this one panicking.
// ---------------------------------------------------------------------------

#[test]
fn a_narrow_span_at_an_extreme_position_is_refused_instead_of_panicking() {
    let parsed = parse_hgvs("NM_001234.1:c.[9223372036854775700del;9223372036854775702del]")
        .expect("the description parses");
    let _ = normalizer().normalize(&parsed);
}

#[test]
fn the_genomic_axis_is_covered_too() {
    // Site 1 only — measured, not assumed, and worth saying because the
    // obvious reading is that a different `frame.delta` reaches the other two
    // additions with a different addend. It does not, three times over:
    // `axis_frame` hard-codes `delta: 0` for `CisKind::Genome`, so
    // `w_hi + frame.delta` adds zero and cannot wrap; `enclosing_exon` returns
    // `None` for that kind before reaching its own addition; and this span runs
    // `10..i64::MAX`, so the width test refuses it before either is reached.
    //
    // What it does pin is that the pad saturation sits on the **shared** path
    // and is not a coding-branch special case: a genomic cis allele carrying an
    // extreme coordinate is refused rather than wrapped. It would not catch a
    // fix applied only to the coding branch of the other two sites — the two
    // exon-clamp tests above are what cover those.
    let normalizer = Normalizer::new(SyntheticBuilder::genomic("ACGTACGTACGTACGTACGT").build());
    let parsed = parse_hgvs("NC_TEST.1:g.[9223372036854775806_9223372036854775807del;10del]")
        .expect("the description parses");
    let _ = normalizer.normalize(&parsed);
}

// ---------------------------------------------------------------------------
// The third site: the exon clamp.
//
// Both inputs below carry a raw member union at the very top of the `i64` range
// that is at most two bases wide, so the padded window passes the width test
// and control reaches the exon clamp. `enclosing_exon`'s `hi + frame.delta`
// runs on that RAW union rather than on the padded window, so the pad's
// saturation does not stand in for it: `NM_001234.1` has `cds_start = 5`, hence
// `frame.delta = 4`, and `i64::MAX + 4` wraps.
//
// Reaching the site is not the same as overflowing at it, and only the second
// is coverage. `a_narrow_span_at_an_extreme_position_…` also *executes* both of
// `enclosing_exon`'s `checked_add`s: its raw union sits at `…700`/`…702`, so
// `crosses_exon_junction` adds `frame.delta` without wrapping, answers "no
// junction crossed" on real arithmetic, and control walks straight into the
// clamp — where `+4` fits too. It exercises the addition without exercising the
// guard. The reproducer above and the genomic case do not get even that far,
// both being refused by the width test first.
//
// So these two inputs are the only ones in the suite whose addition at that
// site OVERFLOWS, and deleting either of them silently uncovers it.
// ---------------------------------------------------------------------------

#[test]
fn two_adjacent_members_at_the_top_of_the_range_are_refused_instead_of_panicking() {
    let parsed = parse_hgvs(EXTREME_ADJACENT_MEMBERS).expect("the description parses");
    let _ = normalizer().normalize(&parsed);
}

#[test]
fn a_lone_member_at_an_extreme_coordinate_is_refused_instead_of_panicking() {
    let parsed = parse_hgvs(EXTREME_LONE_MEMBER).expect("the description parses");
    let _ = normalizer().normalize(&parsed);
}

// ---------------------------------------------------------------------------
// The refusal must be narrow: ordinary coordinates still derive.
//
// Saturating arithmetic is invisible for every input that does not overflow,
// but a fix written as an eager magnitude check would not be — so pin that a
// normal cis allele on the same reference still normalizes.
// ---------------------------------------------------------------------------
#[test]
fn an_ordinary_cis_allele_on_the_same_reference_still_normalizes() {
    let parsed = parse_hgvs("NM_001234.1:c.[4del;8del]").expect("the description parses");
    normalizer()
        .normalize(&parsed)
        .expect("an ordinary two-member allele must still normalize");
}

// ---------------------------------------------------------------------------
// Why the two remaining unchecked additions upstream of the pad are safe.
//
// `crosses_exon_junction`'s comment used to end "nothing on this path adds
// unchecked", which is false as written: `collect_canonical_edits` tests an
// insertion anchor with `if e != s + 1`, and `edit_span_union` reads an `Ins`
// footprint as `(*gap, *gap + 1)`. Both run before the pad and neither is
// checked.
//
// The panic class is nonetheless closed, and the thing that closes it is the
// PARSER, not arithmetic — so it is pinned here rather than left as prose that
// a later parser change could silently invalidate. Each refusal below removes
// one way for an `Insertion` to arrive with `s == i64::MAX`; together they
// leave `s <= i64::MAX - 1`, which is what makes `s + 1` in range at both
// sites. The final case is the positive control: the largest anchor that DOES
// parse, so the bound is exercised rather than merely asserted.
// ---------------------------------------------------------------------------

#[test]
fn the_parser_is_what_keeps_an_insertion_anchor_off_i64_max() {
    // A coordinate above `i64::MAX` is not a number the grammar accepts, so
    // `s` can never exceed it in the first place.
    parse_hgvs("NM_001234.1:c.9223372036854775807_9223372036854775808insA")
        .expect_err("a coordinate past i64::MAX must not parse");

    // With the ceiling in place, `s == i64::MAX` forces `e <= s`. `e == s` is a
    // single-position insertion, refused by `DNA/insertion.md:95-101`…
    parse_hgvs("NM_001234.1:c.9223372036854775807_9223372036854775807insA")
        .expect_err("a single-position insertion anchor must be refused");

    // …and `e < s` is a reversed range, refused before that. Checked on the
    // coding axis and again on a genomic one, since the two take different
    // parser paths and only one of them needs to hold for the argument to leak.
    parse_hgvs("NM_001234.1:c.9223372036854775807_10insA")
        .expect_err("a reversed insertion anchor must be refused");
    parse_hgvs("NC_TEST.1:g.9223372036854775807_9223372036854775807insA")
        .expect_err("a single-position insertion anchor must be refused on g. too");

    // The positive control, and the bound itself: `i64::MAX - 1` and `i64::MAX`
    // are adjacent, so this is the highest anchor HGVS admits. It parses, which
    // is what makes the two additions' worst case `(i64::MAX - 1) + 1`.
    //
    // Deliberately NOT normalized. Every test above that hands an extreme
    // coordinate to `normalize` panics under `FERRO_ASSERT_SEQUENCE=1` — see
    // the note on the fourth site below — and there is no reason to add a fifth
    // for a parser claim that `parse_hgvs` alone settles.
    parse_hgvs("NM_001234.1:c.9223372036854775806_9223372036854775807insA")
        .expect("the highest legal insertion anchor must still parse");
}

// ---------------------------------------------------------------------------
// A FOURTH site, on a different entry point, and NOT closed by this change.
//
// `convert::mapper::cds_to_tx_exon_aware` returns `start_tx + offset` unchecked
// on both of its early-return branches, and an extreme `c.` coordinate reaches
// it through `hgvs_to_spdi` rather than through normalization — so nothing on
// the derivation path this module is about stands in front of it.
//
// Measured, not inferred. Of the four seam-oracle flags only
// `FERRO_ASSERT_SEQUENCE=1` converts the output, and under it exactly the four
// tests above that carry an extreme coordinate into `normalize` panic at
// `src/convert/mapper.rs` with `attempt to add with overflow`; the other three
// flags are green. CI's `test-oracle` job does not set that flag (`sweeps` and
// the nightly do, over selections that exclude this module), which is why the
// gate is green and why the site would otherwise stay invisible.
//
// Left to a separate change on purpose: a different module, a different entry
// point, and a refusal there is a behaviour question for the conversion API
// rather than for this derivation. Recorded here with its reproducer so the
// next person does not re-derive it — which is the failure the `crosses_exon_
// junction` comment's own "tracked separately" caused the first time.
// ---------------------------------------------------------------------------
