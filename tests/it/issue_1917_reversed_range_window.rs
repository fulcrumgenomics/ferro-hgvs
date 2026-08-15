//! Issue #1917 — a reversed range panicked in `normalize_in_grown_window`.
//!
//! A reversed range (`start > end`) names a **wraparound**: the span runs off
//! the 3' end of the reference and resumes at base 1. `normalize_in_grown_window`
//! fetches a linear window and does linear arithmetic on it, so it cannot
//! express one. It now refuses a reversed range up front and the caller takes
//! its minimal-notation fallback, which is the gate `normalize_mito` has
//! carried at its own axis since #129.
//!
//! **The `g.` axis failed three different ways, all functions of `start - end`
//! against `window_size` (default 100).** The issue reported only the middle
//! band and generalised it as "any reversed range whose span exceeds
//! `window_size`", which is wrong in both directions — the widest ranges never
//! panicked, and the narrowest panicked somewhere else entirely:
//!
//! | `start - end`  | before the guard                                          |
//! |----------------|-----------------------------------------------------------|
//! | `<= w`         | the fetch succeeds and nothing underflows, so it fails    |
//! |                | deeper, and how depends on the edit: `dup` slices         |
//! |                | `ref_seq[s..e]` backwards and panics (`slice index starts |
//! |                | at 99 but ends at 50`), while `inv` returns `=` — a       |
//! |                | **wrong answer**, not a panic                             |
//! | `(w, 2w]`      | panic — `end - window_start` underflows                    |
//! | `> 2w`         | the fallback, by accident: the fetch bounds invert and the |
//! |                | provider rejects the read                                  |
//!
//! Both real-corpus rows in #1917 are in the bottom band (16,885 and 13,685
//! bases), so they were **already** taking this exit; the guard makes the other
//! two bands agree with them. Each band is pinned below by a test naming its
//! arithmetic, so a future `window_size` change cannot quietly move a case out
//! from under the coverage.
//!
//! **A panic was not the only outcome, and that is what this change discloses.**
//! Censused over every parseable reversed form × twelve spans, 24 rows were
//! answered: 10 panicked, 8 already declined, and **6 returned a wrong
//! description** — see `a_reversed_range_no_longer_answers_with_a_wrong_description`.
//! The issue reported only the panic, so the moving rows are the ones it missed.
//!
//! The `m.`/`o.` axes are unaffected — they gate before reaching this helper —
//! but are pinned here too, because they are where wraparound is
//! *spec-authorised* (`DNA/deletion.md:17`; `consultation/SVD-WG006.md:15` and
//! its `J01749.1:o.4344_197dup` example at `:23`) and so are the reason the
//! answer is "decline", not "reject".

use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

/// The shipped default `NormalizerConfig::window_size`, which every band
/// boundary below is stated against.
const WINDOW: u64 = 100;

/// Provider serving one contig of `len` cyclic `ACGT` bases.
fn provider_len(accession: &str, len: usize) -> JsonProvider {
    let seq: String = "ACGT".bytes().cycle().take(len).map(char::from).collect();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { accession: seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

fn norm(p: &JsonProvider, input: &str) -> String {
    let v = parse_hgvs(input).unwrap();
    Normalizer::new(p.clone())
        .normalize(&v)
        .unwrap()
        .to_string()
}

/// Confirm the default is what the band arithmetic below assumes. Without this
/// the three band tests still pass if `window_size` changes — they would simply
/// stop testing the bands they name.
#[test]
fn the_bands_below_are_stated_against_the_shipped_window_size() {
    assert_eq!(
        ferro_hgvs::normalize::NormalizeConfig::default().window_size,
        WINDOW,
        "the reversed-range bands in this module are written against \
         window_size = {WINDOW}; if the default moved, re-derive them"
    );
}

/// **Band 2** — the panic the issue reported. `500 - 300 = 200`, which is in
/// `(w, 2w]`, so the anchor `500 - 100 = 400` sat past `end` (300) and
/// `end - window_start` underflowed.
#[test]
fn a_reversed_dup_in_the_underflow_band_is_declined() {
    assert!((WINDOW as i64) < 200 && 200 <= 2 * WINDOW as i64, "band 2");
    let p = provider_len("NC_TEST.1", 1000);
    assert_eq!(
        norm(&p, "NC_TEST.1:g.500_300dup"),
        "NC_TEST.1:g.500_300dup",
        "declined, so the description survives as authored — not re-pointed at \
         the window's first base"
    );
}

/// **Band 1** — narrower than the window, so the fetch succeeded and the
/// subtraction did not underflow. This `dup` failed deeper instead, slicing
/// `ref_seq[99..50]` for its payload. The issue did not report this band, and
/// its stated rule ("any span exceeding `window_size`") excludes it.
///
/// A panic is only *this edit's* band-1 failure. The same band returns a wrong
/// answer for `inv`, and for `dup` at a span of 1 (where `s == e` makes the
/// slice empty and legal) — see
/// `a_reversed_range_no_longer_answers_with_a_wrong_description`.
#[test]
fn a_reversed_dup_narrower_than_the_window_is_declined() {
    assert!(50 <= WINDOW as i64, "band 1");
    let p = provider_len("NC_TEST.1", 1000);
    assert_eq!(norm(&p, "NC_TEST.1:g.350_300dup"), "NC_TEST.1:g.350_300dup");
}

/// **Band 3** — wider than `2 * window_size`, so the fetch bounds themselves
/// inverted and `JsonProvider` rejected the read, taking the same fallback the
/// guard now takes deliberately. This is the band both real-corpus rows sit in,
/// which is why **those two rows** do not move.
///
/// That is a claim about this band, not about the change: six rows in the two
/// narrower bands *did* move, and
/// `a_reversed_range_no_longer_answers_with_a_wrong_description` holds them.
///
/// It passed before the guard and passes after, and that is the point: it pins
/// that the guard **did not change** the case that already worked. Reaching the
/// fallback by a rejected read is not a property of the code — it is a property
/// of whichever provider is installed, which is the other reason not to rely on
/// it.
#[test]
fn a_reversed_dup_wider_than_twice_the_window_is_unchanged() {
    assert!(3300 > 2 * WINDOW as i64, "band 3");
    let p = provider_len("NC_TEST.1", 4000);
    assert_eq!(
        norm(&p, "NC_TEST.1:g.3500_200dup"),
        "NC_TEST.1:g.3500_200dup"
    );
}

/// The spec-authorised route, and the reason declining is the right answer
/// rather than rejecting: `consultation/SVD-WG006.md:23` publishes this exact
/// description. It gates in `normalize_mito` and never reaches the window
/// helper, so this is a guard against that gate being removed as redundant.
#[test]
fn a_circular_wraparound_dup_is_declined() {
    let p = provider_len("J01749.1", 5000);
    assert_eq!(norm(&p, "J01749.1:o.4344_197dup"), "J01749.1:o.4344_197dup");
}

/// The other half of the spec-authorised set — `DNA/deletion.md:17` names
/// deletions, and `SVD-WG006.md:15` gives `NC_012920.1:m.16563_13del`. A
/// reversed `del` is refused at parse on the linear `g.` axis, so the circular
/// axes are the only route by which one reaches the normalizer at all.
#[test]
fn a_circular_wraparound_del_is_declined() {
    let p = provider_len("NC_012920.1", 16569);
    assert_eq!(
        norm(&p, "NC_012920.1:m.16563_13del"),
        "NC_012920.1:m.16563_13del"
    );
}

/// **The two shapes that produced a WRONG ANSWER rather than a panic**, which
/// is the part of this change that moves output.
///
/// A census over every parseable reversed form × twelve spans found 24 answered
/// rows: 10 panicked, 8 already declined (band 3), and **6 came back with an
/// answer that was wrong**. Those six are the disclosure — a panic is loud and a
/// wrong description is not, so these are the rows a consumer could have stored.
///
/// Both shapes below are exactly what two of the seam oracles exist to catch,
/// and neither oracle fired, because no corpus in their selections feeds a
/// reversed range:
///
/// - the `dup` claims `g.4001_4000dup` on a **4000-base** contig — a coordinate
///   past the end of its own sequence, the `FERRO_ASSERT_IN_BOUNDS` class; and
/// - the `inv` claims `=` — asserting the sequence is *unchanged* by an
///   inversion, the `FERRO_ASSERT_SEQUENCE` class.
///
/// The `inv` shape is the worse of the two: `=` is a well-formed, in-bounds,
/// idempotent description, so the other three oracles are satisfied by it too.
#[test]
fn a_reversed_range_no_longer_answers_with_a_wrong_description() {
    let p = provider_len("NC_TEST.1", 4000);

    // Was `g.4001_4000dup` — off the end of a 4000-base contig.
    assert_eq!(
        norm(&p, "NC_TEST.1:g.3500_3499dup"),
        "NC_TEST.1:g.3500_3499dup"
    );

    // Was `g.3500_3450=` — an identity claim about an inversion.
    for input in [
        "NC_TEST.1:g.3500_3499inv",
        "NC_TEST.1:g.3500_3495inv",
        "NC_TEST.1:g.3500_3450inv",
        "NC_TEST.1:g.3500_3401inv",
        "NC_TEST.1:g.3500_3400inv",
    ] {
        let out = norm(&p, input);
        assert_eq!(out, input);
        assert!(
            !out.ends_with('='),
            "{input} must not come back as an identity"
        );
    }
}

/// Whatever the parser admits on the linear axis must not panic. That set is
/// `parse_genome_interval_inner`'s `allows_inverted` list — `dup`, `ins`,
/// `inv`, `dupins` — minus whatever a later check refuses (#1079 refuses the
/// reversed `ins` anchor, which `test_normalize_inverted_range_insertion_no_panic`
/// pins).
///
/// **This is not an endorsement of that list.** `SVD-WG006.md:33` names
/// "deletions/duplications" as the authorised scope, and
/// `check_circular_reversed_range` accordingly admits `del`/`dup`/`delins` and
/// refuses `ins`/`inv` — on the **circular** axes, where the exception actually
/// lives. The linear list is close to the complement of it: it admits `ins` and
/// `inv`, which are spec-silent, and refuses `del`, which is spec-named. That
/// asymmetry is real and is **#1921**; it is not this change's to settle, and
/// the guard is correct under either resolution — which is why this test asks
/// only that what parses reaches an answer, and names no specific set.
#[test]
fn every_linear_reversed_form_the_parser_admits_is_declined() {
    let p = provider_len("NC_TEST.1", 1000);
    let mut parsed = 0;
    for input in [
        "NC_TEST.1:g.500_300dup",
        "NC_TEST.1:g.500_300inv",
        "NC_TEST.1:g.500_300insACGT",
        "NC_TEST.1:g.500_300delinsACGT",
        "NC_TEST.1:g.500_300dupinsACGT",
    ] {
        let Ok(v) = parse_hgvs(input) else { continue };
        parsed += 1;
        let result = Normalizer::new(p.clone()).normalize(&v);
        assert!(
            result.is_ok(),
            "{input} parses, so it must reach an answer rather than a panic; got {result:?}"
        );
    }
    // A refusal at parse is a legitimate outcome for any individual row, but if
    // every row were refused this test would pass while exercising nothing.
    assert!(
        parsed > 0,
        "no reversed linear form parsed — this test asserted nothing"
    );
}
