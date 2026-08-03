//! Issue #1270 — the del→repeat path answered the codon-frame gate on the
//! *input* span where ins→repeat answers it on the tract the repeat occupies.
//!
//! `repeated.md` restricts repeat notation in a coding context to units that
//! are a multiple of three. #1185/#1210 established that the restriction is
//! about where the repeat **lands**, not where the input edit sat, and wired
//! `insertion_to_repeat` to re-ask for the tract. `deletion_to_repeat` was
//! still passed `gate.input_span_is_coding()` alone, so the two paths could
//! answer differently about one and the same tract.
//!
//! ## The issue filed this as unverified; it is reachable
//!
//! "I have not constructed a case that reaches it, so this is an asymmetry with
//! a plausible consequence rather than a demonstrated defect." It is
//! demonstrated. On a `T` run spanning `cds_end`, measured before the fix:
//!
//! ```text
//! NM_TEST.1:c.5_6insTT  ->  NM_TEST.1:c.4_*5T[16]   tract-aware: repeat emitted
//! NM_TEST.1:c.5_6del    ->  NM_TEST.1:c.*4_*5del    input-span:  repeat suppressed
//! ```
//!
//! Same tract, same unit (`T`, length 1, not a multiple of three), opposite
//! answers. The deletion 3'-shifts clean out of the CDS — the `c.`-axis
//! shuffle does not clamp it at `cds_end` here — so by the time the repeat is
//! spelled the tract is not coding, and the carve-out applies to it exactly as
//! it does to the insertion.
//!
//! The fix is the one the issue suggests: give `deletion_to_repeat` the same
//! `gate.cds_bounds()` treatment, and move its gate check to after the tract is
//! extended so there is a tract to ask about.
//!
//! Reference-free (`MockProvider` via `SyntheticBuilder`), so these run on PR
//! CI rather than only in the manifest-backed nightly.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// ```text
/// tx      1  2  3 | 4  5  6 | 7 ................ 20 | 21 ....... 26
///         G  C  A | A  T  G | T x 14                | G C G C G C
/// region  5'UTR   | CDS 4..15 .....| 3'UTR 16..26 ..|
/// ```
///
/// The `T` run spans tx 7..20 and so crosses `cds_end = 15`: its 5' half is
/// coding and its 3' half is not. That is the whole point of the fixture — a
/// tract whose own verdict differs from that of an input span inside it.
const CORE: &str = "GCAATGTTTTTTTTTTTTTTGCGCGC";
const CDS: (u64, u64) = (4, 15);

/// A tract of the same shape with a **3-base** unit, where the gate does not
/// block either way. Used to show the change is scoped to the disagreement.
const CORE_UNIT3: &str = "GCAATGCAGCAGCAGCAGCAGCAGGCGC";

fn normalize(core: &str, cds: (u64, u64), input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    let provider: MockProvider = SyntheticBuilder::cds(core, cds.0, cds.1, Strand::Plus).build();
    Normalizer::new(provider)
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string()
}

/// The defect proper: a CDS-resident deletion inside a tract that crosses
/// `cds_end` must reach repeat notation, because the tract it lands on is not
/// wholly coding.
///
/// The reference run is 14 `T`, so deleting `n` of them leaves `14 - n`.
#[test]
fn a_deletion_whose_tract_leaves_the_cds_reaches_repeat_notation() {
    for (input, expected) in [
        ("NM_TEST.1:c.5_6del", "NM_TEST.1:c.4_*5T[12]"),
        ("NM_TEST.1:c.4_5del", "NM_TEST.1:c.4_*5T[12]"),
        ("NM_TEST.1:c.10_11del", "NM_TEST.1:c.4_*5T[12]"),
        ("NM_TEST.1:c.5_7del", "NM_TEST.1:c.4_*5T[11]"),
        ("NM_TEST.1:c.4_9del", "NM_TEST.1:c.4_*5T[8]"),
    ] {
        assert_eq!(
            normalize(CORE, CDS, input),
            expected,
            "`{input}` must be spelled as a repeat over the tract it lands on"
        );
    }
}

/// The asymmetry the issue is about, asserted directly: the two paths must
/// agree about one tract.
///
/// Stated as a relation between the two rather than as two independent pinned
/// strings, so it keeps meaning if the canonical spelling ever moves — the
/// contract is that del and ins answer the *same* question, and the counts
/// bracket the reference's 14 copies from either side.
#[test]
fn the_deletion_and_insertion_paths_agree_about_one_tract() {
    let deleted = normalize(CORE, CDS, "NM_TEST.1:c.5_6del");
    let inserted = normalize(CORE, CDS, "NM_TEST.1:c.5_6insTT");

    // Split at the *location/edit* boundary, not on the unit letter: the
    // accession `NM_TEST.1` contains a `T`, so `split('T')` returns `"NM_"` for
    // both sides and the comparison is vacuously true. (It was written that way
    // first; the self-review caught it.)
    let tract_of = |s: &str| {
        s.rsplit_once("T[")
            .map(|(location, _)| location.to_string())
            .unwrap_or_else(|| panic!("expected a `T[...]` repeat spelling; got `{s}`"))
    };
    assert_eq!(
        tract_of(&deleted),
        tract_of(&inserted),
        "both must name the same tract; got `{deleted}` and `{inserted}`"
    );
    assert_eq!(
        tract_of(&deleted),
        "NM_TEST.1:c.4_*5",
        "and it must be the tract that crosses cds_end"
    );
    assert_eq!(
        deleted, "NM_TEST.1:c.4_*5T[12]",
        "14 reference copies less 2"
    );
    assert_eq!(
        inserted, "NM_TEST.1:c.4_*5T[16]",
        "14 reference copies plus 2"
    );
}

/// A tract that stays wholly inside the CDS must still be gated. Without this
/// the change would read as "the codon rule was removed" rather than "it is
/// answered about the right span".
///
/// Same reference and same input as the first test — only `cds_end` moves, from
/// 15 to 24, which brings the whole `T` run inside the CDS. That isolates the
/// variable to the tract's own verdict: nothing about the edit or the unit
/// changes, and the repeat must disappear.
#[test]
fn a_tract_wholly_inside_the_cds_is_still_gated() {
    // CDS 4..24 puts the entire `T` run (tx 7..20) inside the CDS, so the
    // tract's own verdict is coding and the 1-base unit must stay blocked.
    let output = normalize(CORE, (4, 24), "NM_TEST.1:c.5_6del");
    assert!(
        !output.contains("T["),
        "a coding tract with a 1-base unit must not be spelled as a repeat; got `{output}`"
    );
    // Pinned as well as negated: "contains no repeat" alone would also be
    // satisfied by the normalizer erroring into some unrelated shape.
    assert_eq!(
        output, "NM_TEST.1:c.16_17del",
        "it must fall back to the plain 3'-shifted deletion"
    );
}

/// A unit that satisfies the codon rule is unaffected in either direction —
/// the gate never blocked it, so the re-ask must not change the answer.
#[test]
fn a_three_base_unit_is_unchanged() {
    for input in ["NM_TEST.1:c.4_9del", "NM_TEST.1:c.7_12del"] {
        assert_eq!(
            normalize(CORE_UNIT3, CDS, input),
            "NM_TEST.1:c.4_*9CAG[4]",
            "`{input}` must keep the repeat spelling it already had"
        );
    }
}

/// The tract-aware verdict must not depend on strand.
///
/// `cds_span` is in transcript coordinates and the provider serves
/// transcript-oriented bases, so a minus-strand record should reach the same
/// answer — but "wrong shift direction or strand handling" is exactly what this
/// arithmetic could get wrong, and every other test here is plus-strand.
/// Asserted as parity so it survives a change in the canonical spelling.
#[test]
fn the_tract_verdict_is_strand_independent() {
    for input in ["NM_TEST.1:c.5_6del", "NM_TEST.1:c.4_9del"] {
        let variant = parse_hgvs(input).expect("input must parse");
        let minus: MockProvider = SyntheticBuilder::cds(CORE, CDS.0, CDS.1, Strand::Minus).build();
        let on_minus = Normalizer::new(minus)
            .normalize(&variant)
            .expect("lenient normalization must not reject")
            .to_string();
        assert_eq!(
            normalize(CORE, CDS, input),
            on_minus,
            "`{input}` must normalize the same on either strand"
        );
    }
}

/// A deletion with no tandem tract around it must stay a plain deletion.
#[test]
fn a_deletion_outside_any_tract_is_unchanged() {
    assert_eq!(
        normalize(CORE, CDS, "NM_TEST.1:c.1_2del"),
        "NM_TEST.1:c.1_2del"
    );
}

/// Normalizing the repeat form again must return it byte for byte.
///
/// This is the property #1210 was actually filed about — the input-span verdict
/// made pass 1 and pass 2 disagree, because pass 2's input was by then
/// UTR-resident — so it is asserted here for the deletion path too rather than
/// left to the global oracle.
#[test]
fn the_repeat_spelling_is_a_fixed_point() {
    let once = normalize(CORE, CDS, "NM_TEST.1:c.5_6del");
    assert_eq!(normalize(CORE, CDS, &once), once, "`{once}` must settle");
}
