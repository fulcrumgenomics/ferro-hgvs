//! Issue #91 — Protein 3'-shifting.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/91>.
//!
//! HGVS general.md: "for all descriptions the most 3' position possible
//! is arbitrarily assigned to have been changed." Until this PR, ferro
//! only applied the 3' rule to nucleic-acid axes (c./g./n./r./m.) and
//! left protein variants in their input position. This file pins the
//! contract that protein deletions and duplications now 3'-shift to
//! their canonical (rightmost) anchor when the reference protein
//! provides equal residues immediately downstream.
//!
//! Per the HGVS spec the 3' rule applies to:
//! - **deletions**: a single-residue or contiguous-range deletion shifts
//!   to its 3'-most equivalent position
//! - **duplications**: a duplication shifts to its 3'-most equivalent
//!   position
//! - **substitutions** are not shifted (point change, no rotation)
//! - **frameshifts**, **extensions**, **identity**, and the various
//!   `Repeat` shapes have their own semantics and are not shuffled
//!
//! Insertions (`p.ins`) are handled by issue #92 (the `ins → dup`
//! canonicalization), which depends on this PR landing first so the
//! "most 3' position" is well-defined for the canonical-form rewrite.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// Build a `MockProvider` carrying a single protein sequence keyed by
/// `NP_TESTPROT.1` that contains a run of three alanines (`AAA`) at
/// protein positions 4-6, plus a duplicate `LE` motif at positions 7-8
/// and 9-10. Used to drive the deletion and duplication 3'-shift tests.
///
/// Protein numbering:
/// ```text
///        1234567890
///   seq: MKVAAALELE
///           ^^^ poly-A at p.4..p.6
///              ^^^^ LE at p.7..p.8, repeated at p.9..p.10
/// ```
fn provider_with_polya_protein() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein("NP_TESTPROT.1", "MKVAAALELE");
    provider
}

/// Single-residue deletion in a homopolymer tract must 3'-shift to the
/// rightmost equivalent position. `p.Ala4del` on `MKVAAALELE` shifts
/// from p.Ala4 → p.Ala6 (the last A of the three-residue run).
#[test]
fn protein_single_deletion_shifts_3prime_in_homopolymer() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Ala4del").expect("parse p.Ala4del");

    let normalized = normalizer.normalize(&variant).expect("normalize p.Ala4del");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Ala6del"),
        "expected 3'-shift to land on p.Ala6del; got {out:?}",
    );
}

/// Two-residue deletion `p.Leu7_Glu8del` of `LE` must 3'-shift one
/// residue to `p.Leu9_Glu10del` because the same `LE` motif sits
/// immediately downstream.
#[test]
fn protein_range_deletion_shifts_3prime_through_repeat() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8del").expect("parse p.Leu7_Glu8del");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8del");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu9_Glu10del"),
        "expected 3'-shift to p.Leu9_Glu10del; got {out:?}",
    );
}

/// Duplication in a homopolymer tract follows the same rule. `p.Ala4dup`
/// (duplication of the A at p.4) must 3'-shift to `p.Ala6dup` — the
/// rightmost equivalent anchor in the AAA run.
#[test]
fn protein_duplication_shifts_3prime_in_homopolymer() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Ala4dup").expect("parse p.Ala4dup");

    let normalized = normalizer.normalize(&variant).expect("normalize p.Ala4dup");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Ala6dup"),
        "expected duplication to 3'-shift to p.Ala6dup; got {out:?}",
    );
}

/// Range duplication `p.Leu7_Glu8dup` of the LE motif must shift to
/// `p.Leu9_Glu10dup` because the same motif sits immediately
/// downstream.
#[test]
fn protein_range_duplication_shifts_3prime_through_repeat() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8dup").expect("parse p.Leu7_Glu8dup");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8dup");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Leu9_Glu10dup"),
        "expected duplication to 3'-shift to p.Leu9_Glu10dup; got {out:?}",
    );
}

/// Negative control: a substitution (point change) is never shifted
/// per HGVS spec. `p.Ala4Gly` stays at p.4 regardless of the downstream
/// `AAA` tract.
#[test]
fn protein_substitution_does_not_shift() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Ala4Gly").expect("parse p.Ala4Gly");

    let normalized = normalizer.normalize(&variant).expect("normalize p.Ala4Gly");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Ala4Gly"),
        "substitution must not 3'-shift; got {out:?}",
    );
}

/// Negative control: a deletion outside any tract has nowhere to
/// shift. `p.Val3del` on `MKVAAALELE` — the V at p.3 is unique
/// (positions 1=M, 2=K, 4=A), so the variant stays at p.3.
#[test]
fn protein_deletion_with_no_3prime_equivalent_does_not_shift() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Val3del").expect("parse p.Val3del");

    let normalized = normalizer.normalize(&variant).expect("normalize p.Val3del");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Val3del"),
        "deletion with no 3' equivalent must not shift; got {out:?}",
    );
}

/// Boundary case: a deletion at p.1 (start codon Met) on the protein
/// fixture. `M` is unique at p.1 (next residues are K, V, …), so the
/// variant stays in place. Pins the s0 = 0 branch of the shuffler.
#[test]
fn protein_deletion_at_first_residue_does_not_panic() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Met1del").expect("parse p.Met1del");

    let normalized = normalizer.normalize(&variant).expect("normalize p.Met1del");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Met1del"),
        "deletion at p.1 with no 3' equivalent must stay at p.1; got {out:?}",
    );
}

/// Boundary case: a deletion ending at the last residue of the protein
/// cannot shift because there is no probe residue past the C-terminus.
/// `p.Glu10del` on `MKVAAALELE` — the loop's `probe >= protein.len()`
/// guard fires immediately and the variant stays at p.10.
#[test]
fn protein_deletion_at_last_residue_does_not_shift() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    let variant = parse_hgvs("NP_TESTPROT.1:p.Glu10del").expect("parse p.Glu10del");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Glu10del");
    let out = format!("{}", normalized);

    assert!(
        out.ends_with(":p.Glu10del"),
        "deletion at C-terminus has no probe residue; must not shift. got {out:?}",
    );
}

/// Provider whose protein sequence ends with the stop codon `*` (Ter).
/// The shuffle's rotation predicate `protein[s0] == protein[probe]`
/// must never match a non-Ter residue against `*`, so a deletion that
/// walks up to the Ter naturally halts one step short.
fn provider_with_polya_protein_terminating_in_ter() -> MockProvider {
    let mut provider = MockProvider::new();
    // Position 1234567
    //          MAAAAA*
    provider.add_protein("NP_TER.1", "MAAAAA*");
    provider
}

/// A deletion inside the poly-A run of `MAAAAA*` must 3'-shift up to
/// the last alanine (p.6), one position short of the terminal `*`.
/// Pins that `*` (Ter) acts as a natural shift barrier without
/// requiring a dedicated guard in the shuffle.
#[test]
fn protein_shuffle_halts_before_terminal_ter() {
    let normalizer = Normalizer::new(provider_with_polya_protein_terminating_in_ter());
    let variant = parse_hgvs("NP_TER.1:p.Ala2del").expect("parse p.Ala2del");

    let normalized = normalizer.normalize(&variant).expect("normalize p.Ala2del");
    let out = format!("{}", normalized);

    // Sequence is M A A A A A *; deletion of A starting at p.2 must
    // shift to the last A (p.6) but not into the Ter at p.7.
    assert!(
        out.ends_with(":p.Ala6del"),
        "shuffle must halt at last A (p.6) — Ter at p.7 is a natural barrier; got {out:?}",
    );
    assert!(
        !out.contains("Ter") && !out.contains("p.7"),
        "shuffle must not produce a Ter-positioned anchor; got {out:?}",
    );
}

/// Long-protein regression: a tiny edit on a long synthetic protein
/// must shuffle correctly without depending on a full-protein read.
/// `M` + 2048×`A` + `*` (2050 residues) exercises the windowed
/// shuffler's geometric-growth path — the initial 256-residue
/// look-ahead is exhausted twice before the walk reaches the
/// terminal `*` barrier. Pins that the optimization preserves
/// the input/output contract on long proteins, not just the 10-AA
/// fixture above.
#[test]
fn protein_shuffle_walks_long_homopolymer_with_bounded_reads() {
    let mut provider = MockProvider::new();
    let long_protein: String = std::iter::once('M')
        .chain(std::iter::repeat_n('A', 2048))
        .chain(std::iter::once('*'))
        .collect();
    provider.add_protein("NP_LONG.1", &long_protein);
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NP_LONG.1:p.Ala2del").expect("parse p.Ala2del");

    let normalized = normalizer.normalize(&variant).expect("normalize p.Ala2del");
    let out = format!("{}", normalized);

    // 2048 As at positions 2..=2049; shuffle must shift to the last
    // A (p.2049), one position short of the terminal `*` at p.2050.
    assert!(
        out.ends_with(":p.Ala2049del"),
        "long-homopolymer shuffle must reach the last A; got {out:?}",
    );
}

/// Regression: when a deletion that carried an explicit, NON-redundant
/// residue list is 3'-shifted, the input residue list is stale because
/// the new window covers a rotated stretch of the protein. The output
/// must NOT preserve the pre-shift letters or it emits a
/// self-contradictory HGVS string (positions saying one thing, `delXYZ`
/// saying another). Per HGVS spec the residue list is optional, so the
/// canonical fix is to drop it when the shift moved the anchor —
/// Display then falls back to the position-only `del` form, which is
/// unambiguous against the shifted interval.
///
/// Redundant residue lists (where the first/last letters match the
/// position AAs) are already stripped pre-shift by
/// `is_redundant_protein_deletion_sequence`. This test pins the
/// remaining case: a sequence whose interior letters differ from the
/// position AAs, which would otherwise survive the shift unchanged.
///
/// Protein `MKVAAALELE`: input `p.Leu7_Glu8delLeuAla` declares window
/// positions Leu7,Glu8 (which validate against the reference) but
/// claims the deleted body is `LeuAla` — the trailing `Ala` does not
/// match end_pos.aa (Glu), so the redundancy stripper leaves the
/// sequence in place. The 3'-shifter then moves the anchor to
/// `p.Leu9_Glu10` (the duplicate LE motif), at which point `delLeuAla`
/// is doubly wrong (window no longer Leu-Ala, and never was).
#[test]
fn protein_range_deletion_with_stale_explicit_residues_drops_them_after_shift() {
    let normalizer = Normalizer::new(provider_with_polya_protein());
    // The parser accepts a `delSEQ` body whose interior letters do not
    // match the position AAs — the body is informational, not
    // structural. `p.Leu7_Glu8delLeuAla` validates Leu7,Glu8 against
    // the reference but the residue list `[Leu,Ala]` is non-redundant
    // (last letter Ala != end_pos Glu), so the pre-shift redundancy
    // stripper leaves it in place.
    let variant =
        parse_hgvs("NP_TESTPROT.1:p.Leu7_Glu8delLeuAla").expect("parse p.Leu7_Glu8delLeuAla");

    let normalized = normalizer
        .normalize(&variant)
        .expect("normalize p.Leu7_Glu8delLeuAla");
    let out = format!("{}", normalized);

    // The shift moves the anchor to p.Leu9_Glu10. The pre-shift
    // residue list `LeuAla` no longer describes the deleted window, so
    // it must be dropped — emit position-only `del`.
    assert!(
        out.ends_with(":p.Leu9_Glu10del"),
        "expected 3'-shift with stale residues dropped to p.Leu9_Glu10del; got {out:?}",
    );
    assert!(
        !out.contains("delLeuAla"),
        "stale pre-shift residue list must not survive the 3'-shift; got {out:?}",
    );
}
