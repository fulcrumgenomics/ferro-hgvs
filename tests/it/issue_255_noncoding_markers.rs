//! Audit for #81 § J2 — round-trip fidelity for the out-of-zone `-N` and `*N`
//! markers through the bare [`parse_hgvs`] entry.
//!
//! # The premise this file was written on is FALSE (#1748)
//!
//! This module used to open *"Per HGVS v21, the `n.` coordinate system uses
//! `-N` upstream of transcript start … `*N` downstream of transcript end"*.
//! That is not in v21. It is the **coding** axis's zone list transplanted onto
//! `n.`, and `background/numbering.md:50`–`:54` says the opposite:
//!
//! - `:52` enumerates the whole axis — "nucleotide numbering is `n.1`, `n.2`,
//!   `n.3`, ..., etc., from the first to the last nucleotide of the reference
//!   sequence". No `*` zone and no negative one.
//! - `:53` grants introns, and only introns. Its subject is "nucleotides **in
//!   introns**", so its "see above" is scoped **by its own sentence** to the
//!   coding section's **introns** bullet (`:33`), not to the **UTR** bullet
//!   (`:28`) where `c.-1`/`c.*1` are defined. The scoping is carried by the
//!   sentence rather than by the link: `:53`'s anchor `#DNAc` targets the whole
//!   coding section. `:38`/`:39` do spell `c.-85+1` and `c.*37+1` inside the
//!   introns bullet, but both are anchored to the CDS — on an axis with no CDS
//!   they have nothing to number from.
//! - `:54` forbids describing "variants in nucleotides beyond the boundaries of
//!   a transcript reference sequence, using that transcript reference
//!   sequence". The coding section's closing bullet `:43`–`:44` is that
//!   sentence verbatim, so `c.*1` and the prohibition coexist on `c.` only
//!   because the UTR is *inside* the transcript. `n.` has no UTR bullet, so
//!   `:54` stands unopposed.
//!
//! `numbering.md:45` records that the sibling proposal — extending the
//! numbering "to specifically mark non-transcribed nucleotides" — "[has] been
//! made but [was] rejected".
//!
//! # WHAT MOVED IN THIS FILE, AND WHAT DELIBERATELY DID NOT
//!
//! The two markers are read the same way and are enforced at **different
//! stages**, on measured cost. This file splits along that line:
//!
//! | section | axis now | why |
//! |---|---|---|
//! | **`*N`** (§2, and the two downstream crossings in §5) | **re-pointed to `c.`** | `n.*N` is refused at parse in **every** mode as `E1003`, so it cannot round-trip through any door. `c.*N` is genuinely legal — the zone is anchored to the stop codon and is inside the transcript — so the assertions keep their subject |
//! | **`-N`** (§1, and the two upstream crossings in §5) | **stays on `n.`** | `n.-N` is refused in **strict only**. The bare [`parse_hgvs`] entry applies no mode at all (#1632), so these still round-trip — and they are now the live guard on five real ClinVar rows |
//!
//! **Do not "finish the job" by re-pointing §1 at `c.` as well.** Those tests
//! are load-bearing in a way §2's were not. Five rows across ferro's committed
//! corpora are `n.-N` descriptions NCBI publishes today —
//! `NR_003051.3:n.-57T>C` and `NR_003051.3:n.-30_-7dup` (RMRP),
//! `LRG_163t1:n.-5delins17` (RMRP), `NR_029595.1:n.-4771G>T` (MIR208A) and
//! `NR_033294.1:n.-6G>A` (SNORD118) — and RMRP promoter variants are the
//! clinically conventional case for the spelling. `n.*N`, by contrast, appears
//! **0** times in the same 103,762 `n.`-axis rows, which is what made refusing
//! it free. The five rows are pinned by name in
//! `issue_1748_noncoding_axis_zones.rs`; these round-trips are the shape guard
//! behind them.
//!
//! # Why these tests stay, and stay green
//!
//! What they pin is **round-trip fidelity, not conformance** — that nothing
//! silently drops a `*`, mis-signs a `-`, or reorders a range endpoint on the
//! way back out. That property is worth pinning on whichever axis the spelling
//! is admitted on, which is why §2 moved rather than being deleted.
//!
//! These tests call the bare [`parse_hgvs`], which applies no `ErrorConfig` at
//! all (#1632) and so sits on none of the three mode arms. That is also the
//! entry every `ferro` CLI subcommand parses through, which is why it is worth
//! pinning in its own right. The strict refusal, the unconditional refusal and
//! the mode schedule are pinned separately, in
//! `issue_1748_noncoding_axis_zones.rs`.
//!
//! **Do not read the `n.-N` round-trips as an endorsement of the spelling.**
//! They assert that ferro reproduces what it was given; they assert nothing
//! about whether what it was given is conformant HGVS. It is not.

use ferro_hgvs::parse_hgvs;

/// The non-coding transcript #255's audit used, kept for the `-N` sections.
const TX: &str = "NR_037639.1";

/// A **coding** transcript for the `*N` sections, which moved off the `n.` axis.
/// CFTR: 6070 nt, CDS `[70, 4513)` 0-based, so `c.4443` is its last CDS base and
/// `c.*1557` its last 3'UTR base — the boundary crossings below are real
/// coordinates rather than synthetic ones.
const CDS_TX: &str = "NM_000492.4";

#[track_caller]
fn assert_round_trips(s: &str) {
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("parse {s:?} failed: {e}"));
    let out = format!("{}", v);
    assert_eq!(out, s, "round-trip mismatch for {s:?}");
}

// =============================================================================
// SECTION 1 — Upstream `-N` marker, still on the `n.` axis
//
// Refused in strict only, so the bare entry still round-trips these. This is
// the half that guards the five real ClinVar rows — see the module docs before
// moving any of it.
// =============================================================================

mod upstream {
    use super::*;

    #[test]
    fn neg_single_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.-1A>G"));
    }

    #[test]
    fn neg_deeper_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.-5A>G"));
    }

    #[test]
    fn neg_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.-100_-50del"));
    }

    #[test]
    fn neg_short_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_-1del"));
    }

    #[test]
    fn neg_range_delins_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_-1delinsACG"));
    }

    #[test]
    fn neg_range_dup_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_-1dup"));
    }

    #[test]
    fn neg_range_inv_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_-1inv"));
    }

    /// The shapes the five ClinVar rows actually take, spelled on the audit's
    /// own transcript so this section keeps testing the marker rather than the
    /// accession. The rows themselves are pinned by name in
    /// `issue_1748_noncoding_axis_zones.rs`.
    #[test]
    fn the_real_world_upstream_shapes_round_trip() {
        assert_round_trips(&format!("{TX}:n.-57T>C"));
        assert_round_trips(&format!("{TX}:n.-30_-7dup"));
        assert_round_trips(&format!("{TX}:n.-6G>A"));
    }
}

// =============================================================================
// SECTION 2 — Downstream `*N` marker, RE-POINTED to the coding axis
//
// `n.*N` no longer parses on any door (#1748), so these assertions moved to the
// axis where `*N` is genuinely legal. What each pins is unchanged: the `*` must
// survive the round trip, and range endpoints must not be reordered.
// =============================================================================

mod downstream {
    use super::*;

    #[test]
    fn star_single_sub_round_trips() {
        assert_round_trips(&format!("{CDS_TX}:c.*1A>G"));
    }

    #[test]
    fn star_deeper_sub_round_trips() {
        assert_round_trips(&format!("{CDS_TX}:c.*100A>G"));
    }

    #[test]
    fn star_range_del_round_trips() {
        assert_round_trips(&format!("{CDS_TX}:c.*1_*5del"));
    }

    #[test]
    fn star_range_delins_round_trips() {
        assert_round_trips(&format!("{CDS_TX}:c.*3_*5delinsACG"));
    }

    #[test]
    fn star_range_dup_round_trips() {
        assert_round_trips(&format!("{CDS_TX}:c.*3_*5dup"));
    }

    #[test]
    fn star_range_inv_round_trips() {
        assert_round_trips(&format!("{CDS_TX}:c.*3_*5inv"));
    }

    /// The refusal that replaced the `n.` half of this section, asserted here so
    /// the move is visible from the file it moved out of.
    #[test]
    fn the_noncoding_axis_no_longer_admits_the_star() {
        assert!(
            parse_hgvs(&format!("{TX}:n.*1A>G")).is_err(),
            "`n.*1` is refused at parse in every mode (#1748)"
        );
    }
}

// =============================================================================
// SECTION 3 — Interior positions (positive controls)
// =============================================================================

mod interior {
    use super::*;

    #[test]
    fn interior_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.100A>G"));
    }

    #[test]
    fn interior_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.1_100del"));
    }
}

// =============================================================================
// SECTION 4 — Intronic offsets `N+M` / `N-M`
//
// `:53` grants these explicitly, so they are untouched by #1748 and stay on the
// `n.` axis. The negative-offset cases matter most: the `-1` in `n.100-1` is a
// distance from position 100, not a zone, and must never be read as one.
// =============================================================================

mod intronic_offsets {
    use super::*;

    #[test]
    fn plus_offset_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.100+1A>G"));
    }

    #[test]
    fn minus_offset_sub_round_trips() {
        assert_round_trips(&format!("{TX}:n.100-1A>G"));
    }

    #[test]
    fn plus_offset_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.100+5_100+10del"));
    }

    #[test]
    fn minus_offset_range_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.100-5_100-1del"));
    }
}

// =============================================================================
// SECTION 5 — Insertions and boundary-crossing ranges
//
// Split by marker like §1/§2: the upstream crossings stay on `n.`, the
// downstream ones moved to `c.`, where the crossing is the real CDS→3'UTR one.
// =============================================================================

mod insertions_and_crossings {
    use super::*;

    #[test]
    fn upstream_to_interior_ins_round_trips() {
        assert_round_trips(&format!("{TX}:n.-1_1insAAA"));
    }

    #[test]
    fn interior_to_downstream_ins_round_trips() {
        assert_round_trips(&format!("{CDS_TX}:c.4443_*1insAAA"));
    }

    /// Boundary-spanning ranges are also covered by issue #253 (J1);
    /// pinned here to keep the J2 marker matrix self-contained.
    #[test]
    fn upstream_to_interior_del_round_trips() {
        assert_round_trips(&format!("{TX}:n.-3_5del"));
    }

    #[test]
    fn interior_to_downstream_del_round_trips() {
        assert_round_trips(&format!("{CDS_TX}:c.4440_*3del"));
    }
}
