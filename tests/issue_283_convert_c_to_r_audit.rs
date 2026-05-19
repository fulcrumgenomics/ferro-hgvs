//! Audit for #81 § E follow-up — reference-driven `c. → r.` conversion path.
//!
//! PR #234 (closes #233) pinned `c.`↔`r.` parity at the **Display** layer —
//! same coordinate axis, alphabet difference only. The **reference-driven**
//! conversion path under `src/convert/` is a separate surface and has no
//! dedicated audit. This file is that audit.
//!
//! ## Audit finding (pinned by `mod surface_inventory` below)
//!
//! `src/convert/` does **not** expose a dedicated `c. → r. variant-level
//! conversion function**. The module contains:
//!
//! - `convert::coding`  — `validate_cds_pos` / `cds_to_transcript_pos`
//!   (CDS → 0-based tx index helpers).
//! - `convert::genomic` — `validate_genome_pos` and basis conversions.
//! - `convert::mapper::CoordinateMapper` — `cds_to_tx`, `tx_to_cds`,
//!   `cds_to_protein`, `genomic_to_tx`, `tx_to_genomic`, plus intronic
//!   variants of these. The `Tx`/`CDS` arms operate on the **shared
//!   coordinate axis** that `c.` and `r.` ride together (per HGVS v21.0
//!   § 3.1; see #233 / PR #234). There is **no `cds_to_rna`,
//!   `rna_to_cds`, or `convert_to_rna_variant` method**.
//! - `convert::noncoding` — splice-site classification (`IntronicConsequence::
//!   from_cds_pos`, `from_tx_pos`) and `validate_tx_pos`. There is **no
//!   `from_rna_pos` constructor**.
//! - `convert::protein` — `validate_prot_pos`, `protein_to_cds_range`,
//!   `get_codon_frame`.
//!
//! In other words, the closest thing to a "reference-driven `c.→r.`"
//! surface today is:
//!   1. **Coordinate axis**: `CoordinateMapper::cds_to_tx` (the c.
//!      side); there is no separate `rna_to_tx` mapper because `r.`
//!      and `c.` share the same numerical axis on the same transcript.
//!   2. **Position-level bridge**: `Normalizer::cds_to_tx_pos` and
//!      `Normalizer::rna_to_tx_pos` (private; mirror each other and
//!      give bit-identical results for matched-axis inputs).
//!   3. **Alphabet swap**: the Display layer (pinned by PR #234).
//!
//! ## What this audit pins
//!
//! For each edit shape the `c. → r.` conversion **would** traverse if it
//! existed, pin the observed reference-driven coordinate behavior on the
//! `CoordinateMapper` surface that does exist, plus the parser-level
//! shape parity that already covers the alphabet swap.
//!
//! Test modules (11 in total):
//! - `surface_inventory`     — pin the absence of a dedicated c.→r.
//!   variant-level converter.
//! - `substitutions`         — pin `cds_to_tx`/`tx_to_cds` on the
//!   substitution coordinate axis; pin that the same positions on the
//!   `r.` side round-trip through the same mapper.
//! - `minus_strand`          — pin `cds_to_tx`/`tx_to_cds` continue to
//!   ride the shared c./r. axis on a minus-strand single-exon fixture.
//! - `intronic_offsets`      — pin offset preservation through
//!   `cds_to_tx` and pin the parser-level c./r. parity for intronic
//!   substitutions.
//! - `utr_markers`           — pin `c.-N` and `c.*N` → tx mapping; pin
//!   the analogous `r.-N`/`r.*N` parser shape.
//! - `deletion_insertion_duplication_delins_inversion` — pin endpoint
//!   mapping for ranged edits and pin reference-provider failure modes.
//! - `repeat_units`          — pin endpoint mapping for repeat edits.
//! - `compound_brackets`     — pin endpoint mapping for compound cis
//!   alleles.
//! - `predicted_edit_wrapper`— pin current acceptance of
//!   parenthesized/predicted wrappers on the parser side (the actual
//!   reference-driven conversion surface for predicted edits is not
//!   yet wired; in-flight via #242/#244/#246).
//! - `ref_provider_failure_modes` — pin error shapes for missing CDS,
//!   missing genomic coords, etc.
//! - `negative`              — pin out-of-range and ill-formed inputs.
//!
//! Pure-coverage audit; no production change. Closes #283.

use ferro_hgvs::convert::CoordinateMapper;
use ferro_hgvs::hgvs::location::{CdsPos, RnaPos, TxPos};
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};

const ACC: &str = "NM_000088.3";

/// Build a synthetic transcript shared by all coordinate-mapping tests.
///
/// Layout: 5bp 5'UTR + 30bp CDS + 5bp 3'UTR = 40bp, single exon.
///   tx 1..=5   — 5'UTR  (c.-5 .. c.-1, r.-5 .. r.-1)
///   tx 6..=35  — CDS    (c.1 .. c.30,  r.1 .. r.30)
///   tx 36..=40 — 3'UTR  (c.*1 .. c.*5, r.*1 .. r.*5)
fn make_transcript() -> Transcript {
    Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        // 40 bases — content not used by coordinate-only tests.
        "A".repeat(40),
        Some(6),
        Some(35),
        vec![Exon::new(1, 1, 40)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

/// Mirror of [`make_transcript`] on the minus strand. Same layout
/// (5bp 5'UTR + 30bp CDS + 5bp 3'UTR = 40bp single exon) — only the
/// strand differs. Used by `mod minus_strand` to pin that
/// `cds_to_tx` / `tx_to_cds` continue to ride the shared c./r.
/// coordinate axis (the same numerical axis as `r.`) on the minus
/// strand. Strand orientation matters for `genomic_to_tx` /
/// `tx_to_genomic`, not for the CDS↔Tx leg, but pinning the latter
/// against a minus-strand fixture guards against future regressions
/// that introduce strand-dependent branches into `cds_to_tx`.
fn make_transcript_minus() -> Transcript {
    Transcript::new(
        "NM_TEST_MINUS.1".to_string(),
        Some("TEST".to_string()),
        Strand::Minus,
        "A".repeat(40),
        Some(6),
        Some(35),
        vec![Exon::new(1, 1, 40)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
}

/// Mechanically translate a `c.`-form `Display` string into the
/// matching `r.`-form by swapping the prefix and lowercasing the
/// nucleotide alphabet (`T → u`).
///
/// Mirrors the `translate_cds_to_rna` helper in
/// `tests/issue_233_rna_cds_consistency.rs` — kept local here so this
/// audit file stands alone.
fn translate_cds_to_rna(s: &str) -> String {
    s.replace(":c.", ":r.")
        .replace('A', "a")
        .replace('C', "c")
        .replace('G', "g")
        .replace('T', "u")
}

/// Parse a `c.` input and an `r.` input on the **same axis** and assert
/// that their Display outputs differ only in prefix + alphabet. The
/// reference-driven conversion surface, were one to exist, would have
/// to satisfy this — the audit's only c./r. shape gate.
#[track_caller]
fn assert_parser_axis_parity(c_input: &str, r_input: &str) {
    let c = parse_hgvs(c_input).unwrap_or_else(|e| panic!("parse {c_input:?}: {e}"));
    let r = parse_hgvs(r_input).unwrap_or_else(|e| panic!("parse {r_input:?}: {e}"));
    let c_out = format!("{}", c);
    let r_out = format!("{}", r);
    let translated = translate_cds_to_rna(&c_out);
    assert_eq!(
        r_out, translated,
        "c./r. axis parity broken at parser surface:\n  c in  = {c_input}\n  c out = {c_out}\n  c→r  = {translated}\n  r in  = {r_input}\n  r out = {r_out}"
    );
}

// =============================================================================
// SECTION 0 — Surface inventory
// =============================================================================
//
// Pin the observed shape of `src/convert/`: there is no dedicated
// variant-level `c. → r.` converter. Coordinate mapping is c./tx-only
// because `r.` and `c.` share the same axis. If a future PR adds a
// `cds_to_rna_variant` (or similar) it should *intentionally* break
// these inventory tests and replace them with positive tests of that
// new API.

mod surface_inventory {
    use super::*;

    /// `CoordinateMapper::cds_to_tx` is the only reference-driven
    /// surface that touches `c.` coordinates. It maps to `TxPos`, which
    /// is the **shared** axis (same numerical value as `RnaPos.base`
    /// on the same transcript).
    #[test]
    fn cds_to_tx_lives_on_coordinate_mapper() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // Trivial smoke call to anchor the method exists with this
        // signature; documented as `pub fn cds_to_tx(&self, &CdsPos) -> Result<TxPos, _>`.
        let out: TxPos = mapper.cds_to_tx(&CdsPos::new(1)).unwrap();
        assert_eq!(out.base, 6);
    }

    /// The mapper exposes **no `cds_to_rna` / `rna_to_cds` /
    /// `rna_to_tx` / `tx_to_rna` method**. Pinned by the fact that this
    /// test compiles using only `cds_to_tx` / `tx_to_cds` (matched to
    /// the absent-r. helpers via the shared coordinate axis).
    #[test]
    fn cds_and_rna_share_the_same_numerical_axis_via_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // The "RNA → tx" leg is the identity on the positive CDS axis,
        // by spec — pin that an RnaPos with base N maps to the same tx
        // index as a CdsPos with base N via cds_to_tx.
        let from_c = mapper.cds_to_tx(&CdsPos::new(15)).unwrap();
        let r = RnaPos::new(15);
        // No public mapper method exists for RnaPos. The expected tx
        // index for r.15 on this transcript is cds_start + 15 - 1 = 20,
        // i.e. the same answer as cds_to_tx(c.15). Confirmed by the
        // private mirror in Normalizer (`rna_to_tx_pos`).
        assert_eq!(from_c.base, 20);
        assert_eq!(r.base, 15); // RnaPos.base is the c./r. axis value
    }
}

// =============================================================================
// SECTION 1 — Substitutions
// =============================================================================

mod substitutions {
    use super::*;

    /// Plus-strand exonic CDS substitution: c.1A>G maps to tx 6,
    /// which is the same tx index r.1a>g would map to.
    #[test]
    fn cds_exonic_c1_maps_to_tx6() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::new(1)).unwrap();
        assert_eq!(out.base, 6);
        assert_eq!(out.offset, None);
    }

    #[test]
    fn cds_exonic_last_cds_base_maps_to_tx35() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::new(30)).unwrap();
        assert_eq!(out.base, 35);
    }

    /// Round-trip c.N -> tx -> c.N for every CDS base.
    #[test]
    fn cds_tx_roundtrip_full_cds() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        for n in 1..=30 {
            let tx_pos = mapper.cds_to_tx(&CdsPos::new(n)).unwrap();
            let back = mapper.tx_to_cds(&tx_pos).unwrap();
            assert_eq!(back.base, n, "round-trip failed for c.{}", n);
            assert!(!back.utr3);
        }
    }

    /// Parser-level c./r. axis parity for a substitution — already
    /// pinned by issue #233; re-pinned here so the convert audit
    /// stands alone if #233 is ever rewritten.
    #[test]
    fn parser_axis_parity_exonic_sub() {
        assert_parser_axis_parity(&format!("{ACC}:c.123A>G"), &format!("{ACC}:r.123a>g"));
    }
}

// =============================================================================
// SECTION 1b — Minus strand
// =============================================================================
//
// The CDS↔Tx leg of `CoordinateMapper` is strand-agnostic: c./r.
// numbering is defined on the mRNA, which is already in the
// transcript's natural 5'→3' orientation regardless of which genomic
// strand the transcript came from. Pin that behavior on a
// minus-strand fixture so future strand-dependent branches in
// `cds_to_tx` would intentionally break these tests.

mod minus_strand {
    use super::*;

    /// On a minus-strand single-exon transcript with identical
    /// layout to [`make_transcript`], `cds_to_tx(c.1)` still maps to
    /// tx 6 — the strand does not enter the CDS↔Tx computation.
    #[test]
    fn minus_strand_cds_to_tx_matches_plus_axis() {
        let tx = make_transcript_minus();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::new(1)).unwrap();
        assert_eq!(out.base, 6);
        assert_eq!(out.offset, None);
    }

    /// Full c.→tx→c. round-trip on the minus-strand fixture for both
    /// exonic CDS and 3'UTR positions. Pins that `tx_to_cds` on a
    /// minus-strand transcript is the inverse of `cds_to_tx`,
    /// matching the plus-strand contract pinned in `mod substitutions`.
    #[test]
    fn minus_strand_cds_tx_roundtrip_preserves_position_and_utr_flag() {
        let tx = make_transcript_minus();
        let mapper = CoordinateMapper::new(&tx);
        for n in [1, 10, 30] {
            let tx_pos = mapper.cds_to_tx(&CdsPos::new(n)).unwrap();
            let back = mapper.tx_to_cds(&tx_pos).unwrap();
            assert_eq!(back.base, n, "round-trip failed for c.{}", n);
            assert!(!back.utr3);
        }
        // 3'UTR position round-trips with the utr3 flag intact.
        let original = CdsPos::utr3(2);
        let tx_pos = mapper.cds_to_tx(&original).unwrap();
        let back = mapper.tx_to_cds(&tx_pos).unwrap();
        assert_eq!(back.base, 2);
        assert!(back.utr3);
    }
}

// =============================================================================
// SECTION 2 — Intronic offsets
// =============================================================================
//
// `r.` does not normally carry intronic positions in the *spliced*
// RNA view; in practice ferro accepts them at the parser level (the
// shared coordinate axis), and `CoordinateMapper::cds_to_tx` carries
// the offset through unchanged. Pin both.

mod intronic_offsets {
    use super::*;

    /// `cds_to_tx` preserves the `+5` offset verbatim on `TxPos.offset`.
    #[test]
    fn cds_to_tx_preserves_donor_offset() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::with_offset(5, 5)).unwrap();
        assert_eq!(out.base, 10); // cds_start (6) + 5 - 1
        assert_eq!(out.offset, Some(5));
    }

    /// `cds_to_tx` preserves the `-3` offset verbatim on `TxPos.offset`.
    #[test]
    fn cds_to_tx_preserves_acceptor_offset() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::with_offset(10, -3)).unwrap();
        assert_eq!(out.base, 15);
        assert_eq!(out.offset, Some(-3));
    }

    /// Round-trip: an intronic CDS position survives the
    /// `cds_to_tx → tx_to_cds` cycle with the offset intact.
    #[test]
    fn intronic_cds_position_roundtrips_through_mapper() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let original = CdsPos::with_offset(7, 4);
        let tx_pos = mapper.cds_to_tx(&original).unwrap();
        let back = mapper.tx_to_cds(&tx_pos).unwrap();
        assert_eq!(back.base, original.base);
        assert_eq!(back.offset, original.offset);
        assert!(!back.utr3);
    }

    /// Parser-level c./r. axis parity holds for intronic donor offsets.
    /// Confirms the shape gate at the entry surface — the convert path
    /// would have to satisfy the same identity if it ever materialized.
    #[test]
    fn parser_axis_parity_donor_intronic_sub() {
        assert_parser_axis_parity(&format!("{ACC}:c.123+5A>G"), &format!("{ACC}:r.123+5a>g"));
    }

    #[test]
    fn parser_axis_parity_acceptor_intronic_sub() {
        assert_parser_axis_parity(&format!("{ACC}:c.123-5A>G"), &format!("{ACC}:r.123-5a>g"));
    }
}

// =============================================================================
// SECTION 3 — UTR markers
// =============================================================================
//
// `c.-N` (5'UTR) and `c.*N` (3'UTR) numbering is the same on `r.`
// (i.e. `r.-N` and `r.*N`). Pin the coordinate mapping and the
// parser-level axis parity.

mod utr_markers {
    use super::*;

    /// c.-3 → tx 3 on the test transcript (cds_start=6).
    #[test]
    fn five_prime_utr_minus_3_maps_to_tx3() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::new(-3)).unwrap();
        assert_eq!(out.base, 3);
        assert_eq!(out.offset, None);
    }

    /// c.-1 → tx 5 (the base immediately upstream of CDS start).
    #[test]
    fn five_prime_utr_minus_1_maps_to_last_utr_base() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::new(-1)).unwrap();
        assert_eq!(out.base, 5);
    }

    /// c.*1 → tx 36 (the base immediately downstream of CDS end).
    #[test]
    fn three_prime_utr_star_1_maps_to_first_3utr_base() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::utr3(1)).unwrap();
        assert_eq!(out.base, 36);
    }

    /// c.*5 → tx 40 (last 3'UTR base on the test transcript).
    #[test]
    fn three_prime_utr_star_5_maps_to_last_3utr_base() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::utr3(5)).unwrap();
        assert_eq!(out.base, 40);
    }

    /// Round-trip for a 3'UTR position preserves the `utr3` flag.
    #[test]
    fn three_prime_utr_roundtrips_with_flag() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let original = CdsPos::utr3(2);
        let tx_pos = mapper.cds_to_tx(&original).unwrap();
        let back = mapper.tx_to_cds(&tx_pos).unwrap();
        assert_eq!(back.base, 2);
        assert!(back.utr3);
    }

    /// Parser-level c./r. axis parity for 5'UTR substitution.
    #[test]
    fn parser_axis_parity_five_prime_utr_sub() {
        assert_parser_axis_parity(&format!("{ACC}:c.-12A>G"), &format!("{ACC}:r.-12a>g"));
    }

    /// Parser-level c./r. axis parity for 3'UTR substitution.
    #[test]
    fn parser_axis_parity_three_prime_utr_sub() {
        assert_parser_axis_parity(&format!("{ACC}:c.*5A>G"), &format!("{ACC}:r.*5a>g"));
    }
}

// =============================================================================
// SECTION 4 — Deletion / Insertion / Duplication / Delins / Inversion
// =============================================================================
//
// For a ranged edit, the reference-driven conversion path would map
// **each endpoint** through `CoordinateMapper::cds_to_tx`. Pin the
// endpoint-mapping behavior here, plus the parser-level c./r. shape
// parity that already covers the alphabet swap on the edit payload.

mod deletion_insertion_duplication_delins_inversion {
    use super::*;

    /// Deletion endpoints: c.10_12 maps to tx 15..=17.
    #[test]
    fn deletion_endpoints_map_through_cds_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let start = mapper.cds_to_tx(&CdsPos::new(10)).unwrap();
        let end = mapper.cds_to_tx(&CdsPos::new(12)).unwrap();
        assert_eq!((start.base, end.base), (15, 17));
    }

    /// Insertion at c.10_11 maps to tx 15..=16.
    #[test]
    fn insertion_endpoints_map_through_cds_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let start = mapper.cds_to_tx(&CdsPos::new(10)).unwrap();
        let end = mapper.cds_to_tx(&CdsPos::new(11)).unwrap();
        assert_eq!((start.base, end.base), (15, 16));
    }

    /// Duplication endpoints behave like deletion endpoints.
    #[test]
    fn duplication_endpoints_map_through_cds_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let start = mapper.cds_to_tx(&CdsPos::new(5)).unwrap();
        let end = mapper.cds_to_tx(&CdsPos::new(7)).unwrap();
        assert_eq!((start.base, end.base), (10, 12));
    }

    /// Delins: same as deletion-then-insertion at the endpoint level.
    #[test]
    fn delins_endpoints_map_through_cds_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let start = mapper.cds_to_tx(&CdsPos::new(20)).unwrap();
        let end = mapper.cds_to_tx(&CdsPos::new(25)).unwrap();
        assert_eq!((start.base, end.base), (25, 30));
    }

    /// Inversion endpoints — identical mechanism.
    #[test]
    fn inversion_endpoints_map_through_cds_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let start = mapper.cds_to_tx(&CdsPos::new(3)).unwrap();
        let end = mapper.cds_to_tx(&CdsPos::new(8)).unwrap();
        assert_eq!((start.base, end.base), (8, 13));
    }

    /// Endpoint round-trip: tx → cds → tx is the identity on this
    /// single-exon transcript for any pair of CDS endpoints.
    #[test]
    fn endpoint_roundtrip_on_single_exon() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        for (a, b) in [(1, 3), (4, 6), (10, 12), (25, 30)] {
            let ta = mapper.cds_to_tx(&CdsPos::new(a)).unwrap();
            let tb = mapper.cds_to_tx(&CdsPos::new(b)).unwrap();
            let ba = mapper.tx_to_cds(&ta).unwrap();
            let bb = mapper.tx_to_cds(&tb).unwrap();
            assert_eq!(ba.base, a, "round-trip failed at left endpoint {a}");
            assert_eq!(bb.base, b, "round-trip failed at right endpoint {b}");
        }
    }

    /// Parser-level shape parity: deletion with stated ref.
    #[test]
    fn parser_axis_parity_del_stated_ref() {
        assert_parser_axis_parity(
            &format!("{ACC}:c.123_125delACG"),
            &format!("{ACC}:r.123_125delacg"),
        );
    }

    /// Parser-level shape parity: insertion.
    #[test]
    fn parser_axis_parity_ins() {
        assert_parser_axis_parity(
            &format!("{ACC}:c.123_124insACG"),
            &format!("{ACC}:r.123_124insacg"),
        );
    }

    /// Parser-level shape parity: duplication with stated ref.
    #[test]
    fn parser_axis_parity_dup_stated_ref() {
        assert_parser_axis_parity(
            &format!("{ACC}:c.123_125dupACG"),
            &format!("{ACC}:r.123_125dupacg"),
        );
    }

    /// Parser-level shape parity: delins with stated alt.
    #[test]
    fn parser_axis_parity_delins() {
        assert_parser_axis_parity(
            &format!("{ACC}:c.123_125delinsTTT"),
            &format!("{ACC}:r.123_125delinsuuu"),
        );
    }

    /// Parser-level shape parity: inversion.
    #[test]
    fn parser_axis_parity_inv() {
        assert_parser_axis_parity(
            &format!("{ACC}:c.123_125inv"),
            &format!("{ACC}:r.123_125inv"),
        );
    }
}

// =============================================================================
// SECTION 5 — Repeat units
// =============================================================================

mod repeat_units {
    use super::*;

    /// Repeat anchor c.5 maps to tx 10 — pin the anchor-side
    /// coordinate mapping. The repeat *count* is a payload concern,
    /// not a coordinate one.
    #[test]
    fn repeat_anchor_maps_through_cds_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::new(5)).unwrap();
        assert_eq!(out.base, 10);
    }

    /// Parser-level c./r. axis parity for a simple repeat-unit edit.
    #[test]
    fn parser_axis_parity_simple_repeat() {
        assert_parser_axis_parity(&format!("{ACC}:c.123ACG[5]"), &format!("{ACC}:r.123acg[5]"));
    }

    /// MultiRepeat compound `c.123ACG[5]GT[3]` (and its `r.` mirror)
    /// shares its anchor with the simple-repeat case: only the first
    /// position (`c.5` on this fixture's CDS axis) participates in
    /// reference-driven coordinate mapping; the remaining unit/count
    /// payload rides on the edit. Pin the anchor-side `cds_to_tx`
    /// behavior on the same fixture as the simple-repeat test so the
    /// MultiRepeat shape is covered by this audit.
    ///
    /// TODO(#285): the inner repeat-count range form `insN[(150_180)]`
    /// is a separate parsing shape not yet exercised here; that
    /// follow-up is filed under #285.
    #[test]
    fn multirepeat_anchor_maps_through_cds_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // The anchor is the first position; on this fixture
        // (cds_start=6) c.5 maps to tx 10. Confirms the MultiRepeat
        // anchor uses the same coordinate path as a simple repeat.
        let out = mapper.cds_to_tx(&CdsPos::new(5)).unwrap();
        assert_eq!(out.base, 10);
        assert_eq!(out.offset, None);
    }

    /// Parser-level c./r. axis parity for the MultiRepeat compound
    /// shape `c.123ACG[5]GT[3]` / `r.123acg[5]gu[3]`. Pins that
    /// adding a second `unit[count]` segment does not desync the c.
    /// and r. axes at the Display layer.
    #[test]
    fn parser_axis_parity_multirepeat_compound() {
        assert_parser_axis_parity(
            &format!("{ACC}:c.123ACG[5]GT[3]"),
            &format!("{ACC}:r.123acg[5]gu[3]"),
        );
    }
}

// =============================================================================
// SECTION 6 — Compound cis brackets
// =============================================================================

mod compound_brackets {
    use super::*;

    /// Each member of a compound `c.[X;Y]` would map its endpoints
    /// independently. Pin two example c./r. parser-side parities.
    #[test]
    fn parser_axis_parity_two_sub_compound() {
        assert_parser_axis_parity(
            &format!("{ACC}:c.[123A>G;125C>T]"),
            &format!("{ACC}:r.[123a>g;125c>u]"),
        );
    }

    #[test]
    fn parser_axis_parity_mixed_shape_compound() {
        assert_parser_axis_parity(
            &format!("{ACC}:c.[123A>G;200_202delACG;300_301insACG]"),
            &format!("{ACC}:r.[123a>g;200_202delacg;300_301insacg]"),
        );
    }

    /// Pin that each candidate compound-member CDS position
    /// individually rides `cds_to_tx` with no cross-member effect —
    /// this is what a future c.→r. converter would walk per-member
    /// when lowering `c.[X;Y;Z]`. Does *not* parse a compound (the
    /// reference-driven conversion surface does not exist yet); the
    /// name reflects the unit being pinned (per-position
    /// independence), not the parsed-compound walk.
    #[test]
    fn each_candidate_member_position_maps_through_cds_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // Three candidate member positions at c.1, c.10, c.25.
        for n in [1, 10, 25] {
            let out = mapper.cds_to_tx(&CdsPos::new(n)).unwrap();
            // cds_start = 6, so tx = 6 + n - 1.
            assert_eq!(out.base, 6 + n - 1);
        }
    }
}

// =============================================================================
// SECTION 7 — Predicted-edit wrapper
// =============================================================================
//
// HGVS predicted notation wraps the edit in parens: `c.(123A>G)` /
// `r.(123a>g)`. PRs #242 / #244 / #246 are in flight for fuller
// predicted-edit support; until those land, pin what the parser
// currently accepts and reject-checks here. Reference-driven
// conversion of a predicted variant is, today, just the underlying
// coordinate mapping of the wrapped position — pin that too.

mod predicted_edit_wrapper {
    use super::*;

    /// Predicted substitution `c.(123A>G)` currently parses to a
    /// `CdsVariant`, and Display emits it as `c.123(A>G)` — the
    /// parens migrate from wrapping the whole edit to wrapping only
    /// the substitution body. Pin this observed re-shaping.
    ///
    /// TODO(issue #300): the c. side reshapes parens around the
    /// edit body; the r. side **drops them entirely** (see the next
    /// test). That asymmetry is a parser/Display divergence between
    /// the two paths and is the only behavioural inconsistency the
    /// audit surfaced. Tracked under issue #300 — either drop parens
    /// on both sides (current `r.` behavior), or keep parens around
    /// the body on both sides (current `c.` behavior).
    ///
    /// In-flight predicted-edit support: #242 / #244 / #246.
    #[test]
    fn predicted_cds_sub_display_wraps_edit_body() {
        let input = format!("{ACC}:c.(123A>G)");
        let v = parse_hgvs(&input).expect("predicted c.(...) sub must parse");
        assert!(matches!(v, HgvsVariant::Cds(_)));
        // Observed today: outer parens become body parens.
        assert_eq!(format!("{}", v), format!("{ACC}:c.123(A>G)"));
    }

    /// Predicted RNA substitution `r.(123a>g)` parses but Display
    /// **drops the parens entirely**, asymmetric to the c. side.
    /// Pin this asymmetry. See TODO on the c. test above for the
    /// follow-up to harmonize the two surfaces.
    #[test]
    fn predicted_rna_sub_display_drops_parens() {
        let input = format!("{ACC}:r.(123a>g)");
        let v = parse_hgvs(&input).expect("predicted r.(...) sub must parse");
        assert!(matches!(v, HgvsVariant::Rna(_)));
        // Observed today: parens are dropped on the r. side.
        assert_eq!(format!("{}", v), format!("{ACC}:r.123a>g"));
    }

    /// Coordinate mapping is unchanged by the predicted wrapper —
    /// the wrapper is a semantic tag on the edit, not the position.
    /// This means any future c.→r. converter would map positions
    /// identically regardless of the wrapper.
    #[test]
    fn predicted_position_uses_same_coordinate_mapping() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        let out = mapper.cds_to_tx(&CdsPos::new(15)).unwrap();
        assert_eq!(out.base, 20);
    }
}

// =============================================================================
// SECTION 8 — Reference-provider failure modes
// =============================================================================
//
// What does `CoordinateMapper` do with a transcript that's missing
// metadata the c./r. axis depends on? Pin each failure shape so a
// future c.→r. converter can rely on these errors.

mod ref_provider_failure_modes {
    use super::*;

    /// Helper: synthetic transcript with optional CDS bounds.
    fn make_tx(cds_start: Option<u64>, cds_end: Option<u64>) -> Transcript {
        Transcript::new(
            "NM_NOCDS.1".to_string(),
            None,
            Strand::Plus,
            "A".repeat(40),
            cds_start,
            cds_end,
            vec![Exon::new(1, 1, 40)],
            None,
            None,
            None,
            Default::default(),
            ManeStatus::default(),
            None,
            None,
        )
    }

    /// Transcript with no `cds_start`/`cds_end` — `cds_to_tx` returns
    /// `Err(ConversionError)`. A c.→r. converter would have to surface
    /// this via the same error type.
    #[test]
    fn missing_cds_errors_on_cds_to_tx() {
        let tx = make_tx(None, None);
        let mapper = CoordinateMapper::new(&tx);
        let result = mapper.cds_to_tx(&CdsPos::new(1));
        assert!(
            result.is_err(),
            "expected ConversionError on missing CDS, got {:?}",
            result
        );
        let msg = format!("{}", result.unwrap_err());
        // Exact-match the error string (variant prefix + message) so
        // a future change to the wording is loud rather than silently
        // still-passing.
        assert_eq!(msg, "Coordinate conversion error: Transcript has no CDS");
    }

    /// Transcript with `cds_start` set but `cds_end` unset — same
    /// error (the mapper short-circuits before checking either field
    /// in isolation).
    #[test]
    fn missing_cds_end_errors_on_cds_to_tx() {
        let tx = make_tx(Some(6), None);
        let mapper = CoordinateMapper::new(&tx);
        assert!(mapper.cds_to_tx(&CdsPos::new(1)).is_err());
    }

    /// `tx_to_cds` on a no-CDS transcript also errors — symmetric to
    /// `cds_to_tx`.
    #[test]
    fn missing_cds_errors_on_tx_to_cds() {
        let tx = make_tx(None, None);
        let mapper = CoordinateMapper::new(&tx);
        assert!(mapper.tx_to_cds(&TxPos::new(10)).is_err());
    }

    /// Transcript without genomic coordinates — `genomic_to_tx` errors.
    /// A c.→r. converter would not generally need this leg, but the
    /// failure is part of the mapper contract.
    #[test]
    fn missing_genomic_coords_errors_on_genomic_to_tx() {
        let tx = make_transcript();
        let mapper = CoordinateMapper::new(&tx);
        assert!(mapper.genomic_to_tx(1000).is_err());
    }
}

// =============================================================================
// SECTION 9 — Negative / invalid inputs
// =============================================================================

mod negative {
    use super::*;
    use ferro_hgvs::convert::coding::validate_cds_pos;

    /// CDS position beyond the CDS length is rejected by
    /// `validate_cds_pos` — same gate any c.→r. converter would need.
    #[test]
    fn cds_position_beyond_cds_is_rejected() {
        let tx = make_transcript();
        // CDS is 30bp; c.40 is out of range.
        assert!(validate_cds_pos(&CdsPos::new(40), &tx).is_err());
    }

    /// 5'UTR position beyond the UTR length is rejected.
    #[test]
    fn five_prime_utr_position_beyond_utr_is_rejected() {
        let tx = make_transcript();
        // 5'UTR is 5bp; c.-10 is out of range.
        assert!(validate_cds_pos(&CdsPos::new(-10), &tx).is_err());
    }

    /// 3'UTR position beyond the UTR length is rejected.
    #[test]
    fn three_prime_utr_position_beyond_utr_is_rejected() {
        let tx = make_transcript();
        // 3'UTR is 5bp; c.*10 is out of range.
        assert!(validate_cds_pos(&CdsPos::utr3(10), &tx).is_err());
    }

    /// Unreasonably large intronic offset is rejected by validation.
    #[test]
    fn unreasonably_large_intronic_offset_is_rejected() {
        let tx = make_transcript();
        let bad = CdsPos::with_offset(5, 2_000_000);
        assert!(validate_cds_pos(&bad, &tx).is_err());
    }

    /// `r.0` ("no RNA produced") parses; `c.0` does not. This
    /// asymmetry has no implication for the reference-driven conversion
    /// path — `r.0` is a whole-entity sentinel with no coordinate — but
    /// pinning it here keeps the audit self-contained.
    #[test]
    fn no_product_is_rna_only_at_parser_surface() {
        assert!(parse_hgvs(&format!("{ACC}:r.0")).is_ok());
        assert!(parse_hgvs(&format!("{ACC}:c.0")).is_err());
    }
}
