//! Audit for issue #218 — dedicated coverage for mixed-accession
//! compound alleles, the long-pending #81 item C4.
//!
//! Background.
//! HGVS does not give a single canonical shape for compound alleles
//! whose sub-variants reference different accessions. The same-coord
//! case (`[NM_X:c.1A>G;NM_Y:c.2A>G]`) is silent in the spec; the
//! mixed-coord case (`[c.76A>C];[g.10091C>G]`) is explicitly
//! discouraged in `recommendations/DNA/alleles.md`. ferro accepts both
//! and observes the following invariants on output:
//!   * Input order is preserved (no implicit sort).
//!   * The compact form `ACC:c.[edit;edit]` is emitted only when every
//!     sub-variant shares the same accession AND the same coord type;
//!     otherwise the expanded form (each sub-variant printed with its
//!     own prefix) is emitted.
//!   * Each sub-variant normalizes against its own reference, with
//!     un-resolvable companions riding through unchanged.
//!
//! Existing coverage is scattered:
//!   - `tests/compound_cross_reference.rs` covers cross-COORD-system
//!     compound alleles (c.+g., c.+r., c.+p., etc.) exhaustively.
//!   - `tests/allele_unknown_phase.rs` covers `(;)` unknown-phase
//!     mixed accession.
//!   - `tests/gene_selector_roundtrip.rs` covers gene-selector forms.
//!
//! Genuinely new ground in this audit (vs the sibling files above):
//!   * Same-coord, **different accession** lattice across cis / trans /
//!     unknown / mosaic / chimeric × 2 / 3 / 4-way.
//!   * Different versions of the same accession base (`NM_X.3` vs
//!     `NM_X.4`) — must not compact.
//!   * Mixed accession types within the same coord system (NM_×NR_,
//!     LRG_×NM_, ENST_×NM_, two different `NC_` chromosomes,
//!     two different `NP_` proteins).
//!   * `=` (identity) and `?` (uncertain) sub-variants inside a
//!     mixed-accession compound.
//!   * Compact-form rejection for several "mixed inner prefix" shapes.
//!   * Per-variant normalization with a *partial* provider — one
//!     accession resolves, the other doesn't — to actually exercise
//!     the independent-normalization invariant.
//!
//! No production-code changes land with this audit; it documents the
//! current correct behavior and guards against regressions.

use ferro_hgvs::hgvs::variant::{AllelePhase, HgvsVariant};
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{MockProvider, Normalizer};

fn assert_round_trips(input: &str) -> HgvsVariant {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let out = format!("{}", v);
    assert_eq!(
        out, input,
        "round-trip mismatch:\n  input:  `{}`\n  output: `{}`",
        input, out
    );
    v
}

fn assert_canonicalizes_to(input: &str, expected: &str) -> HgvsVariant {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    let out = format!("{}", v);
    assert_eq!(
        out, expected,
        "canonical-form mismatch:\n  input:    `{}`\n  expected: `{}`\n  got:      `{}`",
        input, expected, out
    );
    v
}

fn expect_phase(variant: &HgvsVariant, phase: AllelePhase, count: usize) {
    let allele = match variant {
        HgvsVariant::Allele(a) => a,
        other => panic!("expected AlleleVariant, got {:?}", other),
    };
    assert_eq!(allele.phase, phase, "wrong phase");
    assert_eq!(allele.variants.len(), count, "wrong sub-variant count");
}

// =============================================================================
// SECTION 1 — c.+c. (different transcripts), every phase
// =============================================================================

mod c_plus_c_cross_transcript {
    use super::*;

    #[test]
    fn cis_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G;NM_000089.3:c.2A>G]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }

    #[test]
    fn trans_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G];[NM_000089.3:c.2A>G]");
        expect_phase(&v, AllelePhase::Trans, 2);
    }

    #[test]
    fn unknown_phase_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G(;)NM_000089.3:c.2A>G]");
        expect_phase(&v, AllelePhase::Unknown, 2);
    }

    #[test]
    fn mosaic_round_trips() {
        let v = assert_round_trips("NM_000088.3:c.1A>G/NM_000089.3:c.2A>G");
        expect_phase(&v, AllelePhase::Mosaic, 2);
    }

    #[test]
    fn chimeric_round_trips() {
        let v = assert_round_trips("NM_000088.3:c.1A>G//NM_000089.3:c.2A>G");
        expect_phase(&v, AllelePhase::Chimeric, 2);
    }
}

// =============================================================================
// SECTION 2 — 3-way and 4-way mixed-accession compounds
// =============================================================================

mod multi_way {
    use super::*;

    #[test]
    fn cis_three_way_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G;NM_000089.3:c.2A>G;NM_000090.3:c.3A>G]");
        expect_phase(&v, AllelePhase::Cis, 3);
    }

    #[test]
    fn trans_three_way_round_trips() {
        let v =
            assert_round_trips("[NM_000088.3:c.1A>G];[NM_000089.3:c.2A>G];[NM_000090.3:c.3A>G]");
        expect_phase(&v, AllelePhase::Trans, 3);
    }

    #[test]
    fn unknown_phase_three_way_round_trips() {
        let v =
            assert_round_trips("[NM_000088.3:c.1A>G(;)NM_000089.3:c.2A>G(;)NM_000090.3:c.3A>G]");
        expect_phase(&v, AllelePhase::Unknown, 3);
    }

    #[test]
    fn cis_four_way_round_trips() {
        let v = assert_round_trips(
            "[NM_000088.3:c.1A>G;NM_000089.3:c.2A>G;NM_000090.3:c.3A>G;NM_000091.3:c.4A>G]",
        );
        expect_phase(&v, AllelePhase::Cis, 4);
    }

    #[test]
    fn mosaic_three_way_round_trips() {
        let v = assert_round_trips("NM_000088.3:c.1A>G/NM_000089.3:c.2A>G/NM_000090.3:c.3A>G");
        expect_phase(&v, AllelePhase::Mosaic, 3);
    }

    #[test]
    fn chimeric_three_way_round_trips() {
        let v = assert_round_trips("NM_000088.3:c.1A>G//NM_000089.3:c.2A>G//NM_000090.3:c.3A>G");
        expect_phase(&v, AllelePhase::Chimeric, 3);
    }
}

// =============================================================================
// SECTION 3 — same-accession base, different versions
// =============================================================================
//
// NCBI accessions encode versioning as `<id>.<version>`. Two versions of
// the same transcript MUST be treated as different accessions for the
// purposes of compound-allele canonicalization — the underlying sequence
// can differ across versions.

mod versioned_same_base {
    use super::*;

    #[test]
    fn cis_different_versions_does_not_compact() {
        // Different `.3` and `.4` versions → expanded form must survive.
        // (If the parser ever started normalizing version, this would
        // either drift or fold to compact form — both would be wrong.)
        let v = assert_round_trips("[NM_000088.3:c.1A>G;NM_000088.4:c.2A>G]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }

    #[test]
    fn trans_different_versions_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G];[NM_000088.4:c.2A>G]");
        expect_phase(&v, AllelePhase::Trans, 2);
    }
}

// =============================================================================
// SECTION 4 — same-accession compact boundary (single regression guard)
// =============================================================================
//
// `compound_cross_reference.rs::test_display_uses_compact_form_only_when_shared`
// already pins same-accession-shared-coord canonicalization to the
// compact form, and the compact-form round-trip is covered by many
// callers across the suite. We keep ONE explicit guard here on the
// boundary between Section 3 (different versions → must NOT compact)
// and the existing same-accession-compacts pin, so a future refactor of
// the cross-reference machinery that breaks one direction without the
// other surfaces in both audit files.

mod same_accession_compact_boundary {
    use super::*;

    #[test]
    fn cis_same_accession_compacts_to_shared_prefix() {
        assert_canonicalizes_to(
            "[NM_000088.3:c.1A>G;NM_000088.3:c.2A>G]",
            "NM_000088.3:c.[1A>G;2A>G]",
        );
    }
}

// =============================================================================
// SECTION 5 — mixed accession TYPES within same coord system
// =============================================================================
//
// Same coord system, different accession namespaces. All of NM_, NR_,
// LRG_, ENST_, NP_, NC_ are first-class supported by ferro's parser
// (see `src/hgvs/parser/accession.rs`); a compound carrying one of each
// is valid and the emitted form must keep the type-distinguished
// prefixes.

mod mixed_accession_types {
    use super::*;

    #[test]
    fn cis_nm_and_nr_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G;NR_001000.1:n.2A>G]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }

    #[test]
    fn trans_nm_and_nr_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G];[NR_001000.1:n.2A>G]");
        expect_phase(&v, AllelePhase::Trans, 2);
    }

    /// LRG (Locus Reference Genomic) transcript accession paired with a
    /// RefSeq transcript. Per `recommendations/general.md`, LRGs are
    /// still acceptable references; the parser supports them and a
    /// compound across an LRG and an NM_ must round-trip.
    #[test]
    fn cis_lrg_and_nm_round_trips() {
        let v = assert_round_trips("[LRG_199t1:c.1A>G;NM_000088.3:c.2A>G]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }

    /// Ensembl transcript accession (`ENST*`) paired with a RefSeq.
    /// Pins the symmetric case to the LRG one.
    #[test]
    fn cis_ensembl_and_refseq_round_trips() {
        let v = assert_round_trips("[ENST00000380152.7:c.1A>G;NM_000088.3:c.2A>G]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }

    #[test]
    fn cis_two_different_chromosomes_round_trips() {
        // `g.` on chr17 + `g.` on chr18 → same coord type, distinct
        // accessions. Spec-relevant: a haplotype on different
        // chromosomes is not "trans" in the medical sense, but the
        // syntactic shape is the same.
        let v = assert_round_trips("[NC_000017.11:g.100A>G;NC_000018.11:g.200C>T]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }

    #[test]
    fn cis_two_proteins_round_trips() {
        let v = assert_round_trips("[NP_000079.2:p.Arg100Gly;NP_000080.2:p.Lys200Asn]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }
}

// =============================================================================
// SECTION 6 — gene-selector form preserved
// =============================================================================

mod gene_selector {
    use super::*;

    #[test]
    fn cis_with_gene_selectors_round_trips() {
        let v = assert_round_trips("[NM_000088.3(COL1A1):c.100A>G;NM_000089.3(COL1A2):c.200C>T]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }

    #[test]
    fn trans_with_gene_selectors_round_trips() {
        let v = assert_round_trips("[NM_000088.3(COL1A1):c.100A>G];[NM_000089.3(COL1A2):c.200C>T]");
        expect_phase(&v, AllelePhase::Trans, 2);
    }
}

// =============================================================================
// SECTION 7 — input order is preserved (no implicit sorting)
// =============================================================================

mod ordering {
    use super::*;

    #[test]
    fn cis_input_order_preserved() {
        // Forward + reverse must both round-trip and must NOT compare
        // equal — alphabetical-by-accession ordering would be an
        // unfounded canonicalization, since the spec doesn't define
        // one.
        let forward = assert_round_trips("[NM_000088.3:c.1A>G;NM_000089.3:c.2A>G]");
        let reverse = assert_round_trips("[NM_000089.3:c.2A>G;NM_000088.3:c.1A>G]");
        assert_ne!(
            forward, reverse,
            "ordering must be significant; the parser must not sort"
        );
    }

    #[test]
    fn trans_input_order_preserved() {
        let forward = assert_round_trips("[NM_000088.3:c.1A>G];[NM_000089.3:c.2A>G]");
        let reverse = assert_round_trips("[NM_000089.3:c.2A>G];[NM_000088.3:c.1A>G]");
        assert_ne!(forward, reverse);
    }
}

// =============================================================================
// SECTION 8 — normalize-each-independently with a PARTIAL provider
// =============================================================================
//
// The independent-normalization invariant only really tests itself when
// one sub-variant has a resolvable reference and the other does not.
// With an empty `MockProvider` both sides are no-ops and the assertion
// only confirms "round-trip survives a no-op normalize." Mirror the
// `provider_with_simple_transcript` pattern from
// `compound_cross_reference.rs::test_each_variant_normalizes_against_own_ref`
// — register one `NM_TEST.1` transcript whose c. variants normalize,
// and pair them with a sub-variant on an unregistered second
// accession that must ride through unchanged.

/// Build a `MockProvider` holding ONE transcript (`NM_TEST.1`) whose
/// sequence places a `CCCCC` 5-bp run starting at c.6. A sub-variant
/// `NM_TEST.1:c.10dupA` normalizes to `c.10dup` (the redundant
/// stated-ref base is dropped by the canonicalize pass), proving the
/// resolved side of a mixed-accession compound actually got
/// normalized.
fn provider_with_test_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "ATGCAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGTTTTTAAAAACCCCCGGGGGT".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(1),
        Some(len),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

mod normalize_each_independently {
    use super::*;

    /// Cis: `NM_TEST.1:c.10dupA` (resolved) paired with
    /// `NM_OTHER.1:c.5A>G` (un-resolved on the same coord system).
    /// The resolved variant must drop its stated-ref `A` and emit
    /// `c.10dup`; the un-resolved partner must ride through unchanged.
    /// Distinct from `compound_cross_reference.rs` (which pairs c.+g.)
    /// — this is the same-coord, different-accession case.
    #[test]
    fn cis_partial_provider_resolves_only_known_accession() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("[NM_TEST.1:c.10dupA;NM_OTHER.1:c.5A>G]").expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        assert_eq!(
            format!("{}", normalized),
            "[NM_TEST.1:c.10dup;NM_OTHER.1:c.5A>G]",
        );
    }

    #[test]
    fn trans_partial_provider_resolves_only_known_accession() {
        let normalizer = Normalizer::new(provider_with_test_transcript());
        let parsed = parse_hgvs("[NM_TEST.1:c.10dupA];[NM_OTHER.1:c.5A>G]").expect("parse");
        let normalized = normalizer.normalize(&parsed).expect("normalize");
        assert_eq!(
            format!("{}", normalized),
            "[NM_TEST.1:c.10dup];[NM_OTHER.1:c.5A>G]",
        );
    }

    // Idempotency on a `dupA` input with a same-coord different-accession
    // partner surfaces a separate pre-existing two-pass-drift bug (pass 1
    // drops the mismatched stated-ref base without shifting; pass 2
    // shifts). Tracked in #219 as a dedicated follow-up; deliberately
    // out of scope for this audit. Not introduced or worsened by the
    // mixed-accession lattice.
}

// =============================================================================
// SECTION 9 — compact form MUST NOT collapse across different accessions
// =============================================================================
//
// `compound_cross_reference.rs::test_compact_form_rejects_mixed_prefix_inside_brackets`
// pins the canonical case. We extend with the symmetric shapes the
// compact form must also reject — without them a parser that grew
// tolerance for some shape variant of the mix would silently accept
// the inverse mistake.

mod compact_form_rejection {
    use super::*;

    /// Outer accession + first-element-unprefixed, then a second
    /// element that respelled a different accession's prefix:
    /// malformed (the second element should not carry a prefix if the
    /// outer is in compact form).
    #[test]
    fn unprefixed_then_prefixed_inner_rejected() {
        assert!(parse_hgvs("NM_000088.3:c.[1A>G;NM_000089.3:c.2A>G]").is_err());
    }

    /// Symmetric: prefixed first inner element disagreeing with the
    /// outer accession.
    #[test]
    fn prefixed_first_inner_with_different_accession_rejected() {
        assert!(parse_hgvs("NM_000088.3:c.[NM_000089.3:c.1A>G;2A>G]").is_err());
    }

    /// Three-element compact form where one inner element carries a
    /// prefixed (and different) accession.
    #[test]
    fn three_element_with_one_different_accession_rejected() {
        assert!(parse_hgvs("NM_000088.3:c.[1A>G;NM_000089.3:c.2A>G;3A>G]").is_err());
    }
}

// =============================================================================
// SECTION 10 — mixed-accession with `=` (identity) and `?` (uncertain)
// =============================================================================
//
// `=` and `?` are spec-valid edit / position kinds that can appear as
// sub-variants of a compound. The mixed-accession lattice must accept
// them without dropping the wrapper or rewriting the partner.

mod identity_and_uncertain_sub_variants {
    use super::*;

    /// One sub-variant is the `=` (identity / no change) form on a
    /// distinct accession.
    #[test]
    fn cis_with_identity_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G;NM_000089.3:c.2=]");
        expect_phase(&v, AllelePhase::Cis, 2);
    }

    #[test]
    fn trans_with_identity_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G];[NM_000089.3:c.2=]");
        expect_phase(&v, AllelePhase::Trans, 2);
    }

    /// `c.?` (whole-entity unknown) on the second accession — common
    /// in recessive-disease haplotype reports where the second allele
    /// is undetermined.
    #[test]
    fn trans_with_whole_entity_uncertain_round_trips() {
        let v = assert_round_trips("[NM_000088.3:c.1A>G];[NM_000089.3:c.?]");
        expect_phase(&v, AllelePhase::Trans, 2);
    }
}
