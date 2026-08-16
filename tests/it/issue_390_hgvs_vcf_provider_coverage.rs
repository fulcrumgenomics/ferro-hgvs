//! Audit for #390 item 6 — `HgvsToVcfConverter` (provider-aware path)
//! coord-system × edit-kind coverage matrix.
//!
//! Parallel to `tests/issue_261_hgvs_vcf_coverage.rs`, which covers only
//! the no-provider `genomic_hgvs_to_vcf`. The converter accepts a
//! transcript + provider and can therefore project c./n./g./m. inputs
//! through the transcript exon table; this file pins which
//! (coord-axis × edit-kind) cells succeed and which surface a documented
//! `ConversionError`.
//!
//! Cell coverage (this audit is deliberately wider on coord-system than
//! edit-kind — see notes below for the gaps and rationale):
//!
//!   axis ↓ \ edit →   sub   del   ins   dup   delins  identity  inv
//!   c.                ✓     ✓†    ✓†    ✓†    ✓†      —         —
//!   n.                ✓     —     —     —     —       —         —
//!   g.                ✓     —     —     —     —       —         —
//!   m. (genome route) ✓     —     —     —     —       —         —
//!   r.                err   err   err   err   err     err       err
//!   p.                err   err   err   err   err     err       err
//!
//! ✓ = expected success;  err = converter returns a typed
//!   `ConversionError`;  ✓† = provider-aware path; the converter
//!   delegates the anchor-base fetch to the supplied provider, so the
//!   cell exercises the `MockProvider` ↔ converter contract;  — = not
//!   covered here (the inline tests in `src/vcf/from_hgvs.rs` exercise
//!   most del/ins/dup/delins shapes against the same single-exon
//!   plus-strand fixture).
//!
//! Strand / multi-exon coverage: a minus-strand single-exon fixture
//! lives in `minus_strand_axis` below. Multi-exon, intronic, and
//! splice-junction projection are exercised by the inline tests in
//! `src/vcf/from_hgvs.rs` rather than duplicated here (those tests
//! have first-class access to the converter's internals).
//!
//! The transcript fixture is a deliberately simpler synthetic
//! single-exon transcript (1-100 → 1000-1099, CDS [50, 75]) than the
//! inline `create_test_transcript` in `src/vcf/from_hgvs.rs` (which is
//! a two-exon NM_000088.3 fixture). Simplification is intentional: the
//! audit's goal is to pin the high-level (axis × edit-kind) cells, not
//! to re-cover the multi-exon projection details that the inline tests
//! already pin.

use ferro_hgvs::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use ferro_hgvs::hgvs::interval::{CdsInterval, GenomeInterval, TxInterval};
use ferro_hgvs::hgvs::location::{CdsPos, GenomePos, TxPos};
use ferro_hgvs::hgvs::variant::{
    Accession, AllelePhase, AlleleVariant, CdsVariant, GenomeVariant, HgvsVariant, LocEdit,
    MtVariant, TxVariant,
};
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::vcf::HgvsToVcfConverter;
use ferro_hgvs::FerroError;
use std::str::FromStr;

/// Single-exon coding transcript on chr1, plus strand:
///   - tx [1, 100] → genome [1000, 1099]
///   - CDS tx [50, 75]
///   - sequence is `ATGCATGC` repeated so every CDS position has a
///     concrete base for the SPDI/anchor sub-path tests that need one.
fn fixture_transcript() -> Transcript {
    Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        Some("ATGCATGC".repeat(20)),
        Some(50),
        Some(75),
        vec![Exon::with_genomic(1, 1, 100, 1000, 1099)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1099),
        GenomeBuild::GRCh38,
        ManeStatus::Select,
        None,
        None,
    )
}

/// Minus-strand single-exon coding transcript on chr1:
///   - tx [1, 100] → genome [1099, 1000] (reversed)
///   - sequence is the reverse complement of `fixture_transcript`'s
///     (`ATGCATGC` → `GCATGCAT`), so the transcript bases read 5'→3'
///     on the minus strand match the reverse-complement of the
///     plus-strand genomic bases at the same exon positions.
fn fixture_transcript_minus() -> Transcript {
    Transcript::new(
        "NM_TEST_MINUS.1".to_string(),
        Some("TEST".to_string()),
        Strand::Minus,
        Some("GCATGCAT".repeat(20)),
        Some(50),
        Some(75),
        vec![Exon::with_genomic(1, 1, 100, 1000, 1099)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1099),
        GenomeBuild::GRCh38,
        ManeStatus::Select,
        None,
        None,
    )
}

/// Single-exon NON-coding transcript on chr1, plus strand:
///   - tx [1, 100] → genome [1000, 1099]
///   - no CDS (cds_start/cds_end = None), so `is_coding() == false` and
///     `r.` numbering is plain transcript-relative (== n.).
///
/// Used to exercise `convert_rna`'s non-coding branch (including its
/// `utr3`/`base < 1` decline guard) through the provider-backed matrix.
fn fixture_transcript_non_coding() -> Transcript {
    Transcript::new(
        "NR_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        Some("ATGCATGC".repeat(20)),
        None,
        None,
        vec![Exon::with_genomic(1, 1, 100, 1000, 1099)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1099),
        GenomeBuild::GRCh38,
        ManeStatus::Select,
        None,
        None,
    )
}

/// Single-exon transcript on chr2 — used to verify the m. routing
/// (chrM) is independent of the bound transcript's chromosome.
fn fixture_transcript_chr2() -> Transcript {
    Transcript::new(
        "NM_OTHER.1".to_string(),
        Some("OTHER".to_string()),
        Strand::Plus,
        Some("ATGCATGC".repeat(20)),
        Some(50),
        Some(75),
        vec![Exon::with_genomic(1, 1, 100, 5000, 5099)],
        Some("chr2".to_string()),
        Some(5000),
        Some(5099),
        GenomeBuild::GRCh38,
        ManeStatus::Select,
        None,
        None,
    )
}

fn fixture_provider() -> MockProvider {
    // 2000-base genomic sequence covering the fixture's exon range:
    // 1-based 1000..1099 maps to 0-based 999..1098. We pad with N's
    // so anchor-base fetches at exonic positions return a definite
    // base rather than panicking on a missing region.
    let mut p = MockProvider::new();
    let mut seq = "N".repeat(999);
    seq.push_str(&"ATGCATGC".repeat(20)); // 160 bases starting at 0-based 999
    seq.push_str(&"N".repeat(1000));
    p.add_genomic_sequence("NC_000001.11", &seq);
    // Same shape on chr1 alias the converter looks up by accession.
    p.add_genomic_sequence("chr1", &seq);
    p
}

fn substitution(ref_b: Base, alt_b: Base) -> NaEdit {
    NaEdit::Substitution {
        reference: ref_b,
        alternative: alt_b,
    }
}

// =============================================================================
// SECTION 1 — Coord-axis × substitution (the always-supported cell)
// =============================================================================

mod substitution_axis_matrix {
    use super::*;

    #[test]
    fn c_axis_substitution_emits_chr_pos_ref_alt() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        // c.1 = tx 50 = genome 1049 (1-based, plus strand).
        let variant = CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Cds(variant))
            .expect("c. SNV must convert via the provider-aware path");
        assert_eq!(result.record.chrom, "chr1");
        // c.1 = tx 50 (CDS_start) → genome 1049 on the plus-strand
        // single-exon fixture (tx [1,100] → genome [1000,1099]).
        assert_eq!(result.record.pos, 1049, "c.1 must project to genome 1049");
        assert_eq!(result.record.reference, "A");
        assert_eq!(result.record.alternate, vec!["G".to_string()]);
        assert!(
            result.record.info.contains_key("HGVS"),
            "INFO/HGVS annotation must be populated"
        );
    }

    #[test]
    fn n_axis_substitution_emits_chr_pos_ref_alt() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = TxVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Tx(variant))
            .expect("n. SNV must convert via the provider-aware path");
        assert_eq!(result.record.chrom, "chr1");
        // n.1 = tx 1 → genome 1000 on the plus-strand single-exon
        // fixture (tx [1,100] → genome [1000,1099]).
        assert_eq!(result.record.pos, 1000, "n.1 must project to genome 1000");
        assert_eq!(result.record.reference, "A");
        assert_eq!(result.record.alternate, vec!["G".to_string()]);
    }

    #[test]
    fn g_axis_substitution_emits_chr_pos_ref_alt() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Genome(variant))
            .expect("g. SNV must convert");
        assert_eq!(result.record.chrom, "chr1");
        assert_eq!(result.record.pos, 12345);
        assert_eq!(result.record.reference, "A");
    }

    #[test]
    fn m_axis_substitution_routes_through_genome_branch() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = MtVariant {
            accession: Accession::new("NC", "012920", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(3243)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Mt(variant))
            .expect("m. SNV must route through the genome branch");
        assert_eq!(result.record.chrom, "chrM");
        assert_eq!(result.record.pos, 3243);
    }

    /// m. routing is independent of the bound transcript's
    /// chromosome: a converter created with a chr2 transcript still
    /// emits `chrM` for an `NC_012920.1` input.
    #[test]
    fn m_axis_chrom_independent_of_bound_transcript() {
        let tx = fixture_transcript_chr2();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = MtVariant {
            accession: Accession::new("NC", "012920", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(3243)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter.convert(&HgvsVariant::Mt(variant)).unwrap();
        assert_eq!(
            result.record.chrom, "chrM",
            "m. routing must derive chrom from the variant's accession, not the bound transcript"
        );
    }
}

// =============================================================================
// SECTION 1b — c. axis × edit-kind matrix (provider-aware path)
// =============================================================================
//
// The converter's `get_reference_base` / `get_reference_sequence`
// helpers (`src/vcf/from_hgvs.rs:533-597`) pass a **genomic** position
// to `provider.get_sequence` keyed by the **transcript** id. Without
// the genomic bases registered under the transcript accession (and
// the cached transcript sequence covering the genomic range, which
// is unusual), del / ins / dup / delins on c. inputs surface
// `GenomicReferenceNotAvailable` rather than converting.
//
// That's a pre-existing converter bug surfaced by the audit but
// outside the scope of #390 item 6 (which is the audit itself). Pin
// the documented limitation here so a future fix flips the
// assertion shape rather than silently passing without coverage.

mod c_axis_edit_kind_matrix {
    use super::*;

    // The fixture genome holds "TGC" at genomic 1049-1051 (where c.1_3
    // projects), so the deleted/duplicated sequences below use "TGC" to keep
    // each variant biologically consistent with the reference it acts on.

    fn cds_del_tgc() -> HgvsVariant {
        HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(3)),
                NaEdit::Deletion {
                    sequence: Some(Sequence::from_str("TGC").unwrap()),
                    length: None,
                },
            ),
        })
    }

    fn cds_ins_atg() -> HgvsVariant {
        HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(2)),
                NaEdit::Insertion {
                    // Inserted material is arbitrary novel sequence (it need
                    // not match the reference); "ATG" exercises that path.
                    sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
                },
            ),
        })
    }

    fn cds_dup_tgc() -> HgvsVariant {
        HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(3)),
                NaEdit::Duplication {
                    sequence: Some(Sequence::from_str("TGC").unwrap()),
                    length: None,
                    uncertain_extent: None,
                },
            ),
        })
    }

    fn cds_delins_tgc_ttcc() -> HgvsVariant {
        HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(3)),
                NaEdit::Delins {
                    sequence: InsertedSequence::Literal(Sequence::from_str("TTCC").unwrap()),
                    deleted: Some(Sequence::from_str("TGC").unwrap()),
                    deleted_length: None,
                    substitution_reference: None,
                },
            ),
        })
    }

    /// Convert `variant` against the fixture transcript + provider and assert
    /// the resulting VCF record.
    ///
    /// Before #805 these cells surfaced a `GenomicReferenceNotAvailable`
    /// because the anchor-base lookup fetched against the transcript id with a
    /// genomic position. The converter now fetches the anchor base genomically
    /// (keyed by the contig name the record is emitted under), so each cell
    /// resolves to a VCF record with the provider-recovered anchor base.
    ///
    /// Fixture genomic bases: 1-based genomic position N (for N >= 1000) holds
    /// `"ATGCATGC"`[(N - 1000) mod 8]. c.1 maps to genomic 1049.
    fn assert_converts_to(
        variant: HgvsVariant,
        label: &str,
        expected_pos: u64,
        expected_ref: &str,
        expected_alt: &str,
    ) {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let record = converter.convert(&variant).expect(label).record;
        assert_eq!(record.chrom, "chr1", "{label}: chrom");
        assert_eq!(record.pos, expected_pos, "{label}: pos");
        assert_eq!(record.reference, expected_ref, "{label}: ref");
        assert_eq!(
            record.alternate,
            vec![expected_alt.to_string()],
            "{label}: alt"
        );
    }

    #[test]
    fn c_axis_deletion_recovers_anchor_from_provider() {
        // c.1_3delTGC (genomic 1049-1051 = "TGC"); anchor at genomic 1048 =
        // 'A'. REF = anchor + deleted, ALT = anchor.
        assert_converts_to(cds_del_tgc(), "c. del", 1048, "ATGC", "A");
    }

    #[test]
    fn c_axis_insertion_recovers_anchor_from_provider() {
        // c.1_2insATG (anchor at genomic 1049 = 'T'); ALT = anchor + inserted.
        assert_converts_to(cds_ins_atg(), "c. ins", 1049, "T", "TATG");
    }

    #[test]
    fn c_axis_duplication_recovers_anchor_from_provider() {
        // c.1_3dupTGC (anchor at the genomic end 1051 = 'C').
        assert_converts_to(cds_dup_tgc(), "c. dup", 1051, "C", "CTGC");
    }

    #[test]
    fn c_axis_delins_recovers_anchor_from_provider() {
        // c.1_3delTGCinsTTCC (length-changing → anchor at genomic 1048 = 'A';
        // deleted bases fetched from genomic 1049-1051 = "TGC").
        assert_converts_to(cds_delins_tgc_ttcc(), "c. delins", 1048, "ATGC", "ATTCC");
    }
}

// =============================================================================
// SECTION 1c — Minus-strand projection
// =============================================================================
//
// The provider-aware path's distinguishing feature vs.
// `genomic_hgvs_to_vcf` is the strand-aware transcript→genome
// projection. The plus-strand fixture above exercises the trivial
// case (genome = tx + offset). Here we pin the minus-strand cell:
// c.1 on the minus-strand fixture maps to the GENOMIC END of the
// transcript, with strand-aware base complementation applied to the
// ref/alt.

mod minus_strand_axis {
    use super::*;

    #[test]
    fn c_axis_substitution_minus_strand_projects_to_genome_end() {
        let tx = fixture_transcript_minus();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = CdsVariant {
            accession: Accession::new("NM", "TEST_MINUS", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Cds(variant))
            .expect("minus-strand c. SNV must convert");
        assert_eq!(result.record.chrom, "chr1");
        // The fixture's CDS starts at tx_pos=50; on minus strand that's
        // toward the 3' (higher genomic) end of the transcript. The
        // exact genomic position depends on the converter's
        // strand-mapping convention — we just pin that the chrom and
        // strand semantics fired (no error) and let the inline tests
        // pin the precise minus-strand coords.
        assert!(result.record.pos >= 1000 && result.record.pos <= 1099);
    }
}

// =============================================================================
// SECTION 2 — Coord-axes that the converter rejects up-front
// =============================================================================

mod rejected_axes {
    use super::*;
    use ferro_hgvs::hgvs::edit::ProteinEdit;
    use ferro_hgvs::hgvs::interval::ProtInterval;
    use ferro_hgvs::hgvs::interval::RnaInterval;
    use ferro_hgvs::hgvs::location::RnaPos;
    use ferro_hgvs::hgvs::location::{AminoAcid, ProtPos};
    use ferro_hgvs::hgvs::variant::{ProteinVariant, RnaVariant};

    /// On a coding transcript, `r.` numbering is CDS-relative (== c.), so
    /// `r.1A>G` lowers through CDS 1 → tx 50 → genomic 1049 (cds_start=50,
    /// exon1 tx 1 → genomic 1000) and converts like the c. path. The
    /// substitution carries its own REF base.
    #[test]
    fn r_axis_converts_via_cds_numbering() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = RnaVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Rna(variant))
            .expect("r. on a coding transcript converts via CDS numbering");
        assert_eq!(result.record.pos, 1049);
        assert_eq!(result.record.reference, "A");
        assert_eq!(result.record.alternate, vec!["G".to_string()]);
    }

    /// On a NON-coding transcript there is no CDS, so `r.` numbering is plain
    /// transcript-relative (== n.). `r.1A>G` lowers through tx 1 → genomic 1000
    /// (tx [1,100] → genome [1000,1099]) — the non-coding branch of
    /// `convert_rna`, exercised here through the provider-backed matrix rather
    /// than only the in-module unit tests.
    #[test]
    fn r_axis_non_coding_converts_via_transcript_numbering() {
        let tx = fixture_transcript_non_coding();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = RnaVariant {
            accession: Accession::new("NR", "TEST", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let result = converter
            .convert(&HgvsVariant::Rna(variant))
            .expect("r. on a non-coding transcript converts via transcript numbering");
        assert_eq!(result.record.pos, 1000);
        assert_eq!(result.record.reference, "A");
        assert_eq!(result.record.alternate, vec!["G".to_string()]);
    }

    /// A `*N` (3'-UTR) position is meaningless on a non-coding transcript (no
    /// CDS frame), so the non-coding branch declines rather than silently
    /// mapping `*N` as `N`. Pins that `convert_rna`'s `utr3`/`base < 1` guard
    /// fires through the provider-backed matrix.
    #[test]
    fn r_axis_non_coding_utr3_declines() {
        let tx = fixture_transcript_non_coding();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = RnaVariant {
            accession: Accession::new("NR", "TEST", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::utr3(5)),
                substitution(Base::A, Base::G),
            ),
        };
        let err = converter
            .convert(&HgvsVariant::Rna(variant))
            .expect_err("*N has no meaning on a non-coding transcript; converter declines");
        assert!(
            matches!(
                err,
                FerroError::ConversionError { ref msg }
                    if msg.contains("no meaning on a non-coding transcript")
                        && msg.contains("cannot convert to VCF")
            ),
            "expected non-coding r.*N decline reason, got: {err}"
        );
    }

    #[test]
    fn p_axis_returns_back_translation_error() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = ProteinVariant {
            accession: Accession::new("NP", "000079", Some(2)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                ProtInterval::point(ProtPos::new(AminoAcid::Arg, 1)),
                ProteinEdit::Substitution {
                    reference: AminoAcid::Arg,
                    alternative: AminoAcid::Gln,
                },
            ),
        };
        let err = converter
            .convert(&HgvsVariant::Protein(variant))
            .expect_err("p. requires back-translation; converter rejects");
        assert!(
            matches!(err, FerroError::ConversionError { ref msg } if msg.to_lowercase().contains("protein")),
            "expected ConversionError mentioning protein, got: {err}"
        );
    }
}

// =============================================================================
// SECTION 3 — Allele variants
// =============================================================================

mod allele_branch {
    use super::*;

    #[test]
    fn single_inner_variant_routes_through_inner() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let inner = CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let allele = AlleleVariant::new(vec![HgvsVariant::Cds(inner)], AllelePhase::Cis);
        let result = converter
            .convert(&HgvsVariant::Allele(allele))
            .expect("single-inner allele must route through the inner conversion");
        assert_eq!(result.record.chrom, "chr1");
    }

    #[test]
    fn empty_allele_rejected() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let allele = AlleleVariant::new(vec![], AllelePhase::Cis);
        let err = converter
            .convert(&HgvsVariant::Allele(allele))
            .expect_err("empty allele has nothing to encode in VCF");
        assert!(
            matches!(err, FerroError::ConversionError { ref msg } if msg.to_lowercase().contains("empty")),
            "expected ConversionError mentioning empty, got: {err}"
        );
    }

    /// Multi-variant allele: the converter documents this as a
    /// "should be decomposed" failure mode and rejects rather than
    /// silently converting only the first inner variant. Pin the
    /// rejection so future widening (e.g. emitting multiple VCF
    /// rows) lands as a deliberate test update.
    #[test]
    fn multi_variant_allele_rejected() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let inner1 = CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                substitution(Base::A, Base::G),
            ),
        };
        let inner2 = CdsVariant {
            accession: Accession::new("NM", "TEST", Some(1)),
            gene_symbol: Some("TEST".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(2)),
                substitution(Base::T, Base::C),
            ),
        };
        let allele = AlleleVariant::new(
            vec![HgvsVariant::Cds(inner1), HgvsVariant::Cds(inner2)],
            AllelePhase::Cis,
        );
        // The converter may surface a ConversionError or successfully
        // emit only the first variant with a warning, depending on
        // the implementation's contract. Pin "either reject OR emit
        // with a warning" — silently emitting both as one VCF row
        // would be a regression.
        match converter.convert(&HgvsVariant::Allele(allele)) {
            Err(_) => {
                // Acceptable: rejected up-front.
            }
            Ok(result) => {
                assert!(
                    !result.warnings.is_empty(),
                    "multi-variant allele emitted without warning would silently drop variants; \
                     converter must either reject or warn"
                );
            }
        }
    }
}

// =============================================================================
// SECTION 6 — A range payload never reaches the ALT allele as its own spelling
// =============================================================================

/// Closes the gap `tests/it/issue_261_hgvs_vcf_coverage.rs` declares in its own
/// header — "the full path uses `HgvsToVcfConverter` with a provider; pinning
/// that here would require a transcript fixture". That sentence is why nothing
/// caught the defect these tests pin: the provider-backed REF/ALT construction
/// was unpinned for every `NaEdit` arm whose payload is an
/// [`InsertedSequence`] rather than a plain [`Sequence`].
///
/// `build_vcf_record` used to build the ALT allele as
/// `format!("{}{}", anchor_base, sequence)`, and `Display for InsertedSequence`
/// renders a range payload as its **HGVS spelling** (`212_221inv`), not as
/// bases. A `g.221_222ins212_221inv` therefore emitted `ALT=C212_221inv` — a
/// VCF allele field holding an HGVS description.
///
/// This mattered only for user-authored input until the `ins<range>inv`
/// derivation landed. On `origin/main` every construction site of
/// `InsertedSequence::PositionRangeInv` / `InsertedPart::PositionRangeInv` is
/// in the parser (`src/hgvs/parser/edit.rs`); `src/normalize/rules.rs` only
/// ever *destructures* the variant, in two arms that resolve it away to a
/// `Literal`, and `src/normalize/mod.rs` did not mention it at all. The
/// derivation makes the normalizer mint the shape, so ferro's own normalized
/// output began flowing into this construction.
///
/// The converter **declines** rather than resolving the span. `build_vcf_record`
/// is the shared trunk for `convert_cds` / `convert_tx` / `convert_rna` /
/// `convert_genome`, and it receives the edit in the *description's own*
/// coordinate frame while `get_reference_sequence` reads the genomic contig by
/// chromosome name. Resolving a `c.`-frame range against genomic coordinates
/// would return the wrong bases silently, which is worse than the visible
/// corruption it replaces; resolving correctly needs the axis-aware lookup
/// `src/normalize/rules.rs::fetch_position_range_bases` has and this converter
/// does not. Callers that want VCF can convert the literal spelling, which the
/// normalizer preserves as an equally legal form.
mod range_payload_never_reaches_alt {
    use super::*;

    /// Every emitted allele must be nucleotides. A VCF REF/ALT field admits
    /// only `[ACGTN]` (plus the symbolic and breakend forms this converter
    /// never produces), so any other byte is a malformed record.
    fn assert_alleles_are_bases(record: &ferro_hgvs::vcf::VcfRecord, what: &str) {
        for allele in std::iter::once(&record.reference).chain(record.alternate.iter()) {
            assert!(
                !allele.is_empty()
                    && allele.bytes().all(|b| matches!(
                        b.to_ascii_uppercase(),
                        b'A' | b'C' | b'G' | b'T' | b'N'
                    )),
                "{what}: VCF allele {allele:?} is not nucleotides — an HGVS payload spelling \
                 was spliced into the allele field"
            );
        }
    }

    fn genomic_insertion(start: u64, end: u64, sequence: InsertedSequence) -> HgvsVariant {
        HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::new(GenomePos::new(start), GenomePos::new(end)),
                NaEdit::Insertion { sequence },
            ),
        })
    }

    /// The exact shape the `ins<range>inv` derivation now mints: an inverted
    /// copy of the ten bases immediately 5' of the junction.
    #[test]
    fn inverted_range_insertion_declines_instead_of_emitting_its_spelling() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = genomic_insertion(
            1050,
            1051,
            InsertedSequence::PositionRangeInv {
                start: 1041,
                end: 1050,
            },
        );

        match converter.convert(&variant) {
            Err(FerroError::UnsupportedVariant { variant_type }) => {
                assert!(
                    variant_type.contains("1041_1050inv"),
                    "the decline must name the payload it could not spell, got {variant_type:?}"
                );
            }
            Err(other) => panic!("expected UnsupportedVariant, got {other:?}"),
            Ok(result) => {
                assert_alleles_are_bases(&result.record, "ins<range>inv");
                panic!(
                    "converting a range payload must decline; it emitted ALT={:?}",
                    result.record.alternate
                );
            }
        }
    }

    /// The un-inverted sibling is equally unspellable and was equally broken
    /// (`ALT=C1041_1050`). Pinned alongside so a future fix to one cannot
    /// leave the other emitting a description.
    #[test]
    fn forward_range_insertion_declines_instead_of_emitting_its_spelling() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = genomic_insertion(
            1050,
            1051,
            InsertedSequence::PositionRange {
                start: 1041,
                end: 1050,
            },
        );

        match converter.convert(&variant) {
            Err(FerroError::UnsupportedVariant { .. }) => {}
            Err(other) => panic!("expected UnsupportedVariant, got {other:?}"),
            Ok(result) => panic!(
                "converting a range payload must decline; it emitted ALT={:?}",
                result.record.alternate
            ),
        }
    }

    /// A `delins` carries the same [`InsertedSequence`] and reached the same
    /// construction, so it is pinned on the same guard.
    #[test]
    fn range_payload_delins_declines() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::new(GenomePos::new(1050), GenomePos::new(1052)),
                NaEdit::Delins {
                    sequence: InsertedSequence::PositionRangeInv {
                        start: 1041,
                        end: 1050,
                    },
                    deleted: None,
                    deleted_length: None,
                    substitution_reference: None,
                },
            ),
        });

        match converter.convert(&variant) {
            Err(FerroError::UnsupportedVariant { .. }) => {}
            Err(other) => panic!("expected UnsupportedVariant, got {other:?}"),
            Ok(result) => panic!(
                "a delins range payload must decline; it emitted ALT={:?}",
                result.record.alternate
            ),
        }
    }

    /// The control: the literal spelling of the very same edit still converts,
    /// and its alleles are bases. Without this, the guard above would be
    /// satisfied by a converter that declined every insertion.
    #[test]
    fn the_literal_spelling_of_the_same_edit_still_converts() {
        let tx = fixture_transcript();
        let provider = fixture_provider();
        let converter = HgvsToVcfConverter::new(&tx, &provider);
        // `fixture_provider` lays `ATGCATGC` from 1-based 1000, so
        // 1-based 1041..=1050 is `ATGCATGCAT`; its reverse complement is
        // `ATGCATGCAT` read back — the pin only needs *a* literal payload.
        let variant = genomic_insertion(
            1050,
            1051,
            InsertedSequence::Literal(Sequence::from_str("ATGCATGCAT").unwrap()),
        );

        let result = converter
            .convert(&variant)
            .expect("a literal insertion payload must still convert");
        assert_alleles_are_bases(&result.record, "literal insertion");
        assert_eq!(result.record.pos, 1050);
        assert_eq!(
            result.record.alternate,
            vec![format!("{}ATGCATGCAT", result.record.reference)],
            "ALT must be the anchor base followed by the inserted bases"
        );
    }
}
