//! Normalizer-independent same/different oracle. Reduces both outputs to
//! reference-anchored SPDI tuples (via the verified-independent
//! `spdi::convert::hgvs_to_spdi`), applies each to a shared reference window,
//! and compares the edited sequences. Shares NO code with `src/normalize/`.
//! See design spec §5.1.1 / §5.1.1a.

use crate::error::FerroError;
use crate::hgvs::interval::{CdsInterval, TxInterval};
use crate::hgvs::location::{CdsPos, TxPos};
use crate::hgvs::variant::HgvsVariant;
use crate::project::VariantProjector;
use crate::reference::provider::ReferenceProvider;
use crate::spdi::convert::hgvs_to_spdi;

/// A reference-anchored edit: SPDI-equivalent (accession, 0-based interbase
/// position, deleted bases, inserted bases).
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct EditTuple {
    pub accession: String,
    pub position: u64,
    pub deleted: String,
    pub inserted: String,
}

/// Reduce a variant to an `EditTuple` via `hgvs_to_spdi`. Normalizer-free.
pub(crate) fn edit_tuple(
    v: &HgvsVariant,
    provider: &(impl ReferenceProvider + ?Sized),
) -> Result<EditTuple, FerroError> {
    let spdi = hgvs_to_spdi(v, provider)?;
    Ok(EditTuple {
        accession: spdi.sequence.clone(),
        position: spdi.position,
        deleted: spdi.deletion.clone(),
        inserted: spdi.insertion.clone(),
    })
}

/// Apply an edit tuple to a reference window and return the edited sequence.
/// `window_start` is the 0-based coordinate of `window[0]`, in the SAME basis
/// as `t.position` (SPDI 0-based interbase). Pure string operation.
pub(crate) fn apply_spdi_to_window(window: &str, window_start: u64, t: &EditTuple) -> String {
    let rel = (t.position - window_start) as usize;
    let bytes = window.as_bytes();
    let del_len = t.deleted.len();
    // Invariant guaranteed by the sole caller (`compare_variants`): the window
    // is fetched to fully contain both edit footprints, so the splice indices
    // are always in bounds. `hgvs_to_spdi` validates deletion bounds against
    // the reference first, so a malformed tuple never reaches here.
    debug_assert!(
        rel + del_len <= bytes.len(),
        "edit footprint {}..{} exceeds window length {}",
        rel,
        rel + del_len,
        bytes.len()
    );
    let mut out = String::with_capacity(window.len() + t.inserted.len());
    out.push_str(std::str::from_utf8(&bytes[..rel]).unwrap());
    out.push_str(&t.inserted);
    out.push_str(std::str::from_utf8(&bytes[rel + del_len..]).unwrap());
    out
}

/// Result of the normalizer-independent same/different comparison.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SameVariant {
    /// Both variants edit the shared reference window to the same sequence.
    Same,
    /// Both reduce to SPDI on the same accession but edit it differently.
    Different,
    /// The two variants reduce to SPDI on *different* accessions, so there is
    /// no shared reference window to compare on — the oracle declines.
    BasisMismatch {
        /// The accession ferro's output reduced to.
        ferro: String,
        /// The accession the other tool's output reduced to.
        other: String,
    },
}

/// Padding added on each side of the edit footprint so the comparison window
/// contains the full ambiguous/rolled region for either representation.
const WINDOW_PAD: u64 = 32;

/// Compare two variants by edited-sequence over a shared reference window,
/// independent of `src/normalize/`.
///
/// Each variant is reduced to a reference-anchored SPDI tuple ([`edit_tuple`],
/// via `hgvs_to_spdi`), and both tuples are applied to the *same* window of
/// reference sequence ([`apply_spdi_to_window`]); the verdict is
/// [`SameVariant::Same`] iff the two edited sequences are byte-identical. This
/// makes the comparison roll-invariant: two shift-equivalent spellings of the
/// same edit (e.g. a duplication written at different positions within a
/// repeat run) produce different SPDI tuples but the same edited sequence, and
/// so compare as `Same` — the case a naive SPDI-tuple-equality check gets
/// wrong.
///
/// When the two tuples land on different accessions there is no shared window
/// to compare on, so the function returns [`SameVariant::BasisMismatch`] rather
/// than an error. Both variants must already be in the same frame (frame
/// resolution is a separate step); here they are compared on their SPDI
/// accession directly.
///
/// # Errors
///
/// Returns [`FerroError`] if either variant cannot be reduced to an SPDI tuple
/// (e.g. an unsupported edit or an accession the provider lacks), or if the
/// reference window cannot be fetched.
pub fn compare_variants(
    ferro: &HgvsVariant,
    other: &HgvsVariant,
    provider: &(impl ReferenceProvider + ?Sized),
) -> Result<SameVariant, FerroError> {
    let tf = edit_tuple(ferro, provider)?;
    let to = edit_tuple(other, provider)?;

    if tf.accession != to.accession {
        return Ok(SameVariant::BasisMismatch {
            ferro: tf.accession,
            other: to.accession,
        });
    }

    // Pad the window by WINDOW_PAD plus the widest edit footprint so the rolled
    // region of either representation is fully contained on both sides.
    let widest = [
        tf.deleted.len(),
        tf.inserted.len(),
        to.deleted.len(),
        to.inserted.len(),
    ]
    .into_iter()
    .max()
    .unwrap_or(0) as u64;
    let pad = WINDOW_PAD + widest;

    let lo = tf.position.min(to.position).saturating_sub(pad);
    let hi =
        (tf.position + tf.deleted.len() as u64).max(to.position + to.deleted.len() as u64) + pad;
    // Clamp the high bound to the sequence length: `get_sequence` errors if
    // `end` runs past the end of the reference. `lo` is already floored at 0 by
    // the saturating subtraction above.
    let hi = hi.min(provider.get_sequence_length(&tf.accession)?);

    // Fetch the shared reference window on the SPDI accession, in the SAME
    // 0-based basis as SPDI positions (SPDI position 0 == `get_sequence` start
    // 0), so `apply_spdi_to_window` splices with no coordinate adjustment.
    let window = provider.get_sequence(&tf.accession, lo, hi)?;

    let obs_ferro = apply_spdi_to_window(&window, lo, &tf);
    let obs_other = apply_spdi_to_window(&window, lo, &to);
    Ok(if obs_ferro == obs_other {
        SameVariant::Same
    } else {
        SameVariant::Different
    })
}

/// Return `v` in a frame the oracle can apply-to-reference.
///
/// Exonic transcript (`c.`/`n.`) and genomic variants pass through
/// unchanged. Intronic `c.`/`n.` variants are projected onto the chromosome
/// (`NC_`) frame via [`VariantProjector::project_to_nc_frame`], because
/// intronic bases are absent from the mature transcript and
/// [`hgvs_to_spdi`] hard-errors on them (§5.1.1a). The oracle only needs the
/// `NC_` pivot, not a re-anchored `NG_`/`LRG_` parent frame — it compares
/// edited sequences over a window (roll-invariant), so a raw, non-parent
/// pivot suffices and works even on a bare `NM_:c.` input that carries no
/// explicit `genomic_context`.
///
/// # Errors
///
/// Propagates [`VariantProjector::project_to_nc_frame`]'s error contract when
/// an intronic variant cannot be pivoted to the `NC_` frame (e.g. the
/// transcript is absent from cdot).
pub fn resolve_frame<P: ReferenceProvider + Clone>(
    v: &HgvsVariant,
    projector: &VariantProjector<P>,
) -> Result<HgvsVariant, FerroError> {
    if variant_is_intronic(v) {
        projector.project_to_nc_frame(v)
    } else {
        Ok(v.clone())
    }
}

/// True if either endpoint of a `c.`/`n.` variant carries a nonzero intronic
/// offset (`CdsPos::is_intronic` / `TxPos::is_intronic`, the same predicate
/// `resolve_cds_to_tx` in `src/spdi/convert.rs` uses to reject intronic input).
/// Genomic and other variant kinds have no offset notation and are never
/// intronic by this definition.
fn variant_is_intronic(v: &HgvsVariant) -> bool {
    match v {
        HgvsVariant::Cds(c) => cds_location_is_intronic(&c.loc_edit.location),
        HgvsVariant::Tx(n) => tx_location_is_intronic(&n.loc_edit.location),
        _ => false,
    }
}

/// True if either resolved endpoint of a CDS interval is intronic. A compound
/// boundary (`UncertainBoundary::Range`, e.g. `(4185+1_4186-1)`) has no single
/// resolved position — `UncertainBoundary::inner` returns `None` for it — so
/// it is conservatively treated as non-intronic here; only bare/parenthesized
/// single positions carry a `CdsPos` to test.
fn cds_location_is_intronic(interval: &CdsInterval) -> bool {
    interval.start.inner().is_some_and(CdsPos::is_intronic)
        || interval.end.inner().is_some_and(CdsPos::is_intronic)
}

/// Transcript-interval (`n.`) sibling of [`cds_location_is_intronic`].
fn tx_location_is_intronic(interval: &TxInterval) -> bool {
    interval.start.inner().is_some_and(TxPos::is_intronic)
        || interval.end.inner().is_some_and(TxPos::is_intronic)
}

/// Test-only fixture support for [`resolve_frame`], factored into its own
/// module so `frame_tests` can build a real intronic `VariantProjector`
/// without duplicating it per test.
#[cfg(test)]
pub(crate) mod frame_tests_support {
    use crate::data::cdot::{CdotMapper, CdotTranscript};
    use crate::data::projection::Projector;
    use crate::project::VariantProjector;
    use crate::reference::mock::MockProvider;
    use crate::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
    use crate::reference::Strand as ProvStrand;
    use std::sync::OnceLock;

    /// Two-exon coding transcript fixture (`NM_INTR.1`) with a real intron
    /// between the exons, so an intronic `c.` offset has genomic bases to
    /// project onto. Mirrors (a deliberate small duplicate of, since the
    /// original is `#[cfg(test)]`-private to its own module)
    /// `project::projector`'s `make_intronic_test_data`:
    ///
    ///   Exon 1: genome [1000, 1010), tx [0, 10)
    ///   Intron: genome [1010, 2000)
    ///   Exon 2: genome [2000, 2010), tx [10, 20)
    pub(crate) fn intronic_fixture() -> (VariantProjector<MockProvider>, MockProvider) {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_INTR.1".to_string(),
            CdotTranscript {
                cds_start_incomplete: false,
                gene_name: Some("INTRGENE".to_string()),
                // A real versioned `NC_` chromosome accession (not a `chr`-style
                // alias), so the arbitration NC pivot in `nc_pivot_context` can
                // anchor bare-`NM_`/`NG_`/`LRG_`-parent inputs to it.
                contig: "NC_000001.11".to_string(),
                strand: ProvStrand::Plus,
                exons: vec![[1000, 1010, 0, 10], [2000, 2010, 10, 20]],
                cds_start: Some(0),
                cds_end: Some(18),
                gene_id: None,
                protein: Some("NP_INTR.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        let projector = Projector::new(cdot);

        let mut provider = MockProvider::new();
        provider.add_transcript(Transcript {
            cds_start_incomplete: false,
            id: "NM_INTR.1".to_string(),
            gene_symbol: Some("INTRGENE".to_string()),
            strand: TxStrand::Plus,
            sequence: Some("ATGCGCAAAGGGTAACCC".to_string()),
            cds_start: Some(1),
            cds_end: Some(18),
            exons: vec![Exon::new(1, 1, 10), Exon::new(2, 11, 20)],
            // Keep chromosome, cdot `contig` (NC_000001.11 above), and the
            // add_genomic_sequence key below all aligned to the same versioned
            // NC_ accession so an apply-to-reference through this projector
            // resolves the sequence (no chr-alias vs NC_ split).
            chromosome: Some("NC_000001.11".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2009),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });
        let exon1 = "ATGCGCAAAG";
        let intron = "N".repeat(990);
        let exon2 = "GGTAACCCNN";
        let prefix = "N".repeat(999);
        let suffix = "N".repeat(100);
        provider.add_genomic_sequence(
            "NC_000001.11",
            format!("{prefix}{exon1}{intron}{exon2}{suffix}"),
        );

        let vp = VariantProjector::new(projector, provider.clone());
        (vp, provider)
    }
}

#[cfg(test)]
mod frame_tests {
    use super::*;
    use crate::hgvs::edit::NaEdit;
    use crate::hgvs::interval::CdsInterval;
    use crate::hgvs::location::CdsPos;
    use crate::hgvs::parser::parse_hgvs;
    use crate::hgvs::variant::{CdsVariant, LocEdit};
    use crate::project::accession::parse_accession;
    use frame_tests_support::intronic_fixture;

    // -- variant_is_intronic: pure, no projector needed ---------------------

    #[test]
    fn intronic_offset_on_start_or_end_is_detected() {
        let upstream_offset = parse_hgvs("NM_TEST.1:c.68-7A>G").unwrap();
        assert!(
            variant_is_intronic(&upstream_offset),
            "c.68-7 carries a nonzero offset"
        );

        let downstream_offset = parse_hgvs("NM_TEST.1:c.4+7A>G").unwrap();
        assert!(
            variant_is_intronic(&downstream_offset),
            "c.4+7 carries a nonzero offset"
        );
    }

    #[test]
    fn exonic_cds_variant_is_not_intronic() {
        let exonic = parse_hgvs("NM_TEST.1:c.100A>G").unwrap();
        assert!(!variant_is_intronic(&exonic), "c.100 has no offset");
    }

    #[test]
    fn genomic_variant_is_never_intronic() {
        let genomic = parse_hgvs("NC_000001.11:g.100A>G").unwrap();
        assert!(
            !variant_is_intronic(&genomic),
            "g. variants carry no offset notation"
        );
    }

    #[test]
    fn asymmetric_range_with_one_intronic_endpoint_is_intronic() {
        // c.10_11+5del: an exonic start (c.10, no offset) and an intronic
        // end (c.11+5, offset +5) — a genuinely asymmetric range. Exercises
        // the `start.is_intronic() || end.is_intronic()` OR-branch on two
        // different endpoints, rather than the single-point variants above
        // where start and end are the same position.
        let asymmetric = parse_hgvs("NM_TEST.1:c.10_11+5del").unwrap();
        assert!(
            variant_is_intronic(&asymmetric),
            "c.10_11+5del has an intronic end even though its start is exonic"
        );
    }

    // -- resolve_frame: exonic passthrough (REQUIRED) ------------------------

    #[test]
    fn resolve_frame_passes_exonic_variant_through_unchanged() {
        let (projector, _provider) = intronic_fixture();
        // Exonic passthrough never touches the projector, so it needs no
        // genomic_context on the accession — only the intronic branch does.
        let exonic = parse_hgvs("NM_INTR.1:c.4A>G").unwrap();
        let resolved = resolve_frame(&exonic, &projector).unwrap();
        assert_eq!(resolved.to_string(), exonic.to_string());
    }

    // -- resolve_frame: intronic -> genomic (PREFERRED, not required) -------
    //
    // Covered end-to-end here rather than left to delegation, since the
    // two-exon `intronic_fixture` above made it cheap: it reuses
    // `project::projector`'s own `make_intronic_test_data` layout (plus
    // strand, exon1 genome [1000,1010)/tx[0,10), intron genome [1010,2000),
    // exon2 genome [2000,2010)/tx[10,20)), which that module's
    // `project_to_genomic_intronic_offset_plus` test already exercises with
    // an identical `c.10+5` input, projecting to `g.1014`.

    #[test]
    fn resolve_frame_projects_intronic_variant_to_genomic() {
        let (projector, _provider) = intronic_fixture();

        // c.10+5: 5 bases into the intron past the exon-1/intron boundary
        // (exon 1 is tx [0,10) / genome [1000,1010)) -> g.1014 on NC_000001.11.
        let cds = CdsVariant {
            accession: parse_accession("NM_INTR.1")
                .with_genomic_context(parse_accession("NC_000001.11")),
            gene_symbol: Some("INTRGENE".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::with_offset(10, 5)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        let intronic = HgvsVariant::Cds(cds);
        assert!(
            variant_is_intronic(&intronic),
            "c.10+5 must be routed through the intronic branch"
        );

        let resolved = resolve_frame(&intronic, &projector).unwrap();
        match resolved {
            HgvsVariant::Genome(g) => {
                assert_eq!(g.accession.to_string(), "NC_000001.11");
                let start = g
                    .loc_edit
                    .location
                    .start
                    .inner()
                    .expect("resolved genomic start should be concrete");
                assert_eq!(start.base, 1014, "expected g.1014 for c.10+5");
            }
            other => panic!("expected a Genome variant, got {other:?}"),
        }
    }

    // -- resolve_frame: intronic under a non-NC_ parent (REGRESSION #7e9220f) --
    //
    // The bare-`NM_` intronic fix (`project_to_nc_frame`) originally only
    // *filled in* a missing `genomic_context` from cdot's contig, leaving an
    // already-present `NG_`/`LRG_` parent untouched — so `project_to_genomic_nc`
    // stamped the short `NG_`/`LRG_` accession onto NC-scale chromosome
    // coordinates, overflowed the parent sequence, and dropped the whole
    // arbitration to `Verdict::Inconclusive`. The pivot must instead REPLACE the
    // parent with the transcript's true chromosome (`NC_`) accession for *every*
    // parent class. The existing `resolve_frame_projects_intronic_variant_to_genomic`
    // above uses an `NC_` parent (the one class that was never broken) and so
    // cannot catch this; these two do.

    /// Shared assertion: `c.10+5del` under `context_parent` must pivot to the
    /// chromosome (`NC_000001.11`) frame at `g.1014`, never staying on the
    /// non-`NC_` parent accession.
    fn assert_pivots_to_nc_g1014(context_parent: &str) {
        let (projector, _provider) = intronic_fixture();
        let cds = CdsVariant {
            accession: parse_accession("NM_INTR.1")
                .with_genomic_context(parse_accession(context_parent)),
            gene_symbol: Some("INTRGENE".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::with_offset(10, 5)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        let intronic = HgvsVariant::Cds(cds);
        assert!(
            variant_is_intronic(&intronic),
            "c.10+5 must route through the intronic branch"
        );
        match resolve_frame(&intronic, &projector).unwrap() {
            HgvsVariant::Genome(g) => {
                assert_eq!(
                    g.accession.to_string(),
                    "NC_000001.11",
                    "pivot must re-anchor a {context_parent} parent to the chromosome (NC_) \
                     accession, not stamp NC coordinates under {context_parent}"
                );
                let start = g
                    .loc_edit
                    .location
                    .start
                    .inner()
                    .expect("resolved genomic start should be concrete");
                assert_eq!(start.base, 1014, "expected g.1014 for c.10+5");
            }
            other => panic!("expected a Genome variant, got {other:?}"),
        }
    }

    #[test]
    fn resolve_frame_projects_ng_parent_intronic_variant_to_nc() {
        // `NG_007107.2(NM_004992.3):c.378-17del`-shaped input (the notation
        // Mutalyzer emits). The `NG_` parent must be replaced by the NC frame.
        assert_pivots_to_nc_g1014("NG_048011.1");
    }

    #[test]
    fn resolve_frame_projects_lrg_parent_intronic_variant_to_nc() {
        // `LRG_<n>(NM_...):c.…` — an `LRG_` genomic parent must likewise pivot
        // onto the chromosome (`NC_`) frame rather than the short `LRG_` record.
        assert_pivots_to_nc_g1014("LRG_1");
    }

    #[test]
    fn resolve_frame_projects_bare_nm_intronic_variant_to_nc() {
        // A bare `NM_:c.` intronic input (no `genomic_context`) still pivots to
        // the chromosome frame via cdot's contig — the case #7e9220f first fixed.
        let (projector, _provider) = intronic_fixture();
        let cds = CdsVariant {
            accession: parse_accession("NM_INTR.1"),
            gene_symbol: Some("INTRGENE".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::with_offset(10, 5)),
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        match resolve_frame(&HgvsVariant::Cds(cds), &projector).unwrap() {
            HgvsVariant::Genome(g) => {
                assert_eq!(g.accession.to_string(), "NC_000001.11");
                assert_eq!(
                    g.loc_edit.location.start.inner().unwrap().base,
                    1014,
                    "expected g.1014 for c.10+5"
                );
            }
            other => panic!("expected a Genome variant, got {other:?}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn apply_spdi_splices_deletion() {
        // window = "AACGTT", 0-based window_start = 10.
        // delete 1 base "C" at position 12 (0-based interbase before index 2).
        let t = EditTuple {
            accession: "X".into(),
            position: 12,
            deleted: "C".into(),
            inserted: String::new(),
        };
        let obs = apply_spdi_to_window("AACGTT", 10, &t);
        assert_eq!(obs, "AAGTT");
    }

    #[test]
    fn apply_spdi_splices_insertion() {
        // insert "TT" at interbase position 12 (after 2 bases of the window).
        let t = EditTuple {
            accession: "X".into(),
            position: 12,
            deleted: String::new(),
            inserted: "TT".into(),
        };
        let obs = apply_spdi_to_window("AACGTT", 10, &t);
        assert_eq!(obs, "AATTCGTT");
    }
}

#[cfg(test)]
mod verdict_tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;
    use crate::reference::mock::MockProvider;

    /// A 105 bp synthetic genomic contig with an A-free left flank
    /// (1-based 1..=50), a 5 bp poly-A run (1-based 51..=55), and an
    /// A-free right flank (1-based 56..=105). The poly-A run lets a
    /// single-base duplication be spelled at several positions inside the
    /// run while denoting the same variant (insert one A → a 6-A run).
    ///
    /// Returns the provider and the exact contig sequence (for auditing).
    fn provider_with_a_run() -> (MockProvider, String) {
        let left: String = "CGT".chars().cycle().take(50).collect();
        let right: String = "GTC".chars().cycle().take(50).collect();
        let seq = format!("{left}AAAAA{right}");
        assert_eq!(seq.len(), 105);
        // 1-based 51..=55 (0-based 50..55) is the poly-A run.
        assert_eq!(&seq[50..55], "AAAAA");
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence("NC_000001.11", &seq);
        (provider, seq)
    }

    #[test]
    fn shift_equivalent_dups_are_same_variant() {
        // Two spellings of "insert one A into the poly-A run" at different
        // positions inside the run. `g.51dup` duplicates the 5'-most A;
        // `g.55dup` duplicates the 3'-most A. Both denote a 6-A run, so the
        // oracle must call them Same — this is the roll-invariance property
        // that a naive SPDI-tuple-equality check would get wrong (the two
        // SPDI tuples are `NC_000001.11:51::A` and `NC_000001.11:55::A`,
        // which are unequal as tuples but equal once applied to reference).
        let (provider, _seq) = provider_with_a_run();
        let ferro = parse_hgvs("NC_000001.11:g.51dup").unwrap();
        let other = parse_hgvs("NC_000001.11:g.55dup").unwrap();
        let verdict = compare_variants(&ferro, &other, &provider).unwrap();
        assert_eq!(verdict, SameVariant::Same);
    }

    #[test]
    fn genuinely_different_edits_are_different() {
        // A dup in the poly-A run vs a substitution in the left flank: these
        // edit the reference to genuinely different sequences.
        let (provider, seq) = provider_with_a_run();
        // 0-based index 9 == 1-based 10; confirm the reference base is 'C'.
        assert_eq!(&seq[9..10], "C");
        let ferro = parse_hgvs("NC_000001.11:g.51dup").unwrap();
        let other = parse_hgvs("NC_000001.11:g.10C>T").unwrap();
        let verdict = compare_variants(&ferro, &other, &provider).unwrap();
        assert_eq!(verdict, SameVariant::Different);
    }

    #[test]
    fn different_accession_is_basis_mismatch() {
        // Two variants on two distinct registered contigs reduce to SPDI
        // tuples with different accessions, so the oracle declines to
        // compare edited sequences and reports BasisMismatch.
        let (mut provider, _seq) = provider_with_a_run();
        provider.add_genomic_sequence("NC_000002.11", "T".repeat(105));
        let ferro = parse_hgvs("NC_000001.11:g.51dup").unwrap();
        let other = parse_hgvs("NC_000002.11:g.51dup").unwrap();
        let verdict = compare_variants(&ferro, &other, &provider).unwrap();
        assert_eq!(
            verdict,
            SameVariant::BasisMismatch {
                ferro: "NC_000001.11".to_string(),
                other: "NC_000002.11".to_string(),
            }
        );
    }

    #[test]
    fn oracle_module_does_not_import_the_normalizer() {
        // Structural honesty guard (design spec §5.1.1): the oracle source
        // must not IMPORT the normalizer module. Build the forbidden import
        // needle at runtime so this assertion's own source does not contain
        // the literal (otherwise the scan matches itself and always fails).
        let src = include_str!("oracle.rs");
        let forbidden_import = format!("use crate::{}", "normalize");
        let forbidden_shuffle = format!("{}::shuffle", "normalize");
        assert!(
            !src.contains(&forbidden_import) && !src.contains(&forbidden_shuffle),
            "oracle.rs must not import src/normalize (honesty property)"
        );
    }
}
