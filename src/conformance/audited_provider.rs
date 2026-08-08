//! A [`ReferenceProvider`] wrapper that records every reference read it
//! **fails**, and counts every one it serves.
//!
//! # Why a hermetic gate needs this
//!
//! It exists because of a measured false negative, not a hypothetical one. A
//! window fixture captured from a normalize pass is exactly as wide as the reads
//! that pass happened to make, so the fixture serves the recorded pass and
//! errors on any pass that reads further. `src/normalize/mod.rs` then turns that
//! error into a *success*:
//!
//! ```text
//! let ref_seq = match seq_result {
//!     Ok(s) => s,
//!     Err(_) => return Ok((HV::Genome(self.canonicalize_genome_variant(variant)), vec![])),
//! };
//! ```
//!
//! — the input is returned unchanged, with an empty warning vector. For a
//! lone-member variant that is byte-identical to a correct answer, so a guard
//! built on such a fixture is **green while the defect is fully present**. That
//! has already happened in this repository: a guard passed on three of its five
//! rows with the defect in full.
//!
//! An under-serving slice is therefore indistinguishable, at the output, from a
//! clean normalization. The remedy is to stop looking only at the output: assert
//! that **no reference read failed** while producing it, and assert it *before*
//! any verdict is computed, so "the fixture could not answer" always outranks
//! "ferro answered differently".
//!
//! # Why a successful-read floor as well as a failure ceiling
//!
//! Zero failures is trivially satisfiable by reading nothing. A refactor that
//! stops consulting the reference — or a corpus edit that leaves every row
//! short-circuiting before the first read — would satisfy the ceiling and prove
//! nothing. [`AuditedProvider::successful_reads`] is the other half: a caller
//! asserts a positive floor alongside the zero ceiling, so the audit reports on
//! a pass that actually happened.
//!
//! # What is audited, and why it is not simply "the `Result`-returning reads"
//!
//! Audited: the reads that source reference **bases** — [`get_transcript`],
//! [`get_transcript_for_variant`], [`get_transcript_for_accession`],
//! [`get_sequence`] and [`get_genomic_sequence`]. Those are exactly the accesses
//! `examples/common/recording.rs` records and a
//! [`WindowFixture`](super::reference_window::WindowFixture) stores, so a
//! failure among them means a missing window and "regenerate the slice" is a
//! remedy that works.
//!
//! Everything else is delegated verbatim, including two `Result`-returning
//! methods, because auditing them would report a legitimate absence as a fixture
//! defect and send the reader to a regeneration that cannot fix it:
//!
//! * `get_protein_sequence` — a window fixture carries no protein data to lose,
//!   and `has_protein_data` is the honest signal a caller checks. Worse,
//!   [`ReferenceProvider::get_protein_length`]'s default implementation
//!   binary-searches the length by probing `get_protein_sequence` and reading
//!   `Err` as "past the end", so a *successful* length query is built out of
//!   ~17 deliberate failures. Auditing that read would turn intended control
//!   flow into a wall of reported defects.
//! * `get_sequence_length` — it sources no bases, and the recorder does not
//!   record length lookups, so a fixture can never be regenerated to answer one
//!   it does not already cover. In-tree callers treat `Err` as "length unknown"
//!   and proceed, and any read that could actually starve normalization goes
//!   through `get_sequence` or `get_genomic_sequence`, which *are* audited.
//!
//! The `has_*` probes and the placement/selector lookups return no `Result` at
//! all and are delegated for the same reason: an existence check that answers
//! "no" is an answer, not a missing window.
//!
//! # Every method is forwarded, deliberately
//!
//! `ReferenceProvider` has many defaulted methods, and a wrapper that omits one
//! silently falls through to the trait default — measuring itself rather than
//! the provider it wraps. That is exactly the `ArcProvider` defect of #726
//! (`resolve_legacy_gene_selector` and `get_transcript_for_accession` were not
//! forwarded), and the trait's own maintainer note calls it out. So this type
//! forwards the **whole** surface, audited or not.
//!
//! [`get_transcript`]: ReferenceProvider::get_transcript
//! [`get_transcript_for_variant`]: ReferenceProvider::get_transcript_for_variant
//! [`get_transcript_for_accession`]: ReferenceProvider::get_transcript_for_accession
//! [`get_sequence`]: ReferenceProvider::get_sequence
//! [`get_genomic_sequence`]: ReferenceProvider::get_genomic_sequence

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use crate::hgvs::variant::Accession;
use crate::reference::provider::GenomicPlacement;
use crate::reference::transcript::Transcript;
use crate::{FerroError, HgvsVariant, ReferenceProvider};

/// The shared tally. Held behind an [`Arc`] so a cloned provider — which every
/// `Normalizer` construction takes — accumulates into one record.
#[derive(Debug, Default)]
struct Audit {
    failures: Mutex<Vec<String>>,
    successes: AtomicUsize,
}

/// Wraps any [`ReferenceProvider`], forwarding the whole trait surface and
/// recording the outcome of every base-sourcing read.
///
/// See the module documentation for what is audited and why the rest is
/// delegated verbatim.
#[derive(Debug)]
pub struct AuditedProvider<P> {
    inner: P,
    audit: Arc<Audit>,
}

// Derived `Clone` would demand `Audit: Clone` and, worse, would be *wrong*: the
// tally must be shared, not copied, or a clone handed to a `Normalizer` records
// into a tally nobody reads.
impl<P: Clone> Clone for AuditedProvider<P> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
            audit: Arc::clone(&self.audit),
        }
    }
}

impl<P> AuditedProvider<P> {
    /// Wrap `inner`, starting with an empty tally.
    pub fn new(inner: P) -> Self {
        Self {
            inner,
            audit: Arc::new(Audit::default()),
        }
    }

    /// The wrapped provider.
    pub fn inner(&self) -> &P {
        &self.inner
    }

    /// Every failed base-sourcing read, as `"<call>: <error>"`, in the order
    /// they happened.
    pub fn failures(&self) -> Vec<String> {
        self.audit.failures.lock().expect("audit lock").clone()
    }

    /// How many base-sourcing reads were **served**.
    ///
    /// Assert a positive floor on this alongside the zero-failure ceiling: an
    /// audit that saw no reads at all reports zero failures and proves nothing.
    pub fn successful_reads(&self) -> usize {
        self.audit.successes.load(Ordering::Relaxed)
    }

    /// Forget everything recorded so far, so one provider can audit several
    /// passes independently.
    pub fn reset(&self) {
        self.audit.failures.lock().expect("audit lock").clear();
        self.audit.successes.store(0, Ordering::Relaxed);
    }

    /// Tally one audited read, passing its outcome through untouched.
    fn record<T>(
        &self,
        what: impl FnOnce() -> String,
        outcome: Result<T, FerroError>,
    ) -> Result<T, FerroError> {
        match &outcome {
            Ok(_) => {
                self.audit.successes.fetch_add(1, Ordering::Relaxed);
            }
            Err(e) => self
                .audit
                .failures
                .lock()
                .expect("audit lock")
                .push(format!("{}: {e}", what())),
        }
        outcome
    }
}

impl<P: ReferenceProvider> ReferenceProvider for AuditedProvider<P> {
    // -- audited: the reads that source reference bases -----------------------

    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        let outcome = self.inner.get_transcript(id);
        self.record(|| format!("get_transcript({id})"), outcome)
    }

    fn get_transcript_for_variant(
        &self,
        variant: &HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
        let outcome = self.inner.get_transcript_for_variant(variant);
        self.record(|| format!("get_transcript_for_variant({variant})"), outcome)
    }

    fn get_transcript_for_accession(
        &self,
        accession: &Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
        let outcome = self.inner.get_transcript_for_accession(accession);
        self.record(
            || format!("get_transcript_for_accession({})", accession.full()),
            outcome,
        )
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        let outcome = self.inner.get_sequence(id, start, end);
        self.record(|| format!("get_sequence({id}, {start}, {end})"), outcome)
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        let outcome = self.inner.get_genomic_sequence(contig, start, end);
        self.record(
            || format!("get_genomic_sequence({contig}, {start}, {end})"),
            outcome,
        )
    }

    // -- delegated verbatim ---------------------------------------------------

    fn genomic_placement(&self, parent: &Accession) -> Option<GenomicPlacement> {
        self.inner.genomic_placement(parent)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        self.inner.genomic_placement_on_build(parent, build)
    }

    fn infer_genome_build(&self, accession: &Accession) -> Option<&'static str> {
        self.inner.infer_genome_build(accession)
    }

    fn resolve_legacy_gene_selector(
        &self,
        selector: &str,
        ng_parent: Option<&Accession>,
    ) -> Option<String> {
        self.inner.resolve_legacy_gene_selector(selector, ng_parent)
    }

    fn sole_hosted_transcript(&self, ng_parent: &Accession) -> Option<String> {
        self.inner.sole_hosted_transcript(ng_parent)
    }

    fn has_transcript(&self, id: &str) -> bool {
        self.inner.has_transcript(id)
    }

    fn has_transcript_version_exact(&self, id: &str) -> bool {
        self.inner.has_transcript_version_exact(id)
    }

    fn has_genomic_data(&self) -> bool {
        self.inner.has_genomic_data()
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.inner.get_protein_sequence(accession, start, end)
    }

    fn get_protein_length(&self, accession: &str) -> Result<u64, FerroError> {
        self.inner.get_protein_length(accession)
    }

    fn has_protein_data(&self) -> bool {
        self.inner.has_protein_data()
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        self.inner.get_sequence_length(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::mock::MockProvider;
    use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand};

    fn provider() -> AuditedProvider<MockProvider> {
        let mut inner = MockProvider::new();
        inner.add_transcript(Transcript::new(
            "NM_TEST.1".to_string(),
            Some("TEST".to_string()),
            Strand::Plus,
            "ACGTACGTACGT".to_string(),
            Some(1),
            Some(12),
            vec![Exon::new(1, 1, 12)],
            None,
            None,
            None,
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        ));
        AuditedProvider::new(inner)
    }

    #[test]
    fn a_served_read_is_counted_and_not_reported() {
        let audited = provider();
        assert_eq!(audited.get_sequence("NM_TEST.1", 0, 4).unwrap(), "ACGT");
        assert_eq!(audited.successful_reads(), 1);
        assert!(audited.failures().is_empty());
    }

    #[test]
    fn a_failed_read_is_reported_and_not_counted_as_served() {
        let audited = provider();
        assert!(audited.get_transcript("NM_ABSENT.1").is_err());
        assert_eq!(audited.successful_reads(), 0);
        let failures = audited.failures();
        assert_eq!(failures.len(), 1);
        assert!(
            failures[0].starts_with("get_transcript(NM_ABSENT.1): "),
            "a failure must name the call that made it: {:?}",
            failures[0]
        );
    }

    #[test]
    fn a_clone_shares_the_tally() {
        // The load-bearing property: `Normalizer::new` takes the provider by
        // value, so every pass runs against a clone. A per-clone tally would
        // report zero failures no matter what the pass did.
        let audited = provider();
        let clone = audited.clone();
        let _ = clone.get_sequence("NM_TEST.1", 0, 4);
        let _ = clone.get_sequence("NM_ABSENT.1", 0, 4);
        assert_eq!(audited.successful_reads(), 1);
        assert_eq!(audited.failures().len(), 1);
    }

    #[test]
    fn reset_clears_both_halves_of_the_tally() {
        let audited = provider();
        let _ = audited.get_sequence("NM_TEST.1", 0, 4);
        let _ = audited.get_sequence("NM_ABSENT.1", 0, 4);
        audited.reset();
        assert_eq!(audited.successful_reads(), 0);
        assert!(audited.failures().is_empty());
    }

    #[test]
    fn a_protein_length_probe_is_not_reported_as_a_missing_window() {
        // `get_protein_length`'s default binary-searches by probing
        // `get_protein_sequence` and reading `Err` as "past the end", so a
        // *successful* query is built out of many deliberate failures. Auditing
        // that read would bury a real fixture defect under ~17 false positives.
        let audited = provider();
        assert!(!audited.has_protein_data());
        assert_eq!(audited.get_protein_length("NP_TEST.1").ok(), Some(0));
        assert!(
            audited.failures().is_empty(),
            "protein probes must not be reported as failed reference reads: {:?}",
            audited.failures()
        );
    }

    #[test]
    fn an_absent_sequence_length_is_not_reported_as_a_missing_window() {
        // A length lookup sources no bases and the recorder never captures one,
        // so no regeneration could satisfy it; reporting it would send the
        // reader to a remedy that cannot work.
        let audited = provider();
        assert!(audited.get_sequence_length("NC_ABSENT.1").is_err());
        assert!(audited.failures().is_empty());
    }

    #[test]
    fn existence_probes_are_forwarded_rather_than_defaulted() {
        // The #726 shape: an un-forwarded method falls through to the trait
        // default, so the wrapper measures itself instead of the provider.
        let audited = provider();
        assert!(audited.has_transcript("NM_TEST.1"));
        assert!(!audited.has_transcript("NM_ABSENT.1"));
        assert!(audited.has_genomic_data() == audited.inner().has_genomic_data());
        assert!(audited.failures().is_empty());
    }
}
