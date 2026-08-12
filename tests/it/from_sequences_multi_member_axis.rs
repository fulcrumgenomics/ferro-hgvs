//! `from_sequences` against the **real-world** multi-member cis alleles — with
//! no reference, so it runs everywhere.
//!
//! `examples/harvest_multi_member_cis.rs` swept ClinVar, CMRG and Paraphase for
//! inputs with more than one member: **592 out of 9 949 738 rows**. Those 592
//! are the entire real-world evidence base for the code path this surface
//! enters, and multi-member is precisely the shape the downstream consumer
//! submits — so this is the axis that most needs to be green on every PR.
//!
//! # Why this used to skip, and why it no longer does
//!
//! The rows name real accessions (`LRG_`, `NC_`, `NM_`), so judging them needed
//! a prepared reference and the axis was manifest-gated. CI has no manifest, so
//! **the only real-world coverage this surface had never ran there** — it
//! reported a skip and exited 0, which is the failure mode the sibling corpus
//! module exists to remove.
//!
//! It did not have to be gated. `from_sequences` reads **no reference** — that
//! is its defining property — so the only thing a provider was supplying was
//! the *window*, and a window is a value. `examples/capture_multi_member_pairs.rs`
//! captures those values once against the prepared reference, and everything
//! below runs against the committed pairs with no provider at all.
//!
//! That is stronger than the `reference-windows.json` pattern used for the
//! projection corpora, not merely smaller: the fixture stores exactly the bytes
//! the function under test consumes, so the capture cannot drift from what the
//! test does with it.
//!
//! The fixture is inventoried in `tests/fixtures/CORPUS_LAYOUT.md`, which records
//! why it is committed at all: `tests/**/*.rs` policy is to generate test data
//! rather than commit it, and this is a named exception with its reason stated
//! at the file rather than only in a commit message.
//!
//! # What the census says, and what it must not be read as
//!
//! Most rows **refuse**, and that is the design envelope rather than a defect
//! rate. Quote the derived count as "how many real-world multi-member alleles
//! this genomic surface can describe", never as a success rate.
//!
//! **The reason is the accession, not the grid — measured, and the opposite of
//! what this paragraph used to say.** It read: "Real submitted multi-member
//! alleles are overwhelmingly structural … and the alignment grid over such a
//! window is orders of magnitude past the budget, so the surface declines and
//! names the knob." Splitting the refusals by message gives **384 accession
//! refusals and 0 grid refusals**, every one of them an `NM_` transcript.
//!
//! The structural rows are real, but they never reach this bucket: the capture
//! drops a row whose span exceeds `MAX_SPAN` (4 096) before the provider is
//! asked, and counts it as `too_wide`. `LRG_542:g.[101177_102434delins36;…]`
//! spans 26 kb and is one of those. So the grid budget explains rows that are
//! **absent from the fixture**, not rows inside it that refuse — which is why
//! the counter is now split ([`Census::refused_accession`] /
//! [`Census::refused_other`]) rather than left as one number with a narrative
//! attached.
//!
//! The counts are not pinned as constants. They move with the harvest and with
//! the capture's own span cap, and a census whose numbers are asserted becomes
//! a thing to re-bless rather than a thing to read. What **is** asserted is that
//! the buckets partition the corpus, that the derived count is non-zero, and
//! that everything which derives is valid and stable — so the absolute claims
//! cannot pass vacuously.

use std::fs;
use std::sync::{Arc, OnceLock};

use serde::Deserialize;

use ferro_hgvs::reference::Transcript;
use ferro_hgvs::{
    from_sequences_detailed, FerroError, FromSequencesOptions, ReferenceProvider, SequencePair,
};

/// Captured by `examples/capture_multi_member_pairs.rs`. Committed, so this
/// module needs no manifest and no provider.
const FIXTURE_PATH: &str = "tests/fixtures/cis/multi_member_sequence_pairs.json.gz";

#[derive(Deserialize)]
struct Pair {
    input: String,
    accession: String,
    position: u64,
    reference: String,
    alternate: String,
    /// The row's SPDI key, computed by the generator against the prepared
    /// reference from the **input** description — the independent oracle.
    canonical_spdi: Option<String>,
}

#[derive(Deserialize)]
struct Captured {
    rows_scanned: usize,
    too_wide: usize,
    unwindowed: usize,
    pairs: Vec<Pair>,
}

fn captured() -> &'static Captured {
    static FIXTURE: OnceLock<Captured> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        // Gzipped to clear the repo's 500 KB pre-commit limit; ~100 KB packed,
        // so it is a plain committed file rather than an out-of-tree release
        // asset — which matters, because a fixture that is absent makes its
        // suite skip green, the exact failure this module was written to
        // remove.
        let file = fs::File::open(FIXTURE_PATH)
            .unwrap_or_else(|e| panic!("failed to open {FIXTURE_PATH}: {e}"));
        let mut text = String::new();
        std::io::Read::read_to_string(&mut flate2::read::GzDecoder::new(file), &mut text)
            .unwrap_or_else(|e| panic!("failed to decompress {FIXTURE_PATH}: {e}"));
        serde_json::from_str(&text)
            .unwrap_or_else(|e| panic!("failed to parse {FIXTURE_PATH}: {e}"))
    })
}

/// A [`ReferenceProvider`] serving exactly one captured window and nothing else.
///
/// This is what makes the denotation oracle reachable without a manifest. The
/// captured `canonical_spdi` was computed at capture time from the **input**
/// description against the prepared reference; to compare a *derived*
/// description against it, the derived side has to be keyed too — and keying
/// goes through `hgvs_to_spdi` and an SPDI splice, which needs a provider.
///
/// The window is the only reference that provider has to serve, because every
/// member of a derived description is inside it by construction. Requests
/// outside it are refused rather than padded with filler: filler would be
/// invented bases, and a key computed against invented bases is worse than no
/// key at all.
///
/// `get_sequence_length` reports the window's own 3' end, which is what lets
/// `apply_to_reference_padded` set `window_is_final` and stop `canonical_spdi`'s
/// widening loop. The consequence is stated where it matters: a change whose 3'
/// roll reaches the window edge is keyed against a shorter runway than the
/// capture had, so [`Census::denotation_edge_bounded`] withholds those rather
/// than reporting a disagreement the fixture caused.
struct WindowProvider {
    accession: String,
    /// 0-based offset of `bases[0]` on `accession`.
    start: u64,
    bases: String,
}

impl WindowProvider {
    fn for_pair(pair: &Pair) -> Self {
        Self {
            accession: pair.accession.clone(),
            // The fixture's `position` is 1-based; providers are 0-based.
            start: pair.position - 1,
            bases: pair.reference.clone(),
        }
    }

    fn end(&self) -> u64 {
        self.start + self.bases.len() as u64
    }

    fn slice(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        if id != self.accession || start < self.start || end > self.end() || end < start {
            return Err(FerroError::ReferenceNotFound { id: id.to_string() });
        }
        let (lo, hi) = ((start - self.start) as usize, (end - self.start) as usize);
        Ok(self.bases[lo..hi].to_string())
    }
}

impl ReferenceProvider for WindowProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        self.slice(id, start, end)
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        self.slice(contig, start, end)
    }

    fn has_genomic_data(&self) -> bool {
        true
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        if id == self.accession {
            Ok(self.end())
        } else {
            Err(FerroError::ReferenceNotFound { id: id.to_string() })
        }
    }
}

/// Where each captured pair ended up. Every pair lands in exactly one bucket.
#[derive(Debug, Default)]
struct Census {
    pairs: usize,
    /// The surface refused because the accession is not genomic — an `NM_`
    /// transcript, in every measured case.
    ///
    /// Split out from a single `refused` counter because that counter carried a
    /// narrative ("the grid budget, overwhelmingly") that the split refutes:
    /// **384 accession, 0 grid**. See the module docs.
    refused_accession: usize,
    /// The surface refused for any other reason — the grid budget, a 5'-edge
    /// anchor. Currently zero, and worth watching: a rise here is a real change
    /// in what the surface can describe, where a rise in `refused_accession` is
    /// only a change in what the harvest contains.
    refused_other: usize,
    /// The pair itself is not one `SequencePair::new` admits.
    ///
    /// **This is a structural zero on this branch, and it is kept so it can stop
    /// being one.** It was written when `validate` admitted only `ACGTN`, to
    /// count the real submitted rows an IUPAC ambiguity code got refused for
    /// (`NM_000518.4:c.[20A>T;249G>Y]` is a ClinVar row). This same change
    /// widened the alphabet to the IUPAC-IUBMB set `general.md:48` admits — the
    /// `alignment-only-symbol-in-a-description` ruling excludes only `X` and `-`
    /// — so the question the bucket described as open was **closed by the
    /// change that introduced the bucket**, and it now counts nothing.
    ///
    /// The doc used to still describe it as an open question, on both this field
    /// and `SequencePair::new`. Corrected in both places rather than deleted: the
    /// bucket keeps the partition total honest, and it is the counter that would
    /// move if the alphabet were ever narrowed again. Read a non-zero here as a
    /// regression in `validate`, not as data about the corpus.
    ambiguity_codes: usize,
    /// A description came back.
    derived: usize,
    /// … of those, the ones whose denotation was compared against the captured
    /// key.
    ///
    /// Counted rather than assumed: a row the generator could not key gives the
    /// comparison nothing to run against, and a skip that reads as a pass is the
    /// failure mode this module exists to remove.
    denotation_compared: usize,
    /// … and the ones withheld because the window is not a wide enough reference
    /// to key against.
    ///
    /// Two ways in. The derived key's own 3' roll reaches the captured window's
    /// edge — the capture keyed the input against the whole contig, so a change
    /// still rolling at the window's end would be keyed against a shorter
    /// runway. Or `canonical_spdi` declines outright, which over this provider
    /// means it asked for bases outside the window.
    ///
    /// Both are artefacts of the fixture rather than disagreements, so they are
    /// withheld — and counted, so the withholding is visible rather than silent.
    /// **Zero on the current capture**, which is the number to read: every
    /// derived row is judged. A rise here is the fixture drifting away from the
    /// windows the oracle needs, not the derivation changing.
    denotation_edge_bounded: usize,
}

fn run() -> (Census, Vec<String>, Vec<String>, Vec<String>, Vec<String>) {
    let options = FromSequencesOptions::default();
    let mut census = Census::default();
    let (mut invalid, mut prohibited, mut unstable, mut redenoted) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());

    for pair in &captured().pairs {
        census.pairs += 1;
        let Ok(sequences) = SequencePair::new(
            pair.accession.clone(),
            pair.position,
            pair.reference.clone(),
            pair.alternate.clone(),
        ) else {
            census.ambiguity_codes += 1;
            continue;
        };

        let derived = match sequences.derive(&options) {
            Ok(derived) => derived,
            Err(refusal) => {
                // Classified on the refusal's own text rather than on the error
                // variant: `InvalidCoordinates` carries the accession refusal
                // *and* the 5'-anchor one, so the variant cannot tell them apart
                // — the same reason `from_sequences_proptest` matches on the
                // anchor clause.
                if refusal
                    .to_string()
                    .contains("names a transcript or protein")
                {
                    census.refused_accession += 1;
                } else {
                    census.refused_other += 1;
                }
                continue;
            }
        };
        census.derived += 1;
        let rendered = derived.variant.to_string();

        if ferro_hgvs::parse_hgvs(&rendered).is_err() && invalid.len() < 5 {
            invalid.push(format!("{} -> {rendered}", pair.input));
        }
        if let Some(clause) = violated_prohibition(&rendered) {
            if prohibited.len() < 5 {
                prohibited.push(format!("{rendered}: {clause}"));
            }
        }
        // Denotation, against the key the generator captured from the *input*
        // against the prepared reference.
        //
        // The derived side is keyed through `spdi::canonical_spdi`, which
        // reaches the bases via `hgvs_to_spdi` and an SPDI splice — deliberately
        // NOT through `apply_edits_to_window`, which is the applier the
        // derivation itself and its internal `verify_round_trip` both use. That
        // separation is the whole point: a fault in the shared applier cancels
        // on both sides of a self-check and cannot cancel here.
        //
        // What this replaces was vacuous. It read
        // `sequences.alternate != pair.alternate`, and `sequences` was built
        // from `pair.alternate.clone()` four statements earlier while
        // `SequencePair::new` stores it verbatim — so the condition was
        // constantly false, `redenoted` could never be pushed, and the
        // `assert!(redenoted.is_empty())` below could never fire. The captured
        // `canonical_spdi` — the entire reason the fixture carries that field —
        // was only ever checked for `is_empty()`. `denotation_compared` still
        // incremented, so the `> 0` guard passed and the dead check read as live.
        if let Some(expected) = pair.canonical_spdi.as_deref().filter(|s| !s.is_empty()) {
            let provider = WindowProvider::for_pair(pair);
            match ferro_hgvs::spdi::canonical_spdi(&derived.variant, &provider) {
                Ok(actual) => {
                    // Withheld, not compared, when the derived key still runs to
                    // the window's 3' edge: the capture keyed against the whole
                    // contig and this provider holds only the window, so the two
                    // rolls had different runways. See
                    // `Census::denotation_edge_bounded`.
                    let key_end = actual.position + actual.deletion.len() as u64;
                    if key_end >= provider.end() {
                        census.denotation_edge_bounded += 1;
                    } else {
                        census.denotation_compared += 1;
                        if actual.to_string() != expected && redenoted.len() < 5 {
                            redenoted.push(format!(
                                "{}: input keyed {expected}, derived {} keys {actual}",
                                pair.input, derived.variant
                            ));
                        }
                    }
                }
                // A derived description the SPDI path cannot key is not evidence
                // about the derivation — the same discipline the sibling corpus
                // applies when either side declines.
                Err(_) => census.denotation_edge_bounded += 1,
            }
        }

        // Idempotence over the same window. The derivation is a pure function of
        // the window, so re-deriving from it must reproduce the answer — and
        // unlike the denotation check below, nothing about this needs a
        // reference.
        if let Ok(again) = from_sequences_detailed(
            &pair.accession,
            pair.position,
            &pair.reference,
            &pair.alternate,
            &options,
        ) {
            if again.variant.to_string() != rendered && unstable.len() < 5 {
                unstable.push(format!("{rendered} -> {}", again.variant));
            }
        }
    }
    (census, invalid, prohibited, unstable, redenoted)
}

/// Prohibitions that make a description invalid outright, not merely
/// non-preferred. Deliberately small; see the sibling corpus module.
fn violated_prohibition(output: &str) -> Option<&'static str> {
    if output.contains(' ') {
        return Some("general.md:96 — spaces are not permitted in any HGVS description");
    }
    if output.contains('X') || output.contains('x') {
        return Some("standards.md:39 — `X` is an alignment-only symbol, not a nucleotide");
    }
    if let Some((_, body)) = output.split_once(":g.") {
        if body.contains('+') || body.contains('-') || body.contains('*') {
            return Some("checklist.md:16 — a `g.` description admits no `+`/`-`/`*` offset");
        }
    }
    None
}

/// **Every derived real-world allele parses, violates no absolute prohibition,
/// denotes what its input denoted, and is a fixed point of its own derivation.**
///
/// The denotation half is checked against an **independent** comparand: the
/// generator computes each row's `canonical_spdi` from the *input* description,
/// against the prepared reference, and commits it. The derived side is keyed
/// through `spdi::canonical_spdi` over a [`WindowProvider`], which reaches the
/// bases via `hgvs_to_spdi` and an SPDI splice rather than through
/// `apply_edits_to_window` — so neither side can agree with the derivation
/// merely because the derivation produced it.
///
/// An earlier version of this module skipped the check, arguing that
/// `from_sequences` verifies its own round trip so the property was already
/// held. That is a self-consistency check — the verifier and the verified are
/// one piece of code, so a bug in the applier passes it — and it is the one
/// construction the repository's test policy names outright.
///
/// **The version after that one was worse, because it looked live.** It gated on
/// `sequences.alternate != pair.alternate`, comparing a field against the value
/// it had just been constructed from, so the condition was constantly false and
/// the `redenoted` assertion below could not fire for any input whatsoever —
/// while `denotation_compared` kept incrementing and the `> 0` guard kept
/// passing. The captured key was read only for `is_empty()`. Three places
/// described that as an independent oracle: this doc, `CORPUS_LAYOUT.md`, and
/// `examples/capture_multi_member_pairs.rs`. All three are corrected.
#[test]
fn every_derived_real_world_allele_is_valid_and_stable() {
    let (census, invalid, prohibited, unstable, redenoted) = run();
    eprintln!("from_sequences multi-member (hermetic): {census:?}");

    assert!(
        census.derived > 0,
        "no real-world row derived, so the assertions below checked nothing: {census:?}"
    );
    assert!(invalid.is_empty(), "unparseable output:\n{invalid:#?}");
    assert!(prohibited.is_empty(), "prohibited output:\n{prohibited:#?}");
    assert!(unstable.is_empty(), "not idempotent:\n{unstable:#?}");
    assert!(redenoted.is_empty(), "denotation changed:\n{redenoted:#?}");
    assert!(
        census.denotation_compared > 0,
        "no row was compared against its independent SPDI key, so the denotation check ran \
         against nothing: {census:?}"
    );
    // A floor on the *comparison*, not merely on its non-emptiness. The
    // `> 0` guard above was already satisfied by the vacuous version this
    // replaced, so it is not on its own evidence that the oracle is doing work.
    assert!(
        census.denotation_compared * 2 > census.derived,
        "fewer than half the derived rows reached the key comparison ({} of {}), so the \
         oracle is mostly withholding rather than judging: {census:?}",
        census.denotation_compared,
        census.derived
    );
}

/// **Every captured pair is accounted for, and the capture accounts for every
/// harvested row.**
///
/// Two levels, because a row can go missing at either. Without the second, the
/// derived count could quietly be computed against a capture that had silently
/// shrunk — the failure mode `CaptureLedger` exists for, applied by hand
/// because this is a test rather than a generator writing an artifact.
#[test]
fn the_census_accounts_for_every_harvested_row() {
    let (census, ..) = run();
    let fixture = captured();

    assert_eq!(
        census.refused_accession + census.refused_other + census.derived + census.ambiguity_codes,
        census.pairs,
        "the buckets do not sum to the captured pairs: {census:?}"
    );
    assert_eq!(
        census.pairs,
        fixture.pairs.len(),
        "the capture did not load in full"
    );
    assert_eq!(
        fixture.pairs.len() + fixture.too_wide + fixture.unwindowed,
        fixture.rows_scanned,
        "the capture dropped rows it did not account for: {} pairs + {} too wide + {} unwindowed \
         != {} scanned",
        fixture.pairs.len(),
        fixture.too_wide,
        fixture.unwindowed,
        fixture.rows_scanned
    );
}

/// **This axis needs no reference**, which is what lets it run in CI at all.
///
/// Pinned rather than left implicit: the module's whole reason for existing is
/// that the manifest gate came off, and a future edit that reaches for a
/// provider would silently re-gate it — reporting a skip and exiting 0, which
/// is exactly the state this replaced.
#[test]
fn the_axis_is_hermetic() {
    let pair = &captured().pairs[0];
    let sequences = SequencePair::new(
        pair.accession.clone(),
        pair.position,
        pair.reference.clone(),
        pair.alternate.clone(),
    )
    .expect("a captured pair is well-formed");
    // No `Normalizer`, no provider, no manifest — the derivation takes only the
    // four values the fixture carries, and at least one of them derives.
    //
    // Asserted rather than merely called: `let _ = derive(..)` would pass even
    // if every input refused, so it pinned nothing at all.
    let _ = sequences;
    assert!(
        captured()
            .pairs
            .iter()
            .filter_map(|p| SequencePair::new(
                p.accession.clone(),
                p.position,
                p.reference.clone(),
                p.alternate.clone(),
            )
            .ok())
            .any(|p| p.derive(&FromSequencesOptions::default()).is_ok()),
        "no captured pair derives without a provider, so this axis is not \
         exercising the hermetic path it claims to pin"
    );
}
