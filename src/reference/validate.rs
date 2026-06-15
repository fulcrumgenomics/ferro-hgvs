//! Validation of transcript reference records (issue #520).
//!
//! Two tiers, both reporting [`TranscriptAnomaly`] tagged with an
//! [`AnomalyConfidence`]:
//!
//! 1. **Offline structural** ([`validate_transcript_record`]) — self-contained
//!    checks over a single [`Transcript`]'s own sequence, CDS coordinates, and
//!    exon alignment; no network, no external record.
//! 2. **Authoritative comparison** ([`validate_against_authoritative`]) —
//!    compares the served transcript's length, CDS bounds, and paired protein
//!    against the canonical RefSeq record for the exact accession version (an
//!    [`AuthoritativeRecord`], parsed from GenBank in
//!    [`crate::reference::authoritative`]). This is the precise tier that
//!    catches the mis-pairing / non-canonical-length cases the offline tier
//!    cannot (e.g. `NM_012459.2` paired with `NP_036591.3`; `NM_000193.2`
//!    served at 4650 nt vs the canonical 1576 nt).
//!
//! ## Offline structural tier
//!
//! Anomalies are tagged:
//!
//! - **`Corruption`** — high-confidence structural problems: the served
//!   sequence is *shorter* than the exon alignment needs (missing bases), the
//!   CDS coordinates are out of range, the CDS length is not a triplet, or the
//!   record fails to load at all.
//! - **`Advisory`** — open-reading-frame heuristics (missing start codon,
//!   missing stop codon, premature internal stop) that have **known biological
//!   false positives**: selenoproteins recode an in-frame `TGA` as
//!   selenocysteine, some transcripts use alternative (`CTG`/`GTG`/…) or absent
//!   start codons, stop-codon readthrough exists, and `partial` RefSeq records
//!   legitimately lack a start or stop. Treat these as "review", not proof.
//!
//! Deliberately out of scope for the **offline** tier (handled by the
//! authoritative tier instead):
//!
//! - A served sequence *longer* than the exon alignment is **not** flagged
//!   offline — a poly-A tail or unaligned 3' bases routinely make the
//!   transcript FASTA longer than the genome-derived alignment. The precise
//!   over-length check (e.g. `NM_000193.2`) lives in
//!   [`validate_against_authoritative`], which compares against the canonical
//!   transcript length.
//! - A sequence internally consistent but mis-paired to the wrong protein
//!   version (e.g. `NM_012459.2`) is invisible offline; the authoritative tier
//!   catches it via the CDS / protein_id comparison.

use crate::backtranslate::codon::{Codon, CodonTable};
use crate::error::FerroError;
use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::{Exon, Transcript};

/// The kind of structural anomaly found in a transcript record.
///
/// `#[non_exhaustive]`: later phases add authoritative-comparison variants, so
/// downstream `match`es must include a wildcard arm.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AnomalyKind {
    /// The served sequence is shorter than the transcript length implied by the
    /// exon alignment (max exon end): the alignment references bases the served
    /// sequence does not contain.
    LengthMismatch,
    /// The CDS coordinates extend beyond the available sequence or past the
    /// exon-alignment extent (or are otherwise degenerate).
    CdsOutOfRange,
    /// The CDS length is not a multiple of three.
    CdsNotTriplet,
    /// The CDS does not begin with a start codon (`ATG`). Advisory — see the
    /// module docs on alternative/absent start codons.
    MissingStartCodon,
    /// The CDS does not end with a stop codon (`TAA`/`TAG`/`TGA`). Advisory.
    MissingStopCodon,
    /// The CDS contains a premature (internal) stop codon. Advisory — an
    /// in-frame `TGA` may be recoded selenocysteine, not a true stop.
    InternalStopCodon,
    /// The record failed to load from the provider for a reason other than
    /// "not found" (e.g. degenerate coordinates rejected during construction).
    LoadError,
    /// The served transcript length disagrees with the authoritative RefSeq
    /// record for the exact accession version.
    AuthoritativeLengthMismatch,
    /// The served CDS coordinates disagree with the authoritative record.
    AuthoritativeCdsMismatch,
    /// The served paired protein accession disagrees with the authoritative
    /// record (e.g. `NM_012459.2` served with `NP_036591.3` instead of the
    /// canonical `NP_036591.2`).
    AuthoritativeProteinIdMismatch,
    /// Translating the served CDS does not reproduce the canonical protein
    /// sequence. Catches base-level non-canonical sequences that the
    /// coordinate/length/protein_id comparisons miss (e.g. a same-length CDS
    /// with point edits). Advisory — the standard codon table cannot model
    /// selenocysteine (`TGA`) recoding or alternative start codons, so it has
    /// the same biological false-positive classes as the offline codon checks
    /// (the obvious selenoprotein readthrough case is suppressed; see
    /// [`validate_translation_against_protein`]).
    TranslationMismatch,
}

/// How much confidence an anomaly carries that the record is genuinely corrupt.
///
/// `#[non_exhaustive]` for the same forward-compatibility reason as
/// [`AnomalyKind`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AnomalyConfidence {
    /// High-confidence structural corruption (length/frame/coords/load).
    Corruption,
    /// A heuristic open-reading-frame check with known biological false
    /// positives (selenocysteine `TGA`, alternative/absent start codons,
    /// readthrough, partial CDS). Advisory only.
    Advisory,
}

/// A single anomaly found while validating a transcript record.
///
/// `#[non_exhaustive]`: later phases may add fields (e.g. a structured
/// position); construct it only within this crate.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct TranscriptAnomaly {
    /// The transcript accession the anomaly was found on.
    pub transcript_id: String,
    /// The kind of anomaly.
    pub kind: AnomalyKind,
    /// How much confidence the anomaly carries (see [`AnomalyConfidence`]).
    pub confidence: AnomalyConfidence,
    /// A human-readable description with the concrete values involved.
    pub detail: String,
}

impl TranscriptAnomaly {
    fn new(id: &str, kind: AnomalyKind, confidence: AnomalyConfidence, detail: String) -> Self {
        Self {
            transcript_id: id.to_string(),
            kind,
            confidence,
            detail,
        }
    }
}

/// Validate a single transcript record's internal structural consistency.
///
/// Returns one [`TranscriptAnomaly`] per problem found (possibly several).
/// A record with no sequence, or no CDS coordinates, simply skips the checks
/// that need that data — coordinate-only and non-coding records are not
/// anomalies in themselves.
pub fn validate_transcript_record(tx: &Transcript) -> Vec<TranscriptAnomaly> {
    let mut anomalies = Vec::new();

    let Some(sequence) = tx.sequence.as_deref() else {
        // Coordinate-only record: nothing base-level to validate.
        return anomalies;
    };
    // Work on bytes throughout: nucleotide data is ASCII, but a stray
    // multi-byte char in a malformed record must not panic an `&str` slice.
    let bytes = sequence.as_bytes();
    let seq_len = bytes.len() as u64;

    // --- Length consistency: served bases vs exon-alignment extent ----------
    // The exon alignment tiles transcript coordinates [1, max_end]
    // contiguously (gaps are rejected upstream during ingestion). A served
    // sequence *shorter* than that extent means the alignment references bases
    // the sequence lacks — unambiguous corruption. A *longer* sequence is the
    // normal poly-A / unaligned-3' case and is intentionally NOT flagged here
    // (see the module docs).
    let exon_extent = tx.exons.iter().map(|e| e.end).max();
    if let Some(extent) = exon_extent {
        if seq_len < extent {
            anomalies.push(TranscriptAnomaly::new(
                &tx.id,
                AnomalyKind::LengthMismatch,
                AnomalyConfidence::Corruption,
                format!(
                    "served sequence is {seq_len} nt, shorter than the {extent} nt \
                     spanned by the exon alignment (bases are missing)"
                ),
            ));
        }
    }

    // --- CDS open-reading-frame sanity --------------------------------------
    // Only a fully-absent CDS `(None, None)` is the legitimate non-coding case.
    // A half-populated pair is structurally incomplete and must be flagged as
    // corruption rather than silently passing as non-coding.
    let (cds_start, cds_end) = match (tx.cds_start, tx.cds_end) {
        (Some(start), Some(end)) => (start, end),
        (None, None) => return anomalies,
        (start, end) => {
            anomalies.push(TranscriptAnomaly::new(
                &tx.id,
                AnomalyKind::CdsOutOfRange,
                AnomalyConfidence::Corruption,
                format!("CDS bounds are incomplete: start={start:?}, end={end:?}"),
            ));
            return anomalies;
        }
    };

    // `cds_start`/`cds_end` are 1-based inclusive transcript coordinates. The
    // CDS must lie within both the served sequence *and* the exon alignment: a
    // CDS end past the exon extent points into the tolerated poly-A/unaligned
    // tail (`exon_extent < cds_end <= seq_len`), i.e. coding coordinates outside
    // the alignment — corruption the `cds_end > seq_len` bound alone misses.
    let cds_past_exon_extent = exon_extent.is_some_and(|extent| cds_end > extent);
    if cds_start == 0 || cds_end < cds_start || cds_end > seq_len || cds_past_exon_extent {
        anomalies.push(TranscriptAnomaly::new(
            &tx.id,
            AnomalyKind::CdsOutOfRange,
            AnomalyConfidence::Corruption,
            format!(
                "CDS {cds_start}..{cds_end} (1-based) is invalid for a {seq_len} nt sequence \
                 or exceeds the exon-alignment extent"
            ),
        ));
        // Without an in-range CDS the codon-level checks cannot run.
        return anomalies;
    }

    // Byte slice is safe: bounds were just checked against the byte length.
    let cds: Vec<u8> = bytes[(cds_start as usize - 1)..(cds_end as usize)].to_ascii_uppercase();

    if !cds.len().is_multiple_of(3) {
        anomalies.push(TranscriptAnomaly::new(
            &tx.id,
            AnomalyKind::CdsNotTriplet,
            AnomalyConfidence::Corruption,
            format!("CDS length {} is not a multiple of 3", cds.len()),
        ));
        // A non-triplet CDS makes codon framing ambiguous; stop here.
        return anomalies;
    }

    let codons: Vec<&[u8]> = cds.chunks(3).collect();

    // Start codon (advisory). Skip ambiguous codons (containing N/IUPAC) — we
    // cannot tell whether they are a real start.
    if let Some(first) = codons.first() {
        if is_acgt_codon(first) && !is_start_codon(first) {
            anomalies.push(TranscriptAnomaly::new(
                &tx.id,
                AnomalyKind::MissingStartCodon,
                AnomalyConfidence::Advisory,
                format!(
                    "CDS does not begin with ATG (starts with {})",
                    String::from_utf8_lossy(first)
                ),
            ));
        }
    }

    // Stop codon (advisory).
    if let Some(last) = codons.last() {
        if is_acgt_codon(last) && !is_stop_codon(last) {
            anomalies.push(TranscriptAnomaly::new(
                &tx.id,
                AnomalyKind::MissingStopCodon,
                AnomalyConfidence::Advisory,
                format!(
                    "CDS does not end with a stop codon (ends with {})",
                    String::from_utf8_lossy(last)
                ),
            ));
        }
    }

    // Internal stop (advisory): an unambiguous stop codon anywhere before the
    // final codon. Ambiguous codons are skipped; a recoded `TGA`
    // (selenocysteine) is still reported here but as Advisory, by design.
    let internal = &codons[..codons.len().saturating_sub(1)];
    if let Some(pos) = internal
        .iter()
        .position(|c| is_acgt_codon(c) && is_stop_codon(c))
    {
        anomalies.push(TranscriptAnomaly::new(
            &tx.id,
            AnomalyKind::InternalStopCodon,
            AnomalyConfidence::Advisory,
            format!(
                "in-frame stop codon at CDS codon {} ({}); may be recoded selenocysteine",
                pos + 1,
                String::from_utf8_lossy(internal[pos])
            ),
        ));
    }

    anomalies
}

/// Validate a set of transcripts pulled from `provider` by accession, returning
/// every anomaly found across them.
///
/// An accession the provider does not have (`ReferenceNotFound`) is skipped —
/// a missing transcript is a different concern (handled by `ferro check`'s
/// file/manifest validation). Any *other* load error is surfaced as a
/// [`AnomalyKind::LoadError`] anomaly rather than silently dropped: a record
/// that fails to even construct (e.g. degenerate cdot coordinates) is exactly
/// the corruption this tool exists to find.
///
/// Callers pass the accession set they care about (e.g. a corpus/used set)
/// rather than scanning the full ~100k-record reference, which is left to a
/// dedicated full-scan path.
pub fn validate_transcripts<P: ReferenceProvider>(
    provider: &P,
    ids: &[String],
) -> Vec<TranscriptAnomaly> {
    let mut anomalies = Vec::new();
    for id in ids {
        match provider.get_transcript(id) {
            Ok(tx) => anomalies.extend(validate_transcript_record(&tx)),
            Err(FerroError::ReferenceNotFound { .. }) => {} // unknown accession: not our concern
            Err(e) => anomalies.push(TranscriptAnomaly::new(
                id,
                AnomalyKind::LoadError,
                AnomalyConfidence::Corruption,
                format!("transcript failed to load: {e}"),
            )),
        }
    }
    anomalies
}

/// Compare a served transcript against the authoritative RefSeq record for its
/// exact accession version, returning a [`TranscriptAnomaly`] for each
/// disagreement.
///
/// Unlike the offline structural checks ([`validate_transcript_record`]), these
/// are precise: the authoritative record is the canonical truth for the exact
/// version, so a disagreement is high-confidence corruption (`Corruption`). All
/// comparisons are in **1-based** coordinates (the served `Transcript`'s CDS is
/// already 1-based; the [`AuthoritativeRecord`] is 1-based by construction).
///
/// Checks:
/// - **length**: served sequence length vs `tx_length`;
/// - **CDS**: served `(cds_start, cds_end)` vs the authoritative bounds;
/// - **protein_id**: served paired protein vs the authoritative `/protein_id`.
///
/// CDS and protein_id are compared only when **both** sides carry the value. A
/// served `None` is deliberately NOT treated as a "non-coding / no-protein vs
/// coding" mismatch: the provider commonly leaves these fields unpopulated for
/// reasons unrelated to corruption (e.g. the cdot-synthesis path sets
/// `protein_id: None`), so flagging a presence difference would be
/// false-positive-prone. A presence-aware check belongs to a richer
/// authoritative ingestion, not this comparison.
pub fn validate_against_authoritative(
    tx: &Transcript,
    auth: &AuthoritativeRecord,
) -> Vec<TranscriptAnomaly> {
    let mut anomalies = Vec::new();
    let id = &tx.id;

    if let Some(seq) = tx.sequence.as_deref() {
        let served = seq.len() as u64;
        if served != auth.tx_length {
            anomalies.push(TranscriptAnomaly::new(
                id,
                AnomalyKind::AuthoritativeLengthMismatch,
                AnomalyConfidence::Corruption,
                format!(
                    "served sequence is {served} nt but the authoritative record \
                     {} is {} nt",
                    auth.accession, auth.tx_length
                ),
            ));
        }
    }

    // Both sides are 1-based inclusive. `Transcript.cds_end` needs no
    // conversion from cdot because cdot's 0-based *exclusive* end is
    // numerically identical to a 1-based inclusive end (see cdot.rs and
    // multi_fasta.rs CDS construction); the authoritative `cds_end` is the
    // GenBank `a..b` value, also 1-based inclusive.
    match (tx.cds_start, tx.cds_end, auth.cds_start, auth.cds_end) {
        (Some(ts), Some(te), Some(as_), Some(ae)) if (ts, te) != (as_, ae) => {
            anomalies.push(TranscriptAnomaly::new(
                id,
                AnomalyKind::AuthoritativeCdsMismatch,
                AnomalyConfidence::Corruption,
                format!(
                    "served CDS {ts}..{te} (1-based) disagrees with the authoritative \
                     {as_}..{ae}"
                ),
            ));
        }
        // A served record with exactly one CDS bound set is structurally
        // incomplete: surface corruption instead of passing clean because the
        // four-`Some` arm failed to match.
        (Some(_), None, _, _) | (None, Some(_), _, _) => {
            anomalies.push(TranscriptAnomaly::new(
                id,
                AnomalyKind::CdsOutOfRange,
                AnomalyConfidence::Corruption,
                "served CDS bounds are incomplete".to_string(),
            ));
        }
        _ => {}
    }

    if let (Some(tp), Some(ap)) = (tx.protein_id.as_deref(), auth.protein_id.as_deref()) {
        if tp != ap {
            anomalies.push(TranscriptAnomaly::new(
                id,
                AnomalyKind::AuthoritativeProteinIdMismatch,
                AnomalyConfidence::Corruption,
                format!("served protein {tp} disagrees with the authoritative {ap}"),
            ));
        }
    }

    anomalies
}

/// Validate every transcript named in `overrides` against its authoritative
/// record, pulling the served transcript from `provider`.
///
/// As in [`validate_transcripts`], an accession the provider does not have is
/// skipped (`ReferenceNotFound`); any other load error surfaces as a
/// [`AnomalyKind::LoadError`] anomaly. The skip is deliberate even though we
/// hold authoritative data for the accession: "we have the canonical record
/// but the reference does not ship this transcript" is an inventory/coverage
/// concern for the prepare/check layer, not a corruption of a served record.
pub fn validate_against_overrides<P: ReferenceProvider>(
    provider: &P,
    overrides: &CanonicalOverrides,
) -> Vec<TranscriptAnomaly> {
    let mut anomalies = Vec::new();
    for (accession, auth) in &overrides.records {
        match provider.get_transcript(accession) {
            Ok(tx) => anomalies.extend(validate_against_authoritative(&tx, auth)),
            Err(FerroError::ReferenceNotFound { .. }) => {}
            Err(e) => anomalies.push(TranscriptAnomaly::new(
                accession,
                AnomalyKind::LoadError,
                AnomalyConfidence::Corruption,
                format!("transcript failed to load: {e}"),
            )),
        }
    }
    anomalies
}

/// A single correction applied to (or required for) a served transcript when
/// reconciling it with the authoritative record.
///
/// `#[non_exhaustive]`: later phases may add correction kinds.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Correction {
    /// The CDS bounds were overridden from the authoritative record.
    CdsCorrected {
        /// The served `(cds_start, cds_end)` before correction.
        from: (Option<u64>, Option<u64>),
        /// The authoritative 1-based inclusive `(cds_start, cds_end)`.
        to: (u64, u64),
    },
    /// The paired protein accession was overridden from the authoritative record.
    ProteinIdCorrected {
        /// The served `protein_id` before correction.
        from: Option<String>,
        /// The authoritative `protein_id`.
        to: String,
    },
    /// The served bases were **replaced** with the authoritative sequence — the
    /// served sequence was the wrong length, or a same-length edit. Only
    /// possible when the override carries the canonical `sequence`.
    SequenceCorrected {
        /// Served sequence length before correction (0 if coordinate-only).
        from: u64,
        /// Authoritative sequence length.
        to: u64,
    },
    /// The served sequence length disagrees with the authoritative record AND
    /// the override carries no canonical `sequence` to replace it with, so
    /// CDS/protein were **not** overridden (the coordinates index a
    /// different-length sequence). Fetch the canonical bases (carry `sequence`
    /// in the overrides) to make this correctable.
    SequenceReingestionRequired {
        /// Served sequence length.
        served: u64,
        /// Authoritative transcript length.
        authoritative: u64,
    },
}

/// Reconcile a served transcript's CDS/protein metadata with the authoritative
/// record for its exact accession version, mutating `tx` in place and returning
/// the corrections applied.
///
/// This fixes the **metadata-mismatch** class (e.g. `NM_012459.2`, whose served
/// `.2` bases are canonical but whose cdot CDS/`protein_id` point at the `.3`
/// short isoform): the CDS bounds and `protein_id` are overridden from the
/// authoritative record so downstream protein prediction translates the
/// canonical reading frame.
///
/// It also fixes the **wrong-sequence** class (e.g. `NM_000193.2` served at
/// 4650 nt vs the canonical 1576) when the override carries the canonical
/// `sequence`: the served bases are replaced ([`Correction::SequenceCorrected`])
/// — this also catches a same-length-but-edited sequence (byte comparison, not
/// just length). When the override carries **no** `sequence` and the served
/// length disagrees with the authoritative length, the coordinates can't be
/// applied safely, so nothing is overridden and a
/// [`Correction::SequenceReingestionRequired`] is reported instead (carry the
/// canonical bases to make it correctable). The served length is the sequence
/// byte length, or, for a coordinate-only record, the exon-sum length.
///
/// CDS and `protein_id` are corrected **together or not at all**: a corrected
/// reading frame paired with a stale protein label (or vice versa) is more
/// misleading than no correction, so the override is applied only when the
/// authoritative record carries the full coding triple (both CDS bounds *and* a
/// `protein_id`).
///
/// # Limitations
///
/// - **Same-length-but-edited** sequences are only caught when the override
///   carries the canonical `sequence` (byte comparison); without it, length
///   equality is the only proxy and a same-length edit passes as canonical.
/// - Replacing the bases does **not** rewrite the exon structure (the
///   authoritative exon layout isn't carried). Protein prediction reads the CDS
///   bases via `cds_start`/`cds_end`, so it is correct; genomic projection on a
///   sequence-corrected transcript may be inconsistent with the stale exons.
/// - This corrects only the served [`Transcript`]. The projection path also
///   reads the **cdot mapper's** parallel `protein`/CDS fields and prefers them,
///   so for full effect the wiring must reconcile the cdot record too (or route
///   projection through the corrected transcript).
/// - **Ordering contract:** call this on a freshly loaded transcript *before*
///   it is cached or translated. Mutating an already-translated transcript
///   would leave a cached reference translation (the projector's
///   `ref_protein_cache`) stale.
///
/// Returns an empty vec when there is no override for `tx.id` or nothing needs
/// changing.
pub fn apply_canonical_overrides(
    tx: &mut Transcript,
    overrides: &CanonicalOverrides,
) -> Vec<Correction> {
    let Some(auth) = overrides.get(&tx.id) else {
        return Vec::new();
    };
    // `overrides` is file-backed: a hand-edited / corrupt entry whose embedded
    // accession disagrees with its map key would rewrite this record with
    // another accession's metadata. Trust the override only when it self-reports
    // the accession we looked it up by.
    if auth.accession != tx.id {
        return Vec::new();
    }
    let mut corrections = Vec::new();

    // Sequence reconciliation:
    // - If the override carries canonical bases, replace a non-canonical served
    //   sequence (wrong-length OR same-length-but-edited) so the authoritative
    //   coordinates apply to the right bases.
    // - Otherwise fall back to the length gate: only proceed when the served
    //   length matches the authoritative length (the coordinates would
    //   otherwise index a different-length sequence); if not, report that the
    //   canonical bases need to be carried/re-ingested and stop.
    match auth.sequence.as_deref() {
        Some(canonical) => {
            // Compare case-insensitively: a soft-masked (lowercase) served
            // sequence that is otherwise canonical is NOT a real edit, so don't
            // report a spurious correction. Only genuine base differences (or a
            // length difference) trigger a replacement.
            let differs = tx
                .sequence
                .as_deref()
                .is_none_or(|s| !s.eq_ignore_ascii_case(canonical));
            if differs {
                let from = tx.sequence.as_deref().map_or(0, |s| s.len() as u64);
                let to = canonical.len() as u64;
                corrections.push(Correction::SequenceCorrected { from, to });
                tx.sequence = Some(canonical.to_string());
                // A length-changing replacement leaves the exon structure (which
                // tiled the old length) stale — e.g. it would trip the offline
                // `LengthMismatch` check (served shorter than the exon extent).
                // The authoritative exon layout isn't carried, and was wrong for
                // a wrong-length record anyway, so collapse to a single exon
                // spanning the corrected transcript (mRNA is contiguous in
                // transcript space). Genomic projection on such a record is
                // already out of scope (see the doc caveat).
                if from != to {
                    tx.exons = vec![Exon::new(1, 1, to)];
                }
            }
        }
        None => {
            // Served length: sequence byte length, else the exon-sum length for
            // a coordinate-only record. The exon-sum branch assumes the 1-based
            // *inclusive* `Transcript` exon convention (`end - start + 1`); it is
            // reached only when `tx.sequence` is `None`, which a sequence-bearing
            // (FASTA-backed or cdot-synthesized) transcript never is — so the
            // synthesize path's 0-based exon coordinates don't flow through here.
            let served_len = tx.sequence.as_deref().map(|s| s.len() as u64).or_else(|| {
                let sum: u64 = tx
                    .exons
                    .iter()
                    .map(|e| e.end.saturating_sub(e.start) + 1)
                    .sum();
                (sum > 0).then_some(sum)
            });
            if let Some(served) = served_len {
                if served != auth.tx_length {
                    corrections.push(Correction::SequenceReingestionRequired {
                        served,
                        authoritative: auth.tx_length,
                    });
                    return corrections;
                }
            }
        }
    }

    // All-or-nothing for the coding triple: only override when the authoritative
    // record is fully coding (CDS bounds + protein_id present).
    let (Some(as_), Some(ae), Some(ap)) =
        (auth.cds_start, auth.cds_end, auth.protein_id.as_deref())
    else {
        return corrections;
    };

    // Reject structurally-impossible CDS bounds before mutating: require
    // `1 <= cds_start <= cds_end <= tx_length`. A malformed override past the
    // sequence length would index out of range and corrupt translation.
    if as_ == 0 || ae < as_ || ae > auth.tx_length {
        return corrections;
    }

    if (tx.cds_start, tx.cds_end) != (Some(as_), Some(ae)) {
        corrections.push(Correction::CdsCorrected {
            from: (tx.cds_start, tx.cds_end),
            to: (as_, ae),
        });
        tx.cds_start = Some(as_);
        tx.cds_end = Some(ae);
    }

    if tx.protein_id.as_deref() != Some(ap) {
        corrections.push(Correction::ProteinIdCorrected {
            from: tx.protein_id.clone(),
            to: ap.to_string(),
        });
        tx.protein_id = Some(ap.to_string());
    }

    corrections
}

/// The protein produced by translating a CDS, plus whether translation stopped
/// on a `TGA` codon specifically. Only `TGA` is recoded as selenocysteine, so
/// the caller needs to distinguish a `TGA` stop (candidate Sec readthrough)
/// from a `TAA`/`TAG` stop (a true nonsense codon) — see
/// [`validate_translation_against_protein`].
struct TranslatedCds {
    /// One-letter amino-acid string, with the terminating stop omitted.
    protein: String,
    /// Whether translation terminated on an in-frame `TGA` (as opposed to
    /// `TAA`/`TAG`, or running off the end of the CDS).
    terminated_by_tga: bool,
}

/// Translate a CDS byte slice to a one-letter amino-acid string, stopping at
/// the first stop codon (the stop is not emitted). Returns `None` if any codon
/// is ambiguous/unparseable (`N`/IUPAC) — one ambiguous codon aborts the whole
/// translation so the caller skips the comparison rather than reporting a
/// spurious mismatch.
///
/// The returned [`TranslatedCds`] also records whether the terminating stop was
/// `TGA` (the only stop recoded as selenocysteine), so the caller can suppress
/// Sec readthrough without also suppressing a `TAA`/`TAG` nonsense edit at the
/// same site.
///
/// Lowercase / soft-masked bases are fine: `Codon::parse_bytes` → `Base::from_char`
/// upper-cases each base, so no pre-upper-casing is needed here (unlike the
/// string-based offline codon checks). A trailing partial (<3 nt) codon is
/// dropped — an in-range CDS is expected to be a triplet (the offline
/// `CdsNotTriplet` check guards that separately).
fn translate_cds_bytes(cds: &[u8]) -> Option<TranslatedCds> {
    let table = CodonTable::standard();
    let mut protein = String::with_capacity(cds.len() / 3);
    for chunk in cds.chunks(3) {
        if chunk.len() < 3 {
            break; // trailing partial codon
        }
        let codon = Codon::parse_bytes(chunk)?;
        if table.is_stop(&codon) {
            return Some(TranslatedCds {
                protein,
                terminated_by_tga: chunk.eq_ignore_ascii_case(b"TGA"),
            });
        }
        protein.push(table.amino_acid_for(&codon)?.to_one_letter());
    }
    Some(TranslatedCds {
        protein,
        terminated_by_tga: false,
    })
}

/// Compare the protein obtained by translating the served CDS against the
/// canonical protein sequence, returning an [`AnomalyKind::TranslationMismatch`]
/// (`Advisory`) anomaly when they differ.
///
/// This is the base-level (strategy-B) check: the coordinate / length /
/// protein_id comparisons in [`validate_against_authoritative`] can all agree
/// while the served *bases* still translate to the wrong protein (a same-length
/// CDS with point edits). Supplying `canonical_protein` requires the
/// authoritative `NP_` sequence — i.e. protein-FASTA ingestion, the data edge
/// not covered here.
///
/// Biological adjustments so canonical RefSeq records don't false-positive:
/// - the first residue is forced to `M` (translation initiation installs Met
///   regardless of the start codon; canonical `NP_` always begins with `M`);
/// - the obvious **selenocysteine readthrough** case is suppressed — the
///   standard table stops at an in-frame `TGA`, so a selenoprotein translates
///   to a prefix ending exactly where the canonical protein has a `U` (Sec);
///   that specific shape returns `None` rather than a mismatch.
///
/// The anomaly is `Advisory`, not `Corruption`: the standard codon table cannot
/// model every recoding/alternative-start case, so a residual mismatch is a
/// review signal, not proof of corruption.
///
/// Returns `None` (no anomaly) when the record is non-coding, has no sequence,
/// has out-of-range CDS bounds, or the CDS contains ambiguous codons that can't
/// be translated.
pub fn validate_translation_against_protein(
    tx: &Transcript,
    canonical_protein: &str,
) -> Option<TranscriptAnomaly> {
    let (cds_start, cds_end) = (tx.cds_start?, tx.cds_end?);
    let seq = tx.sequence.as_deref()?;
    let bytes = seq.as_bytes();
    if cds_start == 0 || cds_end < cds_start || cds_end as usize > bytes.len() {
        return None;
    }
    let TranslatedCds {
        mut protein,
        terminated_by_tga,
    } = translate_cds_bytes(&bytes[(cds_start as usize - 1)..(cds_end as usize)])?;
    // Translation initiation installs Met regardless of the start codon, so the
    // canonical NP_ always begins with M; mirror that to avoid a position-1
    // false positive on alternative (CTG/GTG/…) start codons.
    if !protein.is_empty() {
        protein.replace_range(0..1, "M");
    }
    let translated = protein;

    // Canonical `NP_` sequences carry no trailing stop; tolerate a stray `*`.
    let canonical = canonical_protein.trim_end_matches('*');
    if translated == canonical {
        return None;
    }

    // Selenocysteine readthrough: the standard table stops at the first
    // in-frame TGA, so a selenoprotein translates to a prefix that ends exactly
    // where the canonical protein has a `U` (Sec). Suppress that expected case,
    // but only when the terminating codon was actually `TGA` — a `TAA`/`TAG`
    // stop at the same position is a genuine nonsense edit, not Sec recoding,
    // and must still be flagged.
    if terminated_by_tga
        && canonical.starts_with(&translated)
        && canonical.as_bytes().get(translated.len()) == Some(&b'U')
    {
        return None;
    }

    Some(TranscriptAnomaly::new(
        &tx.id,
        AnomalyKind::TranslationMismatch,
        AnomalyConfidence::Advisory,
        format!(
            "translated served CDS ({} aa) does not match the canonical protein ({} aa)",
            translated.len(),
            canonical.len()
        ),
    ))
}

/// A coding transcript whose CDS starts with **neither `ATG` nor a recognized
/// alternative initiator** (`CTG`/`GTG`/`TTG`) in the served FASTA (issue #629).
///
/// This is the data-side signal that the cdot CDS coordinates and the
/// transcript FASTA disagree (the cdot annotation release and the mRNA FASTA
/// release are different sequence revisions): `cds_start` has drifted off the
/// true start, so the reading frame is wrong and the reference protein
/// mistranslates. #625 makes the *runtime* projection robust to this (it
/// declines protein prediction); this scan surfaces the *underlying data*
/// inconsistency at `ferro check` time so a bad reference set is caught up
/// front rather than silently losing protein-consequence support.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdsStartInconsistency {
    /// The transcript accession (with version), e.g. `NM_000425.3`.
    pub transcript_id: String,
    /// The first CDS codon actually observed in the served FASTA, upper-cased,
    /// e.g. `AGG` — a codon that is not a recognized translation initiator.
    pub observed_start_codon: String,
}

/// Summary of a [`scan_cds_start_codons`] pass over a set of transcripts (#629).
///
/// Only the **start codon** is examined — deliberately *not* internal stops:
/// selenoproteins legitimately recode an in-frame `TGA` as selenocysteine, so a
/// premature-internal-stop scan would false-positive on every selenoprotein.
///
/// A non-`ATG` start is **not** automatically an inconsistency. The near-cognate
/// initiation codons `CTG`/`GTG`/`TTG` are recognized translation start codons
/// (RefSeq annotates them; the reading frame is correct and the initiator Met is
/// still installed), so they are counted as [`alternative_start`](Self::alternative_start)
/// rather than flagged — they do **not** cause the mistranslation #625/#629 are
/// about. Only a first codon that is neither `ATG` nor a recognized alternative
/// initiator (e.g. `AGC`, `GAT`, `CAG`) signals that `cds_start` has drifted off
/// the true start — the data-side root cause behind #625. This mirrors the
/// canonical-start set already used by the GFF loader's `W-LOAD-201` check.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CdsStartCodonReport {
    /// Number of coding transcripts whose start codon could be examined
    /// (CDS present, sequence available, first codon in range and unambiguous).
    pub coding_examined: usize,
    /// Coding transcripts whose first CDS codon is neither `ATG` nor a
    /// recognized alternative initiator, and that are not on the allowlist —
    /// the genuine `cds_start`-off-the-start signal. Sorted by `transcript_id`.
    pub inconsistent: Vec<CdsStartInconsistency>,
    /// Coding transcripts starting with a recognized alternative initiation
    /// codon (`CTG`/`GTG`/`TTG`). Legitimate; counted for transparency, not
    /// flagged.
    pub alternative_start: usize,
    /// Coding transcripts that *would* have been flagged but were suppressed by
    /// the user-supplied allowlist (known-legitimate non-initiator starts).
    pub allowlisted: usize,
    /// Coding transcripts whose first codon contained an ambiguity code
    /// (`N`/IUPAC) and so could not be classified.
    pub ambiguous_skipped: usize,
    /// Coding transcripts that could not be examined because they carry no
    /// served sequence, or `cds_start` indexes outside the served bases. These
    /// are a *different* corruption class (handled by the structural checks in
    /// [`validate_transcript_record`]) and are only counted here, not flagged.
    pub unresolved: usize,
    /// Transcripts that failed to load from the provider for a reason other
    /// than "not found" (e.g. degenerate coordinates rejected during
    /// construction). Counted, not dropped — mirroring `validate_transcripts`'
    /// `LoadError` accounting; an unknown accession (`ReferenceNotFound`) is
    /// skipped silently as genuinely out of scope.
    pub load_errors: usize,
}

impl CdsStartCodonReport {
    /// Whether the scan found any (non-allowlisted, non-alternative) start-codon
    /// inconsistency.
    pub fn has_inconsistencies(&self) -> bool {
        !self.inconsistent.is_empty()
    }
}

/// Strip a trailing `.<version>` from an accession, returning the base
/// accession (e.g. `NM_000425.3` → `NM_000425`). Returns the input unchanged
/// when there is no `.`-delimited numeric version suffix.
fn versionless_accession(accession: &str) -> &str {
    match accession.rsplit_once('.') {
        // Only treat a purely-numeric suffix as a version (RefSeq accessions
        // never embed a `.` other than the version separator, but guard anyway
        // so e.g. `LRG_199t1` is left intact).
        Some((base, ver)) if !ver.is_empty() && ver.bytes().all(|b| b.is_ascii_digit()) => base,
        _ => accession,
    }
}

/// Scan the named coding transcripts for CDS **start-codon** consistency (#629).
///
/// For every accession in `ids` that the provider can serve as a *coding*
/// transcript (CDS bounds populated), the first CDS codon — `seq[cds_start-1 ..
/// cds_start+2]` in the served FASTA, where `cds_start` is the 1-based cdot CDS
/// start — is classified. `ATG` is consistent; the recognized alternative
/// initiators `CTG`/`GTG`/`TTG` are counted as `alternative_start`; any other
/// unambiguous codon is reported as a [`CdsStartInconsistency`] unless its
/// accession (with **or** without the version suffix) appears in `allowlist`.
///
/// Accessions the provider does not have, non-coding records, ambiguous first
/// codons (`N`/IUPAC), and records whose CDS start indexes outside the served
/// sequence are not flagged — see [`CdsStartCodonReport`] for how each is
/// accounted. Internal stop codons are intentionally **not** examined (see
/// [`CdsStartCodonReport`]).
pub fn scan_cds_start_codons<P, I, S>(
    provider: &P,
    ids: I,
    allowlist: &std::collections::BTreeSet<String>,
) -> CdsStartCodonReport
where
    P: ReferenceProvider,
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    let allowlisted =
        |id: &str| allowlist.contains(id) || allowlist.contains(versionless_accession(id));

    let mut report = CdsStartCodonReport::default();
    for id in ids {
        let id = id.as_ref();
        let tx = match provider.get_transcript(id) {
            Ok(tx) => tx,
            // Unknown accession: genuinely not this scan's concern (it is the
            // file/manifest layer's job to flag a missing transcript).
            Err(FerroError::ReferenceNotFound { .. }) => continue,
            // A transcript we expected to serve (callers pass coding accessions
            // present in the FASTA) failed to construct — count it rather than
            // dropping it, mirroring `validate_transcripts`' LoadError handling.
            Err(_) => {
                report.load_errors += 1;
                continue;
            }
        };
        // Non-coding records have no CDS to validate.
        let (Some(cds_start), Some(_)) = (tx.cds_start, tx.cds_end) else {
            continue;
        };
        // `cds_start` is 1-based inclusive. Need the served bases to look at the
        // actual start codon; a coordinate-only record cannot be examined here.
        let Some(seq) = tx.sequence.as_deref() else {
            report.unresolved += 1;
            continue;
        };
        let bytes = seq.as_bytes();
        let start = cds_start as usize;
        // Guard: `cds_start >= 1` and the three-base codon must lie within the
        // served sequence. An out-of-range CDS start is a structural corruption
        // class handled elsewhere; only count it here.
        if start == 0 || start + 2 > bytes.len() {
            report.unresolved += 1;
            continue;
        }
        let codon = bytes[start - 1..start + 2].to_ascii_uppercase();
        if !is_acgt_codon(&codon) {
            report.ambiguous_skipped += 1;
            continue;
        }
        report.coding_examined += 1;
        if is_start_codon(&codon) {
            continue; // ATG — consistent
        }
        if is_alternative_start_codon(&codon) {
            // CTG/GTG/TTG: recognized near-cognate initiator; frame is intact.
            report.alternative_start += 1;
            continue;
        }
        if allowlisted(id) {
            report.allowlisted += 1;
            continue;
        }
        report.inconsistent.push(CdsStartInconsistency {
            transcript_id: id.to_string(),
            observed_start_codon: String::from_utf8_lossy(&codon).into_owned(),
        });
    }
    report
        .inconsistent
        .sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));
    report
}

/// Load a CDS start-codon allowlist from a newline-delimited text file (#629).
///
/// One accession per line; blank lines and `#`-prefixed comment lines are
/// ignored, and surrounding whitespace is trimmed. Either a versioned
/// (`NM_000425.3`) or versionless (`NM_000425`) accession is accepted —
/// [`scan_cds_start_codons`] matches against both forms, so a versionless entry
/// allows every version of that accession (alternative start codons are stable
/// across minor revisions).
pub fn load_cds_allowlist(
    path: &std::path::Path,
) -> Result<std::collections::BTreeSet<String>, FerroError> {
    let contents = std::fs::read_to_string(path).map_err(|e| FerroError::Io {
        msg: format!("failed to read CDS allowlist {}: {e}", path.display()),
    })?;
    Ok(contents
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_string)
        .collect())
}

/// Whether `codon` is the canonical start codon `ATG` (case-insensitive).
fn is_start_codon(codon: &[u8]) -> bool {
    codon.eq_ignore_ascii_case(b"ATG")
}

/// Whether `codon` is a recognized **alternative** (near-cognate) translation
/// initiation codon — `CTG`, `GTG`, or `TTG` (case-insensitive). RefSeq
/// annotates transcripts with these starts; the reading frame is correct and
/// the initiator Met is still installed, so they are legitimate non-`ATG`
/// starts, not the off-the-start drift #629 looks for. This matches the
/// canonical-start set used by the GFF loader's `W-LOAD-201` check.
fn is_alternative_start_codon(codon: &[u8]) -> bool {
    codon.eq_ignore_ascii_case(b"CTG")
        || codon.eq_ignore_ascii_case(b"GTG")
        || codon.eq_ignore_ascii_case(b"TTG")
}

/// Whether `codon` is one of the three standard stop codons.
fn is_stop_codon(codon: &[u8]) -> bool {
    codon.eq_ignore_ascii_case(b"TAA")
        || codon.eq_ignore_ascii_case(b"TAG")
        || codon.eq_ignore_ascii_case(b"TGA")
}

/// Whether every base in `codon` is an unambiguous A/C/G/T (already
/// upper-cased). Codons containing `N` or other IUPAC ambiguity codes are not
/// classifiable as start/stop and are skipped by the codon checks.
fn is_acgt_codon(codon: &[u8]) -> bool {
    codon.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::Exon;

    /// Build a coding transcript with a single exon spanning the whole
    /// sequence and the given 1-based inclusive CDS bounds.
    fn coding_tx(seq: &str, cds_start: u64, cds_end: u64) -> Transcript {
        Transcript {
            id: "NM_TEST.1".to_string(),
            sequence: Some(seq.to_string()),
            cds_start: Some(cds_start),
            cds_end: Some(cds_end),
            exons: vec![Exon::new(1, 1, seq.len() as u64)],
            ..Default::default()
        }
    }

    fn kinds(anomalies: &[TranscriptAnomaly]) -> Vec<AnomalyKind> {
        anomalies.iter().map(|a| a.kind).collect()
    }

    #[test]
    fn clean_coding_transcript_has_no_anomalies() {
        // ATG | AAA | TAA — start, one sense codon, stop. 9 nt, exon [1,9].
        assert!(validate_transcript_record(&coding_tx("ATGAAATAA", 1, 9)).is_empty());
    }

    #[test]
    fn lowercase_clean_transcript_has_no_anomalies() {
        // Soft-masked / lowercase bases must be upper-cased before codon checks.
        assert!(validate_transcript_record(&coding_tx("atgaaataa", 1, 9)).is_empty());
    }

    #[test]
    fn served_shorter_than_alignment_is_corruption() {
        // Non-coding record so only the length check runs. Exon extent 9 but
        // the served sequence is 6 nt — the alignment references missing bases.
        let tx = Transcript {
            id: "NR_TEST.1".to_string(),
            sequence: Some("ACGTAC".to_string()), // 6 nt
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 9)],
            ..Default::default()
        };
        let found = validate_transcript_record(&tx);
        assert_eq!(kinds(&found), vec![AnomalyKind::LengthMismatch]);
        assert_eq!(found[0].confidence, AnomalyConfidence::Corruption);
    }

    #[test]
    fn served_longer_than_alignment_is_not_flagged_polya() {
        // A served sequence LONGER than the exon-alignment extent is normal
        // (poly-A / unaligned 3' bases) and must NOT be flagged — only a
        // shorter-than-alignment sequence is corruption.
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.sequence = Some(format!("ATGAAATAA{}", "A".repeat(200))); // 209 nt, exon extent 9
        assert!(
            validate_transcript_record(&tx).is_empty(),
            "served-longer (poly-A) must not be a LengthMismatch"
        );
    }

    #[test]
    fn cds_not_a_triplet_is_corruption() {
        // CDS 1..8 = 8 nt, not a multiple of 3.
        let found = validate_transcript_record(&coding_tx("ATGAAATAA", 1, 8));
        assert_eq!(kinds(&found), vec![AnomalyKind::CdsNotTriplet]);
        assert_eq!(found[0].confidence, AnomalyConfidence::Corruption);
    }

    #[test]
    fn missing_start_codon_is_advisory() {
        // CDS begins with CTG, not ATG.
        let found = validate_transcript_record(&coding_tx("CTGAAATAA", 1, 9));
        assert_eq!(kinds(&found), vec![AnomalyKind::MissingStartCodon]);
        assert_eq!(found[0].confidence, AnomalyConfidence::Advisory);
    }

    #[test]
    fn missing_stop_codon_is_advisory() {
        // CDS ends with AAA, not a stop.
        let found = validate_transcript_record(&coding_tx("ATGAAAAAA", 1, 9));
        assert_eq!(kinds(&found), vec![AnomalyKind::MissingStopCodon]);
        assert_eq!(found[0].confidence, AnomalyConfidence::Advisory);
    }

    #[test]
    fn internal_stop_codon_is_advisory() {
        // ATG | TAA | TAA — in-frame stop at codon 2 (advisory: could be Sec).
        let found = validate_transcript_record(&coding_tx("ATGTAATAA", 1, 9));
        assert_eq!(kinds(&found), vec![AnomalyKind::InternalStopCodon]);
        assert_eq!(found[0].confidence, AnomalyConfidence::Advisory);
        assert!(found[0].detail.contains("codon 2"));
    }

    #[test]
    fn ambiguous_codons_are_not_flagged() {
        // ATN start and NNN internal/terminal codons are not classifiable —
        // they must not produce missing-start / missing-stop / internal-stop.
        let tx = coding_tx("ATNNNNNNN", 1, 9);
        assert!(
            validate_transcript_record(&tx).is_empty(),
            "ambiguous (N) codons must be skipped, not flagged"
        );
    }

    #[test]
    fn cds_out_of_range_is_flagged_without_panicking() {
        // CDS end past the end of the sequence.
        let found = validate_transcript_record(&coding_tx("ATGAAATAA", 1, 99));
        assert_eq!(kinds(&found), vec![AnomalyKind::CdsOutOfRange]);
        assert_eq!(found[0].confidence, AnomalyConfidence::Corruption);
    }

    #[test]
    fn half_populated_cds_bounds_are_corruption() {
        // Only (None, None) is the legitimate non-coding case. A record with
        // exactly one CDS bound set is structurally incomplete and must be
        // flagged rather than silently treated as non-coding.
        for (start, end) in [(Some(1), None), (None, Some(9))] {
            let tx = Transcript {
                id: "NM_TEST.1".to_string(),
                sequence: Some("ATGAAATAA".to_string()),
                cds_start: start,
                cds_end: end,
                exons: vec![Exon::new(1, 1, 9)],
                ..Default::default()
            };
            let found = validate_transcript_record(&tx);
            assert_eq!(
                kinds(&found),
                vec![AnomalyKind::CdsOutOfRange],
                "start={start:?} end={end:?}"
            );
            assert_eq!(found[0].confidence, AnomalyConfidence::Corruption);
        }
    }

    #[test]
    fn cds_extending_past_exon_alignment_is_corruption() {
        // The served sequence is legitimately longer than the exon-alignment
        // extent (poly-A / unaligned 3' tail, which is NOT flagged on its own),
        // but the CDS runs past the exon extent into that tail — i.e. coding
        // coordinates outside the exon alignment. That is corruption, even
        // though `cds_end` stays within the served sequence, so the plain
        // `cds_end > seq_len` check alone would miss it.
        let mut tx = coding_tx("ATGAAATAA", 1, 12); // CDS end 12 > exon extent 9
        tx.sequence = Some(format!("ATGAAATAA{}", "A".repeat(51))); // 60 nt; exon extent stays 9
        let found = validate_transcript_record(&tx);
        assert_eq!(kinds(&found), vec![AnomalyKind::CdsOutOfRange]);
        assert_eq!(found[0].confidence, AnomalyConfidence::Corruption);
    }

    #[test]
    fn cds_ending_exactly_at_exon_extent_in_polya_is_clean() {
        // Boundary guard: a CDS ending exactly at the exon extent within a much
        // longer (poly-A) served sequence must NOT be flagged — `cds_end` equal
        // to the extent is in range, only strictly past it is corruption.
        let mut tx = coding_tx("ATGAAATAA", 1, 9); // CDS end 9 == exon extent 9
        tx.sequence = Some(format!("ATGAAATAA{}", "A".repeat(200))); // 209 nt; extent 9
        assert!(
            validate_transcript_record(&tx).is_empty(),
            "CDS ending at the exon extent (boundary) must not be flagged"
        );
    }

    #[test]
    fn non_ascii_sequence_does_not_panic() {
        // A stray multi-byte char must not panic the CDS byte slice. (The
        // codon will simply not match ACGT and be skipped.)
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.sequence = Some("ATGé\u{00e9}TAA".to_string()); // multi-byte é in the middle
        tx.cds_start = Some(1);
        tx.cds_end = Some(tx.sequence.as_ref().unwrap().len() as u64);
        // Must not panic; we don't assert specific anomalies, only no panic.
        let _ = validate_transcript_record(&tx);
    }

    #[test]
    fn coordinate_only_record_is_not_an_anomaly() {
        let tx = Transcript {
            id: "NM_TEST.1".to_string(),
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            ..Default::default()
        };
        assert!(validate_transcript_record(&tx).is_empty());
    }

    #[test]
    fn noncoding_record_runs_only_the_length_check() {
        // No CDS, sequence length matches the exon extent → no anomaly.
        let tx = Transcript {
            id: "NR_TEST.1".to_string(),
            sequence: Some("ACGTACGT".to_string()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 8)],
            ..Default::default()
        };
        assert!(validate_transcript_record(&tx).is_empty());
    }

    #[test]
    fn validate_transcripts_flags_only_the_anomalous_member() {
        use crate::reference::mock::MockProvider;
        let mut provider = MockProvider::new();
        // Clean: ATG AAA TAA, exon [1,9].
        provider.add_transcript(coding_tx("ATGAAATAA", 1, 9));
        // Anomalous: internal stop at codon 2, distinct id.
        let mut bad = coding_tx("ATGTAATAA", 1, 9);
        bad.id = "NM_BAD.1".to_string();
        provider.add_transcript(bad);

        let found = validate_transcripts(
            &provider,
            &[
                "NM_TEST.1".to_string(),
                "NM_BAD.1".to_string(),
                "NM_MISSING.1".to_string(), // not present → skipped, no panic
            ],
        );
        assert_eq!(
            found.len(),
            1,
            "only the anomalous transcript should be flagged"
        );
        assert_eq!(found[0].transcript_id, "NM_BAD.1");
        assert_eq!(found[0].kind, AnomalyKind::InternalStopCodon);
    }

    /// Provider whose `get_transcript` always fails with a non-not-found error,
    /// to exercise the `LoadError` surfacing path.
    struct FailingProvider;
    impl ReferenceProvider for FailingProvider {
        fn get_transcript(&self, _id: &str) -> Result<std::sync::Arc<Transcript>, FerroError> {
            Err(FerroError::InvalidCoordinates {
                msg: "degenerate coordinates".to_string(),
            })
        }
        fn get_sequence(&self, _id: &str, _start: u64, _end: u64) -> Result<String, FerroError> {
            Err(FerroError::InvalidCoordinates {
                msg: "n/a".to_string(),
            })
        }
    }

    #[test]
    fn validate_transcripts_surfaces_non_not_found_load_errors() {
        let found = validate_transcripts(&FailingProvider, &["NM_BROKEN.1".to_string()]);
        assert_eq!(kinds(&found), vec![AnomalyKind::LoadError]);
        assert_eq!(found[0].transcript_id, "NM_BROKEN.1");
        assert_eq!(found[0].confidence, AnomalyConfidence::Corruption);
    }

    // --- authoritative comparison (phase 2) ---------------------------------

    fn auth(
        accession: &str,
        tx_length: u64,
        cds: Option<(u64, u64)>,
        protein: Option<&str>,
    ) -> AuthoritativeRecord {
        AuthoritativeRecord {
            accession: accession.to_string(),
            tx_length,
            cds_start: cds.map(|(s, _)| s),
            cds_end: cds.map(|(_, e)| e),
            protein_id: protein.map(str::to_string),
            sequence: None,
        }
    }

    #[test]
    fn authoritative_match_yields_no_anomalies() {
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.protein_id = Some("NP_X.1".to_string());
        let a = auth("NM_TEST.1", 9, Some((1, 9)), Some("NP_X.1"));
        assert!(validate_against_authoritative(&tx, &a).is_empty());
    }

    #[test]
    fn authoritative_protein_and_cds_mismatch_are_corruption() {
        // NM_012459.2 shape: served pairs the short isoform (NP_036591.3, CDS
        // 34..285); the authoritative record is NP_036591.2, CDS 31..327.
        let s = "A".repeat(327);
        let mut tx = coding_tx(&s, 34, 285);
        tx.protein_id = Some("NP_036591.3".to_string());
        let a = auth("NM_012459.2", 327, Some((31, 327)), Some("NP_036591.2"));
        let found = validate_against_authoritative(&tx, &a);
        let ks = kinds(&found);
        assert!(
            ks.contains(&AnomalyKind::AuthoritativeCdsMismatch),
            "got {ks:?}"
        );
        assert!(
            ks.contains(&AnomalyKind::AuthoritativeProteinIdMismatch),
            "got {ks:?}"
        );
        assert!(found
            .iter()
            .all(|x| x.confidence == AnomalyConfidence::Corruption));
    }

    #[test]
    fn authoritative_length_mismatch_is_flagged() {
        // NM_000193.2 shape: served 4650 nt vs authoritative 1576 nt.
        let s = "A".repeat(4650);
        let mut tx = coding_tx(&s, 342, 1730);
        tx.protein_id = Some("NP_000184.1".to_string());
        let a = auth("NM_000193.2", 1576, Some((152, 1540)), Some("NP_000184.1"));
        let ks = kinds(&validate_against_authoritative(&tx, &a));
        assert!(
            ks.contains(&AnomalyKind::AuthoritativeLengthMismatch),
            "got {ks:?}"
        );
        assert!(
            ks.contains(&AnomalyKind::AuthoritativeCdsMismatch),
            "got {ks:?}"
        );
    }

    #[test]
    fn validate_against_overrides_pulls_served_transcript_from_provider() {
        use crate::reference::mock::MockProvider;
        let mut provider = MockProvider::new();
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.protein_id = Some("NP_WRONG.9".to_string());
        provider.add_transcript(tx);

        let mut ov = CanonicalOverrides::default();
        ov.insert(auth("NM_TEST.1", 9, Some((1, 9)), Some("NP_RIGHT.1")));
        ov.insert(auth("NM_ABSENT.1", 100, Some((1, 99)), Some("NP_ABSENT.1"))); // skipped

        let found = validate_against_overrides(&provider, &ov);
        assert_eq!(
            kinds(&found),
            vec![AnomalyKind::AuthoritativeProteinIdMismatch]
        );
        assert_eq!(found[0].transcript_id, "NM_TEST.1");
    }

    #[test]
    fn missing_served_protein_or_cds_is_not_a_presence_mismatch() {
        // A served record that leaves protein_id/CDS unpopulated must NOT be
        // flagged against a coding authoritative record (a provider `None`
        // means "not populated", not "non-coding").
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.protein_id = None;
        tx.cds_start = None;
        tx.cds_end = None;
        let a = auth("NM_TEST.1", 9, Some((1, 9)), Some("NP_X.1"));
        assert!(
            validate_against_authoritative(&tx, &a).is_empty(),
            "presence-only differences must not be flagged"
        );
    }

    #[test]
    fn coordinate_only_served_record_skips_length_check() {
        // No served sequence → length check skipped despite a differing length.
        let tx = Transcript {
            id: "NM_TEST.1".to_string(),
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_X.1".to_string()),
            exons: vec![Exon::new(1, 1, 9)],
            ..Default::default()
        };
        let a = auth("NM_TEST.1", 9999, Some((1, 9)), Some("NP_X.1"));
        assert!(validate_against_authoritative(&tx, &a).is_empty());
    }

    #[test]
    fn validate_against_overrides_surfaces_non_not_found_load_errors() {
        let mut ov = CanonicalOverrides::default();
        ov.insert(auth("NM_BROKEN.1", 9, Some((1, 9)), Some("NP_X.1")));
        let found = validate_against_overrides(&FailingProvider, &ov);
        assert_eq!(kinds(&found), vec![AnomalyKind::LoadError]);
        assert_eq!(found[0].transcript_id, "NM_BROKEN.1");
    }

    #[test]
    fn authoritative_half_populated_served_cds_is_corruption() {
        // A served record with only one CDS bound set is incomplete; the
        // authoritative comparison must surface corruption rather than passing
        // clean because the four-`Some` `if let` failed to match.
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.cds_end = None;
        let a = auth("NM_TEST.1", 9, Some((1, 9)), Some("NP_X.1"));
        let ks = kinds(&validate_against_authoritative(&tx, &a));
        assert!(ks.contains(&AnomalyKind::CdsOutOfRange), "got {ks:?}");
    }

    // --- correction (phase 3) -----------------------------------------------

    #[test]
    fn correction_overrides_cds_and_protein_when_sequence_is_canonical_length() {
        // NM_012459.2 shape: canonical-length bases but cdot points at the
        // short isoform. Metadata is corrected in place.
        let s = "A".repeat(327);
        let mut tx = coding_tx(&s, 34, 285);
        tx.id = "NM_012459.2".to_string();
        tx.protein_id = Some("NP_036591.3".to_string());

        let mut ov = CanonicalOverrides::default();
        ov.insert(auth(
            "NM_012459.2",
            327,
            Some((31, 327)),
            Some("NP_036591.2"),
        ));

        let corrections = apply_canonical_overrides(&mut tx, &ov);
        assert!(corrections.contains(&Correction::CdsCorrected {
            from: (Some(34), Some(285)),
            to: (31, 327),
        }));
        assert!(corrections.contains(&Correction::ProteinIdCorrected {
            from: Some("NP_036591.3".to_string()),
            to: "NP_036591.2".to_string(),
        }));
        // Mutation applied.
        assert_eq!((tx.cds_start, tx.cds_end), (Some(31), Some(327)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_036591.2"));
    }

    #[test]
    fn correction_requires_reingestion_when_sequence_length_is_wrong() {
        // NM_000193.2 shape: served bases are the wrong length, so authoritative
        // coordinates do not apply — no override, re-ingestion reported.
        let s = "A".repeat(4650);
        let mut tx = coding_tx(&s, 342, 1730);
        tx.id = "NM_000193.2".to_string();
        tx.protein_id = Some("NP_000184.1".to_string());

        let mut ov = CanonicalOverrides::default();
        ov.insert(auth(
            "NM_000193.2",
            1576,
            Some((152, 1540)),
            Some("NP_000184.1"),
        ));

        let corrections = apply_canonical_overrides(&mut tx, &ov);
        assert_eq!(
            corrections,
            vec![Correction::SequenceReingestionRequired {
                served: 4650,
                authoritative: 1576,
            }]
        );
        // Untouched: coordinates for the canonical 1576 nt would be wrong here.
        assert_eq!((tx.cds_start, tx.cds_end), (Some(342), Some(1730)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_000184.1"));
    }

    #[test]
    fn correction_replaces_wrong_length_sequence_when_canonical_present() {
        // Wrong-sequence class WITH canonical bases carried: the served 12 nt is
        // replaced by the canonical 9 nt, then CDS/protein applied.
        let mut tx = coding_tx("ATGAAATAAGGG", 1, 12);
        tx.id = "NM_X.2".to_string();
        tx.protein_id = Some("NP_OLD.1".to_string());
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_X.2".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_NEW.2".to_string()),
            sequence: Some("ATGAAATAA".to_string()),
        });
        let corrections = apply_canonical_overrides(&mut tx, &ov);
        assert!(corrections
            .iter()
            .any(|c| matches!(c, Correction::SequenceCorrected { from: 12, to: 9 })));
        assert_eq!(tx.sequence.as_deref(), Some("ATGAAATAA"));
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_NEW.2"));
    }

    #[test]
    fn correction_replaces_same_length_edited_sequence() {
        // Same length (9) but a point edit vs canonical → replaced (byte
        // comparison, not just length).
        let mut tx = coding_tx("ATGCAATAA", 1, 9);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_TEST.1".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_X.1".to_string()),
            sequence: Some("ATGAAATAA".to_string()),
        });
        let corrections = apply_canonical_overrides(&mut tx, &ov);
        assert!(corrections
            .iter()
            .any(|c| matches!(c, Correction::SequenceCorrected { from: 9, to: 9 })));
        assert_eq!(tx.sequence.as_deref(), Some("ATGAAATAA"));
    }

    #[test]
    fn sequence_corrected_transcript_has_no_self_inflicted_length_mismatch() {
        // Served 12 nt (exons tile [1,12]); canonical is 9 nt. After correction
        // the sequence is 9 nt; the offline length check must NOT flag a
        // LengthMismatch (the exons are collapsed to the corrected length).
        let mut tx = coding_tx("ATGAAATAAGGG", 1, 9);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_TEST.1".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_X.1".to_string()),
            sequence: Some("ATGAAATAA".to_string()),
        });
        apply_canonical_overrides(&mut tx, &ov);
        assert_eq!(tx.sequence.as_deref(), Some("ATGAAATAA"));
        let anomalies = validate_transcript_record(&tx);
        assert!(
            !anomalies
                .iter()
                .any(|a| a.kind == AnomalyKind::LengthMismatch),
            "sequence-corrected transcript must not self-trip LengthMismatch: {anomalies:?}"
        );
    }

    #[test]
    fn lowercase_served_sequence_is_not_spuriously_corrected() {
        // Served bases are canonical except for case (soft-masked); this must
        // NOT be reported as a SequenceCorrected.
        let mut tx = coding_tx("atgaaataa", 1, 9);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_TEST.1".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_X.1".to_string()),
            sequence: Some("ATGAAATAA".to_string()),
        });
        let corrections = apply_canonical_overrides(&mut tx, &ov);
        assert!(
            !corrections
                .iter()
                .any(|c| matches!(c, Correction::SequenceCorrected { .. })),
            "case-only difference must not be a SequenceCorrected: {corrections:?}"
        );
    }

    #[test]
    fn correction_is_noop_when_already_canonical_or_no_override() {
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.protein_id = Some("NP_X.1".to_string());

        // No override entry for this accession.
        assert!(apply_canonical_overrides(&mut tx, &CanonicalOverrides::default()).is_empty());

        // Override that already matches → no corrections.
        let mut ov = CanonicalOverrides::default();
        ov.insert(auth("NM_TEST.1", 9, Some((1, 9)), Some("NP_X.1")));
        assert!(apply_canonical_overrides(&mut tx, &ov).is_empty());
    }

    #[test]
    fn correction_applies_to_coordinate_only_record_via_exon_sum_length() {
        // No sequence; the exon-sum length (9) matches the authoritative length.
        let mut tx = Transcript {
            id: "NM_TEST.1".to_string(),
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OLD.1".to_string()),
            exons: vec![Exon::new(1, 1, 9)],
            ..Default::default()
        };
        let mut ov = CanonicalOverrides::default();
        ov.insert(auth("NM_TEST.1", 9, Some((4, 9)), Some("NP_NEW.2")));
        let corrections = apply_canonical_overrides(&mut tx, &ov);
        assert!(corrections
            .iter()
            .any(|c| matches!(c, Correction::CdsCorrected { .. })));
        assert_eq!((tx.cds_start, tx.cds_end), (Some(4), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_NEW.2"));
    }

    #[test]
    fn correction_requires_reingestion_for_coordinate_only_length_mismatch() {
        // No sequence; the exon-sum length (9) disagrees with the authoritative
        // length (100) → re-ingestion required, record untouched.
        let mut tx = Transcript {
            id: "NM_TEST.1".to_string(),
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OLD.1".to_string()),
            exons: vec![Exon::new(1, 1, 9)],
            ..Default::default()
        };
        let mut ov = CanonicalOverrides::default();
        ov.insert(auth("NM_TEST.1", 100, Some((1, 99)), Some("NP_NEW.2")));
        let corrections = apply_canonical_overrides(&mut tx, &ov);
        assert_eq!(
            corrections,
            vec![Correction::SequenceReingestionRequired {
                served: 9,
                authoritative: 100,
            }]
        );
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
    }

    #[test]
    fn correction_skips_partial_authoritative_coding_triple() {
        // Authoritative has CDS but no protein_id → all-or-nothing: nothing is
        // overridden (a corrected CDS with a stale protein label is misleading).
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.protein_id = Some("NP_X.1".to_string());
        let mut ov = CanonicalOverrides::default();
        ov.insert(auth("NM_TEST.1", 9, Some((4, 9)), None));
        assert!(apply_canonical_overrides(&mut tx, &ov).is_empty());
        assert_eq!(
            (tx.cds_start, tx.cds_end),
            (Some(1), Some(9)),
            "no half-correction"
        );
    }

    #[test]
    fn correction_rejects_override_with_mismatched_accession() {
        // A file-backed override whose map key disagrees with its embedded
        // accession (a hand-edited / corrupt overrides file) must not rewrite
        // the served record with another accession's metadata.
        let json = r#"{
            "schema_version": 1,
            "records": {
                "NM_TEST.1": {
                    "accession": "NM_OTHER.9",
                    "tx_length": 9,
                    "cds_start": 4,
                    "cds_end": 9,
                    "protein_id": "NP_OTHER.9"
                }
            }
        }"#;
        let ov = CanonicalOverrides::from_json(json).unwrap();
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.protein_id = Some("NP_X.1".to_string());
        assert!(
            apply_canonical_overrides(&mut tx, &ov).is_empty(),
            "mismatched-accession override must be rejected"
        );
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_X.1"));
    }

    #[test]
    fn correction_rejects_out_of_range_authoritative_cds_bounds() {
        // A malformed override whose CDS bounds fall outside `1..=tx_length`
        // (here `cds_end` 99 > tx_length 9) must not be applied — installing it
        // would index past the sequence and corrupt downstream translation.
        let mut tx = coding_tx("ATGAAATAA", 1, 9);
        tx.protein_id = Some("NP_X.1".to_string());
        let mut ov = CanonicalOverrides::default();
        ov.insert(auth("NM_TEST.1", 9, Some((1, 99)), Some("NP_NEW.2")));
        assert!(
            apply_canonical_overrides(&mut tx, &ov).is_empty(),
            "out-of-range CDS override must be rejected"
        );
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_X.1"));
    }

    // --- translation comparison (phase 4) -----------------------------------

    #[test]
    fn translation_matches_canonical_protein() {
        // ATG|TGG|TAA → "MW" + stop.
        let tx = coding_tx("ATGTGGTAA", 1, 9);
        assert!(validate_translation_against_protein(&tx, "MW").is_none());
    }

    #[test]
    fn translation_tolerates_trailing_stop_in_canonical() {
        let tx = coding_tx("ATGTGGTAA", 1, 9);
        assert!(validate_translation_against_protein(&tx, "MW*").is_none());
    }

    #[test]
    fn translation_mismatch_is_advisory() {
        // Served translates to "MW" but the canonical protein is "MV" — a real
        // residue divergence (not a prefix), reported as Advisory.
        let tx = coding_tx("ATGTGGTAA", 1, 9);
        let found = validate_translation_against_protein(&tx, "MV").expect("mismatch");
        assert_eq!(found.kind, AnomalyKind::TranslationMismatch);
        assert_eq!(found.confidence, AnomalyConfidence::Advisory);
    }

    #[test]
    fn translation_suppresses_selenocysteine_readthrough() {
        // ATG|TGA|GGG|TAA: standard table stops at the in-frame TGA → "M".
        // Canonical "MUG" has Sec (U) at the recoded position → suppressed.
        let tx = coding_tx("ATGTGAGGGTAA", 1, 12);
        assert!(
            validate_translation_against_protein(&tx, "MUG").is_none(),
            "a selenoprotein readthrough must not be flagged"
        );
    }

    #[test]
    fn translation_flags_non_tga_stop_before_selenocysteine() {
        // ATG|TAA|GGG|TAA: the standard table stops at the in-frame TAA → "M",
        // the same prefix the canonical selenoprotein "MUG" produces. But only
        // TGA is recoded as selenocysteine; a TAA/TAG here is a genuine
        // nonsense edit at the Sec site and must NOT be suppressed.
        let tx = coding_tx("ATGTAAGGGTAA", 1, 12);
        let found = validate_translation_against_protein(&tx, "MUG")
            .expect("a non-TGA stop at a Sec site must be flagged, not suppressed");
        assert_eq!(found.kind, AnomalyKind::TranslationMismatch);
        assert_eq!(found.confidence, AnomalyConfidence::Advisory);
    }

    #[test]
    fn translation_forces_initiator_met_for_alternative_start() {
        // CTG|TGG|TAA translates to "LW"; forcing the initiator Met gives "MW",
        // matching the canonical protein — no false position-1 mismatch.
        let tx = coding_tx("CTGTGGTAA", 1, 9);
        assert!(validate_translation_against_protein(&tx, "MW").is_none());
    }

    #[test]
    fn translation_skipped_for_ambiguous_codons() {
        // An `N`-containing codon can't be translated → skip (not a mismatch).
        let tx = coding_tx("ATGNNGTAA", 1, 9);
        assert!(validate_translation_against_protein(&tx, "MW").is_none());
    }

    #[test]
    fn translation_skipped_for_noncoding_record() {
        let mut tx = coding_tx("ATGTGGTAA", 1, 9);
        tx.cds_start = None;
        tx.cds_end = None;
        assert!(validate_translation_against_protein(&tx, "MW").is_none());
    }

    // --- CDS start-codon consistency scan (#629) ----------------------------

    use std::collections::BTreeSet;

    fn allowlist(entries: &[&str]) -> BTreeSet<String> {
        entries.iter().map(|s| s.to_string()).collect()
    }

    /// Provider seeded with three coding transcripts: one ATG-clean, one with a
    /// non-ATG start (the #629 inconsistency), one non-coding, plus an
    /// ambiguous-start coding transcript.
    fn scan_provider() -> crate::reference::mock::MockProvider {
        use crate::reference::mock::MockProvider;
        let mut p = MockProvider::new();

        let mut clean = coding_tx("ATGAAATAA", 1, 9);
        clean.id = "NM_000001.1".to_string();
        p.add_transcript(clean);

        // cds_start lands on AGG, not ATG — the L1CAM/TPM3-shaped #629 case.
        let mut bad = coding_tx("AGGAAATAA", 1, 9);
        bad.id = "NM_000425.3".to_string();
        p.add_transcript(bad);

        let mut noncoding = coding_tx("ACGTACGT", 1, 9);
        noncoding.id = "NR_000002.1".to_string();
        noncoding.cds_start = None;
        noncoding.cds_end = None;
        p.add_transcript(noncoding);

        // First codon contains an ambiguity code → unclassifiable, must skip.
        let mut ambiguous = coding_tx("ANGAAATAA", 1, 9);
        ambiguous.id = "NM_000003.1".to_string();
        p.add_transcript(ambiguous);

        // Recognized alternative initiation codon (CTG) — legitimate, must not
        // be flagged as an inconsistency.
        let mut alt_start = coding_tx("CTGAAATAA", 1, 9);
        alt_start.id = "NM_000004.1".to_string();
        p.add_transcript(alt_start);

        p
    }

    #[test]
    fn scan_flags_only_the_non_initiator_coding_transcript() {
        let provider = scan_provider();
        let ids = [
            "NM_000001.1".to_string(),
            "NM_000425.3".to_string(),
            "NR_000002.1".to_string(),
            "NM_000003.1".to_string(),
            "NM_000004.1".to_string(),
            "NM_ABSENT.1".to_string(), // provider doesn't have it → skipped
        ];
        let report = scan_cds_start_codons(&provider, &ids, &BTreeSet::new());

        // Examined = the unambiguous coding transcripts (clean ATG + bad AGG +
        // alt CTG); the non-coding and ambiguous ones do not count.
        assert_eq!(report.coding_examined, 3);
        assert_eq!(report.alternative_start, 1); // CTG
        assert_eq!(report.ambiguous_skipped, 1);
        assert_eq!(report.allowlisted, 0);
        assert_eq!(report.inconsistent.len(), 1);
        assert_eq!(report.inconsistent[0].transcript_id, "NM_000425.3");
        assert_eq!(report.inconsistent[0].observed_start_codon, "AGG");
        assert!(report.has_inconsistencies());
    }

    #[test]
    fn scan_recognized_alternative_start_codons_are_not_inconsistent() {
        use crate::reference::mock::MockProvider;
        let mut p = MockProvider::new();
        for (i, codon) in ["CTG", "GTG", "TTG"].iter().enumerate() {
            let mut tx = coding_tx(&format!("{codon}AAATAA"), 1, 9);
            tx.id = format!("NM_00010{i}.1");
            p.add_transcript(tx);
        }
        let ids = ["NM_000100.1", "NM_000101.1", "NM_000102.1"].map(String::from);
        let report = scan_cds_start_codons(&p, &ids, &BTreeSet::new());
        assert_eq!(report.coding_examined, 3);
        assert_eq!(report.alternative_start, 3);
        assert!(!report.has_inconsistencies());
    }

    #[test]
    fn scan_lowercase_start_codon_is_consistent() {
        use crate::reference::mock::MockProvider;
        let mut p = MockProvider::new();
        let mut soft_masked = coding_tx("atgaaataa", 1, 9);
        soft_masked.id = "NM_000009.1".to_string();
        p.add_transcript(soft_masked);
        let report = scan_cds_start_codons(&p, ["NM_000009.1".to_string()], &BTreeSet::new());
        assert_eq!(report.coding_examined, 1);
        assert!(!report.has_inconsistencies());
    }

    #[test]
    fn scan_allowlist_suppresses_by_exact_and_versionless_accession() {
        let provider = scan_provider();
        let ids = ["NM_000425.3".to_string()];

        // Exact, versioned match.
        let exact = scan_cds_start_codons(&provider, &ids, &allowlist(&["NM_000425.3"]));
        assert!(!exact.has_inconsistencies());
        assert_eq!(exact.allowlisted, 1);

        // Versionless match allows every version of the accession.
        let versionless = scan_cds_start_codons(&provider, &ids, &allowlist(&["NM_000425"]));
        assert!(!versionless.has_inconsistencies());
        assert_eq!(versionless.allowlisted, 1);

        // A different accession on the allowlist does not suppress it.
        let other = scan_cds_start_codons(&provider, &ids, &allowlist(&["NM_999999"]));
        assert!(other.has_inconsistencies());
    }

    #[test]
    fn scan_counts_out_of_range_cds_start_as_unresolved_not_inconsistent() {
        use crate::reference::mock::MockProvider;
        let mut p = MockProvider::new();
        // cds_start past the end of the served sequence: a structural class
        // handled elsewhere; the start-codon scan only counts it as unresolved.
        let mut oor = coding_tx("ATGAAATAA", 1, 9);
        oor.id = "NM_000010.1".to_string();
        oor.cds_start = Some(99);
        p.add_transcript(oor);
        let report = scan_cds_start_codons(&p, ["NM_000010.1".to_string()], &BTreeSet::new());
        assert_eq!(report.unresolved, 1);
        assert_eq!(report.coding_examined, 0);
        assert!(!report.has_inconsistencies());
    }

    #[test]
    fn scan_counts_coordinate_only_coding_record_as_unresolved() {
        use crate::reference::mock::MockProvider;
        let mut p = MockProvider::new();
        // Coding (CDS bounds set) but no served sequence: cannot examine the
        // start codon, so it must be counted as unresolved, not flagged.
        let mut coord_only = coding_tx("ATGAAATAA", 1, 9);
        coord_only.id = "NM_000011.1".to_string();
        coord_only.sequence = None;
        p.add_transcript(coord_only);
        let report = scan_cds_start_codons(&p, ["NM_000011.1".to_string()], &BTreeSet::new());
        assert_eq!(report.unresolved, 1);
        assert_eq!(report.coding_examined, 0);
        assert!(!report.has_inconsistencies());
    }

    #[test]
    fn scan_counts_non_not_found_load_errors() {
        // A provider error other than ReferenceNotFound must be counted, not
        // silently skipped (mirrors `validate_transcripts`' LoadError path).
        let report = scan_cds_start_codons(
            &FailingProvider,
            ["NM_BROKEN.1".to_string()],
            &BTreeSet::new(),
        );
        assert_eq!(report.load_errors, 1);
        assert_eq!(report.coding_examined, 0);
        assert!(!report.has_inconsistencies());
    }

    #[test]
    fn scan_inconsistent_list_is_sorted_by_accession() {
        use crate::reference::mock::MockProvider;
        let mut p = MockProvider::new();
        for (idx, id) in ["NM_000300.1", "NM_000100.1", "NM_000200.1"]
            .iter()
            .enumerate()
        {
            let mut tx = coding_tx("AGGAAATAA", 1, 9); // all non-ATG
            tx.id = id.to_string();
            // vary nothing else; idx only to avoid unused warnings
            let _ = idx;
            p.add_transcript(tx);
        }
        let ids = [
            "NM_000300.1".to_string(),
            "NM_000100.1".to_string(),
            "NM_000200.1".to_string(),
        ];
        let report = scan_cds_start_codons(&p, &ids, &BTreeSet::new());
        let got: Vec<&str> = report
            .inconsistent
            .iter()
            .map(|i| i.transcript_id.as_str())
            .collect();
        assert_eq!(got, vec!["NM_000100.1", "NM_000200.1", "NM_000300.1"]);
    }

    #[test]
    fn versionless_accession_strips_only_numeric_version() {
        assert_eq!(versionless_accession("NM_000425.3"), "NM_000425");
        assert_eq!(versionless_accession("NM_000425"), "NM_000425");
        // A non-numeric suffix (e.g. an LRG transcript tag) is left intact.
        assert_eq!(versionless_accession("LRG_199t1"), "LRG_199t1");
    }

    #[test]
    fn load_cds_allowlist_ignores_comments_and_blanks() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "# legitimate alternative start codons").unwrap();
        writeln!(f, "NM_000425").unwrap();
        writeln!(f).unwrap();
        writeln!(f, "  NM_152263.2  ").unwrap();
        f.flush().unwrap();
        let set = load_cds_allowlist(f.path()).unwrap();
        assert_eq!(set.len(), 2);
        assert!(set.contains("NM_000425"));
        assert!(set.contains("NM_152263.2"));
    }
}
