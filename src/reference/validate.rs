//! Offline structural validation of transcript reference records.
//!
//! These checks are **self-contained**: they look only at a single
//! [`Transcript`]'s own sequence, CDS coordinates, and exon alignment, with no
//! network access and no authoritative external record. They catch a class of
//! reference-data corruption where the served bases and the alignment/CDS
//! metadata disagree, or the CDS does not translate as a sane open reading
//! frame.
//!
//! Scope (issue #520, phase 0): this is the cheap, offline tier. Anomalies are
//! tagged with an [`AnomalyConfidence`]:
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
//! Deliberately out of scope here:
//!
//! - A served sequence *longer* than the alignment is **not** flagged — a
//!   poly-A tail or unaligned 3' bases routinely make the transcript FASTA
//!   longer than the genome-derived exon alignment. Detecting an
//!   implausibly-long sequence precisely (e.g. issue #505's `NM_000193.2`,
//!   served ~3× too long) needs the authoritative transcript length, which is
//!   the authoritative-comparison phase, not this offline tier.
//! - A sequence that is internally consistent but non-canonical or mis-paired
//!   to the wrong protein version (e.g. `NM_012459.2`) — also a later phase,
//!   since it requires the authoritative RefSeq record for the exact version.

use crate::error::FerroError;
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::Transcript;

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
    /// The CDS coordinates extend beyond the available sequence (or are
    /// otherwise degenerate).
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
    if let Some(exon_extent) = tx.exons.iter().map(|e| e.end).max() {
        if seq_len < exon_extent {
            anomalies.push(TranscriptAnomaly::new(
                &tx.id,
                AnomalyKind::LengthMismatch,
                AnomalyConfidence::Corruption,
                format!(
                    "served sequence is {seq_len} nt, shorter than the {exon_extent} nt \
                     spanned by the exon alignment (bases are missing)"
                ),
            ));
        }
    }

    // --- CDS open-reading-frame sanity --------------------------------------
    let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) else {
        // Non-coding (no CDS): only the length check applies.
        return anomalies;
    };

    // `cds_start`/`cds_end` are 1-based inclusive transcript coordinates.
    if cds_start == 0 || cds_end < cds_start || cds_end > seq_len {
        anomalies.push(TranscriptAnomaly::new(
            &tx.id,
            AnomalyKind::CdsOutOfRange,
            AnomalyConfidence::Corruption,
            format!("CDS {cds_start}..{cds_end} (1-based) is invalid for a {seq_len} nt sequence"),
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

/// Whether `codon` is the canonical start codon `ATG` (case-insensitive).
fn is_start_codon(codon: &[u8]) -> bool {
    codon.eq_ignore_ascii_case(b"ATG")
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
        fn get_transcript(&self, _id: &str) -> Result<Transcript, FerroError> {
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
}
