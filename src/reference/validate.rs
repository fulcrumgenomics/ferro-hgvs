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

use crate::error::FerroError;
use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
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
    let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) else {
        // Non-coding (no CDS): only the length check applies.
        return anomalies;
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
    if let (Some(ts), Some(te), Some(as_), Some(ae)) =
        (tx.cds_start, tx.cds_end, auth.cds_start, auth.cds_end)
    {
        if (ts, te) != (as_, ae) {
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
}
