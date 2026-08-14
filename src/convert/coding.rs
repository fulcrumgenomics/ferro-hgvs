//! CDS (c.) coordinate handling
//!
//! # Coordinate System
//!
//! | Position Type | Basis | Notes |
//! |---------------|-------|-------|
//! | `CdsPos.base` | 1-based | Positive for CDS, negative for 5'UTR |
//! | `cds_start`, `cds_end` | 1-based | From transcript metadata |
//! | Return values | 0-based | For array indexing |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::error::FerroError;
use crate::hgvs::location::CdsPos;
use crate::reference::transcript::Transcript;

/// Validate a CDS position against a transcript
pub fn validate_cds_pos(pos: &CdsPos, transcript: &Transcript) -> Result<(), FerroError> {
    let cds_start = transcript
        .cds_start
        .ok_or_else(|| FerroError::ConversionError {
            msg: "Transcript has no CDS".to_string(),
        })?;
    let cds_end = transcript
        .cds_end
        .ok_or_else(|| FerroError::ConversionError {
            msg: "Transcript has no CDS end".to_string(),
        })?;

    let cds_length = (cds_end - cds_start + 1) as i64;

    if pos.utr3 {
        // 3' UTR position
        let utr3_length = transcript.sequence_length() - cds_end;
        if pos.base < 1 || pos.base > utr3_length as i64 {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "3' UTR position *{} is out of range (max *{})",
                    pos.base, utr3_length
                ),
            });
        }
    } else if pos.base < 0 {
        // 5' UTR position
        let utr5_length = (cds_start - 1) as i64;
        // `unsigned_abs` for the same reason as the offset checks below:
        // `i64::MIN.abs()` overflows, and `CdsPos`'s fields are public, so a
        // library caller can hand this a magnitude no parse can produce
        // (`parse_cds_pos` negates a `digit1`, whose minimum is
        // `-i64::MAX`). `utr5_length` is non-negative by construction, so
        // `max(0)` only makes the `u64` comparison well-typed (#1767).
        if pos.base.unsigned_abs() > utr5_length.max(0) as u64 {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "5' UTR position {} is out of range (min -{})",
                    pos.base, utr5_length
                ),
            });
        }
    } else if pos.base > cds_length {
        // CDS position out of range
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "CDS position {} is out of range (max {})",
                pos.base, cds_length
            ),
        });
    }

    // An unknown offset (`c.100+?` / `c.100-?`) is not a distance, so it is not
    // a coordinate this can validate — say so rather than reporting it as an
    // out-of-range number (#1767).
    if pos.has_unknown_offset() {
        return Err(FerroError::InvalidCoordinates {
            msg: "Unknown intronic offset (`+?` / `-?`) denotes no transcript coordinate"
                .to_string(),
        });
    }

    // Validate intronic offset if present
    if let Some(offset) = pos.offset {
        // `unsigned_abs` because `i64::MIN.abs()` overflows (#1767).
        if offset.unsigned_abs() > 1_000_000 {
            // Sanity check for very large offsets
            return Err(FerroError::InvalidCoordinates {
                msg: format!("Intronic offset {} is unreasonably large", offset),
            });
        }
    }

    Ok(())
}

/// Convert a CDS position to a 0-based transcript position
pub fn cds_to_transcript_pos(pos: &CdsPos, transcript: &Transcript) -> Result<u64, FerroError> {
    validate_cds_pos(pos, transcript)?;

    let cds_start = transcript.cds_start.unwrap(); // Already validated
    let cds_end = transcript.cds_end.unwrap();

    let tx_pos = if pos.utr3 {
        cds_end + pos.base as u64
    } else if pos.base < 0 {
        // 5'UTR: HGVS numbering skips c.0 (c.-1 is the base immediately
        // upstream of c.1), so c.-N maps to tx position cds_start - N.
        // Issue #97.
        (cds_start as i64 + pos.base) as u64
    } else if pos.base == 0 {
        // c.0 is not a valid HGVS position; preserve the legacy mapping
        // (last 5'UTR base, equivalent to c.-1) rather than failing.
        cds_start.saturating_sub(1)
    } else {
        cds_start + pos.base as u64 - 1
    };

    Ok(tx_pos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            cds_start_incomplete: false,
            id: "NM_TEST.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: Some("AAAAATGCCCAAAGGGTTTTAAAAAA".to_string()), // 26 bases
            cds_start: Some(6),
            cds_end: Some(20),
            exons: vec![Exon::new(1, 1, 26)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_cds_to_transcript_pos_5utr_and_legacy_c0() {
        // Pin the issue #97 fix and the legacy c.0 mapping in the public
        // helper, mirroring Normalizer::cds_to_tx_pos.
        //
        // 5'UTR (issue #97): HGVS skips c.0, so c.-N maps to tx position
        // cds_start - N — i.e. (cds_start as i64 + pos.base) as u64.
        let tx = make_test_transcript();
        let cds_start = tx.cds_start.unwrap(); // 6

        let result = cds_to_transcript_pos(&CdsPos::new(-1), &tx).unwrap();
        assert_eq!(result, (cds_start as i64 + (-1)) as u64);

        // Legacy c.0 mapping: treat as the last 5'UTR base (== c.-1) via
        // cds_start.saturating_sub(1) rather than failing.
        let result = cds_to_transcript_pos(&CdsPos::new(0), &tx).unwrap();
        assert_eq!(result, cds_start.saturating_sub(1));

        // Defensive: cds_start == 0 must not underflow. saturating_sub
        // pins the result at 0 instead of wrapping to u64::MAX.
        let tx_zero_cds = Transcript {
            cds_start_incomplete: false,
            id: "NM_TEST_ZERO.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: Some("ATGCCCAAAGGG".to_string()),
            cds_start: Some(0),
            cds_end: Some(5),
            exons: vec![Exon::new(1, 1, 12)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        let result = cds_to_transcript_pos(&CdsPos::new(0), &tx_zero_cds).unwrap();
        assert_eq!(result, tx_zero_cds.cds_start.unwrap().saturating_sub(1));
        assert_eq!(result, 0);
    }

    /// Boundaries of `validate_cds_pos`, one assertion per edge: both sides
    /// of the 3'UTR guard (`*0` err, `*1` ok, `*6` ok, `*7` err), both sides
    /// of the 5'UTR guard (`-5` ok, `-6` err), both sides of the CDS guard
    /// (`15` ok, `16` err), and the offset sanity limit on both signs — the
    /// guard is `offset.unsigned_abs() > 1_000_000`, pinned here at the
    /// ±1_000_000 (ok) and ±1_000_001 (err) pairs, so dropping the
    /// magnitude check would otherwise sail through unpinned on the
    /// negative side.
    ///
    /// One offset class is deliberately NOT pinned here: the uncertain-offset
    /// sentinels `-?` / `+?` (`OFFSET_UNKNOWN_NEGATIVE` /
    /// `OFFSET_UNKNOWN_POSITIVE` = `i64::MIN` / `i64::MAX`, declared in
    /// `src/hgvs/parser/position.rs`).
    ///
    /// This paragraph previously recorded that class as an OPEN production
    /// gap — that `i64::MIN.abs()` overflows, so `validate_cds_pos` panicked
    /// in debug on `c.100-?` instead of returning an error (issue #1087).
    /// That is no longer true: `validate_cds_pos` declines the sentinels up
    /// front via `pos.has_unknown_offset()`, and its magnitude checks use
    /// `unsigned_abs()` throughout, so `c.100-?` and `c.100+?` both return
    /// `Err` and neither overflows. That class is pinned by
    /// `unknown_offset_sentinels_are_rejected` below.
    #[test]
    fn validate_cds_pos_boundaries() {
        let tx = make_test_transcript();

        assert!(
            validate_cds_pos(&CdsPos::new(1), &tx).is_ok(),
            "1 is the first CDS base"
        );
        assert!(
            validate_cds_pos(&CdsPos::new(20), &tx).is_err(),
            "20 is out of CDS range"
        );
        assert!(
            validate_cds_pos(&CdsPos::new(-3), &tx).is_ok(),
            "-3 is within the 5'UTR"
        );
        assert!(
            validate_cds_pos(&CdsPos::new(-10), &tx).is_err(),
            "-10 is out of 5'UTR range"
        );
        assert!(
            validate_cds_pos(&CdsPos::utr3(3), &tx).is_ok(),
            "*3 is within the 3'UTR"
        );
        assert!(
            validate_cds_pos(&CdsPos::utr3(10), &tx).is_err(),
            "*10 is out of 3'UTR range"
        );

        // 3'UTR runs *1..*6 inclusive.
        assert!(
            validate_cds_pos(&CdsPos::utr3(0), &tx).is_err(),
            "*0 is not a valid 3'UTR position"
        );
        assert!(
            validate_cds_pos(&CdsPos::utr3(1), &tx).is_ok(),
            "*1 is the first 3'UTR base"
        );
        assert!(
            validate_cds_pos(&CdsPos::utr3(6), &tx).is_ok(),
            "*6 is the last 3'UTR base"
        );
        assert!(
            validate_cds_pos(&CdsPos::utr3(7), &tx).is_err(),
            "*7 is past the 3'UTR"
        );

        // 5'UTR runs -1..-5 inclusive.
        assert!(
            validate_cds_pos(&CdsPos::new(-5), &tx).is_ok(),
            "-5 is the last 5'UTR base"
        );
        assert!(
            validate_cds_pos(&CdsPos::new(-6), &tx).is_err(),
            "-6 is past the 5'UTR"
        );

        // CDS runs 1..15 inclusive.
        assert!(
            validate_cds_pos(&CdsPos::new(15), &tx).is_ok(),
            "15 is the last CDS base"
        );
        assert!(
            validate_cds_pos(&CdsPos::new(16), &tx).is_err(),
            "16 is past the CDS"
        );

        // The offset sanity limit is `offset.unsigned_abs() > 1_000_000`, so
        // 1_000_000 is legal on both the positive and negative side.
        let at_limit = CdsPos::with_offset(1, 1_000_000);
        assert!(
            validate_cds_pos(&at_limit, &tx).is_ok(),
            "offset 1_000_000 is at the limit"
        );

        let past_limit = CdsPos::with_offset(1, 1_000_001);
        assert!(
            validate_cds_pos(&past_limit, &tx).is_err(),
            "offset 1_000_001 is past it"
        );

        let at_negative_limit = CdsPos::with_offset(1, -1_000_000);
        assert!(
            validate_cds_pos(&at_negative_limit, &tx).is_ok(),
            "offset -1_000_000 is at the limit"
        );

        let past_negative_limit = CdsPos::with_offset(1, -1_000_001);
        assert!(
            validate_cds_pos(&past_negative_limit, &tx).is_err(),
            "offset -1_000_001 is past it"
        );
    }

    /// Exact tx positions for each branch of `cds_to_transcript_pos`.
    #[test]
    fn cds_to_transcript_pos_arithmetic() {
        let tx = make_test_transcript();

        // CDS: cds_start + base - 1.
        assert_eq!(cds_to_transcript_pos(&CdsPos::new(1), &tx).unwrap(), 6);
        assert_eq!(cds_to_transcript_pos(&CdsPos::new(2), &tx).unwrap(), 7);
        assert_eq!(cds_to_transcript_pos(&CdsPos::new(15), &tx).unwrap(), 20);

        // 3'UTR: cds_end + base.
        assert_eq!(cds_to_transcript_pos(&CdsPos::utr3(1), &tx).unwrap(), 21);
        assert_eq!(cds_to_transcript_pos(&CdsPos::utr3(3), &tx).unwrap(), 23);
    }

    /// The uncertain-offset sentinels `-?` / `+?` (`OFFSET_UNKNOWN_NEGATIVE` /
    /// `OFFSET_UNKNOWN_POSITIVE` = `i64::MIN` / `i64::MAX`) are not distances,
    /// so `validate_cds_pos` must decline them rather than measure them.
    ///
    /// Two failure modes are pinned at once:
    ///
    /// 1. **Panic.** `i64::MIN.abs()` overflows, so reaching this with a plain
    ///    `.abs()` panics in debug on `c.1-?`. Reaching the assertion at all
    ///    proves the magnitude check uses `unsigned_abs()`.
    /// 2. **Refusal for the wrong reason.** This is the one that makes the
    ///    assertion below check the MESSAGE rather than just `is_err()`.
    ///    Both sentinels have enormous magnitudes, so the
    ///    `unsigned_abs() > 1_000_000` check rejects them on its own — which
    ///    means deleting the `has_unknown_offset()` decline entirely leaves
    ///    an `is_err()`-only assertion **green**. Verified by sabotage: with
    ///    the guard removed, an `is_err()` version of this test still passed.
    ///    A sentinel refused as "unreasonably large" has been silently
    ///    reclassified as a measured distance that happens to be too big,
    ///    which is the opposite of what #1767 decided.
    #[test]
    fn unknown_offset_sentinels_are_rejected() {
        use crate::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};

        let tx = make_test_transcript();

        for (label, offset) in [
            ("c.1-?", OFFSET_UNKNOWN_NEGATIVE),
            ("c.1+?", OFFSET_UNKNOWN_POSITIVE),
        ] {
            let err = validate_cds_pos(&CdsPos::with_offset(1, offset), &tx)
                .expect_err("must be declined, not measured");
            let msg = err.to_string();
            assert!(
                msg.contains("denotes no transcript coordinate"),
                "{label} must be refused AS A SENTINEL; got: {msg}"
            );
            assert!(
                !msg.contains("unreasonably large"),
                "{label} must not be refused as a mere magnitude; got: {msg}"
            );
        }
    }
}
