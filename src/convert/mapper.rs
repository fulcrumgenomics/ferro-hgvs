//! Coordinate mapper for converting between coordinate systems
//!
//! # Coordinate Systems
//!
//! This module handles conversions between multiple coordinate systems:
//!
//! | System | Basis | Notes |
//! |--------|-------|-------|
//! | Genomic | 1-based | Closed `[start, end]` exon intervals ([`crate::reference::transcript::Exon`]); `GenomicLocation.position` is an HGVS `g.` number |
//! | Transcript (tx) | 1-based | `TxPos.base` is 1-based |
//! | CDS (c.) | 1-based | `CdsPos.base` is 1-based, negative for 5'UTR |
//! | Protein (p.) | 1-based | `ProtPos.number` is 1-based |
//!
//! The genomic row is 1-based **inclusive on both ends** because this mapper only
//! ever reads a materialised [`crate::reference::transcript::Transcript`], never raw
//! cdot JSON: `genomic_to_tx` tests `g_start <= pos <= g_end` against
//! `Exon::genomic_start`/`genomic_end`, which are 1-based inclusive by construction
//! (whichever loader built them has already converted). Raw cdot's 0-based half-open
//! genomic bounds are a different layer entirely and never reach this module — see
//! [`crate::coords`] for that layer and [`crate::data::cdot`] for the in-memory one.
//!
//! ## Key conversions:
//! - CDS → Tx: `cds_to_tx()` - accounts for CDS start/end and UTR regions
//! - Tx → CDS: `tx_to_cds()` - inverse of above
//! - Genomic → Tx: `genomic_to_tx()` - uses exon boundaries
//! - Tx → Genomic: `tx_to_genomic()` - accounts for strand
//! - CDS → Protein: `cds_to_protein()` - divides by 3 for codon position
//!
//! ## Intronic positions:
//! CDS positions can have intronic offsets (e.g., c.100+5, c.200-10).
//! Use `*_with_intron()` methods to handle these cases.
//!
//! # The two frames, and which one `c.`/`n.` lives in
//!
//! The conversions above split into two frames that must not borrow each
//! other's arithmetic. Getting that wrong is #1619.
//!
//! **The sequence frame — `cds_to_tx` / `tx_to_cds`.** A `c.` or `n.` position
//! names a base of the transcript reference sequence and nothing else.
//! `docs/background/numbering.md:21` anchors `c.1` on "the **`A`** of the `ATG`
//! translation initiation (start) codon"; `:52` numbers `n.` "from the first to
//! the last nucleotide of the reference sequence"; and `:40` states that a
//! coding DNA reference sequence "**does not contain** intron or 5' and 3' gene
//! flanking sequences". So the axis is a plain count of the transcript's own
//! bases, and `c.N` is `cds_start + N - 1` on the flat sequence — the same
//! bytes `ReferenceProvider::get_sequence` serves and the normalizer indexes.
//! **These two methods therefore do not read the exon table at all.**
//!
//! Two shapes have each tried to make them read it, and both were wrong:
//!
//! - a CIGAR insertion, i.e. transcript bases absent from the genome. `c.`
//!   counts them, so no adjustment (#944 — `NM_015120.4` `c.87` is `n.198`,
//!   not `n.201`).
//! - a transcript-coordinate gap, i.e. transcript bases that align to no exon.
//!   `c.` counts those too (#1619 — `NM_033517.1` carries a real 39-base cdot
//!   hole, and RefSeq annotates its CDS as the contiguous `1..5196` across it,
//!   matching the 1731-aa `NP_277052.1`). Walking the exon list moved every
//!   `c.` position on that accession, so `hgvs_to_spdi` and the normalizer
//!   named different bases of the same transcript.
//!
//! **The genome frame — `genomic_to_tx` / `tx_to_genomic` and everything built
//! on them.** Mapping between a contig and a transcript genuinely is exon- and
//! CIGAR-aware, and stays so. That is where the exon table belongs, and where a
//! gap correctly means "this transcript base aligns to nothing" rather than
//! "renumber the positions after it".
//!
//! The composed conversions (`genomic_to_cds`, `cds_to_genomic`) are one step
//! of each: exon-aware across genome↔tx, flat across tx↔CDS.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::error::FerroError;
use crate::hgvs::location::{AminoAcid, CdsPos, ProtPos, TxPos, CDS_BASE_UNKNOWN};
use crate::reference::transcript::{IntronPosition, Strand, Transcript};

/// Maps coordinates between different systems for a transcript
pub struct CoordinateMapper<'a> {
    transcript: &'a Transcript,
}

impl<'a> CoordinateMapper<'a> {
    /// Create a new mapper for a transcript
    pub fn new(transcript: &'a Transcript) -> Self {
        Self { transcript }
    }

    /// Convert CDS position to transcript position
    ///
    /// CDS positions are relative to the start of the coding sequence,
    /// while transcript positions are absolute on the transcript.
    ///
    /// **This is a sequence-frame conversion and it does not consult the exon
    /// table** — see the module's "two frames" note. `c.`/`n.` numbering counts
    /// the transcript's own bases, so
    /// `c.N` is `cds_start + N - 1` on the flat transcript, with no adjustment
    /// for exon geometry, for a CIGAR insertion (#944), or for a
    /// transcript-coordinate gap in the genome alignment (#1619).
    ///
    /// One input where this is **not** interchangeable with the normalizer's
    /// own sequence-frame copy (`Normalizer::cds_to_tx_pos`): a `c.-N` further
    /// upstream than the transcript is long. This method returns a `TxPos` with
    /// a *negative* `base` (which [`Self::fold_non_intronic_utr_offset`] relies
    /// on), where the normalizer's copy errors.
    ///
    /// `c.0` used to be a second such input, and is the divergence #1772 was
    /// opened on: it fell into the `base < 1` branch and answered `cds_start`,
    /// where the normalizer continued the formula and answered `cds_start - 1`.
    /// It no longer reaches that branch. `CDS_BASE_UNKNOWN` is the integer `0`,
    /// so a `base == 0` anchor — whether the `c.?` sentinel, a `*0`, a
    /// terminus marker, or one cis-collapse built — is refused outright by the
    /// no-base guard below rather than computed on. #1772's question is
    /// therefore moot *here*; it is still live for the normalizer's own copy.
    ///
    /// # Errors
    ///
    /// Beyond a transcript with no CDS, this declines two kinds of input that
    /// are not coordinates at all rather than computing on them:
    ///
    /// - a base of `0`, which names no nucleotide on this 1-based axis. That
    ///   integer is shared by three distinct shapes — the `c.?` unknown
    ///   sentinel, a `c.*0` in the 3'UTR zone, and a `pter`/`qter`/`cen`
    ///   terminus marker — and all three are refused, together with a marker
    ///   carrying any base at all;
    /// - a base whose conversion **overflows** `i64`, which is #1690.
    pub fn cds_to_tx(&self, pos: &CdsPos) -> Result<TxPos, FerroError> {
        // A CDS base of 0 names no nucleotide. `background/numbering.md:31`
        // states it inside the definition of the `c.` axis itself, immediately
        // after the `-1` and `*1` lines: "there is no nucleotide `c.0`". That
        // is an existence claim, so there is no transcript position to answer
        // with, and the doctrine is already written on the sibling predicate
        // `CdsPos::has_unknown_offset` for an unknown *offset*: the sentinel is
        // not a distance, and "no coordinate can be derived from them at all".
        //
        // The `*` zone is settled by the line above it, `:30`: nucleotides 3'
        // of the stop "are marked with a `*` (asterisk) and numbered `c.*1`,
        // `c.*2`, `c.*3`, etc., going further downstream". The enumeration
        // starts at `*1`, so `*0` is outside it — and the coordinate it was
        // being answered with, `cds_end`, is the last CDS base, which already
        // has a plain `c.` spelling.
        //
        // **Keyed on `base == 0`, not on `is_unknown()`, and that is the whole
        // point of this guard.** `is_unknown()` is
        // `base == CDS_BASE_UNKNOWN && !utr3 && special.is_none()`, so it
        // recognises only ONE of the three shapes that carry the integer 0, and
        // the two it misses are the more reachable ones:
        //
        //   c.?              base 0, utr3=false, special=None   is_unknown
        //   c.*0             base 0, utr3=TRUE                  NOT is_unknown
        //   c.pter/qter/cen  base 0, special=Some               NOT is_unknown
        //
        // Measured on `NM_TEST.1` (cds_start 6, cds_end 35) before this widened:
        // `c.*0` answered tx 35, which is `c.30`'s own answer, and each of
        // `c.pter`, `c.qter` and `c.cen` answered tx 6, which is `c.1`'s. That
        // collapse survives to the entry point: through `hgvs_to_spdi` on
        // `JsonProvider::with_test_data`, `c.pterdel`, `c.qterdel`, `c.cendel`
        // and `c.1del` all produced the bit-identical triple
        // `NM_001234.1:4:A:`. `pter` and `qter` are opposite ends of a
        // chromosome, so answering them the same coordinate — and the same
        // coordinate as an ordinary `c.1` — is not a near miss.
        //
        // This is also what `Normalizer::cds_to_tx_pos` already says it does
        // (#1799, `src/normalize/mod.rs`): "Base 0 names no nucleotide whatever
        // put it there, so one predicate covers both", citing this method as
        // the mapper side of the same seam. Until now that citation was false.
        //
        // A terminus marker is refused on its SHAPE as well as on its base. It
        // is a genomic landmark, not a displacement from the CDS start, so no
        // `c.` arithmetic applies to it whatever integer accompanies it — the
        // constructors all set 0 today, and keying on the shape too means a
        // later one that does not cannot silently re-open this. The markers are
        // resolved, where they can be resolved at all, by
        // `project_cds_terminus_to_parent` against a named genomic parent,
        // which is a different frame from this one.
        if pos.base == CDS_BASE_UNKNOWN || pos.special.is_some() {
            // The rendering usually does NOT read `c.0`, and the message says
            // so rather than asserting a spelling the reader cannot see:
            // `CDS_BASE_UNKNOWN` *is* 0, so the unknown sentinel prints `c.?`,
            // a marker prints as the marker, and only `c.*0` prints its zero.
            let shape = if pos.special.is_some() {
                "it is a telomere/centromere marker, which names a genomic landmark rather than \
                 a position on this transcript's CDS axis"
            } else if pos.utr3 {
                "the 3'UTR zone is numbered from `*1`, so `*0` names no nucleotide"
            } else {
                "the CDS base is unknown"
            };
            return Err(FerroError::ConversionError {
                msg: format!("cannot resolve c.{pos} to a transcript position: {shape}"),
            });
        }

        let cds_start = self
            .transcript
            .cds_start
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no CDS".to_string(),
            })?;

        let cds_end = self
            .transcript
            .cds_end
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no CDS end".to_string(),
            })?;

        // Every branch below adds an unbounded `i64` (nothing bounds a parsed
        // `c.` coordinate — `parse_hgvs` holds no provider) to a transcript
        // coordinate, so every branch can overflow. Refusing is the answer
        // #1658 already took at the three sibling sites in `normalize/merge.rs`:
        // a coordinate that cannot be represented is not one this transcript
        // has, and a wrapped one is worse than an error because release builds
        // then compute against a negative position in silence.
        let overflowed = || FerroError::ConversionError {
            msg: format!(
                "cannot resolve c.{pos} to a transcript position on {}: the coordinate is out \
                 of range",
                self.transcript.id
            ),
        };
        let cds_start = i64::try_from(cds_start).map_err(|_| overflowed())?;
        let cds_end = i64::try_from(cds_end).map_err(|_| overflowed())?;

        let tx_base = if pos.utr3 {
            // 3' UTR: *N = cds_end + N
            cds_end.checked_add(pos.base)
        } else if pos.base < 1 {
            // 5' UTR: -N = cds_start - N (N is positive in HGVS, stored negative),
            // e.g. c.-3 is 3 bases before the CDS start.
            cds_start.checked_add(pos.base)
        } else {
            // Normal CDS position: c.1 = cds_start, c.2 = cds_start + 1, …
            cds_start
                .checked_add(pos.base)
                .and_then(|b| b.checked_sub(1))
        }
        .ok_or_else(overflowed)?;

        // `downstream: false` is a decision, not a default. `CdsPos.utr3`
        // (`c.*N`) is the 3'UTR of a *coding* transcript, which lies **inside**
        // the transcript — `c.*1` is the base after the stop codon, still
        // transcribed — whereas `TxPos.downstream` (`n.*N`) marks a nucleotide
        // past the transcript's last base. The two `*` markers are not the same
        // marker, so carrying `utr3` across into `downstream` would be wrong;
        // the arms above have already resolved `c.*N` to its plain `n.` base.
        // Pinned by `cds_to_tx_keeps_a_coding_three_prime_utr_position_in_transcript`.
        Ok(TxPos {
            base: tx_base,
            offset: pos.offset,
            downstream: false,
        })
    }

    /// Convert transcript position to CDS position
    ///
    /// The inverse of [`Self::cds_to_tx`], and on the same flat sequence axis:
    /// no exon walk, no CIGAR adjustment. The one input on which the round trip
    /// does not close is `c.0` — it converts to `cds_start` and comes back as
    /// `c.1`, because HGVS has no `c.0` for it to come back as.
    ///
    /// # Refuses `n.*N`
    ///
    /// [`TxPos::downstream`] marks the `n.*N` notation, which names a nucleotide
    /// **past the transcript's last base**. There is no CDS position to convert
    /// it to, so this refuses rather than answering:
    ///
    /// - `background/numbering.md:52` gives the whole of `n.` numbering as
    ///   "`n.1`, `n.2`, `n.3`, ..., etc., from the first to the last nucleotide
    ///   of the reference sequence" — no `*` zone exists on this axis;
    /// - `:54` states that "it is **not** allowed to describe variants in
    ///   nucleotides beyond the boundaries of a transcript reference sequence,
    ///   using that transcript reference sequence".
    ///
    /// Before this refusal existed the flag was simply unread: classification
    /// compared `base` against the CDS bounds, so `n.*5` on a transcript with a
    /// 5-base 5'UTR was answered as `c.-1` — a different nucleotide, on the
    /// wrong side of the CDS, with no diagnostic. That was already known and
    /// worked around one layer up, where `project_*_direct` rejects `*N` inputs
    /// "because `tx_to_cds` ... never reads this flag"; the check now lives with
    /// the contract it belongs to, and that caller's earlier, better-worded
    /// refusal still runs first.
    pub fn tx_to_cds(&self, pos: &TxPos) -> Result<CdsPos, FerroError> {
        if pos.downstream {
            return Err(FerroError::ConversionError {
                msg: format!(
                    "n.*{} lies beyond the transcript's last base and has no CDS position: \
                     the n. axis numbers only n.1..n.<length> (background/numbering.md:52) \
                     and a variant outside a transcript's boundaries may not be described \
                     against that transcript (background/numbering.md:54)",
                    pos.base
                ),
            });
        }

        let cds_start = self
            .transcript
            .cds_start
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no CDS".to_string(),
            })? as i64;
        let cds_end = self
            .transcript
            .cds_end
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no CDS end".to_string(),
            })? as i64;

        let base = pos.base;

        if base < cds_start {
            // 5' UTR: `c.-N`, N bases 5' of the CDS start.
            Ok(CdsPos {
                base: base - cds_start, // negative for 5' UTR
                offset: pos.offset,
                utr3: false,
                special: None,
            })
        } else if base > cds_end {
            // 3' UTR: `c.*N`, N bases 3' of the last CDS base.
            Ok(CdsPos {
                base: base - cds_end,
                offset: pos.offset,
                utr3: true,
                special: None,
            })
        } else {
            // Within the CDS.
            Ok(CdsPos {
                base: base - cds_start + 1,
                offset: pos.offset,
                utr3: false,
                special: None,
            })
        }
    }

    /// Convert CDS position to protein position
    ///
    /// Protein positions are (CDS position + 2) / 3 (rounded up)
    pub fn cds_to_protein(&self, pos: &CdsPos) -> Result<ProtPos, FerroError> {
        if pos.base < 1 || pos.utr3 {
            return Err(FerroError::ConversionError {
                msg: "Cannot convert UTR position to protein".to_string(),
            });
        }

        if pos.is_intronic() {
            return Err(FerroError::ConversionError {
                msg: "Cannot convert intronic position to protein".to_string(),
            });
        }

        // Protein position is 1-indexed, and each codon is 3 bases
        let aa_number = ((pos.base - 1) / 3 + 1) as u64;

        // Get the amino acid at this position if we have the sequence
        // For now, use Xaa (unknown)
        let aa = self.get_amino_acid_at(aa_number).unwrap_or(AminoAcid::Xaa);

        Ok(ProtPos::new(aa, aa_number))
    }

    /// Get the amino acid at a protein position (if sequence is available)
    ///
    /// Returns `None` whenever no codon can be read — including when the
    /// transcript's own CDS metadata does not name one. `Transcript` stores
    /// `cds_start` exactly as its provider supplied it and validates nothing
    /// (`validate.rs`'s CDS check is out-of-band and report-only), so `0` is
    /// representable, and `0` is a 1-based coordinate naming no transcript
    /// base. `Transcript::utr5_length` defends the same field the same way.
    fn get_amino_acid_at(&self, position: u64) -> Option<AminoAcid> {
        let cds_start = self.transcript.cds_start?;
        let cds_end = self.transcript.cds_end?;

        // Get codon start position (0-based in sequence).
        //
        // Both subtrahends are 1-based, so `0` is not an off-by-one to absorb
        // on either — it names no base and no residue, and there is nothing to
        // read. Unchecked, this panicked in *both* profiles and the release
        // failure moved: debug aborted here on `attempt to subtract with
        // overflow`, while release wrapped to `usize::MAX` and aborted on the
        // slice below, reporting a string index instead of the bad metadata.
        let cds_start_0based = cds_start.checked_sub(1)?;
        let codon_offset = position.checked_sub(1)?.checked_mul(3)?;
        let codon_start = usize::try_from(cds_start_0based.checked_add(codon_offset)?).ok()?;
        let codon_end = codon_start.checked_add(3)?;

        // Need cached transcript bases to translate; coordinate-only
        // transcripts return `None`.
        let seq = self.transcript.sequence.as_deref()?;
        if codon_end > cds_end as usize || codon_end > seq.len() {
            return None;
        }

        let codon = &seq[codon_start..codon_end];
        Self::translate_codon(codon)
    }

    /// Translate a codon to an amino acid.
    ///
    /// Total: every one of the 22 match arms (21 named plus the `_`
    /// fall-through) returns `Some`, so this never returns `None` despite
    /// the `Option` signature. The `_` arm is what makes it total: it
    /// coerces *every* unrecognised codon — including `N`-containing ones
    /// from real reference sequences — to `Some(Xaa)`, indistinguishable
    /// from a genuine unknown residue.
    fn translate_codon(codon: &str) -> Option<AminoAcid> {
        match codon.to_uppercase().as_str() {
            "TTT" | "TTC" => Some(AminoAcid::Phe),
            "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => Some(AminoAcid::Leu),
            "ATT" | "ATC" | "ATA" => Some(AminoAcid::Ile),
            "ATG" => Some(AminoAcid::Met),
            "GTT" | "GTC" | "GTA" | "GTG" => Some(AminoAcid::Val),
            "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => Some(AminoAcid::Ser),
            "CCT" | "CCC" | "CCA" | "CCG" => Some(AminoAcid::Pro),
            "ACT" | "ACC" | "ACA" | "ACG" => Some(AminoAcid::Thr),
            "GCT" | "GCC" | "GCA" | "GCG" => Some(AminoAcid::Ala),
            "TAT" | "TAC" => Some(AminoAcid::Tyr),
            "TAA" | "TAG" | "TGA" => Some(AminoAcid::Ter),
            "CAT" | "CAC" => Some(AminoAcid::His),
            "CAA" | "CAG" => Some(AminoAcid::Gln),
            "AAT" | "AAC" => Some(AminoAcid::Asn),
            "AAA" | "AAG" => Some(AminoAcid::Lys),
            "GAT" | "GAC" => Some(AminoAcid::Asp),
            "GAA" | "GAG" => Some(AminoAcid::Glu),
            "TGT" | "TGC" => Some(AminoAcid::Cys),
            "TGG" => Some(AminoAcid::Trp),
            "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => Some(AminoAcid::Arg),
            "GGT" | "GGC" | "GGA" | "GGG" => Some(AminoAcid::Gly),
            _ => Some(AminoAcid::Xaa), // Unknown
        }
    }

    /// Convert genomic position to transcript position
    ///
    /// Returns None if the genomic position is in an intron.
    /// The transcript must have genomic coordinates set on its exons.
    ///
    /// **This is a genome-frame conversion and it does consult the exon
    /// table** — see the module's "two frames" note, and contrast
    /// [`Self::cds_to_tx`], which must not. The two answer different questions
    /// and therefore *must* disagree on a transcript whose alignment carries a
    /// transcript-coordinate gap; `tests/it/axis_frame_disagreement.rs` pins
    /// that, so collapsing them onto one implementation fails loudly rather
    /// than reintroducing #1619.
    pub fn genomic_to_tx(&self, genomic_pos: u64) -> Result<Option<TxPos>, FerroError> {
        // Check if transcript has genomic coordinates
        if !self.transcript.has_genomic_coords() {
            return Err(FerroError::ConversionError {
                msg: "Transcript does not have genomic coordinates".to_string(),
            });
        }

        // Find which exon contains this genomic position
        for exon in &self.transcript.exons {
            let (g_start, g_end) = match (exon.genomic_start, exon.genomic_end) {
                (Some(s), Some(e)) => (s, e),
                _ => continue,
            };

            // Check if position is within this exon's genomic coordinates
            if genomic_pos >= g_start && genomic_pos <= g_end {
                // Calculate offset within exon based on strand
                let tx_pos = match self.transcript.strand {
                    Strand::Plus => {
                        // Plus strand: genomic position increases with transcript position
                        let offset_in_exon = genomic_pos - g_start;
                        exon.start + offset_in_exon
                    }
                    Strand::Minus => {
                        // Minus strand: genomic position decreases with transcript position
                        let offset_in_exon = g_end - genomic_pos;
                        exon.start + offset_in_exon
                    }
                    Strand::Unknown => return Ok(None),
                };

                return Ok(Some(TxPos::new(tx_pos as i64)));
            }
        }

        // Position is in an intron (between exons)
        Ok(None)
    }

    /// Convert genomic position to transcript position with intronic offset support
    ///
    /// Unlike `genomic_to_tx`, this method returns a TxPos with an offset for
    /// intronic positions instead of returning None.
    pub fn genomic_to_tx_with_intron(&self, genomic_pos: u64) -> Result<TxPos, FerroError> {
        // First check if it's exonic
        if let Some(tx_pos) = self.genomic_to_tx(genomic_pos)? {
            return Ok(tx_pos);
        }

        // It's intronic - find the intron and calculate offset
        if let Some((_intron, intron_pos)) = self.transcript.find_intron_at_genomic(genomic_pos) {
            return Ok(TxPos::with_offset(
                intron_pos.tx_boundary_pos as i64,
                intron_pos.offset,
            ));
        }

        Err(FerroError::ConversionError {
            msg: format!(
                "Genomic position {} is outside transcript bounds",
                genomic_pos
            ),
        })
    }

    /// Get intronic position information for a genomic position
    ///
    /// Returns None if the position is exonic or outside the transcript.
    pub fn get_intron_position(&self, genomic_pos: u64) -> Option<IntronPosition> {
        self.transcript
            .find_intron_at_genomic(genomic_pos)
            .map(|(_, pos)| pos)
    }

    /// Convert transcript position to genomic position
    ///
    /// Returns None if the transcript position is not covered by exons with genomic coords.
    ///
    /// **Genome frame**, like [`Self::genomic_to_tx`] and unlike
    /// [`Self::tx_to_cds`]: the alignment answer of `None` means "this
    /// transcript base aligns to nothing", which is a fact about the alignment
    /// and not about the transcript's own numbering — that base still has a
    /// `c.`/`n.` number.
    ///
    /// `None` is not exclusively that answer, though, so do not read it as
    /// "unaligned" without ruling the other two out: it is also returned for a
    /// non-positive `tx_pos.base`, and for a `Strand::Unknown` transcript,
    /// whose bases are all reported unaligned however complete its exon table
    /// is. A caller using `tx_to_genomic(..).is_some()` as its definition of
    /// "aligned" — `tests/it/axis_frame_disagreement.rs`'s genome-frame walk
    /// does — is relying on the transcript having a known strand.
    pub fn tx_to_genomic(&self, tx_pos: &TxPos) -> Result<Option<u64>, FerroError> {
        // Check if transcript has genomic coordinates
        if !self.transcript.has_genomic_coords() {
            return Err(FerroError::ConversionError {
                msg: "Transcript does not have genomic coordinates".to_string(),
            });
        }

        // Find which exon contains this transcript position
        // Only positive transcript positions can be in exons
        if tx_pos.base < 1 {
            return Ok(None);
        }
        let tx_base = tx_pos.base as u64;

        for exon in &self.transcript.exons {
            if tx_base >= exon.start && tx_base <= exon.end {
                let (g_start, g_end) = match (exon.genomic_start, exon.genomic_end) {
                    (Some(s), Some(e)) => (s, e),
                    _ => continue,
                };

                // Calculate offset within exon
                let offset_in_exon = tx_base - exon.start;

                let genomic_pos = match self.transcript.strand {
                    Strand::Plus => {
                        // Plus strand: transcript position increases with genomic position
                        g_start + offset_in_exon
                    }
                    Strand::Minus => {
                        // Minus strand: transcript position increases as genomic position decreases
                        g_end - offset_in_exon
                    }
                    Strand::Unknown => return Ok(None),
                };

                return Ok(Some(genomic_pos));
            }
        }

        // Position not found in exons
        Ok(None)
    }

    /// Convert genomic position to CDS position
    ///
    /// Returns None if the position is intronic.
    pub fn genomic_to_cds(&self, genomic_pos: u64) -> Result<Option<CdsPos>, FerroError> {
        // First convert to transcript position
        let tx_pos = self.genomic_to_tx(genomic_pos)?;

        match tx_pos {
            Some(pos) => Ok(Some(self.tx_to_cds(&pos)?)),
            None => Ok(None),
        }
    }

    /// Convert genomic position to CDS position with intronic offset support
    ///
    /// Unlike `genomic_to_cds`, this method returns a CdsPos with an offset for
    /// intronic positions instead of returning None.
    pub fn genomic_to_cds_with_intron(&self, genomic_pos: u64) -> Result<CdsPos, FerroError> {
        // Get tx position with intron support
        let tx_pos = self.genomic_to_tx_with_intron(genomic_pos)?;

        // Convert to CDS, preserving the offset
        self.tx_to_cds(&tx_pos)
    }

    /// Convert CDS position to genomic position
    pub fn cds_to_genomic(&self, cds_pos: &CdsPos) -> Result<Option<u64>, FerroError> {
        // First convert to transcript position
        let tx_pos = self.cds_to_tx(cds_pos)?;
        self.tx_to_genomic(&tx_pos)
    }

    /// Get the chromosome name for this transcript
    pub fn chromosome(&self) -> Option<&str> {
        self.transcript.chromosome.as_deref()
    }

    /// Convert CDS position with intronic offset to genomic position
    ///
    /// For intronic variants like c.100+5 or c.200-10, this calculates the
    /// genomic position by:
    /// 1. Converting the CDS base position to transcript position
    /// 2. Using the intron mapping to find the genomic position
    ///
    /// # Returns
    /// The genomic position, or an error if the transcript lacks genomic coordinates
    pub fn cds_to_genomic_with_intron(&self, cds_pos: &CdsPos) -> Result<u64, FerroError> {
        // First convert to transcript position
        let tx_pos = self.cds_to_tx(cds_pos)?;

        // If not intronic, use the standard method
        if !cds_pos.is_intronic() {
            return self
                .tx_to_genomic(&tx_pos)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("Position {} not found in exons", tx_pos.base),
                });
        }

        // For intronic positions, use the intronic_to_genomic method
        let offset = cds_pos.offset.ok_or_else(|| FerroError::ConversionError {
            msg: "Expected intronic offset".to_string(),
        })?;

        self.transcript
            .intronic_to_genomic(tx_pos.base as u64, offset)
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Could not convert intronic position {}+{} to genomic",
                    tx_pos.base, offset
                ),
            })
    }

    /// Fold a UTR offset that has *no enclosing intron* into a plain
    /// (non-intronic) CDS coordinate.
    ///
    /// A `+offset`/`-offset` is only meaningful when an intron spans the
    /// boundary. In the UTRs that intron can be absent for two reasons:
    /// - the boundary is the transcript terminus, so the offset points into the
    ///   contiguous genome past the final/first exon (`c.*824+10` 3' of the last
    ///   exon → `c.*834`); or
    /// - the transcript model simply has no intron there (e.g. a single-exon
    ///   supplemental transcript), so the offset is a spurious decoration.
    ///
    /// In both cases the offset is not intronic and folds linearly into the
    /// base coordinate, matching the canonical interpretation and letting the
    /// standard (non-intronic) machinery map and normalize the position instead
    /// of the intron-window shuffle that has no intron to anchor on (#760).
    ///
    /// # Returns
    /// The folded position, or `None` when it should keep its existing
    /// handling: no offset, a position in the CDS proper (where a missing
    /// intron means a genuinely invalid description, not an extended
    /// coordinate), or a genuine enclosing intron.
    pub fn fold_non_intronic_utr_offset(&self, pos: &CdsPos) -> Option<CdsPos> {
        let offset = pos.offset?;
        if offset == 0 {
            return None;
        }

        // Only UTR offsets fold. A mid-CDS offset with no intron is invalid, not
        // an extended coordinate, so it keeps erroring rather than folding.
        let in_utr = pos.utr3 || pos.base < 0;
        if !in_utr {
            return None;
        }

        let tx = self.cds_to_tx(pos).ok()?;

        // A negative tx base sits upstream of transcript base 1 (a deep 5'UTR
        // offset, e.g. `c.-N` with N past the transcript start), so no interior
        // intron can enclose it. Skip the intron probe — which is `u64`-indexed
        // and cannot represent a pre-transcript base — and fold linearly below.
        // A genuine enclosing intron only applies for a non-negative tx base.
        if tx.base >= 0 {
            let tx_base = tx.base as u64;
            // A genuine intron encloses the offset → leave it to the intronic path.
            if self
                .transcript
                .find_intron_at_tx_boundary(tx_base, offset)
                .is_some()
            {
                return None;
            }
        }

        // Transcript coordinates run 5'->3' regardless of strand, so a signed
        // offset shifts the boundary directly; the (contiguous) sequence past
        // the boundary maps back through the standard CDS conversion, which
        // accepts coordinates beyond the transcript termini.
        let folded_tx = tx.base.checked_add(offset)?;
        self.tx_to_cds(&TxPos::new(folded_tx)).ok()
    }

    /// Convert transcript position with intronic offset to genomic position
    pub fn tx_to_genomic_with_intron(&self, tx_pos: &TxPos) -> Result<u64, FerroError> {
        // If not intronic, use the standard method
        if tx_pos.offset.is_none() {
            return self
                .tx_to_genomic(tx_pos)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("Position {} not found in exons", tx_pos.base),
                });
        }

        let offset = tx_pos.offset.unwrap();

        self.transcript
            .intronic_to_genomic(tx_pos.base as u64, offset)
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Could not convert intronic position {}+{} to genomic",
                    tx_pos.base, offset
                ),
            })
    }

    /// Convert genomic position to CDS position with intronic offset support
    pub fn genomic_to_cds_intronic(&self, genomic_pos: u64) -> Result<CdsPos, FerroError> {
        // First check if it's exonic
        if let Some(cds_pos) = self.genomic_to_cds(genomic_pos)? {
            return Ok(cds_pos);
        }

        // It's intronic - get the transcript boundary and offset
        let (tx_boundary, offset) = self
            .transcript
            .genomic_to_intronic(genomic_pos)
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!("Genomic position {} is not within transcript", genomic_pos),
            })?;

        // Convert transcript boundary to CDS
        let cds_boundary = self.tx_to_cds(&TxPos::new(tx_boundary as i64))?;

        Ok(CdsPos {
            base: cds_boundary.base,
            offset: Some(offset),
            utr3: cds_boundary.utr3,
            special: None,
        })
    }
}

/// Result of genomic coordinate lookup
#[derive(Debug, Clone, PartialEq)]
pub struct GenomicLocation {
    /// Chromosome name
    pub chromosome: String,
    /// Genomic position (1-based, HGVS g. format)
    pub position: u64,
    /// Strand orientation
    pub strand: Strand,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::CigarOp;
    use crate::hgvs::location::SpecialPosition;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            cds_start_incomplete: false,
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            // 5'UTR (5bp) + CDS (30bp) + 3'UTR (3bp) = 38bp.
            // NOTE: cds_start = 6 is the `T` of the visible `ATG` (whose `A` is at
            // 0-based 4), so the first codon is seq[5..8] == "TGC" and p.1 is Cys,
            // not Met. Use `make_codon_sweep_transcript` for amino-acid assertions.
            sequence: Some("AAAAATGCCCAAAGGGTTTAGGCCCAAAGGGTTATAAA".to_string()),
            cds_start: Some(6),
            cds_end: Some(35),
            exons: vec![Exon::new(1, 1, 38)],
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

    /// One codon per `translate_codon` match arm, 21 codons, all 21 amino acids
    /// distinct. Three properties are load-bearing — do not "tidy" them:
    ///
    /// 1. All 21 amino acids differ, so a mutated `codon_start` reads a codon
    ///    that translates differently rather than coinciding.
    /// 2. The 5'UTR's first codon is `CCC` (Pro), never `ATG`. In
    ///    `get_amino_acid_at`'s `codon_start` arithmetic
    ///    (`cds_start as usize - 1 + (position as usize - 1) * 3`), a
    ///    `+ -> *` mutant degenerates `codon_start` to 0 at p.1 and reads
    ///    seq[0..3]; a Met-translating 5'UTR would let it survive.
    /// 3. The 3'UTR translates (`GGG` -> Gly). In `get_amino_acid_at`'s
    ///    past-CDS guard (`codon_end > cds_end as usize || codon_end >
    ///    seq.len()`), a `|| -> &&` mutant stops it firing at p.22 and reads
    ///    seq[68..71]; Gly there is what distinguishes it from the `Xaa` the
    ///    original returns.
    ///
    /// Layout: 5'UTR 5bp (cds_start = 6), CDS 63bp (cds_end = 68), 3'UTR 6bp,
    /// total 74bp.
    fn make_codon_sweep_transcript() -> Transcript {
        Transcript {
            id: "NM_CODONS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            sequence: Some(
                "CCCCC".to_string()
                    + "ATGTTTCTGATTGTGTCACCCACAGCCTAT"
                    + "CATCAAAATAAAGATGAATGTTGGCGCGGA"
                    + "TAA"
                    + "GGGTTT",
            ),
            cds_start: Some(6),
            cds_end: Some(68),
            exons: vec![Exon::new(1, 1, 74)],
            ..Default::default()
        }
    }

    /// A deliberately malformed transcript whose `cds_end` exceeds its sequence
    /// length. This is the ONLY way to reach `get_amino_acid_at`'s second guard
    /// clause: under `||`, `codon_end > seq.len()` is unreachable whenever
    /// `codon_end > cds_end` has already fired.
    ///
    /// Sequence is 12bp, cds_start = 1, cds_end = 100.
    ///
    /// Safe to reuse with `validate_cds_pos`: its 3'UTR guard computes
    /// `transcript.sequence_length().checked_sub(cds_end)`, so `12 - 100` here
    /// declines with a malformed-transcript `ConversionError` rather than
    /// underflowing and panicking in debug (#1909).
    fn make_cds_end_past_sequence_transcript() -> Transcript {
        Transcript {
            id: "NM_SHORTSEQ.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            sequence: Some("ATGTTTCTGATT".to_string()),
            cds_start: Some(1),
            cds_end: Some(100),
            exons: vec![Exon::new(1, 1, 12)],
            ..Default::default()
        }
    }

    /// Pins the amino acid at every codon of a 21-arm sweep.
    ///
    /// `test_cds_to_protein` asserts only `.number`, leaving the amino acid
    /// itself unchecked: `get_amino_acid_at` returning `None` is invisible
    /// there because `cds_to_protein` swallows it with `.unwrap_or(Xaa)`,
    /// and every `translate_codon` arm falls through to `_ => Some(Xaa)`
    /// when deleted, so a mutated arm and a correct one are indistinguishable
    /// by number alone. This test checks `.aa` directly so both are caught.
    #[test]
    fn cds_to_protein_pins_every_codon_arm() {
        let tx = make_codon_sweep_transcript();
        let mapper = CoordinateMapper::new(&tx);

        let expected = [
            AminoAcid::Met,
            AminoAcid::Phe,
            AminoAcid::Leu,
            AminoAcid::Ile,
            AminoAcid::Val,
            AminoAcid::Ser,
            AminoAcid::Pro,
            AminoAcid::Thr,
            AminoAcid::Ala,
            AminoAcid::Tyr,
            AminoAcid::His,
            AminoAcid::Gln,
            AminoAcid::Asn,
            AminoAcid::Lys,
            AminoAcid::Asp,
            AminoAcid::Glu,
            AminoAcid::Cys,
            AminoAcid::Trp,
            AminoAcid::Arg,
            AminoAcid::Gly,
            AminoAcid::Ter,
        ];

        for (index, aa) in expected.iter().enumerate() {
            let protein_number = (index + 1) as u64;
            let cds_base = (index as i64) * 3 + 1;
            let got = mapper.cds_to_protein(&CdsPos::new(cds_base)).unwrap();
            assert_eq!(got.number, protein_number, "p.{} number", protein_number);
            assert_eq!(
                got.aa, *aa,
                "p.{} amino acid (c.{})",
                protein_number, cds_base
            );
        }
    }

    /// p.21 is the last codon (`codon_end == cds_end`, 68) and must translate;
    /// p.22 runs past the CDS and must yield `Xaa`.
    ///
    /// The p.22 case is what kills a `|| -> &&` mutant of `get_amino_acid_at`'s
    /// guard (the `codon_end > cds_end as usize || codon_end > seq.len()`
    /// early return in `get_amino_acid_at`): with `&&` the guard does
    /// not fire, the function reads the 3'UTR's `GGG`, and `cds_to_protein`
    /// returns Gly instead of Xaa.
    ///
    /// The two indirect assertions above observe the guard only through
    /// `cds_to_protein`'s `.unwrap_or(AminoAcid::Xaa)`, which — exactly as
    /// `get_amino_acid_at_guards_a_cds_end_past_the_sequence`'s doc comment
    /// diagnoses for its own case — collapses five distinct causes of
    /// `None` into one `Xaa` value: a missing `cds_start`, a missing
    /// `cds_end`, a missing cached `sequence`, `codon_end > cds_end` (the
    /// primary guard clause), and `codon_end > seq.len()` (the secondary
    /// one). `translate_codon` is not a sixth cause — it is total (every
    /// match arm, including `_`, yields `Some`, coercing any unrecognised
    /// codon to `Xaa` rather than `None`), so it can never contribute a
    /// `None` here. That test unmasks the *secondary*,
    /// malformed-transcript-only guard clause (`codon_end > seq.len()`) by
    /// calling `get_amino_acid_at` directly; this test exercises the
    /// *primary* clause (`codon_end > cds_end as usize`), the one that fires
    /// on a well-formed transcript, so it gets the same direct treatment
    /// here.
    #[test]
    fn cds_to_protein_stops_at_the_cds_end() {
        let tx = make_codon_sweep_transcript();
        let mapper = CoordinateMapper::new(&tx);

        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(61)).unwrap().aa,
            AminoAcid::Ter,
            "c.61 is the last codon and must still translate"
        );
        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(64)).unwrap().aa,
            AminoAcid::Xaa,
            "c.64 is past cds_end and must not read into the 3'UTR"
        );

        // Direct, unmasked: `get_amino_acid_at` itself must draw the p.21/p.22
        // line at `codon_end > cds_end`, not just leave `cds_to_protein` to
        // paper over it with `Xaa`.
        assert_eq!(
            mapper.get_amino_acid_at(21),
            Some(AminoAcid::Ter),
            "get_amino_acid_at(21) must translate unmasked, codon_end 68 == cds_end 68"
        );
        assert_eq!(
            mapper.get_amino_acid_at(22),
            None,
            "get_amino_acid_at(22) must return None unmasked, codon_end 71 > cds_end 68"
        );
    }

    /// Exercises `get_amino_acid_at`'s SECOND guard clause
    /// (`codon_end > seq.len()`), which is unreachable on a well-formed
    /// transcript because the first clause fires first.
    ///
    /// cds_end = 100 but the sequence is 12bp, so:
    ///   p.1 codon_end = 3  (< 12)  -> translates
    ///   p.4 codon_end = 12 (== 12) -> translates, the boundary
    ///   p.5 codon_end = 15 (> 12)  -> Xaa
    ///
    /// The assertions below observe the guard only through `cds_to_protein`,
    /// which coerces every `None` from `get_amino_acid_at` — regardless of
    /// which of at least five distinct causes produced it — into `Xaa` via
    /// `cds_to_protein`'s `.unwrap_or(AminoAcid::Xaa)` fallback. That masks
    /// the guard: the assertion cannot fail for the *right* reason. The
    /// direct call below observes `get_amino_acid_at` unmasked, alongside
    /// the existing indirect ones.
    #[test]
    fn get_amino_acid_at_guards_a_cds_end_past_the_sequence() {
        let tx = make_cds_end_past_sequence_transcript();
        let mapper = CoordinateMapper::new(&tx);

        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(1)).unwrap().aa,
            AminoAcid::Met,
            "codon_end 3 is well inside a 12bp sequence"
        );
        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(10)).unwrap().aa,
            AminoAcid::Ile,
            "codon_end 12 == seq.len() is the last readable codon"
        );
        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(13)).unwrap().aa,
            AminoAcid::Xaa,
            "codon_end 15 > seq.len() must not read out of bounds"
        );

        // Direct, unmasked: `get_amino_acid_at` itself must return `None`
        // past the end (p.5, codon_end 15 > seq.len() 12), not just leave
        // `cds_to_protein` to paper over it with `Xaa`.
        assert_eq!(
            mapper.get_amino_acid_at(5),
            None,
            "get_amino_acid_at(5) must return None unmasked, codon_end 15 > seq.len() 12"
        );
    }

    #[test]
    fn test_cds_to_tx_normal() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // Pins `make_test_transcript`'s length (5'UTR 5bp + CDS 30bp + 3'UTR
        // 3bp = 38bp) against the hardcoded `38` literal below — that catches
        // the sequence literal changing size, not the comment drifting from
        // it (both copies once said `= 40bp` while the sequence stayed 38bp,
        // and this assertion would not have caught that). The fixture is
        // duplicated with the same data but different code in
        // `tests/it/convert_tests.rs` (a `Transcript { .. }` literal here vs.
        // `Transcript::new(..)` there). It is also not sequence-specific:
        // `sequence_length()` falls back to summing exon lengths when
        // `sequence` is `None`, and this fixture's `Exon::new(1, 1, 38)`
        // sums to 38 too, so dropping the sequence would still pass.
        assert_eq!(tx.sequence_length(), 38, "make_test_transcript is 38bp");

        // c.1 should be tx position 6
        let result = mapper.cds_to_tx(&CdsPos::new(1)).unwrap();
        assert_eq!(result.base, 6);
    }

    #[test]
    fn test_cds_to_tx_5utr() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // c.-3 should be tx position 3
        let result = mapper.cds_to_tx(&CdsPos::new(-3)).unwrap();
        assert_eq!(result.base, 3);
    }

    #[test]
    fn test_cds_to_tx_3utr() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // c.*2 should be tx position 37 (cds_end=35 + 2)
        let result = mapper.cds_to_tx(&CdsPos::utr3(2)).unwrap();
        assert_eq!(result.base, 37);
    }

    #[test]
    fn test_tx_to_cds_normal() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // tx position 10 should be c.5 (10 - 6 + 1)
        let result = mapper.tx_to_cds(&TxPos::new(10)).unwrap();
        assert_eq!(result.base, 5);
        assert!(!result.utr3);
    }

    #[test]
    fn test_tx_to_cds_5utr() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // tx position 3 should be c.-3
        let result = mapper.tx_to_cds(&TxPos::new(3)).unwrap();
        assert_eq!(result.base, -3);
    }

    #[test]
    fn test_tx_to_cds_3utr() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // tx position 37 should be c.*2
        let result = mapper.tx_to_cds(&TxPos::new(37)).unwrap();
        assert_eq!(result.base, 2);
        assert!(result.utr3);
    }

    /// `c.*3+5` on a single-exon transcript whose 3'UTR ends at `*3` has no
    /// enclosing intron — it is the extended 3'UTR coordinate `c.*8`. Folding
    /// lets the standard (non-intronic) machinery map it. Regression for #760.
    #[test]
    fn fold_three_prime_utr_offset_at_terminus() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // cds_end=35, single exon ends at tx 38 → *3 is the final transcript base.
        let pos = CdsPos {
            base: 3,
            offset: Some(5),
            utr3: true,
            special: None,
        };
        let folded = mapper
            .fold_non_intronic_utr_offset(&pos)
            .expect("should fold");
        assert_eq!(folded.base, 8);
        assert!(folded.utr3);
        assert_eq!(folded.offset, None);
    }

    /// A 3'UTR offset *interior* to a single-exon transcript (no intron there at
    /// all, the degenerate supplemental-transcript shape behind the #760 corpus
    /// row) still folds: `c.*2+1` → `c.*3`.
    #[test]
    fn fold_three_prime_utr_offset_interior() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // *2 maps to tx 37, interior to the single exon (which ends at tx 38).
        let pos = CdsPos {
            base: 2,
            offset: Some(1),
            utr3: true,
            special: None,
        };
        let folded = mapper
            .fold_non_intronic_utr_offset(&pos)
            .expect("should fold");
        assert_eq!(folded.base, 3);
        assert!(folded.utr3);
        assert_eq!(folded.offset, None);
    }

    /// `c.-5-3` anchored at the first transcript base extends 5' past the
    /// terminus to `c.-8`.
    #[test]
    fn fold_five_prime_utr_offset_at_terminus() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // cds_start=6 → c.-5 maps to tx 1, the first transcript base.
        let pos = CdsPos {
            base: -5,
            offset: Some(-3),
            utr3: false,
            special: None,
        };
        let folded = mapper
            .fold_non_intronic_utr_offset(&pos)
            .expect("should fold");
        assert_eq!(folded.base, -8);
        assert!(!folded.utr3);
        assert_eq!(folded.offset, None);
    }

    /// A deep 5'UTR offset whose anchor base already maps upstream of
    /// transcript base 1 (a negative tx base) must still fold linearly rather
    /// than short-circuit into the intronic path: `c.-7-2` (tx -1) → `c.-9`.
    /// Regression for the negative-tx-base 5'UTR edge of #760.
    #[test]
    fn fold_five_prime_utr_offset_negative_tx_base() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // cds_start=6 → c.-7 maps to tx -1, upstream of transcript base 1.
        let pos = CdsPos {
            base: -7,
            offset: Some(-2),
            utr3: false,
            special: None,
        };
        let folded = mapper
            .fold_non_intronic_utr_offset(&pos)
            .expect("should fold");
        assert_eq!(folded.base, -9);
        assert!(!folded.utr3);
        assert_eq!(folded.offset, None);
    }

    /// A mid-CDS offset with no intron is a genuinely invalid description, not
    /// an extended coordinate — it must not be folded.
    #[test]
    fn fold_offset_declines_cds_position() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // c.10+5 maps to tx 15, inside the CDS — not a UTR offset.
        let pos = CdsPos {
            base: 10,
            offset: Some(5),
            utr3: false,
            special: None,
        };
        assert!(mapper.fold_non_intronic_utr_offset(&pos).is_none());
    }

    /// A non-intronic position (no offset) is never folded.
    #[test]
    fn fold_offset_declines_non_intronic_position() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        assert!(mapper
            .fold_non_intronic_utr_offset(&CdsPos::utr3(2))
            .is_none());
    }

    #[test]
    fn test_cds_to_protein() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // c.1, c.2, c.3 should all map to p.1
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(1)).unwrap().number, 1);
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(2)).unwrap().number, 1);
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(3)).unwrap().number, 1);

        // c.4, c.5, c.6 should all map to p.2
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(4)).unwrap().number, 2);
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(5)).unwrap().number, 2);
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(6)).unwrap().number, 2);
    }

    #[test]
    fn test_translate_codon() {
        assert_eq!(
            CoordinateMapper::translate_codon("ATG"),
            Some(AminoAcid::Met)
        );
        assert_eq!(
            CoordinateMapper::translate_codon("TAA"),
            Some(AminoAcid::Ter)
        );
        assert_eq!(
            CoordinateMapper::translate_codon("GGG"),
            Some(AminoAcid::Gly)
        );
        // An unrecognised codon falls through to the `_ => Some(Xaa)` arm.
        assert_eq!(
            CoordinateMapper::translate_codon("NNN"),
            Some(AminoAcid::Xaa)
        );
    }

    fn make_genomic_transcript_plus() -> Transcript {
        // Plus strand transcript with 3 exons
        // Exon 1: tx 1-10, genomic 1000-1009
        // Exon 2: tx 11-20, genomic 2000-2009
        // Exon 3: tx 21-30, genomic 3000-3009
        Transcript {
            cds_start_incomplete: false,
            id: "NM_GENOMIC.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("A".repeat(30)),
            cds_start: Some(5),
            cds_end: Some(25),
            exons: vec![
                Exon::with_genomic(1, 1, 10, 1000, 1009),
                Exon::with_genomic(2, 11, 20, 2000, 2009),
                Exon::with_genomic(3, 21, 30, 3000, 3009),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(3009),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    fn make_genomic_transcript_minus() -> Transcript {
        // Minus strand transcript with 3 exons
        // Exon 1: tx 1-10, genomic 3009-3000 (reversed)
        // Exon 2: tx 11-20, genomic 2009-2000 (reversed)
        // Exon 3: tx 21-30, genomic 1009-1000 (reversed)
        Transcript {
            cds_start_incomplete: false,
            id: "NM_GENOMIC_MINUS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Minus,
            sequence: Some("A".repeat(30)),
            cds_start: Some(5),
            cds_end: Some(25),
            exons: vec![
                Exon::with_genomic(1, 1, 10, 3000, 3009),
                Exon::with_genomic(2, 11, 20, 2000, 2009),
                Exon::with_genomic(3, 21, 30, 1000, 1009),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(3009),
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
    fn test_genomic_to_tx_plus_strand() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // Exon 1: genomic 1000 -> tx 1
        assert_eq!(mapper.genomic_to_tx(1000).unwrap(), Some(TxPos::new(1)));
        // Exon 1: genomic 1005 -> tx 6
        assert_eq!(mapper.genomic_to_tx(1005).unwrap(), Some(TxPos::new(6)));
        // Exon 2: genomic 2000 -> tx 11
        assert_eq!(mapper.genomic_to_tx(2000).unwrap(), Some(TxPos::new(11)));
        // Exon 3: genomic 3009 -> tx 30
        assert_eq!(mapper.genomic_to_tx(3009).unwrap(), Some(TxPos::new(30)));
    }

    #[test]
    fn test_genomic_to_tx_minus_strand() {
        let tx = make_genomic_transcript_minus();
        let mapper = CoordinateMapper::new(&tx);

        // Exon 1: genomic 3009 -> tx 1 (minus strand, higher genomic = earlier tx)
        assert_eq!(mapper.genomic_to_tx(3009).unwrap(), Some(TxPos::new(1)));
        // Exon 1: genomic 3000 -> tx 10
        assert_eq!(mapper.genomic_to_tx(3000).unwrap(), Some(TxPos::new(10)));
        // Exon 2: genomic 2009 -> tx 11
        assert_eq!(mapper.genomic_to_tx(2009).unwrap(), Some(TxPos::new(11)));
    }

    #[test]
    fn test_genomic_to_tx_intronic() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // Position 1500 is between exon 1 (1009) and exon 2 (2000) - intronic
        assert_eq!(mapper.genomic_to_tx(1500).unwrap(), None);
    }

    #[test]
    fn test_tx_to_genomic_plus_strand() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // tx 1 -> genomic 1000
        assert_eq!(mapper.tx_to_genomic(&TxPos::new(1)).unwrap(), Some(1000));
        // tx 15 -> genomic 2004
        assert_eq!(mapper.tx_to_genomic(&TxPos::new(15)).unwrap(), Some(2004));
    }

    #[test]
    fn test_tx_to_genomic_minus_strand() {
        let tx = make_genomic_transcript_minus();
        let mapper = CoordinateMapper::new(&tx);

        // tx 1 -> genomic 3009 (minus strand)
        assert_eq!(mapper.tx_to_genomic(&TxPos::new(1)).unwrap(), Some(3009));
        // tx 10 -> genomic 3000
        assert_eq!(mapper.tx_to_genomic(&TxPos::new(10)).unwrap(), Some(3000));
    }

    #[test]
    fn test_genomic_to_cds() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // genomic 1004 -> tx 5 -> c.1 (cds_start = 5)
        let cds_pos = mapper.genomic_to_cds(1004).unwrap().unwrap();
        assert_eq!(cds_pos.base, 1);
        assert!(!cds_pos.utr3);
    }

    #[test]
    fn test_cds_to_genomic() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // c.1 -> tx 5 -> genomic 1004
        assert_eq!(mapper.cds_to_genomic(&CdsPos::new(1)).unwrap(), Some(1004));
    }

    #[test]
    fn test_no_genomic_coords_error() {
        let tx = make_test_transcript(); // No genomic coords
        let mapper = CoordinateMapper::new(&tx);

        assert!(mapper.genomic_to_tx(1000).is_err());
        assert!(mapper.tx_to_genomic(&TxPos::new(1)).is_err());
    }

    /// A transcript whose exons do not tile transcript space: transcript
    /// positions 1, 95 and 194 fall in no exon.
    ///
    /// Real cdot data does carry such holes — 58 of 474,818 GRCh38 multi-exon
    /// builds (0.012%), the smallest of them 23 bases and none of them one base
    /// — so the shape is not itself malformed. The **sizes** here are synthetic:
    /// one base per junction, which is a shape upstream never produces. That is
    /// deliberate and is now the point of the helper. Under the flat sequence
    /// axis the gap size cannot affect a `c.`/`n.` answer at all, so a hole of
    /// the least plausible size is the sharpest available check that nothing
    /// reads the exon table on this path.
    ///
    /// It does **not** describe cdot's *encoding*. This helper previously
    /// claimed to "simulate how cdot encodes transcripts with virtual intron
    /// positions"; cdot has no such encoding. Its `cds_start_i`/`cds_end_i` are
    /// 1-based inclusive coordinates on the *spliced* transcript, and a hole
    /// means transcript bases that align to no genomic exon — not extra
    /// coordinates that `c.` numbering skips.
    ///
    /// The false claim was load-bearing: it was the stated justification for the
    /// one-base-per-junction holes that
    /// `tests/fixtures/sequences/normalization_transcripts.json` carried in every
    /// multi-exon record, which drifted `NM_000492.4:c.1520` ten bases 3' (#1619).
    /// `tests/it/normalization_transcripts_exon_contract.rs` now guards that
    /// fixture; the three tests below use this helper to pin that the walk those
    /// holes used to trigger is gone.
    fn make_transcript_with_gaps() -> Transcript {
        // Transcript with 3 exons and gaps between them
        // Exon 1: tx 2-94 (93 bases)
        // Gap at tx 95
        // Exon 2: tx 96-193 (98 bases), CDS starts at tx 114
        // Gap at tx 194
        // Exon 3: tx 195-247 (53 bases)
        Transcript {
            cds_start_incomplete: false,
            id: "NM_GAPS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("A".repeat(244)), // Not used for coordinate tests
            cds_start: Some(114),            // 1-based CDS start
            cds_end: Some(247),              // 1-based CDS end
            exons: vec![
                Exon::new(1, 2, 94),    // Exon 1: tx 2-94
                Exon::new(2, 96, 193),  // Exon 2: tx 96-193 (gap at 95)
                Exon::new(3, 195, 247), // Exon 3: tx 195-247 (gap at 194)
            ],
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

    /// #1619: a `c.` position is `cds_start + N - 1` on the flat transcript.
    /// The exon holes at tx 95 and 194 must change nothing — `c.81` is tx 194,
    /// the very base that lies in a hole, because that base is in the
    /// transcript sequence and `c.` numbering counts it.
    #[test]
    fn cds_to_tx_ignores_a_transcript_coordinate_gap() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);
        // cds_start = 114.
        for (cds, expected_tx) in [
            (1, 114),
            (10, 123),
            (80, 193),
            (81, 194),
            (82, 195),
            (133, 246),
        ] {
            assert_eq!(
                mapper.cds_to_tx(&CdsPos::new(cds)).unwrap().base,
                expected_tx,
                "c.{cds} must be the flat offset; the exon walk answered {}",
                expected_tx + 1,
            );
        }
    }

    /// The inverse, on the same axis. tx 194 is a real transcript base and so
    /// has a `c.` number (`c.81`); it is only the *genome* that has nothing
    /// there.
    #[test]
    fn tx_to_cds_ignores_a_transcript_coordinate_gap() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);
        for (tx_base, expected_cds) in [
            (114, 1),
            (123, 10),
            (193, 80),
            (194, 81),
            (195, 82),
            (247, 134),
        ] {
            assert_eq!(
                mapper.tx_to_cds(&TxPos::new(tx_base)).unwrap().base,
                expected_cds,
                "n.{tx_base} must be the flat CDS offset",
            );
        }
    }

    #[test]
    fn cds_tx_round_trip_across_a_transcript_coordinate_gap() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

        for cds_pos in [1, 10, 80, 81, 82, 100, 133] {
            let tx_pos = mapper.cds_to_tx(&CdsPos::new(cds_pos)).unwrap();
            let back = mapper.tx_to_cds(&tx_pos).unwrap();
            assert_eq!(
                back.base, cds_pos,
                "Roundtrip failed for CDS {}: got tx {}, back to CDS {}",
                cds_pos, tx_pos.base, back.base
            );
        }
    }

    fn make_cigar_insertion_transcript() -> Transcript {
        Transcript {
            cds_start_incomplete: false,
            id: "NM_INS.1".to_string(),
            gene_symbol: Some("INS".to_string()),
            strand: Strand::Plus,
            sequence: Some("A".repeat(40)),
            cds_start: Some(1),
            cds_end: Some(40),
            exons: vec![Exon::new(1, 1, 40)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            // M5 I3 M32: 3 inserted transcript bases at tx 6,7,8.
            exon_cigars: vec![Some(vec![
                CigarOp::Match(5),
                CigarOp::Insertion(3),
                CigarOp::Match(32),
            ])],
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn cds_to_tx_is_transcript_native_across_cigar_insertion() {
        // #944: c./n. numbering counts inserted transcript bases, so CDS->tx is
        // naive with NO CIGAR adjustment. cds_start=1 => c.N -> tx N (1-based).
        let tx = make_cigar_insertion_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // Before the insertion — unaffected.
        assert_eq!(mapper.cds_to_tx(&CdsPos::new(3)).unwrap().base, 3);
        // After the 3bp insertion — must NOT shift by +3 (bug produced 13).
        assert_eq!(mapper.cds_to_tx(&CdsPos::new(10)).unwrap().base, 10);
    }

    #[test]
    fn cds_tx_round_trip_across_cigar_insertion() {
        // tx_to_cds must invert cds_to_tx with no CIGAR adjustment (#944).
        let tx = make_cigar_insertion_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // Round-trip a post-insertion position: c.10 -> tx 10 -> c.10.
        let tx_pos = mapper.cds_to_tx(&CdsPos::new(10)).unwrap();
        assert_eq!(tx_pos.base, 10);
        assert_eq!(mapper.tx_to_cds(&tx_pos).unwrap().base, 10);
        // Direct n->c: tx 10 -> c.10 (bug produced c.7).
        assert_eq!(mapper.tx_to_cds(&TxPos::new(10)).unwrap().base, 10);
    }

    /// `chromosome()` had no assertion anywhere, so `-> None`, `-> Some("")`
    /// and `-> Some("xyzzy")` all survived.
    #[test]
    fn chromosome_reports_the_transcript_contig() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);
        assert_eq!(mapper.chromosome(), Some("chr1"));
    }

    /// `get_intron_position` has one production caller (in `src/vcf/to_hgvs.rs`)
    /// and no test reached it, so the whole-function `-> None` survived.
    ///
    /// Exon 1 ends at genomic 1009 (tx 10); genomic 1012 is 3 bases into
    /// intron 1.
    #[test]
    fn get_intron_position_reports_an_intronic_offset() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        let intron = mapper
            .get_intron_position(1012)
            .expect("genomic 1012 lies in intron 1");
        assert_eq!(intron.intron_number, 1, "genomic 1012 lies in intron 1");
        assert_eq!(intron.offset, 3, "3 bases past the exon 1 boundary");
        assert_eq!(intron.tx_boundary_pos, 10, "exon 1 ends at tx 10");

        // An exonic position is not intronic.
        assert!(mapper.get_intron_position(1005).is_none());
    }

    /// `fold_non_intronic_utr_offset` folds only 3'UTR and negative-base
    /// (5'UTR) positions. `base == 0` is neither of those — it is
    /// `CDS_BASE_UNKNOWN` (in `src/hgvs/location.rs`), the `c.?`
    /// unknown-position sentinel, not an arithmetic edge of the UTR check.
    /// It must still decline to fold, and doing so also happens to pin
    /// `pos.base < 0` against a `<= 0` mutant, since 0 is where the two
    /// diverge. Base 0 is the interesting input because this gate and
    /// `cds_to_tx` reach the same verdict by different routes: this gate
    /// classifies 5'UTR as `pos.base < 0`, so base 0 is simply not a UTR
    /// position to it, while `cds_to_tx` refuses base 0 outright as a value
    /// that names no nucleotide. Neither computes on it, and this test pins
    /// that the gate declines on its own terms rather than by borrowing
    /// `cds_to_tx`'s refusal.
    #[test]
    fn fold_offset_declines_unknown_position() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        let pos = CdsPos::unknown(Some(5));

        // The decline happens at the `in_utr` gate, three lines before
        // `cds_to_tx` is even called: `pos.utr3 || pos.base < 0` is
        // `false || (0 < 0)` = `false` for base 0 (`CDS_BASE_UNKNOWN`), so
        // `fold_non_intronic_utr_offset` returns `None` right there and the
        // `c.?` sentinel never reaches the conversion at all.
        // (`cds_to_tx`'s own refusal of this sentinel is pinned separately by
        // `cds_to_tx_refuses_an_unknown_position` below.)
        assert!(mapper.fold_non_intronic_utr_offset(&pos).is_none());
    }

    // -----------------------------------------------------------------------
    // Unvalidated CDS metadata, and coordinates that do not name a base.
    //
    // Everything below guards arithmetic that this file used to perform on a
    // value it had not established was a coordinate at all. The three shapes
    // are independent and each reaches a different expression, so they are
    // pinned separately rather than folded into one table-driven test.
    // -----------------------------------------------------------------------

    /// A transcript whose `cds_start` is `0` — the value `Transcript` accepts
    /// and never validates.
    ///
    /// `Transcript::new` takes `cds_start: Option<u64>` verbatim, and nothing on
    /// the construction path rejects `0`: `validate.rs`'s CDS check is
    /// out-of-band and report-only, so a provider serving corrupt or truncated
    /// CDS metadata reaches the mapper with it. `Transcript::utr5_length`
    /// already defends the same field with `saturating_sub(1)`, which is the
    /// in-tree precedent that `0` is representable here rather than a shape
    /// only a test can build.
    ///
    /// `0` is not merely out of range: `cds_start` is a **1-based** transcript
    /// coordinate, so `0` names no base of the transcript at all.
    fn make_zero_cds_start_transcript() -> Transcript {
        Transcript {
            cds_start_incomplete: false,
            id: "NM_ZEROCDS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGCCCAAAGGGTTTAGG".to_string()),
            cds_start: Some(0),
            cds_end: Some(18),
            exons: vec![Exon::new(1, 1, 18)],
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

    /// `get_amino_acid_at` computed its codon offset as
    /// `cds_start as usize - 1 + (position as usize - 1) * 3`, subtracting from
    /// an unvalidated `cds_start` in `usize`, where subtraction below zero is
    /// not representable.
    ///
    /// It panics in **both** build profiles, which is why this test is worth
    /// running under `--release` as well: debug aborts on the subtraction
    /// itself (`attempt to subtract with overflow`), and release wraps it to
    /// `usize::MAX` and aborts three lines later on the slice
    /// (`start byte index 18446744073709551615 is out of bounds`). The release
    /// failure is the worse of the two because it moves — the reported location
    /// is the string index, which names neither the bad input nor the
    /// subtraction that produced it.
    ///
    /// Reached through the ordinary public path `cds_to_protein`, which is what
    /// `python.rs`'s `c_to_p` calls.
    #[test]
    fn cds_to_protein_declines_a_zero_cds_start_instead_of_underflowing() {
        let tx = make_zero_cds_start_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // The protein *number* is arithmetic on the `c.` position alone and is
        // still derivable, so the conversion succeeds. What must not happen is
        // that resolving its residue underflows: with no usable `cds_start`
        // there is no codon to read, which is the `None` this method already
        // returns for a coordinate-only transcript.
        let prot = mapper
            .cds_to_protein(&CdsPos::new(1))
            .expect("the protein number does not depend on cds_start being usable");
        assert_eq!(prot.number, 1);
        assert_eq!(
            prot.aa,
            AminoAcid::Xaa,
            "a cds_start of 0 names no transcript base, so no residue is resolvable",
        );
    }

    /// The same guard on the second subtrahend. `position` is a 1-based protein
    /// number, so `0` names no residue; it is private and today unreachable
    /// from `cds_to_protein` (which refuses `base < 1` first), and it is pinned
    /// anyway because the two subtractions are one expression and a caller
    /// added later would not know that.
    #[test]
    fn a_zero_protein_number_resolves_no_residue() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        assert_eq!(mapper.get_amino_acid_at(0), None);
        // The neighbouring live value still resolves, so the guard is not just
        // returning `None` for everything: `NM_TEST.1` has cds_start = 6 and
        // sequence `AAAAATGCCC…`, so codon 1 is `TGC`.
        assert_eq!(mapper.get_amino_acid_at(1), Some(AminoAcid::Cys));
    }

    /// `c.?` declares the position unknown, so no transcript coordinate can be
    /// derived from it. `cds_to_tx` answered `Ok(cds_start)` — `c.1`'s own
    /// coordinate — because `CDS_BASE_UNKNOWN` is the integer `0`, and `0` took
    /// the `base < 1` 5'UTR branch with a zero displacement.
    ///
    /// The doctrine is already written in this tree, twice over. On the sibling
    /// predicate `CdsPos::has_unknown_offset` (`src/hgvs/location.rs`), for an
    /// unknown *offset*: "no coordinate can be derived from them at all", and
    /// "callers doing coordinate arithmetic MUST check this first". And at two
    /// sites in `src/project/projector.rs`, which reject `is_unknown()` before
    /// calling down here precisely because — in that code's own words —
    /// "`base <= 0` would misclassify it as 5'UTR and return `Ok`". That is
    /// this defect, described by its callers; the guard belongs at the site
    /// they were working around.
    #[test]
    fn cds_to_tx_refuses_an_unknown_position() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        assert!(
            mapper.cds_to_tx(&CdsPos::unknown(None)).is_err(),
            "c.? must not be answered with c.1's transcript coordinate",
        );
        // An offset does not make the base known: `c.?+5` is still anchored on
        // a position that names nothing.
        assert!(mapper.cds_to_tx(&CdsPos::unknown(Some(5))).is_err());
    }

    /// The neighbours on both sides of the refused value still convert, so the
    /// guard is `base == 0` and not a widened 5'UTR refusal. `NM_TEST.1` has
    /// cds_start = 6.
    #[test]
    fn cds_to_tx_still_answers_the_positions_either_side_of_the_sentinel() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        assert_eq!(mapper.cds_to_tx(&CdsPos::new(1)).unwrap().base, 6);
        assert_eq!(mapper.cds_to_tx(&CdsPos::new(-1)).unwrap().base, 5);
        assert_eq!(mapper.cds_to_tx(&CdsPos::utr3(1)).unwrap().base, 36);
    }

    /// #1690: `cds_to_tx`'s addition is unchecked, so a `c.` coordinate near
    /// `i64::MAX` panics with `attempt to add with overflow` instead of being
    /// refused — in debug, and wrapping to a nonsense negative coordinate in
    /// release.
    ///
    /// This is the same defect class as the two above and lives in the same
    /// expression, so it is guarded here rather than left behind. The issue's
    /// own reproducer is three tests in
    /// `tests/it/issue_1487_canonical_window_overflow.rs`, which reach it
    /// through `hgvs_to_spdi` under `FERRO_ASSERT_SEQUENCE=1`; the unit test
    /// below pins the refusal unconditionally, in every profile and with no
    /// oracle armed.
    ///
    /// Note the issue names `cds_to_tx_exon_aware` at `mapper.rs:113`. That
    /// function no longer exists — #1735 deleted the exon walk from this axis —
    /// and the overflow moved to the flat conversion that replaced it, which
    /// the recorded panic locates at `mapper.rs:114`.
    #[test]
    fn cds_to_tx_refuses_an_extreme_coordinate_instead_of_overflowing() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // The CDS body and the 3'UTR both add a positive base to a positive
        // transcript coordinate, so both overflow at the top of the range.
        assert!(mapper.cds_to_tx(&CdsPos::new(i64::MAX)).is_err());
        assert!(mapper.cds_to_tx(&CdsPos::utr3(i64::MAX)).is_err());

        // The 5'UTR branch is the one that CANNOT overflow, and saying so is
        // worth more than a third assertion: `cds_start` is a non-negative
        // `i64` and the branch is taken only for `base < 1`, so the sum lies in
        // `i64::MIN..=cds_start` by construction. `c.-9223372036854775808`
        // therefore converts, to a coordinate no transcript has — which is the
        // callers' question, and they already ask it
        // (`noncoding_from_coding` refuses `n.` below 1, `resolve_cds_to_tx`
        // bounds both endpoints against the transcript length). The
        // `checked_add` on that branch is kept for symmetry and is deliberately
        // not claimed as tested.
        assert!(mapper.cds_to_tx(&CdsPos::new(i64::MIN)).unwrap().base < 0);
    }

    // -----------------------------------------------------------------------
    // The other two shapes that carry CDS base 0.
    //
    // `is_unknown()` is `base == 0 && !utr3 && special.is_none()`, so a guard
    // written on it recognises exactly one of the three shapes that carry the
    // integer `0`. Each of the two it misses gets its own test, so a future
    // regression names WHICH arm re-opened rather than reporting "the base-0
    // guard broke".
    //
    // Each asserts the COLLAPSE, not merely that some error occurs: it pairs
    // the no-base input with the ordinary coordinate that used to receive the
    // same answer, and requires the two to be told apart. A test that only
    // asserted `is_err()` would pass against a guard that had accidentally
    // widened to refuse the ordinary neighbour too.
    // -----------------------------------------------------------------------

    /// `c.*0` (base `0`, `utr3 = true`).
    ///
    /// `is_unknown()` is false for it — the `utr3` conjunct fails — so it took
    /// the 3'UTR arm with a zero displacement and came back as `cds_end`, which
    /// is the transcript coordinate of the LAST CDS BASE. On `NM_TEST.1`
    /// (cds_start 6, cds_end 35) that made `c.*0` and `c.30` both answer tx 35.
    ///
    /// The `*` zone is numbered from `*1`: `background/numbering.md:30` gives
    /// the nucleotides 3' of the stop codon as "numbered `c.*1`, `c.*2`,
    /// `c.*3`, etc., going further downstream", and `:31` adds "there is no
    /// nucleotide `c.0`". So `*0` is outside the enumeration, and answering it
    /// with a nucleotide that already has a plain `c.` spelling is worse than
    /// refusing.
    ///
    /// Reachable only by a **Rust** caller constructing a `CdsPos` — the fields
    /// are `pub`. Not from a string: `parse_hgvs` refuses `c.*0del` outright
    /// (measured: `Failed to parse variant … code: Verify`). And not from the
    /// Python bindings either, which is a change of reach that arrived with the
    /// rebase rather than with this guard: #1741 (`b3568679`, on `main` but not
    /// on this branch's original base) added `reject_zero_base`, and `c_to_g`,
    /// `c_to_p`, `c_to_n` and `n_to_c` all call it before converting. `c_to_n`
    /// passes its own `utr3` through and `reject_zero_base` keys on the base
    /// alone, so `c_to_n(tx, 0, utr3=True)` — this shape — is refused at the
    /// boundary. Do not restate the older "Rust API and the Python
    /// `c_to_n`/`c_to_g` entry points" reach for any base-0 shape; it was true
    /// when written and #1741 closed the Python half.
    #[test]
    fn cds_to_tx_refuses_a_zero_three_prime_utr_position() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        let star_zero = CdsPos {
            base: 0,
            offset: None,
            utr3: true,
            special: None,
        };
        assert!(
            !star_zero.is_unknown(),
            "the shape under test is precisely the one `is_unknown()` does not see",
        );

        // The collapse: `c.30` is a real coordinate and must keep its answer,
        // and `c.*0` must no longer be given the same one.
        let last_cds_base = mapper
            .cds_to_tx(&CdsPos::new(30))
            .expect("c.30 is the last CDS base and must still convert");
        assert_eq!(last_cds_base.base, 35);
        assert!(
            mapper.cds_to_tx(&star_zero).is_err(),
            "c.*0 names no nucleotide and must not be answered with c.30's \
             transcript coordinate ({})",
            last_cds_base.base,
        );

        // The refusal is about the zero, not about the 3'UTR zone: `c.*1`, the
        // first nucleotide the zone actually numbers, still converts.
        assert_eq!(mapper.cds_to_tx(&CdsPos::utr3(1)).unwrap().base, 36);
    }

    /// `c.pter` / `c.qter` / `c.cen` (base `0`, `special = Some(..)`).
    ///
    /// `is_unknown()` is false for these — the `special.is_none()` conjunct
    /// fails — so they took the `base < 1` 5'UTR arm with a zero displacement
    /// and came back as `cds_start`, which is `c.1`'s own answer. On
    /// `NM_TEST.1` all three markers and `c.1` answered tx 6.
    ///
    /// `pter` and `qter` are opposite ends of a chromosome, so the defect is
    /// not a near miss in either direction: the two markers agreed with each
    /// other as well as with an ordinary coding position. End to end through
    /// `hgvs_to_spdi` the three descriptions `c.pterdel`, `c.qterdel` and
    /// `c.1del` produced one bit-identical SPDI triple; that half is pinned in
    /// `tests/it/issue_1487_canonical_window_overflow.rs`, since it needs a
    /// provider.
    ///
    /// A marker is a *genomic* landmark. Where it can be resolved at all it is
    /// resolved by `project_cds_terminus_to_parent` against a named genomic
    /// parent, which declines outright without one — so a marker arriving here,
    /// on the transcript's own flat axis, has no coordinate to be given.
    #[test]
    fn cds_to_tx_refuses_a_terminus_marker() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        let cds_first_base = mapper
            .cds_to_tx(&CdsPos::new(1))
            .expect("c.1 must still convert");
        assert_eq!(cds_first_base.base, 6);

        for marker in [CdsPos::pter(), CdsPos::qter(), CdsPos::cen()] {
            assert!(
                !marker.is_unknown(),
                "c.{marker} is precisely a shape `is_unknown()` does not see",
            );
            assert!(
                mapper.cds_to_tx(&marker).is_err(),
                "c.{marker} names no position on this transcript's CDS axis and \
                 must not be answered with c.1's transcript coordinate ({})",
                cds_first_base.base,
            );
        }

        // The marker is refused on its shape, not only on its base, so a
        // marker carrying some other integer cannot slip past. The public
        // constructors all set `0` today; keying on the shape as well means a
        // later one that does not cannot silently re-open this.
        let marker_with_a_base = CdsPos {
            base: 7,
            offset: None,
            utr3: false,
            special: Some(SpecialPosition::Qter),
        };
        assert!(
            mapper.cds_to_tx(&marker_with_a_base).is_err(),
            "a terminus marker is not a displacement from the CDS start, whatever \
             base accompanies it",
        );
    }
}

#[cfg(test)]
mod intronic_debug_tests {
    use super::*;
    use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
    use std::sync::OnceLock;

    /// Create a transcript mimicking NM_003742.4 exon structure (first 14 exons)
    /// Based on annotation file which uses 0-based half-open coordinates:
    /// CDS: start=127, end=4093
    /// Exon 12: [1324, 1435)
    /// Exon 13: [1435, 1561)
    /// Exon 14: [1561, 1765)
    fn create_nm003742_like_transcript() -> Transcript {
        // After multi_fasta.rs conversion (0-based to 1-based):
        // start = original_start + 1, end = original_end (half-open becomes inclusive)
        let exons = vec![
            Exon {
                number: 1,
                start: 1,
                end: 100,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 2,
                start: 101,
                end: 203,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 3,
                start: 204,
                end: 225,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 4,
                start: 226,
                end: 277,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 5,
                start: 278,
                end: 516,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 6,
                start: 517,
                end: 604,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 7,
                start: 605,
                end: 738,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 8,
                start: 739,
                end: 910,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 9,
                start: 911,
                end: 1035,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 10,
                start: 1036,
                end: 1210,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 11,
                start: 1211,
                end: 1324,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 12,
                start: 1325,
                end: 1435,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 13,
                start: 1436,
                end: 1561,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 14,
                start: 1562,
                end: 1765,
                genomic_start: None,
                genomic_end: None,
            },
        ];

        Transcript {
            cds_start_incomplete: false,
            id: "NM_003742.4".to_string(),
            gene_symbol: Some("ABCB11".to_string()),
            strand: Strand::Minus,
            sequence: None,       // Not needed for coordinate tests
            cds_start: Some(128), // 127 + 1 = 128 (1-based)
            cds_end: Some(4093),  // Same (half-open end = inclusive 1-based)
            exons,
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn debug_nm003742_intronic() {
        let transcript = create_nm003742_like_transcript();

        eprintln!("\n=== NM_003742.4-like Transcript ===");
        eprintln!("cds_start: {:?}", transcript.cds_start);
        eprintln!("cds_end: {:?}", transcript.cds_end);
        eprintln!("num exons: {}", transcript.exons.len());

        // Check for gaps
        let mut sorted_exons: Vec<_> = transcript.exons.iter().collect();
        sorted_exons.sort_by_key(|e| e.start);

        eprintln!("\nExons 11-14:");
        for e in &sorted_exons[10..14] {
            eprintln!("  Exon {}: tx {}..{}", e.number, e.start, e.end);
        }

        let has_gaps = sorted_exons.windows(2).any(|w| w[0].end + 1 != w[1].start);
        eprintln!("\nHas gaps in tx coords: {}", has_gaps);

        // Calculate expected tx for c.1435
        let cds_start = transcript.cds_start.unwrap();
        let expected_tx = cds_start as i64 + 1434;
        eprintln!("\nFor c.1435:");
        eprintln!("  cds_start = {}", cds_start);
        eprintln!("  expected tx = {} + 1434 = {}", cds_start, expected_tx);

        // Use the mapper
        let mapper = CoordinateMapper::new(&transcript);
        let cds_pos = CdsPos {
            base: 1435,
            offset: None,
            utr3: false,
            special: None,
        };
        let tx_result = mapper.cds_to_tx(&cds_pos);
        eprintln!("  actual tx from mapper: {:?}", tx_result);

        // Verify the calculation
        assert!(tx_result.is_ok());
        let tx_pos = tx_result.unwrap();

        // c.1435 should map to tx 128 + 1434 = 1562
        assert_eq!(
            tx_pos.base, 1562,
            "c.1435 should map to tx 1562, but got tx {}",
            tx_pos.base
        );
    }

    /// A minimal cdot document spelled in cdot's own RAW conventions: genomic
    /// bounds 0-based half-open (`alt_start` inclusive, `alt_end` exclusive)
    /// and cDNA bounds 1-based **inclusive**.
    ///
    /// | exon | raw genomic    | raw tx (1-based incl) | expected internal tx |
    /// |------|----------------|-----------------------|----------------------|
    /// | 1    | `[5000, 5100)` | `1..100`              | `[0, 100)`           |
    /// | 2    | `[6000, 6050)` | `101..150`            | `[100, 150)`         |
    ///
    /// Two exons deliberately. The second is what makes the *half-open* half of
    /// the claim testable at all: under 0-based half-open storage exon 2's
    /// `tx_start` equals exon 1's `tx_end` exactly, where 1-based-inclusive
    /// storage would leave it one higher.
    ///
    /// `start_codon`/`stop_codon` are deliberately ABSENT. `from_genome_build`
    /// copies those through verbatim when they are present, so a fixture
    /// carrying them would bypass the conversion under test and assert nothing
    /// about the coordinate basis. Omitting them forces the genomic-CDS
    /// fallback — the arm that has to reconcile raw genomic bounds against the
    /// already-converted exon table, which is where #742's off-by-one lived.
    const RAW_CDOT_BASIS_FIXTURE: &str = r#"{
        "transcripts": {
            "NM_BASIS_TEST.1": {
                "gene_name": "BASISTEST",
                "genome_builds": {
                    "GRCh38": {
                        "contig": "NC_000001.11",
                        "strand": "+",
                        "cds_start": 5010,
                        "cds_end": 6041,
                        "exons": [
                            [5000, 5100, 1, 1, 100, "M100"],
                            [6000, 6050, 2, 101, 150, "M50"]
                        ]
                    }
                }
            }
        }
    }"#;

    /// cdot tx coordinates are stored 0-based half-open (#742).
    ///
    /// On load, `from_genome_build` converts cdot's raw 1-based-inclusive tx
    /// bounds into the engine's HGVS convention: tx 0-based half-open, the
    /// first exon starting at 0 rather than 1. Violating that convention is a
    /// guaranteed panic on the first exon, so the guard must not be able to
    /// stop running quietly.
    ///
    /// This runs on every machine, against a programmatically-built fixture.
    /// Its predecessor asserted only when an operator-created
    /// `benchmark-output/` symlink was present, and returned silently
    /// otherwise — i.e. it asserted nothing on essentially every run,
    /// including CI. The real-data cross-check lives on in
    /// `cdot_tx_basis_cross_checked_against_real_cdot` below, which is an
    /// addition to this test rather than a substitute for it.
    #[test]
    fn cdot_tx_coordinates_are_0_based_half_open() {
        use crate::data::cdot::CdotMapper;

        let cdot = CdotMapper::from_reader_with_build(RAW_CDOT_BASIS_FIXTURE.as_bytes(), "GRCh38")
            .expect("basis fixture should parse");
        let tx = cdot
            .get_transcript("NM_BASIS_TEST.1")
            .expect("basis fixture defines NM_BASIS_TEST.1");

        assert_eq!(tx.exons.len(), 2, "fixture should yield both exons");

        // 0-based: raw tx 1..100 (1-based inclusive) becomes [0, 100).
        let (first_start, first_end) = (tx.exons[0][2], tx.exons[0][3]);
        assert_eq!(
            (first_start, first_end),
            (0, 100),
            "exon 1 raw tx 1..100 (1-based inclusive) should convert to 0-based \
             half-open [0, 100), got [{first_start}, {first_end})"
        );

        // 0-based: raw tx 101..150 becomes [100, 150).
        let (second_start, second_end) = (tx.exons[1][2], tx.exons[1][3]);
        assert_eq!(
            (second_start, second_end),
            (100, 150),
            "exon 2 raw tx 101..150 (1-based inclusive) should convert to 0-based \
             half-open [100, 150), got [{second_start}, {second_end})"
        );

        // Half-open: adjoining exons meet at one shared value. Under
        // 1-based-inclusive storage exon 2 would start at exon 1's end + 1.
        assert_eq!(
            second_start, first_end,
            "half-open tx bounds should adjoin: exon 2 tx_start ({second_start}) must equal \
             exon 1 tx_end ({first_end}); a difference of 1 means the bounds are stored \
             1-based inclusive"
        );

        // CDS start is 0-based inclusive in transcript space. The fixture's
        // 5'UTR is 10 bases (raw genomic CDS start 5010 against an exon opening
        // at raw genomic 5000), so c.1 sits at tx offset 10 — it would be 11 if
        // the CDS bound were stored 1-based.
        assert_eq!(
            tx.cds_start,
            Some(10),
            "cds_start should be the 0-based tx offset of c.1 (10, after a 10-base 5'UTR)"
        );

        // `cds_end` is deliberately not asserted here. The genomic-CDS fallback
        // applies a further `+ 1` to the mapped end whose correctness depends on
        // whether cdot's raw genomic `cds_end` is inclusive or exclusive — the
        // struct doc says exclusive while the arithmetic implies inclusive. That
        // is a separate question from the coordinate basis this test pins, and
        // pinning today's value here would cement whichever reading is wrong.
    }

    /// Why the real-cdot lookup did or did not produce a path.
    ///
    /// The three non-`Found` arms are kept apart on purpose: they call for
    /// three different messages, and collapsing them into one `None` is how a
    /// skip notice comes to state a condition that does not hold. Telling an
    /// operator whose `benchmark-output/` symlink dangles or whose reference
    /// directory is unreadable to "symlink benchmark-output/ to a prepared
    /// reference" is wrong advice — they already did.
    enum CdotLookup {
        /// A `cdot-*.refseq.GRCh38.json` was found at this path.
        Found(std::path::PathBuf),
        /// The directory does not exist. The default-checkout and CI case, and
        /// also what a dangling symlink reports (both surface as `NotFound`).
        Absent,
        /// The directory exists but could not be read — permissions, or an I/O
        /// error. A wired-up-but-broken reference, not an absent one.
        Unreadable(std::io::Error),
        /// The directory was read and holds no matching cdot JSON.
        NoMatch,
    }

    /// Locate a RefSeq GRCh38 cdot JSON under `dir` without pinning a cdot
    /// release.
    ///
    /// Pattern-matching rather than naming a version is the point: the caller's
    /// predecessor hardcoded `cdot-0.2.32.refseq.GRCh38.json`, so regenerating
    /// the prepared reference at any newer cdot release (0.2.33 is current)
    /// would have returned it to never running, silently.
    ///
    /// When several releases sit side by side the pick is the
    /// **lexicographically last** filename, which is deterministic but is *not*
    /// "the newest release" — `cdot-0.2.9.…` sorts after `cdot-0.2.32.…`.
    /// Any real cdot answers the question this test asks, so the ordering only
    /// has to be stable, not semantic.
    fn locate_refseq_grch38_cdot(dir: &std::path::Path) -> CdotLookup {
        let entries = match std::fs::read_dir(dir) {
            Ok(entries) => entries,
            Err(err) if err.kind() == std::io::ErrorKind::NotFound => return CdotLookup::Absent,
            Err(err) => return CdotLookup::Unreadable(err),
        };
        let mut matches: Vec<std::path::PathBuf> = entries
            .filter_map(Result::ok)
            .map(|entry| entry.path())
            .filter(|path| {
                path.file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| {
                        name.starts_with("cdot-") && name.ends_with(".refseq.GRCh38.json")
                    })
            })
            .collect();
        matches.sort();
        match matches.pop() {
            Some(path) => CdotLookup::Found(path),
            None => CdotLookup::NoMatch,
        }
    }

    /// Cross-check the 0-based basis against a real prepared reference.
    ///
    /// Kept as an *addition* to `cdot_tx_coordinates_are_0_based_half_open`,
    /// never as a substitute. It is not a gate — `benchmark-output/` is an
    /// operator-created symlink to a prepared reference, absent in a default
    /// checkout and in CI — so it names the specific condition it hit and
    /// returns.
    ///
    /// **The skip notice does not reach a passing run's console.** Both
    /// `cargo test` and `cargo nextest run` capture a test's stdout/stderr and
    /// replay it only for failures, so the `eprintln!` below is visible only
    /// under `--nocapture` / `--no-capture` (or nextest's
    /// `--success-output immediate`). Do not read it as making the skip
    /// self-announcing: what stops this test from being the silent no-op its
    /// predecessor was is that the basis is now guarded unconditionally by
    /// `cdot_tx_coordinates_are_0_based_half_open`, which needs no reference at
    /// all. Turning the skip itself into a hard failure would need an opt-in
    /// switch (an env var an operator sets when they expect the reference to be
    /// wired), which is deliberately not added here.
    ///
    /// **Read what this does and does not cover.** `CdotMapper::load` prefers a
    /// sibling `.rkyv` cache over the JSON, and a prepared reference normally
    /// ships one, so on a *warm* reference this deserializes coordinates that
    /// were converted at prepare time rather than re-running
    /// `from_genome_build`. It therefore guards the prepared *artifact's*
    /// stored basis — a mis-prepared or stale reference — and NOT the
    /// conversion. Measured, not assumed: inverting the conversion to emit
    /// 1-based leaves this test green while
    /// `cdot_tx_coordinates_are_0_based_half_open` fails.
    ///
    /// That warm-reference reading is a description of the common case, not a
    /// guarantee: `load` falls back to parsing the JSON whenever the sibling
    /// `.rkyv` is missing *or older than the JSON* (its `cache_is_fresh`
    /// check), and the fallback then writes a fresh multi-hundred-MB `.rkyv`
    /// into the operator's reference directory as a side effect. So on a cold
    /// or stale reference this test does pay the ~489 MiB parse, and re-running
    /// it is cheap only from the second run on. That is the reason not to
    /// *force* the JSON arm, and it is also why the fixture above — which
    /// covers the conversion on every machine, in microseconds — is the
    /// primary guard rather than this one.
    ///
    /// The fixture's own faithfulness was checked against this same file rather
    /// than assumed: real cdot spells raw exon tx bounds 1-based **inclusive**
    /// (adjacent exons read `…4521]`, `[4522…`) with the first exon's raw
    /// `tx_start` at 1 — exactly the spelling `RAW_CDOT_BASIS_FIXTURE` uses.
    ///
    /// It asserts a structural invariant over a sample of real transcripts
    /// rather than one accession's `cds_start`. The pinned
    /// `NM_003742.4 cds_start == 127` it replaces was live external data that
    /// upstream can revise, which would have failed this test for a reason
    /// having nothing to do with the coordinate basis.
    ///
    /// Exon-to-exon tx contiguity is deliberately NOT asserted here: real cdot
    /// carries genuine transcript-coordinate holes (58 of 474,818 GRCh38
    /// multi-exon builds, 23–2718 bases), so contiguity is a property of the
    /// fixture's geometry and not of cdot at large. `tx_start == 0` on the
    /// first exon is universal, and is the discriminating one anyway — it reads
    /// 1 under 1-based storage.
    #[test]
    fn cdot_tx_basis_cross_checked_against_real_cdot() {
        use crate::data::cdot::CdotMapper;
        use std::path::PathBuf;

        let cdot_dir = PathBuf::from("benchmark-output/cdot");
        // Each arm states the condition it actually hit. A single "no cdot file
        // found" message would misdiagnose the two broken-setup cases as an
        // absent reference and hand the operator a fix they have already applied.
        const SKIP_PREAMBLE: &str = "cdot_tx_basis_cross_checked_against_real_cdot: skipping — \
             this is a cross-check only; the basis itself is guarded unconditionally by \
             cdot_tx_coordinates_are_0_based_half_open.";
        let cdot_path = match locate_refseq_grch38_cdot(&cdot_dir) {
            CdotLookup::Found(path) => path,
            CdotLookup::Absent => {
                eprintln!(
                    "{SKIP_PREAMBLE} {} does not exist (a dangling benchmark-output symlink \
                     reports the same). To run it, symlink benchmark-output/ to a prepared \
                     reference.",
                    cdot_dir.display()
                );
                return;
            }
            CdotLookup::Unreadable(err) => {
                eprintln!(
                    "{SKIP_PREAMBLE} {} exists but could not be read: {err}. The reference is \
                     wired up and BROKEN, not absent — re-linking benchmark-output/ will not \
                     fix this.",
                    cdot_dir.display()
                );
                return;
            }
            CdotLookup::NoMatch => {
                eprintln!(
                    "{SKIP_PREAMBLE} {} was read but holds no cdot-*.refseq.GRCh38.json. A \
                     prepared reference should carry one; check that prepare completed and \
                     that the RefSeq GRCh38 cdot was not stored gzipped-only.",
                    cdot_dir.display()
                );
                return;
            }
        };

        let cdot = CdotMapper::load(&cdot_path).unwrap_or_else(|err| {
            panic!("failed to load real cdot at {}: {err}", cdot_path.display())
        });

        // A sample is enough: the invariant is per-transcript, and loading the
        // full RefSeq cdot already dominates this test's cost.
        //
        // Sorted before sampling. `transcripts_on_contig` hands back
        // `contig_index` order, which on the JSON arm is `HashMap` iteration
        // order and therefore differs from process to process — an unsorted
        // `take` would compare a different 200 transcripts on every run, so a
        // failure would not reproduce and the covered set would be unknowable.
        const SAMPLE_SIZE: usize = 200;
        let mut accessions = cdot.transcripts_on_contig("NC_000001.11");
        accessions.sort_unstable();
        let attempted = accessions.len().min(SAMPLE_SIZE);
        let mut checked = 0usize;
        let mut unresolved: Vec<&str> = Vec::new();
        let mut exonless: Vec<&str> = Vec::new();
        for accession in accessions.iter().take(SAMPLE_SIZE) {
            let Some(tx) = cdot.get_transcript(accession) else {
                unresolved.push(*accession);
                continue;
            };
            let Some(first_exon) = tx.exons.first() else {
                exonless.push(*accession);
                continue;
            };
            assert_eq!(
                first_exon[2],
                0,
                "{accession}: first exon tx_start should be 0 (0-based half-open) in {}, got {}",
                cdot_path.display(),
                first_exon[2]
            );
            checked += 1;
        }

        // A zero denominator would let this report success having compared
        // nothing — the same failure mode the silent `return` had.
        assert!(
            checked > 0,
            "no transcripts checked from {} — expected transcripts on NC_000001.11, \
             so this run proved nothing",
            cdot_path.display()
        );

        // ...and a near-zero denominator would let it report success having
        // compared almost nothing, which `checked > 0` alone does not exclude:
        // 199 of 200 could fall through the two `continue` arms and the test
        // would still pass. Both arms are structurally unreachable — every
        // accession `transcripts_on_contig` yields is a key of the very
        // `transcripts` map `get_transcript` looks in first, and every cdot
        // genome build deserializes with a non-empty exon list — so a drop here
        // is information about the reference, not noise to absorb into the
        // count. Account for it rather than letting it thin the sample
        // silently.
        assert_eq!(
            checked,
            attempted,
            "compared {checked} of {attempted} sampled transcripts from {}: {} did not resolve \
             via get_transcript ({unresolved:?}) and {} carried an empty exon table \
             ({exonless:?}) — an unaccounted drop makes the comparison count meaningless",
            cdot_path.display(),
            unresolved.len(),
            exonless.len()
        );
    }
}
