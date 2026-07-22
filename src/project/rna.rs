//! Predicted RNA consequence surface: `c.`/`n.` → `r.(…)`.
//!
//! Numbering is CDS-relative (identical to `c.`); the 3′ shift and
//! range-reference resolution are performed in transcript-sequence space.
//!
//! Intronic positions have no `r.` representation of their own: an RNA
//! reference sequence is the *spliced* transcript, and the spec is explicit
//! that it therefore "can … **not be used** to describe variants affecting
//! these sequences" (`background/numbering.md`). An edit with an intronic
//! endpoint is only rendered when its exonic span is itself the consequence —
//! a whole-exon skip, which `RNA/splicing.md` spells with plain exonic
//! positions. Otherwise the axis declines rather than re-anchor onto the
//! neighbouring exonic base, which exists on the mature RNA and would silently
//! describe a different nucleotide. See `exonic_span`.
//! RNA rendering (T→U, lowercase, predicted `( )` wrapper) is handled by the
//! existing `RnaVariant` Display path; this module only builds the value.

use crate::convert::mapper::CoordinateMapper;
use crate::hgvs::edit::{InsertedPart, InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::{CdsInterval, RnaInterval};
use crate::hgvs::location::{CdsPos, RnaPos, TxPos};
use crate::hgvs::variant::{AlleleVariant, HgvsVariant, LocEdit, RnaVariant};
use crate::normalize::boundary::Boundaries;
use crate::normalize::config::ShuffleDirection;
use crate::normalize::shuffle::shuffle;
use crate::reference::transcript::Transcript;
use crate::sequence::reverse_complement;

/// Map a CDS position to an RNA position (identical numbering; drop `special`).
fn cds_pos_to_rna(p: &CdsPos) -> Option<RnaPos> {
    if p.special.is_some() {
        return None; // pter/qter not representable on the r. axis
    }
    Some(RnaPos {
        base: p.base,
        offset: p.offset,
        utr3: p.utr3,
    })
}

/// Map a CDS interval to an RNA interval (certain boundaries only).
fn cds_interval_to_rna(iv: &CdsInterval) -> Option<RnaInterval> {
    let start = cds_pos_to_rna(iv.start.inner()?)?;
    let end = cds_pos_to_rna(iv.end.inner()?)?;
    Some(RnaInterval::new(start, end))
}

/// 3′-shift a single deletion interval along the spliced transcript sequence,
/// in CDS coordinates. Exonic deletions only; the intronic/special/utr3 guard is
/// a safety fallback (`exonic_span` already stripped or declined every intronic
/// endpoint before this runs) — such intervals are returned unchanged.
fn shift_deletion_3prime(iv: &CdsInterval, transcript: &Transcript) -> Option<CdsInterval> {
    let start_cds = iv.start.inner()?;
    let end_cds = iv.end.inner()?;
    if start_cds.is_intronic()
        || end_cds.is_intronic()
        || start_cds.special.is_some()
        || end_cds.special.is_some()
        || start_cds.utr3
        || end_cds.utr3
    {
        return Some(iv.clone());
    }
    let seq = transcript.sequence.as_ref()?;
    let mapper = CoordinateMapper::new(transcript);
    let tx_start = mapper.cds_to_tx(start_cds).ok()?.base; // 1-based
    let tx_end = mapper.cds_to_tx(end_cds).ok()?.base; // 1-based
    if tx_start < 1 || tx_end < tx_start {
        return Some(iv.clone());
    }
    let start0 = (tx_start - 1) as u64; // 0-based inclusive
    let end0 = tx_end as u64; // 0-based exclusive
    let boundaries = Boundaries::new(0, seq.len() as u64);
    let r = shuffle(
        seq.as_bytes(),
        &[],
        start0,
        end0,
        &boundaries,
        ShuffleDirection::ThreePrime,
    );
    let new_tx_start = TxPos::new((r.start + 1) as i64); // back to 1-based
    let new_tx_end = TxPos::new(r.end as i64);
    let new_start = mapper.tx_to_cds(&new_tx_start).ok()?;
    let new_end = mapper.tx_to_cds(&new_tx_end).ok()?;
    Some(CdsInterval::new(new_start, new_end))
}

/// Reduce a CDS interval to the exonic span the `r.` axis may name, declining
/// (`None`) when no such span exists.
///
/// A non-zero intronic offset (`c.3675-45`, `c.831+2`) names a nucleotide that
/// is spliced out of the mature RNA, so it has no `r.` position of its own.
/// Dropping the offset re-anchors the description onto the exon-boundary base —
/// which *does* exist on the RNA — so the output stays well-formed while
/// describing a different nucleotide. That is only acceptable when the exonic
/// span is itself the thing being described, which is true for exactly one
/// shape:
///
/// **Whole exon(s).** A span reaching 5′ into the intron *preceding* its start
/// exon and 3′ into the intron *following* its end exon removes complete exons,
/// and the spec spells that consequence with plain exonic positions —
/// `RNA/splicing.md`: "`NC_000023.11(NM_004006.2):r.650_831del` … the sequence
/// from nucleotide `r.650` to `r.831` (exon 8) is deleted from the transcript".
/// HGVS intron numbering makes the test purely local: a **negative** offset is
/// anchored on the first base of the downstream exon and a **positive** offset
/// on the last base of the upstream exon (`background/numbering.md`), so
/// `start.offset < 0 && end.offset > 0` means `start.base ..= end.base` is
/// exactly a whole number of exons.
///
/// Everything else declines: a partial-exon edit with one intronic endpoint
/// (`c.3675-45_3692del` — #1086 D5) has no exon-sized RNA consequence to name,
/// and a purely intronic edit (`c.100+5_100+10del`, `c.831+2T>A`) touches no
/// exonic base at all. For those, `background/numbering.md` is explicit that an
/// RNA reference "can therefore **not be used** to describe variants affecting
/// these sequences"; the spec's remedy is a genome-anchored `r.` carrying the
/// offsets, which this surface does not emit.
///
/// Note this is *not* a general "make `r.` agree with `c.`" rule: the 3′-shift
/// difference between the two axes is spec-mandated (`RNA/deletion.md`) and is
/// deliberately left alone — only exonic spans reach the shift, exactly as
/// before.
fn exonic_span(iv: &CdsInterval) -> Option<CdsInterval> {
    let start_in = iv.start.inner()?;
    let end_in = iv.end.inner()?;
    if start_in.is_intronic() || end_in.is_intronic() {
        let spans_whole_exons =
            start_in.offset.is_some_and(|o| o < 0) && end_in.offset.is_some_and(|o| o > 0);
        if !spans_whole_exons {
            return None;
        }
    }
    // A `Some(0)` offset is exonic but would render as a spurious `+0`; the
    // whole-exon offsets above are deliberately dropped (that is the exon-span
    // reduction). Either way the output is exon-anchored.
    let strip = |p: &CdsPos| CdsPos {
        base: p.base,
        offset: None,
        utr3: p.utr3,
        special: p.special,
    };
    Some(CdsInterval::new(strip(start_in), strip(end_in)))
}

/// Resolve a CDS position range to literal DNA bases against `transcript.sequence`.
///
/// `start`/`end` are 1-based inclusive CDS coordinates; `invert` reverse-complements
/// the slice (for `…inv`). Mirrors `fetch_position_range_bases` (rules.rs): the
/// transcript index is `cds_start + pos - 1`. Returns DNA (uppercase) bases; the
/// RNA T→U + lowercase transform is applied later by Display.
fn resolve_range_bases(
    start: u64,
    end: u64,
    invert: bool,
    transcript: &Transcript,
) -> Option<Sequence> {
    let cds_start = transcript.cds_start?;
    let seq = transcript.sequence.as_ref()?;
    let tx_start = cds_start.checked_add(start.checked_sub(1)?)?;
    let tx_end = cds_start.checked_add(end.checked_sub(1)?)?;
    // `cds_start` is documented 1-based but not type-enforced; guard the `- 1`
    // against underflow (mirrors `fetch_position_range_bases`' cds_start==0 reject).
    let lo = (tx_start.checked_sub(1)?) as usize;
    let bases = seq.get(lo..(tx_end as usize))?.to_string();
    let bases = if invert {
        reverse_complement(&bases)
    } else {
        bases
    };
    bases.parse::<Sequence>().ok()
}

/// Resolve a position-range insert against the transcript sequence, in CDS
/// coordinates, to a literal `InsertedSequence` (DNA bases; Display applies T→U +
/// lowercase). Handles the bare `PositionRange`/`PositionRangeInv` variants and the
/// single-part `Complex([PositionRange | PositionRangeInv])` shape the parser
/// produces for `delins191_194inv` — both collapse to a single `Literal`. Inserts
/// that aren't a resolvable same-reference position range pass through unchanged.
fn resolve_range_insert(
    ins: &InsertedSequence,
    transcript: &Transcript,
) -> Option<InsertedSequence> {
    let (start, end, invert) = match ins {
        InsertedSequence::PositionRange { start, end } => (*start, *end, false),
        InsertedSequence::PositionRangeInv { start, end } => (*start, *end, true),
        InsertedSequence::Complex(parts) => match parts.as_slice() {
            [InsertedPart::PositionRange { start, end }] => (*start, *end, false),
            [InsertedPart::PositionRangeInv { start, end }] => (*start, *end, true),
            _ => return Some(ins.clone()),
        },
        other => return Some(other.clone()),
    };
    Some(InsertedSequence::Literal(resolve_range_bases(
        start, end, invert, transcript,
    )?))
}

/// Predict the RNA form from the normalized coding (`c.`) variant.
///
/// Returns `None` (axis "not available") for anything not representable as a
/// concrete exonic RNA edit. Never panics, never fabricates output.
pub(crate) fn predict_rna(coding: &HgvsVariant, transcript: &Transcript) -> Option<HgvsVariant> {
    // Cis/compound allele: predict each member, then strip its predicted `( )`
    // wrapper to a plain (Certain) RNA edit so the allele renders the spec-correct
    // `r.[(a;b)]` — one inside-bracket wrapper from the allele's `uncertain` flag —
    // not `r.[((a);(b))]`. (mutalyzer emits the non-spec outside-bracket `r.([a;b])`;
    // ferro intentionally diverges per recommendations/RNA/alleles.md.)
    if let HgvsVariant::Allele(allele) = coding {
        let mut rna_members = Vec::with_capacity(allele.variants.len());
        for member in &allele.variants {
            let predicted = predict_rna(member, transcript)?;
            let HgvsVariant::Rna(r) = predicted else {
                return None;
            };
            let certain = RnaVariant {
                accession: r.accession,
                gene_symbol: r.gene_symbol,
                loc_edit: LocEdit::new(r.loc_edit.location, r.loc_edit.edit.into_inner()?),
            };
            rna_members.push(HgvsVariant::Rna(certain));
        }
        return Some(HgvsVariant::Allele(AlleleVariant::new_uncertain(
            rna_members,
            allele.phase,
        )));
    }
    // Require a transcript sequence (later tasks read it for shifting and
    // range-reference resolution); absent → axis not available.
    transcript.sequence.as_ref()?;
    let cds = match coding {
        HgvsVariant::Cds(c) => c,
        _ => return None, // non-c. inputs unsupported (alleles handled above)
    };
    let edit: &NaEdit = cds.loc_edit.edit.inner()?;
    // Reduce to the exonic span the r. axis may name (declining an intronic
    // endpoint that has none), then 3′-shift deletions along the spliced sequence.
    let exonic = exonic_span(&cds.loc_edit.location)?;
    let location = match edit {
        NaEdit::Deletion { .. } => shift_deletion_3prime(&exonic, transcript)?,
        _ => exonic,
    };
    let rna_interval = cds_interval_to_rna(&location)?;
    // Resolve position-range/inversion inserts against the transcript sequence to
    // literal bases so Display renders RNA bases rather than position notation.
    let edit_out: NaEdit = match edit {
        NaEdit::Delins {
            sequence,
            deleted,
            deleted_length,
            substitution_reference,
        } => NaEdit::Delins {
            sequence: resolve_range_insert(sequence, transcript)?,
            deleted: deleted.clone(),
            deleted_length: *deleted_length,
            substitution_reference: substitution_reference.clone(),
        },
        NaEdit::Insertion { sequence } => NaEdit::Insertion {
            sequence: resolve_range_insert(sequence, transcript)?,
        },
        other => other.clone(),
    };
    let rna_variant = RnaVariant {
        accession: cds.accession.clone(),
        gene_symbol: cds.gene_symbol.clone(),
        loc_edit: LocEdit::new_predicted(rna_interval, edit_out),
    };
    Some(HgvsVariant::Rna(rna_variant))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;
    use crate::reference::transcript::{Exon, Strand, Transcript};

    /// A 300 nt plus-strand transcript, CDS = whole transcript (cds_start = 1),
    /// single exon. Base 169 is `T`, used by the substitution test. Bases 141/142
    /// are distinct from their neighbours so `c.141_142del` is already 3′-most and
    /// does not shuffle (the rest of the transcript is a homopolymeric `A` run).
    fn coding_tx() -> Transcript {
        let mut seq = vec![b'A'; 300];
        seq[168] = b'T'; // c.169 (1-based) -> index 168
        seq[139] = b'G'; // c.140 -> blocks a left/right run through the deletion
        seq[140] = b'G'; // c.141 (deleted) -> distinct from c.143 ('A') => no 3′ shift
        Transcript {
            id: "NM_000001.1".to_string(),
            strand: Strand::Plus,
            sequence: Some(String::from_utf8(seq).unwrap()),
            cds_start: Some(1),
            cds_end: Some(300),
            exons: vec![Exon::new(1, 1, 300)],
            ..Default::default()
        }
    }

    fn predict(input: &str, tx: &Transcript) -> String {
        let v = parse_hgvs(input).expect("parse");
        predict_rna(&v, tx)
            .expect("predict_rna returned None")
            .to_string()
    }

    #[test]
    fn non_coding_input_returns_none() {
        // A genomic (g.) input is not a c./n. coding variant — the surface is not
        // available for it (this is also what makes the projector's genome-path
        // wiring safe: `.rna` stays `None` for non-`Cds` inputs).
        let tx = coding_tx();
        let v = parse_hgvs("NC_000001.11:g.169A>T").expect("parse");
        assert!(predict_rna(&v, &tx).is_none());
    }

    #[test]
    fn missing_transcript_sequence_returns_none() {
        // Without a transcript sequence the 3′ shift / range resolution cannot run,
        // so the axis is "not available" rather than a fabricated output.
        let mut tx = coding_tx();
        tx.sequence = None;
        let v = parse_hgvs("NM_000001.1:c.169T>A").expect("parse");
        assert!(predict_rna(&v, &tx).is_none());
    }

    #[test]
    fn substitution_renders_rna_predicted() {
        let tx = coding_tx();
        assert_eq!(
            predict("NM_000001.1:c.169T>A", &tx),
            "NM_000001.1:r.(169u>a)"
        );
    }

    #[test]
    fn deletion_renders_rna() {
        let tx = coding_tx();
        assert_eq!(
            predict("NM_000001.1:c.141_142del", &tx),
            "NM_000001.1:r.(141_142del)"
        );
    }

    /// Plus-strand transcript with an `AAA` run at c.567..569 (indices 566..568)
    /// so `c.568del` 3′-normalizes to `c.569del`.
    fn homopolymer_tx() -> Transcript {
        let mut seq = vec![b'C'; 600];
        seq[566] = b'A'; // c.567
        seq[567] = b'A'; // c.568
        seq[568] = b'A'; // c.569
        Transcript {
            id: "NM_000002.1".to_string(),
            strand: Strand::Plus,
            sequence: Some(String::from_utf8(seq).unwrap()),
            cds_start: Some(1),
            cds_end: Some(600),
            exons: vec![Exon::new(1, 1, 600)],
            ..Default::default()
        }
    }

    #[test]
    fn deletion_three_prime_shifts_in_transcript_frame() {
        let tx = homopolymer_tx();
        assert_eq!(
            predict("NM_000002.1:c.568del", &tx),
            "NM_000002.1:r.(569del)"
        );
    }

    #[test]
    fn insertion_literal_renders_rna() {
        let tx = coding_tx();
        assert_eq!(
            predict("NM_000001.1:c.169_170insATA", &tx),
            "NM_000001.1:r.(169_170insaua)"
        );
    }

    #[test]
    fn delins_literal_renders_rna() {
        let tx = coding_tx();
        assert_eq!(
            predict("NM_000001.1:c.5_6delinsTAG", &tx),
            "NM_000001.1:r.(5_6delinsuag)"
        );
    }

    fn utr_tx() -> Transcript {
        Transcript {
            id: "NM_000003.1".to_string(),
            strand: Strand::Plus,
            sequence: Some("G".repeat(400)),
            cds_start: Some(1),
            cds_end: Some(297),
            exons: vec![Exon::new(1, 1, 400)],
            ..Default::default()
        }
    }

    #[test]
    fn utr3_position_preserved() {
        let tx = utr_tx();
        assert_eq!(
            predict("NM_000003.1:c.297_*1del", &tx),
            "NM_000003.1:r.(297_*1del)"
        );
    }

    #[test]
    fn whole_exon_span_reduces_to_its_exonic_positions() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, LocEdit};

        // `c.677-18_704+65del` brackets whole exon(s): the start offset is
        // negative (anchored on the exon's FIRST base) and the end offset is
        // positive (anchored on the exon's LAST base), so c.677..704 is exactly
        // the skipped exon. RNA/splicing.md spells that consequence with plain
        // exonic positions (`r.650_831del`, "exon 8 is deleted from the
        // transcript"), so the offsets are dropped here by design — this is the
        // one intronic shape #1086 keeps.
        let mut tx = coding_tx();
        let mut s = vec![b'A'; 900];
        s[704] = b'C'; // c.705 differs from the deleted `A` run end => no 3' shift
        tx.sequence = Some(String::from_utf8(s).unwrap());
        tx.cds_end = Some(900);
        tx.exons = vec![Exon::new(1, 1, 900)];

        let iv = CdsInterval::new(
            CdsPos::with_offset(677, -18), // c.677-18
            CdsPos::with_offset(704, 65),  // c.704+65
        );
        let cds = CdsVariant {
            accession: Accession::new("NM", "000001", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                iv,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        let out = predict_rna(&HgvsVariant::Cds(cds), &tx).expect("some");
        assert_eq!(out.to_string(), "NM_000001.1:r.(677_704del)");
    }

    /// Multi-exon whole-exon skip (#1111): a base range whose start offset is
    /// negative (anchored on the first base of one exon) and end offset is
    /// positive (anchored on the last base of a *later* exon) brackets a whole
    /// number of exons. `exonic_span` keys purely off the offset *signs* — the
    /// exon count is not a distinct code path — so a span across two
    /// consecutive exons reduces to plain exonic positions exactly like the
    /// single-exon case above. This pins that data point directly.
    #[test]
    fn multi_exon_whole_exon_span_reduces_to_its_exonic_positions() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, LocEdit};

        // Three contiguous exons in transcript space: exon 1 = tx 1..100,
        // exon 2 = tx 101..200, exon 3 = tx 201..300. CDS = whole transcript
        // (cds_start = 1), so c.pos == tx.pos. `c.101-5_300+5del` anchors its
        // start on the first base of exon 2 (negative offset) and its end on
        // the last base of exon 3 (positive offset) — bracketing exons 2 and 3
        // in full. The deletion ends at the transcript's last base, so there is
        // no room for a 3′ shift and the offsets drop to `r.(101_300del)`.
        let seq = vec![b'A'; 300];
        let tx = Transcript {
            id: "NM_000005.1".to_string(),
            strand: Strand::Plus,
            sequence: Some(String::from_utf8(seq).unwrap()),
            cds_start: Some(1),
            cds_end: Some(300),
            exons: vec![
                Exon::new(1, 1, 100),
                Exon::new(2, 101, 200),
                Exon::new(3, 201, 300),
            ],
            ..Default::default()
        };
        let iv = CdsInterval::new(
            CdsPos::with_offset(101, -5), // c.101-5 (first base of exon 2)
            CdsPos::with_offset(300, 5),  // c.300+5 (last base of exon 3)
        );
        let cds = CdsVariant {
            accession: Accession::new("NM", "000005", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                iv,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        let out = predict_rna(&HgvsVariant::Cds(cds), &tx).expect("some");
        assert_eq!(out.to_string(), "NM_000005.1:r.(101_300del)");
    }

    #[test]
    fn partial_exon_span_with_intronic_endpoints_declines() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, LocEdit};

        // Both endpoints intronic, but bracketed the WRONG way round:
        // `c.677+18_704-65del` starts in the intron *following* exon-base 677
        // and ends in the intron *preceding* exon-base 704, so it clips the
        // flanking exons rather than removing whole ones. There is no exon-sized
        // RNA consequence to name — decline (#1086).
        let tx = coding_tx();
        let iv = CdsInterval::new(CdsPos::with_offset(150, 18), CdsPos::with_offset(180, -65));
        let cds = CdsVariant {
            accession: Accession::new("NM", "000001", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                iv,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        assert!(predict_rna(&HgvsVariant::Cds(cds), &tx).is_none());
    }

    #[test]
    fn single_intronic_endpoint_declines() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, LocEdit};

        // #1086 D5, the reported shape: only the *start* is intronic
        // (`c.150-45_180del`), so the span clips into an exon rather than
        // removing whole ones. Stripping the offset would name r.150 — an exonic
        // base that exists, but is not `c.150-45`.
        let tx = coding_tx();
        let iv = CdsInterval::new(CdsPos::with_offset(150, -45), CdsPos::new(180));
        let cds = CdsVariant {
            accession: Accession::new("NM", "000001", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                iv,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        assert!(predict_rna(&HgvsVariant::Cds(cds), &tx).is_none());
    }

    #[test]
    fn zero_offset_endpoint_is_exonic_and_renders() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, LocEdit};

        // A `Some(0)` offset is exonic, not intronic: it must still render (and
        // without a spurious `+0`). Pins that the #1086 guard keys off
        // `is_intronic()` rather than `offset.is_some()`.
        let tx = coding_tx();
        let mut start = CdsPos::new(141);
        start.offset = Some(0);
        let iv = CdsInterval::new(start, CdsPos::new(142));
        let cds = CdsVariant {
            accession: Accession::new("NM", "000001", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                iv,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        let out = predict_rna(&HgvsVariant::Cds(cds), &tx).expect("some");
        assert_eq!(out.to_string(), "NM_000001.1:r.(141_142del)");
    }

    #[test]
    fn range_ref_inv_insert_resolves_to_bases() {
        // CDS = whole transcript; put a known 4-mer at c.191..194 (idx 190..193).
        let mut seq = vec![b'A'; 300];
        seq[190] = b'A';
        seq[191] = b'C';
        seq[192] = b'G';
        seq[193] = b'T'; // c.191..194 = "ACGT"; revcomp = "ACGT" -> rna "acgu"
        let tx = Transcript {
            id: "NM_000004.1".to_string(),
            strand: Strand::Plus,
            sequence: Some(String::from_utf8(seq).unwrap()),
            cds_start: Some(1),
            cds_end: Some(300),
            exons: vec![Exon::new(1, 1, 300)],
            ..Default::default()
        };
        assert_eq!(
            predict("NM_000004.1:c.206_209delins191_194inv", &tx),
            "NM_000004.1:r.(206_209delinsacgu)"
        );
    }

    /// A minus-strand single-exon transcript, CDS = whole transcript. `predict_rna`
    /// operates purely in spliced transcript-sequence space (it reads
    /// `transcript.sequence` directly and never projects to the genome), so a
    /// minus-strand transcript behaves identically to a plus-strand one for the same
    /// spliced sequence. These tests pin that invariant: strand orientation does not
    /// change the RNA prediction. Base 169 is `T` (for the substitution); bases
    /// 141/142 are distinct from their neighbours so `c.141_142del` does not shuffle.
    fn minus_strand_tx() -> Transcript {
        let mut seq = vec![b'A'; 300];
        seq[168] = b'T'; // c.169 (1-based) -> index 168
        seq[139] = b'G'; // c.140 -> blocks a run through the deletion
        seq[140] = b'G'; // c.141 (deleted) -> distinct from c.143 ('A') => no 3' shift
        Transcript {
            id: "NM_000006.1".to_string(),
            strand: Strand::Minus,
            sequence: Some(String::from_utf8(seq).unwrap()),
            cds_start: Some(1),
            cds_end: Some(300),
            exons: vec![Exon::new(1, 1, 300)],
            ..Default::default()
        }
    }

    #[test]
    fn minus_strand_substitution_renders_rna() {
        let tx = minus_strand_tx();
        assert_eq!(
            predict("NM_000006.1:c.169T>A", &tx),
            "NM_000006.1:r.(169u>a)"
        );
    }

    #[test]
    fn minus_strand_deletion_renders_rna() {
        let tx = minus_strand_tx();
        assert_eq!(
            predict("NM_000006.1:c.141_142del", &tx),
            "NM_000006.1:r.(141_142del)"
        );
    }

    /// A two-exon plus-strand transcript whose exons abut in transcript space
    /// (exon 1 = tx 1..150, exon 2 = tx 151..300), CDS = whole transcript. Bases
    /// c.149..152 (indices 148..151) are the distinct 4-mer `CGTA` flanked by the
    /// homopolymeric `A` run, so a junction-spanning `c.149_152del` is already
    /// 3'-most and the coordinate projection across the exon boundary is exercised
    /// without a confounding shuffle.
    fn multi_exon_tx() -> Transcript {
        let mut seq = vec![b'A'; 300];
        seq[148] = b'C'; // c.149 (exon 1, last two exonic bases)
        seq[149] = b'G'; // c.150
        seq[150] = b'T'; // c.151 (exon 2, first two exonic bases)
        seq[151] = b'A'; // c.152
        Transcript {
            id: "NM_000007.1".to_string(),
            strand: Strand::Plus,
            sequence: Some(String::from_utf8(seq).unwrap()),
            cds_start: Some(1),
            cds_end: Some(300),
            exons: vec![Exon::new(1, 1, 150), Exon::new(2, 151, 300)],
            ..Default::default()
        }
    }

    #[test]
    fn multi_exon_deletion_spans_exon_junction() {
        let tx = multi_exon_tx();
        // The deletion straddles the exon 1 / exon 2 boundary (c.150 | c.151);
        // coordinate projection through the mapper must keep it intact on the r. axis.
        assert_eq!(
            predict("NM_000007.1:c.149_152del", &tx),
            "NM_000007.1:r.(149_152del)"
        );
    }

    #[test]
    fn purely_intronic_edit_returns_none() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, LocEdit};

        // An intron-only deletion whose endpoints share an exon anchor —
        // c.100+5_100+10del — has no RNA consequence. Stripping its offsets would
        // otherwise collapse it to a spurious single-base c.100_100del; predict_rna
        // must decline (return None) instead. Covered by `exonic_span` since #1086
        // (offset `+5` is not the whole-exon `-`/`+` bracketing shape).
        let tx = coding_tx();
        let start = CdsPos::with_offset(100, 5); // c.100+5
        let end = CdsPos::with_offset(100, 10); // c.100+10
        let iv = CdsInterval::new(start, end);
        let cds = CdsVariant {
            accession: Accession::new("NM", "000001", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                iv,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        assert!(predict_rna(&HgvsVariant::Cds(cds), &tx).is_none());
    }

    #[test]
    fn purely_intronic_edit_across_adjacent_anchors_returns_none() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, LocEdit};

        // A deletion lying wholly within a single intron whose endpoints anchor to
        // *adjacent* exonic bases — c.100+5_101-3del. The span never touches an
        // exonic base, so it has no RNA consequence. Stripping the offsets would
        // otherwise fabricate a spurious two-base exonic c.100_101del; predict_rna
        // must decline (return None) instead. Covered by `exonic_span` since #1086.
        let tx = coding_tx();
        let start = CdsPos::with_offset(100, 5); // c.100+5
        let end = CdsPos::with_offset(101, -3); // c.101-3
        let iv = CdsInterval::new(start, end);
        let cds = CdsVariant {
            accession: Accession::new("NM", "000001", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                iv,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ),
        };
        assert!(predict_rna(&HgvsVariant::Cds(cds), &tx).is_none());
    }

    #[test]
    fn cis_allele_emits_spec_correct_inside_bracket_form() {
        // SPEC NOTE: HGVS recommends the predicted `( )` INSIDE the brackets —
        // `r.[(a;b)]` — per assets/hgvs-nomenclature/docs/recommendations/RNA/
        // alleles.md:32. mutalyzer emits the outside-bracket `r.([a;b])`, which
        // is non-spec and which ferro intentionally does NOT produce. ferro
        // therefore diverges from mutalyzer here by being more spec-correct; the
        // conformance row is annotated as an `improvement` (Task 9).
        let mut seq = vec![b'A'; 300];
        seq[273] = b'G'; // c.274
        seq[277] = b'A'; // c.278
        let tx = Transcript {
            id: "NM_000005.1".to_string(),
            strand: Strand::Plus,
            sequence: Some(String::from_utf8(seq).unwrap()),
            cds_start: Some(1),
            cds_end: Some(300),
            exons: vec![Exon::new(1, 1, 300)],
            ..Default::default()
        };
        assert_eq!(
            predict("NM_000005.1:c.[274G>T;278A>G]", &tx),
            "NM_000005.1:r.[(274g>u;278a>g)]"
        );
    }
}
