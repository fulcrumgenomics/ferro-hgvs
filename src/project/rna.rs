//! Predicted RNA consequence surface: `c.`/`n.` → `r.(…)`.
//!
//! Numbering is CDS-relative (identical to `c.`); the 3′ shift, intron-clamp,
//! and range-reference resolution are performed in transcript-sequence space.
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
/// a safety fallback (intronic offsets are already stripped by `clamp_intronic`
/// before this runs) — such intervals are returned unchanged.
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

/// Strip non-zero intronic offsets from a CDS interval, keeping the exon-anchor
/// `base`. Collapses an intron-spanning edit to its exonic span on the RNA
/// (introns are not present in the spliced transcript). Must run before any
/// 3′ shift, which is undefined on intronic-offset positions.
fn clamp_intronic(iv: &CdsInterval) -> Option<CdsInterval> {
    // Collapsing to the exonic span: the output is always exon-anchored, so the
    // intronic offset is dropped unconditionally (an exonic input already has
    // `None`; a `Some(0)` boundary offset would otherwise render as `+0`).
    let strip = |p: &CdsPos| CdsPos {
        base: p.base,
        offset: None,
        utr3: p.utr3,
        special: p.special,
    };
    let start_in = iv.start.inner()?;
    let end_in = iv.end.inner()?;
    let start = strip(start_in);
    let end = strip(end_in);
    // A purely intronic edit (both endpoints carry a non-zero intronic offset)
    // that lies wholly within a single intron has no exonic component and hence no
    // RNA consequence — decline rather than fabricate a spurious exonic span. This
    // covers two shapes:
    //   * shared exon anchor — `c.100+5_100+10del` collapses to `c.100_100del`;
    //   * adjacent exon anchors — `c.100+5_101-3del` (the span sits in the lone
    //     intron between bases 100 and 101, which are consecutive exonic positions)
    //     collapses to `c.100_101del`.
    // Both endpoints intronic with `end.base - start.base <= 1` means no exonic base
    // falls inside the span.
    if end.base.saturating_sub(start.base) <= 1 && start_in.is_intronic() && end_in.is_intronic() {
        return None;
    }
    Some(CdsInterval::new(start, end))
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
    // Clamp intronic offsets to the exonic span first, then 3′-shift deletions
    // along the spliced sequence (shifting an intronic-offset edit is undefined).
    let clamped = clamp_intronic(&cds.loc_edit.location)?;
    let location = match edit {
        NaEdit::Deletion { .. } => shift_deletion_3prime(&clamped, transcript)?,
        _ => clamped,
    };
    let rna_interval = cds_interval_to_rna(&location)?;
    // Resolve position-range/inversion inserts against the transcript sequence to
    // literal bases so Display renders RNA bases rather than position notation.
    let edit_out: NaEdit = match edit {
        NaEdit::Delins {
            sequence,
            deleted,
            deleted_length,
        } => NaEdit::Delins {
            sequence: resolve_range_insert(sequence, transcript)?,
            deleted: deleted.clone(),
            deleted_length: *deleted_length,
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
    fn intronic_endpoints_clamp_to_exonic_span() {
        use crate::hgvs::edit::NaEdit;
        use crate::hgvs::interval::CdsInterval;
        use crate::hgvs::location::CdsPos;
        use crate::hgvs::variant::{Accession, CdsVariant, LocEdit};

        // sequence present, cds_start = 1; needs >=704 nt in-bounds. Make
        // position 705 distinct from the deleted run so the post-clamp 3′ shift
        // leaves `677_704del` in place rather than shuffling rightward.
        let mut tx = coding_tx();
        let mut s = vec![b'A'; 900];
        s[704] = b'C'; // 0-based 704 == c.705 differs from the deleted `A` run end
        tx.sequence = Some(String::from_utf8(s).unwrap());
        tx.cds_end = Some(900);
        tx.exons = vec![Exon::new(1, 1, 900)];

        let start = CdsPos::with_offset(677, -18); // c.677-18
        let end = CdsPos::with_offset(704, 65); // c.704+65
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
        let out = predict_rna(&HgvsVariant::Cds(cds), &tx).expect("some");
        assert_eq!(out.to_string(), "NM_000001.1:r.(677_704del)");
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
        // c.100+5_100+10del — has no RNA consequence. Clamping its offsets would
        // otherwise collapse it to a spurious single-base c.100_100del; predict_rna
        // must decline (return None) instead.
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
        // must decline (return None) instead.
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
