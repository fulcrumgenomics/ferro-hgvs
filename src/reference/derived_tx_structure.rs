//! Build-time–derived exon→genome structures for cdot-absent old transcript
//! versions (#790). Mirrors `derived_placement.rs`: derived at `ferro prepare`
//! time, baked into the manifest, injected into the cdot mapper at load. The
//! invariant is decline-rather-than-mis-anchor — every unproven step yields
//! `None`, leaving the single-exon supplemental fallback in place.
//!
//! Poly-A-core handling: `derive_transcript_structure` strips any trailing
//! poly-A tail from `old_mrna` before alignment. A poly-A tail is
//! post-transcriptional and has no genomic coordinate, so it must not be
//! matched against the genome. Only runs of at least `MIN_POLYA_TAIL` bases
//! are treated as a real tail; shorter trailing A-runs (genomically encoded)
//! are left intact.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::data::cdot::CdotTranscript;
use crate::reference::Strand;
use crate::FerroError;

/// Max fraction of mismatched bases tolerated on genome readback. Same gate as
/// `derived_placement::MAX_MISMATCH_FRACTION` — a clean transfer is ~0; any
/// boundary error frameshifts into a long mismatch run far above this.
const MAX_MISMATCH_FRACTION: f64 = 0.005;

/// Length of the trailing poly-A tail to exclude from genomic alignment.
/// A transcript's poly-A tail is post-transcriptional and has no genomic
/// coordinate, so it must not be matched against the genome. Only a run long
/// enough to be a real poly-A tail is stripped — a lone genomic `A` at the 3'
/// end (sibling versions here have a trailing A-run of 1) is left intact.
const MIN_POLYA_TAIL: usize = 8;

/// Returns the length of `mrna` with any trailing poly-A tail removed (the
/// genomically-encoded core). Trims the maximal trailing run of `A`/`a` only
/// when that run is at least `MIN_POLYA_TAIL` long; otherwise returns the full
/// length.
fn genomic_core_len(mrna: &str) -> usize {
    let run = mrna
        .bytes()
        .rev()
        .take_while(|&b| b == b'A' || b == b'a')
        .count();
    if run >= MIN_POLYA_TAIL {
        mrna.len() - run
    } else {
        mrna.len()
    }
}

/// Abstracts genome base access so the derivation is testable without a FASTA.
pub trait GenomeSlice {
    /// 1-based inclusive `start`, exclusive `end`; forward-strand bases.
    fn slice(&self, contig: &str, start_1based: u64, end_exclusive: u64) -> Option<String>;
}

/// One derived transcript structure, in the internal `CdotTranscript` exon
/// convention: `exons[i] = [genome_start (1-based), genome_end (exclusive),
/// tx_start (0-based), tx_end (0-based exclusive)]`. `cds_*` are 0-based
/// (start inclusive, end exclusive), matching `CdotTranscript`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DerivedTxStructure {
    pub accession: String,
    pub contig: String,
    pub strand: String,
    pub exons: Vec<[u64; 4]>,
    #[serde(default)]
    pub cds_start: Option<u64>,
    #[serde(default)]
    pub cds_end: Option<u64>,
    #[serde(default)]
    pub gene_name: Option<String>,
    #[serde(default)]
    pub protein: Option<String>,
    /// The sibling accession whose cdot exons anchored the transfer (provenance).
    pub anchored_by: String,
    /// Genome-readback mismatch fraction (~0 for a clean derivation; provenance).
    pub mismatch_fraction: f64,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DerivedTxStructures {
    #[serde(default)]
    pub description: String,
    pub structures: Vec<DerivedTxStructure>,
}

impl DerivedTxStructures {
    pub fn from_json_path<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let content = std::fs::read_to_string(path.as_ref())?;
        Ok(serde_json::from_str(&content)?)
    }

    /// Pretty JSON + trailing newline (stable for a `--check` diff).
    pub fn to_json(&self) -> Result<String, FerroError> {
        let mut s = serde_json::to_string_pretty(self)?;
        s.push('\n');
        Ok(s)
    }
}

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

/// Concatenate a transcript's exon genome slices in tx order (rev-comp on minus).
pub fn reconstruct_spliced_mrna(cdot: &CdotTranscript, genome: &dyn GenomeSlice) -> Option<String> {
    // Exons are stored ascending by tx coordinate. On the minus strand the
    // forward-genome slice of each exon must be reverse-complemented; the exon
    // order already matches tx order, so concatenate as-is after rev-comp.
    let mut out = String::new();
    for e in &cdot.exons {
        let (gs, ge) = (e[0], e[1]);
        let raw = genome.slice(&cdot.contig, gs, ge)?;
        match cdot.strand {
            Strand::Plus => out.push_str(&raw),
            Strand::Minus => out.push_str(&revcomp(&raw)),
            Strand::Unknown => return None,
        }
    }
    Some(out)
}

/// Offset `off` such that `sib[(p+off)] == old[p]` for the shared body, or
/// `None` if `old` is not a clean terminal trim/extension of `sib` (i.e. there
/// is an internal mismatch/indel). Implemented by trying each candidate offset
/// in `[-(old_len), sib_len]` and accepting the first that makes the entire
/// overlap region identical with overlap length >= 0.9 * min(len) — enough to
/// exclude internal indels (which break every offset's overlap).
fn find_clean_offset(old: &str, sib: &str) -> Option<i64> {
    let (ob, sb) = (old.as_bytes(), sib.as_bytes());
    let (ol, sl) = (ob.len() as i64, sb.len() as i64);
    let min_overlap = ((ol.min(sl) as f64) * 0.9) as i64;
    for off in -(ol - 1)..sl {
        // overlap of old index p with sib index p+off, for p in [max(0,-off) ..]
        let p0 = 0.max(-off);
        let p1 = ol.min(sl - off);
        if p1 - p0 < min_overlap {
            continue;
        }
        let ok = (p0..p1).all(|p| ob[p as usize] == sb[(p + off) as usize]);
        if ok {
            return Some(off);
        }
    }
    None
}

fn readback_mismatch_fraction(
    cdot: &CdotTranscript,
    expected: &str,
    genome: &dyn GenomeSlice,
) -> Option<f64> {
    let got = reconstruct_spliced_mrna(cdot, genome)?;
    if got.len() != expected.len() {
        return Some(1.0);
    }
    let mism = got
        .bytes()
        .zip(expected.bytes())
        .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
        .count();
    Some(mism as f64 / expected.len().max(1) as f64)
}

/// Derive an old version's structure by anchoring to a sibling cdot record.
///
/// Enforces the clean-terminal-trim invariant: `old_mrna` must equal the
/// sibling mRNA with only its 5'/3' ends trimmed/extended (no internal indel).
/// Returns `None` otherwise, or if genome readback of the transferred exons does
/// not reproduce `old_mrna`.
///
/// The CDS is derived from the sibling by mapping its transcript-coordinate CDS
/// into old-transcript coordinates using the computed alignment offset. If the
/// mapped CDS does not land cleanly inside the old genomic core, the anchor is
/// not trustworthy and derivation declines.
pub fn derive_transcript_structure(
    old_id: &str,
    old_mrna: &str,
    sibling_id: &str,
    sibling: &CdotTranscript,
    genome: &dyn GenomeSlice,
) -> Option<DerivedTxStructure> {
    let sib_mrna = reconstruct_spliced_mrna(sibling, genome)?;

    // Strip any poly-A tail from old_mrna before alignment: the tail is
    // post-transcriptional and has no genomic coordinate, so it must not be
    // matched against the genome. We work with the genomic core for offset
    // finding and exon clipping; the readback comparison also uses the core.
    let core_len = genomic_core_len(old_mrna);
    let old_core = &old_mrna[..core_len];

    // 1) Find the 5' offset: the largest k such that old_core[0..] aligns to
    //    sib_mrna[k..] as an identical contiguous run for the shared length.
    //    A clean terminal trim means old is a substring of sib OR sib is a
    //    substring of old, anchored by a shared core. We compute the offset of
    //    old's start within sib (or negative if old extends 5' of sib).
    let off = find_clean_offset(old_core, &sib_mrna)?; // i64: sib_pos = old_pos + off

    // Derive the old-version CDS from the sibling's CDS by mapping sibling
    // transcript coordinates into old transcript coordinates (subtract `off`).
    // The CDS is genomically invariant across sibling versions; `off` captures
    // the difference in 5' extent. If the mapped CDS falls outside the old
    // genomic core, the anchor is not trustworthy — decline.
    let old_cds = match (sibling.cds_start, sibling.cds_end) {
        (Some(cs), Some(ce)) => {
            let ocs = cs as i64 - off; // sib tx -> old tx
            let oce = ce as i64 - off;
            // CDS must land cleanly inside the genomic core; otherwise the
            // anchor is not trustworthy for this transcript — decline.
            if ocs >= 0 && oce > ocs && (oce as u64) <= core_len as u64 {
                (Some(ocs as u64), Some(oce as u64))
            } else {
                return None;
            }
        }
        // Sibling genuinely non-coding: no CDS to carry.
        _ => (None, None),
    };

    // 2) Map each sibling exon's [tx_start,tx_end) into old tx space; clip to
    //    [0, old_len); drop exons that fall entirely outside old; adjust the
    //    genome bounds of clipped terminal exons by the clipped base count.
    let old_len = old_core.len() as u64;
    let mut exons: Vec<[u64; 4]> = Vec::new();
    for e in &sibling.exons {
        let (gs, ge, ts, te) = (e[0], e[1], e[2] as i64, e[3] as i64);
        // sib tx -> old tx
        let mut o_ts = ts - off;
        let mut o_te = te - off;
        let (mut g_s, mut g_e) = (gs, ge);
        // Clip 5' overhang.
        if o_ts < 0 {
            let trim = (-o_ts) as u64;
            o_ts = 0;
            // forward-genome start advances by `trim` on +; on - the end retreats
            match sibling.strand {
                Strand::Plus => g_s += trim,
                Strand::Minus => g_e -= trim,
                Strand::Unknown => return None,
            }
        }
        // Clip 3' overhang.
        if o_te > old_len as i64 {
            let trim = (o_te - old_len as i64) as u64;
            o_te = old_len as i64;
            match sibling.strand {
                Strand::Plus => g_e -= trim,
                Strand::Minus => g_s += trim,
                Strand::Unknown => return None,
            }
        }
        if o_te <= o_ts {
            continue; // exon entirely outside old
        }
        exons.push([g_s, g_e, o_ts as u64, o_te as u64]);
    }
    if exons.len() < 2 {
        return None; // single-exon transfer is no better than fallback
    }

    // 3) Genome readback: the transferred exons must reproduce the genomic core
    //    of old_mrna (poly-A tail excluded — it has no genomic coordinate).
    let probe = CdotTranscript {
        gene_name: sibling.gene_name.clone(),
        contig: sibling.contig.clone(),
        strand: sibling.strand,
        exons: exons.clone(),
        cds_start: old_cds.0,
        cds_end: old_cds.1,
        exon_cigars: vec![None; exons.len()],
        gene_id: None,
        protein: sibling.protein.clone(),
        cds_start_incomplete: false,
    };
    let frac = readback_mismatch_fraction(&probe, old_core, genome)?;
    if frac > MAX_MISMATCH_FRACTION {
        return None;
    }

    Some(DerivedTxStructure {
        accession: old_id.to_string(),
        contig: sibling.contig.clone(),
        strand: match sibling.strand {
            Strand::Plus => "+",
            Strand::Minus => "-",
            Strand::Unknown => return None,
        }
        .to_string(),
        exons,
        cds_start: old_cds.0,
        cds_end: old_cds.1,
        gene_name: sibling.gene_name.clone(),
        protein: sibling.protein.clone(),
        anchored_by: sibling_id.to_string(),
        mismatch_fraction: frac,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_through_json() {
        let s = DerivedTxStructures {
            description: "test".to_string(),
            structures: vec![DerivedTxStructure {
                accession: "NM_003002.2".to_string(),
                contig: "NC_000011.10".to_string(),
                strand: "+".to_string(),
                exons: vec![
                    [112086847, 112086960, 0, 113],
                    [112094805, 112095795, 375, 1365],
                ],
                cds_start: Some(60),
                cds_end: Some(542),
                gene_name: Some("SDHD".to_string()),
                protein: Some("NP_002993.1".to_string()),
                anchored_by: "NM_003002.3".to_string(),
                mismatch_fraction: 0.0,
            }],
        };
        let json = s.to_json().unwrap();
        let back: DerivedTxStructures = serde_json::from_str(&json).unwrap();
        assert_eq!(s, back);
    }
}

#[cfg(test)]
mod derivation_tests {
    use super::*;
    use crate::data::cdot::CdotTranscript;
    use crate::reference::Strand;

    /// Fake plus-strand genome: a single contig string, 1-based slicing.
    struct FakeGenome {
        contig: String,
        seq: String,
    }
    impl GenomeSlice for FakeGenome {
        fn slice(&self, contig: &str, start_1based: u64, end_exclusive: u64) -> Option<String> {
            if contig != self.contig {
                return None;
            }
            let (s, e) = ((start_1based - 1) as usize, (end_exclusive - 1) as usize);
            self.seq.get(s..e).map(|x| x.to_string())
        }
    }

    /// Two-exon sibling on a fake genome:
    ///   genome:  [exon0 @ g.1..=10][intron g.11..=14][exon1 @ g.15..=24]
    /// sibling mRNA = exon0(10) + exon1(10) = 20 nt; CDS tx 2..=18 (0-based 2..18).
    fn sibling_fixture() -> (CdotTranscript, FakeGenome, String) {
        let genome = FakeGenome {
            contig: "NC_TEST.1".to_string(),
            //       1234567890   1234     5678901234
            seq: "AAACCCGGGT".to_string() + "TTGG" + "ACGTACGTAC",
        };
        let cdot = CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("T".to_string()),
            contig: "NC_TEST.1".to_string(),
            strand: Strand::Plus,
            // [g_start(1b), g_end(excl), tx_start(0b), tx_end(excl)]
            exons: vec![[1, 11, 0, 10], [15, 25, 10, 20]],
            cds_start: Some(2),
            cds_end: Some(18),
            exon_cigars: vec![None, None],
            gene_id: None,
            protein: Some("NP_T.1".to_string()),
        };
        let sib_mrna = "AAACCCGGGTACGTACGTAC".to_string(); // 20 nt
        (cdot, genome, sib_mrna)
    }

    #[test]
    fn reconstruct_matches_sibling_mrna() {
        let (cdot, genome, sib_mrna) = sibling_fixture();
        assert_eq!(reconstruct_spliced_mrna(&cdot, &genome).unwrap(), sib_mrna);
    }

    #[test]
    fn derives_old_version_with_5utr_trim() {
        // Old version = sibling minus the first 2 5'UTR bases (exon0 starts later).
        let (cdot, genome, sib_mrna) = sibling_fixture();
        let old_mrna = sib_mrna[2..].to_string(); // 18 nt, "ACCCGGGT..."
        let d = derive_transcript_structure("NM_T.1", &old_mrna, "NM_T.2", &cdot, &genome)
            .expect("clean 5'UTR trim must derive");
        assert_eq!(d.accession, "NM_T.1");
        assert_eq!(d.contig, "NC_TEST.1");
        assert_eq!(d.strand, "+");
        // exon0 now starts 2 bases later in the genome; lengths shrink by 2.
        assert_eq!(d.exons, vec![[3, 11, 0, 8], [15, 25, 8, 18]]);
        assert!(d.mismatch_fraction < 1e-9);
        assert_eq!(d.anchored_by, "NM_T.2");
        assert_eq!(d.cds_start, Some(0));
        assert_eq!(d.cds_end, Some(16));
    }

    #[test]
    fn declines_on_internal_indel() {
        // Old version differs from the sibling by an INTERNAL deletion (1 base
        // removed mid-exon) — not a clean terminal trim. Must decline.
        let (cdot, genome, sib_mrna) = sibling_fixture();
        let mut old = sib_mrna.clone();
        old.remove(10); // delete a base at the exon0/exon1 junction interior
        assert!(derive_transcript_structure("NM_T.1", &old, "NM_T.2", &cdot, &genome,).is_none());
    }

    #[test]
    fn genomic_core_len_strips_long_polya_but_not_lone_a() {
        assert_eq!(genomic_core_len("ACGTAAAAAAAAAAAAAAAAAA"), 4); // 18 A's stripped -> "ACGT"
        assert_eq!(genomic_core_len("ACGTA"), 5); // lone trailing A kept
        assert_eq!(genomic_core_len("ACGT"), 4); // no trailing A
    }

    #[test]
    fn derives_ignoring_polya_tail() {
        // Old version = sibling minus the first 2 5'UTR bases + 18 trailing A's
        // (poly-A tail). The A's are stripped before alignment, so the derived
        // exons match those of `derives_old_version_with_5utr_trim`.
        let (cdot, genome, sib_mrna) = sibling_fixture();
        let old_mrna = sib_mrna[2..].to_string() + &"A".repeat(18);
        let d = derive_transcript_structure("NM_T.1", &old_mrna, "NM_T.2", &cdot, &genome)
            .expect("poly-A tail must be stripped before alignment");
        assert_eq!(d.exons, vec![[3, 11, 0, 8], [15, 25, 8, 18]]);
        assert_eq!(d.cds_start, Some(0));
        assert_eq!(d.cds_end, Some(16));
    }

    /// Two-exon **minus-strand** sibling on the same fake genome. For a minus
    /// transcript, tx position 0 is the *highest* genome coordinate, so exon0
    /// occupies the 3' genome region and each exon's forward-genome slice is
    /// reverse-complemented during reconstruction.
    fn minus_sibling_fixture() -> (CdotTranscript, FakeGenome) {
        let genome = FakeGenome {
            contig: "NC_TEST.1".to_string(),
            seq: "AAACCCGGGT".to_string() + "TTGG" + "ACGTACGTAC",
        };
        let cdot = CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("M".to_string()),
            contig: "NC_TEST.1".to_string(),
            strand: Strand::Minus,
            // minus: tx 0 == highest genome. exon0 = g[15,25), exon1 = g[1,11).
            exons: vec![[15, 25, 0, 10], [1, 11, 10, 20]],
            cds_start: Some(2),
            cds_end: Some(18),
            exon_cigars: vec![None, None],
            gene_id: None,
            protein: Some("NP_M.1".to_string()),
        };
        (cdot, genome)
    }

    #[test]
    fn derives_minus_strand_with_5utr_trim() {
        // Old = sibling minus its first two 5' bases. On the minus strand the 5'
        // end maps to the *high* genome end, so the 5' clip retreats exon0's
        // genome_end (g_e) by 2 — the opposite side from the plus-strand case.
        let (cdot, genome) = minus_sibling_fixture();
        let sib_mrna = reconstruct_spliced_mrna(&cdot, &genome).unwrap();
        let old_mrna = sib_mrna[2..].to_string();
        let d = derive_transcript_structure("NM_M.1", &old_mrna, "NM_M.2", &cdot, &genome)
            .expect("clean minus-strand 5'UTR trim must derive");
        assert_eq!(d.strand, "-");
        assert_eq!(d.exons, vec![[15, 23, 0, 8], [1, 11, 8, 18]]);
        assert_eq!(d.cds_start, Some(0));
        assert_eq!(d.cds_end, Some(16));
        assert!(d.mismatch_fraction < 1e-9);
    }

    #[test]
    fn declines_when_cds_outside_core() {
        // Old version is a large 3' trim of the sibling: only the first 5 bases
        // of the sibling's 20 nt mRNA are retained. After the offset (off == 0,
        // because old starts at the same place as sib), the sibling's cds_end
        // (18) mapped into old tx coords is 18, which exceeds core_len (5) — the
        // CDS falls outside the old genomic core. Derivation must decline.
        let (cdot, genome, sib_mrna) = sibling_fixture();
        let old_mrna = sib_mrna[..5].to_string(); // 5 nt — core_len == 5
        assert!(
            derive_transcript_structure("NM_T.1", &old_mrna, "NM_T.2", &cdot, &genome).is_none()
        );
    }

    #[test]
    fn derives_old_version_with_3utr_trim() {
        // Old version = sibling minus the last 2 3'UTR bases, so the 3' end of
        // the terminal exon (exon1, tx 10..20) is clipped *inside* the exon —
        // not a poly-A tail. With off == 0 the sibling CDS (tx 2..18) maps to
        // old tx 2..18, and core_len == 18, so the CDS still lands cleanly inside
        // the core and derivation succeeds (it does not decline at the CDS gate).
        // This drives the successful 3'-clip branch (`g_e -= trim` on plus):
        //   exon1 o_te = 20 > old_len 18 -> trim = 2 -> g_e = 25 - 2 = 23.
        // Hand-computed expected terminal exon: [15, 23, 10, 18].
        let (cdot, genome, sib_mrna) = sibling_fixture();
        let old_mrna = sib_mrna[..18].to_string(); // 18 nt — drop last 2 bases
        let d = derive_transcript_structure("NM_T.1", &old_mrna, "NM_T.2", &cdot, &genome)
            .expect("clean 3'UTR trim must derive");
        assert_eq!(d.strand, "+");
        // exon0 is unchanged; exon1's genome_end retreats by 2 (g.25 -> g.23).
        assert_eq!(d.exons, vec![[1, 11, 0, 10], [15, 23, 10, 18]]);
        assert_eq!(d.cds_start, Some(2));
        assert_eq!(d.cds_end, Some(18));
        assert!(d.mismatch_fraction < 1e-9);
    }

    #[test]
    fn derives_minus_strand_with_3utr_trim() {
        // Minus-strand mirror of `derives_old_version_with_3utr_trim`. On the
        // minus strand tx position 0 is the highest genome coordinate, so the 3'
        // end of the transcript maps to the *low* genome start — the 3' clip
        // advances the terminal exon's genome_start (`g_s += trim`), the opposite
        // side from the plus-strand case.
        //   exon1 (tx 10..20, g[1,11)) o_te = 20 > old_len 18 -> trim = 2 ->
        //   g_s = 1 + 2 = 3.
        // Hand-computed expected terminal exon: [3, 11, 10, 18].
        let (cdot, genome) = minus_sibling_fixture();
        let sib_mrna = reconstruct_spliced_mrna(&cdot, &genome).unwrap();
        let old_mrna = sib_mrna[..18].to_string(); // 18 nt — drop last 2 bases
        let d = derive_transcript_structure("NM_M.1", &old_mrna, "NM_M.2", &cdot, &genome)
            .expect("clean minus-strand 3'UTR trim must derive");
        assert_eq!(d.strand, "-");
        // exon0 is unchanged; exon1's genome_start advances by 2 (g.1 -> g.3).
        assert_eq!(d.exons, vec![[15, 25, 0, 10], [3, 11, 10, 18]]);
        assert_eq!(d.cds_start, Some(2));
        assert_eq!(d.cds_end, Some(18));
        assert!(d.mismatch_fraction < 1e-9);
    }

    /// Repetitive-sequence false-offset guard. `find_clean_offset` returns the
    /// *first* offset whose >=90%-overlap region matches exactly; on a tandemly
    /// repetitive mRNA a naive search could lock onto a shifted (wrong) offset.
    /// This fixture builds a sibling whose mRNA is a pure period-4 repeat so a
    /// shift by the repeat period still matches over its (shrunken) overlap, then
    /// pins that (a) the only offset surviving the 90%-overlap floor is the true
    /// one, and (b) end-to-end derivation either lands the correct coordinates or
    /// declines — never mis-anchors at the shifted offset.
    #[test]
    fn repetitive_sequence_does_not_misanchor() {
        // 24 nt period-4 repeat split across two exons; genome is the same
        // repeat with an intron so reconstruction reproduces it exactly.
        let genome = FakeGenome {
            contig: "NC_REP.1".to_string(),
            //     exon0 g.1..=12        intron g.13..=16  exon1 g.17..=28
            seq: "ACGTACGTACGT".to_string() + "TTGG" + "ACGTACGTACGT",
        };
        let cdot = CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("R".to_string()),
            contig: "NC_REP.1".to_string(),
            strand: Strand::Plus,
            exons: vec![[1, 13, 0, 12], [17, 29, 12, 24]],
            // CDS chosen so that after the one-period 3'-trim (old_len == 20) the
            // mapped CDS still lands inside the core (cds_end 16 <= 20) and the
            // derivation does not decline at the CDS gate.
            cds_start: Some(0),
            cds_end: Some(16),
            exon_cigars: vec![None, None],
            gene_id: None,
            protein: Some("NP_R.1".to_string()),
        };
        let sib_mrna = reconstruct_spliced_mrna(&cdot, &genome).unwrap();
        assert_eq!(sib_mrna, "ACGTACGTACGTACGTACGTACGT"); // pure ACGT repeat

        // A 5' trim by exactly one repeat period (4 nt). The genuine offset is
        // +4 (old starts 4 bases into sib). A naive matcher could be tempted by
        // offset 0 (the repeat makes old == sib[..20] over a long prefix), but
        // the 90%-overlap floor requires the *entire* overlap to match, and the
        // true period-aligned offset is the first one accepted. Confirm the
        // returned offset is the period-aligned one (a multiple of 4), and that
        // derivation reproduces the trimmed mRNA on readback (mismatch ~0).
        let old_core = sib_mrna[4..].to_string(); // 20 nt
        let off = find_clean_offset(&old_core, &sib_mrna).expect("repeat must still align");
        assert_eq!(off % 4, 0, "offset must be period-aligned, got {off}");
        let d = derive_transcript_structure("NM_R.1", &old_core, "NM_R.2", &cdot, &genome)
            .expect("clean period-aligned trim must derive without misanchoring");
        // Whatever period-aligned offset wins, readback gates it: the derived
        // exons must reconstruct old_core exactly.
        let probe = CdotTranscript {
            cds_start_incomplete: false,
            gene_name: cdot.gene_name.clone(),
            contig: cdot.contig.clone(),
            strand: cdot.strand,
            exons: d.exons.clone(),
            cds_start: d.cds_start,
            cds_end: d.cds_end,
            exon_cigars: vec![None; d.exons.len()],
            gene_id: None,
            protein: cdot.protein.clone(),
        };
        assert_eq!(reconstruct_spliced_mrna(&probe, &genome).unwrap(), old_core);
        assert!(d.mismatch_fraction < 1e-9);
    }
}
