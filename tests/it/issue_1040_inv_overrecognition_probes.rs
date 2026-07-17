//! Issue #1040: probes for `inv` over-recognition (reverse-complement window
//! carved out of a `delins`), pinning ferro's post-#1036 (v0.8.1) behavior.
//!
//! Background
//! ----------
//! #1040 reports that ferro splits a *contiguous* multi-nucleotide change to
//! expose an internal reverse-complement sub-window as an `inv`, where the
//! whole run should stay a single `delins`. That mechanism was the #1034 bug,
//! fixed by #1036: `rules::decompose_delins` now types a maximal contiguous
//! mismatch run as `inv` only when the *entire* run is a reverse complement
//! (`DNA/inversion.md`; `general.md:34`). A reverse-complement *sub-run* of a
//! longer contiguous change is no longer carved out.
//!
//! The literal anonymized example in #1040 is a single 10-nt contiguous
//! substitution run whose whole-run reverse complement does NOT equal the alt,
//! so on v0.8.1 it already normalizes to one `delins` — the spec-correct form
//! mutalyzer produces. This suite pins that, and probes the surrounding
//! boundary of the `inv` typer so any *future* over-recognition regression is
//! caught and so the true residual divergence surface is documented.
//!
//! Two `inv`-producing code paths exist, and both are probed here:
//!   * `decompose_delins` — length-preserving compound runs (the #1040 path).
//!   * `canonicalize_delins` step 4d — a single span; types `inv` only when the
//!     *shared-affix-trimmed* span is a whole reverse complement. Its trimmed
//!     ref and alt must be equal length, so a genuinely length-changing delins
//!     can never become an `inv`, and `shorten_inversion`'s outer-pair peeling
//!     is a provable no-op post-trim (a peelable outer pair would already have
//!     been trimmed as a shared affix). Length-changing revcomp probes are
//!     therefore expected to stay `delins`.
//!
//! Spec basis for the pinned outputs (`DNA/inversion.md`):
//!   * an inversion is "more than one nucleotide" that is the reverse
//!     complement — so a 2-nt reverse complement is a *valid* `inv`;
//!   * "two variants separated by one or more nucleotides should be described
//!     individually and not as a delins" — the delins exception is codon-only
//!     (two variants separated by one nucleotide affecting one amino acid), so
//!     on a genomic axis reverse-complement runs separated by an unchanged base
//!     are individual `inv`s, not one merged `delins`.
//!
//! Every assertion below is ferro's actual v0.8.1 output and is HGVS-spec
//! compliant. The remaining ferro-vs-mutalyzer `inv` divergence candidate is
//! the 2-nt inversion (P3/P4/P5): ferro emits `inv` per the spec's ">1 nt"
//! rule, whereas some tools render a 2-nt reverse complement as a `delins`.
//! That is a spec-vs-tool typing choice, not the #1040 sub-run carving, which
//! is gone.
//!
//! 3'-end code path
//! ----------------
//! The related #1041 investigation found that a variant within `window_size`
//! (100 bp) of a short contig's 3' end took a *different* code path on v0.8.1:
//! the fetch window overran EOF, the read errored, and the variant fell into
//! the minimal-notation fallback that runs neither 3' shifting nor the
//! `delins -> inv/sub/dup` canonicalization (fixed for #1041 by #1042's window
//! clamp). Because the real ferro-vs-mutalyzer comparison sites variants on
//! short `NG_` contigs, the near-3'-end path is the one that actually runs
//! there — so the last section pins the #1040 example on a short contig too.
//! Near the 3' end the fallback simply leaves a `delins` as a `delins`, which
//! is the correct answer for #1040 (the fallback cannot *carve* an `inv`); the
//! result is a single `delins` on both v0.8.1 and the #1042 branch.

use crate::common::synthetic::{normalize_to_string, SyntheticBuilder};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

const G: &str = "NC_TEST.1";
const C: &str = "NM_TEST.1";

/// Build a genomic provider whose core sits at 1-based position 257 and
/// normalize `input` against it.
fn g(core: &str, input: &str) -> String {
    normalize_to_string(SyntheticBuilder::genomic(core).build(), input)
}

/// Build a CDS provider (whole core is the CDS) and normalize `input`.
fn c(core: &str, input: &str) -> String {
    let cds_end = core.len() as u64;
    normalize_to_string(
        SyntheticBuilder::cds(core, 1, cds_end, Strand::Plus).build(),
        input,
    )
}

// ---------------------------------------------------------------------------
// The #1040 report: a contiguous run must stay a single `delins` (FIXED).
// ---------------------------------------------------------------------------

#[test]
fn issue_1040_literal_ten_nt_contiguous_run_stays_delins() {
    // ref[257..266] = CAGTGACTAG (10 explicit substitutions, all mismatch).
    //   revcomp(CAGTGACTAG) = CTAGTCACTG != TGTCACGACT  => NOT a whole-run inv.
    //   Interior sub-window 259_262 (GTGA -> TCAC) IS a revcomp but must not be
    //   carved out of the contiguous run.
    // Spec-correct (and what mutalyzer produces): one delins.
    assert_eq!(
        g(
            "CAGTGACTAG",
            &format!(
                "{G}:g.[257C>T;258A>G;259G>T;260T>C;261G>A;262A>C;263C>G;264T>A;265A>C;266G>T]"
            ),
        ),
        format!("{G}:g.257_266delinsTGTCACGACT"),
    );
}

#[test]
fn issue_1040_same_run_as_delins_input_stays_delins() {
    // The identical change expressed as a delins rather than an allele of subs.
    assert_eq!(
        g("CAGTGACTAG", &format!("{G}:g.257_266delinsTGTCACGACT")),
        format!("{G}:g.257_266delinsTGTCACGACT"),
    );
}

#[test]
fn issue_1034_minimal_contiguous_run_stays_delins() {
    // ref CTG -> alt ACA. revcomp(CTG)=CAG != ACA; interior TG->CA is a revcomp
    // sub-run that must not be carved. The #1034 minimal case, both input forms.
    assert_eq!(
        g("CTG", &format!("{G}:g.257_259delinsACA")),
        format!("{G}:g.257_259delinsACA"),
    );
    assert_eq!(
        g("CTG", &format!("{G}:g.[257C>A;258T>C;259G>A]")),
        format!("{G}:g.257_259delinsACA"),
    );
}

// ---------------------------------------------------------------------------
// Guard: a genuine whole-run reverse complement STILL emits `inv`.
// (Under-recognizing these would be the opposite defect.)
// ---------------------------------------------------------------------------

#[test]
fn whole_run_reverse_complement_emits_inv() {
    // 3-nt whole-run revcomp: GCT -> AGC (= revcomp(GCT)).
    assert_eq!(
        g("GCT", &format!("{G}:g.257_259delinsAGC")),
        format!("{G}:g.257_259inv"),
    );
    // 6-nt whole-run revcomp: ACGTGC -> GCACGT (= revcomp(ACGTGC)).
    assert_eq!(
        g("ACGTGC", &format!("{G}:g.257_262delinsGCACGT")),
        format!("{G}:g.257_262inv"),
    );
}

// ---------------------------------------------------------------------------
// Boundary probes: `canonicalize_delins` step 4d + `shorten_inversion`.
// ---------------------------------------------------------------------------

#[test]
fn shared_affix_trimming_reveals_inner_inv() {
    // ACGAGT -> ACTCGT: shared prefix "AC", shared suffix "GT"; the trimmed
    // core GA -> TC is a whole revcomp, so the inv is the trimmed 259_260 span.
    // This exercises the trimming path; the inv is spec-correct (whole trimmed
    // span is the reverse complement).
    assert_eq!(
        g("ACGAGT", &format!("{G}:g.257_262delinsACTCGT")),
        format!("{G}:g.259_260inv"),
    );
}

#[test]
fn length_changing_revcomp_delins_stays_delins() {
    // A revcomp-looking change whose ref/alt differ in length can NOT be typed
    // `inv` (4d requires equal trimmed lengths). ref GCT (257_259) -> alt AG:
    //   revcomp(GCT) = AGC, and the alt "AG" is a truncated (length-changing)
    //   reverse complement. No shared affix reduces it to a sub/del/ins/dup, so
    //   it stays a delins rather than being over-typed as an `inv`.
    assert_eq!(
        g("GCT", &format!("{G}:g.257_259delinsAG")),
        format!("{G}:g.257_259delinsAG"),
    );
}

// ---------------------------------------------------------------------------
// 2-nt inversions: the remaining ferro-vs-tool divergence candidate.
// Per DNA/inversion.md an inversion is ">1 nt", so a 2-nt reverse complement
// is a valid `inv`. ferro emits `inv`; some tools render it as a `delins`.
// ---------------------------------------------------------------------------

#[test]
fn two_nt_reverse_complement_emits_inv() {
    // TG -> CA (= revcomp(TG)).
    assert_eq!(
        g("TG", &format!("{G}:g.257_258delinsCA")),
        format!("{G}:g.257_258inv"),
    );
    // GG -> CC (= revcomp(GG)).
    assert_eq!(
        g("GG", &format!("{G}:g.257_258delinsCC")),
        format!("{G}:g.257_258inv"),
    );
}

// ---------------------------------------------------------------------------
// Reverse-complement runs separated by an unchanged base: individual `inv`s on
// a genomic axis (general.md:34 — the delins exception is codon-only).
// ---------------------------------------------------------------------------

#[test]
fn revcomp_runs_separated_by_identity_are_individual_on_genomic_axis() {
    // ref TGACA -> alt CAATG: 257_258 (TG->CA) and 260_261 (CA->TG) are each a
    // revcomp, separated by the unchanged base 259 (A). On a genomic axis they
    // are described individually, not merged into one delins.
    assert_eq!(
        g("TGACA", &format!("{G}:g.257_261delinsCAATG")),
        format!("{G}:g.[257_258inv;260_261inv]"),
    );
}

#[test]
fn revcomp_runs_in_distinct_codons_are_individual_on_cds_axis() {
    // Same change on a CDS axis: the two invs fall in *different* codons
    // (1_2 in codon 1, 4_5 in codon 2), so the codon-frame delins exception
    // (two variants separated by one nt affecting one amino acid) does not
    // apply and they stay individual.
    assert_eq!(
        c("TGACA", &format!("{C}:c.1_5delinsCAATG")),
        format!("{C}:c.[1_2inv;4_5inv]"),
    );
}

// ---------------------------------------------------------------------------
// #1040 near a short contig's 3' end — the path the real comparison exercises.
// ---------------------------------------------------------------------------

/// Genome-capable provider: one `contig` of `len` cyclic ACGT bytes with
/// `payload` written 1-based at `pos1` (mirrors the #1041 repro harness).
fn provider_len(contig: &str, len: usize, pos1: usize, payload: &str) -> JsonProvider {
    let mut bases: Vec<u8> = "ACGT".bytes().cycle().take(len).collect();
    for (i, b) in payload.bytes().enumerate() {
        bases[pos1 - 1 + i] = b;
    }
    let seq = String::from_utf8(bases).unwrap();
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

fn norm(p: &JsonProvider, input: &str) -> String {
    Normalizer::new(p.clone())
        .normalize(&parse_hgvs(input).unwrap())
        .unwrap()
        .to_string()
}

#[test]
fn issue_1040_contiguous_run_near_3prime_end_stays_delins() {
    // ref[200..209] = CAGTGACTAG on a 250 bp contig: the variant end (209) plus
    // window_size (100) = 309 overruns the contig, so this takes the near-3'-end
    // path that the real comparison exercises. It must still be a single delins,
    // not a carved `[..delins;..inv;..delins]`. Holds on v0.8.1 and post-#1042.
    let near = provider_len("c", 250, 200, "CAGTGACTAG");
    assert_eq!(
        norm(&near, "c:g.200_209delinsTGTCACGACT"),
        "c:g.200_209delinsTGTCACGACT",
    );
}

#[test]
fn issue_1040_compound_allele_near_3prime_end_stays_delins() {
    // The reported input shape (a compound of adjacent substitutions), sited
    // near the 3' end. It collapses to the same single delins, not a carve.
    let near = provider_len("c", 250, 200, "CAGTGACTAG");
    assert_eq!(
        norm(
            &near,
            "c:g.[200C>T;201A>G;202G>T;203T>C;204G>A;205A>C;206C>G;207T>A;208A>C;209G>T]",
        ),
        "c:g.200_209delinsTGTCACGACT",
    );
}
