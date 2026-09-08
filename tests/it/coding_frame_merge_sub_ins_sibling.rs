//! The `[sub; unchanged; ins]` sibling of the coding-axis codon-frame merge.
//!
//! `coding_frame_merge_axis_asymmetry.rs` pins the `[ins; unchanged; ins]` shape
//! (`c.{p}delinsATC`), and the ruled arm's `TandemDupRun` guard restores it to the
//! shipped `origin/main` (`CanonicalCoalesced`) form. That guard keys on the merge's
//! OUTPUT SIGNATURE — a span-one delins whose alt interior holds the reference base
//! at `ref_start`. That signature is specific to two pure insertions around one
//! retained base.
//!
//! `coalesce_coding_frame_separation` merges ANY pair one unchanged base apart within
//! one codon, not just two insertions. A substitution one base from an insertion
//! (`c.[10G>A;11_12insAC]`) is the same merge on a different shape: the retained base
//! sits at `ref_end`, not `ref_start`, so the original guard's `ref[ref_start]` test
//! misses it and the step-8 `PlusK` peel re-splits the merge into `[sub; dup]`.
//!
//! Measured: on `origin/main`'s shipped `CanonicalCoalesced` arm the merge is kept;
//! the ruled flip split it. This is the same "restore the shipped status quo" mandate
//! as the `[ins; unchanged; ins]` guard (GitHub #1744, a rule-2 preference departure
//! per #1725; `codon-exception-vs-coincidence-carve-out-precedence` rules codon-first),
//! not a new adjudication. Fully hermetic: a `MockProvider`, no `FERRO_MANIFEST`.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

const CODING_TX: &str = "NM_SIB.1";
const GENOMIC_CONTIG: &str = "NC_SIB.1";

/// A single-exon `c.` transcript over `core` with `CDS_START = 1`, so `c.p` is the
/// `p`-th base of `core`. Chosen so codon 4 spans `c.10_12`.
fn coding_provider(core: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    let len = core.len() as u64;
    let pad = "ACGT".repeat(64);
    let g_start = pad.len() as u64 + 1;
    let g_end = pad.len() as u64 + len;
    let transcript = Transcript::new(
        CODING_TX.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        core.to_string(),
        Some(1),
        Some(len),
        vec![Exon::with_genomic(1, 1, len, g_start, g_end)],
        Some(GENOMIC_CONTIG.to_string()),
        Some(g_start),
        Some(g_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    provider.add_genomic_sequence(GENOMIC_CONTIG, format!("{pad}{core}{pad}"));
    provider.add_transcript(transcript);
    provider
}

fn normalize(provider: MockProvider, input: &str) -> String {
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::default());
    let parsed = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    normalizer
        .normalize(&parsed)
        .unwrap_or_else(|e| panic!("normalize {input}: {e}"))
        .to_string()
}

/// The reachable `[sub; unchanged; ins]` case. `core` positions 10/11/12 = `G`/`C`/`A`
/// (codon 4); `c.10G>A` one base from `c.11_12insAC`, retained base `c.11` (`C`).
///
/// On the shipped `CanonicalCoalesced` arm this is the lone `c.10delinsACA`. The ruled
/// flip split it into `c.[10G>A;11_12dup]`; this asserts the merge is kept.
const SPLITTING_CORE: &str = "ATGCCTGAAGCAACGTACGTACGT";
const SPLITTING_INPUT: &str = "NM_SIB.1:c.[10G>A;11_12insAC]";
const MERGED_FORM: &str = "NM_SIB.1:c.10delinsACA";

#[test]
fn a_sub_one_base_from_an_insertion_keeps_the_codon_frame_merge() {
    assert_eq!(
        normalize(coding_provider(SPLITTING_CORE).clone(), SPLITTING_INPUT),
        MERGED_FORM,
        "the [sub; unchanged; ins] codon-frame merge must be kept, matching the shipped \
         CanonicalCoalesced arm — not re-split into [sub; dup] by the PlusK peel",
    );
}

#[test]
fn the_merged_form_is_self_confluent() {
    assert_eq!(
        normalize(coding_provider(SPLITTING_CORE), MERGED_FORM),
        MERGED_FORM,
        "the lone delins spelling must be a fixed point",
    );
}
