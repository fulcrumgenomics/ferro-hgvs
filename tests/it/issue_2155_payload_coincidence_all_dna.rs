//! #2155 — a contiguous payload-coincidence change collapses to one spanning
//! `delins` on every cis-reachable DNA axis (`c./g./m./n.`) and stays split on
//! `r.` (RNA is out of the DNA carve-out's jurisdiction). See the ruling
//! `delins-payload-coincidence-carve-out-is-coding-dna-scoped` (superseded to
//! all-DNA by this change) and its siblings.

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

const PAD: usize = 256;
// Block at 1-based 10_17 = CTTAGTTA -> AAACAAAC (equal length, distance 7,
// two interior A's forced-unchanged in all 3 minimal alignments).
const CORE: &str = "AGAACCCCCCTTAGTTAAGAACAAAAGCAACAATCTTCGTGGTCCTGG";
const CONTIG: &str = "NC_TEST.1";
const SINGLE_DELINS_PAYLOAD: &str = "delinsAAACAAAC";
const FRAGMENTED_TAIL: &str = "[10_12delinsAA;14_16delinsCAA;17_18insC]";

fn padded() -> String {
    let pad = "ACGT".repeat(PAD / 4);
    format!("{pad}{CORE}{pad}")
}

fn tx(acc: &str, cds: Option<(u64, u64)>) -> Transcript {
    let n = CORE.len() as u64;
    let (gs, ge) = (PAD as u64 + 1, PAD as u64 + n);
    Transcript::new(
        acc.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        CORE.to_string(),
        cds.map(|c| c.0),
        cds.map(|c| c.1),
        vec![Exon::with_genomic(1, 1, n, gs, ge)],
        Some(CONTIG.to_string()),
        Some(gs),
        Some(ge),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    )
}

fn norm(p: MockProvider, input: &str) -> String {
    let nz = Normalizer::with_config(p, NormalizeConfig::default());
    nz.normalize(&parse_hgvs(input).unwrap())
        .unwrap()
        .to_string()
}

fn with_contig() -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence(CONTIG, padded());
    p
}

/// The `DNA/delins.md:47` payload-coincidence carve-out, widened from
/// coding-DNA-only to every cis-reachable DNA axis. `c.`, `n.` and `g.` all
/// collapse the fragmented split into one spanning `delins`; `r.` — outside the
/// carve-out's jurisdiction, since a `DNA/` clause cannot scope the RNA axis —
/// stays fragmented.
#[test]
fn dna_axes_collapse_the_contiguous_change_rna_does_not() {
    // c.: CDS_START = 1 so c.p == the transcript's p-th base.
    let mut pc = with_contig();
    pc.add_transcript(tx("NM_TEST.1", Some((1, 48))));
    assert_eq!(
        norm(pc, &format!("NM_TEST.1:c.10_17{SINGLE_DELINS_PAYLOAD}")),
        "NM_TEST.1:c.10_17delinsAAACAAAC",
    );
    // n.: same transcript, no CDS -> collapses under the widened scope.
    let mut pn = with_contig();
    pn.add_transcript(tx("NR_TEST.1", None));
    assert_eq!(
        norm(pn, &format!("NR_TEST.1:n.10_17{SINGLE_DELINS_PAYLOAD}")),
        "NR_TEST.1:n.10_17delinsAAACAAAC",
    );
    // g.: block at genomic 266_273.
    assert_eq!(
        norm(
            with_contig(),
            &format!("{CONTIG}:g.266_273{SINGLE_DELINS_PAYLOAD}")
        ),
        "NC_TEST.1:g.266_273delinsAAACAAAC",
    );
    // m.: same synthetic contig, addressed on the mitochondrial axis. `m.` is
    // selected purely by the coordinate letter (`GenomeKind::Mt`, shared with
    // `g.`'s allele grammar) and carries no accession requirement, so the same
    // hermetic contig serves it -> collapses under the widened scope exactly
    // like `g.`.
    assert_eq!(
        norm(
            with_contig(),
            &format!("{CONTIG}:m.266_273{SINGLE_DELINS_PAYLOAD}")
        ),
        "NC_TEST.1:m.266_273delinsAAACAAAC",
    );
    // r.: RNA stays split — jurisdiction, DNA clause cannot scope it.
    let mut pr = with_contig();
    pr.add_transcript(tx("NR_RNA.1", None));
    let r_out = norm(
        pr,
        &format!("NR_RNA.1:r.10_17{}", SINGLE_DELINS_PAYLOAD.to_lowercase()),
    );
    // Exact pin: `r.` is outside `DNA/delins.md:47`'s jurisdiction, so the
    // authored members stay individual. A `contains(';')` fallback would pass for
    // any multi-member output — a wrong span or payload included — which is the
    // shape this row exists to catch. Re-pin only with a named ruling.
    assert_eq!(
        r_out,
        format!("NR_RNA.1:r.{}", FRAGMENTED_TAIL.to_lowercase()),
        "r. must stay fragmented",
    );
}
