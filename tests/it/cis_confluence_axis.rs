//! Confluence over designed multi-member cis alleles (#1443).
//!
//! # The gap
//!
//! `canonicalize_from_sequence` is gated on `members.len() > 1`, so only a
//! multi-member cis allele reaches the partitioner at all. Part 1 of #1443
//! harvested every such input from the bulk corpora and found **650 rows out of
//! 9 949 738** — real corpora are observed variants, where two changes rarely
//! land close enough to share a block. Consumers submitting *designed* alleles
//! have the opposite distribution, so the partitioner is effectively untested
//! and no real corpus can fix that.
//!
//! `examples/generate_cis_confluence_corpus.rs` fills it. Each class it emits is
//! one synthetic reference, one denoted sequence, and N ≥ 2 distinct spellings
//! verified — by an applier independent of the normalizer — to denote that
//! sequence. This module normalizes every spelling in every class and asks the
//! confluence question: **does a class reach exactly one output?**
//!
//! # This is a census, not a pass/fail gate
//!
//! It does not converge today, and a permanently red suite would be worth
//! nothing. So the numbers below are **pinned baselines**: the target is total
//! convergence, and each pinned figure must only ever move in the direction of
//! it. A change that raises a divergence count is a regression and must be
//! explained; a change that lowers one should re-bless the number in the same
//! commit.
//!
//! Read the pins together with `CLAUDE.md`'s note on representation stability.
//! Confluence and stability are different properties, and a fix for the first
//! moves shipped strings — so a commit that lowers a divergence count here still
//! owes the release its `dump_normalized_corpus` measurement.
//!
//! # Sequence preservation is the harder half
//!
//! A class that converged on the *wrong* sequence would be worse than one that
//! diverged, so every output is also applied back to the reference — through
//! `hgvs_to_spdi`, the same way `tests/it/common/synthetic.rs`'s
//! `assert_padded_preserving` does it — and compared with the class's denoted
//! sequence. That count is asserted at zero.
//!
//! # Corpus
//!
//! The corpus is a generated artifact (gitignored, like the spec fixtures) and
//! is regenerated on demand through `common::fixture_gen`, so a fresh checkout
//! just works. Its parameters are what the pins are measured over: regenerating
//! with different `--seeds` re-rolls every number here.
//!
//! # The pins were once per oracle configuration, and no longer are (#1454)
//!
//! Each direction used to carry two pinned censuses, selected by an
//! `expected_census` helper: one for a plain run and one for
//! `FERRO_ASSERT_IDEMPOTENT=1`. Two classes normalized to a non-fixed-point
//! output (#1454), so the oracle turned them into panics that were counted
//! `declined` rather than two distinct outputs, and one pin could only ever be
//! green in one of CI's two jobs.
//!
//! #1454 is fixed, `declined` is 0 in both configurations, and the two censuses
//! are now byte-identical — so there is one constant per direction again.
//!
//! **The collapse did not land where this module predicted, and the difference
//! is the useful part.** The old text expected `converged` to rise by two. It
//! rises by *one* per direction (3': 6632 -> 6633; 5': 6627 -> 6628), and it
//! *falls* by one against the oracle pins (6634 -> 6633; 6629 -> 6628).
//!
//! Measured per class rather than inferred from the net figure — a `+1` is also
//! what "converged three, broke two" looks like, and the two must be told apart.
//! Diffing every class's output set across the fix: **exactly one class changes
//! convergence status, and none regresses.**
//!
//! ```text
//! s04-c-m2-sep3-p4-rot4   the pure #1454 class
//!   before  [13_15delinsTAA;16_19delinsTAAT] | [13_15delinsTAA;16_19inv]
//!   after   [13_15delinsTAA;16_19inv]                            CONVERGED
//!
//! s00-c-m2-sep5-p1-rot4   the issue's own reproducer
//!   before  [10_12delinsTAA;13_15delinsTAT]  | [9dup;15del]
//!   after   [10_12delinsTAA;13_15inv]        | [9dup;15del]      still splits
//! ```
//!
//! The second class had the #1454 defect too, and the fix corrects it — that
//! member is a fixed point now. The class nonetheless still splits, because its
//! competing output is a different **partition** of the same variant rather than
//! a different typing of one member, and no amount of inversion typing merges
//! `[9dup;15del]` with a two-member delins. That residual belongs to the
//! #1235 / #1419-#1421 partition-non-uniqueness family and is deliberately left
//! alone here: choosing between those two partitions moves shipped
//! representations, which is its own measurement.
//!
//! So the `+2` was never reachable by fixing #1454. It was an artifact: under
//! the oracle, `s00-c-m2-sep5-p1-rot4`'s delins spelling *panicked*, leaving
//! `[9dup;15del]` as the only surviving outcome, and a class with one outcome
//! scores as `converged`. The oracle was crediting the defect with a convergence
//! that was really a suppressed output.

use std::collections::BTreeSet;
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

use serde::Deserialize;

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

use crate::common::cis_apply_oracle::apply_with;
use crate::common::fixture_gen;

const CORPUS_RELATIVE_PATH: &str = "tests/fixtures/cis/cis_confluence_corpus.json";

/// Kept in step with `examples/generate_cis_confluence_corpus.rs`. The test
/// rebuilds the providers rather than the corpus carrying 512 bases of padding
/// per class, so the two must agree on the padding, the accessions and the CDS.
const PAD_OFFSET: usize = 256;
const GENOMIC_CONTIG: &str = "NC_TEST.1";
const TX_ACCESSION: &str = "NM_TEST.1";
const TX_CONTIG: &str = "chr_synth";
const CDS_START: u64 = 1;
const CDS_END: u64 = 63;

// ---------------------------------------------------------------------------
// The pinned baseline
// ---------------------------------------------------------------------------

/// The 3'-direction census, pinned. See the module docs: every divergence
/// figure is a baseline that must only ever go down, and `converged` a floor
/// that must only ever go up.
///
/// **58.8% of designed classes converge** (6 633 of 11 272). The 4 639 that do
/// not are the measurement this module exists to make: nothing in the repo
/// could state that number before, because only 650 real rows reach the
/// partitioner at all — and `multi_member_cis_axis` measures 82% convergence
/// over the respellable ones of those, which is what a corpus of far-apart
/// members looks like.
///
/// Note what the remainder is *not*: `declined` and `sequence_changed` are both
/// zero, so every divergence is two well-formed, sequence-preserving outputs
/// for one variant. That is the confluence defect itself, not a parse failure
/// or a corrupted sequence wearing its clothes.
///
/// `converged` was 6 632 before #1454; see the module docs for which single
/// class moved and why the other #1454 class could not.
const THREE_PRIME: Census = Census {
    classes: 11_272,
    spellings: 47_392,
    declined: 0,
    converged: 6_633,
    split_two: 4_505,
    split_three: 115,
    split_more: 19,
    underdetermined: 0,
    sequence_changed: 0,
};

/// The 5'-direction census, pinned. `--direction 5prime` is a supported public
/// option, and confluence is a property of the normalizer rather than of one
/// shuffle direction, so it is measured in full rather than spot-checked.
///
/// It lands within five classes of the 3' figure (6 628 against 6 633), which
/// is worth reading as evidence: the divergence is a property of how the
/// partitioner splits a block, not of which end of an ambiguous run the shuffle
/// walks to. A fix that moved only one of these two numbers would be treating a
/// symptom.
const FIVE_PRIME: Census = Census {
    classes: 11_272,
    spellings: 47_392,
    declined: 0,
    converged: 6_628,
    split_two: 4_530,
    split_three: 105,
    split_more: 9,
    underdetermined: 0,
    sequence_changed: 0,
};

// ---------------------------------------------------------------------------
// Corpus
// ---------------------------------------------------------------------------

/// One confluence class. A structural mirror of the generator's own `Class`;
/// an integration test cannot import a type out of an example, so the two are
/// kept in step by this file failing to deserialize if they drift.
#[derive(Deserialize)]
struct Class {
    id: String,
    axis: String,
    core: String,
    denoted: String,
    members: usize,
    separation: usize,
    spellings: Vec<String>,
}

#[derive(Deserialize)]
struct Corpus {
    classes: Vec<Class>,
}

fn corpus_path() -> PathBuf {
    fixture_gen::fixture_path(CORPUS_RELATIVE_PATH)
}

fn corpus() -> &'static Corpus {
    static CORPUS: OnceLock<Corpus> = OnceLock::new();
    CORPUS.get_or_init(|| {
        let path = corpus_path();
        fixture_gen::ensure_generated_example_fixture(
            &path,
            "generate_cis_confluence_corpus",
            ".cis_confluence_corpus",
            "cis confluence corpus",
            || {},
        );
        let text = fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()));
        serde_json::from_str(&text)
            .unwrap_or_else(|e| panic!("failed to parse {}: {e}", path.display()))
    })
}

// ---------------------------------------------------------------------------
// Providers — the same construction the generator used
// ---------------------------------------------------------------------------

fn padded(core: &str) -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{core}{pad}")
}

/// The provider a class is drawn against, and the full sequence its coordinates
/// address.
fn reference_for(class: &Class) -> (MockProvider, String) {
    if class.axis == "g" {
        let mut provider = MockProvider::new();
        provider.add_genomic_sequence(GENOMIC_CONTIG, padded(&class.core));
        return (provider, padded(&class.core));
    }
    let mut provider = MockProvider::new();
    let tx_len = class.core.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
    let transcript = Transcript::new(
        TX_ACCESSION.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        class.core.clone(),
        Some(CDS_START),
        Some(CDS_END),
        vec![exon],
        Some(TX_CONTIG.to_string()),
        Some(g_start),
        Some(g_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    provider.add_genomic_sequence(TX_CONTIG, padded(&class.core));
    provider.add_transcript(transcript);
    (provider, class.core.clone())
}

/// The sequence a class denotes, on the axis's own reference.
fn expected_sequence(class: &Class) -> String {
    if class.axis == "g" {
        padded(&class.denoted)
    } else {
        class.denoted.clone()
    }
}

// ---------------------------------------------------------------------------
// Measurement
// ---------------------------------------------------------------------------

/// What one direction's run found. Every field is pinned; see the module docs
/// for which way each is allowed to move.
#[derive(Debug, Default, PartialEq, Eq)]
struct Census {
    /// Classes measured.
    classes: usize,
    /// Spellings normalized across all classes.
    spellings: usize,
    /// Spellings normalization declined (or panicked on). Not folded into the
    /// distinct-output count: a decline is a different failure from a
    /// divergence and conflating them would let one hide the other.
    declined: usize,
    /// Classes whose surviving spellings all reached **one** output. The goal
    /// is for this to equal `classes`.
    converged: usize,
    /// Classes that split into exactly two distinct outputs.
    split_two: usize,
    /// … exactly three.
    split_three: usize,
    /// … four or more.
    split_more: usize,
    /// Classes left with fewer than two spellings that normalized at all, so
    /// confluence is not decidable for them. Counted rather than silently
    /// skipped: a rise here would hollow the axis out while every other number
    /// improved.
    underdetermined: usize,
    /// Outputs whose applied sequence differs from the class's. Asserted at
    /// zero — a class that converges on the wrong sequence is worse than one
    /// that diverges.
    sequence_changed: usize,
}

/// A worst-offender line, for the failure message when a pin moves.
struct Divergence {
    id: String,
    outputs: Vec<String>,
}

fn measure(direction: ShuffleDirection) -> (Census, Vec<Divergence>) {
    let classes = &corpus().classes;
    let mut census = Census::default();
    let mut worst: Vec<Divergence> = Vec::new();

    // Classes are sorted by an id that starts with the core index and the axis,
    // so every class sharing a reference is contiguous. Building one provider
    // and one normalizer per group instead of per class is the difference
    // between this running in seconds and in minutes.
    let mut start = 0usize;
    while start < classes.len() {
        let mut end = start + 1;
        while end < classes.len()
            && classes[end].axis == classes[start].axis
            && classes[end].core == classes[start].core
        {
            end += 1;
        }
        let (provider, reference) = reference_for(&classes[start]);
        let normalizer = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::default().with_direction(direction),
        );

        for class in &classes[start..end] {
            census.classes += 1;
            let expected = expected_sequence(class);
            let mut outputs: BTreeSet<String> = BTreeSet::new();
            let mut normalized_spellings = 0usize;
            for spelling in &class.spellings {
                census.spellings += 1;
                let Ok(variant) = parse_hgvs(spelling) else {
                    census.declined += 1;
                    continue;
                };
                let normalized = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    normalizer.normalize(&variant)
                }));
                let output = match normalized {
                    Ok(Ok(value)) => value.to_string(),
                    _ => {
                        census.declined += 1;
                        continue;
                    }
                };
                normalized_spellings += 1;
                if apply_with(&provider, &reference, &output).as_deref() != Some(expected.as_str())
                {
                    census.sequence_changed += 1;
                }
                outputs.insert(output);
            }

            if normalized_spellings < 2 {
                census.underdetermined += 1;
                continue;
            }
            match outputs.len() {
                1 => census.converged += 1,
                2 => census.split_two += 1,
                3 => census.split_three += 1,
                _ => census.split_more += 1,
            }
            if outputs.len() > 1 && worst.len() < 10 {
                worst.push(Divergence {
                    id: class.id.clone(),
                    outputs: outputs.into_iter().collect(),
                });
            }
        }
        start = end;
    }

    (census, worst)
}

fn report(direction: &str, census: &Census, worst: &[Divergence]) -> String {
    let mut out = format!(
        "cis confluence ({direction}): {} classes, {} spellings, {} declined\n  \
         converged: {}\n  split 2: {}\n  split 3: {}\n  split 4+: {}\n  \
         underdetermined: {}\n  sequence changed: {}\n",
        census.classes,
        census.spellings,
        census.declined,
        census.converged,
        census.split_two,
        census.split_three,
        census.split_more,
        census.underdetermined,
        census.sequence_changed,
    );
    for divergence in worst {
        out.push_str(&format!(
            "  {} -> {:?}\n",
            divergence.id, divergence.outputs
        ));
    }
    out
}

/// Assert one direction's census against its pin, printing the measured numbers
/// either way so a moved pin can be re-blessed from the test output.
fn assert_census(direction: ShuffleDirection, label: &str, pinned: &Census) {
    let (measured, worst) = measure(direction);
    println!("{}", report(label, &measured, &worst));
    assert_eq!(
        measured.sequence_changed, 0,
        "{label}: {} normalized outputs no longer denote their class's sequence — a class that \
         converges on the wrong sequence is worse than one that diverges",
        measured.sequence_changed
    );
    assert_eq!(
        &measured,
        pinned,
        "{label}: the confluence census moved. Every divergence figure must only ever go DOWN \
         and `converged` only ever UP; if this change lowers one, re-bless the pin in the same \
         commit and say so in the PR. Measured:\n{}",
        report(label, &measured, &worst)
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// The corpus is only worth measuring if it is actually dense in the shape it
/// claims to cover, so its own contract is checked before anything is read out
/// of it. A generator change that quietly emitted a hundred classes would
/// otherwise show up as an improved convergence rate.
#[test]
fn the_corpus_is_a_dense_set_of_real_confluence_classes() {
    let corpus = corpus();
    assert!(
        corpus.classes.len() >= 5_000,
        "the corpus collapsed to {} classes; it is meant to be thousands",
        corpus.classes.len()
    );
    for class in &corpus.classes {
        assert!(
            class.spellings.len() >= 2,
            "{} is a singleton, which is not a confluence class",
            class.id
        );
        assert_eq!(
            class.spellings.iter().collect::<BTreeSet<_>>().len(),
            class.spellings.len(),
            "{} repeats a spelling",
            class.id
        );
        assert_ne!(
            class.denoted, class.core,
            "{} denotes its own reference, so it describes no variant",
            class.id
        );
    }
    // Both axes, every member count and every separation the generator
    // enumerates must be present — a parameter contributing nothing would make
    // its share of the census meaningless.
    for axis in ["g", "c"] {
        assert!(corpus.classes.iter().any(|c| c.axis == axis));
    }
    for members in [2, 3, 4] {
        assert!(corpus.classes.iter().any(|c| c.members == members));
    }
    for separation in [0, 1, 2, 3, 5, 8] {
        assert!(corpus.classes.iter().any(|c| c.separation == separation));
    }
}

#[test]
fn three_prime_confluence_census() {
    assert_census(ShuffleDirection::ThreePrime, "3prime", &THREE_PRIME);
}

#[test]
fn five_prime_confluence_census() {
    assert_census(ShuffleDirection::FivePrime, "5prime", &FIVE_PRIME);
}
