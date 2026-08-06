//! Confluence over designed multi-member cis alleles on the `n.` and `r.` axes.
//!
//! # The gap this closes
//!
//! [`super::cis_confluence_axis`] censuses the same generator's `g.` and `c.`
//! output. Those were the two axes #1443 shipped, and they leave the other half
//! of `Normalizer::is_splittable_single_member`'s gate unmeasured: that gate
//! admits `HV::Genome` and `HV::Mt` only, so widening it to `c.`/`n.`/`r.` moves
//! the `n.` and `r.` axes with no committed census to move *against*. The same
//! is true of `FERRO_PARTITION`: the partitioner is axis-blind, but "axis-blind
//! in the source" and "measured on that axis" are different claims, and only the
//! second is evidence.
//!
//! So this module is the `n.`/`r.` half, deliberately kept as a separate corpus
//! and a separate set of pins rather than folded into the `g,c` ones. Folding
//! them in would have re-rolled every number `cis_confluence_axis` pins, which
//! is a change to an existing baseline dressed up as new coverage.
//!
//! # Everything else is the `g,c` module's contract
//!
//! Same generator, same class shape, same census fields, same
//! must-only-move-one-way discipline, same zero-tolerance on
//! `sequence_changed`. Read [`super::cis_confluence_axis`]'s module docs for
//! all of it; only what differs is restated here.
//!
//! # What differs
//!
//! * **The reference.** Both axes are drawn against one `NR_` transcript with
//!   **no CDS**, where `g.`/`c.` use a padded contig and a CDS-bearing `NM_`.
//! * **The alphabet.** An `r.` spelling is lowercase over `acgu`. It is only a
//!   rendering difference: `hgvs_to_spdi` rewrites `u` to `T` before any
//!   sequence comparison, so the stored sequence is DNA on both axes.
//! * **The census is per axis.** `g,c` pins one figure across both of its axes;
//!   here each axis is pinned separately, because the question the axis gate
//!   asks is whether `n.` and `r.` behave *like each other* and like `c.`. A
//!   combined number could hide one axis regressing while the other improved.
//! * **One pin per axis and direction, not two.** `g,c` carries a second census
//!   per direction for `FERRO_ASSERT_IDEMPOTENT=1`, because two of its `c.`
//!   classes normalize to a non-fixed-point output (#1454) and the oracle turns
//!   those into `declined` panics. Neither axis here has such a class: measured
//!   under all three oracles the four censuses below are unchanged, so there is
//!   no second configuration to pin and no `expected_census` selector. If a
//!   future change gives `n.` or `r.` a non-fixed-point output, these pins go
//!   red in CI's `test-oracle` job while staying green in `Test` — which is the
//!   signal to add the split rather than to re-bless one number.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

use serde::Deserialize;

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

use crate::common::cis_apply_oracle::apply_with;
use crate::common::fixture_gen;

const CORPUS_RELATIVE_PATH: &str = "tests/fixtures/cis/cis_confluence_nr_corpus.json";

/// Kept in step with `examples/generate_cis_confluence_corpus.rs`, exactly as
/// [`super::cis_confluence_axis`] keeps its own copies: the test rebuilds the
/// provider rather than the corpus carrying 512 bases of padding per class.
const PAD_OFFSET: usize = 256;
const NR_ACCESSION: &str = "NR_TEST.1";
const TX_CONTIG: &str = "chr_synth";

// ---------------------------------------------------------------------------
// The pinned baseline
// ---------------------------------------------------------------------------

/// The `n.` 3'-direction census, pinned. See the module docs: every divergence
/// figure is a baseline that must only ever go down, and `converged` a floor
/// that must only ever go up.
///
/// **46.6% of designed `n.` classes converge** (2 628 of 5 636), against the
/// 58.8% `cis_confluence_axis` pins across its combined `g,c` corpus. The
/// combined figure is the thing to distrust: measured per axis on that same
/// corpus, `g.` converges 4 000 of 5 636 (71.0%) and `c.` 2 632 of 5 636
/// (46.7%), so 58.8% is the average of two populations that do not resemble
/// each other.
///
/// Lined up, the four axes are **`g.` 71.0%, `c.` 46.7%, `n.` 46.6%, `r.`
/// 46.6%** — one outlier and a tight cluster, and `g.` is precisely the axis
/// `Normalizer::is_splittable_single_member` admits. That is what the axis
/// gate is worth, stated as a number rather than as a mechanism: the corpus
/// carries a lone spanning `delins` spelling for every class, and on `g.` it is
/// re-derived to agree with the bracketed spellings while on `c.`/`n.`/`r.` it
/// stays whole and diverges from them.
///
/// Note what the remainder is *not*: `declined` and `sequence_changed` are both
/// zero, so every divergence is two well-formed, sequence-preserving outputs for
/// one variant — the confluence defect itself, not a parse failure or a
/// corrupted sequence wearing its clothes.
const N_THREE_PRIME: Census = Census {
    classes: 5_636,
    spellings: 23_696,
    declined: 0,
    converged: 2_628,
    split_two: 2_941,
    split_three: 58,
    split_more: 9,
    underdetermined: 0,
    sequence_changed: 0,
};

/// The `r.` 3'-direction census, pinned. See [`N_THREE_PRIME`].
const R_THREE_PRIME: Census = Census {
    classes: 5_636,
    spellings: 23_696,
    declined: 0,
    converged: 2_627,
    split_two: 2_941,
    split_three: 59,
    split_more: 9,
    underdetermined: 0,
    sequence_changed: 0,
};

/// The `n.` 5'-direction census, pinned. `--direction 5prime` is a supported
/// public option and confluence is a property of the normalizer rather than of
/// one shuffle direction, so it is measured in full rather than spot-checked.
const N_FIVE_PRIME: Census = Census {
    classes: 5_636,
    spellings: 23_696,
    declined: 0,
    converged: 2_625,
    split_two: 2_954,
    split_three: 53,
    split_more: 4,
    underdetermined: 0,
    sequence_changed: 0,
};

/// The `r.` 5'-direction census, pinned. See [`N_FIVE_PRIME`].
const R_FIVE_PRIME: Census = Census {
    classes: 5_636,
    spellings: 23_696,
    declined: 0,
    converged: 2_624,
    split_two: 2_954,
    split_three: 54,
    split_more: 4,
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
            &["--axes", "n,r"],
            ".cis_confluence_nr_corpus",
            "cis confluence n./r. corpus",
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
///
/// One construction serves both axes: an `NR_` record stores DNA whichever way
/// a spelling renders it, so an `r.` class differs from its `n.` twin only in
/// the spellings, never in the reference. Mirrors the generator's
/// `noncoding_provider`.
fn reference_for(class: &Class) -> (MockProvider, String) {
    let mut provider = MockProvider::new();
    let tx_len = class.core.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
    let transcript = Transcript::new(
        NR_ACCESSION.to_string(),
        Some("SYNTH_NR".to_string()),
        Strand::Plus,
        class.core.clone(),
        None,
        None,
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

// ---------------------------------------------------------------------------
// Measurement
// ---------------------------------------------------------------------------

/// What one axis in one direction found. Field-for-field the same shape as
/// [`super::cis_confluence_axis`]'s `Census`, so the two sets of numbers are
/// directly comparable; see that module for what each field means and which way
/// it is allowed to move.
#[derive(Debug, Default, PartialEq, Eq)]
struct Census {
    classes: usize,
    spellings: usize,
    declined: usize,
    converged: usize,
    split_two: usize,
    split_three: usize,
    split_more: usize,
    underdetermined: usize,
    sequence_changed: usize,
}

/// A worst-offender line, for the failure message when a pin moves.
struct Divergence {
    id: String,
    outputs: Vec<String>,
}

fn measure(direction: ShuffleDirection, axis: &str) -> (Census, Vec<Divergence>) {
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
        if classes[start].axis != axis {
            start = end;
            continue;
        }
        let (provider, reference) = reference_for(&classes[start]);
        let normalizer = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::default().with_direction(direction),
        );

        for class in &classes[start..end] {
            census.classes += 1;
            let expected = class.denoted.clone();
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

fn report(label: &str, census: &Census, worst: &[Divergence]) -> String {
    let mut out = format!(
        "cis confluence ({label}): {} classes, {} spellings, {} declined\n  \
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

/// Assert one axis-and-direction census against its pin, printing the measured
/// numbers either way so a moved pin can be re-blessed from the test output.
fn assert_census(direction: ShuffleDirection, axis: &str, label: &str, pinned: &Census) {
    let (measured, worst) = measure(direction, axis);
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
fn the_nr_corpus_is_a_dense_set_of_real_confluence_classes() {
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
    // This corpus must be the `n,r` one and nothing else. A stale `g,c` artifact
    // installed under this path would make every census below measure zero
    // classes, which `assert_census` does catch — it compares the whole `Census`
    // for equality rather than only checking the one-way movement its failure
    // message describes, so `classes: 0` cannot pass as "went down". What this
    // guard buys is therefore one failure that names the wrong artifact, instead
    // of four opaque census mismatches to work backwards from.
    assert!(
        corpus
            .classes
            .iter()
            .all(|c| c.axis == "n" || c.axis == "r"),
        "the corpus carries an axis other than n./r.; it is the wrong artifact"
    );
    for axis in ["n", "r"] {
        assert!(corpus.classes.iter().any(|c| c.axis == axis));
    }
    for members in [2, 3, 4] {
        assert!(corpus.classes.iter().any(|c| c.members == members));
    }
    for separation in [0, 1, 2, 3, 5, 8] {
        assert!(corpus.classes.iter().any(|c| c.separation == separation));
    }
}

/// `r.` spellings must actually be RNA-rendered, or the axis is `n.` twice over
/// and its census is a duplicate rather than a measurement.
///
/// The two per-spelling checks at the bottom are both **negative** — no
/// uppercase, no `t` — and negatives alone cannot tell an RNA rendering from a
/// lowercasing that never mapped `t` to `u`: that carries no `t` and no capital
/// either, and it is precisely the "`n.` twice over" outcome above. Nor does
/// looking for a `u` fix it, because `dup` supplies one in every spelling that
/// names a duplication.
///
/// So the positive half is asserted against the `n.` twin, which is the actual
/// contract: an `r.` class differs from the `n.` class with the same id only in
/// how its bases are spelled, never in its positions, its members or its
/// reference. `mapped_a_base` then pins that the `t` -> `u` arm was really
/// taken, so a corpus whose payloads happened to be `t`-free could not pass this
/// vacuously.
#[test]
fn the_rna_axis_is_spelled_in_the_rna_alphabet() {
    let corpus = corpus();
    let rna: Vec<&Class> = corpus.classes.iter().filter(|c| c.axis == "r").collect();
    assert!(!rna.is_empty(), "no r. classes to check");

    let by_id: BTreeMap<&str, &Class> = corpus
        .classes
        .iter()
        .filter(|c| c.axis == "n")
        .map(|c| (c.id.as_str(), c))
        .collect();
    let mut mapped_a_base = false;
    for class in &rna {
        let twin_id = class.id.replacen("-r-", "-n-", 1);
        let twin = by_id.get(twin_id.as_str()).unwrap_or_else(|| {
            panic!(
                "{} has no n. twin at {twin_id}, so the two axes are not the same corpus",
                class.id
            )
        });
        assert_eq!(
            class.spellings.len(),
            twin.spellings.len(),
            "{} and {twin_id} are one design on two axes but carry different \
             spelling counts",
            class.id
        );
        for (rendered, dna) in class.spellings.iter().zip(&twin.spellings) {
            let (_, dna_body) = dna
                .split_once(":n.")
                .unwrap_or_else(|| panic!("{dna} is not an n. description"));
            if dna_body.contains('T') {
                mapped_a_base = true;
            }
            // Safe to fold the whole body: positions are digits and every edit
            // keyword (`del`, `dup`, `ins`, `delins`) is already lowercase and
            // `t`-free, so only bases are touched.
            let expected: String = dna_body
                .to_ascii_lowercase()
                .chars()
                .map(|base| if base == 't' { 'u' } else { base })
                .collect();
            let (_, rna_body) = rendered
                .split_once(":r.")
                .unwrap_or_else(|| panic!("{rendered} is not an r. description"));
            assert_eq!(
                rna_body, expected,
                "{rendered} is not {dna} with its bases spelled in acgu"
            );
        }
    }
    assert!(
        mapped_a_base,
        "no n. spelling carried a `T`, so the t -> u mapping was never exercised \
         and these r. classes could be their n. twins with only the case changed"
    );

    for class in rna {
        for spelling in &class.spellings {
            let (_, description) = spelling
                .split_once(":r.")
                .unwrap_or_else(|| panic!("{spelling} is not an r. description"));
            assert!(
                !description.chars().any(|c| c.is_ascii_uppercase()),
                "{spelling} carries an uppercase base, so it is not RNA-rendered"
            );
            assert!(
                !description.contains('t'),
                "{spelling} carries a `t`, which RNA spells `u`"
            );
        }
    }
}

#[test]
fn noncoding_three_prime_confluence_census() {
    assert_census(
        ShuffleDirection::ThreePrime,
        "n",
        "n. 3prime",
        &N_THREE_PRIME,
    );
}

#[test]
fn rna_three_prime_confluence_census() {
    assert_census(
        ShuffleDirection::ThreePrime,
        "r",
        "r. 3prime",
        &R_THREE_PRIME,
    );
}

#[test]
fn noncoding_five_prime_confluence_census() {
    assert_census(ShuffleDirection::FivePrime, "n", "n. 5prime", &N_FIVE_PRIME);
}

#[test]
fn rna_five_prime_confluence_census() {
    assert_census(ShuffleDirection::FivePrime, "r", "r. 5prime", &R_FIVE_PRIME);
}
