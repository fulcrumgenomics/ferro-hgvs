//! Dump the **divergent** confluence classes measured by
//! `tests/it/cis_confluence_axis.rs`, one record per class, and census them by
//! *shape*.
//!
//! # Why
//!
//! The axis test pins a census — 11 272 classes, 6 629 converged in the 3'
//! direction — and prints ten worst offenders. That is enough to detect a
//! regression and nowhere near enough to *adjudicate* the 4 643 that diverge.
//! Adjudicating them one row at a time is not feasible and is not the right
//! unit anyway: the prior campaign
//! (`campaigns/2026-08-06-separation-rule-adjudication.md`) collapsed 1 957 rows
//! to ~32 shapes with the top nine covering 91%, and the same collapse works
//! here.
//!
//! So this emits, per divergent class, every spelling and what it normalized
//! to, plus a **shape key** derived from the disagreement itself rather than
//! from the class's design parameters:
//!
//! ```text
//! 1×delins | 2×[delins,delins]
//! ```
//!
//! meaning one spelling reached a single spanning `delins` and another reached
//! a two-member allele of two `delins`. Shapes are what a spec clause can be
//! applied to; individual rows are not.
//!
//! # Geometry, and the arithmetic trap it exists to avoid
//!
//! An insertion consumes **no** reference position. Each member is modelled as
//! a closed interval `[lo, hi]` over the reference, with an insertion left
//! *empty* (`lo = A + 1`, `hi = A`), so that
//!
//! ```text
//! separation = next.lo - prev.hi - 1
//! ```
//!
//! holds uniformly across every edit kind. Treating an insertion's `A_B` anchor
//! as a consumed two-base span double-counts and has already invalidated one
//! published distribution.
//!
//! The intervals come from `hgvs_to_spdi`, which applies the description to the
//! reference — the same normalizer-independent oracle the corpus generator used
//! to establish ground truth. Nothing here asks the normalizer what a member
//! covers.
//!
//! # Usage
//!
//! ```text
//! cargo run --features dev --example dump_confluence_divergences -- --stats
//! cargo run --features dev --example dump_confluence_divergences -- --out /tmp/div.json
//! cargo run --features dev --example dump_confluence_divergences -- --direction 5prime --stats
//! ```

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use serde::{Deserialize, Serialize};

use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

/// Kept in step with `examples/generate_cis_confluence_corpus.rs` and
/// `tests/it/cis_confluence_axis.rs`; all three must agree on the padding, the
/// accessions and the CDS or the coordinates mean different things.
const PAD_OFFSET: usize = 256;
const GENOMIC_CONTIG: &str = "NC_TEST.1";
const TX_ACCESSION: &str = "NM_TEST.1";
const TX_CONTIG: &str = "chr_synth";
const CDS_START: u64 = 1;
const CDS_END: u64 = 63;

const DEFAULT_CORPUS: &str = "tests/fixtures/cis/cis_confluence_corpus.json";

#[derive(Parser, Debug)]
#[command(about = "Dump and census the divergent cis confluence classes")]
struct Cli {
    /// The corpus written by `generate_cis_confluence_corpus`.
    #[arg(long, default_value = DEFAULT_CORPUS)]
    corpus: PathBuf,
    /// Write the per-class records here as JSON.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Print the shape census only.
    #[arg(long)]
    stats: bool,
    /// `3prime` (default) or `5prime`.
    #[arg(long, default_value = "3prime")]
    direction: String,
    /// How many shapes to list in the census.
    #[arg(long, default_value_t = 40)]
    top: usize,
    /// How many exemplar classes to print per shape.
    #[arg(long, default_value_t = 0)]
    examples: usize,
}

// ---------------------------------------------------------------------------
// Corpus
// ---------------------------------------------------------------------------

#[derive(Deserialize)]
struct Class {
    id: String,
    axis: String,
    core: String,
    denoted: String,
    members: usize,
    separation: usize,
    payload: usize,
    pattern: String,
    block_start: usize,
    reference_block: String,
    alt_block: String,
    spellings: Vec<String>,
}

#[derive(Deserialize)]
struct Corpus {
    classes: Vec<Class>,
}

/// The confluence question asked one level up: a **variant** is one reference
/// and one denoted sequence, and several classes can denote the same one by
/// different designs. `cis_confluence_axis` compares outputs only *within* a
/// class, so a variant two classes describe differently is invisible to it.
///
/// This matters beyond bookkeeping, because the design parameter the separation
/// rule keys on is exactly what differs between such classes: two members one
/// nucleotide apart and four members zero apart can denote the same sequence.
/// `(axis, reference core, denoted sequence)` — one *variant*, regardless of how
/// many designs reach it.
type VariantKey = (String, String, String);

/// Every output any class of one variant reached, whether every one of those
/// classes converged internally, and how many classes there were.
type VariantOutputs = (BTreeSet<String>, bool, usize);

/// Everything the pass declined to measure.
///
/// Printed unconditionally by [`render_census`], and every field is expected to
/// be **zero**: the spellings come from the corpus generator and the outputs
/// from a successful `normalize`. That is the argument for printing them, not
/// for omitting them. A silently-dropped spelling leaves its class with fewer
/// surviving outputs, and a class reduced to one survivor is then recorded as
/// *converged* — so an uncounted drop does not merely shrink the sample, it
/// moves a class from the divergent side of the census to the converged side.
/// A printed `0` turns that into evidence; silence is indistinguishable from a
/// run where it was not zero.
#[derive(Default)]
struct Drops {
    /// Spellings `parse_hgvs` rejected.
    unparsed_spellings: usize,
    /// Spellings normalization errored or panicked on.
    unnormalized_spellings: usize,
    /// Divergent classes excluded because some output's member spans could not
    /// be resolved, so the class's geometry was never measured. Excluded rather
    /// than classified: see [`members_of`].
    unresolved_classes: usize,
}

struct CrossClass {
    /// Distinct `(axis, reference, denoted)` variants in the corpus.
    variants: usize,
    /// Variants reachable from more than one class.
    multi_class: usize,
    /// Multi-class variants whose classes do not agree on a single output.
    multi_class_divergent: usize,
    /// … of those, the ones where **every** class converged internally, so the
    /// disagreement is visible only across classes.
    hidden_by_the_class_boundary: usize,
}

// ---------------------------------------------------------------------------
// Records
// ---------------------------------------------------------------------------

/// One normalized output, with the geometry the adjudication turns on.
#[derive(Serialize, Clone)]
struct Output {
    text: String,
    /// Number of cis members.
    arity: usize,
    /// Per-member edit kind. Ordered by [`Output::spans`], not by written
    /// order — the two are one member order so that `kinds[i]` and `spans[i]`
    /// always describe the same member.
    kinds: Vec<String>,
    /// Per-member closed reference interval `[lo, hi]`, 0-based over the core,
    /// ascending. An insertion is empty (`lo = hi + 1`), so a separation
    /// computed as `next.lo - prev.hi - 1` is uniform across kinds.
    spans: Vec<(i64, i64)>,
    /// Gaps between consecutive members, from the spans above. Empty for a
    /// single-member output.
    separations: Vec<i64>,
    /// Which of the class's spellings reached this output.
    from: Vec<String>,
}

/// One divergent class.
#[derive(Serialize)]
struct Record {
    id: String,
    axis: String,
    /// The design's own parameters, so a shape can be cross-tabbed against them.
    design_members: usize,
    design_separation: usize,
    design_payload: usize,
    design_pattern: String,
    core: String,
    denoted: String,
    block_start: usize,
    reference_block: String,
    alt_block: String,
    /// The fine shape key: the sorted arity/kind signature of the distinct
    /// outputs. Useful for exemplars; too fine to adjudicate against (663
    /// distinct values over 4 643 classes).
    shape: String,
    /// The coarse shape key the adjudication is written against:
    /// `<relation>/<gap class>`. See [`Record::relation`] and
    /// [`Record::gap_class`].
    family: String,
    /// `typing` when every output partitions the reference identically and only
    /// the edit kinds differ; `partition` when the member boundaries differ.
    disagreement: String,
    /// How the outputs' arities relate: `spanning-vs-split` when one output is
    /// a single member and another is not; `split-vs-split` when they are all
    /// multi-member but of different arity; `same-arity` when the arities agree
    /// and only the boundaries or the kinds move — which includes the case
    /// where every output is a single spanning member.
    relation: String,
    /// Every distinct output arity, ascending.
    arities: Vec<usize>,
    /// The interior gaps of the **most split** output — the one the merge
    /// question is asked about. Every gap, not only the first: checking only
    /// the first under-counted a prior distribution by 54 per arm.
    split_gaps: Vec<i64>,
    /// Bucketed from [`Record::split_gaps`]: `all-0` (every interior gap is
    /// adjacency), `has-0` (some but not all), `min-1` (smallest gap is one
    /// nucleotide), `min-2+` (every gap is two or more). This is the quantity
    /// `general.md:34-35` and `DNA/delins.md:16-18` key on.
    gap_class: String,
    /// True when every output's arity equals the member count of every spelling
    /// that reached it — i.e. the normalizer preserved the input's partition
    /// rather than re-deriving one. This is the mechanism behind the divergence
    /// and it is worth stating separately from the shape.
    input_partition_preserving: bool,
    outputs: Vec<Output>,
}

// ---------------------------------------------------------------------------
// Providers — the same construction the generator and the axis test use
// ---------------------------------------------------------------------------

fn padded(core: &str) -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{core}{pad}")
}

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

fn axis_offset(axis: &str) -> i64 {
    if axis == "g" {
        PAD_OFFSET as i64
    } else {
        0
    }
}

// ---------------------------------------------------------------------------
// Shape
// ---------------------------------------------------------------------------

/// Classify one rendered member by its edit kind. Textual because the question
/// asked is which *notation* the normalizer chose, which is exactly what the
/// downstream consumer stores.
fn kind_of(member: &str) -> &'static str {
    let body = member.rsplit(':').next().unwrap_or(member);
    if body.contains("delins") {
        "delins"
    } else if body.contains('>') {
        "sub"
    } else if body.ends_with("inv") {
        "inv"
    } else if body.ends_with("dup") {
        "dup"
    } else if body.ends_with("del") {
        "del"
    } else if body.contains("ins") {
        "ins"
    } else if body.ends_with('=') {
        "identity"
    } else if body.contains('[') {
        "repeat"
    } else {
        "other"
    }
}

/// Split a normalized description into its cis members' rendered text, dropping
/// the accession and axis prefix.
fn member_texts(text: &str) -> Vec<String> {
    let body = match text.split_once(':') {
        Some((_, rest)) => rest,
        None => text,
    };
    let body = match body.split_once('.') {
        Some((_, rest)) => rest,
        None => body,
    };
    let body = body.trim();
    if let Some(inner) = body.strip_prefix('[').and_then(|b| b.strip_suffix(']')) {
        inner.split(';').map(|s| s.trim().to_string()).collect()
    } else {
        vec![body.to_string()]
    }
}

/// Each member's rendered text paired with its closed reference interval, from
/// `hgvs_to_spdi`, in **one** member order: ascending by interval. An insertion
/// is returned empty (`lo = hi + 1`) so a separation is uniform across kinds.
///
/// The single order is the point. `kinds` and `spans` used to be derived
/// separately — the former from [`member_texts`], which preserves the order the
/// members are *written* in, the latter sorted by coordinate — and then paired
/// positionally, so for any allele whose members are not already written in
/// coordinate order `kinds[i]` and `spans[i]` described different members. The
/// decided ruling's rationale enumerates adjacency **kind-pairs**, which is
/// exactly that join, so the mismatch reached a record rather than only a
/// fixture.
///
/// Returns `None` — never a partial answer — when any member's span cannot be
/// resolved. An empty `Vec` would flow on into `separations` and
/// [`classify_gaps`], which reports `no-gaps` for it, so a multi-member output
/// whose geometry was never measured would be filed as having no interior gaps.
fn members_of(
    provider: &MockProvider,
    text: &str,
    offset: i64,
) -> Option<Vec<(String, (i64, i64))>> {
    let members: Vec<HgvsVariant> = match parse_hgvs(text).ok()? {
        HgvsVariant::Allele(allele) => allele.variants.clone(),
        single => vec![single],
    };
    let texts = member_texts(text);
    if texts.len() != members.len() {
        return None;
    }
    let mut pairs = Vec::with_capacity(members.len());
    for (member, member_text) in members.iter().zip(texts) {
        let triple = hgvs_to_spdi(member, provider).ok()?;
        let lo = i64::try_from(triple.position).ok()? - offset;
        let hi = lo + triple.deletion.len() as i64 - 1;
        pairs.push((member_text, (lo, hi)));
    }
    pairs.sort_by_key(|(_, span)| *span);
    Some(pairs)
}

/// `next.lo - prev.hi - 1` over consecutive members. Uniform because an
/// insertion's interval is empty.
fn separations(spans: &[(i64, i64)]) -> Vec<i64> {
    spans.windows(2).map(|w| w[1].0 - w[0].1 - 1).collect()
}

/// Bucket a split output's interior gaps into the classes the separation rule
/// distinguishes. **Every** gap is considered — a codon or adjacency test
/// applied only to the first gap under-counts.
fn classify_gaps(gaps: &[i64]) -> String {
    if gaps.is_empty() {
        return "no-gaps".to_string();
    }
    let min = gaps.iter().copied().min().unwrap_or(0);
    if gaps.iter().all(|g| *g <= 0) {
        "all-0".to_string()
    } else if min <= 0 {
        "has-0".to_string()
    } else if min == 1 {
        "min-1".to_string()
    } else {
        "min-2+".to_string()
    }
}

/// How a class's distinct output arities relate. `arities` must be **sorted and
/// deduplicated**, and this is only asked of a class with two or more distinct
/// outputs.
///
/// `spanning-vs-split` asserts that one output is a single spanning member and
/// another is split, so it needs both a `1` and something else — hence
/// `contains(&1)` rather than `first() == Some(&1)`. The latter filed a class
/// whose outputs are *all* single-member (`arities == [1]`, two spellings of one
/// spanning member) as `spanning-vs-split`, asserting a split that is not there,
/// and made `same-arity` unreachable at arity 1. `spanning-vs-split/min-2+` is
/// the headline family and the module's central mechanism claim is stated over
/// it, so a class that does not belong to the claim must not be counted in it.
fn relation_of(arities: &[usize]) -> &'static str {
    if arities.len() > 1 && arities.contains(&1) {
        "spanning-vs-split"
    } else if arities.len() > 1 {
        "split-vs-split"
    } else {
        "same-arity"
    }
}

fn signature(arity: usize, kinds: &[String]) -> String {
    if arity == 1 {
        format!("1x{}", kinds[0])
    } else {
        format!("{arity}x[{}]", kinds.join(","))
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

fn main() -> ExitCode {
    let cli = Cli::parse();
    let direction = match cli.direction.as_str() {
        "3prime" => ShuffleDirection::ThreePrime,
        "5prime" => ShuffleDirection::FivePrime,
        other => {
            eprintln!("error: --direction must be 3prime or 5prime, got {other}");
            return ExitCode::FAILURE;
        }
    };

    let text = match std::fs::read_to_string(&cli.corpus) {
        Ok(text) => text,
        Err(e) => {
            eprintln!(
                "error: reading {}: {e}\nGenerate it with:\n  cargo run --features dev \
                 --example generate_cis_confluence_corpus",
                cli.corpus.display()
            );
            return ExitCode::FAILURE;
        }
    };
    let corpus: Corpus = match serde_json::from_str(&text) {
        Ok(corpus) => corpus,
        Err(e) => {
            eprintln!("error: parsing {}: {e}", cli.corpus.display());
            return ExitCode::FAILURE;
        }
    };

    // Zeroed immediately before `collect`, so the ratio below describes this
    // corpus rather than whatever normalization the set-up above happened to do.
    ferro_hgvs::normalize::dev_partitioners::reset_preserve_census();
    let (records, cross, drops) = collect(&corpus, direction);
    print!(
        "{}",
        render_census(&corpus, &records, &cross, &drops, cli.top, cli.examples)
    );

    // The partition-preserving arm's own hit rate, printed unconditionally
    // because it is the denominator for every confluence figure above: each
    // decline falls back to `partition_block`, i.e. to a re-derived partition,
    // so a non-zero decline rate means the numbers describe a BLEND of two
    // partitioners rather than the partition model. That distinction was
    // reported wrong once, from an uninstrumented run, at 0% when it was 33.4%.
    //
    // `kept + declined` is the number of blocks that reached the arm, which is
    // not the class count — one class has many spellings and some never reach
    // canonicalization at all.
    //
    // A zero here cannot be a *disabled instrument*, which is the one other
    // thing it could plausibly mean and the reading the census doc on
    // `preserve_census` exists to rule out. `record_preserve_outcome` is a no-op
    // without the `dev` feature, but this example cannot be built without it:
    // `dev_partitioners` is `#[cfg(feature = "dev")]`, so the two calls below do
    // not exist otherwise, and `Cargo.toml` marks this target
    // `required-features = ["dev"]`, so cargo declines the target by name rather
    // than building it against a stubbed counter. There is no configuration in
    // which this branch reports a silenced instrument as an empty corpus.
    let (kept, declined) = ferro_hgvs::normalize::dev_partitioners::preserve_census();
    let total = kept + declined;
    if total == 0 {
        println!(
            "\npartition-preserving arm: NOT REACHED (0 invocations) — \
             this corpus measures nothing about the partition rule"
        );
    } else {
        println!(
            "\npartition-preserving arm: {kept} kept, {declined} declined of {total} \
             ({:.1}% declined, each falling back to a re-derived partition)",
            100.0 * declined as f64 / total as f64,
        );
    }

    if cli.stats {
        return ExitCode::SUCCESS;
    }
    let Some(path) = cli.out else {
        return ExitCode::SUCCESS;
    };
    let mut rendered = match serde_json::to_string_pretty(&records) {
        Ok(text) => text,
        Err(e) => {
            eprintln!("error: serializing: {e}");
            return ExitCode::FAILURE;
        }
    };
    rendered.push('\n');
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        if let Err(e) = std::fs::create_dir_all(parent) {
            eprintln!("error: creating {}: {e}", parent.display());
            return ExitCode::FAILURE;
        }
    }
    if let Err(e) = std::fs::write(&path, rendered) {
        eprintln!("error: writing {}: {e}", path.display());
        return ExitCode::FAILURE;
    }
    println!(
        "wrote {} divergent classes to {}",
        records.len(),
        path.display()
    );
    ExitCode::SUCCESS
}

fn collect(corpus: &Corpus, direction: ShuffleDirection) -> (Vec<Record>, CrossClass, Drops) {
    let classes = &corpus.classes;
    let mut records = Vec::new();
    let mut drops = Drops::default();
    // Keyed by `(axis, reference, denoted)` — one entry per *variant*, holding
    // every output any of its classes reached and whether each class converged.
    let mut variants: BTreeMap<VariantKey, VariantOutputs> = BTreeMap::new();
    // Classes sharing a reference are contiguous by id, so one provider per
    // group rather than per class — the same optimization the axis test makes.
    let mut start = 0usize;
    while start < classes.len() {
        let mut end = start + 1;
        while end < classes.len()
            && classes[end].axis == classes[start].axis
            && classes[end].core == classes[start].core
        {
            end += 1;
        }
        let (provider, _reference) = reference_for(&classes[start]);
        let normalizer = Normalizer::with_config(
            provider.clone(),
            NormalizeConfig::default().with_direction(direction),
        );
        let offset = axis_offset(&classes[start].axis);

        for class in &classes[start..end] {
            let mut grouped: BTreeMap<String, Vec<String>> = BTreeMap::new();
            for spelling in &class.spellings {
                let Ok(variant) = parse_hgvs(spelling) else {
                    drops.unparsed_spellings += 1;
                    continue;
                };
                let normalized = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    normalizer.normalize(&variant)
                }));
                let Ok(Ok(value)) = normalized else {
                    drops.unnormalized_spellings += 1;
                    continue;
                };
                grouped
                    .entry(value.to_string())
                    .or_default()
                    .push(spelling.clone());
            }
            let entry = variants
                .entry((
                    class.axis.clone(),
                    class.core.clone(),
                    class.denoted.clone(),
                ))
                .or_insert_with(|| (BTreeSet::new(), true, 0));
            entry.0.extend(grouped.keys().cloned());
            entry.1 &= grouped.len() <= 1;
            entry.2 += 1;

            if grouped.len() < 2 {
                continue;
            }

            let mut outputs = Vec::with_capacity(grouped.len());
            let mut unresolved = false;
            for (text, from) in grouped {
                let Some(members) = members_of(&provider, &text, offset) else {
                    unresolved = true;
                    break;
                };
                let kinds: Vec<String> = members
                    .iter()
                    .map(|(member, _)| kind_of(member).to_string())
                    .collect();
                let spans: Vec<(i64, i64)> = members.iter().map(|(_, span)| *span).collect();
                let separations = separations(&spans);
                outputs.push(Output {
                    arity: kinds.len(),
                    kinds,
                    spans,
                    separations,
                    from,
                    text,
                });
            }
            // Excluded, not classified. Every downstream quantity — the gap
            // class, the partition comparison behind `disagreement`, the
            // adjacency kind-pairs — is read off the geometry, so a class one
            // of whose outputs could not be measured has no honest shape.
            if unresolved {
                drops.unresolved_classes += 1;
                continue;
            }

            let mut signatures: Vec<String> = outputs
                .iter()
                .map(|o| signature(o.arity, &o.kinds))
                .collect();
            signatures.sort();
            let shape = signatures.join(" | ");

            let partitions: BTreeSet<Vec<(i64, i64)>> =
                outputs.iter().map(|o| o.spans.clone()).collect();
            let disagreement = if partitions.len() == 1 {
                "typing".to_string()
            } else {
                "partition".to_string()
            };

            let mut arities: Vec<usize> = outputs.iter().map(|o| o.arity).collect();
            arities.sort_unstable();
            arities.dedup();
            let relation = relation_of(&arities).to_string();

            let split_gaps = outputs
                .iter()
                .max_by_key(|o| o.arity)
                .map(|o| o.separations.clone())
                .unwrap_or_default();
            let gap_class = classify_gaps(&split_gaps);
            let family = format!("{relation}/{gap_class}");
            let input_partition_preserving = outputs
                .iter()
                .all(|o| o.from.iter().all(|s| member_texts(s).len() == o.arity));

            records.push(Record {
                id: class.id.clone(),
                axis: class.axis.clone(),
                design_members: class.members,
                design_separation: class.separation,
                design_payload: class.payload,
                design_pattern: class.pattern.clone(),
                core: class.core.clone(),
                denoted: class.denoted.clone(),
                block_start: class.block_start,
                reference_block: class.reference_block.clone(),
                alt_block: class.alt_block.clone(),
                shape,
                family,
                disagreement,
                relation,
                arities,
                split_gaps,
                gap_class,
                input_partition_preserving,
                outputs,
            });
        }
        start = end;
    }

    let multi: Vec<_> = variants.values().filter(|v| v.2 > 1).collect();
    let cross = CrossClass {
        variants: variants.len(),
        multi_class: multi.len(),
        multi_class_divergent: multi.iter().filter(|v| v.0.len() > 1).count(),
        hidden_by_the_class_boundary: multi.iter().filter(|v| v.0.len() > 1 && v.1).count(),
    };
    (records, cross, drops)
}

fn render_census(
    corpus: &Corpus,
    records: &[Record],
    cross: &CrossClass,
    drops: &Drops,
    top: usize,
    examples: usize,
) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "classes: {}  divergent: {}",
        corpus.classes.len(),
        records.len()
    );
    let _ = writeln!(
        out,
        "  dropped (expected 0 everywhere): {} spellings unparsed, {} spellings unnormalized, \
         {} divergent classes excluded for unresolved member spans",
        drops.unparsed_spellings, drops.unnormalized_spellings, drops.unresolved_classes
    );

    let _ = writeln!(
        out,
        "  cross-class: {} distinct (axis, reference, denoted) variants; {} reachable from more \
         than one class, of which {} do not agree on one output — {} of those have every class \
         converging internally, so the class boundary is the only thing hiding them",
        cross.variants,
        cross.multi_class,
        cross.multi_class_divergent,
        cross.hidden_by_the_class_boundary
    );

    let mut by_disagreement: BTreeMap<&str, usize> = BTreeMap::new();
    let mut by_gap_class: BTreeMap<&str, usize> = BTreeMap::new();
    let mut by_axis: BTreeMap<&str, usize> = BTreeMap::new();
    let mut by_design_sep: BTreeMap<usize, usize> = BTreeMap::new();
    let mut by_shape: BTreeMap<&str, usize> = BTreeMap::new();
    let mut by_family: BTreeMap<&str, usize> = BTreeMap::new();
    let mut by_family_axis: BTreeMap<String, usize> = BTreeMap::new();
    let mut preserving = 0usize;
    for record in records {
        *by_disagreement
            .entry(record.disagreement.as_str())
            .or_default() += 1;
        *by_gap_class.entry(record.gap_class.as_str()).or_default() += 1;
        *by_axis.entry(record.axis.as_str()).or_default() += 1;
        *by_design_sep.entry(record.design_separation).or_default() += 1;
        *by_shape.entry(record.shape.as_str()).or_default() += 1;
        *by_family.entry(record.family.as_str()).or_default() += 1;
        *by_family_axis
            .entry(format!("{} [{}]", record.family, record.axis))
            .or_default() += 1;
        if record.input_partition_preserving {
            preserving += 1;
        }
    }
    let _ = writeln!(
        out,
        "  input-partition-preserving: {preserving} ({:.1}%)",
        100.0 * preserving as f64 / records.len().max(1) as f64
    );

    let section = |title: &str, rows: Vec<(String, usize)>, out: &mut String| {
        let _ = writeln!(out, "  {title}:");
        for (key, count) in rows {
            let _ = writeln!(
                out,
                "    {key}: {count} ({:.1}%)",
                100.0 * count as f64 / records.len().max(1) as f64
            );
        }
    };
    section(
        "by disagreement",
        by_disagreement
            .iter()
            .map(|(k, v)| ((*k).to_string(), *v))
            .collect(),
        &mut out,
    );
    section(
        "by axis",
        by_axis
            .iter()
            .map(|(k, v)| ((*k).to_string(), *v))
            .collect(),
        &mut out,
    );
    section(
        "by gap class of the most-split output",
        by_gap_class
            .iter()
            .map(|(k, v)| ((*k).to_string(), *v))
            .collect(),
        &mut out,
    );
    section(
        "by design separation",
        by_design_sep
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect(),
        &mut out,
    );
    let mut families: Vec<(&str, usize)> = by_family.into_iter().collect();
    families.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
    let mut cumulative = 0usize;
    let _ = writeln!(
        out,
        "  FAMILIES (the adjudication unit): {}",
        families.len()
    );
    for (family, count) in &families {
        cumulative += count;
        let _ = writeln!(
            out,
            "    {count:6} ({:5.1}%, cum {:5.1}%)  {family}",
            100.0 * *count as f64 / records.len().max(1) as f64,
            100.0 * cumulative as f64 / records.len().max(1) as f64,
        );
    }
    section(
        "family x axis",
        by_family_axis.into_iter().collect(),
        &mut out,
    );

    let mut shapes: Vec<(&str, usize)> = by_shape.into_iter().collect();
    shapes.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
    let _ = writeln!(out, "  shapes: {} distinct", shapes.len());
    let mut cumulative = 0usize;
    for (shape, count) in shapes.iter().take(top) {
        cumulative += count;
        let _ = writeln!(
            out,
            "    {count:6} ({:5.1}%, cum {:5.1}%)  {shape}",
            100.0 * *count as f64 / records.len().max(1) as f64,
            100.0 * cumulative as f64 / records.len().max(1) as f64,
        );
        for record in records.iter().filter(|r| r.shape == *shape).take(examples) {
            let _ = writeln!(out, "        {} [{}]", record.id, record.disagreement);
            for output in &record.outputs {
                let _ = writeln!(out, "          {}   <- {:?}", output.text, output.from);
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The separation arithmetic is the one thing here that has already
    /// produced published numbers that were wrong, so it is pinned directly:
    /// an insertion consumes no reference position, so its interval is empty
    /// and `next.lo - prev.hi - 1` stays uniform across kinds.
    #[test]
    fn an_insertion_consumes_no_reference_position() {
        // A deletion of [10, 12] followed by an insertion at interbase 15,
        // i.e. between reference bases 14 and 15: the empty interval is
        // (15, 14), so the gap is 15 - 12 - 1 = 2.
        let spans = [(10, 12), (15, 14)];
        assert_eq!(separations(&spans), vec![2]);
        // Two adjacent members leave a gap of zero.
        assert_eq!(separations(&[(10, 12), (13, 15)]), vec![0]);
        // An insertion flush against the end of a deletion: interbase 13, so
        // the empty interval is (13, 12) and the gap is 13 - 12 - 1 = 0.
        assert_eq!(separations(&[(10, 12), (13, 12)]), vec![0]);
        // Every interior gap is reported, not only the first.
        assert_eq!(separations(&[(1, 1), (3, 3), (9, 9)]), vec![1, 5]);
    }

    #[test]
    fn members_and_kinds_are_read_off_the_rendered_form() {
        assert_eq!(
            member_texts("NM_TEST.1:c.[10_12delinsAA;15dup]"),
            vec!["10_12delinsAA".to_string(), "15dup".to_string()]
        );
        assert_eq!(member_texts("NC_TEST.1:g.266_268del"), vec!["266_268del"]);
        assert_eq!(kind_of("10_12delinsAA"), "delins");
        assert_eq!(kind_of("266_268del"), "del");
        assert_eq!(kind_of("266_267insAA"), "ins");
        assert_eq!(kind_of("266dup"), "dup");
        assert_eq!(kind_of("266A>G"), "sub");
        assert_eq!(kind_of("266_270inv"), "inv");
    }

    /// `spanning-vs-split` names the headline family, over which the module's
    /// central mechanism claim is stated, so it must not absorb classes that do
    /// not have a split at all.
    #[test]
    fn an_all_spanning_class_is_not_filed_as_a_spanning_vs_split() {
        // The regression: `arities == [1]` is two spellings reaching two
        // *different* single spanning members. Nothing is split.
        assert_eq!(relation_of(&[1]), "same-arity");
        assert_eq!(relation_of(&[3]), "same-arity");
        // A `1` present alongside anything else is the real shape.
        assert_eq!(relation_of(&[1, 2]), "spanning-vs-split");
        assert_eq!(relation_of(&[1, 2, 4]), "spanning-vs-split");
        // Differing arities with no lone spanning member.
        assert_eq!(relation_of(&[2, 3]), "split-vs-split");
    }

    /// An output whose geometry could not be measured must not reach
    /// [`classify_gaps`]: an empty span list is indistinguishable there from a
    /// genuinely single-member output, which is how an unmeasured multi-member
    /// output used to leave the `has-0` population.
    #[test]
    fn no_gaps_is_what_an_unmeasured_output_would_have_looked_like() {
        assert_eq!(classify_gaps(&[]), "no-gaps");
        assert_eq!(classify_gaps(&[0, 4]), "has-0");
        assert_eq!(classify_gaps(&[3, 4]), "min-2+");
        assert_eq!(classify_gaps(&[1, 4]), "min-1");
    }

    #[test]
    fn a_signature_names_the_arity_and_the_kinds_in_order() {
        assert_eq!(signature(1, &["delins".to_string()]), "1xdelins");
        assert_eq!(
            signature(2, &["del".to_string(), "ins".to_string()]),
            "2x[del,ins]"
        );
    }
}
