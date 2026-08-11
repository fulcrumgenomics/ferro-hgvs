//! Dump normalized output over a deterministic synthetic corpus, and diff two
//! dumps — the measurement a representation-moving PR has to carry (#1441).
//!
//! Representation stability is a shipped guarantee: the downstream consumer keys
//! read counts on the normalized HGVS string, so a PR that moves an output has to
//! say which forms moved, in which direction, over how many rows, and whether the
//! affected inputs were previously **rejected** (no stored string, so free) or
//! previously **accepted** (a migration). Five PRs made that measurement by hand
//! and each rebuilt the harness; this is that harness, committed.
//!
//! ## Usage
//!
//! ```text
//! # at the base revision
//! cargo run --features dev --example dump_normalized_corpus -- --out /tmp/before.tsv
//! # at the candidate revision
//! cargo run --features dev --example dump_normalized_corpus -- --out /tmp/after.tsv
//! # the table for the PR body
//! cargo run --features dev --example dump_normalized_corpus -- \
//!     --compare /tmp/before.tsv --against /tmp/after.tsv
//! ```
//!
//! ## Comparing across a revision that changed the corpus
//!
//! `compare` refuses two dumps that do not cover the same rows, so a revision that
//! adds a shape family cannot be diffed against a dump taken with the older
//! generator. Copy **this file** into the base checkout and dump from there: it
//! depends only on the public API, so it builds against either revision and both
//! dumps then describe the same corpus.
//!
//! ```text
//! git worktree add --detach ../base <base-sha>
//! cp examples/dump_normalized_corpus.rs ../base/examples/
//! (cd ../base && cargo run --features dev --example dump_normalized_corpus -- --out /tmp/before.tsv)
//! cargo run --features dev --example dump_normalized_corpus -- --out /tmp/after.tsv
//! ```
//!
//! ## What this is NOT for
//!
//! Prepared references. `ferro normalize --input <file> --reference <dir> --format
//! tsv` already dumps a real corpus against real reference data, and duplicating it
//! here would mean two things to keep in step. What that path *cannot* do is drive a
//! **synthetic** reference — `--reference` requires a `manifest.json` and the
//! JSON-reference constructor is Python-only — and synthetic references are what
//! densely cover the cis-allele shape families where representation churn actually
//! happens. That is the gap this fills.
//!
//! ## Three traps, the first two hit while measuring #1401
//!
//! 1. **A row's identity is `(reference, axis, direction, input)` — never the input
//!    alone.** Every synthetic reference reuses one accession, so keying on the
//!    input string silently merges rows drawn against different sequences: 2,160
//!    rows collapsed to 135 and the reported movement came out at 6 instead of 498.
//!    Nothing failed. `the_corpus_key_is_unique` pins this.
//! 2. **A within-shape-family rate is not a library-wide rate.** This corpus is
//!    deliberately enriched for the shapes that churn, so its percentages are far
//!    higher than a real corpus would give — #1401 measures 23% here. Quote it as
//!    "of the affected shape family", never as a repo-wide figure.
//! 3. **A zero measures only the families the corpus builds** (#1456). Before the
//!    conflict families below were added, this harness reported `0 of 18,432` for
//!    #1403, #1445 and #1451 in a row — none of which the corpus could express, and
//!    two of which move behaviour for thousands of inputs by their own enumerations.
//!    `compare` now says so in the zero case rather than leaving it to be inferred.
//!    If your change is scoped to a shape no family builds, add the family; do not
//!    quote the zero.
//! 4. **…and only at the scale it builds them** (#1460). The families are drawn
//!    against 20-mers, so a change gated on *length* rather than on shape sees the
//!    same path on both sides no matter how many families exist. #1403 caps a
//!    partition at `MAX_SPLIT_BLOCK` (1024) and measured a guaranteed zero here
//!    even after the conflict families landed. `long_corpus_sequences` is the
//!    answer, and it is deliberately narrow: two cores, one family, 16 rows. Adding
//!    scale is expensive in a way adding a family is not — check the dump cost
//!    before crossing a kilobase core with anything.
//!
//! 5. **…and only at the sequence structure it builds them from** (`#1517`). The
//!    families are drawn *on* the random cores, so any property of the reference
//!    itself is whatever the xorshift produced. Where a reverse-complement block
//!    coincides with its own reference is exactly such a property, and the one
//!    family that built such a block —`delins_hiding_an_inversion` — emits its
//!    pieces adjacent, so nothing in the corpus separated two multi-column runs of
//!    an inversion and a gate keyed on that measured `0 of 78,028`.
//!    `separated_revcomp_runs` is the answer, and like `long_block_inversion` it
//!    brings its own designed references because a random core cannot be asked to
//!    have a coincidence pattern.
//!
//! ## One family knows its own ground truth, and that buys three oracles
//!
//! Everything above measures **movement** — two dumps, one diff — and never
//! whether an output is *right*. Thirteen of the fifteen families cannot be asked:
//! they are sets of descriptions with no record of what the description was meant
//! to denote. `separated_revcomp_runs` is built the other way round, choosing the
//! alternate first and deriving the reference around it, so `(reference,
//! alternate)` is exact for every row. The test module turns that into three real
//! oracles — sequence preservation, confluence, and separation-rule conformance —
//! and two of the three are currently red against `main`. Read their doc comments
//! before quoting them; each records what it measured and which half of its
//! finding is a settled defect and which is an open decision.
//!
//! ## The corpus is deliberately independent of `tests/it/common/`
//!
//! An example cannot reach test helpers, so the reference construction below is its
//! own. That is a feature rather than a compromise: this corpus has to stay
//! **stable across revisions** so that two dumps are comparable, and it must not
//! move because someone refactored `SyntheticBuilder`. It does not need to agree
//! with the test helpers — only with itself, over time. Changing the generator or
//! the seed count re-rolls the corpus and makes old dumps incomparable; treat it
//! the way `sweep_sequences` is treated.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;

use ferro_hgvs::normalize::{NormalizeConfig, ShuffleDirection};
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::reference::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};

/// 256 bases of period-4 `ACGT`, so the base immediately 5' of a core is `T` and
/// the one immediately 3' of it is `A`. A core starting with something other than
/// `A` and ending with something other than `T` therefore cannot extend the pad's
/// own rotation, which bounds a repeat tract to the core.
const PAD_OFFSET: usize = 256;

/// The genomic contig a `g.` row is drawn against.
const GENOMIC_CONTIG: &str = "NC_TEST.1";
/// The transcript a `c.` row is drawn against.
const TX_ACCESSION: &str = "NM_TEST.1";
/// The contig the synthetic transcript is placed on.
const TX_CONTIG: &str = "chr_synth";

/// The CDS of a synthetic coding reference: 1-based inclusive.
///
/// A multiple of three, so the CDS is codon-complete, and deliberately short of the
/// 20-base core so that **five** 3'UTR bases (`c.*1`..`c.*5`) remain. Several shapes
/// then straddle `CDS_END`, which is the region #1398, #1418 and both halves of
/// #1426 lived in — a corpus that stops at the junction reports a confident zero for
/// every one of them.
const CDS_START: u64 = 1;
const CDS_END: u64 = 15;

/// The multi-exon transcript a `cx.` row is drawn against (#1478).
///
/// A *second* coding reference rather than a change to the one above. Re-pointing
/// `CDS_START`/`CDS_END` or the exon layout in place would re-roll every existing
/// `c.` row and make prior dumps incomparable for families that have nothing to do
/// with splicing; an added axis leaves them byte-identical, exactly as an added
/// family does.
const TX_MULTI_ACCESSION: &str = "NM_TESTX.1";
/// The contig the multi-exon synthetic transcript is placed on.
const TX_MULTI_CONTIG: &str = "chr_synth_x";

/// CDS of the multi-exon reference, 1-based inclusive over transcript coordinates.
///
/// `CDS_START` is **4**, not 1, which is the whole point: with a CDS starting at the
/// first base there is no 5'UTR, so `c.-N` is unreachable and every 5'UTR shape
/// measures as a confident zero. Three bases before it give `c.-3`..`c.-1`. The CDS
/// is 12 bases (4..=15), still a multiple of three so it stays codon-complete, and
/// five 3'UTR bases remain as on the single-exon reference.
const CDS_START_MULTI: u64 = 4;
const CDS_END_MULTI: u64 = 15;

/// Exon blocks of the multi-exon transcript, as 0-based half-open ranges into the
/// 20-base core.
///
/// Three exons, so there are **two** exon/exon junctions — at transcript 7/8 and
/// 14/15. `general.md:44` exempts deletions and duplications around such a junction
/// from the 3' rule, and #1450 is that the per-member pipeline honours the exemption
/// while the sequence-first cis path does not. A single-exon reference has no
/// junction at all, which is why that defect cannot appear in this corpus today.
const EXON_SPANS: &[(usize, usize)] = &[(0, 7), (7, 14), (14, 20)];

/// Intronic bases between consecutive exons in the genomic layout.
///
/// Long enough that an exon's genomic span cannot abut its neighbour's, so a shift
/// that runs off an exon end lands in intron rather than silently in the next exon.
const INTRON_LEN: usize = 10;

/// The dump header, written on emit and required verbatim on read.
const HEADER: &str = "reference\taxis\tdirection\tfamily\tinput\toutput\twas_fixed_point\n";

/// Column count of a dump row, checked when reading a dump back.
const COLUMNS: usize = 7;

#[derive(Parser, Debug)]
#[command(about = "Dump normalized output over a synthetic corpus, or diff two dumps (#1441)")]
struct Cli {
    /// Write the dump here (default: stdout). Dump mode only.
    #[arg(long, conflicts_with_all = ["compare", "against"])]
    out: Option<PathBuf>,
    /// Seed count for the reference corpus. Two sequences per seed (an `AT` and an
    /// `ACGT` alphabet). Prefix-stable: a smaller count is a strict prefix.
    /// Dump mode only — in diff mode the corpus comes from the dumps being compared,
    /// so accepting it here would silently ignore it and invite the belief that a
    /// comparison had been re-scoped.
    #[arg(long, default_value_t = 24, conflicts_with_all = ["compare", "against"])]
    seeds: u32,
    /// Diff mode: the baseline dump.
    #[arg(long, requires = "against")]
    compare: Option<PathBuf>,
    /// Diff mode: the candidate dump.
    #[arg(long, requires = "compare")]
    against: Option<PathBuf>,
    /// Also check, per row, whether the normalized output denotes the **same
    /// bases** as its input, and report the rows where it does not (#1514).
    ///
    /// Dump mode only. Findings go to stderr and the dump itself is unchanged,
    /// so a verified run stays byte-comparable with every existing baseline.
    #[arg(long, conflicts_with_all = ["compare", "against"])]
    verify_spdi: bool,
}

/// What [`verify_row`] concluded about one row.
///
/// Three outcomes, not two, and the third is the one that matters for reading
/// the summary: `canonical_spdi` declines an allele whose members overlap, and
/// on this corpus that is thousands of rows. Folding those into "different"
/// buried 45 real findings under ~6,700 rows of noise the first time this
/// check was written by hand.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SpdiVerdict {
    /// Input and output denote the same bases.
    Same,
    /// They denote different bases — the finding this flag exists for.
    Different,
    /// One side could not be applied, so the question has no answer here.
    /// Not a failure.
    Unverifiable,
}

/// Whether `row`'s output denotes the same bases as its input.
///
/// Compared through [`Normalizer::canonical_spdi`], which derives its key from
/// the bases a description *produces* rather than from how it is written, and
/// maximally 3'-shifts it. Two spellings of one edit therefore key identically
/// by construction, so a difference here is a difference in meaning and not in
/// notation — which is exactly what a row-movement count cannot tell you.
///
/// Wrapped in `catch_unwind` for the same reason [`normalize_one`] is: this runs
/// over every row of a corpus built to contain the shapes that break things, and
/// one panicking row must not take the whole dump with it.
///
/// The provider is built from [`Row::core`] and **not** from `Row::reference`
/// (#1624). Those two are the same string for thirteen of the fifteen families,
/// which is what made reading the wrong one survive review: the two families
/// whose `reference` column is a *label* — `long_block_inversion` and
/// `separated_revcomp_runs` — were verified against the label itself. On the
/// `cx` axis that panicked, and on `g.`/`c.` it did something worse, reporting
/// confident `SPDI-MISMATCH` findings whose keys were read out of the letters of
/// `revcomp_sep0_at0`.
fn verify_row(row: &Row) -> (SpdiVerdict, String, String) {
    let declined = |what: &str| format!("<{what}>");
    let (Ok(input), Ok(output)) = (parse_hgvs(&row.input), parse_hgvs(&row.output)) else {
        return (
            SpdiVerdict::Unverifiable,
            declined("unparseable"),
            declined("unparseable"),
        );
    };
    let direction = match row.direction {
        "5prime" => ShuffleDirection::FivePrime,
        _ => ShuffleDirection::ThreePrime,
    };
    let normalizer = Normalizer::with_config(
        provider_for(row.axis, &row.core),
        NormalizeConfig::default().with_direction(direction),
    );
    let key = |v: &_| match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        normalizer.canonical_spdi(v)
    })) {
        Ok(Ok(k)) => Some(k.to_string()),
        Ok(Err(_)) => None,
        Err(_) => None,
    };
    match (key(&input), key(&output)) {
        (Some(a), Some(b)) if a == b => (SpdiVerdict::Same, a, b),
        (Some(a), Some(b)) => (SpdiVerdict::Different, a, b),
        (a, b) => (
            SpdiVerdict::Unverifiable,
            a.unwrap_or_else(|| declined("declined")),
            b.unwrap_or_else(|| declined("declined")),
        ),
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    // This is the harness a bake-off runs through, so it is the entry point that
    // can least afford to serve `live` under a candidate's name. In a release
    // build the library falls safe rather than aborting, and the two silences
    // compound: the dump would come back identical to the `live` dump, and
    // `report_partition_declines` would print nothing at all, because a `live`
    // run attempts no sequence-first partition. Refuse here instead.
    if let Some(message) = ferro_hgvs::normalize::partition_switch_startup_error() {
        eprintln!("error: {message}");
        return ExitCode::FAILURE;
    }
    if let (Some(before), Some(after)) = (&cli.compare, &cli.against) {
        return match compare(before, after) {
            Ok(report) => {
                print!("{report}");
                ExitCode::SUCCESS
            }
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::FAILURE
            }
        };
    }

    let rows = dump(cli.seeds);
    report_partition_declines();
    let mut out = String::new();
    out.push_str(HEADER);
    for row in &rows {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.reference,
            row.axis,
            row.direction,
            row.family,
            row.input,
            row.output,
            row.was_fixed_point
        );
    }
    if cli.verify_spdi {
        let (mut same, mut different, mut unverifiable) = (0usize, 0usize, 0usize);
        for row in &rows {
            let (verdict, key_in, key_out) = verify_row(row);
            match verdict {
                SpdiVerdict::Same => same += 1,
                SpdiVerdict::Unverifiable => unverifiable += 1,
                SpdiVerdict::Different => {
                    different += 1;
                    // Stable, greppable, and on stderr: the recipe for "did this
                    // change make anything newly wrong" is to run this on both
                    // revisions and diff the two reports. Putting it in the dump
                    // instead would add a column, and the reader exact-matches
                    // the header — every committed baseline would stop loading.
                    eprintln!(
                        "SPDI-MISMATCH\t{}\t{} {}\t{}\n    in : {}  => {}\n    out: {}  => {}",
                        row.family,
                        row.axis,
                        row.direction,
                        row.reference,
                        row.input,
                        key_in,
                        row.output,
                        key_out
                    );
                }
            }
        }
        eprintln!(
            "verified {} rows: {different} denote different bases, {same} agree, \
             {unverifiable} unverifiable (member overlap or unparseable)",
            rows.len()
        );
    }

    match &cli.out {
        Some(path) => {
            if let Err(e) = fs::write(path, &out) {
                eprintln!("error: writing {}: {e}", path.display());
                return ExitCode::FAILURE;
            }
            eprintln!("wrote {} rows to {}", rows.len(), path.display());
        }
        None => print!("{out}"),
    }
    ExitCode::SUCCESS
}

/// Say, on stderr, how much of this dump a sequence-first arm actually answered.
///
/// A `FERRO_PARTITION=canonical` dump that comes back identical to the `live`
/// dump reads as *"the candidate changes nothing"*. It can equally mean the
/// candidate declined every block and `partition_block` produced both columns,
/// and no output distinguishes those — which is what makes the line below part
/// of the measurement rather than decoration.
///
/// Printed **unconditionally when the arm ran at all**, including the zero:
/// `0 declined of N` is the reading that licenses quoting the diff, and it is
/// only worth anything if the same line would have shown a non-zero.
fn report_partition_declines() {
    let counts = ferro_hgvs::normalize::partition_decline_counts();
    if counts.attempted == 0 {
        // The selected arm was `live`, which has no decline path, so there is no
        // census to report and a "0 of 0" line would read as a clean bill of
        // health for an arm that did not run.
        //
        // "Selected" is the load-bearing word: the separate
        // `FERRO_SEQFIRST_SHADOW` audit asks a sequence-first partitioner on
        // every block whatever `FERRO_PARTITION` says, and those attempts are
        // deliberately outside this census — they cannot reach the emitted
        // pieces, and the audit reports its own outcomes through `log::debug!`.
        return;
    }
    eprintln!(
        "partitioner: {} of {} blocks declined and were served by `partition_block` \
         (the shipped rule) instead",
        counts.declined, counts.attempted
    );
}

/// One measured row. `was_fixed_point` records whether the *input* was already its
/// own normalized form on the revision that produced this dump — which is what
/// makes the cheap/expensive split derivable from the baseline dump alone, without
/// re-running the old code.
struct Row {
    /// The dump's `reference` column, and part of the row's identity. For
    /// thirteen of the fifteen families this **is** the reference sequence; for
    /// the two that bring their own designed references it is a label, because a
    /// kilobase core repeated on every row is not a column anyone can read.
    reference: String,
    /// The reference **sequence** the row was drawn against.
    ///
    /// Equal to `reference` except on `long_block_inversion` and
    /// `separated_revcomp_runs`, where that column is a label. Deliberately
    /// **not** emitted: the dump's format is exact-matched on read, so adding a
    /// column would stop every committed baseline from loading, and a row's
    /// identity does not need the sequence when the label already determines it.
    ///
    /// It exists because a consumer that needs to *re-derive* the row — which is
    /// what `--verify-spdi` does — cannot get the sequence back out of a label,
    /// and reading the label as one is not an error either provider reports
    /// (#1624).
    core: String,
    axis: &'static str,
    direction: &'static str,
    family: &'static str,
    input: String,
    output: String,
    was_fixed_point: bool,
}

/// The shape families, each drawn from a defect this repository has actually had.
/// Adding one is fine; renaming or reordering re-rolls nothing, but changing what a
/// family *emits* makes old dumps incomparable for that family.
const FAMILIES: &[(&str, &str)] = &[
    (
        "split_vs_spanning_delins",
        "#1232 / #1401 — two delins separated by one unchanged base",
    ),
    (
        "dup_plus_sub",
        "#999 — the split buys a higher-priority type, so it must stay split",
    ),
    (
        "adjacent_junction_ins",
        "#1235 — two insertions at neighbouring junctions",
    ),
    (
        "dup_plus_ins",
        "#1320 / #1323 — a dup and an insertion sharing a junction",
    ),
    (
        "del_plus_sub",
        "#1232's own example — a deletion and a separated substitution",
    ),
    // The four below are the *conflicting* shapes (#1456). Every family above
    // pairs members that occupy disjoint territory, which is why this corpus
    // reported `0 of 18,432` for #1403, #1445 and #1451 in a row — three changes
    // that move behaviour for thousands of inputs by their own enumerations. The
    // two detectors in `normalize::overlap` between them define the shapes: one
    // keys on coincident (and, since #1451, merely intersecting) spans, the other
    // on a junction anchored at or inside a span.
    (
        "coincident_bounds",
        "#1307 — a dup and a substitution on one and the same base",
    ),
    (
        "junction_interior_to_span",
        "#486 — an insertion at a junction interior to a span edit",
    ),
    (
        "partial_overlap_spans",
        "#1451 — two spans that intersect without either containing the other",
    ),
    (
        "nested_spans",
        "#1451 — two spans sharing a bound, one containing the other",
    ),
    (
        "delins_hiding_an_inversion",
        "#1454 — a spanning delins whose trailing piece is an exact reverse complement",
    ),
    // The four below close the corpus over its own output vocabulary: it *emitted*
    // constructs it never fed back in, so a whole class of input went unmeasured.
    // #1454 was one consequence — the corpus produced `inv` outputs but never an
    // input whose derivation had to *decide* on `inv` — and it was found by accident
    // rather than by the corpus, so there was no reason to think it was the only one.
    //
    // Measured on `5a19d1d2`, over 33,792 rows:
    //
    // | construct        | outputs | inputs |
    // |------------------|--------:|-------:|
    // | `inv`            |   3,113 |  3,072 |
    // | repeat `N[k]`    |   1,558 |      0 |
    // | `=`              |     152 |      0 |
    // | 3+ member allele |     170 |      0 |
    //
    // `inv` is the interesting row and the reason `inv_member` is narrower than the
    // other three. #1459 closed the raw count when it added `junction_interior_to_span`
    // — but every one of those 3,072 inputs pairs the `inv` with an insertion that
    // *overlaps* it, and none is lone. So the corpus could only ask what happens to an
    // `inv` already in conflict, never what happens to a well-formed one.
    (
        "inv_member",
        "an `inv` as an input, lone and beside a non-overlapping sibling — #1459 feeds `inv` only in a conflicting pairing",
    ),
    (
        "three_member_allele",
        "a three-member allele as an input — emitted on 170 rows, consumed on none",
    ),
    (
        "identity_member",
        "a no-change `=` member as an input — emitted on 152 rows, consumed on none",
    ),
    (
        "repeat_expansion",
        "a ranged repeat as an input — emitted on 1,558 rows, consumed on none",
    ),
];

/// The family drawn against the **long** cores, and the one shape in this file
/// whose point is its size rather than its geometry (#1460).
///
/// Kept out of `FAMILIES` on purpose. The nine families there are crossed with
/// every short core, and crossing a kilobase core with all of them would multiply
/// the dump cost while adding no coverage: `MAX_SPLIT_BLOCK` gates on the *length*
/// of the block being partitioned, not on how its members are arranged.
const LONG_FAMILY: (&str, &str) = (
    "long_block_inversion",
    "#1403 — a near-palindromic block straddling the MAX_SPLIT_BLOCK cap",
);

/// The long cores, as `(label, sequence)`. Two lengths that straddle the cap the
/// normalizer applies at 1024 bases (`merge::MAX_SPLIT_BLOCK`): one block that
/// fits under it and one that does not, which is the pair a change to the cap moves
/// between. The exact boundary is the one `partition_block`'s own comment records —
/// 1024 confluent, 1026 not.
///
/// The **label** is what lands in the dump's `reference` column, not the sequence:
/// a kilobase core would otherwise be repeated verbatim on every row of the dump.
/// Labels are as good a row identity as the sequence — they only have to be unique
/// and stable — and the short-core rows are untouched, so no existing key moves.
fn long_corpus_sequences() -> Vec<(String, String)> {
    [1024usize, 1100]
        .into_iter()
        .map(|len| (format!("nearpalindrome_{len}"), near_palindromic_core(len)))
        .collect()
}

/// A near-palindromic core: an exact reverse-complement palindrome with two
/// positions perturbed, so inverting the whole core changes exactly **four** bases
/// — the perturbed pair and their two partners.
///
/// The perturbation is what makes the shape usable. An exact palindrome inverts to
/// itself, so the inversion denotes nothing and normalizes to `=`; a random core
/// inverts to something differing at ~75% of its positions, so the equivalent
/// spelled-out allele would carry hundreds of members. Four differences is what
/// gives a block of this size a second spelling short enough to write down, which
/// is what makes it a confluence question at all.
fn near_palindromic_core(len: usize) -> String {
    // An odd length would silently produce a `len - 1` core while the caller's
    // label (`nearpalindrome_{len}`) still claimed `len` — and the label is the
    // dump's `reference` column, so the corpus would misreport the one property
    // this family exists to vary.
    assert!(
        len.is_multiple_of(2),
        "near-palindromic core length must be even, got {len}"
    );
    let half = len / 2;
    let mut state = (len as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let mut seq: Vec<u8> = (0..half)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize]
        })
        .collect();
    // Mirror it into an exact palindrome, then break it at two positions that are
    // not each other's partners — each perturbation costs two differences under
    // reverse complement, so two of them give the four the shape is named for.
    let mirror: Vec<u8> = (0..half).rev().map(|i| complement(seq[i])).collect();
    seq.extend(mirror);
    for offset in [10, half + 7] {
        seq[offset] = complement(seq[offset]);
    }
    String::from_utf8(seq).expect("ACGT is valid UTF-8")
}

/// The family drawn against the **reverse-complement designs**, and the one shape
/// in this file whose point is the *distance between* its changed runs (#1454's
/// shape, at a separation this corpus could not previously reach).
///
/// Kept out of [`FAMILIES`] for the same reason [`LONG_FAMILY`] is: it cannot be
/// drawn on an arbitrary core. Whether a reverse-complement block coincides with
/// its own reference at a given column is a property of the *sequence*, not of the
/// description, so a family crossed with the random 20-mers gets whatever
/// coincidence structure the xorshift happened to produce — which is why
/// `delins_hiding_an_inversion`, the only other family that builds a
/// reverse-complement block, emits its two pieces at **separation 0** and at no
/// other separation. A gate that engages only at separation >= 2 with multi-column
/// pieces therefore fires on no row of this corpus, and a change to one measures a
/// blast radius of zero for structural reasons rather than for safe ones.
///
/// The designs below fix that by building the *reference* around the intended
/// coincidence pattern instead of reading it off a random core — the same trade
/// `long_corpus_sequences` makes for scale.
const REVCOMP_FAMILY: (&str, &str) = (
    "separated_revcomp_runs",
    "a reverse-complement block whose changed runs are multi-column and separated \
     by 0, 1, 2, 4 or 8 unchanged columns",
);

/// Reference columns in each changed run of a `separated_revcomp_runs` design.
///
/// **Two, never one**, and that is the whole point of the family. A one-column run
/// is a substitution, and the separation rules for substitutions are already dense
/// in this corpus (`three_member_allele`, `dup_plus_sub`, `del_plus_sub`). What no
/// family builds is a *multi-column* changed run — a piece that has to be typed
/// `delins` or `inv` rather than `>` — at a separation greater than zero.
const REVCOMP_RUN_WIDTH: usize = 2;

/// Unchanged reference columns between consecutive changed runs.
///
/// `0` is the merge case (`general.md:34`, `DNA/delins.md:16-18` — adjacent changes
/// are one description), `1` is the codon carve-out (`general.md:35` — merged only
/// when both changes fall inside one codon, which needs a reading frame and so can
/// never apply on `g.`), and `2`/`4`/`8` are the split cases. `8` is also
/// `COALESCE_MAX_SEPARATION`, so the sweep sits on that boundary rather than
/// stopping short of it.
const REVCOMP_SEPARATIONS: &[usize] = &[0, 1, 2, 4, 8];

/// Offsets into the core at which a design's block is planted.
///
/// Three, because the *axis* geometry is what varies with the offset and not the
/// block: at 0 the block starts in the `cx` reference's 5'UTR, at 8 it runs past
/// `CDS_END` into the 3'UTR, and every offset crosses at least one of the two
/// `EXON_SPANS` junctions. The block itself is identical at all three, so a
/// difference between them is a difference the axis made.
const REVCOMP_SPAN_OFFSETS: &[usize] = &[0, 4, 8];

/// Length of a `separated_revcomp_runs` core.
///
/// The same 20 as [`corpus_sequences`] deliberately: `CDS_START`/`CDS_END` and
/// `EXON_SPANS` are both tuned to a 20-base core, so matching it is what lets these
/// designs run on all three axes rather than being restricted to `g.`/`c.` the way
/// the kilobase cores are.
const REVCOMP_CORE_LEN: usize = 20;

/// One `separated_revcomp_runs` design: a reference, a block inside it, and the
/// alternate that block's inversion produces.
///
/// Both sides are recorded because this family is the only one in the file that
/// knows its own **ground truth**. Every other family is a set of descriptions and
/// nothing more, so the harness can ask whether an output *moved* but never whether
/// it is *right*. Here the generator built `alternate` itself, which is what makes
/// the three oracles in the test module possible at all.
struct RevcompDesign {
    /// The dump's `reference` column. A label rather than the core, matching
    /// [`long_corpus_sequences`], so a row names the design it came from.
    label: String,
    /// The 20-base reference sequence.
    core: String,
    /// 0-based offset of the block within `core`.
    offset: usize,
    /// Unchanged columns between consecutive changed runs.
    separation: usize,
    /// The block's reference bases: `core[offset..offset + span.len()]`.
    span: String,
    /// `revcomp(span)` — what the block reads after the inversion. This is the
    /// **intended alternate**, and the apply oracle checks the normalized output
    /// against it rather than against the input.
    alternate: String,
}

/// The block of a design at `separation`, as literal reference bases.
///
/// Written out rather than generated, because each one is a small hand
/// construction that a reader should be able to check by eye — and
/// `the_designed_blocks_have_the_coincidence_structure_they_claim` re-derives the
/// run structure from the string, so a typo fails rather than silently producing a
/// different shape. `AATGCACA` (separation 4) is the worked case: it reverse
/// complements to `TGTGCATT`, coinciding at its four interior columns and leaving
/// two 2-base changes.
///
/// **Why the odd separation needs three runs.** A column of a block coincides with
/// its own reverse complement iff `span[i] == complement(span[L-1-i])`, and
/// complement is an involution — so column `i` coincides exactly when its mirror
/// `L-1-i` does. The changed/unchanged pattern of *any* whole-block inversion is
/// therefore a palindrome, and an odd-length block's centre column is its own
/// mirror and can never coincide. Two changed runs separated by an odd gap would
/// need that gap to straddle the centre of an even-length block (impossible: the
/// gap is centred, so its length has the parity of the block) or to sit at the
/// centre of an odd-length block (impossible: the centre column is always changed).
/// So separation 1 is unreachable with two runs and reachable with three, whose
/// two gaps are mirror images of each other. That is a fact about reverse
/// complements, not a limitation of these five strings.
fn revcomp_span(separation: usize) -> &'static str {
    match separation {
        0 => "ACCA",
        1 => "AACTTGAA",
        2 => "AACGAA",
        4 => "AATGCACA",
        8 => "AACGTATACGAA",
        other => unreachable!("no design for separation {other}"),
    }
}

/// The changed runs of a design at `separation`, as `(start, length)` block-local.
///
/// Derived from the parameters rather than from the sequence, because at separation
/// **0** the two runs are adjacent and the sequence alone cannot tell them apart
/// from one run of four — and telling them apart is the entire content of the
/// merge case. `the_designed_blocks_have_the_coincidence_structure_they_claim`
/// pins the two views against each other: the union of these runs must be exactly
/// the set of columns the reverse complement actually changes.
fn revcomp_runs(separation: usize) -> Vec<(usize, usize)> {
    let count = if separation.is_multiple_of(2) { 2 } else { 3 };
    let stride = REVCOMP_RUN_WIDTH + separation;
    (0..count)
        .map(|n| (n * stride, REVCOMP_RUN_WIDTH))
        .collect()
}

/// The full design set: every separation at every offset.
fn revcomp_designs() -> Vec<RevcompDesign> {
    let mut designs = Vec::new();
    for &separation in REVCOMP_SEPARATIONS {
        let span = revcomp_span(separation);
        for &offset in REVCOMP_SPAN_OFFSETS {
            // A design whose block ran off the end of the core would silently
            // produce a truncated span while its label still named the separation
            // it no longer has.
            assert!(
                offset + span.len() <= REVCOMP_CORE_LEN,
                "separation-{separation} block does not fit at offset {offset}"
            );
            designs.push(RevcompDesign {
                label: format!("revcomp_sep{separation}_at{offset}"),
                core: revcomp_core(offset, span),
                offset,
                separation,
                span: span.to_string(),
                alternate: revcomp(span),
            });
        }
    }
    designs
}

/// A design's core: filler everywhere, with `span` planted at `offset`.
///
/// The filler is period-4 `GCTA`, one rotation out of phase with the `ACGT` pad
/// [`padded`] wraps a `g.` core in. That is deliberate — an in-phase filler would
/// let the pad's own rotation continue straight through the flank, so a tract
/// ending at the block edge could extend into the pad and a shift's stopping point
/// would become a property of the padding rather than of the block.
fn revcomp_core(offset: usize, span: &str) -> String {
    let mut core: Vec<u8> = (0..REVCOMP_CORE_LEN).map(|i| b"GCTA"[i % 4]).collect();
    core[offset..offset + span.len()].copy_from_slice(span.as_bytes());
    String::from_utf8(core).expect("GCTA and the designed spans are valid UTF-8")
}

/// The three spellings of one design's variant, on `axis`.
///
/// All three denote the same bases by construction, so they are a **confluence
/// class**: `the_three_spellings_of_a_design_converge` requires one output from
/// all three, and the apply oracle requires that output to produce `alternate`.
///
/// The split spelling is the load-bearing one. It writes each changed run as its
/// own `delins` member, which is the form the spec mandates at separation >= 2
/// (`general.md:34`) and forbids at separation 0 (`DNA/delins.md:16-18`) — so
/// across the five separations the same spelling is required, carved out, and
/// forbidden in turn, and the normalizer has to tell those cases apart.
fn revcomp_inputs_for(axis: &str, design: &RevcompDesign) -> Vec<String> {
    let prefix = prefix_for(axis);
    let p = |i: usize| hgvs_pos(axis, design.offset + i);
    let last = design.span.len() - 1;
    let members: Vec<String> = revcomp_runs(design.separation)
        .iter()
        .map(|&(start, len)| {
            format!(
                "{}_{}delins{}",
                p(start),
                p(start + len - 1),
                &design.alternate[start..start + len]
            )
        })
        .collect();
    vec![
        format!("{prefix}{}_{}inv", p(0), p(last)),
        format!("{prefix}{}_{}delins{}", p(0), p(last), design.alternate),
        format!("{prefix}[{}]", members.join(";")),
    ]
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        // Not a catch-all: a silent fallback would map an unexpected byte to `C`
        // and hand `the_long_cores_are_near_palindromes` a core that looks
        // near-palindromic because the complement lied, not because the sequence
        // is. The cores are built from `ACGT` here, so this is unreachable —
        // which is the reason to say so rather than to absorb it.
        other => panic!("complement of non-ACGT byte {other:#04x}"),
    }
}

/// The two spellings of one whole-core inversion: the `inv` itself, and the same
/// change written as the individual substitutions it makes. On a near-palindromic
/// core those are the same variant, so a normalizer that is confluent has to bring
/// them to one form — and whether it can depends on whether the block was capped.
fn long_inputs_for(axis: &str, core: &str) -> Vec<String> {
    let bytes = core.as_bytes();
    let prefix = prefix_for(axis);
    let inverted: Vec<u8> = bytes.iter().rev().map(|b| complement(*b)).collect();
    let members: Vec<String> = (0..bytes.len())
        .filter(|&i| bytes[i] != inverted[i])
        .map(|i| {
            format!(
                "{}{}>{}",
                hgvs_pos(axis, i),
                bytes[i] as char,
                inverted[i] as char
            )
        })
        .collect();
    let whole = format!(
        "{prefix}{}_{}inv",
        hgvs_pos(axis, 0),
        hgvs_pos(axis, bytes.len() - 1)
    );
    // A core whose inversion changes nothing would emit an empty allele, which is
    // not a description. `the_long_cores_are_near_palindromes` pins that the
    // construction above never gets there.
    if members.is_empty() {
        return vec![whole];
    }
    vec![whole, format!("{prefix}[{}]", members.join(";"))]
}

/// Offset of the first tandem three-base repeat in `core`, if it has one.
///
/// Three bases rather than one, because `repeated.md:21` permits a repeat description
/// on a `c.` reference **only** for units whose length is a multiple of three — a
/// mononucleotide repeat there is spec-invalid (`c.2686A[10]` is called out by name),
/// and a corpus that emits invalid input measures nothing useful. Three also lets the
/// ranged spelling be used throughout, which is the form the pending Community
/// Consultation proposal in that file's NOTE would leave as the only legal one.
///
/// Only 21 of the 48 corpus sequences contain such a tract, so `repeat_expansion`
/// emits nothing for the rest. That is deliberate: synthesising a tract would mean
/// describing a repeat the reference does not have.
fn first_tandem_triplet(core: &str) -> Option<usize> {
    let b = core.as_bytes();
    (0..b.len().saturating_sub(5)).find(|&i| b[i..i + 3] == b[i + 3..i + 6])
}

/// A three-base payload for the leading piece of `delins_hiding_an_inversion`,
/// chosen so the piece is a genuine, non-inverting change of `span`.
///
/// Two degeneracies have to be avoided, and both silently weaken the family rather
/// than failing it. A payload equal to `span` makes the piece a no-op, leaving a
/// one-piece shape that no longer tests a *split*. A payload equal to
/// `revcomp(span)` makes the leading piece an inversion too, so the case no longer
/// isolates a single mis-typed piece. Three candidates suffice: one span rules out
/// at most two of them, and `the_head_payload_is_never_a_no_op_or_an_inversion`
/// checks that exhaustively over all 64 spans rather than trusting the argument.
fn head_payload(span: &str) -> &'static str {
    let inverted = revcomp(span);
    ["GGG", "CCC", "ACA"]
        .into_iter()
        .find(|c| *c != span && *c != inverted)
        .expect("three candidates cannot all collide with one three-base span")
}

/// Reverse complement of `s`, which must be over `ACGT`.
///
/// Used to build the `delins_hiding_an_inversion` family, where the point is that a
/// piece's payload equals the reverse complement of its own span and so should be
/// typed `inv` rather than `delins`.
fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => unreachable!("corpus sequences are over ACGT, found {other}"),
        })
        .collect()
}

fn dump(seeds: u32) -> Vec<Row> {
    let mut rows = Vec::new();
    for (label, core) in long_corpus_sequences() {
        for (axis, direction, dir_label) in axes_and_directions() {
            // The multi-exon axis is skipped for the long cores, and that is a
            // correctness fix rather than a coverage choice.
            //
            // `EXON_SPANS` is a fixed 7/7/6 split of a **20-base** core, chosen
            // so the junctions land inside the shapes the short corpus builds.
            // `multi_exon_provider` hands `Transcript::new` the whole core
            // regardless, so on a 1024-base core the transcript declares 1024
            // bases while its exon table maps 20 and `spliced_genomic` emits
            // only those 20 — every position past 20 is declared and unmapped.
            // Rows drawn against that provider are not wrong-looking (all eight
            // came out as fixed points) but they are not measuring anything
            // either, and a future change touching genomic projection would get
            // nonsense there that reads as real coverage.
            //
            // Scaling `EXON_SPANS` to the core instead would move the junction
            // geometry the short corpus is tuned around. The long cores exist
            // for `long_block_inversion`, which is about block size and not
            // exon structure, so they lose nothing here.
            if axis == "cx" {
                continue;
            }
            let normalizer = normalizer_for(axis, &core, direction);
            for input in long_inputs_for(axis, &core) {
                let output = normalize_through(&normalizer, &input);
                rows.push(Row {
                    reference: label.clone(),
                    core: core.clone(),
                    axis,
                    direction: dir_label,
                    family: LONG_FAMILY.0,
                    was_fixed_point: output == input,
                    input,
                    output,
                });
            }
        }
    }
    // The reverse-complement designs are drawn against their own references for the
    // same reason the long cores are: the property they vary — where a
    // reverse-complement block coincides with its own reference — is a property of
    // the sequence, so it cannot be crossed with the random cores below.
    for design in revcomp_designs() {
        for (axis, direction, dir_label) in axes_and_directions() {
            let normalizer = normalizer_for(axis, &design.core, direction);
            for input in revcomp_inputs_for(axis, &design) {
                let output = normalize_through(&normalizer, &input);
                rows.push(Row {
                    reference: design.label.clone(),
                    core: design.core.clone(),
                    axis,
                    direction: dir_label,
                    family: REVCOMP_FAMILY.0,
                    was_fixed_point: output == input,
                    input,
                    output,
                });
            }
        }
    }
    for core in corpus_sequences(seeds) {
        for (axis, direction, dir_label) in axes_and_directions() {
            // One normalizer per (axis, core, direction) cell, reused across every
            // family and every input in it — see `normalizer_for`.
            let normalizer = normalizer_for(axis, &core, direction);
            for (family, _) in FAMILIES {
                for input in inputs_for(family, axis, &core) {
                    let output = normalize_through(&normalizer, &input);
                    rows.push(Row {
                        reference: core.clone(),
                        core: core.clone(),
                        axis,
                        direction: dir_label,
                        family,
                        was_fixed_point: output == input,
                        input,
                        output,
                    });
                }
            }
        }
    }
    rows
}

/// The axis/direction matrix every core is drawn against, with each direction's
/// dump label. Shared by both loops in [`dump`] so the two cannot disagree about a
/// label — a dump that spelled one direction differently in one loop would compare
/// as an entirely separate set of rows.
fn axes_and_directions() -> Vec<(&'static str, ShuffleDirection, &'static str)> {
    [
        ("g", ShuffleDirection::ThreePrime),
        ("g", ShuffleDirection::FivePrime),
        ("c", ShuffleDirection::ThreePrime),
        ("c", ShuffleDirection::FivePrime),
        ("cx", ShuffleDirection::ThreePrime),
        ("cx", ShuffleDirection::FivePrime),
    ]
    .into_iter()
    .map(|(axis, direction)| {
        let label = match direction {
            ShuffleDirection::ThreePrime => "3prime",
            ShuffleDirection::FivePrime => "5prime",
            // `#[non_exhaustive]`: a new direction must be added to the matrix
            // above deliberately, not silently mislabelled in a dump that a later
            // revision will be diffed against.
            other => unreachable!("unhandled shuffle direction {other:?}"),
        };
        (axis, direction, label)
    })
    .collect()
}

/// Deterministic 20-mers, two per seed. Same xorshift as `sweep_sequences` so the
/// two corpora explore comparable sequence space; kept here rather than shared
/// because an example cannot reach `tests/`.
fn corpus_sequences(seeds: u32) -> Vec<String> {
    let mut sequences = Vec::with_capacity(2 * seeds as usize);
    for seed in 0..seeds {
        for alphabet in [b"AT".as_slice(), b"ACGT".as_slice()] {
            let mut state = u64::from(seed).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            sequences.push(
                (0..20)
                    .map(|_| {
                        state ^= state << 13;
                        state ^= state >> 7;
                        state ^= state << 17;
                        alphabet[(state % alphabet.len() as u64) as usize] as char
                    })
                    .collect(),
            );
        }
    }
    sequences
}

/// 1-based position of core offset `i` on the axis in question. A `g.` row sits in
/// a padded contig; a `c.` row addresses the transcript directly, so the core *is*
/// the transcript and offset 0 is position 1.
fn pos(axis: &str, i: usize) -> usize {
    match axis {
        "g" => PAD_OFFSET + 1 + i,
        _ => 1 + i,
    }
}

/// Render core offset `i` as an HGVS position on `axis`.
///
/// A `c.` position past `CDS_END` is in the 3'UTR and is spelled `*N`, not `N`.
/// Getting this wrong is not cosmetic: the corpus would either emit descriptions
/// that mean the wrong base or, as it did before this was fixed, stop short of the
/// junction entirely and report a confident zero for changes that only affect
/// junction-crossing shapes — which is where #1398, #1418 and both halves of #1426
/// lived.
fn hgvs_pos(axis: &str, i: usize) -> String {
    if axis == "g" {
        return pos(axis, i).to_string();
    }
    let p = pos(axis, i) as u64;
    let (cds_start, cds_end) = cds_bounds(axis);
    if p < cds_start {
        // 5'UTR. Unreachable on the single-exon reference, where `CDS_START` is 1
        // and no transcript position precedes the CDS — which is why `c.-N` never
        // appeared in this corpus before #1478.
        format!("-{}", cds_start - p)
    } else if p <= cds_end {
        (p - cds_start + 1).to_string()
    } else {
        format!("*{}", p - cds_end)
    }
}

/// CDS bounds of the coding reference `axis` is drawn against.
///
/// The general form is a no-op for `c.`: with `CDS_START == 1` the middle arm of
/// `hgvs_pos` reduces to `p`, and the 5'UTR arm is unreachable, so every existing
/// `c.` row spells exactly as it did before this function existed.
fn cds_bounds(axis: &str) -> (u64, u64) {
    if axis == "cx" {
        (CDS_START_MULTI, CDS_END_MULTI)
    } else {
        (CDS_START, CDS_END)
    }
}

/// The `<accession>:<axis letter>.` prefix a description on `axis` carries.
///
/// The corpus's axis token and the HGVS axis letter are different things. `cx` is a
/// corpus dimension — the multi-exon coding reference — whose descriptions are
/// spelled `c.` like any other coding variant. Interpolating the token directly, as
/// both call sites used to, would emit `NM_TESTX.1:cx.`, which is not HGVS and
/// would be recorded as a parse error rather than a measurement.
fn prefix_for(axis: &str) -> String {
    match axis {
        "g" => format!("{GENOMIC_CONTIG}:g."),
        "cx" => format!("{TX_MULTI_ACCESSION}:c."),
        _ => format!("{TX_ACCESSION}:c."),
    }
}

fn inputs_for(family: &str, axis: &str, core: &str) -> Vec<String> {
    let bytes = core.as_bytes();
    let prefix = prefix_for(axis);
    let mut out = Vec::new();
    // The widest shape reads five consecutive offsets (`s`..`s+4`), so `s` runs to
    // `len - 5`. It deliberately runs that far rather than stopping short: on the
    // coding axis the last shapes straddle `CDS_END` and are spelled with `c.*N`,
    // which is the junction region four of this campaign's five defects lived in.
    let last = bytes.len().saturating_sub(4);
    for s in 0..last {
        let p = |i: usize| hgvs_pos(axis, i);
        let base = |i: usize| bytes[i] as char;
        let other = |i: usize| if base(i) == 'A' { 'C' } else { 'A' };
        match family {
            // Two delins separated by exactly one unchanged base, plus the
            // spanning spelling of the same variant. The pair is the confluence
            // question; the spanning form should be a fixed point either way.
            "split_vs_spanning_delins" => {
                let gap = base(s + 2);
                out.push(format!(
                    "{prefix}[{}_{}delinsGGG;{}_{}delinsGG]",
                    p(s),
                    p(s + 1),
                    p(s + 3),
                    p(s + 4)
                ));
                out.push(format!("{prefix}{}_{}delinsGGG{gap}GG", p(s), p(s + 4)));
            }
            "dup_plus_sub" => out.push(format!(
                "{prefix}[{}_{}dup;{}{}>{}]",
                p(s),
                p(s + 1),
                p(s + 3),
                base(s + 3),
                other(s + 3)
            )),
            "adjacent_junction_ins" => out.push(format!(
                "{prefix}[{}_{}insAA;{}_{}insTT]",
                p(s),
                p(s + 1),
                p(s + 1),
                p(s + 2)
            )),
            "dup_plus_ins" => out.push(format!(
                "{prefix}[{}_{}dup;{}_{}insCC]",
                p(s),
                p(s + 2),
                p(s + 3),
                p(s + 4)
            )),
            "del_plus_sub" => out.push(format!(
                "{prefix}[{}_{}del;{}{}>{}]",
                p(s),
                p(s + 1),
                p(s + 3),
                base(s + 3),
                other(s + 3)
            )),
            // A one-base dup and a substitution on that same base: the shape
            // `detect_overlap_conflicts` keys on, and the one #1307 produced out
            // of bounds (`g.[24dup;24C>G]` -> `g.[24C>G;24_25insC]`).
            "coincident_bounds" => out.push(format!(
                "{prefix}[{}dup;{}{}>{}]",
                p(s),
                p(s),
                base(s),
                other(s)
            )),
            // An insertion anchored at the junction *inside* another member's
            // span — mutalyzer `EOVERLAP`, #486. Both spellings the second
            // detector has to see: the span deleted, and the span inverted.
            "junction_interior_to_span" => {
                out.push(format!(
                    "{prefix}[{}_{}del;{}_{}insCC]",
                    p(s),
                    p(s + 1),
                    p(s),
                    p(s + 1)
                ));
                out.push(format!(
                    "{prefix}[{}_{}inv;{}_{}insCC]",
                    p(s),
                    p(s + 1),
                    p(s),
                    p(s + 1)
                ));
            }
            // Two spans that intersect without coinciding — accepted by strict
            // mode until #1451, which reported 10,499 of 25,848 such pairs as
            // denoting nothing while this corpus measured zero movement.
            "partial_overlap_spans" => out.push(format!(
                "{prefix}[{}_{}del;{}_{}del]",
                p(s),
                p(s + 2),
                p(s + 1),
                p(s + 3)
            )),
            // Nested spans sharing their 5' bound. Kept apart from the partial
            // overlap above because containment and mere intersection are
            // different questions to a detector that keys on bounds.
            "nested_spans" => out.push(format!(
                "{prefix}[{}_{}del;{}_{}del]",
                p(s),
                p(s + 1),
                p(s),
                p(s + 3)
            )),
            // A length-preserving spanning delins over `s..s+5`, whose payload is a
            // non-inverting three-base change of `s..s+2` followed by the exact
            // reverse complement of `s+3..s+5`. The natural split is therefore
            // `[s_s+2delins<head>; s+3_s+5inv]` — the trailing piece *is* an
            // inversion, and #1454 is that ferro spells it `delins` on the first
            // pass and only re-types it to `inv` on a second, so the first output is
            // not a fixed point.
            //
            // Both spellings are emitted: the spanning form, and the pre-split form
            // with the trailing piece written as `delins`. The latter is #1454's
            // pass-1 output fed back in, which is the shape every other family here
            // fails to build — measured on `ee0e37ac`, all 8,706 `delins` members
            // the other five families emit have a payload that is *not* the reverse
            // complement of its span, which is why the corpus reported zero
            // non-idempotent rows while #1454 was open.
            //
            // **Three bases per piece, and the width is load-bearing.** A two-base
            // head with a three-base tail reaches zero non-idempotent rows: the
            // trailing piece gets typed `inv` on the first pass and the defect never
            // appears. Three-and-three reaches 30. An odd-length span also can never
            // equal its own reverse complement, so neither piece degenerates to a
            // no-op the way a two-base piece does on `AT`, `TA`, `CG` or `GC`.
            //
            // This shape reads six offsets where every other family reads five, so
            // it is bounded here rather than by widening the shared `last` — moving
            // that would change what all five existing families emit and make every
            // dump taken before this commit incomparable.
            "delins_hiding_an_inversion" => {
                if s + 6 > bytes.len() {
                    continue;
                }
                let head = head_payload(&core[s..s + 3]);
                let tail = revcomp(&core[s + 3..s + 6]);
                out.push(format!("{prefix}{}_{}delins{head}{tail}", p(s), p(s + 5)));
                out.push(format!(
                    "{prefix}[{}_{}delins{head};{}_{}delins{tail}]",
                    p(s),
                    p(s + 2),
                    p(s + 3),
                    p(s + 5)
                ));
            }
            // Three bases, never two: an even-length inversion is a no-op on a
            // palindromic span (`AT`, `TA`, `CG`, `GC`), which would make the input
            // an identity rather than an inversion and quietly weaken the family.
            // An odd-length span can never equal its own reverse complement.
            "inv_member" => {
                out.push(format!("{prefix}{}_{}inv", p(s), p(s + 2)));
                out.push(format!(
                    "{prefix}[{}_{}inv;{}{}>{}]",
                    p(s),
                    p(s + 2),
                    p(s + 4),
                    base(s + 4),
                    other(s + 4)
                ));
            }
            // Three substitutions, each separated from the next by one unchanged
            // base, so the allele is well-formed rather than a set of adjacent
            // changes that should have been written as one edit.
            "three_member_allele" => out.push(format!(
                "{prefix}[{}{}>{};{}{}>{};{}{}>{}]",
                p(s),
                base(s),
                other(s),
                p(s + 2),
                base(s + 2),
                other(s + 2),
                p(s + 4),
                base(s + 4),
                other(s + 4)
            )),
            // `=` is a member the normalizer *writes* when a piece cancels, so it is
            // a shape a consumer can hold and feed back. Paired with a real change so
            // the allele is not wholly vacuous.
            "identity_member" => out.push(format!(
                "{prefix}[{}_{}=;{}{}>{}]",
                p(s),
                p(s + 1),
                p(s + 3),
                base(s + 3),
                other(s + 3)
            )),
            // A genuine two-copy reference tract expanded to three, in the ranged
            // spelling. Emits nothing for a core with no tandem triplet — see
            // `first_tandem_triplet`.
            "repeat_expansion" => {
                if s != 0 {
                    continue;
                }
                if let Some(t) = first_tandem_triplet(core) {
                    let unit = &core[t..t + 3];
                    out.push(format!("{prefix}{}_{}{unit}[3]", p(t), p(t + 5)));
                    out.push(format!("{prefix}{}_{}{unit}[1]", p(t), p(t + 5)));
                }
            }
            _ => unreachable!("unknown family {family}"),
        }
    }
    out
}

/// Normalize, reporting a decline or a panic as data rather than aborting the dump.
/// A row that errors on one revision and succeeds on the other is exactly the
/// cheap-vs-expensive distinction the diff needs to see, so it must not be dropped.
/// The provider a row is drawn against, selected by axis label.
///
/// Factored out of [`normalize_one`] so the `--verify-spdi` pass reads the same
/// reference the row was normalized against. Building a second, subtly different
/// provider there would make the verification answer a question about the wrong
/// sequence — and it would do so silently, since both would still produce keys.
fn provider_for(axis: &str, core: &str) -> MockProvider {
    match axis {
        "g" => genomic_provider(core),
        "cx" => multi_exon_provider(core),
        _ => coding_provider(core),
    }
}

/// The normalizer one `(axis, core, direction)` cell of the matrix is dumped
/// through.
///
/// Split out of [`normalize_one`] so [`dump`] can build it **once per cell**
/// instead of once per row. `provider_for` is not cheap: every call re-pads the
/// core with 512 bases (`padded`), and the coding and multi-exon axes also build a
/// `Transcript` — on the 1024-base long cores that is a ~1.5 kB copy per input, and
/// the corpus has tens of thousands of rows. Nothing about it varies with the row.
fn normalizer_for(axis: &str, core: &str, direction: ShuffleDirection) -> Normalizer<MockProvider> {
    Normalizer::with_config(
        provider_for(axis, core),
        NormalizeConfig::default().with_direction(direction),
    )
}

/// One row's output, through a normalizer the caller already holds.
///
/// The sentinel vocabulary (`<parse-error>` / `<declined>` / `<panic>`) is what
/// `the_corpus_emits_a_block_past_the_split_cap` reads, so it lives here rather
/// than at each call site.
fn normalize_through(normalizer: &Normalizer<MockProvider>, input: &str) -> String {
    let Ok(variant) = parse_hgvs(input) else {
        return "<parse-error>".to_string();
    };
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        normalizer.normalize(&variant)
    })) {
        Ok(Ok(v)) => v.to_string(),
        Ok(Err(_)) => "<declined>".to_string(),
        Err(_) => "<panic>".to_string(),
    }
}

/// [`normalize_through`], building the normalizer for one row.
///
/// `#[cfg(test)]` because that is now the whole of its use: the callers that
/// genuinely normalize a *single* description are the unit tests below and the
/// idempotency re-check, where there is no matrix cell to amortise a normalizer
/// over. [`dump`] holds one per cell via [`normalizer_for`].
#[cfg(test)]
fn normalize_one(axis: &str, core: &str, input: &str, direction: ShuffleDirection) -> String {
    normalize_through(&normalizer_for(axis, core, direction), input)
}

fn padded(core: &str) -> String {
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    format!("{pad}{core}{pad}")
}

fn genomic_provider(core: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence(GENOMIC_CONTIG, padded(core));
    provider
}

/// Genomic sequence for the multi-exon transcript: the core's exon blocks laid out
/// with `INTRON_LEN` bases of filler between them, padded on both sides.
///
/// The filler continues the pad's period-4 `ACGT` rotation rather than using a
/// distinctive motif, so an intron cannot accidentally extend a repeat tract that
/// ends at an exon boundary — which would make a shift's stopping point a property
/// of the filler rather than of the exon edge.
fn spliced_genomic(core: &str) -> String {
    // `EXON_SPANS` maps a fixed-length core, and *both* directions of a mismatch
    // are defects rather than data (#1624). A shorter core slices off the end —
    // which is how a 16-byte label reached here and took the whole `--verify-spdi`
    // pass with it, reported as `end byte index 20 is out of bounds`, a message
    // that names neither the exon table nor the caller. A longer one is the
    // silent direction: `multi_exon_provider` hands `Transcript::new` the whole
    // core, so the transcript would declare bases the exon table does not map and
    // the rows drawn against it would read as coverage while measuring nothing.
    // That is why `dump` skips the `cx` axis for the long cores; this says so
    // where the assumption actually lives.
    let mapped = EXON_SPANS.last().expect("the exon table is not empty").1;
    assert_eq!(
        core.len(),
        mapped,
        "the multi-exon axis maps a {mapped}-base core; got {} bases ({core})",
        core.len()
    );
    let pad = "ACGT".repeat(PAD_OFFSET / 4);
    let mut out = String::with_capacity(2 * PAD_OFFSET + core.len() + 2 * INTRON_LEN);
    out.push_str(&pad);
    for (n, (from, to)) in EXON_SPANS.iter().enumerate() {
        if n > 0 {
            out.push_str(&"ACGT".repeat(INTRON_LEN.div_ceil(4))[..INTRON_LEN]);
        }
        out.push_str(&core[*from..*to]);
    }
    out.push_str(&pad);
    out
}

/// Genomic start (1-based) of exon `n` under the `spliced_genomic` layout.
fn exon_genomic_start(n: usize) -> u64 {
    let exonic_before: usize = EXON_SPANS[..n].iter().map(|(f, t)| t - f).sum();
    (PAD_OFFSET + 1 + exonic_before + n * INTRON_LEN) as u64
}

/// The multi-exon coding reference (#1478): three exons, `CDS_START_MULTI` past the
/// transcript start so a 5'UTR exists, on the plus strand.
fn multi_exon_provider(core: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    let exons: Vec<Exon> = EXON_SPANS
        .iter()
        .enumerate()
        .map(|(n, (from, to))| {
            let g_start = exon_genomic_start(n);
            Exon::with_genomic(
                n as u32 + 1,
                *from as u64 + 1,
                *to as u64,
                g_start,
                g_start + (to - from) as u64 - 1,
            )
        })
        .collect();
    let g_start = exon_genomic_start(0);
    let g_end = exon_genomic_start(EXON_SPANS.len() - 1)
        + (EXON_SPANS[EXON_SPANS.len() - 1].1 - EXON_SPANS[EXON_SPANS.len() - 1].0) as u64
        - 1;
    let transcript = Transcript::new(
        TX_MULTI_ACCESSION.to_string(),
        Some("SYNTHX".to_string()),
        Strand::Plus,
        core.to_string(),
        Some(CDS_START_MULTI),
        Some(CDS_END_MULTI),
        exons,
        Some(TX_MULTI_CONTIG.to_string()),
        Some(g_start),
        Some(g_end),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    provider.add_genomic_sequence(TX_MULTI_CONTIG, spliced_genomic(core));
    provider.add_transcript(transcript);
    provider
}

fn coding_provider(core: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    let tx_len = core.len() as u64;
    let g_start = PAD_OFFSET as u64 + 1;
    let g_end = PAD_OFFSET as u64 + tx_len;
    let exon = Exon::with_genomic(1, 1, tx_len, g_start, g_end);
    let transcript = Transcript::new(
        TX_ACCESSION.to_string(),
        Some("SYNTH".to_string()),
        Strand::Plus,
        core.to_string(),
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
    provider.add_genomic_sequence(TX_CONTIG, padded(core));
    provider.add_transcript(transcript);
    provider
}

/// A row's identity. Deliberately a tuple: see trap 1 in the module docs.
type Key = (String, String, String, String);

fn read_dump(path: &PathBuf) -> Result<BTreeMap<Key, (String, bool, String)>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("reading {}: {e}", path.display()))?;
    let mut rows = BTreeMap::new();
    for (n, line) in text.lines().enumerate() {
        if n == 0 {
            // Exact match, not a prefix. A dump whose columns were reordered or
            // renamed would otherwise be read positionally under the old meanings —
            // silently swapping, say, `output` and `was_fixed_point` and inverting
            // every migration verdict in the report.
            if line != HEADER.trim_end() {
                return Err(format!(
                    "{}:1: unexpected header\n  found:    {line}\n  expected: {}",
                    path.display(),
                    HEADER.trim_end()
                ));
            }
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != COLUMNS {
            return Err(format!(
                "{}:{}: expected {COLUMNS} columns, found {}",
                path.display(),
                n + 1,
                f.len()
            ));
        }
        // Strict, because this column decides the expensive/cheap split — the single
        // most consequential number in the report. Treating anything that is not
        // "true" as false would turn a corrupted dump into an understated migration,
        // which is the direction that gets a change waved through.
        let was_fixed_point = match f[6] {
            "true" => true,
            "false" => false,
            other => {
                return Err(format!(
                    "{}:{}: was_fixed_point must be `true` or `false`, found `{other}`",
                    path.display(),
                    n + 1
                ))
            }
        };
        let key = (
            f[0].to_string(),
            f[1].to_string(),
            f[2].to_string(),
            f[4].to_string(),
        );
        if rows
            .insert(
                key.clone(),
                (f[5].to_string(), was_fixed_point, f[3].to_string()),
            )
            .is_some()
        {
            return Err(format!(
                "{}:{}: duplicate key {key:?} — the dump is not uniquely keyed",
                path.display(),
                n + 1
            ));
        }
    }
    Ok(rows)
}

fn compare(before: &PathBuf, after: &PathBuf) -> Result<String, String> {
    let a = read_dump(before)?;
    let b = read_dump(after)?;

    let only_before = a.keys().filter(|k| !b.contains_key(*k)).count();
    let only_after = b.keys().filter(|k| !a.contains_key(*k)).count();
    if only_before + only_after > 0 {
        return Err(format!(
            "the two dumps do not cover the same corpus ({only_before} rows only in the \
             baseline, {only_after} only in the candidate). Re-dump both with the same \
             --seeds; a changed generator or seed count makes them incomparable."
        ));
    }

    let mut moved: Vec<(&Key, &str, &str, &str)> = Vec::new();
    for (key, (old, _, family)) in &a {
        let (new, _, new_family) = &b[key];
        // The family is part of what the row *means*, not just a label: the
        // per-family table below is read as "this change touched only X". If the two
        // dumps disagree about which family a key belongs to, the family definitions
        // moved between revisions and that table would attribute movement to the
        // wrong shape. Refuse rather than silently report the baseline's view.
        if family != new_family {
            return Err(format!(
                "the two dumps disagree on the shape family for {key:?} \
                 (baseline `{family}`, candidate `{new_family}`). The family \
                 definitions changed between revisions, so a per-family breakdown \
                 would be misattributed; re-dump both at the same definitions."
            ));
        }
        if old != new {
            moved.push((key, old, new, family));
        }
    }

    let mut by_family: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    for (_, _, family) in a.values() {
        by_family.entry(family.as_str()).or_default().0 += 1;
    }
    for (_, _, _, family) in &moved {
        by_family.entry(family).or_default().1 += 1;
    }

    // Expensive == the input was its own normalized form on the BASELINE, so a
    // consumer has that string stored. Cheap == it was not, so nothing stored it.
    let expensive = moved.iter().filter(|(key, _, _, _)| a[*key].1).count();

    let mut r = String::new();
    let _ = writeln!(r, "## Representation stability\n");
    let _ = writeln!(r, "| | rows |");
    let _ = writeln!(r, "|---|---|");
    let _ = writeln!(r, "| total | {} |", a.len());
    let pct = if a.is_empty() {
        0.0
    } else {
        100.0 * moved.len() as f64 / a.len() as f64
    };
    let _ = writeln!(r, "| **moved** | **{} ({pct:.1}%)** |", moved.len());
    let _ = writeln!(
        r,
        "| of which previously **accepted** (a migration) | {expensive} |"
    );
    let _ = writeln!(
        r,
        "| of which previously **not** a fixed point (free) | {} |",
        moved.len() - expensive
    );
    let _ = writeln!(r, "\n### By shape family\n");
    let _ = writeln!(r, "| family | rows | moved |");
    let _ = writeln!(r, "|---|---|---|");
    for (family, (total, m)) in &by_family {
        let _ = writeln!(r, "| `{family}` | {total} | {m} |");
    }
    if moved.is_empty() {
        // A zero is the answer most likely to be quoted, and the one most likely
        // to be wrong: it is a statement about the families above and nothing
        // else. Three PRs in a row (#1403, #1445, #1451) reported `0 of 18,432`
        // from a corpus that built no allele of the shape they changed, and each
        // had to hand-write this caveat into its own PR body (#1456). Say it here
        // instead, and name what was actually covered.
        let covered: Vec<&str> = by_family.keys().copied().collect();
        let _ = writeln!(
            r,
            "\n**No row moved.** That is a measured zero for the {} shape famil{} above and says \
             nothing about a shape this corpus does not build: {}. If the change under test is \
             scoped to a shape absent from that list, this harness has not measured it — add the \
             family rather than quoting the zero.",
            covered.len(),
            if covered.len() == 1 { "y" } else { "ies" },
            covered
                .iter()
                .map(|f| format!("`{f}`"))
                .collect::<Vec<_>>()
                .join(", ")
        );
    } else {
        let _ = writeln!(r, "\n### Sample moved rows\n\n```text");
        for (key, old, new, _) in moved.iter().take(10) {
            let _ = writeln!(
                r,
                "[{} {}] {}\n   before: {old}\n   after : {new}",
                key.1, key.2, key.3
            );
        }
        let _ = writeln!(r, "```");
    }
    let _ = writeln!(
        r,
        "\nThis corpus is enriched for the shapes that churn, so the percentage above is a \
         rate **within these shape families**, not a library-wide rate."
    );
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A row of one of the thirteen sequence-keyed families, where the
    /// `reference` column and the sequence are the same string.
    fn row(axis: &'static str, reference: &str, input: &str, output: &str) -> Row {
        Row {
            reference: reference.to_string(),
            core: reference.to_string(),
            axis,
            direction: "3prime",
            family: "test",
            input: input.to_string(),
            output: output.to_string(),
            was_fixed_point: false,
        }
    }

    /// All three verdicts, on real rows drawn from the corpus's own cores.
    ///
    /// The one that earns the test is `Unverifiable`. `canonical_spdi` declines an
    /// allele whose members overlap, and this corpus is deliberately full of those
    /// — 20,526 of 78,028 rows. Classifying a decline as `Different` would report
    /// twenty thousand findings and bury the eighty real ones, which is exactly
    /// what the first hand-written cut of this check did.
    #[test]
    fn verify_row_separates_disagreement_from_unanswerable() {
        // Same bases, different spelling: a 3'-shifted duplication.
        let (verdict, key_in, key_out) = verify_row(&row(
            "c",
            "TTTTTTTTTAATATATTTTA",
            "NM_TEST.1:c.2_3dup",
            "NM_TEST.1:c.8_9dup",
        ));
        assert!(matches!(verdict, SpdiVerdict::Same));
        assert_eq!(key_in, key_out, "one edit, two spellings, one key");

        // Different bases. The #1513 shape: two insertions at neighbouring
        // junctions, concatenated as though they shared one.
        let (verdict, key_in, key_out) = verify_row(&row(
            "c",
            "TAAAATTATATTTATTATTT",
            "NM_TEST.1:c.[1_2insAA;2_3insTT]",
            "NM_TEST.1:c.2_3insTTAA",
        ));
        assert!(
            matches!(verdict, SpdiVerdict::Different),
            "{key_in} vs {key_out} denote different bases"
        );
        assert_ne!(key_in, key_out);

        // Unanswerable: the members overlap, so there is no single resulting
        // sequence to key on — on either side.
        let (verdict, _, _) = verify_row(&row(
            "c",
            "GCGCTAGTCTCGCCCTGTTA",
            "NM_TEST.1:c.[11_12del;11_12inv]",
            "NM_TEST.1:c.[11_12del;11_12inv]",
        ));
        assert!(matches!(verdict, SpdiVerdict::Unverifiable));

        // An unparseable side is unanswerable too, not a disagreement. `<declined>`
        // and `<panic>` are values `normalize_one` really writes into a dump.
        let (verdict, _, _) = verify_row(&row(
            "c",
            "GCGCTAGTCTCGCCCTGTTA",
            "NM_TEST.1:c.2_3dup",
            "<declined>",
        ));
        assert!(matches!(verdict, SpdiVerdict::Unverifiable));
    }

    /// The verification reads the same reference the row was normalized against.
    ///
    /// `provider_for` is shared with `normalize_one` for this reason: a second,
    /// subtly different provider would answer about the wrong sequence and would
    /// do it silently, since both still produce keys.
    #[test]
    fn verification_uses_the_rows_own_axis_provider() {
        // The multi-exon reference numbers its CDS from transcript 4, so `c.1` is
        // a different base there than on the single-exon one. Same descriptor,
        // different provider, different key.
        let core = "GCGCTAGTCTCGCCCTGTTA";
        let single = verify_row(&row("c", core, "NM_TEST.1:c.1_2dup", "NM_TEST.1:c.1_2dup"));
        let multi = verify_row(&row(
            "cx",
            core,
            "NM_TESTX.1:c.1_2dup",
            "NM_TESTX.1:c.1_2dup",
        ));
        assert!(matches!(single.0, SpdiVerdict::Same));
        assert!(matches!(multi.0, SpdiVerdict::Same));
        assert_ne!(
            single.1, multi.1,
            "the two axes place `c.1` at different transcript positions"
        );
    }

    /// Every row carries the reference **sequence** it was drawn against, which
    /// for two families is not what its `reference` column says (#1624).
    ///
    /// The two views agree on thirteen of the fifteen families, and that is
    /// precisely what made the verify pass's confusion of them survive review —
    /// a bug that only manifests on `long_block_inversion` and
    /// `separated_revcomp_runs` is a bug in 286 rows out of 78,028. So this pins
    /// the distinction itself, on both sides: where the column is a label the
    /// sequence must be the design's, and where it is not the two must be equal.
    #[test]
    fn every_row_carries_the_sequence_it_was_drawn_against() {
        let labelled: BTreeMap<String, String> = long_corpus_sequences()
            .into_iter()
            .chain(revcomp_designs().into_iter().map(|d| (d.label, d.core)))
            .collect();
        let mut label_keyed = 0usize;
        for row in dump(2) {
            // A label is not over `ACGT` — `revcomp_sep0_at0` and
            // `nearpalindrome_1024` both fail this — so it is the one predicate
            // that separates the two kinds of string without consulting the map.
            assert!(
                !row.core.is_empty() && row.core.bytes().all(|b| b"ACGT".contains(&b)),
                "[{}] {} carries something that is not a sequence: {}",
                row.family,
                row.reference,
                row.core
            );
            match labelled.get(&row.reference) {
                Some(core) => {
                    assert_eq!(
                        &row.core, core,
                        "[{}] {} carries the wrong design's sequence",
                        row.family, row.reference
                    );
                    label_keyed += 1;
                }
                None => assert_eq!(
                    row.core, row.reference,
                    "[{}] a sequence-keyed family's two views must be the same string",
                    row.family
                ),
            }
        }
        assert!(
            label_keyed > 0,
            "the label-keyed families must be in the corpus, or this test measures nothing"
        );
    }

    /// A row whose `reference` column is a label is verified against its
    /// **sequence** (#1624).
    ///
    /// `dump(0)` draws no random cores, so every row it returns belongs to one of
    /// the two families that bring their own designed references — exactly the
    /// population `--verify-spdi` could not handle. Before the fix this does not
    /// merely fail: the `cx` rows *panic* inside `spliced_genomic`, because every
    /// `separated_revcomp_runs` label is 16 bytes while the exon table maps 20,
    /// and `verify_row` builds its provider outside its own `catch_unwind`. That
    /// is why one row took the whole pass down instead of being reported as
    /// `Unverifiable`.
    ///
    /// Restricted to `separated_revcomp_runs`. The long cores are the other
    /// label-keyed family, their wiring is pinned by
    /// `every_row_carries_the_sequence_it_was_drawn_against`, and keying two
    /// canonical SPDIs off a kilobase near-palindrome per row is a cost this
    /// test does not need to pay to make its point.
    #[test]
    fn a_label_keyed_row_is_verified_against_its_core() {
        let rows: Vec<Row> = dump(0)
            .into_iter()
            .filter(|row| row.family == REVCOMP_FAMILY.0)
            .collect();
        assert!(!rows.is_empty(), "the revcomp designs must produce rows");
        let mut same = 0usize;
        for row in &rows {
            assert_ne!(
                row.reference, row.core,
                "this family keys its rows by label, so the two views must differ"
            );
            if matches!(verify_row(row).0, SpdiVerdict::Same) {
                same += 1;
            }
        }
        // Not `== rows.len()`: this family is where #1517's separation question
        // lives, and two of its three oracles are red against `main` today. What
        // this asserts is that the verification *ran and answered*, which is the
        // thing that was impossible; whether every row agrees is a different
        // finding, owned by the oracles in this module.
        assert!(
            same > 0,
            "no row verified as denoting the same bases — the pass answered nothing"
        );
    }

    /// Trap 1 from the module docs, pinned. Keying a row on the input string alone
    /// merges rows drawn against different reference sequences — which is how a
    /// measurement of 498 moved rows was reported as 6.
    #[test]
    fn the_corpus_key_is_unique() {
        let rows = dump(3);
        assert!(!rows.is_empty(), "the corpus must not be empty");
        let mut seen = std::collections::HashSet::new();
        for row in &rows {
            let key = (
                row.reference.clone(),
                row.axis,
                row.direction,
                row.input.clone(),
            );
            assert!(
                seen.insert(key.clone()),
                "duplicate corpus key {key:?} — two rows would overwrite each other in a diff"
            );
        }
        // And the input alone is NOT unique, which is precisely why the key is a
        // tuple. If this ever stops holding the trap is gone, but so is the reason
        // to trust that a future refactor kept the tuple for a reason.
        let distinct_inputs = rows
            .iter()
            .map(|r| r.input.as_str())
            .collect::<std::collections::HashSet<_>>()
            .len();
        assert!(
            distinct_inputs < rows.len(),
            "inputs repeat across references, so keying on the input alone would merge rows"
        );
    }

    /// The geometric relationship two members of a cis allele can stand in. An
    /// insertion occupies no reference base, so it is a **junction** between two
    /// bases and never a span; every other edit covers the bases it names.
    #[derive(Debug, Clone, Copy, PartialEq)]
    enum Footprint {
        Span(u64, u64),
        Junction(u64, u64),
    }

    /// Read each member's footprint off a genomic allele by parsing what the
    /// generator emitted. Deliberately derived from the description rather than
    /// from the family table: the census below has to judge what the corpus
    /// actually builds, so that renaming or adding a family cannot satisfy it
    /// vacuously. `None` for anything that is not a two-member `g.` allele.
    fn member_footprints(input: &str) -> Option<Vec<Footprint>> {
        use ferro_hgvs::hgvs::edit::NaEdit;
        use ferro_hgvs::HgvsVariant;
        let HgvsVariant::Allele(allele) = parse_hgvs(input).ok()? else {
            return None;
        };
        let mut out = Vec::new();
        for member in &allele.variants {
            let HgvsVariant::Genome(g) = member else {
                return None;
            };
            let start = g.loc_edit.location.start.inner()?.base;
            let end = g.loc_edit.location.end.inner()?.base;
            out.push(match g.loc_edit.edit.inner() {
                Some(NaEdit::Insertion { .. }) => Footprint::Junction(start, end),
                _ => Footprint::Span(start, end),
            });
        }
        Some(out)
    }

    /// The corpus must build alleles whose members **conflict** (#1456).
    ///
    /// Before this test the five shape families were all disjoint-member pairs, so
    /// the harness reported `0 of 18,432` for three consecutive PRs — #1403, #1445
    /// and #1451 — each of which moves behaviour for thousands of inputs by its own
    /// enumeration. A zero was indistinguishable from blindness, and the stability
    /// disclosure the repository requires was satisfied on paper while telling the
    /// reader nothing.
    #[test]
    fn the_corpus_emits_alleles_whose_members_conflict() {
        let (mut coincident, mut nested, mut partial, mut interior) = (0, 0, 0, 0);
        for core in corpus_sequences(2) {
            for (family, _) in FAMILIES {
                for input in inputs_for(family, "g", &core) {
                    let Some(footprints) = member_footprints(&input) else {
                        continue;
                    };
                    let &[a, b] = &footprints[..] else { continue };
                    match (a, b) {
                        (Footprint::Span(s1, e1), Footprint::Span(s2, e2)) => {
                            if (s1, e1) == (s2, e2) {
                                coincident += 1;
                            } else if (s1 <= s2 && e2 <= e1) || (s2 <= s1 && e1 <= e2) {
                                nested += 1;
                            } else if s1 <= e2 && s2 <= e1 {
                                partial += 1;
                            }
                        }
                        (Footprint::Junction(j0, j1), Footprint::Span(s, e))
                        | (Footprint::Span(s, e), Footprint::Junction(j0, j1))
                            if s <= j0 && j1 <= e =>
                        {
                            interior += 1;
                        }
                        _ => {}
                    }
                }
            }
        }
        assert!(
            coincident > 0,
            "no allele has two members on coincident bounds (`[Ndup;NX>Y]`)"
        );
        assert!(
            interior > 0,
            "no allele puts a junction interior to a span (`[a_bdel;a_binsN]`)"
        );
        assert!(
            partial > 0,
            "no allele has partially-overlapping spans (`[10_14del;12_16del]`)"
        );
        assert!(
            nested > 0,
            "no allele has nested spans (`[10_14del;10_16del]`)"
        );
    }

    /// The longest reference span any member of `input` names, in bases. A
    /// junction insertion spans nothing; every other edit spans what it names.
    fn longest_span(input: &str) -> u64 {
        use ferro_hgvs::hgvs::edit::NaEdit;
        use ferro_hgvs::HgvsVariant;
        fn of(variant: &HgvsVariant) -> u64 {
            match variant {
                HgvsVariant::Allele(a) => a.variants.iter().map(of).max().unwrap_or(0),
                HgvsVariant::Genome(g) => {
                    if matches!(g.loc_edit.edit.inner(), Some(NaEdit::Insertion { .. })) {
                        return 0;
                    }
                    match (
                        g.loc_edit.location.start.inner(),
                        g.loc_edit.location.end.inner(),
                    ) {
                        (Some(s), Some(e)) => e.base.saturating_sub(s.base) + 1,
                        _ => 0,
                    }
                }
                _ => 0,
            }
        }
        parse_hgvs(input).map(|v| of(&v)).unwrap_or(0)
    }

    /// The corpus must build a block past the normalizer's split cap (#1460).
    ///
    /// The cap is a *length* gate — `merge.rs` returns the un-partitioned whole
    /// block once either side exceeds `MAX_SPLIT_BLOCK` — so a corpus of 20-mers
    /// takes the same path on both sides of any change to it, by construction.
    /// That is why #1403 measured `0 of 18,432` here and why the zero was
    /// guaranteed rather than informative. Fixing #1456 added the shapes that were
    /// missing; this adds the scale.
    #[test]
    fn the_corpus_emits_a_block_past_the_split_cap() {
        // `MAX_SPLIT_BLOCK` is crate-private (`src/normalize/merge.rs`), so
        // the bound is restated rather than imported. If it ever moves, this test
        // is measuring the wrong threshold and should be updated with it.
        const SPLIT_CAP: u64 = 1024;
        let rows = dump(1);
        let genomic: Vec<&Row> = rows.iter().filter(|row| row.axis == "g").collect();
        let longest = genomic
            .iter()
            .map(|row| longest_span(&row.input))
            .max()
            .unwrap_or(0);
        assert!(
            longest > SPLIT_CAP,
            "the longest block the corpus builds is {longest} bases, under the {SPLIT_CAP}-base \
             split cap — so every row takes the same path on both sides of a change to it"
        );

        // Building the row is not measuring it. Every sentinel below is recorded
        // in the `output` column like any other string, so a corpus whose long
        // rows all panicked would satisfy the span assertion above and still
        // measure nothing at this scale — the same "a zero is not a result"
        // confusion (#1456) that the shape families were added to end, one level
        // down. At least one past-cap row has to have actually normalized.
        const SENTINELS: [&str; 3] = ["<parse-error>", "<declined>", "<panic>"];
        let past_cap: Vec<&&Row> = genomic
            .iter()
            .filter(|row| longest_span(&row.input) > SPLIT_CAP)
            .collect();
        let measured = past_cap
            .iter()
            .filter(|row| !SENTINELS.contains(&row.output.as_str()))
            .count();
        assert!(
            measured > 0,
            "all {} rows past the {SPLIT_CAP}-base cap normalized to a sentinel, so the corpus \
             builds the scale but measures nothing at it: {:?}",
            past_cap.len(),
            past_cap.iter().map(|r| &r.output).collect::<Vec<_>>()
        );
    }

    /// The long cores are *near* palindromes: inverting one changes exactly four
    /// bases (#1460).
    ///
    /// Both halves of that matter and both are one edit away. An **exact**
    /// palindrome inverts to itself, so the inversion denotes nothing, normalizes
    /// to `=`, and the family measures a pair of descriptions that were never a
    /// confluence question. A core that is not a palindrome at all inverts to
    /// something differing at most of its positions, so the equivalent allele
    /// carries hundreds of members and the pair stops being a spelling anyone
    /// would write. Four is the number that keeps it a real question.
    #[test]
    fn the_long_cores_are_near_palindromes() {
        for (label, core) in long_corpus_sequences() {
            let bytes = core.as_bytes();
            let inverted: Vec<u8> = bytes.iter().rev().map(|b| complement(*b)).collect();
            let differing = (0..bytes.len())
                .filter(|&i| bytes[i] != inverted[i])
                .count();
            assert_eq!(
                differing, 4,
                "{label} inverts with {differing} base changes, not 4 — the two spellings \
                 this family compares are only a confluence question at a handful"
            );
        }
    }

    /// A reduced seed count must be a strict prefix of a larger one, so a quick
    /// local run and a full one are comparable rather than merely similar.
    #[test]
    fn the_sequence_corpus_is_prefix_stable() {
        let few = corpus_sequences(3);
        let many = corpus_sequences(9);
        assert_eq!(few, many[..few.len()]);
        assert_eq!(few.len(), 6, "two sequences per seed");
    }

    /// Every family emits something on every axis, so a family cannot silently
    /// contribute zero rows and make its "0 moved" look like evidence.
    #[test]
    fn every_family_emits_rows_on_both_axes() {
        let core = &corpus_sequences(1)[0];
        for (family, _) in FAMILIES {
            for axis in ["g", "c", "cx"] {
                assert!(
                    !inputs_for(family, axis, core).is_empty(),
                    "family {family} emitted no {axis}. rows"
                );
            }
        }
    }

    /// The corpus feeds back every construct it emits.
    ///
    /// Measured on `5a19d1d2`, before the four closing families existed, the corpus
    /// *produced* repeat notation on 1,558 output rows, a no-change `=` on 152 and
    /// three-member alleles on 170 — while feeding **none** of them in as an input.
    /// #1454 was one consequence of that asymmetry, and it was found by accident
    /// rather than by the corpus.
    ///
    /// `inv` is checked here too even though #1459 already feeds it, because all
    /// 3,072 of those inputs pair it with an overlapping insertion. A count alone
    /// would call that construct covered; it is the *lone* `inv` this asserts.
    ///
    /// This is the guard against the asymmetry reappearing, and it works in two
    /// halves. The four literal assertions pin that each new family exists at all.
    /// The loop after them is the closure check proper: it reads the constructs out
    /// of the **dump's own outputs** and requires an input counterpart for each, so
    /// a family added later that emits one of these without feeding it back fails
    /// here. What it still cannot see is a construct outside the vocabulary below —
    /// that vocabulary is a hand audit, and widening it is the manual step.
    #[test]
    fn the_corpus_is_closed_over_its_output_vocabulary() {
        let core = &corpus_sequences(1)[0];
        let inputs: Vec<String> = FAMILIES
            .iter()
            .flat_map(|(family, _)| {
                ["g", "c"]
                    .into_iter()
                    .flat_map(move |axis| inputs_for(family, axis, core))
            })
            .collect();
        let seen = |pred: &dyn Fn(&str) -> bool| inputs.iter().any(|i| pred(i));

        // A *lone* `inv`, not merely any `inv`. #1459's `junction_interior_to_span`
        // emits 3,072 inv inputs and every one pairs the inversion with an
        // overlapping insertion, so a bare `contains("inv")` here would pass with
        // `inv_member` deleted — asserting nothing this family provides.
        assert!(
            seen(&|i| i.contains("inv") && !i.contains(';')),
            "no input exercises a lone `inv`; only conflicting pairings are covered"
        );
        assert!(
            seen(&|i| i.contains("[3]") || i.contains("[1]")),
            "no input exercises repeat notation"
        );
        // `=` can sit at either position in the allele, so match both delimiters
        // rather than only the closing bracket.
        assert!(
            seen(&|i| i.contains("=;") || i.contains("=]")),
            "no input exercises `=`"
        );
        assert!(
            seen(&|i| i.matches(';').count() >= 2),
            "no input exercises a three-member allele"
        );

        // The four assertions above pin that the new families exist. They do not,
        // on their own, make the test live up to its name: the constructs are
        // named here as literals, so a family that started *emitting* one without
        // feeding it back would not move them. That is the asymmetry this module
        // is about, one level up — a hardcoded audit going stale the same way the
        // original one did.
        //
        // So close the loop against the dump itself: for each construct in the
        // vocabulary, if any normalized **output** exhibits it, some **input**
        // must too. This is what fails when a later family widens the output
        // vocabulary and leaves the input side behind.
        /// A named construct predicate, applied to both sides of the dump.
        type Construct<'a> = (&'a str, &'a dyn Fn(&str) -> bool);

        let rows = dump(1);
        let vocabulary: &[Construct] = &[
            ("lone `inv`", &|s: &str| {
                s.contains("inv") && !s.contains(';')
            }),
            ("repeat notation", &|s: &str| {
                s.contains("[3]") || s.contains("[1]")
            }),
            ("`=`", &|s: &str| s.contains("=;") || s.contains("=]")),
            ("three-member allele", &|s: &str| {
                s.matches(';').count() >= 2
            }),
        ];
        for (name, exhibits) in vocabulary {
            let in_output = rows.iter().any(|row| exhibits(&row.output));
            let in_input = rows.iter().any(|row| exhibits(&row.input));
            assert!(
                !in_output || in_input,
                "the corpus emits {name} as an output but never feeds one in, so it cannot \
                 measure what normalization does to that form — add a family that does"
            );
        }
    }

    /// `repeat_expansion` describes a tract the reference actually has.
    ///
    /// The alternative — anchoring a repeat on a single copy — is spec-questionable
    /// (`repeated.md` defines a repeat as a unit "present several times") and would
    /// make the family assert something about the reference that is not true. The
    /// consequence is partial coverage: only 21 of the 48 corpus sequences contain a
    /// tandem triplet, so the family is silent on the rest. That trade is deliberate.
    #[test]
    fn the_repeat_family_only_describes_a_tract_the_reference_has() {
        let mut covered = 0;
        for core in corpus_sequences(24) {
            let rows = inputs_for("repeat_expansion", "g", &core);
            match first_tandem_triplet(&core) {
                None => assert!(
                    rows.is_empty(),
                    "emitted a repeat for {core}, which has no tandem triplet"
                ),
                Some(t) => {
                    covered += 1;
                    assert_eq!(rows.len(), 2, "expected both counts for {core}");
                    let unit = &core[t..t + 3];
                    assert_eq!(
                        &core[t + 3..t + 6],
                        unit,
                        "{unit} at {t} in {core} is not actually tandem"
                    );
                    assert!(
                        rows[0].ends_with(&format!("{unit}[3]")),
                        "unexpected spelling {}",
                        rows[0]
                    );
                }
            }
        }
        assert_eq!(
            covered, 21,
            "tandem-triplet coverage moved; the generator or seed count changed"
        );
    }

    /// `delins_hiding_an_inversion` emits both spellings, with the trailing piece an
    /// exact reverse complement of its own span.
    ///
    /// Pinned as exact strings rather than as a predicate: a predicate over "some
    /// member is an inversion" is satisfied by a corpus that emits the spanning form
    /// alone, which is precisely the degenerate case this family exists to avoid.
    ///
    /// The head payload differs between the two cores below, and that is the point.
    /// For `TTT` the first candidate `GGG` is neither the span nor its reverse
    /// complement (`AAA`), so it is taken. For `CCC` the span rules out `CCC` and its
    /// reverse complement rules out `GGG`, so `head_payload` falls through to `ACA`.
    #[test]
    fn the_inversion_family_emits_both_spellings_with_a_reverse_complement_tail() {
        let homopolymer = "TTTTTTTTTAATATATTTTA";
        assert_eq!(
            inputs_for("delins_hiding_an_inversion", "g", homopolymer)[..2],
            [
                "NC_TEST.1:g.257_262delinsGGGAAA".to_string(),
                "NC_TEST.1:g.[257_259delinsGGG;260_262delinsAAA]".to_string(),
            ]
        );

        let four_letter = "CCCCCCCCTGACGTATCCTA";
        assert_eq!(
            inputs_for("delins_hiding_an_inversion", "g", four_letter)[..2],
            [
                "NC_TEST.1:g.257_262delinsACAGGG".to_string(),
                "NC_TEST.1:g.[257_259delinsACA;260_262delinsGGG]".to_string(),
            ]
        );
    }

    /// The pieces are three bases each, and that width is load-bearing rather than
    /// arbitrary.
    ///
    /// Measured on `ee0e37ac`: a two-base head with a three-base tail produces **0**
    /// non-idempotent rows across the whole corpus, because the trailing piece is
    /// typed `inv` on the first pass. Widening the head to three bases produces
    /// **30**, all of them #1454's shape — a piece spelled `delins` on pass 1 and
    /// re-typed `inv` on pass 2. So a future edit that narrows either piece silently
    /// returns the corpus to being blind to the defect it was added for.
    #[test]
    fn the_inversion_family_uses_three_base_pieces() {
        let core = &corpus_sequences(1)[0];
        let split = inputs_for("delins_hiding_an_inversion", "g", core)
            .into_iter()
            .find(|s| s.contains('['))
            .expect("the family emits a pre-split spelling");
        let (head, tail) = split
            .split_once(';')
            .expect("the pre-split spelling has two members");
        for (member, payload) in [(head, "head"), (tail, "tail")] {
            let emitted = member
                .rsplit_once("delins")
                .expect("each member is a delins")
                .1
                .trim_end_matches(']');
            assert_eq!(
                emitted.len(),
                3,
                "the {payload} piece must be three bases, found {emitted:?} in {split}"
            );
        }
    }

    /// `head_payload` never degenerates: it is neither a no-op nor a second
    /// inversion. Exhaustive over all 64 three-base spans, so the guarantee does not
    /// rest on the corpus happening not to draw a colliding one.
    #[test]
    fn the_head_payload_is_never_a_no_op_or_an_inversion() {
        for a in "ACGT".chars() {
            for b in "ACGT".chars() {
                for c in "ACGT".chars() {
                    let span: String = [a, b, c].into_iter().collect();
                    let payload = head_payload(&span);
                    assert_ne!(payload, span, "head payload is a no-op on {span}");
                    assert_ne!(
                        payload,
                        revcomp(&span),
                        "head payload is itself an inversion on {span}"
                    );
                }
            }
        }
    }

    /// An odd-length span can never equal its own reverse complement, which is why
    /// the tail piece is three bases and not two: a two-base tail degenerates to a
    /// no-op on every palindromic pair (`AT`, `TA`, `CG`, `GC`).
    #[test]
    fn a_three_base_span_never_equals_its_own_reverse_complement() {
        for a in "ACGT".chars() {
            for b in "ACGT".chars() {
                for c in "ACGT".chars() {
                    let span: String = [a, b, c].into_iter().collect();
                    assert_ne!(revcomp(&span), span, "{span} is its own reverse complement");
                }
            }
        }
        for pair in ["AT", "TA", "CG", "GC"] {
            assert_eq!(revcomp(pair), pair, "{pair} should be palindromic");
        }
    }

    /// The coding corpus must reach past `CDS_END` and spell those positions `c.*N`.
    ///
    /// Pinned because the first version of this file stopped exactly at `CDS_END`, so
    /// it emitted no junction-crossing row at all — and would have reported a
    /// confident "0 moved" for #1417, #1425, #1427 and #1428, every one of which is a
    /// CDS/3'UTR junction defect.
    #[test]
    fn the_coding_corpus_crosses_into_the_three_prime_utr() {
        let core = &corpus_sequences(1)[0];
        let mut utr_rows = 0;
        let mut straddling_rows = 0;
        for (family, _) in FAMILIES {
            for input in inputs_for(family, "c", core) {
                if input.contains(":c.*") || input.contains(";*") || input.contains("_*") {
                    utr_rows += 1;
                    // A shape that names both a CDS position and a `*` position is
                    // the junction-crossing case specifically.
                    if input.contains("_*") || input.contains(";*") {
                        straddling_rows += 1;
                    }
                }
            }
        }
        assert!(
            utr_rows > 0,
            "the coding corpus emitted no 3'UTR (`c.*N`) row; CDS_END={CDS_END} against a \
             {}-base core leaves no room, so junction defects would measure as zero",
            core.len()
        );
        assert!(
            straddling_rows > 0,
            "the coding corpus emitted 3'UTR rows but none straddling CDS_END, which is \
             the shape four of this campaign's five defects lived in"
        );
    }

    /// The multi-exon axis reaches the 5'UTR, and the single-exon axis cannot.
    ///
    /// Both halves matter. The first is the gap #1478 closes; the second is why it
    /// needed a second reference rather than a wider loop over the existing one —
    /// with `CDS_START == 1` there is no transcript position before the CDS, so no
    /// family, however written, can emit a `c.-N` row against it.
    #[test]
    fn only_the_multi_exon_axis_reaches_the_five_prime_utr() {
        let core = &corpus_sequences(1)[0];
        let names_five_prime_utr = |axis: &str| {
            FAMILIES.iter().any(|(family, _)| {
                inputs_for(family, axis, core)
                    .iter()
                    .any(|i| i.contains(".-") || i.contains(";-") || i.contains("_-"))
            })
        };
        assert!(
            names_five_prime_utr("cx"),
            "the multi-exon axis emitted no 5'UTR (`c.-N`) row, so #1478's gap is not closed"
        );
        assert!(
            !names_five_prime_utr("c"),
            "the single-exon axis emitted a 5'UTR row, which CDS_START=1 makes impossible; \
             the axis plumbing has drifted"
        );
    }

    /// A `c.` position before the CDS is spelled `-N`, and the general form is a
    /// no-op on the single-exon reference.
    ///
    /// The second assertion is the guard that this change did not move any existing
    /// row: with `CDS_START == 1` the generalised `hgvs_pos` must still spell
    /// transcript position `n` as plain `n`.
    #[test]
    fn coding_positions_before_the_cds_are_spelled_with_a_minus() {
        assert_eq!(hgvs_pos("cx", 0), "-3");
        assert_eq!(hgvs_pos("cx", 1), "-2");
        assert_eq!(hgvs_pos("cx", 2), "-1");
        assert_eq!(hgvs_pos("cx", 3), "1");
        assert_eq!(hgvs_pos("cx", CDS_END_MULTI as usize - 1), "12");
        assert_eq!(hgvs_pos("cx", CDS_END_MULTI as usize), "*1");

        for i in 0..(CDS_END as usize) {
            assert_eq!(
                hgvs_pos("c", i),
                (i + 1).to_string(),
                "generalising hgvs_pos moved a single-exon coding position"
            );
        }
    }

    /// The multi-exon reference really is spliced: three exons, two junctions, and
    /// genomic spans separated by intron rather than abutting.
    ///
    /// Without the gap the exons would be contiguous in genomic space and the
    /// transcript would be a single-exon one wearing three exon records — which
    /// would look right in the dump and still fail to express #1450.
    #[test]
    fn the_multi_exon_reference_is_actually_spliced() {
        assert!(
            EXON_SPANS.len() >= 3,
            "need at least two exon/exon junctions"
        );
        let core = &corpus_sequences(1)[0];
        let genomic = spliced_genomic(core);

        for (n, (from, to)) in EXON_SPANS.iter().copied().enumerate() {
            let g = exon_genomic_start(n) as usize;
            assert_eq!(
                &genomic[g - 1..g - 1 + (to - from)],
                &core[from..to],
                "exon {n} is not at the genomic offset its Exon record claims"
            );
            if n + 1 < EXON_SPANS.len() {
                let end = g + (to - from) - 1;
                assert_eq!(
                    exon_genomic_start(n + 1) as usize - end - 1,
                    INTRON_LEN,
                    "exons {n} and {} abut; there is no intron between them",
                    n + 1
                );
            }
        }

        let spliced: String = EXON_SPANS.iter().map(|(f, t)| &core[*f..*t]).collect();
        assert_eq!(
            spliced, *core,
            "the exons do not reconstitute the transcript"
        );
    }

    /// A `c.` position past `CDS_END` is `*N`, and one inside it is plain.
    #[test]
    fn coding_positions_past_the_cds_are_spelled_with_a_star() {
        // Offset i maps to 1-based position i+1.
        assert_eq!(hgvs_pos("c", 0), "1");
        assert_eq!(hgvs_pos("c", (CDS_END - 1) as usize), CDS_END.to_string());
        assert_eq!(hgvs_pos("c", CDS_END as usize), "*1");
        assert_eq!(hgvs_pos("c", CDS_END as usize + 1), "*2");
        // The genomic axis has no such notion and is offset by the pad.
        assert_eq!(hgvs_pos("g", 0), (PAD_OFFSET + 1).to_string());
    }

    /// The header is matched exactly, so reordered or renamed columns cannot be read
    /// positionally under their old meanings.
    #[test]
    fn the_reader_rejects_a_wrong_header() {
        let path = fixture("bad-header.tsv", "reference\taxis\tdirection\tfamily\tinput\twas_fixed_point\toutput\nAC\tg\t3prime\tf\tX\tfalse\tY\n");
        let err = read_dump(&path).expect_err("a reordered header must be rejected");
        assert!(err.contains("unexpected header"), "unexpected error: {err}");
    }

    /// `was_fixed_point` is parsed strictly. Anything else would understate the
    /// migration, which is the direction that gets a change waved through.
    #[test]
    fn the_reader_rejects_a_malformed_boolean() {
        let path = fixture(
            "bad-bool.tsv",
            &format!("{HEADER}AC\tg\t3prime\tf\tX\tY\tyes\n"),
        );
        let err = read_dump(&path).expect_err("a non-boolean must be rejected");
        assert!(
            err.contains("must be `true` or `false`"),
            "unexpected error: {err}"
        );
    }

    /// Two dumps that disagree about a key's shape family are refused, rather than
    /// producing a per-family table that attributes movement to the wrong shape.
    #[test]
    fn comparing_dumps_that_disagree_on_a_family_is_refused() {
        let before = fixture(
            "fam-before.tsv",
            &format!("{HEADER}AC\tg\t3prime\tfamily_one\tX\tY\ttrue\n"),
        );
        let after = fixture(
            "fam-after.tsv",
            &format!("{HEADER}AC\tg\t3prime\tfamily_two\tX\tY\ttrue\n"),
        );
        let err = compare(&before, &after).expect_err("a family mismatch must be refused");
        assert!(
            err.contains("disagree on the shape family"),
            "unexpected error: {err}"
        );
    }

    /// A zero has to say what it covered (#1456).
    ///
    /// `0 moved` reads as "this change is safe", but it is only ever a statement
    /// about the families the corpus builds — three PRs in a row quoted a zero from
    /// a corpus that could not express the shape they changed. Naming the covered
    /// families in the zero case is what lets a reader tell the two apart.
    #[test]
    fn a_report_with_no_movement_names_the_families_it_covered() {
        let dump = fixture(
            "zero-scope.tsv",
            &format!(
                "{HEADER}AC\tg\t3prime\tfamily_one\tX\tY\ttrue\n\
                 AC\tg\t3prime\tfamily_two\tP\tQ\ttrue\n"
            ),
        );
        let report = compare(&dump, &dump).expect("a self-comparison must succeed");
        assert!(
            report.contains("No row moved"),
            "a zero must be stated, not left to be read off the table:\n{report}"
        );
        assert!(
            report.contains("does not build"),
            "a zero must say it covers only the families in the corpus:\n{report}"
        );
        for family in ["family_one", "family_two"] {
            assert!(
                report.contains(&format!("`{family}`")),
                "the zero caveat must name {family}:\n{report}"
            );
        }
    }

    /// …and a report that *did* find movement must not carry that claim. The
    /// discriminating half: printing the caveat unconditionally would satisfy the
    /// test above while telling a reader of a moving diff that nothing moved.
    #[test]
    fn a_report_with_movement_does_not_claim_a_zero() {
        let before = fixture(
            "moved-before.tsv",
            &format!("{HEADER}AC\tg\t3prime\tfamily_one\tX\tY\ttrue\n"),
        );
        let after = fixture(
            "moved-after.tsv",
            &format!("{HEADER}AC\tg\t3prime\tfamily_one\tX\tZ\ttrue\n"),
        );
        let report = compare(&before, &after).expect("comparison must succeed");
        assert!(
            !report.contains("No row moved"),
            "a report with a moved row must not claim a zero:\n{report}"
        );
    }

    // ----------------------------------------------------------------------
    // `separated_revcomp_runs` — the family, and the three oracles it enables
    // ----------------------------------------------------------------------
    //
    // Everything above measures *movement*: two dumps, one diff, "did this change
    // the output". Nothing above ever asks whether an output is **right**, and for
    // thirteen of the fifteen families it cannot — they are sets of descriptions,
    // with no record of what the description was meant to denote.
    //
    // `separated_revcomp_runs` is built the other way round: the generator picks
    // the alternate first (`revcomp(span)`) and derives the reference around it, so
    // `(reference, alternate)` is known exactly for every row. That makes three
    // real oracles available for free, and each of them can catch a class of defect
    // a movement count cannot see:
    //
    // 1. **apply** — the output denotes the intended bases. #1416, #1453 and #1431
    //    were each an output that denoted the wrong sequence or none at all, and
    //    each was found by hand.
    // 2. **confluence** — the three spellings of one design agree, and their common
    //    answer is a fixed point.
    // 3. **separation** — the members obey `general.md:34`/`:35`, checked against
    //    the sequence rather than against a pinned string.

    /// Run `f` over every row `separated_revcomp_runs` contributes to a dump.
    ///
    /// Deliberately re-derives the rows rather than filtering [`dump`]: the oracles
    /// need the *design* behind each row — its span, its intended alternate, its run
    /// layout — and a `Row` carries only a label.
    fn for_each_revcomp_row(
        mut f: impl FnMut(&RevcompDesign, &'static str, &'static str, &str, &str),
    ) {
        for design in revcomp_designs() {
            for (axis, direction, dir_label) in axes_and_directions() {
                for input in revcomp_inputs_for(axis, &design) {
                    let output = normalize_one(axis, &design.core, &input, direction);
                    f(&design, axis, dir_label, &input, &output);
                }
            }
        }
    }

    /// The designed blocks really do have the run structure their parameters claim.
    ///
    /// Two views of the same thing are pinned against each other. [`revcomp_runs`]
    /// derives the runs from `(width, separation)`; this re-derives them from the
    /// **sequence**, by comparing the block with its own reverse complement column
    /// by column. A typo in one of the five literals in [`revcomp_span`] moves the
    /// second view and not the first, so it fails here rather than silently
    /// contributing a design at a separation nobody asked for.
    ///
    /// Separation 0 is the one case where the two views legitimately differ: two
    /// adjacent 2-column runs are indistinguishable *as a sequence* from one
    /// 4-column run, which is exactly what "adjacent changes are one description"
    /// means. So the sequence view is checked against the **union** of the runs,
    /// and the count is checked only where a gap separates them.
    #[test]
    fn the_designed_blocks_have_the_coincidence_structure_they_claim() {
        for &separation in REVCOMP_SEPARATIONS {
            let span = revcomp_span(separation);
            let alternate = revcomp(span);
            let runs = revcomp_runs(separation);

            let changed: Vec<usize> = (0..span.len())
                .filter(|&i| span.as_bytes()[i] != alternate.as_bytes()[i])
                .collect();
            let designed: Vec<usize> = runs
                .iter()
                .flat_map(|&(start, len)| start..start + len)
                .collect();
            assert_eq!(
                changed, designed,
                "separation-{separation} block {span} inverts to {alternate}, which changes \
                 columns {changed:?} — not the {designed:?} its run layout claims"
            );
            assert!(
                runs.iter().all(|&(_, len)| len >= REVCOMP_RUN_WIDTH),
                "separation-{separation} has a run narrower than {REVCOMP_RUN_WIDTH} columns, \
                 so it is a substitution shape and not the multi-column one this family adds"
            );
            for pair in runs.windows(2) {
                let [(start, len), (next, _)] = [pair[0], pair[1]];
                assert_eq!(
                    next - (start + len),
                    separation,
                    "runs of the separation-{separation} block are {} columns apart",
                    next - (start + len)
                );
            }
            assert_eq!(
                runs.last().map(|&(s, l)| s + l),
                Some(span.len()),
                "the separation-{separation} block has unchanged columns outside its outermost \
                 runs, so its stated span is not minimal"
            );
        }
    }

    /// The worked case is in the corpus, verbatim.
    ///
    /// `AATGCACA` reverse complements to `TGTGCATT`: a true inversion that coincides
    /// with its own reference at its four interior columns, leaving two 2-base
    /// changes separated by four unchanged ones. It is the shape the family exists
    /// for, so it is pinned as a literal rather than left to be implied by the
    /// parameters.
    #[test]
    fn the_worked_case_is_the_separation_four_design() {
        assert_eq!(revcomp_span(4), "AATGCACA");
        assert_eq!(revcomp("AATGCACA"), "TGTGCATT");
        assert_eq!(revcomp_runs(4), vec![(0, 2), (6, 2)]);
        // …and it reaches the dump, with both flanking runs multi-column.
        let design = revcomp_designs()
            .into_iter()
            .find(|d| d.separation == 4 && d.offset == 0)
            .expect("the separation-4 design exists");
        assert!(design.core.starts_with("AATGCACA"));
        assert_eq!(design.alternate, "TGTGCATT");
        assert_eq!(
            revcomp_inputs_for("g", &design),
            vec![
                "NC_TEST.1:g.257_264inv".to_string(),
                "NC_TEST.1:g.257_264delinsTGTGCATT".to_string(),
                "NC_TEST.1:g.[257_258delinsTG;263_264delinsTT]".to_string(),
            ]
        );
    }

    /// An odd separation needs three runs, and that is a fact about reverse
    /// complements rather than a quirk of the five literals.
    ///
    /// Exhaustive over every block of length <= 8 with exactly two changed runs:
    /// none of them has an odd gap. If this ever finds one, the argument in
    /// [`revcomp_span`]'s docs is wrong and the three-run layout it justifies is
    /// unnecessary.
    #[test]
    fn two_changed_runs_can_never_be_an_odd_number_of_columns_apart() {
        let mut odd_gaps_seen: Vec<String> = Vec::new();
        let mut spans: Vec<Vec<u8>> = vec![Vec::new()];
        for _ in 0..8 {
            spans = spans
                .iter()
                .flat_map(|span| {
                    b"ACGT".iter().map(|base| {
                        let mut wider = span.clone();
                        wider.push(*base);
                        wider
                    })
                })
                .collect();
            for span in &spans {
                let text = std::str::from_utf8(span).expect("ACGT is valid UTF-8");
                let alternate = revcomp(text);
                let changed: Vec<bool> = span
                    .iter()
                    .zip(alternate.as_bytes())
                    .map(|(r, a)| r != a)
                    .collect();
                let runs = changed_runs(&changed);
                if runs.len() != 2 {
                    continue;
                }
                let gap = runs[1].0 - (runs[0].0 + runs[0].1);
                if !gap.is_multiple_of(2) {
                    odd_gaps_seen.push(format!("{text} -> {alternate} (gap {gap})"));
                }
            }
        }
        assert!(
            odd_gaps_seen.is_empty(),
            "a two-run block with an odd gap exists, so `revcomp_span`'s parity argument — \
             and the three-run layout it justifies — is wrong: {odd_gaps_seen:?}"
        );
    }

    /// Maximal runs of `true`, as `(offset, length)` — edges included.
    ///
    /// The sibling of [`interior_runs`], which drops the edge runs on purpose. A
    /// *changed* run at the edge of a block is ordinary (a minimal block starts and
    /// ends changed), whereas an *unchanged* run at the edge is an untrimmed flank,
    /// which is why the two functions differ.
    fn changed_runs(flags: &[bool]) -> Vec<(usize, usize)> {
        let mut runs = Vec::new();
        let mut i = 0;
        while i < flags.len() {
            if !flags[i] {
                i += 1;
                continue;
            }
            let start = i;
            while i < flags.len() && flags[i] {
                i += 1;
            }
            runs.push((start, i - start));
        }
        runs
    }

    /// Gaps between **multi-column** changed runs of every reverse-complement block
    /// `input` names, in unchanged reference columns.
    ///
    /// A member counts as a reverse-complement block when its inserted bases are the
    /// reverse complement of the bases it replaces — however it is spelled, so an
    /// `inv` and the `delins` that says the same thing both qualify. Only gaps whose
    /// *both* flanking runs are at least [`REVCOMP_RUN_WIDTH`] columns are reported,
    /// which is precisely the shape a gate keyed on "multi-column pieces at a
    /// separation" engages on.
    ///
    /// A block with one changed run contributes nothing, so the separation-0 design
    /// does not appear here — adjacency *is* one run, which is the content of the
    /// merge rule rather than a gap in the measurement.
    fn separated_revcomp_gaps(input: &str, core: &str) -> Vec<usize> {
        let Ok(variant) = parse_hgvs(input) else {
            return Vec::new();
        };
        let Ok(facts) = member_facts_of(&variant, &genomic_provider(core)) else {
            return Vec::new();
        };
        let mut gaps = Vec::new();
        for member in &facts {
            if member.reference.is_empty() || member.alternate != revcomp(&member.reference) {
                continue;
            }
            let changed: Vec<bool> = member
                .reference
                .bytes()
                .zip(member.alternate.bytes())
                .map(|(r, a)| r != a)
                .collect();
            gaps.extend(changed_runs(&changed).windows(2).filter_map(|pair| {
                let ((start, len), (next, next_len)) = (pair[0], pair[1]);
                (len >= REVCOMP_RUN_WIDTH && next_len >= REVCOMP_RUN_WIDTH)
                    .then(|| next - (start + len))
            }));
        }
        gaps
    }

    /// No other family builds a reverse-complement block whose multi-column changed
    /// runs are separated at all.
    ///
    /// This is the structural claim the family exists to fix, asserted rather than
    /// argued. `delins_hiding_an_inversion` is the only other family that builds a
    /// reverse-complement piece and its pieces are three columns wide and adjacent;
    /// `long_block_inversion` builds a kilobase one, but its perturbations are
    /// *single*-column runs. So a gate keyed on multi-column pieces at a separation
    /// engaged on no row of this corpus, and a change to one measured `0 of 78,028`
    /// for structural reasons rather than safe ones.
    ///
    /// The second half is what keeps the first honest: the new family really does
    /// reach every separation it claims. Separation 0 is excluded there because it
    /// is not a gap — two adjacent runs are one run — and that is the merge case,
    /// not a hole.
    #[test]
    fn only_this_family_separates_the_runs_of_a_reverse_complement_block() {
        let core = &corpus_sequences(1)[0];
        for (family, _) in FAMILIES {
            for input in inputs_for(family, "g", core) {
                assert!(
                    separated_revcomp_gaps(&input, core).is_empty(),
                    "family {family} already builds a separated reverse-complement block \
                     ({input}), so the structural gap this family was added for has moved"
                );
            }
        }
        for (_, long_core) in long_corpus_sequences() {
            for input in long_inputs_for("g", &long_core) {
                assert!(
                    separated_revcomp_gaps(&input, &long_core).is_empty(),
                    "{} already builds a separated reverse-complement block",
                    LONG_FAMILY.0
                );
            }
        }

        let reached: std::collections::BTreeSet<usize> = revcomp_designs()
            .iter()
            .flat_map(|design| {
                revcomp_inputs_for("g", design)
                    .into_iter()
                    .flat_map(|input| separated_revcomp_gaps(&input, &design.core))
            })
            .collect();
        assert_eq!(
            reached,
            REVCOMP_SEPARATIONS
                .iter()
                .copied()
                .filter(|separation| *separation > 0)
                .collect::<std::collections::BTreeSet<_>>(),
            "the family does not reach every separation it claims"
        );
        // …and separation 0 is present as the one shape that has no gap at all.
        let merged = revcomp_designs()
            .into_iter()
            .find(|design| design.separation == 0)
            .expect("the separation-0 design exists");
        let changed: Vec<bool> = merged
            .span
            .bytes()
            .zip(merged.alternate.bytes())
            .map(|(r, a)| r != a)
            .collect();
        assert_eq!(
            changed_runs(&changed),
            vec![(0, 2 * REVCOMP_RUN_WIDTH)],
            "the separation-0 design must present as one run, which is what adjacency means"
        );
    }

    // -- Oracle 1: sequence preservation -----------------------------------

    /// Where a design's core sits inside the sequence `accession` serves.
    ///
    /// `None` for the multi-exon *contig*, where the core is interleaved with
    /// introns and so is not one contiguous run — a single offset cannot describe
    /// it, and guessing one would compare the wrong bases while still producing a
    /// verdict.
    fn reference_frame(accession: &str, core: &str) -> Option<(String, usize)> {
        match accession {
            GENOMIC_CONTIG | TX_CONTIG => Some((padded(core), PAD_OFFSET)),
            TX_ACCESSION | TX_MULTI_ACCESSION => Some((core.to_string(), 0)),
            _ => None,
        }
    }

    /// What the apply oracle concluded about one row.
    #[derive(Debug)]
    enum ApplyVerdict {
        /// Applying the normalized output to the reference reproduces the alternate
        /// the generator intended.
        Preserved,
        /// It reproduces something else — the defect class this oracle exists for.
        Wrong { got: String, want: String },
        /// It could not be applied, so the question has no answer here. Not a
        /// failure: `apply_to_reference` declines an allele whose members overlap,
        /// and a normalized output can be a sentinel like `<declined>`.
        Unapplicable(String),
    }

    /// Apply `output` to the design's reference and compare with `design.alternate`.
    ///
    /// This is the oracle a movement count cannot be: it asks what the output
    /// *means*, against a ground truth the generator wrote down before the
    /// normalizer ever ran. Comparing the output against its own input instead —
    /// which is what `--verify-spdi` does — cannot catch a spelling that was already
    /// wrong on the way in, and cannot say which of two disagreeing spellings is
    /// the wrong one.
    fn apply_verdict(axis: &str, design: &RevcompDesign, output: &str) -> ApplyVerdict {
        let Ok(variant) = parse_hgvs(output) else {
            return ApplyVerdict::Unapplicable(format!("output does not parse: {output}"));
        };
        let provider = provider_for(axis, &design.core);
        let applied = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            ferro_hgvs::spdi::apply_to_reference(&variant, &provider)
        })) {
            Ok(Ok(applied)) => applied,
            Ok(Err(e)) => return ApplyVerdict::Unapplicable(e.to_string()),
            Err(_) => return ApplyVerdict::Unapplicable("panicked while applying".to_string()),
        };
        let Some((sequence, base)) = reference_frame(&applied.accession, &design.core) else {
            return ApplyVerdict::Unapplicable(format!(
                "no contiguous frame for {}",
                applied.accession
            ));
        };
        let start = applied.start as usize;
        let stop = start + applied.reference.len();
        // A window that does not match the bases this frame serves means the two
        // sides are talking about different sequences, and any verdict drawn from
        // them would be about the wrong locus. Say so rather than compare.
        if stop > sequence.len() || sequence[start..stop] != applied.reference {
            return ApplyVerdict::Unapplicable(format!(
                "the applied window {start}..{stop} does not match the reference frame"
            ));
        }
        let got = format!(
            "{}{}{}",
            &sequence[..start],
            applied.resulting,
            &sequence[stop..]
        );
        let block = base + design.offset;
        let want = format!(
            "{}{}{}",
            &sequence[..block],
            design.alternate,
            &sequence[block + design.span.len()..]
        );
        if got == want {
            ApplyVerdict::Preserved
        } else {
            ApplyVerdict::Wrong { got, want }
        }
    }

    /// **Oracle 1.** Every normalized output still denotes the bases the generator
    /// intended.
    ///
    /// An output that denotes a different sequence — or none — is invisible to a
    /// movement count, to the idempotency oracle (a wrong form can be a perfectly
    /// stable fixed point) and to the re-parse oracle (a wrong form parses fine).
    /// It has produced at least three filed defects here (#1416, #1453, #1431),
    /// each found by hand. This is the check that sees it.
    ///
    /// **Green on the first run**: 270 of 270 rows applied and every one produced
    /// the alternate the generator built; nothing was declined. That is a real
    /// result rather than a vacuous one — the `preserved > 0` assertion below is
    /// what stops a run where everything declined from reading as a clean bill of
    /// health.
    #[test]
    fn every_normalized_output_denotes_the_intended_alternate() {
        let mut preserved = 0usize;
        let mut declined: Vec<String> = Vec::new();
        let mut wrong: Vec<String> = Vec::new();
        for_each_revcomp_row(|design, axis, direction, input, output| {
            match apply_verdict(axis, design, output) {
                ApplyVerdict::Preserved => preserved += 1,
                // Kept, not merely counted: a run where *everything* declined
                // would otherwise report a clean bill of health, and the reason
                // is the only thing that distinguishes the two.
                ApplyVerdict::Unapplicable(reason) => declined.push(reason),
                ApplyVerdict::Wrong { got, want } => wrong.push(format!(
                    "[{} {} {}] {input}\n      -> {output}\n      got  {got}\n      want {want}",
                    design.label, axis, direction
                )),
            }
        });
        assert!(
            wrong.is_empty(),
            "{} normalized outputs denote bases other than the alternate the generator built \
             ({preserved} preserved, {} unapplicable):\n{}",
            wrong.len(),
            declined.len(),
            wrong.join("\n")
        );
        // A corpus where nothing could be applied would pass the assertion above
        // while measuring nothing — the same "a zero is not a result" confusion the
        // shape families were added to end, one level down.
        assert!(
            preserved > 0,
            "no row could be applied at all, so this oracle measured nothing; the first \
             declines were: {:?}",
            declined.iter().take(5).collect::<Vec<_>>()
        );
    }

    // -- Oracle 2: confluence ----------------------------------------------

    /// **Oracle 2.** The three spellings of one design normalize to one string, and
    /// that string is a fixed point.
    ///
    /// Both halves are needed and neither implies the other. Agreement alone is
    /// satisfied by three spellings that all normalize to a form the normalizer
    /// would move again on the next pass; a fixed point alone is satisfied by three
    /// *different* stable answers. Asserting the pair is what makes this a
    /// confluence property rather than a table of expected strings — nothing here
    /// says which of the three forms should win.
    /// **Currently red, deliberately.** Measured on the first run of this oracle,
    /// over the 90 confluence classes the family builds (15 designs x 3 axes x 2
    /// directions): **38 do not converge**, and **0 outputs fail the fixed-point
    /// half**. So the failure is entirely disagreement between spellings, never
    /// drift within one.
    ///
    /// The 38 fall out by axis, and the split is the finding:
    ///
    /// | axis | non-confluent classes of 30 |
    /// |---|---:|
    /// | `g.` | **0** |
    /// | `c.` (single exon) | 14 |
    /// | `cx.` (three exons, 5'UTR) | 24 |
    ///
    /// A reverse-complement block converges on the genomic axis and stops
    /// converging on a coding reference, and every failing class crosses a CDS
    /// boundary. Smallest reproducer — the separation-0 design at offset 0, whose
    /// block is four bases and whose three spellings are three *distinct fixed
    /// points*:
    ///
    /// ```text
    /// NM_TESTX.1:c.-3_1inv                       ->  NM_TESTX.1:c.-3_1inv
    /// NM_TESTX.1:c.-3_1delinsTGGT                ->  NM_TESTX.1:c.-3_1delinsTGGT
    /// NM_TESTX.1:c.[-3_-2delinsTG;-1_1delinsGT]  ->  unchanged
    /// ```
    ///
    /// The same three spellings on the single-exon coding reference all converge to
    /// `NM_TEST.1:c.1_4inv`, so the block is not the problem — the 5'UTR crossing
    /// is. At separation 1 the same thing happens one level down:
    /// `NM_TEST.1:c.[9_10delinsTT;12_13delinsAA;15_*1delinsTT]` re-types the first
    /// two members `inv` and leaves the third, the one straddling `CDS_END`, as a
    /// `delins`.
    ///
    /// Part of this is #1517 (an `inv` spelling and a split spelling that are both
    /// fixed points), but the axis dependence is not: #1517 is axis-neutral by its
    /// own scope note, and `g.` is clean here. That half is unfiled.
    ///
    /// **Not weakened to pass.** Both halves are asserted together on purpose; a
    /// version that only checked the fixed-point half would be green and would be
    /// measuring nothing this family was added for.
    #[test]
    #[ignore = "red on first run: 38 of 90 confluence classes disagree, all on coding \
                axes crossing a CDS boundary; see this test's docs"]
    fn the_three_spellings_of_a_design_converge_to_one_fixed_point() {
        let mut split: Vec<String> = Vec::new();
        let mut drifting: Vec<String> = Vec::new();
        for design in revcomp_designs() {
            for (axis, direction, dir_label) in axes_and_directions() {
                let inputs = revcomp_inputs_for(axis, &design);
                let outputs: Vec<String> = inputs
                    .iter()
                    .map(|input| normalize_one(axis, &design.core, input, direction))
                    .collect();
                let distinct: std::collections::BTreeSet<&str> =
                    outputs.iter().map(String::as_str).collect();
                if distinct.len() > 1 {
                    split.push(format!(
                        "[{} {axis} {dir_label}] {} distinct outputs:\n{}",
                        design.label,
                        distinct.len(),
                        inputs
                            .iter()
                            .zip(&outputs)
                            .map(|(i, o)| format!("      {i}\n        -> {o}"))
                            .collect::<Vec<_>>()
                            .join("\n")
                    ));
                }
                for output in &outputs {
                    let again = normalize_one(axis, &design.core, output, direction);
                    if &again != output {
                        drifting.push(format!(
                            "[{} {axis} {dir_label}] {output}\n        -> {again}",
                            design.label
                        ));
                    }
                }
            }
        }
        assert!(
            split.is_empty() && drifting.is_empty(),
            "{} designs do not converge and {} outputs are not fixed points\n\
             --- non-confluent ---\n{}\n--- not a fixed point ---\n{}",
            split.len(),
            drifting.len(),
            split.join("\n"),
            drifting.join("\n")
        );
    }

    // -- Oracle 3: separation-rule conformance ------------------------------

    /// One member of a normalized output, in its accession's own 0-based frame.
    ///
    /// **An insertion consumes no reference position.** Modelling every member as
    /// the closed interval `[lo, hi]` and leaving an insertion's interval *empty*
    /// (`lo = A + 1`, `hi = A`) is what makes the single formula
    /// `sep = next.lo - prev.hi - 1` hold for every pair. Treating `A_B` as a
    /// two-base consumed span instead double-counts the junction and shifts every
    /// separation in the distribution by one — which invalidates the whole census,
    /// not just the insertion rows.
    struct MemberFacts {
        lo: i64,
        hi: i64,
        reference: String,
        alternate: String,
    }

    /// Read every member of `variant` as a [`MemberFacts`].
    ///
    /// Via SPDI rather than by decoding each `NaEdit` by hand: SPDI is already the
    /// `(position, deleted, inserted)` triple this needs, it resolves the short
    /// forms that name no bases (`del`, `dup`, `inv`) against the reference, and it
    /// puts every axis in one frame — so a `c.` allele and a `g.` allele are
    /// measured by the same arithmetic instead of by two decoders that can drift.
    fn member_facts_of(
        variant: &ferro_hgvs::HgvsVariant,
        provider: &MockProvider,
    ) -> Result<Vec<MemberFacts>, String> {
        use ferro_hgvs::HgvsVariant;
        let members: Vec<&HgvsVariant> = match variant {
            HgvsVariant::Allele(allele) => allele.variants.iter().collect(),
            other => vec![other],
        };
        let mut facts = Vec::with_capacity(members.len());
        for member in members {
            let spdi = ferro_hgvs::spdi::hgvs_to_spdi(member, provider)
                .map_err(|e| format!("{member}: {e}"))?;
            let lo = spdi.position as i64;
            facts.push(MemberFacts {
                lo,
                hi: lo + spdi.deletion.len() as i64 - 1,
                reference: spdi.deletion,
                alternate: spdi.insertion,
            });
        }
        facts.sort_by_key(|f| (f.lo, f.hi));
        Ok(facts)
    }

    /// Which reference columns of `reference` are matched in **every** minimal
    /// alignment of `reference` against `alternate`.
    ///
    /// A `delins` of unequal length has many minimal alignments, so "unchanged
    /// interior" is alignment-dependent and a naive column-by-column comparison
    /// invents violations that the description does not commit to. The rule used
    /// here is the conservative one: a column counts as unchanged only when **no**
    /// minimal-cost path consumes it — neither by deleting it nor by substituting
    /// it. Forward and backward Levenshtein DP give exactly that, since a path
    /// through column `i` costs `prefix(i, j) + step + suffix(i + 1, j')`, and the
    /// column is free to be changed iff some `(j, step)` reaches the optimum.
    ///
    /// An insertion is deliberately not a way to "consume" a reference column: it
    /// advances only the alternate, so it leaves the column matched.
    ///
    /// The direction of the conservatism matters. This can **under**-report — a
    /// column changed in every alignment a human would draw may still be reachable
    /// by some exotic minimal path — and it must never **over**-report, because a
    /// violation claimed here is a claim that the description is wrong.
    fn forced_unchanged(reference: &[u8], alternate: &[u8]) -> Vec<bool> {
        let (n, m) = (reference.len(), alternate.len());
        let idx = |i: usize, j: usize| i * (m + 1) + j;

        let mut prefix = vec![0usize; (n + 1) * (m + 1)];
        for i in 0..=n {
            prefix[idx(i, 0)] = i;
        }
        for j in 0..=m {
            prefix[idx(0, j)] = j;
        }
        for i in 1..=n {
            for j in 1..=m {
                let cost = usize::from(reference[i - 1] != alternate[j - 1]);
                prefix[idx(i, j)] = (prefix[idx(i - 1, j)] + 1)
                    .min(prefix[idx(i, j - 1)] + 1)
                    .min(prefix[idx(i - 1, j - 1)] + cost);
            }
        }

        let mut suffix = vec![0usize; (n + 1) * (m + 1)];
        for i in 0..=n {
            suffix[idx(i, m)] = n - i;
        }
        for j in 0..=m {
            suffix[idx(n, j)] = m - j;
        }
        for i in (0..n).rev() {
            for j in (0..m).rev() {
                let cost = usize::from(reference[i] != alternate[j]);
                suffix[idx(i, j)] = (suffix[idx(i + 1, j)] + 1)
                    .min(suffix[idx(i, j + 1)] + 1)
                    .min(suffix[idx(i + 1, j + 1)] + cost);
            }
        }

        let distance = prefix[idx(n, m)];
        (0..n)
            .map(|i| {
                let changeable = (0..=m).any(|j| {
                    let before = prefix[idx(i, j)];
                    // Delete reference[i].
                    if before + 1 + suffix[idx(i + 1, j)] == distance {
                        return true;
                    }
                    // Substitute reference[i] for alternate[j].
                    j < m
                        && reference[i] != alternate[j]
                        && before + 1 + suffix[idx(i + 1, j + 1)] == distance
                });
                !changeable
            })
            .collect()
    }

    /// Maximal runs of `true` with a `false` on both sides, as `(offset, length)`.
    ///
    /// Interior only: a leading or trailing run of unchanged columns is a
    /// description that failed to trim its own flanks, which is a different defect
    /// with a different rule, and folding the two together would make each
    /// unreadable.
    fn interior_runs(flags: &[bool]) -> Vec<(usize, usize)> {
        let mut runs = Vec::new();
        let mut i = 0;
        while i < flags.len() {
            if !flags[i] {
                i += 1;
                continue;
            }
            let start = i;
            while i < flags.len() && flags[i] {
                i += 1;
            }
            if start > 0 && i < flags.len() {
                runs.push((start, i - start));
            }
        }
        runs
    }

    /// The codon a reference column falls in, or `None` when it has no frame.
    ///
    /// `general.md:35`'s carve-out is about a single amino acid, so it needs a real
    /// reading frame: never on `g.`, and never in a UTR of a coding reference. A
    /// column with no codon cannot license a merge across an unchanged base, which
    /// is the half of the rule that is easy to lose.
    fn codon_index(axis: &str, accession: &str, column: i64) -> Option<i64> {
        if accession != TX_ACCESSION && accession != TX_MULTI_ACCESSION {
            return None;
        }
        let (cds_start, cds_end) = cds_bounds(axis);
        let transcript_position = column + 1;
        if transcript_position < cds_start as i64 || transcript_position > cds_end as i64 {
            return None;
        }
        Some((transcript_position - cds_start as i64) / 3)
    }

    /// A way one description can merge more than the spec allows.
    ///
    /// All three are **over**-merges, which is the direction the rules constrain:
    /// `general.md:34` compels one description for adjacent changes and separate
    /// ones past a single unchanged base, and `general.md:35` is the one carve-out.
    /// Splitting further than required is not a violation of either.
    #[derive(Debug)]
    enum SeparationFinding {
        /// Two members with no unchanged reference column between them
        /// (`general.md:34`, `DNA/delins.md:16-18`).
        AdjacentMembersLeftSplit { at: i64 },
        /// One member spans two or more forced-unchanged columns
        /// (`general.md:34`).
        MergedAcrossUnchangedInterior {
            at: i64,
            length: usize,
            inversion: bool,
        },
        /// One member spans a single forced-unchanged column with no codon to
        /// license it (`general.md:35`).
        MergedAcrossOneColumnWithoutACodon { at: i64, inversion: bool },
    }

    impl SeparationFinding {
        /// How the offending member reads, which is what decides whether the
        /// finding is settled or contested.
        ///
        /// A member whose payload is the exact reverse complement of its own span
        /// is an **inversion**, and `general.md:56` ranks inversion above the
        /// residual `delins` — which is the argument #1517 makes for preferring
        /// `inv` over the split. `general.md:34` says "and **not** as a delins",
        /// so it does not literally reach an `inv`. A member that is *not* an
        /// inversion has no such defence: it is a delins, and the clause names it.
        fn spelling(inversion: bool) -> &'static str {
            if inversion {
                "member, itself an exact inversion (the contested #1517 class),"
            } else {
                "`delins` member"
            }
        }

        /// Whether this finding's offending member is an exact inversion.
        fn is_inversion(&self) -> bool {
            match self {
                Self::AdjacentMembersLeftSplit { .. } => false,
                Self::MergedAcrossUnchangedInterior { inversion, .. }
                | Self::MergedAcrossOneColumnWithoutACodon { inversion, .. } => *inversion,
            }
        }

        /// The finding as a line a reader can act on: what was violated, where,
        /// and which clause says so.
        ///
        /// Spelled out rather than left to `{:?}` because the column is the
        /// reproducer — a finding without its 0-based reference column cannot be
        /// checked by hand against the block it came from.
        fn describe(&self) -> String {
            match self {
                Self::AdjacentMembersLeftSplit { at } => format!(
                    "general.md:34 — members meet at reference column {at} with no unchanged \
                     base between them, so they are one description"
                ),
                Self::MergedAcrossUnchangedInterior {
                    at,
                    length,
                    inversion,
                } => format!(
                    "general.md:34 — one {} spans {length} unchanged reference columns from \
                     {at}, which no minimal alignment can consume",
                    Self::spelling(*inversion)
                ),
                Self::MergedAcrossOneColumnWithoutACodon { at, inversion } => format!(
                    "general.md:35 — one {} spans the unchanged reference column {at} with \
                     no codon shared by the changes flanking it",
                    Self::spelling(*inversion)
                ),
            }
        }
    }

    /// Judge one normalized output against the separation rules.
    ///
    /// `Err` means the output could not be read as a set of reference footprints at
    /// all — a sentinel, an unparseable string, an edit SPDI declines — and is
    /// counted as unevaluated rather than silently passed.
    fn separation_findings(
        axis: &str,
        core: &str,
        output: &str,
    ) -> Result<Vec<SeparationFinding>, String> {
        let variant = parse_hgvs(output).map_err(|e| format!("{output}: {e}"))?;
        let provider = provider_for(axis, core);
        let named = match &variant {
            ferro_hgvs::HgvsVariant::Allele(allele) => allele.variants.first().unwrap_or(&variant),
            other => other,
        };
        let accession = named
            .accession()
            .map(ToString::to_string)
            .unwrap_or_default();
        let facts = member_facts_of(&variant, &provider)?;
        let mut findings = Vec::new();
        for pair in facts.windows(2) {
            if pair[1].lo - pair[0].hi - 1 == 0 {
                findings.push(SeparationFinding::AdjacentMembersLeftSplit { at: pair[1].lo });
            }
        }
        for member in &facts {
            let inversion =
                !member.reference.is_empty() && member.alternate == revcomp(&member.reference);
            let flags = forced_unchanged(member.reference.as_bytes(), member.alternate.as_bytes());
            for (offset, length) in interior_runs(&flags) {
                let at = member.lo + offset as i64;
                if length >= 2 {
                    findings.push(SeparationFinding::MergedAcrossUnchangedInterior {
                        at,
                        length,
                        inversion,
                    });
                    continue;
                }
                // A single unchanged column is merged across only when the changes
                // on either side of it sit in one codon. Both flanks are read, not
                // just one: a merge is licensed by the amino acid the pair shares,
                // so a pair straddling a codon boundary is not licensed even though
                // each half has a frame.
                let before = codon_index(axis, &accession, at - 1);
                let after = codon_index(axis, &accession, at + 1);
                if before.is_none() || before != after {
                    findings.push(SeparationFinding::MergedAcrossOneColumnWithoutACodon {
                        at,
                        inversion,
                    });
                }
            }
        }
        Ok(findings)
    }

    /// The forced-unchanged DP behaves, on cases whose answers are computable by
    /// hand — including the two that make it non-trivial.
    #[test]
    fn the_forced_unchanged_dp_only_claims_columns_no_minimal_alignment_moves() {
        let flags = |r: &str, a: &str| forced_unchanged(r.as_bytes(), a.as_bytes());

        // The worked case: equal lengths, so the only minimal alignment is the
        // column-by-column one and the four interior columns are all forced.
        assert_eq!(
            flags("AATGCACA", "TGTGCATT"),
            vec![false, false, true, true, true, true, false, false]
        );
        assert_eq!(interior_runs(&flags("AATGCACA", "TGTGCATT")), vec![(2, 4)]);

        // Unequal lengths, where a naive column comparison invents a violation.
        // `AAA` -> `AA` deletes any one of the three A's, so no column is forced.
        assert!(flags("AAA", "AA").iter().all(|forced| !forced));

        // …and where the shift genuinely cannot reach the middle: `CAG` -> `TAG`
        // has one minimal alignment (substitute the C), so `A` and `G` are forced,
        // but neither is *interior* — there is no changed column after them.
        assert_eq!(flags("CAG", "TAG"), vec![false, true, true]);
        assert!(interior_runs(&flags("CAG", "TAG")).is_empty());

        // An insertion does not consume a reference column, so the reference is
        // wholly matched and no interior run exists.
        assert!(flags("AC", "AGC").iter().all(|forced| *forced));

        // A genuine spanning delins over an unchanged single base: `CAG` -> `TAT`
        // has cost 2, and the middle `A` cannot be consumed by any 2-cost path.
        assert_eq!(flags("CAG", "TAT"), vec![false, true, false]);
        assert_eq!(interior_runs(&flags("CAG", "TAT")), vec![(1, 1)]);
    }

    /// The member model puts an insertion on an empty interval, so one separation
    /// formula covers every pair.
    #[test]
    fn an_insertion_consumes_no_reference_column() {
        let core = &corpus_sequences(1)[0];
        let provider = genomic_provider(core);
        let variant = parse_hgvs("NC_TEST.1:g.[257_258del;258_259insTT]").expect("parses");
        let facts = member_facts_of(&variant, &provider).expect("both members convert");
        assert_eq!(facts.len(), 2);
        // The deletion consumes 0-based 256..=257; the insertion sits at the
        // junction after 257 and consumes nothing, so its interval is empty.
        assert_eq!((facts[0].lo, facts[0].hi), (256, 257));
        assert_eq!((facts[1].lo, facts[1].hi), (258, 257));
        // …and the uniform formula reads that junction as adjacent, not as a
        // two-base span that would put the two members one column apart.
        assert_eq!(facts[1].lo - facts[0].hi - 1, 0);
    }

    /// **Oracle 3.** No normalized output merges across more than the spec allows.
    ///
    /// Checked against the sequence, not against a pinned string: because
    /// `(reference, alternate)` is known by construction, every member's own
    /// `(deleted, inserted)` pair is recoverable and the rule can be evaluated
    /// directly. `general.md:34` — changes separated by one or more unchanged
    /// nucleotides are described individually, and adjacent ones as a single
    /// description. `general.md:35` — the one carve-out, a single unchanged
    /// nucleotide inside one codon, which needs a reading frame and so cannot apply
    /// on `g.` or in a UTR. `DNA/delins.md:16-18` for the merged form itself.
    /// **Currently red, deliberately.** Measured on the first run of this oracle,
    /// over the 270 rows the family contributes: **138 of 270 outputs violate**,
    /// 132 are clean, 0 unevaluated. The 176 findings partition into two classes
    /// that are not the same claim, and the split is very lopsided:
    ///
    /// | class | findings | status |
    /// |---|---:|---|
    /// | the offending member is itself an **exact reverse complement** | 174 | **contested — #1517** |
    /// | two adjacent members, neither an inversion, left split | 2 | settled; `general.md:34` |
    ///
    /// The contested class is #1517's own subject, read in the opposite direction.
    /// `general.md:34` says two separated variants are described individually "and
    /// **not** as a `delins`" — it does not name `inv`, and `general.md:56` ranks
    /// inversion above the residual `delins`, which is exactly #1517's argument for
    /// preferring `NM_004006.2:c.76_83inv` over the split. #1517 records the
    /// counter-reading this oracle encodes as an available reading, not a settled
    /// one. So these 174 are a decision the repository has not taken rather than a
    /// defect this oracle found — and the classification is by the member's
    /// *content*, not its spelling, so a whole-block reverse complement written
    /// `delins` lands here too.
    ///
    /// The settled class is two rows, and it is new. The separation-0 design at
    /// offset 0 on `cx`, in both directions:
    ///
    /// ```text
    /// core  ACCAGCTAGCTAGCTAGCTA   (CDS_START_MULTI = 4, so c.-3 == transcript 1)
    /// NM_TESTX.1:c.[-3_-2delinsTG;-1_1delinsGT]  ->  unchanged
    /// NM_TEST.1:c.[1_2delinsTG;3_4delinsGT]      ->  NM_TEST.1:c.1_4inv
    /// ```
    ///
    /// Neither member is an inversion on its own (`AC` -> `TG`, `CA` -> `GT`); only
    /// their union is, so `general.md:56` offers no defence and `general.md:34`'s
    /// plainest case applies — two adjacent changes are one description. The same
    /// input merges on the single-exon coding reference and does not on the
    /// multi-exon one, purely because the block starts in the 5'UTR.
    ///
    /// That axis dependence runs through the contested class too and is worth
    /// separating from #1517, which its own scope note calls axis-neutral:
    /// `NM_TEST.1:c.9_*5delinsTTCGTATACGTT` holds eight unchanged interior columns
    /// inside one member, while the identical block on `g.` comes out as
    /// `g.[265_266inv;275_276inv]`.
    ///
    /// **Not weakened to pass.** Narrowing this to the settled class alone would
    /// still be red, and narrowing it to zero would make it documentation.
    #[test]
    #[ignore = "red on first run: 2 settled violations (a CDS-boundary split that the \
                single-exon axis merges, unfiled) and 174 contested ones (#1517); see \
                this test's docs"]
    fn no_normalized_output_merges_across_an_unchanged_interior() {
        let (mut clean, mut unevaluated) = (0usize, 0usize);
        let (mut contested, mut settled) = (0usize, 0usize);
        let mut violations: Vec<String> = Vec::new();
        for_each_revcomp_row(|design, axis, direction, input, output| {
            match separation_findings(axis, &design.core, output) {
                Err(_) => unevaluated += 1,
                Ok(findings) if findings.is_empty() => clean += 1,
                Ok(findings) => {
                    for finding in &findings {
                        if finding.is_inversion() {
                            contested += 1;
                        } else {
                            settled += 1;
                        }
                    }
                    violations.push(format!(
                        "[{} {} {}] {input}\n      -> {output}\n{}",
                        design.label,
                        axis,
                        direction,
                        findings
                            .iter()
                            .map(|f| format!("      {}", f.describe()))
                            .collect::<Vec<_>>()
                            .join("\n")
                    ));
                }
            }
        });
        assert!(
            violations.is_empty(),
            "{} of {} normalized outputs violate the separation rules ({clean} clean, \
             {unevaluated} unevaluated); {settled} findings name a `delins` member and are \
             settled, {contested} name an exact inversion and are #1517's open \
             question:\n{}",
            violations.len(),
            violations.len() + clean + unevaluated,
            violations.join("\n")
        );
        assert!(
            clean > 0,
            "no row could be evaluated at all, so this oracle measured nothing"
        );
    }

    fn fixture(name: &str, contents: &str) -> PathBuf {
        let dir = std::env::temp_dir().join("ferro-dump-corpus-test");
        let _ = fs::create_dir_all(&dir);
        let path = dir.join(name);
        fs::write(&path, contents).expect("write fixture");
        path
    }

    /// A dump round-trips through the reader, and the reader rejects a duplicate
    /// key rather than silently keeping one of the two rows.
    #[test]
    fn the_reader_rejects_a_duplicate_key() {
        let dir = std::env::temp_dir().join("ferro-dump-corpus-test");
        let _ = fs::create_dir_all(&dir);
        let path = dir.join("dup.tsv");
        fs::write(
            &path,
            format!("{HEADER}AC\tg\t3prime\tf\tX\tY\tfalse\nAC\tg\t3prime\tf\tX\tZ\tfalse\n"),
        )
        .expect("write fixture");
        let err = read_dump(&path).expect_err("a duplicate key must be rejected");
        assert!(err.contains("duplicate key"), "unexpected error: {err}");
    }

    /// Comparing a dump with itself reports no movement — the degenerate case that
    /// would otherwise make every table look alarming.
    #[test]
    fn comparing_a_dump_with_itself_reports_no_movement() {
        let dir = std::env::temp_dir().join("ferro-dump-corpus-test");
        let _ = fs::create_dir_all(&dir);
        let path = dir.join("self.tsv");
        let rows = dump(2);
        let mut out = String::from(HEADER);
        for row in &rows {
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                row.reference,
                row.axis,
                row.direction,
                row.family,
                row.input,
                row.output,
                row.was_fixed_point
            );
        }
        fs::write(&path, &out).expect("write dump");
        let report = compare(&path, &path).expect("self-comparison must succeed");
        assert!(
            report.contains("**moved** | **0 (0.0%)**"),
            "a dump compared with itself must report zero movement:\n{report}"
        );
    }
}
