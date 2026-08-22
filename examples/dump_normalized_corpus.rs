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
//!    same path on both sides no matter how many families exist. #1403 capped a
//!    partition at `MAX_SPLIT_BLOCK` — **1024 at the time**, 4096 since #1899
//!    derived it from `MAX_CANONICAL_WINDOW` — and measured a guaranteed zero here
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
//!    **Several corpus totals appear in this file and every one of them is
//!    correct — do not "reconcile" them.** Each is the corpus as it stood when
//!    some measurement was taken, so each is anchored to what changed it:
//!
//!    - `78,028` — before `separated_revcomp_runs` existed, and the only
//!      denominator a measurement taken before that family was added can honestly
//!      carry;
//!    - `78,298` — after it, and `78,298 - 78,028 = 270` is exactly that family's
//!      row count. The two differ by a digit transposition as well as by 270,
//!      which is a coincidence and has already been read once as a typo in one of
//!      them;
//!    - `85,642` — after #1606 added the protein axis. It is the denominator of
//!      item 6 below, and stays attached to that measurement;
//!    - `86,398` — after #1752 added `repeat_beside_a_sibling` (756 rows, in two
//!      steps: 378 when the family landed, 378 more when it was made to emit
//!      shrinking repeats as well as growing ones, which is why `86,020` appears
//!      in that PR's own disclosure as a superseded figure).
//!
//!    When you quote a figure, quote the denominator it was measured against.
//!
//!    **This list is a history, not a census, and deliberately does not say which
//!    entry is current** — the current total is whatever the generator prints
//!    today, and asking the file is the wrong way to find out:
//!
//!    ```text
//!    cargo run --release --features dev --example dump_normalized_corpus -- --out /tmp/now.tsv
//!    ```
//!
//!    That framing is the #1947 correction, and it is the second one this
//!    paragraph has needed. It first said there were exactly **two** totals while
//!    item 6 eight lines below already quoted the third — the more misleading
//!    half of that error, since a reader who trusts an explicit "there are
//!    exactly two, here they are" has no reason to look for a third. Corrected to
//!    three, it was stale again **the next day**, when #1752 added a family: the
//!    entry marked "the corpus **now**" had stopped being now, and nothing about
//!    it said so. Counting the entries and naming a current one are both claims
//!    that go stale on a change that has no reason to touch this comment, which
//!    is the same defect as the family counts below — except that a corpus total
//!    cannot be derived at compile time, so the fix here is to stop making the
//!    claim rather than to import it. The *method* this paragraph states was
//!    right throughout; its own bookkeeping is what kept drifting.
//!
//! 6. **…and only on the molecule types it builds** (#1606). Every family above is
//!    spelled from nucleotides, so until the protein axis was added this corpus
//!    could not emit one `p.` row and a protein-scoped change measured a structural
//!    `0 of 78,298`. #1606 is the live instance: it moves an equal-length protein
//!    `delins` onto its members and had to *report* the zero and name the generator
//!    as its cause, because there was no measurement to quote. With the axis in
//!    place the same change measures **1,584 of 85,642 rows, all in
//!    `protein_equal_length_delins` and none anywhere else**.
//!
//!    Note what these instances have in common and why the list keeps growing:
//!    fixing one blindness does not reveal the next. Geometry (#1456) was only
//!    findable once the conflict families existed, scale (#1460) once geometry was
//!    fixed, transcript geometry (#1478) once scale was — and molecule type was
//!    invisible behind all three, because a corpus with three axes looks well
//!    covered. `compare` names the families it covered (#1459), which catches a
//!    missing *family*; it cannot tell you the axis you never built.
//!
//!    **That chain is illustrative, not an inventory, and this comment states no
//!    count of it on purpose.** It said "the four instances" until #1947, by which
//!    time there were at least six — sequence structure (#1517) and reversed
//!    ranges (#1917) are both missing from the chain above, and #1752 added a
//!    seventh from inside this very file, its `repeat_beside_a_sibling` family
//!    having first shipped able to build only *growing* repeats. A count of
//!    blindnesses is itself a number that goes stale every time the thing it
//!    counts happens, which is the failure this whole section is about. The
//!    maintained list is `CLAUDE.md`'s, under "Assert the property. Measure the
//!    count. Never let a count BE the property".
//!
//! ## One family knows its own ground truth, and that buys three oracles
//!
//! Everything above measures **movement** — two dumps, one diff — and never
//! whether an output is *right*. Eighteen of the twenty-two families cannot be asked:
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

use ferro_hgvs::normalize::{NormalizeConfig, ShuffleDirection, MAX_CANONICAL_BLOCK};
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

/// The protein reference a `p.` row is drawn against (#1606).
///
/// A bare `NP_` accession with no nucleotide counterpart, which is deliberate: it
/// takes `MockProvider::add_protein` and nothing else, so the protein axis needs
/// neither a transcript nor a CDS and cannot be perturbed by the exon geometry the
/// `c.`/`cx.` references carry.
const PROTEIN_ACCESSION: &str = "NP_TEST.1";

/// Residues protein cores are built from, as `(one-letter, three-letter)`.
///
/// `Met` leads the table because it is residue 1 of every core and is **never
/// drawn** into the body — cores draw from `[1..]`. Two reasons, and the second is
/// the one that would silently weaken the corpus. A protein starts at `Met`, so a
/// core that did not would not be a protein; and an edit reaching residue 1 emits
/// the illegal start-loss spelling filed as #1607, so a `Met` in the body would let
/// a shift *arrive* at residue 1 and turn a representation measurement into a
/// rediscovery of that bug. Keeping `Met` out of the body means the only residue-1
/// interaction is one the families never construct.
const RESIDUE_CODES: &[(char, &str)] = &[
    ('M', "Met"),
    ('A', "Ala"),
    ('G', "Gly"),
    ('S', "Ser"),
    ('L', "Leu"),
    ('K', "Lys"),
    ('V', "Val"),
];

/// Residues in a protein core, `Met` included.
///
/// Twenty, matching the DNA cores, so the two halves of the corpus report row
/// counts on the same scale and a per-core figure means the same thing on either.
const PROTEIN_LEN: usize = 20;

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
/// (#1624). Those two are the same string for eighteen of the twenty-two families,
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
    /// eighteen of the twenty-two families this **is** the reference sequence; for
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
    // **The fifth instance of the corpus-blindness class** (#1749), after member
    // geometry (#1456), scale (#1460), transcript geometry (#1478) and molecule
    // type (#1606). `repeat_expansion` above feeds a repeat, but always a *lone*
    // one — so no row in this corpus has ever placed a `repeat` beside a sibling
    // it could collide with, and every question about a repeat's write footprint
    // was structurally unaskable here.
    //
    // That is exactly the geometry #1749 changes: a `repeat` may rewrite the
    // tract it spans *and* lands any expansion at the junction 3' of it, and the
    // three passes in `normalize::overlap` previously gave three different
    // answers for it. Reporting `0 moved` off a corpus that cannot build the
    // shape would have been a claim about the corpus.
    (
        "repeat_beside_a_sibling",
        "#1749 — a ranged repeat, GROWING and SHRINKING, paired with a sibling that intersects \
         its tract or its junction",
    ),
    // Another instance of the same corpus-blindness class, and the one #1610 hit.
    // No ordinal is claimed for it: the header above and the `repeat_beside_a_sibling`
    // comment already number this class differently, and a third restatement is how
    // that drift got started.
    //
    // Every lone `delins` this file built before it is either **equal**-length
    // (`delins_hiding_an_inversion`, `separated_revcomp_runs`) or a **net
    // insertion** (`split_vs_spanning_delins`, whose spanning spelling is a
    // six-base payload over a five-base span). #1610 is scoped to net
    // *deletions*, so it measured a guaranteed structural zero: not one row of
    // the corpus could reach the rule, whatever it did.
    //
    // Note it is the *direction* that was missing, not the shape — the corpus
    // has built lone spanning `delins` rows since #1401. A family list read for
    // shapes looks complete here, which is why the zero would have been quoted.
    (
        "lone_net_deletion_delins",
        "#1610 — a lone minimal net-deletion delins whose payload partly coincides with the reference",
    ),
    (
        "authored_repeat_beside_a_sibling",
        "#1946 — an authored repeat that survives as one member of a two-member allele, \
         which is the only shape that reaches `lower_repeat_edits`",
    ),
    // Another instance of the corpus-blindness class, and the one #2192 hit. Every
    // family above pairs members so that no single member is a **multi-column
    // contiguous run** that the solid-run collapse would coalesce: the multi-member
    // ones are lone substitutions (`three_member_allele`), disjoint spans
    // (`nested_spans`, `partial_overlap_spans`) or a span beside a point change
    // (`del_plus_sub`, `dup_plus_sub`). So no row in this corpus has ever placed a
    // coalescible run beside a **distant** cis member, which is exactly the geometry
    // #2192 is about: a second member stretches the whole-variant hull across a gap,
    // the changed-column set goes non-contiguous, and `coalesce_solid_run`'s
    // contiguity gate declines — the #2174 collapse reaching only single-run
    // alleles. Measuring `0 moved` off a corpus that cannot build run-count >= 2
    // would have been a claim about the corpus.
    (
        "run_beside_a_distant_member",
        "#2192 — a contiguous coalescible run (adjacent changes) beside a distant cis member \
         the whole-hull collapse cannot reach; the run coalesces per-run once a sibling \
         makes the hull non-contiguous",
    ),
];

/// The shape families drawn against the **protein** cores (#1606).
///
/// Kept out of `FAMILIES` because they are not crossed with the DNA axes: a `p.`
/// description is spelled from residues, so nothing above can be evaluated on a
/// protein core and nothing here on a nucleotide one. The split is the same one
/// `LONG_FAMILY` and `REVCOMP_FAMILY` make, and it has the same consequence — every
/// existing row is byte-identical, because an added axis, like an added family,
/// only appends.
///
/// **This is the fourth instance of the corpus-blindness class**, after member
/// geometry (#1456), scale (#1460) and transcript geometry (#1478). The corpus
/// built no protein row at all, so a change scoped to the protein axis measured a
/// structural `0 of 78,298` — indistinguishable in the output from a change that
/// genuinely moved nothing. #1606 is the live instance: it moves an equal-length
/// `p.` delins onto its members, and had to report the zero and name the generator
/// as the reason, because there was no measurement to quote.
const PROTEIN_FAMILIES: &[(&str, &str)] = &[
    (
        "protein_shift_del",
        "#91 — a deletion inside a residue run, which shifts to the C-terminal end",
    ),
    (
        "protein_shift_dup",
        "#91 / #92 — a duplication inside a residue run",
    ),
    (
        "protein_ins_becomes_dup",
        "#92 — an insertion whose payload copies its 5' neighbour, so it re-types as a dup",
    ),
    (
        "protein_equal_length_delins",
        "#1606 — an equal-length delins whose interior residues are unchanged",
    ),
    (
        "protein_cis_separated",
        "#1232 / #1401 on the protein axis — two members separated by unchanged residues",
    ),
];

/// The family drawn against the **long** cores, and the one shape in this file
/// whose point is its size rather than its geometry (#1460).
///
/// Kept out of `FAMILIES` on purpose. The eighteen families there are crossed with
/// every short core, and crossing a kilobase core with all of them would multiply
/// the dump cost while adding no coverage: `MAX_SPLIT_BLOCK` gates on the *length*
/// of the block being partitioned, not on how its members are arranged.
const LONG_FAMILY: (&str, &str) = (
    "long_block_inversion",
    "#1403 — a near-palindromic block straddling the canonical block limit",
);

/// The long cores, as `(label, sequence)`: one block the canonicalizer accepts
/// and one it refuses, which is the pair a change to the limit moves between.
///
/// # Derived from the constant, never restated (#1925)
///
/// These used to be the literals `[1024, 1100]`, straddling `MAX_SPLIT_BLOCK`
/// when that was `1024`. #1899 derived that cap from `MAX_CANONICAL_WINDOW`,
/// taking it to **4096**, and both cores landed on the same side of it — so the
/// one family whose point is *size* stopped measuring size, and the guard below
/// went on passing because it restated `1024` as a literal of its own.
///
/// Two things changed as a result, and the second matters more than the first.
/// The lengths are now computed from [`MAX_CANONICAL_BLOCK`], so the pair cannot
/// be left behind again. And they straddle the limit that **actually fires on
/// the shipping path** — the padded-window test — rather than `MAX_SPLIT_BLOCK`,
/// whose only functional reader also demands `reference.len() != result.len()`
/// and so was never reachable by this family's equal-length `inv` at any size.
///
/// # The gate measures the CHANGED INTERVAL, not the core length
///
/// This is the trap, and the first attempt at this fix fell into it: cores of
/// exactly `MAX_CANONICAL_BLOCK` and `+2` produced **byte-identical behaviour**,
/// because the interval the window is built around runs from the first differing
/// base to the last, and [`near_palindromic_core`] puts those
/// [`EDGE_PERTURBATION`] bases in from each end. So a core of `n` presents a
/// changed interval of `n - 2 * EDGE_PERTURBATION`, and cores straddling the
/// limit give intervals that are both comfortably under it.
///
/// Measured before and after: at cores `[3840, 3842]` both rows normalized
/// identically (`g.257_4096inv` and `g.257_4098inv` each collapsing to four
/// substitutions); the corpus straddled nothing while looking like it did. The
/// cores below are offset so the **intervals** land at the limit and one base
/// past it.
///
/// The margin is `2`, not `1`: [`near_palindromic_core`] refuses an odd length,
/// and the limit itself is even.
///
/// The **label** is what lands in the dump's `reference` column, not the
/// sequence: a kilobase core would otherwise be repeated verbatim on every row.
/// Labels are as good a row identity as the sequence — they only have to be
/// unique and stable.
fn long_corpus_sequences() -> Vec<(String, String)> {
    let at_limit = MAX_CANONICAL_BLOCK as usize + 2 * EDGE_PERTURBATION;
    [at_limit, at_limit + 2]
        .into_iter()
        .map(|len| (format!("nearpalindrome_{len}"), near_palindromic_core(len)))
        .collect()
}

/// How far in from each end [`near_palindromic_core`] breaks the palindrome.
///
/// Load-bearing for [`long_corpus_sequences`], not a free parameter: it sets the
/// gap between a core's length and the width of the changed interval the
/// normalizer's window is actually built around (`len - 2 * EDGE_PERTURBATION`),
/// which is what decides whether the long cores straddle the gate or merely look
/// as though they do.
const EDGE_PERTURBATION: usize = 10;

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
    for offset in [EDGE_PERTURBATION, half + 7] {
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

/// The family drawn against the **tandem-tract designs**, and the one shape in this
/// file whose point is that a *presentational* spelling and a sibling are in the
/// same allele (#1946).
///
/// Kept out of [`FAMILIES`] for the same reason [`REVCOMP_FAMILY`] is: it cannot be
/// drawn on an arbitrary core. A tandem tract of a chosen unit, at a chosen offset,
/// bounded on both sides so it cannot extend, is a property of the *sequence* —
/// `repeat_expansion`, the only other family that describes a tract, reads whatever
/// tract the xorshift happened to produce (`first_tandem_triplet`, 21 of 48 cores)
/// and emits a **lone** member over it at one offset. So the corpus could ask what
/// normalization does to a repeat with no sibling and never what it does to one with.
///
/// # The tax this family exists to measure
///
/// `general.md:56` prioritisation picks `dup` over `ins` and a copy count over `dup`,
/// and ferro applies it **per member**, inside `normalize_na_edit`, before the stages
/// that must compute on the result. Four passes in `src/normalize/merge.rs` exist
/// only to undo that decision once sibling context arrives —
/// `lower_repeat_edits`, `demote_repeats_spanning_siblings`,
/// `demote_coincident_tract_repeats` and `respell_colliding_duplications` — and every
/// one of them is reached only when a repeat or a duplication shares an allele with a
/// sibling. The six shapes below are one apiece plus the two halves of the pair that
/// separates an **authored** repeat from a **ferro-minted** one, which is the
/// distinction a change to *where* the mint happens moves.
///
/// # The zero this was written to close, and what closed it first
///
/// On `f7d177b5`, over the then-85,642-row corpus: **2,098 rows carry a repeat beside
/// a sibling in the OUTPUT and 0 carry one in the INPUT.** The minted half was covered
/// and the authored half was structural.
///
/// **That zero is no longer the base's, and quoting it as though it were would repeat
/// the mistake this family exists to prevent.** Re-measured on `5567412f`, over the
/// 95,614-row corpus of 23 shape families: the input side reads **756**, every one of
/// them from `repeat_beside_a_sibling`, which #1752 added independently while this
/// branch was open. So the structural zero was closed before this branch's base, and
/// what this family adds beside it is a sibling placed *outside* the tract rather than
/// intersecting it — the two are measured against each other in `inputs_for`'s note on
/// this family, which is why both earn a place.
///
/// `the_corpus_feeds_a_repeat_beside_a_sibling` is the guard, and it asserts the
/// property — that some row feeds one in — rather than either figure.
const TRACT_FAMILY: (&str, &str) = (
    "repeat_beside_a_sequence_sibling",
    "#1946 — a repeat or duplication over a tandem tract, sharing its allele with a \
     sibling inside or beside the tract",
);

/// Length of a `repeat_beside_a_sequence_sibling` core.
///
/// The same 20 as [`REVCOMP_CORE_LEN`] and [`corpus_sequences`], and for the same
/// reason: `CDS_START`/`CDS_END` and `EXON_SPANS` are tuned to a 20-base core, so
/// matching it is what lets these designs run on all three DNA axes.
const TRACT_CORE_LEN: usize = 20;

/// Repeat units the designs are built from.
///
/// **Three bases, never one or two, and that is a spec constraint rather than a
/// preference.** `repeated.md:21` permits a repeat description on a `c.` reference
/// only for units whose length is a multiple of three — it calls `c.2686A[10]` out by
/// name — so a mononucleotide or dinucleotide unit would make every `c.`/`cx.` row of
/// the two *authored* repeat shapes spec-invalid, and a corpus that emits invalid
/// input measures nothing. `first_tandem_triplet` is three bases for exactly this
/// reason.
///
/// Two units rather than one because the tract's own internal structure is a
/// different question from its length: `ACG` has three distinct bases, so no rotation
/// of it tiles anything but itself, while `AAT` carries a homopolymer that the
/// rotation search in `rules::insertion_to_duplication` has to get past.
const TRACT_UNITS: &[&str] = &["ACG", "AAT"];

/// Copies of the unit **present in the reference**.
///
/// Three and four, so that a shape spelling `unit[copies - 1]` is a genuine shrink to
/// two or three copies and never degenerates to a single copy — `repeated.md` defines
/// a repeat as a unit "present several times", so `unit[1]` would describe something
/// the notation does not mean.
const TRACT_COPIES: &[usize] = &[3, 4];

/// Offsets into the core at which a tract is planted.
///
/// Two, chosen for what the *axis* does at them rather than for the tract: at 1 the
/// tract starts inside the `cx` reference's 5'UTR (`CDS_START_MULTI` is 4), and at 5
/// the four-copy tract runs past `CDS_END` into the 3'UTR. Both cross at least one
/// `EXON_SPANS` junction. The tract is identical at both, so a difference between
/// them is a difference the axis made — the same argument [`REVCOMP_SPAN_OFFSETS`]
/// makes.
///
/// Never 0: shape `two_junctions_grow_one_tract` inserts at the tract's **5'**
/// junction, which needs a base before the tract to anchor against.
const TRACT_OFFSETS: &[usize] = &[1, 5];

/// One `repeat_beside_a_sequence_sibling` design: a reference and the tandem tract
/// planted in it.
struct TractDesign {
    /// The dump's `reference` column. A label rather than the core, matching
    /// [`RevcompDesign`] and [`long_corpus_sequences`].
    label: String,
    /// The 20-base reference sequence.
    core: String,
    /// The repeated unit, three bases.
    unit: &'static str,
    /// Copies of `unit` present in `core`.
    copies: usize,
    /// 0-based offset of the tract's first base.
    start: usize,
    /// 0-based offset of the tract's last base — `start + 3 * copies - 1`.
    end: usize,
}

/// The full design set: every unit at every copy count at every offset.
fn tract_designs() -> Vec<TractDesign> {
    let mut designs = Vec::new();
    for &unit in TRACT_UNITS {
        for &copies in TRACT_COPIES {
            for &start in TRACT_OFFSETS {
                let end = start + unit.len() * copies - 1;
                // The widest shape reads two offsets past the tract (a sibling at
                // `end + 2`), so a design that did not leave room would silently
                // emit a description addressing a base that is not there.
                assert!(
                    end + 2 < TRACT_CORE_LEN,
                    "{unit}x{copies} at {start} leaves no room for a sibling 3' of the tract"
                );
                designs.push(TractDesign {
                    label: format!("tract_{unit}x{copies}_at{start}"),
                    core: tract_core(unit, copies, start, TRACT_CORE_LEN),
                    unit,
                    copies,
                    start,
                    end,
                });
            }
        }
    }
    designs
}

/// A design's core: filler everywhere, with `copies` copies of `unit` planted at
/// `start`, and the two tract boundaries repaired so the tract cannot extend.
///
/// The filler is period-4 `GCTA`, the same one [`revcomp_core`] uses and for the same
/// reason — one rotation out of phase with the `ACGT` pad [`padded`] wraps a `g.` core
/// in, so the pad's rotation cannot continue through the flank.
///
/// **The repair is the load-bearing part.** Filler alone does not bound a tract: a
/// `GCTA` flank ending in `A` sits immediately 3' of an `ACGACG` tract and continues
/// it, so the real tract would be longer than the design claims and "inside the tract"
/// and "beside it" would stop being different places. Both boundary positions are
/// therefore forced away from the tract's periodic extension —
/// `the_designed_tracts_are_exactly_as_long_as_they_claim` re-derives that from the
/// finished string rather than trusting this construction.
fn tract_core(unit: &str, copies: usize, start: usize, core_len: usize) -> String {
    let mut core: Vec<u8> = (0..core_len).map(|i| b"GCTA"[i % 4]).collect();
    let tract = unit.repeat(copies);
    core[start..start + tract.len()].copy_from_slice(tract.as_bytes());
    let unit_bytes = unit.as_bytes();
    // 5' boundary: the tract's backward extension is the unit's last base.
    core[start - 1] = flank_base(&core, start - 1, unit_bytes[unit_bytes.len() - 1]);
    // 3' boundary: its forward extension is the unit's first base.
    let after = start + tract.len();
    if after < core_len {
        core[after] = flank_base(&core, after, unit_bytes[0]);
    }
    String::from_utf8(core).expect("GCTA and the designed units are valid UTF-8")
}

/// A base for the flank position `at`, avoiding `forbidden` and the two degeneracies
/// that would weaken a design silently.
///
/// `forbidden` is the base that would let the tract extend through this position.
/// Beyond that the choice must not equal the flank base on the far side — a two-base
/// run in the flank is a tract of its own, and a 3' shift out of the designed tract
/// could then keep going — and it must respect [`padded`]'s two edge rules: a core
/// starting with `A` or ending with `T` extends the `ACGT` pad's own rotation.
fn flank_base(core: &[u8], at: usize, forbidden: u8) -> u8 {
    let neighbour = if at == 0 { None } else { Some(core[at - 1]) };
    *b"GCTA"
        .iter()
        .find(|&&candidate| {
            candidate != forbidden
                && Some(candidate) != neighbour
                && !(at == 0 && candidate == b'A')
                && !(at == core.len() - 1 && candidate == b'T')
        })
        .expect("four candidates cannot all be excluded by three constraints")
}

/// The six spellings one tract design contributes on `axis`.
///
/// Each names one pass of the tax in [`TRACT_FAMILY`]'s doc, and the first four are
/// deliberately two matched pairs — the same geometry authored and minted, and the
/// same repeat with its sibling inside the tract and beside it. A change that moves
/// *when* a presentational form is chosen moves one member of a pair and not the
/// other; a change that moves the choice itself moves both.
///
/// The sibling is a substitution throughout, so the only thing varying between
/// shapes is where it sits and what the repeat member is. A richer sibling would be a
/// second variable in a family that already has six shapes and eight references.
fn tract_inputs_for(axis: &str, design: &TractDesign) -> Vec<String> {
    let prefix = prefix_for(axis);
    let bytes = design.core.as_bytes();
    let p = |i: usize| hgvs_pos(axis, i);
    let sub = |i: usize| {
        let base = bytes[i] as char;
        let other = if base == 'A' { 'C' } else { 'A' };
        format!("{}{base}>{other}", p(i))
    };
    let (unit, start, end) = (design.unit, design.start, design.end);
    let span = format!("{}_{}", p(start), p(end));
    // One unchanged base between the tract and the sibling. Zero would be adjacent
    // and is the merge case, which `separated_revcomp_runs` already sweeps; one is
    // the separation at which the two members must stay two.
    let beside = sub(end + 2);
    let inside = sub(start + 1);
    vec![
        // `lower_repeat_edits`: an authored repeat inside a cis allele reaches it
        // without ferro having minted anything, which is why that pass is
        // relocatable to the input boundary and not deletable.
        format!("{prefix}[{span}{unit}[{}];{beside}]", design.copies + 1),
        // `demote_repeats_spanning_siblings`, authored: the repeat's tract span
        // covers the sibling's base.
        format!("{prefix}[{span}{unit}[{}];{inside}]", design.copies + 1),
        // The same pass reached through its `RepeatSource::Removed` arm — a repeat
        // that shrank its tract re-spells as a deletion, not as a duplication.
        format!("{prefix}[{span}{unit}[{}];{inside}]", design.copies - 1),
        // `demote_repeats_spanning_siblings`, minted: an insertion of one unit at
        // the tract's 3' junction is the same variant as the first shape, spelled
        // so that ferro rather than the author chooses the repeat.
        format!("{prefix}[{}_{}ins{unit};{inside}]", p(end), p(end + 1)),
        // Two members that each grow the SAME tract, one at each of its
        // junctions — the input shape `demote_coincident_tract_repeats` (#1316)
        // describes.
        //
        // **It does not reach that pass, and the comment says so rather than
        // claiming a reach it does not have.** Measured by instrumenting the
        // pass's own group loop and dumping the corpus: it builds a candidate
        // group **7,413** times over 85,930 rows and **every one is of size 1**,
        // so the `at.len() < 2` guard always continues and the pass never
        // rewrites anything. This shape misses because the sequence-first cis
        // collapse merges the two junction insertions into one member *before*
        // any per-member mint — `g.[261_262insACG;273_274insACG]` comes out
        // `g.262_273ACG[6]`, one member, so there is never a pair to demote. A
        // second candidate was probed and also missed for its own reason: a
        // tract-shrinking deletion beside a tract-growing insertion cancels to
        // `g.274_275=`.
        //
        // Reaching the pass needs the cis collapse to DECLINE while both members
        // still mint as repeats over one span, which is #1946's item 1 — a
        // non-literal payload dropping the whole allele — and is not reachable
        // from an authored input while `try_expand_genome_ins` expands those at
        // the boundary. The shape is kept because it is the geometry the pass
        // documents and it pins the collapse that pre-empts it; what is not kept
        // is the claim that it exercises the pass.
        format!(
            "{prefix}[{}_{}ins{unit};{}_{}ins{unit}]",
            p(start - 1),
            p(start),
            p(end),
            p(end + 1)
        ),
        // `respell_colliding_duplications`: a duplication whose span holds the
        // sibling's base. A `dup` claims no bases but does name a span, which is
        // the discrepancy that pass repairs.
        format!("{prefix}[{span}dup;{inside}]"),
    ]
}

/// The family drawn against a **long transcript with a long tandem tract** (#1946).
///
/// Every other family in this file is drawn on a 20-base core, and that turns out to
/// hide the one question a render-stage relocation of the repeat mint has to answer.
///
/// # The blindness this closes, measured
///
/// The three repeat-minting rules are handed a `ref_seq` whose provenance differs by
/// axis: on `c.`/`r.`/`n.` it is the **entire spliced transcript**, with no window and
/// therefore no edge, while on `g.` it is a window (`NormalizeConfig::window_size`,
/// 100 each side) that *grows* when a shuffle runs to its edge. Tract discovery walks
/// outward bounded only by `ref_seq.len()`, so a stage that substituted one fixed
/// padded window would find different tract extents near an edge on the transcript
/// axes. That is a behavioural difference, not plumbing.
///
/// Instrumented over the corpus **as it stood before this family**, the question was
/// unanswerable in both directions and neither was visible from the numbers:
///
/// - **Transcript axes: `ref_seq.len()` was 20, always** — the whole synthetic
///   transcript. Any window a relocation would substitute (100, or `CANONICAL_PAD`'s
///   128) is 5-13x *larger* than the entire reference, so the window is the transcript
///   and the two strategies are the same thing by construction.
/// - **Genomic axis: 0 of 3,361 insertion firings had a tract touching the window
///   edge.** The growing-window loop — the entire reason the genomic path differs —
///   never fired for a tract.
///
/// And the reassuring half was equally corpus-shaped: the largest distance any tract
/// walk reached from its anchor was **11 bases**, which is a property of 20-base cores
/// and 12-base designed tracts, not of the algorithm. Real references carry
/// homopolymers and STRs of hundreds of bases.
///
/// # The design decisions, made rather than defaulted
///
/// **Tract lengths straddle both window constants, and the family asserts it.**
/// `the_long_tracts_bracket_the_normalizer_window` requires the set to hold a tract
/// shorter than `window_size` and one longer than twice it, so a change to either the
/// copy counts or the constant fails rather than silently producing a family whose
/// tracts all sit on one side. The five lengths are 96, 102, 126, 129 and 300 bases —
/// two either side of 100, two either side of 128, and one far beyond both.
///
/// **The `c.` axis is not optional here; it is the axis the question is about.** It
/// is where `ref_seq` is the whole transcript, so a long *core* is what makes a long
/// *transcript*. `cx` is skipped for the same reason [`long_corpus_sequences`] skips
/// it: `EXON_SPANS` is a fixed 7/7/6 split of a 20-base core and means nothing here.
///
/// **These sit well below `MAX_CANONICAL_WINDOW` (4096), deliberately.** The longest
/// tract is 300 bases in an 800-base core, so the canonical-window cap and
/// `MAX_SPLIT_BLOCK` never engage and cannot confound the reading. Block-size
/// behaviour at the cap is `long_block_inversion`'s question and is measured there;
/// this family must not answer two questions at once.
const LONG_TRACT_FAMILY: (&str, &str) = (
    "long_tract_window_provenance",
    "#1946 — a tandem tract longer than a normalizer window, on a transcript served whole",
);

/// The repeated unit. Three bases, for [`TRACT_UNITS`]' reason: `repeated.md:21`
/// admits a repeat on a `c.` reference only for units whose length is a multiple of
/// three, and `c.` is the axis this family exists to reach.
const LONG_TRACT_UNIT: &str = "ACG";

/// Copies of the unit in the reference, chosen so the tract lengths bracket both
/// window constants — 96, 102, 126, 129 and 300 bases.
const LONG_TRACT_COPIES: &[usize] = &[32, 34, 42, 43, 100];

/// Length of a `long_tract_window_provenance` core.
///
/// Comfortably more than twice the longest tract, so every design has real flank on
/// both sides and no tract is bounded by the core rather than by its own periodicity.
const LONG_TRACT_CORE_LEN: usize = 800;

/// Offset of the tract within the core. Far enough in that the 5' flank alone exceeds
/// both window constants, so a window centred on the tract's 5' junction is bounded by
/// the window and not by the start of the reference.
const LONG_TRACT_START: usize = 200;

/// The design set: one per copy count.
fn long_tract_designs() -> Vec<TractDesign> {
    LONG_TRACT_COPIES
        .iter()
        .map(|&copies| {
            let end = LONG_TRACT_START + LONG_TRACT_UNIT.len() * copies - 1;
            assert!(
                end + 2 < LONG_TRACT_CORE_LEN,
                "a {copies}-copy tract leaves no 3' flank"
            );
            TractDesign {
                label: format!("longtract_{}x{copies}", LONG_TRACT_UNIT),
                core: tract_core(
                    LONG_TRACT_UNIT,
                    copies,
                    LONG_TRACT_START,
                    LONG_TRACT_CORE_LEN,
                ),
                unit: LONG_TRACT_UNIT,
                copies,
                start: LONG_TRACT_START,
                end,
            }
        })
        .collect()
}

/// The two spellings one long-tract design contributes on `axis`.
///
/// Both are insertions at a tract junction of a payload holding **two** copies of the
/// unit, which is the shape whose mint (`rules::insertion_to_repeat`) walks the tract
/// to discover its extent — the walk this family exists to make longer than a window.
/// One at each junction, so the walk is exercised in both directions: an insertion at
/// the 3' junction reaches the full tract length 5'-wards, and vice versa.
///
/// # Two copies, not one, and that is a hard requirement rather than a preference
///
/// `insertion_to_repeat` opens with
///
/// ```text
/// let added_copies = (inserted_seq.len() / base_unit.len()) as u64;
/// if added_copies < 2 {
///     // Single-copy addition is a duplication, not repeat notation.
///     return None;
/// }
/// ```
///
/// so a payload of one unit can never mint a repeat however long the tract is — it
/// takes the `ins`->`dup` path instead. The first cut of this family inserted a single
/// unit and was measured **declining on every one of its 40 rows**, producing a
/// populated family that exercised none of the machinery it was built for. It read as
/// a result (`g.550_552dup`, a perfectly reasonable output) rather than as a miss,
/// which is the failure this file exists to prevent one level up.
///
/// Deliberately lone members rather than cis alleles. The question is about what
/// reference a *per-member* mint is handed, and a sibling would add cis machinery that
/// has nothing to do with it.
fn long_tract_inputs_for(axis: &str, design: &TractDesign) -> Vec<String> {
    let prefix = prefix_for(axis);
    let p = |i: usize| hgvs_pos(axis, i);
    let payload = design.unit.repeat(2);
    vec![
        format!(
            "{prefix}{}_{}ins{payload}",
            p(design.end),
            p(design.end + 1)
        ),
        format!(
            "{prefix}{}_{}ins{payload}",
            p(design.start - 1),
            p(design.start)
        ),
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
            // regardless, so on a long core the transcript declares the core's
            // full length while its exon table maps 20 and `spliced_genomic` emits
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
    // The long-tract designs. `cx` is skipped for the reason `long_corpus_sequences`
    // skips it — `EXON_SPANS` is a fixed 7/7/6 split of a 20-base core — and the `c.`
    // axis is the point of the family rather than an extra: it is where `ref_seq` is
    // the whole transcript rather than a window (#1946).
    for design in long_tract_designs() {
        for (axis, direction, dir_label) in axes_and_directions() {
            if axis == "cx" {
                continue;
            }
            let normalizer = normalizer_for(axis, &design.core, direction);
            for input in long_tract_inputs_for(axis, &design) {
                let output = normalize_through(&normalizer, &input);
                rows.push(Row {
                    reference: design.label.clone(),
                    core: design.core.clone(),
                    axis,
                    direction: dir_label,
                    family: LONG_TRACT_FAMILY.0,
                    was_fixed_point: output == input,
                    input,
                    output,
                });
            }
        }
    }
    // The tandem-tract designs, drawn against their own references for the same
    // reason the two loops above are: a tract of a chosen unit at a chosen offset,
    // bounded so it cannot extend, is a property of the sequence and cannot be
    // asked of a random core (#1946).
    for design in tract_designs() {
        for (axis, direction, dir_label) in axes_and_directions() {
            let normalizer = normalizer_for(axis, &design.core, direction);
            for input in tract_inputs_for(axis, &design) {
                let output = normalize_through(&normalizer, &input);
                rows.push(Row {
                    reference: design.label.clone(),
                    core: design.core.clone(),
                    axis,
                    direction: dir_label,
                    family: TRACT_FAMILY.0,
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
    // The protein axis, appended last so every row above keeps its position in the
    // dump (#1606).
    //
    // **One direction, not two, and that is measured rather than assumed.** Every
    // loop above crosses its cores with `axes_and_directions`, which carries a 3'
    // and a 5' cell. The protein normalization path never reads
    // `NormalizeConfig::direction` — a `p.` row comes out identical under both — so
    // dumping the second direction would double the protein half of the corpus with
    // byte-identical rows and make a reader believe the axis had been measured under
    // 5' shifting when nothing had. `the_protein_axis_is_direction_invariant` pins
    // that, and it is a live guard rather than a comment: if the protein path ever
    // starts honouring the direction, it fails and says to add the cell here.
    for peptide in protein_sequences(seeds) {
        let normalizer = Normalizer::with_config(
            protein_provider(&peptide),
            NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
        );
        for (family, _) in PROTEIN_FAMILIES {
            for input in protein_inputs_for(family, &peptide) {
                let output = normalize_through(&normalizer, &input);
                rows.push(Row {
                    reference: peptide.clone(),
                    // The protein axis brings no designed references, so the
                    // `reference` column is the peptide itself rather than a
                    // label — `core` is therefore equal to it, which is the
                    // ordinary case the field documents (#1624/#1625).
                    core: peptide.clone(),
                    axis: "p",
                    direction: "3prime",
                    family,
                    was_fixed_point: output == input,
                    input,
                    output,
                });
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

/// Deterministic 20-residue peptides, two per seed, each starting with `Met`.
///
/// Two per seed mirrors [`corpus_sequences`] — one drawn from a two-residue
/// alphabet and one from the full table — so the protein half of the corpus scales
/// with `--seeds` exactly as the DNA half does and stays prefix-stable.
///
/// **Runs are planted, not hoped for.** The residues are drawn in runs of one to
/// three rather than one at a time, and that is the whole design. Every family
/// below is about *shifting*, and a shift only has somewhere to go inside a tract
/// of equal residues; drawing residue-by-residue from a 6-letter alphabet gives a
/// run of 3 about once per 36 positions, so most cores would be fixed points for
/// the two shift families and the corpus would report a near-zero it had generated
/// rather than measured. This is the same lesson as `#1517` — a random core cannot
/// be *asked* to have a structural property, so the property is built in.
fn protein_sequences(seeds: u32) -> Vec<String> {
    let drawn = &RESIDUE_CODES[1..];
    let mut peptides = Vec::with_capacity(2 * seeds as usize);
    for seed in 0..seeds {
        for alphabet in [&drawn[..2], drawn] {
            let mut state = u64::from(seed).wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            let mut peptide = String::with_capacity(PROTEIN_LEN);
            peptide.push('M');
            while peptide.len() < PROTEIN_LEN {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let residue = alphabet[(state % alphabet.len() as u64) as usize].0;
                let run = 1 + (state >> 32) % 3;
                for _ in 0..run {
                    if peptide.len() == PROTEIN_LEN {
                        break;
                    }
                    peptide.push(residue);
                }
            }
            peptides.push(peptide);
        }
    }
    peptides
}

/// Three-letter code of the residue at 1-based `position` of `peptide`.
///
/// Panics on a residue outside [`RESIDUE_CODES`] rather than rendering something
/// unparseable: cores are generated from that table, so a miss means the generator
/// and the table have drifted apart, and a description built from a wrong code
/// would be recorded as a parse error — a corpus row that looks like a measurement
/// and is really a bug in the corpus.
fn residue_at(peptide: &str, position: usize) -> &'static str {
    let one = peptide
        .as_bytes()
        .get(position - 1)
        .map(|b| *b as char)
        .unwrap_or_else(|| panic!("position {position} is past the end of a {PROTEIN_LEN}-mer"));
    RESIDUE_CODES
        .iter()
        .find(|(letter, _)| *letter == one)
        .map(|(_, three)| *three)
        .unwrap_or_else(|| panic!("residue {one} is not in RESIDUE_CODES"))
}

/// A three-letter residue that is **not** the one at `position`, for building an
/// edit whose payload genuinely changes the reference.
///
/// A payload equal to the reference residue would make the edit a no-op that
/// normalizes to `=`, which is a different family's question. The table has seven
/// entries, so a differing one always exists.
fn residue_other_than(peptide: &str, position: usize) -> &'static str {
    let same = residue_at(peptide, position);
    RESIDUE_CODES
        .iter()
        .skip(1)
        .map(|(_, three)| *three)
        .find(|three| *three != same)
        .expect("the residue table holds more than one entry")
}

/// The inputs one protein `family` contributes for one peptide core.
///
/// Every position is read **off the peptide** rather than assumed, because
/// normalization checks the spelled residue against the reference and rejects a
/// mismatch (`Amino acid mismatch at position N`). A hardcoded residue would turn
/// the whole family into rejected rows that still look like a populated corpus.
///
/// All spans start at residue 2 or later. Residue 1 is `Met`, and an edit that
/// reaches it emits the illegal start-loss spelling of #1607 — a defect worth its
/// own guard, not worth silently seeding thousands of corpus rows with.
fn protein_inputs_for(family: &str, peptide: &str) -> Vec<String> {
    let prefix = format!("{PROTEIN_ACCESSION}:p.");
    let at = |p: usize| format!("{}{}", residue_at(peptide, p), p);
    let mut out = Vec::new();
    match family {
        "protein_shift_del" => {
            for start in 2..PROTEIN_LEN {
                out.push(format!("{prefix}{}del", at(start)));
                out.push(format!("{prefix}{}_{}del", at(start), at(start + 1)));
            }
        }
        "protein_shift_dup" => {
            for start in 2..PROTEIN_LEN {
                out.push(format!("{prefix}{}dup", at(start)));
                out.push(format!("{prefix}{}_{}dup", at(start), at(start + 1)));
            }
        }
        "protein_ins_becomes_dup" => {
            // The payload copies the residue 5' of the junction, so the insertion
            // denotes a duplication and the question is whether it is re-typed as
            // one. `general.md:56` ranks duplication above insertion.
            for start in 2..PROTEIN_LEN {
                out.push(format!(
                    "{prefix}{}_{}ins{}",
                    at(start),
                    at(start + 1),
                    residue_at(peptide, start)
                ));
            }
        }
        "protein_equal_length_delins" => {
            // Payload length equals span length and the interior residues repeat
            // the reference, so the span is unchanged in the middle and changed at
            // both ends — the shape #1606 splits. Widths 3 and 4 give one and two
            // unchanged interior residues respectively.
            // `start` is bounded per width rather than once for the widest, so a
            // width-3 span reaches the C-terminal residue instead of stopping one
            // short of it. Deriving the bound this way also means there is no
            // unreachable `end > PROTEIN_LEN` guard to mislead a later reader about
            // whether spans can run off the end — they cannot, by construction.
            for width in [3usize, 4] {
                for start in 2..=(PROTEIN_LEN - width + 1) {
                    let end = start + width - 1;
                    let mut payload = String::new();
                    payload.push_str(residue_other_than(peptide, start));
                    for interior in (start + 1)..end {
                        payload.push_str(residue_at(peptide, interior));
                    }
                    payload.push_str(residue_other_than(peptide, end));
                    out.push(format!("{prefix}{}_{}delins{payload}", at(start), at(end)));
                }
            }
        }
        "protein_cis_separated" => {
            // Two members with two unchanged residues between them: the protein
            // analogue of `split_vs_spanning_delins`, which is where the DNA axes
            // see most of their churn.
            for start in 2..(PROTEIN_LEN - 3) {
                let far = start + 3;
                out.push(format!(
                    "{prefix}[{}{};{}{}]",
                    at(start),
                    residue_other_than(peptide, start),
                    at(far),
                    residue_other_than(peptide, far)
                ));
                out.push(format!(
                    "{prefix}[{}del;{}{}]",
                    at(start),
                    at(far),
                    residue_other_than(peptide, far)
                ));
            }
        }
        other => unreachable!("unknown protein family {other}"),
    }
    out
}

/// The provider a protein row is normalized through.
///
/// `add_protein` and nothing else — no transcript, no CDS, no contig. That is why
/// the protein axis costs so much less per row than the coding axes, which rebuild
/// a `Transcript` per cell (see [`normalizer_for`]).
fn protein_provider(peptide: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_protein(PROTEIN_ACCESSION, peptide);
    provider
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
            // A ranged repeat paired with a sibling, in the three geometries the
            // unified write footprint distinguishes: a deletion intersecting the
            // tract, an insertion at a junction strictly interior to it, and a
            // `dup` writing at the same junction the repeat's expansion lands on.
            // Emits nothing for a core with no tandem triplet — see
            // `first_tandem_triplet`.
            //
            // **Both directions of the repeat, which is the property the
            // footprint keys on.** `repeat_footprint` answers `at_junction(end)`
            // when the repeat grows and `spanning(removed_from, end, false)`
            // when it shrinks, and those are different branches with different
            // collision rules — so a family emitting only one of them would make
            // any `0 moved` over the other a structural zero. The tract is
            // **6** bases with a **3**-base unit, so `[3]` (9 bases) always
            // grows and `[1]` (3 bases) always shrinks: the count is chosen to
            // straddle the tract length rather than to look plausible. Each
            // sibling geometry is therefore emitted twice, once per direction.
            "repeat_beside_a_sibling" => {
                if s != 0 {
                    continue;
                }
                if let Some(t) = first_tandem_triplet(core) {
                    let unit = &core[t..t + 3];
                    // The repeat's tract is `t ..= t + 5`. Growing, its junction
                    // is `t + 5`; shrinking, it rewrites `t + 3 ..= t + 5` and
                    // leaves `t ..= t + 2` as untouched reference.
                    for count in [3, 1] {
                        out.push(format!(
                            "{prefix}[{}_{}{unit}[{count}];{}_{}del]",
                            p(t),
                            p(t + 5),
                            p(t + 4),
                            p(t + 7)
                        ));
                        out.push(format!(
                            "{prefix}[{}_{}{unit}[{count}];{}_{}insCC]",
                            p(t),
                            p(t + 5),
                            p(t + 2),
                            p(t + 3)
                        ));
                        out.push(format!(
                            "{prefix}[{}_{}{unit}[{count}];{}_{}dup]",
                            p(t),
                            p(t + 5),
                            p(t + 4),
                            p(t + 5)
                        ));
                    }
                }
            }
            // A lone spanning `delins` over four reference bases whose payload is
            // three: two novel bases followed by the reference base at `s + 2`, so
            // the payload's last base *coincides* with a reference base the span
            // still holds. That coincidence is the only interior match in the
            // block, and it is what an aligner cuts at — `DNA/delins.md:46`'s
            // construction, at its smallest.
            //
            // **Net deletion is the whole point of the family** (see `FAMILIES`).
            // Four against three, never four against four: an equal-length block
            // has a unique column correspondence and is a different regime, and a
            // net insertion is outside #1610's direction scope. Both are already
            // covered — by `delins_hiding_an_inversion` and by
            // `split_vs_spanning_delins` respectively.
            //
            // The two-member spelling of the same variant is emitted beside it, so
            // the family measures confluence as well as movement: the split writes
            // the deletion at `s + 3` explicitly rather than letting the derivation
            // place it.
            "lone_net_deletion_delins" => {
                let head = format!("{}{}", other(s), other(s + 1));
                let coincident = base(s + 2);
                out.push(format!(
                    "{prefix}{}_{}delins{head}{coincident}",
                    p(s),
                    p(s + 3)
                ));
                out.push(format!(
                    "{prefix}[{}_{}delins{head};{}del]",
                    p(s),
                    p(s + 1),
                    p(s + 3)
                ));
            }
            // An authored repeat that keeps its sibling — the shape `lower_repeat_edits`
            // actually serves, and the one `repeat_beside_a_sequence_sibling` does not
            // reach.
            //
            // # Why this exists alongside `repeat_beside_a_sibling`, measured
            //
            // #1752 added `repeat_beside_a_sibling` (#1749) independently and in this
            // same slot, citing the same structural zero. The two look redundant and
            // are not, and the difference is the sibling's **placement**: that family
            // puts a sibling that *intersects* the tract or its junction, this one puts
            // a substitution two bases **outside** it.
            //
            // Both reach `lower_repeat_edits` — ablating it moves 34 of that family's
            // 756 rows and 24 of this one's 240 — so reach alone does not separate
            // them. What separates them is what `--verify-spdi` reports:
            //
            // ```text
            // repeat_beside_a_sibling   c.[8_13ATA[3];12_15del] -> c.[12_15del;14_*1dup]
            // this family               g.[257_262TTT[1];264T>A] -> g.[256_265T[7];264T>A]
            // ```
            //
            // The first input is **already overlapping** — the deletion at `12_15`
            // intersects the tract at `8_13` — so it is malformed in and malformed out.
            // The second input is **disjoint**: the tract is `257_262` and the
            // substitution is at `264`, two bases clear of it. Normalization *creates*
            // the overlap, emitting a repeat whose span contains the sibling's base.
            //
            // Clean in, malformed out is a different and worse defect than malformed in,
            // malformed out, and it is the one #1983 records. This family is the only
            // one that builds it. Delete it and the reproducer goes with it.
            //
            // # This family found a live defect, and `--verify-spdi` is how
            //
            // Five of its rows report `SPDI-MISMATCH`: the normalized output denotes
            // different bases from the input. Four are `g.` and one `c.`, and the `g.`
            // ones are legible by eye —
            //
            // ```text
            // g.[257_262TTT[1];264T>A]  ->  g.[256_265T[7];264T>A]
            // g.[257_262CCC[1];264C>A]  ->  g.[257_264C[5];264C>A]
            // ```
            //
            // In both, the emitted repeat's tract span **contains the sibling's base**
            // (264 is inside 256_265, and inside 257_264), so the allele is
            // overlapping and malformed. That is exactly the shape
            // `demote_repeats_spanning_siblings` exists to prevent.
            //
            // It cannot prevent this one. That pass decides what to demote by matching
            // on the member's **pre-normalization** edit kind, and its arms are
            // `Deletion | Duplication | Insertion | Delins`, with `_ => None` and a
            // `continue`. An **authored** repeat was already a `NaEdit::Repeat` before
            // normalization, so it falls through and is skipped — the pass covers only
            // the repeats ferro itself mints. That is #1946's thesis as a measurement
            // rather than an argument: a per-shape undo cannot cover a shape it did
            // not create.
            //
            // Kept as a red row rather than removed. The family's value is that it
            // reaches this, and deleting the shape would restore the silence it was
            // added to break.
            //
            // # Why this is on the random cores and not on the designed ones
            //
            // Measured, by instrumenting `merge::canonicalize_from_sequence`. Of the 288
            // rows the designed tract family builds, **144 author a repeat inside a cis
            // allele and not one reaches `lower_repeat_edits`.** Fifty-two of them get as
            // far as `merge::canonicalize_from_sequence`, and **every one arrives with
            // `variants.len() == 1`** — the allele has already collapsed to a single
            // member, so the lone-repeat lockout (`merge.rs:3186`) refuses it. That
            // lockout is deliberate and load-bearing: ablating it moves 7 of 85,930 rows,
            // all `repeat_expansion`, in exactly the direction its comment predicts
            // (`g.257_262AAA[1]` goes from `g.257_263A[4]` to `g.257_259del`).
            //
            // So the gating property is neither unit length nor sibling placement, both
            // of which were tested and ruled out — a two-copy tract on a designed core
            // still reaches zero. It is **member survival**: the repeat has to still have
            // a sibling at the moment the sequence-first pass runs. The designed tracts
            // are 9 and 12 bases growing by 3, and the derived block swallows the
            // sibling; a six-base tract on a random core does not.
            //
            // This shape keeps two members and does reach the pass — measured at 19
            // firings of `lower_repeat_edits` over the 120 `g.`-shape rows it adds. It is
            // the geometry #1946's own reproducer used.
            //
            // Emitted once per core (`s != 0`), like `repeat_expansion`, and only where
            // the core has a tandem triplet with two bases to spare 3' of it — 20 of the
            // 48 cores, one fewer than `repeat_expansion`'s 21, and
            // `the_authored_repeat_family_covers_the_cores_with_room` accounts for the
            // one it drops rather than leaving it silent.
            "authored_repeat_beside_a_sibling" => {
                if s != 0 {
                    continue;
                }
                if let Some(t) = first_tandem_triplet(core) {
                    if t + 8 >= bytes.len() {
                        continue;
                    }
                    let unit = &core[t..t + 3];
                    let sibling = format!("{}{}>{}", p(t + 7), base(t + 7), other(t + 7));
                    out.push(format!(
                        "{prefix}[{}_{}{unit}[3];{sibling}]",
                        p(t),
                        p(t + 5)
                    ));
                    out.push(format!(
                        "{prefix}[{}_{}{unit}[1];{sibling}]",
                        p(t),
                        p(t + 5)
                    ));
                }
            }
            // A contiguous coalescible run beside a DISTANT cis member (#2192).
            // The run is a four-base LEFT ROTATION spelled in its fragmented
            // `[del; ins]` form: deleting `core[a]` and re-inserting it 3' of the
            // window denotes `core[a..a+4]` rotated one place, an equal-length,
            // all-columns-changed block that `coalesce_solid_run` collapses to one
            // `delins`. A distant substitution at the far end of the core stretches
            // the whole-variant hull across a wide gap, so baseline's whole-hull
            // collapse sees a non-contiguous changed-column set and declines —
            // keeping the rotation fragmented — while the per-run pass coalesces it.
            //
            // Built at a FIXED four-base window rather than the sliding `s..s+4`
            // one: the left rotation changes every column only when no cyclically
            // adjacent pair of the window is equal (`ATAT` qualifies with two
            // distinct bases, `AATA` does not), and the distant member needs a gap
            // wider than the five-base window can express. Emits nothing for a core
            // with no such window — reported honestly rather than as a zero (see the
            // family's comment in `FAMILIES`).
            "run_beside_a_distant_member" => {
                if s != 0 {
                    continue;
                }
                // A left rotation of `core[a..a+4]` changes column `j` iff
                // `core[a+j] != core[a+(j+1)%4]`, so every column changes exactly
                // when no cyclically adjacent pair is equal.
                let clean_rotation =
                    |a: usize| (0..4).all(|j| bytes[a + j] != bytes[a + (j + 1) % 4]);
                // Room for the four-base window and a distant member past a gap.
                if let Some(a) = (0..bytes.len().saturating_sub(6)).find(|&a| clean_rotation(a)) {
                    let dist = bytes.len() - 1;
                    out.push(format!(
                        "{prefix}[{}del;{}_{}ins{};{}{}>{}]",
                        p(a),
                        p(a + 3),
                        p(a + 4),
                        base(a),
                        p(dist),
                        base(dist),
                        other(dist)
                    ));
                    // The already-coalesced spelling of the same variant beside the
                    // same distant member. On `main` this ALSO fragments — the
                    // whole-hull collapse re-derives it as `[del; dup]` — so it is a
                    // second moving row, not a fixed point; its value is that the fix
                    // keeps an input already written as one `delins` coalesced rather
                    // than re-fragmenting it.
                    let rot: String = (1..4)
                        .map(|i| base(a + i))
                        .chain(std::iter::once(base(a)))
                        .collect();
                    out.push(format!(
                        "{prefix}[{}_{}delins{rot};{}{}>{}]",
                        p(a),
                        p(a + 3),
                        p(dist),
                        base(dist),
                        other(dist)
                    ));
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
/// `Transcript` — on the long cores that is a multi-kilobase copy per input, and
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

    /// A row of one of the eighteen sequence-keyed families, where the
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
    /// — 20,516 of 78,298 rows. Classifying a decline as `Different` would report
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
    /// The two views agree on eighteen of the twenty-two families, and that is
    /// precisely what made the verify pass's confusion of them survive review —
    /// a bug that only manifests on `long_block_inversion` and
    /// `separated_revcomp_runs` is a bug in 286 rows out of 78,298. So this pins
    /// the distinction itself, on both sides: where the column is a label the
    /// sequence must be the design's, and where it is not the two must be equal.
    #[test]
    fn every_row_carries_the_sequence_it_was_drawn_against() {
        let labelled: BTreeMap<String, String> = long_corpus_sequences()
            .into_iter()
            .chain(revcomp_designs().into_iter().map(|d| (d.label, d.core)))
            .chain(tract_designs().into_iter().map(|d| (d.label, d.core)))
            .chain(long_tract_designs().into_iter().map(|d| (d.label, d.core)))
            .collect();
        let mut label_keyed = 0usize;
        for row in dump(2) {
            // A label is not over `ACGT` — `revcomp_sep0_at0` and any
            // `nearpalindrome_<n>` both fail this — so it is the one predicate
            // that separates the two kinds of string without consulting the map.
            //
            // The alphabet is per-axis, and that is load-bearing rather than a
            // widening: a `p.` row is drawn against a **peptide**, so testing it
            // for `ACGT` would fail every protein row on the strength of the
            // corpus having grown a molecule type. Exempting them instead would
            // be worse — it would leave the protein half of the corpus with no
            // check that `core` holds a sequence at all — so each axis is
            // checked against the alphabet it is actually drawn from.
            let (alphabet_ok, alphabet) = if row.axis == "p" {
                (
                    row.core
                        .chars()
                        .all(|c| RESIDUE_CODES.iter().any(|(one, _)| *one == c)),
                    "one-letter residue codes",
                )
            } else {
                (row.core.bytes().all(|b| b"ACGT".contains(&b)), "ACGT")
            };
            assert!(
                !row.core.is_empty() && alphabet_ok,
                "[{}] {} carries something that is not a {} sequence: {}",
                row.family,
                row.reference,
                alphabet,
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

    /// The corpus must build a block past the normalizer's length gate (#1460).
    ///
    /// The gate is a *length* one — `canonicalize_from_sequence` pads the changed
    /// interval and refuses once the window exceeds `MAX_CANONICAL_WINDOW` — so a
    /// corpus of 20-mers takes the same path on both sides of any change to it, by
    /// construction. That is why #1403 measured `0 of 18,432` here and why the zero
    /// was guaranteed rather than informative. Fixing #1456 added the shapes that
    /// were missing; this adds the scale.
    ///
    /// # The bound is IMPORTED, and that is the whole point (#1925)
    ///
    /// This test used to open `const SPLIT_CAP: u64 = 1024;` with a comment saying
    /// to update it by hand if the constant moved. The constant moved — #1899 took
    /// `MAX_SPLIT_BLOCK` to 4096 — and this kept passing on `1100 > 1024` while
    /// both long cores sat below the real gate. **A guard that restates the value
    /// it guards cannot observe that value changing**, and a comment asking a
    /// future reader to notice is not a mitigation; it is the defect written down.
    ///
    /// It also now names the gate that fires on the shipping path.
    /// `MAX_SPLIT_BLOCK`'s only functional reader additionally requires
    /// `reference.len() != result.len()`, and this family's whole-core `inv` is
    /// equal-length, so it could never have tripped that cap at any size — the
    /// assertion was measuring a threshold this corpus cannot reach even when the
    /// number is right.
    ///
    /// # …and past the gate the family is UNIFORMLY `inv` (#1703)
    ///
    /// This section used to assert the opposite: that the two long cores land on
    /// **opposite sides** of the gate, one re-derived into substitutions below it
    /// and one kept as `inv` above it — a straddle. That straddle is gone, and it
    /// was removed by the ruling this corpus now measures.
    /// `rulings[whole-span-reverse-complement-types-as-inv]` (`DNA/inversion.md:5`,
    /// 2026-08-13, #1703) types an **exact** whole-span reverse complement as one
    /// `inv` UNIFORMLY, with no length or window bound: below the gate
    /// `coalesce_inversion_runs`' ungated whole-span check types it `inv` before
    /// any cut, and above the gate `canonicalize_from_sequence` declines on width
    /// and the per-member pipeline keeps the single `inv`. Both sides now agree,
    /// which is the whole point of the ruling. A near-palindrome core is an exact
    /// whole-span reverse complement, so this family can no longer straddle the
    /// gate by verdict, and the `2 * EDGE_PERTURBATION` offset in
    /// [`long_corpus_sequences`] no longer changes the inv-vs-subs answer (it stays
    /// `inv` either way). The third assertion below therefore pins the uniform
    /// verdict — every measured past-cap core is kept as `inv` — rather than the
    /// straddle, which the ruling made unexpressible for this family.
    #[test]
    fn the_corpus_emits_a_block_past_the_split_cap() {
        use std::collections::BTreeSet;

        let split_cap = MAX_CANONICAL_BLOCK as u64;
        let rows = dump(1);
        let genomic: Vec<&Row> = rows.iter().filter(|row| row.axis == "g").collect();
        let longest = genomic
            .iter()
            .map(|row| longest_span(&row.input))
            .max()
            .unwrap_or(0);
        assert!(
            longest > split_cap,
            "the longest block the corpus builds is {longest} bases, under the {split_cap}-base \
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
            .filter(|row| longest_span(&row.input) > split_cap)
            .collect();
        let measured: Vec<&Row> = past_cap
            .iter()
            .filter(|row| !SENTINELS.contains(&row.output.as_str()))
            .map(|row| **row)
            .collect();
        assert!(
            !measured.is_empty(),
            "all {} rows past the {split_cap}-base cap normalized to a sentinel, so the corpus \
             builds the scale but measures nothing at it: {:?}",
            past_cap.len(),
            past_cap.iter().map(|r| &r.output).collect::<Vec<_>>()
        );

        // Post-#1703 the verdict is UNIFORM, not a straddle: an exact whole-span
        // reverse complement is typed `inv` on both sides of the gate, so every
        // measured past-cap near-palindrome core keeps its `inv`. The verdict is
        // read off the output's spelling, which is how the dump reports it — a
        // core kept as `inv` ends its output with `inv` (trimmed to the changed
        // interval, but an `inv` still); a core re-derived into substitutions
        // would not. Per the ruling, none is re-derived.
        let refused: BTreeSet<&str> = measured
            .iter()
            .filter(|row| row.output.ends_with("inv"))
            .map(|row| row.reference.as_str())
            .collect();
        let derived: BTreeSet<&str> = measured
            .iter()
            .filter(|row| !row.output.ends_with("inv"))
            .map(|row| row.reference.as_str())
            .collect();
        assert!(
            !refused.is_empty() && derived.is_empty(),
            "per `rulings[whole-span-reverse-complement-types-as-inv]` every measured past-cap \
             core is an exact whole-span reverse complement and must keep its `inv`; kept as \
             `inv`: {refused:?}, but re-derived into substitutions (must be empty): {derived:?} — {:?}",
            measured
                .iter()
                .map(|row| (&row.reference, &row.input, &row.output))
                .collect::<Vec<_>>()
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

    /// Every protein family emits something, so none can contribute a silent zero.
    ///
    /// The sibling of `every_family_emits_rows_on_both_axes`, and the reason it is
    /// separate is that the protein families are not crossed with the DNA axes —
    /// there is one axis here, so there is one loop.
    #[test]
    fn every_protein_family_emits_rows() {
        let peptide = &protein_sequences(1)[0];
        for (family, _) in PROTEIN_FAMILIES {
            assert!(
                !protein_inputs_for(family, peptide).is_empty(),
                "protein family {family} emitted no rows"
            );
        }
    }

    /// The protein corpus scales and truncates like the DNA one.
    #[test]
    fn the_protein_corpus_is_prefix_stable() {
        let few = protein_sequences(3);
        let many = protein_sequences(9);
        assert_eq!(few, many[..few.len()]);
        assert_eq!(few.len(), 6, "two peptides per seed");
        assert!(
            few.iter()
                .all(|p| p.len() == PROTEIN_LEN && p.starts_with('M')),
            "every core is a {PROTEIN_LEN}-mer starting at Met"
        );
    }

    /// **The pin that justifies dumping the protein axis at one direction.**
    ///
    /// Every other loop in [`dump`] crosses its cores with `axes_and_directions`,
    /// which carries a 3' and a 5' cell. The protein path does not read
    /// `NormalizeConfig::direction` at all, so the second cell would be a
    /// byte-identical copy of the first — 7,344 more rows asserting nothing, and
    /// worse, a dump that *looks* like it measured the axis under 5' shifting.
    ///
    /// This is a live guard, not a comment restating an assumption. If protein
    /// normalization ever starts honouring the direction, this fails, and the fix
    /// is to add the 5' cell to the protein loop rather than to relax the test.
    ///
    /// Measured over every family on two cores rather than argued from the source,
    /// because "the code does not mention `direction`" is exactly the kind of claim
    /// that stops being true without anyone noticing.
    #[test]
    fn the_protein_axis_is_direction_invariant() {
        let mut compared = 0;
        for peptide in protein_sequences(2) {
            let provider = protein_provider(&peptide);
            let three = Normalizer::with_config(
                provider.clone(),
                NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
            );
            let five = Normalizer::with_config(
                provider,
                NormalizeConfig::default().with_direction(ShuffleDirection::FivePrime),
            );
            for (family, _) in PROTEIN_FAMILIES {
                for input in protein_inputs_for(family, &peptide) {
                    let a = normalize_through(&three, &input);
                    let b = normalize_through(&five, &input);
                    assert_eq!(
                        a, b,
                        "{input} differs by shuffle direction on the protein axis — \
                         the protein loop in `dump` must now carry both directions"
                    );
                    compared += 1;
                }
            }
        }
        assert!(
            compared > 500,
            "only {compared} rows compared — the families stopped emitting"
        );
    }

    /// No protein row is a parse error, a decline or a panic.
    ///
    /// This is the guard on the thing most likely to go quietly wrong here.
    /// Normalization checks each spelled residue against the reference and rejects
    /// a mismatch, so one wrong three-letter code turns a whole family into
    /// `<declined>` rows — a populated-looking corpus measuring nothing. It would
    /// not fail any other test: the rows exist, the key is unique, and the family
    /// is non-vacuous.
    #[test]
    fn no_protein_input_is_rejected() {
        for peptide in protein_sequences(2) {
            let normalizer = Normalizer::with_config(
                protein_provider(&peptide),
                NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
            );
            for (family, _) in PROTEIN_FAMILIES {
                for input in protein_inputs_for(family, &peptide) {
                    let output = normalize_through(&normalizer, &input);
                    assert!(
                        !output.starts_with('<'),
                        "{family}: {input} on {peptide} gave {output}"
                    );
                }
            }
        }
    }

    /// The cores carry residue runs long enough for a shift to move through.
    ///
    /// The generator plants runs rather than drawing residue-by-residue, and this
    /// is what pins that. Without it a later "simplification" of
    /// [`protein_sequences`] to a per-residue draw would leave the two shift
    /// families almost entirely fixed points, and the corpus would report a
    /// near-zero it had generated rather than measured — the corpus-blindness class
    /// one level down.
    #[test]
    fn protein_cores_carry_runs_a_shift_can_move_through() {
        for peptide in protein_sequences(4) {
            let bytes = peptide.as_bytes();
            let longest = bytes
                .iter()
                .enumerate()
                .map(|(i, b)| bytes[i..].iter().take_while(|other| *other == b).count())
                .max()
                .unwrap_or(0);
            assert!(
                longest >= 3,
                "{peptide}'s longest residue run is {longest}, so the shift families \
                 have almost nowhere to move on it"
            );
        }
    }

    /// The measured floor behind the bound above, so the bound is not folklore.
    ///
    /// Over the whole 48-core default corpus the shortest longest-run is **3** and
    /// the longest is 16. The `>= 3` bound is therefore the real floor rather than a
    /// round number picked to pass: a first cut asserted `>= 2`, which every core
    /// clears by a margin and which a regression to per-residue drawing would very
    /// likely clear too, since a two-residue alphabet throws adjacent pairs by
    /// chance. Asserting the distribution here keeps the two numbers honest — if
    /// the generator changes, this fails with what it actually produced.
    #[test]
    fn the_planted_runs_are_longer_than_chance_would_give() {
        let longest_run = |peptide: &str| -> usize {
            let bytes = peptide.as_bytes();
            bytes
                .iter()
                .enumerate()
                .map(|(i, b)| bytes[i..].iter().take_while(|other| *other == b).count())
                .max()
                .unwrap_or(0)
        };
        let cores = protein_sequences(24);
        assert_eq!(cores.len(), 48, "the default corpus is 24 seeds");
        let floor = cores
            .iter()
            .map(|p| longest_run(p))
            .min()
            .expect("non-empty");
        let ceiling = cores
            .iter()
            .map(|p| longest_run(p))
            .max()
            .expect("non-empty");
        assert_eq!(floor, 3, "measured floor moved");
        assert_eq!(ceiling, 16, "measured ceiling moved");
    }

    /// No protein family constructs an edit reaching residue 1.
    ///
    /// Residue 1 is `Met`, and an edit that reaches it emits the illegal start-loss
    /// spelling filed as #1607. Seeding the corpus with that would mean thousands
    /// of rows whose movement is a known bug rather than a representation choice.
    #[test]
    fn no_protein_family_touches_residue_one() {
        for peptide in protein_sequences(2) {
            for (family, _) in PROTEIN_FAMILIES {
                for input in protein_inputs_for(family, &peptide) {
                    assert!(
                        !input.contains("Met1"),
                        "{family} built {input}, which reaches residue 1"
                    );
                }
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

    /// Whether `s` spells a repeat.
    ///
    /// HGVS writes a repeat as `<unit>[<count>]`, so its bracket is preceded by a
    /// base; an allele's own bracket is preceded by the axis dot. That one
    /// character separates the two without parsing, which is what lets the checks
    /// below run over raw dump columns.
    fn holds_repeat_notation(s: &str) -> bool {
        ["A[", "C[", "G[", "T["]
            .iter()
            .any(|token| s.contains(token))
    }

    /// The corpus feeds a repeat **beside a sibling**, not only as a lone variant
    /// (#1946).
    ///
    /// `the_corpus_is_closed_over_its_output_vocabulary` above closes the corpus
    /// over the *constructs* it emits, and `repeat_expansion` is what closed the
    /// repeat row of that table. It closes nothing about **geometry**: every one of
    /// that family's inputs is a lone member, so the corpus asks what normalization
    /// does to a repeat that has no sibling and never what it does to one that has.
    ///
    /// That distinction is the whole subject of #1946. Four passes in
    /// `src/normalize/merge.rs` exist only to undo a repeat or a duplication a
    /// sibling-blind stage minted — `lower_repeat_edits`,
    /// `demote_repeats_spanning_siblings`, `demote_coincident_tract_repeats` and
    /// `respell_colliding_duplications` — and every one of them is reached only when
    /// a repeat and a sibling are in the same allele. Measured on `f7d177b5` over
    /// the then-85,642-row corpus: **2,098 rows carry a repeat beside a sibling in
    /// the OUTPUT and 0 carry one in the INPUT.** So the minted half was covered and
    /// the authored half was a structural zero — a change to where the mint happens
    /// could report `0 moved` for the authored path and mean only that the corpus
    /// could not build it.
    ///
    /// **Re-measured on `5567412f`, the input side reads 756 and the output side
    /// 2,410**, over 95,614 rows and 23 shape families. The zero was closed by
    /// #1752's `repeat_beside_a_sibling`, which supplies all 756. Both figures are
    /// recorded with the revision they were taken on precisely because neither is
    /// what this test checks.
    ///
    /// # Why this asserts a property and pins neither number
    ///
    /// Pinning either would make the count *be* the property, and this guard's whole
    /// subject is a count that was true and then silently stopped being so. What is
    /// asserted is that the corpus mints one **and** feeds one, with the minted side
    /// asserted non-zero first so `0 of 0` cannot pass as a result.
    ///
    /// The predicate is deliberately `contains(';')` plus a repeat token, and not a
    /// parse of the allele body. A regex bounded by the first `]` reads
    /// `g.[257_262TTT[3];261_264del]` as ending at the `]` that closes `TTT[3]`,
    /// finds no `;` in what it captured, and reports **0** on all 756 rows — the
    /// stale figure above, reproduced by an instrument rather than by the corpus.
    #[test]
    fn the_corpus_feeds_a_repeat_beside_a_sibling() {
        let beside_a_sibling = |s: &str| s.contains(';') && holds_repeat_notation(s);
        let rows = dump(2);
        let minted = rows
            .iter()
            .filter(|row| beside_a_sibling(&row.output))
            .count();
        let fed = rows
            .iter()
            .filter(|row| beside_a_sibling(&row.input))
            .count();
        assert!(
            minted > 0,
            "no row mints a repeat beside a sibling, so this test measures nothing"
        );
        assert!(
            fed > 0,
            "the corpus mints a repeat beside a sibling on {minted} rows and feeds one in on \
             none, so a change to where that spelling is minted measures a structural zero \
             on the authored path — add a family that feeds one in"
        );
    }

    /// Every designed tract is exactly as long as its design claims.
    ///
    /// Re-derived from the finished string rather than from the parameters that built
    /// it, because the degeneracy this guards against is invisible in the parameters:
    /// a `GCTA` flank base immediately 3' of an `ACGACG` tract *continues* it, so the
    /// real tract runs one base longer than `end` says and the two sibling placements
    /// `tract_inputs_for` builds — one inside the tract, one beside it — quietly
    /// become the same place. `tract_core` repairs both boundaries; this checks the
    /// repair worked, on all eight designs.
    ///
    /// It also pins the two [`padded`] edge rules, which the repair can reach: with a
    /// tract at offset 1 the repaired 5' boundary **is** offset 0, and a core starting
    /// with `A` extends the `ACGT` pad's own rotation into the flank.
    #[test]
    fn the_designed_tracts_are_exactly_as_long_as_they_claim() {
        let designs = tract_designs();
        assert_eq!(
            designs.len(),
            TRACT_UNITS.len() * TRACT_COPIES.len() * TRACT_OFFSETS.len(),
            "the design matrix lost a dimension"
        );
        for design in &designs {
            let core = design.core.as_bytes();
            assert_eq!(
                core.len(),
                TRACT_CORE_LEN,
                "[{}] wrong length",
                design.label
            );
            assert!(
                core.iter().all(|b| b"ACGT".contains(b)),
                "[{}] not over ACGT: {}",
                design.label,
                design.core
            );
            assert_eq!(
                &design.core[design.start..=design.end],
                design.unit.repeat(design.copies),
                "[{}] the tract is not what the design says",
                design.label
            );
            // Maximality, in both directions: the base immediately 5' of the tract
            // must not be the unit's last, and the one immediately 3' must not be
            // its first.
            let unit = design.unit.as_bytes();
            assert_ne!(
                core[design.start - 1],
                unit[unit.len() - 1],
                "[{}] the tract extends 5' of its stated start",
                design.label
            );
            assert_ne!(
                core[design.end + 1],
                unit[0],
                "[{}] the tract extends 3' of its stated end",
                design.label
            );
            // The designed tract is the only one: two consecutive copies of the unit
            // occur nowhere the design does not claim. Without this, "beside the
            // tract" could land inside a second tract the filler built by accident.
            let tandem = design.unit.repeat(2);
            for i in 0..=core.len() - tandem.len() {
                if design.core[i..i + tandem.len()] == tandem {
                    assert!(
                        i >= design.start && i + tandem.len() <= design.end + 1,
                        "[{}] a second {} tract at offset {i}",
                        design.label,
                        design.unit
                    );
                }
            }
            assert_ne!(core[0], b'A', "[{}] extends the pad 5'", design.label);
            assert_ne!(
                core[TRACT_CORE_LEN - 1],
                b'T',
                "[{}] extends the pad 3'",
                design.label
            );
        }
    }

    /// Rows this gate must not normalize, because normalizing them aborts the
    /// process under the `Test oracle` job's own flags.
    ///
    /// **These are normalizer defects that this corpus revealed, not flakes and not
    /// rows the gate is allowed to be indifferent about.** Each is filed, each is
    /// labelled `high`/`bug`, and each names the oracle it trips:
    ///
    /// - **#2036** — `FERRO_ASSERT_IDEMPOTENT`. Normalizing a long tract is not a
    ///   fixed point: every pass walks three more bases 5' and increments the copy
    ///   count (`g.558_559insACGACG` -> `g.460_558ACG[35]` -> `g.457_558ACG[36]`).
    ///   Found by `long_tract_window_provenance` on its first armed run, which is
    ///   the window-provenance defect that family was built to find.
    /// - **#2037** — `FERRO_ASSERT_IN_BOUNDS`. Combining two members at `c.1` shifts
    ///   the insertion point 5' off the coordinate space
    ///   (`c.[1dup;1T>A]` -> `c.-1_1insA`, naming position 0 of a 20-base transcript).
    ///
    /// **Why an exclusion here rather than a weakened oracle.** The oracles panic
    /// *inside* `normalize`, so a row that trips one takes the whole test binary with
    /// it — and this gate re-normalizes each row **itself**, outside the `catch_unwind`
    /// that [`normalize_through`] wraps `dump`'s own pass in. So the choice is not
    /// between checking and not checking these rows; it is between the gate answering
    /// its own question and the gate aborting on someone else's defect. Suppressing the
    /// oracle, or dropping the rows from the corpus, would instead destroy the evidence
    /// the two issues rest on.
    ///
    /// The list is **shrink-only**: every entry must still be produced by the corpus,
    /// asserted below, so a row that stops being generated — or a defect that gets
    /// fixed and re-admits its row — fails here instead of rotting into a blanket
    /// exemption.
    const ORACLE_TRIPPING_INPUTS: &[(&str, &str)] = &[
        ("NC_TEST.1:g.456_457insACGACG", "#2036"),
        ("NC_TEST.1:g.558_559insACGACG", "#2036"),
        ("NC_TEST.1:g.582_583insACGACG", "#2036"),
        ("NC_TEST.1:g.585_586insACGACG", "#2036"),
        ("NC_TEST.1:g.756_757insACGACG", "#2036"),
        ("NM_TEST.1:c.[1dup;1T>A]", "#2037"),
        ("NM_TEST.1:c.[1dup;1C>A]", "#2037"),
    ];

    /// The reference a render-time mint would rebuild, against the one the
    /// per-member pipeline was actually handed (#1946).
    ///
    /// **This is the gate on relocating the repeat mint.** `render::mint_reference_for`
    /// has to reconstruct, from a settled member and a provider, what
    /// `normalize_na_edit` was given — and `ref_seq` is not one thing: the whole
    /// spliced transcript on `c.`/`r.`/`n.`, a window on `g.`/`m.`, and a genomic
    /// window on the intronic and junction-crossing paths.
    ///
    /// # The comparison is per call site, and the first revision of it was vacuous
    ///
    /// The recorder hands back **every** `ref_seq` a normalization used: one per
    /// member of a cis allele, one per growth attempt, one per
    /// `FERRO_ASSERT_IDEMPOTENT` verification pass, and one per recursive re-entry
    /// (`normalize_na_edit` calls itself from eight arms, always with the same
    /// `ref_seq`). This gate's first revision asked whether **some** entry in that
    /// union matched the rebuilt reference — `seen.iter().any(..)` — and reported
    /// `3903 / 3903` (re-measured on this base; it shipped as `3865/3865` against an
    /// older one). That number could not have read much else: on any row with a
    /// transcript-axis member the whole-transcript entry is somewhere in the union,
    /// so `any` matched for the same reason on every row. Turning it into `all` read
    /// `3761 / 3903` — **142 members diverge, every one on the multi-exon `cx`
    /// axis** — and that is not the tightening it looks like either: different
    /// members of a legitimate multi-member allele genuinely see different windows,
    /// so `all` over the union goes red on correct behaviour.
    ///
    /// Those 142 are the whole story, and per-site selection resolves them rather
    /// than hiding them. Every one sits on a row that records **both** a `Cds`
    /// call handed the 20-base spliced transcript **and** a `BoundarySpanning` call
    /// handed a ~200-base genomic window — `NM_TESTX.1:c.6_11inv` records
    /// `[Cds@d0=20, BoundarySpanning@d0=206, Cds@d0=20, BoundarySpanning@d0=206]`,
    /// because `normalize_cds` tries the spliced frame and then re-normalizes across
    /// the exon junction in genomic coordinates. The mint's whole-transcript rebuild
    /// agrees with the `Cds` call byte-for-byte; the 206 belongs to a site the mint
    /// does not model at all, and is counted below rather than matched.
    ///
    /// Neither is the right question. The right question is per **call site**, which
    /// is why [`NaEditReference`] carries one: for each of the two shapes
    /// `mint_reference_for` claims to model, compare the set of references **it**
    /// builds against the set the pipeline was handed **at the sites that shape
    /// models**, at depth 0. Set equality in both directions, so an extra reference
    /// on either side fails — which is what makes a mis-attributed entry (a
    /// junction-crossing window filed under `Cds`) fail rather than pass.
    ///
    /// # What is asserted
    ///
    /// 1. **Every recorded call is attributed.** A depth-0 call at
    ///    [`NaEditCallSite::Unattributed`] means a normalizer entry point reaches
    ///    `normalize_na_edit` without a site guard, and every conclusion below is
    ///    then about a subset nobody delimited. Asserted zero.
    /// 2. **The recorder did not truncate.** A capped buffer that silently drops
    ///    entries looks exactly like a set that never contained them.
    /// 3. **The whole-transcript sets are equal.** The reference is served entire on
    ///    `c.`/`r.`/`n.`, so there is no window to choose and no reason for a rebuilt
    ///    reference to differ by a byte from the one the pipeline indexed — nor for
    ///    the pipeline to have used one the mint does not build. That is a property,
    ///    statable without any measured number. Asserted over the rows where it is
    ///    answerable, with the answerable-row count asserted non-zero so `0 of 0`
    ///    cannot pass as a result.
    /// 4. **One whole transcript per row.** The set is asserted to be a *singleton*
    ///    wherever it is non-empty, which is the "nothing to choose" claim stated
    ///    directly rather than assumed.
    ///
    /// # What is reported, and why it is not asserted
    ///
    /// **The genomic frame.** It cannot be byte-equal by construction: the pipeline's
    /// window is centred on the *input* position and a render stage only has the
    /// settled one, and the two do not even share a left edge —
    /// `normalize_in_grown_window` fetches from `start - w` where `mint_reference_for`
    /// fetches from `start - 1 - w`. So the question there is containment, and
    /// pinning today's ratio would make the count *be* the property.
    ///
    /// **The frames the mint does not model.** A `c.` member whose span crosses an
    /// exon/exon junction, or which sits in an intron, is normalized against a
    /// *genomic* window by `normalize_boundary_spanning_cds` / `normalize_intronic_cds`
    /// — `mint_reference_for` returns the spliced transcript for it, because that is
    /// what its `CisKind` says. Those rows are counted under
    /// `rows_using_an_unmodelled_frame` rather than being folded into a match, and a
    /// row whose transcript-side calls were **only** at such a site is excluded from
    /// assertion 3 (there is no whole-transcript call to compare against) and counted
    /// under `rows_with_no_modelled_call`. That is the next piece of work, and naming
    /// it as a population is the difference between an open gap and a hidden one.
    ///
    /// # Measured over `dump(1)`
    ///
    /// Read the run's own `MINT-REFERENCE` line rather than this table; it is here so
    /// a reader knows the order of magnitude and which buckets exist.
    ///
    /// | population | result |
    /// |---|---|
    /// | rows where the whole-transcript sets are compared | see the printed line — **asserted equal** |
    /// | genomic-frame members at least as wide as the widest witnessed window | reported |
    /// | members with no rebuildable reference | reported |
    /// | rows using a frame `mint_reference_for` does not model | reported |
    /// | held out (#2036, #2037) | 7 rows |
    ///
    /// **Every asserted figure is identical with the three `FERRO_ASSERT_*` flags set
    /// and unset, measured both ways** — same 2209 rows compared, same zero
    /// mismatches, same `1586 / 2344`, same zero unattributed calls. One *reported*
    /// counter does move: `rows_using_an_unmodelled_frame` reads **102** unarmed and
    /// **103** armed, because `FERRO_ASSERT_IDEMPOTENT`'s verification pass is a
    /// second normalization whose own depth-0 calls are recorded too, and on one row
    /// it reaches a junction-crossing site the first pass did not. That is a fact
    /// about the corpus rather than about the mint, and it is stated rather than
    /// smoothed over: a counter over a *union* is exactly what stays flag-sensitive,
    /// which is why the assertions are over per-site sets instead.
    ///
    /// That distinction is worth stating because a still-earlier revision did
    /// **not** have it: it compared against the *narrowest* reference in the union,
    /// and `FERRO_ASSERT_IDEMPOTENT`'s verification pass contributes narrower entries,
    /// so its genomic figure read `2344/2344` in `Test oracle` and `2274/2344` in
    /// `Test` (measured against that revision's own base — do not carry those two
    /// numbers onto this one) — a gate whose number improves when you arm the oracles
    /// is measuring the oracles. Selecting by call site and depth removes that
    /// dependence at the
    /// source: the verification pass's calls are depth-0 calls of a second
    /// normalization, and they land in the same sets as the first pass's, which for
    /// the transcript frame is the *same* whole transcript and for the genomic frame
    /// is folded into the same maximum.
    ///
    /// Seven rows are held out; see [`ORACLE_TRIPPING_INPUTS`], which names the issue
    /// behind each.
    #[test]
    fn the_render_time_reference_matches_what_the_pipeline_was_given() {
        use ferro_hgvs::normalize::{
            mint_reference_for, na_edit_references_overflowed, reference_digest,
            take_na_edit_references, MintFrame, NaEditCallSite, NaEditReferenceRecording,
            MINT_WINDOW,
        };
        use std::collections::BTreeSet;

        /// The sites at which the pipeline is handed a whole transcript. These are
        /// exactly the sites [`MintFrame::WholeTranscript`] models.
        const WHOLE_TRANSCRIPT_SITES: &[NaEditCallSite] =
            &[NaEditCallSite::Cds, NaEditCallSite::Tx, NaEditCallSite::Rna];
        /// The sites at which the pipeline is handed a `g.`/`m.` window. These are
        /// exactly the sites [`MintFrame::GenomicWindow`] models.
        const GENOMIC_WINDOW_SITES: &[NaEditCallSite] =
            &[NaEditCallSite::GrownWindow, NaEditCallSite::UnshiftedWindow];
        /// Sites the mint does not model at all: a genomic window reached from a
        /// transcript-axis description.
        const UNMODELLED_SITES: &[NaEditCallSite] =
            &[NaEditCallSite::Intronic, NaEditCallSite::BoundarySpanning];

        let mut tx_rows_compared = 0usize;
        let mut tx_row_mismatches: Vec<String> = Vec::new();
        let mut rows_with_no_modelled_call = 0usize;
        let mut rows_using_an_unmodelled_frame = 0usize;
        let mut rows_with_a_declined_member = 0usize;
        let mut unattributed_calls = 0usize;
        let (mut g_total, mut g_contained, mut g_unwitnessed) = (0usize, 0usize, 0usize);
        let mut no_reference = 0usize;
        let mut held_out_seen: BTreeSet<&str> = BTreeSet::new();

        // Built BEFORE arming: `dump` normalizes the whole corpus itself, and those
        // normalizations are not the ones this gate is asking about. Arming first
        // fills the buffer with them, which the cap assertion below catches.
        let rows = dump(1);

        // Armed once for the whole walk. Recording is off by default, so nothing
        // outside this test pays for it or accumulates into it.
        let _recording = NaEditReferenceRecording::arm();

        for row in rows {
            if row.axis == "p" {
                continue;
            }
            // Held out BEFORE normalizing: the oracles abort inside `normalize`, and
            // this gate's own call is not wrapped in `catch_unwind`.
            if let Some((held, _)) = ORACLE_TRIPPING_INPUTS
                .iter()
                .find(|(input, _)| *input == row.input)
            {
                held_out_seen.insert(*held);
                continue;
            }
            let provider = provider_for(row.axis, &row.core);
            let direction = match row.direction {
                "5prime" => ShuffleDirection::FivePrime,
                _ => ShuffleDirection::ThreePrime,
            };
            let normalizer = Normalizer::with_config(
                provider_for(row.axis, &row.core),
                NormalizeConfig::default().with_direction(direction),
            );
            let Ok(input) = parse_hgvs(&row.input) else {
                continue;
            };
            let _ = take_na_edit_references();
            let normalized = normalizer.normalize(&input);
            // Drained unconditionally, and *before* the `else { continue }`. A
            // refused normalization still records the references it was handed
            // before it refused, and leaving them in the buffer bleeds them into
            // the next row's set — which the cap assertion below caught, and which
            // the union-based revision of this gate could not have.
            let seen = take_na_edit_references();
            let Ok(output) = normalized else {
                continue;
            };
            if seen.is_empty() {
                continue;
            }

            // Only depth-0 calls: the eight recursive arms re-enter with the same
            // `ref_seq`, so counting them would weight one entry point by how many
            // rewrites it happened to take.
            let entry_calls = || seen.iter().filter(|e| e.depth == 0);
            unattributed_calls += entry_calls()
                .filter(|e| e.site == NaEditCallSite::Unattributed)
                .count();
            let witnessed_tx: BTreeSet<(usize, u64)> = entry_calls()
                .filter(|e| WHOLE_TRANSCRIPT_SITES.contains(&e.site))
                .map(|e| (e.len, e.digest))
                .collect();
            let widest_genomic = entry_calls()
                .filter(|e| GENOMIC_WINDOW_SITES.contains(&e.site))
                .map(|e| e.len)
                .max();
            let unmodelled: BTreeSet<NaEditCallSite> = entry_calls()
                .filter(|e| UNMODELLED_SITES.contains(&e.site))
                .map(|e| e.site)
                .collect();
            if !unmodelled.is_empty() {
                rows_using_an_unmodelled_frame += 1;
            }

            let members: Vec<_> = match &output {
                ferro_hgvs::hgvs::variant::HgvsVariant::Allele(a) => a.variants.clone(),
                other => vec![other.clone()],
            };
            let mut built_tx: BTreeSet<(usize, u64)> = BTreeSet::new();
            let mut declined_here = false;
            for member in &members {
                let Some(built) = mint_reference_for(member, &provider, MINT_WINDOW) else {
                    no_reference += 1;
                    declined_here = true;
                    continue;
                };
                let digest = reference_digest(&built.bases);
                match built.frame {
                    MintFrame::WholeTranscript => {
                        built_tx.insert((built.bases.len(), digest));
                    }
                    MintFrame::GenomicWindow => {
                        g_total += 1;
                        // The pipeline's window is centred elsewhere and starts one
                        // base 3' of this one, so identity is not the question;
                        // whether the rebuilt window is at least as wide as the
                        // **widest** the pipeline used at a site this frame models
                        // is. Widest, because `insertion_to_repeat` walks outward
                        // bounded by `ref_seq.len()` — a mint is only safe if it is
                        // at least as wide as every window the pipeline used.
                        match widest_genomic {
                            Some(widest) if built.bases.len() >= widest => g_contained += 1,
                            Some(_) => {}
                            // The mint built a genomic window for a member the
                            // pipeline never handed one to. Counted, not passed.
                            None => g_unwitnessed += 1,
                        }
                    }
                }
            }
            if declined_here {
                rows_with_a_declined_member += 1;
            }

            // The whole-transcript comparison is answerable only when the pipeline
            // made a whole-transcript call and every member was minted. A row whose
            // transcript-side work happened entirely in a frame the mint does not
            // model has nothing to compare against, and saying so is the point.
            if witnessed_tx.is_empty() && !built_tx.is_empty() {
                rows_with_no_modelled_call += 1;
            } else if !witnessed_tx.is_empty() && !declined_here {
                tx_rows_compared += 1;
                if built_tx != witnessed_tx && tx_row_mismatches.len() < 12 {
                    tx_row_mismatches.push(format!(
                        "axis={} input={} built={:?} witnessed={:?} unmodelled_sites={:?}",
                        row.axis,
                        row.input,
                        built_tx.iter().map(|(l, _)| *l).collect::<Vec<_>>(),
                        witnessed_tx.iter().map(|(l, _)| *l).collect::<Vec<_>>(),
                        unmodelled,
                    ));
                } else if built_tx != witnessed_tx {
                    tx_row_mismatches.push(String::from("<further mismatches elided>"));
                }
                assert!(
                    witnessed_tx.len() == 1,
                    "the reference is served whole on `c.`/`r.`/`n.`, so one \
                     normalization cannot have been handed two different \
                     whole-transcript references: {} saw {:?}",
                    row.input,
                    witnessed_tx.iter().map(|(l, _)| *l).collect::<Vec<_>>(),
                );
            }
        }

        eprintln!(
            "MINT-REFERENCE whole-transcript sets equal on {tx_rows_compared} rows \
             ({} mismatching); genomic {g_contained}/{g_total} at least as wide \
             ({g_unwitnessed} unwitnessed); {no_reference} members had no rebuildable \
             reference across {rows_with_a_declined_member} rows; \
             {rows_using_an_unmodelled_frame} rows used a frame the mint does not model, \
             {rows_with_no_modelled_call} of them with no modelled call at all; \
             {unattributed_calls} unattributed calls; {} rows held out",
            tx_row_mismatches.len(),
            ORACLE_TRIPPING_INPUTS.len()
        );

        // Shrink-only: a held-out row the corpus no longer builds must be removed from
        // the list deliberately, not left as a silent exemption.
        let expected: BTreeSet<&str> = ORACLE_TRIPPING_INPUTS.iter().map(|(i, _)| *i).collect();
        assert_eq!(
            held_out_seen,
            expected,
            "the held-out list has drifted from the corpus; entries never produced: {:?}",
            expected.difference(&held_out_seen).collect::<Vec<_>>()
        );

        assert!(
            !na_edit_references_overflowed(),
            "the reference recorder hit its cap, so every set above is a truncation \
             of the real one and a missing entry is indistinguishable from a mismatch"
        );
        assert_eq!(
            unattributed_calls, 0,
            "a normalizer entry point reached `normalize_na_edit` without a \
             `NaEditSiteGuard`, so the site selection below is over a population \
             nobody delimited"
        );
        assert!(
            tx_rows_compared > 0,
            "no row reached the whole-transcript comparison, so the equality \
             assertion below is vacuous"
        );
        assert!(
            tx_row_mismatches.is_empty(),
            "a render-time mint must rebuild exactly the whole-transcript references \
             the pipeline was handed at the sites this frame models: the reference is \
             served whole on `c.`/`r.`/`n.`, so there is no window to choose and \
             nothing to disagree about. {} of {tx_rows_compared} rows disagree:\n{}",
            tx_row_mismatches.len(),
            tx_row_mismatches.join("\n"),
        );
    }

    /// The long tracts bracket the normalizer's own window, on both sides.
    ///
    /// **This is the property the family exists for, so it is asserted rather than
    /// assumed.** A set of tracts that all sat on one side of the window would produce
    /// a populated-looking family that could not answer the question — the same shape
    /// as every corpus blindness this file records, one level down.
    ///
    /// `window_size` is read from `NormalizeConfig::default()` rather than restated,
    /// per #1925: a guard that hardcodes the value it guards cannot observe that value
    /// changing. The *other* constant in play, `merge::CANONICAL_PAD` (128), is
    /// `pub(crate)` and an example is a separate crate, so it cannot be imported at
    /// all. Rather than restate `128` — which is exactly the antipattern #1925 names —
    /// the upper bound is expressed as **twice** `window_size`, which exceeds 128 for
    /// any window at or above 64 and so brackets both constants without naming the
    /// unreachable one.
    #[test]
    fn the_long_tracts_bracket_the_normalizer_window() {
        let window = NormalizeConfig::default().window_size as usize;
        let lengths: Vec<usize> = long_tract_designs()
            .iter()
            .map(|d| d.end - d.start + 1)
            .collect();
        assert_eq!(lengths.len(), LONG_TRACT_COPIES.len());
        assert!(
            lengths.iter().any(|&l| l < window),
            "no tract is shorter than the {window}-base window, so the family cannot \
             show a walk that stays inside one: {lengths:?}"
        );
        assert!(
            lengths.iter().any(|&l| l > 2 * window),
            "no tract exceeds twice the {window}-base window, so the family cannot \
             show a walk that leaves one: {lengths:?}"
        );
        // Both sides of the pad as well, without naming it: a tract between the window
        // and twice it, and one beyond.
        assert!(
            lengths.iter().any(|&l| l > window && l < 2 * window),
            "no tract lands between the window and twice it: {lengths:?}"
        );
    }

    /// The long-tract designs are maximal and well formed, like their short siblings.
    ///
    /// Separate from `the_designed_tracts_are_exactly_as_long_as_they_claim` because
    /// that one pins the eight short designs' own count and offsets; this checks the
    /// same *invariants* on a different design set. The boundary repair in
    /// `tract_core` is shared, so a regression in it must fail on both.
    #[test]
    fn the_long_tract_designs_are_maximal() {
        let designs = long_tract_designs();
        assert_eq!(designs.len(), LONG_TRACT_COPIES.len());
        for design in &designs {
            let core = design.core.as_bytes();
            assert_eq!(core.len(), LONG_TRACT_CORE_LEN, "[{}]", design.label);
            assert!(
                core.iter().all(|b| b"ACGT".contains(b)),
                "[{}] not over ACGT",
                design.label
            );
            assert_eq!(
                &design.core[design.start..=design.end],
                design.unit.repeat(design.copies),
                "[{}] the tract is not what the design says",
                design.label
            );
            let unit = design.unit.as_bytes();
            assert_ne!(
                core[design.start - 1],
                unit[unit.len() - 1],
                "[{}] the tract extends 5' of its stated start",
                design.label
            );
            assert_ne!(
                core[design.end + 1],
                unit[0],
                "[{}] the tract extends 3' of its stated end",
                design.label
            );
            assert_ne!(core[0], b'A', "[{}] extends the pad 5'", design.label);
            assert_ne!(
                core[LONG_TRACT_CORE_LEN - 1],
                b'T',
                "[{}] extends the pad 3'",
                design.label
            );
        }
    }

    /// The long-tract family emits both junctions on both DNA axes it runs on.
    ///
    /// `cx` is deliberately absent — asserted, so that wiring it in later has to be a
    /// decision rather than an accident that would draw a 20-base exon table over an
    /// 800-base core.
    #[test]
    fn the_long_tract_family_emits_both_junctions() {
        for design in long_tract_designs() {
            for axis in ["g", "c"] {
                let inputs = long_tract_inputs_for(axis, &design);
                assert_eq!(inputs.len(), 2, "[{} {axis}]", design.label);
                let payload = design.unit.repeat(2);
                assert!(
                    inputs
                        .iter()
                        .all(|i| i.ends_with(&payload) && !i.contains(';')),
                    "[{} {axis}] both shapes are lone insertions of TWO units — one \
                     unit can never mint a repeat: {inputs:?}",
                    design.label
                );
                assert_ne!(inputs[0], inputs[1], "[{} {axis}]", design.label);
            }
        }
        let cx_rows = dump(0)
            .into_iter()
            .filter(|r| r.family == LONG_TRACT_FAMILY.0 && r.axis == "cx")
            .count();
        assert_eq!(
            cx_rows, 0,
            "the long tracts must not be drawn on the cx axis"
        );
    }

    /// The tract family emits all six shapes, on every axis, for every design.
    ///
    /// The sibling of `every_family_emits_rows_on_both_axes`, and separate for the
    /// same reason `every_protein_family_emits_rows` is: this family is not crossed
    /// with the random cores, so it is not in [`FAMILIES`] and that loop cannot see it.
    ///
    /// The count is asserted rather than mere non-emptiness because the six shapes
    /// are the whole content of the family — each names one pass of the tax #1946 is
    /// about — and a shape lost to a bad offset would leave a family that still looks
    /// populated.
    #[test]
    fn the_tract_family_emits_every_shape_on_every_axis() {
        for design in tract_designs() {
            for axis in ["g", "c", "cx"] {
                let inputs = tract_inputs_for(axis, &design);
                assert_eq!(
                    inputs.len(),
                    6,
                    "[{} {axis}] expected six shapes, got {inputs:?}",
                    design.label
                );
                assert!(
                    inputs.iter().all(|i| i.contains(';')),
                    "[{} {axis}] every shape pairs a member with a sibling",
                    design.label
                );
                // Two authored repeats, one authored dup, and three shapes whose
                // presentational form ferro has to choose for itself.
                assert_eq!(
                    inputs.iter().filter(|i| holds_repeat_notation(i)).count(),
                    3,
                    "[{} {axis}] expected three authored repeats",
                    design.label
                );
                assert_eq!(
                    inputs.iter().filter(|i| i.contains("dup")).count(),
                    1,
                    "[{} {axis}] expected one authored duplication",
                    design.label
                );
            }
        }
    }

    /// Every input the tract family emits parses.
    ///
    /// The family's two authored-repeat shapes are the reason this is worth its own
    /// test rather than being left to the dump: `repeated.md:21` admits a repeat on a
    /// `c.` reference only for units whose length is a multiple of three, so a unit
    /// of any other length would make a third of this family spec-invalid on two of
    /// its three axes — and an unparseable input is recorded as a row, not as a
    /// failure, so the corpus would keep looking populated while measuring nothing.
    #[test]
    fn no_tract_input_is_rejected() {
        for design in tract_designs() {
            for axis in ["g", "c", "cx"] {
                for input in tract_inputs_for(axis, &design) {
                    assert!(
                        parse_hgvs(&input).is_ok(),
                        "[{}] {input} does not parse",
                        design.label
                    );
                }
            }
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

    /// `repeat_beside_a_sibling` varies the property `repeat_footprint` keys on:
    /// the repeat must both **grow** and **shrink**.
    ///
    /// # Why this needs a test rather than a comment
    ///
    /// The family was added (#1749) to close a corpus blindness, and shipped
    /// closing half of it: every row was `{unit}[3]` over a 6-base tract with a
    /// 3-base unit, so `3 x 3 = 9 >= 6` and the repeat could only ever grow.
    /// `repeat_footprint`'s **shrinking** branch — `spanning(removed_from, end,
    /// false)`, the one the change is proudest of — was structurally unreachable
    /// from the corpus, so any `0 moved` over it was a claim about the corpus.
    /// That is exactly the failure the family exists to prevent, one level in.
    ///
    /// The check is arithmetic against the tract, not a count of spellings: it
    /// asserts `unit.len() * count` lands on **both** sides of the tract length,
    /// so widening the tract or the unit without re-picking the counts fails
    /// here rather than silently collapsing to one direction again.
    #[test]
    fn the_repeat_sibling_family_varies_the_repeat_direction() {
        let mut covered = 0;
        for core in corpus_sequences(24) {
            let rows = inputs_for("repeat_beside_a_sibling", "g", &core);
            let Some(t) = first_tandem_triplet(&core) else {
                assert!(rows.is_empty(), "emitted a repeat sibling row for {core}");
                continue;
            };
            covered += 1;
            let unit = &core[t..t + 3];
            // The tract the family anchors on is `t ..= t + 5`, six bases.
            const TRACT: usize = 6;
            let mut grows = 0;
            let mut shrinks = 0;
            for row in &rows {
                let count: usize = row
                    .split(&format!("{unit}["))
                    .nth(1)
                    .and_then(|rest| rest.split(']').next())
                    .expect("every row spells the repeat as `<unit>[<count>]`")
                    .parse()
                    .expect("the count is an integer");
                if unit.len() * count >= TRACT {
                    grows += 1;
                } else {
                    shrinks += 1;
                }
            }
            assert!(
                grows > 0 && shrinks > 0,
                "{core}: the family emitted {grows} growing and {shrinks} shrinking \
                 repeats — one direction of `repeat_footprint` is unreachable from \
                 the corpus, so a zero over it would be structural"
            );
        }
        assert_eq!(
            covered, 21,
            "tandem-triplet coverage moved; the generator or seed count changed"
        );
    }

    /// `authored_repeat_beside_a_sibling` accounts for the cores it drops.
    ///
    /// It needs a tandem triplet **and** two bases 3' of it to put the sibling on, so
    /// it covers one fewer core than `repeat_expansion`. Naming both numbers is the
    /// point: a family that silently covered fewer cores than its neighbour would
    /// look like the same population while measuring a smaller one, and the
    /// difference — a single core whose tract starts at offset 14 — is exactly the
    /// kind of drop the generator doctrine says must be counted rather than absorbed.
    #[test]
    fn the_authored_repeat_family_covers_the_cores_with_room() {
        let (mut with_tract, mut covered) = (0, 0);
        for core in corpus_sequences(24) {
            let rows = inputs_for("authored_repeat_beside_a_sibling", "g", &core);
            let Some(t) = first_tandem_triplet(&core) else {
                assert!(
                    rows.is_empty(),
                    "emitted for {core}, which has no tandem triplet"
                );
                continue;
            };
            with_tract += 1;
            if t + 8 >= core.len() {
                assert!(
                    rows.is_empty(),
                    "emitted for {core}, whose tract at {t} leaves no room for a sibling"
                );
                continue;
            }
            covered += 1;
            assert_eq!(rows.len(), 2, "expected both counts for {core}");
            for row in &rows {
                assert!(
                    row.contains(';'),
                    "the whole point is that the repeat keeps a sibling: {row}"
                );
            }
        }
        assert_eq!(with_tract, 21, "tandem-triplet coverage moved");
        assert_eq!(
            covered, 20,
            "coverage moved; exactly one core has a tract too close to the 3' end"
        );
    }

    /// `run_beside_a_distant_member` accounts for the cores it drops.
    ///
    /// The arm emits a pair of rows only for a core carrying a clean four-base left
    /// rotation with room for a distant member — a window `core[a..a+4]` with no
    /// cyclically adjacent pair equal, `a` in `0..len-6`. It emits **nothing** for a
    /// core with no such window (homopolymer-heavy cores have none). Nothing else
    /// records that drop, so without this guard a generator or seed change that
    /// shrank the covered set would leave the family measuring a smaller population
    /// while its `0 moved` still read as a fixed-point measurement — the exact
    /// structural-zero trap the generator doctrine names. Both numbers are pinned so
    /// a shift in either direction fails here rather than passing silently.
    #[test]
    fn the_distant_member_family_accounts_for_the_cores_it_drops() {
        let has_clean_rotation = |core: &str| {
            let b = core.as_bytes();
            (0..b.len().saturating_sub(6)).any(|a| (0..4).all(|j| b[a + j] != b[a + (j + 1) % 4]))
        };
        let (mut emitted, mut dropped) = (0, 0);
        for core in corpus_sequences(24) {
            let rows = inputs_for("run_beside_a_distant_member", "g", &core);
            if has_clean_rotation(&core) {
                emitted += 1;
                assert_eq!(
                    rows.len(),
                    2,
                    "{core} has a clean rotation window but did not emit both rows"
                );
                for row in &rows {
                    // Each row carries the distant member (`;` … `>`) beside the run,
                    // which is the whole geometry the family exists to measure.
                    assert!(
                        row.contains(';') && row.contains('>'),
                        "the distant member is what makes this family the #2192 shape: {row}"
                    );
                }
            } else {
                dropped += 1;
                assert!(
                    rows.is_empty(),
                    "{core} has no clean rotation window yet emitted {rows:?}"
                );
            }
        }
        assert_eq!(
            emitted, 43,
            "clean-rotation coverage moved; the generator or seed count changed"
        );
        assert_eq!(
            dropped, 5,
            "the dropped-core count moved; a core gained or lost a clean rotation window"
        );
        assert_eq!(
            emitted + dropped,
            corpus_sequences(24).len(),
            "every core is either covered or accounted as dropped"
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
    // eighteen of the twenty-two families it cannot — they are sets of descriptions,
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

    /// English for the family counts this file quotes.
    ///
    /// **Composed rather than tabulated, and that is the point.** This used to be a
    /// flat table ending at `twenty`, so the run that took `FAMILIES` past twenty
    /// families died with `no word for 21; extend WORDS in number_word` — a panic
    /// naming a helper, from a test whose job is to name *stale prose sites*. The
    /// table's end was a second, invisible boundary sitting beside the real one, and
    /// extending it by four entries would only have moved the cliff to twenty-five.
    ///
    /// Composing tens and units removes the boundary for any count this corpus can
    /// plausibly reach, so the only thing that can now fail here is the thing the
    /// caller is actually testing.
    fn number_word(n: usize) -> String {
        const UNITS: &[&str] = &[
            "zero",
            "one",
            "two",
            "three",
            "four",
            "five",
            "six",
            "seven",
            "eight",
            "nine",
            "ten",
            "eleven",
            "twelve",
            "thirteen",
            "fourteen",
            "fifteen",
            "sixteen",
            "seventeen",
            "eighteen",
            "nineteen",
        ];
        const TENS: &[&str] = &[
            "", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety",
        ];
        match n {
            0..=19 => UNITS[n].to_string(),
            20..=99 => match (n / 10, n % 10) {
                (t, 0) => TENS[t].to_string(),
                (t, u) => format!("{}-{}", TENS[t], UNITS[u]),
            },
            _ => panic!("no word for {n}; number_word covers 0..=99"),
        }
    }

    /// Strip Markdown and punctuation so a word inside `**bold**` or before a
    /// comma still compares equal.
    fn bare(word: &str) -> String {
        word.trim_matches(|c: char| !c.is_alphanumeric() && c != '-')
            .to_ascii_lowercase()
    }

    /// This file's prose must state the family counts this file actually has.
    ///
    /// # Why this is a guard rather than a correction
    ///
    /// The header is where a reader decides whether their change is measurable,
    /// and it drifted (#1947): seven sites gave the sequence-keyed count as 13
    /// (five of them pairing it with a DNA total of 15), one gave it as 9, while
    /// `FAMILIES` held 14 and the DNA total was 16. Every one of those numbers
    /// understated the corpus, which is the direction that makes a reader trust a
    /// measurement more than they should.
    ///
    /// The stale phrasings are given as digits above **on purpose**: the scan
    /// below reads comment text, so spelling them out here would make this
    /// comment its own first failure. `claude_md_adjudication_tables` records the
    /// same hazard and the same answer — keep the prose clear of the pattern
    /// rather than exempting the file, which would blind the scan to it for good.
    ///
    /// Correcting the seven sites fixes today and guarantees a repeat: the counts
    /// live in a `const` and the prose restates them, which is precisely the
    /// "a count restated instead of imported" shape `CLAUDE.md` names. So the
    /// counts are **derived here** — from `FAMILIES.len()`, from
    /// `PROTEIN_FAMILIES.len()`, and from the singleton `const … _FAMILY`
    /// declarations counted out of the source — and the prose is asserted against
    /// them. Adding a family now fails this test with the sites to fix, instead of
    /// silently making the header wrong again.
    ///
    /// Checked against a deliberate sabotage rather than assumed: appending one
    /// entry to `FAMILIES` turns this red and names every stale site.
    #[test]
    fn the_header_states_the_family_counts_this_file_actually_has() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("dump_normalized_corpus.rs");
        let file = fs::read_to_string(&path).expect("read this generator's own source");

        // Comment lines only. The claims being checked are prose, and the scan
        // would otherwise read this test's own string literals as stale sites —
        // it did, on the first draft.
        let src: String = file
            .lines()
            .filter(|line| line.trim_start().starts_with("//"))
            .collect::<Vec<_>>()
            .join("\n");

        // Two families bring their own designed reference and so are keyed by
        // label rather than by sequence. Counted from the source rather than
        // written down, so adding a third is not silently missed.
        let designed_reference = file
            .lines()
            .filter(|line| line.starts_with("const ") && line.contains("_FAMILY: (&str, &str)"))
            .count();
        assert!(
            designed_reference >= 2,
            "expected at least LONG_FAMILY and REVCOMP_FAMILY; found {designed_reference}"
        );

        let sequence_keyed = number_word(FAMILIES.len());
        let dna_total = number_word(FAMILIES.len() + designed_reference);
        let mut wrong: Vec<String> = Vec::new();

        // How many sites each anchor actually reached. Every anchor below is a
        // `match_indices` loop, so a phrasing that stops appearing does not fail
        // — the loop body simply never runs, `wrong` stays empty, and the test
        // passes having verified nothing. Two of the three hang on a **single
        // sentence each**, and these are 80-column doc comments, so a rewrap is
        // enough: anchor 1 reads across whitespace from `" of the "`, and a
        // rewrap that ends a line there makes the next token the `//!` marker,
        // so the site is skipped in silence.
        //
        // That is the exact shape this whole test is about — a check that keys
        // on a property it never confirms is present — so the counts are
        // asserted below rather than trusted.
        let (mut n1, mut n2, mut n3) = (0usize, 0usize, 0usize);

        // First anchor: the two-number form, pairing the sequence-keyed count
        // with the DNA total.
        for (at, _) in src.match_indices(" of the ") {
            let mut after = src[at + " of the ".len()..].split_whitespace();
            let (Some(total), Some(noun)) = (after.next(), after.next()) else {
                continue;
            };
            if bare(noun) != "families" {
                continue;
            }
            // Counted here rather than at the top of the loop: `" of the "`
            // occurs throughout this file in ordinary prose, and only the sites
            // that survive the `families` test are ones this anchor judges.
            n1 += 1;
            let part = src[..at].split_whitespace().next_back().unwrap_or("");
            if bare(part) != sequence_keyed || bare(total) != dna_total {
                wrong.push(format!(
                    "  \"{} of the {} families\" — expected \"{sequence_keyed} of the \
                     {dna_total} families\"",
                    bare(part),
                    bare(total)
                ));
            }
        }

        // Second anchor: the form in the LONG_FAMILY note, which counts
        // `FAMILIES` alone. Phrased here without quoting the pattern, since this
        // scan reads comment text.
        for (at, _) in src.match_indices(" families there") {
            n2 += 1;
            let word = src[..at].split_whitespace().next_back().unwrap_or("");
            if bare(word) != sequence_keyed {
                wrong.push(format!(
                    "  \"{} families there\" — expected \"{sequence_keyed} families there\"",
                    bare(word)
                ));
            }
        }

        // Third anchor: the form that qualifies the count as sequence-keyed.
        for (at, _) in src.match_indices(" sequence-keyed families") {
            n3 += 1;
            let word = src[..at].split_whitespace().next_back().unwrap_or("");
            if bare(word) != sequence_keyed {
                wrong.push(format!(
                    "  \"{} sequence-keyed families\" — expected \"{sequence_keyed} \
                     sequence-keyed families\"",
                    bare(word)
                ));
            }
        }

        // Asserted BEFORE `wrong`, because an empty `wrong` means two opposite
        // things and only this tells them apart: every site agreed, or there
        // were no sites. Non-zero is the minimum that closes the silent-disarm
        // hole; pinning the exact counts (5/1/1) would additionally catch a new
        // unguarded restatement being *added*, and is deliberately not done
        // here — it would fail on every legitimate edit that adds a sentence.
        assert!(
            n1 > 0 && n2 > 0 && n3 > 0,
            "the count-bearing phrasings this guard anchors on are gone \
             ({n1}/{n2}/{n3} sites) — the prose was reworded and this scan now checks nothing"
        );

        assert!(
            wrong.is_empty(),
            "this file's prose states family counts it does not have.\n\
             derived: FAMILIES = {} ({sequence_keyed}), designed-reference singletons = {}, \
             DNA total = {} ({dna_total})\n\
             anchor sites reached: {n1}/{n2}/{n3}\n\
             stale sites:\n{}",
            FAMILIES.len(),
            designed_reference,
            FAMILIES.len() + designed_reference,
            wrong.join("\n")
        );
    }

    /// `PROTEIN_FAMILIES` is deliberately **not** anchored, and this records why.
    ///
    /// The guard above derives three quantities from `const`s and compares each
    /// against the prose. `PROTEIN_FAMILIES.len()` used to be derived beside them
    /// and then interpolated into the failure message only — never compared —
    /// which reads as a fourth guarded count and is not one. Adding a sixth
    /// protein family would have moved nothing that test can see.
    ///
    /// The reason is that the header states no protein-family count to anchor:
    /// `PROTEIN_FAMILIES` is described qualitatively (which axes it covers, why
    /// it is kept out of `FAMILIES`) and never enumerated in words, so there is
    /// no phrasing to key on and inventing one to guard would be writing the
    /// claim in order to check it.
    ///
    /// So this test asserts the *absence* rather than leaving the gap implicit.
    /// If someone later writes a protein count into the header, this fails and
    /// the answer is to add a fourth anchor beside the other three — not to
    /// delete this.
    #[test]
    fn the_header_states_no_protein_family_count_to_anchor() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("dump_normalized_corpus.rs");
        let file = fs::read_to_string(&path).expect("read this generator's own source");
        let src: String = file
            .lines()
            .filter(|line| line.trim_start().starts_with("//"))
            .collect::<Vec<_>>()
            .join("\n");

        let stated = number_word(PROTEIN_FAMILIES.len());
        let phrasings = [
            format!("{stated} protein families"),
            format!("{stated} protein shape families"),
        ];
        for phrasing in &phrasings {
            assert!(
                !src.contains(phrasing.as_str()),
                "the header now states a protein-family count (\"{phrasing}\"). \
                 `the_header_states_the_family_counts_this_file_actually_has` does not \
                 anchor one, so it would go stale unwatched — add a fourth anchor there \
                 rather than deleting this test."
            );
        }
    }
}
