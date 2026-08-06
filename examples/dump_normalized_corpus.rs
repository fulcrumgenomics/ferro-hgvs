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
}

fn main() -> ExitCode {
    let cli = Cli::parse();
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

/// One measured row. `was_fixed_point` records whether the *input* was already its
/// own normalized form on the revision that produced this dump — which is what
/// makes the cheap/expensive split derivable from the baseline dump alone, without
/// re-running the old code.
struct Row {
    reference: String,
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
    let acc = if axis == "g" {
        GENOMIC_CONTIG
    } else {
        TX_ACCESSION
    };
    let prefix = format!("{acc}:{axis}.");
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
            for input in long_inputs_for(axis, &core) {
                let output = normalize_one(axis, &core, &input, direction);
                rows.push(Row {
                    reference: label.clone(),
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
    for core in corpus_sequences(seeds) {
        for (axis, direction, dir_label) in axes_and_directions() {
            for (family, _) in FAMILIES {
                for input in inputs_for(family, axis, &core) {
                    let output = normalize_one(axis, &core, &input, direction);
                    rows.push(Row {
                        reference: core.clone(),
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
    if p <= CDS_END {
        p.to_string()
    } else {
        format!("*{}", p - CDS_END)
    }
}

fn inputs_for(family: &str, axis: &str, core: &str) -> Vec<String> {
    let bytes = core.as_bytes();
    let acc = if axis == "g" {
        GENOMIC_CONTIG
    } else {
        TX_ACCESSION
    };
    let prefix = format!("{acc}:{axis}.");
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
            _ => unreachable!("unknown family {family}"),
        }
    }
    out
}

/// Normalize, reporting a decline or a panic as data rather than aborting the dump.
/// A row that errors on one revision and succeeds on the other is exactly the
/// cheap-vs-expensive distinction the diff needs to see, so it must not be dropped.
fn normalize_one(axis: &str, core: &str, input: &str, direction: ShuffleDirection) -> String {
    let Ok(variant) = parse_hgvs(input) else {
        return "<parse-error>".to_string();
    };
    let provider = match axis {
        "g" => genomic_provider(core),
        _ => coding_provider(core),
    };
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(direction),
    );
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        normalizer.normalize(&variant)
    })) {
        Ok(Ok(v)) => v.to_string(),
        Ok(Err(_)) => "<declined>".to_string(),
        Err(_) => "<panic>".to_string(),
    }
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
            for axis in ["g", "c"] {
                assert!(
                    !inputs_for(family, axis, core).is_empty(),
                    "family {family} emitted no {axis}. rows"
                );
            }
        }
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
