//! Dump all three block partitions side by side, one row per block.
//!
//! The bake-off harness for ferro's three partitioners. Each input row is a
//! `(reference_block, alternate_block)` pair; each output row is what `live`,
//! `shadow` and `canonical` each make of it, together with the block's
//! Levenshtein distance and each arm's changed-column count so distance-minimality
//! is checkable per arm straight from the dump.
//!
//! The three rules are not variations on one another:
//!
//! * **live** — `partition_block`: a single-gap alignment search plus two narrow
//!   escapes. The shipped rule.
//! * **shadow** — `partition_block_sequence_first`: cut at the alignment steps
//!   common to *every* minimal alignment, then merge runs separated by fewer than
//!   `--min-separation` unchanged bases. `mutalyzer-algebra`'s `local_supremal`.
//! * **canonical** — `partition_block_canonical`: the minimal alignment with the
//!   fewest members, 3'-most among ties. `mutalyzer-algebra`'s `canonical`.
//!
//! # Usage
//!
//! ```text
//! cargo run --features dev --example dump_partitions -- [--min-separation N] [PATH]
//! ```
//!
//! `PATH` defaults to `-` (standard input). Input is tab-separated
//! `reference<TAB>alternate`, one block pair per line; empty lines and lines
//! beginning with `#` are skipped, and either block may be empty. Bases are
//! compared byte-exactly and are **not** case-folded — pass the blocks in the
//! same case the normalizer would, which is upper case.
//!
//! `--min-separation` applies to the `shadow` arm only; the other two take no
//! such parameter. It defaults to the coding one-amino-acid exception (2). Pass
//! `1` to measure the arm as a genomic axis would run it.
//!
//! # Output
//!
//! One header row, then one tab-separated row per input block:
//!
//! ```text
//! ref  alt  distance  live  live_cols  shadow  shadow_cols  canonical  canonical_cols
//! ```
//!
//! A partition renders as its members joined by `|`, each as
//! `ref_start:ref_end/alt` over 0-based half-open reference offsets into the
//! block — so `0:1/|3:4/` is two deletions and `2:2/CAC` is one insertion. An
//! empty partition (an unchanged block) renders as `-`.
//!
//! Two sentinels replace a partition rather than a member list. `DECLINED` is an
//! arm that legitimately returned nothing — the `shadow` arm does on some blocks
//! — and `PANIC` is an arm that panicked, which is caught per arm so one bad
//! block cannot abandon the dump. Their column counts render as `-`.

use ferro_hgvs::normalize::dev_partitioners::{self, DevPiece, DEFAULT_MIN_SEPARATION};
use std::io::{BufRead, BufWriter, Write};
use std::panic::{catch_unwind, AssertUnwindSafe};

/// What one partitioner made of one block.
enum Outcome {
    Partitioned(Vec<DevPiece>),
    /// The arm returned `None`.
    Declined,
    /// The arm panicked. Recorded rather than propagated, so a dump over a large
    /// corpus reports the bad block instead of stopping at it.
    Panicked,
}

impl Outcome {
    /// Run one arm, converting a panic into [`Outcome::Panicked`].
    fn of<F: FnOnce() -> Option<Vec<DevPiece>>>(arm: F) -> Self {
        match catch_unwind(AssertUnwindSafe(arm)) {
            Ok(Some(pieces)) => Outcome::Partitioned(pieces),
            Ok(None) => Outcome::Declined,
            Err(_) => Outcome::Panicked,
        }
    }

    /// The partition as one deterministic field.
    fn render(&self) -> String {
        match self {
            Outcome::Partitioned(pieces) if pieces.is_empty() => "-".to_string(),
            Outcome::Partitioned(pieces) => pieces
                .iter()
                .map(|piece| {
                    format!(
                        "{}:{}/{}",
                        piece.ref_start,
                        piece.ref_end,
                        String::from_utf8_lossy(&piece.alt)
                    )
                })
                .collect::<Vec<_>>()
                .join("|"),
            Outcome::Declined => "DECLINED".to_string(),
            Outcome::Panicked => "PANIC".to_string(),
        }
    }

    /// Changed columns the partition claims, or `-` when there is no partition.
    fn columns(&self) -> String {
        match self {
            Outcome::Partitioned(pieces) => dev_partitioners::changed_columns(pieces).to_string(),
            Outcome::Declined | Outcome::Panicked => "-".to_string(),
        }
    }
}

/// Command-line arguments: an optional input path, the shadow arm's separation
/// threshold, and which axis the live arm is asked about.
struct Args {
    path: String,
    min_separation: u32,
    /// Whether the live arm may apply `DNA/delins.md:44-47`'s
    /// payload-coincidence carve-out, which the operator ruling
    /// `delins-payload-coincidence-carve-out-is-coding-dna-scoped` scopes to
    /// `c.` and nothing else.
    ///
    /// A TSV row is two bare blocks and names no axis, so the default is the
    /// frameless reading — `general.md:34` governs and a split across one
    /// unchanged base stands. `--coding-dna` asks the other question, which is
    /// the only way to dump both sides of the ruling from one corpus.
    coding_dna: bool,
}

fn parse_args() -> Result<Args, String> {
    // Tracked separately from `path`'s default so a second positional argument
    // is an error rather than a silent overwrite: `dump_partitions a.tsv b.tsv`
    // would otherwise read `b.tsv` and never mention `a.tsv`, which on a corpus
    // run looks like the first file simply contained nothing of interest.
    let mut path: Option<String> = None;
    let mut min_separation = DEFAULT_MIN_SEPARATION;
    let mut coding_dna = false;
    let mut argv = std::env::args().skip(1);
    while let Some(arg) = argv.next() {
        match arg.as_str() {
            "--min-separation" => {
                let value = argv
                    .next()
                    .ok_or_else(|| "--min-separation needs a value".to_string())?;
                min_separation = value
                    .parse()
                    .map_err(|_| format!("--min-separation: {value} is not a number"))?;
            }
            "--coding-dna" => coding_dna = true,
            "-h" | "--help" => {
                return Err(
                    "usage: dump_partitions [--min-separation N] [--coding-dna] [PATH]".to_string(),
                )
            }
            other if other.starts_with("--") => return Err(format!("unknown flag {other}")),
            other if path.is_some() => {
                return Err(format!(
                    "unexpected second PATH {other}; dump_partitions reads one file (or `-`)"
                ))
            }
            other => path = Some(other.to_string()),
        }
    }
    Ok(Args {
        path: path.unwrap_or_else(|| "-".to_string()),
        min_separation,
        coding_dna,
    })
}

fn main() {
    let args = match parse_args() {
        Ok(args) => args,
        Err(message) => {
            eprintln!("{message}");
            std::process::exit(2);
        }
    };

    // One terse line per panicking block instead of the default hook's payload
    // and backtrace: a dump over a large corpus is a table, and a hook that
    // writes paragraphs to stderr for every bad row buries the rows that matter.
    // The block itself is already identified by the `PANIC` cell in the table.
    std::panic::set_hook(Box::new(|info| {
        eprintln!(
            "panic: {}",
            info.location().map_or(String::new(), |l| l.to_string())
        );
    }));

    let input: Box<dyn BufRead> = if args.path == "-" {
        Box::new(std::io::stdin().lock())
    } else {
        match std::fs::File::open(&args.path) {
            Ok(file) => Box::new(std::io::BufReader::new(file)),
            Err(error) => {
                eprintln!("cannot read {}: {error}", args.path);
                std::process::exit(2);
            }
        }
    };

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    writeln!(
        out,
        "ref\talt\tdistance\tlive\tlive_cols\tshadow\tshadow_cols\tcanonical\tcanonical_cols"
    )
    .expect("write header");

    for line in input.lines() {
        // Exit rather than panic: a read or UTF-8 error part-way through a corpus
        // leaves a truncated TSV, and a panic's backtrace buries that fact under
        // noise. Status 2 distinguishes it from the usage errors above (1).
        let line = match line {
            Ok(line) => line,
            Err(error) => {
                eprintln!("cannot read input: {error}");
                std::process::exit(2);
            }
        };
        // `is_empty`, not `trim().is_empty()`: a lone tab is the legitimate
        // empty-to-empty block pair, and trimming would silently drop it.
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((reference, alt)) = line.split_once('\t') else {
            eprintln!("skipping line without a tab: {line}");
            continue;
        };
        let (reference, alt) = (reference.as_bytes(), alt.as_bytes());

        // `live` cannot decline, but it can panic, so it goes through the same
        // wrapper with a `Some` of its own.
        let live = Outcome::of(|| Some(dev_partitioners::live(reference, alt, args.coding_dna)));
        let shadow = Outcome::of(|| dev_partitioners::shadow(reference, alt, args.min_separation));
        let canonical = Outcome::of(|| dev_partitioners::canonical(reference, alt));
        let distance = match catch_unwind(AssertUnwindSafe(|| {
            dev_partitioners::edit_distance(reference, alt)
        })) {
            Ok(distance) => distance.to_string(),
            Err(_) => "-".to_string(),
        };

        writeln!(
            out,
            "{}\t{}\t{distance}\t{}\t{}\t{}\t{}\t{}\t{}",
            String::from_utf8_lossy(reference),
            String::from_utf8_lossy(alt),
            live.render(),
            live.columns(),
            shadow.render(),
            shadow.columns(),
            canonical.render(),
            canonical.columns(),
        )
        .expect("write row");
    }
    out.flush().expect("flush output");
}
