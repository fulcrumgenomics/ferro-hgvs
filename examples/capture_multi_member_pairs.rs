//! Capture the multi-member cis alleles as reference/alternate pairs, so the
//! `from_sequences` axis over them runs **without a reference**.
//!
//! # Why this generator exists
//!
//! `tests/fixtures/cis/multi_member_cis_alleles.json` holds the 592 real-world
//! multi-member alleles harvested from ClinVar, CMRG and Paraphase — the entire
//! real-world evidence base for the code path `from_sequences` enters, and the
//! shape the downstream consumer actually submits. Judging them needed a
//! prepared reference, so the axis was manifest-gated and **skipped in CI**:
//! the only real-world coverage this surface has was invisible on every PR.
//!
//! It does not have to be. `from_sequences` reads no reference — that is its
//! defining property — so the only thing the axis needs from a provider is the
//! *window*, and a window is a value that can be committed. This captures those
//! values once, against the prepared reference, and the axis then runs anywhere.
//!
//! That is a smaller fixture than the `reference-windows.json` pattern used for
//! the projection corpora, and a stronger one: it stores exactly the bytes the
//! function under test consumes, so nothing about the capture can drift from
//! what the test does with it.
//!
//! ```sh
//! cargo run --features dev --example capture_multi_member_pairs -- \
//!   --manifest <manifest.json>
//! ```

use ferro_hgvs::conformance::completeness::CaptureLedger;
use ferro_hgvs::conformance::error_mode_stamp::ErrorModeStamp;
use ferro_hgvs::{parse_hgvs, MultiFastaProvider, Normalizer};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// Flank asked of `to_sequences`, matching `merge::CANONICAL_PAD` and every
/// sibling module in this family.
const PAD: u64 = 128;

/// Widest allele span captured, in bases.
///
/// The grid budget refuses well below this, so a wider row could only ever be
/// stored to be refused — at a cost of its whole window on disk. Rows past it
/// are counted rather than dropped silently; see [`Captured::too_wide`].
const MAX_SPAN: u64 = 4_096;

#[derive(Deserialize)]
struct Row {
    input: String,
}

#[derive(Deserialize)]
struct Fixture {
    rows: Vec<Row>,
}

#[derive(Serialize)]
struct Pair {
    input: String,
    accession: String,
    position: u64,
    reference: String,
    alternate: String,
    /// The row's encoding-invariant SPDI key, computed **here**, against the
    /// prepared reference, from the *input* description.
    ///
    /// This is the axis's independent **comparand**. Without it the test could
    /// only compare a derivation against `from_sequences`'s own internal
    /// round-trip check — the verifier and the verified being one piece of code,
    /// which is a self-consistency check that a bug in the applier passes.
    /// Captured from the input rather than from any derived form, so nothing
    /// about it can agree with a derivation by construction.
    ///
    /// **Capturing it is half the job; the consumer has to actually compare
    /// against it.** For a while it did not — `from_sequences_multi_member_axis`
    /// gated the comparison on a condition that was constantly false and read
    /// this field only for `is_empty()`, so the census reported rows as
    /// "compared" that nothing had judged. The consumer now keys the derived
    /// side through `spdi::canonical_spdi` over a window-only provider and
    /// compares the two keys as strings. If that changes again, this doc is
    /// wrong again: it describes a property of the pair, not of the file.
    ///
    /// `None` where `canonical_spdi` declines the row (an edit SPDI cannot
    /// carry, members that overlap); the test counts those rather than
    /// treating a missing key as a pass.
    canonical_spdi: Option<String>,
}

#[derive(Serialize)]
struct Captured {
    description: &'static str,
    generator: &'static str,
    source: &'static str,
    /// The error mode the capture ran under.
    ///
    /// Both captured values are mode-dependent, so a file without this says
    /// nothing about what it is comparable to. `to_sequences` and
    /// `canonical_spdi` each return `Err` on a row the mode rejects, and both
    /// failures are counted rather than stored — so a stricter mode moves rows
    /// from `pairs` into `unwindowed`, and a laxer one moves them back. Two
    /// captures taken under different modes would differ in their row *set*
    /// while looking like the same fixture.
    normalized_under: ErrorModeStamp,
    /// Rows in the source fixture.
    rows_scanned: usize,
    /// Rows whose span exceeds [`MAX_SPAN`]; never offered to the provider.
    too_wide: usize,
    /// Rows the provider could not serve a window for — an edit SPDI cannot
    /// carry, members on different accessions, a span past `MAX_APPLY_WINDOW`.
    unwindowed: usize,
    pairs: Vec<Pair>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let mut manifest = None;
    while let Some(arg) = args.next() {
        if arg == "--manifest" {
            manifest = args.next();
        }
    }
    let manifest = manifest.ok_or("usage: --manifest <manifest.json>")?;
    let provider = Arc::new(MultiFastaProvider::from_manifest(&manifest)?);
    let normalizer = Normalizer::new(Arc::clone(&provider));

    let text = std::fs::read_to_string("tests/fixtures/cis/multi_member_cis_alleles.json")?;
    let fixture: Fixture = serde_json::from_str(&text)?;

    let mut ledger = CaptureLedger::new("multi-member windows");
    let (mut too_wide, mut unwindowed) = (0usize, 0usize);
    let mut pairs = Vec::new();

    for row in &fixture.rows {
        let Some(variant) = ledger.record(row.input.as_str(), parse_hgvs(&row.input)) else {
            continue;
        };
        if declared_span(&row.input).is_none_or(|span| span > MAX_SPAN) {
            too_wide += 1;
            continue;
        }
        match normalizer.to_sequences(&variant, PAD) {
            Ok(pair) => pairs.push(Pair {
                input: row.input.clone(),
                canonical_spdi: normalizer
                    .canonical_spdi(&variant)
                    .ok()
                    .map(|s| s.to_string()),
                accession: pair.accession,
                position: pair.position,
                reference: pair.reference,
                alternate: pair.alternate,
            }),
            Err(_) => unwindowed += 1,
        }
    }
    let counts = ledger.finish()?;
    eprintln!("parsed {} of {} rows", counts.succeeded, counts.attempted);

    let captured = Captured {
        description: "Reference/alternate windows for the real-world multi-member cis alleles, \
                      captured against a prepared reference so the from_sequences axis over them \
                      needs no provider. from_sequences reads no reference, so the window is the \
                      only thing it needs and a window is a value.",
        generator:
            "cargo run --features dev --example capture_multi_member_pairs -- --manifest <m>",
        source: "tests/fixtures/cis/multi_member_cis_alleles.json",
        normalized_under: ErrorModeStamp::of(&normalizer.config().error_config),
        rows_scanned: fixture.rows.len(),
        too_wide,
        unwindowed,
        pairs,
    };
    eprintln!(
        "captured {} pairs ({} too wide, {} unwindowed)",
        captured.pairs.len(),
        captured.too_wide,
        captured.unwindowed
    );
    // Gzipped: the pretty JSON is ~740 KB, past the repo's 500 KB pre-commit
    // limit, and compresses to ~100 KB — small enough to commit plainly, so no
    // LFS pointer and no test that skips green when the fixture is absent.
    let json = serde_json::to_string_pretty(&captured)? + "\n";
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::best());
    std::io::Write::write_all(&mut encoder, json.as_bytes())?;
    std::fs::write(
        "tests/fixtures/cis/multi_member_sequence_pairs.json.gz",
        encoder.finish()?,
    )?;
    Ok(())
}

/// The row's total span in bases, read off its own coordinates.
///
/// Derived from the string rather than by fetching, because this is the
/// pre-filter that decides whether to spend a reference read at all.
fn declared_span(input: &str) -> Option<u64> {
    let body = input.split_once(':')?.1;
    let digits: Vec<u64> = body
        .split(|c: char| !c.is_ascii_digit())
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse().ok())
        .collect();
    Some(digits.iter().max()? - digits.iter().min()?)
}
