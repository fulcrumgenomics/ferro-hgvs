//! Issue #1970 — the `u16` alignment cost grid refuses instead of wrapping, and
//! the refusal is reachable from the ordinary public API.
//!
//! `AlignmentDag`'s cost grids narrowed to `u16`. The bound that makes that safe
//! is `n + m <= u16::MAX`, and it was enforced by a `debug_assert` — which
//! `[profile.release]` compiles out, since it sets neither `debug-assertions`
//! nor `overflow-checks`. Past the bound the shipped library therefore did not
//! panic: `edit_grid`'s boundary row truncated `j as u16` and the
//! `prefix + suffix` sums wrapped, so `on_optimal_path` was computed from wrong
//! costs and the emitted description came off a wrong DAG.
//!
//! Every gate in front of `AlignmentDag::build` bounds **cells** —
//! `(n + 1) * (m + 1)` — never `n + m`, and a `1 x 70,000` grid clears the
//! widest of them by four orders of magnitude. A long literal `ins` payload is
//! exactly that shape, and it arrives through `parse_hgvs` +
//! `Normalizer::normalize` with no feature flag and no `FERRO_PARTITION` arm.
//!
//! The bound itself is pinned from both sides as a unit test next to the code,
//! `align::tests::the_span_bound_is_enforced_from_both_sides`. What this module
//! adds is the half that cannot be seen from there: that the size is **reached**
//! from the public API, that the refusal is taken there rather than the wrap,
//! and that taking it moves no output.

use ferro_hgvs::normalize::partition_decline_counts;
use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

/// Deterministic pseudo-random bases.
///
/// Not `"ACGT"` cycled: a periodic payload is repeat-lowered long before the
/// aligner sees it (`MAX_CANONICAL_WINDOW` caps the synthesized payload there),
/// so a cyclic 70,000-mer normalizes to a 21-character repeat description and
/// never builds a grid at all. Measured while writing this test — the first
/// version of it passed on the unfixed code for exactly that reason.
fn prng(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed | 1;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

/// A genome-capable provider holding one 4,000-base pseudo-random contig.
fn provider() -> JsonProvider {
    let seq = String::from_utf8(prng(4_000, 99)).expect("ACGT is UTF-8");
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { "c": seq },
    });
    let mut file = tempfile::NamedTempFile::new().expect("temp file");
    file.write_all(doc.to_string().as_bytes())
        .expect("write provider");
    JsonProvider::from_json(file.path()).expect("load provider")
}

/// `c:g.1000_1001ins<payload_len pseudo-random bases>`.
fn insertion(payload_len: usize) -> String {
    let payload = String::from_utf8(prng(payload_len, 7)).expect("ACGT is UTF-8");
    format!("c:g.1000_1001ins{payload}")
}

fn normalize(provider: &JsonProvider, input: &str) -> String {
    let variant = parse_hgvs(input).expect("the input parses");
    Normalizer::new(provider.clone())
        .normalize(&variant)
        .expect("the input normalizes")
        .to_string()
}

/// The payload length at which `n + m` first exceeds the `u16` grid.
///
/// The reference block is trimmed to nothing by the time the aligner sees a pure
/// insertion, so `n + m` is the payload length: 65,535 is the last accepted size
/// and 65,536 the first refused one.
const FIRST_REFUSED_PAYLOAD: usize = u16::MAX as usize + 1;

/// A payload one base under the bound must still be aligned — the refusal must
/// not have been made conservative enough to swallow inputs that are fine.
///
/// This is the *under* side of the boundary, and it is the side a careless fix
/// gets wrong: the `debug_assert` this replaced demanded `2 * (n + m)` fit, so
/// it would have refused everything past 32,767 — half the representable range,
/// on inputs whose arithmetic never overflows.
#[test]
fn a_payload_just_under_the_span_bound_is_still_aligned() {
    let provider = provider();
    let before = partition_decline_counts();
    let out = normalize(&provider, &insertion(FIRST_REFUSED_PAYLOAD - 1));
    let after = partition_decline_counts();

    assert!(
        after.attempted > before.attempted,
        "the sequence-first partitioner must have been asked at all — \
         a zero here would make the test vacuous ({before:?} -> {after:?})"
    );
    assert_eq!(
        after.declined,
        before.declined,
        "n + m = {} is inside the u16 grid and must NOT be declined ({before:?} -> {after:?})",
        FIRST_REFUSED_PAYLOAD - 1
    );
    assert_eq!(
        out.len(),
        FIRST_REFUSED_PAYLOAD - 1 + 16,
        "the whole payload must survive into the description"
    );
}

/// A payload one base over the bound must be **declined by the aligner**, and
/// the decline must be invisible in the output.
///
/// The assertion is on the refusal itself — the `declined` census — not on the
/// emitted string, which is the convenient side effect. `partition_block_*`
/// answers `None` and `partition_block_for_rule` falls back, so an input past
/// the bound still gets a description; that is what makes the output a bad
/// subject and the counter a good one.
#[test]
fn a_payload_just_over_the_span_bound_is_declined_by_the_aligner() {
    let provider = provider();
    let before = partition_decline_counts();
    let out = normalize(&provider, &insertion(FIRST_REFUSED_PAYLOAD));
    let after = partition_decline_counts();

    assert!(
        after.attempted > before.attempted,
        "the sequence-first partitioner must have been asked at all ({before:?} -> {after:?})"
    );
    assert!(
        after.declined > before.declined,
        "n + m = {} is past the u16 grid and must be DECLINED, not computed from \
         wrapped costs ({before:?} -> {after:?})",
        FIRST_REFUSED_PAYLOAD
    );
    assert_eq!(
        out.len(),
        FIRST_REFUSED_PAYLOAD + 16,
        "the fallback partitioner must still describe the whole payload"
    );
}

/// The size #1970 reported, end to end, with the description pinned exactly.
///
/// The pin is what makes `Representation-Change: none` a measurement rather than
/// an argument: this string was captured by running the identical test against
/// the PR's own pre-narrowing parent (`92cd3f18`, `u32` grids), where the whole
/// 70,000-base grid is computed rather than declined. Both arms emit it
/// byte-for-byte, so refusing the grid costs this class nothing.
///
/// **This arm is a stability pin, not a wrap detector, and the difference is
/// worth stating.** Measured under sabotage (refusal removed, `debug-assertions`
/// and `overflow-checks` off): the wrap at 70,000 happens and this test still
/// passes. A pure insertion trims the reference block to nothing, so the grid is
/// a single row with exactly one path through it, and `prefix + suffix == total`
/// holds under wrapping arithmetic just as it does without — the truncation is
/// benign at *this* shape. What refuses it is
/// `a_payload_just_over_the_span_bound_is_declined_by_the_aligner`, and what
/// pins the bound is the unit test in `align.rs`; those two are the guards that
/// go red. This one says the refusal is free.
#[test]
fn the_reported_seventy_kilobase_insertion_is_unchanged_by_the_narrowing() {
    let provider = provider();
    let out = normalize(&provider, &insertion(70_000));

    assert_eq!(out.len(), 70_016, "output length, measured on both arms");
    assert!(
        out.starts_with("c:g.") && out.contains("ins"),
        "still a genomic insertion: {}...",
        &out[..out.len().min(32)]
    );

    // The whole 70,016-character string, compared without printing it.
    let mut hash = 0xcbf29ce484222325u64;
    for byte in out.as_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x100000001b3);
    }
    assert_eq!(
        format!("{hash:016x}"),
        "e837728b1be3a7c4",
        "the emitted description must be the one the u32 grids produced"
    );
}
