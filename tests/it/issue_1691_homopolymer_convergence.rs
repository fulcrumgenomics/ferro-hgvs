//! Issue #1691 — `Normalizer::normalize` never reaches a fixed point inside a
//! homopolymer longer than the reference fetch window.
//!
//! Root cause: `normalize_genome` fetches a `±window_size` (default 100 bp)
//! window around the variant and hands `normalize_na_edit` a
//! `Boundaries::new(0, ref_seq.len())`. That right bound is *where the fetch
//! stopped*, not a property of the contig, so inside a long homopolymer the 3'
//! shuffle runs to the window edge and reports the position it happened to
//! reach. Re-normalizing re-centres the window on that position and the shuffle
//! advances another `window_size` bases:
//!
//! ```text
//! c:g.65del -> c:g.165del -> c:g.265del -> c:g.365del -> ...   (step == 100)
//! ```
//!
//! **`MAX_SEQUENCE_COMPARE_WINDOW` (100_000) is a red herring**, named in the
//! issue only because the PR that surfaced this built a homopolymer that long.
//! The bound that truncates the shuffle is `NormalizeConfig::window_size`, three
//! orders of magnitude smaller — and the step size in the walk above is exactly
//! that constant, which is the measurement that identifies it.
//!
//! **A short tract does not reproduce it, and that near-miss is worth pinning.**
//! `canonicalize_from_sequence` alternates the sequence-first pass (which fetches
//! its own `±CANONICAL_PAD = 128` window) with `normalize_core` up to
//! `MAX_PASSES = 4` times, so a tract the alternation can cross end-to-end is
//! rescued: a 200-base and a 600-base run both converge on unfixed code. The
//! defect starts at the first tract the four alternations cannot span, which is
//! **between 700 and 800 bases** — measured on `origin/main` by counting calls to
//! a fixed point from `c:g.65del`:
//!
//! ```text
//! run=600  1 call    run=800  2 calls    run=1000  4 calls
//! run=700  1 call    run=900  3 calls    run=2000 12 calls
//! ```
//!
//! `TRACT` below is 2000 rather than 800 so it sits comfortably past that
//! boundary rather than on it, and rather than the 100_000-plus the issue
//! reports. (This module previously bracketed the onset at "between 600 and
//! 2000", which is true but loose; the 600 in
//! `issue_1691_short_tracts_are_unaffected` is still correctly chosen — it is
//! inside the converging range, not on its edge.)
//!
//! The fix grows the fetch window geometrically until the result rests strictly
//! inside it, capped by `MAX_SHUFFLE_FETCH_WINDOW` (64 KiB). Past the cap the
//! **shift** is refused outright rather than truncated, because a truncated shift
//! only changes the walk's step size — see `normalize_in_grown_window`. The
//! **edit type** is not refused with it: `canonicalize_without_shifting` still
//! re-types the edit against the reference, because that verdict is decided by
//! the bases flanking the edit and so survives the cap.

use ferro_hgvs::{parse_hgvs, JsonProvider, Normalizer};
use std::io::Write;

/// Long enough that `canonicalize_from_sequence`'s four alternations cannot span
/// it, so the walk is visible. See the module comment for the measurement.
const TRACT: usize = 2000;

/// Flanking bases either side of the tract. A multiple of 4 so the cyclic
/// `ACGT` lead ends on `T` — never on the tract's own base, which would extend
/// the run and move its 3' end.
const FLANK: usize = 64;

/// Genome-capable provider serving one contig: `flank` bases of cyclic `ACGT`,
/// then a run of `run` `A`s, then `flank` bases of cyclic `CGTA`.
///
/// The trailing flank starts at `C` rather than `A` for the same reason the
/// leading one ends at `T`: the tract's 3' end must be exactly `flank + run`.
fn homopolymer_provider(contig: &str, flank: usize, run: usize) -> JsonProvider {
    let lead: String = "ACGT".chars().cycle().take(flank).collect();
    let tract: String = std::iter::repeat_n('A', run).collect();
    let trail: String = "CGTA".chars().cycle().take(flank).collect();
    genomic_provider(contig, &format!("{lead}{tract}{trail}"))
}

fn genomic_provider(contig: &str, sequence: &str) -> JsonProvider {
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { contig: sequence },
    });
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(file.path()).unwrap()
}

fn norm(provider: &JsonProvider, input: &str) -> String {
    let variant = parse_hgvs(input).unwrap();
    Normalizer::new(provider.clone())
        .normalize(&variant)
        .unwrap()
        .to_string()
}

/// The defect as the issue states it: repeated normalization must reach a fixed
/// point instead of advancing one window per call.
///
/// Stated as the *chain* rather than as a single equality, so a failure prints
/// the walk and its step size — which is the evidence that identifies the bound
/// being hit.
#[test]
fn issue_1691_normalize_reaches_a_fixed_point_inside_a_long_homopolymer() {
    let provider = homopolymer_provider("c", FLANK, TRACT);

    let mut chain = vec![format!("c:g.{}del", FLANK + 1)];
    for _ in 0..8 {
        let next = norm(&provider, chain.last().unwrap());
        if &next == chain.last().unwrap() {
            break;
        }
        chain.push(next);
    }

    assert!(
        chain.len() <= 2,
        "normalize did not converge; it walked {chain:?}",
    );
}

/// …and the fixed point it reaches is the 3'-most position in the tract, which
/// is what the 3' rule asks for. A test that only asserted convergence would
/// pass on a fix that simply stopped shifting.
#[test]
fn issue_1691_the_fixed_point_is_the_three_prime_end_of_the_tract() {
    let provider = homopolymer_provider("c", FLANK, TRACT);
    let three_prime_most = format!("c:g.{}del", FLANK + TRACT);

    for authored in [FLANK + 1, FLANK + 2, FLANK + TRACT / 2, FLANK + TRACT] {
        assert_eq!(
            norm(&provider, &format!("c:g.{authored}del")),
            three_prime_most,
            "c:g.{authored}del did not reach the 3' end of the tract",
        );
    }
}

/// The insertion and duplication shuffles are driven by the same window bound
/// and walk the same way; both must land on a fixed point too.
#[test]
fn issue_1691_insertions_and_duplications_converge_too() {
    let provider = homopolymer_provider("c", FLANK, TRACT);

    for input in [
        format!("c:g.{}_{}insA", FLANK + 1, FLANK + 2),
        format!("c:g.{}dup", FLANK + 1),
    ] {
        let once = norm(&provider, &input);
        let twice = norm(&provider, &once);
        assert_eq!(once, twice, "{input} did not converge ({once} -> {twice})");
    }
}

/// The 3' growth must stop at the contig end rather than reading past it — the
/// #1041 failure mode, which a naive "just fetch more" fix would reintroduce.
#[test]
fn issue_1691_growth_stops_at_the_contig_end() {
    // No trailing flank: the tract is the last thing on the contig.
    let lead: String = "ACGT".chars().cycle().take(FLANK).collect();
    let tract: String = std::iter::repeat_n('A', TRACT).collect();
    let provider = genomic_provider("c", &format!("{lead}{tract}"));

    let last = FLANK + TRACT;
    assert_eq!(
        norm(&provider, &format!("c:g.{}del", FLANK + 1)),
        format!("c:g.{last}del"),
    );
    assert_eq!(
        norm(&provider, &format!("c:g.{last}del")),
        format!("c:g.{last}del"),
    );
}

/// A tract past `MAX_SHUFFLE_FETCH_WINDOW` (64 KiB) exhausts the growth cap.
/// The answer there is a **refusal to shift**, not a shift truncated at the cap
/// — a truncated shift would merely change the walk's step size from 100 to
/// 65_536 and the description would still never settle.
///
/// This is the case the issue's own 200_000-base reproducer falls into, so it is
/// the one that decides whether `FERRO_ASSERT_IDEMPOTENT` stays quiet on it.
#[test]
fn issue_1691_a_tract_past_the_growth_cap_refuses_to_shift_and_is_stable() {
    // 70_000 > 64 KiB, so the cap is reached with contig left to read.
    let provider = homopolymer_provider("c", FLANK, 70_000);

    let authored = format!("c:g.{}del", FLANK + 1);
    let once = norm(&provider, &authored);
    assert_eq!(
        once, authored,
        "past the growth cap the shift must be refused, not truncated",
    );
    assert_eq!(norm(&provider, &once), once, "the refusal must be stable");
}

/// **Refusing the shift must not also refuse the re-typing.** `c:g.65_66insA`
/// inserts an `A` directly 3' of an `A`, which is a tandem duplication:
/// `DNA/duplication.md:18` — "when a variant can be described as a duplication,
/// it **must** be described as a duplication and not as, e.g., an insertion" —
/// and `DNA/insertion.md:17` renders the insertion spelling as
/// `<code class="invalid">`. That verdict is decided by the bases immediately
/// flanking the insertion point, not by how far the tract runs, so it is
/// available at the growth cap and must still be given there.
///
/// Refusing it traded rule 1 (conformant output, README:176-179, "absolute —
/// never traded") for idempotency, which is not one of the four output rules at
/// all. The capped answer has to be **both**: a `dup`, and a fixed point.
///
/// Both axes, because the growth loop is shared and a per-axis fix would be the
/// drift `clamp_fetch_end_to_contig` already exists to prevent.
#[test]
fn issue_1691_past_the_growth_cap_an_insertion_is_still_typed_as_a_duplication() {
    // 70_000 > 64 KiB, so the cap is reached with contig left to read.
    for (accession, axis) in [("c", "g"), ("NC_012920.1", "m")] {
        let provider = homopolymer_provider(accession, FLANK, 70_000);
        let authored = format!("{accession}:{axis}.{}_{}insA", FLANK + 1, FLANK + 2);

        let once = norm(&provider, &authored);
        assert!(
            once.ends_with("dup"),
            "{authored} is a tandem duplication; past the growth cap it was \
             emitted as {once}, the spelling `DNA/insertion.md:17` marks invalid",
        );
        assert_eq!(
            norm(&provider, &once),
            once,
            "the re-typed answer must still be a fixed point ({authored} -> {once})",
        );
    }
}

/// The rescue path this defect hid behind must keep working: a tract the
/// sequence-first alternation can span end-to-end converged before the fix and
/// must still converge after it.
#[test]
fn issue_1691_short_tracts_are_unaffected() {
    for run in [8_usize, 60, 200, 600] {
        let provider = homopolymer_provider("c", FLANK, run);
        let three_prime_most = format!("c:g.{}del", FLANK + run);
        assert_eq!(
            norm(&provider, &format!("c:g.{}del", FLANK + 1)),
            three_prime_most,
            "a {run}-base tract regressed",
        );
    }
}

/// The same bound sits on the mitochondrial path, which shares the fetch with
/// `normalize_genome`. mtDNA carries no run anywhere near this long, so this
/// pins the shared helper rather than a real `m.` input.
#[test]
fn issue_1691_the_mitochondrial_path_converges_too() {
    let lead: String = "ACGT".chars().cycle().take(FLANK).collect();
    let tract: String = std::iter::repeat_n('A', TRACT).collect();
    let trail: String = "CGTA".chars().cycle().take(FLANK).collect();
    let provider = genomic_provider("NC_012920.1", &format!("{lead}{tract}{trail}"));

    let once = norm(&provider, &format!("NC_012920.1:m.{}del", FLANK + 1));
    assert_eq!(
        norm(&provider, &once),
        once,
        "the mitochondrial path did not converge",
    );
}
