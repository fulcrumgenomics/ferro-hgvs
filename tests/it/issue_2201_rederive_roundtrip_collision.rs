//! Regression + exploration coverage for issue #2201.
//!
//! `Normalizer::rederive` regressed in v0.17.1: for certain clustered
//! multi-edit genomic (`g.`) inputs it derives a description it then cannot
//! re-apply to its own window, raising
//!
//! ```text
//! Coordinate conversion error: derived description <...> could not be
//! re-applied to its window
//! ```
//!
//! Bisected to `0db1ec61` (#2194/#2199), the `coalesce_by_run` two-license
//! rewrite. The shape all eleven reported variants share is a run of adjacent
//! single-base deletions immediately followed by one net-positive insertion; the
//! new coalesce collapses that run into a form that places a `dup` and an `ins`
//! on the SAME junction (e.g. `342_343dup;343_344insCGGAGG`). Both write the
//! reference's `after[343]` slot, so `apply_edits_to_window`'s coincident-slot
//! guard refuses the combination and the round-trip check fails.
//!
//! Authored test-first (RED), then made green by the `merge_coincident_insertions`
//! fix in this same PR; they now stand as regression guards. Layout:
//!
//! * **Part A** — the exact reported geometry (must round-trip).
//! * **Part B** — 1-D sweeps mapping the RED/GREEN boundary around it.
//! * **Tier 1** — the trigger reverse-engineered (`T1a`/`T1c`), a 2-D grid,
//!   both shuffle directions, idempotency.
//! * **Tier 2** — confluence, net-sign symmetry, window-perturbation stability,
//!   a structured window corpus.
//! * **Tier 3** — a neighbourhood fuzzer proven able to fail, and an opt-in
//!   corpus replay.
//! * **Tier 4** — shape coverage (#2206) across the regression's sibling
//!   geometries the del-run+ins tests do not reach: near-miss anchors, a `dup`
//!   5' member, a del-run-then-`dup`, and the `m.` axis. These pin the derived
//!   SPELLING, which `verify_round_trip` (guaranteeing only the bases) cannot.
//!
//! Each exploratory test prints a per-cell RED/GREEN map before asserting, so a
//! single run shows exactly where the boundary sits.

use ferro_hgvs::{parse_hgvs, FromSequencesOptions, JsonProvider, Normalizer, ShuffleDirection};
use std::io::Write;

// ===========================================================================
// Constants
// ===========================================================================

/// The synthetic contig from the issue's standalone repro. Fully synthetic — no
/// real reference bases (verified in the issue). Bases of interest:
/// `335..=342 = TTTCCTTC`, `343 = T`, `344 = T`.
const TEMPLATE: &str = "AGGAAAAAAAAAAGGAAGAAGAAAGGAAGGCAGGGACACGACAGCAAAACACAAGGAAACAAACGCAACACACACAAGACACACACAGCAAGGAAACAGAAAGACACGCAGAACAGAAGAAGAAGAACGCACACGAACAAGGAAAAAACACCCAGGACACCACGAAAAACAGGGAAAGAAAAAACAAAAAAACAAGGAACCACAAGGAGGAGACACAACAACAAGAAGGGACGACGAACAAGAAAAACAGAACAACGCAAAAACAACCACAAAGGCGAAGACCACAAAAAGAAACAGAAAACGAGAACCCGAAACAGAGACAACAAAGAAACACTTTCCTTCTTCTAAGCACAAAGAAGCAGAAACAGACAAGGGAAAGAAAAAACCGAAAAAACACCAAAAAACCCCAAACCCAAAAAACACAAGAAAAACAGAAGAAGGAACAAGACCAAGGCAAGAAAAGCAAAGAAGAAAAAGAAGAAACAAGAAAAAGCCGACACGAAAAACAACAACAAAAAAACAAAAAAGAAAAAGCCAAAGAAGACCACAAACACGCGGAAAACACACACAGCAACGAAAAAACACACCAAACAAAAAACCAAAAAAAAGACCAAGAACAAACACCAAGACAACGACCAGGAGCAGAACGAAAAAGAAAACAACAAAAACAGGCGAGCACACGAACGCAGACACAAAGAGAAAAGGAAAAGAAAAAGAGAAAACAGACACACAAACAAGAAAAACAGAAAAACACACAACAAAGAAAACACCGAAAAGGAAGAGAAAAAGAAGAACAGCCAGACGAGAACCAAAAACAACCGAAACAACAGGAAGGCGCAGGAAAAAAACCAGAACAAAAAAAAGAAAGGGAAGAAGAAACAGAAGAGAAAACAGAGCAAAAGGAGGCGCAGAGCGAGAAACAACACAAACAAACACAAGAGAAACAAGAAGCGCAAGGAGAGCAAACAACAACACAAAAAGACAAGCAAAAAAGCAAAAAAGGACAACACAAGAGCACCACACGGGGAGAGAAACAAGAGGAACACGAAACGCAGCACCAGCAACAAGAAACAGAGACGACAAACAAGAAACCGAGAAAGAAAAAGCCAAGAGAGACGCGCCACCACATGAAAAAAGCACAGCAAGAACAAACGAACAGAGAACGAAGAAGAAGAGAAAAGAAAGGACCAAGAAAAGGCAGAAGCAAAAAGCAACAAACGAAAGAAACCGCGGAACAGACCAAAACACAGCAAAGGCGACAGCAGAGGACA";

/// The 30-base insertion payload from the tightest reported case (#6/#3/#4/…).
/// `PAYLOAD30 == "CCCGAACAAAGAAA" + TRIGGER_TAIL`.
const PAYLOAD30: &str = "CCCGAACAAAGAAATCACGTCTCTCGGAGG";

/// The 3' portion of `PAYLOAD30` that actually drives the collision, isolated by
/// reverse-engineering (see `t1a`). The 14-base prefix becomes the leading
/// `336delins…`; this 16-base tail mis-aligns against the reference's 3' context
/// to expose the `342_343dup;343_344ins` collision.
const TRIGGER_TAIL: &str = "TCACGTCTCTCGGAGG";

/// The interior window (1-based, inclusive) used by the `from_sequences`-level
/// tests: ±60 bases of flank around the `335..=344` hotspot, wide enough that no
/// member rests on a window edge.
const HOTSPOT_WINDOW: (u64, u64) = (275, 404);

/// The eleven synthetic variants from the issue thread — each collided before
/// the fix and must now round-trip.
const REPORTED_VARIANTS: [&str; 11] = [
    "NC_TEST.1:g.[301A>C;335delT;336delT;337delT;338delC;339delC;340delT;341delT;342delC;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG]",
    "NC_TEST.1:g.[335del;336del;337del;338del;339del;340del;341del;342del;343_344insCCAGAACAAAGAAAAAACGTCTCTCGGAGG;1030A>C;1046C>A;1085C>A]",
    "NC_TEST.1:g.[335del;336del;337del;338del;339del;340del;341del;342del;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG;489A>C;491A>C]",
    "NC_TEST.1:g.[335del;336del;337del;338del;339del;340del;341del;342del;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG;622C>T]",
    "NC_TEST.1:g.[335del;336del;337del;338del;339del;340del;341del;342del;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG;833A>C;914C>A]",
    "NC_TEST.1:g.[335delT;336delT;337delT;338delC;339delC;340delT;341delT;342delC;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG]",
    "NC_TEST.1:g.[1123_1124insACAC;1124C>A;1129T>A;1223A>C]",
    "NC_TEST.1:g.[11A>C;276C>A;334_335insCCCCGAACAAAGAAATCACGTCTCTCGGAGG;336del;337del;338del;339del;340del;341del;342del;343del;344del]",
    "NC_TEST.1:g.[69del;70del;315C>G;324C>T;335del;336del;337del;338del;339del;340del;341del;342del;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG]",
    "NC_TEST.1:g.[74A>C;335del;336del;337del;338del;339del;340del;341del;342del;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG]",
    "NC_TEST.1:g.[204A>C;335del;336del;337del;338del;339del;340del;341del;342del;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG]",
];

// ===========================================================================
// Outcome classification
// ===========================================================================

/// How a single `rederive` / `from_sequences` attempt turned out, for the
/// RED/GREEN maps.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Outcome {
    /// Rederived to a description (GREEN). Carries the derived string.
    Ok(String),
    /// The #2201 round-trip failure (RED): "could not be re-applied to its
    /// window".
    RoundTrip(String),
    /// Any other error (parse, coordinate, refusal). Not the #2201 shape.
    Other(String),
}

impl Outcome {
    /// Classify a `Result<HgvsVariant, FerroError>`-shaped outcome from its
    /// `Ok` string or error message.
    fn classify(result: Result<String, String>) -> Self {
        match result {
            Ok(s) => Outcome::Ok(s),
            Err(m) if m.contains("could not be re-applied to its window") => Outcome::RoundTrip(m),
            Err(m) => Outcome::Other(m),
        }
    }
    fn is_ok(&self) -> bool {
        matches!(self, Outcome::Ok(_))
    }
    fn is_collision(&self) -> bool {
        matches!(self, Outcome::RoundTrip(_))
    }
    fn tag(&self) -> &'static str {
        match self {
            Outcome::Ok(_) => "GREEN ok",
            Outcome::RoundTrip(_) => "RED   roundtrip",
            Outcome::Other(_) => "?     other-err",
        }
    }
    fn detail(&self) -> &str {
        match self {
            Outcome::Ok(s) | Outcome::RoundTrip(s) | Outcome::Other(s) => s,
        }
    }
}

// ===========================================================================
// Helpers
// ===========================================================================

/// A genome-capable provider over one contig named `NC_TEST.1`.
///
/// Genomic name on purpose: `from_sequences` (which `rederive` calls) emits
/// `g.` and refuses a transcript/protein accession.
fn provider(seq: &str) -> JsonProvider {
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { "NC_TEST.1": seq },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    JsonProvider::from_json(f.path()).unwrap()
}

/// 1-based inclusive slice of the TEMPLATE, as an owned string.
fn tmpl(lo1: u64, hi1: u64) -> String {
    TEMPLATE[(lo1 as usize - 1)..(hi1 as usize)].to_string()
}

/// Non-homologous filler of length `n` (a repeating tetramer).
fn filler(n: usize) -> String {
    const F: &str = "GACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACT";
    F[..n].to_string()
}

/// Apply `dels` (1-based positions) and one insertion after `ins_after`
/// (1-based) to `reference`, returning the alternate string.
fn splice(reference: &str, dels: &[u64], ins_after: u64, payload: &str) -> String {
    debug_assert!(
        !dels.contains(&ins_after),
        "ins_after {ins_after} is a deleted position; `splice` skips deleted positions \
         before the anchor check, so the payload would be silently dropped and the case \
         would degenerate to a plain deletion"
    );
    let mut alt = String::new();
    for (i, b) in reference.bytes().enumerate() {
        let pos = i as u64 + 1;
        if dels.contains(&pos) {
            continue;
        }
        alt.push(b as char);
        if pos == ins_after {
            alt.push_str(payload);
        }
    }
    alt
}

/// Build a description: a run of `k` base-less deletions at `start..=start+k-1`,
/// then an insertion of `payload` at `(start+k+gap)_(start+k+gap+1)`.
///
/// `gap == 0` is the reported adjacent shape: the insertion's left anchor is the
/// first base past the deletion run. For the tightest reported case #6
/// (`start=335, k=8`) this yields dels `335..=342` and `343_344ins…`, matching
/// the issue byte-for-byte. Each unit of `gap` moves the insertion one base
/// further 3', past that many retained reference bases.
fn del_run_then_insertion(start: u64, k: u64, gap: u64, payload: &str) -> String {
    let mut members: Vec<String> = (0..k).map(|i| format!("{}del", start + i)).collect();
    let ins_lo = start + k + gap; // first base past the run (+ gap retained bases)
    members.push(format!("{ins_lo}_{}ins{payload}", ins_lo + 1));
    format!("NC_TEST.1:g.[{}]", members.join(";"))
}

/// Run `rederive` on one description under an explicit shuffle direction and
/// classify the result.
fn rederive_dir(
    nz: &Normalizer<JsonProvider>,
    desc: &str,
    recommended_form: bool,
    dir: ShuffleDirection,
) -> Outcome {
    let variant = match parse_hgvs(desc) {
        Ok(v) => v,
        Err(e) => return Outcome::Other(format!("parse: {e}")),
    };
    let options = FromSequencesOptions::default().with_direction(dir);
    Outcome::classify(
        nz.rederive(&variant, &options, recommended_form)
            .map(|r| r.to_string())
            .map_err(|e| e.to_string()),
    )
}

/// [`rederive_dir`] with the default (3') shuffle direction.
fn rederive(nz: &Normalizer<JsonProvider>, desc: &str, recommended_form: bool) -> Outcome {
    rederive_dir(nz, desc, recommended_form, ShuffleDirection::ThreePrime)
}

/// Classify a derivation straight from reference/alternate bases through the
/// public `from_sequences` — the cheap substrate (no provider, no HGVS parsing).
/// It runs `verify_round_trip` internally, so an `Ok` is denoted-correct by
/// construction and a `RoundTrip` is the #2201 collision.
fn window_outcome(
    position: u64,
    reference: &str,
    alternate: &str,
    dir: ShuffleDirection,
) -> Outcome {
    let opts = FromSequencesOptions::default().with_direction(dir);
    Outcome::classify(
        ferro_hgvs::from_sequences("NC_TEST.1", position, reference, alternate, &opts)
            .map(|v| v.to_string())
            .map_err(|e| e.to_string()),
    )
}

/// Build `(reference_window, alternate_window)` from a TEMPLATE window
/// `[win_lo..=win_hi]` (1-based): delete the 1-based positions in `dels`, and
/// insert `payload` immediately after 1-based `ins_after`. `position` for
/// `from_sequences` is `win_lo`.
fn window_ref_alt(
    win_lo: u64,
    win_hi: u64,
    dels: &[u64],
    ins_after: u64,
    payload: &str,
) -> (String, String) {
    let reference = tmpl(win_lo, win_hi);
    let rel_dels: Vec<u64> = dels.iter().map(|d| d - win_lo + 1).collect();
    let alt = splice(&reference, &rel_dels, ins_after - win_lo + 1, payload);
    (reference, alt)
}

/// A tiny deterministic xorshift64 PRNG. Deterministic on purpose: `Math.random`
/// equivalents are banned in this repo's harnessed contexts, and a fixed seed
/// makes any counterexample the fuzzer finds reproducible.
struct Rng(u64);
impl Rng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4)]
    }
    fn dna(&mut self, n: usize) -> String {
        (0..n).map(|_| self.base() as char).collect()
    }
}

// ===========================================================================
// Part A — the exact reported geometry (must round-trip)
// ===========================================================================

/// The tightest reported case (#6): eight adjacent deletions `335..=342`
/// immediately followed by a 30-base insertion at `343_344`. The minimal
/// reported reproducer.
#[test]
fn tightest_del_run_then_big_insertion_rederives() {
    let nz = Normalizer::new(provider(TEMPLATE));
    let desc = "NC_TEST.1:g.[335delT;336delT;337delT;338delC;339delC;340delT;341delT;342delC;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG]";
    let got = rederive(&nz, desc, false);
    // Pin the exact derived members, not merely that they re-apply.
    // `merge_coincident_insertions` concatenates coincident payloads in list
    // order, and that order is trusted rather than verified (see the module doc),
    // with `verify_round_trip` the only backstop — so a wrong concatenation that
    // still re-applied would pass an `is_ok()` check. Captured from the current
    // run; if a deliberate derivation change moves it, update the string.
    let Outcome::Ok(derived) = got else {
        panic!("{}: {}", got.tag(), got.detail());
    };
    assert_eq!(
        derived, "NC_TEST.1:g.[336delinsCCCGAACAAAGAAA;338_339insA;340T>G;343_344insCTCGGAGG]",
        "the #2201 derivation for the tightest reported case is pinned"
    );
}

/// The same geometry with base-less `del`s (the other spelling half the reported
/// cases use) rederives identically — the deleted-base annotation is not what
/// triggers the collision.
#[test]
fn tightest_del_run_baseless_dels_rederives() {
    let nz = Normalizer::new(provider(TEMPLATE));
    let desc = del_run_then_insertion(335, 8, 0, PAYLOAD30);
    let got = rederive(&nz, &desc, false);
    assert!(got.is_ok(), "{}: {}", got.tag(), got.detail());
}

/// All eleven reported variants must rederive, under both `recommended_form`
/// values (the issue reports both fail). Prints a per-case RED/GREEN map.
#[test]
fn all_eleven_reported_variants_rederive() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== Part A: eleven reported variants ===");
    let mut failures = 0;
    for (i, v) in REPORTED_VARIANTS.iter().enumerate() {
        for form in [false, true] {
            let got = rederive(&nz, v, form);
            failures += usize::from(!got.is_ok());
            eprintln!(
                "  #{:<2} recommended_form={:<5} {}\n      -> {}",
                i + 1,
                form,
                got.tag(),
                got.detail()
            );
        }
    }
    assert_eq!(
        failures, 0,
        "{failures} of 22 (11 variants × 2 forms) failed to rederive"
    );
}

/// The reversed member order (#8): a net-positive insertion FIRST, then the
/// deletion run 3' of it. Same collision, different arrival order.
#[test]
fn insertion_then_del_run_rederives() {
    let nz = Normalizer::new(provider(TEMPLATE));
    let desc = "NC_TEST.1:g.[334_335insCCCCGAACAAAGAAATCACGTCTCTCGGAGG;336del;337del;338del;339del;340del;341del;342del;343del;344del]";
    let got = rederive(&nz, desc, false);
    assert!(got.is_ok(), "{}: {}", got.tag(), got.detail());
}

/// The short-payload reported case (#7), at a different locus with a 4-base
/// insertion — confirms the family is not specific to the 335 hotspot or a
/// 30-mer payload.
#[test]
fn short_payload_distant_locus_rederives() {
    let nz = Normalizer::new(provider(TEMPLATE));
    let desc = "NC_TEST.1:g.[1123_1124insACAC;1124C>A;1129T>A;1223A>C]";
    let got = rederive(&nz, desc, false);
    assert!(got.is_ok(), "{}: {}", got.tag(), got.detail());
}

// ===========================================================================
// Part B — 1-D sweeps mapping the RED/GREEN boundary
// ===========================================================================

/// Sweep the deletion-run length with the payload fixed.
#[test]
fn sweep_deletion_run_length() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== Part B1: deletion-run length (start=335, gap=0, payload=30mer) ===");
    let mut failures = 0;
    for k in 1..=12u64 {
        let got = rederive(&nz, &del_run_then_insertion(335, k, 0, PAYLOAD30), false);
        failures += usize::from(!got.is_ok());
        eprintln!("  k={:<2} {}  {}", k, got.tag(), got.detail());
    }
    assert_eq!(
        failures, 0,
        "{failures}/12 deletion-run lengths failed to rederive"
    );
}

/// Sweep the insertion payload length (prefixes of `PAYLOAD30`) with an
/// 8-deletion run fixed.
#[test]
fn sweep_insertion_payload_length() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== Part B2: payload length (start=335, k=8 dels, gap=0) ===");
    let mut failures = 0;
    for len in 1..=PAYLOAD30.len() {
        let got = rederive(
            &nz,
            &del_run_then_insertion(335, 8, 0, &PAYLOAD30[..len]),
            false,
        );
        failures += usize::from(!got.is_ok());
        eprintln!("  payload_len={:<2} {}  {}", len, got.tag(), got.detail());
    }
    assert_eq!(
        failures,
        0,
        "{failures}/{} payload lengths failed to rederive",
        PAYLOAD30.len()
    );
}

/// Sweep the gap (retained bases) between the deletion run and the insertion.
/// `gap==0` is the reported adjacent shape.
#[test]
fn sweep_gap_between_run_and_insertion() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== Part B3: gap between run and insertion (start=335, k=8, payload=30mer) ===");
    let mut failures = 0;
    for gap in 0..=5u64 {
        let got = rederive(&nz, &del_run_then_insertion(335, 8, gap, PAYLOAD30), false);
        failures += usize::from(!got.is_ok());
        eprintln!("  gap={:<2} {}  {}", gap, got.tag(), got.detail());
    }
    assert_eq!(failures, 0, "{failures}/6 gaps failed to rederive");
}

/// Sweep the payload CONTENT at fixed length — discriminates "any net-positive
/// insertion beside a run" from "only payloads that induce a dup".
#[test]
fn sweep_payload_content() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== Part B4: payload content (start=335, k=8, gap=0, len=30) ===");
    let mut failures = 0;
    for (name, payload) in [
        ("reported", PAYLOAD30),
        ("all-A", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        ("all-G", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
        ("alt-AC", "ACACACACACACACACACACACACACACAC"),
        ("alt-GT", "GTGTGTGTGTGTGTGTGTGTGTGTGTGTGT"),
    ] {
        assert_eq!(payload.len(), 30, "payload {name} must be 30 bases");
        let got = rederive(&nz, &del_run_then_insertion(335, 8, 0, payload), false);
        failures += usize::from(!got.is_ok());
        eprintln!("  {:<9} {}  {}", name, got.tag(), got.detail());
    }
    assert_eq!(
        failures, 0,
        "{failures}/5 payload contents failed to rederive"
    );
}

/// Sweep the deletion-run start position across the contig — tells whether the
/// collision is a property of the 335 locus's bases or of the geometry itself.
#[test]
fn sweep_run_start_position() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== Part B5: run start position (k=8, gap=0, payload=30mer) ===");
    let mut failures = 0;
    for start in [100u64, 200, 335, 500, 700, 900, 1100] {
        let got = rederive(&nz, &del_run_then_insertion(start, 8, 0, PAYLOAD30), false);
        failures += usize::from(!got.is_ok());
        eprintln!("  start={:<5} {}  {}", start, got.tag(), got.detail());
    }
    assert_eq!(
        failures, 0,
        "{failures}/7 start positions failed to rederive"
    );
}

/// Controls: the deletion run alone and the insertion alone must both rederive.
/// If either is RED, the collision is not specific to the run+insertion
/// combination and the diagnosis is wrong.
#[test]
fn controls_run_alone_and_insertion_alone_rederive() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== Part B6: controls (isolate the combination) ===");
    let run_only =
        "NC_TEST.1:g.[335del;336del;337del;338del;339del;340del;341del;342del]".to_string();
    let ins_only = format!("NC_TEST.1:g.343_344ins{PAYLOAD30}");
    let mut failures = 0;
    for desc in [run_only, ins_only] {
        let got = rederive(&nz, &desc, false);
        failures += usize::from(!got.is_ok());
        eprintln!("  {}  {}\n    {}", got.tag(), desc, got.detail());
    }
    assert_eq!(failures, 0, "a control geometry failed to rederive");
}

// ===========================================================================
// Tier 1 — the trigger reverse-engineered, plus grid / directions / idempotency
// ===========================================================================

/// The trigger is a payload-*tail* × reference-flank coincidence, and it is
/// PREFIX-INDEPENDENT.
///
/// Two hypotheses were tried and *refuted* first (kept here, since a refuted
/// hypothesis recurs):
///   1. "payload copies its adjacent flank bases" — constructing payload = `h`
///      copied flank bases + filler, swept h=0..12, produced **0** collisions.
///   2. a fully-synthetic reconstruction of the hotspot's immediate 5'/3' flanks
///      (`…GAAACAC[core]TCTAAGCAC…`) also does not collide — so the trigger needs
///      the wider TEMPLATE context, not just the local window.
///
/// What DOES reproduce, reverse-engineered with `from_sequences` as an oracle: at
/// the hotspot, **any** filler prefix followed by `TRIGGER_TAIL` collides — the
/// prefix is arbitrary (it just becomes the leading `delins`). This sweeps five
/// unrelated prefixes to show prefix-independence, and sweeps the tail length to
/// find the ~15-base threshold.
#[test]
fn t1a_trigger_tail_is_prefix_independent() {
    let nz = Normalizer::new(provider(TEMPLATE));
    let mut failures = 0;

    eprintln!("\n=== T1a: prefix-independence (14-base prefix + TRIGGER_TAIL, at hotspot) ===");
    for prefix in [
        "AAAAAAAAAAAAAA",
        "CCCCCCCCCCCCCC",
        "GGGGGGGGGGGGGG",
        "GACTGACTGACTGA",
        "CGCGCGCGCGCGCG",
    ] {
        let payload = format!("{prefix}{TRIGGER_TAIL}");
        let got = rederive(&nz, &del_run_then_insertion(335, 8, 0, &payload), false);
        failures += usize::from(!got.is_ok());
        eprintln!("  prefix={prefix}  {}  {}", got.tag(), got.detail());
    }

    eprintln!("\n=== T1a: tail-length threshold (A-prefix padded to total 30) ===");
    // Suffixes of TRIGGER_TAIL, longest first; the collision needs ~15 bases of
    // tail. Only the full-length tail is a hard failure (that is the reproducer);
    // the shorter cells are characterized, not asserted.
    for taillen in (8..=TRIGGER_TAIL.len()).rev() {
        let tail = &TRIGGER_TAIL[TRIGGER_TAIL.len() - taillen..];
        let payload = format!("{}{tail}", "A".repeat(30 - taillen));
        let got = rederive(&nz, &del_run_then_insertion(335, 8, 0, &payload), false);
        if taillen == TRIGGER_TAIL.len() {
            failures += usize::from(!got.is_ok());
        }
        eprintln!("  taillen={:<2} {}  {}", taillen, got.tag(), got.detail());
    }

    assert_eq!(
        failures, 0,
        "{failures} prefix-independent / full-tail cells failed to rederive"
    );
}

/// The k × payload-length 2-D grid. The 1-D sweeps each pass through one corner;
/// a red region defined by an interaction of both axes is a curve those slices
/// can only touch. Prints a compact grid so the shape is visible.
#[test]
fn t1b_run_length_by_payload_length_grid() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== T1b: k (rows) × payload_len (cols) grid; '.'=GREEN 'X'=RED '?'=other ===");
    let plens = [6usize, 12, 18, 24, 30];
    eprint!("   k\\len ");
    for p in plens {
        eprint!("{p:>4}");
    }
    eprintln!();
    let mut failures = 0;
    for k in 1..=12u64 {
        eprint!("   {k:>5}  ");
        for &plen in &plens {
            let got = rederive(
                &nz,
                &del_run_then_insertion(335, k, 0, &PAYLOAD30[..plen]),
                false,
            );
            failures += usize::from(got.is_collision());
            eprint!(
                "{:>4}",
                match got {
                    Outcome::Ok(_) => '.',
                    Outcome::RoundTrip(_) => 'X',
                    Outcome::Other(_) => '?',
                }
            );
        }
        eprintln!();
    }
    assert_eq!(failures, 0, "{failures} grid cells hit the collision");
}

/// Repeat context, two findings that must be kept apart.
///
/// **(a) The trigger SURVIVES a tandem-repeat prefix (would collide without the
/// fix).** Since the prefix is arbitrary (see `t1a`), a `CAG×n` or `AC×n` repeat prefix
/// followed by `TRIGGER_TAIL` collides at the hotspot exactly as a homopolymer
/// prefix does — so an indel whose inserted sequence begins with a short tandem
/// expansion and ends in the trigger tail is in the failing family.
///
/// **(b) A pure tandem-repeat REFERENCE resolves cleanly (control, GREEN).** A
/// clean STR is not the "natural habitat" of this bug: deleting whole units of a
/// `CAG` tract and inserting `n` extra copies always resolves to a single clean
/// `dup`, never the collision — the trigger needs the asymmetric non-repeat 3'
/// context a uniform tract does not provide. Kept as the evidence that clean STR
/// contexts are safe, so the (a)/(b) distinction is not re-litigated.
///
/// Only the (a) cases are asserted to rederive — they are the ones the fix
/// rescues; without it they collide.
#[test]
fn t1c_repeat_context() {
    eprintln!("\n=== T1c(a): repeat-unit prefix + TRIGGER_TAIL at hotspot (reproduces) ===");
    let nz = Normalizer::new(provider(TEMPLATE));
    let mut failures = 0;
    for (unit, n) in [("CAG", 5usize), ("CAG", 8), ("AC", 8), ("AC", 11)] {
        let payload = format!("{}{TRIGGER_TAIL}", unit.repeat(n));
        let got = rederive(&nz, &del_run_then_insertion(335, 8, 0, &payload), false);
        failures += usize::from(!got.is_ok());
        eprintln!("  ({unit})x{n} + tail  {}  {}", got.tag(), got.detail());
    }

    eprintln!(
        "\n=== T1c(b): pure CAG-tract reference — clean tandem, does NOT collide (control) ==="
    );
    let unit = "CAG";
    let seq = format!("{}{}{}", filler(30), unit.repeat(20), filler(30));
    let nz_str = Normalizer::new(provider(&seq));
    let tract_lo = 31u64; // tract occupies 31 ..= 30 + 3*20
    let mut control_collisions = 0;
    for n in 1..=8u64 {
        let (del_lo, del_hi) = (tract_lo, tract_lo + 5); // two 3-base units
        let payload = unit.repeat(n as usize);
        let desc = format!(
            "NC_TEST.1:g.[{del_lo}_{del_hi}del;{del_hi}_{}ins{payload}]",
            del_hi + 1
        );
        let got = rederive(&nz_str, &desc, false);
        control_collisions += usize::from(got.is_collision());
        eprintln!("  n={:<2} {}  {}", n, got.tag(), got.detail());
    }

    assert_eq!(
        control_collisions, 0,
        "control: a clean STR case unexpectedly collided"
    );
    assert_eq!(
        failures, 0,
        "{failures}/4 repeat-prefix reproducers failed to rederive"
    );
}

/// Both shuffle directions across all eleven reported variants. The derived
/// collision may sit at a different junction (or vanish) under 5' shuffling; the
/// issue reports failure regardless, so both must round-trip.
#[test]
fn t1d_both_shuffle_directions() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== T1d: both shuffle directions ===");
    let mut failures = 0;
    for (i, v) in REPORTED_VARIANTS.iter().enumerate() {
        for dir in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            let got = rederive_dir(&nz, v, false, dir);
            failures += usize::from(!got.is_ok());
            eprintln!("  #{:<2} {:?}  {}", i + 1, dir, got.tag());
        }
    }
    assert_eq!(
        failures, 0,
        "{failures} of 22 (11 × 2 directions) failed to rederive"
    );
}

/// Idempotency over the run-length sweep. For every cell that rederives,
/// `rederive` of its own output must be a fixed point — checked on the GREEN
/// cells too, since a cell next to the red region can round-trip once yet not be
/// stable. RED cells cannot be checked and are counted as failures.
#[test]
fn t1e_idempotency_over_run_length_sweep() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== T1e: idempotency (rederive(rederive(x)) == rederive(x)) ===");
    let mut failures = 0;
    for k in 1..=12u64 {
        match rederive(&nz, &del_run_then_insertion(335, k, 0, PAYLOAD30), false) {
            Outcome::Ok(once) => match rederive(&nz, &once, false) {
                Outcome::Ok(twice) if twice == once => {
                    eprintln!("  k={:<2} GREEN fixed-point  {once}", k)
                }
                other => {
                    failures += 1;
                    eprintln!(
                        "  k={:<2} RED   not-a-fixed-point\n    once:  {once}\n    twice: {}",
                        k,
                        other.detail()
                    );
                }
            },
            not_ok => {
                failures += 1;
                eprintln!(
                    "  k={:<2} RED   first rederive failed: {}",
                    k,
                    not_ok.detail()
                );
            }
        }
    }
    assert_eq!(
        failures, 0,
        "{failures}/12 run-length cells were not idempotent fixed points"
    );
}

// ===========================================================================
// Tier 2 — confluence, net-sign symmetry, window stability, structured corpus
// ===========================================================================

/// Confluence for the tightest geometry. Three independent spellings of ONE
/// underlying variant (member-wise del-run + ins, a range-delete + ins, and a
/// single spanning delins) must all rederive to the same canonical string.
/// Convergence is a stronger property than "no error".
#[test]
fn t2a_confluence_of_spellings() {
    let nz = Normalizer::new(provider(TEMPLATE));
    // Region 335..=344 is ref "TTTCCTTCTT"; the variant deletes 335..=342, keeps
    // 343 (T) and 344 (T), and inserts PAYLOAD30 after 343 — so the spanning
    // delins over 335..=344 replaces those ten bases with "T"+PAYLOAD30+"T".
    let spellings = [
        del_run_then_insertion(335, 8, 0, PAYLOAD30),
        format!("NC_TEST.1:g.[335_342del;343_344ins{PAYLOAD30}]"),
        format!("NC_TEST.1:g.335_344delinsT{PAYLOAD30}T"),
    ];
    eprintln!("\n=== T2a: confluence of three spellings ===");
    let outs: Vec<Outcome> = spellings
        .iter()
        .map(|s| {
            let got = rederive(&nz, s, false);
            eprintln!("  {}  {}\n    -> {}", got.tag(), s, got.detail());
            got
        })
        .collect();
    assert!(
        outs.iter().all(Outcome::is_ok),
        "a spelling failed to rederive"
    );
    let first = outs[0].detail();
    assert!(
        outs.iter().all(|o| o.detail() == first),
        "spellings did not converge on one canonical form"
    );
}

/// Net-sign symmetry. The reported geometry is net-positive (del-run + bigger
/// insertion). Explore the mirror (single-base insertions + a bigger deletion,
/// net-negative) and the net-zero case, to see whether the defect family is
/// imbalance-direction-symmetric.
#[test]
fn t2b_net_sign_symmetry() {
    let nz = Normalizer::new(provider(TEMPLATE));
    eprintln!("\n=== T2b: net-sign symmetry ===");
    // Net-negative mirror: single-base insertions at 335_336..342_343, then a
    // larger deletion 344..361.
    let mut mirror: Vec<String> = (335..=342u64)
        .map(|p| format!("{p}_{}insA", p + 1))
        .collect();
    mirror.push("344_361del".to_string());
    let net_negative = format!("NC_TEST.1:g.[{}]", mirror.join(";"));
    // Net-zero: eight deletions + an eight-base insertion.
    let net_zero = del_run_then_insertion(335, 8, 0, &PAYLOAD30[..8]);

    let mut failures = 0;
    for (label, desc) in [("net-negative", &net_negative), ("net-zero", &net_zero)] {
        let got = rederive(&nz, desc, false);
        failures += usize::from(!got.is_ok());
        eprintln!("  {:<12} {}  {}", label, got.tag(), got.detail());
    }
    assert_eq!(
        failures, 0,
        "{failures}/2 net-sign geometries failed to rederive"
    );
}

/// Window-perturbation metamorphism. Derive the tightest geometry from windows
/// of increasing flank around the hotspot; the canonical output must be
/// derivable and identical regardless of where the window boundary falls. A
/// derivation that depends on the window boundary is a second defect class the
/// guard-trip alone cannot see.
#[test]
fn t2c_window_perturbation_stability() {
    eprintln!("\n=== T2c: window-perturbation stability (hotspot 335..344) ===");
    let dels: Vec<u64> = (335..=342).collect();
    let outs: Vec<Outcome> = [12u64, 20, 30, 45, 60]
        .into_iter()
        .map(|pad| {
            let (win_lo, win_hi) = (335 - pad, 344 + pad);
            let (reference, alternate) = window_ref_alt(win_lo, win_hi, &dels, 343, PAYLOAD30);
            let got = window_outcome(win_lo, &reference, &alternate, ShuffleDirection::ThreePrime);
            eprintln!("  pad={:<3} {}  {}", pad, got.tag(), got.detail());
            got
        })
        .collect();
    assert!(outs.iter().all(Outcome::is_ok), "a window failed to derive");
    let first = outs[0].detail();
    assert!(
        outs.iter().all(|o| o.detail() == first),
        "derived canonical form was not stable across window widths"
    );
}

/// Structured window-level corpus. Not exhaustive over all DNA (that space is
/// astronomical); a bounded product of run lengths × payload lengths × two gaps
/// at the hotspot, driven directly through `from_sequences`. Asserts the
/// invariant the round-trip guard enforces: no derivation may collide. Reports
/// how many cells were exercised so a zero is not mistaken for an empty corpus.
#[test]
fn t2d_structured_window_corpus() {
    eprintln!("\n=== T2d: structured window corpus (from_sequences) ===");
    let (win_lo, win_hi) = HOTSPOT_WINDOW;
    let (mut collisions, mut others, mut exercised) = (0, 0, 0);
    for k in 1..=12u64 {
        let dels: Vec<u64> = (335..335 + k).collect();
        for plen in [4usize, 8, 12, 16, 20, 24, 28, 30] {
            for gap in 0..=2u64 {
                // first base past the deletion run (+ gap retained bases), matching
                // `del_run_then_insertion`'s `start + k + gap` anchor
                let ins_after = 335 + k + gap;
                let (reference, alternate) =
                    window_ref_alt(win_lo, win_hi, &dels, ins_after, &PAYLOAD30[..plen]);
                exercised += 1;
                // The contract is that every structured cell DERIVES: the
                // collision is one way to fail it, a refusal / WouldNotRender /
                // coordinate error (`Outcome::Other`) is another, and counting
                // only collisions would let the latter pass silently.
                let outcome =
                    window_outcome(win_lo, &reference, &alternate, ShuffleDirection::ThreePrime);
                if outcome.is_collision() {
                    collisions += 1;
                } else if !outcome.is_ok() {
                    others += 1;
                }
            }
        }
    }
    eprintln!("  exercised {exercised} cells, {collisions} collisions, {others} other errors");
    assert_eq!(
        (collisions, others),
        (0, 0),
        "of {exercised} structured-corpus cells, {collisions} collided and {others} otherwise failed to derive"
    );
}

// ===========================================================================
// Tier 3 — a fuzzer proven able to fail, and an opt-in corpus replay
// ===========================================================================

/// Neighbourhood fuzzer, seeded from the known trigger.
///
/// A uniform-random / flank-homology fuzzer is BLIND here — an earlier version
/// found **0 collisions in 4000 cases**, which refuted the "accidental flank
/// homology" hypothesis rather than proving safety. Per this repo's doctrine a
/// guard must be *shown able to fail*, so this fuzzer explores the NEIGHBOURHOOD
/// of the reported reproducer: it starts from the hotspot with `PAYLOAD30` (known
/// to collide) and applies small random perturbations — payload point-mutations,
/// ±1 length, run length 7–9, gap 0–1. The base case collides without the fix, so
/// this test fails without it and passes with it, and the perturbations map how
/// wide the colliding region was. The pure-random search is re-run and its
/// (expected zero) rate reported, so the blindness is recorded, not hidden.
#[test]
fn t3a_neighborhood_fuzzer_seeded_from_the_trigger() {
    const CASES: usize = 4000;
    let (win_lo, win_hi) = HOTSPOT_WINDOW;

    // --- pure-random control: recorded as a BLIND result, not coverage. ---
    let mut rng = Rng(0x2201_c0ffee);
    let mut random_collisions = 0;
    for _ in 0..CASES {
        let n = 80 + rng.below(60);
        let reference = rng.dna(n);
        let k = 1 + rng.below(12) as u64;
        let start = 21 + rng.below(n - 60) as u64;
        let gap = rng.below(3) as u64;
        let payload_len = 4 + rng.below(28);
        let payload = rng.dna(payload_len);
        let dels: Vec<u64> = (start..start + k).collect();
        let alt = splice(&reference, &dels, start + k + gap, &payload);
        random_collisions += usize::from(
            window_outcome(1, &reference, &alt, ShuffleDirection::ThreePrime).is_collision(),
        );
    }

    // --- neighbourhood search around the known reproducer. ---
    let mut rng = Rng(0xdead_2201);
    let mut collisions = 0;
    let mut others = 0;
    let mut smallest: Option<String> = None;
    for case in 0..CASES {
        // case 0 is the exact reported reproducer; the rest perturb it.
        let (k, gap) = if case == 0 {
            (8, 0)
        } else {
            (7 + rng.below(3) as u64, rng.below(2) as u64)
        };
        let mut payload: Vec<u8> = PAYLOAD30.bytes().collect();
        if case != 0 {
            match rng.below(3) {
                0 if payload.len() > 20 => {
                    payload.pop();
                }
                1 => payload.push(rng.base()),
                _ => {}
            }
            for _ in 0..rng.below(4) {
                let i = rng.below(payload.len());
                payload[i] = rng.base();
            }
        }
        let payload = String::from_utf8(payload).unwrap();
        let dels: Vec<u64> = (335..335 + k).collect();
        let (reference, alternate) = window_ref_alt(win_lo, win_hi, &dels, 335 + k + gap, &payload);
        let outcome = window_outcome(win_lo, &reference, &alternate, ShuffleDirection::ThreePrime);
        if outcome.is_collision() {
            collisions += 1;
            if smallest.as_ref().is_none_or(|s| payload.len() < s.len()) {
                smallest = Some(payload);
            }
        } else if !outcome.is_ok() {
            // A neighbourhood cell that fails to derive for any other reason is
            // still a failure of the "the whole neighbourhood round-trips"
            // contract, not just a collision.
            others += 1;
        }
    }

    eprintln!("\n=== T3a: neighbourhood fuzzer ===");
    eprintln!("  pure-random control: {random_collisions}/{CASES} collisions (blind — expected 0)");
    eprintln!("  neighbourhood: {collisions}/{CASES} collisions, {others} other errors");
    if let Some(payload) = &smallest {
        eprintln!(
            "  smallest colliding payload (len={}): {payload}",
            payload.len()
        );
    }
    // Seeded with a known collider, so a green run means the whole neighbourhood
    // round-trips — the post-fix state. A cell that fails any other way
    // (`Outcome::Other`) breaks that contract too, so both counts must be zero.
    assert_eq!(
        (collisions, others),
        (0, 0),
        "neighbourhood fuzzer found {collisions}/{CASES} colliding and {others}/{CASES} otherwise-failing derivations"
    );
}

/// Opt-in corpus replay, the standing guard for the NEXT occurrence at real-data
/// scale. Set `FERRO_2201_REPLAY` to a file of newline-delimited `NC_TEST.1:g.…`
/// descriptions (against the synthetic TEMPLATE) to replay them through
/// `rederive` and assert none collide. Skips (green) when unset — the large
/// real-world indel corpora (ClinVar/CMRG/Paraphase) require the prepared
/// reference and are exercised in the manifest-backed nightly, not here.
#[test]
fn t3b_optin_corpus_replay() {
    let Ok(path) = std::env::var("FERRO_2201_REPLAY") else {
        eprintln!("\n=== T3b: skipped (FERRO_2201_REPLAY unset) ===");
        return;
    };
    let body = std::fs::read_to_string(&path).expect("read FERRO_2201_REPLAY corpus");
    let nz = Normalizer::new(provider(TEMPLATE));
    let (mut collisions, mut checked) = (0, 0);
    for line in body.lines().map(str::trim).filter(|l| !l.is_empty()) {
        checked += 1;
        if let Outcome::RoundTrip(m) = rederive(&nz, line, false) {
            collisions += 1;
            eprintln!("  RED {line}\n    -> {m}");
        }
    }
    eprintln!("\n=== T3b: replayed {checked} descriptions from {path} ===");
    assert!(checked > 0, "FERRO_2201_REPLAY corpus was empty");
    assert_eq!(
        collisions, 0,
        "{collisions}/{checked} replayed descriptions collided"
    );
}

// ===========================================================================
// Tier 4 — #2206 shape coverage
//
// `verify_round_trip` (src/normalize/from_sequences.rs) runs on EVERY derivation
// and guarantees the derived BASES equal the alternate — so an `Ok` is
// denoted-correct and `merge_coincident_insertions`' concatenation order is
// proven every run. What it cannot see is SPELLING/SHAPE: two spellings that
// denote the same bases are indistinguishable to it. These pin the SHAPE of the
// merge across the regression's sibling geometries, which the del-run+ins tests
// above never exercise. Strings captured from the current run; a deliberate
// derivation change moves them.
// ===========================================================================

/// Adjacent-but-not-coincident insertion anchors are NOT over-merged (#2206).
///
/// `merge_coincident_insertions` folds two pure insertions sharing ONE interbase.
/// The near-miss it must leave alone is two insertions one base apart:
/// `del_run_then_insertion(335, 8, gap)` with `gap >= 1` anchors the insertion
/// past the run, and the derivation places insertions at adjacent but DISTINCT
/// interbases (e.g. `339_340ins…;340_341ins…`). They must stay separate — a fold
/// keyed on adjacency rather than coincidence would wrongly merge them.
#[test]
fn t4a_near_miss_adjacent_anchors_are_not_over_merged() {
    let nz = Normalizer::new(provider(TEMPLATE));
    // gap == 1: the insertion sits one base 3' of the deletion run.
    let got = rederive(&nz, &del_run_then_insertion(335, 8, 1, PAYLOAD30), false);
    let Outcome::Ok(derived) = got else {
        panic!(
            "gap=1 near-miss must derive: {}: {}",
            got.tag(),
            got.detail()
        );
    };
    assert_eq!(
        derived,
        "NC_TEST.1:g.[337T>C;339_340insGAACAAAGAAA;340_341insCACG;343_344insC;345_346insGGAGGC]",
        "adjacent-but-not-coincident anchors must stay separate, not merge"
    );
    // gap == 2: one base further; still distinct anchors, still not merged.
    let got = rederive(&nz, &del_run_then_insertion(335, 8, 2, PAYLOAD30), false);
    let Outcome::Ok(derived) = got else {
        panic!(
            "gap=2 near-miss must derive: {}: {}",
            got.tag(),
            got.detail()
        );
    };
    assert_eq!(
        derived,
        "NC_TEST.1:g.[337T>C;339_340insCGAACAAAGAAA;340_341insCACG;343_344insC;345_346insGGAGG]",
        "adjacent-but-not-coincident anchors must stay separate, not merge"
    );
}

/// A `dup` as the 5' member (not a deletion run) rederives (#2206).
///
/// The regression touched a `dup` OR `del` 5' member; the tests above cover only
/// del runs. Here the 5' member is a duplication, carried through the same
/// coalesce/merge path, and it survives the merge as a `dup`.
#[test]
fn t4b_dup_as_five_prime_member_rederives() {
    let nz = Normalizer::new(provider(TEMPLATE));
    let desc = "NC_TEST.1:g.[335_336dup;337del;338del;339del;340del;341del;342del;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG]";
    let got = rederive(&nz, desc, false);
    let Outcome::Ok(derived) = got else {
        panic!("dup-5' member must derive: {}: {}", got.tag(), got.detail());
    };
    assert_eq!(
        derived,
        "NC_TEST.1:g.[336_337dup;339_340insCGAACAAAGAAA;340_341insCACG;343_344insCTCGGAGG]",
        "a dup 5' member survives the merge as a dup"
    );
}

/// A deletion run followed by a `dup` (not an insertion) rederives (#2206).
///
/// The other sibling: a `dup` on the 3' side of the run. Both a single-base and a
/// two-base 3' dup collapse to one net deletion. A dup that OVERLAPS the run
/// (`342_343dup`) is a self-cancelling allele the parser rejects, so it cannot be
/// an input — noted so the absence is deliberate, not an oversight.
#[test]
fn t4c_del_run_then_dup_rederives() {
    let nz = Normalizer::new(provider(TEMPLATE));
    let run = "335del;336del;337del;338del;339del;340del;341del;342del";
    for (dup, expected) in [
        ("343dup", "NC_TEST.1:g.339_345del"),
        ("343_344dup", "NC_TEST.1:g.338_343del"),
    ] {
        let desc = format!("NC_TEST.1:g.[{run};{dup}]");
        let got = rederive(&nz, &desc, false);
        let Outcome::Ok(derived) = got else {
            panic!("{dup}: must derive: {}: {}", got.tag(), got.detail());
        };
        assert_eq!(derived, expected, "del-run then {dup}");
    }
}

/// The mitochondrial (`m.`) axis derives through the same merge path (#2206).
///
/// `merge_coincident_insertions` is axis-agnostic, and `from_sequences` emits
/// `m.` when the accession is mitochondrial (`NC_012920`). The tightest reported
/// geometry, on `NC_012920.1`, derives to the `m.`-prefixed counterpart of the
/// `g.` case #6 — same members, `m.` axis.
#[test]
fn t4d_mitochondrial_axis_rederives() {
    let doc = serde_json::json!({
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [],
        "genomic_sequences": { "NC_012920.1": TEMPLATE },
    });
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(doc.to_string().as_bytes()).unwrap();
    let nz = Normalizer::new(JsonProvider::from_json(f.path()).unwrap());
    let desc = "NC_012920.1:m.[335del;336del;337del;338del;339del;340del;341del;342del;343_344insCCCGAACAAAGAAATCACGTCTCTCGGAGG]";
    let variant = parse_hgvs(desc).expect("parse m. description");
    let derived = nz
        .rederive(&variant, &FromSequencesOptions::default(), false)
        .expect("m. axis must derive")
        .to_string();
    assert_eq!(
        derived, "NC_012920.1:m.[336delinsCCCGAACAAAGAAA;338_339insA;340T>G;343_344insCTCGGAGG]",
        "the m. axis derives the m.-prefixed counterpart of the tightest g. case"
    );
}
