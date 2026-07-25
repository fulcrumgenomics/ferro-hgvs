//! Property-based normalization idempotency test over reference-backed
//! synthetic fixtures (issue #1157 follow-up campaign).
//!
//! The pre-existing idempotency coverage is reference-free (empty provider, so
//! nothing shifts) or substitution-only (subs never shift). This fuzzes random
//! **repeat-biased** cores (homopolymers + short tandems, where del/dup/ins
//! actually 3'-shift) across `g.`, `c.`, `n.`, and `r.` — the last on a coding
//! *and* a non-coding transcript, which are two different axes — each
//! transcript axis on both strands, and every edit kind (`del`, `dup`, `ins`,
//! `delins`, `inv`, `con`), asserting `norm(norm(x)) == norm(x)`. `c.` likewise
//! appears twice: once where the transcript *is* the CDS, and once on a
//! UTR-bearing transcript, where alone a shift can cross `cds_end` into `*N`.
//!
//! `case_strategy_covers_every_system_and_edit_kind` is the non-vacuity guard:
//! it samples the strategy and fails if any (coordinate system x edit kind)
//! cell is empty, so the property above cannot quietly stop exercising an axis.
//!
//! Scope: **normalization only**. Nothing here goes through `VariantProjector`,
//! which reads the `r.` axis as transcript-relative rather than CDS-relative
//! (issue #1177) and would fail for that unrelated reason.
//!
//! Three properties are fuzzed:
//!   * `norm(norm(x)) == norm(x)` (idempotency), covering the full edit/coord
//!     matrix in the **3' (default, HGVS-canonical)** direction;
//!   * the same idempotency property in the **5'** direction — the mirror of the
//!     3' property, enabled once the campaign's pre-existing 5' bugs were fixed
//!     (see the 5' scope note below);
//!   * multiple spellings of ONE homopolymer-boundary haplotype normalize to ONE
//!     form (confluence), in **both** directions. Idempotency is necessary but
//!     NOT sufficient — a mode can be idempotent yet non-confluent (each spelling
//!     its own fixed point), which the confluence property catches. The 5'
//!     boundary/homopolymer shuffle was reworked to a coordinate-identity clamp
//!     keyed on the post-shuffle edit (see `src/normalize/mod.rs` CDS-start/-end
//!     clamps + the dup-routing recursion in `normalize_na_edit`) so those
//!     spellings now converge.
//!
//! 5' scope note: the general 5' idempotency fuzz is now enabled. Reaching green
//! required fixing three separate pre-existing 5' bugs the fuzzer surfaced, all
//! outside the boundary rework: the boundary dup/insertion→delins family, an
//! insertion→dup step that selected the wrong tract for heterogeneous content
//! (#7), a tandem-repeat `dup` whose `unit[N]` phase rotated on re-normalization
//! (#8, `normalize_repeat` now direction-aware), and the ins→dup sibling of #8 —
//! a single-copy dup over a dinucleotide tract that under-shifted by one when the
//! run continued off-phase 5' (fixed by mirroring the 3' branch's tract
//! alignment). Each has an explicit regression test alongside this fuzz.
//!
//! Cases are generated **valid by construction** (edit endpoints derived from
//! the core length), so there is no skip path at all: an unparseable render or a
//! failing first-pass normalize is a bug in the generator or the normalizer and
//! fails the test rather than silently shrinking the exercised set. Failure
//! seeds are recorded by proptest under
//! `tests/proptest-regressions/normalize_idempotency_proptest.txt` (proptest
//! anchors that directory at the crate's `tests/`, not at this module's dir).

use crate::common::synthetic::{padded, SyntheticBuilder, PAD_OFFSET};
use ferro_hgvs::reference::transcript::Strand;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};
use proptest::prelude::*;
use proptest::test_runner::{Config as ProptestConfig, TestCaseError};

fn normalizer_for(core: &str, sys: Sys, dir: ShuffleDirection) -> Normalizer<MockProvider> {
    Normalizer::with_config(
        provider_for(core, sys),
        NormalizeConfig::default().with_direction(dir),
    )
}

fn both_directions() -> impl Strategy<Value = ShuffleDirection> {
    prop::sample::select(vec![
        ShuffleDirection::ThreePrime,
        ShuffleDirection::FivePrime,
    ])
}

/// Coordinate system under test.
///
/// Every variant is arranged so that **core position `p` is rendered as
/// position `p`** on its own axis (genomic positions additionally carry
/// [`PAD_OFFSET`]), which is what lets one generator emit valid coordinates for
/// all of them from a single core length.
///
/// The `r.` axis appears twice on purpose, because its meaning depends on the
/// transcript it sits on (HGVS `background/numbering.md` L58/L61, issue #469,
/// pinned by `issue_291_rna_axis_convention`):
///
/// - on a **coding** transcript `r.` is **CDS-relative** — literally the same
///   axis as `c.`, so `r.N` and `c.N` name the same base;
/// - on a **non-coding** transcript there is no CDS, so `r.` coincides with
///   `n.` (transcript-relative).
///
/// The `r.` systems therefore share fixtures with their same-axis twins:
/// [`Sys::RnaNoncodingPlus`] uses the very `NR_TEST.1` provider
/// [`Sys::NoncodingPlus`] uses, and [`Sys::RnaCodingPlus`] uses a coding
/// `NM_TEST.1` provider whose `c.` spelling would be the same axis.
#[derive(Debug, Clone, Copy)]
enum Sys {
    /// `g.` on the padded contig `NC_TEST.1`.
    Genomic,
    /// `c.` on `NM_TEST.1` where the transcript *is* the CDS
    /// (`cds_start == 1`), so the CDS translation is the identity and both
    /// transcript boundaries are exercised.
    CdsPlus,
    /// Minus-strand twin of [`Sys::CdsPlus`].
    CdsMinus,
    /// `n.` on the non-coding `NR_TEST.1` — transcript-relative.
    NoncodingPlus,
    /// Minus-strand twin of [`Sys::NoncodingPlus`].
    NoncodingMinus,
    /// `r.` on a **coding** `NM_TEST.1` that carries a real
    /// [`PAD_OFFSET`]-base 5'UTR *and* 3'UTR, so `cds_start == PAD_OFFSET + 1`:
    /// CDS-relative, the same axis as `c.`. The 5'UTR is what gives this teeth
    /// — with `cds_start == 1` the CDS-relative and transcript-relative
    /// readings coincide and the axis would go untested — and the 3'UTR lets
    /// a 3'-shift cross `cds_end` into `*N`, which does happen here and is
    /// idempotent on this axis.
    RnaCodingPlus,
    /// Minus-strand twin of [`Sys::RnaCodingPlus`].
    RnaCodingMinus,
    /// `r.` on the **non-coding** `NR_TEST.1`: no CDS, so transcript-relative,
    /// coinciding with [`Sys::NoncodingPlus`].
    RnaNoncodingPlus,
    /// Minus-strand twin of [`Sys::RnaNoncodingPlus`].
    RnaNoncodingMinus,
    /// `c.` on the same 5'UTR- *and* 3'UTR-bearing `NM_TEST.1` fixture
    /// [`Sys::RnaCodingPlus`] uses — the `c.` spelling of that identical axis.
    ///
    /// Distinct from [`Sys::CdsPlus`], where `cds_start == 1` leaves no UTR on
    /// either side: only here can a 3'-shift run off `cds_end` into `*N`. That
    /// shape was non-idempotent until #1185/#1192 (fixed by #1189) and #1209
    /// (fixed by #1211), which is why these two systems were held out until now.
    CdsUtr3Plus,
    /// Minus-strand twin of [`Sys::CdsUtr3Plus`].
    CdsUtr3Minus,
}

/// Every coordinate system the generator draws from. Also drives the
/// non-vacuity histogram in
/// `case_strategy_covers_every_system_and_edit_kind`, so a system added here
/// without the generator actually emitting it fails that test.
///
/// `CdsUtr3*` — `c.` on a transcript with a 3'UTR — was held out of this list
/// while a 3'-shift crossing `cds_end` into `*N` was non-idempotent there. Two
/// distinct defects caused that, both now fixed: the `dup`/`delins`/`con`
/// family (#1185, fixed by #1189) and an insertion resting at `cds_end` that
/// the #387 clamp rewrote instead of letting it shift on (#1209, fixed by
/// #1211). Neither was a harness artifact — `r.` performed the identical shift
/// over the byte-identical reference in one pass throughout, which is what
/// localised both to the `c.` path. Minimal repros stay pinned in
/// `tests/it/cds_utr3_crossing_shift_idempotency.rs` and
/// `tests/it/issue_1209_cds_end_insertion_shift.rs`.
const ALL_SYSTEMS: [Sys; 11] = [
    Sys::Genomic,
    Sys::CdsPlus,
    Sys::CdsMinus,
    Sys::NoncodingPlus,
    Sys::NoncodingMinus,
    Sys::RnaCodingPlus,
    Sys::RnaCodingMinus,
    Sys::RnaNoncodingPlus,
    Sys::RnaNoncodingMinus,
    Sys::CdsUtr3Plus,
    Sys::CdsUtr3Minus,
];

impl Sys {
    /// Accession and axis prefix this system renders under.
    fn accession_and_prefix(self) -> (&'static str, &'static str) {
        match self {
            Sys::Genomic => ("NC_TEST.1", "g."),
            Sys::CdsPlus | Sys::CdsMinus | Sys::CdsUtr3Plus | Sys::CdsUtr3Minus => {
                ("NM_TEST.1", "c.")
            }
            Sys::RnaCodingPlus | Sys::RnaCodingMinus => ("NM_TEST.1", "r."),
            Sys::NoncodingPlus | Sys::NoncodingMinus => ("NR_TEST.1", "n."),
            Sys::RnaNoncodingPlus | Sys::RnaNoncodingMinus => ("NR_TEST.1", "r."),
        }
    }

    /// True for the `r.` axes, whose literal base sequences must be written in
    /// the RNA alphabet (`acgu`) rather than DNA.
    fn is_rna(self) -> bool {
        matches!(
            self,
            Sys::RnaCodingPlus
                | Sys::RnaCodingMinus
                | Sys::RnaNoncodingPlus
                | Sys::RnaNoncodingMinus
        )
    }
}

#[derive(Debug, Clone)]
struct Case {
    core: String,
    sys: Sys,
    hgvs: String,
}

fn dna(n: std::ops::RangeInclusive<usize>) -> impl Strategy<Value = String> {
    prop::collection::vec(prop::sample::select(vec!['A', 'C', 'G', 'T']), n)
        .prop_map(|v| v.into_iter().collect())
}

/// Repeat-biased core: 2-5 runs, each a 1-3bp unit repeated 1-5 times, so
/// homopolymers and short tandems dominate. Length (of the body) clamped to
/// 6..=40.
///
/// The body is wrapped in homopolymer breaker flanks. `SyntheticBuilder` pads
/// the genomic contig with a perfect `ACGT…` tandem, and its base immediately
/// 5' of the core is `T`, immediately 3' is `A`. Without a breaker, a tract at
/// the core's edge *continues* that pad tandem, creating a 256bp+ repeat an edit
/// can shift into — but only ~100bp per pass (the normalizer window), which is a
/// window-limited-tandem artifact of the harness, not a shuffle-logic bug.
///
/// The flanks are `CCCC…` / `…GGGG` (was a single `G`…`C`). The single-base
/// version only broke a *homopolymer* run against the pad (`T`+`G` is no run,
/// `C`+`A` is no run) and did NOT stop a **rotational** multi-base unit from
/// cycling straight through it: a body ending `…CAG` plus the old `C` flank
/// gives `…CAGC`, and the pad's leading `A` continues the `CAG` cycle. A
/// debugging session lost real time to exactly that — a "3' non-idempotency"
/// that turned out to be the normalizer correctly following a tract which had
/// genuinely leaked into the padding.
///
/// A 4-base homopolymer flank closes it: for a phase walk to cross `GGGG` with a
/// unit of length <= 3, every rotation it visits there must be `G`, i.e. the
/// minimal unit is the homopolymer `G` — which then stops at the pad's leading
/// `A`. Symmetrically `CCCC` stops against the pad's trailing `T`. So no tract
/// of any unit length this strategy can build reaches the padding, and a failure
/// reported against a generated core is always about the normalizer.
///
/// (`repeat_case_strategy` below needs the same property and gets it a different
/// way — it knows its tract's unit, so it derives the breaker from that unit
/// directly rather than relying on a fixed flank.)
///
/// Provenance, since the narrower single-base `G`/`C` flanks this replaced are
/// referenced elsewhere: those let a 2-base rotation continue a few periods into
/// the pad, and that leak is what originally surfaced the `cds_end`-crossing
/// non-idempotency. That case is now pinned explicitly by
/// `tests/it/cds_utr3_crossing_shift_idempotency.rs`, so widening the flanks
/// does not lose it.
fn core_strategy() -> impl Strategy<Value = String> {
    prop::collection::vec((dna(1..=3), 1usize..=5), 2..=5)
        .prop_map(|runs| {
            runs.into_iter()
                .flat_map(|(u, n)| u.repeat(n).chars().collect::<Vec<_>>())
                .collect::<String>()
        })
        .prop_filter("body length 6..=40", |s| s.len() >= 6 && s.len() <= 40)
        .prop_map(|body| format!("CCCC{body}GGGG"))
}

/// A core carrying ONE cleanly-bounded tandem tract, plus a `unit[N]` input
/// spelled against it — the input class `case_strategy` cannot reach (it only
/// emits `del`/`dup`/`ins`/`delins`), which is how the 5' repeat del/dup
/// non-idempotency hid.
///
/// The tract is `unit` repeated `ref_copies` times, flanked on both sides by a
/// run of a single **breaker** base chosen so it can continue the unit's
/// rotation on neither side: `breaker != unit.first()` blocks the 3' phase walk
/// (`three_prime_align_tract` slides while `ref[end] == rotated[0]`) and
/// `breaker != unit.last()` blocks the 5' one (`five_prime_align_tract` slides
/// while `ref[start-1] == rotated.last()`). A 1-3bp unit constrains at most 2 of
/// the 4 bases, so a breaker always exists.
///
/// This matters more than the usual "break the pad run" guard: `SyntheticBuilder`
/// pads with `ACGT…`, and a *rotational* (multi-base) unit can keep cycling into
/// that pad even where a homopolymer run would stop — e.g. a core ending `…CAGC`
/// continues the `CAG` cycle into the pad's leading `A`. Flanking by unit-derived
/// exclusion is what keeps the tract's extent equal to what the test intends.
fn repeat_case_strategy() -> impl Strategy<Value = Case> {
    const FLANK: usize = 3;
    (dna(1..=3), 2usize..=6)
        .prop_flat_map(|(unit, ref_copies)| {
            let ub = unit.as_bytes();
            let (first, last) = (ub[0] as char, ub[ub.len() - 1] as char);
            let breakers: Vec<char> = ['A', 'C', 'G', 'T']
                .into_iter()
                .filter(|&b| b != first && b != last)
                .collect();
            (
                Just(unit),
                Just(ref_copies),
                prop::sample::select(breakers),
                prop::sample::select(vec![Sys::Genomic, Sys::CdsPlus, Sys::CdsMinus]),
            )
        })
        .prop_flat_map(|(unit, ref_copies, breaker, sys)| {
            // count 0..=ref_copies+3 spans every normalize_repeat arm:
            // contraction to zero / with survivors, identity, +1 (dup), expansion.
            (
                Just(unit),
                Just(ref_copies),
                Just(breaker),
                Just(sys),
                0u64..=(ref_copies as u64 + 3),
            )
        })
        .prop_map(|(unit, ref_copies, breaker, sys, count)| {
            let flank: String = std::iter::repeat_n(breaker, FLANK).collect();
            let core = format!("{flank}{}{flank}", unit.repeat(ref_copies));
            // 1-based position of the tract's first base within the core.
            let tract_start_in_core = FLANK + 1;
            let (prefix, acc) = match sys {
                Sys::Genomic => ("g.", "NC_TEST.1"),
                _ => ("c.", "NM_TEST.1"),
            };
            let pos = match sys {
                Sys::Genomic => PAD_OFFSET as usize + tract_start_in_core,
                _ => tract_start_in_core, // cds_start == 1 => c.p maps to core p
            };
            Case {
                core,
                sys,
                hgvs: format!("{acc}:{prefix}{pos}{unit}[{count}]"),
            }
        })
}

/// Re-express a DNA literal in the RNA alphabet used by the `r.` axis
/// (lowercase, `t` -> `u`). Feeding real `u` bases through the `r.` cases also
/// exercises the U->T rewrite normalization applies before comparing an `r.`
/// edit against the DNA-stored transcript (#736).
fn to_rna(dna: &str) -> String {
    dna.chars()
        .map(|c| match c {
            'T' | 't' => 'u',
            other => other.to_ascii_lowercase(),
        })
        .collect()
}

fn case_strategy() -> impl Strategy<Value = Case> {
    core_strategy().prop_flat_map(move |core| {
        let len = core.len();
        let sys = prop::sample::select(ALL_SYSTEMS.to_vec());
        // kind: 0=del 1=dup 2=ins 3=delins 4=inv 5=con
        (
            sys,
            0u8..=5,
            1usize..=len,
            0usize..=3,
            dna(1..=3),
            // `con` donor interval, drawn from the same core so the donated
            // bases are as short as an `ins` payload (1..=3) and the resulting
            // delins cannot 3'-shift any further than the `ins`/`delins` cases
            // already can.
            1usize..=len,
            0usize..=2,
        )
            .prop_map(
                move |(sys, kind, start, span, ins, donor_start, donor_span)| {
                    let end = (start + span).min(len);
                    let (s, e) = (start.min(end), start.max(end));
                    // Every axis renders core position `p` as `p`; only the
                    // genomic contig carries the padding offset. For the
                    // 5'UTR-bearing coding fixture this holds because
                    // `cds_start == PAD_OFFSET + 1`, so `c.p`/`r.p` translate
                    // to transcript `PAD_OFFSET + p` = core position `p`.
                    let render_pos = |p: usize| -> u64 {
                        match sys {
                            Sys::Genomic => (PAD_OFFSET as usize + p) as u64,
                            _ => p as u64,
                        }
                    };
                    let (acc, prefix) = sys.accession_and_prefix();
                    let bases = if sys.is_rna() {
                        to_rna(&ins)
                    } else {
                        ins.clone()
                    };
                    let range = |s: usize, e: usize| {
                        if s == e {
                            format!("{}", render_pos(s))
                        } else {
                            format!("{}_{}", render_pos(s), render_pos(e))
                        }
                    };
                    let body = match kind {
                        0 => format!("{}del", range(s, e)),
                        1 => format!("{}dup", range(s, e)),
                        2 => {
                            // insertion sits between p and p+1; keep p in 1..len
                            let p = s.min(len - 1);
                            format!("{}_{}ins{}", render_pos(p), render_pos(p + 1), bases)
                        }
                        3 => format!("{}delins{}", range(s, e), bases),
                        4 => {
                            // "an inversion covers more than one nucleotide by
                            // definition" (HGVS DNA/inversion.md:16) and the
                            // parser enforces it, so widen a degenerate
                            // interval to two bases instead of emitting an
                            // input the generator knows is invalid.
                            let (s, e) = match (s, e) {
                                (s, e) if e > s => (s, e),
                                (s, _) if s < len => (s, s + 1),
                                (s, _) => (s - 1, s),
                            };
                            format!("{}inv", range(s, e))
                        }
                        _ => {
                            // Always render the donor as an explicit `S_E`
                            // pair, even when `S == E`: a bare `conN` parses
                            // as an opaque cross-reference (`delins[N]`)
                            // rather than a same-reference position range.
                            let ds = donor_start;
                            let de = (ds + donor_span).min(len);
                            format!("{}con{}_{}", range(s, e), render_pos(ds), render_pos(de))
                        }
                    };
                    Case {
                        core: core.clone(),
                        sys,
                        hgvs: format!("{acc}:{prefix}{body}"),
                    }
                },
            )
    })
}

fn provider_for(core: &str, sys: Sys) -> MockProvider {
    let len = core.len() as u64;
    match sys {
        Sys::Genomic => SyntheticBuilder::genomic(core).build(),
        Sys::CdsPlus => SyntheticBuilder::cds(core, 1, len, Strand::Plus).build(),
        Sys::CdsMinus => SyntheticBuilder::cds(core, 1, len, Strand::Minus).build(),
        // A coding transcript of `PAD + core + PAD` whose CDS covers exactly
        // the core, giving it both a 5'UTR and a 3'UTR.
        // `cds_start == PAD_OFFSET + 1` makes the CDS translation a real offset
        // (unlike `CdsPlus`/`CdsMinus`, where it is the identity), which is what
        // makes "`r.` is CDS-relative" an observable claim; the padding doubles
        // as the repeat-tract breaker on both sides.
        // `CdsUtr3*` is the `c.` spelling of that same axis over the same
        // fixture, so it shares the provider outright.
        Sys::RnaCodingPlus | Sys::CdsUtr3Plus => SyntheticBuilder::cds(
            &padded(core),
            PAD_OFFSET + 1,
            PAD_OFFSET + len,
            Strand::Plus,
        )
        .build(),
        Sys::RnaCodingMinus | Sys::CdsUtr3Minus => SyntheticBuilder::cds(
            &padded(core),
            PAD_OFFSET + 1,
            PAD_OFFSET + len,
            Strand::Minus,
        )
        .build(),
        // `n.` and `r.` on a non-coding transcript are likewise the same axis
        // (no CDS to be relative to), so they share the `NR_TEST.1` fixture.
        Sys::NoncodingPlus | Sys::RnaNoncodingPlus => {
            SyntheticBuilder::noncoding(core, Strand::Plus).build()
        }
        Sys::NoncodingMinus | Sys::RnaNoncodingMinus => {
            SyntheticBuilder::noncoding(core, Strand::Minus).build()
        }
    }
}

/// Edit-kind label recovered from a *rendered* case, so the histogram below
/// observes what the generator actually emitted rather than a self-reported
/// tag. Order matters: `delins` contains both `del` and `ins`.
fn edit_kind_label(hgvs: &str) -> &'static str {
    for kind in ["delins", "del", "dup", "inv", "con", "ins"] {
        if hgvs.contains(kind) {
            return kind;
        }
    }
    "?"
}

/// Non-vacuity guard for [`case_strategy`].
///
/// A passing idempotency property proves nothing about a coordinate system or
/// edit kind the strategy never emits, and a `prop_filter` or a mis-ordered
/// `match` arm can silently drop a whole bucket without any test going red.
/// This samples the strategy directly and asserts every
/// (coordinate system x edit kind) cell is actually populated.
#[test]
fn case_strategy_covers_every_system_and_edit_kind() {
    use proptest::strategy::{Strategy as _, ValueTree as _};
    use proptest::test_runner::TestRunner;
    use std::collections::BTreeMap;

    const SAMPLES: usize = 4000;
    const KINDS: [&str; 6] = ["del", "dup", "ins", "delins", "inv", "con"];
    // Derived from the enum rather than spelled out, so adding a `Sys` variant
    // to `ALL_SYSTEMS` without teaching the generator to emit it fails here.
    let systems: Vec<String> = ALL_SYSTEMS.iter().map(|s| format!("{s:?}")).collect();

    let strategy = case_strategy();
    let mut runner = TestRunner::deterministic();
    let mut histogram: BTreeMap<(String, &'static str), usize> = BTreeMap::new();
    for _ in 0..SAMPLES {
        let case = strategy
            .new_tree(&mut runner)
            .expect("case_strategy must always produce a value")
            .current();
        *histogram
            .entry((format!("{:?}", case.sys), edit_kind_label(&case.hgvs)))
            .or_default() += 1;
    }

    let mut report = String::new();
    let mut missing = Vec::new();
    for sys in &systems {
        for kind in KINDS {
            let n = histogram.get(&(sys.clone(), kind)).copied().unwrap_or(0);
            report.push_str(&format!("  {sys:16} {kind:8} {n:5}\n"));
            if n == 0 {
                missing.push(format!("{sys}/{kind}"));
            }
        }
    }
    println!("case_strategy coverage over {SAMPLES} samples:\n{report}");
    assert!(
        missing.is_empty(),
        "case_strategy never generated these (system, edit kind) cells: {missing:?}\n\
         full histogram over {SAMPLES} samples:\n{report}"
    );
}

/// `norm(norm(x)) == norm(x)` for one case in one direction, and the normalized
/// output must always re-parse — a normalizer that emits invalid HGVS (e.g. a
/// degenerate single-position `c.1insAA`) fails the re-parse step.
///
/// The render and first normalize fail LOUDLY, not by silent skip: cases are
/// valid by construction and the provider supplies the very core the coordinates
/// were derived from, so a parse/normalize error is a generator or normalizer
/// bug — and a skip branch here would let the whole property go vacuous unnoticed
/// (measured: 0 parse errors, 0 first-pass normalize errors over 4000 cases).
fn check_idempotent(case: &Case, dir: ShuffleDirection) -> Result<(), TestCaseError> {
    let norm = normalizer_for(&case.core, case.sys, dir);
    let v = parse_hgvs(&case.hgvs)
        .unwrap_or_else(|e| panic!("generator produced unparseable HGVS {}: {e}", case.hgvs));
    let n1 = norm
        .normalize(&v)
        .unwrap_or_else(|e| panic!("first normalize failed: {}: {e}", case.hgvs));
    let f1 = n1.to_string();
    // The normalized output MUST re-parse and re-normalize cleanly; an error
    // here (parse or normalize) is a bug.
    let reparsed = parse_hgvs(&f1).unwrap_or_else(|e| {
        panic!(
            "normalized output does not re-parse: {} -> {f1}: {e}\n  (dir={dir:?}, core={})",
            case.hgvs, case.core
        )
    });
    let n2 = norm
        .normalize(&reparsed)
        .unwrap_or_else(|e| panic!("re-normalize failed: {} -> {f1}: {e}", case.hgvs));
    prop_assert_eq!(
        &f1,
        &n2.to_string(),
        "NOT IDEMPOTENT\n  core={} sys={:?} dir={:?}\n  input={}\n  once ={}\n  twice={}",
        case.core,
        case.sys,
        dir,
        case.hgvs,
        f1,
        n2
    );
    Ok(())
}

proptest! {
    // 12k cases rather than 4k: the case space is 11 coordinate systems x 6 edit
    // kinds = 66 cells (it was 3 x 4 = 12), so 4k would leave only ~61 cases per
    // cell. 12k restores ~182 per cell and still runs in ~1s.
    //
    // `max_local_rejects` is raised off its 65_536 default because
    // `core_strategy`'s `body length 6..=40` filter rejects a steady ~6.7% of
    // draws, so the *absolute* default cap is reached partway through any run
    // over ~900k cases and aborts it with "Too many local rejects" — which looks
    // like a failure but is purely the harness's own quota. 1M rejects covers
    // ~14M cases while still aborting promptly if the generator ever degenerates
    // into rejecting nearly everything. Both values are overridden by
    // `PROPTEST_CASES` / `PROPTEST_MAX_LOCAL_REJECTS` (the `proptest!` macro
    // re-applies the env vars over this config), so a soak is just
    // `PROPTEST_CASES=1000000 cargo nextest run …`. 1M cases verified clean.
    #![proptest_config(ProptestConfig {
        cases: 12_000,
        max_local_rejects: 1_000_000,
        ..ProptestConfig::default()
    })]

    /// 3' (default, HGVS-canonical) direction, fully general edit content.
    #[test]
    fn normalization_is_idempotent_three_prime(case in case_strategy()) {
        check_idempotent(&case, ShuffleDirection::ThreePrime)?;
    }

    /// 5' direction, fully general edit content — the mirror of the 3' property.
    /// Enabled once the three pre-existing 5' bugs the campaign chased down
    /// (boundary dup/insertion→delins, insertion→dup tract selection, and the
    /// tandem-repeat `dup` unit rotation) were fixed.
    #[test]
    fn normalization_is_idempotent_five_prime(case in case_strategy()) {
        check_idempotent(&case, ShuffleDirection::FivePrime)?;
    }

    /// Repeat-notation (`unit[N]`) **inputs**, both directions. `case_strategy`
    /// never emits these, so `normalize_repeat`'s own result arms were only ever
    /// reached as an *output* of a del/dup/ins — never re-entered as an input.
    /// That gap hid the 5' del/dup arms naming the 3'-most copy (non-idempotent).
    #[test]
    fn repeat_input_is_idempotent(case in repeat_case_strategy(), dir in both_directions()) {
        check_idempotent(&case, dir)?;
    }
}

/// Explicit pin for the #1169 multi-base 5'-shift rotation bug (`insCA` emitted
/// unrotated at the shifted position — a *different haplotype*, not merely a
/// non-canonical spelling).
///
/// Its only coverage was the `bf332dd0…` seed in
/// `tests/proptest-regressions/normalize_idempotency_proptest.txt`, recorded as
/// `core: "GTGTGTGTGTGCTCC"`. **Seeds are replayed through the *current*
/// strategy, not stored cases** — and this PR rewrites `core_strategy`'s flanks
/// to `CCCC…GGGG`, so no core it can now generate equals that string. The seed
/// therefore silently stops reproducing this case. Pinning it as a concrete test
/// makes the coverage immune to any future strategy change.
///
/// The odd 11-base `GTGTGTGTGTG` run is the ambiguous-phase shape that made the
/// rotation matter: inserting `TG` has two equally-valid spellings, and the 5'
/// path must rotate the inserted bases as it slides rather than carry them
/// verbatim.
#[test]
fn issue_1169_multibase_five_prime_rotation_is_idempotent() {
    let case = Case {
        core: "GTGTGTGTGTGCTCC".to_string(),
        sys: Sys::Genomic,
        hgvs: "NC_TEST.1:g.267_268insTG".to_string(),
    };
    check_idempotent(&case, ShuffleDirection::FivePrime)
        .expect("#1169 regression: 5' multi-base insertion rotation must be a fixed point");
}

/// A homopolymer confluence case: a run of `base` of length `run_len`, embedded
/// in the CDS, with two distinct HGVS spellings of the SAME haplotype (inserting
/// `unit_len` more copies of `base` into the run at two different positions).
/// Both spellings describe an identical resulting sequence, so a confluent
/// normalizer must map them to the same normal form.
#[derive(Debug, Clone)]
struct ConfluenceCase {
    core: String,
    sys: Sys,
    a: String,
    b: String,
}

/// Build a pure-homopolymer or boundary-adjacent run and two equivalent
/// insertion spellings within it. Using `cds_start == 1`/`cds_end == len`
/// exercises the CDS-start and CDS-end boundary clamps (the confluence-critical
/// region), where all spellings must collapse to one `c.<boundary>delins…`.
fn confluence_case_strategy() -> impl Strategy<Value = ConfluenceCase> {
    (
        prop::sample::select(vec!['A', 'C', 'G', 'T']),
        4usize..=12, // run length
        1usize..=3,  // number of extra copies inserted
        prop::sample::select(vec![Sys::CdsPlus, Sys::CdsMinus]),
    )
        .prop_flat_map(|(base, run_len, extra, sys)| {
            // The whole CDS is the homopolymer run, so c.k maps to core k and
            // both transcript boundaries participate.
            let core: String = std::iter::repeat_n(base, run_len).collect();
            let unit: String = std::iter::repeat_n(base, extra).collect();
            // Insertion sits between positions p and p+1 (1-based), p in 1..run_len.
            // Any two distinct p give the same haplotype (run grows by `extra`).
            let core_c = core.clone();
            let unit_c = unit.clone();
            // `p2` is drawn as a non-zero offset from `p1` modulo the position
            // count, so the two spellings are always DISTINCT. Sampling both
            // independently would let `p1 == p2` through, and comparing a
            // normalized form against itself is trivially equal — a silent
            // weakening of the confluence property on those cases.
            let slots = run_len - 1; // positions are 1..run_len
            (1usize..run_len, 1usize..slots).prop_map(move |(p1, offset)| {
                let p2 = 1 + ((p1 - 1 + offset) % slots);
                let acc = "NM_TEST.1";
                let spell = |p: usize| format!("{acc}:c.{}_{}ins{}", p, p + 1, unit_c);
                ConfluenceCase {
                    core: core_c.clone(),
                    sys,
                    a: spell(p1),
                    b: spell(p2),
                }
            })
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(4000))]

    /// Two spellings of one haplotype (insert the same unit at two positions in
    /// a homopolymer run) must normalize to the SAME form, in both directions.
    /// This catches non-confluence that idempotency alone cannot: before the
    /// boundary-clamp rework the 5' direction mapped `c.1_2insAA`,
    /// `c.5_6insAA`, … to distinct (each individually idempotent)
    /// `c.1delinsAAA` / `c.1_3delinsAAAAA` / … forms.
    #[test]
    fn normalization_is_confluent(case in confluence_case_strategy(), dir in both_directions()) {
        let norm = normalizer_for(&case.core, case.sys, dir);
        // No skip branches: the generator only emits well-formed spellings, so a
        // parse or normalize failure is a BUG, not a case to pass over. Returning
        // `Ok(())` here would let the whole property go vacuous unnoticed — the
        // exact hazard `check_idempotent` documents and avoids.
        let va = parse_hgvs(&case.a)
            .unwrap_or_else(|e| panic!("generator produced unparseable HGVS {}: {e}", case.a));
        let vb = parse_hgvs(&case.b)
            .unwrap_or_else(|e| panic!("generator produced unparseable HGVS {}: {e}", case.b));
        let na = norm
            .normalize(&va)
            .unwrap_or_else(|e| panic!("normalize failed: {}: {e}", case.a));
        let nb = norm
            .normalize(&vb)
            .unwrap_or_else(|e| panic!("normalize failed: {}: {e}", case.b));
        prop_assert_eq!(
            na.to_string(),
            nb.to_string(),
            "NOT CONFLUENT\n  core={} sys={:?} dir={:?}\n  a={} b={}",
            case.core,
            case.sys,
            dir,
            case.a,
            case.b
        );
    }
}
