//! Driver for the exhaustive, non-redundant HGVS spec test enumeration.
//!
//! **What a green run here does and does not prove.** Every assertion in this
//! file is either a *rejection/repair* assertion (input side) or a
//! *well-formedness* assertion over the string ferro emits. None of it is a
//! correctness oracle: "is this legal HGVS?" is not "is this the right answer".
//! The correctness oracles are the differential corpora
//! (`mutalyzer_normalize_tests.rs`, `biocommons_normalize_tests.rs`,
//! `hgvs_rs_projection_tests.rs`). A green enumeration means ferro's output
//! parses and obeys the MUST-level shape rules — nothing more.
//!
//! The fixture is generated and gitignored; regenerate it with
//! `cargo run --features dev --bin generate_spec_enumeration`.
//!
//! **Which assertions here are real oracles (#1272).** Because the fixture is
//! regenerated from the code under test before every CI run,
//! [`enumeration_replays_recorded_behavior`] compares ferro against itself and
//! can only fire on a stale local artifact. The assertions that can actually
//! catch a behaviour change are the ones backed by something *committed*:
//!
//! - [`DIVERGENCE_BUDGET`] and [`PASSING_CENSUS`] — per-status row counts, which
//!   between them total every row.
//! - [`passing_spec_mandated_musts_stay_passing`] — replays against `expected`,
//!   which comes from the committed overrides file, not from ferro.
//!
//! Read a drift in the replay test as "decide whether this change was
//! deliberate", never as "regenerate".
//!
//! Rows whose recorded behaviour diverges from the spec are **not** failures
//! here. They are classified (`repair-diverges`, `false-acceptance`,
//! `invariant-violation-must`, …) and counted against a committed budget, so
//! the suite stays green while the divergence is recorded and any regression
//! or improvement shows up as a budget mismatch.

use std::collections::BTreeMap;
use std::path::PathBuf;

use ferro_hgvs::conformance::reference_window::WindowProvider;
use ferro_hgvs::conformance::spec_projection;
use ferro_hgvs::conformance::{Expectation, NormativeLevel, Status};
use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::project::VariantProjector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use serde::Deserialize;

/// Committed hermetic reference slice backing the `project-*` dimensions.
/// Anchored on `CARGO_MANIFEST_DIR` (via the shared helper) so the path
/// resolves regardless of the test's working directory. Unlike the generated
/// fixtures this slice is committed, so it is returned directly — no
/// `ensure_generated_fixture` regeneration step.
fn projection_windows_path() -> PathBuf {
    crate::common::fixture_gen::fixture_path("tests/fixtures/grammar/spec_enumeration_windows.json")
}

#[derive(Debug, Deserialize)]
struct Enumeration {
    rows: Vec<Row>,
}

#[derive(Debug, Deserialize)]
struct Row {
    id: String,
    dimension: String,
    operation: String,
    error_mode: String,
    target: String,
    expectation: Expectation,
    normative_level: NormativeLevel,
    #[serde(default)]
    expected: Option<String>,
    observed: String,
    status: Status,
    spec_citation: String,
}

fn load() -> Enumeration {
    crate::common::spec_enumeration::ensure_spec_enumeration();
    let path = crate::common::spec_enumeration::spec_enumeration_path();
    let text = std::fs::read_to_string(&path).expect("read enumeration fixture");
    serde_json::from_str(&text).expect("parse enumeration fixture")
}

/// Everything a replay needs: the hermetic normalizer, plus the projector built
/// from the committed reference slice (absent only on a checkout that has not
/// been given the slice, in which case projection rows are skipped rather than
/// failed).
struct Replayer {
    normalizer: Normalizer<MockProvider>,
    projector: Option<VariantProjector<WindowProvider>>,
}

impl Replayer {
    fn new() -> Self {
        // The slice is committed, so its absence means only a checkout that has
        // not materialised it — skip the `project-*` rows then. But if the file
        // *is* present, any load failure (malformed JSON, a missing sidecar
        // FASTA) is a corrupt fixture: fail loudly rather than silently drop all
        // projection coverage, which a blanket `.ok()` would have done.
        let windows_path = projection_windows_path();
        let projector = if windows_path.exists() {
            let windows = spec_projection::load_slice(&windows_path).unwrap_or_else(|e| {
                panic!(
                    "committed projection slice {} is present but failed to load \
                     (corrupt fixture?): {e}",
                    windows_path.display()
                )
            });
            let cdot = CdotMapper::from_transcripts(windows.transcripts.iter());
            Some(VariantProjector::new(
                Projector::new(cdot),
                windows.to_provider(),
            ))
        } else {
            None
        };
        Replayer {
            normalizer: Normalizer::new(MockProvider::new()),
            projector,
        }
    }
}

/// Re-run a row's operation against live ferro.
///
/// Returns `None` for rows whose operation cannot be replayed from the fixture
/// alone (`invariant-check` needs the generator's catalog, which lives with the
/// generator so the rules and their spec citations stay in one place).
fn replay(row: &Row, ctx: &Replayer) -> Option<String> {
    let normalizer = &ctx.normalizer;
    if let Some(axis) = row.operation.strip_prefix("project-") {
        let code = axis.chars().next()?;
        let projector = ctx.projector.as_ref()?;
        let variant = parse_hgvs(&row.target).ok()?;
        let pass = spec_projection::project_all_axes(projector, &variant);
        return pass.axes.get(&code).map(|r| r.as_observed());
    }
    match row.operation.as_str() {
        "reject" | "normalize" => Some(match parse_hgvs(&row.target) {
            Err(e) => format!("parse error: {e}"),
            Ok(v) => match normalizer.normalize_with_diagnostics(&v) {
                Err(e) => format!("normalize error: {e}"),
                Ok(n) => format!("{}", n.result),
            },
        }),
        "parse" if row.dimension == "grammar-form" => Some(match parse_hgvs(&row.target) {
            Ok(v) => match v.coordinate_axis() {
                Some(ax) => format!("axis={}", ax.code()),
                None => "axis=none".to_string(),
            },
            Err(e) => format!("parse error: {e}"),
        }),
        "parse" => {
            let cfg = match row.error_mode.as_str() {
                "strict" => ErrorConfig::strict(),
                "lenient" => ErrorConfig::lenient(),
                "silent" => ErrorConfig::silent(),
                _ => return None,
            };
            Some(match parse_hgvs_with_config(&row.target, cfg) {
                Err(e) => format!("parse error: {e}"),
                Ok(r) => {
                    let mut codes: Vec<String> = r
                        .warnings
                        .iter()
                        .map(|w| format!("{:?}", w.error_type))
                        .collect();
                    codes.sort();
                    codes.dedup();
                    if codes.is_empty() {
                        format!("{}", r.result)
                    } else {
                        format!("{} warnings={}", r.result, codes.join(","))
                    }
                }
            })
        }
        _ => None,
    }
}

/// Every replayable row must reproduce its recorded behaviour exactly. This is
/// the drift lock: it catches a change in ferro that the generator recorded but
/// nothing else asserts.
#[test]
fn enumeration_replays_recorded_behavior() {
    let fx = load();
    let ctx = Replayer::new();
    let mut diffs = Vec::new();
    let mut replayed = 0usize;

    for row in &fx.rows {
        let Some(observed) = replay(row, &ctx) else {
            continue;
        };
        replayed += 1;
        if observed != row.observed {
            diffs.push(format!(
                "  id       : {}\n    target   : {}\n    recorded : {}\n    observed : {}",
                row.id, row.target, row.observed, observed
            ));
        }
    }

    eprintln!(
        "spec enumeration: replayed {replayed} of {} rows",
        fx.rows.len()
    );
    assert!(replayed > 300, "replayed too few rows ({replayed})");
    assert!(
        diffs.is_empty(),
        "{} enumeration row(s) drifted.\n\n\
         This can only mean your *local* fixture is stale: it is gitignored and \
         regenerated before every CI run, so on CI the recorded and observed sides \
         always agree and this assertion cannot fire. It is a stale-artifact detector, \
         not a regression detector — see the module docs and #1272.\n\n\
         So do NOT reach for the regenerate command first. A drift here means your \
         code changed behaviour since the artifact was written, and the question is \
         whether that was deliberate. The guards that actually judge it are the \
         committed ones — DIVERGENCE_BUDGET and PASSING_CENSUS — so run the full \
         module after regenerating and see whether they move:\n  \
         cargo run --features dev --bin generate_spec_enumeration\n\n{}",
        diffs.len(),
        diffs
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// Every row must carry usable provenance: a spec citation pinned to the
/// submodule SHA. The `normative_level` and `expectation` vocabularies are now
/// enforced structurally — an unknown value fails deserialization by naming
/// itself — so the only cross-field invariant left to assert is that a pinned
/// baseline never masquerades as a MUST-level spec expectation.
#[test]
fn every_row_carries_provenance() {
    let fx = load();
    for row in &fx.rows {
        assert!(
            row.spec_citation.contains('@'),
            "row {} has no SHA-pinned citation: {:?}",
            row.id,
            row.spec_citation
        );
        // A pinned baseline must never masquerade as a spec expectation.
        if row.expectation == Expectation::PinnedBaseline {
            assert_ne!(
                row.normative_level,
                NormativeLevel::Must,
                "row {} pins current behaviour but claims MUST-level force",
                row.id
            );
        }
    }
}

/// Divergence budget.
///
/// Each entry is the number of rows currently in a status that records ferro
/// diverging from the spec. The suite stays green, but any change — a
/// regression *or* a fix — trips this and must be re-blessed deliberately.
///
/// Regenerate the numbers with:
///   `cargo run --features dev --bin generate_spec_enumeration -- --census`
/// and read `by_status` in the generated fixture.
const DIVERGENCE_BUDGET: &[(Status, usize)] = &[
    // Spec-forbidden strings ferro parses and renders anyway.
    (Status::FalseAcceptance, 6),
    // Spec names a canonical replacement; ferro accepts the bad string and
    // renders something else. These are the genuine violations.
    (Status::RepairDiverges, 4),
    // Spec names a canonical replacement; ferro rejects the bad string instead
    // of repairing it. Spec-conformant — rejection is always permitted.
    (Status::RejectedNotRepaired, 13),
    // Repairs that need real reference bases; not assertable hermetically.
    (Status::RequiresReference, 10),
    // MUST-level output invariants violated by a string ferro emits.
    (Status::InvariantViolationMust, 0),
    // SHOULD-level (advisory) output invariants. Never a hard failure.
    (Status::InvariantViolationShould, 2),
    // syntax.yaml examples that do not parse into their declared axis.
    (Status::FormAxisDiverges, 3),
    // Projections whose rendered form differs from the one the spec states.
    (Status::ProjectionDiverges, 0),
    // Single-member inputs that project to a multi-member allele — the form
    // `DNA/delins.md:42` calls "not correct". See
    // `Status::ProjectionSplitsSingleMember`; added by #1272, which is also what
    // makes the LRG_199 delins regression visible.
    //
    // 10 -> 9 in #1271: that LRG_199 row is the one this status was created to
    // surface, and extending `separations_are_meaningful` to net deletions is
    // what stops it splitting. The row it vacates reappears in
    // `PASSING_CENSUS`'s `ProjectionPinned` (1167 -> 1168), so the total at that
    // time (2184 rows) was unchanged — nothing left the enumeration, one row
    // moved from a divergence to a pass. The total is 2172 as of #1498, which
    // did drop rows; see `ModeDivergencePinned` below.
    //
    // **#1235's axis widening leaves this constant where it was.** Do not read
    // the note below as a re-bless: the branch measured the four rows *in*
    // (9 -> 13) at the widening commit and measured them back *out* (13 -> 9)
    // once the rest of the branch landed, so this value and its mirror in
    // `PASSING_CENSUS` (`ProjectionPinned`, 1168) are byte-identical to `main`.
    // The analysis is kept because the four rows are the class this status
    // exists for and the shape will recur.
    //
    // All four are one input, `NM_004006.2:c.76_83inv` (`general.md:110`),
    // reaching `project-c`, `project-n`, `project-r` and `project-p`:
    //
    //   c.76_83inv -> c.[76_77delinsTG;82_83delinsTT]
    //   p.(Asn26_Gln28delinsCysAlaLeu) -> p.[(Asn26Cys);(Gln28Leu)]
    //
    // Read the status name and not only its headline citation: these enter by
    // *shape* (a single member projecting to an allele), not by the divergence
    // `delins.md:42` names — that passage is about two variants **one** base
    // apart, and these are four apart, where `general.md:34` plainly asks for
    // individual descriptions.
    //
    // What the four rows really record is that the transcript axes now answer as
    // the genomic axis already did. `c.76_83` is `AATGCACA`, whose reverse
    // complement `TGTGCATT` coincides at four of its eight columns, so the
    // partition finds two runs 4 bases apart and `coalesce_whole_block_inversion`
    // declines to rejoin them: 50% unchanged is above the ~25% a random reverse
    // complement leaves, so `changed_columns_dominate_the_span` reads it as a
    // near-palindrome carrying two independent changes rather than one inversion.
    // None of that reads the axis — the block is equal-length, so neither
    // `separations_are_meaningful` nor `coalesce_coding_frame_separation` is even
    // consulted — so `g.` has always produced this split for this block. Widening
    // the gate did not create the behaviour; it let the transcript axes see it.
    //
    // Had the four rows stayed, they would have been kept rather than chased:
    // making them pass means loosening `coalesce_whole_block_inversion`, which
    // is axis-neutral and shipped, so it would move genomic representations for
    // a class this PR measures nowhere. That is its own change with its own
    // blast radius, not a tail of this one.
    //
    // **Read the pair, never either constant alone.** The two deltas are equal
    // and opposite by construction (`ProjectionSplitsSingleMember` down four,
    // `ProjectionPinned` up four), so a row that vanished from one without
    // appearing in the other would mean something left the enumeration rather
    // than moved between statuses — and the 2172-row total would be the only
    // thing that showed it.
    //
    // **And measure against a regenerated artifact.** All of the above was
    // invisible until this branch was rebased onto a `main` carrying #1535 and
    // #1537: the enumeration artifact is gitignored, so a local file generated
    // against the pre-rebase tree had these two committed guards comparing
    // against a stale recording, and the movement in both directions went
    // unseen. That is exactly why these are the guards that judge behaviour and
    // `enumeration_replays_recorded_behavior` is not.
    //
    // **#1672 holds this at 9; it moves no row.** Two of the nine are the
    // `LRG_199t1` rows `c.145_147delinsTGG` (`delins.md:37`) and
    // `c.235_237delinsTAT` (`delins.md:19`), whose derived `g.` axis renders the
    // split `delins.md:42` calls "not correct". #1664 reads that as a defect and
    // asks for the genomic axis to be merged to match the coding one, which
    // would take this constant to 7.
    //
    // It is declined. `:42` is a conditional whose second conjunct is "together
    // affecting one amino acid", and that cannot be stated on a genomic
    // reference — `general.md:23-31` makes the prefix a claim about the **type
    // of reference sequence**, so `LRG_199:g.…` is genomic however gene-scoped
    // its accession. With `:42` not firing, `DNA/delins.md:17` governs unopposed
    // and asks for exactly these individual descriptions. The two rows are
    // therefore conformant output, and this status is where conformant-but-split
    // rows belong.
    //
    // See `rulings[projection-codon-exception-is-decided-by-the-rendered-axis]`,
    // which is what #1672 contributes here: the ruling that KEEPS this at 9.
    // Declining the merge is also what lets `axis_genomic_idempotent` assert
    // every projected genomic axis is a fixed point with no exemption list at
    // all — a merged genomic string is not re-derivable from its own reference.
    //
    // **Driving this to zero is not the goal, and neither is driving it to 7.** A
    // change that takes it below 9 has merged something `DNA/delins.md:17` wanted
    // split. Read it with its mirror in `PASSING_CENSUS` (`ProjectionPinned`,
    // 1168): the two constants move in equal and opposite steps, so a row that
    // vanished from one without appearing in the other left the enumeration
    // rather than moving between statuses. The other seven are the shape doing
    // its job — `NM_000797.3:c.812_829delins908_925` over four axes and
    // `c.235_238delinsTAGT` over three, whose members end up two or more bases
    // apart, where `general.md:34` asks for the split.
    //
    // # 9 -> 12 (#1649), and why the paragraph directly above does not forbid it
    //
    // The warning above is one-directional: it forbids going **below** 9, on the
    // ground that a row leaving this status has been merged where
    // `DNA/delins.md:17` wanted it split. This move is the other way. Three rows
    // enter, and they are three projections of **one input**:
    //
    // ```text
    // project-c/NM_004006.3:r.123_127delinsag
    // project-n/NM_004006.3:r.123_127delinsag
    // project-r/NM_004006.3:r.123_127delinsag
    // ```
    //
    // `PASSING_CENSUS`'s `ProjectionPinned` falls by the same three (1168 ->
    // 1165), so the pair still moves in equal and opposite steps and no row left
    // the enumeration — which is the check that paragraph asks for, applied.
    //
    // **OPERATOR RULING, 2026-08-11 — the three rows are accepted, and they are
    // not one thing.** They were measured through `spdi::compare_denoted_sequences`
    // against the real reference and denote the input's own bases, so this is a
    // representation move and not a correctness regression.
    //
    // - The `n.` row is a **fix**. #1682's payload-coincidence carve-out is
    //   `c.`-scoped by `rulings[delins-payload-coincidence-carve-out-is-coding-dna-scoped]`,
    //   so `general.md:34` governs the `n.` axis unopposed and the split is what
    //   it asks for. This row was previously merged and should not have been.
    // - The `r.` row is **unruled**. The same record puts `r.` out of the
    //   carve-out in both directions — a DNA document has no jurisdiction over
    //   the RNA axis — but `RNA/delins.md` states no counterpart clause either
    //   way, so nothing selects between the two legal forms here.
    // - The `c.` row is the **residue**, and the only one. Its real question is
    //   whether `coalesce_payload_alignment_split` should be promoted onto the
    //   shipping path: the pass is gated on `PartitionRule::CanonicalCoalesced`
    //   while the shipped rule is `Live`. That is a separate change which #1649
    //   neither makes nor blocks, and it is deliberately not made here.
    // # 12 -> 21 (#1835, the partition default flip)
    //
    // NINE rows enter, and `PASSING_CENSUS`'s `ProjectionPinned` falls by exactly
    // nine (1167 -> 1158), so the pair still moves in equal and opposite steps
    // and no row left the enumeration. **That was measured rather than inferred**:
    // the enumeration was regenerated under `FERRO_PARTITION=live` and diffed
    // against the default arm's by row id — same 2 172 rows on both sides, no id
    // present in one and absent from the other, and exactly nine status changes,
    // every one of them `projection-pinned -> projection-splits-single-member`.
    // That is the check the paragraph above asks for, applied.
    //
    // The nine are **four inputs across their axes**, not nine independent rows:
    //
    // ```text
    // LRG_199t1:c.992_1004delinsAC                    project-{c,g,n,r}   (4)
    // NM_004006.2:r.2623_2803delins2804_2949          project-{c,n,r}     (3)
    // LRG_199t1:c.850_901delinsTTCCTCGATGCCTG         project-g           (1)
    // LRG_199t1:c.9002_9009delinsTTT                  project-g           (1)
    // ```
    //
    // **The direction is the licensed one.** The warning above forbids going
    // BELOW 9, on the ground that a row leaving this status has been merged where
    // `DNA/delins.md:17` wanted it split. This is the other way: nine rows that
    // were merged are now described individually.
    //
    // The two `project-g` rows are the axis scope arriving.
    // `delins-payload-coincidence-carve-out-is-coding-dna-scoped` (decided) puts
    // `DNA/delins.md:47` off the genomic axis, so `general.md:34` governs there
    // unopposed. Note `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` is
    // `delins.md:44-47`'s OWN worked example: its `g.` projection splits here
    // while its `c.` form merges, which is the axis scope doing exactly what it
    // says and is worth reading as the clearest single demonstration of it.
    //
    // The `LRG_199t1:c.992_1004delinsAC` family reaches
    // `c.[992_993del;995_997del;999_1004del]` — three PURE deletions, so
    // `delins-recommendation-reach-when-the-input-arrives-split` (decided) keeps
    // `:47` away from it on every axis: nothing was inserted, so `:46`'s
    // re-alignment mechanism cannot have occurred. That record cites this very
    // derivation by name.
    (Status::ProjectionSplitsSingleMember, 21),
];

/// Census of the **conformant** statuses — the counting half of
/// [`KNOWN_PASSING_STATUSES`].
///
/// **Why this exists (#1272).** `divergence_budget_is_unchanged` iterates only
/// [`DIVERGENCE_BUDGET`], so every status in `KNOWN_PASSING_STATUSES` was
/// counted by nothing at all. That is 2125 of the enumeration's 2172 rows, free
/// to move among themselves with no signal — a projection could go from
/// rendering a string to erroring and no test would notice, because
/// `every_observed_status_is_accounted_for` only rejects status names it has
/// never seen, not changes in how many rows carry each one.
///
/// Together with `DIVERGENCE_BUDGET` this makes the census total, so every row
/// is counted by exactly one committed table.
///
/// **Known residual.** These are counts, not row-id sets, so a *compensating*
/// swap — one row entering a status while another leaves it — is invisible.
/// Committing 2136 row ids would close that, at the cost of re-creating the
/// merge-conflict magnet that got the fixture decommitted in the first place.
/// Counts catch every non-compensating move, which is the realistic regression.
const PASSING_CENSUS: &[(Status, usize)] = &[
    (Status::CorrectlyRejected, 8),
    (Status::Repaired, 6),
    (Status::FormAxisOk, 120),
    (Status::FormAxisPinned, 0),
    (Status::Preserved, 6),
    // #1264 moved 19 rows between three statuses. Stated gross, because the net
    // deltas alone do not determine them — the three nets are linearly dependent
    // (any `Error -> Pinned` count satisfies all three if the other two absorb
    // it), so a reader given only the nets cannot recover what actually moved:
    //
    //   14 rows  ProjectionPinned      -> ProjectionUnavailablePinned
    //    4 rows  ProjectionErrorPinned -> ProjectionUnavailablePinned
    //    1 row   ProjectionErrorPinned -> ProjectionPinned
    //
    // Net: ProjectionPinned -13, ProjectionErrorPinned -5,
    // ProjectionUnavailablePinned +18. The deltas cancel, so the total at that
    // time (2184 rows) across this census and `DIVERGENCE_BUDGET` was
    // unchanged — nothing was dropped from the enumeration, only moved. Every row that moved to
    // `Unavailable` was a projection that had no business rendering:
    //
    // * `r.spl` / `r.spl?` / `r.0` are RNA-level *effects*, not sequence
    //   changes, so they have no g./c./n./p. counterpart. The projector had no
    //   `NaEdit::Splice` handling at all and carried the edit token onto the
    //   target axis verbatim, emitting `NC_000023.11:g.spl`, `…:c.1spl` and
    //   `…:n.245spl` — none of which ferro can re-parse. It also invented a
    //   protein consequence, `NP_003997.1:p.(Met1?)`, for a variant that states
    //   no coding change. These now decline with a reason.
    // * `project-g/LRG_199t1:r.1149_1150insn[100]` rendered
    //   `LRG_199:g.699646_700297insNNN…`, an insertion anchored on two bases 652
    //   apart because they straddle an intron — the exact string reported in
    //   #1264, and one ferro's own parser rejects per `DNA/insertion.md:15`.
    //
    // `r.0` on `NM_004006.1` accounts for all five `Error` departures, and is
    // the reason declining per axis rather than per projection matters. It
    // previously failed the *whole* projection with "Transcript has no CDS", so
    // all five of its rows were errors. Now the one axis that can carry the
    // statement renders it — `project-r/…` moves to `Pinned` as `r.(0)` — while
    // the four that cannot (`project-c`, `project-g`, `project-n`, `project-p`)
    // move to `Unavailable` with a reason. One row improved to rendered, not two.
    //
    // The four `project-r/…r.spl` rows are deliberately unmoved: the `rna` axis
    // is documented as the predicted form and already rendered `r.(spl)`, so the
    // replacement path preserves that wrapper rather than quietly changing it.
    // 1167 -> 1168 in #1271: the LRG_199 `delins.md:44` row that used to land in
    // `ProjectionSplitsSingleMember` now projects as the single member the spec
    // states. See that status's note in `DIVERGENCE_BUDGET`.
    // #1235's axis widening leaves this at 1168. The four `c.76_83inv` rows went
    // out to `ProjectionSplitsSingleMember` at the widening commit (1168 -> 1164)
    // and came back once the rest of that branch landed (1164 -> 1168), so the
    // net is zero and the total is again unchanged. The analysis of those rows,
    // and the reason the two constants must always be read as a pair, are with
    // `ProjectionSplitsSingleMember` in `DIVERGENCE_BUDGET`.
    // #1672 leaves this at 1168 and moves no row. Merging the derived genomic
    // axis of the two `LRG_199t1` rows — what #1664 asks for — would take this
    // to 1170 and `ProjectionSplitsSingleMember` to 7; the ruling declines it,
    // so both constants stand. See that status's note in `DIVERGENCE_BUDGET`
    // for why, and read the two as a pair either way.
    // 1168 -> 1165 (#1649). The three rows are the `project-{c,n,r}` axes of
    // `NM_004006.3:r.123_127delinsag`, moving to `ProjectionSplitsSingleMember`
    // — equal and opposite, as this pair must always be read. The adjudication
    // is recorded with that status in `DIVERGENCE_BUDGET`: sequence-preserving
    // in all three, a fix on `n.`, unruled on `r.`, and one residue on `c.`.
    //
    // 1165 -> 1167 (#1704), and it is an improvement rather than a
    // re-partition: `project-g/c.1704del` and `project-g/c.1704dup` used to read
    // `unavailable: no g. representation for this variant` and now render
    // `NC_000023.11:g.32573745del` / `dup`. `NM_004006.2:c.1704del` 3'-shifts
    // across the exon/intron junction, and while its `c.` answer was spelled on a
    // BARE `NM_004006.2` there was no genomic reference in the description for the
    // genomic axis to derive from. Re-parenting it onto
    // `NC_000023.11(NM_004006.2)` — what `checklist.md:20` requires anyway —
    // restores the anchor, so the projection resolves. `ProjectionUnavailablePinned`
    // takes the matching -2, below. The two changes touch disjoint rows, but the
    // figure was re-measured on the rebase rather than composed from the deltas.
    //
    // Four further rows change their rendered string without changing status, and
    // are a representation change rather than a capability one: `project-c` and
    // `project-n` for the same two inputs now read
    // `NC_000023.11(NM_004006.2):c.1704+1del` / `:n.1948+1del` (and the `dup`
    // pair). Those are the same accession repair, on the axes that were already
    // rendering.
    // #1835: 1167 -> 1158. The nine rows that left are exactly the nine that
    // entered `DIVERGENCE_BUDGET`'s `ProjectionSplitsSingleMember`, verified by a
    // row-id diff of the enumeration regenerated under both arms rather than by
    // the counts agreeing — see that constant for the list and the licensing.
    //
    // 1158 -> 1154 (#1870), with the matching -6 and +10 immediately below. All
    // ten moved rows are drawn from the four `NM_004006.1:c.…` inputs across the
    // five projection axes — a pool of twenty, of which the other ten were
    // already reporting the error and so did not move — and `NM_004006.1` is a
    // transcript version whose CDS the
    // committed reference slice does not resolve. Under the decided
    // `rulings[c-description-against-an-unresolvable-cds-is-refused]` a `c.`
    // description against such a record is refused, so the axes that used to
    // render one now report the same `has no CDS start` conversion error the
    // other axes of the same rows were ALREADY reporting before this change —
    // which is why the three deltas sum to zero and `DIVERGENCE_BUDGET` is
    // untouched. Read all three as one move. The -4 is disjoint from #1835's
    // -9 above — measured on the rebase by regenerating the enumeration, not
    // composed from the two branches' deltas.
    //
    // Worth knowing, and recorded in the ruling record too: `refseq.md:145`
    // publishes `NM_004006.1:c.5697del` as its own worked minus-strand
    // example, so a spec example is in the moving set. That is a statement
    // about the reference this corpus is built on, not about the description —
    // and `project` refused those rows already.
    (Status::ProjectionPinned, 1154),
    // #1704: 487 -> 485, the two `project-g/c.1704{del,dup}` rows that stopped
    // being unavailable. See the note above; read the pair together, since a move
    // between these two statuses is invisible in either number alone.
    // #1870: 485 -> 479. See `ProjectionPinned` above — same ten-row move.
    (Status::ProjectionUnavailablePinned, 479),
    // #1870: 210 -> 220, the +10 receiving the two deltas above.
    (Status::ProjectionErrorPinned, 220),
    // 132 -> 120 (#1498). The 12 rows are all LRG — `LRG_199:c.357+1G>A`,
    // `LRG_199:g.954966C>T`, `LRG_199:g.981731G>A` and `LRG_476:g.4950_39800=`,
    // one per mode — and **no row entered**, which is what says this is the LRG
    // population and not a wider move. Diffed row-by-row against `origin/main`
    // rather than inferred from the net, since a -12 is also what "13 left, 1
    // arrived" looks like.
    //
    // They did not become conformant: with the bare-`LRG_<N>` version check
    // fixed they now reach projection, where the committed reference slice
    // cannot serve `LRG_199` / `LRG_163t1`, so they are excluded as a fixture
    // gap. That is why the enumeration total drops 2184 -> 2172 rather than
    // staying flat. The generator names the absent transcripts on every run, so
    // the exclusion is reported rather than silent.
    (Status::ModeDivergencePinned, 120),
];

/// Does a projected string match the form the spec states? Mirrors
/// `generate_spec_enumeration`'s rule: a fully-qualified spec value is compared
/// verbatim, an accession-less one against the coordinate part only.
fn projection_matches(observed: &str, expected: &str) -> bool {
    if expected.contains(':') {
        observed == expected
    } else {
        observed.rsplit(':').next() == Some(expected)
    }
}

#[test]
fn divergence_budget_is_unchanged() {
    let fx = load();
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for row in &fx.rows {
        *counts.entry(row.status.as_str()).or_default() += 1;
    }
    let mut mismatches = Vec::new();
    for (status, budget) in DIVERGENCE_BUDGET {
        let actual = counts.get(status.as_str()).copied().unwrap_or(0);
        if actual != *budget {
            mismatches.push(format!("  {status}: budget {budget}, actual {actual}"));
        }
    }
    assert!(
        mismatches.is_empty(),
        "spec-divergence budget changed. If ferro improved, lower the budget; if it \
         regressed, this is a real defect. Update DIVERGENCE_BUDGET deliberately.\n{}",
        mismatches.join("\n")
    );
}

/// The conformant half of the census. See [`PASSING_CENSUS`] for why counting
/// these matters — before #1272 they were counted by nothing.
#[test]
fn passing_census_is_unchanged() {
    let fx = load();
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for row in &fx.rows {
        *counts.entry(row.status.as_str()).or_default() += 1;
    }
    let mut mismatches = Vec::new();
    for (status, expected) in PASSING_CENSUS {
        let actual = counts.get(status.as_str()).copied().unwrap_or(0);
        if actual != *expected {
            mismatches.push(format!("  {status}: census {expected}, actual {actual}"));
        }
    }
    assert!(
        mismatches.is_empty(),
        "the conformant-status census changed. A row moved between statuses without \
         changing the divergence budget — e.g. a projection that used to render a \
         string now errors. Decide whether that is an improvement or a regression, \
         then update PASSING_CENSUS deliberately (#1272).\n{}",
        mismatches.join("\n")
    );
}

/// Every row must be counted by exactly one of the two committed tables.
///
/// This is what makes the census *total*. Without it a status could be dropped
/// from both tables and silently stop being counted, which is the failure mode
/// #1272 is about — the guard has to notice its own coverage shrinking.
#[test]
fn every_status_is_counted_by_exactly_one_table() {
    let fx = load();
    let budgeted: usize = DIVERGENCE_BUDGET.iter().map(|(_, n)| n).sum();
    let passing: usize = PASSING_CENSUS.iter().map(|(_, n)| n).sum();
    assert_eq!(
        budgeted + passing,
        fx.rows.len(),
        "the two census tables cover {} rows but the enumeration has {}. Every row must \
         be counted exactly once: add the missing status to DIVERGENCE_BUDGET (if it \
         records a spec divergence) or PASSING_CENSUS (if it is conformant).",
        budgeted + passing,
        fx.rows.len()
    );

    for (status, _) in PASSING_CENSUS {
        assert!(
            KNOWN_PASSING_STATUSES.contains(status),
            "{status} is counted as conformant but is not in KNOWN_PASSING_STATUSES"
        );
    }
    for status in KNOWN_PASSING_STATUSES {
        assert!(
            PASSING_CENSUS.iter().any(|(s, _)| s == status),
            "{status} is allowlisted as conformant but PASSING_CENSUS does not count it"
        );
    }
}

/// Known passing/expected statuses that are deliberately **not** part of the
/// divergence budget: every row in one of these is a conformant outcome, not a
/// recorded spec divergence, so `divergence_budget_is_unchanged` rightly ignores
/// them. Keeping the list explicit is the point of issue #1107 — a *new* status
/// the generator starts emitting is then neither budgeted nor allowlisted, so
/// `every_observed_status_is_accounted_for` fails and names it, turning "a new
/// status category appeared" into a deliberate, reviewed decision instead of a
/// silent pass. Regenerate the roster of observed statuses with:
///   `cargo run --features dev --bin generate_spec_enumeration -- --census`
/// and read `by_status` in the generated fixture.
const KNOWN_PASSING_STATUSES: &[Status] = &[
    // Spec-forbidden string that ferro rejected outright (parse error) — the
    // conformant counterpart of the budgeted `false-acceptance`.
    Status::CorrectlyRejected,
    // Bad string repaired to the canonical form the spec names — the ideal
    // outcome of a repair.
    Status::Repaired,
    // Grammar example that parses into the coordinate axis the spec declares.
    Status::FormAxisOk,
    // Grammar example whose axis the spec does not state; ferro's parsed axis is
    // pinned as a baseline (not a spec matter). Emittable but currently 0 rows.
    Status::FormAxisPinned,
    // Projection that matches the form the spec states — the conformant
    // counterpart of the budgeted `projection-diverges`.
    Status::Preserved,
    // Projection with no spec-stated expectation, pinned as a baseline: ferro
    // rendered a form, declined/unavailable, or errored. None is a spec
    // divergence — the spec mandates no expectation for these.
    Status::ProjectionPinned,
    Status::ProjectionUnavailablePinned,
    Status::ProjectionErrorPinned,
    // Error-mode outcome pinned as ferro policy: the spec says nothing about
    // ferro's error modes, so a per-mode divergence is never a spec expectation.
    Status::ModeDivergencePinned,
];

/// A status is accounted for if it is either budgeted (a tracked divergence in
/// `DIVERGENCE_BUDGET`) or allowlisted (a known passing outcome in
/// `KNOWN_PASSING_STATUSES`).
fn status_is_accounted_for(status: &Status) -> bool {
    DIVERGENCE_BUDGET.iter().any(|(s, _)| s == status) || KNOWN_PASSING_STATUSES.contains(status)
}

/// Returns the distinct observed statuses that are neither budgeted (a tracked
/// spec divergence in `DIVERGENCE_BUDGET`) nor allowlisted (a known passing
/// outcome in `KNOWN_PASSING_STATUSES`), sorted for a stable message. Such a
/// status is a new, unreviewed category — `divergence_budget_is_unchanged` is
/// blind to it, so this is what `every_observed_status_is_accounted_for` gates
/// on.
fn unaccounted_statuses<'a>(statuses: impl IntoIterator<Item = &'a Status>) -> Vec<&'a Status> {
    let mut unknown: Vec<&Status> = statuses
        .into_iter()
        .filter(|status| !status_is_accounted_for(status))
        .collect();
    unknown.sort_unstable();
    unknown.dedup();
    unknown
}

#[test]
fn unaccounted_statuses_flags_a_new_status() {
    // A budgeted status is accounted for; a status neither budgeted nor
    // allowlisted (e.g. a future `projection panicked` defect indicator) must
    // be surfaced by name so it becomes a deliberate decision, not a silent
    // pass. This holds independent of the allowlist's contents.
    let observed = [
        Status::FalseAcceptance,
        Status::from_wire("projection panicked"),
    ];
    let expected = Status::Unknown("projection panicked".to_string());
    assert_eq!(unaccounted_statuses(&observed), vec![&expected]);
}

#[test]
fn budget_and_allowlist_are_disjoint() {
    // `DIVERGENCE_BUDGET` = tracked divergences, `KNOWN_PASSING_STATUSES` =
    // conformant outcomes: the two carve the status vocabulary into disjoint
    // halves. A status in both would blur that intent (and be double-classified),
    // so guard the split explicitly.
    let overlap: Vec<&Status> = KNOWN_PASSING_STATUSES
        .iter()
        .filter(|passing| {
            DIVERGENCE_BUDGET
                .iter()
                .any(|(budgeted, _)| budgeted == *passing)
        })
        .collect();
    assert!(
        overlap.is_empty(),
        "a status is in both DIVERGENCE_BUDGET and KNOWN_PASSING_STATUSES: {overlap:?}"
    );
}

#[test]
fn every_observed_status_is_accounted_for() {
    let fx = load();
    let unknown = unaccounted_statuses(fx.rows.iter().map(|row| &row.status));
    assert!(
        unknown.is_empty(),
        "the enumeration emitted status(es) that are neither in DIVERGENCE_BUDGET \
         nor KNOWN_PASSING_STATUSES: {unknown:?}. A new status category appeared — \
         classify it deliberately: budget it in DIVERGENCE_BUDGET if it records a \
         spec divergence, or add it to KNOWN_PASSING_STATUSES if it is a conformant \
         outcome. (`divergence_budget_is_unchanged` alone would not have caught this.)"
    );
}

/// The acceptance gate from the conformance design: the output-invariant
/// catalog must produce **zero** MUST-level violations over the outputs the
/// repo already blesses. Every MUST violation must trace back to a
/// `false-acceptance` row — a string the spec forbids that ferro accepted —
/// never to a blessed `preserved` output. A violation on a blessed output means
/// the *rule* is wrong, not ferro.
#[test]
fn invariant_catalog_has_no_false_positives_on_blessed_output() {
    let fx = load();
    let offenders: Vec<&Row> = fx
        .rows
        .iter()
        .filter(|r| r.status == Status::InvariantViolationMust)
        .filter(|r| {
            // The generator records the source row's pinned status in `note`;
            // a blessed output is one whose source row is `preserved`.
            parse_hgvs(&r.target).is_ok() && r.expected.as_deref() == Some("no violation")
        })
        .filter(|r| !r.id.contains("output-invariant/"))
        .collect();
    assert!(
        offenders.is_empty(),
        "invariant catalog fired outside the output-invariant dimension: {:?}",
        offenders.iter().map(|r| &r.id).collect::<Vec<_>>()
    );

    // Cross-check the same claim from the other side: no MUST-level invariant
    // row may cite a rule the design lists as advisory.
    for row in fx.rows.iter().filter(|r| r.dimension == "output-invariant") {
        if row.id.contains("/A1-") || row.id.contains("/A2-") {
            assert_eq!(
                row.normative_level,
                NormativeLevel::Should,
                "advisory rule {} must never be MUST-level",
                row.id
            );
        }
    }
}

/// Spec-mandated MUST rows that currently pass must keep passing. This is the
/// ratchet that turns the enumeration into a regression net rather than a
/// snapshot.
#[test]
fn passing_spec_mandated_musts_stay_passing() {
    let fx = load();
    let ctx = Replayer::new();
    let mut regressions = Vec::new();
    let mut checked = 0usize;

    for row in &fx.rows {
        if row.expectation != Expectation::SpecMandated
            || row.normative_level != NormativeLevel::Must
        {
            continue;
        }
        let passing = matches!(
            row.status,
            Status::CorrectlyRejected
                | Status::Repaired
                | Status::RejectedNotRepaired
                | Status::FormAxisOk
                // A projection that already matches the form the spec states.
                | Status::Preserved
        );
        if !passing {
            continue;
        }
        let Some(observed) = replay(row, &ctx) else {
            continue;
        };
        checked += 1;
        let ok = match row.status {
            Status::CorrectlyRejected | Status::RejectedNotRepaired => {
                observed.starts_with("parse error") || observed.starts_with("normalize error")
            }
            Status::Repaired => Some(&observed) == row.expected.as_ref(),
            Status::FormAxisOk => Some(&observed) == row.expected.as_ref(),
            // The spec writes a projected form with or without its accession;
            // an accession-less statement is compared against the coordinate
            // part only, mirroring the generator (see `projection_matches`).
            Status::Preserved => row
                .expected
                .as_ref()
                .is_some_and(|e| projection_matches(&observed, e)),
            _ => true,
        };
        if !ok {
            regressions.push(format!(
                "  {} ({}): expected {:?}, observed {:?}",
                row.id, row.spec_citation, row.expected, observed
            ));
        }
    }

    eprintln!("spec enumeration: {checked} passing spec-mandated MUST rows re-checked");
    assert!(
        regressions.is_empty(),
        "{} spec-mandated MUST assertion(s) regressed:\n{}",
        regressions.len(),
        regressions.join("\n")
    );
}

/// Projecting a result **back onto its own axis** must be a fixed point.
///
/// The one metamorphic relation in this file, and the only assertion here whose
/// oracle is neither the fixture nor a committed count: it compares ferro's
/// projection pipeline against *itself run twice*, which a regeneration cannot
/// reconcile. A projection that does not converge on its own axis is wrong
/// however it was recorded.
///
/// Restricted to same-axis rows on purpose. Re-projecting a *cross*-axis result
/// asks a different question and the answer is a legitimate decline — a `g.`
/// string carries no transcript frame to return through, so `project-c` of a
/// `project-g` result reports "not applicable" by design, not as a defect.
///
/// Two relations that look adjacent were measured and rejected rather than
/// left as follow-ons:
///
/// - **`r.` equals `c.` with `t`→`u`.** False. The axes legitimately differ:
///   `r.` drops intronic offsets, because a spliced RNA has no introns
///   (`c.94-23_188+33=` projects to `r.(94_188=)`), and its numbering base
///   differs (`c.3921del` → `r.(3922del)`). 15 of 271 pairs differ for those
///   reasons. Asserting the relation would pin a false invariant.
/// - **`c` → `g` → `c` round-trip identity.** Not expressible: the return leg
///   is the decline described above.
#[test]
fn a_projection_is_a_fixed_point_on_its_own_axis() {
    let fx = load();
    let ctx = Replayer::new();
    let Some(projector) = ctx.projector.as_ref() else {
        eprintln!("projection slice absent — skipping");
        return;
    };

    let mut divergent = Vec::new();
    let mut checked = 0usize;

    for row in &fx.rows {
        let Some(axis) = row.operation.strip_prefix("project-") else {
            continue;
        };
        let Some(code) = axis.chars().next() else {
            continue;
        };
        // Same-axis rows only: the input's own coordinate prefix must match.
        let core = row.target.rsplit(':').next().unwrap_or(&row.target);
        if !core.starts_with(&format!("{axis}.")) {
            continue;
        }
        let Some(once) = replay(row, &ctx) else {
            continue;
        };
        if !is_rendered(&once) {
            continue;
        }
        let Ok(reparsed) = parse_hgvs(&once) else {
            continue;
        };
        let twice = spec_projection::project_all_axes(projector, &reparsed)
            .axes
            .get(&code)
            .map(|r| r.as_observed());
        let Some(twice) = twice else { continue };
        if !is_rendered(&twice) {
            continue;
        }
        checked += 1;
        if twice != once {
            divergent.push(format!(
                "  {}\n    first pass : {once}\n    second pass: {twice}",
                row.id
            ));
        }
    }

    eprintln!("spec enumeration: {checked} same-axis projections re-projected");
    assert!(
        checked > 300,
        "re-projected too few rows ({checked}) — the same-axis filter or the \
         committed slice changed shape"
    );
    assert!(
        divergent.is_empty(),
        "{} projection(s) are not a fixed point on their own axis. Unlike the \
         replay test, nothing about this is recorded in the fixture — a \
         regeneration cannot make it pass (#1272).\n{}",
        divergent.len(),
        divergent.join("\n")
    );
}

/// Did the projector actually render a description, as opposed to declining?
///
/// `AxisResult`'s non-rendered arms stringify to a prose reason rather than an
/// HGVS string, and all three are legitimate outcomes here, so they are skipped
/// rather than compared.
fn is_rendered(observed: &str) -> bool {
    !observed.starts_with("error")
        && !observed.starts_with("unavailable")
        && !observed.starts_with("not applicable")
}
