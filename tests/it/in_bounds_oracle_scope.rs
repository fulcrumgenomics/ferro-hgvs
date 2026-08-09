//! What `FERRO_ASSERT_IN_BOUNDS` can and cannot see, pinned against the corpus's
//! largest defect class.
//!
//! # The question this file settles
//!
//! `spec_conformance_axis` measures `outputs_leaving_the_transcript` at **371**
//! on the 3' direction. Run the same corpus with `FERRO_ASSERT_IN_BOUNDS=1`
//! armed and the oracle fires on **none** of them — only the idempotency oracle
//! fires, on 7 rows. Two readings were available and they have very different
//! consequences: either the in-bounds oracle is *blind* to the largest class the
//! corpus found, or the corpus counts something the oracle deliberately exempts.
//!
//! **It is the exemption, and it is documented.** The oracle asks one question —
//! *does the position this description names exist on its own sequence?* — and an
//! intronic offset names no position on that sequence at all. There is no length
//! for it to be compared against, so there is nothing for the predicate to
//! answer. `merge::first_out_of_bounds_coordinate`'s doc comment says so in as
//! many words, under "What it deliberately does not check":
//!
//! > **Intronic offsets, unknown, uncertain and special positions** (`pter`,
//! > `qter`) yield no plain axis position, so they are skipped — see
//! > [`readable_endpoints`] … An intronic offset is outside the transcript by
//! > construction and has no length to be compared with.
//!
//! # Where the skip happens, exactly
//!
//! `Normalizer::assert_in_bounds` (`src/normalize/mod.rs:1933`) →
//! `merge::first_out_of_bounds_coordinate` (`merge.rs:9231`) →
//! `member_out_of_bounds_coordinate` (`:9246`) → `readable_endpoints` (`:9290`) →
//! `simple_cds_pos` (`:1141`) / `simple_tx_pos` (`:1171`) / `simple_rna_pos`
//! (`:1201`). Each of those three readers opens with
//! `if pos.is_intronic() { return None; }` (`merge.rs:1146`, `:1176`, `:1206`),
//! and an endpoint that reads as `None` contributes no `(region, axis)` pair, so
//! the `find_map` never evaluates a comparison for it.
//!
//! # So the useful statement is a scope statement, not a bug report
//!
//! Both instruments are working; they measure different properties, and neither
//! subsumes the other:
//!
//! | | asks | answer for `NM_TEST.1:c.10+2del` |
//! |---|---|---|
//! | `FERRO_ASSERT_IN_BOUNDS` | does the named position exist on the served sequence? | unanswerable — skipped |
//! | corpus `outputs_leaving_the_transcript` | did the output name an intronic position its input did not? | yes — counted |
//!
//! **The consequence worth committing is the negative one:** a zero from
//! `FERRO_ASSERT_IN_BOUNDS` over the corpus is a claim about the oracle's scope,
//! not about the outputs, and the oracle must never be cited as coverage for the
//! transcript-exit class. That is the same shape as this repository's recorded
//! blindness family (#1456 / #1460 / #1478) — a zero that reads as safety when it
//! means the instrument could not build the question.
//!
//! Nor do the other two cover it: the intronic output is a fixed point, so
//! `FERRO_ASSERT_IDEMPOTENT` is silent, and `parse_hgvs` accepts an intronic
//! position on a bare `NM_` (pinned by
//! `spec_corpus_regressions::a_bare_transcript_accession_accepts_an_intronic_position`),
//! so `FERRO_ASSERT_REPARSE` is silent too. Both halves are measured in
//! `defect_371_transcript_exit::the_intronic_output_is_a_fixed_point_so_no_seam_oracle_sees_it`.
//! **No seam oracle covers the corpus's largest defect class**, and closing that
//! would mean a fourth predicate, not a widening of this one.
//!
//! # Why the predicate is not called directly here
//!
//! `first_out_of_bounds_coordinate` is `pub(crate)`, so an integration test
//! cannot call it. Its *positive* behaviour is already pinned by its own unit
//! tests in `src/normalize/merge.rs`
//! (`the_known_out_of_range_shapes_are_flagged`,
//! `a_reversed_circular_o_range_is_judged_only_on_its_endpoints`,
//! `an_authored_overrun_is_visible_through_a_special_position`). This file pins
//! the complementary half two ways: **structurally**, by asserting the output's
//! endpoints satisfy exactly the `is_intronic()` condition those readers
//! short-circuit on; and **behaviourally**, by normalizing the same input in a
//! child process with all three oracles armed and observing that it exits
//! cleanly.

use std::process::Command;

use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer, ShuffleDirection};

use crate::defect_371_transcript_exit::{junction_provider, CODING};
use ferro_hgvs::reference::transcript::Strand;

/// Environment variable the harness sets to tell a re-spawned copy of this test
/// binary that it is the child. Its absence is what makes the probe below a
/// no-op in an ordinary suite run.
const CHILD_MARKER: &str = "FERRO_371_ORACLE_CHILD";

/// The exit exemplar, as `(input, output)`.
///
/// The same shape `defect_371_transcript_exit` dissects: a four-base run flush
/// against exon 1's 3' end, with an intron that continues it by two.
const INPUT: &str = "c.7del";
const EXIT_OUTPUT: &str = "c.10+2del";

/// The fixture the exemplar is normalized against: exon 1 ends `…AAAA`, the
/// intron opens `AA`, everything else is `C`.
fn exit_provider(strand: Strand) -> ferro_hgvs::MockProvider {
    let filler: String = std::iter::repeat_n('C', 20).collect();
    let mut exon1 = filler.clone();
    exon1.replace_range(16..20, "AAAA");
    let mut intron = String::from("AA");
    while intron.len() < 30 {
        intron.push(if intron.len() % 2 == 0 { 'C' } else { 'G' });
    }
    junction_provider(
        strand,
        [exon1.as_str(), filler.as_str(), filler.as_str()],
        &intron,
        true,
    )
}

fn normalize_3prime(provider: &ferro_hgvs::MockProvider, input: &str) -> String {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("{input} must parse: {e}"));
    Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    )
    .normalize(&variant)
    .unwrap_or_else(|e| panic!("{input} must normalize: {e}"))
    .to_string()
}

/// Both endpoints of a parsed `c.` description, as `(is_intronic, offset)`.
fn cds_endpoint_offsets(descriptor: &str) -> Vec<(bool, Option<i64>)> {
    let parsed = parse_hgvs(descriptor).expect("the descriptor parses");
    let HgvsVariant::Cds(cds) = parsed else {
        panic!("{descriptor} is not a `c.` description");
    };
    [&cds.loc_edit.location.start, &cds.loc_edit.location.end]
        .into_iter()
        .filter_map(|boundary| boundary.inner())
        .map(|pos| (pos.is_intronic(), pos.offset))
        .collect()
}

// ---------------------------------------------------------------------------
// The structural pin: the endpoint the oracle would read is skipped by name
// ---------------------------------------------------------------------------

/// **Question.** Does the transcript-exit output name a coordinate
/// `first_out_of_bounds_coordinate` can read?
///
/// **No.** Every endpoint it carries is intronic, which is precisely the
/// condition `simple_cds_pos` short-circuits on at `merge.rs:1146` before any
/// bound is consulted. The oracle's silence therefore carries no information
/// about the output: it is the silence of a predicate that was handed nothing to
/// compare.
///
/// **Authority.** `src/normalize/merge.rs:9216-9220`, the "What it deliberately
/// does not check" list on `first_out_of_bounds_coordinate`: "*Intronic offsets,
/// unknown, uncertain and special positions (`pter`, `qter`) yield no plain axis
/// position, so they are skipped … An intronic offset is outside the transcript
/// by construction and has no length to be compared with.*" This is an
/// **adjudicated-deviation**-style record in the `KNOWN_DIVERGENT_INPUTS`
/// tradition: it pins the exemption *as* an exemption, so that a later widening
/// of the oracle has to move this assertion rather than silently subsuming it.
#[test]
fn the_transcript_exit_output_carries_only_endpoints_the_oracle_skips() {
    for strand in [Strand::Plus, Strand::Minus] {
        let output = normalize_3prime(&exit_provider(strand), &format!("{CODING}:{INPUT}"));
        assert_eq!(
            output,
            format!("{CODING}:{EXIT_OUTPUT}"),
            "the fixture must still produce the exit for this test to mean anything ({strand:?})"
        );

        let endpoints = cds_endpoint_offsets(&output);
        assert!(
            !endpoints.is_empty(),
            "VACUOUS — the output carries no readable endpoint at all"
        );
        for (intronic, offset) in &endpoints {
            assert!(
                *intronic,
                "every endpoint of `{output}` must be intronic, which is the condition \
                 `simple_cds_pos` (merge.rs:1146) returns `None` on; got offset {offset:?}"
            );
        }
    }
}

/// **Question.** Is the exemption *needed*, or could the oracle compare the
/// intronic coordinate against the transcript length and be done?
///
/// **It is needed, and the arithmetic shows why.** The input's endpoint sits at
/// `c.10` — transcript position 20 of a 60-base transcript — and the output's at
/// `c.10+2`, which is *the same* transcript position carrying an offset of 2.
/// There is no third number: the offset counts genomic bases the transcript does
/// not contain, so `10 + 2` is not a coordinate on any sequence the provider
/// serves and comparing it against 60 would be comparing two different axes.
///
/// That is the whole content of "has no length to be compared with", stated as an
/// assertion rather than as prose — and it is why closing this gap needs a fourth
/// predicate ("did the output leave the transcript?") rather than a wider bound
/// on this one.
#[test]
fn an_intronic_offset_has_no_axis_position_to_compare_against_a_length() {
    let provider = exit_provider(Strand::Plus);
    let output = normalize_3prime(&provider, &format!("{CODING}:{INPUT}"));

    // The input is plainly on-axis; the output is plainly not.
    for (intronic, offset) in cds_endpoint_offsets(&format!("{CODING}:{INPUT}")) {
        assert!(
            !intronic,
            "the input's endpoint is exonic, offset {offset:?}"
        );
    }
    let exit = cds_endpoint_offsets(&output);
    assert_eq!(
        exit.iter().map(|(_, offset)| *offset).collect::<Vec<_>>(),
        vec![Some(2), Some(2)],
        "the exit is `c.10+2`, i.e. transcript position 20 plus an offset of 2 — and the \
         offset counts bases the transcript does not contain"
    );

    // And the base the offset hangs off is comfortably inside the transcript, so
    // no bound the oracle could tighten would ever reach this row.
    use ferro_hgvs::reference::ReferenceProvider;
    let length = provider
        .get_sequence_length(CODING)
        .expect("the transcript reports a length");
    assert_eq!(length, 60, "the fixture transcript is 60 bases");
}

// ---------------------------------------------------------------------------
// The behavioural pin: armed, and silent
// ---------------------------------------------------------------------------

/// The child-side probe. A no-op unless [`CHILD_MARKER`] is set, so an ordinary
/// suite run pays nothing and the ordinary run cannot be what proves the point.
///
/// When it *is* the child it first checks that the arming actually arrived —
/// without that check, "the child exited 0" would be indistinguishable from "the
/// environment never reached the child", which is the flattering-zero failure
/// mode this repository keeps recording.
#[test]
fn oracle_probe_normalizes_the_transcript_exit_under_arming() {
    if std::env::var_os(CHILD_MARKER).is_none() {
        return;
    }
    for flag in [
        "FERRO_ASSERT_IN_BOUNDS",
        "FERRO_ASSERT_IDEMPOTENT",
        "FERRO_ASSERT_REPARSE",
    ] {
        assert_eq!(
            std::env::var(flag).ok().as_deref(),
            Some("1"),
            "{flag} must be armed in the child, or its silence proves nothing"
        );
    }
    for strand in [Strand::Plus, Strand::Minus] {
        let output = normalize_3prime(&exit_provider(strand), &format!("{CODING}:{INPUT}"));
        assert_eq!(output, format!("{CODING}:{EXIT_OUTPUT}"));
        // Re-normalizing is what makes the idempotency oracle's silence a
        // measurement rather than an absence of opportunity.
        assert_eq!(normalize_3prime(&exit_provider(strand), &output), output);
    }
}

/// Run [`oracle_probe_normalizes_the_transcript_exit_under_arming`] in a fresh
/// copy of this test binary, with `armed` deciding whether the three seam oracle
/// flags are set. Returns whether the child exited successfully.
fn run_probe_child(armed: bool) -> bool {
    let mut command = Command::new(std::env::current_exe().expect("the test binary has a path"));
    command
        .args([
            "--exact",
            "in_bounds_oracle_scope::oracle_probe_normalizes_the_transcript_exit_under_arming",
            "--nocapture",
            "--test-threads",
            "1",
        ])
        .env(CHILD_MARKER, "1");
    for flag in [
        "FERRO_ASSERT_IN_BOUNDS",
        "FERRO_ASSERT_IDEMPOTENT",
        "FERRO_ASSERT_REPARSE",
    ] {
        if armed {
            command.env(flag, "1");
        } else {
            command.env_remove(flag);
        }
    }
    command.status().expect("the test binary re-runs").success()
}

/// **Question.** With all three seam oracles armed, does normalizing the
/// transcript-exit exemplar fire any of them?
///
/// **No — the child exits cleanly.** Which, on its own, would be worth nothing:
/// an oracle that was never armed is silent for free. So this is stated as a pair
/// with the control below, which runs the identical child *without* the flags and
/// requires it to fail. The two together say the harness can tell armed from
/// unarmed, and that armed is silent here.
///
/// What the pair does **not** claim is that the oracle machinery could fire on
/// anything at all — that is pinned where it belongs, on the predicate's own unit
/// tests in `src/normalize/merge.rs`. The reason it cannot fire *here* is the
/// structural one above, which is the stronger of the two statements.
#[test]
fn the_in_bounds_oracle_does_not_fire_on_an_output_that_left_the_transcript() {
    assert!(
        run_probe_child(true),
        "the armed child failed. If it panicked with FERRO_ASSERT_IN_BOUNDS, the oracle's \
         scope has WIDENED to cover the transcript-exit class — delete this file's exemption \
         pins and re-bless `spec_conformance_axis`'s census in the same commit."
    );
}

/// The control for the test above: the same child, unarmed, must fail.
///
/// Its assertion is the arming check inside the probe, so a failure here means
/// the environment is not reaching the child — and therefore that the armed run's
/// success is vacuous. Written as a separate test so the two outcomes are
/// distinguishable in a CI log rather than collapsed into one message.
#[test]
fn the_unarmed_control_child_fails_so_the_armed_run_is_not_vacuous() {
    assert!(
        !run_probe_child(false),
        "the unarmed child SUCCEEDED, so the probe cannot tell armed from unarmed and the \
         companion test proves nothing"
    );
}
