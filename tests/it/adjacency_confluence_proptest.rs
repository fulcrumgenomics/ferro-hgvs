//! Randomised confluence and soundness over cis alleles that contain an
//! insertion next to another member.
//!
//! # Relationship to `cis_allele_confluence_proptest`
//!
//! That module is **not** touched by this one, deliberately. Its checked-in
//! seeds in `tests/proptest-regressions/cis_allele_confluence_proptest.txt` are
//! RNG seeds, not values: they are replayed *through* its strategies, so editing
//! how one of those strategies consumes randomness regenerates different cases
//! and silently retires the guards those seeds represent. Every strategy here is
//! new and lives in this file, so this module's seeds and that module's seeds
//! cannot interfere.
//!
//! # What is generated, and why each knob exists
//!
//! A generator that cannot vary the property a finding keys on reports a
//! structural zero — a clean run that proves nothing. `dump_normalized_corpus`
//! was blind three separate ways for exactly this reason (member geometry,
//! scale, transcript geometry), each invisible until the previous was fixed. So
//! the knobs below are chosen against the specific things this family has
//! broken on, and [`generator_reaches_every_arrangement`] asserts the generator
//! still reaches each of them rather than trusting that it does:
//!
//! - **junction placement** relative to the sibling span: before, at the 5'
//!   edge, strictly interior, at the 3' edge, after. The edge/interior boundary
//!   is the whole question.
//! - **payload against the flanking bases**: a payload equal to its neighbour
//!   3'-shifts, which can carry a junction that started at an edge *into* a
//!   sibling's span. A generator drawing payloads independently of the reference
//!   almost never produces that.
//! - **sibling edit kind**, including the two that keep their bases (`inv`,
//!   `delins`), the one that removes them (`del`, the #1406 exemption), and the
//!   one that writes at a junction of its own (`dup`).
//! - **separation** 0, 1 and 2, since the merge rules change at 1.
//! - **member count** 2 and 3, since order defects need three to distinguish
//!   from a stable sort.
//!
//! # The properties
//!
//! Soundness (an accepted output denotes the input's sequence) is checked
//! against `cis_apply_oracle::apply_with`, which converts each member to SPDI
//! and splices the reference itself rather than going through the normalizer.
//! Checking the normalizer against itself would pass for any self-consistent
//! defect.
//!
//! # …and the soundness half is currently OFF — read this before citing the
//! module as coverage
//!
//! [`an_accepted_allele_denotes_its_input_and_reparses`] is the only property
//! here that consults an oracle outside the normalizer, and it is `#[ignore]`d
//! against **#2013**. So as shipped this module runs three *self-consistency*
//! properties — order-independence, idempotence and the edge-adjacency verdict
//! — and **zero** soundness ones. The paragraph above describes the check the
//! module is built around, not a check that executes today.
//!
//! The ignore is pre-existing behaviour, not something this module introduced:
//! an insertion with an adjacent `dup` and a nearby `del` normalizes to a
//! different sequence, and it reproduces byte-for-byte on `origin/main`.
//! Un-ignoring the property is #2013's acceptance criterion, and doing so is
//! what restores the guarantee this section opens with.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::{parse_hgvs, HgvsVariant, NormalizeConfig, Normalizer};
use proptest::prelude::*;
use proptest::test_runner::Config as ProptestConfig;

use crate::common::cis_apply_oracle::apply_with;
use crate::common::hg38_window::{base_at, local_desc, provider, HG38_WINDOW, LOCAL_CONTIG};

/// Local positions the generator may touch.
///
/// Kept well inside the window so the normalizer's 100-base shift window never
/// runs off an end — an edit that hits the boundary would fail for a reason
/// that has nothing to do with adjacency.
const LO: u64 = 200;
const HI: u64 = 400;

/// Where the insertion's junction sits relative to the sibling span `[s, e]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Placement {
    /// Junction strictly 5' of the span.
    Before,
    /// Junction at `s - 1`: abutting the span's 5' edge.
    FivePrimeEdge,
    /// Junction strictly inside the span.
    Interior,
    /// Junction at `e`: abutting the span's 3' edge, and where a `dup` writes.
    ThreePrimeEdge,
    /// Junction strictly 3' of the span.
    After,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum SiblingKind {
    Sub,
    Del,
    Inv,
    Delins,
    Dup,
}

/// One generated allele, kept structured so the properties can ask about its
/// geometry rather than re-parsing the string they are testing.
#[derive(Debug, Clone)]
struct AdjacencyCase {
    /// Insertion junction: between `gap` and `gap + 1`.
    gap: u64,
    payload: String,
    span_start: u64,
    span_end: u64,
    kind: SiblingKind,
    placement: Placement,
    /// An optional third member, 3' of everything else.
    third: Option<(u64, u64)>,
}

impl AdjacencyCase {
    fn sibling_member(&self) -> String {
        let (s, e) = (self.span_start, self.span_end);
        match self.kind {
            SiblingKind::Sub => {
                let alt = match base_at(s) {
                    b'A' => 'C',
                    b'C' => 'A',
                    b'G' => 'T',
                    _ => 'G',
                };
                format!("{s}{}>{alt}", base_at(s) as char)
            }
            SiblingKind::Del if s == e => format!("{s}del"),
            SiblingKind::Del => format!("{s}_{e}del"),
            SiblingKind::Inv => format!("{s}_{e}inv"),
            SiblingKind::Delins if s == e => format!("{s}delinsAA"),
            SiblingKind::Delins => format!("{s}_{e}delinsAA"),
            SiblingKind::Dup if s == e => format!("{s}dup"),
            SiblingKind::Dup => format!("{s}_{e}dup"),
        }
    }

    /// Members in ascending coordinate order, which is the order an allele is
    /// conventionally written in.
    fn members(&self) -> Vec<String> {
        let ins = format!("{}_{}ins{}", self.gap, self.gap + 1, self.payload);
        let sib = self.sibling_member();
        let mut ordered: Vec<(u64, String)> = vec![(self.gap, ins), (self.span_start, sib)];
        if let Some((s, e)) = self.third {
            ordered.push((
                s,
                if s == e {
                    format!("{s}del")
                } else {
                    format!("{s}_{e}del")
                },
            ));
        }
        ordered.sort_by_key(|(k, _)| *k);
        ordered.into_iter().map(|(_, m)| m).collect()
    }

    fn body(&self) -> String {
        format!("[{}]", self.members().join(";"))
    }
}

/// A payload that is sometimes a copy of the base 3' of the junction, so that
/// 3'-shifting actually happens for a meaningful fraction of cases.
fn payload_strategy(gap: u64) -> impl Strategy<Value = String> {
    let neighbour = (base_at(gap + 1) as char).to_string();
    prop_oneof![
        2 => Just(neighbour.clone()),
        2 => Just(format!("{neighbour}{neighbour}")),
        3 => "[ACGT]{1,4}".prop_map(|s| s),
    ]
}

fn case_strategy() -> impl Strategy<Value = AdjacencyCase> {
    (
        LO..HI,
        1u64..7,
        prop_oneof![
            Just(SiblingKind::Sub),
            Just(SiblingKind::Del),
            Just(SiblingKind::Inv),
            Just(SiblingKind::Delins),
            Just(SiblingKind::Dup),
        ],
        prop_oneof![
            Just(Placement::Before),
            Just(Placement::FivePrimeEdge),
            Just(Placement::Interior),
            Just(Placement::ThreePrimeEdge),
            Just(Placement::After),
        ],
        0u64..3,
        prop::option::of(1u64..4),
    )
        .prop_flat_map(
            move |(span_start, len, kind, placement, separation, third_len)| {
                // A substitution is one base by construction; everything else may
                // be a range. An `inv` needs at least two bases to mean anything.
                let len = match kind {
                    SiblingKind::Sub => 1,
                    SiblingKind::Inv => len.max(2),
                    _ => len,
                };
                let span_end = span_start + len - 1;
                let gap = match placement {
                    Placement::Before => span_start.saturating_sub(separation + 2),
                    Placement::FivePrimeEdge => span_start - 1,
                    // Only meaningful when the span has an interior junction.
                    Placement::Interior => span_start + (len / 2).max(1) - 1,
                    Placement::ThreePrimeEdge => span_end,
                    Placement::After => span_end + separation + 1,
                };
                let third = third_len.map(|l| {
                    let s = span_end.max(gap) + 4;
                    (s, s + l - 1)
                });
                payload_strategy(gap).prop_map(move |payload| AdjacencyCase {
                    gap,
                    payload,
                    span_start,
                    span_end,
                    kind,
                    placement,
                    third,
                })
            },
        )
        .prop_filter("stay inside the window", |c| {
            c.gap >= LO
                && c.span_start >= LO
                && c.span_end < HI
                && c.third.map(|(_, e)| e < HI).unwrap_or(true)
                && c.gap + 1 < HI
        })
}

fn normalize(body: &str) -> Result<String, String> {
    let input = local_desc(body);
    let variant: HgvsVariant = parse_hgvs(&input).map_err(|e| format!("parse: {e}"))?;
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    .map(|v| v.to_string())
    .map_err(|e| format!("{e}"))
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 768, ..ProptestConfig::default() })]

    /// An accepted allele's normalized form denotes the same sequence as its
    /// input, and re-parses.
    ///
    /// Both halves matter and they fail differently: a sequence change is a
    /// wrong answer, an unparseable output is an unusable one. The seam defect
    /// in `insertion_adjacency_defects` is the second kind, and would pass a
    /// sequence-only property.
    #[test]
    fn an_accepted_allele_denotes_its_input_and_reparses(case in case_strategy()) {
        let body = case.body();
        let Ok(out) = normalize(&body) else { return Ok(()); };

        prop_assert!(
            parse_hgvs(&out).is_ok(),
            "{body} normalized to {out}, which does not re-parse"
        );

        let p = provider();
        let before = apply_with(&p, HG38_WINDOW, &local_desc(&body));
        let after = apply_with(&p, HG38_WINDOW, &out);
        // `None` before means the input denotes nothing (an overlapping
        // description); the normalizer is entitled to refuse or repair those, so
        // only cases the oracle can read are compared.
        if let (Some(b), Some(a)) = (before, after) {
            prop_assert_eq!(&b, &a, "{} changed the sequence: {} -> {}", body, b, a);
        }
    }

    /// Member order changes neither the verdict nor the normalized form.
    ///
    /// Compares the *verdict* (accepted or refused) and, when accepted, the
    /// emitted string — not the refusal text. The diagnostic wording does vary
    /// with authoring order, which is a real but separate and much smaller
    /// problem; it is pinned deterministically as
    /// `insertion_adjacency_defects::defect_overlap_diagnostic_depends_on_member_order`
    /// rather than being allowed to mask this property, whose failures are the
    /// serious kind.
    #[test]
    fn member_order_is_not_part_of_what_an_allele_denotes(case in case_strategy()) {
        let members = case.members();
        let forward = format!("[{}]", members.join(";"));
        let reversed = {
            let mut m = members.clone();
            m.reverse();
            format!("[{}]", m.join(";"))
        };
        let f = normalize(&forward);
        let r = normalize(&reversed);
        prop_assert_eq!(
            f.is_ok(),
            r.is_ok(),
            "verdict depends on member order for {}: {:?} vs {:?}",
            forward,
            f,
            r
        );
        if let (Ok(a), Ok(b)) = (&f, &r) {
            prop_assert_eq!(a, b, "member order leaked into the form for {}", forward);
        }
    }

    /// Normalizing twice is normalizing once.
    #[test]
    fn normalization_is_idempotent(case in case_strategy()) {
        let body = case.body();
        let Ok(once) = normalize(&body) else { return Ok(()); };
        let stripped = once
            .strip_prefix(&format!("{LOCAL_CONTIG}:g."))
            .unwrap_or(&once)
            .to_string();
        let Ok(twice) = normalize(&stripped) else {
            return Err(proptest::test_runner::TestCaseError::fail(format!(
                "{body} normalized to {once}, which then failed to normalize"
            )));
        };
        prop_assert_eq!(&once, &twice, "{} is not idempotent", body);
    }

    /// A junction at a span's **edge** is never refused as an overlap.
    ///
    /// The direct statement of the reported gap, as a property: edge-adjacency
    /// is not overlap, whatever the sibling's kind — and with **no** exception,
    /// which is what this change makes true.
    ///
    /// # Both edges, including the one a `dup` writes at
    ///
    /// `ThreePrimeEdge` is in the set deliberately, and it is the geometry this
    /// change is about. It used to be absent, which made the whole property
    /// silent on the 3' edge for *every* kind — `sub`, `del`, `inv` and `delins`
    /// as well as `dup` — and left the carve-out it was paired with unreachable:
    /// with `ThreePrimeEdge` outside `edge`, `!edge` was already true there, so
    /// a `dup_write_collision` disjunct could never change the outcome. A
    /// documented exception that cannot fire is not an exception, it is dead
    /// code standing where coverage is supposed to be.
    ///
    /// The exception is gone rather than repaired, because the ruling this
    /// change ships removes the thing it described: a `dup` sharing its write
    /// junction with one insertion is no longer "a genuine collision under
    /// current policy", it is the accepted case
    /// (`insertion_adjacency_corpus::insertion_at_a_duplications_write_junction_composes_with_it`).
    /// Two junction writers of the *same* kind still conflict, but the generator
    /// emits one insertion beside one span edit, so it cannot build that shape
    /// and no carve-out is needed for it here.
    #[test]
    fn an_edge_adjacent_junction_is_not_treated_as_an_overlap(case in case_strategy()) {
        let edge = matches!(
            case.placement,
            Placement::Before
                | Placement::FivePrimeEdge
                | Placement::ThreePrimeEdge
                | Placement::After
        );
        if !edge {
            return Ok(());
        }
        let body = case.body();
        if let Err(e) = normalize(&body) {
            prop_assert!(
                !(e.contains("coincident cis-allele edits") || e.contains("W5002")),
                "{} was refused as an overlap, but its junction is at an edge: {}",
                body,
                e
            );
        }
    }
}

/// The generator actually reaches every arrangement it claims to.
///
/// Without this, narrowing a strategy later turns the properties above into a
/// structural zero — they would still pass, over a case space that no longer
/// contains the interesting geometry. Asserting reachability is cheap; finding
/// out the hard way is not.
#[test]
fn generator_reaches_every_arrangement() {
    use proptest::strategy::ValueTree;
    use proptest::test_runner::TestRunner;
    use std::collections::HashSet;

    let mut runner = TestRunner::deterministic();
    let strategy = case_strategy();

    let mut placements: HashSet<Placement> = HashSet::new();
    let mut kinds: HashSet<SiblingKind> = HashSet::new();
    let mut three_member = false;
    let mut shifting_payload = false;

    for _ in 0..4000 {
        let case = strategy
            .new_tree(&mut runner)
            .expect("strategy produces a value")
            .current();
        placements.insert(case.placement);
        kinds.insert(case.kind);
        three_member |= case.third.is_some();
        shifting_payload |= case.payload.as_bytes()[0] == base_at(case.gap + 1);
    }

    assert_eq!(
        placements.len(),
        5,
        "not every junction placement is reached: {placements:?}"
    );
    assert_eq!(
        kinds.len(),
        5,
        "not every sibling kind is reached: {kinds:?}"
    );
    assert!(three_member, "three-member alleles are never generated");
    assert!(
        shifting_payload,
        "no payload ever matches its 3' neighbour, so no case can 3'-shift into a sibling"
    );
}
