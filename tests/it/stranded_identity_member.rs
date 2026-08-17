//! PINNED DEFECT: a *derived* identity member that no sibling overlaps, so
//! `drop_identity_members_covered_by_siblings` cannot reach it — leaving one
//! variant with two normal forms.
//!
//! # The shape
//!
//! ```text
//! reference  CGCGCGCGCAATCGCGCGCAATCGCG        10 = A, 11 = A, 12 = T
//!
//! in    TEMPLATE:g.[10_11insA;11del;12_13insAA]
//! out   TEMPLATE:g.[11=;12_13insAA]
//! but   TEMPLATE:g.12_13insAA        ->  TEMPLATE:g.12_13insAA
//! ```
//!
//! `10_11insA` inserts an `A` into the `AA` tract and `11del` takes one back
//! out, so between them they change nothing and collapse to `11=`. The third
//! member is untouched and sits two bases away. The result is a two-member
//! allele whose first member asserts nothing and whose second, spelled alone,
//! is already a fixed point — **two normalized forms for one variant**.
//!
//! # Why the existing rule does not fire, and why that is not a bug in it
//!
//! `merge::drop_identity_members_covered_by_siblings` (#1297, widened to
//! partial overlap by #1416) removes an identity that a sibling's span
//! *overlaps* under `blocks_sibling_shift`. Neither half holds here:
//!
//! - **No overlap.** `11=` and `12_13insAA` are disjoint. The rule scans every
//!   sibling, not just the adjacent one, so this is a reach limit rather than an
//!   adjacency bug — a claim worth stating precisely, because "it only looks at
//!   neighbours" is the natural guess and it is wrong.
//! - **No absorber left to overlap.** The member that cancelled `11=` is
//!   `10_11insA`, and it no longer exists as a member — it was consumed in the
//!   merge. The rule is written against a surviving absorber that grew *over*
//!   the identity's bases; here the pair annihilated each other and the residue
//!   was stranded beside an unrelated third member.
//! - Even had they overlapped, `12_13insAA` is an **insertion**, which both
//!   `claims_reference_bases` and `blocks_sibling_shift` deliberately exclude:
//!   "an insertion's span *is* the gap it occupies, so `g.264_265insA` nominally
//!   covers position 264 while changing nothing there."
//!
//! # Measured
//!
//! Over a 10,720-input three-member sweep (`ins` / middle / `ins` at sliding
//! anchors, both shuffle directions) on the reference above. The denominator is
//! `sweep_inputs()`'s own output — 134 surviving `(p, gap_a, gap_b)` triples ×
//! 4 payloads × 4 payloads × 5 middle spellings — and
//! [`every_stranded_identity_is_non_confluent_with_its_surviving_member`] pins
//! it exactly, so this number cannot drift away from the generator.
//!
//! **Two populations, with two different denominators.** The pinning test
//! passes `|members| members == 2` to [`census`], so it never sees an output
//! with three or more members — and on this corpus *every* self-cancel lives in
//! one. Reading a single "of those" chain across the two columns below
//! therefore chains figures the code never computes together; both are stated,
//! and both are pinned by test
//! ([`the_wider_multi_member_census_is_what_the_header_says`]):
//!
//! | identity members in… | outputs with **≥ 2** members | outputs with **exactly 2** — what the pinning test visits |
//! |---|---|---|
//! | total | 3,352 | **132** |
//! | of those, **self-cancels** (the input member is a no-op on its own) | 3,220 | **0** |
//! | of those, absorbed and *overlapping* a sibling (the rule fires) | 0 | 0 |
//! | of those, absorbed and **disjoint** — this defect | **132** | **132** |
//! | ... non-confluent with the surviving sibling spelled alone | 132 of 132 | **132 of 132** |
//!
//! The stranded split is 65 in the 5' direction and 67 in the 3', so this is
//! **not** a direction disagreement: both directions produce it, and both
//! produce the same string.
//!
//! ## The "overlapping" zero is ~85% structural
//!
//! It reads as "the existing rule already handles every overlapping case", and
//! for most of the population it means "the rule was never eligible to fire".
//! Of the 132 stranded rows the surviving sibling is an **insertion in 112**
//! and a **duplication in 20** (pinned below). `merge::blocks_sibling_shift` is
//! `claims_reference_bases` plus `Duplication`/`DupIns`, and
//! `claims_reference_bases` excludes `Insertion` **deliberately** — an
//! insertion's span *is* the gap it occupies. So for 112 of 132 (85%) no span
//! arrangement could have made the rule fire. Only the other 20 are a genuine
//! reach limit: a `dup` *does* block sibling shift, and there the rule is
//! stopped by disjointness alone.
//!
//! **What this generator can and cannot vary.** `sweep_inputs` hardcodes both
//! outer members as `ins` and varies only the anchor positions, the two
//! insertion payloads and the five middle spellings. The surviving sibling is
//! therefore always an insertion or a duplication derived from one; the corpus
//! cannot build a stranded identity beside a `del`, `delins` or `inv` sibling
//! at all. A zero on the "overlapping" row is a claim about this corpus, not
//! about the rule.
//!
//! ## `stranded` is a floor, not a count
//!
//! The census's `overlaps` predicate is purely positional — it compares the
//! rendered spans — so it reads `12_13insAA` as covering positions 12 and 13.
//! Ferro does not, for the reason directly above. The census therefore treats
//! *more* siblings as overlapping than the real rule does, and every row it
//! skips on that basis would in fact be stranded. Measured, that skip costs
//! **0** rows here (the "overlapping" row is empty on both columns), so 132 is
//! currently exact — but it is a lower bound by construction, and a generator
//! that produced an insertion flush against an identity would separate the two.
//!
//! ## The two-member sweep, committed rather than described
//!
//! [`two_member_sweep_inputs`] builds 900 inputs — the same anchors, payloads
//! and middle spellings with the second insertion dropped — and
//! [`a_two_member_sweep_strands_nothing_and_the_zero_is_structural`] pins the
//! result: **279** identity members, **all 279** self-cancels, **0** absorbed
//! and therefore **0** stranded. The zero is structural rather than reassuring:
//! with one sibling, an identity can only be the echo of an input member that
//! was already a no-op, and there is no third member left over for a residue to
//! be stranded beside. That is why the pinning sweep is three-member.
//!
//! (An earlier revision of this header quoted 5,694 identity members from an
//! uncommitted two-member sweep. No reconstruction reproduced that figure; it
//! is replaced by the committed generator above, whose numbers a reader can
//! re-run.)
//!
//! # What is NOT decided here
//!
//! Whether ferro should drop a derived identity that overlaps nothing. The
//! existing rule's doc states the opposing case deliberately: an identity
//! "says a position was examined and found unchanged. That is real information
//! when it stands alone, so `g.[1002=;1005del]` keeps both." That reasoning is
//! about an **authored** identity; nothing in the input here asserted `11=`.
//! Distinguishing the two is a representation decision with a real trailer, so
//! this file measures and pins rather than fixes.
//!
//! # Why this module is named in `SEQUENCE_ORACLE_EXCLUDE`
//!
//! The output above — an insertion at interbase 11|12 beside an identity
//! claiming 11..12 — **denotes no sequence**, which is #1281's shape and which
//! `FERRO_ASSERT_SEQUENCE` fires on. Two of this module's four tests do:
//! `every_stranded_identity_is_non_confluent_with_its_surviving_member` and
//! `the_wider_multi_member_census_is_what_the_header_says`. That is not a false
//! positive — the oracle is right, and this module pins the defect it is right
//! about, so a test that pins a defect and an oracle that aborts on it cannot
//! both run. Exactly the class `ci.yml`'s `ORACLE_EXCLUDE` documents.
//!
//! So when #1815 armed that flag in `test-oracle`, this module was named in
//! `ci.yml`'s `SEQUENCE_ORACLE_EXCLUDE` alongside the issue that retires the
//! row: **#1655**. It still runs there under the other three oracles, in a
//! second un-partitioned step of the same job, and it still runs in full in the
//! plain `test` job — so arming the fourth oracle cost this file nothing.
//!
//! **Fixing #1655 means deleting this module's row from that filter**; rehoming
//! the sweep the way `ORACLE_EXCLUDE`'s pinned-defect modules are homed is the
//! other way out. Do not silence the fire at the seam — that hollows out the
//! oracle, and this module's own measurement is the evidence the fire is
//! correct.

use std::collections::BTreeMap;

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

const REF: &str = "CGCGCGCGCAATCGCGCGCAATCGCG";

const DIRECTIONS: [ShuffleDirection; 2] =
    [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime];

fn provider() -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("TEMPLATE", REF.to_string());
    p
}

fn norm(input: &str, dir: ShuffleDirection) -> Option<String> {
    let parsed = parse_hgvs(input).ok()?;
    Normalizer::with_config(provider(), NormalizeConfig::default().with_direction(dir))
        .normalize(&parsed)
        .ok()
        .map(|v| v.to_string())
}

/// The `start..=end` a rendered member names, read off its text.
fn span_of(member: &str) -> Option<(u64, u64)> {
    let digits: String = member
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '_')
        .collect();
    let mut parts = digits.split('_').filter(|s| !s.is_empty());
    let start: u64 = parts.next()?.parse().ok()?;
    let end: u64 = parts.next().and_then(|s| s.parse().ok()).unwrap_or(start);
    Some((start, end))
}

/// The span `member` occupies once normalized **on its own**, but only when that
/// alone-normalization is an identity.
///
/// `None` when the member does not normalize, or normalizes to something that
/// changes the sequence — so a `Some` answer is precisely "this member is a
/// no-op, and here is where its no-op lands".
fn solo_identity_span(member: &str, dir: ShuffleDirection) -> Option<(u64, u64)> {
    let solo = norm(&format!("TEMPLATE:g.{member}"), dir)?;
    let rendered = solo.split_once("g.").map(|(_, m)| m)?;
    if !rendered.ends_with('=') {
        return None;
    }
    span_of(rendered)
}

fn members_of(rendered: &str) -> Vec<String> {
    let inner = rendered
        .split_once('[')
        .and_then(|(_, r)| r.rsplit_once(']'))
        .map(|(m, _)| m)
        .unwrap_or("");
    inner.split(';').map(str::to_string).collect()
}

/// The reported instance, spelled out, in both directions.
///
/// Kept separate from the census below so the defect has one exact string a
/// reader can run, and so a fix can be verified without reading a counter.
#[test]
fn the_reported_instance_strands_an_identity_beside_an_unrelated_member() {
    for dir in DIRECTIONS {
        assert_eq!(
            norm("TEMPLATE:g.[10_11insA;11del;12_13insAA]", dir).as_deref(),
            Some("TEMPLATE:g.[11=;12_13insAA]"),
            "PINNED DEFECT ({dir:?}) — `10_11insA` and `11del` annihilate, and the \
             `11=` residue is left beside a disjoint third member"
        );
        // The same variant, spelled as the one member that actually changes
        // anything, is already a fixed point. Two normal forms, one variant.
        assert_eq!(
            norm("TEMPLATE:g.12_13insAA", dir).as_deref(),
            Some("TEMPLATE:g.12_13insAA"),
            "the surviving member alone is a fixed point, so the pair above is \
             non-confluent"
        );
    }
}

/// The five middle spellings both sweeps use, at anchor `q`.
///
/// Shared so the two generators differ in exactly one axis — the presence of a
/// second insertion — which is what makes the two-member sweep a control for
/// the three-member one rather than a different experiment.
fn middle_spellings(q: u64) -> [String; 5] {
    [
        format!("{q}del"),
        format!("{q}_{}del", q + 1),
        format!("{q}dup"),
        format!("{q}_{}inv", q + 1),
        format!("{q}delinsA"),
    ]
}

/// The three-member family that produces the shape, swept.
///
/// A two-member generator cannot build it — every identity it produces is a
/// self-cancel — so a zero from one is structural. That is not asserted from
/// prose: [`two_member_sweep_inputs`] is committed beside this one and
/// [`a_two_member_sweep_strands_nothing_and_the_zero_is_structural`] measures
/// it. It is the reason this sweep is three-member.
fn sweep_inputs() -> Vec<String> {
    let mut out = Vec::new();
    let n = REF.len() as u64;
    for p in 5..(n - 6) {
        for gap_a in 1..=3u64 {
            for gap_b in 1..=3u64 {
                let q = p + gap_a;
                let r = q + gap_b;
                if r + 1 >= n {
                    continue;
                }
                for x in ["A", "AA", "AT", "T"] {
                    for y in ["A", "AA", "AT", "T"] {
                        for middle in middle_spellings(q) {
                            out.push(format!(
                                "TEMPLATE:g.[{p}_{}ins{x};{middle};{r}_{}ins{y}]",
                                p + 1,
                                r + 1
                            ));
                        }
                    }
                }
            }
        }
    }
    out
}

/// The same family with the **second insertion dropped** — the control.
///
/// Identical anchors, payloads and middle spellings to [`sweep_inputs`], one
/// member shorter. It exists so the module header's "a two-member sweep is
/// structurally blind" is a measurement a reader can re-run rather than a
/// recollection, which is what it used to be.
fn two_member_sweep_inputs() -> Vec<String> {
    let mut out = Vec::new();
    let n = REF.len() as u64;
    for p in 5..(n - 6) {
        for gap_a in 1..=3u64 {
            let q = p + gap_a;
            if q + 1 >= n {
                continue;
            }
            for x in ["A", "AA", "AT", "T"] {
                for middle in middle_spellings(q) {
                    out.push(format!("TEMPLATE:g.[{p}_{}ins{x};{middle}]", p + 1));
                }
            }
        }
    }
    out
}

/// The coarse edit kind a rendered member names.
///
/// Probed longest-first: `delins` contains `ins`, so a naive scan would call
/// every `delins` an insertion.
fn edit_kind(member: &str) -> &'static str {
    for kind in ["delins", "dup", "inv", "del", "ins"] {
        if member.contains(kind) {
            return kind;
        }
    }
    if member.ends_with('=') {
        "identity"
    } else {
        "other"
    }
}

/// Every identity member a sweep's outputs contain, classified.
///
/// Each field is a count over one traversal, so no two of them can be chained
/// across different populations by accident — which is exactly the mistake the
/// module header's census used to invite.
#[derive(Debug, Default)]
struct Census {
    /// Identity members seen in the visited outputs.
    identity_members: usize,
    /// Of those, honest echoes of an input member that is a no-op on its own.
    self_cancels: usize,
    /// Of those, absorbed and positionally overlapping a surviving sibling.
    /// The existing rule's territory — see the module header on why this being
    /// zero is largely structural.
    absorbed_overlapping: usize,
    /// Of those, absorbed and **disjoint** — this defect. One rendered row each.
    stranded: Vec<String>,
    /// Of `stranded`, how many are non-confluent with their surviving member
    /// spelled alone.
    non_confluent: usize,
    /// `stranded`, split per entry of [`DIRECTIONS`], in that order.
    stranded_by_direction: [usize; DIRECTIONS.len()],
    /// Of `stranded`, the surviving sibling's edit kind. Asserted to partition
    /// `stranded` by its callers, so a corpus that grew a third kind is a
    /// failure rather than a silently-dropped column.
    surviving_kinds: BTreeMap<&'static str, usize>,
}

/// Classify every identity member the outputs of `inputs` contain.
///
/// `visit` selects which outputs are examined, by member count — `|n| n == 2`
/// reproduces the gate the pinning test applies, `|n| n >= 2` the wider
/// population the module header also states.
fn census(inputs: &[String], visit: fn(usize) -> bool) -> Census {
    let mut out_census = Census::default();
    for (d, dir) in DIRECTIONS.into_iter().enumerate() {
        for input in inputs {
            let Some(out) = norm(input, dir) else {
                continue;
            };
            let out_members = members_of(&out);
            if !visit(out_members.len()) {
                continue;
            }
            let in_members = members_of(input);
            for (i, member) in out_members.iter().enumerate() {
                if !member.ends_with('=') {
                    continue;
                }
                let Some(id_span) = span_of(member) else {
                    continue;
                };
                out_census.identity_members += 1;
                // Self-cancel: *this* identity is the honest echo of an input
                // member that is already a no-op on its own, landing on the same
                // span under the same direction. Such an identity asserts
                // something the input asserted; it is not residue left behind by
                // a merge.
                //
                // The predicate has to name the identity under examination.
                // The unqualified form — "does *any* input member normalize
                // alone to an identity" — asks a different question, and would
                // discard a genuine residue whenever some unrelated member of
                // the same allele happened to be a no-op.
                //
                // That is latent rather than live today, and the distinction is
                // measured, not assumed: under the `|n| n == 2` gate the
                // unqualified form fires **0** times, and `stranded` is **132**
                // under the span-matched form, under the unqualified form, and
                // with no self-cancel check at all. The guard order is what
                // makes it latent — the gate runs first, and every
                // self-cancelling identity on this corpus lives in an output
                // with more than two members. A future change that let a
                // three-member input collapse to two while one input member
                // stayed a no-op would start dropping stranded rows silently
                // under the unqualified form, so it is written the narrow way.
                if in_members
                    .iter()
                    .any(|m| solo_identity_span(m, dir) == Some(id_span))
                {
                    out_census.self_cancels += 1;
                    continue;
                }
                let overlaps = out_members.iter().enumerate().any(|(j, s)| {
                    j != i && span_of(s).is_some_and(|(ss, se)| ss <= id_span.1 && se >= id_span.0)
                });
                if overlaps {
                    out_census.absorbed_overlapping += 1;
                    continue;
                }
                out_census
                    .stranded
                    .push(format!("{dir:?} {input} -> {out}"));
                out_census.stranded_by_direction[d] += 1;

                // A two-member output whose members are BOTH identities has no
                // surviving member to compare against, so it is neither
                // stranded nor not — the census would have to classify it
                // explicitly. Measured **0** times on either sweep, but nothing
                // structurally forbids it, so it panics with its own row rather
                // than through a bare `.unwrap()` whose message names only the
                // `Option`.
                let Some(sibling) = out_members.iter().find(|m| !m.ends_with('=')) else {
                    panic!(
                        "every member of {out:?} (from {input:?}) is an identity, so there \
                         is no surviving member to test confluence against. This shape was \
                         measured 0 times; classify it explicitly rather than widening the \
                         search."
                    );
                };
                *out_census
                    .surviving_kinds
                    .entry(edit_kind(sibling))
                    .or_default() += 1;
                if norm(&format!("TEMPLATE:g.{sibling}"), dir).as_deref() != Some(out.as_str()) {
                    out_census.non_confluent += 1;
                }
            }
        }
    }
    out_census
}

/// **Every** stranded identity this family produces is non-confluent with the
/// single-member spelling of what actually changed. That is the claim worth
/// pinning: the residue is never merely cosmetic.
///
/// The stranded set is asserted non-empty, so a generator change that stops
/// building the shape fails loudly rather than reporting a reassuring zero.
/// Emptiness has **two** causes, though, and the message says so: a fix and a
/// generator regression look identical from inside this test, and only
/// [`the_reported_instance_strands_an_identity_beside_an_unrelated_member`] —
/// which reads no generator — tells them apart.
///
/// Every figure it asserts is the **exactly-two-member** column of the module
/// header's table. The wider column is pinned separately by
/// [`the_wider_multi_member_census_is_what_the_header_says`], because the gate
/// below means this test never observes the population that column counts.
///
/// The input count is pinned **exactly** rather than as a floor, because the
/// module header quotes it as the denominator of every figure in its census. A
/// prose number and a `>` assertion cannot hold each other honest — the floor
/// passes over a wide band of generator edits, so the header would go stale
/// silently. Equality is safe here: `sweep_inputs` is a pure function of `REF`
/// and three constant loop bounds, with no ordering or environment dependence.
#[test]
fn every_stranded_identity_is_non_confluent_with_its_surviving_member() {
    let all = sweep_inputs();
    assert_eq!(
        all.len(),
        10_720,
        "the sweep must build exactly the family the module header describes; \
         if this moved, update the header's denominator in the same change"
    );

    // The gate is what makes the right-hand column of the header's table the
    // population this test measures. Every figure asserted below is that
    // column, and none of them may be read against the left-hand one.
    let c = census(&all, |members| members == 2);

    assert!(
        !c.stranded.is_empty(),
        "the sweep found no stranded identity at all. TWO causes produce this, \
         and they are opposites:\n\
         (a) THE DEFECT WAS FIXED — a rule now reaches a derived identity that \
             overlaps nothing. Delete this pinning and record the ruling.\n\
         (b) A GENERATOR CHANGE stopped building the shape, so the zero is \
             structural and proves nothing.\n\
         `the_reported_instance_strands_an_identity_beside_an_unrelated_member` \
         separates them: it is generator-independent, so it goes RED under (a) \
         and stays GREEN under (b)."
    );
    assert_eq!(
        c.identity_members, 132,
        "the header's right-hand column: identity members in two-member outputs"
    );
    assert_eq!(
        c.self_cancels, 0,
        "no self-cancel survives the two-member gate on this corpus; the \
         header's 3,220 belongs to the ≥2-member column and must not be chained \
         onto this one"
    );
    assert_eq!(
        c.absorbed_overlapping, 0,
        "the existing rule's territory is empty here — and see the module \
         header: 112 of the 132 carry an `ins` sibling, which \
         `blocks_sibling_shift` excludes by construction, so most of this zero \
         is structural rather than a statement about overlap"
    );
    assert_eq!(c.stranded.len(), 132, "PINNED DEFECT — the stranded set");
    assert_eq!(
        c.stranded_by_direction,
        [67, 65],
        "3' then 5', matching `DIRECTIONS` — both directions produce it, so \
         this is not a direction disagreement"
    );
    assert_eq!(
        c.surviving_kinds,
        BTreeMap::from([("ins", 112), ("dup", 20)]),
        "the evidence for the module header's structural reading of the \
         `absorbed_overlapping` zero; a third kind here means the generator \
         gained an axis and the header's 85% is stale"
    );
    assert_eq!(
        c.non_confluent,
        c.stranded.len(),
        "PINNED DEFECT — every stranded identity must (still) be non-confluent \
         with its surviving member spelled alone; {} of {} are",
        c.non_confluent,
        c.stranded.len()
    );
}

/// The wider population the module header's left-hand column states, pinned so
/// the prose cannot drift from the code.
///
/// It is deliberately a **separate** test rather than more assertions on the
/// one above: the two columns come from two traversals with two different
/// gates, and keeping them apart is what stops the next reader chaining "3,352
/// … of those 3,220 … of those 132" through a population the pinning test
/// never visits.
#[test]
fn the_wider_multi_member_census_is_what_the_header_says() {
    let c = census(&sweep_inputs(), |members| members >= 2);
    assert_eq!(
        c.identity_members, 3_352,
        "identity members, ≥2-member outputs"
    );
    assert_eq!(c.self_cancels, 3_220, "self-cancels, ≥2-member outputs");
    assert_eq!(c.absorbed_overlapping, 0, "absorbed and overlapping");
    assert_eq!(
        c.stranded.len(),
        132,
        "the stranded set is the same 132 either way — widening the gate adds \
         only self-cancels, which is the point"
    );
    assert_eq!(
        c.self_cancels + c.absorbed_overlapping + c.stranded.len(),
        c.identity_members,
        "the three classes must partition the identity members, or a row is \
         being silently dropped"
    );
}

/// The control: a two-member sweep strands nothing, and the zero is
/// **structural** — recorded as such rather than reported as reassurance.
///
/// With one sibling an identity can only be the echo of an input member that
/// was already a no-op, so there is never a third member for a residue to be
/// stranded beside. The assertions below say that in the only way that stays
/// honest: the identity-member count is asserted **non-zero** first, so a
/// generator that stopped producing identities at all cannot pass this test by
/// producing `0 of 0`.
#[test]
fn a_two_member_sweep_strands_nothing_and_the_zero_is_structural() {
    let all = two_member_sweep_inputs();
    assert_eq!(
        all.len(),
        900,
        "45 surviving `(p, gap_a)` pairs × 4 payloads × 5 middle spellings"
    );

    let c = census(&all, |members| members >= 2);
    assert_eq!(
        c.identity_members, 279,
        "the denominator, asserted before the zero so `0 of 0` cannot pass"
    );
    assert_eq!(
        c.self_cancels, c.identity_members,
        "every identity a two-member sweep produces is a self-cancel"
    );
    assert_eq!(c.absorbed_overlapping, 0, "nothing is absorbed");
    assert!(
        c.stranded.is_empty(),
        "a two-member sweep cannot build the shape; if this ever finds one, the \
         module header's claim that the zero is structural is wrong and the \
         three-member sweep is no longer the only generator that matters: {:?}",
        c.stranded
    );
}
