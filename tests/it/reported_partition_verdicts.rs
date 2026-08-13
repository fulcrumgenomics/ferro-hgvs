//! Per-spelling verdicts for the externally reported partition family
//! (#1419, #1420, #1421).
//!
//! # What these rows are, and why they were not already runnable
//!
//! Three issues report the same defect from a downstream pipeline: one variant,
//! written two ways, normalizes to two different strings. Each issue body states
//! four things about each of its rows, in prose:
//!
//! 1. what ferro prints for **each** of the two spellings,
//! 2. which of the two the HGVS recommendations make canonical, and on which
//!    line,
//! 3. that `EquivalenceChecker` still calls the two one variant — at
//!    `SequenceMatch` where the two spellings still print different strings,
//!    and at `NormalizedMatch` for the three `1419-r*` pairs that #1649
//!    converged. Either way it is what makes the divergence a *representation*
//!    problem rather than a correctness one, and what leaves consumers a
//!    working fallback,
//! 4. the reference bases the argument in (2) rests on.
//!
//! Only (1)'s *pairwise disagreement* was guarded. `reported_confluence_pairs`
//! pins that the two spellings must converge, and deliberately declines to pin
//! *which* string they converge on — pinning a winner there would freeze an open
//! product decision (#1235) by accident. That is the right call for a ratchet on
//! a future fix, and it leaves the four claims above untested: today's answers,
//! the spec argument for them, the equivalence fallback, and the reference.
//!
//! This module tests those. It is the *present-tense* half — what ferro does
//! now and why that is or is not what the recommendations ask for — where
//! `reported_confluence_pairs` is the *future-tense* half, ratcheting the fix.
//! Neither subsumes the other and they assert disjoint things.
//!
//! # Convention: characterization, not `#[ignore]`
//!
//! Twelve of the eighteen rows below print something the cited recommendation
//! argues against ([`OPEN_GAPS`], which records why the count rose from nine).
//! They are pinned **as they behave today**, with the wanted
//! answer written into the row, rather than being expressed as an `#[ignore]`d
//! test of the wanted answer.
//!
//! That is deliberate and it is the whole point. An `#[ignore]`d test never
//! runs, so it cannot notice the day the answer moves — and for this family the
//! *move* is the event that matters. A normalized string is the key a downstream
//! consumer stores read counts against, so any change to one of these rows is a
//! migration for somebody, whether it is a fix or a regression. Pinned, it
//! surfaces in a diff and has to be argued for in review. Ignored, it happens
//! silently and is discovered downstream.
//!
//! [`OPEN_GAPS`] counts the disagreeing rows so the count cannot drift away from
//! the table: closing one means editing a row *and* the census, deliberately.
//!
//! # The pair model has three states, and all three are sayable
//!
//! A reported pair reaches its [`Row::wanted`] form from **both** spellings,
//! from **one**, or from **neither**. That is [`PairState`], census-pinned per
//! pair in [`PAIR_STATES`] and measured by
//! [`each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names`].
//!
//! It used to have two, and the missing one was the state that means *fixed*.
//! The old `each_pair_reaches_its_canonical_form_from_exactly_one_spelling`
//! asserted `a.verdict != b.verdict` — "a pair is one canonical spelling and one
//! gap" — which is `OneReaches` mistaken for a structural fact, and a
//! `PAIRS_NO_SPELLING_REACHES` exemption list was bolted on beside it to carry
//! the pairs where both rows are `Gap`. `BothReach` had no expression at all,
//! so the remedy that test's own failure message prescribed — *"flip this row's
//! verdict to Canonical, lower OPEN_GAPS"* — fired the `assert_ne!` twenty lines
//! above it. Two independent changes were blocked on exactly that, and neither
//! was a fix to this module: #1616's half of it, and an `authority`-only relabel
//! of the `1421-n*` rows (#1802).
//!
//! **Making the third state sayable is deliberately not an auto-ratchet.** A
//! pair changing state is still a **four-place** edit — this pair's
//! [`PAIR_STATES`] entry, [`PAIR_STATE_CENSUS`], the affected rows' `verdict`
//! fields and [`OPEN_GAPS`] — that has to name the clause that carried it, for
//! the reason [`OPEN_GAPS`] gives at length: a row can reach its `wanted`
//! string by *re-derivation* rather than by a licensed fix. Nothing here moves
//! a census on its own — [`the_pair_state_census_holds`] ties [`PAIR_STATES`]
//! to the `verdict` column, to [`PAIR_STATE_CENSUS`] and, arithmetically, to
//! [`OPEN_GAPS`], so a state that moves without all four being edited fails.
//!
//! **And it is a claim about reaching `wanted`, never about convergence.**
//! Whether the two spellings agree *with each other* is `reported_confluence_pairs`'
//! subject and is counted there. The two came apart at #1649 and stayed apart:
//! all three of today's `NeitherReaches` pairs have converged, on a string that
//! is still not their `wanted` one.
//!
//! # The reference
//!
//! All three issues ran on one synthetic 125 nt contig, reproduced here verbatim
//! so these rows measure the bases they were reported against.
//! [`reference_bases_are_what_the_reports_state`] asserts that premise before
//! anything else; without it every row below could be measuring a different
//! locus and still look green.
//!
//! # Where the wanted form's authority comes from
//!
//! Every row carries a [`Row::wanted`] string and an [`Authority`] recording
//! what backs it. `wanted` is usually the output *its issue asks for, as the
//! issue states it* — but not always, and the exception matters: #1419's six
//! rows carry the **decided chain's** answer instead, the spanning delins, since
//! `canonical-form-choice-when-both-legal` settled that shape after those issues
//! were filed. So read `wanted` as "the form this row is measured against",
//! never as "what the filer asked for". The
//! project's precedence policy is **spec-explicit > Mutalyzer > our judgement**,
//! and applying it row by row does not give one answer for the whole family:
//!
//! * **#1420 is [`Authority::SpecExplicit`].** `delins.md:16` (two or more
//!   *consecutive* changed nucleotides are one delins) settles v4;
//!   `delins.md:17` (variants separated by one or more nucleotides are described
//!   individually) settles v2 and v3. `delins.md:44-47` does not reach those two:
//!   they replace a span with a payload of the **same length**, so reference and
//!   payload correspond column for column and no alignment is ever chosen. Their
//!   one unchanged interior column is simply unchanged — not a "part of the
//!   inserted sequence" that had to be *found* to align.
//! * **#1421 is [`Authority::SpecSelfConflicting`], for the same reason #1419
//!   is.** An earlier revision labelled it `SpecExplicit`, on the grounds that
//!   `delins.md:44-47` "is conditioned on the span losing length". It is not.
//!   `:47` states its condition as *"parts of the inserted sequence `align` with
//!   the reference sequence"* and says nothing whatever about length; the
//!   `c.850_901` example happens to be a net deletion, and that incidental
//!   arithmetic was mistaken for the rule.
//!
//!   Length does bear on it, but the other way round. Every #1421 row replaces
//!   **5** reference bases with an **11**-base payload (net +6), so there is no
//!   column-for-column correspondence to read off, and the wanted split exists
//!   *only* because the two-base unchanged interior reappears in the payload —
//!   `AC` at 30-31 for n1, `TG` at 33-34 for n2, `CA` at 35-36 for n3 — with zero
//!   end-column identity on either side of any of them. Recovering that split
//!   requires choosing an alignment, which is precisely the shape `:44-47` says
//!   to describe as one delins, against `:17` which says to split. The spec
//!   therefore answers both ways, the policy falls through to Mutalyzer, and
//!   Mutalyzer picks the merged delins — so it does not pick the issue's form
//!   either. #1421's wanted forms are recorded, not asserted as targets.
//! * **#1419 is [`Authority::SpecSelfConflicting`], so its wanted form is
//!   recorded and not asserted as the target.** The issue rests on
//!   `general.md:56` (prioritisation: substitution over deletion) and
//!   `delins.md:17`. But `delins.md:44-47` addresses this exact shape — a payload
//!   whose bases align with the reference, giving "an alternative description
//!   like `c.[850_869del;874_881del;887_897del;901_902insG]`" — and says outright
//!   that **the delins format is recommended**. #1419's own `[19_30del;33T>G]` is
//!   a del-plus-sub arising from precisely that alignment coincidence, so the
//!   spec points both ways and the policy falls through to Mutalyzer.
//!
//! # What Mutalyzer answers, measured rather than assumed
//!
//! Mutalyzer cannot resolve the synthetic `TEMPLATE` accession, so its answer was
//! not transcribed for these rows. It was established by reproducing each *shape*
//! on a public accession whose bases were read back out of Mutalyzer's own
//! `ESEQUENCEMISMATCH` detail (`NG_012337.1:g.5001_5120`), then normalizing both
//! spellings through `https://mutalyzer.nl/api/normalize/`. Measured 2026-08-06:
//!
//! | shape reproduced | both spellings converge on |
//! |---|---|
//! | #1419's `[del;del]` vs spanning delins | the **spanning delins** |
//! | #1420 v2/v3's unchanged interior column | the **merged delins** |
//! | #1421's net insertion across an unchanged interior | the **merged delins** |
//!
//! Two things follow, and they point in opposite directions. Mutalyzer is
//! **confluent** on every one of these shapes, which is the property #1430 asks
//! ferro to adopt. And Mutalyzer's chosen representative is the **merged** one in
//! all three, which is what `delins.md:17` forbids for #1420 — so there the spec
//! wins outright and Mutalyzer is not consulted, exactly as the policy orders it.
//! For #1419 and #1421, where the spec conflicts with itself, Mutalyzer is the
//! tie-break the policy names, and it does **not** pick either issue's form.
//!
//! Mutalyzer's merge is also not a rule one could adopt even if the policy
//! allowed it. Two substitutions two nucleotides apart on that accession merge
//! for `T>A`/`C>A` and for `T>G`/`C>T`, and stay individual for `T>C`/`C>G` — the
//! decision follows its aligner's view of the alternate bases, not a separation
//! rule.
//!
//! # The collision this module does not resolve
//!
//! #1430 names Mutalyzer's behaviour as *"not guaranteed to be spec compliant
//! (they take some shortcuts that merge deletion-insertions)"*. That is the same
//! delins-merging the paragraph above uses to adjudicate #1419 — so where this
//! module follows Mutalyzer onto a merged delins, **the issue author's stated view
//! points the other way**. #1430's own three-stage proposal reconciles them
//! (confluence from the Mutalyzer-like derivation in stage 2, spec compliance
//! applied on top in stage 3), and `merge::coalesce`'s doc comment already cites
//! #1430 for that split. But the reconciliation is a design, not a decision that
//! has been taken, and taking it is not this module's to make. #1419's rows are
//! therefore left as characterization pins with the disagreement written down.
//!
//! # Both shuffle directions are pinned
//!
//! [`every_reported_spelling_still_normalizes_as_the_reports_recorded`] pins the
//! 3' answer the reports observed; [`Row::five_prime`] and
//! [`every_reported_spelling_normalizes_as_recorded_under_five_prime`] pin the 5'
//! answer, which nothing pinned before. Fifteen of the eighteen agree across
//! directions; the three that do not are named in [`FIVE_PRIME_MOVERS`].
//!
//! **Two of those three are `Gap` rows that the 3' pass returns verbatim and
//! the 5' pass does not** — see
//! [`every_gap_row_is_returned_exactly_as_authored`] for why that distinction
//! carries weight. The third, `1419-r3/span`, is a `Gap` row that **neither**
//! direction returns verbatim: #1649 moved it off `[19T>G;22_33del]` onto the
//! two-deletion form its `/cis` sibling prints, and the 5' pass then shifts
//! that form's leading deletion one base left. Its pair is
//! [`PairState::NeitherReaches`] in [`PAIR_STATES`], which is what carries the
//! retention exemption that follows.

use crate::common::cis_apply_oracle::{normalize, normalize_in, provider};
use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::ShuffleDirection;

/// The 125 nt contig every row runs against. Position 1 is index 0.
///
/// Synthetic — it names no real sequence. Reproduced from the reports so these
/// rows measure the same bases the reports did.
const TEMPLATE: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

/// Whether ferro's answer for one spelling is the form its issue argues for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    /// Ferro prints the form the cited recommendation makes canonical.
    Canonical,
    /// Ferro prints something else, and the row records what was wanted.
    Gap,
}

/// How many of a pair's two spellings reach the form the pair is measured
/// against.
///
/// The closed enumeration the two-state model lacked — see the module docs for
/// what it replaces and why. Every combination of the two rows' [`Verdict`]s
/// names a variant, so there is no state this declines to express and no
/// exemption list beside it.
///
/// **This is a statement about [`Row::wanted`], not about convergence.**
/// `reported_confluence_pairs` owns "do the two spellings agree"; this module
/// owns "does it reach `wanted`". The two have come apart — every
/// [`PairState::NeitherReaches`] pair below has converged — and keeping them in
/// separate modules is what makes that visible rather than confusing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PairState {
    /// Both spellings print it.
    ///
    /// This is what a **fixed** pair looks like: the canonical form has stopped
    /// being a function of how the variant was written, which is the defect
    /// #1419/#1420/#1421 report. No pair is in this state today; it is a state
    /// rather than a failure so that arriving in it is expressible.
    ///
    /// Arriving here is still not automatically progress — see [`OPEN_GAPS`] on
    /// re-derivation. The state is sayable; banking it as a fix still requires
    /// naming the clause.
    BothReach,
    /// Exactly one spelling prints it.
    ///
    /// The reported defect as filed, and the state six of the nine pairs are
    /// in. It is what the retired `assert_ne!(a.verdict, b.verdict)` read as a
    /// structural fact about every pair.
    OneReaches,
    /// Neither spelling prints it.
    ///
    /// Absorbs the former `PAIRS_NO_SPELLING_REACHES` list, whose doc recorded
    /// why it had to exist — and that reasoning is now this variant's
    /// definition. Applying the decided chain moved #1419's `wanted` onto the
    /// spanning delins, and ferro prints that for neither spelling: the cis
    /// spelling is returned verbatim as its own two-deletion form, and the span
    /// spelling is *split* into the del-plus-sub form. So both rows are `Gap`,
    /// and one of them is a `Gap` row that does **not** retain its input.
    ///
    /// Two consequences the rest of the module keys on, both inherited
    /// unchanged. There is no canonical sibling whose output a gap row's
    /// `wanted` could be compared against, so
    /// [`the_wanted_form_of_every_gap_row_is_its_siblings_output`] has nothing
    /// to say here. And retention is not the property to assert for the whole
    /// pair, so [`every_gap_row_is_returned_exactly_as_authored`] pins the
    /// measured output instead — for the `/cis` half, which does still retain
    /// its input, that is a strictly weaker pin taken deliberately rather than a
    /// dropped one.
    ///
    /// Leaving this state means a spelling started reaching the wanted form. A
    /// pair merely agreeing with itself is not that: all three of today's have
    /// converged as of #1649 (`reported_confluence_pairs`'
    /// `CONVERGING_PAIRS_THREE_PRIME` and `_FIVE_PRIME` are both **3**, not 0,
    /// and [`every_reported_pair_is_still_one_variant_by_equivalence`] pins
    /// `NormalizedMatch` for these three against the denotational rung for the
    /// other six) and all three stay here, because the string they converge on
    /// is not the spanning delins the decided chain requires. That ratchet
    /// belongs in the module that owns it, and it is there.
    NeitherReaches,
}

impl PairState {
    /// The state a pair's two pinned verdicts put it in.
    ///
    /// Total over the verdict pair, which is the property the two-state model
    /// lacked: there is no combination of verdicts this declines to name, so a
    /// pair cannot reach a shape the model has to be extended to describe.
    fn from_verdicts(a: Verdict, b: Verdict) -> Self {
        match (a, b) {
            (Verdict::Canonical, Verdict::Canonical) => Self::BothReach,
            (Verdict::Gap, Verdict::Gap) => Self::NeitherReaches,
            _ => Self::OneReaches,
        }
    }

    /// How many of a pair's two rows are `Gap` in this state.
    ///
    /// The arithmetic bridge between the per-pair census and [`OPEN_GAPS`],
    /// which counts rows. Asserted in [`the_pair_state_census_holds`] so a pair
    /// cannot change state while the row count stands still, or the other way
    /// round.
    fn gap_rows(self) -> usize {
        match self {
            Self::BothReach => 0,
            Self::OneReaches => 1,
            Self::NeitherReaches => 2,
        }
    }
}

/// What backs the [`Row::wanted`] form, under the project's precedence policy
/// (**spec-explicit > Mutalyzer > our judgement**). See the module docs for the
/// per-issue reasoning and for the Mutalyzer measurement it rests on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Authority {
    /// The spec answers this shape outright at the line the issue cites, and
    /// nothing else in the spec answers it the other way. Mutalyzer is not
    /// consulted, and the wanted form is the target a fix should reach.
    SpecExplicit,
    /// The spec answers this shape *both* ways, so the policy falls through to
    /// Mutalyzer.
    ///
    /// **The two families under this authority differ, so do not read one rule
    /// off it.** For #1421 the wanted form is recorded as what the issue asked
    /// for and is not asserted as the target — Mutalyzer converges that shape on
    /// the merged delins, so it does not pick the issue's form. For #1419 the
    /// wanted form is the decided chain's answer rather than the issue's, and
    /// the tie-break is never reached.
    SpecSelfConflicting,
}

/// One spelling of one reported variant.
struct Row {
    /// `<issue>-<row>`, matching the label used in `reported_confluence_pairs`
    /// so the two modules' rows can be lined up by eye.
    label: &'static str,
    /// The spelling as authored in the issue.
    input: &'static str,
    /// What ferro prints for it today under the shipped 3' direction. Measured,
    /// never transcribed.
    output: &'static str,
    /// What ferro prints for it under `ShuffleDirection::FivePrime`. Measured.
    /// Equal to [`Row::output`] for fifteen of the eighteen rows; the three
    /// that differ are listed in [`FIVE_PRIME_MOVERS`].
    five_prime: &'static str,
    verdict: Verdict,
    /// The output **the issue asks for**, in the issue's own terms — the string
    /// its "canonical"/"expected" column names for this spelling. For a
    /// [`Verdict::Canonical`] row this equals [`Row::output`].
    ///
    /// For a [`Verdict::Gap`] row it is the sibling spelling's output **only in
    /// a [`PairState::OneReaches`] pair**, where that is the family's defining
    /// shape and is asserted rather than assumed by
    /// [`the_wanted_form_of_every_gap_row_is_its_siblings_output`]. In a
    /// [`PairState::NeitherReaches`] pair both rows are `Gap`, there is no
    /// canonical sibling, and `wanted` is the decided chain's spanning delins —
    /// which is why that test skips those pairs rather than asserting an
    /// identity that does not hold.
    wanted: &'static str,
    /// What backs [`Row::wanted`]. Both spellings of a pair share one value:
    /// the authority is a property of the shape, not of how it was written.
    authority: Authority,
    /// Why that is or is not the canonical form, citing the line of the pinned
    /// spec checkout the issue cites.
    argument: &'static str,
}

/// The eighteen spellings — two per reported pair.
///
/// Both spellings of a pair are listed rather than only the losing one, because
/// the defect is not "this string is wrong": both strings are well-formed and
/// both denote the same sequence. The defect is that the canonical form is
/// reachable from one spelling and not from the other, and you cannot state that
/// without both halves. See
/// [`each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names`].
const REPORTED_ROWS: &[Row] = &[
    // -- #1419 -- a net deletion spelled as two deletions, versus as one span.
    //
    // Neither spelling is a delins across an unchanged base, so no local rule
    // is violated by either; the issue's own framing is that this is pure
    // partition non-uniqueness. The canonical argument is prioritisation.
    Row {
        label: "1419-r1/cis",
        input: "TEMPLATE:g.[19_23del;27_33del]",
        output: "TEMPLATE:g.[19_23del;27_33del]",
        five_prime: "TEMPLATE:g.[19_23del;27_33del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.19_33delinsCGG",
        authority: Authority::SpecSelfConflicting,
        argument: "CORRECTED. This row wanted `[19_30del;33T>G]`, the form \
                   #1419 asks for, on the ground that general.md:56 \
                   prioritises (1) substitution over (2) deletion so the \
                   change at 33 should be exposed rather than absorbed. \
                   The decided chain, applied. `adjudication-precedence-order` records that general.md:56 \"ranks single-variant TYPE LABELS FOR ONE SPAN — it never ranks a multi-member allele against a spanning description\", and that an earlier attempt to use it that way was refuted on exactly that ground, so :56 cannot settle a merge-versus-split question at all. :56 was the whole of the argument this row used to carry. `canonical-form-choice-when-both-legal` supplies what to do instead: derive from the resulting sequence, then apply every explicit spec tie-break. :56 does not apply; delins.md:46-47 does, this being precisely the alignment-coincidence shape it constructs and then declines — `The \"delins\" format is recommended`. So the wanted form is the SPANNING delins. Ferro prints \
                   neither form, so the verdict stays Gap.",
    },
    Row {
        label: "1419-r1/span",
        input: "TEMPLATE:g.19_33delinsCGG",
        // Measured, and MOVED by #1649's two-deletion alignment:
        // `[19_30del;33T>G]` -> `[19_23del;27_33del]`, in both directions. The
        // splitter can now express *deletion, retained reference, deletion*, so
        // this span no longer has to be flattened onto the del-plus-sub form.
        // What that changes is WHICH split ferro prints — not `wanted`, which
        // is still the spanning delins on delins.md:46-47, and not `verdict`,
        // which is still Gap because ferro still does not print it. The pair
        // now CONVERGES: this is byte-identical to `1419-r1/cis`'s output, so
        // `CONVERGING_PAIRS_THREE_PRIME`/`_FIVE_PRIME` are 3 rather than 0 and
        // the equivalence level is NormalizedMatch. A migration for anyone
        // storing the old string — see the PR's Representation-Change trailer.
        output: "TEMPLATE:g.[19_23del;27_33del]",
        five_prime: "TEMPLATE:g.[19_23del;27_33del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.19_33delinsCGG",
        authority: Authority::SpecSelfConflicting,
        argument: "CORRECTED, and the verdict flipped Canonical -> Gap. This \
                   row was marked Canonical for printing `[19_30del;33T>G]`, \
                   argued from general.md:56 plus a general.md:41 reading that \
                   placed the RETAINED bases as 3' as possible — which is the \
                   3' rule inverted, since :41 assigns the most 3' position to \
                   have been CHANGED. Both grounds are withdrawn. \
                   The decided chain, applied. `adjudication-precedence-order` records that general.md:56 \"ranks single-variant TYPE LABELS FOR ONE SPAN — it never ranks a multi-member allele against a spanning description\", and that an earlier attempt to use it that way was refuted on exactly that ground, so :56 cannot settle a merge-versus-split question at all. :56 was the whole of the argument this row used to carry. `canonical-form-choice-when-both-legal` supplies what to do instead: derive from the resulting sequence, then apply every explicit spec tie-break. :56 does not apply; delins.md:46-47 does, this being precisely the alignment-coincidence shape it constructs and then declines — `The \"delins\" format is recommended`. So the wanted form is the SPANNING delins. That is this row's own \
                   input, so the wanted form is what it arrived as — and ferro \
                   splits it to `[19_23del;27_33del]` instead, which is why this \
                   is a Gap row that does NOT retain its input. Ferro's own \
                   coalesce pass encodes the same recommendation on the spec's \
                   worked example (see \
                   `merge::the_coalesce_pass_reaches_the_spec_worked_example`) \
                   and does not reach this block.",
    },
    Row {
        label: "1419-r2/cis",
        input: "TEMPLATE:g.[19_22del;26_36del]",
        output: "TEMPLATE:g.[19_22del;26_36del]",
        five_prime: "TEMPLATE:g.[19_22del;26_36del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.19_36delinsGCG",
        authority: Authority::SpecSelfConflicting,
        argument: "CORRECTED, one locus over from `1419-r1/cis` and for the \
                   same reason: the general.md:56 ground it wanted \
                   `[19_33del;36A>G]` on is refuted for merge-versus-split, \
                   and the decided chain lands on the spanning delins. See \
                   `1419-r1/cis` for the argument in full.",
    },
    Row {
        label: "1419-r2/span",
        input: "TEMPLATE:g.19_36delinsGCG",
        // Measured, and MOVED by #1649 — see `1419-r1/span`, same mechanism
        // and same consequences: `[19_33del;36A>G]` -> `[19_22del;26_36del]`
        // in both directions, which is `1419-r2/cis`'s output, so this pair
        // converges too. `wanted` and `verdict` are untouched.
        output: "TEMPLATE:g.[19_22del;26_36del]",
        five_prime: "TEMPLATE:g.[19_22del;26_36del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.19_36delinsGCG",
        authority: Authority::SpecSelfConflicting,
        argument: "CORRECTED, verdict flipped Canonical -> Gap, for the same \
                   reason as `1419-r1/span` one locus over. The wanted form is \
                   this row's own input; ferro splits it instead.",
    },
    Row {
        label: "1419-r3/cis",
        input: "TEMPLATE:g.[19_24del;28_33del]",
        output: "TEMPLATE:g.[19_24del;28_33del]",
        // The 5' pass shifts the leading deletion one base left: 19-24 and
        // 18-23 delete the same content out of the run `C T G A T G C` at
        // 18-24, so the denoted sequence is unchanged (asserted by
        // `reported_confluence_pairs::no_reported_pair_normalizes_to_a_\
        // different_sequence`, which sweeps both directions).
        five_prime: "TEMPLATE:g.[18_23del;28_33del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.19_33delinsGGA",
        authority: Authority::SpecSelfConflicting,
        // TRAP, see OPEN_GAPS. The wanted form here is ALSO what a
        // re-partitioning of the block from the sequence produces: two
        // deletions three unchanged nucleotides apart (25, 26, 27) coming back
        // as a substitution plus a deletion. So this row can close for the
        // wrong reason, and closing it for the wrong reason looks identical to
        // closing it for the right one. Do not flip the verdict without naming
        // the clause that licensed the merge.
        argument: "CORRECTED. This row wanted `[19T>G;22_33del]` on the same \
                   general.md:56 ground, with the exposed substitution landing \
                   at the 5' end rather than the 3' end — which is why #1419 \
                   argues no per-end rule fixes the family. That observation \
                   survives the correction: it is an argument against a \
                   per-end rule, not for :56. The ground is refuted for \
                   merge-versus-split all the same, and the decided chain \
                   lands on the spanning delins. See `1419-r1/cis`.",
    },
    Row {
        label: "1419-r3/span",
        input: "TEMPLATE:g.19_33delinsGGA",
        // Measured, and MOVED by #1649 — see `1419-r1/span`. Two differences
        // from its two siblings, both measured. It lands on
        // `[19_24del;28_33del]` rather than keeping the `[19T>G;22_33del]`
        // substitution-plus-deletion form, which is `1419-r3/cis`'s output, so
        // this pair converges as well. And it is now a FIVE_PRIME_MOVER: the 5'
        // pass shifts the leading deletion one base left exactly as it does for
        // `1419-r3/cis`, so the two directions no longer agree here. The
        // argument below still says ferro splits this row rather than printing
        // its `wanted` spanning delins, which remains true of both directions.
        output: "TEMPLATE:g.[19_24del;28_33del]",
        five_prime: "TEMPLATE:g.[18_23del;28_33del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.19_33delinsGGA",
        authority: Authority::SpecSelfConflicting,
        argument: "CORRECTED, verdict flipped Canonical -> Gap, for the same \
                   reason as `1419-r1/span`. The wanted form is this row's own \
                   input; ferro splits it to `[19_24del;28_33del]` instead, and \
                   to `[18_23del;28_33del]` under 5' — it is a FIVE_PRIME_MOVER \
                   too now, alongside its `1419-r3/cis` sibling.",
    },
    // -- #1420 -- one row per member type. v2/v3 must reduce, v4 must coalesce:
    // the corrections point in opposite directions, which is the issue's
    // argument that no per-member-type rule reconciles them.
    Row {
        label: "1420-v2/cis",
        input: "TEMPLATE:g.[37dup;41del]",
        // CLOSED by #1616: the reported spelling now re-derives to its wanted
        // form, and to the same string its `/span` sibling prints, so the pair
        // converges. The `dup` is gone, which is also why the 5' anchor note
        // this row used to carry is gone: there is no duplication left for the
        // 5' pass to shift, and both directions print one string.
        output: "TEMPLATE:g.[38T>A;40_41delinsTG]",
        five_prime: "TEMPLATE:g.[38T>A;40_41delinsTG]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[38T>A;40_41delinsTG]",
        authority: Authority::SpecExplicit,
        // THE TRAP THIS ROW CARRIED IS RESOLVED, not silently dropped — operator
        // ruling, 2026-08-13. The note here read: "this row's members are three
        // unchanged nucleotides apart (38, 39, 40), and the wanted form is what
        // a re-partitioning across them produces — so the row can converge on
        // #1420's own ask by violating `general.md:34`."
        //
        // That reads the separation off THE INPUT'S SPELLING (members at 37 and
        // 41). `rulings[separation-is-a-property-of-the-spelling-not-of-the-\
        // variant]` (decided) holds it is read off the partition re-derived from
        // the resulting sequence. Read that way the output is TWO MEMBERS AT
        // SEPARATION ONE, described individually across the unchanged 39 —
        // exactly what `:34` asks for. **Nothing merged**, so the trap's
        // premise does not hold here. The clause carrying the move is
        // `rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`, whose
        // scope paragraph forbids licensing merges — respected, for the same
        // reason. The companion guard in `reported_confluence_pairs` was deleted
        // in this change; its block comment carries the full reasoning.
        argument: "Wanted `[38T>A;40_41delinsTG]`. Reference 38-41 is `TTGC` \
                   and the result is `ATTG`, so 38 and 40-41 change and 39 \
                   does not (asserted in `reported_spans_change_the_columns_\
                   the_reports_state`). general.md:56 prioritises \
                   (1) substitution over (4) duplication, so the change at 38 \
                   must not be spelled as a `dup`; delins.md:17 keeps the two \
                   runs individual across the unchanged 39. Mutalyzer merges \
                   this shape instead, but the spec is explicit, so under \
                   spec-explicit > Mutalyzer it does not get a vote.",
    },
    Row {
        label: "1420-v2/span",
        input: "TEMPLATE:g.38_41delinsATTG",
        output: "TEMPLATE:g.[38T>A;40_41delinsTG]",
        five_prime: "TEMPLATE:g.[38T>A;40_41delinsTG]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[38T>A;40_41delinsTG]",
        authority: Authority::SpecExplicit,
        argument: "The form #1420 names canonical for v2: substitution exposed \
                   at 38, and the span separated at the unchanged 39 per \
                   delins.md:17.",
    },
    Row {
        label: "1420-v3/cis",
        input: "TEMPLATE:g.[36_37insC;40del]",
        output: "TEMPLATE:g.[36_37insC;40del]",
        five_prime: "TEMPLATE:g.[36_37insC;40del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.[37_38delinsCA;40G>T]",
        authority: Authority::SpecExplicit,
        argument: "Wanted `[37_38delinsCA;40G>T]`. general.md:56 prioritises \
                   (1) substitution over (5) insertion, so the substitution at \
                   40 must be exposed rather than left inside an `ins`+`del` \
                   pair.",
    },
    Row {
        label: "1420-v3/span",
        input: "TEMPLATE:g.37_40delinsCATT",
        output: "TEMPLATE:g.[37_38delinsCA;40G>T]",
        five_prime: "TEMPLATE:g.[37_38delinsCA;40G>T]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[37_38delinsCA;40G>T]",
        authority: Authority::SpecExplicit,
        argument: "The form #1420 names canonical for v3: reference 37-40 is \
                   `ATTG` against `CATT`, so 37-38 change consecutively \
                   (delins.md:16), 39 does not, and 40 is a substitution \
                   separated from them (delins.md:17).",
    },
    Row {
        label: "1420-v4/cis",
        input: "TEMPLATE:g.[21delinsGC;24del]",
        output: "TEMPLATE:g.[21delinsGC;24del]",
        five_prime: "TEMPLATE:g.[21delinsGC;24del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.21_24delinsGCTG",
        authority: Authority::SpecExplicit,
        argument: "Wanted `21_24delinsGCTG` — the opposite direction to v2/v3. \
                   Reference 21-24 is `ATGC` against `GCTG`, so all four \
                   positions change, and delins.md:16 says changes involving \
                   two or more CONSECUTIVE nucleotides are one delins. The \
                   split spelling asserts 22-23 are unchanged; they are not.",
    },
    Row {
        label: "1420-v4/span",
        input: "TEMPLATE:g.21_24delinsGCTG",
        output: "TEMPLATE:g.21_24delinsGCTG",
        five_prime: "TEMPLATE:g.21_24delinsGCTG",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.21_24delinsGCTG",
        authority: Authority::SpecExplicit,
        argument: "The form #1420 names canonical for v4: one delins over four \
                   consecutive changed positions, per delins.md:16.",
    },
    // -- #1421 -- net insertions. Here the SPLIT spelling is the canonical one,
    // the reverse of #1419/#1420, so a fix that simply preferred whichever form
    // has fewer members would get this family backwards.
    Row {
                // PINNED TO THE MEASURED STRINGS by #1616 — NOT a ruling on the form.
        // `rulings[canonical-form-choice-when-both-legal]` (decided) licenses
        // deriving from the resulting sequence and emitting what falls out, so a
        // form neither spelling authored is not wrong on its face. Whether THIS
        // form is right is a `:47`/`:34` question that is deliberately NOT
        // answered here; what is fixed is that the movement is no longer
        // unpinned. Disclosed in the PR's `Representation-Change:` trailer.
        //
        // NOTE FOR THE #1421 RELABEL: this touches `output`/`five_prime` only.
        // The `authority` field is untouched on every row.
        label: "1421-n1/split",
        input: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        output: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        five_prime: "TEMPLATE:g.[29delinsAACACAT;32_33delinsTG]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        authority: Authority::SpecSelfConflicting,
        argument: "The form #1421 names canonical: the substitution at 29 and \
                   the delins at 32-33 are separated by two unchanged \
                   nucleotides, and delins.md:17 says variants separated by one \
                   or more nucleotides are described individually and not as a \
                   delins. Recorded, not targeted: delins.md:44-47 answers the \
                   other way on this shape (see the sibling span row).",
    },
    Row {
                // PINNED TO THE MEASURED STRINGS by #1616 — NOT a ruling on the form.
        // `rulings[canonical-form-choice-when-both-legal]` (decided) licenses
        // deriving from the resulting sequence and emitting what falls out, so a
        // form neither spelling authored is not wrong on its face. Whether THIS
        // form is right is a `:47`/`:34` question that is deliberately NOT
        // answered here; what is fixed is that the movement is no longer
        // unpinned. Disclosed in the PR's `Representation-Change:` trailer.
        //
        // NOTE FOR THE #1421 RELABEL: this touches `output`/`five_prime` only.
        // The `authority` field is untouched on every row.
        label: "1421-n1/span",
        input: "TEMPLATE:g.29_33delinsAACACATACTG",
        output: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        five_prime: "TEMPLATE:g.[29delinsAACACAT;32_33delinsTG]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        authority: Authority::SpecSelfConflicting,
        argument: "Wanted `[29C>A;32_33delinsACATACTG]`. This span merges \
                   across the unchanged 30-31, which delins.md:17 forbids. The \
                   delins.md:18 exception (two variants separated by ONE \
                   nucleotide together affecting one amino acid) does not \
                   reach it: the separation is two nucleotides, and a `g.` \
                   description has no amino acid. See \
                   `the_1421_spans_separate_by_two_nucleotides_not_one`. \
                   delins.md:44-47 DOES reach it, though, which is why this row \
                   is SpecSelfConflicting: 5 reference bases (CACGT) become an \
                   11-base payload, so there is no column-for-column reading and \
                   the split is recoverable only because the unchanged 30-31 \
                   `AC` reappears in the payload — with neither end column \
                   matching. That is the alignment `:47` says to write as one \
                   delins. The spec answers both ways, so the policy falls \
                   through to Mutalyzer, which merges and so does not pick this \
                   wanted form either.",
    },
    Row {
                // PINNED TO THE MEASURED STRINGS by #1616 — NOT a ruling on the form.
        // `rulings[canonical-form-choice-when-both-legal]` (decided) licenses
        // deriving from the resulting sequence and emitting what falls out, so a
        // form neither spelling authored is not wrong on its face. Whether THIS
        // form is right is a `:47`/`:34` question that is deliberately NOT
        // answered here; what is fixed is that the movement is no longer
        // unpinned. Disclosed in the PR's `Representation-Change:` trailer.
        //
        // NOTE FOR THE #1421 RELABEL: this touches `output`/`five_prime` only.
        // The `authority` field is untouched on every row.
        label: "1421-n2/split",
        input: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        output: "TEMPLATE:g.[32G>T;35C>G;38_39insCGACAT]",
        five_prime: "TEMPLATE:g.[32G>T;35C>G;36_37insATCGAC]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        authority: Authority::SpecSelfConflicting,
        argument: "The form #1421 names canonical for row 2; unchanged \
                   interior 33-34. Recorded, not targeted — delins.md:44-47 \
                   answers the other way (see the sibling span row).",
    },
    Row {
                // PINNED TO THE MEASURED STRINGS by #1616 — NOT a ruling on the form.
        // `rulings[canonical-form-choice-when-both-legal]` (decided) licenses
        // deriving from the resulting sequence and emitting what falls out, so a
        // form neither spelling authored is not wrong on its face. Whether THIS
        // form is right is a `:47`/`:34` question that is deliberately NOT
        // answered here; what is fixed is that the movement is no longer
        // unpinned. Disclosed in the PR's `Representation-Change:` trailer.
        //
        // NOTE FOR THE #1421 RELABEL: this touches `output`/`five_prime` only.
        // The `authority` field is untouched on every row.
        label: "1421-n2/span",
        input: "TEMPLATE:g.32_36delinsTTGGAATCGAC",
        output: "TEMPLATE:g.[32G>T;35C>G;38_39insCGACAT]",
        five_prime: "TEMPLATE:g.[32G>T;35C>G;36_37insATCGAC]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        authority: Authority::SpecSelfConflicting,
        argument: "Wanted `[32G>T;35_36delinsGAATCGAC]` per delins.md:17, but \
                   delins.md:44-47 answers the other way: ref 32-36 (GTGCA) \
                   becomes an 11-base payload, and the split survives only \
                   because the unchanged 33-34 `TG` reappears in it, neither end \
                   column matching. SpecSelfConflicting for that reason.",
    },
    Row {
                // PINNED TO THE MEASURED STRINGS by #1616 — NOT a ruling on the form.
        // `rulings[canonical-form-choice-when-both-legal]` (decided) licenses
        // deriving from the resulting sequence and emitting what falls out, so a
        // form neither spelling authored is not wrong on its face. Whether THIS
        // form is right is a `:47`/`:34` question that is deliberately NOT
        // answered here; what is fixed is that the movement is no longer
        // unpinned. Disclosed in the PR's `Representation-Change:` trailer.
        //
        // NOTE FOR THE #1421 RELABEL: this touches `output`/`five_prime` only.
        // The `authority` field is untouched on every row.
        label: "1421-n3/split",
        input: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        output: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        five_prime: "TEMPLATE:g.[34G>T;35_36insACCTTT;37_38delinsCG]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        authority: Authority::SpecSelfConflicting,
        argument: "The form #1421 names canonical for row 3; unchanged \
                   interior 35-36. Recorded, not targeted — delins.md:44-47 \
                   answers the other way (see the sibling span row).",
    },
    Row {
                // PINNED TO THE MEASURED STRINGS by #1616 — NOT a ruling on the form.
        // `rulings[canonical-form-choice-when-both-legal]` (decided) licenses
        // deriving from the resulting sequence and emitting what falls out, so a
        // form neither spelling authored is not wrong on its face. Whether THIS
        // form is right is a `:47`/`:34` question that is deliberately NOT
        // answered here; what is fixed is that the movement is no longer
        // unpinned. Disclosed in the PR's `Representation-Change:` trailer.
        //
        // NOTE FOR THE #1421 RELABEL: this touches `output`/`five_prime` only.
        // The `authority` field is untouched on every row.
        label: "1421-n3/span",
        input: "TEMPLATE:g.34_38delinsTCACCTTTACG",
        output: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        five_prime: "TEMPLATE:g.[34G>T;35_36insACCTTT;37_38delinsCG]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        authority: Authority::SpecSelfConflicting,
        argument: "Wanted `[34G>T;37_38delinsCCTTTACG]` per delins.md:17, but \
                   delins.md:44-47 answers the other way: ref 34-38 (GCAAT) \
                   becomes an 11-base payload, and the split survives only \
                   because the unchanged 35-36 `CA` reappears in it, neither end \
                   column matching. SpecSelfConflicting for that reason.",
    },
];

/// The three rows whose 5' answer differs from their 3' answer.
///
/// Named rather than left implicit because all three are `Gap` rows, and
/// [`every_gap_row_is_returned_exactly_as_authored`] rests on a gap row being
/// returned *untouched* — which is what makes it a second canonical form rather
/// than an under-applied pass. That argument is a 3'-direction argument, and
/// these are where it stops holding.
///
/// **`1419-r3/span` was added by #1649**, which moved that row off
/// `[19T>G;22_33del]` onto the same two-deletion form its `/cis` sibling
/// prints — and the 5' pass shifts that form's leading deletion one base left,
/// which the substitution-plus-deletion form gave it no way to do. So the row
/// joins its own `/cis` sibling here rather than arriving for a new reason.
///
/// **`1420-v2/cis` was removed by #1616.** It was here because the 5' pass moved
/// its duplication's anchor one base left; the row no longer prints a `dup` at
/// all, so the two directions now print one string and there is nothing to move.
/// **Six #1421 rows joined by #1616, and they are the measured set rather than a
/// chosen one.** Under 5' all three #1421 pairs land on a form that is neither
/// spelling's `wanted` — `1421-n2` does so in BOTH directions. This constant is
/// asserted as a set by
/// [`only_the_named_rows_answer_differently_in_the_two_directions`] precisely so
/// the movers cannot drift silently behind the per-row pins, so it is updated
/// with them. **No judgement is implied about the forms themselves** — see the
/// `PINNED TO THE MEASURED STRINGS` notes on the rows.
const FIVE_PRIME_MOVERS: &[&str] = &[
    "1419-r3/cis",
    "1419-r3/span",
    "1421-n1/split",
    "1421-n1/span",
    "1421-n2/split",
    "1421-n2/span",
    "1421-n3/split",
    "1421-n3/span",
];

/// How many of [`REPORTED_ROWS`] print something their issue argues against.
///
/// Eleven of eighteen. It was nine — exactly one per pair — then #1419's three
/// pairs began contributing **two** each, because neither of their spellings
/// prints the wanted form any more — [`PairState::NeitherReaches`], which is
/// where that combination is now expressed — taking it to twelve. #1616 then
/// closed `1420-v2`, which is the `12 -> 11` section below. See
/// [`each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names`].
///
/// **This may only go down.** Up means a change has moved a row *away* from the
/// cited recommendation — with one exception, taken here and named so it cannot
/// be taken again silently: the rise from nine to twelve is a **correction of
/// the target**, not a regression of the output. #1419's six `wanted` fields
/// rested entirely on `general.md:56`, which `adjudication-precedence-order`
/// records as refuted for merge-versus-split, and the decided chain lands them
/// on the spanning `delins` instead. Ferro prints neither the old target nor
/// the new one, so three rows that were `Canonical` against the withdrawn
/// target are `Gap` against the standing one.
///
/// **When this rose to twelve, no output moved** — the rise was entirely in
/// `wanted`, which is an adjudication, and in `verdict`, which is derived from
/// the two.
///
/// **That is no longer true of the file as a whole, and the paragraph above is
/// scoped to the rise rather than describing the current state.** #1649 moved
/// six fields: `output` and `five_prime` on all three `1419-r*/span` rows. The
/// count did not move with them and that is the point — `wanted`, `verdict`
/// and `authority` are untouched on every row, so all three remain `Gap`
/// against the same standing target for the same reason, and `OPEN_GAPS` stays
/// twelve. A moved output with a still-unreached target is a representation
/// change, not a closure. Both directions are measured by
/// [`every_reported_spelling_still_normalizes_as_the_reports_recorded`] and
/// [`every_reported_spelling_normalizes_as_recorded_under_five_prime`].
///
/// # Going down is not automatically progress either
///
/// A gap row closes when its output becomes its `wanted` string. There are two
/// ways for that to happen and they are not equivalent:
///
///  * a **licensed move** — a pass that implements the clause the issue cites,
///    which is the fix; or
///  * a **re-derivation** — the block re-partitioned from the sequence across
///    bases the input left unchanged, landing on a form *neither* spelling
///    asserted and which happens to equal the issue's ask.
///
/// The second has been observed, on an experimental partitioner arm. It was
/// recorded against two rows: `1419-r3/cis` (`g.[19_24del;28_33del]` ->
/// `g.[19T>G;22_33del]`) and `1420-v2/cis` (`g.[37dup;41del]` ->
/// `g.[38T>A;40_41delinsTG]`).
///
/// # CORRECTED by #1616 — this paragraph called those `general.md:34` violations
///
/// It described both as "three unchanged nucleotides apart" and therefore as
/// `general.md:34` violations. **That characterisation is ruled wrong** and it
/// fails in exactly the way the guard deleted with it did: 38, 39 and 40 are
/// unchanged only in the INPUT'S SPELLING (members written at 37 and 41), and
/// `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]`
/// (decided) reads separation off the partition re-derived from the resulting
/// sequence. Read that way `g.[38T>A;40_41delinsTG]` is two members at
/// separation one, described individually — what `:34` asks for. **Nothing
/// merged.**
///
/// This doc also named
/// `reported_confluence_pairs::the_1420_v2_pair_does_not_converge_by_re_derivation`
/// as the pin holding the hazard. **That test no longer exists** — it was
/// deleted by the same operator ruling (2026-08-13); its block comment in that
/// file carries the argument.
///
/// The DICHOTOMY above still stands and is still the thing to apply. What is
/// withdrawn is the claim that these two rows instantiate its second horn. Note
/// the hazard is a row reaching **its wanted string** by re-derivation, which is
/// how a defect gets banked as a fix; a row landing on a form that is *not*
/// `wanted` banks nothing and is not this hazard.
///
/// So lowering this constant requires naming, in the PR, **which clause carried
/// the move** — not merely that the row now prints its wanted form. A row that
/// converged with no clause behind it is a re-derivation and the count must not
/// move.
///
/// **[`PAIR_STATES`] does not weaken any of this, and was written not to.** It
/// makes `BothReach` *expressible*, which is a different thing from making it
/// *automatic*: the state is pinned per pair, so a pair arriving there is an
/// edit to that table, this constant, and the affected rows' verdicts, argued
/// for together. [`the_pair_state_census_holds`] asserts the arithmetic between
/// the two censuses — 0, 1 or 2 gap rows per pair by state, summing to this
/// number — so neither can be moved to accommodate the other in silence.
///
/// # 12 -> 11 (#1616): `1420-v2/cis`, and which clause carried it
///
/// **`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`** (decided)
/// deleted the input-relative weight bound, so the block is re-derived from the
/// resulting sequence rather than handed back. The stored spelling that moved is
/// `TEMPLATE:g.[37dup;41del]` -> `TEMPLATE:g.[38T>A;40_41delinsTG]`, in **both**
/// directions.
///
/// **And the re-derivation half of the dichotomy above does not apply**, which
/// is the part that has to be argued rather than asserted. It was recorded here
/// as the live hazard: the wanted form is "what a re-partitioning across bases
/// the input left unchanged produces". That reads separation off the INPUT'S
/// SPELLING; `rulings[separation-is-a-property-of-the-spelling-not-of-the-variant]`
/// (decided) reads it off the re-derived partition, where the output is two
/// members at separation one described individually — `general.md:34` satisfied,
/// nothing merged. The destination is independently spec-correct on this row's
/// own committed `argument`: `general.md:56` ranks substitution above
/// duplication so 38 may not be a `dup`, and `delins.md:17` keeps the runs
/// individual across the unchanged 39.
///
/// Operator ruling, 2026-08-13. The companion guard
/// `reported_confluence_pairs::the_1420_v2_pair_does_not_converge_by_re_derivation`
/// was deleted in the same change — it pinned this string as forbidden on the
/// input-spelling reading, and two guards in one PR were giving opposite
/// instructions for one event.
const OPEN_GAPS: usize = 11;

/// Every reported pair's [`PairState`], pinned per pair and in the table's own
/// order.
///
/// **Every pair is listed, including the five in the unremarkable state.** That
/// is the difference between this and the `PAIRS_NO_SPELLING_REACHES` list it
/// replaces, which named only the exceptional pairs and so could say nothing
/// about the others: a pair leaving `OneReaches` was a failure of an assertion
/// that had assumed it could not happen, rather than a census moving. Here every
/// transition in the space is the same kind of event — one row of this table
/// changing — and there is no state left that the model has to be extended to
/// describe.
///
/// **Pinned, so nothing ratchets on its own.** Moving an entry is a deliberate
/// edit that must be accompanied by the affected rows' `verdict` fields,
/// [`PAIR_STATE_CENSUS`] and [`OPEN_GAPS`] — [`the_pair_state_census_holds`]
/// fails on any of the three lagging — and, per [`OPEN_GAPS`], by the **named
/// clause** that carried the move. A pair that converged with no clause behind
/// it is a re-derivation and the census must not move.
///
/// Today: three `NeitherReaches` (#1419's, for the reason on that variant),
/// five `OneReaches`, and one `BothReach` — `1420-v2`, whose two spellings both
/// print the wanted form as of #1616. That entry is not a ratchet: it is derived
/// from the pair's two `verdict` fields, which moved in the same change under
/// the clause [`OPEN_GAPS`]'s `12 -> 11 (#1616)` section names
/// (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`).
const PAIR_STATES: &[(&str, PairState)] = &[
    ("1419-r1", PairState::NeitherReaches),
    ("1419-r2", PairState::NeitherReaches),
    ("1419-r3", PairState::NeitherReaches),
    ("1420-v2", PairState::BothReach),
    ("1420-v3", PairState::OneReaches),
    ("1420-v4", PairState::OneReaches),
    ("1421-n1", PairState::OneReaches),
    ("1421-n2", PairState::OneReaches),
    ("1421-n3", PairState::OneReaches),
];

/// How many pairs sit in each state: `(BothReach, OneReaches, NeitherReaches)`.
///
/// The family-wide total, asserted alongside the per-pair entries rather than
/// instead of them. Per-pair alone, two compensating edits — one pair moved into
/// a state and another out of it — leave every individual assertion satisfied
/// while the family's shape has changed; the tuple is what makes that a
/// three-number diff a reviewer has to argue for.
///
/// It is one tuple rather than three constants so the three numbers cannot be
/// edited apart. What actually forces a **new** [`PairState`] variant to be
/// accounted for is the exhaustive `match` over it in
/// [`the_pair_state_census_holds`] and in [`PairState::gap_rows`]: neither
/// compiles until the new variant is given an arm, and the census arm has
/// nowhere to put its count until this tuple grows too.
const PAIR_STATE_CENSUS: (usize, usize, usize) = (1, 5, 3);

/// The state [`PAIR_STATES`] pins for one pair.
///
/// Panics rather than defaulting when a pair is unlisted: a reported pair with
/// no census entry is a table that has fallen out of step, and a silent default
/// would pick one of the three answers for it.
fn pinned_state(pair: &str) -> PairState {
    PAIR_STATES
        .iter()
        .find(|(name, _)| *name == pair)
        .map(|(_, state)| *state)
        .unwrap_or_else(|| panic!("no PAIR_STATES entry for `{pair}`"))
}

/// The `<issue>-<row>` half of an `<issue>-<row>/<spelling>` label.
fn pair_of(label: &'static str) -> &'static str {
    label
        .split('/')
        .next()
        .expect("labels are `<issue>-<row>/<spelling>`")
}

/// [`REPORTED_ROWS`] read two at a time, checked to really be nine pairs.
///
/// `chunks(2)` on its own trusts the table's row order, and that trust is not
/// free. Reorder [`REPORTED_ROWS`] so a spelling sits beside a *different*
/// pair's spelling and
/// [`each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names`]
/// still passes — measured, not assumed: every assertion in it compares a row
/// against its own pinned output or its own `wanted`, so a mispaired neighbour
/// changes nothing it looks at. Only
/// [`every_reported_pair_is_still_one_variant_by_equivalence`] notices, and it
/// reports the mispairing as `NotEquivalent` — that is, as a regression in the
/// equivalence fallback, which is the wrong diagnosis and sends a reader to the
/// wrong file.
///
/// [`the_pair_state_census_holds`] adds a second, cheaper tell for the same
/// error: it zips [`PAIR_STATES`] against these pairs and checks the labels line
/// up, so a reorder that survives the split above is reported as the two tables
/// disagreeing rather than as an equivalence regression.
///
/// Splitting the label costs nothing and makes the table error itself the
/// failure.
fn reported_pairs() -> impl Iterator<Item = (&'static Row, &'static Row)> {
    REPORTED_ROWS.chunks(2).map(|chunk| {
        let [a, b] = chunk else {
            panic!("REPORTED_ROWS must hold both spellings of every pair")
        };
        assert_eq!(
            pair_of(a.label),
            pair_of(b.label),
            "{} and {} are not two spellings of one reported pair; \
             REPORTED_ROWS has fallen out of pair order",
            a.label,
            b.label
        );
        (a, b)
    })
}

/// Normalize against [`TEMPLATE`] in the shipped default 3' direction.
///
/// 3' only, and on purpose: this module pins the strings the reports observed,
/// and the reports ran through the default. `reported_confluence_pairs` is the
/// module that sweeps both directions, because *convergence* has to hold in
/// both for a fix to count, whereas the canonical string itself is
/// direction-dependent by construction (general.md:41).
fn normalized(input: &str) -> String {
    normalize(TEMPLATE, input)
}

/// The premise: the bases the reports' arguments rest on are the bases here.
///
/// Asserted first and separately. Every canonical argument in
/// [`REPORTED_ROWS`] is a claim about which reference positions change, so if
/// the contig above were transcribed wrong, the rows would still pass — they
/// only compare ferro against itself — while every `argument` field silently
/// became fiction.
#[test]
fn reference_bases_are_what_the_reports_state() {
    assert_eq!(TEMPLATE.len(), 125, "contig length");
    for (first, last, expected, source) in [
        (21usize, 24usize, "ATGC", "#1420 v4: reference 21-24"),
        (29, 33, "CACGT", "#1421: positions 29-33"),
        (37, 40, "ATTG", "#1420 v3: reference 37-40"),
        (38, 41, "TTGC", "#1420 v2: reference 38-41"),
    ] {
        assert_eq!(
            &TEMPLATE[first - 1..last],
            expected,
            "{source} should be {expected}"
        );
    }
}

/// Every spelling still prints what the reports recorded for it.
///
/// The eighteen pins. A failure here is a representation change on a shape a
/// downstream consumer keys read counts on, so it is meant to be argued for in
/// review rather than re-blessed: the message names the verdict and the wanted
/// form so a reviewer can tell a fix from a regression without reopening the
/// issue.
#[test]
fn every_reported_spelling_still_normalizes_as_the_reports_recorded() {
    for row in REPORTED_ROWS {
        assert_eq!(
            normalized(row.input),
            row.output,
            "{} moved.\n  input:   {}\n  pinned:  {}\n  verdict: {:?}\n  {}",
            row.label,
            row.input,
            row.output,
            row.verdict,
            row.argument
        );
    }
}

/// Every spelling's 5' answer, pinned as hard as its 3' one.
///
/// Nothing pinned this before: `reported_confluence_pairs` sweeps both
/// directions but only *counts* convergence, and this module pinned 3' alone.
/// A count cannot tell a reader that flipping the direction moves two of the
/// eighteen rows, and moving a normalized string is a migration for whoever
/// stored it whichever direction produced it.
#[test]
fn every_reported_spelling_normalizes_as_recorded_under_five_prime() {
    for row in REPORTED_ROWS {
        assert_eq!(
            normalize_in(TEMPLATE, row.input, ShuffleDirection::FivePrime),
            row.five_prime,
            "{} moved under FivePrime.\n  input:   {}\n  pinned:  {}\n  \
             (its ThreePrime answer is {})",
            row.label,
            row.input,
            row.five_prime,
            row.output
        );
    }
}

/// Exactly two rows answer differently in the two directions, and they are the
/// two [`FIVE_PRIME_MOVERS`] names.
///
/// Asserted separately from the pins above so the *set* cannot drift silently:
/// with only the per-row pins, editing a `five_prime` string to match a new
/// behaviour would keep the suite green while quietly changing which rows are
/// direction-sensitive — and direction-sensitivity is the property
/// [`every_gap_row_is_returned_exactly_as_authored`]'s argument depends on.
#[test]
fn only_the_named_rows_answer_differently_in_the_two_directions() {
    let movers: Vec<&str> = REPORTED_ROWS
        .iter()
        .filter(|row| row.output != row.five_prime)
        .map(|row| row.label)
        .collect();
    assert_eq!(
        movers, FIVE_PRIME_MOVERS,
        "the set of direction-sensitive rows changed; if a row joined it, the \
         5' pass now moves a spelling the 3' pass returns verbatim, which is a \
         representation change under a supported option"
    );
}

/// Every gap row's `wanted` form is the one its sibling spelling already prints.
///
/// This is the family's defining shape stated as a checked fact rather than as
/// prose in an `argument` field. Each issue names an expected output; that
/// output is, in every one of the nine pairs, exactly what ferro produces from
/// the *other* spelling of the same variant. That is what makes the defect
/// spelling-dependence rather than a missing rule — the canonical form is
/// already reachable, just not from here.
///
/// Written as a comparison against the sibling's pinned `output` rather than
/// against a recomputed value, so a change that moved *both* spellings onto some
/// third string would fail here instead of quietly redefining "wanted".
#[test]
fn the_wanted_form_of_every_gap_row_is_its_siblings_output() {
    for (a, b) in reported_pairs() {
        // Exhaustive on purpose: the claim is statable in exactly one of the
        // three states, and adding a fourth must force that decision rather
        // than fall through a `contains` check.
        match pinned_state(pair_of(a.label)) {
            // One canonical spelling and one gap — the shape this test is about.
            PairState::OneReaches => {}
            // Both spellings reach it, so there is no gap row to make the claim
            // about. Vacuous, not skipped:
            // `each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names`
            // is what measures such a pair.
            PairState::BothReach => continue,
            // Neither spelling reaches the wanted form: both rows are `Gap`, so
            // there is no canonical sibling whose output could be compared. Not
            // a convergence claim — see `PairState::NeitherReaches`.
            PairState::NeitherReaches => continue,
        }
        let (canonical, gap) = match a.verdict {
            Verdict::Canonical => (a, b),
            Verdict::Gap => (b, a),
        };
        assert_eq!(
            gap.wanted, canonical.output,
            "{}: the issue asks for `{}`, but the sibling spelling {} prints \
             `{}`. Either the wanted form was transcribed wrong or the family \
             has stopped being a pure partition disagreement.",
            gap.label, gap.wanted, canonical.label, canonical.output
        );
        assert_ne!(
            gap.output, gap.wanted,
            "{}: this row is marked Gap but already prints the wanted form; \
             flip its verdict to Canonical and lower OPEN_GAPS",
            gap.label
        );
    }
}

/// Every canonical row already prints the form its issue asks for.
///
/// The other half of the same claim, and the one that makes `wanted` load-bearing
/// rather than decorative: for seven of the eighteen spellings ferro and the
/// reporter agree outright — the other eleven are [`OPEN_GAPS`] — and that
/// agreement is what a fix must not disturb.
#[test]
fn every_canonical_row_already_prints_its_wanted_form() {
    for row in REPORTED_ROWS {
        if row.verdict != Verdict::Canonical {
            continue;
        }
        assert_eq!(
            row.output, row.wanted,
            "{}: marked Canonical but prints something other than the form its \
             issue names",
            row.label
        );
    }
}

/// Both spellings of a pair carry the same [`Authority`], and the split across
/// the family is the one the module docs argue for.
///
/// Pinned because the authority is what decides whether a row's `wanted` is a
/// *target* or merely a *record*, and that is the single most consequential
/// judgement in this file. Three pairs (#1420's) rest on a spec line nothing
/// contradicts; six (#1419's and #1421's) rest on a spec that answers both
/// ways — but the two self-conflicting families no longer sit the same way
/// under the policy, so stating them together is what made this sentence wrong.
///
/// **#1421's three pairs** still carry the issue's split form as `wanted`,
/// recorded rather than targeted: the policy falls through to Mutalyzer, and the
/// measurement in this module's docs has Mutalyzer converging that shape on the
/// **merged** delins, so it does not pick the issue's form.
///
/// **#1419's three pairs no longer carry the issue's form at all.** The decided
/// chain — `adjudication-precedence-order` withdrawing `general.md:56`, then
/// `canonical-form-choice-when-both-legal` with `delins.md:46-47` — moved all
/// six `wanted` strings onto the spanning delins, so the Mutalyzer tie-break is
/// never reached there. Note the same measurement has Mutalyzer converging that
/// shape on the spanning delins too, i.e. **agreeing** with `wanted` — the
/// opposite of "picks neither". Ferro prints that form for neither spelling,
/// which is why both rows of each `1419-r*` pair are `Gap`.
#[test]
fn the_spec_authority_census_holds() {
    for (a, b) in reported_pairs() {
        assert_eq!(
            a.authority, b.authority,
            "{} and {} are two spellings of one shape, so they cannot have \
             different authorities",
            a.label, b.label
        );
    }

    let explicit = REPORTED_ROWS
        .iter()
        .filter(|row| row.authority == Authority::SpecExplicit)
        .count();
    let conflicting = REPORTED_ROWS
        .iter()
        .filter(|row| row.authority == Authority::SpecSelfConflicting)
        .count();
    assert_eq!(
        (explicit, conflicting),
        (6, 12),
        "the authority split moved. Six rows (#1420's) rest on delins.md:16/17 \
         with nothing in the spec against them — they are length-neutral, so \
         reference and payload correspond column for column and delins.md:44-47 \
         never engages. Twelve (#1419's six and #1421's six) replace a span with \
         a payload of a different length, so their wanted split is recoverable \
         only by choosing an alignment, which is what :44-47 says to write as one \
         delins while :17 says to split. Moving a row between the two is a \
         decision about what a fix is allowed to target, so it belongs in a PR \
         description."
    );

    // Every self-conflicting row belongs to #1419 or #1421, and no other row
    // does. Stated as a set rather than left to the counts above, which two
    // compensating edits could satisfy while relabelling the wrong rows.
    let conflicting_labels: Vec<&str> = REPORTED_ROWS
        .iter()
        .filter(|row| row.authority == Authority::SpecSelfConflicting)
        .map(|row| row.label)
        .collect();
    assert!(
        conflicting_labels
            .iter()
            .all(|l| l.starts_with("1419-") || l.starts_with("1421-")),
        "only #1419's and #1421's rows have a spec that answers both ways; got \
         {conflicting_labels:?}"
    );
}

/// Each reported pair reaches its `wanted` form from exactly the spellings its
/// pinned [`PairState`] names — **measured**, not read off the `verdict` column.
///
/// This is the defect in one assertion, and it is sharper than "the two
/// spellings disagree". Both spellings denote one variant; the canonical form
/// is a property of the variant, so it should be reachable from either. For six
/// of the nine pairs it is reachable from exactly one, which makes reachability
/// a function of how the variant was *written* — the spelling-dependence
/// root-caused on #1419 to a weight bound whose threshold is the input's own
/// spelling.
///
/// The reach set is measured per row against that row's **own** `wanted`, never
/// against the sibling's `output`, so the claim reads the same in all three
/// states. Where a canonical sibling exists the two strings are the same one —
/// that identity is what
/// [`the_wanted_form_of_every_gap_row_is_its_siblings_output`] pins — and where
/// none exists ([`PairState::NeitherReaches`]) only the `wanted` formulation
/// says anything at all. The sibling-relative form is kept as a second,
/// differently-derived assertion in the `OneReaches` arm rather than dropped.
///
/// This test is written so that **every** transition is loud and none is
/// automatic. It replaces an `assert_ne!(a.verdict, b.verdict)` that could not
/// be satisfied by a fixed pair: the remedy its own failure message prescribed
/// was to make both rows `Canonical`, which that assertion then rejected. A pair
/// reaching its wanted form from both spellings now fails *this* test, naming
/// [`PairState::BothReach`] as the entry to move it to — a census edit that has
/// to be argued for, not a green run.
#[test]
fn each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names() {
    for (a, b) in reported_pairs() {
        let pair = pair_of(a.label);
        let pinned = pinned_state(pair);

        let reaching: Vec<&str> = [a, b]
            .into_iter()
            .filter(|row| normalized(row.input) == row.wanted)
            .map(|row| row.label)
            .collect();
        let measured = match reaching.as_slice() {
            [] => PairState::NeitherReaches,
            [_] => PairState::OneReaches,
            _ => PairState::BothReach,
        };

        assert_eq!(
            measured, pinned,
            "{pair}: pinned {pinned:?} in PAIR_STATES, measures {measured:?} \
             (reached by {reaching:?}).\n  \
             A pair changing state is a deliberate multi-place edit: move this \
             pair's PAIR_STATES entry, adjust PAIR_STATE_CENSUS, flip the \
             affected rows' verdicts, and set OPEN_GAPS to match. BothReach is \
             the fixed state and is sayable — but say in the PR WHICH CLAUSE \
             carried the move. A pair that arrived there by re-derivation \
             rather than by a licensed fix must not be banked as a fix; see \
             OPEN_GAPS."
        );

        // Per row, so the `verdict` column cannot drift away from what the row
        // measures while the pair-level count stays right — two rows swapping
        // verdicts leaves `measured` unchanged.
        for row in [a, b] {
            assert_eq!(
                normalized(row.input) == row.wanted,
                row.verdict == Verdict::Canonical,
                "{}: marked {:?}, but reaching `{}` is {}. A row is Canonical \
                 exactly when it prints the form it is measured against.",
                row.label,
                row.verdict,
                row.wanted,
                normalized(row.input) == row.wanted
            );
        }

        if pinned == PairState::OneReaches {
            let (canonical, gap) = match a.verdict {
                Verdict::Canonical => (a, b),
                Verdict::Gap => (b, a),
            };

            assert_eq!(
                normalized(canonical.input),
                canonical.output,
                "{}: the canonical spelling stopped reaching the canonical form",
                canonical.label
            );
            // The sibling-relative half of the same claim, kept because it is
            // derived through a different string: `canonical.output` rather
            // than `gap.wanted`. If the two ever come apart, the pin above
            // fires here and not only in the test that owns that identity.
            assert_ne!(
                normalized(gap.input),
                canonical.output,
                "{} now reaches `{}` too, so the pair has converged on the form \
                 its issue argues for. That is the fix, and it is now a STATE \
                 rather than a contradiction: move {} to PairState::BothReach, \
                 flip this row's verdict to Canonical, lower OPEN_GAPS, and say \
                 in the PR which stored spelling moved and which clause carried \
                 it.",
                gap.label,
                canonical.output,
                pair
            );
        }
    }
}

/// The pair-state census, so no pair can change state in silence.
///
/// Three claims, each of which a different mistake trips. The **per-pair**
/// entries tie [`PAIR_STATES`] to the `verdict` column, so a verdict edited
/// without its census entry fails. The **family tuple** [`PAIR_STATE_CENSUS`]
/// catches two compensating per-pair edits, which every individual entry would
/// otherwise accept. And the **arithmetic tie to [`OPEN_GAPS`]** — 0, 1 or 2 gap
/// rows per pair by state — is derived from the pair states where
/// [`the_open_gap_census_holds`] derives the same number by counting rows, so
/// the two censuses have to agree without either being computed from the other.
///
/// Together with
/// [`each_pair_reaches_its_wanted_form_from_the_spellings_its_state_names`],
/// which measures the same states against ferro's actual output, that closes the
/// chain: measurement, verdict column, per-pair census, family census, row
/// count. Nothing in it moves on its own.
#[test]
fn the_pair_state_census_holds() {
    assert_eq!(
        PAIR_STATES.len() * 2,
        REPORTED_ROWS.len(),
        "PAIR_STATES must carry exactly one entry per reported pair"
    );

    let (mut both, mut one, mut neither) = (0usize, 0usize, 0usize);
    let mut gap_rows = 0usize;

    for ((name, pinned), (a, b)) in PAIR_STATES.iter().zip(reported_pairs()) {
        assert_eq!(
            *name,
            pair_of(a.label),
            "PAIR_STATES and REPORTED_ROWS list the reported pairs in different \
             orders, so this row's state is being read against the wrong pair"
        );

        let derived = PairState::from_verdicts(a.verdict, b.verdict);
        assert_eq!(
            derived, *pinned,
            "{name}: pinned {pinned:?}, but its two verdicts ({:?}, {:?}) make \
             it {derived:?}. Whichever is wrong, both change together — and if \
             a spelling really did start reaching its wanted form, name the \
             clause that carried it (see OPEN_GAPS).",
            a.verdict, b.verdict
        );

        match pinned {
            PairState::BothReach => both += 1,
            PairState::OneReaches => one += 1,
            PairState::NeitherReaches => neither += 1,
        }
        gap_rows += pinned.gap_rows();
    }

    assert_eq!(
        (both, one, neither),
        PAIR_STATE_CENSUS,
        "the family's shape moved: (BothReach, OneReaches, NeitherReaches) is \
         now ({both}, {one}, {neither}). A rise in BothReach is the family \
         being fixed and must name the clause; a fall is a regression."
    );

    assert_eq!(
        gap_rows, OPEN_GAPS,
        "the two censuses disagree: the pair states imply {gap_rows} gap rows \
         and OPEN_GAPS says {OPEN_GAPS}. They are derived differently on \
         purpose — this one from PAIR_STATES, the_open_gap_census_holds by \
         counting rows — so one was moved without the other."
    );
}

/// Every gap row is returned exactly as it was authored.
///
/// Deliberately *not* the same claim as `reported_confluence_pairs`'
/// `every_reported_output_is_a_fixed_point`, which normalizes twice and checks
/// the second pass changes nothing. That would also be satisfied by a row whose
/// first pass moved it somewhere new and then settled.
///
/// This asserts the stronger thing the reports actually observed: the input is
/// *already* a fixed point — ferro does not touch it. That distinction is what
/// rules out the benign reading of this family. A non-canonical spelling that
/// moved, even to some third string, would be an under-applied pass. One that is
/// returned untouched is a second canonical form, and a second canonical form is
/// what splits a consumer's counts across two keys.
///
/// **This is a 3'-direction claim and does not generalize.** Three gap rows do
/// move under 5' ([`FIVE_PRIME_MOVERS`]), so under that option the argument
/// above reaches fewer rows than it does under 3' — which is exactly why the 5'
/// answers are now pinned per row rather than left to a convergence count.
///
/// **And it does not reach every gap row even under 3'.** There are eleven gap
/// rows now ([`OPEN_GAPS`]), and the three `1419-r*/span` rows are gap rows that
/// ferro *splits* rather than returns — see [`PairState::NeitherReaches`], whose
/// pairs this test exempts from the retention claim while still pinning each
/// row's measured output.
#[test]
fn every_gap_row_is_returned_exactly_as_authored() {
    for row in REPORTED_ROWS {
        if row.verdict != Verdict::Gap {
            continue;
        }
        if pinned_state(pair_of(row.label)) == PairState::NeitherReaches {
            // A gap row in one of these pairs may also MOVE: the span spelling
            // leaves its input for the del-plus-sub form, which is not the
            // wanted one. Retention is therefore not the property to assert for
            // the whole pair — but exempting the pair must not silently drop the
            // pin from the `/cis` half, which DOES still retain its input. So
            // pin the MEASURED output, which holds either way, and keep the
            // wanted-form check that would close the pair.
            //
            // The exemption is keyed on the pinned state rather than on a
            // separate list, but it is the same set of pairs and the same
            // weaker pin — see `PairState::NeitherReaches`.
            let output = normalized(row.input);
            assert_eq!(
                output, row.output,
                "{}: measured output moved off its pinned value (`{}` -> `{}`)",
                row.label, row.output, output
            );
            assert_ne!(
                output, row.wanted,
                "{}: now reaches its wanted form; move this pair's PAIR_STATES \
                 entry off NeitherReaches and name the clause that carried it",
                row.label
            );
            continue;
        }
        let output = normalized(row.input);
        assert_eq!(
            output, row.input,
            "{}: the reported spelling is no longer retained verbatim. If it \
             moved to the canonical form the family is fixed; if it moved \
             somewhere else this is a new representation nobody asked for.",
            row.label
        );
    }
}

/// `EquivalenceChecker` still calls each reported pair one variant.
///
/// All three issues state this, and it is the load-bearing mitigation rather
/// than a detail: it is what makes the divergence a representation problem
/// instead of a correctness one, and it is the fallback the issues tell
/// consumers to compare on while the family is open. If it regressed, the
/// reports' "use `EquivalenceChecker` instead" advice would quietly stop
/// working — and nothing else in the suite covers these shapes.
///
/// The exact level, not merely `is_equivalent()`, and it is **not one level for
/// every pair**: a pair whose two spellings normalize to one string is
/// `NormalizedMatch`, and a pair whose spellings still differ stops at a
/// **denotational** rung. Pinning the stronger level where it holds is what
/// stops a convergence from being lost, and pinning the weaker one where it
/// holds is what stops a convergence from being claimed.
///
/// The denotational rung is `CrossAxisSequenceMatch` rather than
/// `SequenceMatch` because these are genomic descriptions, which determine
/// exactly one axis — the genome — so apply-equality on it *is* the whole
/// equivalence relation and there is no second axis left to disagree. (A `c.`
/// pair also determines the genome through its exon alignment, and there the
/// two rungs come apart: `LRG_199t1:c.3921dup` and `c.3922dup` agree on the
/// transcript and not on the genome, so they stop at `SequenceMatch`.)
#[test]
fn every_reported_pair_is_still_one_variant_by_equivalence() {
    let checker = EquivalenceChecker::new(provider(TEMPLATE));
    for (a, b) in reported_pairs() {
        // Derived from the rows' own measured outputs rather than from a second
        // list: the level IS whether the two spellings print one string, so a
        // hand-maintained set of pair names could disagree with the pins above
        // and would then be measuring itself.
        //
        // The three `1419-r*` pairs are `NormalizedMatch` as of #1649, which
        // gave the splitter the two-deletion shape that lets each span spelling
        // reach its `/cis` sibling's form. They stay `PairState::NeitherReaches`:
        // converging on one string is not reaching `wanted`, and that census
        // tracks the target.
        //
        // Every other pair takes the denotational rung. `PAIR_STATES` says
        // nothing about convergence in either direction — a `NeitherReaches`
        // pair may or may not converge, and today all three do — so there is no
        // exemption to take here either.
        let expected = if a.output == b.output {
            EquivalenceLevel::NormalizedMatch
        } else {
            EquivalenceLevel::CrossAxisSequenceMatch
        };
        let left = parse_hgvs(a.input).unwrap_or_else(|e| panic!("{}: {e}", a.label));
        let right = parse_hgvs(b.input).unwrap_or_else(|e| panic!("{}: {e}", b.label));
        let result = checker
            .check(&left, &right)
            .unwrap_or_else(|e| panic!("{} vs {}: {e}", a.label, b.label));
        assert_eq!(
            result.level, expected,
            "{} vs {}: expected {expected:?}, got {:?}. If this is now \
             NormalizedMatch the pair has converged and the pinned rows are \
             stale; anything else is a regression in the fallback the reports \
             direct consumers to.",
            a.label, b.label, result.level
        );
    }
}

/// `(label, first changed position, reference bases, replacement)` for the
/// three #1420 rows whose replacement is the same length as the span it
/// replaces.
///
/// Equal length is what makes the changed/unchanged column pattern well defined
/// without an alignment, which is why only #1420 appears here: #1419's rows are
/// net deletions and #1421's are net insertions, so their columns do not line
/// up one-to-one and the pattern is a choice rather than a fact.
const EQUAL_LENGTH_SPANS: &[(&str, usize, &str, &str)] = &[
    ("1420-v2", 38, "TTGC", "ATTG"),
    ("1420-v3", 37, "ATTG", "CATT"),
    ("1420-v4", 21, "ATGC", "GCTG"),
];

/// The #1420 arguments are checked against the reference, not taken on trust.
///
/// Each row's canonical form follows from *which columns change*: v2 and v3
/// have an unchanged interior column, so delins.md:17 separates them; v4 has
/// none, so delins.md:16 coalesces it. That distinction is the reason v4's
/// correction points the opposite way to v2's and v3's, and it is the one part
/// of #1420's reasoning that is a fact about the sequence rather than a reading
/// of the recommendations — so it is computed here rather than asserted in a
/// comment.
#[test]
fn reported_spans_change_the_columns_the_reports_state() {
    // `X` = changed, `.` = unchanged, one character per position in the span.
    let expected: &[(&str, &str)] = &[
        ("1420-v2", "X.XX"),
        ("1420-v3", "XX.X"),
        ("1420-v4", "XXXX"),
    ];
    // `zip` stops at the shorter side, so a span added to one table and not the
    // other would be skipped rather than reported.
    assert_eq!(
        EQUAL_LENGTH_SPANS.len(),
        expected.len(),
        "every equal-length span needs a column pattern"
    );

    for ((label, first, reference, replacement), (echo, pattern)) in
        EQUAL_LENGTH_SPANS.iter().zip(expected)
    {
        assert_eq!(label, echo, "the two tables fell out of order");
        assert_eq!(
            reference.len(),
            replacement.len(),
            "{label}: this table is only meaningful for equal-length spans"
        );
        assert_eq!(
            &TEMPLATE[first - 1..first - 1 + reference.len()],
            *reference,
            "{label}: reference bases at {first} are not what the report states"
        );

        let observed: String = reference
            .bytes()
            .zip(replacement.bytes())
            .map(|(r, a)| if r == a { '.' } else { 'X' })
            .collect();
        assert_eq!(
            &observed, pattern,
            "{label}: {reference} -> {replacement} changes columns {observed}, \
             not {pattern} as the report states"
        );
    }

    // The consequence, stated so it is not left implicit: two rows have an
    // unchanged interior column and one does not, which is exactly why #1420's
    // three corrections cannot share a per-member-type rule.
    let separable = expected
        .iter()
        .filter(|(_, pattern)| pattern.trim_matches('X').contains('.'))
        .count();
    assert_eq!(
        separable, 2,
        "expected v2 and v3 to be separable under delins.md:17 and v4 to be \
         one run under delins.md:16"
    );
}

/// `(label, last position of the leading substitution, first position of the
/// trailing delins)` for #1421's three split spellings.
const ISSUE_1421_MEMBER_BOUNDS: &[(&str, usize, usize)] = &[
    ("1421-n1", 29, 32),
    ("1421-n2", 32, 35),
    ("1421-n3", 34, 37),
];

/// #1421's rows fall under delins.md:17, not under the delins.md:18 exception.
///
/// Worth its own test because the exception is the only reading under which the
/// retained spans would be correct, and it turns on a number: it covers two
/// variants separated by **one** nucleotide that together affect one amino
/// acid. Every #1421 row is separated by two, and is a `g.` description with no
/// amino acid, so it fails the exception twice over.
///
/// The interior positions are unchanged by construction — they lie strictly
/// between the split spelling's two members and are named by neither, so a cis
/// allele leaves them alone. What needs checking is the width.
///
/// The two member positions are read back out of the corresponding row in
/// [`REPORTED_ROWS`] before the arithmetic runs. Without that the test would be
/// `32 - 29 - 1 == 2` over its own constants — true of the table rather than of
/// the reported variants, and still green if somebody moved a row's coordinates.
#[test]
fn the_1421_spans_separate_by_two_nucleotides_not_one() {
    for (label, substitution, delins_start) in ISSUE_1421_MEMBER_BOUNDS {
        let split_label = format!("{label}/split");
        let row = REPORTED_ROWS
            .iter()
            .find(|row| row.label == split_label)
            .unwrap_or_else(|| panic!("no row labelled `{split_label}`"));

        // The split spelling is `<acc>:g.[<sub-position><ref>><alt>;<delins-start>_...]`,
        // so these two anchors pin the bounds below to the row's own coordinates.
        assert!(
            row.input.contains(&format!("g.[{substitution}")),
            "{label}: the bounds table says the substitution is at \
             {substitution}, but the row's first member is `{}`",
            row.input
        );
        assert!(
            row.input.contains(&format!(";{delins_start}_")),
            "{label}: the bounds table says the delins starts at \
             {delins_start}, but the row's second member is `{}`",
            row.input
        );

        let separation = delins_start - substitution - 1;
        assert_eq!(
            separation, 2,
            "{label}: members at {substitution} and {delins_start} are \
             separated by {separation} nt; the delins.md:18 exception covers \
             one, so a separation of two puts this row squarely under \
             delins.md:17"
        );
    }
}

/// The two reported-family modules describe the same nine pairs, on the same
/// bases.
///
/// Each module holds its own copy of the contig and of all eighteen spellings,
/// and the [`Row::label`] field promises the rows line up with
/// `reported_confluence_pairs`. Nothing enforced that promise, so editing one
/// table would leave the two modules silently measuring different variants while
/// both stayed green — this one pinning verdicts for one set of spellings and
/// `reported_confluence_pairs` ratcheting convergence for another. Since a fix
/// is expected to move rows in *both* modules together, the two falling out of
/// step is a live hazard rather than a theoretical one.
#[test]
fn the_two_reported_modules_describe_the_same_pairs() {
    assert_eq!(
        TEMPLATE,
        crate::reported_confluence_pairs::TEMPLATE,
        "the two modules' contigs have diverged, so their rows are no longer \
         about the same bases"
    );

    let pairs = crate::reported_confluence_pairs::REPORTED_PAIRS;
    assert_eq!(
        pairs.len() * 2,
        REPORTED_ROWS.len(),
        "one module has gained or lost a reported pair"
    );

    for ((label, a, b), (row_a, row_b)) in pairs.iter().zip(reported_pairs()) {
        assert_eq!(
            *label,
            pair_of(row_a.label),
            "the two modules list the reported pairs in different orders"
        );
        assert_eq!(
            [*a, *b],
            [row_a.input, row_b.input],
            "{label}: the two modules disagree about how this pair is spelled"
        );
    }
}

/// The census, so the count and the table cannot drift apart.
#[test]
fn the_open_gap_census_holds() {
    let gaps = REPORTED_ROWS
        .iter()
        .filter(|row| row.verdict == Verdict::Gap)
        .count();
    assert_eq!(
        gaps, OPEN_GAPS,
        "the number of rows disagreeing with their cited recommendation moved \
         (now {gaps}). Down is progress: lower OPEN_GAPS and name the moved \
         representation in the PR. Up means a change pushed a row away from \
         the spec."
    );
    assert_eq!(
        REPORTED_ROWS.len(),
        18,
        "REPORTED_ROWS must stay two spellings per reported pair"
    );
}
