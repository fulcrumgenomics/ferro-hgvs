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
//! 3. that `EquivalenceChecker` still calls the two one variant
//!    (`SequenceMatch`) — which is what makes the divergence a *representation*
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
//! Nine of the eighteen rows below print something the cited recommendation
//! argues against. They are pinned **as they behave today**, with the wanted
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
//! Every row carries a [`Row::wanted`] string — the output *its issue asks for*,
//! as the issue states it — and an [`Authority`] recording what backs it. The
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
//! answer, which nothing pinned before. Sixteen of the eighteen agree across
//! directions; the two that do not are named in
//! [`FIVE_PRIME_MOVERS`], and both are `Gap` rows that the 3' pass returns
//! verbatim and the 5' pass does not — see
//! [`every_gap_row_is_returned_exactly_as_authored`] for why that distinction
//! carries weight.

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
    /// Mutalyzer — which does not pick the issue's form either. The wanted form
    /// is recorded as what the issue asked for, not asserted as the target.
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
    /// Equal to [`Row::output`] for sixteen of the eighteen rows; the two that
    /// differ are listed in [`FIVE_PRIME_MOVERS`].
    five_prime: &'static str,
    verdict: Verdict,
    /// The output **the issue asks for**, in the issue's own terms — the string
    /// its "canonical"/"expected" column names for this spelling. For a
    /// [`Verdict::Canonical`] row this equals [`Row::output`]; for a
    /// [`Verdict::Gap`] row it is the sibling spelling's output, which is the
    /// family's defining shape and is asserted rather than assumed by
    /// [`the_wanted_form_of_every_gap_row_is_its_siblings_output`].
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
/// [`each_pair_reaches_its_canonical_form_from_exactly_one_spelling`].
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
        wanted: "TEMPLATE:g.[19_30del;33T>G]",
        authority: Authority::SpecSelfConflicting,
        argument: "Wanted `[19_30del;33T>G]`. general.md:56 prioritises \
                   (1) substitution over (2) deletion, so the change at 33 \
                   should be exposed as a substitution rather than absorbed \
                   into a deletion. The cis spelling hides it; the equivalent \
                   span reaches it (see `1419-r1/span`), so the partition is \
                   reachable and simply is not reached from here. \
                   BUT the spec also answers this shape the other way: the \
                   wanted form is a deletion plus a substitution arising from \
                   payload bases that align with the reference, which is the \
                   case delins.md:44-47 covers and calls out — `The \"delins\" \
                   format is recommended`. So the spec conflicts with itself \
                   and this row is recorded, not targeted.",
    },
    Row {
        label: "1419-r1/span",
        input: "TEMPLATE:g.19_33delinsCGG",
        output: "TEMPLATE:g.[19_30del;33T>G]",
        five_prime: "TEMPLATE:g.[19_30del;33T>G]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[19_30del;33T>G]",
        authority: Authority::SpecSelfConflicting,
        argument: "The form #1419 names canonical: re-derived from the \
                   resulting sequence, exposing the substitution at 33 per \
                   general.md:56, with the retained bases placed as 3' as \
                   possible per general.md:41. Note ferro is *leaving* the \
                   spanning delins to reach it, and delins.md:44-47 recommends \
                   the spanning delins for exactly this alignment coincidence \
                   — a tension ferro's own coalesce pass encodes (see \
                   `merge::the_coalesce_pass_reaches_the_spec_worked_example`, \
                   which restores the single delins on the spec's own example).",
    },
    Row {
        label: "1419-r2/cis",
        input: "TEMPLATE:g.[19_22del;26_36del]",
        output: "TEMPLATE:g.[19_22del;26_36del]",
        five_prime: "TEMPLATE:g.[19_22del;26_36del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.[19_33del;36A>G]",
        authority: Authority::SpecSelfConflicting,
        argument: "Wanted `[19_33del;36A>G]`, for the same general.md:56 \
                   reason as `1419-r1/cis`, one locus over — and subject to \
                   the same delins.md:44-47 counter-reading.",
    },
    Row {
        label: "1419-r2/span",
        input: "TEMPLATE:g.19_36delinsGCG",
        output: "TEMPLATE:g.[19_33del;36A>G]",
        five_prime: "TEMPLATE:g.[19_33del;36A>G]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[19_33del;36A>G]",
        authority: Authority::SpecSelfConflicting,
        argument: "The form #1419 names canonical for row 2.",
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
        wanted: "TEMPLATE:g.[19T>G;22_33del]",
        authority: Authority::SpecSelfConflicting,
        // TRAP, see OPEN_GAPS. The wanted form here is ALSO what a
        // re-partitioning of the block from the sequence produces: two
        // deletions three unchanged nucleotides apart (25, 26, 27) coming back
        // as a substitution plus a deletion. So this row can close for the
        // wrong reason, and closing it for the wrong reason looks identical to
        // closing it for the right one. Do not flip the verdict without naming
        // the clause that licensed the merge.
        argument: "Wanted `[19T>G;22_33del]`. Same general.md:56 reason, but \
                   note the exposed substitution lands at the 5' end here \
                   rather than the 3' end — which is why #1419 argues no \
                   per-end rule fixes the family. Same delins.md:44-47 \
                   counter-reading as the other two rows.",
    },
    Row {
        label: "1419-r3/span",
        input: "TEMPLATE:g.19_33delinsGGA",
        output: "TEMPLATE:g.[19T>G;22_33del]",
        five_prime: "TEMPLATE:g.[19T>G;22_33del]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[19T>G;22_33del]",
        authority: Authority::SpecSelfConflicting,
        argument: "The form #1419 names canonical for row 3.",
    },
    // -- #1420 -- one row per member type. v2/v3 must reduce, v4 must coalesce:
    // the corrections point in opposite directions, which is the issue's
    // argument that no per-member-type rule reconciles them.
    Row {
        label: "1420-v2/cis",
        input: "TEMPLATE:g.[37dup;41del]",
        output: "TEMPLATE:g.[37dup;41del]",
        // The 5' pass moves the duplication's anchor one base left, into the
        // `AA` at 36-37. `dup` is 3'-anchored by general.md:41 in the shipped
        // direction, so this is the one row where the reported spelling is not
        // even a fixed point once the direction is flipped.
        five_prime: "TEMPLATE:g.[36dup;41del]",
        verdict: Verdict::Gap,
        wanted: "TEMPLATE:g.[38T>A;40_41delinsTG]",
        authority: Authority::SpecExplicit,
        // TRAP, see OPEN_GAPS. This row's members are three unchanged
        // nucleotides apart (38, 39, 40), and the wanted form is what a
        // re-partitioning across them produces — so the row can converge on
        // #1420's own ask by violating `general.md:34`. Observed on an
        // experimental partitioner arm, where it raised
        // `reported_confluence_pairs`' 5' census from 0 to 1 on this row alone.
        // `the_1420_v2_pair_does_not_converge_by_re_derivation` pins the string
        // as forbidden; read it before flipping this verdict.
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
        label: "1421-n1/split",
        input: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        output: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        five_prime: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
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
        label: "1421-n1/span",
        input: "TEMPLATE:g.29_33delinsAACACATACTG",
        output: "TEMPLATE:g.29_33delinsAACACATACTG",
        five_prime: "TEMPLATE:g.29_33delinsAACACATACTG",
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
        label: "1421-n2/split",
        input: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        output: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        five_prime: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        authority: Authority::SpecSelfConflicting,
        argument: "The form #1421 names canonical for row 2; unchanged \
                   interior 33-34. Recorded, not targeted — delins.md:44-47 \
                   answers the other way (see the sibling span row).",
    },
    Row {
        label: "1421-n2/span",
        input: "TEMPLATE:g.32_36delinsTTGGAATCGAC",
        output: "TEMPLATE:g.32_36delinsTTGGAATCGAC",
        five_prime: "TEMPLATE:g.32_36delinsTTGGAATCGAC",
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
        label: "1421-n3/split",
        input: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        output: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        five_prime: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        verdict: Verdict::Canonical,
        wanted: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        authority: Authority::SpecSelfConflicting,
        argument: "The form #1421 names canonical for row 3; unchanged \
                   interior 35-36. Recorded, not targeted — delins.md:44-47 \
                   answers the other way (see the sibling span row).",
    },
    Row {
        label: "1421-n3/span",
        input: "TEMPLATE:g.34_38delinsTCACCTTTACG",
        output: "TEMPLATE:g.34_38delinsTCACCTTTACG",
        five_prime: "TEMPLATE:g.34_38delinsTCACCTTTACG",
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

/// The two rows whose 5' answer differs from their 3' answer.
///
/// Named rather than left implicit because both are `Gap` rows, and
/// [`every_gap_row_is_returned_exactly_as_authored`] rests on a gap row being
/// returned *untouched* — which is what makes it a second canonical form rather
/// than an under-applied pass. That argument is a 3'-direction argument, and
/// these two are where it stops holding.
const FIVE_PRIME_MOVERS: &[&str] = &["1419-r3/cis", "1420-v2/cis"];

/// How many of [`REPORTED_ROWS`] print something their issue argues against.
///
/// Nine of eighteen — exactly one per reported pair, which is the family's
/// shape rather than a coincidence; see
/// [`each_pair_reaches_its_canonical_form_from_exactly_one_spelling`].
///
/// **This may only go down.** Up means a change has moved a row *away* from the
/// cited recommendation.
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
/// The second has been observed, on an experimental partitioner arm, on exactly
/// the two rows below whose `wanted` form is reachable that way: `1419-r3/cis`
/// (`g.[19_24del;28_33del]` -> `g.[19T>G;22_33del]`, two deletions three
/// unchanged nucleotides apart coming back as a substitution plus a deletion)
/// and `1420-v2/cis` (`g.[37dup;41del]` -> `g.[38T>A;40_41delinsTG]`, likewise
/// three unchanged nucleotides apart). Both are `general.md:34` violations that
/// land on an issue's wanted string, which is precisely how a defect gets banked
/// as a fix.
///
/// So lowering this constant requires naming, in the PR, **which clause carried
/// the move** — not merely that the row now prints its wanted form. A row that
/// converged with no clause behind it is a re-derivation and the count must not
/// move. `reported_confluence_pairs::the_1420_v2_pair_does_not_converge_by_re_derivation`
/// pins the `1420-v2` half of that as a forbidden string.
const OPEN_GAPS: usize = 9;

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
/// [`each_pair_reaches_its_canonical_form_from_exactly_one_spelling`] still
/// passes — measured, not assumed: every assertion in it compares a row against
/// its own pinned output, so a mispaired neighbour changes nothing it looks at.
/// Only [`every_reported_pair_is_still_one_variant_by_equivalence`] notices, and
/// it reports the mispairing as `NotEquivalent` — that is, as a regression in
/// the equivalence fallback, which is the wrong diagnosis and sends a reader to
/// the wrong file.
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
/// rather than decorative: for nine of the eighteen spellings ferro and the
/// reporter agree outright, and that agreement is what a fix must not disturb.
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
/// ways, where the precedence policy hands the tie-break to Mutalyzer and
/// Mutalyzer picks neither ferro's answer nor the issue's.
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

/// Each reported pair reaches its canonical form from exactly one spelling.
///
/// This is the defect in one assertion, and it is sharper than "the two
/// spellings disagree". Both spellings denote one variant; the canonical form
/// is a property of the variant, so it should be reachable from either. It is
/// reachable from exactly one, which makes reachability a function of how the
/// variant was *written* — the spelling-dependence root-caused on #1419 to a
/// weight bound whose threshold is the input's own spelling.
///
/// While this test passes, the family is unfixed. It is written to fail the day
/// any pair reaches its canonical form from both spellings, so that progress is
/// as loud as regression.
#[test]
fn each_pair_reaches_its_canonical_form_from_exactly_one_spelling() {
    for (a, b) in reported_pairs() {
        assert_ne!(
            a.verdict, b.verdict,
            "{} and {} have the same verdict; a pair is one canonical spelling \
             and one gap",
            a.label, b.label
        );

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
        assert_ne!(
            normalized(gap.input),
            canonical.output,
            "{} now reaches `{}` too, so the pair has converged on the form \
             its issue argues for. That is the fix: flip this row's verdict to \
             Canonical, lower OPEN_GAPS, and say in the PR which stored \
             spelling moved.",
            gap.label,
            canonical.output
        );
    }
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
/// **This is a 3'-direction claim and does not generalize.** Two of the nine gap
/// rows do move under 5' ([`FIVE_PRIME_MOVERS`]), so under that option the
/// argument above holds for seven of nine, not nine of nine — which is exactly
/// why the 5' answers are now pinned per row rather than left to a convergence
/// count.
#[test]
fn every_gap_row_is_returned_exactly_as_authored() {
    for row in REPORTED_ROWS {
        if row.verdict != Verdict::Gap {
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
/// `SequenceMatch` exactly, not merely `is_equivalent()`: the two spellings
/// normalize to different strings here, so anything stronger would mean the
/// pair had converged and the rows above are stale.
#[test]
fn every_reported_pair_is_still_one_variant_by_equivalence() {
    let checker = EquivalenceChecker::new(provider(TEMPLATE));
    for (a, b) in reported_pairs() {
        let left = parse_hgvs(a.input).unwrap_or_else(|e| panic!("{}: {e}", a.label));
        let right = parse_hgvs(b.input).unwrap_or_else(|e| panic!("{}: {e}", b.label));
        let result = checker
            .check(&left, &right)
            .unwrap_or_else(|e| panic!("{} vs {}: {e}", a.label, b.label));
        assert_eq!(
            result.level,
            EquivalenceLevel::SequenceMatch,
            "{} vs {}: expected SequenceMatch, got {:?}. If this is now \
             NormalizedMatch the pair has converged and the pinned rows are \
             stale; anything else is a regression in the fallback the reports \
             direct consumers to.",
            a.label,
            b.label,
            result.level
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
