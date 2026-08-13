//! No emitted pure insertion may be a tandem duplication (`DNA/insertion.md:5`).
//!
//! # The question, the ruling, the authority
//!
//! **Question.** `DNA/duplication.md:18` — "when a variant can be described as a
//! duplication, it **must** be described as a duplication and not as, e.g., an
//! insertion" — is one of the spec's genuine MUSTs. Does it bind on each *member*
//! of a cis allele, so that a partition exposing a describable duplication is
//! required; or does it bind on each *change* ferro derives, ranking the label
//! only?
//!
//! **Ruling.** The label, not the partition. The ledger record is
//! `duplication-must-ranks-the-label-not-the-partition`; read it there, not here,
//! because a rule written in two places drifts.
//!
//! **Authority**, both grounds:
//!
//! - `docs/recommendations/DNA/duplication.md:17` — "by definition, duplication
//!   may only be used when the additional copy is **directly 3'-flanking** of the
//!   original copy (a "tandem duplication")". `:18` is its sub-bullet, so the
//!   MUST ranks a label inside a scope `:17` fixes.
//! - `docs/background/glossary.md:310-311` — a variant is "a difference between a
//!   reference sequence and an observed sequence". `:18`'s grammatical subject is
//!   *a variant*, so it tests a difference, not a member of a chosen spelling.
//!
//! # Why this file exists at all
//!
//! Because the ruling has a falsifier, and until this file the falsifier had only
//! ever been measured *indirectly*.
//!
//! The reading above is safe only if a piece ferro spells `ins` genuinely is not a
//! duplication. `DNA/insertion.md:5` is definitional, not preferential —
//! "Insertion: a sequence change where, compared to the reference sequence, one or
//! more nucleotides are inserted **and** where the insertion is not a copy of a
//! sequence immediately 5'." A change that *is* such a copy is by definition not an
//! insertion, so an `ins`-labelled output of that shape would be a **rank-1
//! conformance defect**, not a merely non-preferred spelling — and the ruling
//! would be wrong as stated.
//!
//! The argument that it does not happen was that `anchor_for_piece` tries
//! `duplication_anchor` before every other typing a pure insertion can take, so a
//! piece a `dup` describes never reaches the insertion arm. That is an argument
//! about call order. This file asserts the property instead:
//!
//! > no emitted pure `ins` member has `payload == reference[junction - payload.len() .. junction]`
//!
//! which is the same predicate `is_tandem_duplication` computes, evaluated on the
//! rendered output rather than on the code path that produced it.
//!
//! # What the sabotage measured, and what it refuted
//!
//! A guard nobody has watched fail is not a guard, so this one was built against
//! a deliberate break. The break that works is **not** the one the argument above
//! names, and that is worth recording because it corrects where `:18` is actually
//! enforced:
//!
//! | forced off | this file |
//! |---|---|
//! | `merge::is_tandem_duplication` (so `duplication_anchor` never fires) | **still passes** |
//! | `rules::insertion_is_duplication` | **still passes** |
//! | `rules::insertion_to_duplication` | **still passes** |
//! | `five_prime_match` in `normalize::mod`'s insertion→duplication resolution | **FAILS** |
//!
//! So on these shapes the operative `:18` enforcement is the **per-member**
//! promotion in `normalize/mod.rs` — the inline `ref_seq[pos - len..pos] ==
//! rotated_seq` test in candidate arm (b) — and not the sequence-first
//! derivation's `duplication_anchor`. Three sabotages that left every output
//! byte-identical are the evidence; `TEMPLATE:g.[74_75insC;75_76insG;90A>T]`
//! still came back as `g.[73_74dup;90A>T]` under each. With `five_prime_match`
//! forced false the suite reports
//! `TEMPLATE:g.[264_265insATA;274del]` carrying payload `ATA` against a 5' flank
//! of `ATA` — the rank-1 shape, named.
//!
//! Do not read the table as a claim that `duplication_anchor` is dead code. It
//! reads as "no case in *this* corpus reaches it", which is a claim about the
//! corpus; what the table does establish is that citing call order in
//! `anchor_for_piece` is not by itself an argument that the property holds.
//!
//! # What the corpus is built to contain
//!
//! A zero here would otherwise be a claim about the corpus rather than about
//! ferro, so the corpus is constructed to make the violation *reachable*: every
//! sweep case authors an insertion whose payload **is** the run of reference bases
//! immediately 5' of its junction. Those are precisely the inputs that must come
//! back labelled `dup`, and each is offered both alone and beside a distant `del`
//! — the third member whose distance decides whether the derivation coalesces the
//! block, which is the geometry the whole D1 question turns on. Scrambled payloads
//! ride along as the negative half: they legitimately stay `ins`, and they are
//! what keeps the examined count non-zero on every arm.
//!
//! [`the_corpus_reaches_both_labels`] pins that construction, so a corpus that
//! stopped producing `dup`s (or stopped producing `ins`es) fails rather than
//! reporting a vacuous pass.
//!
//! # It guards whichever partition arm the process is running
//!
//! Nothing here reads or sets `FERRO_PARTITION`. The invariant is a property of
//! the emitted description, so it holds — or does not — on every arm, and running
//! the suite under `FERRO_PARTITION=canonical-coalesced` re-evaluates it against
//! the coalescing partitioner without any change to this file. That matters
//! because the coalesced arm is where the D1 rows actually lose their `dup`s.

use ferro_hgvs::hgvs::edit::{InsertedSequence, NaEdit};
use ferro_hgvs::hgvs::uncertainty::Mu;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::parse_hgvs;

use crate::common::cis_apply_oracle::{normalize, sweep_sequences};

/// Sequences drawn from the shared cis corpus generator.
///
/// Fixed, and deliberately **not** routed through the `FERRO_SWEEP_SEEDS` knob.
/// That knob exists to shrink the three exhaustive sweeps, which dominate a
/// local run at ~80s each; this file's whole corpus normalizes in well under a
/// second, so there is nothing to shrink. Routing through it would also oblige
/// `SWEEP_FILTER` in `ci.yml` to name this module —
/// `sweep_filter_invariant::every_seed_knob_consumer_is_named_in_the_sweep_filter`
/// enforces exactly that — which would move the file into the `sweeps` job and
/// make an invariant that should hold everywhere run at a reduced corpus in
/// every other job.
const CORPUS_SEEDS: u32 = 8;

/// `("ACGT" x 64) + core + ("ACGT" x 64)`, so `core[i]` is at 1-based position
/// `257 + i`. The same padding convention as `issue_1261_cis_member_order`, and
/// for the same reason: it keeps every shift well clear of a contig boundary,
/// where the boundary-`delins` rewrites would answer a different question.
fn padded(core: &str) -> String {
    format!("{}{}{}", "ACGT".repeat(64), core, "ACGT".repeat(64))
}

/// A pure-insertion member of a rendered description.
#[derive(Debug)]
struct PureInsertion {
    /// 1-based position of the reference base immediately 5' of the insertion
    /// point — the left half of the flanking pair `DNA/insertion.md:15` requires
    /// ("`positions flanking` should contain **two flanking nucleotides**, e.g.,
    /// `123_124`, not `123_125`").
    junction: u64,
    /// The literal inserted bases, as ASCII.
    payload: Vec<u8>,
}

/// Collect one member's pure insertion, if it is one.
///
/// Deliberately narrow. Only a genomic member, a `Certain` edit, a `Literal`
/// inserted sequence and a plain flanking position pair are decidable: a `dup`
/// renders as `NaEdit::Duplication` and never reaches here, and an uncertain or
/// non-flanking anchor is not a shape the predicate can be evaluated against.
/// Skipping those is safe in the direction that matters — the check can only
/// examine fewer members, never different ones — and
/// [`the_corpus_reaches_both_labels`] is what stops "fewer" from becoming
/// "none".
fn push_pure_insertion(variant: &HgvsVariant, into: &mut Vec<PureInsertion>) {
    let HgvsVariant::Genome(genome) = variant else {
        return;
    };
    let Some(NaEdit::Insertion {
        sequence: InsertedSequence::Literal(sequence),
    }) = genome.loc_edit.edit.inner()
    else {
        return;
    };
    let (Some(Mu::Certain(start)), Some(Mu::Certain(end))) = (
        genome.loc_edit.location.start.as_single(),
        genome.loc_edit.location.end.as_single(),
    ) else {
        return;
    };
    if start.special.is_some() || end.special.is_some() {
        return;
    }
    if end.base.checked_sub(start.base) != Some(1) {
        return;
    }
    into.push(PureInsertion {
        junction: start.base,
        payload: sequence.bases().iter().map(|b| b.to_char() as u8).collect(),
    });
}

/// Whether any member of `description` is rendered as a `dup`.
///
/// The mirror image of [`pure_insertions`], and parsed for the same reason: the
/// question is which **label** ferro applied, which is a property of the AST.
fn renders_a_duplication(description: &str) -> bool {
    fn is_duplication(variant: &HgvsVariant) -> bool {
        let edit = match variant {
            HgvsVariant::Genome(g) => g.loc_edit.edit.inner(),
            _ => return false,
        };
        matches!(edit, Some(NaEdit::Duplication { .. }))
    }

    let variant = parse_hgvs(description).unwrap_or_else(|error| {
        panic!("ferro emitted `{description}`, which it cannot read back: {error}")
    });
    match &variant {
        HgvsVariant::Allele(allele) => allele.variants.iter().any(is_duplication),
        lone => is_duplication(lone),
    }
}

/// Every pure-insertion member of `description`.
fn pure_insertions(description: &str) -> Vec<PureInsertion> {
    let variant = parse_hgvs(description).unwrap_or_else(|error| {
        panic!("ferro emitted `{description}`, which it cannot read back: {error}")
    });
    let mut found = Vec::new();
    match &variant {
        HgvsVariant::Allele(allele) => {
            for member in &allele.variants {
                push_pure_insertion(member, &mut found);
            }
        }
        lone => push_pure_insertion(lone, &mut found),
    }
    found
}

/// The invariant, on one normalized output. Returns how many pure insertions it
/// examined, so the caller can refuse a vacuous sweep.
fn assert_no_emitted_insertion_is_a_tandem_copy(seq: &str, input: &str, output: &str) -> usize {
    let bases = seq.as_bytes();
    let insertions = pure_insertions(output);
    for insertion in &insertions {
        let junction = insertion.junction as usize;
        // Nothing lies 5' of base 1, so a payload longer than the prefix cannot
        // be a copy of it. `checked_sub` is the same guard
        // `is_tandem_duplication` uses for the same reason.
        let Some(source_start) = junction.checked_sub(insertion.payload.len()) else {
            continue;
        };
        let flank = &bases[source_start..junction];
        assert!(
            !flank.eq_ignore_ascii_case(&insertion.payload),
            "`{input}` normalized to `{output}`, whose insertion at junction \
             {junction}|{} carries payload `{}` — a copy of the reference bases \
             `{}` immediately 5' of it. `DNA/insertion.md:5` defines an insertion \
             as one that is NOT such a copy, and `DNA/duplication.md:18` requires \
             a change that is one to be labelled `dup`. This is a rank-1 \
             conformance defect, and it falsifies the ruling record \
             `duplication-must-ranks-the-label-not-the-partition`, which reads \
             `:18` as ranking the label on the strength of this never happening.",
            junction + 1,
            String::from_utf8_lossy(&insertion.payload),
            String::from_utf8_lossy(flank),
        );
    }
    insertions.len()
}

/// The corpus: `(reference, input)` pairs, each authored to put a tandem-copy
/// payload within the derivation's reach.
///
/// For every sweep 20-mer, at every junction across the core and for payload
/// lengths 1..=3, two inputs carry the tandem payload (alone, and beside a `del`
/// 10 nt clear) and two carry a scrambled one. The `del` is the third member
/// whose presence or absence decides whether the block is coalesced — R1–R4 of
/// the D1 row set are exactly that geometry.
fn corpus() -> Vec<(String, String)> {
    let mut cases = Vec::new();
    for core in sweep_sequences(CORPUS_SEEDS) {
        let seq = padded(&core);
        let bases = seq.as_bytes();
        // The core occupies 1-based 257..=276. Stay inside it, leaving room for
        // the 3-base 5' flank and for the third member 10 nt to the 3' side.
        for junction in 260..=270usize {
            for length in 1..=3usize {
                let tandem = String::from_utf8(bases[junction - length..junction].to_vec())
                    .expect("padded reference is ASCII");
                // A payload that is NOT a copy of the 5' flank, so the emitted
                // member legitimately stays `ins`. Rotating each base one step
                // around the alphabet cannot reproduce the flank it came from.
                let scrambled: String = tandem
                    .chars()
                    .map(|base| match base {
                        'A' => 'C',
                        'C' => 'G',
                        'G' => 'T',
                        _ => 'A',
                    })
                    .collect();
                for payload in [&tandem, &scrambled] {
                    cases.push((
                        seq.clone(),
                        format!("TEMPLATE:g.{junction}_{}ins{payload}", junction + 1),
                    ));
                    cases.push((
                        seq.clone(),
                        format!(
                            "TEMPLATE:g.[{junction}_{}ins{payload};{}del]",
                            junction + 1,
                            junction + 10
                        ),
                    ));
                }
            }
        }
    }

    // The two D1 loci that are reproducible on `TEMPLATE`, added verbatim.
    //
    // R1/R2 — `issue_1261_cis_member_order` and its `issue_1320` twin. Authored
    // entirely as insertions; the `dup`s in the pinned value are ferro's own
    // per-member re-spelling, and the `268del` is what blocks the derivation.
    let r1 = padded("CAGTATGCAGGCAA");
    for input in [
        "TEMPLATE:g.[258_259insA;259_260insAG;268del]",
        "TEMPLATE:g.[258_259insA;259_260insAG]",
    ] {
        cases.push((r1.clone(), input.to_string()));
    }
    // R5/R6 — `cis_junction_crossing_shift`'s `DUP_RUN`, the locus
    // `contiguous-insertion-split-by-a-blocked-derivation` was decided on. Here
    // the input does carry a `dup`.
    let dup_run = "TTAATATATAATAATATTAT".to_string();
    for input in [
        "TEMPLATE:g.[4_5insC;5_6dup;15del]",
        "TEMPLATE:g.[3_4insACT;15del]",
        "TEMPLATE:g.[4_5insCTA;15del]",
    ] {
        cases.push((dup_run.clone(), input.to_string()));
    }
    cases
}

/// **The invariant.** Over the whole corpus, no emitted pure `ins` member repeats
/// the reference bases immediately 5' of its own insertion point.
///
/// A failure here is not a test to re-bless. It is a `DNA/insertion.md:5`
/// violation and it falsifies the ruling; the panic message says so.
#[test]
fn no_emitted_insertion_repeats_its_five_prime_flank() {
    let mut examined = 0usize;
    for (seq, input) in corpus() {
        let output = normalize(&seq, &input);
        examined += assert_no_emitted_insertion_is_a_tandem_copy(&seq, &input, &output);
    }
    // A corpus zero is a claim about the corpus. The sweep authors scrambled
    // payloads specifically so that pure insertions survive to be checked; if
    // none reach the assertion, the assertion measured nothing.
    assert!(
        examined > 0,
        "the corpus produced no pure insertion at all, so the invariant was never \
         evaluated — fix the corpus rather than reading this as a pass"
    );
}

/// The corpus is not vacuous **in both directions**: it reaches inputs ferro
/// labels `dup`, and inputs it labels `ins`.
///
/// The first half is what makes the invariant above meaningful. Its assertion can
/// only fail on a member that came back as `ins`, so a corpus in which nothing was
/// ever a candidate duplication passes it for free. The tandem-payload cases are
/// candidates by construction — and this row proves the construction survived
/// normalization rather than being shifted or merged into something else.
#[test]
fn the_corpus_reaches_both_labels() {
    let mut with_dup = 0usize;
    let mut with_ins = 0usize;
    for (seq, input) in corpus() {
        let output = normalize(&seq, &input);
        // Parsed, not `output.contains("dup")`. A substring scan over a
        // rendered description answers a different question than the one this
        // row asks — it would count an accession or an inserted-range payload
        // that happens to spell those three letters, and it is the sibling of
        // the `ins` half two lines down, which *is* parsed. Reading one half of
        // one claim off the text and the other off the AST is how the two come
        // to disagree without either being obviously wrong.
        if renders_a_duplication(&output) {
            with_dup += 1;
        }
        if !pure_insertions(&output).is_empty() {
            with_ins += 1;
        }
    }
    assert!(
        with_dup > 0,
        "no corpus case normalized to anything carrying a `dup`, so the tandem \
         payloads never became duplication candidates and \
         `no_emitted_insertion_repeats_its_five_prime_flank` is vacuous"
    );
    assert!(
        with_ins > 0,
        "no corpus case normalized to a pure insertion, so the invariant has \
         nothing to evaluate"
    );
}
