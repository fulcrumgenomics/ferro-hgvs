//! The authored-inversion sweep: schema, generator parameters, and outcome
//! classification.
//!
//! # What the sweep is
//!
//! 2,075 authored `inv` descriptions on one real transcript — `NM_004006.2`,
//! the DMD coding sequence the HGVS recommendations use for their own worked
//! inversion examples — enumerated as
//! `NM_004006.2:c.{start}_{start + len - 1}inv` for every `start` in
//! `101, 108, …, 2999` and every `len` in `{4, 6, 8, 12, 16}`.
//!
//! The parameters live here rather than in the test so that the generator that
//! captures the reference slice and the gate that replays it cannot enumerate
//! different corpora — the same single-source-of-truth arrangement as
//! [`super::case_harvest`] and [`super::spec_worked_examples`].
//!
//! # What is pinned, and which kind of record each pin is
//!
//! Four layers, and they are **not** redundant — each catches something the
//! others cannot:
//!
//! 1. **Exact outputs for all 2,075 rows** (`cases.tsv`). Characterization
//!    pins. They are the only layer that notices the 3'-shift arithmetic
//!    moving *within* the set of correct answers: `c.101_108inv ->
//!    c.102_107inv` silently becoming `-> c.103_106inv` keeps the same outcome
//!    bucket and the same census total — so the other layers can stay green
//!    while a string a downstream consumer stores has moved.
//! 2. **The sequence oracle**, [`member_edits`] + [`apply_member_edits`]. The
//!    output's own members, applied to the committed bases, must produce the
//!    same sequence as replacing the authored span with its reverse complement.
//!    Neither side calls the normalizer, so this is not a self-consistency
//!    check — and it is indifferent to the *shape* of the answer, which is what
//!    lets it survive a repartitioning.
//! 3. **The shape invariant**, [`SweepOutcome`]. An authored `inv` normalizes
//!    to an `inv` (possibly 3'-shifted or narrowed), to a no-op `=` when the
//!    span equals its own reverse complement, or to allele members that are
//!    each on this transcript and each length-preserving. It never silently
//!    becomes a lone substitution, deletion or duplication, and no member may
//!    change the sequence's length. This is what makes a careless re-blessing
//!    of layer 1 still fail.
//! 4. **The census** — the four bucket counts, as named constants in the gate.
//!    A bulk re-bless that moves the *distribution* is visible in a diff and has
//!    to be argued for.
//!
//! Layer 3 was narrower once, and the narrowing was wrong. It read "never split
//! into allele members", which was true of the tree this corpus was authored
//! against and stopped being true when #1484 widened the sequence-derived axis
//! gate to `c.`/`n.`/`r.`. See [`SweepOutcome::Repartitioned`].
//!
//! Only three rows are adjudications rather than characterization: the spec
//! publishes `NM_004006.2:c.5657_5660inv` and `NM_004006.2:c.4145_4160inv`
//! outright (`DNA/inversion.md:27-34`), and the `delins` spelling of the first
//! converges on it. The gate states which is which; a reader should not have to
//! guess.
//!
//! # The authority for the invariant
//!
//! > Inversion: a sequence change where, compared to a reference sequence,
//! > **more than one nucleotide** replacing the original sequence is the
//! > reverse complement of the original sequence.
//! >   — `assets/hgvs-nomenclature/docs/recommendations/DNA/inversion.md:5`
//!
//! An authored `inv` asserts exactly that relation between a reference span and
//! its replacement. Normalization may move where the span sits (the 3'rule,
//! `DNA/inversion.md:17`), may retype it when the replacement turns out to
//! equal the reference (a palindromic span is its own reverse complement, so
//! nothing changed), and may re-derive the block as several members — but
//! whatever it emits must still *denote* that span replaced by its reverse
//! complement. That is the sequence oracle, and it is the reading of `:5` that
//! does not have to be revised every time the partitioner changes.

use std::path::Path;

use crate::error_handling::ErrorConfig;
use crate::hgvs::edit::{InsertedSequence, NaEdit};
use crate::hgvs::interval::UncertainBoundary;
use crate::hgvs::location::CdsPos;
use crate::hgvs::variant::Accession;
use crate::{HgvsVariant, NormalizeConfig, ShuffleDirection};

/// Repo-relative path to the pinned outputs.
pub const CASES_PATH: &str = "tests/fixtures/inversion-sweep/cases.tsv";
/// Repo-relative path to the committed hermetic reference slice.
pub const WINDOWS_PATH: &str = "tests/fixtures/inversion-sweep/reference-windows.json";
/// Repo-relative path to the vendored HGVS spec checkout citations resolve
/// against.
pub const SPEC_DIR: &str = "assets/hgvs-nomenclature";

/// The transcript the sweep is authored on: the coding DNA reference sequence
/// `DNA/inversion.md:30-34` uses for its own worked inversion examples.
pub const TRANSCRIPT: &str = "NM_004006.2";

/// First `c.` start position of the sweep.
pub const START_FIRST: u64 = 101;
/// Exclusive upper bound on the `c.` start position.
pub const START_LIMIT: u64 = 3000;
/// Stride between consecutive start positions. Deliberately coprime with 2, 3
/// and 4 so the sweep does not sit on one codon phase.
pub const START_STEP: u64 = 7;
/// Inverted span lengths, in nucleotides. All even, all `> 1`: an odd-length
/// span has a middle base that is its own reverse complement only when it is
/// absent, and `DNA/inversion.md:15-16` forbids a one-nucleotide inversion
/// outright.
pub const LENGTHS: [u64; 5] = [4, 6, 8, 12, 16];

/// The three rows that are **adjudications** rather than characterization
/// pins, and whose expected outputs therefore live in the gate's source next to
/// the verbatim clause that settles them (`DNA/inversion.md:27-34`).
///
/// They sit outside [`sweep_inputs`] — `c.5657` and `c.4145` are past the
/// sweep's `start < 3000` bound, and the `delins` row is not an authored `inv`
/// at all — but the reference slice must still serve their reads, so the
/// generator runs them through the same recording pass.
pub const ADJUDICATED_INPUTS: [&str; 3] = [
    "NM_004006.2:c.5657_5660inv",
    "NM_004006.2:c.4145_4160inv",
    "NM_004006.2:c.5657_5660delinsTCAG",
];

/// One authored inversion: a `c.` start and a span length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SweepInput {
    /// 1-based `c.` position of the first inverted nucleotide.
    pub start: u64,
    /// Number of nucleotides inverted.
    pub length: u64,
}

impl SweepInput {
    /// 1-based `c.` position of the last inverted nucleotide.
    pub fn end(&self) -> u64 {
        self.start + self.length - 1
    }

    /// The HGVS description handed to the normalizer.
    pub fn description(&self) -> String {
        format!("{TRANSCRIPT}:c.{}_{}inv", self.start, self.end())
    }
}

/// The whole sweep, in the order the fixture stores it: start-major, then
/// length in [`LENGTHS`] order.
///
/// Deterministic and total — the fixture's input column is asserted equal to
/// this, so the pins and the generator cannot drift onto different corpora.
pub fn sweep_inputs() -> Vec<SweepInput> {
    let mut inputs = Vec::new();
    let mut start = START_FIRST;
    while start < START_LIMIT {
        for length in LENGTHS {
            inputs.push(SweepInput { start, length });
        }
        start += START_STEP;
    }
    inputs
}

/// The normalize settings the generator records and the gate replays:
/// `ferro normalize --error-mode lenient` with 3'-shifting.
///
/// Shared for the same reason [`sweep_inputs`] is. Centralizing the *corpus*
/// while leaving the *configuration* duplicated would let one side change
/// shuffle direction or error mode while every pin still compared equal in
/// shape — the gate would replay a configuration the generator never recorded,
/// and the pins would be silently about something else.
pub fn sweep_normalize_config() -> (ErrorConfig, NormalizeConfig) {
    let error_config = ErrorConfig::lenient();
    let normalize_config =
        NormalizeConfig::for_entry_point(ShuffleDirection::ThreePrime, error_config.clone());
    (error_config, normalize_config)
}

/// What normalization did to one authored inversion.
///
/// The first four are the conformant outcomes and the census buckets; the last
/// three are the ways the invariant can break, kept apart so a failure names
/// *which* way rather than "not an inversion".
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SweepOutcome {
    /// Returned character-for-character unchanged: already the 3'-most
    /// spelling.
    Unchanged,
    /// Still a lone `inv`, but at a different span — 3'-shifted or re-typed
    /// onto a narrower block.
    ShiftedInversion,
    /// A no-op `=`. The authored span equals its own reverse complement, so
    /// inverting it changes nothing.
    PalindromicNoOp,
    /// Repartitioned into allele members, carrying the member count.
    ///
    /// **Conformant, and it did not used to be.** This corpus was authored
    /// against a tree where the `c.` axis did not run the sequence-derived
    /// partitioner, so every row stayed a lone edit and the invariant was
    /// written as "never split into members". #1484 widened that gate to
    /// `c.`/`n.`/`r.`, and 155 of these 2,075 rows repartitioned; #1706's
    /// post-hoc run scan then absorbed 62 of them, leaving **93** — the count
    /// pinned as `CENSUS_REPARTITIONED` in `tests/it/inversion_sweep.rs`, which
    /// is where to read it rather than from this prose. That split was once
    /// read as the *decided* direction for the substitution case, on
    /// `rulings[inversion-vs-two-delins-76-83]`'s note that "#1230's
    /// substitution case is untouched and still splits". **That reading is
    /// superseded**: `rulings[whole-span-reverse-complement-types-as-inv]`
    /// (2026-08-13) types a whole-span reverse complement `inv` uniformly, so
    /// all 93 — every one of them the all-substitutions case after #1706 — are
    /// decided *against* and are expected to stop repartitioning once
    /// #1703/#1541/#1575 implement it.
    ///
    /// So a split is no longer evidence of a defect by itself, and the weight
    /// moves onto the two checks that *are* still falsifiable: no member may
    /// change the sequence's length (below), and the members must jointly
    /// reproduce the inverted span
    /// ([`apply_member_edits`], asserted per row in the gate).
    Repartitioned(usize),
    /// A member changes the sequence's length. Inverting a span replaces it
    /// with a block of the same length, so a `del`, `dup`, `ins`, `repeat` or
    /// unbalanced `delins` member asserts something `DNA/inversion.md:5` does
    /// not. Carries the offending member.
    MemberChangesLength(String),
    /// Retyped to some other lone edit — a substitution, deletion, duplication,
    /// delins. Carries the rendered edit kind.
    Retyped(String),
    /// The output is not a `c.` description on this transcript at all.
    Unexpected(String),
}

impl SweepOutcome {
    /// Whether this outcome satisfies the invariant.
    pub fn is_conformant(&self) -> bool {
        matches!(
            self,
            SweepOutcome::Unchanged
                | SweepOutcome::ShiftedInversion
                | SweepOutcome::PalindromicNoOp
                | SweepOutcome::Repartitioned(_)
        )
    }

    /// Short label for census reporting.
    pub fn label(&self) -> &'static str {
        match self {
            SweepOutcome::Unchanged => "unchanged",
            SweepOutcome::ShiftedInversion => "shifted-still-inv",
            SweepOutcome::PalindromicNoOp => "palindromic-no-op",
            SweepOutcome::Repartitioned(_) => "repartitioned",
            SweepOutcome::MemberChangesLength(_) => "MEMBER-CHANGES-LENGTH",
            SweepOutcome::Retyped(_) => "RETYPED",
            SweepOutcome::Unexpected(_) => "UNEXPECTED",
        }
    }
}

/// Whether `accession` names [`TRANSCRIPT`], **with or without** the genomic
/// wrapper an intronic position requires.
///
/// The wrapper is not a re-anchor, and telling the two apart is the whole point
/// of this helper (#1704). `checklist.md:20` makes `NC_000023.11(NM_004006.2)`
/// mandatory once a description names an intronic position, and four rows of this
/// sweep shift into an intron — so ferro renders them compound, and a bare string
/// comparison read that as "re-anchored onto another transcript", the exact
/// retype-in-disguise the check exists to catch. It is not one: the inner
/// accession is unchanged, and `Accession::without_genomic_context` is what says
/// so structurally rather than by string surgery.
///
/// What is deliberately NOT checked here is that the wrapper only appears where
/// an intronic position needs it. A gratuitous wrapper would be a representation
/// change like any other, and `every_authored_inversion_produces_its_pinned_output`
/// pins all 2,075 strings exactly — so it would fail there, with the string in
/// hand, rather than being silently absorbed by a looser predicate.
fn names_the_sweep_transcript(accession: &Accession) -> bool {
    accession.clone().without_genomic_context().full() == TRANSCRIPT
}

/// Classify a normalized result against the authored input it came from.
///
/// `rendered` is the output string and `input` the description handed in;
/// [`SweepOutcome::Unchanged`] is decided by string equality with `input`, and
/// everything else by the output's shape.
pub fn classify(input: &str, normalized: &HgvsVariant, rendered: &str) -> SweepOutcome {
    if let HgvsVariant::Allele(allele) = normalized {
        for member in &allele.variants {
            if let Some(outcome) = member_violation(member, rendered) {
                return outcome;
            }
        }
        return SweepOutcome::Repartitioned(allele.variants.len());
    }
    let HgvsVariant::Cds(cds) = normalized else {
        return SweepOutcome::Unexpected(format!(
            "output is not a c. description on {TRANSCRIPT}: {rendered}"
        ));
    };
    // The accession is checked, not assumed: without this an output re-anchored
    // onto another transcript — or onto another *version* of this one — would be
    // classified `Unchanged` or `ShiftedInversion` and counted as conformant,
    // which is precisely the retype-in-disguise the invariant exists to catch.
    if !names_the_sweep_transcript(&cds.accession) {
        return SweepOutcome::Unexpected(format!(
            "output is a c. description on {} rather than on {TRANSCRIPT}: {rendered}",
            cds.accession.full()
        ));
    }
    let Some(edit) = cds.loc_edit.edit.inner() else {
        return SweepOutcome::Unexpected(format!("output carries no edit: {rendered}"));
    };
    match edit {
        NaEdit::Inversion { .. } => {
            if rendered == input {
                SweepOutcome::Unchanged
            } else {
                SweepOutcome::ShiftedInversion
            }
        }
        NaEdit::Identity { .. } => SweepOutcome::PalindromicNoOp,
        other => SweepOutcome::Retyped(edit_kind(other).to_string()),
    }
}

/// Why `member` is not an admissible piece of a repartitioned inversion, if it
/// is not.
///
/// Two grounds, both from `DNA/inversion.md:5`. The member must sit on the same
/// coding sequence — a member re-anchored elsewhere is describing a different
/// reference — and it must replace bases with the same number of bases, since
/// the reverse complement of a span is exactly as long as the span.
fn member_violation(member: &HgvsVariant, rendered: &str) -> Option<SweepOutcome> {
    let HgvsVariant::Cds(cds) = member else {
        return Some(SweepOutcome::Unexpected(format!(
            "a member is not a c. description on {TRANSCRIPT}: {rendered}"
        )));
    };
    if !names_the_sweep_transcript(&cds.accession) {
        return Some(SweepOutcome::Unexpected(format!(
            "a member is on {} rather than on {TRANSCRIPT}: {rendered}",
            cds.accession.full()
        )));
    }
    let edit = cds.loc_edit.edit.inner()?;
    match edit {
        // A substitution replaces one base with one base, an inversion a block
        // with its own reverse complement, and `=` nothing at all.
        NaEdit::Substitution { .. } | NaEdit::SubstitutionNoRef { .. } => None,
        NaEdit::Inversion { .. } | NaEdit::Identity { .. } => None,
        // A delins may or may not preserve length, so it is measured rather
        // than assumed — when it can be. An endpoint carrying an intron offset
        // has no `c.`-axis length, so the balance is left to the sequence
        // oracle, which the gate runs over every row whose span lies in one
        // exon.
        NaEdit::Delins { sequence, .. } => {
            let span = plain_span(cds)?;
            let inserted = literal_insertion(sequence)?;
            (inserted.len() as u64 != span).then(|| {
                SweepOutcome::MemberChangesLength(format!(
                    "a delins member replaces {span} base(s) with {}: {rendered}",
                    inserted.len()
                ))
            })
        }
        other => Some(SweepOutcome::MemberChangesLength(format!(
            "a member is a {}, which cannot preserve the span's length: {rendered}",
            edit_kind(other)
        ))),
    }
}

/// The 1-based inclusive `c.` position a boundary names, or `None` when it is
/// not a plain coding position — an intron offset, a `*`/`-` UTR position, a
/// `pter`-style marker, an uncertainty range or a `?`.
///
/// `None` means "this cannot be read on the coding axis", never "zero". Every
/// caller treats it as a refusal to judge rather than as a value.
fn plain_position(boundary: &UncertainBoundary<CdsPos>) -> Option<i64> {
    let UncertainBoundary::Single(mu) = boundary else {
        return None;
    };
    let pos = mu.inner()?;
    if pos.offset.is_some() || pos.utr3 || pos.special.is_some() || pos.base < 1 {
        return None;
    }
    Some(pos.base)
}

/// The number of coding bases a `c.` interval spans, when both endpoints are
/// plain coding positions.
fn plain_span(cds: &crate::hgvs::variant::CdsVariant) -> Option<u64> {
    let start = plain_position(&cds.loc_edit.location.start)?;
    let end = plain_position(&cds.loc_edit.location.end)?;
    u64::try_from(end - start + 1).ok()
}

/// The literal bases an insertion carries, when it is literal at all — a
/// `ins10`, `insN[15]` or multi-part payload has no single spelling to compare.
fn literal_insertion(sequence: &InsertedSequence) -> Option<String> {
    match sequence {
        InsertedSequence::Literal(bases) => Some(bases.to_string()),
        _ => None,
    }
}

/// One member of a normalized output, reduced to the coding-axis block it
/// replaces and what it replaces it with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MemberEdit {
    /// 1-based inclusive `c.` position of the first replaced base.
    pub start: u64,
    /// 1-based inclusive `c.` position of the last replaced base.
    pub end: u64,
    /// The bases that replace `[start, end]`.
    pub replacement: MemberReplacement,
}

/// What a [`MemberEdit`] puts in place of the block it names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MemberReplacement {
    /// Literal bases, from a substitution or a `delins`.
    Literal(String),
    /// The reverse complement of the block itself, from an `inv`.
    ReverseComplement,
    /// Nothing changes, from a `=`.
    Unchanged,
}

/// Reduce a normalized output to its members' coding-axis edits, or say why it
/// cannot be read on that axis.
///
/// This is the front half of the sequence oracle: it converts a description
/// into edits that [`apply_member_edits`] can apply to the stored transcript
/// bases **without going through the normalizer**, which is what stops the
/// check being circular. It deliberately refuses rather than guesses — an
/// intron offset, a UTR position or a non-literal payload means the coding axis
/// is not where this member lives, and a caller that treated a refusal as "no
/// change" would turn an unreadable row into a passing one.
pub fn member_edits(normalized: &HgvsVariant) -> Result<Vec<MemberEdit>, String> {
    let members: Vec<&HgvsVariant> = match normalized {
        HgvsVariant::Allele(allele) => allele.variants.iter().collect(),
        single => vec![single],
    };
    let mut edits = Vec::with_capacity(members.len());
    for member in members {
        let HgvsVariant::Cds(cds) = member else {
            return Err(format!("member is not a c. description: {member}"));
        };
        let start = plain_position(&cds.loc_edit.location.start)
            .ok_or_else(|| format!("member start is not a plain coding position: {member}"))?;
        let end = plain_position(&cds.loc_edit.location.end)
            .ok_or_else(|| format!("member end is not a plain coding position: {member}"))?;
        let edit = cds
            .loc_edit
            .edit
            .inner()
            .ok_or_else(|| format!("member carries no edit: {member}"))?;
        let replacement = match edit {
            NaEdit::Substitution { alternative, .. } => {
                MemberReplacement::Literal(alternative.to_string())
            }
            NaEdit::SubstitutionNoRef { alternative } => {
                MemberReplacement::Literal(alternative.to_string())
            }
            NaEdit::Inversion { .. } => MemberReplacement::ReverseComplement,
            NaEdit::Identity { .. } => MemberReplacement::Unchanged,
            NaEdit::Delins { sequence, .. } => MemberReplacement::Literal(
                literal_insertion(sequence)
                    .ok_or_else(|| format!("member's insertion is not literal: {member}"))?,
            ),
            other => {
                return Err(format!(
                    "member is a {}, which the coding-axis oracle does not apply: {member}",
                    edit_kind(other)
                ))
            }
        };
        edits.push(MemberEdit {
            start: u64::try_from(start)
                .map_err(|_| format!("member start is negative: {member}"))?,
            end: u64::try_from(end).map_err(|_| format!("member end is negative: {member}"))?,
            replacement,
        });
    }
    Ok(edits)
}

/// Apply `edits` to `bases` and return the resulting sequence.
///
/// `bases` is the stored transcript sequence and `cds_start` its 1-based
/// transcript position of `c.1`, so a `c.` position `n` sits at
/// `bases[cds_start + n - 2]`. Applied 3'→5' so an earlier splice never moves a later one's
/// coordinates, and overlapping members are refused rather than silently
/// resolved by application order — two members claiming one base denote no
/// single sequence, which is a defect in the output rather than a detail of the
/// oracle.
pub fn apply_member_edits(
    bases: &str,
    cds_start: u64,
    edits: &[MemberEdit],
) -> Result<String, String> {
    let mut ordered: Vec<&MemberEdit> = edits.iter().collect();
    ordered.sort_by_key(|edit| std::cmp::Reverse(edit.start));
    let mut result = bases.as_bytes().to_vec();
    let mut claimed_from = u64::MAX;
    for edit in ordered {
        if edit.end < edit.start {
            return Err(format!("member spans c.{}_{}", edit.start, edit.end));
        }
        if edit.end >= claimed_from {
            return Err(format!(
                "member c.{}_{} overlaps the one 3' of it, so the output denotes no single \
                 sequence",
                edit.start, edit.end
            ));
        }
        claimed_from = edit.start;
        let lo = usize::try_from(cds_start + edit.start - 2)
            .map_err(|_| format!("c.{} is out of range", edit.start))?;
        let hi = usize::try_from(cds_start + edit.end - 1)
            .map_err(|_| format!("c.{} is out of range", edit.end))?;
        if hi > result.len() || lo >= hi {
            return Err(format!(
                "c.{}_{} maps to [{lo}, {hi}), outside the {} stored bases",
                edit.start,
                edit.end,
                result.len()
            ));
        }
        let block = String::from_utf8(result[lo..hi].to_vec())
            .map_err(|e| format!("stored bases are not ASCII: {e}"))?;
        let replacement = match &edit.replacement {
            MemberReplacement::Literal(bases) => bases.to_ascii_uppercase(),
            MemberReplacement::ReverseComplement => reverse_complement(&block),
            MemberReplacement::Unchanged => block,
        };
        result.splice(lo..hi, replacement.bytes());
    }
    String::from_utf8(result).map_err(|e| format!("edited bases are not ASCII: {e}"))
}

/// A stable, human-readable name for the edit an output was retyped to.
fn edit_kind(edit: &NaEdit) -> &'static str {
    match edit {
        NaEdit::Substitution { .. } | NaEdit::SubstitutionNoRef { .. } => "substitution",
        NaEdit::Deletion { .. } | NaEdit::NPaddedDeletion { .. } => "deletion",
        NaEdit::Insertion { .. } | NaEdit::BreakpointInsertion { .. } => "insertion",
        NaEdit::Delins { .. } => "delins",
        NaEdit::Duplication { .. } => "duplication",
        NaEdit::DupIns { .. } => "dupins",
        NaEdit::Repeat { .. } | NaEdit::MultiRepeat { .. } => "repeat",
        NaEdit::Conversion { .. } => "conversion",
        NaEdit::CopyNumber { .. } => "copy-number",
        _ => "other",
    }
}

/// The reverse complement of an ASCII nucleotide string, uppercased.
///
/// Used to decide whether a span is its own reverse complement, which is the
/// only thing that makes a [`SweepOutcome::PalindromicNoOp`] correct. It is the
/// *independent* oracle for that question — it answers to
/// `DNA/inversion.md:5` rather than to ferro's own output — so it complements
/// the whole IUPAC alphabet rather than only `ACGT`: passing `R` through
/// uncomplemented would make an ambiguity-carrying span look self-complementary
/// when it is not, and the disagreement would surface as a wrong verdict with
/// nothing pointing at the cause.
///
/// `N`, `S` and `W` are their own complements. A byte outside the IUPAC
/// alphabet — which a real transcript record does not carry — is passed through
/// unchanged and is therefore treated as self-complementary; that is stated
/// here rather than left to be discovered, since the only way one can appear is
/// a corrupted fixture.
pub fn reverse_complement(bases: &str) -> String {
    bases
        .bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'U' => 'A',
            b'R' => 'Y',
            b'Y' => 'R',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'V' => 'B',
            b'D' => 'H',
            b'H' => 'D',
            // `S` (G/C), `W` (A/T) and `N` are their own complements.
            other => other as char,
        })
        .collect()
}

/// One pinned row: the authored description and the exact string ferro
/// produces for it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PinnedCase {
    pub input: String,
    pub output: String,
}

/// The committed pin table, loaded from the two-column TSV.
#[derive(Debug, Clone, Default)]
pub struct PinnedCases {
    pub rows: Vec<PinnedCase>,
}

impl PinnedCases {
    /// Parse the TSV: `#` comments and blank lines ignored, then one
    /// `input<TAB>output` row each.
    ///
    /// The error is a plain `String` rather than a [`FerroError`](crate::FerroError)
    /// because every way this can fail is a hand-edit of a generated fixture,
    /// and the remedy is a regeneration command the caller names — not a
    /// condition any library code should be matching on.
    pub fn from_tsv(content: &str) -> Result<Self, String> {
        let mut rows = Vec::new();
        for (number, line) in content.lines().enumerate() {
            let line = line.trim_end();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let Some((input, output)) = line.split_once('\t') else {
                return Err(format!(
                    "{CASES_PATH}:{}: expected `input<TAB>output`, got {line:?}",
                    number + 1
                ));
            };
            if input.is_empty() || output.is_empty() {
                return Err(format!(
                    "{CASES_PATH}:{}: empty column in {line:?}",
                    number + 1
                ));
            }
            rows.push(PinnedCase {
                input: input.to_string(),
                output: output.to_string(),
            });
        }
        Ok(Self { rows })
    }

    /// Load from a TSV file.
    pub fn from_tsv_path(path: &Path) -> Result<Self, String> {
        let content =
            std::fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?;
        Self::from_tsv(&content)
    }

    /// Render the TSV, header comment included, with a trailing newline so a
    /// byte-comparing `--check` is stable.
    pub fn to_tsv(&self, header: &str) -> String {
        let mut out = String::new();
        for line in header.lines() {
            out.push_str("# ");
            out.push_str(line);
            out.push('\n');
        }
        for row in &self.rows {
            out.push_str(&row.input);
            out.push('\t');
            out.push_str(&row.output);
            out.push('\n');
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sweep_enumerates_every_start_length_pair_exactly_once() {
        let inputs = sweep_inputs();
        // 415 starts (101..3000 step 7) x 5 lengths. The last start is
        // `101 + 7 * 414 = 2999`, the largest one under `START_LIMIT`.
        assert_eq!(inputs.len(), 2075);
        let mut seen: Vec<(u64, u64)> = inputs.iter().map(|i| (i.start, i.length)).collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), inputs.len(), "the sweep repeats a (start, len)");
        assert_eq!(
            inputs[0],
            SweepInput {
                start: 101,
                length: 4
            }
        );
        assert_eq!(
            *inputs.last().unwrap(),
            SweepInput {
                start: 2999,
                length: 16
            }
        );
        assert!(
            inputs.iter().all(|input| input.start < START_LIMIT),
            "a start escaped the sweep's exclusive upper bound"
        );
    }

    #[test]
    fn a_description_spans_exactly_its_length() {
        let input = SweepInput {
            start: 101,
            length: 8,
        };
        assert_eq!(input.end(), 108);
        assert_eq!(input.description(), "NM_004006.2:c.101_108inv");
    }

    #[test]
    fn every_sweep_length_is_a_legal_inversion_width() {
        // `DNA/inversion.md:15-16`: "the region inverted (`positions_inverted`)
        // contains **more than one nucleotide**".
        assert!(LENGTHS.iter().all(|&len| len > 1));
    }

    #[test]
    fn a_palindromic_span_is_its_own_reverse_complement() {
        assert_eq!(reverse_complement("CTGA"), "TCAG");
        assert_eq!(reverse_complement("GGCC"), "GGCC");
        assert_eq!(reverse_complement("acgt"), "ACGT");
    }

    #[test]
    fn the_reverse_complement_covers_the_iupac_alphabet() {
        // `R` (A/G) complements to `Y` (C/T) and back; passing either through
        // uncomplemented — as an ACGT-only table does — would make `RY` look
        // like `YR` and an ambiguity-carrying span be judged on wrong bases.
        assert_eq!(reverse_complement("R"), "Y");
        assert_eq!(reverse_complement("Y"), "R");
        assert_eq!(reverse_complement("RY"), "RY");
        assert_eq!(reverse_complement("KMBVDH"), "DHBVKM");
        // Self-complementary: `N` (any), `S` (G/C), `W` (A/T).
        assert_eq!(reverse_complement("N"), "N");
        assert_eq!(reverse_complement("S"), "S");
        assert_eq!(reverse_complement("W"), "W");
    }

    #[test]
    fn an_output_on_another_accession_is_unexpected_rather_than_unchanged() {
        // The shape is a lone `inv` and the string is character-for-character
        // what was handed in — everything except the accession says
        // `Unchanged`. Only the accession check separates them.
        let input = "NM_004006.2:c.101_104inv";
        let ours = crate::parse_hgvs(input).expect("parse");
        assert_eq!(classify(input, &ours, input), SweepOutcome::Unchanged);

        let elsewhere = "NM_000088.3:c.101_104inv";
        let other = crate::parse_hgvs(elsewhere).expect("parse");
        assert!(matches!(
            classify(input, &other, elsewhere),
            SweepOutcome::Unexpected(_)
        ));

        // A different *version* of the same transcript is a different reference
        // sequence, so it is rejected on the same grounds.
        let other_version = "NM_004006.3:c.101_104inv";
        let versioned = crate::parse_hgvs(other_version).expect("parse");
        assert!(matches!(
            classify(input, &versioned, other_version),
            SweepOutcome::Unexpected(_)
        ));
    }

    /// The whole point of the oracle: a repartition and the lone `inv` it
    /// replaced must denote the same bases.
    #[test]
    fn a_repartition_and_the_inversion_it_replaced_denote_one_sequence() {
        // c.2_5 is `TGGC`, whose reverse complement is `GCCA` — so inverting it
        // is the same four substitutions, and the oracle says so from the bases
        // rather than from ferro's classification of either spelling.
        let bases = "ATGGCCTTAG";
        let inversion = [MemberEdit {
            start: 2,
            end: 5,
            replacement: MemberReplacement::ReverseComplement,
        }];
        let split: Vec<MemberEdit> = [(2, "G"), (3, "C"), (4, "C"), (5, "A")]
            .into_iter()
            .map(|(at, base)| MemberEdit {
                start: at,
                end: at,
                replacement: MemberReplacement::Literal(base.to_string()),
            })
            .collect();

        let inverted = apply_member_edits(bases, 1, &inversion).unwrap();
        assert_eq!(inverted, "AGCCACTTAG");
        assert_eq!(apply_member_edits(bases, 1, &split).unwrap(), inverted);
    }

    #[test]
    fn the_oracle_reads_c_positions_through_the_cds_start() {
        // `c.1` is the fourth stored base, so `c.1_2` is `GG`, not `AT`.
        let bases = "ATGGCCTTAG";
        let edits = [MemberEdit {
            start: 1,
            end: 2,
            replacement: MemberReplacement::ReverseComplement,
        }];
        assert_eq!(apply_member_edits(bases, 3, &edits).unwrap(), "ATCCCCTTAG");
    }

    #[test]
    fn the_oracle_refuses_overlapping_members_rather_than_ordering_them() {
        // Two members claiming one base denote no single sequence. Applying
        // them in some order would invent one, and the invented answer would
        // then be compared against the inversion and could even match.
        let edits = [
            MemberEdit {
                start: 2,
                end: 4,
                replacement: MemberReplacement::Literal("AAA".to_string()),
            },
            MemberEdit {
                start: 4,
                end: 5,
                replacement: MemberReplacement::Literal("CC".to_string()),
            },
        ];
        let error = apply_member_edits("ATGGCCTTAG", 1, &edits).unwrap_err();
        assert!(error.contains("overlaps"), "{error}");
    }

    #[test]
    fn member_edits_are_read_from_a_repartitioned_output() {
        let parsed = crate::parse_hgvs("NM_004006.2:c.[101A>T;104A>T]").expect("parse");
        let edits = member_edits(&parsed).expect("readable on the coding axis");
        assert_eq!(
            edits,
            vec![
                MemberEdit {
                    start: 101,
                    end: 101,
                    replacement: MemberReplacement::Literal("T".to_string())
                },
                MemberEdit {
                    start: 104,
                    end: 104,
                    replacement: MemberReplacement::Literal("T".to_string())
                },
            ]
        );
    }

    #[test]
    fn member_edits_refuses_an_intronic_endpoint_instead_of_guessing() {
        // `c.1332-1` is in an intron, so the coding axis cannot address it and
        // the row has to be excluded from the sequence check rather than
        // silently read as `c.1332`.
        let parsed = crate::parse_hgvs("NM_004006.2:c.1322_1332-1inv").expect("parse");
        let error = member_edits(&parsed).unwrap_err();
        assert!(error.contains("not a plain coding position"), "{error}");
    }

    #[test]
    fn a_repartition_is_conformant_and_a_length_changing_member_is_not() {
        let input = "NM_004006.2:c.101_104inv";
        for split in [
            "NM_004006.2:c.[101A>T;104A>T]",
            "NM_004006.2:c.[101A>T;102_103inv;104A>T]",
        ] {
            let parsed = crate::parse_hgvs(split).expect("parse");
            assert_eq!(
                classify(input, &parsed, split),
                SweepOutcome::Repartitioned(if split.contains("inv") { 3 } else { 2 })
            );
        }

        // A deletion member removes a base, and the reverse complement of a
        // span is exactly as long as the span — so this cannot be a spelling of
        // the authored inversion whatever else is true of it.
        let deletion = "NM_004006.2:c.[101A>T;104del]";
        let parsed = crate::parse_hgvs(deletion).expect("parse");
        assert!(matches!(
            classify(input, &parsed, deletion),
            SweepOutcome::MemberChangesLength(_)
        ));

        // Same for an unbalanced delins: three bases out, two in.
        let unbalanced = "NM_004006.2:c.[101_103delinsTT;104A>T]";
        let parsed = crate::parse_hgvs(unbalanced).expect("parse");
        assert!(matches!(
            classify(input, &parsed, unbalanced),
            SweepOutcome::MemberChangesLength(_)
        ));
    }

    /// The direction [`names_the_sweep_transcript`] was added for, which its two
    /// neighbours do not cover: the genomic wrapper is **accepted**.
    ///
    /// `an_output_on_another_accession_is_unexpected_rather_than_unchanged` and
    /// `a_member_on_another_accession_is_unexpected` both test the reject side,
    /// and both would still pass against the bare `accession.full() != TRANSCRIPT`
    /// comparison this helper replaced. Only this test distinguishes them, so
    /// without it the helper is untested in the only direction it changed.
    ///
    /// The claim is that a wrapper is not a re-anchor (#1704): `checklist.md:20`
    /// makes `NC_000023.11(NM_004006.2)` mandatory once a description names an
    /// intronic position, and four rows of this sweep shift into an intron. Under
    /// a bare string comparison those read as "re-anchored onto another
    /// transcript" — the retype-in-disguise `classify` exists to catch — and they
    /// are not one: the inner accession is unchanged.
    ///
    /// Both call sites are exercised, because they are separate comparisons and a
    /// fix applied to one only is the likely regression.
    ///
    /// Note the wrapped row lands in `ShiftedInversion`, not `Unchanged`:
    /// `Unchanged` is a **string-identity** bucket (`rendered == input`) and the
    /// wrapper changes the string. That is the correct bucket — the description
    /// really did move — and it is a conformant one, which is the claim. Reading
    /// `Unchanged` as "same variant" is the mistake this note exists to head off.
    #[test]
    fn a_genomic_wrapper_is_not_a_re_anchor() {
        // The whole-description call site, in `classify`.
        let input = "NM_004006.2:c.101_104inv";
        let wrapped = "NC_000023.11(NM_004006.2):c.101_104inv";
        let parsed = crate::parse_hgvs(wrapped).expect("parse");
        assert_eq!(
            classify(input, &parsed, wrapped),
            SweepOutcome::ShiftedInversion,
            "the wrapper is the reference `checklist.md:20` requires, not a different \
             transcript — the inner accession is unchanged, so this is a conformant \
             re-rendering rather than `Unexpected`"
        );

        // The per-member call site, in `member_violation`.
        let split = "NC_000023.11(NM_004006.2):c.[101A>T;104A>T]";
        let parsed = crate::parse_hgvs(split).expect("parse");
        assert!(
            classify(input, &parsed, split).is_conformant(),
            "a wrapped member is on the sweep transcript too, and `member_violation` is a \
             second comparison that has to agree with the first"
        );

        // The control: a wrapper around a *different* transcript is still
        // rejected, so the helper is not simply ignoring the inner accession.
        let elsewhere = "NC_000023.11(NM_000088.3):c.101_104inv";
        let parsed = crate::parse_hgvs(elsewhere).expect("parse");
        assert!(
            matches!(
                classify(input, &parsed, elsewhere),
                SweepOutcome::Unexpected(_)
            ),
            "stripping the wrapper must not strip the check — the inner accession still decides"
        );
    }

    #[test]
    fn a_member_on_another_accession_is_unexpected() {
        let mixed = "NM_004006.2:c.[101A>T;104A>T]";
        let parsed = crate::parse_hgvs(mixed).expect("parse");
        assert!(classify(mixed, &parsed, mixed).is_conformant());

        let elsewhere = "NM_000088.3:c.[101A>T;104A>T]";
        let parsed = crate::parse_hgvs(elsewhere).expect("parse");
        assert!(matches!(
            classify(mixed, &parsed, elsewhere),
            SweepOutcome::Unexpected(_)
        ));
    }

    #[test]
    fn the_shared_normalize_config_is_lenient_and_three_prime() {
        let (errors, normalize) = sweep_normalize_config();
        assert_eq!(errors.mode, ErrorConfig::lenient().mode);
        assert_eq!(normalize.shuffle_direction, ShuffleDirection::ThreePrime);
    }

    #[test]
    fn the_tsv_round_trips_and_rejects_a_malformed_row() {
        let cases = PinnedCases {
            rows: vec![PinnedCase {
                input: "NM_004006.2:c.101_104inv".to_string(),
                output: "NM_004006.2:c.101_104inv".to_string(),
            }],
        };
        let rendered = cases.to_tsv("generated\ndo not hand-edit");
        assert!(rendered.starts_with("# generated\n# do not hand-edit\n"));
        assert_eq!(PinnedCases::from_tsv(&rendered).unwrap().rows, cases.rows);

        assert!(PinnedCases::from_tsv("no-tab-here").is_err());
        assert!(PinnedCases::from_tsv("input\t").is_err());
    }

    #[test]
    fn only_the_four_conformant_outcomes_are_conformant() {
        assert!(SweepOutcome::Unchanged.is_conformant());
        assert!(SweepOutcome::ShiftedInversion.is_conformant());
        assert!(SweepOutcome::PalindromicNoOp.is_conformant());
        assert!(SweepOutcome::Repartitioned(2).is_conformant());
        assert!(!SweepOutcome::MemberChangesLength("x".to_string()).is_conformant());
        assert!(!SweepOutcome::Retyped("substitution".to_string()).is_conformant());
        assert!(!SweepOutcome::Unexpected("x".to_string()).is_conformant());
    }
}
