//! The authored-inversion sweep: schema, generator parameters, and outcome
//! classification.
//!
//! # What the sweep is
//!
//! 2,075 authored `inv` descriptions on one real transcript — `NM_004006.2`,
//! the DMD coding sequence the HGVS recommendations use for their own worked
//! inversion examples — enumerated as
//! `NM_004006.2:c.{start}_{start + len - 1}inv` for every `start` in
//! `101, 108, …, 2996` and every `len` in `{4, 6, 8, 12, 16}`.
//!
//! The parameters live here rather than in the test so that the generator that
//! captures the reference slice and the gate that replays it cannot enumerate
//! different corpora — the same single-source-of-truth arrangement as
//! [`super::case_harvest`] and [`super::spec_worked_examples`].
//!
//! # What is pinned, and which kind of record each pin is
//!
//! Three layers, and they are **not** redundant — each catches something the
//! others cannot:
//!
//! 1. **Exact outputs for all 2,075 rows** (`cases.tsv`). Characterization
//!    pins. They are the only layer that notices the 3'-shift arithmetic
//!    moving: `c.101_108inv -> c.102_107inv` silently becoming
//!    `-> c.103_106inv` keeps the same outcome bucket, keeps the census total,
//!    and never splits — so layers 2 and 3 stay green while a string a
//!    downstream consumer stores has moved.
//! 2. **The invariant**, [`SweepOutcome`]. An authored `inv` normalizes either
//!    to an `inv` (possibly 3'-shifted or re-typed) or to a no-op `=` when the
//!    span equals its own reverse complement. It is never split into allele
//!    members and never silently becomes a substitution, deletion or
//!    duplication. This is what makes a careless re-blessing of layer 1 still
//!    fail.
//! 3. **The census** — the three bucket counts, as named constants in the gate.
//!    A bulk re-bless that moves the *distribution* is visible in a diff and has
//!    to be argued for.
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
//! `DNA/inversion.md:17`) and may retype it when the replacement turns out to
//! equal the reference (a palindromic span is its own reverse complement, so
//! nothing changed) — but it may not conclude that the replacement is something
//! other than the reverse complement of that span, which is what a split into
//! members or a retype to `del`/`dup`/`sub` would assert.

use std::path::Path;

use crate::hgvs::edit::NaEdit;
use crate::HgvsVariant;

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

/// What normalization did to one authored inversion.
///
/// The first three are the conformant outcomes and the census buckets; the last
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
    /// Split into allele members. `DNA/inversion.md:5` defines the whole span's
    /// replacement as its reverse complement; a split asserts something else.
    SplitIntoMembers(usize),
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
        )
    }

    /// Short label for census reporting.
    pub fn label(&self) -> &'static str {
        match self {
            SweepOutcome::Unchanged => "unchanged",
            SweepOutcome::ShiftedInversion => "shifted-still-inv",
            SweepOutcome::PalindromicNoOp => "palindromic-no-op",
            SweepOutcome::SplitIntoMembers(_) => "SPLIT-INTO-MEMBERS",
            SweepOutcome::Retyped(_) => "RETYPED",
            SweepOutcome::Unexpected(_) => "UNEXPECTED",
        }
    }
}

/// Classify a normalized result against the authored input it came from.
///
/// `rendered` is the output string and `input` the description handed in;
/// [`SweepOutcome::Unchanged`] is decided by string equality with `input`, and
/// everything else by the output's shape.
pub fn classify(input: &str, normalized: &HgvsVariant, rendered: &str) -> SweepOutcome {
    if let HgvsVariant::Allele(allele) = normalized {
        return SweepOutcome::SplitIntoMembers(allele.variants.len());
    }
    let HgvsVariant::Cds(cds) = normalized else {
        return SweepOutcome::Unexpected(format!(
            "output is not a c. description on {TRANSCRIPT}: {rendered}"
        ));
    };
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

/// The reverse complement of an ASCII DNA string, uppercased.
///
/// Used to decide whether a span is its own reverse complement, which is the
/// only thing that makes a [`SweepOutcome::PalindromicNoOp`] correct.
pub fn reverse_complement(bases: &str) -> String {
    bases
        .bytes()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
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
        // 415 starts (101..3000 step 7) x 5 lengths.
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
                start: 2996,
                length: 16
            }
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
    fn only_the_three_conformant_outcomes_are_conformant() {
        assert!(SweepOutcome::Unchanged.is_conformant());
        assert!(SweepOutcome::ShiftedInversion.is_conformant());
        assert!(SweepOutcome::PalindromicNoOp.is_conformant());
        assert!(!SweepOutcome::SplitIntoMembers(2).is_conformant());
        assert!(!SweepOutcome::Retyped("substitution".to_string()).is_conformant());
        assert!(!SweepOutcome::Unexpected("x".to_string()).is_conformant());
    }
}
