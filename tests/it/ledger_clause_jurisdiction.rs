//! A ruling may not speak to an axis whose molecule directory it does not cite.
//!
//! # The failure mode this exists for
//!
//! `assets/hgvs-nomenclature` organises the recommendations by **molecule**:
//! `docs/recommendations/DNA/` covers the `g.`, `c.`, `n.` and `m.` axes,
//! `RNA/` covers `r.`, `protein/` covers `p.`. That directory is the
//! jurisdiction of every clause inside it — a `DNA/` clause is not authority for
//! a claim about `r.`, however exactly it is quoted and however precisely its
//! `file:line` resolves.
//!
//! The existing ledger guards all pass such a record. The spec-fixture
//! generator checks that each citation's quote is really at that `file:line`;
//! `clause_ruling_index.rs` checks that every cited clause is reachable from a
//! clause string; `ruling_citation_currency.rs` checks that prose does not
//! contradict a record's status. A citation can satisfy all three and still not
//! be authority for what the record says, because none of them compares the
//! **axes a record rules on** against the **molecules it cites**.
//!
//! Two instances were found on one day, which is what argued for a guard rather
//! than a third correction:
//!
//! * `delins-codon-carve-out-gap-one` listed `r.142_144delinsugg` and
//!   `r.[142c>u;144a>g]` in `applies_to` while its structured `clauses` were
//!   `DNA/delins.md:18`, `general.md:39` and `consultation/SVD-WG010.md:5`. It
//!   named `RNA/delins.md:18` and `RNA/substitution.md:18` in its rationale
//!   prose and cited neither, so a ruling that is `decided` rested its RNA reach
//!   on an uncited clause. The fix landed with this guard, and is a citation
//!   promotion: the two clauses the rationale already named became structured
//!   entries, and no verdict moved.
//! * A record proposed on PR #1682 cited six clauses, all under `DNA/` or
//!   `general.md`, and ruled in its rationale on `r.` over a coding transcript
//!   and on a non-coding `r.`. `RNA/delins.md` has no counterpart to the
//!   `DNA/delins.md:47` clause that record leans on. That PR is fixing its own
//!   record; this guard is what stops the next one.
//!
//! # What "speaks to" means, and why the rationale is read
//!
//! An axis prefix in an `applies_to` entry is the record's own statement of what
//! it rules on and is the strongest signal available. It is **not** sufficient:
//! the second instance above has a single `applies_to` entry, a `c.` string, and
//! its whole RNA problem lives in the rationale. So axis tokens in the rationale
//! count too — read only inside backticked spans, because that is how this ledger
//! writes an axis and because unbackticked prose is where the false positives
//! are.
//!
//! # What "cites" means, and where `general.md` fits
//!
//! This is the part that was decided from measurement rather than from taste,
//! because the two obvious readings both fail. Counting `general.md` as
//! universal authority makes the guard **toothless**: both known defects cite
//! `general.md`, so that rule flags neither. Counting it for nothing flags
//! records that rule entirely under molecule-agnostic authority
//! (`background/standards.md`, `recommendations/checklist.md`) and mention a
//! `c.` example in passing. Neither is usable.
//!
//! **The counts that argument was made on are not restated here**, because they
//! were measured on a 16-record ledger with both defects present and neither
//! condition survives: this PR fixes the first defect and the second lives on a
//! branch. A count quoted from a base nobody can rebuild is the failure mode
//! `clause_ruling_index.rs` documents at length. What is enforceable is that
//! **each half changes the verdict on the ledger as it stands**, and that is
//! what [`the_scan_is_not_vacuous`] asserts — drop either half and it fails.
//!
//! The rule has two halves:
//!
//! 1. **Jurisdiction attaches to the clause a record rules under.** If the
//!    `governing` clause is molecule-agnostic — `general.md`, `checklist.md`,
//!    `style.md`, anything under `background/`, `consultation/` or `versions/` —
//!    the record carries no molecule-scoped authority to exceed and is
//!    unconstrained. That is where `general.md` gets to be universal, and it is
//!    what keeps the records above out of the flag list — today
//!    `alignment-only-symbol-in-a-description`,
//!    `bare-transcript-intronic-position` and
//!    `absolute-prohibition-enforcement-stage`. A record with no
//!    `governing` clause at all is constrained only if it reaches into a
//!    molecule directory anyway, so that an `undecided` record — which the
//!    generator forbids from naming a governing clause — is not exempt by
//!    construction.
//! 2. **A molecule-agnostic clause carries the jurisdiction its own quoted text
//!    enumerates, and no more.** `general.md:44` reads "**exception**:
//!    deletions/duplications around exon/exon junctions using **c.**, **r.** or
//!    **n.** reference sequences" — it names the RNA axis, so it *is* RNA
//!    authority, and `exon-junction-dup-converge-from-the-far-side` scoping its
//!    ruling onto `r.` is properly cited rather than a defect. `general.md:34`,
//!    `:39` and `:43` name no axis at all, so they extend nothing. Without this
//!    half the guard reports that record too, which measurement says is wrong.
//!
//! # What this does NOT claim
//!
//! * **Not that a cited clause supports the ruling.** It compares directories,
//!   not arguments. A record can cite `RNA/repeated.md` for a delins question
//!   and pass. Whether the clause says what the record needs is review's job.
//! * **Not that the axes it reads are all the axes a record touches.** An axis
//!   discussed in unbackticked prose, or implied ("on a coding transcript"), is
//!   invisible here. This is a floor: it catches a record that *writes* an axis
//!   token it has no directory for.
//! * **Not that a quote is byte-exact.** The generator's citation check is a
//!   whitespace-collapsed substring match over the cited line range, so a green
//!   build is not evidence of exactness — check a quote against the spec file if
//!   a claim rests on it. Half 2 above reads the quote as the ledger stores it,
//!   which is the text a reader sees.
//! * **Not that an axis token is a claim rather than a disclaimer.** A record
//!   that writes `` `p.` `` only to say it does *not* rule on the protein axis
//!   reads identically to one that does. One ledger record does exactly that and
//!   is exempted by name in [`SCOPE_DISCLAIMERS`]; the exemption is
//!   shrink-only, so it cannot grow into a way of silencing the guard.
//! * **Not that `protein/` is exercised by the ledger.** No record cites a
//!   clause under it, so that arm of the mapping is exercised only by the
//!   synthetic cases at the bottom of this file. Said plainly rather than left
//!   for someone to discover, since a rule with an unexercised arm reads as
//!   coverage.

use std::collections::{BTreeMap, BTreeSet};

use crate::common::rulings::{self, Record};

/// A molecule, and the directory that is its jurisdiction.
///
/// Prefixes rather than exact paths: every file under `DNA/` is DNA authority,
/// and the spec adds pages without asking.
const MOLECULE_DIRECTORIES: &[(&str, &str)] = &[
    ("DNA", "docs/recommendations/DNA/"),
    ("RNA", "docs/recommendations/RNA/"),
    ("protein", "docs/recommendations/protein/"),
];

/// Axis prefix -> the molecule whose recommendations govern it.
///
/// `o.` is the circular-genome axis and `m.` the mitochondrial one; both are DNA
/// descriptions and both live under `DNA/`.
const AXIS_MOLECULES: &[(char, &str)] = &[
    ('g', "DNA"),
    ('c', "DNA"),
    ('n', "DNA"),
    ('m', "DNA"),
    ('o', "DNA"),
    ('r', "RNA"),
    ('p', "protein"),
];

/// The molecule an axis prefix belongs to, if the character is an axis at all.
fn molecule_of_axis(prefix: char) -> Option<&'static str> {
    AXIS_MOLECULES
        .iter()
        .find(|(axis, _)| *axis == prefix)
        .map(|(_, molecule)| *molecule)
}

/// The molecule a clause's directory places it under, or `None` when the clause
/// is molecule-agnostic.
///
/// Takes the citation as the ledger spells it, line locator and all, because
/// only the leading path matters.
fn molecule_of_clause(clause: &str) -> Option<&'static str> {
    MOLECULE_DIRECTORIES
        .iter()
        .find(|(_, directory)| clause.starts_with(directory))
        .map(|(molecule, _)| *molecule)
}

/// Every axis prefix appearing in `text`, as `<letter>.`.
///
/// The preceding character must not be alphanumeric, `_` or `.`, which is what
/// keeps two classes of non-axis text out. **Spec filenames**: `inversion.md`
/// and `substitution.md` both end in `n` before the dot, so without the check
/// every citation of either would read as a claim about the `n.` axis.
/// **Abbreviations**: `e.g.` ends in `g` before a dot, and the recommendations
/// use it constantly. Measured over the ledger, six `<dot><axis>.` occurrences:
/// five inside citation quotes, one in a rationale (outside backticks, so
/// unread).
///
/// ```text
/// python3 -c "import json,re;d=json.load(open('tests/fixtures/grammar/hgvs_spec_normalization_overrides.json'));\
/// print(sum(len(re.findall(r'(?<=\.)[gcnmorp]\.',c['quote'])) for r in d['rulings'] for c in r['clauses']))"
/// ```
///
/// The second class is the one that matters for correctness rather than noise,
/// because it over-reads in the *unsafe* direction. `molecules_with_authority`
/// reads the quote of a molecule-agnostic clause, so an `e.g.` in a
/// `general.md` or `checklist.md` quote would hand the record DNA authority it
/// was never granted and quietly stop the guard catching a DNA overreach. **No
/// agnostic quote carries one today** — all five sit on clauses under a
/// molecule directory, whose quotes are never read — so this moves no verdict
/// and is a latent hazard closed rather than a defect fixed.
///
/// [`filenames_are_not_read_as_axes`] pins both classes.
fn axes_in(text: &str) -> BTreeSet<char> {
    let bytes: Vec<char> = text.chars().collect();
    let mut found = BTreeSet::new();
    for (index, window) in bytes.windows(2).enumerate() {
        let (candidate, dot) = (window[0], window[1]);
        if dot != '.' || molecule_of_axis(candidate).is_none() {
            continue;
        }
        let preceded_by_word = index
            .checked_sub(1)
            .map(|before| {
                let previous = bytes[before];
                previous.is_alphanumeric() || previous == '_' || previous == '.'
            })
            .unwrap_or(false);
        if !preceded_by_word {
            found.insert(candidate);
        }
    }
    found
}

/// Every axis prefix appearing inside a backticked span of `text`.
///
/// Scoped to backticks because that is how this ledger spells an axis, and
/// because unbackticked rationale prose is discursive enough that a bare
/// `<letter>.` there is far more often a sentence ending than an axis.
fn axes_in_backticks(text: &str) -> BTreeSet<char> {
    let mut found = BTreeSet::new();
    let mut inside = false;
    for span in text.split('`') {
        if inside {
            found.extend(axes_in(span));
        }
        inside = !inside;
    }
    found
}

/// One reason a record was read as speaking to a molecule, for the diagnostic.
type Evidence = BTreeMap<&'static str, Vec<String>>;

/// The molecules a record speaks to, with why.
fn molecules_spoken_to(record: &Record) -> Evidence {
    let mut spoken: Evidence = BTreeMap::new();
    for description in &record.applies_to {
        for axis in axes_in(description) {
            if let Some(molecule) = molecule_of_axis(axis) {
                spoken
                    .entry(molecule)
                    .or_default()
                    .push(format!("`applies_to` entry {description:?}"));
            }
        }
    }
    for axis in axes_in_backticks(&record.rationale) {
        if let Some(molecule) = molecule_of_axis(axis) {
            spoken
                .entry(molecule)
                .or_default()
                .push(format!("a backticked `{axis}.` in the rationale"));
        }
    }
    spoken
}

/// The molecules a record has authority for, with why.
///
/// Both halves of the rule's second point: a clause under a molecule directory
/// supplies that molecule, and a molecule-agnostic clause supplies exactly the
/// molecules its own quote names by axis token.
fn molecules_with_authority(record: &Record) -> Evidence {
    let mut authority: Evidence = BTreeMap::new();
    for citation in &record.citations {
        match molecule_of_clause(&citation.clause) {
            Some(molecule) => authority
                .entry(molecule)
                .or_default()
                .push(format!("cites {}", citation.clause)),
            None => {
                for axis in axes_in(&citation.quote) {
                    if let Some(molecule) = molecule_of_axis(axis) {
                        authority.entry(molecule).or_default().push(format!(
                            "{} quotes the `{axis}.` axis by name",
                            citation.clause
                        ));
                    }
                }
            }
        }
    }
    authority
}

/// Whether a record's authority is molecule-scoped at all.
///
/// A record ruling under a molecule-agnostic clause is unconstrained; so is one
/// that names no governing clause and reaches into no molecule directory. See
/// half 1 of the rule in this module's docs.
fn is_molecule_scoped(record: &Record) -> bool {
    match &record.governing {
        Some(governing) => molecule_of_clause(governing).is_some(),
        None => record
            .citations
            .iter()
            .any(|c| molecule_of_clause(&c.clause).is_some()),
    }
}

/// One molecule a record speaks to without authority for it.
#[derive(Debug)]
struct Overreach {
    id: String,
    molecule: &'static str,
    spoken_because: Vec<String>,
    authority_for: Vec<&'static str>,
}

/// Every overreach in `records`, in ledger order.
fn overreaches(records: &[Record]) -> Vec<Overreach> {
    let mut found = Vec::new();
    for record in records {
        if !is_molecule_scoped(record) {
            continue;
        }
        let authority = molecules_with_authority(record);
        for (molecule, why) in molecules_spoken_to(record) {
            if !authority.contains_key(molecule) {
                found.push(Overreach {
                    id: record.id.clone(),
                    molecule,
                    spoken_because: why,
                    authority_for: authority.keys().copied().collect(),
                });
            }
        }
    }
    found
}

/// Records that name a molecule in prose **only to disclaim it**, by id,
/// molecule, and the sentence that does the disclaiming.
///
/// The rule reads an axis token as a claim because prose cannot be parsed for
/// intent, and a record that writes `` `p.` `` to say "not this axis" is
/// indistinguishable from one that rules on it. Rather than weaken the
/// tokenizer — which would silently drop the bare-label form across the whole
/// ledger, and that form is how the known defects could recur — the exceptions
/// are named here with their reason.
///
/// [`every_scope_disclaimer_is_still_needed`] makes the list **shrink-only**: a
/// row that stops being needed fails the test and must be deleted, so this
/// cannot become a way of silencing the guard. Adding a row is a legitimate
/// answer; doing it without a quoted sentence is not.
const SCOPE_DISCLAIMERS: &[(&str, &str, &str)] = &[
    (
        "projection-codon-exception-is-decided-by-the-rendered-axis",
        "protein",
        "its `SCOPE, STATED NARROWLY` paragraph reads \"The `p.` axis is the protein consequence \
         itself.\" — the record's only `p.`, and a statement that the protein axis is out of scope \
         rather than a ruling on it. The same paragraph disclaims `r.` and `n.` in the same breath; \
         `r.` needs no row only because `general.md:44`'s quote names that axis",
    ),
    (
        "delins-merge-vs-individual-gap-two-or-more",
        "RNA",
        "its 2026-08-17 (#2155) cross-reference paragraph reads \"That record originally scoped \
         the axis question to `c.` alone (2026-08-11) and was WIDENED to every DNA axis — \
         `c./g./m./n.`; `r.` still out — by operator ruling 2026-08-17 (#2155)\" — a statement \
         that the RNA axis stays out of the (delegated) axis question, not a ruling on it",
    ),
    (
        "delins-payload-coincidence-carve-out-is-coding-dna-scoped",
        "RNA",
        "its widen paragraph reads \"`r.` stays OUT, on the same jurisdiction ground the \
         2026-08-11 text already gives … a `DNA/` document has no authority over the RNA axis \
         regardless of how widely its DNA-side scope is drawn, and `RNA/delins.md` states no \
         `:47` counterpart of its own\" — a disclaimer, not a ruling on `r.`. The 2026-08-11 text \
         it supersedes made the identical disclaimer (\"It rules on the coding DNA axis alone. It \
         makes no ruling about the RNA axis…\") and its own WHAT DOES NOT MOVE paragraph repeats \
         it a third time (\"`r.` moves nothing: `RefShape::all()` has no RNA shape\")",
    ),
    (
        "unequal-length-block-a-placed-gap-is-not-a-separation",
        "RNA",
        "its 2026-08-17 (#2155) supersession paragraph reads \"The only axis this rule still \
         declines is `r.`, on the same jurisdiction ground stated throughout this ledger's \
         DNA-scoped records — a `DNA/` clause cannot scope the RNA axis regardless of how widely \
         its DNA-side scope is drawn, and `RNA/delins.md` states no `:47` counterpart of its \
         own\" — a disclaimer, not a ruling on `r.`",
    ),
];

/// The disclaiming sentence for `(id, molecule)`, when the pair is listed.
fn scope_disclaimer(id: &str, molecule: &str) -> Option<&'static str> {
    SCOPE_DISCLAIMERS
        .iter()
        .find(|(record, listed, _)| *record == id && *listed == molecule)
        .map(|(_, _, why)| *why)
}

// --------------------------------------------------------------------------
// The guard.
// --------------------------------------------------------------------------

/// **The assertion this file exists to make.**
#[test]
fn no_record_speaks_to_a_molecule_it_cannot_cite() {
    let found: Vec<Overreach> = overreaches(&rulings::records())
        .into_iter()
        .filter(|o| scope_disclaimer(&o.id, o.molecule).is_none())
        .collect();
    let report: Vec<String> = found
        .iter()
        .map(|o| {
            format!(
                "\n  {} rules on the {} axes ({}) but its authority reaches only {:?}",
                o.id,
                o.molecule,
                o.spoken_because.join("; "),
                o.authority_for
            )
        })
        .collect();
    assert!(
        found.is_empty(),
        "a ruling names an axis whose molecule directory it does not cite:{}\n\n\
         `assets/hgvs-nomenclature/docs/recommendations/` is organised by molecule, and that \
         directory is the jurisdiction of every clause inside it — a `DNA/` clause cannot scope a \
         claim about `r.`. Cite a clause under the molecule's own directory, or narrow what the \
         record says. If the authority really is a molecule-agnostic clause that names the axis \
         itself (as `general.md:44` does), cite that clause and quote the part that names it.",
        report.join("")
    );
}

/// The scan reads what it claims to, and both carve-outs are load-bearing.
///
/// Modelled on `clause_ruling_index.rs`'s `the_index_is_not_vacuous`. Each
/// assertion below pins a population that a broken reader would empty: a ledger
/// that failed to yield axes, a molecule map that matched no directory, or a
/// carve-out that no record exercises. A rule whose exemptions are dead code
/// reads as nuance and provides none.
#[test]
fn the_scan_is_not_vacuous() {
    let records = rulings::records();
    assert!(
        records.len() >= 10,
        "the ledger reader returned {} records",
        records.len()
    );

    let speaking: Vec<&Record> = records
        .iter()
        .filter(|r| !molecules_spoken_to(r).is_empty())
        .collect();
    assert!(
        speaking.len() >= 8,
        "only {} records were read as speaking to any molecule — the axis tokenizer is broken, \
         and every jurisdiction check below would pass vacuously",
        speaking.len()
    );

    let molecules_reached: BTreeSet<&str> = records
        .iter()
        .flat_map(|r| molecules_with_authority(r).into_keys())
        .collect();
    assert!(
        molecules_reached.len() >= 2,
        "citations resolved to only {molecules_reached:?} — with one molecule in play the \
         comparison cannot fail"
    );

    // The population the guard actually judges. Without this, an
    // `is_molecule_scoped` that answered `false` for everything would make
    // `no_record_speaks_to_a_molecule_it_cannot_cite` pass over an empty set
    // while every other assertion here still held — including the half-1 pin
    // below, which such a break would satisfy *harder*.
    let constrained: Vec<&String> = records
        .iter()
        .filter(|r| is_molecule_scoped(r) && !molecules_spoken_to(r).is_empty())
        .map(|r| &r.id)
        .collect();
    assert!(
        constrained.len() >= 8,
        "only {} records are both molecule-scoped and read as speaking to a molecule, so the \
         guard is judging almost nothing: {constrained:?}",
        constrained.len()
    );

    // Half 1: at least one record rules under molecule-agnostic authority.
    let unscoped: Vec<&String> = records
        .iter()
        .filter(|r| !is_molecule_scoped(r))
        .map(|r| &r.id)
        .collect();
    assert!(
        !unscoped.is_empty(),
        "no record is exempted by its governing clause, so half 1 of the rule is dead code. It \
         exists because three records rule under `checklist.md` or `background/standards.md` and \
         mention a `c.` example; if the ledger no longer has one, delete the carve-out rather \
         than leaving it unexercised"
    );

    // Half 2: at least one molecule is supplied by an axis-naming agnostic quote.
    let by_quote: Vec<(&String, &'static str)> = records
        .iter()
        .flat_map(|record| {
            let direct: BTreeSet<&str> = record
                .citations
                .iter()
                .filter_map(|c| molecule_of_clause(&c.clause))
                .collect();
            molecules_with_authority(record)
                .into_keys()
                .filter(move |m| !direct.contains(m))
                .map(move |m| (&record.id, m))
        })
        .collect();
    assert!(
        !by_quote.is_empty(),
        "no molecule-agnostic clause supplies a molecule through its own quote, so half 2 of the \
         rule is dead code. `general.md:44` naming `**c.**, **r.** or **n.**` is the case it was \
         written for"
    );
}

/// Every [`SCOPE_DISCLAIMERS`] row still exempts something, so the list cannot
/// outlive the records it was written for.
///
/// Shrink-only in the same sense as `LEDGER_EXEMPT` in
/// `generator_completeness.rs`: a row whose record has been deleted, re-cited,
/// or reworded so the token is gone fails here and must be removed. That is the
/// only thing standing between a named exception and a general mute button.
#[test]
fn every_scope_disclaimer_is_still_needed() {
    let raw = overreaches(&rulings::records());
    for (id, molecule, why) in SCOPE_DISCLAIMERS {
        assert!(
            raw.iter().any(|o| o.id == *id && o.molecule == *molecule),
            "`{id}` is no longer read as speaking to {molecule} without authority, so its \
             disclaimer row is dead. Delete it rather than leaving an unexercised exemption. It \
             was there because {why}"
        );
    }
}

/// The `r.` half of `delins-codon-carve-out-gap-one` is cited, not only narrated.
///
/// A named pin on top of the general guard, because this record is the reason
/// the guard exists and a general assertion says nothing about which record it
/// was passing for. The record is `decided`, so its RNA rows rest on structured
/// authority or on nothing.
#[test]
fn the_codon_carve_out_record_cites_rna_for_its_rna_rows() {
    let records = rulings::records();
    let record = records
        .iter()
        .find(|r| r.id == "delins-codon-carve-out-gap-one")
        .expect("the ledger still carries this ruling");

    assert!(
        record
            .applies_to
            .iter()
            .any(|description| description.starts_with("r.")),
        "this pin assumes the record still declares RNA descriptions; if that changed, the RNA \
         citation requirement below may no longer be the right assertion"
    );

    let rna: Vec<&str> = record
        .citations
        .iter()
        .map(|c| c.clause.as_str())
        .filter(|clause| molecule_of_clause(clause) == Some("RNA"))
        .collect();
    assert!(
        !rna.is_empty(),
        "the ruling declares `r.` descriptions in `applies_to` and cites no clause under \
         `docs/recommendations/RNA/`"
    );
}

// --------------------------------------------------------------------------
// The rule, on synthetic records.
//
// Plain `#[test]`, not a `#[cfg(test)] mod tests` — this tree is an integration
// binary and compiles without `cfg(test)`, so a gated module would never run.
//
// The tests above read the committed ledger, so none of them can exercise a
// record the ledger does not contain. These drive the rule directly, which is
// the only way to pin that it *fires* and that each carve-out is doing what its
// justification says.
// --------------------------------------------------------------------------

/// Build a record from its parts, so a case reads as the shape it is testing.
fn record(
    id: &str,
    governing: Option<&str>,
    clauses: &[(&str, &str)],
    applies_to: &[&str],
    rationale: &str,
) -> Record {
    Record {
        id: id.to_string(),
        status: "decided".to_string(),
        // Neither field is read by the jurisdiction rule, which compares a
        // record's cited clauses against the axes its `applies_to` names. They
        // are filled rather than defaulted because `Record` has no `Default`:
        // a record assembled from parts must state every part it carries.
        question: format!("A synthetic record for {id}."),
        // Not read by the jurisdiction rule either; see the comment above.
        summary: None,
        equivalence_classes: Vec::new(),
        guard: rulings::Guard::Declined("a synthetic record, enforced by this file".to_string()),
        // Not a house choice: every case here names a governing clause, which is
        // the shape a house choice may not take.
        house_choice: None,
        rationale: rationale.to_string(),
        applies_to: applies_to.iter().map(|s| s.to_string()).collect(),
        citations: clauses
            .iter()
            .map(|(clause, quote)| rulings::Citation {
                clause: clause.to_string(),
                quote: quote.to_string(),
                role: if Some(*clause) == governing {
                    rulings::Role::Governing
                } else {
                    rulings::Role::Cited
                },
            })
            .collect(),
        governing: governing.map(str::to_string),
    }
}

const DNA_DELINS: (&str, &str) = (
    "docs/recommendations/DNA/delins.md:18",
    "**exception**: two variants separated by one nucleotide",
);
const RNA_DELINS: (&str, &str) = (
    "docs/recommendations/RNA/delins.md:18",
    "**exception**: two variants separated by one nucleotide",
);

/// The shape of both defects: a DNA-governed ruling that rules on `r.`.
#[test]
fn a_dna_governed_ruling_may_not_rule_on_rna() {
    let found = overreaches(&[record(
        "a-dna-ruling-reaching-rna",
        Some(DNA_DELINS.0),
        &[DNA_DELINS],
        &["r.142_144delinsugg"],
        "It also settles the coding `r.` case.",
    )]);
    assert_eq!(
        found.len(),
        1,
        "expected exactly one overreach, got {found:?}"
    );
    assert_eq!(found[0].molecule, "RNA");
}

/// The fix, and the negative control for the case above: citing the RNA clause
/// clears it, with nothing else changed.
#[test]
fn citing_the_rna_clause_clears_it() {
    let found = overreaches(&[record(
        "a-dna-ruling-citing-rna",
        Some(DNA_DELINS.0),
        &[DNA_DELINS, RNA_DELINS],
        &["r.142_144delinsugg"],
        "It also settles the coding `r.` case.",
    )]);
    assert!(found.is_empty(), "{found:?}");
}

/// The rationale is read, not only `applies_to` — the second defect's whole
/// problem lives there, its single `applies_to` entry being a `c.` string.
#[test]
fn an_rna_claim_made_only_in_the_rationale_is_caught() {
    let found = overreaches(&[record(
        "a-ruling-whose-rna-claim-is-prose",
        Some(DNA_DELINS.0),
        &[DNA_DELINS],
        &["LRG_199t1:c.850_901delinsTTCCTCGATGCCTG"],
        "On a coding transcript the same reasoning gives `r.850_901delinsuuccucg`.",
    )]);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0].molecule, "RNA");
}

/// Half 1: a record ruling under a molecule-agnostic clause is unconstrained.
///
/// Without this the guard reports `alignment-only-symbol-in-a-description`,
/// `bare-transcript-intronic-position` and `absolute-prohibition-enforcement-stage`,
/// all three of which rule under `checklist.md` or `background/standards.md`
/// and use a `c.` string only as an example.
#[test]
fn a_ruling_under_molecule_agnostic_authority_is_unconstrained() {
    let checklist = (
        "docs/recommendations/checklist.md:20",
        "a description should be based on a reference sequence",
    );
    let found = overreaches(&[record(
        "a-ruling-under-the-checklist",
        Some(checklist.0),
        &[checklist],
        &[],
        "Strict input hygiene refuses `NM_TEST.1:c.20+2del`.",
    )]);
    assert!(found.is_empty(), "{found:?}");
}

/// Half 2: a molecule-agnostic clause whose quote names an axis supplies that
/// molecule. This is `general.md:44`, and it is why
/// `exon-junction-dup-converge-from-the-far-side` scoping onto `r.` is cited
/// rather than an overreach.
#[test]
fn an_agnostic_clause_that_names_an_axis_supplies_that_molecule() {
    let general = (
        "docs/recommendations/general.md:44",
        "**exception**: deletions/duplications around exon/exon junctions using **c.**, **r.** \
         or **n.** reference sequences",
    );
    let duplication = (
        "docs/recommendations/DNA/duplication.md:26",
        "the variant is described as `LRG_199t1:c.3921dup`",
    );
    let found = overreaches(&[record(
        "a-junction-ruling-scoped-onto-rna",
        Some(duplication.0),
        &[duplication, general],
        &[],
        "SCOPE. Deletions and duplications, on the `c.` and `n.` axes, and on `r.`.",
    )]);
    assert!(found.is_empty(), "{found:?}");
}

/// …and the control for it: an agnostic clause that names no axis extends
/// nothing. `general.md:34` and `:43` are the two the second defect leaned on.
#[test]
fn an_agnostic_clause_that_names_no_axis_extends_nothing() {
    let general = (
        "docs/recommendations/general.md:34",
        "two variants separated by one or more nucleotides should be described individually",
    );
    let found = overreaches(&[record(
        "a-ruling-leaning-on-silent-general-prose",
        Some(DNA_DELINS.0),
        &[DNA_DELINS, general],
        &[],
        "And on a non-coding `r.` the same follows.",
    )]);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0].molecule, "RNA");
}

/// An `undecided` record cannot name a governing clause, so it must not be
/// exempt merely for having none.
#[test]
fn a_record_with_no_governing_clause_is_still_scoped_by_what_it_cites() {
    let found = overreaches(&[record(
        "an-undecided-record-reaching-rna",
        None,
        &[DNA_DELINS],
        &["r.142_144delinsugg"],
        "Unsettled.",
    )]);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0].molecule, "RNA");
}

/// The `protein/` arm, which no ledger record currently exercises.
#[test]
fn the_protein_arm_behaves_like_the_others() {
    let protein = (
        "docs/recommendations/protein/delins.md:21",
        "two variants affecting the same amino acid",
    );
    let reaching = overreaches(&[record(
        "a-dna-ruling-reaching-protein",
        Some(DNA_DELINS.0),
        &[DNA_DELINS],
        &["p.(Arg48Trp)"],
        "",
    )]);
    assert_eq!(reaching.len(), 1, "{reaching:?}");
    assert_eq!(reaching[0].molecule, "protein");

    let cited = overreaches(&[record(
        "a-dna-ruling-citing-protein",
        Some(DNA_DELINS.0),
        &[DNA_DELINS, protein],
        &["p.(Arg48Trp)"],
        "",
    )]);
    assert!(cited.is_empty(), "{cited:?}");
}

/// The tokenizer's own guard: neither a spec filename nor an abbreviation is an
/// axis.
///
/// `inversion.md` and `substitution.md` end in `n` before the dot and
/// `delins.md` in `s`; without the preceding-character check, every citation of
/// the first two would read as a claim about the `n.` axis and the ledger would
/// be one large false positive.
///
/// `e.g.` is the second class and the one that could weaken rather than merely
/// noise up the guard — see [`axes_in`]. A quote is read for authority only
/// when the clause it belongs to is molecule-agnostic, which is the case this
/// must never get wrong, so the pinned string is a real ledger quote rather
/// than a paraphrase of one.
#[test]
fn filenames_are_not_read_as_axes() {
    assert!(axes_in("docs/recommendations/DNA/inversion.md:20").is_empty());
    assert!(axes_in("docs/recommendations/DNA/substitution.md:32").is_empty());
    assert!(axes_in("docs/recommendations/DNA/delins.md:47").is_empty());
    assert!(axes_in("general.md:34").is_empty());

    // Abbreviations. The dot before the axis letter is what disqualifies it.
    // The first is `codon-carve-out-shape-restriction`'s quote of
    // `DNA/duplication.md:18`, verbatim from the ledger.
    assert!(axes_in(
        "when a variant can be described as a duplication, it **must** be described as a \
         duplication and not as, e.g., an insertion"
    )
    .is_empty());
    assert!(axes_in("e.g. two variants").is_empty());

    // And the positives it must still see, including one at the start of the
    // string, one after a `:` accession separator, and one after a `/`, which
    // is how the spec's own prose writes a pair of axes.
    assert_eq!(axes_in("r.142_144delinsugg"), BTreeSet::from(['r']));
    assert_eq!(axes_in("LRG_199t1:c.3921dup"), BTreeSet::from(['c']));
    assert_eq!(
        axes_in("the `c.` and `n.` axes"),
        BTreeSet::from(['c', 'n'])
    );
    assert_eq!(axes_in("the c./r. pair"), BTreeSet::from(['c', 'r']));
}

/// Backtick scoping: an axis token outside backticks is not read.
#[test]
fn only_backticked_axes_are_read_from_a_rationale() {
    assert!(axes_in_backticks("the RNA axis is r. throughout").is_empty());
    assert_eq!(
        axes_in_backticks("the RNA axis is `r.` throughout"),
        BTreeSet::from(['r'])
    );
    // An unterminated backtick opens a span that runs to the end of the text;
    // pinned so the split-based reader's behaviour is a decision, not a
    // surprise. Over-reading here is the safe direction — it can only add a
    // molecule the record must then cite.
    assert_eq!(
        axes_in_backticks("trailing `r.142del"),
        BTreeSet::from(['r'])
    );
}
