//! Adjudication-shaped prose in `src/` and `tests/` must point at a ruling
//! record — a shrink-only ratchet against the thing that has already happened
//! once.
//!
//! # The failure this exists for, measured
//!
//! One implementer's choice — that a reference base is unchanged iff **every**
//! minimum-edit-distance alignment matches it — had no home in the ledger,
//! because a `decided` record was required to name a governing spec clause and
//! this choice has none. So it was written down where it was needed instead, and
//! by 2026-08-14 it existed as **six** independent prose restatements across
//! `src/` and `tests/`. They drifted apart, three of them into the flatly false
//! generalisation that an equal-length block's column correspondence is unique;
//! reconciling them cost several days of adjudication. An HGVS
//! committee-perspective review of the same rule called it "an implementer's
//! choice the recommendations neither require nor forbid" and objected to its
//! being presented as compliance.
//!
//! The ledger now has a home for such a choice — `house_choice`, see
//! `generate_spec_fixture`'s `overrides::HouseChoice` — and the known sites have
//! been pointed at their record. A one-time cleanup regrows. This is the ratchet
//! that makes the seventh restatement cost a deliberate act.
//!
//! # What it checks
//!
//! Every contiguous comment block under `src/` and `tests/` whose text matches
//! one of [`ADJUDICATION_PHRASES`] must **name a ledger record id** somewhere in
//! that same block, or the block's file must carry a row in [`REGROWTH_EXEMPT`].
//!
//! The exemption list is **shrink-only**, on the pattern `LEDGER_EXEMPT` in
//! `tests/it/generator_completeness.rs` already established here: a row whose
//! site no longer flags fails the test and must be deleted, so the list cannot
//! rot into a general excuse for files nobody reads.
//!
//! # What it is NOT, because the obvious design does not survive measurement
//!
//! The natural way to build this is the inverse of `clause_ruling_index.rs`:
//! flag a comment that cites a clause the ledger governs and names no record
//! within N lines. **That predicate was measured over this tree and is
//! unusable** — 707 comment lines in 144 files, and 517 lines carrying a
//! choice/precedence verb with no nearby record. Most are the legitimate habit
//! of quoting the clause a piece of code implements. A guard whose first run
//! flags 700 lines is argued with once and disabled, and then there is no guard.
//!
//! The allowlist form does not transfer at that scale either, and the reason is
//! worth keeping: `LEDGER_EXEMPT` in `generator_completeness.rs` works because
//! it is ~15 rows keyed on a **file path**, which is a stable identifier. A
//! clause-keyed scan would need hundreds of rows keyed on a **site** —
//! `file:line` churns on every edit above it, `(file, sentence)` churns exactly
//! when you reword, and bare `file` is too coarse (one row on `src/normalize/merge.rs`
//! would cover ~124 flagged lines and hide every distinct site inside one
//! exemption).
//!
//! So what transfers from `generator_completeness.rs` is the **ratchet
//! discipline** — adding a row is a legitimate, reviewed act; silence is not —
//! and not the predicate. This file keeps the discipline and replaces the
//! predicate with a narrow phrase list, measured below. Its first run over the
//! tree flagged **three** blocks in total, of which one is exempt. The
//! coarse-key hazard is closed rather than accepted: an exemption row may cover
//! **one** flagged block, and a second block in the same file fails the test
//! (see [`adjudication_shaped_prose_points_at_a_ruling_record_or_is_allowlisted`]).
//!
//! # Read the reach honestly — this is a floor, not a proof
//!
//! Stating the limits is not a disclaimer, it is the operating manual. A guard
//! whose reach is overstated gets cited for coverage it does not provide, which
//! is the same defect as the prose it polices.
//!
//! * **It cannot reach the clause-citing class at all**, and that class is real:
//!   an audit of this tree on 2026-08-14 found fifteen sites where prose decides
//!   something, and **nine** of them are invisible to every existing guard. The
//!   707 lines above are not covered here and are not covered anywhere. This
//!   file closes one shape, not the problem.
//!
//! * **It is a substring scan over comment text.** It cannot tell an
//!   adjudication from a sentence that happens to contain one of the phrases,
//!   and it cannot see an adjudication written in words nobody has used yet.
//! * **[`ADJUDICATION_PHRASES`] is deliberately narrow and was chosen by
//!   measurement, not by taste.** Every candidate was run over the whole tree
//!   first. Phrases that read adjudicative but are ordinary here — "is
//!   canonical" (26 blocks), "the canonical form is" (24), "must be described
//!   as" (4, all of them quoting `DNA/duplication.md:18`), "the spec forbids"
//!   (10), "is a fact about the" (8) — were **rejected**, because a guard with
//!   dozens of false positives is argued with once and disabled the next time.
//!   What survived is the family that states an *ontological* claim about
//!   representation: that something is a fact rather than a choice, or that a
//!   column correspondence is unique. That family had **three** hits in the
//!   whole tree, which is what makes it usable.
//! * **Comments only.** Prose inside an `assert!`/`panic!` message is not
//!   scanned, and that is a known miss with a named instance:
//!   `reported_confluence_pairs.rs`'s assertion message used to *instruct* a
//!   future maintainer to delete an assertion, which is adjudication in the most
//!   consequential possible place. Reaching it needs the string literals
//!   attributed to their enclosing item, which is parsing rather than reading.
//! * **The pointer test is "does this block name a record id", not "does it name
//!   the *right* record".** Relevance is a judgement; membership is checkable.
//! * **Quoting a claim in order to withdraw it still flags.** That is why
//!   `issue_1284_transcript_axis_collision.rs` is on the exemption list rather
//!   than being a bug in the detector — it quotes the withdrawn claim verbatim,
//!   which is the right way to record a withdrawal and reads identically to
//!   asserting it.
//!
//! # Siblings, and the gap between them
//!
//! * `ruling_citation_currency.rs` — prose that cites a record and contradicts
//!   its status. It only sees prose that *already* cites one.
//! * `clause_ruling_index.rs` — clause to record, for the reader who never
//!   looked the record up. Its own docs name the limit: "the person who never
//!   looked the record up writes no citation for it to check."
//!
//! Every existing ledger guard — those two plus `claude_md_adjudication_tables`,
//! `normalization_contract_doc`, `ledger_clause_jurisdiction`,
//! `ruling_guard_field` and `normalization_ruleset_page` — fires only on prose
//! that already names a record, or on the ledger's own contents. **The class that
//! cost this project a week names no record and cites no clause**: the falsified
//! equal-length premise was a bare assertion. This file is the only one that can
//! see that shape, which is the whole of its claim to exist, and the reason its
//! key is a phrase rather than a citation.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use crate::common::rulings;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The phrases that mark a comment block as deciding what the correct output
/// **is**, rather than explaining what the code does.
///
/// Matched case-insensitively against the block's text with its comment markers
/// stripped and its line breaks collapsed to single spaces, so a phrase that
/// wraps across two comment lines still matches.
///
/// The selection rule, which is the part worth keeping: a phrase belongs here
/// only if it asserts that a representational question has an answer
/// independent of anyone's choice — "a fact rather than a choice", "the
/// correspondence is unique" — or announces a choice in the first person. Both
/// shapes are claims a ruling record should be carrying. A phrase that merely
/// names an outcome ("the canonical form is `g.262dup`") explains a case and is
/// not an adjudication, which is why the obvious candidates are absent; see the
/// module docs for the measured counts that got them rejected.
const ADJUDICATION_PHRASES: &[&str] = &[
    "a fact rather than a choice",
    "a choice rather than a fact",
    "a fact and not a choice",
    "a choice and not a fact",
    "a fact, not a choice",
    "a choice, not a fact",
    "is a fact about the variant",
    "are a fact about the variant",
    "correspondence is unique",
    "well defined without an alignment",
    "well-defined without an alignment",
    "we choose",
    "we rule that",
    "is the correct output",
    "is the correct form",
    "are the correct output",
];

/// Sites that match [`ADJUDICATION_PHRASES`], name no record, and are
/// deliberately left as they are. `(repo-relative path, why)`.
///
/// **Shrink-only.** A row whose file no longer flags is reported as stale and
/// must be deleted; the list may not grow without the growth being a reviewed
/// line of this file. Seeded on 2026-08-14 with the tree as it stood after the
/// known sites were repointed at their record.
///
/// Keyed by **file**, not by line, so an unrelated edit above a site does not
/// invalidate the row — a line-keyed list would go stale on every rebase and be
/// deleted for being noisy rather than for being wrong.
const REGROWTH_EXEMPT: &[(&str, &str)] = &[(
    "tests/it/issue_1284_transcript_axis_collision.rs",
    "quotes the withdrawn claim verbatim in order to withdraw it — the block \
     reads \"an intermediate revision of this branch read that as …\" and then \
     states the refutation. Recording a withdrawal by reproducing what was \
     withdrawn is the right shape, and it is indistinguishable to a substring \
     scan from asserting it",
)];

/// Every `.rs` file under `src/` and `tests/`, repo-relative.
///
/// `examples/` is deliberately out of scope: a generator's comments describe a
/// measurement, and the two directories this file polices are the ones whose
/// prose a reader reaches for when asking what ferro's output *should* be.
fn scanned_sources() -> Vec<PathBuf> {
    let root = repo_root();
    let mut files = Vec::new();
    for directory in ["src", "tests"] {
        collect_rust_sources(&root.join(directory), &root, &mut files);
    }
    files.sort();
    assert!(
        files.len() > 100,
        "the source scan found only {} files, which cannot be right for this tree — a scan that \
         silently finds nothing passes vacuously",
        files.len()
    );
    files
}

fn collect_rust_sources(directory: &Path, root: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = std::fs::read_dir(directory) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_rust_sources(&path, root, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(
                path.strip_prefix(root)
                    .expect("a scanned path is under the repo root")
                    .to_path_buf(),
            );
        }
    }
}

/// One contiguous run of comment lines: its 1-based first line, and its text
/// with the markers stripped and the line breaks collapsed.
struct CommentBlock {
    line: usize,
    text: String,
}

/// Split a file into its contiguous comment blocks.
///
/// The block — rather than the line — is the unit of judgement, and that is the
/// opposite choice from `ruling_citation_currency.rs`, deliberately. That file
/// reads a *status word* next to a citation, where a two-line window measurably
/// picks up the neighbouring record's word and reports a contradiction that does
/// not exist. Here the question is whether a paragraph that decides something
/// also points at a record, and the citation is almost never on the same
/// physical line as the claim — line scoping would flag every correctly-cited
/// paragraph in the tree.
fn comment_blocks(source: &str) -> Vec<CommentBlock> {
    let lines: Vec<&str> = source.lines().collect();
    let mut blocks = Vec::new();
    let mut index = 0;
    while index < lines.len() {
        if !lines[index].trim_start().starts_with("//") {
            index += 1;
            continue;
        }
        let start = index;
        let mut text = String::new();
        while index < lines.len() && lines[index].trim_start().starts_with("//") {
            let stripped = lines[index]
                .trim_start()
                .trim_start_matches('/')
                .trim_start_matches('!')
                .trim();
            text.push_str(stripped);
            text.push(' ');
            index += 1;
        }
        blocks.push(CommentBlock {
            line: start + 1,
            text,
        });
    }
    blocks
}

/// The phrases a block matches, if any.
fn adjudication_phrases_in(text: &str) -> Vec<&'static str> {
    let lowered = text.to_ascii_lowercase();
    ADJUDICATION_PHRASES
        .iter()
        .filter(|phrase| lowered.contains(**phrase))
        .copied()
        .collect()
}

/// Whether a block points at the ledger.
///
/// Membership, not relevance: the block must contain the id of a record the
/// ledger actually has. A bare `rulings[...]` with a made-up id does not count,
/// which is the same standard `ruling_citation_currency.rs` holds citations to.
fn names_a_record(text: &str, ids: &BTreeSet<String>) -> bool {
    ids.iter().any(|id| text.contains(id.as_str()))
}

/// Every flagged site in the tree: `(path, line, phrases)`.
fn flagged_sites(ids: &BTreeSet<String>) -> Vec<(PathBuf, usize, Vec<&'static str>)> {
    let root = repo_root();
    let mut sites = Vec::new();
    for relative in scanned_sources() {
        // This file states the phrases in order to look for them, so it matches
        // itself on every one. Excluding it by path is the only honest option —
        // any cleverer rule (skip `const` blocks, skip the module docs) would be
        // a hole another file could climb through.
        if relative == Path::new(SELF_RELATIVE_PATH) {
            continue;
        }
        let source = std::fs::read_to_string(root.join(&relative))
            .unwrap_or_else(|e| panic!("read {}: {e}", relative.display()));
        for block in comment_blocks(&source) {
            let phrases = adjudication_phrases_in(&block.text);
            if phrases.is_empty() || names_a_record(&block.text, ids) {
                continue;
            }
            sites.push((relative.clone(), block.line, phrases));
        }
    }
    sites
}

/// This file, so the scan can skip its own statement of the phrases.
const SELF_RELATIVE_PATH: &str = "tests/it/adjudication_prose_regrowth.rs";

// --------------------------------------------------------------------------
// Tests.
// --------------------------------------------------------------------------

/// The ratchet. Adjudication-shaped prose points at a record, or its file is
/// named here with a reason.
#[test]
fn adjudication_shaped_prose_points_at_a_ruling_record_or_is_allowlisted() {
    let exempt: BTreeMap<&str, &str> = REGROWTH_EXEMPT.iter().copied().collect();
    assert_eq!(
        exempt.len(),
        REGROWTH_EXEMPT.len(),
        "REGROWTH_EXEMPT has a duplicate path; one row would shadow another while the count \
         still looked right"
    );
    for (path, reason) in REGROWTH_EXEMPT {
        assert!(
            !reason.trim().is_empty(),
            "REGROWTH_EXEMPT row {path:?} states no reason; an exemption that does not say why \
             is the silence this file exists to make expensive"
        );
        assert!(
            repo_root().join(path).is_file(),
            "REGROWTH_EXEMPT names {path:?}, which is not a file in this tree"
        );
    }

    let ids: BTreeSet<String> = rulings::records().into_iter().map(|r| r.id).collect();
    let sites = flagged_sites(&ids);

    let mut covered: BTreeMap<&str, Vec<String>> = BTreeMap::new();
    let mut undeclared = Vec::new();
    for (path, line, phrases) in &sites {
        let key = path.display().to_string();
        match exempt.get_key_value(key.as_str()) {
            Some((path, _)) => {
                covered
                    .entry(path)
                    .or_default()
                    .push(format!("{key}:{line}  {phrases:?}"));
            }
            None => undeclared.push(format!("{key}:{line}  {phrases:?}")),
        }
    }

    // A row is keyed by FILE, which is the only key that survives editing. The
    // cost of that key is that one row can silently cover several distinct
    // sites — measured elsewhere in this tree as the reason a clause-keyed
    // version of this guard is unusable, since a single row on
    // `src/normalize/merge.rs` would cover ~124 flagged lines. Cap it at one:
    // a second flagged block in an exempt file has not been reviewed by whoever
    // wrote the row, and an exemption that grows silently is the shape this file
    // exists to make expensive.
    let overloaded: Vec<String> = covered
        .iter()
        .filter(|(_, sites)| sites.len() > 1)
        .map(|(path, sites)| {
            format!(
                "{path} covers {} sites:\n    {}",
                sites.len(),
                sites.join("\n    ")
            )
        })
        .collect();
    assert!(
        overloaded.is_empty(),
        "a REGROWTH_EXEMPT row may cover at most one flagged block; these cover more, so the \
         row's stated reason cannot have been written about all of them. Fix the new site, or \
         rewrite the row's reason to cover both and split the phrase list if they differ:\n  {}",
        overloaded.join("\n  ")
    );

    // Shrink-only, reported first so a PR that repoints a site is told to delete
    // its row rather than being told twice about the same file.
    let stale: Vec<&str> = exempt
        .keys()
        .filter(|key| !covered.contains_key(**key))
        .copied()
        .collect();
    assert!(
        stale.is_empty(),
        "REGROWTH_EXEMPT is shrink-only and these rows no longer describe the tree — the site \
         now cites a record, or the prose is gone. Delete them: {stale:?}"
    );

    assert!(
        undeclared.is_empty(),
        "prose here decides what the correct representation is and points at no ruling record.\n\
         \n\
         Record the decision instead of restating it. If the recommendations settle it, add a \
         `rulings` record naming the governing clause; if they are silent, add one with a \
         `house_choice` under `README.md` rule 5's silent limb or rule 6 — see \
         `generate_spec_fixture`'s `overrides::HouseChoice`. Then cite the record's id here \
         rather than repeating its content, because a rule written in several places is how \
         this project's rulings have drifted apart before.\n\
         \n\
         If this is a false positive — the phrase is incidental, or the block quotes a claim in \
         order to withdraw it — add the file to REGROWTH_EXEMPT with the reason.\n\
         \n  {}",
        undeclared.join("\n  ")
    );
}

/// The detector fires on the shape it was built for, and stays silent on the
/// shapes that look like it.
///
/// Pinned against synthetic text rather than against the tree, so it keeps
/// meaning something once the tree's own instances are cleaned up — the same
/// reason `issue_1615_denoted_sequence_oracle.rs` pins recorded triples instead
/// of re-normalizing.
#[test]
fn the_detector_fires_on_an_ontological_claim_and_not_on_an_explained_outcome() {
    let ids: BTreeSet<String> = ["a-test-record".to_string()].into_iter().collect();
    let judge = |text: &str| {
        let blocks = comment_blocks(text);
        assert_eq!(blocks.len(), 1, "the fixture is one block");
        let phrases = adjudication_phrases_in(&blocks[0].text);
        (!phrases.is_empty(), names_a_record(&blocks[0].text, &ids))
    };

    // The disease: an ontological claim about representation, cited to nothing.
    assert_eq!(
        judge("// For an equal-length block the column correspondence is unique,\n// so which columns changed is a fact rather than a choice.\n"),
        (true, false)
    );
    // The same claim, pointed at a record. Flagged and covered.
    //
    // The fixture id is spelled bare rather than in the `rulings[…]` form, and
    // the word that would make it a citation is kept off the line, so that
    // `ruling_citation_currency.rs`'s scan — which reads *this* file like any
    // other — does not report a synthetic id as a citation of a record the
    // ledger does not have. `common/rulings.rs` keeps its own fixture out of
    // that scan's way the same way, and for the same reason.
    assert_eq!(
        judge("// Which columns changed is a fact rather than a choice only\n// under the decided `a-test-record`.\n"),
        (true, true)
    );
    // The phrase wrapping across two comment lines still matches, which is why
    // the block rather than the line is the unit.
    assert_eq!(
        judge("// … so which columns changed is a fact\n// rather than a choice.\n"),
        (true, false)
    );

    // Explaining an outcome is not adjudicating one. These are the shapes that
    // made the broader phrase lists unusable; if any of them starts flagging,
    // the list has been widened past what this file can carry.
    for benign in [
        "// The bug returned the unshifted `g.258dup`; the canonical form is `g.262dup`.\n",
        "// `duplication.md:18` — \"when a variant can be described as a duplication, it must be\n// described as a duplication\" — is a MUST.\n",
        "// `delinsATG` is canonical; must NOT be matched.\n",
        "// Nothing in the spec ranks the two, so no rule is being broken.\n",
        "// `ATG` is the only `Met` codon, so this is not a choice.\n",
    ] {
        assert!(!judge(benign).0, "false positive on: {benign}");
    }
}

/// `//!` and `///` and `//` all read as comments, and a non-comment line ends
/// the block.
///
/// The second half is the load-bearing one: without it a whole file would be one
/// block, and a single citation anywhere in it would launder every claim in it.
#[test]
fn a_comment_block_ends_at_the_first_non_comment_line() {
    let blocks = comment_blocks(
        "//! module doc\nuse std::fmt;\n\n/// item doc\nfn f() {}\n    // indented\n    let x = 1;\n",
    );
    let starts: Vec<usize> = blocks.iter().map(|b| b.line).collect();
    assert_eq!(starts, vec![1, 4, 6]);
    assert_eq!(blocks[0].text.trim(), "module doc");
    assert_eq!(blocks[2].text.trim(), "indented");
}
