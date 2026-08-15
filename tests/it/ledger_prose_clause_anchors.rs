//! A prose citation of a spec clause must name the line that clause is actually
//! on — not merely a line that exists.
//!
//! # The failure this closes
//!
//! The ledger cites spec clauses twice, in two registers. The **structured**
//! `clauses` array carries a `file:line` plus the quoted text, and
//! `generate_spec_fixture`'s `verify_citation` fails the build unless that quote
//! is still on those lines in the pinned spec checkout. The **prose** —
//! `question`, `rationale`, `guard.none`, and the `by_input` notes — cites the
//! same clauses in running text, carries no quote, and until this file was
//! checked by nothing at all.
//!
//! So a submodule bump that moves a clause fails loudly on one register and is
//! silent on the other. That happened: #1793 bumped `6f85311` -> `565b973`,
//! where upstream deleted a duplicated line at old `general.md:22` and every old
//! line from 23 on shifted by one. The structured array was repaired because the
//! build demanded it; the prose was repaired in three records and left in the
//! other seventeen, so **51 prose citations of `general.md:34` went on naming a
//! line that now holds the codon *exception* while meaning the separation rule
//! one line above it** (#2009).
//!
//! (#2009 states 45. That was measured at `8881a54e`; the ledger grew before this
//! branch, and at this commit's own base the count is 51. Every figure below is
//! measured at this commit, and each says which ledger it was taken over.)
//!
//! That is worse than an off-by-one, and it is why this is worth a guard rather
//! than a one-off sweep. The two clauses argue opposite ways: `:33` splits
//! ("described individually and **not** as a delins") and `:34` merges
//! ("**exception**: … should be described as a delins"). A reader following a
//! citation landed on the clause that says the reverse of what the sentence
//! citing it claimed.
//!
//! # Why a line-keyed table cannot do this, and what does
//!
//! The obvious guard — "every prose citation must name a line the structured
//! array also names" — does not fire on this defect, and it is worth recording
//! why so it is not proposed again. `general.md:34` *is* a structured citation:
//! two records cite it legitimately, meaning the exception. A rule keyed on
//! `file:line` cannot tell a prose `:34` that means the exception from one that
//! means the separation rule, because a line number carries no meaning. Only a
//! **phrase** does.
//!
//! So the check reads the meaning off the citing sentence. [`ANCHORS`] names a
//! handful of clause identities. Each carries a **signature** — a phrase used to
//! *locate* the clause in the spec checkout, so the expected line number is
//! derived from the spec on every run and never restated here — and a set of
//! **triggers**, phrases the ledger's own prose uses when it is talking about
//! that clause. For each prose citation of `general.md`, the sentence containing
//! it is examined; if exactly one identity triggers, the citation must name that
//! identity's located line.
//!
//! Measured: over the ledger at the commit before this file, that fires on **27**
//! citations across 20 record/field pairs. Over the repaired ledger it fires on
//! none, checking 29.
//!
//! # How much this actually covers: 29 of 117, and 10 of the 53 that mattered
//!
//! **This guard is a floor, not a proof, and the fraction is small enough that it
//! must be stated as a number rather than as a caveat.** Measured over the
//! repaired ledger at this commit, by the algorithm below:
//!
//! | population | n |
//! |---|---:|
//! | qualified `general.md:N` prose citations | **117** |
//! | of those, **judged** by the trigger rule | **29** (25%) |
//! | skipped — the sentence carries no trigger | **84** |
//! | skipped — the sentence carries two citations (2 sentences) | **4** |
//! | citations of `general.md:33`, the separation rule #2009 was about | **53** |
//! | of those 53, judged | **10** |
//! | of those 53, in zero-trigger sentences and so unreachable | **43** |
//! | citations this PR repaired in total (100 qualified + 81 bare) | **181** |
//! | of those 181, judged here | **29** (16%) |
//!
//! So **the defect this file exists for can be reintroduced in 43 of the 53
//! places it lived, without reddening anything.** That was verified, not
//! inferred: reverting `general.md:33` -> `:34` in the zero-trigger sentence
//! "Probes 1-4 show it splitting at separation ONE and TWO on a genomic axis,
//! which is what `general.md:33` asks for."
//! (`adjudication-precedence-order.rationale`) leaves all five checks green.
//! Reverting the *same* citation in a trigger-bearing sentence fires. The
//! difference is the sentence, not the citation.
//!
//! **Do not read the 29 as coverage.** Read it as: the sentences that say what
//! clause they mean are checked, and most sentences do not say.
//!
//! ## Widening the triggers was measured and is NOT the fix
//!
//! The obvious next move is to widen `separation-rule`'s triggers until the 43
//! are reachable. Measured, one candidate phrase at a time, over the repaired
//! ledger — where every citation is known correct, so **any fire is a false
//! failure**:
//!
//! | extra trigger | judged | false failures |
//! |---|---:|---:|
//! | `"separation"` | 46 | **1** |
//! | `"separation rule"` | 31 | **1** |
//! | `"separated by"` / `"two variants separated"` | 29 / 28 | **1** |
//! | `"governs"` | 31 | 0 |
//! | `"asks for"`, `"unchanged nucleotide"` | 30 | 0 |
//! | `"not as a delins"`, `"splits"`, `"separation of two or more"` | 29 | 0 |
//!
//! The only candidates that reach a useful number of the 43 are the ones that
//! fire falsely, and the ones that fire cleanly buy one or two citations while
//! being generic words (`"governs"` is clean here only because no sentence
//! happens to say "`general.md:34` governs" today — that is a fact about the
//! corpus, not about the trigger). The false failure `"separation"` produces is
//! instructive: it is on a sentence *about the bump itself* ("the separation rule
//! moved 34 -> 33") which cites `general.md:22`. A phrase that names a clause and
//! a phrase that discusses one are not distinguishable by substring.
//!
//! The same limit already has a live instance in the shipped table, worth knowing
//! before it lands on someone: `:55`, `:56` and `:57` are one bullet and its two
//! sub-bullets, **all three about prioritisation**, so inserting the word
//! "prioritisation" beside the already-correct `general.md:56` in
//! `duplication-must-ranks-the-label-not-the-partition` makes this guard demand
//! `:55` — a false failure on correct prose. The failure text names the right
//! escape ("the trigger is too broad"), and the trigger stays narrow for exactly
//! this reason: a trigger that never fires costs a skipped citation, one that
//! fires wrongly costs a correct citation being "repaired" into a wrong one.
//!
//! # What it does NOT catch, stated rather than left to be discovered
//!
//! - **Citations whose sentence names no trigger.** Most of the ledger's
//!   citations are bare cross-references ("`general.md:33` still governs"), and
//!   for those the prose says nothing that could identify the clause. They are
//!   counted as skipped and checked only by [`every_prose_citation_resolves`]'s
//!   weaker in-range / non-blank rule below. On the repaired ledger **84 of the
//!   117** are in this class.
//! - **Sentences carrying two or more `general.md` citations.** A trigger cannot
//!   be attributed to one of two citations, so both are skipped. **Two such
//!   sentences on the repaired ledger, carrying 4 citations between them.**
//! - **Sentences that trigger two identities.** Same reason, skipped. **None on
//!   the repaired ledger**, so this arm is currently untested by the corpus.
//! - **Bare `` `:33` `` shorthands.** The ledger cites a clause by bare line
//!   number once its file is established by an earlier sentence. Those carry no
//!   file, so nothing here can resolve them, and #2009 had to repair **81** of
//!   them by hand. This is the largest hole — it is structural rather than a
//!   matter of trigger breadth, since the scanner keys on the literal
//!   `general.md:` — and closing it needs an attribution convention in the prose,
//!   not a cleverer matcher.
//! - **Files other than `general.md`.** `general.md` is the only spec file whose
//!   line numbers moved in #1793, and it is the only one whose basename resolves
//!   unambiguously inside `docs/` — prose writes `delins.md:17` for a file that
//!   exists three times (`DNA/`, `RNA/`, `protein/`). Widening this is a
//!   reasonable follow-up; guessing which `delins.md` a sentence means is not.
//! - **Whether a record's argument is *right*.** This checks that a citation
//!   points at the clause its own sentence is about. It has nothing to say about
//!   whether that clause supports the claim.
//!
//! # Why the expected line is located rather than written down
//!
//! Restating `33` here would make this a change detector for the number, in the
//! exact way the repository `CLAUDE.md` warns about under "Assert the property.
//! Measure the count." — and it would be the *same* number that just went stale,
//! kept in a second place. The signature phrases are located in the spec on every
//! run, so the next bump moves both the expectation and the check together, and a
//! signature that stops resolving fails
//! [`every_anchor_signature_locates_exactly_one_line`] by name instead of
//! silently disarming everything below it.

use std::collections::BTreeSet;
use std::path::PathBuf;

use crate::common::rulings::{self, Guard, LEDGER_RELATIVE_PATH};

/// The vendored spec checkout, relative to the crate root — the same default
/// `generate_spec_fixture` takes for `--spec-dir`.
const SPEC_DIR: &str = "assets/hgvs-nomenclature";

/// The one spec file whose lines moved in #1793, and the only one whose
/// basename is unambiguous inside `docs/`.
const SPEC_FILE: &str = "docs/recommendations/general.md";

/// A clause the ledger argues over, in three parts.
struct Anchor {
    /// A name for the clause, used only in failure messages.
    identity: &'static str,
    /// A phrase that occurs on exactly one line of [`SPEC_FILE`]. The line it is
    /// found on IS the expectation; nothing here restates a line number.
    signature: &'static str,
    /// Phrases the ledger's prose uses when the sentence is about this clause.
    /// Deliberately narrow: a trigger that fires on ordinary prose costs a false
    /// failure, while one that never fires costs only a skipped citation.
    triggers: &'static [&'static str],
}

/// The clauses `general.md` prose actually litigates.
///
/// Six, not an attempt at completeness. Each earns its place by being cited in
/// prose several times *and* by having a phrase the citing sentences reliably
/// use. `separation-rule` and `codon-exception` are the pair #2009 is about —
/// adjacent lines that argue opposite ways — and the rest are the clauses whose
/// prose citations the same bump moved.
const ANCHORS: &[Anchor] = &[
    Anchor {
        identity: "separation-rule",
        signature: r#"should be described individually and **not** as a "delins""#,
        triggers: &[
            "described individually",
            "separated by one or more nucleotide",
        ],
    },
    Anchor {
        identity: "codon-exception",
        signature: "**exception**: two variants separated by one nucleotide, together affecting \
                    one amino acid",
        triggers: &[
            "codon exception",
            "together affecting one amino acid",
            "separated by one nucleotide",
        ],
    },
    Anchor {
        identity: "prioritisation",
        signature: "the preferred description is: (1) substitution, (2) deletion, (3) inversion",
        triggers: &[
            "prioritisation",
            "ranks substitution",
            "ranks single-variant",
            "(1) substitution",
            "ranks the type of one description",
        ],
    },
    Anchor {
        identity: "self-replacement-prohibition",
        signature: "descriptions removing part of a reference sequence and replacing it with part \
                    of the same sequence are not allowed",
        triggers: &["removing part of a reference sequence"],
    },
    Anchor {
        identity: "junction-exception",
        signature: "**exception**: deletions/duplications around exon/exon junctions",
        triggers: &["exon/exon junction"],
    },
    Anchor {
        identity: "three-prime-universality",
        signature: "the 3'rule applies to ALL descriptions",
        triggers: &["to all descriptions (genome"],
    },
];

/// Floor for how many citations the trigger rule must actually judge.
///
/// Measured at **29** on the commit that added this file. The floor sits at
/// less than half of that, so ordinary prose edits cannot trip it while a
/// matcher that stopped recognising sentences, citations or triggers still
/// would — the same reasoning `ruling_citation_currency.rs`'s
/// `STATUS_CLAIMS_FLOOR` gives.
const ANCHORED_CITATIONS_FLOOR: usize = 12;

/// Floor for how many prose citations exist at all, for the weaker check.
///
/// Measured at **117** on the commit that added this file. A run that resolved
/// none of them would pass every assertion below by looking at nothing.
const PROSE_CITATIONS_FLOOR: usize = 40;

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The spec file's lines, 0-indexed, as read from the pinned checkout.
fn spec_lines() -> Vec<String> {
    let path = crate_root().join(SPEC_DIR).join(SPEC_FILE);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "read {}: {e} — the spec submodule is probably not initialised. Run\n    \
             git -c protocol.file.allow=always submodule update --init {SPEC_DIR}",
            path.display()
        )
    });
    text.lines().map(str::to_string).collect()
}

/// Whitespace-collapsed, lower-cased, for substring comparison.
///
/// The same normalisation `generate_spec_fixture::verify_citation` applies, plus
/// case folding: the ledger's prose shouts (`DESCRIBED INDIVIDUALLY`) where the
/// spec does not.
fn collapse(text: &str) -> String {
    text.split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .to_lowercase()
}

/// The 1-based lines of [`SPEC_FILE`] on which `signature` occurs.
fn lines_matching(lines: &[String], signature: &str) -> Vec<usize> {
    let needle = collapse(signature);
    lines
        .iter()
        .enumerate()
        .filter(|(_, line)| collapse(line).contains(&needle))
        .map(|(index, _)| index + 1)
        .collect()
}

/// One prose field of the ledger, with enough provenance to name it in a
/// failure.
struct ProseField {
    owner: String,
    field: &'static str,
    text: String,
}

/// Every prose field of the ledger that can carry a clause citation.
///
/// The records come through `common::rulings`, so "what a record is" keeps its
/// single definition. `by_input` has no such reader and is parsed here — it is
/// a flat map of notes, and the one note that cites `general.md` cites it the
/// same way a rationale does.
fn prose_fields() -> Vec<ProseField> {
    let mut fields = Vec::new();
    for record in rulings::records() {
        fields.push(ProseField {
            owner: record.id.clone(),
            field: "question",
            text: record.question.clone(),
        });
        fields.push(ProseField {
            owner: record.id.clone(),
            field: "rationale",
            text: record.rationale.clone(),
        });
        if let Guard::Declined(reason) = &record.guard {
            fields.push(ProseField {
                owner: record.id.clone(),
                field: "guard.none",
                text: reason.clone(),
            });
        }
    }

    let path = rulings::ledger_path();
    let text =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let document: serde_json::Value =
        serde_json::from_str(&text).unwrap_or_else(|e| panic!("parse {}: {e}", path.display()));
    let by_input = document
        .get("by_input")
        .and_then(|v| v.as_object())
        .unwrap_or_else(|| panic!("{LEDGER_RELATIVE_PATH} has no `by_input` object"));
    for (input, entry) in by_input {
        if let Some(note) = entry.get("note").and_then(|v| v.as_str()) {
            fields.push(ProseField {
                owner: input.clone(),
                field: "by_input.note",
                text: note.to_string(),
            });
        }
    }
    fields
}

/// A `general.md:<line>` or `general.md:<first>-<last>` citation found in prose.
struct Citation {
    first: usize,
    last: usize,
    /// The citation as written, for the failure message.
    spelling: String,
}

/// Every `general.md:<line>` citation in `text`, in order.
///
/// Written by hand rather than with a regex crate because the grammar is three
/// tokens wide and the test tree carries no regex dependency.
fn citations_in(text: &str) -> Vec<Citation> {
    const MARKER: &str = "general.md:";
    let bytes = text.as_bytes();
    let mut found = Vec::new();
    let mut from = 0usize;
    while let Some(offset) = text[from..].find(MARKER) {
        let start = from + offset + MARKER.len();
        let digits: String = text[start..]
            .chars()
            .take_while(char::is_ascii_digit)
            .collect();
        from = start.max(from + offset + 1);
        if digits.is_empty() {
            continue;
        }
        let after = start + digits.len();
        let mut last = digits.clone();
        if bytes.get(after) == Some(&b'-') {
            let tail: String = text[after + 1..]
                .chars()
                .take_while(char::is_ascii_digit)
                .collect();
            if !tail.is_empty() {
                last = tail;
                from = after + 1 + last.len();
            }
        } else {
            from = after;
        }
        let first: usize = digits.parse().expect("digits parse");
        let last: usize = last.parse().expect("digits parse");
        found.push(Citation {
            spelling: if first == last {
                format!("{MARKER}{first}")
            } else {
                format!("{MARKER}{first}-{last}")
            },
            first,
            last,
        });
    }
    found
}

/// Split prose into sentences.
///
/// Deliberately crude — a break after `.`, `?` or `!` followed by whitespace,
/// and at every newline. Over-splitting only shrinks the window a trigger has to
/// appear in, which costs skipped citations rather than wrong verdicts;
/// under-splitting is what lets a neighbouring clause's trigger be attributed to
/// the wrong citation, which costs a false failure. When the two error modes are
/// not symmetric, take the one that fails safe.
fn sentences(text: &str) -> Vec<&str> {
    let bytes = text.as_bytes();
    let mut out = Vec::new();
    let mut start = 0usize;
    for (index, byte) in bytes.iter().enumerate() {
        let terminator = matches!(byte, b'.' | b'?' | b'!')
            && bytes
                .get(index + 1)
                .is_some_and(|next| next.is_ascii_whitespace());
        let newline = *byte == b'\n';
        if terminator || newline {
            let end = if newline { index } else { index + 1 };
            if text.is_char_boundary(start) && text.is_char_boundary(end) && start < end {
                out.push(&text[start..end]);
            }
            start = end;
        }
    }
    if start < text.len() && text.is_char_boundary(start) {
        out.push(&text[start..]);
    }
    out
}

/// Every anchor's signature must locate exactly one line of the spec file.
///
/// This is the load-bearing precondition for everything below: the expectation
/// is *derived* from the spec, so a signature that resolves to zero lines (the
/// clause was reworded) or to several (the clause was duplicated — which is
/// precisely what old `general.md:17`/`:22` were before this bump) would
/// silently disarm the check it feeds. It fails here, by name, instead.
#[test]
fn every_anchor_signature_locates_exactly_one_line() {
    let lines = spec_lines();
    assert!(
        lines.len() > 100,
        "{SPEC_FILE} has only {} lines — the spec checkout is not what this expects",
        lines.len()
    );
    let mut located = Vec::new();
    for anchor in ANCHORS {
        let hits = lines_matching(&lines, anchor.signature);
        assert_eq!(
            hits.len(),
            1,
            "the `{}` anchor's signature matches {} lines of {SPEC_FILE} ({hits:?}); it must \
             match exactly one, since the line it matches IS the expectation. Re-point the \
             signature at the clause's current wording.",
            anchor.identity,
            hits.len()
        );
        located.push((anchor.identity, hits[0]));
    }
    let distinct: BTreeSet<usize> = located.iter().map(|(_, line)| *line).collect();
    assert_eq!(
        distinct.len(),
        located.len(),
        "two anchors located the same line of {SPEC_FILE}: {located:?} — one of the signatures \
         is not distinctive enough to identify its clause"
    );
}

/// A citation's triggers must not overlap another anchor's, or a sentence about
/// one clause would be judged against another.
///
/// Checked against the anchors themselves rather than against the ledger: this
/// is a property of the table, and finding it via a failing ledger citation
/// would report the drift at the wrong place.
#[test]
fn no_trigger_belongs_to_two_anchors() {
    // Identities are the key the self-comparison below skips on, so they must be
    // distinct or an anchor would be compared against a namesake and excused.
    // `std::ptr::eq` would be the obvious skip and is wrong here: `ANCHORS` is a
    // `const`, so it is inlined at each use site and the two loops are not
    // guaranteed to walk one allocation. It passes today only because rustc
    // deduplicates them.
    let identities: BTreeSet<&str> = ANCHORS.iter().map(|anchor| anchor.identity).collect();
    assert_eq!(
        identities.len(),
        ANCHORS.len(),
        "two anchors share an identity; the self-comparison below skips on it, so a duplicate \
         would excuse a real collision"
    );

    let mut collisions = Vec::new();
    for anchor in ANCHORS {
        for trigger in anchor.triggers {
            assert!(
                trigger.chars().all(|c| !c.is_ascii_uppercase()),
                "the `{}` anchor's trigger {trigger:?} contains upper case; triggers are matched \
                 against lower-cased prose and would never fire",
                anchor.identity
            );
            for other in ANCHORS {
                if other.identity == anchor.identity {
                    continue;
                }
                if other.triggers.iter().any(|t| t.contains(trigger)) {
                    collisions.push(format!(
                        "`{}`'s trigger {trigger:?} is contained in a trigger of `{}`",
                        anchor.identity, other.identity
                    ));
                }
            }
        }
    }
    assert!(
        collisions.is_empty(),
        "triggers must identify one clause each:\n  {}",
        collisions.join("\n  ")
    );
}

/// Every prose citation of `general.md` must name a line that exists and holds
/// text.
///
/// The weak, universal half — it judges the citations the trigger rule below
/// skips. It is genuinely weak: it would not have caught #2009's central defect,
/// because `general.md:34` is a real, non-blank line. What it does catch is a
/// citation that has fallen off the end of the file or onto one of the blank
/// lines that separate the recommendation bullets, which is what a citation
/// looks like once a clause has been *deleted* rather than moved. Over the
/// pre-repair ledger it fired on `general.md:96`, whose target is blank after
/// the bump.
#[test]
fn every_prose_citation_resolves() {
    let lines = spec_lines();
    let mut faults = Vec::new();
    let mut seen = 0usize;
    for field in prose_fields() {
        for citation in citations_in(&field.text) {
            seen += 1;
            if citation.last > lines.len() {
                faults.push(format!(
                    "{}.{} cites {}, past the end of {SPEC_FILE} ({} lines)",
                    field.owner,
                    field.field,
                    citation.spelling,
                    lines.len()
                ));
                continue;
            }
            if lines[citation.first - 1..citation.last]
                .iter()
                .all(|line| line.trim().is_empty())
            {
                faults.push(format!(
                    "{}.{} cites {}, which is blank in {SPEC_FILE}",
                    field.owner, field.field, citation.spelling
                ));
            }
        }
    }
    assert!(
        faults.is_empty(),
        "these prose citations in {LEDGER_RELATIVE_PATH} do not resolve against the pinned spec \
         checkout. A submodule bump moves prose citations exactly as it moves the structured \
         `clauses` array, and only the latter fails the build on its own:\n  {}",
        faults.join("\n  ")
    );
    assert!(
        seen >= PROSE_CITATIONS_FLOOR,
        "only {seen} prose citations of {SPEC_FILE} were found in {LEDGER_RELATIVE_PATH}, below \
         the floor of {PROSE_CITATIONS_FLOOR} — `citations_in` is probably broken, and this \
         check would pass by reading nothing"
    );
}

/// A prose citation whose sentence identifies the clause must name that clause's
/// line.
///
/// This is the check #2009 exists for. See the module docs for the reasoning and
/// for the four classes it deliberately skips.
#[test]
fn a_citation_names_the_clause_its_sentence_is_about() {
    let lines = spec_lines();
    let located: Vec<(&Anchor, usize)> = ANCHORS
        .iter()
        .map(|anchor| {
            let hits = lines_matching(&lines, anchor.signature);
            assert_eq!(
                hits.len(),
                1,
                "the `{}` anchor's signature does not locate one line — see \
                 `every_anchor_signature_locates_exactly_one_line`, which reports this properly",
                anchor.identity
            );
            (anchor, hits[0])
        })
        .collect();

    let mut faults = Vec::new();
    let mut checked = 0usize;
    for field in prose_fields() {
        for sentence in sentences(&field.text) {
            let citations = citations_in(sentence);
            if citations.is_empty() {
                continue;
            }
            let lowered = sentence.to_lowercase();
            let triggered: Vec<&(&Anchor, usize)> = located
                .iter()
                .filter(|(anchor, _)| anchor.triggers.iter().any(|t| lowered.contains(t)))
                .collect();
            // Zero triggers: the sentence says nothing that identifies a clause.
            // Two or more: it talks about two, and nothing attributes either to
            // this citation. Two or more citations: the trigger cannot be
            // attributed to one of them. All three are skips, not verdicts.
            if triggered.len() != 1 || citations.len() != 1 {
                continue;
            }
            let (anchor, expected) = triggered[0];
            let citation = &citations[0];
            checked += 1;
            if citation.first <= *expected && *expected <= citation.last {
                continue;
            }
            faults.push(format!(
                "{}.{} cites {} in a sentence about the `{}` clause, which is at \
                 {SPEC_FILE}:{expected}:\n      {}",
                field.owner,
                field.field,
                citation.spelling,
                anchor.identity,
                sentence.trim()
            ));
        }
    }

    assert!(
        faults.is_empty(),
        "{} prose citation(s) in {LEDGER_RELATIVE_PATH} name a different line from the clause \
         their own sentence is about:\n  {}\n\n\
         The line number is what moves when the spec submodule is bumped; the sentence is not. \
         Re-point the citation at the line the clause is on now — do NOT re-word the sentence to \
         match a stale number, and do NOT relax the anchor. If the sentence really is about a \
         clause other than the one that triggered, the trigger is too broad and belongs in \
         `ANCHORS` with a narrower phrase.",
        faults.len(),
        faults.join("\n  ")
    );

    assert!(
        checked >= ANCHORED_CITATIONS_FLOOR,
        "only {checked} citations were judged, below the floor of {ANCHORED_CITATIONS_FLOOR} — \
         the sentence splitter, the citation scanner or the triggers stopped matching, and this \
         check would pass by judging nothing"
    );
}

/// The two scanners, pinned against literals.
///
/// Both are hand-rolled, and both have a failure mode that looks like success:
/// a citation scanner that finds nothing and a sentence splitter that returns
/// one giant sentence each leave the checks above passing while judging much
/// less than they claim.
#[test]
fn the_scanners_read_citations_and_sentence_boundaries() {
    let found = citations_in(
        "`general.md:33` governs, `general.md:35-38` is dead text, and `delins.md:17` agrees; \
         `general.md:` names no line.",
    );
    let spellings: Vec<&str> = found.iter().map(|c| c.spelling.as_str()).collect();
    assert_eq!(
        spellings,
        ["general.md:33", "general.md:35-38"],
        "the scanner must read the single-line and range forms and skip a bare `general.md:`"
    );
    assert_eq!((found[1].first, found[1].last), (35, 38));

    let split = sentences("One. Two? Three!\nFour — with `c.76_83inv` and e.g. a period.");
    assert!(
        split.len() >= 5,
        "the splitter returned {split:?}; a splitter that returns one sentence would let a \
         neighbouring clause's trigger be attributed to the wrong citation"
    );
    assert!(
        split.iter().any(|s| s.contains("c.76_83inv")),
        "an HGVS description's internal period must not lose the text around it: {split:?}"
    );
}
