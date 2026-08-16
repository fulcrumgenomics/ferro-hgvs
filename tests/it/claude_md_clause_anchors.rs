//! A `general.md:N` clause citation in the repository `CLAUDE.md` must name the
//! line that clause is actually on in the pinned spec checkout — not merely a
//! line that exists.
//!
//! # The failure this closes
//!
//! This is the `CLAUDE.md` sibling of [`ledger_prose_clause_anchors`]. The
//! adjudication ledger cites spec clauses in two registers, and both are guarded:
//! the structured `clauses` array by `generate_spec_fixture`'s `verify_citation`,
//! and the ledger's own prose by [`ledger_prose_clause_anchors`]. `CLAUDE.md`
//! cites the same clauses — in its prose and, densely, in the decided/open
//! ruling tables — and until this file was guarded by nothing.
//!
//! So the #1793 submodule bump (`6f85311` -> `565b973`), which deleted a
//! duplicated line at old `general.md:22` and shifted every old line from 23 on
//! down by one, moved every clause `CLAUDE.md` cites while leaving the numbers
//! written beside them untouched. #2009/#2014 swept the ledger; #2057 is the same
//! defect in `CLAUDE.md`, which that sweep did not reach. The sharpest instance:
//! the two clauses argue opposite ways — `:33` splits ("described individually
//! and **not** as a delins") and `:34` merges ("**exception**: … should be
//! described as a delins") — so a reader following a stale `general.md:34`
//! citation of the separation rule lands on the clause that says the reverse.
//!
//! # Why a line-keyed table cannot do this, and what does
//!
//! The obvious guard — "every citation must name a line the ledger's structured
//! array also names" — cannot fire on this defect. `general.md:34` *is* a real,
//! non-blank, legitimately-cited line: it is the codon exception, and `CLAUDE.md`
//! cites it correctly for that. A line number carries no meaning, so only a
//! **phrase** can tell a `:34` that means the exception from a stale `:34` that
//! means the separation rule one line above it.
//!
//! So the check reads the meaning off the citing sentence, exactly as
//! [`ledger_prose_clause_anchors`] does. [`ANCHORS`] names the clauses `CLAUDE.md`
//! litigates. Each carries a **signature** — a phrase found on exactly one line
//! of the spec checkout, so the expected line number is *derived from the spec on
//! every run* and never restated here — and a set of **triggers**, phrases
//! `CLAUDE.md`'s own sentences use when they are about that clause. For each
//! `general.md:N` citation, the sentence containing it is examined; if exactly one
//! identity triggers and the sentence carries exactly one citation, that citation
//! must name (or, for a range, contain) that identity's located line.
//!
//! # What it does NOT catch, stated rather than left to be discovered
//!
//! - **Bare `` `:N` `` shorthands.** `CLAUDE.md` cites a clause by bare line
//!   number once its file is established by an earlier `general.md:` token in the
//!   same row. The scanner keys on the literal `general.md:`, so those are not
//!   seen — the same structural hole [`ledger_prose_clause_anchors`] documents,
//!   and #2015's territory across the wider corpus. The #2057 fix re-points the
//!   bare shorthands that share a row with a qualified anchor by hand, so no
//!   edited row is internally contradictory; nothing here enforces that.
//! - **Files other than `general.md`.** It is the only spec file whose lines
//!   moved in #1793 and the only one whose basename resolves unambiguously inside
//!   `docs/` (`delins.md` exists under `DNA/`, `RNA/` and `protein/`).
//! - **Citations whose sentence names no trigger.** Most citations in the tables
//!   are bare cross-references; their sentence says nothing that identifies the
//!   clause, so they are skipped and covered only by the weaker
//!   [`every_qualified_citation_resolves`] rule.
//! - **Whether a record's argument is *right*.** This checks that a citation
//!   points at the clause its own sentence is about, nothing more.
//!
//! # Why the expected line is located rather than written down
//!
//! Restating `33` here would make this a change detector for the number, in the
//! exact way `CLAUDE.md` warns about under "Assert the property, measure the
//! count" — and it would be the same number that just went stale, kept in a
//! second place. The signatures are located in the spec on every run, so the next
//! bump moves both the expectation and the check together, and a signature that
//! stops resolving fails [`every_anchor_signature_locates_exactly_one_line`] by
//! name instead of silently disarming everything below it.

use std::collections::BTreeSet;
use std::path::PathBuf;

/// The vendored spec checkout, relative to the crate root — the same default
/// `generate_spec_fixture` takes for `--spec-dir`.
const SPEC_DIR: &str = "assets/hgvs-nomenclature";

/// The one spec file whose lines moved in #1793, and the only one whose basename
/// is unambiguous inside `docs/`.
const SPEC_FILE: &str = "docs/recommendations/general.md";

/// The document this guard is about, relative to the crate root.
const CLAUDE_MD: &str = "CLAUDE.md";

/// A clause `CLAUDE.md` argues over, in three parts.
struct Anchor {
    /// A name for the clause, used only in failure messages.
    identity: &'static str,
    /// A phrase that occurs on exactly one line of [`SPEC_FILE`]. The line it is
    /// found on IS the expectation; nothing here restates a line number.
    signature: &'static str,
    /// Phrases `CLAUDE.md`'s prose uses when the sentence is about this clause.
    /// Deliberately narrow: a trigger that fires on ordinary prose costs a false
    /// failure, while one that never fires costs only a skipped citation.
    triggers: &'static [&'static str],
}

/// The clauses `CLAUDE.md`'s prose actually litigates.
///
/// `separation-rule` and `codon-exception` are the adjacent pair #2057 is about;
/// the rest are the other `general.md` clauses whose `CLAUDE.md` citations the
/// same bump moved. Each earns its place by being cited in a sentence that names
/// it with a phrase no other clause's sentence uses.
const ANCHORS: &[Anchor] = &[
    Anchor {
        identity: "separation-rule",
        signature: r#"should be described individually and **not** as a "delins""#,
        triggers: &["separation of two or more", "the excluded alternative"],
    },
    Anchor {
        identity: "codon-exception",
        signature: "**exception**: two variants separated by one nucleotide, together affecting \
                    one amino acid",
        triggers: &[
            "together affecting one amino acid",
            "codon exception",
            "separated by one nucleotide",
        ],
    },
    Anchor {
        identity: "prioritisation",
        signature: "the preferred description is: (1) substitution, (2) deletion, (3) inversion",
        triggers: &["ranks substitution", "above inversion"],
    },
    Anchor {
        identity: "self-replacement-prohibition",
        signature: "descriptions removing part of a reference sequence and replacing it with part \
                    of the same sequence are not allowed",
        triggers: &["are not allowed"],
    },
    Anchor {
        identity: "three-prime-universality",
        signature: "the 3'rule applies to ALL descriptions",
        triggers: &["one variant's several descriptions", "all descriptions"],
    },
    Anchor {
        identity: "dna-level-symbols",
        signature: "nucleotides in CAPITALS using",
        triggers: &["admits only iupac-iubmb symbols"],
    },
    Anchor {
        identity: "rna-level-symbols",
        signature: "nucleotides in lower case using",
        triggers: &["description's symbols"],
    },
    Anchor {
        identity: "protein-re-derivation",
        signature: "protein variant descriptions should be derived from comparing",
        triggers: &["strongest statement of the method"],
    },
    Anchor {
        identity: "svd-wg-forward-note",
        signature: "the SVD-WG is preparing a proposal to modify this recommendation",
        triggers: &["reads as forward guidance"],
    },
];

/// Floor for how many citations the trigger rule must actually judge.
///
/// Measured at 12 on the commit that added this file. The floor sits below half
/// of that, so ordinary prose edits cannot trip it while a matcher that stopped
/// recognising sentences, citations or triggers still would — the same reasoning
/// [`ledger_prose_clause_anchors`]'s `ANCHORED_CITATIONS_FLOOR` gives.
const ANCHORED_CITATIONS_FLOOR: usize = 6;

/// Floor for how many `general.md:N` citations exist at all, for the weak check.
///
/// Measured at 32 on the commit that added this file. A run that resolved none of
/// them would pass every assertion below by looking at nothing.
const QUALIFIED_CITATIONS_FLOOR: usize = 20;

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The spec file's lines, as read from the pinned checkout.
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

/// The whole of [`CLAUDE_MD`].
fn claude_md() -> String {
    let path = crate_root().join(CLAUDE_MD);
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// Whitespace-collapsed, lower-cased, for substring comparison.
///
/// The same normalisation `generate_spec_fixture::verify_citation` applies, plus
/// case folding: `CLAUDE.md`'s prose shouts (`ALL descriptions`) where the spec
/// does not.
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

/// A `general.md:<line>` or `general.md:<first>-<last>` citation found in prose.
struct Citation {
    first: usize,
    last: usize,
    /// The citation as written, for the failure message.
    spelling: String,
}

/// Every `general.md:<line>` citation in `text`, in order.
///
/// Hand-rolled rather than a regex because the grammar is three tokens wide and
/// the test tree carries no regex dependency — the same scanner
/// [`ledger_prose_clause_anchors`] uses, kept self-contained so this guard and
/// its floors do not perturb that file's measured constants.
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
/// Deliberately crude — a break after `.`, `?` or `!` followed by whitespace, and
/// at every newline. Over-splitting only shrinks the window a trigger has to
/// appear in, which costs skipped citations rather than wrong verdicts;
/// under-splitting lets a neighbouring clause's trigger be attributed to the
/// wrong citation, which costs a false failure. When the error modes are not
/// symmetric, take the one that fails safe. In `CLAUDE.md`'s single-line table
/// rows the newline break is inert and `. ` does the work.
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
/// The load-bearing precondition for everything below: the expectation is
/// *derived* from the spec, so a signature that resolves to zero lines (the
/// clause was reworded) or to several (the clause was duplicated) would silently
/// disarm the check it feeds. It fails here, by name, instead.
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
        "two anchors located the same line of {SPEC_FILE}: {located:?} — one of the signatures is \
         not distinctive enough to identify its clause"
    );
}

/// A citation's triggers must not overlap another anchor's, or a sentence about
/// one clause would be judged against another.
#[test]
fn no_trigger_belongs_to_two_anchors() {
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

/// Every `general.md:N` citation in `CLAUDE.md` must name a line that exists and
/// holds text.
///
/// The weak, universal half — it judges the citations the trigger rule below
/// skips. It would not have caught #2057's central defect (`general.md:34` is a
/// real, non-blank line), but it catches a citation that has fallen off the end
/// of the file or onto one of the blank lines separating the recommendation
/// bullets, which is what a citation looks like once a clause has been *deleted*.
#[test]
fn every_qualified_citation_resolves() {
    let lines = spec_lines();
    let text = claude_md();
    let mut faults = Vec::new();
    let mut seen = 0usize;
    for sentence in sentences(&text) {
        for citation in citations_in(sentence) {
            seen += 1;
            if citation.last > lines.len() {
                faults.push(format!(
                    "{} cites {}, past the end of {SPEC_FILE} ({} lines)",
                    CLAUDE_MD,
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
                    "{} cites {}, which is blank in {SPEC_FILE}",
                    CLAUDE_MD, citation.spelling
                ));
            }
        }
    }
    assert!(
        faults.is_empty(),
        "these `general.md` citations in {CLAUDE_MD} do not resolve against the pinned spec \
         checkout. A submodule bump moves prose citations exactly as it moves the ledger's \
         structured `clauses` array, and only the latter fails the build on its own:\n  {}",
        faults.join("\n  ")
    );
    assert!(
        seen >= QUALIFIED_CITATIONS_FLOOR,
        "only {seen} `general.md:N` citations were found in {CLAUDE_MD}, below the floor of \
         {QUALIFIED_CITATIONS_FLOOR} — `citations_in` is probably broken, and this check would \
         pass by reading nothing"
    );
}

/// A `general.md` citation whose sentence identifies the clause must name that
/// clause's line.
///
/// This is the check #2057 exists for. See the module docs for the reasoning and
/// for the classes it deliberately skips.
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

    let text = claude_md();
    let mut faults = Vec::new();
    let mut checked = 0usize;
    for sentence in sentences(&text) {
        let citations = citations_in(sentence);
        if citations.is_empty() {
            continue;
        }
        let lowered = sentence.to_lowercase();
        let triggered: Vec<&(&Anchor, usize)> = located
            .iter()
            .filter(|(anchor, _)| anchor.triggers.iter().any(|t| lowered.contains(t)))
            .collect();
        // Zero triggers: the sentence says nothing that identifies a clause. Two
        // or more: it talks about two, and nothing attributes either to this
        // citation. Two or more citations: the trigger cannot be attributed to
        // one of them. All three are skips, not verdicts.
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
            "{CLAUDE_MD} cites {} in a sentence about the `{}` clause, which is at \
             {SPEC_FILE}:{expected}:\n      {}",
            citation.spelling,
            anchor.identity,
            sentence.trim()
        ));
    }

    assert!(
        faults.is_empty(),
        "{} `general.md` citation(s) in {CLAUDE_MD} name a different line from the clause their \
         own sentence is about:\n  {}\n\n\
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
/// Both are hand-rolled, and both have a failure mode that looks like success: a
/// citation scanner that finds nothing and a sentence splitter that returns one
/// giant sentence each leave the checks above passing while judging much less
/// than they claim.
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
