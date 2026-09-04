//! A clause citation in the project's prose docs — `general.md:N`, or a
//! path-qualified `DNA/delins.md:N` and its siblings — must name a line the
//! clause is actually on in the pinned spec checkout, not merely a line that
//! exists.
//!
//! # The failure this closes
//!
//! This is the docs sibling of [`ledger_prose_clause_anchors`], which guards the
//! ruling ledger's own prose the same way. The #1793 submodule bump
//! (`6f85311` -> `565b973`) deleted a duplicated line at old `general.md:22` and
//! shifted every later line down by one, moving every clause the project cites
//! while leaving the numbers written beside them untouched. The ledger sweep
//! (#2009/#2014) reached the ledger; the CLAUDE.md guards (#2057) reached
//! CLAUDE.md. When the CLAUDE.md guards were retired, the clause citations they
//! protected had moved into `docs/READING_THE_SPEC.md`, `docs/TESTING.md` and
//! `CONTRIBUTING.md`, where nothing scanned them — this file closes that gap.
//!
//! # Why a line-keyed check, and how it stays honest
//!
//! Each clause carries a **signature**: a phrase found on exactly one line of
//! its spec file, so the expected line number is *derived from the spec on every
//! run* and never restated here. The set of lines the signatures resolve to is
//! the set of clause lines the project litigates. Every citation to one of those
//! files must name one of those lines; a bump that shifts a clause moves the
//! derived line out from under a stale citation and reddens this test.
//!
//! The docs cite these clauses mostly by enumeration (`general.md:33`, `:34`,
//! `:55`; `DNA/delins.md:16`, `:17`, `:18`, `:47`) rather than one clause per
//! argumentative sentence, so a sentence-trigger check has almost nothing to
//! judge. The line-set check covers every citation instead, which is the right
//! shape for an enumeration and still catches the uniform line shift that
//! motivated the original guard. It is deliberately blind to a pure swap of two
//! adjacent clauses (both lines stay litigated); that is not the failure a
//! submodule bump produces.
//!
//! # Which citations are line-keyed, and which are not
//!
//! A citation is line-keyed only when its file is **unambiguous**. `general.md`
//! is unique by name, so a bare `general.md:N` is safe. A bare `delins.md:N`
//! would be ambiguous across `DNA/`, `RNA/` and `protein/`, so it is *not*
//! keyed — but the docs cite these clauses **path-qualified** (`DNA/delins.md`,
//! `DNA/inversion.md`, `DNA/duplication.md`, `background/basics.md`), which names
//! the file exactly and is what makes them keyable here. A bare, unqualified
//! `delins.md:N` is left unchecked by design; qualify it to bring it under the
//! guard. This also means a bare `:N` shorthand is only keyed when its
//! `<file>.md:` establisher sits on the **same source line** (attribution resets
//! at each newline), so keep a bare `:N` on the line of the token it inherits,
//! or spell the qualified form out when a wrap would separate them.
//!
//! The files line-keyed here are `general.md`, `DNA/delins.md`,
//! `DNA/inversion.md`, `DNA/duplication.md` and `background/basics.md` (see
//! [`SPEC_FILES`]). The docs also cite other spec files — `style.md`,
//! `RNA/adjoined_transcript.md`, `DNA/other.md`, `DNA/repeated.md`,
//! `RNA/repeated.md`, `background/refseq.md`, and the `consultation/` proposals —
//! whose line numbers are **not** machine-verified. Those citations must be
//! re-checked by hand after a submodule bump; add a `SpecFile` with a unique
//! signature to bring one under the guard.

use std::path::PathBuf;

/// The pinned spec checkout, relative to the crate root.
const SPEC_DIR: &str = "assets/hgvs-nomenclature";

/// A phrase unique to one line of a spec file. The line it is found on IS the
/// expectation; nothing here restates a line number.
struct Signature {
    /// A short name for the clause, for the failure message only.
    identity: &'static str,
    /// The phrase, in the spec's own words. Matched whitespace-collapsed and
    /// case-insensitively, so wrapping and case do not matter.
    phrase: &'static str,
}

/// A spec file whose citations this test line-keys, with the phrase signatures
/// for the clauses the docs cite in it.
struct SpecFile {
    /// The citation token the docs use. `general.md` is cited bare (as
    /// `general.md:N`); the `DNA/` files are cited path-qualified (as
    /// `DNA/delins.md:N`). A `<token>:N` citation — and a bare `:N` following it
    /// on the same line — is checked against this file. The match is exact, so a
    /// bare `delins.md:N` does not resolve to `DNA/delins.md` and is left
    /// unchecked, which is the ambiguity guard.
    token: &'static str,
    /// Path under [`SPEC_DIR`] to the file the signatures are found in.
    path: &'static str,
    /// Phrases, each unique to one line of `path`.
    signatures: &'static [Signature],
}

/// `general.md` clauses. These are the same signatures the CLAUDE.md guard used,
/// for the clauses the docs actually cite by line.
const GENERAL_SIGNATURES: &[Signature] = &[
    Signature {
        identity: "separation-rule",
        phrase: r#"should be described individually and **not** as a "delins""#,
    },
    Signature {
        identity: "codon-exception",
        phrase: "**exception**: two variants separated by one nucleotide, together affecting \
                 one amino acid",
    },
    Signature {
        identity: "prioritisation",
        phrase: "the preferred description is: (1) substitution, (2) deletion, (3) inversion",
    },
    Signature {
        identity: "three-prime-rule",
        phrase: "the most 3' position possible of the reference sequence is arbitrarily assigned",
    },
    Signature {
        identity: "self-replacement-prohibition",
        phrase: "descriptions removing part of a reference sequence and replacing it with part \
                 of the same sequence are not allowed",
    },
    Signature {
        identity: "svd-wg-forward-note",
        phrase: "the SVD-WG is preparing a proposal to modify this recommendation",
    },
];

/// `DNA/delins.md` clauses. The docs cite `:16` (consecutive-nucleotides-are-
/// delins), `:17` (individual-not-delins), `:18` (the codon exception), the
/// worked example at `:44` and its `:46` "alternative description" align-note and
/// `:47` "delins recommended" note, plus the two trailing `!!! note` Q&As the
/// docs walk through: the separated-variants note (reached by `:79-84`) and the
/// BRCA1 modified-answer note (cited both as the `:86-89` range and as `:89`, the
/// line recording that the two-member spelling was removed).
const DELINS_SIGNATURES: &[Signature] = &[
    Signature {
        identity: "delins-consecutive-are-delins",
        phrase: "changes involving two or more consecutive nucleotides are described as \
                 deletion/insertion",
    },
    Signature {
        identity: "delins-separation-rule",
        phrase: "two variants separated by one or more nucleotides should be described \
                 individually",
    },
    Signature {
        identity: "delins-worked-example",
        phrase: "c.850_901delinsTTCCTCGATGCCTG",
    },
    Signature {
        identity: "delins-codon-exception",
        phrase: "**exception**: two variants separated by one nucleotide, together affecting \
                 one amino acid",
    },
    Signature {
        identity: "delins-alternative-description",
        phrase: r#"parts of the inserted sequence "align" with the reference sequence"#,
    },
    Signature {
        identity: "delins-format-recommended",
        phrase: r#"**The "delins" format is recommended**"#,
    },
    Signature {
        identity: "delins-provenance-qa",
        phrase: "the two variants may have been reported",
    },
    Signature {
        identity: "delins-brca1-modified-answer",
        phrase: "the answer was modified",
    },
];

/// `DNA/inversion.md` clauses. The docs cite `:20`, inversion's own copy of the
/// individual-not-delins rule.
const INVERSION_SIGNATURES: &[Signature] = &[Signature {
    identity: "inversion-separation-rule",
    phrase: "two variants separated by one or more nucleotides should be described individually",
}];

/// `DNA/duplication.md` clauses. The docs cite `:18` as read-strong for its
/// "**must** be described as a duplication" wording, and `:86`'s note marking
/// the worked example above it as part of the (undecided) SVD-WG003 proposal.
const DUPLICATION_SIGNATURES: &[Signature] = &[
    Signature {
        identity: "duplication-must-be-duplication",
        phrase: "when a variant can be described as a duplication, it **must** be described as a \
                 duplication",
    },
    Signature {
        identity: "duplication-svd-wg003-note",
        phrase: "proposal SVD-WG003 (undecided)",
    },
];

/// `background/basics.md` clauses. The docs cite `:38` for the spec's stated
/// design values ("stable, meaningful, memorable, unequivocal"), which the
/// "stability is first" argument leans on and repeats, so it is line-keyed.
const BASICS_SIGNATURES: &[Signature] = &[Signature {
    identity: "basics-design-values",
    phrase: "designed to be **stable**, **meaningful**, **memorable**, and **unequivocal**",
}];

/// Every spec file whose citations are line-keyed.
const SPEC_FILES: &[SpecFile] = &[
    SpecFile {
        token: "general.md",
        path: "docs/recommendations/general.md",
        signatures: GENERAL_SIGNATURES,
    },
    SpecFile {
        token: "DNA/delins.md",
        path: "docs/recommendations/DNA/delins.md",
        signatures: DELINS_SIGNATURES,
    },
    SpecFile {
        token: "DNA/inversion.md",
        path: "docs/recommendations/DNA/inversion.md",
        signatures: INVERSION_SIGNATURES,
    },
    SpecFile {
        token: "DNA/duplication.md",
        path: "docs/recommendations/DNA/duplication.md",
        signatures: DUPLICATION_SIGNATURES,
    },
    SpecFile {
        token: "background/basics.md",
        path: "docs/background/basics.md",
        signatures: BASICS_SIGNATURES,
    },
];

/// The prose docs that cite spec clauses by line. `CLAUDE.md` carries none today
/// but is scanned so a citation added to it later is line-keyed from the start.
const DOCS: &[&str] = &[
    "docs/READING_THE_SPEC.md",
    "CONTRIBUTING.md",
    "docs/TESTING.md",
    "CLAUDE.md",
];

/// Floor for how many clause citations this test actually checks.
///
/// The docs cite these tracked spec files many times over, and the count grows
/// as the docs cite more clauses, so an exact figure is deliberately not pinned
/// here — it would go stale on every prose edit that adds a citation. The floor
/// sits well below the actual count, so ordinary edits cannot trip it while a
/// citation scanner that stopped finding citations (a parser regression) would.
const CITATIONS_FLOOR: usize = 12;

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// A spec file's lines, as read from the pinned checkout.
fn spec_lines(spec_path: &str) -> Vec<String> {
    let path = crate_root().join(SPEC_DIR).join(spec_path);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "read {}: {e} — the spec submodule is probably not initialised. Run\n    \
             git -c protocol.file.allow=always submodule update --init {SPEC_DIR}",
            path.display()
        )
    });
    text.lines().map(str::to_string).collect()
}

/// Whitespace-collapsed and lower-cased, for substring comparison — the same
/// normalisation `ledger_prose_clause_anchors` applies, so a signature written
/// in the spec's own words matches regardless of wrapping or case.
fn collapse(text: &str) -> String {
    text.split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .to_lowercase()
}

/// The 1-based lines of `lines` on which `phrase` occurs.
fn lines_matching(lines: &[String], phrase: &str) -> Vec<usize> {
    let needle = collapse(phrase);
    lines
        .iter()
        .enumerate()
        .filter(|(_, line)| collapse(line).contains(&needle))
        .map(|(index, _)| index + 1)
        .collect()
}

/// The set of line numbers the clauses of `spec_file` currently occupy, derived
/// from its signatures on every run. Panics if a signature no longer names
/// exactly one line — that is a broken anchor, not a citation problem.
fn litigated_lines(spec_file: &SpecFile) -> Vec<usize> {
    let lines = spec_lines(spec_file.path);
    spec_file
        .signatures
        .iter()
        .map(|signature| {
            let hits = lines_matching(&lines, signature.phrase);
            assert_eq!(
                hits.len(),
                1,
                "the `{}` signature must occur on exactly one line of {}, found {hits:?}. The \
                 spec text moved or the signature drifted; update the signature to a phrase \
                 unique to that clause.",
                signature.identity,
                spec_file.path
            );
            hits[0]
        })
        .collect()
}

/// A clause citation found in a doc.
struct Citation {
    /// Index into [`SPEC_FILES`] of the file this citation names.
    file: usize,
    first: usize,
    last: usize,
    /// The citation as written, for the failure message.
    spelling: String,
    /// The doc it was found in, and its 1-based line.
    source: String,
}

/// Parse a run of digits, optionally `first-last`, into `(first, last)`.
fn parse_range(digits: &str) -> Option<(usize, usize)> {
    let mut parts = digits.splitn(2, '-');
    let first: usize = parts.next()?.parse().ok()?;
    let last: usize = match parts.next() {
        Some(tail) => tail.parse().ok()?,
        None => first,
    };
    Some((first, last))
}

/// The index into [`SPEC_FILES`] whose token exactly matches `file`, if any.
/// Exact match is the ambiguity guard: a bare `delins.md` does not match the
/// path-qualified `DNA/delins.md` token, so it is left unchecked.
fn tracked_file(file: &str) -> Option<usize> {
    SPEC_FILES.iter().position(|f| f.token == file)
}

/// Every tracked-file clause line cited in `text`, qualified
/// (`DNA/delins.md:N`) or bare (`:N` after a `<file>.md:` establisher). Bare
/// shorthands are attributed to the most recent tracked `<file>.md:` token on
/// the SAME line, and the attribution resets at each newline, so a bare `:N`
/// never borrows a file from an unrelated line. An untracked establisher
/// (`style.md:9`) resets the attribution too, so a bare `:N` after it is not
/// mis-keyed. Citations are read from inline-code spans (backtick-delimited),
/// which is how the docs write them.
fn citations_in(text: &str, source_name: &str) -> Vec<Citation> {
    let mut found = Vec::new();
    for (row, line) in text.lines().enumerate() {
        let mut current: Option<usize> = None;
        for (index, span) in line.split('`').enumerate() {
            // Even indices are outside backticks; only inline code can be a citation.
            if index % 2 == 0 {
                continue;
            }
            let record = |found: &mut Vec<Citation>, file: usize, digits: &str, spelling: &str| {
                if let Some((first, last)) = parse_range(digits) {
                    found.push(Citation {
                        file,
                        first,
                        last,
                        spelling: spelling.to_string(),
                        source: format!("{source_name}:{}", row + 1),
                    });
                }
            };
            if let Some(rest) = span.strip_prefix(':') {
                if let Some(file) = current {
                    record(&mut found, file, rest, span);
                }
            } else if let Some(marker) = span.find(".md:") {
                let file = &span[..marker + ".md".len()];
                current = tracked_file(file);
                if let Some(index) = current {
                    record(&mut found, index, &span[marker + ".md:".len()..], span);
                }
            }
            // Any other inline code (e.g. `dup`) leaves the current file intact:
            // it neither is a citation nor establishes a new file context.
        }
    }
    found
}

/// Read every doc's tracked-file citations.
fn all_citations() -> Vec<Citation> {
    let root = crate_root();
    let mut all = Vec::new();
    for doc in DOCS {
        let path = root.join(doc);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        all.extend(citations_in(&text, doc));
    }
    all
}

#[test]
fn every_signature_locates_exactly_one_line() {
    // `litigated_lines` asserts uniqueness; calling it for every file is the check.
    let total: usize = SPEC_FILES
        .iter()
        .map(|spec_file| litigated_lines(spec_file).len())
        .sum();
    let expected: usize = SPEC_FILES.iter().map(|f| f.signatures.len()).sum();
    assert_eq!(total, expected, "every signature must resolve to a line");
}

#[test]
fn every_doc_clause_citation_names_a_litigated_line() {
    let litigated: Vec<Vec<usize>> = SPEC_FILES.iter().map(litigated_lines).collect();
    let citations = all_citations();

    assert!(
        citations.len() >= CITATIONS_FLOOR,
        "only {} clause citations found across {DOCS:?}; expected at least {CITATIONS_FLOOR}. \
         The citation scanner probably stopped matching — a real drop to zero would mean the \
         docs no longer cite the spec by line.",
        citations.len()
    );

    let mut stale = Vec::new();
    for citation in &citations {
        let lines = &litigated[citation.file];
        // A single line must be litigated; a range must contain a litigated line.
        let names_a_clause = (citation.first..=citation.last).any(|line| lines.contains(&line));
        if !names_a_clause {
            stale.push(format!(
                "{} cites `{}`, which is not a line any litigated clause of {} is on \
                 (current clause lines: {lines:?})",
                citation.source, citation.spelling, SPEC_FILES[citation.file].path
            ));
        }
    }

    assert!(
        stale.is_empty(),
        "these doc citations no longer name a clause line in the pinned spec — the submodule \
         almost certainly moved under them:\n  {}",
        stale.join("\n  ")
    );
}
