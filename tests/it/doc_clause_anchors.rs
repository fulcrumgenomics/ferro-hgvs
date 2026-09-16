//! A `general.md:N` clause citation in the project's prose docs must name a
//! line the clause is actually on in the pinned spec checkout — not merely a
//! line that exists.
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
//! protected had moved into `docs/READING_THE_SPEC.md` and `CONTRIBUTING.md`,
//! where nothing scanned them — this file closes that gap.
//!
//! # Why a line-keyed check, and how it stays honest
//!
//! Each clause carries a **signature**: a phrase found on exactly one line of
//! the spec checkout, so the expected line number is *derived from the spec on
//! every run* and never restated here. The set of lines the signatures resolve
//! to is the set of clause lines the project litigates. Every `general.md:N`
//! citation in the docs must name one of those lines; a bump that shifts a
//! clause moves the derived line out from under a stale citation and reddens
//! this test.
//!
//! Unlike the CLAUDE.md guard this succeeds, the docs cite `general.md` clauses
//! mostly by enumeration (`general.md:33`, `:34`, `:55`, `:57`) rather than one
//! clause per argumentative sentence, so a sentence-trigger check has almost
//! nothing to judge. The line-set check covers every citation instead, which is
//! the right shape for an enumeration and still catches the uniform line shift
//! that motivated the original guard. It is deliberately blind to a pure swap of
//! two adjacent clauses (both lines stay litigated); that is not the failure a
//! submodule bump produces.
//!
//! Only `general.md` is checked, exactly as the CLAUDE.md guard did: a bare
//! `delins.md:` is ambiguous across `DNA/`, `RNA/` and `protein/`, so those
//! citations are not line-keyed here.

use std::path::PathBuf;

/// The pinned spec checkout, relative to the crate root.
const SPEC_DIR: &str = "assets/hgvs-nomenclature";
/// The one spec file whose citations are unambiguous by name.
const SPEC_FILE: &str = "docs/recommendations/general.md";

/// The prose docs that cite `general.md:N` clauses. `docs/TESTING.md` and
/// `CLAUDE.md` carry none today but are scanned so a citation added to them
/// later is line-keyed from the start.
const DOCS: &[&str] = &[
    "docs/READING_THE_SPEC.md",
    "CONTRIBUTING.md",
    "docs/TESTING.md",
    "CLAUDE.md",
];

/// A phrase that occurs on exactly one line of [`SPEC_FILE`]. The line it is
/// found on IS the expectation; nothing here restates a line number. These are
/// the same signatures [`ledger_prose_clause_anchors`] uses, for the clauses the
/// docs actually cite.
const SIGNATURES: &[(&str, &str)] = &[
    (
        "separation-rule",
        r#"should be described individually and **not** as a "delins""#,
    ),
    (
        "codon-exception",
        "**exception**: two variants separated by one nucleotide, together affecting \
         one amino acid",
    ),
    (
        "prioritisation",
        "the preferred description is: (1) substitution, (2) deletion, (3) inversion",
    ),
    (
        "self-replacement-prohibition",
        "descriptions removing part of a reference sequence and replacing it with part \
         of the same sequence are not allowed",
    ),
    (
        "svd-wg-forward-note",
        "the SVD-WG is preparing a proposal to modify this recommendation",
    ),
];

/// Floor for how many `general.md:N` citations this test actually checks.
///
/// Measured at 7 on the commit that added this file (the `general.md:33/:34/:55/
/// :57` enumeration and the `:57` and `:35-38` citations in
/// `docs/READING_THE_SPEC.md`, plus the `general.md:33` example in
/// `CONTRIBUTING.md`). The floor sits below that, so ordinary prose edits cannot
/// trip it while a citation scanner that stopped finding citations would.
const CITATIONS_FLOOR: usize = 6;

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

/// Whitespace-collapsed and lower-cased, for substring comparison — the same
/// normalisation `ledger_prose_clause_anchors` applies, so a signature written
/// in the spec's own words matches regardless of wrapping or case.
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

/// A `general.md:N` or `general.md:N-M` citation found in a doc.
struct Citation {
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

/// Every `general.md` line cited in `text`, qualified (`general.md:N`) or bare
/// (`:N` after a `general.md:` establisher). Bare shorthands are attributed to
/// the most recent `<file>.md:` token on the SAME line, and the attribution
/// resets at each newline, so a bare `:N` never borrows a file from an unrelated
/// line. Citations are read from inline-code spans (backtick-delimited), which
/// is how the docs write them.
fn general_md_citations(text: &str, source_name: &str) -> Vec<Citation> {
    let mut found = Vec::new();
    for (row, line) in text.lines().enumerate() {
        let mut current_is_general = false;
        for (index, span) in line.split('`').enumerate() {
            // Even indices are outside backticks; only inline code can be a citation.
            if index % 2 == 0 {
                continue;
            }
            let record = |found: &mut Vec<Citation>, digits: &str, spelling: &str| {
                if let Some((first, last)) = parse_range(digits) {
                    found.push(Citation {
                        first,
                        last,
                        spelling: spelling.to_string(),
                        source: format!("{source_name}:{}", row + 1),
                    });
                }
            };
            if let Some(rest) = span.strip_prefix(':') {
                if current_is_general {
                    record(&mut found, rest, span);
                }
            } else if let Some(marker) = span.find(".md:") {
                let file = &span[..marker + ".md".len()];
                current_is_general = file.ends_with("general.md");
                if current_is_general {
                    record(&mut found, &span[marker + ".md:".len()..], span);
                }
            }
            // Any other inline code (e.g. `dup`) leaves the current file intact:
            // it neither is a citation nor establishes a new file context.
        }
    }
    found
}

/// Read every doc's `general.md` citations.
fn all_citations() -> Vec<Citation> {
    let root = crate_root();
    let mut all = Vec::new();
    for doc in DOCS {
        let path = root.join(doc);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        all.extend(general_md_citations(&text, doc));
    }
    all
}

/// The set of `general.md` line numbers the project's clauses currently occupy,
/// derived from the signatures on every run. Panics if a signature no longer
/// names exactly one line — that is a broken anchor, not a citation problem.
fn litigated_lines() -> Vec<usize> {
    let lines = spec_lines();
    SIGNATURES
        .iter()
        .map(|(identity, signature)| {
            let hits = lines_matching(&lines, signature);
            assert_eq!(
                hits.len(),
                1,
                "the `{identity}` signature must occur on exactly one line of {SPEC_FILE}, \
                 found {hits:?}. The spec text moved or the signature drifted; update the \
                 signature to a phrase unique to that clause."
            );
            hits[0]
        })
        .collect()
}

#[test]
fn every_signature_locates_exactly_one_line() {
    // `litigated_lines` asserts uniqueness; calling it is the check.
    let lines = litigated_lines();
    assert_eq!(
        lines.len(),
        SIGNATURES.len(),
        "every signature must resolve to a line"
    );
}

#[test]
fn every_doc_general_md_citation_names_a_litigated_line() {
    let litigated = litigated_lines();
    let citations = all_citations();

    assert!(
        citations.len() >= CITATIONS_FLOOR,
        "only {} `general.md` citations found across {DOCS:?}; expected at least \
         {CITATIONS_FLOOR}. The citation scanner probably stopped matching — a real drop \
         to zero would mean the docs no longer cite the spec by line.",
        citations.len()
    );

    let mut stale = Vec::new();
    for citation in &citations {
        // A single line must be litigated; a range must contain a litigated line.
        let names_a_clause = (citation.first..=citation.last).any(|line| litigated.contains(&line));
        if !names_a_clause {
            stale.push(format!(
                "{} cites `{}`, which is not a line any litigated clause is on \
                 (current clause lines: {litigated:?})",
                citation.source, citation.spelling
            ));
        }
    }

    assert!(
        stale.is_empty(),
        "these doc citations no longer name a clause line in the pinned spec — the \
         submodule almost certainly moved under them:\n  {}",
        stale.join("\n  ")
    );
}
