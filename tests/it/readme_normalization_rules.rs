//! Drift guard for the **Normalization rules** section of `README.md`.
//!
//! That section is the project's stated contract for what the normalizer's output
//! is allowed to be: four rules about the output, three about the gaps. It is cited
//! by rule number — from PR descriptions, from ruling records, and from
//! `CONTRIBUTING.md`, whose "Declaring a representation change" section is rule 7's
//! mechanism. A silent renumbering or a dropped rule turns every one of those
//! citations into a pointer at the wrong text, and prose carries no compiler.
//!
//! This follows the shape of `tests/python/test_representation_change_trailer.py`,
//! which pins the decline vocabulary across `release-plz.toml`, `CONTRIBUTING.md`
//! and `scripts/check_representation_change.py` for the same reason: three copies of
//! one statement drift, and the prose copy drifts first.
//!
//! Deliberately **not** a markdown parser. It asserts that seven known rule openers
//! are present, exactly once each, in ascending order, and each under the subsection
//! that classifies it — enough that an edit to the ruleset has to be deliberate, and
//! little enough that rewording the body of a rule does not fail the build.
//!
//! Every assertion is scoped to the section it is about rather than to the whole
//! file — [`section`] for Markdown, [`banner_block`] for the `// Record N —`
//! blocks of an adjudication test module — and the two cross-links are matched as
//! whole Markdown tokens rather than as bare `#anchor` fragments. Both are the
//! difference between "this string exists in a long document" and "the ruleset
//! states this" — a guard satisfiable from a code block or an unrelated section is
//! one that passes while the thing it guards is broken.
//!
//! That applies to the *pinned* side of a cross-claim as much as to the README
//! side, which is why [`banner_block`] exists: a whole-file `contains` over a
//! 500-line Rust module is the same weak question in a different language.
//!
//! Scoping is applied on **both** sides of each claim, which is a strictly stronger
//! test than scoping one side. Pinning the rules' order without binding each rule to
//! its half lets rule 4 be reclassified from the output contract to the procedure
//! with every assertion still green; asserting CONTRIBUTING.md's backlink against
//! the whole file lets it migrate out of the section a reader arriving from rule 7
//! actually lands in.

use std::path::PathBuf;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read(relative: &str) -> String {
    let path = repo_root().join(relative);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
}

/// The Markdown section introduced by `heading`, up to the next heading at the
/// same or a shallower level (or end of input).
///
/// Every assertion below runs against a slice rather than the whole file,
/// because a whole-file search answers a weaker question than the one being
/// asked. `README.md` is long and full of HGVS prose and fenced examples, so a
/// rule opener or an organising heading appearing *anywhere* — a later section,
/// a code block, a changelog snippet — would satisfy a whole-file `find` while
/// the ruleset itself was incomplete or renumbered. That is the failure this
/// guard exists to catch, so it must not be satisfiable from outside the
/// section it guards.
///
/// The match is **line-anchored** (the whole trimmed line must equal `heading`)
/// rather than a substring `find`. A substring search for `## Normalization
/// rules` also matches inside a `### Normalization rules` heading, since `##`
/// is a prefix of `###` — so a deeper heading added later would silently move
/// where the slice begins.
fn section<'a>(text: &'a str, heading: &str) -> &'a str {
    let level = heading.chars().take_while(|c| *c == '#').count();
    assert!(level > 0, "`{heading}` is not a Markdown heading");

    let mut start = None;
    let mut end = text.len();
    let mut offset = 0usize;
    for line in text.split_inclusive('\n') {
        if start.is_none() {
            if line.trim_end() == heading {
                start = Some(offset);
            }
        } else {
            // A heading at the same or a shallower level closes this section.
            let depth = line.chars().take_while(|c| *c == '#').count();
            if depth > 0 && depth <= level && line[depth..].starts_with(' ') {
                end = offset;
                break;
            }
        }
        offset += line.len();
    }

    let start = start.unwrap_or_else(|| panic!("no `{heading}` heading found"));
    &text[start..end]
}

/// The body of a `// Record N — …` banner block in an adjudication test module,
/// up to the next banner rule (or end of input).
///
/// The Rust counterpart of [`section`], and it exists for the same reason: the
/// pinned side of a cross-claim has to be scoped too, or the claim degrades to
/// "this token appears somewhere in a 500-line file". `cis_confluence_adjudication.rs`
/// carries four records separated by
///
/// ```text
/// // ---------------------------------------------------------------------------
/// // Record 1 — settled: separation is a property of the spelling
/// // ---------------------------------------------------------------------------
/// ```
///
/// so the banner is the structure to key on. `title` is matched line-anchored,
/// exactly as `section` matches a heading, and the block ends at the next rule
/// line — which is the *opening* rule of the following record, since a title's
/// own closing rule is consumed first.
///
/// **What this does not buy**, stated because overclaiming a guard is worse than
/// not having one: it scopes to the record, not to an executed assertion. A
/// token sitting only in a doc comment inside the right record still passes.
/// What it rules out is the failure that actually happens here — the case
/// migrating to a different record, or surviving only as commented-out text
/// elsewhere in the file — while the README goes on citing it.
fn banner_block<'a>(text: &'a str, title: &str) -> &'a str {
    /// The prefix a banner rule line starts with; long enough not to match an
    /// ordinary `// --` comment.
    const RULE: &str = "// -----";

    // `(offset, line)` pairs, so a slice can be taken from any line boundary.
    let lines: Vec<(usize, &str)> = {
        let mut offset = 0usize;
        text.split_inclusive('\n')
            .map(|line| {
                let at = offset;
                offset += line.len();
                (at, line)
            })
            .collect()
    };

    let title_index = lines
        .iter()
        .position(|(_, line)| line.trim_end() == title)
        .unwrap_or_else(|| panic!("no `{title}` banner found"));

    // The title's own closing rule is the first rule at or after it; the block
    // is everything from just past that rule to the next rule line.
    let is_rule = |index: usize| lines[index].1.trim_start().starts_with(RULE);
    let closing = (title_index + 1..lines.len())
        .find(|index| is_rule(*index))
        .unwrap_or_else(|| panic!("the `{title}` banner has no closing rule line"));
    let start = lines[closing].0 + lines[closing].1.len();
    let end = (closing + 1..lines.len())
        .find(|index| is_rule(*index))
        .map_or(text.len(), |index| lines[index].0);

    &text[start..end]
}

/// The section heading, and the subsection headings that organise it.
///
/// Listed in the order they must appear: the two halves of the ruleset first
/// (`The output contract` holds rules 1-4, `The procedure` holds 5-7), then the
/// explanatory material.
const HEADINGS: &[&str] = &[
    "## Normalization rules",
    "### The output contract",
    "### The procedure",
    "### Why 2 and 3 are best effort, and 1 and 4 are not",
    "### A worked example of reading force from prose",
    "### What rule 3 excludes",
    "### Known limitation",
];

/// The coordinates the README's rule 3 example is stated on.
///
/// The example exists in two places by necessity — the README states the rule and
/// `cis_confluence_adjudication.rs` pins the behaviour — so it is exactly the shape
/// this repository drifts on. Neither copy is the authority for the other's *prose*;
/// what must not diverge is the case itself, since a README example on coordinates
/// no test covers is a claim nothing checks.
const RULE_3_EXAMPLE_TOKENS: &[&str] = &[
    "NC_000001.11",
    "ATGAGGGGCCACTGT",
    "g.[1001009del;1001010del;1001013del]",
    "g.[1001009del;1001011del;1001013del]",
    "g.[1001009_1001010del;1001013del]",
];

/// The banner block of `cis_confluence_adjudication.rs` that owns the rule 3
/// example, so the pinned side of the cross-claim is scoped like the README side.
const RULE_3_PINNED_RECORD: &str =
    "// Record 1 — settled: separation is a property of the spelling";

/// The seven rules, by number and name, exactly as the README opens each one.
///
/// The name is part of the guard because the numbers alone are what a renumbering
/// preserves: swapping two rules keeps `1.`..`7.` intact and silently re-points
/// every "rule 5" citation in the repository.
const RULE_OPENERS: &[&str] = &[
    "1. **Conformant.**",
    "2. **Recommended form.**",
    "3. **Confluent.**",
    "4. **Deterministic.**",
    "5. **Where the spec is silent, ambiguous, or self-contradictory:**",
    "6. **Among multiple conformant forms:**",
    "7. **Disclosure.**",
];

#[test]
fn readme_states_all_seven_normalization_rules() {
    let readme = read("README.md");
    let rules = section(&readme, "## Normalization rules");
    let mut cursor = 0usize;
    for opener in RULE_OPENERS {
        assert_eq!(
            rules.matches(opener).count(),
            1,
            "the `## Normalization rules` section must state `{opener}` exactly once; \
             a rule was dropped, duplicated, or reworded"
        );
        let at = rules
            .find(opener)
            .expect("checked to occur exactly once just above");
        assert!(
            at >= cursor,
            "the `## Normalization rules` section states `{opener}` out of order; \
             the rules are cited by number, so reordering them re-points every \
             existing citation"
        );
        cursor = at;
    }
}

/// Each rule must sit under the subsection that classifies it.
///
/// [`readme_states_all_seven_normalization_rules`] pins the openers' *order*,
/// which is a strictly weaker claim than the one the module doc makes: it walks
/// the rules with one cursor and the headings with another, and the two never
/// meet. Rule 4 could move below `### The procedure` — reclassifying an
/// **Absolute** output-contract guarantee as procedure — and both walks would
/// still pass, because the relative order of every opener and every heading is
/// unchanged. The half a rule lives in is what says whether it is a promise
/// about output or a step in handling a gap, so it is the part worth pinning.
#[test]
fn each_rule_sits_under_the_heading_that_classifies_it() {
    const OUTPUT_CONTRACT: usize = 4; // rules 1-4; the rest belong to the procedure

    let readme = read("README.md");
    let rules = section(&readme, "## Normalization rules");
    let output_contract = section(rules, "### The output contract");
    let procedure = section(rules, "### The procedure");

    let (contract_rules, procedure_rules) = RULE_OPENERS.split_at(OUTPUT_CONTRACT);
    for (half, openers) in [
        ("### The output contract", contract_rules),
        ("### The procedure", procedure_rules),
    ] {
        let slice = if half == "### The output contract" {
            output_contract
        } else {
            procedure
        };
        for opener in openers {
            assert_eq!(
                slice.matches(opener).count(),
                1,
                "`{opener}` must be stated exactly once under `{half}`; a rule that \
                 crosses between the two halves changes what it claims — the output \
                 contract binds ferro's output, the procedure only binds how gaps \
                 are handled"
            );
        }
    }
}

#[test]
fn readme_keeps_the_headings_that_organise_the_rules() {
    let readme = read("README.md");
    let rules = section(&readme, "## Normalization rules");
    let mut cursor = 0usize;
    for heading in HEADINGS {
        // Exactly once, not merely present: a duplicated organising heading
        // splits the ruleset into two places, and `find` would happily anchor on
        // whichever came first.
        assert_eq!(
            rules.matches(heading).count(),
            1,
            "the `## Normalization rules` section must carry the heading `{heading}` \
             exactly once"
        );
        let at = rules
            .find(heading)
            .expect("checked to occur exactly once just above");
        assert!(
            at >= cursor,
            "the `## Normalization rules` section states `{heading}` out of order"
        );
        cursor = at;
    }
}

#[test]
fn the_ruleset_and_its_disclosure_mechanism_point_at_each_other() {
    // `CONTRIBUTING.md`'s trailer section is rule 7's mechanism, and rule 7 is that
    // section's reason for existing. Either pointer going stale leaves a reader on
    // one half of one statement.
    let readme = read("README.md");
    let contributing = read("CONTRIBUTING.md");

    // Whole Markdown link tokens, not bare fragments: `contains("…#anchor")`
    // also passes on prose that merely names the anchor, or on a stale reference
    // left behind after the link itself was deleted.
    //
    // Both sides are asserted **inside** the section the link belongs to. The
    // README side must sit in the ruleset, since rule 7 is what it belongs to;
    // the CONTRIBUTING side must sit in `Declaring a representation change`,
    // since that section is rule 7's mechanism. A whole-file `contains` on
    // either side passes while the backlink has migrated somewhere the reader
    // arriving from the other document will never reach — which is precisely
    // the drift these two pointers exist to prevent.
    //
    // `section()` panics with a clear message when its heading is absent, so
    // each call below doubles as the "the anchor still exists" assertion; the
    // link targets are the two headings themselves.
    let rules = section(&readme, "## Normalization rules");
    let declaring = section(&contributing, "### Declaring a representation change");

    assert!(
        rules.contains("[CONTRIBUTING.md](CONTRIBUTING.md#declaring-a-representation-change)"),
        "the `## Normalization rules` section must link to CONTRIBUTING.md's \
         `Declaring a representation change` section, which is rule 7's mechanism"
    );
    assert!(
        declaring.contains("[normalization rules](README.md#normalization-rules)"),
        "CONTRIBUTING.md's `Declaring a representation change` section must link \
         back to the README's normalization rules, so the trailer's purpose is \
         reachable from the trailer's instructions"
    );
}

/// The ledger record that names `README.md` as canonical must find the section
/// it names.
///
/// The `adjudication-precedence-order` record was rewritten from a restatement
/// of the ruleset into a **pointer** at it. A pointer has a failure mode a
/// restatement does not: it can dangle. Nothing asserted that it resolved —
/// `readme_normalization_rules.rs` guards the README against itself and knows
/// nothing of the ledger, and `ruling_citation_currency.rs` checks the ledger's
/// citations of *spec clauses*, not of this repository's own documents.
///
/// The gap was live, and named at the time it was live: while the ruleset PR was
/// unmerged, this record's opening claim — that the ruleset "lives in
/// `README.md`" — was false on `main`, and the required merge order was carried
/// only by a PR comment. The guard could not be written on that branch because
/// it would have been red until the README section landed. It has, so this is
/// that guard, and from here the order is enforced by the suite rather than
/// remembered.
///
/// Deliberately checks that the record *is still a pointer*, not only that the
/// pointer resolves. "Do not restate the rules" is the substance of the ruling,
/// and a record that quietly grew a copy of them would satisfy a
/// resolves-only check while recreating the drift the ruling exists to stop.
#[test]
fn the_ledgers_pointer_at_the_readme_ruleset_resolves() {
    let ledger_path =
        repo_root().join("tests/fixtures/grammar/hgvs_spec_normalization_overrides.json");
    let text = std::fs::read_to_string(&ledger_path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", ledger_path.display()));
    let ledger: serde_json::Value =
        serde_json::from_str(&text).expect("the overrides ledger is valid JSON");

    let record = ledger["rulings"]
        .as_array()
        .expect("`rulings` is an array")
        .iter()
        .find(|r| r["id"] == "adjudication-precedence-order")
        .expect(
            "the `adjudication-precedence-order` record must exist; it is what names the README \
             ruleset as canonical",
        );
    let rationale = record["rationale"].as_str().expect("a string rationale");

    assert!(
        rationale.contains("README.md"),
        "`adjudication-precedence-order` no longer names `README.md`. If the canonical ruleset \
         moved, move this guard with it; if the record went back to restating the rules, that is \
         the drift the ruling forbids"
    );

    // `section` panics with the heading it could not find, so this call IS the
    // "the pointer resolves" assertion — the same construction
    // `the_ruleset_and_its_disclosure_mechanism_point_at_each_other` relies on.
    let readme = read("README.md");
    let rules = section(&readme, "## Normalization rules");
    assert!(
        !rules.trim().is_empty(),
        "the `## Normalization rules` section exists but is empty, so the ledger's pointer \
         resolves to nothing"
    );

    // Still a pointer, not a restatement. The rules are numbered `1.`..`7.` in
    // the README; a record that had grown its own copy would carry several of
    // their bolded names.
    let restated: Vec<&str> = [
        "**Confluent.**",
        "**Deterministic.**",
        "**Idempotent.**",
        "**Disclosed.**",
    ]
    .into_iter()
    .filter(|name| rationale.contains(name))
    .collect();
    assert!(
        restated.is_empty(),
        "`adjudication-precedence-order` has started restating the ruleset it is supposed to \
         point at ({restated:?}). The ruling is that the rules are stated in exactly one place"
    );
}

/// The README's rule 3 example must still be the case the adjudication test pins.
///
/// [`the_ruleset_and_its_disclosure_mechanism_point_at_each_other`] guards a
/// cross-*link*; this guards a cross-*claim*, which is the weaker link of the two.
/// A link breaks loudly when its anchor moves. An example does not: the README can
/// go on quoting coordinates long after the test stopped covering them, or the test
/// can be re-pinned onto a different run, and both documents still read as though
/// they agree. That is the failure this repository has recorded against itself
/// repeatedly — one fact in several places, drifting apart quietly.
///
/// Deliberately token-based rather than a diff. The two texts *should* differ: one
/// states a rule to a reader, the other exercises a normalizer. Only the case they
/// are both about — the reference window, the two spellings, and the decided output
/// — has to stay shared.
///
/// **Both sides are scoped**, which is what the module doc above claims of every
/// assertion here and what this test would otherwise have been the one exception
/// to. The README side is bound to `### What rule 3 excludes` by [`section`]; the
/// pinned side is bound to record 1's banner block by [`banner_block`], since a
/// whole-file `contains` over a 500-line module answers the weaker question
/// [`section`]'s own doc argues against — it would stay green with the case moved
/// into an unrelated record, or left behind as commented-out text.
#[test]
fn the_rule_3_example_is_the_case_the_adjudication_test_pins() {
    let readme = read("README.md");
    let rules = section(&readme, "## Normalization rules");
    let example = section(rules, "### What rule 3 excludes");
    let pinned_file = read("tests/it/cis_confluence_adjudication.rs");
    let pinned = banner_block(&pinned_file, RULE_3_PINNED_RECORD);

    for token in RULE_3_EXAMPLE_TOKENS {
        assert!(
            example.contains(token),
            "the README's rule 3 example no longer states `{token}`; if the example moved to \
             different coordinates, move `RULE_3_EXAMPLE_TOKENS` and the pinned test with it"
        );
        assert!(
            pinned.contains(token),
            "the README's rule 3 example states `{token}`, which \
             `tests/it/cis_confluence_adjudication.rs` no longer covers under \
             `{RULE_3_PINNED_RECORD}`. The README would be quoting a case nothing checks — or \
             the case moved to another record, in which case move `RULE_3_PINNED_RECORD` too"
        );
    }

    // The link is what sends a reader from the claim to its evidence, so it is
    // part of the claim rather than a courtesy.
    assert!(
        example.contains(
            "[`tests/it/cis_confluence_adjudication.rs`](tests/it/cis_confluence_adjudication.rs)"
        ),
        "the rule 3 example must link to the test that pins it"
    );
}
