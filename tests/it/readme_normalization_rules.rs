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
//! Every assertion is scoped to the section it is about ([`section`]) rather than
//! to the whole file, and the two cross-links are matched as whole Markdown tokens
//! rather than as bare `#anchor` fragments. Both are the difference between "this
//! string exists in a long document" and "the ruleset states this" — a guard
//! satisfiable from a code block or an unrelated section is one that passes while
//! the thing it guards is broken.
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
    "### Known limitation",
];

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
