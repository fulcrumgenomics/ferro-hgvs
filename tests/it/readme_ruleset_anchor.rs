//! Guard for #2143: the numbered normalization ruleset must be hosted in `README.md`.
//!
//! House choices are cited repo-wide as "`README.md` rule N" — by `HouseRule::label`
//! (rendered into `docs/NORMALIZATION_CONTRACT.md`), by the ruling ledger's rationale
//! prose, and by ~60 code and doc sites. #2121 relocated the ruleset text into
//! `docs/src/reference/normalization-rules.md` and slimmed the README, which left every
//! "`README.md` rule N" citation pointing at a file that no longer stated the rules — the
//! dangling citation #2143 reported. The fix re-anchors a compact numbered list in the
//! README; this test pins that anchor so a later slim cannot silently re-break the
//! citations.
//!
//! `include_str!` resolves relative to this source file, so it finds the repo-root
//! `README.md` under both the `it` binary and the `ferro-hgvs-soak-tests` binary (which
//! `#[path]`-includes this module from a different `CARGO_MANIFEST_DIR`).

const README: &str = include_str!("../../README.md");

/// Every rule number cited as "`README.md` rule N" resolves to a numbered item here.
#[test]
fn readme_hosts_the_numbered_normalization_ruleset() {
    const HEADING: &str = "## Normalization ruleset";
    let start = README.find(HEADING).unwrap_or_else(|| {
        panic!(
            "README.md is missing the `## Normalization ruleset` section — house choices are \
             cited repo-wide as `README.md` rule N, so the ruleset must be hosted here (#2143)."
        )
    });
    // Restrict the rule-presence checks to this section — from its heading to the next `##`
    // heading (or EOF) — so an unrelated numbered bold item elsewhere in the README cannot
    // satisfy them after the canonical ruleset loses a rule.
    let body = &README[start..];
    let end = body[HEADING.len()..]
        .find("\n## ")
        .map_or(body.len(), |offset| offset + HEADING.len());
    let section = &body[..end];
    for n in 1..=7u8 {
        // Assert the property (the rule is present as a list item), not an exact rendering.
        let item = format!("\n{n}. **");
        assert!(
            section.contains(&item),
            "README.md is missing normalization rule {n} (looked for a `{n}. **…` list item \
             in the `## Normalization ruleset` section). If the ruleset moved out of the README, \
             re-anchor it here or the ~60 repo-wide `README.md` rule N citations dangle (#2143)."
        );
    }
}
