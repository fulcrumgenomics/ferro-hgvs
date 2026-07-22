//! MUST-level rejection of a frameshift anchored at an **unchanged** amino
//! acid — `p.His150HisfsTer10` / `p.His150Hisfs*10` — issue #1079.
//!
//! # Spec
//!
//! `recommendations/protein/frameshift.md:47-49`:
//!
//! > `p.Gln151Thrfs*9` (not `p.His150Hisfs*10`) — […] Since frameshift variants
//! > start with the **first amino acid changed**, the description
//! > `p.His150HisfsTer10` (or `p.His150Hisfs*10`) is not correct.
//!
//! A frameshift must anchor at the first residue the shifted reading frame
//! actually changes. `p.His150His…` names `His150` as both the reference and
//! the new residue — an unchanged anchor — so the description is spec-invalid.
//!
//! # Rejection, not repair
//!
//! The correct form (`p.Gln151Thrfs*9`) names the first *changed* residue and
//! the length of the shifted tail; recovering it needs the shifted protein
//! sequence, which the description does not carry. So — like the immediate-stop
//! (`fsTer1`) and initiator-loss rules — this is a hard rejection, not a
//! rewrite.
//!
//! # Scope
//!
//! Only a frameshift that **states** its new anchor residue and makes it equal
//! the reference is rejected. The short form `p.Arg97fs` (no stated residue)
//! and any genuine change (`p.Gln151ThrfsTer9`) are untouched.

use ferro_hgvs::parse_hgvs;

#[test]
fn rejects_unchanged_anchor_frameshift() {
    // `His150` → `His150` is no change; the anchor must be the first *changed*
    // residue (frameshift.md:47-49). The `Ter#` spelling currently parses and
    // must now be rejected.
    for input in [
        "NP_003997.1:p.His150HisfsTer10",
        "NP_003997.1:p.His150HisfsTer?",
        // 1-letter spelling: `AminoAcid` is a normalized enum, so `H` == `His`.
        "NP_003997.1:p.H150HfsTer10",
        // Predicted `( )` wraps the edit; the anchor is still unchanged.
        "NP_003997.1:p.(His150HisfsTer10)",
        // `FrameshiftAlternatives` where *every* alternative equals the
        // reference — no reading is a change, so it is an unchanged anchor.
        "NP_003997.1:p.His150(His^His)fsTer23",
    ] {
        assert!(
            parse_hgvs(input).is_err(),
            "frameshift.md:47-49 forbids the unchanged anchor in {input:?}"
        );
    }
}

#[test]
fn rejects_uncertain_position_unchanged_anchor() {
    // The anchor validator must read the position through `.inner()`, not gate
    // on `Mu::Certain`, so an *uncertain-position* unchanged anchor is caught
    // too rather than silently slipping through.
    assert!(
        parse_hgvs("NP_003997.1:p.(His150)HisfsTer10").is_err(),
        "an uncertain-position unchanged anchor is still spec-invalid"
    );
}

#[test]
fn rejection_names_the_spec_rule() {
    let msg = parse_hgvs("NP_003997.1:p.His150HisfsTer10")
        .unwrap_err()
        .to_string();
    assert!(
        msg.contains("frameshift.md:47-49"),
        "the diagnostic should cite the frameshift-anchor rule; got: {msg}"
    );
}

#[test]
fn changed_anchor_frameshift_still_parses() {
    // A genuine change at the anchor (`Gln151` → `Thr`) is a legal frameshift.
    for input in [
        "NP_003997.1:p.Gln151ThrfsTer9",
        "NP_003997.1:p.Arg97ProfsTer23",
        "NP_003997.1:p.Tyr4ValfsTer2",
        // `FrameshiftAlternatives` with at least one changed reading is a valid
        // frameshift — only an all-unchanged set is rejected.
        "NP_003997.1:p.His150(His^Ala)fsTer23",
        "NP_003997.1:p.Gly719(Ala^Ser)fsTer23",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "{input:?} changes the anchor residue and must parse"
        );
    }
}

#[test]
fn short_form_frameshift_still_parses() {
    // The short form names no new residue, so it makes no unchanged-anchor
    // claim and must not be rejected (frameshift.md short format).
    for input in ["NP_003997.1:p.Arg97fs", "NP_003997.1:p.His150fs"] {
        assert!(
            parse_hgvs(input).is_ok(),
            "{input:?} states no anchor residue and must parse"
        );
    }
}
