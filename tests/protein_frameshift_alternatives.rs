//! Frameshift with alternative new residues (#544).
//!
//! `p.Gly719(Ala^Ser)fsTer23` — a frameshift at Gly719 whose new residue is
//! Ala OR Ser (the `^` alternatives are wrapped in `(...)`), terminating at
//! Ter23. HGVS general.md.

use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn protein_frameshift_alternatives_round_trips() {
    let s = "NP_003997.1:p.Gly719(Ala^Ser)fsTer23";
    let v = parse_hgvs(s).unwrap_or_else(|e| panic!("must parse `{s}`: {e}"));
    assert!(
        matches!(v, HgvsVariant::Protein(_)),
        "expected Protein for `{s}`, got {v:?}"
    );
    assert_eq!(format!("{v}"), s, "round-trip for `{s}`");
}
