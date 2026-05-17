//! Audit for issue #228 — `src/vcf/annotate.rs::determine_consequence`
//! misclassifies C-terminal protein extensions as `Nonsense` because
//! the protein-edit dispatcher uses `hgvs.contains("Ter")` as the
//! Nonsense heuristic, which fires on the `Ter315` start position of a
//! stop-loss extension (and, post-#225, also on the `extTer10` count).
//!
//! Spec: an extension at the stop codon is a `stop_lost`
//! (`Consequence::StopLoss`), not a `stop_gained` (`Consequence::Nonsense`).
//! See SO:0001587 (stop_gained) vs SO:0001578 (stop_lost).
//!
//! The principled fix is to drive the classification off the
//! `ProteinEdit` enum variants directly, not off substring matches on
//! the `Debug`-formatted edit or the HGVS string.

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::transcript::ManeStatus;
use ferro_hgvs::vcf::{determine_consequence, Consequence, HgvsAnnotation};

fn protein_annotation(input: &str) -> HgvsAnnotation {
    let variant =
        parse_hgvs(input).unwrap_or_else(|e| panic!("parse failed for `{}`: {}", input, e));
    HgvsAnnotation {
        variant,
        gene_symbol: None,
        transcript_accession: Some("NP_TEST".to_string()),
        is_coding: true,
        is_intronic: false,
        mane_status: ManeStatus::None,
        intron_number: None,
        intron_position: None,
        intronic_consequence: None,
    }
}

// =============================================================================
// SECTION 1 — Extension at the stop codon → StopLoss (the original bug)
// =============================================================================

mod extension_is_stop_loss {
    use super::*;

    #[test]
    fn c_terminal_extension_with_ter_canonical_classifies_as_stop_loss() {
        let ann = protein_annotation("NP_003997.1:p.Ter315TyrextTer10");
        assert_eq!(
            determine_consequence(&ann),
            Consequence::StopLoss,
            "C-terminal extension must classify as StopLoss, not Nonsense",
        );
    }

    /// Same shape, star canonical form (accepted on input, even though
    /// Display canonicalizes to `Ter`). The dispatcher must classify by
    /// the parsed enum, not the raw HGVS text.
    #[test]
    fn c_terminal_extension_with_star_input_classifies_as_stop_loss() {
        let ann = protein_annotation("NP_003997.1:p.Ter315Tyrext*10");
        assert_eq!(determine_consequence(&ann), Consequence::StopLoss);
    }

    /// Unknown extension count `Ter?` — classify by enum variant, not
    /// by whether the HGVS contains `Ter`.
    #[test]
    fn c_terminal_extension_unknown_count_classifies_as_stop_loss() {
        let ann = protein_annotation("NP_003997.1:p.Ter315TyrextTer?");
        assert_eq!(determine_consequence(&ann), Consequence::StopLoss);
    }

    /// N-terminal extension (Met loss) — start position is `Met1`, no
    /// `Ter` substring at all. Still must classify as StopLoss because
    /// Extension semantics are about gain-of-translated-sequence on
    /// either terminus.
    #[test]
    fn n_terminal_extension_classifies_as_stop_loss() {
        let ann = protein_annotation("NP_003997.1:p.Met1Valext-12");
        assert_eq!(determine_consequence(&ann), Consequence::StopLoss);
    }

    /// Predicted parens around extension — same classification.
    #[test]
    fn predicted_c_terminal_extension_classifies_as_stop_loss() {
        let ann = protein_annotation("NP_003997.1:p.(Ter315TyrextTer10)");
        assert_eq!(determine_consequence(&ann), Consequence::StopLoss);
    }
}

// =============================================================================
// SECTION 2 — Nonsense classification regression guards
// =============================================================================
//
// The fix must NOT regress nonsense classification — substitutions to
// `Ter` are the canonical nonsense case and must still surface as
// `Consequence::Nonsense`.

mod nonsense_classification {
    use super::*;

    #[test]
    fn substitution_to_ter_classifies_as_nonsense() {
        let ann = protein_annotation("NP_003997.1:p.Tyr4Ter");
        assert_eq!(determine_consequence(&ann), Consequence::Nonsense);
    }

    /// Star-form input — same classification.
    #[test]
    fn substitution_to_star_input_classifies_as_nonsense() {
        let ann = protein_annotation("NP_003997.1:p.Tyr4*");
        assert_eq!(determine_consequence(&ann), Consequence::Nonsense);
    }

    /// One-letter star input — same classification.
    #[test]
    fn substitution_one_letter_star_classifies_as_nonsense() {
        let ann = protein_annotation("NP_003997.1:p.Y4*");
        assert_eq!(determine_consequence(&ann), Consequence::Nonsense);
    }

    /// Predicted parens around nonsense — same classification.
    #[test]
    fn predicted_substitution_to_ter_classifies_as_nonsense() {
        let ann = protein_annotation("NP_003997.1:p.(Tyr4Ter)");
        assert_eq!(determine_consequence(&ann), Consequence::Nonsense);
    }
}

// =============================================================================
// SECTION 3 — Frameshift classification regression guards
// =============================================================================
//
// Frameshift is checked first today; the rewrite must preserve that.
// Frameshifts whose Display contains `Ter` (`p.Arg97ProfsTer23`) must
// not silently collapse to Nonsense.

mod frameshift_classification {
    use super::*;

    #[test]
    fn frameshift_with_ter_in_display_classifies_as_frameshift() {
        let ann = protein_annotation("NP_003997.1:p.Arg97ProfsTer23");
        assert_eq!(determine_consequence(&ann), Consequence::Frameshift);
    }

    #[test]
    fn frameshift_short_form_classifies_as_frameshift() {
        let ann = protein_annotation("NP_003997.1:p.Arg97fs");
        assert_eq!(determine_consequence(&ann), Consequence::Frameshift);
    }

    #[test]
    fn predicted_frameshift_classifies_as_frameshift() {
        let ann = protein_annotation("NP_003997.1:p.(Arg97ProfsTer23)");
        assert_eq!(determine_consequence(&ann), Consequence::Frameshift);
    }
}

// =============================================================================
// SECTION 4 — Other protein edit kinds (regression guards)
// =============================================================================

mod other_edit_kinds {
    use super::*;

    /// Missense substitution (non-Ter alternative).
    #[test]
    fn non_ter_substitution_classifies_as_missense() {
        let ann = protein_annotation("NP_003997.1:p.Arg97Trp");
        assert_eq!(determine_consequence(&ann), Consequence::Missense);
    }

    #[test]
    fn deletion_classifies_as_inframe_deletion() {
        let ann = protein_annotation("NP_003997.1:p.Lys100del");
        assert_eq!(determine_consequence(&ann), Consequence::InframeDeletion);
    }

    #[test]
    fn insertion_classifies_as_inframe_insertion() {
        let ann = protein_annotation("NP_003997.1:p.Lys100_Arg101insGlySer");
        assert_eq!(determine_consequence(&ann), Consequence::InframeInsertion);
    }
}
