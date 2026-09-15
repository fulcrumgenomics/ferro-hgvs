//! A Mutalyzer-style, multi-axis view of a variant — computed entirely by ferro.
//!
//! Mutalyzer's `/normalize` returns, for one input, a normalized description plus
//! the equivalent descriptions across coordinate systems (and any messages). This
//! module bundles ferro's own [`normalize`](crate::Normalizer::normalize) output
//! and its [`VariantProjection`] axes into that familiar shape so callers used to
//! Mutalyzer can get the whole picture from one call.
//!
//! **Interface-compatible, not output-compatible.** The *values* here are ferro's
//! own spec-grounded normalization and projection, so they can differ from what
//! Mutalyzer would emit — that is expected, not a bug. Only the *shape* (one
//! object carrying the normalized form + the equivalent axes + warnings) is
//! Mutalyzer-familiar. Mutalyzer is not a spec oracle for ferro; see the
//! `adjudication-precedence-order` ruling.

use std::fmt;

use crate::project::result::VariantProjection;

/// A warning attached to a [`MutalyzerResult`]: ferro's own code plus a
/// human-readable message.
///
/// The `code` is ferro's warning-code vocabulary (e.g. `REFSEQ_MISMATCH`), not
/// Mutalyzer's — a compat *interface* passes ferro's diagnostics through rather
/// than inventing a mapping to Mutalyzer's code set.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct MutalyzerWarning {
    pub code: String,
    pub message: String,
}

/// A variant resolved across coordinate axes, presented in a Mutalyzer-familiar
/// shape and computed entirely by ferro.
///
/// Each axis is `Option<String>` — the rendered HGVS description on that axis, or
/// `None` when it is not derivable from the input and available reference data
/// (never an error; see [`VariantProjection`]). `normalized` is the input's own
/// axis, normalized; it is set by the constructor that also runs normalization
/// (see [`crate::project::VariantProjector`]'s `mutalyzer` method) and is left
/// `None` by [`MutalyzerResult::from_projection`] alone.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct MutalyzerResult {
    /// The input description, echoed verbatim.
    pub input: String,
    /// The input normalized on its own coordinate axis, when available.
    pub normalized: Option<String>,
    /// The genomic (`g.`) description.
    pub genomic: Option<String>,
    /// The coding (`c.`) description.
    pub coding: Option<String>,
    /// The non-coding transcript (`n.`) description.
    pub noncoding: Option<String>,
    /// The predicted RNA consequence (`r.`).
    pub rna: Option<String>,
    /// The predicted protein consequence (`p.`).
    pub protein: Option<String>,
    /// The gene symbol, when known.
    pub gene_symbol: Option<String>,
    /// The transcript the axes were resolved against.
    pub transcript_id: String,
    /// Diagnostics emitted while normalizing the input (ferro's codes).
    pub warnings: Vec<MutalyzerWarning>,
}

impl MutalyzerResult {
    /// Build a result from an input description and a [`VariantProjection`].
    ///
    /// Fills the axis fields, gene symbol, transcript id and warnings from the
    /// projection. `normalized` is left `None` — the caller that ran
    /// normalization sets it via [`MutalyzerResult::with_normalized`], since a
    /// projection does not carry the input's own-axis normalized form as a
    /// distinct field.
    pub fn from_projection(
        input: impl Into<String>,
        projection: &VariantProjection,
        protein_style: crate::hgvs::location::ProteinRenderStyle,
    ) -> Self {
        let render = |axis: &Option<crate::hgvs::variant::HgvsVariant>| {
            axis.as_ref().map(|variant| variant.to_string())
        };
        let mut warnings: Vec<MutalyzerWarning> = projection
            .normalization_warnings
            .iter()
            .map(|warning| MutalyzerWarning {
                code: warning.code().to_string(),
                message: warning.to_string(),
            })
            .collect();
        // Fold each declined axis's reason in, so `ferro mutalyzer` reports the
        // same "why is this axis `-`" diagnostics that `ferro project --axis g`
        // already surfaces for the same input. `VariantProjection` explains an
        // absent axis through `axis_decline_reasons` (e.g. a missing genomic
        // parent); without this the bundle would show a bare `-` and no reason,
        // contradicting the module doc's "(and any messages)" contract.
        for (axis, reason) in [
            ("genomic", &projection.axis_decline_reasons.genomic),
            ("coding", &projection.axis_decline_reasons.coding),
            ("noncoding", &projection.axis_decline_reasons.noncoding),
            ("rna", &projection.axis_decline_reasons.rna),
            ("protein", &projection.axis_decline_reasons.protein),
        ] {
            if let Some(reason) = reason {
                warnings.push(MutalyzerWarning {
                    code: format!("{axis}_axis_declined"),
                    message: reason.clone(),
                });
            }
        }
        MutalyzerResult {
            input: input.into(),
            normalized: None,
            genomic: render(&projection.genomic),
            coding: render(&projection.coding),
            noncoding: render(&projection.noncoding),
            rna: render(&projection.rna),
            // Protein is the one axis whose rendering depends on the caller's
            // style (`protein_stop` / `amino_acid_code`); the g./c./n./r. axes
            // have no styled form. Render it through `protein_style` so this
            // bundle honours a projector's configured style like every other
            // protein surface (#1050), rather than defaulting to Ter/three.
            protein: projection
                .protein
                .as_ref()
                .map(|variant| variant.to_styled_string(protein_style)),
            gene_symbol: projection.gene_symbol.clone(),
            transcript_id: projection.transcript_id.clone(),
            warnings,
        }
    }

    /// Set the `normalized` (input's own axis) description, returning `self` for
    /// chaining.
    #[must_use]
    pub fn with_normalized(mut self, normalized: Option<String>) -> Self {
        self.normalized = normalized;
        self
    }
}

impl fmt::Display for MutalyzerResult {
    /// A plain, aligned, Mutalyzer-familiar block. Absent axes are shown as `-`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let show = |axis: &Option<String>| axis.clone().unwrap_or_else(|| "-".to_string());
        writeln!(f, "Input:       {}", self.input)?;
        writeln!(f, "Normalized:  {}", show(&self.normalized))?;
        writeln!(f, "Genomic:     {}", show(&self.genomic))?;
        writeln!(f, "Coding:      {}", show(&self.coding))?;
        writeln!(f, "Noncoding:   {}", show(&self.noncoding))?;
        writeln!(f, "RNA:         {}", show(&self.rna))?;
        writeln!(f, "Protein:     {}", show(&self.protein))?;
        if let Some(symbol) = &self.gene_symbol {
            writeln!(f, "Gene:        {symbol}")?;
        }
        if self.warnings.is_empty() {
            write!(f, "Warnings:    (none)")
        } else {
            // No trailing newline, so callers can add a consistent record
            // separator regardless of whether warnings are present.
            write!(f, "Warnings:")?;
            for warning in &self.warnings {
                write!(f, "\n  [{}] {}", warning.code, warning.message)?;
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::project::result::{AxisDeclineReasons, VariantProjection};

    /// A projection with only the axes we set; everything else defaulted.
    fn projection_with(
        genomic: Option<&str>,
        coding: Option<&str>,
        protein: Option<&str>,
    ) -> VariantProjection {
        let parse = |s: Option<&str>| s.map(|d| crate::parse_hgvs(d).expect("valid test HGVS"));
        VariantProjection {
            genomic: parse(genomic),
            coding: parse(coding),
            noncoding: None,
            protein: parse(protein),
            rna: None,
            transcript_id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            is_frameshift: false,
            is_intronic: false,
            is_utr: false,
            affects_init: false,
            normalization_warnings: Vec::new(),
            axis_decline_reasons: AxisDeclineReasons::default(),
        }
    }

    #[test]
    fn from_projection_renders_each_present_axis_as_a_string() {
        let projection = projection_with(
            Some("NC_000017.11:g.50201978C>A"),
            Some("NM_000088.3:c.589G>T"),
            Some("NP_000079.2:p.Gly197Cys"),
        );
        let result = MutalyzerResult::from_projection(
            "NM_000088.3:c.589G>T",
            &projection,
            crate::hgvs::location::ProteinRenderStyle::default(),
        );
        assert_eq!(result.input, "NM_000088.3:c.589G>T");
        assert_eq!(
            result.genomic.as_deref(),
            Some("NC_000017.11:g.50201978C>A")
        );
        assert_eq!(result.coding.as_deref(), Some("NM_000088.3:c.589G>T"));
        assert_eq!(result.protein.as_deref(), Some("NP_000079.2:p.Gly197Cys"));
        assert_eq!(result.transcript_id, "NM_000088.3");
        assert_eq!(result.gene_symbol.as_deref(), Some("COL1A1"));
        // Not set by `from_projection` alone.
        assert_eq!(result.normalized, None);
    }

    #[test]
    fn an_absent_axis_is_none_and_renders_as_a_dash() {
        let projection = projection_with(None, Some("NM_000088.3:c.589G>T"), None);
        let result = MutalyzerResult::from_projection(
            "NM_000088.3:c.589G>T",
            &projection,
            crate::hgvs::location::ProteinRenderStyle::default(),
        );
        assert_eq!(result.genomic, None);
        assert_eq!(result.rna, None);
        assert!(result.to_string().contains("Genomic:     -"));
    }

    #[test]
    fn with_normalized_sets_the_normalized_axis() {
        let projection = projection_with(None, Some("NM_000088.3:c.589G>T"), None);
        let result = MutalyzerResult::from_projection(
            "NM_000088.3:c.589G>T",
            &projection,
            crate::hgvs::location::ProteinRenderStyle::default(),
        )
        .with_normalized(Some("NM_000088.3:c.589G>T".to_string()));
        assert_eq!(result.normalized.as_deref(), Some("NM_000088.3:c.589G>T"));
        assert!(result
            .to_string()
            .contains("Normalized:  NM_000088.3:c.589G>T"));
    }

    #[test]
    fn serializes_with_the_documented_json_field_names() {
        let projection = projection_with(
            Some("NC_000017.11:g.50201978C>A"),
            Some("NM_000088.3:c.589G>T"),
            Some("NP_000079.2:p.Gly197Cys"),
        );
        let result = MutalyzerResult::from_projection(
            "NM_000088.3:c.589G>T",
            &projection,
            crate::hgvs::location::ProteinRenderStyle::default(),
        )
        .with_normalized(Some("NM_000088.3:c.589G>T".to_string()));
        let value: serde_json::Value =
            serde_json::to_value(&result).expect("MutalyzerResult serializes");
        let object = value.as_object().expect("serializes to a JSON object");
        // The CLI's `--format json` output and the guide (docs/src/guide/
        // mutalyzer-view.md) document this schema; a renamed field silently
        // breaks downstream parsers, so pin the exact key set.
        let mut keys: Vec<&str> = object.keys().map(String::as_str).collect();
        keys.sort_unstable();
        assert_eq!(
            keys,
            [
                "coding",
                "gene_symbol",
                "genomic",
                "input",
                "noncoding",
                "normalized",
                "protein",
                "rna",
                "transcript_id",
                "warnings",
            ]
        );
    }

    #[test]
    fn axis_decline_reasons_are_folded_into_warnings() {
        // A declined axis (here genomic, e.g. no genomic parent) carries its
        // reason on the projection; the bundle must report it so `ferro
        // mutalyzer` and `ferro project --axis g` give the same diagnostic.
        let mut projection = projection_with(None, Some("NM_000088.3:c.589G>T"), None);
        projection.axis_decline_reasons.genomic =
            Some("no genomic parent placement for a bare NM_ input".to_string());
        let result = MutalyzerResult::from_projection(
            "NM_000088.3:c.589G>T",
            &projection,
            crate::hgvs::location::ProteinRenderStyle::default(),
        );
        assert_eq!(result.genomic, None);
        let declined = result
            .warnings
            .iter()
            .find(|w| w.code == "genomic_axis_declined")
            .expect("a declined genomic axis surfaces its reason as a warning");
        assert_eq!(
            declined.message,
            "no genomic parent placement for a bare NM_ input"
        );
        // The reason also shows in the text block's Warnings section.
        assert!(result.to_string().contains("genomic_axis_declined"));
    }

    #[test]
    fn protein_axis_honours_the_render_style() {
        use crate::hgvs::location::{AaCode, ProteinRenderStyle, TerStyle};
        let projection = projection_with(None, None, Some("NP_000079.2:p.Gly197Cys"));

        // Default style renders three-letter codes.
        let default = MutalyzerResult::from_projection(
            "NM_000088.3:c.589G>T",
            &projection,
            ProteinRenderStyle::default(),
        );
        assert_eq!(default.protein.as_deref(), Some("NP_000079.2:p.Gly197Cys"));

        // A one-letter style must reach the protein axis (was previously rendered
        // with the default `to_string()` regardless of the caller's style).
        let one = ProteinRenderStyle {
            stop: TerStyle::Ter,
            aa_code: AaCode::One,
        };
        let styled = MutalyzerResult::from_projection("NM_000088.3:c.589G>T", &projection, one);
        assert_eq!(styled.protein.as_deref(), Some("NP_000079.2:p.G197C"));
    }
}
