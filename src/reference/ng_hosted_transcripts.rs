//! Per-`NG_`-version hosted-transcript map for NG_-parent-aware legacy
//! `GENE_v001` selector resolution (#792). Built from the `NG_` GenBank mRNA
//! features (gene + transcript_id); a selector resolves only when the exact
//! parent version hosts exactly one transcript for the gene.

use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NgHostedTranscripts {
    schema_version: u32,
    /// NG_ accession.version → gene (upper) → hosted transcript ids (feature order, deduped).
    records: FxHashMap<String, FxHashMap<String, Vec<String>>>,
}

impl NgHostedTranscripts {
    pub fn from_records(
        records: impl IntoIterator<Item = (String, Vec<(String, String)>)>,
    ) -> Self {
        let mut out: FxHashMap<String, FxHashMap<String, Vec<String>>> = FxHashMap::default();
        for (ng, pairs) in records {
            let per_gene = out.entry(ng).or_default();
            for (gene, tx) in pairs {
                let list = per_gene.entry(gene).or_default();
                if !list.contains(&tx) {
                    list.push(tx);
                }
            }
        }
        Self {
            schema_version: 1,
            records: out,
        }
    }

    /// The single transcript `ng_acc_version` hosts for `gene_upper`, or `None`
    /// when absent or ambiguous (≥2 hosted).
    pub fn hosted_unique(&self, ng_acc_version: &str, gene_upper: &str) -> Option<&str> {
        let list = self.records.get(ng_acc_version)?.get(gene_upper)?;
        match list.as_slice() {
            [one] => Some(one.as_str()),
            _ => None,
        }
    }

    /// The single transcript this NG_ hosts across ALL its genes, or None when
    /// absent, multi-gene, or the sole gene hosts >1 transcript. For synthesizing
    /// the selector on a bare-NG_ c. input where the mapping is unambiguous (#923).
    pub fn sole_hosted(&self, ng_acc_version: &str) -> Option<&str> {
        let per_gene = self.records.get(ng_acc_version)?;
        let mut genes = per_gene.values();
        let only_gene = genes.next()?;
        if genes.next().is_some() {
            return None;
        }
        match only_gene.as_slice() {
            [one] => Some(one.as_str()),
            _ => None,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn to_json(&self) -> serde_json::Result<String> {
        serde_json::to_string_pretty(self)
    }

    pub fn from_json(s: &str) -> serde_json::Result<Self> {
        let v: Self = serde_json::from_str(s)?;
        if v.schema_version != 1 {
            return Err(serde::de::Error::custom(format!(
                "unsupported ng_hosted_transcripts schema_version {}",
                v.schema_version
            )));
        }
        Ok(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hosted_unique_resolves_single_and_declines_multi_and_missing() {
        let nh = NgHostedTranscripts::from_records([
            (
                "NG_012337.1".to_string(),
                vec![("TIMM8B".to_string(), "NM_012459.2".to_string())],
            ),
            (
                "NG_000000.1".to_string(),
                vec![
                    ("MULTI".to_string(), "NM_1.1".to_string()),
                    ("MULTI".to_string(), "NM_2.1".to_string()),
                ],
            ),
        ]);
        // exact version + single hosted → resolves
        assert_eq!(
            nh.hosted_unique("NG_012337.1", "TIMM8B"),
            Some("NM_012459.2")
        );
        // multi-hosted → declines
        assert_eq!(nh.hosted_unique("NG_000000.1", "MULTI"), None);
        // wrong version / unknown gene → declines
        assert_eq!(nh.hosted_unique("NG_012337.3", "TIMM8B"), None);
        assert_eq!(nh.hosted_unique("NG_012337.1", "SDHD"), None);
    }

    #[test]
    fn sole_hosted_resolves_single_gene_single_tx_and_declines_otherwise() {
        let nh = NgHostedTranscripts::from_records([
            // Single gene hosting exactly one transcript → resolves.
            (
                "NG_008939.1".to_string(),
                vec![("SLC26A4".to_string(), "NM_000441.1".to_string())],
            ),
            // Two genes → ambiguous → None.
            (
                "NG_MULTI_GENE.1".to_string(),
                vec![
                    ("GENEA".to_string(), "NM_1.1".to_string()),
                    ("GENEB".to_string(), "NM_2.1".to_string()),
                ],
            ),
            // Sole gene hosting two transcripts → ambiguous → None.
            (
                "NG_MULTI_TX.1".to_string(),
                vec![
                    ("GENEC".to_string(), "NM_3.1".to_string()),
                    ("GENEC".to_string(), "NM_4.1".to_string()),
                ],
            ),
        ]);
        assert_eq!(nh.sole_hosted("NG_008939.1"), Some("NM_000441.1"));
        assert_eq!(nh.sole_hosted("NG_MULTI_GENE.1"), None);
        assert_eq!(nh.sole_hosted("NG_MULTI_TX.1"), None);
        // Absent NG_ → None.
        assert_eq!(nh.sole_hosted("NG_999999.9"), None);
    }

    #[test]
    fn json_round_trips() {
        let nh = NgHostedTranscripts::from_records([(
            "NG_012337.1".to_string(),
            vec![("TIMM8B".to_string(), "NM_012459.2".to_string())],
        )]);
        let back = NgHostedTranscripts::from_json(&nh.to_json().unwrap()).unwrap();
        assert_eq!(
            back.hosted_unique("NG_012337.1", "TIMM8B"),
            Some("NM_012459.2")
        );
    }

    #[test]
    fn from_json_rejects_unknown_schema_version() {
        let json = r#"{"schema_version":2,"records":{}}"#;
        let err = NgHostedTranscripts::from_json(json).unwrap_err();
        assert!(
            err.to_string()
                .contains("unsupported ng_hosted_transcripts schema_version 2"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn from_json_accepts_schema_version_1() {
        let json = r#"{"schema_version":1,"records":{"NG_012337.1":{"TIMM8B":["NM_012459.2"]}}}"#;
        let nh = NgHostedTranscripts::from_json(json).unwrap();
        assert_eq!(
            nh.hosted_unique("NG_012337.1", "TIMM8B"),
            Some("NM_012459.2")
        );
    }
}
