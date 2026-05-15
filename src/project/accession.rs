//! Shared `Accession` parsing for transcript and protein identifiers.

use crate::hgvs::variant::Accession;

/// Parse an HGVS-style accession string like "NM_000088.3" or "NP_000079.2".
///
/// Accepts:
/// - "PREFIX_NUMBER.VERSION" (e.g. `NM_000088.3` → prefix `NM`, number `000088`, version `Some(3)`)
/// - "PREFIX_NUMBER" (no version)
/// - Anything else (e.g. Ensembl-style `ENSP00000256509`): pass through as the prefix with empty number.
pub(crate) fn parse_accession(s: &str) -> Accession {
    if let Some((prefix_num, version)) = s.rsplit_once('.') {
        if let Ok(v) = version.parse::<u32>() {
            if let Some((prefix, number)) = prefix_num.split_once('_') {
                return Accession::new(prefix, number, Some(v));
            }
        }
    }
    if let Some((prefix, number)) = s.split_once('_') {
        return Accession::new(prefix, number, None);
    }
    Accession::new(s, "", None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_nm_with_version() {
        let acc = parse_accession("NM_000088.3");
        assert_eq!(acc.to_string(), "NM_000088.3");
    }

    #[test]
    fn parse_np_with_version() {
        let acc = parse_accession("NP_000079.2");
        assert_eq!(acc.to_string(), "NP_000079.2");
    }

    #[test]
    fn parse_without_version() {
        let acc = parse_accession("NM_000088");
        assert_eq!(acc.to_string(), "NM_000088");
    }
}
