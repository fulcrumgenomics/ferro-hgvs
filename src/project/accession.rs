//! Shared `Accession` parsing for transcript and protein identifiers.

use crate::hgvs::variant::Accession;

/// Parse an HGVS-style accession string like "NM_000088.3" or "ENSP00000256509.1".
///
/// Accepts:
/// - "PREFIX_NUMBER.VERSION" (e.g. `NM_000088.3` → prefix `NM`, number `000088`, version `Some(3)`)
/// - "PREFIX_NUMBER" (no version)
/// - Ensembl-style "PREFIXNUMBER.VERSION" or "PREFIXNUMBER" with no underscore
///   (e.g. `ENSP00000256509.1` → prefix `ENSP`, number `00000256509`, version `Some(1)`,
///   `ensembl_style = true` so `Display` does not synthesize an underscore).
/// - Anything else: pass through as the prefix with empty number.
pub(crate) fn parse_accession(s: &str) -> Accession {
    if let Some((prefix_num, version)) = s.rsplit_once('.') {
        if let Ok(v) = version.parse::<u32>() {
            if let Some((prefix, number)) = prefix_num.split_once('_') {
                return Accession::new(prefix, number, Some(v));
            }
            if let Some((prefix, number)) = split_alpha_digits(prefix_num) {
                return Accession::with_style(prefix, number, Some(v), true);
            }
            // No underscore and no alpha→digit split: treat the entire
            // prefix_num as a no-underscore prefix so Display emits
            // `<prefix_num>.<v>` without synthesizing a trailing underscore.
            return Accession::with_style(prefix_num, "", Some(v), true);
        }
    }
    if let Some((prefix, number)) = s.split_once('_') {
        return Accession::new(prefix, number, None);
    }
    if let Some((prefix, number)) = split_alpha_digits(s) {
        return Accession::with_style(prefix, number, None, true);
    }
    // Same rationale as above: no underscore and no alpha→digit split, so
    // use `ensembl_style = true` to avoid a `<s>_` Display rendering.
    Accession::with_style(s, "", None, true)
}

/// Split an accession at the first ASCII digit, returning `(alpha_prefix,
/// digits_and_rest)`. Returns `None` if the input does not start with an
/// alphabetic character followed by at least one digit.
///
/// This is the no-underscore (Ensembl-style) split path; it avoids the
/// `format!("{}_{}", prefix, number)` rendering that `Accession::full` uses
/// when `ensembl_style = false`.
fn split_alpha_digits(s: &str) -> Option<(&str, &str)> {
    let digit_start = s.find(|c: char| c.is_ascii_digit())?;
    if digit_start == 0 {
        return None;
    }
    let prefix = &s[..digit_start];
    if !prefix.chars().all(|c| c.is_ascii_alphabetic()) {
        return None;
    }
    Some((prefix, &s[digit_start..]))
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

    #[test]
    fn split_alpha_digits_requires_ascii_alpha_prefix() {
        // Strictly alphabetic prefix is accepted (Ensembl-style path).
        assert_eq!(
            split_alpha_digits("ENSP00000256509"),
            Some(("ENSP", "00000256509"))
        );
        // Non-alphabetic characters before the first digit are rejected so we
        // do not misclassify arbitrary IDs as Ensembl-style.
        assert_eq!(split_alpha_digits("A B1"), None);
        assert_eq!(split_alpha_digits("A-B1"), None);
        // A digit at position 0 is already rejected.
        assert_eq!(split_alpha_digits("1ABC"), None);
    }
}
