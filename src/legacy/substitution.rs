//! Old substitution syntax parser.
//!
//! Handles deprecated substitution formats:
//! - `c.100_102>ATG` → `c.100_102delinsATG` (multi-base substitution with arrow)
//! - `c.100insA` → `c.100_101insA` (single position insertion)

/// Old substitution format detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OldSubstitutionFormat {
    /// Arrow instead of delins: `c.100_102>ATG`
    ArrowSubstitution,
    /// Single position insertion: `c.100insA`
    SinglePositionInsertion,
}

/// Parse old substitution format.
pub fn parse_old_substitution(input: &str) -> Option<(OldSubstitutionFormat, String)> {
    convert_old_substitution(input).map(|s| (OldSubstitutionFormat::ArrowSubstitution, s))
}

/// Convert old substitution syntax to modern format.
///
/// Handles:
/// - `c.100_102>ATG` → `c.100_102delinsATG`
/// - `g.12345_12350>ACGTAC` → `g.12345_12350delinsACGTAC`
pub fn convert_old_substitution(input: &str) -> Option<String> {
    // Look for pattern: prefix + interval + > + sequence
    // Must have underscore (interval) followed by > (not => or other)

    // Find the coordinate type prefix (c., g., n., etc.)
    let coord_start = input.find(['c', 'g', 'n', 'm', 'o'])?;
    let after_type = &input[coord_start..];

    // Must have a . after the type
    if after_type.chars().nth(1).is_none_or(|c| c != '.') {
        return None;
    }

    // Look for interval pattern: number_number
    let underscore_pos = after_type.find('_')?;
    let arrow_pos = after_type.find('>')?;

    // Arrow must come after underscore
    if arrow_pos <= underscore_pos {
        return None;
    }

    // Check that what's between underscore and arrow is a number (possibly with offset)
    let between = &after_type[underscore_pos + 1..arrow_pos];
    if between.is_empty() || !between.chars().next()?.is_ascii_digit() {
        return None;
    }

    // Check it's not already a delins (>del or similar)
    if after_type[..arrow_pos].contains("del") {
        return None;
    }

    // Get the sequence after >
    let sequence = &after_type[arrow_pos + 1..];
    if sequence.is_empty() || !sequence.chars().all(|c| "ACGTacgtNn".contains(c)) {
        return None;
    }

    // Build the converted string
    let prefix = &input[..coord_start];
    let before_arrow = &after_type[..arrow_pos];

    Some(format!("{}{}delins{}", prefix, before_arrow, sequence))
}

/// Convert single position insertion to modern two-position format.
///
/// `c.100insA` → `c.100_101insA`
#[allow(dead_code)]
pub fn convert_single_position_insertion(input: &str) -> Option<String> {
    // Look for pattern: prefix + position + ins + sequence (no underscore in position)

    let ins_pos = input.find("ins")?;
    let before_ins = &input[..ins_pos];

    // Find the coordinate type marker (c., g., n., etc.)
    let coord_marker = before_ins.rfind(['c', 'g', 'n', 'm', 'o'])?;

    // Check the character after the coordinate type is a period
    if before_ins.get(coord_marker + 1..coord_marker + 2) != Some(".") {
        return None;
    }

    // Get the position part (after "c." or similar)
    let pos_part = &before_ins[coord_marker + 2..];

    // Check there's no underscore in the position part (which would indicate interval)
    if pos_part.contains('_') {
        return None; // Already has two positions
    }

    // Find the position number
    let pos_start = before_ins.rfind(|c: char| !c.is_ascii_digit())?;
    let position: u64 = before_ins[pos_start + 1..].parse().ok()?;

    // Build the converted string
    let prefix = &before_ins[..=pos_start];
    let after_ins = &input[ins_pos..];

    Some(format!(
        "{}{}_{}{}",
        prefix,
        position,
        position + 1,
        after_ins
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_old_substitution_simple() {
        let result = convert_old_substitution("NM_000088.3:c.100_102>ATG");
        assert_eq!(result, Some("NM_000088.3:c.100_102delinsATG".to_string()));
    }

    #[test]
    fn test_old_substitution_genomic() {
        let result = convert_old_substitution("NC_000001.11:g.12345_12350>ACGTAC");
        assert_eq!(
            result,
            Some("NC_000001.11:g.12345_12350delinsACGTAC".to_string())
        );
    }

    #[test]
    fn test_old_substitution_lowercase() {
        let result = convert_old_substitution("NM_000088.3:c.100_102>atg");
        assert_eq!(result, Some("NM_000088.3:c.100_102delinsatg".to_string()));
    }

    #[test]
    fn test_old_substitution_not_interval() {
        // Single position with > is still a normal substitution
        let result = convert_old_substitution("NM_000088.3:c.100A>G");
        assert_eq!(result, None);
    }

    #[test]
    fn test_old_substitution_already_delins() {
        // Already correct format
        let result = convert_old_substitution("NM_000088.3:c.100_102delinsATG");
        assert_eq!(result, None);
    }

    #[test]
    fn test_single_position_insertion() {
        let result = convert_single_position_insertion("NM_000088.3:c.100insA");
        assert_eq!(result, Some("NM_000088.3:c.100_101insA".to_string()));
    }

    #[test]
    fn test_single_position_insertion_multi() {
        let result = convert_single_position_insertion("c.500insACGT");
        assert_eq!(result, Some("c.500_501insACGT".to_string()));
    }

    #[test]
    fn test_already_two_position_insertion() {
        // Already has two positions
        let result = convert_single_position_insertion("c.100_101insA");
        assert_eq!(result, None);
    }
}
