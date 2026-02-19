//! Fuzz target for the main HGVS parser entry point
//!
//! This target feeds arbitrary byte strings to the parser to find crashes,
//! panics, or memory issues.

#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Try to convert bytes to a string
    if let Ok(input) = std::str::from_utf8(data) {
        // Don't fuzz extremely long inputs - diminishing returns
        if input.len() > 1000 {
            return;
        }

        // The parser should never panic or crash on any input
        let _ = ferro_hgvs::parse_hgvs(input);
    }
});
