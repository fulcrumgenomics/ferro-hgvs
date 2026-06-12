//! Golden/contract tests for the curated tool-support matrix.
#![cfg(feature = "dev")]

use ferro_hgvs::tool_support::{Matrix, Normalize, Validate};

fn matrix() -> Matrix {
    Matrix::load("docs/tool_support_matrix.json").expect("load matrix")
}

#[test]
fn matrix_invariants_hold() {
    matrix().validate_invariants().expect("invariants");
}

#[test]
fn adjudicated_cells_are_correct() {
    let m = matrix();
    let p = m.find_row("coordinate_types", "p").unwrap();
    assert_eq!(p.support["ferro"].normalize, Normalize::Full);
    assert_eq!(p.support["hgvs-rs"].normalize, Normalize::Errors);
    assert_eq!(p.support["hgvs-rs"].validate, Validate::Permissive);
    assert_eq!(
        p.support["biocommons"].normalize,
        Normalize::UnsupportedByDesign
    );
    assert_eq!(p.support["mutalyzer"].normalize, Normalize::NetworkOnly);

    let intronic = m.find_row("coordinate_types", "intronic").unwrap();
    assert_eq!(intronic.support["ferro"].normalize, Normalize::Full);
    assert_eq!(intronic.support["biocommons"].normalize, Normalize::No);
    assert_eq!(
        intronic.support["hgvs-rs"].normalize,
        Normalize::UnsupportedByDesign
    );

    let r = m.find_row("coordinate_types", "r").unwrap();
    assert_eq!(r.support["hgvs-rs"].normalize, Normalize::Full);

    let conv = m.find_row("variant_types", "conversion").unwrap();
    assert_eq!(conv.support["ferro"].normalize, Normalize::Full);
}

#[test]
fn normalization_view_renders_expected() {
    let m = matrix();
    let out = m
        .render_markdown_view("normalization_capabilities")
        .unwrap();
    let expected = "\
| Pattern Type | ferro | mutalyzer | biocommons | hgvs-rs |
|--------------|:-----:|:---------:|:----------:|:-------:|
| Genomic (g.) | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) exonic | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) intronic | ✓ | ✓** | ✗ | ✗ |
| Non-coding (n.) | ✓ | ✓ | ✓ | ✓ |
| RNA (r.) | ✓ | ✓ | ✓ | ✓ |
| Protein (p.) | ✓ | Net* | ✗ | ✗ |
";
    // Compare the table portion (ignore the trailing footnote legend, which is
    // asserted separately to avoid coupling to exact footnote prose).
    let table: String = out.lines().take(8).collect::<Vec<_>>().join("\n");
    assert_eq!(table, expected.trim_end());
    // Canonical marker order (sorted keys: `net` < `rewrite`): * = net, ** = rewrite.
    assert!(out.contains("\\* mutalyzer protein normalization requires network access"));
    assert!(out.contains("\\*\\* mutalyzer intronic support is enabled by default"));
}

#[test]
fn newly_verified_cells() {
    let m = matrix();

    let inversion = m.find_row("variant_types", "inversion").unwrap();
    assert_eq!(inversion.support["hgvs-rs"].normalize, Normalize::Full);

    let conversion = m.find_row("variant_types", "conversion").unwrap();
    assert_eq!(
        conversion.support["mutalyzer"].normalize,
        Normalize::UnsupportedByDesign
    );

    let repeat = m.find_row("variant_types", "repeat").unwrap();
    assert_eq!(repeat.support["biocommons"].validate, Validate::No);

    let enst = m.find_row("reference_types", "ENST").unwrap();
    assert_eq!(enst.support["ferro"].validate, Validate::Yes);
    assert_eq!(enst.support["mutalyzer"].normalize, Normalize::Full);
    assert_eq!(enst.support["biocommons"].normalize, Normalize::No);
}
