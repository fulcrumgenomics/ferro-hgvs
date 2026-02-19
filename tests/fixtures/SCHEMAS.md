# Test Fixture Schemas

This document describes the JSON schema formats used by the test fixtures in ferro-hgvs.

## Directory Structure

```
tests/fixtures/
├── grammar/            # Parsing/grammar test cases
│   ├── biocommons.json
│   └── parsing.json
├── normalization/      # Normalization test cases
│   ├── mutalyzer.json
│   └── normalization.json
├── transcripts/        # Mock transcript data
│   └── mock_transcripts.json
├── validation/         # Validation test cases
│   └── variantvalidator.json
└── SCHEMAS.md          # This file
```

---

## Transcript Schema

**File:** `transcripts/mock_transcripts.json`

Mock transcript data for unit and integration tests.

```json
[
  {
    "id": "string",              // Accession ID (e.g., "NM_000088.3")
    "gene_symbol": "string",     // Gene symbol (e.g., "COL1A1")
    "strand": "string",          // "+" or "-"
    "sequence": "string",        // Full transcript sequence (A, C, G, T)
    "cds_start": number | null,  // 1-based CDS start position (null for non-coding)
    "cds_end": number | null,    // 1-based CDS end position (null for non-coding)
    "exons": [                   // Array of exon definitions
      {
        "number": number,        // Exon number (1-based)
        "start": number,         // 1-based start position in transcript
        "end": number            // 1-based end position in transcript
      }
    ]
  }
]
```

### Example

```json
{
  "id": "NM_000088.3",
  "gene_symbol": "COL1A1",
  "strand": "+",
  "sequence": "ATGCCCCCC...",
  "cds_start": 1,
  "cds_end": 87,
  "exons": [
    { "number": 1, "start": 1, "end": 30 },
    { "number": 2, "start": 31, "end": 60 }
  ]
}
```

---

## Parsing Test Schema

**File:** `grammar/parsing.json`

Test cases for HGVS variant parsing.

```json
{
  "parsing": [
    {
      "input": "string",         // HGVS variant string to parse
      "valid": boolean,          // Whether parsing should succeed
      "type": "string",          // Expected variant type (e.g., "GenomeVariant", "CdsVariant")
      "description": "string"    // Human-readable test description
    }
  ]
}
```

### Variant Types

| Type | Description |
|------|-------------|
| `GenomeVariant` | Genomic variant (g.) |
| `CdsVariant` | Coding DNA variant (c.) |
| `TxVariant` | Non-coding transcript variant (n.) |
| `RnaVariant` | RNA variant (r.) |
| `ProteinVariant` | Protein variant (p.) |
| `MtVariant` | Mitochondrial variant (m.) |

---

## Biocommons Test Schema

**File:** `grammar/biocommons.json`

Test cases from the biocommons/hgvs repository.

```json
{
  "description": "string",       // File description
  "source": "string",            // Source URL
  "tests": [
    {
      "input": "string",         // HGVS variant string
      "valid": boolean,          // Whether parsing should succeed
      "type": "string",          // Expected type (abbreviated: "Cds", "Genome", etc.)
      "source": "string"         // Test source identifier
    }
  ]
}
```

### Type Abbreviations (Biocommons)

| Abbreviation | Full Type |
|--------------|-----------|
| `Cds` | CdsVariant |
| `Genome` | GenomeVariant |
| `Tx` | TxVariant |
| `Rna` | RnaVariant |
| `Protein` | ProteinVariant |
| `Mt` | MtVariant |

---

## Mutalyzer Test Schema

**File:** `normalization/mutalyzer.json`

Test cases from the mutalyzer/hgvs-parser repository.

```json
{
  "description": "string",
  "source": "string",
  "parsing_tests": [
    {
      "input": "string",         // HGVS variant string
      "valid": boolean,          // Whether parsing should succeed
      "type": "string",          // Expected type (abbreviated)
      "supported": boolean,      // Whether this syntax is supported by ferro-hgvs
      "source": "string"         // Source identifier
    }
  ],
  "normalization_tests": [       // Optional: normalization tests
    {
      "input": "string",         // Input HGVS variant
      "normalized": "string",    // Expected normalized form
      "description": "string"    // Test description
    }
  ]
}
```

---

## VariantValidator Test Schema

**File:** `validation/variantvalidator.json`

Test cases from VariantValidator.

```json
{
  "description": "string",
  "source": "string",
  "note": "string",              // Additional notes
  "test_cases": [
    {
      "input": "string",         // HGVS variant string
      "type": "string",          // Expected type (abbreviated)
      "valid": boolean,          // Whether parsing should succeed
      "description": "string"    // Human-readable description (e.g., gene name, variant effect)
    }
  ]
}
```

---

## Usage Notes

### Adding New Test Cases

1. Add test cases to the appropriate fixture file
2. Ensure the schema matches this documentation
3. Run tests to verify: `cargo test`

### Type Consistency

Some fixture files use abbreviated type names (`Cds`, `Genome`) while others use full names (`CdsVariant`, `GenomeVariant`). The test code handles both formats.

### Non-coding Transcripts

For non-coding transcripts (e.g., NR accessions), set `cds_start` and `cds_end` to `null`.

### Validation

- All positions are 1-based (matching HGVS convention)
- Sequences should only contain valid nucleotide characters: A, C, G, T (or U for RNA)
- Exons should be non-overlapping and cover the full transcript
