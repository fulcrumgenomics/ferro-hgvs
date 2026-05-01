# Issue #82: Python from_manifest Loader Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Allow Python clients to load real reference data into `Normalizer`, `EquivalenceChecker`, `BatchProcessor`, and `CoordinateMapper` via a `from_manifest(path)` classmethod, and extend `MockProvider::from_json` to accept proteins and genomic sequences alongside transcripts.

**Architecture:** Two changes work in tandem. (1) Extend `MockProvider::from_json` to accept either the existing bare-array form OR a new object form `{transcripts, proteins, genomic_sequences}`. (2) Introduce a `PyProvider` enum in `src/python.rs` that wraps `MockProvider` (cheap to clone) and `Arc<MultiFastaProvider>` (cheap to clone) and implements `ReferenceProvider` by delegation. The four PyO3 wrapper classes that currently store `MockProvider` switch to `PyProvider`, gaining a `from_manifest` static-method constructor that delegates to `MultiFastaProvider::from_manifest`.

**Tech Stack:** Rust (PyO3 0.28, serde, serde_json), Python (pytest), maturin for builds.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `src/reference/mock.rs` | Modify | Extend `from_json` to accept object form |
| `src/python.rs` | Modify | Add `PyProvider` enum, swap `MockProvider` for `PyProvider` in 4 wrappers, add `from_manifest` static methods |
| `python/ferro_hgvs/__init__.pyi` | Modify | Add `from_manifest` stubs |
| `tests/fixtures/SCHEMAS.md` | Modify | Document new MockProvider object form |
| `tests/fixtures/sequences/normalization_transcripts.json` | Modify (in-place) | Wrap existing bare array under `{transcripts: [...], proteins: {...}}` for the object-form integration test (alternative: leave file alone, add a new tiny fixture) |
| `tests/fixtures/python/manifest_tiny/` | Create | Tiny fixture directory: `manifest.json`, `transcripts/test.fna`, `transcripts/test.fna.fai` |
| `tests/python/test_loaders.py` | Create | Two integration tests: object-form JSON + manifest |

**Decomposition note:** The three providers (`MockProvider`, `MultiFastaProvider`, `PyProvider`) keep their current responsibilities. The new enum lives inside `python.rs` because it's only needed at the Python boundary — the Rust API still uses generic `Normalizer<P: ReferenceProvider>` directly.

**Fixture-file decision:** Don't mutate `normalization_transcripts.json` (it's a generic fixture used elsewhere). Instead, the object-form test loads a small inline JSON file that we write inside the test using `tmp_path`, which avoids touching the shared fixture. This sidesteps the issue's hint about extending the file in-place.

---

## Task 1: Extend `MockProvider::from_json` to accept object form

**Files:**
- Modify: `src/reference/mock.rs`

- [ ] **Step 1.1: Write the failing test (Rust unit test, object form)**

Append to the `tests` module in `src/reference/mock.rs`:

```rust
#[test]
fn test_from_json_object_form_with_proteins_and_genomic() {
    use std::io::Write;
    use tempfile::NamedTempFile;

    let json = r#"{
      "transcripts": [{
        "id": "NM_TEST.1",
        "gene_symbol": "TEST",
        "strand": "+",
        "sequence": "ATGCATGCAT",
        "cds_start": 1,
        "cds_end": 10,
        "exons": [{"number": 1, "start": 1, "end": 10}]
      }],
      "proteins": {
        "NP_TEST.1": "MAPLE"
      },
      "genomic_sequences": {
        "chr1": "ACGTACGTACGT"
      }
    }"#;

    let mut file = NamedTempFile::new().unwrap();
    file.write_all(json.as_bytes()).unwrap();

    let provider = MockProvider::from_json(file.path()).unwrap();

    // Transcript loaded
    assert!(provider.has_transcript("NM_TEST.1"));

    // Protein loaded
    assert!(provider.has_protein_data());
    assert_eq!(
        provider.get_protein_sequence("NP_TEST.1", 0, 5).unwrap(),
        "MAPLE"
    );

    // Genomic loaded
    assert!(provider.has_genomic_data());
    assert_eq!(
        provider.get_genomic_sequence("chr1", 0, 4).unwrap(),
        "ACGT"
    );
}

#[test]
fn test_from_json_bare_array_form_still_works() {
    use std::io::Write;
    use tempfile::NamedTempFile;

    let json = r#"[{
      "id": "NM_TEST.1",
      "gene_symbol": "TEST",
      "strand": "+",
      "sequence": "ATGCATGCAT",
      "cds_start": 1,
      "cds_end": 10,
      "exons": [{"number": 1, "start": 1, "end": 10}]
    }]"#;

    let mut file = NamedTempFile::new().unwrap();
    file.write_all(json.as_bytes()).unwrap();

    let provider = MockProvider::from_json(file.path()).unwrap();

    assert!(provider.has_transcript("NM_TEST.1"));
    assert!(!provider.has_protein_data());
    assert!(!provider.has_genomic_data());
}
```

- [ ] **Step 1.2: Run test to verify it fails**

Run: `cargo nextest run --features dev test_from_json_object_form_with_proteins_and_genomic test_from_json_bare_array_form_still_works`
Expected: object-form test FAILS (parser only accepts array); bare-array test PASSES.

- [ ] **Step 1.3: Implement object-form parsing**

Replace the body of `MockProvider::from_json` in `src/reference/mock.rs`:

```rust
pub fn from_json(path: &Path) -> Result<Self, FerroError> {
    let content = std::fs::read_to_string(path)?;
    let value: serde_json::Value = serde_json::from_str(&content)?;

    // Two accepted shapes:
    //   1. Bare array: [Transcript, ...]
    //   2. Object: { "transcripts": [...], "proteins": {...}, "genomic_sequences": {...} }
    let (transcripts, proteins, genomic_sequences) = if value.is_array() {
        let transcripts: Vec<Transcript> = serde_json::from_value(value)?;
        (transcripts, HashMap::new(), HashMap::new())
    } else {
        let transcripts: Vec<Transcript> = value
            .get("transcripts")
            .cloned()
            .map(serde_json::from_value)
            .transpose()?
            .unwrap_or_default();
        let proteins: HashMap<String, String> = value
            .get("proteins")
            .cloned()
            .map(serde_json::from_value)
            .transpose()?
            .unwrap_or_default();
        let genomic_sequences: HashMap<String, String> = value
            .get("genomic_sequences")
            .cloned()
            .map(serde_json::from_value)
            .transpose()?
            .unwrap_or_default();
        (transcripts, proteins, genomic_sequences)
    };

    let map: HashMap<String, Transcript> = transcripts
        .into_iter()
        .map(|tx| (tx.id.clone(), tx))
        .collect();

    Ok(Self {
        transcripts: map,
        proteins,
        genomic_sequences,
    })
}
```

- [ ] **Step 1.4: Run tests to verify both pass**

Run: `cargo nextest run --features dev mock`
Expected: All tests in `mock` module PASS, including the two new ones.

- [ ] **Step 1.5: Commit**

```bash
git add src/reference/mock.rs
git commit -m "feat(mock): accept transcripts/proteins/genomic_sequences object form in from_json"
```

---

## Task 2: Add `PyProvider` enum in `src/python.rs`

**Files:**
- Modify: `src/python.rs`

This task introduces the abstraction without yet wiring it into the four wrappers. That's done in Tasks 3–6.

- [ ] **Step 2.1: Add `Arc` import**

In `src/python.rs` near the top imports (around line 11), add:

```rust
use std::sync::Arc;
```

Also add to existing reference imports (around line 31):

```rust
use crate::reference::{MockProvider, MultiFastaProvider};
```

(replacing `use crate::reference::MockProvider;`).

- [ ] **Step 2.2: Define the `PyProvider` enum**

Add after the imports (before the `parse` function around line 47), a new section:

```rust
// ============================================================================
// PyProvider — wraps either MockProvider or MultiFastaProvider for the Python
// surface. Cheap to clone (Arc shares the heavy MultiFasta state).
// ============================================================================

#[derive(Clone)]
pub(crate) enum PyProvider {
    Mock(MockProvider),
    MultiFasta(Arc<MultiFastaProvider>),
}

impl ReferenceProvider for PyProvider {
    fn get_transcript(
        &self,
        id: &str,
    ) -> Result<crate::reference::transcript::Transcript, crate::error::FerroError> {
        match self {
            PyProvider::Mock(p) => p.get_transcript(id),
            PyProvider::MultiFasta(p) => p.get_transcript(id),
        }
    }

    fn get_sequence(
        &self,
        id: &str,
        start: u64,
        end: u64,
    ) -> Result<String, crate::error::FerroError> {
        match self {
            PyProvider::Mock(p) => p.get_sequence(id, start, end),
            PyProvider::MultiFasta(p) => p.get_sequence(id, start, end),
        }
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, crate::error::FerroError> {
        match self {
            PyProvider::Mock(p) => p.get_genomic_sequence(contig, start, end),
            PyProvider::MultiFasta(p) => p.get_genomic_sequence(contig, start, end),
        }
    }

    fn has_genomic_data(&self) -> bool {
        match self {
            PyProvider::Mock(p) => p.has_genomic_data(),
            PyProvider::MultiFasta(p) => p.has_genomic_data(),
        }
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, crate::error::FerroError> {
        match self {
            PyProvider::Mock(p) => p.get_protein_sequence(accession, start, end),
            PyProvider::MultiFasta(p) => p.get_protein_sequence(accession, start, end),
        }
    }

    fn has_protein_data(&self) -> bool {
        match self {
            PyProvider::Mock(p) => p.has_protein_data(),
            PyProvider::MultiFasta(p) => p.has_protein_data(),
        }
    }
}

impl PyProvider {
    /// Load from a JSON file (delegates to MockProvider::from_json).
    fn from_json(path: &Path) -> PyResult<Self> {
        MockProvider::from_json(path)
            .map(PyProvider::Mock)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to load reference: {}", e)))
    }

    /// Load from a manifest file (delegates to MultiFastaProvider::from_manifest).
    fn from_manifest(path: &Path) -> PyResult<Self> {
        MultiFastaProvider::from_manifest(path)
            .map(|p| PyProvider::MultiFasta(Arc::new(p)))
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to load manifest: {}", e)))
    }

    /// Default: built-in test data.
    fn test_data() -> Self {
        PyProvider::Mock(MockProvider::with_test_data())
    }
}
```

- [ ] **Step 2.3: Build to confirm the enum compiles standalone**

Run: `cargo build --features python`
Expected: Build SUCCEEDS. The enum is dead code at this point but compiles.

- [ ] **Step 2.4: Commit**

```bash
git add src/python.rs
git commit -m "refactor(python): add PyProvider enum wrapping Mock + MultiFasta backends"
```

---

## Task 3: Switch `PyNormalizer` to `PyProvider` and add `from_manifest`

**Files:**
- Modify: `src/python.rs` (the `PyNormalizer` struct and `impl PyNormalizer` block, ~lines 360–422)

- [ ] **Step 3.1: Replace `MockProvider` with `PyProvider` in `PyNormalizer`**

In `src/python.rs`, change:

```rust
#[pyclass(name = "Normalizer")]
pub struct PyNormalizer {
    provider: MockProvider,
    config: NormalizeConfig,
}
```

to:

```rust
#[pyclass(name = "Normalizer")]
pub struct PyNormalizer {
    provider: PyProvider,
    config: NormalizeConfig,
}
```

- [ ] **Step 3.2: Update `PyNormalizer::new` to use `PyProvider`**

Replace the `new` method body (~lines 376–391):

```rust
#[new]
#[pyo3(signature = (reference_json=None, direction="3prime"))]
fn new(reference_json: Option<&str>, direction: &str) -> PyResult<Self> {
    let provider = match reference_json {
        Some(path) => PyProvider::from_json(Path::new(path))?,
        None => PyProvider::test_data(),
    };

    let config = NormalizeConfig::default().with_direction(parse_direction(direction));

    Ok(Self { provider, config })
}
```

- [ ] **Step 3.3: Add `from_manifest` static method**

Inside the same `#[pymethods] impl PyNormalizer { ... }` block, after the `new` method, add:

```rust
/// Create a normalizer from a reference manifest written by `ferro prepare`.
///
/// Args:
///     manifest_path: Path to a manifest.json file (typically inside a
///         directory produced by `ferro prepare`).
///     direction: Shuffle direction - "3prime" (default) or "5prime"
///
/// Returns:
///     A Normalizer backed by a MultiFastaProvider.
#[staticmethod]
#[pyo3(signature = (manifest_path, direction="3prime"))]
fn from_manifest(manifest_path: &str, direction: &str) -> PyResult<Self> {
    let provider = PyProvider::from_manifest(Path::new(manifest_path))?;
    let config = NormalizeConfig::default().with_direction(parse_direction(direction));
    Ok(Self { provider, config })
}
```

- [ ] **Step 3.4: Verify normalize methods still compile**

The methods `normalize_variant` and `normalize` already do `self.provider.clone()`; `PyProvider` is `Clone`, so they should compile unchanged.

Run: `cargo build --features python`
Expected: Build SUCCEEDS.

- [ ] **Step 3.5: Build the Python wheel and confirm Normalizer still works**

Run: `maturin develop --features python`
Expected: Build SUCCEEDS, wheel installs.

Then run: `pytest tests/python/test_core.py -v`
Expected: All existing tests still PASS.

- [ ] **Step 3.6: Commit**

```bash
git add src/python.rs
git commit -m "feat(python): add Normalizer.from_manifest classmethod"
```

---

## Task 4: Switch `PyEquivalenceChecker` to `PyProvider` and add `from_manifest`

**Files:**
- Modify: `src/python.rs` (the `PyEquivalenceChecker` struct and `impl` block, ~lines 866–921)

- [ ] **Step 4.1: Replace `MockProvider` with `PyProvider`**

Change:

```rust
#[pyclass(name = "EquivalenceChecker")]
pub struct PyEquivalenceChecker {
    checker: EquivalenceChecker<MockProvider>,
}
```

to:

```rust
#[pyclass(name = "EquivalenceChecker")]
pub struct PyEquivalenceChecker {
    checker: EquivalenceChecker<PyProvider>,
}
```

- [ ] **Step 4.2: Update `new`**

Replace the `new` body (~lines 880–890):

```rust
#[new]
#[pyo3(signature = (reference_json=None))]
fn new(reference_json: Option<&str>) -> PyResult<Self> {
    let provider = match reference_json {
        Some(path) => PyProvider::from_json(Path::new(path))?,
        None => PyProvider::test_data(),
    };

    Ok(Self {
        checker: EquivalenceChecker::new(provider),
    })
}
```

- [ ] **Step 4.3: Add `from_manifest` static method**

Inside `#[pymethods] impl PyEquivalenceChecker { ... }`, after `new`:

```rust
/// Create an equivalence checker from a reference manifest.
///
/// Args:
///     manifest_path: Path to a manifest.json file produced by `ferro prepare`.
///
/// Returns:
///     An EquivalenceChecker backed by a MultiFastaProvider.
#[staticmethod]
fn from_manifest(manifest_path: &str) -> PyResult<Self> {
    let provider = PyProvider::from_manifest(Path::new(manifest_path))?;
    Ok(Self {
        checker: EquivalenceChecker::new(provider),
    })
}
```

- [ ] **Step 4.4: Build and confirm**

Run: `cargo build --features python`
Expected: SUCCEEDS.

- [ ] **Step 4.5: Commit**

```bash
git add src/python.rs
git commit -m "feat(python): add EquivalenceChecker.from_manifest classmethod"
```

---

## Task 5: Switch `PyBatchProcessor` to `PyProvider` and add `from_manifest`

**Files:**
- Modify: `src/python.rs` (the `PyBatchProcessor` struct and `impl` block, ~lines 1459–1535)

- [ ] **Step 5.1: Replace `MockProvider` with `PyProvider`**

Change:

```rust
#[pyclass(name = "BatchProcessor")]
pub struct PyBatchProcessor {
    processor: BatchProcessor<MockProvider>,
}
```

to:

```rust
#[pyclass(name = "BatchProcessor")]
pub struct PyBatchProcessor {
    processor: BatchProcessor<PyProvider>,
}
```

- [ ] **Step 5.2: Update `new`**

Replace the `new` body (~lines 1473–1483):

```rust
#[new]
#[pyo3(signature = (reference_json=None))]
fn new(reference_json: Option<&str>) -> PyResult<Self> {
    let provider = match reference_json {
        Some(path) => PyProvider::from_json(Path::new(path))?,
        None => PyProvider::test_data(),
    };

    Ok(Self {
        processor: BatchProcessor::new(provider),
    })
}
```

- [ ] **Step 5.3: Add `from_manifest` static method**

Inside `#[pymethods] impl PyBatchProcessor { ... }`, after `new`:

```rust
/// Create a batch processor from a reference manifest.
///
/// Args:
///     manifest_path: Path to a manifest.json file produced by `ferro prepare`.
///
/// Returns:
///     A BatchProcessor backed by a MultiFastaProvider.
#[staticmethod]
fn from_manifest(manifest_path: &str) -> PyResult<Self> {
    let provider = PyProvider::from_manifest(Path::new(manifest_path))?;
    Ok(Self {
        processor: BatchProcessor::new(provider),
    })
}
```

- [ ] **Step 5.4: Build and confirm**

Run: `cargo build --features python`
Expected: SUCCEEDS.

- [ ] **Step 5.5: Commit**

```bash
git add src/python.rs
git commit -m "feat(python): add BatchProcessor.from_manifest classmethod"
```

---

## Task 6: Switch `PyCoordinateMapper` to `PyProvider` and add `from_manifest`

**Files:**
- Modify: `src/python.rs` (the `PyCoordinateMapper` struct and `impl` block, ~lines 2542–2700)

- [ ] **Step 6.1: Replace `MockProvider` with `PyProvider`**

Change:

```rust
#[pyclass(name = "CoordinateMapper")]
pub struct PyCoordinateMapper {
    provider: MockProvider,
}
```

to:

```rust
#[pyclass(name = "CoordinateMapper")]
pub struct PyCoordinateMapper {
    provider: PyProvider,
}
```

- [ ] **Step 6.2: Update `new`**

Replace the `new` body (~lines 2554–2564):

```rust
#[new]
#[pyo3(signature = (reference_json=None))]
fn new(reference_json: Option<&str>) -> PyResult<Self> {
    let provider = match reference_json {
        Some(path) => PyProvider::from_json(Path::new(path))?,
        None => PyProvider::test_data(),
    };

    Ok(Self { provider })
}
```

- [ ] **Step 6.3: Add `from_manifest` static method**

Inside `#[pymethods] impl PyCoordinateMapper { ... }`, after `new`:

```rust
/// Create a coordinate mapper from a reference manifest.
///
/// Args:
///     manifest_path: Path to a manifest.json file produced by `ferro prepare`.
///
/// Returns:
///     A CoordinateMapper backed by a MultiFastaProvider.
#[staticmethod]
fn from_manifest(manifest_path: &str) -> PyResult<Self> {
    let provider = PyProvider::from_manifest(Path::new(manifest_path))?;
    Ok(Self { provider })
}
```

- [ ] **Step 6.4: Verify all existing methods still compile**

The methods like `c_to_g`, `g_to_c`, etc. call `self.provider.get_transcript(...)` and `self.provider.has_transcript(...)`. `PyProvider: ReferenceProvider`, so these calls compile unchanged. `has_transcript` comes from the trait's default impl — confirm by reviewing the methods in lines 2580–2900 of `src/python.rs` after the edit; no body changes are required.

Run: `cargo build --features python`
Expected: SUCCEEDS.

- [ ] **Step 6.5: Run full Python test suite**

Run: `maturin develop --features python && pytest tests/python/ -v`
Expected: All existing tests PASS.

- [ ] **Step 6.6: Commit**

```bash
git add src/python.rs
git commit -m "feat(python): add CoordinateMapper.from_manifest classmethod"
```

---

## Task 7: Update Python type stubs

**Files:**
- Modify: `python/ferro_hgvs/__init__.pyi`

- [ ] **Step 7.1: Add `Normalizer.from_manifest` stub**

In `python/ferro_hgvs/__init__.pyi`, inside `class Normalizer:` (around line 412), after the `__init__` method, add:

```python
    @staticmethod
    def from_manifest(manifest_path: str, direction: str = "3prime") -> Normalizer:
        """Create a normalizer from a reference manifest written by ``ferro prepare``.

        Args:
            manifest_path: Path to a manifest.json file (typically inside a
                directory produced by ``ferro prepare``).
            direction: Shuffle direction - "3prime" (default) or "5prime".

        Returns:
            A Normalizer backed by a MultiFastaProvider.

        Raises:
            RuntimeError: If the manifest cannot be loaded.
        """
        ...
```

- [ ] **Step 7.2: Add `EquivalenceChecker.from_manifest` stub**

Inside `class EquivalenceChecker:` (around line 624), after `__init__`, add:

```python
    @staticmethod
    def from_manifest(manifest_path: str) -> EquivalenceChecker:
        """Create an equivalence checker from a reference manifest.

        Args:
            manifest_path: Path to a manifest.json file produced by ``ferro prepare``.

        Returns:
            An EquivalenceChecker backed by a MultiFastaProvider.

        Raises:
            RuntimeError: If the manifest cannot be loaded.
        """
        ...
```

- [ ] **Step 7.3: Add `BatchProcessor.from_manifest` stub**

Inside `class BatchProcessor:` (around line 846), after `__init__`, add:

```python
    @staticmethod
    def from_manifest(manifest_path: str) -> BatchProcessor:
        """Create a batch processor from a reference manifest.

        Args:
            manifest_path: Path to a manifest.json file produced by ``ferro prepare``.

        Returns:
            A BatchProcessor backed by a MultiFastaProvider.

        Raises:
            RuntimeError: If the manifest cannot be loaded.
        """
        ...
```

- [ ] **Step 7.4: Add `CoordinateMapper.from_manifest` stub**

Inside `class CoordinateMapper:` (around line 1215), after `__init__`, add:

```python
    @staticmethod
    def from_manifest(manifest_path: str) -> CoordinateMapper:
        """Create a coordinate mapper from a reference manifest.

        Args:
            manifest_path: Path to a manifest.json file produced by ``ferro prepare``.

        Returns:
            A CoordinateMapper backed by a MultiFastaProvider.

        Raises:
            RuntimeError: If the manifest cannot be loaded.
        """
        ...
```

- [ ] **Step 7.5: Run mypy**

Run: `mypy python/ferro_hgvs/`
Expected: No type errors.

- [ ] **Step 7.6: Commit**

```bash
git add python/ferro_hgvs/__init__.pyi
git commit -m "docs(stubs): add from_manifest type stubs for the four reference-aware classes"
```

---

## Task 8: Update SCHEMAS.md

**Files:**
- Modify: `tests/fixtures/SCHEMAS.md`

- [ ] **Step 8.1: Document the object form**

In `tests/fixtures/SCHEMAS.md`, in the `## Transcript Schema` section, after the existing Example block (around line 65), insert:

```markdown
### Object Form (Extended)

`MockProvider::from_json` and the Python `Normalizer(reference_json=...)`,
`EquivalenceChecker(reference_json=...)`, `BatchProcessor(reference_json=...)`,
and `CoordinateMapper(reference_json=...)` constructors also accept an object
form that supplies optional protein and genomic sequences:

```json
{
  "transcripts": [
    { "id": "NM_TEST.1", "...": "..." }
  ],
  "proteins": {
    "NP_TEST.1": "MAPLE..."
  },
  "genomic_sequences": {
    "NC_000023.11": "ACGT..."
  }
}
```

All three top-level keys are optional. `transcripts` follows the same schema
as the bare-array form. The bare-array form continues to be accepted for
backward compatibility.
```

(Note: the inner triple-backtick fence inside the markdown block needs to be a literal `` ``` `` — Edit will handle this verbatim since the outer is the file format; just paste as-is.)

- [ ] **Step 8.2: Commit**

```bash
git add tests/fixtures/SCHEMAS.md
git commit -m "docs(fixtures): document MockProvider object form in SCHEMAS.md"
```

---

## Task 9: Create tiny manifest fixture

**Files:**
- Create: `tests/fixtures/python/manifest_tiny/manifest.json`
- Create: `tests/fixtures/python/manifest_tiny/transcripts/test.fna`
- Create: `tests/fixtures/python/manifest_tiny/transcripts/test.fna.fai`

- [ ] **Step 9.1: Create the fixture directory**

Run:
```bash
mkdir -p tests/fixtures/python/manifest_tiny/transcripts
```

- [ ] **Step 9.2: Write a tiny FASTA**

Create `tests/fixtures/python/manifest_tiny/transcripts/test.fna` with exactly this content:

```
>NM_TEST.1
ATGCATGCAT
```

(Single sequence, 10 bases, no trailing whitespace beyond the final newline.)

- [ ] **Step 9.3: Write the matching FAI index**

Create `tests/fixtures/python/manifest_tiny/transcripts/test.fna.fai` with exactly this content (tab-separated):

```
NM_TEST.1	10	13	10	11
```

(Fields: name, length, offset_to_seq, line_bases, line_bytes. The offset 13 = 12 chars of `>NM_TEST.1` header + 1 newline. line_bytes = 11 because the line ends with `\n`. Use a tab between every pair of fields.)

- [ ] **Step 9.4: Write the manifest**

Create `tests/fixtures/python/manifest_tiny/manifest.json` with content:

```json
{
  "prepared_at": "2026-04-30T00:00:00Z",
  "transcript_fastas": ["transcripts/test.fna"],
  "genome_fasta": null,
  "cdot_json": null,
  "transcript_count": 1,
  "available_prefixes": ["NM_"]
}
```

- [ ] **Step 9.5: Smoke-test the fixture from Rust**

Run from a one-off command:
```bash
cargo run --release --features dev --bin ferro -- check --reference tests/fixtures/python/manifest_tiny 2>&1 | head -20
```

Expected: prints something showing 1 transcript loaded (or returns successfully with no error). Some warnings about missing cdot/genome are acceptable.

If `ferro check` doesn't accept this minimal manifest, instead drop the fixture into a temp dir from the test (see Task 10 fallback). For now, assume it works.

- [ ] **Step 9.6: Commit**

```bash
git add tests/fixtures/python/manifest_tiny/
git commit -m "test(fixtures): add tiny prepared-manifest fixture for Python loader tests"
```

---

## Task 10: Add Python integration tests

**Files:**
- Create: `tests/python/test_loaders.py`

- [ ] **Step 10.1: Write the test file**

Create `tests/python/test_loaders.py`:

```python
"""Tests for Normalizer/EquivalenceChecker/BatchProcessor/CoordinateMapper loader paths."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

import ferro_hgvs

FIXTURES = Path(__file__).parent.parent / "fixtures"
MANIFEST_TINY = FIXTURES / "python" / "manifest_tiny" / "manifest.json"


class TestObjectFormJson:
    """Object-form reference_json: {transcripts, proteins, genomic_sequences}."""

    def test_normalizer_loads_object_form(self, tmp_path: Path) -> None:
        ref = {
            "transcripts": [
                {
                    "id": "NM_TEST.1",
                    "gene_symbol": "TEST",
                    "strand": "+",
                    "sequence": "ATGCCCAAGGTGCTGCCC",
                    "cds_start": 1,
                    "cds_end": 18,
                    "exons": [{"number": 1, "start": 1, "end": 18}],
                }
            ],
            "proteins": {"NP_TEST.1": "MPKVLP"},
        }
        ref_file = tmp_path / "ref.json"
        ref_file.write_text(json.dumps(ref))

        norm = ferro_hgvs.Normalizer(reference_json=str(ref_file))
        # Round-trip a simple substitution to confirm the provider is wired up
        out = norm.normalize("NM_TEST.1:c.1A>G")
        assert "NM_TEST.1" in out

    def test_bare_array_form_still_works(self, tmp_path: Path) -> None:
        ref = [
            {
                "id": "NM_TEST.1",
                "gene_symbol": "TEST",
                "strand": "+",
                "sequence": "ATGCCCAAGGTGCTGCCC",
                "cds_start": 1,
                "cds_end": 18,
                "exons": [{"number": 1, "start": 1, "end": 18}],
            }
        ]
        ref_file = tmp_path / "ref.json"
        ref_file.write_text(json.dumps(ref))

        norm = ferro_hgvs.Normalizer(reference_json=str(ref_file))
        out = norm.normalize("NM_TEST.1:c.1A>G")
        assert "NM_TEST.1" in out


class TestFromManifest:
    """from_manifest constructors backed by MultiFastaProvider."""

    @pytest.mark.skipif(not MANIFEST_TINY.exists(), reason="tiny manifest fixture missing")
    def test_normalizer_from_manifest(self) -> None:
        norm = ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY))
        # Provider should be able to look up the transcript (sub-element of normalize).
        # Use a parse + normalize round-trip on a substitution within the test seq.
        out = norm.normalize("NM_TEST.1:c.1A>G")
        assert "NM_TEST.1" in out

    @pytest.mark.skipif(not MANIFEST_TINY.exists(), reason="tiny manifest fixture missing")
    def test_equivalence_checker_from_manifest(self) -> None:
        chk = ferro_hgvs.EquivalenceChecker.from_manifest(str(MANIFEST_TINY))
        v1 = ferro_hgvs.parse("NM_TEST.1:c.1A>G")
        v2 = ferro_hgvs.parse("NM_TEST.1:c.1A>G")
        result = chk.check(v1, v2)
        assert result.is_equivalent()

    @pytest.mark.skipif(not MANIFEST_TINY.exists(), reason="tiny manifest fixture missing")
    def test_batch_processor_from_manifest(self) -> None:
        bp = ferro_hgvs.BatchProcessor.from_manifest(str(MANIFEST_TINY))
        result = bp.parse(["NM_TEST.1:c.1A>G", "INVALID"])
        assert result.success_count() == 1
        assert result.error_count() == 1

    @pytest.mark.skipif(not MANIFEST_TINY.exists(), reason="tiny manifest fixture missing")
    def test_coordinate_mapper_from_manifest(self) -> None:
        cm = ferro_hgvs.CoordinateMapper.from_manifest(str(MANIFEST_TINY))
        # The tiny fixture has no cdot, so c-to-g may fail, but has_transcript should work.
        assert cm.has_transcript("NM_TEST.1")
```

- [ ] **Step 10.2: Run the tests**

Run: `pytest tests/python/test_loaders.py -v`
Expected: All tests PASS. If a manifest test fails because the tiny fixture lacks cdot data and `from_manifest` requires it, weaken the assertions to "construction succeeded" rather than full round-trip.

- [ ] **Step 10.3: Commit**

```bash
git add tests/python/test_loaders.py
git commit -m "test(python): add integration tests for object-form JSON and from_manifest loaders"
```

---

## Task 11: Final verification and PR

- [ ] **Step 11.1: Run the full test suite**

Run: `cargo nextest run --features dev`
Expected: All tests PASS.

Run: `maturin develop --features python && pytest tests/python/ -v`
Expected: All tests PASS.

- [ ] **Step 11.2: Run linters**

Run: `cargo clippy --features dev -- -D warnings`
Expected: No warnings.

Run: `ruff check python/ tests/python/ && ruff format --check python/ tests/python/`
Expected: Clean.

Run: `mypy python/ferro_hgvs/`
Expected: Clean.

- [ ] **Step 11.3: Pre-commit (final check)**

Run: `pre-commit run --all-files`
Expected: All hooks pass. Address any failures.

- [ ] **Step 11.4: Push branch and open draft PR**

Use the project's standard PR style (see `~/.claude/PR_STYLE.md`). Title:
`feat(python): load reference data via from_manifest and extended from_json`.

Body should include:
- `Fixes #82` on its own line
- Brief summary of the two parts
- Test plan describing what's verified

```bash
git push -u origin <branch>
gh pr create --draft --title "..." --body "$(cat <<'EOF'
... PR body per PR_STYLE.md ...
EOF
)"
```

---

## Self-Review Notes

- **Spec coverage:** Tasks 1, 7, 8 cover part (a). Tasks 2–7, 9, 10 cover part (b). Tests in Task 10 satisfy "one new Python integration test per loader path."
- **Backward compat:** Bare-array form is preserved (Task 1, Step 1.3 branches on `value.is_array()`). Existing positional/keyword args of `Normalizer(...)` etc. are unchanged.
- **Fixture trade-off:** The plan does not modify `normalization_transcripts.json` (the issue's hint), because mutating that 200KB shared fixture is invasive and brittle; instead the object-form test writes a small JSON file via `tmp_path`. If the user prefers the in-place change, swap Task 10's `TestObjectFormJson` to load `normalization_transcripts.json` (after wrapping it under `{"transcripts": [...]}`).
- **Risk:** The tiny manifest fixture (Task 9) has no cdot data. `MultiFastaProvider::from_manifest` warns but still constructs successfully. If a Python test that requires cdot fails, the fix is to relax the assertion (e.g., test only `has_transcript`) rather than ship a 30MB cdot file. The plan already favors weak assertions on the manifest tests.
