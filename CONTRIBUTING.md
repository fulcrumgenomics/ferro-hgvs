# Contributing to ferro-hgvs

Thank you for considering contributing to ferro-hgvs!

## Getting Started

### Prerequisites

- Rust (stable)
- Git
- [uv](https://docs.astral.sh/uv/) (for Python development)

### Setup

```bash
git clone https://github.com/fulcrumgenomics/ferro-hgvs.git
cd ferro-hgvs
cargo build
cargo test --features dev
```

### Python Bindings Setup

```bash
uv sync --group dev
uv run maturin develop --features python
uv run pytest
```

After modifying Python dependencies in `pyproject.toml`, run `uv lock` and commit
the updated `uv.lock`. CI uses `--locked` and will fail if the lockfile is out of
sync with `pyproject.toml`.

## Development Workflow

### Making Changes

1. Create a branch: `git checkout -b feature/your-feature-name`
2. Write your code and add tests
3. Ensure all tests pass: `cargo nextest run --features dev` (or `cargo test --features dev`)
4. Run lints: `cargo clippy --features dev -- -D warnings && cargo fmt --check`

### Commit Messages

Follow conventional commit format:

```
type(scope): description
```

Types: `feat`, `fix`, `docs`, `refactor`, `test`, `chore`

Examples:
```
feat(parser): add support for repeat variants
fix(normalize): correct boundary detection for UTR regions
```

### Submitting a Pull Request

1. Push your branch and open a PR on GitHub
2. Fill out the PR template
3. Wait for CI to pass
4. Request review from maintainers

## Testing

```bash
cargo nextest run --features dev          # All tests (preferred)
cargo test --features dev                 # Alternative
cargo nextest run -E 'test(test_name)'    # Specific test
cargo test -- --nocapture                 # With output
cargo bench                              # Benchmarks
```

## Code Style

- Format with `cargo fmt`
- Lint with `cargo clippy --features dev -- -D warnings`
- Document all public APIs with doc comments

## Updating the v21.0 Normalization Fixture

`tests/fixtures/grammar/hgvs_spec_normalization.json` pins ferro's current
`normalize()` output for every variant string in the HGVS v21.0 spec
(vendored at `assets/hgvs-nomenclature/`). The companion test
`tests/hgvs_spec_normalization_tests.rs` fails any time a row's observed
output drifts from the recorded `current`.

If your PR changes any normalization output:

1. Make sure the vendored spec submodule is checked out, then regenerate
   the fixture:

   ```bash
   git submodule update --init assets/hgvs-nomenclature
   cargo run --features dev --example generate_spec_fixture
   ```

2. Inspect the diff. For each changed row, verify the new `current` against
   the v21.0 spec text under `assets/hgvs-nomenclature/docs/recommendations/`.
   When ferro's new output now matches the spec's canonical form, that row
   should also have `current == spec_expected` after regen.

3. If a row's spec-canonical form differs from the input string (a "pair" — e.g.
   spec input `c.79GC>TT` is canonicalized to `c.79_80delinsTT`), record the
   divergence in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`
   keyed on the input. Override entry shape:

   ```jsonc
   {
     "by_input": {
       "<exact input string>": {
         "status": "pair",                          // optional: compliant | pair | reject | needs-reference
         "spec_expected": "<spec's canonical form>", // optional: defaults to input
         "todo": "<https://… link>"                 // optional: defaults to a #83 link when current != spec_expected
       }
     }
   }
   ```

   Override keys must match a real fixture input — typos are caught by the
   generator.

4. CI verifies byte-identical regeneration via:

   ```bash
   cargo run --features dev --example generate_spec_fixture -- --check
   ```

5. To bump the spec to a newer upstream version, update the submodule pointer
   under `assets/hgvs-nomenclature/` (e.g., `git -C assets/hgvs-nomenclature
   checkout <new-tag>`) and regenerate the fixture.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
