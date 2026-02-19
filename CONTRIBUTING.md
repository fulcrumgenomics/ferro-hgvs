# Contributing to ferro-hgvs

Thank you for considering contributing to ferro-hgvs!

## Getting Started

### Prerequisites

- Rust (stable)
- Git

### Setup

```bash
git clone https://github.com/fulcrumgenomics/ferro-hgvs.git
cd ferro-hgvs
cargo build
cargo test --features dev
```

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

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
