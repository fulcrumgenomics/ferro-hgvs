# Error Handling

ferro-hgvs provides configurable error handling with three modes:

| Mode | Behavior |
|------|----------|
| `strict` | Reject non-conformant input (default) |
| `lenient` | Auto-correct with warnings |
| `silent` | Auto-correct silently |

```bash
# Use lenient mode to auto-correct common issues
ferro parse --error-mode lenient "p.val600glu"  # Corrects to p.Val600Glu

# Ignore specific warnings
ferro parse --ignore W1001,W2001 "p.val600glu"

# Get help on any error/warning code
ferro explain W1001
ferro explain --list
```

## Configuration File

Create `.ferro.toml` in your project directory:

```toml
[error-handling]
mode = "lenient"
ignore = ["W1001", "W2001"]  # Silently correct these
reject = ["W3003"]           # Always reject these
```
