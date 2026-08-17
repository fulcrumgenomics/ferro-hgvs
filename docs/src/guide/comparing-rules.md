# Comparing normalization rules (`FERRO_PARTITION`)

> **Unstable.** `FERRO_PARTITION` is an evaluation switch, **not** a supported feature. It is not covered by semantic versioning, its values may change, and it is expected to be removed once the normalization rule is settled. Do not depend on it in production pipelines.

Normalization cuts each changed region of sequence into allele members. `FERRO_PARTITION` selects which rule does that cutting, so a candidate rule can be measured against the shipped one over a real corpus before anything changes for users.

| Value | Rule |
|-------|------|
| unset / empty / `canonical-coalesced` | **The shipped rule**, and what every normal invocation uses. `canonical`, plus the `delins.md:44-47` merge: a split whose payload realigns as one block is re-spelled as a single `delins`. Applied after the downstream passes rather than at partition time, so it cuts identically to `canonical` and differs only in what survives. |
| `live` | The rule shipped up to and including v0.14.0: a single-gap alignment search plus two narrow escapes. **No longer the default** — set it by name to reproduce pre-flip output. |
| `shadow` | Cut only at alignment steps common to *every* minimal alignment. |
| `canonical` | The member-count-minimal minimal alignment, without the merge above. |

With the variable unset — or set to the empty string — output is byte-identical to a build with no switch at all.

> **The default moved in v0.15.0**, from `live` to `canonical-coalesced`. If you are comparing against a stored corpus normalized by v0.14.0 or earlier, the like-for-like arm is now `FERRO_PARTITION=live`, not the unset one. The change is disclosed as a representation change in `CHANGELOG.md`; the ruling behind it is that a description is derived from the **resulting sequence** rather than preserved from the input's spelling, which `partition_block` — `live`'s cutter — cannot do.

## Running an A/B comparison

Run the same input twice and diff the normalized descriptions. In `tsv` format the columns are `line, input, normalized, changed, status, detail`, so column 3 is the normalized string:

```bash
# 1. a baseline arm, named explicitly. `live` is the pre-v0.15.0 rule, which is
#    the one to use when the question is "what moved for my stored corpus";
#    name the arm rather than relying on unset, which is now the NEW rule.
FERRO_PARTITION=live ferro normalize --input variants.txt --reference /path/to/reference \
  --format tsv --error-mode lenient -j 10 > shipped.tsv

# 2. the candidate rule — here, the shipped default
ferro normalize --input variants.txt --reference /path/to/reference \
  --format tsv --error-mode lenient -j 10 > candidate.tsv

# 3. what moved
diff <(cut -f3 shipped.tsv) <(cut -f3 candidate.tsv) | head

# 4. how many moved, counting only rows that succeeded on both sides
#    (columns 3/5 are `normalized`/`status`; a row that failed carries no
#     normalized string, so comparing column 3 alone would score two different
#     failures as identical)
paste <(cut -f3,5 shipped.tsv) <(cut -f3,5 candidate.tsv) \
  | awk -F'\t' '$2 == "ok" && $4 == "ok" && $1 != $3' | wc -l

# and, separately, rows whose status itself changed
paste <(cut -f5 shipped.tsv) <(cut -f5 candidate.tsv) | awk -F'\t' '$1 != $2' | wc -l
```

This works on a **stock release build** — the switch is not behind a build feature, so no special binary or wheel is needed.

## Two traps worth knowing

**A misspelled value is refused, loudly.** `FERRO_PARTITION=canonicl` makes `ferro` exit with an error before it reads any input, naming the value you gave and the arms this build has:

```
FERRO_PARTITION="canonicl" is not a partitioner this build has. This build's arms
are: live, shadow, canonical, canonical-coalesced. Refusing rather than falling
back to `live`, because a bake-off served the shipped rule under a candidate's
name reports that the candidate changes nothing.
```

Naming *this build's* arms is the point: a value that exists on some other branch, or that used to exist, is reported as absent here rather than quietly answered as `live`.

The refusal comes from the CLI, not from the normalizer. Up to and including v0.14.0 it was a panic raised deep inside normalization, which meant a development-only switch could abort any process — a long-running service, or a Python caller, across the FFI boundary — that merely happened to have the variable set. A release build of the *library* now falls safe to `live` for a value that names no arm and keeps the refusal for its caller to report. Every binary and every example in this repository reports it — `ferro`, `ferro-web`, `ferro-benchmark`, both spec generators, and all 23 declared examples, which are the bake-off harnesses, the window extractors and the artifact generators — and that is pinned by a test whose denominator is `Cargo.toml`'s own `[[bin]]` and `[[example]]` tables rather than by this sentence, so a target added later cannot quietly skip it. The one class of cargo target still outside that obligation is `[[bench]]`: criterion's `criterion_main!` generates the `main`, so there is no hand-written entry point to put the call in. If you embed ferro as a library and run bake-offs through it, read `ferro_hgvs::normalize::partition_switch_startup_error()` at startup and do the same.

Builds up to and including v0.13.1 fell back to `live` instead, and produced a clean, empty diff that read as "the candidate changes nothing". The fallback emitted a warning through the `log` facade, but the `ferro` CLI installs no logger, so that warning reached no stream and `RUST_LOG` could not surface it — there was no signal at all. If you are on an older build, treat every empty diff as unproven.

A positive control on an input **known** to differ is still worth running, because it catches the other way a comparison can be vacuous: a variable that never reached the process at all (a lost `export`, a `sudo` that scrubbed the environment, a container that did not forward it). That case is indistinguishable from `unset`, which is legitimately `live`, so no amount of validation inside ferro can catch it. On your own corpus a zero remains ambiguous between "the switch is not taking effect" and "this corpus has no affected variants".

These three inputs run against the built-in test data, so they need no `--reference` and no prepared reference directory. Between them they separate all four arms — every pair of arms disagrees on at least one row, so the control tells you *which* arm you got, not merely that something changed:

```bash
printf 'NM_001234.1:c.[5_6insAC;9del]\nNM_001234.1:c.2_6delinsGA\nNM_001234.1:c.4_10delinsAC\n' > control.txt

for arm in live shadow canonical canonical-coalesced; do
  echo "== $arm"
  if FERRO_PARTITION=$arm ferro normalize --input control.txt --format tsv \
       --error-mode lenient > "control.$arm.tsv"; then
    tail -n +2 "control.$arm.tsv" | cut -f3
  else
    echo "   FAILED (exit $?) -- ferro did not produce this arm's column"
  fi
done
```

The status check is not boilerplate. An unrecognised arm now aborts the process, and a pipeline reports the exit status of its *last* command — so `ferro … | tail | cut` reports `cut`'s success and the loop prints an empty column under the arm's heading. An empty column and an aborted run look identical, which is the same "a broken measurement reads as a result" failure this whole section is about. Redirecting first and testing the status makes the abort say so.

| input | `live` | `shadow` | `canonical` | `canonical-coalesced` (= unset) |
|---|---|---|---|---|
| `c.[5_6insAC;9del]` | `c.[5_6insA;7_9delinsCAA]` | `c.[5_6insAC;11del]` | `c.[5_6insAC;11del]` | `c.[5_6insAC;11del]` |
| `c.2_6delinsGA` | `c.2_6delinsGA` | `c.2_6delinsGA` | `c.[2del;4_6delinsA]` | `c.2_6delinsGA` |
| `c.4_10delinsAC` | `c.4_10delinsAC` | `c.4_10delinsAC` | `c.[4C>A;6_10del]` | `c.[4C>A;6_10del]` |

**All six pairs of arms are separated**, which is what makes this tell you *which* arm you got rather than merely that something changed: row 1 separates `live` from the other three, row 2 isolates `canonical` (it is the only arm that does not merge that block), and row 3 separates `{live, shadow}` from `{canonical, canonical-coalesced}`. Rows 1 and 3 together separate `live` from `canonical-coalesced`; row 3 alone separates `shadow` from it.

Because `canonical-coalesced` is now the default, running with the variable **unset** must reproduce that last column exactly. If it reproduces the `live` column instead, you are on a pre-v0.15.0 build.

If the arm you selected does not produce its column above, the variable is not reaching ferro and any comparison you run is meaningless. Only once the control behaves is a zero on your own corpus informative.

> **This control has now been wrong twice, so check it rather than trusting it.** The first version offered `c.[2del;9del]`, `c.[3del;9del]` and `c.[2del;9dup]` and claimed `canonical` answers `c.[2del;33del]`; on those three inputs no arm differed from `live` at all. Its replacement — `c.[5_6insAC;9del]`, `c.[2del;5del]`, `c.[2del;9del]` — was **also** wrong in four of its twelve cells when re-measured: `live` row 1 read `c.[5_6insA;7_9delinsCAA]` and not `c.6_9delinsACCAA`, and rows 2 and 3 separated *nothing*, every arm answering `c.[2del;6del]` and `c.[2del;11del]`. Both failures share one cause: `NM_001234.1` is a G homopolymer from `c.9` to `c.33`, so deletion pairs inside it shuffle to a common form on every arm instead of partitioning differently. A discriminating row has to be a `delins` whose payload re-aligns, which is what the three above are. The table is measured, not composed.

**From Python, it must be set before the first normalization.** The value is read once per process and cached, so:

```python
import os
os.environ["FERRO_PARTITION"] = "canonical"   # must precede the first normalize call
import ferro_hgvs
```

Setting it after any variant has been normalized silently does nothing, and it cannot be changed within a running process. Comparing two rules therefore means two separate processes (or two CLI runs, as above), not two calls in one script.
