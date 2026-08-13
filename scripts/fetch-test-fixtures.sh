#!/usr/bin/env bash
#
# Fetch the large test fixtures that are hosted as GitHub release assets rather
# than committed to the repository.
#
# The repository used to track these four corpora (156 MB) with Git LFS. LFS
# bandwidth is metered even on a public repository and ten CI jobs fetched them
# per run, which exhausted the account's LFS budget and turned every `lfs: true`
# job red at checkout. Release-asset downloads are not metered, so they live on
# a release instead and are fetched on demand.
#
# Usage:
#
#   scripts/fetch-test-fixtures.sh            # fetch anything missing or stale
#   scripts/fetch-test-fixtures.sh --force    # re-download even if verified
#   scripts/fetch-test-fixtures.sh --verify   # verify what is on disk; fetch nothing
#
# Every download is verified against `tests/fixtures/CHECKSUMS.sha256`. A failed
# download, a missing asset, or a digest mismatch is a HARD failure: the suites
# that read these fixtures skip green when the files are absent, so a fetch that
# failed quietly would delete the ClinVar, CMRG and Paraphase coverage while
# making CI look faster.
#
# No authentication is used or required — the assets are public, so this works
# on a forked-PR runner with no token.

set -euo pipefail

readonly REPO="${FERRO_FIXTURE_REPO:-fulcrumgenomics/ferro-hgvs}"
readonly TAG="${FERRO_FIXTURE_TAG:-test-fixtures-v1}"
readonly BASE_URL="https://github.com/${REPO}/releases/download/${TAG}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
readonly repo_root
readonly CHECKSUMS="${repo_root}/tests/fixtures/CHECKSUMS.sha256"

mode="fetch"
case "${1:-}" in
    --force) mode="force" ;;
    --verify) mode="verify" ;;
    --help | -h)
        sed -n '3,26p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
        exit 0
        ;;
    "") ;;
    *)
        echo "error: unknown argument '$1' (expected --force, --verify or --help)" >&2
        exit 2
        ;;
esac
readonly mode

# `sha256sum` on Linux, `shasum -a 256` on macOS. Resolved once rather than per
# file so an environment with neither fails immediately and by name.
if command -v sha256sum >/dev/null 2>&1; then
    sha256_of() { sha256sum "$1" | cut -d' ' -f1; }
elif command -v shasum >/dev/null 2>&1; then
    sha256_of() { shasum -a 256 "$1" | cut -d' ' -f1; }
else
    echo "error: neither 'sha256sum' nor 'shasum' is available; cannot verify fixtures" >&2
    exit 1
fi

if [[ ! -f "${CHECKSUMS}" ]]; then
    echo "error: checksum manifest not found: ${CHECKSUMS}" >&2
    exit 1
fi

# Parse the manifest into parallel arrays. Comment and blank lines are skipped;
# everything else must be `<64-hex digest><spaces><repo-relative path>`, and a
# line that is not is an error rather than something to ignore quietly.
declare -a digests=() paths=()
while read -r digest path; do
    [[ -z "${digest}" || "${digest}" == \#* ]] && continue
    if [[ ! "${digest}" =~ ^[0-9a-f]{64}$ || -z "${path}" ]]; then
        echo "error: malformed line in ${CHECKSUMS}: '${digest} ${path}'" >&2
        exit 1
    fi
    digests+=("${digest}")
    paths+=("${path}")
done <"${CHECKSUMS}"

if [[ ${#paths[@]} -eq 0 ]]; then
    echo "error: ${CHECKSUMS} lists no fixtures — refusing to report success over an empty set" >&2
    exit 1
fi

echo "Fixtures: ${#paths[@]} file(s) from ${REPO}@${TAG}"

failures=0
for i in "${!paths[@]}"; do
    rel="${paths[$i]}"
    want="${digests[$i]}"
    dest="${repo_root}/${rel}"
    asset="$(basename "${rel}")"

    if [[ -f "${dest}" && "${mode}" != "force" ]]; then
        got="$(sha256_of "${dest}")"
        if [[ "${got}" == "${want}" ]]; then
            echo "  ok       ${rel}"
            continue
        fi
        if [[ "${mode}" == "verify" ]]; then
            echo "  MISMATCH ${rel}" >&2
            echo "           expected ${want}" >&2
            echo "           actual   ${got}" >&2
            failures=$((failures + 1))
            continue
        fi
        echo "  stale    ${rel} (digest mismatch; re-downloading)"
    elif [[ "${mode}" == "verify" ]]; then
        echo "  MISSING  ${rel}" >&2
        failures=$((failures + 1))
        continue
    fi

    mkdir -p "$(dirname "${dest}")"
    # Download to a temporary file so a failed or truncated transfer can never
    # leave a partial file where the suites would read it as the real fixture.
    tmp="${dest}.download.$$"
    # shellcheck disable=SC2064  # expand ${tmp} now, at trap-installation time
    trap "rm -f '${tmp}'" EXIT
    echo "  fetch    ${rel}"
    if ! curl --fail --location --show-error --silent \
        --retry 3 --retry-delay 2 --retry-connrefused \
        --output "${tmp}" "${BASE_URL}/${asset}"; then
        echo "error: failed to download ${asset} from ${BASE_URL}" >&2
        rm -f "${tmp}"
        failures=$((failures + 1))
        continue
    fi

    got="$(sha256_of "${tmp}")"
    if [[ "${got}" != "${want}" ]]; then
        echo "error: digest mismatch for ${asset}" >&2
        echo "       expected ${want}" >&2
        echo "       actual   ${got}" >&2
        echo "       The release asset does not match ${CHECKSUMS}. The asset may have" >&2
        echo "       been replaced, or the download was corrupted. NOT installing it." >&2
        rm -f "${tmp}"
        failures=$((failures + 1))
        continue
    fi

    mv "${tmp}" "${dest}"
    trap - EXIT
    echo "  ok       ${rel}"
done

if [[ ${failures} -gt 0 ]]; then
    echo >&2
    echo "error: ${failures} of ${#paths[@]} test fixture(s) could not be obtained and verified." >&2
    echo "       This is fatal on purpose: the bulk-fixture and exhaustive suites SKIP" >&2
    echo "       GREEN when their fixtures are absent, so continuing would silently drop" >&2
    echo "       the ClinVar, CMRG and Paraphase coverage instead of failing." >&2
    exit 1
fi

echo "All ${#paths[@]} fixture(s) present and verified."
