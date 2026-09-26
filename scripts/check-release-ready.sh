#!/usr/bin/env bash
# Verify Reader checks succeeded for the exact private release commit.

set -euo pipefail

ref="${1:-HEAD}"
# `^{commit}` peels annotated tags to the commit they point at; git rev-parse
# on a bare tag name otherwise returns the tag object's own SHA, which never
# has a CI/audit run against it.
sha="$(git rev-parse "${ref}^{commit}")"

echo "Checking release readiness for $ref ($sha)..."

check_workflow() {
    local workflow="$1"
    local run_json
    run_json="$(gh run list --repo pgarrett-scripps/opentfraw -w "$workflow" -c "$sha" --json status,conclusion,url -L 1)"

    if [[ "$run_json" == "[]" ]]; then
        echo "FAIL: no run of $workflow found for commit $sha."
        return 1
    fi

    local status conclusion url
    status="$(jq -r '.[0].status' <<<"$run_json")"
    conclusion="$(jq -r '.[0].conclusion' <<<"$run_json")"
    url="$(jq -r '.[0].url' <<<"$run_json")"

    if [[ "$status" != "completed" ]]; then
        echo "FAIL: latest $workflow run for $sha is not completed (status: $status). $url"
        return 1
    fi

    if [[ "$conclusion" != "success" ]]; then
        echo "FAIL: latest $workflow run for $sha did not succeed (conclusion: $conclusion). $url"
        return 1
    fi

    echo "OK: $workflow passed for $sha. $url"
    return 0
}

ok=1

check_workflow reader.yml || ok=0

if [[ "$ok" -ne 1 ]]; then
    exit 1
fi

echo "Release ready: Reader checks is green for $sha."
