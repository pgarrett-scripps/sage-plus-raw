#!/usr/bin/env bash
# Run every CI check that can run on one Linux machine, in one command.
#
#   scripts/local-gate.sh
#
# Each step is labelled with the CI job it stands for, and the summary at the end
# says which of the ten jobs this covered and which it could not. The exit status
# is 0 only when every step that ran passed; every step runs even after one
# fails, so a single run names everything that is wrong.
#
# Two jobs cannot run here at all: `test (macos-latest)` and `test
# (windows-latest)` need those operating systems. A third, `fuzz smoke`, runs its
# three libFuzzer targets under a pinned nightly with cargo-fuzz and
# AddressSanitizer; this script compiles the fuzz workspace against its own
# lockfile instead, which is what catches a source change that breaks the
# harness. See CONTRIBUTING.md.
#
# Builds go to ~/.cache by default, on the root drive, not into the repository.
set -uo pipefail

cd "$(dirname "$0")/.."
repo="$PWD"

export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$HOME/.cache/sage-plus-raw-target}"
export CARGO_TERM_COLOR="${CARGO_TERM_COLOR:-always}"
export CARGO_INCREMENTAL=0
mkdir -p "$CARGO_TARGET_DIR"

# Keep these in sync with .github/workflows/ci.yml.
rust_msrv=1.85.0
fuzz_nightly=nightly-2026-09-26
fuzz_manifest=crates/sage-plus-raw/fuzz/Cargo.toml
python="${PYTHON:-python3}"

labels=()
jobs=()
results=()

run() {
    local label="$1" job="$2"
    shift 2
    printf '\n=== %s (CI: %s) ===\n' "$label" "$job"
    printf '    %s\n' "$*"
    if "$@"; then
        labels+=("$label")
        jobs+=("$job")
        results+=(pass)
    else
        printf '!!! %s FAILED\n' "$label"
        labels+=("$label")
        jobs+=("$job")
        results+=(FAIL)
    fi
}

skip() {
    local label="$1" job="$2" why="$3"
    printf '\n=== %s (CI: %s) === skipped: %s\n' "$label" "$job" "$why"
    labels+=("$label")
    jobs+=("$job")
    results+=("skip: $why")
}

# rustfmt: the crate workspace, as CI does.
run "fmt, crate" rustfmt cargo fmt --all --check

# The fuzz workspace is formatted too, which CI does not check: its targets are
# edited by hand and `cargo fmt --all` from the root does not reach them.
run "fmt, fuzz workspace (beyond CI)" rustfmt \
    cargo fmt --all --check --manifest-path "$fuzz_manifest"

# clippy: every target of the crate, denying warnings.
run "clippy" clippy \
    cargo clippy --locked -p sage-plus-raw --all-targets -- -D warnings

# test: the inline unit and doc tests first, then the gated public-fixture
# suite, exactly as the ubuntu job runs it.
run "test, unit and doc" "test (ubuntu-latest)" \
    cargo test --locked -p sage-plus-raw

# The fixtures are pinned by SHA-256 and verified, so this is a no-op once they
# are present; it needs the network the first time.
run "public fixtures fetched" "test (ubuntu-latest)" \
    "$python" scripts/fetch_public_fixtures.py
run "reference arrays prepared" "test (ubuntu-latest)" \
    "$python" scripts/prepare_reference_arrays.py
run "test, public fixtures" "test (ubuntu-latest)" \
    env THERMO_PUBLIC_CORPUS=corpus/public \
    cargo test --locked -p sage-plus-raw --test public_fixtures -- --ignored

# msrv: the published crate must build on the oldest supported compiler.
if rustup toolchain list 2>/dev/null | grep -q "^$rust_msrv"; then
    run "msrv build" "MSRV ($rust_msrv)" \
        cargo "+$rust_msrv" build --locked -p sage-plus-raw --all-targets
else
    skip "msrv build" "MSRV ($rust_msrv)" \
        "install it with: rustup toolchain install $rust_msrv"
fi

# docs: a broken intra-doc link fails CI.
run "docs" docs env RUSTDOCFLAGS="-D warnings" \
    cargo doc --locked -p sage-plus-raw --no-deps

# deny: advisories, licences, bans and sources, the four checks CI names.
if command -v cargo-deny >/dev/null 2>&1; then
    # cargo-deny takes --locked before the subcommand, unlike cargo itself, so
    # `cargo deny check --locked` is not a valid command line.
    run "deny" "cargo deny" cargo deny check advisories licenses bans sources
else
    skip "deny" "cargo deny" "install it with: cargo install cargo-deny --locked"
fi

# package: the licence files the published crate carries must match the
# repository root, and the dry run must build the packaged sources.
check_licences() {
    cmp LICENSE crates/sage-plus-raw/LICENSE || return 1
    cmp NOTICE crates/sage-plus-raw/NOTICE || return 1
    echo "LICENSE and NOTICE match crates/sage-plus-raw/"
}

package_files="$CARGO_TARGET_DIR/package-files.txt"
check_package() {
    cargo package --locked -p sage-plus-raw --list >"$package_files" || return 1
    # `include` in crates/sage-plus-raw/Cargo.toml decides this list: the two
    # licence files, the README it copies in from the root, the sources, the
    # examples and the integration tests. CHANGELOG.md is not packaged.
    local required
    for required in LICENSE NOTICE README.md Cargo.toml src/lib.rs; do
        grep -qx "$required" "$package_files" || {
            echo "missing $required from package"
            return 1
        }
    done
    if grep -E '^(fuzz|benchmarks|scripts|corpus|\.github)/' "$package_files"; then
        echo "development-only files in package"
        return 1
    fi
    wc -l <"$package_files" | xargs printf '%s files in the package\n'
}

# `cargo publish --dry-run` refuses an unclean tree, which a local run often has.
dirty_flag=()
dirty_note=""
if [ -n "$(git -C "$repo" status --porcelain)" ]; then
    dirty_flag=(--allow-dirty)
    dirty_note=" (uncommitted changes: packaged with --allow-dirty)"
fi

run "packaged licences" "publish dry-run" check_licences
run "package contents" "publish dry-run" check_package
run "publish dry run" "publish dry-run" \
    cargo publish --locked -p sage-plus-raw --dry-run "${dirty_flag[@]}"

# fuzz smoke, as far as it goes here: resolve the fuzz workspace against its own
# lockfile and compile every target. A source change that breaks the harness
# fails here instead of in the one CI job the other nine do not stand in for.
run "fuzz lockfile resolves" "fuzz smoke" \
    cargo metadata --locked --manifest-path "$fuzz_manifest" --format-version 1
run "fuzz workspace compiles" "fuzz smoke" \
    env CARGO_TARGET_DIR="$CARGO_TARGET_DIR/fuzz" \
    cargo check --locked --all-targets --manifest-path "$fuzz_manifest"
run "fuzz lock unchanged" "fuzz smoke" \
    git diff --exit-code -- crates/sage-plus-raw/fuzz/Cargo.lock
run "fuzz corpus seeds" "fuzz smoke" "$python" scripts/seed_fuzz.py

if command -v cargo-fuzz >/dev/null 2>&1; then
    skip "fuzz target runs" "fuzz smoke" \
        "cargo-fuzz is installed; CI runs the three targets under $fuzz_nightly with AddressSanitizer"
else
    skip "fuzz target runs" "fuzz smoke" \
        "needs cargo-fuzz 0.13.1 and $fuzz_nightly; CI runs the three targets"
fi

skip "test, macOS" "test (macos-latest)" "no macOS host"
skip "test, Windows" "test (windows-latest)" "no Windows host"

printf '\n===== summary =====\n'
failed=0
for i in "${!labels[@]}"; do
    case "${results[$i]}" in
    pass) printf '  pass  %-42s %s\n' "${labels[$i]}" "[${jobs[$i]}]" ;;
    FAIL)
        printf '  FAIL  %-42s %s\n' "${labels[$i]}" "[${jobs[$i]}]"
        failed=1
        ;;
    *) printf '  skip  %-42s %s  %s\n' "${labels[$i]}" "[${jobs[$i]}]" "${results[$i]#skip: }" ;;
    esac
done

printf '\nCI jobs fully covered here: rustfmt, clippy, test (ubuntu-latest), MSRV (%s), docs, cargo deny, publish dry-run\n' "$rust_msrv"
printf 'CI jobs partly covered:     fuzz smoke (resolves, compiles and seeds; does not fuzz)\n'
printf 'CI jobs not covered:        test (macos-latest), test (windows-latest)\n'
printf 'Gated tests needing private data stay #[ignore]d here and in CI: THERMO_LARGE_RAW,\n'
printf '  OPENTFRAW_LUMOS_PROFILE_RAW, OPENTFRAW_LOCK_MASS_*, OPENTFRAW_SOURCE_CID_* (CORPUS.md).\n'
printf 'Target directory: %s%s\n' "$CARGO_TARGET_DIR" "$dirty_note"

if [ "$failed" -ne 0 ]; then
    printf '\nLOCAL_GATE_FAILED\n'
    exit 1
fi
printf '\nLOCAL_GATE_PASSED\n'
