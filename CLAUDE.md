# CLAUDE.md

A pure-Rust reader for Thermo `.raw` files: no .NET, no vendor DLLs, `unsafe_code = "forbid"`,
dependencies only `thiserror` and `roxmltree`, one workspace member `crates/sage-plus-raw` (library
`sage_plus_raw`). A fork of [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw); upstream's Python
crate, mzML writer and docs site are deleted here and stay deleted. It feeds
[sage-plus](https://github.com/pgarrett-scripps/sage-plus), which consumes it **from crates.io**,
not from this folder: nothing changed here reaches sage-plus until a release is published.

## The public surface sage-plus pins

`RawFileReader::open_path`, `try_iter_spectra`, `SpectrumRecord`, `PrecursorInfo`, `Error`. Changing
any of those signatures, or removing a public item or output column, breaks sage-plus's build;
additions are fine. `Error` is not `#[non_exhaustive]`, so a new variant is breaking (`ErrorKind`
is). Ask before any removal or signature change. sage-plus reads activation **only** from the
rendered filter string (`sage/src/spectrum.rs`, `AcquisitionGroup::from_thermo_filter`), so a change
in `scan_filter.rs` changes sage-plus's behaviour: keep existing filter strings byte-identical
unless that is the task.

## Module layout and invariants

- Byte I/O: `binary_reader`, `bytes`. Untrusted lengths are capped
  (`Error::AllocationTooLarge`); keep it that way.
- File metadata: `header`, `raw_file_info`, `run_header`, `sample_info`, `seq_row`, `audit_tag`,
  `error_log`, `device`, `instrument_method`, `compound_file` (MS-CFB method container).
- Trailer schema: `generic_data`, `metadata`. Per-scan decode: `scan_index`, `scan_event`,
  `event_stream`, `scan_data`, `scan_format`. Extraction: `spectrum`. Rendering: `scan_filter`.
- `reader.rs` (~1.5k lines) is the only hard file: the schema-search and event-layout inference
  heuristics live there. Read it before touching event decoding.
- `event_stream::decode` reads v66 events from their self-describing layout and the fixed-size
  tables are only a per-file fallback; every scan above MS1 has a precursor; raw reaction words
  are preserved, not normalised; a malformed file returns an `Error`, never a panic
  (`clippy::unwrap_used` is warned outside tests).

## Evidence rule

A decoded byte needs public corpus evidence, cited at the point of use (`CONTRIBUTING.md`,
`CORPUS.md`). No vendor code, decompiled implementations or SDK headers; prefer synthetic fixtures
for malformed input. Never commit private acquisitions, credentials or paths that identify research
participants.

## Running things

MSRV is **1.85** and stays there; upstream's 1.88 APIs (`slice::as_chunks`, `usize::is_multiple_of`)
are off limits. Several agents share this machine, so every workspace cargo build, test or clippy
goes through the scheduler:

```sh
ajs submit --title "sage-plus-raw local gate" --description "what this gates" --meta est=15m \
  --cpu 6 --mem 12G --max-runtime 45m \
  -- bash -lc '/home/patrick-garrett/Repos/sage-plus-raw/scripts/local-gate.sh'
```

`scripts/local-gate.sh` is the check to run before a commit, not a hand-written subset: it runs
every one of the ten CI jobs that a single Linux machine can run, in the form CI runs it, labels
each step with the job it stands for, and ends by printing what it covered, what it covered partly
and what it could not. It keeps going after a failure, so one run names everything that is wrong;
`LOCAL_GATE_PASSED` on the last line is the thing to look for. A subset is how a change can reach
CI with a fuzz harness that no longer compiles: nine jobs pass and `fuzz smoke` fails, because the
ad-hoc command covered the crate and not `crates/sage-plus-raw/fuzz`, which is its own workspace
with its own `Cargo.lock` that CI resolves with `--locked`, so a version bump must touch it too.
Only `cargo fmt` and a single-crate `cargo check -p sage-plus-raw` run directly.

Releases follow `RELEASING.md`: the tag push publishes to crates.io from CI, and only Ty approves
one.

## Fixtures

Plain `cargo test` runs the inline unit tests only. The integration tests are `#[ignore]` and gated
on `THERMO_PUBLIC_CORPUS`, `THERMO_LARGE_RAW`, `OPENTFRAW_LUMOS_PROFILE_RAW`, the three
`OPENTFRAW_LOCK_MASS_*` (`HFX`, `PLUS`, `EXPLORIS`) and the two `OPENTFRAW_SOURCE_CID_*` (`RAW`,
`PLUS_RAW`) variables. Fetch the public ones with `python3 scripts/fetch_public_fixtures.py`, pinned
by SHA-256 in `scripts/public-fixtures.json`; they are test inputs, not ours to rewrite.

## Upstream fixes

`upstream` is Sigilweaver/OpenTFRaw. Do not merge the branch: it conflicts in ~8 library sources and
in every file this fork deleted. Cherry-pick one commit at a time with the rename applied, keeping
upstream authorship, adding `(cherry picked from commit <sha>)` and saying why it needed adapting:

```sh
git format-patch -1 --stdout <sha> | sed 's|crates/opentfraw/|crates/sage-plus-raw/|g' | git am -3
```

Deferred by decision, do not pick: `eb79dff`, `b2418f2`, `e5f4aff` (API-breaking, waiting on a
0.2.0), the openmassspec-core 3.0 series, and `d86f574` (MSRV 1.88). `audits/sage-plus-raw.md` in
the agent-teams sageplus project classifies every upstream commit.
