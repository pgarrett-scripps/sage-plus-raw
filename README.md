# sage-plus-raw

Fork of Sigilweaver OpenTFRaw by Nathan Riley and contributors; not affiliated with Sigilweaver or Thermo Fisher Scientific.

Native Rust Thermo RAW reading for [Sage Plus](https://github.com/pgarrett-scripps/sage-plus), without .NET or vendor DLLs. The package is `sage-plus-raw` and the library is `sage_plus_raw`. It is not published to crates.io. Depend on it by Git revision:

```toml
[dependencies]
sage-plus-raw = { git = "https://github.com/pgarrett-scripps/opentfraw", rev = "FULL_REVIEWED_COMMIT" }
```

Upstream copyright and Apache-2.0 licensing are retained. See [LICENSE](LICENSE), [NOTICE](NOTICE), [ATTRIBUTION.md](ATTRIBUTION.md) and [CONTRIBUTORS.md](CONTRIBUTORS.md).

## Reading

Open with `RawFileReader::open_path`. Use `try_iter_spectra` for a stream of `Result<SpectrumRecord>` values, or `try_extract_spectrum` for one zero-based scan index. `try_extract_spectrum_arrays` returns full-precision `f64` arrays. `scan_metadata` reads metadata without peaks.

```sh
cargo run --locked -p sage-plus-raw --example strict_read -- acquisition.raw
```

- Scan events are decoded from their self-describing layout (counted reactions, mass ranges, calibrators, source fragmentation lists and, in v66, a label string). The layout is accepted only when the stream ends exactly. Otherwise the reader falls back to the fixed-size layouts that opentfraw 1.4.1 used, so every file 1.4.1 opens still opens.
- Every scan above MS1 has a precursor. For multiplexed (MSX) and multi-stage (EThcD, SPS-MS3) events it describes the first reaction; every reaction stays in `raw.scan_events`.
- `master_scan_number` comes from the "Master Scan Number" trailer, or else the nearest earlier scan one MS level lower. "Master Index" is a 0/1 flag and is never used.
- The isolation target comes from "Target M/Z" or the first reaction. "MS2 Isolation Offset" is an offset, exposed by `ScanParams::isolation_offset_mz`.
- A peak with a non-finite or non-positive m/z, or an intensity that is not finite or overflows `f32`, is dropped and counted in `dropped_peaks`. It does not fail the scan.
- Unreadable scans return errors with operation, scan and byte-offset context (`Error::kind`, `Error::context`, `Error::root_cause`).
- Activation is reported only when qualified (v66 reaction tags 11 = HCD, 9 = ETD). Nondependent MS2 is not evidence of DIA; `is_dia` is a legacy heuristic.

## Validation

```sh
cargo fmt --all --check
cargo clippy --locked -p sage-plus-raw --all-targets -- -D warnings
cargo test --locked -p sage-plus-raw
python3 scripts/fetch_public_fixtures.py
python3 scripts/prepare_reference_arrays.py
THERMO_PUBLIC_CORPUS=corpus/public cargo test --locked -p sage-plus-raw --test public_fixtures -- --ignored
```

The public fixtures are ProteoWizard RAW files and reference mzML, hash-checked on download. The larger evaluation corpus is described in [CORPUS.md](CORPUS.md) and fetched with `scripts/fetch_corpus.py`. Private acquisitions are never committed.

Fuzz targets live in `crates/sage-plus-raw/fuzz`. Seed them with `python scripts/seed_fuzz.py`, then run for example `cargo +nightly fuzz run strict_decode -- -max_total_time=20` from `crates/sage-plus-raw`.

`cargo run --release --locked -p sage-plus-raw --example reader_benchmark -- acquisition.raw 3 42` prints JSON timings for open, sequential and shuffled reads.
