# sage-plus-raw

[![CI](https://github.com/pgarrett-scripps/sage-plus-raw/actions/workflows/ci.yml/badge.svg)](https://github.com/pgarrett-scripps/sage-plus-raw/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/sage-plus-raw.svg)](https://crates.io/crates/sage-plus-raw)
[![docs.rs](https://img.shields.io/docsrs/sage-plus-raw)](https://docs.rs/sage-plus-raw)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

Fork of Sigilweaver OpenTFRaw by Nathan Riley and contributors; not affiliated with Sigilweaver or Thermo Fisher Scientific.

Native Rust Thermo RAW reading for [Sage Plus](https://github.com/pgarrett-scripps/sage-plus), without .NET or vendor DLLs. The package is `sage-plus-raw` and the library is `sage_plus_raw`. The minimum supported Rust version is 1.85.

Upstream copyright and Apache-2.0 licensing are retained. See [LICENSE](LICENSE), [NOTICE](NOTICE), [ATTRIBUTION.md](ATTRIBUTION.md) and [CONTRIBUTORS.md](CONTRIBUTORS.md).

## Install

```sh
cargo add sage-plus-raw
```

or in `Cargo.toml`:

```toml
[dependencies]
sage-plus-raw = "0.2"
```

## Usage

```rust,no_run
use sage_plus_raw::{try_iter_spectra, RawFileReader};
use std::{fs::File, io::BufReader};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let raw = RawFileReader::open_path("acquisition.raw")?;
    let mut source = BufReader::new(File::open("acquisition.raw")?);
    for spectrum in try_iter_spectra(&raw, &mut source, false) {
        let spectrum = spectrum?;
        println!("scan {}: {} peaks", spectrum.scan_number, spectrum.mz.len());
    }
    Ok(())
}
```

## Reading

Open with `RawFileReader::open_path`. Use `try_iter_spectra` for a stream of `Result<SpectrumRecord>` values, or `try_extract_spectrum` for one zero-based scan index. `try_extract_spectrum_arrays` returns full-precision `f64` arrays. `scan_metadata` reads metadata without peaks.

A reader keeps only what decoding reads of a file's per-scan trailer: the ten values in `RawFileReader::scan_trailers`, 136 bytes a scan. To read any other trailer or instrument-log field through `scan_params` or `status_log_entry`, open the file with `RawFileReader::open_path_with(path, OpenOptions::with_records())`, which keeps the records themselves (on a 136k-scan Orbitrap file, 550 MiB and 10 million allocations more). Without them those accessors say so rather than report a missing field.

```sh
cargo run --locked -p sage-plus-raw --example strict_read -- acquisition.raw
```

- Scan events are decoded from their self-describing layout (counted reactions, mass ranges, calibrators, source fragmentation lists and, in v66, a label string). The layout is accepted only when the stream ends exactly. Otherwise the reader falls back to the fixed-size layouts that opentfraw 1.4.1 used, so every file 1.4.1 opens still opens.
- Every scan above MS1 has a precursor. For multiplexed (MSX) and multi-stage (EThcD, SPS-MS3) events it describes the first reaction; every reaction stays in `raw.scan_events`.
- `master_scan_number` comes from the "Master Scan Number" trailer, or else the nearest earlier scan one MS level lower. "Master Index" is a 0/1 flag and is never used.
- The isolation target comes from "Target M/Z" or the first reaction. "MS2 Isolation Offset" is an offset, exposed by `ScanParams::isolation_offset_mz`.
- A peak with a non-finite or non-positive m/z, or an intensity that is not finite or overflows `f32`, is dropped and counted in `dropped_peaks`. It does not fail the scan.
- The instrument status log is a time series: `RawFileReader::status_log_record` gives the last record written at or before a scan's start time, and `status_log_error` says why a log did not decode. Its records need `OpenOptions::with_records`.
- Unreadable scans return errors with operation, scan and byte-offset context (`Error::kind`, `Error::context`, `Error::root_cause`).
- Activation is reported only when qualified (v66 reaction tags 11 = HCD, 9 = ETD). Nondependent MS2 is not evidence of DIA; `is_dia` is a legacy heuristic.

## Examples

Six programs in `crates/sage-plus-raw/examples`, each taking a path to a `.raw` file. Run one with `cargo run --release --locked -p sage-plus-raw --example <name> -- <args>`.

- `dump <file.raw> [--max-scans N]` prints the file header, instrument, sequence row, run header, the first scan index entries and scan events, scan 1's parameters, the error log and the first scans' peaks: the quickest look at what a file contains.
- `strict_read <file.raw>` streams every spectrum and stops at the first decoding error, printing scan number, MS level and peak count per scan.
- `review_metadata <file.raw>` prints the trailer and scan-parameter addresses, the MS-level histogram and the first events' reactions and calibration coefficients.
- `reader_benchmark <file.raw> [REPEATS] [SEED] [--profile]` prints JSON timings for opening the file and for sequential and shuffled scan reads.
- `reuse_buffers <file.raw> [REPEATS] [SCANS] [--profile]` prints JSON timings for fresh-allocation against reused-buffer decoding (`try_extract_spectrum` against `try_extract_spectrum_into`), the pattern a downstream reader loop uses.
- `single_handle <file.raw> [REPEATS] [--profile]` prints JSON timings for two handles on the same path against one owned handle (`try_iter_spectra` against `iter_spectra_from_path`).

`reuse_buffers` and `single_handle` each assert that the two paths they compare decoded the same number of peaks, so they are a consistency check as well as a benchmark. The three timing examples take `--profile` as their last argument to decode profile spectra instead of centroids.

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
