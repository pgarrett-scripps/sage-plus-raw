# opentfraw-reader

Native Rust Thermo RAW reading without .NET or vendor DLLs. This independently maintained fork of [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw) preserves upstream Apache-2.0 history and adds qualified v66 decoding, scientific metadata corrections and strict extraction. The launch candidate is `opentfraw-reader 0.1.0-rc.2`.

This is publication preparation only. The repository remains private, no public crate has been published, and the manifests block registry publishing. The distinct crate name `opentfraw-reader` was unregistered when checked on 2026-09-26. That check does not reserve the name.

For current use, pin an authenticated Git dependency to a reviewed full commit. Use the Cargo dependency alias `opentfraw` to preserve existing Rust imports:

```toml
[dependencies]
opentfraw = { package = "opentfraw-reader", git = "https://github.com/pgarrett-scripps/opentfraw", rev = "FULL_REVIEWED_COMMIT" }
```

After an explicitly authorized public release, the planned registry dependency is:

```toml
[dependencies]
opentfraw = { package = "opentfraw-reader", version = "=0.1.0-rc.2" }
```

The registry example is not available yet. The upstream `opentfraw` package is a separate project. No fork docs.rs page or hosted documentation site is claimed to exist. Generate API documentation locally with `cargo doc --locked -p opentfraw-reader --no-deps --open`.

## Reading and exporting

Open with `RawFileReader::open_path`. Use `try_extract_spectrum` for one zero-based scan index or `try_iter_spectra` for a stream of `Result<SpectrumRecord>` values. Both propagate unreadable arrays and validate emitted peak values. Empty valid spectra remain valid. The low-level `read_scan_peaks` API remains available. Run the checked streaming example with `cargo run --locked -p opentfraw-reader --example strict_read -- acquisition.raw`.

`scan_metadata` reads metadata without reading peaks. Precursor targets and selected monoisotopic masses remain separate. Centroid arrays report centroid mode even when their source acquisition was profile mode. Profile extraction requires a supported calibration shape, or explicitly identified ion-trap m/z coordinates.

The legacy `extract_spectrum` and `iter_spectra` functions remain available for compatibility. They discard errors and must not be used for lossless conversion.

`write_mzml`, `write_indexed_mzml` and the matching `OpenTfRawSource` methods propagate input and output failures. Discard output on any failure because it may be partial. The external `openmassspec_core::SpectrumSource` trait cannot return errors, so users calling its iterator or canonical writer directly must call `OpenTfRawSource::finish` afterward. The adapter stops at its first failure.

## Qualified scope and limits

- Qualified v66 event layouts include mixed 232-byte and 288-byte records and the documented four-reaction MSX layout.
- Independent public HCD and ETD references distinguish v66 reaction tags 11 and 9. The preamble activation byte alone cannot distinguish these acquisitions. Other tags remain unknown. Reconstructed filters include activation only when qualified.
- Multiple reactions remain available in `scan_events`. A singular precursor is unavailable for multiplexed scans and MS3 or higher.
- Convenience mzML export rejects MSX, SPS, MS3 or higher and unknown activation. Its shared writer would otherwise invent CID for unknown activation.
- Nondependent MS2 is not sufficient evidence of DIA. The legacy `is_dia` boolean is only that heuristic, not an acquisition classification.
- Collision energy is exported only from labeled trailer values. Unqualified reaction energy units are not guessed. SRM polarity and isolation width remain unknown when not sourced.
- Parent scan IDs outside the current run or referring to the current or a later scan are omitted. Full cropped-file parent identity reconstruction is not implemented.
- Unsupported SPS and source-CID fixture layouts return errors. Other upstream format revisions remain available but are not comprehensively qualified by this private release.

See [REVIEW.md](REVIEW.md) for evidence and remaining limits. Original upstream format documentation is in `docs/docs`, and may describe behavior superseded by this private branch.

## Validation

```sh
cargo fmt --all --check
cargo clippy --locked -p opentfraw-reader --all-targets -- -D warnings
cargo test --locked -p opentfraw-reader
python3 scripts/fetch_public_fixtures.py
THERMO_PUBLIC_CORPUS=corpus/public cargo test --locked -p opentfraw-reader --test public_fixtures -- --ignored
cargo package --locked -p opentfraw-reader
```

CI runs native checks, public fixture regressions, Python binding compilation and package verification on Linux, macOS and Windows. Rust 1.85 is checked separately. The optional Python binding retains the import `opentfraw` and uses a private candidate distribution identity `opentfraw-reader`. Python wheels are not published or included in this launch scope.

Private large acquisitions are never committed. The local regression additionally exercises the documented 5,445-scan acquisition. See [RELEASING.md](RELEASING.md) for safe registry preflight and the later launch gate.

## Strict arrays and contextual errors

`try_extract_spectrum_arrays` returns `SpectrumArrays` with `f64` m/z and intensity values, metadata and the effective scan mode. It preserves full decoded profile precision. `try_extract_spectrum` and `try_iter_spectra` return the existing `SpectrumRecord` with `f32` intensities and reject nonfinite narrowing. The lossy `extract_spectrum` and `iter_spectra` compatibility helpers are deprecated.

`Error::kind()` distinguishes `Io`, `Corrupt`, `Unsupported`, `InvalidInput` and `ResourceLimit` (reserved for future explicit budgets). Input-length allocation guards are corruption. `Error::context()` identifies the operation and available input path, scan index, scan number and packet offset. Nested `std::error::Error::source()` chains preserve the original cause. Packet offsets are absolute file positions. Inner decoder offsets may be relative to the packet. `Error::root_cause()` permits typed leaf inspection without parsing text. Source and sink failures stop exports. Call `OpenTfRawSource::finish()` after using its infallible external trait iterator.

Python `peaks`, `profile` and `scan` use strict decoding. `scan` labels its emitted centroid arrays as centroid, while the filter string retains the recorded acquisition mode. Errors map to `OSError`, `ValueError` for corruption, `NotImplementedError` for unsupported modes, `IndexError` for invalid scans, and `MemoryError` for explicit resource budgets if introduced. Profile output remains `float64`.

## Qualification and performance

CI checks independent ProteoWizard centroid and profile arrays for four pinned public acquisitions. These qualify existing layouts rather than additional instrument families. Profile comparisons remove reference-only zero padding and use the reference arrays' declared 32-bit precision. See [REVIEW.md](REVIEW.md) for support limits.

Run `cargo run --release --locked -p opentfraw-reader --example reader_benchmark -- acquisition.raw 3 42 --profile` for JSON timings of repeated open, sequential and deterministically shuffled reads. Omit `--profile` for centroids. The executable reports peak counts and fails on any decoding error. It does not log the acquisition path. Compare the same input, build flags, machine and cache conditions. Opening warms the filesystem cache before scan timings. Measure peak RSS externally if needed.

The bounded fuzz CI job uses public verified seeds, three targets, address sanitizer and a fixed seed. Each target runs for 20 seconds with a 1 GiB RSS cap and 5-second input timeout. It uploads no corpus or crash artifacts. To reproduce, run `python scripts/fetch_public_fixtures.py`, `python scripts/seed_fuzz.py`, then `cargo +nightly-2026-09-26 fuzz run strict_decode -- -max_total_time=20 -timeout=5 -rss_limit_mb=1024 -max_len=1048576 -seed=42` from `crates/opentfraw`. Bounded smoke fuzzing is regression coverage, not exhaustive assurance.

FAIMS compensation voltage is emitted only when the trailer explicitly enables FAIMS and the voltage is finite. Enabled zero voltage is retained. Raw trailer values remain inspectable when FAIMS is disabled or its state is unknown.
