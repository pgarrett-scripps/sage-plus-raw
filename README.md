# opentfraw-reader

Native Rust Thermo RAW reading without .NET or vendor DLLs. This independently maintained fork of [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw) preserves upstream Apache-2.0 history and adds qualified v66 decoding, scientific metadata corrections and strict extraction. The launch candidate is `opentfraw-reader 0.1.0-rc.1`.

This is publication preparation only. The repository remains private, no public crate has been published, and the manifests block registry publishing. The distinct crate name `opentfraw-reader` was unregistered when checked on 2026-09-26. That check does not reserve the name.

For current use, pin an authenticated Git dependency to a reviewed full commit. Use the Cargo dependency alias `opentfraw` to preserve existing Rust imports:

```toml
[dependencies]
opentfraw = { package = "opentfraw-reader", git = "https://github.com/pgarrett-scripps/opentfraw", rev = "FULL_REVIEWED_COMMIT" }
```

After an explicitly authorized public release, the planned registry dependency is:

```toml
[dependencies]
opentfraw = { package = "opentfraw-reader", version = "=0.1.0-rc.1" }
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
