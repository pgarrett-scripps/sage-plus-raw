# Thermo reader refactor

This change separates responsibilities that previously accumulated in `reader.rs` and `mzml.rs`. It preserves the public API, package identity, qualification limits and publishing guards.

## Responsibility boundaries

- `binary_reader.rs` owns bounded byte reads, seek bookkeeping, string decoding and their existing regression tests. Format decoders now depend directly on these helpers instead of the file orchestrator.
- `reader.rs` owns opening RAW files, locating their structures and reading scan packets. Binary layout decisions and scan-index behavior are unchanged.
- `metadata.rs` owns typed views of trailer and instrument-status records. These accessors no longer share a file with allocation guards and scan packet I/O.
- `spectrum.rs` owns metadata derivation, strict spectrum validation and compatibility extraction helpers. It has no writer dependency. The same extraction functions and record types remain available through both the crate root and `mzml` re-exports.
- `mzml.rs` owns format translation, chromatograms, adapter error handling and export gates. Unsupported acquisition checks still happen before convenience writers begin output.

The private positional spectrum tuple is now `DecodedSpectrum`, with named metadata, peak-array and emitted-mode fields. This makes the handoff from decoding to XML translation explicit. Only this handoff and its decoder are visible within the crate. No new public module or generic framework was introduced.

Profile extraction consumes already validated pairs in one split pass. It no longer repeats the positive-mass filter for each output array. The validation still runs before splitting, and neither the calibration formula nor point ordering changes. The unused nominal-mode argument was removed from the private array resolver because emitted mode depends on the arrays actually returned.

## Validation

Existing binary allocation and malformed-input regressions moved with their implementations. The typed-accessor coverage check now reads the metadata and spectrum modules, preserving its original coverage requirement. Public fixture regressions exercise the same scientific results and export failures through the existing API paths.

Run the normal native tests, public fixtures, Clippy, Python binding compilation, strict Rust documentation build, package inventory check, package verification and disposable registry preflight. Consumer compatibility and before/after fixture comparisons supplement these checks. No claim of broader acquisition support or performance improvement is made by this structural refactor.

The before/after differential against commit `7998a20` matched across eight Thermo acquisitions. It compared metadata, events, chromatograms, bitwise peak arrays, strict errors, and regular and indexed mzML bytes or error results. The large acquisition exercised 10,890 centroid and profile spectrum attempts and returned 3,238,475 peaks across those modes without a difference.

The unchanged consumer API also passed the isolated MSConvert suite with 111 tests passing and 10 expected skips, plus the Sage strict-reader suite with 10 tests passing and one expected skip. Native validation passed 117 unit tests, six documentation tests, six public-fixture regressions and the local large-acquisition regression. Rust 1.85, Python binding compilation, strict Rust documentation, Clippy and package verification passed.
