# Thermo reader review and qualification

Public-launch preparation renames the package to `opentfraw-reader 0.1.0-rc.1` while retaining the Rust library name `opentfraw`. No public release has occurred. The qualification below was completed for the preceding private release and remains the supported scope.

Reviewed upstream OpenTFRaw 1.6.0 at `e300d34`, the private initial implementation at `5602402`, and consumer usage in MSConvert RS. The private release is `1.6.1-private.1`. This is a scoped parser and scientific metadata qualification, not a claim that every RAW revision has been validated.

## Fixed findings

| Priority | Finding | Resolution |
| --- | --- | --- |
| P1 | Average v66 event sizes shift metadata in mixed 232-byte and 288-byte streams. | Qualified per-record decoding validates fields and exact stream boundaries. |
| P1 | Four-reaction MSX events lose three reactions. | Preserve all four qualified reaction records and reject singular-precursor export. |
| P1 | Selected monoisotopic mass replaces the actual isolation target. | Keep explicit or reaction isolation targets separate from selected masses. |
| P1 | The byte-24 activation lookup calls the independent ETD acquisition HCD. | Decode the qualified v66 reaction tags, 9 for ETD and 11 for HCD. Other tags remain unknown. |
| P1 | Extraction and mzML writing silently skip corrupt peak arrays. | Add strict extraction and iteration, stop the adapter on first failure, and propagate writer errors. |
| P1 | Canonical writer defaults unknown activation to CID. | Reject mzML export when activation is unknown, before output begins. |
| P1 | SPS, MS3 and multiplexed precursors can be flattened into one precursor. | Reject convenience export for these modes and omit singular MS3 metadata. |
| P1 | Centroid arrays can be labeled profile and missing FT calibration can be treated as m/z. | Derive mode from returned arrays and require supported profile calibration. Reject nonfinite or invalid peak arrays. |
| P1 | Missing scan events can default to MS1 during strict extraction. | Strict extraction requires a recognized MS level for packet-format scans. |
| P2 | Several instrument CV IDs identify unrelated instruments. SRM chromatograms use the selected-ion-current term. | Correct all literal CV accession/name pairs against PSI-MS 4.2.2. Add instrument and chromatogram regressions. |
| P2 | SRM reports unsourced positive polarity and width 0.7. | Return unknown instead of invented values. |
| P2 | Reaction energy is labeled eV without unit evidence and cropped files retain dangling parent IDs. | Export labeled trailer energy only and omit parent IDs outside the run or not preceding the scan. |

## Activation evidence

The independent `BSA-FT-ETD.raw` and `BSA-FT-HCD.raw` acquisitions both have RAW version 66 and identical 136-byte preambles, including byte 24 equal to 4. Their first reaction tags differ: 9 for ETD and 11 for HCD. The corresponding vendor reference mzML files use `MS:1000598` and `MS:1000422`. All four MSX HCD reactions and the offset-isolation HCD fixture also use tag 11. This rules out a universal preamble-only lookup.

`ScanEvent::qualified_activation` accepts these tags for v66 events with the 136-byte preamble and a zero second reaction word. Every reaction must agree. Unknown or mixed tags return unknown. `ScanEventPreamble::activation` returns unknown because that byte alone is insufficient. Raw bytes and reaction words remain available.

The public references are pinned to [ProteoWizard commit 699dd489](https://github.com/ProteoWizard/pwiz/tree/699dd48953f81e4b14a02e5ae018a3f13c3b3b54/pwiz/data/vendor_readers/Thermo/Reader_Thermo_Test.data). The fixture downloader verifies each SHA-256 from `scripts/public-fixtures.json` before testing. Instrument and chromatogram terms were checked against [the official PSI-MS ontology](https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo), version 4.2.2 dated 2026-09-18.

## Verification

- Native formatting, Clippy with warnings denied, unit tests and documentation tests.
- Rust 1.85 native compilation and optional Python binding compilation.
- Six explicitly invoked public-fixture regressions covering paired activation references, centroid labeling, corrupt source errors, adapter error capture, profile calibration, MSX preservation, offset isolation, lossy export guards and instrument terms.
- Linux, macOS and Windows CI explicitly downloads and runs those public regressions. An absent fixture fails the explicitly invoked suite.
- Strict local extraction reads all 5,445 scans and 470,235 stored centroid peaks from the private large acquisition without loss. The acquisition has 3,445 MS1 scans and 2,000 MS2 scans.
- Earlier byte-level qualification showed all 5,445 event records and 470,235 peaks match the repaired MSConvert RS baseline.
- Private acquisition files are never committed. The large source SHA-256 is `b8f8fe0986a03388ac582bd8a49812a8ad8d78d3cec70bce01c345774a13e3b0`.

Local large-fixture invocation:

```sh
THERMO_LARGE_RAW=/path/to/large-vendor.raw \
THERMO_PUBLIC_CORPUS=/path/to/pwiz-thermo \
cargo test --locked -p opentfraw-reader --test local_metadata -- --ignored
```

## Remaining limits

- The legacy `extract_spectrum` and `iter_spectra` APIs still discard errors for compatibility. Scientific processing must use `try_extract_spectrum`, `try_iter_spectra`, or low-level Result-returning methods.
- The external `SpectrumSource` trait is infallible. Direct trait consumers must check `OpenTfRawSource::finish` after full iteration or canonical writing. Strict adapter writer methods perform that check automatically. Partial output must be discarded on error.
- SPS and MS3 precursor chains and cropped original-parent identities need a richer output model. Current convenience writers reject these unsupported modes. The public SPS and source-CID layouts fail during parsing.
- Nondependent MS2 alone does not establish DIA. `is_dia` remains a documented compatibility heuristic.
- Only documented v66 layouts and HCD/ETD reaction tags are qualified. Other uniform layouts retain upstream heuristics and require additional independent fixtures. Strict profile extraction rejects unavailable or unsupported calibration rather than treating frequencies as masses.
- Activation modes beyond qualified HCD and ETD remain unknown. Collision energies with unqualified units remain available only through raw fields.
- Python bindings are compile-checked but this release does not distribute Python wheels.

## Consumer adoption and provenance

Pin the private Git revision or release tag. Use the strict APIs and preserve unsupported-mode gates during consumer integration. Consumer migration is separate from parser release qualification.

Upstream history and Apache-2.0 attribution are retained. The upstream remote remains `https://github.com/Sigilweaver/OpenTFRaw`.
