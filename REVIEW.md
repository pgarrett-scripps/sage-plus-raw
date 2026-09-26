# Thermo reader review and changes

Reviewed upstream OpenTFRaw 1.6.0 at `e300d34`, together with MSConvert RS at `e6b9936`. This is a scoped parser and metadata review, not a claim that every RAW revision has been qualified.

## Fixed findings

| Priority | Finding | Resolution |
| --- | --- | --- |
| P1 | The v66 event reader guesses a uniform body from the average stream length. Mixed 232-byte and 288-byte events then shift later metadata and calibration. | Decode the qualified layout from each record. Validate every field and the exact final stream boundary. Unsupported variable layouts return an error. |
| P1 | The public FT-HCD-MSX fixture contains four 56-byte precursor reactions. The fallback parser retains only the first. | Decode all reactions in the qualified multiplexed layout. Validate reaction padding, mass range, calibration count and tail. |
| P1 | A monoisotopic selected mass can replace the actual isolation target in the public metadata API. | Prefer the explicit isolation target and event reaction target. A selected mass alone never establishes an isolation target. |
| P2 | SRM metadata always reports a 0.7 m/z isolation width, without source evidence. | Return unknown width instead of fabricating a measurement. |
| P1 | The mzML convenience writers can serialize just one reaction from a multiplexed scan. | Reject multiplexed export before writing any output. All raw reactions remain available on `scan_events`. The singular metadata precursor is unavailable for these scans. |

The default Cargo workspace target is the native Rust library. Python bindings remain an explicit optional workspace target. No .NET runtime or vendor DLL is introduced.

## Verification

- Formatting and Clippy with warnings denied pass.
- 116 unit tests and 5 documentation tests pass. The upstream conformance harness passes but skips its real decode when its optional corpus is absent.
- The explicitly invoked local regression passes on the 5,445-scan Thermo acquisition and public ProteoWizard fixtures.
- All 5,445 event records and 470,235 stored centroid peaks exactly match the MSConvert RS baseline using upstream OpenTFRaw plus its qualified v66 repair.
- The large acquisition has 3,445 MS1 scans and 2,000 MS2 scans.
- The independent MSX reference targets are 522.987607, 698.067188, 502.978512 and 768.099021. All four reactions now decode with width 5.0 and energy 27.0.
- The source-CID three-scan acquisition returns an unsupported-layout error instead of shifted event fields.
- Tests cover mixed event lengths, truncation, extra bytes, unsupported record fields, multiple reactions and target-versus-selected-mass semantics.

Local fixture invocation:

```sh
THERMO_LARGE_RAW=/path/to/large-vendor.raw \
THERMO_PUBLIC_CORPUS=/path/to/pwiz-thermo \
cargo test --locked -p opentfraw --test local_metadata -- --ignored
```

Private instrument acquisitions are not committed. The large source hash is `b8f8fe0986a03388ac582bd8a49812a8ad8d78d3cec70bce01c345774a13e3b0`.

## Remaining review findings and limits

- P1: Activation-code interpretation disagrees with the independent BSA-FT-ETD reference. Correct decoding needs revision and instrument-specific evidence. This branch does not invent a replacement activation value.
- P1: SPS/MS3 precursor chains and original parent identities for cropped files need a richer output model. MSConvert RS already rejects these acquisition modes. These guards must remain during integration.
- P1: The legacy spectrum and writer APIs can silently skip unreadable peak arrays. A follow-up should expose a strict error-returning extraction API and make the writer propagate it.
- P2: Nondependent MS2 does not by itself establish DIA. The upstream `is_dia` heuristic remains unsuitable as an acquisition classification guarantee.
- Only the documented v66 layouts are repaired. Other uniform layouts still use upstream heuristics and require broader fixture qualification.
- The multi-precursor mzML guard covers the convenience writer functions. Direct users of the upstream vendor-neutral adapter must honor the absent singular precursor or avoid unsupported acquisition modes.
- Linux validation is complete for this change. Windows and macOS execution have not been run locally.

## Integration

The package name remains `opentfraw`. Pin the private Git revision when adopting it in MSConvert RS. Keep the application's unsupported-mode gates and validate again before removing its compatibility repair. No consumer repository was changed by this reader work.

Upstream history and Apache-2.0 attribution are retained. The upstream remote remains `https://github.com/Sigilweaver/OpenTFRaw`.
