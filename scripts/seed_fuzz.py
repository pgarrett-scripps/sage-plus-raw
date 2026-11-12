"""Seed fuzzers with verified public acquisitions so mutations reach decoders.

Every target also gets malformed copies of those acquisitions, derived from
them here rather than committed: the fixtures are well-formed files, so a
fuzzer seeded with them alone has to mutate its own way to the error paths
(`Error::AllocationTooLarge`, `UnexpectedEof`, `CorruptData`), and the
project prefers synthetic input for malformed cases anyway. The derivations
are deterministic, so a corpus is the same on every machine and in CI.
"""
import hashlib
import json
from pathlib import Path

# crates/sage-plus-raw/src/header.rs: `FileHeader::SIZE`, the file header at 0x0.
HEADER_SIZE = 1356
# A declared count or size no file can honour, which is what the allocation
# cap and the bounded reads exist for.
IMPLAUSIBLE = (0xFFFF_FFFF).to_bytes(4, "little")
TARGETS = ["open", "read_scan_peaks", "strict_decode"]


def truncations(data):
    """Copies that end inside a structure the reader has to bound-check."""
    return {
        # Inside the file header, and one field short of its end.
        "trunc-header-mid": data[: HEADER_SIZE // 2],
        "trunc-header-end": data[: HEADER_SIZE - 4],
        # Inside the method and run-header structures that follow it.
        "trunc-quarter": data[: len(data) // 4],
        # Inside the scan data, part way through a packet.
        "trunc-mid-packet": data[: len(data) * 3 // 5],
        # With the trailing scan index and trailer streams cut short.
        "trunc-tail": data[: max(0, len(data) - 1024)],
    }


def overwritten_fields(data):
    """Copies with a declared version, count or size field made implausible."""
    seeds = {}
    for name, offset, value in (
        # The file header's version word (crates/sage-plus-raw/src/header.rs:
        # magic u16, an 18-byte signature, four unused words, then version).
        ("version-unsupported", 36, (99).to_bytes(4, "little")),
        # The character count of the first Pascal string of the SeqRow that
        # follows the file header (crates/sage-plus-raw/src/seq_row.rs: its
        # InjectionData preamble is 64 bytes).
        ("count-first-string", HEADER_SIZE + 64, IMPLAUSIBLE),
    ):
        if offset < 0 or offset + 4 > len(data):
            continue
        seeds[name] = data[:offset] + value + data[offset + 4 :]
    return seeds


root = Path(__file__).resolve().parent.parent
manifest = json.loads((root / "scripts/public-fixtures.json").read_text())
written = 0
for name, expected in manifest["files"].items():
    if not name.endswith(".raw"):
        continue
    data = (root / "corpus/public" / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"Fixture digest mismatch: {name}")
    seeds = {name: data}
    stem = Path(name).stem
    for tag, seed in {**truncations(data), **overwritten_fields(data)}.items():
        seeds[f"{stem}.{tag}.raw"] = seed
    for target in TARGETS:
        destination = root / "crates/sage-plus-raw/fuzz/corpus" / target
        destination.mkdir(parents=True, exist_ok=True)
        for seed_name, seed in seeds.items():
            (destination / seed_name).write_bytes(seed)
            written += 1
print(f"seeded {written} files across {len(TARGETS)} targets")
