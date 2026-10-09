"""Seed fuzzers with verified public acquisitions so mutations reach decoders."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent.parent
manifest = json.loads((root / "scripts/public-fixtures.json").read_text())
for name, expected in manifest["files"].items():
    if not name.endswith(".raw"):
        continue
    data = (root / "corpus/public" / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"Fixture digest mismatch: {name}")
    for target in ["open", "read_scan_peaks", "strict_decode"]:
        destination = root / "crates/sage-plus-raw/fuzz/corpus" / target
        destination.mkdir(parents=True, exist_ok=True)
        (destination / name).write_bytes(data)
