"""Fetch small public regression acquisitions at an immutable upstream revision."""
import hashlib
import json
from pathlib import Path
import urllib.request

root = Path(__file__).resolve().parent.parent
manifest = json.loads((root / "scripts/public-fixtures.json").read_text())
destination = root / "corpus/public"
destination.mkdir(parents=True, exist_ok=True)
base = (
    "https://raw.githubusercontent.com/ProteoWizard/pwiz/"
    + manifest["revision"]
    + "/pwiz/data/vendor_readers/Thermo/Reader_Thermo_Test.data/"
)
for name, expected in manifest["files"].items():
    path = destination / name
    if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == expected:
        continue
    with urllib.request.urlopen(base + name, timeout=60) as response:
        data = response.read()
    if hashlib.sha256(data).hexdigest() != expected:
        raise RuntimeError(f"Fixture digest mismatch: {name}")
    path.write_bytes(data)
    print(f"Verified {name}")
