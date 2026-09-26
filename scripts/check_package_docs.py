"""Check that native packages retain their license, provenance and support scope."""
from pathlib import Path
import subprocess
import tomllib

ROOT = Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "crates/opentfraw"
DOCUMENTS = ("LICENSE", "NOTICE", "ATTRIBUTION.md", "CONTRIBUTORS.md", "SECURITY.md", "REVIEW.md", "RELEASING.md")


def main():
    for name in DOCUMENTS:
        if (ROOT / name).read_bytes() != (PACKAGE / name).read_bytes():
            raise SystemExit(f"Packaged document differs from workspace document: {name}")
    for directory in ("opentfraw", "opentfraw-py"):
        manifest = tomllib.loads((ROOT / "crates" / directory / "Cargo.toml").read_text())
        if manifest["package"].get("publish") is not False:
            raise SystemExit(f"Private publishing guard missing: {directory}")
    files = subprocess.check_output(
        ["cargo", "package", "--locked", "--offline", "--allow-dirty", "--list", "-p", "opentfraw-reader"],
        cwd=ROOT,
        text=True,
    ).splitlines()
    missing = set(DOCUMENTS + ("README.md",)) - set(files)
    if missing:
        raise SystemExit(f"Required package documents missing: {sorted(missing)}")
    print("Package inventory includes license, attribution, support scope and release instructions")


if __name__ == "__main__":
    main()
