"""Validate registry packaging in a disposable committed copy without publishing."""
import hashlib
import io
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parent.parent
RELATIVE_MANIFEST = Path("crates/opentfraw/Cargo.toml")


def main():
    source_manifest = ROOT / RELATIVE_MANIFEST
    original = source_manifest.read_bytes()
    manifest = tomllib.loads(original.decode())
    if manifest["package"].get("publish") is not False:
        raise SystemExit("Refusing preflight without the real private publishing guard")
    if subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=normal"], cwd=ROOT
    ).strip():
        raise SystemExit("Commit the candidate first, preflight requires a clean checkout")
    snapshot = subprocess.check_output(["git", "archive", "--format=tar", "HEAD"], cwd=ROOT)
    try:
        with tempfile.TemporaryDirectory(prefix="opentfraw-registry-preflight-") as temporary:
            target = Path(temporary)
            with tarfile.open(fileobj=io.BytesIO(snapshot)) as archive:
                archive.extractall(target, filter="data")
            disposable_manifest = target / RELATIVE_MANIFEST
            copied = disposable_manifest.read_text()
            if copied.count("publish = false") != 1:
                raise SystemExit("Expected exactly one publishing guard in exported native manifest")
            disposable_manifest.write_text(copied.replace("publish = false", 'publish = ["crates-io"]', 1))
            environment = {
                key: value for key, value in os.environ.items()
                if key != "CARGO_REGISTRY_TOKEN"
                and not (key.startswith("CARGO_REGISTRIES_") and key.endswith("_TOKEN"))
            }
            environment["CARGO_TARGET_DIR"] = str(target / "target")
            subprocess.run(
                ["cargo", "publish", "--registry", "crates-io", "--dry-run", "--locked", "--allow-dirty", "-p", "opentfraw-reader"],
                cwd=target,
                env=environment,
                check=True,
            )
            archives = list((target / "target/package").glob("opentfraw-reader-*.crate"))
            if len(archives) != 1:
                raise SystemExit("Expected one verified registry candidate archive")
            package = archives[0]
            print(f"Verified disposable archive: {package.name}")
            print(f"SHA-256: {hashlib.sha256(package.read_bytes()).hexdigest()}")
            print("Dry-run completed, no upload or source-manifest change performed")
    finally:
        if source_manifest.read_bytes() != original:
            raise SystemExit("Real source manifest changed during registry preflight")


if __name__ == "__main__":
    main()
