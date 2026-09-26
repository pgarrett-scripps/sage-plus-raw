# Private releases

This repository distributes private Git tags and verified `.crate` archives through private GitHub releases. It does not publish to crates.io or PyPI.

1. Set the shared version in `Cargo.toml` to a distinct private version, such as `1.6.1-private.1`, and update `Cargo.lock` and `CHANGELOG.md`.
2. Run formatting, Clippy, native tests, public fixture tests, Python compilation and `cargo package --locked -p opentfraw`.
3. Commit and push. Run `scripts/check-release-ready.sh` after Reader checks succeeds for that exact commit. Its matrix includes Linux, macOS, Windows and the minimum Rust version.
4. Run `cargo package --locked -p opentfraw` from the clean checkout. This embeds the committed source revision in `.cargo_vcs_info.json`.
5. Compute a SHA-256 checksum for `target/package/opentfraw-VERSION.crate` and retain it with the archive.
6. Create an annotated `vVERSION` tag at the verified commit and upload the archive and checksum to a private GitHub release. Include the supported-mode limits from `REVIEW.md` in the release notes.
7. Verify the release assets and repository visibility. Pin consumers to the reviewed tag or full commit hash.

`reader.yml` does not publish artifacts automatically. Public registry names and upstream Python packaging instructions do not grant registry publishing authorization.
