# Releasing sage-plus-raw

Releases are cut by pushing a `vX.Y.Z` tag. The [release workflow](.github/workflows/release.yml) then:

1. checks that the tag equals `v` + the `sage-plus-raw` version in `Cargo.toml`, and that `CHANGELOG.md` has a non-empty `## [X.Y.Z]` section;
2. runs clippy, the unit tests and the public-fixture tests on Linux, macOS and Windows;
3. runs `cargo publish` (skipped if that version is already on crates.io);
4. creates the GitHub release `vX.Y.Z`, with that version's CHANGELOG section as the notes (skipped if it already exists).

The crates.io upload comes before the GitHub release because it cannot be undone. If the workflow fails after publishing, re-run it: both steps are skipped when already done.

## Steps

1. Pick the version (Semantic Versioning; while below 1.0, a breaking API change bumps the minor version).
2. On a branch, set `version` under `[workspace.package]` in the root `Cargo.toml`, then run `cargo check` so `Cargo.lock` updates.
3. In `CHANGELOG.md`, rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD` and add a new empty `## [Unreleased]` above it. Check the notes with `sh scripts/changelog_section.sh X.Y.Z`.
4. Run the local checks from `CONTRIBUTING.md` plus `cargo publish --dry-run -p sage-plus-raw`, open a pull request, and merge it once CI is green.
5. Tag the merge commit on `main` and push the tag:

   ```sh
   git switch main && git pull
   git tag -a vX.Y.Z -m "sage-plus-raw X.Y.Z"
   git push origin vX.Y.Z
   ```

6. Watch the **Release** workflow. When it finishes, check <https://crates.io/crates/sage-plus-raw>, the GitHub release, and the docs.rs build at <https://docs.rs/sage-plus-raw>.

Only `CHANGELOG.md` sections above the "Upstream OpenTFRaw history" heading are used for release notes, so an upstream OpenTFRaw entry with the same number is never picked up by mistake.

## crates.io credentials

The workflow first tries [crates.io Trusted Publishing](https://crates.io/docs/trusted-publishing) through `rust-lang/crates-io-auth-action@v1`, which exchanges the job's GitHub OIDC token for a short-lived crates.io token. No secret is stored. If that fails, it falls back to the `CARGO_REGISTRY_TOKEN` repository secret.

Trusted Publishing can only be configured for a crate that already exists, so the first release needs a token:

1. **First release.** Create a crates.io API token (Account Settings, API Tokens) with the `publish-new` and `publish-update` scopes, limited to the `sage-plus-raw` crate name if possible, and a short expiry. Save it as the `CARGO_REGISTRY_TOKEN` repository secret (Settings, Secrets and variables, Actions). Then release as above.
2. **After the first release.** On crates.io, open the crate's Settings, Trusted Publishing, and add a GitHub publisher: owner `pgarrett-scripps`, repository `sage-plus-raw`, workflow filename `release.yml`, environment left empty.
3. Delete the `CARGO_REGISTRY_TOKEN` secret and revoke the token on crates.io. Later releases use Trusted Publishing only.

Optionally, restrict publishing further by requiring "Trusted Publishing only" in the crate's crates.io settings.

## If something goes wrong

- **Tag does not match the version.** Delete the tag (`git push origin :refs/tags/vX.Y.Z` and `git tag -d vX.Y.Z`), fix `Cargo.toml` on `main`, and tag again. Nothing was published.
- **A broken version was published.** crates.io versions cannot be deleted or replaced. Run `cargo yank --version X.Y.Z sage-plus-raw` and release a fixed patch version.
