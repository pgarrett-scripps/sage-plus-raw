# Public launch readiness

Prepared on 2026-09-26. This is a release preparation record, not authorization to publish. The GitHub repository remains private and the real package manifest retains `publish = false`.

## Publication boundary

No public repository, registry upload, public documentation deployment or package-name reservation was performed. A package name appearing available does not reserve it. Existing private release tags and assets remain unchanged.

The final launch requires an explicit decision to make the repository public and publish the package. After that decision, follow `RELEASING.md`, verify the intended crates.io account and its email, enable its narrowly scoped publishing credential, recheck name availability and rerun the release checks. Never commit credentials. The registry performs server-side checks that a local dry run cannot prove.

## History and retained GitHub content

Gitleaks 8.30.1 was obtained from its official release and its archive checksum verified. A local scan of all reachable Git refs found no secret matches. This is a scanner result, not proof that every possible secret is absent. Upstream commit history, authorship and license notices are intentionally retained.

The existing private GitHub release and issue metadata were inspected and scanned without secret findings. Existing release assets contain code archives, checksums and validation reports. No raw acquisition was attached. There were no uploaded Actions artifacts at the time of inspection. Historical CI logs may contain ordinary checkout paths and public fixture URLs.

Changing visibility exposes repository history and retained GitHub content. Repeat this review immediately before launch if new commits, releases, issues or artifacts have been added. Keep private experimental acquisitions and user-provided data outside Git and release uploads.

## Supporting references

- [Cargo publication requirements and dry runs](https://doc.rust-lang.org/cargo/reference/publishing.html)
- [GitHub repository visibility changes](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/managing-repository-settings/setting-repository-visibility)
- [docs.rs build metadata](https://docs.rs/about/metadata)

## Thermo scope

The intended registry name is `opentfraw-reader`, separate from upstream `opentfraw`. Existing callers can retain the Rust import name through the documented Cargo dependency alias. The initial public candidate has its own version, independent of upstream releases and the earlier private tag.

The reviewed history contained no committed RAW or mzML acquisition files. Public regression inputs are fetched from checksum-pinned ProteoWizard URLs and are not redistributed inside the crate. Local qualification results record counts and hashes, not private acquisition contents.

Preserve the documented supported-layout limits and use strict Result-returning extraction for scientific processing. A public package release does not qualify additional acquisition modes or make a consumer migration scientifically equivalent.
