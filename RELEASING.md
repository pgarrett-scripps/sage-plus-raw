# Public launch preparation

The candidate package is `opentfraw-reader 0.1.0-rc.1`. The Rust library name remains `opentfraw`. Existing private tags and artifacts retain their original package identity and version.

Nothing in these instructions authorizes publication or a visibility change. Real manifests retain `publish = false`, and Python retains its private classifier. There is no automatic publishing or website deployment workflow.

## Safe preflight

From a clean, committed checkout:

```sh
python3 scripts/registry_preflight.py
```

The script verifies the real publishing guards, exports the committed repository into a temporary directory, and changes the native crate publishing guard only in that disposable copy. It runs `cargo publish --dry-run --locked` for the candidate package. Cargo packages and verifies the crate but does not upload it. The temporary copy is deleted afterward, and the real manifests are checked again. No token is requested or written.

The native `.crate` must contain LICENSE, NOTICE, ATTRIBUTION.md, CONTRIBUTORS.md, SECURITY.md, REVIEW.md and RELEASING.md alongside README and source. `scripts/check_package_docs.py` checks that those copies match the repository documents. To build a reviewable local archive without registry access:

```sh
cargo package --locked --offline -p opentfraw-reader
```

Only use `--offline` when the required dependencies and index entries are already cached. Package verification from a clean checkout records the commit in `.cargo_vcs_info.json`. Retain the archive SHA-256 with qualification evidence. CI runs package verification, API documentation and the isolated registry preflight without a registry token.

## Launch gate requiring later authorization

Before any public launch, review the exact committed candidate and its CI results. Recheck crate-name availability because an availability check does not reserve a name. Confirm GitHub private vulnerability reporting and a working maintainer contact route. Review source, history, fixtures, licenses and dependency advisories for public distribution. Confirm the support limits in README and REVIEW.

The account owner must then explicitly authorize public repository visibility and registry publication. Only after that authorization should maintainers remove the native publishing guard in a reviewed commit, configure registry ownership and a scoped authentication mechanism, publish the approved version, and verify its registry metadata and generated documentation. Website deployment and Python publication are separate decisions requiring their own review and authorization.

Do not use the upstream `opentfraw` crates.io or PyPI identities or assume access to upstream publisher credentials. Do not modify existing private release tags. A future public publishing workflow should use least-privilege permissions and an explicit approval environment.
